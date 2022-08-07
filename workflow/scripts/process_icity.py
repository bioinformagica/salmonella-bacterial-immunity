#!/usr/bin/env python3
# Imports
# native
import os
import sys
import re
from collections import Counter
from collections import defaultdict
import random
from itertools import chain
import multiprocessing

# external
import numpy as np
import modin.pandas as pd
from modin.config import Engine
Engine.put("ray")

os.environ["MODIN_CPUS"] = str(int(multiprocessing.cpu_count() * 0.40))

def eprint(*args):
    print(*args, file=sys.stderr)

def get_sequential_groups(vals: iter) -> list:
    """ This functions take as input a SORTED iterables with UNIQUE values
        and returns a list with all groups that are sequential labaled by a
        interger index.
        Ex:
        $ get_sequential_groups([1528, 1529, 1530, 1531, 1559, 1560, 1561, 1562, 1565, 1566])
        [1, 1, 1, 1, 2, 2, 2, 2, 3, 3]
    """
    groups_index = 1
    groups = []
    stack = []
    for i in vals:

        if not stack:
            stack.append(i)
            groups.append(groups_index)
            continue

        if i == stack[-1] + 1:
            stack.append(i)
            groups.append(groups_index)

        else:
            groups_index += 1
            groups.append(groups_index)
            stack = [i]

    return groups


def check_neighborhood(all_items: iter, range_to_check: int = 10) -> list:
    if not any(all_items):
        for i in all_items:
            yield None
        return

    all_indexes = tuple(enumerate(all_items))
    baits_location = sorted(set(i for i, x in all_indexes if x))

    for indx, _ in all_indexes:

        if indx in baits_location:
            yield 'is_bait'
            continue

        if min(map(lambda x: abs(indx - x), baits_location)) <= range_to_check:
            yield 'is_close'
            continue

        yield None


def get_distance_from_baits(df) -> tuple:
    bait_coords = set(df.query('bait_assign == "is_bait"').start.values.tolist())

    for start_coord, bait_assign in df.values.tolist():
        if not bait_assign:
            yield None
            continue

        if start_coord in bait_coords:
            yield 0
            continue

        yield min(map(lambda x: abs(start_coord - x), bait_coords))

def get_name_percentage(info_list):
    total_info = len(info_list)
    name_counts = Counter(re.findall('(?<=Name\=).+?(?=\;)' ,';'.join(info_list)))
    total_count = sum(name_counts.values())

    results = []
    for name, count in name_counts.items():
        results.append('{}={}%'.format(name, round((count/total_info)*100, 2)))

    if not total_count == total_info:
        results.append('Hypothetical protein={}%'.format(round(((total_info - total_count)/total_info)*100, 2)))

    return ';'.join(results)


def get_name_systems(info_list):
    name_counts = Counter(re.findall('(?<=system\=).+?(?=\;)', ';'.join(info_list)))

    results = []

    if not name_counts:
        return None

    for name, count in name_counts.items():
        results.append('{}'.format(name))

    return ';'.join(results)


def get_edges(topology: iter):
    for indx, cluster_id in enumerate(topology):
        try:
            yield (cluster_id, topology[indx + 1])
        except IndexError:
            pass


def assign_loci_id(bait_assign_iter):
    current_loci_index = 1
    is_inside_loci = False

    for description in bait_assign_iter:

        if not description and is_inside_loci:
            current_loci_index += 1

        if not description:
            is_inside_loci = False
            yield None
            continue

        is_inside_loci = True
        yield 'loci_{}'.format(str(current_loci_index))


def most_common_name(annotations: iter):
    for annotation in annotations:
        if not annotation:
            yield 'NA'
            continue

        yield sorted(map(lambda x: x.split('='), annotation.split(';')), key = lambda x: float(x[1].replace('%', '')), reverse=True)[0][0]


def main(all_cds_location, dfs_prediction, cluster_file):
    # Cytoscape variables
    max_diameter_size = 230
    min_diameter_size = 5
    max_edge_size = 60
    min_edge_size = 5
    max_border_size = 60
    min_border_size = 1
    min_bait_count_to_be_bait = .2
    default_color = '#d4ecff'

    ###### Merging bait prediction with cds info #######
    eprint('Reading cds and bait tables ...')
    cds_df = (pd.read_csv(all_cds_location, header=None, sep='\t')
              .rename(columns = dict(zip(range(6), ['locus_tag', 'contig_id', 'start', 'end', 'strand', 'info'])))
              .assign(contig_index = lambda df: df.contig_id.apply(lambda contig_id: int(contig_id.split('_')[1])))
              .sort_values(['contig_id', 'contig_index', 'start'])
              .drop('contig_index', axis=1)
              .assign(csd_lenght = lambda df: df.end - df.start )
              .merge(
                  pd.read_csv(dfs_prediction, header=None, sep=' ').rename(columns={0: 'locus_tag', 1: 'dfs_prediction'}),
                  on=['locus_tag'],
                  how='left'
              )
              .merge(
                  pd.read_csv(cluster_file, header=None, sep='\t').rename(columns={0: 'seeds', 1: 'locus_tag'}),
                  on=['locus_tag'],
                  how='inner'
              )
              .groupby('contig_id', as_index=False)
              .apply(lambda gdf: (gdf
                                  .assign(bait_assign = lambda x: tuple(check_neighborhood(x.dfs_prediction.fillna(0).values)))
                                  .assign(distance_to_bait = lambda x: tuple(get_distance_from_baits(x.loc[:, ['start', 'bait_assign']])))
                                  .assign(loci_id = lambda x: tuple(assign_loci_id(x.bait_assign)))
                                  ))
              )

    eprint('Writing cds table ...')
    cds_df.to_csv('cds_df.tsv', sep='\t', index=False)

    eprint('Calculating seeds df ...')
    seeds_df = (cds_df.groupby('seeds')
                .agg({
                    'seeds'           : 'count',
                    'bait_assign'     : [
                        lambda x: x.where(lambda x: x != None).dropna().shape[0],
                        lambda x: x.where(lambda x: x == 'is_bait').dropna().shape[0],
                    ],
                    'locus_tag'       : lambda x: ','.join(x.values.tolist()),
                    'distance_to_bait': 'mean',
                    'info'            : lambda x: get_name_percentage(x.values),
                    'dfs_prediction'  : lambda x: get_name_systems(x.fillna('NA')),
                    'csd_lenght'      : 'mean',
                })
                .droplevel(0, axis=1))

    seeds_df.columns = ['members_count', 'members_close_to_baits', 'is_bait_count', 'locus_tags', 'mean_distance_to_bait', 'members_names', 'dfs_names', 'cds_lenght_mean']


    seeds_df = (seeds_df
                .query('members_close_to_baits > 0')
                .assign(icity = lambda df: df.members_close_to_baits / df.members_count)
                .sort_values('mean_distance_to_bait', ascending=False)
                [['members_count', 'icity', 'mean_distance_to_bait', 'members_close_to_baits', 'is_bait_count', 'locus_tags', 'members_names', 'dfs_names', 'cds_lenght_mean']]
                .assign(is_bait_perc = lambda df: df.is_bait_count / df.members_count)
                [['members_count', 'icity', 'mean_distance_to_bait', 'members_close_to_baits', 'locus_tags', 'members_names', 'is_bait_count', 'is_bait_perc', 'dfs_names', 'cds_lenght_mean']]
                .reset_index())


    dfs_dict = (cds_df
                .groupby(['contig_id', 'loci_id'])
                .apply(
                    lambda df: df.assign(group_systems = tuple(assign_loci_id(df.fillna(0).dfs_prediction)))
                )
                .query('group_systems.notna()')
                .groupby(['contig_id', 'loci_id', 'group_systems'])
                .agg({'seeds': lambda x: sorted(tuple(x.values))})
                .droplevel('group_systems')
                .reset_index()
                .groupby(['contig_id', 'loci_id'])
                .agg({'seeds': lambda x: x.values.tolist()})).seeds.to_dict()

    dfs_count = len(set(map(tuple, chain.from_iterable(dfs_dict.values()))))

    seeds_df = seeds_df.merge(
        (cds_df.loc[:,['seeds', 'contig_id', 'loci_id']]
         .query('loci_id.notna()')
         .drop_duplicates()
         .assign(tup_query = lambda df: tuple(map(lambda x: tuple(x), df[['contig_id', 'loci_id']].values.tolist())))
         .assign(diversity_count = lambda df: df.tup_query.apply(lambda x:  tuple(map(lambda y: tuple(y), dfs_dict.get(x, None))) ))
         .groupby('seeds')
         .agg({'diversity_count': lambda x: len(set(chain.from_iterable(x.values)))})
         .assign(diversity_score = lambda df: df.diversity_count / dfs_count)),
        on=['seeds'],
        how='inner'
    )[['seeds', 'members_count', 'icity', 'diversity_score', 'mean_distance_to_bait', 'members_close_to_baits', 'locus_tags', 'members_names', 'is_bait_count', 'is_bait_perc', 'dfs_names', 'cds_lenght_mean']]

    eprint('Writing seeds df ...')
    seeds_df.to_csv('seeds_df.tsv', sep='\t', index=False)

    eprint('Calculating nodes df ...')
    max_members_count = seeds_df.members_count.max()
    color_pallet = defaultdict(lambda: "#%06x" % random.randint(0, 0xFFFFFF))
    nodes = (seeds_df
             .drop('locus_tags', axis=1)
             .rename(columns={'seeds': 'id'})
             .assign(node_size = lambda df: df.members_count.apply(lambda x: int((max_diameter_size * x)/max_members_count)))
             .assign(node_shape = lambda df: df.is_bait_perc.apply(lambda x: 'TRIANGLE' if x > min_bait_count_to_be_bait else 'ELLIPSE'))
             .assign(dfs_label = lambda df: df.dfs_names.apply(lambda x: None if not x else x.upper().split(';')[-1].split('_')[0]))
             .assign(node_color = lambda df: df.dfs_label.apply(lambda x: default_color if not x else color_pallet[x]))
             .assign(node_color = lambda df: tuple(map(lambda x: default_color if x == 'ELLIPSE' else x[1], df[['node_shape', 'node_color']].values.tolist())))
             .assign(border_size = lambda df: df.icity.apply(lambda x: int(x * max_border_size)))
             .assign(border_size = lambda df: tuple(map(lambda x: min_border_size if x[0] == 'TRIANGLE' else x[1], df[['node_shape', 'border_size']].values.tolist()))))

    eprint('Writing nodes df ...')
    nodes.to_csv('nodes_df.tsv', sep='\t', index=False)

    eprint('Calculating edges df ...')
    edges = (cds_df
             .loc[:,['contig_id','seeds', 'bait_assign']]
             .assign(loci_id = lambda df: tuple(assign_loci_id(df.bait_assign)))
             .query('bait_assign.notna()')
             .groupby(['contig_id', 'loci_id'])
             .agg({
                 'seeds' : lambda series: tuple(get_edges(series.values))
             })
             .explode('seeds')
             .reset_index()
             .loc[:,'seeds']
             .apply(lambda edge: sorted(edge, key=lambda seed_id: int(seed_id.split('_')[1])))
             .to_frame()
             .assign(source = lambda df: df.seeds.apply(lambda x: x[0]))
             .assign(target = lambda df: df.seeds.apply(lambda x: x[1]))
             .drop('seeds', axis=1)
             .groupby(['source', 'target'], as_index=False)
             .size()
             .sort_values('size', ascending=False)
             .rename(columns={'size': 'weight'})
             .reset_index(drop=True))


    max_edge_count = edges.weight.max()
    edges = edges.assign(edge_size = lambda df: df.weight.apply(lambda x: max(min_edge_size, int((max_edge_size * x)/max_edge_count))))

    eprint('Writing edges df ...')
    edges.to_csv('edges_df.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main(*sys.argv[1:])
