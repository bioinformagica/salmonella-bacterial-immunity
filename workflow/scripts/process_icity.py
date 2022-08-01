#!/usr/bin/env python3
# Imports
# native
import os
import sys
import re
from collections import Counter
import multiprocessing

# external
import numpy as np
import modin.pandas as pd
from modin.config import Engine
Engine.put("ray")

os.environ["MODIN_CPUS"] = str(int(multiprocessing.cpu_count() * 0.75))


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
    name_counts = Counter(re.findall('(?<=system\=).+?(?=\;)' ,';'.join(info_list)))

    results = []
    for name, count in name_counts.items():
        results.append('{}'.format(name))

    return ';'.join(results)


def main(all_cds_location, dfs_prediction, cluster_file, output_fname):
    eprint('Reading prokka and dfs prediction ...')
    cds_df = (pd.read_csv(all_cds_location, header=None, sep='\t')
              .rename(columns = dict(zip(range(6), ['locus_tag', 'contig_id', 'start', 'end', 'strand', 'info'])))
              .assign(contig_index = lambda df: df.contig_id.apply(lambda contig_id: int(contig_id.split('_')[1])))
              .sort_values(['contig_id', 'contig_index', 'start'])
              .drop('contig_index', axis=1)
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
                                  ))
              )

    eprint('Calculating metrics ...')
    seeds_df = (cds_df.groupby('seeds')
                .agg({
                    'seeds'           : 'count',
                    'bait_assign'     : [
                        lambda x: x.dropna().shape[0],
                        lambda x: x.where(lambda x: x == 'is_bait').dropna().shape[0]
                    ],
                    'locus_tag'       : lambda x: ','.join(x.values.tolist()),
                    'distance_to_bait': 'mean',
                    'info'            : lambda x: get_name_percentage(x.values),
                    'dfs_prediction'  : lambda x: get_name_systems(x.fillna('NA'))
                })
                .droplevel(0, axis=1))

    seeds_df.columns = ['members_count', 'members_close_to_baits', 'dfs_count', 'locus_tags', 'mean_distance_to_bait', 'members_names', 'dfs_names']

    seeds_df = (seeds_df
                .query('members_close_to_baits > 0')
                .assign(icity = lambda df: df.members_close_to_baits / df.members_count)
                .sort_values('mean_distance_to_bait', ascending=False)
                .assign(diversity_score = lambda df: df.dfs_count / df.dfs_count.sum())
                [['members_count', 'icity', 'diversity_score', 'mean_distance_to_bait', 'members_close_to_baits', 'dfs_count', 'locus_tags', 'members_names', 'dfs_names']]
                .assign(is_known = lambda df: df.members_count == df.dfs_count)
                .sort_values(['icity', 'diversity_score', 'mean_distance_to_bait'], ascending=[False, False, True])
                [['members_count', 'icity', 'diversity_score', 'mean_distance_to_bait', 'members_close_to_baits', 'locus_tags', 'members_names', 'dfs_names', 'is_known']]
                .reset_index()
                )

    seeds_df.to_csv(output_fname, sep='\t', index=False)

    eprint('All jobs finished !')

if __name__ == '__main__':
    main(*sys.argv[1:])
