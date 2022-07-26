#!/usr/bin/env python3
# native
from collections import defaultdict
from collections import Counter
from itertools import repeat
from operator import itemgetter
import sys

# external
import pandas as pd
import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt

def get_percentage_of_member_count(df):
    results = []
    for total, count_dict in df.loc[:,['members_count', 'members_names']].values.tolist():
        record = []
        for name, count in count_dict.items():
            record.append('{}:{}%'.format(name, round((count/total) * 100, 2)))
        results.append(','.join(record))
    return results

def get_bait_distances(main_df):
    temp = []
    for i, df in main_df.query('casset_id.notna()').groupby(['contig_id', 'casset_id']):
        try:
            if i[0] == '2052597.6_10':
                temp.append(df.assign(closest_bait_distance = lambda df: df.start.apply(lambda query_start: min(map(
                     lambda bait_start: abs(query_start - bait_start),
                     df.query('dfs_prediction == dfs_prediction').start.values.tolist()
                )))))
        except:
            pass
    return pd.concat(
        map(
            lambda gdf: pd.DataFrame(
                data=[(gdf[0], gdf[1].closest_bait_distance.mean(),  gdf[1].closest_bait_distance.std(), gdf[1].closest_bait_distance.std()/gdf[1].closest_bait_distance.mean())],
                columns=('seeds', 'mean_closest_bait_distance', 'std_closest_bait_distance', 'cv_closest_bait_distance')
            ),
            pd.concat(temp)[['seeds', 'closest_bait_distance']].groupby('seeds')
        )
    )

def read_prokka(prokka_gff, seeds_with_close_dfs, casset_info):
    return (pd.read_csv(prokka_gff, sep='\t', header=None)
            .rename(columns=dict(zip(range(7), ['contig_id', 'locus_tag', 'start', 'end', 'strand', 'info', 'dfs_prediction'])))
            .merge(
                pd.read_csv(seeds_with_close_dfs, header=None, sep='\t').rename(columns=dict(zip(range(2), ['seeds', 'members']))),
                left_on=['locus_tag'],
                right_on=['members'],
                how='left'
            )
            .drop('members', axis=1)
            .assign(contig = lambda df: df.contig_id.apply(lambda x: x.split('_')[0]))
            .assign(contig_order = lambda df: df.contig_id.apply(lambda x: int(x.split('_')[1])))
            .sort_values(['contig', 'contig_order', 'start', 'end'])
            .drop(['contig', 'contig_order'], axis=1)
            .merge(
                (pd.read_csv(casset_info, sep='\t', header=None)[1]
                    .str.split(':', expand=True)
                    .iloc[:,[1,4]]
                    .sort_values(4)
                    .rename(columns={ 1 :'casset_id', 4 :'locus_tag'})),
                on=['locus_tag'],
                how='left'
            )
            )

def eprint(*args):
    print(*args, file=sys.stderr)

def main(prokka_gff, seeds_with_close_dfs, casset_info, output_fname):
    eprint('Starting ...')
    eprint('Reading prokka annotation and dfs prediction ...')
    prokka_annotation_with_dfs_prediction_df = read_prokka(prokka_gff, seeds_with_close_dfs, casset_info)

    eprint('Finished table read.')
    eprint('Calculation metrics ...')
    (prokka_annotation_with_dfs_prediction_df
     .loc[:,['seeds']]
     .drop_duplicates()
     .reset_index(drop=True)
     .assign(members_count = lambda df: df.seeds.apply(lambda seed: prokka_annotation_with_dfs_prediction_df.query('seeds == @seed').shape[0]))
     .assign(members_names = lambda df: df.seeds.apply(
         lambda seed: Counter(prokka_annotation_with_dfs_prediction_df.query('seeds == @seed')['info'].apply(
             lambda info: dict(map(lambda x: x.split('='), info.split(';'))).get('Name', 'Hypothetical Protein')
         ).values.tolist())
     ))
     .assign(members_names = lambda df: get_percentage_of_member_count(df))
     .assign(icity = lambda df: df.seeds.apply(
         lambda seed: prokka_annotation_with_dfs_prediction_df.query('seeds == @seed').assign(icity = lambda df1: df1.casset_id.dropna().shape[0]/df1.casset_id.shape[0]).icity.values.tolist()[0]
     ))
     .assign(members_close_to_baits = lambda df: df.seeds.apply(
         lambda seed: prokka_annotation_with_dfs_prediction_df.query('seeds == @seed').casset_id.dropna().shape[0]
     ))
     .merge(
         get_bait_distances(prokka_annotation_with_dfs_prediction_df),
         on='seeds',
         how='left'
     )
     .assign(diversity_of_close_systems = lambda df: df.seeds.apply(lambda seed: Counter(prokka_annotation_with_dfs_prediction_df.query('seeds == @seed and casset_id.notna()')[['contig_id', 'casset_id']]
                                                                                         .merge(
                                                                                             prokka_annotation_with_dfs_prediction_df.query('casset_id.notna()')[['contig_id','casset_id', 'dfs_prediction']].dropna(),
                                                                                             on=['contig_id','casset_id'],
                                                                                             how='inner').dfs_prediction.values.tolist())))
     .assign(diversity_score = lambda df: df.diversity_of_close_systems.apply(lambda x: sum(x.values())))
     .assign(all_dfs = prokka_annotation_with_dfs_prediction_df.dfs_prediction.dropna().shape[0])
     .assign(diversity_score = lambda df: tuple(map(lambda tup: tup[0]/tup[1], zip(df.diversity_score.values.tolist(), df.all_dfs.values.tolist()))))
     .assign(diversity_of_close_systems = lambda df: df.diversity_of_close_systems.apply(lambda x: str(dict(x))).replace('{}', np.nan))
     .drop('all_dfs', axis=1)
     ).fillna('NA').to_csv(output_fname, index=False, sep='\t')

    eprint('All jobs finished !')

if __name__ == '__main__':
    main(*sys.argv[1:])
