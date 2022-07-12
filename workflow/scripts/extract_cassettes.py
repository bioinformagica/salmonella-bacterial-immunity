#!/usr/bin/env python3
"""
Extract cassettes from prokka gbk files with the
Padloc and DefesenFinder output files.

Usage:
    extract_cassettes.py ( --gbk_file=PATH ) ( --padloc_table=PATH ) ( --defensefinder_table=PATH )
                         ( --output_path=PATH ) [ --n_genes=INT ]

Options:
    --gbk_file=PATH             Path to the prokka gbk file.
    --padloc_table=PATH         Output padloc table.
    --defensefinder_table=PATH  Output defensefinder table.
    --output_path=PATH          Path to write the output files.
    --n_genes=INT               Number of genes to get from up and down stream
                                of the interest genes a n_gene = 10 will generate
                                cassettes with 21 genes [default: 10].
"""
import os
import sys
from pathlib import Path
from typing import List
from collections import namedtuple
import itertools
import logging
logger = logging.getLogger(os.path.basename(__file__).replace('.py', ''))
from docopt import docopt
import pandas as pd
from Bio import SeqIO

def read_padloc(file_path: str, separator: str, columns_to_keep: List[str], renamed_names: List[str]):
    return pd.read_csv(file_path, sep=separator) \
             .loc[:, columns_to_keep] \
             .rename(columns=dict(zip(columns_to_keep, renamed_names)))

def read_defensefinder(file_path: str, separator: str, columns_to_keep: List[str], renamed_names: List[str]):
    df = pd.read_csv(file_path, sep=separator)

    temp = df.iloc[:, [1,2,5,7]] \
             .assign(protein_in_syst=df['protein_in_syst'].str.split(','), name_of_profiles_in_sys=df['name_of_profiles_in_sys'].str.split(',')) \
             .explode('protein_in_syst') \
             .explode('name_of_profiles_in_sys')

    a = temp.iloc[:,[0,1,2]].explode('protein_in_syst').drop_duplicates().reset_index(drop=True)
    b = temp.iloc[:,[0,1,3]].explode('name_of_profiles_in_sys').drop_duplicates().reset_index(drop=True)

    return pd.concat([a, b], axis=1) \
             .iloc[:,[2,5,1]] \
             .reset_index(drop=True) \
             .rename(columns=dict(zip(columns_to_keep, renamed_names)))


def parser_gbk_file(gbk_file: str) -> list:
    contig = namedtuple('Contig', 'id features')
    # Parsing gbk file
    contigs = []
    with open(gbk_file) as f:
        for record in SeqIO.parse(f, 'genbank'):
            features = []
            for feature in filter(lambda fet: fet.type == 'CDS', record.features):
                features.append(feature)
            contigs.append(contig(record.id, features))
    return contigs

def get_cassettes_locus_tags(all_loci: str, interest_loci: str, n_genes: int) -> list:
    cassettes_locus_tags = []
    for interest_locus in interest_loci:
        interest_locus_index = all_loci.index(interest_locus)
        left_loci = all_loci[:interest_locus_index]
        if len(left_loci) > n_genes:
            left_loci = all_loci[ interest_locus_index - n_genes : interest_locus_index ]
        right_loci = all_loci[ interest_locus_index : interest_locus_index + n_genes + 1 ]
        cassettes_locus_tags.append(left_loci + right_loci)
    return cassettes_locus_tags

def merge_overlapping_intervals(intervals: tuple) -> tuple:
    stack = []

    for interval in sorted(intervals, key=lambda tup: tup[0]):

        if not stack:
            stack.append(interval)
            continue

        query = stack[-1]

        if interval[0] <= query[1]:
            new_end = max(query[1], interval[1])
            stack[-1] = (query[0], new_end)
            continue

        stack.append(interval)

    return stack


def merge_overlapping_cassettes(cassettes: tuple) -> tuple:
    get_start_end_from_cassette = lambda cassette: (cassette[0].location.start.position, cassette[-1].location.end.position)
    remove_duplicated_loci = lambda merged_cassettes: sorted(set(itertools.chain.from_iterable(merged_cassettes)), key = lambda feature: feature.location.start.position)

    cassettes_ranges = tuple(map(lambda cassette: get_start_end_from_cassette(cassette), cassettes))
    merged_intervals = merge_overlapping_intervals(cassettes_ranges)

    merged_cassettes_loci = remove_duplicated_loci(cassettes)

    # if there is only on interval we can just return all
    # loci after remove duplications
    if len(merged_intervals) == 1:
        return [merged_cassettes_loci]

    merged_cassettes = []
    for interval in merged_intervals:
        cassette = list(filter(lambda feature: feature.location.start.position in range(*interval), merged_cassettes_loci))
        merged_cassettes.append(cassette)

    return merged_cassettes

def extract_cassettes_from_gbk(gbk_file: str, defence_genes_locus_tags: list, n_genes: int, output_path: str):
    contigs = parser_gbk_file(gbk_file)
    strand = { -1 : '-', 1 : '+' }

    records = []
    for contig in contigs:
        contig_locus_tags = tuple(map(lambda fet: fet.qualifiers.get('locus_tag')[0], contig.features))
        interest_loci = tuple(filter(lambda locus_tag: locus_tag in defence_genes_locus_tags, contig_locus_tags))
        if not interest_loci:
            continue

        cassettes_locus_tags = get_cassettes_locus_tags(all_loci=contig_locus_tags, interest_loci=interest_loci, n_genes=n_genes)
        cassettes = tuple(map(lambda cassette: tuple(filter(lambda fet: fet.qualifiers.get('locus_tag')[0] in cassette, contig.features)), cassettes_locus_tags))
        cassettes = merge_overlapping_cassettes(cassettes)

        logger.info('Contig {} have {} cassette(s).'.format(contig.id, len(cassettes)))

        cassette_counter = 1
        for cassette in cassettes:
            for feature in cassette:
                translation = feature.qualifiers.get('translation', [None])[0]
                locus_tag = feature.qualifiers.get('locus_tag')[0]
                if translation:
                    records.append(
                        '>Contig_{}:Cassette_{}:GOI_{}:Strand_{}:{}\n{}\n'.format(
                            contig.id,
                            cassette_counter,
                            str(locus_tag in defence_genes_locus_tags),
                            strand[feature.strand],
                            locus_tag,
                            translation,
                        )
                    )
            cassette_counter += 1

    output_file_name = os.path.join(output_path, 'Cassettes.faa')
    logger.info('Writing cassettes to file {}.'.format(output_file_name))
    with open(output_file_name, 'w') as f:
        f.write(''.join(records))



def main(*args, **kwargs) -> None:
    # Logging setup
    logging.basicConfig(
        level=logging.DEBUG,
        datefmt="%Y-%m-%d %H:%M",
        format="[%(name)s][%(asctime)s][%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(),
        ]
    )

    logger.info('ARGS: {}'.format(kwargs))


    merged_format = ['locus_tags', 'names', 'systems']

    logger.info('Reading padloc table ...')
    # Reading feature tables
    padloc_df = read_padloc(
        file_path=kwargs['--padloc_table'],
        separator=',',
        columns_to_keep=['target.name', 'protein.name', 'system'],
        renamed_names=merged_format,
    ).assign(tool = 'padloc')

    logger.info('Reading defensefinder table ...')
    defensefinder_df = read_defensefinder(
        file_path=kwargs['--defensefinder_table'],
        separator='\t',
        columns_to_keep=['protein_in_syst', 'name_of_profiles_in_sys', 'subtype'],
        renamed_names=merged_format,
    ).assign(tool = 'defensefinder')

    logger.info('Merging tables ...')
    # Joining results
    merged_features_df = pd.concat([padloc_df, defensefinder_df]).reset_index(drop=True)

    logger.info('Writing merge csv ...')
    merged_features_df.to_csv(os.path.join(kwargs['--output_path'], 'merged_defense_systems_prediction.csv'), index=False, sep=',')

    # Read gbk file
    logger.info('Starting reading the gbk file ...')
    extract_cassettes_from_gbk(
        gbk_file=kwargs['--gbk_file'],
        defence_genes_locus_tags=sorted(merged_features_df.locus_tags.unique().tolist(), key=lambda x: int(x.split('_')[1])), # this must be sorted
        n_genes=int(kwargs['--n_genes']),
        output_path=kwargs['--output_path'],
    )

    logger.info('Finished all jobs !')






if __name__ == '__main__':
    main(**docopt(__doc__))
