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

def extract_cassettes_from_gbk(gbk_file: str, defence_genes_locus_tags: set, n_genes: int, output_path: str):
    contig = namedtuple('Contig', 'id features defence_genes')
    protein = namedtuple('Protein', 'locus_tag feature')
    cassette = namedtuple('Cassette', 'locus_tags start end')
    contigs_with_defence_genes = {}

    is_overlapping = lambda a, b: b.start > a.start and b.start < a.end
    sort_locus_tags = lambda locus_tag: int(locus_tag.split('_')[1])
    fasta_records_to_write = []
    output_file_name = os.path.join(output_path, 'Cassettes.faa')


    logger.info('Retriving locus_tags from each contig ...')
    for record in SeqIO.parse(gbk_file, 'genbank'):
        contig_id = record.id
        contig_features = []
        for feature in record.features:
            if feature.type == 'CDS':
                try:
                    contig_features.append(protein(feature.qualifiers['locus_tag'][0], feature))
                except KeyError as e:
                    pass

        # check if contig has a defense system protein
        check_locustags = tuple(map(lambda protein: protein.locus_tag, contig_features))
        defence_genes = list(map(lambda x: x in defence_genes_locus_tags, check_locustags))
        defence_genes = [ locus_tag for locus_tag in check_locustags if locus_tag in defence_genes_locus_tags ]

        if not defence_genes:
            continue

        defence_genes = sorted(defence_genes, key=sort_locus_tags)

        contigs_with_defence_genes[contig_id] = contig(contig_id, contig_features, defence_genes)


    logger.info('There are {} contigs with defence systems.'.format(len(contigs_with_defence_genes)))
    for contig in contigs_with_defence_genes.values():
        cassettes = []
        all_locus_tags = tuple(map(lambda protein: protein.locus_tag, contig.features))
        for defence_gene in contig.defence_genes:
            gene_index = all_locus_tags.index(defence_gene)
            if gene_index < n_genes:
                cassette_locus_tags = all_locus_tags[:gene_index + n_genes + 1 ]
            else:
                cassette_locus_tags = all_locus_tags[gene_index - n_genes : gene_index + n_genes + 1 ]

            first_locus_tag, last_locus_tag = cassette_locus_tags[0], cassette_locus_tags[-1]
            first_locus_tag_position, last_locus_tag_position = None, None

            proteins = iter(contig.features)
            protein = next(proteins, None)

            while not all((first_locus_tag_position, last_locus_tag_position)) or protein:
                if protein.locus_tag == first_locus_tag:
                    first_locus_tag_position = protein.feature.location.start.position

                if protein.locus_tag == last_locus_tag:
                    last_locus_tag_position = protein.feature.location.end.position

                protein = next(proteins, None)

            cassettes.append(cassette(cassette_locus_tags, first_locus_tag_position, last_locus_tag_position))

        if len(cassettes) > 1:
            # Merging overlapping cassettes
            cassettes.sort(key = lambda cassette: cassette.start)
            cassettes_overlaps_merge = []
            cassettes_overlaps_merge.append(cassettes[0])

            for i in range(1, len(cassettes)):
                pop_element = cassettes_overlaps_merge.pop()
                if is_overlapping(pop_element, cassettes[i]):
                    new_locus_tags = sorted(list(set([*pop_element.locus_tags, *cassettes[i].locus_tags])), key=sort_locus_tags)
                    new_start, new_end = pop_element.start, max(pop_element.end, cassettes[i].end)
                    new_element = cassette(new_locus_tags, new_start, new_end)
                    cassettes_overlaps_merge.append(new_element)
                else:
                    cassettes_overlaps_merge.append(pop_element)
                    cassettes_overlaps_merge.append(cassettes[i])
        else:
            cassettes_overlaps_merge = cassettes

        logger.info('Contig {} has {} cassette(s).'.format(contig.id, len(cassettes_overlaps_merge)))

        # Write cassettes
        cassette_counter = 1
        for i in cassettes_overlaps_merge:
            for protein in contig.features:
                if protein.locus_tag in i.locus_tags:
                    try:
                        fasta_records_to_write.append('>Contig_{}:Cassette_{}:{}\n{}\n'.format(
                            contig.id,
                            cassette_counter,
                            protein.locus_tag,
                            protein.feature.qualifiers['translation'][0]
                        ))
                    except KeyError as e:
                        pass
            cassette_counter += 1

    logger.info('Writing cassettes to file {}.'.format(output_file_name))
    with open(output_file_name, 'w') as f:
        f.write(''.join(fasta_records_to_write))




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
        defence_genes_locus_tags=set(merged_features_df.locus_tags),
        n_genes=int(kwargs['--n_genes']),
        output_path=kwargs['--output_path'],
    )

    logger.info('Finished all jobs !')






if __name__ == '__main__':
    main(**docopt(__doc__))
