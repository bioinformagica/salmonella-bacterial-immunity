#!/usr/bin/env python3
"""
Extract cassettes from prokka gbk files with the
Padloc and DefesenFinder output files.

Usage:
    extract_cassettes.py ( --prokka_dir=PATH ) ( --prodigal_gois=PATH )
                         [ --logfile=PATH ] [ --n_genes=INT ]
                         [ --sample_input=INT ] [ --threads=INT ]

Options:
    --prokka_dir=PATH           Path to the prokka gbk files.
    --prodigal_gois=PATH        Table with prodigal GOIs.
    --logfile=PATH              Path to store log [default: ./extract_cassettes.log].
    --n_genes=INT               Number of genes to get from up and down stream
                                of the interest genes a n_gene = 10 will generate
                                cassettes with 21 genes [default: 10].
    --sample_input=INT          Only run max N genome cassette extration.
    --threads=INT               Numbers of threads to use [default: 1].
"""
import os
import sys
from pathlib import Path
from typing import List
from collections import namedtuple
import itertools
import subprocess
from multiprocessing import Pool
from functools import partial
import random
random.seed(42)
import gzip
import logging
logger = logging.getLogger(os.path.basename(__file__).replace('.py', ''))
from docopt import docopt
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

def run_gnufind(query_dir, match_string, ftype='f'):
    cmd = [
        'find',
        query_dir,
        '-type', ftype,
        '-name', match_string
    ]
    found_files = subprocess.run(cmd, check=True, stdout=subprocess.PIPE).stdout.decode('UTF-8').splitlines()
    assert found_files, 'No files found ! Dir {} with prefix {}.'.format(query_dir, match_string)
    logger.info('Found {} in dir {} with pattern {}.'.format(len(found_files), query_dir, match_string))
    return found_files


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

def check_cds_overlap(feature, coord):
    feature_start = feature.location.start.position
    feature_end = feature.location.end.position
    feature_strand = feature.location.strand
    return max(feature_start, coord[0]) < min(feature_end, coord[1]) and feature_strand == coord[2]


def extract_cassettes_from_gbk(gbk_file: str, prodigal_gois_df, n_genes: int, no_clobber=True):
    logger.info('Starting extraction from file: {}.'.format(gbk_file))

    output_file_name = gbk_file.replace('.gbk', '.cassettes.faa.gz')
    if Path(output_file_name).exists():
        logger.info('File already exists, skiping: {}.'.format(output_file_name))
        return


    goi = namedtuple('Goi', 'prediction feature')
    wanted_contigs = prodigal_gois_df.contig_id.unique().tolist()
    contigs = dict(map(lambda contig: (contig.id, contig), tuple(filter(lambda contig: contig.id in wanted_contigs, parser_gbk_file(gbk_file)))))
    strand = { -1 : '-', 1 : '+' }
    records = []

    # Merging PRODIGAL and prokka annotation
    for contig, df in sorted(prodigal_gois_df.groupby('contig_id'), key=lambda tup: int(tup[0].split('_')[-1])):
        gois = []
        prodigal_GOI_cds_coords = df.loc[:,['start', 'end', 'strand', 'prediction']].values.tolist()
        for feature in contigs[contig].features:
            if not prodigal_GOI_cds_coords:
                break
            for idx, coord in enumerate(prodigal_GOI_cds_coords):
                if check_cds_overlap(feature, coord):
                    gois.append(goi(coord[-1], feature))
                    _ = prodigal_GOI_cds_coords.pop(idx)
                    break

        all_loci = tuple(map(lambda feature: feature.qualifiers.get('locus_tag')[0], contigs[contig].features))
        goi_loci = tuple(map(lambda goi: goi.feature.qualifiers.get('locus_tag')[0], gois))
        cassettes_locus_tags = get_cassettes_locus_tags(all_loci=all_loci, interest_loci=goi_loci, n_genes=n_genes)


        cassettes = tuple(map(lambda cassette: tuple(filter(lambda fet: fet.qualifiers.get('locus_tag')[0] in cassette, contigs[contig].features)), cassettes_locus_tags))
        cassettes = merge_overlapping_cassettes(cassettes)

        logger.info('Contig {} have {} cassette(s).'.format(contigs[contig].id, len(cassettes)))

        cassette_counter = 1
        for cassette in cassettes:
            for feature in cassette:
                translation = feature.qualifiers.get('translation', [None])[0]
                locus_tag = feature.qualifiers.get('locus_tag')[0]
                if translation:
                    records.append(
                        '>Contig_{}:Cassette_{}:GOI_{}:Strand_{}:{}\n{}\n'.format(
                            contigs[contig].id,
                            cassette_counter,
                            str(locus_tag in goi_loci),
                            strand[feature.strand],
                            locus_tag,
                            translation,
                        )
                    )
            cassette_counter += 1

    with gzip.open(output_file_name, 'wb') as f:
        f.write(''.join(records).encode())

    logger.info('Finished extraction: {}.'.format(output_file_name))



def main(*args, **kwargs) -> None:
    # Logging setup
    logging.basicConfig(
        level=logging.DEBUG,
        datefmt="%Y-%m-%d %H:%M",
        format="[%(name)s][%(asctime)s][%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(kwargs['--logfile']),
            # logging.StreamHandler(),
        ]
    )

    logger.info('ARGS: {}'.format(kwargs))

    logger.info('Listing prokka gbk files ...')
    gbk_file_cache = Path(kwargs['--prokka_dir']).absolute().parent / 'prokka_dir_cache.txt'
    if gbk_file_cache.exists():
        logger.info('prokka_dir cache exists, reading it ...')
        with open(gbk_file_cache) as f:
            gbk_files_list = f.read().splitlines()
    else:
        logger.info('prokka_dir cache do not exists, creating it ...')
        gbk_files_list = run_gnufind(kwargs['--prokka_dir'], '*_out.gbk' )
        with open(gbk_file_cache, 'w') as o:
            o.write('\n'.join(gbk_files_list))

    gbk_files: dict = dict(tuple(map(
        lambda file_obj: (file_obj.parent.name.replace('_out', ''), str(file_obj)),
        map(lambda file_path: Path(file_path), gbk_files_list)
    )))

    logger.info('Reading prodigal GOIs ...')
    prodigal_GOIs_df = (pd.read_csv(kwargs['--prodigal_gois'], sep=';', dtype={'genome_id': str})
                         .groupby('genome_id').filter(lambda df: bool(gbk_files.get(df.name)))
                         .assign(contig_id = lambda df: df.locus_tag.str.split('_').apply(lambda x: '_'.join(x[:2])) ))


    extract_cassettes_from_gbk_part = partial(extract_cassettes_from_gbk, n_genes=int(kwargs['--n_genes']))

    jobs = tuple(map(lambda tup: (gbk_files[tup[0]], tup[1]), prodigal_GOIs_df.groupby('genome_id')))

    if kwargs['--sample_input'] and len(jobs) > int(kwargs['--sample_input']):
        jobs = random.sample(jobs, int(kwargs['--sample_input']))

    logger.info('Starting {} parallel jobs with {} threads'.format(len(jobs), kwargs['--threads']))
    with Pool(processes=int(kwargs['--threads'])) as p:
        _ = p.starmap(extract_cassettes_from_gbk_part, jobs)

    logger.info('Finished all jobs !')

if __name__ == '__main__':
    main(**docopt(__doc__))
