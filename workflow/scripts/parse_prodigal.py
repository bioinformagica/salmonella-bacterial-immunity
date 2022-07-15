#!/usr/bin/env python
"""
Mag GOIs to prodigal headers

Usage:
    parse_prodigal.py ( --prodigal_headers_csv=PATH ) ( --padloc_dir=PATH )
                      ( --defensefinder_dir=PATH ) ( --output=PATH )

Options:
    --prodigal_headers_csv=PATH File with prodigal headers.
    --padloc_dir=PATH           Path to padloc predictions.
    --defensefinder_dir=PATH    Path to defesenfinder predictions.
    --output=PATH               Name of the output_table.
"""

# Native modules
from pathlib import Path
import sys
import re
import subprocess
import itertools
import logging
logger = logging.getLogger(Path(__file__).stem)

# External modules
from docopt import docopt
import pandas as pd

def prodigal_csv_to_df(records: tuple):
    columns = ('locus_tag', 'start', 'end', 'strand', 'id')
    return (pd.DataFrame(data=records, columns=columns)
             .assign(id = lambda df: df.id.str.replace('ID=', ''))
             .assign(genome_id = lambda df: df.locus_tag.str.split('_', expand=True).iloc[:, 0] ))

def run_gnufind(query_dir, match_string, ftype='f'):
    cmd = [
        'find', '-L',
        query_dir,
        '-type', ftype,
        '-name', match_string
    ]
    found_files = subprocess.run(cmd, check=True, stdout=subprocess.PIPE).stdout.decode('UTF-8').splitlines()
    assert found_files, 'No files found ! Dir {} with prefix {}.'.format(query_dir, match_string)
    logger.info('Found {} in dir {} with pattern {}.'.format(len(found_files), query_dir, match_string))
    return tuple(map(lambda x: str(Path(x).absolute()), found_files))

def get_defensefinder_df(defensefinder_output_files):
    return pd.concat(map(lambda tsv: pd.read_csv(tsv, sep='\t'), defensefinder_output_files)) \
            .loc[:,['protein_in_syst', 'name_of_profiles_in_sys', 'type']] \
            .assign(protein_in_syst = lambda df: df.protein_in_syst.str.split(',')) \
            .assign(name_of_profiles_in_sys = lambda df: df.name_of_profiles_in_sys.str.split(',')) \
            .explode(['protein_in_syst', 'name_of_profiles_in_sys']).reset_index(drop=True)

def get_padloc_df(padloc_output_files):
    return pd.concat(map(lambda csv: pd.read_csv(csv, sep=','), padloc_output_files))

def main(**kwargs):

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

    logger.info('Reading defensefinder tables ...')

    defensefinder_files_cache = Path(kwargs['--output']).absolute().parent / 'defensefinder_tables_paths.txt'
    padloc_files_cache = Path(kwargs['--output']).absolute().parent / 'padloc_tables_paths.txt'

    if defensefinder_files_cache.exists():
        logger.info('Defensefinder cache exists reading it ...')
        with open(defensefinder_files_cache) as f:
            defensefinder_files = f.read().splitlines()
    else:
        logger.info('Defensefinder cache  do not exists creating ...')
        defensefinder_files = run_gnufind(kwargs['--defensefinder_dir'], '*defense_finder_systems.tsv')
        with open(defensefinder_files_cache, 'w') as f:
            f.write('\n'.join(defensefinder_files))

    if padloc_files_cache.exists():
        logger.info('Defensefinder cache exists reading it ...')
        with open(padloc_files_cache) as f:
            padloc_files = f.read().splitlines()
    else:
        logger.info('Defensefinder cache  do not exists creating ...')
        padloc_files = run_gnufind(kwargs['--padloc_dir'], '*.fasta_padloc.csv')
        with open(padloc_files_cache, 'w') as f:
            f.write('\n'.join(padloc_files))

    defensefinder_map = get_defensefinder_df(defensefinder_files) \
        .reset_index(drop=True) \
        .rename(columns=dict(zip(['protein_in_syst', 'name_of_profiles_in_sys', 'type'], ['locus_tag', 'name', 'system']))) \
        .groupby('locus_tag', as_index=False) \
        .apply(lambda df: df.head(1)) \
        .assign(annotation = lambda df: tuple(map(lambda x: ','.join(tuple(map('='.join, zip(*x)))), zip(itertools.repeat(('name', 'system')), df.loc[:,['name', 'system']].values.tolist()))) ) \
        .assign(annotation = lambda df: 'defensefinder:' + df['annotation']) \
        .drop(['name', 'system'], axis=1) \
        .set_index('locus_tag', drop=True).annotation.to_dict()


    logger.info('Reading padloc tables ...')

    padloc_map = get_padloc_df(padloc_files).loc[:,['target.name', 'protein.name','system']] \
        .rename(columns=dict(zip(['target.name', 'protein.name','system'], ['locus_tag', 'name', 'system']))) \
        .groupby('locus_tag', as_index=False) \
        .apply(lambda df: df.head(1)) \
        .assign(annotation = lambda df: tuple(map(lambda x: ','.join(tuple(map('='.join, zip(*x)))), zip(itertools.repeat(('name', 'system')), df.loc[:,['name', 'system']].values.tolist()))) ) \
        .assign(annotation = lambda df: 'padloc:' + df['annotation']) \
        .drop(['name', 'system'], axis=1) \
        .set_index('locus_tag', drop=True).annotation.to_dict()

    logger.info('Reading prodigal headers ...')

    records_cache = Path(kwargs['--output']).absolute().parent / 'records_prodigal_cache.txt'

    if records_cache.exists():
        logger.info('Records_Prodigal cache exists reading it ...')
        with open(records_cache) as f:
            records = tuple(map(lambda x: tuple(x.split(';')), f.read().splitlines()))
    else:
        logger.info('Records_Prodigal cache  do not exists creating ...')
        records = []
        check_prodigal_record = lambda record: re.match('^[\d\w]+\.[\d\w]+\_\d+\_\d+\;\d+\;\d+\;\-?\d\;ID\=\d+\_\d+$', record)
        with open(kwargs['--prodigal_headers_csv']) as f:
            for line in f:
                splitted = line.strip().split(';')
                if any(map(lambda loci_dir: loci_dir.get(splitted[0]), [defensefinder_map, padloc_map])) and check_prodigal_record(line):
                    records.append(tuple(splitted))
        with open(records_cache, 'w') as o:
            o.write('\n'.join(tuple(map(lambda x: ';'.join(x), records))))

    (prodigal_csv_to_df(records)
     .assign(prediction = lambda df: df.locus_tag.apply(lambda locus_tag: ';'.join(map(lambda prediction_map: prediction_map.get(locus_tag, 'NA'), (padloc_map, defensefinder_map)))))
     .query('prediction != "NA;NA"').to_csv(kwargs['--output'], sep=';', index=False))

    logger.info('All finished !')

if __name__ == '__main__':
    main(**docopt(__doc__))
