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
import re
import subprocess
import itertools
import logging
logger = logging.getLogger(Path(__file__).stem)

# External modules
from docopt import docopt
import pandas as pd

def prodigal_csv_to_df(file_path: str):
    check_prodigal_record = lambda record: re.match('^[\d\w]+\.[\d\w]+\_\d+\_\d+\;\d+\;\d+\;\-?\d\;ID\=\d+\_\d+$', record)
    with open(file_path) as f:
        columns = ('locus_tag', 'start', 'end', 'strand', 'id')
        return pd.DataFrame(data=map(lambda x: x.replace('ID=', '').split(';'), filter(check_prodigal_record, f.read().splitlines())), columns=columns) \
                .assign(genome_id = lambda df: df.locus_tag.str.split('_', expand=True).iloc[:, 0] )

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

    defensefinder_map = get_defensefinder_df(run_gnufind(kwargs['--defensefinder_dir'], '*defense_finder_systems.tsv')) \
        .reset_index(drop=True) \
        .rename(columns=dict(zip(['protein_in_syst', 'name_of_profiles_in_sys', 'type'], ['locus_tag', 'name', 'system']))) \
        .groupby('locus_tag', as_index=False) \
        .apply(lambda df: df.head(1)) \
        .assign(annotation = lambda df: tuple(map(lambda x: ','.join(tuple(map('='.join, zip(*x)))), zip(itertools.repeat(('name', 'system')), df.loc[:,['name', 'system']].values.tolist()))) ) \
        .assign(annotation = lambda df: 'defensefinder:' + df['annotation']) \
        .drop(['name', 'system'], axis=1) \
        .set_index('locus_tag', drop=True).annotation.to_dict()


    logger.info('Reading padloc tables ...')

    padloc_map = get_padloc_df(run_gnufind(kwargs['--padloc_dir'], '*.fasta_padloc.csv')).loc[:,['target.name', 'protein.name','system']] \
        .rename(columns=dict(zip(['target.name', 'protein.name','system'], ['locus_tag', 'name', 'system']))) \
        .groupby('locus_tag', as_index=False) \
        .apply(lambda df: df.head(1)) \
        .assign(annotation = lambda df: tuple(map(lambda x: ','.join(tuple(map('='.join, zip(*x)))), zip(itertools.repeat(('name', 'system')), df.loc[:,['name', 'system']].values.tolist()))) ) \
        .assign(annotation = lambda df: 'padloc:' + df['annotation']) \
        .drop(['name', 'system'], axis=1) \
        .set_index('locus_tag', drop=True).annotation.to_dict()

    logger.info('Reading prodigal headers ...')
    prodigal_csv_to_df(kwargs['--prodigal_headers_csv']) \
        .assign(prediction = lambda df: df.locus_tag.apply(lambda locus_tag: ';'.join(map(lambda prediction_map: prediction_map.get(locus_tag, 'NA'), (padloc_map, defensefinder_map))))) \
        .query('prediction != "NA;NA"').to_csv(kwargs['--output'], sep=';', index=False)

    logger.info('All finished !')

if __name__ == '__main__':
    main(**docopt(__doc__))
