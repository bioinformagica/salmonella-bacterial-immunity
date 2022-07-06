#!/usr/bin/env python3
"""
Extract cassettes from prokka gbk files with the
Padloc and DefesenFinder output files.

Usage:
    extract_cassettes.py ( --gbk_file=PATH ) ( --padloc_table=PATH ) ( --defensefinder_table=PATH )
                         ( --output_path=PATH ) [ --prefix=STR ]

Options:
    --gbk_file=PATH             Path to the prokka gbk file.
    --padloc_table=PATH         Output padloc table.
    --defensefinder_table=PATH  Output defensefinder table.
    --output_path=PATH          Path to write the output files.
    --prefix=STR String         to prefix output files.

"""
import os
import sys
from pathlib import Path
from typing import List
from collections import defaultdict
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

    a = temp.iloc[:,[0,1,2]].explode('protein_in_syst').drop_duplicates()
    b = temp.iloc[:,[0,1,3]].explode('name_of_profiles_in_sys').drop_duplicates()

    return pd.concat([a, b], axis=1) \
             .iloc[:,[2,5,1]] \
             .reset_index(drop=True) \
             .rename(columns=dict(zip(columns_to_keep, renamed_names)))


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


    merged_format = ['locus_tag', 'name', 'system']

    # Reading feature tables
    padloc_df = read_padloc(
        file_path=kwargs['--padloc_table'],
        separator=',',
        columns_to_keep=['target.name', 'protein.name', 'system'],
        renamed_names=merged_format,
    ).assign(tool = 'padloc')

    defensefinder_df = read_defensefinder(
        file_path=kwargs['--defensefinder_table'],
        separator='\t',
        columns_to_keep=['protein_in_syst', 'name_of_profiles_in_sys', 'subtype'],
        renamed_names=merged_format,
    ).assign(tool = 'defensefinder')

    # Joining results
    merged_features_df = pd.concat([padloc_df, defensefinder_df]).reset_index(drop=True)

    # Read gbk file


if __name__ == '__main__':
    main(**docopt(__doc__))
