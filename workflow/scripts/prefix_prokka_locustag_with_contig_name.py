#!/usr/bin/env python3

import sys
from Bio import SeqIO

def usage() -> None:
    help_massage = """
Description:

    This script reads a prokka gbk file and prints a faa file with
    locustags prefixed by the contig name to be used by defencefinder.

Usage:

    python3 {} <gbk_file_path>

    """.format(__file__)
    print(help_massage)

def main(gbk_prokka_file: str) -> None:
    with open(gbk_prokka_file, 'r') as f:
        for record in SeqIO.parse(f, 'genbank'):
            for feature in record.features:
                if feature.type == 'CDS':
                    try:
                        print('>{}_{}\n{}'.format(
                            record.id,
                            feature.qualifiers['locus_tag'][0],
                            feature.qualifiers['translation'][0],
                        ))
                    except KeyError as e:
                        pass

if __name__ == '__main__':

    if not len(sys.argv[1:]) == 1:
        usage()
        sys.exit(1)
    try:
        main(*sys.argv[1:])
    except BrokenPipeError:
        sys.exit(0)
