# [[file:../../main.org::*Get Gene Neighborhood][Get Gene Neighborhood:1]]
rule get_ppanggolin_input_list:
    input:
        gffs_files = input_table.gff_file_paths,
    output:
        organisms_fasta_list = temp(temp_dir_path / 'input.ppanggolin.tsv'),
    run:
        # Write gffs ids and paths to a file to use as input in ppanggolin.
        input_table[['genome_ids', 'gff_file_paths']].to_csv(
            output.organisms_fasta_list,
            sep='\t',
            index=False,
            header=None
        )

rule run_ppanggolin:
    input:
        get_ppanggolin_input_list.output.organisms_fasta_list,
    output:
        pangenome_dir = directory(outputs_dir_path / 'ppanggolin'),
    conda:
        ENVS_DIR_PATH / 'ppanggolin_env.yaml',
    threads:
        get_cores_perc(1)
    shell:
        'ppanggolin all'
        ' --anno {input}'
        ' --output {output.pangenome_dir}'
        ' --tmpdir {TEMP_DIR_PATH}'
        ' --cpu {threads}'
# Get Gene Neighborhood:1 ends here
