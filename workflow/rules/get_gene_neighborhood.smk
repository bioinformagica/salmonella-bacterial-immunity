# [[file:../../main.org::*Find defense genes with padloc][Find defense genes with padloc:1]]
# rule padloc_search_defense_genes:
#     input:
#         gff_file = join_path(config['prokka_dir'], '{genome_id}_out', '{genome_id}_out.gff'),
#         faa_file = join_path(config['prokka_dir'], '{genome_id}_out', '{genome_id}_out.faa')
#     output:
#         padloc_out_dir = directory(join_path(results_dir, 'padloc', '{genome_id}')),
#         padloc_csv = join_path(results_dir, 'padloc', '{genome_id}', '{genome_id}_out_padloc.csv'),
#     params:
#         gff_nofasta_file = join_path(config['prokka_dir'], '{genome_id}_out', '{genome_id}_out_nofasta.gff'),
#     threads:
#         get_cores_perc(0.1)
#     conda:
#         '../envs/padloc_env.yaml'
#     shell:
#         "sed '/^##FASTA/Q' {input.gff_file} > {params.gff_nofasta_file} && "
#         'mkdir -p {output.padloc_out_dir} && '
#         'padloc --faa {input.faa_file} --gff {params.gff_nofasta_file} --outdir {output.padloc_out_dir} && '
#         'rm -v {params.gff_nofasta_file}'
# Find defense genes with padloc:1 ends here

# [[file:../../main.org::*Find defense genes with defense finder][Find defense genes with defense finder:1]]
# rule defensefinder_search_defense_genes:
#     input:
#         faa_file = join_path(config['prokka_dir'], '{genome_id}_out', '{genome_id}_out.faa'),
#     output:
#         defensefinder_out_dir = directory(join_path(results_dir, 'defensefinder', '{genome_id}')),
#         defensefinder_tsv = join_path(results_dir, 'defensefinder', '{genome_id}', 'defense_finder_systems.tsv'),
#     params:
#         ,**config['params']['defensefinder'],
#     threads:
#         get_cores_perc(0.1)
#     conda:
#         '../envs/defensefinder_env.yaml'
#     shell:
#         'defense-finder run '
#         '--db-type {params.db_type} '
#         '--out-dir {output.defensefinder_out_dir} '
#         '--workers {threads} '
#         '{input.faa_file} '
# Find defense genes with defense finder:1 ends here

# [[file:../../main.org::*Prodigal: find GOIs][Prodigal: find GOIs:1]]
rule add_GOIs_to_prodigal_headers:
    input:
        prodigal_headers = join_path(results_dir, 'prodigal_annotation.csv'),
        padloc_dir = join_path(results_dir, 'Padloc_csv'),
        defensefinder_dir = join_path(results_dir, 'DefenseFinder'),
        parser = join_path(scripts_dir, 'parse_prodigal.py'),
    output:
        prodigal_with_GOIs = join_path(results_dir, 'prodigal_with_GOIs.csv')
    conda:
        '../envs/misc_env.yaml'
    shell:
        'python3 {input.parser}'
        ' --prodigal_headers_csv {input.prodigal_headers}'
        ' --padloc_dir {input.padloc_dir}'
        ' --defensefinder_dir {input.defensefinder_dir}'
        ' --output {output.prodigal_with_GOIs}'
# Prodigal: find GOIs:1 ends here

# [[file:../../main.org::*Get gene cassettes][Get gene cassettes:1]]
rule get_gene_cassettes:
    input:
        prokka_dir = join_path(results_dir, 'Prokka'),
        prodigal_with_GOIs = join_path(results_dir, 'prodigal_with_GOIs.csv'),
        script = join_path(scripts_dir, 'extract_cassettes.py'),
    output:
        finish_signal = join_path(results_dir, 'finish_signal_cassette_extraction.txt')
    params:
        **config['params']['get_gene_cassettes'],
        log_dir = join_path(snakefile_path, '..', 'logs'),
    conda:
        '../envs/misc_env.yaml'
    threads:
        get_cores_perc(1)
    shell:
       'python3 /home/hugo/projects/salmonella-bacterial-immunity/workflow/scripts/extract_cassettes.py'
       ' --prokka_dir results/Prokka/'
       ' --prodigal_gois results/prodigal_with_GOIs.csv'
       ' --threads {threads}'
       ' --logfile  $( echo {params.log_dir}/{rule}_$(date +%Y_%m_%d_-_%H_%M_%S).log )'
       ' --sample_input 100 &&'
       ' touch {output.finish_signal}'
# Get gene cassettes:1 ends here

# [[file:../../main.org::*Combine faa files][Combine faa files:1]]
rule merge_cassette_proteins:
    input:
        cassettes_faa = expand(join_path(results_dir, 'cassettes', 'extracted', '{genome_id}', 'Cassettes.faa.gz'), genome_id=GENOME_IDS),
    output:
        concatenated_faa_file = join_path(results_dir, 'cassettes', 'clusters', 'merged_cassette_proteins.faa.gz')
    threads:
        1
    shell:
        'zcat {input.cassettes_faa} > {output.concatenated_faa_file}'
# Combine faa files:1 ends here

# [[file:../../main.org::*Cluster cassette genes][Cluster cassette genes:1]]
rule mmseqs2_cluster_proteins:
    input:
        concatenated_faa_file = join_path(results_dir, 'cassettes', 'clusters', 'merged_cassette_proteins.faa.gz')
    output:
        mmseqs2_rep_gz = join_path(results_dir, 'cassettes', 'clusters', 'mmseqs2_rep_seq.fasta.gz'),
        mmseqs2_tsv_gz = join_path(results_dir, 'cassettes', 'clusters', 'mmseqs2_cluster.tsv.gz'),
    params:
        **config['params']['mmseqs2'],
        mmseqs2_prefix = join_path(results_dir, 'cassettes', 'clusters', 'mmseqs2'),
        mmseqs2_fasta = join_path(results_dir, 'cassettes', 'clusters', 'mmseqs2_all_seqs.fasta'),
        mmseqs2_tmp = join_path(results_dir, 'cassettes', 'clusters', 'tmp'),
        mmseqs2_rep = join_path(results_dir, 'cassettes', 'clusters', 'mmseqs2_rep_seq.fasta'),
        mmseqs2_tsv = join_path(results_dir, 'cassettes', 'clusters', 'mmseqs2_cluster.tsv'),
    conda:
        '../envs/mmseqs2_env.yaml'
    threads:
        get_cores_perc(1)
    shell:
        'mmseqs easy-cluster {input.concatenated_faa_file} {params.mmseqs2_prefix} '
        '{params.mmseqs2_tmp} --threads {threads} '
        '-c {params.min_cov} --min-seq-id {params.min_ide} && '
        'gzip -v {params.mmseqs2_rep} {params.mmseqs2_tsv} && '
        'rm -rfv {params.mmseqs2_fasta} {params.mmseqs2_tmp}'
# Cluster cassette genes:1 ends here

# [[file:../../main.org::*Create network dataframes][Create network dataframes:1]]
# rule create_network_dataframes:
#     input:
#         concatenated_faa_file = join_path(results_dir, 'cassettes', 'clusters', 'merged_cassette_proteins.faa.gz'),
#         mmseqs2_tsv = join_path(results_dir, 'cassettes', 'clusters', 'mmseqs2_cluster.tsv'),
#         get_network_scrip =
#     output:
#     params:
#     conda:
#         '../envs/mmseqs2_env.yaml'
#     threads:
#         get_cores_perc(1)
#     shell:
# Create network dataframes:1 ends here
