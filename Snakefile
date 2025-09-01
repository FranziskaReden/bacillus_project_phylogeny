# Snakemake workflow for Bacillus project
# Copyright (c) 2025 Franziska Reden
# This workflow is licensed under the GNU General Public License v3.0 (GPL-3.0)
# See LICENSE file for details.

configfile: "config.yaml"

rule all:
    input:
        expand("bacillus_assemblies/prokka/{genome}/{genome}.gff", genome=config["genome"]),
        expand("alignments/faa/muscle/{marker}.faa.afa", marker=config["marker_genes"])

rule annotate:
    input: "bacillus_assemblies/{genome}/{genome}.fasta"
    output: "bacillus_assemblies/prokka/{genome}/{genome}.gff"
    conda: "environments/prokka_envs.yaml"
    shell:"""
        prokka --outdir bacillus_assemblies/prokka/{wildcards.genome} --prefix {wildcards.genome} {input} --force"
    """

rule get_gtdb:
    output: 
        "gtdb/bac120_marker_genes_all.tar.gz",
        "gtdb/bac120_marker_genes_reps.tar.gz",
        "bacillus_assemblies/prokka/hmm_PGAP.tsv",
        "gtdb/bac120_taxonomy.tsv",
        "gtdb/bac120.tree"
    conda: "environments/prokka_envs.yaml"
    script: "scripts/retrieve_gtdb_data.py"

rule retrieve_marker_genes:
    input: 
        "gtdb/bac120_marker_genes_all.tar.gz",
        "gtdb/bac120_marker_genes_reps.tar.gz",
        "bacillus_assemblies/prokka/hmm_PGAP.tsv",
        "gtdb/bac120_taxonomy.tsv"
    output: "gtdb/bac120_taxonomy_expanded.tsv"
    conda: "environments/prokka_envs.yaml"
    script: "scripts/retrieve_marker_genes.py"

rule compare_marker_genes:
    input: "gtdb/bac120_taxonomy_expanded.tsv"
    output: "bacillus_assemblies/prokka/{genome}/{genome}_closest_genomes.tsv",
        "bacillus_assemblies/prokka/{genome}/{genome}_genes_identity.tsv"
    conda: "environments/prokka_envs.yaml"
    shell: """
        python scripts/compare_marker_genes.py -a {wildcards.genome}
    """

rule choose_genomes:
    input: 
        expand("bacillus_assemblies/prokka/{genome}/{genome}_closest_genomes.tsv", genome=config["genome"]),
        expand("bacillus_assemblies/prokka/{genome}/{genome}_genes_identity.tsv", genome=config["genome"])
    output: 
        "marker_genes/best_hits.tsv",
        expand("alignments/faa/{marker}.faa", marker=config["marker_genes"]),
        expand("alignments/fna/{marker}.fna", marker=config["marker_genes"])
    conda: "environments/prokka_envs.yaml"
    script: "scripts/choose_genomes.py"

rule create_aa_alignment:
    input: "alignments/faa/{marker}.faa"
    output: 
        alignments="alignments/faa/muscle/{marker}.faa.efa",
        confidence="alignments/faa/muscle/{marker}.faa.cefa",
        best_alignment="alignments/faa/muscle/{marker}.faa.afa"
    conda: "environments/prokka_envs.yaml"
    threads: 2
    shell: """
        muscle5 -align {input} -stratified -output {output.alignments} -threads 2 > {output.alignments}.log 2>&1
        muscle5 -addconfseq {output.alignments} -output {output.confidence} >> {output.alignments}.log 2>&1
        muscle5 -maxcc {output.alignments} -output {output.best_alignment} >> {output.alignments}.log 2>&1
    """

rule overlay_alignment:
    input: "alignments/faa/muscle/{marker}.faa.afa"
    output: "alignments/fna/muscle/{marker}.fna.afa"
    conda: "environments/prokka_envs.yaml"
    shell: """
        python scripts/overlap_alignment.py -m {wildcards.marker}
    """

rule find_sub_model:
    input: "alignments/fna/muscle/{marker}.fna.afa"
    output: "alignments/fna/muscle/{marker}.iqtree"
    threads: 10
    shell: """
        BASENAME=$(basename "{input}" .fna.afa)
        bin/iqtree2 -s {input} --prefix modeltest/"$BASENAME" -m MF -T {threads} --seqtype DNA
    """

rule ml_tree:
    input: 
        model_file=model_file="alignments/fna/muscle/{marker}.iqtree"
        alignment_file="alignments/fna/muscle/{marker}.fna.afa"
    output: "alignments/fna/muscle/{marker}.mlsearch"
    threads: 10 
    shell: """
        MODEL=$(grep "Best-fit model according to BIC:" {input.model_file} | awk '{print $NF}')
        BASENAME=$(basename "{input.alignment_file}" .fna.afa)

        for i in $(seq 1 50); do
            bin/iqtree2 -s {input.alignment_file} -m "$MODEL" -T AUTO --seqtype DNA --keep-ident --prefix model_selection/"$BASENAME"/"$BASENAME"_"$i" >> "$BASENAME"/"$BASENAME".iqlog 2>&1
        done;
    """