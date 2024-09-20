"""
Author: Susheel Bhanu BUSI & Amy Thorpe
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-07-01]
Run: snakemake -s workflow/rules/metabolishmm.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: To run metabolishmm on MAGs
"""

MAGS=glob_wildcards(os.path.join(RESULTS_DIR, "bins/finalbins/{mags}.fa")).mags


############################################
rule metabolishmm:
    input:
        os.path.join(RESULTS_DIR, "metabolishmm/summarize_metabolism_out")  
    output:
        touch("status/metabolishmm.done")


############################################
localrules: 


############################################
# rules to run metabolisHMM
rule metadata:
    input:
        stat=os.path.join(RESULTS_DIR, "bins/gtdbtk_final/gtdbtk.bac120.summary.tsv")
    output:
        os.path.join(RESULTS_DIR, "data/metadata")
    shell:
        """cat {input} | awk '{{print $1\",\"$2}}' | sed '/^user/d' | sed 's/.fasta.contigs//g' | sed 's/.fasta_sub.contigs//g' | awk '!visited[$0]++' > {output}"""

rule get_data:
    output:
        directory("curated_markers"),
        dummy="status/metabolishmm_DB.done"
    message:
        "downloading the marker DB for metabolisHMM"
    shell:
        "wget https://github.com/elizabethmcd/metabolisHMM/releases/download/v1.9/metabolisHMM_markers_v1.9.tgz && tar -xzvf metabolisHMM_markers_v1.9.tgz && "
        "touch {output.dummy}"

rule prep_bins:
    input:
        bins=os.path.join(RESULTS_DIR, "bins/finalbins"),
        dummy="status/metabolishmm_DB.done"
    output:
        bins=directory(os.path.join(RESULTS_DIR, "mags/fna_mags")),
        dummy="status/bin_prep.done"
    log:
        os.path.join(RESULTS_DIR, "logs/metabolishmm_bin_prep.log")
    message:
        "Changing bins from .FA to .FNA"
    shell:
        "(date && mkdir -p {output.bins} &&"
        """for file in {input.bins}/*.fa; do name=$(basename $file .fa).fna; cp -v $file {output.bins}/$name; done && """
        "touch {output.dummy} && date) &> >(tee {log})"

rule run_metabolishmm:
    input:
        bins=os.path.join(RESULTS_DIR, "mags/fna_mags"),
        meta=rules.metadata.output,
        dummy="status/bin_prep.done"
    output:
        sum=directory(os.path.join(RESULTS_DIR, "metabolishmm/summarize_metabolism_out")),
        indiv=directory(os.path.join(RESULTS_DIR, "metabolishmm/individual_metabolism_out"))
    conda:
        os.path.join(ENV_DIR, "metabolishmm_new.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/metabolishmm.log")
    message:
        "Running MetabolisHMM on all MAGs"
    shell:
        "(date && "
        "summarize-metabolism --input {input.bins} --output {output.sum} --metadata {input.meta} --summary summarize_metabolism.csv --heatmap summarize_metabolism.pdf --aggregate ON --plotting ON && "
        "summarize-metabolism --input {input.bins} --output {output.indiv} --metadata {input.meta} --summary individual_metabolism.csv --heatmap individual_metabolism.pdf --plotting ON && "
        "date) &> >(tee {log})"
