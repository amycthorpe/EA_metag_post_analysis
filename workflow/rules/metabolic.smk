"""
Author: Susheel Bhanu BUSI & Amy Thorpe
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-09-23]
Run: snakemake -s workflow/rules/metabolic.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: To run metabolic on MAGs
"""

############################################
rule metabolic:
    input:
        os.path.join(RESULTS_DIR, "metabolic/450_samples_out")
    output:
        touch("status/metabolic.done")


############################################
localrules:


############################################
rule metabolic_reads:
    input:
        reads=os.path.join(DATA_DIR, "reads")
    output:
        reads=os.path.join(RESULTS_DIR, "metabolic/reads.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/metabolic/reads.log")
    conda:
        "METABOLIC_v4.0"
    message:
        "Setting up reads file for metabolic"
    shell:
        "(date && for file in {input.reads}/*; do echo $file; done | paste -d, - - > {output.reads} && "
        "head {output.reads} && date) &> >(tee {log})"

rule run_metabolic:
    input:
        bin=os.path.join(RESULTS_DIR, "bins/finalbins"),
        reads=os.path.join(RESULTS_DIR, "metabolic/reads.txt")
    output:
        metabolic=directory(os.path.join(RESULTS_DIR, "metabolic/450_samples_out"))
    log:
        os.path.join(RESULTS_DIR, "logs/metabolic/run_metabolic.log")
    conda:
        "METABOLIC_v4.0"
    params:
        metabolic=config["metabolic"]["path"]
    message:
        "Running metabolic"
    shell:
        "(date && perl {params.metabolic}/METABOLIC-C.pl -in-gn {input.bin} -r {input.reads} -o {output.metabolic} && "
        "date) &> >(tee {log})"
