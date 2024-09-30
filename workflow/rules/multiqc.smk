"""
Author: Susheel Bhanu BUSI & Amy Thorpe
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-09-23]
Run: snakemake -s workflow/rules/multiqc.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: To run multiqc on all samples
"""

############################################
rule multiqc:
    input:
        os.path.join(RESULTS_DIR, "multiqc/multiqc_report.html")
    output:
        touch("status/multiqc.done")


############################################
localrules:


############################################
rule multiqc_fastqc:
    input:
        os.path.join(RESULTS_DIR, "multiqc_input/")
    output:
        html=os.path.join(RESULTS_DIR, "multiqc/multiqc_report.html"),
        stat=os.path.join(RESULTS_DIR, "multiqc/multiqc_data/multiqc_fastqc.txt"),
    log:
        os.path.join(RESULTS_DIR, "multiqc/multiqc.log")
    threads:
        config["multiqc"]["threads"]
    conda:
        os.path.join(ENV_DIR, "multiqc.yaml")
    message:
        "MultiQC (FastQC)"
    shell:
        "multiqc --interactive -p -f -m fastqc -o $(dirname {output.html}) {input} &> {log}"
