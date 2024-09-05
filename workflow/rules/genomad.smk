"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-02-27]
Run: snakemake -s workflow/rules/genomad.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To classify viruses and plasmids using geNomad
"""


############################################
rule genomad:
    input:
        expand(os.path.join(RESULTS_DIR, "annotation/genomad/{sid}_summary"), sid=SAMPLES.index)
    output:
        touch("status/genomad.done")


############################################
# localrules: download_genomad


############################################
# Annotating the viruses and plasmids in assemblies
rule download_genomad:
    output:
#        db=os.path.join(DB_DIR, "genomad/genomad_db/version.txt"),
        dummy=os.path.join(DB_DIR, "genomad/db_download.done")
    conda:
        "genomad" # os.path.join(ENV_DIR, "genomad.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/genomad.download_db.log") 
    message:
        "Downloading the geNomad database"
    shell:
        """
        (date
        if [ -e "$(dirname {output.dummy})/genomad_db/version.txt" ]; then
            echo "geNomad database 'genomad_db' already exists."
            touch {output.dummy}
        else
            echo "downloading geNomad database."
            genomad download-database $(dirname {output.dummy})
            touch {output.dummy}
        fi
        date) &> >(tee {log})
        """

rule run_genomad:
    input:
        contig=os.path.join(DATA_DIR, "assembly/{sid}/{sid}.fasta"),
#        db=os.path.join(DB_DIR, "genomad/genomad_db/version.txt")
        dummy=os.path.join(DB_DIR, "genomad/db_download.done")
    output:
        directory(os.path.join(RESULTS_DIR, "annotation/genomad/{sid}_summary")),
    conda:
        "genomad" # os.path.join(ENV_DIR, "genomad.yaml")
    threads:
        config["genomad"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/{sid}_genomad.log")
    params:
        split=config["genomad"]["splits"]
#        db=os.path.join(DB_DIR, "genomad/genomad_db/version.txt"),
    message:
        "Running geNomad on {wildcards.sid} assembly"
    shell:
        "(date && genomad end-to-end --cleanup --splits {params.split} {input.contig} $(dirname {output}) $(dirname {input.dummy})/genomad_db && date) &> >(tee {log})"
