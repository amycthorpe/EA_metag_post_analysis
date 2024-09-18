"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/dereplicate.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: Post-processing bins
"""

MAGS=glob_wildcards(os.path.join(DATA_DIR, "mags/{mags}.fa")).mags

############################################
rule bin_taxqual:
    input:
        os.path.join(RESULTS_DIR, "bins/gtdbtk_final"),
        os.path.join(RESULTS_DIR, "bins/checkm2/quality_report.tsv") 
    output:
        touch("status/bin_taxqual.done")


############################################
# localrules: 


############################################
# Dereplicating the bins

rule checkm:
    input:
        os.path.join(DATA_DIR, "mags/")
    output:
        os.path.join(RESULTS_DIR, "bins/pre_drep_checkm.tsv")
    conda:
        os.path.join(ENV_DIR, "checkm.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/dasTool_checkm.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/dasTool_checkm.err.log")
    threads:
        config["checkm"]["threads"]
    message:
        "Running CheckM on each DasTool output"
    shell:
        """
        (date && checkm lineage_wf -t {threads} -f {output} --tab_table -x fa {input} $(dirname {output})/dastool_checkm && date) 2> {log.err} > {log.out}
        """

rule drep_prepare:
    input:
        check=rules.checkm.output
    output:
        final=os.path.join(RESULTS_DIR, "bins/checkmbeforedrep.tsv"),
        temp=temp(os.path.join(RESULTS_DIR, "bins/checkm_temp.tsv"))
    shell:
        """
        awk '{{print $1".fa,"$13","$14}}' {input.check} | perl -pe "s/Bin.*\n//g" > {output.temp} && 
        (echo "genome,completeness,contamination" && cat {output.temp}) > {output.final}
        """

def symlink_relative(files, input_dir, output_dir):
    """create symlink with and adjust for relative path"""

    input_dir_rel=os.path.relpath(input_dir, output_dir)
    
    for f in files:
        source_path = os.path.join(input_dir_rel, f)
        target_path = os.path.join(output_dir, f)
        os.symlink(source_path, target_path)
        
# rule get_all_bins:
#    input:
#        bins=os.path.join(RESULTS_DIR, "bins/checkmbeforedrep.tsv")
#    output:
#        directory(os.path.join(RESULTS_DIR,"bins/binsbeforedrep")),
#    run:
#        os.mkdir(output[0])
#        bin_folder = os.path.join(os.path.dirname(input.bins), 'das_DASTool_bins')
#        fasta_files = [f for f in os.listdir(bin_folder) if f.endswith(".fa")]
#        symlink_relative(fasta_files, bin_folder, output[0])
           
rule drep:
    input:
        check=rules.drep_prepare.output.final,
        bins=os.path.join(DATA_DIR, "mags/")
    output:
        temp=directory(os.path.join(RESULTS_DIR, "bins/drep/dereplicated_genomes")),
        final=directory(os.path.join(RESULTS_DIR, "bins/finalbins"))
    conda:
        os.path.join(ENV_DIR, "drep.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/drep/drep.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/drep/drep.err.log")
    params:
        comp=config["drep"]["comp"],
        cont=config["drep"]["cont"]
    threads:
        config["drep"]["threads"]
    message:
        "Running dRep on all the bins"
    shell:
        "(date && "
        "dRep dereplicate $(dirname {output.temp}) -p {threads} -comp {params.comp} -con {params.cont} --genomeInfo {input.check} -g {input.bins}/*fa && "
        "cp -r {output.temp} {output.final} && date) 2> {log.err} > {log.out}"

# GTDBTK taxonomy
rule gtdbtk:
    input:
        rules.drep.output.final
    output:
        directory(os.path.join(RESULTS_DIR, "bins/gtdbtk_final"))
    log:
        os.path.join(RESULTS_DIR, "logs/gtdbtk.log")
    conda:
        os.path.join(ENV_DIR, "gtdbtk.yaml")
    params:
        config["gtdbtk"]["path"]
    threads:
        config["gtdbtk"]["threads"]
    message:
        "Running GTDB on MAGs"
    shell:
        "(date && "
        "export GTDBTK_DATA_PATH={params} && gtdbtk classify_wf --cpus {threads} -x fa --genome_dir {input} --out_dir {output} && "
        "date) &> >(tee {log})"

# install checkm database
rule checkm_db:
    output:
        os.path.join(DB_DIR, "CheckM2_database/uniref100.KO.1.dmnd")
    log:
        os.path.join(RESULTS_DIR, "logs/checkm2_db.log")
    conda:
        "/home/amytho/miniconda3/envs/checkm2"
#        os.path.join(ENV_DIR, "checkm.yaml")
    message:
        "Downloading the checkm2 database"
    shell:
        "(date && checkm2 database --download --path $(dirname $(dirname {output})) && date) &> >(tee {log})"
   
# Checking bin quality
rule checkm_final:
    input:
        drep=rules.drep.output.final,
        db=rules.checkm_db.output[0]
    output:
        tsv=os.path.join(RESULTS_DIR, "bins/checkm2/quality_report.tsv")
    conda:
        "/home/amytho/miniconda3/envs/checkm2"  # /hdd0/susbus/tools/conda_envs/checkm2" # "checkm2"
#        os.path.join(ENV_DIR, "checkm.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/checkm/checkm.out.log")
    threads:
        config["checkm"]["threads"]
    params:
        ext=config["checkm"]["extension"],
        db=os.path.join(DB_DIR, "CheckM2_database/uniref100.KO.1.dmnd")
       # checkm2=os.path.join(SUBMODULES, "checkm2/bin/checkm2")
    message:
        "Running Final Checkm on dereplicated output"
    shell:
        "(date && mkdir -p $(dirname $(dirname {output.tsv}))/checkm2_tmp && chmod -R 777 $(dirname $(dirname {output.tsv}))/checkm2_tmp && "
        "export CHECKM2DB={params.db} && "
        "checkm2 predict --tmpdir $(dirname $(dirname {output.tsv}))/checkm2_tmp  --threads {threads} -x {params.ext} --input {input.drep} --output-directory $(dirname {output.tsv}) --force && "
        "date) &> >(tee {log})"
