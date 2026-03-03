# Run the MIDAS species module to generate species profiles.
rule profileSpeciesAbundances:
    input:
        # r1=lambda wildcards:
        #    dict(df.groupby['sample'])))[wildcards.sample]['fastq1'],
        # r2=lambda wildcards:
        #    dict(df.groupby['sample'])))[wildcards.sample]['fastq2']
        r1=lambda wc: fastq_map[wc.sample]["fastq1"],
        r2=lambda wc: fastq_map[wc.sample]["fastq2"],
    output:
        profile="workflow/out/midasOutput/{sample}/species/species_profile.txt"
    threads: config['maxCPUs']
    # conda:
    #    "MIDASpython2"
    # use conda activate MIDASpython2 before running snakemake (or in the SLURM job)
    shell:
        r"""
        export PYTHONPATH=/dfs7/xuek5-lab/aespelet/tools/MIDAS:$PYTHONPATH
        /dfs7/xuek5-lab/aespelet/tools/MIDAS/scripts/run_midas.py species workflow/out/midasOutput/{wildcards.sample} \
        -1 {input.r1} -2 {input.r2} -t {threads}
        """

# Annotate each species profile with the sample name for further processing.
rule annotateSpeciesAbundancesBySubject:
    input:
        "workflow/out/midasOutput/{sample}/species/species_profile.txt"
    output:
        "workflow/out/midasOutput/{sample}/species/species_profile_subject.txt"
    shell:
        "sed 's/$/\t{wildcards.sample}/g' {input} > {output}"

# Concatenate all species abundance profiles into a single file.
rule concatenateSpeciesAbundances:
   input:
       expand(
           "workflow/out/midasOutput/{sample}/species/species_profile_subject.txt",
           sample=samples
       )
   output:
       "workflow/out/midasOutput/species/species_profile_all.txt"
   shell:
       r"""
       mkdir -p $(dirname {output})
       # print a single header, then append all non-header rows from all files
       printf "species_id\tcount_reads\tcoverage\trelative_abundance\tsample\n" > {output}
       awk 'FNR>1' {input} | awk '$3 > 0' >> {output}
       """

