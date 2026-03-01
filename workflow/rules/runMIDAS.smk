# Run the MIDAS species module to generate species profiles.
rule profileSpeciesAbundances:
    input:
        r1=join(config["filterdir"],"{sample}-filtered.1.fastq.gz"),
        r2=join(config["filterdir"],"{sample}-filtered.2.fastq.gz")
    output:
        profile="workflow/out/midasOutput/{sample}/species/species_profile.txt"
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/MIDASpython2-no-builds.yml"
    shell:
        """
        run_midas.py species workflow/out/midasOutput/{wildcards.sample} \
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
        expand("workflow/out/midasOutput/HouseholdTransmission-Stool-{sample}/species/species_profile_subject.txt",sample=samples)
    output:
        "workflow/out/midasOutput/species/species_profile_all.txt"
    shell:
        # Concatenate all species profiles without the file headers.
        "cat <( tail -n +2 {input} ) | grep -v '==>' | sed '/^$/d' | awk '$3 > 0' |"
        # Add the file header back and include the new "sample" field.
        "sed '1 i\species_id\tcount_reads\tcoverage\trelative_abundance\tsample'> {output}"

# Identify the abundant species to analyze for SNPs in each subject.
rule identifyAbundantSpeciesBySubject:
    input:
        lambda wildcards: expand("workflow/out/midasOutput/HouseholdTransmission-Stool-{subject}-{sample}/species/species_profile.txt",
            subject=wildcards.subject,
            sample=dict(tuple(df.groupby(['subject'])))[wildcards.subject]['timepoint'].tolist())
    params:
        minCoverage=config["runMIDAS_speciesMinCoverage"]
    output:
        "workflow/out/midasOutput/species/abundantSpecies_{subject}.txt"
    shell:
        # Retain species with coverage greater than the minimum coverage threshold.
        "cat <( tail -n +2 {input} ) | grep -v '==>' | sed '/^$/d' | awk '$3 > {params.minCoverage}"
            "' | cut -f1 | sort | uniq > {output}"
