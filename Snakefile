# This is a snakefile for analysing RenSeq variants.
# First, pacbio hifi reads are trimmed with cutadapt.
# Then, the trimmed reads are aligned to the reference genome with minimap2.
# The aligned reads are sorted and indexed with samtools.
# Variants are called with gatk HaplotypeCaller.

# Import the config file
configfile: "config.yaml"

# Read the samples file into a pandas dataframe
samples = pd.read_table(config["samples"], header=0).set_index(["sample"], drop=False)

# function to get the sample name from the wildcards
def get_reads(wildcards):
    return samples["Reads"][wildcards.sample]

rule all:
    input:
        expand("aligned_reads/{sample}.bam", sample=samples["sample"])

rule cutadapt:
    input:
        get_reads
    params:
        fiveprime=config["fiveprime"],
        threeprime=config["threeprime"]
    output:
        fq=temp("trimmed_reads/{sample}.fq"),
        intermediate=temp("trimmed_reads/{sample}_intermediate.fq")
    log:
        "logs/cutadapt/{sample}.log"
    threads:
        8
    conda:
        "envs/variant_analysis.yaml"
    resources:
        mem_mb=4000
    shell:
        """(cutadapt -j {threads} -g ^{params.fiveprime} -o {output.intermediate} {input}) 2> {log}
        (cutadapt -j {threads} -a {params.threeprime}$ -o {output.fq} {output.intermediate}) 2>> {log}"""

rule minimap2:
    input:
        fq="trimmed_reads/{sample}.fq",
        reference=config["reference"]
    output:
        bam="aligned_reads/{sample}.bam"
    log:
        "logs/minimap2/{sample}.log"
    threads:
        8
    conda:
        "envs/variant_analysis.yaml"
    resources:
        mem_mb=8000
    shell:
        """minimap2 -ax map-hifi -t {threads} {input.reference} {input.fq} | samtools sort -@ {threads} -o {output.bam} -"""

rule gatk:
    input:
        bam="aligned_reads/{sample}.bam",
        reference=config["reference"]
    output:
        vcf="variants/{sample}.vcf"
    log:
        "logs/gatk/{sample}.log"
    threads:
        8
    conda:
        "envs/variant_analysis.yaml"
    resources:
        mem_mb=8000
    shell:
        """gatk HaplotypeCaller -R {input.reference} -I {input.bam} -O {output.vcf}"""

# This rule is for annotating the variants with snpEff
# the snpEff database will need to be generated manually beforehand
# hardcoded for now
rule snpeff:
    input:
        vcf="variants/{sample}.vcf"
    output:
        vcf="annotated_variants/{sample}.vcf"
    log:
        "logs/snpeff/{sample}.log"
    threads:
        4
    resources:
        mem_mb=4000
    conda:
        "envs/variant_analysis.yaml"
    shell:
        """snpEff solanum_verrucosum {input.vcf} > {output.vcf} 2> {log}"""

# snpsift interval on all nlr genes
rule snpsift:
    input:
        vcf="annotated_variants/{sample}.vcf"
    output:
        vcf="annotated_variants/{sample}_nlr.vcf"
    params:
        nlr=config["nlr"]
    log:
        "logs/snpsift/{sample}.log"
    threads:
        4
    resources:
        mem_mb=4000
    conda:
        "envs/variant_analysis.yaml"
    shell:
        """cat {input.vcf} | SnpSift intervals {params.nlr} > {output.vcf} 2> {log}"""

# now make a table of the variants using snpsift extractfields.
# this will be a tab delimited file with the following columns:
# CHROM POS REF ALT QUAL FILTER AC AF AN DP
rule snpsift_table:
    input:
        vcf="annotated_variants/{sample}_nlr.vcf"
    output:
        table="tables/{sample}_nlr.tsv"
    log:
        "logs/snpsift_table/{sample}.log"
    threads:
        4
    resources:
        mem_mb=4000
    conda:
        "envs/variant_analysis.yaml"
    shell:
        """SnpSift extractFields {input.vcf} CHROM POS REF ALT QUAL FILTER AC AF AN DP > {output.table} 2> {log}"""

