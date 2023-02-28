# This is a snakefile for analysing RenSeq variants.
# First, pacbio hifi reads are trimmed with cutadapt.
# Then, the trimmed reads are aligned to the reference genome with minimap2.
# The aligned reads are sorted and indexed with samtools.
# Variants are called with gatk HaplotypeCaller.

import pandas as pd

# Import the config file
configfile: "config.yaml"

# Read the samples file into a pandas dataframe
samples = pd.read_csv(config["samples"], header=0).set_index(["sample"], drop=False)

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
        mem_mb=16000
    shell:
        """minimap2 -ax map-hifi -t {threads} {input.reference} {input.fq} | samtools sort -@ {threads} -o {output.bam} -"""
