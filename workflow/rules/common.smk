import numpy as np
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

report: "../report/workflow.rst"

container: "mambaorg/micromamba:alpine"

###### Config file and sample sheets #####
configfile: "config/config.yaml"

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], dtype=str).set_index(
    ["Sample", "Unit"], drop=False
)

samples = pd.read_csv("config/samples.csv", dtype=str).fillna(value="")
samples.set_index(["Sample", "Unit"], drop=False,inplace=True)

samples.index = samples.index.set_levels(
    [i.astype(str) for i in samples.index.levels]
)  # enforce str in index
validate(samples, schema="../schemas/samples.schema.yaml")
# validate(samples, schema="workflow/schemas/samples.schema.yaml")

##### Wildcard constraints #####
units=samples.index.get_level_values(1)
wildcard_constraints:
    v="snvs|indels",
    s="|".join(samples.index.get_level_values(0)),
    u="|".join(units)
# "" if units[0] else 
##### Helper functions #####

# contigs in reference genohime
def get_contigs():
    return pd.read_table(config["genome"].get("fai"), header=None, usecols=[0], dtype=str)

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = samples.loc[(wildcards.s, wildcards.u), ["fq1", "fq2"]]
    if config["fastq"].get("pe"):
        return [fastqs.fq1, fastqs.fq2]
    return [fastqs.fq1]

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{s}{u}\tSM:{s}\tLB:{lb}\tPL:{pl}'".format(
        s=wildcards.s,
        u=wildcards.u,
        lb=config["fastq"].get("omics"),
        pl=config["fastq"].get("platform")
    )

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if config["fastq"].get("pe"):
        # paired-end sample
        return expand(
            "results/trimmed/{s}{u}.{group}.fastq.gz",
            group=[1, 2],**wildcards
        )
    # single end sample
    return "results/trimmed/{s}{u}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/recal/{s}{u}.bam",
        s=wildcards.s,
        u=samples.loc[wildcards.s,:].Unit,
    )


def get_interval_padding(config=config,default=""):
    padding = config["processing"].get("region-padding")
    if padding:
        return "--interval-padding {}".format(padding)
    return default


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "results/mapped/{s}{u}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "results/dedup/{s}{u}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f


def get_snpeff_reference():
    return "{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"])


def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.v == "snvs" else "INDEL"
    )


def get_filter(wildcards):
    return {"snv-hard-filter": config["filtering"]["hard"][wildcards.v]}
