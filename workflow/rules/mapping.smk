rule fastp_se:
    input:
        sample=get_fastq,
    output:
        trimmed=temp("results/trimmed/{s}{u}.fastq.gz"),
        html="report/{s}{u}.fastp.html",
        json="report/{s}{u}.fastp.json",
    log:
        "logs/fastp/{s}{u}.log",
    threads: 8,
    wrapper:
        config["warpper_mirror"]+"bio/fastp"

rule fastp_pe:
    input:
        sample=get_fastq,
    output:
        trimmed=[temp("results/trimmed/{s}{u}.1.fastq.gz"), temp("results/trimmed/{s}{u}.2.fastq.gz")],
        html="report/{s}{u}.fastp.html",
        json="report/{s}{u}.fastp.json",
    log:
        "logs/fastp/{s}{u}.log",
    threads: 8
    wrapper:
        config["warpper_mirror"]+"bio/fastp"

rule fastp_pe_wo_trimming:
    input:
        sample=get_fastq,
    output:
        html="report/pe_wo_trimming/{s}{u}.fastp.html",
        json="report/pe_wo_trimming/{s}{u}.fastp.json",
    log:
        "logs/fastp/pe_wo_trimming/{s}{u}.log",
    threads: 8,
    wrapper:
        config["warpper_mirror"]+"bio/fastp"

rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=multiext(config["genome"]["ref"]+".64", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "results/mapped/{s}{u}.sorted.bam",
    log:
        "logs/bwa_mem/{s}{u}.log",
    params:
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
    threads: 16,
    wrapper:
        config["warpper_mirror"]+"bio/bwa/mem"

rule mark_duplicates_spark:
    input:
        "results/mapped/{s}{u}.sorted.bam",
    output:
        bam=temp("results/dedup/{s}{u}.bam"),
        metrics="results/dedup/{s}{u}.metrics.txt",
    log:
        "logs/dedup/{s}{u}.log",
    params:
        extra="--remove-sequencing-duplicates",  # optional
        java_opts="",  # optional
        #spark_runner="",  # optional, local by default
        #spark_master="",  # optional
        #spark_extra="", # optional
    resources:
        # Memory needs to be at least 471859200 for Spark, so 589824000 when
        # accounting for default JVM overhead of 20%. We round round to 650M.
        mem_mb=lambda wildcards, input: max([input.size_mb * 0.25, 650]),
    threads: 16
    wrapper:
         config["warpper_mirror"]+"bio/gatk/markduplicatesspark"

rule gatk_baserecalibrator:
    input:
        bam="results/dedup/{s}{u}.bam",
        ref=config["genome"]["ref"],
        dict=config["genome"]["dict"],
        known=config["genome"]["dbsnp"],  # optional known sites - single or a list
    output:
        recal_table="results/recal/{s}{u}.grp",
    log:
        "logs/gatk/baserecalibrator/{s}{u}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        config["warpper_mirror"]+"bio/gatk/baserecalibrator"

rule gatk_applybqsr:
    input:
        bam="results/dedup/{s}{u}.bam",
        ref=config["genome"]["ref"],
        dict=config["genome"]["dict"],
        recal_table="results/recal/{s}{u}.grp",
    output:
        bam=protected("results/recal/{s}{u}.bam"),
    log:
        "logs/gatk/gatk_applybqsr/{s}{u}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
        embed_ref=True,  # embed the reference in cram output
    resources:
        mem_mb=1024,
    wrapper:
        config["warpper_mirror"]+"bio/gatk/applybqsr"

rule samtools_index:
    input:
        "results/recal/{s}{u}.bam",
    output:
        "results/recal/{s}{u}.bam.bai",
    log:
        "logs/samtools_index/{s}{u}.log",
    params:
        extra="",  # optional params string
    threads: 16  # This value - 1 will be sent to -@
    wrapper:
        config["warpper_mirror"]+"bio/samtools/index"
