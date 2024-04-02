if "restrict-regions" in config["processing"]:
    rule compose_regions:
        input:
            config["processing"]["restrict-regions"],
        output:
            "results/called/{contig}.regions.bed",
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"

rule haplotype_caller_gvcf:
    input:
        bam=get_sample_bams, # single or list of bam files
        ref=config["genome"]["ref"],
        known=config["genome"]["dbsnp"],
        intervals=(
            "results/called/{contig}.regions.bed"
            if config["processing"].get("restrict-regions")
            else ""
        ),
    output:
        gvcf=protected("results/called/{s}.{contig}.g.vcf.gz"),
    # 	bam="{sample}.assemb_haplo.bam",
    log:
        "logs/gatk/haplotypecaller/{s}.{contig}.log",
    params:
        extra=get_interval_padding(),
    threads: 16,
    resources:
        mem_mb=2048,
    wrapper:
        config["warpper_mirror"]+"bio/gatk/haplotypecaller"

rule combinegvcfs:
    input:
        gvcfs=expand(
            "results/called/{s}.{{contig}}.g.vcf.gz", s=samples.index.get_level_values(0)
        ),
        ref=config["genome"]["ref"],
    output:
        gvcf="results/called/all.{contig}.g.vcf.gz",
    log:
        "logs/gatk/combinegvcfs.{contig}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        config["warpper_mirror"]+"bio/gatk/combinegvcfs"

rule genotype_gvcfs:
    input:
        gvcf="results/called/all.{contig}.g.vcf.gz",  # combined gvcf over multiple samples
	# N.B. gvcf or genomicsdb must be specified
	# in the latter case, this is a GenomicsDB data store
        ref=config["genome"]["ref"],
    output:
        vcf=temp("results/genotyped/all.{contig}.vcf.gz"),
    log:
        "logs/gatk/genotypegvcfs.{contig}.log",
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"], # optional
        java_opts="", # optional
    resources:
        mem_mb=1024
    wrapper:
        config["warpper_mirror"]+"bio/gatk/genotypegvcfs"

rule merge_vcfs:
    input:
        vcfs=lambda w: expand(
            "results/genotyped/all.{contig}.vcf.gz", contig=get_contigs()
        ),
    output:
        "results/genotyped/all.vcf.gz",
    log:
        "logs/picard/mergevcfs.log",
    params:
        extra="",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        config["warpper_mirror"]+"bio/picard/mergevcfs"
