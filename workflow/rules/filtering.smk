rule gatk_select:
    input:
        vcf="results/genotyped/all.vcf.gz",
        ref=config["genome"]["ref"],
    output:
        vcf=temp("results/filtered/all.{v}.vcf.gz"),
    log:
        "logs/gatk/selectvariants/{v}.log",
    params:
        extra=get_vartype_arg,  # optional filter arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        config["warpper_mirror"]+"bio/gatk/selectvariants"

rule gatk_filter:
    input:
        vcf="results/filtered/all.{v}.vcf.gz",
        ref=config["genome"]["ref"],
#        intervals="targets.bed",
    output:
        vcf=temp("results/filtered/all.{v}.hardfiltered.vcf.gz"),
    log:
        "logs/gatk/variantfiltration/{v}.log",
    params:
        filters=get_filter,
        extra="",  # optional arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        config["warpper_mirror"]+"bio/gatk/variantfiltration"

rule haplotype_caller:
    input:
        vcf="results/filtered/all.{v}.vcf.gz",
        ref=config["genome"]["ref"],
        fai=config["genome"]["fai"],
        dict=config["genome"]["dict"],
        mills=config["genome"]["mills"],
        mills_idx=config["genome"]["mills_idx"],
        omni=config["genome"]["omni"],
        omni_idx=config["genome"]["omni_idx"],
        g1k=config["genome"]["g1k"],
        g1k_idx=config["genome"]["g1k_idx"],
        dbsnp=config["genome"]["dbsnp"],
        dbsnp_idx=config["genome"]["dbsnp_idx"],
    output:
        vcf=temp("results/filtered/all.{v}.recalibrated.vcf.gz"),
        idx=temp("results/filtered/all.{v}.recalibrated.vcf.gz.idx"),
        tranches=temp("results/filtered/all.{v}.recalibrated.vcf.gz.tranches"),
    log:
        "logs/gatk/variantrecalibrator/{v}.log",
    params:
        mode="BOTH",  # set mode, must be either SNP, INDEL or BOTH
        resources={
            "mills": {"known": False, "training": True, "truth": True, "prior": 15.0},
            "omni": {"known": False, "training": True, "truth": False, "prior": 12.0},
            "g1k": {"known": False, "training": True, "truth": False, "prior": 10.0},
            "dbsnp": {"known": True, "training": False, "truth": False, "prior": 2.0},
        },
        annotation=["MQ", "QD", "SB"],
        extra=config["params"]["gatk"]["VariantRecalibrator"],
    threads: 16
    resources:
        mem_mb=1024,
    wrapper:
         config["warpper_mirror"]+"bio/gatk/variantrecalibrator"

rule merge_calls:
    input:
        vcfs=expand(
            "results/filtered/all.{v}.{filtertype}.vcf.gz",
            v=["snvs", "indels"],
            filtertype="recalibrated"
            if config["filtering"]["vqsr"]
            else "hardfiltered",
        ),
    output:
        "results/filtered/all.vcf.gz",
    log:
        "logs/picard/merge-filtered.log",
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