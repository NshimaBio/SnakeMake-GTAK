# rule fastp_multiqc:
#     input:
#         expand("report/{s}{u}.fastp.json", s=samples.index.get_level_values(0),u=samples.index.get_level_values(1)),
#     output:
#         "report/fastp_multiqc.html",
#     params:
#         extra="",  # Optional: extra parameters for multiqc.
#     log:
#         "logs/fastp/fastp_multiqc.log",
#     wrapper:
#         config["warpper_mirror"]+"bio/multiqc"
rule samtools_stats:
    input:
        bam="results/recal/{s}{u}.bam",
        # bed="",              #Optional input, specify target regions
    output:
        "results/qc/samtools_stats/{s}{u}.txt",
    params:
        extra="",  # Optional: extra arguments.
        region="",  # Optional: region string.
    log:
        "logs/samtools-stats/{s}{u}.log",
    wrapper:
        config["warpper_mirror"]+"bio/samtools/stats"

rule multiqc:
    input:
        expand(
            [
                "results/qc/samtools_stats/{t.Sample}{t.Unit}.txt",
                "report/{t.Sample}{t.Unit}.fastp.json",
                "results/qc/dedup/{t.Sample}{t.Unit}.metrics.txt",
            ],
            t=samples.itertuples(),
        ),
    output:
        report(
            "report/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "logs/mapping_multiqc.log",
    wrapper:
        config["warpper_mirror"]+"bio/multiqc"