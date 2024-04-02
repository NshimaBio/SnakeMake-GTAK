rule annotate_variants:
    input:
        calls="results/filtered/all.vcf.gz",  # .vcf, .vcf.gz or .bcf
        cache="resources/vep/cache",  # can be omitted if fasta and gff are specified
        plugins="resources/vep/plugins",
        # optionally add reference genome fasta
        fasta=config["genome"]["ref"],
        fai=config["genome"]["fai"], # fasta index
        gff=config["genome"]["gff"],
        csi=config["genome"]["gff_idx"], # tabix index
        # add mandatory aux-files required by some plugins if not present in the VEP plugin directory specified above.
        # aux files must be defined as following: "<plugin> = /path/to/file" where plugin must be in lowercase
        # revel = path/to/revel_scores.tsv.gz
    output:
        calls=report("results/annotated/all.vcf.gz",
            caption="../report/vcf.rst",
            category="Calls",
        ),  # .vcf, .vcf.gz or .bcf
        stats=report("results/stats/all.stats.html",
            caption="../report/stats.rst",
            category="Calls",
        ),
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["params"]["vep"]["plugins"],
        extra=config["params"]["vep"]["extra"],  # optional: extra arguments
    log:
        "logs/vep/annotate.log",
    threads: 8
    wrapper:
        config["warpper_mirror"]+"v1.25.0/bio/vep/annotate"
