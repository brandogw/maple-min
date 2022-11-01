include: "rules/common.smk"

rule all:
    input:
        expand("data/processed/mutation_stats/{sample}",sample = samples["Sample"]),

# Modules
include: "rules/aligning.smk"
include: "rules/analysis.smk"
