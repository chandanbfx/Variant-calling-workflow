configfile: "config/config.yaml"


import csv

samples = {}

with open("config/samples.tsv") as tsvfile:
    reader = csv.DictReader(tsvfile, delimiter="\t")
    for row in reader:
        sample_id = row["sample"]
        samples[sample_id] = {
            "fq1": row["fq1"],
            "fq2": row["fq2"]
        }
        
#all rule modules in EXECUTION ORDER
include: "rules/common.smk"
include: "rules/preprocessing.smk"
include: "rules/alignment.smk"
include: "rules/post_alignment.smk"
include: "rules/variant_calling.smk"
include: "rules/filtering.smk"
include: "rules/annotation.smk"
include: "rules/utilities.smk"

rule all:
    input:
        "variants/Allvcf_stats.txt",
        "variants/filtered_stats.txt",
        rawVcf="variants/vcf.tar.gz",
        filteredVcf="variants/vcf_filtered.tar.gz",
        dbSNP_Ann="annotation/vcf_dbSNP.tar.gz",
        clinvar_Ann="annotation/vcf_clinVar.tar.gz"
