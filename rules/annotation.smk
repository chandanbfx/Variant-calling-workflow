"""
Variant annotation rules
"""

rule snpeff_annotation:
    """
    Annotate variants with SnpEff.
    """
    input:
        vcf="variants/filtered/all_filtered.vcf.gz"
    output:
        ann_vcf="annotation/all.snpEffAnn.vcf",
        stats="annotation/snpEff_summary.html",
        genes="annotation/all.summary.genes.txt"
    resources:
        mem_mb=config["siftRam"],
        disk_mb=config["siftDisk"]
    threads: config.get("siftThreads", 8)
    conda: "../envs/annotation.yaml"
    log: "logs/snpeff/annotation.log"
    shell:
        "snpEff -Xmx8G -s {output.stats} GRCh37.75 {input.vcf} > {output.ann_vcf}"

rule snpsift_dbSNP:
    """
    Add dbSNP annotations using SnpSift.
    """
    input:
        vcf="annotation/all.snpEffAnn.vcf",
        db=config["dbSNP"]
    output:
        "annotation/all.dbSnpANN.vcf"
    resources:
        mem_mb=config["siftRam"],
        disk_mb=config["siftDisk"]
    threads: config.get("siftThreads", 8)
    conda: "../envs/annotation.yaml"
    log: "logs/snpsift/dbnsfp.log"
    shell:
        "SnpSift dbSNP -v -db {input.db} {input.vcf} > {output}"

rule snpSift_clinvar:
	input:
		vari="annotation/all.snpEffAnn.vcf",
		clinvar=config["clinVar"],
		index=config["clinVar_index"]
	output:
		Anno="annotation/all.clinvarANN.vcf"
	conda:
		"environment.yaml"
	log:
		"logs/sift_clinvar/all.log"
	resources:  mem_mb = config["siftRam"], disk_mb = config["siftDisk"]
	threads: config.get("siftThreads")
	shell:
		"SnpSift -Xmx16G annotate -V {input.clinvar} {input.vari} > {output}"

