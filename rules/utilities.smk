"""
Utility rules for VCF processing and statistics generation
"""

rule generate_stats:
    """
    Generate final statistics reports for raw and filtered VCFs
    """
    input:
        raw_vcf="variants/vcf/all.vcf.gz",
        filtered_vcf="variants/filtered/all_filtered.vcf.gz",
        annotated_vcf="annotation/all.dbSnpANN.vcf"
    output:
        all_stats="variants/Allvcf_stats.txt",       # Must match rule all
        filtered_stats="variants/filtered_stats.txt"  # Must match rule all
    shell:
        """
        rtg vcfstats {input.filtered_vcf} > {output.filtered_stats}
        rtg vcfstats {input.annotated_vcf} > {output.all_stats}
        """

rule split_raw_vcfs:
    """
    Split raw joint-called VCF into per-sample VCFs
    """
    input:
        vcf="variants/vcf/all.vcf.gz",
        index="variants/vcf/all.vcf.gz.tbi"
    output:
        tar="variants/vcf.tar.gz",
        touchfile=touch("variants/rawvcfFile/check.txt")
    params:
        outdir="variants/rawvcfFile"
    conda:
        "../envs/variant_calling.yaml"
    script:
        "../scripts/split_vcf.py"

rule split_filtered_vcfs:
    """
    Split filtered VCF into per-sample VCFs
    """
    input:
        vcf="variants/filtered/all_filtered.vcf.gz",
        index="variants/filtered/all_filtered.vcf.gz.tbi"
    output:
        tar="variants/vcf_filtered.tar.gz",
        touchfile=touch("variants/vcfFiltered/check.txt")
    params:
        outdir="variants/vcfFiltered"
    conda:
        "../envs/variant_calling.yaml"
    script:
        "../scripts/split_vcf.py"

rule split_annotated_vcfs:
    """
    Split annotated VCF into per-sample VCFs with dbSNP annotations
    """
    input:
        vcf="annotation/all.dbSnpANN.vcf"
    output:
        tar="annotation/vcf_dbSNP.tar.gz",
        touchfile=touch("annotation/dbSNP/check.txt")
    params:
        outdir="annotation/dbSNP"
    conda:
        "../envs/annotation.yaml"
    script:
        "../scripts/split_vcf.py"

rule split_clinvar_vcfs:
    """
    Split annotated VCF into per-sample VCFs with ClinVar annotations
    """
    input:
        vcf="annotation/all.clinvarANN.vcf"
    output:
        tar="annotation/vcf_clinVar.tar.gz",
        touchfile=touch("annotation/clinVar/check.txt")
    params:
        outdir="annotation/clinVar"
    conda:
        "../envs/annotation.yaml"
    script:
        "../scripts/split_vcf.py"
