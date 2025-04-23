"""
Variant Calling Rules - GATK Best Practices
"""

include: "../rules/common.smk"

rule haplotype_caller:
    """
    Generate GVCFs for each sample using GATK HaplotypeCaller
    """
    input:
        fa=config["reference"],
        bam="sorted_reads/{sample}.markdup.bam",
        interval=config.get("bedfile", "")
    output:
        vcf="variants/{sample}.gvcf.gz",
        index="variants/{sample}.gvcf.gz.tbi"
    params:
        extra_args="--ERC GVCF",
        sample_name=lambda w: config["samples"].get(w.sample, {}).get("name", w.sample) if "samples" in config 
               else samples_df.loc[w.sample, "sample_name"] if "sample_name" in samples_df.columns 
               else w.sample,
        mem_gb=lambda wildcards, resources: int(resources.mem_mb) // 1024,
        interval_arg=lambda wildcards, input: f"-L {input.interval}" if input.interval else ""
    resources:
        mem_mb=config["gatkRam"],
        disk_mb=config["gatkDisk"]
    conda:
        "../envs/variant_calling.yaml"
    log:
        "logs/haplotypecaller/{sample}.log"
    threads: config.get("gatkThreads", 4)
    shell:
        """
        gatk --java-options "-Xmx{params.mem_gb}G" HaplotypeCaller \
        -I {input.bam} \
        -R {input.fa} \
        {params.interval_arg} \
        -O {output.vcf} \
        {params.extra_args}
        """

rule combine_gvcfs:
    """
    Combine multiple GVCFs into one
    """
    input:
        fa=config["reference"],
        gvcfs=expand("variants/{sample}.gvcf.gz", sample=SAMPLES),
        indices=expand("variants/{sample}.gvcf.gz.tbi", sample=SAMPLES)
    output:
        combined="variants/combined.g.vcf.gz",
        index="variants/combined.g.vcf.gz.tbi"
    params:
        mem_gb=lambda wildcards, resources: int(resources.mem_mb) // 1024,
        gvcf_args=lambda wildcards, input: " ".join(f"-V {vcf}" for vcf in input.gvcfs)
    resources:
        mem_mb=config.get("gatkRam", 64000)
    conda:
        "../envs/variant_calling.yaml"
    log:
        "logs/combinegvcfs/combined.log"
    shell:
        """
        gatk --java-options "-Xmx{params.mem_gb}G" CombineGVCFs \
        -R {input.fa} \
        {params.gvcf_args} \
        -O {output.combined}
        """

rule genotype_gvcfs:
    """
    Joint genotype GVCFs
    """
    input:
        fa=config["reference"],
        gvcf="variants/combined.g.vcf.gz",
        index="variants/combined.g.vcf.gz.tbi"
    output:
        vcf="variants/vcf/all.vcf.gz",
        index="variants/vcf/all.vcf.gz.tbi"
    params:
        mem_gb=lambda wildcards, resources: int(resources.mem_mb) // 1024,
    resources:
        mem_mb=config.get("gatkRam", 64000)
    conda:
        "../envs/variant_calling.yaml"
    log:
        "logs/genotypegvcfs/joint.log"
    shell:
        """
        gatk --java-options "-Xmx{params.mem_gb}G" GenotypeGVCFs \
        -R {input.fa} \
        -V {input.gvcf} \
        -O {output.vcf}
        """
