"""
Read alignment and BAM processing rules
"""

rule bwa_index:
    input:
        fa=REF_FA
    output:
        index=touch("trimmed/check.txt")
    resources:
        mem_mb=config["bwaRam"],
        disk_mb=config["bwaDisk"]
    conda: "../envs/alignment.yaml"
    shell:
        """
        bwa index {input.fa}
        """

rule bwa_map:
    """
    Align reads using BWA-MEM.
    """
    input:
        forw="trimmed/{sample}_1_clean.fastq.gz",
        rev="trimmed/{sample}_2_clean.fastq.gz",
        fa=REF_FA,
        index=touch("trimmed/check.txt")
    output:
        temp("mapped_reads/{sample}.bam")
    resources:
        mem_mb=config["bwaRam"],
        disk_mb=config["bwaDisk"]
    threads: config.get("bwathreads", 16)
    conda: "../envs/alignment.yaml"
    log: "logs/bwa/{sample}.log"
    shell:
        "bwa mem -t {threads} {input.fa} {input.forw} {input.rev} | "
        "samtools view -Sb - > {output}"

rule sort_bam:
    """
    Sort BAM file by coordinate.
    """
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.sorted.bam"
    resources:
        mem_mb=130000,
        disk_mb=100000
    threads: 14
    conda: "../envs/alignment.yaml"
    log: "logs/sort/{sample}.log"
    shell:
        "gatk --java-options -Xmx56G SortSam "
        "-I {input} -O {output} -SO coordinate"
