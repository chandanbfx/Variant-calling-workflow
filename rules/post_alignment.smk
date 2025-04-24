#genearte summary of alignment metrics 
rule AlignmentSummary:
    input:
        fa=REF_FA,
        bam="sorted_reads/{sample}.sorted.bam"
    output:
        "mapped_reads/{sample}.AlignmentSummary.txt",
    resources: mem_mb=config["qcRam"], disk_mb=config["qcDisk"]
    conda:
        "environment.yaml"
    log:
        "logs/gatkAlignment/{sample}.log"
    threads: 8
    shell:
        "gatk --java-options -Xmx{resources.mem_mb}M CollectAlignmentSummaryMetrics -I {input.bam} -R {input.fa} -O {output}"
rule mark_duplicates:
    input:
        bam="sorted_reads/{sample}.sorted.bam",
        metrices="mapped_reads/{sample}.AlignmentSummary.txt"
    output:
        bam="sorted_reads/{sample}.markdup.bam",
        metrics="sorted_reads/{sample}.markduplicate.txt"
    resources:
        mem_mb=config["gatkRam"],
        disk_mb=config["gatkDisk"]
    shell:
        "gatk MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics}"

rule base_recalibration:
    input:
        bam="sorted_reads/{sample}.markdup.bam",
        known=config["known_Sites"]
    output:
        table="AddRG/{sample}.recal.table"
    shell:
        "gatk BaseRecalibrator -I {input.bam} -R {config[reference]} --known-sites {input.known} -O {output.table}"

rule apply_recalibration:
    input:
        fa=REF_FA,
        bam="AddRG/{sample}.markdup.bam",
        table="AddRG/{sample}.recal.table"
    output:
        "AddRG/{sample}.recalibrated.bam"
    resources:  mem_mb = config["gatkRam"], disk_mb = config["gatkDisk"]
    threads: config.get("gatkThreads")
    conda:
        "environment.yaml"
    log:
        "logs/gatkApplyBQSR/{sample}.log"
    shell:
        "gatk --java-options -Xmx{resources.mem_mb}M ApplyBQSR -I {input.bam} -R {input.fa} --bqsr-recal-file {input.table} -O {output}"
