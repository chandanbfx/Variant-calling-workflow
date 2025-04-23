"""
Quality control and read preprocessing rules
"""

rule Trim:
    """
    Trim low-quality bases using bbduk.
    """
    input:
        #forward=lambda wildcards: samples.loc[wildcards.sample, "fq1"],
        #reversed=lambda wildcards: samples.loc[wildcards.sample, "fq2"]
        fq1=lambda wildcards: samples[wildcards.sample]["fq1"],
        fq2=lambda wildcards: samples[wildcards.sample]["fq2"],
    output:
        forw="trimmed/{sample}_1_Trim.fastq.gz",
        reve="trimmed/{sample}_2_Trim.fastq.gz"
    params:
        trimq=config["phred_cutoff"]
    resources:
        mem_mb=config["qcRam"],
        disk_mb=config["qcDisk"]
    threads: config.get("qcThreads", 8)
    conda: "../envs/preprocessing.yaml"
    log: "logs/trim/{sample}.log"
    shell:
        "bbduk.sh in1={input.fq1} in2={input.fq2} "
        "qtrim=rl trimq={params.trimq} "
        "out1={output.forw} out2={output.reve}"

rule fastqc:
    """
    Run FastQC on trimmed reads.
    """
    input:
        forw="trimmed/{sample}_1_Trim.fastq.gz",
        reve="trimmed/{sample}_2_Trim.fastq.gz"
    output:
        html=["trimmed/{sample}_1_Trim_fastqc.html",
              "trimmed/{sample}_2_Trim_fastqc.html"],
        zip=["trimmed/{sample}_1_Trim_fastqc.zip",
             "trimmed/{sample}_2_Trim_fastqc.zip"]
    resources:
        mem_mb=config["qcRam"],
        disk_mb=config["qcDisk"]
    threads: 1
    conda: "../envs/preprocessing.yaml"
    log: "logs/fastqc/{sample}.log"
    shell:
        "fastqc {input.forw} {input.reve} -o trimmed/"

rule trimAdapter:
    """
    Trim adapter using bbduk.
    """
    input:
        forw="trimmed/{sample}_1_Trim.fastq.gz",
        reve="trimmed/{sample}_2_Trim.fastq.gz",
    output:
        forwar="trimmed/{sample}_1_clean.fastq.gz",
        revers="trimmed/{sample}_2_clean.fastq.gz"
    resources:  mem_mb = config["qcRam"], disk_mb = config["qcDisk"]
    threads: config.get("qcThreads")
    params:
        adapter=config["adapter"]
    conda:
        "environment.yaml"
    log:
        "logs/Adap/{sample}.log"
    shell:
        "bbduk.sh in1={input.forw} in2={input.reve} out1={output.forwar} out2={output.revers} ref={params.adapter} ktrim=r k=21 mink=11 hdist=2 tpe tbo"
