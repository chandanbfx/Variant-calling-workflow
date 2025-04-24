#VQSR in indel mode for indels quality score recalibration.
rule vqsr_indel:
    input:
        call="variants/vcf/all.vcf.gz",
        index="variants/vcf/all.vcf.gz.tbi",
        Mills=config["mills"],
        Mills_index=config["mills_index"],
        Axiom=config["Axiom"],
        Axiom_index=config["Axiom_index"],
        dbsnp=config["known_Sites"],
        dbsnp_index=config["known_Sites_index"],
    output:
        calls="variants/vcf/all_indel.recal",
        calls_index="variants/vcf/all_indel.recal.idx",
        tranch="variants/vcf/all_indels.tranches"
    params:
        tranche=config["indel_tranche_value"],
        annotation=config["indel_annotation"],
        max_gaussians=config["indel_max_gaussians"]	
    resources:  mem_mb = config["gatkRam"], disk_mb = config["gatkDisk"]
    threads: config.get("gatkThreads")
    conda:
        "variant_calling.yaml"
    log:
        "logs/vqsr_indel/all.log"
    shell:
        "gatk --java-options '-Xmx24G -Xms24G' VariantRecalibrator -V {input.call} --trust-all-polymorphic {params.tranche} {params.annotation} -mode INDEL --max-gaussians {params.max_gaussians} -resource:mills,known=false,training=true,truth=true,prior=12 {input.Mills} -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {input.Axiom} -resource:dbsnp,known=true,training=false,truth=false,prior=2 {input.dbsnp} -O {output.calls} --tranches-file {output.tranch}"

rule apply_vqsr_indel:
    input:
        call="variants/vcf/all.vcf.gz",
        index="variants/vcf/all.vcf.gz.tbi",
        recal="variants/vcf/all_indel.recal",
        recal_index="variants/vcf/all_indel.recal.idx",
        tranch="variants/vcf/all_indels.tranches"
    output:
        filtered="variants/filtered/all_indelfiltered.vcf.gz",
        filtered_index="variants/filtered/all_indelfiltered.vcf.gz.tbi",
    resources:  mem_mb = config["gatkRam"], disk_mb = config["gatkDisk"]
    threads: config.get("gatkThreads")
    params:
        indelTranches=config["indel_tranches"]
    conda:
        "variant_calling.yaml"
    log:
        "logs/apply_vqsr_indel/all.log"
    shell:
        "gatk --java-options '-Xmx5g -Xms5g' ApplyVQSR -V {input.call} --recal-file {input.recal} --tranches-file {input.tranch} --truth-sensitivity-filter-level {params.indelTranches} --create-output-variant-index true -mode INDEL -O {output.filtered}"

#VQSR in snp mode for snp quality score recalibration
rule vqsr_snp:
    input:
        vcf="variants/filtered/all_indelfiltered.vcf.gz",
        hapmap=config["hapmap"],
        hapmap_index=config["hapmap_index"],
        omni=config["1000G_omni"],
        omni_index=config["omni_index"],
        phase_1000=config["1000G_phase"],
        phase_1000_index=config["1000G_phase_index"],
        dbSNP=config["known_Sites"],
        dbSNP_index=config["known_Sites_index"]       
    output:
        recal="variants/vcf/all_snp.recal",
        tranches="variants/vcf/all_snp.tranches"
    conda:
        "variant_calling.yaml"
    shell:
        "gatk VariantRecalibrator -V {input.vcf} -mode SNP "
        "-resource:hapmap,known=false,training=true,truth=true,prior=15 {input.hapmap} "
        "-resource:1000G_phase,known=false,training=true,truth=true,prior=15 {input.phase_1000} "
        "-resource:omni,known=false,training=true,truth=true,prior=15 {input.omni} "
        "-resource:dbSNP,known=false,training=true,truth=true,prior=15 {input.dbSNP} "
        "-O {output.recal} --tranches-file {output.tranches}"
rule apply_vqsr_snp:
    input:
        vcf="variants/filtered/all_indelfiltered.vcf.gz",
        recal="variants/vcf/all_snp.recal",
        tranches="variants/vcf/all_snp.tranches"
    output:
        "variants/filtered/all_filtered.vcf.gz",
        "variants/filtered/all_filtered.vcf.gz.tbi"
    conda:
        "variant_calling.yaml"
    shell:
        "gatk ApplyVQSR -V {input.vcf} --recal-file {input.recal} "
        "--tranches-file {input.tranches} -O {output}"
