# Germline Variant Calling Pipeline (Snakemake)

A modular and reproducible variant calling pipeline built using [Snakemake](https://snakemake.readthedocs.io/).
Designed for short-read NGS data and tested on **chr1** of the human genome.

---

## Features

- Adapter trimming and quality control (BBDuk)
- Read alignment using BWA
- Variant discovery using GATK
- SNP/INDEL filtration and VQSR
- Annotation using dbSNP, ClinVar, and others
- Modular Snakemake rules and conda-managed environments

---

## Project Structure

```
.
├── Snakefile                 # Main workflow entry point
├── config/                  # Configuration files (samples, params)
│   ├── config.yaml
│   └── samples.tsv
├── envs/                    # Conda environment definitions
│   ├── preprocessing.yaml
│   ├── alignment.yaml
│   ├── variant_calling.yaml
│   ├── annotation.yaml
│   └── environment.yaml
├── rules/                   # Modular rule files (e.g. preprocessing, alignment)
├── scripts/                 # Custom scripts used in workflow
├── resources/               # BED files, adapter files
├── reference/               # chr1.fa and related reference indices
├── database/                # dbSNP, ClinVar, VQSR resources (not pushed to GitHub)
├── data/                    # Raw FASTQ files (ignored by Git)
├── results/                 # Output files and VCFs (ignored by Git)
├── logs/                    # Log files (ignored by Git)
└── LICENSE
```

---

## Setup Instructions

### 1. Clone the Repository

```bash
git clone https://github.com/chandanbfx/germline-variant-calling-pipeline.git
cd germline-variant-calling-pipeline
```

### 2. Install Mamba (Optional but Recommended)

```bash
conda install -n base -c conda-forge mamba
```

### 3. Run a Dry Test

```bash
snakemake -n --use-conda --cores 4
```

### 4. Run the Pipeline

```bash
snakemake --use-conda --cores 8
```

---

## Configuration

Located in `config/config.yaml`, this file controls:
- Sample sheet: `config/samples.tsv`
- Adapter file: `resources/adapters.fa`
- Reference genome: `reference/chr1.fa`
- BED file: `resources/chr1_only.bed`
- Annotation and VQSR files: under `database/`

---

## External Databases (Not Included)

The following large files must be placed in `database/` manually:

- **dbSNP**: `00-common_all.vcf.gz`
- **ClinVar**: `clinvar_20250421.vcf.gz`
- **VQSR**: Mills, HapMap, 1000G, Axiom

You can download these from:
- [Broad GATK Resource Bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/)
- [dbSNP FTP](https://ftp.ncbi.nih.gov/snp/)
- [ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/)

---

---

# External Link for
## Reference Files

You can download `chr1.fa` from the following link:
[Download chr1.fa from Google Drive] (https://drive.google.com/file/d/10ao_pvFXoG58Bpi2wusxXLQgB7BlPdW-/view?usp=sharing)

---

## Test Dataset

For lightweight testing, this project uses:
- Paired-end FASTQs aligned to `chr1`
- Subset BED file for `chr1`

---

## License

MIT License

---

## Acknowledgements

- [Snakemake](https://snakemake.readthedocs.io/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [GATK](https://gatk.broadinstitute.org/)
- [SnpEff/SnpSift](http://snpeff.sourceforge.net/)
- [BBMap](https://sourceforge.net/projects/bbmap/)
- [FastQc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
