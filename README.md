## Overview

This pipeline automates the processing of VCF (Variant Call Format) files through functional annotation and multi-tier filtering to identify clinically relevant variants. It is designed with clinical laboratory best practices in mind, including comprehensive logging, audit trails, and standardized quality control metrics.

## Features

- Automated variant annotation using SnpEff for reliable effect prediction

- Tiered filtering to categorize variants by functional impact

- Clinical relevance focus on protein-affecting variants

- Detailed, auditable reports with logs, HTML summaries, and statistics

- Gene-level outputs for downstream pathway or enrichment analysis

## Dependencies

### Required Tools

- `bcftools` (v1.15+) - VCF manipulation and filtering
- `tabix` - VCF indexing
- `Java` (v8+) - Required for SnpEff
- `SnpEff` - Variant annotation tool

### Installation
#### 1. Install bcftools and tabix
**Ubuntu/Debian:**
```
sudo apt-get update
sudo apt-get install -y bcftools tabix
```
**macOS (with Homebrew):**

```
brew install bcftools htslib
```
**Conda (any OS):**
```
conda install -c bioconda bcftools htslib
```
#### 2. Install Java
**Ubuntu/Debian:**
```
sudo apt-get install -y default-jre
```
**macOS:**
```
brew install openjdk
```
#### 3. Install SnpEff
```
# Create tools directory
mkdir -p ~/tools
cd ~/tools

# Download SnpEff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff

# Test installation
java -jar snpEff.jar -version

# Download human genome database (GRCh38)
java -jar snpEff.jar download -v GRCh38.99
```

## Directory Structure
```
vcf_annotation_pipeline/
├── vcf_annotation_filter.sh    # Main pipeline script
├── data/                        # Input VCF files
│   └── sample.vcf.gz
└── results/                     # Output files (auto-created)
    ├── sample_annotated.vcf.gz
    ├── sample_high_impact.vcf.gz
    ├── sample_moderate_impact.vcf.gz
    ├── sample_coding.vcf.gz
    ├── sample_genes.txt
    ├── sample_snpeff_summary.html
    ├── sample_summary.txt
    └── sample_pipeline.log
```

## Usage
### Basic Usage
```
# Make script executable
chmod +x vcf_annotation_filter.sh

# Run pipeline
./vcf_annotation_filter.sh -i data/sample.vcf.gz -o my_sample
```

### Command Line Options
```
Required Arguments:
  -i INPUT_VCF        Input VCF file (can be gzipped)
  -o OUTPUT_PREFIX    Output file prefix (without directory)

Optional Arguments:
  -r RESULTS_DIR      Results directory (default: results)
  -s SNPEFF_DIR       Path to SnpEff directory (default: ~/tools/snpEff)
  -v                  Verbose mode
  -h                  Show help message
```

### Example Commands
```
# Basic usage with default settings
./vcf_annotation_filter.sh -i data/patient001.vcf.gz -o patient001

# Specify custom SnpEff location
./vcf_annotation_filter.sh -i data/sample.vcf.gz -o sample -s /opt/snpeff

# Use custom output directory
./vcf_annotation_filter.sh -i data/sample.vcf.gz -o sample -r custom_results

# Verbose mode for debugging
./vcf_annotation_filter.sh -i data/sample.vcf.gz -o sample -v
```

## Output Files

- `{PREFIX}_annotated.vcf.gz` – Full annotated VCF with SnpEff predictions (ANN field includes gene, transcript, effect, impact). Use for comprehensive review or submissions.

- `{PREFIX}_high_impact.vcf.gz` – Stop-gain/loss, frameshift, and splice-site variants. Use for priority clinical review.

- `{PREFIX}_moderate_impact.vcf.gz` – Missense and in-frame indel variants. Use for secondary review and further validation.

- `{PREFIX}_coding.vcf.gz` – All coding-region variants (HIGH + MODERATE). Use for gene panel or protein-coding studies.

- `{PREFIX}_snpeff_summary.html / .csv` – Interactive and machine-readable summaries with variant distributions and impact statistics.

- `{PREFIX}_genes.txt` – Alphabetical list of affected genes for pathway or enrichment analysis.

- `{PREFIX}_summary.txt` – Overall processing summary with QC metrics and filtering results.

- `{PREFIX}_pipeline.log` – Timestamped log for reproducibility and troubleshooting.

## Interpreting Results

- **HIGH Impact (🔴)** – Likely loss-of-function (nonsense, frameshift, splice-site). Review first; often pathogenic.

- **MODERATE Impact (🟡)** – Protein-altering (missense, in-frame indels). Review with supporting evidence.

- **LOW/MODIFIER (🟢)** – Synonymous or intronic; generally benign.

### Suggested Workflow:

1. Review HIGH impact variants in relevant genes; confirm with clinical databases (ClinVar, OMIM).

2. Filter rare variants (gnomAD_AF < 0.01) using bcftools.

3. Evaluate MODERATE impact variants with conservation and in silico predictors.

4. Validate key findings (Sanger, segregation, or functional assays).

## Example: 1000 Genomes Data
### Download Test Data
```
# Create data directory
mkdir -p data

# Download chromosome 22 from 1000 Genomes
cd data
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi

# Create test subset
bcftools view 1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz | head -2500 | bgzip > test_sample.vcf.gz
tabix -p vcf test_sample.vcf.gz
cd ..

# Run pipeline
./vcf_annotation_filter.sh -i data/test_sample.vcf.gz -o test_1000g
```

## Future Enhancements

- Integration with ClinVar annotations
- gnomAD population frequency filtering
- Support for multi-sample VCFs
- Variant prioritization scoring
- IGV batch script generation
- ACMG classification framework

## Acknowledgments

- 1000 Genomes Project for publicly available test data
- SnpEff team for the excellent annotation tool
- Samtools/bcftools developers for robust VCF processing tools




