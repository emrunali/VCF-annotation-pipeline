#!/bin/bash

################################################################################
# VCF Annotation and Filtering Pipeline
# Version: 1.0.0
#
# Description:
#   Automated pipeline for annotating VCF files with SnpEff and filtering variants
#   based on functional impact and clinical significance.
#   Designed for learning clinical bioinformatics workflows.
#
# Dependencies:
#   - bcftools (v1.15+)
#   - tabix
#   - Java (for SnpEff)
#   - SnpEff (installed locally)
#
# Usage:
#   ./vcf_annotation_filter.sh -i data/input.vcf.gz -o results/sample_name [-s /path/to/snpeff]
#
################################################################################

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# Default parameters
SNPEFF_DIR="$HOME/tools/snpEff"
GENOME_DB=""  # Will be auto-detected
OUTPUT_PREFIX=""
INPUT_VCF=""
LOG_FILE=""
VERBOSE=0
RESULTS_DIR="results"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

################################################################################
# Functions
################################################################################

log_message() {
    local level=$1
    shift
    local message="$@"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    # Only write to log file if LOG_FILE is set and not empty
    if [ -n "${LOG_FILE}" ]; then
        echo "[${timestamp}] [${level}] ${message}" | tee -a "${LOG_FILE}"
    else
        echo "[${timestamp}] [${level}] ${message}"
    fi
    
    case $level in
        ERROR)
            echo -e "${RED}[ERROR]${NC} ${message}" >&2
            ;;
        WARNING)
            echo -e "${YELLOW}[WARNING]${NC} ${message}" >&2
            ;;
        SUCCESS)
            echo -e "${GREEN}[SUCCESS]${NC} ${message}"
            ;;
    esac
}

usage() {
    cat << EOF
VCF Annotation and Filtering Pipeline (SnpEff)

Usage: $0 -i INPUT_VCF -o OUTPUT_PREFIX [OPTIONS]

Required Arguments:
    -i INPUT_VCF        Input VCF file (can be gzipped)
                        Recommended: place in data/ directory
    -o OUTPUT_PREFIX    Output file prefix (without directory)
                        Results will be saved to results/ directory

Optional Arguments:
    -r RESULTS_DIR     Results directory (default: results)
    -s SNPEFF_DIR      Path to SnpEff directory (default: ~/tools/snpEff)
    -v                 Verbose mode
    -h                 Show this help message

Example:
    # Recommended directory structure:
    # ./vcf_annotation_filter.sh
    # ./data/sample.vcf.gz
    # ./results/ (created automatically)
    
    $0 -i data/sample.vcf.gz -o my_sample

Output Files (in results/ directory):
    - {PREFIX}_annotated.vcf.gz       : Fully annotated VCF
    - {PREFIX}_high_impact.vcf.gz     : High impact variants
    - {PREFIX}_moderate_impact.vcf.gz : Moderate impact variants
    - {PREFIX}_coding.vcf.gz          : All coding variants
    - {PREFIX}_snpeff_summary.html    : SnpEff statistics (open in browser)
    - {PREFIX}_summary.txt            : Processing summary
    - {PREFIX}_pipeline.log           : Detailed log file
    - {PREFIX}_genes.txt              : List of affected genes

EOF
    exit 1
}

setup_directories() {
    log_message "INFO" "Setting up directory structure..."
    
    # Create results directory if it doesn't exist
    if [ ! -d "${RESULTS_DIR}" ]; then
        mkdir -p "${RESULTS_DIR}"
        log_message "INFO" "Created results directory: ${RESULTS_DIR}"
    fi
    
    log_message "SUCCESS" "Directory structure ready"
}

check_dependencies() {
    log_message "INFO" "Checking dependencies..."
    
    local missing_deps=0
    
    # Check bcftools and tabix
    for cmd in bcftools tabix java; do
        if ! command -v $cmd &> /dev/null; then
            log_message "ERROR" "Required tool not found: $cmd"
            missing_deps=1
        else
            if [ "$cmd" = "java" ]; then
                local version=$(java -version 2>&1 | head -n1)
            else
                local version=$($cmd --version 2>&1 | head -n1 || echo "unknown")
            fi
            log_message "INFO" "Found $cmd: $version"
        fi
    done
    
    # Check SnpEff
    if [ ! -f "${SNPEFF_DIR}/snpEff.jar" ]; then
        log_message "ERROR" "SnpEff not found at ${SNPEFF_DIR}"
        log_message "ERROR" "Please install SnpEff or specify path with -s flag"
        missing_deps=1
    else
        log_message "INFO" "Found SnpEff at ${SNPEFF_DIR}"
        
        # Auto-detect available GRCh38 database
        log_message "INFO" "Detecting available genome databases..."
        if [ -d "${SNPEFF_DIR}/data" ]; then
            local available_dbs=$(ls "${SNPEFF_DIR}/data" | grep -i "grch38\|hg38" | head -1)
            if [ -n "$available_dbs" ]; then
                GENOME_DB="$available_dbs"
                log_message "INFO" "Found genome database: ${GENOME_DB}"
            else
                log_message "WARNING" "No GRCh38/hg38 database found in ${SNPEFF_DIR}/data"
                log_message "INFO" "Available databases:"
                ls "${SNPEFF_DIR}/data" | head -5 | while read db; do
                    log_message "INFO" "  - $db"
                done
                log_message "ERROR" "Please download a database with: java -jar ${SNPEFF_DIR}/snpEff.jar download GRCh38.99"
                missing_deps=1
            fi
        fi
    fi
    
    if [ $missing_deps -eq 1 ]; then
        log_message "ERROR" "Missing required dependencies."
        exit 1
    fi
    
    log_message "SUCCESS" "All dependencies found"
}

validate_input() {
    log_message "INFO" "Validating input file: ${INPUT_VCF}"
    
    if [ ! -f "${INPUT_VCF}" ]; then
        log_message "ERROR" "Input file not found: ${INPUT_VCF}"
        exit 1
    fi
    
    # Check if file is gzipped
    if [[ "${INPUT_VCF}" == *.gz ]]; then
        if ! gunzip -t "${INPUT_VCF}" 2>/dev/null; then
            log_message "ERROR" "Input file is corrupted or not a valid gzip file"
            exit 1
        fi
    fi
    
    # Validate VCF format
    if ! bcftools view -h "${INPUT_VCF}" &> /dev/null; then
        log_message "ERROR" "Invalid VCF format"
        exit 1
    fi
    
    # Count variants
    local variant_count=$(bcftools view -H "${INPUT_VCF}" | wc -l)
    log_message "INFO" "Input VCF contains ${variant_count} variants"
    
    log_message "SUCCESS" "Input validation complete"
}

annotate_with_snpeff() {
    local input=$1
    local output=$2
    
    log_message "INFO" "Starting SnpEff annotation..."
    log_message "INFO" "Using genome database: ${GENOME_DB}"
    log_message "INFO" "This may take several minutes depending on variant count..."
    
    # Run SnpEff - redirect stderr to log, stdout to output
    java -Xmx4g -jar "${SNPEFF_DIR}/snpEff.jar" \
        -v "${GENOME_DB}" \
        -stats "${RESULTS_DIR}/$(basename ${OUTPUT_PREFIX})_snpeff_summary.html" \
        -csvStats "${RESULTS_DIR}/$(basename ${OUTPUT_PREFIX})_snpeff_summary.csv" \
        "${input}" \
        > "${output}" \
        2>> "${LOG_FILE}"
    
    local exit_code=$?
    
    if [ $exit_code -eq 0 ] && [ -f "${output}" ] && [ -s "${output}" ]; then
        log_message "SUCCESS" "SnpEff annotation complete"
        
        # Count annotated variants
        local annotated_count=$(grep -v "^#" "${output}" | wc -l)
        log_message "INFO" "Successfully annotated ${annotated_count} variants"
    else
        log_message "ERROR" "SnpEff annotation failed with exit code ${exit_code}"
        if [ -f "${output}" ]; then
            log_message "ERROR" "Output file size: $(wc -c < ${output}) bytes"
        fi
        exit 1
    fi
    
    # Compress and index
    bgzip -f "${output}"
    tabix -p vcf "${output}.gz"
    log_message "INFO" "Output compressed and indexed"
}

filter_high_impact() {
    local input=$1
    local output=$2
    
    log_message "INFO" "Filtering for HIGH impact variants..."
    
    # Filter for HIGH impact variants from SnpEff annotation
    bcftools view -i 'ANN ~ "HIGH"' \
        -O z \
        -o "${output}" \
        "${input}"
    
    tabix -p vcf "${output}"
    
    local filtered_count=$(bcftools view -H "${output}" | wc -l)
    log_message "SUCCESS" "High impact filtering complete: ${filtered_count} variants retained"
}

filter_moderate_impact() {
    local input=$1
    local output=$2
    
    log_message "INFO" "Filtering for MODERATE impact variants..."
    
    # Filter for MODERATE impact variants
    bcftools view -i 'ANN ~ "MODERATE"' \
        -O z \
        -o "${output}" \
        "${input}"
    
    tabix -p vcf "${output}"
    
    local filtered_count=$(bcftools view -H "${output}" | wc -l)
    log_message "SUCCESS" "Moderate impact filtering complete: ${filtered_count} variants retained"
}

filter_coding_variants() {
    local input=$1
    local output=$2
    
    log_message "INFO" "Filtering for coding region variants..."
    
    # Filter for coding variants (missense, nonsense, frameshift, splice, etc.)
    bcftools view -i 'ANN ~ "missense_variant" || ANN ~ "stop_gained" || ANN ~ "stop_lost" || ANN ~ "frameshift_variant" || ANN ~ "splice" || ANN ~ "start_lost" || ANN ~ "inframe"' \
        -O z \
        -o "${output}" \
        "${input}"
    
    tabix -p vcf "${output}"
    
    local filtered_count=$(bcftools view -H "${output}" | wc -l)
    log_message "SUCCESS" "Coding variant filtering complete: ${filtered_count} variants retained"
}

extract_gene_list() {
    local input=$1
    local output=$2
    
    log_message "INFO" "Extracting affected genes list..."
    
    # Extract gene names from ANN field
    bcftools query -f '%ANN\n' "${input}" | \
    grep -oP '\|[^\|]+\|' | \
    tr -d '|' | \
    sort -u | \
    grep -v "^$" > "${output}" 2>/dev/null || touch "${output}"
    
    local gene_count=$(wc -l < "${output}" 2>/dev/null || echo 0)
    log_message "INFO" "Found ${gene_count} unique affected genes"
}

generate_summary() {
    local input_vcf=$1
    local annotated_vcf=$2
    local high_impact_vcf=$3
    local moderate_impact_vcf=$4
    local coding_vcf=$5
    local summary_file=$6
    
    log_message "INFO" "Generating summary report..."
    
    local input_count=$(bcftools view -H "${input_vcf}" | wc -l)
    local total_variants=$(bcftools view -H "${annotated_vcf}" | wc -l)
    local high_impact_variants=$(bcftools view -H "${high_impact_vcf}" | wc -l)
    local moderate_impact_variants=$(bcftools view -H "${moderate_impact_vcf}" | wc -l)
    local coding_variants=$(bcftools view -H "${coding_vcf}" | wc -l)
    
    # Calculate percentages
    local high_pct=0
    local mod_pct=0
    local coding_pct=0
    if [ ${total_variants} -gt 0 ]; then
        high_pct=$(awk "BEGIN {printf \"%.1f\", (${high_impact_variants}/${total_variants})*100}")
        mod_pct=$(awk "BEGIN {printf \"%.1f\", (${moderate_impact_variants}/${total_variants})*100}")
        coding_pct=$(awk "BEGIN {printf \"%.1f\", (${coding_variants}/${total_variants})*100}")
    fi
    
    cat > "${summary_file}" << EOF
================================================================================
VCF Annotation and Filtering Pipeline - Summary Report
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')
Input File: ${INPUT_VCF}
Output Prefix: $(basename ${OUTPUT_PREFIX})
Results Directory: ${RESULTS_DIR}

Annotation Tool: SnpEff
Reference Assembly: ${GENOME_DB}
SnpEff Location: ${SNPEFF_DIR}

Variant Statistics:
-------------------
Input Variants: ${input_count}
Total Annotated Variants: ${total_variants}

Impact Distribution:
--------------------
HIGH Impact Variants: ${high_impact_variants} (${high_pct}%)
  - Stop gained/lost, frameshift, splice site disruptions
  - Likely loss of function or significant protein disruption

MODERATE Impact Variants: ${moderate_impact_variants} (${mod_pct}%)
  - Missense variants, in-frame indels
  - May affect protein function

Coding Variants (All): ${coding_variants} (${coding_pct}%)
  - All variants affecting protein coding sequence

Output Files:
-------------
1. $(basename ${OUTPUT_PREFIX})_annotated.vcf.gz
   - All variants with complete SnpEff functional annotations
   - Includes impact predictions, gene names, transcript IDs
   
2. $(basename ${OUTPUT_PREFIX})_high_impact.vcf.gz
   - HIGH impact variants only (priority for clinical review)
   - Stop codons, frameshifts, splice site disruptions
   
3. $(basename ${OUTPUT_PREFIX})_moderate_impact.vcf.gz
   - MODERATE impact variants (missense, in-frame indels)
   - Require additional evidence for pathogenicity assessment
   
4. $(basename ${OUTPUT_PREFIX})_coding.vcf.gz
   - All coding region variants combined
   - Useful for gene panel analysis
   
5. $(basename ${OUTPUT_PREFIX})_snpeff_summary.html
   - Interactive HTML report with charts and statistics
   - Open in web browser for detailed visualization
   
6. $(basename ${OUTPUT_PREFIX})_snpeff_summary.csv
   - Machine-readable statistics for downstream analysis
   
7. $(basename ${OUTPUT_PREFIX})_genes.txt
   - List of all genes affected by variants
   - One gene per line, useful for pathway analysis
   
8. $(basename ${OUTPUT_PREFIX})_pipeline.log
   - Complete processing log with timestamps
   - Audit trail for clinical laboratory compliance

Tools Used:
-----------
- SnpEff $(java -jar ${SNPEFF_DIR}/snpEff.jar -version 2>&1 | head -1 | awk '{print $3}') for functional annotation
- bcftools $(bcftools --version | head -1 | awk '{print $2}') for VCF manipulation
- Functional impact predictions based on Sequence Ontology terms

Clinical Relevance:
-------------------
HIGH impact variants are prioritized for clinical review as they typically
result in protein truncation or loss of function. MODERATE impact variants
(missense, in-frame indels) require further assessment based on:
  - Conservation scores (phyloP, phastCons)
  - Protein domain location and function
  - Population frequency (gnomAD MAF < 0.01)
  - Clinical databases (ClinVar pathogenic/likely pathogenic)
  - In silico predictions (SIFT, PolyPhen-2, CADD)

Recommended Next Steps:
-----------------------
1. Review HIGH impact variants in disease-relevant genes
2. Filter by population frequency using gnomAD (MAF < 0.01)
3. Check ClinVar for known pathogenic variants
4. Assess MODERATE impact variants in candidate genes
5. Prioritize based on inheritance pattern and phenotype
6. Validate clinically significant variants via Sanger sequencing

Pipeline Version: 1.0.0
Execution Time: Complete log available in $(basename ${OUTPUT_PREFIX})_pipeline.log
================================================================================
EOF
    
    log_message "SUCCESS" "Summary report generated: ${summary_file}"
    cat "${summary_file}"
}

################################################################################
# Main Pipeline
################################################################################

# Parse command line arguments
while getopts "i:o:r:s:vh" opt; do
    case $opt in
        i) INPUT_VCF="$OPTARG" ;;
        o) OUTPUT_PREFIX="$OPTARG" ;;
        r) RESULTS_DIR="$OPTARG" ;;
        s) SNPEFF_DIR="$OPTARG" ;;
        v) VERBOSE=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required arguments
if [ -z "${INPUT_VCF}" ] || [ -z "${OUTPUT_PREFIX}" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Prepend results directory to output prefix if not already included
if [[ ! "${OUTPUT_PREFIX}" =~ ^${RESULTS_DIR}/ ]]; then
    OUTPUT_PREFIX="${RESULTS_DIR}/$(basename ${OUTPUT_PREFIX})"
fi

# Set up log file BEFORE calling any functions that use it
LOG_FILE="${OUTPUT_PREFIX}_pipeline.log"

# Create results directory first if it doesn't exist
mkdir -p "${RESULTS_DIR}"

# Create/truncate log file
: > "${LOG_FILE}"

# Pipeline start
log_message "INFO" "=========================================="
log_message "INFO" "VCF Annotation & Filtering Pipeline v1.0.0"
log_message "INFO" "Using SnpEff for Functional Annotation"
log_message "INFO" "=========================================="

# Run pipeline steps
setup_directories
check_dependencies
validate_input

# Define output files
ANNOTATED_VCF="${OUTPUT_PREFIX}_annotated.vcf"
HIGH_IMPACT_VCF="${OUTPUT_PREFIX}_high_impact.vcf.gz"
MODERATE_IMPACT_VCF="${OUTPUT_PREFIX}_moderate_impact.vcf.gz"
CODING_VCF="${OUTPUT_PREFIX}_coding.vcf.gz"
GENES_FILE="${OUTPUT_PREFIX}_genes.txt"
SUMMARY_FILE="${OUTPUT_PREFIX}_summary.txt"

# Execute pipeline
annotate_with_snpeff "${INPUT_VCF}" "${ANNOTATED_VCF}"
filter_high_impact "${ANNOTATED_VCF}.gz" "${HIGH_IMPACT_VCF}"
filter_moderate_impact "${ANNOTATED_VCF}.gz" "${MODERATE_IMPACT_VCF}"
filter_coding_variants "${ANNOTATED_VCF}.gz" "${CODING_VCF}"
extract_gene_list "${ANNOTATED_VCF}.gz" "${GENES_FILE}"
generate_summary "${INPUT_VCF}" "${ANNOTATED_VCF}.gz" "${HIGH_IMPACT_VCF}" "${MODERATE_IMPACT_VCF}" "${CODING_VCF}" "${SUMMARY_FILE}"

log_message "SUCCESS" "Pipeline complete! All outputs generated successfully."
log_message "INFO" "Results saved to: ${RESULTS_DIR}/"
log_message "INFO" "Check ${OUTPUT_PREFIX}_pipeline.log for detailed processing information"
log_message "INFO" "View ${RESULTS_DIR}/$(basename ${OUTPUT_PREFIX})_snpeff_summary.html for detailed annotation statistics"

exit 0
