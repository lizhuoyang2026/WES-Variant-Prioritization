# WES Clinical Pipeline for Inborn Errors of Immunity (IEI)

## Overview
This repository contains an automated, end-to-end Whole Exome Sequencing (WES) computational pipeline developed to identify and prioritize pathogenic rare variants in patients with Inborn Errors of Immunity (IEI). 

The pipeline is optimized for High-Performance Computing (HPC) environments using the SLURM workload manager. It takes raw sequencing data (FASTQ) through alignment, variant calling, comprehensive clinical annotation, and targeted gene panel filtering to output a highly curated list of candidate actionable variants.

## Pipeline Architecture
The workflow is executed via a unified Bash shell script (`wes_pipeline.sh`) divided into two primary phases:

### Phase 1: Data Preparation
1. **Read Alignment:** Raw paired-end reads (`.fastq.gz`) are aligned to the human reference genome (hg38) using `BWA-MEM`. 
2. **Sorting & Indexing:** Alignments are piped directly into `samtools` for coordinate sorting and indexing.
3. **Duplicate Marking:** PCR and optical duplicates are identified and removed using `GATK MarkDuplicates` to prevent artificial amplification biases.

### Phase 2: Variant Discovery & Prioritization
4. **Variant Calling:** Germline short variants (SNPs and Indels) are identified using `GATK HaplotypeCaller` to generate a raw VCF.
5. **Clinical Annotation:** The raw VCF is annotated using `ANNOVAR`. The pipeline integrates multiple databases to assess population frequency, pathogenicity, and functional impact:
   - `refGene` (Gene-based annotation)
   - `clinvar` (Clinical significance)
   - `1000g2015aug_all` & `exac03` (Population allele frequencies)
   - `dbnsfp47a` (In silico prediction scores)
   - `avsnp150` (dbSNP identifiers)
6. **Clinical IEI Filtering:** The heavily annotated variant list is cross-referenced against a custom, curated Inborn Errors of Immunity gene panel (`clean_iuis_panel.txt`). Variants outside this panel are filtered out, resulting in a clinically actionable candidate list.

### Phase 3: Downstream Analysis & Visualization (R/RMarkdown)
7. **Automated Reporting (`R_Analysis&Visualization.Rmd`):** A comprehensive R script processes the annotated outputs to perform:
   - **Genotype Quality Control:** Evaluates read depth (DP) and allele fractions (AF) directly from the VCF.
   - **Constraint Metrics Integration:** Maps candidate genes against gnomAD constraint metrics (pLI vs LOEUF) to assess Loss-of-Function intolerance.
   - **Data Visualization:** Generates Nature-style, publication-ready composite figures (SVG/PDF/TIFF) detailing functional categories, allele frequencies, chromosome distributions, and sequential filtering waterfalls.
   - **Final Prioritization:** Outputs a dynamically generated PDF report and a final consolidated CSV table mapping variants to specific IEI inheritance patterns.

## Environment & Dependencies

### 1. HPC Environment (Phase 1 & 2)
This script is designed for a Linux-based HPC cluster running the SLURM scheduler. 
* **Environment Manager:** Conda 
* **Core Bioinformatics Tools Required:**
  * BWA
  * Samtools
  * GATK4 (Genome Analysis Toolkit)
  * ANNOVAR
  * Perl, Bash, and standard UNIX utilities (`grep`, `head`)

### 2. R Environment (Phase 3)
* **Language:** R (v4.2+)
* **Core Packages:** `tidyverse`, `vcfR`, `ggplot2`, `patchwork`, `ggrepel`, `kableExtra`
* **Graphics Engines:** `svglite`, `ragg`, `showtext`

## Usage

### 1. Configuration
Before running, ensure the directory variables at the top of the script (`PROJECT_DIR`, `REF`, `ANNOVAR_DIR`, `DATA_DIR`) point to your specific HPC file paths. 

Verify that your IEI gene panel (`clean_iuis_panel.txt`) and raw FASTQ files (`SAMPLE_1_1.fastq.gz`, `SAMPLE_1_2.fastq.gz`) are placed in the correct directories.

### 2. Execution
Submit the batch script to the SLURM scheduler:
```bash
sbatch wes_pipeline.sh
```

### 3. Generate the RMarkdown Report
Once the bash pipeline finishes, compile the downstream analysis report by running the RMarkdown script in your local or server R environment, supplying the target sample ID:

```R
rmarkdown::render("R_Analysis&Visualization.Rmd", params = list(sample_id = "SAMPLE_1"))
```
