#!/bin/bash
#SBATCH --job-name=WES_Pipeline
#SBATCH --partition=cpu-epyc-genoa
#SBATCH --qos=long
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@example.com
#SBATCH --output=/path/to/project/logs/wes_pipeline_%j.log
#SBATCH --error=/path/to/project/logs/wes_pipeline_%j.err

echo "Starting End-to-End Clinical WES Pipeline..."
date

# --- WAKE UP THE SOFTWARE (CONDA) ---
# Ensure miniconda path matches your cluster environment
source /path/to/miniconda3/bin/activate fyp_env

# --- 1. SET UP VARIABLES ---
PROJECT_DIR="/path/to/project"
REF="$PROJECT_DIR/reference/hg38.fasta"
ANNOVAR_DIR="$PROJECT_DIR/annovar"
PANEL="$PROJECT_DIR/clean_iuis_panel.txt"

# Ensure these match the exact filenames you uploaded to your data folder
DATA_DIR="/path/to/data"
FASTQ_R1="$DATA_DIR/SAMPLE_1_1.fastq.gz"
FASTQ_R2="$DATA_DIR/SAMPLE_1_2.fastq.gz"

RG_INFO="@RG\tID:SAMPLE_1\tSM:SAMPLE_1\tPL:ILLUMINA"

ALIGNED_BAM="$PROJECT_DIR/sample1_aligned.bam"
DEDUP_BAM="$PROJECT_DIR/sample1_dedup.bam"
RAW_VCF="$PROJECT_DIR/sample1_raw_variants.vcf"

# --- PHASE 1: DATA PREPARATION ---

echo "Step 1: BWA-MEM Alignment and Samtools Sorting..."
bwa mem -t 16 -R "$RG_INFO" $REF $FASTQ_R1 $FASTQ_R2 | samtools sort -@ 8 -o $ALIGNED_BAM -

echo "Step 2: GATK MarkDuplicates..."
gatk MarkDuplicates \
    -I $ALIGNED_BAM \
    -O $DEDUP_BAM \
    -M $PROJECT_DIR/sample1_marked_dup_metrics.txt \
    --REMOVE_DUPLICATES true

echo "Step 3: Samtools Indexing..."
samtools index $DEDUP_BAM

# --- PHASE 2: VARIANT DISCOVERY ---

echo "Step 4: GATK HaplotypeCaller..."
gatk HaplotypeCaller \
    -R $REF \
    -I $DEDUP_BAM \
    -O $RAW_VCF

echo "Step 5: ANNOVAR Annotation Engine..."
perl $ANNOVAR_DIR/table_annovar.pl $RAW_VCF $ANNOVAR_DIR/humandb/ \
    -buildver hg38 \
    -out $PROJECT_DIR/sample1_annotated \
    -remove -protocol refGene,clinvar,avsnp150,1000g2015aug_all,exac03,dbnsfp47a \
    -operation g,f,f,f,f,f \
    -nastring . \
    -vcfinput

echo "Step 6: Clinical IEI Filtering..."
head -n 1 $PROJECT_DIR/sample1_annotated.hg38_multianno.txt > $PROJECT_DIR/sample1_iei_candidates.txt
grep -F -w -f $PANEL $PROJECT_DIR/sample1_annotated.hg38_multianno.txt >> $PROJECT_DIR/sample1_iei_candidates.txt

echo "Pipeline Complete!"
date