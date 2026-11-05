#Quality control for 16S and ITS amplicon samples. Includes adaptor/quality filtering with trimmomatic,
#as well as phix (phix174 reference) and human contamination (hg19 reference) removal with bbduk. 

#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4
#SBATCH --account=def-nagissa
#SBATCH --job-name=Illumina_read_cleaning
#SBATCH --array=1-23%10

#NOTE: Threads were kept consistent at 4 to match quality control with metatranscriptome samples. Threading can alter the initial organization and subsequent read recovery
#during QC. To improve reproducibility of our results, the thread count was kept consistent for the entire QC pipeline across both metagenome and metatranscriptome samples.

#Some flags to stop the whole pipeline if an error occurs
set -euo pipefail
set -x
threads=$SLURM_CPUS_PER_TASK

#Set sample_list to loop through
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 16S_samplenames.txt)

#Load modules in Nibi cluster (Compute Canada)
module load StdEnv/2020
module load trimmomatic
module load bbmap
module load fastqc

echo "${sample}"

#Making all the directories + define location of adaptor file for trimmomatic
indir=/project/6105702/CMO_2024_Metagenome/goordial/second_round_of_sequencing/goordial_air_snowfall_second_sequencing_round/16S
outdir=/project/6105702/CMO_2024_Metagenome/goordial/second_round_of_sequencing/goordial_air_snowfall_second_sequencing_round/16S/clean
mkdir -p ${outdir}/clean_fastq/readcounts
mkdir -p ${outdir}/clean_fastq/normalized
mkdir -p ${outdir}/stats

adaptor=/project/6105702/2025_cmo_metagenome/general_scripts/NexteraPE-PE.fa

# Create output folder for logs
logdir=${outdir}/illumina_clean_logs
mkdir -p ${logdir}

# Initialize read count log
readcount_dir=${outdir}/clean_fastq/readcounts
mkdir -p ${readcount_dir}

readcount_log=${readcount_dir}/readcounts_${sample}.txt
echo "Step,Forward_reads,Reverse_reads" > $readcount_log

#Set the “current” clean files to the initial clumpify output
r1=$(ls ${indir}/*${sample}*R1*.fastq.gz)
r2=$(ls ${indir}/*${sample}*R2*.fastq.gz)
current_R1=${outdir}/clean_fastq/${sample}_clean_R1.fastq.gz
current_R2=${outdir}/clean_fastq/${sample}_clean_R2.fastq.gz


# Log starting read count per file
echo "raw,$(zcat $r1 | wc -l | awk '{print $1/4}'),$(zcat $r2 | wc -l | awk '{print $1/4}')" >> $readcount_log

#Trimmomatic
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $SLURM_CPUS_PER_TASK \
  ${r1} ${r2} \
  ${current_R1} ${current_R1}_unpaired \
  ${current_R2} ${current_R2}_unpaired \
  ILLUMINACLIP:${adaptor}:2:30:10:8:True SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36

# Log read counts after fastp
echo "fastp,$(zcat $current_R1 | wc -l | awk '{print $1/4}'),$(zcat $current_R2 | wc -l | awk '{print $1/4}')" >> $readcount_log

# Remove phix contamination with bbduk
tmp_R1=${outdir}/clean_fastq/${sample}_tmp_R1.fastq.gz
tmp_R2=${outdir}/clean_fastq/${sample}_tmp_R2.fastq.gz

mkdir -p ${outdir}/stats/phixstats

bbduk.sh \
in1=$current_R1 \
in2=$current_R2 \
out=$tmp_R1 \
out2=$tmp_R2 \
ref=$(dirname $(which bbduk.sh))/resources/phix174_ill.ref.fa.gz \
k=31 \
hdist=1 \
stats=${outdir}/stats/phixstats/${sample}_phix_stats.txt \
ow=t \
threads=$SLURM_CPUS_PER_TASK \
-Xmx40g

mv $tmp_R1 $current_R1
mv $tmp_R2 $current_R2

# Log read counts after phix removal
echo "phix_removal,$(zcat $current_R1 | wc -l | awk '{print $1/4}'),$(zcat $current_R2 | wc -l | awk '{print $1/4}')" >> $readcount_log

# Remove human contamination with bbduk
tmp_R1=${outdir}/clean_fastq/${sample}_tmp_R1.fastq.gz
tmp_R2=${outdir}/clean_fastq/${sample}_tmp_R2.fastq.gz

mkdir -p ${outdir}/human_contam

removehuman.sh \
path=/project/6105702/CMO_2024_Metagenome/reference_contam/hg19_index/ \
in1=$current_R1 \
in2=$current_R2 \
outu1=$tmp_R1 \
outu2=$tmp_R2 \
outm1=${outdir}/human_contam/${sample}_humancon_R1.fastq.gz \
outm2=${outdir}/human_contam/${sample}_humancon_R2.fastq.gz \
ow=t \
threads=$SLURM_CPUS_PER_TASK \
-Xmx40g

mv $tmp_R1 $current_R1
mv $tmp_R2 $current_R2

# Log read counts after human removal
echo "human_removal,$(zcat $current_R1 | wc -l | awk '{print $1/4}'),$(zcat $current_R2 | wc -l | awk '{print $1/4}')" >> $readcount_log" >> ${readcount_log}
