#!/usr/bin/env bash


# ================================================================================
# ================================================================================

#       Parsing input parameters and providing necessary error messages

# ================================================================================
# ================================================================================

# Resetting parameters that are input to this script
unset runNum
unset samp
unset cores

# Help message for usage
usage()
{
    echo "usage: align_sort_split.sh -r run_number -s sample_name -c available_cores"
    echo "    -r run_number: Run number, used for folder naming."
    echo "    -s sample_name: Sample name, used for file naming."
    echo "    -c available_cores: Number of cores available for multithreaded " \
         "processing in SAMtools and bowtie2. "
}

# Parsing input parameters
while getopts ":r:s:c:h" opt; do
  case $opt in
    r)
      export runNum=$OPTARG
      ;;
    s)
      export samp=$OPTARG
      ;;
    c)
      export cores=$OPTARG
      ;;
    h)
      usage
      exit 0
      ;;
    \?)
      echo "ERROR: Invalid option for -$OPTARG" >&2
      usage
      exit 1
      ;;
    :)
      echo "ERROR: Option -$OPTARG requires an argument." >&2
      usage
      exit 1
      ;;
  esac
done

# Checking that all options are provided
if [[ -z "$runNum" || -z "$samp" || -z "$cores" ]]
then
    echo "ERROR: Options -r, -s, and -c are required." >&2
    usage
    exit 1
fi






# ================================================================================
# ================================================================================

#       Alignment to mouse genome, mm9

# ================================================================================
# ================================================================================

cd /lustre1/lan/musGBS/fastq/run_${runNum}/separated

module load bowtie2/latest

export btBuild=/lustre1/lan/musGenome/bt2_mm9/bt2_mm9
export outDir=/lustre1/lan/musGBS/bam_mapped/run_${runNum}

export read1_list=`ls -m ${samp}*_1.fastq.gz | tr -d ' \n'`
export read2_list=`ls -m ${samp}*_2.fastq.gz | tr -d ' \n'`

bowtie2 -p ${cores} -x ${btBuild} \
    -1 ${read1_list} \
    -2 ${read2_list} \
    --rg-id ${samp}_${runNum} \
    --rg SM:${samp} \
    --rg PL:ILLUMINA \
    --rg LB:${runNum} \
    -S ${outDir}/${samp}.sam &> ${outDir}/align_${samp}.log



# ================================================================================
# ================================================================================

#       Info for SAM flags

# ================================================================================
# ================================================================================

# https://broadinstitute.github.io/picard/explain-flags.html
# Flag        Description
# 1           the read is paired in sequencing
# 2           the read is mapped in a proper pair
# 4           the query sequence itself is unmapped
# 8           the mate is unmapped
# 16          read reverse strand
# 32          mate reverse strand
# 64          the read is the first read in a pair
# 128         the read is the second read in a pair
# 256         the alignment is not primary
# 512         the read fails platform/vendor quality checks
# 1024        the read is either a PCR or an optical duplicate
# 2048        supplementary alignment





# ================================================================================
# ================================================================================

#       Convert SAM to sorted BAM

# ================================================================================
# ================================================================================

module load samtools/latest
export PATH=$PATH:/usr/local/apps/bedtools/2.25.0/bin

cd ${outDir}

# Convert to BAM and remove SAM
samtools view -bh -@ $(expr ${cores} - 1) -o ${samp}.bam ${samp}.sam
rm ${samp}.sam
# Sort BAM, remove unsorted BAM, index sorted BAM
samtools sort -o ${samp}_sorted.bam -T ${samp}_ss -@ ${cores} \
    ${samp}.bam
mv ${samp}_sorted.bam ${samp}.bam
samtools index -b ${samp}.bam


# ================================================================================
# ================================================================================

#       Create mapped BAM files

# ================================================================================
# ================================================================================

# Filter FOR reads that were mapped, as were their mates
samtools view -F 4 -F 8 -bh -@ $(expr ${cores} - 1) -o ${samp}_bothMapped.bam \
    ${samp}.bam
# Filter FOR reads on chrX that had mate unmapped but were mapped themselves
samtools view -F 4 -f 8 -bh -@ $(expr ${cores} - 1) -o ${samp}_Xpartial.bam \
    ${samp}.bam chrX

# Merge them w/o messing with header and index final BAM
# (No need to re-sort bc merging maintains position-sorting)
samtools merge -cp ${samp}_mapped.bam \
    ${samp}_bothMapped.bam ${samp}_Xpartial.bam
samtools index -b ${samp}_mapped.bam
# Remove temporary files
rm ${samp}_bothMapped.bam ${samp}_Xpartial.bam





# ================================================================================
# ================================================================================

#       Create unmapped FASTQ files

# ================================================================================
# ================================================================================

# Filter for reads that were unmapped AND their mates were unmapped
samtools view -f 12 -bh -@ $(expr ${cores} - 1) -o ${samp}_allUnmapped.bam \
    ${samp}.bam
# Filter for reads on chrX that had mate mapped but were unmapped themselves
samtools view -f 4 -F 8 -bh -@ $(expr ${cores} - 1) -o ${samp}_Xpart1.bam \
    ${samp}.bam chrX
# Filter for reads on chrX that were mapped themselves but had mate unmapped
samtools view -f 8 -F 4 -bh -@ $(expr ${cores} - 1) -o ${samp}_Xpart2.bam \
    ${samp}.bam chrX


# For all BAMs above, sort by name
for f in ${samp}_allUnmapped.bam ${samp}_Xpart?.bam
do
    samtools sort -o ${f/.bam/_sorted.bam} -n -T ${f/.bam/_s} -@ ${cores} ${f}
    mv ${f/.bam/_sorted.bam} ${f}
done

# Merge them w/o messing with header (maintains name-sorting)
samtools merge -ncp ${samp}_unmapped.bam \
    ${samp}_allUnmapped.bam ${samp}_Xpart?.bam

# Convert BAM file to fastq files
bedtools bamtofastq -i ${samp}_unmapped.bam \
                      -fq ${samp}_1.fastq \
                      -fq2 ${samp}_2.fastq

# Combine them, sort by name, and gzip
cat ${samp}_?.fastq | \
    paste - - - - | \
    sort -k1,1 -t " " | \
    tr "\t" "\n" | \
    gzip > ${samp}.fastq.gz

# Remove temporary files
rm ${samp}_allUnmapped.bam ${samp}_Xpart?.bam \
    ${samp}_unmapped.bam ${samp}_?.fastq





# ================================================================================
# ================================================================================

#       Moving files

# ================================================================================
# ================================================================================

# Move files other than mapped BAMs out of here, first checking that all output
# directories exist

for f in "logs" "../bam/run_${runNum}" "../fastq_unmapped/run_${runNum}"
do
    mkdir -p "$f"
done

mv *.log ./logs/
mv ${samp}.bam ${samp}.bam.bai ../bam/run_${runNum}/
mv ${samp}.fastq.gz ../fastq_unmapped/run_${runNum}/
