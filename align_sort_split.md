# Initial GBS data processing

## Lucas Nell, June 2016

# Common objects to all following scripts

All following sections will use one or more of these objects.

```bash
export runNum=132571
export cores=4
export samp=A10_712
```



# Align

Align to mm9 genome using `bowtie2`. The options starting with `--rg` specify aspects
of the read group, which would allow us to skip adding read groups using PICARD.
Thus this would save time if specified in pipelines utilizing GATK program(s).

```bash
cd /lustre1/lan/musGBS/fastq/run_${runNum}/separated

module load bowtie2/latest

export btBuild=/lustre1/lan/musGenome/bt2_mm9/bt2_mm9
export outDir=/lustre1/lan/musGBS/bam_mapped/run_${runNum}

if [ ! -d "${outDir}" ]; then mkdir ${outDir}; fi

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
```



# SAM -> sorted BAM

```bash
module load samtools/latest
export PATH=$PATH:/usr/local/apps/bedtools/2.25.0/bin

cd /lustre1/lan/musGBS/bam_mapped/run_${runNum}

# Convert to BAM and remove SAM
samtools view -bh -@ $(expr ${cores} - 1) -o ${samp}.bam ${samp}.sam
rm ${samp}.sam
# Sort BAM, remove unsorted BAM, index sorted BAM
samtools sort -o ${samp}_sorted.bam -T ${samp}_ss -@ ${cores} \
    ${samp}.bam
mv ${samp}_sorted.bam ${samp}.bam
samtools index -b ${samp}.bam
```



# Create mapped BAM files

Mapped BAM files were those filtered for reads that met one of the following criteria:

- Both the focal read and its mate were mapped
- The focal read mapped to chrX, but its mate was unmapped

These files are used in the aligned-Stacks and the Bolnick pipelines.
I'm still unsure which pipeline is best because the second run of sequencing is not yet
done, so I haven't been able to finish all pipelines with a full data set.

The below code does the filtering, producing a final mapped BAM for this sample and run.

```bash
module load samtools/latest
export PATH=$PATH:/usr/local/apps/bedtools/2.25.0/bin

cd /lustre1/lan/musGBS/bam_mapped/run_${runNum}

# Filter FOR reads that were mapped, as were their mates
samtools view -F 4 -F 8 -bh -@ $(expr ${cores} - 1) -o ${samp}_bothMapped.bam \
    ${samp}.bam
# Filter FOR reads on chrX that had mate unmapped but were mapped themselves
samtools view -F 4 -f 8 -bh -@ $(expr ${cores} - 1) -o ${samp}_Xpartial.bam \
    ${samp}.bam chrX

# Merge them w/o messing with header, and index final BAM
# (No need to re-sort bc merging maintains position-sorting)
samtools merge -cp ${samp}_mapped.bam \
    ${samp}_bothMapped.bam ${samp}_Xpartial.bam
samtools index -b ${samp}_mapped.bam
# Remove temporary files
rm ${samp}_bothMapped.bam ${samp}_Xpartial.bam
```



# Create unmapped FASTQ files

The unmapped FASTQ files contain reads that match one of the following criteria:

- The focal read and their mates were both unmapped
- The focal read or its mate mapped to chrX, but one was unmapped

These files are used to create de novo stacks in Stacks, which we can then use to create
genetic maps to find stacks that show linkage with stacks on the X chromosome. This way,
we can find stacks on the X chromosome that are in areas not yet on mouse genome mm9.

To extract these reads, I had to filter BAM files using `samtools`, sort filtered BAMs by
name, then merge all filtered BAMs to one BAM. I then converted this merged BAM to FASTQ
files, separated by read #, using `bedtools bamtofastq`. I could then `gzip` these
FASTQ files to reduce file sizes.

Filtering code:

```bash
module load samtools/latest
export PATH=$PATH:/usr/local/apps/bedtools/2.25.0/bin

cd /lustre1/lan/musGBS/bam_mapped/run_${runNum}

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
```


# Moving files

The below code moves all files into their more permanent location on the cluster, first
making sure that all output locations exist.

```bash
for f in "logs" "../bam/run_${runNum}" "../fastq_unmapped/run_${runNum}"
do
    mkdir -p "$f"
done

mv *.log ./logs/
mv ${samp}.bam ${samp}.bam.bai ../bam/run_${runNum}/
mv ${samp}.fastq.gz ../fastq_unmapped/run_${runNum}/
```
