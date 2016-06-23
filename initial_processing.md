# Initial GBS data processing

## Lucas Nell, June 2016

# Combining FASTQ from sequencing facility

Two libraries of 96 samples, using paired-end genotyping by sequencing (GBS) using the
restriction enzyme *ApeKI*. Libraries are from runs 132571 and ... **FILL IN LATER**.

Files from run 132571 (located in
`/Volumes/MW_18TB/NextGen_RawData/Mouse_GBS_CASTxDOM_F2`):
```bash
GBS-2_S1_L001_R1_001.fastq.gz	GBS-2_S1_L003_R1_001.fastq.gz
GBS-2_S1_L001_R2_001.fastq.gz	GBS-2_S1_L003_R2_001.fastq.gz
GBS-2_S1_L002_R1_001.fastq.gz	GBS-2_S1_L004_R1_001.fastq.gz
GBS-2_S1_L002_R2_001.fastq.gz	GBS-2_S1_L004_R2_001.fastq.gz
```

These were combined by read # using the following code:
```bash
cat *_R1_*.fastq.gz > musGBS_132571_1.fastq.gz
cat *_R2_*.fastq.gz > musGBS_132571_2.fastq.gz
```

These resulting files are located in `/lustre1/lan/musGBS/fastq/run_132571`.



# Demultiplex

The above files were then separated by individual using Stacks' `process_radtags`,
resulting in 192 files per run. For run 132571, the code for this would be...
```bash
export runNum=132571

cd /lustre1/lan/musGBS/fastq/run_${runNum}

module load stacks/1.40

export barcodes=barcodes_${runNum}.txt

mkdir -p separated

process_radtags \
    -1 musGBS_${runNum}_1.fastq.gz \
    -2 musGBS_${runNum}_2.fastq.gz \
    -b ${barcodes} \
    -o ./separated/ \
    --inline_inline \
    -e apeKI \
    -r -c -q \
    -i gzfastq &> demultiplex_${runNum}.log
```

The output folder (`./separated/`) will contain a bunch of files ending with
`*.rem.?.fq.gz` (`?` will be `1` or `2`). The following from the
[Stacks manual](http://catchenlab.life.illinois.edu/stacks/manual/#clean)
explains what these files are:

> The `process_radtags` program wants to keep the reads in *phase*, so that the first
> read in the `sample_XXX.1.fq` file is the mate of the first read in the
> `sample_XXX.2.fq` file. Likewise for the second pair of reads being the second
> record in each of the two files and so on. When one read in a pair is discarded due
> to low quality or a missing restriction enzyme cut site, the remaining read can't
> simply be output to the `sample_XXX.1.fq` or `sample_XXX.2.fq` files as it would
> cause the remaining reads to fall out of phase. Instead, this read is considered a
> remainder read and is output into the `sample_XXX.rem.1.fq` file if the paired-end
> was discarded, or the `sample_XXX.rem.2.fq` file if the single-end was discarded.

I relegated these files to a folder called `rem_files`.


# Fix file and read names

File renaming was only done to conform to other fastq files in my directory. This is
unnecessary if you don't mind the *.1.fq.gz *.2.fq.gz naming from Stacks.

Because `process_radtags` adds "_2" to the end of read names in the R2 files and
"_1" at the end of R1 files, I'm going change the R1 files to match the read names in
the R1 files.

```bash
export runNum=132571

cd /lustre1/lan/musGBS/fastq/run_${runNum}/separated

for f in *.fq.gz
do
    # Fix file names
    g=`echo ${f/%.fq.gz/} | sed 's/.1$/_1/g; s/.2$/_2/g'`.fastq.gz
    mv $f $g
    # Rename reads in *_2 files
    if [[ ${g} == *"_2.fastq.gz"* ]]
    then
        gunzip -c ${g} | sed 's/_2$/_1/g' - | gzip > ${g/_2.fastq/_2_tmp.fastq}
        mv ${g/_2.fastq/_2_tmp.fastq} ${g}
    fi
done
```
