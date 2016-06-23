# Bioinformatics Pipeline, edited from Dan Bolnick, 2015

## Lucas Nell, June 2016

# Data

See `initial_processing.md` for how I made the mapped BAM files.

# Call SNPs using `samtools`
I use `samtools mpileup` piped to `bcftools call` for SNP calling.
The Bolnick pipeline was a little outdated, so I had to use
[this guide](https://github.com/samtools/bcftools/wiki/HOWTOs)
for an update on calling SNPs using this method.

```bash
module load samtools/latest
module load vcftools/0.1.12b

export ref=/lustre1/lan/musGenome/fasta/mm9.fa
cd /lustre1/lan/musGBS/bam_mapped
export bamFiles=`ls *.bam`

samtools mpileup -t DP,SP -uf ${ref} ${bamFiles} | bcftools call -mv -V indels \
    -O v -o ../vcf/musGBS.vcf &> ../vcf/musGBS_call.log
```



# Filter using `vcftools` and output to 012 files

### meanDP filter
I then used `vcftools` to filter for "...sites with mean depth values (over all included
individuals) greater than or equal to the '--min-meanDP' value and less than or equal
to the '--max-meanDP' value" ([`vcftools`
manual](https://vcftools.github.io/man_latest.html)).
Per Dan Bolnick's pipeline, I chose a max of 75 and min of 1.


### More filtering and 012 files
I next filtered for sites with depths ≤ 100 and ≥ 4, 6, and 8, then for
completeness of 0.5, 0.8, and 0.9 to create 9 new VCF files in total.

Then I made 012 files from these filtered VCF files.
Below is a description of 012 files.

> This option outputs the genotypes as a large matrix. Three files are produced. The
> first, with suffix ".012", contains the genotypes of each individual on a separate
> line. Genotypes are represented as 0, 1 and 2, where the number represent that number
> of non-reference alleles. Missing genotypes are represented by -1. The second file,
> with suffix ".012.indv" details the individuals included in the main file. The third
> file, with suffix ".012.pos" details the site locations included in the main file.

(from the [`vcftools` manual](https://vcftools.github.io/man_latest.html))

*Note:* I included the lines starting with `echo` in the code below to make sure the log
file has identification for which stderr/stdout goes with which filter.


### Code

(It is assumed that you've already loaded necessary modules from using the code above.)

```bash
cd /lustre1/lan/musGBS/vcf

# -----------
# meanDP filter
# -----------
vcftools --vcf musGBS.vcf --out musGBS_1-75 \
    --min-meanDP 1 --max-meanDP 75 --recode &> musGBS_1-75.log

# -----------
# More filtering and 012 files
# -----------
export logFile=musGBS_1-75_filt012.log

echo -e "=============   Started `date "+%H:%M:%S %d-%h-%Y"`   =============\n\n\n" \
    > ${logFile}

for dp in 4 6 8
do
    echo -e "~~~~~~~~~~~   --minDP ${dp}   ~~~~~~~~~~~\n" >> ${logFile}
    vcftools --vcf musGBS_1-75.recode.vcf --recode --maxDP 100 --minDP ${dp} \
        --out musGBS_1-75_DP${dp} &>> ${logFile}
    echo -e "\n\n\n" >> ${logFile}
    for comp in 0.5 0.8 0.9
    do
        echo -e "~~~~~~~~~~~   --minDP ${dp} --max-missing ${comp}   ~~~~~~~~~~~\n" \
            >> ${logFile}
        vcftools --vcf musGBS_1-75_DP${dp}.recode.vcf --recode \
            --out musGBS_1-75_DP${dp}_comp${comp} --max-missing ${comp} &>> ${logFile}
        echo -e "~~~~~~~~~~~   --minDP ${dp} --max-missing ${comp} --012   ~~~~~~~~~~~\n" \
            >> ${logFile}
        vcftools --vcf musGBS_1-75_DP${dp}_comp${comp}.recode.vcf --012 \
            --out musGBS_1-75_DP${dp}_comp${comp} &>> ${logFile}
    done
done

echo -e "\n\n\n=============   Ended `date "+%H:%M:%S %d-%h-%Y"`   =============" \
    >> ${logFile}
```




# Combined code

A full submission script for this pipeline:

```bash
#PBS -S /bin/bash
#PBS -q batch
#PBS -N musGBS_Bpipe
#PBS -l nodes=1:ppn=1:mwnode
#PBS -l mem=2gb
#PBS -l walltime=48:00:00
#PBS -M lucnell@gmail.com
#PBS -m ae
#PBS -j oe

module load samtools/latest
module load vcftools/0.1.12b

# -----------
# Call SNPs
# -----------
export ref=/lustre1/lan/musGenome/fasta/mm9.fa
cd /lustre1/lan/musGBS/bam_mapped
export bamFiles=`ls *.bam`

samtools mpileup -t DP,SP -uf ${ref} ${bamFiles} | bcftools call -mv -V indels \
    -O v -o ../vcf/musGBS.vcf &> ../vcf/musGBS_call.log


cd /lustre1/lan/musGBS/vcf

# -----------
# meanDP filter
# -----------
vcftools --vcf musGBS.vcf --out musGBS_1-75 \
    --min-meanDP 1 --max-meanDP 75 --recode &> musGBS_1-75.log

# -----------
# More filtering and 012 files
# -----------
export logFile=musGBS_1-75_filt012.log

echo -e "=============   Started `date "+%H:%M:%S %d-%h-%Y"`   =============\n\n\n" \
    > ${logFile}

for dp in 4 6 8
do
    echo -e "~~~~~~~~~~~   --minDP ${dp}   ~~~~~~~~~~~\n" >> ${logFile}
    vcftools --vcf musGBS_1-75.recode.vcf --recode --maxDP 100 --minDP ${dp} \
        --out musGBS_1-75_DP${dp} &>> ${logFile}
    echo -e "\n\n\n" >> ${logFile}
    for comp in 0.5 0.8 0.9
    do
        echo -e "~~~~~~~~~~~   --minDP ${dp} --max-missing ${comp}   ~~~~~~~~~~~\n" \
            >> ${logFile}
        vcftools --vcf musGBS_1-75_DP${dp}.recode.vcf --recode \
            --out musGBS_1-75_DP${dp}_comp${comp} --max-missing ${comp} &>> ${logFile}
        echo -e "~~~~~~~~~~~   --minDP ${dp} --max-missing ${comp} --012   ~~~~~~~~~~~\n" \
            >> ${logFile}
        vcftools --vcf musGBS_1-75_DP${dp}_comp${comp}.recode.vcf --012 \
            --out musGBS_1-75_DP${dp}_comp${comp} &>> ${logFile}
    done
done

echo -e "\n\n\n=============   Ended `date "+%H:%M:%S %d-%h-%Y"`   =============" \
    >> ${logFile}
```
