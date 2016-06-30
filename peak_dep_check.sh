#!/usr/bin/env bash


# ================================================================================
# ================================================================================

#       Parsing options and providing necessary error messages

# ================================================================================
# ================================================================================

# Resetting parameters that are required for this script
unset inFile outFile # (outFile isn't required, but is dependent on inFile)
# Setting defaults for others
export mapq=false
export genomeFile=/lustre1/lan/musGenome/genome/mm9.genome

# Help message for usage
usage()
{
    echo "usage: peak_dep_check.sh -i input_file [-o output_file] [-q] [-g genome_file]"
    echo "    input_file: Input BAM file to check read depth in."
    echo "    output_file: Name of output depth file, NOT including extension. Defaults" \
              "to input filename with '.bed' extension instead of '.bam'."
    echo "    [-q]: Filter by MAPQ >= 20? Takes 'true' or 'false'. Defaults to 'false'."
    echo "    [genome_file]: Tab-delimited text file with chromosome names in col 1, " \
              "length of chromosomes in col 2. " \
              "Defaults to '/lustre1/lan/musGenome/genome/mm9.genome'."
}

# Parsing options
while getopts ":i:o:qg:h" opt
do
    case $opt in
        i)
            export inFile=$OPTARG
            ;;
        o)
            export outFile=$OPTARG
            ;;
        q)
            export mapq=true
            ;;
        g)
            export genomeFile=$OPTARG
            ;;
        h)
            usage
            exit 0
            ;;
        \?)
            echo "ERROR: Invalid argument for -$OPTARG" >&2
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

# Checking that all required parameters are provided
if [[ -z "$inFile" ]]
then
    echo "ERROR: Option -i is required." >&2
    usage
    exit 1
fi

if [[ -z "$outFile" ]]
then
    export outFile=${inFile/%.bam/.bed}
fi



# ================================================================================
# ================================================================================

#

# ================================================================================
# ================================================================================

export PATH=$PATH:/usr/local/apps/bedtools/latest/bin
module load samtools/latest

if [ ${mapq} = "true" ]
then
    samtools view -bh -q 20 ${inFile} | \
        bedtools genomecov -bga -split -ibam stdin -g ${genomeFile} > ${outFile}.bed
else
    bedtools genomecov -bga -split -ibam ${inFile} -g ${genomeFile} > ${outFile}.bed
fi
