#!/usr/bin/env bash


# ================================================================================
# ================================================================================

#       Parsing options and providing necessary error messages

# ================================================================================
# ================================================================================

# Resetting parameters that are input to this script
unset inFiles outFile fileFormat


# Help message for usage
usage()
{
    echo "usage: merge_runs.sh -i input_files -o output_file -f file_format"
    echo "    input_files: List of input BAM or gzipped FASTQ files to merge into 1, " \
              "comma-separated (no spaces)."
    echo "    output_file: File name of merged file."
    echo "    file_format: File format for input files (should be the same for all). " \
              "Takes 'bam' or 'gzfastq'."
}

# Parsing options
while getopts ":i:o:f:h" opt
do
    case $opt in
        i)
            export inFiles=$OPTARG
            ;;
        o)
            export outFile=$OPTARG
            ;;
        f)
            export fileFormat=$OPTARG
            if [[ "$fileFormat" != "bam" && "$fileFormat" != "gzfastq" ]]
            then
                echo "ERROR: Option -f only takes 'bam' or 'gzfastq'." >&2
                usage
                exit 1
            fi
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
if [[ -z "$inFiles" || -z "$outFile" || -z "$fileFormat" ]]
then
    echo "ERROR: Options -i, -o, and -f are required." >&2
    usage
    exit 1
fi


# ------------
# Checking that filename extensions match what's provided
# ------------
# Setting parameters for checks
if [ "${fileFormat}" = "bam" ]
then
    i=1
    ext="bam"
else
    i=2
    ext="fastq gz"
fi

# For input files
for file in ${inFiles//,/\ }
do
    g=(${file//./\ })
    if [ "${g[@]:(-i)}" != "${ext}" ]
    then
        echo "ERROR: Input filename extension(s) do not match the file format " \
             "provided." >&2
        exit 1
    fi
done

# For output file
g=(${outFile//./\ })
if [ "${g[@]:(-i)}" != "${ext}" ]
then
    echo "ERROR: Output filename extension does not match the file format " \
         "provided." >&2
    exit 1
fi



# ================================================================================
# ================================================================================

#       Create files

# ================================================================================
# ================================================================================

# Splitting $outFile into containing folder and filename
if [[ "${outFile}" == *"/"* ]]
then
    tmp=(${outFile//\//\ })
    export outFileName=${tmp[${#tmp[@]}-1]}
    unset tmp[${#tmp[@]}-1]
    export outPath=/`echo ${tmp[@]} | tr ' ' '/'`
else
    export outPath='.'
    export outFileName=${outFile}
fi



mkdir -p ${outPath}/logs


if [ "${fileFormat}" = "bam" ]
then
    export javMem=2
    module load java/latest
    module load samtools/latest
    module load picard/2.4.1

    java -Xmx${javMem}g \
        -jar picard.jar \
        MergeSamFiles \
        # Below will create something like `I=one.bam I=two.bam ... I=N.bam`
        I=${inFiles/,/\ I=} \
        O=${outFile} \
        &> ${outPath}/logs/${outFileName/%.bam/}.log

    samtools index -b ${outFile} &>> ${outPath}/logs/${outFileName/%.bam/}.log

else
    # Combine *.fastq.gz, sort by name, and re-gzip
    gunzip -c ${inFiles//,/\ } | \
        paste - - - - | \
        sort -k1,1 -t " " | \
        tr "\t" "\n" | \
        gzip \
        1> ${outFile} \
        2> ${outPath}/logs/${outFileName/%.fastq.gz/}.log
fi
