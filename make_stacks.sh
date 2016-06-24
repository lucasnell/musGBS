#!/usr/bin/env bash


# ================================================================================
# ================================================================================

#       Parsing options and providing necessary error messages

# ================================================================================
# ================================================================================

# Resetting parameters that are input to this script
unset outPath inFile mapped


# Help message for usage
usage()
{
    echo "usage: make_stacks.sh -i input_file -o output_directory -m mapped"
    echo "    input_file: Input BAM or gzipped FASTQ file."
    echo "    output_directory: Folder in which to put all stacks files."
    echo "    mapped: Are these mapped (BAM) files? Takes 'true' or 'false'."
}

# Parsing options
while getopts ":r:i:h" opt; do
  case $opt in
    o)
      export outPath=$OPTARG
      ;;
    i)
      export inFile=$OPTARG
      ;;
    m)
      export mapped=$OPTARG
      if [[ "$mapped" != "true" && "$mapped" != "false" ]]
      then
          echo "ERROR: Option -m only takes 'true' or 'false'." >&2
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

# Checking that all options are provided
if [[ -z "$outPath" || -z "$inFile" || -z "$mapped" ]]
then
    echo "ERROR: Options -o, -i, and -m are required." >&2
    usage
    exit 1
fi




# ================================================================================
# ================================================================================

#       Create stacks

# ================================================================================
# ================================================================================

module load stacks/1.40

mkdir -p ${outPath}/logs

# Get integer id from numeric chars after 1st underscore
# (e.g., file H7_2013_132571.bam would be 2013)
tmp=(`echo ${inFile} | sed 's/\.bam//g; s/\.fastq//g; s/\.gz//g; s/_/\ /g'`)
export int_id=${tmp[1]}

if [ "${mapped}" = "true" ]
then
    pstacks \
        -t bam \
        -f ${inFile} \
        -o ${outPath} \
        -i ${int_id} \
        &> ${outPath}/logs/${inFile/%.bam/}.log
else
    ustacks \
        -t gzfastq \
        -f ${inFile} \
        -o ${outPath} \
        -i ${int_id} \
        &> ${outPath}/logs/${inFile/%.fastq.gz/}.log
fi
