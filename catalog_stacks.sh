#!/usr/bin/env bash



# THIS FILE HASN'T BEEN EDITED AT ALL, AND THE BELOW SECTION IS FROM make_stacks.sh




# ================================================================================
# ================================================================================

#       Parsing options and providing necessary error messages

# ================================================================================
# ================================================================================

# Resetting parameters that are input to this script
unset inStacks batchNum outPath
export threads=1


# Help message for usage
usage()
{
    echo "usage: catalog_stacks.sh -i input_stacks -b batch_number -o output_directory " \
              "[-t threads]"
    echo "    input_stacks: Input stacks (filenames without `.tags.tsv`, " \
              "`.snps.tsv`, or `.alleles.tsv` extensions); multiples should be " \
              "separated by commas."
    echo "    batch_number: The number for this batch. Must be an integer."
    echo "    output_directory: Folder in which to put all stacks files."
    echo "    threads: Number of threads to use. Defaults to 1."
}

# Parsing options
while getopts ":i:b:o:t:h" opt; do
  case $opt in
    i)
      export inStacks=$OPTARG
      ;;
    b)
      export batchNum=$OPTARG
      ;;
    o)
      export outPath=$OPTARG
      ;;
    t)
      export threads=$OPTARG
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

# Checking that all options are provided
if [[ -z "$inStacks" || -z "$batchNum" || -z "$outPath" ]]
then
    echo "ERROR: Options -i, -b, and -o are required." >&2
    usage
    exit 1
fi




# ================================================================================
# ================================================================================

#       Running cstacks

# ================================================================================
# ================================================================================


stacks_tmp=(${inStacks//,/\ })
export inStacks=${stacks_tmp[@]/#/' -s '}

module load stacks/1.40

cstacks \
    -b ${batchNum} \
    ${inStacks} \
    -o ${outPath} \
    -p ${threads} &> ${outPath}/catalog.log
