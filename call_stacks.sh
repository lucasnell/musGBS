#!/usr/bin/env bash




# ================================================================================
# ================================================================================

#       Parsing options and providing necessary error messages

# ================================================================================
# ================================================================================

# Resetting parameters that are input to this script
unset inStacks catFile outPath
export threads=1


# Help message for usage
usage()
{
    echo "usage: call_stacks.sh -i input_stacks -c catalog_file -o output_directory " \
              "[-t threads]"
    echo "    input_stacks: Input stacks (filenames without `.tags.tsv`, " \
              "`.snps.tsv`, or `.alleles.tsv` extensions); multiples should be " \
              "separated by commas."
    echo "    catalog_file: Stacks catalog to call from; this will be the filename " \
              "excluding `.catalog` and the extensions also excluded in `input_stacks`." \
              " It should end in an underscored followed by an integer."
    echo "    output_directory: Folder in which to put all stacks files."
    echo "    threads: Number of threads to use. Defaults to 1."
}

# Parsing options
while getopts ":i:c:o:t:h" opt; do
  case $opt in
    i)
      export inStacks=$OPTARG
      ;;
    c)
      export catFile=$OPTARG
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
if [[ -z "$inStacks" || -z "$catFile" || -z "$outPath" ]]
then
    echo "ERROR: Options -i, -c, and -o are required." >&2
    usage
    exit 1
fi




# ================================================================================
# ================================================================================

#       Running sstacks

# ================================================================================
# ================================================================================


stacks_tmp=(${inStacks//,/\ })
export inStacks=${stacks_tmp[@]/#/' -s '}
# Getting the batch id from the catalog name
tmp=(${catFile//_/\ })
export batch_id=${tmp[@]:(-1)}

module load stacks/1.40

sstacks \
    -b ${batch_id} \
    -c ${catFile} \
    ${inStacks} \
    -o ${outPath} \
    -p ${threads} &> ${outPath}/search.log
