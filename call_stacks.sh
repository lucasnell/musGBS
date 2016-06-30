#!/usr/bin/env bash



# THIS FILE HASN'T BEEN EDITED AT ALL, AND THE BELOW SECTION IS FROM make_stacks.sh




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
