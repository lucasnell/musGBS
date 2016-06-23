#!/usr/bin/env bash


# Example usage:
# ./initial_processing.sh -r 132571 -i /lustre1/lan/musGBS/fastq/original/run_${RUN}



# ================================================================================
# ================================================================================

#       Parsing options and providing necessary error messages

# ================================================================================
# ================================================================================

# Resetting parameters that are input to this script
unset runNum
unset inDir


# Help message for usage
usage()
{
    echo "usage: initial_processing.sh -r run_number -i input_directory"
    echo "    -r run_number: Run number, used for folder naming."
    echo "    -i input_directory: Directory to find FASTQ files for processing. It's " \
         "assumed that all files in this directory are to be processed."
}

# Parsing options
while getopts ":r:i:h" opt; do
  case $opt in
    r)
      export runNum=$OPTARG
      ;;
    i)
      export inDir=$OPTARG
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
if [[ -z "$runNum" || -z "$inDir" ]]
then
    echo "ERROR: Options -r and -i are required." >&2
    usage
    exit 1
fi



# ================================================================================
# ================================================================================

#       Combine FASTQ files by read #

# ================================================================================
# ================================================================================

export fastqDir=/lustre1/lan/musGBS/fastq/run_${runNum}

mkdir -p ${fastqDir}

cd ${inDir}

cat *_R1_*.fastq.gz > ${fastqDir}/musGBS_${runNum}_1.fastq.gz
cat *_R2_*.fastq.gz > ${fastqDir}/musGBS_${runNum}_2.fastq.gz




# ================================================================================
# ================================================================================

#       Separate by individual using `process_radtags`

# ================================================================================
# ================================================================================

cd ${fastqDir}

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


cd ./separated

mkdir -p rem_files
mv *.rem.* ./rem_files/




# ================================================================================
# ================================================================================

#       Fixing file and read names

# ================================================================================
# ================================================================================


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
