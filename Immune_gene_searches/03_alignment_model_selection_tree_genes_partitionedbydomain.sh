#!/usr/bin/bash

##################################################################
# ensure correct script usage

start_=$(date +%s)

# activates python3.8 virtual environment
if [[ -d ~/python38venv ]]; then
	source ~/python38venv/bin/activate
fi

if [ $# -ne 1 ]; then
  echo "Usage: $0 folder_with_domain_files"
  exit 1
fi

##################################################################
# check input
seq_dir=$1

#check if it is a directory
if ! [[ -d $seq_dir ]]; then
  echo "$seq_dir is not a directory. Please provide a directory"
  exit 1
# if last character is "/", remove it
elif [[ "${seq_dir: -1}" == "/" ]]; then
  seq_dir=${seq_dir:0: -1}
fi

#check if files inside the directory are fasta
for file in $(find $seq_dir -type f); do
  first_line=$(head -1 "${file}")
  first_seqlines=$(head -10 "${file}" | grep -v "^>")
  if ! [[ "${first_line}" =~ ^">".* ]]; then
    echo "First line of $file does not start with >. Provide only fasta files."
    exit 1
  elif [[ $(echo "${first_seqlines}" | grep -i -e [RNDQEHILKMFPSWYVBZX] | wc -l) -eq 0 ]]; then
    echo "!!!WARNING !!! your database $file seems to contain nucleotides and not amino acids."
  fi
done


# send stdout and stderr to logfile
scriptname=$0
logfile="$(basename $scriptname .sh)-$(basename $seq_dir)-logfile.log"
touch $logfile
exec > >(tee "$logfile") 2>&1
printf "\nLOGS: stdout and stderr will be saved to $logfile. If file already exists, it will be overwritten.\n"

stop_=$(date +%s)
time=$(($stop_-$start_))
printf "\n\nTime needed:${time} s\n\n"
################################################
start_=$(date +%s)
# align sequences in each file (domain) with MAFFT G-INS-i (global, iterative refinement with wsp and consistency score)
# save them to output directory
printf "Aligning sequences in $seq_dir with MAFFT\n"
alignments_dir="${seq_dir}_aligned"
success=0
fails=0
mkdir $alignments_dir
for input_file in $(find $seq_dir -type f -name '*.fa*'); do
  output_file="${alignments_dir}/$(basename $input_file)"
  echo "mafft --globalpair --maxiterate 1000 $input_file > $output_file"
  mafft --globalpair --maxiterate 1000 $input_file > $output_file
  status=$?
  if [[ $status -eq 0 ]]; then
    let success++
  else
    let fails++
  fi
done

printf "${success} successful alignments, ${fails} failed."

stop_=$(date +%s)
time=$(($stop_-$start_))
printf "\n\nTime needed:${time} s\n\n"

start_=$(date +%s)

# concatenate all of the alignments with AMAS (relaxed phylip format) + obtain initial partition file (nexus format)
# find script path
echo "Concatenating alignment files in $alignment_dir with AMAS"
path=${0%/*}
tree_dir="${seq_dir}_tree"
if ! [[ -d $tree_dir ]]; then
  mkdir $tree_dir
fi

conc_file="${tree_dir}/$(basename $seq_dir)_concatenated"
amas_out=$(python ${path}/amas/amas/AMAS.py concat -f fasta -d aa -i $(find $alignments_dir -type f) -u phylip -t "${conc_file}.phy" --part-format nexus -p "${conc_file}_partition.txt")
echo "${amas_out}"

status=$?
if [[ $status -eq 0 ]]; then
  echo "Successfully concatenated files and saved them to ${conc_file}.phy. Partition file saved to ${conc_file}_partition.txt"
else
  echo "Something went wrong concantenating the alignment files. Try again"
  exit 1
fi

echo ""

stop_=$(date +%s)
time=$(($stop_-$start_))
printf "\nTime needed:${time} s\n\n"

# search for best partition and best model with iqtree2 (or?)
# separate step or together with previous one: ML tree with iqtree2 ??
start_=$(date +%s)

printf "Running iqtree2\n\n"

# run iqtree2 for concatenated file and partition file from AMAS
# -m TESTMERGE: choose right partition scheme and model resembling PartitionFinder
  # (by just considering the invariable site and Gamma rate heterogeneity)
  # -mset to test only a subset of models to exclude clade-specific models and mitochondrial/plastidial
  # -AIC uses Akaike's information criterion instead of Bayesian -> is this the best option??
# -B 1000: ultrafast bootstrap (resamples the sites within partitions), with 1000 replicates.
  # number yield by ultrafast bootstrap is to interpret as % prob that a clade is correct
# -T AUTO: determine automatically number of CPU threads to be used
iqtree2 -s "${conc_file}.phy" -p "${conc_file}_partition.txt" -m TESTMERGE -mset Dayhoff,DCMut,LG,JTT,WAG,VT,Poisson,PMB,Q.pfam -B 1000 -T AUTO

stop_=$(date +%s)
time=$(($stop_-$start_))
printf "\nTime needed:${time} s\n"
