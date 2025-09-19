#!/usr/bin/bash

##################################################################
# ensure correct script usage

# activates python3.8 virtual environment
if [[ -d ~/python38venv ]]; then
	source ~/python38venv/bin/activate
fi

if [ $# -eq 0 ]; then
  echo "Usage: $0 followed by directories containing hmmer output domtables from which to extract sequences or domains."
  exit 1
fi

start_=$(date +%s)

path=${0%/*}
folders=$@

# determine input files
domtables=()
i=0
for folder in ${folders[@]}; do
  for domt in $(find $folder -type f -name '*domtable.out'); do
    domtables[$i]=$domt
    let i++
  done
done

# extract domain sequences from database
while true; do
  read -p "Extract full 'gene' sequences [S] or domain sequences only [D]? " seq_type
  read -p "Extract domain sequences from databases: above reporting threshold (E-value<10) [R] or above inclusion threshold (E-value <0.01) [I]? " answer
  read -p "Enter path of directory where sequences will be saved: " output_dir
  echo ""

  # find script path
  path=${0%/*}

  # extract domains over REPORTING threshold
  if [[ $answer == "R" ]]; then
    dom_r=$(python ${path}/extract_domains.py -R -$seq_type -o "$output_dir" "${domtables[@]}")
    echo "${dom_r}"
    # python script to extract fasta sequences from database and save to file
    break

  #extract domains over INCLUSION threshold
  elif [[ $answer == "I" ]]; then
    dom_r=$(python ${path}/extract_domains.py -I -$seq_type -o "$output_dir" "${domtables[@]}")
    echo "${dom_r}"
    break
  fi
done

stop_=$(date +%s)
time=$(($stop_-$start_))
printf "\nTime needed:${time} s\n"
