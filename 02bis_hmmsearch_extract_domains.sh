#!/usr/bin/bash

##################################################################
# ensure correct script usage
start_=$(date +%s)

# activates python3.8 virtual environment
if [[ -d ~/python38venv ]]; then
	source ~/python38venv/bin/activate
fi

if [ $# -lt 2 ]; then
  echo "Usage: $0 hmm_profiles_directory fasta1 fasta2 fasta3. Use ./path/to/directory_name"
  exit 1
fi

#first argument is directory with hmm profiles
prof_dir="${1}"
# if last character is "/", remove it
if [[ "${prof_dir: -1}" == "/" ]]; then
  prof_dir=${prof_dir:0: -1}
fi
shift 1
#all other arguments are files to search
fastafiles=$@

#check if argument provided is a directory
if ! [[ -d "$prof_dir" ]]; then
  echo "$prof_dir is not a directory. Please provide a single directory.  Use . if current directory or ./directory_name"
  exit 1
elif ! [[ ${prof_dir:0:2} == "./" ]]; then
  if [[ ${prof_dir:0:2} == "/" ]]; then
    prof_dir=".$prof_dir"
  else
    prof_dir="./$prof_dir"
  fi
  echo "Using $prof_dir as path for profile directory"
else
  for file in $fastafiles; do
    first_line=$(head -1 "${file}")
    first_seqlines=$(head -10 "${file}" | grep -v "^>")
    if ! [[ "${first_line}" =~ ^">".* ]]; then
      echo "First line of $file does not start with >. Provide fasta file."
      exit 1
    elif [[ $(echo "${first_seqlines}" | grep -i -e [RNDQEHILKMFPSWYVBZX] | wc -l) -eq 0 ]]; then
      echo "!!!WARNING !!! your database $file seems to contain nucleotides and not amino acids."
    fi
  done
fi

###########################################################
# send stdout and stderr to logfile
scriptname=$0
logfile="$(basename $scriptname .sh)-$(basename $prof_dir)-logfile.log"
i=1
while true; do
  if [[ -f "$logfile" ]]; then
    logfile="${logfile%.log}$i.log"
    let i++
  else
    break
  fi
done

touch $logfile
exec > >(tee "$logfile") 2>&1
printf "\n¨LOGS: stdout and stderr will be saved to $logfile. If file already exists, it will be overwritten.\n\n"

stop_=$(date +%s)
time=$(($stop_-$start_))
printf "\nTime needed:${time} s\n"

##########################################################
#search fastafiles with hmm profiles
start_=$(date +%s)

#success and fail counters
search_success=0
search_fail=0

printf "\n\n"

j=0
domtables=()

for fasta in ${fastafiles[@]}; do
  base=$(basename $fasta)
  hmmsearch_dir="${prof_dir%_hmm_profiles}_${base%%.*}_hmm_search"
  if ! [[ -d $hmmsearch_dir ]]; then
    mkdir $hmmsearch_dir
  fi
  for hmm_profile in $(find $prof_dir -type f); do
    echo "Searching $hmm_profile in $fasta"
    dom="${hmmsearch_dir}/$(basename ${hmm_profile%.*})-${base%%.*}"
    hmmsearch -o ${dom}.out --domtblout ${dom}_domtable.out ${hmm_profile} ${fasta} \
    && printf "Successfully searched ${fasta} with ${hmm_profile}.\n\nOutput saved to ${dom}.out\nTable output saved to ${dom}_table.out\nDomain table output saved to ${dom}_domtable.out\n"

    #catch exit code for hmmbuild and count success and failures
    status_hmmsearch=$?
    if [[ $status_hmmsearch -eq 0 ]]; then
      let search_success++
    else
      let search_fail++
    fi

    # save domtable filenames to array
    domtables[$j]="${dom}_domtable.out"
    let j++

  printf "\n"
  done
printf "\n"
done


#report successes and failures
printf "$search_success successful searches, $search_fail searches failed.\n"

stop_=$(date +%s)
time=$(($stop_-$start_))
printf "\nTime needed:${time} s\n"

################################################
# extract gene or domain sequences from database
while true; do
  read -p "Extract full 'gene' sequences [S] or domain sequences only [D]?" seq_type
  read -p "Extract domain sequences from databases: above reporting threshold (E-value<10) [R] or above inclusion threshold (E-value <0.01) [I]? " answer
  read -p "Enter path of directory where sequences will be saved: " output_dir

  # write what will be done to stdout
  if [[$seq_type == "S" ]]; then
    echo -n "Extracting full sequences "
  elif [[$seq_type == "R" ]]; then
    echo -n "Extracting domains "
  fi
  
  if [[$answer == "I" ]]; then
    echo -n "above inclusion threshold (E-value <0.01)."
  elif [[$answer == "R" ]]; then
    echo -n "above reporting threshold (E-value<10)."
  fi
  
  echo " OUtput will be saved in $output_dir"

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
