#!/usr/bin/bash

##################################################################
# ensure correct script usage
start_=$(date +%s)

# activates python3.8 virtual environment
if [[ -d ~/python38venv ]]; then
	source ~/python38venv/bin/activate
fi

if [ $# -lt 2 ]; then
  echo "Usage: $0 directory_with_alignments proteome_to_search1 proteome2 proteome3. Use . if current directory or ./path/to/directory_name"
  exit 1
fi

#first argument is directory with alignments
al_dir="${1}"
# if last character is "/", remove it
if [[ "${al_dir: -1}" == "/" ]]; then
  al_dir=${al_dir:0: -1}
fi
shift 1
# following arguments are proteomes to search
databases=$@

#check if directory is provided for alignments and if databases are proteomes in fasta format
if ! [[ -d "$al_dir" ]]; then
  echo "$al_dir is not a directory. Please provide a single directory.  Use . if current directory or ./directory_name"
  exit 1
else
  for file in $databases; do
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
logfile="$(basename $scriptname .sh)-$(basename $al_dir)-logfile.log"
touch $logfile
exec > >(tee "$logfile") 2>&1
printf "\n¨LOGS: stdout and stderr will be saved to $logfile. If file already exists, it will be overwritten.\n\n"

stop_=$(date +%s)
time=$(($stop_-$start_))
printf "\nTime needed:${time} s\n"

#########################################################
#build hmm profiles for all files in the directory
#should work with alignments of aa, RNA or DNA sequences (ONLY TESTED FOR AA SO FAR)
start_=$(date +%s)

#success and fails counters
build_success=0
build_fail=0
# create folder to store hmm profiles if it does not exist
# removes _domain_alignments and an eventual / from the name of the directory and adds _hmm_profiles
hmm_profile_dir="${al_dir%_domain_alignments}_hmm_profiles"
if ! [[ -d $hmm_profile_dir ]]; then
  mkdir $hmm_profile_dir
fi
#execute hmmbuild for each file in the directory
for i in $(find $al_dir -type f); do
  echo "Building hmm profile for ${i}"
  hmm_profile_file="${hmm_profile_dir}/$(basename ${i%.*}).hmm"
  hmmbuild "$hmm_profile_file" ${i} && echo "Built hmm profile for ${i} and saved it to $hmm_profile_file"

  #catch exit code for hmmbuild and count success and failures
  status=$?
  if [[ $status -eq 0 ]]; then
    let build_success++
  else
    let build_fail++
  fi

  printf "\n"
done

###################################################

#report successes and failures
printf "$build_success hmm profiles successfully built, $build_fail failed.\n"

stop_=$(date +%s)
time=$(($stop_-$start_))
printf "\nTime needed:${time} s\n"

# ask for user input to continue
while true; do
  read -p "Do you wish to continue with hmmsearch on the databases ${databases[@]}? [Y/N] " answer
  if [[ $answer == "N" ]]; then
    exit 0
  elif [[ $answer == "Y" ]]; then
    break
  fi
done

################################################
#search translated transcriptome/protein database
start_=$(date +%s)

#success and fail counters
search_success=0
search_fail=0

printf "\n\n"

j=0
domtables=()

for database in ${databases[@]}; do
  base_db=$(basename $database)
  hmmsearch_dir="${hmm_profile_dir%_hmm_profiles}_${base_db%%.*}_hmm_search"
  if ! [[ -d $hmmsearch_dir ]]; then
    mkdir $hmmsearch_dir
  fi
  for hmm_profile in $(find $hmm_profile_dir -type f); do
    echo "Searching $hmm_profile in $database"
    dom="${hmmsearch_dir}/$(basename ${hmm_profile%.*})-${base_db%%.*}"
    hmmsearch -o ${dom}.out --domtblout ${dom}_domtable.out ${hmm_profile} ${database} \
    && printf "Successfully searched ${database} with ${hmm_profile}.\n\nOutput saved to ${dom}.out\nTable output saved to ${dom}_table.out\nDomain table output saved to ${dom}_domtable.out\n"

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
  read -p "Extract (domain) sequences from databases: above reporting threshold (E-value<10) [R] or above inclusion threshold (E-value <0.01) [I]? " answer
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
