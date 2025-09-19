#!/usr/bin/bash

##################################################################
# ensure correct script usage

# activates python3.8 virtual environment
if [[ -d ~/python38venv ]]; then
	source ~/python38venv/bin/activate
fi

if [ $# -ne 1 ]; then
  echo "Usage: $0 CD_search_hitdata"
  exit 1
fi

#Extract names and accession numbers for domains from a CD search
hitdata=$1

##################################################################
# check input file

#check if input file is correct
if ! grep -q "#Batch CD-search tool" $hitdata; then
  echo "Input file should be Batch CD-search output." 
  exit 1
fi
if ! grep -q "#datatype.*hits.*" $hitdata; then
  echo "Data type should be hits"
  exit 1
fi

if grep -q "%" $hitdata; then
  echo "WARNING! There is a % symbol in your query names. The output would not be reliable given how the table is parsed.!"
  exit 1
fi

# send stdout and stderr to logfile
scriptname=$0
logfile="$(basename $scriptname .sh)-$(basename $hitdata .txt)-logfile.log"
touch $logfile
exec > >(tee "$logfile") 2>&1
printf "\nLOGS: stdout and stderr will be saved to $logfile. If file already exists, it will be overwritten.\n\n"

#####################################################
# From CDD gene search hitdata, extract domain alignments for each domain

# create file with domain names and accession numbers
output="$(basename $hitdata .txt)_domains_names_ANs.txt"

  # write column headers to file
grep "^Query.*" ${hitdata} | sed 's/\t/%/g' | cut -d % -f 8,9 | sed 's/%/\t/g' > "temporary${output}"

  # write ANs and domain names for hit type: specific. If duplicates, remove them
grep -v "^#" "${hitdata}" | grep "specific" | grep -v "non-specific" | sed 's/\t/%/g' | cut -d % -f 8,9 | sed 's/%%*/\t/g' >> "temporary${output}"
awk '!array[$0]++' temporary${output} > ${output}
rm temporary${output}
echo "Extracted all domain names and accession numbers from ${hitdata} and saved them to ${output} without repetitions:"
printf "\n$(cat $output)\n\n"
echo "To avoid downloading alignments for some domains, manually remove the corresponding line from file $output"

while true; do
  read -p "Write 'Done' when ready to proceed. " answer

  if [ $answer == "Done" ]; then
    break
  fi
done

printf "\n\nAlignments for the following domains will be extracted:\n\n$(cat $output)\n\n"

# Create variable with only accession numbers
acc_nums=$(grep -v "Accession" ${output} | cut -f 1)

# Download alignments for accession numbers in acc_nums and save them in directory_to_make
# python CDD_search.py $acc_nums
script_path=${scriptname%/*}
hitdata_b=$(basename $hitdata)
directory_to_make="${hitdata_b%.txt}_domain_alignments"

  #call python script - CDD API downloading alignments
printf "$(python "${script_path}/CDD_search.py" "${acc_nums}" "$directory_to_make")"
printf "\n"

exit 0
