import sys
import pandas as pd
import re
import os
import argparse
import time
from pprintpp import pprint

# constants: evalue thresholds
# default from hmmer - both for single domain and full sequence
INCL_THR=float(0.01)
REP_THR=float(10.0)

############################################################################
_start_=time.time()

# parse input arguments and options
parser = argparse.ArgumentParser(description="Extract domains from hmmer domtable files")
# R and I are mutually exclusive (above Reporting or Inclusion threshold)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-R", action='store_true')
group.add_argument("-I", action='store_true')
# D and S are mutually exclusive
group_type = parser.add_mutually_exclusive_group(required=True)
group_type.add_argument("-D", action='store_true')
group_type.add_argument("-S", action='store_true')
# -o to specify output name
parser.add_argument("-o", action='store', dest='main_output_dir')
# all the remanining arguments go in the domtable_list
parser.add_argument('domtable_list',nargs='*')
#parse_args imports and uses sys
args = parser.parse_args()

main_output_dir = args.main_output_dir

############################################################################
# Extract information from domtables
# initialize success counter
success_n=0

# initialize dictionary with keys and lists as values to save info from domtable to dataframe
hmmsearch_dict={"domain_query":[],"seq_hit":[],"database":[],"e-value_seq":[],"e-value_dom":[],"ali_from":[],"ali_to":[],"env_from":[],"env_to":[]}
for file in args.domtable_list:

  # intialize variable to check if only one database per file
  already_database=False

  with open(file) as f:
    for line in f:

      # check info lines
      if re.match("#",line):

        # if program is not hmmsearch, exit with error status
        if re.match("# Program:",line):
          program=line.split(":",1)[1].strip()
          if program != "hmmsearch":
            print(f"ERROR: Program used should be hmmsearch, no {program}. Exiting python script.")
            print(f"Info extraction successful for {success_n} files, missing {len(args.domtable_list)-success_n} files")
            sys.exit(1)

        # if current directory is not the one where hmmsearch was launched, exit with error status
        if re.match("# Current dir:",line):
          curr_dir=line.split(":",1)[1].strip()
          if os.getcwd() != curr_dir:
            print(f"ERROR: Current directory to extract from {file} should be {curr_dir}! Exiting python program")
            print(f"Info extraction successful for {success_n} files, missing {len(args.domtable_list)-success_n} files")
            sys.exit(1)

        # determine database for file
        if re.match("# Target file:",line):
          if already_database == True:
            print(f"ERROR: multiple target files found in {file}. Exiting python program")
            print(f"Info extraction successful for {success_n} files, missing {len(domtable_list)-success_n} files")
            sys.exit(1)
          else:
            database=line.split(":",1)[1].strip()
            already_database=True

    #go back to start of the file
    f.seek(0,0)
    # extract info from domtable lines (not starting with #)
    for line in [ l for l in f if re.match('#',l) == None ]:
      info=line.strip().split()
      hmmsearch_dict["domain_query"].append(info[3])
      hmmsearch_dict["seq_hit"].append(info[0])
      hmmsearch_dict["database"].append(database)
      hmmsearch_dict["e-value_seq"].append(float(info[6]))
      hmmsearch_dict["e-value_dom"].append(float(info[12])) #i-Evalue
      hmmsearch_dict["ali_from"].append(int(info[17])-1) # -1 to have pythonic coordinates (start from 0)
      hmmsearch_dict["ali_to"].append(int(info[18])-1) # -1 to have pythonic coordinates (start from 0) 
      hmmsearch_dict["env_from"].append(int(info[19])-1) # -1 to have pythonic coordinates (start from 0)
      hmmsearch_dict["env_to"].append(int(info[20])-1) # -1 to have pythonic coordinates (start from 0) 

  success_n += 1


length_dict_list=len(hmmsearch_dict["domain_query"])
for key in hmmsearch_dict.keys():
  if len(hmmsearch_dict[key]) != length_dict_list:
    sys.exit(f"Could not match sequence information reliably. List of query seq contains {length_firs_list} domains, while {hmmsearch_dict[key]} contains {len(hmmsearch_dict[key])}")

print(f"Successfully extracted information from {success_n} files out of {len(args.domtable_list)}")

_end_=time.time()
time_used=_end_ - _start_
print(f"\nTime needed: {time_used} s\n")


#############################################################
# determine sequences to extract (above threshold - Reporting or inclusion based on preference)
_start_=time.time()

#if not already present, mkdir in current folder ( . ) containing folders with sequences
if os.path.isdir(main_output_dir) == False:
  os.mkdir(main_output_dir)

#create new dictionaries with only sequences above threshold
# above inclusion or above reporting threshold
if args.I == True and args.R == False:
  print("Selecting sequences above inclusion threshold")
  #set E-value threshold to inclusion threshold
  E_thr=INCL_THR
elif args.I == False and args.R == True:
  print("Selecting sequences above reporting threshold")
  #set E-value threshold to reporting threshold
  E_thr=REP_THR
else:
  sys.exit("Evalue threshold option error")

###############################
# if domains to extract
if args.D == True and args.S == False:
# create a dictionary to store only query name, hit name, start and end of domain envelope
# use domain E-value to compare threshold
  dict_abovethr = {"domain_query":[],"seq_hit":[],"env_from":[],"env_to":[], "database":[]}
  for i in range(length_dict_list):
    if hmmsearch_dict["e-value_dom"][i] < E_thr:
      dict_abovethr["domain_query"].append(hmmsearch_dict["domain_query"][i])
      dict_abovethr["seq_hit"].append(hmmsearch_dict["seq_hit"][i])
      dict_abovethr["database"].append(hmmsearch_dict["database"][i])
      dict_abovethr["env_from"].append(hmmsearch_dict["env_from"][i])
      dict_abovethr["env_to"].append(hmmsearch_dict["env_to"][i])
# if sequences to extract
elif args.D == False and args.S == True:
  dict_abovethr = {"domain_query":[],"seq_hit":[], "database":[]}
  for i in range(length_dict_list):
    if hmmsearch_dict["e-value_seq"][i] < E_thr:
      dict_abovethr["domain_query"].append(hmmsearch_dict["domain_query"][i])
      dict_abovethr["seq_hit"].append(hmmsearch_dict["seq_hit"][i])
      dict_abovethr["database"].append(hmmsearch_dict["database"][i])

length_abvt_dict_list=len(dict_abovethr["domain_query"])
for key in dict_abovethr.keys():
  if len(dict_abovethr[key]) != length_abvt_dict_list:
    sys.exit(f"Could not match sequence information reliably. List of query seq contains {length_abvt_dict_list} domains, while {dict_abovethr[key]} contains {len(dict_abovethr[key])}")


print(f"{length_abvt_dict_list} sequences kept")

##############################################################
# Extract sequences and save them to files

# for each database present in dict_abovethr (one time only)
set_databases=set(dict_abovethr["database"])
domfile_names=set()

database={}
for db in set_databases:
  print(f"Extracting info from database {db}")
  # create a temporary dictionary with all sequences from each database - eliminating breaklines if present
  
  # based on Han Genome2OR
  # open file and intialize the variables
  # Flag is to mark when a sequence is 'registered' in the dict and a new sequence starts
  with open(db) as proteome:
    seqname,seqline,flag="","",False
    # read line by line
    line = proteome.readline()
    while line:
      # encountering a new defline
      if line[0] == '>':
        #  save the last sequence in the dictionary
        if flag:
          database[seqname]=seqline
          # then set flag back to False and seqname and seqline back to empty strings
          flag=False
        seqline=""
        # name of the next sequence from current defline w/o > and \n at the end
        seqname = line[1:]
        if seqname[-1:] == "\n":
          seqname=seqname[:-1]
      # if the line is not a defline
      else:
        # add its content to current seqline
        flag=True
        seqline += line.strip()
      line=proteome.readline()
    # process last line (no > after it)
    if seqname in database.keys():
      print(f"WARNING!!!! {seqname} duplicate in two databases")
    database[seqname]=seqline

print("Extracting domains/sequences")
  
# Create and/or open files to append domain sequences
for j in range(length_abvt_dict_list):

  domain=dict_abovethr["domain_query"][j]
  seq_name=dict_abovethr["seq_hit"][j]
  if "env_from" in dict_abovethr.keys() and "env_to" in dict_abovethr.keys():
    fr=dict_abovethr["env_from"][j]
    to=dict_abovethr["env_to"][j]+1
    seqseq=database[seq_name][fr:to]
  else:
    seqseq=database[seq_name]
  
  # if extracting domains: output files with domain names
  dom_file=main_output_dir + "/" + domain + ".fasta"
  
  # if extracting sequences: only one output file 
  seq_file=main_output_dir + "/"+"hmmsearch_sequences.fasta"
  
  # set output
  if args.D == False and args.S == True:
    out_file = seq_file
  elif args.D == True and args.S == False:
    out_file = dom_file
    domfile_names.add(dom_file)
    
  # append sequence to file if not already present
  seq_present = set()
  with open(out_file,'a+') as fastafile:
    # move to beginning of file and read every line
    fastafile.seek(0,0)
    lines = fastafile.readlines()
    for c in range(len(lines)):
      # add every sequence as defline%sequence [1:] to remove '>', [:-1] to remove '\n'
      if lines[c][0] == ">":
        seq_present.add(lines[c][1:-1]+"%"+lines[c+1][:-1])

    # if sequence not present in set seq_present, append it to file
    item = seq_name+"%"+seqseq
    if item not in seq_present:
      fastafile.write(">"+seq_name+"\n")
      fastafile.write(seqseq+"\n")

if args.D == False and args.S == True:
  print(f"Sequences appended to file {seq_file}")
elif args.D == True and args.S == False:
  print(f"Sequences appended to files: {domfile_names}")
  print(f"Files named after domains, sequences keep name from database")

_end_=time.time()
time_used=_end_ - _start_
print(f"\nTime needed: {time_used} s\n")
