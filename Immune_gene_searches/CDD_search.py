import sys
import urllib.request
import lxml
from bs4 import BeautifulSoup
from bs4.element import Comment
import os
import re

# python 3.8
#If uncorrect number of arguments, exit and print the correct usage
if len(sys.argv) != 3:
  sys.exit("Not executed python script. Usage python CDD_search_API.py CDD_accession_number_list path")


def main():
  #create directory to store output. If it exists already, exit the program
  directory=sys.argv[2]
  try:
    os.mkdir(directory)
  except FileExistsError:
    sys.exit(f"!!! ERROR !!! Directory {directory} already exists. Please check and rename it.")
    
  #store cdd accession numbers to a list
  ANs=sys.argv[1].split("\n")
  
  #extract alignment for all accession numbers and save them in different files
  for AN in ANs:
    AN=str(AN).rstrip('\n')
    
    #if ANs contains whitespaces (i.e. input file had more than one column), return error
    if re.search("\s", AN):
      os.rmdir(directory)
      sys.exit("Input file should be a single column of accession numbers")

    url="https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid="+AN+"&seqout=1&maxaln=-1"
    print(f"Retrieving alignment from {url}")

    #get alignments and write only the visible text to a new file
    try:
      handle = urllib.request.urlopen(url)
      soup = text_from_html(handle)
      
    except:
      print(f"!!! WARNING! Something went wrong downloading the alignment for {AN} !!!")
    #if invalid input from CDD give a warning and print invalid input message from CDD in the corresponding file
    if re.search("Invalid input", str(soup)):
      print(f"!!! WARNING! Invalid input for {AN} !!!")
    
    # remove breaklines within sequence
      #separate fastas from each other (defline+sequence)
    seqs=soup.split(">")
      #remove first empty element
    seqs.pop(0)
    seq_list=[]
      #split defline from sequence. Create a list of defline, sequence, defline, sequence etc
    for seq in seqs:
      seq_two = seq.split("\n",maxsplit=1)
      seq_list.append(seq_two[0])
      seq_list.append(seq_two[1].replace("\n",""))

    print(f"Saving to {directory}/{AN}.fasta")
    with open(f"{directory}/{AN}.fasta","w") as writer:
      for i in range(len(seq_list)):
        if i%2 == 0:
          writer.write(">"+seq_list[i]+"\n")
        else:
          writer.write(seq_list[i]+"\n")


#from https://stackoverflow.com/questions/1936466/beautifulsoup-grab-visible-webpage-text
#not so necessary as soup only has one paragraph 
def tag_visible(element):
    if element.parent.name in ['style', 'script', 'head', 'title', 'meta', '[document]']:
        return False
    elif isinstance(element, Comment):
        return False
    return True


def text_from_html(body):
    soup = BeautifulSoup(body, 'lxml')
    texts = soup.findAll(text=True)
    visible_texts = filter(tag_visible, texts)
    return u"".join(t.strip() for t in visible_texts)


if __name__ == "__main__":
  main()
