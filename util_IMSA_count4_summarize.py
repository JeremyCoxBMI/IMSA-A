# -*- coding: utf-8 -*-

# input line: 
# <program> <species or genus>

# 2015.12.23 UPDATE:  Taxa output file has an additional column in output, the program is changed to accomodate this
# EXAMPLE python util_IMSA_count4_summarize.py firstTaxon.IMSA+A_4count.txt


import sys
import os
import string
#import queue

if len(sys.argv) != 2:
  print "Syntax is %s <string>" % sys.argv[0]

name = sys.argv[1] 

os.system("ls *"+name+"* > temp.in")

inList = open("temp.in",'r')

outFile = open(name.split(".")[0]+".summarize.txt",'w')

fileNames = []
dataDicts = [] # list of dictionaries
taxaTOname = {}

numColumns = 0
clade = ""

for line in inList:
  line = string.replace(line,"\n","")
  fileNames.append(line)
  tempin = open(line,'r')
  countD = {}
  skipFirst = True
  for line2 in tempin:
      line2 = string.replace(line2,'\t\n',"")
      line2 = string.replace(line2,'\n',"")
      #print(">> "+line2+"  ::  ",line2.split('\t') )

      #UPDATE 12/8/2015
      numColumns = len(line2.split('\t'))
      #print "detected number columns: ", numColumns

      if not skipFirst:
          if numColumns == 6:
            (taxaKey, taxaName, count1, count2, count3, count4) = line2.split('\t')
          elif numColumns == 7:
            (taxaKey, taxaName, clade, count1, count2, count3, count4) = line2.split('\t')
          elif numColumns == 8:
            (taxaKey, taxaName, clade, count1, count2, count3, count4, junk) = line2.split('\t')
          if taxaKey not in taxaTOname:
              taxaTOname[taxaKey] = (taxaName, clade)
          countD[taxaKey] = count1+"\t"+count2+"\t"+count3+"\t"+count4
      else:
        skipFirst = False
  dataDicts.append(countD)

# OUTPUT Header Line
outFile.write("Taxa ID\tTaxa Name\t")
if numColumns == 7:
    outFile.write("clade\t")
for file in fileNames:
    outFile.write(file+"\tUnique hits\tPartial hits\tPartial sum\t")
outFile.write("\n")

# OUTPUT Aligned files lines
for taxaKey in taxaTOname:
    (taxaName, clade) = taxaTOname[taxaKey]
    if numColumns == 6:
        outFile.write("%s\t%s\t" % (taxaKey, taxaName) )
    elif numColumns == 7:
        outFile.write("%s\t%s\t%s\t" % (taxaKey, taxaName, clade) )

    for td in dataDicts:
        if taxaKey in td:
            outFile.write(td[taxaKey]+"\t")
        else:
	  outFile.write("\t\t\t\t")
    outFile.write("\n")

outFile.close()