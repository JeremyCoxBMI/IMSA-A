# -*- coding: utf-8 -*-

# input line: 
# <program> <species or genus>

# 2015.12.23 UPDATE:  Taxa output file has an additional column in output, the program is changed to accomodate this
# EXAMPLE python util_IMSA_count4_summarize.py firstTaxon.IMSA+A_4count.txt
# 2017.03.02 UPDATE:  Fixing a bug related to new IMSA+A_4count reports, adding an IMSA+A count only feature


import sys
import os
import string
#import queue

if len(sys.argv) != 2:
  print "Syntax is %s <string>" % sys.argv[0]

name = sys.argv[1] 

os.system("ls *"+name+"*4count.txt > temp.in")
inList = open("temp.in",'r')

#inList = os.listdir(".")

outFile = open(name.split(".")[0]+".summarize.txt",'w')
outFile2 = open(name.split(".")[0]+".summarize.IMSA+A.count.txt",'w')

fileNames = []
dataDicts = [] # list of dictionaries
taxaTOname = {}

numColumns = 0
clade = ""

for line in inList:
  line = string.replace(line,"\n","")
  #if line
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
            (taxaKey, taxaName, junk, clade, count1, count2, count3, count4) = line2.split('\t')
          if taxaKey not in taxaTOname:
              taxaTOname[taxaKey] = (taxaName, clade)
          countD[taxaKey] = count1+"\t"+count2+"\t"+count3+"\t"+count4
      else:
        skipFirst = False
  dataDicts.append(countD)

# OUTPUT Header Line
outFile.write("Taxa ID\tTaxa Name\t")
outFile2.write("Taxa ID\tTaxa Name")
if numColumns >= 7:
    outFile.write("Kingdom\t")
    outFile2.write("\tKingdom")
for file in fileNames:
    outFile.write(file+"\tUnique hits\tPartial hits\tPartial sum\t")
    outFile2.write("\t"+file)
outFile2.write("\tsum")
outFile.write("\n")
outFile2.write("\n")

# OUTPUT Aligned files lines
for taxaKey in taxaTOname:
    (taxaName, clade) = taxaTOname[taxaKey]
    if numColumns == 6:
        outFile.write("%s\t%s\t" % (taxaKey, taxaName) )
        outFile2.write("\t%s\t%s" % (taxaKey, taxaName) )
    elif numColumns >= 7:
        outFile.write("%s\t%s\t%s\t" % (taxaKey, taxaName, clade) )
        outFile2.write("%s\t%s\t%s" % (taxaKey, taxaName, clade) )


    sum=0
    for td in dataDicts:
        if taxaKey in td:
            string4 = td[taxaKey]
            uniq = string4.split()[1]
            outFile.write(string4+"\t")
            outFile2.write("\t"+uniq)
            sum+=int(uniq)
        else:
            outFile.write("\t\t\t\t")
            outFile2.write("\t")
    outFile2.write("\t"+str(sum))
    outFile.write("\n")
    outFile2.write("\n")

outFile.close()
outFile2.close()