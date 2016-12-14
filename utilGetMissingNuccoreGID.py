# -*- coding: utf-8 -*-

import os
import string
import re
import sys

from systemSettings import *

#https://www.ncbi.nlm.nih.gov/books/NBK179288/
#defined in systemSettings.py
command = efetchPath+''' -db nuccore -id "XXX" -mode xml > test.out'''


def extractTaxa( aLine ):
  right = aLine.split(">")[1]
  left = right.split("<")[0]
  return int(left)


if __name__== "__main__":

  #os.system("module unload perl")
  #os.system("module load perl/5.22.1")

  if len(sys.argv) != 2:
      print "Takes one argument:  IMSA pipeline's output of unidentified contigs"
      print "Example:  A50_sub_assemble.fa.tax_unidentifiedTaxaAlignments.txt"
      print ""

  outF = open(EXTRA_BLAST_TAX_DB,'a')

  inf = open(sys.argv[1])

  gis = {}
  for line in inf:
    splits = line.split()
    gid=splits[1].split("|")[1]
    gis[gid]=1

  print "there are %d unique gid's" % len(gis.keys())
  #updated to not repeat gid lookup for multiple cases
  for gid in gis:
    os.system(string.replace(command, "XXX", gid))
    innit = open("test.out")
    taxonFound = False
    taxaID=""
    for line2 in innit:
      if re.search("taxon",line2) != None:
        taxonFound = True
      if taxonFound and re.search("Object-id_id",line2) != None:
        taxaID = extractTaxa( line2 )
        break
    innit.close()
    if taxaID != "":
        outF.write( gid+"\t"+str(taxaID)+"\n" )
        print >> sys.stderr, gid+"\t"+str(taxaID)
    else:
        print "  WARNING ", gid, " was not found"
