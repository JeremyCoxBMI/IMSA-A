'''
Copyright (C) 2016  Jeremy Wayne Cox    jeremy.cox@cchmc.org

This code is distributed under Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0) License.

For a human readable interpretation of the license, please see
http://creativecommons.org/licenses/by-nc-sa/4.0/

NOTICE
This code is a modification and extension of the program Arron-IMSA
distributed by sourceforge.net.
https://sourceforge.net/projects/arron-imsa/

The IMSA license file and notice of copyright is included with this release.
'''


###################################################
# Utility program: allows user to quickly create a fasta file containing only the sequences that created unique hits
# to a certain Taxon ID.
#
###################################################


import sys


def main(argv=None):
    if argv is None:
        argv = sys.argv

    if (len(argv) != 3):
        print "\tProper usage:\n\tpython "+argv[0]+" <assembled fasta file> <TaxonID>"

        exit(0)


    #uniq = open(argv[3])
    #fasta = open(argv[2])
    clade = argv[1].strip()     #   I am paranoid about whitespace

    taxonID = int(argv[2])

    listQueryKeep = {}

    for line in open(argv[1]+".tax_uniqueAlignments.txt"):
        splits = line.strip().split('\t')
        if int(splits[2]) == taxonID:
            listQueryKeep[ splits[0] ] = 1

    write = False
    for line in open(argv[1]):
        if line[0] == ">":
            short = line[1:].strip().split()[0]  #inchworm compatible
            if short in listQueryKeep:
                write = True
                line = ">"+short
            else:
                write = False
        if write:
            print line.strip()


##############################################
if __name__ == "__main__":
    sys.exit(main(None))