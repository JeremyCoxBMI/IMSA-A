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
# at a chosen clade level.
#
###################################################


import sys


def main(argv=None):
    if argv is None:
        argv = sys.argv

    if (len(argv) != 5):
        print "\tProper usage:\n\tpython "+argv[0]+" <clade selection> <fasta file> <unique alignments file> <output file>"
        print "\t\tclade selection choices: firstTaxa, species, genus, or family"
        exit(0)

    outF = open(argv[4],'w')
    #uniq = open(argv[3])
    #fasta = open(argv[2])
    clade = argv[1].strip()     #   I am paranoid about whitespace

    listQueryKeep = {}

    for line in open(argv[3]):
        splits = line.split('\t')
        if splits[1] == clade:
            listQueryKeep[ splits[0] ] = splits[2].strip()

    write = False
    for line in open(argv[2]):
        if line[0] == ">":
            short = line[1:].strip().split()[0]  #inchworm compatible
            if short in listQueryKeep:
                write = True
                line = ">"+short+" "+clade+":taxonID:"+listQueryKeep[short]+"\n"
            else:
                write = False
        if write:
            outF.write(line)

    outF.close()

##############################################
if __name__ == "__main__":
    sys.exit(main(None))