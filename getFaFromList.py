#!/usr/local/bin/python
#
# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)


import sys
from getopt import getopt
import os
from utils import sequenceUtils

HELP_STRING = """Given a list (one title per line) and a fasta input, creates a fasta file of the records

    -l   list
    -f   fasta input
    -o   output file name
    -q   input is fastq instead of fasta
"""

def main(argv=None):

    if argv is None:
        argv = sys.argv

    listInput = None
    fastaInput = None
    outputName = None
    isFastq = False

    try:
        optlist, args = getopt(argv[1:], "hl:f:o:q")
    except:
        print "Illegal option"
        print ""
        print HELP_STRING
        sys.exit(1)

    for (opt, opt_arg) in optlist:
        if opt == "-h":
            print ""
            print HELP_STRING
            sys.exit(1)
        elif opt == "-l":
            listInput = opt_arg
        elif opt == "-f":
            fastaInput = opt_arg
        elif opt == "-o":
            outputName = opt_arg
        elif opt == "-q":
            isFastq = True

    if not listInput or not fastaInput or not outputName:
        print "You must include the list (-l), the fasta initial input (-f) and the output file name (-o)."
        sys.exit(1)


    foundQueries = {}
    for line in open(listInput):
        query = line.strip()
        foundQueries[query] = 1
        
    print "There are %s unique items in your list" % (len(foundQueries))
    #allfound = foundQueries.keys()
    #allfound.sort()
    #for k in allfound:
    #    print "|%s|" % (k)

    out = open(outputName, "w")
    fastaFileHandle = open(fastaInput)
    countFound = 0

    if not isFastq:
        for (title, sequence) in sequenceUtils.FastaIterator(fastaInput):
            found = False
            for k,v in foundQueries.iteritems():
                if title.find(k) >= 0:
                    out.write(">%s\n%s\n" % (title, sequence))
                    countFound += 1
                    found = k

            if found:
                del foundQueries[found]
    else:
        for (title, sequence, quality) in sequenceUtils.FastqIterator(fastaInput):
            found = False
            for k,v in foundQueries.iteritems():
                if title.find(k) >= 0:
                    out.write(">%s\n%s\n" % (title, sequence))
                    countFound += 1
                    found = k

            if found:
                del foundQueries[found]


    print "There were %s items from the list that were found." % (countFound)
    print "The following items were in the list but not found in the input fasta file:"
    for k,v in foundQueries.iteritems():
        print k


##############################################
if __name__ == "__main__":
    sys.exit(main(None))
                


