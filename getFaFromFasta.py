#!/usr/local/bin/python
#
# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)



import sys
from getopt import getopt

import sequenceUtils

HELP_STRING = """Given blast results and the initial fasta input, creates a fasta file of the blast results

    -b   blast results
    -f   fasta initial input
    -o   output file name
    -e   eval threshold (optional)
"""

def main(argv=None):

    if argv is None:
        argv = sys.argv

    blastResults = None
    fastaInput = None
    outputName = None
    evalThreshold = None

    try:
        optlist, args = getopt(argv[1:], "hb:f:o:e:")
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
        elif opt == "-b":
            blastResults = opt_arg
        elif opt == "-f":
            fastaInput = opt_arg
        elif opt == "-o":
            outputName = opt_arg
        elif opt == "-e":
            evalThreshold = float(opt_arg)

    if not blastResults or not fastaInput or not outputName:
        print "You must include the blast results (-b), the fasta initial input (-f) and the output file name (-o)."
        sys.exit(1)


    foundQueries = {}
    for line in open(blastResults):
        query = line.split()[0]
        eval = float(line.split()[10])
        foundQueries[query] = eval
        
    print "There are %s blast results" % (len(foundQueries))
    allfound = foundQueries.keys()
    allfound.sort()
    for k in allfound:
        print "|%s|" % (k)

    out = open(outputName, "w")
    fastaFileHandle = open(fastaInput)
    countFound = 0
    countPassed = 0
    for (title, sequence) in sequenceUtils.FastaIterator(fastaInput):
        found = None
        for k,v in foundQueries.iteritems():
            if title.find(k) >= 0:
                if not evalThreshold or v <= evalThreshold:
                    out.write(">%s\n%s\n" % (title, sequence))
                    countPassed += 1
                countFound += 1
                found = k
        if found:
            del foundQueries[k]

    print "%s of the blast items were found and %s passed the eval threshold" % (countFound, countPassed)
    print "The following items were in the blast results but not found in the input fasta file:"
    for k,v in foundQueries.iteritems():
        print k


##############################################
if __name__ == "__main__":
    sys.exit(main(None))
                


