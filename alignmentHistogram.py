#!/usr/local/bin/python
#
# Given blast results and a bin size, creates the data to make a histogram of the
# number of reads blasting to each bin.  Ignores any pair information.

# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)

import sys
from getopt import getopt
import os


HELP_STRING = """
Given blast results and a bin size, creates the data to make a histogram of the
number of reads blasting to each bin.  Ignores any pair information.

Required:
    -i   Result file (default is blast format, but sam okay with -s)

Optional:
    -b   Bin size
    -t   Target filter.  Only results to this target will be counted.
          For example, using this flag lets you use results to nt.
    -e   Eval filter.  Only results at this evalue or less will be counted.
    -f   Bit score filter.  Only results with this bit score or higher will be used.
            Can be combined with the eval filter.
    -s   Alignment is in SAM format instead of blast
"""

DEFAULT_BIN_SIZE = 1000

def main(argv=None):

    if argv is None:
        argv = sys.argv

    blastResults = None
    binSize = DEFAULT_BIN_SIZE
    targetFilter = None
    evalFilter = None
    bitScoreFilter = None
    isSamFormat = False

    try:
        optlist, args = getopt(argv[1:], "hi:st:e:b:f:")
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
            binSize = int(opt_arg)
        elif opt == "-i":
            blastResults = opt_arg
        elif opt == "-t":
            targetFilter = opt_arg
        elif opt == "-e":
            evalFilter = float(opt_arg)
        elif opt == "-f":
            bitScoreFilter = float(opt_arg)
        elif opt == "-s":
            isSamFormat = True

    if not blastResults:
        print "You must specify results using '-i'"
        print
        print HELP_STRING
        sys.exit(1)

    alignmentHistogram(blastResults, binSize, targetFilter, evalFilter, bitScoreFilter, isSamFormat)

def alignmentHistogram(blastResults, binSize, targetFilter, evalFilter, bitScoreFilter, isSamFormat):

    binCounts = {}

    numAlignments = 0
    numAlignmentsPassTarget = 0
    numAlignmentsSaved = 0
    
    for line in open(blastResults):
        pieces = line.split()
        if isSamFormat:
            if line.startswith("@"):
                continue
            target = pieces[2]
            eval = 0
            bitScore = int(pieces[4])
            bin = int(pieces[3]) / binSize
            if target != "*":
                numAlignments += 1
                
        else:
            #print pieces
            target = pieces[1]
            eval = float(pieces[10])
            bitScore = float(pieces[11])
            subStart = int(pieces[8])
            subEnd = int(pieces[9])
            if subStart < subEnd:
                bin = subStart / binSize
            else:
                bin = subEnd / binSize
                
        if targetFilter:
            if target != targetFilter:
                continue
        numAlignmentsPassTarget += 1
        
        if evalFilter:
            if eval > evalFilter:
                continue
        if bitScoreFilter:
            if bitScore < bitScoreFilter:
                continue

        if not binCounts.has_key(bin):
            binCounts[bin] = 0
        binCounts[bin] += 1
        numAlignmentsSaved += 1
        
    allkeys = binCounts.keys()
    if len(allkeys) == 0:
        print "There are no results!  If you used -t, your target wasn't found."
        sys.exit(0)
        
    maxBin = max(allkeys)

    print "Bin\tIndex\tCount"
    #sumValues = 0
    for i in range(0, maxBin+1):
        if binCounts.has_key(i):
            print "%s\t%s\t%s" % (i, i*binSize, binCounts[i])
            #sumValues += binCounts[i]
        else:
            print "%s\t%s\t0" % (i, i*binSize)
    #print "Sum is %s" % (sumValues)
    
    print "\nThere were %s alignments, %s passed the target filter and %s were saved.  There are %s reads in the histogram." % (
        numAlignments, numAlignmentsPassTarget, numAlignmentsSaved, sum(binCounts.values()))



##############################################
if __name__ == "__main__":
    sys.exit(main(None))
                


