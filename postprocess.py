#!/usr/local/bin/python
#
# This script is designed to do the typical processing done after the IMSA pipeline completes.
# This script takes in the blast vs. NT, filters to best hits, and runs the taxonomy report.

# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)


import sys
from getopt import getopt
import os

import pipelineUtils
import reportPairConcordance
import sequenceUtils
import blastUtils

HELP_STRING = """Given the results of the filtered reads blasted against NT, this script
does the best hit blast filter and then runs the taxonomy report.

Required Inputs:
    -b    Blast output file

Optional Inputs:
    -f    Fasta input file.  If this is given then a fasta file of all viral results will be created.
"""

def main(argv=None):

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt(argv[1:], "hb:f:")
    except:
        print "Illegal option"
        print ""
        print HELP_STRING
        print ""
        e = sys.exc_info()[1]
        print e
        sys.exit(1)

    if len(optlist) == 0:
        print ""
        print HELP_STRING
        sys.exit(1)

    blastResults = None
    fastaInput = None
    for (opt, opt_arg) in optlist:
        if opt == "-h":
            print ""
            print HELP_STRING
            sys.exit(1)

        elif opt == "-b":
            blastResults = opt_arg

        elif opt == "-f":
            fastaInput = opt_arg
                                                                
    if not blastResults:
        print "You must specify a blast result file (use -b)."
        print
        print HELP_STRING
        sys.exit(1)

    doPostProcess(blastResults, fastaInput)

def doPostProcess(blastResults, fastaInput):

    rootName, extension = os.path.splitext(blastResults)
    bestHitsName = rootName + ".bestHits.bln"
    dp = rootName + ".bubble_"
    reportName = rootName + ".tax_"

    # do best hits
    blastUtils.keepAllBestHits(blastResults, bestHitsName)

    # report taxonomy
    pipelineUtils.reportTaxonomy(bestHitsName, printTargets=True, outputPrefix=reportName, dotPrefix=dp)

    if fastaInput:
        pipelineUtils.getFastaForTaxonomy([10239], bestHitsName, fastaInput, rootName+"_viral.fa")
        pipelineUtils.getFastaForTaxonomy([10292], bestHitsName, fastaInput, rootName+"_herp.fa")
        pipelineUtils.getFastaForTaxonomy([11632], bestHitsName, fastaInput, rootName+"_retro.fa")

##############################################
if __name__ == "__main__":
    sys.exit(main(None))

            
