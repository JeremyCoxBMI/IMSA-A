#!/usr/local/bin/python
#
# This script just runs the "reportTaxonomy" piece of the IMSA pipeline

# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)

import sys
from getopt import getopt
import os

import pipelineUtils

HELP_STRING = """Runs the reportTaxonomy piece of the IMSA pipeline.

Required Inputs:
    -b    Blast output file
    -d    Optional.  Prefix for the bubble diagrams.  If none is given then bubble diagrams are not generated.
"""

def main(argv=None):

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt(argv[1:], "hb:d:")
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
    bubblePrefix = None
    
    for (opt, opt_arg) in optlist:
        if opt == "-h":
            print ""
            print HELP_STRING
            sys.exit(1)

        elif opt == "-b":
            blastResults = opt_arg

        elif opt == "-d":
            bubblePrefix = opt_arg
                                                                
    if not blastResults:
        print "You must specify a blast result file (use -b)."
        print
        print HELP_STRING
        sys.exit(1)

    # report taxonomy
    rootName, extension = os.path.splitext(blastResults)
    reportName = rootName + ".tax_"
    pipelineUtils.reportTaxonomy(blastResults, printTargets=True, outputPrefix=reportName, dotPrefix=bubblePrefix)


##############################################
if __name__ == "__main__":
    sys.exit(main(None))

            
