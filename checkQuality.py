#!/usr/local/bin/python
#
# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)

import sys
from getopt import getopt
import os
from imsa import pipelineUtils

HELP_STRING = """Runs the quality filter part of imsa, just to check numbers.

    -i    Input fastq 1
    -j    Input fastq 2
    -n    name of output file (default = log listener default)
    -b    number bases
    -t    Threshold (default = 15)
"""

def main(argv=None):

    if argv is None:
        argv = sys.argv


    inputFastq1 = None
    inputFastq2 = None
    logName = None
    numBases = 3
    threshold = 15
    
    try:
        optlist, args = getopt(argv[1:], "hi:j:n:t:")
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
        elif opt == "-i":
            inputFastq1 = opt_arg
        elif opt == "-j":
            inputFastq2 = opt_arg
        elif opt == "-n":
            logName = opt_arg
        elif opt == "-b":
            numBases = int(opt_arg)
        elif opt == "-t":
            threshold = int(opt_arg)
            
    if not inputFastq1 or not inputFastq2:
        print "You must include the inputFastq1 (-i) and inputFastq2 (-j)"
        sys.exit(1)

    listener = pipelineUtils.getLogListener(logName)
    pipelineUtils.pairedQualityFilter(inputFastq1, inputFastq2, None, None, numBases, threshold, offset=33, keepMixed=False,
                        delimiter="#0/", mylistener=listener)
    
            



##############################################
if __name__ == "__main__":
    sys.exit(main(None))
                


