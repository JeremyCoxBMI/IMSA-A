#!/usr/local/bin/python
#
# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)


import sys
from getopt import getopt
import os
from glob import glob

import pipelineUtils

HELP_STRING = """Given taxonomy IDs and blast results connecting a gi with a read, pulls all the reads from the given taxonomies into
a fasta file.  Runs the pipelineUtils method getFastaForTaxonomy.

    -i  id(s) to pull from the fasta file.  For multiple, should be comma-separated without spaces (i.e. 2,5,10)
    -b  blast input file
    -f  fasta input file
    -o  output fasta file
    
"""

def main(argv=None):

    if argv is None:
        argv = sys.argv

    ids = None
    intIds = None
    blastfileGlob = None
    fastafileGlob = None
    outputfile = None

    try:
        optlist, args = getopt(argv[1:], "hb:f:o:i:")
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
            ids = opt_arg.split(",")
            intIds = [int(i) for i in ids]
        elif opt == "-b":
            blastfileGlob = glob(opt_arg)
        elif opt == "-f":
            fastafileGlob = glob(opt_arg)
        elif opt == "-o":
            outputfile = opt_arg

    if not blastfileGlob or not fastafileGlob or not outputfile or not intIds:
        print "You must include the ids to retrieve (-i), blast file (-b), input fasta file (-f) and the output filename (-o)"
        print blastfileGlob, fastafileGlob, outputfile, intIds
        sys.exit(1)

    #print intIds
    pipelineUtils.getFastaForTaxonomy(intIds, blastfileGlob, fastafileGlob, outputfile)
            



##############################################
if __name__ == "__main__":
    sys.exit(main(None))
                


