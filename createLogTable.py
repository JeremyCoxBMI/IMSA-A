#!/usr/local/bin/python
#
# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)

import sys
from getopt import getopt
import os
from glob import glob
from imsa import pipelineUtils

HELP_STRING = """Creates the log table from all the log files that match the input file glob.
Do not forget to include the glob string in quotes, otherwise it won't work.

   -i   log table glob (string with wildcards)
    
"""

def main(argv=None):

    if argv is None:
        argv = sys.argv

    logTableGlob = None

    try:
        optlist, args = getopt(argv[1:], "hi:")
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
            logTableGlob = glob(opt_arg)

    if logTableGlob == None:
        print "You must include the path (with wildcards) for log files (-i)."
        sys.exit(1)

    if len(logTableGlob) == 0:
        print "Your input path did not include any files.  Make sure you use quotes and try again."
        sys.exit(1)

    for f in logTableGlob:
        print "="*60
        print f
        print "-"*20
        pipelineUtils.createLogTable(f)
        


##############################################
if __name__ == "__main__":
    sys.exit(main(None))
                


