#!/usr/local/bin/python
#
# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)

import sys
from getopt import getopt
import os
import eutilsWrapper
import pipelineUtils

HELP_STRING = """Adds a division column (the second column) to a species (or genus) cluster table.

   -i   cluster table
   -o   output file name
"""

def main(argv=None):

    if argv is None:
        argv = sys.argv

    clusterTable = None
    outputName = None

    try:
        optlist, args = getopt(argv[1:], "hi:o:")
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
            clusterTable = opt_arg
        elif opt == "-o":
            outputName = opt_arg

    if not clusterTable or not outputName:
        print "You must include the cluster table file (-i) and the output file name (-o)."
        sys.exit(1)

    out = open(outputName, "w")

    addDivisionColumn(clusterTable, outputName)

def addDivisionColumn(clusterTable, outputName):

    # first get all the species (or genus, etc) ids
    ids = {}
    for line in open(clusterTable):
        name = line.split()[0]
        try:
            id = int(name.split("_")[0])
            ids[id] = 1
        except:
            print "Unable to get tax id for '%s'.  Skipping." % (name)

    # then look up the taxonomy for each id
    print "Looking up taxonomy records for each id"
    taxRecords = eutilsWrapper.getTaxa(ids.keys())

    # then go back through the file and print the new second column
    out = open(outputName, "w")
    count = 0
    print "Writing output file"
    for line in open(clusterTable):
        pieces = line.split()        
        count += 1
        if count == 1:
            out.write("%s\t%s\t%s\n" % (pieces[0], "Division", "\t".join(pieces[1:])))
            continue

        id = int(pieces[0].split("_")[0])
        try:
            division = taxRecords[id].getDivision()
        except:
            raise
        
        out.write("%s\t%s\t%s\n" % (pieces[0], division, "\t".join(pieces[1:])))

    print "Wrote %s lines to the output file." % (count)

##############################################
if __name__ == "__main__":
    sys.exit(main(None))
                


