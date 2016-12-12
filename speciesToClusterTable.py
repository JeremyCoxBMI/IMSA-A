#!/usr/local/bin/python
#

# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)

import sys
from getopt import getopt
import os
from glob import glob
from imsa import pipelineUtils

HELP_STRING = """Creates a cluster table from the information in the species files.

   -i   species file glob (string with wildcards)
   -o   output file name
   -t   threshold to include a species (default = 1)
    
"""

def main(argv=None):

    if argv is None:
        argv = sys.argv

    speciesFileGlob = None
    outputName = None
    threshold = 1.0

    try:
        optlist, args = getopt(argv[1:], "hi:o:t:")
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
            speciesFileGlob = glob(opt_arg)
        elif opt == "-o":
            outputName = opt_arg
        elif opt == "-t":
            threshold = float(opt_arg)

    if not speciesFileGlob:
        print "You must include the path (with wildcards) for species files (-i)."
        sys.exit(1)

    if not outputName:
        print "You must include an output name (-o)"
        sys.exit(1)

    out = open(outputName, "w")
    speciesToUse = {}
    speciesPerSample = {}

    # first figure out the species to keep
    for f in speciesFileGlob:
        print "%s_%s" % (f.split("_")[0], f.split("_")[1].split("/")[0])
        for line in open(f):
            if line.split()[-1] == "Count":
                continue
            #print line.split("\t")
            [id, name, count] = line.split("\t")
            count = float(count)

            if count >= threshold:
                speciesToUse[id] = name

    # now, create the table.  Build it in memory because the samples are on the columns
    for f in speciesFileGlob:
        for line in open(f):
            if line.split()[-1] == "Count":
                continue

            [id, name, count] = line.split("\t")
            if speciesToUse.has_key(id):
                if not speciesPerSample.has_key(id):
                    speciesPerSample[id] = {}
                speciesPerSample[id][f] = count

    out.write("SpeciesName")
    for f in speciesFileGlob:
        out.write("\t")
        out.write("%s_%s" % (f.split("_")[0], f.split("_")[1].split("/")[0]))
    out.write("\n")
    
    for (speciesId, allcounts) in speciesPerSample.iteritems():
        out.write(speciesId)
        out.write("_")
        out.write(speciesToUse[speciesId].replace(" ", "_"))
        for f in speciesFileGlob:
            out.write("\t")
            if allcounts.has_key(f):
                out.write(allcounts[f].strip())
            else:
                out.write("0")
        out.write("\n")
        


##############################################
if __name__ == "__main__":
    sys.exit(main(None))
                


