#!/usr/local/bin/python
#
# Given a list of reads or a fasta file along with a list of files, counts all the occurences of
# the reads in the input list of files.  If the inputs are blast files, there is also the option
# of printing out the hits to an output file.

# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)



import sys
from getopt import getopt
from glob import glob

import sequenceUtils

HELP_STRING = """Given a list of reads or a fasta file along with a list of files, counts all the occurences of
the reads in the input list of files.  If the inputs are blast files, there is also the option
of printing out the hits to an output file.  

    -r    Input read file
    -s    Search file pattern
    -o    Whether to create an output file for each file matching the search pattern with the hits.  Default=False
"""

def main(argv=None):

    if argv is None:
        argv = sys.argv

    inputReads = None
    inputFilePattern = None
    outputFiles = False

    try:
        optlist, args = getopt(argv[1:], "hr:s:o:")
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
        elif opt == "-r":
            inputReads = opt_arg
        elif opt == "-s":
            inputFilePattern = opt_arg
        elif opt == "-o":
            outputFiles = opt_arg.startswith("t") or opt_arg.startswith("T")

    if not inputReads:
        print "You must use '-r' to enter a file of read names to search for."
        sys.exit(1)
    if not inputFilePattern:
        print "You must use '-s' to enter a pattern for the files to search through (will be used in 'glob')."
        sys.exit(1)


    inputFiles = glob(inputFilePattern)
    if len(inputFiles) == 0:
        print "The file pattern '%s' did not match any files." % (inputFilePattern)
    else:
        print "The file pattern '%s' matched the following files:" % (inputFilePattern)
        print "\n".join(inputFiles)

    print
    print
    findReads(inputReads, inputFiles, outputFiles)

def findReads(inputReads, inputFiles, outputFiles):

    ir = open(inputReads)
    line = ir.readline()
    isFasta = line.startswith(">")
    ir.close()

    reads = {}
    if isFasta:
        for (title, sequence) in sequenceUtils.FastaIterator(inputReads):
            reads[title] = 1
    else:
        print "The input file is not fasta.  It must be a blast file or a list of read titles..."
        for line in open(inputReads):
            query = line.split()[0]
            reads[query] = 1

    print "countFound\tcountTotal\tFile\tisFasta"    
    for f in inputFiles:
        countFound = 0
        countTotal = 0
        test = open(f)
        line = test.readline()
        isFasta = line.startswith(">")

        if isFasta:
            for (title, sequence) in sequenceUtils.FastaIterator(f):
                countTotal += 1
                if reads.has_key(title):
                    countFound += 1
        else:
            out = None
            if outputFiles:
                out = open(f+".found.txt", "w")
            for line in open(f):
                countTotal += 1
                query = line.split()[0]
                if reads.has_key(query):
                    countFound += 1
                    if out:
                        out.write(line)

        print "%s\t%s\t%s\t%s" % (countFound, countTotal, f, isFasta)

##############################################
if __name__ == "__main__":
    sys.exit(main(None))
                


