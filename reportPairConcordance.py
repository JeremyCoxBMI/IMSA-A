#!/usr/local/bin/python
#
# Given a genome fasta file and a blast file reporting pair read alignments to the genome fasta file,
# reports the pair concordance (i.e. the number of pairs facing each other, average distance between, etc)
# across the length of the genome.

# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)


import sys
from getopt import getopt

import sequenceUtils

HELP_STRING = """Given a blast file reporting read pair alignments to the genome
fasta file, this script reports the pair concordance across the lengths of the genome (i.e. the
number of pairs facing each other, average distance between pairs, etc).

Required Inputs:
      -b    blast output file
      -g    genome size
      -r    read size

Optional Inputs:
      -a    Start Limit.  Only bins happening after this index will be counted
      -z    Stop Limit.  Only bins happening before this index will be counted
      -t    Target in the blast file.  All hits to other targets will be ignored
      -d    Delimiter.  Default = #0/
      -s    Bin size.  Default = 100 bp
      -f    Fasta file of the original read sequences.  If this is given then two files will be created,
               one for potential integration sites off the start of the target sequence and one for
               potential integraion sites off the end.

"""

DEFAULT_DELIMITER = "#0/"
DEFAULT_BIN_SIZE = 100

def main(argv=None):

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt(argv[1:], "hb:g:r:t:d:s:a:z:f:")
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
    genomeSize = None
    readSize = None
    delimiter = DEFAULT_DELIMITER
    targetFilter = None
    binSize = DEFAULT_BIN_SIZE
    startLimit = None
    stopLimit = None
    fastaFile = None

    for (opt, opt_arg) in optlist:
        if opt == "-h":
            print ""
            print HELP_STRING
            sys.exit(1)

        elif opt == "-b":
            blastResults = opt_arg
        elif opt == "-g":
            genomeSize = int(opt_arg)
        elif opt == "-r":
            readSize = int(opt_arg)
        elif opt == "-d":
            delimiter = opt_arg
        elif opt == "-t":
            targetFilter = opt_arg
        elif opt == "-s":
            binSize = int(opt_arg)
        elif opt == "-a":
            startLimit = int(opt_arg)
        elif opt == "-z":
            stopLimit = int(opt_arg)
        elif opt == "-f":
            fastaFile = opt_arg
        
    if not blastResults:
        print "You must specify a blast input file using -b."
        print ""
        print HELP_STRING
        sys.exit(1)

    if not readSize or not genomeSize:
        print "You must specify the genome and read sizes with -g and -r respectively"
        print ""
        print HELP_STRING
        sys.exit(1)

    reportPairConcordance(blastResults, genomeSize, readSize, delimiter, targetFilter, binSize,
                          startLimit, stopLimit, fastaFile)


def reportPairConcordance(blastResults, genomeSize, readSize, delimiter, targetFilter,
                          binSize, startLimit=None, stopLimit=None, fastaFile=None):

    fAlignments = {}
    rAlignments = {}
    for line in open(blastResults):
        pieces = line.split()
        query = pieces[0]
        target = pieces[1]
        if targetFilter:
            if target != targetFilter:
                continue
        qStart = int(pieces[6])
        qEnd = int(pieces[7])
        sStart = int(pieces[8])
        sEnd = int(pieces[9])
        eval = float(pieces[10])
        bitscore = float(pieces[11])

        isForward = sEnd > sStart
        if isForward:
            readStart = sStart - qStart + 1
            fAlignments[query] = readStart
        else:
            readStart = sStart + qStart - 1
            rAlignments[query] = readStart


    # assume possible integration reads are within 500 of the end of the genome.  In the first
    # bin, they are reverse reads with no forward match.  In the last bin they are
    # forward reads with no reverse match.  For now, take note of the name of the
    # 'missing' read.  Then get them from the fastaFile later
    missingStart = {}
    missingEnd = {}

    perBin = {}
    count = 0
    for (query, readStart) in fAlignments.iteritems():
        bin = int(readStart) / int(binSize)
        #print bin, readStart, binSize
        #count += 1
        #if count > 100:
        #    break

        if not perBin.has_key(bin):
            perBin[bin] = {"discordant":0, "discordantR":0, "cordant":0, "circle":0, "sizes":[], "notfound":0, "circlesize":[], "rnotfound":0}
        delimStart = query.find(delimiter)
        if delimStart < 0:
            raise Exception("Unable to find delimiter '%s' in query '%s'" % (delimiter, query))
        rootName, curPair = query.split(delimiter)
        if curPair == "1":
            otherName = rootName + delimiter + "2"
        elif curPair == "2":
            otherName = rootName + delimiter + "1"
        else:
            raise Exception("Expected the post-delimter options of 1 and 2, but found '%s'" % (curPair))

        if fAlignments.has_key(otherName):
            perBin[bin]["discordant"] += 1
        elif not rAlignments.has_key(otherName):
            perBin[bin]["notfound"] += 1
            if readStart > (genomeSize - 500) and fastaFile:
                missingEnd[otherName] = 0
        else:
            otherStart = rAlignments[otherName]
            if otherStart < readStart:
                perBin[bin]["circle"] += 1
                perBin[bin]["circlesize"].append( (genomeSize-(readStart+readSize)) + (otherStart-readSize))
            else:
                ampSize = otherStart - readStart
                perBin[bin]["cordant"] += 1
                perBin[bin]["sizes"].append(ampSize)


    # for the reverse reads, we've already counted pairs so we only need to count the notfound ones
    for (query, readStart) in rAlignments.iteritems():
        bin = readStart / binSize
        if not perBin.has_key(bin):
            perBin[bin] = {"discordant":0, "discordantR":0, "cordant":0, "circle":0, "sizes":[], "notfound":0, "circlesize":[], "rnotfound":0}
        delimStart = query.find(delimiter)
        if delimStart < 0:
            raise Exception("Unable to find delimiter '%s' in query '%s'" % (delimiter, query))
        rootName, curPair = query.split(delimiter)
        if curPair == "1":
            otherName = rootName + delimiter + "2"
        elif curPair == "2":
            otherName = rootName + delimiter + "1"
        else:
            raise Exception("Expected the post-delimter options of 1 and 2, but found '%s'" % (curPair))

        if rAlignments.has_key(otherName):
            perBin[bin]["discordantR"] += 1
        elif not fAlignments.has_key(otherName):
            perBin[bin]["rnotfound"] += 1
            if readStart < 500 and fastaFile:
                missingStart[otherName] = 0


    print "Bin_Start\tNum_Not_Found_Fwd\tNum_Not_Found_Rev\tNum_Discordant_Fwd\tNum_Discordant_Rev\tNum_Good\tAmplicon_Min\tAmplicon_Max\tAmplicon_Mean\tNum_Circles\tCircle_Min\tCircle_Max\tCircle_Mean"
    for bin in perBin.keys():
        binStart = bin * binSize
        useBin = True
        if startLimit and binStart < startLimit:
            useBin = False
        if stopLimit and (binStart + binSize) > stopLimit:
            useBin = False
        if not useBin:
            continue
        try:
            insertMin = min(perBin[bin]["sizes"])
            insertMax = max(perBin[bin]["sizes"])
            insertMean = sum(perBin[bin]["sizes"]) / float(len(perBin[bin]["sizes"]))
        except ValueError:
            insertMin = insertMax = insertMean = 0.0
        try:
            circleMin = min(perBin[bin]["circlesize"])
            circleMax = max(perBin[bin]["circlesize"])
            circleMean = sum(perBin[bin]["circlesize"]) / float(len(perBin[bin]["circlesize"]))
        except ValueError:
            circleMin = circleMax = circleMean = 0.0

        print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%.2f\t" % (binStart, perBin[bin]["notfound"], perBin[bin]["rnotfound"],
                                         perBin[bin]["discordant"], perBin[bin]["discordantR"], perBin[bin]["cordant"],
                                                                            insertMin, insertMax, insertMean,
                                         perBin[bin]["circle"], circleMin, circleMax, circleMean)



    
    if fastaFile:
        print
        print "Creating insertion fasta files..."
        print "There are %s missingStart and %s missingEnd reads" % (len(missingStart), len(missingEnd))
        outStart = open(fastaFile + ".startInsert.fa", "w")
        outEnd = open(fastaFile + ".endInsert.fa", "w")
        numStart = numEnd = 0
        for (title, record) in sequenceUtils.FastaIterator(fastaFile):
            if missingStart.has_key(title):
                numStart += 1
                outStart.write(">%s\n%s\n" % (title, record))
            elif missingEnd.has_key(title):
                numEnd += 1
                outEnd.write(">%s\n%s\n" % (title, record))
        print "Found %s start and %s end reads" % (numStart, numEnd)

##############################################
if __name__ == "__main__":
    sys.exit(main(None))
                
