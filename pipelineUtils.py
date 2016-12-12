#!/usr/local/bin/python
#
# This module contains helper functions used in the Arron Lab Solexa Pipeline.
# These functions are fairly high level pieces of the pipeline.  Smaller
# helper functions are in sequenceUtils.py (assumed to be
# in the user's PYTHONPATH).


# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)


FILE_TYPE_FASTA = "fasta"
FILE_TYPE_FASTQ = "fastq"
FILE_TYPE_BLAST = "blast" # -m 8
FILE_TYPE_BLAT = "blat" # psl

LOG_PREFIX = "log"
SUMMARY_PREFIX = "summary"

# general python
import getpass
import time
import os
import commands
import inspect
import datetime
import operator
import random
import subprocess

# IMSA
import listener
import config
import sequenceUtils
import eutilsWrapper
import arronGrid as ag
import config

def getLogListener(filename=None):
    """Creates a log file and returns the open file handle."""
    if not filename:
        filename = "%s_%s_%s.txt" % (LOG_PREFIX, getpass.getuser(), time.strftime("%y_%m_%d_%H_%M_%S"))
    l = listener.getListener(5, filename)
    return l

def getArgString(frame):
    paramString = ""
    args, x, y, values = inspect.getargvalues(frame)
    for i in args:
        paramString += "\n                                         %s = %s" % (i, values[i])
    return paramString    
    

def getFileType(inputName):
    """Determine what type of file the given file is.  For now this is a simple check
    of the first few lines to determine structure.  If necessary we can adjust this 
    function so that it also checks that the whole file is valid."""
    firstLines = []
    i = open(inputName)
    try:
        while len(firstLines) < 4:
            firstLines.append(i.readline())
    except IOError:
        raise Exception("Unable to read at least 4 lines from '%s'." % (inputName))
        
    if firstLines[0].startswith(">"):
        return FILE_TYPE_FASTA
    elif firstLines[0].startswith("@"):
        if firstLines[2].startswith("+"):
            return FILE_TYPE_FASTQ
        else:
            raise Exception("Unclear file format for '%s'.  The first line starts with '@', but the third line doesn't start with '+'." % inputName)
    elif firstLines[0].startswith("psLayout"):
        return FILE_TYPE_BLAT
    else:
        # figure out the number of pieces and assume accordingly
        pieces = firstLines[0].split("\t")
        if len(pieces) == 21:
            return FILE_TYPE_BLAT
        if len(pieces) == 12:
            return FILE_TYPE_BLAST
        
    raise Exception("I don't know what type of file '%s' is.\nThe first four lines are:\n\n%s" % (inputName, "".join(firstLines)))
        

def convertToFasta(inputFastq, outputFasta, mylistener):
    if mylistener:
        mylistener.reportImportantInfo("Converting '%s' (fastq) into '%s' (fasta)" % (inputFastq, outputFasta))
    out = open(outputFasta, "w")
    numRecords = 0
    for (titleStr, seqStr, qualityStr) in sequenceUtils.FastqIterator(inputFastq):
        out.write(">%s\n%s\n" % (titleStr, seqStr))
        numRecords += 1

    mylistener.reportInfo("Converted %s records from fastq to fasta." % numRecords)

def removeNs(inputFasta, outputFasta, n_threshold, rep_threshold, mylistener):
    """Removes reads from the file that contain more than the threshold number of N's (doesn't have
    to be continuous)."""
    if mylistener:
        mylistener.reportImportantInfo("Removing reads with more than %s N's and reads with %s or more of the same base in a row.  Variable values:%s" % (n_threshold, rep_threshold, getArgString(inspect.currentframe())))

    NUM_REPEAT = rep_threshold
    repeats = ["A"*NUM_REPEAT, "C"*NUM_REPEAT, "G"*NUM_REPEAT, "T"*NUM_REPEAT]
    out = open(outputFasta, "w")
    countGood = 0
    countBadRepeat = 0
    countBadN = 0
    countBadSVA =0
    for (title, sequence) in sequenceUtils.FastaIterator(inputFasta):
        if sequence.count("N") <= n_threshold:
            isBad = False
            for r in repeats:
                if r in sequence:
                    isBad = True
                    countBadRepeat += 1
                    break
            if not isBad:
                countGood += 1
                out.write(">%s\n%s\n" % (title, sequence))
        else:
            countBadN += 1

    if listener:
        mylistener.reportInfo("There were %s reads with more than %s N's, %s reads with repeats, and %s reads that passed the filter." % (
        countBadN, n_threshold, countBadRepeat, countGood))


def qualityFilter(inputFastq, outputFastq, numBases, threshold, offset=33, mylistener=None):
    """Removes the low quality reads from inputFastq, writing to outputFastq.  Low quality reads are defined as reads
    where numBases or more bases have a quality score at or below threshold.  Offset is the ascii offset to convert
    the quality letter into a score."""

    if mylistener:
        mylistener.reportImportantInfo("Removing low-quality read pairs with more than %s bases with quality scores less than %s.  Variable values:%s" % (numBases,
                                       threshold,
                                       getArgString(inspect.currentframe())))


    out = open(outputFastq, "w")
    numGood = 0
    numBad = 0
    numPrimer1 = 0
    numPrimer2 = 0
    for (title, sequence, qualityStr) in sequenceUtils.FastqIterator(inputFastq):
        if readIsGood(qualityStr, numBases, threshold, offset, sequence):
            #if sequence.find("ATCGGAAGAGCGGTTCA") >= 0:
            #    numPrimer1 += 1
            #elif sequence.find("GCAGGAATGCCGAGACC") >= 0:
            #    numPrimer2 += 1
            #else:
            numGood += 1
            out.write("@%s\n%s\n+%s\n%s\n" % (title, sequence, title, qualityStr))
        else:
            #print "Found a bad one!"
            #print title
            #print qualityStr
            #print qualityInts
            numBad += 1

    print "There were %s good and %s bad sequences." % (numGood, numBad)

    if mylistener:
        mylistener.reportInfo("There were %s good reads, %s bad reads." % (numGood, numBad))

def pairedQualityFilter(inputFastq1, inputFastq2, outputFastq1, outputFastq2, numBases, threshold, offset=33, keepMixed=False,
                        delimiter="#0/", mylistener=None):
    """A pair-aware quality filter.  Reads are determined to be low quality if more than 'numBases' bases
    have a quality score of 'threshold' or less.  If both reads in a pair are good then the reads are put in
    the output files.  If both are bad, then they are tossed.  If the results are mixed (one good, one bad), then
    they are both kept if keepMixed is True, otherwise they are both tossed.
    The code assumes the read pairs are in the same order in inputFastq1 and inputFastq2, and that order will be
    maintained in the output files.
    """

    if mylistener:
        mylistener.reportImportantInfo("Removing low-quality read pairs with more than %s bases with quality scores less than %s.  Variable values:%s" % (numBases,
        threshold,
        getArgString(inspect.currentframe())))

    if outputFastq1:
        out1 = open(outputFastq1, "w")
        out2 = open(outputFastq2, "w")
    numAllGood = 0
    numMixed = 0
    numAllBad = 0
    numPrimer1first = 0 # number of sequence1 that have the first half of the primer
    numPrimer1second = 0 # number of sequence1 that have the second half of the primer
    numPrimer2first = 0
    numPrimer2second = 0
    
    it2 = sequenceUtils.FastqIterator(inputFastq2)
    for (title1, sequence1, qualityStr1) in sequenceUtils.FastqIterator(inputFastq1):
        (title2, sequence2, qualityStr2) = it2.next()

        title1 = title1.split()[0]
        title2 = title2.split()[0]

        if title1.split(delimiter)[0] != title2.split(delimiter)[0]:
            raise Exception("InputFastq1 and InputFastq2 appear to have reads in different orders.  '%s' != '%s'" % (
                title1, title2))

        good1 = readIsGood(qualityStr1, numBases, threshold, offset, sequence1)
        good2 = readIsGood(qualityStr2, numBases, threshold, offset, sequence2)

        #if sequence1:
        #    if sequence1.find("ATCGGAAGAGCGGTTCA") >= 0:
        #        numPrimer1first += 1
        #    elif sequence1.find("GCAGGAATGCCGAGACC") >= 0:
        #        numPrimer1second += 1

        #if sequence2:
        #    if sequence2.find("ATCGGAAGAGCGGTTCA") >= 0:
        #        numPrimer2first += 1
        #    elif sequence2.find("GCAGGAATGCCGAGACC") >= 0:
        #        numPrimer2second += 1


        if good1 and good2:
            numAllGood += 1
            if outputFastq1:
                out1.write("@%s\n%s\n+%s\n%s\n" % (title1, sequence1, title1, qualityStr1))
                out2.write("@%s\n%s\n+%s\n%s\n" % (title2, sequence2, title2, qualityStr2))
        elif not good1 and not good2:
            numAllBad += 1
        else:
            numMixed += 1
            if keepMixed and outputFastq1:
                out1.write("@%s\n%s\n+%s\n%s\n" % (title1, sequence1, title1, qualityStr1))
                out2.write("@%s\n%s\n+%s\n%s\n" % (title2, sequence2, title2, qualityStr2))

    mylistener.reportInfo("There were %s good pairs, %s bad pairs and %s mixed pairs." % (numAllGood, numAllBad, numMixed))
    #mylistener.reportInfo("Looking at the first read, there were %s with the first half of the primer and %s with the second half." % (numPrimer1first, numPrimer1second))
    #mylistener.reportInfo("Looking at the second read, there were %s with the first half of the primer and %s with the second half." % (numPrimer2first, numPrimer2second))

    if mylistener:
        if keepMixed:
            mylistener.reportInfo("The mixed pairs were kept.")
        else:
            mylistener.reportInfo("The mixed pairs were not kept.")

        

def readIsGood(qualityStr, numBases, threshold, offset, sequence=None):
    #if sequence:
    #    if sequence.find("ATCGGAAGAGCGGTTCA") >= 0:
    #        return False
    #    elif sequence.find("GCAGGAATGCCGAGACC") >= 0:
    #        return False
    
    qualityInts = map(lambda x: ord(x)-offset, qualityStr)
    numBelow = 0
    for x in qualityInts:
        if x > 45 or x < -5:
            raise Exception("You might have the offset wrong.  The current offset is %s.  There is a quality score of %s.\nqualstr: %s\nqualint: %s" % (offset, x, qualityStr, qualityInts))
            #print qualityStr
            #print qualityInts
        if x <= threshold:
            numBelow += 1

    return numBelow < numBases


def combineFiles(files, outputName):
    out = open(outputName, "w")
    for f in files:
        for line in open(f):
            out.write(line)


def runBowtie(inputFastq1, inputFastq2, targetDB, outputBowtie, bowtieParams, tmpDir, nodesToUse,
              bowtieLogFile, mylistener):

    if mylistener:
        mylistener.reportImportantInfo("Running Bowtie.  Variable values:%s" % (getArgString(inspect.currentframe())))

    if inputFastq2 == "None":    
        inputFastq2 = None

    if inputFastq2:
        cmd = "%s --quiet --sam --sam-nohead %s -p %s %s -1 %s -2 %s %s/%s" % (config.PATH_TO_BOWTIE,
                                                        bowtieParams, nodesToUse, targetDB,
                                                        inputFastq1, inputFastq2, tmpDir, outputBowtie)
    else:
        cmd = "%s --quiet --sam --sam-nohead %s -p %s %s %s %s/%s" % (config.PATH_TO_BOWTIE,
                                                        bowtieParams, nodesToUse, targetDB, inputFastq1,
                                                        tmpDir, outputBowtie)

    print "Running command:", cmd
    try:
        p = subprocess.Popen([cmd], shell=True)
        p.wait()
    except:
        raise  


def runBowtie2(inputFastq1, inputFastq2, targetDB, outputBowtie, bowtieParams, tmpDir, nodesToUse,
              outputUnaligned, bowtieLogFile, mylistener):

    if mylistener:
        mylistener.reportImportantInfo("Running Bowtie.  Variable values:%s" % (getArgString(inspect.currentframe())))

    if inputFastq2 == "None":    
        inputFastq2 = None

    if inputFastq2:
        cmd = "%s %s -p %s --un-conc %s -x %s -1 %s -2 %s -S %s/%s" % (config.PATH_TO_BOWTIE2,
                                                        bowtieParams, nodesToUse, outputUnaligned, targetDB,
                                                        inputFastq1, inputFastq2, tmpDir, outputBowtie)
    else:
        cmd = "%s %s -p %s --un %s -x %s -U %s -S %s/%s" % (config.PATH_TO_BOWTIE2,
                                                        bowtieParams, nodesToUse, outputUnaligned, targetDB,
                                                        inputFastq1, inputFastq2, tmpDir, outputBowtie)


    print "Running command:", cmd
    try:
        p = subprocess.Popen([cmd], shell=True)
        p.wait()
    except:
        raise  



def filterBowtieSingleEnd(bowtieResults, inputFasta, outputFasta, mylistener=None):

    if mylistener:
        mylistener.reportImportantInfo("Filtering Bowtie output.  Variable values:%s" % (getArgString(inspect.currentframe())))
            
    # read through the bowtieResults in batches
    # assume bowtie results will be in the same order as the fasta.  Everything in the bowtie will be in the fasta, but
    #   not the reverse
    count = 0
    lastRead = None
    savedMatches = {}
    br = open(bowtieResults)
    line = br.readline()
    input = sequenceUtils.FastaIterator(inputFasta)
    out = open(outputFasta, "w")
    countHits = 0
    countMisses = 0
    while len(line) > 2:
        # go until we've got enough results for a batch.  Keep numbers in X * 1,000,000 to make code easy to read.
        while count < (20 * 1000000) and len(line) > 2:
            readname = line.split("\t")[0]
            savedMatches[readname] = 0
            lastRead = readname
            line = br.readline()
            count +=1

        print "read a batch of length %s" % (len(savedMatches))

        # now, read through input fasta for each saved match
        while True:
            (title, sequence) = input.next()
            if savedMatches.has_key(title):
                countHits += 1
            else:
                countMisses += 1
                out.write(">%s\n%s\n" % (title, sequence))
            if title == lastRead:
                break

        # ok, clear and prep for next batch
        savedMatches = {}
        count = 0
        print "  Done with a batch.  So far there are %s hits and %s misses." % (countHits, countMisses)

    print "Done!  There were %s hits and %s misses." % (countHits, countMisses)

    if (mylistener):
        mylistener.reportInfo("In bowtieFilter, there were %s hits and %s misses" % (countHits, countMisses))


def filterUsingBowtie(inputFastq1, inputFastq2, bowtieResults, inputType=FILE_TYPE_FASTQ,
                    delimiter=None, pairedHits=None, pairedMisses=None, printHits1=None, printHits2=None,
                      printMisses1=None, printMisses2=None,
                    divideHalf=False, mylistener=None):
    """Given an input fasta and a set of blat results, filters the fasta so only those fasta
    records without hits (as defined by params) are in the outputFasta.
    """

    if mylistener:
        mylistener.reportImportantInfo("Filtering Bowtie output.  Variable values:%s" % (getArgString(inspect.currentframe())))

    countHits, countMisses = doFilterUsingBowtie(inputFastq1, inputFastq2, bowtieResults, inputType, 
                    delimiter, pairedHits, pairedMisses, printHits1, printHits2, printMisses1,
                    printMisses2, mylistener)

    if mylistener:
        mylistener.reportInfo("In bowtieFilter, there were %s hits and %s misses" % (countHits, countMisses))

    #print countHits, countMisses

    return (0)


def doFilterUsingBowtie(inputFastq1, inputFastq2, bowtieResults, inputType=FILE_TYPE_FASTQ, 
                    delimiter=None, pairedHits=None, pairedMisses=None, printHits1=None, printHits2=None,
                    printMisses1=None, printMisses2=None, mylistener=None):
    """Actually does the work."""
    if inputType != FILE_TYPE_FASTQ and inputType != FILE_TYPE_FASTA:
        raise Exception("The file type for doFilterUsingBowtie is expected to be FASTA or FASTQ.  File type was %s" % inputType)

    if inputFastq2 == "None":
        inputFastq2 = None
    if delimiter == "None":
        delimiter = None
        
    #if mylistener:
    #    mylistener.reportInfo("BOWTIE ALIGNMENT IS UNPAIRED!")
    #delimiter = None
        
    outHits1 = outHits2 = outMisses1 = outMisses2 = None
    if printHits1:
        outHits1 = open(printHits1, "w")
    if printHits2:
        outHits2 = open(printHits2, "w")
    if printMisses1:
        outMisses1 = open(printMisses1, "w")
    if printMisses2:
        outMisses2 = open(printMisses2, "w")

    bowtieResultsFH = open(bowtieResults)
    # PLAN:  read in a large number (10 million) hits, making sure to end on a paired end properly.  Then process these results.
    #        repeat, as necessary, until the whole bowtieResults file is read.  This depends on the results being in SAM format,
    #        where both hits and misses are reported.  Also, we have to re-add the "/1" and "/2" onto the read name.  For now
    #        assume it's "/1" and "/2"
    line = bowtieResultsFH.readline()
    totalCountHits = totalCountMisses = 0


    if delimiter:
        d = True
    else:
        d = False
    print "Delimiter is '%s' which is %s" % (delimiter, d)
    
    while len(line) > 2:
        hits = {}
        misses = {}
        count = 1
        #  keeping the count as 100*1000000 to make it easy to read how many millions)
        while count < (50 * 1000000) and len(line) > 2:
            if count % (10 * 1000000) == 0:
                print count, time.strftime("%y_%m_%d_%H_%M_%S")
            fields = line.split("\t")
            query=fields[0].split()[0]
            flag = int(fields[1])
            pairMapped = flag & 2 == 2
            readUnmapped = flag & 4 == 4
            mateUnmapped = flag & 8 == 8
            pairedRead = flag & 1 == 1
            isFirstRead = flag & 64 == 64
            isSecondRead = flag & 128 == 128

            # FIX ME!!!  WE NEED THE DELIMITER AND THE SAM OUTPUT TO BE MUCH LESS HACK-Y
            if query[-2:] == "#0":
                query = query[:-2]
            if delimiter:
                if pairedRead and isFirstRead:
                    query = query + delimiter + "1"
                elif pairedRead and isSecondRead:
                    query = query + delimiter + "2"

            if inputFastq2 != None:
                if pairedHits:
                    if not readUnmapped and not mateUnmapped:
                        hits[query] = None
                    else:
                        misses[query] = None
                else:
                    if readUnmapped and mateUnmapped:
                        misses[query] = None
                    else:
                        hits[query] = None
            else:
                if readUnmapped:
                    misses[query] = None
                else:
                    hits[query] = None


            line = bowtieResultsFH.readline()
            count += 1

        if inputType == FILE_TYPE_FASTQ:
            countHits1, countMisses1 = iterateBowtieResults(inputFastq1, delimiter, hits, misses, outHits1, outMisses1)
            if inputFastq2:
                countHits2, countMisses2 = iterateBowtieResults(inputFastq2, delimiter, hits, misses, outHits2, outMisses2)
            else:
                countHits2 = countMisses2 = 0
        else:
            countHits1, countMisses1 = iterateBowtieResultsFasta(inputFastq1, delimiter, hits, misses, outHits1, outMisses1)
            if inputFastq2:
                countHits2, countMisses2 = iterateBowtieResultsFasta(inputFastq2, delimiter, hits, misses, outHits2, outMisses2)
            else:
                countHits2 = countMisses2 = 0

        #print countHits1, countHits2, countMisses1, countMisses2

        totalCountHits += countHits1+countHits2
        totalCountMisses += countMisses1+countMisses2

    return totalCountHits, totalCountMisses

def iterateBowtieResults(inputFastq, delimiter, hits, misses, outHits, outMisses):
    countHits = 0
    countMisses = 0
    totalCount = 0
    #print "in iterate", len(hits), len(misses)
    #print hits.keys()[:5]
    for (title, sequence, quality) in sequenceUtils.FastqIterator(inputFastq):
        read = title.split()[0]
        totalCount += 1
        #if totalCount < 10:
        #    print "in iterate:", read
        if read in hits:
            countHits += 1
            if outHits:
                outHits.write("@%s\n%s\n+%s\n%s\n" % (read, sequence, read, quality))
        if read in misses:
            countMisses += 1
            if outMisses:
                outMisses.write("@%s\n%s\n+%s\n%s\n" % (read, sequence, read, quality))

    return countHits, countMisses


def iterateBowtieResultsFasta(inputFasta, delimiter, hits, misses, outHits, outMisses):
    countHits = 0
    countMisses = 0
    for (title, sequence) in sequenceUtils.FastaIterator(inputFasta):
        read = title.strip()
        if read in hits:
            countHits += 1
            if outHits:
                outHits.write(">%s\n%s\n" % (read, sequence))
        if read in misses:
            countMisses += 1
            if outMisses:
                outMisses.write(">%s\n%s\n" % (read, sequence))

    return countHits, countMisses


def getNumFromTitle(title):
    # works for
    #   HWUSI-EAS053R_0021_FC:1:1:4834:1047#0/1
    #   HWI-EAS318:5:1:10:21#0/1:N
    #print title, title.split(":"), title.split(":")[3]

    
    return int(title.split(":")[3])
    #return int(title.split(":")[5])
        
def determineDivisionPoints(readFile, numDivisions=2, isFasta=True, binSize=1):
    """Given a fasta file of reads, returns the best number to use to divide the file in half based
    on the read name.  For example, in the read 'HWUSI-EAS053R_0026_FC709NYAAXX:7:1:1809:1239#0/1'
    the number 1809 is used."""

    bins = {}
    #count = 0
    maxNum = 0
    if isFasta:
        it = sequenceUtils.FastaIterator(readFile)
    else:
        it = sequenceUtils.FastqIterator(readFile)

    for record in it:
        title = record[0]
        num = getNumFromTitle(title)
        if num > maxNum:
            maxNum = num
        bin = num / binSize
        if not bins.has_key(bin):
            bins[bin] = 0
        bins[bin] += 1

    sizeOfEachDivision = sum(bins.values()) / numDivisions
    maxBin = max(bins.keys())

    totalSinceLastDivision = 0
    divisionPoints = []
    for i in range(maxBin):
        if bins.has_key(i):
            totalSinceLastDivision += bins[i]
            if totalSinceLastDivision > sizeOfEachDivision:
                divisionPoints.append(i * binSize)
                totalSinceLastDivision = 0

    divisionPoints.append(maxNum)
    
    if len(divisionPoints) != (numDivisions):
        for k,v in bins.iteritems():
            print k, v
        raise Exception("Unable to find division points.  The getNumFromTitle function may not be parsing the read names properly.  SizeOfEachDivision=%s, maxBin=%s, divisionPoints=%s" % (sizeOfEachDivision, maxBin, divisionPoints))

    return divisionPoints

                 
def divideFile(inputFile, divisionPoints, fileType='fasta'):
    """Divides fasta, fastq, blast, and blat files by the divisionPoints set.
    The accepted file type are 'fasta', 'fastq', 'blast', and 'blat'."""

    if fileType not in ['fasta', 'fastq', 'blast', 'blat']:
        raise Exception("In divideFile, the fileType must be one of 'fasta', 'fastq', 'blast', or 'blat'.  The type '%s' is not recognized." % fileType)

    rootName, extension = os.path.splitext(inputFile)
    outNames = []
    outFiles = []
    for i in range(len(divisionPoints)):
        outName = rootName+"_%s." % (i + 1)+extension
        outNames.append(outName)
        outFiles.append(open(outName, "w"))

    countsPerDivision = [0] * len(divisionPoints)
    if fileType == 'fasta':
        it = sequenceUtils.FastaIterator(inputFile)
    elif fileType == 'fastq':
        it = sequenceUtils.FastqItereator(inputFile)
    else:
        it = open(inputFile)
    for record in it:
        # skip blanks
        if len(record) < 3:
            continue
        if fileType == 'blat' and record.startswith("psLayout"):
            continue
        if fileType == 'blat' and record.find("match") >= 0:
            continue
        if fileType == 'blat' and record.startswith("-------"):
            continue
        try:
            if fileType == 'fasta' or fileType == 'fastq':
                title = record[0]
                sequence = record[1]
            elif fileType == 'blast':
                title = record.split()[0]
            elif fileType == 'blat':
                title = record.split()[9]
            else:
                raise Exception("Unrecognized file type '%s'" % (fileType))
            num = getNumFromTitle(title)
        except:
            print "Error parsing line:"
            print record
            raise

        written = False
        for i in range(len(divisionPoints)):
            if num <= divisionPoints[i]:
                written = True
                countsPerDivision[i] += 1
                if fileType == 'fasta':
                    outFiles[i].write(">%s\n%s\n" % (title, sequence))
                elif fileType == 'fastq':
                    quality = record[2]
                    outFiles[i].write("@%s\n%s\n+%s\n%s\n)" % (title, sequence, title, quality))
                else:
                    outFiles[i].write(record)
                break
        if not written:
            raise Exception("The title '%s' was above the maximum division point.  Record:" % (title, record))
                 
    print "Done dividing file '%s'.  The counts per division were: %s" % (inputFile, countsPerDivision)

    return outNames
            
def runBlat(inputFasta, targetDB, outputBlat, params, oocFile, tmpDir, nodesToUse, blatLogFile, mylistener):
    """Runs gridBlat to run the blat if USE_GRID is on, otherwise runs the blat straight"""
    if mylistener:
        mylistener.reportImportantInfo("Running blat.  Variable values:%s" % (getArgString(inspect.currentframe())))

    # check things we always want to check
    if not os.path.exists(inputFasta):
        raise Exception("The input file for BLAST '%s' does not exist" % (inputFasta))

    # submit to grid or straight
    if config.USE_GRID:
        runGridBlat(inputFasta, targetDB, outputBlat, params, oocFile, tmpDir, nodesToUse, blatLogFile, mylistener)
    else:
        runStraightBlat(inputFasta, targetDB, outputBlat, params, oocFile, mylistener)

def runGridBlat(inputFasta, targetDB, outputBlat, params, ooc, tmpDir, nodesToUse, blatLogFile, mylistener):
    """Runs the blat on the SGE.
    """
    
    openLog = open(blatLogFile, "w")
    jobname = "s%s_%s" % (inputFasta[:5], random.randint(0,1000))
    openLog.write("Creating job %s\n" % (jobname))

    openLog.write("Starting blat job %s of %s against %s with params '%s'\n" % (jobname, inputFasta, targetDB, params))
    openLog.flush()
    if ooc and params.find("ooc")<0:
        params += " -ooc=%s " % (ooc)
    ag.startBlatJob(jobname, targetDB, inputFasta, outputBlat, params)

    # now monitor it until it finishes
    timeSinceUpdate = 0
    MINOR_WAIT = 300   # 5 minutes between "is it still running?" checks
    MAJOR_WAIT = 6000  # print status to the log file every hour
    while(True):
        time.sleep(MINOR_WAIT)
        if ag.isJobRunning(jobname):
            timeSinceUpdate += MINOR_WAIT
        else:
            break

        if timeSinceUpdate >= MAJOR_WAIT:
            sequenceUtils.blastOrBlatAnalysis(inputFasta, outputBlat, isBlast=False, logFile=openLog)
            openLog.flush()
            timeSinceUpdate = 0

    openLog.write("The job is no longer running...\n")
    percentComplete = sequenceUtils.blastOrBlatAnalysis(inputFasta, outputBlat, isBlast=False, logFile=openLog)

    if mylistener:
        mylistener.reportImportantInfo("The blat job %s is no longer running.  It appears to be %.2f%% complete." % (jobname, percentComplete))
    
    if percentComplete < 80.0:
        raise Exception("There appears to be an error.  The blat is no longer running, but blastOrBlatAnalysis determined the blat was only %.2f%% complete.  It may be that there weren't very many hits, or the cluster job may have failed." % (percentComplete))
    

def runStraightBlat(inputFasta, targetDB, outputBlat, params, oocFile, mylistener):
    """Runs the blat command directly (without a grid)."""
    if mylistener:
        mylistener.reportImportantInfo("Running blat.  Variable values:%s" % (getArgString(inspect.currentframe())))

    if oocFile and params.find("ooc")<0:
        params += " -ooc=%s " % (oocFile)
              
    cmd = "%s %s %s %s %s" % (config.PATH_TO_BLAT, targetDB, inputFasta, params, outputBlat)

    print "Running command:", cmd
    try:
        p = subprocess.Popen([cmd], shell=True)
        p.wait()
    except:
        raise

    percentComplete = sequenceUtils.blastOrBlatAnalysis(inputFasta, outputBlat, isBlast=False, logFile=None)

    if mylistener:
        mylistener.reportImportantInfo("The blat job is no longer running.  It appears to be %.2f%% complete." % (percentComplete))

def filterUsingBlat(inputFasta, blatResults, minPercent=None, minCoverage=None, minTotalPercent=None,
                    delimiter=None, pairedHits=None, pairedMisses=None, printHits=None, printMisses=None,
                    divideHalf=False, mylistener=None):
    """Given an input fasta and a set of blat results, filters the fasta so only those fasta
    records without hits (as defined by params) are in the outputFasta.

    Code is adapted from Peter's blatFilter, but adapted to handle single end or paired end and to allow
    pairedHits and pairedMisses
    """

    if mylistener:
        mylistener.reportImportantInfo("Filtering Blat output.  Variable values:%s" % (getArgString(inspect.currentframe())))

    if not divideHalf:
        countHits, countMisses = doFilterUsingBlat(inputFasta, blatResults, minPercent, minCoverage, minTotalPercent,
                    delimiter, pairedHits, pairedMisses, printHits, printMisses,
                    mylistener)

        if mylistener:
            mylistener.reportInfo("In blatFilter, there were %s hits and %s misses" % (countHits, countMisses))
    else:
        numDivisions = 4
        divPoints = determineDivisionPoints(inputFasta, numDivisions)
        #print "DivNum:", divNum
        fastaFiles = divideFile(inputFasta, divPoints, 'fasta')
        blatFiles = divideFile(blatResults, divPoints, 'blat')
        printHitsFiles = []
        printMissesFiles = []
        sumHits = 0
        sumMisses = 0
        for i in range(numDivisions):
            printHitsFile = printMissesFile = None
            if printHits:
                hitsRoot, ext = os.path.splitext(printHits)
                printHitFile = hitsRoot + "_%s" % (i+1) + ext
                printHitsFiles.append(printHitFile)
            if printMisses:
                missesRoot, ext = os.path.splitext(printMisses)
                printMissesFile = missesRoot + "_%s" % (i+1) + ext
                printMissesFiles.append(printMissesFile)

            countHits, countMisses = doFilterUsingBlat(fastaFiles[i], blatFiles[i], minPercent, minCoverage, minTotalPercent,
                    delimiter, pairedHits, pairedMisses, printHitsFile, printMissesFile,
                    mylistener)
            sumHits += countHits
            sumMisses += countMisses

            if mylistener:
                mylistener.reportInfo("For part %s of the blat filtering, there were %s hits and %s misses" % (i+1, countHits, countMisses))


        if printMisses:
            combineFiles ( printMissesFiles, printMisses)
        if printHits:
            combineFiles ( printHitsFiles, printHits)

        if mylistener:
            mylistener.reportInfo("In blatFilter, there were %s hits and %s misses" % (sumHits, sumMisses))

    return (0)



def doFilterUsingBlat(inputFasta, blatResults, minPercent=None, minCoverage=None, minTotalPercent=None,
                    delimiter=None, pairedHits=None, pairedMisses=None, printHits=None, printMisses=None,
                    mylistener=None):
    """Actually does the work."""

    if not minPercent:
        minPercent = 0
    if not minCoverage:
        minCoverage = 0
    if not minTotalPercent:
        minTotalPercent = 0

    if delimiter != None:
        pairedEndSuffixes = []

    hits = {}
    for line in open(blatResults):
        fields = line.split("\t")
        if len(fields) != 21:
            print "Found a header line:", line.strip()
            continue
        if len(line) < 10 or line.startswith("pslLayout") or line.startswith("match"):
            print "Found a header line:", line.strip()
            continue
                    
        query=fields[9]

        if delimiter != None and len(pairedEndSuffixes) != 2:
            parts = query.split(delimiter)
            if len(pairedEndSuffixes) == 0:
                #print "Split %s into %s using %s.  Appending %s to suffixes." % (query, parts, delimiter, parts[-1])
                pairedEndSuffixes.append(parts[-1])
            elif len(pairedEndSuffixes) == 1:
                if parts[-1] not in pairedEndSuffixes:
                    #print "Split %s into %s.  Appending %s to suffixes." % (query, parts, parts[-1])
                    pairedEndSuffixes.append(parts[-1])
        
        if query not in hits:
            totalMatches = int(fields[0]) + int(fields[2])
            sumBlockSizes = sum(map(float,filter(lambda x: x!= "", fields[18].split(","))))
            lenQuery = int(fields[10])
            #print "Done with the calculations.  totalMatches=%s, sumBlockSizes=%s, lenQuery=%s" % (totalMatches, sumBlockSizes, lenQuery)
            #print "sumBlockSizes/lenQuery=%s, minCoverage=%s, totalMatches/sumBlockSizes=%s, minPercent=%s, totalMatches/lenQuery=%s, minTotalPercent=%s" % (sumBlockSizes/lenQuery, minCoverage, totalMatches/sumBlockSizes, minPercent, float(totalMatches)/lenQuery, minTotalPercent)
            if sumBlockSizes/lenQuery >= minCoverage and totalMatches/sumBlockSizes >= minPercent and float(totalMatches)/lenQuery >= minTotalPercent:
                hits[query] = None
                
    print "Found %s hits in the blat results" % (len(hits))

    if len(hits) == 0:
        raise Exception("The blat PSL file had no hits.  Stopping IMSA execution.")

    if printHits:
        outHits = open(printHits, "w")
    if printMisses:
        outMisses = open(printMisses, "w")

    countHits = 0
    countMisses = 0
    #print
    #print
    for (title, sequence) in sequenceUtils.FastaIterator(inputFasta):
        read = title.strip().split()[0]
        #if countHits + countMisses < 10:
        #    print "'%s'" % (read)
        if delimiter != None:
            parts = read.split(delimiter)
            if parts[-1] == pairedEndSuffixes[0]:
                pairedEnd = parts[0] + delimiter + pairedEndSuffixes[1]
            else:
                pairedEnd = parts[0] + delimiter + pairedEndSuffixes[0]

            if pairedHits:
                if read in hits or pairedEnd in hits:
                    countHits += 1
                    if printHits:
                        outHits.write(">%s\n%s\n" % (read, sequence))
                else:
                    countMisses += 1
                    if printMisses:
                        outMisses.write(">%s\n%s\n" % (read, sequence))
            elif pairedMisses:
                if read in hits and pairedEnd in hits:
                    countHits += 1
                    if printHits:
                        outHits.write(">%s\n%s\n" % (read, sequence))
                else:
                    countMisses += 1
                    if printMisses:
                        outMisses.write(">%s\n%s\n" % (read, sequence))
        else:
            if read in hits:
                countHits += 1
                if printHits:
                    outHits.write(">%s\n%s\n" % (read, sequence))
            else:
                countMisses += 1
                if printMisses:
                    outMisses.write(">%s\n%s\n" % (read, sequence))

    return countHits, countMisses


def runBlast(inputFasta, targetDB, outputBlast, params, tmpDir, nodesToUse, blastLogFile, mylistener):
    """Runs gridBlast to run the blast if USE_GRID is on, otherwise runs the blast straight"""
    if mylistener:
        mylistener.reportImportantInfo("Running blast.  Variable values:%s" % (getArgString(inspect.currentframe())))

    # check things we always want to check
    if not os.path.exists(inputFasta):
        raise Exception("The input file for BLAST '%s' does not exist" % (inputFasta))
                                
    # submit grid or straight
    if config.USE_GRID:
        runGridBlast(inputFasta, targetDB, outputBlast, params, tmpDir, nodesToUse, blastLogFile, mylistener)
    else:
        runStraightBlast(inputFasta, targetDB, outputBlast, params, mylistener)


def runGridBlast(inputFasta, targetDB, outputBlast, params, tmpDir, nodesToUse, logFile, mylistener):
    """Uses the new ArronGrid to do the blast submission to the cluster!"""

    openLog = open(logFile, "w")
    jobname = "s%s_%s" % (inputFasta[:5], random.randint(0,1000))
    openLog.write("Creating job %s\n" % (jobname))

    if params.find("num_threads") < 0 and nodesToUse > 1:
        params += " -num_threads %s" % (nodesToUse)

    openLog.write("Starting job %s of %s against %s with params '%s'\n" % (jobname, inputFasta, targetDB, params))
    openLog.flush()
    ag.startBlastJob(jobname, targetDB, inputFasta, outputBlast, params)

    # now monitor it until it finishes
    timeSinceUpdate = 0
    MINOR_WAIT = 300   # 5 minutes between "is it still running?" checks
    MAJOR_WAIT = 6000  # print status to the log file every hour
    while(True):
        time.sleep(MINOR_WAIT)
        if ag.isJobRunning(jobname):
            timeSinceUpdate += MINOR_WAIT
        else:
            break

        if timeSinceUpdate >= MAJOR_WAIT:
            sequenceUtils.blastOrBlatAnalysis(inputFasta, outputBlast, isBlast=True, logFile=openLog)
            openLog.flush()
            timeSinceUpdate = 0

    openLog.write("The job is no longer running...\n")
    percentComplete = sequenceUtils.blastOrBlatAnalysis(inputFasta, outputBlast, isBlast=True, logFile=openLog)

    if mylistener:
        mylistener.reportImportantInfo("The blast job %s is no longer running.  It appears to be %.2f%% complete." % (jobname, percentComplete))
    
    if percentComplete < 80.0:
        raise Exception("There appears to be an error.  The blast is no longer running, but blastOrBlatAnalysis determined the blast was only %.2f%% complete.  It may be that there weren't very many hits, or the cluster job may have failed." % (percentComplete))
    

def runStraightBlast(inputFasta, targetDB, outputBlast, params, mylistener):
    if mylistener:
        mylistener.reportImportantInfo("Running blast.  Variable values:%s" % (getArgString(inspect.currentframe())))
              
    cmd = "%s -outfmt 6 %s -db %s -query %s -out %s" % (config.PATH_TO_BLASTN, params, targetDB, inputFasta, outputBlast)

    print "Running command:", cmd
    try:
        p = subprocess.Popen([cmd], shell=True)
        p.wait()
    except:
        raise

    percentComplete = sequenceUtils.blastOrBlatAnalysis(inputFasta, outputBlast, isBlast=True, logFile=None)

    if mylistener:
        mylistener.reportImportantInfo("The blast job is no longer running.  It appears to be %.2f%% complete." % (percentComplete))



def filterUsingBlast(inputFasta, blastResults, minPercent=None, minLength=None, minNumIds=None,
                     maxEval=None, minBitScore=None, delimiter='#0/', pairedHits=False, pairedMisses=True,
                     printHits=None, printMisses=None, mylistener=None):
    """Submits a job for doBlastFilter to the grid, using the doBlastFilter.py script in this directory.
    Roundabout way to do it, but the filter is memory intensive so it should be done on the grid.
    """

    if mylistener:
        mylistener.reportImportantInfo("Filtering Blast output '%s'.  Variable values:%s" % (blastResults,
                                                                     getArgString(inspect.currentframe())))



    doBlastFilter(inputFasta, blastResults, minPercent, minLength, minNumIds,
                  maxEval, minBitScore, delimiter, pairedHits, pairedMisses,
                  printHits, printMisses, mylistener)


    #cmd = "python %s/%s/doBlastFilter.py -f %s -b %s -d '%s'" % (constants.SRC_DIRECTORY,
    #    constants.PIPELINE_DIRECTORY,
    #    inputFasta, blastResults, delimiter)

    #if minPercent:
    #    cmd += " -p %s" % minPercent

    #if minLength:
    #    cmd += " -l %s" % minLength
    #if minNumIds:
    #    cmd += " -n %s" % minNumIds
    #if maxEval:
    #    cmd += " -e %s" % maxEval
    #if minBitScore:
    #    cmd += " -s %s" % minBitScore
    #if printHits:
    #    cmd += " -o %s" % printHits
    #elif printMisses:
    #    cmd += " -m %s" % printMisses

    #if pairedHits:
    #    cmd += " -u False"
    #else:
    #    cmd += " -u True"

    #print "Running command:", cmd
    #try:
    #    p = subprocess.Popen([cmd], shell=True)
    #    p.wait()
    #except:
    #    raise  


    # then read through the output to find the number of hits to report.  Wow.  Very roundabout!
    #fastaSize = 0
    #for (title, sequence) in sequenceUtils.FastaIterator(inputFasta):
    #    fastaSize += 1
    #numHits = 0
    #numMisses = 0
    #if printHits:
    #    for (title, sequence) in sequenceUtils.FastaIterator(printHits):
    #        numHits += 1
    #    numMisses = fastaSize - numHits
    #elif printMisses:
    #    for (title, sequence) in sequenceUtils.FastaIterator(printMisses):
    #        numMisses +=1
    #    numHits = fastaSize - numMisses

    if mylistener:
        mylistener.reportInfo("PAIRED BLAST FILTERING IS TURNED OFF!")
        #mylistener.reportInfo("In blastFilter, there were %s hits and %s misses." % (numHits, numMisses))
    return (0)



def doBlastFilter(inputFasta, blastResults, minPercent=None, minLength=None, minNumIds=None,
                  maxEval=None, minBitScore=None, delimiter='#0/', pairedHits=False, pairedMisses=True,
                  printHits=None, printMisses=None, mylistener=None):
    """Given an input fasta and a set of blast results, filters the fasta, putting the hits in the file
    defined by printHits (if not none) and the misses in the file defined by printMisses (if not none).

    Code is adapted from blastFilter.py by Peter
    """

    #print getArgString(inspect.currentframe())

    delimiter = None
    
    if delimiter != None:
        pairedEndSuffixes = []

    if not printHits and not printMisses:
        raise Exception("You must specify an output, either as printHits or printMisses, in filterUsingBlast")

    
    curm8 = open(blastResults)
    line = curm8.readline()
    hits = {}
    while line:
        if not line.startswith("#") and len(line.split("\t")) == 12:
            isHit = False
            fields = line.split("\t")
            query = fields[0]

            if ((minPercent == None or float(fields[2]) >= minPercent) and (minLength == None or int(fields[3]) >= minLength) and \
                (minNumIds == None or ((float(fields[3])-float(fields[4])-float(fields[5]))>=minNumIds)) and \
                (maxEval == None or float(fields[10]) < maxEval) and (minBitScore == None or float(fields[11]) > minBitScore)):
                isHit = True

            if delimiter != None and len(pairedEndSuffixes) != 2:
                parts = query.split(delimiter)
                if len(pairedEndSuffixes) == 0:
                    pairedEndSuffixes.append(parts[-1])
                elif len(pairedEndSuffixes) == 1:
                    if parts[-1] not in pairedEndSuffixes:
                        pairedEndSuffixes.append(parts[-1])

            if isHit:
                hits[query] = None

        line = curm8.readline()

    curm8.close()
    pairedEndSuffixes = ['1', '2']

    fastaIn = sequenceUtils.FastaIterator(inputFasta)
    if printHits:
        queryHits = open(printHits, "w")
    if printMisses:
        queryMisses = open(printMisses, "w")

    countHits = 0
    countMisses = 0
    for (title, sequence) in fastaIn:
        read = title.strip()
        if delimiter != None:
            parts = read.split(delimiter)
            if parts[-1] == pairedEndSuffixes[0]:
                pairedEnd = parts[0] + delimiter + pairedEndSuffixes[1]
            else:
                pairedEnd = parts[0] + delimiter + pairedEndSuffixes[0]

            if pairedHits:
                if read in hits or pairedEnd in hits:
                    countHits += 1
                    if printHits:
                        queryHits.write(">%s\n%s\n" % (read, sequence))
                else:
                    countMisses += 1
                    if printMisses:
                        queryMisses.write(">%s\n%s\n" % (read, sequence))
            elif pairedMisses:
                if read in hits and pairedEnd in hits:
                    countHits += 1
                    if printHits:
                        queryHits.write(">%s\n%s\n" % (read, sequence))
                else:
                    countMisses += 1
                    if printMisses:
                        queryMisses.write(">%s\n%s\n" % (read, sequence))
        else:
            if read in hits:
                countHits += 1
                if printHits:
                    queryHits.write(">%s\n%s\n" % (read, sequence))
            else:
                countMisses += 1
                if printMisses:
                    queryMisses.write(">%s\n%s\n" % (read, sequence))

    if mylistener:
        mylistener.reportInfo("In blastFilter, there were %s hits and %s misses." % (countHits, countMisses))
    return (0)



def runLZWfilter(inputFasta, outputFasta, lzwThreshold, mylistener):

    if mylistener:
        mylistener.reportImportantInfo("Running LZW filter.  Variable values:%s" % (getArgString(inspect.currentframe())))

    lzwThreshold = float(lzwThreshold)
    out = open(outputFasta, "w")
    countPasses = 0
    countFails = 0
    for (title, sequence) in sequenceUtils.FastaIterator(inputFasta):
        curRatio = sequenceUtils.LZWratio(sequence)
        if curRatio >= lzwThreshold:
            countPasses +=1
            out.write(">%s\n%s\n" % (title, sequence))
        else:
            countFails += 1

    if mylistener:
        mylistener.reportImportantInfo("LZW filter at level %.2f had %s passes and %s fails." % (lzwThreshold, countPasses, countFails))

    

def convertBackToFastq(originalFastq, inputFasta, outputFastq, mylistener):
    if mylistener:
        mylistener.reportImportantInfo("Using '%s' to convert '%s' to fastq file '%s'" % (originalFastq, inputFasta, outputFastq))
    open(outputFasta, "w")

def createLogTable(logFile, outputName=None):
    """Given a log file, creates a table of the hits/misses from blast/blat."""
    out = None
    if outputName:
        out = open(outputName, "w")

    curOp = None
    curStart = None
    curDB = None
    curParams = []

    for line in open(logFile):
        if line.find("Running ") >= 0:
            curParams = []
            pieces = line.split()
            curOp = pieces[3][:-1]
            curStart = convertTime(pieces[0], pieces[1])
        elif line.find("minPercent")>=0 or line.find("minLength")>=0 or line.find("minNumIds")>=0 or line.find("maxEval")>=0 or line.find("minBitScore")>=0 or line.find("minCoverage")>=0 or line.find("minTotalPercent")>=0:
            if not line.find("None")>=0:
                curParams.append(line.strip())
        elif line.find("targetDB")>= 0:
            curDB=line.split()[2].split("/")[-1]
        elif line.find("In blastFilter")>=0 or line.find("In blatFilter")>0 or line.find("In bowtieFilter")>0:
            pieces = line.split()
            numHits = pieces[6]
            numMisses = pieces[9]
            finishTime = convertTime(pieces[0], pieces[1])
            try:
                runTime = finishTime-curStart
            except:
                runTime = "NA"
            try:
                percentHits = float(numHits) * 100 / ( float(numHits) + float(numMisses))
            except ZeroDivisionError:
                percentHits = 0.0
            towrite = "%s\t%s\t%s\t%s\t%s\t%.2f\t%s" % (curOp, curDB, runTime, numHits, numMisses,
                                                        percentHits, ", ".join(curParams))
            if out:
                out.write(towrite)
                out.write("\n")
            else:
                print towrite
            curOp = None
            curDB = None
            curStart = None
            curParams = []
        
def convertTime(yearStr, timeStr):
    #print yearStr, yearStr.split("//"), yearStr.split("/")
    year, month, day = yearStr.split("/")
    hour, minute, second = timeStr.split(":")[:3]
    return datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))


class taxNode:

    def __init__ (self, sname, id, rank):
        self.sname = sname
        self.id = id
        self.rank = rank
        self.count = 0
        self.totalCount = 0

    def __repr__ (self):

        width=3.0
        if self.totalCount > 0:
            ratio = float(self.count) / self.totalCount
            if ratio > config.MAX_NODE_SIZE_CUTOFF:
                width = config.MAX_NODE_SIZE
            else:
                width = (ratio * (1/config.MAX_NODE_SIZE_CUTOFF) * (config.MAX_NODE_SIZE-config.MIN_NODE_SIZE)) + config.MIN_NODE_SIZE

            #print "%s\t\t%.3f\t\t%.3f" % (self.totalCount, ratio, width)

        height = width / config.WIDTH_TO_HEIGHT
        fontcolor = "black"
        return '%s [color=%s,fixedsize=true,width=%.1f,height=%.1f,shape=ellipse,style=bold,fontcolor=%s,label="%s\\n%.0f"];' % (
            self.id, config.RANK_TO_COLOR[self.rank], width, height, fontcolor, self.sname, self.count)


class taxEdge:

    def __init__ (self, name1, name2):
        self.name1 = name1
        self.name2 = name2

    def __repr__ (self):

        return "%s -> %s;" % (self.name1, self.name2)

def reportTaxonomy(blastInput, printTargets=False, outputPrefix=None, dotPrefix=None, dotLimit=1.0):
    """Given a blast input file, reports the taxonomy lineage for the hits.  By default the report
    will be printed to standard output, but if an outputPrefix is given then output files, one each
    for species, genus, family and division, will be created.
    The Blast input file should have already been run through a Best Blast filter if you want to only
    count unique hits.  Can have the extra counts column added, as per blastUtils.keepAllBestHits."""
    
    targets = readBlastIntoTargetDict(blastInput)

    if printTargets:
        if outputPrefix:
            outTargets = open(outputPrefix+"targets.txt", "w")

        sortedTargets = sorted(targets.iteritems(), key=operator.itemgetter(1), reverse=True)
        for (k, v) in sortedTargets:
            if outputPrefix:
                outTargets.write("%s=%s\n" % (k, v))
            else:
                print k, "=", v
        if outputPrefix:
            outTargets.write("\n\n")
        else:
            print
            
    # change full target to gi
    gis = {}
    for k in targets.keys():
        gi = k.split("|")[1]
        gis[gi] = targets[k]

    taxonomy = {}
    print "Getting taxonomy ids from local file for %s gis" % (len(gis))
    for line in open(config.BLAST_TAX_DB):
        (gi, taxid) = line.split()
        if gis.has_key(gi):
            taxonomy[gi] = int(taxid)

    print "Done getting taxonomy ids, got %s ids" % (len(taxonomy))
    
    countTaxIds = {}
    for (gi, count) in gis.iteritems():
        try:
            taxid = taxonomy[gi]
            if not countTaxIds.has_key(taxid):
                countTaxIds[taxid] = 0
            countTaxIds[taxid] += count
        except KeyError:
            print "Unable to find gi %s in taxonomy results.  Ignoring this gi" % (gi)

    # get the scientific name for each taxa
    print "Accessing NCBI for %s taxonomies" % len(countTaxIds)
    names = eutilsWrapper.getTaxa(countTaxIds.keys())
    print "Got %s taxonomies from NCBI" % (len(names))

    # Build all the pretty reports
    nodes = {}
    edges = {}


    notFound = 0
    for (taxId, count) in countTaxIds.iteritems():
        try:
            fullTax = names[taxId]
        except:
            notFound += 1
            print "Unable to find taxId %s in NCBI results.  (%s not found so far)" % (taxId, notFound)
            continue

        species = fullTax.getSpecies()
        genus = fullTax.getGenus()
        family = fullTax.getFamily()
        division = fullTax.getDivision()

        if not nodes.has_key(division):
            nodes[division] = {}
            nodes[division][division] = taxNode(division, division, "division")
        if not edges.has_key(division):
            edges[division] = {}

        nodes[division][division].count += count

        if species and genus:
            if not edges[division].has_key( (genus.taxId,species.taxId)):
                edges[division][(genus.taxId,species.taxId)] = taxEdge(genus.taxId, species.taxId)
        if genus and family:
            if not edges[division].has_key( (family.taxId, genus.taxId)):
                edges[division][(family.taxId,genus.taxId)] = taxEdge(family.taxId, genus.taxId)
        if family:
            if not edges[division].has_key( (division, family.taxId)):
                edges[division][(division,family.taxId)] = taxEdge(division, family.taxId)

        if species:
            if not nodes[division].has_key(species.taxId):
                nodes[division][species.taxId] = taxNode(species.sname, species.taxId, "species")
            nodes[division][species.taxId].count += count

        if genus:
            if not nodes[division].has_key(genus.taxId):
                nodes[division][genus.taxId] = taxNode(genus.sname, genus.taxId, "genus")
            nodes[division][genus.taxId].count += count

        if family:
            if not nodes[division].has_key(family.taxId):
                nodes[division][family.taxId] = taxNode(family.sname, family.taxId, "family")
            nodes[division][family.taxId].count += count

    # organize the data for easy text report
    speciesD = {}
    genusD = {}
    familyD = {}
    divisionD = {}

    for (division, nodeSet) in nodes.iteritems():
        divisionD[(division,division)] = nodes[division][division].count
        for (taxId, node) in nodeSet.iteritems():
            if node.rank == "species":
                speciesD[(taxId, node.sname)] = node.count
            elif node.rank == "genus":
                genusD[(taxId, node.sname)] = node.count
            elif node.rank == "family":
                familyD[(taxId, node.sname)] = node.count

    if outputPrefix:
        outSpecies = open(outputPrefix+"species.txt", "w")
        outGenus = open(outputPrefix+"genus.txt", "w")
        outFamily = open(outputPrefix+"family.txt", "w")
        outDivision = open(outputPrefix+"division.txt", "w")

        outSpecies.write("%s\t%s\t%s\n" % ("Species ID", "Scientific Name", "Count"))
        outSpecies.write(getTaxDictionaryStr(speciesD))
        outGenus.write("%s\t%s\t%s\n" % ("Genus ID", "Scientific Name", "Count"))
        outGenus.write(getTaxDictionaryStr(genusD))
        outFamily.write("%s\t%s\t%s\n" % ("Family ID", "Scientific Name", "Count"))
        outFamily.write(getTaxDictionaryStr(familyD))
        outDivision.write("%s\t%s\t%s\n" % ("Division ID", "Scientific Name", "Count"))
        outDivision.write(getTaxDictionaryStr(divisionD))
    else:
        print ("%s\t%s\t%s" % ("Species ID", "Scientific Name", "Count"))
        print (getTaxDictionaryStr(speciesD))
        print ("%s\t%s\t%s" % ("Genus ID", "Scientific Name", "Count"))
        print (getTaxDictionaryStr(genusD))
        print ("%s\t%s\t%s" % ("Family ID", "Scientific Name", "Count"))
        print (getTaxDictionaryStr(familyD))
        print ("%s\t%s\t%s" % ("Division ID", "Scientific Name", "Count"))
        print (getTaxDictionaryStr(divisionD))


    if dotPrefix:
        for division in nodes.keys():
            makeDotFile(nodes[division], edges[division], dotPrefix+"_"+division+".txt", dotLimit, nodes[division][division].count)


def makeDotFile(nodeDict, edgesDict, outputName, dotLimit, totalCount):

    out = open(outputName, "w")
    out.write("digraph TaxonomyReport{\n\n")
    includedNodes = {}
    for (taxId, node) in nodeDict.iteritems():
        if node.count >= dotLimit:
            includedNodes[taxId] = 1
            node.totalCount = totalCount
            out.write("%s\n" % node)
    for (taxId, edge) in edgesDict.iteritems():
        if includedNodes.has_key(edge.name1) and includedNodes.has_key(edge.name2):
            out.write("%s\n" % edge)

    out.write("overlap=false\n}\n")

def readBlastIntoTargetDict(blastInput, keepQueries=False):
    """Reads the blast to create a dictionary where the key is the target and the value
    is the count of the number of hits for that target (1 if it's unique, could be 0.5 if the
    query hits two targets equally and the file has been run through the keepAllBlastHits
    filter to add the counts column.).  If keepQueries is true than the queries that hit
    the target are saved instead of the count value."""
    
    targets = {}
    for line in open(blastInput):
        pieces = line.split()
        if len(pieces) < 12 or len(pieces) > 13:
            print "Ignoring line, it has %s pieces:" % len(pieces)
            print line
            continue
        elif len(pieces) == 12:
            hitVal = 1
        elif len(pieces) == 13:
            hitVal = float(pieces[12])
            
        target = pieces[1]
        if not targets.has_key(target):
            if keepQueries:
                targets[target] = []
            else:
                targets[target] = 0
        if keepQueries:
            targets[target].append(pieces[0])
        else:
            targets[target] += hitVal

    return targets
        
def getTaxDictionaryStr(d):
    retval = ""
    sortedTax = sorted(d.iteritems(), key=operator.itemgetter(1), reverse=True)
    #print sortedTax
    
    for ((taxId, sname), count) in sortedTax:
        retval += "%s\t%s\t%s\n" % (taxId, sname, count)

    return retval

def getFastaForTaxonomy(taxIdList, blastInputGlob, fastaInputGlob,  fastaOutput, printLineageOfHits=False):
    """Given a fasta file and the blast results of that fasta file against nt (or another NCBI database
    with gi-formatted headings), this function will pull all the reads that match against any of the
    tax ids in the taxIdList and put them in the output file.  The tax id can be at any level (species,
    genus, etc).

    Hint:  The viral superkingdom taxId is 10239, so use that taxid to get all viruses.
    """

    allTargets = {}
    for blastInput in blastInputGlob:
        targets = readBlastIntoTargetDict(blastInput, keepQueries=True)
        for k,v in targets.iteritems():
            if not allTargets.has_key(k):
                allTargets[k] = []
            allTargets[k].extend(v)

    # change full target to gi
    gis = {}
    for k in targets.keys():
        gi = k.split("|")[1]
        gis[gi] = targets[k]

    taxonomy = {}
    usedTaxonomies = set()
    print "Getting taxonomy ids from local file for %s gis" % (len(gis))
    for line in open(BLAST_TAX_DB):
        (gi, taxid) = line.split()
        if gis.has_key(gi):
            usedTaxonomies.add(taxid)
            taxonomy[gi] = int(taxid)

    print "Done getting taxonomy ids.  Found %s gis, resulting in %s taxonomy ids" % (len(taxonomy), len(usedTaxonomies))

    # get the scientific name for each taxa
    print "Accessing NCBI for %s taxonomies" % len(usedTaxonomies)
    names = eutilsWrapper.getTaxa(usedTaxonomies)
    print "Got %s taxonomies from NCBI" % (len(names))

    print "Determining which reads to keep"
    queriesToKeep = set()
    #for q in queriesToKeep:
    #    print q
        
    for target in targets.keys():
        gi = target.split("|")[1]
        try:
            taxid = taxonomy[gi]
        except KeyError:
            print "Unable to find taxonomy for %s.  Skipping." % (gi)
            continue
        try:
            fullTax = names[taxid]
        except KeyError:
            print "Unable to find full tax for %s.  Skipping." % (taxid)
            continue

        for taxObj in fullTax.lineage:
            if taxObj.taxId in taxIdList:
                if printLineageOfHits:
                    print "found %s reads matching %s" % (len(targets[target]), target)
                    for taxObj in fullTax.lineage:
                        print "       %s: %s (%s)" % (taxObj.taxId, taxObj.sname, taxObj.rank)
                queriesToKeep = queriesToKeep.union(gis[gi])

    print "Keeping %s reads" % (len(queriesToKeep))
    for x in queriesToKeep:
        print "    %s" % (x)

    out = open(fastaOutput, "w")
    for fastaInput in fastaInputGlob:
        print "Reading fasta %s" % (fastaInput)
        countFound = 0
        for (title, record) in sequenceUtils.FastaIterator(fastaInput):
            if title in queriesToKeep:
                out.write(">%s\n%s\n" % (title, record))
                countFound += 1

        print "  Found %s reads in that fasta" % (countFound)
