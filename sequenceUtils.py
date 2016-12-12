from types import StringTypes, StringType
#import ctypes
import string


# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)

# Rather than do this through AOS, just do it ourselves.  Maybe a bit slower in Python, but
# not that much slower
#AOS_LOADED = False
#try:
#    __aos = ctypes.cdll.LoadLibrary("libaos.so.1.1.1")
#    __aos.energy.restype=ctypes.c_float
#    AOS_LOADED = True
#except Exception as inst:
#    print "AOS Loading error:", inst
#    print "Unable to load AOS.  The only functionality that uses AOS is computing LZW ratios."

#====================================================================================
#  General Utils
#====================================================================================

def median(list):
    """Returns the median of a list.  Makes a copy before sorting so the original list is not changed."""
    if len(list) == 0:
        return 0
    
    listCopy = list[:]
    listCopy.sort()
    if len(listCopy) % 2 == 0:
        return (listCopy[ (len(list)-1) / 2] + listCopy[ len(list)/2 ]) / 2.0
    else:
        return listCopy[ (len(list)/2) ]
    
def mean(list):
    """Returns the mean of the list."""
    
    if len(list) == 0:
        return None
    
    try:
        total = 0
        numItems = 0
        # have to write our own sum to ignore None values
        for x in list:
            if x:
                total += x
                numItems += 1
    except TypeError:
        raise Exception("Unable to take the sum of the list in mathUtils.mean for list: %s" % (list))
    
    return total / float(numItems)

def FastqIterator(fh):
    """return an iterator of Records found in file handle, fh.
    Original code from Kael and Dale.
    """
    def readTotitle(fh, titleChar):
        """returns a tuple ([lines before the next title line], next tile line)
        """
        preLines = []
        while True:
            l = fh.readline().strip()
            if l.startswith(titleChar):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)

    if type(fh) in StringTypes:
        fh = file(fh)
    
    preLines,nextTitleLine =readTotitle(fh,'@')

    while nextTitleLine != None:
        seqTitle = nextTitleLine[1:].rstrip()
        preLines,nextTitleLine=readTotitle(fh,'+')
        qualTitle = nextTitleLine[1:].rstrip()
        if len(qualTitle.strip()) > 0 and seqTitle != qualTitle:
            print seqTitle
            print preLines
            print qualTitle
            raise Exception("Error in parsing: @title sequence entry must be immediately followed by corresponding +title quality entry.")
        seqLines = preLines
        qualLines = []
        for i in range(len(seqLines)): # Quality characters should be the same length as the sequence
            qualLines.append( fh.readline().strip() )
        
        preLines,nextTitleLine=readTotitle(fh,'@')

        yield (seqTitle, ''.join(seqLines), ''.join(qualLines))


def FastaIterator(fh):
    """return an iterator of Records found in file handle, fh.
        Original code from Kael and Dale.
    """
    def readTotitle(fh):
        """returns a tuple ([lines before the next title line], next tile line)
        """
        preLines = []
        while True:
            l = fh.readline().strip()
            if l.startswith('>'):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)


    if type(fh) in StringTypes:
        fh = file(fh)
    
    preLines,nextTitleLine =readTotitle(fh)

    while nextTitleLine != None:
        title = nextTitleLine[1:].rstrip()
        preLines,nextTitleLine=readTotitle(fh)
        yield (title,''.join(map(lambda x: x.rstrip(),preLines)))


def trimFasta(inputFasta, outputFasta, trimFront=0, trimRear=0):
    if trimFront == 0 and trimRear == 0:
        raise Exception("You must specify a trim value for either the front or the rear.")

    out = open(outputFasta, "w")
    for (title, record) in FastaIterator(inputFasta):
        if trimFront > 0:
            trimmed = record[trimFront:]
        else:
            trimmed = record
        if trimRear > 0:
            trimmed = trimmed[:-trimRear]
        out.write(">%s\n%s\n" % (title, trimmed))


def trimFastq(inputFastq, outputFastq, trimFront=0, trimRear=0):
    if trimFront == 0 and trimRear == 0:
        raise Exception("You must specify a trim value for either trimFront or trimRear")

    out = open(outputFastq, "w")
    for (title, sequence, quality) in FastqIterator(inputFastq):
        if trimFront > 0:
            trimSeq = sequence[trimFront:]
            trimQual = quality[trimFront:]
        else:
            trimSeq = sequence
            trimQual = quality
        if trimRear > 0:
            trimSeq = trimSeq[:-trimRear]
            trimQual = trimQual[:-trimRear]
        if len(trimSeq) < 1:
            raise Exception("In trimFastq, the length of the sequence is zero.  Original length=%s, trimFront=%s, trimRear=%s" % (len(sequence), trimFront, trimRear))
        out.write("@%s\n%s\n+%s\n%s\n" % (title, trimSeq, title, trimQual))


def reverse(seq):
    """return a copy of the reversed seq
    """
    rv = list(seq)
    rv.reverse()
    return ''.join(rv)

def reverseComplement(seq):
    """return the reverse complement of a sequence.
    Original code from Kael and Dale."""
    seq=seq.upper()
    # complement
    compl = complement(seq)
    # reverse
    return compl[::-1]

def complement(seq,transl=None):
    """Return complement of seq.
    Original code from Kael and Dale.
    """
    transl = string.maketrans('aAcCgGtTnNxX-\t\n ','tTgGcCaAnNxX-\t\n ')
    compl = seq.translate(transl)
    return compl

def aosEnergy(s1,s2):
    """Return melting energy for association os s1/s2.
    """

    if type(s1) != StringType or type(s2) != StringType:
        raise Exception("s1  and s2 must be strings")

    if len(s2) != len(s1):
        raise Exception("s1 and s2 are different lengths.")

    s2=reverse(s2)

    return __aos.energy(ctypes.c_char_p(s1),ctypes.c_char_p(s2),ctypes.c_int(0))


def LZWsize(seq):
    return len(LZWcompress(seq))

def LZWratio(seq):
    if len(seq) == 0:
        raise Exception("Divide by zero error: length of sequence in LZW ratio is 0.")
    return float(LZWsize(seq)) / len(seq)


# compress and decompress are LZW implementations copied from http://rosettacode.org/wiki/LZW_compression
def LZWcompress(uncompressed):
    """Compress a string to a list of output symbols."""
 
    # Build the dictionary.
    dict_size = 256
    dictionary = dict((chr(i), chr(i)) for i in xrange(dict_size))
    # in Python 3: dictionary = {chr(i): chr(i) for i in range(dict_size)}
 
    w = ""
    result = []
    for c in uncompressed:
        wc = w + c
        if wc in dictionary:
            w = wc
        else:
            result.append(dictionary[w])
            # Add wc to the dictionary.
            dictionary[wc] = dict_size
            dict_size += 1
            w = c
 
    # Output the code for w.
    if w:
        result.append(dictionary[w])
    return result
 
 
def LZWdecompress(compressed):
    """Decompress a list of output ks to a string."""
 
    # Build the dictionary.
    dict_size = 256
    dictionary = dict((chr(i), chr(i)) for i in xrange(dict_size))
    # in Python 3: dictionary = {chr(i): chr(i) for i in range(dict_size)}
 
    w = result = compressed.pop(0)
    for k in compressed:
        if k in dictionary:
            entry = dictionary[k]
        elif k == dict_size:
            entry = w + w[0]
        else:
            raise ValueError('Bad compressed k: %s' % k)
        result += entry
 
        # Add w+entry[0] to the dictionary.
        dictionary[dict_size] = w + entry[0]
        dict_size += 1
 
        w = entry
    return result


def sam2Fastq(inputSam, outputFastq, outputPaired=False, checkPairs=False,
              fixReverse=False, printStats=False, printStatsDuring=False):
    """Converts a SAM file into a FASTQ file.
    If fixReverse is true then reads that have been reverse complemented in
    the alignment will be returned to their original sequence.
    If printStats is true then some
    basic statistics are printed about the number of matching pairs, etc."""

    if outputPaired:
        out1 = open(outputFastq+"_1.txt", "w")
        out2 = open(outputFastq+"_2.txt", "w")
    else:
        out = open(outputFastq+".txt", "w")

    countProper = 0
    countThisNotOther = 0
    countOtherNotThis = 0
    countNeither = 0
    countBadBitwiseFlag = 0
    countSortedPair = 0
    countBadPair = 0

    countTotal = 0
    countPaired = 0
    countFirstRead = 0
    countSecondRead = 0

    lastReadName = None
    # whether the last read is the first (as opposed to second) of a pair
    lastReadisFirst = None
    lastReadLine = None
    
    if printStatsDuring:
        print "\t".join(["total", "paired", "firstRead", "secondRead", "proper", "thisNotOther", "otherNotThis", "neither"])

    for line in open(inputSam):
        if line.startswith("@"):
            continue
        if len(line) < 2:
            continue

        try:
            [qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual] = line.split()[:11]
            flag = int(flag)
        except:
            print "Unable to parse this line, it doesn't appear to be SAM format:"
            print line
            print line.split()
            continue

        pairedRead = flag & 1 == 1
        if pairedRead:
            isFirstRead = flag & 64 == 64
            isSecondRead = flag & 128 == 128
            if (isFirstRead and isSecondRead):
                raise Exception("A read cannot be the first read AND the second read.  See %s" % (qname))
            if (not isFirstRead and not isSecondRead):
                raise Exception("A read cannot be paired but neither the first read nor second read. See %s." % (qname))

        # check pairs BEFORE adding \1 and \2 to the names
        okToPrint = False
        if checkPairs:
            if not pairedRead:
                countBadPair += 1
                print ("Error!  The paired flag is false")
            else:
                if lastReadName:
                    if qname != lastReadName:
                        countBadPair += 1
                        print "ERROR! The last read was %s and this read is %s" % (lastReadName, qname)
                        lastReadName = qname
                        lastReadisFirst = isFirstRead
                        lastReadLine = None
                    elif not ( (isFirstRead and not lastReadIsFirst) or (not isFirstRead and lastReadIsFirst)):
                        countBadPair += 1
                        print "ERROR!  The last read pair value is %s and isFirstRead=%s" % (lastReadisFirst, isFirstRead)
                        lastReadName = qname
                        lastRaedIsFirst = isFirstRead
                        lastReadLine = None
                    else:
                        #print "OK TO PRINT"
                        okToPrint = True
                else:
                    #print "Setting lastReadName because blank"
                    lastReadName = qname
                    lastReadIsFirst = isFirstRead
                    lastReadLine = ""
                    
        if pairedRead:
            if isFirstRead:
                qname = qname + "/1"
            elif isSecondRead:
                qname = qname + "/2"
            else:
                raise Exception("Paired Read, yet neither first read nor second read is true.")


        if fixReverse or printStats:

            pairMapped = flag & 2 == 2
            readUnmapped = flag & 4 == 4
            mateUnmapped = flag & 8 == 8
            readReversed = flag & 16 == 16
            mateReversed = flag & 32 == 32

            if pairMapped and readUnmapped:
                #raise Exception("Disagreement on the bitwise flag!  Flag=%s, which says pairMapped and readUnmapped!" % (flag))
                #print "Bad read %s, Disagreement on the bitwise flag!  Flag=%s, which says pairMapped and readUnmapped!" % (qname, flag)
                countBadBitwiseFlag += 1
                
            if pairMapped and mateUnmapped:
                #raise Exception("Disagreement on the bitwise flag!  Flag=%s, which says pairMapped and mateUnmapped!" % (flag))
                #print "Bad read %s, Disagreement on the bitwise flag!  Flag=%s, which says pairMapped and readUnmapped!" % (qname, flag)
                countBadBitwiseFlag += 1

            if fixReverse and readReversed:
                seq = reverseComplement(seq)
                qual = reverse(qual)

            if printStats:
                countTotal += 1
                if pairMapped:
                    countProper += 1
                if readUnmapped and mateUnmapped:
                    countNeither += 1
                if readUnmapped and not mateUnmapped:
                    countOtherNotThis += 1
                if not readUnmapped and mateUnmapped:
                    countThisNotOther += 1
                if isFirstRead:
                    countFirstRead += 1
                if isSecondRead:
                    countSecondRead += 1
                if pairedRead:
                    countPaired += 1
                if printStatsDuring and countTotal % 1000000==0:
                    print "\t".join([str(x) for x in [countTotal, countPaired, countFirstRead, countSecondRead,
                                                      countProper, countThisNotOther, countOtherNotThis, countNeither,
                                                      countBadBitwiseFlag]])
                                        


        if pairedRead and not checkPairs:
            if firstRead:
                out1.write("@%s\n%s\n+\n%s\n" % (qname, seq, qual))
            else:
                out2.write("@%s\n%s\n+\n%s\n" % (qname, seq, qual))
        elif pairedRead and checkPairs:
            if okToPrint:
                #print "PRINTING!"
                countSortedPair += 1
                if isFirstRead:
                    out1.write("@%s\n%s\n+\n%s\n" % (qname, seq, qual))
                    out2.write(lastReadLine)
                else:
                    out2.write("@%s\n%s\n+\n%s\n" % (qname, seq, qual))
                    out1.write(lastReadLine)
                lastReadName = None
                lastReadisFirst = None
                lastReadLine = None
            else:
                #print "Skipping Print"
                lastReadLine = "@%s\n%s\n+\n%s\n" % (qname, seq, qual)
                    
        else:
            out.write("@%s\n%s\n+\n%s\n" % (qname, seq, qual))

    if lastReadName:
        countBadPair += 1

    if printStats: 
        print "There were %s reads.  Of these, %s were paired reads, with %s first reads and %s second reads." % (countTotal, countPaired, countFirstRead, countSecondRead)
        print "There were %s paired alignments, %s where one aligned and the other didn't (and %s the other way) and %s where neither read aligned." % (countProper,
                                                                                                             countThisNotOther, countOtherNotThis, countNeither)
        print "There were %s bad bitwise flag, %s sorted pairs and %s bad pairs" % (countBadBitwiseFlag, countSortedPair, countBadPair)

def fastqToOneLine(inputFastq, outputOneLine):

    out = open(outputOneLine, "w")

    count = 0
    countRecords = 0
    for line in open(inputFastq):
        line = line.strip()
        count = count + 1

        out.write(line)
        if count % 4 == 0:
            countRecords += 1
            out.write("\n")
        else:
            out.write("\t")

    print "There were %s records in %s lines (check %s = %s)." % (countRecords, count, countRecords, count/4.0)

def oneLineToFastq(inputOneLine, outputFastq):

    out = open(outputFastq, "w")

    count = 0
    for line in open(inputOneLine):
        pieces = line.split()
        out.write("\n".join(pieces))
        out.write("\n")
        count += 1

    print "There were %s records converted." % (count)

def subtractFastas(fasta1, fasta2,  fasta3):
    """Creates a file fasta3 containing the records in fasta1 but not fasta2, i.e. fasta1 - fasta2 = fasta3."""

    d2 = {}
    for (title, sequence) in FastaIterator(fasta2):
        d2[title] = 1

    countFasta1 = 0
    countFasta2 = len(d2)
    countFasta3 = 0
    out = open(fasta3, "w")
    for (title, sequence) in FastaIterator(fasta1):
        countFasta1 += 1
        if not d2.has_key(title):
            countFasta3 += 1
            out.write(">%s\n%s\n" % (title, sequence))
            

    print "Created %s.  %s - %s = %s records in %s" % (fasta3, countFasta1, countFasta2, countFasta3, fasta3)

def graphQualityPerPosition(inputFastq):
    """Creates a histogram counting the quality value per position in the input fastq file."""

    histD = {}

    count = 0
    for (titleStr, seqStr, qualityStr) in FastqIterator(inputFastq):
        count += 1
        if count < 2000000:
            continue
        if count > 3000000:
            break

        qInts = convertQualityStr(qualityStr)
        for i in range(len(qInts)):
            q = qInts[i]
            if q < 0 or q > 40:
               raise Exception("Invalid quality value %s at position %s of %s" % (q, i, qualityStr))
            if not histD.has_key(i):
                histD[i] = [0]*41
            histD[i][q] += 1

    print "Histogram of quality score per position"
    allk = histD.keys()
    allk.sort()
    for k in allk:
        print "%s|" % k, "|".join(str(x) for x in histD[k])

def convertQualityStr(strVersion, offset=64):
    """Converts a quality string to a list of integers.
    an offset of 33 is Phred style and an offet of 64 is Solexa style"""
    return map( lambda x: ord(x)-offset, strVersion )


def blastOrBlatAnalysis(inputFasta, outputAlignment, isBlast=True, logFile=None):
    dResults = {}
    countResults = 0
    lastRead = ""
    for line in open(outputAlignment):
        if len(line) < 10 or line.startswith("pslLayout") or line.startswith("match"):
            continue
        countResults += 1
        try:
            if isBlast:
                title = line.split()[0]
            else:
                title = line.split()[9]
        except IndexError:
            #ignore this line, it may be malformatted because the file is partial
            continue
        if not dResults.has_key(title):
            dResults[title] = 0
        dResults[title] += 1
        lastRead = title

    lineCount = 0
    lastReadCount = 0
    if len(lastRead) > 2:
        for (title, sequence) in FastaIterator(inputFasta):
            lineCount += 1

            if title == lastRead:
                lastReadCount = lineCount

    try:
        percentHit = ( (len(dResults) * 100.0) / lastReadCount )
    except ZeroDivisionError:
        percentHit = 0

    try:
        percentComplete = ((lastReadCount*100.0)/lineCount)
    except ZeroDivisionError:
        percentComplete = 0.0
        
    msg = "The last alignment result was found on line %s of %s (%.2f%% complete).  Approximately %.2f%% of the queries are hits." % (
        lastReadCount, lineCount, percentComplete, percentHit)
    if logFile:
        logFile.write(msg)
        logFile.write("\n")
    else:
        print msg

    return percentComplete
