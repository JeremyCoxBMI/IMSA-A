#!/usr/local/bin/python
#
# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)

import sequenceUtils

def keepAllBestHits(blastInput, blastOutput, ties=True, addHitColumn=True):
    """Filters a blast input file, writing only the best hit for each query input.  If there is a tie for the
    best hit, all the best hits are written if ties is True, otherwise only the first is kept.
    If addHitColumn is true then a column is added with the hit count for the line -- if this hit is the only
    hit for the query the value is 1, otherwise the value is 1/numHits, in other words 0.5 if there are two
    best hits, 0.33 if there are three, etc.
    Reads the whole file into memory so cannot be used for large blast results."""

    queryDict = {}

    for line in open(blastInput):
        pieces = line.split()
        if len(pieces) != 12 and len(pieces) != 13:
            print "Ignoring line, it has %s pieces:" % len(pieces)
            print line
            continue

        query = pieces[0]
        target = pieces[1]
        bitScore = float(pieces[11])
        if not queryDict.has_key(query):
            queryDict[query] = (bitScore, [pieces[:12]])
        else:
            if queryDict[query][0] < bitScore:
                queryDict[query] = (bitScore, [pieces[:12]])
            elif ties and queryDict[query][0] == bitScore:
                # only save if the target is different
                found = False
                for previousPieceList in queryDict[query][1]:
                    if target == previousPieceList[1]:
                        found = True
                        break
                if not found:
                    queryDict[query][1].append(pieces[:12])

    out = open(blastOutput, "w")
    for (query, (bitScore, piecesList)) in queryDict.iteritems():
        eachHitVal = 1.0 / len(piecesList)
        for oneLinePieces in piecesList:
            if not addHitColumn:
                out.write("\t".join(oneLinePieces))
                out.write("\n")
            else:
                out.write("\t".join(oneLinePieces))
                out.write("\t%.3f\n" % eachHitVal)


def returnRecordTitles(blastDictFasta, blastResults):
    """For each hit in the blast results, returns the full title of the record in the blastDictFasta."""

    d = {}
    for (title, seq) in sequenceUtils.FastaIterator(blastDictFasta):
        d[title] = seq

    for line in open(blastResults):
        (read, rec) = line.split()[:2]
        found = False
        for (k,v) in d.iteritems():
            if k.find(rec) >= 0:
                print "%s\t%s\t%s" % (read, rec, k)
                found = True
                break
        if not found:
            print "%s\t%s\tNOT FOUND" % (read, rec)
                
