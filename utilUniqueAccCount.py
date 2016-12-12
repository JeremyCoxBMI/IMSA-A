__author__ = 'COX1KB'

import sys

from postprocesscount4acc import parseAccession


def extract_gi( name ):
    splits = name.split("|")
    for k in range(len(splits)-1):
        if splits[k] == "gi":
            return splits[k+1]

    #if name[0:2] == "NW":
    #    return name
    #return ""

    return name

def parseAccession(text):
    acc = "0"
    j = 0
    if text[0] == "g" and text[1] == "i":  ##old school format
        for i in text.split("|"):
            if i == "ref":
                acc = text.split("|")[j + 1].split(".")[0]  # exclude .suffix

                break
            j += 1
    else:  ##new format (I presume... )
        acc = text.split(".")[0]

    return acc



def AddTuple (a, b):
    return (a[0]+b[0], a[1]+b[1])

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "Proper usage is "
        print "  python ", sys.argv[0], " <besthits bln file>"
        print "  it is your option to redirect output to a file"
        exit()

    uniq = {}
    partial = {}
    gis = {}

    for line in open(sys.argv[1], "r"):
        splits = line.split()
        match_pc = float(splits[-1])  #match PerCentage
        #gi = extract_gi(splits[1])
        gi = parseAccession(splits[1])

        if len(gi) > 0:
            gis[gi] = 1
            if match_pc == 1:
                if not gi in uniq:
                    uniq[gi] = 0
                uniq[gi] += 1
            else:
                if not gi in partial:
                    partial[gi] = (0,0)
                partial[gi] = AddTuple(partial[gi], (1,match_pc))
        else:
            print "Error in processing >> ", line

    numUnique = 0
    print "gi\ttotal\tunique count\tpartial count\tpartial sum"
    for gi in gis:
        if gi in uniq:
            uc = uniq[gi]
        else:
            uc = 0
        if gi in partial:
            pc = partial[gi][0]
            ps = partial[gi][1]
        else:
            pc = 0
            ps = 0
        if uc > 0:
            numUnique += 1
        tot = uc + ps
        print gi, "\t", tot, "\t", uc, "\t", pc, "\t", ps

    print >> sys.stderr, sys.argv[1], " contains ", numUnique, " gi with unique alignments"