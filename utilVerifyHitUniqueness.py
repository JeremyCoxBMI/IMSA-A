import sys


# Program takes output from utilCompareHitQualityHost and looks at overlap with non-organism hits

def toHash(aList):
    #print aList
    #print "YO"
    return aList[1]+aList[8]



def calculate_overlap(a1, a2, b1, b2):
    if (a1 < b1 and a2 <= b2 ):
        return a2 - b1
    elif (a1 < b1 and b2 <= a2):
        return b2 - b1
    elif (b1 <= a1 and a2 <= b2):
        return a2 - a1
    elif (b1 <= a1 and b2 <= a2):
        return b2 - a1


def overlap(a1, a2, b1, b2):
    result = False

    if  (a1 <= b2 and b2 <= a2) \
        or  \
        (a1 <= b1 and b1 <= a2 ) :
        result = True

    return result

def printResults(dict_host, dict_organisms):

    #print "printResults size hosts\t", len(dict_host), "\tsize organisms\t", len(dict_organisms)
    removeMe = {}
    for organism in dict_organisms:
        for host in dict_host:
            o = dict_organisms[organism]
            h = dict_host[host]
            if overlap(o[0],o[1],h[0], h[1]):
                #print "o: ", o
                z = calculate_overlap(o[0],o[1],h[0],h[1])
                hash = toHash(o[2])
                if hash in removeMe:
                    if removeMe[hash] < z:
                        removeMe[hash] = z
                else:
                    removeMe[hash] = z
                #print o, "\teliminated"
                break

    for org in dict_organisms:
        #print "org: ", dict_organisms[org]
        if toHash(dict_organisms[org][2]) in removeMe:
            length = int(dict_organisms[org][0].split("_")[-1])
            overlap = float(100)*(removeMe[toHash(dict_organisms[org][2]])/length
            print "not {0:.2}\t".format(overlap)+"\t".join(dict_organisms[org][2]).replace("\n","")
        else:
            print "unique\t"+"\t".join(dict_organisms[org][2]).replace("\n","")




if __name__ == "__main__":

    compareHitHostFile = sys.argv[1]
    # uniqueAlignmentsFile = sys.argv[2]
    #
    # unique = {}
    # for line in open(uniqueAlignmentsFile):
    #   splits = line.split("\t")
    #   name=splits[0]
    #   clade=splits[1]
    #   taxon=int(splits[2])
    #
    #   if clade=="species":
    #     unique[name] = taxon



    currQuery = ""

    dict_organisms = {}
    dict_host = {}

    outSet = [ "unique by overlap","query", "reference", "species taxon", "species name","matches", "align length", "percent identity",
               "true pident by query length", "query start", "query end", "note" , "ref start", "ref end", "bitscore"]
    print "\t".join(outSet)

    #firstLine = True
    k=1
    for line in open(compareHitHostFile):

        if k > 2:
            splits = line.split("\t")
            if currQuery != splits[0]:
                #write answer
                #print "X"
                printResults(dict_host, dict_organisms)

                #reset loop
                currQuery = splits[0]
                dict_host = {}
                dict_organisms = {}

            #add alignment to data structure
            #print "Y"
            #taxon = unique[name]
            taxon = int(splits[2])
            query_start = int(splits[8])
            query_end = int(splits[9])

            if (taxon == 9606 or taxon == 40674):
                dict_host[splits[1]] = (query_start, query_end, splits)
            else:
                dict_organisms[splits[1]] = (query_start, query_end, splits)
        k += 1
