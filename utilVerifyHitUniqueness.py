import sys
from mylibrary import *

# Program takes output from utilCompareHitQualityHost and looks at overlap with non-organism hits

lib = MyLibrary()

#percent overlap human
outFile_0percent = open(sys.argv[1]+".0percent.bln","w")
outFile_10percent = open(sys.argv[1]+".10percent.bln","w")
outFile_20percent = open(sys.argv[1]+".20percent.bln","w")


def printResults(dict_host, dict_organisms):


        if len(dict_organisms) ==0 and len(dict_host) == 0:
            #mybreakdebug = 1
            return

        # if len(dict_organisms) > 0:
        #     key = dict_organisms.keys()[0]
            # if dict_organisms[key][2][0].split("_")[-1]=="227":
            #     mybreakdebug = 1
        # if len(dict_host) > 0:
        #     key = dict_host.keys()[0]
        #     if dict_host[key][2][0].split("_")[-1]=="227":
        #         mybreakdebug = 1

        #print "printResults size hosts\t", len(dict_host), "\tsize organisms\t", len(dict_organisms)
        removeMe = {}
        for organism in dict_organisms:
            hash = dict_organisms[organism][2][0]
            z = 0.0
            for host in dict_host:
                o = dict_organisms[organism]
                h = dict_host[host]
                #if (lib.overlap(o[0],o[1],h[0], h[1])):
                    #print "o: ", o
                z = lib.calculate_overlap(o[0],o[1],h[0],h[1])
                #hash = lib.toHash(dict_organisms[organism][2])
                if hash in removeMe and removeMe[hash] > z:
                        removeMe[hash] = z
                # if hash in removeMe:
                #     if removeMe[hash] < z:
                #         removeMe[hash] = z
                # else:
                #     removeMe[hash] = z
                #print o, "\teliminated"
                break

        for org in dict_organisms:
            #print "org: ", dict_organisms[org]
            prevLine = "\t".join(dict_organisms[org][2]).replace("\n","")
            # print >> sys.stderr,dict_organisms[org][2][0]
            # print >> sys.stderr,""
            # print >> sys.stderr,""
            length = int(dict_organisms[org][2][0].split("_")[-1])
            #hash = lib.toHash(dict_organisms[org][2])
            if dict_organisms[org][2][0] in removeMe:
                overlap = float(100)*(removeMe[  dict_organisms[org][2][0]  ])/length
            else:
                overlap =0.0

            if overlap > 0:
                print "not\t{0:.2f}\t".format(overlap) + prevLine
            else:
                print "unique\t0.00\t"+prevLine

            if overlap == 0.0:
                outFile_0percent.write("\t".join(dict_organisms[org][2]))
                outFile_10percent.write("\t".join(dict_organisms[org][2]))
                outFile_20percent.write("\t".join(dict_organisms[org][2]))
            elif overlap <= 10.0:
                outFile_10percent.write("\t".join(dict_organisms[org][2]))
                outFile_20percent.write("\t".join(dict_organisms[org][2]))
            elif overlap <= 20.0:
                outFile_20percent.write("\t".join(dict_organisms[org][2]))

            # if overlap > 0.0:   #removeMe[dict_organisms[org][2][0]] > 0:
            #     print "not\t{0:.2f}\t".format(overlap) + prevLine
            # else:
            #     print "unique\t0.00\t"+prevLine
                ###




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

    outSet = [ "unique","human overlap %","query", "reference", "species taxon", "species name","matches", "align length", "percent identity",
               "true pident by query length", "query start", "query end", "note" , "ref start", "ref end", "bitscore"]
    print "\t".join(outSet)

    #firstLine = True
    k=1
    for line in open(compareHitHostFile):

        if k > 2:
            splits = line.split("\t")
            # if splits[0].split("_")[-1]=="227":
            #     mybreakdebug = 1

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

            #encode special name to make each one unique
            if (taxon == 9606 or taxon == 40674):
                dict_host[splits[0]+splits[1]+str(query_start)+str(query_end)] = (query_start, query_end, splits)
            else:
                dict_organisms[splits[0]+splits[1]+str(query_start)+str(query_end)]= (query_start, query_end, splits)
        k += 1
