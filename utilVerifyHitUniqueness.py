import sys
from mylibrary import *

# Program takes output from utilCompareHitQualityHost and looks at overlap with non-organism hits

lib = MyLibrary()



def printResults(dict_host, dict_organisms, dict_0, dict_10, dict_20):


        # if len(dict_organisms) ==0 and len(dict_host) == 0:
        #     #mybreakdebug = 1
        #     return

        removeMe = {}

        #calculate overlaps with host sequences for each Microbiome Sequence
        for organism in dict_organisms:
            o = dict_organisms[organism]
            hash = o[2][0]

            #overlap minimum is 0
            removeMe[hash] = 0.0
            for host in dict_host:
                h = dict_host[host]
                z = lib.calculate_overlap(o[0],o[1],h[0],h[1])
                #maximum overlap
                if removeMe[hash] < z:
                    removeMe[hash] = z
                    # print >> sys.stderr, "REMOVING\t" + "\t".join(o[2]).replace("\n","")
                    # print >> sys.stderr, "query: [",str(o[0]),", ",str(o[1]),"]\thost [",str(h[0]),", ",str(h[1]),"]"

            prevLine = "\t".join(o[2]).replace("\n","")
            length = int(o[2][0].split("_")[-1]) #queryName

            if removeMe[o[2][0]] > 0: #queryName
                overlap = float(100)*(removeMe[  o[2][0]  ])/length
                print "not\t{0:.2f}\t".format(overlap) + prevLine
            else:
                overlap =0.0
                print "not\t{0:.2f}\t".format(overlap) + prevLine

            # if overlap > 0:
            #     print "not\t{0:.2f}\t".format(overlap) + prevLine
            # else:
            #     print "unique\t0.00\t"+prevLine

            if overlap == 0.0:
                dict_0[ (o[2][0], o[2][1]) ] = 1
                # outFile_0percent.write("\t".join(o[2]))
                # outFile_10percent.write("\t".join(o[2]))
                # outFile_20percent.write("\t".join(o[2]))
            elif overlap <= 10.0:
                dict_10[ (o[2][0], o[2][1]) ] = 1
                #outFile_10percent.write("\t".join(o[2]))
                #outFile_20percent.write("\t".join(o[2]))
            elif overlap <= 20.0:
                dict_20[ (o[2][0], o[2][1]) ] = 1
                #outFile_20percent.write("\t".join(o[2]))

if __name__ == "__main__":

    compareHitHostFile = sys.argv[1]
    bln_source = sys.argv[2]
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
    dict_0 = {}
    dict_10 = {}
    dict_20 = {}

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
                printResults(dict_host, dict_organisms,dict_0, dict_10, dict_20)

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

    #percent overlap human
    outFile_0percent = open("0percent.bln","w")
    outFile_10percent = open("10percent.bln","w")
    outFile_20percent = open("20percent.bln","w")


    for line in open(bln_source):
        if len(line) > 10:   #skip empty lines
            splits = line.split()
            coor = (splits[0], splits[1])

            if coor in dict_0:
                outFile_0percent.write(line)
                outFile_10percent.write(line)
                outFile_20percent.write(line)
            elif coor in dict_10:
                outFile_10percent.write(line)
                outFile_20percent.write(line)
            elif coor in dict_20:
                outFile_20percent.write(line)