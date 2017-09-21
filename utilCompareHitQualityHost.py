# -*- coding: utf-8 -*-

## SEVERAL INFORMATION PIECES NEEDED
#  to rank sequence quality
#     length of original sequence



# arg 1: file of blast hit outputs
# arg 2: file of uniqueAlignments
# arg 3: file of original sequences
import copy

import sys
from systemSettings import *
from postprocesscount4 import *
from postprocesscount4acc import *


def processReferences2Taxon(refToTaxon):



    print >> sys.stderr, "       : Accessing NCBI for %d taxonomies" % len(taxID_to_lookup)
    #V2: getTaxa for GI alignments and the fungiDB alignments
    names = eutilsWrapper.getTaxa( taxID_to_lookup.keys() )   #lists of integers
    print >> sys.stderr, "       : Got %s taxonomies from NCBI" % (len(names))


if __name__ == "__main__":


    print sys.argv

    names=buildNames()

    seq_len = {}
    unique = {}

    name=""
    bases = 0

    blastInput = sys.argv[1]
    uniqueAlignmentsFile = sys.argv[2]
    fastA_ref_File = sys.argv[3]
    outputFile = sys.argv[4]


    bases = 0
    #aligned sequence length
    for line in open(fastA_ref_File):
        if line[0]==">":
            seq_len[name] = bases

            #reset
            bases = 0
            name = line.split()[0][1:]
        else:
            bases = bases + len(line)-1  #no endl
    #write last
    seq_len[name] = bases


    for line in open(uniqueAlignmentsFile):
      splits = line.split("\t")
      name=splits[0]
      clade=splits[1]
      taxon=int(splits[2])

      if clade=="species":
        unique[name] = taxon

    #print seq_len.keys()
    #exit()


    query_to_results = {}
    refToTaxon = {}

    for line in open(blastInput):
        splits = line.split("\t")
        query = splits[0]
        ref = splits[1]
        pident = splits[2]
        num_mismatch=int(splits[4])
	#print query
        length = seq_len[query]
        query_start = int(splits[6])
        query_end = int(splits[7])
        aligned_len = query_end - query_start + 1

        ref_start = int(splits[8])
        ref_end = int(splits[9])
        bitscore=float(splits[11])

        true_pident = float(aligned_len - num_mismatch) / length

        if not query in query_to_results:
            query_to_results[query] = []   #empty list
        query_to_results[query].append(  (query, ref, aligned_len, length, pident, true_pident, query_start, query_end,ref_start,ref_end, bitscore)  )

        #add reference to be analyzed
        refToTaxon[ref]=1

    #process refToTaxon so it conatins firstTaxa hits

    print >> sys.stderr, 'STEP 110: build fungiLookup'
    fungiLookup = {}  # lookup table from file
    print >> sys.stderr, '        : using lookup file ', fungiLookupFile
    for line in open(fungiLookupFile, 'r'):
        splits = line.split('\t')
        key = splits[0]
        speciesID = int(splits[1])
        fungiLookup[key] = speciesID

    # STEP 120:
    print >> sys.stderr, 'STEP 120: extract accession number targets, SID_to_lookup from file'

    #contains integer taxID to lookup in NCBI -- will contain FungiDB species and regular database taxID
    taxID_to_lookup = {}

    accs = {}
    extractLookupInfo(accs, taxID_to_lookup, fungiLookup, blastInput)

    # STEP 140: build taxonomy dictionary from gis -- parse local file
    print >> sys.stderr, 'STEP 140: build taxonomy dictionary from accs -- parse local file     LONG STEP'
    print >> sys.stderr, '       : checking for pickle'

    import cPickle as pickle
    if not os.path.exists(blastInput+".taxonomy.pickle"):
        print >> sys.stderr, '       : no pickle to load'
        taxonomy = {}
        print >> sys.stderr, "       : Getting taxonomy ids from local file for %s accs" % (len(accs))

        #short circuit if no need to parse file
        if len(accs)>0:
            print >> sys.stderr, '       :       : searching %s' % ACC_BLAST_TAX_DB
            for line in open(ACC_BLAST_TAX_DB):
                (acc, accv, taxid, gi) = line.split()
                if accs.has_key(acc):
                    taxonomy[acc] = int(taxid)
            print >> sys.stderr, '       :       : searching %s' % ACC_EXTRA_BLAST_TAX_DB
            for line in open(ACC_EXTRA_BLAST_TAX_DB):
                (acc, accv, taxid, gi) = line.split('\t')
                if accs.has_key(acc):
                    taxonomy[acc] = int(taxid)

        print >> sys.stderr, '       : creating pickle on disk'
        pickle.dump( taxonomy, open( blastInput+".taxonomy.pickle", "wb" ) )
    else:
        print >> sys.stderr, '       : loading pickle'
        taxonomy = pickle.load(  open( blastInput+".taxonomy.pickle", "rb" )  )

    print >> sys.stderr, "       : Done getting taxonomy ids, got %s ids" % (len(taxonomy))


     # STEP 150: build taxID_to_lookup from taxonomy
    # could be combined in step 140, but 1 pickle is enough
    print >> sys.stderr, 'STEP 150: increase taxID_to_lookup from taxonomy list'
    for k in taxonomy:
        taxID_to_lookup[ taxonomy[k] ] = 1  # writing integers as keys

    #gis contains gi keys as strings
    #taxonomy contains taxid (int) coded to gi keys (string)

    # STEP 160: get the scientific name for each taxa from DB
    print >> sys.stderr, 'STEP 160: get the scientific name for each taxa from DB'

    print >> sys.stderr, "       : Accessing NCBI for %d taxonomies" % len(taxID_to_lookup)
    #V2: getTaxa for GI alignments and the fungiDB alignments
    names = eutilsWrapper.getTaxa( taxID_to_lookup.keys() )   #lists of integers
    print >> sys.stderr, "       : Got %s taxonomies from NCBI" % (len(names))


    #seq_len = {}
    #unique = {}





    # STEP 170: build allHits from blastInput: mapping query, to list of species ID that are hits

    print >> sys.stderr, 'STEP 170: build allHits from blastInput: mapping query, to list of species ID that are hits'

    outF = open(outputFile, "w")
    currQuery = ""
    hitList = []
    first = True
    taxaCount = {}
    speciesCount = {}
    genusCount = {}
    familyCount = {}


    #temporary list to store values looked up via eUtils
    lookups = {}

    outSet = [ "query", "reference", "species taxon", "species name","matches", "align length", "percent identity",
               "true pident by query length", "query start", "query end", "note" , "ref start", "ref end", "bit score"]
    print "\t".join(outSet)

    for line in open(blastInput, "r"):
        #loop control
        splits = line.split('\t')

        #what about excluding human?
        #print unique[splits[0]]
        #print splits[1]


        speciesName = "not found"

        #Need a filter to not lookup human sequences
        # if splits[0] in unique and unique[splits[0]] != 9606 and unique[splits[0]] != 40674:  #filters human?  need parameteR?
        #     lookup = False
        # human = False
        # #Identify species
        # if human:
        #     speciesID = 40674
        #     speciesName = "Mammal"
        # else:
        #analysis of line - do every time, even if find end
        k = splits[1]
        acc = parseAccession(k, fungiLookup)

        fullTax = 0
        if acc == "0" and k in fungiLookup:  #gi not found, cuz new fungiDB reference sequence name
            speciesID = fungiLookup[k]
            # if not speciesID == 40674:
            #     (taxaID, speciesID, genusID, familyID) = lookupFungi( k, fungiLookup, names)
        elif acc in taxonomy:  #2015-07-24
            taxID = taxonomy[acc]
            #taxaID = taxID
            if taxID in names:
                fullTax = names[taxID]
                taxaID = taxID
                (speciesID, genusID, familyID) = fullTax2IDS( fullTax )
            elif taxID in lookups.keys():
                #new condition: if I looked up and got a merged entry, I saved it and this looks it up
                (taxID, fullTax) = lookups[taxID]
                (speciesID, genusID, familyID) = fullTax2IDS( fullTax )
                taxaID = taxID
            else:
                # attempt to lookup missing values
                print >> sys.stderr, "taxID %d was not found in the NCBI lookup; it is being looked up again, possible merged record" % taxID
                #print "taxID %d was not found in the NCBI lookup; it is being looked up again, possible merged record" % taxID
                single =  eutilsWrapper.getTaxa( [taxID] )
                if len(single)==1:
                    newTaxaID = single.keys()[0]
                    lookups[taxID] = (newTaxaID, single[newTaxaID])
                    print >> sys.stderr, "taxID %d was found in the NCBI lookup as taxID %d" % (taxID, int(newTaxaID))
                    #print "taxID %d was found in the NCBI lookup as taxID %d" % (taxID, int(newTaxaID))
                    fullTax = single[newTaxaID]
                    (speciesID, genusID, familyID) = fullTax2IDS( fullTax )
                else:
                    print >> sys.stderr, "taxID %d was STILL not found in the NCBI lookup; it is being ignored" % taxID
                    #missing.write("taxaID\t"+str(taxID)+"\n")
                    taxaID = -1
                    speciesID = -1
                    genusID = -1
                    familyID = -1
                    #fullTax = 0   #using default
        else:
            print "not in taxonomy DB\t"+line
            taxaID = -1
            speciesID = -1
            genusID = -1
            familyID = -1
            #fullTax = 0    #using default



        if not fullTax == 0:
            speciesName = fullTax.getSpecies().sname
        else:
            # if speciesID == 40674:
            #     speciesName = "mammal"
            if speciesID in names:
                speciesName = names[speciesID]
            else:
                speciesName = "not found"

        queryList = query_to_results[ splits[0] ]
        mytuple = ()
        for x in range(len(queryList)):
            if splits[1] == queryList[x][1]: #reference name matches
                mytuple = copy.deepcopy(queryList[x])
                #queryList.remove(x)  #guarantees same reference sequence scores not reported for multiple hits
                del queryList[x]
                break
        #(query, ref, aligned_len, length, pident, true_pident) = mytuple
        (query, ref, aligned_len, length, pident, true_pident, query_start, query_end,ref_start,ref_end, bitscore) = mytuple
        true_pident = int(100*true_pident)

        if splits[0] in unique and unique[splits[0]] == speciesID:  #winner
            outSet = [ splits[0], splits[1], str(speciesID), str(speciesName), str(aligned_len), str(length), str(pident),
                       str(true_pident), str(query_start), str(query_end), "<- query||ref -->", str(ref_start),
                       str(ref_end), str(bitscore)+"\t*"]
            print "\t".join(outSet)
        else:
            outSet = [ splits[0], splits[1], str(speciesID), str(speciesName), str(aligned_len), str(length), str(pident),
                       str(true_pident), str(query_start), str(query_end), "<- query||ref -->", str(ref_start),
                       str(ref_end), str(bitscore)]
            print "\t".join(outSet)
        #else:
            #print all sequences?
        # else:
        #     queryList = query_to_results[ splits[0] ]
        #     mytuple = ()
        #     for x in range(len(queryList)):
        #         if splits[1] == queryList[x][1]: #reference name matches
        #             mytuple = copy.deepcopy(queryList[x])
        #             #queryList.remove(x)  #guarantees same reference sequence scores not reported for multiple hits
        #             del queryList[x]
        #             break
        #     (query, ref, aligned_len, length, pident, true_pident, query_start, query_end,ref_start,ref_end, bitscore) = mytuple
        #     true_pident = int(100*true_pident)
        #     outSet = [ splits[0], splits[1], str(9606), "human", str(aligned_len), str(length), str(pident),
        #                str(true_pident), str(query_start), str(query_end), "<- query||ref -->", str(ref_start),
        #                str(ref_end), str(bitscore)]
        #     print "\t".join(outSet)
            # if splits[0] in unique:  #by if statement above, mammal or human
            #     outSet = [ splits[0], splits[1], str(9606), "human", str(aligned_len), str(length), str(pident),
            #            str(true_pident), str(query_start), str(query_end), "<- query||ref -->", str(ref_start),
            #            str(ref_end), str(bitscore)]
            #     print "\t".join(outSet)
            # else:
            #
            #     outSet = [ splits[0], splits[1], str(speciesID), speciesName, str(aligned_len), str(length), str(pident),
            #            str(true_pident), str(query_start), str(query_end), "<- query||ref -->", str(ref_start),
            #            str(ref_end), str(bitscore)+"\t*"]
            #     print "\t".join(outSet)

    #end for

