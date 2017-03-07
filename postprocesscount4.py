#!/usr/local/bin/python
'''
Copyright (C) 2016  Jeremy Wayne Cox    jeremy.cox@cchmc.org

This code is distributed under Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0) License.

For a human readable interpretation of the license, please see
http://creativecommons.org/licenses/by-nc-sa/4.0/

NOTICE
This code is a modification and extension of the program Arron-IMSA
distributed by sourceforge.net.
https://sourceforge.net/projects/arron-imsa/

The IMSA license file and notice of copyright is included with this release.
'''


"""
JWC (Jeremy.Cox@cchmc.org) modifying this code to do different postprocessing.

THIS SCRIPT CONTAINS POSTPROCESSOR UPGRADES ASSOCIATED WITH IMSA+A
However, the postprocessor can be used with regular IMSA pipeline results.

MODIFICATIONS

* Counting Scheme

Original IMSA tabulated 'total count' = hits+partial sums,
IMSA+A tracks total count (float), but also counts
unique hits(int), partial hits(int), and partial hit sum for the partial hits (float).

Essentially, this modification required substituting a tuple for a float throughout the code.
4count results are filtered to include only results with unique hits > 0.
This required tweaking multiple libraries, which are now embedded here.
Functions, classes, and variables which were tweaked were renamed with suffix JWC.

Two files are created
    instead of
        FILE.report.txt
    two files
        FILE.report.IMSA_count.txt
        FILE.report.IMSA+A_4count.txt   #4 counts, filtered to include unique count > 0

* Virus accomodations

To help with accurate classification of viruses, a new first taxon report is created.
(FILE.firstTaxon.IMSA_count.txt and FILE.firstTaxon.IMSA+A_4count.txt).

* divisionID removed

divisionIDs were removed from analysis.  This only reports firstTaxon, species, genus, family.

_________________________
 Additional New Features
*************************

- Can recognize species by reference sequence names not containing a gi number, when those numbers are added to file
referred to in systemSettings.fungiLookUp
- Counting method modified to recalculate uniqueness at every clade.  This is very similar to LCA binning by MEGAN.
- Much more detailed messages during execution, directed to STDERR.
- Optimized to read blast besthits input line by line, instead of reading entire file all at once to allow for larger files.
- Saves computationally intensive data objects to pickles to support resuming / recalculating.

If the blast alignments file changes, all these pickles must be deleted for accurate computation.

Generates a new report:   FILE.tax_unidentifiedTaxaAlignments.txt
These missing entries can be looked up and added to systemSettings.EXTRA_BLAST_TAX_DB, which IMSA will parse to find names.
    Note that if systemSettings.EXTRA_BLAST_TAX_DB file is updated, you must delete .taxonomy.pickle when re-running.

Handles unknowns gracefully, reporting taxon ID = -1
    As previously mentioned, there is built-in support to help eliminate this with USER intervention.
    However, IMSA-A now has a way to handle the error gracefully and report gracefully, so that the program does not crash.


"""

# This script is designed to do the typical processing done after the IMSA pipeline completes.
# This script takes in the blast vs. NT, filters to best hits, and runs the taxonomy report.

# Paths to local data files

import config
#systemSettings overrides config settings if necessary
from systemSettings import *
import sys
from getopt import getopt
import os
import pipelineUtils
import blastUtils
import eutilsWrapper

HELP_STRING = """Given the results of the filtered reads blasted against reference database, this script
does the best hit blast filter and then runs the taxonomy reports.  Uses Python 2.7.

Required Inputs:
    -b    Blast output file

Optional Inputs:
    -f    this option is deprecated for IMSA+A 4 count.  (It may still function)
          This program now outputs file with extension .firstTaxon.IMSA+A_4count.txt for better virus classification.

REQUIRES PYTHON 2.7
"""

#new function
def fullTax2IDS( fullTax ):
    """ Converts fullTax to taxon ID's, considering error checking.
        (Original IMSA did not have error checking)
    :param fullTax: INPUT: eutilsWrapper.Taxonomy object
    :return:        (int firstTaxonID, int speciesID, int genusID, int familyID)
    """
    species = fullTax.getSpecies()
    if species:
        speciesID = species.taxId
    else:
        speciesID = -1

    #V2 new code -- this change was made in v2, now in v3 it is incorporated to a function
    genus = fullTax.getGenus()
    if genus:
        genusID = genus.taxId
    else:
        genusID = -1

    family = fullTax.getFamily()
    if family:
        familyID = family.taxId
    else:
        familyID = -1

    return (speciesID, genusID, familyID)

# helper function to use local augmented reference sequence names for lookup
def lookupFungi( k, fungiLookup, names):
    """ takes a reference sequence name and lookups relevant taxonomy
        -1 indicates not found
    :param k:           (string) INPUT: name of ref sequence to lookup
    :param fungiLookup: (dict) INPUT: mapping ref sequence names (non-gi) to speciesID
    :param names:       (dict) INPUT: names database from etuils lookup: mapping taxID (int) to
                        eutilsWrapper.Taxonomy object
    :return:            (int taxaID, int speciesID, int genusID, int familyID)
    """
    speciesID = fungiLookup[k]
    taxaID = speciesID
    fullTax = names[speciesID]
    genus = fullTax.getGenus()
    if genus:
        genusID = genus.taxId
    else:
        genusID = -1
    family = fullTax.getFamily()
    if family:
        familyID = family.taxId
    else:
        familyID = -1

    return (taxaID, speciesID, genusID, familyID)

# new function
def getTaxa( sID, names, nodes ):
    """ Finds genus ID and family ID given the species ID from the NCBI taxonomy tree
        -1 represents not found, which happens a lot with viruses
    :param sID:     integer INPUT: species ID
    :param names:   NCBI names database mapping: int taxaID to taxaName (str)
    :param nodes:   NCBI nodes database mapping:int taxaID to (parentID, level)
    :return:        (int, int) OUTPUT: (genus ID, family ID)
    """
    genusID = -1
    familyID = -1

    cur = sID
    while cur != 1: #while not root
        if cur in names:
            if nodes[cur][1] == "genus":
                genusID = cur
            if nodes[cur][1] == "family":
                familyID = cur
                break
            cur = nodes[cur][0]
        else:
            #genusID and familyID can no longer be found
            break

    return (genusID, familyID)

# new function to add immutable tuples.  Note that tuple[0] += tuple2[0] does not work
def addToTuple(first, second):
    return ( first[0] + second[0], first[1] + second[1], first[2] + second[2], first[3] + second[3] )

# original IMSA code modified to handle 4-tuple
class taxNodeJWC:
    def __init__(self, sname, id, rank):
        self.sname = sname
        self.id = id
        self.rank = rank
        self.count = (0, 0, 0, 0)
        self.totalCount = (0, 0, 0, 0)

    def __repr__(self):

        width = 3.0
        if self.totalCount[0] > 0:
            ratio = float(self.count[0]) / self.totalCount[0]
            if ratio > config.MAX_NODE_SIZE_CUTOFF:
                width = config.MAX_NODE_SIZE
            else:
                width = (ratio * (1 / config.MAX_NODE_SIZE_CUTOFF) * (
                    config.MAX_NODE_SIZE - config.MIN_NODE_SIZE)) + config.MIN_NODE_SIZE

                #print "%s\t\t%.3f\t\t%.3f" % (self.totalCount, ratio, width)

        height = width / config.WIDTH_TO_HEIGHT
        fontcolor = "black"
        return '%s [color=%s,fixedsize=true,width=%.1f,height=%.1f,shape=ellipse,style=bold,fontcolor=%s,label="%s\\n%.0f"];' % (
            self.id, config.RANK_TO_COLOR[self.rank], width, height, fontcolor, self.sname, self.count[0])


# original IMSA code modified to handle 4-tuple
def getTaxDictionaryStrJWC(d):
    retval = ""
    sortedTax = sorted(d.iteritems(), reverse=True)
    #print sortedTax

    for ((taxId, sname), count) in sortedTax:
        retval += "%s\t%s\t%s\n" % (taxId, sname, catV(count))

    return retval

# Helper function
# converts 4-tuple to a string
def catV(v):
    #v is a 4-tuple
    result = str(v[0]) + '\t' + str(v[1]) + '\t' + str(v[2]) + '\t' + str(v[3]) + '\t'
    return result

# global variable
blastResults = ""

# unchanged IMSA code
def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt(argv[1:], "hb:f:")
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
    fastaInput = None
    for (opt, opt_arg) in optlist:
        if opt == "-h":
            print ""
            print HELP_STRING
            sys.exit(1)

        elif opt == "-b":
            blastResults = opt_arg

        elif opt == "-f":
            fastaInput = opt_arg

    if not blastResults:
        print "You must specify a blast result file (use -b)."
        print
        print HELP_STRING
        sys.exit(1)

    doPostProcess(blastResults, fastaInput)

# unchanged IMSA code
def doPostProcess(blastResults, fastaInput):
    rootName, extension = os.path.splitext(blastResults)
    bestHitsName = rootName + ".bestHits.bln"
    bestHitsName2 = rootName + ".besthits.bln"
    dp = rootName + ".bubble_"
    reportName = rootName + ".tax_"

    # do best hits
    print >> sys.stderr, 'STEP 5: calculate besthits file'
    print >> sys.stderr, '      : besthits filename: >', bestHitsName, "<"
    print >> sys.stderr, '      : besthits filename: >', bestHitsName2, "<"

    #updated 2015-08-04 so to stop stupid bug from not seeing ".besthits.bln" files, like when copied from DOS
    if os.path.exists(bestHitsName2):
        bestHitsName = bestHitsName2
        print >> sys.stderr, '      : found besthits filename: >', bestHitsName2, "<"

    if not os.path.exists(bestHitsName):
        blastUtils.keepAllBestHits(blastResults, bestHitsName)
    else:
        print >> sys.stderr, '<skip>: besthits file already exists'

    # report taxonomy
    reportTaxonomyJWC(bestHitsName, printTargets=True, outputPrefix=reportName, dotPrefix=dp)

    if fastaInput:
        pipelineUtils.getFastaForTaxonomy([10239], bestHitsName, fastaInput, rootName + "_viral.fa")
        pipelineUtils.getFastaForTaxonomy([10292], bestHitsName, fastaInput, rootName + "_herp.fa")
        pipelineUtils.getFastaForTaxonomy([11632], bestHitsName, fastaInput, rootName + "_retro.fa")

# new function
def extractLookupInfo(gis, sID_to_lookup, fungiLookup, blastInput):
    """ This function pulls gi numbers out of blast input file.
        Gi numbers and speciesID for non-gi reference sequences are returned.
    :param gis:             (dict) OUTPUT: list of gi numbers (as text) found in blastinput file
    :param sID_to_lookup:   (dict) OUTPUT: list of speciesID we need to lookup in NCBI
    :param fungiLookup:     (dict) INPUT: mapping ref sequence names (non-gi) to speciesID
    :param blastInput:      (string) INPUT: name of file containing blast hits
    :return:
    """
    input = open(blastInput,"r")

    for line in input:
        #parse blast lines for reference sequence names
        splits = line.split('\t')

        # search for gi number in seq name
        gi = 0 #checking for gi number
        j=0
        for i in splits[1].split("|"):
            if i == "gi":
                gi = splits[1].split("|")[j + 1]
                break
            j += 1
        if gi == 0:  #gi not found, try lookup in additional seqname/speciesID list (fungiLookup)
            if splits[1] in fungiLookup:
                sID = fungiLookup[ splits[1] ]
                sID_to_lookup[sID] = 1
        else:
            gis[gi] = 1

#new function
def findingUnknownAlignments( taxonomy, blastInput, filenameOUT):
    """ If a blast file and looks for any gi codes, which were not found, then writes to a file.
        Can be used by USER to update the systemSettings.EXTRA_BLAST_TAX_DB file.
    :param taxonomy:    (dict) INPUT: maps gi numbers found in database to their taxon ID
                                        used here as a list of gi numbers found
    :param blastInput:  (string) INPUT: name of input blast file
    :param filenameOUT: (string) INPUT: name of file to write entries to
    :return:
    """

    outF = open(filenameOUT,"w")
    for line in open(blastInput,"r"):
        #parse blast lines for reference sequence gi numbers
        splits = line.split('\t')

        gi = 0 #checking for gi number
        j=0
        for i in splits[1].split("|"):
            if i == "gi":
                gi = splits[1].split("|")[j + 1]
                break
            j += 1

        #report on gi number not found
        if gi != 0 and not gi in taxonomy:
            #gi can not be non-gi fungiDB code
            #if not found earlier in taxadump files (thus not in taxonomy)
            outF.write(line)
    outF.close()


# new function
def processHitList(taxaCount, speciesCount, genusCount, familyCount, hitList, currQuery, uniqueLog):
    """ Counting at each clade varies by the uniqueness criterion.  Therefore, this function calculates hits
        at multiple clades given the hits in hitlist.
    :param taxaCount:   INPUT/OUTPUT: dict mapping taxonID to 4tuple count.  This function adds to it.
    :param speciesCount: INPUT/OUTPUT: dict mapping taxonID to 4tuple count.  This function adds to it.
    :param genusCount:  INPUT/OUTPUT: dict mapping taxonID to 4tuple count.  This function adds to it.
    :param familyCount: INPUT/OUTPUT: dict mapping taxonID to 4tuple count.  This function adds to it.
    :param hitList:     INPUT: list of hits to taxons as tuples (taxaID, speciesID, genusID, familyID)
    :param currQuery:   INPUT: (string) containing the query name to which this hitList corresponds
    :param uniqueLog:   INPUT: open file for writing log messages
    :return:
    """

    if len(hitList) == 1: #unique hit
        result = (1,1,0,0)

        taxaID = hitList[0][0]
        if not taxaCount.has_key(taxaID):
            taxaCount[taxaID] = (0,0,0,0)
        taxaCount[taxaID] = addToTuple( taxaCount[taxaID], result )
        uniqueLog.write(currQuery+"\t"+"firstTaxa\t"+str(taxaID)+"\n")

        speciesID = hitList[0][1]
        if not speciesCount.has_key(speciesID):
            speciesCount[speciesID] = (0,0,0,0)
        speciesCount[speciesID] = addToTuple( speciesCount[speciesID], result )
        uniqueLog.write(currQuery+"\t"+"species\t"+str(speciesID)+"\n")

        genusID = hitList[0][2]
        if not genusCount.has_key(genusID):
            genusCount[genusID] = (0,0,0,0)
        genusCount[genusID] = addToTuple( genusCount[genusID], result )
        uniqueLog.write(currQuery+"\t"+"genus\t"+str(genusID)+"\n")

        familyID = hitList[0][3]
        if not familyCount.has_key(familyID):
            familyCount[familyID] = (0,0,0,0)
        familyCount[familyID] = addToTuple( familyCount[familyID], result )
        uniqueLog.write(currQuery+"\t"+"family\t"+str(familyID)+"\n")
    else:   # not a unique hit

        taxaD = {}
        speciesD = {}
        genusD = {}
        familyD = {}

        #count entries at every clade
        for entry in hitList:
            if entry[0] in taxaD:
                taxaD[entry[0]] += 1
            else:
                taxaD[entry[0]] = 1

            if entry[1] in speciesD:
                speciesD[entry[1]] += 1
            else:
                speciesD[entry[1]] = 1

            if entry[2] in genusD:
                genusD[entry[2]] += 1
            else:
                genusD[entry[2]] = 1

            if entry[3] in familyD:
                familyD[entry[3]] += 1
            else:
                familyD[entry[3]] = 1

        #determine if unique at each clade or not, count accordingly
        for taxaID in taxaD:
            partial = float(taxaD[taxaID]) / len(hitList)
            if partial == 1: #unique hit
                result = (1, 1, 0, 0)
                uniqueLog.write(currQuery+"\tfirstTaxa\t"+str(taxaID)+"\n")
            else:
                result = (partial, 0, 1, partial)
            # write result to the count object
            if not taxaCount.has_key(taxaID):
                taxaCount[taxaID] = (0,0,0,0)
            taxaCount[taxaID] = addToTuple( taxaCount[taxaID], result )
            #if taxaID == -1:
            #    print "This is a debug trap for debugging."

        for speciesID in speciesD:
            partial = float(speciesD[speciesID]) / len(hitList)
            if partial == 1: #unique hit
                result = (1, 1, 0, 0)
                uniqueLog.write(currQuery+"\tspecies\t"+str(speciesID)+"\n")
            else:
                result = (partial, 0, 1, partial)
            if not speciesCount.has_key(speciesID):
                speciesCount[speciesID] = (0,0,0,0)
            speciesCount[speciesID] = addToTuple( speciesCount[speciesID], result )

        for genusID in genusD:
            partial = float(genusD[genusID]) / len(hitList)
            if partial == 1: #unique hit
                result = (1, 1, 0, 0)
                uniqueLog.write(currQuery+"\t"+"genus\t"+str(genusID)+"\n")
            else:
                result = (partial, 0, 1, partial)
            if not genusCount.has_key(genusID):
                genusCount[genusID] = (0,0,0,0)
            genusCount[genusID] = addToTuple( genusCount[genusID], result )

        for familyID in familyD:
            partial = float(familyD[familyID]) / len(hitList)
            if partial == 1: #unique hit
                result = (1, 1, 0, 0)
                uniqueLog.write(currQuery+"\t"+"family\t"+str(familyID)+"\n")
            else:
                result = (partial, 0, 1, partial)
            if not familyCount.has_key(familyID):
                familyCount[familyID] = (0,0,0,0)
            familyCount[familyID] = addToTuple( familyCount[familyID], result )
    # end   if len(hitList) == 1:
    return


# this is a helper function used to determine Kingdom of groups
# without hard coding the Kingdom names, this is difficult; kingdom or superkingdom?
# let's try good old simple Kingdom for the time being

def findKingdom( taxaID, nodes):
    curr = int(taxaID)
    if curr in nodes:
        while nodes[curr][1] != "kingdom" and nodes[curr][1] != "superkingdom" :
            curr = nodes[curr][0]
    else:
        curr = 1    #default determination is root

    return curr

    # if (curr == -1): #ERROR code
    #     return -1
    # while curr != 1:
    #         #caps required
    #     if names[curr] == "Viruses" or names[curr] == "Bacteria" or names[curr] == "Fungi":
    #         break
    #     if names[curr] == "Stramenopiles":
    #         print >> sys.stderr, "Stramenopiles found %d %s, counting as FUNGI " % (int(taxaID),names[int(taxaID)])
    #         curr = 4751 #change to Fungi
    #         break
    #     curr = nodes[curr][0]
    #
    # if curr == 1:
    #     return -1
    # return names[curr]



# original IMSA code modified to handle 4-tuple and new counting scheme and upgrades
def fillOutCountFiles(KEYS, familyCount, taxaNames, taxaNodes, outFile, outFile2, category, savedKingdoms):
    outFile.write("%s\t%s\t%s\t%s\t%s\n" %
                  (category+" ID", "Scientific Name", "IMSA+A count", "Kingdom", "Total\tUnique clade hits\tPartial clade hits\tPartial clade sum"))

    outFile2.write("%s\t%s\t%s\n" %
               (category+" ID", "Scientific Name", "IMSA count"))

    for taxID in KEYS:
        if taxID == -1 or taxID == "-1":
            sname = -1
            # count = taxaCount[taxID]
        elif taxID in taxaNames:
            sname = taxaNames[taxID]
        else:
            sname = -1
        count = familyCount[taxID]

        if taxID in taxaNodes:
            kingdom = findKingdom(taxID, taxaNodes)
        elif taxID in savedKingdoms:
            kingdom = savedKingdoms[taxID]
        else:
            kingdom = 1


        if (count[1] > 0):  # if unique count > 0
            myLine = "%s\t%s\t%s\t%s\t%s\n" % (str(taxID), sname, count[1], taxaNames[kingdom], catV(count))
            outFile.write(myLine)
        myLine = "%s\t%s\t%s\n" % (str(taxID), sname, str(count[0]))
        outFile2.write(myLine)

def prettyEntry(entry):    # 4 tuple: score, taxon ID, scientific name, kingdom
    return '{:15}\t{:35}\t{:15}\t{:15}\n'.format(entry[1],entry[2],str(entry[0]),entry[3])


def processOutputFileHumanReadable(outBig, filename):
    firstLine = True
    listEntries = []
    header = []
    for line in open(filename):
        splits = line.split("\t")
        if firstLine:
            header = splits[0:4]
            firstLine = False
        else:
            if int(splits[2]) > 0:
                listEntries.append( (int(splits[2]), splits[0], splits[1][:33], splits[3]) )

    listEntries.sort(reverse=True)
    outBig.write('{:15}\t{:35}\t{:15}\t{:15}\n'.format(header[0], header[1],header[2], header[3])+"\n")

    for entry in listEntries:
        if entry[2] == 1921421:
            print "debugging trap"
        if entry[1] == -1 or entry[1] == "-1":
            entry = (entry[0], "none", entry[2], entry[3])
        if entry[2] == -1 or entry[2] == "-1":
            entry = (entry[0], entry[1], "not found in DB", entry[3])
        outBig.write(prettyEntry(entry))

    return len(listEntries)

def reportTaxonomyJWC(blastInput, printTargets=False, outputPrefix=None, dotPrefix=None, dotLimit=1.0):
    """Given a blast input file, reports the taxonomy lineage for the hits.  By default the report
    will be printed to standard output, but if an outputPrefix is given then output files, one each
    for species, genus, family and division, will be created.
    The Blast input file should have already been run through a Best Blast filter if you want to only
    count unique hits.  Can have the extra counts column added, as per blastUtils.keepAllBestHits."""

    '''
    fungiLookup - dictionary mapping
        seqID (string) - the reference sequence name used for looking up blast hits
        to
        speciesID (string) - taxon ID number as text
    gis - dictionary mapping                (used as non-repeating list)
        gid (string) - the gid, the reference number used for looking up blast hits (as text)
        to 1
    taxID_to_lookup - dictionary mapping    (used as non-repeating list)
        taxID (int) - taxon ID to lookup
        to 1
    taxonomy - dictionary mapping
        gid (string)
        to
        taxID (int)
    names
        taxID (int)
        to
        eutilsWrapper.Taxonomy object
    4tuple
        (float, int, int, float) representing (IMSA count, unique count, partial count, partial sum)

    taxaCount = {}
    speciesCount = {}
    genusCount = {}
    familyCount = {}
        map
        taxID (int)
        to
        4tuple

    lookups -- dictionary to store eUtil lookup results
        taxID (int)
        to
          _variable_               _class_
        (newTaxaID,             (int)
        single[newTaxaID])      eutilsWrapper.Taxonomy
    '''

    # STEP 110: build fungiLookup
    print >> sys.stderr, 'STEP 110: build fungiLookup'
    fungiLookup = {}  # lookup table from file
    print >> sys.stderr, '        : using lookup file ', fungiLookupFile
    for line in open(fungiLookupFile, 'r'):
        splits = line.split('\t')
        key = splits[0]
        speciesID = int(splits[1])
        fungiLookup[key] = speciesID

    # STEP 120:
    print >> sys.stderr, 'STEP 120: extract gi targets, SID_to_lookup from file'

    #contains integer taxID to lookup in NCBI -- will contain FungiDB species and regular database taxID
    taxID_to_lookup = {}

    # dict containing gi's to lookup, requiring NCBI lookup
    gis = {}

    extractLookupInfo(gis, taxID_to_lookup, fungiLookup, blastInput)
        # gis is indexed by string numbers
        # taxID_to_lookup is indexed by integer numbers

    # STEP 130: write gis table to file
    print >> sys.stderr, 'STEP 130: write gis table to file'
    outGis = open(outputPrefix + "GIS.txt", "w")
    for key in gis:
        outGis.write("%s\n" % (key))
    outGis.close()

    # STEP 140: build taxonomy dictionary from gis -- parse local file
    print >> sys.stderr, 'STEP 140: build taxonomy dictionary from gis -- parse local file     LONG STEP'
    print >> sys.stderr, '       : checking for pickle'

    import cPickle as pickle
    if not os.path.exists(blastInput+".taxonomy.pickle"):
        print >> sys.stderr, '       : no pickle to load'
        taxonomy = {}
        print >> sys.stderr, "       : Getting taxonomy ids from local file for %s gis" % (len(gis))

        #short circuit if no need to parse file
        if len(gis)>0:
            print >> sys.stderr, '       :       : searching %s' % BLAST_TAX_DB
            for line in open(BLAST_TAX_DB):
                (gi, taxid) = line.split()
                if gis.has_key(gi):
                    taxonomy[gi] = int(taxid)
            print >> sys.stderr, '       :       : searching %s' % EXTRA_BLAST_TAX_DB
            for line in open(EXTRA_BLAST_TAX_DB):
                (gi, taxid) = line.split()
                if gis.has_key(gi):
                    taxonomy[gi] = int(taxid)

        print >> sys.stderr, '       : creating pickle on disk'
        pickle.dump( taxonomy, open( blastInput+".taxonomy.pickle", "wb" ) )
    else:
        print >> sys.stderr, '       : loading pickle'
        taxonomy = pickle.load(  open( blastInput+".taxonomy.pickle", "rb" )  )

    print >> sys.stderr, "       : Done getting taxonomy ids, got %s ids" % (len(taxonomy))


    # STEP 145: write Alignments to gi nubmers, which were not found
    print >> sys.stderr, 'STEP 145: write Alignments to gi numbers, which were not found'
    findingUnknownAlignments( taxonomy, blastInput, outputPrefix + "unidentifiedTaxaAlignments.txt")


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

    # STEP 165: write names to file
    print >> sys.stderr, 'STEP 165: write names to file'
    OPP = open(outputPrefix + "NAMES.txt", "w")
    for name in names:
        OPP.write("\n" + str(name) + "\n" + str(names[name]))
    OPP.close()


    # STEP 168: output taxa reports (LCA binning version)
    print >> sys.stderr, 'STEP 168: import ncbi gi database'
    print >> sys.stderr, '       : build taxaNames'
    taxaNames = buildNames()
    print >> sys.stderr, '       : build taxaNodes'
    (taxaNodes, levels) = buildNodes()


    # STEP 170: build allHits from blastInput: mapping query, to list of species ID that are hits

    print >> sys.stderr, 'STEP 170: build allHits from blastInput: mapping query, to list of species ID that are hits'

    currQuery = ""
    hitList = []
    first = True
    taxaCount = {}
    speciesCount = {}
    genusCount = {}
    familyCount = {}

    missing = open(outputPrefix + "unresolved_gi.txt",'w')
    uniqLog = open(outputPrefix + "uniqueAlignments.txt",'w')

    #temporary list to store values looked up via eUtils
    lookups = {}

    for line in open(blastInput, "r"):
        #loop control
        splits = line.split()

        if splits[0] != currQuery:
            #process stored output
            if first:
                first = False
                currQuery = splits[0]
            else:
                #process output stored in hitList
                processHitList(taxaCount, speciesCount, genusCount, familyCount, hitList, currQuery, uniqLog)

                #reset loop variables
                hitList = []
                currQuery = splits[0]

        #analysis of line - do every time, even if find end
        gi = 0
        j = 0
        k = splits[1]
        for i in k.split("|"):
            if i == "gi":
                gi = k.split("|")[j + 1]
                break
            j += 1
        #taxaID=-1
        if gi == 0 and k in fungiLookup:  #gi not found, cuz new fungiDB reference sequence name
            (taxaID, speciesID, genusID, familyID) = lookupFungi( k, fungiLookup, names)
        elif gi in taxonomy:  #2015-07-24
            taxID = taxonomy[gi]
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
                    missing.write("taxaID\t"+str(taxID)+"\n")
                    taxaID = -1
                    speciesID = -1
                    genusID = -1
                    familyID = -1
        else:
            missing.write("not in taxonomy DB\t"+line)
            taxaID = -1
            speciesID = -1
            genusID = -1
            familyID = -1
        #end if gi == 0:
        hitList.append( (taxaID, speciesID, genusID, familyID) )
    #end for
    processHitList(taxaCount, speciesCount, genusCount, familyCount, hitList, currQuery, uniqLog)
    missing.close()

    # STEP 180: output taxa reports (LCA binning version)
    print >> sys.stderr, 'STEP 180: output taxa reports (LCA binning version)'
    print >> sys.stderr, '       : writing reports'

    # 2016-12-08 need to save Kingdom information and pass forward
    savedKingdom = dict()

    filename = outputPrefix + "firstTaxon.IMSA+A_4count.txt"
    outFile = open(filename, "w")
    outFile.write( "%s\t%s\t%s\t%s\t%s\n" %
                ("Taxa ID", "Scientific Name", "IMSA+A count", "Kingdom","Total\tUnique clade hits\tPartial clade hits\tPartial clade sum"))
    KEYS=taxaCount.keys()
    KEYS.sort()
    numFirstTaxon = 0
    for taxID in KEYS:
        kingdom = 1
        if taxID == -1 or taxID == "-1":
            sname = -1
            count = taxaCount[taxID]
        elif taxID in taxaNames:    #if found in NCBI local database at systemSettings.PATH
            sname = taxaNames[taxID]
            count = taxaCount[taxID]
            kingdom = findKingdom(taxID,taxaNodes)
        else:
            # debug 2016-01-15 added step to prevent unknown/merged record to taxa report
            print "*&*&* Calling eUtils with taxID: ", taxID
            single = eutilsWrapper.getTaxa( [taxID] )
            print "*&*&* result of eUtils : ", single
            print "*&*&* result of eUtils single[0] : ", single[ single.keys()[0] ]
            if len(single)==1:
                count = taxaCount[taxID]
                try:
                    ID = int( single[ single.keys()[0] ].getSpecies().taxId )
                except:
                    print "*&*& could not look up taxaID ", taxID
                    ID = int(taxID)
                taxID = ID
                #print "*&*&* fullTax is ", single[ single.keys()[0] ]
                print "*&*&* taxID is ", taxID
                #sname = taxaNames[ taxID ]     #bug 12/06/2016
                sname = single[ single.keys()[0] ].getSpecies().sname
                print "*&*&* sname is ", sname
                kingdom = findKingdom(single[ single.keys()[0] ].getPhylum().taxId, taxaNodes)
                if single[single.keys()[0]].getSpecies():
                    savedKingdom[single[ single.keys()[0]].getSpecies().taxId ] = kingdom
                if single[single.keys()[0]].getGenus():
                    savedKingdom[single[single.keys()[0]].getGenus().taxId] = kingdom
                if single[single.keys()[0]].getFamily():
                    savedKingdom[single[single.keys()[0]].getFamily().taxId] = kingdom
            else:
                sname = -1
                count = taxaCount[taxID]

        if taxID in taxaNodes:
            level = taxaNodes[taxID][1]
        else:
            level = "unknown"
        #myLine = "%s\t%s\t%s\t%s\n" % (str(taxID), sname, level, catV(count))
        myLine = "%s\t%s\t%s\t%s\t%s\n" % (str(taxID), sname, count[1], taxaNames[kingdom], catV(count))
        if int(count[1]) > 0:
            numFirstTaxon += 1
        outFile.write(myLine)
    outFile.close()

    filename = outputPrefix + "species.IMSA+A_4count.txt"
    outFile = open(filename, "w")
    filename2 = outputPrefix + "species.IMSA_count.txt"
    outFile2 = open(filename2, "w")
    KEYS=speciesCount.keys()
    KEYS.sort()
    fillOutCountFiles(KEYS, speciesCount, taxaNames, taxaNodes, outFile, outFile2, "Species", savedKingdom)
    outFile.close()
    outFile2.close()


    filename = outputPrefix + "genus.IMSA+A_4count.txt"
    outFile = open(filename, "w")
    filename2 = outputPrefix + "genus.IMSA_count.txt"
    outFile2 = open(filename2, "w")
    KEYS=genusCount.keys()
    KEYS.sort()
    fillOutCountFiles(KEYS, genusCount, taxaNames, taxaNodes, outFile, outFile2, "Genus", savedKingdom)
    outFile.close()
    outFile2.close()


    filename = outputPrefix + "family.IMSA+A_4count.txt"
    outFile = open(filename, "w")
    filename2 = outputPrefix + "family.IMSA_count.txt"
    outFile2 = open(filename2, "w")

    KEYS=familyCount.keys()
    KEYS.sort()

    fillOutCountFiles(KEYS, familyCount, taxaNames, taxaNodes, outFile, outFile2, "Family", savedKingdom)

    outFile.close()
    outFile2.close()


    print >> sys.stderr, 'STEP 190: output human readable report'
    # We can read the files we have written, sort, and spit back out with proper justification
    # make a function ?
    # we can check file size etc etc and look for errors
    numberTaxa = dict()
    outBig = open(outputPrefix + ".IMSA+A.HUMAN_READABLE_REPORT.txt","w")

    numSpecies = 0
    numGenus = 0
    numFamily = 0
    for key in speciesCount.keys():
        if int(speciesCount[key][1]) > 0:
            numSpecies += 1

    for key in genusCount.keys():
        if int(genusCount[key][1]) > 0:
            numGenus += 1

    for key in familyCount.keys():
        if int(familyCount[key][1]) > 0:
            numFamily += 1

    outBig.write("IMSA+A metataxonomics report\n\nSUMMARY\n")
    outBig.write("{:10}".format((str(numFirstTaxon)))+"Total unique lowest taxa identified\n")
    outBig.write("{:10}".format(str(numSpecies)) + "Total unique Species identified\n")
    outBig.write("{:10}".format(str(numGenus)) + "Total unique Genera identified\n")
    outBig.write("{:10}".format(str(numFamily)) + "Total unique Families identified\n")
    outBig.write("\nCONTENT\n")
    outBig.write("Section A - List of identified lowest taxon (for Virus detection)\n")
    outBig.write("Section B - List of identified species\n")
    outBig.write("Section C - List of identified genera (recommended)\n")
    outBig.write("Section D - List of identified families\n")
    outBig.write("Section E - Errors from analysis requiring human intervention\n")
    outBig.write("Section F - List of output files and their brief description\n")

    outBig.write("\nSection A\n\n")
    outBig.write("Taxon at the lowest clade associated with the reference sequence.\n"+
                 "This report is most useful for detecting viruses, which can be omitted in other reports.\n")
    numberTaxa["firstTaxon"] = processOutputFileHumanReadable( outBig, outputPrefix + "firstTaxon.IMSA+A_4count.txt" )

    outBig.write("\nSection B\n\n")
    outBig.write("List of identified species.\n\n")
    numberTaxa["species"] = processOutputFileHumanReadable( outBig, outputPrefix + "species.IMSA+A_4count.txt" )


    outBig.write("\nSection C\n\n")
    outBig.write("List of identified genera.\n\n")
    numberTaxa["genus"] = processOutputFileHumanReadable( outBig, outputPrefix + "genus.IMSA+A_4count.txt" )

    outBig.write("\nSection D\n\n")
    outBig.write("List of identified families.\nThis report can be useful when the sequenced organism(s) in the sample is not in the reference database.\n")
    outBig.write("Reviewing results at higher clade levels may prevent the identification of many closely related organisms.\n\n")
    numberTaxa["family"] = processOutputFileHumanReadable(outBig, outputPrefix + "family.IMSA+A_4count.txt")

    # outBig.write("\n\nNumber of taxa found by clade\n")
    # clades = ["firstTaxon","species","genus","family"]
    #
    # for key in clades:
    #     outBig.write("{:15}".format(numberTaxa[key])+ " " + key + "\n")


    outBig.write("\nSection E\n\nErrors from analysis requiring human intervention\n\n")
    f = open(outputPrefix + "unidentifiedTaxaAlignments.txt")
    data = f.read()
    f.close()
    outBig.write("Sequence alignment conversion to taxa:\n")
    if len(data) != 0:
        outBig.write("WARNING!\n"+
                     "There were reference sequences which could not be converted to taxa.\n\n"+
                     outputPrefix + "unidentifiedTaxaAlignments.txt" + " file contains blast alignments to reference sequence names.\n\n"+
                     "These names must be manually looked up, input into appropriate files, and then the program re-run.\n"+
                     "Please refer to documentation for detailed instructions to resolve this error, section 'Additional Output file'.\n\n")
    else:
        outBig.write("SUCCESSFUL\n\n")

    outBig.write("Reference sequence names conversion to taxa:\n")
    f = open(outputPrefix + "unresolved_gi.txt")
    data = f.read()
    f.close()

    if len(data) != 0:
        outBig.write("WARNING!\n"+
                     "There were reference sequence GI numbers which could not be converted to taxa.\n\n" +
                     outputPrefix + "unresolved_gi.txt"+ " contains a list of GI numbers.\n\n"+
                     "These need to be manually looked up, input, and then process re-run\n"+
                     "Please refer to documentation for detailed instructions to resolve this error, section 'Additional Output file'.\n")
    else:
        outBig.write("SUCCESSFUL\n")


    outBig.write("\nSection F\n\n")
    outBig.write("\nList of output files and their brief description:\n")
    outBig.write("\t" + outputPrefix +          "species.IMSA_count.txt - Original IMSA report" + "\n")
    outBig.write("\t" + outputPrefix +          "genus.IMSA_count.txt   - Original IMSA report" + "\n")
    outBig.write("\t" + outputPrefix +          "family.IMSA_count.txt  - Original IMSA report"+"\n")
    outBig.write("\t" + outputPrefix +          "firstTaxon.IMSA+A_4count.txt - IMSA+A detailed counts; use for further analysis" + "\n")
    outBig.write("\t" + outputPrefix +          "species.IMSA+A_4count.txt    - IMSA+A detailed counts; use for further analysis" + "\n")
    outBig.write("\t" + outputPrefix +          "genus.IMSA+A_4count.txt      - IMSA+A detailed counts; use for further analysis" + "\n")
    outBig.write("\t" + outputPrefix +          "family.IMSA+A_4count.txt     - IMSA+A detailed counts; use for further analysis" + "\n")
    outBig.write("\t" + outputPrefix +          "NAMES.txt         - list of taxonomies with names looked up via NCBI" + "\n")
    outBig.write("\t" + outputPrefix +          "GIS.txt           - list of GI numbers (or other sequence names) found" + "\n")
    outBig.write("\t" + outputPrefix +          "uniqueAlignments.txt           - query names resulting in unique hits across all clades and corresponding Taxon ID" + "\n")
    outBig.write("\t" + outputPrefix +          "unidentifiedTaxaAlignments.txt - serious errors where reference sequence could not be converted to a taxon ID recorded here" + "\n")
    outBig.write("\t" + outputPrefix +          "unresolved_gi.txt              - serious errors where reference sequence could not be converted to a taxon ID recorded here" + "\n")
    outBig.write("\t" + outputPrefix[:-4] + "taxonomy.pickle       - intermediate python binary file.  Delete if your blast alignments change." + "\n")
    outBig.write("\t" + outputPrefix[:-4] + "bestHits.bln          - intermediate blast alignments, which are used for final counts.  Delete if your blast alignments change." + "\n")
    outBig.write("\t" + outputPrefix[:-4] + "bestHits.bln.pickle   - intermediate python binary file.  Delete if your blast alignments change." + "\n")

    outBig.close()

# original IMSA code modified to handle 4-tuple
def readBlastIntoTargetDictJWC(blastInput, keepQueries=False):
    """Reads the blast to create a dictionary where the key is the target and the value
    is the count of the number of hits for that target (1 if it's unique, could be 0.5 if the
    query hits two targets equally and the file has been run through the keepAllBlastHits
    filter to add the counts column.).  If keepQueries is true than the queries that hit
    the target are saved instead of the count value."""

    targets = {}
    k=0
    for line in open(blastInput):
        k+=1
        if k%100000 == 0:
            print >> sys.stderr, "Processing ", float(k)/1000000, " million"
        pieces = line.split()
        if len(pieces) < 12 or len(pieces) > 13:
            print "Ignoring line, it has %s pieces:" % len(pieces)
            print line
            continue
        elif len(pieces) == 12:
            hitVal = 1
        elif len(pieces) == 13:
            hitVal = float(pieces[12])

        #2015-07-10 : return query and hit
        target = (pieces[0], pieces[1])
        if not targets.has_key(target):
            if keepQueries:
                targets[target] = []
            else:
                targets[target] = (0, 0, 0, 0)
        """
        new tuple definition
        totalsum, # unique hits, # partial hits, partial sum
        """

        if keepQueries:
            targets[target].append(pieces[0])
        else:
            if hitVal == 1:
                targets[target] = (targets[target][0] + hitVal,
                                   targets[target][1] + 1,
                                   targets[target][2],
                                   targets[target][3])
            else:
                targets[target] = (targets[target][0] + hitVal,
                                   targets[target][1],
                                   targets[target][2] + 1,
                                   targets[target][3] + hitVal)

    return targets


##############################################
if __name__ == "__main__":
    sys.exit(main(None))

            
