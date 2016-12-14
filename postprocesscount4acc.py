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
from postprocess import doPostProcess

"""
JWC (Jeremy.Cox@cchmc.org) modifying this code to do different postprocessing.

THIS SCRIPT IS POSTPROCESSOR UPGRADES ASSOCIATED WITH IMSA+A
However, the postprocessor can be used with regular IMSA pipeline results.

This code is modified version of postprocesscount4.py.  The previous used GI numbers.  This script now uses Accession numbers instead.

Note that code base is kept the same, only some functions are overriden to use accession numbers.

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

#trying to keep code in one location; overrides included here
from postprocesscount4 import *


HELP_STRING = """Given the results of the filtered reads blasted against reference database, this script
does the best hit blast filter and then runs the taxonomy reports.  Uses Python 2.7.

Required Inputs:
    -b    Blast output file

Optional Inputs:
    -f    this option is deprecated for IMSA+A 4 count.  (It may still function)
          This program now outputs file with extension .firstTaxon.IMSA+A_4count.txt for better virus classification.

REQUIRES PYTHON 2.7
"""


# global variable
blastResults = ""

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


#new function 12-06-2016
def parseAccession(text, fungiLookup):
    acc = "0"
    j = 0
    if text in fungiLookup:
        acc = "0"
    elif text[0] == "g" and text[1] == "i":  ##old school format
        for i in text.split("|"):
            if i == "ref":
                acc = text.split("|")[j + 1].split(".")[0]  # exclude .suffix

                break
            j += 1
    else:  ##new format (I presume... )
        acc = text.split(".")[0]



    return acc

# new function
def extractLookupInfo(accs, sID_to_lookup, fungiLookup, blastInput):
    """ This function pulls gi numbers out of blast input file.
        Gi numbers and speciesID for non-gi reference sequences are returned.
    :param accs:            (dict) OUTPUT: list of accession numbers (as text) found in blastinput file
    :param sID_to_lookup:   (dict) OUTPUT: list of speciesID we need to lookup in NCBI
    :param fungiLookup:     (dict) INPUT: mapping ref sequence names (non-gi) to speciesID
    :param blastInput:      (string) INPUT: name of file containing blast hits
    :return:
    """
    input = open(blastInput,"r")

    for line in input:
        #parse blast lines for reference sequence names
        splits = line.split('\t')
        acc = parseAccession(splits[1], fungiLookup)
        if acc == "0":  # not found, try lookup in additional seqname/speciesID list (fungiLookup)
            if splits[1] in fungiLookup:
                sID = fungiLookup[splits[1]]
                sID_to_lookup[sID] = 1
        else:
            accs[acc] = 1


#new function
def findingUnknownAlignments( taxonomy, blastInput, filenameOUT, fungiLookup):
    """ If a blast file and looks for any gi codes, which were not found, then writes to a file.
        Can be used by USER to update the systemSettings.EXTRA_BLAST_TAX_DB file.
    :param taxonomy:    (dict) INPUT: maps accession numbers found in database to their taxon ID
                                            used here as a list of accession numbers found
    :param blastInput:  (string) INPUT: name of input blast file
    :param filenameOUT: (string) INPUT: name of file to write entries to
    :return:
    """

    outF = open(filenameOUT,"w")
    for line in open(blastInput,"r"):
        #parse blast lines for reference sequence gi numbers
        splits = line.split('\t')

        acc = parseAccession(splits[1], fungiLookup)

        #report on gi number not found
        if acc != "0" and not acc in taxonomy:
            #acc can not be non-acc fungiDB code, because parseAccession checked; would have assigned a "0"
            #if not found earlier in taxadump files (thus not in taxonomy)
            outF.write(line)
    outF.close()



# original IMSA code modified to handle 4-tuple and new counting scheme and upgrades
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
    accs - dictionary mapping                (used as non-repeating list)
        accesion numbers (without.version suffix) (string) - referenceces we are seeking to look up

    taxID_to_lookup - dictionary mapping    (used as non-repeating list)
        taxID (int) - taxon ID to lookup
        to 1
    taxonomy - dictionary mapping
        ###gid (string)
        accession number (without .version suffix) (string)
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
    print >> sys.stderr, 'STEP 120: extract accession number targets, SID_to_lookup from file'

    #contains integer taxID to lookup in NCBI -- will contain FungiDB species and regular database taxID
    taxID_to_lookup = {}

    # dict containing gi's to lookup, requiring NCBI lookup
    #gis = {}   #removed for
    accs = {}

    #extractLookupInfo(gis, taxID_to_lookup, fungiLookup, blastInput)
        # gis is indexed by string numbers
        # taxID_to_lookup is indexed by integer numbers

    extractLookupInfo(accs, taxID_to_lookup, fungiLookup, blastInput)
        # gis is indexed by string numbers
        # taxID_to_lookup is indexed by integer numbers



    # STEP 130: write gis table to file
    print >> sys.stderr, 'STEP 130: write accession numbers table to file'
    outGis = open(outputPrefix + "ACCS.txt", "w")
    for key in accs:
        outGis.write("%s\n" % (key))
    outGis.close()

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
                (acc, accv, taxid, gi) = line.split()
                if accs.has_key(acc):
                    taxonomy[acc] = int(taxid)

        print >> sys.stderr, '       : creating pickle on disk'
        pickle.dump( taxonomy, open( blastInput+".taxonomy.pickle", "wb" ) )
    else:
        print >> sys.stderr, '       : loading pickle'
        taxonomy = pickle.load(  open( blastInput+".taxonomy.pickle", "rb" )  )

    print >> sys.stderr, "       : Done getting taxonomy ids, got %s ids" % (len(taxonomy))


    # STEP 145: write Alignments to gi nubmers, which were not found
    print >> sys.stderr, 'STEP 145: write Alignments to accession numbers, which were not found'
    findingUnknownAlignments( taxonomy, blastInput, outputPrefix + "unidentifiedTaxaAlignments.txt", fungiLookup)


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
    print >> sys.stderr, 'STEP 168: import ncbi taxonomy database'
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

    missing = open(outputPrefix + "unresolved_acc.txt",'w')
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
        k = splits[1]
        acc = parseAccession(k, fungiLookup)

        if acc == "0" and k in fungiLookup:  #gi not found, cuz new fungiDB reference sequence name
            (taxaID, speciesID, genusID, familyID) = lookupFungi( k, fungiLookup, names)
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
    missing.close()

    # STEP 180: output taxa reports (LCA binning version)
    print >> sys.stderr, 'STEP 180: output taxa reports (LCA binning version)'
    print >> sys.stderr, '       : writing reports'

    # 2016-12-08 need to save Kingdom information and pass forward
    savedKingdom = dict()

    filename = outputPrefix + "firstTaxon.IMSA+A_4count.txt"
    outFile = open(filename, "w")
    outFile.write( "%s\t%s\t%s\n" %
                ("Taxa ID", "Scientific Name", "Clade Level\tTotal\tUnique clade hits\tPartial clade hits\tPartial clade sum") )
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
            single =  eutilsWrapper.getTaxa( [taxID] )
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
    KEYS = speciesCount.keys()
    KEYS.sort()
    fillOutCountFiles(KEYS, speciesCount, taxaNames, taxaNodes, outFile, outFile2, "Species", savedKingdom)
    outFile.close()
    outFile2.close()

    filename = outputPrefix + "genus.IMSA+A_4count.txt"
    outFile = open(filename, "w")
    filename2 = outputPrefix + "genus.IMSA_count.txt"
    outFile2 = open(filename2, "w")
    KEYS = genusCount.keys()
    KEYS.sort()
    fillOutCountFiles(KEYS, genusCount, taxaNames, taxaNodes, outFile, outFile2, "Genus", savedKingdom)
    outFile.close()
    outFile2.close()

    filename = outputPrefix + "family.IMSA+A_4count.txt"
    outFile = open(filename, "w")
    filename2 = outputPrefix + "family.IMSA_count.txt"
    outFile2 = open(filename2, "w")

    KEYS = familyCount.keys()
    KEYS.sort()

    fillOutCountFiles(KEYS, familyCount, taxaNames, taxaNodes, outFile, outFile2, "Family", savedKingdom)

    outFile.close()
    outFile2.close()

    print >> sys.stderr, 'STEP 190: output human readable report'
    # We can read the files we have written, sort, and spit back out with proper justification
    # make a function ?
    # we can check file size etc etc and look for errors
    numberTaxa = dict()
    outBig = open(outputPrefix + ".IMSA+A.HUMAN_READABLE_REPORT.txt", "w")

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
    outBig.write("{:10}".format((str(numFirstTaxon))) + "Total unique lowest taxa identified\n")
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
    outBig.write("Taxon at the lowest clade associated with the reference sequence.\n" +
                 "This report is most useful for detecting viruses, which can be omitted in other reports.\n")
    numberTaxa["firstTaxon"] = processOutputFileHumanReadable(outBig, outputPrefix + "firstTaxon.IMSA+A_4count.txt")

    outBig.write("\nSection B\n\n")
    outBig.write("List of identified species.\n\n")
    numberTaxa["species"] = processOutputFileHumanReadable(outBig, outputPrefix + "species.IMSA+A_4count.txt")

    outBig.write("\nSection C\n\n")
    outBig.write("List of identified genera.\n\n")
    numberTaxa["genus"] = processOutputFileHumanReadable(outBig, outputPrefix + "genus.IMSA+A_4count.txt")

    outBig.write("\nSection D\n\n")
    outBig.write(
        "List of identified families.\nThis report can be useful when the sequenced organism(s) in the sample is not in the reference database.\n")
    outBig.write(
        "Reviewing results at higher clade levels may prevent the identification of many closely related organisms.\n\n")
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
        outBig.write("WARNING!\n" +
                     "There were reference sequences which could not be converted to taxa.\n\n" +
                     outputPrefix + "unidentifiedTaxaAlignments.txt" + " file contains blast alignments to reference sequence names.\n\n" +
                     "These names must be manually looked up, input into appropriate files, and then the program re-run.\n" +
                     "Please refer to documentation for detailed instructions to resolve this error, section 'Additional Output file'.\n\n")
    else:
        outBig.write("SUCCESSFUL\n\n")


    f = open(outputPrefix + "unresolved_acc.txt")
    data = f.read()
    f.close()
    outBig.write("Reference sequence names conversion to taxa:\n")
    if len(data) != 0:
        outBig.write("WARNING!\n" +
                     "There were reference sequence GI numbers which could not be converted to taxa.\n\n" +
                     outputPrefix + "unresolved_acc.txt" + " contains a list of Accession numbers.\n\n" +
                     "These need to be manually looked up, input, and then process re-run\n" +
                     "Please refer to documentation for detailed instructions to resolve this error, section 'Additional Output file'.\n")
    else:
        outBig.write("SUCCESSFUL\n")

    outBig.write("\nSection F\n\n")
    outBig.write("\nList of output files and their brief description:\n")
    outBig.write("\t" + outputPrefix + "species.IMSA_count.txt - Original IMSA report" + "\n")
    outBig.write("\t" + outputPrefix + "genus.IMSA_count.txt   - Original IMSA report" + "\n")
    outBig.write("\t" + outputPrefix + "family.IMSA_count.txt  - Original IMSA report" + "\n")
    outBig.write(
        "\t" + outputPrefix + "firstTaxon.IMSA+A_4count.txt - IMSA+A detailed counts; use for further analysis" + "\n")
    outBig.write(
        "\t" + outputPrefix + "species.IMSA+A_4count.txt    - IMSA+A detailed counts; use for further analysis" + "\n")
    outBig.write(
        "\t" + outputPrefix + "genus.IMSA+A_4count.txt      - IMSA+A detailed counts; use for further analysis" + "\n")
    outBig.write(
        "\t" + outputPrefix + "family.IMSA+A_4count.txt     - IMSA+A detailed counts; use for further analysis" + "\n")
    outBig.write("\t" + outputPrefix + "NAMES.txt         - list of taxonomies with names looked up via NCBI" + "\n")
    outBig.write("\t" + outputPrefix + "ACCS.txt          - list of Accession numbers (or other sequence names) found" + "\n")
    outBig.write(
        "\t" + outputPrefix + "uniqueAlignments.txt           - query names resulting in unique hits across all clades and corresponding Taxon ID" + "\n")
    outBig.write(
        "\t" + outputPrefix + "unidentifiedTaxaAlignments.txt - serious errors where reference sequence could not be converted to a taxon ID recorded here" + "\n")
    outBig.write(
        "\t" + outputPrefix + "unresolved_acc.txt             - serious errors where reference sequence could not be converted to a taxon ID recorded here" + "\n")
    outBig.write("\t" + outputPrefix[
                        :-4] + "taxonomy.pickle       - intermediate python binary file.  Delete if your blast alignments change." + "\n")
    outBig.write("\t" + outputPrefix[
                        :-4] + "bestHits.bln          - intermediate blast alignments, which are used for final counts.  Delete if your blast alignments change." + "\n")
    outBig.write("\t" + outputPrefix[
                        :-4] + "bestHits.bln.pickle   - intermediate python binary file.  Delete if your blast alignments change." + "\n")

    outBig.close()

# # original IMSA code modified to handle 4-tuple
# def readBlastIntoTargetDictJWC(blastInput, keepQueries=False):
#     """Reads the blast to create a dictionary where the key is the target and the value
#     is the count of the number of hits for that target (1 if it's unique, could be 0.5 if the
#     query hits two targets equally and the file has been run through the keepAllBlastHits
#     filter to add the counts column.).  If keepQueries is true than the queries that hit
#     the target are saved instead of the count value."""
#
#     targets = {}
#     k=0
#     for line in open(blastInput):
#         k+=1
#         if k%100000 == 0:
#             print >> sys.stderr, "Processing ", float(k)/1000000, " million"
#         pieces = line.split()
#         if len(pieces) < 12 or len(pieces) > 13:
#             print "Ignoring line, it has %s pieces:" % len(pieces)
#             print line
#             continue
#         elif len(pieces) == 12:
#             hitVal = 1
#         elif len(pieces) == 13:
#             hitVal = float(pieces[12])
#
#         #2015-07-10 : return query and hit
#         target = (pieces[0], pieces[1])
#         if not targets.has_key(target):
#             if keepQueries:
#                 targets[target] = []
#             else:
#                 targets[target] = (0, 0, 0, 0)
#         """
#         new tuple definition
#         totalsum, # unique hits, # partial hits, partial sum
#         """
#
#         if keepQueries:
#             targets[target].append(pieces[0])
#         else:
#             if hitVal == 1:
#                 targets[target] = (targets[target][0] + hitVal,
#                                    targets[target][1] + 1,
#                                    targets[target][2],
#                                    targets[target][3])
#             else:
#                 targets[target] = (targets[target][0] + hitVal,
#                                    targets[target][1],
#                                    targets[target][2] + 1,
#                                    targets[target][3] + hitVal)
#
#     return targets


##############################################
if __name__ == "__main__":
    sys.exit(main(None))

            
