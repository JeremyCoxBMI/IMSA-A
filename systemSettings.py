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

# LOCAL DATABASES
PATH_TO_OASES="/path/to/oases/"
OASES_SCRIPT=PATH_TO_OASES+"/oases/scripts/oases_pipeline.py -m 19 -M 33 -o singleEnd -d "

# NCBI database
PATH="/data/ncbi/gi_taxid/"
BLAST_TAX_DB = PATH+"gi_taxid_nucl.dmp"
NAMES_BLAST_TAX_DB = PATH+"names.dmp"
NODES_BLAST_TAX_DB = PATH+"nodes.dmp"
EXTRA_BLAST_TAX_DB = PATH+"extra_JWC_gi_taxid_nucl.dmp"

# custom database index; augments lookup capability to include select genomes from fungiDB
# this file establishes taxonomy lookup
fungiLookupFile = "/home/username/imsa_v2/fungiDB.seqNames.2.speciesID"

def buildNames():
    inFile = open(NAMES_BLAST_TAX_DB)

    result = dict()

    for line in inFile:
        line = line.replace("\t","")
        line = line.replace("|\n","")
        splits = line.split("|")
        if splits[-1] == 'scientific name':
            key = int( line.split("|")[0] )
            name = line.split("|")[1]
            result[key] = name

    return result

def buildNodes():
    inFile = open(NODES_BLAST_TAX_DB)

    result = dict()
    levels = dict()

    for line in inFile:
        #print line
        line = line.replace("\t","")
        line = line.replace("|\n","")
        splits = line.split("|")
        #print splits
        taxaID = int(splits[0])
        parentID = int(splits[1])
        level = splits[2]
        result[taxaID] = (parentID, level)
        levels[level]=1

    return (result, levels)

def buildReverseNames():
    inFile = open(PATH+"names.dmp")

    result = dict()


    for line in inFile:
        line = line.replace("\t","")
        line = line.replace("|\n","")
        splits = line.split("|")
        if splits[-1] == 'scientific name':
            key = int( line.split("|")[0] )
            name = line.split("|")[1]
            result[name] = key

    return result

def buildReverseAndForwardNames():
    inFile = open(NAMES_BLAST_TAX_DB)

    result = dict()
    reverse = dict()

    for line in inFile:
        line = line.replace("\t","")
        line = line.replace("|\n","")
        splits = line.split("|")
        if splits[-1] == 'scientific name':
            key = int( line.split("|")[0] )
            name = line.split("|")[1]
            result[key] = name
            reverse[name]=key

    return (result, reverse)