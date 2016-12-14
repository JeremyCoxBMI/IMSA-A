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

#GID based alignment database (old / reverse compatible)
PATH="/mnt/Dshare/gi_taxid/"
BLAST_TAX_DB = PATH+"gi_taxid_nucl.dmp"
NAMES_BLAST_TAX_DB = PATH+"names.dmp"
NODES_BLAST_TAX_DB = PATH+"nodes.dmp"

# This file allows you to manually define Accesion numbers to Taxon ID conversion, if it somehow was left out of
# the database dump, and ISMA then throws an error.
# Add to this tab-delimited file lines matching the format of the BLAST_TAX_DB file:
# 2 columns: GI number, taxonID
# ACTION ITEM: you need to create an empty file until there is data to put in this file
EXTRA_BLAST_TAX_DB = PATH+"extra_gi_taxid_nucl.dmp"             #note, you may need to create an empty file

#ACCESSION number based alignments database (new)
ACC_PATH="/mnt/Dshare/acc_taxid/"
ACC_BLAST_TAX_DB = ACC_PATH+"nucl_gb.accession2taxid"
ACC_NAMES_BLAST_TAX_DB = PATH+"names.dmp"
ACC_NODES_BLAST_TAX_DB = PATH+"nodes.dmp"

# This file allows you to manually define Accesion numbers to Taxon ID conversion, if it somehow was left out of
# the database dump, and ISMA then throws an error.
# Add to this tab-delimited file lines matching the format of the ACC_BLAST_TAX_DB file:
# 4 columns: Accession Number, Accession.Version, taxonID, GI number
# Note the GI number can be a blank column, but the column needs to be present in the file
# ACTION ITEM: you need to create an empty file until there is data to put in this file
ACC_EXTRA_BLAST_TAX_DB = ACC_PATH+"extra_acc_taxid_nucl.dmp"

# custom database index; augments lookup capability to include select genomes from fungiDB
# this file establishes taxonomy lookup
# User can add more codes here, format is tab-delimited, 1 line per entry
# 2 columns: sequenceName, taxon ID
fungiLookupFile = "/home/osboxes/imsa.develop/imsa-a/fungiDB.seqNames.2.speciesID"

# For running nuccore utilities, eDirect must be installed
# https://www.ncbi.nlm.nih.gov/books/NBK179288/
# ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.zip
efetchPath="~/EDirect/edirect/efetch"

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