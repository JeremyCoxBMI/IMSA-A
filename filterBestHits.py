#!/usr/local/bin/python

# This script filters the best hit blast results (the file with the score on the end of
# each line) to only include scores above a certain threshold.  This is done in an
# intelligent way, looking at the score within a species, so hits are not filtered if
# they hit multiple records of the same species.

# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)


import sys
import eutilsWrapper
import config

try:
    bestHitName = sys.argv[1]
except:
    print "Usage: python filterBestHits.py bestHitBlastResults [threshold] [outputName]"
    sys.exit(1)
    
try:
    threshold = float(sys.argv[2])
except:
    print "No thresold entered.  Using the default value of 0.05"
    threshold = 0.05

try:
    outputName = sys.argv[3]
except:
    outputName = "%s.filtered%s.bln" % (bestHitName[:-4], threshold)
    print "No output name entered.  Using: %s" % (outputName)


# old way: filter just based on the score at the end of the line
# new way: look up the tax id using the gi then sum across tax id before filtering

gis = {}
for line in open(bestHitName):
    target = line.split()[1].split("|")[1]
    gis[target] = 1

taxonomy = {}
print "Getting taxonomy ids from local file for %s gis" % (len(gis))
for line in open(config.BLAST_TAX_DB):
    (gi, taxid) = line.split()
    if gis.has_key(gi):
        taxonomy[gi] = int(taxid)

fullTax = eutilsWrapper.getTaxa(taxonomy.values())
                                    
print "filtering blast results based on tax id"
countUsed = 0
countSkipped = 0
out = open(outputName, "w")
currentQuery = None
currentResults = []
notFoundTax = {}
countFoundTax = 0
for line in open(bestHitName):
    query = line.split()[0]

    if currentQuery == query:
        currentResults.append(line)
    else:
        #print
        #print
        allTaxScores = {}
        for r in currentResults:
            score = float(r.split()[-1])
            try:
                gi = r.split()[1].split("|")[1]
                tax = taxonomy[gi]
                species = fullTax[tax].getSpecies().taxId
                countFoundTax += 1
            except KeyError:
                tax = 0
                #print "Unable to find taxonomy for gi %s" % (gi)
                if not notFoundTax.has_key(gi):
                    notFoundTax[gi] = 0
                notFoundTax[gi] += 1
                species = 0
            except:
                species = tax
                
            if not allTaxScores.has_key(species):
                allTaxScores[species] = 0
            if species != 0:
                allTaxScores[species] += score

            #print gi, tax, species, score

        #print allTaxScores

        for r in currentResults:
            try:
                gi = r.split()[1].split("|")[1]
                tax = taxonomy[gi]
                species = fullTax[tax].getSpecies().taxId
            except KeyError:
                tax = 0
                species = 0
            except:
                species = tax
            
            if allTaxScores[species] >= threshold:
                out.write(r)
                countUsed += 1
                #print "Used", gi, tax
            else:     
                countSkipped += 1
                #print "Filtered", gi, tax

        currentQuery = query
        currentResults = [line]

print "Found %s taxonomy results." % (countFoundTax)
print "The following gis didn't have known taxonomies (number of times the gi came up also reported.)"
for k,v in notFoundTax.iteritems():
    print k, v 
print "Done.  Filtered at %s, Kept %s records but skipped %s" % (threshold, countUsed, countSkipped)
