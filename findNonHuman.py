#!/usr/local/bin/python
#
# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)


import sys
import sequenceUtils as su

import config

# read through the blast file to get all the gis
args = sys.argv

try:
    inputBlast = sys.argv[1]
except:
    print "Error reading arguments!"
    print "Usage:  python findNonHuman.py blastResults"
    sys.exit(1)

gis = {}
for line in open(inputBlast):
    pieces = line.split()
    gis[pieces[1].split("|")[1]] = 1
                            

# get the taxonomy for each gi
taxonomy = {}
print "Getting taxonomy ids from local file for %s gis" % (len(gis))
for line in open(config.BLAST_TAX_DB):
    (gi, taxid) = line.split()
    if gis.has_key(gi):
        taxonomy[gi] = int(taxid)

# read through the blast again.  For all hits to one read, see if there is 1 human (or close).
# If not, print read name.
humanSpecies = [9606, 9598, 9544, 9597, 9601]
hpvGuess = [333760, 337041, 10566]
lastRead = ""
currentSpecies = []
isHuman = False
readsNotHuman = {}
countHpv = 0
for line in open(inputBlast):
    pieces = line.split()
    if pieces[0] != lastRead:
        # okay, interpret what we've got
        if not isHuman: 
            print "Read %s only found in species %s" % (lastRead, currentSpecies)
            readsNotHuman[lastRead] = 1
            for sp in hpvGuess:
                if sp in currentSpecies:
                    countHpv += 1
                    break

        # set up for the next round
        lastRead = pieces[0]
        currentSpecies = []
        isHuman = False

    gi = pieces[1].split("|")[1]
    try:
        tax = taxonomy[gi]
        currentSpecies.append(tax)
        if tax in humanSpecies:
            isHuman = True
    except:
        print "Unable to find a species for %s for read %s." % (pieces[1], pieces[0])
        currentSpecies.append(None)

print "There were %s non-human reads." % (len(readsNotHuman))
print "The non-human reads are:"
for k,v in readsNotHuman.iteritems():
    print k

print "There were %s HPV non-human reads" % (countHpv)
