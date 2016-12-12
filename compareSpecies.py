#!/usr/local/bin/python
#
# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)

import sys

args = sys.argv

try:
    species1 = sys.argv[1]
    species2 = sys.argv[2]
except:
    print "Error reading arguments!"
    print "Usage:  python compareSpecies.py file1 file2 [threshold]"
    sys.exit(1)

threshold = 0.05
try:
    threshold = float(sys.argv[3])
except:
    print "The threshold is the default value of %s" % (threshold)

# note:  some of the names have changed.  The names from file1 will be used.
idToName = {}
count1 = {}
count2 = {}

for line in open(species2):
    if line.find("Count") > 0:
        continue

    pieces = line.split()
    idToName[pieces[0]] = "_".join(pieces[1:-1])
    count2[pieces[0]] = float(pieces[-1])

for line in open(species1):
    if line.find("Count") > 0:
        continue

    pieces = line.split()
    idToName[pieces[0]] = "_".join(pieces[1:-1])
    count1[pieces[0]] = float(pieces[-1])

print "The ids with different values are:"
print "ID\tName\tCount1\tCount2"
countSame = 0
countDiff = 0
for id, name in idToName.iteritems():
    c1 = c2 = 0
    if count1.has_key(id):
        c1 = count1[id]
    if count2.has_key(id):
        c2 = count2[id]

    if c1 != c2 and (c1 > threshold or c2 > threshold):
        countDiff += 1
        print "%s\t%s\t%s\t%s" % (id, name, c1, c2)
    else:
        countSame += 1

print "There were %s that were the same and %s that were different" % (countSame, countDiff)
