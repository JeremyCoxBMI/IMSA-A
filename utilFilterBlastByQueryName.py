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


###################################################
# Utility program:
#   Take a list of query names and filter the blast results to se
#
#  EXAMPLE:
#     grep -P "\s469$" 4A8F_*uniqueAlignments.txt | cut -f1 > actinobacter.query.unique.hits.txt
#     python utilFilterBlastByQueryName.py actinobacter.query.unique.hits.txt 4A8F_sub_assemble.fa.bln  4A8F_sub_assemble.fa.bln.filter.469.bln
#
###################################################

import sys



if __name__ == "__main__":
    hits2filter = sys.argv[1]
    blast6 = sys.argv[2]
    output = sys.argv[3]

    alignNames = dict()

    for line in open(hits2filter):
        alignNames[line.strip()] = 1

    outF = open(output,"w")

    write = False
    for line in open(blast6):
        if line.strip().split("\t")[0] in alignNames:
            outF.write(line)