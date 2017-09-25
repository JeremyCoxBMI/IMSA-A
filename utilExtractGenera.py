
from systemSettings import *
import sys

STOP = ["family","order","class","phylum","kingdom","superkingdom"]

if __name__ == "__main__":


    print >> sys.stderr, "Importing names and nodes"
    (nodes, levels) = buildNodes()
    names = buildNames()

    print >> sys.stderr, "reading input file ", sys.argv[1]
    print >> sys.stderr, "writing output file ", sys.argv[2]

    acc_to_lookup = {}

    print >> sys.stderr, "Parsing Input into Memory"
    for line in open(sys.argv[1]):
        if line.split("/")[0][0] != 'n':  #skip empty lines
            file = line.split("/")[10]
            accession = "_".join(file.split("_")[0:2]).split(".")[0]
            acc_to_lookup[accession] = line

    print >> sys.stderr, "Converting ACC to TaxonID"
    acc_to_taxid = {}
    for line in open(ACC_BLAST_TAX_DB):
        (acc, accv, taxid, gi) = line.split()
        acc = acc.replace("_","")
        if acc_to_lookup.has_key(acc):
            acc_to_taxid[acc] = int(taxid)

    print >> sys.stderr, "Writing output: filename to species and genus"
    outF = open(sys.argv[2],'w')
    outF.write("file\tspecies ID\tspecies Name\tgenus ID\tgenus name\n")

    for key in acc_to_lookup.keys():

        taxid = acc_to_lookup[key]

        if key in acc_to_taxid:
            file = acc_to_taxid[key]
        else:
            file = key+"\tERROR: no file name"

        speciesID = -1
        genusID = -1

        while(True):
            if taxid in nodes:
                (parent, level) = nodes[taxid]
                if level == "species":
                    speciesID = taxid
                elif level == "genus":
                    genusID = taxid
                    break
                elif level in STOP:
                    break
                taxid = parent
            else:
                break

        if speciesID == -1:
            speciesName = "unknown"
        else:
            speciesName = names[speciesID]

        if genusID == -1:
            genusName = "unknown"
        else:
            genusName = names[genusID]

        outF.write(file+"\t"+str(speciesID)+"\t"+speciesName+"\t"+str(genusID)+"\t"+str(genusName)+"\n")