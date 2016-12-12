#!/usr/local/bin/python
#
# This module wraps access to the NCBI eutils.

# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)

#IMSA+A update 2016-12
#import urllib
import urllib2 as urllib

from xml.dom import minidom
import re
import StringIO

######################
# NCBI URLs
EUTILS_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ESUMMARY_URL = EUTILS_URL + "esummary.fcgi?"
EFETCH_URL = EUTILS_URL + "efetch.fcgi?"
ESEARCH_URL = EUTILS_URL + "esearch.fcgi?usehistory=y&"

BATCH_SIZE=250


class Taxonomy:
    def __init__(self, taxonNode):
        self.lineage = []
        #self.division = None
        #self.taxId = None
        #self.sname = None
        taxId = None
        division = None
        sname = None
        
        for cn in taxonNode.childNodes:
            if cn.nodeName=="TaxId":
                taxId = int(cn.childNodes[0].nodeValue)
            elif cn.nodeName=="ScientificName":
                sname = str(cn.childNodes[0].nodeValue)
            elif cn.nodeName=="Division":
                division = str(cn.childNodes[0].nodeValue)
            elif cn.nodeName=="Rank":
                rank = cn.childNodes[0].nodeValue
            elif cn.nodeName=="LineageEx":
                for ccn in cn.childNodes:
                    #print "ccn", ccn.nodeName
                    cRank=None
                    for cccn in ccn.childNodes:
                        #print "cccn", cccn.nodeName
                        if cccn.nodeName == "Rank":
                            cRank = str(cccn.childNodes[0].nodeValue)
                        elif cccn.nodeName == "ScientificName":
                            cSname = str(cccn.childNodes[0].nodeValue)
                        elif cccn.nodeName == "TaxId":
                            cTaxId = int(cccn.childNodes[0].nodeValue)
                    if not cRank:
                        continue
                    lineageTaxon = OneTaxon(cTaxId, cSname, cRank)
                    self.lineage.append(lineageTaxon)

        # make sure to append the current rank to the lineage
        if rank and taxId and sname:
            currentTaxon = OneTaxon(taxId, sname, rank)
            self.lineage.append(currentTaxon)

        self.thisTaxon = OneTaxon(taxId, sname, division)

    def getSuperKingdom(self):
        return self._tryKey("superkingdom")

    def _getLineageItem(self, key):
        """Retrieves the items with the rank equal to key from the
        lineage list.  Just returns the first one -- this function
        should be called only for items guaranteed to have only one,
        like species, not 'no rank'"""
        for t in self.lineage:
            if t.rank == key:
                return t
        return None

    def isVirus(self):
        #sk = self.getSuperKingdom()
        #for (id,val) in sk:
        #    if val == "Viruses":
        #        return True
        #return False
        return self.division == "Viruses"

    def isBacteria(self):
        #sk = self.getSuperKingdom()
        #for (id, val) in sk:
        #    if val == "Bacteria":
        #        return True
        #return False
        return self.division == "Bacteria"

    def getDivision(self):
        return self.thisTaxon.rank

    def getPhylum(self):
        return self._getLineageItem("phylum")

    def getClass(self):
        return self._getLineageItem("class")
    
    def getOrder(self):
        return self._getLineageItem("order")

    def getFamily(self):
        return self._getLineageItem("family")

    def getGenus(self):
        return self._getLineageItem("genus")

    def getSpecies(self):
        return self._getLineageItem("species")

    def getSubSpeces(self):
        return self._getLineageItem("subspecies")

    def getSubFamily(self):
        return self._getLineageItem("subfamily")


class OneTaxon:
    """This class represents one node in the lineageEx tree."""
    def __init__(self, taxId, sname, rank):
        self.taxId = taxId
        self.sname = sname
        self.rank = rank

    def __repr__(self):
        return "{{taxid=%s, sname='%s', rank=%s}}" % (self.taxId, self.sname, self.rank)


def batchList(iterable,batchSize=None,batchCount=None):
    """ return a list of lists representing batches of items from
    iterable. Either batchSize or batchCount should be an integer.
    """
    if batchCount == None and batchSize == None:
        raise ValueError , "only batchSize or batchCount may be specified"
    elif batchCount != None and batchSize != None:
        raise ValueError , "only batchSize or batchCount may be specified"

    rv=[]
    try:
        inList=list(iterable)
    except TypeError:
        # got non sequence type
        return [[iterable]]
        
    if batchSize != None:
        batchSize = int(batchSize)
    else:
        if len(inList)%batchCount==0:
            batchSize=len(inList)/batchCount
        else:
            batchSize=(len(inList)/batchCount) +1

    for start in range(0,len(inList),batchSize):
        rv.append(inList[start:start+batchSize])
    return rv


def parseMultiXML(iFile):

    XML_RE="\<\?xml.+\?\>"
    DOC_RE="\<\!DOCTYPE.+\>"

    timeToYield = False

    xmlHeader = ""
    docHeader = ""
    currentXML = ""
    for line in iFile:
        if re.search(XML_RE, line):
            xmlHeader = line
            timeToYield = True
        elif re.search(DOC_RE, line):
            docHeader = line
            timeToYield = True
        else:
            currentXML += line
            #print currentXML
            timeToYield = False

        if timeToYield:
            timeToYield = False
            if len(currentXML.strip()) == 0:
                continue
            else:
                yield xmlHeader + docHeader + currentXML
                currentXML = ""

    if len(currentXML) > 0:
        yield xmlHeader + docHeader + currentXML


def efetch(url, args, debug=False):
    """Runs the urlopen command and returns the stream."""
    
    url = url + '&'.join([urllib.quote(k)+'='+urllib.quote(str(v)) for k,v in args.items() if v!= None])

    if debug:
        print "Contacting NCBI.  Url:"
        print url

    output = urllib.urlopen(url)
    return StringIO.StringIO(output.read())

def getTaxa(ids):
    """Returns a dictionary where the key is the id and the value is the
    Taxon class object for each id in the list that can be found
    using the NCBI EUTILS."""
    d = {}
    idLists = batchList(ids, batchSize=BATCH_SIZE)
    #print idLists
    for idList in idLists:
        ids = _getThrottledTaxa(idList)
        #print len(ids)
        d.update(ids)

    return d

def _getThrottledTaxa(ids):
    args = {}
    args["db"] = "taxonomy"
    args["mode"] = "xml"
    args["id"] = ",".join([str(x) for x in ids])
    args["report"] = "brief"

    sioBuffer = efetch(EFETCH_URL, args)
    #print sioBuffer.read()
    #sioBuffer.seek(0)

    docs = []
    for xml in parseMultiXML(sioBuffer):
        try:
            doc = minidom.parse(StringIO.StringIO(xml))
        except:
            print "Unable to parse:"
            print xml
        #print doc
        #doc.seek(0)
        docs.append(doc)

    d = {}
    #print "There are %s docs in _getThrottledTaxa and there were %s ids" % (len(docs), len(ids))
    for doc in docs:
        try:
            ts = doc.getElementsByTagName("TaxaSet")[0]
        except IndexError:
            print "XML Doc returned from NCBI was empty."
            print doc
            continue
        
        for t in ts.childNodes:
            if t.nodeName == "Taxon":
                x = Taxonomy(t)
                #print thisTaxon.taxId, thisTaxon.division
                d[x.thisTaxon.taxId] = x
        
    return d

def testTaxa(inputFile):
    """Given the xml record from EUtils taxonomy db (brief), create a taxon class."""

    infile = open(inputFile)

    docs = []
    for xml in parseMultiXML(infile):
        doc = minidom.parse(StringIO.StringIO(xml))
        docs.append(doc)

    d = {}
    for doc in docs:
        t = Taxonomy(doc)
        d[t.thisTaxon.taxId] = t

    return d


def getTaxIdFromGi (ids):
    """Returns a dictionary where the key is the id and the value is the
    Taxonomy id for each id in the list that can be found
    using the NCBI EUTILS."""
    d = {}
    idLists = batchList(ids, batchSize=BATCH_SIZE)
    #print idLists
    for idList in idLists:
        ids = _getThrottledTaxIds(idList)
        #print len(ids)
        d.update(ids)

    return d

def _getThrottledTaxIds(ids):
    args = {}
    args["db"] = "nuccore"
    args["retmode"] = "xml"
    args["id"] = ",".join([str(x) for x in ids])
    args["rettype"] = "fasta"

    sioBuffer = efetch(EFETCH_URL, args)
    #print sioBuffer.read()
    #sioBuffer.seek(0)

    docs = []
    for xml in parseMultiXML(sioBuffer):
        try:
            doc = minidom.parse(StringIO.StringIO(xml))
        except:
            print "Unable to parse:"
            print xml
        #print doc
        #doc.seek(0)
        docs.append(doc)

    d = {}
    for doc in docs:
        try:
            ts = doc.getElementsByTagName("TSeqSet")[0]
        except IndexError:
            print "XML Doc returned from NCBI was empty."
            print doc
            continue
        
        for t in ts.childNodes:
            if t.nodeName == "TSeq":
                gi = taxId = None
                for cn in t.childNodes:
                    if cn.nodeName == "TSeq_taxid":
                        #print cn.nodeName
                        #print cn.childNodes[0].nodeValue
                        taxId = int(cn.childNodes[0].nodeValue)
                    elif cn.nodeName == "TSeq_gi":
                        gi = int(cn.childNodes[0].nodeValue)

                if gi and taxId:
                    d[gi] = taxId
                else:
                    print "Error reading.  gi=%s and taxId=%s" % (gi, taxId)
        
    return d


def getFastaFromGi (ids):
    """Returns a dictionary where the key is the id and the value is the
    Taxonomy id for each id in the list that can be found
    using the NCBI EUTILS."""
    d = {}
    idLists = batchList(ids, batchSize=BATCH_SIZE)
    #print idLists
    for idList in idLists:
        ids = _getThrottledFasta(idList)
        #print len(ids)
        d.update(ids)

    return d

def _getThrottledFasta(ids):
    args = {}
    args["db"] = "nuccore"
    args["retmode"] = "xml"
    args["id"] = ",".join([str(x) for x in ids])
    args["rettype"] = "fasta"

    sioBuffer = efetch(EFETCH_URL, args)
    #print sioBuffer.read()
    #sioBuffer.seek(0)

    docs = []
    for xml in parseMultiXML(sioBuffer):
        try:
            doc = minidom.parse(StringIO.StringIO(xml))
        except:
            print "Unable to parse:"
            print xml
        #print doc
        #doc.seek(0)
        docs.append(doc)

    d = {}
    for doc in docs:
        try:
            ts = doc.getElementsByTagName("TSeqSet")[0]
        except IndexError:
            print "XML Doc returned from NCBI was empty."
            print doc
            continue
        
        for t in ts.childNodes:
            if t.nodeName == "TSeq":
                accver = seq = type = None
                for cn in t.childNodes:
                    if cn.nodeName == "TSeq_seqtype":
                        type = cn.getAttribute("value")
                    elif cn.nodeName == "TSeq_sequence":
                        #print cn.nodeName
                        #print cn.childNodes[0].nodeValue
                        seq = cn.childNodes[0].nodeValue
                    elif cn.nodeName == "TSeq_accver":
                        accver = cn.childNodes[0].nodeValue

                if accver and seq and type=="nucleotide":
                    d[accver] = seq
                else:
                    #print "Error reading.  accver=%s and type=%s and seq=%s" % (accver, type, seq)
                    if type=="nucleotide":
                        print "Error reading.  accver=%s and type=%s and seq=%s" % (accver, type, seq)
    return d
