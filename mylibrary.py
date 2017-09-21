def toHash(aList):
    #print aList
    #print "YO"
    return aList[1]+aList[8]

def overlap(a1, a2, b1, b2):
    result = False

    if  (a1<=b2 and b2<=a2) or (a1<=b1 and b1<=a2):
        result = True

    return result

def calculate_overlap(a1, a2, b1, b2):
    if (a1 < b1 and a2 <= b2 ):
        return a2 - b1
    elif (a1 < b1 and b2 <= a2):
        return b2 - b1
    elif (b1 <= a1 and a2 <= b2):
        return a2 - a1
    elif (b1 <= a1 and b2 <= a2):
        return b2 - a1
    #return 0

