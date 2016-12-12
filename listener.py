#====================================================================================
#  Listener Functions
#====================================================================================
# Code borrowed from EPredict code written by me in 2010.  mtd.

# This module is a slightly fancy log file utility.  At the beginning of the code, a listener
# is created with a specific verbosity level, between 'silent' and 'verbose'.  All subsequent code then
# posts messages without knowing the initial user's preferred verbosity and the listener only logs those
# messages that meet or exceed the user's preferred verbosity level.

# NOTE: I'm currently avoiding Python's AbstractBaseClass because it is
# new in Python 2.6 and I'm not sure I want to require Pythin 2.6. 
# BUT, that is the right object-oriented choice here (or an interface,
# but Python doesn't support interfaces).

# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)

import time 

class Listener:
    """This Listener class is the 'Silent Listener' (or 'Null Listener') instance 
    and provides the base class for the rest of the listener classes.
    
    Levels:
    Silent (this class) = Level 0
    Error Only = Level 1 -- reportError only
    Important Only ("normal") = Level 2 -- reportError, reportImportantInfo
    Talkative = Level 3 -- reportError, reportImportantInfo, reportWarning
    Verbose = Level 4 -- reportError, reportImportantInfo, reportWarning, reportInfo
    Debugging = Level 5 -- reportError, reportImportantInfo, reportWarning, reportInfo, reportUnimportantInfo
    """
    def __init__(self, outputFileName=None):
        if outputFileName:
            self.out = open(outputFileName, "w")
        else:
            self.out = None
    
    def printStr(self, msg, level=None):
        """This class is responsible for printing a string to the output location.  Although it
        is not used by this class, since this is the Silent Listener, it is used by all classes
        inheriting from this class.
        If you need to change from a command line interface and want to stop printing to the command
        line or an error text file, change this code.
        """

        levelIndicator = ""
        if level == 1:
            levelIndicator = "*****"
        elif level == 2:
            levelIndicator = "    *"
        elif level == 3:
            levelIndicator = "         *"
        elif level == 4:
            levelIndicator = "         +"
        elif level == 5:
            levelIndicator = "         -"

        fullmsg = "%s: %s%s\n" % (time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), levelIndicator, msg)
        
        if self.out:
            self.out.write(fullmsg)
            self.out.flush()
        else:
            print fullmsg
    
    def reportError(self, str):
        """To be used when E-Predict has found an error condition but can find a way
        around it -- i.e.  not a situation when an Exception must be thrown, but a
        situation when the user really, really should be informed that something 
        went wrong.
        """
        pass

    def reportImportantInfo(self, str):
        """To be used for info that the user will almost always want to see.  For example,
        this is used to print the number of genomes and oligos found in the E-Matrix and
        the list of arrays analized.
        """
        pass
    
    def reportWarning(self, str):
        """To be used when E-Predict encounters a possibly dangerous situation. 
        For example, if there are oligos on the array that are not found in the E-Matrix,
        they should be printed as warnings."""
        pass
    
    def reportInfo(self, str):
        """To be used for less vital information that the user may want to see.
        """
        pass
    
    def reportUnimportantInfo(self, str):
        """This option should only be used for debugging.  For example, this would print
        all the oligos on the array that are also found in the E-Matrix and it would print
        all the user weights and their values at each iteration of the algorithm."""
        pass
    
    
class Level1Listener(Listener):
    
    def reportError(self, str):
        self.printStr(str, 1)
        
class Level2Listener(Listener):
    
    def reportError(self, str):
        self.printStr(str, 1)
        
    def reportImportantInfo(self, str):
        self.printStr(str, 2)
        
class Level3Listener(Listener):
    
    def reportError(self, str):
        self.printStr(str, 1)
        
    def reportImportantInfo(self, str):
        self.printStr(str, 2)
        
    def reportWarning(self, str):
        self.printStr(str, 3)
    
class Level4Listener(Listener):
    
    def reportError(self, str):
        self.printStr(str, 1)
        
    def reportImportantInfo(self, str):
        self.printStr(str, 2)
        
    def reportWarning(self, str):
        self.printStr(str, 3)
        
    def reportInfo(self, str):
        self.printStr(str, 4)
    
class Level5Listener(Listener):
    
    def reportError(self, str):
        self.printStr(str, 1)
        
    def reportImportantInfo(self, str):
        self.printStr(str, 2)
        
    def reportWarning(self, str):
        self.printStr(str, 3)
        
    def reportInfo(self, str):
        self.printStr(str, 4)
        
    def reportUnimportantInfo(self, str):
        self.printStr(str, 5)
        
def getListener(verbosity, outputFilename=None):
    
    if verbosity == 0:
        return Listener(outputFilename)
    elif verbosity == 1:
        return Level1Listener(outputFilename)
    elif verbosity == 2:
        return Level2Listener(outputFilename)
    elif verbosity == 3:
        return Level3Listener(outputFilename)
    elif verbosity == 4:
        return Level4Listener(outputFilename)
    elif verbosity == 5:
        return Level5Listener(outputFilename)
    else:
        raise epredictErrors.BadInput("Invalid verbosity level '%s'.  Verbosity levels must be between 0 and 5." % (verbosity))
