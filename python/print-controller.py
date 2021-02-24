from io import TextIOWrapper

class printController():

    '''Provides a way to control interactive output
    =============================================
    This class provides a way to uniformly control
    the amount of real-time output from executing
    functions.  The user can set a "print-level"
    in an object created by this class, then 
    put use the "print" method to cause output.
    The "print" method has two arguments.  The
    first is a "level" argument.  This allows
    the user to set the degree of output
    desired.  The second is a string which will
    be output.  The user can also direct the
    output to the terminal, to the file and the 
    terminal, or to the file alone. Users can 
    open different file objects to customize 
    behavior for different parts of a program.
    '''
     
    
    def __init__(self, outputFile=None, termOutput=True):

        '''Creates a printController object.
        ====================================
        Arguments and defaults are:
        Defaults are
        printLevel = 0.  Must be integer >=0
        printLevel = 0 means no output
        outputFile = None.
        outputFIle must be open, writable, file object or None
        termOutput = True. Controls output to terminal
        termOutput must be boolean.
        '''

        # I construct three self variables with dummy values
        # then use setter methods to set their true values.
        # That lets me perform type and value checks.
        self.termOutput = True
        self.setTermOutput(termOutput)
        self.outFile = None
        self.setOutputFile(outputFile)
        # Default is no output.
        # If you want output you must use setPrintLevel set self.printLevel > 0
        # We set this ourselves, so no need to type and value check.
        self.printLevel = 0 

    def setPrintLevel(self, newLevel):
        '''Sets print level,  Arg must be integer >=0'''
        if newLevel is not integer:
            raise TypeError("printController.setPrintLevel: Print level must be integer")
        # This isn't strictly necesssary, but allowing for
        # negative print levels would make things confusing.
        if newLevel < 0:
            raise ValueError("printController.setPrintLevel: "
                                 + "Print level must be >= 0")
        self.printLevel = newLevel
        return self.prtlvl

    def getPrintLeval(self):
        return self.printLevel

    def setTermOutput(self, trueOrFalse):
        '''Sets term parameter. argv must be boolean.'''
        if trueOrFalse is not bool:
            raise TypeError("printController.setTermOutput: "
                                + "Terminal control arg must be 'True' or 'False'")
        self.termOutput = trueOrFalse
        return self.termOutput

    def getTermOutput(self):
        return self.termOutput

    def setOutputFile(self, newFileObject):
        '''Sets internal output file to an open, writable file object or None'''
        
        if newFileObject is None:
            self.outFile = None
        else:
            # Argument must be an open, writable file object.
            if (not newFileObject is TextIOWrapper):
                raise TypeErrof("printController.setOutputFile: "
                                +"Argument of setOutputFile must be a file object.")
            # We only get here if newFileObject is really a file
            if newFileObject.closed or (not newFileObject.writable()):
                # ValueError may not be the best choice, but it's reasonable
                raise ValueError("printController.setOutputFile: "
                                + "Argument of setOuputFile must be open and writable.")
            # We get here only if newFileObject is open and writable
            self.outFile = newFileObject
        return self.outFile

    def getOutputFile(self):
        return self.outFile

    def output(self, prtlvl, text):
        '''Provides text to be printed at a specified printLevel.'''
        
        # prtlvl gives the verbosity level at which this is to be printed.
        # self.printLevel is the verbosity level youi want.
        # If the verbosity level you want is greater than or equal to the
        # the text's verbosity level, then the text is printed.
        # (This comment needs to be clearer.)
        if not prtlvl is integer:
            raise TypeErrof("printController.output: "
                            + "Print level must be an integer.").
        # We print only if self.printLevel > 0 because 0 turns off all printing
        # and only if the level we want is >= level of this text.
        if (self.printLevel > 0) and (self.printLevel >= prtlvl):
            if not text = string:
                raise TypeError("printController.output: "
                                    + "Argument to output must be a string.")
            if self.printLevel >= prtlvl:
                if self.termOutput:
                    print(text)
                if self.outFile:
                    self.outFIile.write(text)
        
        return text

   


    
