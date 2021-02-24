class nashConstants():
  
    '''This file contains constants used througout "Compact Numerical Methods"
     ======================================================================
    '''

    def __init__(self):
        self.big = 1.0E+35    # a very large number
        self.Maxconst = 25    # Maximum number of constants in data record
        self.Maxobs = 100     # Maximum number of observations in data record
        self.Maxparm  = 25    # Maximum number of parameters to adjust
        self.Maxvars = 10     # Maximum number of variables in data record
        self.acctol = 0.0001  # acceptable point tolerance for minimisation codes
        self.maxm = 20        # Maximum number or rows in a matrix
        self.maxn = 20        # Maximum number of columns in a matrix
        self.maxmn = 40       # maxn+maxm, the number of rows in a working array
        self.maxsym = 210     # maximum number of elements of a symmetric matrix
              # which need to be stored = maxm * (maxm + 1)/2 
        self.reltest = 10.0   # a relative size used to check equality of numbers.
              # Numbers x and y are considered equal if the
              # floating-point representation of reltest+x equals
              # that of reltest+y.
        self.stepredn = 0.2   # factor to reduce stepsize in line search
        self.yearwrit = 1990  # year in which file was written

        self.banner = ""    # program name and description
        self.confile = ""   # file for output of console image
        self.confname = ""  # a name for confile
        self.dfile = ""     # file for output of console image
        self.dfname = ""    # a name for confile
        self.infile = ""    # an input file (keyboard image)
        self.infname = ""   # a name for the input file
        self.con = ""       # the console file

  

    


#   str2 = string[2];
#   rmatrix = array[1..maxm, 1..maxn] of real; {a real matrix}
#   wmatrix = array[1..maxmn, 1..maxn] of real; {a working array, formed
#                   as one real matrix stacked on another}
#   smatvec = array[1..maxsym] of real; {a vector to store a symmetric matrix
#               as the row-wise expansion of its lower triangle}
#   rvector = array[1..maxm] of real;  {a real vector. We will use vectors
#               of m elements always. While this is NOT space efficient,
#               it simplifies program codes.}
#   cgmethodtype= (Fletcher_Reeves,Polak_Ribiere,Beale_Sorenson);
#     {three possible forms of the conjugate gradients updating formulae}
#   probdata = record
#           m     : integer; {number of observations}
#           nvar  : integer; {number of variables}
#           nconst: integer; {number of constants}
#           vconst: array[1..Maxconst] of real;
#           Ydata : array[1..Maxobs, 1..Maxvars] of real;
#           nlls  : boolean; {true if problem is nonlinear least squares}
#         end;
# {

#   The following variables allow us to keep a copy of all screen
#   information in a file for some of the codes.  Pascal requires a
#   variable (confile in this case) for the file itself.  The string
#   variable confname is used for the name of the file.  Similar variables
#   allow problem data to be read from the file dfile named dfname.
