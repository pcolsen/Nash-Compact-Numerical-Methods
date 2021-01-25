# Nash-Compact-Numerical-Methods

This is a project to convert the Pacal Routines in John Nash's book "Compact Numerical Methods for Computers: Linear Algebra and FUnction Minimization" into Python,
as well as to bring up-to-date some of the existing implementations of the algorithms in other programming languages.

Our goal is to provide a good set of modest but usable numerical methods in Python source for people who can't use NumPy or SciPy because there's no libraries for their platform, their platform is too small, or because the work in an environment where binary libraries can't be imported. Note that these algorithms are generally not the ones chosen for mainstream use. Without some effort at tuning, they will be slower and possibly provide less accurate solutions than codes from more extensive libraries. On the other hand, they are relatively short, hopefully easier to follow, and almost certainly more adaptable to new computing environments.

We plan to do this in two stages.  First we plan to translate the code from Pascal into Python so that the interfaces conform as closesly as possible the original Pascal so that users can refer to Professor Nash's original text.  Second we plan to bring the interfases as close as possible to those of NumPy and SciPy.

Note: The original Pascal code (available as the pascal collection on netlib.org) is not ready to use.  It was written to run under TurboPascal version 5 running on MS-DOS. A very large proportion of the code is present to allow for the deficiencies of the operating system and to allow for general problems to be solved. In the present collection, we have stripped away all but essential code to illustrate the algorithms, which for Pascal are run under the cross-platform Free Pascal Compiler (https://www.freepascal.org/). 

In the process of examining the Pascal code, Professor Nash has started to add Fortran, BASIC and R implementations of some of the algorithms, as well as provide some documentation to explain the choices made in both the algorithms and their particular implementations. We hope this will be instructive for those using the algorithms in new environments. We welcome comments and collaborations.
