C----------------------------------------------------------------------+---
C  GENTORS.INC - Include variables for GENTORS.F
C----------------------------------------------------------------------+---

CVAX      IMPLICIT NONE

      INTEGER       ICUT, MXCELL, MXLOG, NAMINO, NDISTR
      PARAMETER    (
     -              ICUT    =    6,
     -              MXCELL  = 7200,
     -              MXLOG   = 2000,
     -              NAMINO  =   20,
     -              NDISTR  =    7
     -             )

      REAL          CUTOFF, EMIN, SIGMA
      PARAMETER    (
     -              CUTOFF  =  0.0,
     -              EMIN    =  1.0E-06,
     -              SIGMA   =  1.0 / 50.0
     -             )

C---- PostScript parameters
      INTEGER       BBOXX1, BBOXX2, BBOXY1, BBOXY2, NTICKX, NTICKY
      PARAMETER    (
     -              BBOXX1 =  30,
     -              BBOXX2 = 550,
     -              BBOXY1 =  50,
     -              BBOXY2 = 780,
     -              NTICKX =   4,
     -              NTICKY =   4
     -             )

      INTEGER       XLEFT3, XSEP3, XWID3, YBOTT3, YHIGH3, YSEP3, YTITGP
      PARAMETER    (
     -              XLEFT3 =  80,
     -              XSEP3  =  30,
     -              XWID3  = 130,
     -              YBOTT3 = 120,
     -              YHIGH3 = 130,
     -              YSEP3  =  70,
     -              YTITGP =  8
     -             )


      CHARACTER*1   CHAIN
      CHARACTER*3   AMINO(NAMINO)
      CHARACTER*4   BCODE, PDBEXT, PDBPRE, RINEXT
      CHARACTER*78  FILPS
      CHARACTER*200 PDBDIR

      COMMON/CHARS/ CHAIN,

     -              AMINO,

     -              BCODE, PDBEXT, PDBPRE, RINEXT,

     -              FILPS,

     -              PDBDIR

      INTEGER       DISTRB, MVALUE, NCOUNT(NAMINO+1), NFILES, NRESID,
     -              PLOPT, TRESID
      COMMON/INTS/  DISTRB, MVALUE, NCOUNT, NFILES, NRESID,
     -              PLOPT, TRESID

      LOGICAL       ENDLST, IFAIL, TWODEE
      COMMON/LOGS/
     -              ENDLST, IFAIL, TWODEE

      REAL          ADJFAC, ALPHA, AVEDIF, MAXACT, MAXDIF, MINACT,
     -              NEWMAX, NEWMIN, PERMAX, PERMIN, PI, RADDEG, RMSDIF,
     -              STEP(2), VALBEG(2), VALEND(2)
      COMMON/REALS/
     -              ADJFAC, ALPHA, AVEDIF, MAXACT, MAXDIF, MINACT,
     -              NEWMAX, NEWMIN, PERMAX, PERMIN, PI, RADDEG, RMSDIF,
     -              STEP, VALBEG, VALEND

C----------------------------------------------------------------------+---
C Parameters
C ----------
C
C ADJFAC - Adjustment factor to normalise each Fourier series approximation
C          to the original distribution
C ALPHA  - Angular step used in calculating the Fourier series
C          approximation
C AMINO(NAMINO) - Three-letter residue code
C ARRAY(MXCELL,NAMINO+1) - Residue-by-residue counts of numbers of
C          points falling in each of the defined cells in the plot. In
C          routine SIPPL, these counts are converted into energy
C          differences
C AVEDIF - Average difference between actual and model distributions
C BBOXX1, BBOXX2, BBOXY1, BBOXY2 - Limits of bounding box round picture on
C          PostScript page
C BCODE  - Brookhaven code, as extracted from the name of the .rin file
C CHAIN  - Chain identifier for current protein structure. If blank, then
C          all chains are to be included
C CUTOFF - Cell occupancy cut-off. Only cells containing more examples
C          in the total distribution than this number are considred when
C          calculating energy differences
C DISTRB - Distribution being generated (1 - NDISTR)
C EMIN   - Minimum value, taken to be zero
C ENDLST - Flag indicating whether the end of the files listed in gentors.lst
C          has been reached
C IARRAY(MXCELL) - Integer versions of values in ARRAY, converted
C          to integers for more compact storage on disk
C ICUT   - Cutoff parameter used in calculation of log odds distribution
C          for each amino acid
C IFAIL  - Failure flag causing program to abort due to some error
C          condition
C INDEX(MXCELL) - Pointers used when sorting cell occupancies in
C          decreasing order of size
C MXCELL - Maximum number of cells stored for each distribution
C MXLOG  - Maximum number of occupancy values for the cells - used for
C          plotting log of frequency of each occupancy vs the occupany
C NAMINO - Number of amino acid types
C NCELL1, NCELL2 - Dimensions of array holding distribution
C NCOUNT(NAMINO+1) - Count of residues of each residue type in the dataset
C NDISTR - Number of distributions to be generated
C NEWMAX, NEWMIN - Maximum and minimum grid values of calculated distrib.
C NFILES - Number of files read in from gentors.lst
C NPARAM - Maximum number of parameters for the Fourier series approximation
C NRESID - Number of residues in current structure
C NTICKX, NTICKY - Numbers of ticks on each little graph on the page
C PDBDIR - Directory holding all the PDB files
C PDBEXT - PDB file extension
C PDBPRE - PDB file prefix
C PERMAX, PERMIN - Maximum and minimum values on the plotted graphs
C PI     - Value of pi
C PLOPT  - Plot option: 1=Plot frequency distributions; 2=Plot energy
C          differences
C RADDEG - Conversion factor from radians to degrees
C RHS(1,N2PARM) - Right-hand sides of normal equations
C RINEXT - .rin file extension
C RMSDIF - Rms difference between the actual distribution and that
C          calculated from the Fourier series approximation
C SIGMA  - Value of weight given to each quanta of information when
C          calculating the Sippl energy differences
C SOLTN(N2PARM) - Solution vector holding solutions to normal equations
C STEP(2) - Size of cell in each direction
C TARRAY(MXCELL) - Total counts of points falling in each of the cells
C TRESID - Total number of residues in the data set
C TWODEE - Flag indicating whether the distribution is a 1D or 2D
C          distribution
C VALBEG(2), VALEND(2) - Start and end limits for the value pair selected
C XLEFT3 - Left margin for 3x3 plot
C XSEP3  - x-distance between graphs in 3x3 plot
C XWID3  - Width of graphs in 3x3 plot
C YBOTT3 - Bottom margin for 3x3 plot
C YHIGH3 - Height of each graph in 3x3 plot
C YSEP3  - y-separation between graphs in 3x3 plot
C YTITGP - y-separation to graph title
C
C---------------------------------------------------------------------------
