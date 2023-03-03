C**************************************************************************
C
C  CALNORM.F - Program to calculate the normalization factors for the
C              dihedral angle distributions in prodata.
C
C              Written by Roman Laskowski, University College, London,
C              May 1996.
C
C              Modified 20 May 1999. (RAL)
C
C--------------------------------------------------------------------------
C
C Compilation and linking (on unix)
C -----------------------
C
C f77 -u -c calnorm.f
C f77 -u -c ps.f
C f77 -o calnorm calnorm.o ps.o
C
C--------------------------------------------------------------------------
C
C     Files
C     -----
C
C  2  prodata        - File holding the torsion angle distributions for
C                      phi-psi and chi1-chi2 combinations. Data generated
C                      by program gentors from a given database of high
C                      resolution protein structures. (Latest version is
C                      163 non-homologous structures from the Apr 93
C                      release of the Brookhaven databank)
C  8  prodata.new    - Output file being a copy of the input prodata file
C                      with all the new normalization factors slooted in
C
C--------------------------------------------------------------------------
C
C     Subroutine calling tree
C     -----------------------
C
C     MAIN    --> INITS   --> RAN2
C
C                 Read in and uncompress the torsion angle distributions
C             --> GETDAT
C             --> SMEAR
C             --> CALCLO
C             --> FINORM  --> RAN2
C             --> UPDATE
C
C----------------------------------------------------------------------+---


      PROGRAM CALNOR

      INCLUDE 'calnorm.inc'

      INTEGER       DISTRB, IRAN, NCELL1, NCELL2
      LOGICAL       NONORM
      REAL          ENERGY(MXCELL*(NAMINO+1)), NOBSER(MXCELL*(NAMINO+1))

C---- Initialise variables
      NONORM = .TRUE.
      CALL INITS(IRAN)
      IF (IFAIL) GO TO 999

C---- Loop for the different types of distributions to be plotted
      DO 100, DISTRB = 1, NDISTR

C----     Read in the torsion angle distributions and uncompress the
C         data
          CALL GETDAT(DISTRB,TWODEE,NOBSER,MXCELL,NCELL1,NCELL2,NAMINO,
     -        VALBEG,VALEND,STEP,NCOUNT,NRMEAN(1,DISTRB),
     -        NRMSTD(1,DISTRB),2,NONORM,IFAIL)
          IF (IFAIL) GO TO 999

C----     Smear the observations over the 2D distributions
          IF (TWODEE) THEN
              CALL SMEAR(NOBSER,ENERGY,NCELL1,NCELL2,NCOUNT,NAMINO)
          ENDIF

C----     Calculate log-odds scores
          CALL CALCLO(NOBSER,ENERGY,NCELL1,NCELL2,NAMINO,NCOUNT)

C----     Calculate normalization factors
          CALL FINORM(IRAN,DISTRB,TWODEE,NOBSER,ENERGY,NCELL1,NCELL2,
     -        NAMINO,NCOUNT,NRMEAN,NRMSTD,NDISTR)

 100  CONTINUE

C---- Copy the prodata file, slotting in the new normalization factors
      CALL UPDATE(NRMEAN,NRMSTD,NDISTR,NAMINO,IFAIL)

 999  CONTINUE
      IF (IFAIL) THEN
          PRINT*, '*** Program calnorm terminated with error'
      ELSE
         PRINT*, '* Program complete'
      ENDIF
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE INITS  -  Initialise variables
C
C----------------------------------------------------------------------+---

      SUBROUTINE INITS(IRAN)

      INCLUDE 'calnorm.inc'

      INTEGER       IRAN
      REAL          RAN2, X

C---- Initialise variables

C---- Initialise the random-number generator
      IRAN = -1
      X = RAN2(IRAN)

C---- Open data file, prodata
      OPEN(UNIT=2, FILE='prodata', STATUS='OLD',
     -     FORM='FORMATTED', ACCESS='SEQUENTIAL',
CVAX     -     CARRIAGECONTROL = 'LIST', READONLY,
     -     ERR=900)

      GO TO 999

900   CONTINUE
      PRINT*, '*** ERROR. Unable to open parameter file, prodata'
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.

999   CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE CALCPE  -  Convert the raw data-counts into energy values
C                        using the inverse Boltzmann distribution method
C                        of Sippl.
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE CALCPE(DISTRB,TWODEE,NOBSER,ENERGY,NCELL,NCELL1,NAMINO,
     -    NCOUNT,SIGMA,RT)

      INTEGER       NCELL, NCELL1, NAMINO

      INTEGER       DISTRB, IAMINO, ICELL, JCELL, NCOUNT(NAMINO + 1)
      LOGICAL       TWODEE
      REAL          ENERGY(NCELL,NCELL,NAMINO + 1), F, G, GRDMAX,
     -              NOBSER(NCELL,NCELL,NAMINO + 1), RT, SIGMA, TERM1,
     -              TERM2

C---- Loop over all the residue types for this distribution
      DO 600, IAMINO = 1, NAMINO

C----     Process only if have data points for this residue's distribution
          IF (NCOUNT(IAMINO).GT.0) THEN

C----         Initialise the maximum energy for this residue-type
              GRDMAX = 0.0

C----         Loop through all the cells for the current residue
              DO 200, JCELL = 1, NCELL
                  DO 100, ICELL = 1, NCELL

C----                 If no data points at all for this cell, set
C                     energy to a high level
                      IF (NOBSER(ICELL,JCELL,NAMINO + 1).EQ.0.0) THEN
                          ENERGY(ICELL,JCELL,IAMINO) = 999.9

C----                 Calculate potential energy for this cell
                      ELSE
                          F = NOBSER(ICELL,JCELL,NAMINO + 1)
     -                        / REAL(NCOUNT(NAMINO + 1))
                          G = NOBSER(ICELL,JCELL,IAMINO)
     -                        / REAL(NCOUNT(IAMINO))
                          TERM1 = LOG(1.0 + NCOUNT(IAMINO) * SIGMA)
                          TERM2 = LOG(1.0 + NCOUNT(IAMINO) * SIGMA
     -                        * G / F)
                          ENERGY(ICELL,JCELL,IAMINO)
     -                        = RT * (TERM1 - TERM2)
                          GRDMAX
     -                        = MAX(ENERGY(ICELL,JCELL,IAMINO),GRDMAX)
                      ENDIF
 100              CONTINUE
 200          CONTINUE

C----         Loop through all the grid points setting all invalid data points
C             to the maximum energy value
              DO 400, JCELL = 1, NCELL
                  DO 300, ICELL = 1, NCELL
                      IF (ENERGY(ICELL,JCELL,IAMINO).GT.GRDMAX) THEN
                          ENERGY(ICELL,JCELL,IAMINO) = GRDMAX
                      ENDIF
 300              CONTINUE
 400          CONTINUE
          ENDIF
 600  CONTINUE

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE FINORM  -  Calculate the normalisation factors for each
C                        distribution by randomly generating a large number
C                        of 'trial' distribution.
C
C----------------------------------------------------------------------+--- 

      SUBROUTINE FINORM(IRAN,DISTRB,TWODEE,NOBSER,ENERGY,NCELL1,NCELL2,
     -    NAMINO,NCOUNT,NRMEAN,NRMSTD,NDISTR)

      INTEGER       NLOOPS, NRESID

CDEBUG
C      PARAMETER    (NLOOPS = 100, NRESID = 100)
      PARAMETER    (NLOOPS = 10, NRESID = 100)
CDEBUG

      INTEGER       NAMINO, NCELL1, NCELL2,  NDISTR

      INTEGER       DISTRB, IAMINO, ICELL, IOBS, IRAN, LOOP, NCELLS,
     -              NCOUNT(NAMINO + 1), NOBS, NPOINT, NSTORE, NVALUE
      LOGICAL       TWODEE
      REAL          ENERGY(NCELL1*NCELL2,NAMINO + 1), MEAN,
     -              NOBSER(NCELL1*NCELL2,NAMINO + 1),
     -              NRMEAN(NAMINO+1,NDISTR), NRMSTD(NAMINO+1,NDISTR),
     -              RAN2, SCORAV, SCORE, SCOSTD, STDEV, X

C---- Initialise variables
      IF (TWODEE) THEN
          NCELLS = NCELL1 * NCELL2
      ELSE
          NCELLS = NCELL1
      ENDIF

C---- Loop over all the residue types for this distribution
      PRINT*, 'DISTRIBUTION', DISTRB
      DO 1000, IAMINO = 1, NAMINO

C----     Process only if have data points for this residue's distribution
          IF (NCOUNT(IAMINO).GT.0) THEN

C----         Initialise counts for this residue type
              MEAN = 0.0
              NVALUE = 0
              STDEV = 0.0

C----         Loop over the appropriate number of times
              DO 800, LOOP = 1, NLOOPS

C----             Initialise the maximum and minimum scores for this
C                 residue-type
                  NSTORE = 0
                  SCORAV = 0.0
                  SCOSTD = 0.0

C----             Loop until have generated all the required points
                  DO 600, NPOINT = 1, NRESID

C----                 Generate a random observation
                      X = RAN2(IRAN)
                      IOBS = INT(X * REAL(NCOUNT(IAMINO)))
                      IF (IOBS.LT.1) IOBS = 1
                      IF (IOBS.GT.NCOUNT(IAMINO)) IOBS = NCOUNT(IAMINO)

C----                 Initialise obervations count
                      NOBS = 0

C----                 Loop over all the cells
                      DO 500, ICELL = 1, NCELLS

C----                     Increment observations count with value in this
C                         cell
                          NOBS = NOBS + NOBSER(ICELL,IAMINO)

C----                     If have hit the number we're after, this is
C                         the cell
                          IF (NOBS.GE.IOBS) THEN

C----                         Get the cell's "energy" value and add to
C                             cumulative mean and st. dev. scores
                              NSTORE = NSTORE + 1
                              SCORE = ENERGY(ICELL,IAMINO)
                              SCORAV = SCORAV + SCORE
                              SCOSTD = SCOSTD + SCORE * SCORE

C----                         Finish here
                              GO TO 510
                          ENDIF
 500                  CONTINUE

C----                 Found the cell we were after
 510                  CONTINUE
 600              CONTINUE

C----             Update accumulators for this residue type
                  MEAN = MEAN + SCORAV
                  STDEV = STDEV + SCOSTD
                  NVALUE = NVALUE + NSTORE

C----             Calculate the mean and standard deviation score
                  IF (NSTORE.GT.0) THEN
                      SCORAV = SCORAV / REAL(NSTORE)
                      SCOSTD = SQRT(SCOSTD / REAL(NSTORE) 
     -                     - SCORAV * SCORAV)
                  ENDIF
 800          CONTINUE

C----         Calculate the mean and standard deviation score for this
C             residue type
              IF (NVALUE.GT.0) THEN
                  MEAN = MEAN / REAL(NVALUE)
                  STDEV = SQRT(STDEV / REAL(NVALUE) - MEAN * MEAN)
              ENDIF
              PRINT*, IAMINO, '. Overall Mean and stdev:',
     -            MEAN, STDEV

C----         Store the mean and standard deviation values
              NRMEAN(IAMINO,DISTRB) = MEAN
              NRMSTD(IAMINO,DISTRB) = STDEV
          ENDIF
 1000 CONTINUE

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  FUNCTION RAN2  -  Returns a uniform random deviate between 0.0
C                    and 1.0. Set IDUM to any negative value to
C                    initialise or reinitialise the sequence.
C                    This is a fast random number generator, its
C                    principal limitation being that it returns one
C                    of only 714,025 possible values.
C                    [Taken from Numerical Recipes and modified]
C
C----------------------------------------------------------------------+--- 

      REAL FUNCTION RAN2(IDUM)

      SAVE

      INTEGER       IA, IC, M
      REAL          RM
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112E-6)

      INTEGER       IDUM, IFF, IR(97), IY, J
      REAL          SEED

      DATA IFF /0/

      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN

C----     Open random seed file and read in seed
          SEED = 0.0
          OPEN (UNIT=13,STATUS='UNKNOWN',FILE='RANSEED.DAT',ERR=20)
          READ (13,*,END=20,ERR=20) SEED
20        CONTINUE
          IDUM = SEED
          IF (IDUM.EQ.0) IDUM = 1
          CLOSE (13)

          IDUM=MOD(IC-IDUM,M)
          DO 11 J=1,97
              IDUM=MOD(IA*IDUM+IC,M)
              IR(J)=IDUM
11        CONTINUE
          IDUM=MOD(IA*IDUM+IC,M)
          IY=IDUM
      ENDIF
      J=1+(97*IY)/M
      IF(J.GT.97.OR.J.LT.1) STOP 'Random number error'
      IY=IR(J)
      RAN2=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM

C---- If this is the first call, use the generated random number to
C     create a new seed and write to disk
      IF (IFF.EQ.0) THEN
          IFF = 1
          SEED = IY * RM * 100000
          OPEN (UNIT=13,STATUS='UNKNOWN',FILE='RANSEED.DAT',ERR=120)
          WRITE (13,*,ERR=120) SEED
120       CONTINUE
          CLOSE (13)
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C     SUBROUTINE UPDATE  -   Update the prodata file to add in the
C                            normalisation scores
C
C----------------------------------------------------------------------+---

      SUBROUTINE UPDATE(NRMEAN,NRMSTD,NDISTR,NAMINO,IFAIL)

      INTEGER       NAMINO, NDISTR

      CHARACTER*80  FNAME, IREC
      INTEGER       IAMINO, IDISTR, LENSTR, LINE
      LOGICAL       ENDFIL, IFAIL, INEW, INORM, IPOS, WRNORM, WROUT
      REAL          NRMEAN(NAMINO+1,NDISTR), NRMSTD(NAMINO+1,NDISTR)

C---- Initialise variables
      IFAIL = .FALSE.
      LINE = 0

C---- Open output file, gentorsN.dat, where N is the distribution number
      FNAME = 'prodata.new'
      OPEN(UNIT=8,FILE=FNAME,STATUS='UNKNOWN',ERR=900,
CVAX     -    CARRIAGECONTROL='LIST',
     -    FORM='FORMATTED',ACCESS='SEQUENTIAL')

C---- Rewind the input prodata file
      REWIND(2)

C---- Initialise flags
      ENDFIL = .FALSE.
      INEW =  .TRUE.
      INORM = .FALSE.

C---- Read through the prodata file, copying out each line, and slotting
C     in the new normalization factors
 100  CONTINUE

C----     Read in the next line from the input file
          LINE = LINE + 1
          READ(2,120,ERR=902,END=800) IREC
 120      FORMAT(A)
          WROUT = .FALSE.
          WRNORM = .FALSE.

C----     If this is the start of a distribution, then get which one
          IF (IREC(1:13).EQ.'Distribution ') THEN
              READ(IREC,160,ERR=904) IDISTR
 160          FORMAT(13X,I1)

C----         Set flags
              INEW =  .TRUE.
              INORM = .FALSE.
              WROUT = .TRUE.

C----     If this is the line containing the normalization factors, then
C         ignore
          ELSE IF (IREC(1:4).EQ.'Norm') THEN

C----         Set flags
              INEW =  .TRUE.
              INORM = .FALSE.

C----     If expecting the distribution for the next amino acid, get
C         which one we have
          ELSE IF (INEW) THEN

C----         Get the amino acid code
              READ(IREC,440,ERR=904) IAMINO
 440          FORMAT(5X,I2)

C----         Check that have a valid amino acid
              IF (IAMINO.LE.0 .OR. IAMINO.GT.NAMINO) GO TO 906

C----         Set flags
              INEW =  .FALSE.
              WROUT = .TRUE.

C----     Otherwise, must have a details line, so write straight out
          ELSE
              WROUT = .TRUE.

C----         Check whether this is the last of the details lines (ie it
C             has an asterisk)
              IPOS = INDEX(IREC,'*')
              IF (IPOS.NE.0) THEN
                  INEW = .TRUE.
                  WRNORM = .TRUE.
              ENDIF
          ENDIF

C----     If record to be written out, then do so
          IF (WROUT) THEN

C----         Write the record
              WRITE(8,120) IREC(1:LENSTR(IREC))
          ENDIF

C----     If writing out the normalisation factors, then do so
          IF (WRNORM) THEN

C----         Write the normalisation factors
              WRITE(8,520) NRMEAN(IAMINO,IDISTR), NRMSTD(IAMINO,IDISTR)
 520          FORMAT('Normalisation: ',2F10.4)
          ENDIF

C---- Loop back for the next record
      GO TO 100

C---- End of file reached
 800  CONTINUE


      GO TO 999

C---- Fatal errors
 900  CONTINUE
      PRINT*, '*** Error opening ', FNAME(1:LENSTR(FNAME)),
     -   'output file. Program aborted.'
      GO TO 990

 902  CONTINUE
      PRINT*, '*** Error reading parameter file, prodata, at',
     -     ' line', LINE
      GO TO 990

904   CONTINUE
      PRINT*, '*** Error reading parameter file, prodata, at',
     -     ' line', LINE
      GO TO 990

 906  CONTINUE
      PRINT*, '*** ERROR. Invalid amino acid code in parameter file,',
     -    ' prodata, at line', LINE, '       :', IAMINO
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
