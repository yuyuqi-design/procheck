C----------------------------------------------------------------------+---
C CALNORM.INC - Include variables for CALNORM.FOR
C----------------------------------------------------------------------+---

CVAX      IMPLICIT NONE

      INTEGER       MCHANG, MCHLEN, MINOBS, MXCELL, MXCOLR, MXRES,
     -              NAMINO, NDISTR, NMAIN, NPLANE, NALP
      PARAMETER    (
     -              MCHANG =     20,
     -              MCHLEN =     11,
     -              MINOBS  =    20,
     -              MXCELL  =  7200,
     -              MXCOLR  =    50,
     -              MXRES   = 20000,
     -              NAMINO  =    20,
     -              NDISTR  =     7,
     -              NMAIN  =  MCHANG + MCHLEN,
     -              NPLANE =      9,
     -              NALP   =  MCHANG + MCHLEN + NPLANE
     -             )

      INTEGER       MXFILE, MXLEND
      PARAMETER    (
     -              MXFILE = 60,
     -              MXLEND =  8
     -             )

      INTEGER       MAXRNG
      PARAMETER    (MAXRNG = 200)

      REAL          MINSHD, RT, SIGMA, VALLOW, VALUPP
      PARAMETER    (
     -              MINSHD =  0.995,
     -              RT     =  0.582,
     -              SIGMA  =  1 / 25.0,
     -              VALLOW =  0.00,
     -              VALUPP = -3.00
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

      REAL          XRAMC1, XRAMC2, YRAMC1, YRAMC2
      PARAMETER    (
     -              XRAMC1 =   98.0,
     -              XRAMC2 =  502.0,
     -              YRAMC1 =  313.0,
     -              YRAMC2 =  709.0
     -             )

      REAL          CHARHT, HWIDX, HWIDY
      PARAMETER    (
     -              CHARHT = (YRAMC2 - YRAMC1) / 400.0,
     -              HWIDX  = (XRAMC2 - XRAMC1) / 200.0,
     -              HWIDY  = (YRAMC2 - YRAMC1) / 200.0
     -             )

      INTEGER       XLEFT1, XWID1, YBOTT1, YHIGH1, YSEP1
      PARAMETER    (
     -              XLEFT1 = 170,
     -              XWID1  = 250,
     -              YBOTT1 = 127,
     -              YHIGH1 = 250,
     -              YSEP1  =  50
     -             )

      INTEGER       XLEFT3, XSEP3, XWID3, YBOTT3, YHIGH3, YSEP3, YTITGP
      PARAMETER    (
     -              XLEFT3 =  74,
     -              XSEP3  =  20,
     -              XWID3  = 140,
     -              YBOTT3 = 150,
     -              YHIGH3 = 130,
     -              YSEP3  =  58,
     -              YTITGP =  10
     -             )

      REAL          MKSZ, TXTSZ
      PARAMETER    (
     -              MKSZ   =  3.0,
     -              TXTSZ  =  6.0
     -             )

      REAL          GSTEP
      PARAMETER    (
     -              GSTEP  =  0.1
     -             )

      CHARACTER*1   CHAIN
      CHARACTER*3   AMINO(NAMINO+1)
      CHARACTER*10  RESID(MXRES)
      CHARACTER*12  COLNAM(MXCOLR)
      CHARACTER*80  BRCODE, FILHTM, FILLAN, FILPS, FILRIN(MXFILE),
     -              FILRNG, FILSDH, FILSUM, PDBFIL, RAMHED
      CHARACTER*(MXLEND) FILID(MXFILE)

      COMMON/CHARS/ CHAIN,

     -              AMINO,

     -              RESID,

     -              COLNAM,

     -              BRCODE, FILHTM, FILLAN, FILPS, FILRIN,
     -              FILRNG, FILSDH, FILSUM, PDBFIL, RAMHED,
     -              FILID

      INTEGER       ACTNUM(MXFILE), BLEN, COLCHI(7), COLRAM(7),
     -              COLRGP(5), FFILE, IEND, IPLOT, ISTART, LABRES,
     -              LENID(MXFILE), MCHNUM(NMAIN), MODEL(MXRES),
     -              MODFRM(MAXRNG), MODNUM(MXFILE), MODTO(MAXRNG),
     -              MRANGE, NALLOW, NAMLEN, NCOLOR, NCORE,
     -	            NCOUNT(NAMINO+1), NENDCH, NFILE, NGENER, NGLY,
     -              NLABEL, NMODEL, NONPG, NOUTSI, NPOINT, NPRO, NRANGE,
     -              NRAMPL, NRESID, NUMAVE(3), PSLEN, RLEN,
     -              SAVRES(MXRES), TOPMOD

      COMMON/INTS/
     -              ACTNUM, BLEN, COLCHI, COLRAM,
     -              COLRGP, FFILE, IEND, IPLOT, ISTART, LABRES,
     -              LENID, MCHNUM, MODEL,
     -              MODFRM, MODNUM, MODTO,
     -              MRANGE, NALLOW, NAMLEN, NCOLOR, NCORE,
     -	            NCOUNT, NENDCH, NFILE, NGENER, NGLY,
     -              NLABEL, NMODEL, NONPG, NOUTSI, NPOINT, NPRO, NRANGE,
     -              NRAMPL, NRESID, NUMAVE, PSLEN, RLEN,
     -              SAVRES, TOPMOD

      LOGICAL       ALLRAM, COMBPS, CRONLY, DISALL, DOPLOT(NDISTR),
     -              ENSEMB, HAVRAN, IFAIL, INCLET, INCOLC, INCOLG,
     -              INCOLR, MWANT(MXFILE), NMR, NOLINE, NOSHAD, PLABEL,
     -              PLTRES(NAMINO + 1), PLOTRM, RBXFIL, RBXMOD,
     -	            RMSPLT, RSELEC, TWODEE, USENGH, WITHAN
      COMMON/LOGS/
     -              ALLRAM, COMBPS, CRONLY, DISALL, DOPLOT,
     -              ENSEMB, HAVRAN, IFAIL, INCLET, INCOLC, INCOLG,
     -              INCOLR, MWANT, NMR, NOLINE, NOSHAD, PLABEL,
     -              PLTRES, PLOTRM, RBXFIL, RBXMOD,
     -	            RMSPLT, RSELEC, TWODEE, USENGH, WITHAN

      REAL          ALLOWP, COREPC, COVAVE(3), GAUMEA, GAUSTD,
     -              GENERP, LIMCHI, LIMGP, MCHMEA(NMAIN),
     -              NRMEAN(NAMINO+1,NDISTR), NRMSTD(NAMINO+1,NDISTR),
     -              OUTSIP, PERMAX, PERMIN, RBXSIZ, RESOL,
     -		    RGB(3,MXCOLR), SCORE(NDISTR,MXRES), SCOVAL(3,MXRES),
     -              STEP(2), VALBEG(2), VALEND(2), VALUE(NDISTR,MXRES)
      COMMON/REALS/ 
     -              ALLOWP, COREPC, COVAVE, GAUMEA, GAUSTD,
     -              GENERP, LIMCHI, LIMGP, MCHMEA,
     -              NRMEAN, NRMSTD,
     -              OUTSIP, PERMAX, PERMIN, RBXSIZ, RESOL,
     -		    RGB, SCORE, SCOVAL,
     -              STEP, VALBEG, VALEND, VALUE

      SAVE /CHARS/, /INTS/, /LOGS/, /REALS/


C----------------------------------------------------------------------+---
C Parameters
C ----------
C
C ACTNUM(MXFILE) - Actual model number corresponding to every entry in
C          MODNUM(MXFILE)
C ALLOWP - Percentage of residues in the allowed regions of the Ramachandran
C          plot
C ALLRAM - Flag indicating whether all the Ramachandran plot are required
C AMINO(NAMINO) - Three-letter residue code
C BBOXX1, BBOXX2, BBOXY1, BBOXY2 - Limits of bounding box round picture on
C          PostScript page
C BLEN   - Length of subheading on all plots
C BRCODE - Name of protein filename, as printed at the top of each plot
C CHAIN  - Chain identifier for current protein structure. If blank, then
C          all chains are to be included
C CHARHT - Character height for plots
C COLCHI(5), COLRAM(7), COLRGP(5) - User-defined colours for the Chi1-Chi2,
C          Gly & Pro, and Ramachandran plots, respectively
C COLNAM(MXCOLR) - Names identifying the user-defined colours
C COMBPS - Flag indicating whether all pages of a given plot are to be
C          combined into a single PostScript file
C COREPC - Percentage of residues in the core regions of the Ramachandran
C          plot
C COVAVE(3) - Average covalent scores for main-chain bond lengths, angles,
C          and overall
C CRONLY - Flag indicating whether only core region of Ramachandran plot
C          to be highlighted
C DISALL - Flag indicating whether the given residue in the disallowed
C          region of the Ramachandran plot
C DOPLOT(2) - Flags indicating which of the plots are required and which
C          not
C ENERGY(MXCELL,NAMINO+1) - Boltzmann energies obtained from the
C          counts stored in NOBSER
C ENSEMB - Flag indicating whether program has been called for an
C          ensemble of structures
C FFILE  - First file in a NMR ensemble that contains a wanted model
C FILHTM - Name of output html file holding summary information in html
C          format
C FILID(MXFILE) - Identifying name of each file in the ensemble
C FILLAN - Name of input .lan file
C FILPS  - Name of output PostScript file
C FILRIN(MXFILE) - Name of input .rin file(s)
C FILRNG - Name of input ranges file
C FILSDH - Name of input .sdh file
C FILSUM - Name of output .sum summary file
C GAUMEAN, GAUSTD - Normalization factors for log-odds scores generated
C          by the Gaussian distributions used for bond lengths and angles
C GENERP - Percentage of residues in the generous regions of the Ramachandran
C          plot
C GSTEP  - Interval size for calculating log-odds on Gaussians for main-
C          chain bond lengths and bond angles
C HAVRAN - Flag indicating whether have a ranges file
C HWIDX, HWIDY - Half-widths of markers for points in Ramachandran plot
C IEND   - End character-position for raw filename
C IFAIL  - Failure flag causing program to abort due to some error
C          condition
C INCLET - Flag determining whether lettering of regions required on
C          Ramachandran plot
C INCOLC, INCOLG, INCOLR - Flags indicating whether PostScript files are
C          required for the Chi1-Chi2, Gly & Pro, and Ramachandran plots
C IPLOT  - Plot numbers for PostScript pages
C ISTART - First character-position for raw filename
C LABRES - Which regions on the Ramachandran plot are to be labelled
C          (0 = disall, 1 = generous, 2 = allowed, 3 = core)
C LENID(MXFILE) - Length of each identifying name, FILID, for each file
C          in the ensemble
C LIMCHI, LIMGP - Limiting score determining which residues are to be
C          highlighted on the Chi1-Chi2 and the Gly & Pro Ramachandran
C          plots, repectively
C MCHANG - Number of main-chain bond angle types
C MCHLEN - Number of main-chain bond length types
C MCHMEA(NMAIN) - Mean value for each type of main-chain bond length and
C          angle
C MCHNUM(NMAIN) - Number of values for each type of main-chain bond
C          length and angle
C MINOBS - Minimum number of observations in an array for calculation of
C          probabilities
C MINSHD - Minimum shading corresponding to a favourable G-factor
C MAXRNG - Maximum number of model- or residue-ranges entered in the
C          ranges file
C MKSZ   - Marker size on plots
C MODEL(MXRES) - Model-number corresponding to residue (for NMR ensembles)
C MODFRM(MAXRNG), MODTO(MAXRNG) - Model-number ranges read in from the
C          ranges file
C MODNUM(MXFILE) - Model numbers (as defined by number of .rin file) of
C          those members of the ensemble that have been selected
C MRANGE - Number of model-number ranges in the ranges file
C MWANT(MXFILE) - Flags indicating which models are required
C MXCELL - Maximum number of cells stored for each distribution
C MXCOLR - Maximum number of RGB colours that can be defined by the user
C MXFILE - Maximum number of files allowed in NMR ensemble
C MXLEND - Maximum length of ensemble file-id that can be stored
C MXRES  - Maximum number of residues that can be stored
C NALLOW - Number of residues in the allowed regions of the Ramachandran plot
C NALP   - Total number of angle/length/plane types held in the .lan file
C NAMINO - Number of amino acid types
C NAMLEN - Length of the PDB code plus chain ID (if relevant)
C NCELL1, NCELL2 - Dimensions of array holding distribution
C NCORE  - Number of residues in the core regions of the Ramachandran plot
C NCOUNT(NAMINO+1) - Count of residues of each residue type in the dataset
C NCTOR1, NCTOR2 - Number of cells in the x- and y-directions for the
C          2D torsion angle distributions
C NDISTR   - Number of options for which distributions have been generated
C NENDCH - Number of end residues having meaningless phi-psi values
C NFILE  - Number of input .rin files
C NGENER - Number of residues in the generous regions of the Ramachandran
C          plot
C NGLY   - Number of glycine residues in the protein
C NLABEL - Number of labelled residues on the distribution plots
C NMAIN  - Number of main-chain bond length and bond angle types
C NMODEL - Number of models wanted
C NMR    - Flag indicating that the ensemble is a set of NMR models coming
C          from a single file
C NOBSER(MXCELL,NAMINO+1) - Residue-by-residue counts of numbers of
C          points falling in each of the defined cells in the plot. In
C          routine CALCPE, these counts are converted into energy
C          differences, stored in array ENERGY
C NOLINE - Flag determing whether lines round regions of Ramachandran plot
C          required
C NONPG  - Number of non-Gly and non-Pro residues in the protein
C NOSHAD - Flag determing whether shading of regions in Ramachandran plot
C          required
C NOUTSI - Number of residues in the disallowed regions of the Ramachandran
C          plot
C NPLANE - Number of types of planar groups held in the .lan file
C NPOINT - Number of points plotted on the distribution plots
C NPRO   - Number of proline residues
C NRAMPL - Number of points plotted on the Ramachandran plot
C NRANGE - Number of residue-number ranges in the ranges file
C NRESID - Number of residues in current structure
C NRMEAN(NAMINO+1,NDISTR) - Normalization factor for each distribution,
C          being the mean protein-by-protein score
C NRMSTD(NAMINO+1,NDISTR) - Standard deviation of NRMEAN, used in
C          normalizing the calculated log-odds scores for each residue
C          and each distribution
C NTICKX, NTICKY - Numbers of ticks on each little graph on the page
C NUMAVE(3) - Number of values contributing to the mean
C OUTSIP - Percentage of residues in the disallowed regions of the
C          Ramachandran plot
C PDBFIL - Input name of original PDB file
C PERMAX, PERMIN - Maximum and minimum values on the plotted graphs
C PLABEL - Flag indicating whether plot filename to be printed on the plot
C          itself
C PLOTRM - Flag signifying whether Ramachandran plot is required
C PLTRES(NAMINO + 1) - Flags indicating which residues are required for
C          plotting
C PSLEN  - Length of name of PostScript file
C RAMHED - Subheading for Ramachandran plot
C RBXFIL - Flag indicating whether data points on Ramachandran plot are to
C          be filled in
C RBXMOD - Flag indicating whether data points on Ramachandran plot are to
C          contain the model number
C RBXSIZ - Relative size of data-points on Ramachandran plot (1.0 corresponds
C          to default size)
C RESFRM(MAXRNG), RESTO(MAXRNG) - Residue-number ranges read in from the
C          ranges file
C RESID(MXRES) - Chain, residue name and sequence id
C RESOL  - Resolution of structure
C 
C RGB(3,MXCOLR) - User-defined RGB colour definitions
C RLEN   - Length of heading, RAMHED
C RMSPLT - Flag indicating whether the ensemble Ramachandran plots are
C          to be split up so that each one is on a separate page.
C RSELEC - Flag indicating whether residue-range selection on
C RT     - Product of Boltzmann's const and temperature. Used when deriving
C          potentials of mean force from the observed atom distributions
C SAVRES(MXRES) - Saved serial number for each residue used (necessary
C          when processing only a single chain)
C SCORE(NDISTR,MXRES) - Residue-by-residue log-odds score for each option
C SCOVAL(3,MXRES) - Residue-by-residue covalent scores for main-chain
C          bond lengths, angles, and overall
C SIGMA  - Weight for a single quantum of information, used when deriving
C          potentials of mean force from the observed atom distributions
C STEP(2) - Size of each cell in each direction
C TOPMOD - Highest model-number in the ensemble
C TWODEE - Flag indicating whether the current distribution is a 2D one
C TXTSZ  - Size of text for disallowed residues
C USENGH - Flag indicating whether the Engh & Huber means or the protein's
C	   own mean values to be used in calculating the covalent quality
C	   scores
C VALBEG(2), VALEND(2) - Start and end limits for the value pair selected
C VALLOW, VALUPP - Lower and upper G-factors defining the colour range
C          (values outside this range get the same colour as at the
C          range boundary
C VALUE(2,MXRES) - Residue-by-residue torsion-angle pairs
C WITHAN - Flag indicating whether file-handles are required for
C          PostScript output filenames
C XRAMC1, XRAMC2, YRAMC1, YRAMC2 - Graph area for main Ramachandran plot
C XLEFT1, XLEFT3 - Left margin for 1x3 and 3x3 plots
C XSEP3  - x-distance between graphs for 3x3 plots
C XWID1, XWID3  - Width of graphs for 1x3 and 3x3 plots
C YBOTT1, YBOTT3 - Bottom margin for 1x3 and 3x3 plots
C YHIGH1, YHIGH3 - Height of each graph for 1x3 and 3x3 plots
C YSEP1, YSEP3  - y-separation between graphs for 1x3 and 3x3 plots
C YTITGP - y-separation to graph title
C
C---------------------------------------------------------------------------
