C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
C
C     --------------------------------------------------------------
      SUBROUTINE SCANFN (IFILE,JFILE)                                     I00010
C                                                                         I00020
      IMPLICIT REAL*8          (V)                                     ! I00030
C                                                                         I00040
C     DRIVER FOR CONVOLVING INSTRUMENTAL SCANNING FUNCTION                I00050
C     WITH SPECTRUM                                                       I00060
C                                                                         I00070
      COMMON S(3850),R1(5000),N1(5000)                                    I00080
C                                                                         I00090
      character*8      XID,       HMOLID,      YID,SCANID        
      real*8               SECANT,       XALTZ 
C                                                                         I00110
      COMMON /HVERSN/  HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,
     *                HVROPR,HVRPST,HVRPLT,HVRTST,HVRUTL,HVRXMR
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       I00120
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND,   I00130
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    I00140
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   I00150
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       I00160
     *               DVSC,XDUM,V1SHFT                                     I00170
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 I00180
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            I00190
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      I00200
      COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,
     *                HWF3,DXF3,NX3,N3MAX
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I00230
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I00240
     *              NLTEFL,LNFIL4,LNGTH4                                  I00250
      COMMON /SCINF/ HWHM,JEMIT,JFN,SAMPLE,SCANID,NPTS,XF(6018)           I00260
      COMMON /FLFORM/ CFORM                                               I00270
      COMMON /RCTSV/ JDUM,SDUM,JFLG,RNJDM,NB,IPC,VLFT,VCNT,VRGT,
     *               WGTL,WGTR
C                                                                         I00280
      CHARACTER*12 BCD,HTRANS,HABSRB,HRADIA                               I00290
      CHARACTER*11 CFORM                                                  I00300
      CHARACTER*8 HSCNID(0:6)                                             I00310
      CHARACTER*8 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *            HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR
      CHARACTER SCNOUT*7,SCNINF*7,CTAPE*4                                 I00320
      LOGICAL OP                                                          I00330
C                                                                         I00340
      DIMENSION FILHDR(2),SUMR(4)                                         I00350
      DIMENSION HWJ(0:6),DXJ(0:6),NJ(0:6),NJMX(0:6),SMPLJ(0:6),           I00360
     *          XSCAL(0:6)                                                I00370
C                                                                         I00380
      EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),                 I00390
     *            (FSCDID(6),ISCHDR) , (FSCDID(12),XSCID),                I00400
     *            (FSCDID(13),XHWHM) , (FSCDID(14),IDABS),                I00410
     *            (FSCDID(16),LAYR1)                                      I00420
C                                                                         I00430
      DATA HTRANS / 'TRANSMISSION'/,HABSRB / ' ABSORPTION '/,             I00440
     *     HRADIA / ' RADIANCE   '/                                       I00450
      DATA SCNOUT / '       '/,SCNINF / 'SCNINTF'/,CTAPE / 'TAPE'/        I00460
C                                                                         I00470
      DATA HSCNID(0) / 'RECTANGL'/,HWJ(0) / 1.         /,                 I00480
     *     DXJ(0) / 0.0  /,NJ(0) / 0    /,NJMX(0) / 0    /,               I00490
     *     SMPLJ(0) / .5 /,XSCAL(0) / 0.          /                       I00500
      DATA HSCNID(1) / 'TRIANGLE'/,HWJ(1) / 2.         /,                 I00510
     *     DXJ(1) / 0.02 /,NJ(1) / 101  /,NJMX(1) / 251  /,               I00520
     *     SMPLJ(1) / 2. /,XSCAL(1) / 0.          /                       I00530
      DATA HSCNID(2) / 'GAUSS   '/,HWJ(2) / 4.         /,                 I00540
     *     DXJ(2) / 0.02 /,NJ(2) / 201  /,NJMX(2) / 251  /,               I00550
     *     SMPLJ(2) / 4. /,XSCAL(2) / 0.          /                       I00560
C                                                                         I00570
C     SINCSQ: 54.18 HALFWIDTHS CORRESPONDS TO 24 ZERO CROSSINGS           I00580
C             PI CORRESPONDS TO X=2.257609141                             I00590
C                                                                         I00600
      DATA HSCNID(3) / 'SINCSQ  '/,HWJ(3) / 54.1826    /,                 I00610
     *     DXJ(3) / 0.02 /,NJ(3) / 2710 /,NJMX(3) / 2760 /,               I00620
     *     SMPLJ(3) / 4. /,XSCAL(3) / 1.391557377 /                       I00630
C                                                                         I00640
C     SINC: 119.33 HALFWIDTHS CORRESPONDS TO 72 ZERO CROSSINGS            I00650
C           PI CORRESPONDS TO X=1.657400255                               I00660
C                                                                         I00670
      DATA HSCNID(4) / 'SINC    '/,HWJ(4) / 119.332818 /,                 I00680
     *     DXJ(4) / 0.02 /,NJ(4) / 5968 /,NJMX(4) / 6018 /,               I00690
     *     SMPLJ(4) / 4. /,XSCAL(4) / 1.89549425  /                       I00700
      DATA HSCNID(5) / 'VRCTCENT'/,HWJ(5) / 1.         /,                 I00680
     *     DXJ(5) / 0.0  /,NJ(5) / 0    /,NJMX(5) / 0    /,               I00690
     *     SMPLJ(5) / .5 /,XSCAL(5) / 0.          /                       I00700
      DATA HSCNID(6) / 'VRCTLEFT'/,HWJ(6) / 1.         /,                 I00680
     *     DXJ(6) / 0.0  /,NJ(6) / 0    /,NJMX(6) / 0    /,               I00690
     *     SMPLJ(6) / .5 /,XSCAL(6) / 0.          /                       I00700
C                                                                         I00710
C----------------------------------------------------------------------   I00720
C                                                                         I00730
C    ADDITIONAL SCANNING FUNCTIONS MAY READILY BE ADDED TO THOSE          I00740
C      CURRENTLY IMPLEMENTED IN THIS VERSION OF LBLRTM:                   I00750
C                                                                         I00760
C    A SHAPE SUBROUTINE FOR THE DESIRED FUNCTION MUST BE CREATED-         I00770
C     THIS SUBROUTINE PRECALCULATES THE FUNCTION FOR SUBSEQUENT           I00780
C      LOOKUP.  SEE FOR EXAMPLE SUBROUTINE SHAPEG FOR THE GAUSSIAN        I00790
C                                                                         I00800
C    THE SHAPE SUBROUTINE SETS UP THE SYMMETRIC FUNCTION IN ARRAY FG      I00810
C     AT EQUAL INCREMENTS OF THE HALFWIDTH, 'DXF'. THE VALUE OF 'DXF'     I00820
C     IS SET IN THIS SUBROUTINE BY THE VALUE OF 'DXJ(?)'                  I00830
C                                                                         I00840
C    A DATA CARD MUST BE CREATED FOR EACH SCANNING FUNCTION DEFINING      I00850
C    THE FOLLOWING QUANTITIES:                                            I00860
C                                                                         I00870
C    HWJ(?)   EXTENT OF THE FUNCTION (BOUND) FROM THE CENTER IN UNITS     I00880
C               OF HALFWIDTH                                              I00890
C                                                                         I00900
C    DXJ(?)   INCREMENT AT WHICH THE FUNCTION IS STORED IN UNITS          I00910
C               OF HALFWIDTH                                              I00920
C                                                                         I00930
C    NJ(?)    THE NUMBER OF POINTS FROM THE CENTER TO THE FUNCTION        I00940
C               BOUND                                                     I00950
C                                                                         I00960
C    NJMAX(?) SIZE OF THE ARRAY IN WHICH THE FUNCTION IS STORED           I00970
C               FUNCTION VALUES BETWEEN NJ AND NJMAX ARE ZERO             I00980
C                                                                         I00990
C    SMPL(?)  DEFAULT VALUE OF THE SAMPLING INCREMENT IN RECIPRICAL       I01000
C               HALFWIDTH UNITS: E.G. A VALUE OF FOUR MEANS THAT THE      I01010
C               OUTPUT SPACING, 'DV', IN WAVENUMBERS WILL BE 1/4 THE      I01020
C               HALFWIDTH VALUE IN WAVENUMBERS, 'HWHM'.                   I01030
C                                                                         I01040
C    XSCAL(?) REQUIRED FOR PERIODIC FUNTIONS. THE VALUE OF THE            I01050
C               FUNCTION ARGUMENT IN RADIANS FOR WHICH THE                I01060
C               FUNCTION VALUE IS 0.5, E.G.                               I01070
C                   SINX/X = 0.5 FOR X = 1.89549425, XSCAL(4)             I01080
C                                                                         I01090
C    CONSIDERATION MUST BE GIVEN TO THE ISSUE OF FUNCTION                 I01100
C      NORMALIZATION FOR FUNCTIONS THAT DO NOT HAVE RAPID                 I01110
C      CONVERGENCE TO ZERO (SINX/X)                                       I01120
C                                                                         I01130
C                                                               SAC       I01140
C                                                                         I01150
C----------------------------------------------------------------------   I01160
C                                                                         I01170
C
C     ASSIGN SCCS VERSION NUMBER TO MODULE 
C
      HVRPST = '$Revision$'
C
      PI = 2.*ASIN(1.)                                                    I01180
C                                                                         I01190
C  SET THE MAXIMIM NUMBER OF AVAILABLE FUNCTIONS:                         I01200
C                                                                         I01210
      NFNMAX = 6                                                          I01220
C                                                                         I01230
C  NLIMF IS ONE MORE THAN THE SIZE OF OUTPUT (CONVOLVED) ARRAY            I01240
C                                                                         I01250
      NLIMF = 2401                                                        I01260
      NREN = 0                                                            I01270
      NSHIFT = 32                                                         I01280
      IFLSAV = 0                                                          I01290
      IPRT = 1                                                            I01300
C                                                                         I01310
   10 CONTINUE                                                            I01320
      SUMOUT = 0.                                                         I01330
      SMIN = 999999.                                                      I01340
      SMAX = -99999.                                                      I01350
      DVOSAV = 0.                                                         I01360
      SUMR(1) = SUMOUT                                                    I01370
      SUMR(2) = SMIN                                                      I01380
      SUMR(3) = SMAX                                                      I01390
      SUMR(4) = DVOSAV                                                    I01400
C                                                                         I01410
      IEOFT = 1                                                           I01420
C                                                                         I01430
C     READ IN CONTROL PARAMETERS: SEE INSTRUCTIONS FO DEFINITIONS         I01440
C                                                                         I01450
      READ (IRD,900,END=80) HWHM,V1,V2,JEMIT,JFN,JVAR,SAMPL,IUNIT,        I01460
     *                      IFILST,NIFILS,JUNIT,NPTS                      I01470
C                                                                         I01480
      IF (HWHM.LE.0.) GO TO 70                                            I01490
C                                                                         I01500
C     JEMIT=-1   SCANFN CONVOLVED WITH ABSORPTION                         I01510
C     JEMIT=0    SCANFN CONVOLVED WITH TRANSMISSION                       I01520
C     JEMIT=1    SCANFN CONVOLVED WITH EMISSION                           I01530
C                                                                         I01540
      JABS = 0                                                            I01550
      IDABS = 0                                                           I01560
      IF (JEMIT.LT.0) THEN                                                I01570
         JABS = 1                                                         I01580
         JEMIT = 0                                                        I01590
         IDABS = -1                                                       I01600
      ENDIF                                                               I01610
      IDABST = IDABS                                                      I01620
C                                                                         I01630
C     JVAR=1 FOR A VARIABLE SLIT FUNCTION (NOT FOR JFN=0)                 I01640
C     THE CODING IN CNVSCN  RESULTS IN HWHM=1./ (VI-V1)**2                I01650
C     HWHM IS CONSTANT FOR EACH PANEL AS PROGRAMMED                       I01660
C                                                                         I01710
      IFN = ABS(JFN)                                                      I01720
      IF (IFN.GT.NFNMAX) THEN                                             I01730
         WRITE (IPR,*)    'SCANF; JFN GT LIMIT'
         STOP             'SCANF; JFN GT LIMIT'
      ENDIF
C                                                                         I01740
      READ (HSCNID(IFN),905) SCANID                                       I01750
C                                                                         I01760
      HWF = HWJ(IFN)                                                      I01770
      DXF = DXJ(IFN)                                                      I01780
      NF = NJ(IFN)                                                        I01790
      NFMAX = NJMX(IFN)                                                   I01800
      SAMPLE = SMPLJ(IFN)                                                 I01810
      XSCALE = XSCAL(IFN)                                                 I01820
C                                                                         I01830
C     CHECK FOR NEGATIVE JFN OR NEGATIVE SAMPL                            I01840
C                                                                         I01850
C     FOR NEGATIVE JFN, USER IS SUPPLYING FIRST ZERO CROSSING FOR THE     I01860
C     PERIODIC FUNCTION IN HWHM.  SET HWHM=(FIRST ZERO)/(PI/XSCALE)       I01870
C
C     For JFN=5,6 user is supplying instrument field of view half angle
C     in degrees in HWHM.  Trap if JFN=-5,-6.
C                                                                         I01910
      IF (JFN.LT.0) THEN                                                  I01920
         JFN = ABS(JFN)                                                   I01930
         IF ((JFN.EQ.3).OR.(JFN.EQ.4)) THEN                               I01940
            HWHM = HWHM/(PI/XSCALE)                                       I01950
         ELSE                                                             I01960
            WRITE (IPR,910) JFN                                           I01970
            STOP 'SCANFN; INVALID JFN'                                    I01980
         ENDIF                                                            I01990
      ENDIF                                                               I02000
C
C     SET DVINT TO DETERMINE IF INTERPOLATION IS NECESSARY
C     - For JFN = 5,6, set DVINT to 1/12 the width of the first box.
C       HWHM should carry the value of the field of view half angle 
C       (in degrees).  This is converted to radians.  The box width
C       formula is
C
C                  width = V1*(1/2 angle FOV)**2/2
C     
C       and the degrees-to-radians formula is
C
C                  rad = deg*3.141592654/180.
C
C     - For JFN not equal to 5 or 6, set DVINT to 1/12 the value of
C       HWHM.  HWHM should carry the true value of the Half Width
C       at Half Maximum of the scanning function at this point.
C
      IF ((JFN.EQ.5).OR.(JFN.EQ.6)) THEN
         DVINT = V1*(HWHM*3.141592654/180.)**2/24
      ELSE
         DVINT = HWHM/12.                                                 I02010
      ENDIF
C
C     - For positive SAMPL, set SAMPLE equal to SAMPL (the number
C       of points per half width).
C     - For negative SAMPL, user is supplying desired DELVO
C       (outgoing spectral spacing).  SAMPLE (the number of sample
C       points per half width) is set such that SAMPLE=HWHM/DELVO
C       (Half Width at Half Max over user input outgoing spectral
C       spacing), and the outgoing spectral spacing DVO will be
C       recalculated using HWHM and SAMPLE below.
C
      IF (SAMPL.LT.0.) THEN
         SAMPLE = HWHM/(-SAMPL)                                           I02020
      ELSEIF (SAMPL.GT.0.) THEN
         SAMPLE = SAMPL                                                   I02030
      ENDIF
C                                                                         I02040
C     SET UP SELECTED SCANNING FUNCTION:                                  I02050
C                                                                         I02060
      IF (JFN.EQ.1) CALL SHAPET (XF)                                      I02070
      IF (JFN.EQ.2) CALL SHAPEG (XF)                                      I02080
      IF (JFN.EQ.3) CALL SINCSQ (XF,XSCALE)                               I02090
      IF (JFN.EQ.4) CALL SINC (XF,XSCALE)                                 I02100
C                                                                         I02110
      IF (IUNIT.LE.0) IUNIT = IFILE                                       I02120
      IFILE = IUNIT                                                       I02130
      IFILST = MAX(IFILST,1)                                              I02140
      IF (NIFILS.LE.0) NIFILS = 99                                        I02150
C                                                                         I02160
C     SKIP TO SELECTED 'FILE'                                             I02170
C                                                                         I02180
      REWIND IFILE                                                        I02190
      IF (IFILST.GT.1) CALL SKIPFL (IFILST-1,IFILE,IEOF)                  I02200
C                                                                         I02210
C     READ FILE HEADER FOR SELECTED 'FILE'                                I02220
C                                                                         I02230
   20 CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)                            I02240
      IF (IEOF.EQ.0) GO TO 10                                             I02250
      IDABS = IDABST                                                      I02260
C                                                                         I02270
      WRITE (IPR,915) XID,(YID(M),M=1,2)                                  I02280
      WRITE (IPR,920) LAYR1,LAYER                                         I02290
      WRITE (IPR,925) SECANT,PAVE,TAVE,DV,V1C,V2C                         I02300
      WRITE (IPR,930) WBROAD,(HMOLID(M),WK(M),M=1,NMOL)                   I02310
C                                                                         I02320
C     CHECK FOR INTERPOLATION AND OPEN OUTPUT FILE IF NECESSARY           I02330
C                                                                         I02340
C     IFILE INTERPOLATED ONTO JFILE                                       I02350
C                                                                         I02360
      IF (JUNIT.LE.0) JUNIT = JFILE                                       I02370
      JFILE = JUNIT                                                       I02380
C                                                                         I02390
C     IF DV NOT FINE ENOUGH, FIRST INTERPOLATE                            I02400
C                                                                         I02410
      IF (DV.GT.DVINT) THEN                                               I02420
         IFLSAV = IFILE                                                   I02430
         JFLSAV = JFILE                                                   I02440
         IEOFSC = 1                                                       I02450
         JFILE = 77                                                       I02460
         INQUIRE (UNIT=JFILE,OPENED=OP)                                   I02470
         IF (OP) CLOSE (JFILE)                                            I02480
         SCNOUT = SCNINF                                                  I02490
         OPEN (JFILE,FILE=SCNOUT,STATUS='UNKNOWN',FORM=CFORM)             I02500
         REWIND JFILE                                                     I02510
         IBUF = 0                                                         I02520
C                                                                         I02530
C     INTERPOLATE:                                                        I02540
C                                                                         I02550
         CALL SCNINT (IFILE,JFILE,DVINT,JEMIT,NPTS,IBUF)                  I02560
C                                                                         I02570
         IEOFSV = IEOFSC                                                  I02580
         WRITE (IPR,935)                                                  I02590
         IFILE = JFILE                                                    I02600
         REWIND IFILE                                                     I02610
         CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)                         I02620
         JFILE = JFLSAV                                                   I02630
      ELSE                                                                I02640
         IFLSAV = 0                                                       I02650
      ENDIF                                                               I02660
      INQUIRE (UNIT=JFILE,OPENED=OP)                                      I02670
      IF (.NOT.OP) THEN                                                   I02680
         WRITE (SCNOUT,940) CTAPE,JFILE                                   I02690
         OPEN (JFILE,FILE=SCNOUT,STATUS='UNKNOWN',FORM=CFORM)             I02700
         REWIND JFILE                                                     I02710
      ENDIF                                                               I02720
C                                                                         I02730
      ISCAN = ISCHDR                                                      I02740
      IF (ISCAN.LE.0.OR.XSCID.EQ.-99.) ISCAN = 0                          I02750
      IF (ISCHDR.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCHDR                   I02760
      ISCHDR = ISCAN+1                                                    I02770
      JTREM = -1                                                          I02780
      IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0                        I02790
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2                        I02800
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1                        I02810
      ISCANT = MOD(ISCAN,1000)                                            I02820
      IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 2                       I02830
      IF (JTREM.LT.0) THEN                                                I02840
         WRITE(IPR,*)   ' SCANF; JTREM LT 0'                              I02840
         STOP           ' SCANF; JTREM LT 0'                                    
      ENDIF                                                                     
      WRITE (IPR,945) SCANID,IFILE,IFILST,NIFILS,JEMIT,JFN,JVAR,JABS      I02850
C                                                                         I02860
C     JTREM=0   SCANFN CONVOLVED WITH EXPONENTIATED                       I02870
C                      ABSORPTION COEFFICIENT                             I02880
C     JTREM=1   SCANFN CONVOLVED WITH EMISSION                            I02890
C     JTREM=2   SCANFN CONVOLVED WITH TRANSMISSION                        I02900
C                                                                         I02910
      DVI = DV                                                            I02920
      DVSAV = DVI                                                         I02930
C
C     Compute output spectral spacing.  For JFN not 5 or 6, at this
C     point HWHM always contains the value of the Half Width at Half
C     Maximum of the scanning function, and SAMPLE always contains
C     the number of points per half width of the scanning function.
C
C     For JFN = 5,6 at this point, HWHM contains the value of the
C     field of view half angle (in degrees), and SAMPLE contains
C     the ratio of the field of view half angle to the specified
C     output spectral spacing (the quotient of HWHM and SAMPLE
C     results in the circuitous calculation of the previously input
C     DVO).
C
      DVO = HWHM/SAMPLE                                                   I02940
      IF (JFN.EQ.0) THEN                                                  I02950
         IRATIO = DVO/DVI+0.5                                             I02960
         DVO = FLOAT(IRATIO)*DVI                                          I02970
         IF (IRATIO.LT.2) THEN                                            I02980
            WRITE (IPR,950)                                               I02990
            GO TO 10                                                      I03000
         ENDIF                                                            I03010
      ENDIF                                                               I03020
C                                                                         I03030
C     BOUND AT THIS POINT IS THE WAVENUMBER VALUE                         I03040
C     OF HALF THE SCANNING FUNCTION                                       I03050
C                                                                         I03060
      BOUND = HWF*HWHM                                                    I03070
      DV = DVO                                                            I03080
      V1C = V1                                                            I03090
      V2C = V2                                                            I03100
      SCIND = JVAR+10*(JFN+10*(JEMIT))                                    I03110
      XSCID = SCIND+0.01                                                  I03120
      XHWHM = HWHM                                                        I03130
      CALL BUFOUT (JFILE,FILHDR(1),NFHDRF)                                I03140
      WRITE (IPR,955) HWHM,BOUND,JFILE,V1,V2,DVO                          I03150
      NBOUND = (2.*HWF)*SAMPLE+0.01                                       I03160
C                                                                         I03170
C     NBOUND IS THE NUMBER OF SPECTRAL VALUES SPANNED                     I03180
C     BY THE FULL SCANNING FUNCTION                                       I03190
C                                                                         I03200
C     RESET BOUND BASED ON NBOUND                                         I03210
C                                                                         I03220
      BOUND = FLOAT(NBOUND)*DVO/2.                                        I03230
      MAXF = NLIMF+2*NBOUND+NSHIFT                                        I03240
C                                                                         I03250
      TIMRDF = 0.                                                         I03260
      TIMCNV = 0.                                                         I03270
      TIMPNL = 0.                                                         I03280
      IEOFSC = 1                                                          I03290
      NLO = NSHIFT+1                                                      I03300
      SUMIN = 0.                                                          I03310
      NHI = NLIMF+NSHIFT-1                                                I03320
      DO 30 I = 1, MAXF                                                   I03330
         N1(I) = 0.                                                       I03340
         R1(I) = 0.                                                       I03350
   30 CONTINUE                                                            I03360
      INIT = 0                                                            I03370
      IDATA = -1                                                          I03380
      IPANEL = -1
      JFLG = -1
      VFT = V1-FLOAT(NSHIFT)*DV                                           I03390
      VBOT = V1-BOUND                                                     I03400
      VTOP = V2+BOUND                                                     I03410
C                                                                         I03420
      IF (JEMIT.EQ.0.AND.IDABS.EQ.0) BCD = HTRANS                         I03430
      IF (JEMIT.EQ.0.AND.IDABS.EQ.-1) BCD = HABSRB                        I03440
      IF (JEMIT.EQ.1) BCD = HRADIA                                        I03450
      IF (NPTS.GT.0) WRITE (IPR,960) BCD                                  I03460
C                                                                         I03470
   40 CALL CPUTIM (TIME0)                                                 I03480
C                                                                         I03490
      IF (IEOFSC.LE.0) GO TO 60                                           I03500
C                                                                         I03510
C     READ DATA TO BE CONVOLVE FROM IFILE  AND PUT INTO ARRAY S           I03520
C                                                                         I03530
      CALL RDSCAN (S,JTREM,IFILE,ISCAN,IPRT)                              I03540
C                                                                         I03550
CPRT  WRITE(IPR,965) IEOFSC,IDATA                                         I03560
C                                                                         I03570
      CALL CPUTIM (TIME)                                                  I03580
      TIMRDF = TIMRDF+TIME-TIME0                                          I03590
C                                                                         I03600
      IF (IEOFSC.LE.0) GO TO 60                                           I03610
C                                                                         I03620
C     SHRKSC MAY SHRINK (COMPRESS) THE DATA; DVI IS MODIFIED ACCORDINGL   I03630
C                                                                         I03640
      IF ((JFN.NE.0).AND.(JFN.NE.5).AND.(JFN.NE.6)) THEN
         CALL SHRKSC (INIT,HWHM)                                          I03650
      ENDIF
C                                                                         I03660
   50 CONTINUE                                                            I03670
C                                                                         I03680
C     PERFORM THE CONVOLUTION OF XF ON S TO GIVE R1                       I03690
C                                                                         I03700
      IF (JFN.EQ.0) THEN                                                  I03710
         CALL CNVRCT (S,HWHM,R1,XF)                                       I03720
      ELSEIF (JFN.EQ.5) THEN
         CALL CNVVRC (S,HWHM,R1,XF)                                       I03720
      ELSEIF (JFN.EQ.6) THEN
         CALL CNVVRL (S,HWHM,R1,XF)                                       I03720
      ELSE                                                                I03730
         CALL CONVSC (S,HWHM,R1,XF)                                       I03740
      ENDIF                                                               I03750
C                                                                         I03760
CPRT  WRITE(IPR,965) IEOFSC,IDATA,IPANEL                                  I03770
C                                                                         I03780
      IF (IPANEL.EQ.0) GO TO 40                                           I03790
C                                                                         I03800
   60 CONTINUE                                                            I03810
C                                                                         I03820
C     OUTPUT PANEL TO JFILE, NPTS VALUES OF R1                            I03830
C                                                                         I03840
      IF (JFN.EQ.0.OR.JFN.EQ.5.OR.JFN.EQ.6) THEN                          I03850
         CALL PNLRCT (R1,JFILE,SUMR,NPTS)                                 I03860
      ELSE                                                                I03870
         CALL PANLSC (R1,JFILE,SUMR,NPTS)                                 I03880
      ENDIF                                                               I03890
C                                                                         I03900
      IF ((ISTOP.NE.1).AND.(IEOFSC.LT.0)) GO TO 60                        I03910
      IF ((ISTOP.NE.1).AND.(IEOFSC.GT.0)) GO TO 50                        I03920
      CALL CPUTIM (TIME)                                                  I03930
      WRITE (IPR,970) TIME,TIMRDF,TIMCNV,TIMPNL                           I03940
      CALL ENDFIL (JFILE)                                                 I03950
C                                                                         I03960
      SUMIN = SUMIN*DVSAV                                                 I03970
C                                                                         I03980
      WRITE (IPR,975) SUMIN                                               I03990
C                                                                         I04000
      IF (IFLSAV.NE.0) THEN                                               I04010
         IFILE = IFLSAV                                                   I04020
         IEOFSC = IEOFSV                                                  I04030
      ENDIF                                                               I04040
      IF (IEOFSC.EQ.1) CALL SKIPFL (1,IFILE,IEOFSC)                       I04050
C                                                                         I04060
      IEOFT = IEOFT+1                                                     I04070
C                                                                         I04080
      SUMOUT = SUMR(1)                                                    I04090
      SMIN = SUMR(2)                                                      I04100
      SMAX = SUMR(3)                                                      I04110
      DVOSAV = SUMR(4)                                                    I04120
C                                                                         I04130
      SUMOUT = SUMOUT*DVOSAV                                              I04140
      WRITE (IPR,980) SUMOUT,SMIN,SMAX                                    I04150
C                                                                         I04160
      IF (IEOFT.LE.NIFILS.AND.IEOFSC.LT.0) GO TO 20                       I04170
C                                                                         I04180
      GO TO 10                                                            I04190
C                                                                         I04200
   70 CONTINUE                                                            I04210
C                                                                         I04220
   80 RETURN                                                              I04230
C                                                                         I04240
  900 FORMAT (3F10.3,3(3X,I2),F10.4,4(3X,I2),I5)                          I04250
  905 FORMAT (A8)                                                         I04260
  910 FORMAT (//,' *****  INVALID VALUE FOR JFN = ',I2,'  *****',/)       I04270
  915 FORMAT ('1',' **SCANFN** ',/,'0',10A8,2X,2(1X,A8,1X))               I04280
  920 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)            I04290
  925 FORMAT ('0 SECANT =',F15.5,/'0 PRESS(MB) =',F12.5/'0 TEMP(K) =',    I04300
     *        F11.2,/'0 DV(CM-1) = ',F12.8,/'0 V1(CM-1) = ',F12.6,/       I04310
     *        '0 V2(CM-1) = ',F12.6)                                      I04320
  930 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//5X,'WBROAD = ',       I04330
     *        1PE10.3,/(5X,A6,' = ',1PE10.3))                             I04340
  935 FORMAT ('0',' **SCANFN** ',/)                                       I04350
  940 FORMAT (A4,I2.2)                                                    I04360
  945 FORMAT ('0','***',A8,'***',//6X,'INPUT FILE NUMBER =',I3,           I04370
     *        ' ,IFILST = ',I5,' ,NIFILS = ',I5,',JEMIT =',I2,            I04380
     *        ' ,JFN =',I2,' ,JVAR =',I2,'  ,JABS =',I2)                  I04390
  950 FORMAT ('0',60X,'****** IRATIO LESS THAN 2, NO SCANFN ******')      I04400
  955 FORMAT (1X,'     HWHM OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,     I04410
     *        5X,'BOUND OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,6X,      I04420
     *        'OUTPUT FILE NUMBER =',I3,',   V1 =',F12.5,',   V2 =',      I04430
     *        F12.5,5X,' DV OUT',F12.8)                                   I04440
  960 FORMAT (///,'0',5X,A12,/)                                           I04450
  965 FORMAT ('0',5X,'IEOFSC =',I3,'  IDATA =',I3,'  IPANEL =',I3,/)      I04460
  970 FORMAT ('0',5X,'TIME =',F7.3,',  READ =',F6.3,',  CONV. =',F7.3,    I04470
     *        ',  PANEL =',F6.3)                                          I04480
  975 FORMAT ('0    SUMIN  =',1P,E16.9)                                   I04490
  980 FORMAT ('0    SUMOUT =',1P,E16.9,'  MIN =',E16.9,'  MAX =',E16.9)   I04500
C                                                                         I04510
      END                                                                 I04520
C
C     --------------------------------------------------------------
C
      SUBROUTINE SCANRD (DVINT,IEMIT)                                     I04530
C                                                                         I04540
      IMPLICIT REAL*8          (V)                                     ! I04550
C                                                                         I04560
C     READ CONTROL CARD FOR SCANNING WITH WEIGHTING FUNCTIONS             I04570
C                                                                         I04580
      COMMON S(3850),R1(5000)                                             I04590
C                                                                         I04600
      character*8      XID,       HMOLID,      YID,SCANID
      real*8               SECANT,       XALTZ 
C                                                                         I04620
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       I04630
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND,   I04640
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    I04650
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   I04660
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       I04670
     *               DVSC,XDUM,V1SHFT                                     I04680
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 I04690
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            I04700
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      I04710
      COMMON /SCSHAP/ HWFS,DXFS,NFS,NFMAXS
      COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,
     *                HWF3,DXF3,NX3,N3MAX
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I04740
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I04750
     *              NLTEFL,LNFIL4,LNGTH4                                  I04760
      COMMON /SCINF/ HWHM,JEMIT,JFN,SAMPLE,SCANID,NPTS,XF(6018)           I04770
      COMMON /FLFORM/ CFORM                                               I04780
C                                                                         I04790
      CHARACTER*8 HSCNID(0:6)                                             I04800
      CHARACTER CFORM*11,TAPE13*6,CTAPE*4                                 I04810
      LOGICAL OP                                                          I04820
C                                                                         I04830
      DIMENSION FILHDR(2)                                                 I04840
      DIMENSION HWJ(0:6),DXJ(0:6),NJ(0:6),NJMX(0:6),SMPLJ(0:6),           I04850
     *          XSCAL(0:6)                                                I04860
C                                                                         I04870
      EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(6),ISCHDR),                I04880
     *            (FSCDID(12),XSCID) , (FSCDID(13),XHWHM),                I04890
     *            (FSCDID(14),IDABS) , (FSCDID(16),LAYR1)                 I04900
C                                                                         I04910
      DATA HSCNID(0) / 'RECTANGL'/,HWJ(0) / 1.         /,                 I04920
     *     DXJ(0) / 0.0  /,NJ(0) / 0    /,NJMX(0) / 0    /,               I04930
     *     SMPLJ(0) / .5 /,XSCAL(0) / 0.          /                       I04940
      DATA HSCNID(1) / 'TRIANGLE'/,HWJ(1) / 2.         /,                 I04950
     *     DXJ(1) / 0.02 /,NJ(1) / 101  /,NJMX(1) / 251  /,               I04960
     *     SMPLJ(1) / 2. /,XSCAL(1) / 0.          /                       I04970
      DATA HSCNID(2) / 'GAUSS   '/,HWJ(2) / 4.         /,                 I04980
     *     DXJ(2) / 0.02 /,NJ(2) / 201  /,NJMX(2) / 251  /,               I04990
     *     SMPLJ(2) / 4. /,XSCAL(2) / 0.          /                       I05000
C                                                                         I05010
C     SINCSQ: 54.18 HALFWIDTHS CORRESPONDS TO 24 ZERO CROSSINGS           I05020
C             PI CORRESPONDS TO X=2.257609141                             I05030
C                                                                         I05040
      DATA HSCNID(3) / 'SINCSQ  '/,HWJ(3) / 54.1826    /,                 I05050
     *     DXJ(3) / 0.02 /,NJ(3) / 2710 /,NJMX(3) / 2760 /,               I05060
     *     SMPLJ(3) / 4. /,XSCAL(3) / 1.391557377 /                       I05070
C                                                                         I05080
C     SINC: 119.33 HALFWIDTHS CORRESPONDS TO 72 ZERO CROSSINGS            I05090
C           PI CORRESPONDS TO X=1.657400255                               I05100
C                                                                         I05110
      DATA HSCNID(4) / 'SINC    '/,HWJ(4) / 119.332818 /,                 I05120
     *     DXJ(4) / 0.02 /,NJ(4) / 5968 /,NJMX(4) / 6018 /,               I05130
     *     SMPLJ(4) / 4. /,XSCAL(4) / 1.89549425  /                       I05140
      DATA HSCNID(5) / 'VRCTCENT'/,HWJ(5) / 1.         /,                 I05120
     *     DXJ(5) / 0.0  /,NJ(5) / 0    /,NJMX(5) / 0    /,               I05130
     *     SMPLJ(5) / .5 /,XSCAL(5) / 0.          /                       I05140
      DATA HSCNID(6) / 'VRCTLEFT'/,HWJ(6) / 1.         /,                 I05120
     *     DXJ(6) / 0.0  /,NJ(6) / 0    /,NJMX(6) / 0    /,               I05130
     *     SMPLJ(6) / .5 /,XSCAL(6) / 0.          /                       I05140
C                                                                         I05150
      DATA TAPE13 / '      '/,CTAPE / 'TAPE'/                             I05160
C                                                                         I05170
      PI = 2.*ASIN(1.)                                                    I05180
C
C  SET THE MAXIMIM NUMBER OF AVAILABLE FUNCTIONS:
C
      NFNMAX = 6
C                                                                         I05190
      NLIMF = 2401                                                        I05200
      NSHIFT = 32                                                         I05210
      READ (IRD,900,END=10) HWHM,V1,V2,JEMIT,JFN,JVAR,SAMPL,NNFILE,NPTS   I05220
C                                                                         I05230
      IF (HWHM.LE.0.) THEN                                                I05240
         WRITE(IPR,*) ' SCANRD * HWHM NEGATIVE '                          I05242
         STOP         ' SCANRD * HWHM NEGATIVE '
      ENDIF                                                                     
C                                                                         I05250
C     JEMIT=-1   SCANFN CONVOLVED WITH ABSORPTION                         I05260
C     JEMIT=0    SCANFN CONVOLVED WITH TRANSMISSION                       I05270
C     JEMIT=1    SCANFN CONVOLVED WITH EMISSION                           I05280
C                                                                         I05290
      JABS = 0                                                            I05300
C                                                                         I05310
C     THE FOLLOWING CARDS HAVE BEEN RETRAINED                             I05320
C     FOR POSSIBLE FUTURE CODE ENHANCEMENTS                               I05330
C                                                                         I05340
CC    IF (JEMIT.LT.0) THEN                                                I05350
CC       JABS=1                                                           I05360
CC       JEMIT=0                                                          I05370
CC    ENDIF                                                               I05380
C                                                                         I05390
C     JVAR=1 FOR A VARIABLE SLIT FUNCTION (NOT FOR JFN=0)                 I05400
C     THE CODING IN CNVSCN  RESULTS IN HWHM=1./ (VI-V1)**2                I05410
C     HWHM IS CONSTANT FOR EACH PANEL AS PROGRAMMED                       I05420
C     FOLLOWING VALUES INITIALIZE FOR RECTANGLE                           I05430
C                                                                         I05440
      IFN = ABS(JFN)                                                      I05450
      IF (IFN.GT.NFNMAX) THEN                                             I05460
         WRITE(IPR,*)' SCANF; JFN GT LIMIT'                               I05462
         STOP        ' SCANF; JFN GT LIMIT'
      ENDIF
C                                                                         I05470
      READ (HSCNID(IFN),905) SCANID                                       I05480
C                                                                         I05490
C     JVAR=1 FOR A VARIABLE SLIT FUNCTION (NOT FOR JFN=0)                 I05500
C     THE CODING IN CNVSCN  RESULTS IN HWHM=1./ (VI-V1)**2                I05510
C     HWHM IS CONSTANT FOR EACH PANEL AS PROGRAMMED                       I05520
C     FOLLOWING VALUES INITIALIZE FOR RECTANGLE                           I05530
C                                                                         I05540
      HWF = HWJ(IFN)                                                      I05550
      DXF = DXJ(IFN)                                                      I05560
      NF = NJ(IFN)                                                        I05570
      NFMAX = NJMX(IFN)                                                   I05580
      SAMPLE = SMPLJ(IFN)                                                 I05590
      XSCALE = XSCAL(IFN)                                                 I05600
C
C     Set values of HWFS, DXFS, NFS, & NFMAXS to HWF, DXF, NF,
C     & NFMAX for use when entering HIRAC1 between SCANRD and
C     SCNMRG.
C
      HWFS = HWF
      DXFS = DXF
      NFS = NFS
      NFMAXS = NFMAX
C                                                                         I05610
C     CHECK FOR NEGATIVE JFN OR NEGATIVE SAMPL                            I05620
C                                                                         I05630
C     FOR NEGATIVE JFN, USER IS SUPPLYING FIRST ZERO CROSSING FOR THE     I05640
C     PERIODIC FUNCTION IN HWHM.  SET HWHM=(FIRST ZERO)/(PI/XSCALE)       I05650
C
C     For JFN=5,6 user is supplying instrument field of view half angle
C     in degrees in HWHM.
C                                                                         I05660
C     FOR NEGATIVE SAMPL, USER IS SUPPLYING DESIRED DELVO.                I05670
C     SET SAMPLE=HWHM/DELVO.                                              I05680
C                                                                         I05690
      IF (JFN.LT.0) THEN                                                  I05700
         JFN = ABS(JFN)                                                   I05710
         IF ((JFN.EQ.3).OR.(JFN.EQ.4)) THEN                               I05720
            HWHM = HWHM/(PI/XSCALE)                                       I05730
         ELSE                                                             I05740
            WRITE (IPR,910) JFN                                           I05750
            STOP 'SCANRD; INVALID JFN'                                    I05760
         ENDIF                                                            I05770
      ENDIF                                                               I05780
C                                                                         I05790
      IF (SAMPL.LT.0.) SAMPLE = HWHM/(-SAMPL)                             I05800
      IF (SAMPL.GT.0.) SAMPLE = SAMPL                                     I05810
C                                                                         I05820
      IF (JFN.EQ.1) CALL SHAPET (XF)                                      I05830
      IF (JFN.EQ.2) CALL SHAPEG (XF)                                      I05840
      IF (JFN.EQ.3) CALL SINCSQ (XF,XSCALE)                               I05850
      IF (JFN.EQ.4) CALL SINC (XF,XSCALE)                                 I05860
C                                                                         I05870
      IF (NNFILE.NE.NFILE.AND.NNFILE.GT.0) THEN                           I05880
         INQUIRE (UNIT=NFILE,OPENED=OP)                                   I05890
         IF (OP) CLOSE (NFILE)                                            I05900
         NFILE = NNFILE                                                   I05910
         INQUIRE (UNIT=NFILE,OPENED=OP)                                   I05920
         IF (.NOT.OP) THEN                                                I05930
            WRITE (TAPE13,915) CTAPE,NFILE                                I05940
            OPEN (NFILE,FILE=TAPE13,STATUS='UNKNOWN',FORM=CFORM)          I05950
            REWIND NFILE                                                  I05960
         ENDIF                                                            I05970
      ENDIF                                                               I05980
C                                                                         I05990
      WRITE (IPR,920) SCANID,JEMIT,JFN,JVAR,SAMPL,NPTS                    I06000
      IEMIT = JEMIT                                                       I06010
C                                                                         I06020
C    BOUND AT THIS POINT IS THE WAVENUMBER VALUE                          I06030
C    OF HALF THE SCANNING FUNCTION                                        I06040
C                                                                         I06050
      DVO = HWHM/SAMPLE                                                   I06060
      DVINT = HWHM/12.                                                    I06070
      BOUND = HWF*HWHM                                                    I06080
      V1C = V1                                                            I06090
      V2C = V2                                                            I06100
      XHWHM = HWHM                                                        I06110
      WRITE (IPR,925) HWHM,BOUND,NFILE,V1,V2                              I06120
      RETURN                                                              I06130
   10 PRINT 930                                                           I06140
      STOP                                                                I06150
C                                                                         I06160
  900 FORMAT (3F10.3,3(3X,I2),F10.4,15X,2I5)                              I06170
  905 FORMAT (A8)                                                         I06180
  910 FORMAT (//,' *****  INVALID VALUE FOR JFN = ',I2,'  *****',/)       I06190
  915 FORMAT (A4,I2.2)                                                    I06200
  920 FORMAT ('1',5X,'SCANRD',5X,'***',A8,'***',/,/,' JEMIT =',I2,        I06210
     *        ' JFN =',I2,' ,JVAR =',I2,' ,SAMPL =',F10.4,'  ,NPTS =',    I06220
     *        I5)                                                         I06230
  925 FORMAT (1X,'     HWHM OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,     I06240
     *        5X,'BOUND OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,6X,      I06250
     *        'OUTPUT FILE NUMBER =',I3,',   V1 =',F12.5,',   V2 =',      I06260
     *        F12.5)                                                      I06270
  930 FORMAT (' END OF FILE TAPE5',/,' (NOTE TAPE10 ALREADY CREATED )')   I06280
C                                                                         I06290
      END                                                                 I06300
C
C     --------------------------------------------------------------
C
      SUBROUTINE SCNINT (IFILE,JFILE,DVINT,JEMIT,NPTS,IBUF)               I06310
C                                                                         I06320
      IMPLICIT REAL*8          (V)                                     ! I06330
C                                                                         I06340
C**********************************************************************   I06350
C                                                                         I06360
C     INTERPOLATION FUNCTION DRIVER FOR WEIGHTING FUNCTIONS               I06370
C                                                                         I06380
C     FOUR-POINT VERSION    (MARCH 1990)                                  I06390
C                                                                         I06400
C**********************************************************************   I06410
C                                                                         I06420
C     THE INPUT DATA WILL BE PUT INTO T(5) = S(1) WITH THE LAST           I06430
C     4 POINTS OF THE PREVIOUS PANEL PUT INTO T(1 TO 4).                  I06440
C     THIS SCHEME PERMITS 6 POINT INTERPOLATION.                          I06450
C                                                                         I06460
C     S IS NOMINALLY 2401 POINTS BUT MAY NEED TO BE EXTENDED BY           I06470
C     2 POINTS TO PERMIT 4 POINT INTERPOLATION UP TO THE LAST             I06480
C     DATA POINT.                                                         I06490
C                                                                         I06500
      COMMON T(2410),R(2401)                                              I06510
      DIMENSION S(2406)                                                   I06520
      EQUIVALENCE (T(5),S(1))                                             I06530
C                                                                         I06540
      character*8      XID,       HMOLID,      YID        
      real*8               SECANT,       XALTZ 
C                                                                         I06560
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       I06570
     *                WK(60),PZL,PZU,TZL,TZU,WN2   ,DV ,V1C,V2C,TBOUND,   I06580
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    I06590
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   I06600
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       I06610
     *               DVSC,XDUM,V1SHFT                                     I06620
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 I06630
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            I06640
      COMMON /INPNL/ V1I,V2I,DVI,NNI                                      I06650
      COMMON /OUTPNL/ V1J,V2J,DVJ,NNJ                                     I06660
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I06670
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I06680
     *              NLTEFL,LNFIL4,LNGTH4                                  I06690
C                                                                         I06700
      DIMENSION FILHDR(2),RSTAT(3)                                        I06710
C                                                                         I06720
      EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),                 I06730
     *            (FSCDID(6),ISCHDR) , (FSCDID(12),XSCID),                I06740
     *            (FSCDID(13),XHWHM) , (FSCDID(14),IDABS),                I06750
     *            (FSCDID(16),LAYR1)                                      I06760
C                                                                         I06770
      CHARACTER*12 BCD,HTRANS,HABSRB,HRADIA                               I06780
C                                                                         I06790
      DATA HTRANS / 'TRANSMISSION'/,HABSRB / ' ABSORPTION '/,             I06800
     *     HRADIA / ' RADIANCE   '/                                       I06810
C                                                                         I06820
C----------------------------------------------------------------------   I06830
C     JEMIT=-1  INTERPOLATE ABSORPTION                                    I06840
C     JEMIT=0   INTERPOLATE TRANSMISSION                                  I06850
C     JEMIT=1   INTERPOLATE EMISSION                                      I06860
C     JEMIT=2   INTERPOLATE OPTICAL DEPTH                                 I06870
C----------------------------------------------------------------------   I06880
C                                                                         I06890
      WRITE (IPR,900)                                                     I06900
      CALL CPUTIM (TIME1)                                                 I06910
      TIMRDF = 0.0                                                        I06920
      TIMCNV = 0.0                                                        I06930
      TIMPNL = 0.0                                                        I06940
C                                                                         I06950
      V1SAV = V1                                                          I06960
      V2SAV = V2                                                          I06970
      DVSAV = DV                                                          I06980
      DVOSAV = 0.                                                         I06990
C                                                                         I07000
      DVO = DVINT                                                         I07010
C                                                                         I07020
C     I4PT = 1 FOR FOUR POINT INTERPOLATION                               I07030
C                                                                         I07040
      I4PT = 1                                                            I07050
      ICNVRT = 1                                                          I07060
      IF (DVO.LE.0.) GO TO 40                                             I07070
C                                                                         I07080
      IF (IBUF.EQ.1) REWIND IFILE                                         I07090
      REWIND JFILE                                                        I07100
C                                                                         I07110
C     BUFFER IN THE FILE HEADER ON UNIT (IFILE)                           I07120
C     BUFFER OUT ON UNIT (JFILE)                                          I07130
C                                                                         I07140
      IF (IBUF.EQ.1) THEN                                                 I07150
         CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)                         I07160
         IF (IEOF.EQ.0) GO TO 30                                          I07170
         JABS = 0                                                         I07180
         IDABS = 0                                                        I07190
         IF (JEMIT.LT.0) THEN                                             I07200
            JABS = 1                                                      I07210
            JEMIT = 0                                                     I07220
            IDABS = -1                                                    I07230
         ENDIF                                                            I07240
      ENDIF                                                               I07250
      V1 = V1C                                                            I07260
      V2 = V2C                                                            I07270
      DVI = DV                                                            I07280
C                                                                         I07290
C   V2 IS ONLY APPROXIMATE                                                I07300
C                                                                         I07310
      NUM = (((V2-V1)/DVO)+0.5)                                           I07320
      V2 = V1+FLOAT(NUM)*DVO                                              I07330
      NUM = NUM+1                                                         I07340
      WRITE (IPR,905) V1,V2,DVO,NUM,JEMIT,I4PT,IFILE,JFILE,NPTS           I07350
C                                                                         I07360
      ISCAN = ISCHDR                                                      I07370
      IF (ISCAN.LE.0.OR.XSCID.EQ.-99.) ISCAN = 0                          I07380
      IF (ISCHDR.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCHDR                   I07390
      ISCHDR = ISCAN+10                                                   I07400
      V1C = V1                                                            I07410
      V2C = V2                                                            I07420
      DV = DVO                                                            I07430
C                                                                         I07440
      SCNID = 100*JEMIT                                                   I07450
      XSCID = SCNID+0.01                                                  I07460
C                                                                         I07470
      CALL BUFOUT (JFILE,FILHDR(1),NFHDRF)                                I07480
C                                                                         I07490
      JTREM = -1                                                          I07500
      IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0                        I07510
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2                        I07520
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.2)) JTREM = 2                        I07530
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1                        I07540
      ISCANT = MOD(ISCAN,1000)                                            I07550
      IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 2                       I07560
      IF (JTREM.LT.0) THEN                                                I07570
         WRITE(IPR,*) ' JTREM.LT.0 AT I07570'                             I07572
         STOP         ' JTREM.LT.0 AT I07570'
      ENDIF                                                                     
      WRITE (IPR,910) IFILE,IEMIT,JEMIT,JTREM,JABS                        I07580
C                                                                         I07590
      IDATA = -1                                                          I07600
C                                                                         I07610
C     NEED TO SAVE LAST IBOUND POINTS OF EACH PANEL TO ATTACH TO NEXT     I07620
C                                                                         I07630
      IBOUND = 4                                                          I07640
C                                                                         I07650
C     VBOT IS LOWEST NEEDED WAVENUMBER, VTOP IS HIGHEST                   I07660
C                                                                         I07670
      BOUND = FLOAT(IBOUND)*DV                                            I07680
      VBOT = V1-BOUND                                                     I07690
      VTOP = V2+BOUND                                                     I07700
C                                                                         I07710
      IF (JEMIT.EQ.0.AND.IDABS.EQ.0) BCD = HTRANS                         I07720
      IF (JEMIT.EQ.0.AND.IDABS.EQ.-1) BCD = HABSRB                        I07730
      IF (JEMIT.EQ.1) BCD = HRADIA                                        I07740
      IF (NPTS.GT.0) WRITE (IPR,915) BCD                                  I07750
C                                                                         I07760
C     ZERO OUT T(1 TO IBOUND)                                             I07770
C                                                                         I07780
      DO 10 II = 1, IBOUND                                                I07790
         T(II) = 0.0                                                      I07800
   10 CONTINUE                                                            I07810
C                                                                         I07820
C     READ FROM IFILE UNTIL THE FIRST REQUIRED POINT IS REACHED           I07830
C     AND LOAD DATA INTO S                                                I07840
C                                                                         I07850
      CALL RDPANL (S,JTREM,IFILE,ISCAN,JEMIT,ICNVRT)                      I07860
      IF (IEOFSC.LE.0) GO TO 20                                           I07870
C                                                                         I07880
C     DO INTERPOLATION                                                    I07890
C                                                                         I07900
      CALL INTERP (IFILE,JFILE,I4PT,IBOUND,NPTS,JTREM,ISCAN,JEMIT,        I07910
     *             RSTAT,ICNVRT)                                          I07920
C                                                                         I07930
      CALL CPUTIM (TIME2)                                                 I07940
      CALL ENDFIL (JFILE)                                                 I07950
C                                                                         I07960
C     WRITE STATISTICS                                                    I07970
C                                                                         I07980
      WRITE (IPR,920) RSTAT(1),RSTAT(2),RSTAT(3)                          I07990
      TIMTOT = TIME2-TIME1                                                I08000
      TIMCNV = TIMTOT-TIMRDF-TIMPNL                                       I08010
      WRITE (IPR,925) TIMTOT,TIMRDF,TIMCNV,TIMPNL                         I08020
C                                                                         I08030
      GO TO 30                                                            I08040
C                                                                         I08050
   20 CONTINUE                                                            I08060
      WRITE (IPR,930) IFILE                                               I08070
C                                                                         I08080
   30 CONTINUE                                                            I08090
      V1 = V1SAV                                                          I08100
      V2 = V2SAV                                                          I08110
      DV = DVSAV                                                          I08120
      RETURN                                                              I08130
C                                                                         I08140
   40 CONTINUE                                                            I08150
      WRITE (IPR,935) DVINT                                               I08160
C                                                                         I08170
      RETURN                                                              I08180
C                                                                         I08190
  900 FORMAT (/,'0***SCNINT***',/)                                        I08200
  905 FORMAT (5X,'V1=',F14.8,' V2=',F14.8,' DVO=',E14.6,' NUM=',I8,/,     I08210
     *        5X,'JEMIT=',I3,' I4PT=',I3,' IUNIT=',I3,' JUNIT=',I3,       I08220
     *        ' NPTS=',I5)                                                I08230
  910 FORMAT (5X,'INPUT FILE NUMBER =',I3,' IEMIT=',I3,' JEMIT=',I3,      I08240
     *        ' JTREM=',I3,' JABS=',I3)                                   I08250
  915 FORMAT (///,'0',5X,A12,/)                                           I08260
  920 FORMAT ('0    SUMOUT =',1P,E16.9,'  MIN =',E16.9,'  MAX =',E16.9)   I08270
  925 FORMAT (/,5X,'SCNINT TIME: TOTAL = ',F8.3,' READ = ',F8.3,          I08280
     *        ' INTERP = ',F8.3,' WRITE = ',F8.3,/)                       I08290
  930 FORMAT (/,5X,'SCNINT- ERROR: EOF ON INPUT UNIT ',I4,                I08300
     *        ' BEFORE V1 WAS REACHED',/)                                 I08310
  935 FORMAT (/,5X,'SCNINT- ERROR: DVINT .LT. ZERO ; DVINT =',F12.4,/)    I08320
C                                                                         I08330
      END                                                                 I08340
C
C     --------------------------------------------------------------
C
      SUBROUTINE SCNMRG (IFILE,JFILE)                                     I08350
C                                                                         I08360
      IMPLICIT REAL*8          (V)                                     ! I08370
C                                                                         I08380
C     DRIVER FOR CONVOLVING INSTRUMENTAL SCANNING FUNCTION                I08390
C     WITH SPECTRUM                                                       I08400
C                                                                         I08410
      COMMON S(3850),R1(5000)                                             I08420
C                                                                         I08430
      character*8      XID,       HMOLID,      YID,SCANID
      real*8               SECANT,       XALTZ 
C                                                                         I08450
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       I08460
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND,   I08470
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    I08480
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   I08490
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       I08500
     *               DVSC,XDUM,V1SHFT                                     I08510
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 I08520
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            I08530
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      I08540
      COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,
     *                HWF3,DXF3,NX3,N3MAX
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I08570
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I08580
     *              NLTEFL,LNFIL4,LNGTH4                                  I08590
      COMMON /SCINF/ HWHM,JEMIT,JFN,SAMPLE,SCANID,NPTS,XF(6018)           I08600
      COMMON /RCTSV/ JJ,SUMJ,JFLG,RNJ,NB,IPC,VLFT,VCNT,VRGT,WGTL,WGTR     I15280
C                                                                         I08610
      DIMENSION FILHDR(2)                                                 I08620
      DIMENSION SUMR(4)                                                   I08630
C                                                                         I08640
      EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),                 I08650
     *            (FSCDID(6),ISCHDR) , (FSCDID(12),XSCID),                I08660
     *            (FSCDID(13),XHWHM) , (FSCDID(14),IDABS),                I08670
     *            (FSCDID(16),LAYR1)                                      I08680
C                                                                         I08690
      CHARACTER*12 BCD,HTRANS,HABSRB,HRADIA                               I08700
C                                                                         I08710
      DATA HTRANS / 'TRANSMISSION'/,HABSRB / ' ABSORPTION '/,             I08720
     *     HRADIA / ' RADIANCE   '/                                       I08730
C                                                                         I08740
C     IUNIT INPUT FILE                                                    I08750
C     JUNIT OUTPUT FILE                                                   I08760
C                                                                         I08770
      IUNIT = IFILE                                                       I08780
      JUNIT = JFILE                                                       I08790
      NREN = 0                                                            I08800
      IPRT = 1                                                            I08810
      IDABS = 0                                                           I08820
      IF (JEMIT.LT.0) THEN                                                I08830
         JABS = 1                                                         I08840
         JEMIT = 0                                                        I08850
         IDABS = -1                                                       I08860
      ENDIF                                                               I08870
      IDABST = IDABS                                                      I08880
      IFILST = 1                                                          I08890
      NIFILS = 9999                                                       I08900
C                                                                         I08910
      SUMOUT = 0.                                                         I08920
      SMIN = 999999.                                                      I08930
      SMAX = -99999.                                                      I08940
      DVOSAV = 0.                                                         I08950
      SUMR(1) = SUMOUT                                                    I08960
      SUMR(2) = SMIN                                                      I08970
      SUMR(3) = SMAX                                                      I08980
      SUMR(4) = DVOSAV                                                    I08990
      NSHIFT = 32
C                                                                         I09000
      REWIND IUNIT                                                        I09010
      CALL BUFIN (IUNIT,IEOF,FILHDR(1),NFHDRF)                            I09020
      IF (IEOF.EQ.0) GO TO 50                                             I09030
C                                                                         I09040
      DVSAV = DV                                                          I09050
      IDABS = IDABST                                                      I09060
C                                                                         I09070
      WRITE (IPR,900) XID,(YID(M),M=1,2)                                  I09080
      WRITE (IPR,905) LAYR1,LAYER                                         I09090
      WRITE (IPR,910) SECANT,PAVE,TAVE,DV,V1C,V2C                         I09100
      WRITE (IPR,915) WBROAD,(HMOLID(M),WK(M),M=1,NMOL)                   I09110
C                                                                         I09120
      ISCAN = ISCHDR                                                      I09130
      IF (ISCAN.LE.0.OR.XSCID.EQ.-99.) ISCAN = 0                          I09140
      IF (ISCHDR.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCHDR                   I09150
      ISCHDR = ISCAN+1                                                    I09160
      JTREM = -1                                                          I09170
      IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0                        I09180
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2                        I09190
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1                        I09200
      ISCANT = MOD(ISCAN,1000)                                            I09210
      IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 2                       I09220
C
      IF (JTREM.LT.0) THEN                                                I09230
         WRITE(IPR,*) ' SCANF; JTREM LT 0'                                I09230
         STOP         ' SCANF; JTREM LT 0'                                I09232
      ENDIF                                                      
C
      WRITE (IPR,920) SCANID,IUNIT,IFILST,NIFILS,JEMIT,JFN,JVAR,JABS      I09240
C                                                                         I09250
C     JTREM=0   SCANFN CONVOLVED WITH EXPONENTIATED                       I09260
C                      ABSORPTION COEFFICIENT                             I09270
C     JTREM=1   SCANFN CONVOLVED WITH EMISSION                            I09280
C     JTREM=2   SCANFN CONVOLVED WITH TRANSMISSION                        I09290
C                                                                         I09300
      DVI = DV                                                            I09310
      DVO = HWHM/SAMPLE                                                   I09320
C                                                                         I09330
C    BOUND AT THIS POINT IS THE WAVENUMBER VALUE                          I09340
C    OF HALF THE SCANNING FUNCTION                                        I09350
C                                                                         I09360
      BOUND = HWF*HWHM                                                    I09370
      DV = DVO                                                            I09380
      V1C = V1                                                            I09390
      V2C = V2                                                            I09400
      SCIND = JVAR+10*(JFN+10*(JEMIT))                                    I09410
      XSCID = SCIND+0.01                                                  I09420
      XHWHM = HWHM                                                        I09430
      CALL BUFOUT (JUNIT,FILHDR(1),NFHDRF)                                I09440
      WRITE (IPR,925) HWHM,BOUND,JUNIT,V1,V2,DVO                          I09450
      NBOUND = (2.*HWF)*SAMPLE+0.01                                       I09460
C                                                                         I09470
C     BOUND AT THIS POINT IS THE WAVENUMBER VALUE OF THE                  I09480
C     FULL SCANNING FUNCTION                                              I09490
C                                                                         I09500
      BOUND = FLOAT(NBOUND)*DVO/2.                                        I09510
      MAXF = NLIMF+2*NBOUND+NSHIFT                                        I09520
C                                                                         I09530
      TIMRDF = 0.                                                         I09540
      TIMCNV = 0.                                                         I09550
      TIMPNL = 0.                                                         I09560
      IEOFSC = 1                                                          I09570
      NLO = NSHIFT+1                                                      I09580
      SUMIN = 0.                                                          I09590
      NHI = NLIMF+NSHIFT-1                                                I09600
      DO 10 I = 1, MAXF                                                   I09610
         R1(I) = 0.                                                       I09630
   10 CONTINUE                                                            I09640
      INIT = 0                                                            I09650
      IDATA = -1                                                          I09660
      IPANEL = -1                                                         I09663
      JFLG = -1                                                           I09667
      VFT = V1-FLOAT(NSHIFT)*DV                                           I09670
      VBOT = V1-BOUND                                                     I09680
      VTOP = V2+BOUND                                                     I09690
C                                                                         I09700
      IF (JEMIT.EQ.0.AND.IDABS.EQ.0) BCD = HTRANS                         I09710
      IF (JEMIT.EQ.0.AND.IDABS.EQ.-1) BCD = HABSRB                        I09720
      IF (JEMIT.EQ.1) BCD = HRADIA                                        I09730
      IF (NPTS.GT.0) WRITE (IPR,930) BCD                                  I09740
   20 CALL CPUTIM (TIME0)                                                 I09750
      IF (IEOFSC.LE.0) GO TO 40                                           I09760
      CALL RDSCAN (S,JTREM,IUNIT,ISCAN,IPRT)                              I09770
C                                                                         I09780
CPRT  WRITE(IPR,935) IEOFSC,IDATA                                         I09790
C                                                                         I09800
      CALL CPUTIM (TIME)                                                  I09810
      TIMRDF = TIMRDF+TIME-TIME0                                          I09820
C                                                                         I09830
      IF (IEOFSC.LE.0) GO TO 40                                           I09840
      IF (JFN.NE.0) CALL SHRKSC (INIT,HWHM)                               I09850
C                                                                         I09860
C     SHRKSC MAY SHRINK (COMPRESS) THE DATA;                              I09870
C     DVI IS MODIFIED ACCORDINGLY                                         I09880
C                                                                         I09890
   30 CONTINUE                                                            I09900
      IF (JFN.EQ.0) THEN                                                  I09910
         CALL CNVRCT (S,HWHM,R1,XF)                                       I09920
      ELSEIF (JFN.EQ.5) THEN
         CALL CNVVRC (S,HWHM,R1,XF)                                       I09920
      ELSEIF (JFN.EQ.6) THEN
         CALL CNVVRL (S,HWHM,R1,XF)                                       I09920
      ELSE                                                                I09930
         CALL CONVSC (S,HWHM,R1,XF)                                       I09940
      ENDIF                                                               I09950
C                                                                         I09960
CPRT  WRITE(IPR,935) IEOFSC,IDATA,IPANEL                                  I09970
C                                                                         I09980
      IF (IPANEL.EQ.0) GO TO 20                                           I09990
C                                                                         I10000
   40 CONTINUE                                                            I10010
      IF (JFN.EQ.0.OR.JFN.EQ.5.OR.JFN.EQ.6) THEN                          I10020
         CALL PNLRCT (R1,JUNIT,SUMR,NPTS)                                 I10030
      ELSE                                                                I10040
         CALL PANLSC (R1,JUNIT,SUMR,NPTS)                                 I10050
      ENDIF                                                               I10060
      IF ((ISTOP.NE.1).AND.(IEOFSC.GT.0)) GO TO 30                        I10070
      CALL CPUTIM (TIME)                                                  I10080
      WRITE (IPR,940) TIME,TIMRDF,TIMCNV,TIMPNL                           I10090
C                                                                         I10100
      SUMIN = SUMIN*DVSAV                                                 I10110
C                                                                         I10120
      WRITE (IPR,945) SUMIN                                               I10130
C                                                                         I10140
      IF (IEOFSC.EQ.1) CALL SKIPFL (1,IUNIT,IEOFSC)                       I10150
C                                                                         I10160
      IEOFT = IEOFT+1                                                     I10170
C                                                                         I10180
C                                                                         I10190
      SUMOUT = SUMR(1)                                                    I10200
      SMIN = SUMR(2)                                                      I10210
      SMAX = SUMR(3)                                                      I10220
      DVOSAV = SUMR(4)                                                    I10230
C                                                                         I10240
      SUMOUT = SUMOUT*DVOSAV                                              I10250
      WRITE (IPR,950) SUMOUT,SMIN,SMAX                                    I10260
C                                                                         I10270
   50 RETURN                                                              I10280
C                                                                         I10290
  900 FORMAT ('0',' **SCNMRG** ',/,'0',10A8,2X,2(1X,A8,1X))               I10300
  905 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)            I10310
  910 FORMAT ('0 SECANT =',F15.5,/'0 PRESS(MB) =',F12.5/'0 TEMP(K) =',    I10320
     *        F11.2,/'0 DV(CM-1) = ',F12.8,/'0 V1(CM-1) = ',F12.6,/,      I10330
     *        '0 V2(CM-1) = ',F12.6)                                      I10340
  915 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//5X,'WBROAD = ',       I10350
     *        1PE10.3,/(5X,A6,' = ',1PE10.3))                             I10360
  920 FORMAT ('0','***',A8,'***',//6X,'INPUT FILE NUMBER =',I3,           I10370
     *        ' ,IFILST = ',I5,' ,NIFILS = ',I5,',JEMIT =',I2,            I10380
     *        ' ,JFN =',I2,' ,JVAR =',I2,'  ,JABS =',I2)                  I10390
  925 FORMAT (1X,'     HWHM OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,     I10400
     *        5X,'BOUND OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,6X,      I10410
     *        'OUTPUT FILE NUMBER =',I3,',   V1 =',F12.5,',   V2 =',      I10420
     *        F12.5,5X,' DV OUT',F12.8)                                   I10430
  930 FORMAT (///,'0',5X,A12,/)                                           I10440
  935 FORMAT ('0',5X,'IEOFSC =',I3,'  IDATA =',I3,'  IPANEL =',I3,/)      I10450
  940 FORMAT ('0',5X,'TIME =',F7.3,',  READ =',F6.3,',  CONV. =',F7.3,    I10460
     *        ',  PANEL =',F6.3)                                          I10470
  945 FORMAT ('0    SUMIN  =',1P,E16.9)                                   I10480
  950 FORMAT ('0    SUMOUT =',1P,E16.9,'  MIN =',E16.9,'  MAX =',E16.9)   I10490
C                                                                         I10500
      END                                                                 I10510
C
C     --------------------------------------------------------------
C
      SUBROUTINE SHRKSC (INIT,HWHM)                                       I10520
C                                                                         I10530
      IMPLICIT REAL*8          (V)                                     ! I10540
C                                                                         I10550
C     THIS SUBROUTINE COMPRESSES (SHRINKS) THE INPUT TO THE CONVOLUTION   I10560
C     ROUTINE FOR THE SCANNING FUNCTION TO ACCELERATE THE CALCULATION     I10570
C                                                                         I10580
      COMMON S(3850),R1(5000),SS(200)                                     I10590
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   I10600
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       I10610
     *               DVSC,XDUM,V1SHFT                                     I10620
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            I10630
      COMMON /RSCAN/ V1I,V2I,DVI,NLIM                                     I10640
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I10650
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I10660
     *              NLTEFL,LNFIL4,LNGTH4                                  I10670
      DIMENSION JRATIO(24)                                                I10680
C                                                                         I10690
      DATA JRATIO / 1,2,3,4,5,6,8,10,12,15,16,20,24,25,30,32,40,48,50,    I10700
     *             60,75,80,100,120 /                                     I10710
C                                                                         I10720
      CALL CPUTIM (TIME0)                                                 I10730
      NLIMS = NLIM                                                        I10740
      IF (NREN.GT.0) THEN                                                 I10750
         DO 10 I = 1, NREN                                                I10760
            S(I) = SS(I)                                                  I10770
   10    CONTINUE                                                         I10780
      ENDIF                                                               I10790
      NREN = 0                                                            I10800
      IF (INIT.EQ.0) THEN                                                 I10810
         DVSC = HWHM/12.                                                  I10820
         IRATSH = DVSC/DVI+0.5                                            I10830
         DO 20 I = 2, 24                                                  I10840
            IF (JRATIO(I).GT.IRATSH) THEN                                 I10850
               IRATSH = JRATIO(I-1)                                       I10860
               GO TO 30                                                   I10870
            ENDIF                                                         I10880
   20    CONTINUE                                                         I10890
   30    IF (IRATSH.GT.JRATIO(24)) IRATSH = JRATIO(24)                    I10900
         IF (IRATSH.LE.1) RETURN                                          I10910
         DVSC = FLOAT(IRATSH)*DVI                                         I10920
         V1SHFT = FLOAT(IRATSH-1)*DVI/2.                                  I10930
         WRITE (IPR,900) IRATSH                                           I10940
         SRATIO = IRATSH                                                  I10950
         IRATM1 = IRATSH-1                                                I10960
         INIT = 1                                                         I10970
      ENDIF                                                               I10980
      IF (IRATSH.LE.1) RETURN                                             I10990
      NREN = NLIM-(NLIM/IRATSH)*IRATSH                                    I11000
C                                                                         I11010
CPRT  WRITE(IPR,905) V1I,V1SHFT,DVSC,NREN                                 I11020
C                                                                         I11030
      V1I = V1I+V1SHFT                                                    I11040
      IMIN = 1                                                            I11050
      IMAX = NLIM-IRATM1-NREN                                             I11060
C                                                                         I11070
      K = 0                                                               I11080
      DO 50 I = IMIN, IMAX, IRATSH                                        I11090
         SUMK = 0.                                                        I11100
         JHI = I+IRATM1                                                   I11110
         K = K+1                                                          I11120
         DO 40 J = I, JHI                                                 I11130
            SUMK = SUMK+S(J)                                              I11140
   40    CONTINUE                                                         I11150
         S(K) = SUMK/SRATIO                                               I11160
   50 CONTINUE                                                            I11170
C                                                                         I11180
      V2I = V1I+DVSC*FLOAT(K-1)                                           I11190
      NLIM = K                                                            I11200
      DVI = DVSC                                                          I11210
      ILO = ((VBOT-V1I)/DVI)+1.5                                          I11220
      ILO = MAX(ILO,1)                                                    I11230
      IHI = ((VTOP-V1I)/DVI)+1.5                                          I11240
      IHI = MIN(IHI,NLIM)                                                 I11250
C                                                                         I11260
CPRT  WRITE(IPR,910) ILO,IHI                                              I11270
C                                                                         I11280
      IF (NREN.GT.0) THEN                                                 I11290
         DO 60 I = 1, NREN                                                I11300
            II = NLIMS-NREN+I                                             I11310
            SS(I) = S(II)                                                 I11320
   60    CONTINUE                                                         I11330
      ENDIF                                                               I11340
      CALL CPUTIM (TIME)                                                  I11350
      TIMCNV = TIMCNV+TIME-TIME0                                          I11360
C                                                                         I11370
      RETURN                                                              I11380
C                                                                         I11390
  900 FORMAT ('   SHRINK RATIO = ',I5)                                    I11400
  905 FORMAT ('   V1I =',F10.3,'  V1SHFT =',F10.3,'  DVSC =',F12.5,       I11410
     C        '   NREN =',I4)                                             I11420
  910 FORMAT ('   ILO =',I4,'  IHI =',I4)                                 I11430
C                                                                         I11440
      END                                                                 I11450
C
C     --------------------------------------------------------------
C
      SUBROUTINE SHAPET (XF)                                              I11460
C                                                                         I11470
C     SUBROUTINE SHAPET SETS UP THE TRIANGULAR SCANNING FUNCTION          I11480
C                                                                         I11490
      COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,
     *                HWF3,DXF3,NX3,N3MAX
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I11520
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I11530
     *              NLTEFL,LNFIL4,LNGTH4                                  I11540
      DIMENSION XF(*)                                                     I11550
C                                                                         I11560
      XTRIAN(X) = 1.-0.5*X                                                I11570
      DO 10 I = 1, NFMAX                                                  I11580
         XF(I) = 0.                                                       I11590
   10 CONTINUE                                                            I11600
      XF(1) = 0.5                                                         I11610
      SUM = XF(1)                                                         I11620
      DO 20 I = 2, NF                                                     I11630
         X = FLOAT(I-1)*DXF                                               I11640
         XF(I) = 0.5*XTRIAN(X)                                            I11650
         SUM = SUM+2.*XF(I)                                               I11660
   20 CONTINUE                                                            I11670
      SUM = SUM*DXF                                                       I11680
C                                                                         I11690
CPRT  WRITE(IPR,900) NF,DXF,SUM                                           I11700
C                                                                         I11710
      RETURN                                                              I11720
C                                                                         I11730
  900 FORMAT ('0',5X,'NF =',I5,',  DXF =',F7.5,',    SUM =',F18.15)       I11740
C                                                                         I11750
      END                                                                 I11760
C
C     --------------------------------------------------------------
C
      SUBROUTINE RDSCAN (S,JTREM,IFILE,ISCAN,IPRT)                        I11770
C                                                                         I11780
      IMPLICIT REAL*8          (V)                                     ! I11790
C                                                                         I11800
C     SUBROUTINE RDSCAN INPUTS PANELS FROM IFILE RESULTING                I11810
C     FROM THE LBLRTM CALCULATION FOR CONVOLUTION                         I11820
C     WITH THE SELECTED SCANNING FUNCTION                                 I11830
C                                                                         I11840
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           I11850
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   I11860
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       I11870
     *               DVSC,XDUM,V1SHFT                                     I11880
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 I11890
      COMMON /RSCAN/ VMIN,VMAX,DVI,NNI                                    I11900
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I11910
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I11920
     *              NLTEFL,LNFIL4,LNGTH4                                  I11930
      DIMENSION DUMMY(2),PNLHDR(2)                                        I11940
      DIMENSION S(*)                                                      I11950
C                                                                         I11960
      EQUIVALENCE (PNLHDR(1),VMIN)                                        I11970
C                                                                         I11980
CPRT  WRITE(IPR,900) VBOT,VTOP                                            I11990
C                                                                         I12000
      IDUM1 = 0                                                           I12010
      IDUM2 = 0                                                           I12020
      ISCANT = MOD(ISCAN,1000)                                            I12030
      IF(JTREM.EQ.0.AND.ISCANT.GE.1) GO TO 60                             I12035
      IF (ISCAN.LT.1) THEN                                                I12040
         IF (JTREM.EQ.1) IDUM1 = 1                                        I12050
         IF (JTREM.EQ.2) IDUM2 = 1                                        I12060
      ENDIF                                                               I12070
   10 CALL BUFIN (IFILE,IEOFSC,PNLHDR(1),NPHDRF)                          I12080
      IF (IEOFSC.LE.0) GO TO 50                                           I12090
      NLOW = NREN+1                                                       I12100
      IF (NREN.LE.0) NLOW = 1                                             I12110
      VMIN = VMIN-(NLOW-1)*DVI                                            I12120
      NNB = NNI                                                           I12130
      NNI = NNI+NLOW-1                                                    I12140
      IF ((IDATA.EQ.-1).AND.(VMIN.GT.VBOT).AND.(IPRT.EQ.1))               I12150
     *     WRITE (IPR,905)                                                I12160
      IDATA = 0                                                           I12170
      IF (VMAX.GE.VBOT) GO TO 20                                          I12180
      IF (IDUM2.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),1)                I12190
      CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)                                I12200
      IF (IDUM1.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),1)                I12210
      GO TO 10                                                            I12220
   20 IF (JTREM.EQ.0) THEN                                                I12230
         CALL BUFIN (IFILE,IEOFSC,S(NLOW),NNB)                            I12240
         DO 30 I = NLOW, NNI                                              I12250
            SI = S(I)                                                     I12260
            S(I) = 1.                                                     I12270
            IF (SI.GT.1.0E-04) THEN                                       I12280
               IF (SI.LT.ARGMIN) THEN                                     I12290
                  S(I) = EXP(-SI)                                         I12300
               ELSE                                                       I12310
                  S(I) = EXPMIN                                           I12320
               ENDIF                                                      I12330
            ELSE                                                          I12340
               S(I) = 1.-SI                                               I12350
            ENDIF                                                         I12360
   30    CONTINUE                                                         I12370
      ELSE                                                                I12380
C                                                                         I12390
         IF (IDUM2.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),1)             I12400
         CALL BUFIN (IFILE,IEOFSC,S(NLOW),NNB)                            I12410
         IF (IDUM1.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),1)             I12420
      ENDIF                                                               I12430
C                                                                         I12440
CPRT  WRITE(IPR,910) VMIN,VMAX,DVI,NLOW,NNI                               I12450
C                                                                         I12460
      IF (JABS.NE.0) THEN                                                 I12470
         DO 40 I = NLOW, NNI                                              I12480
            S(I) = 1.-S(I)                                                I12490
   40    CONTINUE                                                         I12500
      ENDIF                                                               I12510
      ILO = 1                                                             I12520
      IHI = NNI                                                           I12530
      DIF = (VMIN-VBOT)/DVI                                               I12540
      IF (DIF.LT.0.) ILO = -DIF+1.5                                       I12550
      IF (VMAX.LE.VTOP) RETURN                                            I12560
      IHI = (VTOP-VMIN)/DVI+1.5                                           I12570
      IDATA = 1                                                           I12580
      RETURN                                                              I12590
   50 IF (IPRT.EQ.1) WRITE (IPR,915)                                      I12600
      RETURN                                                              I12610
C                                                                         I12620
   60 WRITE(IPR,920) JTREM,ISCAN                                          I12622
      RETURN                                                              I12625
C                                                                         I12627
  900 FORMAT ('0',/,'0   READING SPECTRUM, VBOT =',F10.3,', VTOP =',      I12630
     *        F10.3)                                                      I12640
  905 FORMAT ('0 ********** FIRST VALUE USED ON IFILE; CHECK IFILE ')     I12650
  910 FORMAT (10X,'VMIN =',F10.3,',  VMAX =',F10.3,',  DVI=',F7.5,',      I12660
     *        NLOW=',I4,',  NNI=',I4)                                     I12670
  915 FORMAT ('0 ********** END OF FILE ENCOUNTERED; CHECK IFILE ')       I12680
  920 FORMAT(' ERROR IN INPUT',/,'  JTREM =',I2,'  ISCAN=',I5)            I12685
C                                                                         I12690
      END                                                                 I12700
C
C     --------------------------------------------------------------
C
      SUBROUTINE CONVSC (S,HWHMV1,R1,XF)                                  I12710
C                                                                         I12720
      IMPLICIT REAL*8          (V)                                     ! I12730
C                                                                         I12740
C     SUBROUTINE CONVSC PERFORMS THE CONVOLUTION WITH THE SELECTED        I12750
C     SCANNING FUNCTION                                                   I12760
C                                                                         I12770
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   I12780
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       I12790
     *               DVSC,XDUM,V1SHFT                                     I12800
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 I12810
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            I12820
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      I12830
      COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,
     *                HWF3,DXF3,NX3,N3MAX
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I12860
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I12870
     *              NLTEFL,LNFIL4,LNGTH4                                  I12880
      DIMENSION S(*),R1(*),XF(*)                                          I12890
C                                                                         I12900
      CALL CPUTIM (TIME0)                                                 I12910
      IF (ILO.GT.IHI) GO TO 60                                            I12920
      RATIO = DVI/DVO                                                     I12930
      DVODX = DVO/DXF                                                     I12940
      HWBND = HWF/DVO                                                     I12950
      ZINT = ((V1I-VFT)/DVO)                                              I12960
      HWHM = HWHMV1                                                       I12970
      ITST = -1                                                           I12980
C                                                                         I12990
      DO 50 I = ILO, IHI                                                  I13000
         IF (S(I).EQ.0.) GO TO 50                                         I13010
         IF (I.LT.ITST) GO TO 20                                          I13020
         ITST = 9999                                                      I13030
         IF (JVAR.EQ.0) GO TO 10                                          I13040
         VI = FLOAT(I-1)*DVI+V1I                                          I13050
         HWHM = HWHMV1*(VI/V1)**2                                         I13060
         ITST = I+IFIX(1./DVI)                                            I13070
   10    CONTINUE                                                         I13080
         ZSLOPE = DVODX/HWHM                                              I13090
         ZBOUND = HWBND*HWHM                                              I13100
         XNORM = DVI/HWHM                                                 I13110
C                                                                         I13120
CPRT     WRITE(IPR,900) VI,HWHM                                           I13130
C                                                                         I13140
   20    CONTINUE                                                         I13150
         ZPEAK = FLOAT(I-1)*RATIO+ZINT                                    I13160
         JMAX = ZPEAK+ZBOUND+1.5                                          I13170
         IF (JMAX.LE.MAXF) GO TO 30                                       I13180
         ILAST = I-1                                                      I13190
         GO TO 60                                                         I13200
C                                                                         I13210
   30    JMIN = ZPEAK-ZBOUND+1.5                                          I13220
         JMIN = MAX(JMIN,1)                                               I13230
         SUMIN = SUMIN+S(I)                                               I13240
         SI = XNORM*S(I)                                                  I13250
         ZF = (FLOAT(JMIN-1)-ZPEAK)*ZSLOPE                                I13260
         DO 40 JF = JMIN, JMAX                                            I13270
            IT = ABS(ZF)+1.5                                              I13280
            R1(JF) = R1(JF)+SI*XF(IT)                                     I13290
            ZF = ZF+ZSLOPE                                                I13300
   40    CONTINUE                                                         I13310
C                                                                         I13320
   50 CONTINUE                                                            I13330
      ILAST = IHI                                                         I13340
      IPANEL = IDATA                                                      I13350
      GO TO 70                                                            I13360
C                                                                         I13370
   60 IPANEL = 1                                                          I13380
   70 CALL CPUTIM (TIME)                                                  I13390
      TIMCNV = TIMCNV+TIME-TIME0                                          I13400
      ILO = ILAST+1                                                       I13410
C                                                                         I13420
      RETURN                                                              I13430
C                                                                         I13440
  900 FORMAT ('0 AVE PANEL WAVENUMBER = ',F12.4,5X,'HWHM = ',F10.5)       I13450
C                                                                         I13460
      END                                                                 I13470
C
C     --------------------------------------------------------------
C
      SUBROUTINE PANLSC (R1,JFILE,SUMR,NPTS)                              I13480
C                                                                         I13490
      IMPLICIT REAL*8          (V)                                     ! I13500
C                                                                         I13510
C     SUBROUTINE PANLSC OUTPUTS THE RESULTS OF THE SCANNING FUNCTION      I13520
C     TO FILE JFILE                                                       I13530
C                                                                         I13540
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I13550
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I13560
     *              NLTEFL,LNFIL4,LNGTH4                                  I13570
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   I13580
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       I13590
     *               DVSC,XDUM,V1SHFT                                     I13600
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 I13610
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            I13620
      COMMON /SPANEL/ V1P,V2P,DV,NLIM                                     I13630
      DIMENSION PNLHDR(2)                                                 I13640
      DIMENSION R1(*),SUMR(*)                                             I13650
C                                                                         I13660
      EQUIVALENCE (PNLHDR(1),V1P)                                         I13670
C                                                                         I13680
      CALL CPUTIM (TIME0)                                                 I13690
C                                                                         I13700
      SUMOUT = SUMR(1)                                                    I13710
      SMIN = SUMR(2)                                                      I13720
      SMAX = SUMR(3)                                                      I13730
      DV = DVO                                                            I13740
      ISTOP = 0                                                           I13750
      NNHI = (V2-VFT)/DV+1.5                                              I13760
      IF (NHI.GE.NNHI) ISTOP = 1                                          I13770
      IF (ISTOP.EQ.1) NHI = NNHI                                          I13780
      NLIM = NHI-NLO+1                                                    I13790
      V1P = VFT+FLOAT(NLO-1)*DV                                           I13800
      V2P = VFT+FLOAT(NHI-1)*DV                                           I13810
C                                                                         I13820
C     V1P IS FIRST FREQ OF PANEL                                          I13830
C     V2P IS LAST  FREQ OF PANEL                                          I13840
C                                                                         I13850
      CALL BUFOUT (JFILE,PNLHDR(1),NPHDRF)                                I13860
      CALL BUFOUT (JFILE,R1(NLO),NLIM)                                    I13870
      VFT = VFT+FLOAT(NLIMF-1)*DV                                         I13880
      IF (NPTS.GT.0) THEN                                                 I13890
         WRITE (IPR,900) V1P,V2P,DVO,NLIM                                 I13900
         WRITE (IPR,905)                                                  I13910
         NNPTS = NPTS                                                     I13920
         IF (NPTS.GT.(NLIM/2)+1) NNPTS = NLIM/2+1                         I13930
         IJLIM = NLIM-NNPTS+1                                             I13940
         DO 10 IJ = 1, NNPTS                                              I13950
            IK = IJ+IJLIM-1                                               I13960
            VI = V1P+FLOAT(IJ-1)*DVO                                      I13970
            VK = V1P+FLOAT(IK-1)*DVO                                      I13980
            JJ = NLO+IJ-1                                                 I13990
            KK = NLO+IK-1                                                 I14000
            WRITE (IPR,910) IJ,VI,R1(JJ),IK,VK,R1(KK)                     I14010
   10    CONTINUE                                                         I14020
      ENDIF                                                               I14030
      NLIMHI = NLIM+NLO-1                                                 I14040
      DO 20 I = NLO, NLIMHI                                               I14050
         SMIN = MIN(SMIN,R1(I))                                           I14060
         SMAX = MAX(SMAX,R1(I))                                           I14070
         SUMOUT = SUMOUT+R1(I)                                            I14080
   20 CONTINUE                                                            I14090
      IF (ISTOP.EQ.1) GO TO 50                                            I14100
      DO 30 J = NLIMF, MAXF                                               I14120
         R1(J-NLIMF+1) = R1(J)                                            I14130
   30 CONTINUE                                                            I14150
      DO 40 J = MAXF-NLIMF+2, MAXF                                        I14160
         R1(J) = 0.                                                       I14170
   40 CONTINUE                                                            I14180
      NLO = NSHIFT+1                                                      I14190
   50 SUMR(1) = SUMOUT                                                    I14200
      SUMR(2) = SMIN                                                      I14210
      SUMR(3) = SMAX                                                      I14220
      SUMR(4) = DVO                                                       I14230
      CALL CPUTIM (TIME)                                                  I14240
      TIMPNL = TIMPNL+TIME-TIME0                                          I14250
C                                                                         I14260
      RETURN                                                              I14270
C                                                                         I14280
  900 FORMAT ('0 V1P =',F12.5,' V2P =',F12.5,' DVOUT =',F12.8,' NLIM ='   I14290
     *   ,I10)                                                            I14300
  905 FORMAT ('0')                                                        I14310
  910 FORMAT (I5,0PF12.5,1PE12.5,I15,0PF12.5,1PE12.5)                     I14320
C                                                                         I14330
      END                                                                 I14340
C
C     --------------------------------------------------------------
C
      SUBROUTINE PNLRCT (R1,JFILE,SUMR,NPTS)                              I14350
      IMPLICIT REAL*8          (V)                                       I14360
C                                                                         I14370
C     SUBROUTINE PNLRCT OUTPUTS THE RESULTS OF THE SCANNING FUNCTION      I14380
C     TO FILE JFILE                                                       I14390
C                                                                         I14400
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I14410
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I14420
     *              NLTEFL,LBL4FL,LNGTH4                                  I14430
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   I14440
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       I14450
     *               DVSC,XDUM,V1SHFT                                     I14460
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 I14470
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            I14480
      COMMON /SPANEL/ V1P,V2P,DV,NLIM                                     I14490
      DIMENSION PNLHDR(2)                                                 I14500
      DIMENSION R1(*),SUMR(*)                                             I14510
      EQUIVALENCE (PNLHDR(1),V1P)                                         I14520
C                                                                         I14530
      CALL CPUTIM (TIME0)                                                 I14540
C                                                                         I14550
      SUMOUT = SUMR(1)                                                    I14560
      SMIN = SUMR(2)                                                      I14570
      SMAX = SUMR(3)                                                      I14580
      DV = DVO                                                            I14590
      ISTOP = 0                                                           I14600
      NNHI = (V2-VFT)/DV+1.5                                              I14610
      IF (NHI.GE.NNHI) ISTOP = 1                                          I14620
      IF (ISTOP.EQ.1) NHI = NNHI                                          I14630
      NLIM = NHI-NLO+1                                                    I14640
C                                                                         I14650
      V1P = VFT+FLOAT(NLO-1)*DV                                           I14660
      V2P = VFT+FLOAT(NHI-1)*DV                                           I14670
C                                                                         I14680
C     V1P IS FIRST FREQ OF PANEL                                          I14690
C     V2P IS LAST  FREQ OF PANEL                                          I14700
C                                                                         I14710
      CALL BUFOUT (JFILE,PNLHDR(1),NPHDRF)                                I14720
      CALL BUFOUT (JFILE,R1(NLO),NLIM)                                    I14730
      VFT = VFT+FLOAT(NLIMF-1)*DV                                         I14740
      IF (NPTS.GT.0) THEN                                                 I14750
         WRITE (IPR,900) V1P,V2P,DVO,NLIM                                 I14760
         WRITE (IPR,905)                                                  I14770
         NNPTS = NPTS                                                     I14780
         IF (NPTS.GT.(NLIM/2)+1) NNPTS = NLIM/2+1                         I14790
         IJLIM = NLIM-NNPTS+1                                             I14800
         DO 10 IJ = 1, NNPTS                                              I14810
            IK = IJ+IJLIM-1                                               I14820
            VI = V1P+FLOAT(IJ-1)*DVO                                      I14830
            VK = V1P+FLOAT(IK-1)*DVO                                      I14840
            JJ = NLO+IJ-1                                                 I14850
            KK = NLO+IK-1                                                 I14860
            WRITE (IPR,910) IJ,VI,R1(JJ),IK,VK,R1(KK)                     I14870
   10    CONTINUE                                                         I14880
      ENDIF                                                               I14890
      NLIMHI = NLIM+NLO-1                                                 I14900
      DO 20 I = NLO, NLIMHI                                               I14910
         SMIN = MIN(SMIN,R1(I))                                           I14920
         SMAX = MAX(SMAX,R1(I))                                           I14930
         SUMOUT = SUMOUT+R1(I)                                            I14940
   20 CONTINUE                                                            I14950
C                                                                         I14960
      IF (ISTOP.EQ.1) GO TO 50                                            I14970
      DO 30 J = NLIMF, MAXF
         R1(J-NLIMF+1) = R1(J)
   30 CONTINUE
      DO 40 J = MAXF-NLIMF+2, MAXF
         R1(J) = 0.
   40 CONTINUE
      NLO = NSHIFT+1
   50 SUMR(1) = SUMOUT                                                    I14980
      IPANEL = -1                                                         I14990
      SUMR(2) = SMIN                                                      I15000
      SUMR(3) = SMAX                                                      I15010
      SUMR(4) = DVO                                                       I15020
      CALL CPUTIM (TIME)                                                  I15030
      TIMPNL = TIMPNL+TIME-TIME0                                          I15040
      RETURN                                                              I15050
C                                                                         I15060
  900    FORMAT('0 V1P =',F12.5,' V2P =',F12.5,' DVOUT =',F12.8,          I15070
     *   ' NLIM =',I10)                                                   I15080
  905    FORMAT('0')
  910    FORMAT(I5,0PF12.5,1PE12.5,I15,0PF12.5,1PE12.5)                   I15100
C                                                                         I15110
      END                                                                 I15120
C
C     --------------------------------------------------------------
C
      SUBROUTINE CNVRCT (S,HWHM,R1,XF)                                    I15130
      IMPLICIT REAL*8          (V)                                       I15140
C                                                                         I15150
C     SUBROUTINE CNVRCT PERFORMS THE CONVOLUTION WITH AN ALTERNATE        I15160
C     RECTANGULAR SCANNING FUNCTION (ADJACENT BOXES OF ONE SIZE,          I15170
C     EQUAL TO 2*HWHM)
C                                                                         I15180
C     THE CONVOLUTION IS A WEIGHTED SUM THAT PROPERLY WEIGHS THE INPUT 
C     POINTS WITH THE FRACTION OF THAT POINT THAT COMPLETELY FALLS WITHIN 
C     THE OUTPUT BOX.  OUTPUT RADIANCE IS THE SUMMED RADIANCE DIVIDED
C     BY THE SUM OF THE WEIGHTS.
C
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   I15190
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       I15200
     *               DVSC,XDUM,V1SHFT                                     I15210
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I15220
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I15230
     *              NLTEFL,LBL4FL,LNGTH4                                  I15240
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 I15250
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            I15260
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      I15270
      COMMON /RCTSV/ JN,SUMJ,JFLG,RNJ,NB,IPC,VLFT,VCNT,VRGT,WGTL,WGTR     I15280
      DIMENSION S(*),R1(*),XF(*)                                          I15290
C                                                                         I15300
C     LBLRTM flags                                                        I15310
C     JFLG = -1:  first time through; increment NB when box is full       I15320
C     JFLG =  0:  subsequent calls:increment NB when box full             I15330
C     JFLG =  1:  out of data, return for more; do not increment NB       I15340
C     IDATA =-1:  first time through; or, need more data                  I15350
C     IDATA = 0:  data present                                            I15360
C     IDATA = 1:  no data present                                         I15370
C     IPANEL=-1:  first time through; or, after panel written             I15380
C     IPANEL= 0:  panel not full                                          I15390
C     IPANEL= 1:  panel is full                                           I15400
C                                                                         I15410
      CALL CPUTIM (TIME0)                                                 I15420
      RATIO = DVO/DVI                                                     I15430
C                                                                         I15440
C     During first call or if entering after writing a panel,             I15450
C     initialize: SUMJ (radiance sum), and
C                 RNJ (accumulator for number of input points in 
C                      current box, i.e. the sum of the weights)
C                 JN (box counter from 1 at VFT)
C     During first call only,
C     initialize: NB (box counter from 1 at V1), 
C                 IPC (output panel counter).
C                                                                         I15470
      IF (IPANEL.EQ.-1) THEN                                              I15480
         SUMJ = 0.                                                        I15500
         RNJ = 0.                                                         I15510
         JN = NLO                                                         I15490
         IF (JFLG.EQ.-1) THEN
            NB = 1
            IPC = 1
         ENDIF
      ENDIF                                                               I15520
C                                                                         I15530
C     Check that number of points in current panel, NNI, is correct.
C
      NNIV2 = (V2I-V1I)/DVI+1.0001                                        I15540
      IF (NNI.GT.NNIV2) NNI = NNIV2                                       I15550
C                                                                         I15560
C     Top of loop over NB boxes                                           I15600
C                                                                         I15610
   10 IF (NLO.LE.NHI) THEN                                                I15620

         VCNT = V1+(NB-1)*DVO
         IF (VCNT.GT.V2) THEN
            IPANEL = 1
            RETURN
         ENDIF
C
         VLFT = VCNT-HWHM
         VRGT = VCNT+HWHM
C
C        Find lbl panel indices for points which fall within current
C        box.
C
         RL = (VLFT-V1I)/DVI+1
         RR = (VRGT-V1I)/DVI+1
C
         IL = INT(RL+0.5)
         IH = INT(RR+0.5)
C
C        Calculate weight for each end point, inner points weighted
C        as 1.  NEP is the number of endpoints in use.
C
         VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
         WGTL = (VLBLR-VLFT)/DVI
         VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
         WGTR = (VRGT-VLBLL)/DVI
         NEP = 2
C                                                                         I15670
C        Set flag if last data point on current input panel reached       I15680
C                                                                         I15690
         IF (IH.GT.NNI) THEN                                              I15700
            IH = NNI                                                      I15710
            JFLG = 1                                                      I15720
C
C        If retrieving next panel while box sum is in progress, then
C        use weight of 1. for temp. right endpoint at IH = NNI = 2400
C        calculate partial sum below, then return.  If only one point
C        is included in this sum (IL = IH = NNI), use weight of 0
C        for right point, and add only left endpoint to sum.
            VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
            WGTR = 1.
            IF (IL.EQ.IH) THEN 
               WGTR = 0.
               NEP = 1
            ENDIF
         ENDIF                                                            I15730
C
C        If returning with new panel to partially summed box, then set
C        weight for temporary left endpoint to 1.  If it's the last
C        point going into the box, then count it as final right endpoint,
C        and use weight of 0 for left point (since IL = IH = 1).
C          
         IF (IL.LE.1) THEN
            IL = 1
            VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
            WGTL = 1.
            IF (IL.GE.IH) THEN 
               WGTL = 0.
               NEP = 1
            ENDIF
         ENDIF
C                                                                         I15740
C        If retrieving next panel while box sum is not progress, then
C        check that left edge of current output box is beyond last data on
C        panel (IL.GT.NNI), if so, go back for new panel without summing
C
         IF (JFLG.EQ.1.AND.IDATA.EQ.0.AND.IL.GT.NNI) THEN
            IPANEL = 0                                                    I15910
            JFLG = 0                                                      I15920
            RETURN                                                        I15930
         ENDIF
C                                                                         I15740
C        If last point on current input panel is reached, and there is    I15750
C        no more data to retrieve, then return                            I15760
C                                                                         I15770
         IF (JFLG.EQ.1.AND.IDATA.EQ.1) RETURN                             I15780
C                                                                         I15790
C        Compute sum for current box number NB, for all points but        I15800
C        the end points
C                                                                         I15810
         DO 20 I = IL+1, IH-1                                             I15820
            SUMJ = SUMJ+S(I)                                              I15830
   20    CONTINUE                                                         I15840
C
C        Add weighted end points to sum
         SUMJ = SUMJ+S(IL)*WGTL
         SUMJ = SUMJ+S(IH)*WGTR
C
C        Define sum of the weights, where all inner points are weighted
C        as 1, and the end points are weighted with the fraction that
C        occurs within the box.
C
         RNJ = RNJ+(IH-IL+1-NEP)+WGTL+WGTR                                I15850
C                                                                         I15860
C        If out of data on current input panel, go back for more;         I15870
C        partial SUMJ, current NB, and JFLG are saved in COMMON RCTSV     I15880
C                                                                         I15890
         IF (JFLG.EQ.1.AND.IDATA.EQ.0) THEN                               I15900
            IPANEL = 0                                                    I15910
            JFLG = 0                                                      I15920
            RETURN                                                        I15930
         ENDIF                                                            I15940
C                                                                         I15950
C        IPANEL=IDATA                                                     I15960
C                                                                         I15970
         SUMIN = SUMIN+SUMJ                                               I15980
C                                                                         I16000
C        Compute average radiance for completed box                       I16010
C                                                                         I16020
         R1(JN) = SUMJ/RNJ                                                I16030
C
         ILPR = IH+1                                                      I16040
C
C        Increment current box counters, initialize SUMJ and RNJ
         JN = JN+1                                                        I16050
         SUMJ = 0.                                                        I16050
         RNJ = 0.                                                         I16050
C
C        Output panel when number of boxes, NB, reaches a multiple of     I16060
C        2400, using then incrementing current output panel number, IPC.
C                                                                         I16070
         IF (NB.EQ.IPC*(NHI-NLO+1)) THEN                                  I16080
            IPANEL = 1                                                    I16090
            IPC = IPC+1
            NB = NB+1
            RETURN                                                        I16100
         ENDIF                                                            I16110
C                                                                         I16150
C        Increment NB
         NB = NB+1

C        Go back to top of loop over NB boxes                             I16160
C                                                                         I16170
         GO TO 10                                                         I16180
      ENDIF                                                               I16190
C                                                                         I16200
      RETURN                                                              I16210
      END                                                                 I16220
C
C     --------------------------------------------------------------
C
      SUBROUTINE CNVVRC (S,AFOV,R1,XF)                                    I15130
      IMPLICIT REAL*8          (V)                                       I15140
C                                                                         I15150
C     SUBROUTINE CNVVRC PERFORMS THE CONVOLUTION WITH A RECTANGULAR       I15160
C     SCANNING FUNCTION OF VARIABLE SIZE, WHERE THE BOX SIZE IS           I15170
C     WAVENUMBER DEPENDENT.  V1, V2, AND DVO ARE USED TO DEFINE THE
C     CENTER OF THE OUTPUT BOXES.  BOXES OVERLAP WHERE NECESSARY TO
C     INSURE A CONSTANT DVO.
C
C     BFOV is used to determine the resolution (box size), which is
C     spectrally variable.
C
C     AFOV is passed in from calls to CNVVRC from SCANFN and SCNMRG
C     as HWHM, since the value of HWHM on Record 8.1 on TAPE5 holds
C     the place of the half angle of the instrument field of view in
C     degrees.
C    
C     BOX WIDTH EQUALS V*B**2/2, AND THE SHIFT EQUALS HALF THE BOX WIDTH
C     V*(1-B**2/4), WHERE B IS THE HALF ANGLE OF THE INSTRUMENT FIELD
C     OF VIEW IN RADIANS.
C                                                                         I15180
C     THE CONVOLUTION IS A WEIGHTED SUM THAT PROPERLY WEIGHS THE INPUT 
C     POINTS WITH THE FRACTION OF THAT POINT THAT COMPLETELY FALLS WITHIN 
C     THE OUTPUT BOX.  OUTPUT RADIANCE IS THE SUMMED RADIANCE DIVIDED
C     BY THE SUM OF THE WEIGHTS.
C
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   I15190
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       I15200
     *               DVSC,XDUM,V1SHFT                                     I15210
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I15220
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I15230
     *              NLTEFL,LBL4FL,LNGTH4                                  I15240
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 I15250
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            I15260
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      I15270
      COMMON /RCTSV/ JN,SUMJ,JFLG,RNJ,NB,IPC,VLFT,VCNT,VRGT,WGTL,WGTR     I15280
      DIMENSION S(*),R1(*),XF(*)                                          I15290
C                                                                         I15300
C     LBLRTM flags                                                        I15310
C     JFLG = -1:  first time through; increment NB when box is full       I15320
C     JFLG =  0:  subsequent calls:increment NB when box full             I15330
C     JFLG =  1:  out of data, return for more; do not increment NB       I15340
C     IDATA =-1:  first time through; or, need more data                  I15350
C     IDATA = 0:  data present                                            I15360
C     IDATA = 1:  no data present                                         I15370
C     IPANEL=-1:  first time through; or, after panel written             I15380
C     IPANEL= 0:  panel not full                                          I15390
C     IPANEL= 1:  panel is full                                           I15400
C                                                                         I15410
      CALL CPUTIM (TIME0)                                                 I15420
C
C     Convert AFOV to BFOV, the half angle of the instrument field of
C     view in radians. (For IRIS-D, AFOV equals 2.5 degrees)

      BFOV = AFOV*3.141592654/180.
C
      RATIO = DVO/DVI                                                     I15430
C                                                                         I15440
C     During first call or if entering after writing a panel,             I15450
C     initialize: SUMJ (radiance sum), and
C                 RNJ (accumulator for number of input points in 
C                      current box, i.e. the sum of the weights)
C                 JN (box counter from 1 at VFT)
C     During first call only,
C     initialize: NB (box counter from 1 at V1), 
C                 IPC (output panel counter).
C                                                                         I15470
      IF (IPANEL.EQ.-1) THEN                                              I15480
         SUMJ = 0.                                                        I15500
         RNJ = 0.                                                         I15510
         JN = NLO                                                         I15490
         IF (JFLG.EQ.-1) THEN
            NB = 1
            IPC = 1
         ENDIF
      ENDIF                                                               I15520
C                                                                         I15530
C     Check that number of points in current panel, NNI, is correct.
C
      NNIV2 = (V2I-V1I)/DVI+1.0001                                        I15540
      IF (NNI.GT.NNIV2) NNI = NNIV2                                       I15550
C
C     Top of loop over NB boxes                                           I15600
C                                                                         I15610
   10 IF (NLO.LE.NHI) THEN                                                I15620

C     For current box find wavenumber at center and left/right edges.
C     For first box, VCNT equals V1.  When current box exceeds V2,
C     then exit.

         VCNT = V1+(NB-1)*DVO
         IF (VCNT.GT.V2) THEN
            IPANEL = 1
            RETURN
         ENDIF

         VLFT = VCNT*(1-BFOV**2/4)
         VRGT = VCNT*(1+BFOV**2/4)

C     Find lbl panel indices for points which fall within current box.

         RL = (VLFT-V1I)/DVI+1
         RR = (VRGT-V1I)/DVI+1
C
         IL = INT(RL+0.5)
         IH = INT(RR+0.5)
C
C     Calculate weight for each end point, inner points weighted as 1.
C     NEP is the number of endpoints in use.
C
         VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
         WGTL = (VLBLR-VLFT)/DVI
         VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
         WGTR = (VRGT-VLBLL)/DVI
         NEP = 2
C
C        Set flag if last data point on current input panel reached       I15680
C                                                                         I15690
         IF (IH.GT.NNI) THEN                                              I15700
            IH = NNI                                                      I15710
            JFLG = 1                                                      I15720
C
C        If retrieving next panel while box sum is in progress, then
C        use weight of 1. for temp. right endpoint at IH = NNI = 2400
C        calculate partial sum below, then return.  If only one point
C        is included in this sum (IL = IH = NNI), use weight of 0
C        for right point, and add only left endpoint to sum.
            VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
            WGTR = 1.
            IF (IL.EQ.IH) THEN 
               WGTR = 0.
               NEP = 1
            ENDIF
         ENDIF                                                            I15730
C
C        If returning with new panel to partially summed box, then set
C        weight for temporary left endpoint to 1.  If it's the last
C        point going into the box, then count it as final right endpoint,
C        and use weight of 0 for left point (since IL = IH = 1).
C
         IF (IL.LE.1) THEN
            IL = 1
            VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
            WGTL = 1.
            IF (IL.GE.IH) THEN 
               WGTL = 0.
               NEP = 1
            ENDIF
         ENDIF
C                                                                         I15740
C        If retrieving next panel while box sum is not progress, then
C        check that left edge of current output box is beyond last data on
C        panel (IL.GT.NNI), if so, go back for new panel without summing
C
         IF (JFLG.EQ.1.AND.IDATA.EQ.0.AND.IL.GT.NNI) THEN
            IPANEL = 0                                                    I15910
            JFLG = 0                                                      I15920
            RETURN                                                        I15930
         ENDIF
C
C        If last point on current input panel is reached, and there is    I15750
C        no more data to retrieve, then return                            I15760
C                                                                         I15770
         IF (JFLG.EQ.1.AND.IDATA.EQ.1) RETURN                             I15780
C                                                                         I15790
C        Compute sum for current box number NB, for all points but        I15800
C        the end points
C                                                                         I15810
         DO 20 I = IL+1, IH-1                                             I15820
            SUMJ = SUMJ+S(I)                                              I15830
   20    CONTINUE                                                         I15840
C
C        Add weighted end points to sum
         SUMJ = SUMJ+S(IL)*WGTL
         SUMJ = SUMJ+S(IH)*WGTR
C
C        Define sum of the weights, where all inner points are weighted
C        as 1, and the end points are weighted with the fraction that
C        occurs within the box.
C
         RNJ = RNJ+(IH-IL+1-NEP)+WGTL+WGTR                                I15850
C                                                                         I15860
C        If out of data on current input panel, go back for more;         I15870
C        partial SUMJ, current NB, and JFLG are saved in COMMON RCTSV     I15880
C                                                                         I15890
         IF (JFLG.EQ.1.AND.IDATA.EQ.0) THEN                               I15900
            IPANEL = 0                                                    I15910
            JFLG = 0                                                      I15920
            RETURN                                                        I15930
         ENDIF                                                            I15940
C                                                                         I15950
C        IPANEL=IDATA                                                     I15960
C                                                                         I15970
         SUMIN = SUMIN+SUMJ                                               I15980
C                                                                         I15990
C                                                                         I16000
C        Compute average radiance for completed box                       I16010
C                                                                         I16020
         R1(JN) = SUMJ/RNJ                                                I16030
C                                                                         I16050
C        Increment current box counters, initialize SUMJ and RNJ
         JN = JN+1                                                        I16050
         SUMJ = 0.                                                        I16050
         RNJ = 0.                                                         I16050
C
C        Output panel when number of boxes, NB, reaches a multiple of     I16060
C        2400, using then incrementing current output panel number, IPC.
C                                                                         I16070
         IF (NB.EQ.IPC*(NHI-NLO+1)) THEN                                  I16080
            IPANEL = 1                                                    I16090
            IPC = IPC+1
            NB = NB+1
            RETURN                                                        I16100
         ENDIF                                                            I16110

C        Increment NB
         NB = NB+1
C                                                                         I16150
C        Go back to top of loop over NB boxes                             I16160
C                                                                         I16170
         GO TO 10                                                         I16180

      ENDIF                                                               I16190
C                                                                         I16200
      RETURN                                                              I16210
      END                                                                 I16220
C
C     --------------------------------------------------------------
C
      SUBROUTINE CNVVRL (S,AFOV,R1,XF)                                    I15130
      IMPLICIT REAL*8          (V)                                       I15140
C                                                                         I15150
C     SUBROUTINE CNVVRL PERFORMS THE CONVOLUTION WITH A RECTANGULAR       I15160
C     SCANNING FUNCTION OF VARIABLE SIZE, WHERE THE BOX SIZE IS           I15170
C     WAVENUMBER DEPENDANT.  V1, V2, AND DVO ARE USED TO DEFINE THE
C     LEFT EDGE OF THE OUTPUT BOXES.  BOXES OVERLAP WHERE NECESSARY
C     TO INSURE A CONSTANT DVO.
C
C     BFOV is used to determine the resolution (box size), which is
C     spectrally variable.
C    
C     AFOV is passed in from calls to CNVVRL from SCANFN and SCNMRG
C     as HWHM, since the value of HWHM on Record 8.1 on TAPE5 holds
C     the place of the half angle of the instrument field of view in
C     degrees.
C    
C     BOX WIDTH EQUALS V*B**2/2, AND THE SHIFT EQUALS HALF THE BOX WIDTH
C     V*(1-B**2/4), WHERE B IS THE HALF ANGLE OF THE INSTRUMENT FIELD
C     OF VIEW IN RADIANS.
C                                                                         I15180
C     THE CONVOLUTION IS A WEIGHTED SUM THAT PROPERLY WEIGHS THE INPUT 
C     POINTS WITH THE FRACTION OF THAT POINT THAT COMPLETELY FALLS WITHIN 
C     THE OUTPUT BOX.  OUTPUT RADIANCE IS THE SUMMED RADIANCE DIVIDED
C     BY THE SUM OF THE WEIGHTS.
C
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   I15190
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       I15200
     *               DVSC,XDUM,V1SHFT                                     I15210
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I15220
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I15230
     *              NLTEFL,LBL4FL,LNGTH4                                  I15240
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 I15250
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            I15260
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      I15270
      COMMON /RCTSV/ JN,SUMJ,JFLG,RNJ,NB,IPC,VLFT,VCNT,VRGT,WGTL,WGTR     I15280
      DIMENSION S(*),R1(*),XF(*)                                          I15290
C                                                                         I15300
C     LBLRTM flags                                                        I15310
C     JFLG = -1:  first time through; increment NB when box is full       I15320
C     JFLG =  0:  subsequent calls:increment NB when box full             I15330
C     JFLG =  1:  out of data, return for more; do not increment NB       I15340
C     IDATA =-1:  first time through; or, need more data                  I15350
C     IDATA = 0:  data present                                            I15360
C     IDATA = 1:  no data present                                         I15370
C     IPANEL=-1:  first time through; or, after panel written             I15380
C     IPANEL= 0:  panel not full                                          I15390
C     IPANEL= 1:  panel is full                                           I15400
C                                                                         I15410
      CALL CPUTIM (TIME0)                                                 I15420
C
C     Convert AFOV to BFOV, the half angle of the instrument field of
C     view in radians. (For IRIS-D, AFOV equals 2.5 degrees)

      BFOV = AFOV*3.141592654/180.

      RATIO = DVO/DVI                                                     I15430
C                                                                         I15440
C     During first call or if entering after writing a panel,             I15450
C     initialize: SUMJ (radiance sum), and
C                 RNJ (accumulator for number of input points in 
C                      current box, i.e. the sum of the weights)
C                 JN (box counter from 1 at VFT)
C     During first call only,
C     initialize: NB (box counter from 1 at V1), 
C                 IPC (output panel counter).
C                                                                         I15470
      IF (IPANEL.EQ.-1) THEN                                              I15480
         SUMJ = 0.                                                        I15500
         RNJ = 0.                                                         I15510
         JN = NLO                                                         I15490
         IF (JFLG.EQ.-1) THEN
            NB = 1
            IPC = 1
         ENDIF
      ENDIF                                                               I15520
C                                                                         I15530
C     Check that number of points in current panel, NNI, is correct.
C
      NNIV2 = (V2I-V1I)/DVI+1.0001                                        I15540
      IF (NNI.GT.NNIV2) NNI = NNIV2                                       I15550
C                                                                         I15560
C     Top of loop over NB boxes                                           I15600
C                                                                         I15610
   10 IF (NLO.LE.NHI) THEN                                                I15620
C
C     For current box find wavenumber at the left and right edges.
C     For first box, VLFT equals V1.  When current box exceeds V2,
C     then exit.
C
         VLFT = V1+(NB-1)*DVO
         IF (VLFT.GT.V2) THEN
            IPANEL = 1
            RETURN
         ENDIF
C
         VCNT = VLFT*(1+BFOV**2/4)
         VRGT = VLFT*(1+BFOV**2/2)
C
C     Find lbl panel indices for points which fall within current box.
C
         RL = (VLFT-V1I)/DVI+1
         RR = (VRGT-V1I)/DVI+1
C
         IL = INT(RL+0.5)
         IH = INT(RR+0.5)
C
C     Calculate weight for each end point, inner points weighted as 1,
C     NEP is the number of endpoints in use.
C
         VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
         WGTL = (VLBLR-VLFT)/DVI
         VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
         WGTR = (VRGT-VLBLL)/DVI
         NEP = 2
C
C        Set flag if last data point on current input panel reached       I15680
C                                                                         I15690
         IF (IH.GT.NNI) THEN                                              I15700
            IH = NNI                                                      I15710
            JFLG = 1                                                      I15720
C
C        If retrieving next panel while box sum is in progress, then
C        use weight of 1. for temp. right endpoint at IH = NNI, 
C        calculate partial sum below, then return.  If only one point
C        is included in this sum (IL = IH = NNI), use weight of 0
C        for right point, and add only left endpoint to sum.
            VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
            WGTR = 1.
            IF (IL.EQ.IH) THEN 
               WGTR = 0.
               NEP = 1
            ENDIF
         ENDIF                                                            I15730
C
C        If returning with new panel to partially summed box, then set
C        weight for temporary left endpoint to 1.  If it's the last
C        point going into the box, then count it as final right endpoint,
C        and use weight of 0 for left point (since IL = IH = 1).
C
         IF (IL.LE.1) THEN
            IL = 1
            VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
            WGTL = 1.
            IF (IL.GE.IH) THEN 
               WGTL = 0.
               NEP = 1
            ENDIF
         ENDIF
C                                                                         I15740
C        If retrieving next panel while box sum is not in progress, then
C        check that left edge of current output box is beyond last data on
C        panel (IL.GT.NNI), if so, go back for new panel without summing
C
         IF (JFLG.EQ.1.AND.IDATA.EQ.0.AND.IL.GT.NNI) THEN
            IPANEL = 0                                                    I15910
            JFLG = 0                                                      I15920
            RETURN                                                        I15930
         ENDIF
C
C        If last point on current input panel is reached, and there is    I15750
C        no more data to retrieve, then return                            I15760
C                                                                         I15770
         IF (JFLG.EQ.1.AND.IDATA.EQ.1) RETURN                             I15780
C                                                                         I15790
C        Compute sum for current box number NB, using a weight of 1.0     I15800
C        for all points but the end points, which use WGTL and WGTR
C                                                                         I15810
         DO 20 I = IL+1, IH-1                                             I15820
            SUMJ = SUMJ+S(I)                                              I15830
   20    CONTINUE                                                         I15840
C
C        Add weighted end points to sum
         SUMJ = SUMJ+S(IL)*WGTL
         SUMJ = SUMJ+S(IH)*WGTR
C
C        Define sum of the weights, where all inner points are weighted
C        as 1, and the end points are weighted with the fraction that
C        occurs within the box.
C
         RNJ = RNJ+(IH-IL+1-NEP)+WGTL+WGTR                                I15850
C                                                                         I15860
C        If out of data on current input panel, go back for more;         I15870
C        partial SUMJ, current NB, and JFLG are saved in COMMON RCTSV     I15880
C                                                                         I15890
         IF (JFLG.EQ.1.AND.IDATA.EQ.0) THEN                               I15900
            IPANEL = 0                                                    I15910
            JFLG = 0                                                      I15920
            RETURN                                                        I15930
         ENDIF                                                            I15940
C                                                                         I15950
C        IPANEL=IDATA                                                     I15960
C                                                                         I15970
         SUMIN = SUMIN+SUMJ                                               I15980
C                                                                         I16000
C        Compute average radiance for completed box                       I16010
C                                                                         I16020
         R1(JN) = SUMJ/RNJ                                                I16030
C
C        Increment current box counters, initialize SUMJ and RNJ
         JN = JN+1                                                        I16040
         SUMJ = 0.                                                        I16040
         RNJ = 0.                                                         I16040
C                                                                         I16050
C        Output panel when number of boxes, NB, reaches a multiple of     I16060
C        2400, using then incrementing current output panel number, IPC.
C                                                                         I16070
         IF (NB.EQ.IPC*(NHI-NLO+1)) THEN                                  I16080
            IPANEL = 1                                                    I16090
            IPC = IPC+1
            NB = NB+1
            RETURN                                                        I16100
         ENDIF                                                            I16110
C
C        Increment NB
         NB = NB+1
C                                                                         I16150
C        Go back to top of loop over NB boxes                             I16160
C                                                                         I16170
         GO TO 10                                                         I16180
C
      ENDIF                                                               I16190
C                                                                         I16200
      RETURN                                                              I16210
      END                                                                 I16220
C
C     --------------------------------------------------------------
C
      SUBROUTINE SINCSQ (XF,XSCALE)                                       I16230
C                                                                         I16240
C     SUBROUTINE SINCSQ SETS UP THE SINCSQ SCANNING FUNCTION              I16250
C                                                                         I16260
      COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,
     *                HWF3,DXF3,NX3,N3MAX
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I16290
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I16300
     *              NLTEFL,LNFIL4,LNGTH4                                  I16310
      DIMENSION XF(*)                                                     I16320
C                                                                         I16330
C     DATA XSCALE / 1.391557377 /                                         I16340
C                                                                         I16350
      XSINC2(X) = (SIN(X)/X)**2                                           I16360
      PI = 2.*ASIN(1.)                                                    I16370
C                                                                         I16380
C     PI CORRESPONDS TO X=2.257609141                                     I16390
C                                                                         I16400
      XNORM = XSCALE/PI                                                   I16410
      DO 10 I = 1, NFMAX                                                  I16420
         XF(I) = 0.                                                       I16430
   10 CONTINUE                                                            I16440
      XF(1) = XNORM                                                       I16450
      SUM = XF(1)                                                         I16460
      DO 20 I = 2, NF                                                     I16470
         X = FLOAT(I-1)*DXF                                               I16480
         XF(I) = XNORM*XSINC2(X*XSCALE)                                   I16490
         SUM = SUM+2.*XF(I)                                               I16500
   20 CONTINUE                                                            I16510
      SUM = SUM*DXF                                                       I16520
C                                                                         I16530
CPRT  WRITE(IPR,900) NF,DXF,SUM                                           I16540
C                                                                         I16550
      RETURN                                                              I16560
C                                                                         I16570
  900 FORMAT ('0',5X,'NF =',I5,',  DXF =',F7.5,',    SUM =',F18.15)       I16580
C                                                                         I16590
      END                                                                 I16600
C
C     --------------------------------------------------------------
C
      SUBROUTINE SINC (XF,XSCALE)                                         I16610
C                                                                         I16620
C     SUBROUTINE SINC SETS UP THE SINC SCANNING FUNCTION                  I16630
C                                                                         I16640
      COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,
     *                HWF3,DXF3,NX3,N3MAX
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         I16670
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        I16680
     *              NLTEFL,LNFIL4,LNGTH4                                  I16690
      DIMENSION XF(*)                                                     I16700
C                                                                         I16710
C     DATA XSCALE / 1.89549425  /                                         I16720
C                                                                         I16730
      XSINC(X) = (SIN(X)/X)                                               I16740
      PI = 2.*ASIN(1.)                                                    I16750
C                                                                         I16760
C     PI CORRESPONDS TO X=1.657400255                                     I16770
C                                                                         I16780
      XNORM = XSCALE/PI                                                   I16790
      DO 10 I = 1, NFMAX                                                  I16800
         XF(I) = 0.                                                       I16810
   10 CONTINUE                                                            I16820
      XF(1) = XNORM                                                       I16830
      SUM = XF(1)                                                         I16840
      DO 20 I = 2, NF                                                     I16850
         X = FLOAT(I-1)*DXF                                               I16860
         XF(I) = XNORM*XSINC(X*XSCALE)                                    I16870
         SUM = SUM+2.*XF(I)                                               I16880
   20 CONTINUE                                                            I16890
      SUM = SUM*DXF                                                       I16900
C                                                                         I16910
CPRT  WRITE(IPR,900) NF,DXF,SUM                                           I16920
C                                                                         I16930
      RETURN                                                              I16940
C                                                                         I16950
  900 FORMAT ('0',5X,'NF =',I5,',  DXF =',F7.5,',    SUM =',F18.15)       I16960
C                                                                         I16970
      END                                                                 I16980
C
C     --------------------------------------------------------------
C
      SUBROUTINE INTRPL (IFILE,JFILE)                                     J00010
C                                                                         J00020
      IMPLICIT REAL*8          (V)                                     ! J00030
C                                                                         J00040
C**********************************************************************   J00050
C                                                                         J00060
C     INTERPOLATION FUNCTION DRIVER: FOUR POINT VERSION                   J00070
C                                                                         J00080
C           A.E.R. INC.           (MARCH 1990)                            J00090
C                                                                         J00100
C**********************************************************************   J00110
C                                                                         J00120
C     THE INPUT DATA WILL BE PUT INTO T(5) WITH THE LAST                  J00130
C     4 POINTS OF THE PREVIOUS PANEL PUT INTO T(1 TO 4).                  J00140
C     THIS SCHEME PERMITS 6 POINT INTERPOLATION.                          J00150
C                                                                         J00160
C     S IS NOMINALLY 2401 POINTS BUT MAY NEED TO BE EXTENDED BY TWO (2)   J00170
C     POINTS TO PERMIT 4 POINT INTERPOLATION UP TO THE LAST DATA POINT.   J00180
C                                                                         J00190
      COMMON T(2410),R(2401)                                              J00200
      DIMENSION S(2406)                                                   J00210
      EQUIVALENCE (S(1),T(5))                                             J00220
C                                                                         J00230
      character*8      XID,       HMOLID,      YID        
      real*8               SECANT,       XALTZ 
C                                                                         J00250
      COMMON /HVERSN/  HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,
     *                HVROPR,HVRPST,HVRPLT,HVRTST,HVRUTL,HVRXMR
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       J00260
     *                WK(60),PZL,PZU,TZL,TZU,WN2   ,DV ,V1C,V2C,TBOUND,   J00270
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    J00280
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   J00290
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       J00300
     *               DVSC,XDUM,V1SHFT                                     J00310
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 J00320
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            J00330
      COMMON /INPNL/ V1I,V2I,DVI,NNI                                      J00340
      COMMON /OUTPNL/ V1J,V2J,DVJ,NNJ                                     J00350
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         J00360
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        J00370
     *              NLTEFL,LNFIL4,LNGTH4                                  J00380
      COMMON /FLFORM/ CFORM                                               J00390
C                                                                         J00400
      DIMENSION FILHDR(2),RSTAT(3)                                        J00410
C                                                                         J00420
      EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),                 J00430
     *            (FSCDID(6),ISCHDR) , (FSCDID(12),XSCID),                J00440
     *            (FSCDID(13),XHWHM) , (FSCDID(14),IDABS),                J00450
     *            (FSCDID(16),LAYR1)                                      J00460
C                                                                         J00470
      CHARACTER*12 BCD,HTRANS,HABSRB,HRADIA                               J00480
      CHARACTER CFORM*11,SCNOUT*6,CTAPE*4                                 J00490
      CHARACTER*8 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *            HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR
      LOGICAL OP                                                          J00500
C                                                                         J00510
      DATA HTRANS / 'TRANSMISSION'/,HABSRB / ' ABSORPTION '/,             J00520
     *     HRADIA / ' RADIANCE   '/                                       J00530
      DATA SCNOUT / '      '/,CTAPE / 'TAPE'/                             J00540
C                                                                         J00550
C----------------------------------------------------------------------   J00560
C     JEMIT=-1  INTERPOLATE ABSORPTION                                    J00570
C     JEMIT=0   INTERPOLATE TRANSMISSION                                  J00580
C     JEMIT=1   INTERPOLATE EMISSION                                      J00590
C     JEMIT=2   INTERPOLATE OPTICAL DEPTH                                 J00600
C----------------------------------------------------------------------   J00610
C                                                                         J00620
C
C     ASSIGN SCCS VERSION NUMBER TO MODULE 
C
      HVRPST = '$Revision$'
C
   10 CONTINUE                                                            J00630
      CALL CPUTIM (TIME1)                                                 J00640
      TIMRDF = 0.0                                                        J00650
      TIMCNV = 0.0                                                        J00660
      TIMPNL = 0.0                                                        J00670
C                                                                         J00680
      READ (IRD,900,END=50) DVO,V1,V2,JEMIT,I4PT,IUNIT,IFILST,NIFILS,     J00690
     *                      JUNIT,NPTS                                    J00700
      IF (DVO.LE.0.) GO TO 40                                             J00710
C                                                                         J00720
C     V2 IS ONLY APPROXIMATE                                              J00730
C                                                                         J00740
      NUM = (((V2-V1)/DVO)+0.5)                                           J00750
      V2 = V1+FLOAT(NUM)*DVO                                              J00760
      NUM = NUM+1                                                         J00770
      WRITE (IPR,905) V1,V2,DVO,NUM,JEMIT,I4PT,IUNIT,IFILST,JUNIT,NPTS    J00780
C                                                                         J00790
C     SET INPUT(IFILE), OUTPUT(JFILE) UNITS.                              J00800
C                                                                         J00810
      IF (IUNIT.LE.0) IUNIT = IFILE                                       J00820
      IFILE = IUNIT                                                       J00830
      INQUIRE (UNIT=IFILE,OPENED=OP)
      IF (.NOT.OP) THEN
         WRITE (SCNOUT,910) CTAPE,IFILE
         OPEN (IFILE,FILE=SCNOUT,STATUS='UNKNOWN',FORM=CFORM)
      ENDIF
      IFILST = MAX(IFILST,1)                                              J00840
      IF (NIFILS.LE.0) NIFILS = 99                                        J00850
      IF (JUNIT.LE.0) JUNIT = JFILE                                       J00860
      JFILE = JUNIT                                                       J00870
      INQUIRE (UNIT=JFILE,OPENED=OP)                                      J00880
      IF (.NOT.OP) THEN                                                   J00890
         WRITE (SCNOUT,910) CTAPE,JFILE                                   J00900
         OPEN (JFILE,FILE=SCNOUT,STATUS='UNKNOWN',FORM=CFORM)             J00910
         REWIND JFILE                                                     J00920
      ENDIF                                                               J00930
C                                                                         J00940
      REWIND IFILE                                                        J00950
      IF (IFILST.GT.1) CALL SKIPFL (IFILST-1,IFILE,IEOF)                  J00960
C                                                                         J00970
C     BUFFER IN THE FILE HEADER ON UNIT (IFILE)                           J00980
C     BUFFER OUT ON UNIT (JFILE)                                          J00990
C                                                                         J01000
      CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)                             J01010
      IF (IEOF.EQ.0) GO TO 10                                             J01020
C                                                                         J01030
      WRITE (IPR,915) XID,(YID(M),M=1,2)                                  J01040
      WRITE (IPR,920) LAYR1,LAYER                                         J01050
      WRITE (IPR,925) SECANT,PAVE,TAVE,DV,V1C,V2C                         J01060
      WRITE (IPR,930) WN2,(HMOLID(M),WK(M),M=1,NMOL)                      J01070
C                                                                         J01080
      JABS = 0                                                            J01090
      IDABS = 0                                                           J01100
      IF (JEMIT.LT.0) THEN                                                J01110
         JABS = 1                                                         J01120
         JEMIT = 0                                                        J01130
         IDABS = -1                                                       J01140
      ENDIF                                                               J01150
C                                                                         J01160
      ISCAN = ISCHDR                                                      J01170
      IF (ISCAN.LE.0.OR.XSCID.EQ.-99.) ISCAN = 0                          J01180
      IF (ISCHDR.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCHDR                   J01190
      ISCHDR = ISCAN+10                                                   J01200
      V1C = V1                                                            J01210
      V2C = V2                                                            J01220
      DV = DVO                                                            J01230
C                                                                         J01240
      SCNID = 100*JEMIT                                                   J01250
      XSCID = SCNID+0.01                                                  J01260
C                                                                         J01270
      CALL BUFOUT (JFILE,FILHDR(1),NFHDRF)                                J01280
C                                                                         J01290
      ICNVRT = 0                                                          J01300
      JTREM = -1                                                          J01310
      IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0                        J01320
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2                        J01330
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.2)) JTREM = 2                        J01340
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1                        J01350
      ISCANT = MOD(ISCAN,1000)                                            J01360
      IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 2                       J01370
      IF (JTREM.LT.0) STOP 71048                                          J01380
      WRITE (IPR,935) IEMIT,JEMIT,JTREM                                   J01390
      WRITE (IPR,940) IFILE,IFILST,NIFILS,JEMIT,JABS                      J01400
C                                                                         J01410
      IDATA = -1                                                          J01420
C                                                                         J01430
C     NEED TO SAVE LAST IBOUND POINTS OF EACH PANEL TO ATTACH TO NEXT     J01440
C                                                                         J01450
      IBOUND = 4                                                          J01460
C                                                                         J01470
C     VBOT IS LOWEST NEEDED WAVENUMBER, VTOP IS HIGHEST                   J01480
C                                                                         J01490
      BOUND = FLOAT(IBOUND)*DV                                            J01500
      VBOT = V1-BOUND                                                     J01510
      VTOP = V2+BOUND                                                     J01520
C                                                                         J01530
      IF (JEMIT.EQ.0.AND.IDABS.EQ.0) BCD = HTRANS                         J01540
      IF (JEMIT.EQ.0.AND.IDABS.EQ.-1) BCD = HABSRB                        J01550
      IF (JEMIT.EQ.1) BCD = HRADIA                                        J01560
      IF (NPTS.GT.0) WRITE (IPR,945) BCD                                  J01570
C                                                                         J01580
C     ZERO OUT T(1 TO IBOUND)                                             J01590
C                                                                         J01600
      DO 20 II = 1, IBOUND                                                J01610
         T(II) = 0.0                                                      J01620
   20 CONTINUE                                                            J01630
C                                                                         J01640
C     READ FROM IFILE UNTIL THE FIRST REQUIRED POINT IS REACHED           J01650
C     AND LOAD DATA INTO S                                                J01660
C                                                                         J01670
      CALL RDPANL (S,JTREM,IFILE,ISCAN,JEMIT,ICNRT)                       J01680
      IF (IEOFSC.LE.0) GO TO 30                                           J01690
C                                                                         J01700
C     DO INTERPOLATION                                                    J01710
C                                                                         J01720
      CALL INTERP (IFILE,JFILE,I4PT,IBOUND,NPTS,JTREM,ISCAN,JEMIT,        J01730
     *             RSTAT,ICNVRT)                                          J01740
C                                                                         J01750
      CALL CPUTIM (TIME2)                                                 J01760
      CALL ENDFIL (JFILE)                                                 J01770
C                                                                         J01780
C     WRITE STATISTICS                                                    J01790
C                                                                         J01800
      WRITE (IPR,950) RSTAT(1),RSTAT(2),RSTAT(3)                          J01810
      TIMTOT = TIME2-TIME1                                                J01820
      TIMCNV = TIMTOT-TIMRDF-TIMPNL                                       J01830
      WRITE (IPR,955) TIMTOT,TIMRDF,TIMCNV,TIMPNL                         J01840
      GO TO 10                                                            J01850
C                                                                         J01860
   30 CONTINUE                                                            J01870
      WRITE (IPR,960) IFILE                                               J01880
C                                                                         J01890
      GO TO 10                                                            J01900
C                                                                         J01910
   40 CONTINUE                                                            J01920
      WRITE (IPR,965)                                                     J01930
C                                                                         J01940
      RETURN                                                              J01950
C                                                                         J01960
   50 CONTINUE                                                            J01970
      WRITE (IPR,970) IRD                                                 J01980
      STOP ' INTRPL'                                                      J01990
C                                                                         J02000
  900 FORMAT (3F10.3,2I5,15X,5I5)                                         J02010
  905 FORMAT (5X,'V1=',F14.8,' V2=',F14.8,' DVO=',E14.6,' NUM=',I8,/5X,   J02020
     *   'JEMIT=',I3,' I4PT=',I3,' IUNIT=',I3,' IFILST=',I3,/5X,          J02030
     *   'JUNIT=',I3,' NPTS=',I5)                                         J02040
  910 FORMAT (A4,I2.2)                                                    J02050
  915 FORMAT (//,' ***INTRPL***',/,'0',10A8,2X,2(1X,A8,1X))               J02060
  920 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)            J02070
  925 FORMAT ('  SECANT =',F15.5,/'  PRESS(MB) =',F12.5/'  TEMP(K) =',    J02080
     *   F11.2,/'  DV(CM-1) = ',F12.8,/,'  V1(CM-1) = ',F12.6,/,          J02090
     *   '  V2(CM-1) = ',F12.6)                                           J02100
  930 FORMAT (/,'  COLUMN DENSITY (MOLECULES/CM**2)',/,5X,'WBROAD = ',    J02110
     *   1PE10.3,/(5X,A6,' = ',1PE10.3))                                  J02120
  935 FORMAT (5X,'IEMIT=',I5,' JEMIT=',I5,' JTREM=',I5)                   J02130
  940 FORMAT (5X,'INPUT FILE NUMBER =',I3,' ,IFILST = ',I3,               J02140
     *   ' ,NIFILS = ',I3,',JEMIT =',I2,' ,JABS =',I2)                    J02150
  945 FORMAT (///,'0',5X,A12,/)                                           J02160
  950 FORMAT (/,5X,'SUMOUT =',1P,E16.9,'  MIN =',E16.9,'  MAX =',E16.9)   J02170
  955 FORMAT (/,5X,'TIME: TOTAL = ',F8.3,' READ = ',F8.3,' INTERP = ',    J02180
     *   F8.3,' WRITE = ',F8.3)                                           J02190
  960 FORMAT (/,5X,'INTRPL- ERROR: EOF ON INPUT UNIT ',I4,                J02200
     *   ' BEFORE V1 WAS REACHED')                                        J02210
  965 FORMAT (/,5X,'END OF INTERPOLATION REQUESTS')                       J02220
  970 FORMAT (/,5X,' INTRP - ERROR: EOF ON STANDARD INPUT, UNIT = ',I4)   J02230
C                                                                         J02240
      END                                                                 J02250
C
C     --------------------------------------------------------------
C
      SUBROUTINE RDPANL (S,JTREM,IFILE,ISCAN,JEMIT,ICNVRT)                J02260
C                                                                         J02270
      IMPLICIT REAL*8          (V)                                     ! J02280
C                                                                         J02290
C     SUBROUTINE RDPANL INPUTS PANELS FROM IFILE RESULTING FROM THE       J02300
C     LBLRTM CALCULATION                                                  J02310
C                                                                         J02320
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           J02330
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   J02340
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       J02350
     *               DVSC,XDUM,V1SHFT                                     J02360
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            J02370
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 J02380
      COMMON /INPNL/ VMIN,VMAX,DVI,NNI                                    J02390
      COMMON /RPANL/ V1P,V2P,DVP,NLIMP                                    J02400
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         J02410
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        J02420
     *              NLTEFL,LNFIL4,LNGTH4                                  J02430
      DIMENSION DUMMY(2),PNLHDR(2)                                        J02440
      DIMENSION S(*)                                                      J02450
C                                                                         J02460
      EQUIVALENCE (PNLHDR(1),V1P)                                         J02470
C                                                                         J02480
C----------------------------------------------------------------------   J02490
C                                                                         J02500
      CALL CPUTIM (TIME1)                                                 J02510
      IDUM1 = 0                                                           J02520
      IDUM2 = 0                                                           J02530
      ISCANT = MOD(ISCAN,1000)                                            J02540
      IF (JTREM.EQ.0.AND.ISCANT.GE.1) GO TO 70                            J02550
      IF (ISCAN.LT.1) THEN                                                J02560
         IF (JTREM.EQ.1) IDUM1 = 1                                        J02570
         IF (JTREM.EQ.2) IDUM2 = 1                                        J02580
      ENDIF                                                               J02590
   10 CALL BUFIN (IFILE,IEOFSC,PNLHDR(1),NPHDRF)                          J02600
      IF (IEOFSC.LE.0) THEN                                               J02610
         WRITE (IPR,900)                                                  J02620
         GO TO 60                                                         J02630
      ELSE                                                                J02640
         VMIN = V1P                                                       J02650
         VMAX = V2P                                                       J02660
         DVI = DVP                                                        J02670
         NNI = NLIMP                                                      J02680
      ENDIF                                                               J02690
C                                                                         J02700
      IF ((IDATA.EQ.-1).AND.(VMIN.GT.VBOT)) WRITE (IPR,905)               J02710
      IDATA = 0                                                           J02720
      IF (VMAX.GE.VBOT) GO TO 20                                          J02730
      IF (IDUM2.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)                J02740
      CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)                                J02750
      IF (IDUM1.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)                J02760
      GO TO 10                                                            J02770
C                                                                         J02780
   20 IF (JTREM.EQ.0) THEN                                                J02790
         CALL BUFIN (IFILE,IEOFSC,S(1),NNI)                               J02800
         IF (JEMIT.NE.2.AND.ICNVRT.EQ.0) THEN                             J02810
            DO 30 I = 1, NNI                                              J02820
               SI = S(I)                                                  J02830
               S(I) = 1.                                                  J02840
               IF (SI.GT.0.) THEN                                         J02850
                  IF (SI.GE.ARGMIN) THEN                                  J02860
                     S(I) = EXPMIN                                        J02870
                  ELSE                                                    J02880
                     S(I) = EXP(-SI)                                      J02890
                  ENDIF                                                   J02900
               ENDIF                                                      J02910
   30       CONTINUE                                                      J02920
         ENDIF                                                            J02930
      ELSE                                                                J02940
C                                                                         J02950
         IF (IDUM2.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)             J02960
         CALL BUFIN (IFILE,IEOFSC,S(1),NNI)                               J02970
         IF (IDUM1.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)             J02980
      ENDIF                                                               J02990
C                                                                         J03000
      IF (JABS.EQ.1.AND.ICNVRT.EQ.0) THEN                                 J03010
         DO 40 I = 1, NNI                                                 J03020
            S(I) = 1.-S(I)                                                J03030
   40    CONTINUE                                                         J03040
      ENDIF                                                               J03050
      IF (JEMIT.EQ.2.AND.ICNVRT.EQ.0) THEN                                J03060
         DO 50 I = 1, NNI                                                 J03070
            S(I) = -ALOG(S(I))                                            J03080
   50    CONTINUE                                                         J03090
      ENDIF                                                               J03100
C                                                                         J03110
      VMIN = VMIN-4.0*DVI                                                 J03120
      NNI = NNI+4                                                         J03130
      ILO = 1                                                             J03140
      IHI = NNI                                                           J03150
      DIF = (VMIN-VBOT)/DVI                                               J03160
      IF (DIF.LT.0.) ILO = -DIF+1.5                                       J03170
      IF (VMAX.GT.VTOP) THEN                                              J03180
         IHI = (VTOP-VMIN)/DVI+1.5                                        J03190
         IDATA = 1                                                        J03200
      ENDIF                                                               J03210
C                                                                         J03220
   60 CALL CPUTIM (TIME2)                                                 J03230
      TIMRDF = TIMRDF+TIME2-TIME1                                         J03240
C                                                                         J03250
      RETURN                                                              J03260
C                                                                         J03270
   70 WRITE (IPR,910) JTREM,ISCAN                                         J03280
C                                                                         J03290
      RETURN                                                              J03300
C                                                                         J03310
  900 FORMAT ('0 ********** END OF FILE ENCOUNTERED; CHECK IFILE ')       J03320
  905 FORMAT ('0 ********** FIRST VALUE USED ON IFILE; CHECK IFILE ')     J03330
  910 FORMAT (' ERROR IN INPUT',/,'  JTREM =',I2,'  ISCAN=',I5)           J03340
C                                                                         J03350
      END                                                                 J03360
C
C     --------------------------------------------------------------
C
      SUBROUTINE OTPANL (R1,JFILE,NPTS)                                   J03370
C                                                                         J03380
      IMPLICIT REAL*8          (V)                                     ! J03390
C                                                                         J03400
C     SUBROUTINE OTPANL OUTPUTS THE RESULTS OF THE INTERPOLATION ON       J03410
C     TO FILE JFILE                                                       J03420
C                                                                         J03430
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         J03440
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        J03450
     *              NLTEFL,LNFIL4,LNGTH4                                  J03460
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            J03470
      COMMON /OUTPNL/ V1P,V2P,DVP,NLIM                                    J03480
      DIMENSION PNLHDR(2),R1(*)                                           J03490
C                                                                         J03500
      EQUIVALENCE (PNLHDR(1),V1P)                                         J03510
C                                                                         J03520
C----------------------------------------------------------------------   J03530
C                                                                         J03540
      CALL CPUTIM (TIME1)                                                 J03550
      IF (NLIM.LE.0) GO TO 20                                             J03560
C                                                                         J03570
      CALL BUFOUT (JFILE,PNLHDR(1),NPHDRF)                                J03580
      CALL BUFOUT (JFILE,R1(1),NLIM)                                      J03590
C                                                                         J03600
      IF (NPTS.GT.0) THEN                                                 J03610
         WRITE (IPR,900) V1P,V2P,DVP,NLIM                                 J03620
         WRITE (IPR,905)                                                  J03630
         NNPTS = NPTS                                                     J03640
         IF (NPTS.GT.(NLIM/2)+1) NNPTS = NLIM/2+1                         J03650
         IJLIM = NLIM-NNPTS+1                                             J03660
         DO 10 IJ = 1, NNPTS                                              J03670
            IK = IJ+IJLIM-1                                               J03680
            VIJ = V1P+FLOAT(IJ-1)*DVP                                     J03690
            VIK = V1P+FLOAT(IK-1)*DVP                                     J03700
            WRITE (IPR,910) IJ,VIJ,R1(IJ),IK,VIK,R1(IK)                   J03710
   10    CONTINUE                                                         J03720
      ENDIF                                                               J03730
C                                                                         J03740
   20 CALL CPUTIM (TIME2)                                                 J03750
      TIMPNL = TIMPNL+TIME2-TIME1                                         J03760
C                                                                         J03770
      RETURN                                                              J03780
C                                                                         J03790
  900 FORMAT ('0 V1P =',F12.5,' V2P =',F12.5,' DVP =',F12.8,' NLIM =',    J03800
     *        I8)                                                         J03810
  905 FORMAT ('0')                                                        J03820
  910 FORMAT (I5,0PF12.5,1PE12.5,I15,0PF12.5,1PE12.5)                     J03830
C                                                                         J03840
      END                                                                 J03850
C
C     --------------------------------------------------------------
C
      SUBROUTINE INTERP (IFILE,JFILE,I4PT,IBOUND,NPTS,JTREM,ISCAN,        J03860
     *                   JEMIT,RSTAT,ICNVRT)                              J03870
C                                                                         J03880
      IMPLICIT REAL*8          (V)                                     ! J03890
C                                                                         J03900
C**********************************************************************   J03910
C     THIS SUBROUTINE INTERPOLATES THE SPECTRAL DATA FROM IFILE, ON       J03920
C     A GRID DEFINED BY V1I,V2I,DVI, AND NNI, ONTO THE GRID FROM          J03930
C     V1 TO V2 WITH A DV OF DVO AND WRITES THE RESULT TO JFILE.           J03940
C     THE INTERPOLATION IS EITHER LINEAR (I4PT = 0) OR 4 POINT (I4PT      J03950
C     = 1).   IBOUND IS THE NUMBER OF POINTS NEEDED FROM THE PREVIOUS     J03960
C     INPUT PANEL, WHILE NPTS IS THE NUMBER OF POINTS TO BE PRINTED       J03970
C     AT THE BEGINNING AND END OF EACH OUTPUT PANEL. JTREM,ISCAN, AND     J03980
C     JEMIT RELATE TO THE INPUT DATA AND ARE NEEDED BY RDPANL.            J03990
C     RSTAT(3) RETURNS THE SUM, MIN, AND MAX OF THE INTERPOLATED          J04000
C     SPECTRUM.                                                           J04010
C**********************************************************************   J04020
C                                                                         J04030
C     THE INPUT DATA WILL BE PUT INTO T(5) WITH THE LAST                  J04040
C     IBOUND POINTS OF THE PREVIOUS PANEL PUT INTO T(1-4)                 J04050
C                                                                         J04060
      COMMON T(2410),R(2401)                                              J04070
      DIMENSION S(2406)                                                   J04080
      EQUIVALENCE (T(5),S(1))                                             J04090
C                                                                         J04100
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   J04110
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       J04120
     *               DVSC,XDUM,V1SHFT                                     J04130
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         J04140
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        J04150
     *              NLTEFL,LNFIL4,LNGTH4                                  J04160
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 J04170
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            J04180
      COMMON /INPNL/ V1I,V2I,DVI,NNI                                      J04190
      COMMON /OUTPNL/ V1J,V2J,DVJ,NNJ                                     J04200
      DIMENSION C1(0:202),C2(0:202),C3(0:202),C4(0:202),RSTAT(3)          J04210
C                                                                         J04220
      DATA NUMCOF / 201 /                                                 J04230
C                                                                         J04240
      CALL CPUTIM (TIME1)                                                 J04250
C                                                                         J04260
C     SET UP FOUR POINT INTERPOLATION COEFICIENTS FOR P FOR 201           J04270
C     POINTS BETWEEN 0 AND 1.0, with an extra point at each end           J04280
C                                                                         J04290
      IF (I4PT.NE.0) THEN                                                 J04300
         XNUMCF = FLOAT(NUMCOF)
         DO 10 IP = 0, NUMCOF+1                                           J04310
            P = (FLOAT(IP)-1.0)/(XNUMCF-1.0)                              J04320
            PP = P**2                                                     J04330
            C1(IP) = -P/2.0*(1-P)**2                                      J04340
            C2(IP) = 1.0-PP*(3.0-2.0*P)+PP/2.0*(1.0-P)                    J04350
            C3(IP) = PP*(3.0-2.0*P)+P/2.0*(1.0-P)**2                      J04360
            C4(IP) = -PP/2*(1.0-P)                                        J04370
   10    CONTINUE                                                         J04380
      ENDIF                                                               J04390
C                                                                         J04400
C     V1 IS REQUESTED LOWER V, VBOT = V1-VBOUND.  VBOUND ALLOWS FOR       J04410
C     4 POINT INTERPOLATION OF THE FIRST DATA POINT.                      J04420
C     THE FIRST PANEL OF DATA IS STORED IN T, STARTING AT T(5)            J04430
C     WITH A CORRESPONDING WAVENUMBER V1I.                                J04440
C     T(1-4) ARE USED TO STORE THE LAST IBOUND POINTS FROM THE            J04450
C     PREVIOUS PANEL, BUT ARE ZEROED OUT FOR THE FIRST PANEL.             J04460
C     THE INTERPOLATED POINTS ARE STORED IN THE ARRAY R.                  J04470
C     THE INDEX I REFERS TO THE INPUT POINTS, J TO THE OUTPUT POINTS.     J04480
C     P VARIES FROM 0 TO 1 AND IS THE FRACTIONAL DISTANCE OF THE          J04490
C     CURRENT OUTPUT WAVENUMBER VJ TO THE NEXT LOWEST INPUT WAVENUMBER    J04500
C     RELATIVE TO THE INPUT DV: P = (VJ-VI(II))/DVI                       J04510
C     INRANG IS 0 IF VJ IS WITHIN THE RANGE OF THE INPUT DATA, -1         J04520
C     IF VJ IS LESS THAN THE INPUT DATA, AND 1 IF IT IS GREATER           J04530
C                                                                         J04540
C     INITIALIZE THE VARIABLES                                            J04550
C                                                                         J04560
      V1J = V1                                                            J04570
      DVJ = DVO                                                           J04580
      VRATIO = DVJ/DVI                                                    J04590
      VJ = V1J                                                            J04600
C                                                                         J04610
      RMIN = 1.0E15                                                       J04620
      RMAX = -1.0                                                         J04630
      RSUM = 0.0                                                          J04640
C                                                                         J04650
C     EXTRAPOLATE DOWN TO V1I-DVI SO THAT THE POINT I=4 IS AVAILABLE      J04660
C     FOR THE FIRST PANEL.  THIS ALLOWS 4 POINT INTERPOLATION BETWEEN     J04670
C     V1I AND V1I+DVI                                                     J04680
C                                                                         J04690
      T(4) = 2.0*T(5)-T(6)                                                J04700
C                                                                         J04710
C     LOOP OVER THE OUTPUT PANELS                                         J04720
C                                                                         J04730
C     IF V1J .LT. V1I, THEN ZERO FILL UP TO V1I.                          J04740
C                                                                         J04750
   20 IF (V1J.LT.V1I) THEN                                                J04760
         J1 = 1                                                           J04770
         J2 = MIN(INT((V1I-V1J)/DVJ)+1,2400)                              J04780
C                                                                         J04790
C     FILL IN                                                             J04800
C                                                                         J04810
         DO 30 J = J1, J2                                                 J04820
            R(J) = 0.0                                                    J04830
   30    CONTINUE                                                         J04840
C                                                                         J04850
            V2J = V1J+DVJ*(J2-1)                                          J04870
            NNJ = J2                                                      J04880
            CALL OTPANL (R,JFILE,NPTS)                                    J04890
            V1J = V2J+DVJ                                                 J04900
            VJ = V1J                                                      J04910
C                                                                         J04950
         GO TO 20                                                         J04960
      ENDIF                                                               J04970
C                                                                         J04980
C     AT THIS POINT, VJ >= V1I                                            J04990
C                                                                         J05000
   40 CONTINUE                                                            J05010
C                                                                         J05020
C     I INDEXES THE LARGEST VI .LE. VJ                                    J05030
C     AND AT THIS POINT SHOULD .GE. 1.                                    J05040
C                                                                         J05050
      I = (VJ-V1I)/DVI+1.00001                                            J05060
      IF (I.LT.1) THEN                                                    J05070
         WRITE (IPR,*) ' INTERP-ERROR: I SHOULD >= 1, IS ',I              J05080
      ENDIF                                                               J05090
      VI = V1I+DVI*FLOAT(I-1)                                             J05100
C                                                                         J05110
C     P IS INCREMENTED BY ADDING DVJ/DVI BUT WILL BE REINITIALIZED        J05120
C     HERE FOR EACH OUTPUT PANEL TO AVOID THE ACCUMULATION OF             J05130
C     TRUNCATION ERRORS                                                   J05140
C                                                                         J05150
      P = (VJ-VI)/DVI                                                     J05160
C                                                                         J05170
      J1 = INT((VJ-V1J)/DVJ+1.001)                                        J05180
      J2 = MIN(INT((V2-V1J)/DVJ+1.001),INT((V2I-DVI-V1J)/DVJ+1.),2400)    J05190
C                                                                         J05200
C     LOOP OVER A SINGLE OUTPUT PANEL                                     J05210
C                                                                         J05220
      IF (I4PT.GT.0) THEN                                                 J05230
C                                                                         J05240
C     4 POINT INTERPOLATION                                               J05250
C                                                                         J05260
         DO 50 J = J1, J2                                                 J05270
C                                                                         J05280
C     PERFORM INTERPOLATION                                               J05290
C                                                                         J05300
            IP = P*XNUMCF+1.00001                                         J05310
            R(J) = C1(IP)*T(I-1)+C2(IP)*T(I)+C3(IP)*T(I+1)+               J05320
     *             C4(IP)*T(I+2)                                          J05330
C                                                                         J05340
C     ACCUMULATE STATISTICS                                               J05350
C                                                                         J05360
            RMIN = MIN(RMIN,R(J))                                         J05370
            RMAX = MAX(RMAX,R(J))                                         J05380
            RSUM = RSUM+R(J)                                              J05390
C                                                                         J05400
C     INCREMENT P AND I                                                   J05410
C                                                                         J05420
            P = P+VRATIO                                                  J05430
            IF (P.GE.1.0) THEN                                            J05440
               I = I+P                                                    J05450
               P = P-FLOAT(INT(P))                                        J05460
            ENDIF                                                         J05470
C                                                                         J05480
   50    CONTINUE                                                         J05490
C                                                                         J05500
      ELSE                                                                J05510
C                                                                         J05520
C     LINEAR INTERPOLATION                                                J05530
C                                                                         J05540
         DO 60 J = J1, J2                                                 J05550
C                                                                         J05560
C     PERFORM INTERPOLATION                                               J05570
C                                                                         J05580
            R(J) = T(I)*(1.0-P)+T(I+1)*P                                  J05590
C                                                                         J05600
C     ACCUMULATE STATISTICS                                               J05610
C                                                                         J05620
            RMIN = MIN(RMIN,R(J))                                         J05630
            RMAX = MAX(RMAX,R(J))                                         J05640
            RSUM = RSUM+R(J)                                              J05650
C                                                                         J05660
C     INCREMENT P AND I                                                   J05670
C                                                                         J05680
            P = P+VRATIO                                                  J05690
            IF (P.GE.1.0) THEN                                            J05700
               I = I+P                                                    J05710
               P = P-FLOAT(INT(P))                                        J05720
            ENDIF                                                         J05730
   60    CONTINUE                                                         J05740
C                                                                         J05750
      ENDIF                                                               J05760
C                                                                         J05770
C     VJ IS THE FREQUENCY OF THE NEXT OUTPUT POINT (NOT THE LAST          J05780
C     POINT IS THE CURRENT PANEL)                                         J05790
C                                                                         J05800
      VJ = V1J+DVJ*J2                                                     J05810
C                                                                         J05820
C     IF THE OUTPUT PANEL IS FULL OR IF V2 REACHED,                       J05830
C     WRITE OUT THE PANEL                                                 J05840
C                                                                         J05850
      IF (J2.GE.2400.OR.VJ.GE.V2) THEN                                    J05860
         NNJ = J2                                                         J05870
         V2J = V1J+DVJ*(J2-1)                                             J05880
         CALL OTPANL (R,JFILE,NPTS)                                       J05890
         V1J = V2J+DVJ                                                    J05900
         J2 = 0                                                           J05910
      ENDIF                                                               J05920
C                                                                         J05930
C     IF REACHED V2, THEN FINISH                                          J05940
C                                                                         J05950
      IF (VJ.GE.V2) GO TO 100                                             J05960
C                                                                         J05970
C     IF THE INPUT FILE REACHED AN EOF, THEN ZERO FILL TO END             J05980
C                                                                         J05990
      IF (IEOFSC.LE.0) GO TO 80                                           J06000
C                                                                         J06010
C     IF THE DATA FROM CURRENT INPUT PANEL IS EXHAUSTED, GET MORE         J06020
C                                                                         J06030
      IF (I.GE.NNI-2) THEN                                                J06040
C                                                                         J06050
C     SHIFT THE LAST IBOUND POINTS DOWN TO T(1-4)                         J06060
C                                                                         J06070
         DO 70 II = 1, IBOUND                                             J06080
            T(II) = T(II+NNI-IBOUND)                                      J06090
   70    CONTINUE                                                         J06100
C                                                                         J06110
C     GET THE NEXT PANEL OF DATA AND RESET I                              J06120
C                                                                         J06130
         CALL RDPANL (S,JTREM,IFILE,ISCAN,JEMIT,ICNVRT)                   J06140
         IF (IEOFSC.LE.0) THEN                                            J06150
C                                                                         J06160
C     IF EOF ON INPUT FILE, THEN EXTRAPOLATE OUT TWO MORE                 J06170
C     POINTS BEYOND I=NNI SO THAT 4 POINT INTERPOLATION CAN               J06180
C     BE PERFORMED UP TO VJ=V2I. (ACTUALLY, ONLY T(NNI+1) NEED            J06190
C     BE EXTRAPOLATED, T(NNI+2) NEED ONLY BE DEFINED.)                    J06200
C     EXTEND THE INPUT PANEL BY ONE POINT  AND LOOP AROUND THE            J06210
C     INTERPOLATION ONE LAST TIME                                         J06220
C                                                                         J06230
            T(NNI+1) = 2.0*T(NNI)-T(NNI-1)                                J06240
            T(NNI+2) = 0.0                                                J06250
            V2I = V2I+DVI                                                 J06260
         ENDIF                                                            J06270
      ENDIF                                                               J06280
C                                                                         J06290
C     LOOP BACK                                                           J06300
C                                                                         J06310
      GO TO 40                                                            J06320
C                                                                         J06330
   80 CONTINUE                                                            J06340
      J1 = J2+1                                                           J06350
      J2 = MIN(INT((V2-V1J)/DVJ+1.0001),2400)                             J06360
C                                                                         J06370
      DO 90 J = J1, J2                                                    J06380
         R(J) = 0.0                                                       J06390
   90 CONTINUE                                                            J06400
      VJ = V1J+DVJ*J2                                                     J06410
C                                                                         J06420
C     IF THE OUTPUT PANEL IS FULL OR IF V2 REACHED,                       J06430
C     WRITE OUT THE PANEL                                                 J06440
C                                                                         J06450
      IF (J2.GE.2400.OR.VJ.GE.V2) THEN                                    J06460
         NNJ = J2                                                         J06470
         V2J = V1J+DVJ*(J2-1)                                             J06480
         CALL OTPANL (R,JFILE,NPTS)                                       J06490
         V1J = V2J+DVJ                                                    J06500
         J2 = 0                                                           J06510
      ENDIF                                                               J06520
C                                                                         J06530
C     IF REACHED V2, THEN FINISH                                          J06540
C                                                                         J06550
      IF (VJ.LT.V2) GO TO 80                                              J06560
C                                                                         J06570
  100 CONTINUE                                                            J06580
C                                                                         J06590
      RSTAT(1) = RSUM*DVJ                                                 J06600
      RSTAT(2) = RMIN                                                     J06610
      RSTAT(3) = RMAX                                                     J06620
C                                                                         J06630
      CALL CPUTIM (TIME2)                                                 J06640
      TIMCNV = TIMCNV+TIME2-TIME1                                         J06650
C                                                                         J06660
      RETURN                                                              J06670
C                                                                         J06680
      END                                                                 J06690
C
C     --------------------------------------------------------------
C
      SUBROUTINE FLTRFN (IFILE)                                           L00010
C                                                                         L00020
C     NFLTPT sets the maximum number of points in the incoming filter
C
      PARAMETER (NFLTPT = 1001)
C
      IMPLICIT REAL*8 (V)
C                                                                         L00040
      COMMON S(2650),R1(3750)                                             L00050
C                                                                         L00060
      character*8      XID,       HMOLID,      YID        
      real*8               SECANT,       XALTZ 
C                                                                         L00080
      COMMON /HVERSN/  HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,
     *                HVROPR,HVRPST,HVRPLT,HVRTST,HVRUTL,HVRXMR
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       L00090
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DVC,V1C,V2C,TBOUND,   L00100
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    L00110
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   L00120
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       L00130
     *               DVSC,XDUM,V1SHFT                                     L00140
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 L00150
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            L00160
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      L00170
      COMMON /COMFLT/ V1F,V2F,DVF,NPTS,NPTF,JEMIT,IUNIT,IFILST,NIFILS,    L00180
     *                HEDDR(9),XF(NFLTPT),SUMFLT                          L00190
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         L00200
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        L00210
     *              NLTEFL,LNFIL4,LNGTH4                                  L00220
      DIMENSION FILHDR(2)                                                 L00230
      CHARACTER*80 CVAR                                                   L00240
      CHARACTER*8 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *            HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR
C                                                                         L00250
      EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),                 L00260
     *            (FSCDID(6),ISCAN) , (FSCDID(7),IPLOT),                  L00270
     *            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                  L00280
     *            (FSCDID(12),SCNID) , (FSCDID(13),HWHM),                 L00290
     *            (FSCDID(16),LAYR1)                                      L00300
C                                                                         L00310
C
C     ASSIGN SCCS VERSION NUMBER TO MODULE 
C
      HVRPST = '$Revision$'
C
      NLIMF = 2401                                                        L00320
      NREN = 0                                                            L00330
      IPRT = 1                                                            L00340
      NSHIFT = 0                                                          L00350
   10 READ (IRD,900) V1F,DVF,NPTF,JEMIT,IUNIT,IFILST,NIFILS,HEDDR         L00360
C
C     Test to ensure NPTF is less than NFLTPT, the maximum number
C     of filter points allowed
C
      IF (NPTF.GT.NFLTPT) THEN
         WRITE(IPR,*) 'FLTRFN: NPTS > NFLTPT limit', NLFTPT
         STOP 'FLTRFN: NPTS > NFLTPT limit'
      ENDIF
C
      JABS = 0                                                            L00370
      IF (JEMIT.GE.0) GO TO 20                                            L00380
      JEMIT = 0                                                           L00390
      JABS = 1                                                            L00400
   20 IEOFT = 1                                                           L00410
      IF (IUNIT.LE.0) IUNIT = IFILE                                       L00420
      IFILE = IUNIT                                                       L00430
      IF (NIFILS.LE.0) NIFILS = 99                                        L00440
C                                                                         L00450
      IF (V1F.LT.0) RETURN                                                L00460

c
c     DVF < 0 option flags V1F value to be the center frequency
c     Check that there NPTF is odd (to ensure a center frequency),
c     save center frequency value, and reset V1F to endpoint value.

      if (DVF.lt.0.) then

         dvf = abs(dvf)

         if (mod((nptf-1),2).ne.0) then
            write(*,*) 'Use of V1F as center frequency requires odd
     *           number of points'
            write(ipr,*) 'Use of V1F as center frequency requires odd
     *           number of points, stopping in FLTFRN'
            stop 'FLTRFN'
         endif

         V1F_center = V1F
         nptf_half = (abs(nptf)-1)/2
         V1F = V1F_center - DVF*float(nptf_half)
         write(ipr,*) ' ``````````````````````````````'
         write(ipr,*) ' Use of V1F as center frequency:'
         write(ipr,*) '   V center = ',V1F_center
         write(ipr,*) '   V1F      = ',V1F
         write(ipr,*) " ''''''''''''''''''''''''''''''"
      endif

C                                                                         L00470
      WRITE (IPR,905)                                                     L00480
      REWIND IFILE                                                        L00490
      IF (IFILST.GT.1) CALL SKIPFL (IFILST-1,IFILE,IEOF)                  L00500
      IEOFSC = 0                                                          L00510
      ISTOP = 0                                                           L00520
C                                                                         L00530
      IF (NPTF.LE.0) GO TO 30                                             L00540
      NPTS = NPTF                                                         L00550
   30 V2F = V1F+DVF*FLOAT(NPTS-1)                                         L00560
      WRITE (IPR,910) V1F,V2F,DVF,NPTF,JEMIT,JABS,IUNIT,IFILST,NIFILS,    L00570
     *                HEDDR                                               L00580
      V1 = V1F                                                            L00590
      V2 = V2F                                                            L00600
      DV = DVF                                                            L00610
      IF (NPTF.LE.0) GO TO 40                                             L00620
      READ (IRD,915) CVAR                                                 L00630
      READ (IRD,CVAR) (XF(I),I=1,NPTS)                                    L00640
      WRITE (IPR,CVAR) (XF(I),I=1,NPTS)                                   L00650
C                                                                         L00660
C     MAKE ADJUSTMENT FOR END POINT CORRECTIONS                           L00670
C                                                                         L00680
      XF(1) = 0.5*XF(1)                                                   L00690
      XF(NPTS) = 0.5*XF(NPTS)                                             L00700
   40 SUMFLT = 0.0                                                        L00710
      DO 50 I = 1, NPTS                                                   L00720
         SUMFLT = SUMFLT+XF(I)                                            L00730
   50 CONTINUE                                                            L00740
      SUMFLT = SUMFLT*DVF                                                 L00750
C                                                                         L00760
   60 RFILTR = 0.0                                                        L00770
      VFT = V1                                                            L00780
      VBOT = V1                                                           L00790
      VTOP = V2                                                           L00800
      TIMRDF = 0.0                                                        L00810
      TIMCNV = 0.0                                                        L00820
      CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)                            L00830
      ISCHDR = ISCAN                                                      L00840
      IF (ISCAN.LE.0.OR.SCNID.EQ.-99.) ISCAN = 0                          L00850
      IF (ISCHDR.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCHDR                   L00860
      IF (MOD(ISCAN,1000).EQ.0) GO TO 70                                  L00870
      JEMSCN = SCNID/100.                                                 L00880
      IF (JEMIT.EQ.JEMSCN) GO TO 70                                       L00890
      WRITE (IPR,920)                                                     L00900
      CALL SKIPFL (1,IFILE,IEOF)                                          L00910
      IF (IEOF.EQ.0) GO TO 10                                             L00920
      GO TO 60                                                            L00930
   70 CONTINUE                                                            L00940
      IF (IEOF.LT.1) GO TO 10                                             L00950
C                                                                         L00960
C     JEMIT=-1 FILTER PASSED OVER ABSORPTION                              L00970
C     JEMIT=0  FILTER PASSED OVER TRANSMISSION                            L00980
C     JEMIT=1  FILTER PASSED OVER EMISSION                                L00990
C                                                                         L01000
      JTREM = -1                                                          L01010
      IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0                        L01020
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2                        L01030
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1                        L01040
      ISCANT = MOD(ISCAN,1000)                                            L01050
      IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 2                       L01060
      IF (JTREM.LT.0) GO TO 10                                            L01070
      WRITE (IPR,925) XID,(YID(M),M=1,2)                                  L01080
      WRITE (IPR,930) LAYR1,LAYER                                         L01090
      WRITE (IPR,935) SECANT,PAVE,TAVE,DVC,V1C,V2C                        L01100
      WRITE (IPR,940) WBROAD,(HMOLID(M),WK(M),M=1,NMOL)                   L01110
      WRITE (IPR,945) V1F,V2F,DVF,NPTF,IEMIT,JEMIT,JABS,IUNIT,IFILST,     L01120
     *   NIFILS,HEDDR                                                     L01130
      IDATA = -1                                                          L01140
   80 CALL CPUTIM (TIMEO)                                                 L01150
      CALL RDSCAN (S,JTREM,IFILE,ISCAN,IPRT)                              L01160
      CALL CPUTIM (TIME)                                                  L01170
      TIMRDF = TIMRDF+TIME-TIMEO                                          L01180
      IF (IEOFSC.NE.1) GO TO 90                                           L01190
      CALL CNVFLT (S,RFILTR,XF)                                           L01200
      IF (IDATA.EQ.1) GO TO 90                                            L01210
      GO TO 80                                                            L01220
   90 CALL CPUTIM (TIME)                                                  L01230
      WRITE (IPR,950) TIME,TIMRDF,TIMCNV                                  L01240
      IF (JEMIT.EQ.1) GO TO 100                                           L01250
      TRNSM = RFILTR                                                      L01260
      RFILTR = RFILTR*DVC/SUMFLT                                          L01270
      IF (JABS.EQ.0) WRITE (IPR,955) RFILTR,SUMFLT,TRNSM                  L01280
      IF (JABS.EQ.1) WRITE (IPR,960) RFILTR,SUMFLT,TRNSM                  L01290
      GO TO 110                                                           L01300
  100 RFILTR = RFILTR*DVC                                                 L01310
      WRITE (IPR,965) RFILTR,SUMFLT                                       L01320
  110 IF (IEOFSC.EQ.1) CALL SKIPFL (1,IFILE,IEOFSC)                       L01330
      IEOFT = IEOFT+1                                                     L01340
      IF (IEOFT.GT.NIFILS) GO TO 10                                       L01350
      IF (IEOFSC.LT.0) GO TO 60                                           L01360
      GO TO 10                                                            L01370
C                                                                         L01380
  900 FORMAT (2F10.4,5I5,8A4,A3)                                          L01390
  905 FORMAT ('1',/'   ***  FILTER ***',8(' ********** '))                L01400
  910 FORMAT ('0   V1F=',F10.4,' V2F=',F10.4,',DVF=',F10.4,',NPTF =',     L01410
     *        I5,/,'0',10X,', JEMIT= ',I2,', JABS= ',I2,                  L01420
     *        ', INPUT FILE= ',I3,' ,IFILST =',I5,' ,NIFILS =',I5,2X,     L01430
     *        8A4,A3)                                                     L01440
  915 FORMAT (A80)                                                        L01450
  920 FORMAT ('0  RESULT FROM SCANNING FUNCTION INCONSISTENT WITH ',      L01460
     *        'FILTER REQUEST')                                           L01470
  925 FORMAT ('0',//,8(' -------  '),/,'0',10A8,2X,2(1X,A8,1X))           L01480
  930 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)            L01490
  935 FORMAT ('0 SECANT =',F15.5,/'0 PRESS(MB) =',F12.5/'0 TEMP(K) =',    L01500
     *        F11.2,/'0 DV(CM-1) = ',F12.8,/'0 V1(CM-1) = ',F12.6,/       L01510
     *        '0 V2(CM-1) = ',F12.6)                                      L01520
  940 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//5X,'WBROAD = ',       L01530
     *        1PE10.3,/(5X,A6,' = ',1PE10.3))                             L01540
  945 FORMAT ('0   V1F=',F10.4,' V2F=',F10.4,',DVF=',F10.4,',NPTF =',     L01550
     *        I5,/'0 , IEMIT= ',I2,', JEMIT= ',I2,', JABS= ',I2,          L01560
     *        ', INPUT FILE= ',I3,' ,IFILST =',I5,' ,NIFILS =',I5,2X,     L01570
     *        8A4,A3)                                                     L01580
  950 FORMAT ('0',5X,'TIME =',F7.3,',  READ =',F6.3,',  CONV. =',F7.3)    L01590
  955 FORMAT ('0  INTEGRATED TRANSMISSION = ',1PE14.5,                    L01600
     *        '  NORMALIZATION OF  THE FILTER = ',E14.5,/                 L01610
     *        '0 UNNORMALIZED INTEGRATED TRANSMISSION =  ',E14.5)         L01620
  960 FORMAT ('0  INTEGRATED ABSORPTION = ',1PE14.5,                      L01630
     *        '  NORMALIZATION OF  THE FILTER = ',E14.5,/                 L01640
     *        '0 UNNORMALIZED INTEGRATED ABSORPTION =    ',E14.5)         L01650
  965 FORMAT ('0 INTEGRATED EMISSION = ',1PE14.5,                         L01660
     *        '  NORMALIZATION OF THE',' FILTER = ',E14.5)                L01670
C                                                                         L01680
      END                                                                 L01690
C
C     --------------------------------------------------------------
C
      SUBROUTINE FLTRRD (IFILE)                                           L01700
C                                                                         L01710
C     NFLTPT sets the maximum number of points in the incoming filter
C
      PARAMETER (NFLTPT = 1001)
C
      IMPLICIT REAL*8           (V) 
C                                                                         L01730
C     READ CONTROL CARD FOR FILTER WITH WEIGHTING FUNCTIONS               L01740
C                                                                         L01750
      COMMON S(2650),R1(3750)                                             L01760
C                                                                         L01770
      character*8      XID,       HMOLID,      YID        
      real*8               SECANT,       XALTZ 
C                                                                         L01790
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       L01800
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DVC,V1C,V2C,TBOUND,   L01810
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    L01820
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   L01830
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       L01840
     *               DVSC,XDUM,V1SHFT                                     L01850
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 L01860
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            L01870
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      L01880
      COMMON /COMFLT/ V1F,V2F,DVF,NPTS,NPTF,JEMIT,IUNIT,IFILST,NIFILS,    L01890
     *                HEDDR(9),XF(NFLTPT),SUMFLT                          L01900
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         L01910
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        L01920
     *              NLTEFL,LNFIL4,LNGTH4                                  L01930
      COMMON /FLFORM/ CFORM                                               L01940
      DIMENSION FILHDR(2)                                                 L01950
C                                                                         L01960
      CHARACTER*80 CVAR                                                   L01970
      CHARACTER CFORM*11,TAPE13*6,CTAPE*4                                 L01980
      LOGICAL OP                                                          L01990
C                                                                         L02000
      EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),                 L02010
     *            (FSCDID(6),ISCAN) , (FSCDID(7),IPLOT),                  L02020
     *            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                  L02030
     *            (FSCDID(12),SCNID) , (FSCDID(13),HWHM),                 L02040
     *            (FSCDID(16),LAYR1)                                      L02050
C                                                                         L02060
      DATA TAPE13 / '      '/,CTAPE / 'TAPE'/                             L02070
C                                                                         L02080
      NLIMF = 2401                                                        L02090
      NSHIFT = 0                                                          L02100
      READ (IRD,900) V1F,DVF,NPTF,JEMIT,NNFILE,HEDDR                      L02110
C                                                                         L02120
      JABS = 0                                                            L02130
      IEOFT = 1                                                           L02140
      IUNIT = IFILE                                                       L02150
      IFILE = IUNIT                                                       L02160
      IFILST = 1                                                          L02170
      NIFILS = 1                                                          L02180
      IF (NNFILE.NE.NFILE.AND.NNFILE.GT.0) THEN                           L02190
         INQUIRE (UNIT=NFILE,OPENED=OP)                                   L02200
         IF (OP) CLOSE (NFILE)                                            L02210
         NFILE = NNFILE                                                   L02220
         INQUIRE (UNIT=NFILE,OPENED=OP)                                   L02230
         IF (.NOT.OP) THEN                                                L02240
            WRITE (TAPE13,905) CTAPE,NFILE                                L02250
            OPEN (NFILE,FILE=TAPE13,STATUS='UNKNOWN',FORM=CFORM)          L02260
            REWIND NFILE                                                  L02270
         ENDIF                                                            L02280
      ENDIF                                                               L02290
C                                                                         L02300
      IF (V1F.LT.0) RETURN                                                L02310
      WRITE (IPR,910)                                                     L02320
      REWIND IFILE                                                        L02330
      IF (IFILST.GT.1) CALL SKIPFL (IFILST-1,IFILE,IEOF)                  L02340
      IEOFSC = 0                                                          L02350
      ISTOP = 0                                                           L02360
C                                                                         L02370
      IF (NPTF.LE.0) GO TO 10                                             L02380
      NPTS = NPTF                                                         L02390
   10 V2F = V1F+DVF*FLOAT(NPTS-1)                                         L02400
      WRITE (IPR,915) V1F,V2F,DVF,NPTF,JEMIT,JABS,IUNIT,IFILST,NIFILS,    L02410
     *                NFILE,HEDDR                                         L02420
      V1 = V1F                                                            L02430
      V2 = V2F                                                            L02440
      DV = DVF                                                            L02450
      IF (NPTF.LE.0) GO TO 20                                             L02460
      READ (IRD,920) CVAR                                                 L02470
      READ (IRD,CVAR) (XF(I),I=1,NPTS)                                    L02480
      WRITE (IPR,CVAR) (XF(I),I=1,NPTS)                                   L02490
C                                                                         L02500
C     MAKE ADJUSTMENT FOR END POINT CORRECTIONS                           L02510
C                                                                         L02520
      XF(1) = 0.5*XF(1)                                                   L02530
      XF(NPTS) = 0.5*XF(NPTS)                                             L02540
   20 SUMFLT = 0.0                                                        L02550
      DO 30 I = 1, NPTS                                                   L02560
         SUMFLT = SUMFLT+XF(I)                                            L02570
   30 CONTINUE                                                            L02580
      SUMFLT = SUMFLT*DVF                                                 L02590
C                                                                         L02600
      RETURN                                                              L02610
C                                                                         L02620
  900 FORMAT (2F10.4,I5,I5,I5,10X,8A4,A3)                                 L02630
  905 FORMAT (A4,I2.2)                                                    L02640
  910 FORMAT ('1',/'   ***  FILTER ***',8(' ********** '))                L02650
  915 FORMAT ('0   V1F=',F10.4,' V2F=',F10.4,',DVF=',F10.4,               L02660
     *        ',NPTF =',I5,/,'0',10X,', JEMIT= ',I2,', JABS= ',           L02670
     *        I2,', INPUT FILE= ',I3,' ,IFILST =',I5,' ,NIFILS =',        L02680
     *        I5,', OUTPUT FILE= ',I3,2X,8A4,A3)                          L02690
  920 FORMAT (A80)                                                        L02700
C                                                                         L02710
      END                                                                 L02720
C
C     --------------------------------------------------------------
C
      SUBROUTINE FLTMRG (IFILE,JFILE)                                     L02730
C                                                                         L02740
C     NFLTPT sets the maximum number of points in the incoming filter
C
      PARAMETER (NFLTPT = 1001)
C
      IMPLICIT REAL*8           (V) 
C                                                                         L02760
C     SUBROUTINE FLTMRG CALCULATES AND OUTPUTS THE RESULTS                L02770
C     OF THE FILTER TO FILE JFILE                                         L02780
C                                                                         L02790
      COMMON S(2650),R1(3750)                                             L02800
C                                                                         L02810
      character*8      XID,       HMOLID,      YID        
      real*8               SECANT,       XALTZ 
C                                                                         L02830
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       L02840
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DVC,V1C,V2C,TBOUND,   L02850
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    L02860
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   L02870
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       L02880
     *               DVSC,XDUM,V1SHFT                                     L02890
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   L02900
     *              AVFIX,LAYMA,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,        L02910
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       L02920
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      L02930
     *              EXTID(10)                                             L02940
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 L02950
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            L02960
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      L02970
      COMMON /COMFLT/ V1F,V2F,DVF,NPTS,NPTF,JEMIT,IUNIT,IFILST,NIFILS,    L02980
     *                HEDDR(9),XF(NFLTPT),SUMFLT                          L02990
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         L03000
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        L03010
     *              NLTEFL,LNFIL4,LNGTH4                                  L03020
      COMMON /SPANEL/ V1P,V2P,DV,NLIM                                     L03030
      DIMENSION FILHDR(2),PNLHDR(2),RF(4)                                 L03040
      EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),                 L03050
     *            (FSCDID(6),ISCNHD) , (FSCDID(7),IPLOT),                 L03060
     *            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                  L03070
     *            (FSCDID(12),SCNID) , (FSCDID(13),HWHM),                 L03080
     *            (FSCDID(16),LAYR1) , (PNLHDR(1),V1P)                    L03090
C                                                                         L03100
      NLIMF = 2401                                                        L03110
      NREN = 0                                                            L03120
      IPRT = 1                                                            L03130
      NSHIFT = 0                                                          L03140
      IUNIT = IFILE                                                       L03150
      REWIND IFILE                                                        L03160
C                                                                         L03170
   10 RFILTR = 0.0                                                        L03180
      VFT = V1                                                            L03190
      VBOT = V1                                                           L03200
      VTOP = V2                                                           L03210
      TIMRDF = 0.0                                                        L03220
      TIMCNV = 0.0                                                        L03230
C                                                                         L03240
      CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)                            L03250
C                                                                         L03260
      ISCAN = ISCNHD                                                      L03270
      IF (ISCAN.LE.0.OR.SCNID.EQ.-99.) ISCAN = 0                          L03280
      IF (ISCNHD.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCNHD                   L03290
      ISCNHD = ISCAN+100                                                  L03300
      IF (MOD(ISCAN,1000).EQ.0) GO TO 20                                  L03310
      JEMSCN = SCNID/100.                                                 L03320
      IF (JEMIT.EQ.JEMSCN) GO TO 20                                       L03330
      WRITE (IPR,900)                                                     L03340
      CALL SKIPFL (1,IFILE,IEOF)                                          L03350
C                                                                         L03360
      IF (IEOF.EQ.0) RETURN                                               L03370
C                                                                         L03380
      GO TO 10                                                            L03390
   20 CONTINUE                                                            L03400
C                                                                         L03410
      IF (IEOF.LT.1) RETURN                                               L03420
C                                                                         L03430
C     JEMIT=-1 FILTER PASSED OVER ABSORPTION                              L03440
C     JEMIT=0  FILTER PASSED OVER TRANSMISSION                            L03450
C     JEMIT=1  FILTER PASSED OVER EMISSION                                L03460
C                                                                         L03470
      JTREM = -1                                                          L03480
      IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0                        L03490
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2                        L03500
      IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1                        L03510
      ISCANT = MOD(ISCAN,1000)                                            L03520
      IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 2                       L03530
C                                                                         L03540
      IF (JTREM.LT.0) RETURN                                              L03550
C                                                                         L03560
      WRITE (IPR,905) XID,(YID(M),M=1,2)                                  L03570
      WRITE (IPR,910) LAYR1,LAYER                                         L03580
      WRITE (IPR,915) SECANT,PAVE,TAVE,DVC,V1C,V2C                        L03590
      WRITE (IPR,920) WBROAD,(HMOLID(M),WK(M),M=1,NMOL)                   L03600
      WRITE (IPR,925) V1F,V2F,DVF,NPTF,IEMIT,JEMIT,JABS,IUNIT,IFILST,     L03610
     *                NIFILS,HEDDR                                        L03620
C                                                                         L03630
      XSCID = 100*JEMIT                                                   L03640
      SCNID = XSCID+0.01                                                  L03650
C                                                                         L03660
      CALL BUFOUT (JFILE,FILHDR(1),NFHDRF)                                L03670
      IDATA = -1                                                          L03680
   30 CALL CPUTIM (TIMEO)                                                 L03690
      CALL RDSCAN (S,JTREM,IFILE,ISCAN,IPRT)                              L03700
      CALL CPUTIM (TIME)                                                  L03710
      TIMRDF = TIMRDF+TIME-TIMEO                                          L03720
      IF (IEOFSC.NE.1) GO TO 40                                           L03730
      CALL CNVFLT (S,RFILTR,XF)                                           L03740
      IF (IDATA.EQ.1) GO TO 40                                            L03750
      GO TO 30                                                            L03760
C                                                                         L03770
   40 CALL CPUTIM (TIME)                                                  L03780
      WRITE (IPR,930) TIME,TIMRDF,TIMCNV                                  L03790
      IF (JEMIT.EQ.1) GO TO 50                                            L03800
      TRNSM = RFILTR                                                      L03810
      RFILTR = RFILTR*DVC/SUMFLT                                          L03820
      RF(1) = RFILTR                                                      L03830
      RF(2) = SUMFLT                                                      L03840
      RF(3) = TRNSM                                                       L03850
      RF(4) = PLAY                                                        L03860
      IF (JABS.EQ.0) WRITE (IPR,935) RFILTR,SUMFLT,TRNSM                  L03870
      IF (JABS.EQ.1) WRITE (IPR,940) RFILTR,SUMFLT,TRNSM                  L03880
C                                                                         L03890
C     V1P IS FIRST FREQ OF PANEL                                          L03900
C     V2P IS LAST  FREQ OF PANEL                                          L03910
C                                                                         L03920
      V1P = V1F                                                           L03930
      V2P = V2F                                                           L03940
      DVP = DVF                                                           L03950
      NLIM = 4                                                            L03960
      CALL BUFOUT (JFILE,PNLHDR(1),NPHDRF)                                L03970
      CALL BUFOUT (JFILE,RF(1),NLIM)                                      L03980
      GO TO 60                                                            L03990
C                                                                         L04000
   50 RFILTR = RFILTR*DVC                                                 L04010
      WRITE (IPR,945) RFILTR,SUMFLT                                       L04020
C                                                                         L04030
   60 IF (IEOFSC.EQ.1) CALL SKIPFL (1,IFILE,IEOFSC)                       L04040
      IEOFT = IEOFT+1                                                     L04050
      IF (IEOFT.GT.NIFILS) RETURN                                         L04060
      IF (IEOFSC.LT.0) GO TO 10                                           L04070
C                                                                         L04080
      RETURN                                                              L04090
C                                                                         L04100
  900 FORMAT ('0  RESULT FROM SCANNING FUNCTION INCONSISTENT WITH ',      L04110
     *        'FILTER REQUEST')                                           L04120
  905 FORMAT ('0',//,8(' -------  '),/,'0',10A8,2X,2(1X,A8,1X))           L04130
  910 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)            L04140
  915 FORMAT ('0 SECANT =',F15.5,/'0 PRESS(MB) =',F12.5/'0 TEMP(K) =',    L04150
     *        F11.2,/'0 DV(CM-1) = ',F12.8,/'0 V1(CM-1) = ',F12.6,/       L04160
     *        '0 V2(CM-1) = ',F12.6)                                      L04170
  920 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//5X,'WBROAD = ',       L04180
     *        1PE10.3,/(5X,A6,' = ',1PE10.3))                             L04190
  925 FORMAT ('0   V1F=',F10.4,' V2F=',F10.4,',DVF=',F10.4,',NPTF =',     L04200
     *        I5,/,'0 , IEMIT= ',I2,', JEMIT= ',I2,', JABS= ',I2,         L04210
     *        ', INPUT FILE= ',I3,' ,IFILST =',I5,' ,NIFILS =',I5,2X,     L04220
     *        8A4,A3)                                                     L04230
  930 FORMAT ('0',5X,'TIME =',F7.3,',  READ =',F6.3,',  CONV. =',F7.3)    L04240
  935 FORMAT ('0  INTEGRATED TRANSMISSION = ',1PE14.5,                    L04250
     *        '  NORMALIZATION OF  THE FILTER = ',E14.5,/                 L04260
     *        '0 UNNORMALIZED INTEGRATED TRANSMISSION =  ',E14.5)         L04270
  940 FORMAT ('0  INTEGRATED ABSORPTION = ',1PE14.5,                      L04280
     *        '  NORMALIZATION OF  THE FILTER = ',E14.5,/                 L04290
     *        '0 UNNORMALIZED INTEGRATED ABSORPTION =    ',E14.5)         L04300
  945 FORMAT ('0 INTEGRATED EMISSION = ',1PE14.5,                         L04310
     *        '  NORMALIZATION OF THE    FILTER = ',E14.5)                L04320
C                                                                         L04330
      END                                                                 L04340
C
C     --------------------------------------------------------------
C
      SUBROUTINE CNVFLT (S,RFILTR,XF)                                     L04350
C                                                                         L04360
C     NFLTPT sets the maximum number of points in the incoming filter
C
      PARAMETER (NFLTPT = 1001)
C
      IMPLICIT REAL*8           (V) 
C                                                                         L04380
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 L04390
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      L04400
      COMMON /COMFLT/ V1F,V2F,DVF,NPTS,NPTF,JEMIT,IUNIT,IFILST,NIFILS,    L04410
     *                HEDDR(9),XFS(NFLTPT),SUMFLT                         L04420
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   L04430
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       L04440
     *               DVSC,XDUM,V1SHFT                                     L04450
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            L04460
C                                                                         L04470
      DIMENSION XF(*),S(*)                                                L04480
C                                                                         L04490
      CALL CPUTIM (TIMEO)                                                 L04500
      IMIN = (V1F-V1I)/DVI+1.5                                            L04510
      IMIN = MAX(IMIN,ILO)                                                L04520
      IMAX = (V2F+V1F-V1I)/DVI+1.5                                        L04530
      IMAX = MIN(IMAX,IHI)                                                L04540
      XIF0 = (V1I-V1F)/DVF+1.5                                            L04550
      XDVIF = DVI/DVF                                                     L04560
      DO 10 I = IMIN, IMAX                                                L04570
         IFL = XIF0+XDVIF*FLOAT(I-1)                                      L04580
         RFILTR = RFILTR+S(I)*XF(IFL)                                     L04590
   10 CONTINUE                                                            L04600
      IF (IMAX.LT.IHI) VFT = VFT+((FLOAT(IHI)-FLOAT(ILO))+1.0)*DVI        L04610
      CALL CPUTIM (TIME)                                                  L04620
      TIMCNV = TIMCNV+TIME-TIMEO                                          L04630
C                                                                         L04640
      RETURN                                                              L04650
C                                                                         L04660
      END                                                                 L04670
C
C     --------------------------------------------------------------
C
      SUBROUTINE FLTPRT (IFILE)                                           L04680
C                                                                         L04690
C     NFLTPT sets the maximum number of points in the incoming filter
C
      PARAMETER (NFLTPT = 1001)
C
      IMPLICIT REAL*8           (V)     
C                                                                         L04710
C     SUBROUTINE FLTPRT READS FROM IFILE AND FORMATS OUT THE RESULTS      L04720
C     OF THE FILTERED WEIGHTING FUNCTION TO IPR                           L04730
C                                                                         L04740
      COMMON S(2650),R1(3750)                                             L04750
C                                                                         L04760
      character*8      XID,       HMOLID,      YID        
      real*8               SECANT,       XALTZ 
C                                                                         L04780
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       L04790
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DVC,V1C,V2C,TBOUND,   L04800
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    L04810
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   L04820
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       L04830
     *               DVSC,XDUM,V1SHFT                                     L04840
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   L04850
     *              AVFIX,LAYMA,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,        L04860
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       L04870
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      L04880
     *              EXTID(10)                                             L04890
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 L04900
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL                            L04910
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      L04920
      COMMON /COMFLT/ V1F,V2F,DVF,NPTS,NPTF,JEMIT,IUNIT,IFILST,NIFILS,    L04930
     *                HEDDR(9),XF(NFLTPT),SUMFLT                          L04940
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         L04950
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        L04960
     *              NLTEFL,LNFIL4,LNGTH4                                  L04970
      COMMON /SPANEL/ V1P,V2P,DV,NLIM                                     L04980
      DIMENSION FILHDR(2),PNLHDR(2),RF(4)                                 L04990
C                                                                         L05000
      EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),                 L05010
     *            (FSCDID(6),ISCAN) , (FSCDID(7),IPLOT),                  L05020
     *            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                  L05030
     *            (FSCDID(12),SCNID) , (FSCDID(13),HWHM),                 L05040
     *            (FSCDID(16),LAYR1) , (PNLHDR(1),V1P)                    L05050
C                                                                         L05060
      NLIMF = 2401                                                        L05070
      NSHIFT = 0                                                          L05080
      IUNIT = IFILE                                                       L05090
      REWIND IFILE                                                        L05100
C                                                                         L05110
      WRITE (IPR,900)                                                     L05120
      OLDTRN = 1.0                                                        L05130
      TIMRDF = 0.0                                                        L05140
      CALL CPUTIM (TIMEO)                                                 L05150
   10 CONTINUE                                                            L05160
      CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)                            L05170
      IF (IEOF.EQ.-99) GO TO 10                                           L05180
      IF (IEOF.EQ.0) GO TO 20                                             L05190
C                                                                         L05200
      CALL BUFIN (IFILE,IEOF,PNLHDR(1),NPHDRF)                            L05210
      CALL BUFIN (IFILE,IEOF,RF(1),NLIM)                                  L05220
C                                                                         L05230
      RFILTR = RF(1)                                                      L05240
      SUMFLT = RF(2)                                                      L05250
      TRNSM = RF(3)                                                       L05260
      PLAY = RF(4)                                                        L05270
      DTAU = OLDTRN-RFILTR                                                L05280
      OLDTRN = RFILTR                                                     L05290
      WRITE (IPR,905) LAYER,PLAY,DTAU,RFILTR,TRNSM                        L05300
C                                                                         L05310
      IF (IEOF.EQ.1) GO TO 10                                             L05320
   20 WRITE (IPR,910) SUMFLT                                              L05330
      CALL CPUTIM (TIME)                                                  L05340
      TIMRDF = TIMRDF+TIME-TIMEO                                          L05350
      WRITE (IPR,915) TIME,TIMRDF                                         L05360
C                                                                         L05370
      RETURN                                                              L05380
C                                                                         L05390
  900 FORMAT ('1',15X,                                                    L05400
     *        'FORMATTED RESULTS OF FILTERED WEIGHTING FUNCTIONS',        L05410
     *        2(/),10X,'LAYER',5X,'PRESSURE',4X,'TRANSMISSION',6X,        L05420
     *        'INTEGRATED',6X,'UNNORMALIZED',/,10X,'  #  ',5X,            L05430
     *        '  (MB)  ',5X,'(N-1) - (N)',5X,'TRANSMISSION',5X,           L05440
     *        '  INT TRANS ',/)                                           L05450
  905 FORMAT (10X,I3,5X,F8.3,5X,1PE13.6,4X,1PE13.6,4X,1PE13.6)            L05460
  910 FORMAT ('0',5X,'NORMALIZATION OF THE FILTER =',1PE13.6)             L05470
  915 FORMAT ('0',5X,'TIME =',F7.3,',  IN FLTPRT =',F7.3)                 L05480
C                                                                         L05490
      END                                                                 L05500
