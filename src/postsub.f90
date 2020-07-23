!     path:      $HeadURL$
!     revision:  $Revision$
!     created:   $Date$
!     presently: %H%  %T%
!
!  --------------------------------------------------------------------------
! |  Copyright ©, Atmospheric and Environmental Research, Inc., 2015         |
! |                                                                          |
! |  All rights reserved. This source code is part of the LBLRTM software    |
! |  and is designed for scientific and research purposes. Atmospheric and   |
! |  Environmental Research, Inc. (AER) grants USER the right to download,   |
! |  install, use and copy this software for scientific and research         |
! |  purposes only. This software may be redistributed as long as this       |
! |  copyright notice is reproduced on any copy made and appropriate         |
! |  acknowledgment is given to AER. This software or any modified version   |
! |  of this software may not be incorporated into proprietary software or   |
! |  commercial software offered for sale without the express written        |
! |  consent of AER.                                                         |
! |                                                                          |
! |  This software is provided as is without any express or implied          |
! |  warranties.                                                             |
! |                       (http://www.rtweb.aer.com/)                        |
!  --------------------------------------------------------------------------
!
SUBROUTINE SCANFN (IFILE,JFILE)
!
   IMPLICIT REAL*8          (V)
!
!     DRIVER FOR CONVOLVING INSTRUMENTAL SCANNING FUNCTION
!     WITH SPECTRUM
!
   COMMON S(3850),R1(5000),N1(5000)
!
   character*8      XID,       HMOLID,      YID,SCANID
   real*8               SECANT,       XALTZ
!
   COMMON /CVRPST/HNAMPST, HVRPST
   COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /RSCAN/ V1I,V2I,DVI,NNI
   COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,             &
   &                HWF3,DXF3,NX3,N3MAX
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /SCINF/ HWHM,JEMIT,JFN,SAMPLE,SCANID,NPTS,XF(6018)
   COMMON /FLFORM/ CFORM
   COMMON /RCTSV/ JDUM,SDUM,JFLG,RNJDM,NB,IPC,VLFT,VCNT,VRGT,        &
   &               WGTL,WGTR
!
   CHARACTER*12 BCD,HTRANS,HABSRB,HRADIA
   CHARACTER*11 CFORM
   CHARACTER*8 HSCNID(0:6)
   CHARACTER*18 HNAMPST,HVRPST
   CHARACTER h_ifil*7,SCNOUT*7,SCNINF*7,CTAPE*4
   LOGICAL OP
!
   DIMENSION FILHDR(2),SUMR(4)
   DIMENSION HWJ(0:6),DXJ(0:6),NJ(0:6),NJMX(0:6),SMPLJ(0:6),         &
   &          XSCAL(0:6)
!
   EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(6),ISCHDR) , (FSCDID(12),XSCID),              &
   &            (FSCDID(13),XHWHM) , (FSCDID(14),IDABS),              &
   &            (FSCDID(16),LAYR1)
!
   DATA I_1/1/, I_1000/1000/
!
   DATA HTRANS / 'TRANSMISSION'/,HABSRB / ' ABSORPTION '/,           &
   &     HRADIA / ' RADIANCE   '/
   DATA SCNOUT / '       '/,SCNINF / 'SCNINTF'/,CTAPE / 'TAPE'/
!
   DATA HSCNID(0) / 'RECTANGL'/,HWJ(0) / 1.         /,               &
   &     DXJ(0) / 0.0  /,NJ(0) / 0    /,NJMX(0) / 0    /,             &
   &     SMPLJ(0) / .5 /,XSCAL(0) / 0.          /
   DATA HSCNID(1) / 'TRIANGLE'/,HWJ(1) / 2.         /,               &
   &     DXJ(1) / 0.02 /,NJ(1) / 101  /,NJMX(1) / 251  /,             &
   &     SMPLJ(1) / 2. /,XSCAL(1) / 0.          /
   DATA HSCNID(2) / 'GAUSS   '/,HWJ(2) / 4.         /,               &
   &     DXJ(2) / 0.02 /,NJ(2) / 201  /,NJMX(2) / 251  /,             &
   &     SMPLJ(2) / 4. /,XSCAL(2) / 0.          /
!
!     SINCSQ: 54.18 HALFWIDTHS CORRESPONDS TO 24 ZERO CROSSINGS
!             PI CORRESPONDS TO X=2.257609141
!
   DATA HSCNID(3) / 'SINCSQ  '/,HWJ(3) / 54.1826    /,               &
   &     DXJ(3) / 0.02 /,NJ(3) / 2710 /,NJMX(3) / 2760 /,             &
   &     SMPLJ(3) / 4. /,XSCAL(3) / 1.391557377 /
!
!     SINC: 119.33 HALFWIDTHS CORRESPONDS TO 72 ZERO CROSSINGS
!           PI CORRESPONDS TO X=1.657400255
!
   DATA HSCNID(4) / 'SINC    '/,HWJ(4) / 119.332818 /,               &
   &     DXJ(4) / 0.02 /,NJ(4) / 5968 /,NJMX(4) / 6018 /,             &
   &     SMPLJ(4) / 4. /,XSCAL(4) / 1.89549425  /
   DATA HSCNID(5) / 'VRCTCENT'/,HWJ(5) / 1.         /,               &
   &     DXJ(5) / 0.0  /,NJ(5) / 0    /,NJMX(5) / 0    /,             &
   &     SMPLJ(5) / .5 /,XSCAL(5) / 0.          /
   DATA HSCNID(6) / 'VRCTLEFT'/,HWJ(6) / 1.         /,               &
   &     DXJ(6) / 0.0  /,NJ(6) / 0    /,NJMX(6) / 0    /,             &
   &     SMPLJ(6) / .5 /,XSCAL(6) / 0.          /
!
!----------------------------------------------------------------------
!
!    ADDITIONAL SCANNING FUNCTIONS MAY READILY BE ADDED TO THOSE
!      CURRENTLY IMPLEMENTED IN THIS VERSION OF LBLRTM:
!
!    A SHAPE SUBROUTINE FOR THE DESIRED FUNCTION MUST BE CREATED-
!     THIS SUBROUTINE PRECALCULATES THE FUNCTION FOR SUBSEQUENT
!      LOOKUP.  SEE FOR EXAMPLE SUBROUTINE SHAPEG FOR THE GAUSSIAN
!
!    THE SHAPE SUBROUTINE SETS UP THE SYMMETRIC FUNCTION IN ARRAY FG
!     AT EQUAL INCREMENTS OF THE HALFWIDTH, 'DXF'. THE VALUE OF 'DXF'
!     IS SET IN THIS SUBROUTINE BY THE VALUE OF 'DXJ(?)'
!
!    A DATA CARD MUST BE CREATED FOR EACH SCANNING FUNCTION DEFINING
!    THE FOLLOWING QUANTITIES:
!
!    HWJ(?)   EXTENT OF THE FUNCTION (BOUND) FROM THE CENTER IN UNITS
!               OF HALFWIDTH
!
!    DXJ(?)   INCREMENT AT WHICH THE FUNCTION IS STORED IN UNITS
!               OF HALFWIDTH
!
!    NJ(?)    THE NUMBER OF POINTS FROM THE CENTER TO THE FUNCTION
!               BOUND
!
!    NJMAX(?) SIZE OF THE ARRAY IN WHICH THE FUNCTION IS STORED
!               FUNCTION VALUES BETWEEN NJ AND NJMAX ARE ZERO
!
!    SMPL(?)  DEFAULT VALUE OF THE SAMPLING INCREMENT IN RECIPRICAL
!               HALFWIDTH UNITS: E.G. A VALUE OF FOUR MEANS THAT THE
!               OUTPUT SPACING, 'DV', IN WAVENUMBERS WILL BE 1/4 THE
!               HALFWIDTH VALUE IN WAVENUMBERS, 'HWHM'.
!
!    XSCAL(?) REQUIRED FOR PERIODIC FUNTIONS. THE VALUE OF THE
!               FUNCTION ARGUMENT IN RADIANS FOR WHICH THE
!               FUNCTION VALUE IS 0.5, E.G.
!                   SINX/X = 0.5 FOR X = 1.89549425, XSCAL(4)
!
!    CONSIDERATION MUST BE GIVEN TO THE ISSUE OF FUNCTION
!      NORMALIZATION FOR FUNCTIONS THAT DO NOT HAVE RAPID
!      CONVERGENCE TO ZERO (SINX/X)
!
!                                                               SAC
!
!----------------------------------------------------------------------
!
!
!     ASSIGN CVS VERSION NUMBER TO MODULE
!
   HVRPST = '$Revision$'
!
   PI = 2.*ASIN(1.)
!
!  SET THE MAXIMIM NUMBER OF AVAILABLE FUNCTIONS:
!
   NFNMAX = 6
!
!  NLIMF IS ONE MORE THAN THE SIZE OF OUTPUT (CONVOLVED) ARRAY
!
   NLIMF = 2401
   NREN = 0
   NSHIFT = 32
   IFLSAV = 0
   IPRT = 1
!
10 CONTINUE
   SUMOUT = 0.
   SMIN = 999999.
   SMAX = -99999.
   DVOSAV = 0.
   SUMR(1) = SUMOUT
   SUMR(2) = SMIN
   SUMR(3) = SMAX
   SUMR(4) = DVOSAV
!
   IEOFT = 1
!
!     READ IN CONTROL PARAMETERS: SEE INSTRUCTIONS FO DEFINITIONS
!
   READ (IRD,900,END=80) HWHM,V1,V2,JEMIT,JFN,JVAR,SAMPL,IUNIT,      &
   &                      IFILST,NIFILS,JUNIT,NPTS
!
   IF (HWHM.LE.0.) GO TO 70
!
!     JEMIT=-1   SCANFN CONVOLVED WITH ABSORPTION
!     JEMIT=0    SCANFN CONVOLVED WITH TRANSMISSION
!     JEMIT=1    SCANFN CONVOLVED WITH EMISSION
!
   JABS = 0
   IDABS = 0
   IF (JEMIT.LT.0) THEN
      JABS = 1
      JEMIT = 0
      IDABS = -1
   ENDIF
   IDABST = IDABS
!
!     JVAR=1 FOR A VARIABLE SLIT FUNCTION (NOT FOR JFN=0)
!     THE CODING IN CNVSCN  RESULTS IN HWHM=1./ (VI-V1)**2
!     HWHM IS CONSTANT FOR EACH PANEL AS PROGRAMMED
!
   IFN = ABS(JFN)
   IF (IFN.GT.NFNMAX) THEN
      WRITE (IPR,*) 'SCANF; JFN GT LIMIT'
      STOP 'SCANF; JFN GT LIMIT'
   ENDIF
!
   READ (HSCNID(IFN),905) SCANID
!
   HWF = HWJ(IFN)
   DXF = DXJ(IFN)
   NF = NJ(IFN)
   NFMAX = NJMX(IFN)
   SAMPLE = SMPLJ(IFN)
   XSCALE = XSCAL(IFN)
!
!     CHECK FOR NEGATIVE JFN OR NEGATIVE SAMPL
!
!     FOR NEGATIVE JFN, USER IS SUPPLYING FIRST ZERO CROSSING FOR THE
!     PERIODIC FUNCTION IN HWHM.  SET HWHM=(FIRST ZERO)/(PI/XSCALE)
!
!     For JFN=5,6 user is supplying instrument field of view half angle
!     in degrees in HWHM.  Trap if JFN=-5,-6.
!
   IF (JFN.LT.0) THEN
      JFN = ABS(JFN)
      IF ((JFN.EQ.3).OR.(JFN.EQ.4)) THEN
         HWHM = HWHM/(PI/XSCALE)
      ELSE
         WRITE (IPR,910) JFN
         STOP 'SCANFN; INVALID JFN'
      ENDIF
   ENDIF
!
!     SET DVINT TO DETERMINE IF INTERPOLATION IS NECESSARY
!     - For JFN = 5,6, set DVINT to 1/12 the width of the first box.
!       HWHM should carry the value of the field of view half angle
!       (in degrees).  This is converted to radians.  The box width
!       formula is
!
!                  width = V1*(1/2 angle FOV)**2/2
!
!       and the degrees-to-radians formula is
!
!                  rad = deg*3.141592654/180.
!
!     - For JFN not equal to 5 or 6, set DVINT to 1/12 the value of
!       HWHM.  HWHM should carry the true value of the Half Width
!       at Half Maximum of the scanning function at this point.
!
   IF ((JFN.EQ.5).OR.(JFN.EQ.6)) THEN
      DVINT = V1*(HWHM*3.141592654/180.)**2/24
   ELSE
      DVINT = HWHM/12.
   ENDIF
!
!     - For positive SAMPL, set SAMPLE equal to SAMPL (the number
!       of points per half width).
!     - For negative SAMPL, user is supplying desired DELVO
!       (outgoing spectral spacing).  SAMPLE (the number of sample
!       points per half width) is set such that SAMPLE=HWHM/DELVO
!       (Half Width at Half Max over user input outgoing spectral
!       spacing), and the outgoing spectral spacing DVO will be
!       recalculated using HWHM and SAMPLE below.
!
   IF (SAMPL.LT.0.) THEN
      SAMPLE = HWHM/(-SAMPL)
   ELSEIF (SAMPL.GT.0.) THEN
      SAMPLE = SAMPL
   ENDIF
!
!     SET UP SELECTED SCANNING FUNCTION:
!
   IF (JFN.EQ.1) CALL SHAPET (XF)
   IF (JFN.EQ.2) CALL SHAPEG (XF)
   IF (JFN.EQ.3) CALL SINCSQ (XF,XSCALE)
   IF (JFN.EQ.4) CALL SINC (XF,XSCALE)
!
   IF (IUNIT.LE.0) IUNIT = IFILE
   IFILE = IUNIT
   IFILST = MAX(IFILST,I_1)
   IF (NIFILS.LE.0) NIFILS = 99
!
!     SKIP TO SELECTED 'FILE'
!
   inquire(ifile,opened=op)
   IF (.NOT.OP) THEN
      WRITE (h_ifil,940) CTAPE,IFILE
      OPEN (IFILE,FILE=h_ifil,STATUS='UNKNOWN',FORM=CFORM)
   ENDIF

   REWIND IFILE
   IF (IFILST.GT.1) CALL SKIPFL (IFILST-1,IFILE,IEOF)
!
!     READ FILE HEADER FOR SELECTED 'FILE'
!
20 CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)
   IF (IEOF.EQ.0) GO TO 10
   IDABS = IDABST
!
   WRITE (IPR,915) XID,(YID(M),M=1,2)
   WRITE (IPR,920) LAYR1,LAYER
   WRITE (IPR,925) SECANT,PAVE,TAVE,DV,V1C,V2C
   WRITE (IPR,930) WBROAD,(HMOLID(M),WK(M),M=1,NMOL)
!
!     CHECK FOR INTERPOLATION AND OPEN OUTPUT FILE IF NECESSARY
!
!     IFILE INTERPOLATED ONTO JFILE
!
   IF (JUNIT.LE.0) JUNIT = JFILE
   JFILE = JUNIT
!
!     IF DV NOT FINE ENOUGH, FIRST INTERPOLATE
!
   IF (DV.GT.DVINT) THEN
      IFLSAV = IFILE
      JFLSAV = JFILE
      IEOFSC = 1
      JFILE = 77
      INQUIRE (UNIT=JFILE,OPENED=OP)
      IF (OP) CLOSE (JFILE)
      SCNOUT = SCNINF
      OPEN (JFILE,FILE=SCNOUT,STATUS='UNKNOWN',FORM=CFORM)
      REWIND JFILE
      IBUF = 0
!
!     INTERPOLATE:
!
      CALL SCNINT (IFILE,JFILE,DVINT,JEMIT,NPTS,IBUF)
!
      IEOFSV = IEOFSC
      WRITE (IPR,935)
      IFILE = JFILE
      REWIND IFILE
      CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)
      JFILE = JFLSAV
   ELSE
      IFLSAV = 0
   ENDIF
   INQUIRE (UNIT=JFILE,OPENED=OP)
   IF (.NOT.OP) THEN
      WRITE (SCNOUT,940) CTAPE,JFILE
      OPEN (JFILE,FILE=SCNOUT,STATUS='UNKNOWN',FORM=CFORM)
      REWIND JFILE
   ENDIF
!
   ISCAN = ISCHDR
   IF (ISCAN.LE.0.OR.XSCID.EQ.-99.) ISCAN = 0
   IF (ISCHDR.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCHDR
   ISCHDR = ISCAN+1
   JTREM = -1
   IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1
   ISCANT = MOD(ISCAN,I_1000)
   IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 2
   IF (JTREM.LT.0) THEN
      WRITE(IPR,*) ' SCANF; JTREM LT 0'
      STOP ' SCANF; JTREM LT 0'
   ENDIF
   WRITE (IPR,945) SCANID,IFILE,IFILST,NIFILS,JEMIT,JFN,JVAR,JABS
!
!     JTREM=0   SCANFN CONVOLVED WITH EXPONENTIATED
!                      ABSORPTION COEFFICIENT
!     JTREM=1   SCANFN CONVOLVED WITH EMISSION
!     JTREM=2   SCANFN CONVOLVED WITH TRANSMISSION
!
   DVI = DV
   DVSAV = DVI
!
!     Compute output spectral spacing.  For JFN not 5 or 6, at this
!     point HWHM always contains the value of the Half Width at Half
!     Maximum of the scanning function, and SAMPLE always contains
!     the number of points per half width of the scanning function.
!
!     For JFN = 5,6 at this point, HWHM contains the value of the
!     field of view half angle (in degrees), and SAMPLE contains
!     the ratio of the field of view half angle to the specified
!     output spectral spacing (the quotient of HWHM and SAMPLE
!     results in the circuitous calculation of the previously input
!     DVO).
!
   DVO = HWHM/SAMPLE
   IF (JFN.EQ.0) THEN
      IRATIO = DVO/DVI+0.5
      DVO = REAL(IRATIO)*DVI
      IF (IRATIO.LT.2) THEN
         WRITE (IPR,950)
         GO TO 10
      ENDIF
   ENDIF
!
!     BOUND AT THIS POINT IS THE WAVENUMBER VALUE
!     OF HALF THE SCANNING FUNCTION
!
   BOUND = HWF*HWHM
   DV = DVO
   V1C = V1
   V2C = V2
   SCIND = JVAR+10*(JFN+10*(JEMIT))
   XSCID = SCIND+0.01
   XHWHM = HWHM
   CALL BUFOUT (JFILE,FILHDR(1),NFHDRF)
   WRITE (IPR,955) HWHM,BOUND,JFILE,V1,V2,DVO
   NBOUND = (2.*HWF)*SAMPLE+0.01
!
!     NBOUND IS THE NUMBER OF SPECTRAL VALUES SPANNED
!     BY THE FULL SCANNING FUNCTION
!
!     RESET BOUND BASED ON NBOUND
!
   BOUND =  REAL(NBOUND)*DVO/2.
   MAXF = NLIMF+2*NBOUND+NSHIFT
!
   TIMRDF = 0.
   TIMCNV = 0.
   TIMPNL = 0.
   IEOFSC = 1
   NLO = NSHIFT+1
   SUMIN = 0.
   NHI = NLIMF+NSHIFT-1
   DO 30 I = 1, MAXF
      N1(I) = 0.
      R1(I) = 0.
30 END DO
   INIT = 0
   IDATA = -1
   IPANEL = -1
   JFLG = -1
   VFT = V1- REAL(NSHIFT)*DV
   VBOT = V1-BOUND
   VTOP = V2+BOUND
!
   IF (JEMIT.EQ.0.AND.IDABS.EQ.0) BCD = HTRANS
   IF (JEMIT.EQ.0.AND.IDABS.EQ.-1) BCD = HABSRB
   IF (JEMIT.EQ.1) BCD = HRADIA
   IF (NPTS.GT.0) WRITE (IPR,960) BCD
!
40 CALL CPUTIM (TIME0)
!
   IF (IEOFSC.LE.0) GO TO 60
!
!     READ DATA TO BE CONVOLVE FROM IFILE  AND PUT INTO ARRAY S
!
   CALL RDSCAN (S,JTREM,IFILE,ISCAN,IPRT)
!
!PRT  WRITE(IPR,965) IEOFSC,IDATA
!
   CALL CPUTIM (TIME)
   TIMRDF = TIMRDF+TIME-TIME0
!
   IF (IEOFSC.LE.0) GO TO 60
!
!     SHRKSC MAY SHRINK (COMPRESS) THE DATA; DVI IS MODIFIED ACCORDINGL
!
   IF ((JFN.NE.0).AND.(JFN.NE.5).AND.(JFN.NE.6)) THEN
      CALL SHRKSC (INIT,HWHM)
   ENDIF
!
50 CONTINUE
!
!     PERFORM THE CONVOLUTION OF XF ON S TO GIVE R1
!
   IF (JFN.EQ.0) THEN
      CALL CNVRCT (S,HWHM,R1,XF)
   ELSEIF (JFN.EQ.5) THEN
      CALL CNVVRC (S,HWHM,R1,XF)
   ELSEIF (JFN.EQ.6) THEN
      CALL CNVVRL (S,HWHM,R1,XF)
   ELSE
      CALL CONVSC (S,HWHM,R1,XF)
   ENDIF
!
!PRT  WRITE(IPR,965) IEOFSC,IDATA,IPANEL
!
   IF (IPANEL.EQ.0) GO TO 40
!
60 CONTINUE
!
!     OUTPUT PANEL TO JFILE, NPTS VALUES OF R1
!
   IF (JFN.EQ.0.OR.JFN.EQ.5.OR.JFN.EQ.6) THEN
      CALL PNLRCT (R1,JFILE,SUMR,NPTS)
   ELSE
      CALL PANLSC (R1,JFILE,SUMR,NPTS)
   ENDIF
!
   IF ((ISTOP.NE.1).AND.(IEOFSC.LT.0)) GO TO 60
   IF ((ISTOP.NE.1).AND.(IEOFSC.GT.0)) GO TO 50
   CALL CPUTIM (TIME)
   WRITE (IPR,970) TIME,TIMRDF,TIMCNV,TIMPNL
   CALL ENDFIL (JFILE)
!
   SUMIN = SUMIN*DVSAV
!
   WRITE (IPR,975) SUMIN
!
   IF (IFLSAV.NE.0) THEN
      IFILE = IFLSAV
      IEOFSC = IEOFSV
   ENDIF
   IF (IEOFSC.EQ.1) CALL SKIPFL (1,IFILE,IEOFSC)
!
   IEOFT = IEOFT+1
!
   SUMOUT = SUMR(1)
   SMIN = SUMR(2)
   SMAX = SUMR(3)
   DVOSAV = SUMR(4)
!
   SUMOUT = SUMOUT*DVOSAV
   WRITE (IPR,980) SUMOUT,SMIN,SMAX
!
   IF (IEOFT.LE.NIFILS.AND.IEOFSC.LT.0) GO TO 20
!
   GO TO 10
!
70 CONTINUE
!
80 RETURN
!
900 FORMAT (3F10.3,3(3X,I2),F10.4,4(3X,I2),I5)
905 FORMAT (A8)
910 FORMAT (//,' *****  INVALID VALUE FOR JFN = ',I2,'  *****',/)
915 FORMAT ('1',' **SCANFN** ',/,'0',10A8,2X,2(1X,A8,1X))
920 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)
925 FORMAT ('0 SECANT =',F15.5,/'0 PRESS(MB) =',F12.5/'0 TEMP(K) =',  &
   &        F11.2,/'0 DV(CM-1) = ',F12.8,/'0 V1(CM-1) = ',F12.6,/     &
   &        '0 V2(CM-1) = ',F12.6)
930 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//5X,'WBROAD = ',     &
   &        1PE10.3,/(5X,A6,' = ',1PE10.3))
935 FORMAT ('0',' **SCANFN** ',/)
940 FORMAT (A4,I2.2)
945 FORMAT ('0','***',A8,'***',//6X,'INPUT FILE NUMBER =',I3,         &
   &        ' ,IFILST = ',I5,' ,NIFILS = ',I5,',JEMIT =',I2,          &
   &        ' ,JFN =',I2,' ,JVAR =',I2,'  ,JABS =',I2)
950 FORMAT ('0',60X,'****** IRATIO LESS THAN 2, NO SCANFN ******')
955 FORMAT (1X,'     HWHM OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,   &
   &        5X,'BOUND OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,6X,    &
   &        'OUTPUT FILE NUMBER =',I3,',   V1 =',F12.5,',   V2 =',    &
   &        F12.5,5X,' DV OUT',F12.8)
960 FORMAT (///,'0',5X,A12,/)
965 FORMAT ('0',5X,'IEOFSC =',I3,'  IDATA =',I3,'  IPANEL =',I3,/)
970 FORMAT ('0',5X,'TIME =',F7.3,',  READ =',F6.3,',  CONV. =',F7.3,  &
   &        ',  PANEL =',F6.3)
975 FORMAT ('0    SUMIN  =',1P,E16.9)
980 FORMAT ('0    SUMOUT =',1P,E16.9,'  MIN =',E16.9,'  MAX =',E16.9)
!
end subroutine SCANFN
!
!     --------------------------------------------------------------
!
SUBROUTINE SCANRD (DVINT,IEMIT,iaj)
!
   !
   IMPLICIT REAL*8          (V)
!
!     READ CONTROL CARD FOR SCANNING WITH WEIGHTING FUNCTIONS
!
   COMMON S(3850),R1(5000)
!
   character*8      XID,       HMOLID,      YID,SCANID
   real*8               SECANT,       XALTZ
!

   common /scanaj/HWHMa,V1a,V2a,JEMITa,JFNa,JVARa,                   &
   &    SAMPLa,NNFILEa,NPTSa


   COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /RSCAN/ V1I,V2I,DVI,NNI
   COMMON /SCSHAP/ HWFS,DXFS,NFS,NFMAXS
   COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,             &
   &                HWF3,DXF3,NX3,N3MAX
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /SCINF/ HWHM,JEMIT,JFN,SAMPLE,SCANID,NPTS,XF(6018)
   COMMON /FLFORM/ CFORM
!
   CHARACTER*8 HSCNID(0:6)
   CHARACTER CFORM*11,TAPE13*6,CTAPE*4
   LOGICAL OP
!
   DIMENSION FILHDR(2)
   DIMENSION HWJ(0:6),DXJ(0:6),NJ(0:6),NJMX(0:6),SMPLJ(0:6),         &
   &          XSCAL(0:6)
!
   EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(6),ISCHDR),              &
   &            (FSCDID(12),XSCID) , (FSCDID(13),XHWHM),              &
   &            (FSCDID(14),IDABS) , (FSCDID(16),LAYR1)
!
   DATA HSCNID(0) / 'RECTANGL'/,HWJ(0) / 1.         /,               &
   &     DXJ(0) / 0.0  /,NJ(0) / 0    /,NJMX(0) / 0    /,             &
   &     SMPLJ(0) / .5 /,XSCAL(0) / 0.          /
   DATA HSCNID(1) / 'TRIANGLE'/,HWJ(1) / 2.         /,               &
   &     DXJ(1) / 0.02 /,NJ(1) / 101  /,NJMX(1) / 251  /,             &
   &     SMPLJ(1) / 2. /,XSCAL(1) / 0.          /
   DATA HSCNID(2) / 'GAUSS   '/,HWJ(2) / 4.         /,               &
   &     DXJ(2) / 0.02 /,NJ(2) / 201  /,NJMX(2) / 251  /,             &
   &     SMPLJ(2) / 4. /,XSCAL(2) / 0.          /
!
!     SINCSQ: 54.18 HALFWIDTHS CORRESPONDS TO 24 ZERO CROSSINGS
!             PI CORRESPONDS TO X=2.257609141
!
   DATA HSCNID(3) / 'SINCSQ  '/,HWJ(3) / 54.1826    /,               &
   &     DXJ(3) / 0.02 /,NJ(3) / 2710 /,NJMX(3) / 2760 /,             &
   &     SMPLJ(3) / 4. /,XSCAL(3) / 1.391557377 /
!
!     SINC: 119.33 HALFWIDTHS CORRESPONDS TO 72 ZERO CROSSINGS
!           PI CORRESPONDS TO X=1.657400255
!
   DATA HSCNID(4) / 'SINC    '/,HWJ(4) / 119.332818 /,               &
   &     DXJ(4) / 0.02 /,NJ(4) / 5968 /,NJMX(4) / 6018 /,             &
   &     SMPLJ(4) / 4. /,XSCAL(4) / 1.89549425  /
   DATA HSCNID(5) / 'VRCTCENT'/,HWJ(5) / 1.         /,               &
   &     DXJ(5) / 0.0  /,NJ(5) / 0    /,NJMX(5) / 0    /,             &
   &     SMPLJ(5) / .5 /,XSCAL(5) / 0.          /
   DATA HSCNID(6) / 'VRCTLEFT'/,HWJ(6) / 1.         /,               &
   &     DXJ(6) / 0.0  /,NJ(6) / 0    /,NJMX(6) / 0    /,             &
   &     SMPLJ(6) / .5 /,XSCAL(6) / 0.          /
!
   DATA TAPE13 / '      '/,CTAPE / 'TAPE'/
!
   PI = 2.*ASIN(1.)
!
!  SET THE MAXIMIM NUMBER OF AVAILABLE FUNCTIONS:
!
   NFNMAX = 6
!
   NLIMF = 2401
   NSHIFT = 32

   if (iaj.eq.0) then
      READ (IRD,900,END=10) HWHM,V1,V2,JEMIT,JFN,JVAR,              &
      &        SAMPL,NNFILE,NPTS
   else
      HWHM=HWHMa
      V1=V1a
      V2=V2a
      JEMIT=JEMITa
      JFN=JFNa
      JVAR=JVARa
      SAMPL=SAMPLa
      NNFILE=NNFILEa
      NPTS=NPTSa
   endif



!
   IF (HWHM.LE.0.) THEN
      WRITE(IPR,*) ' SCANRD * HWHM NEGATIVE '
      STOP ' SCANRD * HWHM NEGATIVE '
   ENDIF
!
!     JEMIT=-1   SCANFN CONVOLVED WITH ABSORPTION
!     JEMIT=0    SCANFN CONVOLVED WITH TRANSMISSION
!     JEMIT=1    SCANFN CONVOLVED WITH EMISSION
!
   JABS = 0
!
!     THE FOLLOWING CARDS HAVE BEEN RETRAINED
!     FOR POSSIBLE FUTURE CODE ENHANCEMENTS
!
!C    IF (JEMIT.LT.0) THEN
!C       JABS=1
!C       JEMIT=0
!C    ENDIF
!
!     JVAR=1 FOR A VARIABLE SLIT FUNCTION (NOT FOR JFN=0)
!     THE CODING IN CNVSCN  RESULTS IN HWHM=1./ (VI-V1)**2
!     HWHM IS CONSTANT FOR EACH PANEL AS PROGRAMMED
!     FOLLOWING VALUES INITIALIZE FOR RECTANGLE
!
   IFN = ABS(JFN)
   IF (IFN.GT.NFNMAX) THEN
      WRITE(IPR,*)' SCANF; JFN GT LIMIT'
      STOP ' SCANF; JFN GT LIMIT'
   ENDIF
!
   READ (HSCNID(IFN),905) SCANID
!
!     JVAR=1 FOR A VARIABLE SLIT FUNCTION (NOT FOR JFN=0)
!     THE CODING IN CNVSCN  RESULTS IN HWHM=1./ (VI-V1)**2
!     HWHM IS CONSTANT FOR EACH PANEL AS PROGRAMMED
!     FOLLOWING VALUES INITIALIZE FOR RECTANGLE
!
   HWF = HWJ(IFN)
   DXF = DXJ(IFN)
   NF = NJ(IFN)
   NFMAX = NJMX(IFN)
   SAMPLE = SMPLJ(IFN)
   XSCALE = XSCAL(IFN)
!
!     Set values of HWFS, DXFS, NFS, & NFMAXS to HWF, DXF, NF,
!     & NFMAX for use when entering HIRAC1 between SCANRD and
!     SCNMRG.
!
   HWFS = HWF
   DXFS = DXF
   NFS = NF
   NFMAXS = NFMAX
!
!     CHECK FOR NEGATIVE JFN OR NEGATIVE SAMPL
!
!     FOR NEGATIVE JFN, USER IS SUPPLYING FIRST ZERO CROSSING FOR THE
!     PERIODIC FUNCTION IN HWHM.  SET HWHM=(FIRST ZERO)/(PI/XSCALE)
!
!     For JFN=5,6 user is supplying instrument field of view half angle
!     in degrees in HWHM.
!
!     FOR NEGATIVE SAMPL, USER IS SUPPLYING DESIRED DELVO.
!     SET SAMPLE=HWHM/DELVO.
!
   IF (JFN.LT.0) THEN
      JFN = ABS(JFN)
      IF ((JFN.EQ.3).OR.(JFN.EQ.4)) THEN
         HWHM = HWHM/(PI/XSCALE)
      ELSE
         WRITE (IPR,910) JFN
         STOP 'SCANRD; INVALID JFN'
      ENDIF
   ENDIF
!
   IF (SAMPL.LT.0.) SAMPLE = HWHM/(-SAMPL)
   IF (SAMPL.GT.0.) SAMPLE = SAMPL
!
   IF (JFN.EQ.1) CALL SHAPET (XF)
   IF (JFN.EQ.2) CALL SHAPEG (XF)
   IF (JFN.EQ.3) CALL SINCSQ (XF,XSCALE)
   IF (JFN.EQ.4) CALL SINC (XF,XSCALE)
!
   IF (NNFILE.NE.NFILE.AND.NNFILE.GT.0) THEN
      INQUIRE (UNIT=NFILE,OPENED=OP)
      IF (OP) CLOSE (NFILE)
      NFILE = NNFILE
      INQUIRE (UNIT=NFILE,OPENED=OP)
      IF (.NOT.OP) THEN
         WRITE (TAPE13,915) CTAPE,NFILE
         OPEN (NFILE,FILE=TAPE13,STATUS='UNKNOWN',FORM=CFORM)
         REWIND NFILE
      ENDIF
   ENDIF
!
   WRITE (IPR,920) SCANID,JEMIT,JFN,JVAR,SAMPL,NPTS
   IEMIT = JEMIT
!
!    BOUND AT THIS POINT IS THE WAVENUMBER VALUE
!    OF HALF THE SCANNING FUNCTION
!
   DVO = HWHM/SAMPLE
   DVINT = HWHM/12.
   BOUND = HWF*HWHM
   V1C = V1
   V2C = V2
   XHWHM = HWHM
   WRITE (IPR,925) HWHM,BOUND,NFILE,V1,V2
   RETURN
10 PRINT 930
   STOP
!
900 FORMAT (3F10.3,3(3X,I2),F10.4,15X,2I5)
905 FORMAT (A8)
910 FORMAT (//,' *****  INVALID VALUE FOR JFN = ',I2,'  *****',/)
915 FORMAT (A4,I2.2)
920 FORMAT ('1',5X,'SCANRD',5X,'***',A8,'***',/,/,' JEMIT =',I2,      &
   &        ' JFN =',I2,' ,JVAR =',I2,' ,SAMPL =',F10.4,'  ,NPTS =',  &
   &        I5)
925 FORMAT (1X,'     HWHM OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,   &
   &        5X,'BOUND OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,6X,    &
   &        'OUTPUT FILE NUMBER =',I3,',   V1 =',F12.5,',   V2 =',    &
   &        F12.5)
930 FORMAT (' END OF FILE TAPE5',/,' (NOTE TAPE10 ALREADY CREATED )')
!
end subroutine SCANRD
!
!     --------------------------------------------------------------
!
subroutine scanrd_aj

! This subroutine is used with the analytic jacobian calculations
! to store the TAPE5 scan information rather than reading it in for
! each layer/level to be scanned.  It also fixes problems with the
! duplicate use of common blocks for scanning the cross-sectional
! molecules to the correct resolution.

   implicit real*8 (v)

   common /scanaj/HWHM,V1,V2,JEMIT,JFN,JVAR,SAMPL,NNFILE,NPTS

   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4

   READ (IRD,900,END=10) HWHM,V1,V2,JEMIT,JFN,JVAR,SAMPL,NNFILE,NPTS
900 FORMAT (3F10.3,3(3X,I2),F10.4,15X,2I5)

   return

10 print 930
930 FORMAT (' END OF FILE TAPE5 in scanrd_aj')
end subroutine scanrd_aj
!
!     --------------------------------------------------------------
!
SUBROUTINE SCNINT (IFILE,JFILE,DVINT,JEMIT,NPTS,IBUF)
!
   !
   IMPLICIT REAL*8          (V)
!
!**********************************************************************
!
!     INTERPOLATION FUNCTION DRIVER FOR WEIGHTING FUNCTIONS
!
!     FOUR-POINT VERSION    (MARCH 1990)
!
!**********************************************************************
!
!     THE INPUT DATA WILL BE PUT INTO T(5) = S(1) WITH THE LAST
!     4 POINTS OF THE PREVIOUS PANEL PUT INTO T(1 TO 4).
!     THIS SCHEME PERMITS 6 POINT INTERPOLATION.
!
!     S IS NOMINALLY 2401 POINTS BUT MAY NEED TO BE EXTENDED BY
!     2 POINTS TO PERMIT 4 POINT INTERPOLATION UP TO THE LAST
!     DATA POINT.
!
   COMMON T(2410),R(2401)
   DIMENSION S(2406)
   EQUIVALENCE (T(5),S(1))
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WN2   ,DV ,V1C,V2C,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /INPNL/ V1I,V2I,DVI,NNI
   COMMON /OUTPNL/ V1J,V2J,DVJ,NNJ
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
   DIMENSION FILHDR(2),RSTAT(3)
!
   EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(6),ISCHDR) , (FSCDID(12),XSCID),              &
   &            (FSCDID(13),XHWHM) , (FSCDID(14),IDABS),              &
   &            (FSCDID(16),LAYR1)
!
   CHARACTER*12 BCD,HTRANS,HABSRB,HRADIA
!
   DATA I_1000/1000/
!
   DATA HTRANS / 'TRANSMISSION'/,HABSRB / ' ABSORPTION '/,           &
   &     HRADIA / ' RADIANCE   '/
!
!----------------------------------------------------------------------
!     JEMIT=-1  INTERPOLATE ABSORPTION
!     JEMIT=0   INTERPOLATE TRANSMISSION
!     JEMIT=1   INTERPOLATE EMISSION
!     JEMIT=2   INTERPOLATE OPTICAL DEPTH
!----------------------------------------------------------------------
!
   WRITE (IPR,900)
   CALL CPUTIM (TIME1)
   TIMRDF = 0.0
   TIMCNV = 0.0
   TIMPNL = 0.0
!
   V1SAV = V1
   V2SAV = V2
   DVSAV = DV
   DVOSAV = 0.
!
   DVO = DVINT
!
!     I4PT = 1 FOR FOUR POINT INTERPOLATION
!
   I4PT = 1
   ICNVRT = 1
   IF (DVO.LE.0.) GO TO 40
!
   IF (IBUF.EQ.1) REWIND IFILE
   REWIND JFILE
!
!     BUFFER IN THE FILE HEADER ON UNIT (IFILE)
!     BUFFER OUT ON UNIT (JFILE)
!
   IF (IBUF.EQ.1) THEN
      CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)
      IF (IEOF.EQ.0) GO TO 30
      JABS = 0
      IDABS = 0
      IF (JEMIT.LT.0) THEN
         JABS = 1
         JEMIT = 0
         IDABS = -1
      ENDIF
   ENDIF
   V1 = V1C
   V2 = V2C
   DVI = DV
!
!   V2 IS ONLY APPROXIMATE
!
   NUM = (((V2-V1)/DVO)+0.5)
   V2 = V1+ REAL(NUM)*DVO
   NUM = NUM+1
   WRITE (IPR,905) V1,V2,DVO,NUM,JEMIT,I4PT,IFILE,JFILE,NPTS
!
   ISCAN = ISCHDR
   IF (ISCAN.LE.0.OR.XSCID.EQ.-99.) ISCAN = 0
   IF (ISCHDR.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCHDR
   ISCHDR = ISCAN+10
   V1C = V1
   V2C = V2
   DV = DVO
!
   SCNID = 100*JEMIT
   XSCID = SCNID+0.01
!
   CALL BUFOUT (JFILE,FILHDR(1),NFHDRF)
!
   JTREM = -1
   IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.2)) JTREM = 2
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1
   ISCANT = MOD(ISCAN,I_1000)
   IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 2
   WRITE (IPR,910) IFILE,IEMIT,JEMIT,JTREM,JABS
   IF (JTREM.LT.0) THEN
      WRITE(IPR,*) ' Invalid JTREM in SCNINT '
      STOP ' Invalid JTREM in SCNINT '
   ENDIF
!
   IDATA = -1
!
!     NEED TO SAVE LAST IBOUND POINTS OF EACH PANEL TO ATTACH TO NEXT
!
   IBOUND = 4
!
!     VBOT IS LOWEST NEEDED WAVENUMBER, VTOP IS HIGHEST
!
   BOUND =  REAL(IBOUND)*DV
   VBOT = V1-BOUND
   VTOP = V2+BOUND
!
   IF (JEMIT.EQ.0.AND.IDABS.EQ.0) BCD = HTRANS
   IF (JEMIT.EQ.0.AND.IDABS.EQ.-1) BCD = HABSRB
   IF (JEMIT.EQ.1) BCD = HRADIA
   IF (NPTS.GT.0) WRITE (IPR,915) BCD
!
!     ZERO OUT T(1 TO IBOUND)
!
   DO 10 II = 1, IBOUND
      T(II) = 0.0
10 END DO
!
!     READ FROM IFILE UNTIL THE FIRST REQUIRED POINT IS REACHED
!     AND LOAD DATA INTO S
!
   CALL RDPANL (S,JTREM,IFILE,ISCAN,JEMIT,ICNVRT)
   IF (IEOFSC.LE.0) GO TO 20
!
!     DO INTERPOLATION
!
   CALL INTERP (IFILE,JFILE,I4PT,IBOUND,NPTS,JTREM,ISCAN,JEMIT,      &
   &             RSTAT,ICNVRT)
!
   CALL CPUTIM (TIME2)
   CALL ENDFIL (JFILE)
!
!     WRITE STATISTICS
!
   WRITE (IPR,920) RSTAT(1),RSTAT(2),RSTAT(3)
   TIMTOT = TIME2-TIME1
   TIMCNV = TIMTOT-TIMRDF-TIMPNL
   WRITE (IPR,925) TIMTOT,TIMRDF,TIMCNV,TIMPNL
!
   GO TO 30
!
20 CONTINUE
   WRITE (IPR,930) IFILE
!
30 CONTINUE
   V1 = V1SAV
   V2 = V2SAV
   DV = DVSAV
   RETURN
!
40 CONTINUE
   WRITE (IPR,935) DVINT
!
   RETURN
!
900 FORMAT (/,'0***SCNINT***',/)
905 FORMAT (5X,'V1=',F14.8,' V2=',F14.8,' DVO=',E14.6,' NUM=',I8,/,   &
   &        5X,'JEMIT=',I3,' I4PT=',I3,' IUNIT=',I3,' JUNIT=',I3,     &
   &        ' NPTS=',I5)
910 FORMAT (5X,'INPUT FILE NUMBER =',I3,' IEMIT=',I3,' JEMIT=',I3,    &
   &        ' JTREM=',I3,' JABS=',I3)
915 FORMAT (///,'0',5X,A12,/)
920 FORMAT ('0    SUMOUT =',1P,E16.9,'  MIN =',E16.9,'  MAX =',E16.9)
925 FORMAT (/,5X,'SCNINT TIME: TOTAL = ',F8.3,' READ = ',F8.3,        &
   &        ' INTERP = ',F8.3,' WRITE = ',F8.3,/)
930 FORMAT (/,5X,'SCNINT- ERROR: EOF ON INPUT UNIT ',I4,              &
   &        ' BEFORE V1 WAS REACHED',/)
935 FORMAT (/,5X,'SCNINT- ERROR: DVINT .LT. ZERO ; DVINT =',F12.4,/)
!
end subroutine SCNINT
!
!     --------------------------------------------------------------
!
SUBROUTINE SCNMRG (IFILE,JFILE)
!
   !
   IMPLICIT REAL*8          (V)
!
!     DRIVER FOR CONVOLVING INSTRUMENTAL SCANNING FUNCTION
!     WITH SPECTRUM
!
   COMMON S(3850),R1(5000)
!
   character*8      XID,       HMOLID,      YID,SCANID
   real*8               SECANT,       XALTZ
!
   COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /RSCAN/ V1I,V2I,DVI,NNI
   COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,             &
   &                HWF3,DXF3,NX3,N3MAX
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /SCINF/ HWHM,JEMIT,JFN,SAMPLE,SCANID,NPTS,XF(6018)
   COMMON /RCTSV/ JJ,SUMJ,JFLG,RNJ,NB,IPC,VLFT,VCNT,VRGT,WGTL,WGTR
!
   DIMENSION FILHDR(2)
   DIMENSION SUMR(4)
!
   EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(6),ISCHDR) , (FSCDID(12),XSCID),              &
   &            (FSCDID(13),XHWHM) , (FSCDID(14),IDABS),              &
   &            (FSCDID(16),LAYR1)
!
   CHARACTER*12 BCD,HTRANS,HABSRB,HRADIA
!
   DATA I_1000/1000/
!
   DATA HTRANS / 'TRANSMISSION'/,HABSRB / ' ABSORPTION '/,           &
   &     HRADIA / ' RADIANCE   '/
!
!     IUNIT INPUT FILE
!     JUNIT OUTPUT FILE
!
   IUNIT = IFILE
   JUNIT = JFILE
   NREN = 0
   IPRT = 1
   IDABS = 0
   IF (JEMIT.LT.0) THEN
      JABS = 1
      JEMIT = 0
      IDABS = -1
   ENDIF
   IDABST = IDABS
   IFILST = 1
   NIFILS = 9999
!
   SUMOUT = 0.
   SMIN = 999999.
   SMAX = -99999.
   DVOSAV = 0.
   SUMR(1) = SUMOUT
   SUMR(2) = SMIN
   SUMR(3) = SMAX
   SUMR(4) = DVOSAV
   NSHIFT = 32
!
   REWIND IUNIT
   CALL BUFIN (IUNIT,IEOF,FILHDR(1),NFHDRF)
   IF (IEOF.EQ.0) GO TO 50
!
   DVSAV = DV
   IDABS = IDABST
!
   WRITE (IPR,900) XID,(YID(M),M=1,2)
   WRITE (IPR,905) LAYR1,LAYER
   WRITE (IPR,910) SECANT,PAVE,TAVE,DV,V1C,V2C
   WRITE (IPR,915) WBROAD,(HMOLID(M),WK(M),M=1,NMOL)
!
   ISCAN = ISCHDR
   IF (ISCAN.LE.0.OR.XSCID.EQ.-99.) ISCAN = 0
   IF (ISCHDR.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCHDR
   ISCHDR = ISCAN+1
   JTREM = -1
   IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1
   ISCANT = MOD(ISCAN,I_1000)
   IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 4
!
   WRITE (IPR,920) SCANID,IUNIT,IFILST,NIFILS,JEMIT,JFN,JVAR,JABS
!
   IF (JTREM.LT.0) THEN
      WRITE(IPR,*) ' SCANF; JTREM LT 0'
      STOP ' SCANF; JTREM LT 0'
   ENDIF
!
!     JTREM=0   SCANFN CONVOLVED WITH EXPONENTIATED
!                      ABSORPTION COEFFICIENT
!     JTREM=1   SCANFN CONVOLVED WITH EMISSION
!     JTREM=2   SCANFN CONVOLVED WITH TRANSMISSION
!
   DVI = DV
   DVO = HWHM/SAMPLE
!
!    BOUND AT THIS POINT IS THE WAVENUMBER VALUE
!    OF HALF THE SCANNING FUNCTION
!
   BOUND = HWF*HWHM
   DV = DVO
   V1C = V1
   V2C = V2
   SCIND = JVAR+10*(JFN+10*(JEMIT))
   XSCID = SCIND+0.01
   XHWHM = HWHM
   CALL BUFOUT (JUNIT,FILHDR(1),NFHDRF)
   WRITE (IPR,925) HWHM,BOUND,JUNIT,V1,V2,DVO
   NBOUND = (2.*HWF)*SAMPLE+0.01
!
!     BOUND AT THIS POINT IS THE WAVENUMBER VALUE OF THE
!     FULL SCANNING FUNCTION
!
   BOUND =  REAL(NBOUND)*DVO/2.
   MAXF = NLIMF+2*NBOUND+NSHIFT
!
   TIMRDF = 0.
   TIMCNV = 0.
   TIMPNL = 0.
   IEOFSC = 1
   NLO = NSHIFT+1
   SUMIN = 0.
   NHI = NLIMF+NSHIFT-1
   DO 10 I = 1, MAXF
      R1(I) = 0.
10 END DO
   INIT = 0
   IDATA = -1
   IPANEL = -1
   JFLG = -1
   VFT = V1- REAL(NSHIFT)*DV
   VBOT = V1-BOUND
   VTOP = V2+BOUND
!
   IF (JEMIT.EQ.0.AND.IDABS.EQ.0) BCD = HTRANS
   IF (JEMIT.EQ.0.AND.IDABS.EQ.-1) BCD = HABSRB
   IF (JEMIT.EQ.1) BCD = HRADIA
   IF (NPTS.GT.0) WRITE (IPR,930) BCD
20 CALL CPUTIM (TIME0)
   IF (IEOFSC.LE.0) GO TO 40
   CALL RDSCAN (S,JTREM,IUNIT,ISCAN,IPRT)
!
!PRT  WRITE(IPR,935) IEOFSC,IDATA
!
   CALL CPUTIM (TIME)
   TIMRDF = TIMRDF+TIME-TIME0
!
   IF (IEOFSC.LE.0) GO TO 40
   IF (JFN.NE.0) CALL SHRKSC (INIT,HWHM)
!
!     SHRKSC MAY SHRINK (COMPRESS) THE DATA;
!     DVI IS MODIFIED ACCORDINGLY
!
30 CONTINUE
   IF (JFN.EQ.0) THEN
      CALL CNVRCT (S,HWHM,R1,XF)
   ELSEIF (JFN.EQ.5) THEN
      CALL CNVVRC (S,HWHM,R1,XF)
   ELSEIF (JFN.EQ.6) THEN
      CALL CNVVRL (S,HWHM,R1,XF)
   ELSE
      CALL CONVSC (S,HWHM,R1,XF)
   ENDIF
!
!PRT  WRITE(IPR,935) IEOFSC,IDATA,IPANEL
!
   IF (IPANEL.EQ.0) GO TO 20
!
40 CONTINUE
   IF (JFN.EQ.0.OR.JFN.EQ.5.OR.JFN.EQ.6) THEN
      CALL PNLRCT (R1,JUNIT,SUMR,NPTS)
   ELSE
      CALL PANLSC (R1,JUNIT,SUMR,NPTS)
   ENDIF
   IF ((ISTOP.NE.1).AND.(IEOFSC.GT.0)) GO TO 30
   CALL CPUTIM (TIME)
   WRITE (IPR,940) TIME,TIMRDF,TIMCNV,TIMPNL
!
   SUMIN = SUMIN*DVSAV
!
   WRITE (IPR,945) SUMIN
!
   IF (IEOFSC.EQ.1) CALL SKIPFL (1,IUNIT,IEOFSC)
!
   IEOFT = IEOFT+1
!
!
   SUMOUT = SUMR(1)
   SMIN = SUMR(2)
   SMAX = SUMR(3)
   DVOSAV = SUMR(4)
!
   SUMOUT = SUMOUT*DVOSAV
   WRITE (IPR,950) SUMOUT,SMIN,SMAX
!
50 RETURN
!
900 FORMAT ('0',' **SCNMRG** ',/,'0',10A8,2X,2(1X,A8,1X))
905 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)
910 FORMAT ('0 SECANT =',F15.5,/'0 PRESS(MB) =',F12.5/'0 TEMP(K) =',  &
   &        F11.2,/'0 DV(CM-1) = ',F12.8,/'0 V1(CM-1) = ',F12.6,/,    &
   &        '0 V2(CM-1) = ',F12.6)
915 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//5X,'WBROAD = ',     &
   &        1PE10.3,/(5X,A6,' = ',1PE10.3))
920 FORMAT ('0','***',A8,'***',//6X,'INPUT FILE NUMBER =',I3,         &
   &        ' ,IFILST = ',I5,' ,NIFILS = ',I5,',JEMIT =',I2,          &
   &        ' ,JFN =',I2,' ,JVAR =',I2,'  ,JABS =',I2)
925 FORMAT (1X,'     HWHM OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,   &
   &        5X,'BOUND OF INSTRUMENT FUNCTION =',F12.8,' CM-1'/,6X,    &
   &        'OUTPUT FILE NUMBER =',I3,',   V1 =',F12.5,',   V2 =',    &
   &        F12.5,5X,' DV OUT',F12.8)
930 FORMAT (///,'0',5X,A12,/)
935 FORMAT ('0',5X,'IEOFSC =',I3,'  IDATA =',I3,'  IPANEL =',I3,/)
940 FORMAT ('0',5X,'TIME =',F7.3,',  READ =',F6.3,',  CONV. =',F7.3,  &
   &        ',  PANEL =',F6.3)
945 FORMAT ('0    SUMIN  =',1P,E16.9)
950 FORMAT ('0    SUMOUT =',1P,E16.9,'  MIN =',E16.9,'  MAX =',E16.9)
!
end subroutine SCNMRG
!
!     --------------------------------------------------------------
!
SUBROUTINE SHRKSC (INIT,HWHM)
!
   !
   IMPLICIT REAL*8          (V)
!
!     THIS SUBROUTINE COMPRESSES (SHRINKS) THE INPUT TO THE CONVOLUTION
!     ROUTINE FOR THE SCANNING FUNCTION TO ACCELERATE THE CALCULATION
!
   COMMON S(3850),R1(5000),SS(200)
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /RSCAN/ V1I,V2I,DVI,NLIM
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   DIMENSION JRATIO(24)
!
   DATA I_1/1/
!
   DATA JRATIO / 1,2,3,4,5,6,8,10,12,15,16,20,24,25,30,32,40,48,50,  &
   &             60,75,80,100,120 /
!
   CALL CPUTIM (TIME0)
   NLIMS = NLIM
   IF (NREN.GT.0) THEN
      DO 10 I = 1, NREN
         S(I) = SS(I)
10    CONTINUE
   ENDIF
   NREN = 0
   IF (INIT.EQ.0) THEN
      DVSC = HWHM/12.
      IRATSH = DVSC/DVI+0.5
      DO 20 I = 2, 24
         IF (JRATIO(I).GT.IRATSH) THEN
            IRATSH = JRATIO(I-1)
            GO TO 30
         ENDIF
20    CONTINUE
30    IF (IRATSH.GT.JRATIO(24)) IRATSH = JRATIO(24)
      IF (IRATSH.LE.1) RETURN
      DVSC = REAL(IRATSH)*DVI
      V1SHFT = REAL(IRATSH-1)*DVI/2.
      WRITE (IPR,900) IRATSH
      SRATIO = IRATSH
      IRATM1 = IRATSH-1
      INIT = 1
   ENDIF
   IF (IRATSH.LE.1) RETURN
   NREN = NLIM-(NLIM/IRATSH)*IRATSH
!
!PRT  WRITE(IPR,905) V1I,V1SHFT,DVSC,NREN
!
   V1I = V1I+V1SHFT
   IMIN = 1
   IMAX = NLIM-IRATM1-NREN
!
   K = 0
   DO 50 I = IMIN, IMAX, IRATSH
      SUMK = 0.
      JHI = I+IRATM1
      K = K+1
      DO 40 J = I, JHI
         SUMK = SUMK+S(J)
40    CONTINUE
      S(K) = SUMK/SRATIO
50 END DO
!
   V2I = V1I+DVSC* REAL(K-1)
   NLIM = K
   DVI = DVSC
   ILO = ((VBOT-V1I)/DVI)+1.5
   ILO = MAX(ILO,I_1)
   IHI = ((VTOP-V1I)/DVI)+1.5
   IHI = MIN(IHI,NLIM)
!
!PRT  WRITE(IPR,910) ILO,IHI
!
   IF (NREN.GT.0) THEN
      DO 60 I = 1, NREN
         II = NLIMS-NREN+I
         SS(I) = S(II)
60    CONTINUE
   ENDIF
   CALL CPUTIM (TIME)
   TIMCNV = TIMCNV+TIME-TIME0
!
   RETURN
!
900 FORMAT ('   SHRINK RATIO = ',I5)
905 FORMAT ('   V1I =',F10.3,'  V1SHFT =',F10.3,'  DVSC =',F12.5,     &
   &        '   NREN =',I4)
910 FORMAT ('   ILO =',I4,'  IHI =',I4)
!
end subroutine SHRKSC
!
!     --------------------------------------------------------------
!
SUBROUTINE SHAPET (XF)
!
!     SUBROUTINE SHAPET SETS UP THE TRIANGULAR SCANNING FUNCTION
!
   COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,             &
   &                HWF3,DXF3,NX3,N3MAX
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   DIMENSION XF(*)
!
   XTRIAN(X) = 1.-0.5*X
   DO 10 I = 1, NFMAX
      XF(I) = 0.
10 END DO
   XF(1) = 0.5
   SUM = XF(1)
   DO 20 I = 2, NF
      X = REAL(I-1)*DXF
      XF(I) = 0.5*XTRIAN(X)
      SUM = SUM+2.*XF(I)
20 END DO
   SUM = SUM*DXF
!
!PRT  WRITE(IPR,900) NF,DXF,SUM
!
   RETURN
!
900 FORMAT ('0',5X,'NF =',I5,',  DXF =',F7.5,',    SUM =',F18.15)
!
end subroutine SHAPET
!
!     --------------------------------------------------------------
!
SUBROUTINE RDSCAN (S,JTREM,IFILE,ISCAN,IPRT)
!
   !
   IMPLICIT REAL*8          (V)
!
!     SUBROUTINE RDSCAN INPUTS PANELS FROM IFILE RESULTING
!     FROM THE LBLRTM CALCULATION FOR CONVOLUTION
!     WITH THE SELECTED SCANNING FUNCTION
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /RSCAN/ VMIN,VMAX,DVI,NNI
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   DIMENSION DUMMY(2),PNLHDR(2)
   DIMENSION S(*)
!
   EQUIVALENCE (PNLHDR(1),VMIN)
!
   DATA I_1000/1000/
!
!PRT  WRITE(IPR,900) VBOT,VTOP
!
   IDUM1 = 0
   IDUM2 = 0
   ISCANT = MOD(ISCAN,I_1000)
   IF(JTREM.EQ.0.AND.ISCANT.GE.1) GO TO 60
   IF (ISCAN.LT.1) THEN
      IF (JTREM.EQ.1) IDUM1 = 1
      IF (JTREM.EQ.2) IDUM2 = 1
   ENDIF
10 CALL BUFIN (IFILE,IEOFSC,PNLHDR(1),NPHDRF)
   IF (IEOFSC.LE.0) GO TO 50
   NLOW = NREN+1
   IF (NREN.LE.0) NLOW = 1
   VMIN = VMIN-(NLOW-1)*DVI
   NNB = NNI
   NNI = NNI+NLOW-1
   IF ((IDATA.EQ.-1).AND.(VMIN.GT.VBOT).AND.(IPRT.EQ.1))             &
   &     WRITE (IPR,905)
   IDATA = 0
   IF (VMAX.GE.VBOT) GO TO 20
   IF (IDUM2.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),1)
   CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)
   IF (IDUM1.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),1)
   GO TO 10
20 IF (JTREM.EQ.0 .OR. JTREM.EQ.4 ) THEN
      CALL BUFIN (IFILE,IEOFSC,S(NLOW),NNB)
      DO 30 I = NLOW, NNI
         SI = S(I)
         S(I) = 1.
         IF (SI.GT.1.0E-04) THEN
            IF (SI.LT.ARGMIN) THEN
               S(I) = EXP(-SI)
            ELSE
               S(I) = EXPMIN
            ENDIF
         ELSE
            S(I) = 1.-SI
         ENDIF
30    CONTINUE
   ELSE
!
      IF (IDUM2.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),1)
      CALL BUFIN (IFILE,IEOFSC,S(NLOW),NNB)
      IF (IDUM1.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),1)
   ENDIF
!
!PRT  WRITE(IPR,910) VMIN,VMAX,DVI,NLOW,NNI
!
   IF (JABS.NE.0) THEN
      DO 40 I = NLOW, NNI
         S(I) = 1.-S(I)
40    CONTINUE
   ENDIF
   ILO = 1
   IHI = NNI
   DIF = (VMIN-VBOT)/DVI
   IF (DIF.LT.0.) ILO = -DIF+1.5
   IF (VMAX.LE.VTOP) RETURN
   IHI = (VTOP-VMIN)/DVI+1.5
   IDATA = 1
   RETURN
50 IF (IPRT.EQ.1) WRITE (IPR,915)
   RETURN
!
60 WRITE(IPR,920) JTREM,ISCAN
   RETURN
!
900 FORMAT ('0',/,'0   READING SPECTRUM, VBOT =',F10.3,', VTOP =',    &
   &        F10.3)
905 FORMAT ('0 ********** FIRST VALUE USED ON IFILE; CHECK IFILE ')
910 FORMAT (10X,'VMIN =',F10.3,',  VMAX =',F10.3,',  DVI=',F7.5,',    &
   &        NLOW=',I4,',  NNI=',I4)
915 FORMAT ('0 ********** END OF FILE ENCOUNTERED; CHECK IFILE ')
920 FORMAT(' ERROR IN INPUT',/,'  JTREM =',I2,'  ISCAN=',I5)
!
end subroutine RDSCAN
!
!     --------------------------------------------------------------
!
SUBROUTINE CONVSC (S,HWHMV1,R1,XF)
!
   !
   IMPLICIT REAL*8          (V)
!
!     SUBROUTINE CONVSC PERFORMS THE CONVOLUTION WITH THE SELECTED
!     SCANNING FUNCTION
!
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /RSCAN/ V1I,V2I,DVI,NNI
   COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,             &
   &                HWF3,DXF3,NX3,N3MAX
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   DIMENSION S(*),R1(*),XF(*)
!
   DATA I_1/1/
!
   CALL CPUTIM (TIME0)
   IF (ILO.GT.IHI) GO TO 60
   RATIO = DVI/DVO
   DVODX = DVO/DXF
   HWBND = HWF/DVO
   ZINT = ((V1I-VFT)/DVO)
   HWHM = HWHMV1
   ITST = -1
!
   DO 50 I = ILO, IHI
      IF (S(I).EQ.0.) GO TO 50
      IF (I.LT.ITST) GO TO 20
      ITST = 9999
      IF (JVAR.EQ.0) GO TO 10
      VI = REAL(I-1)*DVI+V1I
      HWHM = HWHMV1*(VI/V1)**2
      ITST = I+ INT(1./DVI)
10    CONTINUE
      ZSLOPE = DVODX/HWHM
      ZBOUND = HWBND*HWHM
      XNORM = DVI/HWHM
!
!PRT     WRITE(IPR,900) VI,HWHM
!
20    CONTINUE
      ZPEAK = REAL(I-1)*RATIO+ZINT
      JMAX = ZPEAK+ZBOUND+1.5
      IF (JMAX.LE.MAXF) GO TO 30
      ILAST = I-1
      GO TO 60
!
30    JMIN = ZPEAK-ZBOUND+1.5
      JMIN = MAX(JMIN,I_1)
      SUMIN = SUMIN+S(I)
      SI = XNORM*S(I)
      ZF = ( REAL(JMIN-1)-ZPEAK)*ZSLOPE
      DO 40 JF = JMIN, JMAX
         IT = ABS(ZF)+1.5
         R1(JF) = R1(JF)+SI*XF(IT)
         ZF = ZF+ZSLOPE
40    CONTINUE
!
50 END DO
   ILAST = IHI
   IPANEL = IDATA
   GO TO 70
!
60 IPANEL = 1
70 CALL CPUTIM (TIME)
   TIMCNV = TIMCNV+TIME-TIME0
   ILO = ILAST+1
!
   RETURN
!
900 FORMAT ('0 AVE PANEL WAVENUMBER = ',F12.4,5X,'HWHM = ',F10.5)
!
end subroutine CONVSC
!
!     --------------------------------------------------------------
!
SUBROUTINE PANLSC (R1,JFILE,SUMR,NPTS)
!
   !
   IMPLICIT REAL*8          (V)
!
!     SUBROUTINE PANLSC OUTPUTS THE RESULTS OF THE SCANNING FUNCTION
!     TO FILE JFILE
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /SPANEL/ V1P,V2P,DV,NLIM
   DIMENSION PNLHDR(2)
   DIMENSION R1(*),SUMR(*)
!
   EQUIVALENCE (PNLHDR(1),V1P)
!
   CALL CPUTIM (TIME0)
!
   SUMOUT = SUMR(1)
   SMIN = SUMR(2)
   SMAX = SUMR(3)
   DV = DVO
   ISTOP = 0
   NNHI = (V2-VFT)/DV+1.5
   IF (NHI.GE.NNHI) ISTOP = 1
   IF (ISTOP.EQ.1) NHI = NNHI
   NLIM = NHI-NLO+1
   V1P = VFT+ REAL(NLO-1)*DV
   V2P = VFT+ REAL(NHI-1)*DV
!
!     V1P IS FIRST FREQ OF PANEL
!     V2P IS LAST  FREQ OF PANEL
!
   CALL BUFOUT (JFILE,PNLHDR(1),NPHDRF)
   CALL BUFOUT (JFILE,R1(NLO),NLIM)
   VFT = VFT+ REAL(NLIMF-1)*DV
   IF (NPTS.GT.0) THEN
      WRITE (IPR,900) V1P,V2P,DVO,NLIM
      WRITE (IPR,905)
      NNPTS = NPTS
      IF (NPTS.GT.(NLIM/2)+1) NNPTS = NLIM/2+1
      IJLIM = NLIM-NNPTS+1
      DO 10 IJ = 1, NNPTS
         IK = IJ+IJLIM-1
         VI = V1P+ REAL(IJ-1)*DVO
         VK = V1P+ REAL(IK-1)*DVO
         JJ = NLO+IJ-1
         KK = NLO+IK-1
         WRITE (IPR,910) IJ,VI,R1(JJ),IK,VK,R1(KK)
10    CONTINUE
   ENDIF
   NLIMHI = NLIM+NLO-1
   DO 20 I = NLO, NLIMHI
      SMIN = MIN(SMIN,R1(I))
      SMAX = MAX(SMAX,R1(I))
      SUMOUT = SUMOUT+R1(I)
20 END DO
   IF (ISTOP.EQ.1) GO TO 50
   DO 30 J = NLIMF, MAXF
      R1(J-NLIMF+1) = R1(J)
30 END DO
   DO 40 J = MAXF-NLIMF+2, MAXF
      R1(J) = 0.
40 END DO
   NLO = NSHIFT+1
50 SUMR(1) = SUMOUT
   SUMR(2) = SMIN
   SUMR(3) = SMAX
   SUMR(4) = DVO
   CALL CPUTIM (TIME)
   TIMPNL = TIMPNL+TIME-TIME0
!
   RETURN
!
900 FORMAT ('0 V1P =',F12.5,' V2P =',F12.5,' DVOUT =',F12.8,' NLIM =' &
   &   ,I10)
905 FORMAT ('0')
910 FORMAT (I5,0PF12.5,1PE12.5,I15,0PF12.5,1PE12.5)
!
end subroutine PANLSC
!
!     --------------------------------------------------------------
!
SUBROUTINE PNLRCT (R1,JFILE,SUMR,NPTS)
   IMPLICIT REAL*8          (V)
!
!     SUBROUTINE PNLRCT OUTPUTS THE RESULTS OF THE SCANNING FUNCTION
!     TO FILE JFILE
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LBL4FL,LNGTH4
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /SPANEL/ V1P,V2P,DV,NLIM
   DIMENSION PNLHDR(2)
   DIMENSION R1(*),SUMR(*)
   EQUIVALENCE (PNLHDR(1),V1P)
!
   CALL CPUTIM (TIME0)
!
   SUMOUT = SUMR(1)
   SMIN = SUMR(2)
   SMAX = SUMR(3)
   DV = DVO
   ISTOP = 0
   NNHI = (V2-VFT)/DV+1.5
   IF (NHI.GE.NNHI) ISTOP = 1
   IF (ISTOP.EQ.1) NHI = NNHI
   NLIM = NHI-NLO+1
!
   V1P = VFT+ REAL(NLO-1)*DV
   V2P = VFT+ REAL(NHI-1)*DV
!
!     V1P IS FIRST FREQ OF PANEL
!     V2P IS LAST  FREQ OF PANEL
!
   CALL BUFOUT (JFILE,PNLHDR(1),NPHDRF)
   CALL BUFOUT (JFILE,R1(NLO),NLIM)
   VFT = VFT+ REAL(NLIMF-1)*DV
   IF (NPTS.GT.0) THEN
      WRITE (IPR,900) V1P,V2P,DVO,NLIM
      WRITE (IPR,905)
      NNPTS = NPTS
      IF (NPTS.GT.(NLIM/2)+1) NNPTS = NLIM/2+1
      IJLIM = NLIM-NNPTS+1
      DO 10 IJ = 1, NNPTS
         IK = IJ+IJLIM-1
         VI = V1P+ REAL(IJ-1)*DVO
         VK = V1P+ REAL(IK-1)*DVO
         JJ = NLO+IJ-1
         KK = NLO+IK-1
         WRITE (IPR,910) IJ,VI,R1(JJ),IK,VK,R1(KK)
10    CONTINUE
   ENDIF
   NLIMHI = NLIM+NLO-1
   DO 20 I = NLO, NLIMHI
      SMIN = MIN(SMIN,R1(I))
      SMAX = MAX(SMAX,R1(I))
      SUMOUT = SUMOUT+R1(I)
20 END DO
!
   IF (ISTOP.EQ.1) GO TO 50
   DO 30 J = NLIMF, MAXF
      R1(J-NLIMF+1) = R1(J)
30 END DO
   DO 40 J = MAXF-NLIMF+2, MAXF
      R1(J) = 0.
40 END DO
   NLO = NSHIFT+1
50 SUMR(1) = SUMOUT
   IPANEL = -1
   SUMR(2) = SMIN
   SUMR(3) = SMAX
   SUMR(4) = DVO
   CALL CPUTIM (TIME)
   TIMPNL = TIMPNL+TIME-TIME0
   RETURN
!
900 FORMAT('0 V1P =',F12.5,' V2P =',F12.5,' DVOUT =',F12.8,        &
   &   ' NLIM =',I10)
905 FORMAT('0')
910 FORMAT(I5,0PF12.5,1PE12.5,I15,0PF12.5,1PE12.5)
!
end subroutine PNLRCT
!
!     --------------------------------------------------------------
!
SUBROUTINE CNVRCT (S,HWHM,R1,XF)
   IMPLICIT REAL*8          (V)
!
!     SUBROUTINE CNVRCT PERFORMS THE CONVOLUTION WITH AN ALTERNATE
!     RECTANGULAR SCANNING FUNCTION (ADJACENT BOXES OF ONE SIZE,
!     EQUAL TO 2*HWHM)
!
!     THE CONVOLUTION IS A WEIGHTED SUM THAT PROPERLY WEIGHS THE INPUT
!     POINTS WITH THE FRACTION OF THAT POINT THAT COMPLETELY FALLS WITHI
!     THE OUTPUT BOX.  OUTPUT RADIANCE IS THE SUMMED RADIANCE DIVIDED
!     BY THE SUM OF THE WEIGHTS.
!
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LBL4FL,LNGTH4
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /RSCAN/ V1I,V2I,DVI,NNI
   COMMON /RCTSV/ JN,SUMJ,JFLG,RNJ,NB,IPC,VLFT,VCNT,VRGT,WGTL,WGTR
   DIMENSION S(*),R1(*),XF(*)
!
!     LBLRTM flags
!     JFLG = -1:  first time through; increment NB when box is full
!     JFLG =  0:  subsequent calls:increment NB when box full
!     JFLG =  1:  out of data, return for more; do not increment NB
!     IDATA =-1:  first time through; or, need more data
!     IDATA = 0:  data present
!     IDATA = 1:  no data present
!     IPANEL=-1:  first time through; or, after panel written
!     IPANEL= 0:  panel not full
!     IPANEL= 1:  panel is full
!
   CALL CPUTIM (TIME0)
   RATIO = DVO/DVI
!
!     During first call or if entering after writing a panel,
!     initialize: SUMJ (radiance sum), and
!                 RNJ (accumulator for number of input points in
!                      current box, i.e. the sum of the weights)
!                 JN (box counter from 1 at VFT)
!     During first call only,
!     initialize: NB (box counter from 1 at V1),
!                 IPC (output panel counter).
!
   IF (IPANEL.EQ.-1) THEN
      SUMJ = 0.
      RNJ = 0.
      JN = NLO
      IF (JFLG.EQ.-1) THEN
         NB = 1
         IPC = 1
      ENDIF
   ENDIF
!
!     Check that number of points in current panel, NNI, is correct.
!
   NNIV2 = (V2I-V1I)/DVI+1.0001
   IF (NNI.GT.NNIV2) NNI = NNIV2
!
!     Top of loop over NB boxes
!
10 IF (NLO.LE.NHI) THEN

      VCNT = V1+(NB-1)*DVO
      IF (VCNT.GT.V2) THEN
         IPANEL = 1
         RETURN
      ENDIF
!
      VLFT = VCNT-HWHM
      VRGT = VCNT+HWHM
!
!        Find lbl panel indices for points which fall within current
!        box.
!
      RL = (VLFT-V1I)/DVI+1
      RR = (VRGT-V1I)/DVI+1
!
      IL = INT(RL+0.5)
      IH = INT(RR+0.5)
!
!        Calculate weight for each end point, inner points weighted
!        as 1.  NEP is the number of endpoints in use.
!
      VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
      WGTL = (VLBLR-VLFT)/DVI
      VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
      WGTR = (VRGT-VLBLL)/DVI
      NEP = 2
!
!        Set flag if last data point on current input panel reached
!
      IF (IH.GT.NNI) THEN
         IH = NNI
         JFLG = 1
!
!        If retrieving next panel while box sum is in progress, then
!        use weight of 1. for temp. right endpoint at IH = NNI = 2400
!        calculate partial sum below, then return.  If only one point
!        is included in this sum (IL = IH = NNI), use weight of 0
!        for right point, and add only left endpoint to sum.
         VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
         WGTR = 1.
         IF (IL.EQ.IH) THEN
            WGTR = 0.
            NEP = 1
         ENDIF
      ENDIF
!
!        If returning with new panel to partially summed box, then set
!        weight for temporary left endpoint to 1.  If it's the last
!        point going into the box, then count it as final right endpoint
!        and use weight of 0 for left point (since IL = IH = 1).
!
      IF (IL.LE.1) THEN
         IL = 1
         VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
         WGTL = 1.
         IF (IL.GE.IH) THEN
            WGTL = 0.
            NEP = 1
         ENDIF
      ENDIF
!
!        If retrieving next panel while box sum is not progress, then
!        check that left edge of current output box is beyond last data
!        panel (IL.GT.NNI), if so, go back for new panel without summing
!
      IF (JFLG.EQ.1.AND.IDATA.EQ.0.AND.IL.GT.NNI) THEN
         IPANEL = 0
         JFLG = 0
         RETURN
      ENDIF
!
!        If last point on current input panel is reached, and there is
!        no more data to retrieve, then return
!
      IF (JFLG.EQ.1.AND.IDATA.EQ.1) RETURN
!
!        Compute sum for current box number NB, for all points but
!        the end points
!
      DO 20 I = IL+1, IH-1
         SUMJ = SUMJ+S(I)
20    CONTINUE
!
!        Add weighted end points to sum
      SUMJ = SUMJ+S(IL)*WGTL
      SUMJ = SUMJ+S(IH)*WGTR
!
!        Define sum of the weights, where all inner points are weighted
!        as 1, and the end points are weighted with the fraction that
!        occurs within the box.
!
      RNJ = RNJ+(IH-IL+1-NEP)+WGTL+WGTR
!
!        If out of data on current input panel, go back for more;
!        partial SUMJ, current NB, and JFLG are saved in COMMON RCTSV
!
      IF (JFLG.EQ.1.AND.IDATA.EQ.0) THEN
         IPANEL = 0
         JFLG = 0
         RETURN
      ENDIF
!
!        IPANEL=IDATA
!
      SUMIN = SUMIN+SUMJ
!
!        Compute average radiance for completed box
!
      R1(JN) = SUMJ/RNJ
!
      ILPR = IH+1
!
!        Increment current box counters, initialize SUMJ and RNJ
      JN = JN+1
      SUMJ = 0.
      RNJ = 0.
!
!        Output panel when number of boxes, NB, reaches a multiple of
!        2400, using then incrementing current output panel number, IPC.
!
      IF (NB.EQ.IPC*(NHI-NLO+1)) THEN
         IPANEL = 1
         IPC = IPC+1
         NB = NB+1
         RETURN
      ENDIF
!
!        Increment NB
      NB = NB+1

!        Go back to top of loop over NB boxes
!
      GO TO 10
   ENDIF
!
   RETURN
end subroutine CNVRCT
!
!     --------------------------------------------------------------
!
SUBROUTINE CNVVRC (S,AFOV,R1,XF)
   IMPLICIT REAL*8          (V)
!
!     SUBROUTINE CNVVRC PERFORMS THE CONVOLUTION WITH A RECTANGULAR
!     SCANNING FUNCTION OF VARIABLE SIZE, WHERE THE BOX SIZE IS
!     WAVENUMBER DEPENDENT.  V1, V2, AND DVO ARE USED TO DEFINE THE
!     CENTER OF THE OUTPUT BOXES.  BOXES OVERLAP WHERE NECESSARY TO
!     INSURE A CONSTANT DVO.
!
!     BFOV is used to determine the resolution (box size), which is
!     spectrally variable.
!
!     AFOV is passed in from calls to CNVVRC from SCANFN and SCNMRG
!     as HWHM, since the value of HWHM on Record 8.1 on TAPE5 holds
!     the place of the half angle of the instrument field of view in
!     degrees.
!
!     BOX WIDTH EQUALS V*B**2/2, AND THE SHIFT EQUALS HALF THE BOX WIDTH
!     V*(1-B**2/4), WHERE B IS THE HALF ANGLE OF THE INSTRUMENT FIELD
!     OF VIEW IN RADIANS.
!
!     THE CONVOLUTION IS A WEIGHTED SUM THAT PROPERLY WEIGHS THE INPUT
!     POINTS WITH THE FRACTION OF THAT POINT THAT COMPLETELY FALLS WITHI
!     THE OUTPUT BOX.  OUTPUT RADIANCE IS THE SUMMED RADIANCE DIVIDED
!     BY THE SUM OF THE WEIGHTS.
!
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LBL4FL,LNGTH4
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /RSCAN/ V1I,V2I,DVI,NNI
   COMMON /RCTSV/ JN,SUMJ,JFLG,RNJ,NB,IPC,VLFT,VCNT,VRGT,WGTL,WGTR
   DIMENSION S(*),R1(*),XF(*)
!
!     LBLRTM flags
!     JFLG = -1:  first time through; increment NB when box is full
!     JFLG =  0:  subsequent calls:increment NB when box full
!     JFLG =  1:  out of data, return for more; do not increment NB
!     IDATA =-1:  first time through; or, need more data
!     IDATA = 0:  data present
!     IDATA = 1:  no data present
!     IPANEL=-1:  first time through; or, after panel written
!     IPANEL= 0:  panel not full
!     IPANEL= 1:  panel is full
!
   CALL CPUTIM (TIME0)
!
!     Convert AFOV to BFOV, the half angle of the instrument field of
!     view in radians. (For IRIS-D, AFOV equals 2.5 degrees)

   BFOV = AFOV*3.141592654/180.
!
   RATIO = DVO/DVI
!
!     During first call or if entering after writing a panel,
!     initialize: SUMJ (radiance sum), and
!                 RNJ (accumulator for number of input points in
!                      current box, i.e. the sum of the weights)
!                 JN (box counter from 1 at VFT)
!     During first call only,
!     initialize: NB (box counter from 1 at V1),
!                 IPC (output panel counter).
!
   IF (IPANEL.EQ.-1) THEN
      SUMJ = 0.
      RNJ = 0.
      JN = NLO
      IF (JFLG.EQ.-1) THEN
         NB = 1
         IPC = 1
      ENDIF
   ENDIF
!
!     Check that number of points in current panel, NNI, is correct.
!
   NNIV2 = (V2I-V1I)/DVI+1.0001
   IF (NNI.GT.NNIV2) NNI = NNIV2
!
!     Top of loop over NB boxes
!
10 IF (NLO.LE.NHI) THEN

!     For current box find wavenumber at center and left/right edges.
!     For first box, VCNT equals V1.  When current box exceeds V2,
!     then exit.

      VCNT = V1+(NB-1)*DVO
      IF (VCNT.GT.V2) THEN
         IPANEL = 1
         RETURN
      ENDIF

      VLFT = VCNT*(1-BFOV**2/4)
      VRGT = VCNT*(1+BFOV**2/4)

!     Find lbl panel indices for points which fall within current box.

      RL = (VLFT-V1I)/DVI+1
      RR = (VRGT-V1I)/DVI+1
!
      IL = INT(RL+0.5)
      IH = INT(RR+0.5)
!
!     Calculate weight for each end point, inner points weighted as 1.
!     NEP is the number of endpoints in use.
!
      VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
      WGTL = (VLBLR-VLFT)/DVI
      VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
      WGTR = (VRGT-VLBLL)/DVI
      NEP = 2
!
!        Set flag if last data point on current input panel reached
!
      IF (IH.GT.NNI) THEN
         IH = NNI
         JFLG = 1
!
!        If retrieving next panel while box sum is in progress, then
!        use weight of 1. for temp. right endpoint at IH = NNI = 2400
!        calculate partial sum below, then return.  If only one point
!        is included in this sum (IL = IH = NNI), use weight of 0
!        for right point, and add only left endpoint to sum.
         VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
         WGTR = 1.
         IF (IL.EQ.IH) THEN
            WGTR = 0.
            NEP = 1
         ENDIF
      ENDIF
!
!        If returning with new panel to partially summed box, then set
!        weight for temporary left endpoint to 1.  If it's the last
!        point going into the box, then count it as final right endpoint
!        and use weight of 0 for left point (since IL = IH = 1).
!
      IF (IL.LE.1) THEN
         IL = 1
         VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
         WGTL = 1.
         IF (IL.GE.IH) THEN
            WGTL = 0.
            NEP = 1
         ENDIF
      ENDIF
!
!        If retrieving next panel while box sum is not progress, then
!        check that left edge of current output box is beyond last data
!        panel (IL.GT.NNI), if so, go back for new panel without summing
!
      IF (JFLG.EQ.1.AND.IDATA.EQ.0.AND.IL.GT.NNI) THEN
         IPANEL = 0
         JFLG = 0
         RETURN
      ENDIF
!
!        If last point on current input panel is reached, and there is
!        no more data to retrieve, then return
!
      IF (JFLG.EQ.1.AND.IDATA.EQ.1) RETURN
!
!        Compute sum for current box number NB, for all points but
!        the end points
!
      DO 20 I = IL+1, IH-1
         SUMJ = SUMJ+S(I)
20    CONTINUE
!
!        Add weighted end points to sum
      SUMJ = SUMJ+S(IL)*WGTL
      SUMJ = SUMJ+S(IH)*WGTR
!
!        Define sum of the weights, where all inner points are weighted
!        as 1, and the end points are weighted with the fraction that
!        occurs within the box.
!
      RNJ = RNJ+(IH-IL+1-NEP)+WGTL+WGTR
!
!        If out of data on current input panel, go back for more;
!        partial SUMJ, current NB, and JFLG are saved in COMMON RCTSV
!
      IF (JFLG.EQ.1.AND.IDATA.EQ.0) THEN
         IPANEL = 0
         JFLG = 0
         RETURN
      ENDIF
!
!        IPANEL=IDATA
!
      SUMIN = SUMIN+SUMJ
!
!
!        Compute average radiance for completed box
!
      R1(JN) = SUMJ/RNJ
!
!        Increment current box counters, initialize SUMJ and RNJ
      JN = JN+1
      SUMJ = 0.
      RNJ = 0.
!
!        Output panel when number of boxes, NB, reaches a multiple of
!        2400, using then incrementing current output panel number, IPC.
!
      IF (NB.EQ.IPC*(NHI-NLO+1)) THEN
         IPANEL = 1
         IPC = IPC+1
         NB = NB+1
         RETURN
      ENDIF

!        Increment NB
      NB = NB+1
!
!        Go back to top of loop over NB boxes
!
      GO TO 10

   ENDIF
!
   RETURN
end subroutine CNVVRC
!
!     --------------------------------------------------------------
!
SUBROUTINE CNVVRL (S,AFOV,R1,XF)
   IMPLICIT REAL*8          (V)
!
!     SUBROUTINE CNVVRL PERFORMS THE CONVOLUTION WITH A RECTANGULAR
!     SCANNING FUNCTION OF VARIABLE SIZE, WHERE THE BOX SIZE IS
!     WAVENUMBER DEPENDANT.  V1, V2, AND DVO ARE USED TO DEFINE THE
!     LEFT EDGE OF THE OUTPUT BOXES.  BOXES OVERLAP WHERE NECESSARY
!     TO INSURE A CONSTANT DVO.
!
!     BFOV is used to determine the resolution (box size), which is
!     spectrally variable.
!
!     AFOV is passed in from calls to CNVVRL from SCANFN and SCNMRG
!     as HWHM, since the value of HWHM on Record 8.1 on TAPE5 holds
!     the place of the half angle of the instrument field of view in
!     degrees.
!
!     BOX WIDTH EQUALS V*B**2/2, AND THE SHIFT EQUALS HALF THE BOX WIDTH
!     V*(1-B**2/4), WHERE B IS THE HALF ANGLE OF THE INSTRUMENT FIELD
!     OF VIEW IN RADIANS.
!
!     THE CONVOLUTION IS A WEIGHTED SUM THAT PROPERLY WEIGHS THE INPUT
!     POINTS WITH THE FRACTION OF THAT POINT THAT COMPLETELY FALLS WITHI
!     THE OUTPUT BOX.  OUTPUT RADIANCE IS THE SUMMED RADIANCE DIVIDED
!     BY THE SUM OF THE WEIGHTS.
!
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LBL4FL,LNGTH4
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /RSCAN/ V1I,V2I,DVI,NNI
   COMMON /RCTSV/ JN,SUMJ,JFLG,RNJ,NB,IPC,VLFT,VCNT,VRGT,WGTL,WGTR
   DIMENSION S(*),R1(*),XF(*)
!
!     LBLRTM flags
!     JFLG = -1:  first time through; increment NB when box is full
!     JFLG =  0:  subsequent calls:increment NB when box full
!     JFLG =  1:  out of data, return for more; do not increment NB
!     IDATA =-1:  first time through; or, need more data
!     IDATA = 0:  data present
!     IDATA = 1:  no data present
!     IPANEL=-1:  first time through; or, after panel written
!     IPANEL= 0:  panel not full
!     IPANEL= 1:  panel is full
!
   CALL CPUTIM (TIME0)
!
!     Convert AFOV to BFOV, the half angle of the instrument field of
!     view in radians. (For IRIS-D, AFOV equals 2.5 degrees)

   BFOV = AFOV*3.141592654/180.

   RATIO = DVO/DVI
!
!     During first call or if entering after writing a panel,
!     initialize: SUMJ (radiance sum), and
!                 RNJ (accumulator for number of input points in
!                      current box, i.e. the sum of the weights)
!                 JN (box counter from 1 at VFT)
!     During first call only,
!     initialize: NB (box counter from 1 at V1),
!                 IPC (output panel counter).
!
   IF (IPANEL.EQ.-1) THEN
      SUMJ = 0.
      RNJ = 0.
      JN = NLO
      IF (JFLG.EQ.-1) THEN
         NB = 1
         IPC = 1
      ENDIF
   ENDIF
!
!     Check that number of points in current panel, NNI, is correct.
!
   NNIV2 = (V2I-V1I)/DVI+1.0001
   IF (NNI.GT.NNIV2) NNI = NNIV2
!
!     Top of loop over NB boxes
!
10 IF (NLO.LE.NHI) THEN
!
!     For current box find wavenumber at the left and right edges.
!     For first box, VLFT equals V1.  When current box exceeds V2,
!     then exit.
!
      VLFT = V1+(NB-1)*DVO
      IF (VLFT.GT.V2) THEN
         IPANEL = 1
         RETURN
      ENDIF
!
      VCNT = VLFT*(1+BFOV**2/4)
      VRGT = VLFT*(1+BFOV**2/2)
!
!     Find lbl panel indices for points which fall within current box.
!
      RL = (VLFT-V1I)/DVI+1
      RR = (VRGT-V1I)/DVI+1
!
      IL = INT(RL+0.5)
      IH = INT(RR+0.5)
!
!     Calculate weight for each end point, inner points weighted as 1,
!     NEP is the number of endpoints in use.
!
      VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
      WGTL = (VLBLR-VLFT)/DVI
      VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
      WGTR = (VRGT-VLBLL)/DVI
      NEP = 2
!
!        Set flag if last data point on current input panel reached
!
      IF (IH.GT.NNI) THEN
         IH = NNI
         JFLG = 1
!
!        If retrieving next panel while box sum is in progress, then
!        use weight of 1. for temp. right endpoint at IH = NNI,
!        calculate partial sum below, then return.  If only one point
!        is included in this sum (IL = IH = NNI), use weight of 0
!        for right point, and add only left endpoint to sum.
         VLBLL = (V1I+(IH-1)*DVI)-DVI/2.
         WGTR = 1.
         IF (IL.EQ.IH) THEN
            WGTR = 0.
            NEP = 1
         ENDIF
      ENDIF
!
!        If returning with new panel to partially summed box, then set
!        weight for temporary left endpoint to 1.  If it's the last
!        point going into the box, then count it as final right endpoint
!        and use weight of 0 for left point (since IL = IH = 1).
!
      IF (IL.LE.1) THEN
         IL = 1
         VLBLR = (V1I+(IL-1)*DVI)+DVI/2.
         WGTL = 1.
         IF (IL.GE.IH) THEN
            WGTL = 0.
            NEP = 1
         ENDIF
      ENDIF
!
!        If retrieving next panel while box sum is not in progress, then
!        check that left edge of current output box is beyond last data
!        panel (IL.GT.NNI), if so, go back for new panel without summing
!
      IF (JFLG.EQ.1.AND.IDATA.EQ.0.AND.IL.GT.NNI) THEN
         IPANEL = 0
         JFLG = 0
         RETURN
      ENDIF
!
!        If last point on current input panel is reached, and there is
!        no more data to retrieve, then return
!
      IF (JFLG.EQ.1.AND.IDATA.EQ.1) RETURN
!
!        Compute sum for current box number NB, using a weight of 1.0
!        for all points but the end points, which use WGTL and WGTR
!
      DO 20 I = IL+1, IH-1
         SUMJ = SUMJ+S(I)
20    CONTINUE
!
!        Add weighted end points to sum
      SUMJ = SUMJ+S(IL)*WGTL
      SUMJ = SUMJ+S(IH)*WGTR
!
!        Define sum of the weights, where all inner points are weighted
!        as 1, and the end points are weighted with the fraction that
!        occurs within the box.
!
      RNJ = RNJ+(IH-IL+1-NEP)+WGTL+WGTR
!
!        If out of data on current input panel, go back for more;
!        partial SUMJ, current NB, and JFLG are saved in COMMON RCTSV
!
      IF (JFLG.EQ.1.AND.IDATA.EQ.0) THEN
         IPANEL = 0
         JFLG = 0
         RETURN
      ENDIF
!
!        IPANEL=IDATA
!
      SUMIN = SUMIN+SUMJ
!
!        Compute average radiance for completed box
!
      R1(JN) = SUMJ/RNJ
!
!        Increment current box counters, initialize SUMJ and RNJ
      JN = JN+1
      SUMJ = 0.
      RNJ = 0.
!
!        Output panel when number of boxes, NB, reaches a multiple of
!        2400, using then incrementing current output panel number, IPC.
!
      IF (NB.EQ.IPC*(NHI-NLO+1)) THEN
         IPANEL = 1
         IPC = IPC+1
         NB = NB+1
         RETURN
      ENDIF
!
!        Increment NB
      NB = NB+1
!
!        Go back to top of loop over NB boxes
!
      GO TO 10
!
   ENDIF
!
   RETURN
end subroutine CNVVRL
!
!     --------------------------------------------------------------
!
SUBROUTINE SINCSQ (XF,XSCALE)
!
!     SUBROUTINE SINCSQ SETS UP THE SINCSQ SCANNING FUNCTION
!
   COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,             &
   &                HWF3,DXF3,NX3,N3MAX
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   DIMENSION XF(*)
!
!     DATA XSCALE / 1.391557377 /
!
   XSINC2(X) = (SIN(X)/X)**2
   PI = 2.*ASIN(1.)
!
!     PI CORRESPONDS TO X=2.257609141
!
   XNORM = XSCALE/PI
   DO 10 I = 1, NFMAX
      XF(I) = 0.
10 END DO
   XF(1) = XNORM
   SUM = XF(1)
   DO 20 I = 2, NF
      X = REAL(I-1)*DXF
      XF(I) = XNORM*XSINC2(X*XSCALE)
      SUM = SUM+2.*XF(I)
20 END DO
   SUM = SUM*DXF
!
!PRT  WRITE(IPR,900) NF,DXF,SUM
!
   RETURN
!
900 FORMAT ('0',5X,'NF =',I5,',  DXF =',F7.5,',    SUM =',F18.15)
!
end subroutine SINCSQ
!
!     --------------------------------------------------------------
!
SUBROUTINE SINC (XF,XSCALE)
!
!     SUBROUTINE SINC SETS UP THE SINC SCANNING FUNCTION
!
   COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,             &
   &                HWF3,DXF3,NX3,N3MAX
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   DIMENSION XF(*)
!
!     DATA XSCALE / 1.89549425  /
!
   XSINC(X) = (SIN(X)/X)
   PI = 2.*ASIN(1.)
!
!     PI CORRESPONDS TO X=1.657400255
!
   XNORM = XSCALE/PI
   DO 10 I = 1, NFMAX
      XF(I) = 0.
10 END DO
   XF(1) = XNORM
   SUM = XF(1)
   DO 20 I = 2, NF
      X = REAL(I-1)*DXF
      XF(I) = XNORM*XSINC(X*XSCALE)
      SUM = SUM+2.*XF(I)
20 END DO
   SUM = SUM*DXF
!
!PRT  WRITE(IPR,900) NF,DXF,SUM
!
   RETURN
!
900 FORMAT ('0',5X,'NF =',I5,',  DXF =',F7.5,',    SUM =',F18.15)
!
end subroutine SINC
!
!     --------------------------------------------------------------
!
SUBROUTINE INTRPL (IFILE,JFILE)
!
   !
   IMPLICIT REAL*8          (V)
!
!**********************************************************************
!
!     INTERPOLATION FUNCTION DRIVER: FOUR POINT VERSION
!
!           A.E.R. INC.           (MARCH 1990)
!
!**********************************************************************
!
!     THE INPUT DATA WILL BE PUT INTO T(5) WITH THE LAST
!     4 POINTS OF THE PREVIOUS PANEL PUT INTO T(1 TO 4).
!     THIS SCHEME PERMITS 6 POINT INTERPOLATION.
!
!     S IS NOMINALLY 2401 POINTS BUT MAY NEED TO BE EXTENDED BY TWO (2)
!     POINTS TO PERMIT 4 POINT INTERPOLATION UP TO THE LAST DATA POINT.
!
   COMMON T(2410),R(2401)
   DIMENSION S(2406)
   EQUIVALENCE (S(1),T(5))
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WN2   ,DV ,V1C,V2C,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /INPNL/ V1I,V2I,DVI,NNI
   COMMON /OUTPNL/ V1J,V2J,DVJ,NNJ
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /FLFORM/ CFORM
!
   DIMENSION FILHDR(2),RSTAT(3)
!
   EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(6),ISCHDR) , (FSCDID(12),XSCID),              &
   &            (FSCDID(13),XHWHM) , (FSCDID(14),IDABS),              &
   &            (FSCDID(16),LAYR1)
!
   CHARACTER*12 BCD,HTRANS,HABSRB,HRADIA
   CHARACTER CFORM*11,SCNOUT*6,CTAPE*4
   LOGICAL OP
!
   DATA I_1/1/, I_1000/1000/
!
   DATA HTRANS / 'TRANSMISSION'/,HABSRB / ' ABSORPTION '/,           &
   &     HRADIA / ' RADIANCE   '/
   DATA SCNOUT / '      '/,CTAPE / 'TAPE'/
!
!----------------------------------------------------------------------
!     JEMIT=-1  INTERPOLATE ABSORPTION
!     JEMIT=0   INTERPOLATE TRANSMISSION
!     JEMIT=1   INTERPOLATE EMISSION
!     JEMIT=2   INTERPOLATE OPTICAL DEPTH
!----------------------------------------------------------------------
!
10 CONTINUE
   CALL CPUTIM (TIME1)
   TIMRDF = 0.0
   TIMCNV = 0.0
   TIMPNL = 0.0
!
   READ (IRD,900,END=50) DVO,V1,V2,JEMIT,I4PT,IUNIT,IFILST,NIFILS,   &
   &                      JUNIT,NPTS
   IF (DVO.LE.0.) GO TO 40
!
!     V2 IS ONLY APPROXIMATE
!
   NUM = (((V2-V1)/DVO)+0.5)
   V2 = V1+ REAL(NUM)*DVO
   NUM = NUM+1
   WRITE (IPR,905) V1,V2,DVO,NUM,JEMIT,I4PT,IUNIT,IFILST,JUNIT,NPTS
!
!     SET INPUT(IFILE), OUTPUT(JFILE) UNITS.
!
   IF (IUNIT.LE.0) IUNIT = IFILE
   IFILE = IUNIT
   INQUIRE (UNIT=IFILE,OPENED=OP)
   IF (.NOT.OP) THEN
      WRITE (SCNOUT,910) CTAPE,IFILE
      OPEN (IFILE,FILE=SCNOUT,STATUS='UNKNOWN',FORM=CFORM)
   ENDIF
   IFILST = MAX(IFILST,I_1)
   IF (NIFILS.LE.0) NIFILS = 99
   IF (JUNIT.LE.0) JUNIT = JFILE
   JFILE = JUNIT
   INQUIRE (UNIT=JFILE,OPENED=OP)
   IF (.NOT.OP) THEN
      WRITE (SCNOUT,910) CTAPE,JFILE
      OPEN (JFILE,FILE=SCNOUT,STATUS='UNKNOWN',FORM=CFORM)
      REWIND JFILE
   ENDIF
!
   REWIND IFILE
   IF (IFILST.GT.1) CALL SKIPFL (IFILST-1,IFILE,IEOF)
!
!     BUFFER IN THE FILE HEADER ON UNIT (IFILE)
!     BUFFER OUT ON UNIT (JFILE)
!
   CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)
   IF (IEOF.EQ.0) GO TO 10
!
   WRITE (IPR,915) XID,(YID(M),M=1,2)
   WRITE (IPR,920) LAYR1,LAYER
   WRITE (IPR,925) SECANT,PAVE,TAVE,DV,V1C,V2C
   WRITE (IPR,930) WN2,(HMOLID(M),WK(M),M=1,NMOL)
!
   JABS = 0
   IDABS = 0
   IF (JEMIT.LT.0) THEN
      JABS = 1
      JEMIT = 0
      IDABS = -1
   ENDIF
!
   ISCAN = ISCHDR
   IF (ISCAN.LE.0.OR.XSCID.EQ.-99.) ISCAN = 0
   IF (ISCHDR.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCHDR
   ISCHDR = ISCAN+10
   V1C = V1
   V2C = V2
   DV = DVO
!
   SCNID = 100*JEMIT
   XSCID = SCNID+0.01
!
   CALL BUFOUT (JFILE,FILHDR(1),NFHDRF)
!
   ICNVRT = 0
   JTREM = -1
   IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.2)) JTREM = 2
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1
   ISCANT = MOD(ISCAN,I_1000)
   IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 2
   WRITE (IPR,935) IEMIT,JEMIT,JTREM
   WRITE (IPR,940) IFILE,IFILST,NIFILS,JEMIT,JABS
   IF (JTREM.LT.0) then
      WRITE(IPR,*) ' Invalid JTREM in INTRPL '
      STOP ' Invalid JTREM in INTRPL '
   endif
!
   IDATA = -1
!
!     NEED TO SAVE LAST IBOUND POINTS OF EACH PANEL TO ATTACH TO NEXT
!
   IBOUND = 4
!
!     VBOT IS LOWEST NEEDED WAVENUMBER, VTOP IS HIGHEST
!
   BOUND = REAL(IBOUND)*DV
   VBOT = V1-BOUND
   VTOP = V2+BOUND
!
   IF (JEMIT.EQ.0.AND.IDABS.EQ.0) BCD = HTRANS
   IF (JEMIT.EQ.0.AND.IDABS.EQ.-1) BCD = HABSRB
   IF (JEMIT.EQ.1) BCD = HRADIA
   IF (NPTS.GT.0) WRITE (IPR,945) BCD
!
!     ZERO OUT T(1 TO IBOUND)
!
   DO 20 II = 1, IBOUND
      T(II) = 0.0
20 END DO
!
!     READ FROM IFILE UNTIL THE FIRST REQUIRED POINT IS REACHED
!     AND LOAD DATA INTO S
!
   CALL RDPANL (S,JTREM,IFILE,ISCAN,JEMIT,ICNvRT)
   IF (IEOFSC.LE.0) GO TO 30
!
!     DO INTERPOLATION
!
   CALL INTERP (IFILE,JFILE,I4PT,IBOUND,NPTS,JTREM,ISCAN,JEMIT,   &
      RSTAT,ICNVRT)
!
   CALL CPUTIM (TIME2)
   CALL ENDFIL (JFILE)
!
!     WRITE STATISTICS
!
   WRITE (IPR,950) RSTAT(1),RSTAT(2),RSTAT(3)
   TIMTOT = TIME2-TIME1
   TIMCNV = TIMTOT-TIMRDF-TIMPNL
   WRITE (IPR,955) TIMTOT,TIMRDF,TIMCNV,TIMPNL
   GO TO 10
!
30 CONTINUE
   WRITE (IPR,960) IFILE
!
   GO TO 10
!
40 CONTINUE
   WRITE (IPR,965)
!
   RETURN
!
50 CONTINUE
   WRITE (IPR,970) IRD
   STOP ' INTRPL'
!
900 FORMAT (3F10.3,2I5,15X,5I5)
905 FORMAT (5X,'V1=',F14.8,' V2=',F14.8,' DVO=',E14.6,' NUM=',I8,/5X, &
   &   'JEMIT=',I3,' I4PT=',I3,' IUNIT=',I3,' IFILST=',I3,/5X,        &
   &   'JUNIT=',I3,' NPTS=',I5)
910 FORMAT (A4,I2.2)
915 FORMAT (//,' ***INTRPL***',/,'0',10A8,2X,2(1X,A8,1X))
920 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)
925 FORMAT ('  SECANT =',F15.5,/'  PRESS(MB) =',F12.5/'  TEMP(K) =',  &
   &   F11.2,/'  DV(CM-1) = ',F12.8,/,'  V1(CM-1) = ',F12.6,/,        &
   &   '  V2(CM-1) = ',F12.6)
930 FORMAT (/,'  COLUMN DENSITY (MOLECULES/CM**2)',/,5X,'WBROAD = ',  &
   &   1PE10.3,/(5X,A6,' = ',1PE10.3))
935 FORMAT (5X,'IEMIT=',I5,' JEMIT=',I5,' JTREM=',I5)
940 FORMAT (5X,'INPUT FILE NUMBER =',I3,' ,IFILST = ',I3,             &
   &   ' ,NIFILS = ',I3,',JEMIT =',I2,' ,JABS =',I2)
945 FORMAT (///,'0',5X,A12,/)
950 FORMAT (/,5X,'SUMOUT =',1P,E16.9,'  MIN =',E16.9,'  MAX =',E16.9)
955 FORMAT (/,5X,'TIME: TOTAL = ',F8.3,' READ = ',F8.3,' INTERP = ',  &
   &   F8.3,' WRITE = ',F8.3)
960 FORMAT (/,5X,'INTRPL- ERROR: EOF ON INPUT UNIT ',I4,              &
   &   ' BEFORE V1 WAS REACHED')
965 FORMAT (/,5X,'END OF INTERPOLATION REQUESTS')
970 FORMAT (/,5X,' INTRP - ERROR: EOF ON STANDARD INPUT, UNIT = ',I4)
!
end subroutine INTRPL
!
!     --------------------------------------------------------------
!
SUBROUTINE RDPANL (S,JTREM,IFILE,ISCAN,JEMIT,ICNVRT)
!
   !
   IMPLICIT REAL*8          (V)
!
!     SUBROUTINE RDPANL INPUTS PANELS FROM IFILE RESULTING FROM THE
!     LBLRTM CALCULATION
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /INPNL/ VMIN,VMAX,DVI,NNI
   COMMON /RPANL/ V1P,V2P,DVP,NLIMP
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   DIMENSION DUMMY(2),PNLHDR(2)
   DIMENSION S(*)
!
   EQUIVALENCE (PNLHDR(1),V1P)
!
   DATA I_1000/1000/
!
!----------------------------------------------------------------------
!
   CALL CPUTIM (TIME1)
   IDUM1 = 0
   IDUM2 = 0
   ISCANT = MOD(ISCAN,I_1000)
   IF (JTREM.EQ.0.AND.ISCANT.GE.1) GO TO 70
   IF (ISCAN.LT.1) THEN
      IF (JTREM.EQ.1) IDUM1 = 1
      IF (JTREM.EQ.2) IDUM2 = 1
   ENDIF
10 CALL BUFIN (IFILE,IEOFSC,PNLHDR(1),NPHDRF)
   IF (IEOFSC.LE.0) THEN
      WRITE (IPR,900)
      GO TO 60
   ELSE
      VMIN = V1P
      VMAX = V2P
      DVI = DVP
      NNI = NLIMP
   ENDIF
!
   IF ((IDATA.EQ.-1).AND.(VMIN.GT.VBOT)) WRITE (IPR,905)
   IDATA = 0
   IF (VMAX.GE.VBOT) GO TO 20
   IF (IDUM2.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)
   CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)
   IF (IDUM1.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)
   GO TO 10
!
20 IF (JTREM.EQ.0) THEN
      CALL BUFIN (IFILE,IEOFSC,S(1),NNI)
      IF (JEMIT.NE.2.AND.ICNVRT.EQ.0) THEN
         DO 30 I = 1, NNI
            SI = S(I)
            S(I) = 1.
            IF (SI.GT.0.) THEN
               IF (SI.GE.ARGMIN) THEN
                  S(I) = EXPMIN
               ELSE
                  S(I) = EXP(-SI)
               ENDIF
            ENDIF
30       CONTINUE
      ENDIF
   ELSE
!
      IF (IDUM2.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)
      CALL BUFIN (IFILE,IEOFSC,S(1),NNI)
      IF (IDUM1.EQ.1) CALL BUFIN (IFILE,IEOFSC,DUMMY(1),2)
   ENDIF
!
   IF (JABS.EQ.1.AND.ICNVRT.EQ.0) THEN
      DO 40 I = 1, NNI
         S(I) = 1.-S(I)
40    CONTINUE
   ENDIF
   IF (JEMIT.EQ.2.AND.ICNVRT.EQ.0) THEN
      DO 50 I = 1, NNI
         S(I) = - LOG(S(I))
50    CONTINUE
   ENDIF
!
   VMIN = VMIN-4.0*DVI
   NNI = NNI+4
   ILO = 1
   IHI = NNI
   DIF = (VMIN-VBOT)/DVI
   IF (DIF.LT.0.) ILO = -DIF+1.5
   IF (VMAX.GT.VTOP) THEN
      IHI = (VTOP-VMIN)/DVI+1.5
      IDATA = 1
   ENDIF
!
60 CALL CPUTIM (TIME2)
   TIMRDF = TIMRDF+TIME2-TIME1
!
   RETURN
!
70 WRITE (IPR,910) JTREM,ISCAN
!
   RETURN
!
900 FORMAT ('0 ********** END OF FILE ENCOUNTERED; CHECK IFILE ')
905 FORMAT ('0 ********** FIRST VALUE USED ON IFILE; CHECK IFILE ')
910 FORMAT (' ERROR IN INPUT',/,'  JTREM =',I2,'  ISCAN=',I5)
!
end subroutine RDPANL
!
!     --------------------------------------------------------------
!
SUBROUTINE OTPANL (R1,JFILE,NPTS)
!
   !
   IMPLICIT REAL*8          (V)
!
!     SUBROUTINE OTPANL OUTPUTS THE RESULTS OF THE INTERPOLATION ON
!     TO FILE JFILE
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /OUTPNL/ V1P,V2P,DVP,NLIM
   DIMENSION PNLHDR(2),R1(*)
!
   EQUIVALENCE (PNLHDR(1),V1P)
!
!----------------------------------------------------------------------
!
   CALL CPUTIM (TIME1)
   IF (NLIM.LE.0) GO TO 20
!
   CALL BUFOUT (JFILE,PNLHDR(1),NPHDRF)
   CALL BUFOUT (JFILE,R1(1),NLIM)
!
   IF (NPTS.GT.0) THEN
      WRITE (IPR,900) V1P,V2P,DVP,NLIM
      WRITE (IPR,905)
      NNPTS = NPTS
      IF (NPTS.GT.(NLIM/2)+1) NNPTS = NLIM/2+1
      IJLIM = NLIM-NNPTS+1
      DO 10 IJ = 1, NNPTS
         IK = IJ+IJLIM-1
         VIJ = V1P+ REAL(IJ-1)*DVP
         VIK = V1P+ REAL(IK-1)*DVP
         WRITE (IPR,910) IJ,VIJ,R1(IJ),IK,VIK,R1(IK)
10    CONTINUE
   ENDIF
!
20 CALL CPUTIM (TIME2)
   TIMPNL = TIMPNL+TIME2-TIME1
!
   RETURN
!
900 FORMAT ('0 V1P =',F12.5,' V2P =',F12.5,' DVP =',F12.8,' NLIM =',  &
   &        I8)
905 FORMAT ('0')
910 FORMAT (I5,0PF12.5,1PE12.5,I15,0PF12.5,1PE12.5)
!
end subroutine OTPANL
!
!     --------------------------------------------------------------
!
SUBROUTINE INTERP (IFILE,JFILE,I4PT,IBOUND,NPTS,JTREM,ISCAN,      &
&                   JEMIT,RSTAT,ICNVRT)
!
   IMPLICIT REAL*8          (V)
   real*8 p
!
!**********************************************************************
!     THIS SUBROUTINE INTERPOLATES THE SPECTRAL DATA FROM IFILE, ON
!     A GRID DEFINED BY V1I,V2I,DVI, AND NNI, ONTO THE GRID FROM
!     V1 TO V2 WITH A DV OF DVO AND WRITES THE RESULT TO JFILE.
!     THE INTERPOLATION IS EITHER LINEAR (I4PT = 0) OR 4 POINT (I4PT
!     = 1).   IBOUND IS THE NUMBER OF POINTS NEEDED FROM THE PREVIOUS
!     INPUT PANEL, WHILE NPTS IS THE NUMBER OF POINTS TO BE PRINTED
!     AT THE BEGINNING AND END OF EACH OUTPUT PANEL. JTREM,ISCAN, AND
!     JEMIT RELATE TO THE INPUT DATA AND ARE NEEDED BY RDPANL.
!     RSTAT(3) RETURNS THE SUM, MIN, AND MAX OF THE INTERPOLATED
!     SPECTRUM.
!**********************************************************************
!
!     THE INPUT DATA WILL BE PUT INTO T(5) WITH THE LAST
!     IBOUND POINTS OF THE PREVIOUS PANEL PUT INTO T(1-4)
!
   COMMON T(2410),R(2401)
   DIMENSION S(2406)
   EQUIVALENCE (T(5),S(1))
!
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /INPNL/ V1I,V2I,DVI,NNI
   COMMON /OUTPNL/ V1J,V2J,DVJ,NNJ
   DIMENSION C1(0:202),C2(0:202),C3(0:202),C4(0:202),RSTAT(3)
!
   DATA NUMCOF / 201 /
   DATA I_2400/2400/
!
   CALL CPUTIM (TIME1)
!
!     SET UP FOUR POINT INTERPOLATION COEFICIENTS FOR P FOR 201
!     POINTS BETWEEN 0 AND 1.0, with an extra point at each end
!
   I_0 = 0
   IF (I4PT.NE.0) THEN
      XNUMCF = REAL(NUMCOF)
      DO 10 IP = I_0, NUMCOF+1
         P = ( REAL(IP)-1.0)/(XNUMCF-1.0)
         PP = P**2
         C1(IP) = -P/2.0*(1-P)**2
         C2(IP) = 1.0-PP*(3.0-2.0*P)+PP/2.0*(1.0-P)
         C3(IP) = PP*(3.0-2.0*P)+P/2.0*(1.0-P)**2
         C4(IP) = -PP/2*(1.0-P)
10    CONTINUE
   ENDIF
!
!     V1 IS REQUESTED LOWER V, VBOT = V1-VBOUND.  VBOUND ALLOWS FOR
!     4 POINT INTERPOLATION OF THE FIRST DATA POINT.
!     THE FIRST PANEL OF DATA IS STORED IN T, STARTING AT T(5)
!     WITH A CORRESPONDING WAVENUMBER V1I.
!     T(1-4) ARE USED TO STORE THE LAST IBOUND POINTS FROM THE
!     PREVIOUS PANEL, BUT ARE ZEROED OUT FOR THE FIRST PANEL.
!     THE INTERPOLATED POINTS ARE STORED IN THE ARRAY R.
!     THE INDEX I REFERS TO THE INPUT POINTS, J TO THE OUTPUT POINTS.
!     P VARIES FROM 0 TO 1 AND IS THE FRACTIONAL DISTANCE OF THE
!     CURRENT OUTPUT WAVENUMBER VJ TO THE NEXT LOWEST INPUT WAVENUMBER
!     RELATIVE TO THE INPUT DV: P = (VJ-VI(II))/DVI
!     INRANG IS 0 IF VJ IS WITHIN THE RANGE OF THE INPUT DATA, -1
!     IF VJ IS LESS THAN THE INPUT DATA, AND 1 IF IT IS GREATER
!
!     INITIALIZE THE VARIABLES
!
   V1J = V1
   DVJ = DVO
   VRATIO = DVJ/DVI
   VJ = V1J
!
   RMIN = 1.0E15
   RMAX = -1.0
   RSUM = 0.0
!
!     EXTRAPOLATE DOWN TO V1I-DVI SO THAT THE POINT I=4 IS AVAILABLE
!     FOR THE FIRST PANEL.  THIS ALLOWS 4 POINT INTERPOLATION BETWEEN
!     V1I AND V1I+DVI
!
   T(4) = 2.0*T(5)-T(6)
!
!     LOOP OVER THE OUTPUT PANELS
!
!     IF V1J .LT. V1I, THEN ZERO FILL UP TO V1I.
!
20 IF (V1J.LT.V1I) THEN
      J1 = 1
      J2 = MIN(INT((V1I-V1J)/DVJ)+1,I_2400)
!
!     FILL IN
!
      DO 30 J = J1, J2
         R(J) = 0.0
30    CONTINUE
!
      V2J = V1J+DVJ*(J2-1)
      NNJ = J2
      CALL OTPANL (R,JFILE,NPTS)
      V1J = V2J+DVJ
      VJ = V1J
!
      GO TO 20
   ENDIF
!
!     AT THIS POINT, VJ >= V1I
!
40 CONTINUE
!
!     I INDEXES THE LARGEST VI .LE. VJ
!     AND AT THIS POINT SHOULD .GE. 1.
!
   I = (VJ-V1I)/DVI+1.00001
   IF (I.LT.1) THEN
      WRITE (IPR,*) ' INTERP-ERROR: I SHOULD >= 1, IS ',I
   ENDIF
   VI = V1I+DVI* REAL(I-1)
!
!     P IS INCREMENTED BY ADDING DVJ/DVI BUT WILL BE REINITIALIZED
!     HERE FOR EACH OUTPUT PANEL TO AVOID THE ACCUMULATION OF
!     TRUNCATION ERRORS
!
   P = (VJ-VI)/DVI
!
   J1 = INT((VJ-V1J)/DVJ+1.001)
   J2 =                                                              &
   &   MIN(INT((V2-V1J)/DVJ+1.001),INT((V2I-DVI-V1J)/DVJ+1.),I_2400)
!
!     LOOP OVER A SINGLE OUTPUT PANEL
!
   IF (I4PT.GT.0) THEN
!
!     4 POINT INTERPOLATION
!
      DO 50 J = J1, J2
!
!     PERFORM INTERPOLATION
!
         IP = P*XNUMCF+1.00001
         R(J) = C1(IP)*T(I-1)+C2(IP)*T(I)+C3(IP)*T(I+1)+ C4(IP)*T(I+ &
            2)
!
!     ACCUMULATE STATISTICS
!
         RMIN = MIN(RMIN,R(J))
         RMAX = MAX(RMAX,R(J))
         RSUM = RSUM+R(J)
!
!     INCREMENT P AND I
!
         P = P+VRATIO
         IF (P.GE.1.0) THEN
            I = I+P
            P = P- REAL(INT(P))
         ENDIF
!
50    CONTINUE
!
   ELSE
!
!     LINEAR INTERPOLATION
!
      DO 60 J = J1, J2
!
!     PERFORM INTERPOLATION
!
         R(J) = T(I)*(1.0-P)+T(I+1)*P
!
!     ACCUMULATE STATISTICS
!
         RMIN = MIN(RMIN,R(J))
         RMAX = MAX(RMAX,R(J))
         RSUM = RSUM+R(J)
!
!     INCREMENT P AND I
!
         P = P+VRATIO
         IF (P.GE.1.0) THEN
            I = I+P
            P = P- REAL(INT(P))
         ENDIF
60    CONTINUE
!
   ENDIF
!
!     VJ IS THE FREQUENCY OF THE NEXT OUTPUT POINT (NOT THE LAST
!     POINT IS THE CURRENT PANEL)
!
   VJ = V1J+DVJ*J2
!
!     IF THE OUTPUT PANEL IS FULL OR IF V2 REACHED,
!     WRITE OUT THE PANEL
!
   IF (J2.GE.2400.OR.VJ.GE.V2) THEN
      NNJ = J2
      V2J = V1J+DVJ*(J2-1)
      CALL OTPANL (R,JFILE,NPTS)
      V1J = V2J+DVJ
      J2 = 0
   ENDIF
!
!     IF REACHED V2, THEN FINISH
!
   IF (VJ.GE.V2) GO TO 100
!
!     IF THE INPUT FILE REACHED AN EOF, THEN ZERO FILL TO END
!
   IF (IEOFSC.LE.0) GO TO 80
!
!     IF THE DATA FROM CURRENT INPUT PANEL IS EXHAUSTED, GET MORE
!
   IF (I.GE.NNI-2) THEN
!
!     SHIFT THE LAST IBOUND POINTS DOWN TO T(1-4)
!
      DO 70 II = 1, IBOUND
         T(II) = T(II+NNI-IBOUND)
70    CONTINUE
!
!     GET THE NEXT PANEL OF DATA AND RESET I
!
      CALL RDPANL (S,JTREM,IFILE,ISCAN,JEMIT,ICNVRT)
      IF (IEOFSC.LE.0) THEN
!
!     IF EOF ON INPUT FILE, THEN EXTRAPOLATE OUT TWO MORE
!     POINTS BEYOND I=NNI SO THAT 4 POINT INTERPOLATION CAN
!     BE PERFORMED UP TO VJ=V2I. (ACTUALLY, ONLY T(NNI+1) NEED
!     BE EXTRAPOLATED, T(NNI+2) NEED ONLY BE DEFINED.)
!     EXTEND THE INPUT PANEL BY ONE POINT  AND LOOP AROUND THE
!     INTERPOLATION ONE LAST TIME
!
         T(NNI+1) = 2.0*T(NNI)-T(NNI-1)
         T(NNI+2) = 0.0
         V2I = V2I+DVI
      ENDIF
   ENDIF
!
!     LOOP BACK
!
   GO TO 40
!
80 CONTINUE
   J1 = J2+1
   J2 = MIN(INT((V2-V1J)/DVJ+1.0001),I_2400)
!
   DO 90 J = J1, J2
      R(J) = 0.0
90 END DO
   VJ = V1J+DVJ*J2
!
!     IF THE OUTPUT PANEL IS FULL OR IF V2 REACHED,
!     WRITE OUT THE PANEL
!
   IF (J2.GE.2400.OR.VJ.GE.V2) THEN
      NNJ = J2
      V2J = V1J+DVJ*(J2-1)
      CALL OTPANL (R,JFILE,NPTS)
      V1J = V2J+DVJ
      J2 = 0
   ENDIF
!
!     IF REACHED V2, THEN FINISH
!
   IF (VJ.LT.V2) GO TO 80
!
100 CONTINUE
!
   RSTAT(1) = RSUM*DVJ
   RSTAT(2) = RMIN
   RSTAT(3) = RMAX
!
   CALL CPUTIM (TIME2)
   TIMCNV = TIMCNV+TIME2-TIME1
!
   RETURN
!
end subroutine INTERP
!
!     --------------------------------------------------------------
!
SUBROUTINE FLTRFN (IFILE)
!
!     NFLTPT sets the maximum number of points in the incoming filter
!
   USE lblparams, ONLY: NFLTPT
!
   IMPLICIT REAL*8 (V)
!
   COMMON S(4650),R1(5750)
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ

   LOGICAL OP
!
   COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DVC,V1C,V2C,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /RSCAN/ V1I,V2I,DVI,NNI
   COMMON /COMFLT/ V1F,V2F,DVF,NPTS,NPTF,JEMIT,IUNIT,IFILST,NIFILS,  &
   &                HEDDR(9),XF(NFLTPT),SUMFLT
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /FLFORM/ CFORM
   DIMENSION FILHDR(2)
   CHARACTER*80 CVAR
   CHARACTER*11 CFORM
   character ctape*4,fltinf*7,fltout*7
   integer*4 itest,itest2
!
   EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(6),ISCAN) , (FSCDID(7),IPLOT),                &
   &            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                &
   &            (FSCDID(12),SCNID) , (FSCDID(13),HWHM),               &
   &            (FSCDID(16),LAYR1)
!
   DATA i_2/2/, I_1000/1000/
!
   DATA FLTINF / '       '/,FLTOUT / '       ' /,                    &
   &     CTAPE / 'TAPE'/
!
   data itest / 0 /, itest2 / 0 /
   save itest
   save itest2
!
   NLIMF = 2401
   NREN = 0
   IPRT = 1
   NSHIFT = 0
10 READ (IRD,900) V1F,DVF,NPTF,JEMIT,IUNIT,IFILST,NIFILS,junit,HEDDR
!
!     Test to ensure NPTF is less than NFLTPT, the maximum number
!     of filter points allowed
!
   IF (NPTF.GT.NFLTPT) THEN
      WRITE(IPR,*) 'FLTRFN: NPTS > NFLTPT limit', NFLTPT
      STOP 'FLTRFN: NPTS > NFLTPT limit'
   ENDIF
!
   JABS = 0
   IF (JEMIT.GE.0) GO TO 20
   JEMIT = 0
   JABS = 1
20 IEOFT = 1
   IF (IUNIT.LE.0) IUNIT = IFILE
   IFILE = IUNIT
   IF (NIFILS.LE.0) NIFILS = 99
!
   IF (V1F.LT.0) RETURN

!
   WRITE (IPR,905)
!
!     DVF < 0 option flags V1F value to be the center frequency
!     Check that there NPTF is odd (to ensure a center frequency),
!     save center frequency value, and reset V1F to endpoint value.

   if (DVF.lt.0.) then

      dvf = abs(dvf)

      if (mod((nptf-1),i_2).ne.0) then
         write(*,*) 'Use of V1F as center frequency requires odd     &
         &           number of points'
         write(ipr,*) 'Use of V1F as center frequency requires odd   &
         &           number of points, stopping in FLTFRN'
         stop 'FLTRFN'
      endif

      V1F_center = V1F
      nptf_half = (abs(nptf)-1)/2
      V1F = V1F_center - DVF* REAL(nptf_half)
      write(ipr,*) ' ``````````````````````````````'
      write(ipr,*) ' Use of V1F as center frequency:'
      write(ipr,*) '   V center = ',V1F_center
      write(ipr,*) '   V1F      = ',V1F
      write(ipr,*) " ''''''''''''''''''''''''''''''"
   endif

   REWIND IFILE
   inquire(ifile,opened=op)
   IF (.NOT.OP) THEN
      WRITE (FLTINF,970) CTAPE,IFILE
      OPEN (IFILE,FILE=FLTINF,STATUS='UNKNOWN',FORM=CFORM)
      REWIND IFILE
   ENDIF

   JFILE = junit
   IF (JFILE.NE.0) THEN
      INQUIRE(JFILE,OPENED=OP)
      IF (ITEST.EQ.0) THEN
         IF (OP) CLOSE(JFILE)
         ITEST = 1
         FLTOUT = 'FLT_OUT'
         OPEN (JFILE,FILE=FLTOUT,STATUS='UNKNOWN')
         REWIND JFILE
      ENDIF
   ENDIF

   IF (IFILST.GT.1) CALL SKIPFL (IFILST-1,IFILE,IEOF)
   IEOFSC = 0
   ISTOP = 0
!
   IF (NPTF.LE.0) GO TO 30
   NPTS = NPTF
30 V2F = V1F+DVF* REAL(NPTS-1)
   WRITE (IPR,910) V1F,V2F,DVF,NPTF,JEMIT,JABS,IUNIT,IFILST,NIFILS,  &
   &                HEDDR
   V1 = V1F
   V2 = V2F
   DV = DVF
   IF (NPTF.LE.0) GO TO 40
   READ (IRD,915) CVAR
   READ (IRD,CVAR) (XF(I),I=1,NPTS)
   WRITE (IPR,CVAR) (XF(I),I=1,NPTS)
!
!     MAKE ADJUSTMENT FOR END POINT CORRECTIONS
!
   XF(1) = 0.5*XF(1)
   XF(NPTS) = 0.5*XF(NPTS)
40 SUMFLT = 0.0
   DO 50 I = 1, NPTS
      SUMFLT = SUMFLT+XF(I)
50 END DO
   SUMFLT = SUMFLT*DVF
!
60 RFILTR = 0.0
   VFT = V1
   VBOT = V1
   VTOP = V2
   TIMRDF = 0.0
   TIMCNV = 0.0
   CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)
   ISCHDR = ISCAN
   IF (ISCAN.LE.0.OR.SCNID.EQ.-99.) ISCAN = 0
   IF (ISCHDR.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCHDR
   IF (MOD(ISCAN,I_1000).EQ.0) GO TO 70
   JEMSCN = SCNID/100.
   IF (JEMIT.EQ.JEMSCN) GO TO 70
   WRITE (*,920)
   WRITE (IPR,920)
   CALL SKIPFL (1,IFILE,IEOF)
   IF (IEOF.EQ.0) GO TO 10
   GO TO 60
70 CONTINUE
   IF (IEOF.LT.1) GO TO 10
!
!     JEMIT=-1 FILTER PASSED OVER ABSORPTION
!     JEMIT=0  FILTER PASSED OVER TRANSMISSION
!     JEMIT=1  FILTER PASSED OVER EMISSION
!
   JTREM = -1
   IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1
   ISCANT = MOD(ISCAN,I_1000)
   IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 2
   IF (JTREM.LT.0) GO TO 10
   WRITE (IPR,925) XID,(YID(M),M=1,2)
   WRITE (IPR,930) LAYR1,LAYER
   WRITE (IPR,935) SECANT,PAVE,TAVE,DVC,V1C,V2C
   WRITE (IPR,940) WBROAD,(HMOLID(M),WK(M),M=1,NMOL)
   WRITE (IPR,945) V1F,V2F,DVF,NPTF,IEMIT,JEMIT,JABS,IUNIT,IFILST,   &
   &   NIFILS,HEDDR
   IF ((JFILE.NE.0).AND.(ITEST2.EQ.0)) THEN
      WRITE (JFILE,925) XID,(YID(M),M=1,2)
      WRITE (JFILE,935) SECANT,PAVE,TAVE,DVC,V1C,V2C
      WRITE (JFILE,940) WBROAD,(HMOLID(M),WK(M),M=1,NMOL)
      WRITE(JFILE,980)
      ITEST2 = 1
   ENDIF
   IDATA = -1
80 CALL CPUTIM (TIMEO)
   CALL RDSCAN (S,JTREM,IFILE,ISCAN,IPRT)
   CALL CPUTIM (TIME)
   TIMRDF = TIMRDF+TIME-TIMEO
   IF (IEOFSC.NE.1) GO TO 90
   CALL CNVFLT (S,RFILTR,XF)
   IF (IDATA.EQ.1) GO TO 90
   GO TO 80
90 CALL CPUTIM (TIME)
   WRITE (IPR,950) TIME,TIMRDF,TIMCNV
   IF (JEMIT.NE.1) THEN
      TRNSM = RFILTR
      RFILTR = RFILTR*DVC/SUMFLT
      IF (JABS.EQ.0) WRITE (IPR,955) RFILTR,SUMFLT,TRNSM
      IF (JABS.EQ.1) WRITE (IPR,960) RFILTR,SUMFLT,TRNSM
      IF (JFILE.NE.0) THEN
         WRITE(JFILE,975) RFILTR
      ENDIF
   ELSE
100   RFILTR = RFILTR*DVC
      WRITE (IPR,965) RFILTR,SUMFLT,RFILTR/SUMFLT
      IF (JFILE.NE.0) THEN
         WRITE(JFILE,975) RFILTR,SUMFLT,RFILTR/SUMFLT
      ENDIF
   ENDIF
110 IF (IEOFSC.EQ.1) CALL SKIPFL (1,IFILE,IEOFSC)
   IEOFT = IEOFT+1
   IF (IEOFT.GT.NIFILS) GO TO 10
   IF (IEOFSC.LT.0) GO TO 60
   GO TO 10
!
900 FORMAT (2F10.4,6I5,8A4,A3)
905 FORMAT ('1',/'   ***  FILTER ***',8(' ********** '))
910 FORMAT ('0   V1F=',F10.4,' V2F=',F10.4,',DVF=',F10.4,',NPTF =',   &
   &        I5,/,'0',10X,', JEMIT= ',I2,', JABS= ',I2,                &
   &        ', INPUT FILE= ',I3,' ,IFILST =',I5,' ,NIFILS =',I5,2X,   &
   &        8A4,A3)
915 FORMAT (A80)
920 FORMAT ('0  RESULT FROM SCANNING FUNCTION INCONSISTENT WITH ',    &
   &        'FILTER REQUEST')
925 FORMAT ('0',//,8(' -------  '),/,'0',10A8,2X,2(1X,A8,1X))
930 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)
935 FORMAT ('0 SECANT =',F15.5,/'0 PRESS(MB) =',F12.5/'0 TEMP(K) =',  &
   &        F11.2,/'0 DV(CM-1) = ',F12.8,/'0 V1(CM-1) = ',F12.6,/     &
   &        '0 V2(CM-1) = ',F12.6)
940 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//5X,'WBROAD = ',     &
   &        1PE10.3,/(5X,A6,' = ',1PE10.3))
945 FORMAT ('0   V1F=',F10.4,' V2F=',F10.4,',DVF=',F10.4,',NPTF =',   &
   &        I5,/'0 , IEMIT= ',I2,', JEMIT= ',I2,', JABS= ',I2,        &
   &        ', INPUT FILE= ',I3,' ,IFILST =',I5,' ,NIFILS =',I5,2X,   &
   &        8A4,A3)
950 FORMAT ('0',5X,'TIME =',F7.3,',  READ =',F6.3,',  CONV. =',F7.3)
955 FORMAT ('0  INTEGRATED TRANSMISSION = ',1PE14.5,                  &
   &        '  NORMALIZATION OF  THE FILTER = ',E14.5,/               &
   &        '0 UNNORMALIZED INTEGRATED TRANSMISSION =  ',E14.5)
960 FORMAT ('0  INTEGRATED ABSORPTION = ',1PE14.5,                    &
   &        '  NORMALIZATION OF  THE FILTER = ',E14.5,/               &
   &        '0 UNNORMALIZED INTEGRATED ABSORPTION =    ',E14.5)
!  965 FORMAT ('0 INTEGRATED EMISSION = ',1PE14.5,
!     *        '  NORMALIZATION OF THE',' FILTER = ',E14.5)
965 FORMAT ('0 INTEGRATED EMISSION = ',1PE14.5,                       &
   &        '  NORMALIZATION OF THE',' FILTER = ',1PE14.5,            &
   &        ' NORM. EMISSION',E14.5)

970 format (a4,i2.2)
975 format (1p,e14.5,1p,e14.5,1p,e14.5,0p)
980 format (' Filter output:')
!
end subroutine FLTRFN
!
!     --------------------------------------------------------------
!
SUBROUTINE FLTRRD
!
!     NFLTPT sets the maximum number of points in the incoming filter
!
   USE lblparams, ONLY: NFLTPT
!
   IMPLICIT REAL*8           (V)
!
!     READ CONTROL CARD FOR FILTER WITH WEIGHTING FUNCTIONS
!
   COMMON S(4650),R1(5750)
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DVC,V1C,V2C,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /RSCAN/ V1I,V2I,DVI,NNI
   COMMON /COMFLT/ V1F,V2F,DVF,NPTS,NPTF,JEMIT,IUNIT,IFILST,NIFILS,  &
   &                HEDDR(9),XF(NFLTPT),SUMFLT
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /FLFORM/ CFORM
   DIMENSION FILHDR(2)
!
   CHARACTER*80 CVAR
   CHARACTER CFORM*11,TAPE13*6,CTAPE*4
   LOGICAL OP
!
   EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(6),ISCAN) , (FSCDID(7),IPLOT),                &
   &            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                &
   &            (FSCDID(12),SCNID) , (FSCDID(13),HWHM),               &
   &            (FSCDID(16),LAYR1)
!
   DATA TAPE13 / '      '/,CTAPE / 'TAPE'/
!
   NLIMF = 2401
   NSHIFT = 0
   READ (IRD,900) V1F,DVF,NPTF,JEMIT,NNFILE,HEDDR
!
   JABS = 0
   IF (NNFILE.NE.NFILE.AND.NNFILE.GT.0) THEN
      INQUIRE (UNIT=NFILE,OPENED=OP)
      IF (OP) CLOSE (NFILE)
      NFILE = NNFILE
      INQUIRE (UNIT=NFILE,OPENED=OP)
      IF (.NOT.OP) THEN
         WRITE (TAPE13,905) CTAPE,NFILE
         OPEN (NFILE,FILE=TAPE13,STATUS='UNKNOWN',FORM=CFORM)
         REWIND NFILE
      ENDIF
   ENDIF
!
   IF (V1F.LT.0) RETURN
   WRITE (IPR,910)
!
   IF (NPTF.LE.0) GO TO 10
   NPTS = NPTF
10 V2F = V1F+DVF* REAL(NPTS-1)
   WRITE (IPR,915) V1F,V2F,DVF,NPTF,JEMIT,JABS,NFILE,HEDDR
   V1 = V1F
   V2 = V2F
   DV = DVF
   IF (NPTF.LE.0) GO TO 20
   READ (IRD,920) CVAR
   READ (IRD,CVAR) (XF(I),I=1,NPTS)
   WRITE (IPR,CVAR) (XF(I),I=1,NPTS)
!
!     MAKE ADJUSTMENT FOR END POINT CORRECTIONS
!
   XF(1) = 0.5*XF(1)
   XF(NPTS) = 0.5*XF(NPTS)
20 SUMFLT = 0.0
   DO 30 I = 1, NPTS
      SUMFLT = SUMFLT+XF(I)
30 END DO
   SUMFLT = SUMFLT*DVF
!
   RETURN
!
900 FORMAT (2F10.4,I5,I5,I5,10X,8A4,A3)
905 FORMAT (A4,I2.2)
910 FORMAT ('1',/'   ***  FILTER ***',8(' ********** '))
915 FORMAT ('0   V1F=',F10.4,' V2F=',F10.4,',DVF=',F10.4,             &
   &        ',NPTF =',I5,/,'0',10X,', JEMIT= ',I2,', JABS= ',I2,      &
   &        ',NIFILS =',I5,', OUTPUT FILE= ',I3,2X,8A4,A3)
920 FORMAT (A80)
!
end subroutine FLTRRD
!
!     --------------------------------------------------------------
!
SUBROUTINE FLTMRG (IFILE,JFILE)
!
!     NFLTPT sets the maximum number of points in the incoming filter
!
   USE lblparams, ONLY: NFLTPT
!
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE FLTMRG CALCULATES AND OUTPUTS THE RESULTS
!     OF THE FILTER TO FILE JFILE
!
   COMMON S(4650),R1(5750)
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DVC,V1C,V2C,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYMA,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID

   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /RSCAN/ V1I,V2I,DVI,NNI
   COMMON /COMFLT/ V1F,V2F,DVF,NPTS,NPTF,JEMIT,IUNIT,IFILST,NIFILS,  &
   &                HEDDR(9),XF(NFLTPT),SUMFLT
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /SPANEL/ V1P,V2P,DV,NLIM
   DIMENSION FILHDR(2),PNLHDR(2),RF(4)
   EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(6),ISCNHD) , (FSCDID(7),IPLOT),               &
   &            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                &
   &            (FSCDID(12),SCNID) , (FSCDID(13),HWHM),               &
   &            (FSCDID(16),LAYR1) , (PNLHDR(1),V1P)
!
   DATA I_1000/1000/
!
   NLIMF = 2401
   NREN = 0
   IPRT = 1
   NSHIFT = 0
   IUNIT = IFILE
   REWIND IFILE
!
10 RFILTR = 0.0
   VFT = V1
   VBOT = V1
   VTOP = V2
   TIMRDF = 0.0
   TIMCNV = 0.0
!
   CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)
!
   ISCAN = ISCNHD
   IF (ISCAN.LE.0.OR.SCNID.EQ.-99.) ISCAN = 0
   IF (ISCNHD.GE.1000.AND.ISCAN.EQ.0) ISCAN = ISCNHD
   ISCNHD = ISCAN+100
   IF (MOD(ISCAN,I_1000).EQ.0) GO TO 20
   JEMSCN = SCNID/100.
   IF (JEMIT.EQ.JEMSCN) GO TO 20
   WRITE (*,900)
   WRITE (IPR,900)
   CALL SKIPFL (1,IFILE,IEOF)
!
   IF (IEOF.EQ.0) RETURN
!
   GO TO 10
20 CONTINUE
!
   IF (IEOF.LT.1) RETURN
!
!     JEMIT=-1 FILTER PASSED OVER ABSORPTION
!     JEMIT=0  FILTER PASSED OVER TRANSMISSION
!     JEMIT=1  FILTER PASSED OVER EMISSION
!
   JTREM = -1
   IF ((IEMIT.EQ.0).AND.(JEMIT.EQ.0)) JTREM = 0
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.0)) JTREM = 2
   IF ((IEMIT.EQ.1).AND.(JEMIT.EQ.1)) JTREM = 1
   ISCANT = MOD(ISCAN,I_1000)
   IF ((ISCANT.GE.1).AND.(JEMIT.EQ.0)) JTREM = 2
!
   IF (JTREM.LT.0) RETURN
!
   WRITE (IPR,905) XID,(YID(M),M=1,2)
   WRITE (IPR,910) LAYR1,LAYER
   WRITE (IPR,915) SECANT,PAVE,TAVE,DVC,V1C,V2C
   WRITE (IPR,920) WBROAD,(HMOLID(M),WK(M),M=1,NMOL)
   WRITE (IPR,925) V1F,V2F,DVF,NPTF,IEMIT,JEMIT,JABS,IUNIT,IFILST,   &
   &                NIFILS,HEDDR
!
   XSCID = 100*JEMIT
   SCNID = XSCID+0.01
!
   CALL BUFOUT (JFILE,FILHDR(1),NFHDRF)
   IDATA = -1
30 CALL CPUTIM (TIMEO)
   CALL RDSCAN (S,JTREM,IFILE,ISCAN,IPRT)
   CALL CPUTIM (TIME)
   TIMRDF = TIMRDF+TIME-TIMEO
   IF (IEOFSC.NE.1) GO TO 40
   CALL CNVFLT (S,RFILTR,XF)
   IF (IDATA.EQ.1) GO TO 40
   GO TO 30
!
40 CALL CPUTIM (TIME)
   WRITE (IPR,930) TIME,TIMRDF,TIMCNV
   IF (JEMIT.EQ.1) GO TO 50
   TRNSM = RFILTR
   RFILTR = RFILTR*DVC/SUMFLT
   RF(1) = RFILTR
   RF(2) = SUMFLT
   RF(3) = TRNSM
   RF(4) = PLAY
   IF (JABS.EQ.0) WRITE (IPR,935) RFILTR,SUMFLT,TRNSM
   IF (JABS.EQ.1) WRITE (IPR,940) RFILTR,SUMFLT,TRNSM
!
!     V1P IS FIRST FREQ OF PANEL
!     V2P IS LAST  FREQ OF PANEL
!
   V1P = V1F
   V2P = V2F
   DVP = DVF
   NLIM = 4
   CALL BUFOUT (JFILE,PNLHDR(1),NPHDRF)
   CALL BUFOUT (JFILE,RF(1),NLIM)
   GO TO 60
!
50 RFILTR = RFILTR*DVC
   WRITE (IPR,945) RFILTR,SUMFLT
!
60 IF (IEOFSC.EQ.1) CALL SKIPFL (1,IFILE,IEOFSC)
   IEOFT = IEOFT+1
   IF (IEOFT.GT.NIFILS) RETURN
   IF (IEOFSC.LT.0) GO TO 10
!
   RETURN
!
900 FORMAT ('0  RESULT FROM SCANNING FUNCTION INCONSISTENT WITH ',    &
   &        'FILTER REQUEST')
905 FORMAT ('0',//,8(' -------  '),/,'0',10A8,2X,2(1X,A8,1X))
910 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)
915 FORMAT ('0 SECANT =',F15.5,/'0 PRESS(MB) =',F12.5/'0 TEMP(K) =',  &
   &        F11.2,/'0 DV(CM-1) = ',F12.8,/'0 V1(CM-1) = ',F12.6,/     &
   &        '0 V2(CM-1) = ',F12.6)
920 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//5X,'WBROAD = ',     &
   &        1PE10.3,/(5X,A6,' = ',1PE10.3))
925 FORMAT ('0   V1F=',F10.4,' V2F=',F10.4,',DVF=',F10.4,',NPTF =',   &
   &        I5,/,'0 , IEMIT= ',I2,', JEMIT= ',I2,', JABS= ',I2,       &
   &        ', INPUT FILE= ',I3,' ,IFILST =',I5,' ,NIFILS =',I5,2X,   &
   &        8A4,A3)
930 FORMAT ('0',5X,'TIME =',F7.3,',  READ =',F6.3,',  CONV. =',F7.3)
935 FORMAT ('0  INTEGRATED TRANSMISSION = ',1PE14.5,                  &
   &        '  NORMALIZATION OF  THE FILTER = ',E14.5,/               &
   &        '0 UNNORMALIZED INTEGRATED TRANSMISSION =  ',E14.5)
940 FORMAT ('0  INTEGRATED ABSORPTION = ',1PE14.5,                    &
   &        '  NORMALIZATION OF  THE FILTER = ',E14.5,/               &
   &        '0 UNNORMALIZED INTEGRATED ABSORPTION =    ',E14.5)
945 FORMAT ('0 INTEGRATED EMISSION = ',1PE14.5,                       &
   &        '  NORMALIZATION OF THE    FILTER = ',E14.5)
!
end subroutine FLTMRG
!
!     --------------------------------------------------------------
!
SUBROUTINE CNVFLT (S,RFILTR,XF)
!
!     NFLTPT sets the maximum number of points in the incoming filter
!
   USE lblparams, ONLY: NFLTPT
!
   IMPLICIT REAL*8           (V)
!
   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /RSCAN/ V1I,V2I,DVI,NNI
   COMMON /COMFLT/ V1F,V2F,DVF,NPTS,NPTF,JEMIT,IUNIT,IFILST,NIFILS,  &
   &                HEDDR(9),XFS(NFLTPT),SUMFLT
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
!
   DIMENSION XF(*),S(*)
!
   CALL CPUTIM (TIMEO)
   IMIN = (V1F-V1I)/DVI+1.0001
!      IMIN = (V1F-V1I)/DVI+1.5
   IMIN = MAX(IMIN,ILO)
   IMAX = (V2F+V1F-V1I)/DVI+1.0001
!      IMAX = (V2F+V1F-V1I)/DVI+1.5
   IMAX = MIN(IMAX,IHI)
   XIF0 = (V1I-V1F)/DVF+1.0001
!      XIF0 = (V1I-V1F)/DVF+1.5
   XDVIF = DVI/DVF
   DO 10 I = IMIN, IMAX
      IFL = XIF0+XDVIF* REAL(I)
!         IFL = XIF0+XDVIF* REAL(I-1)
!
!        Linearly interpolate filter function XF to avoid
!        discontinuities in output spectrum
!
      v1s = v1i+(i-1)*dvi
      v2s = v1i+i*dvi
      vxf1 = v1f + (ifl-1)*dvf
      vxf2 = v1f + (ifl)*dvf
      p = (vxf2-v2s)/(vxf2-vxf1)
      RFILTR = RFILTR+S(I)*(p*XF(IFL)+(1.-p)*xf(IFL+1))
10 END DO
   IF (IMAX.LT.IHI) VFT = VFT+(( REAL(IHI)- REAL(ILO))+1.0)*DVI
   CALL CPUTIM (TIME)
   TIMCNV = TIMCNV+TIME-TIMEO
!
   RETURN
!
end subroutine CNVFLT
!
!     --------------------------------------------------------------
!
SUBROUTINE FLTPRT (IFILE)
!
!     NFLTPT sets the maximum number of points in the incoming filter
!
   USE lblparams, ONLY: NFLTPT
!
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE FLTPRT READS FROM IFILE AND FORMATS OUT THE RESULTS
!     OF THE FILTERED WEIGHTING FUNCTION TO IPR
!
   COMMON S(4650),R1(5750)
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DVC,V1C,V2C,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI, &
   &               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,     &
   &               DVSC,XDUM,V1SHFT
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYMA,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID

   COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /RSCAN/ V1I,V2I,DVI,NNI
   COMMON /COMFLT/ V1F,V2F,DVF,NPTS,NPTF,JEMIT,IUNIT,IFILST,NIFILS,  &
   &                HEDDR(9),XF(NFLTPT),SUMFLT
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /SPANEL/ V1P,V2P,DV,NLIM
   DIMENSION FILHDR(2),PNLHDR(2),RF(4)
!
   EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(6),ISCAN) , (FSCDID(7),IPLOT),                &
   &            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                &
   &            (FSCDID(12),SCNID) , (FSCDID(13),HWHM),               &
   &            (FSCDID(16),LAYR1) , (PNLHDR(1),V1P)
!
   NLIMF = 2401
   NSHIFT = 0
   IUNIT = IFILE
   REWIND IFILE
!
   WRITE (IPR,900)
   OLDTRN = 1.0
   TIMRDF = 0.0
   CALL CPUTIM (TIMEO)
10 CONTINUE
   CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)
   IF (IEOF.EQ.-99) GO TO 10
   IF (IEOF.EQ.0) GO TO 20
!
   CALL BUFIN (IFILE,IEOF,PNLHDR(1),NPHDRF)
   CALL BUFIN (IFILE,IEOF,RF(1),NLIM)
!
   RFILTR = RF(1)
   SUMFLT = RF(2)
   TRNSM = RF(3)
   PLAY = RF(4)
   DTAU = OLDTRN-RFILTR
   OLDTRN = RFILTR
   WRITE (IPR,905) LAYER,PLAY,DTAU,RFILTR,TRNSM
!
   IF (IEOF.EQ.1) GO TO 10
20 WRITE (IPR,910) SUMFLT
   CALL CPUTIM (TIME)
   TIMRDF = TIMRDF+TIME-TIMEO
   WRITE (IPR,915) TIME,TIMRDF
!
   RETURN
!
900 FORMAT ('1',15X,                                                  &
   &        'FORMATTED RESULTS OF FILTERED WEIGHTING FUNCTIONS',      &
   &        2(/),10X,'LAYER',5X,'PRESSURE',4X,'TRANSMISSION',6X,      &
   &        'INTEGRATED',6X,'UNNORMALIZED',/,10X,'  #  ',5X,          &
   &        '  (MB)  ',5X,'(N-1) - (N)',5X,'TRANSMISSION',5X,         &
   &        '  INT TRANS ',/)
905 FORMAT (10X,I3,5X,F8.3,5X,1PE13.6,4X,1PE13.6,4X,1PE13.6)
910 FORMAT ('0',5X,'NORMALIZATION OF THE FILTER =',1PE13.6)
915 FORMAT ('0',5X,'TIME =',F7.3,',  IN FLTPRT =',F7.3)
!
end subroutine FLTPRT
