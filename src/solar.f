C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
C
C     ----------------------------------------------------------------
C
      SUBROUTINE SOLINT(IFILE,LFILE,NPTS,INFLAG,IOTFLG,JULDAT)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      
C               LAST MODIFICATION:    1 November 1995
C                                                                      
C                  IMPLEMENTATION:    P.D. Brown
C                                                                      
C             ALGORITHM REVISIONS:    S.A. Clough
C                                     P.D. Brown
C                                                                      
C                                                                      
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.      
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139       
C                                                                      
C----------------------------------------------------------------------
C                                                                      
C               WORK SUPPORTED BY:    THE ARM PROGRAM                  
C                                     OFFICE OF ENERGY RESEARCH        
C                                     DEPARTMENT OF ENERGY             
C                                                                      
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      PARAMETER (NSOL=2000001)
C
      IMPLICIT REAL*8           (V)
C
C     ------------------------------------------------------------
C     SUBROUTINE SOLINT interpolates solar radiances from the binary
C     file SOLAR.RAD.  The following are input and output options:
C
C       INFLAG = 0   => input transmittance from TAPE12 (default,
C                       where TAPE12 includes the monochromatic radiance
C                       and transmittance).
C              = 1   => input optical depths from TAPE12 and
C                       convert to transmittance.
C              = 2   => input R1, T1, T2, r1,  where
C
C                           _
C                           |
C                Observer |-O-|
C                           |
C                           -      \|/ Sun
C                          /|\    --O--
C                           |      /|\
C                           |   S2,/  /
C                           |   T2/  /
C                      R1,T1|    / |/_ R3
C                           |   /  /
C                           |  /  /
C                         r1||/_|/_r3
C            Ground    ---------------------------
C                     ///////////////////////////
C
C              = 3   => input transmittance from TAPE12 (CHARTS-type
C                       output)
C
C
C       IOTFLG = 0   => attenuate w/transmittance & output (default).
C              = 1   => attenuate and add to radiance from TAPE12
C                       (requires INFLAG = 1).
C              = 2   => Calculate solar contribution Rs = S2*T2*r1*T1
C                       and add to thermal contribution R1.
C
C     Output radiance goes to TAPE13.
C
C     ------------------------------------------------------------
C
      COMMON /MANE/ P0,TEMP0,NLAYER,DDUM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,
     *              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,
     *              EXTID(10)
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2
C
      character*8      XID,       HMOLID,      YID        
      real*8               SECANT,       XALTZ 
c
      character*8      XIDS,       HMLIDS,       YIDS        
      real*8                SECNTS,       XALTZS 
c
C
      COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,
     *               EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
      COMMON /SLHDR/ XIDS(10),SECNTS,PAVES,TAVES,HMLIDS(60),XALTZS(4),
     *               WKS(60),PZLS,PZUS,TZLS,TZUS,WBRODS,DVS,V1S,V2S,
     *               TBONDS,
     *               EMSIVS,FSCDDS(17),NMOLS,LAYERS,YI1S,YIDS(10),
     *               LSTWDS
      COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
      COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4
      COMMON /CVRSOL/ HVRSOL
C
C     ----------------------------------------------------------------
C     Parameter and common blocks for direct input of emission and
C     reflection function values
C
      PARAMETER (NMAXCO=4040)
      COMMON /EMSFIN/ V1EMIS,V2EMIS,DVEMIS,NLIMEM,ZEMIS(NMAXCO)
      COMMON /RFLTIN/ V1RFLT,V2RFLT,DVRFLT,NLIMRF,ZRFLT(NMAXCO)
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP
C     ----------------------------------------------------------------
C
      DIMENSION XFILHD(2),XSOLHD(2),PNLHDR(2),OPNLHD(2)
      DIMENSION A1(0:100),A2(0:100),A3(0:100),A4(0:100)
      DIMENSION A1T2(0:100),A2T2(0:100),A3T2(0:100),A4T2(0:100)
      DIMENSION A1RF(0:100),A2RF(0:100),A3RF(0:100),A4RF(0:100)
      DIMENSION TRAO(2),TRAN(2410)
      DIMENSION RADO(2),RADN(2410)
      DIMENSION OPTO(2),OPTN(2410)
C
      DIMENSION SOLAR(-1:NSOL)
      DIMENSION TRAN2(-1:NSOL)
      DIMENSION RAD2(-1:NSOL)
      DIMENSION XRFLT(-1:NSOL)
      DIMENSION SOLRAD(2410)
C
      CHARACTER*40 CYID
      CHARACTER*15 HVRSOL
C
      EQUIVALENCE (XSOLHD(1),XIDS(1))
      EQUIVALENCE (XFILHD(1),XID(1)),(PNLHDR(1),V1P),
     *            (OPNLHD(1),V1PO)
      EQUIVALENCE (TRAN(1),TRAO(1)),(RADN(1),RADO(1)),
     *            (OPTN(1),OPTO(1)),
     *            (FSCDID(4),IAERSL),(FSCDID(5),IEMIT),
     *            (FSCDID(6),ISCHDR),(FSCDID(7),IPLOT),
     *            (FSCDID(8),IPATHL),(FSCDID(12),XSCID),
     *            (FSCDID(16),LAYR1)
C
C     ************************************************************
C     ****** THIS PROGRAM DOES MERGE FOR SOLAR RADIANCE AND ******
C     ****** TRANMITTANCE USING FOUR POINT INTERPOLATION    ******
C     ************************************************************
C
C
C     ASSIGN CVS VERSION NUMBER TO MODULE 
C
      HVRSOL = '$Revision$'
C
C     -------------------
C
C     Open file SOLAR.RAD
C
      ISOLFL = 19
      OPEN(UNIT=ISOLFL,FILE='SOLAR.RAD',FORM='UNFORMATTED',
     *     STATUS='OLD')
C
      CALL CPUTIM (TIME)
      WRITE (IPR,900) TIME
      NPANLS = 0
      TIMRD = 0.0
      TIMOT = 0.0

C     Calculate Earth distance to sun given Julian Day JULDAT. Used to
C     scale solar source function. Formula taken from "Atmospheric Radiative
C     Transfer", J. Lenoble, 1993.

c     Test validity of JULDAT

      if ((juldat.lt.0).or.(juldat.gt.366)) then
         write(*,*) 'JULDAT = ',juldat,' is out of range 0-366.'
         write(ipr,*) 'JULDAT = ',juldat,' is out of range 0-366.'
         stop 'Stopped in SOLINT'
      endif

c     If JULDAT = 0 , then set XJUL_SCALE to 1

      if (juldat.eq.0) then
         XJUL_SCALE = 1.0
         write(ipr,*) 'JULDAT = 0, no scaling of solar source function'
      else
         theta = 2*pi*(float(JULDAT)-1.)/365.
         XJUL_SCALE = 1.00011 + 0.034221*cos(theta) + 
     *        1.28E-3*sin(theta) + 7.19E-4*cos(2.*theta) + 
     *        7.7E-5*sin(2.*theta)
         write(ipr,*) 'JULDAT = ',JULDAT,
     *        ', scale factor for solar source function = ',XJUL_SCALE
      endif


C
C     FOR AEROSOL RUNS, MOVE YID (IFILE) INTO YID (LFILE)
C
C     Read file header of solar radiance file and determine dv ratio
C
      IF (IAERSL.GT.0) WRITE (CYID,'(5A8)') (YID(I),I=3,7)
      CALL BUFIN (ISOLFL,KEOF,XFILHD(1),NFHDRF)
cpdb      IF (IAERSL.GT.0) READ (CYID,'(5A8)') (YID(I),I=3,7)
cpdb      IF (JPATHL.GE.1) IPATHL = JPATHL
      DVK = DV
C
C     Read in file header of transmittance/optical depth file
C
      CALL BUFIN (IFILE,LEOF,XFILHD(1),NFHDRF)
      DVL = DV
      DVMIN = DVL
C
      ATYPE = 9.999E09
      IF (DVK.EQ.DVL) ATYPE = 0.
      IF (DVL.GT.DVK) ATYPE = DVK/(DVL-DVK)+0.5
      IF (DVL.LT.DVK) ATYPE = -DVL/(DVK-DVL)-0.5
C
C     IF (ATYPE .GT. 0) STOP  ' SOLINT; ATYPE GT 0 '
C
C
C     Write file information out to TAPE6
C     If INFLAG = 2, then also open files and read headers
C     for downward transmittance from solar path and solar
C     reflectance at the ground.
C
      WRITE (IPR,905) ISOLFL,IFILE,LFILE,ATYPE,INFLAG,IOTFLG
      IF (INFLAG.EQ.0) THEN
         WRITE(IPR,920) IFILE
      ELSEIF (INFLAG.EQ.1) THEN
         WRITE(IPR,925) IFILE
      ELSEIF (INFLAG.EQ.2) THEN
         ISLTRN = 20
         OPEN(UNIT=ISLTRN,FILE='SOL.PATH.T2',FORM='UNFORMATTED',
     *        STATUS='OLD')
         CALL BUFIN (ISLTRN,ITEOF,XSOLHD(1),NFHDRF)
         DVT2 = DVS
         ISLRFL = 21
         OPEN(UNIT=ISLRFL,FILE='SOL.REFLECTANCE',FORM='UNFORMATTED',
     *        STATUS='OLD')
         CALL BUFIN (ISLRFL,IREOF,XSOLHD(1),NFHDRF)
         DVRF = DVS
C
C        **************************************************************
C
C        Read Record 1.4 from TAPE5
C
         READ (IRD,970,END=80) TMPBND,(BNDEMI(IBND),IBND=1,3),
     *                         (BNDRFL(IBND),IBND=1,3)
C
         BNDTST = ABS(BNDRFL(1))+ABS(BNDRFL(2))+ABS(BNDRFL(3))
         IF (BNDTST.NE.0.) IBPROP = 1
C
C        If BNDEMI(1) < 0, read in coefficients from file 'EMISSION'
C        If BNDEMI(1) > 0, check to see if emissivity is reasonable
C
C        UNIT ICOEF used for input files
C
         ICOEF = 14
C         
         IF (BNDEMI(1).LT.0) THEN
            OPEN (UNIT=ICOEF,FILE='EMISSION',STATUS='OLD')
            CALL READEM(ICOEF)
            CLOSE (ICOEF)
         ELSE
            XVMID = (V1+V2)/2.
            EMITST = BNDEMI(1)+BNDEMI(2)*XVMID+BNDEMI(3)*XVMID*XVMID
            IF (EMITST.LT.0..OR.EMITST.GT.1.) THEN
               WRITE (IPR,975) XVMID,EMITST
               STOP 'BNDEMI'
            ENDIF
         ENDIF
C
C        ------------------------------------------------------------
C        *** NOTE: REFLECTION FUNCTION NOT CURRENTLY INCORPORATED ***
C        *** INTO SOLAR RADIATIVE TRANSFER CALCULATION            ***
C
C        If BNDRFL(1) < 0, read in coefficients from file 'REFLECTION'
C        If BNDRFL(1) > 0, check to see if reflectivity is reasonable
C     
         IF (BNDRFL(1).LT.0) THEN
            OPEN (UNIT=ICOEF,FILE='REFLECTION',STATUS='OLD')
            CALL READRF(ICOEF)
            CLOSE (ICOEF)
         ELSE
            REFTST = BNDRFL(1)+BNDRFL(2)*XVMID+BNDRFL(3)*XVMID*XVMID
            IF (REFTST.LT.0..OR.REFTST.GT.1.) THEN
               WRITE (IPR,980) XVMID,REFTST
               STOP 'BNDRFL'
            ENDIF
            DVEMIS = MIN(DVL,DVK,DVT2,DVRF)
         ENDIF
C        ------------------------------------------------------------
C
C        TBOUND is the boundary temperature. TBOUND=0. for no boundary
C        EMISIV is the boundary emissivity
C        Set default for EMISIV
C
         EMITST = ABS(BNDEMI(1))+ABS(BNDEMI(2))+ABS(BNDEMI(3))
         IF ((TMPBND.GT.0.).AND.(EMITST.EQ.0.)) BNDEMI(1) = 1.
         EMISIV = BNDEMI(1)
         TBOUND = TMPBND
         XKTBND = TBOUND/RADCN2
         WRITE (IPR,985) V1,V2,TBOUND,(BNDEMI(IBND),IBND=1,3),
     *                   (BNDRFL(IBND),IBND=1,3)
C
C     **************************************************************
C
C        Determine the minimum and maximum DV of all files and reset
C        ATYPE (-1. is only a flag for nonzero ATYPE)
C     
         DVMIN = MIN(DVL,DVK,DVT2,DVRF,DVEMIS)
         DVMAX = MAX(DVL,DVK,DVT2,DVRF,DVEMIS)
c         IF (DVMAX.EQ.DVMIN) THEN
c            ATYPE = 0.0
c         ELSE
            ATYPE = -1.
c         ENDIF
         WRITE(IPR,927) IFILE,ISLTRN,ISLRFL
         WRITE(IPR,941) DVT2,DVRF
C
      ELSEIF (INFLAG.EQ.3) THEN
         WRITE(IPR,926) IFILE
      ENDIF
      WRITE(IPR,906) DVL,DVK
C
C     Test for INFLAG and IOTFLG compatibility
C
      IF (IOTFLG.EQ.0) THEN
         WRITE(IPR,930) LFILE
      ELSEIF (IOTFLG.EQ.1) THEN
         IF (INFLAG.EQ.1) STOP 'ERROR: INFLAG=1, IOTFLG=1'
         WRITE(IPR,935) LFILE
      ELSEIF (IOTFLG.EQ.2) THEN
         IF (INFLAG.NE.2) STOP 'ERROR: INFLAG not 2, IOTFLG=2'
         WRITE(IPR,940) LFILE
      ENDIF
      IEMIT = 1
      SECANT = 0.
      DV = DVL
C
C     Set XSCID and ISCHDR to values which allow for potential 
C     interpolation of output file
C
      XSCID = 100.01
      ISCHDR = 10
C
C     Output file header
C
      CALL BUFOUT (LFILE,XFILHD(1),NFHDRF)
C
      IF (ATYPE.EQ.0.) THEN
C
C     ======================================
C        1/1 ratio only
C     ======================================
C
         WRITE(IPR,950) DVMIN
         IPANEM = 0
   30    CONTINUE
         IPANEM = IPANEM+1
         CALL CPUTIM (TIMSL1)
         CALL SOLIN (V1P,V2P,DVP,NLIM,ISOLFL,SOLAR(1),LSEOF)
         CALL CPUTIM (TIMSL2)
         TIMRD = TIMRD+TIMSL2-TIMSL1
         IF (LSEOF.LE.0) GO TO 110
C
C
C        If INFLAG = 0, then read radiance and transmittance
C        If INFLAG = 1, then read optical depth
C        If INFLAG = 2, then read radiance and transmittance
C                       and call SOLIN to read in r1 and T2
C        If INFLAG = 3, then read transmittance
C
         IF (INFLAG.EQ.0) THEN
            CALL SOLIN2 (V1PO,V2PO,DVPO,NLIMO,IFILE,RADO(1),
     *           TRAO(1),LEOF)
         ELSEIF (INFLAG.EQ.1) THEN
            CALL SOLIN (V1PO,V2PO,DVPO,NLIMO,IFILE,OPTO(1),
     *           LEOF)
            DO 35 I = 1,NLIMO
               TRAN(I) = EXP(-OPTN(I))
 35         CONTINUE
         ELSEIF (INFLAG.EQ.2) THEN
            CALL SOLIN2 (V1PO,V2PO,DVPO,NLIMO,IFILE,RADO(1),
     *           TRAO(1),LEOF)
            CALL SOLIN2 (V1T2,V2T2,DVT2,NLIMT2,ISLTRN,RAD2(1),
     *           TRAN2(1),LSEOF)
            CALL SOLIN (V1RF,V2RF,DVRF,NLIMRF,ISLRFL,XRFLT(1),
     *           LSEOF)
            IF (ABS(V1T2-V1RF).GT.0.001) THEN 
               WRITE(IPR,*) 'SOLINT: PANELS DO NOT MATCH:'
               WRITE(IPR,*) '  V1T2 = ',V1T2,'  V1RF = ',V1RF
               STOP 'SOLINT: PANELS DO NOT MATCH: SEE TAPE6'
            ENDIF
         ELSEIF (INFLAG.EQ.3) THEN
            CALL SOLIN (V1PO,V2PO,DVPO,NLIMO,IFILE,TRAO(1),
     *           LEOF)
         ENDIF
         CALL CPUTIM (TIMSL3)
         TIMRD = TIMRD+TIMSL3-TIMSL2
C
C        If IOTFLG = 0, then calculate attenuated solar radiance
C        If IOTFLG = 1, then calculate attenuated solar radiance
C                       plus atmospheric radiance
C        If IOTFLG = 2, then calculate attenuated solar radiance
C                       through the reflected atmosphere plus
C                       atmospheric radiance
C
C        Solar irradiance is input from SOLAR.RAD (W/m2 cm-1).
C        Convert to radiance units (W/cm2 sr cm-1) by multiplying
C        by 1/6.8e-5.

         conv_ster = 1./(1.e4*6.8e-5)

C
C        Combine XJUL_SCALE and conv_ster into one scale factor SCAL_FAC

         SCAL_FAC = conv_ster*XJUL_SCALE
C
         IF (IOTFLG.EQ.0) THEN
            DO 40 I = 1, NLIM
               SOLRAD(I) = SOLAR(I)*SCAL_FAC*TRAN(I)
 40         CONTINUE
         ELSEIF (IOTFLG.EQ.1) THEN
            DO 41 I = 1, NLIM
               SOLRAD(I) = SOLAR(I)*SCAL_FAC*TRAN(I)+RADN(I)
 41         CONTINUE
         ELSEIF (IOTFLG.EQ.2) THEN
            IF (TBOUND.EQ.0) THEN
               DO 42 I = 1, NLIM
                  SOLRAD(I) = SOLAR(I)*SCAL_FAC*TRAN2(I)*XRFLT(I)*
     *                 TRAN(I)+RADN(I)
 42            CONTINUE
            ELSE
               DO 43 I = 1, NLIM
                  EMDUM = -1.
                  BBDUM = -1.
                  VBND = V1PO+(I-1)*DVPO
                  ZEMSV = EMISFN(VBND,DVPO,VINEM,EMDEL,EMDUM)
                  BBND = BBFN(VBND,DVPO,VBND,XKTBND,VINEW,BBDEL,BBDUM)
                  SOLRAD(I) = (SOLAR(I)*SCAL_FAC*TRAN2(I)*XRFLT(I)
     *                 +ZEMSV*BBND)*TRAN(I)+RADN(I)
 43            CONTINUE
            ENDIF
         ENDIF
C
         CALL CPUTIM (TIMSL2)
         CALL SOLOUT(V1PO,V2PO,DVPO,NLIMO,SOLRAD,LFILE,NPTS,NPANLS)
         CALL CPUTIM (TIMSL3)
         TIMOT = TIMOT+TIMSL3-TIMSL2
         GO TO 30
C
      ENDIF
C
C     ======================================
C     All ratios except 1/1
C     ======================================
C
      WRITE(IPR,951) DVMIN
      DO 50 JP = 0,100
         APG = JP
         P = 0.01*APG
C
C        The following are the constants for the Lagrange
C        4 point interpolation
C
         A1(JP) = -P*(P-1.0)*(P-2.0)/6.0
         A2(JP) = (P**2-1.0)*(P-2.0)*0.5
         A3(JP) = -P*(P+1.0)*(P-2.0)*0.5
         A4(JP) = P*(P**2-1.0)/6.0
         A1T2(JP) = A1(JP)
         A2T2(JP) = A2(JP)
         A3T2(JP) = A3(JP)
         A4T2(JP) = A4(JP)
         A1RF(JP) = A1(JP)
         A2RF(JP) = A2(JP)
         A3RF(JP) = A3(JP)
         A4RF(JP) = A4(JP)
   50 CONTINUE
C
C     *** Beginning of loop that does merge  ***
C
      NPE = 0
      V1P = 0.0
      V2P = 0.0
      DVP = 0.0
      SOLAR(0) = 0.0
      LSEOF = 1
      NPT2 = 0
      V1T2 = 0.0
      V2T2 = 0.0
      DVT2 = 0.0
      TRAN2(0) = 0.0
      LT2EOF = 1
      NPRF = 0
      V1RF = 0.0
      V2RF = 0.0
      DVRF = 0.0
      XRFLT(0) = 0.0
      LRFEOF = 1
      V1PO = 0.0
      V2PO = 0.0
      DVPO = 0.0
C
C     ============================================================
C
      IPANEM = 0
   60 CONTINUE
      IPANEM = IPANEM+1
C     
C     If INFLAG = 0, then read radiance and transmittance
C     If INFLAG = 1, then read optical depth
C     If INFLAG = 2, then read radiance and transmittance
C                    and call SOLIN to read in r1 and T2
C     If INFLAG = 3, then read transmittance
C
      IF (INFLAG.EQ.0) THEN
         CALL SOLIN2 (V1PO,V2PO,DVPO,NLIMO,IFILE,RADO(1),
     *        TRAO(1),LEOF)
         IF (LEOF.LE.0) GO TO 110
      ELSEIF (INFLAG.EQ.1) THEN
         CALL SOLIN (V1PO,V2PO,DVPO,NLIMO,IFILE,OPTO(1),
     *        LEOF)
         IF (LEOF.LE.0) GO TO 110
         DO 65 I = 1,NLIMO
            TRAN(I) = EXP(-OPTN(I))
 65      CONTINUE
      ELSEIF (INFLAG.EQ.2) THEN
         CALL SOLIN2 (V1PO,V2PO,DVPO,NLIMO,IFILE,RADO(1),
     *        TRAO(1),LEOF)
         IF (LEOF.LE.0) GO TO 110
C
C        -----------------------------------------------------
C        TRAN2 and RAD2 read in
C
         IF (V2T2.LE.V2PO+DVT2.AND.LT2EOF.GT.0) THEN
 66         CALL CPUTIM(TIMSL2)
            CALL SOLIN2 (V1T2,V2T2,DVT2,NLIMT2,ISLTRN,RAD2(npt2+1),
     *           TRAN2(npt2+1),LT2EOF)
            CALL CPUTIM(TIMSL3)
            TIMRD = TIMRD+TIMSL3-TIMSL2
            IF (LT2EOF.LE.0) GO TO 67
            V1T2 = V1T2-FLOAT(NPT2)*DVT2
 1          NPT2 = NLIMT2+NPT2
            IF (V2T2.LE.V2PO+DVT2) GO TO 66
         ENDIF
C
C        Zero point of first panel
C
 67      IF (TRAN2(0).EQ.0.0) THEN
            TRAN2(-1) = TRAN2(1)
            TRAN2(0) = TRAN2(1)
         ENDIF
C
C        End point of last panel
C
         IF (V2T2+DVT2.GE.V2) THEN
            TRAN2(NPT2+1) = TRAN2(NPT2)
            TRAN2(NPT2+2) = TRAN2(NPT2)
         ENDIF
C
C        -----------------------------------------------------
C        XRFLT read in
C
         IF (V2RF.LE.V2PO+DVRF.AND.LRFEOF.GT.0) THEN
 68         CALL CPUTIM(TIMSL2)
            CALL SOLIN (V1RF,V2RF,DVRF,NLIMRF,ISLRFL,XRFLT(nprf+1),
     *           LRFEOF)
            CALL CPUTIM(TIMSL3)
            TIMRD = TIMRD+TIMSL3-TIMSL2
            IF (LRFEOF.LE.0) GO TO 69
            V1RF = V1RF-FLOAT(NPRF)*DVRF 
            NPRF = NLIMRF+NPRF
            IF (V2RF.LE.V2PO+DVRF) GO TO 68
         ENDIF
C
C        Zero point of first panel
C
 69      IF (XRFLT(0).EQ.0.0) THEN
            XRFLT(-1) = XRFLT(1)
            XRFLT(0) = XRFLT(1)
         ENDIF
C
C        End point of last panel
C
         IF (V2RF+DVRF.GE.V2) THEN
            XRFLT(NPRF+1) = XRFLT(NPRF)
            XRFLT(NPRF+2) = XRFLT(NPRF)
         ENDIF
      ELSEIF (INFLAG.EQ.3) THEN
         CALL SOLIN (V1PO,V2PO,DVPO,NLIMO,IFILE,TRAO(1),
     *        LEOF)
         IF (LEOF.LE.0) GO TO 110
      ENDIF
C     -----------------------------------------------------
C
      CALL CPUTIM(TIMSL3)
      TIMRD = TIMRD+TIMSL3-TIMSL2
      II = 1
C
C     Buffer in panels from solar radiance file
C
      IF (V2P.LE.V2PO+DVP .AND.LSEOF.GT.0) THEN
   70    CALL CPUTIM(TIMSL2)
         CALL SOLIN(V1P,V2P,DVP,NLIM,ISOLFL,SOLAR(NPE+1),LSEOF)
         CALL CPUTIM(TIMSL3)
         TIMRD = TIMRD+TIMSL3-TIMSL2
         IF (LSEOF.LE.0) GO TO 80
         V1P = V1P-FLOAT(NPE)*DVP
         NPE = NLIM+NPE
         IF (V2P.LE.V2PO+DVP) GO TO 70
      ENDIF
C
C     Zero point of first panel
C
   80 IF (SOLAR(0).EQ.0.0) THEN
         SOLAR(-1) = SOLAR(1)
         SOLAR(0) = SOLAR(1)
      ENDIF
C
C     End point of last panel
C
      IF (V2P+DVP.GE.V2) THEN
         SOLAR(NPE+1) = SOLAR(NPE)
         SOLAR(NPE+2) = SOLAR(NPE)
      ENDIF
C
C     -----------------------------------------------------
C
C     NPL is the location of first element on arrays RADO and TRAO
C
      NPL = 1
      NPLT2 = 1
      NPLRF = 1
C
C     Set spectral spacing ratios (uses DVMIN instead of DVL, for
C     coding where it is not assumed that RADN & TRAN have the
C     smallest spacing).
C
      RATDVR = DVMIN/DVK
      IF (INFLAG.EQ.2) THEN
         RTDVT2 = DVMIN/DVT2
         RTDVRF = DVMIN/DVRF
      ENDIF
C
C     Test for ratios of 1 (no need for interpolation).  Reset
C     interpolation coefficients if necessary.
C
      IF (RATDVR.EQ.1) THEN
         DO 82 JP = 0,100
            A1(JP) = 0.0
            A2(JP) = 1.0
            A3(JP) = 0.0
            A4(JP) = 0.0
 82      CONTINUE
      ENDIF
      IF (INFLAG.EQ.2) THEN
         IF (RTDVT2.EQ.1) THEN
            DO 84 JP = 0,100
               A1T2(JP) = 0.0
               A2T2(JP) = 1.0
               A3T2(JP) = 0.0
               A4T2(JP) = 0.0
 84         CONTINUE
         ENDIF
         IF (RTDVRF.EQ.1) THEN
            DO 86 JP = 0,100
               A1RF(JP) = 0.0
               A2RF(JP) = 1.0
               A3RF(JP) = 0.0
               A4RF(JP) = 0.0
 86         CONTINUE
         ENDIF
      ENDIF
C
C     -----------------------------------------------------
C
C     FJJ is offset by 2. (for rounding purposes)
C
      FJ1DIF = (V1PO-V1P)/DVP+1.+2.
      IF (INFLAG.EQ.2) THEN
         FJDFT2 = (V1PO-V1T2)/DVT2+1.+2.
         FJDFRF = (V1PO-V1RF)/DVRF+1.+2.
      ENDIF
C
C     ============================================================
C
C     ***** Beginning of loop that does merge  *****
C
C     If IOTFLG = 0, then calculate attenuated solar radiance
C     If IOTFLG = 1, then calculate attenuated solar radiance
C                    plus atmospheric radiance
C     If IOTFLG = 2, then calculate attenuated solar radiance
C                    through the reflected atmosphere plus
C                    atmospheric radiance
C
C     Solar irradiance is input from SOLAR.RAD (W/m2 cm-1).
C     Convert to radiance units (W/cm2 sr cm-1) by multiplying
C     by 1/6.8e-5.

      conv_ster = 1./(1.e4*6.8e-5)

C
C     Combine XJUL_SCALE and conv_ster into one scale factor SCAL_FAC
      
      SCAL_FAC = conv_ster*XJUL_SCALE
C

      IF (IOTFLG.EQ.0) THEN
         DO 90 II = 1, NLIMO
            FJJ = FJ1DIF+RATDVR*FLOAT(II-1)
            JJ = IFIX(FJJ)-2
            JP = (FJJ-FLOAT(JJ))*100.-199.5
            SOLRAD(II) = (A1(JP)*SOLAR(JJ-1)+A2(JP)*SOLAR(JJ)+
     *           A3(JP)*SOLAR(JJ+1)+A4(JP)*SOLAR(JJ+2))*SCAL_FAC*
     *           TRAN(II)
c
 90      CONTINUE
      ELSEIF (IOTFLG.EQ.1) THEN
         DO 91 II = 1, NLIMO
            FJJ = FJ1DIF+RATDVR*FLOAT(II-1)
            JJ = IFIX(FJJ)-2
            JP = (FJJ-FLOAT(JJ))*100.-199.5
            SOLRAD(II) = (A1(JP)*SOLAR(JJ-1)+A2(JP)*SOLAR(JJ)+
     *           A3(JP)*SOLAR(JJ+1)+A4(JP)*SOLAR(JJ+2))*SCAL_FAC*
     *           TRAN(II)+RADN(II)
 91      CONTINUE
      ELSEIF (IOTFLG.EQ.2) THEN
         IF (TBOUND.EQ.0.) THEN
            DO 92 II = 1, NLIMO
               FJJ = FJ1DIF+RATDVR*FLOAT(II-1)
               FJJT2 = FJDFT2+RTDVT2*FLOAT(II-1)
               FJJRF = FJDFRF+RTDVRF*FLOAT(II-1)
               JJ = IFIX(FJJ)-2
               JJT2 = IFIX(FJJT2)-2
               JJRF = IFIX(FJJRF)-2
               JP = (FJJ-FLOAT(JJ))*100.-199.5
               JPT2 = (FJJT2-FLOAT(JJT2))*100.-199.5
               JPRF = (FJJRF-FLOAT(JJRF))*100.-199.5
               ZSOL = (A1(JP)*SOLAR(JJ-1)+A2(JP)*SOLAR(JJ)+
     *              A3(JP)*SOLAR(JJ+1)+A4(JP)*SOLAR(JJ+2))*SCAL_FAC
               ZTR2 = (A1T2(JPT2)*TRAN2(JJT2-1)+
     *              A2T2(JPT2)*TRAN2(JJT2)+
     *              A3T2(JPT2)*TRAN2(JJT2+1)+
     *              A4T2(JPT2)*TRAN2(JJT2+2))
               ZRFL = (A1RF(JPRF)*XRFLT(JJRF-1)+
     *              A2RF(JPRF)*XRFLT(JJRF)+
     *              A3RF(JPRF)*XRFLT(JJRF+1)+
     *              A4RF(JPRF)*XRFLT(JJRF+2))
               SOLRAD(II) = ZSOL*ZTR2*ZRFL*TRAN(II)+RADN(II)
 92            CONTINUE
         ELSE
            DO 93 II = 1, NLIMO
               FJJ = FJ1DIF+RATDVR*FLOAT(II-1)
               FJJT2 = FJDFT2+RTDVT2*FLOAT(II-1)
               FJJRF = FJDFRF+RTDVRF*FLOAT(II-1)
               JJ = IFIX(FJJ)-2
               JJT2 = IFIX(FJJT2)-2
               JJRF = IFIX(FJJRF)-2
               JP = (FJJ-FLOAT(JJ))*100.-199.5
               JPT2 = (FJJT2-FLOAT(JJT2))*100.-199.5
               JPRF = (FJJRF-FLOAT(JJRF))*100.-199.5
               ZSOL = (A1(JP)*SOLAR(JJ-1)+A2(JP)*SOLAR(JJ)+
     *              A3(JP)*SOLAR(JJ+1)+A4(JP)*SOLAR(JJ+2))*SCAL_FAC
               ZTR2 = (A1T2(JPT2)*TRAN2(JJT2-1)+
     *              A2T2(JPT2)*TRAN2(JJT2)+
     *              A3T2(JPT2)*TRAN2(JJT2+1)+
     *              A4T2(JPT2)*TRAN2(JJT2+2))
               ZRFL = (A1RF(JPRF)*XRFLT(JJRF-1)+
     *              A2RF(JPRF)*XRFLT(JJRF)+
     *              A3RF(JPRF)*XRFLT(JJRF+1)+
     *              A4RF(JPRF)*XRFLT(JJRF+2))
               EMDUM = -1.
               BBDUM = -1.
               VBND = V1PO+(II-1)*DVPO
               ZEMSV = EMISFN(VBND,DVPO,VINEM,EMDEL,EMDUM)
               BBND = BBFN(VBND,DVPO,VBND,XKTBND,VINEW,BBDEL,BBDUM)
               SOLRAD(II) = (ZSOL*ZTR2*ZRFL+ZEMSV*BBND)*
     *              TRAN(II)+RADN(II)
 93         CONTINUE
         ENDIF
      ENDIF
C
      NPL = JJ-1
      NPLT2 = JJT2-1
      NPLRF = JJRF-1
C
      CALL CPUTIM (TIMSL1)
C
C     ============================================================
C
C     Output attenuated radiance
C
      CALL SOLOUT(V1PO,V2PO,DVPO,NLIMO,SOLRAD,LFILE,NPTS,NPANLS)
      CALL CPUTIM (TIMSL2)
      TIMOT = TIMOT+TIMSL2-TIMSL1
C
C     ============================================================
C
C     Reset element locations for each array
C
C     ============================================================
C     NPL is now location of first element in the array SOLAR to
C     be used for next pass.
C
      IPL = -2
      DO 100 NL = NPL, NPE
         IPL = IPL+1
         SOLAR(IPL) = SOLAR(NL)
  100 CONTINUE
C
      V1P = V1P+FLOAT(NPL+1)*DVP
      NPE = IPL
C
      IF (IOTFLG.EQ.2) THEN
C
C        ------------------------------------------------------------
C        NPLT2 is now location of first element in the array TRAN2 to
C        be used for next pass.
C
         IPL = -2
         DO 102 NL = NPLT2, NPT2
            IPL = IPL+1
            TRAN2(IPL) = TRAN2(NL)
 102     CONTINUE
C     
         V1T2 = V1T2+FLOAT(NPLT2+1)*DVT2
         NPT2 = IPL
C     
C        -------------------------------------------------------------
C        NPLRF is now location of first element in the array XRFLT to
C        be used for next pass.
C
         IPL = -2
         DO 104 NL = NPLRF, NPRF
            IPL = IPL+1
            XRFLT(IPL) = XRFLT(NL)
 104     CONTINUE
C     
         V1RF = V1RF+FLOAT(NPLRF+1)*DVRF
         NPRF = IPL
      ENDIF
C     ============================================================
C
      GO TO 60
  110 CONTINUE
C
      CALL CPUTIM (TIME1)
      TIM = TIME1-TIME
      WRITE (IPR,910) TIME1,TIM,TIMRD,TIMOT
C
      RETURN
C
  900 FORMAT ('0 THE TIME AT THE START OF SOLINT IS ',F12.3)
  905 FORMAT ('0 FILE ',I5,' MERGED WITH FILE ',I5,' ONTO FILE',
     *        I5,'  WITH XTYPE=',G15.5,/,'0 INFLAG = ',I5,4X,
     *        'IOTFLG = ',I5)
 906  FORMAT ('0          Thermal spectrum spacing = ',E10.5,/,
     *        '0   Solar radiance spectral spacing = ',E10.5,/)
  910 FORMAT ('0 THE TIME AT THE END OF SOLINT IS ',F12.3/F12.3,
     *        ' SECS WERE REQUIRED FOR THIS SOLAR MERGE',F12.3,
     *        ' - READ - ',F12.3,' - SOLOUT - ',F12.3)
 920  FORMAT ('0 Radiance and Transmittance read in from unit',I5)
 925  FORMAT ('0 Optical Depths read in from unit',I5)
 926  FORMAT ('0 Transmittance read in from unit',I5)
 927  FORMAT ('0 Thermal upwelling Radiance & Transmittance',
     *        ' read in from unit',I5,/,
     *        '  Solar reflectance function read in from unit',
     *        I5,/,
     *        '  Thermal downwelling Radiance & Transmittance',
     *        ' read in from unit',I5,/)
 930  FORMAT ('0 Attenuated solar radiance output to unit',I5,/)
 935  FORMAT ('0 Attenuated solar radiance + atmospheric radiance',
     *        1x,'output to unit',I5,/)
 940  FORMAT ('0 Attenuated solar radiance + atmospheric radiance',
     *        1x,'(including effects of reflection) output to unit',
     *        I5,/)
 941  FORMAT ('0       Thermal downwelling spacing = ',E10.5,/,
     *        '0       Reflection function spacing = ',E10.5)
 950  FORMAT ('0 No interpolation needed: using DV = ',E10.5,//)
 951  FORMAT ('0 Interpolating to spectral spacing = ',E10.5,//)
 970  FORMAT (7E10.3)
 975  FORMAT ('0 FOR VNU = ',F10.3,' THE EMISSIVITY = ',E10.3,
     *        ' AND IS NOT BOUNDED BY (0.,1.) ')
 980  FORMAT ('0 FOR VNU = ',F10.3,' THE REFLECTIVITY = ',E10.3,
     *        ' AND IS NOT BOUNDED BY (0.,1.) ')
 985  FORMAT (5(/),'0*********** BOUNDARY PROPERTIES ***********',/,
     *        '0 V1(CM-1) = ',F12.4,/,'0 V2(CM-1) = ',F12.4,/,
     *        '0 TBOUND   = ',F12.4,5X,'BOUNDARY EMISSIVITY   = ',
     *        3(1PE11.3),/,'0',29X,'BOUNDARY REFLECTIVITY = ',
     *        3(1PE11.3))

C
      END
C
C     ----------------------------------------------------------------
C
      SUBROUTINE SOLIN (V1P,V2P,DVP,NLIM,KFILE,XARRAY,KEOF)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      
C               LAST MODIFICATION:    1 November 1995
C                                                                      
C                  IMPLEMENTATION:    P.D. Brown
C                                                                      
C             ALGORITHM REVISIONS:    S.A. Clough
C                                     P.D. Brown
C                                                                      
C                                                                      
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.      
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139       
C                                                                      
C----------------------------------------------------------------------
C                                                                      
C               WORK SUPPORTED BY:    THE ARM PROGRAM                  
C                                     OFFICE OF ENERGY RESEARCH        
C                                     DEPARTMENT OF ENERGY             
C                                                                      
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8           (V)
C
C     SUBROUTINE SOLIN inputs files for use with solar radiation
C     calculations for interpolation in SOLINT.  Reads files with
C     one record per panel.
C
      character*8      XID,       HMOLID,      YID       
      real*8               SECANT,       XALTZ 
C
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,
     *                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4
      COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF
C
      DIMENSION PNLHDR(2),XARRAY(*)
C
      EQUIVALENCE (PNLHDR(1),V1PBF)
C
      CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
      IF (KEOF.LE.0) RETURN
      CALL BUFIN (KFILE,KEOF,XARRAY(1),NLIMBF)
C
      V1P = V1PBF
      V2P = V2PBF
      DVP = DVPBF
      NLIM = NLIMBF
C
      RETURN
C
      END
C
C     ----------------------------------------------------------------
C
      SUBROUTINE SOLIN2 (V1P,V2P,DVP,NLIM,KFILE,XARAY1,XARAY2,KEOF)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      
C               LAST MODIFICATION:    1 November 1995
C                                                                      
C                  IMPLEMENTATION:    P.D. Brown
C                                                                      
C             ALGORITHM REVISIONS:    S.A. Clough
C                                     P.D. Brown
C                                                                      
C                                                                      
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.      
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139       
C                                                                      
C----------------------------------------------------------------------
C                                                                      
C               WORK SUPPORTED BY:    THE ARM PROGRAM                  
C                                     OFFICE OF ENERGY RESEARCH        
C                                     DEPARTMENT OF ENERGY             
C                                                                      
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8           (V)
C
C     SUBROUTINE SOLIN inputs files for use with solar radiation
C     calculations for interpolation in SOLINT.  Reads files with
C     two records per panel.
C
      character*8      XID,       HMOLID,      YID       
      real*8               SECANT,       XALTZ 
C
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,
     *                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4
      COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF
C
      DIMENSION PNLHDR(2),XARAY1(*),XARAY2(2)
C
      EQUIVALENCE (PNLHDR(1),V1PBF)
C
      CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
      IF (KEOF.LE.0) RETURN
      CALL BUFIN (KFILE,KEOF,XARAY1(1),NLIMBF)
      CALL BUFIN (KFILE,KEOF,XARAY2(1),NLIMBF)
C
      V1P = V1PBF
      V2P = V2PBF
      DVP = DVPBF
      NLIM = NLIMBF
C
      RETURN
C
      END
C
C     ----------------------------------------------------------------
C
      SUBROUTINE SOLOUT (V1P,V2P,DVP,NLIM,SOLRAD,LFILE,NPTS,NPANLS)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      
C               LAST MODIFICATION:    3 April 1994          
C                                                                      
C                  IMPLEMENTATION:    P.D. Brown
C                                                                      
C             ALGORITHM REVISIONS:    S.A. Clough
C                                     P.D. Brown
C                                                                      
C                                                                      
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.      
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139       
C                                                                      
C----------------------------------------------------------------------
C                                                                      
C               WORK SUPPORTED BY:    THE ARM PROGRAM                  
C                                     OFFICE OF ENERGY RESEARCH        
C                                     DEPARTMENT OF ENERGY             
C                                                                      
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8           (V)
C
C     SUBROUTINE SOLOUT OUTPUTS ATTENUATED SOLAR RADIANCE (INTERPOLATED)
C     TO LFILE
C
      COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4
      DIMENSION PNLHDR(2)
      DIMENSION SOLRAD(*)
C
      EQUIVALENCE (PNLHDR(1),V1PBF)
C
      REAL SOLRAD
C
      NPANLS = NPANLS+1
      V1PBF = V1P
      V2PBF = V2P
      DVPBF = DVP
      NLIMBF = NLIM
C
      CALL BUFOUT (LFILE,PNLHDR(1),NPHDRF)
      CALL BUFOUT (LFILE,SOLRAD(1),NLIMBF)
C
      IF (NPTS.GT.0) THEN
         IF (NPANLS.EQ.1) WRITE (IPR,900)
         WRITE (IPR,905)
         NNPTS = NPTS
         IF (NPTS.GT.(NLIMBF/2)+1) NNPTS = (NLIMBF/2)+1
         JEND = NLIMBF-NNPTS+1
         DO 10 J = 1, NNPTS
            VJ = V1PBF+FLOAT(J-1)*DVPBF
            K = J+JEND-1
            VK = V1PBF+FLOAT(K-1)*DVPBF
            WRITE (IPR,910) J,VJ,SOLRAD(J),K,VK,SOLRAD(K)
   10    CONTINUE
      ENDIF
C
      RETURN
C
  900 FORMAT ('0 ','LOCATION  WAVENUMBER',4X,'RADIANCE',13X,
     *        'LOCATION   WAVENUMBER',4X,
     *        'RADIANCE')
  905 FORMAT (' ')
  910 FORMAT (I8,2X,F12.6,1P,E15.7,0P,9X,I8,2X,F12.6,1P,E15.7,0P)
C
      END

