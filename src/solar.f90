!     revision:  $Revision$
!     created:   $Date$
!     presently: %H%  %T%
!
!  --------------------------------------------------------------------------
! |  Copyright �, Atmospheric and Environmental Research, Inc., 2015         |
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
SUBROUTINE SOLINT(IFILE,LFILE,NPTS,INFLAG,IOTFLG,JULDAT,ISOLVAR,  &
   SCON,SOLCYCFRAC,SOLVAR)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    24 May 2002
!
!                  IMPLEMENTATION:    P.D. Brown
!
!             ALGORITHM REVISIONS:    S.A. Clough
!                                     P.D. Brown
!                                     M.W. Shephard
!
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     131 Hartwell Ave, Lexington MA 02421
!
!----------------------------------------------------------------------
!
!               WORK SUPPORTED BY:    THE ARM PROGRAM
!                                     OFFICE OF ENERGY RESEARCH
!                                     DEPARTMENT OF ENERGY
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   USE phys_consts, ONLY: pi, radcn2
! mji - This array expanded to accommodate spectral extension to 86500 cm-1
!      PARAMETER (NSOL=2000001)
   PARAMETER (NSOL=2700001)
!
   IMPLICIT REAL*8           (V)
!
!     ------------------------------------------------------------
!     SUBROUTINE SOLINT interpolates solar radiances from the binary
!     file SOLAR.RAD.  The following are input and output options:
!
!       INFLAG = 0   => input transmittance from TAPE12 (default,
!                       where TAPE12 includes the monochromatic radiance
!                       and transmittance).
!              = 1   => input optical depths from TAPE12 and
!                       convert to transmittance.
!              = 2   => input R1, T1, T2, r1,  where
!
!                           _
!                           |
!                Observer |-O-|
!                           |
!                           -      \|/ Sun
!                          /|\    --O--
!                           |      /|\
!                           |   S2,/  /
!                           |   T2/  /
!                      R1,T1|    / |/_ R3
!                           |   /  /
!                           |  /  /
!                         r1||/_|/_r3
!            Ground    ---------------------------
!                     ///////////////////////////
!
!              = 3   => input transmittance from TAPE12 (CHARTS-type
!                       output)
!
!
!       IOTFLG = 0   => attenuate w/transmittance & output (default).
!              = 1   => attenuate and add to radiance from TAPE12
!                       (requires INFLAG = 1).
!              = 2   => Calculate solar contribution Rs = S2*T2*r1*T1
!                       and add to thermal contribution R1.
!
!
!        ISOLVAR = -1 => uses solar source file with a single component:
!                       no temporal variability; assumes Kurucz extraterrestrial
!                       solar irradiance, which yields a solar constant of 1368.22 Wm-2,
!                       unless scaled by SCON.
!
!        ISOLVAR =  0 => uses solar source file with a single component:
!                       no temporal variability; assumes NRLSSI2 extraterrestrial
!                       solar irradiance, which yields a solar constant of 1360.85 Wm-2,
!                       (for the spectral range 100-50000 cm-1 with quiet sun, facular
!                       and sunspot contributions fixed to the mean of
!                       Solar Cycles 13-24 and averaged over the mean solar cycle),
!                       unless scaled by SCON.

!        SCON   = 0.0 => no scaling of internal solar irradiance
!               > 0.0 => ISOLVAR = -1 or 0
!                        Total solar irradiance is scaled to SCON
!               > 0.0 => ISOLVAR = 1
!                        integral of total solar irradiance averaged over solar cycle
!                        is scaled to SCON
!
!        ISOLVAR =  1 => uses solar source file with multiple components;
!                       facular brightening and sunspot blocking amplitudes
!                       are by default determined by the fraction of the
!                       way into the solar cycle (see SOLCYCFRAC) or can
!                       be scaled (see SOLVAR)
!        ISOLVAR =  2 => uses solar source file with multiple components;
!                       facular brightening and sunspot blocking amplitudes
!                       are determined by the Mg and SB indeces (see SOLVAR)
!
!        SOLCYCFRAC   Solar cycle fraction (0-1); fraction of the way through the mean 11-year
!                     cycle with 0 and 1 defined as the minimum phase of the solar cycle
!                     (ISOLVAR=1 only)

!       SOLVAR        Solar variability scaling factors or indices (ISOLVAR=1,2 only)
!                     ISOLVAR = 1 =>
!                      SOLVAR(1)    Facular (Mg) index amplitude scale factor
!                      SOLVAR(2)    Sunspot (SB) index amplitude scale factor

!                      ISOLVAR = 2 =>
!                      SOLVAR(1)    Facular (Mg) index as defined in the NRLSSI2 model;
!                                   used for modeling specific solar activity
!                      SOLVAR(2)    Sunspot (SB) index as defined in the NRLSSI2 model;
!                                   used for modeling specific solar activity
!
!     Output radiance goes to TAPE13.
!
!
!     ------------------------------------------------------------
!
   COMMON /MANE/ P0,TEMP0,NLAYER,DDUM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   character*8      XIDS,       HMLIDS,       YIDS
   real*8                SECNTS,       XALTZS

!
!
   COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),      &
   &               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,  &
   &               EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /SLHDR/ XIDS(10),SECNTS,PAVES,TAVES,HMLIDS(60),XALTZS(4),  &
   &               WKS(60),PZLS,PZUS,TZLS,TZUS,WBRODS,DVS,V1S,V2S,    &
   &               TBONDS,                                            &
   &               EMSIVS,FSCDDS(17),NMOLS,LAYERS,YI1S,YIDS(10),      &
   &               LSTWDS
   COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
   COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /CVRSOL/ HNAMSOL,HVRSOL
!
!     ----------------------------------------------------------------
!     Parameter and common blocks for direct input of emission and
!     reflection function values
!
   PARAMETER (NMAXCO=4040)
   COMMON /EMSFIN/ V1EMIS,V2EMIS,DVEMIS,NLIMEM,ZEMIS(NMAXCO)
   COMMON /RFLTIN/ V1RFLT,V2RFLT,DVRFLT,NLIMRF,ZRFLT(NMAXCO)
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,pad_3,&
   &                angle_path,secant_diffuse,secant_path,diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!     ----------------------------------------------------------------
!
   DIMENSION XFILHD(2),XSOLHD(2),PNLHDR(2),OPNLHD(2)
   DIMENSION A1(0:100),A2(0:100),A3(0:100),A4(0:100)
   DIMENSION A1T2(0:100),A2T2(0:100),A3T2(0:100),A4T2(0:100)
   DIMENSION A1RF(0:100),A2RF(0:100),A3RF(0:100),A4RF(0:100)
   DIMENSION TRAO(2),TRAN(2410)
   DIMENSION RADO(2),RADN(2410)
   DIMENSION OPTO(2),OPTN(2410)
!
   real          solvar(2),svar(3)

   DIMENSION SOLAR(-1:NSOL)
   DIMENSION SOLAR3(2400,3)
   DIMENSION TRAN2(-1:NSOL)
   DIMENSION RAD2(-1:NSOL)
   DIMENSION XRFLT(-1:NSOL)
   DIMENSION SOLRAD(2410)
!
   CHARACTER*40 CYID
   CHARACTER*18 HNAMSOL,HVRSOL
!
   EQUIVALENCE (XSOLHD(1),XIDS(1))
   EQUIVALENCE (XFILHD(1),XID(1)),(PNLHDR(1),V1P),                   &
   &            (OPNLHD(1),V1PO)
   EQUIVALENCE (TRAN(1),TRAO(1)),(RADN(1),RADO(1)),                  &
   &            (OPTN(1),OPTO(1)),                                    &
   &            (FSCDID(4),IAERSL),(FSCDID(5),IEMIT),                 &
   &            (FSCDID(6),ISCHDR),(FSCDID(7),IPLOT),                 &
   &            (FSCDID(8),IPATHL),(FSCDID(12),XSCID),                &
   &            (FSCDID(16),LAYR1)
   equivalence (dv_sol,xhdr_s(218)),                                 &
   &            (v1_sol,xhdr_s(219)),(v2_sol,xhdr_s(221))
!
   real*4 xhdr_s(265),dv_sol
!
!     ************************************************************
!     ****** THIS PROGRAM DOES MERGE FOR SOLAR RADIANCE AND ******
!     ****** TRANMITTANCE USING FOUR POINT INTERPOLATION    ******
!     ************************************************************
!
!
!     ASSIGN CVS VERSION NUMBER TO MODULE
!
   HVRSOL = '$Revision$'
!
!     -------------------
!
!     Open file SOLAR.RAD
!
   ISOLFL = 19
   OPEN(UNIT=ISOLFL,FILE='SOLAR.RAD',FORM='UNFORMATTED',             &
   &     STATUS='OLD')
!
!     Note that the file SOLAR.RAD is always single precision. Provision
!     has been made to deal with the case in which the current program i
!     double precision.
!
   CALL CPUTIM (TIME)
   WRITE (IPR,900) TIME
   NPANLS = 0
   TIMRD = 0.0
   TIMOT = 0.0

!    Calculate solar scaling factors

   call scale_solar (isolvar,scon,solcycfrac,solvar,svar)
!
!     Calculate Earth distance to sun given Julian Day JULDAT. Used to
!     scale solar source function. Formula taken from "Atmospheric Radia
!     Transfer", J. Lenoble, 1993.

!     Test validity of JULDAT

   if ((juldat.lt.0).or.(juldat.gt.366)) then
      write(*,*) 'JULDAT = ',juldat,' is out of range 0-366.'
      write(ipr,*) 'JULDAT = ',juldat,' is out of range 0-366.'
      stop 'Stopped in SOLINT'
   endif

!     If JULDAT = 0 , then set XJUL_SCALE to 1

   if (juldat.eq.0) then
      XJUL_SCALE = 1.0
      write(ipr,*) 'JULDAT = 0, no scaling of solar source function'
   else
      theta = 2*pi*( REAL(JULDAT)-1.)/365.
      XJUL_SCALE = 1.00011 + 0.034221*cos(theta) +                   &
      &        1.28E-3*sin(theta) + 7.19E-4*cos(2.*theta) +              &
      &        7.7E-5*sin(2.*theta)
      write(ipr,*) 'JULDAT = ',JULDAT,                               &
      &        ', scale factor for solar source function = ',XJUL_SCALE
   endif

!     Combine Julian day scaling with solar source scaling when there is no solar variability
   if (isolvar.le.0) then
      xjul_scale = xjul_scale*solvar(1)
   endif

!     FOR AEROSOL RUNS, MOVE YID (IFILE) INTO YID (LFILE)
!
!     Read file header of solar radiance file and determine dv ratio
!
   IF (IAERSL.GT.0) WRITE (CYID,'(5A8)') (YID(I),I=3,7)
!
   read (isolfl) (xhdr_s(i),i=1,264)
!
   dv = dv_sol
   v1 = v1_sol
   v2 = v2_sol
!
!pdb      IF (IAERSL.GT.0) READ (CYID,'(5A8)') (YID(I),I=3,7)
!pdb      IF (JPATHL.GE.1) IPATHL = JPATHL
   DVK = DV
!
!     Read in file header of transmittance/optical depth file
!
   CALL BUFIN (IFILE,LEOF,XFILHD(1),NFHDRF)
   DVL = DV
   DVMIN = DVL
!
   ATYPE = 9.999E09
   IF (DVK.EQ.DVL) ATYPE = 0.
   IF (DVL.GT.DVK) ATYPE = DVK/(DVL-DVK)+0.5
   IF (DVL.LT.DVK) ATYPE = -DVL/(DVK-DVL)-0.5
!
!     IF (ATYPE .GT. 0) STOP  ' SOLINT; ATYPE GT 0 '
!
!
!     Write file information out to TAPE6
!     If INFLAG = 2, then also open files and read headers
!     for downward transmittance from solar path and solar
!     reflectance at the ground.
!
   WRITE (IPR,905) ISOLFL,IFILE,LFILE,ATYPE,INFLAG,IOTFLG
   IF (INFLAG.EQ.0) THEN
      WRITE(IPR,920) IFILE
   ELSEIF (INFLAG.EQ.1) THEN
      WRITE(IPR,925) IFILE
   ELSEIF (INFLAG.EQ.2) THEN
      ISLTRN = 20
      OPEN(UNIT=ISLTRN,FILE='SOL.PATH.T2',FORM='UNFORMATTED',        &
         STATUS='OLD')
      CALL BUFIN (ISLTRN,ITEOF,XSOLHD(1),NFHDRF)
      DVT2 = DVS
      ISLRFL = 21
      OPEN(UNIT=ISLRFL,FILE='SOL.REFLECTANCE',FORM='UNFORMATTED',    &
         STATUS='OLD')
      CALL BUFIN (ISLRFL,IREOF,XSOLHD(1),NFHDRF)
      DVRF = DVS
!
!        **************************************************************
!
!        Read Record 1.4 from TAPE5
!
      READ (IRD,970,END=80) TMPBND,(BNDEMI(IBND),IBND=1,3), (BNDRFL( &
         IBND),IBND=1,3), surf_refl
!
      BNDTST = ABS(BNDRFL(1))+ABS(BNDRFL(2))+ABS(BNDRFL(3))
      IF (BNDTST.NE.0.) IBPROP = 1
!
!        If BNDEMI(1) < 0, read in coefficients from file 'EMISSIVITY'
!        If BNDEMI(1) > 0, check to see if emissivity is reasonable
!
!        UNIT ICOEF used for input files
!
      ICOEF = 14
!
      IF (BNDEMI(1).LT.0) THEN
         OPEN (UNIT=ICOEF,FILE='EMISSIVITY',STATUS='OLD')
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
!
!        ------------------------------------------------------------
!        *** NOTE: REFLECTION FUNCTION NOT CURRENTLY INCORPORATED ***
!        *** INTO SOLAR RADIATIVE TRANSFER CALCULATION            ***
!
!        If BNDRFL(1) < 0, read in coefficients from file 'REFLECTIVITY'
!        If BNDRFL(1) > 0, check to see if reflectivity is reasonable
!
      IF (BNDRFL(1).LT.0) THEN
         OPEN (UNIT=ICOEF,FILE='REFLECTIVITY',STATUS='OLD')
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
!        ------------------------------------------------------------
!
!        TBOUND is the boundary temperature. TBOUND=0. for no boundary
!        EMISIV is the boundary emissivity
!        Set default for EMISIV
!
      EMITST = ABS(BNDEMI(1))+ABS(BNDEMI(2))+ABS(BNDEMI(3))
      IF ((TMPBND.GT.0.).AND.(EMITST.EQ.0.)) BNDEMI(1) = 1.
      EMISIV = BNDEMI(1)
      TBOUND = TMPBND
      XKTBND = TBOUND/RADCN2
      WRITE (IPR,985) V1,V2,TBOUND,(BNDEMI(IBND),IBND=1,3), (BNDRFL( &
         IBND),IBND=1,3), surf_refl
!
!     **************************************************************
!
!        Determine the minimum and maximum DV of all files and reset
!        ATYPE (-1. is only a flag for nonzero ATYPE)
!
      DVMIN = MIN(DVL,DVK,DVT2,DVRF,DVEMIS)
      DVMAX = MAX(DVL,DVK,DVT2,DVRF,DVEMIS)
!         IF (DVMAX.EQ.DVMIN) THEN
!            ATYPE = 0.0
!         ELSE
      ATYPE = -1.
!         ENDIF
      WRITE(IPR,927) IFILE,ISLTRN,ISLRFL
      WRITE(IPR,941) DVT2,DVRF
!
   ELSEIF (INFLAG.EQ.3) THEN
      WRITE(IPR,926) IFILE
   ENDIF
   WRITE(IPR,906) DVL,DVK
!
!     Test for INFLAG and IOTFLG compatibility
!
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
!
!     Set XSCID and ISCHDR to values which allow for potential
!     interpolation of output file
!
   XSCID = 100.01
   ISCHDR = 10
!
!     Output file header
!
   CALL BUFOUT (LFILE,XFILHD(1),NFHDRF)
!
   IF (ATYPE.EQ.0.) THEN
!
!     ======================================
!        1/1 ratio only
!     ======================================
!
      WRITE(IPR,950) DVMIN
      IPANEM = 0
30    CONTINUE
      IPANEM = IPANEM+1
      CALL CPUTIM (TIMSL1)
      if (isolvar.le.0) then
         CALL SOLIN_sgl (V1P,V2P,DVP,NLIM,ISOLFL,SOLAR(1),LSEOF)
      else
         CALL SOLIN_tri (V1P,V2P,DVP,NLIM,ISOLFL,SOLAR3,LSEOF)
         call comb_solar (solar3,solar(1),nlim,isolvar,svar)
      end if

      CALL CPUTIM (TIMSL2)
      TIMRD = TIMRD+TIMSL2-TIMSL1
      IF (LSEOF.LE.0) GO TO 110
!
!
!        If INFLAG = 0, then read radiance and transmittance
!        If INFLAG = 1, then read optical depth
!        If INFLAG = 2, then read radiance and transmittance
!                       and call SOLIN to read in r1 and T2
!        If INFLAG = 3, then read transmittance

!
      IF (INFLAG.EQ.0) THEN
         CALL SOLIN2 (V1PO,V2PO,DVPO,NLIMO,IFILE,RADO(1), TRAO(1),   &
            LEOF)
      ELSEIF (INFLAG.EQ.1) THEN
         CALL SOLIN (V1PO,V2PO,DVPO,NLIMO,IFILE,OPTO(1), LEOF)
         DO 35 I = 1,NLIMO
            TRAN(I) = EXP(-OPTN(I))
35       CONTINUE
      ELSEIF (INFLAG.EQ.2) THEN
         CALL SOLIN2 (V1PO,V2PO,DVPO,NLIMO,IFILE,RADO(1), TRAO(1),   &
            LEOF)
         CALL SOLIN2 (V1T2,V2T2,DVT2,NLIMT2,ISLTRN,RAD2(1), TRAN2(1),&
            LSEOF)
         CALL SOLIN (V1RF,V2RF,DVRF,NLIMRF,ISLRFL,XRFLT(1), LSEOF)
         IF (ABS(V1T2-V1RF).GT.0.001) THEN
            WRITE(IPR,*) 'SOLINT: PANELS DO NOT MATCH:'
            WRITE(IPR,*) '  V1T2 = ',V1T2,'  V1RF = ',V1RF
            STOP 'SOLINT: PANELS DO NOT MATCH: SEE TAPE6'
         ENDIF
      ELSEIF (INFLAG.EQ.3) THEN
         CALL SOLIN (V1PO,V2PO,DVPO,NLIMO,IFILE,TRAO(1), LEOF)
      ENDIF
      CALL CPUTIM (TIMSL3)
      TIMRD = TIMRD+TIMSL3-TIMSL2
!
!        If IOTFLG = 0, then calculate attenuated solar radiance
!        If IOTFLG = 1, then calculate attenuated solar radiance
!                       plus atmospheric radiance
!        If IOTFLG = 2, then calculate attenuated solar radiance
!                       through the reflected atmosphere plus
!                       atmospheric radiance
!
!        Solar irradiance is input from SOLAR.RAD (W/m2 cm-1).
!        Convert to radiance units (W/cm2 sr cm-1) by multiplying
!        by 1/6.8e-5.

      conv_ster = 1./(1.e4*6.8e-5)

!
!        Combine XJUL_SCALE and conv_ster into one scale factor SCAL_FAC

      SCAL_FAC = conv_ster*XJUL_SCALE
!
      IF (IOTFLG.EQ.0) THEN
         DO 40 I = 1, NLIM
            SOLRAD(I) = SOLAR(I)*SCAL_FAC*TRAN(I)
40       CONTINUE
      ELSEIF (IOTFLG.EQ.1) THEN
         DO 41 I = 1, NLIM
            SOLRAD(I) = SOLAR(I)*SCAL_FAC*TRAN(I)+RADN(I)
41       CONTINUE
      ELSEIF (IOTFLG.EQ.2) THEN
         IF (TBOUND.EQ.0) THEN
            DO 42 I = 1, NLIM
               SOLRAD(I) = SOLAR(I)*SCAL_FAC*TRAN2(I)*XRFLT(I)*      &
                  TRAN(I)+RADN(I)
42          CONTINUE
         ELSE
            DO 43 I = 1, NLIM
               EMDUM = -1.
               BBDUM = -1.
               VBND = V1PO+(I-1)*DVPO
               ZEMSV = EMISFN(VBND,DVPO,VINEM,EMDEL,EMDUM)
               BBND = PLANCK(VBND, XKTBND)
               SOLRAD(I) = (SOLAR(I)*SCAL_FAC*TRAN2(I)*XRFLT(I)      &
                  +ZEMSV*BBND)*TRAN(I)+RADN(I)
43          CONTINUE
         ENDIF
      ENDIF
!
      CALL CPUTIM (TIMSL2)
      CALL SOLOUT(V1PO,V2PO,DVPO,NLIMO,SOLRAD,LFILE,NPTS,NPANLS)
      CALL CPUTIM (TIMSL3)
      TIMOT = TIMOT+TIMSL3-TIMSL2
      GO TO 30
!
   ENDIF
!
!     ======================================
!     All ratios except 1/1
!     ======================================
!
   WRITE(IPR,951) DVMIN
   DO 50 JP = 0,100
      APG = JP
      P = 0.01*APG
!
!        The following are the constants for the Lagrange
!        4 point interpolation
!
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
50 END DO
!
!     *** Beginning of loop that does merge  ***
!
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
!
!     ============================================================
!
   IPANEM = 0
60 CONTINUE
   IPANEM = IPANEM+1
!
!     If INFLAG = 0, then read radiance and transmittance
!     If INFLAG = 1, then read optical depth
!     If INFLAG = 2, then read radiance and transmittance
!                    and call SOLIN to read in r1 and T2
!     If INFLAG = 3, then read transmittance
!
   IF (INFLAG.EQ.0) THEN
      CALL SOLIN2 (V1PO,V2PO,DVPO,NLIMO,IFILE,RADO(1), TRAO(1),LEOF)
      IF (LEOF.LE.0) GO TO 110
   ELSEIF (INFLAG.EQ.1) THEN
      CALL SOLIN (V1PO,V2PO,DVPO,NLIMO,IFILE,OPTO(1), LEOF)
      IF (LEOF.LE.0) GO TO 110
      DO 65 I = 1,NLIMO
         TRAN(I) = EXP(-OPTN(I))
65    CONTINUE
   ELSEIF (INFLAG.EQ.2) THEN
      CALL SOLIN2 (V1PO,V2PO,DVPO,NLIMO,IFILE,RADO(1), TRAO(1),LEOF)
      IF (LEOF.LE.0) GO TO 110
!
!        -----------------------------------------------------
!        TRAN2 and RAD2 read in
!
      IF (V2T2.LE.V2PO+DVT2.AND.LT2EOF.GT.0) THEN
66       CALL CPUTIM(TIMSL2)
         CALL SOLIN2 (V1T2,V2T2,DVT2,NLIMT2,ISLTRN,RAD2(npt2+1),     &
            TRAN2(npt2+1),LT2EOF)
         CALL CPUTIM(TIMSL3)
         TIMRD = TIMRD+TIMSL3-TIMSL2
         IF (LT2EOF.LE.0) GO TO 67
         V1T2 = V1T2- REAL(NPT2)*DVT2
1        NPT2 = NLIMT2+NPT2
         IF (V2T2.LE.V2PO+DVT2) GO TO 66
      ENDIF
!
!        Zero point of first panel
!
67    IF (TRAN2(0).EQ.0.0) THEN
         TRAN2(-1) = TRAN2(1)
         TRAN2(0) = TRAN2(1)
      ENDIF
!
!        End point of last panel
!
      IF (V2T2+DVT2.GE.V2) THEN
         TRAN2(NPT2+1) = TRAN2(NPT2)
         TRAN2(NPT2+2) = TRAN2(NPT2)
      ENDIF
!
!        -----------------------------------------------------
!        XRFLT read in
!
      IF (V2RF.LE.V2PO+DVRF.AND.LRFEOF.GT.0) THEN
68       CALL CPUTIM(TIMSL2)
         CALL SOLIN (V1RF,V2RF,DVRF,NLIMRF,ISLRFL,XRFLT(nprf+1),     &
            LRFEOF)
         CALL CPUTIM(TIMSL3)
         TIMRD = TIMRD+TIMSL3-TIMSL2
         IF (LRFEOF.LE.0) GO TO 69
         V1RF = V1RF- REAL(NPRF)*DVRF
         NPRF = NLIMRF+NPRF
         IF (V2RF.LE.V2PO+DVRF) GO TO 68
      ENDIF
!
!        Zero point of first panel
!
69    IF (XRFLT(0).EQ.0.0) THEN
         XRFLT(-1) = XRFLT(1)
         XRFLT(0) = XRFLT(1)
      ENDIF
!
!        End point of last panel
!
      IF (V2RF+DVRF.GE.V2) THEN
         XRFLT(NPRF+1) = XRFLT(NPRF)
         XRFLT(NPRF+2) = XRFLT(NPRF)
      ENDIF
   ELSEIF (INFLAG.EQ.3) THEN
      CALL SOLIN (V1PO,V2PO,DVPO,NLIMO,IFILE,TRAO(1), LEOF)
      IF (LEOF.LE.0) GO TO 110
   ENDIF
!     -----------------------------------------------------
!
   CALL CPUTIM(TIMSL3)
   TIMRD = TIMRD+TIMSL3-TIMSL2
   II = 1
!
!     Buffer in panels from solar radiance file
!
   IF (V2P.LE.V2PO+DVP .AND.LSEOF.GT.0) THEN
70    CALL CPUTIM(TIMSL2)
      if (isolvar.le.0) then
         CALL SOLIN_sgl (V1P,V2P,DVP,NLIM,ISOLFL,SOLAR(NPE+1),LSEOF)
      else
         CALL SOLIN_tri (V1P,V2P,DVP,NLIM,ISOLFL,SOLAR3,LSEOF)
         call comb_solar (solar3,solar(npe+1),nlim,isolvar,svar)
      endif

      CALL CPUTIM(TIMSL3)
      TIMRD = TIMRD+TIMSL3-TIMSL2
      IF (LSEOF.LE.0) GO TO 80
      V1P = V1P- REAL(NPE)*DVP
      NPE = NLIM+NPE
      IF (V2P.LE.V2PO+DVP) GO TO 70
   ENDIF
!
!     Zero point of first panel
!
80 IF (SOLAR(0).EQ.0.0) THEN
      SOLAR(-1) = SOLAR(1)
      SOLAR(0) = SOLAR(1)
   ENDIF
!
!     End point of last panel
!
   IF (V2P+DVP.GE.V2) THEN
      SOLAR(NPE+1) = SOLAR(NPE)
      SOLAR(NPE+2) = SOLAR(NPE)
   ENDIF
!
!     -----------------------------------------------------
!
!     NPL is the location of first element on arrays RADO and TRAO
!
   NPL = 1
   NPLT2 = 1
   NPLRF = 1
!
!     Set spectral spacing ratios (uses DVMIN instead of DVL, for
!     coding where it is not assumed that RADN & TRAN have the
!     smallest spacing).
!
   RATDVR = DVMIN/DVK
   IF (INFLAG.EQ.2) THEN
      RTDVT2 = DVMIN/DVT2
      RTDVRF = DVMIN/DVRF
   ENDIF
!
!     Test for ratios of 1 (no need for interpolation).  Reset
!     interpolation coefficients if necessary.
!
   IF (RATDVR.EQ.1) THEN
      DO 82 JP = 0,100
         A1(JP) = 0.0
         A2(JP) = 1.0
         A3(JP) = 0.0
         A4(JP) = 0.0
82    CONTINUE
   ENDIF
   IF (INFLAG.EQ.2) THEN
      IF (RTDVT2.EQ.1) THEN
         DO 84 JP = 0,100
            A1T2(JP) = 0.0
            A2T2(JP) = 1.0
            A3T2(JP) = 0.0
            A4T2(JP) = 0.0
84       CONTINUE
      ENDIF
      IF (RTDVRF.EQ.1) THEN
         DO 86 JP = 0,100
            A1RF(JP) = 0.0
            A2RF(JP) = 1.0
            A3RF(JP) = 0.0
            A4RF(JP) = 0.0
86       CONTINUE
      ENDIF
   ENDIF
!
!     -----------------------------------------------------
!
!     FJJ is offset by 2. (for rounding purposes)
!
   FJ1DIF = (V1PO-V1P)/DVP+1.+2.
   IF (INFLAG.EQ.2) THEN
      FJDFT2 = (V1PO-V1T2)/DVT2+1.+2.
      FJDFRF = (V1PO-V1RF)/DVRF+1.+2.
   ENDIF
!
!     ============================================================
!
!     ***** Beginning of loop that does merge  *****
!
!     If IOTFLG = 0, then calculate attenuated solar radiance
!     If IOTFLG = 1, then calculate attenuated solar radiance
!                    plus atmospheric radiance
!     If IOTFLG = 2, then calculate attenuated solar radiance
!                    through the reflected atmosphere plus
!                    atmospheric radiance
!
!     Solar irradiance is input from SOLAR.RAD (W/m2 cm-1).
!     Convert to radiance units (W/cm2 sr cm-1) by multiplying
!     by 1/6.8e-5.

   conv_ster = 1./(1.e4*6.8e-5)

!
!     Combine XJUL_SCALE and conv_ster into one scale factor SCAL_FAC

   SCAL_FAC = conv_ster*XJUL_SCALE
!

   IF (IOTFLG.EQ.0) THEN
      DO 90 II = 1, NLIMO
         FJJ = FJ1DIF+RATDVR* REAL(II-1)
         JJ = INT(FJJ)-2
         JP = (FJJ- REAL(JJ))*100.-199.5
         SOLRAD(II) = (A1(JP)*SOLAR(JJ-1)+A2(JP)*SOLAR(JJ)+ A3(JP)*  &
            SOLAR(JJ+1)+A4(JP)*SOLAR(JJ+2))*SCAL_FAC* TRAN(II)
!
90    CONTINUE
   ELSEIF (IOTFLG.EQ.1) THEN
      DO 91 II = 1, NLIMO
         FJJ = FJ1DIF+RATDVR* REAL(II-1)
         JJ = INT(FJJ)-2
         JP = (FJJ- REAL(JJ))*100.-199.5
         SOLRAD(II) = (A1(JP)*SOLAR(JJ-1)+A2(JP)*SOLAR(JJ)+ A3(JP)*  &
            SOLAR(JJ+1)+A4(JP)*SOLAR(JJ+2))*SCAL_FAC* TRAN(II)+RADN(II)
91    CONTINUE
   ELSEIF (IOTFLG.EQ.2) THEN
      IF (TBOUND.EQ.0.) THEN
         DO 92 II = 1, NLIMO
            FJJ = FJ1DIF+RATDVR* REAL(II-1)
            FJJT2 = FJDFT2+RTDVT2* REAL(II-1)
            FJJRF = FJDFRF+RTDVRF* REAL(II-1)
            JJ = INT(FJJ)-2
            JJT2 = INT(FJJT2)-2
            JJRF = INT(FJJRF)-2
            JP = (FJJ- REAL(JJ))*100.-199.5
            JPT2 = (FJJT2- REAL(JJT2))*100.-199.5
            JPRF = (FJJRF- REAL(JJRF))*100.-199.5
            ZSOL = (A1(JP)*SOLAR(JJ-1)+A2(JP)*SOLAR(JJ)+ A3(JP)*     &
               SOLAR(JJ+1)+A4(JP)*SOLAR(JJ+2))*SCAL_FAC
            ZTR2 = (A1T2(JPT2)*TRAN2(JJT2-1)+ A2T2(JPT2)*TRAN2(JJT2)+&
               A3T2(JPT2)*TRAN2(JJT2+1)+ A4T2(JPT2)*TRAN2(JJT2+2))
            ZRFL = (A1RF(JPRF)*XRFLT(JJRF-1)+ A2RF(JPRF)*XRFLT(JJRF)+&
               A3RF(JPRF)*XRFLT(JJRF+1)+ A4RF(JPRF)*XRFLT(JJRF+2))
            SOLRAD(II) = ZSOL*ZTR2*ZRFL*TRAN(II)+RADN(II)
92       CONTINUE
      ELSE
         DO 93 II = 1, NLIMO
            FJJ = FJ1DIF+RATDVR* REAL(II-1)
            FJJT2 = FJDFT2+RTDVT2* REAL(II-1)
            FJJRF = FJDFRF+RTDVRF* REAL(II-1)
            JJ = INT(FJJ)-2
            JJT2 = INT(FJJT2)-2
            JJRF = INT(FJJRF)-2
            JP = (FJJ- REAL(JJ))*100.-199.5
            JPT2 = (FJJT2- REAL(JJT2))*100.-199.5
            JPRF = (FJJRF- REAL(JJRF))*100.-199.5
            ZSOL = (A1(JP)*SOLAR(JJ-1)+A2(JP)*SOLAR(JJ)+ A3(JP)*     &
               SOLAR(JJ+1)+A4(JP)*SOLAR(JJ+2))*SCAL_FAC
            ZTR2 = (A1T2(JPT2)*TRAN2(JJT2-1)+ A2T2(JPT2)*TRAN2(JJT2)+&
               A3T2(JPT2)*TRAN2(JJT2+1)+ A4T2(JPT2)*TRAN2(JJT2+2))
            ZRFL = (A1RF(JPRF)*XRFLT(JJRF-1)+ A2RF(JPRF)*XRFLT(JJRF)+&
               A3RF(JPRF)*XRFLT(JJRF+1)+ A4RF(JPRF)*XRFLT(JJRF+2))
            EMDUM = -1.
            VBND = V1PO+(II-1)*DVPO
            ZEMSV = EMISFN(VBND,DVPO,VINEM,EMDEL,EMDUM)
            BBND = PLANCK(VBND, XKTBND)
            SOLRAD(II) = (ZSOL*ZTR2*ZRFL+ZEMSV*BBND)* TRAN(II)+RADN( &
               II)
93       CONTINUE
      ENDIF
   ENDIF
!
   NPL = JJ-1
   NPLT2 = JJT2-1
   NPLRF = JJRF-1
!
   CALL CPUTIM (TIMSL1)
!
!     ============================================================
!
!     Output attenuated radiance
!
   CALL SOLOUT(V1PO,V2PO,DVPO,NLIMO,SOLRAD,LFILE,NPTS,NPANLS)
   CALL CPUTIM (TIMSL2)
   TIMOT = TIMOT+TIMSL2-TIMSL1
!
!     ============================================================
!
!     Reset element locations for each array
!
!     ============================================================
!     NPL is now location of first element in the array SOLAR to
!     be used for next pass.
!
   IPL = -2
   DO 100 NL = NPL, NPE
      IPL = IPL+1
      SOLAR(IPL) = SOLAR(NL)
100 END DO
!
   V1P = V1P+ REAL(NPL+1)*DVP
   NPE = IPL
!
   IF (IOTFLG.EQ.2) THEN
!
!        ------------------------------------------------------------
!        NPLT2 is now location of first element in the array TRAN2 to
!        be used for next pass.
!
      IPL = -2
      DO 102 NL = NPLT2, NPT2
         IPL = IPL+1
         TRAN2(IPL) = TRAN2(NL)
102   CONTINUE
!
      V1T2 = V1T2+ REAL(NPLT2+1)*DVT2
      NPT2 = IPL
!
!        -------------------------------------------------------------
!        NPLRF is now location of first element in the array XRFLT to
!        be used for next pass.
!
      IPL = -2
      DO 104 NL = NPLRF, NPRF
         IPL = IPL+1
         XRFLT(IPL) = XRFLT(NL)
104   CONTINUE
!
      V1RF = V1RF+ REAL(NPLRF+1)*DVRF
      NPRF = IPL
   ENDIF
!     ============================================================
!
   GO TO 60
110 CONTINUE
!
   CALL CPUTIM (TIME1)
   TIM = TIME1-TIME
   WRITE (IPR,910) TIME1,TIM,TIMRD,TIMOT
!
   RETURN
!
900 FORMAT ('0 THE TIME AT THE START OF SOLINT IS ',F12.3)
905 FORMAT ('0 FILE ',I5,' MERGED WITH FILE ',I5,' ONTO FILE',        &
   &        I5,'  WITH XTYPE=',G15.5,/,'0 INFLAG = ',I5,4X,           &
   &        'IOTFLG = ',I5)
906 FORMAT ('0          Thermal spectrum spacing = ',E10.5,/,         &
   &        '0   Solar radiance spectral spacing = ',E10.5,/)
910 FORMAT ('0 THE TIME AT THE END OF SOLINT IS ',F12.3/F12.3,        &
   &        ' SECS WERE REQUIRED FOR THIS SOLAR MERGE',F12.3,         &
   &        ' - READ - ',F12.3,' - SOLOUT - ',F12.3)
920 FORMAT ('0 Radiance and Transmittance read in from unit',I5)
925 FORMAT ('0 Optical Depths read in from unit',I5)
926 FORMAT ('0 Transmittance read in from unit',I5)
927 FORMAT ('0 Thermal upwelling Radiance & Transmittance',           &
   &        ' read in from unit',I5,/,                                &
   &        '  Solar reflectance function read in from unit',         &
   &        I5,/,                                                     &
   &        '  Thermal downwelling Radiance & Transmittance',         &
   &        ' read in from unit',I5,/)
930 FORMAT ('0 Attenuated solar radiance output to unit',I5,/)
935 FORMAT ('0 Attenuated solar radiance + atmospheric radiance',     &
   &        1x,'output to unit',I5,/)
940 FORMAT ('0 Attenuated solar radiance + atmospheric radiance',     &
   &        1x,'(including effects of reflection) output to unit',    &
   &        I5,/)
941 FORMAT ('0       Thermal downwelling spacing = ',E10.5,/,         &
   &        '0       Reflection function spacing = ',E10.5)
950 FORMAT ('0 No interpolation needed: using DV = ',E10.5,//)
951 FORMAT ('0 Interpolating to spectral spacing = ',E10.5,//)
970 FORMAT (7E10.3,4X,A1)
975 FORMAT ('0 FOR VNU = ',F10.3,' THE EMISSIVITY = ',E10.3,          &
   &        ' AND IS NOT BOUNDED BY (0.,1.) ')
980 FORMAT ('0 FOR VNU = ',F10.3,' THE REFLECTIVITY = ',E10.3,        &
   &        ' AND IS NOT BOUNDED BY (0.,1.) ')
985 FORMAT (5(/),'0*********** BOUNDARY PROPERTIES ***********',/,    &
   &        '0 V1(CM-1) = ',F12.4,/,'0 V2(CM-1) = ',F12.4,/,          &
   &        '0 TBOUND   = ',F12.4,5X,'BOUNDARY EMISSIVITY   = ',      &
   &        3(1PE11.3),/,'0',29X,'BOUNDARY REFLECTIVITY = ',          &
   &        3(1PE11.3), ' SURFACE REFLECTIVITY = ',4X,1A)

!
end subroutine SOLINT
!
!     ----------------------------------------------------------------
!
SUBROUTINE SOLIN (V1P,V2P,DVP,NLIM,KFILE,XARRAY,KEOF)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    1 November 1995
!
!                  IMPLEMENTATION:    P.D. Brown
!
!             ALGORITHM REVISIONS:    S.A. Clough
!                                     P.D. Brown
!
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!----------------------------------------------------------------------
!
!               WORK SUPPORTED BY:    THE ARM PROGRAM
!                                     OFFICE OF ENERGY RESEARCH
!                                     DEPARTMENT OF ENERGY
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE SOLIN inputs files for use with solar radiation
!     calculations for interpolation in SOLINT.  Reads files with
!     one record per panel.
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /BUFPNL_s/ V1PBF,V2PBF,DVPBF,NLIMBF
!
   DIMENSION PNLHDR(2),XARRAY(*)
!
   EQUIVALENCE (PNLHDR(1),V1PBF)
!
   !real*4 dvpbf,pnlhdr
   !integer*4 nlimbf
!
   CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
   IF (KEOF.LE.0) RETURN
   CALL BUFIN (KFILE,KEOF,XARRAY(1),NLIMBF)
!
   V1P = V1PBF
   V2P = V2PBF
   DVP = DVPBF
   NLIM = NLIMBF
!
   RETURN
!
end subroutine SOLIN
!
!     ----------------------------------------------------------------
!
SUBROUTINE SOLIN_sgl (V1P,V2P,DVP,NLIM,KFILE,XARRAY,KEOF)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    1 November 1995
!
!                  IMPLEMENTATION:    P.D. Brown
!
!             ALGORITHM REVISIONS:    S.A. Clough
!                                     P.D. Brown
!
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!----------------------------------------------------------------------
!
!               WORK SUPPORTED BY:    THE ARM PROGRAM
!                                     OFFICE OF ENERGY RESEARCH
!                                     DEPARTMENT OF ENERGY
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE SOLIN_sgl inputs files for use with solar radiation
!     calculations for interpolation in SOLINT.  Reads files with
!     one record per panel. SOLIN_sgl reads single precision files
!
   COMMON /BUFPNL_s/ V1PBF,V2PBF,dvpbf,nlimbf
!
   DIMENSION PNLHDR(2),XARRAY(*),xarray_s(2410)
!
   EQUIVALENCE (PNLHDR(1),V1PBF)
!
   real*4 dvpbf,pnlhdr,xarray_s

   integer*4 kfil_s,keof_s,nphdr_s,nphdrf,nlimbf

   data nphdrf / 6 /

   kfil_s = KFILE
   keof_s = KEOF

   CALL BUFIN_sgl (kfil_s,keof_s,pnlhdr(1),nphdrf)

   KEOF = keof_s
   IF (KEOF.LE.0) RETURN

   CALL BUFIN_sgl (kfil_s,keof_s,xarray_s(1),nlimbf)
!
   KEOF = keof_s

   V1P = V1PBF
   V2P = V2PBF
   DVP = dvpbf
   NLIM = nlimbf
!
!     The variable XARRAY (either single or double, depending on the
!     complile option specified, is set equal to the real*4 variable xar

   do 10 i=1,nlimbf
      XARRAY(i) = xarray_s(i)
10 continue
!
   RETURN
!
end subroutine SOLIN_sgl
!
!     ----------------------------------------------------------------
SUBROUTINE SOLIN_tri (V1P,V2P,DVP,NLIM,KFILE,XARRAY,KEOF)
!
!
!                  Written  :         January 2017
!
!                  IMPLEMENTATION:    K. E. Cady-Pereira

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE SOLIN_tri inputs files for use with solar radiation
!     calculations for interpolation in SOLINT.  Reads files with
!     three records per panel. SOLIN_tri reads single precision files
!
   COMMON /BUFPNL_s/ V1PBF,V2PBF,dvpbf,nlimbf
!
   DIMENSION PNLHDR(2),XARRAY(2400,3)
!
   EQUIVALENCE (PNLHDR(1),V1PBF)
!
   real*4 dvpbf,pnlhdr,xarray_s(2410)

   integer*4 kfil_s,keof_s,nphdr_s,nphdrf,nlimbf

   data nphdrf / 6 /

   kfil_s = KFILE
   keof_s = KEOF

   CALL BUFIN_sgl (kfil_s,keof_s,pnlhdr(1),nphdrf)

   KEOF = keof_s
   IF (KEOF.LE.0) RETURN

!     The variable XARRAY (either single or double, depending on the
!     complile option specified, is set equal to the real*4 variable xar

   do ivar=1,3
      CALL BUFIN_sgl (kfil_s,keof_s,xarray_s(1),nlimbf)
      xarray(1:nlimbf,ivar) = xarray_s
      !if (ivar.eq.1) then
      !   do il=1,nlimbf
      !      print *, xarray(il,ivar)
      !   end do
      !endif
   end do
!
   KEOF = keof_s

   V1P = V1PBF
   V2P = V2PBF
   DVP = dvpbf
   NLIM = nlimbf
!
   RETURN
!
end subroutine SOLIN_tri
!
!
SUBROUTINE SOLIN2 (V1P,V2P,DVP,NLIM,KFILE,XARAY1,XARAY2,KEOF)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    1 November 1995
!
!                  IMPLEMENTATION:    P.D. Brown
!
!             ALGORITHM REVISIONS:    S.A. Clough
!                                     P.D. Brown
!
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!----------------------------------------------------------------------
!
!               WORK SUPPORTED BY:    THE ARM PROGRAM
!                                     OFFICE OF ENERGY RESEARCH
!                                     DEPARTMENT OF ENERGY
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE SOLIN inputs files for use with solar radiation
!     calculations for interpolation in SOLINT.  Reads files with
!     two records per panel.
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /BUFPNL_s/ V1PBF,V2PBF,DVPBF,NLIMBF
!
   DIMENSION PNLHDR(2),XARAY1(*),XARAY2(2)
!
   EQUIVALENCE (PNLHDR(1),V1PBF)
!
   !real*4 dvpbf,pnlhdr
   !integer*4 nlimbf
!
   CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
   IF (KEOF.LE.0) RETURN
   CALL BUFIN (KFILE,KEOF,XARAY1(1),NLIMBF)
   CALL BUFIN (KFILE,KEOF,XARAY2(1),NLIMBF)
!
   V1P = V1PBF
   V2P = V2PBF
   DVP = DVPBF
   NLIM = NLIMBF
!
   RETURN
!
end subroutine SOLIN2
!     ----------------------------------------------------------------
!
SUBROUTINE scale_solar (isolvar,scon,solcycfrac,solvar,svar)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!                         Written:    January 2017
!
!                  IMPLEMENTATION:    K. E. Cady-Pereira
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! Calculates the scaling factors  of the solar term from the options
! selected by the user
!

   use solar_cycle

   dimension solvar(2)   ! input solar variability scale factors or indices
   real indsolvar(2)     ! calculated solar index amplitude scale factors or indices
   real svar(3)          ! final quiet, facular and sunspot term scale factors

   real nsfm1_inv, nsfm1_inv_2
   real intfrac

   indsolvar(:) = 1.0
   svar(:) = 1.0

   if (isolvar.eq.1) then
      indsolvar=solvar

! Adjust amplitude scaling to be 1.0 at solar min,to be the requested indsolvar at solar max,
!  and to vary between those values at other solcycfrac.
      if (indsolvar(1).ne.1.0.or.indsolvar(2).ne.1.0) then
         if (solcycfrac.ge.0.0.and.solcycfrac.lt.solcyc_min) then
            wgt = (solcycfrac+1.0-solcyc_max)/(1.0+solcyc_min-solcyc_max)
            indsolvar(:) = indsolvar(:) + wgt * (1.0-indsolvar(:))
         endif
         if (solcycfrac.ge.solcyc_min.and.solcycfrac.le.solcyc_max) then
            wgt = (solcycfrac-solcyc_min)/(solcyc_max-solcyc_min)
            indsolvar(:) = 1.0 + wgt * (indsolvar(:)-1.0)
         endif
         if (solcycfrac.gt.solcyc_max.and.solcycfrac.le.1.0) then
            wgt = (solcycfrac-solcyc_max)/(1.0+solcyc_min-solcyc_max)
            indsolvar(:) = indsolvar(:) + wgt * (1.0-indsolvar(:))
         endif
      endif

!   Interpolate svar_f_0 and svar_s_0 from lookup tables using provided solar cycle fraction
!   Lookup tables points are located at the middle of each month, except for first and last points,
!   which correspond to the first and last days of the 11-year solar cycle.
!   Initial half interval (1)
      if (solcycfrac .le. 0.0) then
         tmp_f_0 = mgavgcyc(1)
         tmp_s_0 = sbavgcyc(1)
!   Final half interval (1)
      elseif (solcycfrac .ge. 1.0) then
         tmp_f_0 = mgavgcyc(nsolfrac)
         tmp_s_0 = sbavgcyc(nsolfrac)
!   Main whole intervals (131)
      else
         nsfm1_inv = 1.0 / (nsolfrac-2)
         nsfm1_inv_2 = nsfm1_inv/2.0
         if (solcycfrac.le.nsfm1_inv_2) then
            isfid = 1
            fraclo = 0.0
            frachi = nsfm1_inv_2
         elseif (solcycfrac.gt.(1.0-nsfm1_inv_2)) then
            isfid = nsolfrac-1
            fraclo = 1.0-nsfm1_inv_2
            frachi = 1.0
         else
            isfid = floor((solcycfrac-nsfm1_inv_2) * (nsolfrac-2)) + 2
            fraclo = (isfid-2) * nsfm1_inv+nsfm1_inv_2
            frachi = (isfid-1) * nsfm1_inv+nsfm1_inv_2
         endif
         intfrac = (solcycfrac - fraclo) / (frachi - fraclo)
         tmp_f_0 = mgavgcyc(isfid) + intfrac * (mgavgcyc(isfid+1) - mgavgcyc(isfid))
         tmp_s_0 = sbavgcyc(isfid) + intfrac * (sbavgcyc(isfid+1) - sbavgcyc(isfid))
      endif
      svar_f_0 = tmp_f_0
      svar_s_0 = tmp_s_0
   endif

   if (isolvar.eq.2) indsolvar(:) = solvar

   if (scon.eq.0.0) then
!   No solar cycle and no solar variability (Kurucz solar irradiance)
      if (isolvar .eq. -1) then
         solvar(:) = 1.0
      endif

!   Mean solar cycle with no solar variability (NRLSSI2 model solar irradiance)
!   Quiet sun, facular, and sunspot terms averaged over the mean solar cycle
!   (defined as average of Solar Cycles 13-24).
      if (isolvar .eq. 0) then
         solvar(:) = 1.0
      endif

!   Mean solar cycle with solar variability (NRLSSI2 model)
!   Facular and sunspot terms interpolated from LUTs to input solar cycle
!   fraction for mean solar cycle. Scalings defined below to convert from
!   averaged Mg and SB terms to Mg and SB terms interpolated here.
!   (Includes optional facular and sunspot amplitude scale factors)
      if (isolvar .eq. 1) then
         svar(2) = indsolvar(1) * (svar_f_0 - Foffset) / (svar_f_avg - Foffset)
         svar(3) = indsolvar(2) * (svar_s_0 - Soffset) / (svar_s_avg - Soffset)
      endif

!   Specific solar cycle with solar variability (NRLSSI2 model)
!   Facular and sunspot index terms input directly to model specific
!   solar cycle.  Scalings defined below to convert from averaged
!   Mg and SB terms to specified Mg and SB terms.
      if (isolvar .eq. 2) then
         svar(2)= (indsolvar(1) - Foffset) / (svar_f_avg - Foffset)
         svar(3) = (indsolvar(2) - Soffset) / (svar_s_avg - Soffset)
      endif

   endif

   if (scon.gt.0) then

!   No solar cycle and no solar variability (Kurucz or NRLSSI2 model solar irradiance)
!   Quiet sun, facular, and sunspot terms averaged over the mean solar cycle
!   (defined as average of Solar Cycles 13-24).
!   Scale from internal solar constant to requested solar constant.
      if (isolvar .eq. -1) then
         solvar(:) = scon/scon_kurucz
      else if (isolvar.eq.0) then
         solvar(:) = scon/scon_nrlssi2
      end if

!   Mean solar cycle with solar variability (NRLSSI2 model)
!   Facular and sunspot terms interpolated from LUTs to input solar cycle
!   fraction for mean solar cycle. Scalings defined below to convert from
!   averaged Mg and SB terms to Mg and SB terms interpolated here.
!   Scale internal solar constant to requested solar constant.
!   (Includes optional facular and sunspot amplitude scale factors)

      if (isolvar .eq. 1) then
         svar(1)= (scon - (indsolvar(1) * Fint + indsolvar(2) * Sint)) / Iint
         svar(2)= indsolvar(1) * (svar_f_0 - Foffset) / (svar_f_avg - Foffset)
         svar(3)= indsolvar(2) * (svar_s_0 - Soffset) / (svar_s_avg - Soffset)
      endif

   endif


   return

end subroutine scale_solar
!     ----------------------------------------------------------------

!
subroutine comb_solar (solar3,solar,nlim,isolvar,svar)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!                         Written:    January 2017
!
!                  IMPLEMENTATION:    K. E. Cady-Pereira
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! Combines the three components (quiet, facular and sunspot contributions) according
! to the options selected by the user
!

   use solar_cycle

   DIMENSION SOLAR(*)
   DIMENSION SOLAR3(2400,3)
   dimension svar(3)

   !if ((isolvar.eq.0.OR.isolvar.eq.-1).AND.scon.eq.0.0) then
   !if (isolvar.eq.0) then
   !   solar(1:nlim) = solar3(1:nlim,1)+solar3(1:nlim,2)+solar3(1:nlim,3)
   !endif

   if (isolvar.gt.0.AND.isolvar.le.2) then
      solar(1:nlim) = solar3(1:nlim,1)*svar(1)     &
      &     +solar3(1:nlim,2)*svar(2)+solar3(1:nlim,3)*svar(3)
   end if
   return

end subroutine comb_solar
!     ----------------------------------------------------------------

SUBROUTINE SOLOUT (V1P,V2P,DVP,NLIM,SOLRAD,LFILE,NPTS,NPANLS)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    3 April 1994
!
!                  IMPLEMENTATION:    P.D. Brown
!
!             ALGORITHM REVISIONS:    S.A. Clough
!                                     P.D. Brown
!
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!----------------------------------------------------------------------
!
!               WORK SUPPORTED BY:    THE ARM PROGRAM
!                                     OFFICE OF ENERGY RESEARCH
!                                     DEPARTMENT OF ENERGY
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE SOLOUT OUTPUTS ATTENUATED SOLAR RADIANCE (INTERPOLATED)
!     TO LFILE
!
   COMMON /BUFPNL_s/ V1PBF,V2PBF,DVPBF,NLIMBF
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   DIMENSION PNLHDR(2)
   DIMENSION SOLRAD(*)
!
   EQUIVALENCE (PNLHDR(1),V1PBF)
!
   REAL SOLRAD
!
   !real*4 dvpbf,pnlhdr
   !integer*4 nlimbf
!
   NPANLS = NPANLS+1
   V1PBF = V1P
   V2PBF = V2P
   DVPBF = DVP
   NLIMBF = NLIM
!
   CALL BUFOUT (LFILE,PNLHDR(1),NPHDRF)
   CALL BUFOUT (LFILE,SOLRAD(1),NLIMBF)
!
   IF (NPTS.GT.0) THEN
      IF (NPANLS.EQ.1) WRITE (IPR,900)
      WRITE (IPR,905)
      NNPTS = NPTS
      IF (NPTS.GT.(NLIMBF/2)+1) NNPTS = (NLIMBF/2)+1
      JEND = NLIMBF-NNPTS+1
      DO 10 J = 1, NNPTS
         VJ = V1PBF+ REAL(J-1)*DVPBF
         K = J+JEND-1
         VK = V1PBF+ REAL(K-1)*DVPBF
         WRITE (IPR,910) J,VJ,SOLRAD(J),K,VK,SOLRAD(K)
10    CONTINUE
   ENDIF
!
   RETURN
!
900 FORMAT ('0 ','LOCATION  WAVENUMBER',4X,'RADIANCE',13X,            &
   &        'LOCATION   WAVENUMBER',4X,                               &
   &        'RADIANCE')
905 FORMAT (' ')
910 FORMAT (I8,2X,F12.6,1P,E15.7,0P,9X,I8,2X,F12.6,1P,E15.7,0P)
!
end subroutine SOLOUT
