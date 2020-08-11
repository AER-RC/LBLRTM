!     path:      $HeadURL$
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
SUBROUTINE LBLATM
!
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
   IMPLICIT REAL*8           (V)
!
!**********************************************************************
!
!     LBLATM IS AN ATMOSPHERIC RAY TRACE PROGRAM.
!     IT CREATES AND FORMATS THE ATMOSPHERIC INPUTS FOR THE AFGL
!     LINE-BY-LINE TRANSMITTANCE/RADIANCE PROGRAM LBLRTM.
!
!     SEE THE COMMENTS IN SUBROUTINE ATMPTH FOR DETAILED INSTRUCTIONS O
!     THE USAGE OF THE ATMOSPHERIC INPUTS.
!
!     The geometry was modified for LBLRTM to reflect changes
!     implemented in MODTRAN to solve problems with inconsistent
!     path parameters.
!     These changes include changing some variables and functions to
!     double precision.
!
!**********************************************************************
!-
!-                      STATEMENT FLAGS
!-
!-    LBLATM HAS BEEN STRUCTURED TO HAVE ENHANCED PORTABILITY UNDER
!-    FORTRAN 77.  TWO FLAGS (COLUMN73) HAVE BEEN USED TO FACILITATE
!-    PROGRAM CONVERSION.
!-
!-   &    IDENTIFIES STATEMENTS REQUIRED FOR WORD SIZE LESS THAN 8 CHAR
!-               ALL STATEMENTS FLAGGED WITH & IN COLUMN 73 HAVE
!-               STARTING IN COLUMN 1. THESE TWO CHARACTERS MUST
!-               BE CHANGED TO BLANKS FOR COMPUTERS WITH WORD SIZE
!-               LESS THAN 8 CHARACTERS.
!-
!-   !    IDENTIFIES STATEMENTS REQUIRED TO DOUBLE PRECISION THE
!-               VARIABLES NEEDED FOR CALCULATIONS WHICH NEED MORE
!-               THAN 32 BITS TO OBTAIN SUFFICIENT ACCURACY (I.E.
!-               THE FREQUENCIES). STATEMENTS FLAGGED WITH ! HAVE
!-               STARTING IN COLUMN 1. THESE TWO CHARACTERS SHOULD BE
!-               CHANGED TO BLANKS FOR COMPUTERS HAVING SINGLE
!-               PRECISION LESS THAN 10 SIGNIFICANT DIGITS.
!-
!-   >    IDENTIFIES STATEMENTS THAT MAY BE USEFUL FOR CONVERSION,
!-               TYPICALLY SYSTEM SPECIFIC CALLS (I.E. DATE, TIME,
!-               CPU TIME, RANDOM NUMBER, ETC.).
!-
!----------------------------------------------------------------------
!
!     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
!     MXLAY IS THE MAXIMUM NUMBER OF OUTPUT LAYERS
!     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
!         STORED IN ZMDL (INPUT)
!     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
!         OBTAINED BY MERGING ZMDL AND ZOUT
!     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
!

!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=47,MXTRAC=22)
!
   COMMON /PATHD/ PBAR(MXLAY),TBAR(MXLAY),AMOUNT(MXMOL,MXLAY),       &
   &               WN2L(MXLAY),DVL(MXLAY),WTOTL(MXLAY),ALBL(MXLAY),   &
   &               ADBL(MXLAY),AVBL(MXLAY),H2OSL(MXLAY),IPATH(MXLAY), &
   &               ITYL(MXLAY),SECNTA(MXLAY),HT1,HT2,ALTZ(0:MXLAY),   &
   &               PZ(0:MXLAY),TZ(0:MXLAY)
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   CHARACTER*8      XID,       HMOLID,      YID
   Real*8               SECANT,       XALTZ
!
   COMMON /CVRATM/ HNAMATM,HVRATM
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WN2   ,DV ,V1 ,V2 ,TBOUND, &
   &           EMISIV,FSCDID(17),nmol_flhdr,LAYER ,YI1,YID(10),LSTWDF
!
   EQUIVALENCE (FSCDID(3),IXSCNT) , (FSCDID(5),IEMIT)
!
   CHARACTER*8      HMOLS
!
   COMMON /HMOLS/ HMOLS(MXMOL),JUNIT(MXMOL),WMOL(MXMOL),JUNITP,      &
   &               JUNITT
   COMMON /HMOLC/ HMOLC(MXMOL)
   CHARACTER*18 HNAMATM,HVRATM
   CHARACTER*8 HMOLC
   character*4 ht1,ht2
!
!     ********************************************************
!
!        NEW DATA FORMAT - GENERIC UNITS
!
!
!     *********************************************************
!
!     IRD, IPR, IPU ARE UNIT NUMBERS FOR INPUT, OUTPUT, PUNCH
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,n_zero,NOP,H1F,H2F,      &
   &                ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH,       &
   &                XVBAR, HMINF,PHIF,IERRF,HSPACE

   COMMON /c_drive/ ref_lat,hobs,ibmax_b,immax_b,                    &
   &                 lvl_1_2,jchar_st(10,2),wm(mxzmd)
!
   character*1 jchar_st
!
   CHARACTER*8      HDATE,HTIME
!
   COMMON /BNDRY/ ZBND(MXFSC),PBND(MXFSC),TBND(MXFSC),ALORNZ(MXFSC), &
   &               ADOPP(MXFSC),AVOIGT(MXFSC)
   COMMON /ZOUTP/ ZOUT(MXLAY),SOUT(MXLAY),RHOSUM(MXLAY),             &
   &               AMTTOT(MXMOL),AMTCUM(MXMOL),ISKIP(MXMOL)
!
   CHARACTER*8      HMOD

   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)

   COMMON /DEAMT/ DENM(MXMOL,MXZMD),DENP(MXMOL,MXPDIM),DRYAIR(MXZMD)


   COMMON /CNTRL/ I1,I2,I3,I4,NBNDL,I6,I7,NBNDF,I9
   COMMON /PCHINF/ MUNITS,CTYPE(MXLAY)
!
   CHARACTER*3 CTYPE
!
!     ASSIGN CVS VERSION NUMBER TO MODULE
!
   HVRATM = '$Revision$'
!
!     IBDIM IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
!     IOUTDM IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
!     IMDIM IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
!         STORED IN ZMDL (INPUT)
!     IPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH OBTAINE
!         BY MERGING ZMDL AND ZOUT
!     KDIM IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
!
   KDIM = MXMOL
   IMDIM = MXZMD
   IOUTDM = MXLAY
   IPDIM = MXPDIM
   IBDIM = MXFSC
!
   CALL LBLDAT(HDATE)
   CALL FTIME (HTIME)
   WRITE (IPR,900) HDATE,HTIME
!
   DO 10 M = 1, MXMOL
      READ (HMOLC(M),905) HMOLS(M)
      HMOLID(M) = HMOLS(M)
10 END DO
!
   IXSECT = IXSCNT/10
   CALL ATMPTH (xid,IEMIT, IXSECT)
!
   SECANT = 1.0
   nmol_flhdr = nmol
!
!     FOR IXSECT = 1, CALL XAMNTS
!

   XV1 = V1-25.
   XV2 = V2+25.
   IF (IXSECT.EQ.1) CALL XAMNTS (XV1,XV2)
!
   RETURN
!
900 FORMAT ('1',20X,'*****PROGRAM LBLATM*****     ',A10,5X,A10,///)
905 FORMAT (A8)
!
end subroutine LBLATM
!
!     ----------------------------------------------------------------
!
SUBROUTINE ATMPTH (xid,IEMIT, IXSECT)
!
!**********************************************************************
!
!
!
!
!                  ATMPTH   (ATMOSPHERIC PATH)
!
!
!
!
!                            WILLIAM O. GALLERY
!                          + GAIL   P.  ANDERSON
!                            FRANCIS X. KNEIZYS
!                            JAMES   H. CHETWYND JR.
!                            SHEPARD A. CLOUGH
!
!
!                           +(POINT OF CONTACT FOR THIS PROGRAM)
!
!                                      AIR FORCE GEOPHYSICS LAB
!                                      OPTICAL PHYSICS DIVISION
!                                      HANSCOM AFB
!                                      BEDFORD, MA.  01731
!                                      617-861-4774
!
!
!                                      REVISED:   JULY 1990
!
!**********************************************************************
!
!
!     USER INSTRUCTIONS:
!
!     ATMPTH CALCULATES THE DENSITY WEIGHTED MEAN TEMPERATURE AND
!     PRESSURE AND THE INTEGRATED ABSORBER AMOUNTS (IN MOLECULES
!     CM-2) FOR EACH LAYER ALONG A PATH THROUGH A LAYERED
!     ATMOSPHERE, INCLUDING THE EFFECTS OF REFRACTION AND THE  EARTH'S
!     CURVATURE.  ATMPTH IS DESIGNED TO PREPARE THE ATMOSPHERIC INPUTS
!     TO THE PROGRAM LBLRTM WHICH DOES A LINE-BY-LINE CALCULATION OF
!     ATMOSPHERIC TRANSMITTANCE OR RADIANCE AND IS DESCRIBED IN
!     REFERENCE (1).  THE CONTROL CARDS REQUIRED TO RUN ATMPTH ARE
!     DESCRIBED LATER IN THESE COMMENTS.  A DETAILED DESCRIPTION
!     OF THE ALGORITHM USED HERE AND A DISCUSSION OF THE EFFECTS OF
!     THE EARTH'S CURVATURE AND REFRACTION ARE GIVEN IN REFERENCE (2).
!
!     THE DEFINITIONS AND USES OF THE PATH PARAMETERS ITYPE, H1, H2,
!     ANGLE, RANGE, BETA, AND LEN ARE DESCRIBED IN REFERENCE (2) AND
!     ARE THE SAME AS IN REFERENCE (4).
!
!     THERE ARE SIX BUILT IN ATMOSPHERIC PROFILES WHICH DEFINE THE
!     PRESSURE, TEMPERATURE, AND MIXING RATIOS OF THE 28 MOLECULAR
!     SPECIES INCLUDING H2O, CO2, O3, N2O, CO, CH4, AND O2 ON THE AFGL
!     ATMOSPHERIC LINE PARAMETERS COMPILATION AT 50 STANDARD
!     ALTITUDES.  THESE MODEL ATMOSPHERES ARE DESCRIBED IN
!     REFERENCE (3).  THE USER MAY ALSO INPUT AN ATMOSPHERIC
!     PROFILE AS DESCRIBED LATER (SEE ALSO THE COMMENTS IN
!     THE SUBROUTINE NSMDL). TWENTY-0NE ADDITIONAL MIXING RATIO PROFILE
!     FOR SPECIES CORRESPONDING TO THE MOLECULES ON THE AFGL TRACE GAS
!     COMPILATION ARE INCLUDED.
!
!     THE PRINCIPAL OUTPUT CONSISTS OF THE INTEGRATED ABSORBER AMOUNTS
!     FOR A SET OF LAYERS TO BE INPUT TO THE LINE-BY-LINE CALCULATION.
!     THE NUMBER OF THESE LAYERS REPRESENTS A TRADEOFF BETWEEN ACCURACY
!     AND COMPUTATIONAL SPEED OF THE LINE-BY-LINE CALCULATION.  THE
!     USER HAS THE OPTION OF INPUTTING HIS OWN SET OF LAYER BOUNDARIES
!     OR OF LETTING THE SUBROUTINE AUTLAY GENERATE THESE LAYERS
!     AUTOMATICALLY.  IF THE USER INPUTS BOUNDARY ALTITUDES,  THEY NEED
!     NOT FALL ON THE ATMOSPHERIC PROFILE BOUNDARIES OR INCLUDE THE
!     PATH ENDPOINTS. IF AUTOMATIC LAYERING IS SELECTED, THE USER MAY
!     SPECIFY THE MAXIMUM HALFWIDTH RATIO ACROSS A LAYER AND THE
!     MAXIMUM TEMPERATURE DIFFERENCE ACROSS A LAYER.
!
!     IT IS DIFFICULT TO SPECIFY APRIORI THE RELATIONSHIP BETWEEN
!     THE NUMBER OF LAYERS AND THE ACCURACY:  THE ACCURACY DEPENDS UPON
!     SUCH FACTORS AS THE SPECTRAL REGION, THE DISTRIBUTION OF THE
!     MOLECULES OF INTEREST, THE PARTICULAR PATH TAKEN, AND WHETHER
!     TRANSMITTANCE OR RADIANCE IS CALCULATED. THE LAYERING CREATED
!     BY THE DEFAULT VALUES OF AVTRAT (1.5) AND TDIFF1 (5.0 K) AND
!     TDIFF2 (8.0 K) SHOULD BE CONSIDERED A POINT OF DEPARTURE FOR
!     SUBSEQUENT CALCULATIONS. THE USER SHOULD THEN EXPERIMENT WITH
!     DIFFERENT LAYERING UNTIL THE RESULTS ARE CONSISTENT WITH
!     HIS ACCURACY REQUIREMENTS.
!
!     TO SAVE COMPUTER TIME IN LBLRTM, THE LAYER AMOUNTS ARE ZEROED
!     OUT WHEN
!         1.  THE CUMULATIVE AMOUNT FOR THAT LAYER AND ABOVE IS LESS
!             THAN 0.1 PERCENT OF THE TOTAL,
!         AND
!         2.  A. TRANSMITTANCE IS CALCUALTED (IEMIT = 0)
!             OR
!             B. RADIANCE IS CALCULATED (IEMIT = 1) AND THE PATH IS
!                LOOKING UP ( IPATH = 3)
!     O2 IS  NOT CONSIDERED IN THIS SCHEME.  IF THE ABSORBER
!     FOR A LAYER FOR ALL THE MOLECULES (EXCEPT O2) ARE ZEROED
!     OUT, THEN THAT LAYER AND THOSE ABOVE ARE ELIMINATED
!
!     TO CALCULATE THE AMOUNTS FOR THE TRACE GASES (MOLECULES 8 THROUGH
!     31) THE USER MUST INCREASE NMOL ON CARD 3.1.
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!     OUTPUT :
!
!     THE PRINTED OUTPUT IS ON FILE IPR (DEFAULT=6). SELECTING
!     NOPRNT=1 SUPRESSES THE PRINTING OF THE ATMOSPHERIC PROFILES
!     AND THE LAYER-BY-LAYER RESULTS FOR THE REFRACTED PATH.
!     IF IPUNCH = 1, THEN THE LBLRTM INPUT DATA IS ALSO PUT ON FILE
!     IPU (DEFAULT=7) AND CONSISTS OF A SINGLE CARD IMAGE GIVING THE
!     NUMBER OF LAYERS LMAX AND A 70 CHARACTER FIELD DESCRIBING THE
!     PROFILE AND THE PATH, FOLLOWED BY TWO (OR MORE) CARD IMAGES FOR
!     EACH OF THE LMAX LAYERS
!
!        CARD 2.1    IFORM,LMAX,NMOL,SECNT0,HMOD (1X,I1,I3,I5,F10.6,3A8)
!             IFORM  = COLUMN AMOUNT FORMAT FLAG
!             LMAX   = NUMBER OF LBLRTM LAYERS, MAY DIFFER FROM
!                      IBMAX DEPENDING ON THE PATH.
!             NMOL   = NUMBER OF MOLECULES SELECTED
!             SECNT0 = EFFECTIVE SECANT (SCALE FACTOR) FOR THE AMOUNTS
!             HMOD   = 24 CHARACTER FIELD.
!
!        CARD 2.1.1  PBAR(L),TBAR(L),SECNTK(L),ITYL(L),IPATH(L),
!                    ALTZ(L-1),PZ(L-1),TZ(L-1),ALTZ(L),PZ(L),TZ(L)
!                  (E15.7,2F10.4,A3,I2,1X,F7.2,F8.3,F7.2,F7.2,F8.3,F7.2)
!             PBAR   =  AVERAGE PRESSURE (MB)
!             TBAR   =  AVERAGE TEMPERATURE (K)
!             SECNTK = SCALE FACTOR FOR COLUMN AMOUNT (DEFAULT=0)
!             ITYL  : OVERRIDES THE LBLRTM INTERNAL CALCULATION FOR
!                     ITYPE, NORMALLY LEFT BLANK
!             IPATH : IF THE PATH DOES NOT GO THROUGH A TANGENT HEIGHT,
!                         IF H1.LT.H2   IPATH = 3
!                         IF H1.GT.H2   IPATH = 1
!                      IF THE PATH GOES THROUGH A TANGENT HEIGHT, THEN
!                         FOR THE LAYERS FROM THE TANGENT HEIGHT TO
!                         MIN(H1,H2),   IPATH = 2
!                         FOR THE LAYERS (IF ANY) FROM MIN(H1,H2)
!                         TO H1,  IPATH = 1
!                         FOR THE LAYERS (IF ANY) FROM MIN(H1,H2)
!                         TO H2,  IPATH = 3
!                      FOR A HORIZONTAL PATH,  IPATH = 0
!             ALTZ(L)    UPPER BOUNDARY ALTITUDE (CURRENT LAYER)
!             ALTZ(L-1)  LOWER BOUNDARY ALTITUDE (FOR FIRST LAYER ONLY)
!             PZ(L)      PRESSURE AT ALTZ(L), MB
!             PZ(L-1)    PRESSURE AT ATLZ(L-1),  (FOR FIRST LAYER ONLY)
!             TZ(L)      TEMPERATURE AT ALTZ(L), DEGREES K
!             TZ(L-1)    TEMPERATURE AT ALTZ(L-1),(FOR FIRST LAYER ONLY
!
!        CARD 2.1.2           (AMOUNT(K,L),K=1,7),WBROADL(L)
!                             (1P8E15.7)
!        CARD 2.1.3           (AMOUNT(K,L),K=8,NMOL)
!                             (1P8E15.7)
!             AMOUNT(K)   COLUMN DENSITIES FOR THE K'TH MOLECULAR
!                         SPECIES (MOLECULES CM-2)
!             WBROADL(L)  COLUMN DENSITY FOR BROADENING GASES
!                         (MOLECULES CM-2)
!
!        CARDS 2.1 ARE REPEATED UNITL LMAX LAYERS ARE SPECIFIED.
!
!----------------------------------------------------------------------
!
!  REFERENCES:
!
! (1) LBLRTM - A USERS' GUIDE (AVAILABLE FROM S.A. CLOUGH AT
!                    THE ABOVE ADDRESS)
!        SEE ALSO:
!          FASCODE - FAST ATMOSPHERIC SIGNATURE CODE
!          (SPECTRAL TRANSMITTANCE AND RADIANCE)
!                                                       AFGL-TR-78-0081
!
! (2) AIR MASS COMPUTER PROGRAM FOR ATMOSPHERIC TRANSMITTANCE/RADIANCE:
!     FSCATM
!        W. O. GALLERY, F. X. KNEIZYS, AND S. A. CLOUGH
!                                                       AFGL-TR-83-0065
!
! (3) AFGL ATMOSPHERIC CONSTITUENT PROFILES (0-120 KM)
!        G. P. ANDERSON, S. A. CLOUGH, F.X. KNEIZYS, J. H. CHETWYND
!        AND E. P. SHETTLE
!                                                       AFGL-TR-86-0110
!
! (4) ATMOSPHERIC TRANSMITTANCE/RADIANCE:
!     COMPUTER CODE LOWTRAN 5
!                                                       AFGL-TR-80-0067
!
!**********************************************************************
!
   USE phys_consts, ONLY: pi, clight, avogad, alosmt, gascon
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
   PARAMETER (NXPBAR=MXLAY*(14+MXMOL)+2,NXZOUT=MXLAY*3+MXMOL*3)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),          &
   &                MXPNL1(MXLAY),MSLAY1,ISFILE,JSFILE,KSFILE,        &
   &                LSFILE,MSFILE,IEFILE,JEFILE,KEFILE
   COMMON /MSCONS/ AIRMSS(MXLAY),TGRND,SEMIS(3),HMINMS,HMAXMS,       &
   &                MSFLAG,MSWIT,IODFIL,MSTGLE
   COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,n_zero,NOP,H1F,H2F,      &
   &                ANGLEF,RANGEF,BETAF,LENF,V1,V2,RO,IPUNCH,XVBAR,   &
   &                HMINF,PHIF,IERRF,HSPACE

   COMMON /c_drive/ ref_lat,hobs,ibmax_b,immax_b,                    &
   &                 lvl_1_2,jchar_st(10,2),wm(mxzmd)
!
   character*1 jchar_st
!
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
!
!     BLANK COMMON FOR ZMDL
!
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
   COMMON /DEAMT/ DENM(MXMOL,MXZMD),DENP(MXMOL,MXPDIM),DRYAIR(MXZMD)
!
   CHARACTER*8      HMOLS
!
   COMMON /HMOLS/ HMOLS(MXMOL),JUNIT(MXMOL),WMOL(MXMOL),JUNITP,      &
   &               JUNITT
   COMMON /HMOLC/ HMOLC(MXMOL)
   CHARACTER*8 HMOLC
!
!     ********************************************************
!
!        NEW DATA FORMAT : GENERIC UNITS
!
!
!     *********************************************************
!
   COMMON /PATHD/ PBAR(MXLAY),TBAR(MXLAY),AMOUNT(MXMOL,MXLAY),       &
   &               WN2L(MXLAY),DVL(MXLAY),WTOTL(MXLAY),ALBL(MXLAY),   &
   &               ADBL(MXLAY),AVBL(MXLAY),H2OSL(MXLAY),IPATH(MXLAY), &
   &               ITYL(MXLAY),SECNTA(MXLAY),HT1,HT2,ALTZ(0:MXLAY),   &
   &               PZ(0:MXLAY),TZ(0:MXLAY)
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID

   COMMON /BNDRY/ ZBND(MXFSC),PBND(MXFSC),TBND(MXFSC),ALORNZ(MXFSC), &
   &              ADOPP(MXFSC),AVOIGT(MXFSC)
   COMMON /ZOUTP/ ZOUT(MXLAY),SOUT(MXLAY),RHOSUM(MXLAY),             &
   &               AMTTOT(MXMOL),AMTCUM(MXMOL),ISKIP(MXMOL)
   COMMON /PCHINF/ MUNITS,CTYPE(MXLAY)
   COMMON /FIXITYL/ IFXTYP
   DIMENSION XPBAR(NXPBAR),XZOUT(NXZOUT),WMT(MXMOL)
   DIMENSION TTMP(2),WVTMP(2),PTMP(2),ZTMP(2)

   COMMON /IADFLG/ NSPCRT,IMRGSAV
   dimension densave(mxzmd)
!
   EQUIVALENCE (PBAR(1),XPBAR(1)) , (ZOUT(1),XZOUT(1))
!
   character*8 xid(10)
   CHARACTER*48 CFORM1,CFORM2
   CHARACTER*8 COTHER
   CHARACTER*7 PAFORM(2)
   CHARACTER*4 HT1HRZ,HT2HRZ,HT1SLT,HT2SLT,PZFORM(5),  ht1,ht2
   CHARACTER*3 CTYPE
   CHARACTER*10 SREF_LAT
   REAL REF_LAT
   REAL DEN_AJ(MXMOL, MXLAY), A_AJ
   INTEGER i_aj
!

   DATA COTHER / 'OTHER   '/
   DATA AVRATS / 1.5 /,TDIF1S / 5.0 /,TDIF2S / 8.0 /
   DATA HT1HRZ / ' AT '/,HT2HRZ / ' KM '/,HT1SLT / ' TO '/,          &
   &     HT2SLT / ' KM '/
   DATA PZFORM / 'F8.6','F8.5','F8.4','F8.3','F8.2'/
   DATA PAFORM / '1PE15.7','  G15.7'/
   DATA CFORM1 / '(1PE15.7,0PF10.2,10X,A3,I2,1X,2(F7.3,F8.3,F7.2))'/
   DATA CFORM2 / '(  G15.7,0PF10.2,10X,A3,I2,23X,(F7.3,F8.3,F7.2))'/
   DATA IERROR / 0 /, IPASS / 0 /
   DATA T296 /296.0/
!
!     IAMT = 1: CALCULATE AMOUNTS, IAMT = 2: DO NOT CALCULATE AMOUNTS
!
   DATA IAMT / 1 /
!
!     AIRMS1 IS ONE AIRMASS OR THE TOTAL AMOUNT FOR A VERTICAL PATH
!     FROM GROUND TO SPACE
!
   DATA AIRMS1 / 2.153E25 /
   SECNT0 = 1.0
!
   DEG = 180.0/PI
!
!     GCAIR IS THE GAS CONSTANT FOR RHO IN MOL CM(-3), P IN MB, AND
!     T IN K
!
   GCAIR = 1.0E-3*GASCON/AVOGAD
!
!     ADCON IS THE CONSTANT FOR THE DOPPLER HALFWIDTH
!
   ADCON = SQRT(2.0* LOG(2.0)*GASCON/CLIGHT**2)
!
!     ZERO OUT COMMON BLOCKS
!
   DO 10 N = 1, MXMOL
      WMT(N) = 0.0
10 END DO
!
!       COMMON /DEAMT/ DRYAIR with BLANK COMMON RFNDXM
!
   DO 20 N = 1, IMDIM
      DRYAIR(N) = 0.
      RFNDXM(N) = 0.
20 END DO
   DO 30 N = 1, IPDIM
      IF (N.LE.IPDIM-2) THEN
         ZPTH(N) = 0.
         PP(N) = 0.
         TP(N) = 0.
         RFNDXP(N) = 0.
         SP(N) = 0.
         PPSUM(N) = 0.
         TPSUM(N) = 0.
         RHOPSM(N) = 0.
      ENDIF
!
!     COMMON /DEAMT/ DENP  with BLANK COMMON AMTP
!
      DO 28 M = 1, KDIM
         DENP(M,N) = 0.
         AMTP(M,N) = 0.
28    CONTINUE
30 END DO
!
!     /PATHD/
!
   DO 40 N = 1, NXPBAR
      XPBAR(N) = 0.0
40 END DO
!
!     /ZOUT/

   DO 50 N = 1, NXZOUT
      XZOUT(N) = 0.0
50 END DO
!
   IF (IREAD.LE.0) THEN
!
!     READ CONTROL CARD 3.1
!

      READ (IRD,900) MODEL,ITYPE,IBMAX_B,n_zero,NOPRNT,NMOL,IPUNCH,  &
         IFXTYP,MUNITS,RE,HSPACE,XVBAR,dumrd,sref_lat

! Set REF_LAT variable based on input character string in SREF_LAT
      IF (SREF_LAT .eq. '          ') THEN
         REF_LAT = 45.0
      ELSE
         READ(SREF_LAT,"(F10.3)") REF_LAT
      ENDIF

!     check to see that a value for dumrd, formerly co2mx,
!     has not been entered. Specified mixing ratios and accumulated
!     column amounts are now handled in lblrtm for all molecules.

      if (dumrd .ne. 0.) then
         write (*,*) ' - a value has been read for co2mx'
         write (*,*) ' - this option has been replaced by a'
         write (*,*) '   similar option for all molecular species'
         write (*,*) ' - see latest instructions  '
         stop ' lblatm.f:  co2mx '
      endif

      IBMAX = ABS(IBMAX_B)
   ENDIF
!
   NOP = NOPRNT
   RO = RE

   IF (NOPRNT.GE.0) THEN
      WRITE (IPR,902)
      WRITE (IPR,904) MODEL,ITYPE,IBMAX,n_zero,NOPRNT,NMOL,IPUNCH,   &
         IFXTYP,MUNITS,RE,HSPACE,XVBAR,REF_LAT
   ENDIF
!
   M = MODEL
   IF (NMOL.EQ.0) NMOL = KMXNOM
   IF (ITYPE.LT.1.OR.ITYPE.GT.3) GO TO 290
   IF (M.LT.0.OR.M.GT.6) GO TO 290
   IF (IBMAX.GT.IBDIM) GO TO 290
   IF (NMOL.GT.KDIM) GO TO 290

   IF (IPUNCH.ge.1 .and. ipass.eq.0) then

!        Tape7 only opened if ipu selected and first time through
      ipass = 1
      OPEN (UNIT=IPU,FILE='TAPE7',STATUS='UNKNOWN')
      write (ipu,905) ipass, xid
   endif
!
   if (ipunch.eq.2 .and. ipass.eq.1) then

!     open these files if derivatives are being calculated
      ipass = 2
      open (97,FILE='AJ_atmosphere', STATUS='UNKNOWN',FORM=          &
         'UNFORMATTED')
      write (97) xid

      if (IXSECT.GE.1) then
         open (20,file='AJ_xs_amnts', status='unknown',form =           &
            'unformatted')
      endif
   endif
!
   IF (RE.NE.0.0) GO TO 60
   RE = 6371.23
   IF (M.EQ.1) RE = 6378.39
   IF (M.EQ.4.OR.M.EQ.5) RE = 6356.91
   RO = RE
60 CONTINUE
!
   IF (HSPACE.EQ.0.) HSPACE = 100.
   IF (XVBAR.LE.0.) THEN
      XVBAR = (V1+V2)/2.
      IF (V2.LT.V1) XVBAR = V1
   ENDIF
!      if (REF_LAT .eq.0) REF_LAT = 45.

!
   IF (NOPRNT.GE.0) THEN
      WRITE (IPR,906)
      WRITE (IPR,904) MODEL,ITYPE,IBMAX,n_zero,NOPRNT,NMOL,IPUNCH,&
         IFXTYP,MUNITS,RE,HSPACE,XVBAR,REF_LAT
   ENDIF
!
   IF (ITYPE.EQ.1) THEN
!
!
!     =>   HORIZONTAL PATH SELECTED
!
!
      IF (NOPRNT.GE.0) WRITE (IPR,908)

!
!        READ IN CONTROL CARD 3.2
!
      IF (IREAD.LE.0) READ (IRD,910) H1F,RANGEF
      RANGE = RANGEF
      ZH = H1F
      H1 = ZH
      H2 = 0.
      H2F = H2
      ANGLE = 0.
      ANGLEF = ANGLE
      BETA = 0.
      BETAF = BETA
      LEN = 0
      LENF = LEN
      IF (NOPRNT.GE.0) WRITE (IPR,912) ZH,RANGE
!
!        SET UP THE ATMOSPHERIC PROFILE
!
      CALL MDLATM (ITYPE,M,IREAD,HSPACE)
!
      IF (IMMAX.EQ.1) THEN
         ZH = ZMDL(1)
         H1F = ZH
         H1 = ZH
         PH = PM(1)
         TH = TM(1)
         RHOBAR = ALOSMT*PH*TZERO/(PZERO*TH)
         DO 70 K = 1, NMOL
            DENP(K,1) = DENM(K,1)
70       CONTINUE
         DENW(1) = DENM(1,1)
         GO TO 110
      ENDIF
!
!        INTERPOLATE ATMOSPHERIC PROFILE DENSITIES TO ZH
!
      DO 80 IM = 2, IMMAX
         IF (ZH.LT.ZMDL(IM)) GO TO 90
80    CONTINUE
      IM = IMMAX
90    CONTINUE
      A = (ZH-ZMDL(IM-1))/(ZMDL(IM)-ZMDL(IM-1))
      CALL EXPINT (PH,PM(IM-1),PM(IM),A)
      TH = TM(IM-1)+(TM(IM)-TM(IM-1))*A
      RHOBAR = ALOSMT*PH*TZERO/(PZERO*TH)
      DO 100 K = 1, NMOL
         CALL EXPINT (DENP(K,1),DENM(K,IM-1),DENM(K,IM),A)
100   CONTINUE
!
110   CONTINUE
      IF (NOPRNT.GE.0) THEN
         WRITE (IPR,914) HMOD,ZH,PH,TH,(HMOLS(K),K=1,NMOL)
         WRITE (IPR,916) RHOBAR,(DENP(K,1),K=1,NMOL)
      ENDIF
!
!        COMPUTE AMOUNTS FOR A HORIZONTAL PATH
!
      DO 120 K = 1, NMOL
         AMOUNT(K,1) = DENP(K,1)*RANGE*1.0E+5
120   CONTINUE
      AMTAIR = RHOBAR*RANGE*1.0E+5
      IF (NOPRNT.GE.0) THEN
         WRITE (IPR,918) HMOD,ZH,PH,TH,RANGE,(HMOLS(K),K=1,NMOL)
         WRITE (IPR,920) AMTAIR,(AMOUNT(K,1),K=1,NMOL)
      ENDIF
      IPATH(1) = 0
      LMAX = 1
      NLAYRS = 1
!
      SUMAMT = 0.
      DO 130 K = 1, NMOL
         SUMAMT = SUMAMT+AMOUNT(K,1)
130   CONTINUE
      WN2L(1) = AMTAIR-SUMAMT
!
      PBAR(1) = PH
      TBAR(1) = TH
      ALTZ(0) = -RANGE
      ZOUT(1) = ZH
      IOUTMX = 1
      SECNTA(1) = 1.
      ALTZ(1) = ZH
      ht1 = ht1hrz
      ht2 = ht2hrz
!
!        > Write atmosphere to TAPE7 (in E15.7 format) <
!
      IF (IPUNCH.ge.1) THEN
         IFORM = 1
         WRITE (IPU,924) IFORM,LMAX,NMOL,SECNT0,HMOD,RANGE,ZH
!
!           -------------------------------------
!           > Write molecular information in    <
!           >  - mixing ratio if MUNITS is 1    <
!           >  - column density if MUNITS is 0  <
!           -------------------------------------
!
         IF (MUNITS.EQ.1) THEN
            DRAIR = WN2L(1)
            DO 135 M = 2,NMOL
               DRAIR = DRAIR + AMOUNT(M,1)
135         CONTINUE
!
!             > If DRAIR is zero, then write out AMOUNT only    <
!             > (since AMOUNT zero => mixing ratio zero)        <
!
            IF (DRAIR.EQ.0 .AND. NOPRNT.GE.0) THEN
               WRITE (IPU,926) PH,TH,IPATH(1),ZH,ZH, (AMOUNT(K,1),&
                  K=1,7),WN2L(1), (AMOUNT(K,1),K=8,NMOL)
            ELSEIF (NOPRNT.GE.0) THEN
               WRITE (IPU,926) PH,TH,IPATH(1),ZH,ZH, (AMOUNT(K,1)/&
                  DRAIR,K=1,7),WN2L(1), (AMOUNT(K,1)/DRAIR,K=8,NMOL)
            ENDIF
         ELSE
!
!             Test to make sure there are no fractional molecular
!             amounts written out (will cause PATH to assume
!             mixing ratio)
!
            DO 137 K=1,NMOL
               IF (AMOUNT(K,1).LT.1.) THEN
                  IF (NOPRNT.GE.0) WRITE(IPR,1000) K,AMOUNT(K,1)
                  AMOUNT(K,1) = 0.0
               ENDIF
137         CONTINUE
!
            IF (NOPRNT.GE.0) WRITE (IPU,926) PH,TH,IPATH(1),ZH,ZH,&
               (AMOUNT(K,1),K=1,7),WN2L(1), (AMOUNT(K,1),K=8,NMOL)
         ENDIF
      ENDIF
!
   ELSE
!
!
!     =>   SLANT PATH SELECTED
!
!
!        ITYPE = 2 OR 3: SLANT PATH THROUGH THE ATMOSPHERE
!
      IF (NOPRNT.GE.0) WRITE (IPR,930) ITYPE
!
!      >  READ IN CONTROL CARD 3.2 CONTAINING SLANT PATH PARAMETERS <
!
      IF (IREAD.LE.0) READ (IRD,932) H1F,H2F,ANGLEF,RANGEF,BETAF, &
         LENF,HOBS
      H1 = H1F
      H2 = H2F
      ANGLE = ANGLEF
      RANGE = RANGEF
      BETA = BETAF
      LEN = LENF
      IF (NOPRNT.GE.0) THEN
         IF (IBMAX_B .LT. 0) THEN
            WRITE (IPR,933) H1,H2,ANGLE,RANGE,BETA,LEN
         ELSE
            WRITE (IPR,934) H1,H2,ANGLE,RANGE,BETA,LEN
         ENDIF
      ENDIF
!
!        > GENERATE OR READ IN LBLRTM BOUNDARY LAYERS <
!
      IF (IBMAX.EQ.0) THEN
!
!           > SELECT AUTOMATIC LAYERING <
!
         IF (IREAD.LE.0) THEN
            READ (IRD,936) AVTRAT,TDIFF1,TDIFF2,ALTD1,ALTD2
            IF (AVTRAT.EQ.0.0) AVTRAT = AVRATS
            IF (TDIFF1.EQ.0.0) TDIFF1 = TDIF1S
            IF (TDIFF2.EQ.0.0) TDIFF2 = TDIF2S
            IF ((ALTD2.LE.0).OR.(ALTD2.LE.ALTD1)) THEN
               ALTD1 = 0.
               ALTD2 = 100.
            ENDIF
            IF (NOPRNT.GE.0) WRITE (IPR,938) AVTRAT,TDIFF1,TDIFF2,&
               ALTD1,ALTD2
            IF (AVTRAT.LE.1.0.OR.TDIFF1.LE.0.0.OR.TDIFF2.LE.0.0)  &
               GO TO 320
         ENDIF
         GO TO 150
!
      ENDIF
!
!        > READ IN LBLRTM BOUNDARY LAYERS <
!
      IF (IREAD.LE.0) THEN
         IF (IBMAX_B .LT. 0) THEN
            READ (IRD,940) (PBND(IB),IB=1,IBMAX)
            IF (NOPRNT.GE.0) WRITE (IPR,943) (IB,PBND(IB),IB=1,   &
               IBMAX)
         ELSE
            READ (IRD,940) (ZBND(IB),IB=1,IBMAX)
            IF (NOPRNT.GE.0) WRITE (IPR,942) (IB,ZBND(IB),IB=1,   &
               IBMAX)
         ENDIF
      ENDIF
!
      IF (IBMAX.EQ.0) GO TO 150

      IF (IBMAX_B .GT. 0) THEN
         DO 140 IB = 2, IBMAX
            IF (ZBND(IB).LE.ZBND(IB-1)) GO TO 300
140      CONTINUE
      ENDIF
      IF (IBMAX_B .LT. 0) THEN
         DO 145 IB = 2,IBMAX
            IF (PBND(IB) .GE. PBND(IB-1)) GO TO 305
145      CONTINUE
      ENDIF
150   CONTINUE

!
!        > SET UP ATMOSPHERIC PROFILE <
!

      CALL MDLATM (ITYPE,M,IREAD,HSPACE)


! INTERPOLATE PBND GRID ONTO ZBND GRID.

! TO ENSURE THAT CALCULATED/INPUT ZMDL'S WILL MATCH CALCULATED USER-LEVE
! ALTITUDES, A COMBINATION OF INTERPOLATION AND HYDROSTATICS ARE USED.
! ZBND = A * F1(P) + (1 - A) * F2(P), WHERE
! F1(P) = INTERPOLATION IN LN(P), F2(P) = HYDROSTATIC CALCULATION

      IF (IBMAX_B .LT. 0) THEN

         ISTART = 2

         DO 160 IP=1,IBMAX
            PTMP(1) = 0.0
            TTMP(1) = 0.0
            WVTMP(1) = 0.0
            ZTMP(1) = 0.0

            PTMP(2) = 0.0
            TTMP(2) = 0.0
            WVTMP(2) = 0.0
            ZTMP(2) = 0.0

            DO 161 LIP=ISTART,IMMAX
               IF (PBND(IP) .GT. PM(LIP)) GO TO 162
161         CONTINUE
            LIP=IMMAX
162         CONTINUE
            IF (PBND(IP) .EQ. PM(LIP-1)) THEN
               ZBND(IP) = ZMDL(LIP-1)
               TBND(IP) = TM(LIP-1)
            ELSE
               IF(PBND(IP) .EQ. PM(LIP)) THEN
                  ZBND(IP) = ZMDL(LIP)
                  TBND(IP) = TM(LIP)
               ELSE

! PERFORM INTERPOLATION IN LN(PM)
                  HIP = (ZMDL(LIP)-ZMDL(LIP-1))/ LOG(PM(LIP)/PM(  &
                     LIP-1))
                  ZINT = ZMDL(LIP-1)+ HIP* LOG(PBND(IP)/PM(LIP-1))

! PERFORM ALTITUDE CALCULATION USING HYDROSTATIC EQUATION
                  PTMP(1) = PM(LIP-1)
                  ZTMP(1) = ZMDL(LIP-1)
                  TTMP(1) = TM(LIP-1)
                  WVTMP(1) = DENW(LIP-1)

                  PTMP(2) = PBND(IP)

                  TIP = (TM(LIP)-TM(LIP-1))/ LOG(PM(LIP)/PM(LIP-1)&
                     )
                  TTMP(2) = TM(LIP-1)+ TIP* LOG(PBND(IP)/PM(LIP-1)&
                     )

                  WVIP = (DENW(LIP)-DENW(LIP-1))/ LOG(PM(LIP)/PM( &
                     LIP-1))
                  WVTMP(2) = DENW(LIP-1) + WVIP* LOG(PBND(IP)/PM( &
                     LIP-1))
                  CALL CMPALT(2,PTMP,TTMP, WVTMP,ZTMP(1),REF_LAT, &
                     ZTMP)

! COMBINE THE INTERPOLATION AND THE HYDROSTATIC CALCULATION

                  RATP = LOG(PBND(IP)/PM(LIP-1))/ LOG(PM(LIP)/PM( &
                     LIP-1))

! Ignore exponential interolation term (zint) if user has input a profile specified on pressure levels
                  if (immax_b.lt.0) then
                     A =0.
                  else
                     A = RATP**3
                  endif

                  ZBND(IP) = A*ZINT + (1-A)*ZTMP(2)
                  TBND(IP) = TTMP(2)
               ENDIF
            ENDIF

            ISTART = LIP

160      CONTINUE

! INTERPOLATE H1, H2 ONTO ALTITUDE GRID
         PTMP(1) = 0.0
         TTMP(1) = 0.0
         WVTMP(1) = 0.0
         ZTMP(1) = 0.0

         PTMP(2) = 0.0
         TTMP(2) = 0.0
         WVTMP(2) = 0.0
         ZTMP(2) = 0.0

         DO 166 LIP = 2,IMMAX
            IF (H1 .GT. PM(LIP)) GO TO 167
166      CONTINUE
         LIP = IMMAX
167      CONTINUE
         IF (H1 .EQ. PM(LIP-1)) THEN
            H1 = ZMDL(LIP-1)
         ELSE
            IF(H1 .EQ. PM(LIP)) THEN
               H1 = ZMDL(LIP)
            ELSE

! PERFORM INTERPOLATION IN LN(PM)
               HIP = (ZMDL(LIP)-ZMDL(LIP-1))/ LOG(PM(LIP)/PM(LIP- &
                  1))
               ZINT = ZMDL(LIP-1)+ HIP* LOG(H1/PM(LIP-1))

! PERFORM ALTITUDE CALCULATION USING HYDROSTATIC EQUATION
               PTMP(1) = PM(LIP-1)
               ZTMP(1) = ZMDL(LIP-1)
               TTMP(1) = TM(LIP-1)
               WVTMP(1) = DENW(LIP-1)

               PTMP(2) = H1

               TIP = (TM(LIP)-TM(LIP-1))/ LOG(PM(LIP)/PM(LIP-1))
               TTMP(2) = TM(LIP-1)+ TIP* LOG(H1/PM(LIP-1))

               WVIP = (DENW(LIP)-DENW(LIP-1))/ LOG(PM(LIP)/PM(LIP-&
                  1))
               WVTMP(2) = DENW(LIP-1) + WVIP* LOG(H1/PM(LIP-1))
               CALL CMPALT(2,PTMP,TTMP, WVTMP,ZTMP(1),REF_LAT,    &
                  ZTMP)

! COMBINE THE INTERPOLATION AND THE HYDROSTATIC CALCULATION

               RATP = LOG(H1/PM(LIP-1))/ LOG(PM(LIP)/PM(LIP-1))
               if (immax_b.lt.0) then
                  A =0.
               else
                  A = RATP**3
               endif


               H1 = A*ZINT + (1-A)*ZTMP(2)
            ENDIF
         ENDIF

         IF (H1 .LT. 0.0) THEN
            PRINT 946, H1,ZTMP(1)
            IF (NOPRNT.GE.0) WRITE (IPR,946) H1,ZTMP(1)
            STOP ' COMPUTED ALTITUDE VALUE OF H1 IS NEGATIVE'
         ENDIF

         PTMP(1) = 0.0
         TTMP(1) = 0.0
         WVTMP(1) = 0.0
         ZTMP(1) = 0.0

         PTMP(2) = 0.0
         TTMP(2) = 0.0
         WVTMP(2) = 0.0
         ZTMP(2) = 0.0

         DO 168 LIP = 2,IMMAX
            IF (H2 .GT. PM(LIP)) GO TO 169
168      CONTINUE
         LIP = IMMAX
169      CONTINUE
         IF (H2 .EQ. PM(LIP-1)) THEN
            H2 = ZMDL(LIP-1)
         ELSE
            IF(H2 .EQ. PM(LIP)) THEN
               H2 = ZMDL(LIP)
            ELSE
! PERFORM INTERPOLATION IN LN(PM)
               HIP = (ZMDL(LIP)-ZMDL(LIP-1))/ LOG(PM(LIP)/PM(LIP- &
                  1))
               ZINT = ZMDL(LIP-1)+ HIP* LOG(H2/PM(LIP-1))

! PERFORM ALTITUDE CALCULATION USING HYDROSTATIC EQUATION
               PTMP(1) = PM(LIP-1)
               ZTMP(1) = ZMDL(LIP-1)
               TTMP(1) = TM(LIP-1)
               WVTMP(1) = DENW(LIP-1)

               PTMP(2) = H2

               TIP = (TM(LIP)-TM(LIP-1))/ LOG(PM(LIP)/PM(LIP-1))
               TTMP(2) = TM(LIP-1)+ TIP* LOG(H2/PM(LIP-1))

               WVIP = (DENW(LIP)-DENW(LIP-1))/ LOG(PM(LIP)/PM(LIP-&
                  1))
               WVTMP(2) = DENW(LIP-1) + WVIP* LOG(H2/PM(LIP-1))
               CALL CMPALT(2,PTMP,TTMP, WVTMP,ZTMP(1),REF_LAT,    &
                  ZTMP)

! COMBINE THE INTERPOLATION AND THE HYDROSTATIC CALCULATION

               RATP = LOG(H2/PM(LIP-1))/ LOG(PM(LIP)/PM(LIP-1))

               if (immax_b.lt.0) then
                  A =0.
               else
                  A = RATP**3
               endif

               H2 = A*ZINT + (1-A)*ZTMP(2)
            ENDIF
         ENDIF

         IF (H2 .LT. 0.0) THEN
            PRINT 946, H2,ZTMP(1)
            IF (NOPRNT.GE.0) WRITE (IPR,946) H2,ZTMP(1)
            STOP ' COMPUTED ALTITUDE VALUE OF H2 IS NEGATIVE'
         ENDIF
      ENDIF

      IERB = 0
      IF (IBMAX.GE.1) THEN
         IF (ZBND(1).LT.ZMDL(1)) THEN
            IERB = 1
            IF (NOPRNT.GE.0) WRITE (IPR,944)
            IF (ABS(ZBND(1)-ZMDL(1)).LE.0.0001) THEN
               ZBND(1) = ZMDL(1)
            ELSE
               PRINT 946,ZBND(1),ZMDL(1)
               IF (NOPRNT.GE.0) WRITE (IPR,946) ZBND(1),ZMDL(1)
               STOP ' BOUNDARIES OUTSIDE OF ATMOS'
            ENDIF
         ENDIF
      ENDIF
!
!
!        > COMPUTE THE REFRACTIVE INDEX PROFILE        <
!        > RFNDXM IS 1.0-INDEX                         <
!        > EQUATION FOR RFNDXM IS FROM LOWTRAN (REF 3) <
!
      IF (NOPRNT.GE.0)                                               &
      &        WRITE(IPR,*) '   - Using LOWTRAN6 refractive index -'
!
      DO 170 IM = 1, IMMAX
         PPH2O = DENM(1,IM)*PZERO*TM(IM)/(TZERO*ALOSMT)
!
!	    Approximation to refraction index (from LOWTRAN5)
!
!           RFNDXM(IM) = ((77.46+0.459E-8*XVBAR**2)*PM(IM)/TM(IM)-
!    *                   (PPH2O/1013.0)*(43.49-0.347E-8*XVBAR**2))*
!    *                   1.0E-6
!
!	    Approximation to refraction index (from LOWTRAN6)
!
         RFNDXM(IM)=((83.42+(185.08/(1.0-(XVBAR/1.14E+5)**2))+    &
            (4.11/(1.0-(XVBAR/6.24E+4)**2)))*(PM(IM)*288.15)/        &
            (1013.25*TM(IM))-(43.49-(XVBAR/1.7E+4)**2)*(PPH2O/       &
            1013.25)) *1.0E-06
170   CONTINUE
!
!        > PRINT ATMOSPHERIC PROFILE <
!
      IF (NOPRNT.GE.0) WRITE (IPR,948) M,HMOD
!
      IF (NOPRNT.EQ.0) THEN
         WRITE (IPR,950) (HMOLS(K),K=1,NMOL)
         WRITE (IPR,952)
      ENDIF
!
      DO 180 IM = 1, IMMAX
!
!        > DENG=DENM(1,IM)*2.989641E-17 <
!
         DENAIR = ALOSMT*(PM(IM)/PZERO)*(TZERO/TM(IM))
         densave(im) = denair
         IF (NOPRNT.EQ.0) WRITE (IPR,954) IM,ZMDL(IM),PM(IM),TM(  &
            IM),RFNDXM(IM), DENAIR,(DENM(K,IM),K=1,NMOL)


180   CONTINUE

      IF (NOPRNT.EQ.0) WRITE (IPR,951) (HMOLS(K),K=1,NMOL)

      DO 188 IM = 1, IMMAX

!     Calculate mixing ratio, using dry air

         dry_air = densave(im)-denm(1,im)
         IF (NOPRNT.EQ.0) WRITE (IPR,954) IM,ZMDL(IM),PM(IM),TM(  &
            IM),RFNDXM(IM), DENsave(im),((DENM(K,IM)/DRY_AIR)*1.e6,K=&
            1,NMOL)

188   continue

!_______________________________________________________________________

!
!        > REDUCE SLANT PATH PARAMETERS TO STANDARD FORM <
!
      CALL FSCGEO (H1,H2,ANGLE,RANGE,BETA,ITYPE,LEN,HMIN,PHI,     &
         IERROR, HOBS)
      IF (IERROR.NE.0) GO TO 310
!
!        > SET UP LBLRTM LAYER BOUNDARIES <
!
      IF (IBMAX.NE.0) GO TO 200
!
!        > AUTOMATIC LAYERING SELECTED <
!
      HMAX = MAX(H1,H2)
      CALL AUTLAY (HMIN,HMAX,XVBAR,AVTRAT,TDIFF1,TDIFF2,ALTD1,    &
         ALTD2, IERROR)
      GO TO 220
200   CONTINUE
!
!        > USER SUPPLIED LAYERING <
!
      if (noprnt .ge. 0) WRITE (IPR,956)
      IF (IBMAX_B .LT. 0) THEN
         DO 205 IB = 1, IBMAX
            CALL HALFWD_P(ZBND(IB),XVBAR,PBND(IB),TBND(IB),       &
               ALORNZ(IB),ADOPP(IB),AVOIGT(IB))
205      CONTINUE
      ELSE
         DO 210 IB = 1, IBMAX
            CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),        &
               ALORNZ(IB),ADOPP(IB),AVOIGT(IB))
210      CONTINUE
      ENDIF
220   CONTINUE
      if (noprnt .ge. 0) WRITE (IPR,958) ALZERO,AVMWT,XVBAR
      DO 230 IB = 1, IBMAX
         ZETA = ALORNZ(IB)/(ALORNZ(IB)+ADOPP(IB))
         RATIO = 0.0
         DTEMP = 0.0
         IF (IB.NE.IBMAX) RATIO = AVOIGT(IB)/AVOIGT(IB+1)
         IF (IB.NE.IBMAX) DTEMP = ABS(TBND(IB)-TBND(IB+1))
         if (noprnt .ge. 0) WRITE (IPR,960) IB,ZBND(IB),PBND(IB),&
            TBND(IB),ALORNZ(IB), ADOPP(IB),ZETA,AVOIGT(IB),RATIO,    &
            DTEMP
230   CONTINUE
      IF (IERROR.NE.0) STOP ' IERROR'
!
!        > CALCULATE THE REFRACTED PATH THROUGH THE ATMOSPHERE <
!
      CALL RFPATH (H1,H2,ANGLE,PHI,LEN,HMIN,IAMT,RANGE,BETA,      &
         BENDNG)
!
!        > PRINT AMOUNTS BY LAYER AND ACCUMULATE TOTALS <
!
      IF (NOPRNT.ge.0) WRITE (IPR,962) (HMOLS(K),K=1,NMOL)
      I2 = IPMAX-1
      AIRTOT = 0.0
      DO 240 K = 1, NMOL
         AMTTOT(K) = 0.0
240   CONTINUE
      HMID = MIN(H1,H2)
      DO 260 I = 1, I2
         FAC = 1.0
         IF (LEN.EQ.1.AND.ZPTH(I+1).LE.HMID) FAC = 2.0
         AMTAIR = RHOPSM(I)*1.0E5
         AIRTOT = AIRTOT+FAC*AMTAIR
         DO 250 K = 1, NMOL
            AMTTOT(K) = AMTTOT(K)+FAC*AMTP(K,I)
250      CONTINUE
         IF (NOPRNT.ge.0) WRITE (IPR,964) I,ZPTH(I),ZPTH(I+1),    &
            AMTAIR, (AMTP(K,I),K=1,NMOL)
260   CONTINUE
      IF (NOPRNT.ge.0) WRITE (IPR,966) H1,H2,AIRTOT,(AMTTOT(K),K= &
         1,NMOL)
!
!        > PRINT SUMMARY OF PATH <
!
      AIRMAS = AIRTOT/AIRMS1
      if (noprnt .ge. 0) WRITE (IPR,968) HMOD,H1,H2,ANGLE,RANGE,  &
         BETA,PHI,HMIN,BENDNG,LEN,AIRMAS
      IF (ITYPE.EQ.3) ITYPE = 2
      H1F = H1
      H2F = H2
      ANGLEF = ANGLE
      PHIF = PHI
      IERRF = IERROR
      LENF = LEN
      HMINF = HMIN
!
!        > CONDENSE THE AMOUNTS INTO THE LBLRTM OUTPUT LAYERS ZOUT, <
!        > WHICH ARE DEFINED BY THE BOUNDARIES ZBND FROM HMIN TO    <
!        > HMAX ALSO, ZERO OUT THE AMOUNT FOR A MOLECULE IF THE     <
!        > CUMULATIVE AMOUNT FOR THAT LAYER AND ABOVE IN LESS THAN  <
!        > 0.1 PERCENT OF THE TOTAL                                 <
!
      CALL FPACK (H1,H2,HMID,LEN,IEMIT,n_zero)
!
!        > OUTPUT THE PROFILE <
!
      LMAX = IOUTMX-1
      IF (NMOL.LE.7) THEN
         if (noprnt.ge.0) WRITE (IPR,970) (HMOLS(K),K=1,NMOL),    &
            COTHER
      ELSE
         if (noprnt.ge.0) WRITE (IPR,970) (HMOLS(K),K=1,7),COTHER,&
            (HMOLS(K),K=8,NMOL)
      ENDIF
      IF (IPUNCH.ge.1) THEN
         IFORM = 1
         WRITE (IPU,972) IFORM,LMAX,NMOL,SECNT0,(HMOD(I),I=1,2),  &
            H1,H2,ANGLE,LEN
      ENDIF

      IF (IPUNCH.eq.2) THEN
         WRITE (97) LMAX,NMOL,SECNT0,(HMOD(I),I=1,2),H1,H2,ANGLE, &
            LEN
      ENDIF

      SUMN2 = 0.
      SUMRS = 0.
      PWTD = 0.
      TWTD = 0.
      WTOT = 0.
      P_h2o = 0.
      T_h2o = 0.
!
!
!        Read from/Write to "IFIXTYPE" file: if IFXTYP = -2, use
!        preset ITYL values; if IFXTYP = 2, calculate and write out
!        ITYL values.

      IF (IFXTYP.eq.-2) then
         open(99,file='IFIXTYPE',status='old', form='formatted')
      elseif (IFXTYP.eq.2) then
         open(99,file='IFIXTYPE',status='unknown', form=          &
            'formatted')
      endif

      DO 280 L = 1, LMAX
         FACTOR = 1.
         IF (IPATH(L).EQ.2) FACTOR = 2.
         SUMWK = 0.
         DO 270 MOL = 1, NMOL
            SUMWK = SUMWK+AMOUNT(MOL,L)
            WMT(MOL) = WMT(MOL)+AMOUNT(MOL,L)*FACTOR
270      CONTINUE
         WTOTL(L) = SUMWK+WN2L(L)
         SUMN2 = SUMN2+WN2L(L)*FACTOR
         SUMRS = SUMRS+RHOSUM(L)*FACTOR
         WTOT = WTOT+WTOTL(L)*FACTOR
         PWTD = PWTD+PBAR(L)*WTOTL(L)*FACTOR
         TWTD = TWTD+TBAR(L)*WTOTL(L)*FACTOR
         P_h2o = P_h2o+PBAR(L)*amount(1,l)*FACTOR
         T_h2o = T_h2o+TBAR(L)*amount(1,l)*FACTOR
!
!           > Determine ITYL(L), if desired (when setting the ratio <
!           > from layer to layer).  Default is ITYL(L) left blank, <
!           >                                                       <
!           >                    CTYPE(L) = '   '                   <
!
         CTYPE(L) = '   '
         IF (IFXTYP.EQ.1) THEN
            FRH2O = AMOUNT(1,L)/WTOTL(L)
            ALFCOR = (PBAR(L)/PZERO)*SQRT(T296/TBAR(L))
            ADBAR = 3.581155E-07*XVBAR*SQRT(TBAR(L)/AVMWT)
            CALL FIXTYP(IEMIT,FRH2O,ALFCOR,OLDDV,L,CTYPE(L))
            read(ctype(L),1100) ityl(l)
         elseif (ifxtyp.eq.2) then
            FRH2O = AMOUNT(1,L)/WTOTL(L)
            ALFCOR = (PBAR(L)/PZERO)*SQRT(T296/TBAR(L))
            ADBAR = 3.581155E-07*XVBAR*SQRT(TBAR(L)/AVMWT)
            CALL FIXTYP(IEMIT,FRH2O,ALFCOR,OLDDV,L,CTYPE(L))
            read(ctype(l),1100) ityl(l)
            write(99,1100) ityl(l)
         elseif (ifxtyp.eq.-2) then
            read(99,1100) ityl(l)
            write(ctype(l),1100) ityl(l)
         ENDIF
!
!
!           > Write atmosphere to TAPE6 in column density <
!
         if (noprnt .ge. 0) then
            IF (NMOL.LE.7) THEN
               WRITE (IPR,974) L,ZOUT(L),ZOUT(L+1),CTYPE(L),IPATH(&
                  L), PBAR(L),TBAR(L),RHOSUM(L), (AMOUNT(K,L),K=1,   &
                  NMOL),WN2L(L)
            ELSE
               WRITE (IPR,976) L,ZOUT(L),ZOUT(L+1),CTYPE(L),IPATH(&
                  L), PBAR(L),TBAR(L),RHOSUM(L), (AMOUNT(K,L),K=1,7),&
                  WN2L(L), (AMOUNT(K,L),K=8,NMOL)
            ENDIF
         endif
!
!           > Write atmosphere to TAPE7 <
!
         IF (IPUNCH.ge.1) THEN
            LTST = L
            IF (L.EQ.1) LTST = 0
            PTST = LOG10(PZ(LTST))
            NPTST = PTST+2
            IF (PTST.LT.0.0) NPTST = 1
            CFORM1(38:41) = PZFORM(NPTST)
            CFORM2(38:41) = PZFORM(NPTST)
            NPTST = 1
            IF (PBAR(L).GE.0.1) NPTST = 2
            CFORM1(2:8) = PAFORM(NPTST)
            CFORM2(2:8) = PAFORM(NPTST)
            IF (L.EQ.1) THEN
               WRITE (IPU,CFORM1) PBAR(L),TBAR(L),CTYPE(L),    &
                  IPATH(L), ALTZ(L-1),PZ(L-1),TZ(L-1), ALTZ(L),   &
                  PZ(L), TZ(L)
            ELSE
               WRITE (IPU,CFORM2) PBAR(L),TBAR(L),CTYPE(L),    &
                  IPATH(L), ALTZ(L), PZ(L), TZ(L)
            ENDIF
!
!           -------------------------------------
!           > Write molecular information in    <
!           >  - mixing ratio if MUNITS is 1    <
!           >  - column density if MUNITS is 0  <
!           -------------------------------------
!
            IF (MUNITS.EQ.1) THEN
               DRAIR = WN2L(L)
               DO 275 M = 2,NMOL
                  DRAIR = DRAIR + AMOUNT(M,L)
275            CONTINUE
!
!                 > If DRAIR is zero, then write out AMOUNT only    <
!                 > (since AMOUNT zero => mixing ratio zero)        <
!
               IF (DRAIR.EQ.0) THEN
                  WRITE (IPU,978) (AMOUNT(K,L),K=1,7),WN2L(L)
                  IF (NMOL.GT.7) WRITE (IPU,978) (AMOUNT(K,L), &
                     K=8,NMOL)
               ELSE
                  WRITE (IPU,978) (AMOUNT(K,L)/DRAIR,K=1,7),   &
                     WN2L(L)
                  IF (NMOL.GT.7) WRITE (IPU,978) (AMOUNT(K,L)/ &
                     DRAIR,K=8,NMOL)
               ENDIF
            ELSE
!
!                 Test to make sure there are no fractional molecular
!                 amounts written out (will cause PATH to assume
!                 mixing ratio)
!
               DO 277 K=1,NMOL
                  IF (AMOUNT(K,L).LT.1.) THEN
                     WRITE(IPR,1000) K,L
                     AMOUNT(K,L) = 0.0
                  ENDIF
277            CONTINUE
!
               WRITE (IPU,978) (AMOUNT(K,L),K=1,7),WN2L(L)
               IF (NMOL.GT.7) WRITE (IPU,978) (AMOUNT(K,L),K=8,&
                  NMOL)
            ENDIF
         ENDIF
280   CONTINUE

!        Close "IFIXTYPE" file, if used
      IF (abs(IFXTYP).eq.2) then
         rewind(99)
         close(99)
      endif
!
!        > Write atmosphere to TAPE6 in mixing ratio <
!
      if (noprnt .ge. 0) then
         IF (NMOL.LE.7) THEN
            WRITE (IPR,973) (HMOLS(K),K=1,NMOL),COTHER
         ELSE
            WRITE (IPR,973) (HMOLS(K),K=1,7),COTHER, (HMOLS(K),&
               K=8,NMOL)
         ENDIF
      endif
!
      DO 285 L = 1, LMAX
         DRAIR = WN2L(L)
         DO 283 M = 2,NMOL
            DRAIR = DRAIR + AMOUNT(M,L)
283      CONTINUE
!
!           > If DRAIR is zero, then write out AMOUNT only    <
!           > (since AMOUNT zero => mixing ratio zero)        <
!
         IF (DRAIR.EQ.0 .and. noprnt.ge.0) THEN
            IF (NMOL.LE.7) THEN
               WRITE (IPR,974) L,ZOUT(L),ZOUT(L+1),CTYPE(L),&
                  IPATH(L), PBAR(L),TBAR(L),RHOSUM(L), (AMOUNT(&
                  K,L),K=1,NMOL),WN2L(L)
            ELSE
               WRITE (IPR,976) L,ZOUT(L),ZOUT(L+1),CTYPE(L),&
                  IPATH(L), PBAR(L),TBAR(L),RHOSUM(L), (AMOUNT(&
                  K,L),K=1,7),WN2L(L), (AMOUNT(K,L),K=8,NMOL)
            ENDIF
         ELSE if (noprnt.ge.0) then
            IF (NMOL.LE.7) THEN
               WRITE (IPR,974) L,ZOUT(L),ZOUT(L+1),CTYPE(L),&
                  IPATH(L), PBAR(L),TBAR(L),RHOSUM(L), (AMOUNT(&
                  K,L)/DRAIR,K=1,NMOL),WN2L(L)
            ELSE
               WRITE (IPR,976) L,ZOUT(L),ZOUT(L+1),CTYPE(L),&
                  IPATH(L), PBAR(L),TBAR(L),RHOSUM(L), (AMOUNT(&
                  K,L)/DRAIR,K=1,7),WN2L(L), (AMOUNT(K,L)/     &
                  DRAIR,K=8,NMOL)
            ENDIF
         ENDIF
285   CONTINUE
!
      PWTD = PWTD/WTOT
      TWTD = TWTD/WTOT
      P_h2o = P_h2o/wmt(1)
      T_h2o = T_h2o/wmt(1)

      L = LMAX
      if (noprnt .ge. 0) then
         WRITE (IPR,979) L,ZOUT(1),ZOUT(L+1),P_h2o,T_h2o
         IF (NMOL.LE.7) THEN
            WRITE (IPR,980) (HMOLS(K),K=1,NMOL),COTHER
            WRITE (IPR,982) L,ZOUT(1),ZOUT(L+1),PWTD,TWTD,     &
               SUMRS, (WMT(K),K=1,NMOL),SUMN2
         ELSE
            WRITE (IPR,980) (HMOLS(K),K=1,7),COTHER, (HMOLS(K),&
               K=8,NMOL)
            WRITE (IPR,984) L,ZOUT(1),ZOUT(L+1),PWTD,TWTD,     &
               SUMRS, (WMT(K),K=1,7),SUMN2,(WMT(K),K=8,NMOL)
         ENDIF
      endif
!
      NLAYRS = LMAX
      HT1 = HT1SLT
      HT2 = HT2SLT
!
!-----------------------------------------------------------
! write to ATMOSPHERE file (97) information needed for level to layer ca

      if (ipunch.eq.2) then
         write (97) ibmax,(pbar(l),tbar(l),l=1,ibmax-1)

!*******************************************************************
! MJA 10-17-2012 Try out fix for AJ bug when default atmospheres are used.
! We shouldn't be printing denm, but rather a density interpolated to
! the same grid as PBND and TBND (i.e., on ZOUT)
!*******************************************************************
         IF (MODEL .EQ. 0) THEN ! User-supplied Atm.
            write (97) (pbnd(l),tbnd(l),(denm(k,l),k=1,nmol),l=1, &
               ibmax)
         ELSE !Default Atm.
            !Need to interpolate density (molec/cm3) for each molecule
            !from the model levels in DENM (corresponding to ZMDL)
            !to the boundaries of the output layering (i.e., that
            !correspond to ZOUT)
            do l = 1, ibmax
               ! Find first ZMDL point that is greater than ZOUT(l)
               do i_aj = 1, immax
                  if (ZMDL(i_aj) .GT. ZOUT(l)) exit
               enddo
               ! Calculate weighting
               A_AJ = (ZOUT(l)-ZMDL(i_aj-1))/(ZMDL(i_aj)-ZMDL(i_aj-1))
               do k = 1, nmol
                  DEN_AJ(k,l) = DENM(k,i_aj-1)* &
                     ((DENM(k,i_aj)/DENM(k,i_aj-1))**A_AJ)
               enddo
            enddo
            write (97) (pbnd(l),tbnd(l),(den_aj(k,l),k=1,nmol),l=1, &
               ibmax)
         ENDIF
         close (97)
      endif
!-----------------------------------------------------------

   ENDIF
!
   RETURN
!
!     ERROR MESSAGES
!
290 WRITE (IPR,986) MODEL,ITYPE,NMOL,IBMAX
!
   STOP ' CARD 3.1'
!
!
300 WRITE (IPR,988) (ZBND(I),I=1,IBMAX)
   PRINT 988,(ZBND(I),I=1,IBMAX)
!
   STOP ' ZBND'

301 PRINT 988,(ZTMP(I),I=2,IBMAX)
   STOP ' USER INPUT LEVELS TOO CLOSE - IBMAX'
!
305 WRITE (IPR,989) (PBND(I),I=1,IBMAX)
   PRINT 989,(PBND(I),I=1,IBMAX)
!
   STOP ' PBND'
!
310 WRITE (IPR,990)
!
   STOP ' ERROR - FSCGEO'
!
320 WRITE (IPR,992) AVTRAT,TDIFF1,TDIFF2
!
   STOP ' AVTRAT,TDIFF'
!
900 FORMAT (7I5,I2,1X,I2,4F10.3,A10)
902 FORMAT (' CONTROL CARD 3.1: MODEL AND OPTIONS ')
904 FORMAT (/,10X,'MODEL   = ',I5,/,10X,'ITYPE   = ',I5,/,10X,        &
   &        'IBMAX   = ',I5,/,10X,'n_zero  = ',I5,/,10X,'NOPRNT  = ', &
   &        I5,/,10X,'NMOL    = ',I5,/,10X,'IPUNCH  = ',I5,/,10X,     &
   &        'IFXTYP  = ',I5,/,10X,'MUNITS  = ',I5,/,10X,'RE      = ', &
   &        F10.3,' KM',/,10X,'HSPACE  = ',F10.3,' KM',/,10X,         &
   &        'VBAR    = ',F10.3,' CM-1',                               &
   &                     /,10X,'REF_LAT = ',F10.3, ' DEG')
905 format('$',i5, 10a8)
906 FORMAT (///,' CONTROL CARD 3.1 PARAMETERS WITH DEFAULTS:')
908 FORMAT (///,' HORIZONTAL PATH SELECTED')
910 FORMAT (F10.3,10X,10X,F10.3)
912 FORMAT (///,' CONTROL CARD 3.2:',//,10X,'Z     = ',F10.3,' KM',/, &
   &       10X,'RANGE = ',F10.3,' KM')
914 FORMAT (///' PRESSURE, TEMPERATURE, AND DENSITIES INTERPOLATED',  &
   &        ' FROM THE FOLLOWING ATMOSPHERIC MODEL: ',//,10X,3A8,//,  &
   &        10X, 'Z     = ',F10.3,' KM',/,10X,'P     = ',F10.3,' MB', &
   &        /,10X,'T     = ',F10.3,' K',//,10X,'DENSITIES :',T26,     &
   &        'AIR',(T30,8A10))
916 FORMAT (T63,'(MOL CM-3)',//,T20,1PE10.3,(T30,8E10.3))
918 FORMAT ('0SINGLE LAYER INPUT TO LBLRTM',//,10X,'MODEL = ',3A8,/,  &
   &        10X,'Z     = ',F10.3,' KM',/,10X,'P     = ',F10.3,' MB',  &
   &        /,10X,'T     = ',F10.3,' K',/,10X,'RANGE = ',F10.3,' KM', &
   &        //,10X,'AMOUNTS (MOL CM-2):',T36,'AIR',(T32,8A10))
920 FORMAT (//,T30,1PE10.2,(T30,8E10.2))
922 FORMAT (A4)
924 FORMAT (1X,I1,I3,I5,F10.6,3A8,' * ',F7.3,' KM PATH AT ',F7.3,     &
   &        ' KM ALT')
926 FORMAT (E15.7,F10.4,10X,I5,1X,F7.3,15X,F7.3,/,(1P8E15.7))
928 FORMAT (//,' MULTIPLE SCATTERING TURNED OFF, HMIN = ',F10.6,      &
   &        ' > HMAXMS = ',F10.6,/)
930 FORMAT (///,' SLANT PATH SELECTED, ITYPE = ',I5)
932 FORMAT (5F10.4,I5,5X,F10.4)
933 FORMAT (///' CONTROL CARD 3.2:  SLANT PATH PARAMETERS',//,10X,    &
   &       'H1      = ',F12.6,' MBAR',/,10X,'H2      = ',F12.6,' MBAR'&
   &        /,10X,'ANGLE   = ',F12.6,' DEG',/,10X,'RANGE   = ',F12.6, &
   &        ' KM',/,10X,'BETA    = ',F12.6,' DEG',/,10X,'LEN     = ', &
   &        I10)
934 FORMAT (///' CONTROL CARD 3.2:  SLANT PATH PARAMETERS',//,10X,    &
   &        'H1      = ',F12.6,' KM',/,10X,'H2      = ',F12.6,' KM',  &
   &        /,10X,'ANGLE   = ',F12.6,' DEG',/,10X,'RANGE   = ',F12.6, &
   &        ' KM',/,10X,'BETA    = ',F12.6,' DEG',/,10X,'LEN     = ', &
   &        I10)
936 FORMAT (5F10.3)
938 FORMAT (///,' AUTOLAYERING SELECTED',//,10X,'AVTRAT    = ',F8.2,  &
   &        /,10X,'TDIFF1    = ',F8.2,/,10X,'TDIFF2    = ',F8.2,/,    &
   &        10X,'ALTD1     = ',F8.2,/,10X,'ALTD2     = ',F8.2)
940 FORMAT (8F10.3)
942 FORMAT (///,' USER DEFINED BOUNDARIES FOR LBLRTM LAYERS',/,10X,   &
   &        'I',4X,'Z (KM)',//,(10X,I4,F15.8))
943 FORMAT (///,' USER DEFINED BOUNDARIES FOR LBLRTM LAYERS',/,10X,  &
   &        'I',4X,'P (MB)',//,(10X,I4,F15.8))
944 FORMAT (' ERROR IN USER INPUT BOUNDARIES ')
946 FORMAT (' BOUNDARIES ARE OUTSIDE THE RANGE OF THE ATMOSPHERE',/,  &
   &        ' BOUNDARY = ',F10.2,' ATMOSPHERE =',F10.2,/,             &
   &        ' RESET BOUNDARY GT THAN ATMOSPHERE')
948 FORMAT ('1ATMOSPHERIC PROFILE SELECTED IS: M = ',I3,5X,3A8)
950 FORMAT (/,T4,'I',T13,'Z',T24,'P',T38,'T',T46,'REFRACT',T73,       &
   &        'DENSITY  (MOLS CM-3)',/,T46,'INDEX-1',/,T12,'(KM)',T23,  &
   &        '(MB)',T37,'(K)',T46,'*1.0E6',T63,'AIR',(T68,8(6X,A9)))
951 FORMAT (/,T4,'I',T13,'Z',T24,'P',T38,'T',T46,'REFRACT',T55,       &
   &        'DENSITY',T70,'MIXING RATIO (BASED UPON DRY AIR) (ppmv)',/&
   &        T46,'INDEX-1',T56,                                        &
   &        '(MOL CM-3)'/,T12,                                        &
   &        '(KM)',T23,                                               &
   &        '(MB)',T37,'(K)',T46,'*1.0E6',T61,'AIR',(T68,8(6X,A9)))
952 FORMAT (/)
954 FORMAT (I4,F11.5,F15.8,F11.5,6P,F11.5,1P,E15.7,(T68,1P,8E15.7))
956 FORMAT (///,' HALFWIDTH INFORMATION ON THE USER SUPPLIED ',       &
   &        'LBLRTM BOUNDARIES',/,' THE FOLLOWING VALUES ARE ',       &
   &        'ASSUMED:')
958 FORMAT (10X,'ALZERO    = ',F9.3,' CM-1 = AVERAGE LORENTZ WIDTH ', &
   &        'AT STP',/,10X,'AVMWT     = ',F8.2,                       &
   &        '       = AVERAGE  MOLECULAR WEIGHT',/,10X,               &
   &        'VBAR      = ',F8.2,'  CM-1 = AVERAGE WAVENUMBER',///,    &
   &        T5,'I',T12,'Z',T22,'P',T32,'T',T39,'LORENTZ',T49,         &
   &        'DOPPLER',T61,'ZETA',T70,'VOIGT',T80,'VOIGT',T90,         &
   &        'TEMP',/,T11,'(KM)',T21,'(MB)',T31,'(K)',T40,'(CM-1)',    &
   &        T50,'(CM-1)',T70,'(CM-1)',T80,'RATIO',T90,'DIFF (K)',/)
960 FORMAT (I5,F10.3,F12.5,F9.2,F9.5,F10.5,F10.3,F10.5,F10.2,F10.1)
962 FORMAT ('1INTEGRATED ABSORBER AMOUNTS BY LAYER',///,T5,           &
   &        'I  LAYER BOUNDARIES',T55,'INTEGRATED AMOUNTS ',          &
   &        '(MOL CM-2)',/,T11,'FROM',T22,'TO',T29,'AIR',T36,         &
   &        8(1X,A8,1X),/,T11,'(KM)',T21,'(KM)',(T37,8A10))
964 FORMAT (I5,2F10.3,1P,E10.3,(T36,1P,8E10.3))
966 FORMAT ('0TOTAL',F9.3,F10.3,1P,E10.3,(T36,1P,8E10.3))
968 FORMAT ('1 SUMMARY OF THE GEOMETRY CALCULATION',//,10X,           &
   &        'MODEL   = ',4X,3A8,/10X,'H1      = ',F12.6,' KM',/,10X,  &
   &        'H2      = ',F12.6,' KM',/,10X,'ANGLE   = ',F12.6,' DEG', &
   &        /,10X,'RANGE   = ',F12.6,' KM',/,10X,'BETA    = ',F12.6,  &
   &        ' DEG',/,10X,'PHI     = ',F12.6,' DEG',/,10X,             &
   &        'HMIN    = ',F12.6,' KM',/,10X,'BENDING = ',F12.6,' DEG', &
   &        /,10X,'LEN     = ',I10,/,10X,'AIRMAS  = ',G12.6,          &
   &        'RELATIVE TO A VERTICAL PATH , GROUND TO SPACE')
970 FORMAT ('0FINAL SET OF LAYERS FOR INPUT TO LBLRTM',/,             &
   &        ' A LAYER AMOUNT MAY BE SET TO ZERO IF THE CUMULATIVE ',  &
   &        'AMOUNT FOR THAT LAYER AND ABOVE IS LESS THAN 0.1 ',      &
   &        'PERCENT',/,' OF THE TOTAL AMOUNT. THIS IS DONE ONLY ',   &
   &        'FOR THE FOLLOWING CASES',/,5X,'1.  IEMIT = 0 ',          &
   &        '(TRANSMITTANCE)',/,5X,'2.  IEMIT = 1 (RADIANCE) AND ',   &
   &        'IPATH = 3 (PATH LOOKING UP)',/,' O2 IS NOT INCLUDED',/,  &
   &        ' IF THE AMOUNTS FOR ALL THE MOLECULES BUT O2 ARE ',      &
   &        'ZEROED, THE REMAINING LAYERS ARE ELIMINATED',///,T13,    &
   &        'LAYER',T23,'I',T25,'I',/,T4,'L',T10,'BOUNDARIES',T23,    &
   &        'T',T25,'P',T31,'PBAR',T40,'TBAR',T70,                    &
   &        'INTEGRATED AMOUNTS (MOLS CM-2)',/,T9,'FROM',T18,         &
   &        'TO',T23,'Y',T25,'T',/,T9,'(KM)',T17,'(KM)',T23,'P',T25,  &
   &        'H',T31,'(MB)',T41,'(K)',T53,'AIR',(T59,8(6X,A9)))
972 FORMAT (1X,I1,I3,I5,F10.6,2A8,' H1=',F8.2,' H2=',F8.2,            &
   &        ' ANG=',F8.3,' LEN=',I2)
973 FORMAT  (///,'1',3X,'------------------------------------',/,     &
   &        T13,'LAYER',T23,'I',T25,'I',/,T4,'L',T10,'BOUNDARIES',    &
   &        T23,'T',T25,'P',T31,'PBAR',T40,'TBAR',                    &
   &        T68,'MOLECULAR MIXING RATIOS BY LAYER',/,T9,'FROM',       &
   &        T18,'TO',T23,'Y',T25,'T',/,T9,'(KM)',T17,'(KM)',T23,'P',  &
   &        T25,'H',T31,'(MB)',T41,'(K)',T53,'AIR',(T59,8(6X,A9)))
974 FORMAT ('0',I3,2F8.3,A3,I2,F11.5,F8.2,1X,1P9E15.7)
976 FORMAT ('0',I3,2F8.3,A3,I2,F11.5,F8.2,1X,1P9E15.7,/,              &
   &            (60X,1P8E15.7))
978 FORMAT (1P8E15.7)
979 FORMAT ('0',/,'0',T4,'L  PATH BOUNDARIES',T28,'P_h2o',T37,'T_h2o' &
   &           ,/,'0',I3,2F8.3,2X,F11.5,F8.2)
980 FORMAT ('0',/,'0',T4,'L  PATH BOUNDARIES',T28,'PBAR',T37,'TBAR',  &
   &        T65,'ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH',/,T9,  &
   &        'FROM',T18,'TO',/,T9,'(KM)',T17,'(KM)',T28,'(MB)',T38,    &
   &        '(K)',T47,'AIR',(T54,8(1X,A9)))
982 FORMAT ('0',I3,2F8.3,2X,F11.5,F8.2,1X,1P,9E10.3)
984 FORMAT ('0',I3,2F8.3,2X,F11.5,F8.2,1X,1P,9E10.3,/,(52X,8E10.3))
986 FORMAT (///,' ERROR IN INPUT, CONTROL CARD 3.1: ONE OF THE ',     &
   &        'PARAMETERS MODEL, ITYPE, NMOL IS OUT OF RANGE',//,10X,   &
   &        'MODEL   = ',I5,/,10X,'ITYPE   = ',I5,/,10X,'NMOL    = ', &
   &        I5,10X,' IBMAX =',I5)
988 FORMAT (///,' ERROR: BOUNDARY ALTITUDES FOR LBLRTM LAYERS ',      &
   &        'ARE NEGATIVE OR NOT IN ASCENDING ORDER',//,5X,' ZBND ',  &
   &        /,(10F10.4))
989 FORMAT (///,' ERROR: BOUNDARY PRESSURES FOR LBLRTM LAYERS ',      &
   &        'ARE POSITIVE OR NOT IN DESCENDING ORDER',//,5X,' PBND ', &
   &        /,(10F10.4))
990 FORMAT ('0ERROR FLAG RETURNED FROM FSCGEO:  AN ERROR OCCURED ',   &
   &        'IN PROCESSING THE SLANT PATH PARAMETERS',/,'0PROGRAM ',  &
   &        'STOP')
992 FORMAT (///,' ERROR: EITHER AVTRAT.LE.1.0 OR TDIFF.LE.0',/,       &
   &        '0PROGRAM STOP  -  AVTRAT = ',E12.6,' TDIFF1 = ',F10.4,   &
   &        ' TDIFF2 = ',F10.4)
1000 FORMAT ('*** WARNING: Zeroing molecule #',i2.2,' amount ',        &
   &        'in layer #',i3.3)
1100 format (i3)
1101 format (a3)
!
end subroutine ATMPTH
!
!     ----------------------------------------------------------------
!
BLOCK DATA ATMCON
!
!     *****************************************************************
!     THIS SUBROUTINE INITIALIZES THE CONSTANTS USED IN THE
!     PROGRAM. CONSTANTS RELATING TO THE ATMOSPHERIC PROFILES ARE STORE
!     IN BLOCK DATA MLATMB.
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
!
   CHARACTER*8      HMOLS
!
   COMMON /HMOLS/ HMOLS(MXMOL),JUNIT(MXMOL),WMOL(MXMOL),JUNITP,      &
   &               JUNITT
   COMMON /HMOLC/ HMOLC(MXMOL)
   COMMON /FIXITYL/ IFXTYP
   CHARACTER*8 HMOLC
!     IFXTYP is the flag for fixing the value of ITYL
   DATA IFXTYP /0/
!
!     IBDIM IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
!     IOUTDM IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
!     IMDIM IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
!         STORED IN ZMDL (INPUT)
!     IPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH OBTAINE
!         BY MERGING ZMDL AND ZOUT
!     KDIM IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
!
   DATA KMXNOM / 7 /
!
!     DELTAS IS THE NOMINAL SLANT PATH INCREMENT IN KM.
!
   DATA DELTAS / 5.0 /
   DATA PZERO / 1013.25 /,TZERO / 273.15 /
!
!     ALZERO IS THE MEAN LORENTZ HALFWIDTH AT PZERO AND 296.0 K.
!     AVMWT IS THE MEAN MOLECULAR WEIGHT USED TO AUTOMATICALLY
!     GENERATE THE LBLRTM BOUNDARIES IN AUTLAY
!
   DATA ALZERO / 0.04 /,AVMWT / 36.0 /
!
!     ORDER OF MOLECULES H2O(1), CO2(2), O3(3), N2O(4), CO(5), CH4(6),
!         O2(7), NO(8), SO2(9), NO2(10), NH3(11), HNO3(12), OH(13),
!         HF(14 ), HCL(15), HBR(16), HI(17), CLO(18), OCS(19), H2CO(20)
!         HOCL(21), N2(22), HCN(23), CH3CL(24), H2O2(25), C2H2(26),
!         C2H6(27), PH3(28), COF2(29), SF6(30), H2S(31), HCOOH(32),
!         HO2(33), O(34), ClONO2(35), NO+(36), HOBr(37), C2H4(38),
!         CH3OH(39), CH3Br(40), CH3CN(41), CF4(42), C4H2(43), HC3N(44),
!         H2(45), CS(46), SO3(47)
!
   DATA HMOLC / '  H2O   ' , '  CO2   ' , '   O3   ' , '  N2O   ' ,  &
   &             '   CO   ' , '  CH4   ' , '   O2   ' , '   NO   ' ,  &
   &             '  SO2   ' , '  NO2   ' , '  NH3   ' , ' HNO3   ' ,  &
   &             '   OH   ' , '   HF   ' , '  HCL   ' , '  HBR   ' ,  &
   &             '   HI   ' , '  CLO   ' , '  OCS   ' , ' H2CO   ' ,  &
   &             ' HOCL   ' , '   N2   ' , '  HCN   ' , ' CH3CL  ' ,  &
   &             ' H2O2   ' , ' C2H2   ' , ' C2H6   ' , '  PH3   ' ,  &
   &             ' COF2   ' , '  SF6   ' , '  H2S   ' , ' HCOOH  ' ,  &
   &             '  HO2   ' , '   O+   ' , ' ClONO2 ' , '   NO+  ' ,  &
   &             '  HOBr  ' , ' C2H4   ' , ' CH3OH  ' , ' CH3Br  ' ,  &
   &             ' CH3CN  ' , '  CF4   ' , ' C4H2   ' , ' HC3N   ' ,  &
   &             '   H2   ' , '   CS   ' , '  SO3   '/
!
!     MOLECULAR WEIGHTS
!
   DATA AMWT /  18.015 ,  44.010 , 47.998 , 44.01 ,                 &
   &              28.011 ,  16.043 , 31.999 , 30.01 ,                 &
   &              64.06  ,  46.01  , 17.03  , 63.01 ,                 &
   &              17.00  ,  20.01  , 36.46  , 80.92 ,                 &
   &             127.91  ,  51.45  , 60.08  , 30.03 ,                 &
   &              52.46  ,  28.014 , 27.03  , 50.49 ,                 &
   &              34.01  ,  26.03  , 30.07  , 34.00 ,                 &
   &              66.01  , 146.05  , 34.08  , 46.03 ,                 &
   &              33.00  ,  15.99  , 98.    , 30.00 ,                 &
   &              97.    ,  28.05  , 32.04  , 94.94 ,                 &
   &              41.05  ,  88.0043, 50.06  , 51.05 ,                 &
   &               2.016 ,  44.08  , 80.066 /
!     approx;
!     MJA, 10/04/2011 fixed ethylene molecular weight (Molecule #38)
end block data ATMCON
!
!     ----------------------------------------------------------------
!
BLOCK DATA MLATMB
!
!     *****************************************************************
!     THIS SUBROUTINE INITIALIZES THE 6 BUILT-IN ATMOSPHERIC PROFILES
!     (FROM 'OPTICAL PROPERTIES OF THE ATMOSPHERE, THIRD EDITION'
!     AFCRL-72-0497 (AD 753 075), 'U.S. STANDARD ATMOSPHERE 1976' AND
!     'SUPPLEMENTS 1966'), PLUS COLLECTED CONSTITUENT PROFILES (REF)
!     AND SETS OTHER CONSTANTS RELATED TO THE ATMOSPHERIC PROFILES
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=47,MXTRAC=22)
   PARAMETER (MXZ50=MXZMD-50)
!
   COMMON /MLATM/ ALT(MXZMD),                                        &
   &     P1(MXZMD),P2(MXZMD),P3(MXZMD),P4(MXZMD),P5(MXZMD),P6(MXZMD), &
   &     T1(MXZMD),T2(MXZMD),T3(MXZMD),T4(MXZMD),T5(MXZMD),T6(MXZMD), &
   &     AMOL11(MXZMD),AMOL12(MXZMD),AMOL13(MXZMD),AMOL14(MXZMD),     &
   &     AMOL15(MXZMD),AMOL16(MXZMD),AMOL17(MXZMD),AMOL18(MXZMD),     &
   &     AMOL21(MXZMD),AMOL22(MXZMD),AMOL23(MXZMD),AMOL24(MXZMD),     &
   &     AMOL25(MXZMD),AMOL26(MXZMD),AMOL27(MXZMD),AMOL28(MXZMD),     &
   &     AMOL31(MXZMD),AMOL32(MXZMD),AMOL33(MXZMD),AMOL34(MXZMD),     &
   &     AMOL35(MXZMD),AMOL36(MXZMD),AMOL37(MXZMD),AMOL38(MXZMD),     &
   &     AMOL41(MXZMD),AMOL42(MXZMD),AMOL43(MXZMD),AMOL44(MXZMD),     &
   &     AMOL45(MXZMD),AMOL46(MXZMD),AMOL47(MXZMD),AMOL48(MXZMD),     &
   &     AMOL51(MXZMD),AMOL52(MXZMD),AMOL53(MXZMD),AMOL54(MXZMD),     &
   &     AMOL55(MXZMD),AMOL56(MXZMD),AMOL57(MXZMD),AMOL58(MXZMD),     &
   &     AMOL61(MXZMD),AMOL62(MXZMD),AMOL63(MXZMD),AMOL64(MXZMD),     &
   &     AMOL65(MXZMD),AMOL66(MXZMD),AMOL67(MXZMD),AMOL68(MXZMD),     &
   &     ZST(MXZMD),PST(MXZMD),TST(MXZMD),AMOLS(MXZMD,MXMOL)
   COMMON /MLATMC/ ATMNAM(6)
   CHARACTER*24 ATMNAM
!
!     COMMON /TRAC/ TRAC(MXZMD,MXTRAC)
!
   COMMON /TRAC/ ANO(MXZMD),SO2(MXZMD),ANO2(MXZMD),ANH3(MXZMD),      &
   &              HNO3(MXZMD),OH(MXZMD),HF(MXZMD),HCL(MXZMD),         &
   &              HBR(MXZMD),HI(MXZMD),CLO(MXZMD),OCS(MXZMD),         &
   &              H2CO(MXZMD),HOCL(MXZMD),AN2(MXZMD),HCN(MXZMD),      &
   &              CH3CL(MXZMD),H2O2(MXZMD),C2H2(MXZMD),C2H6(MXZMD),   &
   &              PH3(MXZMD), COF2(MXZMD), SF6(MXZMD), H2S(MXZMD),    &
   &              HCOOH(MXZMD), HO2(MXZMD), O(MXZMD), CLONO2(MXZMD),  &
   &           NOPLUS(MXZMD), HOBr(MXZMD), C2H4(MXZMD), CH3OH(MXZMD), &
   &              CH3BR(MXZMD), CH3CN(MXZMD), CF4(MXZMD), C4H2(MXZMD),&
   &              HC3N(MXZMD), H2(MXZMD), CS(MXZMD), SO3(MXZMD),      &
   &              TDUM(MXZMD)
!
   DATA ATMNAM(1) / 'TROPICAL                '/
   DATA ATMNAM(2) / 'MIDLATITUDE SUMMER      '/
   DATA ATMNAM(3) / 'MIDLATITUDE WINTER      '/
   DATA ATMNAM(4) / 'SUBARCTIC SUMMER        '/
   DATA ATMNAM(5) / 'SUBARCTIC WINTER        '/
   DATA ATMNAM(6) / 'U. S. STANDARD,  1976   '/
!
!     DATA ALT (KM) /
!
   DATA ALT /       0.0,       1.0,       2.0,       3.0,       4.0, &
   &                 5.0,       6.0,       7.0,       8.0,       9.0, &
   &                10.0,      11.0,      12.0,      13.0,      14.0, &
   &                15.0,      16.0,      17.0,      18.0,      19.0, &
   &                20.0,      21.0,      22.0,      23.0,      24.0, &
   &                25.0,      27.5,      30.0,      32.5,      35.0, &
   &                37.5,      40.0,      42.5,      45.0,      47.5, &
   &                50.0,      55.0,      60.0,      65.0,      70.0, &
   &                75.0,      80.0,      85.0,      90.0,      95.0, &
   &               100.0,     105.0,     110.0,     115.0,     120.0, &
   &               MXZ50*0.0 /
!
!     DATA PRESSURE /
!
   DATA P1 /  1.013E+03, 9.040E+02, 8.050E+02, 7.150E+02, 6.330E+02, &
   &           5.590E+02, 4.920E+02, 4.320E+02, 3.780E+02, 3.290E+02, &
   &           2.860E+02, 2.470E+02, 2.130E+02, 1.820E+02, 1.560E+02, &
   &           1.320E+02, 1.110E+02, 9.370E+01, 7.890E+01, 6.660E+01, &
   &           5.650E+01, 4.800E+01, 4.090E+01, 3.500E+01, 3.000E+01, &
   &           2.570E+01, 1.763E+01, 1.220E+01, 8.520E+00, 6.000E+00, &
   &           4.260E+00, 3.050E+00, 2.200E+00, 1.590E+00, 1.160E+00, &
   &           8.540E-01, 4.560E-01, 2.390E-01, 1.210E-01, 5.800E-02, &
   &           2.600E-02, 1.100E-02, 4.400E-03, 1.720E-03, 6.880E-04, &
   &           2.890E-04, 1.300E-04, 6.470E-05, 3.600E-05, 2.250E-05, &
   &           MXZ50*0.0 /
!
   DATA P2 /  1.013E+03, 9.020E+02, 8.020E+02, 7.100E+02, 6.280E+02, &
   &           5.540E+02, 4.870E+02, 4.260E+02, 3.720E+02, 3.240E+02, &
   &           2.810E+02, 2.430E+02, 2.090E+02, 1.790E+02, 1.530E+02, &
   &           1.300E+02, 1.110E+02, 9.500E+01, 8.120E+01, 6.950E+01, &
   &           5.950E+01, 5.100E+01, 4.370E+01, 3.760E+01, 3.220E+01, &
   &           2.770E+01, 1.907E+01, 1.320E+01, 9.300E+00, 6.520E+00, &
   &           4.640E+00, 3.330E+00, 2.410E+00, 1.760E+00, 1.290E+00, &
   &           9.510E-01, 5.150E-01, 2.720E-01, 1.390E-01, 6.700E-02, &
   &           3.000E-02, 1.200E-02, 4.480E-03, 1.640E-03, 6.250E-04, &
   &           2.580E-04, 1.170E-04, 6.110E-05, 3.560E-05, 2.270E-05, &
   &           MXZ50*0.0 /
!
   DATA P3 /  1.018E+03, 8.973E+02, 7.897E+02, 6.938E+02, 6.081E+02, &
   &           5.313E+02, 4.627E+02, 4.016E+02, 3.473E+02, 2.993E+02, &
   &           2.568E+02, 2.199E+02, 1.882E+02, 1.611E+02, 1.378E+02, &
   &           1.178E+02, 1.007E+02, 8.610E+01, 7.360E+01, 6.280E+01, &
   &           5.370E+01, 4.580E+01, 3.910E+01, 3.340E+01, 2.860E+01, &
   &           2.440E+01, 1.646E+01, 1.110E+01, 7.560E+00, 5.180E+00, &
   &           3.600E+00, 2.530E+00, 1.800E+00, 1.290E+00, 9.400E-01, &
   &           6.830E-01, 3.620E-01, 1.880E-01, 9.500E-02, 4.700E-02, &
   &           2.220E-02, 1.030E-02, 4.560E-03, 1.980E-03, 8.770E-04, &
   &           4.074E-04, 2.000E-04, 1.057E-04, 5.980E-05, 3.600E-05, &
   &           MXZ50*0.0 /
!
   DATA P4 /  1.010E+03, 8.960E+02, 7.929E+02, 7.000E+02, 6.160E+02, &
   &           5.410E+02, 4.740E+02, 4.130E+02, 3.590E+02, 3.108E+02, &
   &           2.677E+02, 2.300E+02, 1.977E+02, 1.700E+02, 1.460E+02, &
   &           1.260E+02, 1.080E+02, 9.280E+01, 7.980E+01, 6.860E+01, &
   &           5.900E+01, 5.070E+01, 4.360E+01, 3.750E+01, 3.228E+01, &
   &           2.780E+01, 1.923E+01, 1.340E+01, 9.400E+00, 6.610E+00, &
   &           4.720E+00, 3.400E+00, 2.480E+00, 1.820E+00, 1.340E+00, &
   &           9.870E-01, 5.370E-01, 2.880E-01, 1.470E-01, 7.100E-02, &
   &           3.200E-02, 1.250E-02, 4.510E-03, 1.610E-03, 6.060E-04, &
   &           2.480E-04, 1.130E-04, 6.000E-05, 3.540E-05, 2.260E-05, &
   &           MXZ50*0.0 /
!
   DATA P5 /  1.013E+03, 8.878E+02, 7.775E+02, 6.798E+02, 5.932E+02, &
   &           5.158E+02, 4.467E+02, 3.853E+02, 3.308E+02, 2.829E+02, &
   &           2.418E+02, 2.067E+02, 1.766E+02, 1.510E+02, 1.291E+02, &
   &           1.103E+02, 9.431E+01, 8.058E+01, 6.882E+01, 5.875E+01, &
   &           5.014E+01, 4.277E+01, 3.647E+01, 3.109E+01, 2.649E+01, &
   &           2.256E+01, 1.513E+01, 1.020E+01, 6.910E+00, 4.701E+00, &
   &           3.230E+00, 2.243E+00, 1.570E+00, 1.113E+00, 7.900E-01, &
   &           5.719E-01, 2.990E-01, 1.550E-01, 7.900E-02, 4.000E-02, &
   &           2.000E-02, 9.660E-03, 4.500E-03, 2.022E-03, 9.070E-04, &
   &           4.230E-04, 2.070E-04, 1.080E-04, 6.000E-05, 3.590E-05, &
   &           MXZ50*0.0 /
!
   DATA P6 /  1.013E+03, 8.988E+02, 7.950E+02, 7.012E+02, 6.166E+02, &
   &           5.405E+02, 4.722E+02, 4.111E+02, 3.565E+02, 3.080E+02, &
   &           2.650E+02, 2.270E+02, 1.940E+02, 1.658E+02, 1.417E+02, &
   &           1.211E+02, 1.035E+02, 8.850E+01, 7.565E+01, 6.467E+01, &
   &           5.529E+01, 4.729E+01, 4.047E+01, 3.467E+01, 2.972E+01, &
   &           2.549E+01, 1.743E+01, 1.197E+01, 8.258E+00, 5.746E+00, &
   &           4.041E+00, 2.871E+00, 2.060E+00, 1.491E+00, 1.090E+00, &
   &           7.978E-01, 4.250E-01, 2.190E-01, 1.090E-01, 5.220E-02, &
   &           2.400E-02, 1.050E-02, 4.460E-03, 1.840E-03, 7.600E-04, &
   &           3.200E-04, 1.450E-04, 7.100E-05, 4.010E-05, 2.540E-05, &
   &           MXZ50*0.0 /
!
!     DATA TEMPERATURE /
!
   DATA T1 /     299.70,    293.70,    287.70,    283.70,    277.00, &
   &              270.30,    263.60,    257.00,    250.30,    243.60, &
   &              237.00,    230.10,    223.60,    217.00,    210.30, &
   &              203.70,    197.00,    194.80,    198.80,    202.70, &
   &              206.70,    210.70,    214.60,    217.00,    219.20, &
   &              221.40,    227.00,    232.30,    237.70,    243.10, &
   &              248.50,    254.00,    259.40,    264.80,    269.60, &
   &              270.20,    263.40,    253.10,    236.00,    218.90, &
   &              201.80,    184.80,    177.10,    177.00,    184.30, &
   &              190.70,    212.00,    241.60,    299.70,    380.00, &
   &              MXZ50*0.0 /
!
   DATA T2 /     294.20,    289.70,    285.20,    279.20,    273.20, &
   &              267.20,    261.20,    254.70,    248.20,    241.70, &
   &              235.30,    228.80,    222.30,    215.80,    215.70, &
   &              215.70,    215.70,    215.70,    216.80,    217.90, &
   &              219.20,    220.40,    221.60,    222.80,    223.90, &
   &              225.10,    228.45,    233.70,    239.00,    245.20, &
   &              251.30,    257.50,    263.70,    269.90,    275.20, &
   &              275.70,    269.30,    257.10,    240.10,    218.10, &
   &              196.10,    174.10,    165.10,    165.00,    178.30, &
   &              190.50,    222.20,    262.40,    316.80,    380.00, &
   &              MXZ50*0.0 /
!
   DATA T3 /     272.20,    268.70,    265.20,    261.70,    255.70, &
   &              249.70,    243.70,    237.70,    231.70,    225.70, &
   &              219.70,    219.20,    218.70,    218.20,    217.70, &
   &              217.20,    216.70,    216.20,    215.70,    215.20, &
   &              215.20,    215.20,    215.20,    215.20,    215.20, &
   &              215.20,    215.50,    217.40,    220.40,    227.90, &
   &              235.50,    243.20,    250.80,    258.50,    265.10, &
   &              265.70,    260.60,    250.80,    240.90,    230.70, &
   &              220.40,    210.10,    199.80,    199.50,    208.30, &
   &              218.60,    237.10,    259.50,    293.00,    333.00, &
   &              MXZ50*0.0 /
!
   DATA T4 /     287.20,    281.70,    276.30,    270.90,    265.50, &
   &              260.10,    253.10,    246.10,    239.20,    232.20, &
   &              225.20,    225.20,    225.20,    225.20,    225.20, &
   &              225.20,    225.20,    225.20,    225.20,    225.20, &
   &              225.20,    225.20,    225.20,    225.20,    226.60, &
   &              228.10,    231.00,    235.10,    240.00,    247.20, &
   &              254.60,    262.10,    269.50,    273.60,    276.20, &
   &              277.20,    274.00,    262.70,    239.70,    216.60, &
   &              193.60,    170.60,    161.70,    161.60,    176.80, &
   &              190.40,    226.00,    270.10,    322.70,    380.00, &
   &              MXZ50*0.0 /
!
   DATA T5 /     257.20,    259.10,    255.90,    252.70,    247.70, &
   &              240.90,    234.10,    227.30,    220.60,    217.20, &
   &              217.20,    217.20,    217.20,    217.20,    217.20, &
   &              217.20,    216.60,    216.00,    215.40,    214.80, &
   &              214.20,    213.60,    213.00,    212.40,    211.80, &
   &              211.20,    213.60,    216.00,    218.50,    222.30, &
   &              228.50,    234.70,    240.80,    247.00,    253.20, &
   &              259.30,    259.10,    250.90,    248.40,    245.40, &
   &              234.70,    223.90,    213.10,    202.30,    211.00, &
   &              218.50,    234.00,    252.60,    288.50,    333.00, &
   &              MXZ50*0.0 /
!
   DATA T6 /     288.20,    281.70,    275.20,    268.70,    262.20, &
   &              255.70,    249.20,    242.70,    236.20,    229.70, &
   &              223.30,    216.80,    216.70,    216.70,    216.70, &
   &              216.70,    216.70,    216.70,    216.70,    216.70, &
   &              216.70,    217.60,    218.60,    219.60,    220.60, &
   &              221.60,    224.00,    226.50,    229.60,    236.50, &
   &              243.40,    250.40,    257.30,    264.20,    270.60, &
   &              270.70,    260.80,    247.00,    233.30,    219.60, &
   &              208.40,    198.60,    188.90,    186.90,    188.40, &
   &              195.10,    208.80,    240.00,    300.00,    360.00, &
   &              MXZ50*0.0 /
!
!     DATA  H2O /
!
   DATA AMOL11 /                                                     &
   &           2.593E+04, 1.949E+04, 1.534E+04, 8.600E+03, 4.441E+03, &
   &           3.346E+03, 2.101E+03, 1.289E+03, 7.637E+02, 4.098E+02, &
   &           1.912E+02, 7.306E+01, 2.905E+01, 9.900E+00, 6.220E+00, &
   &           4.000E+00, 3.000E+00, 2.900E+00, 2.750E+00, 2.600E+00, &
   &           2.600E+00, 2.650E+00, 2.800E+00, 2.900E+00, 3.200E+00, &
   &           3.250E+00, 3.600E+00, 4.000E+00, 4.300E+00, 4.600E+00, &
   &           4.900E+00, 5.200E+00, 5.500E+00, 5.700E+00, 5.900E+00, &
   &           6.000E+00, 6.000E+00, 6.000E+00, 5.400E+00, 4.500E+00, &
   &           3.300E+00, 2.100E+00, 1.300E+00, 8.500E-01, 5.400E-01, &
   &           4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01, &
   &           MXZ50*0.0 /
!
   DATA AMOL21 /                                                     &
   &           1.876E+04, 1.378E+04, 9.680E+03, 5.984E+03, 3.813E+03, &
   &           2.225E+03, 1.510E+03, 1.020E+03, 6.464E+02, 4.129E+02, &
   &           2.472E+02, 9.556E+01, 2.944E+01, 8.000E+00, 5.000E+00, &
   &           3.400E+00, 3.300E+00, 3.200E+00, 3.150E+00, 3.200E+00, &
   &           3.300E+00, 3.450E+00, 3.600E+00, 3.850E+00, 4.000E+00, &
   &           4.200E+00, 4.450E+00, 4.700E+00, 4.850E+00, 4.950E+00, &
   &           5.000E+00, 5.100E+00, 5.300E+00, 5.450E+00, 5.500E+00, &
   &           5.500E+00, 5.350E+00, 5.000E+00, 4.400E+00, 3.700E+00, &
   &           2.950E+00, 2.100E+00, 1.330E+00, 8.500E-01, 5.400E-01, &
   &           4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01, &
   &           MXZ50*0.0 /
!
   DATA AMOL31 /                                                     &
   &           4.316E+03, 3.454E+03, 2.788E+03, 2.088E+03, 1.280E+03, &
   &           8.241E+02, 5.103E+02, 2.321E+02, 1.077E+02, 5.566E+01, &
   &           2.960E+01, 1.000E+01, 6.000E+00, 5.000E+00, 4.800E+00, &
   &           4.700E+00, 4.600E+00, 4.500E+00, 4.500E+00, 4.500E+00, &
   &           4.500E+00, 4.500E+00, 4.530E+00, 4.550E+00, 4.600E+00, &
   &           4.650E+00, 4.700E+00, 4.750E+00, 4.800E+00, 4.850E+00, &
   &           4.900E+00, 4.950E+00, 5.000E+00, 5.000E+00, 5.000E+00, &
   &           4.950E+00, 4.850E+00, 4.500E+00, 4.000E+00, 3.300E+00, &
   &           2.700E+00, 2.000E+00, 1.330E+00, 8.500E-01, 5.400E-01, &
   &           4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01, &
   &           MXZ50*0.0 /
!
   DATA AMOL41 /                                                     &
   &           1.194E+04, 8.701E+03, 6.750E+03, 4.820E+03, 3.380E+03, &
   &           2.218E+03, 1.330E+03, 7.971E+02, 3.996E+02, 1.300E+02, &
   &           4.240E+01, 1.330E+01, 6.000E+00, 4.450E+00, 4.000E+00, &
   &           4.000E+00, 4.000E+00, 4.050E+00, 4.300E+00, 4.500E+00, &
   &           4.600E+00, 4.700E+00, 4.800E+00, 4.830E+00, 4.850E+00, &
   &           4.900E+00, 4.950E+00, 5.000E+00, 5.000E+00, 5.000E+00, &
   &           5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00, &
   &           4.950E+00, 4.850E+00, 4.500E+00, 4.000E+00, 3.300E+00, &
   &           2.700E+00, 2.000E+00, 1.330E+00, 8.500E-01, 5.400E-01, &
   &           4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01, &
   &           MXZ50*0.0 /
!
   DATA AMOL51 /                                                     &
   &           1.405E+03, 1.615E+03, 1.427E+03, 1.166E+03, 7.898E+02, &
   &           4.309E+02, 2.369E+02, 1.470E+02, 3.384E+01, 2.976E+01, &
   &           2.000E+01, 1.000E+01, 6.000E+00, 4.450E+00, 4.500E+00, &
   &           4.550E+00, 4.600E+00, 4.650E+00, 4.700E+00, 4.750E+00, &
   &           4.800E+00, 4.850E+00, 4.900E+00, 4.950E+00, 5.000E+00, &
   &           5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00, &
   &           5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00, &
   &           4.950E+00, 4.850E+00, 4.500E+00, 4.000E+00, 3.300E+00, &
   &           2.700E+00, 2.000E+00, 1.330E+00, 8.500E-01, 5.400E-01, &
   &           4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01, &
   &           MXZ50*0.0 /
!
   DATA AMOL61 /                                                     &
   &           7.745E+03, 6.071E+03, 4.631E+03, 3.182E+03, 2.158E+03, &
   &           1.397E+03, 9.254E+02, 5.720E+02, 3.667E+02, 1.583E+02, &
   &           6.996E+01, 3.613E+01, 1.906E+01, 1.085E+01, 5.927E+00, &
   &           5.000E+00, 3.950E+00, 3.850E+00, 3.825E+00, 3.850E+00, &
   &           3.900E+00, 3.975E+00, 4.065E+00, 4.200E+00, 4.300E+00, &
   &           4.425E+00, 4.575E+00, 4.725E+00, 4.825E+00, 4.900E+00, &
   &           4.950E+00, 5.025E+00, 5.150E+00, 5.225E+00, 5.250E+00, &
   &           5.225E+00, 5.100E+00, 4.750E+00, 4.200E+00, 3.500E+00, &
   &           2.825E+00, 2.050E+00, 1.330E+00, 8.500E-01, 5.400E-01, &
   &           4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01, &
   &           MXZ50*0.0 /
!
!     DATA CO2 /
!
   DATA AMOL12 /                                                     &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02, &
   &           1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01, &
   &           MXZ50*0.0 /
!
   DATA AMOL22 /                                                     &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02, &
   &           1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01, &
   &           MXZ50*0.0 /
!
   DATA AMOL32 /                                                     &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02, &
   &           1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01, &
   &           MXZ50*0.0 /
!
   DATA AMOL42 /                                                     &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02, &
   &           1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01, &
   &           MXZ50*0.0 /
!
   DATA AMOL52 /                                                     &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02, &
   &           1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01, &
   &           MXZ50*0.0 /
!
   DATA AMOL62 /                                                     &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, &
   &           3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02, &
   &           1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01, &
   &           MXZ50*0.0 /
!
!     DATA OZONE /
!
   DATA AMOL13 /                                                     &
   &           2.869E-02, 3.150E-02, 3.342E-02, 3.504E-02, 3.561E-02, &
   &           3.767E-02, 3.989E-02, 4.223E-02, 4.471E-02, 5.000E-02, &
   &           5.595E-02, 6.613E-02, 7.815E-02, 9.289E-02, 1.050E-01, &
   &           1.256E-01, 1.444E-01, 2.500E-01, 5.000E-01, 9.500E-01, &
   &           1.400E+00, 1.800E+00, 2.400E+00, 3.400E+00, 4.300E+00, &
   &           5.400E+00, 7.800E+00, 9.300E+00, 9.850E+00, 9.700E+00, &
   &           8.800E+00, 7.500E+00, 5.900E+00, 4.500E+00, 3.450E+00, &
   &           2.800E+00, 1.800E+00, 1.100E+00, 6.500E-01, 3.000E-01, &
   &           1.800E-01, 3.300E-01, 5.000E-01, 5.200E-01, 5.000E-01, &
   &           4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04, &
   &           MXZ50*0.0 /
!
   DATA AMOL23 /                                                     &
   &           3.017E-02, 3.337E-02, 3.694E-02, 4.222E-02, 4.821E-02, &
   &           5.512E-02, 6.408E-02, 7.764E-02, 9.126E-02, 1.111E-01, &
   &           1.304E-01, 1.793E-01, 2.230E-01, 3.000E-01, 4.400E-01, &
   &           5.000E-01, 6.000E-01, 7.000E-01, 1.000E+00, 1.500E+00, &
   &           2.000E+00, 2.400E+00, 2.900E+00, 3.400E+00, 4.000E+00, &
   &           4.800E+00, 6.000E+00, 7.000E+00, 8.100E+00, 8.900E+00, &
   &           8.700E+00, 7.550E+00, 5.900E+00, 4.500E+00, 3.500E+00, &
   &           2.800E+00, 1.800E+00, 1.300E+00, 8.000E-01, 4.000E-01, &
   &           1.900E-01, 2.000E-01, 5.700E-01, 7.500E-01, 7.000E-01, &
   &           4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04, &
   &           MXZ50*0.0 /
!
   DATA AMOL33 /                                                     &
   &           2.778E-02, 2.800E-02, 2.849E-02, 3.200E-02, 3.567E-02, &
   &           4.720E-02, 5.837E-02, 7.891E-02, 1.039E-01, 1.567E-01, &
   &           2.370E-01, 3.624E-01, 5.232E-01, 7.036E-01, 8.000E-01, &
   &           9.000E-01, 1.100E+00, 1.400E+00, 1.800E+00, 2.300E+00, &
   &           2.900E+00, 3.500E+00, 3.900E+00, 4.300E+00, 4.700E+00, &
   &           5.100E+00, 5.600E+00, 6.100E+00, 6.800E+00, 7.100E+00, &
   &           7.200E+00, 6.900E+00, 5.900E+00, 4.600E+00, 3.700E+00, &
   &           2.750E+00, 1.700E+00, 1.000E-00, 5.500E-01, 3.200E-01, &
   &           2.500E-01, 2.300E-01, 5.500E-01, 8.000E-01, 8.000E-01, &
   &           4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04, &
   &           MXZ50*0.0 /
!
   DATA AMOL43 /                                                     &
   &           2.412E-02, 2.940E-02, 3.379E-02, 3.887E-02, 4.478E-02, &
   &           5.328E-02, 6.564E-02, 7.738E-02, 9.114E-02, 1.420E-01, &
   &           1.890E-01, 3.050E-01, 4.100E-01, 5.000E-01, 6.000E-01, &
   &           7.000E-01, 8.500E-01, 1.000E+00, 1.300E+00, 1.700E+00, &
   &           2.100E+00, 2.700E+00, 3.300E+00, 3.700E+00, 4.200E+00, &
   &           4.500E+00, 5.300E+00, 5.700E+00, 6.900E+00, 7.700E+00, &
   &           7.800E+00, 7.000E+00, 5.400E+00, 4.200E+00, 3.200E+00, &
   &           2.500E+00, 1.700E+00, 1.200E+00, 8.000E-01, 4.000E-01, &
   &           2.000E-01, 1.800E-01, 6.500E-01, 9.000E-01, 8.000E-01, &
   &           4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04, &
   &           MXZ50*0.0 /
!
   DATA AMOL53 /                                                     &
   &           1.802E-02, 2.072E-02, 2.336E-02, 2.767E-02, 3.253E-02, &
   &           3.801E-02, 4.446E-02, 7.252E-02, 1.040E-01, 2.100E-01, &
   &           3.000E-01, 3.500E-01, 4.000E-01, 6.500E-01, 9.000E-01, &
   &           1.200E+00, 1.500E+00, 1.900E+00, 2.450E+00, 3.100E+00, &
   &           3.700E+00, 4.000E+00, 4.200E+00, 4.500E+00, 4.600E+00, &
   &           4.700E+00, 4.900E+00, 5.400E+00, 5.900E+00, 6.200E+00, &
   &           6.250E+00, 5.900E+00, 5.100E+00, 4.100E+00, 3.000E+00, &
   &           2.600E+00, 1.600E+00, 9.500E-01, 6.500E-01, 5.000E-01, &
   &           3.300E-01, 1.300E-01, 7.500E-01, 8.000E-01, 8.000E-01, &
   &           4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04, &
   &           MXZ50*0.0 /
!
   DATA AMOL63 /                                                     &
   &           2.660E-02, 2.931E-02, 3.237E-02, 3.318E-02, 3.387E-02, &
   &           3.768E-02, 4.112E-02, 5.009E-02, 5.966E-02, 9.168E-02, &
   &           1.313E-01, 2.149E-01, 3.095E-01, 3.846E-01, 5.030E-01, &
   &           6.505E-01, 8.701E-01, 1.187E+00, 1.587E+00, 2.030E+00, &
   &           2.579E+00, 3.028E+00, 3.647E+00, 4.168E+00, 4.627E+00, &
   &           5.118E+00, 5.803E+00, 6.553E+00, 7.373E+00, 7.837E+00, &
   &           7.800E+00, 7.300E+00, 6.200E+00, 5.250E+00, 4.100E+00, &
   &           3.100E+00, 1.800E+00, 1.100E+00, 7.000E-01, 3.000E-01, &
   &           2.500E-01, 3.000E-01, 5.000E-01, 7.000E-01, 7.000E-01, &
   &           4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04, &
   &           MXZ50*0.0 /
!
!     DATA  N2O /
!
   DATA AMOL14 /                                                     &
   &           3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, &
   &           3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01, &
   &           3.179E-01, 3.140E-01, 3.095E-01, 3.048E-01, 2.999E-01, &
   &           2.944E-01, 2.877E-01, 2.783E-01, 2.671E-01, 2.527E-01, &
   &           2.365E-01, 2.194E-01, 2.051E-01, 1.967E-01, 1.875E-01, &
   &           1.756E-01, 1.588E-01, 1.416E-01, 1.165E-01, 9.275E-02, &
   &           6.693E-02, 4.513E-02, 2.751E-02, 1.591E-02, 9.378E-03, &
   &           4.752E-03, 3.000E-03, 2.065E-03, 1.507E-03, 1.149E-03, &
   &           8.890E-04, 7.056E-04, 5.716E-04, 4.708E-04, 3.932E-04, &
   &           3.323E-04, 2.837E-04, 2.443E-04, 2.120E-04, 1.851E-04, &
   &           MXZ50*0.0 /
!
   DATA AMOL24 /                                                     &
   &           3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, &
   &           3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01, 3.163E-01, &
   &           3.096E-01, 2.989E-01, 2.936E-01, 2.860E-01, 2.800E-01, &
   &           2.724E-01, 2.611E-01, 2.421E-01, 2.174E-01, 1.843E-01, &
   &           1.607E-01, 1.323E-01, 1.146E-01, 1.035E-01, 9.622E-02, &
   &           8.958E-02, 8.006E-02, 6.698E-02, 4.958E-02, 3.695E-02, &
   &           2.519E-02, 1.736E-02, 1.158E-02, 7.665E-03, 5.321E-03, &
   &           3.215E-03, 2.030E-03, 1.397E-03, 1.020E-03, 7.772E-04, &
   &           6.257E-04, 5.166E-04, 4.352E-04, 3.727E-04, 3.237E-04, &
   &           2.844E-04, 2.524E-04, 2.260E-04, 2.039E-04, 1.851E-04, &
   &           MXZ50*0.0 /
!
   DATA AMOL34 /                                                     &
   &           3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, &
   &           3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01, 3.163E-01, &
   &           3.096E-01, 2.989E-01, 2.936E-01, 2.860E-01, 2.800E-01, &
   &           2.724E-01, 2.611E-01, 2.421E-01, 2.174E-01, 1.843E-01, &
   &           1.621E-01, 1.362E-01, 1.230E-01, 1.124E-01, 1.048E-01, &
   &           9.661E-02, 8.693E-02, 7.524E-02, 6.126E-02, 5.116E-02, &
   &           3.968E-02, 2.995E-02, 2.080E-02, 1.311E-02, 8.071E-03, &
   &           4.164E-03, 2.629E-03, 1.809E-03, 1.321E-03, 1.007E-03, &
   &           7.883E-04, 6.333E-04, 5.194E-04, 4.333E-04, 3.666E-04, &
   &           3.140E-04, 2.717E-04, 2.373E-04, 2.089E-04, 1.851E-04, &
   &           MXZ50*0.0 /
!
   DATA AMOL44 /                                                     &
   &           3.100E-01, 3.100E-01, 3.100E-01, 3.100E-01, 3.079E-01, &
   &           3.024E-01, 2.906E-01, 2.822E-01, 2.759E-01, 2.703E-01, &
   &           2.651E-01, 2.600E-01, 2.549E-01, 2.494E-01, 2.433E-01, &
   &           2.355E-01, 2.282E-01, 2.179E-01, 2.035E-01, 1.817E-01, &
   &           1.567E-01, 1.350E-01, 1.218E-01, 1.102E-01, 9.893E-02, &
   &           8.775E-02, 7.327E-02, 5.941E-02, 4.154E-02, 3.032E-02, &
   &           1.949E-02, 1.274E-02, 9.001E-03, 6.286E-03, 4.558E-03, &
   &           2.795E-03, 1.765E-03, 1.214E-03, 8.866E-04, 6.756E-04, &
   &           5.538E-04, 4.649E-04, 3.979E-04, 3.459E-04, 3.047E-04, &
   &           2.713E-04, 2.439E-04, 2.210E-04, 2.017E-04, 1.851E-04, &
   &           MXZ50*0.0 /
!
   DATA AMOL54 /                                                     &
   &           3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, &
   &           3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01, 3.163E-01, &
   &           3.096E-01, 2.989E-01, 2.936E-01, 2.860E-01, 2.800E-01, &
   &           2.724E-01, 2.611E-01, 2.421E-01, 2.174E-01, 1.843E-01, &
   &           1.621E-01, 1.362E-01, 1.230E-01, 1.122E-01, 1.043E-01, &
   &           9.570E-02, 8.598E-02, 7.314E-02, 5.710E-02, 4.670E-02, &
   &           3.439E-02, 2.471E-02, 1.631E-02, 1.066E-02, 7.064E-03, &
   &           3.972E-03, 2.508E-03, 1.726E-03, 1.260E-03, 9.602E-04, &
   &           7.554E-04, 6.097E-04, 5.024E-04, 4.210E-04, 3.579E-04, &
   &           3.080E-04, 2.678E-04, 2.350E-04, 2.079E-04, 1.851E-04, &
   &           MXZ50*0.0 /
!
   DATA AMOL64 /                                                     &
   &           3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, &
   &           3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01, &
   &           3.179E-01, 3.140E-01, 3.095E-01, 3.048E-01, 2.999E-01, &
   &           2.944E-01, 2.877E-01, 2.783E-01, 2.671E-01, 2.527E-01, &
   &           2.365E-01, 2.194E-01, 2.051E-01, 1.967E-01, 1.875E-01, &
   &           1.756E-01, 1.588E-01, 1.416E-01, 1.165E-01, 9.275E-02, &
   &           6.693E-02, 4.513E-02, 2.751E-02, 1.591E-02, 9.378E-03, &
   &           4.752E-03, 3.000E-03, 2.065E-03, 1.507E-03, 1.149E-03, &
   &           8.890E-04, 7.056E-04, 5.716E-04, 4.708E-04, 3.932E-04, &
   &           3.323E-04, 2.837E-04, 2.443E-04, 2.120E-04, 1.851E-04, &
   &           MXZ50*0.0 /
!
!     DATA CO /
!
   DATA AMOL15 /                                                     &
   &           1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01, &
   &           1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01, &
   &           9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02, &
   &           3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02, &
   &           1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02, &
   &           1.521E-02, 1.722E-02, 1.995E-02, 2.266E-02, 2.487E-02, &
   &           2.738E-02, 3.098E-02, 3.510E-02, 3.987E-02, 4.482E-02, &
   &           5.092E-02, 5.985E-02, 6.960E-02, 9.188E-02, 1.938E-01, &
   &           5.688E-01, 1.549E+00, 3.849E+00, 6.590E+00, 1.044E+01, &
   &           1.705E+01, 2.471E+01, 3.358E+01, 4.148E+01, 5.000E+01, &
   &           MXZ50*0.0 /
!
   DATA AMOL25 /                                                     &
   &           1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01, &
   &           1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01, &
   &           9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02, &
   &           3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02, &
   &           1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02, &
   &           1.521E-02, 1.722E-02, 1.995E-02, 2.266E-02, 2.487E-02, &
   &           2.716E-02, 2.962E-02, 3.138E-02, 3.307E-02, 3.487E-02, &
   &           3.645E-02, 3.923E-02, 4.673E-02, 6.404E-02, 1.177E-01, &
   &           2.935E-01, 6.815E-01, 1.465E+00, 2.849E+00, 5.166E+00, &
   &           1.008E+01, 1.865E+01, 2.863E+01, 3.890E+01, 5.000E+01, &
   &           MXZ50*0.0 /
!
   DATA AMOL35 /                                                     &
   &           1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01, &
   &           1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01, &
   &           9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02, &
   &           3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02, &
   &           1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02, &
   &           1.498E-02, 1.598E-02, 1.710E-02, 1.850E-02, 1.997E-02, &
   &           2.147E-02, 2.331E-02, 2.622E-02, 3.057E-02, 3.803E-02, &
   &           6.245E-02, 1.480E-01, 2.926E-01, 5.586E-01, 1.078E+00, &
   &           1.897E+00, 2.960E+00, 4.526E+00, 6.862E+00, 1.054E+01, &
   &           1.709E+01, 2.473E+01, 3.359E+01, 4.149E+01, 5.000E+01, &
   &           MXZ50*0.0 /
!
   DATA AMOL45 /                                                     &
   &           1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01, &
   &           1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01, &
   &           9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02, &
   &           3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02, &
   &           1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02, &
   &           1.510E-02, 1.649E-02, 1.808E-02, 1.997E-02, 2.183E-02, &
   &           2.343E-02, 2.496E-02, 2.647E-02, 2.809E-02, 2.999E-02, &
   &           3.220E-02, 3.650E-02, 4.589E-02, 6.375E-02, 1.176E-01, &
   &           3.033E-01, 7.894E-01, 1.823E+00, 3.402E+00, 5.916E+00, &
   &           1.043E+01, 1.881E+01, 2.869E+01, 3.892E+01, 5.000E+01, &
   &           MXZ50*0.0 /
!
   DATA AMOL55 /                                                     &
   &           1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01, &
   &           1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01, &
   &           9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02, &
   &           3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02, &
   &           1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02, &
   &           1.521E-02, 1.722E-02, 2.037E-02, 2.486E-02, 3.168E-02, &
   &           4.429E-02, 6.472E-02, 1.041E-01, 1.507E-01, 2.163E-01, &
   &           3.141E-01, 4.842E-01, 7.147E-01, 1.067E+00, 1.516E+00, &
   &           2.166E+00, 3.060E+00, 4.564E+00, 6.877E+00, 1.055E+01, &
   &           1.710E+01, 2.473E+01, 3.359E+01, 4.149E+01, 5.000E+01, &
   &           MXZ50*0.0 /
!
   DATA AMOL65 /                                                     &
   &           1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01, &
   &           1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01, &
   &           9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02, &
   &           3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02, &
   &           1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02, &
   &           1.498E-02, 1.598E-02, 1.710E-02, 1.850E-02, 2.009E-02, &
   &           2.220E-02, 2.497E-02, 2.824E-02, 3.241E-02, 3.717E-02, &
   &           4.597E-02, 6.639E-02, 1.073E-01, 1.862E-01, 3.059E-01, &
   &           6.375E-01, 1.497E+00, 3.239E+00, 5.843E+00, 1.013E+01, &
   &           1.692E+01, 2.467E+01, 3.356E+01, 4.148E+01, 5.000E+01, &
   &           MXZ50*0.0 /
!
!     DATA  CH4 /
!
   DATA AMOL16 /                                                     &
   &           1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, &
   &           1.700E+00, 1.700E+00, 1.699E+00, 1.697E+00, 1.693E+00, &
   &           1.685E+00, 1.675E+00, 1.662E+00, 1.645E+00, 1.626E+00, &
   &           1.605E+00, 1.582E+00, 1.553E+00, 1.521E+00, 1.480E+00, &
   &           1.424E+00, 1.355E+00, 1.272E+00, 1.191E+00, 1.118E+00, &
   &           1.055E+00, 9.870E-01, 9.136E-01, 8.300E-01, 7.460E-01, &
   &           6.618E-01, 5.638E-01, 4.614E-01, 3.631E-01, 2.773E-01, &
   &           2.100E-01, 1.651E-01, 1.500E-01, 1.500E-01, 1.500E-01, &
   &           1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01, &
   &           1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02, &
   &           MXZ50*0.0 /
!
   DATA AMOL26 /                                                     &
   &           1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.697E+00, &
   &           1.687E+00, 1.672E+00, 1.649E+00, 1.629E+00, 1.615E+00, &
   &           1.579E+00, 1.542E+00, 1.508E+00, 1.479E+00, 1.451E+00, &
   &           1.422E+00, 1.390E+00, 1.356E+00, 1.323E+00, 1.281E+00, &
   &           1.224E+00, 1.154E+00, 1.066E+00, 9.730E-01, 8.800E-01, &
   &           7.888E-01, 7.046E-01, 6.315E-01, 5.592E-01, 5.008E-01, &
   &           4.453E-01, 3.916E-01, 3.389E-01, 2.873E-01, 2.384E-01, &
   &           1.944E-01, 1.574E-01, 1.500E-01, 1.500E-01, 1.500E-01, &
   &           1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01, &
   &           1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02, &
   &           MXZ50*0.0 /
!
   DATA AMOL36 /                                                     &
   &           1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.697E+00, &
   &           1.687E+00, 1.672E+00, 1.649E+00, 1.629E+00, 1.615E+00, &
   &           1.579E+00, 1.542E+00, 1.508E+00, 1.479E+00, 1.451E+00, &
   &           1.422E+00, 1.390E+00, 1.356E+00, 1.323E+00, 1.281E+00, &
   &           1.224E+00, 1.154E+00, 1.066E+00, 9.730E-01, 8.800E-01, &
   &           7.931E-01, 7.130E-01, 6.438E-01, 5.746E-01, 5.050E-01, &
   &           4.481E-01, 3.931E-01, 3.395E-01, 2.876E-01, 2.386E-01, &
   &           1.944E-01, 1.574E-01, 1.500E-01, 1.500E-01, 1.500E-01, &
   &           1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01, &
   &           1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02, &
   &           MXZ50*0.0 /
!
   DATA AMOL46 /                                                     &
   &           1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.697E+00, &
   &           1.687E+00, 1.672E+00, 1.649E+00, 1.629E+00, 1.615E+00, &
   &           1.579E+00, 1.542E+00, 1.506E+00, 1.471E+00, 1.434E+00, &
   &           1.389E+00, 1.342E+00, 1.290E+00, 1.230E+00, 1.157E+00, &
   &           1.072E+00, 9.903E-01, 9.170E-01, 8.574E-01, 8.013E-01, &
   &           7.477E-01, 6.956E-01, 6.442E-01, 5.888E-01, 5.240E-01, &
   &           4.506E-01, 3.708E-01, 2.992E-01, 2.445E-01, 2.000E-01, &
   &           1.660E-01, 1.500E-01, 1.500E-01, 1.500E-01, 1.500E-01, &
   &           1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01, &
   &           1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02, &
   &           MXZ50*0.0 /
!
   DATA AMOL56 /                                                     &
   &           1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.697E+00, &
   &           1.687E+00, 1.672E+00, 1.649E+00, 1.629E+00, 1.615E+00, &
   &           1.579E+00, 1.542E+00, 1.506E+00, 1.471E+00, 1.434E+00, &
   &           1.389E+00, 1.342E+00, 1.290E+00, 1.230E+00, 1.161E+00, &
   &           1.084E+00, 1.014E+00, 9.561E-01, 9.009E-01, 8.479E-01, &
   &           7.961E-01, 7.449E-01, 6.941E-01, 6.434E-01, 5.883E-01, &
   &           5.238E-01, 4.505E-01, 3.708E-01, 3.004E-01, 2.453E-01, &
   &           1.980E-01, 1.590E-01, 1.500E-01, 1.500E-01, 1.500E-01, &
   &           1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01, &
   &           1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02, &
   &           MXZ50*0.0 /
!
   DATA AMOL66 /                                                     &
   &           1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, &
   &           1.700E+00, 1.700E+00, 1.699E+00, 1.697E+00, 1.693E+00, &
   &           1.685E+00, 1.675E+00, 1.662E+00, 1.645E+00, 1.626E+00, &
   &           1.605E+00, 1.582E+00, 1.553E+00, 1.521E+00, 1.480E+00, &
   &           1.424E+00, 1.355E+00, 1.272E+00, 1.191E+00, 1.118E+00, &
   &           1.055E+00, 9.870E-01, 9.136E-01, 8.300E-01, 7.460E-01, &
   &           6.618E-01, 5.638E-01, 4.614E-01, 3.631E-01, 2.773E-01, &
   &           2.100E-01, 1.650E-01, 1.500E-01, 1.500E-01, 1.500E-01, &
   &           1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01, &
   &           1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02, &
   &           MXZ50*0.0 /
!
!     DATA O2 /
!
   DATA AMOL17 /                                                     &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05, &
   &           1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04, &
   &           MXZ50*0.0 /
!
   DATA AMOL27 /                                                     &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05, &
   &           1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04, &
   &           MXZ50*0.0 /
!
   DATA AMOL37 /                                                     &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05, &
   &           1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04, &
   &           MXZ50*0.0 /
!
   DATA AMOL47 /                                                     &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05, &
   &           1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04, &
   &           MXZ50*0.0 /
!
   DATA AMOL57 /                                                     &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05, &
   &           1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04, &
   &           MXZ50*0.0 /
!
   DATA AMOL67 /                                                     &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, &
   &           2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05, &
   &           1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04, &
   &           MXZ50*0.0 /
!
!     DATA DENSITY /
!
   DATA AMOL18 /                                                     &
   &           2.450E+19, 2.231E+19, 2.028E+19, 1.827E+19, 1.656E+19, &
   &           1.499E+19, 1.353E+19, 1.218E+19, 1.095E+19, 9.789E+18, &
   &           8.747E+18, 7.780E+18, 6.904E+18, 6.079E+18, 5.377E+18, &
   &           4.697E+18, 4.084E+18, 3.486E+18, 2.877E+18, 2.381E+18, &
   &           1.981E+18, 1.651E+18, 1.381E+18, 1.169E+18, 9.920E+17, &
   &           8.413E+17, 5.629E+17, 3.807E+17, 2.598E+17, 1.789E+17, &
   &           1.243E+17, 8.703E+16, 6.147E+16, 4.352E+16, 3.119E+16, &
   &           2.291E+16, 1.255E+16, 6.844E+15, 3.716E+15, 1.920E+15, &
   &           9.338E+14, 4.314E+14, 1.801E+14, 7.043E+13, 2.706E+13, &
   &           1.098E+13, 4.445E+12, 1.941E+12, 8.706E+11, 4.225E+11, &
   &           MXZ50*0.0 /
!
   DATA AMOL28 /                                                     &
   &           2.496E+19, 2.257E+19, 2.038E+19, 1.843E+19, 1.666E+19, &
   &           1.503E+19, 1.351E+19, 1.212E+19, 1.086E+19, 9.716E+18, &
   &           8.656E+18, 7.698E+18, 6.814E+18, 6.012E+18, 5.141E+18, &
   &           4.368E+18, 3.730E+18, 3.192E+18, 2.715E+18, 2.312E+18, &
   &           1.967E+18, 1.677E+18, 1.429E+18, 1.223E+18, 1.042E+18, &
   &           8.919E+17, 6.050E+17, 4.094E+17, 2.820E+17, 1.927E+17, &
   &           1.338E+17, 9.373E+16, 6.624E+16, 4.726E+16, 3.398E+16, &
   &           2.500E+16, 1.386E+16, 7.668E+15, 4.196E+15, 2.227E+15, &
   &           1.109E+15, 4.996E+14, 1.967E+14, 7.204E+13, 2.541E+13, &
   &           9.816E+12, 3.816E+12, 1.688E+12, 8.145E+11, 4.330E+11, &
   &           MXZ50*0.0 /
!
   DATA AMOL38 /                                                     &
   &           2.711E+19, 2.420E+19, 2.158E+19, 1.922E+19, 1.724E+19, &
   &           1.542E+19, 1.376E+19, 1.225E+19, 1.086E+19, 9.612E+18, &
   &           8.472E+18, 7.271E+18, 6.237E+18, 5.351E+18, 4.588E+18, &
   &           3.931E+18, 3.368E+18, 2.886E+18, 2.473E+18, 2.115E+18, &
   &           1.809E+18, 1.543E+18, 1.317E+18, 1.125E+18, 9.633E+17, &
   &           8.218E+17, 5.536E+17, 3.701E+17, 2.486E+17, 1.647E+17, &
   &           1.108E+17, 7.540E+16, 5.202E+16, 3.617E+16, 2.570E+16, &
   &           1.863E+16, 1.007E+16, 5.433E+15, 2.858E+15, 1.477E+15, &
   &           7.301E+14, 3.553E+14, 1.654E+14, 7.194E+13, 3.052E+13, &
   &           1.351E+13, 6.114E+12, 2.952E+12, 1.479E+12, 7.836E+11, &
   &           MXZ50*0.0 /
!
   DATA AMOL48 /                                                     &
   &           2.549E+19, 2.305E+19, 2.080E+19, 1.873E+19, 1.682E+19, &
   &           1.508E+19, 1.357E+19, 1.216E+19, 1.088E+19, 9.701E+18, &
   &           8.616E+18, 7.402E+18, 6.363E+18, 5.471E+18, 4.699E+18, &
   &           4.055E+18, 3.476E+18, 2.987E+18, 2.568E+18, 2.208E+18, &
   &           1.899E+18, 1.632E+18, 1.403E+18, 1.207E+18, 1.033E+18, &
   &           8.834E+17, 6.034E+17, 4.131E+17, 2.839E+17, 1.938E+17, &
   &           1.344E+17, 9.402E+16, 6.670E+16, 4.821E+16, 3.516E+16, &
   &           2.581E+16, 1.421E+16, 7.946E+15, 4.445E+15, 2.376E+15, &
   &           1.198E+15, 5.311E+14, 2.022E+14, 7.221E+13, 2.484E+13, &
   &           9.441E+12, 3.624E+12, 1.610E+12, 7.951E+11, 4.311E+11, &
   &           MXZ50*0.0 /
!
   DATA AMOL58 /                                                     &
   &           2.855E+19, 2.484E+19, 2.202E+19, 1.950E+19, 1.736E+19, &
   &           1.552E+19, 1.383E+19, 1.229E+19, 1.087E+19, 9.440E+18, &
   &           8.069E+18, 6.898E+18, 5.893E+18, 5.039E+18, 4.308E+18, &
   &           3.681E+18, 3.156E+18, 2.704E+18, 2.316E+18, 1.982E+18, &
   &           1.697E+18, 1.451E+18, 1.241E+18, 1.061E+18, 9.065E+17, &
   &           7.742E+17, 5.134E+17, 3.423E+17, 2.292E+17, 1.533E+17, &
   &           1.025E+17, 6.927E+16, 4.726E+16, 3.266E+16, 2.261E+16, &
   &           1.599E+16, 8.364E+15, 4.478E+15, 2.305E+15, 1.181E+15, &
   &           6.176E+14, 3.127E+14, 1.531E+14, 7.244E+13, 3.116E+13, &
   &           1.403E+13, 6.412E+12, 3.099E+12, 1.507E+12, 7.814E+11, &
   &           MXZ50*0.0 /
!
   DATA AMOL68 /                                                     &
   &           2.548E+19, 2.313E+19, 2.094E+19, 1.891E+19, 1.704E+19, &
   &           1.532E+19, 1.373E+19, 1.228E+19, 1.094E+19, 9.719E+18, &
   &           8.602E+18, 7.589E+18, 6.489E+18, 5.546E+18, 4.739E+18, &
   &           4.050E+18, 3.462E+18, 2.960E+18, 2.530E+18, 2.163E+18, &
   &           1.849E+18, 1.575E+18, 1.342E+18, 1.144E+18, 9.765E+17, &
   &           8.337E+17, 5.640E+17, 3.830E+17, 2.524E+17, 1.761E+17, &
   &           1.238E+17, 8.310E+16, 5.803E+16, 4.090E+16, 2.920E+16, &
   &           2.136E+16, 1.181E+16, 6.426E+15, 3.386E+15, 1.723E+15, &
   &           8.347E+14, 3.832E+14, 1.711E+14, 7.136E+13, 2.924E+13, &
   &           1.189E+13, 5.033E+12, 2.144E+12, 9.688E+11, 5.114E+11, &
   &           MXZ50*0.0 /
!
   DATA ANO /  3.00E-04,  3.00E-04,  3.00E-04,  3.00E-04,  3.00E-04, &
   &            3.00E-04,  3.00E-04,  3.00E-04,  3.00E-04,  3.00E-04, &
   &            3.00E-04,  3.00E-04,  3.00E-04,  2.99E-04,  2.95E-04, &
   &            2.83E-04,  2.68E-04,  2.52E-04,  2.40E-04,  2.44E-04, &
   &            2.55E-04,  2.77E-04,  3.07E-04,  3.60E-04,  4.51E-04, &
   &            6.85E-04,  1.28E-03,  2.45E-03,  4.53E-03,  7.14E-03, &
   &            9.34E-03,  1.12E-02,  1.19E-02,  1.17E-02,  1.10E-02, &
   &            1.03E-02,  1.01E-02,  1.01E-02,  1.03E-02,  1.15E-02, &
   &            1.61E-02,  2.68E-02,  7.01E-02,  2.13E-01,  7.12E-01, &
   &            2.08E+00,  4.50E+00,  7.98E+00,  1.00E+01,  1.00E+01, &
   &            MXZ50*0.0 /
!
   DATA SO2 /  3.00E-04,  2.74E-04,  2.36E-04,  1.90E-04,  1.46E-04, &
   &            1.18E-04,  9.71E-05,  8.30E-05,  7.21E-05,  6.56E-05, &
   &            6.08E-05,  5.79E-05,  5.60E-05,  5.59E-05,  5.64E-05, &
   &            5.75E-05,  5.75E-05,  5.37E-05,  4.78E-05,  3.97E-05, &
   &            3.19E-05,  2.67E-05,  2.28E-05,  2.07E-05,  1.90E-05, &
   &            1.75E-05,  1.54E-05,  1.34E-05,  1.21E-05,  1.16E-05, &
   &            1.21E-05,  1.36E-05,  1.65E-05,  2.10E-05,  2.77E-05, &
   &            3.56E-05,  4.59E-05,  5.15E-05,  5.11E-05,  4.32E-05, &
   &            2.83E-05,  1.33E-05,  5.56E-06,  2.24E-06,  8.96E-07, &
   &            3.58E-07,  1.43E-07,  5.73E-08,  2.29E-08,  9.17E-09, &
   &            MXZ50*0.0 /
!
   DATA ANO2 / 2.30E-05,  2.30E-05,  2.30E-05,  2.30E-05,  2.30E-05, &
   &            2.30E-05,  2.30E-05,  2.30E-05,  2.30E-05,  2.32E-05, &
   &            2.38E-05,  2.62E-05,  3.15E-05,  4.45E-05,  7.48E-05, &
   &            1.71E-04,  3.19E-04,  5.19E-04,  7.71E-04,  1.06E-03, &
   &            1.39E-03,  1.76E-03,  2.16E-03,  2.58E-03,  3.06E-03, &
   &            3.74E-03,  4.81E-03,  6.16E-03,  7.21E-03,  7.28E-03, &
   &            6.26E-03,  4.03E-03,  2.17E-03,  1.15E-03,  6.66E-04, &
   &            4.43E-04,  3.39E-04,  2.85E-04,  2.53E-04,  2.31E-04, &
   &            2.15E-04,  2.02E-04,  1.92E-04,  1.83E-04,  1.76E-04, &
   &            1.70E-04,  1.64E-04,  1.59E-04,  1.55E-04,  1.51E-04, &
   &            MXZ50*0.0 /
!
   DATA ANH3 / 5.00E-04,  5.00E-04,  4.63E-04,  3.80E-04,  2.88E-04, &
   &            2.04E-04,  1.46E-04,  9.88E-05,  6.48E-05,  3.77E-05, &
   &            2.03E-05,  1.09E-05,  6.30E-06,  3.12E-06,  1.11E-06, &
   &            4.47E-07,  2.11E-07,  1.10E-07,  6.70E-08,  3.97E-08, &
   &            2.41E-08,  1.92E-08,  1.72E-08,  1.59E-08,  1.44E-08, &
   &            1.23E-08,  9.37E-09,  6.35E-09,  3.68E-09,  1.82E-09, &
   &            9.26E-10,  2.94E-10,  8.72E-11,  2.98E-11,  1.30E-11, &
   &            7.13E-12,  4.80E-12,  3.66E-12,  3.00E-12,  2.57E-12, &
   &            2.27E-12,  2.04E-12,  1.85E-12,  1.71E-12,  1.59E-12, &
   &            1.48E-12,  1.40E-12,  1.32E-12,  1.25E-12,  1.19E-12, &
   &            MXZ50*0.0 /
!
   DATA HNO3 / 5.00E-05,  5.96E-05,  6.93E-05,  7.91E-05,  8.87E-05, &
   &            9.75E-05,  1.11E-04,  1.26E-04,  1.39E-04,  1.53E-04, &
   &            1.74E-04,  2.02E-04,  2.41E-04,  2.76E-04,  3.33E-04, &
   &            4.52E-04,  7.37E-04,  1.31E-03,  2.11E-03,  3.17E-03, &
   &            4.20E-03,  4.94E-03,  5.46E-03,  5.74E-03,  5.84E-03, &
   &            5.61E-03,  4.82E-03,  3.74E-03,  2.59E-03,  1.64E-03, &
   &            9.68E-04,  5.33E-04,  2.52E-04,  1.21E-04,  7.70E-05, &
   &            5.55E-05,  4.45E-05,  3.84E-05,  3.49E-05,  3.27E-05, &
   &            3.12E-05,  3.01E-05,  2.92E-05,  2.84E-05,  2.78E-05, &
   &            2.73E-05,  2.68E-05,  2.64E-05,  2.60E-05,  2.57E-05, &
   &            MXZ50*0.0 /
!
   DATA OH /   4.40E-08,  4.40E-08,  4.40E-08,  4.40E-08,  4.40E-08, &
   &            4.40E-08,  4.40E-08,  4.41E-08,  4.45E-08,  4.56E-08, &
   &            4.68E-08,  4.80E-08,  4.94E-08,  5.19E-08,  5.65E-08, &
   &            6.75E-08,  8.25E-08,  1.04E-07,  1.30E-07,  1.64E-07, &
   &            2.16E-07,  3.40E-07,  5.09E-07,  7.59E-07,  1.16E-06, &
   &            2.18E-06,  5.00E-06,  1.17E-05,  3.40E-05,  8.35E-05, &
   &            1.70E-04,  2.85E-04,  4.06E-04,  5.11E-04,  5.79E-04, &
   &            6.75E-04,  9.53E-04,  1.76E-03,  3.74E-03,  7.19E-03, &
   &            1.12E-02,  1.13E-02,  6.10E-03,  1.51E-03,  2.42E-04, &
   &            4.47E-05,  1.77E-05,  1.19E-05,  1.35E-05,  2.20E-05, &
   &            MXZ50*0.0 /
!
   DATA HF /   1.00E-08,  1.00E-08,  1.23E-08,  1.97E-08,  3.18E-08, &
   &            5.63E-08,  9.18E-08,  1.53E-07,  2.41E-07,  4.04E-07, &
   &            6.57E-07,  1.20E-06,  1.96E-06,  3.12E-06,  4.62E-06, &
   &            7.09E-06,  1.05E-05,  1.69E-05,  2.57E-05,  4.02E-05, &
   &            5.77E-05,  7.77E-05,  9.90E-05,  1.23E-04,  1.50E-04, &
   &            1.82E-04,  2.30E-04,  2.83E-04,  3.20E-04,  3.48E-04, &
   &            3.72E-04,  3.95E-04,  4.10E-04,  4.21E-04,  4.24E-04, &
   &            4.25E-04,  4.25E-04,  4.25E-04,  4.25E-04,  4.25E-04, &
   &            4.25E-04,  4.25E-04,  4.25E-04,  4.25E-04,  4.25E-04, &
   &            4.25E-04,  4.25E-04,  4.25E-04,  4.25E-04,  4.25E-04, &
   &            MXZ50*0.0 /
!
   DATA HCL /  1.00E-03,  7.49E-04,  5.61E-04,  4.22E-04,  3.19E-04, &
   &            2.39E-04,  1.79E-04,  1.32E-04,  9.96E-05,  7.48E-05, &
   &            5.68E-05,  4.59E-05,  4.36E-05,  6.51E-05,  1.01E-04, &
   &            1.63E-04,  2.37E-04,  3.13E-04,  3.85E-04,  4.42E-04, &
   &            4.89E-04,  5.22E-04,  5.49E-04,  5.75E-04,  6.04E-04, &
   &            6.51E-04,  7.51E-04,  9.88E-04,  1.28E-03,  1.57E-03, &
   &            1.69E-03,  1.74E-03,  1.76E-03,  1.79E-03,  1.80E-03, &
   &            1.80E-03,  1.80E-03,  1.80E-03,  1.80E-03,  1.80E-03, &
   &            1.80E-03,  1.80E-03,  1.80E-03,  1.80E-03,  1.80E-03, &
   &            1.80E-03,  1.80E-03,  1.80E-03,  1.80E-03,  1.80E-03, &
   &            MXZ50*0.0 /
!
   DATA HBR /  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06, &
   &            1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06, &
   &            1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06, &
   &            1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06, &
   &            1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06, &
   &            1.71E-06,  1.76E-06,  1.90E-06,  2.26E-06,  2.82E-06, &
   &            3.69E-06,  4.91E-06,  6.13E-06,  6.85E-06,  7.08E-06, &
   &            7.14E-06,  7.15E-06,  7.15E-06,  7.15E-06,  7.15E-06, &
   &            7.15E-06,  7.15E-06,  7.15E-06,  7.15E-06,  7.15E-06, &
   &            7.15E-06,  7.15E-06,  7.15E-06,  7.15E-06,  7.15E-06, &
   &            MXZ50*0.0 /
!
   DATA HI /   3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06, &
   &            3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06, &
   &            3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06, &
   &            3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06, &
   &            3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06, &
   &            3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06, &
   &            3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06, &
   &            3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06, &
   &            3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06, &
   &            3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06, &
   &            MXZ50*0.0 /
!
   DATA CLO /  1.00E-08,  1.00E-08,  1.00E-08,  1.00E-08,  1.00E-08, &
   &            1.00E-08,  1.00E-08,  1.00E-08,  1.01E-08,  1.05E-08, &
   &            1.21E-08,  1.87E-08,  3.18E-08,  5.61E-08,  9.99E-08, &
   &            1.78E-07,  3.16E-07,  5.65E-07,  1.04E-06,  2.04E-06, &
   &            4.64E-06,  8.15E-06,  1.07E-05,  1.52E-05,  2.24E-05, &
   &            3.97E-05,  8.48E-05,  1.85E-04,  3.57E-04,  5.08E-04, &
   &            6.07E-04,  5.95E-04,  4.33E-04,  2.51E-04,  1.56E-04, &
   &            1.04E-04,  7.69E-05,  6.30E-05,  5.52E-05,  5.04E-05, &
   &            4.72E-05,  4.49E-05,  4.30E-05,  4.16E-05,  4.03E-05, &
   &            3.93E-05,  3.83E-05,  3.75E-05,  3.68E-05,  3.61E-05, &
   &            MXZ50*0.0 /
!
   DATA OCS /  6.00E-04,  5.90E-04,  5.80E-04,  5.70E-04,  5.62E-04, &
   &            5.55E-04,  5.48E-04,  5.40E-04,  5.32E-04,  5.25E-04, &
   &            5.18E-04,  5.09E-04,  4.98E-04,  4.82E-04,  4.60E-04, &
   &            4.26E-04,  3.88E-04,  3.48E-04,  3.09E-04,  2.74E-04, &
   &            2.41E-04,  2.14E-04,  1.88E-04,  1.64E-04,  1.37E-04, &
   &            1.08E-04,  6.70E-05,  2.96E-05,  1.21E-05,  4.31E-06, &
   &            1.60E-06,  6.71E-07,  4.35E-07,  3.34E-07,  2.80E-07, &
   &            2.47E-07,  2.28E-07,  2.16E-07,  2.08E-07,  2.03E-07, &
   &            1.98E-07,  1.95E-07,  1.92E-07,  1.89E-07,  1.87E-07, &
   &            1.85E-07,  1.83E-07,  1.81E-07,  1.80E-07,  1.78E-07, &
   &            MXZ50*0.0 /
!
   DATA H2CO / 2.40E-03,  1.07E-03,  4.04E-04,  2.27E-04,  1.40E-04, &
   &            1.00E-04,  7.44E-05,  6.04E-05,  5.01E-05,  4.22E-05, &
   &            3.63E-05,  3.43E-05,  3.39E-05,  3.50E-05,  3.62E-05, &
   &            3.62E-05,  3.58E-05,  3.50E-05,  3.42E-05,  3.39E-05, &
   &            3.43E-05,  3.68E-05,  4.03E-05,  4.50E-05,  5.06E-05, &
   &            5.82E-05,  7.21E-05,  8.73E-05,  1.01E-04,  1.11E-04, &
   &            1.13E-04,  1.03E-04,  7.95E-05,  4.82E-05,  1.63E-05, &
   &            5.10E-06,  2.00E-06,  1.05E-06,  6.86E-07,  5.14E-07, &
   &            4.16E-07,  3.53E-07,  3.09E-07,  2.76E-07,  2.50E-07, &
   &            2.30E-07,  2.13E-07,  1.98E-07,  1.86E-07,  1.75E-07, &
   &            MXZ50*0.0 /
!
   DATA HOCL / 7.70E-06,  1.06E-05,  1.22E-05,  1.14E-05,  9.80E-06, &
   &            8.01E-06,  6.42E-06,  5.42E-06,  4.70E-06,  4.41E-06, &
   &            4.34E-06,  4.65E-06,  5.01E-06,  5.22E-06,  5.60E-06, &
   &            6.86E-06,  8.77E-06,  1.20E-05,  1.63E-05,  2.26E-05, &
   &            3.07E-05,  4.29E-05,  5.76E-05,  7.65E-05,  9.92E-05, &
   &            1.31E-04,  1.84E-04,  2.45E-04,  2.96E-04,  3.21E-04, &
   &            3.04E-04,  2.48E-04,  1.64E-04,  9.74E-05,  4.92E-05, &
   &            2.53E-05,  1.50E-05,  1.05E-05,  8.34E-06,  7.11E-06, &
   &            6.33E-06,  5.78E-06,  5.37E-06,  5.05E-06,  4.78E-06, &
   &            4.56E-06,  4.37E-06,  4.21E-06,  4.06E-06,  3.93E-06, &
   &            MXZ50*0.0 /
!
   DATA AN2 /  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05, &
   &            7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05, &
   &            7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05, &
   &            7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05, &
   &            7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05, &
   &            7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05, &
   &            7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05, &
   &            7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05, &
   &            7.81E+05,  7.81E+05,  7.81E+05,  7.80E+05,  7.79E+05, &
   &            7.77E+05,  7.74E+05,  7.70E+05,  7.65E+05,  7.60E+05, &
   &            MXZ50*0.0 /
!
   DATA HCN /  1.70E-04,  1.65E-04,  1.63E-04,  1.61E-04,  1.60E-04, &
   &            1.60E-04,  1.60E-04,  1.60E-04,  1.60E-04,  1.60E-04, &
   &            1.60E-04,  1.60E-04,  1.60E-04,  1.59E-04,  1.57E-04, &
   &            1.55E-04,  1.52E-04,  1.49E-04,  1.45E-04,  1.41E-04, &
   &            1.37E-04,  1.34E-04,  1.30E-04,  1.25E-04,  1.19E-04, &
   &            1.13E-04,  1.05E-04,  9.73E-05,  9.04E-05,  8.46E-05, &
   &            8.02E-05,  7.63E-05,  7.30E-05,  7.00E-05,  6.70E-05, &
   &            6.43E-05,  6.21E-05,  6.02E-05,  5.88E-05,  5.75E-05, &
   &            5.62E-05,  5.50E-05,  5.37E-05,  5.25E-05,  5.12E-05, &
   &            5.00E-05,  4.87E-05,  4.75E-05,  4.62E-05,  4.50E-05, &
   &            MXZ50*0.0 /
!
   DATA CH3CL/ 7.00E-04,  6.70E-04,  6.43E-04,  6.22E-04,  6.07E-04, &
   &            6.02E-04,  6.00E-04,  6.00E-04,  5.98E-04,  5.94E-04, &
   &            5.88E-04,  5.79E-04,  5.66E-04,  5.48E-04,  5.28E-04, &
   &            5.03E-04,  4.77E-04,  4.49E-04,  4.21E-04,  3.95E-04, &
   &            3.69E-04,  3.43E-04,  3.17E-04,  2.86E-04,  2.48E-04, &
   &            1.91E-04,  1.10E-04,  4.72E-05,  1.79E-05,  7.35E-06, &
   &            3.03E-06,  1.32E-06,  8.69E-07,  6.68E-07,  5.60E-07, &
   &            4.94E-07,  4.56E-07,  4.32E-07,  4.17E-07,  4.05E-07, &
   &            3.96E-07,  3.89E-07,  3.83E-07,  3.78E-07,  3.73E-07, &
   &            3.69E-07,  3.66E-07,  3.62E-07,  3.59E-07,  3.56E-07, &
   &            MXZ50*0.0 /
!
   DATA H2O2 / 2.00E-04,  1.95E-04,  1.92E-04,  1.89E-04,  1.84E-04, &
   &            1.77E-04,  1.66E-04,  1.49E-04,  1.23E-04,  9.09E-05, &
   &            5.79E-05,  3.43E-05,  1.95E-05,  1.08E-05,  6.59E-06, &
   &            4.20E-06,  2.94E-06,  2.30E-06,  2.24E-06,  2.68E-06, &
   &            3.68E-06,  5.62E-06,  1.03E-05,  1.97E-05,  3.70E-05, &
   &            6.20E-05,  1.03E-04,  1.36E-04,  1.36E-04,  1.13E-04, &
   &            8.51E-05,  6.37E-05,  5.17E-05,  4.44E-05,  3.80E-05, &
   &            3.48E-05,  3.62E-05,  5.25E-05,  1.26E-04,  3.77E-04, &
   &            1.12E-03,  2.00E-03,  1.68E-03,  4.31E-04,  4.98E-05, &
   &            6.76E-06,  8.38E-07,  9.56E-08,  1.00E-08,  1.00E-09, &
   &            MXZ50*0.0 /
!
   DATA C2H2 / 3.00E-04,  1.72E-04,  9.57E-05,  6.74E-05,  5.07E-05, &
   &            3.99E-05,  3.19E-05,  2.80E-05,  2.55E-05,  2.40E-05, &
   &            2.27E-05,  2.08E-05,  1.76E-05,  1.23E-05,  7.32E-06, &
   &            4.52E-06,  2.59E-06,  1.55E-06,  8.63E-07,  5.30E-07, &
   &            3.10E-07,  1.89E-07,  1.04E-07,  5.75E-08,  2.23E-08, &
   &            8.51E-09,  4.09E-09,  2.52E-09,  1.86E-09,  1.52E-09, &
   &            1.32E-09,  1.18E-09,  1.08E-09,  9.97E-10,  9.34E-10, &
   &            8.83E-10,  8.43E-10,  8.10E-10,  7.83E-10,  7.60E-10, &
   &            7.40E-10,  7.23E-10,  7.07E-10,  6.94E-10,  6.81E-10, &
   &            6.70E-10,  6.59E-10,  6.49E-10,  6.40E-10,  6.32E-10, &
   &            MXZ50*0.0 /
!
   DATA C2H6 / 2.00E-03,  2.00E-03,  2.00E-03,  2.00E-03,  1.98E-03, &
   &            1.95E-03,  1.90E-03,  1.85E-03,  1.79E-03,  1.72E-03, &
   &            1.58E-03,  1.30E-03,  9.86E-04,  7.22E-04,  4.96E-04, &
   &            3.35E-04,  2.14E-04,  1.49E-04,  1.05E-04,  7.96E-05, &
   &            6.01E-05,  4.57E-05,  3.40E-05,  2.60E-05,  1.89E-05, &
   &            1.22E-05,  5.74E-06,  2.14E-06,  8.49E-07,  3.42E-07, &
   &            1.34E-07,  5.39E-08,  2.25E-08,  1.04E-08,  6.57E-09, &
   &            4.74E-09,  3.79E-09,  3.28E-09,  2.98E-09,  2.79E-09, &
   &            2.66E-09,  2.56E-09,  2.49E-09,  2.43E-09,  2.37E-09, &
   &            2.33E-09,  2.29E-09,  2.25E-09,  2.22E-09,  2.19E-09, &
   &            MXZ50*0.0 /
!
!    Every molecule after here as a small but non-zero default profile
   DATA PH3 /  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA COF2 / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA SF6  / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA H2S  / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA HCOOH /1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA HO2  / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA O    / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA CLONO2  / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA NOPLUS  / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA HOBR  / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA C2H4 / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA CH3OH / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA CH3BR / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA CH3CN / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA CF4 / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA C4H2 / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA HC3N / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA H2   / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA CS   / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /

   DATA SO3  / 1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14, &
   &            MXZ50*0.0 /
!
end block data MLATMB
!
!     ----------------------------------------------------------------
!
SUBROUTINE MDLATM (ITYPE,MDL,IREAD,HSPACE)
!
!     *****************************************************************
!     THIS SUBROUTINE LOADS ONE OF THE 6 BUILT IN ATMOSPHERIC PROFILES
!     OR CALLS NSMDL TO READ IN A USER SUPPLIED PROFILE.
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),          &
   &                MXPNL1(MXLAY),MSLAY1,ISFILE,JSFILE,KSFILE,        &
   &                LSFILE,MSFILE,IEFILE,JEFILE,KEFILE
   COMMON /MSCONS/ AIRMSS(MXLAY),TGRND,SEMIS(3),HMINMS,HMAXMS,       &
   &                MSFLAG,MSWIT,IODFIL,MSTGLE

   COMMON /c_drive/ ref_lat,hobs,ibmax_b,immax_b,                    &
   &                 lvl_1_2,jchar_st(10,2),wm(mxzmd)
! common block for layer-to-level analytical jacobians
   common /dlaydlev/ilevdx,imoldx,iupdwn,                            &
   &    dxdL(mxlay,0:mxmol),dxdU(mxlay,0:mxmol)
!
   character*1 jchar_st
!
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
   COMMON /MLATM/ ALT(MXZMD),PMDL(MXZMD,6),TMDL(MXZMD,6),            &
   &               AMOL(MXZMD,8,6),ZST(MXZMD),PST(MXZMD),TST(MXZMD),  &
   &               AMOLS(MXZMD,MXMOL)
   COMMON /MLATMC/ ATMNAM(6)
   CHARACTER*24 ATMNAM
   COMMON /TRAC/ TRAC(MXZMD,MXTRAC)
   COMMON /DEAMT/ DENM(MXMOL,MXZMD),DENP(MXMOL,MXPDIM),DRYAIR(MXZMD)
!
!     ZMDL BLANK COMMON ALTITUDES FOR LBLRTM BOUNDRIES
!     ZMAX /PARMTR/ HIGHEST LBLRTM ALT
!     ZMIN /PARMTR/ LOWEST LBLRTM ALT
!     ZPTH BLANK COMMON
!     ZST /MLATM/ ORIGINAL LBLRTM ALTITUDES
!
   IF (MDL.EQ.0) GO TO 40
   IF (MDL.GE.1) IMMAX = 50
   DO 30 I = 1, IMMAX
      ZMDL(I) = ALT(I)
      PM(I) = PMDL(I,MDL)
      TM(I) = TMDL(I,MDL)
!
!        > Calculate water density and subtract from <
!        > total density to obtain dry air density  <
!
      DENM(1,I) = AMOL(I,1,MDL)*AMOL(I,8,MDL)*1.0E-6
      DRYAIR(I) = AMOL(I,8,MDL) - DENM(1,I)
      DO 10 K = 1, 7
         IF (K.GT.NMOL) GO TO 10
         DENM(K,I) = AMOL(I,K,MDL)*1.0E-6*DRYAIR(I)
10    CONTINUE

      DENW(I) = DENM(1,I)
      DO 20 K = 8, 28
         IF (K.GT.NMOL) GO TO 30
         ITR = K-7
!
!           < TRAC is the trace constituent information, >
!           < obtained from LBLLOW                       >
!
         DENM(K,I) = TRAC(I,ITR)*1.0E-6*DRYAIR(I)
20    CONTINUE
!
30 END DO
!
   READ (ATMNAM(MDL),900) (HMOD(L),L=1,3)
   GO TO 50
!
40 CALL NSMDL (ITYPE,MDL)
!
   if (imoldx.eq.-99) then
      if (immax.ne.ibmax) then
         write(ipr,*) 'Error in Atmosphere Specification:'
         write(ipr,*) '   Desired levels must match input grid'
         write(ipr,*) '   for analytic jacobian calculation'
         stop 'error in level grid:  see TAPE6'
      endif
   endif

50 ZMIN = ZMDL(1)
!
   DO 70 I = 1, IMMAX
      ZST(I) = ZMDL(I)
      PST(I) = PM(I)
      TST(I) = TM(I)
      IF (HSPACE+0.001.GT.ZMDL(I)) ISPACE = I
      DO 60 M = 1, NMOL
         AMOLS(I,M) = DENM(M,I)
60    CONTINUE
70 END DO
!
   IMMAX = ISPACE
   ZMAX = ZMDL(IMMAX)
!
   IMLOW = IMMAX
!
   RETURN
!
900 FORMAT (3A8)
!
end subroutine MDLATM
!
!     ----------------------------------------------------------------
!
SUBROUTINE NSMDL (ITYPE,MDL)
!
!     *****************************************************************
!
!
!     NOTES TO USER:
!
!     THIS SUBROUTINE IS FOR READING IN AN ATMOSPHERIC PROFILE
!     CORRESPONDING TO MODEL = 0.  THE PROFILE IS READ IN AFTER
!     CONTROL CARD 3.4
!
!     CARD 3.4    IMMAX,(HMOD(I),I=1,3)
!                   (I5,3A8)
!
!             IMMAX  NUMBER OF BOUNDARIES FOR THE PROFILE
!
!             HMOD   A 24 CHARACTER HEADER DESCRIBING THE PROFILE
!
!     SEE DETAILS IN RDUNIT ON CARDS 3.5 AND 3.6.1 ... 3.6.N
!
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4

   COMMON /c_drive/ ref_lat,hobs,ibmax_b,immax_b,                    &
   &                 lvl_1_2,jchar_st(10,2),wm(mxzmd)
!
   character*1 jchar_st
!
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
   COMMON /DEAMT/ DENM(MXMOL,MXZMD),DENP(MXMOL,MXPDIM),DRYAIR(MXZMD)
!
   CHARACTER*8      HMOLS
!
   COMMON /HMOLS/ HMOLS(MXMOL),JUNIT(MXMOL),WMOL(MXMOL),JUNITP,      &
   &              JUNITT
!
!     ***********************************************************
!
!     COMMON BLOCK FOR GENERIC MOLECULAR DATA INPUT
!
!     ************************************************************
!
   if (noprnt .ge.0)  WRITE (IPR,900)
   READ (IRD,905) IMMAX_B,HMOD
!
   IMMAX = ABS(IMMAX_B)
   IMLOW = IMMAX
   if (noprnt .ge.0) WRITE (IPR,910) IMMAX,HMOD
   IF (IMMAX.GT.IMDIM) GO TO 30
!
   DO 20 IM = 1, IMMAX
!
!     READ IN GENERIC UNITS FOR USER MODEL
!
      CALL RDUNIT (IMMAX_B,IM,ZMDL(IM),PM(IM),TM(IM),NMOL,NOPRNT)
!
!     CONVERSION OF GENERIC UNITS TO DENSITIES FOR LBLRTM RUNS
!
      CALL CONVRT (PM(IM),TM(IM),JUNIT,WMOL,IM,NMOL,NOPRNT)
!
      DENW(IM) = DENM(1,IM)
20 END DO
!

   IF (IMMAX_B .LT. 0) THEN
      CALL CMPALT (IMMAX,PM,TM,DENW, ZMDL(1),REF_LAT,ZMDL)
   ENDIF

   DO 25 IM = 2,IMMAX
      IF (ZMDL(IM) .LE. ZMDL(IM-1)) GO TO 35
25 CONTINUE

   RETURN
!
30 CONTINUE
   if (noprnt .ge.0) WRITE (IPR,915) IMMAX,IMDIM
!
   STOP ' LEVEL ERROR IN NSMDL '
!
35 CONTINUE

   if (noprnt .ge.0) WRITE (IPR,920) im,im+1,ZMDL(IM),ZMDL(IM+1)

   STOP 'INPUT ALTITUDES NOT IN ASCENDING ORDER'

900 FORMAT (///,' READING IN USER SUPPLIED MODEL ATMOSPHERE')
905 FORMAT (I5,3A8)
910 FORMAT (//,10X,'IMMAX    = ',I5,/,10X,'PROFILE = ',3A8)
915 FORMAT (/,' NUMBER OF PROFILE LEVELS IMMAX = ',I5,                &
   &        ' EXCEEDS THE MAXIMUM ALLOWED = ',I5)
920 FORMAT (///,' ERROR: INPUT ALTITUDES FOR LBLRTM LAYERS ',         &
   &        'ARE NOT IN ASCENDING ORDER',//,5X,                       &
   &        ' ZMDL AT GRID PT I,I+1 = ',i5,i5,/,(2F10.4))
!
end subroutine NSMDL
!
!     ----------------------------------------------------------------
!
SUBROUTINE HEADPR (IPR,NOPRNT)
!
!     SUBROUTINE TO WRITE HEADER INFORMATION FOR MODEL  0
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!
   CHARACTER*8      HMOLS
!
!
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /HMOLS/ HMOLS(MXMOL),JUNIT(MXMOL),WMOL(MXMOL),JUNITP,      &
   &               JUNITT
!
   WRITE (IPR,900)
   WRITE (IPR,905)
   WRITE (IPR,910) (I,HMOLS(I),I=1,MXMOL)
   WRITE (IPR,915)
!
   RETURN
!
900 FORMAT (/,'  THE USER HAS ELECTED TO PROVIDE THE REQUIRED',/,     &
   &        '  MODEL ATMOSPHERE SPECIFICATIONS.',/,/,                 &
   &        '  SEE DOCUMENTATION OR "SUBROUTINE RDUNIT" FOR ',/,      &
   &        '  ADDITIONAL INFORMATION.',//)
905 FORMAT ('  USER OPTIONS FOR PRESSURE AND TEMPERATURE ',//,        &
   &        '               JCHAR   JUNIT ',//,                       &
   &        '    PRESSURE " ",A      10    PRESSURE IN (MB)',/,       &
   &        '                 B      11       "     "  (ATM)',/,      &
   &        '                 C      12       "     "  (TORR)',/,     &
   &        '                1-6    1-6    DEFAULT TO SPECIFIED',     &
   &        ' MODEL ATMOSPHERE',//,                                   &
   &        '    TEMP     " ",A      10    AMBIENT TEMP IN DEG(K)',/, &
   &        '                 B      11       "     "   "   " (C)',/, &
   &        '                1-6    1-6    DEFAULT TO SPECIFIED',     &
   &        ' MODEL ATMOSPHERE',//)
910 FORMAT (/,' AVAILABLE     ',7('(',I2,')',A8),/,' MOL. SPECIES',   &
   &        (T16,7('(',I2,')',A8)))
915 FORMAT (/,'  POTENTIAL CHOICE OF UNITS FOR ABOVE SPECIES',/,      &
   &        ' JCHAR = " ",A    - VOLUME MIXING RATIO (PPMV)',/,       &
   &        '       = B        - NUMBER DENSITY (CM-3)',/,            &
   &        '       = C        - MASS MIXING RATIO (GM/KG)',/,        &
   &        '       = D        - MASS DENSITY (GM M-3)',/,            &
   &        '       = E        - PARTIAL PRESSURE (MB)',/,            &
   &        '       = F        - DEW POINT TEMP (K) * H2O ONLY *',/,  &
   &        '       = G        - DEW POINT TEMP (C) * H2O ONLY *',/,  &
   &        '       = H        - RELATIVE HUMIDITY (PERCENT) ',       &
   &                             '*H2O ONLY*',/,                      &
   &        '       = I        - AVAILABLE FOR USER DEFINITION',/,    &
   &        '       = 1-6      - DEFAULT TO SPECIFIED MODEL ',        &
   &        'ATMOSPHERE',/,' JCHAR MUST BE LESS THAN "J"',/)
!
end subroutine HEADPR
!
!     ----------------------------------------------------------------
!
SUBROUTINE RDUNIT (IMMAX_B,IM,ZMDL,PM,TM,NMOL,NOPRNT)
!
!     *******************************************************
!
!       SUBROUTINE DESIGNED TO READ NEW MOLECULAR DATA INPUT
!        PARAMETERS - JCHAR = INPUT KEY (SEE BELOW)
!                     WMOL  = INPUT VALUE FOR LAYER
!
!     ***  ROUTINE ALSO ACCEPTS VARIABLE UNITS ON PRESS AND TEMP
!     ***  THE ASSOCIATED 'JUNIT' DEFINITIONS ARE CONTAINED IN
!               JUNITP, AND JUNITT
!          SEE INPUT KEY BELOW
!
!
!       NMOL = NUMBER OF MOLECULAR SPECIES TO BE CONSIDERED
!               (ORDER IS THAT OF AFGL LINE PARAMETER TAPE)
!
!     FOR MOLECULAR SPECIES ONLY
!
!       JCHAR   JUNIT
!
!     " ",A      10    VOLUME MIXING RATIO (PPMV)
!         B      11    NUMBER DENSITY (CM-3)
!         C      12    MASS MIXING RATIO (GM(K)/KG(AIR))
!         D      13    MASS DENSITY (GM M-3)
!         E      14    PARTIAL PRESSURE (MB)
!         F      15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY
!         G      16     "    "     "  (TD IN T(C)) - H2O ONLY
!         H      17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY
!         I      18    AVAILABLE FOR USER DEFINITION
!        1-6    1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE
!                                                (SEE KEY BELOW)
!
!     ****************************************************************
!     ****************************************************************
!
!     ***** OTHER 'JCHAR' SPECIFICATIONS - JCHARP,JCHART
!
!       JCHAR   JUNIT
!
!      " ",A     10    PRESSURE IN (MB)
!          B     11       "     "  (ATM)
!          C     12       "     "  (TORR)
!         1-6   1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE
!
!      " ",A     10    AMBIENT TEMPERATURE IN DEG(K)
!          B     11       "         "       "  " (C)
!          C     12       "         "       "  " (F)
!         1-6   1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE
!
!     ***** DEFINITION OF "DEFAULT" CHOICES FOR PROFILE SELECTION *****
!
!      FOR THE USER WHO WISHES TO ENTER ONLY SELECTED ORIGINAL
!      VERTICAL PROFILES AND WANTS STANDARD ATMOSPHERE SPECIFICATIONS
!      FOR THE OTHERS, THE FOLLOWING OPTION IS AVAILABLE
!
!     *** JCHAR(P,T OR K) MUST = 1-6 (AS ABOVE)
!
!      FOR MOLECULES 8-35, ONLY US STD PROFILES ARE AVIALABLE
!      THEREFORE, WHEN  'JCHAR(K) = 1-5', JCHAR(K) WILL BE RESET TO 6
!
!     *************************************************************
!     *************************************************************
!
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
   CHARACTER*8      HMOLS
!
   CHARACTER*8      HMOD
   COMMON /CMN/HMOD(3),                                              &
   &       ZMDL_st(MXZMD),PM_st(MXZMD),TM_st(MXZMD),RFNDXM_st(MXZMD), &
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
   COMMON /HMOLS/ HMOLS(MXMOL),JUNIT(MXMOL),WMOL(MXMOL),JUNITP,      &
   &               JUNITT
   CHARACTER*1 JCHAR,JCHARP,JCHART,JLONG
!
   COMMON /MCHAR/ JCHAR(MXMOL),JCHARP,JCHART,JLONG

   COMMON /c_drive/ ref_lat,hobs,ibmax_b,immax_dum,                  &
   &                 lvl_1_2,jchar_st(10,2),wm(mxzmd)
!
   character*1 jchar_st
!
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
!
!     ***********************************************************
!
!       COMMON BLOCK FOR GENERIC MOLECULAR DATA INPUT
!
!
!     ************************************************************
!
   DIMENSION JOLD(MXMOL),KUNIT(MXMOL)
!
   DATA JOLD / MXMOL*99 /
   DATA C1 / 18.9766 /,C2 / -14.9595 /,C3 / -2.4388 /
!
   IF (IM.EQ.0) CALL HEADPR (IPR,NOPRNT)
!
!     *********************************************************
!
!     INPUT READ FOR 'MODEL = 0", I.E. USER-SUPPLIED VERITCAL
!
!     **********************************************************
!
   READ (IRD,900) ZMDL,PM,TM,JCHARP,JCHART,JLONG,                    &
   &               (JCHAR(K),K=1,MXMOL)
   ISAME = 0
   JUNITP = JOU(JCHARP)
   JUNITT = JOU(JCHART)
   DO 10 K = 1, NMOL
      JUNIT(K) = JOU(JCHAR(K))
!
!    TEST TO SEE IF INPUT UNITS HAVE CHANGED FROM PREVIOUS READ
!
      IF (JOLD(K).NE.JUNIT(K)) ISAME = 1
      KUNIT(K) = JUNIT(K)+1
10 END DO
!
!     Read in moleclar information at E15.8 format for flag JLONG='L'
   IF (JLONG.EQ.'L') THEN
      READ (IRD,906) (WMOL(K),K=1,NMOL)
   ELSEIF (JLONG.EQ.' ') THEN
      READ (IRD,905) (WMOL(K),K=1,NMOL)
   ELSE
      WRITE(*,*) 'INVALID VALUE FOR JLONG ON RECORD 3.5: ',JLONG
      STOP 'RDUNIT'
   ENDIF
   IF (IM.EQ.0) WRITE (IPR,910)
!
   if (noprnt .ge. 0) then

      IF (JLONG.EQ.'L') THEN
         WRITE (IPR,916) IM,ZMDL,JCHARP,PM,JCHART,TM, (K,JCHAR(K),WMOL( &
            K),K=1,NMOL)
      ELSE
         WRITE (IPR,915) IM,ZMDL,JCHARP,PM,JCHART,TM, (K,JCHAR(K),WMOL( &
            K),K=1,NMOL)
      ENDIF
   endif

   DO 20 I = 1, NMOL
      JOLD(I) = JUNIT(I)
20 END DO
   CALL CHECK (PM,JUNITP,1)
   CALL CHECK (TM,JUNITT,2)
!
!     SUBROUTINE DEFAULT DEFINES WMOL  FOR JCHAR  1-6
!
   IF (IMMAX_B .LT. 0) THEN
      CALL DEFALT_P (PM,TM)
   ELSE
      CALL DEFALT (ZMDL,PM,TM)
   ENDIF

   RETURN
!
900 FORMAT (3E10.3,5X,2A1,1X,A1,1X,47A1)
905 FORMAT (8E10.3)
906 FORMAT (8E15.8)
910 FORMAT (//,'  ECHO INPUT PARAMETERS FOR USER PROVIDED MODEL',/,   &
   &        '0   (P : UNIT)=   ',5X,'(T : UNIT)=   ',5X,              &
   &        '(MOLECULE NUMBER : UNIT)=   ')
915 FORMAT ('0',I4,1X,'(ALT:KM)=',F7.3,4X,'(P:',A1,')=',G11.5,4X,     &
   &        '(T:',A1,')=',F8.3,/,(5X,7(' (',I2,':',A1,')=',1PE10.3)))
916 FORMAT ('0',I4,1X,'(ALT:KM)=',F7.3,4X,'(P:',A1,')=',G11.5,4X,     &
   &        '(T:',A1,')=',F8.3,/,(5X,7(' (',I2,':',A1,')=',1PE15.8)))
!
end subroutine RDUNIT
FUNCTION JOU (CHAR)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
   CHARACTER*1 CHAR,HOLVEC(22)
   DIMENSION INDX1(22)
!
   DATA (HOLVEC(I),I=1,22) /                                         &
   &                '1','2','3','4','5','6','0','0','0','0',' ','A',  &
   &                'B','C','D','E','F','G','H','I','J','K'/
   DATA (INDX1(I),I=1,22) /                                          &
   &                  1,  2,  3,  4,  5,  6,  0,  0,  0,  0, 10, 10,  &
   &                 11, 12, 13, 14, 15, 16, 17, 18, 19, 20/
!
   INDX = 0
   DO 10 I = 1, 22
      IF (HOLVEC(I).NE.CHAR) GO TO 10
      INDX = INDX1(I)
      GO TO 20
10 END DO
20 IF (INDX.EQ.0) THEN
      WRITE (IPR,900) CHAR
      STOP ' JOU: BAD PARAM '
   ENDIF
   JOU = INDX
!
   RETURN
!
900 FORMAT ('0 INVALID PARAMETER :',2X,A1)
!
end function JOU
!
!     ----------------------------------------------------------------
!
SUBROUTINE CHECK (A,IA,KEY)
!
!      UNITS CONVERSION FOR P AND T
!
!     A = P OR T     AND  IA =JUNITP(I.E. MB,ATM,TORR)
!                            =JUNITT(I.E. DEG K OR C)
!                            =JUNITR(I.E. KM,M,OR CM)
!
   DATA PMB / 1013.25 /,PTORR / 760. /,DEGK / 273.15 /
!
   IF (IA.LE.10) RETURN
!
   GO TO (10,20,30) KEY
!
!     PRESSURE CONVERSIONS
!
10 IF (IA.EQ.11) THEN
      A = A*PMB
      RETURN
   ELSEIF (IA.EQ.12) THEN
      A = A*PMB/PTORR
      RETURN
   ELSE
      STOP ' CHECK(P)'
   ENDIF
!
!     TEMPERATURE COMVERSIONS
!
20 IF (IA.LE.11) THEN
      A = A+DEGK
      RETURN
   ELSE
      STOP ' CHECK(T)'
   ENDIF
!
!      RANGE CONVERSIONS
!
30 IF (IA.EQ.11) THEN
      A = A/1.E3
      RETURN
   ELSEIF (IA.EQ.12) THEN
      A = A/1.E5
      RETURN
   ELSE
      STOP ' CHECK(R)'
   ENDIF
!
end subroutine CHECK
!
!     ----------------------------------------------------------------
!
SUBROUTINE DEFALT (Z,P,T)
!
!     *****************************************************************
!
!     THIS SUBROUTINE LOADS ONE OF THE 6 BUILT IN ATMOSPHERIC PROFILES
!     FROM WHICH IT WILL INTERPOLATE "DEFAULT" VALUES FOR ALTITUDE "Z"
!
!
!      ***  THIS SUBROUTINE IS CALLED BY "RDUNIT" WHICH
!      ***  READS USER SUPPLIED INPUT PROFILES OR SINGLE VALUES
!      ***  UNDER "MODEL = 0     " SPECIFICATIONS
!
!      *** SEE DOCUMENTATION FOR CLARIFICATION ***
!
!     SUBROUTINE "DEFALT"IS TRIGGERRED WHENEVER ANY ONE OF
!     THE INPUT PARAMETERS JCHARP, JCART, (JCHAR(K),K=1,NMOL) IS = 1-6
!
!     FOR SIMPLICITY, ALL INTERPOLATIONS ARE DONE AT ONE TIME BECAUSE
!     THE LAGRANGE WEIGHTS (4PT), BASED ON (ALT-Z), REMAIN UNCHANGED
!
!     JCHARP,JCHART AND JCHAR(K) FOR K<8 ALLOW MODEL-DEPENDENT CHOICES
!
!                   JCHAR=JUNIT
!
!                        1       CHOOSES TROPICAL
!                        2         "     MID-LATITUDE SUMMER
!                        3         "     MID-LATITUDE WINTER
!                        4         "     HIGH-LAT SUMMER
!                        5         "     HIGH-LAT WINTER
!                        6         "     US STANDARD
!
!
!     JUNIT(K) FOR K>7 CHOOSES FROM THE SINGLE TRACE CONSTITUENT
!        PROFILES, ALL APPRORIATE FOR THE US STD ATMOSPHERE
!
!     ***  NOTE ***  T<0 WILL ALSO PRINT OUT A MESSAGE INDICATING
!     ***  A POSSIBLE MISAPPLICATION OF TEMPERATURE UNITS, (K) VS (C)
!
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
!
   CHARACTER*8      HMOLS
!
   COMMON /HMOLS/ HMOLS(MXMOL),JUNIT(MXMOL),WMOL(MXMOL),JUNITP,      &
   &               JUNITT
   COMMON /MLATM/ ALT(MXZMD),PMATM(MXZMD,6),TMATM(MXZMD,6),          &
   &               AMOL(MXZMD,8,6),ZST(MXZMD),PST(MXZMD),TST(MXZMD),  &
   &               AMOLS(MXZMD,MXMOL)
   COMMON /MLATMC/ ATMNAM(6)
   CHARACTER*24 ATMNAM
   COMMON /TRAC/ TRAC(MXZMD,MXTRAC)
!
!     *** 4PT INTERPOLATION FUNCTION
!
   VAL(A1,A2,A3,A4,X1,X2,X3,X4) = A1*X1+A2*X2+A3*X3+A4*X4
!
   ILOWER = 0
   IUPPER = 0
   IM50 = 50
   DO 10 IM = 2, IM50
      I2 = IM
      IF (ALT(IM).GE.Z) GO TO 20
10 END DO
   I2 = IM50
20 I1 = I2-1
   I0 = I2-2
   I3 = I2+1
   IF (I0.LT.1) GO TO 30
   IF (I3.GT.IM50) GO TO 40
!
   GO TO 60
!
!     LOWER ENDPOINT CORRECTION
!
30 CONTINUE
   ILOWER = 1
   I0 = I1
   I1 = I2
   I2 = I3
   I3 = I3+1
   GO TO 60
!
!     UPPER ENDPOINT CORRECTION
!
40 CONTINUE
   IUPPER = 1
   IF (Z.GT.ALT(IM50)) GO TO 50
   I3 = I2
   I2 = I1
   I1 = I0
   I0 = I1-1
   GO TO 60
!
!      UPPER ENDPOINT EXTRAPOLATION
!
50 CONTINUE
   Z0 = ALT(I0)
   Z1 = ALT(I1)
   Z2 = ALT(I2)
   Z3 = Z2+2.*(Z-Z2)
   IUPPER = 2
   WRITE (IPR,900) Z
!
   if (JUNITP.le.6.OR.JUNITT.LE.6) STOP 'DEFAULT Z'
   do k=1,nmol
      IF (JUNIT(K).le.6)  STOP 'DEFAULT Z'
   end do
!
!     LAGRANGE CONTINUATION
!
60 CONTINUE
!
!     LAGRANGE COEF DETERMINATION
!
   Z1 = ALT(I1)
   Z2 = ALT(I2)
   Z0 = ALT(I0)
   Z3 = ALT(I3)
   DEN1 = (Z0-Z1)*(Z0-Z2)*(Z0-Z3)
   DEN2 = (Z1-Z2)*(Z1-Z3)*(Z1-Z0)
   DEN3 = (Z2-Z3)*(Z2-Z0)*(Z2-Z1)
   DEN4 = (Z3-Z0)*(Z3-Z1)*(Z3-Z2)
   A1 = ((Z-Z1)*(Z-Z2)*(Z-Z3))/DEN1
   A2 = ((Z-Z2)*(Z-Z3)*(Z-Z0))/DEN2
   A3 = ((Z-Z3)*(Z-Z0)*(Z-Z1))/DEN3
   A4 = ((Z-Z0)*(Z-Z1)*(Z-Z2))/DEN4
!
!     TEST INPUT PARAMETERS (JUNIT'S) SEQUENTIALLY FOR TRIGGER
!      I.E.  JUNIT(P,T,K) = 1-6
!
   IF (JUNITP.GT.6) GO TO 70
   MATM = JUNITP
!
!     WRITE (IPR,60) Z,MATM
!
   X1 =  LOG(PMATM(I0,MATM))
   X2 =  LOG(PMATM(I1,MATM))
   X3 =  LOG(PMATM(I2,MATM))
   X4 =  LOG(PMATM(I3,MATM))
   IF (IUPPER.EQ.2) X4 = X3+2*(X3-X2)
   P = VAL(A1,A2,A3,A4,X1,X2,X3,X4)
   P = EXP(P)
70 IF (JUNITT.GT.6) GO TO 80
   MATM = JUNITT
!
!     WRITE (IPR,65) Z,MATM
!
   X1 = TMATM(I0,MATM)
   X2 = TMATM(I1,MATM)
   X3 = TMATM(I2,MATM)
   X4 = TMATM(I3,MATM)
   T = VAL(A1,A2,A3,A4,X1,X2,X3,X4)
80 DO 110 K = 1, NMOL
      IF (JUNIT(K).GT.6) GO TO 110
!
      IF (K.GT.7) GO TO 90
      MATM = JUNIT(K)
!
!     WRITE (IPR,70) K,HMOLS(K),Z,MATM
!
      X1 = AMOL(I0,K,MATM)
      X2 = AMOL(I1,K,MATM)
      X3 = AMOL(I2,K,MATM)
      X4 = AMOL(I3,K,MATM)
      GO TO 100
90    ITR = K-7
      MATM = 6
!
!     WRITE (IPR,70) K,HMOLS(K),Z,MATM
!
      X1 = TRAC(I0,ITR)
      X2 = TRAC(I1,ITR)
      X3 = TRAC(I2,ITR)
      X4 = TRAC(I3,ITR)
!
100   WMOL(K) = VAL(A1,A2,A3,A4,X1,X2,X3,X4)
      JUNIT(K) = 10
110 END DO
!
   RETURN
!
900 FORMAT (/,'   *** Z IS GREATER THAN 120 KM ***, Z = ',F10.3)
!
end subroutine DEFALT
!
!     ----------------------------------------------------------------
!

SUBROUTINE DEFALT_P (P,T)
!
!     *****************************************************************
!
!     THIS SUBROUTINE LOADS ONE OF THE 6 BUILT IN ATMOSPHERIC PROFILES
!     FROM WHICH IT WILL INTERPOLATE "DEFAULT" VALUES FOR PRESSURE "P"
!
!
!      ***  THIS SUBROUTINE IS CALLED BY "RDUNIT" WHICH
!      ***  READS USER SUPPLIED INPUT PROFILES OR SINGLE VALUES
!      ***  UNDER "MODEL = 0     " SPECIFICATIONS
!
!      *** SEE DOCUMENTATION FOR CLARIFICATION ***
!
!     SUBROUTINE "DEFALT"IS TRIGGERRED WHENEVER ANY ONE OF
!     THE INPUT PARAMETERS JCHARP, JCART, (JCHAR(K),K=1,NMOL) IS = 1-6
!
!     FOR SIMPLICITY, ALL INTERPOLATIONS ARE DONE AT ONE TIME BECAUSE
!     THE LAGRANGE WEIGHTS (4PT), BASED ON (ALT-Z), REMAIN UNCHANGED
!
!     JCHARP,JCHART AND JCHAR(K) FOR K<8 ALLOW MODEL-DEPENDENT CHOICES
!
!                   JCHAR=JUNIT
!
!                        1       CHOOSES TROPICAL
!                        2         "     MID-LATITUDE SUMMER
!                        3         "     MID-LATITUDE WINTER
!                        4         "     HIGH-LAT SUMMER
!                        5         "     HIGH-LAT WINTER
!                        6         "     US STANDARD
!
!
!     JUNIT(K) FOR K>7 CHOOSES FROM THE SINGLE TRACE CONSTITUENT
!        PROFILES, ALL APPRORIATE FOR THE US STD ATMOSPHERE
!
!     ***  NOTE ***  T<0 WILL ALSO PRINT OUT A MESSAGE INDICATING
!     ***  A POSSIBLE MISAPPLICATION OF TEMPERATURE UNITS, (K) VS (C)
!
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
!
   CHARACTER*8      HMOLS
!
   COMMON /HMOLS/ HMOLS(MXMOL),JUNIT(MXMOL),WMOL(MXMOL),JUNITP,      &
   &               JUNITT
   COMMON /MLATM/ ALT(MXZMD),PMATM(MXZMD,6),TMATM(MXZMD,6),          &
   &               AMOL(MXZMD,8,6),ZST(MXZMD),PST(MXZMD),TST(MXZMD),  &
   &               AMOLS(MXZMD,MXMOL)
   COMMON /MLATMC/ ATMNAM(6)
   CHARACTER*24 ATMNAM
   COMMON /TRAC/ TRAC(MXZMD,MXTRAC)
!
!     *** 4PT INTERPOLATION FUNCTION
!
   VAL(A1,A2,A3,A4,X1,X2,X3,X4) = A1*X1+A2*X2+A3*X3+A4*X4
!
   XLOG_P =  LOG(P)
!
   DO 200 J_MDL=1,6

      ILOWER = 0
      IUPPER = 0
      LVL_50 = 50
      DO 10 LVL = 2, LVL_50
         I2 = LVL
         IF (P .GE. PMATM(LVL,J_MDL)) GO TO 20
10    CONTINUE
      I2 = LVL_50
20    I1 = I2-1
      I0 = I2-2
      I3 = I2+1
      IF (I0.LT.1) GO TO 30
      IF (I3.GT.LVL_50) GO TO 40
!
      GO TO 60
!
!     LOWER ENDPOINT CORRECTION
!
30    CONTINUE
      ILOWER = 1
      I0 = I1
      I1 = I2
      I2 = I3
      I3 = I3+1
      GO TO 60
!
!     UPPER ENDPOINT CORRECTION
!
40    CONTINUE
      IUPPER = 1
      IF (P .LE. PMATM(LVL_50,J_MDL)) GO TO 50
      I3 = I2
      I2 = I1
      I1 = I0
      I0 = I1-1
      GO TO 60
!
!      UPPER ENDPOINT EXTRAPOLATION
!
50    CONTINUE
      P_0 =  LOG(PMATM(I0,J_MDL))
      P_1 =  LOG(PMATM(I1,J_MDL))
      P_2 =  LOG(PMATM(I2,J_MDL))
      P_3 = P_2+2.*(XLOG_P-P_2)
      IUPPER = 2
      WRITE (IPR,900) P
!
      STOP 'DEFAULT P'
!
!     LAGRANGE CONTINUATION
!
60    CONTINUE
!
!     LAGRANGE COEF DETERMINATION
!
      P_0 =  LOG(PMATM(I0,J_MDL))
      P_1 =  LOG(PMATM(I1,J_MDL))
      P_2 =  LOG(PMATM(I2,J_MDL))
      P_3 =  LOG(PMATM(I3,J_MDL))
      DEN1 = (P_0-P_1)*(P_0-P_2)*(P_0-P_3)
      DEN2 = (P_1-P_2)*(P_1-P_3)*(P_1-P_0)
      DEN3 = (P_2-P_3)*(P_2-P_0)*(P_2-P_1)
      DEN4 = (P_3-P_0)*(P_3-P_1)*(P_3-P_2)
      A1 = ((XLOG_P-P_1)*(XLOG_P-P_2)*(XLOG_P-P_3))/DEN1
      A2 = ((XLOG_P-P_2)*(XLOG_P-P_3)*(XLOG_P-P_0))/DEN2
      A3 = ((XLOG_P-P_3)*(XLOG_P-P_0)*(XLOG_P-P_1))/DEN3
      A4 = ((XLOG_P-P_0)*(XLOG_P-P_1)*(XLOG_P-P_2))/DEN4
!
!     TEST INPUT PARAMETERS (JUNIT'S) SEQUENTIALLY FOR TRIGGER
!      I.E.  JUNIT(P,T,K) = 1-6
!
!     FOR THIS VERSION OF THE SUBROUTINE DRIVEN BY PRESSURE P
!     JUNITP IS THE MODEL ATMOSPHERES TO BE USED FOR THE ALTITUDE
!
70    IF (JUNITT.GT.6 .OR. JUNITT.NE.J_MDL) GO TO 80
      MATM = JUNITT
!
!     WRITE (IPR,65) P_,MATM
!
      X1 = TMATM(I0,MATM)
      X2 = TMATM(I1,MATM)
      X3 = TMATM(I2,MATM)
      X4 = TMATM(I3,MATM)
      T = VAL(A1,A2,A3,A4,X1,X2,X3,X4)

80    DO 110 K = 1, NMOL
         IF (JUNIT(K).GT.6 .OR. JUNIT(K).NE.J_MDL) GO TO 110
!
         IF (K.GT.7) GO TO 90
         MATM = JUNIT(K)
!
!     WRITE (IPR,70) K,HMOLS(K),P_,MATM
!
         X1 = AMOL(I0,K,MATM)
         X2 = AMOL(I1,K,MATM)
         X3 = AMOL(I2,K,MATM)
         X4 = AMOL(I3,K,MATM)
         GO TO 100
90       ITR = K-7
         MATM = 6
!
!     WRITE (IPR,70) K,HMOLS(K),P_,MATM
!
         X1 = TRAC(I0,ITR)
         X2 = TRAC(I1,ITR)
         X3 = TRAC(I2,ITR)
         X4 = TRAC(I3,ITR)

100      WMOL(K) = VAL(A1,A2,A3,A4,X1,X2,X3,X4)
         JUNIT(K) = 10
!
110   CONTINUE
!
200 CONTINUE

   RETURN
!
900 FORMAT (/,'   *** P IS GREATER THAN P(120 KM)  ***, P = ',        &
   &     1PE10.4)
!
end subroutine DEFALT_P
!
!     ----------------------------------------------------------------
!

SUBROUTINE CONVRT (P,T,JUNIT,WMOL,IM,NMOL,NOPRNT)
!
!*************************************************************
!
!        WRITTEN APR, 1985 TO ACCOMMODATE 'JCHAR' DEFINITIONS FOR
!        UNIFORM DATA INPUT -
!
!      JCHAR    JUNIT
!
!    " ",A       10    VOLUME MIXING RATIO (PPMV)
!        B       11    NUMBER DENSITY (CM-3)
!        C       12    MASS MIXING RATIO (GM(K)/KG(AIR))
!        D       13    MASS DENSITY (GM M-3)
!        E       14    PARTIAL PRESSURE (MB)
!        F       15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY
!        G       16     "    "     "  (TD IN T(C)) - H2O ONLY
!        H       17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY
!        I       18    AVAILABLE FOR USER DEFINITION
!        J       19    REQUEST DEFAULT TO SPECIFIED MODEL ATMOSPHERE
!
!***************************************************************
!
   USE phys_consts, ONLY: avogad, alosmt
   USE planet_consts, ONLY: airmwt
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
   COMMON /DEAMT/ DENM(MXMOL,MXZMD),DENP(MXMOL,MXPDIM),DRYAIR(MXZMD)
!
!
   INTEGER JUNIT(MXMOL)
   DIMENSION WMOL(MXMOL)
   DATA C1 / 18.9766 /,C2 / -14.9595 /,C3 / -2.4388 /
!
   RHOAIR = ALOSMT*(P/PZERO)*(TZERO/T)
   A = TZERO/T
!
!     Get water vapor density
!
   CALL WATVAP (P,T,JUNIT(1),WMOL(1),DENM(1,IM),NOPRNT)
!
!     Determine density of dry air
!
   DRYAIR(IM) = RHOAIR - DENM(1,IM)
!
!     Loop through other molecules
!
   DO 70 K=2,NMOL
      B = AVOGAD/AMWT(K)
      R = AIRMWT/AMWT(K)
      DENM(K,IM) = 0.0

      IF (JUNIT(K).GT.10) GO TO 20
!
!     GIVEN VOL. MIXING RATIO
!
      DENM(K,IM) = WMOL(K)*DRYAIR(IM)*1.E-6
      GO TO 70
20    IF (JUNIT(K).NE.11) GO TO 30
!
!     GIVEN NUMBER DENSITY (CM-3)
!
      DENM(K,IM) = WMOL(K)
      GO TO 70
30    CONTINUE
      IF (JUNIT(K).NE.12) GO TO 40
!
!     GIVEN MASS MIXING RATIO (GM KG-1)
!
      DENM(K,IM) = R*WMOL(K)*1.0E-3*DRYAIR(IM)
      GO TO 70
40    CONTINUE
      IF (JUNIT(K).NE.13) GO TO 50
!
!     GIVEN MASS DENSITY (GM M-3)
!
      DENM(K,IM) = B*WMOL(K)*1.0E-6
      GO TO 70
50    CONTINUE
      IF (JUNIT(K).NE.14) GO TO 60
!
!     GIVEN PARTIAL PRESSURE (MB)
!
      DENM(K,IM) = ALOSMT*(WMOL(K)/PZERO)*(TZERO/T)
      GO TO 70
60    CONTINUE
!
!     JUNIT(18) available for user definition here
!
!
      IF (JUNIT(K).GT.14) THEN
         WRITE (IPR,900) K,JUNIT(K)
         STOP ' CONVRT '
      ENDIF
!
70 END DO
!
900 FORMAT (/,'   **** ERROR IN CONVRT ****, JUNIT(',I5,') = ',I5)
!
   RETURN
!
end subroutine CONVRT
!
!     ----------------------------------------------------------------
!
SUBROUTINE WATVAP (P,T,JUNIT,WMOL,DENNUM,NOPRNT)
!
!**********************************************************************
!
!        WRITTEN APR, 1985 TO ACCOMMODATE 'JCHAR' DEFINITIONS FOR
!        UNIFORM DATA INPUT -
!
!     JCHAR    JUNIT
!
!    " ",A       10    VOLUME MIXING RATIO (PPMV)
!        B       11    NUMBER DENSITY (CM-3)
!        C       12    MASS MIXING RATIO (GM(K)/KG(AIR))
!        D       13    MASS DENSITY (GM M-3)
!        E       14    PARTIAL PRESSURE (MB)
!        F       15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY
!        G       16     "    "     "  (TD IN T(C)) - H2O ONLY
!        H       17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY
!        I       18    AVAILABLE FOR USER DEFINITION
!        J       19    REQUEST DEFAULT TO SPECIFIED MODEL ATMOSPHERE
!
!     THIS SUBROUTINE COMPUTES THE WATERVAPOR NUMBER DENSITY (MOL CM-3)
!     GIVE HUMIDITY  # TD = DEW POINT TEMP(K,C), RH = RELATIVE
!     (PERCENT), PPH2O = WATER VAPOR PARTIAL PRESSURE (MB), DENH2O =
!     WATER VAPOR MASS DENSITY (GM M-3),AMSMIX = MASS MIXING RATIO
!     (GM/KG).
!                     THE FUNCTION DENSAT FOR THE SATURATION
!     WATER VAPOR DENSITY OVER WATER IS ACCURATE TO BETTER THAN 1
!     PERCENT FROM -50 TO +50 DEG C. (SEE THE LOWTRAN3 OR 5 REPORT)
!
!       'JUNIT' GOVERNS CHOICE OF UNITS -
!
!**********************************************************************
!
   USE phys_consts, ONLY: avogad, alosmt
   USE planet_consts, ONLY: airmwt
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
!
   DATA C1 / 18.9766 /,C2 / -14.9595 /,C3 / -2.4388 /
!
   DENSAT(ATEMP) = ATEMP*B*EXP(C1+C2*ATEMP+C3*ATEMP**2)*1.0E-6
!
   RHOAIR = ALOSMT*(P/PZERO)*(TZERO/T)
   A = TZERO/T
   B = AVOGAD/AMWT(1)
   R = AIRMWT/AMWT(1)
   IF (JUNIT.NE.10) GO TO 10
!
!     GIVEN VOL. MIXING RATIO

!     Convert using density of dry air.

   WMOL = WMOL*1.E-06
   DENNUM = (WMOL/(1.+WMOL))*RHOAIR
   GO TO 90
10 IF (JUNIT.NE.11) GO TO 20
!
!     GIVEN NUMBER DENSITY (CM-3)
!
   DENNUM = WMOL
   GO TO 90
20 CONTINUE
   IF (JUNIT.NE.12) GO TO 30
!
!     GIVEN MASS MIXING RATIO (GM KG-1)

!     Convert using density of dry air.  The following quadratic is
!
   WMOL = WMOL*R*1.0E-3
   DENNUM = (WMOL/(1.+WMOL))*RHOAIR
   GO TO 90
30 CONTINUE
   IF (JUNIT.NE.13) GO TO 40
!
!     GIVEN MASS DENSITY (GM M-3)
!
   DENNUM = B*WMOL*1.0E-6
   GO TO 90
40 CONTINUE
   IF (JUNIT.NE.14) GO TO 50
!
!     GIVEN WATER VAPOR PARTIAL PRESSURE (MB)
!
   DENNUM = ALOSMT*(WMOL/PZERO)*(TZERO/T)
   GO TO 90
50 CONTINUE
   IF (JUNIT.NE.15) GO TO 60
!
!     GIVEN DEWPOINT (DEG K)
!
   ATD = TZERO/(WMOL)
   DENNUM = DENSAT(ATD)*(WMOL)/T
   GO TO 90
60 CONTINUE
   IF (JUNIT.NE.16) GO TO 70
!
!     GIVEN DEWPOINT (DEG C)
!
   ATD = TZERO/(TZERO+WMOL)
   DENNUM = DENSAT(ATD)*(TZERO+WMOL)/T
   GO TO 90
70 CONTINUE
   IF (JUNIT.NE.17) GO TO 80
!
!     GIVEN RELATIVE HUMIDITY (PERCENT)
!
   DENNUM = DENSAT(A)*(WMOL/100.0)
   GO TO 90
80 WRITE (IPR,900) JUNIT
   STOP 'JUNIT'
90 CONTINUE
   DENST = DENSAT(A)
   RHP = 100.0*(DENNUM/DENST)
   IF (NOPRNT .ge. 0) WRITE (IPR,905) RHP
   IF (RHP.LE.100.0) GO TO 100
   if (noprnt .ge. 0) WRITE (IPR,910) RHP
100 CONTINUE
!
   RETURN
!
900 FORMAT (/,'  **** ERROR IN WATVAP ****, JUNIT = ',I5)
905 FORMAT (8X,'RH = ',F6.2)
910 FORMAT (/,' ****** WARNING (FROM WATVAP) # RELATIVE HUMIDTY = ',  &
   &        G10.3,' IS GREATER THAN 100 PERCENT')
!
end subroutine WATVAP
!
!     ----------------------------------------------------------------
!
SUBROUTINE FSCGEO (H1,H2,ANGLE,RANGE,BETA,ITYPE,LEN,HMIN,PHI,     &
&                   IERROR,HOBS)
!
!     -------------------------------------------------------------
!     This routine was modified for LBLRTM to reflect changes
!     implemented in MODTRAN to solve problems with inconsistent
!     path parameters.
!     It was also modified to eliminate GOTO statements in order to
!     make the program easier to understand.
!     These changes were obtained from H. Snell (March, 1996).
!     -------------------------------------------------------------
!
!     *****************************************************************
!     FSCGEO INTERPRETS THE ALLOWABLE COMBINATIONS OF INPUT PATH
!     PARAMETERS INTO THE STANDARD SET H1,H2,ANGLE,PHI,HMIN, AND LEN.
!     THE ALLOWABLE COMBINATIONS OF INPUT PARAMETERS ARE- FOR ITYPE = 2
!     (SLANT PATH H1 TO H2) A. H1, H2, AND ANGLE, B. H1, ANGLE, AND
!     RANGE, C. H1, H2, AND RANGE, D. H1, H2, AND BETA -
!     FOR ITYPE = 3 (SLANT PATH H1 TO SPACE, H2 = ZMAX(=100 KM,M=1 TO 6
!     A. H1 AND ANGLE, B. H1 AND HMIN (INPUT AS H2).
!     THE SUBROUTINE ALSO DETECTS BAD INPUT (IMPOSSIBLE GEOMETRY) AND
!     ITYPE = 2 CASES WHICH INTERSECT THE EARTH, AND RETURNS THESE
!     CASES WITH ERROR FLAGS.
!     THE SUBROUTINE FNDHMN IS CALLED TO CALCULATE HMIN, THE MINIMUM
!     HEIGHT ALONG THE PATH, AND PHI, THE ZENITH ANGLE AT H2, USING THE
!     ATMOSPHERIC PROFILE STORED IN /MDATA/
!     *****************************************************************
!
   USE phys_consts, ONLY: pi
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
!
   ITER = 0
!
!     Check for error
!
   IF ((ITYPE.NE.3).AND.(ITYPE.NE.2)) GOTO 90
!
   IF (ITYPE.EQ.3) THEN
!
!     Slant path to space
!     NOTE: If both HMIN and ANGLE are zero, then ANGLE is
!           assumed specified
!
      IF (H2.EQ.0) THEN
!
!             Case 3A: H1,SPACE,ANGLE
!
         WRITE (IPR,900)
         H2 = ZMAX
         CALL FNDHMN (H1,ANGLE,H2,LEN,ITER,HMIN,PHI,IERROR)

      ELSE
!
!             Case 3B: H1,HMIN,SPACE
!
         WRITE (IPR,905)
         HMIN = H2
         H2 = ZMAX
         IF (H1.LT.HMIN) GO TO 80
         CALL FNDHMN (HMIN,90.0,H1,LEN,ITER,HMIN,ANGLE,IERROR)
         CALL FNDHMN (HMIN,90.0,H2,LEN,ITER,HMIN,PHI,IERROR)
         IF (HMIN.LT.H1) LEN = 1
      ENDIF
   ENDIF
!
   IF (ITYPE.EQ.2) THEN
!
!       Assign the variable ISELCT to the following cases
!       (depending on input parameters):
!
!       -----------------------------------------------
!       H1   H2   ANGLE  RANGE  BETA  =>   CASE  ISELCT
!       -----------------------------------------------
!       X    X      X                       2A     21
!       X           X      X                2B     22
!       X    X             X                2C     23
!       X    X                   X          2D     24
!       -----------------------------------------------
!
      IF (RANGE.GT.0.0) THEN
!
!           Must be Case 2B or Case 2C
!
         IF (H2.GT.0.0) THEN
!
!              Case 2C
!
            ISELCT=23
         ELSEIF (ANGLE.EQ.0.0) THEN
            WRITE(IPR,1000)
            WRITE(*,1000)
            ISELCT=23
         ELSE
!
!              Case 2B
!
            ISELCT=22
         ENDIF
      ELSEIF (BETA.GT.0.0) THEN
!
!           Case 2D (beta cannot be zero)
!
         ISELCT=24
      ELSE
!
!           Case 2A, since RANGE and BETA are both zero
!
         ISELCT=21
      ENDIF
!
      IF (ISELCT.EQ.21) THEN
!
!           Case 2A: H1, H2, ANGLE
!
         if (noprnt .ge.0) WRITE (IPR,910)
         IF (H1.GE.H2.AND.ANGLE.LE.90.0) GO TO 110
         IF (H1.EQ.0.0.AND.ANGLE.GT.90.0) GO TO 120
         IF (H2.LT.H1.AND.ANGLE.GT.90.0) WRITE (IPR,915) LEN
         H2ST = H2
         CALL FNDHMN (H1,ANGLE,H2,LEN,ITER,HMIN,PHI,IERROR)
         IF (H2.NE.H2ST) GO TO 120
      ENDIF
!
      IF (ISELCT.EQ.22) THEN
!
!           Case 2B: H1, ANGLE, RANGE
!           Assume refraction
!
         if (noprnt .ge.0) WRITE (IPR,920)
         CALL NEWH2(H1,H2,ANGLE,RANGE,BETA,LEN,HMIN,PHI)
      ENDIF
!
      IF (ISELCT.EQ.23) THEN
!
!           Case 2C: H1, H2, RANGE
!
         if (noprnt .ge.0) WRITE (IPR,930)
         IF (ABS(H1-H2).GT.RANGE) GO TO 100
         R1 = H1+RE
         R2 = H2+RE
!
         ZARG2 = (H1**2-H2**2+RANGE**2+2.0*RE*(H1-H2)) /(2.0*R1*     &
            RANGE)
         ERARG2 = ABS(ZARG2)-1.0
         IF ((ERARG2.LE.1.0E-6).AND.(ERARG2.GE.0.0)) THEN
            IF (ZARG2.LT.0.0) THEN
               ZARG2 = -1.0
            ELSE
               ZARG2 = 1.0
            ENDIF
         ENDIF
         ANGLE = 180.0-ACOS(ZARG2)*DEG
         ZARG3 = (H2**2-H1**2+RANGE**2+2*RE*(H2-H1))/(2.0*R2*RANGE)
         ERARG3 = ABS(ZARG3)-1.0
         IF ((ERARG3.LE.1.0E-6).AND.(ERARG3.GE.0.0)) THEN
            IF (ZARG3.LT.0.0) THEN
               ZARG3 = -1.0
            ELSE
               ZARG3 = 1.0
            ENDIF
         ENDIF
         PHI = 180.0-ACOS(ZARG3)*DEG
         BETA = PHI+ANGLE-180.
!
         IF (RANGE.GT.2.0.AND.BETA.GT.0) THEN
            CALL FDBETA (H1,H2,BETA,ANGLE,PHI,LEN,HMIN,IERROR)
         ELSE
            LEN = 0
            IF (ANGLE.GT.90.0.AND.PHI.GT.90.0) LEN = 1
            CALL FNDHMN (H1,ANGLE,H2,LEN,ITER,HMIN,PHI,IERROR)
         ENDIF
      ENDIF
!
      IF (ISELCT.EQ.24) THEN
!
!        Case 2D: H1, H2, BETA
!
         CALL FDBETA (H1,H2,BETA,ANGLE,PHI,LEN,HMIN,IERROR)
      ENDIF
   ENDIF
!
!     End of allowed cases
!
!     Test IERROR and recheck LEN
!
   IF (IERROR.NE.0) RETURN
   LEN = 0
   IF (HMIN.LT.  MIN(H1,H2)) LEN = 1
!
!     Reduce path endpoints above ZMAX to ZMAX
!
   IF (HMIN.GE.ZMAX) GO TO 130
   IF (H1.GT.ZMAX.OR.H2.GT.ZMAX) CALL REDUCE (H1,H2,ANGLE,PHI,ITER)
!
!     At this point the following parameters are defined-
!         H1,H2,ANGLE,PHI,HMIN,LEN
!
!     Calculate sin(PHI) and sin(ANGLE) and output
!
   radconv = 2.*pi/360.
   sinphi = sin(radconv*phi)
   sinangle = sin(radconv*angle)
   if (noprnt .ge. 0) WRITE (IPR,935)                                &
   &                   H1,H2,ANGLE,sinangle,PHI,sinphi,HMIN,LEN

!
!     Calculate and output geometry from satellite above 120km.
!     Subtract from 180 degrees to correctly place angle in the
!     3rd quadrant.
!
   if (hobs.gt.0.) then
      if (h2.gt.h1) then
         h_toa = h2
         sintoa = sinphi
         toa_ang = phi
      else
         h_toa = h1
         sintoa = sinangle
         toa_ang = angle
      endif
      sintoa_sat = ((re+h_toa)/(re+hobs))*sintoa
      toa_sat = 180. - asin(sintoa_sat)/radconv
      sintoa_sat = sin(radconv*toa_sat)
      diffangle = toa_sat - toa_ang
      WRITE (IPR,937) hobs,toa_sat,sintoa_sat,diffangle
   endif


!
   RETURN
!
!     Error messages
!
80 CONTINUE
   WRITE (IPR,940) H1,HMIN
   GO TO 140
90 WRITE (IPR,945) ITYPE,ITYPE
   GO TO 140
100 WRITE (IPR,950) H1,H2,RANGE
   GO TO 140
110 CONTINUE
   WRITE (IPR,955) H1,H2,ANGLE
   GO TO 140
120 WRITE (IPR,960)
   GO TO 140
130 WRITE (IPR,965) ZMAX,H1,H2,HMIN
140 IERROR = 1
!
   RETURN
!
900 FORMAT (//,' CASE 3A: GIVEN H1,H2=SPACE,ANGLE')
905 FORMAT (//,' CASE 3B: GIVEN H1, HMIN, H2=SPACE')
910 FORMAT (//,' CASE 2A: GIVEN H1, H2, ANGLE')
915 FORMAT (//,' EITHER A SHORT PATH (LEN=0) OR A LONG PATH ',        &
   &        'THROUGH A TANGENT HEIGHT (LEN=1) IS POSSIBLE: LEN = ',   &
   &        I3)
920 FORMAT (//,' CASE 2B:, GIVEN H1, ANGLE, RANGE',//,10X,            &
   &        'NOTE: H2 IS COMPUTED FROM H1, ANGLE, AND RANGE ',        &
   &        'ASSUMING REFRACTION')
925 FORMAT (//,10X,'CALCULATED H2 IS LESS THAN ZERO:',/,10X,          &
   &        'RESET H2 = 0.0 AND RANGE = ',F10.3)
930 FORMAT (//,' CASE 2C: GIVEN H1, H2, RANGE',//,10X,                &
   &        'NOTE: ANGLE IS COMPUTED FROM H1, H2, AND RANGE ',        &
   &        'ASSUMING NO REFRACTION')
935 FORMAT (///,' SLANT PATH PARAMETERS IN STANDARD FORM',/           &
   &        /,10X,'H1         = ',F12.6,' KM',                        &
   &        /,10X,'H2         = ',F12.6,' KM',                        &
   &        /,10X,'ANGLE      = ',F12.6,' DEG',                       &
   &        /,10X,'sin(ANGLE) = ',F12.6,                              &
   &        /,10X,'PHI        = ',F12.6,' DEG',                       &
   &        /,10X,'sin(PHI)   = ',F12.6,                              &
   &        /,10X,'HMIN       = ',F12.6,' KM',                        &
   &        /,10X,'LEN        = ',I10)
937 FORMAT (///,' SLANT PATH PARAMETERS AT SATELLITE',/               &
   &        /,10X,'H_SAT        = ',F12.6,' KM',                      &
   &        /,10X,'PHI_SAT      = ',F12.6,' DEG'                      &
   &        /,10X,'sin(PHI_SAT) = ',F12.6,                            &
   &        /,10X,'PHI_SAT-PHI  = ',F12.6,' DEG')
940 FORMAT ('0FSCGEO: CASE 3B (H1,HMIN,SPACE): ERROR IN INPUT DATA',  &
   &        //,10X,'H1 = ',F12.6,'    IS LESS THAN HMIN = ',F12.6)
945 FORMAT ('0FSCGEO: ERROR IN INPUT DATA, ITYPE NOT EQUAL TO ',      &
   &        ' 2, OR 3.   ITYPE = ',I10,E23.14)
950 FORMAT ('0FSCGEO: CASE 2C (H1,H2,RANGE): ERROR IN INPUT DATA',    &
   &        //,10X,'ABS(H1-H2) GT RANGE;  H1 = ',F12.6,'    H2 = ',   &
   &        F12.6,'    RANGE = ',F12.6)
955 FORMAT ('0FSCGEO: CASE 2A (H1,H2,ANGLE): ERROR IN INPUT DATA',    &
   &        //,10X,'H1 = ',F12.6,'    IS GREATER THAN OR EQUAL TO',   &
   &        ' H2 = ',F12.6,/,10X,'AND ANGLE = ',F12.6,'    IS LESS',  &
   &        ' THAN OR EQUAL TO 90.0')
960 FORMAT ('0FSCGEO: ITYPE = 2: SLANT PATH INTERSECTS THE EARTH',    &
   &        ' AND CANNOT REACH H2')
965 FORMAT (' FSCGEO:  THE ENTIRE PATH LIES ABOVE THE TOP ZMAX ',     &
   &        'OF THE ATMOSPHERIC PROFILE',//,10X,'ZMAX = ',G12.6,5X,   &
   &        '  H1 = ',G12.6,5X,'  H2 = ',G12.6,'  HMIN = ',G12.6)
!
1000 FORMAT (/3X, 'Ambiguous Inputs:',/3X,'H1 and RANGE are both > 0', &
   &    /3X,'but H2 and ANGLE = 0',//3X,'Path could be 2B or 2C',     &
   &    //5X,'will assume 2C',//3X,                                   &
   &    'change in FSCGEO if 2B is desired')
end subroutine FSCGEO
!
!     ----------------------------------------------------------------
!
SUBROUTINE REDUCE (H1,H2,ANGLE,PHI,ITER)
!
!     *****************************************************************
!     ZMAX IS THE HIGHEST LEVEL IN THE ATMOSPHERIC PROFILE STORED IN
!     COMMON /MDATA/.  IF H1 AND/OR H2 ARE GREATER THAN ZMAX, THIS
!     SUBROUTINE REDUCES THEM TO ZMAX AND RESETS ANGLE AND/OR PHI
!     AS NECESSARY. THIS REDUCTION IS NECESSARY,FOR EXAMPLE FOR
!     SATELLITE ALTITUDES, BECAUSE (1) THE DENSITY PROFILES ARE
!     POORLY DEFINED ABOVE ZMAX AND (2) THE CALCULATION TIME FOR
!     PATHS ABOVE ZMAX CAN BE EXCESSIVE ( EG. FOR GEOSYNCRONOUS
!     ALTITUDES)
!     *****************************************************************
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
!
   IF (H1.LE.ZMAX.AND.H2.LE.ZMAX) RETURN
   CALL FINDSH (H1,SH,GAMMA)
   CPATH = ANDEX(H1,SH,GAMMA)*(RE+H1)*SIN(ANGLE/DEG)
   CALL FINDSH (ZMAX,SH,GAMMA)
   CZMAX = ANDEX(ZMAX,SH,GAMMA)*(RE+ZMAX)
   ANGMAX = 180.0-ASIN(CPATH/CZMAX)*DEG
   IF (H1.LE.ZMAX) GO TO 10
   H1 = ZMAX
   ANGLE = ANGMAX
10 CONTINUE
   IF (H2.LE.ZMAX) GO TO 20
   H2 = ZMAX
   PHI = ANGMAX
20 CONTINUE
   IF (ITER.EQ.0) WRITE (IPR,900) ZMAX,ANGMAX
!
   RETURN
!
900 FORMAT (///,' FROM SUBROUTINE REDUCE : ',/,10X,'ONE OR BOTH OF',  &
   &        ' H1 AND H2 ARE ABOVE THE TOP OF THE ATMOSPHERIC ',       &
   &        'PROFILE ZMAX = ',F10.3,'  AND HAVE BEEN RESET TO ZMAX.', &
   &        /,10X,'ANGLE AND/OR PHI HAVE ALSO BEEN RESET TO THE ',    &
   &        'ZENITH ANGLE AT ZMAX = ',F10.3,' DEG')
!
end subroutine REDUCE
!
!     ----------------------------------------------------------------
!
SUBROUTINE FDBETA (H1,H2,BETAS,ANGLE,PHI,LEN,HMIN,IERROR)
!
!     *****************************************************************
!     GIVEN H1,H2,AND BETA (THE EARTH CENTERED ANGLE) THIS SUBROUTINE
!     CALCULATES THE INITIAL ZENITH ANGLE AT H1 THROUGH AN ITERATIVE
!     PROCEDURE
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   REAL*8           RA,RB,SG,ANGLE1,ANGLE2,BETA,DBETA
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON /BNDRY/ ZBND(MXFSC),PBND(MXFSC),TBND(MXFSC),ALORNZ(MXFSC), &
   &               ADOPP(MXFSC),AVOIGT(MXFSC)
!
   DATA TOLRNC / 5.0E-3 /,ITERMX / 10 /,BETD / 0.04 /
   DATA ZER / 0. /
!
   BETA = BETAS
   IFLAG = 0
   IF (H1.LE.H2) THEN
      IORDER = 1
      HA = H1
      HB = H2
   ELSE
      IORDER = -1
      HA = H2
      HB = H1
   ENDIF
!
!     IF AUTOLAYERING SELECTED(IBMAX = 0) THEN SET UP DUMMY
!     LBLRTM OUTPUT LAYERS
!
   IBMSAV = IBMAX
   IF (IBMAX.EQ.0) THEN
      IBMAX = 2
      ZBND(1) = ZMIN
      ZBND(2) = ZMAX
   ENDIF
!
!     SET PARAMETER TO SUPRESS CALCULATION OF AMOUNTS
!
   IAMTB = 2
!
!     GUESS AT ANGLE, INTEGRATE TO FIND BETA, TEST FOR
!     CONVERGENCE, AND ITERATE
!     FIRST GUESS AT ANGLE: USE THE GEOMETRIC SOLUTION (NO REFRACTION)
!
   WRITE (IPR,900)
   ITER = 0
   RA = RE+HA
   RB = RE+HB
   SG = SQRT((HA-HB)**2+4.0*RA*RB*(SIN(BETA/(2.0*DEG)))**2)
   ANGLE1 = 180.0-ACOS((HA**2-HB**2+2.0*RE*(HA-HB)+SG**2)            &
   &         /(2.0*RA*SG))*DEG
   HMIN = HA
   IF (ANGLE1.GT.90.0) HMIN = RA*SIN(ANGLE1/DEG)-RE
   HMING = HMIN
   ANGLS1 = ANGLE1
   CALL FNDHMN (HA,ANGLS1,HB,LEN,ITER,HMIN,PHI,IERROR)
   LEN = 0
   IF (HMIN.LT.HA) LEN = 1
   CALL RFPATH (HA,HB,ANGLS1,PHI,LEN,HMIN,IAMTB,RANGE,BETA1,BENDNG)
   WRITE (IPR,905) ITER,ANGLS1,BETA,ZER,SG,HMING,ZER,ZER
!
!     OBTAIN DERIVATIVE
!
   SG = SQRT((HA-HB)**2+4.0*RA*RB*(SIN((BETA+BETD)/(2.0*DEG)))**2)
   ANGLEP = 180.0-ACOS((HA**2-HB**2+2.0*RE*(HA-HB)+SG**2)            &
   &         /(2.0*RA*SG))*DEG
   DANG = ANGLE1-ANGLEP
   IF (HMIN.LT.0.0) THEN
      IFLAG = 1
      HMIN = 0.0
      CALL FNDHMN (HMIN,90.0,HA,LEN,ITER,HMIN,ANGLS1,IERROR)
   ENDIF
   ITER = 1
   LEN = 0
   IF (ANGLE1.GT.90.0) LEN = 1
   CALL FNDHMN (HA,ANGLS1,HB,LEN,ITER,HMIN,PHI,IERROR)
   LEN = 0
   IF (HMIN.LT.HA) LEN = 1
   CALL RFPATH (HA,HB,ANGLS1,PHI,LEN,HMIN,IAMTB,RANGE,BETA1,BENDNG)
   DBETA = BETA-BETA1
   WRITE (IPR,905) ITER,ANGLS1,BETA1,DBETA,RANGE,HMIN,PHI,BENDNG
   IF (IFLAG.EQ.1.AND.BETA1.LT.BETA) GO TO 90
50 CONTINUE
   ANGLEP = ANGLE1-DANG
   LEN = 0
   IF (ANGLEP.GT.90.0) LEN = 1
   CALL FNDHMN (HA,ANGLEP,HB,LEN,ITER,HMIN,PHI,IERROR)
   LEN = 0
   IF (HMIN.LT.HA) LEN = 1
   CALL RFPATH (HA,HB,ANGLEP,PHI,LEN,HMIN,IAMTB,RANGE,BETAP,BENDNG)
   IF (ABS(BETA1-BETAP).LT.TOLRNC) GO TO 60
   ITER = ITER+1
   DC = BETAP-BETA1
   DERIV = -DC/BETD
   ANGLE2 = ANGLE1+(ANGLE1-ANGLEP)*(BETA-BETA1)/(BETA1-BETAP)
   ANGLS2 = ANGLE2
   LEN = 0
   IF (ANGLE2.GT.90.0) LEN = 1
   CALL FNDHMN (HA,ANGLS2,HB,LEN,ITER,HMIN,PHI,IERROR)
   LEN = 0
   IF (HMIN.LT.HA) LEN = 1
   CALL RFPATH (HA,HB,ANGLS2,PHI,LEN,HMIN,IAMTB,RANGE,BETA2,BENDNG)
   DBETA = BETA-BETA2
   WRITE (IPR,905) ITER,ANGLS2,BETA2,DBETA,RANGE,HMIN,PHI,BENDNG
   IF (BETA2.LT.BETA.AND.HMIN.LT.0.0) GO TO 90
   ANGLE1 = ANGLE2
   ANGLS1 = ANGLE1
   BETA1 = BETA2
   IF (ABS(BETA-BETA2).LT.TOLRNC) GO TO 70
   IF (ITER.GT.ITERMX) GO TO 100
   GO TO 50
60 ANGLE2 = ANGLEP
   ANGLS2 = ANGLE2
   BETA = BETAP
70 CONTINUE
   IF (HMIN.LT.0.0) GO TO 90
!
!     CONVERGED TO A SOLUTION
!
   ANGLE = ANGLE2
   BETA = BETA2
!
!     ASSIGN ANGLE AND PHI TO PROPER H1 AND H2
!
   IF (IORDER.NE.1) THEN
      TEMP = PHI
      PHI = ANGLE
      ANGLE = TEMP
   ENDIF
   IBMAX = IBMSAV
   BETAS = BETA
!
   RETURN
!
!     ERROR MESSAGES
!
90 CONTINUE
   WRITE (IPR,910)
   GO TO 110
100 CONTINUE
   WRITE (IPR,915) H1,H2,BETA,ITER,ANGLE1,BETA1,ANGLE2,BETA2
!
110 IERROR = 1
!
   RETURN
!
900 FORMAT (///,' CASE 2D: GIVEN H1, H2,  BETA:',//,                  &
   &        ' ITERATE AROUND ANGLE UNTIL BETA CONVERGES',//,          &
   &        ' ITER    ANGLE',T21,'BETA',T30,'DBETA',T40,'RANGE',      &
   &        T51,'HMIN',T61,'PHI',T70,'BENDING',/,T10,'(DEG)',T21,     &
   &        '(DEG)',T30,'(DEG)',T41,'(KM)',T51,'(KM)',T60,'(DEG)',    &
   &        T71,'(DEG)',/)
905 FORMAT (I5,3F10.4,2F10.3,2F10.4)
910 FORMAT ('0FDBETA, CASE 2D(H1,H2,BETA): REFRACTED TANGENT ',       &
   &        'HEIGHT IS LESS THAN ZERO-PATH INTERSECTS THE EARTH',     &
   &        //,10X,'BETA IS TOO LARGE FOR THIS H1 AND H2')
915 FORMAT ('0FDBETA, CASE 2D (H1,H2,BETA): SOLUTION DID NOT ',       &
   &        ' CONVERGE',//,10X,'H1 = ',F12.6,'    H2 = ',F12.6,       &
   &        '    BETA = ',F12.6,'    ITERATIONS = ',I4,//,10X,        &
   &        'LAST THREE ITERATIONS ',//,(10X,'ANGLE = ',F15.9,        &
   &        '    BETA = ',F15.9))
!
end subroutine FDBETA
!
!     ----------------------------------------------------------------
!
SUBROUTINE FNDHMN (H1,ANGLE,H2,LEN,ITER,HMIN,PHI,IERROR)
!
!     *****************************************************************
!     THIS SUBROUTINE CALCULATES THE MINIMUM ALTITUDE HMIN ALONG
!     THE REFRACTED PATH AND THE FINAL ZENITH ANGLE PHI.
!     THE PARAMETER LEN INDICATES WHETHER THE PATH GOES THROUGH
!     A TANGENT HEIGHT (LEN=1) OR NOT (LEN=0).  IF ANGLE > 90 AND
!     H1 > H2, THEN LEN CAN EITHER BE 1 OR 0, AND THE CHOICE IS
!     LEFT TO THE USER.
!     THE (INDEX OF REFRACTION - 1.0) IS MODELED AS AN EXPONENTIAL
!     BETWEEN THE LAYER BOUNDARIES, WITH A SCALE HEIGHT SH AND AN
!     AMOUNT AT THE GROUND GAMMA.
!     CPATH IS THE REFRACTIVE CONSTANT FOR THIS PATH AND
!     EQUALS  INDEX(H1)*(RE+H1)*SIN(ANGLE).
!     *****************************************************************
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
!
   REAL*8           CPATH,CRFRCT,ANDEXD,SH,GAMMA,CT1,CTP,            &
   &                 CH2,CMIN
!
   DATA DH / 0.2 /,ETA / 5.0E-7 /
!
!>    ETA MAY BE TOO SMALL FOR SOME COMPUTERS. TRY 1.0E-7 FOR 32 BIT
!>    WORD MACHINES
!
   CRFRCT(H) = (RE+H)*ANDEXD(H,SH,GAMMA)
   N = 0
   CALL FNDSHD (H1,SH,GAMMA)
   CPATH = CRFRCT(H1)*SIN(ANGLE/DEG)
   CALL FNDSHD (H2,SH,GAMMA)
   CH2 = CRFRCT(H2)
   IF (ABS(CPATH/CH2).GT.1.0) GO TO 70
   IF (ANGLE.LE.90.0) THEN
      LEN = 0
      HMIN = H1
      GO TO 60
   ENDIF
   IF (H1.LE.H2) LEN = 1
   IF (LEN.NE.1) THEN
      LEN = 0
      HMIN = H2
      GO TO 60
   ENDIF
!
!     LONG PATH THROUGH A TANGENT HEIGHT.
!     SOLVE ITERATIVELY FOR THE TANGENT HEIGHT HT.
!     HT IS THE HEIGHT FOR WHICH  INDEX(HT)*(RE+HT) = CPATH.
!
   CALL FNDSHD (0.0,SH,GAMMA)
   CMIN = CRFRCT(0.0)
!
!     FOR BETA CASES (ITER>0), ALLOW FOR HT < 0.0
!
   IF (ITER.EQ.0.AND.CPATH.LT.CMIN) GO TO 50
   HT1 = H1*SIN(ANGLE/DEG)+(SIN(ANGLE/DEG)-1.0)*RE
!
!     ITERATE TO FIND HT
!
30 CONTINUE
   N = N+1
   CALL FNDSHD (HT1,SH,GAMMA)
   CT1 = CRFRCT(HT1)
   IF (ABS((CPATH-CT1)/CPATH).LT.ETA) GO TO 40
   IF (N.GT.15) GO TO 80
   HTP = HT1-DH
   CALL FNDSHD (HTP,SH,GAMMA)
   CTP = CRFRCT(HTP)
   DERIV=(CT1-CTP)/DH
   HT1=HT1+(CPATH-CT1)/DERIV
   GO TO 30
40 CONTINUE
   HMIN=HT1
   GO TO 60
50 CONTINUE
!
!     TANGENT PATH INTERSECTS EARTH
!
   H2 = 0.0
   HMIN = 0.0
   LEN = 0
   CH2 = CMIN
   WRITE (IPR,900) H1,ANGLE
60 CONTINUE
!
!     CALCULATE THE ZENITH ANGLE PHI AT H2
!
   PHI = ASIN(CPATH/CH2)*DEG
   IF (ANGLE.LE.90.0.OR.LEN.EQ.1) PHI = 180.0-PHI
!
   RETURN
!
!     H2 LT TANGENT HEIGHT FOR THIS H1 AND ANGLE
!
70 CONTINUE
   WRITE (IPR,905)
   IERROR = 2
!
   RETURN
!
80 CONTINUE
   DC = CPATH-CT1
   WRITE (IPR,910) N,CPATH,CT1,DC,HT1
!
   STOP ' FNDHMN '
!
900 FORMAT (///,' TANGENT PATH WITH H1 = ',F10.3,' AND ANGLE = ',     &
   &        F10.3,' INTERSECTS THE EARTH',//,10X,                     &
   &        'H2 HAS BEEN RESET TO 0.0 AND LEN TO 0')
905 FORMAT ('0H2 IS LESS THAN THE TANGENT HEIGHT FOR THIS PATH ',     &
   &        'AND CANNOT BE REACHED')
910 FORMAT (///,'0FROM SUBROUTINE FNDHMN :',//,10X,                   &
   &        'THE PROCEEDURE TO FIND THE TANGENT HEIGHT DID NOT ',     &
   &        'CONVERG AFTER ',I3,'  ITERATIONS',//,10X,'CPATH   = ',   &
   &        F12.5,' KM',//,10X,'CT1     = ',F12.5,' KM',//,10X,       &
   &        'DC      = ',E12.3,' KM',//,10X,'HT1     = ',F12.5,' KM')
!
end subroutine FNDHMN
!
!     ----------------------------------------------------------------
!
SUBROUTINE FINDSH (H,SH,GAMMA)
!
!     *****************************************************************
!     GIVEN AN ALTITUDE H, THIS SUBROUTINE FINDS THE LAYER BOUNDARIES
!     Z(I1) AND Z(I2) WHICH CONTAIN H,  THEN CALCULATES THE SCALE
!     HEIGHT (SH) AND THE VALUE AT THE GROUND (GAMMA+1) FOR THE
!     REFRACTIVITY (INDEX OF REFRACTION -1)
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD) ,      &
   &       AMTP(MXMOL,MXPDIM)
!
   DO 10 IM = 2, IMMAX
      I2 = IM
      IF (ZMDL(IM).GE.H) GO TO 20
10 END DO
   I2 = IMMAX
20 CONTINUE
   I1 = I2-1
   CALL SCALHT (ZMDL(I1),ZMDL(I2),RFNDXM(I1),RFNDXM(I2),SH,GAMMA)
!
   RETURN
!
end subroutine FINDSH
!
!     ----------------------------------------------------------------
!
SUBROUTINE SCALHT (Z1,Z2,RFNDX1,RFNDX2,SH,GAMMA)
!
!     *****************************************************************
!     THIS SUBROUTINE CALCULATES THE SCALE HEIGHT SH OF THE (INDEX OF
!     REFRACTION-1.0) FROM THE VALUES OF THE INDEX AT THE ALTITUDES Z1
!     AND Z2 ( Z1 < Z2). IT ALSO CALCULATES THE EXTRAPOLATED VALUE
!     GAMMA OF THE (INDEX-1.0) AT Z = 0.0
!     *****************************************************************
!
   RF1 = RFNDX1+1.0E-20
   RF2 = RFNDX2+1.0E-20
   RATIO = RF1/RF2
   IF (ABS(RATIO-1.0).LT.1.0E-05) GO TO 10
   SH = (Z2-Z1)/ LOG(RATIO)
   GAMMA = RF1*(RF2/RF1)**(-Z1/(Z2-Z1))
   GO TO 20
10 CONTINUE
!
!     THE VARIATION IN THE INDEX OF REFRACTION WITH HEIGHT IS
!     INSIGNIFICANT OR ZERO
!
   SH = 0.0
   GAMMA = RFNDX1
20 CONTINUE
!
   RETURN
!
end subroutine SCALHT
FUNCTION ANDEX (H,SH,GAMMA)
!
!     *****************************************************************
!     COMPUTES THE INDEX OF REFRACTION AT HEIGHT H, SH IS THE
!     SCALE HEIGHT, GAMMA IS THE VALUE AT H=0 OF THE REFRACTIVITY =
!     INDEX-1
!     *****************************************************************
!
   IF (SH.EQ.0.0) THEN
      ANDEX = 1.0+GAMMA
   ELSE
      ANDEX = 1.0+GAMMA*EXP(-H/SH)
   ENDIF
!
   RETURN
!
end function ANDEX
FUNCTION RADREF (H,SH,GAMMA)
!
!     *****************************************************************
!     COMPUTES THE RADIUS OF CURVATURE OF THE REFRACTED RAY FOR
!     A HORIZONTAL PATH.  RADREF = ANDEX/ D(ANDEX)/D(RADIUS)
!     *****************************************************************
!
   DATA BIGNUM / 1.0E36 /
!
   IF (SH.EQ.0.0) GO TO 10
   RADREF = SH*(1.0+EXP(H/SH)/GAMMA)
!
   RETURN
!
10 RADREF = BIGNUM
!
   RETURN
!
end function RADREF
!
!     ----------------------------------------------------------------
!
SUBROUTINE RFPATH (H1,H2,ANGLE,PHI,LEN,HMIN,IAMT,RANGE,BETA,      &
&                   BENDNG)
!
!     -------------------------------------------------------------
!     This routine was modified for LBLRTM to reflect changes
!     implemented in MODTRAN to solve problems with inconsistent
!     path parameters.
!     It was also modified to eliminate GOTO statements in order to
!     make the program easier to understand.
!     These changes were obtained from H. Snell (March, 1996).
!     -------------------------------------------------------------
!
!     *****************************************************************
!     THIS SUBROUTINE TRACES THE REFRACTED RAY FROM H1 WITH AN
!     INITIAL ZENITH ANGLE ANGLE TO H2 WHERE THE ZENITH ANGLE IS PHI,
!     AND CALCULATES THE ABSORBER AMOUNTS (IF IAMT.EQ.1) ALONG
!     THE PATH.  IT STARTS FROM THE LOWEST POINT ALONG THE PATH
!     (THE TANGENT HEIGHT HMIN IF LEN = 1 OR HA = MIN(H1,H2) IF LEN = 0
!     AND PROCEEDS TO THE HIGHEST POINT.  BETA AND RANGE ARE THE
!     EARTH CENTERED ANGLE AND THE TOTAL DISTANCE RESPECTIVELY
!     FOR THE REFRACTED PATH FROM H1 TO H2, AND BENDNG IS THE TOTAL
!     BENDING ALONG THE PATH
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
   REAL*8           DS,DBEND,S,SINAI,COSAI,CPATH,ANDEXD,SH,GAMMA
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
   CHARACTER*2 HLOW(2)
   DATA HLOW / 'H1','H2'/
   DATA I_2/2/
!
!     REORDER H1 AND H2 TO HA AND HB (HA .LE. HB)
!
   IF (H1.LE.H2) THEN
      IORDER = 1
      HA = H1
      HB = H2
      ANGLEA = ANGLE
   ELSE
      IORDER = -1
      HA = H2
      HB = H1
      ANGLEA = PHI
   ENDIF
!
!     MERGE THE ATMOSPHERIC PROFILE STORED IN ZMDL WITH H1,H2,(HMIN) AN
!     THE BOUNDARIES ZBND
!
   CALL AMERGE (H1,H2,HMIN,LEN)
   IF (IAMT.EQ.1.AND.NOPRNT.ge.0) WRITE (IPR,900)
!
!     CALCULATE CPATH SEPERATELY FOR LEN = 0,1
!
   IF (LEN.EQ.0) THEN
      CALL FNDSHD (HA,SH,GAMMA)
      CPATH = (RE+HA)*ANDEXD(HA,SH,GAMMA)*SIN(ANGLEA/DEG)
   ELSE
      CALL FNDSHD (HMIN,SH,GAMMA)
      CPATH = (RE+HMIN)*ANDEXD(HMIN,SH,GAMMA)
   ENDIF
!
   BETA = 0.0
   S = 0.0
   BENDNG = 0.0
   IF (LEN.EQ.1) THEN
!
!     TANGENT PATH
!
      IF (IORDER.EQ.-1) THEN
         IHLOW = 2
      ELSE
         IHLOW = 1
      ENDIF
      IF (IAMT.EQ.1.AND.NOPRNT.ge.0) WRITE (IPR,905) HLOW(IHLOW)
      SINAI = 1.0
      COSAI = 0.0
      THETA = 90.0
   ELSE
!
!     SHORT PATH
!
!     ANGLEA IS THE ZENITH ANGLE AT HA IN DEG
!     SINAI IS SIN OF THE INCIDENCE ANGLE
!     COSAI IS CARRIED SEPERATELY TO AVOID A PRECISION PROBLEM
!     WHEN SINAI IS CLOSE TO 1.0
!
      THETA = ANGLEA
      IF (ANGLEA.LE.45.0) THEN
         SINAI = SIN(ANGLEA/DEG)
         COSAI = -COS(ANGLEA/DEG)
      ELSE
         SINAI = COS((90.0-ANGLEA)/DEG)
         COSAI = -SIN((90.0-ANGLEA)/DEG)
      ENDIF
      IF (IORDER.EQ.-1) THEN
         IHLOW = 2
      ELSE
         IHLOW = 1
      ENDIF
      IHIGH = MOD(IHLOW,I_2)+1
      IF (IAMT.EQ.1.AND.NOPRNT.ge.0) WRITE (IPR,910) HLOW(IHLOW),    &
         HLOW(IHIGH)
   ENDIF
!
!     LOOP OVER THE LAYERS
!
   J2 = IPMAX-1
   DO 100 J = 1, J2
      CALL SCLHTD (ZPTH(J),ZPTH(J+1),RFNDXP(J),RFNDXP(J+1),SH,GAMMA)
      CALL ALAYER (J,SINAI,COSAI,CPATH,SH,GAMMA,IAMT,DS,DBEND)
      DBEND = DBEND*DEG
      PHI = ASIN(SINAI)*DEG
      DBETA = THETA-PHI+DBEND
      PHI = 180.0-PHI
      S = S+DS
      BENDNG = BENDNG+DBEND
      BETA = BETA+DBETA
      IF (IAMT.EQ.1) THEN
         PBAR = PPSUM(J)/RHOPSM(J)
         TBAR = TPSUM(J)/RHOPSM(J)
         RHOBAR = RHOPSM(J)/DS
         IF (NOPRNT.ge.0) WRITE (IPR,915) J,ZPTH(J),ZPTH(J+1),       &
            THETA,DS,S,DBETA,BETA,PHI,DBEND,BENDNG,PBAR, TBAR,RHOBAR
      ENDIF
      THETA = 180.0-PHI
!
      IF (LEN.EQ.1) THEN
!
!            For tangent paths, double the quantities BENDNG,BETA,
!            and S for the symmetric part of the path
!
         IF ((J+1).EQ.IPHMID) THEN
            BENDNG = 2.0*BENDNG
            BETA = 2.0*BETA
            S = 2.0*S
            IF (IAMT.EQ.1.AND.NOPRNT.ge.0) WRITE (IPR,920) S,BETA,   &
               BENDNG
            IF (IPHMID.NE.IPMAX) THEN
               IF (IORDER.EQ.-1) THEN
                  IHLOW = 2
               ELSE
                  IHLOW = 1
               ENDIF
               IHIGH = MOD(IHLOW,I_2)+1
               IF (IAMT.EQ.1.AND.NOPRNT.ge.0) WRITE (IPR,910) HLOW(  &
                  IHLOW),HLOW(IHIGH)
            ENDIF
         ENDIF
      ENDIF
100 END DO
   IF (IORDER.EQ.-1) PHI = ANGLEA
   RANGE = S
!
   RETURN
!
900 FORMAT ('1CALCULATION OF THE REFRACTED PATH THROUGH THE ',        &
   &        'ATMOSPHERE',///,T5,'I',T14,'ALTITUDE',T30,'THETA',T38,   &
   &        'DRANGE',T47,'RANGE',T57,'DBETA',T65,'BETA',T76,'PHI',    &
   &        T84,'DBEND',T91,'BENDING',T102,'PBAR',T111,'TBAR',T119,   &
   &        'RHOBAR',/,T11,'FROM',T22,'TO',/,T11,'(KM)',T21,'(KM)',   &
   &        T30,'(DEG)',T39,'(KM)',T48,'(KM)',T57,'(DEG)',T65,        &
   &        '(DEG)',T75,'(DEG)',T84,'(DEG)',T92,'(DEG)',T102,'(MB)',  &
   &        T112,'(K)',T117,'(MOL CM-3)',/)
905 FORMAT (' ',T10,'TANGENT',T20,A2,/,T10,'HEIGHT',/)
910 FORMAT (' ',T14,A2,' TO ',A2,/)
915 FORMAT (' ',I4,2F10.3,10F9.3,1PE9.2)
920 FORMAT ('0',T10,'DOUBLE RANGE, BETA, BENDING',/,T10,              &
   &        'FOR SYMMETRIC PART OF PATH',T44,F9.3,T62,F9.3,T89,       &
   &        F9.3,/)
!
end subroutine RFPATH
!
!     ----------------------------------------------------------------
!
SUBROUTINE AMERGE (H1,H2,HMIN,LEN)
!
!     *****************************************************************
!     AMERGE CREATES A SET OF LAYER BOUNDARIES ZOUT WHICH INCLUDES
!     HMIN, (HMID), HMAX AND ALL OF ZBND BETWEEN HMIN AND HAMX.
!     ZOUT DEFINES THE LAYERS FOR THE LBLRTM CALCULATION.
!     ZOUT IS THEN MERGED WITH THE ATMOSPHERIC PROFILE IN ZMDL INTO ZPT
!     INTERPOLATING TO THE LEVELS ZOUT WHEN NECESSARY.  THE RAY
!     TRACE IS CALCULATED USING THE PROFILE IN ZPTH.
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
   COMMON /DEAMT/ DENM(MXMOL,MXZMD),DENP(MXMOL,MXPDIM),DRYAIR(MXZMD)
   COMMON /BNDRY/ ZBND(MXFSC),PBND(MXFSC),TBND(MXFSC),ALORNZ(MXFSC), &
   &               ADOPP(MXFSC),AVOIGT(MXFSC)
   COMMON /ZOUTP/ ZOUT(MXLAY),SOUT(MXLAY),RHOSUM(MXLAY),             &
   &               AMTTOT(MXMOL),AMTCUM(MXMOL),ISKIP(MXMOL)
   DIMENSION ZH(3)
!
   DATA TOL / 5.E-4 /
!
   DATA I_2/2/

!
!     HMID .EQ. MINIMUM OF H1, H2
!
   HMID =   MIN(H1,H2)
   HMAX =   MAX(H1,H2)
   IHMAX = 2
   ZH(1) = HMIN
   IF (LEN.EQ.0) THEN
      ZH(2) = HMAX
   ELSE
      ZH(2) = HMID
      IF (ABS(H1-H2).LT.TOL) H1 = H2
      IF (H1.NE.H2) THEN
         IHMAX = 3
         ZH(3) = HMAX
      ENDIF
   ENDIF
!
!     MERGE ZH AND ZBND BETWEEN ZH(1) AND ZH(IHMAX) TO CREAT ZOUT
!
   ZOUT(1) = ZH(1)
   DO 30 I1 = 1, IBMAX
      IF (ABS(ZBND(I1)-ZH(1)).LT.TOL) ZH(1) = ZBND(I1)
      IF (ZBND(I1).GT.ZH(1)) GO TO 40
30 END DO
   I1 = IBMAX
40 CONTINUE
!
!     ZBND(I1) IS SMALLEST ZBND .GT. ZH(1)
!
   IOUT = 1
   IB = I1
   IH = 2
50 CONTINUE
   IOUT = IOUT+1
   IF (IB.GT.IBMAX) GO TO 60
   IF (ABS(ZBND(IB)-ZH(IH)).LT.TOL) ZH(IH) = ZBND(IB)
   IF (ZBND(IB).LT.ZH(IH)) GO TO 70
   IF (ZBND(IB).EQ.ZH(IH)) IB = IB+1
!
!     INSERT ZH(IH)
!
60 CONTINUE
   ZOUT(IOUT) = ZH(IH)
   IH = IH+1
   IF (IH.GT.IHMAX) GO TO 80
   GO TO 50
!
!     INSERT ZBND(IB)
!
70 CONTINUE
   ZOUT(IOUT) = ZBND(IB)
   IB = IB+1
   GO TO 50
!
80 CONTINUE
   IOUTMX = IOUT
!
!     NOW MERGE ZOUT AND ZMDL INTO ZPTH (FROM ZOUT(1) TO ZOUT(IOUTMX))
!     AND INTERPOLATE PRESSURE, TEMPERATURE, AND DENSITY WHEN
!     NECESSARY
!
!     FIND SMALLEST ZMDL .GT. HMIN
!
   DO 90 IM = 1, IMMAX
      IF (ZMDL(IM).GE.HMIN) GO TO 100
90 END DO
   WRITE (IPR,900) HMIN
   STOP ' AMERGE - HMIN '
100 CONTINUE
   IPHMID = 0
   IP = 0
   IOUT = 1
110 CONTINUE
   IP = IP+1
   IF (IP.GT.IPDIM) THEN
      WRITE (IPR,905) IPDIM
      STOP ' AMERGE - IPDIM '
   ENDIF
   IF (IM.GT.IMMAX) GO TO 130
   IF (ABS(ZOUT(IOUT)-ZMDL(IM)).LT.TOL) ZMDL(IM) = ZOUT(IOUT)
   IF (ZOUT(IOUT).LT.ZMDL(IM)) GO TO 130
   IF (ZOUT(IOUT).EQ.ZMDL(IM)) IOUT = IOUT+1
!
!     INSERT ZMDL(IM)
!
   ZPTH(IP) = ZMDL(IM)
   PP(IP) = PM(IM)
   TP(IP) = TM(IM)
   RFNDXP(IP) = RFNDXM(IM)
   DO 120 K = 1, NMOL
      DENP(K,IP) = DENM(K,IM)
120 END DO
   IM = IM+1
   IF (ABS(ZPTH(IP)-HMID).LT.TOL) HMID = ZPTH(IP)
   IF (ZPTH(IP).EQ.HMID) IPHMID = IP
   IF (ABS(ZPTH(IP)-ZOUT(IOUTMX)).LT.TOL) ZOUT(IOUTMX) = ZPTH(IP)
   IF (ZPTH(IP).EQ.ZOUT(IOUTMX)) GO TO 150
   GO TO 110
!
!     INSERT LEVEL FROM ZOUT(IOUT) AND INTERPOLATE
!
130 CONTINUE
   ZPTH(IP) = ZOUT(IOUT)
   JM = IM
   JM = MAX(JM,I_2)
   A = (ZOUT(IOUT)-ZMDL(JM-1))/(ZMDL(JM)-ZMDL(JM-1))
   CALL EXPINT (PP(IP),PM(JM-1),PM(JM),A)
   TP(IP) = TM(JM-1)+(TM(JM)-TM(JM-1))*A
   CALL EXPINT (RFNDXP(IP),RFNDXM(JM-1),RFNDXM(JM),A)
   DO 140 K = 1, NMOL
      CALL EXPINT (DENP(K,IP),DENM(K,JM-1),DENM(K,JM),A)
140 END DO
   IF (ABS(ZPTH(IP)-HMID).LT.TOL) ZPTH(IP) = HMID
   IF (ZPTH(IP).EQ.HMID) IPHMID = IP
   IOUT = IOUT+1
   IF (ABS(ZPTH(IP)-ZOUT(IOUTMX)).LT.TOL) ZPTH(IP) = ZOUT(IOUTMX)
   IF (ZPTH(IP).EQ.ZOUT(IOUTMX)) GO TO 150
   GO TO 110
150 CONTINUE
   IPMAX = IP
!
   RETURN
!
900 FORMAT ('0FROM AMERGE- ATMOSPHERIC PROFILE IN ZMDL DOES NOT',     &
   &        ' EXTEND UP TO HMIN = ',E12.5)
905 FORMAT ('0FROM AMERGE- MERGING THE ATMOSPHERIC PROFILE AND THE ', &
   &        'LBLRTM BOUNDARIES INTO ZPTH(IPDIM) EXCEEDS THE ',        &
   &        'DIMENSION IPDIM = ',I5)
!
end subroutine AMERGE
!
!     ----------------------------------------------------------------
!
SUBROUTINE ALAYER (J,SINAI,COSAI,CPATH,SH,GAMMA,IAMT,S,BEND)
!
!     -------------------------------------------------------------
!     This routine was modified for LBLRTM to reflect changes
!     implemented in MODTRAN to solve problems with inconsistent
!     path parameters.
!     It was also modified to eliminate GOTO statements in order to
!     make the program easier to understand.
!     These changes were obtained from H. Snell (March, 1996).
!     -------------------------------------------------------------
!
!     *****************************************************************
!     THIS SUBROUTINE TRACES THE OPTICAL RAY THROUGH ONE LAYER FROM
!     Z1 TO Z2 AND IF IAMT.NE.2 CALCULATES THE INTEGRATED ABSORBER
!     AMOUNTS FOR THE LAYER. SINAI IS THE SIN OF THE INITIAL INCIDENCE
!     ANGLE (= 180 - ZENITH ANGLE). COSAI IS CARRIED SEPERATELY TO
!     AVOID A PRECISION PROBLEM NEAR SINAI = 1. CPATH IS THE CONSTANT
!     OF REFRACTION FOR THE PATH = INDEX*RADIUS*SINAI, SH AND GAMMA ARE
!     THE SCALE HEIGHT AND THE AMOUNT AT THE GROUND FOR THE REFRACTIVIT
!     (= 1-INDEX OF REFRACTION), S IS THE REFRACTED PATH LENGTH THROUGH
!     THE LAYER, BETA IS THE EARTH CENTERED ANGLE, AND BEND IS THE
!     BENDING THROUGH THE LAYER. IAMT CONTROLS WHETHER AMOUNTS ARE
!     CALCULATED OR NOT.
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
   COMMON /DEAMT/ DENM(MXMOL,MXZMD),DENP(MXMOL,MXPDIM),DRYAIR(MXZMD)
   DIMENSION HDEN(MXMOL),DENA(MXMOL),DENB(MXMOL)
!
   REAL*8           S,BEND,DS,DBEND,W1,W2,W3,DSDX1,DSDX2,DSDX3,      &
   &       DBNDX1,DBNDX2,DBNDX3,R1,R2,R3,X1,X2,X3,RATIO1,RATIO2,      &
   &       RATIO3,SINAI1,SINAI2,SINAI3,COSAI1,COSAI2,COSAI3,Y1,Y3,    &
   &       CPATH,DX,DH,SINAI,COSAI,D31,D32,D21,DHMIN,                 &
   &       SH,GAMMA,ANDEXD,RADRFD
!
   DATA EPSILN / 1.0E-5 /
!
!     INITIALIZE VARIABLES FOR THE CALCULATION OF THE PATH
!
   N = 0
   Z1 = ZPTH(J)
   Z2 = ZPTH(J+1)
   H1 = Z1
   R1 = RE+H1
   DHMIN = DELTAS**2/(2.0*R1)
   SINAI1 = SINAI
   COSAI1 = COSAI
   IF ((1.0-SINAI).LT.EPSILN)                                        &
   &     Y1 = COSAI1**2/2.0+COSAI1**4/8.0+COSAI1**6*3.0/48.0
   Y3 = 0.0
   X1 = -R1*COSAI1
   RATIO1 = R1/RADRFD(H1,SH,GAMMA)
   DSDX1 = 1.0/(1.0-RATIO1*SINAI1**2)
   DBNDX1 = DSDX1*SINAI1*RATIO1/R1
   S = 0.0
   BEND = 0.0
   IF (IAMT.NE.2) THEN
!
!         Initialize the variables for the calculation of the
!         absorber amounts
!
      PA = PP(J)
      PB = PP(J+1)
      IF (PB.EQ.PA) THEN
         WRITE(*,*) PB
         STOP 'LBLATM: PRESSURES IN ADJOINING LAYERS MUST DIFFER'
      ENDIF
      TA = TP(J)
      TB = TP(J+1)
      RHOA = PA/(GCAIR*TA)
      RHOB = PB/(GCAIR*TB)
      DZ = ZPTH(J+1)-ZPTH(J)
      HP = -DZ/ LOG(PB/PA)
      IF (ABS(RHOB/RHOA-1.0).GE.EPSILN) THEN
         HRHO = -DZ/ LOG(RHOB/RHOA)
      ELSE
         HRHO = 1.0E30
      ENDIF
      DO 40 K = 1, NMOL
         DENA(K) = DENP(K,J)
         DENB(K) = DENP(K,J+1)
         IF ((DENA(K).EQ.0.0.OR.DENB(K).EQ.0.0).OR. (ABS(1.0-DENA(K)/&
            DENB(K)).LE.EPSILN)) THEN
!
!                 Use linear interpolation
!
            HDEN(K) = 0.0
         ELSE
!
!                 Use exponential interpolation
!
            HDEN(K) = -DZ/ LOG(DENB(K)/DENA(K))
         ENDIF
40    CONTINUE
   ENDIF
!
!     LOOP THROUGH PATH
!     INTEGRATE PATH QUANTITIES USING QUADRATIC INTEGRATION WITH
!     UNEQUALLY SPACED POINTS
!
60 CONTINUE
   N = N+1
   DH = -DELTAS*COSAI1
   DH = MAX(DH,DHMIN)
   H3 = H1+DH
   IF (H3.GT.Z2) H3 = Z2
   DH = H3-H1
   R3 = RE+H3
   H2 = H1+DH/2.0
   R2 = RE+H2
   SINAI2 = CPATH/(ANDEXD(H2,SH,GAMMA)*R2)
   SINAI3 = CPATH/(ANDEXD(H3,SH,GAMMA)*R3)
   RATIO2 = R2/RADRFD(H2,SH,GAMMA)
   RATIO3 = R3/RADRFD(H3,SH,GAMMA)
   IF ((1.0-SINAI2).LE.EPSILN) THEN
!
!        Near a tangent height, COSAI = -SQRT(1-SINAI**2) loses
!        precision. use the following algorithm to get COSAI.
!
      Y3 = Y1+(SINAI1*(1.0-RATIO1)/R1+4.0*SINAI2*(1.0-RATIO2)/R2+    &
         SINAI3*(1.0-RATIO3)/R3)*DH/6.0
      COSAI3 = -SQRT(2.0*Y3-Y3**2)
      X3 = -R3*COSAI3
      DX = X3-X1
      W1 = 0.5*DX
      W2 = 0.0
      W3 = 0.5*DX
   ELSE
      COSAI2 = -SQRT(1.0-SINAI2**2)
      COSAI3 = -SQRT(1.0-SINAI3**2)
      X2 = -R2*COSAI2
      X3 = -R3*COSAI3
!
!        Calculate weights
!
      D31 = X3-X1
      D32 = X3-X2
      D21 = X2-X1
      IF (D32.EQ.0.0.OR.D21.EQ.0.0) THEN
         W1 = 0.5*D31
         W2 = 0.0
         W3 = 0.5*D31
      ELSE
         W1 = (2.0-D32/D21)*D31/6.0
         W2 = D31**3/(D32*D21*6.0)
         W3 = (2.0-D21/D32)*D31/6.0
      ENDIF
   ENDIF
   DSDX2 = 1.0/(1.0-RATIO2*SINAI2**2)
   DSDX3 = 1.0/(1.0-RATIO3*SINAI3**2)
   DBNDX2 = DSDX2*SINAI2*RATIO2/R2
   DBNDX3 = DSDX3*SINAI3*RATIO3/R3
!
!     INTEGRATE
!
   DS = W1*DSDX1+W2*DSDX2+W3*DSDX3
   S = S+DS
   DBEND = W1*DBNDX1+W2*DBNDX2+W3*DBNDX3
   BEND = BEND+DBEND
   IF (IAMT.NE.2) THEN
!
!         Calculate amounts
!
      DSDZ = DS/DH
      PB = PA*EXP(-DH/HP)
      RHOB = RHOA*EXP(-DH/HRHO)
      IF ((DH/HRHO).GE.EPSILN) THEN
         PPSUM(J) = PPSUM(J)+DSDZ*(HP/(1.0+HP/HRHO)) *(PA*RHOA-PB*   &
            RHOB)
         TPSUM(J) = TPSUM(J)+DSDZ*HP*(PA-PB)/GCAIR
         RHOPSM(J) = RHOPSM(J)+DSDZ*HRHO*(RHOA-RHOB)
      ELSE
         PPSUM(J) = PPSUM(J)+0.5*DS*(PA*RHOA+PB*RHOB)
         TPSUM(J) = TPSUM(J)+0.5*DS*(PA+PB)/GCAIR
         RHOPSM(J) = RHOPSM(J)+0.5*DS*(RHOA+RHOB)
      ENDIF
      DO 130 K = 1, NMOL
         IF ((HDEN(K).EQ.0.0).OR. (ABS(DH/HDEN(K)).LT.EPSILN)) THEN
!
!                 Linear interpolation
!                 1.0E05 factor converts units km to cm
!
            DENB(K)=DENP(K,J)+(DENP(K,J+1)-DENP(K,J))*(H3-Z1)/DZ
            AMTP(K,J) = AMTP(K,J)+0.5*(DENA(K)+DENB(K))*DS*1.0E5
         ELSE
!
!                 Exponential interpolation
!
            DENB(K) = DENP(K,J)*EXP(-(H3-Z1)/HDEN(K))
            AMTP(K,J) = AMTP(K,J)+DSDZ*HDEN(K) *(DENA(K)-DENB(K))*   &
               1.0E5
         ENDIF
130   CONTINUE
      PA = PB
      RHOA = RHOB
      DO 140 K = 1, NMOL
         DENA(K) = DENB(K)
140   CONTINUE
   ENDIF
!
   IF (H3.LT.Z2) THEN
      H1 = H3
      R1 = R3
      SINAI1 = SINAI3
      RATIO1 = RATIO3
      Y1 = Y3
      COSAI1 = COSAI3
      X1 = X3
      DSDX1 = DSDX3
      DBNDX1 = DBNDX3
   ELSE
      SINAI = SINAI3
      COSAI = COSAI3
      SP(J) = S
      RETURN
   ENDIF
!
   GO TO 60
!
end subroutine ALAYER
!
!     ----------------------------------------------------------------
!
SUBROUTINE AUTLAY (HMIN,HMAX,XVBAR,AVTRAT,TDIFF1,TDIFF2,ALTD1,    &
&                   ALTD2,IERROR)
!
!     *****************************************************************
!     THIS SUBROUTINE AUTOMATICALLY SELECTS A SET OF LBLRTM BOUNDARY
!     LEVELS WHICH SATISFY THE FOLLOWING TWO TESTS:
!          1. THE RATIO OF THE VOIGT HALFWIDTHS BETWEEN BOUNDARIES
!             IS LESS THAN OR EQUAL TO AVTRAT, AND
!          2. THE TEMPERATURE DIFFERENCE BETWEEN BOUNDARIES IS
!             LESS THAN OR EQUAL TO TDIFF
!     TDIFF VARIES FROM TDIFF1 AT HMIN TO TDIFF2 AT HMAX,
!     WITH EXPONENTIAL INTERPOLATION BETWEEN
!     THESE BOUNDARIES ARE ROUNDED DOWN TO THE NEAREST TENTH KM
!     NOTE THAT THESE TESTS APPLY TO THE LAYER BOUNDARIES
!     NOT TO THE AVERAGE VALUES FROM ONE LAYER TO THE NEXT.
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
   COMMON /BNDRY/ ZBND(MXFSC),PBND(MXFSC),TBND(MXFSC),ALORNZ(MXFSC), &
   &               ADOPP(MXFSC),AVOIGT(MXFSC)
   DIMENSION AVTM(MXZMD)
!
!     FUNCTION ZROUND ROUNDS THE ALTITUDE Z DOWN TO THE
!     NEAREST TENTH KM
!
   ZROUND(ZX) = 0.1* REAL( INT(10.0*ZX))
   HMIN = MAX(HMIN,ZMDL(1))
!
   DO 10 IM = 2, IMMAX
      IHMIN = IM
      IF (ZMDL(IM).GT.HMIN) GO TO 20
10 END DO
20 CONTINUE
   HTOP = HMAX
   HTOP = MIN(HTOP,ZMAX)
   IM = IHMIN-1
   ZZ = ZMDL(IM)
   CALL HALFWD (ZZ,XVBAR,P,T,AL,AD,AVTM(IM))
   IB = 1
   ZBND(IB) = HMIN
   IM = IHMIN
   CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB),         &
   &             ADOPP(IB),AVOIGT(IB))
!
!     BEGIN IM LOOP
!
30 CONTINUE
   IB = IB+1
   IF (IB.GT.IBDIM) GO TO 90
   IBM1 = IB-1
   TMIN = TBND(IBM1)
   TMAX = TBND(IBM1)
   IND = 0
!
!     BEGIN IB LOOP
!
40 CONTINUE
   IPASS = 0
   ZBND(IB) = ZMDL(IM)
   ZBNDTI = ZMDL(IM)
   IF (ZBND(IB).GE.HTOP) ZBND(IB) = HTOP
   CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB),         &
   &             ADOPP(IB),AVOIGT(IB))
   AVTM(IM) = AVOIGT(IB)
!
!     TEST THE RATIO OF THE VOIGT WIDTHS AGAINST AVTRAT
!
   IF ((AVOIGT(IB-1)/AVOIGT(IB)).LT.AVTRAT) GO TO 50
!
!     ZMDL(IM) FAILS THE HALFWIDTH RATIO TEST
!
   IPASS = 1
   AVOIGT(IB) = AVOIGT(IB-1)/AVTRAT
   X = AVTM(IM)/AVTM(IM-1)
   ALOGX = 1.-X
   IF (ABS(ALOGX).LT.0.001) THEN
      ZBND(IB) = (ZMDL(IM)+ZMDL(IM-1))/2.
      GO TO 50
   ELSE
      ALOGX = LOG(X)
   ENDIF
   Y = AVOIGT(IB)/AVTM(IM-1)
   ALOGY = 1.-Y
   IF (ABS(ALOGY).GT.0.001) ALOGY =  LOG(Y)
   ZBND(IB) = ZMDL(IM-1)+(ZMDL(IM)-ZMDL(IM-1))*ALOGY/ALOGX
50 CONTINUE
!
!     TEST THE TEMPERATURE DIFFERENCE AGAINST TDIFF
!
   FAC = (ZBND(IB-1)-ALTD1)/(ALTD2-ALTD1)
   CALL EXPINT (TDIFF,TDIFF1,TDIFF2,FAC)
   IF (TM(IM).GT.TMAX) THEN
      IND = 1
      TMAX = TM(IM)
   ENDIF
   IF (TM(IM).LT.TMIN) THEN
      IND = 2
      TMIN = TM(IM)
   ENDIF
   IF (TMAX-TMIN.LE.TDIFF) GO TO 60
   IF (IND.EQ.1) TBND(IB) = TMIN+TDIFF
   IF (IND.EQ.2) TBND(IB) = TMAX-TDIFF
!
!     ZBND(IB) FAILS THE TEMPERATURE DIFFERENCE TEST
!
   IPASS = 2
   IF (ABS(TM(IM)-TM(IM-1)).LT.0.0001) THEN
      ZBNDTI = (ZMDL(IM)+ZMDL(IM-1))/2.
   ELSE
      ZBNDTI = ZMDL(IM-1)+(ZMDL(IM)-ZMDL(IM-1))* (TBND(IB)-TM(IM-1))/&
         (TM(IM)-TM(IM-1))
   ENDIF
60 CONTINUE
   IF (ZBNDTI.LT.ZBND(IB)) ZBND(IB) = ZBNDTI
   IF (ZBND(IB).GE.HTOP) THEN
      ZBND(IB) = HTOP
      IF (ZBND(IB)-ZBND(IB-1).LE.0.1) THEN
         IB = IB-1
         ZBND(IB) = HTOP
         CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB),   &
            ADOPP(IB),AVOIGT(IB))
      ENDIF
      GO TO 80
   ENDIF
   IF (IPASS.NE.0) GO TO 70
!
!     BOTH HALFWIDTH AND TEMPERATURE TEST PASS FOR ZBND(IB) = ZMDL(IM),
!     NOW TRY ZBND(IB) = ZMDL(IM+1)
!
   IM = IM+1
   GO TO 40
70 CONTINUE
!
!     ONE OF THE TESTS FAILED AND A NEW BOUNDRY ZBND WAS PRODUCED
!
   ZBND(IB) = ZROUND(ZBND(IB))
   CALL HALFWD (ZBND(IB),XVBAR,PBND(IB),TBND(IB),ALORNZ(IB),         &
   &             ADOPP(IB),AVOIGT(IB))
   GO TO 30
80 CONTINUE
   IBMAX = IB
   WRITE (IPR,900) AVTRAT,TDIFF1,HMIN,TDIFF2,HMAX
!
   RETURN
!
90 CONTINUE
   WRITE (IPR,905) IBDIM
   IBMAX = IBDIM
   IERROR = 5
!
   RETURN
!
900 FORMAT (///,                                                      &
   &        ' LBLRTM LAYER BOUNDARIES PRODUCED BY THE AUTOMATIC ',    &
   &        'LAYERING ROUTINE AUTLAY',/,' THE USER SHOULD EXAMINE ',  &
   &        'THESE BOUNDARIES AND MODIFY THEM IF APPROPRIATE',/,      &
   &        ' THE FOLLOWING PARAMETERS ARE USED:',//,10X,             &
   &        'AVTRAT    = ',F8.2,'       = MAX RATIO OF VOIGT WIDTHS', &
   &        /,10X,'TDIFF1    = ',F8.2,'       = MAX TEMP DIFF AT ',   &
   &        F4.0,' KM',/10X,'TDIFF2    = ',F8.2,                      &
   &        '       = MAX TEMP DIFF AT ',F4.0,' KM')
905 FORMAT (///,' ERROR IN AUTLAY:',/,5X,'THE NUMBER OF ',            &
   &        'GENERATED LAYER BOUNDARIES EXCEEDS THE DIMENSION IBDIM', &
   &        ' OF THE ARRAY ZBND.  IBDIM = ',I5,/,5X,'PROBABLE CAUSE', &
   &        ': EITHER AVTRAT AND/OF TDIFF ARE TOO SMALL',/,5X,        &
   &        'THE GENERATED LAYERS FOLLOW')
!
end subroutine AUTLAY
!
!     ----------------------------------------------------------------
!
SUBROUTINE HALFWD_P (Z,XVBAR,P,T,ALORNZ,ADOPP,AVOIGT)
!
!     *****************************************************************
!     GIVEN AN PRESSURE AND TEMP. AND AVERAGE WAVENUMBER VBAR, THIS
!     SUBROUTINE
!     CALCULATES THE LORENTZ, THE DOPPLER, AND THE VOIGT HALFWIDTHS
!     (AT HALFHEIGHT) ALORNZ, ADOPP, AND AVOIGT RESPECTIVELY FOR
!     THE ALTITUDE Z
!     AN AVERAGE LORENTZ WIDTH ALZERO AND AN AVERAGE MOLECULAR
!     WEIGHT AVMWT ARE ASSUMED
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
!     FUNCTIONS
!     ALZERO IS AT 1013.25 MB AND 296.0 K
!
   ALPHAL(P,T) = ALZERO*(P/PZERO)*SQRT(296.0/T)
   ALPHAD(T,V) = ADCON*V*SQRT(T/AVMWT)
   ALPHAV(AL,AD) = 0.5*(AL+SQRT(AL**2+4.0*AD**2))
!
   ALORNZ = ALPHAL(P,T)
   ADOPP = ALPHAD(T,XVBAR)
   AVOIGT = ALPHAV(ALORNZ,ADOPP)
!
   RETURN
!
end subroutine HALFWD_P
!
!     ----------------------------------------------------------------
!
!
!     ----------------------------------------------------------------
!
SUBROUTINE HALFWD (Z,XVBAR,P,T,ALORNZ,ADOPP,AVOIGT)
!
!     *****************************************************************
!     GIVEN AN ALTITUDE Z AND AN AVERAGE WAVENUMBER VBAR, THIS
!     SUBROUTINE INTERPOLATES P AND T FROM THE PROFILE IN ZMDL  AND
!     CALCULATES THE LORENTZ, THE DOPPLER, AND THE VOIGT HALFWIDTHS
!     (AT HALFHEIGHT) ALORNZ, ADOPP, AND AVOIGT RESPECTIVELY FOR
!     THE ALTITUDE Z
!     AN AVERAGE LORENTZ WIDTH ALZERO AND AN AVERAGE MOLECULAR
!     WEIGHT AVMWT ARE ASSUMED
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
!     FUNCTIONS
!     ALZERO IS AT 1013.25 MB AND 296.0 K
!
   ALPHAL(P,T) = ALZERO*(P/PZERO)*SQRT(296.0/T)
   ALPHAD(T,V) = ADCON*V*SQRT(T/AVMWT)
   ALPHAV(AL,AD) = 0.5*(AL+SQRT(AL**2+4.0*AD**2))
!
   DO 10 I2 = 2, IMMAX
      IM = I2
      IF (ZMDL(IM).GE.Z) GO TO 20
10 END DO
   IM = IMMAX
20 CONTINUE
   FAC = (Z-ZMDL(IM-1))/(ZMDL(IM)-ZMDL(IM-1))
   CALL EXPINT (P,PM(IM-1),PM(IM),FAC)
   T = TM(IM-1)+(TM(IM)-TM(IM-1))*FAC
   ALORNZ = ALPHAL(P,T)
   ADOPP = ALPHAD(T,XVBAR)
   AVOIGT = ALPHAV(ALORNZ,ADOPP)
!
   RETURN
!
end subroutine HALFWD
!
!     ----------------------------------------------------------------
!
SUBROUTINE FPACK (H1,H2,HMID,LEN,IEMIT,n_zero)
!
!     *****************************************************************
!     FPACK TAKES THE AMOUNTS STORED IN THE LAYERS DEFINED BY ZPTH AND
!     PACKS THEM INTO THE LAYERS DEFINED BY ZOUT.  IT ALSO ZEROS OUT
!     LAYER AMOUNTS IF THE AMOUNT FOR THAT LAYER AND ABOVE IS LESS
!     THAN 0.1 PERCENT OF THE TOTAL FOR THAT MOLECULE if THE
!     n_zero OPTION IS SELECTED.
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
   IMPLICIT REAL*8           (V)
!
   CHARACTER*8      XID,       HMOLID,      YID
   Real*8               SECANT,       XALTZ
!
!
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
   COMMON /DEAMT/ DENM(MXMOL,MXZMD),DENP(MXMOL,MXPDIM),DRYAIR(MXZMD)
   COMMON /PATHD/ PBAR(MXLAY),TBAR(MXLAY),AMOUNT(MXMOL,MXLAY),       &
   &               WN2L(MXLAY),DVL(MXLAY),WTOTL(MXLAY),ALBL(MXLAY),   &
   &               ADBL(MXLAY),AVBL(MXLAY),H2OSL(MXLAY),IPATH(MXLAY), &
   &               ITYL(MXLAY),SECNTA(MXLAY),HT1,HT2,ALTZ(0:MXLAY),   &
   &               PZ(0:MXLAY),TZ(0:MXLAY)
   COMMON /ZOUTP/ ZOUT(MXLAY),SOUT(MXLAY),RHOSUM(MXLAY),             &
   &               AMTTOT(MXMOL),AMTCUM(MXMOL),ISKIP(MXMOL)
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WN2   ,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NDUM,LAYER ,YI1,YID(10),LSTWDF
!
   character*4 ht1,ht2

   I2 = IPMAX-1
   IOUT = 1
   PZ(0) = PP(1)
   TZ(0) = TP(1)
!
!     If entry in TAPE5 for TBOUND < 0, use TZ(O) as boundary
!     temperature
!
   IF (TBOUND.LT.0.) TBOUND = TZ(0)
!
   DO 20 IP = 1, I2
      PBAR(IOUT) = PBAR(IOUT)+PPSUM(IP)
      TBAR(IOUT) = TBAR(IOUT)+TPSUM(IP)
      RHOSUM(IOUT) = RHOSUM(IOUT)+RHOPSM(IP)
      SOUT(IOUT) = SOUT(IOUT)+SP(IP)
      DO 10 K = 1, NMOL
         AMOUNT(K,IOUT) = AMOUNT(K,IOUT)+AMTP(K,IP)
10    CONTINUE
      IF (ZPTH(IP+1).EQ.ZOUT(IOUT+1)) THEN
         PZ(IOUT) = PP(IP+1)
         TZ(IOUT) = TP(IP+1)
         IOUT = IOUT+1
      ENDIF
20 END DO
   IF (IOUT.NE.IOUTMX) GO TO 110
!
!     CALCULATE THE DENSITY WEIGHTED PRESSURE AND TEMPERATURE AND
!     ZERO OUT LAYER AMOUNTS AFTER 99.9 PERCENT OF THE TOTAL
!
   iskip(7) = 0
!
   DO 30 K = 1, NMOL
      AMTCUM(K) = 0.0
      ISKIP(K) = 0
      IF (AMTTOT(K).EQ.0.0) ISKIP(K) = 1
30 END DO
   L2 = IOUTMX-1
   LMAX = L2
   DO 90 L = 1, L2
      PBAR(L) = PBAR(L)/RHOSUM(L)
      TBAR(L) = TBAR(L)/RHOSUM(L)
!
!     ADJUST RHOSUM FOR THE PATH LENGTH IN CM NOT KM
!
      RHOSUM(L) = RHOSUM(L)*1.0E+5
!
      SUMAMT = 0.
      DO 40 K = 1, NMOL
         SUMAMT = SUMAMT+AMOUNT(K,L)
40    CONTINUE
      WN2L(L) = RHOSUM(L)-SUMAMT
!
!     CALCULATE 'EFFECTIVE SECANT' SECNTA
!
      SECNTA(L) = SOUT(L)/(ZOUT(L+1)-ZOUT(L))
      IF (L.EQ.1) ALTZ(0) = ZOUT(1)
      ALTZ(L) = ZOUT(L+1)
!
!     SET  IPATH
!
      IF (LEN.EQ.1) GO TO 50
      IF (H1.LT.H2) IPATH(L) = 3
      IF (H1.GT.H2) IPATH(L) = 1
      GO TO 60
50    CONTINUE
      IF (ZOUT(L).LT.HMID) IPATH(L) = 2
      IF (ZOUT(L).GE.HMID.AND.H1.GT.H2) IPATH(L) = 1
      IF (ZOUT(L).GE.HMID.AND.H1.LT.H2) IPATH(L) = 3
60    CONTINUE
!
!     TEST FOR ZEROING OF AMOUNTS
!
      ISKPT = 0
      nmol_max = nmol
      IF (ISKIP(7).EQ.1) nmol_max = nmol - 1
      FAC = 1.0
      IF (IPATH(L).EQ.2) FAC = 2.0
!
      DO 80 K = 1, NMOL

         IF (n_zero.ne.2) go to 70

         IF (ISKIP(K).NE.1) THEN
            IF (K.EQ.7 .OR. (IEMIT.EQ.1.AND.IPATH(L).NE.3)) GO TO 70
            IF (((AMTTOT(K)-AMTCUM(K))/AMTTOT(K)).GT.0.001) GO TO 70
         ENDIF
!
!     ZERO OUT THIS AMOUNT
!
         ISKIP(K) = 1
         AMOUNT(K,L) = 0.0
         ISKPT = ISKPT+1
!
!     IF ALL BUT O2 ARE ZEROED, ELIMINATE ALL HIGHER LAYERS
!
         IF (ISKPT.GE.(NMOL_max)) GO TO 100
70       CONTINUE
         AMTCUM(K) = AMTCUM(K)+FAC*AMOUNT(K,L)
80    CONTINUE
      LMAX = L
90 END DO
100 CONTINUE
   IOUTMX = LMAX+1
!
   RETURN
!
110 WRITE (IPR,900) IOUT,IOUTMX
!
   STOP ' ERROR FPACK '
!
900 FORMAT ('0FROM FPACK-  ERROR, IOUT = ',I5,'  DOES NOT MATCH ',    &
   &        'IOUTMX = ',I5)
!
end subroutine FPACK
!
!     ----------------------------------------------------------------
!
!
SUBROUTINE FIXTYP(IEMIT,FRH2O,ALFCOR,OLDDV,L,CINP)
!
!     *****************************************************************
!     This subroutine calculates ITYL, the ITYPE (ratio of DV from
!     one layer to the next) for each layer for output to TAPE7, if
!     desired (IFXTYP = 1).
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   CHARACTER*3 CINP
!
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID

   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
!
   DATA I_2/2/
!
   DV = 0.
!
!     Correct for water self broadening
!
   H2OSLF = (1.-FRH2O+5.*FRH2O)
   ALBAR = ALZERO*ALFCOR*H2OSLF
!
   AVBAR = 0.5*(ALBAR+SQRT(ALBAR*ALBAR+4.*ADBAR*ADBAR))
!
   DV = AVBAR/SAMPLE
!
   TYPE = 0.
   ITYPE = 99
!
!     DV is assumed to be less than 1
!     Set DV to 3 significant figures
!
   IF (L.EQ.1) THEN
      ISCAL = LOG10(DV)-3.
      SCAL = 10.**ISCAL
      IDV = (DV/SCAL)+0.5
!
!        Set IDV to be even
!
      IF (MOD(IDV,I_2).GT.0) IDV = IDV+1
      DV = SCAL* REAL(IDV)
   ELSE
      TYPE = OLDDV/DV
      TYPMAX = 2.5
      IF (TYPE.GT.TYPMAX) THEN
         IPROB = 1
         ISTOP = 1
      ELSEIF (TYPE.GE.1.2) THEN
!
!           TYPE is between 1.2 and TYPMAX
!
         DV = OLDDV
         ITYPE = 1./(TYPE-1.)+0.5
         IF (ITYPE.EQ.3) ITYPE = 2
         DV = OLDDV* REAL(ITYPE)/ REAL(ITYPE+1)
      ELSEIF (TYPE.GE.0.8) THEN
!
!           TYPE is between 0.8 and 1.2 (set to 1.0)
!
         DV = OLDDV
         ITYPE = 0
      ELSE
!
!           TYPE is less than 0.8
!
         DV = OLDDV
         ITYPE = 0
         IF (IEMIT.NE.1) THEN
            ITYPE = TYPE/(1.-TYPE)+0.5
            DV = DV* REAL(ITYPE+1)/ REAL(ITYPE)
            ITYPE = -ITYPE
         ENDIF
      ENDIF
   ENDIF
!
   OLDDV = DV
!
   WRITE(CINP,900) ITYPE
!
   RETURN
!
900 FORMAT(I3)
!
end subroutine FIXTYP
!
!     ----------------------------------------------------------------
!
SUBROUTINE XAMNTS (XV1,XV2)
!
!     *****************************************************************
!     THIS SUBROUTINE GENERATES THE ABSORBER AMOUNTS FOR THE SELECTED
!     HEAVY MOLECULES FOR WHICH "CROSS-SECTION" SPECTRAL DATA IS
!     AVAILABLE.  THE USER SELECTS THE DESIRED MOLECULES BY
!     USING THE CHEMICAL FORMULA, E.G. "CF2CL2" OR BY AN ALIAS, IN
!     THIS CASE, "F12".  THEREAFTER, THE MOLECULES ARE IDENTIFIED BY
!     AN INDEX.  THE USER MAY EITHER SELECT A STANDARD PROFILE OR READ
!     IN ONE.
!
!                             A.E.R. INC.     (AUGUST 1990)
!    *****************************************************************
!
   USE phys_consts, ONLY: alosmt
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC, MX_XS
!      PARAMETER (MXFSC=600,MXLAY=MXFSC+3,MXZMD=6000,                    &
!     &     MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,mx_xs=38,MXTRAC=22)
!
   IMPLICIT REAL*8           (V)
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
!     IFIL CARRIES FILE INFORMATION
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
!     LAMCHN CARRIES HARDWARE SPECIFIC PARAMETERS
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
   COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,n_zero,NOP,H1F,H2F,      &
   &                ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH,XVBAR, &
   &                HMINF,PHIF,IERRF,HSPACE
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON /MLATM/ ALT(MXZMD),PMDL(MXZMD,6),TMDL(MXZMD,6),            &
   &               AMOL(MXZMD,8,6),ZST(MXZMD),PST(MXZMD),             &
   &               TST(MXZMD),AMOLS(MXZMD,MXMOL)
   COMMON /DEAMT/ DENM(MXMOL,MXZMD),DENP(MXMOL,MXPDIM),DRYAIR(MXZMD)
!
!     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES
!     FOR THE CROSS-SECTION MOLECULES.
!     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES
!
   CHARACTER*10 XSFILE,XSNAME,ALIAS
   COMMON /XSECTF/ XSFILE(6,5,MX_XS),XSNAME(MX_XS),ALIAS(4,MX_XS)
   COMMON /XSECTR/ V1FX(5,MX_XS),V2FX(5,MX_XS),DVFX(5,MX_XS),        &
   &                WXM(MX_XS),NTEMPF(5,MX_XS),NSPECR(MX_XS),         &
   &                IXFORM(5,MX_XS),XSMASS(MX_XS),XDOPLR(5,MX_XS),    &
   &                NUMXS,IXSBIN
!
!     AMOLX(L,I)=MIXING RATIO (PPMV) OF THE I'TH MOLECULE FOR THE L'TH
!     LEVEL, ALTX(L)= ALTITUDE OF THE L'TH LEVEL, LAYXMX LEVELS MAX
!
   COMMON /MLATMX/ LAYXMX,ALTX(MXZMD),AMOLX(MXZMD,MX_XS)
   COMMON /PATHD/ PBAR(MXLAY),TBAR(MXLAY),AMOUNT(MXMOL,MXLAY),       &
   &               WN2L(MXLAY),DVL(MXLAY),WTOTL(MXLAY),ALBL(MXLAY),   &
   &               ADBL(MXLAY),AVBL(MXLAY),H2OSL(MXLAY),IPATH(MXLAY), &
   &               ITYL(MXLAY),SECNTA(MXLAY),HT1,HT2,ALTZ(0:MXLAY),   &
   &               PZ(0:MXLAY),TZ(0:MXLAY)
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
!     IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE
!     MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES
!     (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR
!     L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD
!     MOLECULES.
!
   COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(MX_XS),XAMNT(MX_XS,MXLAY)
   COMMON /ZOUTP/ ZOUT(MXLAY),SOUT(MXLAY),RHOSUM(MXLAY),             &
   &               AMTTOT(MXMOL),AMTCUM(MXMOL),ISKIP(MXMOL)
   COMMON /PCHINF/ MUNITS,CTYPE(MXLAY)
!
   DIMENSION XAMNTT(MX_XS)
!
   CHARACTER*48 CFORM1,CFORM2
   CHARACTER*10 HOTHER
   CHARACTER*7 PAFORM(2)
   CHARACTER*4 PZFORM(5),ht1,ht2
   CHARACTER*3 CTYPE
!
   DATA HOTHER / ' OTHER    '/
   DATA PZFORM / 'F8.6','F8.5','F8.4','F8.3','F8.2'/
   DATA PAFORM / '1PE15.7','  G15.7'/
   DATA CFORM1 / '(1PE15.7,0PF10.2,10X,A3,I2,1X,2(F7.3,F8.3,F7.2))'/
   DATA CFORM2 / '(  G15.7,0PF10.2,10X,A3,I2,23X,(F7.3,F8.3,F7.2))'/
!
   WRITE (IPR,900)
!
   NOZSAV = n_zero
   n_zero = 1
   IOMXSV = IOUTMX
!
!     READ IN THE NUMBER OF MOLECULES IXMOLS, AND THE FLAG IPRFL
!     INDICATING WHETHER A STANDARD PROFILE (0) OR A USER-INPUT PROFILE
!     (1) WILL BE USED.
!
   READ (IRD,905) IXMOLS,IPRFL,IXSBIN
   WRITE (IPR,910) IXMOLS,IPRFL
!
   IF (IPRFL.EQ.0) THEN
      WRITE (IPR,915)
   ELSEIF (IPRFL.EQ.1) THEN
      WRITE (IPR,920)
   ELSE
      WRITE (IPR,925)
      STOP 'STOPPED IN XAMNTS'
   ENDIF
!
!     READ IN DESIRED 'CROSS SECTIONS'
!
   CALL XSREAD (XV1,XV2)
!
   WRITE (IPR,930) (I,XSNAME(I),I=1,IXMOLS)
!
!     CALL XPROFL TO GENERATE THE DENSITY PROFILES OF THE CROSS-SECTION
!     MOLECULES.  THE PROFILES OF THE X-MOLECULES WILL BE STORED IN
!     DENM(J,I), J=1,IXMOLS, I=1,IMMAX, AT THE LEVELS ZMDL(IMMAX)
!
   CALL XPROFL (IPRFL)
!
!     SET NMOL = IXMOLS FOR THE AMOUNT CALCULATION, BUT RESET IT AFTER.
!     ALSO SET NOPRNT TO 1 TO SUPRESS PRINTING IN THE RAYTRACE.
!
   NMOLSV = NMOL
   NMOL = IXMOLS
   NOPRSV = NOPRNT
   NOPRNT = 1
!
   DO 10 I = 1, IOUTDM
      DO 8 K = 1, IXMAX
         XAMNT(K,I) = 0.0
8     CONTINUE
10 END DO
!
!     GET THE STANDARD-FORM SLANT PATH PARAMETERS H1, H2, ANGLE, PHI,
!     HMIN, AND LEN FROM /ADRIVE/ AS H1F, H2F, ANGLEF, PHIF, HMINF
!     AND LENF. USE THESE AS THE INPUTS TO THE RAYTRACE SUBROUTINE
!     RFPATH TO CALCULATE THE ABSORBER AMOUNTS.
!
   IF (ITYPE.EQ.1) THEN
!
!     =>  HORIZONTAL PATH
!
!        > GET NUMBER DENSITIES OF X-MOLECULES AT H1F <
!
      IF (IMMAX.EQ.1) THEN
!
!           > ITYPE = 1, HOMOGENOUS PATH <
!
         PH = PM(1)
         TH = TM(1)
         DO 20 K = 1, IXMOLS
            DENP(K,1) = DENM(K,1)
20       CONTINUE
!
      ELSE
!
!           > INTERPOLATE NUMBER DENSITIES TO H1F <
!
         ZH = H1F
         DO 30 L = 1, IMMAX
            IF (ZH.LT.ZMDL(L)) GO TO 40
30       CONTINUE
         L = IMMAX
40       CONTINUE
         A = (ZH-ZMDL(L-1))/(ZMDL(L)-ZMDL(L-1))
         CALL EXPINT (PH,PM(L-1),PM(L),A)
         TH = TM(L-1)+(TM(L)-TM(L-1))*A
         DO 50 K = 1, IXMOLS
            CALL EXPINT (DENP(K,1),DENM(K,L-1),DENM(K,L),A)
50       CONTINUE
      ENDIF
!
!     > CALCULATE PATH AMOUNTS <
!
      DO 60 K = 1, IXMOLS
         XAMNT(K,1) = DENP(K,1)*RANGEF*1.0E+5
60    CONTINUE
      RANGE = RANGEF
!
      LMAX = NLAYRS
      IOUTMX = LMAX+1
!
   ELSE
!
!     => SLANT PATH
!
!        > ZERO OUT ARRAYS <
!
      DO 70 N = 1, IPDIM
         IF (N.LE.IPDIM-2) THEN
            ZPTH(N) = 0.0
            PP(N) = 0.0
            TP(N) = 0.0
            RFNDXP(N) = 0.0
            SP(N) = 0.0
            PPSUM(N) = 0.0
            TPSUM(N) = 0.0
            RHOPSM(N) = 0.0
         ENDIF
         DO 68 M = 1, KDIM
            DENP(M,N) = 0.0
            AMTP(M,N) = 0.0
68       CONTINUE
70    CONTINUE
!
!        > CALCULATE THE REFRACTIVITY <
!
      WRITE(IPR,*) '   - Using LOWTRAN6 refractive index -'
!
      DO 80 IM = 1, IMMAX
         PPH2O = AMOLS(IM,1)*PZERO*TM(IM)/(TZERO*ALOSMT)
!
!	    Approximation to refraction index (from LOWTRAN5)
!
!           RFNDXM(IM) = ((77.46+0.459E-8*XVBAR**2)*PM(IM)/TM(IM)-
!    *                   (PPH2O/1013.0)*(43.49-0.347E-8*XVBAR**2))*
!    *                   1.0E-6
!
!	    Approximation to refraction index (from LOWTRAN6)
!
         RFNDXM(IM)=((83.42+(185.08/(1.0-(XVBAR/1.14E+5)**2))+       &
            (4.11/(1.0-(XVBAR/6.24E+4)**2)))*(PM(IM)*288.15)/ (1013.25* &
            TM(IM))-(43.49-(XVBAR/1.7E+4)**2)*(PPH2O/1013.25)) *1.0E-06
80    CONTINUE
      CALL RFPATH (H1F,H2F,ANGLEF,PHIF,LENF,HMINF,1,RANGE,BETA,      &
         BENDNG)
!
!        > CROSS-SECTION ABSORBER AMOUNTS ARE NOW IN AMTP(J,I).   <
!        > CONDENSE THE AMOUNTS INTO THE LAYERS DEFINDED BY ZOUT. <
!
      I2 = IPMAX-1
      IOUT = 1
      DO 100 IP = 1, I2
!
         DO 90 K = 1, IXMOLS
            XAMNT(K,IOUT) = XAMNT(K,IOUT)+AMTP(K,IP)
90       CONTINUE
         IF (ZPTH(IP+1).EQ.ZOUT(IOUT+1)) IOUT = IOUT+1
!
100   CONTINUE
!
      IF (IOUT.NE.IOUTMX) THEN
         WRITE (IPR,935) IOUT,IOUTMX
         STOP 'STOPPED IN XAMNTS, IOUT .NE. IOUTMX'
      ENDIF
!
      IOUTMX = IOMXSV
      LMAX = IOUTMX-1
!
   ENDIF
!
!     CROSS-SECTION AMOUNTS ARE NOW IN XAMNT. PRINT THEM OUT.
!                      (in E15.7 format)
!
   IF (IPUNCH.ge.1) THEN
      IFRMX = 1
      WRITE (IPU,940) IXMOLS,IXSBIN
      WRITE (IPU,945) (XSNAME(K),K=1,7),HOTHER,(XSNAME(K),K=8,NMOL)
      IF (ITYPE.EQ.1) THEN
         WRITE (IPU,950) IFRMX,LMAX,NMOL,SECNT0,HMOD,RANGE,ZH
!
!           -------------------------------------
!           > Write molecular information in    <
!           >  - mixing ratio if MUNITS is 1    <
!           >  - column density if MUNITS is 0  <
!           -------------------------------------
!
         IF (MUNITS.EQ.1) THEN
            DRAIR = WN2L(1)
            DO 105 M = 2,NMOLSV
               DRAIR = DRAIR + AMOUNT(M,1)
105         CONTINUE
!
!              > If DRAIR is zero, then write out XAMNT only    <
!              > (since XAMNT zero => mixing ratio zero)        <
!
            IF (DRAIR.EQ.0) THEN
               WRITE (IPU,955) PH,TH,IPATH(1),ZH,ZH, (XAMNT(K,1),K=1,&
                  7),WN2L(1), (XAMNT(K,1),K=8,NMOL)
            ELSE
               WRITE (IPU,955) PH,TH,IPATH(1),ZH,ZH, (XAMNT(K,1)/    &
                  DRAIR,K=1,7),WN2L(1), (XAMNT(K,1)/DRAIR,K=8,NMOL)
            ENDIF
         ELSE
!
!              Test to make sure there are no fractional molecular
!              amounts written out (will cause PATH to assume
!              mixing ratio)
!
            DO 107 K=1,NMOL
               IF (XAMNT(K,1).LT.1.) THEN
                  WRITE(IPR,1000) K,1
                  XAMNT(K,1) = 0.0
               ENDIF
107         CONTINUE
!
            WRITE (IPU,955) PH,TH,IPATH(1),ZH,ZH, (XAMNT(K,1),K=1,7),&
               WN2L(1), (XAMNT(K,1),K=8,NMOL)
         ENDIF
      ELSE
         WRITE (IPU,960) IFRMX,LMAX,NMOL,SECNT0,(HMOD(I),I=1,2),     &
            H1F,H2F,ANGLE,LENF
      ENDIF
   ENDIF
!
   WRITE (IPR,965) (XSNAME(I),I=1,IXMOLS)
!
   DO 110 K = 1, IXMOLS
      XAMNTT(K) = 0.0
110 END DO
!
   DO 130 L = 1, NLAYRS
!
!        > Write atmosphere to TAPE6 in column density <
!
      IF (ITYPE.EQ.1) THEN
         WRITE (IPR,970) L,ZOUT(L),ZOUT(L),(XAMNT(K,L),K=1,IXMOLS)
      ELSE
         WRITE (IPR,970) L,ZOUT(L),ZOUT(L+1),(XAMNT(K,L),K=1,IXMOLS)
      ENDIF
      DO 120 K = 1, IXMOLS
         FAC = 1.0
         IF (IPATH(L).EQ.2) FAC = 2.0
         XAMNTT(K) = XAMNTT(K)+FAC*XAMNT(K,L)
120   CONTINUE
!
      IF (IPUNCH.ge.1.AND.ITYPE.NE.1) THEN
         LTST = L
         IF (L.EQ.1) LTST = 0
         PTST = LOG10(PZ(LTST))
         NPTST = PTST+2
         IF (PTST.LT.0.0) NPTST = 1
         CFORM1(38:41) = PZFORM(NPTST)
         CFORM2(38:41) = PZFORM(NPTST)
         NPTST = 1
         IF (PBAR(L).GE.0.1) NPTST = 2
         CFORM1(2:8) = PAFORM(NPTST)
         CFORM2(2:8) = PAFORM(NPTST)
         IF (L.EQ.1) THEN
            WRITE (IPU,CFORM1) PBAR(L),TBAR(L),CTYPE(L),IPATH(L),    &
               ALTZ(L-1),PZ(L-1),TZ(L-1),ALTZ(L), PZ(L),TZ(L)
         ELSE
            WRITE (IPU,CFORM2) PBAR(L),TBAR(L),CTYPE(L),IPATH(L),    &
               ALTZ(L),PZ(L),TZ(L)
         ENDIF
!
!           -------------------------------------
!           > Write molecular information in    <
!           >  - mixing ratio if MUNITS is 1    <
!           >  - column density if MUNITS is 0  <
!           -------------------------------------
!
         IF (MUNITS.EQ.1) THEN
            DRAIR = WN2L(L)
            DO 125 M = 2,NMOLSV
               DRAIR = DRAIR + AMOUNT(M,L)
125         CONTINUE
!
!              > If DRAIR is zero, then write out XAMNT only    <
!              > (since XAMNT zero => mixing ratio zero)        <
!
            IF (DRAIR.EQ.0) THEN
               WRITE (IPU,975) (XAMNT(K,L),K=1,7),WN2L(L)
               IF (NMOL.GT.7) WRITE (IPU,975) (XAMNT(K,L),K=8,NMOL)
            ELSE
               WRITE (IPU,975) (XAMNT(K,L)/DRAIR,K=1,7),WN2L(L)
               IF (NMOL.GT.7) WRITE (IPU,975) (XAMNT(K,L)/DRAIR,K=8, &
                  NMOL)
            ENDIF
         ELSE
!
!              Test to make sure there are no fractional molecular
!              amounts written out (will cause PATH to assume
!              mixing ratio)
!
            DO 127 K=1,NMOL
               IF (XAMNT(K,L).LT.1.) THEN
                  WRITE(IPR,1000) K,L
                  XAMNT(K,L) = 0.0
               ENDIF
127         CONTINUE
!
            WRITE (IPU,975) (XAMNT(K,L),K=1,7),WN2L(L)
            IF (NMOL.GT.7) WRITE (IPU,975) (XAMNT(K,L),K=8,NMOL)
         ENDIF
      ENDIF
!
130 END DO
!
!_______________________________________________________________________
!     write cross sections to ifil_xs = 20 for derivatives

   if (ipunch.eq.2) then

      open (20,file='AJ_xs_amnts',                                   &
      &                        status='unknown',form = 'unformatted')
      write (20) IXMAX,IXMOLS,                                       &
      &        ( IXINDX(mol),(XAMNT(mol,l),l=1,nlayrs),mol=1,ixmols )

      close (20)

   endif
!_______________________________________________________________________
!
!     > Write atmosphere to TAPE6 in mixing ratio <
!
   WRITE(IPR,973)
   DO 135 L = 1, NLAYRS
      DRAIR = WN2L(L)
      DO 133 M = 2,NMOLSV
         DRAIR = DRAIR + AMOUNT(M,L)
133   CONTINUE
!
!        > If DRAIR is zero, then write out XAMNT only    <
!        > (since XAMNT zero => mixing ratio zero)        <
!
      IF (DRAIR.EQ.0) THEN
         IF (ITYPE.EQ.1) THEN
            WRITE (IPR,970) L,ZOUT(L),ZOUT(L), (XAMNT(K,L),K=1,      &
               IXMOLS)
         ELSE
            WRITE (IPR,970) L,ZOUT(L),ZOUT(L+1), (XAMNT(K,L),K=1,    &
               IXMOLS)
         ENDIF
      ELSE
         IF (ITYPE.EQ.1) THEN
            WRITE (IPR,970) L,ZOUT(L),ZOUT(L), (XAMNT(K,L)/DRAIR,K=1,&
               IXMOLS)
         ELSE
            WRITE (IPR,970) L,ZOUT(L),ZOUT(L+1), (XAMNT(K,L)/DRAIR,K=&
               1,IXMOLS)
         ENDIF
      ENDIF
135 END DO
!
   WRITE (IPR,980) (XAMNTT(K),K=1,IXMOLS)
!
!     DONE
!
   n_zero = NOZSAV
   NMOL = NMOLSV
   NOPRNT = NOPRSV
!
   RETURN
!
900 FORMAT ('1***** XAMNTS: ABSORBER AMOUNTS FOR THE CROSS-',         &
   &        'SECTION MOLECULES *****')
905 FORMAT (3I5)
910 FORMAT (/,'     IXMOLS      ISTD',/,2I10)
915 FORMAT (/,'     USER INPUT PROFILE SELECTED')
920 FORMAT (/,'     STANDARD PROFILE SELECTED')
925 FORMAT (/,'  ERROR: IPRFL IS NOT 0 OR 1, STOP')
930 FORMAT (/,'  THE CROSS-SECTION MOLECULES SELECTED ARE: ',/,/,     &
   &        (5X,I5,3X,A))
935 FORMAT (/,'  XAMNTS: ERROR- IOUT = ',I5,                          &
   &        '  DOES NOT MATCH IOUTMX = ',I5)
940 FORMAT (2(I5,5X),' THE FOLLOWING CROSS-SECTIONS WERE SELECTED:')
945 FORMAT (8A10)
950 FORMAT (1X,I1,I3,I5,F10.6,3A8,' * ',F7.3,' KM PATH AT ',F7.3,     &
   &        ' KM ALT')
955 FORMAT (E15.7,F10.4,10X,I5,1X,F7.3,15X,F7.3,/,(1P8E15.7))
960 FORMAT (1X,I1,I3,I5,F10.6,2A8,' H1=',F8.3,' H2=',F8.3,' ANG=',    &
   &        F8.3,' LEN=',I2)
965 FORMAT (//,'  LAYER AMOUNTS FOR THE CROSS-SECTION MOLECULES',//,  &
   &        '           LAYER          AMOUNTS (MOLS/CM2)',/,         &
   &        '   L    FROM     TO ',/,'        (KM)    (KM)',4X,8A10,  &
   &        /,25X,8A10)
970 FORMAT (1X,I3,2F8.3,3X,1P8E10.3,/,23X,1P8E15.7)
973 FORMAT ('1',3X,'------------------------------------',/,          &
   &            3X,'  MOLECULAR MIXING RATIOS BY LAYER',/)
975 FORMAT (1P8E15.7)
980 FORMAT (//,1X,'TOTAL AMOUNT FOR PATH ',1P8E15.7)
1000 FORMAT ('*** WARNING: Zeroing molecule #',i2.2,' amount ',        &
   &        'in layer #',i3.3)
!
end subroutine XAMNTS
!
!     ----------------------------------------------------------------
!
SUBROUTINE XPROFL (IPRFL)
!
!     *****************************************************************
!     THIS SUBROUTINE GENERATES THE DENSITY PROFILES OF THE CROSS-
!     SECTION MOLECULES.  IT STORES THE PROFILES IN THE ARRAY DENM IN
!     /DEAMT/ AT THE ALTITUDES ZMDL, WHICH ARE THE SAME ALTITUDES THAT
!     THE PROFILES OF THE MOLECULAR AMOUNTS ARE DEFINED ON.  (NOTE: THE
!     ACTUAL ALTITUDES USED ARE FROM ZST WHICH IS A COPY OF ZMDL.)
!     IPRFL IS A FLAG INDICATING THAT THE STANDARD PROFILES (0) OR A
!     USER-INPUT PROFILE (1) IS TO BE USED.
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC, MX_XS
!      PARAMETER (MXFSC=600,MXLAY=MXFSC+3,MXZMD=6000,                    &
!     &     MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,mx_xs=38,MXTRAC=22)
!
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)

!     IFIL CARRIES FILE INFORMATION
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
!     LAMCHN CARRIES HARDWARE SPECIFIC PARAMETERS
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN

   COMMON /c_drive/ ref_lat,hobs,ibmax_b,immax_b,                    &
   &                 lvl_1_2,jchar_st(10,2),wm(mxzmd)
!
   character*1 jchar_st
!
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
!
!     IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE
!     MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES
!     (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR
!     L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD
!     MOLECULES.
!
   COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(MX_XS),XAMNT(MX_XS,MXLAY)
   COMMON /MLATM/ ALT(MXZMD),PMDL(MXZMD,6),TMDL(MXZMD,6),            &
   &               AMOL(MXZMD,8,6),ZST(MXZMD),PST(MXZMD),             &
   &               TST(MXZMD),AMOLS(MXZMD,MXMOL)
   COMMON /DEAMT/ DENM(MXMOL,MXZMD),DENP(MXMOL,MXPDIM),DRYAIR(MXZMD)
!
!     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES
!     FOR THE CROSS-SECTION MOLECULES.
!     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES
!
   CHARACTER*10 XSFILE,XSNAME,ALIAS
   COMMON /XSECTF/ XSFILE(6,5,MX_XS),XSNAME(MX_XS),ALIAS(4,MX_XS)
!
!     AMOLX(L,I)=MIXING RATIO (PPMV) OF THE I'TH MOLECULE FOR THE L'TH
!     LEVEL, ALTX(L)= ALTITUDE OF THE L'TH LEVEL, LAYXMX LEVELS MAX
!
   COMMON /MLATMX/ LAYXMX,ALTX(MXZMD),AMOLX(MXZMD,MX_XS)

   DIMENSION ZX(MXZMD),DTMP(MX_XS,MXZMD),DENX(MX_XS,MXZMD)
   DIMENSION PX(MXZMD)
   DIMENSION ZTMP(2),PTMP(2),TTMP(2),WVTMP(2)
   CHARACTER JCHAR(MX_XS,MXZMD)*1,XTITLE*50
!
!     LOAD THE PROFILES OF ALTITUDE, PRESSURE, AND TEMPERATURE THAT
!     WERE USED TO CALCULATE THE MOLECULAR AMOUNTS BACK INTO THE
!     ARRAYS ZMDL, PM, AND TM FROM THE ARRAYS ZST, PST, AND TST
!
   DO 10 L = 1, IMMAX
      ZMDL(L) = ZST(L)
      PM(L) = PST(L)
      TM(L) = TST(L)
10 END DO
!
   IF (IPRFL.GT.0) THEN
!
!     A STANDARD PROFILE FOR X-MOLECULES DENSITY PROFILES HAS BEEN
!     SELECTED. THE PROFILES OF VOLUME MIXING RATIO ARE IN AMOLX
!     STORED AT THE LEVELS ALTX. LOAD THE ALTITUDES INTO ZX AND
!     DENX RESPECTIVELY.
!
      LAYX = LAYXMX
      DO 30 L = 1, LAYX
         ZX(L) = ALTX(L)
         DO 20 K = 1, IXMOLS
            DENX(K,L) = AMOLX(L,IXINDX(K))
20       CONTINUE
30    CONTINUE
!
   ELSE
!
!     A USER-INPUT PROFILE HAS BEEN SELECTED. READ IN THE PROFILES
!     AND INTERPOLATE THEM TO THE LEVELS ZMDL.
!
!     READ IN THE PROFILES. NOTE THAT ZORP CAN BE EITHER ALTITUDE
!     OR PRESSURE, DEPENDING UPON THE VALUE OF IZORP: 0 FOR
!     ALTITUDE, 1 FOR PRESSURE.
!
      WRITE (IPR,900)
      READ (IRD,905) LAYX,IZORP,XTITLE
      WRITE (IPR,910) LAYX,IZORP,XTITLE
      IF (LAYX.GT.LAYXMX) THEN
         WRITE (IPR,915) LAYXMX
         STOP 'STOPPED IN XPROFL'
      ENDIF
!
      WRITE (IPR,920) (K,K=1,IXMOLS)
!
      IF (IZORP .EQ. 0) THEN
         DO 50 L = 1, LAYX
            READ (IRD,925) ZX(L),(JCHAR(I,L),I=1,IXMOLS)
            WRITE (IPR,930) ZX(L),(JCHAR(I,L),I=1,IXMOLS)
!
            READ (IRD,935) (DTMP(K,L),K=1,IXMOLS)
            WRITE (IPR,940) (DTMP(K,L),K=1,IXMOLS)
50       CONTINUE
      ELSE
         DO 60 L = 1, LAYX
            READ (IRD,925) PX(L),(JCHAR(I,L),I=1,IXMOLS)
            WRITE (IPR,930) PX(L),(JCHAR(I,L),I=1,IXMOLS)
!
            READ (IRD,935) (DTMP(K,L),K=1,IXMOLS)
            WRITE (IPR,940) (DTMP(K,L),K=1,IXMOLS)
60       CONTINUE
      ENDIF
!
      IF (IZORP .EQ. 1) THEN
!
! INTERPOLATE PX GRID ONTO ZX GRID.

! TO ENSURE THAT CALCULATED/INPUT ZMDL'S WILL MATCH CALCULATED USER-LEVE
! ALTITUDES, A COMBINATION OF INTERPOLATION AND HYDROSTATICS ARE USED.
! ZBND = A * F1(P) + (1 - A) * F2(P), WHERE
! F1(P) = INTERPOLATION IN LN(P), F2(P) = HYDROSTATIC CALCULATION
         ISTART = 2

         DO 160 IP=1,LAYX
            PTMP(1) = 0.0
            TTMP(1) = 0.0
            WVTMP(1) = 0.0
            ZTMP(1) = 0.0

            PTMP(2) = 0.0
            TTMP(2) = 0.0
            WVTMP(2) = 0.0
            ZTMP(2) = 0.0

            DO 161 LIP=ISTART,IMMAX
               IF (PX(IP) .GT. PM(LIP)) GO TO 162
161         CONTINUE
            LIP=IMMAX
162         CONTINUE

            IF (PX(IP) .EQ. PM(LIP-1)) THEN
               ZX(IP) = ZMDL(LIP-1)
            ELSE

               IF(PX(IP) .EQ. PM(LIP)) THEN
                  ZX(IP) = ZMDL(LIP)
               ELSE

!     PERFORM INTERPOLATION IN LN(PM)
                  HIP = (ZMDL(LIP)-ZMDL(LIP-1))/ LOG(PM(LIP)/PM(LIP- &
                     1))
                  ZINT = ZMDL(LIP-1)+ HIP* LOG(PX(IP)/PM(LIP-1))

!     PERFORM ALTITUDE CALCULATION USING HYDROSTATIC EQUATION
                  PTMP(1) = PM(LIP-1)
                  ZTMP(1) = ZMDL(LIP-1)
                  TTMP(1) = TM(LIP-1)
                  WVTMP(1) = DENW(LIP-1)

                  PTMP(2) = PX(IP)

                  TIP = (TM(LIP)-TM(LIP-1))/ LOG(PM(LIP)/PM(LIP-1))
                  TTMP(2) = TM(LIP-1)+ TIP* LOG(PX(IP)/PM(LIP-1))

                  WVIP = (DENW(LIP)-DENW(LIP-1))/ LOG(PM(LIP)/PM(LIP-&
                     1))
                  WVTMP(2) = DENW(LIP-1) + WVIP* LOG(PX(IP)/PM(LIP-1)&
                     )

                  CALL CMPALT(2,PTMP,TTMP, WVTMP,ZTMP(1),REF_LAT,    &
                     ZTMP)
!     COMBINE THE INTERPOLATION AND THE HYDROSTATIC CALCULATION

                  RATP = LOG(PX(IP)/PM(LIP-1))/ LOG(PM(LIP)/PM(LIP-1))
                  if (immax_b.lt.0) then
                     A =0.
                  else
                     A = RATP**3
                  endif

                  ZX(IP) = A*ZINT + (1-A)*ZTMP(2)

               ENDIF
            ENDIF

            IF (IP .NE. 1) THEN
               IF (ZX(IP).LE.ZX(IP-1)) GO TO 300
            ENDIF

            ISTART = LIP
            CALL XTRACT (IP,DTMP,JCHAR,ZX(IP))
            DO 40 K = 1, IXMOLS
               DENX(K,IP) = DTMP(K,IP)
40          CONTINUE

160      CONTINUE
!
      ELSE
!            stop 'do 171'
!
!	   izorp is zero to get to this
!
         DO 171 L=1,LAYX
            CALL XTRACT (L,DTMP,JCHAR,ZX(L))
            DO 41 K = 1, IXMOLS
               DENX(K,L) = DTMP(K,L)
41          CONTINUE
171      CONTINUE
      ENDIF
!
   ENDIF
!
!     INTERPOLATE THE DENSITY PROFILE DENX DEFINED ON ZX TO DENM
!     DEFINED ON ZMDL, THEN CONVERT MIXING RATIO TO NUMBER DENSITY.
!
   CALL XINTRP (ZX,DENX,LAYX,IXMOLS)
!
   RETURN

! ERROR MESSAGES
300 WRITE(IPR,988) (ZX(I),I=1,LAYX)
   PRINT 988,(ZX(I),I=1,IP)

   STOP 'ZX IN XPROFL'
!
900 FORMAT (/,' READING IN A PROFILE FOR THE CROSS-SECTION',          &
   &        ' MOLECULES')
905 FORMAT (2I5,A)
910 FORMAT (/,'  LAYERS = ',I5,/,'     IZORP = ',I5,                  &
   &        '  (0 FOR ALTITUDE, 1 FOR PRESSURE)',/,                   &
   &        '  TITLE  = ',A50)
915 FORMAT (/,'  XPROFL: ERROR- LAYX > LAYXMX = ',I4)
920 FORMAT (/,'      Z OR P     JCHAR',/,'  ',                        &
   &        8('    DENX(',I2,')'))
925 FORMAT (F10.3,5X,38A1)
930 FORMAT (2X,F10.3,5X,38A1)
935 FORMAT (8E10.3)
940 FORMAT (2X,1p,8E12.3)
988 FORMAT (///,' ERROR: BOUNDARY ALTITUDES FOR CROSS_SECTION LEVELS',&
   &        'ARE NEGATIVE OR NOT IN ASCENDING ORDER',//,5X,' ZX ',    &
   &        /,(10F10.4))
!
end subroutine XPROFL
!
!     ----------------------------------------------------------------
!
SUBROUTINE XTRACT (ILEV,DTMP,JCHAR,Z)
!
!     *****************************************************************
!     FOR EACH MOLECULE K FOR WHICH JCHAR(K,ILEV) IS '1', THIS SUBROUTIN
!     INTERPOLATES THE MIXING RATIO DTMP(K,ILEV) AT THE ALTITUDE Z
!     FROM THE STANDARD PROFILE IN AMOLX ON THE ALTITUDE GRID ALTX.
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC, MX_XS
!      PARAMETER (MXFSC=600,MXLAY=MXFSC+3,MXZMD=6000,                    &
!     &     MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,mx_xs=38,MXTRAC=22)
!
!     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES
!     FOR THE CROSS-SECTION MOLECULES.
!     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES
!
   CHARACTER*10 XSFILE,XSNAME,ALIAS
   COMMON /XSECTF/ XSFILE(6,5,MX_XS),XSNAME(MX_XS),ALIAS(4,MX_XS)
!
!     AMOLX(L,I)=MIXING RATIO (PPMV) OF THE I'TH MOLECULE FOR THE L'TH
!     LEVEL, ALTX(L)= ALTITUDE OF THE L'TH LEVEL, LAYXMX LEVELS MAX
!
   COMMON /MLATMX/ LAYXMX,ALTX(MXZMD),AMOLX(MXZMD,MX_XS)
!
!     IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE
!     MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES
!     (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR
!     L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD
!     MOLECULES.
!
   COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(MX_XS),XAMNT(MX_XS,MXLAY)
!
   DIMENSION DTMP(MX_XS,MXZMD)
   CHARACTER*1 JCHAR(MX_XS,MXZMD)
!
!     FIND SMALLEST ALTX(L) GT Z
!
   DO 10 L = 2, LAYXMX
      IF (Z.LE.ALTX(L)) GO TO 20
10 END DO
   L = LAYXMX
!
20 CONTINUE
!
   DO 30 K = 1, IXMOLS
      IF (JCHAR(K,ILEV).EQ.'1') THEN
!
!     INTERPOLATE MIXING RATIO FROM STANDARD PROFILE
!
         A = (Z-ALTX(L-1))/(ALTX(L)-ALTX(L-1))
         CALL EXPINT (DTMP(K,ILEV),AMOLX(L,IXINDX(K)), AMOLX(L-1,    &
            IXINDX(K)),A)
      ENDIF
30 END DO
!
   RETURN
!
end subroutine XTRACT
!
!     ----------------------------------------------------------------
!
SUBROUTINE XINTRP (ZX,DENX,LAYX,IXMOLS)
!
!     *****************************************************************
!     THIS SUBROUTINE INTERPLOLATES THE PROFILE DENX ON THE ALTITUDE
!     GRID ZX INTO DENM ON THE GRID ZMDL.  EXPONENTIAL INTERPOLATION
!     IS USED.
!     *****************************************************************
!
   USE phys_consts, ONLY: alosmt
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC, MX_XS
!      PARAMETER (MXFSC=600,MXLAY=MXFSC+3,MXZMD=6000,                    &
!     &     MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,mx_xs=38,MXTRAC=22)
!
!     IFIL CARRIES FILE INFORMATION
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /CNSTATM/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL)
!
!     LAMCHN CARRIES HARDWARE SPECIFIC PARAMETERS
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
   COMMON /MLATM/ ALT(MXZMD),PMDL(MXZMD,6),TMDL(MXZMD,6),            &
   &               AMOL(MXZMD,8,6),ZST(MXZMD),PST(MXZMD),             &
   &               TST(MXZMD),AMOLS(MXZMD,MXMOL)
   COMMON /DEAMT/ DENM(MXMOL,MXZMD),DENP(MXMOL,MXPDIM),DRYAIR(MXZMD)
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
!
   DIMENSION ZX(MXZMD),DENX(MX_XS,MXZMD)
!
   LX = 2
   DO 30 L = 1, IMMAX
!
!        > FIND THE SMALLEST ZX GE ZMDL(L) <
!
10    CONTINUE
      IF (ZMDL(L).LE.ZX(LX).OR.LX.EQ.LAYX) THEN
         A = (ZMDL(L)-ZX(LX-1))/(ZX(LX)-ZX(LX-1))
         IF (A.LT.0.0 .OR. A.GT.1.0) WRITE (IPR,900)
!
!           > IF DRYAIR FOR LAYER NOT CALCULATED PREVIOUSLY <
!           > (USING NORMAL MOLECULES), THEN CALCULATE THE  <
!           > NUMBER DENSITY OF AIR                         <
!
         IF (DRYAIR(L).EQ.0.) DRYAIR(L) = ALOSMT*(PM(L)/PZERO)/(TM(L)&
            /TZERO)
!
         DO 20 K = 1, IXMOLS
            CALL EXPINT (DENM(K,L),DENX(K,LX-1),DENX(K,LX),A)
!
!              > CONVERT MIXING RATIO (PPMV) TO NUMBER DENSITY <
!
            DENM(K,L) = DRYAIR(L)*DENM(K,L)*1.0E-6
20       CONTINUE
         GO TO 30
      ELSE
         LX = LX+1
      ENDIF
      GO TO 10
!
30 END DO
!
   RETURN
!
900 FORMAT (//,'  XINTPL: CAUTION- EXTRAPOLATING X-SECTION PROFILE')
!
end subroutine XINTRP
!
! -------------------------------------------------------------------
!
BLOCK DATA XMLATM
!
!     *****************************************************************
!     THIS BLOCK DATA SUBROUTINE INITIALIZES THE STANDARD PROFILES
!     FOR THE "CROSS-SECTION" MOLECULES, THAT IS, THE MOLECULES FOR
!     WHICH THE SPECTRAL DATA IS IN THE FORM OF CROSS-SECTIONS
!     (ABSORPTION COEFFICIENTS) INSTEAD OF LINE PARAMETERS.
!     THE PROFILES OF VOLUME MIXING RATIOS GIVEN HERE ARE FROM:
!
!                ?????????????????????????????????
!
!     *****************************************************************
!
!     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES
!     FOR THE CROSS-SECTION MOLECULES.
!
!     AMOLX(L,I)=MIXING RATIO (PPMV) OF THE I'TH MOLECULE FOR THE L'TH
!     LEVEL, ALTX(L)= ALTITUDE OF THE L'TH LEVEL, LAYXMX LEVELS MAX
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC, MX_XS
!      PARAMETER (MXFSC=600,MXLAY=MXFSC+3,MXZMD=6000,                    &
!     &     MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,mx_xs=38,MXTRAC=22)
   PARAMETER (MXZ50=MXZMD-50)
!
   COMMON /MLATMX/ LAYXMX,ALTX(MXZMD),                               &
   &                AMOL1(MXZMD),  AMOL2(MXZMD),  AMOL3(MXZMD),       &
   &                AMOL4(MXZMD),  AMOL5(MXZMD),  AMOL6(MXZMD),       &
   &                AMOL7(MXZMD),  AMOL8(MXZMD),  AMOL9(MXZMD),       &
   &                AMOL10(MXZMD), AMOL11(MXZMD), AMOL12(MXZMD),      &
   &                AMOL13(MXZMD), AMOL14(MXZMD), AMOL15(MXZMD),      &
   &                AMOL16(MXZMD), AMOL17(MXZMD), AMOL18(MXZMD),      &
   &                AMOL19(MXZMD), AMOL20(MXZMD), AMOL21(MXZMD),      &
   &                AMOL22(MXZMD), AMOL23(MXZMD), AMOL24(MXZMD),      &
   &                AMOL25(MXZMD), AMOL26(MXZMD), AMOL27(MXZMD),      &
   &                AMOL28(MXZMD), AMOL29(MXZMD), AMOL30(MXZMD),      &
   &                AMOL31(MXZMD), AMOL32(MXZMD), AMOL33(MXZMD),      &
   &                AMOL34(MXZMD), AMOL35(MXZMD), AMOL36(MXZMD),      &
   &                AMOL37(MXZMD), AMOL38(MXZMD), AMOL39(MXZMD),      &
   &                AMOL40(MXZMD) 
!
   DATA LAYXMX / 6000 /
!
   DATA ALTX /                                                       &
   &       0.0,       1.0,       2.0,       3.0,       4.0,           &
   &       5.0,       6.0,       7.0,       8.0,       9.0,           &
   &      10.0,      11.0,      12.0,      13.0,      14.0,           &
   &      15.0,      16.0,      17.0,      18.0,      19.0,           &
   &      20.0,      21.0,      22.0,      23.0,      24.0,           &
   &      25.0,      27.5,      30.0,      32.5,      35.0,           &
   &      37.5,      40.0,      42.5,      45.0,      47.5,           &
   &      50.0,      55.0,      60.0,      65.0,      70.0,           &
   &      75.0,      80.0,      85.0,      90.0,      95.0,           &
   &     100.0,     105.0,     110.0,     115.0,     120.0,           &
   &     MXZ50*0.0/
!
!     DATA AMOL1 / CLONO2 /
!
   DATA AMOL1 /                                                      &
   &  4.737E-06, 4.162E-06, 3.587E-06, 2.891E-06, 2.195E-06,          &
   &  1.717E-06, 1.238E-06, 9.775E-07, 7.170E-07, 6.231E-07,          &
   &  5.292E-07, 5.313E-07, 5.334E-07, 1.451E-06, 2.368E-06,          &
   &  1.037E-05, 1.837E-05, 5.571E-05, 9.304E-05, 1.748E-04,          &
   &  2.566E-04, 3.681E-04, 4.796E-04, 5.910E-04, 7.024E-04,          &
   &  7.724E-04, 8.587E-04, 7.428E-04, 4.585E-04, 2.005E-04,          &
   &  5.867E-05, 8.818E-06, 1.319E-06, 1.610E-07, 1.889E-08,          &
   &  1.855E-09, 7.032E-11, 2.870E-12, 2.174E-13, 3.025E-14,          &
   &  3.257E-15, 2.634E-17, 3.313E-20, 2.134E-23, 1.366E-25,          &
   &  4.128E-28, 3.433E-30, 0.       , 0.       , 0.       ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL2 / HNO4 /
!
   DATA AMOL2 /                                                      &
   &  8.851E-07, 2.031E-06, 3.177E-06, 7.444E-06, 1.171E-05,          &
   &  1.962E-05, 2.752E-05, 3.284E-05, 3.816E-05, 3.576E-05,          &
   &  3.336E-05, 2.928E-05, 2.519E-05, 2.814E-05, 3.109E-05,          &
   &  4.406E-05, 5.703E-05, 8.127E-05, 1.055E-04, 1.344E-04,          &
   &  1.632E-04, 1.946E-04, 2.260E-04, 2.541E-04, 2.822E-04,          &
   &  3.010E-04, 3.016E-04, 2.175E-04, 1.151E-04, 4.626E-05,          &
   &  1.448E-05, 3.333E-06, 8.857E-07, 2.244E-07, 5.992E-08,          &
   &  1.669E-08, 2.490E-09, 3.659E-10, 7.913E-11, 3.168E-11,          &
   &  7.075E-12, 2.676E-13, 3.296E-15, 1.497E-18, 2.227E-21,          &
   &  1.918E-24, 2.507E-26, 3.283E-28, 1.605E-29, 2.026E-30,          &
   &  MXZ50*0.0/
!
!     DATA AMOL3 / CHCL2F /
!
   DATA AMOL3 /                                                      &
   &  50*-99.                                              ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL4 / CCL4 /
!
   DATA AMOL4 /                                                      &
   &  1.300E-04, 1.300E-04, 1.299E-04, 1.299E-04, 1.298E-04,          &
   &  1.297E-04, 1.296E-04, 1.295E-04, 1.294E-04, 1.293E-04,          &
   &  1.292E-04, 1.289E-04, 1.285E-04, 1.266E-04, 1.247E-04,          &
   &  1.187E-04, 1.127E-04, 1.026E-04, 9.256E-05, 8.037E-05,          &
   &  6.817E-05, 5.611E-05, 4.405E-05, 3.395E-05, 2.385E-05,          &
   &  1.701E-05, 5.027E-06, 8.202E-07, 1.204E-07, 1.304E-08,          &
   &  1.050E-09, 4.864E-11, 5.081E-12, 5.372E-13, 5.548E-14,          &
   &  5.688E-15, 2.281E-16, 5.092E-18, 1.699E-19, 3.184E-21,          &
   &  9.600E-23, 1.638E-24, 4.605E-26, 6.985E-28, 1.743E-29,          &
   &  2.224E-31, 4.283E-33, 0.       , 0.       , 0.       ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL5 / CCL3F /
!
   DATA AMOL5 /                                                      &
   &  1.400E-04, 1.400E-04, 1.399E-04, 1.399E-04, 1.398E-04,          &
   &  1.397E-04, 1.396E-04, 1.396E-04, 1.395E-04, 1.394E-04,          &
   &  1.392E-04, 1.389E-04, 1.386E-04, 1.368E-04, 1.349E-04,          &
   &  1.292E-04, 1.234E-04, 1.138E-04, 1.041E-04, 9.216E-05,          &
   &  8.021E-05, 6.799E-05, 5.576E-05, 4.480E-05, 3.384E-05,          &
   &  2.550E-05, 9.634E-06, 2.441E-06, 5.553E-07, 1.024E-07,          &
   &  1.581E-08, 1.939E-09, 3.811E-10, 7.716E-11, 1.585E-11,          &
   &  3.658E-12, 4.173E-13, 3.465E-14, 3.353E-15, 2.383E-16,          &
   &  2.084E-17, 1.346E-18, 1.080E-19, 6.099E-21, 4.246E-22,          &
   &  1.923E-23, 1.110E-24, 5.158E-26, 3.393E-27, 3.738E-28,          &
   &  MXZ50*0.0/
!
!     DATA AMOL6 / CCL2F2 /
!
   DATA AMOL6 /                                                      &
   &  2.400E-04, 2.400E-04, 2.399E-04, 2.399E-04, 2.398E-04,          &
   &  2.398E-04, 2.397E-04, 2.396E-04, 2.395E-04, 2.394E-04,          &
   &  2.393E-04, 2.390E-04, 2.387E-04, 2.370E-04, 2.353E-04,          &
   &  2.300E-04, 2.247E-04, 2.157E-04, 2.066E-04, 1.952E-04,          &
   &  1.838E-04, 1.712E-04, 1.585E-04, 1.452E-04, 1.319E-04,          &
   &  1.183E-04, 8.552E-05, 5.683E-05, 3.498E-05, 2.013E-05,          &
   &  1.111E-05, 6.014E-06, 3.446E-06, 1.998E-06, 1.181E-06,          &
   &  7.687E-07, 3.876E-07, 1.818E-07, 8.265E-08, 3.432E-08,          &
   &  1.380E-08, 4.984E-09, 1.704E-09, 4.917E-10, 1.272E-10,          &
   &  2.351E-11, 3.640E-12, 4.251E-13, 4.981E-14, 8.792E-15,          &
   &  MXZ50*0.0/
!
!     DATA AMOL7 / C2CL2F4 /
!
   DATA AMOL7 /                                                      &
   &  1.200E-05, 1.200E-05, 1.200E-05, 1.200E-05, 1.199E-05,          &
   &  1.199E-05, 1.199E-05, 1.199E-05, 1.198E-05, 1.198E-05,          &
   &  1.198E-05, 1.197E-05, 1.196E-05, 1.191E-05, 1.185E-05,          &
   &  1.167E-05, 1.149E-05, 1.120E-05, 1.090E-05, 1.053E-05,          &
   &  1.015E-05, 9.731E-06, 9.311E-06, 8.865E-06, 8.419E-06,          &
   &  7.949E-06, 6.770E-06, 5.620E-06, 4.547E-06, 3.623E-06,          &
   &  2.884E-06, 2.315E-06, 1.906E-06, 1.600E-06, 1.375E-06,          &
   &  1.231E-06, 1.037E-06, 8.645E-07, 7.140E-07, 5.799E-07,          &
   &  4.610E-07, 3.530E-07, 2.524E-07, 1.588E-07, 7.585E-08,          &
   &  2.131E-08, 3.107E-09, 2.089E-10, 1.084E-11, 8.968E-13,          &
   &  MXZ50*0.0/
!
!     DATA AMOL8 / C2CL3F3 /
!
   DATA AMOL8 /                                                      &
   &  1.900E-05, 1.900E-05, 1.899E-05, 1.899E-05, 1.898E-05,          &
   &  1.898E-05, 1.897E-05, 1.896E-05, 1.895E-05, 1.894E-05,          &
   &  1.893E-05, 1.890E-05, 1.887E-05, 1.871E-05, 1.854E-05,          &
   &  1.803E-05, 1.751E-05, 1.664E-05, 1.576E-05, 1.466E-05,          &
   &  1.356E-05, 1.236E-05, 1.116E-05, 9.931E-06, 8.702E-06,          &
   &  7.515E-06, 4.787E-06, 2.678E-06, 1.366E-06, 6.370E-07,          &
   &  2.820E-07, 1.222E-07, 5.930E-08, 2.949E-08, 1.507E-08,          &
   &  8.617E-09, 3.550E-09, 1.304E-09, 4.610E-10, 1.427E-10,          &
   &  4.322E-11, 1.131E-11, 2.861E-12, 5.798E-13, 1.059E-13,          &
   &  1.182E-14, 9.372E-16, 3.545E-17, 1.223E-18, 6.979E-20,          &
   &  MXZ50*0.0/
!
!     DATA AMOL9 / N2O5 /
!
   DATA AMOL9 /                                                      &
   &  1.312E-10, 4.065E-10, 6.818E-10, 5.329E-08, 1.059E-07,          &
   &  1.177E-06, 2.248E-06, 2.435E-06, 2.622E-06, 2.526E-06,          &
   &  2.430E-06, 2.714E-06, 2.998E-06, 7.354E-06, 1.171E-05,          &
   &  3.638E-05, 6.105E-05, 1.157E-04, 1.703E-04, 2.471E-04,          &
   &  3.239E-04, 4.204E-04, 5.168E-04, 6.318E-04, 7.468E-04,          &
   &  8.576E-04, 9.888E-04, 7.845E-04, 4.140E-04, 1.556E-04,          &
   &  4.229E-05, 7.489E-06, 1.426E-06, 2.195E-07, 3.706E-08,          &
   &  6.586E-09, 4.858E-10, 7.919E-12, 1.913E-13, 2.626E-15,          &
   &  3.692E-17, 5.125E-19, 2.169E-21, 6.096E-25, 6.336E-28,          &
   &  9.855E-32, 0.       , 0.       , 0.       , 0.       ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL10 / HNO3 /
!
   DATA AMOL10 /                                                     &
   &  5.738E-05, 6.671E-05, 7.603E-05, 8.176E-05, 8.748E-05,          &
   &  9.153E-05, 9.558E-05, 9.914E-05, 1.027E-04, 1.111E-04,          &
   &  1.195E-04, 1.431E-04, 1.667E-04, 3.217E-04, 4.766E-04,          &
   &  9.273E-04, 1.378E-03, 2.070E-03, 2.762E-03, 3.514E-03,          &
   &  4.266E-03, 4.891E-03, 5.516E-03, 5.858E-03, 6.200E-03,          &
   &  6.170E-03, 5.684E-03, 4.611E-03, 3.245E-03, 1.978E-03,          &
   &  1.015E-03, 3.855E-04, 1.252E-04, 3.480E-05, 9.533E-06,          &
   &  2.792E-06, 5.898E-07, 1.885E-07, 4.912E-08, 1.021E-08,          &
   &  2.233E-09, 1.122E-09, 3.566E-11, 3.213E-14, 7.770E-17,          &
   &  9.752E-20, 1.129E-21, 2.151E-23, 1.720E-24, 2.813E-25,          &
   &  MXZ50*0.0/
!
!     DATA AMOL11 / CF4 /
!
   DATA AMOL11 /                                                     &
   &  7.000E-05, 7.000E-05, 7.000E-05, 7.000E-05, 7.000E-05,          &
   &  7.000E-05, 7.000E-05, 7.000E-05, 7.000E-05, 7.000E-05,          &
   &  7.000E-05, 7.000E-05, 7.000E-05, 7.000E-05, 7.000E-05,          &
   &  7.000E-05, 7.000E-05, 7.000E-05, 7.000E-05, 7.000E-05,          &
   &  7.000E-05, 7.000E-05, 7.000E-05, 7.000E-05, 7.000E-05,          &
   &  7.000E-05, 7.000E-05, 6.906E-05, 6.711E-05, 6.547E-05,          &
   &  6.500E-05, 6.500E-05, 6.500E-05, 6.500E-05, 6.275E-05,          &
   &  5.700E-05, 5.262E-05, 5.000E-05, 4.641E-05, 3.989E-05,          &
   &  3.492E-05, 3.081E-05, 3.023E-05, 3.000E-05, 3.000E-05,          &
   &  3.000E-05, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,          &
   &  MXZ50*0.0/
!
!     DATA AMOL12 / CHCLF2 /  ????
!
   DATA AMOL12 /                                                     &
   &  9.350E-05, 9.350E-05, 9.350E-05, 9.350E-05, 9.350E-05,          &
   &  9.350E-05, 9.350E-05, 9.350E-05, 9.350E-05, 9.350E-05,          &
   &  9.350E-05, 9.350E-05, 9.350E-05, 9.350E-05, 9.350E-05,          &
   &  9.350E-05, 9.350E-05, 9.341E-05, 9.331E-05, 9.321E-05,          &
   &  9.310E-05, 9.297E-05, 9.285E-05, 9.282E-05, 9.280E-05,          &
   &  9.265E-05, 9.138E-05, 8.837E-05, 8.245E-05, 7.061E-05,          &
   &  5.868E-05, 4.835E-05, 4.288E-05, 3.561E-05, 2.736E-05,          &
   &  2.040E-05, 1.509E-05, 9.775E-06, 1.700E-06, 0.000E+00,          &
   &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,          &
   &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,          &
   &  MXZ50*0.0/
!
!     DATA AMOL13 / CCLF3 /
!
   DATA AMOL13 /                                                     &
   &  5.000E-06, 5.000E-06, 5.000E-06, 5.000E-06, 5.000E-06,          &
   &  5.000E-06, 5.000E-06, 5.000E-06, 5.000E-06, 5.000E-06,          &
   &  5.000E-06, 5.000E-06, 5.000E-06, 4.580E-06, 4.200E-06,          &
   &  4.050E-06, 3.900E-06, 3.750E-06, 3.600E-06, 3.450E-06,          &
   &  3.300E-06, 3.200E-06, 3.100E-06, 3.050E-06, 3.000E-06,          &
   &  2.900E-06, 2.570E-06, 2.200E-06, 2.100E-06, 2.050E-06,          &
   &  1.920E-06, 1.800E-06, 1.670E-06, 1.550E-06, 1.420E-06,          &
   &  1.300E-06, 1.080E-06, 8.970E-07, 7.460E-07, 6.200E-07,          &
   &  5.150E-07, 1.000E-07, 1.000E-08, 1.000E-09, 1.000E-10,          &
   &  1.000E-12, 1.000E-14, 1.000E-15, 1.000E-15, 1.000E-15,          &
   &  MXZ50*0.0/
!
!     DATA AMOL14 / C2CLF5 /
!
   DATA AMOL14 /                                                     &
   &  4.000E-06, 4.000E-06, 4.000E-06, 4.000E-06, 4.000E-06,          &
   &  4.000E-06, 4.000E-06, 4.000E-06, 4.000E-06, 4.000E-06,          &
   &  4.000E-06, 4.000E-06, 4.000E-06, 3.740E-06, 3.500E-06,          &
   &  3.240E-06, 3.000E-06, 2.790E-06, 2.600E-06, 2.450E-06,          &
   &  2.300E-06, 2.090E-06, 1.900E-06, 1.740E-06, 1.600E-06,          &
   &  1.500E-06, 1.320E-06, 1.200E-06, 1.070E-06, 9.490E-07,          &
   &  8.240E-07, 7.000E-07, 5.920E-07, 5.140E-07, 4.600E-07,          &
   &  4.200E-07, 3.530E-07, 2.970E-07, 2.500E-07, 2.110E-07,          &
   &  1.770E-07, 1.000E-09, 1.000E-11, 1.000E-13, 1.000E-15,          &
   &  1.000E-15, 1.000E-15, 1.000E-15, 1.000E-15, 1.000E-15,          &
   &  MXZ50*0.0/
!
!     DATA AMOL15 / NO2 /
!
   DATA AMOL15 /                                                     &
   &  2.30E-05,  2.30E-05,  2.30E-05,  2.30E-05,  2.30E-05,           &
   &  2.30E-05,  2.30E-05,  2.30E-05,  2.30E-05,  2.32E-05,           &
   &  2.38E-05,  2.62E-05,  3.15E-05,  4.45E-05,  7.48E-05,           &
   &  1.71E-04,  3.19E-04,  5.19E-04,  7.71E-04,  1.06E-03,           &
   &  1.39E-03,  1.76E-03,  2.16E-03,  2.58E-03,  3.06E-03,           &
   &  3.74E-03,  4.81E-03,  6.16E-03,  7.21E-03,  7.28E-03,           &
   &  6.26E-03,  4.03E-03,  2.17E-03,  1.15E-03,  6.66E-04,           &
   &  4.43E-04,  3.39E-04,  2.85E-04,  2.53E-04,  2.31E-04,           &
   &  2.15E-04,  2.02E-04,  1.92E-04,  1.83E-04,  1.76E-04,           &
   &  1.70E-04,  1.64E-04,  1.59E-04,  1.55E-04,  1.51E-04,           &
   &  MXZ50*0.0 /
!
!     DATA AMOL16 / PAN /  !GEOS-Chem profile, malvarad@aer.com
!     01/01/2006-01/16/2006 12 S 172 W, 2x25, 47L GEOS5
   DATA AMOL16 /                                                     &
   &  6.594E-08, 1.467E-07, 8.605E-07, 2.253E-06, 3.506E-06,          &
   &  5.169E-06, 1.088E-05, 1.620E-05, 1.892E-05, 2.213E-05,          &
   &  2.758E-05, 2.672E-05, 2.464E-05, 3.670E-05, 4.520E-05,          &
   &  4.595E-05, 3.741E-05, 2.511E-05, 1.568E-05, 8.231E-06,          &
   &  2.787E-06, 8.456E-07, 1.467E-07, 2.545E-08, 1.671E-08,          &
   &  1.128E-08, 0.       , 0.       , 0.       , 0.       ,          &
   &  0.       , 0.       , 0.       , 0.       , 0.       ,          &
   &  0.       , 0.       , 0.       , 0.       , 0.       ,          &
   &  0.       , 0.       , 0.       , 0.       , 0.       ,          &
   &  0.       , 0.       , 0.       , 0.       , 0.       ,          &
   &  MXZ50*0.0/

!     DATA AMOL17 / ACET (CH3COCH3) / 2 wk average GEOS-Chem profile
!     01/01/2006-01/16/2006, 58 N 152 W, 2x25, 47L GEOS5
!     malvarad@aer.com
   DATA AMOL17 /                                                     &
   &  4.770E-04, 5.405E-04, 6.139E-04, 7.226E-04, 8.225E-04,          &
   &  9.284E-04, 9.847E-04, 1.012E-03, 1.033E-03, 9.969E-04,          &
   &  8.522E-04, 5.500E-04, 2.464E-04, 1.033E-04, 5.134E-05,          &
   &  2.252E-05, 1.196E-05, 6.204E-06, 3.677E-06, 2.185E-06,          &
   &  1.145E-06, 5.874E-07, 2.714E-07, 1.254E-07, 1.025E-07,          &
   &  8.470E-08, 0.       , 0.       , 0.       , 0.       ,          &
   &  0.       , 0.       , 0.       , 0.       , 0.       ,          &
   &  0.       , 0.       , 0.       , 0.       , 0.       ,          &
   &  0.       , 0.       , 0.       , 0.       , 0.       ,          &
   &  0.       , 0.       , 0.       , 0.       , 0.       ,          &
   &  MXZ50*0.0/


!     DATA AMOL18 / CH3CN / !ARCTAS-B background profile, malvarad@aer.com
   DATA AMOL18 /                                                     &
   &  1.000E-04, 1.000E-04, 1.000E-04, 1.000E-04, 1.000E-04,          &
   &  1.000E-04, 1.000E-04, 1.000E-04, 1.000E-04, 1.000E-04,          &
   &  8.000E-05, 6.000E-05, 4.000E-05, 2.000E-05, 0.000E-04,          &
   &  0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04,          &
   &  0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04,          &
   &  0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04,          &
   &  0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04,          &
   &  0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04,          &
   &  0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04,          &
   &  0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04, 0.000E-04,          &
   &  MXZ50*0.0/
!
!     DATA AMOL19 / CHF2CF3 /
!
   DATA AMOL19 /                                                     &
   &  2.000E-06, 1.983E-06, 1.968E-06, 1.951E-06, 1.934E-06,          &
   &  1.919E-06, 1.902E-06, 1.887E-06, 1.870E-06, 1.852E-06,          &
   &  1.830E-06, 1.823E-06, 1.815E-06, 1.790E-06, 1.765E-06,          &
   &  1.726E-06, 1.688E-06, 1.631E-06, 1.575E-06, 1.492E-06,          &
   &  1.413E-06, 1.333E-06, 1.258E-06, 1.191E-06, 1.127E-06,          &
   &  1.076E-06, 9.556E-07, 8.289E-07, 7.052E-07, 5.940E-07,          &
   &  4.940E-07, 4.133E-07, 3.562E-07, 3.101E-07, 2.902E-07,          &
   &  2.839E-07, 2.679E-07, 2.655E-07, 2.680E-07, 2.702E-07,          &
   &  2.725E-07, 2.725E-07, 2.725E-07, 2.702E-07, 2.655E-07,          &
   &  2.609E-07, 2.609E-07, 2.609E-07, 2.609E-07, 2.609E-07,          &
   &  MXZ50*0.0 /
!
!     DATA AMOL20 / CFH2CF3 /
!
   DATA AMOL20 /                                                     &
   &  9.659E-06, 9.578E-06, 9.504E-06, 9.423E-06, 9.342E-06,          &
   &  9.268E-06, 9.186E-06, 9.113E-06, 9.031E-06, 8.943E-06,          &
   &  8.839E-06, 8.803E-06, 8.766E-06, 8.643E-06, 8.522E-06,          &
   &  8.335E-06, 8.153E-06, 7.875E-06, 7.606E-06, 7.204E-06,          &
   &  6.824E-06, 6.439E-06, 6.077E-06, 5.750E-06, 5.442E-06,          &
   &  5.196E-06, 4.615E-06, 4.003E-06, 3.406E-06, 2.869E-06,          &
   &  2.386E-06, 1.996E-06, 1.720E-06, 1.497E-06, 1.402E-06,          &
   &  1.371E-06, 1.294E-06, 1.282E-06, 1.294E-06, 1.305E-06,          &
   &  1.316E-06, 1.316E-06, 1.316E-06, 1.305E-06, 1.282E-06,          &
   &  1.260E-06, 1.260E-06, 1.260E-06, 1.260E-06, 1.260E-06,          &
   &  MXZ50*0.0 /
!
!     DATA AMOL21 / CF3CH3 /
!
   DATA AMOL21 /                                                     &
   &  2.000E-06, 1.983E-06, 1.968E-06, 1.951E-06, 1.934E-06,          &
   &  1.919E-06, 1.902E-06, 1.887E-06, 1.870E-06, 1.852E-06,          &
   &  1.830E-06, 1.823E-06, 1.815E-06, 1.790E-06, 1.765E-06,          &
   &  1.726E-06, 1.688E-06, 1.631E-06, 1.575E-06, 1.492E-06,          &
   &  1.413E-06, 1.333E-06, 1.258E-06, 1.191E-06, 1.127E-06,          &
   &  1.076E-06, 9.556E-07, 8.289E-07, 7.052E-07, 5.940E-07,          &
   &  4.940E-07, 4.133E-07, 3.562E-07, 3.101E-07, 2.902E-07,          &
   &  2.839E-07, 2.679E-07, 2.655E-07, 2.680E-07, 2.702E-07,          &
   &  2.725E-07, 2.725E-07, 2.725E-07, 2.702E-07, 2.655E-07,          &
   &  2.609E-07, 2.609E-07, 2.609E-07, 2.609E-07, 2.609E-07,          &
   &  MXZ50*0.0/
!
!     DATA AMOL22 / CH3CHF2 /
!
   DATA AMOL22 /                                                     &
   &  2.000E-06, 1.983E-06, 1.968E-06, 1.951E-06, 1.934E-06,          &
   &  1.919E-06, 1.902E-06, 1.887E-06, 1.870E-06, 1.852E-06,          &
   &  1.830E-06, 1.823E-06, 1.815E-06, 1.790E-06, 1.765E-06,          &
   &  1.726E-06, 1.688E-06, 1.631E-06, 1.575E-06, 1.492E-06,          &
   &  1.413E-06, 1.333E-06, 1.258E-06, 1.191E-06, 1.127E-06,          &
   &  1.076E-06, 9.556E-07, 8.289E-07, 7.052E-07, 5.940E-07,          &
   &  4.940E-07, 4.133E-07, 3.562E-07, 3.101E-07, 2.902E-07,          &
   &  2.839E-07, 2.679E-07, 2.655E-07, 2.680E-07, 2.702E-07,          &
   &  2.725E-07, 2.725E-07, 2.725E-07, 2.702E-07, 2.655E-07,          &
   &  2.609E-07, 2.609E-07, 2.609E-07, 2.609E-07, 2.609E-07,          &
   &  MXZ50*0.0/
!
!     DATA AMOL23 / CH2F2 /
!
   DATA AMOL23 /                                                     &
   &    3.854e-06,  3.803e-06,  3.755e-06,  3.723e-06,  3.707e-06,          &
   &    3.697e-06,  3.689e-06,  3.681e-06,  3.670e-06,  3.655e-06,          &
   &    3.629e-06,  3.584e-06,  3.534e-06,  3.491e-06,  3.442e-06,          &
   &    3.365e-06,  3.237e-06,  3.041e-06,  2.789e-06,  2.526e-06,          &
   &    2.294e-06,  2.079e-06,  1.906e-06,  1.772e-06,  1.670e-06,          &
   &    1.590e-06,  1.433e-06,  1.273e-06,  1.085e-06,  8.952e-07,          &
   &    7.474e-07,  6.506e-07,  5.913e-07,  5.419e-07,  4.862e-07,          &
   &    4.298e-07,  3.569e-07,  3.123e-07,  2.704e-07,  2.200e-07,          &
   &    1.640e-07,  1.233e-07,  1.014e-07,  8.941e-08,  8.403e-08,          &
   &    7.907e-08,  7.479e-08,  7.113e-08,  6.833e-08,  6.617e-08,          &
   &     MXZ50*0.0/
!     DATA AMOL24 / CCl2FCH3 /
!
!
   DATA AMOL24 /                                                     &
   &  1.120E-05, 1.111E-05, 1.102E-05, 1.093E-05, 1.083E-05,          &
   &  1.075E-05, 1.065E-05, 1.057E-05, 1.047E-05, 1.037E-05,          &
   &  1.025E-05, 1.021E-05, 1.016E-05, 1.002E-05, 9.882E-06,          &
   &  9.665E-06, 9.454E-06, 9.132E-06, 8.820E-06, 8.354E-06,          &
   &  7.913E-06, 7.466E-06, 7.046E-06, 6.667E-06, 6.310E-06,          &
   &  6.025E-06, 5.351E-06, 4.642E-06, 3.949E-06, 3.327E-06,          &
   &  2.767E-06, 2.315E-06, 1.995E-06, 1.737E-06, 1.626E-06,          &
   &  1.590E-06, 1.500E-06, 1.487E-06, 1.500E-06, 1.513E-06,          &
   &  1.526E-06, 1.526E-06, 1.526E-06, 1.513E-06, 1.487E-06,          &
   &  1.462E-06, 1.462E-06, 1.462E-06, 1.462E-06, 1.462E-06,          &
   &  MXZ50*0.0/
!
!     DATA AMOL25 / CH3CClF2 /
!
   DATA AMOL25 /                                                     &
   &  1.077E-05, 1.068E-05, 1.059E-05, 1.050E-05, 1.041E-05,          &
   &  1.033E-05, 1.024E-05, 1.016E-05, 1.007E-05, 9.967E-06,          &
   &  9.852E-06, 9.811E-06, 9.770E-06, 9.633E-06, 9.498E-06,          &
   &  9.290E-06, 9.087E-06, 8.777E-06, 8.478E-06, 8.030E-06,          &
   &  7.606E-06, 7.176E-06, 6.773E-06, 6.408E-06, 6.065E-06,          &
   &  5.791E-06, 5.144E-06, 4.462E-06, 3.796E-06, 3.198E-06,          &
   &  2.660E-06, 2.225E-06, 1.917E-06, 1.669E-06, 1.563E-06,          &
   &  1.528E-06, 1.442E-06, 1.429E-06, 1.442E-06, 1.454E-06,          &
   &  1.467E-06, 1.467E-06, 1.467E-06, 1.454E-06, 1.429E-06,          &
   &  1.405E-06, 1.405E-06, 1.405E-06, 1.405E-06, 1.405E-06,          &
   &  MXZ50*0.0/
!
!     DATA AMOL26 / ?????? /
!
   DATA AMOL26 /                                                     &
   &  50*-99.                                              ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL27 / CHCl2CF3 /
!
   DATA AMOL27 /                                                     &
   &  50*-99.                                              ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL28 / CHClFCF3 /
!
   DATA AMOL28 /                                                     &
   &  50*-99.                                              ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL29 / CHCl2C2F5 /
!
   DATA AMOL29 /                                                     &
   &  50*-99.                                              ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL30 / C2HCl2F3CF2 /
!
   DATA AMOL30 /                                                     &
   &  50*-99.                                              ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL31 / SO2 /
!
   DATA AMOL31 /                                                     &
   &  3.00E-04,  2.74E-04,  2.36E-04,  1.90E-04,  1.46E-04,           &
   &  1.18E-04,  9.71E-05,  8.30E-05,  7.21E-05,  6.56E-05,           &
   &  6.08E-05,  5.79E-05,  5.60E-05,  5.59E-05,  5.64E-05,           &
   &  5.75E-05,  5.75E-05,  5.37E-05,  4.78E-05,  3.97E-05,           &
   &  3.19E-05,  2.67E-05,  2.28E-05,  2.07E-05,  1.90E-05,           &
   &  1.75E-05,  1.54E-05,  1.34E-05,  1.21E-05,  1.16E-05,           &
   &  1.21E-05,  1.36E-05,  1.65E-05,  2.10E-05,  2.77E-05,           &
   &  3.56E-05,  4.59E-05,  5.15E-05,  5.11E-05,  4.32E-05,           &
   &  2.83E-05,  1.33E-05,  5.56E-06,  2.24E-06,  8.96E-07,           &
   &  3.58E-07,  1.43E-07,  5.73E-08,  2.29E-08,  9.17E-09,           &
   &  MXZ50*0.0 /
!
!     DATA AMOL32 / ISOP /
!
   DATA AMOL32 /                                                     &
   &  1.01E-04,  2.85E-05,  6.77E-06,  1.88E-06,  1.18E-06,           &
   &  1.15E-06,  1.04E-06,  1.09E-06,  1.27E-06,  9.39E-07,           &
   &  6.96E-07,  4.94E-07,  3.41E-07,  2.15E-07,  1.07E-07,           &
   &  3.13E-08,  4.64E-09,  2.72E-10,  1.26E-11,  1.08E-13,           &
   &  4.90E-14,  7.79E-17,  3.94E-17,  7.86E-19,  8.24E-21,           &
   &  2.21E-21,  3.53E-22,  3.50E-22,  3.46E-22,  3.42E-22,           &
   &  3.39E-22,  3.36E-22,  3.34E-22,  3.31E-22,  3.29E-22,           &
   &  3.27E-22,  3.22E-22,  2.45E-22,  1.89E-22,  1.65E-22,           &
   &  0.00E+00,  0.00E+00,  0.00E+00,  0.00E+00,  0.00E+00,           &
   &  0.00E+00,  0.00E+00,  0.00E+00,  0.00E+00,  0.00E+00,           &
   &  MXZ50*0.0  /
!
!
!     DATA AMOL33 / CHF3 /
!
   DATA AMOL33 /                                                     &
   &    2.144e-05,  2.140e-05,  2.136e-05,  2.134e-05,  2.132e-05,          &
   &    2.131e-05,  2.130e-05,  2.129e-05,  2.128e-05,  2.126e-05,          &
   &    2.122e-05,  2.114e-05,  2.106e-05,  2.098e-05,  2.088e-05,          &
   &    2.074e-05,  2.050e-05,  2.014e-05,  1.967e-05,  1.918e-05,          &
   &    1.874e-05,  1.834e-05,  1.801e-05,  1.777e-05,  1.759e-05,          &
   &    1.746e-05,  1.723e-05,  1.701e-05,  1.675e-05,  1.652e-05,          &
   &    1.635e-05,  1.624e-05,  1.619e-05,  1.613e-05,  1.604e-05,          &
   &    1.593e-05,  1.575e-05,  1.560e-05,  1.544e-05,  1.521e-05,          &
   &    1.489e-05,  1.456e-05,  1.434e-05,  1.422e-05,  1.417e-05,          &
   &    1.412e-05,  1.408e-05,  1.404e-05,  1.401e-05,  1.398e-05,          &
   &     MXZ50*0.0/
!
!     DATA AMOL34 / BrO /
!
   DATA AMOL34 /                                                     &
   &    0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  1.741e-11,          &
   &    4.149e-11,  6.558e-11,  8.966e-11,  1.281e-10,  1.853e-10,          &
   &    2.908e-10,  4.351e-10,  6.323e-10,  8.120e-10,  1.013e-09,          &
   &    1.211e-09,  1.416e-09,  1.657e-09,  1.898e-09,  2.174e-09,          &
   &    2.543e-09,  3.049e-09,  3.654e-09,  4.460e-09,  5.329e-09,          &
   &    6.482e-09,  1.106e-08,  1.855e-08,  2.969e-08,  4.489e-08,          &
   &    6.543e-08,  8.670e-08,  1.016e-07,  1.028e-07,  9.041e-08,          &
   &    7.027e-08,  2.802e-08,  0.000e+00,  0.000e+00,  0.000e+00,          &
   &    0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,          &
   &    0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,          &
   &  MXZ50*0.0/
!
!     DATA AMOL35 / ?????? /
!
   DATA AMOL35 /                                                     &
   &  50*-99.                                              ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL36 / ?????? /
!
   DATA AMOL36 /                                                     &
   &  50*-99.                                              ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL37 / ?????? /
!
   DATA AMOL37 /                                                     &
   &  50*-99.                                              ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL38 / ?????? /
!
   DATA AMOL38 /                                                     &
   &  50*-99.                                              ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL39 / ?????? /
!
   DATA AMOL39 /                                                     &
   &  50*-99.                                              ,          &
   &  MXZ50*0.0/
!
!     DATA AMOL39 / ?????? /
!
   DATA AMOL40 /                                                     &
   &  50*-99.                                              ,          &
   &  MXZ50*0.0/
   !
end block data XMLATM
!
! -------------------------------------------------------------------
!
SUBROUTINE NEWH2(H1,H2,ANGLE,RANGE,BETA,LEN,HTAN,PHI)
!
!     Changed for LBLRTM to correct geometry problems
!
!     THIS ROUTINE DETERMINES H2,BETA, TANGENT HEIGHT AND LEN.
!     ADOPTED FROM THE MODTRAN2 GEOMETRY PACKAGE
!
!     INPUTS ARE: H1, ZENTIH ANGLE (ANGLE) AND RANGE.
!     LEN = 1 IF THE PATH GOES THROUGH HTAN.
!
!     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO FASE01
!     MXLAY IS THE MAXIMUM NUMBER OF OUTPUT LAYERS
!     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
!         STORED IN ZMDL (INPUT)
!     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
!         OBTAINED BY MERGING ZMDL AND ZOUT
!     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600,MXLAY=MXFSC+3,MXZMD=6000,                    &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
!
!     BLANK COMMON FOR ZMDL
!
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
   REAL*8           CPATH,CPJ,CPJ1,SH,GAMMA,ANDEXD,CRFRCT,RE2
!
   CRFRCT(H)=(RE2+H)*ANDEXD(H,SH,GAMMA)
!
   RE2=RE
!     COMPUTE CPATH OR PATH CONSTANT
   CALL FNDSHD(H1,SH,GAMMA)
   CPATH = CRFRCT(H1)*SIN(ANGLE/DEG)
!
!     ANGLE = 90 at H1 implies that H1 = tangent height
!
   IF (ANGLE.EQ.90.0) THEN
      HTAN=H1
   ELSE
      DO 100 J=1,IMMAX
         IF (H1.GE.ZMDL(J)) JMAX=J
100   CONTINUE
      JMAX=JMAX+1
      ZJ1=ZMDL(JMAX)
      CPJ1=CRFRCT(ZJ1)
      HTAN=-1.0
      DO 200 J=JMAX,1,-1
         IF (HTAN.LT.0.0) THEN
            IF (J.EQ.1) THEN
               HTAN=0.0
            ELSE
               CPJ=CPJ1
               ZJ=ZJ1
               ZJ1=ZMDL(J-1)
               CPJ1=CRFRCT(ZJ1)
               IF ((CPATH.LE.CPJ).AND.(CPATH.GE.CPJ1)) THEN
                  HTAN=RTBIS(ZJ1,CPJ1,ZJ,CPJ,CPATH)
               ENDIF
            ENDIF
         ENDIF
200   CONTINUE
   ENDIF
!
!     Find H2, BETA AND LEN
!
   CALL FNDPTH(CPATH,H1,HTAN,H2,RANGE,BETA,LEN,ANGLE,PHI)
!
!     Ensure LEN is not reset in FSCGEO if direct path
   IF (LEN.EQ.0) HTAN=H2
!
!     IF (ANGLE .LE. 90.0) HTAN CARRIES HMIN NOT HTAN
   IF (ANGLE .LE. 90.0) HTAN = MIN(H1,H2)
!
   RETURN
!
end subroutine NEWH2
!
! ----------------------------------------------------------------
!
FUNCTION RTBIS(X1,CX1,X2,CX2,CPATH)
!
!     THIS FUNCTION FINDS THE ROOT OF
!            FUNC(X) = X*REFRACTIVE INDEX - CPA
!
!     THE ROOT IS ACTUALLY THE TANGENT HEIGHT, BETWEEN X1 AND X2.
!     THIS ROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS, ET AL.
!
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL

   REAL*8           CX1,CX2,CPATH,F,FMID,SH,GAMMA,ANDEXD
   DATA XACC/1E-5/
   PARAMETER (JMAX=40)
!
   FMID=CX2-CPATH
   F=CX1-CPATH
   IF(F*FMID.GE.0.) STOP 'ROOT MUST BE BRACKETED FOR BISECTION.'
   IF(F.LT.0.)THEN
      RTBIS=X1
      DX=X2-X1
   ELSE
      RTBIS=X2
      DX=X1-X2
   ENDIF
   DO 11 J=1,JMAX
      DX=DX*.5
      XMID=RTBIS+DX
      CALL FNDSHD(XMID,SH,GAMMA)
      FMID=ANDEXD(XMID,SH,GAMMA)*(XMID+RE)-CPATH
      IF(FMID.LE.0.)RTBIS=XMID
      IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
11 END DO
!
!     COMES HERE IF UNABLE TO SOLVE.
!
   IF (ABS(CX2) .LT. ABS(CX1)) THEN
      RTBIS = X2
   ELSE
      RTBIS = X1
   ENDIF
   RETURN
end function RTBIS
!
! ----------------------------------------------------------------
!
SUBROUTINE FNDPTH(CPATH,H1,HTAN,H2,RANGEI,BETA,LEN,ANGLE,PHI)
!
!     THIS ROUTINE DETERMINES H2, BETA AND LEN.
!     INPUTS ARE H1, HTAN (TANGENT HEIGHT), RANGE (RANGEI) AND
!     THE PATH CONSTANT, CPATH.
!     RANGEO IS THE OUTPUT RANGE WHICH SHOULD EQUAL THE INPUT RANGE.
!
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
!
   REAL*8           SAVE,STHETA,CAPRJ,PNTGRN,CTHETA,CTHET1,DX,       &
   &     DRNG,DBETA,R,DIFF,CPATH,ANDEXD,SH,GAMMA,RX,RATIO,RPLDR
!
   DATA DR/0.005/
!
   IF (RANGEI .LT. DR) STOP 'STOPPED IN FNDPTH'
!
!     (RANGEI .LT. DR) SHOULD NOT HAPPEN; SO THIS CHECK IS REDUNDANT.
!
   RANGEO = 0
   DO 200 I = 1, 2
!
      IF (ANGLE .LE. 90.0000 .AND. I .EQ. 1) GO TO 200
!
!        IF (ANGLE .LE. 90.0000) THE PATH DOES NOT GO THROUGH HTAN.
!        IF (ANGLE .LE. 90.0000) THE I = 1 CALCULATION SHOULD NOT BE DON
!        IF (ANGLE .LE. 90.0000) FOR I = 2, R1 = H1
!
      IF (I .EQ. 1) THEN
         R1 = H1
         R2 = HTAN
      ELSEIF (I .EQ. 2) THEN
         IF (HTAN .LT. 0.001 .AND. ANGLE .GT. 90) GO TO 200
!
!           IF (HTAN APPROXIMATELY 0) THEN YOU ARE ABOUT TO HIT THE EART
!
         R2 = ZMAX
         IF (ANGLE .LE. 90.0000) THEN
            R1 = H1
         ELSE
            R1 =HTAN
         ENDIF
      ENDIF
      IF (R2 .LT. R1) THEN
         DZ = -DR
      ELSE
         DZ = DR
      ENDIF
!

      z = r1
      DO 100 while (z.lt.r2)
         Z2=Z
         R=Z+RE
         CALL FNDSHD(Z2,SH,GAMMA)
         RX=ANDEXD(Z2,SH,GAMMA)
         STHETA = CPATH/(RX*R)
         IF (STHETA .GT. 1.0) STHETA = 1.
         IF (STHETA .LT.-1.0) STHETA =-1.
         SAVE = STHETA
         CTHETA = SQRT(1.0-STHETA**2)
         IF (R1 .GT. R2) CTHETA = -CTHETA
!
!           IF (R1 .GT. R2) THEN CTHETA IS NEGATIVE BECAUSE THETA .GT. 9
!
         RATIO=-(RX*SH)/(RX-1.0)
         CAPRJ = -R/RATIO
         PNTGRN = 1.0/(1.0-CAPRJ*STHETA*STHETA)
         RPLDR = R+DZ
         Z2 = Z+DZ
         CALL FNDSHD(Z2,SH,GAMMA)
         RX=ANDEXD(Z2,SH,GAMMA)
         STHETA = CPATH/(RX*RPLDR)
         CTHET1 = CTHETA
         CTHETA = SQRT(1.0-STHETA**2)
         IF (R1 .GT. R2) CTHETA = -CTHETA
         DX=CTHETA*DZ+(CTHETA-CTHET1)*R
         DRNG = PNTGRN*DX
         RANGEO = RANGEO + DRNG
!
         DBETA = (((SAVE+STHETA)*0.5) * (PNTGRN*DX)) / (Z-0.5*DZ+RE)
         BETA = BETA+DBETA
         IF (RANGEO .GE. RANGEI) THEN
            DIFF = (RANGEI-(RANGEO-DRNG))
            H2 = Z + (DZ/DRNG)*DIFF
            BETA = BETA*DEG
            IF (I .EQ. 2) THEN
               LEN = 1
               IF (ANGLE .LE. 90.0000) LEN = 0
               IF (H2 .LT. HTAN) THEN
!
!                    THIS WILL BE THE CASE IF I = 2, AND YOU HAVE
!                    GONE THROUGH THE R-LOOP BARELY (ONLY) ONCE.
!
                  H2 = HTAN
                  LEN = 0
               ENDIF
            ELSE
               LEN = 0
            ENDIF
!
!              CORRECTION FOR VERY SHORT PATHS; HERE IT IS ABOUT 5 KM
!
            IF (RANGEI .LT. 5.0 .AND. RANGEO/RANGEI .GT. 1.05) THEN
!
!                 CALCULATE BETA BY STARIGHT LINE GEOMETRY.
!
               PERP = SIN(ANGLE/DEG)*RANGEI
               BASE = COS(ANGLE/DEG)*RANGEI + RE+H1
               BETA = ATAN(PERP/BASE)*DEG
               RANGEO = RANGEI
!
!                 H2 = BASE - RE
!
               H2 = COS(ANGLE/DEG)*RANGEI+H1
            ENDIF
            PHI = 180.0 - ACOS(CTHETA)*DEG
            RETURN
         ENDIF
         z=z+dz
100   CONTINUE
200 END DO
!
!     COMES HERE IF YOU HAVE REACHED ZMAX, BUT YOUR RANGEI IS STILL
!     NOT EQUAL TO OUTPUT VALUE.
!     IN THIS CASE DO THE FOLLOWING.
!
   RANGEI = RANGEO
   H2 = ZMAX
   IF (ANGLE .LE. 90) THEN
      LEN = 0
   ELSE
      LEN = 1
   ENDIF
   IF (HTAN .LT. 0.001 .AND. ANGLE .GT. 90) THEN
!
!        YOU HAVE HIT THE EARTH IF YOU ARE AT THIS POINT OF THE CODE
!
      LEN = 0
      H2 = 0
   ENDIF
   BETA = BETA*DEG
   PHI = 180.0 - ACOS(CTHETA)*DEG
!
   RETURN
end subroutine FNDPTH
!
!     ----------------------------------------------------------------
!
SUBROUTINE FNDSHD (H,SH,GAMMA)
!
!     Double precision version of FINDSH - needed for improved geometry
!
!     *****************************************************************
!     GIVEN AN ALTITUDE H, THIS SUBROUTINE FINDS THE LAYER BOUNDARIES
!     Z(I1) AND Z(I2) WHICH CONTAIN H,  THEN CALCULATES THE SCALE
!     HEIGHT (SH) AND THE VALUE AT THE GROUND (GAMMA+1) FOR THE
!     REFRACTIVITY (INDEX OF REFRACTION -1)
!     *****************************************************************
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600,MXLAY=MXFSC+3,MXZMD=6000,                    &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
   COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),AVH(16),TX(16),W(16)
   COMMON WPATH(IM2,16),TBBY(IM2)
   COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),AVX2(47),AWCCON(5)
!
   REAL*8                SH,GAMMA
   CHARACTER*8      HMOD
!
   COMMON /CMN/HMOD(3),ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),RFNDXM(MXZMD),&
   &       ZPTH(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),  &
   &       TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),       &
   &       AMTP(MXMOL,MXPDIM)
!
   DO 10 IM = 2, IMMAX
      I2 = IM
      IF (ZMDL(IM).GE.H) GO TO 20
10 END DO
   I2 = IMMAX
20 CONTINUE
   I1 = I2-1
   CALL SCLHTD (ZMDL(I1),ZMDL(I2),RFNDXM(I1),RFNDXM(I2),SH,GAMMA)
!
   RETURN
!
end subroutine FNDSHD
!
!     ----------------------------------------------------------------
!
SUBROUTINE SCLHTD (Z1,Z2,RFNDX1,RFNDX2,SH,GAMMA)
!
!     Double precision version of SCALHT - needed for improved geometry
!
!     *****************************************************************
!     THIS SUBROUTINE CALCULATES THE SCALE HEIGHT SH OF THE (INDEX OF
!     REFRACTION-1.0) FROM THE VALUES OF THE INDEX AT THE ALTITUDES Z1
!     AND Z2 ( Z1 < Z2). IT ALSO CALCULATES THE EXTRAPOLATED VALUE
!     GAMMA OF THE (INDEX-1.0) AT Z = 0.0
!     *****************************************************************
!
   REAL*8           SH,GAMMA
!
   RF1 = RFNDX1+1.0E-20
   RF2 = RFNDX2+1.0E-20
   RATIO = RF1/RF2
   IF (ABS(RATIO-1.0).LT.1.0E-05) GO TO 10
   SH = (Z2-Z1)/ LOG(RATIO)
   GAMMA = RF1*(RF2/RF1)**(-Z1/(Z2-Z1))
   GO TO 20
10 CONTINUE
!
!     THE VARIATION IN THE INDEX OF REFRACTION WITH HEIGHT IS
!     INSIGNIFICANT OR ZERO
!
   SH = 0.0
   GAMMA = RFNDX1
20 CONTINUE
!
   RETURN
!
end subroutine SCLHTD
!
! ----------------------------------------------------------------
!
FUNCTION ANDEXD (H,SH,GAMMA)
!
!     Double precision version of ANDEX - needed for improved geometry
!
!     *****************************************************************
!     COMPUTES THE INDEX OF REFRACTION AT HEIGHT H, SH IS THE
!     SCALE HEIGHT, GAMMA IS THE VALUE AT H=0 OF THE REFRACTIVITY =
!     INDEX-1
!     *****************************************************************
!
   REAL*8     andexd, SH,GAMMA
!
   IF (SH.EQ.0.0) THEN
      ANDEXD = 1.0+GAMMA
   ELSE
      ANDEXD = 1.0+GAMMA*EXP(-H/SH)
   ENDIF
!
   RETURN
!
end function ANDEXD
!
! ----------------------------------------------------------------
!
FUNCTION RADRFD (H,SH,GAMMA)
!
!     Double precision version of RADREF - needed for improved geometry
!
!     *****************************************************************
!     COMPUTES THE RADIUS OF CURVATURE OF THE REFRACTED RAY FOR
!     A HORIZONTAL PATH.  RADREF = ANDEX/ D(ANDEX)/D(RADIUS)
!     *****************************************************************
!
   REAL*8     radrfd, SH,GAMMA,BIGNUM
   DATA BIGNUM / 1.0D36 /
!
   IF (SH.EQ.0.0) GO TO 10
   RADRFD = SH*(1.0+EXP(H/SH)/GAMMA)
!
   RETURN
!
10 RADRFD = BIGNUM
!
   RETURN
!
end function RADRFD
!

SUBROUTINE CMPALT(ILVL,PM,TM,DENW,REF_Z,REF_LAT,ZMDL)

!**************************************************************
!     AUTHOR: TONY CLOUGH, JENNIFER DELAMERE, JOHN WARDEN
!             JANUARY 2001
!     PROGRAM TO CALCULATE ALTITUDE LEVEL (ZMDL) GIVEN
!     PRESSURE (PM), TEMPERATURE (TM) AND THE NUMBER DENSITY
!     OF WATER VAPOR (DENW) USING THE HYDROSTATIC EQUATION
!
!     INPUT:
!      A) PRESSURE (MBAR)
!      B) TEMPERATURE (KELVIN)
!      C) NUMBER DENSITY OF WATER VAPOR
!
!     OUTPUT:
!      IDEAL GAS LAW: P.E.CIDDOR (1996), Refractive index of
!      air: New equations for the visible and near infrared,
!      Applied Optics, 35(9), 1566-1573.
!**************************************************************

   USE phys_consts, ONLY: boltz, gascon
   USE planet_consts, ONLY: xmass_dry, grav_const
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)

   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4

   COMMON /PARMTR/ DEG,GCAIR,RE,DELTAS,ZMIN,ZMAX,NOPRNT,IMMAX,       &
   &                IMDIM,IBMAX,IBDIM,IOUTMX,IOUTDM,IPMAX,            &
   &                IPHMID,IPDIM,KDIM,KMXNOM,NMOL
!
   REAL PM(MXZMD),TM(MXZMD),DENW(MXZMD),ZMDL(MXZMD)
   REAL H2O_MIXRAT(MXZMD),COMP_FACTOR(MXZMD),ZTEMP(MXZMD)

   REAL Y
   REAL CHI0
   REAL T0,DT
   REAL C1,C2,C3
   REAL A, B, ALPHA
!      REAL BTZ
   REAL XINT_TOT

   DATA CA0/1.58123E-6/,CA1/-2.9331E-8/,CA2/1.1043E-10/
   DATA CB0/5.707E-6/,CB1/-2.051E-8/
   DATA CC0/1.9898E-4/,CC1/-2.376E-6/
   DATA CD/1.83E-11/,CE/-0.0765E-8/

   DATA XMASS_H2O/0.018015/

! CALCULATE GRAVITY AT REFERENCE LATITUDE AT SURFACE

   G0 = GRAV_CONST(REF_LAT)

! CALCULATE THE NUMBER DENSITY OF TOTAL AIR MOLECULES [MOLEC/CM^3]
! CALCULATE THE COMPRESSIBILITY FACTOR (COMP_FAC) FOR THE
! IDEAL GAS LAW
   XMASS_RATIO = XMASS_H2O/XMASS_DRY
   DO 10 J=1,ILVL
      DT = TM(J) - 273.15
      TOTAL_AIR = PM(J)*1.0E+3/(BOLTZ*TM(J))
      DRY_AIR = TOTAL_AIR - DENW(J)
      H2O_MIXRAT(J) = DENW(J)/DRY_AIR
      CHIM = XMASS_RATIO*H2O_MIXRAT(J)
      COMP_FACTOR(J) = 1. - (PM(J)*100/TM(J))* (CA0 + CA1*DT + CA2*  &
         DT**2 + (CB0 + CB1*DT)*CHIM + (CC0 + CC1*DT)*CHIM**2) +        &
         (CD + CE*CHIM**2)*(PM(J)*100./TM(J))**2
10 END DO

! CONVERT REFERENCE ALTITUDE TO METERS

   ZTEMP(1) = REF_Z*1000.0
   ZMDL(1) = REF_Z

   DO 20 I=1, ILVL - 1
      GAVE = G0*(RE/(RE+ZTEMP(I)/1000.0))**2
      Y = LOG(PM(I+1)/PM(I))

      IF (Y .NE. 0.0) THEN
         CHI0 = H2O_MIXRAT(I)
         DCHI = (H2O_MIXRAT(I+1)-H2O_MIXRAT(I))/Y

         T0 = TM(I)
         DT = (TM(I+1) - TM(I))/Y

         C1 = T0 + T0*CHI0
         C2 = T0*DCHI + DT*CHI0 + DT
         C3 = DT*DCHI

         B = 1 + XMASS_RATIO*CHI0
         A = XMASS_RATIO*DCHI
         ALPHA = A/B

         IF ( ABS(ALPHA*Y) .GE. 0.01) THEN
            write(ipr,*) 'LAYER ',I, &
               ' THICKER THAN IDEAL FOR ALTITUDE CALCULATION'
         ENDIF

         XINT_TOT = C1*Y + 0.5*(C2-C1*ALPHA)*Y**2 + 0.3333*(C3-C2*   &
            ALPHA+C1*ALPHA**2)*Y**3

         XINT_TOT = -XINT_TOT*(GASCON*1.0E-7)/(XMASS_DRY*GAVE*B)

         ZTEMP(I+1) = ZTEMP(I) + XINT_TOT*COMP_FACTOR(I)
         ZMDL(I+1) = ZTEMP(I+1)/1000.
      ELSE
         ZTEMP(I+1) = ZMDL(I)*1000.0
         ZMDL(I+1) = ZMDL(I)
      ENDIF
20 CONTINUE
   RETURN

end subroutine CMPALT

