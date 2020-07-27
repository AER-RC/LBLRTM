!     path:      $HeadURL$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
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
SUBROUTINE PLTLBL (IENDPL)
!
   IMPLICIT REAL*8           (V)
!
   COMMON XX(2450),YY(2450)
!
   character*8      XID,       HMOL  ,      YID,PROGID
   real*8               SEC   ,       XALTZ
!
   COMMON /CVRPLT/ HNAMPLT,HVRPLT
   COMMON /PLTHDR/ XID(10),SEC,P0,T0,HMOL(60),XALTZ(4),              &
   &                W(60),PZL,PZU,TZL,TZU,WBROAD,DVT,V1V,V2V,TBOUND,  &
   &                EMISIV,FSCDID(17),NMOL,NLAYER,YID1,YID(10),LSTWDF
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
   COMMON /POINTS/ XXI,YI
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /FLFORM/ CFORM
   DIMENSION PROGID(3),PLTHDR(2),PNLHDR(2),DUM(2),JCFN(0:3)
!
   CHARACTER*8 FSTAT,NSTAT,OSTAT
   CHARACTER*8 COPT,CDIF,CRAT
   CHARACTER*18 HNAMPLT,HVRPLT

   CHARACTER*11 CFORM,BFORM,FFORM
   CHARACTER*25 TAPEJJ,TAPELL,TAPEMM,TAPE29,TAPEST(99),              &
   &             JFILEN,LFILEN,MFILEN,JOUTNM,CFILEN(3),BLNKNM
   CHARACTER CEX*2,CEXST*2,CTAPE*4,HOTHER*6,CPRGID*60
!
   LOGICAL OP,EX,IFIRST
!
   EQUIVALENCE (PLTHDR(1),XID(1)) , (PNLHDR(1),V1P),                 &
   &            (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),              &
   &            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),              &
   &            (FSCDID(5),IEMIT) , (FSCDID(6),ISCAN),                &
   &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),               &
   &            (FSCDID(9),JRAD) , (FSCDID(10),ITEST),                &
   &            (FSCDID(11),IMRG) , (FSCDID(12),SCNID),               &
   &            (FSCDID(13),HWHM) , (FSCDID(14),IDABS),               &
   &            (FSCDID(15),IATM) , (FSCDID(16),LAYR1),               &
   &            (FSCDID(17),NLAYFS)
   EQUIVALENCE (CFILEN(1),JFILEN) , (CFILEN(2),LFILEN),              &
   &            (CFILEN(3),MFILEN)
!
   DATA JCFN/1,1,3,3/
!
   DATA HOTHER / ' OTHER'/
   DATA CDIF / 'DIFRENCE'/,CRAT / ' RATIO  '/
   DATA CEXST/'EX'/,CTAPE / 'TAPE'/
   DATA BFORM/'FORMATTED  '/
   DATA NSTAT/'NEW'/,OSTAT/'OLD'/
   DATA TAPEST/99*'                         '/
   DATA BLNKNM/'                         '/
!
!                  PROGRAM TO PLOT LBLRTM RESULTS.
!
!     CARD 1: NAME, PHONE EXTENSION, ID
!       FORMAT (A30)
!
!     *********************************************
!     *** CARDS 2A AND 3A FOR PLOT              ***
!     *** CARD 2B FOR FILE DIFFERENCE OR RATIO  ***
!     *********************************************
!
!     CARD 2A: V1,V2,XSIZE,DELV,NUMSBX,NOENDX,LFILE,LSKIPF,SCALE,
!              IOPT,I4P,IXDEC
!       FORMAT (4F10.4,4I5,F10.3,I2,I3,I5)
!
!
!       V1 IS THE INITIAL WAVENUMBER OF THE PLOT,
!       V2 IS THE FINAL WAVENUMBER OF THE PLOT
!            (USE INTEGRAL VALUES OF CM-1)
!          A NEGATIVE NUMBER FOR V1 WILL TERMINATE PLOTTING SEQUENCE.
!
!       XSIZE IS THE NUMBER OF INCHES OF THE X-AXIS.
!
!       DELV IS THE NUMBER OF WAVENUMBERS (CM-1) PER MAJOR DIVISION.
!
!       NUMSBX IS THE NUMBER OF SUBDIVISIONS PER MAJOR DIVISION OF
!          X-AXIS, I.E., 1 GIVES NO SUBDIVISION TIC-MARKS,
!                        2 GIVES ONE SUBDIVISION MARK
!                          (IN OTHER WORDS 2 COMPARTMENTS
!                           PER MAJOR DIVISION),
!                        ETC.
!
!          ** FOR IOPT = 1, NUMSBX BECOMES THE QUANTITY JCAL WHICH  **
!          ** CONTROLS THE USE OF SYMBOLS ON THE OVERLAYED PLOT:    **
!          **                                                       **
!          **        JCAL = 0 WILL PRODUCE A LINE PLOT WITHOUT      **
!          **                 SYMBOLS                               **
!          **        JCAL > 0 WILL PRODUCE A LINE PLOT WITH SYMBOLS **
!          **                 AT EVERY JCALTH POINT                 **
!          **        JCAL < 0 WILL PRODUCE A POINT PLOT WITH A      **
!          **                 SYMBOL AT EVERY JCALTH POINT          **
!
!       NOENDX CONTROLS THE NUMBERS AT EITHER END OF THE X-AXIS, I.E.,
!            1 SUPRESSES THE NUMBERS AT EITHER END OF THE AXIS,
!            2 SUPRESSES THE BEGINNING NUMBER,
!            3 THE ENDING NUMBER, AND
!            0 LEAVES BOTH.
!
!          ** FOR IOPT=1, NOENDX BECOMES THE QUANTITY LCAL WHICH    **
!          ** IS A NUMBER DESCRIBING THE SYMBOL TO BE USED ON THE   **
!          ** OVERLAYED PLOT.  FOR A LIST OF SYMBOLS SEE THE AFGL   **
!          ** CYBER PLOTTER HANDOUT, APPENDIX B.                    **
!
!       LFILE IS THE TAPE NUMBER OF FILE TO BE READ FROM LBLRTM.
!
!       LSKIPF IS THE NUMBER OF FILES TO BE SKIPPED IN TAPE LFILE
!            (NUMBER OF FILE TO BE PLOTTED WILL BE LSKIPF + 1).
!
!       SCALE ENABLES ONE TO ENLARGE OR REDUCE A PLOT.
!
!       IOPT = 0 FOR PLOT,            = 1 FOR PLOT OVERLAY,
!            = 2 FOR FILE DIFFERENCE, = 3 FOR FILE RATIO.
!
!       I4P = 0 FOR LINEAR CONNECTION OF POINTS,
!           = 1 FOR FOUR POINT INTERP.
!
!       IXDEC IS NUMBER OF FIGURES AFTER DECIMAL POINT ON X-AXIS.
!
!
!     *** CARD 3A   (REQUIRED FOR CARD 2A, OTHERWISE OMIT) ***
!
!     CARD 3A: YMIN,YMAX,YSIZE,DELY,NUMSBY,NOENDY,IDEC,JEMIT,JPLOT,
!                                   LOGPLT,JHDR,JDUMMY,JOUT,JPLTFL
!         FORMAT (2G10.4,2F10.3,6I5,2(I2,I3))
!
!
!       YMIN IS ORIGIN VALUE OF Y-AXIS,
!       YMAX IS VALUE AT TOP OF Y-AXIS
!            (THESE VALUES SHOULD BE EXPONENTS WHEN LOG PLOT IS
!             DESIRED, WHERE THE DIFFERENCE YMAX-YMIN WILL COINCIDE
!             TO THE NUMBER OF DECADES SPANNED).
!
!          ** FOR IOPT = 1, YMIN DETERMINES THE Y-AXIS OFFSET OF   **
!          ** THE OVERLAYED PLOT WITH RESPECT TO THE ORIGINAL PLOT **
!
!       YSIZE IS THE NUMBER OF INCHES OF THE Y-AXIS.
!
!       DELY IS THE NUMBER OF UNITS OF Y VALUES PER MAJOR DIVISION
!            (= 1. WHEN LOG PLOT IS DESIRED).
!
!       NUMSBY IS THE NUMBER OF SUBDIVISIONS PER MAJOR DIVISION
!            OF Y-AXIS (SAME DEFINITION AS FOR NUMSBX).
!
!       NOENDY CONTROLS THE NUMBERS AT EITHER END OF Y-AXIS
!            (SAME DEFINITION AS FOR NOENDX).
!
!       IDEC IS NUMBER OF FIGURES AFTER DECIMAL POINT ON LINEAR Y-AXIS.
!
!       JEMIT = 0 FOR ABSORPTION, = 1 FOR RADIANCE.
!
!       JPLOT = 0 FOR TRANSMISSION, = 1 FOR OPTICAL DEPTH (JEMIT=0)
!       JPLOT = 0 FOR WATTS, = 1 FOR TEMPERATURE PLOT (JEMIT=1).
!       JPLOT = 2 ,JEMIT = 0 DECIBELS SCALE
!
!       LOGPLT = 0 FOR LINEAR Y-AXIS, = 1 FOR LOG SCALE.
!
!       JHDR = 0 FOR HEADER PANEL TO BE PLOTTED
!       JHDR = 1 FOR HEADER PANEL NOT PLOTTED
!
!       JDUMMY - NOT CURRENTLY USED
!
!       JOUT = 0 FOR PLOT TO SCREEN
!       JOUT = 1 FOR PLOT SAVED TO FILE
!       JOUT = 2 FOR PLOT TO SCREEN AND SAVED TO FILE
!       JOUT = 3 FOR PLOT PRINTED (ASCII) TO FILE
!       JOUT = 4 FOR PLOT TO SCREEN AND PRINTED (ASCII) TO FILE
!
!       JPLTFL IS THE FILE TO WHICH THE PLOT IS SAVED OR PRINTED
!          IF JOUT .GE. 1      (DEFAULT = TAPE29)
!
!
!      ** REPEAT CARD 2A & 3A, OR CARD 2B **
!
!
!     CARD 2B: V1,V2,JFILE,JSKIPF,LFILE,LSKIPF,IOPT,MFILE
!       FORMAT (2F10.4,20X,4I5,10X,I2,3X,I5)
!
!       SUBROUTINE FILOPT WILL DIFFERENCE OR RATIO
!         TWO LBLRTM OUTPUT FILES.
!
!              FILES MUST BE 1) UNFORMATTED LBLRTM FILES
!
!                            2) SINGLE QUANTITY FILES
!                                  I.E. CONTAIN ONLY ONE
!                                       OF THE FOLLOWING:
!
!                                        A) OPTICAL DEPTHS
!                                        B) TRANSMITTANCE
!                                        C) RADIANCE
!                                        D) TEMPERATURE
!
!                  *** NOTE: STANDARD LBLRTM OUTPUT FILES USUALLY
!                            CONTAIN BOTH TRANSMITTANCE AND RADIANCE.
!                            THE USER MUST EITHER SCAN OR PLOT THE
!                            DESIRED QUANTITY TO A FILE FOR USE AS
!                            INPUT TO THIS ROUTINE.
!
!       V1 IS THE INITIAL WAVENUMBER OF THE DIFFERENCE OR RATIO
!       V2 IS THE FINAL WAVENUMBER OF THE DIFFERENCE OR RATIO
!
!            A NEGATIVE NUMBER FOR V1 WILL TERMINATE PLOTTING SEQUENCE.
!
!       JFILE IS THE TAPE NUMBER OF FILE TO BE READ FROM LBLRTM.
!                    (NO DEFAULT)
!
!       JSKIPF IS THE NUMBER OF FILES TO BE SKIPPED IN TAPE JFILE
!            (NUMBER OF FILE TO BE USED WILL BE JSKIPF + 1).
!
!       LFILE IS THE TAPE NUMBER OF FILE TO BE READ FROM LBLRTM.
!
!       LSKIPF IS THE NUMBER OF FILES TO BE SKIPPED IN TAPE LFILE
!            (NUMBER OF FILE TO BE USED WILL BE LSKIPF + 1).
!
!       IOPT = 0 FOR PLOT,
!            = 1 FOR PLOT OVERLAY,
!            = 2 FOR FILE DIFFERENCE,
!            = 3 FOR FILE RATIO.
!
!          ***  DIFFERENCE WRITTEN ON MFILE IS (JFILE - LFILE)  ***
!          ***     RATIO WRITTEN ON MFILE IS (JFILE/LFILE)      ***
!
!       MFILE IS THE TAPE NUMBER OF FILE FOR DIFFERENCE/RATIO OUTPUT.
!
!       ** REPEAT CARDS 2A & 3A, OR CARD 2B **
!
!**********************************************************************
!
   DATA I_0/0/, I_1/1/, I_3/3/, I_6/6/, I_10/10/, I_100/100/
!
!     ASSIGN SCCS VERSION NUMBER TO MODULE
!
   HVRPLT =  '$Revision$'
!
   IENDPL = 0
   JHDR = 0
   JOUT = 0
   JPLTFL = 0
   JDUMMY = 0
   YMIN = 0.
   YMAX = 0.
   YSIZ = 0.
   DELY = 0.
   NUMSBY = 0.
   NOENDY = 0.
   IDEC = 0
   JEMIT = 0
   JPLOT = 0
   LOGPLT = 0
   NPLT = 0
   XLOGCN = - LOG10(EXP(1.))
   X3 = 0.
   YPL = 0.
   READ (IRD,900) CPRGID,CEX
   IEXTRC=0
   IF (CEX.EQ.CEXST) IEXTRC = 1
!
10 WRITE (IPR,902) CPRGID
   IFIRST=.TRUE.
   JFILEN = BLNKNM
   LFILEN = BLNKNM
   MFILEN = BLNKNM
!
   READ (IRD,905) V1,V2,XSIZ,DELV,NUMSBX,NOENDX,LFILE,LSKIPF,        &
   &               SCALE,IOPT,I4P,IXDEC
!
   IF (V1.LT.0) GO TO 240
!
   IF (IEXTRC.EQ.1) THEN
      READ (IRD,906) (CFILEN(J),J=1,JCFN(IOPT))
   ENDIF
!
   WRITE (TAPELL,930) CTAPE,LFILE
   IF (IOPT.LE.1) LFILEN = JFILEN
   IF (LFILEN.NE.BLNKNM) THEN
      CALL CLJUST(LFILEN,25)
      TAPELL = LFILEN
   ENDIF
   INQUIRE (UNIT=LFILE,OPENED=OP)
   IF (OP .AND. TAPEST(LFILE).NE.TAPELL) THEN
      CLOSE (LFILE)
      OP=.FALSE.
   ENDIF
   IF (.NOT.OP) THEN
      INQUIRE (FILE=TAPELL,EXIST=EX)
      FSTAT=NSTAT
      IF (EX) FSTAT=OSTAT
      OPEN (LFILE,FILE=TAPELL,STATUS=FSTAT,FORM=CFORM)
      TAPEST(LFILE) = TAPELL
   ENDIF
!
!     FOR IOPT = 0 - PLOT
!     FOR IOPT = 1 - OVERLAY PLOT
!     FOR IOPT = 2 - DIFFERENCE VALUES
!     FOR IOPT = 3 - RATIO VALUES
!
   IF (IOPT.GT.1) THEN
!sac     IF (V1.LT.0) GO TO 240
      JFILE = NUMSBX
      JSKIPF = NOENDX
      MFILE = IXDEC
      IRDOPT = 0
      CALL CLJUST(JFILEN,25)
      CALL CLJUST(LFILEN,25)
      CALL CLJUST(MFILEN,25)
      IF (JFILEN.EQ.BLNKNM) IRDOPT = IRDOPT+1
      IF (LFILEN.EQ.BLNKNM) IRDOPT = IRDOPT+1
      IF (IRDOPT.EQ.1) THEN
         WRITE(IPR,907) IRDOPT,JFILEN,LFILEN
         STOP ' IRDOPT '
      ELSE
         WRITE (TAPEJJ,930) CTAPE,JFILE
         IF (JFILEN.NE.BLNKNM) TAPEJJ = JFILEN
         INQUIRE (UNIT=JFILE,OPENED=OP)
         IF (OP .AND. TAPEST(JFILE).NE.TAPEJJ) THEN
            CLOSE (JFILE)
            OP=.FALSE.
         ENDIF
         IF (.NOT.OP) THEN
            INQUIRE (FILE=TAPEJJ,EXIST=EX)
            FSTAT=NSTAT
            IF (EX) FSTAT=OSTAT
            OPEN (JFILE,FILE=TAPEJJ,STATUS=FSTAT,FORM=CFORM)
            TAPEST(JFILE) = TAPELL
         ENDIF
         WRITE (TAPEMM,930) CTAPE,MFILE
         IF (MFILEN.NE.BLNKNM) TAPEMM = MFILEN
         INQUIRE (UNIT=MFILE,OPENED=OP)
         IF (OP .AND. TAPEST(MFILE).NE.TAPEMM) THEN
            CLOSE (MFILE)
            OP=.FALSE.
         ENDIF
         IF (.NOT.OP) THEN
            INQUIRE (FILE=TAPEMM,EXIST=EX)
            FSTAT=NSTAT
            IF (EX) FSTAT=OSTAT
            OPEN (MFILE,FILE=TAPEMM,STATUS=FSTAT,FORM=CFORM)
            TAPEST(MFILE) = TAPEMM
         ENDIF
      ENDIF
      CALL FILOPT (V1,V2,JFILE,JSKIPF,LFILE,LSKIPF,MFILE,LENGTH,     &
         IOPT)
      GO TO 10
   ENDIF
   IF ((V1.LT.0.).AND.(DELV.LT.0.)) THEN
      DELV = -DELV
   ELSE
      IF (V1.LT.0.) GO TO 240
   ENDIF
!
   JOUTNM = BLNKNM
   READ (IRD,910) YMINR,YMAX,YSIZ,DELY,NUMSBY,NOENDY,IDEC,JEMIT,     &
   &                  JPLOT,LOGPLT,JHDR,JDUMMY,JOUT,JPLTFL
   IF (IEXTRC.EQ.1) THEN
      READ (IRD,911) JOUTM
   ENDIF
   IF (JOUT.LT.0.OR.JOUT.GT.4) THEN
      WRITE (IPR,915) JOUT
      JOUT = 0
   ENDIF
   IFLAGJ = 1
   NOAXES = 0
   IF (IOPT.EQ.1) THEN
      IF (NPLT.EQ.0) THEN
         WRITE (IPR,920)
         GO TO 10
      ENDIF
      NOAXES = 1
      IFLAGJ = 2
      JHDR = 1
      V1 = V1STOR
      V2 = V2STOR
      XSIZ = XSIZES
      SCALE = SCALES
      YMIN = YMINST-YMINR
      YMAX = YMAXST-YMINR
      WRITE (IPR,927) YMINR,NUMSBX,NOENDX
   ELSE
      V1STOR = V1
      V2STOR = V2
      XSIZES = XSIZ
      SCALES = SCALE
      YMIN = YMINR
      YMINST = YMIN
      YMAXST = YMAX
   ENDIF
   IF (JOUT.EQ.1.OR.JOUT.EQ.3) THEN
      IFLAGJ = 0
   ELSE
      NPLT = NPLT+1
   ENDIF
   IF (NPLT.GT.1.AND.IFLAGJ.EQ.1) CALL PLOT (XSIZE+X3,-YPL,-3)
   IF (NPLT.GE.1.AND.IFLAGJ.GE.1) IENDPL = 1
   IF (SCALE.LT.0.01) SCALE = 1.
   IF (JOUT.GE.1) THEN
      IF (JPLTFL.EQ.0) JPLTFL = 29
      WRITE (TAPE29,930) CTAPE,JPLTFL
      IF (JOUTNM.NE.BLNKNM) THEN
         CALL CLJUST(JOUTNM,25)
         TAPE29 = JOUTNM
      ENDIF
!**      CALL CLJUST(JOUTNM,25)
!**      IF (JOUTNM.NE.BLNKNM) TAPE29 = JOUTNM
      INQUIRE (UNIT=JPLTFL,OPENED=OP)
      IF (OP .AND. TAPEST(JPLTFL).NE.TAPE29) THEN
         CLOSE (JPLTFL)
         OP=.FALSE.
      ENDIF
      IF (.NOT.OP) THEN
         INQUIRE (FILE=TAPE29,EXIST=EX)
         FSTAT=NSTAT
         IF (EX) FSTAT=OSTAT
         FFORM=CFORM
         IF (JOUT.GT.2) FFORM=BFORM
         OPEN (JPLTFL,FILE=TAPE29,STATUS=FSTAT,FORM=FFORM)
         TAPEST(JPLTFL) = TAPE29
      ENDIF
   ENDIF
   WRITE (IPR,932)
   WRITE (IPR,935) V1,V2,XSIZ,DELV,NUMSBX,NOENDX,LFILE,LSKIPF,SCALE, &
   &                IOPT,I4P,IXDEC
   WRITE (IPR,937)
   WRITE (IPR,940) YMIN,YMAX,YSIZ,DELY,NUMSBY,NOENDY,IDEC,JEMIT,     &
   &                JPLOT,LOGPLT,JHDR,JOUT,JPLTFL
   XSIZE = XSIZ*SCALE
   IF (JOUT.NE.1.AND.JOUT.NE.3.AND.NOAXES.EQ.0) THEN
      READ (CPRGID,945) PROGID
      CALL PLTID3 (PROGID,XSIZE+20.,11.0,1.0)
   ENDIF
!
   REWIND LFILE
   CALL SKIPFL (LSKIPF,LFILE,IEOF)
   IEOF = LSKIPF
   CALL BUFIN (LFILE,LEOF,PLTHDR(1),NFHDRF)
!
   ICNTNM = MOD(IXSCNT,I_10)
   IXSECT = IXSCNT/10
!
   WRITE (COPT,950) XID(10)
   IFUNCT = 0
   IF (COPT.EQ.CDIF) IFUNCT = 1
   IF (COPT.EQ.CRAT) IFUNCT = 2
   IF (ISCAN.GE.1000) THEN
      ISCANT = MOD(ISCAN,I_100)
      IF (ISCANT.LE.0.OR.SCNID.EQ.-99.) ISCAN = ISCAN-ISCANT
   ELSE
      IF (ISCAN.LE.0.OR.SCNID.EQ.-99.) ISCAN = 0
   ENDIF
!
   YSIZE = YSIZ*SCALE
   YPL = (10.0-YSIZE)/2.
   IF (YPL.GT.0.) GO TO 20
   WRITE (IPR,955)
   YSIZE = 10.
20 HGT = 0.140*SCALE
   IF (NOAXES.EQ.0) DX = (V2-V1)/XSIZE
   DIVLNX = DELV/DX
   NUMDVX = (V2-V1)/DELV+0.01
   SFY = 1.
   IF (JPLOT.EQ.0.AND.JEMIT.EQ.1.AND.IFUNCT.NE.2) THEN
!
!     PLANCK BLACK BODY (INPUT TEMPERATURE - PLOT RADIANCE)
!
      IF (YMAX.GE.2.0) THEN
         CALL BBSCLE
      ELSE
         IF (LOGPLT.GT.0) DELY = 1.
         NUMDVY = (YMAX-YMIN)/DELY+0.01
      ENDIF
      IF (LOGPLT.EQ.0) THEN
         NMAX = LOG10(ABS(YMAX))
         SFY = 10.**(-NMAX+1)
      ENDIF
   ELSE
      IF (LOGPLT.GT.0) DELY = 1.
      NUMDVY = (YMAX-YMIN)/DELY+0.01
   ENDIF
   IF (NOAXES.EQ.0) DY = (YMAX-YMIN)/YSIZE
   DIVLNY = DELY/DY
   WRITE (IPR,937)
   WRITE (IPR,940) YMIN,YMAX,YSIZ,DELY,NUMSBY,NOENDY,IDEC,JEMIT,     &
   &                JPLOT,LOGPLT,JHDR,JOUT,JPLTFL
   IF (IXDEC.LT.1) IXDEC = -1
   IF (IDEC.LT.1) IDEC = -1
!
   IF (JEMIT.EQ.2.AND.IEMIT.EQ.1.AND.ISCAN.GE.1) THEN
      JEMIT = 1
      IEMIT = 2
   ENDIF
   IF (JOUT.GE.1) THEN
      IF (JOUT.LE.2) THEN
         V1S = V1V
         V2S = V2V
         DVS = DVT
         ISCANS = ISCAN
         IEMITS = IEMIT
         IF (JEMIT.EQ.1.AND.JPLOT.EQ.1.AND.IEMIT.GE.1) IEMIT = 2
         V1V = V1
         V2V = V2
         ISCAN = ISCAN+1000
         IF (I4P.EQ.1) DVT = DVT/4.
         CALL BUFOUT (JPLTFL,PLTHDR(1),NFHDRF)
         V1V = V1S
         V2V = V2S
         DVT = DVS
         ISCAN = ISCANS
         IEMIT = IEMITS
      ELSE
         WRITE (JPLTFL,960) XID,(YID(M),M=1,2)
         WRITE (JPLTFL,962) LAYR1,NLAYER
         WRITE (JPLTFL,965) SEC,P0,T0,DVT,V1V,V2V
         WRITE (JPLTFL,967) HOTHER,WBROAD,(HMOL(I),W(I),I=1,NMOL)
         WRITE (JPLTFL,970) IHIRAC,ILBLF4,ICNTNM,IXSECT,IAERSL,      &
            IEMIT,ISCAN,IPLOT,IPATHL,JRAD,SCNID, HWHM
         IF (JEMIT.EQ.0.AND.JPLOT.EQ.0) WRITE (JPLTFL,972)
         IF (JEMIT.EQ.0.AND.JPLOT.EQ.1) WRITE (JPLTFL,975)
         IF (JEMIT.EQ.1.AND.JPLOT.EQ.0) WRITE (JPLTFL,977)
         IF (JEMIT.EQ.1.AND.JPLOT.EQ.1) WRITE (JPLTFL,980)
      ENDIF
   ENDIF
   WRITE (IPR,982) XID,(YID(M),M=1,2)
   WRITE (IPR,962) LAYR1,NLAYER
   WRITE (IPR,965) SEC,P0,T0,DVT,V1V,V2V
   WRITE (IPR,967) HOTHER,WBROAD,(HMOL(I),W(I),I=1,NMOL)
   WRITE (IPR,970) IHIRAC,ILBLF4,ICNTNM,IXSECT,IAERSL,IEMIT,ISCAN,   &
   &                IPLOT,IPATHL,JRAD,SCNID,HWHM
   WRITE (IPR,982)
   ISCANT = MOD(ISCAN,I_100)
   IF (ISCANT.EQ.0) GO TO 30
   JEMSCN = SCNID/100.
   IF (JEMIT.EQ.JEMSCN) GO TO 30
   WRITE (IPR,985)
   GO TO 10
30 CONTINUE
   IF (JOUT.NE.1.AND.JOUT.NE.3.AND.NOAXES.EQ.0) THEN
      IF (JHDR.EQ.0) THEN
         CALL HEADER
      ELSE
         CALL PLOT (1.0,YPL,-3)
      ENDIF
   ENDIF
!
   CALL AXES (IGO,IPR,IREJ,SFY)
   IF (IREJ.EQ.1) GO TO 230
!
   J = 1
   IF (I4P.EQ.1) J = 2
   TIMLOP = 0.
   TIMLIN = 0.
   ISCANT = MOD(ISCAN,I_100)
   IF (ISCANT.EQ.0) GO TO 40
   IF (IGO.EQ.7.OR.IGO.EQ.8) THEN
      WRITE (IPR,985)
      GO TO 230
   ENDIF
   IF (IGO.GE.25.AND.IGO.LE.30) THEN
      WRITE (IPR,985)
      GO TO 230
   ENDIF
   IF (IGO.NE.12.AND.MOD(IGO,I_3).EQ.3) THEN
      WRITE (IPR,985)
      GO TO 230
   ENDIF
!
40 CALL BUFIN (LFILE,LEOF,PNLHDR(1),NPHDRF)
   IF (LEOF.LE.0) GO TO 220
   IF (V2P.GE.V1) GO TO 50
   CALL BUFIN (LFILE,LEOF,DUM(1),2)
   IF (ISCAN.GE.1) GO TO 40
   IF (MOD(IGO+1,I_3).NE.0) GO TO 40
   CALL BUFIN (LFILE,LEOF,DUM(1),2)
   GO TO 40
50 NLO = J
   NHI = NLIM+J-1
   WRITE (IPR,987) V1P,V2P,DV,NLIM,NLO,NHI
   IF (ISCAN.GE.1) GO TO 60
   IF (MOD(IGO+1,I_3).NE.0) GO TO 60
   IF (MOD(IGO+1,I_6).EQ.0) GO TO 60
   CALL BUFIN (LFILE,LEOF,DUM(1),2)
60 CALL BUFIN (LFILE,LEOF,Y(NLO),NLIM)
   IF (ISCAN.GE.1) GO TO 70
   IF (MOD(IGO+1,I_3).NE.0) GO TO 70
   IF (MOD(IGO+1,I_6).NE.0) GO TO 70
   CALL BUFIN (LFILE,LEOF,DUM(1),2)
!
70 JMIN = J
   JMAX = J
   NEWPTS = 0
   NST = NLO
   IF (V1P.LT.V1) NST = (V1-V1P)/DV+0.5+ REAL(NLO)
   NND = NHI
   IF (V2P.GT.V2) NND = (V2-V1P)/DV+0.5+ REAL(NLO)
   CALL CPUTIM (TIME0)
   DO 80 I = NST, NND
      XXI = V1P+DV* REAL(I-NLO)
!
!     YI=Y(I)
!   J RANGES FROM NLO TO NLO-1+NEWPTS (=NPTS).
!
      XX(J) = XXI
      J = J+1
80 END DO
   J = JMIN
!
   ISCANT = MOD(ISCAN,I_100)
!
   IF (ISCANT.EQ.0)                                                  &
   &    GO TO (100,90,210,210,90,210,90,110,210,210,120,90,160,170,   &
   &           210,210,210,210,140,130,210,210,130,210,130,150,210,   &
   &           210,210,210,180,190,210,210,210,210) IGO
   IF (ISCANT.GE.1)                                                  &
   &    GO TO (90,90,210,210,90,210,210,210,210,210,120,90,170,170,   &
   &           210,210,210,210,130,130,210,210,130,210,210,210,210,   &
   &           210,210,210,190,190,210,210,210,210) IGO
90 CALL LINT (NST,NND,J)
   GO TO 200
100 CALL EXPT (NST,NND,J)
   GO TO 200
110 CALL XNTLOG (NST,NND,J)
   GO TO 200
120 CALL TEMPFN (NST,NND,J)
   GO TO 200
130 CALL TENLOG (NST,NND,J)
   GO TO 200
140 CALL LINSTC (NST,NND,J,XLOGCN)
   GO TO 200
150 CALL XLOGLN (NST,NND,J)
   GO TO 200
160 CALL DBOD (NST,NND,J)
   GO TO 200
170 CALL DBTR (NST,NND,J)
   GO TO 200
180 CALL LDBOPD (NST,NND,J)
   GO TO 200
190 CALL LDBTR (NST,NND,J)
200 NEWPTS = J-NLO
   CALL MNMX (NST,NND,JMIN,JMAX)
   IF (NEWPTS.EQ.0) GO TO 210
   NPTS = NLO-1+NEWPTS
   WRITE (IPR,990) NPTS,XX(JMIN),YY(JMIN),XX(JMAX),YY(JMAX)
   IF (IFIRST) THEN
      YYMIN=YY(JMIN)
      XXMIN=XX(JMIN)
      YYMAX=YY(JMAX)
      XXMAX=XX(JMAX)
      IFIRST=.FALSE.
   ELSE
      IF (YYMIN.GT.YY(JMIN)) THEN
         YYMIN=YY(JMIN)
         XXMIN=XX(JMIN)
      ENDIF
      IF (YYMAX.LT.YY(JMAX)) THEN
         YYMAX=YY(JMAX)
         XXMAX=XX(JMAX)
      ENDIF
   ENDIF
   CALL CPUTIM (TIME1)
   TIMLOP = TIMLOP+TIME1-TIME0
   IF (I4P.EQ.0) CALL FSCLIN (NLO,NPTS,XX,YY)
   IF (I4P.EQ.1) CALL FPLINE (NLO,NPTS)
   CALL CPUTIM (TIME2)
   TIMLIN = TIMLIN+TIME2-TIME1
   J = 2
   IF (I4P.EQ.1) J = 4
210 IF (V2P.LT.V2) GO TO 40
!
220 X3 = 3.
   WRITE (IPR,992) XXMIN,YYMIN,XXMAX,YYMAX
   WRITE (IPR,995) TIMLOP,TIMLIN
230 CONTINUE
   IF (JOUT.EQ.1.OR.JOUT.EQ.2) CALL ENDFIL (JPLTFL)
   IF (JOUT.GE.1) CLOSE (JPLTFL)
   GO TO 10
!
240 RETURN
!
900 FORMAT (A60,18X,A2)
902 FORMAT ('1'/30X,A60)
905 FORMAT (4F10.4,4I5,F10.3,I2,I3,I5)
906 FORMAT (A25)
907 FORMAT (/,' ILLEGAL VALUE FOR IRDOPT =',I2,/,' INDICATES THAT ',  &
   &        'ONLY FILE OF THE TWO REQUIRED FILES WAS NAMED. ',/       &
   &        ' BOTH FILES MUST EITHER BE NAMED, OR BOTH NOT NAMED.'/,  &
   &        ' NAME OF JFILE =',A25,' NAME OF LFILE =',A25,/)
910 FORMAT (2G10.4,2G10.3,6I5,2(I2,I3))
911 FORMAT (A25)
915 FORMAT ('0',10X,'INVALID VALUE FOR JOUT = ',I5,', RESET TO ZERO')
920 FORMAT (//,'  *** OVERLAY FEATURE REQUIRES A STANDARD PLOT *** ', &
   &        /,'  *** ON WHICH THE LINE CAN BE OVERLAYED ||    *** ')
927 FORMAT (//,'  --- OVERLAY PLOT WILL HAVE AN OFFSET OF ',G10.3,/,  &
   &        '      WITH JCAL = ',I5,' AND LCAL = ',I5,/)
930 FORMAT (A4,I2.2,'                   ')
932 FORMAT ('0'/8X,'V1        V2     XSIZE      DELV  SBX  NDX LFIL', &
   &        ' SKPF     SCALE IOPT  I4P XDEC')
935 FORMAT (4F10.4,4I5,F10.3,3I5)
937 FORMAT ('0       YMIN       YMAX      YSIZE       DELY  SBY  ',   &
   &        'NDY IDEC  JEM JPLT  LOG JHDR JOUT PFIL')
940 FORMAT (1X,1p,2G11.4,2G11.3,0p,9I5)
945 FORMAT (3A10)
950 FORMAT (A8)
955 FORMAT ('     YSIZE TOO LARGE, HAS BEEN RESET TO 10 INCHES')
960 FORMAT ('1'/10X,10A8,2X,2(1X,A8,1X))
962 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)
965 FORMAT ('0',10X,'SECANT =',F15.5,/,'0',10X,'PRESS(MB) =',F12.5,/, &
   &        '0',10X,'TEMP =',F8.2,/,'0',10X,'DV =',F12.8,' CM-1',/,   &
   &        '0',10X,'V1 =',F12.6,' CM-1',/,'0',10X,'V2 =',F12.6,      &
   &        ' CM-1',/,'0',10X)
967 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//(5X,A6,' = ',       &
   &        1PE10.3))
970 FORMAT ('0',/4X,'IHIRAC    ILBLF4    ICNTNM    IXSECT    IAERSL', &
   &        '     IEMIT     ISCAN     IPLOT    IPATHL      JRAD    ', &
   &        'JFNVAR      HWHM'/10I10,F10.0,F10.5)
972 FORMAT ('0',5X,'WAVENUMBER',7X,' TRANSMISSION',/)
975 FORMAT ('0',5X,'WAVENUMBER',7X,'OPTICAL DEPTH',/)
977 FORMAT ('0',5X,'WAVENUMBER',7X,'   RADIANCE  ',/)
980 FORMAT ('0',5X,'WAVENUMBER',7X,' TEMPERATURE ',/)
982 FORMAT ('0'/10X,10A8,2X,2(1X,A8,1X))
985 FORMAT ('     RESULT FROM SCANNING FUNCTION INCONSISTENT WITH',   &
   &        ' PLOT REQUEST')
987 FORMAT (10X,' * PANEL *',F10.3,F10.3,F20.6,3I10)
990 FORMAT ('0',9X,'NO. OF PTS =',I5,',    FREQ AT MIN =',F10.3,      &
   &        '  MIN VALUE =',1PE11.3,',    FREQ AT MAX =',0PF10.3,     &
   &        '  MAX VALUE =',1PE11.3/)
992 FORMAT ('0',6X,'**** PLOT TOTALS ****    FREQ AT MIN =',F10.3,    &
   &        '  MIN VALUE =',1PE11.3,',    FREQ AT MAX =',0PF10.3,     &
   &        '  MAX VALUE =',1PE11.3/)
995 FORMAT (10X,'TIME FOR LOOP =',F8.3,' SECONDS, TIME FOR LINE =',   &
   &        F8.3)
!
end subroutine PLTLBL
SUBROUTINE FILOPT (V1,V2,JFILE,JSKIPF,LFILE,LSKIPF,MFILE,LENGTH,  &
&                   IOPT)
!
   IMPLICIT REAL*8           (V)
!
!----------------------------------------------------------------------
!             R.D. WORSHAM
!             ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.
!             JANUARY 1990
!----------------------------------------------------------------------
!
!     THIS SUBROUTINE WILL DIFFERENCE OR RATIO TWO
!          LBLRTM OUTPUT OR PLOT FILES
!
!          FILES MUST BE    1) UNFORMATTED LBLRTM FILES
!
!                           2) SINGLE QUANTITY FILES
!                              I.E. CONTAIN ONLY ONE OF THE FOLLOWING:
!
!                                    A) OPTICAL DEPTHS
!                                    B) TRANSMITTANCE
!                                    C) RADIANCE
!                                    D) TEMPERATURE
!
!                           3) THE DIFFERENCE IN DV BETWEEN THE FILES
!                              MUST BE LESS THAN 1.E-8
!
!               *** NOTE: STANDARD LBLRTM OUTPUT FILES USUALLY
!                         CONTAIN BOTH TRANSMITTANCE AND RADIANCE
!                         THE USER MUST CREATE EITHER SCANNED OR
!                         PLOT FILES FOR USE AS INPUT TO THIS ROUTINE
!
!----------------------------------------------------------------------
!
   COMMON XX(2450),YY(2450)
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
   character*8      XID,    HMOL,      YID,HDATE,HTIME
   real*8               SEC,     XALTZ
!
   COMMON /PLTHDR/ XID(10),SEC,P0  ,T0  ,HMOL(60),XALTZ(4),          &
   &                W(60),PZL,PZU,TZL,TZU,WBROAD,DVT,V1V,V2V,TBOUND,  &
   &                EMISIV,FSCDID(17),NMOL,NLAYER,YID1,YID(10),LSTWDF
   COMMON /JYCOM/ V1PJ,V2PJ,DVJ,NLIMJ,YJ(2502)
   COMMON /LYCOM/ V1PL,V2PL,DVL,NLIML,YL(2502)
   DIMENSION PLTHDR(2),PNLHDJ(2),PNLHDL(2),COPT(2)
!
   EQUIVALENCE (PLTHDR(1),XID(1)) , (PNLHDJ(1),V1PJ),                &
   &            (PNLHDL(1),V1PL) , (FSCDID(1),IHIRAC),                &
   &            (FSCDID(2),ILBLF4) , (FSCDID(3),IXSCNT),              &
   &            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(6),ISCAN) , (FSCDID(7),IPLOT),                &
   &            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                &
   &            (FSCDID(10),ITEST) , (FSCDID(11),IMRG),               &
   &            (FSCDID(12),SCNID) , (FSCDID(13),HWHM),               &
   &            (FSCDID(14),IDABS) , (FSCDID(15),IATM),               &
   &            (FSCDID(16),LAYR1) , (FSCDID(17),NLAYFS),             &
   &            (YID(1),HDATE) , (YID(2),HTIME)
!
   LOGICAL IFIRST
   CHARACTER COPT*8,HOTHER*6
!
   DATA COPT / 'DIFRENCE',' RATIO  '/
   DATA HOTHER / ' OTHER'/
!
   DATA I_0/0/, I_1/1/, I_3/3/, I_6/6/, I_10/10/, I_100/100/
!
   REWIND JFILE
   REWIND LFILE
   CALL SKIPFL (JSKIPF,JFILE,IJEOF)
   CALL SKIPFL (LSKIPF,LFILE,ILEOF)
   IJEOF = JSKIPF
   ILEOF = LSKIPF
   WRITE (IPR,900) JFILE
   CALL BUFIN (JFILE,JEOF,PLTHDR(1),NFHDRF)
!
   ICNTNM = MOD(IXSCNT,I_10)
   IXSECT = IXSCNT/10
!
   WRITE (IPR,905) XID,(YID(M),M=1,2)
   WRITE (IPR,910) LAYR1,NLAYER
   WRITE (IPR,915) SEC,P0,T0,DVT,V1V,V2V
   WRITE (IPR,920) HOTHER,WBROAD,(HMOL(I),W(I),I=1,NMOL)
   ISCAN = MAX(ISCAN,I_0)
   WRITE (IPR,925) IHIRAC,ILBLF4,ICNTNM,IXSECT,IAERSL,IEMIT,ISCAN,   &
   &                IPLOT,IPATHL,JRAD,SCNID,HWHM
   WRITE (IPR,905)
!
   WRITE (IPR,930) LFILE
   CALL BUFIN (LFILE,LEOF,PLTHDR(1),NFHDRF)
!
   ICNTNM = MOD(IXSCNT,I_10)
   IXSECT = IXSCNT/10
!
   READ (COPT(IOPT-1),935) XID(10)
   CALL BUFOUT (MFILE,PLTHDR(1),NFHDRF)
   WRITE (IPR,905) XID,(YID(M),M=1,2)
   WRITE (IPR,910) LAYR1,NLAYER
   WRITE (IPR,915) SEC,P0,T0,DVT,V1V,V2V
   WRITE (IPR,920) HOTHER,WBROAD,(HMOL(I),W(I),I=1,NMOL)
   ISCAN = MAX(ISCAN,I_0)
   WRITE (IPR,925) IHIRAC,ILBLF4,ICNTNM,IXSECT,IAERSL,IEMIT,ISCAN,   &
   &                IPLOT,IPATHL,JRAD,SCNID,HWHM
   WRITE (IPR,905)
   IF (IOPT.EQ.2) WRITE (IPR,940)
   IF (IOPT.EQ.3) WRITE (IPR,945)
   J = 0
   L = 0
   M = 0
   LIMOUT = 2400
   IFIRST = .TRUE.
   NLOW = 1
   NLIMJ = 0
   NLIML = 0
!
10 IF (J.GE.NLIMJ) THEN
      CALL BUFIN (JFILE,JEOF,PNLHDJ(1),NPHDRF)
      WRITE (IPR,950) V1PJ,V2PJ,DVJ,NLIMJ,JFILE,JEOF
      IF (JEOF.LE.0.AND.M.NE.0) GO TO 50
      IF (JEOF.LE.0) GO TO 60
      CALL BUFIN (JFILE,JEOF,YJ(NLOW),NLIMJ)
      J = 0
   ENDIF
   IF (L.GE.NLIML) THEN
      CALL BUFIN (LFILE,LEOF,PNLHDL(1),NPHDRF)
      WRITE (IPR,950) V1PL,V2PL,DVL,NLIML,LFILE,LEOF
      IF (LEOF.LE.0.AND.M.NE.0) GO TO 50
      IF (LEOF.LE.0) GO TO 60
      CALL BUFIN (LFILE,LEOF,YL(NLOW),NLIML)
      L = 0
   ENDIF
!
!     CHECK DV TO INSURE SAME POINT SPACING
!
   IF (IFIRST) THEN
      V1TST = DVJ*1.E-3
      DVTST = V1TST/2400.
      DVDIF = DVJ-DVL
      IDVTST = 0
      IF (ABS(DVDIF).GT.DVTST) THEN
         IDVTST = 1
         WRITE (IPR,955) DVJ,DVL,DVTST
      ELSE
!
!     IF DV IS OKAY, CHECK STARTING POINTS FOR BEATING
!
         NTEST = (V1PJ-V1PL)/DVJ
         V1NEW = V1PL+ REAL(NTEST)*DVJ
         V1DIF = V1NEW-V1PJ
         IF (ABS(V1DIF).GT.V1TST) THEN
            IDVTST = 1
            WRITE (IPR,960) V1PJ,V1PL,V1TST
         ENDIF
      ENDIF
      IF (IDVTST.EQ.1) STOP ' FILOPT - PANELS DO NOT MATCH '
   ENDIF
   IF (V1PJ.NE.V1PL) THEN
      IF (V1PJ.GT.V1PL) THEN
         VDIF = (V1PJ-V1PL)/DVJ+0.5
         L = VDIF
         IF (L.GE.NLIML) THEN
            NLIML = 0
            GO TO 10
         ENDIF
      ELSE
         VDIF = (V1PL-V1PJ)/DVJ+0.5
         J = VDIF
         IF (J.GE.NLIMJ) THEN
            NLIMJ = 0
            GO TO 10
         ENDIF
      ENDIF
   ENDIF
   IF (V2PJ.GE.V1) GO TO 20
   NLIMJ = 0
20 IF (V2PL.GE.V1.AND.NLIMJ.NE.0) GO TO 30
   NLIML = 0
   GO TO 10
!
30 NST = 0
   IF (V1PJ.LT.V1) NST = (V1-V1PJ)/DVJ+0.5+ REAL(NLOW)
   J = J+NST
   L = L+NST
40 CONTINUE
   J = J+1
   IF (J.GT.NLIMJ) GO TO 10
   L = L+1
   IF (L.GT.NLIML) THEN
      J = J-1
      GO TO 10
   ENDIF
   M = M+1
!
   XX(M) = V1PJ+DVJ* REAL(J-1)
!
!     FOR IOPT = 2 - DIFFERENCE VALUES
!     FOR IOPT = 3 - RATIO VALUES
!
   IF (IOPT.EQ.2) THEN
      YY(M) = YJ(J)-YL(L)
   ELSE
      IF (YL(L).NE.0) THEN
         YY(M) = YJ(J)/YL(L)
      ELSE
         YY(M) = 0.0
      ENDIF
   ENDIF
   IF (IFIRST) THEN
      YMIN = YY(M)
      YMAX = YY(M)
      IFIRST = .FALSE.
   ELSE
      YMIN = MIN(YMIN,YY(M))
      YMAX = MAX(YMAX,YY(M))
   ENDIF
   IF (XX(M).GE.V2) GO TO 50
   IF (J.EQ.NLIMJ) GO TO 50
   IF (M.LE.LIMOUT) GO TO 40
!
50 NPTS = M
   V1PJS = V1PJ
   V2PJS = V2PJ
   NLIMJS = NLIMJ
   V1PJ = XX(1)
   V2PJ = XX(NPTS)
   NLIMJ = NPTS
   CALL BUFOUT (MFILE,PNLHDJ(1),NPHDRF)
   CALL BUFOUT (MFILE,YY(NLOW),NPTS)
   M = 0
   IF (V2PJ.LT.V2) THEN
      V1PJ = V1PJS
      V2PJ = V2PJS
      NLIMJ = NLIMJS
      GO TO 40
   ENDIF
!
60 WRITE (IPR,965) YMIN,YMAX
   CALL ENDFIL (MFILE)
!
   RETURN
!
900 FORMAT (/,'0',4X,'FILOPT ***** READING FROM FILE',I3,' *****')
905 FORMAT ('0',/,10X,10A8,2X,2(1X,A8,1X))
910 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)
915 FORMAT ('0',10X,'SECANT =',F15.5/'0',10X,'PRESS(MB) =',F12.5,/,   &
   &        '0',10X,'TEMP =',F8.2/'0',10X,'DV =',F12.8,' CM-1',/,     &
   &        '0',10X,'V1 =',F12.6,' CM-1'/'0',10X,'V2 =',F12.6,        &
   &        ' CM-1',/,'0',10X)
920 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//(5X,A6,' = ',       &
   &        1PE10.3))
925 FORMAT ('0',/4X,'IHIRAC    ILBLF4    ICNTNM    IXSECT    ',       &
   &        'IAERSL     IEMIT     ISCAN     IPLOT    IPATHL      ',   &
   &        'JRAD    JFNVAR      HWHM',/,10(I10),F10.0,F10.5)
930 FORMAT (/,'1',4X,'FILOPT ***** READING FROM FILE',I3,' *****')
935 FORMAT (A8)
940 FORMAT (2(/),5X,'-----------------------------------------',/,    &
   &        5X,'OPT = 1 * (FILE 1)-(FILE 2) ** DIFFERENCE',/,         &
   &        5X,'-----------------------------------------')
945 FORMAT (2(/),5X,'------------------------------------',/,         &
   &        5X,'OPT = 2 * (FILE 1)/(FILE 2) ** RATIO',/,              &
   &        5X,'------------------------------------')
950 FORMAT (/,5X,'* PANEL *   V1P = ',F10.3,' V2P = ',F10.3,          &
   &        ' DVP = ',F13.6,' NLIM = ',I5,' IFILE = ',I3,             &
   &        ' IEOF = ',I4)
955 FORMAT ('  PANELS DO NOT MATCH FOR INPUT FILES ',/,'  DVJ = ',    &
   &        F20.10,'  DVL = ',F20.10,' DVTST = ',F14.10)
960 FORMAT ('  PANELS DO NOT MATCH FOR INPUT FILES ',/,' V1PJ = ',    &
   &        F20.10,' V1PL = ',F20.10,' V1TST = ',F14.10)
965 FORMAT (/,6X,'MINIMUM Y VALUE = ',1PE13.6,5X,                     &
   &        'MAXIMUM Y VALUE = ',1PE13.6)
!
end subroutine FILOPT
SUBROUTINE HEADER
!
   IMPLICIT REAL*8           (V)
!
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES
!
   character*8      XID,       HMOL,        YID
   real*8               SEC,          XALTZ
!
   COMMON /PLTHDR/ XID(10),SEC,P0  ,T0  ,HMOL(60),XALTZ(4),          &
   &                W(60),PZL,PZU,TZL,TZU,WBROAD,DVT,V1V,V2V,TBOUND,  &
   &                EMISIV,FSCDID(17),NMOL,NLAYER,YI1,YID(10),LSTWDF
!
   CHARACTER HTEN*2,HSE*5,HPR*11,HTP*10,HDV*11,HV1*11,HV2*11
   CHARACTER HSC*13,HVA*11,HHW*11,HLR*11,H4C*10,HAMNTS*19,HOTHER*7
   CHARACTER*10 HSLIT(0:4)
   CHARACTER XID1*72,XID2*30
!
   EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),              &
   &            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),              &
   &            (FSCDID(5),IEMIT) , (FSCDID(6),ISCAN),                &
   &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),               &
   &            (FSCDID(9),JRAD) , (FSCDID(10),ITEST),                &
   &            (FSCDID(11),IMRG) , (FSCDID(12),SCNID),               &
   &            (FSCDID(13),HWHM) , (FSCDID(14),IDABS),               &
   &            (FSCDID(15),IATM) , (FSCDID(16),LAYR1),               &
   &            (FSCDID(17),NLAYFS) , (YID(1),HDATE),                 &
   &            (YID(2),HTIME) , (YI1,IMULT)
!
   DATA HSLIT / ' RECTANGLE',' TRIANGLE ','  GAUSS   ','  SINCSQ  ', &
   &             '   SINC   '/
   DATA XP1 / 0.1 /, XP15 / 0.15 /, X1 / 1.0 /, X2 / 2.0 /,          &
   &     X26 / 2.6 /, X27 / 2.7 /, X3 / 3.0 /, X6 / 6.0 /,            &
   &     X7 / 7.0 /, X8 / 8.0 /, X8P3 / 8.3 /, X11 / 11. /,           &
   &     YP3 / 0.3 /, YSHF / 0.45 /, Y10 / 10.0 /
   DATA HTEN / '10'/,                     HSE / 'SEC ='/,            &
   &     HPR / 'PRESS (MB)='/,             HTP / 'TEMP (K) ='/,       &
   &     HDV / 'DV(CM-1) = '/,             HV1 / 'V1(CM-1) = '/,      &
   &     HV2 / 'V2(CM-1) = '/,             HSC / 'SCAN FUNCTION'/,    &
   &     HVA / ' (VARIABLE)'/,             HHW / 'HWHM(CM-1)='/,      &
   &     HLR / 'NLAYERS  = '/,             H4C / 'H4CXAEISRM'/,       &
   &     HAMNTS / 'AMOUNTS (MOL/CM**2)'/,  HOTHER / ' OTHER '/
!
   DATA I_10/10/, I_100/100/
!
   ICNTNM = MOD(IXSCNT,I_10)
   IXSECT = IXSCNT/10
!
   WRITE (XID1,'( 9(A8   ))') (XID(I),I=1,9)
   WRITE (XID2,'(3(A10  ))') XID(10),(YID(I),I=1,2)
   YT = Y10-YSHF
   CALL SYMBOL (0.0,Y10,XP15,XID1,0.0,72)
   CALL SYMBOL (0.0,YT,XP15,XID2,0.0,30)
   YT = YT-2.*YSHF
   CALL SYMBOL (X1,YT,XP15,HSE,0.0,5)
   RSEC = SEC
   CALL NUMBER (X2,YT,XP15,RSEC,0.0,1)
   YT = YT-YSHF
   CALL SYMBOL (X1,YT,XP15,HPR,0.0,11)
   CALL NUMBER (X3,YT,XP15,P0,0.0,3)
   YT = YT-YSHF
   CALL SYMBOL (X1,YT,XP15,HTP,0.0,10)
   CALL NUMBER (X27,YT,XP15,T0,0.0,3)
   YT = YT-YSHF
   CALL SYMBOL (X1,YT,XP15,HDV,0.0,11)
   CALL NUMBER (X26,YT,XP15,DVT,0.0,5)
   YT = YT-YSHF
   CALL SYMBOL (X1,YT,XP15,HV1,0.0,11)
   RV1V = V1V
   CALL NUMBER (X26,YT,XP15,RV1V,0.0,2)
   YT = YT-YSHF
   CALL SYMBOL (X1,YT,XP15,HV2,0.0,11)
   RV2V = V2V
   CALL NUMBER (X26,YT,XP15,RV2V,0.0,2)
   ISCANT = MOD(ISCAN,I_100)
   IF (ISCANT.EQ.0) GO TO 20
   JFNVAR = SCNID
   JEMIT = JFNVAR/100
   JFN = (JFNVAR-100*JEMIT)/10
   JVAR = JFNVAR-100*JEMIT-10*JFN
   YT = YT-YSHF
   CALL SYMBOL (X1,YT,XP15,HSC,0.0,13)
   IF (JVAR.EQ.0) GO TO 10
   YT = YT-YP3
   CALL SYMBOL (X1,YT,XP15,HVA,0.0,11)
10 YT = YT-YP3
!
   CALL SYMBOL (X1,YT,XP15,HSLIT(JFN),0.0,10)
   YT = YT-YP3
   CALL SYMBOL (X1,YT,XP15,HHW,0.0,11)
!
   CALL NUMBER (X27,YT,XP15,HWHM,0.0,5)
20 YT = YT-YSHF
   CALL SYMBOL (X1,YT,XP15,HLR,0.0,11)
   XLAYER = NLAYER
   CALL NUMBER (X26,YT,XP15,XLAYER,0.0,-1)
   YT = YT-YSHF
   CALL SYMBOL (X1,YT,XP15,H4C,0.0,10)
   YT = YT-YP3
   JRAD = JRAD+4
   CONTRL = JRAD+10*(ISCANT+100*(IEMIT+10*(IAERSL+10*(IXSECT+10*     &
   &                 (ICNTNM+10*(ILBLF4+10*IHIRAC))))))
   CONTRL = CONTRL*10+IMULT
   CALL NUMBER (X1,YT,XP15,CONTRL,0.0,-1)
   YT = YT-YP3
   IF (IATM.NE.0) CALL LBH1H2 (YT,X1,YID)
   YT = YT-YP3
   IF (IAERSL.NE.0) CALL LBAERS (YT,X1,YID)
   YT = Y10-2.*YSHF
   CALL SYMBOL (X6,YT,XP15,HAMNTS,0.0,19)
   YT = YT-YSHF
   YTE = YT+XP1
   RW = 0.
   HIP = 0.
   IF (WBROAD.LT.1.E-23) GO TO 30
   IP =  LOG10(WBROAD)
   HIP = IP
   RW = WBROAD/10.**IP
30 CALL SYMBOL (X6,YT,XP15,HOTHER,0.0,7)
   CALL NUMBER (X7,YT,XP15,RW,0.0,3)
   CALL SYMBOL (X8,YT,XP15,HTEN,0.0,2)
   CALL NUMBER (X8P3,YTE,XP1,HIP,0.0,-1)
   DO 40 M = 1, NMOL
      NM = M
      IF (W(M).LE.0.) GO TO 40
      YT = YT-YSHF
      IF (YT.LT.1.0) GO TO 50
      YTE = YT+XP1
      RW = 0.
      HIP = 0.
      IP = LOG10(W(M))
      HIP = IP
      RW = W(M)/10.**IP
      CALL SYMBOL (X6,YT,XP15,HMOL(M),0.0,6)
      CALL NUMBER (X7,YT,XP15,RW,0.0,3)
      CALL SYMBOL (X8,YT,XP15,HTEN,0.0,2)
      CALL NUMBER (X8P3,YTE,XP1,HIP,0.0,-1)
40 END DO
   CALL PLOT (X11,0.0,-3)
   CALL PLOT (1.0,YPL,-3)
!
   RETURN
!
50 YT = Y10-2.*YSHF
   X12 = X6+X6
   X13 = X12+X1
   X14 = X13+X1
   X14P3 = X14+0.3
   CALL SYMBOL (X12,YT,XP15,HAMNTS,0.0,19)
   DO 60 M = NM, NMOL
      IF (W(M).LE.0.) GO TO 60
      YT = YT-YSHF
      YTE = YT+XP1
      RW = 0.
      HIP = 0.
      IP = LOG10(W(M))
      HIP = IP
      RW = W(M)/10.**IP
      CALL SYMBOL (X12,YT,XP15,HMOL(M),0.0,6)
      CALL NUMBER (X13,YT,XP15,RW,0.0,3)
      CALL SYMBOL (X14,YT,XP15,HTEN,0.0,2)
      CALL NUMBER (X14P3,YTE,XP1,HIP,0.0,-1)
60 END DO
   XFN = X14+X3
   CALL PLOT (XFN,0.0,-3)
   CALL PLOT (1.0,YPL,-3)
!
   RETURN
!
end subroutine HEADER
SUBROUTINE AXES (IGO,LOUT,IREJ,SFY)
!
   IMPLICIT REAL*8           (V)
!
   character*8      XID,    HMOL,      YID,HDATE,HTIME
   real*8               SEC,     XALTZ
!
   COMMON /PLTHDR/ XID(10),SEC,P0  ,T0  ,HMOL(60),XALTZ(4),          &
   &                W(60),PZL,PZU,TZL,TZU,WBROAD,DVT,V1V,V2V,TBOUND,  &
   &                EMISIV,FSCDID(17),NMOL,NLAYER,YID1,YID(10),LSTWDF
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES
!
   CHARACTER YIDC*43,ASTR*3
!
   CHARACTER TITL1*12,TITL2*13,TITL3*50,TITL4*11,XWAVEN*17,XWAVEL*23
   CHARACTER BLK*2,HSCALE*10,GIGAHZ*15,TITL5*12,TITL6*9,TITL1D*25,   &
   &          TITL1R*20,TITL2D*26,TITL2R*21,TITL3D*34,TITL3R*29,      &
   &          TITL4D*24,TITL4R*19,TITL5D*23,TITL5R*18,TITL6D*21,      &
   &          TITL6R*16
!
   EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),              &
   &            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),              &
   &            (FSCDID(5),IEMIT) , (FSCDID(6),ISCAN),                &
   &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),               &
   &            (FSCDID(9),JRAD) , (FSCDID(10),ITEST),                &
   &            (FSCDID(11),IMRG) , (FSCDID(12),SCNID),               &
   &            (FSCDID(13),HWHM) , (FSCDID(14),IDABS),               &
   &            (FSCDID(15),IATM) , (FSCDID(16),LAYR1),               &
   &            (FSCDID(17),NLAYFS)
!
   DATA NX / -17 /,LX / 23 /,NT1 / 12 /,NT2 / 13 /,NT3 / 37 /,       &
   &     NT4 / 11 /,NGH / 15 /,NT6 / 9 /,NT1D / 25 /,NT1R / 20 /,     &
   &     NT2D / 26 /,NT2R / 21 /,NT3D / 21 /,NT3R / 16 /,NT4D / 24 /, &
   &     NT4R / 19 /,NT5D / 23 /,NT5R / 18 /,NT6D / 21 /,NT6R / 16 /
   DATA XWAVEN / 'WAVENUMBER (CM-1)'/,                               &
   &     XWAVEL / 'WAVELENGTH (MICROMETER)'/,                         &
   &     GIGAHZ / 'FREQUENCY (GHZ)'/,                                 &
   &     BLK / ' '/,                                                  &
   &     TITL1 / 'TRANSMISSION'/,                                     &
   &     TITL2 / 'OPTICAL DEPTH'/,                                    &
   &     TITL3 / ' RADIANCE WATTS /  (CM**2*STER*CM-1)  * '/,         &
   &     TITL4 / 'TEMPERATURE'/,                                      &
   &     TITL5 / ' ABSORPTION '/,                                     &
   &     TITL6 / ' DECIBELS'/
   DATA TITL1D / ' TRANSMISSION DIFFERENCE '/,                       &
   &     TITL1R / ' TRANSMISSION RATIO '/,                            &
   &     TITL2D / ' OPTICAL DEPTH DIFFERENCE '/,                      &
   &     TITL2R / ' OPTICAL DEPTH RATIO '/,                           &
   &     TITL3D / ' RADIANCE DIFFERENCE  * '/,                        &
   &     TITL3R / ' RADIANCE RATIO  * '/,                             &
   &     TITL4D / ' TEMPERATURE DIFFERENCE '/,                        &
   &     TITL4R / ' TEMPERATURE RATIO '/,                             &
   &     TITL5D / ' ABSORPTION DIFFERENCE '/,                         &
   &     TITL5R / ' ABSORPTION RATIO '/,                              &
   &     TITL6D / ' DECIBELS DIFFERENCE '/,                           &
   &     TITL6R / ' DECIBELS RATIO '/
!
   DATA ASTR / '*  '/
   DATA XP1 / 0.1 /,X11 / 11.0 /,YP5 / 0.5 /
   DATA CL / 29.9792458 /
!
   IREJ = 0
   IGO = IEMIT+3*JEMIT+6*JPLOT+18*LOGPLT+1
   IF (JOUT.EQ.1.OR.JOUT.EQ.3.OR.NOAXES.EQ.1) THEN
      GO TO (10,10,140,140,10,140,10,10,140,140,10,10,10,10,140,     &
         140,140,140,10,10,140,140,10,140,10,10,140,140,140, 140,10,10, &
         140,140,140,140) IGO
10    RETURN
   ENDIF
!
   RV1 = V1
   CALL AXISL (0.0,0.0,XWAVEN,NX,NUMDVX,DIVLNX,NUMSBX,RV1,DELV,      &
   &            IXDEC,0.,HGT,1,0,NOENDX,0,0)
   IF (V1.LT.100.) GO TO 40
   CALL AXISL (0.0,YSIZE,BLK,1,NUMDVX,DIVLNX,NUMSBX,RV1,DELV,        &
   &            IXDEC,0.,HGT,0,0,NOENDX,1,0)
!
!   CREATE WAVELENGTH SCALE (MICRONS) OF APPROXIMATELY FIVE VALUES:
!
   IF (YSIZE.GT.8.) GO TO 60
   YSZSH = YSIZE+YP5
   YINCH = YSZSH+HGT
   CALL AXISL (0.0,YSZSH,XWAVEL,LX,1,XSIZE,1,RV1,DELV,-1,0.,HGT,0,0, &
   &            1,0,0)
   TOP = 10000./V1
   STEP = (TOP-10000./V2)/5.
   NDP = - LOG10(STEP)+1.
   ISTEP = STEP*(10.**NDP)
   IF (ISTEP.GE.6) ISTEP = 10+10*((ISTEP-5)/10)
   STEP =  REAL(ISTEP)/(10.**NDP)
   IF (NDP.LE.0) NDP = -1
   TOP = TOP-(int(TOP/STEP))*step
20 BOT = 10000./TOP
   IF (BOT.GT.V2) GO TO 60
   XINCH = (BOT-V1)/DX
   IF (XINCH.GT.XSIZE) GO TO 30
   CALL SYMBOL (XINCH,YSZSH,XP1,'3',0.0,-1)
!
!     ('3' AS FOURTH ARGUMENT REPRESENTS PLUS SYMBOL FOR TIC MARK)
!
   DIGITS = 0.
   IF (TOP.GE.10.) DIGITS =  LOG10(TOP)
   DIGITS = AINT(DIGITS+1.0E-12)
   SIZNUM = HGT*(DIGITS+ REAL(NDP)+1.7)
   XPOS = XINCH-0.5*SIZNUM
   IF (XPOS.GT.XSIZE) GO TO 30
   CALL NUMBER (XPOS,YINCH,HGT,TOP,0.0,NDP)
30 TOP = TOP-STEP
   GO TO 20
!
!     CREATE FREQUENCY SCALE IN GIGAHERTZ
!
40 RV1 = V1
   CALL AXISL (0.0,YSIZE,GIGAHZ,NGH,1,XSIZE,1,RV1,DELV,-1,0.,HGT,0,  &
   &            0,1,0,0)
   IF (DELV.LT.1./CL) GO TO 60
   HLFHGT = 0.5*HGT
   YPOS = YSIZE+HGT+HLFHGT
   DO 50 I = 1, NUMDVX
      IG = (V1+ REAL(I-1)*DELV)*CL
      G = IG
      IF (IG.NE.0) G = IG+1
      XG = ((G/CL)-V1)/DX
      IF (XG.GT.XSIZE) GO TO 50
      CALL PLOT (XG,YSIZE,3)
      CALL PLOT (XG,YSIZE+HLFHGT,2)
      DIGITS = 0.
      IF (G.GE.10.) DIGITS = LOG10(G)
      DIGITS = AINT(DIGITS+1.0E-12)
      SIZNUM = HGT*(DIGITS+0.7)
      XPOS = XG-0.5*SIZNUM
      IF (XPOS.GT.XSIZE) GO TO 50
      CALL NUMBER (XPOS,YPOS,HGT,G,0.0,-1)
50 END DO
!
60 GO TO (70,70,140,140,100,140,90,90,140,140,130,130,110,110,       &
   &       140,140,140,140,70,70,140,140,100,140,90,90,140,140,       &
   &       140,140,120,120,140,140,140,140) IGO
!
70 IF (IDABS.LT.0) GO TO 80
   IF (IGO.LT.13) THEN
      IF (IFUNCT.EQ.0) CALL PLTLIN (TITL1,BLK,NT1,1.0)
      IF (IFUNCT.EQ.1) CALL PLTLIN (TITL1D,BLK,NT1D,1.0)
      IF (IFUNCT.EQ.2) CALL PLTLIN (TITL1R,BLK,NT1R,1.0)
   ELSE
      IF (IFUNCT.EQ.0) CALL PLTLOG (TITL1,BLK,NT1)
      IF (IFUNCT.EQ.1) CALL PLTLOG (TITL1D,BLK,NT1D)
      IF (IFUNCT.EQ.2) CALL PLTLOG (TITL1R,BLK,NT1R)
   ENDIF
   GO TO 150
!
80 IF (IGO.LT.13) THEN
      IF (IFUNCT.EQ.0) CALL PLTLIN (TITL5,BLK,NT1,1.0)
      IF (IFUNCT.EQ.1) CALL PLTLIN (TITL5D,BLK,NT1D,1.0)
      IF (IFUNCT.EQ.2) CALL PLTLIN (TITL5R,BLK,NT1R,1.0)
   ELSE
      IF (IFUNCT.EQ.0) CALL PLTLOG (TITL5,BLK,NT1)
      IF (IFUNCT.EQ.1) CALL PLTLOG (TITL5D,BLK,NT1D)
      IF (IFUNCT.EQ.2) CALL PLTLOG (TITL5R,BLK,NT1R)
   ENDIF
   GO TO 150
!
90 IF (IGO.LT.13) THEN
      IF (IFUNCT.EQ.0) CALL PLTLIN (TITL2,BLK,NT2,1.0)
      IF (IFUNCT.EQ.1) CALL PLTLIN (TITL2D,BLK,NT2D,1.0)
      IF (IFUNCT.EQ.2) CALL PLTLIN (TITL2R,BLK,NT2R,1.0)
   ELSE
      IF (IFUNCT.EQ.0) CALL PLTLOG (TITL2,BLK,NT2)
      IF (IFUNCT.EQ.1) CALL PLTLOG (TITL2D,BLK,NT2D)
      IF (IFUNCT.EQ.2) CALL PLTLOG (TITL2R,BLK,NT2R)
   ENDIF
   GO TO 150
!
100 IF (IGO.GE.13) THEN
      IF (IFUNCT.EQ.0) CALL PLTLOG (TITL3,BLK,NT3)
      IF (IFUNCT.EQ.1) CALL PLTLOG (TITL3D,BLK,NT3D)
      IF (IFUNCT.EQ.2) CALL PLTLOG (TITL3R,BLK,NT3R)
      GO TO 150
   ENDIF
   WRITE (HSCALE,'(1PE10.0)') SFY
!
   NT = 50
   NTD = 34
   NTR = 29
   TITL3(41:50) = HSCALE
   TITL3D(25:34) = HSCALE
   TITL3R(20:29) = HSCALE
   IF (IFUNCT.EQ.0) CALL PLTLIN (TITL3,BLK,NT,SFY)
   IF (IFUNCT.EQ.1) CALL PLTLIN (TITL3D,BLK,NTD,SFY)
   IF (SFY.EQ.1.) NTR = NT3R
   IF (IFUNCT.EQ.2) CALL PLTLIN (TITL3R,BLK,NTR,SFY)
   GO TO 150
!
110 CALL PLTLIN (TITL6,BLK,NT6,1.0)
   GO TO 150
120 CALL PLTLOG (TITL6,BLK,NT6)
   GO TO 150
130 IF (IFUNCT.EQ.0) CALL PLTEMP (TITL4,BLK,NT4)
   IF (IFUNCT.EQ.1) CALL PLTEMP (TITL4D,BLK,NT4D)
   IF (IFUNCT.EQ.2) CALL PLTEMP (TITL4R,BLK,NT4R)
   GO TO 150
140 IREJ = 1
   WRITE (LOUT,900) IGO
!
   RETURN
!
!     WRITE DATE AND TIME TO BOTTOM OF PLOT
!
150 YBIAS = -3.75*1.25*HGT-0.25
   XBIAS = YBIAS-0.5
   CALL LBLDAT(HDATE)
   CALL FTIME (HTIME)
   WRITE (YIDC,'(2(A10),A3,2(A10))') HDATE,HTIME,ASTR,(YID(I),I=1,2)
   CALL SYMBOL (XBIAS,YBIAS,0.10,YIDC,0.0,43)
!
   RETURN
!
900 FORMAT (10X,'PLOT REQUEST INCOMPATIBLE WITH LBLRTM, IGO =',I2)
!
end subroutine AXES
SUBROUTINE BBSCLE
!
   USE phys_consts, ONLY: radcn2
   IMPLICIT REAL*8           (V)
!
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES
!
!
!     PLANCK  BLACK BODY
!
   XKTMX = YMAX/RADCN2
   YMAX1 = PLANCK(V2, XKTMX)
   IF (YMIN.GT.0.) THEN
      XKTMN = YMIN/RADCN2
      YMIN1 = PLANCK(V1, XKTMN)
      RATYF = 0.
   ELSE
      YMIN1 = 0.
      RATYF = YMIN/YMAX
   ENDIF
   NS =  LOG10(YMAX1)
   IF (YMIN1.GT.0.) THEN
      NB = LOG10(YMIN1)
   ELSE
      NB = NS-5
   ENDIF
!
!     LINEAR RADIANCE
!
   IF (LOGPLT.EQ.0) THEN
      YMAX = 10.**NS
      YMIN = 10.**NB
      ZMAX = YMAX1/YMAX
      ZMIN = YMIN1/YMIN
      IF (ZMAX.GT.0.1.AND.ZMAX.LE.0.2) YMAX = .2*YMAX
      IF (ZMAX.GT.0.2.AND.ZMAX.LE.0.5) YMAX = .5*YMAX
!
!        IF(ZMAX.GT.0.5.AND.ZMAX.LE.1.0) YMAX=YMAX
!
      IF (ZMIN.GT.0.1.AND.ZMIN.LE.0.2) YMIN = .1*YMIN
      IF (ZMIN.GT.0.2.AND.ZMIN.LE.0.5) YMIN = .2*YMIN
      IF (ZMIN.GT.0.5.AND.ZMIN.LE.1.0) YMIN = .5*YMIN
      IF (RATYF.NE.0..AND.IFUNCT.EQ.1) YMIN = RATYF*YMAX
      DELY = (YMAX-YMIN)/10.
      NUMDVY = 10
   ELSE
!
!     LOG RADIANCE
!
      YMAX = NS
      IF (NB.EQ.NS) NB = NS-1
      YMIN = NB
      NUMDVY = NS-NB
      DELY = 1.
   ENDIF
!
   RETURN
!
end subroutine BBSCLE
SUBROUTINE PLTLIN (YTITLE,BLK,NY,SFY)
!
   IMPLICIT REAL*8           (V)
!
   character YTITLE,blk
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES
!
   CALL AXISL (0.0,0.0,YTITLE,NY,NUMDVY,DIVLNY,NUMSBY,YMIN*SFY,      &
   &            DELY*SFY,IDEC,90.0,HGT,1,1,NOENDY,0,0)
   CALL AXISL (XSIZE,0.0,BLK,-1,NUMDVY,DIVLNY,NUMSBY,YMIN*SFY,       &
   &            DELY*SFY,IDEC,90.0,HGT,1,1,NOENDY,1,0)
!
   RETURN
!
end subroutine PLTLIN
SUBROUTINE PLTEMP (YTITLE,BLK,NY)
!
   IMPLICIT REAL*8           (V)
!
   character YTITLE,blk
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES
!
   CALL AXISL (0.0,0.0,YTITLE,NY,NUMDVY,DIVLNY,NUMSBY,YMIN,DELY,     &
   &            IDEC,90.0,HGT,1,1,NOENDY,0,0)
   IF (YMIN.LE.0.0) THEN
      CALL AXISL (XSIZE,0.0,BLK,-1,NUMDVY,DIVLNY,NUMSBY,YMIN,DELY,   &
         IDEC,90.0,HGT,1,1,NOENDY,1,0)
   ELSE
      CALL AX2
   ENDIF
!
   RETURN
!
end subroutine PLTEMP
SUBROUTINE PLTLOG (YTITLE,BLK,NY)
!
   IMPLICIT REAL*8           (V)
!
   character YTITLE,blk
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES
!
   CALL AXLOG (0.0,0.0,YTITLE,NY,NUMDVY,DIVLNY,YMIN,90.0,HGT,1,1,    &
   &            NOENDY,0,0)
   CALL AXLOG (XSIZE,0.0,BLK,-1,NUMDVY,DIVLNY,YMIN,90.0,HGT,1,1,     &
   &            NOENDY,1,0)
!
   RETURN
!
end subroutine PLTLOG
SUBROUTINE AX2
!
   USE phys_consts, ONLY: radcn1, radcn2
   IMPLICIT REAL*8           (V)
!
   CHARACTER TITL3*50,HTEN*2
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES
!
   DATA HTEN / '10'/,                                                &
   &     TITL3 / ' RADIANCE WATTS /  (CM**2*STER*CM-1)  * '/,         &
   &     NT3 / 37 /
!
   BB(X,TAVE) = (RADCN1*X**3)/(EXP(RADCN2*X/TAVE)-1.0)
   TEM(X,BBY) = RADCN2*X/ LOG(RADCN1*X**3/BBY+1.0)
!
   ISKIP = 0
   CALL PLOT (XSIZE,0.0,3)
   CALL PLOT (XSIZE,0.0,2)
   CALL PLOT (XSIZE,YSIZE,2)
   CALL PLOT (XSIZE,YSIZE,3)
   X2 = XSIZE+HGT
   X3 = XSIZE+3.*HGT
   YOLD = YSIZE+1.001*HGT
   HGTE = HGT*(2./3.)
   S2 = V2
   BB1 = BB(S2,YMAX)
   IB1 =  LOG10(BB1)
   B1 = 10.**IB1
10 DO 30 I = 1, 5
      FN = 2*(5-I)
      IF (I.EQ.5) FN = 1.
      BC = B1*FN
      TMP = TEM(S2,BC)
      IF (TMP.GT.YMAX) GO TO 30
      YIN = (TMP-YMIN)/DY
      IF (YIN.LT.HGT) GO TO 40
      DELY = YOLD-YIN
      IF (DELY.LE.HGT) GO TO 20
      IF ((ISKIP.EQ.1).AND.(I.NE.5)) GO TO 30
      YOLD = YIN
      CALL SYMBOL (XSIZE,YIN,HGT,'3',0.0,-1)
      IF (I.NE.5) GO TO 30
      YIN = YIN-HGT/2.
      CALL SYMBOL (X2,YIN,HGT,HTEN,0.0,2)
      HIP = IB1
      YINE = YIN+HGTE
      CALL NUMBER (X3,YINE,HGTE,HIP,0.0,-1)
      GO TO 30
20    ISKIP = 1
30 END DO
   B1 = B1/10.
   IB1 = IB1-1
   GO TO 10
40 YLT = NT3*HGT*.333
   YHV = YSIZE/2.
   YIN = YHV-YLT
   X4 = X3+4.5*HGT
   CALL SYMBOL (X4,YIN,HGT,TITL3,90.0,NT3)
!
   RETURN
!
end subroutine AX2
SUBROUTINE FPLINE (NLO,NPTS)
!
   IMPLICIT REAL*8           (V)
!
   COMMON XX(2450),YY(2450)
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELX,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
   DIMENSION X(2450)
!
   EQUIVALENCE (X(1),XX(1))
!
!     FOUR POINT LAGRANGE INTERPOLATION:
!
   IF (NLO.EQ.4) GO TO 10
!
!     INITIALIZE FIRST PANEL:
!
   YY(1) = YY(2)
   DV1 = DV/4.
   DV2 = 2.*DV1
   DV3 = 3.*DV1
   X00 = -7./128.
   X01 = 105./128.
   X02 = 35./128.
   X03 = -5./128.
   X10 = -1./16.
   X11 = 9./16.
   IF (NPTS.LT.4) GO TO 50
!
10 XXI = XX(NLO)
   ILO = 2
20 IHI = ILO+600
   IF (IHI.GE.NPTS) IHI = NPTS-1
   IMAX = IHI-1
   IP0 = 1
   DO 30 I = ILO, IMAX
      XI = XXI+DV* REAL(I-NLO)
      X(IP0) = XI
      X(IP0+1) = XI+DV1
      X(IP0+2) = XI+DV2
      X(IP0+3) = XI+DV3
      Y(IP0) = YY(I)
      Y(IP0+1) = X00*YY(I-1)+X01*YY(I)+X02*YY(I+1)+X03*YY(I+2)
      Y(IP0+2) = X10*(YY(I-1)+YY(I+2))+X11*(YY(I)+YY(I+1))
      Y(IP0+3) = X03*YY(I-1)+X02*YY(I)+X01*YY(I+1)+X00*YY(I+2)
      IP0 = IP0+4
30 END DO
   X(IP0) = XI+DV
   Y(IP0) = YY(IHI)
   J1 = 1
   CALL FSCLIN (J1,IP0,X,Y)
   IF (IHI+1.EQ.NPTS) GO TO 40
   ILO = IHI
   GO TO 20
!
!     PRESERVE LAST THREE POINTS OF PANEL FOR SUBSEQUENT PANEL:
!
40 YY(1) = YY(IMAX)
   YY(2) = YY(IHI)
   YY(3) = YY(NPTS)
!
   RETURN
!
!     CASE IN WHICH INITIAL PANEL HAS LESS THAN THREE POINTS:
!
50 YY(3) = YY(NPTS)
   YY(2) = YY(NPTS-1)
   YY(1) = YY(2)
!
   RETURN
!
end subroutine FPLINE
SUBROUTINE LINT (NST,NND,J)
!
   IMPLICIT REAL*8           (V)
!
   COMMON XX(2450),YY(2450)
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
!
   DO 10 I = NST, NND
      YY(J) = Y(I)
      J = J+1
10 END DO
!
   RETURN
!
end subroutine LINT
SUBROUTINE EXPT (NST,NND,J)
!
   IMPLICIT REAL*8           (V)
!
   COMMON XX(2450),YY(2450)
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
!
   DO 10 I = NST, NND
      YY(J) = 0.
      IF (Y(I).lt.20.) YY(J) = EXP(-Y(I))
      J = J+1
10 END DO
!
   RETURN
!
end subroutine EXPT
SUBROUTINE XNTLOG (NST,NND,J)
!
   IMPLICIT REAL*8           (V)
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON XX(2450),YY(2450)
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
!
   DO 10 I = NST, NND
      IF (Y(I).LE.0.) Y(I) = EXPMIN
      YY(J) = - LOG(Y(I))
      J = J+1
10 END DO
!
   RETURN
!
end subroutine XNTLOG
SUBROUTINE TEMPFN (NST,NND,J)
!
   USE phys_consts, ONLY: radcn1, radcn2
   IMPLICIT REAL*8           (V)
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON XX(2450),YY(2450)
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
!
   DO 10 I = NST, NND
      YY(J) = 0.
      IF (Y(I).GE.EXPMIN) THEN
         X = RADCN1*(XX(J)**3)/Y(I)
         IF (X.GE.EXPMIN) THEN
            Z = LOG(X+1.)
            YY(J) = RADCN2*(XX(J)/Z)
         ENDIF
      ENDIF
      J = J+1
10 END DO
!
   RETURN
!
end subroutine TEMPFN
SUBROUTINE TENLOG (NST,NND,J)
!
   IMPLICIT REAL*8           (V)
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON XX(2450),YY(2450)
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
!
   DO 10 I = NST, NND
      IF (Y(I).LE.0.) Y(I) = EXPMIN
      YY(J) = LOG10(Y(I))
      J = J+1
10 END DO
!
   RETURN
!
end subroutine TENLOG
SUBROUTINE DBOD (NST,NND,J)
!
   IMPLICIT REAL*8           (V)
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON XX(2450),YY(2450)
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
!
   DATA CON / 4.3429448 /
!
   DO 10 I = NST, NND
      YY(J) = Y(I)*CON
      J = J+1
10 END DO
!
   RETURN
!
end subroutine DBOD
SUBROUTINE DBTR (NST,NND,J)
!
   IMPLICIT REAL*8           (V)
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON XX(2450),YY(2450)
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
!
   DATA CON / 4.3429448 /
!
   DO 10 I = NST, NND
      IF (Y(I).LE.0.) Y(I) = EXPMIN
      YY(J) = -CON* LOG(Y(I))
      J = J+1
10 END DO
!
   RETURN
!
end subroutine DBTR
SUBROUTINE LDBOPD (NST,NND,J)
!
   IMPLICIT REAL*8           (V)
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON XX(2450),YY(2450)
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
!
   DATA CON / 4.3429448 /
!
   DO 10 I = NST, NND
      IF (Y(I).LE.0.) Y(I) = EXPMIN
      YY(J) = LOG10(Y(I)*CON)
      J = J+1
10 END DO
!
   RETURN
!
end subroutine LDBOPD
SUBROUTINE LDBTR (NST,NND,J)
!
   IMPLICIT REAL*8           (V)
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON XX(2450),YY(2450)
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
!
   DATA CON / 4.3429448 /
!
   DO 10 I = NST, NND
      IF (Y(I).LE.0.) Y(I) = EXPMIN
      YY(J) = -CON* LOG(Y(I))
      IF (YY(J).LE.0.) YY(J) = EXPMIN
      YY(J) = LOG10(YY(J))
      J = J+1
10 END DO
!
   RETURN
!
end subroutine LDBTR
SUBROUTINE LINSTC (NST,NND,J,XLOGCN)
!
   IMPLICIT REAL*8           (V)
!
   COMMON XX(2450),YY(2450)
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
!
   DO 10 I = NST, NND
      YY(J) = Y(I)*XLOGCN
      J = J+1
10 END DO
!
   RETURN
!
end subroutine LINSTC
SUBROUTINE XLOGLN (NST,NND,J)
!
   IMPLICIT REAL*8           (V)
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON XX(2450),YY(2450)
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
!
   DO 10 I = NST, NND
      CAY = Y(I)
      IF (CAY.EQ.0) CAY = EXPMIN
      CAY = - LOG(CAY)
      IF (CAY.EQ.0) CAY = EXPMIN
      YY(J) = LOG10(CAY)
      J = J+1
10 END DO
!
   RETURN
!
end subroutine XLOGLN
SUBROUTINE MNMX (NST,NND,JMIN,JMAX)
!
   IMPLICIT REAL*8           (V)
!
   COMMON XX(2450),YY(2450)
   COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)
   J = JMIN
   DO 10 I = NST, NND
      IF (YY(J).LT.YY(JMIN)) JMIN = J
      IF (YY(J).GT.YY(JMAX)) JMAX = J
      J = J+1
10 END DO
   RETURN
end subroutine MNMX
SUBROUTINE FSCLIN (NLO,NPTS,XX,YY)
!
   IMPLICIT REAL*8           (V)
!
   DIMENSION XX(*),YY(*),PNLHDR(2)
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /JPLTFL/ V1P,V2P,DVP,NLIM
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELX,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES
!
   EQUIVALENCE (PNLHDR(1),V1P) , (JPLT,NUMSBX) , (LPLT,NOENDX)
   save dvp_save
!
   YYSTOR = YY(NPTS)
!

   IF (JOUT.GE.1) THEN
      V1P = XX(NLO)
      V2P = XX(NPTS)
      if (npts .eq. nlo) then
         dvp = dvp_save
      else
         DVP = (V2P-V1P)/ REAL(NPTS-NLO)
         dvp_save = dvp
      endif
      NLOW = NLO-1
      NLIM = NPTS-NLOW
      IF (JOUT.LE.2) THEN
         CALL BUFOUT (JPLTFL,PNLHDR(1),NPHDRF)
         CALL BUFOUT (JPLTFL,YY(NLO),NPTS)
      ELSE
         DO 10 II = NLO, NPTS
            VOUT = V1P+ REAL(II-NLO)*DVP
            WRITE (JPLTFL,900) VOUT,YY(II)
10       CONTINUE
      ENDIF
      IF (JOUT.EQ.1.OR.JOUT.EQ.3) THEN
         XX(1) = XX(NPTS)
         YY(1) = YYSTOR
         RETURN
      ENDIF
   ENDIF
!
!   TRUNCATE POINTS OUTSIDE PLOTTING LIMITS:
!
   DO 20 J = 1, NPTS
      YY(J) = MAX(YY(J),YMIN)
      YY(J) = MIN(YY(J),YMAX)
20 END DO
!
   XX(NPTS+1) = V1
   YY(NPTS+1) = YMIN
   XX(NPTS+2) = DX
   YY(NPTS+2) = DY
   IF (NOAXES.EQ.1) THEN
      JCAL = JPLT
      LCAL = LPLT
   ELSE
      JCAL = 0
      LCAL = 0
   ENDIF
!
   CALL LINE (XX,YY,NPTS,1,JCAL,LCAL)
!
!     USE ABOVE CALL TO LINE FOR STANDARD CALCOMP SYSTEM.
!     CALL BELOW MAKES USE OF MODIFIED CALCOMP ROUTINE
!     AT GEOPHYSICS DIRECTORATE OF THE PHILLIPS LABORATIORY.
!
!>    CALL LINE(XX,YY,NPTS,1,JCAL,LCAL,V1,DX,YMIN,DY,0.08)
!
   XX(1) = XX(NPTS)
   YY(1) = YYSTOR
!
   RETURN
!
900 FORMAT (3X,F15.8,4X,1p,G19.8E3)
!
end subroutine FSCLIN
SUBROUTINE AXISL (X,Y,BCD,N,NUMDIV,DIVLEN,NUMSUB,BEGNUM,DELNUM,   &
&                  NUMDEC,THETA,HEIGHT,NRPT,NTURN,NOEND,LSUPR,     &
&                  LTURN)
!
!   MODIFIED VERSION OF AXISL
!
!   WRITTEN BY RICHARD L. TAYLOR   RADC/ET   EEC   NOVEMBER 1980
!
!   ROUTINE TO PLOT A LABELLED LINEAR AXIS
!
!
!   X AND Y ARE THE STARTING COORDINATES OF THE AXIS RELATIVE TO THE
!
!        CURRENT ORIGIN
!
!   BCD IS THE LABEL OF THE AXIS EXPRESSED AS A HOLLERITH CONSTANT
!
!   N IS THE NUMBER OF CHARACTERS IN THE LABEL
!
!        NEGATIVE N PLACES THE LABEL ON THE CLOCKWISE SIDE OF THE AXIS
!
!        POSITIVE N PLACES THE LABEL ON THE COUNTERCLOCKWISE SIDE
!
!   NUMDIV IS THE NUMBER OF MAJOR DIVISIONS
!
!   DIVLEN IS THE LENGTH IN INCHES OF A MAJOR DIVISION
!
!   NUMSUB IS THE NUMBER OF MINOR DIVISIONS PER MAJOR DIVISION
!
!        1 GIVES NO SUBDIVISION TICS, 2 GIVES ONE SUBDIVISION TIC, ETC.
!
!   BEGNUM IS THE NUMBER FOR THE BEGINNING OF THE AXIS
!
!   DELNUM IS THE DELTA NUMBER FOR A MAJOR DIVISION
!
!   NUMDEC IS THE NUMBER OF DECIMAL PLACES DESIRED
!
!        NUMDEC EQUAL TO -1 SUPPRESSES THE DECIMAL POINT
!
!   THETA IS THE ANGLE OF THE AXIS IN DEGREES  (0.0 FOR X, 90.0 FOR Y)
!
!   HEIGHT IS THE HEIGHT OF THE NUMBERS IN INCHES
!
!   NRPT IS THE REPEAT FACTOR FOR THE SCALE NUMBERS (USUALLY INTEGER 1)
!
!        WHEN NRPT IS ZERO THE SCALE NUMBERS WILL BE SUPPRESSED;
!
!        WHEN NRPT = 2, EVERY 2ND SCALE NUMBER WILL BE PRODUCED; ETC.
!
!   NTURN EQUAL TO 1 TURNS THE AXIS NUMBERS BY 90 DEGREES CLOCKWISE,
!        -1 TURNS NUMBERS BY 90 DEGREES COUNTERCLOCKWISE, 0 FOR NO TURN
!
!   NOEND EQUAL TO 1 SUPPRESSES THE NUMBERS AT EITHER END OF THE AXIS,
!
!        2 SUPPRESSES THE BEGINNING NUMBER, 3 THE ENDING NUMBER
!
!   LSUPR EQUAL TO 1 SUPPRESSES THE LABEL
!    LTURN NOT USED
!
   USE phys_consts, ONLY: pi
   COMMON /TITLOC/ XPOS,YPOS
   character BCD(*)
!
   THETA1 = THETA-90.*NTURN
   ANGLE = (PI/180.)*THETA
   SINANG = SIN(ANGLE)
   COSANG = COS(ANGLE)
   SIGNAX =  REAL(ISIGN(1,N))
   SIZMAJ = 0.25*HEIGHT+0.05
   OFFST = HEIGHT*1.5
   DXMAJ = -SIZMAJ*SINANG*SIGNAX
   DYMAJ = SIZMAJ*COSANG*SIGNAX
   DXMIN = 0.5*DXMAJ
   DYMIN = 0.5*DYMAJ
   NSUB = NUMSUB
   IF (NUMSUB.LT.1) NSUB = 1
   SUBDIV = DIVLEN/ REAL(NSUB)
   BCDSIZ = 1.25*HEIGHT
   YBIAS = (-0.50+SIGN(1.25,SIGNAX))*HEIGHT+DYMAJ
   NABS = IABS(N)
   BCDLEN = ( REAL(NABS)-0.4)*BCDSIZ
   S = DIVLEN* REAL(NUMDIV)
   DIVCOS = DIVLEN*COSANG
   DIVSIN = DIVLEN*SINANG
   SIZMAX = HEIGHT
!
!   DRAW DIVISION NUMBERS
!
   NDIV = 0
10 DIGITS = 0.0
   XTIC = X+DIVCOS* REAL(NDIV)
   YTIC = Y+DIVSIN* REAL(NDIV)
   IF (NRPT.EQ.0) GO TO 30
   NSUPR = NDIV-(NDIV/NRPT)*NRPT
   IF (NSUPR.NE.0) GO TO 30
   IF ((NOEND.EQ.1.OR.NOEND.EQ.2).AND.NDIV.EQ.0) GO TO 30
   IF ((NOEND.EQ.1.OR.NOEND.EQ.3).AND.NDIV.EQ.NUMDIV) GO TO 30
   DIVNUM = BEGNUM+DELNUM* REAL(NDIV)
   IF (ABS(DIVNUM).GE.10.0) DIGITS =  LOG10(ABS(DIVNUM))
   DIGITS = AINT(DIGITS+1.0E-12)
   IF (DIVNUM.LT.0.0) DIGITS = DIGITS+1.0
   SIZNUM = (DIGITS+ REAL(NUMDEC)+1.7)*HEIGHT
   XBIAS = -0.5*SIZNUM
   XBIAS1 = 0.
   YBIAS1 = 0.
   IF (NTURN.EQ.0) GO TO 20
   YBIAS1 = YBIAS-SIZNUM-OFFST
   IF (N.LT.0) YBIAS1 = YBIAS+OFFST
   XBIAS1 = XBIAS+HEIGHT*0.5
20 XPOS = XTIC-YBIAS*SINANG+XBIAS*COSANG+YBIAS1*SINANG-XBIAS1*COSANG
   YPOS = YTIC+YBIAS*COSANG+XBIAS*SINANG-XBIAS1*SINANG-YBIAS1*COSANG
   CALL NUMBER (XPOS,YPOS,HEIGHT,DIVNUM,THETA1,NUMDEC)
   SIZMAX = MAX(SIZMAX,SIZNUM)
30 CONTINUE
!
!   DRAW TIC MARKS
!
   CALL PLOT (XTIC,YTIC,3)
   CALL PLOT (XTIC+DXMAJ,YTIC+DYMAJ,2)
   IF (NDIV.EQ.NUMDIV) GO TO 60
   IF (NUMSUB.LE.1) GO TO 50
   DO 40 J = 2, NUMSUB
      SUBLEN = SUBDIV* REAL(J-1)
      XSTIC = XTIC+SUBLEN*COSANG
      YSTIC = YTIC+SUBLEN*SINANG
      CALL PLOT (XSTIC+DXMIN,YSTIC+DYMIN,3)
      CALL PLOT (XSTIC,YSTIC,2)
40 END DO
50 NDIV = NDIV+1
   GO TO 10
!
!   DRAW AXIS
!
60 CALL PLOT (XTIC,YTIC,3)
   CALL PLOT (X,Y,2)
!
!   DRAW LABEL
!
   IF (LSUPR.EQ.1.OR.NABS.EQ.0) RETURN
   XBIAS1 = 0
   YBIAS1 = -SIZNUM-OFFST
   IF (N.LT.0) YBIAS1 = -YBIAS1
   IF (NTURN.EQ.0) YBIAS1 = 0
   XBIAS = 0.5*(S-BCDLEN)
   YBIAS = (-0.50+SIGN(3.25,SIGNAX))*BCDSIZ
   XPOS = X-YBIAS*SINANG+XBIAS*COSANG+YBIAS1*SINANG-XBIAS1*COSANG
   YPOS = Y+YBIAS*COSANG+XBIAS*SINANG-XBIAS1*SINANG-YBIAS1*COSANG
   CALL SYMBOL (XPOS,YPOS,BCDSIZ,BCD,THETA,NABS)
!
   RETURN
!
end subroutine AXISL
SUBROUTINE AXLOG (X,Y,BCD,N,NUMCYC,CYCLEN,BEGEXP,THETA,HEIGHT,    &
&                  NRPT,NTURN,NOEND,LSUPR,LTURN)
!
!
!   WRITTEN BY RICHARD L. TAYLOR   RADC/ET   EEC   NOVEMBER 1980
!
!
!
!   ROUTINE TO PLOT A LABELLED LOGARITHMIC AXIS
!
!
!
!
!   X AND Y ARE THE STARTING COORDINATES OF THE AXIS RELATIVE TO THE
!
!        CURRENT ORIGIN
!
!   BCD IS THE LABEL OF THE AXIS EXPRESSED AS A HOLLERITH CONSTANT
!
!   N IS THE NUMBER OF CHARACTERS IN THE LABEL
!
!        NEGATIVE N PLACES THE LABEL ON THE CLOCKWISE SIDE OF THE AXIS
!
!        POSITIVE N PLACES THE LABEL ON THE COUNTERCLOCKWISE SIDE
!
!   NUMCYC IS THE NUMBER OF CYCLES DESIRED
!
!   CYCLEN IS THE LENGTH OF ONE CYCLE IN INCHES
!
!   BEGEXP IS THE EXPONENT FOR THE BEGINNING OF THE AXIS
!
!   THETA IS THE ANGLE OF THE AXIS IN DEGREES  (0.0 FOR X, 90.0 FOR Y)
!
!   HEIGHT IS THE HEIGHT IN INCHES OF THE TENS
!
!   NRPT IS THE REPEAT FACTOR FOR THE SCALE NUMBERS (USUALLY INTEGER 1)
!
!        WHEN NRPT IS ZERO THE SCALE NUMBERS WILL BE SUPPRESSED;
!
!        WHEN NRPT = 2, EVERY 2ND SCALE NUMBER WILL BE PRODUCED;
!
!        WHEN NRPT = 3, EVERY 3RD SCALE NUMBER WILL BE PRODUCED; ETC.
!
!    NTURN EQUAL TO 1 TURNS THE AXIS NUMBERS BY 90 DEGREES CLOCKWISE,
!                  -1 TURNS NUMBERS BY 90 DEGREES COUNTERCLOCKWISE,
!                   0 FOR NO TURN
!
!   NOEND EQUAL TO 1 SUPPRESSES THE NUMBERS AT EITHER END OF THE AXIS,
!
!        NOEND EQUAL TO 2 SUPPRESSES ONLY THE STARTING NUMBER, AND
!
!        NOEND EQUAL TO 3 SUPPRESSES ONLY THE ENDING NUMBER
!
!   LSUPR EQUAL TO 1 SUPPRESSES THE LABEL
!
!    LTURN EQUAL TO 1 TURNS THE LABEL BY 90 DEGREES CLOCKWISE,
!         -1 TURNS LABEL BY 90 DEGREES COUNTERCLOCKWISE, 0 FOR NO TURN
!
!
   USE phys_consts, ONLY: pi
   COMMON /TITLOC/ XPOS,YPOS
   DIMENSION        SUBDIV(10),DIVLOG(8)
   character BCD(*)
   CHARACTER HTEN*2
!
   DATA DIVLOG / 0.301029995664, 0.477121254720, 0.602059991328,     &
   &              0.698970004336, 0.778151250384, 0.845098040014,     &
   &              0.903089986992, 0.954242509439 /
   DATA HTEN / '10'/
!
   THETA1 = THETA-90.*NTURN
   ANGLE = (PI/180.)*THETA
   SINANG = SIN(ANGLE)
   COSANG = COS(ANGLE)
   SIGNAX =  REAL(ISIGN(1,N))
   SIZMAJ = 0.25*HEIGHT+0.05
   OFFST = HEIGHT*1.5
   DXMAJ = -SIZMAJ*SINANG*SIGNAX
   DYMAJ = SIZMAJ*COSANG*SIGNAX
   DXMIN = 0.5*DXMAJ
   DYMIN = 0.5*DYMAJ
   BCDSIZ = 1.25*HEIGHT
   ENLARG = 1.5
   EXPSIZ = 0.60*HEIGHT*ENLARG
   NABS = IABS(N)
   BCDLEN = ( REAL(NABS)-0.4)*BCDSIZ
   S = CYCLEN* REAL(NUMCYC)
   NUMTIC = 1.
   if (cyclen .lt.1.) NUMTIC = 2- INT(cyclen)
   NUMLOG = 8/NUMTIC
   XBIAS = 1.85*HEIGHT*ENLARG
   YBIAS = 0.70*HEIGHT*ENLARG
   EXPBX = -YBIAS*SINANG+XBIAS*COSANG
   EXPBY = YBIAS*COSANG+XBIAS*SINANG
   IF (NTURN.EQ.0) GO TO 10
   EXPBX = XBIAS*SINANG+YBIAS*COSANG
   EXPBY = YBIAS*SINANG-XBIAS*COSANG
10 DO 20 I = 2, 9
      SUBDIV(I) = DIVLOG(I-1)*CYCLEN
20 END DO
   NNUMB = NUMCYC+1
   SIZMAX = EXPSIZ
   EXP = BEGEXP
   DO 30 I = 1, NNUMB
      DIGITS = 0.0
      IF (ABS(EXP).GE.10.0) DIGITS = LOG10(ABS(EXP))
      DIGITS = AINT(DIGITS+1.0E-12)+0.7
      IF (EXP.LT.0.0) DIGITS = DIGITS+1.0
      SIZNUM = DIGITS*EXPSIZ
      SIZMAX = MAX(SIZMAX,SIZNUM)
      EXP = EXP+1.0
30 END DO
   SIZNUM = SIZMAX+2.0*HEIGHT
!
!   DRAW CYCLE NUMBERS AND EXPONENTS
!
   NCYCLE = 0
   EXP = BEGEXP
   XBIAS = -0.7*HEIGHT
   YBIAS = (-0.75+SIGN(2.00,SIGNAX))*HEIGHT
   XBIAS1 = 0.
   YBIAS1 = 0.
   IF (NTURN.EQ.0) GO TO 40
   XBIAS1 = XBIAS+HEIGHT*0.5
   YBIAS1 = YBIAS1-OFFST-SIZNUM
   IF (N.LT.0) YBIAS1 = YBIAS+OFFST
40 TENBX = -YBIAS*SINANG+XBIAS*COSANG+YBIAS1*SINANG-XBIAS1*COSANG
   TENBY = YBIAS*COSANG+XBIAS*SINANG-YBIAS1*COSANG-XBIAS1*SINANG
50 XTIC = X+ REAL(NCYCLE)*CYCLEN*COSANG
   YTIC = Y+ REAL(NCYCLE)*CYCLEN*SINANG
   IF (NRPT.EQ.0) GO TO 60
   NSUPR = NCYCLE-(NCYCLE/NRPT)*NRPT
   IF (NSUPR.NE.0) GO TO 60
   IF ((NOEND.EQ.1.OR.NOEND.EQ.2).AND.NCYCLE.EQ.0) GO TO 60
   IF ((NOEND.EQ.1.OR.NOEND.EQ.3).AND.NCYCLE.EQ.NUMCYC) GO TO 60
   XPOS = XTIC+TENBX
   YPOS = YTIC+TENBY
   CALL SYMBOL (XPOS,YPOS,HEIGHT,HTEN,THETA1,2)
   CALL NUMBER ((XPOS+EXPBX),(YPOS+EXPBY),EXPSIZ,EXP,THETA1,-1)
60 CONTINUE
!
!   DRAW TIC MARKS
!
   CALL PLOT (XTIC,YTIC,3)
   CALL PLOT (XTIC+DXMAJ,YTIC+DYMAJ,2)
   IF (NCYCLE.EQ.NUMCYC) GO TO 90
   IF (NRPT.LT.0) GO TO 80
   DO 70 ILOG = 1, NUMLOG
      I = ILOG*NUMTIC+1/NUMTIC
      XLOG = XTIC+SUBDIV(I)*COSANG
      YLOG = YTIC+SUBDIV(I)*SINANG
      CALL PLOT (XLOG+DXMIN,YLOG+DYMIN,3)
      CALL PLOT (XLOG,YLOG,2)
70 END DO
80 NCYCLE = NCYCLE+1
   EXP = EXP+1.0
   GO TO 50
!
!   DRAW AXIS
!
90 CALL PLOT (XTIC,YTIC,3)
   CALL PLOT (X,Y,2)
!
!   DRAW LABEL
!
   IF (LSUPR.EQ.1.OR.NABS.EQ.0) RETURN
   XBIAS = 0.5*(S-BCDLEN)
   YBIAS = (-0.50+SIGN(3.25,SIGNAX))*BCDSIZ
   THETA2 = THETA-90.*LTURN
   XBIAS2 = 0.
   OFFST = HEIGHT*2.5
   YBIAS2 = -SIZNUM-OFFST
   IF (N.LT.0) YBIAS2 = OFFST
   IF (LTURN.EQ.0) GO TO 100
   XBIAS2 = XBIAS-0.5*(S-HEIGHT)
100 XPOS = X-YBIAS*SINANG+XBIAS*COSANG+YBIAS2*SINANG-XBIAS2*COSANG
   YPOS = Y+YBIAS*COSANG+XBIAS*SINANG-YBIAS2*COSANG-XBIAS2*SINANG
   CALL SYMBOL (XPOS,YPOS,BCDSIZ,BCD,THETA2,NABS)
!
   RETURN
!
end subroutine AXLOG
SUBROUTINE LBH1H2 (YT,XT,YID)
!
!     LABEL H1,H2,ANGLE
!
   DIMENSION YID(10)
!
   CHARACTER*8      YID
!
   CHARACTER CARD3*30,ACRD3*30
   CHARACTER GOUT*64
!
   DATA CARD3 / '     H1,       H2,       ANGL'/
!
   WRITE (GOUT,900) (YID(I),I=3,7)
   READ (GOUT,905) H1,H2,ANGLE
   WRITE (GOUT,910) H1,H2,ANGLE
   READ (GOUT,915) ACRD3
   YT = YT-0.2
   CALL SYMBOL (XT,YT,0.12,CARD3,0.0,30)
   YT = YT-0.2
   CALL SYMBOL (XT,YT,0.12,ACRD3,0.0,30)
!
   RETURN
!
900 FORMAT (8A8)
905 FORMAT (3F6.3)
910 FORMAT (3F10.3)
915 FORMAT (A30)
!
end subroutine LBH1H2
SUBROUTINE LBAERS (YT,XT,YID)
!
!     LABEL AERSOL
!
   DIMENSION YID(10)
!
   CHARACTER*8      YID
!
   CHARACTER CARD2*85,ACRD2*85
   CHARACTER GOUT*85
!
   DATA CARD2/           'IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,   VIS,&
   &      WSS,      WHH,      RAINRT,   GNDALT'/
!
   WRITE (GOUT,900) (YID(I),I=3,7)
   READ (GOUT,905) IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,IVIS,IWSS,    &
   &                IWHH,IRAINR,IGNDAL
   VIS = IVIS/10
   WSS = IWSS/10
   WHH = IWHH/10
   RAINRT = IRAINR/10
   GNDALT = IGNDAL/10
   WRITE (GOUT,910) IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH, &
   &                 RAINRT,GNDALT
   READ (GOUT,915) ACRD2
   YT = YT-0.2
   CALL SYMBOL (XT,YT,0.12,CARD2,0.0,85)
   YT = YT-0.2
   CALL SYMBOL (XT,YT,0.12,ACRD2,0.0,85)
!
   RETURN
!
900 FORMAT (8A8)
905 FORMAT (18X,4I1,I2,I1,I4,3I3,I2)
910 FORMAT (2I5,I6,I7,2I5,5F10.3)
915 FORMAT (A85)
!
end subroutine LBAERS
SUBROUTINE YDIH1 (H1,H2,ANGLE,YID)
!
!     PUT H1 H2 ANGLE INTO YID
!
   DIMENSION YID(10)
!
   CHARACTER*8      YID
!
   CHARACTER GOUT*64,HOUT*40
   WRITE (HOUT,900) (YID(I),I=3,7)
   IF (H1.LT.100.AND.H2.LT.100.) THEN
      WRITE (GOUT(1:18),905) H1,H2,ANGLE
   ELSE
      WRITE (GOUT(1:18),910) H1,H2,ANGLE
   ENDIF
   GOUT(19:40) = HOUT(19:40)
   READ (GOUT,900) (YID(I),I=3,7)
!
   RETURN
!
900 FORMAT (8A8)
905 FORMAT (2F6.3,F6.2)
910 FORMAT (3F6.2)
!
end subroutine YDIH1
