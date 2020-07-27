!     path:      $HeadURL$
!     revision:  $Revision$
!     created:   $Date$
!     presently: %H%  %T%
!
!  --------------------------------------------------------------------------
! |  Copyright ��, Atmospheric and Environmental Research, Inc., 2015         |
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
SUBROUTINE XMERGE (NPTS,LFILE,MFILE,JPATHL)
!
   IMPLICIT REAL*8           (V)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
!     XMERGE CALL ABSMRG,EMINIT,RADMRG
!
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYER,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   COMMON /CVRXMR/ HNAMXMR,HVRXMR
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYHDR,YI1,YID(10),LSTWDF
!
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,      &
   &    pad_3,angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!
   COMMON /XME/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(5000)
!
   EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),              &
   &            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),              &
   &            (FSCDID(5),IEMIT) , (FSCDID(7),IPLOT),                &
   &            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                &
   &            (FSCDID(11),IMRG) , (FSCDID(16),LAYR1),               &
   &            (FSCDID(17),NLAYHD)
!
   CHARACTER*18 HNAMXMR,HVRXMR
!
!     ASSIGN CVS VERSION NUMBER TO MODULE
!
   HVRXMR = '$Revision$'
!
   IOD = 0
!
   IEXFIL = 20
   IAFIL = 14
!
!     --------------------------------------------------------------
!
!     Special case of merging optical depths from multiple files to
!     multiple files, for use when calculating radiance derivatives.
!     The subroutine call tree is as follows:
!
!        LBLSUB called XLAYER when IHIRAC+IATM+IMRG > 0
!        XLAYER called XMERGE when IMRG   = 10
!        XMERGE calls  ABSMRG when IMRG   = 10
!
   IF (IMRG.EQ.10) THEN
      CALL ABSMRG (NPTS,LFILE,MFILE,JPATHL)
      CALL ENDFIL (MFILE)
      RETURN
   ENDIF
!
!     --------------------------------------------------------------
!
!     WHEN IAERSL EQUALS 2 CALL ADARSL TO ADD ABSORPTION AND SCATTERING
!     TO COMMON BLOCKS FOR USE IN A TRANSMITTANCE RUN
!
!     --------------------------------------------------------------
!
!     IEMIT = 0  =>  Optical depth calculation only
!
   IF (IEMIT.EQ.0) THEN
      IF (LAYER.EQ.1) THEN
         CALL ABINIT (NPTS,MFILE,JPATHL)
      ELSE
         WRITE (IPR,900) XID,(YID(M),M=1,2)
         CALL ABSMRG (NPTS,LFILE,MFILE,JPATHL)
      ENDIF
   ELSE
!
!     --------------------------------------------------------------
!
!     IEMIT > 0 TO REACH THIS STATEMENT
!
      WRITE (IPR,900) XID,(YID(M),M=1,2)
      IF (IAERSL.GE.1 .and. iaersl.ne.5) then
         IF (LAYER.EQ.1) REWIND IEXFIL
         CALL GETEXT (IEXFIL,LAYER,IEMIT)
      ENDIF
!
      TBND = 0.
!
!        -----------------------------------------------------------
!
!        IEMIT = 1  =>  Radiance and Transmittance calculated
!
      IF (IEMIT.EQ.1) THEN
         IF (IMRG.eq.36 .or. IMRG.eq.46) THEN
            IF (LAYER.EQ.1) THEN
               TBND = TMPBND
               CALL FLINIT (NPTS,MFILE,JPATHL,TBND, 1)
            ELSE
               CALL FLUXUP (NPTS,LFILE,MFILE,JPATHL,TBND)
            ENDIF
         ELSE
            IF (LAYER.EQ.lh1) THEN
               IF (IPATHL.EQ.1) TBND = TMPBND
               CALL EMINIT (NPTS,MFILE,JPATHL,TBND)
            ELSE
               IF (IPATHL.EQ.3.AND.LAYER.EQ.LH2) TBND = TMPBND
               CALL RADMRG (NPTS,LFILE,MFILE,JPATHL,TBND)
            ENDIF
         ENDIF
      ENDIF
!
!        -----------------------------------------------------------
!
!        IEMIT = 3  =>  Radiance, Transmittance and Radiance
!                       derivative calculated
!                       ipathl=1 = downlooking (upwelling)
!                       ipathl=3 = uplooking (downwelling)
      IF (IEMIT.EQ.3) THEN
         IF (LAYER.EQ.1) THEN
            IF (IPATHL.EQ.1) TBND = TMPBND
            CALL EMADL1 (NPTS,MFILE,JPATHL,TBND)
         ELSE
            IF (IPATHL.EQ.3.AND.LAYER.EQ.LH2) TBND = TMPBND
            CALL EMADMG (NPTS,LFILE,MFILE,JPATHL,TBND)
         ENDIF
      ENDIF
   ENDIF

   RETURN
!
900 FORMAT (///,1X,10A8,2X,2(1X,A8,1X))
!
END SUBROUTINE XMERGE
!
!     ----------------------------------------------------------------
!
SUBROUTINE XMERGI (NPTS,LFILE,MFILE,JPATHL)
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
      MXMOL, MXTRAC, MX_XS
   IMPLICIT REAL*8           (V)
!
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
!     XMERGI CALL ABINIT,ABSINT,ABSOUT
!
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYER,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID

   COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),          &
   &                MSPNL1(MXLAY),                                    &
   &                MSLAY1,ISFILE,JSFILE,KSFILE,LSFILE,MSFILE,IEFILE, &
   &                JEFILE,KEFILE
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYHDR,YI1,YID(10),LSTWDF
!
   COMMON /XMI/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(4819)
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,      &
   &    pad_3,angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!
   EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),              &
   &            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),              &
   &            (FSCDID(5),IEMIT) , (FSCDID(7),IPLOT),                &
   &            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                &
   &            (FSCDID(11),IMRG) , (FSCDID(16),LAYR1),               &
   &            (FSCDID(17),NLAYHD)
!
   IF (IEMIT.EQ.1) GO TO 10
   IF (LAYER.EQ.LH1.AND.IANT.NE.-1) THEN
      CALL ABINIT (NPTS,MFILE,JPATHL)
   ELSE
!
      WRITE (IPR,900) XID,(YID(M),M=1,2)
      CALL ABSINT (NPTS,LFILE,MFILE,JPATHL)
   ENDIF
!
   GO TO 20
!
!     IEMIT = 1 TO REACH THIS STATEMENT
!
10 CONTINUE
   WRITE (IPR,900) XID,(YID(M),M=1,2)
   IF (IAERSL.GE.1 .and. iaersl.ne.5) then
      IF (LAYER.EQ.1) REWIND IEXFIL
      CALL GETEXT (IEXFIL,LAYER,IEMIT)
   ENDIF
!
   TBND = 0.
!
   IF (IMRG.NE.35.AND.IMRG.NE.45) THEN
      IF (LAYER.EQ.LH1.AND.IANT.NE.-1) THEN
         IF (JPATHL.EQ.1.AND.LAYER.EQ.1) TBND = TMPBND
         CALL EMINIT (NPTS,MFILE,JPATHL,TBND)
      ELSE
         IF (JPATHL.EQ.3.AND.LAYER.EQ.LH2) TBND = TMPBND
         CALL RADINT (NPTS,LFILE,MFILE,JPATHL,TBND)
      ENDIF
   ELSE
      IF (LAYER.EQ.LH1.AND.IANT.NE.-1) THEN
         TBND = TMPBND
         CALL FLINIT (NPTS,MFILE,JPATHL,TBND,-1)
      elseif (layer.eq.1) then
         CALL FLUXDN (NPTS,LFILE,MFILE,JPATHL,TBND,1)
      else
         CALL FLUXDN (NPTS,LFILE,MFILE,JPATHL,TBND,-1)
      ENDIF
   ENDIF
!
20 CONTINUE
!
   RETURN
!
900 FORMAT (///,1X,10A8,2X,2(1X,A8,1X))
!
END SUBROUTINE XMERGI
!
!     ----------------------------------------------------------------
!
SUBROUTINE ABINIT (NPTS,MFILE,JPATHL)
!
   IMPLICIT REAL*8           (V)
!
   COMMON ODLAY(-2:2407)
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /ABSHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /ABSPNL/ V1P,V2P,DVP,NLIM,NSHFT,NPNTS
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
   DIMENSION XFILHD(2),PNLHDR(2)
   DIMENSION ODLAYR(2)
!
   CHARACTER*40 CEXT,CYID
!
   EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P)
   EQUIVALENCE (ODLAY(1),ODLAYR(1)) , (FSCDID(4),IAERSL),            &
   &            (FSCDID(5),IEMIT) , (FSCDID(7),IPLOT),                &
   &            (FSCDID(8),IPATHL) , (FSCDID(16),LAYR1)
!
!
!     ***********************************************************
!     ****** THIS SUBROUTINE INITALIZES MERGE FOR OPTICAL  ******
!     ****** DEPTHS                                        ******
!     ***********************************************************
!
   CALL CPUTIM (TIME)
   WRITE (IPR,900) TIME
   NPANLS = 0
!
   CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)
   IF (JPATHL.GE.1) IPATHL = JPATHL
   PLAY = PAVE
   TLAY = TAVE
!
!     FOR AEROSOL RUNS, MOVE EXTID INTO YID
!
   IF (iaersl.ge.1 .and. iaersl.ne.5) THEN
      WRITE (CEXT,'(10A4)') EXTID
      WRITE (CYID,'(5A8)') (YID(I),I=3,7)
      CYID(19:40) = CEXT(19:40)
      READ (CYID,'(5A8)') (YID(I),I=3,7)
   ENDIF
!
   IEMIT = 0
   FACT = 1.
   IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.
   DO 10 MOL = 1, NMOL
      WK(MOL) = WK(MOL)*FACT
10 END DO
   WBROAD = WBROAD*FACT
   LAYR1 = LAYER
   WRITE (IPR,905) LAYR1,KFILE,MFILE
!
   CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)
   DVXM = DV
!
20 CONTINUE
!
   CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
   IF (KEOF.LE.0) GO TO 40
   CALL BUFIN (KFILE,KEOF,ODLAYR(1),NLIM)
!
   IF (FACT.EQ.2.) THEN
      DO 30 I = 1, NLIM
         ODLAYR(I) = ODLAYR(I)+ODLAYR(I)
30    CONTINUE
   ENDIF
!
   CALL ABSOUT (V1P,V2P,DVP,NLIM,1,MFILE,NPTS,ODLAYR,NPANLS)
   GO TO 20
!
40 CONTINUE
!
   CALL CPUTIM (TIME1)
   TIM = TIME1-TIME
   WRITE (IPR,910) TIME1,TIM
!
   RETURN
!
900 FORMAT ('0 THE TIME AT THE START OF ABINIT IS ',F12.3)
905 FORMAT ('0 INITIAL LAYER',I5,/,'0 FILE ',I5,                      &
   &        ' INITIALIZED ONTO FILE',I5)
910 FORMAT ('0 THE TIME AT THE END OF ABINIT IS ',F12.3/F12.3,        &
   &        ' SECS WERE REQUIRED FOR THIS MERGE ')
!
END SUBROUTINE ABINIT
!
!     ---------------------------------------------------------------
!
SUBROUTINE OPNMRG(LFILE,PATH1,L1,FM1,PATH2,L2,FM2,MFILE,          &
&                  PATH3,FM3)
!
!     This subroutine opens file for merging in ABSMRG, when IMRG = 10
!
   LOGICAL OP
   CHARACTER*57 FILE1,FILE2,FILE3
   CHARACTER*55 PATH1,PATH2,PATH3
   CHARACTER*11 CFORM
   CHARACTER*10 FM1,FM2,FM3
!
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID

   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
!           123456789-123456789-123456789-123456789-123456789-1234567
   DATA FILE1 /                                                      &
   &     '                                                         '/,&
   &     FILE2 /                                                      &
   &     '                                                         '/,&
   &     FILE3 /                                                      &
   &     '                                                         '/
   DATA CFORM / 'UNFORMATTED' /
!
   INQUIRE (UNIT=LFILE,OPENED=OP)
   IF (OP) CLOSE (LFILE)
   WRITE(FILE1,FM1) PATH1,L1
   OPEN(UNIT=LFILE,FILE=FILE1,FORM=CFORM,STATUS='OLD')
!
   IF (L2.NE.L1) THEN
      INQUIRE (UNIT=KFILE,OPENED=OP)
      IF (OP) CLOSE (KFILE)
      WRITE(FILE2,FM2) PATH2,L2
      OPEN(UNIT=KFILE,FILE=FILE2,FORM=CFORM,STATUS='OLD')
   ELSE
      FILE2 = 'NO FILE USED'
   ENDIF
!
   INQUIRE (UNIT=MFILE,OPENED=OP)
   IF (OP) CLOSE (MFILE)
   WRITE(FILE3,FM3) PATH3,L2
   OPEN(UNIT=MFILE,FILE=FILE3,FORM=CFORM,STATUS='UNKNOWN')
!
!     Write procedure
!
   WRITE(IPR,900) FILE1,FILE2,FILE3
!
   RETURN
!
900 FORMAT ('     Merged file:  ',A57,/,                              &
   &        '       With file:  ',A57,/,                              &
   &        '  To obtain file:  ',A57,/)
!
END SUBROUTINE OPNMRG
!
!     ----------------------------------------------------------------
!
SUBROUTINE ABSMRG (NPTS,LFILE,MFILE,JPATHL)
!
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE ABSMRG MERGES ABSORPTION VALUES FROM KFILE INTO LFILE
!
   COMMON R1(2410),OLDR1(2410)
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /ABSHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
   COMMON /ABSPNL/ V1P,V2P,DVP,NLIM,NSHFT,NPNTS
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   DIMENSION SAVOR1(5),A1(10),A2(10),A3(10),A4(10)
   DIMENSION WKSAV(35)
   DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)
!
   CHARACTER*40 CYID
!
   EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                 &
   &            (OPNLHD(1),V1PO) , (FSCDID(4),IAERSL),                &
   &            (FSCDID(5),IEMIT) , (FSCDID(8),IPATHL),               &
   &            (FSCDID(16),LAYR1)
!
   CALL CPUTIM (TIME)
   NPANLS = 0
   IF (NOPR.EQ.0) WRITE (IPR,900) TIME
   CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)
   DVL = DV
   LAY1SV = LAYR1
   PL = PAVE
   TL = TAVE
   WTOTL = 0.
   DO 10 MOL = 1, NMOL
      WTOTL = WTOTL+WK(MOL)
      WKSAV(MOL) = WK(MOL)
10 END DO
   WTOTL = WTOTL+WBROAD
   WN2SAV = WBROAD
!
!     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)
!
   IF (iaersl.ge.1 .and. iaersl.ne.5)                                &
   &                 WRITE (CYID,'(5A8)') (YID(I),I=3,7)
   CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)
   IF (iaersl.ge.1 .and. iaersl.ne.5)                                &
   &                 READ (CYID,'(5A8)') (YID(I),I=3,7)
   IF (JPATHL.GE.1) IPATHL = JPATHL
   PLAY = PAVE
   TLAY = TAVE
   DVK = DV
   LAYR1 = LAY1SV
   FACT = 1.
   IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.
   IF (DVL.EQ.DVK) ITYPE = 0
   IF (DVL.GT.DVK) ITYPE = DVK/(DVL-DVK)+0.5
   IF (DVL.LT.DVK) ITYPE = -INT(DVL/(DVK-DVL)+0.5)
   SAVOR1(4) = 0.
!
!     ITYPE .LT. 0  IF DV(K-1) IS LESS THAN DV(K)
!
   IF (ITYPE.LT.0) STOP ' ABSMRG: ITYPE LT 0 '
   ITYPE = IABS(ITYPE)
   WTOTK = 0.
   DO 20 MOL = 1, NMOL
      WTOTK = WTOTK+FACT*WK(MOL)
      WK(MOL) = FACT*WK(MOL)+WKSAV(MOL)
20 END DO
   WTOTK = WTOTK+FACT*WBROAD
   PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)
   TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)
   WBROAD = FACT*WBROAD+WN2SAV
   SECANT = 0.
   IF (NOPR.EQ.0) WRITE (IPR,905) LAYR1,LAYER,KFILE,LFILE,MFILE
   IEMIT = 0
!
!     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES
!
   CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)
   DVXM = DV
   DO 30 K = 1, 5
      SAVOR1(K) = 0.
30 END DO
   ATYPE = ITYPE
   AP = 1.0/(ATYPE+1.0)
   IF (ITYPE.NE.0) GO TO 80
!
!     1/1 RATIO ONLY
!
40 CONTINUE
   CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)
   IF (LEOF.LE.0) GO TO 50
   CALL BUFIN (LFILE,LEOF,OLDR1(1),NLIMO)
50 CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
   IF (KEOF.LE.0) GO TO 250
   CALL BUFIN (KFILE,KEOF,R1(1),NLIM)
!
   IF (FACT.EQ.1) THEN
      DO 60 KOD = 1, NLIM
         R1(KOD) = R1(KOD)+OLDR1(KOD)
60    CONTINUE
!
   ELSE
      DO 70 KOD = 1, NLIM
         R1(KOD) = R1(KOD)+R1(KOD)+OLDR1(KOD)
70    CONTINUE
   ENDIF
!
   CALL ABSOUT (V1P,V2P,DVP,NLIM,1,MFILE,NPTS,R1,NPANLS)
!
   GO TO 40
!
!     ALL RATIOS EXCEPT 1/1
!
80 LL = ITYPE+1
   DO 90 JPG = 1, ITYPE
      APG = JPG
      P = 1.0-(AP*APG)
!
!    THE FOLLOWING ARE THE CONSTANTS FOR THE LAGRANGE 4 POINT
!    INTERPOLATION.
!
      A1(JPG) = -P*(P-1.0)*(P-2.0)/6.0
      A2(JPG) = (P**2-1.0)*(P-2.0)*0.5
      A3(JPG) = -P*(P+1.0)*(P-2.0)*0.5
      A4(JPG) = P*(P**2-1.0)/6.0
90 END DO
!
!     ********  BEGINNING OF LOOP THAT DOES INTERPOLATION  *********
!
   CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)
   CALL BUFIN (LFILE,LEOF,OLDR1(1),NLIMO)
   MAXLF = NLIMO
   NVS = 1
   CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
   CALL BUFIN (KFILE,KEOF,R1(1),NLIM)
   IF (KEOF.LE.0) GO TO 250
   II = 1
   DIF = DVP*0.01
   IF (ABS(V1PO-V1P).LT.DIF) GO TO 120
!
!     V1P  <  V1PO  LASER OPTION
!
100 NVS = NVS+1
   V1PN = V1PO+DVPO*(NVS-1)
   IF (ABS(V1PN-V1P).LT.DIF) GO TO 120
   IF (V1PN.LT.V1P) GO TO 100
!
110 II = II+1
   V1PP = V1P+DVP*(II-1)
   IF (ABS(V1PN-V1PP).LT.DIF) GO TO 120
   IF (V1PP.LT.V1PN) GO TO 110
!
   GO TO 100
120 R1(II) = FACT*R1(II)+OLDR1(NVS)
   V1PN = V1PO+DVPO*(NVS-1)
   V1PP = V1P+DVP*(II-1)
   II = II+1
130 JJ = 1
!
   DO 240 JPG = 1, LL
      IF (JPG.EQ.LL) GO TO 140
      IF (NVS.EQ.1) GO TO 150
      GO TO 170
140   IF (FACT.EQ.1.) THEN
         R1(II) = R1(II)+OLDR1(NVS)
      ELSE
         R1(II) = R1(II)+R1(II)+OLDR1(NVS)
      ENDIF
      V1PN = V1PO+DVPO*(NVS-1)
      V1PP = V1P+DVP*(II-1)
      GO TO 190
150   IF (SAVOR1(4).EQ.0.0) GO TO 160
      OLDR1Y = SAVOR1(4)
      GO TO 180
160   OLDR1Y = OLDR1(1)
      GO TO 180
170   OLDR1Y = OLDR1(NVS-1)
180   OLDR1I = A1(JJ)*OLDR1Y+A2(JJ)*OLDR1(NVS)+A3(JJ)*OLDR1(NVS+1)+  &
         A4(JJ)*OLDR1(NVS+2)
      IF (FACT.EQ.1.) THEN
         R1(II) = R1(II)+OLDR1I
      ELSE
         R1(II) = R1(II)+R1(II)+OLDR1I
      ENDIF
190   NVS = NVS+1
      IF (NVS.LE.MAXLF-2) GO TO 200
      SAVOR1(1) = OLDR1(NVS-1)
      SAVOR1(2) = OLDR1(NVS)
      SAVOR1(3) = OLDR1(NVS+1)
      SAVOR1(4) = OLDR1(NVS-2)
      CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)
      IF (LEOF.LE.0) GO TO 210
      MAXLF = NLIMO+3
      CALL BUFIN (LFILE,LEOF,OLDR1(4),NLIMO)
      OLDR1(1) = SAVOR1(1)
      OLDR1(2) = SAVOR1(2)
      OLDR1(3) = SAVOR1(3)
      NVS = 2
200   II = II+1
      JJ = JJ+1
      IF (II.GT.NLIM) GO TO 230
      GO TO 240
210   II = II+1
      AVRG = (SAVOR1(3)+SAVOR1(2))*0.5
220   R1(II) = FACT*R1(II)+AVRG
      II = II+1
      IF (II.LE.NLIM) GO TO 220
!
!     WRITE OUTPUT FILE
!
230   CALL ABSOUT (V1P,V2P,DVP,NLIM,1,MFILE,NPTS,R1,NPANLS)
!
      CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
      IF (KEOF.LE.0) GO TO 250
      CALL BUFIN (KFILE,KEOF,R1(1),NLIM)
      II = 1
240 END DO
   NVS = NVS-1
   GO TO 130
250 CONTINUE
!
   CALL CPUTIM (TIME1)
   TIM = TIME1-TIME
   IF (NOPR.EQ.0) WRITE (IPR,910) TIME1,TIM
   RETURN
!
900 FORMAT ('0 THE TIME AT THE START OF ABSMRG IS ',F12.3)
905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,   &
   &        ' MERGED WITH FILE ',I5,' ONTO FILE',I5)
910 FORMAT (' THE TIME AT THE END OF ABSMRG IS',F12.3/F12.3,          &
   &        ' SECS. WERE REQUIRED FOR THIS ADDITION')
!
END SUBROUTINE ABSMRG
!
!     ----------------------------------------------------------------
!
SUBROUTINE ABSINT (NPTS,LFILE,MFILE,JPATHL)
!
   IMPLICIT REAL*8           (V)
!
   COMMON NEWOD(2410),ODLAY(-2:2407)
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /ABSHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
   COMMON /ABSPNL/ V1P,V2P,DVP,NLIM,NSHFT,NPNTS
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
   DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)
   DIMENSION A1(0:100),A2(0:100),A3(0:100),A4(0:100)
   DIMENSION OLDOD(2),ODLAYR(2)
   DIMENSION WKSAV(35)
!
   CHARACTER*40 CYID
   REAL NEWOD
!
   EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                 &
   &            (OPNLHD(1),V1PO)
   EQUIVALENCE (NEWOD(1),OLDOD(1)) , (ODLAY(1),ODLAYR(1)),           &
   &            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),               &
   &            (FSCDID(16),LAYR1)
!
!     ***********************************************************
!     ****** THIS SUBROUTINE DOES LAYER MERGE FOR OPTICAL  ******
!     ****** DEPTHS USING FOUR POINT GENERAL INTERPOLATION ******
!     ***********************************************************
!
   CALL CPUTIM (TIME)
   WRITE (IPR,900) TIME
   NPANLS = 0
!
   CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)
   DVL = DV
   LAY1SV = LAYR1
   PL = PAVE
   TL = TAVE
   WTOTL = 0.
   DO 10 MOL = 1, NMOL
      WTOTL = WTOTL+WK(MOL)
      WKSAV(MOL) = WK(MOL)
10 END DO
   WTOTL = WTOTL+WBROAD
   WN2SAV = WBROAD
!
!     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)
!
   IF (iaersl.ge.1 .and. iaersl.ne.5)                                &
   &                 WRITE (CYID,'(5A8)') (YID(I),I=3,7)
   CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)
   IF (iaersl.ge.1 .and. iaersl.ne.5)                                &
   &                 READ (CYID,'(5A8)') (YID(I),I=3,7)
   IF (JPATHL.GE.1) IPATHL = JPATHL
   PLAY = PAVE
   TLAY = TAVE
   DVK = DV
   LAYR1 = LAY1SV
   FACT = 1.
!
!     IF(IPATHL.EQ.2 .AND. IANT.EQ.0) FACT=2.
!
   IF (IPATHL.EQ.2.AND.IANT.EQ.0) STOP ' ABSINT: FACT=2.  '
   ATYPE = 9.999E09
   IF (DVK.EQ.DVL) ATYPE = 0.
   IF (DVL.GT.DVK) ATYPE = DVK/(DVL-DVK)+0.5
   IF (DVL.LT.DVK) ATYPE = -DVL/(DVK-DVL)-0.5
   IF (ATYPE.GT.0) STOP ' ABSINT; ATYPE GT 0 '
   WTOTK = 0.
   WRITE (IPR,905) LAYR1,LAYER,KFILE,LFILE,MFILE,ATYPE
   IEMIT = 0
   DO 20 MOL = 1, NMOL
      WTOTK = WTOTK+WK(MOL)*FACT
      WK(MOL) = WK(MOL)*FACT+WKSAV(MOL)
20 END DO
   WTOTK = WTOTK+WBROAD*FACT
   WBROAD = WBROAD*FACT+WN2SAV
   PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)
   TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)
   SECANT = 0.
   DV = DVL
!
!     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES
!
   CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)
   DVXM = DV
!
   IF (ATYPE.EQ.0.) THEN
!
!     1/1 RATIO ONLY
!
30    CONTINUE
!
      CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
      IF (KEOF.LE.0) GO TO 90
      CALL BUFIN (KFILE,KEOF,ODLAYR(1),NLIM)
!
      CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)
      CALL BUFIN (LFILE,LEOF,OLDOD(1),NLIMO)
!
      DO 40 I = 1, NLIM
         NEWOD(I) = ODLAYR(I)+OLDOD(I)
40    CONTINUE
      CALL ABSOUT (V1PO,V2PO,DVPO,NLIMO,1,MFILE,NPTS,NEWOD,NPANLS)
      GO TO 30
!
   ENDIF
!
!     ALL RATIOS EXCEPT 1/1
!
   DO 50 JP = 0, 100
      APG = JP
      P = 0.01*APG
!
!     THE FOLLOW ARE THE CONSTANTS FOR THE LAGRANGE 4 POINT
!     INTERPOLATION
!
      A1(JP) = -P*(P-1.0)*(P-2.0)/6.0
      A2(JP) = (P**2-1.0)*(P-2.0)*0.5
      A3(JP) = -P*(P+1.0)*(P-2.0)*0.5
      A4(JP) = P*(P**2-1.0)/6.0
50 END DO
!
   CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
   IF (KEOF.LE.0) GO TO 90
   CALL BUFIN (KFILE,KEOF,ODLAYR(1),NLIM)
!
   ODLAY(-2) = ODLAY(1)
   ODLAY(-1) = ODLAY(1)
   ODLAY(0) = ODLAY(1)
!
   RATDV = DVL/DVK
!
60 CONTINUE
!
   CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)
   IF (LEOF.LE.0) GO TO 90
   CALL BUFIN (LFILE,LEOF,OLDOD(1),NLIMO)
!
!     FJJ IS OFFSET BY 2. FOR ROUNDING PURPOSES
!
   FJ1DIF = (V1PO-V1P)/DVP+1.+2.
!
!     ***** BEGINNING OF LOOP THAT DOES MERGE  *****
!
   DO 80 II = 1, NLIMO

!
70    CONTINUE
!
      FJJ = FJ1DIF+RATDV* REAL(II-1)
      JJ = INT(FJJ)-2
!
      IF (JJ+2.GT.NLIM) THEN
         ODLAY(-2) = ODLAY(NLIM-2)
         ODLAY(-1) = ODLAY(NLIM-1)
         ODLAY(0) = ODLAY(NLIM)
         V1PST = V1P
         V2PST = V2P
         NLIMST = NLIM
!
         CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
!
         IF (KEOF.LE.0) THEN
            V1P = V1PST
            DVP = DVK
            V2P = V2PST+2.*DVP
            NLIM = NLIMST+2
            ODLAY(NLIM-1) = ODLAY(NLIM-2)
            ODLAY(NLIM) = ODLAY(NLIM-2)
         ELSE
            CALL BUFIN (KFILE,KEOF,ODLAYR(1),NLIM)
         ENDIF
!
         FJ1DIF = (V1PO-V1P)/DVP+1.+2.
         GO TO 70
      ENDIF
!
!     JP = (FJJ- REAL(JJ))*100. + 0.5 - 200.
!
      JP = (FJJ- REAL(JJ))*100.-199.5
      IF (JP.GT.100) THEN
         WRITE (IPR,910) JP,JJ,NLIM
         STOP
      ENDIF
!
!     INTERPOLATE THE OLD TRANSMISSION
!
      ODLAYI = A1(JP)*ODLAY(JJ-1)+A2(JP)*ODLAY(JJ)+ A3(JP)*ODLAY(JJ+ &
         1)+A4(JP)*ODLAY(JJ+2)
      IF (ODLAYI.LT.0.) ODLAYI = 0.
!
      NEWOD(II) = ODLAYI+OLDOD(II)
!
80 END DO
!
   CALL ABSOUT (V1PO,V2PO,DVPO,NLIMO,1,MFILE,NPTS,NEWOD,NPANLS)
!
   GO TO 60
!
90 CONTINUE
!
   CALL CPUTIM (TIME1)
   TIM = TIME1-TIME
   WRITE (IPR,915) TIME1,TIM
!
   RETURN
!
900 FORMAT ('0 THE TIME AT THE START OF ABSINT IS ',F12.3)
905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,   &
   &        ' MERGED WITH FILE ',I5,' ONTO FILE',I5,'  WITH XTYPE=',  &
   &        G15.5)
910 FORMAT ('0 JP, JJ, NLIM ',3I6)
915 FORMAT ('0 THE TIME AT THE END OF ABSINT IS ',F12.3/F12.3,        &
   &        ' SECS WERE REQUIRED FOR THIS MERGE ')
!
END SUBROUTINE ABSINT
!
!     ----------------------------------------------------------------
!
SUBROUTINE ABSOUT (V1PO,V2PO,DVPO,NLIMO,JLO,MFILE,NPTS,R1,NPANLS)
!
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE ABSOUT OUPUTS THE MERGED RESULT (R1) ONTO MFILE
!
   COMMON /ABSPNI/ V1P,V2P,DVP,NLIM
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   DIMENSION PNLHDR(2),R1(*)
!
   EQUIVALENCE (PNLHDR(1),V1P)
!
   V1P = V1PO
   V2P = V2PO
   DVP = DVPO
   NLIM = NLIMO
!
   NPANLS = NPANLS+1
   CALL BUFOUT (MFILE,PNLHDR(1),NPHDRF)
   CALL BUFOUT (MFILE,R1(JLO),NLIM)
   IF (NPTS.LE.0) GO TO 20
   IF (NPANLS.EQ.1) WRITE (IPR,900)
   WRITE (IPR,905)
   NNPTS = NPTS
   IF (NPTS.GT.(NLIM/2)+1) NNPTS = (NLIM/2)+1
   JHILIM = JLO+NLIM-NNPTS
   DO 10 I = 1, NNPTS
      J = JLO+I-1
      K = JHILIM+I-1
      VJ = V1P+ REAL(J-JLO)*DVP
      VK = V1P+ REAL(K-JLO)*DVP
      WRITE (IPR,910) J,VJ,R1(J),K,VK,R1(K)
10 END DO
20 CONTINUE
!
   RETURN
!
900 FORMAT ('0 ','LOCATION  WAVENUMBER',2X,'OPT DPTH',27X,            &
   &        'LOCATION   WAVENUMBER',2X,'OPT DPTH')
905 FORMAT (' ')
910 FORMAT (I8,2X,F12.6,1P,E15.7,0P,20X,I8,2X,F12.6,1P,E15.7)
!
END SUBROUTINE ABSOUT
!
!     ----------------------------------------------------------------
FUNCTION PLANCK (VI,XKT)
!
   USE phys_consts, ONLY: radcn1
   IMPLICIT REAL*8           (V)
!
!     FUNCTION BBFN CALCULATES BLACK BODY FN FOR WAVENUMBER VALUE VI
!     AND CALCULATES THE WAVENUMBER VALUE (VINEW) FOR NEXT BBFN CALC.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    01 JANUARY 2020
!
!                  IMPLEMENTATION:    
!                                     I.   POLONSKY
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
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  if xkt small, e.g. xkt = 0., trap and return
! to ensure smoothnenss in single precision make all calculation in double precision   
   VPLANCK  = 0d0
   if (XKT > 0.0) then
      VKT = real(XKT, kind=8)
      VIOKT = VI/VKT
      IF (VIOKT < 1d-2) THEN
         VPLANCK = (VI**2)*VKT/(1.+0.5*VIOKT)
      ELSEIF (VIOKT.LE.80d0) THEN
         VEXPNEG = EXP(-VIOKT)
         VPLANCK = (VI**3)*VEXPNEG/(1d0-VEXPNEG)
      ENDIF
   ENDIF
   PLANCK = VPLANCK*RADCN1
   RETURN
!
END FUNCTION PLANCK

FUNCTION PLANCK_DT (VI,XKT, BBVAL)
!
   USE phys_consts, ONLY: radcn2
   IMPLICIT REAL*8           (V)
!
!     FUNCTION PLANCK_DT calculates the derivative of the black body fn
!     analytically for wavenumber value VI
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!                  IMPLEMENTATION:    P.D. Brown, I. Polonsky
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
!
!
!     Incoming:
!                BBVAL      Planck function at VI
!                VI         Beginning wavenumber value for which BBdT is
!                XKT        Temperature in units of wavenumbers
!     Outgoing:
!                PLANCK_DT  Derivative of Planck function at VI

!
!     Using the exact function,
!
!       BBFN{prime} = BBFN*(planck*clight*gnu/(boltz*t*t))*
!                     1/[1-exp(-planck*clight*gnu/(boltz*t))]
!     where
!                   planck*clight/boltz = RADCN2
!               boltz*t/(planck*clight) = XKT
!             planck*clight*gnu/boltz*t = XVIOKT
!
!     and we can solve easily for t:
!                                     t = XKT*RADCN2.
!
! to ensure smoothnenss in single precision make all calculation in double precision   
  
   VPLANCK_DT = 0d0
   IF (XKT.GT.0.0) THEN
      VKT = real(XKT, kind=8)
      VIOKT = VI/VKT
      VTKELV  = VKT*real(RADCN2, kind=8)
      IF (VIOKT.LE.1d-2) THEN
         VPLANCK_DT = 1d0/(-VTKELV*(1d0+5d-1*VIOKT))
      ELSEIF (VIOKT.LE.80.0) THEN
         VPLANCK_DT = VIOKT/(VTKELV*(1d0-EXP(-VIOKT)))
      ENDIF
   ENDIF
   PLANCK_DT = VPLANCK_DT*BBVAL
   RETURN
END FUNCTION PLANCK_DT

!______________________________________________________________________
!     ----------------------------------------------------------------
!
FUNCTION BBFN (VI,DVI,V2I,XKT,VINEW,BBDEL,BBLAST)
!
   USE phys_consts, ONLY: radcn1
   IMPLICIT REAL*8           (V)
!
!     FUNCTION BBFN CALCULATES BLACK BODY FN FOR WAVENUMBER VALUE VI
!     AND CALCULATES THE WAVENUMBER VALUE (VINEW) FOR NEXT BBFN CALC.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    23 AUGUST 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
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
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   DATA I_1/1/
!
   DATA FACTOR / 0.00003 /
!
!     if xkt small, e.g. xkt = 0., trap and return

   if (xkt.lt.1.e-6) then
      bblast = 0.
      bbfn   = 0.

      bbdel  = 0.
      vinew  = 6.0E+5
      return
   endif

   XVI = VI
   XVIOKT = XVI/XKT
   EXPNEG = EXP(-XVIOKT)
   GNU2 = XVI*XVI
   BG2  = XVIOKT*XVIOKT
!
!     Initialize BBLAST for BBLAST negative
!
   IF (BBLAST.LT.0.) THEN
      IF (XKT.GT.0.0) THEN
         IF (XVIOKT.LE.0.01) THEN
            BBLAST = RADCN1*(XVI**2)*XKT/(1.+0.5*XVIOKT)
         ELSEIF (XVIOKT.LE.80.0) THEN
            BBLAST = RADCN1*(XVI**3)/(EXP(XVIOKT)-1.)
         ELSE
            BBLAST = 0.
         ENDIF
      ELSE
         BBLAST = 0.
      ENDIF
   ENDIF
!
!     SET BBFN EQUAL TO BLACK BODY FUNCTION AT VI
!
!     BBLAST IS BBFN(VI) FOR EACH SUBSEQUENT CALL
!
   BBFN = BBLAST
!
   INTVLS = 1
   DELTAV2 = V2I - VI
   IF (XKT.GT.0.0) THEN
!
      IF (XVIOKT.LE.0.01) THEN
         IF (VINEW.GE.0.0) THEN
            XDELT = (GNU2 * (4.+4.*XVIOKT + BG2))/ (10.*BG2 - 24.*   &
               XVIOKT + 8.)
            DELTAV = SQRT(ABS(FACTOR*XDELT))
            IF (DELTAV .GT. DELTAV2) DELTAV = DELTAV2
            INTVLS = ((DELTAV)/DVI) + 1.001
            INTVLS = MAX(INTVLS,I_1)
            VINEW = VI+DVI* REAL(INTVLS)
         ELSE
            VINEW = ABS(VINEW)
            INTVLS = ((VINEW-VI)/DVI) + 1.001
         ENDIF
         XVINEW = VINEW
!
         BBNEXT = RADCN1*(XVINEW**2)*XKT/(1.+0.5*XVINEW/XKT)
      ELSEIF (XVIOKT.LE.80.0) THEN
         IF (VINEW.GE.0.0) THEN
            FRONT = XVIOKT/(1.-EXPNEG)
            BOX = 3.0 - FRONT
            DELT2C = (1./GNU2)*(2.*BOX-FRONT*(1.+BOX-FRONT*EXPNEG))
            DELTAV = SQRT(ABS(FACTOR/DELT2C))
            IF (DELTAV .GT. DELTAV2) DELTAV = DELTAV2
            INTVLS = ((DELTAV)/DVI) + 1.001
            INTVLS = MAX(INTVLS,I_1)
            VINEW = VI+DVI* REAL(INTVLS)
         ELSE
            VINEW = ABS(VINEW)
!
!              The following IF test added for cases where
!              XVIOKT > 80 on one call, and XVIOKT < 80 on
!              the next call (numerical artifact causing
!              the change over the XVIOKT = 80 boundary)
!
            IF (VINEW.EQ.6.0E+05) THEN
               BBNEXT=0.
               BBDEL = (BBNEXT-BBFN)/ REAL(INTVLS)
               BBLAST = BBNEXT
               RETURN
            ENDIF
            INTVLS = ((VINEW-VI)/DVI) + 1.001
            INTVLS = MAX(INTVLS,I_1)
         ENDIF
         XVINEW = VINEW
!
         BBNEXT = RADCN1*(XVINEW**3)/(EXP(XVINEW/XKT)-1.)
      ELSE
         BBNEXT = 0.
         VINEW = 6.0E+5
      ENDIF
   ELSE
      BBNEXT = 0.
      VINEW = 6.0E+5
   ENDIF
!
   BBDEL = (BBNEXT-BBFN)/ REAL(INTVLS)
   BBLAST = BBNEXT
!
   RETURN
!
END FUNCTION BBFN
!
!______________________________________________________________________

FUNCTION  BBDTFN(BBVAL,VI,DVI,V2I,XKT,VDnew,BBdTdel,BBDTLAST)
!
   USE phys_consts, ONLY: radcn2
   IMPLICIT REAL*8           (V)
!
!     FUNCTION bbdTfn calculates the derivative of the black body fn
!     analytically for wavenumber value VI
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!                  IMPLEMENTATION:    P.D. Brown
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
!
   DATA I_1/1/
!
   DATA FACTOR / 0.003 /
!
!     Incoming:
!                BBVAL      Planck function at VI
!                VI         Beginning wavenumber value for which BBdT is
!                DVI        Wavenumber grid on which wil BBdT will be ap
!                V2I        Ending wavenumber value for which BBdT is re
!                XKT        Temperature in units of wavenumbers
!                VDnew      Wavenumber value for next calculation of fun
!                           Also Used as flag:
!                                         negative    specifies upper wa
!                                         positive    upper wavenumber v
!                BBdTlast   Last value of BBdT at wavenumber value VDdel
!                           Also used as a flag: negative indicates firs
!     Outgoing:
!                BBdTlast   Set to a nominal value for first trip throug
!                VDdel      Upper wavenumber value (value at which BBdT
!                BBdTdel    Incremental change of BBdT for every DVI inc
!                VI         Beginning wavenumber value for BBdT
!                BBdT       Derivative of Planck function at VI

   XVI = VI
   XVIOKT = XVI/XKT
   EXPNEG = EXP(-XVIOKT)
   GNU2 = XVI*XVI
   BG2  = XVIOKT*XVIOKT
   TKELV  = XKT*RADCN2
!
!     If first call, initialize BBDTLAST
!
!     For linear approximation,
!
!       BBFN{prime} = -BBFN/(t*(1 + 0.5*(planck*clight*gnu/boltz*t)))
!
!     Using the exact function,
!
!       BBFN{prime} = BBFN*(planck*clight*gnu/(boltz*t*t))*
!                     1/[1-exp(-planck*clight*gnu/(boltz*t))]
!     where
!                   planck*clight/boltz = RADCN2
!               boltz*t/(planck*clight) = XKT
!             planck*clight*gnu/boltz*t = XVIOKT
!
!     and we can solve easily for t:
!                                     t = XKT*RADCN2.
!
   IF (BBDTLAST.LT.0.) THEN
      IF (XKT.GT.0.0) THEN
         IF (XVIOKT.LE.0.01) THEN
            BBDTLAST = BBVAL/(-TKELV*(1.+0.5*XVIOKT))
         ELSEIF (XVIOKT.LE.80.0) THEN
            BBDTLAST = BBVAL*XVIOKT/(TKELV*(1-EXPNEG))
         ELSE
            BBDTLAST = 0.
         ENDIF
      ELSE
         BBDTLAST = 0.
      ENDIF
   ENDIF
!
!     Set BBAD equal to black body function derivative
!
!     BBDTLAST is BBAD(VI) for each subsequent call
!
   BBDTFN = BBDTLAST
!
   INTVLS = 1
   DELTAV2 = V2I - VI
   IF (XKT.GT.0.0) THEN
!
      IF (XVIOKT.LE.0.01) THEN
         IF (VDNEW.GE.0.0) THEN
            XDELT = (GNU2 * (4.+4.*XVIOKT + BG2))/ (10.*BG2 - 24.*   &
               XVIOKT + 8.)
            DELTAV = SQRT(ABS(FACTOR*XDELT))
            IF (DELTAV .GT. DELTAV2) DELTAV = DELTAV2
            INTVLS = (DELTAV)/DVI
            INTVLS = MAX(INTVLS,I_1)
            VDNEW = VI+DVI* REAL(INTVLS)
         ELSE
            VDNEW = ABS(VDNEW)
            INTVLS = (VDNEW-VI)/DVI
            INTVLS = MAX(INTVLS,I_1)
         ENDIF
         XVDNEW = VDNEW
         BBDTFN = BBVAL/(-TKELV*(1.+0.5*XVDNEW/XKT))
      ELSEIF (XVIOKT.LE.80.0) THEN
         IF (VDNEW.GE.0.0) THEN
            FRONT = XVIOKT/(1.-EXPNEG)
            BOX = 3.- FRONT
            DELT2C = (1./GNU2)*(2.*BOX-FRONT*(1.+BOX-FRONT*EXPNEG))
            DELTAV = SQRT(ABS(FACTOR/DELT2C))
            IF (DELTAV .GT. DELTAV2) DELTAV = DELTAV2
            INTVLS = (DELTAV)/DVI
            INTVLS = MAX(INTVLS,I_1)
            VDNEW = VI+DVI* REAL(INTVLS)
         ELSE
            VDNEW = ABS(VDNEW)
!
!              The following IF test added for cases where
!              XVIOKT > 80 on one call, and XVIOKT < 80 on
!              the next call (numerical artifact causing
!              the change over the XVIOKT = 80 boundary)
!
            IF (VDNEW.EQ.6.0E+05) THEN
               BBDTFN=0.
               VDNEW = VDNEW-DVI+0.00001
               BBdTdel = (BBDTFN-BBDTLAST)/ REAL(INTVLS)
               BBDTLAST = BBDTFN
               RETURN
            ENDIF
            INTVLS = (VDNEW-VI)/DVI
            INTVLS = MAX(INTVLS,I_1)
         ENDIF
         XVDNEW = VDNEW
         BBDTFN = BBVAL*(XVDNEW/XKT)/(TKELV*(1-EXP(-XVDNEW/XKT)))
      ELSE
         BBDTFN = 0.
         VDNEW = 6.0E+5
      ENDIF
   ELSE
      BBDTFN = 0.
      VDNEW = 6.0E+5
   ENDIF
!
   BBdTdel = (BBDTFN-BBDTLAST)/ REAL(INTVLS)
!
   VDNEW = VDNEW-DVI+0.00001
   BBDTLAST = BBDTFN
!
   RETURN
!
END FUNCTION BBDTFN
!
!     ----------------------------------------------------------------
!
FUNCTION EMISFN (VI,DVI,VINEM,EMDEL,EMLAST)
!
   USE lblparams, ONLY: NMAXCO
   IMPLICIT REAL*8           (V)
!
!     FUNCTION EMISFN CALCULATES BOUNDARY EMISSIVITY FOR WAVE NUMBER
!     VALUE CORRESPONDING TO VI AND VINEM, AND THEN CALCULATES THE
!     LINEAR CHANGE BETWEEN THE EMISSIVITY VALUES AT VI AND VINEM
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    23 AUGUST 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
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
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     ----------------------------------------------------------------
!     Parameter and common block for direct input of emission function
!     values
!
   !PARAMETER (NMAXCO=4040)
   COMMON /EMSFIN/ V1EMIS,V2EMIS,DVEMIS,NLIMEM,ZEMIS(NMAXCO)
!     ----------------------------------------------------------------
!
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,      &
   &    pad_3,angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!
   EQUIVALENCE (BNDEMI(1),A) , (BNDEMI(2),B) , (BNDEMI(3),C)
!
   DATA I_1/1/
!
   DATA FACTOR / 0.001 /
!
!     ***************************************************
!     Check for A < 0.  If so, use input values read in from file
!     "EMISSION"
!     ***************************************************
!
   IF (A.LT.0.) THEN
!
!        Determine elements of EMISSION function to use with
!        input frequency
!
      NELMNT = INT((VI-V1EMIS)/DVEMIS)
!
!        Test for bounds on EMISSION function
!
      IF ((NELMNT.LT.0).OR.(NELMNT.GE.NLIMEM-1)) THEN
         WRITE(*,*) 'Frequency range of calculation exceeded',       &
            ' emissivity input.'
         WRITE(*,*) ' VI = ',VI,' V1EMIS = ',V1EMIS,' V2EMIS = ',    &
            V2EMIS
         STOP 'ERROR IN EMISFN'
      ENDIF
!
!        Interpolate to obtain appropriate EMISSION value
!
      V1A = V1EMIS+DVEMIS*NELMNT
      V1B = V1EMIS+DVEMIS*(NELMNT+1)
      CALL LINTCO(V1A,ZEMIS(NELMNT+1),V1B,ZEMIS(NELMNT+2),VI,ZINT,   &
         ZDEL)
      EMISFN = ZINT
      VINEM = V1B
      EMDEL = ZDEL*DVI
      EMLAST = ZEMIS(NELMNT+1)
      RETURN
!
   ENDIF
!
!     ***************************************************
!     The following uses a quadratic formula for emission
!     ***************************************************
!
!     CHECK FOR CONSTANT E (INDEPENDENT OF VI)
!     IF CONSTANT RETURN LARGE VALUE FOR VINEM
!
   IF (B.EQ.0..AND.C.EQ.0.) THEN
      EMISFN = A
      VINEM = 9.99E+9
      EMDEL = 0.0
      EMLAST = EMISFN
      RETURN
   ENDIF
!
   XVI = VI
   IF (EMLAST.LT.0.) THEN
      EMLAST = A+B*XVI+C*XVI*XVI
   ENDIF
!
!     SET EMISFN EQUAL TO EMISSIVITY AT VI
!
!     EMLAST IS EMISFN(VI) FOR EACH SUBSEQUENT CALL
!
   EMISFN = EMLAST
!
   IF (VINEM.GE.0.0) THEN
      XVNEXT = XVI+FACTOR/ABS((B+2.*C*XVI))
      XVNEXT = MIN(XVNEXT,(XVI+DVI*2400))
      INTVLS = (XVNEXT-XVI)/DVI
      INTVLS = MAX(INTVLS,I_1)
      XVNEXT = XVI+DVI* REAL(INTVLS)
   ELSE
      XVNEXT = ABS(VINEM)
      INTVLS = (XVNEXT-XVI)/DVI
      INTVLS = MAX(INTVLS,I_1)
   ENDIF
!
   EMNEXT = A+B*XVNEXT+C*XVNEXT*XVNEXT
!
   EMDEL = (EMNEXT-EMISFN)/ REAL(INTVLS)
!
   VINEM = XVNEXT
   EMLAST = EMNEXT
!
   RETURN
!
END FUNCTION EMISFN
!
!     ----------------------------------------------------------------
!
FUNCTION REFLFN (VI,DVI,VINRF,RFDEL,RFLAST)
!
   USE lblparams, ONLY: NMAXCO
   IMPLICIT REAL*8           (V)
!
!     FUNCTION REFLFN CALCULATES BOUNDARY REFLECTIVITY FOR WAVE NUMBER
!     VALUE CORRESPONDING TO VI AND VINRF, AND THEN CALCULATES THE
!     LINEAR CHANGE BETWEEN THE REFLECTIVITY VALUES AT VI AND VINRF
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    23 AUGUST 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
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
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     ----------------------------------------------------------------
!     Parameter and common block for direct input of reflection
!     function values
!
!      PARAMETER (NMAXCO=4040)
   COMMON /RFLTIN/ V1RFLT,V2RFLT,DVRFLT,NLIMRF,ZRFLT(NMAXCO)
!     ----------------------------------------------------------------
!
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,      &
   &    pad_3,angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!
   EQUIVALENCE (BNDRFL(1),A) , (BNDRFL(2),B) , (BNDRFL(3),C)
!
   DATA FACTOR / 0.001 /
!
   DATA I_1/1/
!
!     ***************************************************
!     Check for A < 0.  If so, use input values read in from file
!     "REFLECTION"
!     ***************************************************
!
   IF (A.LT.0.) THEN
!
!        Determine elements of REFLECTION function to use with
!        input frequency
!
      NELMNT = INT((VI-V1RFLT)/DVRFLT)
!
!        Test for bounds on REFLECTION function
!
      IF ((NELMNT.LT.0).OR.(NELMNT.GE.NLIMRF-1)) THEN
         WRITE(*,*) 'Frequency range of calculation exceeded',       &
            ' reflectivity input.'
         WRITE(*,*) ' VI = ',VI,' V1RFLT = ',V1RFLT,' V2RFLT = ',    &
            V2RFLT
         STOP 'ERROR IN REFLFN'
      ENDIF
!
!        Interpolate to obtain appropriate reflection value
!
      V1A = V1RFLT+DVRFLT*NELMNT
      V1B = V1RFLT+DVRFLT*(NELMNT+1)
      CALL LINTCO(V1A,ZRFLT(NELMNT+1),V1B,ZRFLT(NELMNT+2),VI,ZINT,   &
         ZDEL)
      REFLFN = ZINT
      VINRF = V1B
      RFDEL = ZDEL*DVI
      RFLAST = ZRFLT(NELMNT+1)
      RETURN
!
   ENDIF
!
!     ***************************************************
!     The following uses a quadratic formula for emission
!     ***************************************************
!
!     CHECK FOR CONSTANT R (INDEPENDENT OF VI)
!     IF CONSTANT RETURN LARGE VALUE FOR VINRF
!
   IF (B.EQ.0..AND.C.EQ.0.) THEN
      REFLFN = A
      VINRF = 9.99E+9
      RFDEL = 0.0
      RFLAST = REFLFN
      RETURN
   ENDIF
!
   XVI = VI
   IF (RFLAST.LT.0.) THEN
      RFLAST = A+B*XVI+C*XVI*XVI
   ENDIF
!
!     SET REFLFN EQUAL TO REFLECTIVITY AT VI
!
!     RFLAST IS REFLFN(VI) FOR EACH SUBSEQUENT CALL
!
   REFLFN = RFLAST
!
   IF (VINRF.GE.0.0) THEN
      XVNEXT = XVI+FACTOR/ABS((B+2.*C*XVI))
      XVNEXT = MIN(XVNEXT,(XVI+DVI*2400))
      INTVLS = (XVNEXT-XVI)/DVI
      INTVLS = MAX(INTVLS,I_1)
      XVNEXT = XVI+DVI* REAL(INTVLS)
   ELSE
      XVNEXT = ABS(VINRF)
      INTVLS = (XVNEXT-XVI)/DVI
      INTVLS = MAX(INTVLS,I_1)
   ENDIF
!
   RFNEXT = A+B*XVNEXT+C*XVNEXT*XVNEXT
!
   RFDEL = (RFNEXT-REFLFN)/ REAL(INTVLS)
!
   VINRF = XVNEXT
   RFLAST = RFNEXT
!
   RETURN
!
END FUNCTION REFLFN
!
!     ----------------------------------------------------------------
!
SUBROUTINE EMIN (V1P,V2P,DVP,NLIM,KFILE,EM,EMB,TR,KEOF,NPANLS)
!
   USE lblparams, ONLY: NN_TBL, dbg, od_lo
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE EMIN INPUTS OPTICAL DEPTH VALUES FROM KFILE AND
!       CALCULATES SOURCE FUNCTION FOR THE LAYER.
!       THIS VERSION WORKS FOR AEROSOLS AND NLTE.
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    03 MARDCH 2006
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
!                                     M.W. SHEPHARD
!
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     131 Hartwell Ave., Lexington, MA. 02421
!
!----------------------------------------------------------------------
!
!               WORK SUPPORTED BY:    THE ARM PROGRAM
!                                     OFFICE OF ENERGY RESEARCH
!                                     DEPARTMENT OF ENERGY
!
!
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF
   COMMON /RMRG/ XKT,XKTA,XKTB,SECNT
!
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,      &
   &    pad_3,angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!
   DIMENSION PNLHDR(2),EM(*),EMB(*),TR(*)
!
   EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4)
   EQUIVALENCE (PNLHDR(1),V1PBF)
   EQUIVALENCE (FSCDID(4),IAERSL)
!
   data itbl_calc/-99/, aa /0.278/
!
   CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
   IF (KEOF.LE.0) RETURN
   CALL BUFIN (KFILE,KEOF,TR(1),NLIMBF)
!
!     TR contains the optical depths at this stage
!
   IF (IHIRAC.EQ.4) CALL BUFIN (KFILE,KEOF,EM(1),NLIMBF)
!
!     EM contains the optical depth corrections for nlte at this stage
!
   IF (NPANLS.LT.1) then
      if (IAERSL.EQ.0 .or. iaersl.eq.5) then
         WRITE (IPR,900)
      else
         WRITE (IPR,905)
      endif
   ENDIF
!
   EXT = 0.
   ADEL = 0.
   RADFN0 = 0.
   RDEL = 0.
   BB = 0.
   BBDEL = 0.
   BBA = 0.
   BBDLA = 0.
   BBB = 0.
   BBDLB = 0.
!
   V1P = V1PBF
   V2P = V2PBF
   DVP = DVPBF
   NLIM = NLIMBF
   VI = V1P-DVP
   VIDV = VI
   VIBB = VI
   VAER = VI
   VDUM = VI
   BBLAST = -1.
   BBLXTA = -2.
   BBLXTB = -3.
   RDLAST = -1.
   BBDUM = -4.
   RDDUM = -1.
   NLIM1 = 0
   NLIM2 = 0
!
   rec_6 = 1./6.
!
! **********************************************************************
!
   IF (IAERSL.EQ.0 .or. iaersl.eq.5) THEN
      IAFBB = -1
   ELSE
      RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDDUM)
      EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
      IF (VITST.LT.VAER) then
         IAFBB = 1
      ELSE
         IAFBB = 2
      ENDIF   
   ENDIF
!
!     - THIS SECTION TREATS THE CASE WHERE THE LAYER CONTRIBUTES
!       TO THE RADIATIVE TRANSFER ONLY ONCE
!
!     - WITH XKTA=0 THIS ALGORITHM REVERTS TO THE ORIGINAL
!
! **********************************************************************
!
   IF (XKTB.LE.0.) THEN
!
!     - THIS SECTION TREATS THE LTE CASE
!
      IF (IHIRAC.NE.4) THEN
         if (dbg(1)) then
            print *,'EMIN::XKTB.LE.0. .AND. IHIRAC /= 4: LTE CASE:  CHECKED'
            dbg(1) = .false.
         END IF
!
         VI = V1P
         bb_dif = 0.
10       NLIM1 = NLIM2+1
!
         IF (IAFBB.EQ.-1) THEN
            NLIM2 = NLIM
         ELSEIF (IAFBB.EQ.1) THEN
            RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
            EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
         ELSEIF (IAFBB.EQ.2) THEN
            EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            VITST = -VIDV
            RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
         ENDIF
!
!
         DO I = NLIM1, NLIM2
            BB = PLANCK(VI,XKT)
            IF (XKTA.GT.0.) THEN
               bb_dif = PLANCK(VI,XKTA)-BB
            ENDIF            
!
            ODVI = TR(I)+EXT*RADFN0
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
            if (odvi .lt. -od_lo) odvi = -od_lo
!
            tr_i = exp(-odvi)
            if (odvi .le. od_lo) then
               f_i = rec_6*odvi
               !
            else
               f_i = 1. - 2.*(tr_i/(tr_i-1.) + 1./odvi   )
               !
            end if
            tr(i) = tr_i
            em(i) = (1.-tr_i) * (bb + bb_dif * f_i)
!
!              Increment interpolation values
!
            EXT = EXT+ADEL
            RADFN0 = RADFN0+RDEL
!
            VI = VI+ DVP
         END DO
!
         IF (NLIM2.LT.NLIM) GO TO 10
      ELSE
!
!     - THIS SECTION TREATS THE NLTE CASE
!
         if (dbg(2)) then
            print *,'EMIN::XKTB.LE.0. .AND. IHIRAC == 4: NLTE CASE: NOT CHECKED'
            dbg(2) = .false.
         ENDIF
         VI = V1P
30       NLIM1 = NLIM2+1
!
         IF (IAFBB.EQ.-1) THEN
            NLIM2 = NLIM
         ELSEIF (IAFBB.EQ.1) THEN
            RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
            EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
         ELSEIF (IAFBB.EQ.2) THEN
            EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            VITST = -VIDV
            RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
         ENDIF
!
         DO I = NLIM1, NLIM2
            BB = PLANCK(VI,XKT)
            IF (XKTA.GT.0.) THEN
               bb_dif = PLANCK(VI,XKTA)-BB
            ENDIF            
!              tr(i) contains the layer optical depths at this stage

            ODVI = TR(I)+EXT*RADFN0
!
!              em(i) contains the ratio differences from
!                                              lte of the state populati
            c_nlte = em(i)
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
            if (odvi .lt. -od_lo) odvi = -od_lo
            tr_i = exp(-odvi)
            if (odvi .le. od_lo) then
               f_i = rec_6*odvi
               !
               abs_i = (odvi   - c_nlte) * (1.- 0.5 * odvi)
            else
               f_i = 1. - 2.*(tr_i/(tr_i-1.)     + 1./odvi   )
               !
               abs_i = (1. - c_nlte/odvi) * (1.-tr_i)
            end if
            tr(i) = tr_i
            em(i) = abs_i * (bb + bb_dif * f_i)
!
!              Increment interpolation values
!
            EXT = EXT+ADEL
            RADFN0 = RADFN0+RDEL
!
            VI = VI + DVP
         END DO
         !
         IF (NLIM2.LT.NLIM) GO TO 30
!
      ENDIF
! --------------------------------------------------------------
   ELSE
! --------------------------------------------------------------
!
!     - THIS SECTION TREATS THE CASE WHERE THE LAYER CONTRIBUTES
!       TO THE RADIATIVE TRANSFER TWICE:
!
!     - FOR TANGENT PATHS AND FOR THE CASE OF THE REFLECTED ATMOSPHERE
!
      IF (IHIRAC.NE.4) THEN
!
!     - THIS SECTION TREATS THE LTE CASE
!
         bb_dif_a = 0.
         bb_dif_b = 0.
         VI = V1P
50       NLIM1 = NLIM2+1
!
         IF (IAFBB.EQ.-1) THEN
            NLIM2 = NLIM
         ELSEIF (IAFBB.EQ.1) THEN
            RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
            EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
         ELSEIF (IAFBB.EQ.2) THEN
            EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            VITST = -VIDV
            RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
         ENDIF
!ccc
!     This calculation  is for specular reflection for the downwelling
!ccc

         if (surf_refl .eq. 's') then
!
            if (dbg(3)) then
               print *,'EMIN::XKTB.LE.0. .AND. IHIRAC /= 4: LTE CASE: specular:  CHECKED'
               dbg(3) = .false.
            ENDIF
            DO I = NLIM1, NLIM2
               BB = PLANCK(VI,XKT)
               IF (XKTA.GT.0.) THEN
                  bb_dif_a = PLANCK(VI,XKTA)-BB
               ENDIF            
               IF (XKTB.GT.0.) THEN
                  bb_dif_b = PLANCK(VI,XKTB)-BB
               ENDIF            
!
               ODVI = TR(I)+EXT*RADFN0
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
               if (odvi .lt. -od_lo) odvi = -od_lo
               tr_i = exp(-odvi)
               if (odvi .le. od_lo) then
                  f_i = rec_6*odvi
                  !
               else
                  f_i = 1. - 2.*(tr_i/(tr_i-1.) + 1./odvi   )
                  !
               end if
               tr(i) = tr_i
               abs_i = (1.-tr_i)
               em(i) = abs_i * (bb + bb_dif_a * f_i)
               emb(i)= abs_i * (bb + bb_dif_b * f_i)
!
!     Increment interpolation values
!
               EXT = EXT+ADEL
               RADFN0 = RADFN0+RDEL
               VI = VI + DVP
            END DO
   !
         elseif (surf_refl .eq. 'l') then
            if (dbg(4)) then
               print *,'EMIN::XKTB.LE.0. .AND. IHIRAC /= 4: LTE CASE: Lambertian:  NOT FIXED'
               dbg(4) = .false.
            endif
!ccc
!     The following calculation is for an approximation to the
!     downwelling flux for application to Lambertian surfaces. The
!     'diffusivity' approximation is used with the assumption that the
!     dowmwelling flux is isotropic and that the surface scatters
!     isotropically.  with a value obtained from the
!     The value of the diffusivity angle corresponds to a secant of 1.67
!     diffuse_fac is the factor that is used to scale the optical depth
!     on the 'observer side' of the path to that required for the 'back
!ccc
            DO I = NLIM1, NLIM2
               BB = PLANCK(VI,XKT)
               IF (XKTA.GT.0.) THEN
                  bb_dif_a = PLANCK(VI,XKTA)-BB
               ENDIF            
               IF (XKTB.GT.0.) THEN
                  bb_dif_b = PLANCK(VI,XKTB)-BB
               ENDIF            
!
               ODVI = TR(I)+EXT*RADFN0
               odvi_d = diffuse_fac * odvi
!
!       for odvi outside the range of the table,  set optical depth to bo
!
               if (odvi .lt. -od_lo) odvi = -od_lo
               tr_i   = exp(-odvi)
               tr_d_i = exp(-odvi_d)
               if (odvi .le. od_lo) then
                  f_i   = rec_6*odvi
                  f_d_i = rec_6*odvi_d
                  !
               else
                  f_i   = 1. - 2.*(tr_i  /(tr_i  -1.) + 1./odvi   )
                  f_d_i = 1. - 2.*(tr_d_i/(tr_d_i-1.) + 1./odvi_d )
                  !
               end if
               TR(i) = tr_i
               em(i)  = (1.-tr_i  ) * (bb + bb_dif_a * f_i)
               emb(i) = (1.-tr_d_i) * (bb + bb_dif_b * f_d_i)
!---
!
!     Increment interpolation values
!
               EXT = EXT+ADEL
               RADFN0 = RADFN0+RDEL
               VI = VI + DVP
            END DO
   !
         else
            WRITE (IPR,906) surf_refl
            STOP 'INVALID SURFACE REFLECTIVITY FLAG'
         endif
!
         IF (NLIM2.LT.NLIM) GO TO 50
!
      ELSE
!
!     - THIS SECTION TREATS THE CASE OF NLTE
!
         bb_dif_a = 0.
         bb_dif_b = 0.
         VI = V1P
70       NLIM1 = NLIM2+1
!
         IF (IAFBB.EQ.-1) THEN
            NLIM2 = NLIM
         ELSEIF (IAFBB.EQ.1) THEN
            RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
            EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
         ELSEIF (IAFBB.EQ.2) THEN
            EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            VITST = -VIDV
            RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
         ENDIF
!
!ccc
!     This calculation  is for specular reflection for the downwelling
!ccc

         if (surf_refl .eq. 's') then
            if (dbg(5)) then
               print *,'EMIN::XKTB.LE.0. .AND. IHIRAC == 4: NLTE CASE: specular: CHECKED'
               dbg(5) = .false.
            ENDIF

            DO I = NLIM1, NLIM2
               BB = PLANCK(VI,XKT)
               IF (XKTA.GT.0.) THEN
                  bb_dif_a = PLANCK(VI,XKTA)-BB
               ENDIF            
               IF (XKTB.GT.0.) THEN
                  bb_dif_b = PLANCK(VI,XKTB)-BB
               ENDIF            
!     tr(i) contains the layer optical depths at this stage

               ODVI = TR(I)+EXT*RADFN0
!
!     em(i) contains the ratio differences from lte of the state populat
!
               c_nlte = em(i)
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
               if (odvi .lt. -od_lo) odvi = -od_lo
               tr_i   = exp(-odvi)
               if (odvi .le. od_lo) then
                  f_i   = rec_6*odvi
                  !
                  abs_i = (odvi   - c_nlte) * (1.- 0.5 * odvi)
               else
                  f_i   = 1. - 2.*(tr_i  /(tr_i  -1.) + 1./odvi   )
                  !
                  abs_i = (1. - c_nlte/odvi) * (1.-tr_i)
               end if
               TR(i) = tr_i
               em(i) = abs_i * (bb + bb_dif_a * f_i)
               emb(i)= abs_i * (bb + bb_dif_b * f_i)

!
!     Increment interpolation values
!
               EXT = EXT+ADEL
               RADFN0 = RADFN0+RDEL
               VI = VI + DVP
            END DO
   
         elseif (surf_refl .eq. 'l') then
            if (dbg(6)) then
               print *,'EMIN::XKTB.LE.0. .AND. IHIRAC == 4: NLTE CASE: Lambertian: NOT CHECKED'
               dbg(6) = .false.
            endif

!ccc
!     The following calculation is for an approximation to the
!     downwelling flux for application to Lambertian surfaces. The
!     'diffusivity' approximation is used with the assumption that the
!     dowmwelling flux is isotropic and that the surface scatters
!     isotropically.  with a value obtained from the
!     The value of the diffusivity angle corresponds to a secant of 1.67
!     diffuse_fac is the factor that is used to scale the optical depth
!     on the 'observer side' of the path to that required for the 'back
!ccc
            DO I = NLIM1, NLIM2
               BB = PLANCK(VI,XKT)
               IF (XKTA.GT.0.) THEN
                  bb_dif_a = PLANCK(VI,XKTA)-BB
               ENDIF            
               IF (XKTB.GT.0.) THEN
                  bb_dif_b = PLANCK(VI,XKTB)-BB
               ENDIF            
!     tr(i) contains the layer optical depths at this stage

               ODVI = TR(I)+EXT*RADFN0
               odvi_d = diffuse_fac * odvi
!
!     em(i) contains the ratio differences from  lte of the state popula
!
               c_nlte = em(i)
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
               if (odvi .lt. -od_lo) odvi = -od_lo
               tr_i   = exp(-odvi)
               tr_d_i = exp(-odvi_d)
               if (odvi .le. od_lo) then
                  f_i   = rec_6*odvi
                  f_d_i = rec_6*odvi_d
                  !
                  abs_i = (odvi   - c_nlte) * (1.- 0.5 * odvi)
               else
                  f_i   = 1. - 2.*(tr_i  /(tr_i  -1.) + 1./odvi   )
                  f_d_i = 1. - 2.*(tr_d_i/(tr_d_i-1.) + 1./odvi_d )
                  !
                  abs_i  = (1. - c_nlte/odvi  ) * (1.-tr_i)
               end if
               TR(i)  = tr_i
               abs_d_i= (1. - c_nlte/odvi_d) * (1.-tr_d_i)
               em(i)  = abs_i   * (bb + bb_dif_a * f_i)
               emb(i) = abs_d_i * (bb + bb_dif_b * f_d_i)
!
!     Increment interpolation values
!
               EXT = EXT+ADEL
               RADFN0 = RADFN0+RDEL
               VI = VI + DVP
            END DO
   
         else

            WRITE (IPR,906) surf_refl
            STOP 'INVALID SURFACE REFLECTIVITY FLAG'

         endif
!
         IF (NLIM2.LT.NLIM) GO TO 70
!
      ENDIF
   ENDIF
!
!     ---------------------------------------------------------------
!
   RETURN
!
900 FORMAT ('0EMISSION AND TRANSMISSION  (MOLECULAR) ')
905 FORMAT ('0EMISSION AND TRANSMISSION (AEROSOLS EFFECTS INCLUDED)')
906 FORMAT (' THE SURFACE REFLECTIVITY FLAG OF ', A1, 'IS NOT VALID')
907 FORMAT (' THE SURFACE REFLECTIVITY FLAG OF: ', A1)
!
END SUBROUTINE EMIN
!
!     ----------------------------------------------------------------
!
SUBROUTINE EMINIT (NPTS,MFILE,JPATHL,TBND)
!
   USE phys_consts, ONLY: radcn2
   use lblparams, ONLY: dbg, od_lo
   IMPLICIT REAL*8           (V)
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    14 AUGUST 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
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
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   COMMON NEWEM(2410),NEWTR(2410)
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYER,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /EMIHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYDUM,YI1,YID(10),LSTWDF
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,pad_3,&
   &    angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!
   COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
!
   common /rddn_pnl/ v1pdn,v2pdn,dvpdn,nlimdn,                       &
   &                  rddn_sfc(2400),tr_sfc(2400)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /RMRG/ XKT,XKTA,XKTB,SECNT
!
   CHARACTER*40 CEXT,CYID
!
   DIMENSION EMLAYB(2410)
   DIMENSION XFILHD(2),OPNLHD(2)
   DIMENSION EMLAYR(2),TRALYR(2)
   dimension xdn_dum(2),tr_dum(2)
!
   EQUIVALENCE (XFILHD(1),XID(1)) , (OPNLHD(1),V1PO)
   EQUIVALENCE (NEWEM(1),EMLAYR(1)) , (NEWTR(1),TRALYR(1)),          &
   &            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),               &
   &            (FSCDID(11),IMRG) , (FSCDID(16),LAYR1)
!
   REAL NEWEM,NEWTR
!
!
! *********************************************************************
! ****  THIS SUBROUTINE COMPUTES THE EMISSION FOR THE FIRST LAYER  ****
! *********************************************************************
!
!     TBND IS THE BOUNDARY BLACK BODY TEMPERATUE
!
!     IPATHL = -1 IS FOR downwelling radiance REFLECTED ATMOSPHERE
!     IPATHL =  0 IS FOR THE HORIZONTAL PATH CASE (HOMOGENEOUS LAYER)
!     IPATHL =  1 IS for upwelling radiance (TO DENSER LAYERS)
!     IPATHL =  2 IS FOR THE SYMMETRIC TANGENT PATH radiance
!     IPATHL =  3 IS for downwelling radiance (TO LESS DENSE LAYERS)
!     IPATHL = 31 is for downwelling radiance (to more dense layers)
!
   CALL CPUTIM (TIME)
!
!      ** NOTE ON IPATHL =2
!           THE TANGENT MERGE ROUTINES ARE DIVIDED INTO ANTERIOR (1ST)
!           AND POSTERIOR (2ND) LAYER CROSSINGS.  THIS RECOGNITION IS
!           TRIGGERED BY THE VALUE OF "IANT"
!
!          IF  IANT.EQ. 1  THEN POSTERIOR MERGE
!          IF  IANT.EQ. 0  THEN NORMAL MERGE
!          IF  IANT.EQ.-1  THEN ANTERIOR MERGE
!
   WRITE (IPR,900) TIME
   NPANLS = 0
   CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)
   IF (JPATHL.GE.1) IPATHL = JPATHL
   PLAY = PAVE
   TLAY = TAVE
!
!     FOR AEROSOL RUNS, MOVE EXTID INTO YID
!
   IF (iaersl.ge.1 .and. iaersl.ne.5) THEN
      WRITE (CEXT,'(10A4)') EXTID
      WRITE (CYID,'(5A8)') (YID(I),I=3,7)
      CYID(19:40) = CEXT(19:40)
      READ (CYID,'(5A8)') (YID(I),I=3,7)
   ENDIF
!
   imrg_sfc = 0

! check to see if this is upwelling case: open 'RDDN_SFC file

   IF (imrg.eq.41 .and. IBPROP.EQ.1.AND.IPATHL.EQ.1) then
      k_rddn = 27

      OPEN(UNIT=k_rddn,FILE='RDDN_SFC',FORM='UNFORMATTED', STATUS=   &
         'OLD',iostat=ios)

      if (ios.gt.0) then
         stop 'not able to open RDDN_SFC'
      endif

      imrg_sfc = 1
!        buffer past file header
      call bufin(k_rddn,keof,xdn_dum(1),1)

   endif
!
!     IF BOUNDARY PROPERTIES ARE SUPPLIED, AND DOWNWARD LOOKING
!     CASE; SET IPATHL TO REFLECTED ATMOSPHERE CASE
!
   IF (imrg.ne.41 .and. IBPROP.EQ.1.AND.IPATHL.EQ.1) IPATHL = -1

   IEMIT = 1
   FACT = 1.
   TIMEM = 0.0
   IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.
   DO 10 MOL = 1, NMOL
      WK(MOL) = WK(MOL)*FACT
10 END DO
   WBROAD = WBROAD*FACT
   LAYR1 = LAYER
   WRITE (IPR,905) LAYR1,LAYER,KFILE,MFILE
   CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)
   DVXM = DV
   XKT = TAVE/RADCN2
   XKTBND = TBND/RADCN2
   IF (IPATHL.EQ.-1) THEN
      XKTA = TZU/RADCN2
      XKTB = TZL/RADCN2
   ENDIF
   IF (IPATHL.EQ.0) THEN
      XKTA = 0.
      XKTB = 0.
   ENDIF
   IF (IPATHL.EQ.1) THEN
      XKTA = TZU/RADCN2
      XKTB = 0.
   ENDIF
   IF (IPATHL.EQ.2) THEN
      XKTA = TZU/RADCN2
      XKTB = TZL/RADCN2
   ENDIF
   IF (IPATHL.EQ.3 .or. ipathl.eq.31) THEN
      XKTA = TZL/RADCN2
      XKTB = 0.
   ENDIF
   WRITE (IPR,910) IPATHL,IANT
!
20 CONTINUE
!
   CALL CPUTIM (TIMEM1)
   CALL EMIN (V1PO,V2PO,DVPO,NLIMO,KFILE,EMLAYR,EMLAYB,TRALYR,    &
      KEOF, NPANLS)
   CALL CPUTIM (TIMEM2)
   TIMEM = TIMEM+TIMEM2-TIMEM1
   IF (KEOF.LE.0) GO TO 80
   VI = V1PO-DVPO
   VIDVBD = VI
   VIDVEM = VI
   VIDVRF = VI
   BBLAST = -1.
   EMLAST = -1.
   RFLAST = -1.
   IF (IPATHL.EQ.2.AND.IANT.EQ.0) THEN
      DO 30 J = 1, NLIMO
         TRJ = TRALYR(J)
         NEWEM(J) = EMLAYR(J)+EMLAYB(J)*TRJ
         TRALYR(J) = TRALYR(J)*TRJ
30    CONTINUE
   ELSEIF ((IPATHL.EQ.1).AND.(TBND.GT.0.)) THEN
!
      if (dbg(7)) then
         print *,'EMINIT (IPATHL.EQ.1).AND.(TBND.GT.0.)::  CHECKED'
         dbg(7) = .false.
      endif
      NLIM1 = 0
      NLIM2 = 0
      EMDUM = 0.
      BBDUM = 0.
      EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMDUM)
      IEMBB = 1
!
      VI = V1PO
40    NLIM1 = NLIM2+1
!
      EMISIV = EMISFN(VI,DVPO,VIDV,EMDEL,EMLAST)
!
      IF (VIDV.GE.9.E+4) THEN
         NLIM2 = NLIMO+1
      ELSE
         NLIM2 = (VIDV-V1PO)/DVPO+1.001
      ENDIF
      NLIM2 = MIN(NLIM2,NLIMO)
!
      DO J = NLIM1, NLIM2
         V=V1PO+ REAL(J-1)*DVPO
         BB = PLANCK(VI,XKTBND)
         NEWEM(J) = EMLAYR(J)+TRALYR(J)*BB*EMISIV
!
!           Increment interpolation value
!
         EMISIV = EMISIV+EMDEL
         VI = VI+ DVPO
      END DO
!
      IF (NLIM2.LT.NLIMO) GO TO 40
!
   ELSEIF (IPATHL.EQ.-1 .or. imrg_sfc.gt.0) then
!
      if (dbg(8)) then
         print *,'EMINIT (IPATHL.EQ.1) .or. imrg_sfc.gt.0:: NOT CHECKED', IPATHL, imrg_sfc
         dbg(8) = .false.
      endif
      NLIM1 = 0
      NLIM2 = 0
      EMDUM = 0.
      RFDUM = 0.
      BBDUM = 0.
      EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMDUM)
      REFLCT = REFLFN(VI,DVPO,VIDVRF,RFDEL,RFDUM)
      IF (VIDVEM.LE.VIDVRF) THEN
         IEMBB = 1
      ELSE
         IEMBB = 2
      ENDIF
      !
      VI = V1PO
60    NLIM1 = NLIM2+1
!
      IF (IEMBB.EQ.1) THEN
         EMISIV = EMISFN(VI,DVPO,VIDV,EMDEL,EMLAST)
         VIDVRF = -VIDV
         REFLCT = REFLFN(VI,DVPO,VIDVRF,RFDEL,RFLAST)
      ELSE
         REFLCT = REFLFN(VI,DVPO,VIDV,RFDEL,RFLAST)
         VIDVEM = -VIDV
         EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMLAST)
      ENDIF
!
      IF (VIDV.GE.9.E+4) THEN
         NLIM2 = NLIMO+1
      ELSE
         NLIM2 = (VIDV-V1PO)/DVPO+1.001
      ENDIF
      NLIM2 = MIN(NLIM2,NLIMO)
!
!     check to see if contribution from reflected precalculated
!     downwelling radiance at the surface (RDDN_SFC) is to be included

      if (imrg_sfc .eq. 0 ) then

         DO J = NLIM1, NLIM2
            V=V1PO+ REAL(J-1)*DVPO
            BB = PLANCK(VI,XKTBND)
            NEWEM(J) = EMLAYR(J)+EMLAYB(J)*REFLCT*TRALYR(J)+ TRALYR( &
               J)*BB*EMISIV
!
!              Increment interpolation value
!
            EMISIV = EMISIV+EMDEL
            REFLCT = REFLCT+RFDEL
            VI = VI+ DVPO
         END DO
         !
      else

         call bufin(k_rddn,keof,xdn_dum(1),nphdrf)
         call bufin(k_rddn,keof,rddn_sfc(1),nlimdn)
         call bufin(k_rddn,keof,tr_dum(1),1)

         DO J = NLIM1, NLIM2
            V=V1PO+ REAL(J-1)*DVPO
            BB = PLANCK(VI,XKTBND)
            NEWEM(J) = EMLAYR(J)+rddn_sfc(j)*REFLCT*TRALYR(J)+       &
               TRALYR(J)*BB*EMISIV
!
!                  Increment interpolation value
!
            EMISIV = EMISIV+EMDEL
            REFLCT = REFLCT+RFDEL
            VI = VI+ DVPO
         END DO

      endif
!
      IF (NLIM2.LT.NLIMO) GO TO 60
!
   ENDIF
   CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,NEWEM,NEWTR,MFILE,NPTS,NPANLS)
   GO TO 20

80 CALL CPUTIM (TIME1)
   TIME = TIME1-TIME
   WRITE (IPR,915) TIME,TIMEM
!
   RETURN
!
900 FORMAT (' TIME AT THE START OF --EMINIT-- ',F10.3)
905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,                &
   &        '0 INPUT FILE =',I5,' OUTPUT FILE =',I5)
910 FORMAT ('0 IPATHL AND IANT',2I5)
915 FORMAT (' TIME REQUIRED FOR --EMINIT-- ',F10.3,                   &
   &        ' --EMIN-- ',F10.3)
!
END SUBROUTINE EMINIT
!
!     ----------------------------------------------------------------
!
SUBROUTINE RADMRG (NPTS,LFILE,MFILE,JPATHL,TBND)
!
   USE phys_consts, ONLY: radcn2
   IMPLICIT REAL*8           (V)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    8 APRIL 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
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
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   COMMON RADN(2410),TRAN(2410),RADO(0:5000)
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),      &
   &               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,  &
   &               EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,      &
   &    pad_3,angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!
   COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
   COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /XME/ v1_pad,v2_pad,dv_pad,n_pad,TRAO(0:5000)
   COMMON /RMRG/ XKT,XKTA,XKTB,SECNT
!
   DIMENSION RADLYB(2410)
   DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)
   DIMENSION A1(10),A2(10),A3(10),A4(10)
   DIMENSION RADLYR(2),TRALYR(2),RADOI(2),TRAOI(2)
   DIMENSION WKSAV(35)
!
   CHARACTER*40 CYID
!
   EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                 &
   &            (OPNLHD(1),V1PO)
   EQUIVALENCE (RADO(1),RADOI(1)) , (TRAO(1),TRAOI(1)),              &
   &            (RADN(1),RADLYR(1)) , (TRAN(1),TRALYR(1)),            &
   &            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),               &
   &            (FSCDID(16),LAYR1)
!
!
!
!
!      ************************************************************
!      ****** THIS SUBROUTINE DOES LAYER MERGE FOR RADIANCE  ******
!      ************************************************************
!
!     IPATHL = -1 IS FOR downwelling radiance REFLECTED ATMOSPHERE
!     IPATHL =  0 IS FOR THE HORIZONTAL PATH CASE (HOMOGENEOUS LAYER)
!     IPATHL =  1 IS for upwelling radiance (TO DENSER LAYERS)
!     IPATHL =  2 IS FOR THE SYMMETRIC TANGENT PATH radiance
!     IPATHL =  3 IS for downwelling radiance (TO LESS DENSE LAYERS)
!     IPATHL = 31 is for downwelling radiance (to more dense layers)
!
!
!      ** NOTE ON IPATHL = 2
!            THE TANGENT MERGE ROUTINES ARE DIVIDED INTO ANTERIOR (1ST)
!            AND POSTERIOR (2ND) LAYER CROSSINGS   THIS RECOGNITION IS
!            TRIGGERED BY THE VALUE OF "IANT"
!
!          IF  IANT.EQ. 1  THEN POSTERIOR MERGE
!          IF  IANT.EQ. 0  THEN NORMAL MERGE
!          IF  IANT.EQ.-1  THEN ANTERIOR MERGE
!
   CALL CPUTIM (TIME)
   WRITE (IPR,900) TIME
   NPANLS = 0
   TIMEM = 0.0
   TIMRD = 0.0
   TIMOT = 0.0
!
   CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)
   LAY1SV = LAYR1
   DVL = DV
   PL = PAVE
   TL = TAVE
   WTOTL = 0.
!
   DO 10 MOL = 1, NMOL
      WTOTL = WTOTL+WK(MOL)
      WKSAV(MOL) = WK(MOL)
10 END DO
!
   WTOTL = WTOTL+WBROAD
   WN2SAV = WBROAD
!
!     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)
!
   IF (iaersl.ge.1 .and. iaersl.ne.5)                                &
   &                 WRITE (CYID,'(5A8)') (YID(I),I=3,7)
   CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)
   IF (iaersl.ge.1 .and. iaersl.ne.5)                                &
   &                 READ (CYID,'(5A8)') (YID(I),I=3,7)
!
   IF (JPATHL.GE.1) IPATHL = JPATHL
   PLAY = PAVE
   TLAY = TAVE
!
!     IF BOUNDARY PROPERTIES ARE SUPPLIED, AND DOWNWARD LOOKING
!     CASE; SET IPATHL TO REFLECTED ATMOSPHERE CASE
!
   IF (IBPROP.EQ.1.AND.IPATHL.EQ.1) IPATHL = -1
   TAVK = TAVE
   DVK = DV
   FACT = 1.
   IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.
!
   IF (DVL.EQ.DVK) THEN
      ITYPE = 0
   ELSEIF (DVL.GT.DVK) THEN
      ITYPE = DVK/(DVL-DVK)+0.5
   ELSE
!
!     DVL.LT.DVK
!
      ITYPE = -INT(DVL/(DVK-DVL)+0.5)
   ENDIF
   IF (ITYPE.LT.0) STOP ' RADMRG; ITYPE LT 0 '
!
   WTOTK = 0.
   LAYR1 = LAY1SV
   WRITE (IPR,905) LAYR1,LAYER,KFILE,LFILE,MFILE
   IEMIT = 1
   DO 20 MOL = 1, NMOL
      WTOTK = WTOTK+WK(MOL)*FACT
      WK(MOL) = WK(MOL)*FACT+WKSAV(MOL)
20 END DO
   WTOTK = WTOTK+WBROAD*FACT
   WBROAD = WBROAD*FACT+WN2SAV
   PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)
   TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)
   SECANT = 0.
!
!     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES
!
   CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)
   DVXM = DV
   XKT = TAVK/RADCN2
!
   WRITE (IPR,910) IPATHL,IANT
!
   IF (IPATHL.EQ.-1) THEN
      XKTA = TZU/RADCN2
      XKTB = TZL/RADCN2
   ELSEIF (IPATHL.EQ.1) THEN
      XKTA = TZU/RADCN2
      XKTB = 0.
   ELSEIF (IPATHL.EQ.2) THEN
      XKTA = TZU/RADCN2
      XKTB = TZL/RADCN2
   ELSEIF (IPATHL.EQ.3 .or. ipathl.eq.31) THEN
      XKTA = TZL/RADCN2
      XKTB = 0.
   ELSE
      STOP ' RADMRG; IPATHL '
   ENDIF
!
   ATYPE = ITYPE
   LL = ITYPE+1
   AP = 1.0/(ATYPE+1.0)
!
!     A1, A2, A3 AND A4 ARE THE CONSTANTS
!     FOR THE LAGRANGE 4 POINT INTERPOLATION
!
   DO 30 JPG = 1, ITYPE
      APG = JPG
      IPL = JPG+1
      P = 1.0-(AP*APG)
      A1(IPL) = -P*(P-1.0)*(P-2.0)/6.0
      A2(IPL) = (P**2-1.0)*(P-2.0)*0.5
      A3(IPL) = -P*(P+1.0)*(P-2.0)*0.5
      A4(IPL) = P*(P**2-1.0)/6.0
30 END DO
!
!     *** BEGINNING OF LOOP THAT DOES MERGE ***
!
   NPE = 0
   RADO(0) = 0.0
   TRAO(0) = 0.0
   V1PO = 0.0
   V2PO = 0.0
   DVPO = 0.0
!
40 CONTINUE
!
   CALL CPUTIM (TIMEM1)
   CALL EMIN (V1P,V2P,DVP,NLIM,KFILE,RADLYR,RADLYB,TRALYR,KEOF,      &
   &           NPANLS)
   CALL CPUTIM (TIMEM2)
   TIMEM = TIMEM+TIMEM2-TIMEM1
   IF (KEOF.LE.0) GO TO 80
   II = 1
!
   IF (V2PO.LE.V2P+DVPO) THEN
50    CALL CPUTIM (TIMEM1)
      CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)
      IF (LEOF.LE.0) GO TO 60
      CALL BUFIN (LFILE,LEOF,RADOI(NPE+1),NLIMO)
      CALL BUFIN (LFILE,LEOF,TRAOI(NPE+1),NLIMO)
      CALL CPUTIM (TIMEM2)
      TIMRD = TIMRD+TIMEM2-TIMEM1
      NPE = NLIMO+NPE
      IF (V2PO.LE.V2P+DVPO) GO TO 50
   ENDIF
!
!     ZERO POINT OF FIRST PANEL
!
60 IF (RADO(0).EQ.0.0.AND.TRAO(0).EQ.0.0) THEN
      RADO(0) = RADO(1)
      TRAO(0) = TRAO(1)
   ENDIF
!
!     END POINT OF LAST PANEL
!
   IF (V2PO+DVPO.GE.V2) THEN
      RADO(NPE+1) = RADO(NPE)
      RADO(NPE+2) = RADO(NPE)
      TRAO(NPE+1) = TRAO(NPE)
      TRAO(NPE+2) = TRAO(NPE)
   ENDIF
!
   NPL = 1
!
!     NPL IS LOCATION OF FIRST ELEMENT ON ARRAYS RADO AND TRAO
!
   ipath_flg = ipathl
   if (ipathl .eq. -1) then
      if (surf_refl .eq. 's') ipath_flg = -10
      if (surf_refl .eq. 'l') ipath_flg = -11
   endif
!
   CALL RADNN (RADN,TRAN,RADO,TRAO,RADLYB,NLIM,V1P,DVP,              &
   &           IPATH_flg,A1,A2,A3,A4,LL,NPL)
!
   CALL CPUTIM (TIMEM1)
!
   IF (TBND.GT.0.) CALL EMBND (V1P,V2P,DVP,NLIM,RADN,TRAN,TBND)
!
   CALL EMOUT (V1P,V2P,DVP,NLIM,RADN,TRAN,MFILE,NPTS,NPANLS)
   CALL CPUTIM (TIMEM2)
   TIMOT = TIMOT+TIMEM2-TIMEM1
!
!     NPL IS NOW LOCATION OF FIRST ELEMENT TO BE USED FOR NEXT PASS
!
   IPL = -1
   DO 70 NL = NPL, NPE
      IPL = IPL+1
      RADO(IPL) = RADO(NL)
      TRAO(IPL) = TRAO(NL)
70 END DO
!
   NPE = IPL
!
   GO TO 40
80 CONTINUE
!
   CALL CPUTIM (TIME1)
   TIM = TIME1-TIME
   WRITE (IPR,915) TIME1,TIM,TIMEM,TIMRD,TIMOT
!
   RETURN
!
!
900 FORMAT ('0 THE TIME AT THE START OF RADMRG IS ',F12.3)
905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,   &
   &        ' MERGED WITH FILE ',I5,' ONTO FILE',I5)
910 FORMAT ('0 IPATHL AND IANT',2I5)
915 FORMAT ('0 THE TIME AT THE END OF RADMRG IS ',F12.3/F12.3,        &
   &        ' SECS WERE REQUIRED FOR THIS MERGE  - EMIN - ',          &
   &        F12.3,' - READ - ',F12.3,' - EMOUT - ',F12.3)
!
END SUBROUTINE RADMRG
!
!     ----------------------------------------------------------------
!
SUBROUTINE RADNN (RADLYR,TRALYR,RADO,TRAO,RADLYB,NLIM,            &
&                  V1P,DVP,IPATH_flg,A1,A2,A3,A4,LL,NPL)
!
   USE lblparams, ONLY : NDIM, ND2
   IMPLICIT REAL*8           (V)
!
!     THIS SUBROUTINE CALCULATES THE NEW RADIANCE AND TRANSMISSION
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    8 APRIL 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!                       ALGORITHM:    R.D. WORSHAM
!                                     S.A. CLOUGH
!                                     J.L. MONCET
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
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,      &
   &    pad_3,angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
! if this changes, make sure it is changed in subroutines in xmerge.f
!      parameter (ndim=2410, nd2=5000)

!
   DIMENSION RADLYR(NDIM),TRALYR(NDIM),RADO(0:ND2),TRAO(0:ND2),      &
   &          RADLYB(NDIM),A1(*),A2(*),A3(*),A4(*)
!
   DATA I_1/1/
!
   dimension tr_diffus(-100:1100), tr_diffus_dif(-100:1100)
!
   save i_trdiffus, tr_diffus, tr_diffus_dif, xnpts
!
   data i_trdiffus/-99/
!
!     lambertian

!ccc
!     The following call sets up a table to obtain the transmittance
!     from the layer to the surface at the diffusivity angle
!     (backside) from the transmittance from the surface to the layer on
!     'observer side' of the path.  The quantity 'diffuse_fac'
!     is the appropriate factor to be used for this purpose.  Note that
!     optical depth for the downwelling radiance has already been
!     properly computed in subroutine emin.
!
!     i_trdiffus       flag to establish if the table needs to be set up
!     tr_diffus        table of diffuse transmittances given the
!                             'observer side' transmittance
!     tr_diffus_dif    first differences in table tr_diffus for interpol
!     xnpoints         number of values in the table
!     diffuse_fac      factor to obtainn oberver side' optical depth fro
!                              'backside' optical depth
!
!ccc

   if (ipath_flg .eq. -11 .and. i_trdiffus.eq.-99)                   &
   &     call tr_diffus_tbl                                           &
   &        (i_trdiffus, tr_diffus, tr_diffus_dif, xnpts, diffuse_fac)
!
!
   LLM1 = LL-1
   LLM1 = MAX(LLM1,I_1)
!
!     LOOPING OVER POINTS WITH SAME WEIGHTS
!
   DO 110 NL = 1, LL
      IPL = (NPL+NL-1)-LLM1
      IF (NL.GT.1) IPL = IPL-1
!
      IF (NL.EQ.1) THEN
!
!     EXACT FREQUENCY - NO INTERPOLATION
!
         IF (IPATH_FLG.EQ.1) THEN
!
            DO 10 I = NL, NLIM, LL
               IPL = IPL+LLM1
               RADLYR(I) = TRALYR(I)*RADO(IPL)+RADLYR(I)
               TRALYR(I) = TRALYR(I)*TRAO(IPL)
10          CONTINUE
!
         ELSEIF (IPATH_FLG.EQ.2) THEN
!
            DO 20 I = NL, NLIM, LL
               IPL = IPL+LLM1
               TRTEMP = TRALYR(I)*TRAO(IPL)
               RADLYR(I) = TRALYR(I)*RADO(IPL)+RADLYR(I)+ RADLYB(I)* &
                  TRTEMP
               TRALYR(I) = TRALYR(I)*TRTEMP
20          CONTINUE
!
         ELSEIF (IPATH_FLG.EQ.3) THEN
!
            DO 30 I = NL, NLIM, LL
               IPL = IPL+LLM1
               RADLYR(I) = RADO(IPL)+RADLYR(I)*TRAO(IPL)
               TRALYR(I) = TRALYR(I)*TRAO(IPL)
30          CONTINUE
!
         ELSEIF (IPATH_FLG.EQ.31) THEN
!
            DO 31 I = NL, NLIM, LL
               IPL = IPL+LLM1
               RADLYR(I) = RADLYR(I) + TRALYR(IPL)* RADO(IPL)
               TRALYR(I) = TRALYR(I)*TRAO(IPL)
31          CONTINUE
!
         ELSEIF (IPATH_FLG.EQ.-10 .or. IPATH_FLG.EQ.-11) THEN
!
            VI = V1P-DVP
            DVI = DVP* REAL(LL)
            VIDVRF = VI
            RFLAST = -1.
            NLIM1 = 0
            NLIM2 = NL-1
!
40          NLIM1 = NLIM2+1
!
            VI = V1P+ REAL(NLIM1-1)*DVP
            REFLCT = REFLFN(VI,DVI,VIDVRF,RFDEL,RFLAST)
!
            IF (VIDVRF.GE.9.E+4) THEN
               NLIM2 = NLIM+1
            ELSE
               NLIM2 = (VIDVRF+DVI-DVP-V1P)/DVP+1.001
            ENDIF
            NLIM2 = MIN(NLIM2,NLIM)
!
!              Test to make sure LL divides evenly into (NLIM2-NLIM1+1)
!
            NRMNDR = MOD(NLIM2-NLIM1,LL)
            IF ((NRMNDR.NE.0).AND.(NLIM2+(LL-NRMNDR).LT.2400)) THEN
               NLIM2 = NLIM2+(LL-NRMNDR)
            ENDIF
!
            IF (IPATH_FLG.EQ.-10) THEN
!
!     specular
!
               DO 50 I = NLIM1, NLIM2, LL
                  IPL = IPL+LLM1
                  TRTEMP = TRALYR(I)*TRAO(IPL)
                  RADLYR(I) = TRALYR(I)*RADO(IPL)+RADLYR(I)+ RADLYB( &
                     I)*TRAO(IPL)*TRTEMP*REFLCT
                  TRALYR(I) = TRTEMP
!
!     Increment interpolation values
!
                  REFLCT = REFLCT+RFDEL
                  ILAST = I
50             CONTINUE
!
            endif
!
            IF (IPATH_FLG.EQ.-11) THEN
!
!     Lambertian
!
!     the effective flux transmittance, diffus_tr, is obtained by interp
!     from the table tr_diffus.  The effective argument is TRAO(IPL).
!     comments at the begiining of this subroutine for further informati
!
               DO 52 I = NLIM1, NLIM2, LL
                  IPL = IPL+LLM1

                  trao_ipl = TRAO(IPL)
                  TRTEMP = TRALYR(I)*trao_ipl

                  itd = xNpts*trao_ipl
                  diff = xNpts*trao_ipl - REAL(itd)
                  diffus_tr = tr_diffus(itd)+diff*tr_diffus_dif(  &
                     itd)
!
                  RADLYR(I) = TRALYR(I)*RADO(IPL)+RADLYR(I)+      &
                     RADLYB(I)*diffus_tr*TRTEMP*REFLCT

                  TRALYR(I) = TRTEMP
!
!           Increment interpolation values
!
                  REFLCT = REFLCT+RFDEL
                  ILAST = I
52             CONTINUE
            endif
!
            IF (NLIM2.LT.NLIM) THEN
               NLIM2 = ILAST + LL-1
               GO TO 40
            ENDIF
!
         ENDIF
!
!     NOT EXACT FREQUENCY - INTERPOLATE RESULT
!
      ELSE
!
         A1N = A1(NL)
         A2N = A2(NL)
         A3N = A3(NL)
         A4N = A4(NL)
!
         IF (IPATH_FLG.EQ.1) THEN
            DO 60 I = NL, NLIM, LL
               IPL = IPL+LLM1
!
!     INTERPOLATE THE OLD RADIANCE
!
               RADLYR(I) = TRALYR(I)*(A1N*RADO(IPL-1)+A2N*RADO(&
                  IPL)+ A3N*RADO(IPL+1)+A4N*RADO(IPL+2))+RADLYR(I)
!
!     INTERPOLATE THE OLD TRANSMISSION
!
               TRALYR(I) = TRALYR(I)*(A1N*TRAO(IPL-1)+A2N*TRAO(&
                  IPL)+ A3N*TRAO(IPL+1)+A4N*TRAO(IPL+2))
60          CONTINUE
!
         ELSEIF (IPATH_FLG.EQ.2) THEN
!
            DO 70 I = NL, NLIM, LL
               IPL = IPL+LLM1
!
!     INTERPOLATE THE OLD TRANSMISSION
!
               TRTEMP = TRALYR(I)*(A1N*TRAO(IPL-1)+A2N*TRAO(   &
                  IPL)+ A3N*TRAO(IPL+1)+A4N*TRAO(IPL+2))
!
!     INTERPOLATE THE OLD RADIANCE
!
               RADLYR(I) = TRALYR(I)*(A1N*RADO(IPL-1)+A2N*RADO(&
                  IPL)+ A3N*RADO(IPL+1)+A4N*RADO(IPL+2))+ RADLYR( &
                  I)+RADLYB(I)*TRTEMP
               TRALYR(I) = TRALYR(I)*TRTEMP
70          CONTINUE
!
         ELSEIF (IPATH_FLG.EQ.3) THEN
!
            DO 80 I = NL, NLIM, LL
               IPL = IPL+LLM1
!
!     INTERPOLATE THE OLD TRANSMISSION
!
               TRAOI = A1N*TRAO(IPL-1)+A2N*TRAO(IPL)+ A3N*TRAO(&
                  IPL+1)+A4N*TRAO(IPL+2)
!
!     INTERPOLATE THE OLD RADIANCE
!
               RADLYR(I) = A1N*RADO(IPL-1)+A2N*RADO(IPL)+      &
                  A3N*RADO(IPL+1)+A4N*RADO(IPL+2)+ RADLYR(I)*     &
                  TRAOI
               TRALYR(I) = TRALYR(I)*TRAOI
80          CONTINUE
!
         ELSEIF (IPATH_FLG.EQ.-10 .or. IPATH_FLG.EQ.-11) THEN
!
            VI = V1P-DVP
            DVI = DVP* REAL(LL)
            VIDVRF = VI
            RFLAST = -1
            NLIM1 = 0
            NLIM2 = NL-1
!
90          NLIM1 = NLIM2+1
!
            VI = V1P+ REAL(NLIM1-1)*DVP
            REFLCT = REFLFN(VI,DVI,VIDVRF,RFDEL,RFLAST)
!
            IF (VIDVRF.GE.9.E+4) THEN
               NLIM2 = NLIM+1
            ELSE
               NLIM2 = (VIDVRF+DVI-DVP-V1P)/DVP+1.001
            ENDIF
            NLIM2 = MIN(NLIM2,NLIM)
!
!              Test to make sure LL divides evenly into (NLIM2-NLIM1+1)
!
            NRMNDR = MOD(NLIM2-NLIM1,LL)
            IF ((NRMNDR.NE.0).AND.(NLIM2+(LL-NRMNDR).LT.2400)) &
               THEN
               NLIM2 = NLIM2+(LL-NRMNDR)
            ENDIF
!
            IF (IPATH_FLG.EQ.-10) THEN
!
!              specular
!
               DO 100 I = NLIM1, NLIM2, LL
                  IPL = IPL+LLM1
!
!     INTERPOLATE THE OLD TRANSMISSION
!
                  TRAOI = A1N*TRAO(IPL-1)+A2N*TRAO(IPL)+       &
                     A3N*TRAO(IPL+1)+A4N*TRAO(IPL+2)
!
                  TRTEMP = TRALYR(I)*TRAOI
!
!     INTERPOLATE THE OLD RADIANCE
!
                  RADLYR(I) = TRALYR(I)* (A1N*RADO(IPL-1)+A2N* &
                     RADO(IPL)+ A3N*RADO(IPL+1)+A4N*RADO(IPL+2))+ &
                     RADLYR(I)+RADLYB(I)*TRAOI*TRTEMP*REFLCT
                  TRALYR(I) = TRTEMP
!
!           Increment interpolation values
!
                  REFLCT = REFLCT+RFDEL
                  ILAST = I
100            CONTINUE
!
            endif
!
            IF (IPATH_FLG.EQ.-11) THEN
!
!              Lambertian
!
!     the effective flux transmittance, diffus_tr, is obtained by interp
!     from the table tr_diffus.  The effective argument is TRAOI.  See t
!     comments at the begiining of this subroutine for further informati
!
               DO 102 I = NLIM1, NLIM2, LL
                  IPL = IPL+LLM1
!
!     INTERPOLATE THE OLD TRANSMISSION
!
                  TRAOI = A1N*TRAO(IPL-1)+A2N*TRAO(IPL)+       &
                     A3N*TRAO(IPL+1)+A4N*TRAO(IPL+2)
!
!
                  itd = xNpts*traoi
                  diff = xNpts*traoi - REAL(itd)
                  diffus_tr = tr_diffus(itd)+diff*             &
                     tr_diffus_dif(itd)
!
                  TRTEMP = TRALYR(I)*TRAOI
!
!     INTERPOLATE THE OLD RADIANCE
!
                  RADLYR(I) = TRALYR(I)* (A1N*RADO(IPL-1)+A2N* &
                     RADO(IPL)+ A3N*RADO(IPL+1)+A4N*RADO(IPL+2))+ &
                     RADLYR(I)+RADLYB(I)*diffus_tr*TRTEMP*REFLCT
                  TRALYR(I) = TRTEMP
!
!           Increment interpolation values
!
                  REFLCT = REFLCT+RFDEL
                  ILAST = I
102            CONTINUE
!
            endif
!
            IF (NLIM2.LT.NLIM) THEN
               NLIM2 = ILAST + LL -1
               GO TO 90
            ENDIF
!
         ENDIF
!
      ENDIF
!
110 END DO
!
   NPL = IPL
!
   RETURN
!
END SUBROUTINE RADNN
!
!     ----------------------------------------------------------------
!
SUBROUTINE tr_diffus_tbl                                          &
&     (i_trdiffus, tr_diffus, tr_diffus_dif, xnpts, diffuse_fac)
!
!ccc
!     This subroutine sets up a table to obtain the transmittance
!     from the layer to the surface at the diffusivity angle
!     (backside) from the transmittance from the surface to the layer on
!     'observer side' of the path.  The quantity 'diffuse_fac'
!     is the appropriate factor to be used for this purpose.  Note that
!     optical depth for the downwelling radiance has already been been
!     properly computed in subroutine emin.
!
!     i_trdiffus       flag to establish if the table needs to be set up
!     tr_diffus        table of diffuse transmittances given the
!                             'observer side' transmittance
!     tr_diffus_dif    first differences in table tr_diffus for interpol
!     xnpoints         number of values in the table
!     diffuse_fac      factor to obtainn oberver side' optical depth fro
!                              'backside' optical depth
!
!ccc
!
   dimension tr_diffus(-100:1100), tr_diffus_dif(-100:1100)
!
   i_trdiffus = +11
   Npts  = 500
   Nxtra = Npts/10
   xNpts =  REAL(Npts)
   trdiff = 1./xNpts

   do 5 i=0,Npts+Nxtra
      tri =  REAL(i)*trdiff
      tr_diffus (i) = tri**diffuse_fac
      tr_diffus_dif (i) = (tri+trdiff)**diffuse_fac-tr_diffus (i)
5  continue
!
   do 6 i=-Nxtra,-1
      tr_diffus (i) = 0.
      tr_diffus_dif (i) = 0.
6  continue

   return

END SUBROUTINE tr_diffus_tbl
!
!     ----------------------------------------------------------------
!
SUBROUTINE RADINT (NPTS,LFILE,MFILE,JPATHL,TBND)
!
   USE phys_consts, ONLY: radcn2
   IMPLICIT REAL*8           (V)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    5 APRIL 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
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
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   COMMON RADN(2410),TRAN(2410),RADLYR(-1:4818)
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),      &
   &               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,  &
   &               EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
   COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /XMI/ v1_pqd,v2_dum,dv_dum,n_dum,TRALYR(-1:4818)
   COMMON /RMRG/ XKT,XKTA,XKTB,SECNT
!
   DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)
   DIMENSION A1(0:100),A2(0:100),A3(0:100),A4(0:100)
   DIMENSION RADO(2),TRAO(2)
   DIMENSION WKSAV(35)
!
   DIMENSION RADLYB(-1:4818)
!
   CHARACTER*40 CYID
!
   EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                 &
   &            (OPNLHD(1),V1PO)
   EQUIVALENCE (RADN(1),RADO(1)) , (TRAN(1),TRAO(1)),                &
   &            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),               &
   &            (FSCDID(16),LAYR1)
!
!     ************************************************************
!     ****** THIS SUBROUTINE DOES LAYER MERGE FOR EMISSION  ******
!     ****** USING FOUR POINT GENERAL INTERPOLATION         ******
!     ************************************************************
!
   CALL CPUTIM (TIME)
   WRITE (IPR,900) TIME
   NPANLS = 0
   TIMEM = 0.0
   TIMRD = 0.0
   TIMTB = 0.0
   TIMOT = 0.0
!
   CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)
   DVL = DV
   LAY1SV = LAYR1
   PL = PAVE
   TL = TAVE
   WTOTL = 0.
   DO 10 MOL = 1, NMOL
      WTOTL = WTOTL+WK(MOL)
      WKSAV(MOL) = WK(MOL)
10 END DO
   WTOTL = WTOTL+WBROAD
   WN2SAV = WBROAD
!
!     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)
!
   IF (iaersl.ge.1 .and. iaersl.ne.5)                                &
   &     WRITE (CYID,'(5A8)') (YID(I),I=3,7)
   CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)
   IF (iaersl.ge.1 .and. iaersl.ne.5)                                &
   &     READ (CYID,'(5A8)') (YID(I),I=3,7)
   IF (JPATHL.GE.1) IPATHL = JPATHL
   PLAY = PAVE
   TLAY = TAVE
   XKT = TAVE/RADCN2
   XKTA = TZU/RADCN2
   XKTB = 0.
   DVK = DV
   LAYR1 = LAY1SV
   FACT = 1.
   IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.
   ATYPE = 9.999E09
   IF (DVK.EQ.DVL) ATYPE = 0.
   IF (DVL.GT.DVK) ATYPE = DVK/(DVL-DVK)+0.5
   IF (DVL.LT.DVK) ATYPE = -DVL/(DVK-DVL)-0.5
!
!     IF (ATYPE .GT. 0) STOP  ' RADINT; ATYPE GT 0 '
!
   WTOTK = 0.
   WRITE (IPR,905) LAYR1,LAYER,KFILE,LFILE,MFILE,ATYPE
   IEMIT = 1
   DO 20 MOL = 1, NMOL
      WTOTK = WTOTK+WK(MOL)*FACT
      WK(MOL) = WK(MOL)*FACT+WKSAV(MOL)
20 END DO
   WTOTK = WTOTK+WBROAD*FACT
   WBROAD = WBROAD*FACT+WN2SAV
   PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)
   TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)
   SECANT = 0.
   DV = DVL
!
!     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES
!
   CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)
   DVXM = DV
!
   IF (ATYPE.EQ.0.) THEN
!
!     1/1 RATIO ONLY
!
30    CONTINUE
      CALL CPUTIM (TIMEM1)
      CALL EMIN (V1P,V2P,DVP,NLIM,KFILE,RADLYR(1),RADLYB(1), TRALYR( &
         1),KEOF,NPANLS)
      CALL CPUTIM (TIMEM2)
      TIMEM = TIMEM+TIMEM2-TIMEM1
      IF (KEOF.LE.0) GO TO 110
      CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)
      CALL BUFIN (LFILE,LEOF,RADO(1),NLIMO)
      CALL BUFIN (LFILE,LEOF,TRAO(1),NLIMO)
      CALL CPUTIM (TIMEM3)
      TIMRD = TIMRD+TIMEM3-TIMEM2
      DO 40 I = 1, NLIM
         RADN(I) = RADO(I)+RADLYR(I)*TRAO(I)
         TRAN(I) = TRALYR(I)*TRAO(I)
40    CONTINUE
      CALL CPUTIM (TIMEM1)
      IF (TBND.GT.0.) CALL EMBND (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,    &
         TBND)
      CALL CPUTIM (TIMEM2)
      TIMTB = TIMTB+TIMEM2-TIMEM1
      CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,MFILE,NPTS,NPANLS)
      CALL CPUTIM (TIMEM3)
      TIMOT = TIMOT+TIMEM3-TIMEM2
      GO TO 30
!
   ENDIF
!
!     ALL RATIOS EXCEPT 1/1
!
   DO 50 JP = 0, 100
      APG = JP
      P = 0.01*APG
!
!     THE FOLLOW ARE THE CONSTANTS FOR THE LAGRANGE 4 POINT
!     INTERPOLATION
!
      A1(JP) = -P*(P-1.0)*(P-2.0)/6.0
      A2(JP) = (P**2-1.0)*(P-2.0)*0.5
      A3(JP) = -P*(P+1.0)*(P-2.0)*0.5
      A4(JP) = P*(P**2-1.0)/6.0
50 END DO
!
!     *** BEGINNING OF LOOP THAT DOES MERGE  ***
!
   NPE = 0
   RADLYR(0) = 0.0
   TRALYR(0) = 0.0
   V1P = 0.0
   V2P = 0.0
   DVP = 0.0
   V1PO = 0.0
   V2PO = 0.0
   DVPO = 0.0
   KEOF = 1
!
60 CONTINUE
!
   CALL CPUTIM (TIMEM1)
   CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)
   IF (LEOF.LE.0) GO TO 110
   CALL BUFIN (LFILE,LEOF,RADO(1),NLIMO)
   CALL BUFIN (LFILE,LEOF,TRAO(1),NLIMO)
   CALL CPUTIM (TIMEM2)
   TIMRD = TIMRD+TIMEM2-TIMEM1
   II = 1
!
   IF (V2P.LE.V2PO+DVP .AND. KEOF.GT.0) THEN
70    CALL CPUTIM (TIMEM2)
      CALL EMIN (V1P,V2P,DVP,NLIM,KFILE,RADLYR(NPE+1),RADLYB(NPE+1), &
         TRALYR(NPE+1),KEOF,NPANLS)
      CALL CPUTIM (TIMEM3)
      TIMEM = TIMEM+TIMEM3-TIMEM2
      IF (KEOF.LE.0) GO TO 80
      V1P = V1P- REAL(NPE)*DVP
      NPE = NLIM+NPE
      IF (V2P.LE.V2PO+DVP) GO TO 70
   ENDIF
!
!     ZERO POINT OF FIRST PANEL
!
80 IF (RADLYR(0).EQ.0.0.AND.TRALYR(0).EQ.0.0) THEN
      TRALYR(-1) = TRALYR(1)
      RADLYR(-1) = RADLYR(1)
      RADLYB(-1) = RADLYB(1)
      TRALYR(0) = TRALYR(1)
      RADLYR(0) = RADLYR(1)
      RADLYB(0) = RADLYB(1)
   ENDIF
!
!     END POINT OF LAST PANEL
!
   IF (V2P+DVP.GE.V2) THEN
      TRALYR(NPE+1) = TRALYR(NPE)
      RADLYR(NPE+1) = RADLYR(NPE)
      RADLYB(NPE+1) = RADLYB(NPE)
      TRALYR(NPE+2) = TRALYR(NPE)
      RADLYR(NPE+2) = RADLYR(NPE)
      RADLYB(NPE+2) = RADLYB(NPE)
   ENDIF
!
!     NPL IS LOCATION OF FIRST ELEMENT ON ARRAYS RADO AND TRAO
!
   NPL = 1
!
   RATDV = DVL/DVK
!
!     FJJ IS OFFSET BY 2. FOR ROUNDING PURPOSES
!
   FJ1DIF = (V1PO-V1P)/DVP+1.+2.
!
!     ***** BEGINNING OF LOOP THAT DOES MERGE  *****
!
   DO 90 II = 1, NLIMO
!
      FJJ = FJ1DIF+RATDV* REAL(II-1)
      JJ = INT(FJJ)-2
!
      JP = (FJJ- REAL(JJ))*100.-199.5
!
!     INTERPOLATE THE OLD EMISSION
!
      RADN(II) = RADO(II)+(A1(JP)*RADLYR(JJ-1)+A2(JP)*RADLYR(JJ)+    &
         A3(JP)*RADLYR(JJ+1)+A4(JP)*RADLYR(JJ+2))*TRAO(II)
!
!     INTERPOLATE THE OLD TRANSMISSION
!
      TRAN(II) = (A1(JP)*TRALYR(JJ-1)+A2(JP)*TRALYR(JJ)+ A3(JP)*     &
         TRALYR(JJ+1)+A4(JP)*TRALYR(JJ+2))*TRAO(II)
!
90 END DO
!
   NPL = JJ-1
!
   CALL CPUTIM (TIMEM1)
   IF (TBND.GT.0.) CALL EMBND (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,TBND)
   CALL CPUTIM (TIMEM2)
   TIMTB = TIMTB+TIMEM2-TIMEM1
   CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,MFILE,NPTS,NPANLS)
   CALL CPUTIM (TIMEM3)
   TIMOT = TIMOT+TIMEM3-TIMEM2
!
!     NPL IS NOW LOCATION OF FIRST ELEMENT TO BE USED FOR NEXT PASS
!
   IPL = -2
   DO 100 NL = NPL, NPE
      IPL = IPL+1
      TRALYR(IPL) = TRALYR(NL)
      RADLYR(IPL) = RADLYR(NL)
      RADLYB(IPL) = RADLYB(NL)
100 END DO
!
   V1P = V1P+ REAL(NPL+1)*DVP
   NPE = IPL
!
   GO TO 60
110 CONTINUE
!
   CALL CPUTIM (TIME1)
   TIM = TIME1-TIME
   WRITE (IPR,910) TIME1,TIM,TIMEM,TIMRD,TIMTB,TIMOT
!
   RETURN
!
900 FORMAT ('0 THE TIME AT THE START OF RADINT IS ',F12.3)
905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,   &
   &        ' MERGED WITH FILE ',I5,' ONTO FILE',I5,'  WITH XTYPE=',  &
   &        G15.5)
910 FORMAT ('0 THE TIME AT THE END OF RADINT IS ',F12.3/F12.3,        &
   &        ' SECS WERE REQUIRED FOR THIS MERGE  - EMIN - ',F12.3,    &
   &        ' - READ - ',F12.3,' - EMBND - ',F12.3,' - EMOUT - ',     &
   &        F12.3)
!
END SUBROUTINE RADINT
!
!     ----------------------------------------------------------------
!
SUBROUTINE EMBND (V1PO,V2PO,DVPO,NLIMO,NEWEM,NEWTR,TBND)
!
   USE phys_consts, ONLY: radcn2
   USE lblparams, ONLY: dbg, od_lo
   IMPLICIT REAL*8           (V)
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    9 APRIL 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
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
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   DIMENSION NEWEM(*),NEWTR(*)
!
   REAL NEWEM,NEWTR
!
   XKTBND = TBND/RADCN2
   VI = V1PO-DVPO
   VIDVBD = VI
   VIDVEM = VI
   BBLAST = -1.
   EMLAST = -1.
   NLIM1 = 0
   NLIM2 = 0
   EMDUM = 0.
   BBDUM = 0.
   EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMDUM)
!
   if (dbg(9)) then
      print *, 'EMBND :: NOT FIXED'
      dbg(9) = .false.
   endif
   VI = V1PO
10 NLIM1 = NLIM2+1
!
   EMISIV = EMISFN(VI,DVPO,VIDV,EMDEL,EMLAST)
!
   IF (VIDV.GE.9.E+4) THEN
      NLIM2 = NLIMO+1
   ELSE
      NLIM2 = (VIDV-V1PO)/DVPO+1.001
   ENDIF
   NLIM2 = MIN(NLIM2,NLIMO)
!
   DO J = NLIM1, NLIM2
      BB = PLANCK(VI,XKTBND)
      NEWEM(J) = NEWEM(J)+NEWTR(J)*EMISIV*BB
!
!        Increment interpolation values
!
      EMISIV = EMISIV+EMDEL
      VI = VI + DVPO
   END DO
!
   IF (NLIM2.LT.NLIMO) GO TO 10
!
   RETURN
!
END SUBROUTINE EMBND
!
!     ----------------------------------------------------------------
!
SUBROUTINE EMOUT (V1P,V2P,DVP,NLIM,NEWEM,NEWTR,MFILE,NPTS,NPANLS)
!
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE EMOUT OUTPUTS MERGED EMISSION AND TRANSMITTANCE RESULT
!     TO MFILE
!
   COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   DIMENSION PNLHDR(2)
   DIMENSION NEWEM(*),NEWTR(*)
!
   EQUIVALENCE (PNLHDR(1),V1PBF)
!
   REAL NEWEM,NEWTR
!
   NPANLS = NPANLS+1
   V1PBF = V1P
   V2PBF = V2P
   DVPBF = DVP
   NLIMBF = NLIM
!
   CALL BUFOUT (MFILE,PNLHDR(1),NPHDRF)
   CALL BUFOUT (MFILE,NEWEM(1),NLIMBF)
   CALL BUFOUT (MFILE,NEWTR(1),NLIMBF)
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
         WRITE (IPR,910) J,VJ,NEWEM(J),NEWTR(J), K,VK,NEWEM(K),NEWTR(&
            K)
10    CONTINUE
   ENDIF
!
   RETURN
!
900 FORMAT ('0 ','LOCATION  WAVENUMBER',2X,'RADIANCE',7X,             &
   &        'TRANSMITTANCE',22X,'LOCATION   WAVENUMBER',2X,           &
   &        'RADIANCE',7X,'TRANSMITTANCE')
905 FORMAT (' ')
910 FORMAT (I8,2X,F12.6,1P2E15.7,20X,I8,2X,0PF12.6,1P2E15.7)
!
END SUBROUTINE EMOUT
!
!     ----------------------------------------------------------------
!
SUBROUTINE EMDM (V1P,V2P,DVP,NLIM,KFILE,EM,EMB,TR,KEOF,NPANLS)
!
   USE lblparams, ONLY: NN_TBL, dbg, od_lo
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE EMDM INPUTS OPTICAL DEPTH VALUES FROM KFILE AND
!       CALCULATES SOURCE FUNCTION FOR THE LAYER.
!       THIS VERSION WORKS FOR AEROSOLS AND NLTE.
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    03 MARDCH 2006
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
!
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     131 Hartwell Ave., Lexington, MA. 02421
!
!----------------------------------------------------------------------
!
!               WORK SUPPORTED BY:    THE ARM PROGRAM
!                                     OFFICE OF ENERGY RESEARCH
!                                     DEPARTMENT OF ENERGY
!
!
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF
   COMMON /RMRG/ XKT,XKTA,XKTB,SECNT
   COMMON /EMDXSV/ BBEFF(2410),BBSAV(2410),BBASAV(2410),             &
   &               BBDSAV(2410),BBDLSAV(2410),fsav(2410),dF_dtau(2410)
!
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,      &
   &    pad_3,angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!
   DIMENSION PNLHDR(2),EM(*),EMB(*),TR(*)
!
   EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4)
   EQUIVALENCE (PNLHDR(1),V1PBF)
   EQUIVALENCE (FSCDID(4),IAERSL)
!
   data itbl_calc/-99/, aa /0.278/
!
   CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
   IF (KEOF.LE.0) RETURN
   CALL BUFIN (KFILE,KEOF,TR(1),NLIMBF)
!
!     TR contains the optical depths at this stage
!
   IF (IHIRAC.EQ.4) CALL BUFIN (KFILE,KEOF,EM(1),NLIMBF)
!
!     EM contains the optical depth corrections for nlte at this stage
!
   IF (NPANLS.LT.1) then
      if (IAERSL.EQ.0 .or. iaersl.eq.5) then
         WRITE (IPR,900)
      else
         WRITE (IPR,905)
      endif
   ENDIF
!
   EXT = 0.
   ADEL = 0.
   RADFN0 = 0.
   RDEL = 0.
   BB = 0.
   BBDEL = 0.
   BBA = 0.
   BBDLA = 0.
   BBB = 0.
   BBDLB = 0.
!
   V1P = V1PBF
   V2P = V2PBF
   DVP = DVPBF
   NLIM = NLIMBF
   VI = V1P-DVP
   VIDV = VI
   VIBB = VI
   VAER = VI
   VDUM = VI
   BBLAST = -1.
   BBLXTA = -2.
   BBLXTB = -3.
   RDLAST = -1.
   BBDUM = -4.
   RDDUM = -1.
   NLIM1 = 0
   NLIM2 = 0
!
   rec_6 = 1./6.
! **********************************************************************
!
   IF (IAERSL.EQ.0 .or. iaersl.eq.5) THEN
      IAFBB = -1
   ELSE
      RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDDUM)
      EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
      IF (VITST.LT.VAER) then
         IAFBB = 1
      ELSE
         IAFBB = 2
      ENDIF   
   ENDIF
!
!     - THIS SECTION TREATS THE CASE WHERE THE LAYER CONTRIBUTES
!       TO THE RADIATIVE TRANSFER ONLY ONCE
!
!     - WITH XKTA=0 THIS ALGORITHM REVERTS TO THE ORIGINAL
!
! **********************************************************************
!
   IF (XKTB.LE.0.) THEN
!
!     - THIS SECTION TREATS THE LTE CASE
!
      IF (IHIRAC.NE.4) THEN
!
         if (dbg(10)) then
            print *, 'EMDM::IHIRAC.NE.4 XKTB.LE.0. ::LTE:: CHECKED'
            dbg(10) = .false.
         ENDIF
         VI = V1P
         bb_dif=0.
10       NLIM1 = NLIM2+1
!
         IF (IAFBB.EQ.-1) THEN
            NLIM2 = NLIM
         ELSEIF (IAFBB.EQ.1) THEN
            RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
            EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
         ELSEIF (IAFBB.EQ.2) THEN
            EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            VITST = -VIDV
            RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
         ENDIF
!
         DO I = NLIM1, NLIM2
            BB = PLANCK(VI,XKT)
            IF (XKTA.GT.0.) THEN
               bb_dif = PLANCK(VI,XKTA)-BB
            ENDIF            
!
            ODVI = TR(I)+EXT*RADFN0
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
            if (odvi .lt. -od_lo) odvi = -od_lo
            tr_i   = exp(-odvi)
            if (odvi .le. od_lo) then
               f_i       = rec_6*odvi
               df_dtau_i = rec_6
               !
            else
               f_i        = 1. - 2.*(tr_i/(tr_i -1.)    + 1./odvi   )
               df_dtau_i  =     -2.*(tr_i/(tr_i -1.)**2 - 1./odvi**2)
               !
            end if
            TR(i)    = tr_i
            bbeff(i) = bb + bb_dif*f_i
            em(i)    = (1.-tr_i) * bbeff(i)
!
!---
! Store quantities for derivative calculation

            fsav(i)    = f_i
            dF_dtau(i) = df_dtau_i
            BBSAV(I)   = BB
            BBASAV(I)  = bb + bb_dif
!---

!              Increment interpolation values
!
            EXT = EXT+ADEL
            RADFN0 = RADFN0+RDEL
            VI = VI + DVP
         END DO
!
         IF (NLIM2.LT.NLIM) GO TO 10
      ELSE
!
!     - THIS SECTION TREATS THE NLTE CASE
!
         bb_dif=0.
         VI = V1P
30       NLIM1 = NLIM2+1
!
         VI = V1P+ REAL(NLIM1-1)*DVP
         IF (IAFBB.EQ.-1) THEN
            NLIM2 = NLIM
         ELSEIF (IAFBB.EQ.1) THEN
            RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
         ELSEIF (IAFBB.EQ.2) THEN
            EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            VITST = -VIDV
            RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
         ENDIF
!
         DO I = NLIM1, NLIM2
            BB = PLANCK(VI,XKT)
            IF (XKTA.GT.0.) THEN
               bb_dif = PLANCK(VI,XKTA)-BB
            ENDIF            
!              tr(i) contains the layer optical depths at this stage

            ODVI = TR(I)+EXT*RADFN0

!              em(i) contains the ratio differences from
!                                              lte of the state populati
            c_nlte = em(i)
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
            if (odvi .lt. -od_lo) odvi = -od_lo
            tr_i   = exp(-odvi)
            if (odvi .le. od_lo) then
               f_i       = rec_6*odvi
               df_dtau_i = rec_6
               !
               abs_i     = (odvi   - c_nlte) * (1.- 0.5 * odvi)
            else
               f_i        = 1. - 2.*(tr_i/(tr_i -1.)    + 1./odvi   )
               df_dtau_i  =     -2.*(tr_i/(tr_i -1.)**2 - 1./odvi**2)
               !
               abs_i  = (1. - c_nlte/odvi  ) * (1.-tr_i)
            end if
            TR(i)    = tr_i
            bbeff(i) = bb + bb_dif * f_i
            em(i)    = abs_i * bbeff(i)

            fsav(i)   = f_i
            dF_dtau(i)= df_dtau_i
            BBSAV(I)  = BB
            BBASAV(I) = bb + bb_dif
!---
!              Increment interpolation values
!
            EXT = EXT+ADEL
            RADFN0 = RADFN0+RDEL
            VI = VI + DVP
         END DO
!
         IF (NLIM2.LT.NLIM) GO TO 30
!
      ENDIF
! --------------------------------------------------------------
   ELSE
! --------------------------------------------------------------
!
!     - THIS SECTION TREATS THE CASE WHERE THE LAYER CONTRIBUTES
!       TO THE RADIATIVE TRANSFER TWICE:
!
!     - FOR TANGENT PATHS AND FOR THE CASE OF THE REFLECTED ATMOSPHERE
!
      IF (IHIRAC.NE.4) THEN
!
!     - THIS SECTION TREATS THE LTE CASE
!
         bb_dif_a = 0.
         bb_dif_b = 0.
         VI = V1P
50       NLIM1 = NLIM2+1
!
         VI = V1P+ REAL(NLIM1-1)*DVP
         IF (IAFBB.EQ.-1) THEN
            NLIM2 = NLIM
         ELSEIF (IAFBB.EQ.1) THEN
            RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
         ELSEIF (IAFBB.EQ.2) THEN
            EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            VITST = -VIDV
            RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
         ENDIF
!ccc
!     This calculation  is for specular reflection for the downwelling
!ccc

         if (surf_refl .eq. 's') then
            if (dbg(12)) then
               print *, 'EMDM::IHIRAC!=4 XKTB>0. ::LTE:: specular :: NOT CHECKED'
               dbg(12) = .false.
            ENDIF
!
            DO I = NLIM1, NLIM2
               BB = PLANCK(VI,XKT)
               IF (XKTA.GT.0.) THEN
                  bb_dif_a = PLANCK(VI,XKTA)-BB
               ENDIF            
               IF (XKTB.GT.0.) THEN
                  bb_dif_b = PLANCK(VI,XKTB)-BB
               ENDIF            
               ODVI = TR(I)+EXT*RADFN0
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
               if (odvi .lt. -od_lo) odvi = -od_lo
               tr_i   = exp(-odvi)
               if (odvi .le. od_lo) then
                  f_i       = rec_6*odvi
                  df_dtau_i = rec_6
                  !
               else
                  f_i        = 1. - 2.*(tr_i/(tr_i -1.)    + 1./odvi   )
                  df_dtau_i  =     -2.*(tr_i/(tr_i -1.)**2 - 1./odvi**2)
                  !
               end if
               TR(i)   = tr_i
               abs_i   = (1.-tr_i)
               bbeff(i)= bb + bb_dif_a * f_i
               em(i)   = abs_i * bbeff(i)
               emb(i)  = abs_i * (bb + bb_dif_b * f_i)
               fsav(i)    = f_i
               dF_dtau(i) = df_dtau_i
!---
! Store BB, BBA, and XX for derivative source term
               BBSAV(I) = BB
               BBASAV(I) = bb + bb_dif_a
!---
!
!     Increment interpolation values
!
               EXT = EXT+ADEL
               RADFN0 = RADFN0+RDEL
               VI = VI + DVP
            END DO
   !
         elseif (surf_refl .eq. 'l') then
            if (dbg(13)) then
               print *, 'EMDM::IHIRAC!=4 XKTB>0. ::LTE:: Lambertian :: NOT CHECKED'
               dbg(13) = .false.
            endif
!ccc
!     The following calculation is for an approximation to the
!     downwelling flux for application to Lambertian surfaces. The
!     'diffusivity' approximation is used with the assumption that the
!     dowmwelling flux is isotropic and that the surface scatters
!     isotropically.  with a value obtained from the
!     The value of the diffusivity angle corresponds to a secant of 1.67
!     diffuse_fac is the factor that is used to scale the optical depth
!     on the 'observer side' of the path to that required for the 'back
!ccc
            DO I = NLIM1, NLIM2
               BB = PLANCK(VI,XKT)
               IF (XKTA.GT.0.) THEN
                  bb_dif_a = PLANCK(VI,XKTA)-BB
               ENDIF            
               IF (XKTB.GT.0.) THEN
                  bb_dif_b = PLANCK(VI,XKTB)-BB
               ENDIF            
!
               ODVI = TR(I)+EXT*RADFN0
               odvi_d = diffuse_fac * odvi
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
               if (odvi .lt. -od_lo) odvi = -od_lo
               tr_i   = exp(-odvi)
               tr_d_i = exp(-odvi_d)
               if (odvi .le. od_lo) then
                  f_i       = rec_6*odvi
                  f_d_i     = rec_6*odvi_d
                  df_dtau_i = rec_6
                  !
               else
                  f_i        = 1. - 2.*(tr_i  /(tr_i   -1.)    + 1./odvi   )
                  f_d_i      = 1. - 2.*(tr_d_i/(tr_d_i -1.)    + 1./odvi_d )
                  df_dtau_i  =     -2.*(tr_i  /(tr_i   -1.)**2 - 1./odvi**2)
                  !
               end if
               TR(i)   = tr_i
               fsav(i) = f_i

               bbeff(i)= (bb + bb_dif_a * f_i)
               em(i)   = (1.-tr_i  ) * bbeff(i)
               emb(i)  = (1.-tr_d_i) * (bb + bb_dif_b * f_d_i) !mja, 10-27-2011
!---
! Store BB, BBA, and XX for derivative source term
               dF_dtau(i) = df_dtau_i
               BBSAV(I) = BB
               BBASAV(I) = bb + bb_dif_a
!---
!
!     Increment interpolation values
!
               EXT = EXT+ADEL
               RADFN0 = RADFN0+RDEL
               VI = VI + DVP
            END DO
   !
         else
            WRITE (IPR,906) surf_refl
            STOP 'EMDM:INVALID SURFACE REFLECTIVITY FLAG'
         endif
!
         IF (NLIM2.LT.NLIM) GO TO 50
!
      ELSE
!
!     - THIS SECTION TREATS THE CASE OF NLTE
!
         bb_dif_a = 0.
         bb_dif_b = 0.
         VI = V1P
70       NLIM1 = NLIM2+1
!
         IF (IAFBB.EQ.-1) THEN
            NLIM2 = NLIM
         ELSEIF (IAFBB.EQ.1) THEN
            RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
         ELSEIF (IAFBB.EQ.2) THEN
            EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            VITST = -VIDV
            RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
         ENDIF
!
!ccc
!     This calculation  is for specular reflection for the downwelling
!ccc

         if (surf_refl .eq. 's') then
            if (dbg(14)) then
               print *, 'EMDM::IHIRAC==4 XKTB>0. ::NLTE:: specular :: NOT CHECKED'
               dbg(14) = .false.
            END IF   
            DO I = NLIM1, NLIM2
               BB = PLANCK(VI,XKT)
               IF (XKTA.GT.0.) THEN
                  bb_dif_a = PLANCK(VI,XKTA)-BB
               ENDIF            
               IF (XKTB.GT.0.) THEN
                  bb_dif_b = PLANCK(VI,XKTB)-BB
               ENDIF            

!     tr(i) contains the layer optical depths at this stage

               ODVI = TR(I)+EXT*RADFN0
!
!     em(i) contains the ratio differences from lte of the state populat
!
               c_nlte = em(i)
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
               if (odvi .lt. -od_lo) odvi = -od_lo
               tr_i   = exp(-odvi)
               if (odvi .le. od_lo) then
                  f_i        = rec_6*odvi
                  df_dtau_i  = rec_6
                  !
                  abs_i   = (odvi   - c_nlte) * (1.- 0.5 * odvi)
               else
                  f_i        = 1. - 2.*(tr_i  /(tr_i   -1.)    + 1./odvi   )
                  df_dtau_i  =     -2.*(tr_i  /(tr_i   -1.)**2 - 1./odvi**2)
                  !
                  abs_i   = (1. - c_nlte/odvi  ) * (1.-tr_i)
               end if
               TR(i)     = tr_i
               bbeff(i)  = bb + bb_dif_a * f_i
               em(i)     = abs_i * bbeff(i)
               emb(i)    = abs_i * (bb + bb_dif_b * f_i)
               fsav(i)   = f_i
               dF_dtau(i)= df_dtau_i

!---
! Store BB, BBA, and XX for derivative source term
               BBSAV(I) = BB
               BBASAV(I) = bb + bb_dif_a
!---
!
!     Increment interpolation values
!
               EXT = EXT+ADEL
               RADFN0 = RADFN0+RDEL
               VI = VI + DVP
            END DO
   
         elseif (surf_refl .eq. 'l') then
            if (dbg(14)) then
               print *, 'EMDM::IHIRAC==4 XKTB>0. ::NLTE:: Lambertian :: NOT CHECKED'
               dbg(14) = .false.
            endif
!ccc
!     The following calculation is for an approximation to the
!     downwelling flux for application to Lambertian surfaces. The
!     'diffusivity' approximation is used with the assumption that the
!     dowmwelling flux is isotropic and that the surface scatters
!     isotropically.  with a value obtained from the
!     The value of the diffusivity angle corresponds to a secant of 1.67
!     diffuse_fac is the factor that is used to scale the optical depth
!     on the 'observer side' of the path to that required for the 'back
!ccc
            DO I = NLIM1, NLIM2
               BB = PLANCK(VI,XKT)
               IF (XKTA.GT.0.) THEN
                  bb_dif_a = PLANCK(VI,XKTA)-BB
               ENDIF            
               IF (XKTB.GT.0.) THEN
                  bb_dif_b = PLANCK(VI,XKTB)-BB
               ENDIF            

!     tr(i) contains the layer optical depths at this stage

               ODVI = TR(I)+EXT*RADFN0
               odvi_d = diffuse_fac * odvi
!
!     em(i) contains the ratio differences from  lte of the state popula
!
               c_nlte = em(i)
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
               if (odvi .lt. -od_lo) odvi = -od_lo
               tr_i   = exp(-odvi)
               tr_d_i = exp(-odvi_d)
               if (odvi .le. od_lo) then
                  f_i       = rec_6*odvi
                  f_d_i     = rec_6*odvi_d
                  df_dtau_i = rec_6
                  !
                  abs_i     = (odvi   - c_nlte) * (1.- 0.5 * odvi)
               else
                  f_i        = 1. - 2.*(tr_i  /(tr_i   -1.)    + 1./odvi   )
                  f_d_i      = 1. - 2.*(tr_d_i/(tr_d_i -1.)    + 1./odvi_d )
                  df_dtau_i  =     -2.*(tr_i  /(tr_i   -1.)**2 - 1./odvi**2)
                  !
                  abs_i   = (1.0 - c_nlte/odvi) * (1.-tr_i)
               end if
               TR(i)   = tr_i
               abs_d_i = (1.0 - c_nlte/odvi) * (1.-tr_d_i)
               bbeff(i)= bb + bb_dif_a * f_i

               em(i)   = abs_i   * bbeff(i)
               emb(i)  = abs_d_i * (bb + bb_dif_b * f_d_i)

               fsav(i)   = f_i
               dF_dtau(i)= dF_dtau_i

!     Store BB, BBA, and XX for derivative source term
               BBSAV(I) = BB
               BBASAV(I) = bb + bb_dif_a
!---
!
!     Increment interpolation values
!
               EXT = EXT+ADEL
               RADFN0 = RADFN0+RDEL
               VI = VI + DVP
            END DO
   
         else

            WRITE (IPR,906) surf_refl
            STOP 'EMDM::INVALID SURFACE REFLECTIVITY FLAG'

         endif
!
         IF (NLIM2.LT.NLIM) GO TO 70
!
      ENDIF
   ENDIF
!
   RETURN
!
900 FORMAT ('0EMISSION AND TRANSMISSION  (MOLECULAR) ')
905 FORMAT ('0EMISSION AND TRANSMISSION (AEROSOLS EFFECTS INCLUDED)')
906 FORMAT (' THE SURFACE REFLECTIVITY FLAG OF ', A1, 'IS NOT VALID')
907 FORMAT (' THE SURFACE REFLECTIVITY FLAG OF: ', A1)
!
END SUBROUTINE EMDM
!
!     ---------------------------------------------------------------
SUBROUTINE EMDT (V1P,V2P,DVP,NLIM,KFILE,EM,EMB,TR,KEOF,NPANLS)
!
   USE lblparams, ONLY: NN_TBL, dbg, od_lo
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE EMDT inputs optical depth values from kfile and
!       calculates source function for the layer.
!
!       This version is used for analytic temperature derivatives.
!       Non-LTE and limb not yet implemented.
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!                  IMPLEMENTATION:    P.D. Brown
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF
   COMMON /RMRG/ XKT,XKTA,XKTB,SECNT
   COMMON /EMDXSV/ BBEFF(2410),BBSAV(2410),BBASAV(2410),             &
   &               BBDSAV(2410),BBDLSAV(2410),fsav(2410),dF_dtau(2410)
!
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,      &
   &    pad_3,angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3

   data itbl_calc/-99/, aa /0.278/
!
   DIMENSION PNLHDR(2),EM(*),EMB(*),TR(*)
!
   EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4)
   EQUIVALENCE (PNLHDR(1),V1PBF)
   EQUIVALENCE (FSCDID(4),IAERSL)
!
   CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
   IF (KEOF.LE.0) RETURN
   CALL BUFIN (KFILE,KEOF,TR(1),NLIMBF)
!
!     TR contains the optical depths at this stage
!
   IF (IHIRAC.EQ.4) CALL BUFIN (KFILE,KEOF,EM(1),NLIMBF)
!
!     EM contains the optical depth corrections for nlte at this stage
!
   IF (NPANLS.LT.1) then
      if (IAERSL.EQ.0 .or. iaersl.eq.5) then
         WRITE (IPR,900)
      else
         WRITE (IPR,905)
      endif
   ENDIF
!
   EXT = 0.
   ADEL = 0.
   RADFN0 = 0.
   RDEL = 0.
   BB = 0.
   BBDEL = 0.
   BBA = 0.
   BBDLA = 0.
   BBB = 0.
   BBDLB = 0.

   BBD = 0.
   BBdTdel = 0.
   bbdTlast = -1.

   bbdT = 0.
   bbdT1= 0.
   bbdT2= -1.
!
   V1P = V1PBF
   V2P = V2PBF
   DVP = DVPBF
   NLIM = NLIMBF
   VI = V1P-DVP
   VIDV = VI
   VIBB = VI
   VAER = VI
   VDUM = VI
   BBLAST = -1.
   BBLXTA = -2.
   BBLXTB = -3.
   RDLAST = -1.
   BBDUM = -4.
   RDDUM = -1.
   NLIM1 = 0
   NLIM2 = 0
!
   ! vi
   vbbdT= v1p-dvp
   ! vi
   vidd = v1p-dvp

   rec_6 = 1./6.
! **********************************************************************
!
   IF (IAERSL.EQ.0 .or.iaersl.eq.5) THEN
      IAFBB = -1
   ELSE
      RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDDUM)
      EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
      IF (VITST.LT.VAER) then
         IAFBB = 1
      ELSE
         IAFBB = 2
      ENDIF   

   ENDIF
!
!     - THIS SECTION TREATS THE CASE WHERE THE LAYER CONTRIBUTES
!       TO THE RADIATIVE TRANSFER ONLY ONCE
!
!     - WITH XKTA=0 THIS ALGORITHM REVERTS TO THE ORIGINAL
!
! **********************************************************************
!
   IF (XKTB.LE.0.) THEN
!
!     - THIS SECTION TREATS THE LTE CASE
!
      IF (IHIRAC.NE.4) THEN
         if (dbg(15)) then
            print *, 'EMDT::IHIRAC!=4 XKTB<=0. :: LTE: CHECKED'
            dbg(15) = .false.
         endif
!
         bb_dif=0.
         VI = V1P
10       NLIM1 = NLIM2+1
!
         IF (IAFBB.EQ.-1) THEN
            NLIM2 = NLIM
         ELSEIF (IAFBB.EQ.1) THEN
            RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
         ELSEIF (IAFBB.EQ.2) THEN
            EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            VITST = -VIDV
            RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
         ENDIF
!
         DO I = NLIM1, NLIM2
            BB = PLANCK(VI,XKT)
            BBD = PLANCK_DT(VI,XKT,BB)
            IF (XKTA.GT.0.) THEN
               bb_dif = PLANCK(VI,XKTA)-BB
               bbdT = PLANCK_DT(VI,XKTA,BB)
            ENDIF            
!
            ODVI = TR(I)+EXT*RADFN0
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
            if (odvi .lt. -od_lo) odvi = -od_lo
            tr_i   = exp(-odvi)
            if (odvi .le. od_lo) then
               f_i        = rec_6*odvi
               df_dtau_i  = rec_6
            else
               f_i        = 1. - 2.*(tr_i  /(tr_i   -1.)    + 1./odvi   )
               df_dtau_i  =     -2.*(tr_i  /(tr_i   -1.)**2 - 1./odvi**2)
            end if
            TR(i)    = tr_i
            bbeff(i) = bb + bb_dif * f_i
            em(i)    = (1.-tr_i) * bbeff(i)
!---
! Store quantities for derivative calculation
!
            fsav(i) = f_i
            dF_dtau(i) = df_dtau_i
            bbdsav(i) = bbd
            bbdlsav(i) = bbdT
            BBSAV(I) = BB
            BBASAV(I) = bb + bb_dif
!---
!
!              Increment interpolation values
!
            EXT = EXT+ADEL
            RADFN0 = RADFN0+RDEL
            VI = VI + DVP
         END DO
!
         IF (NLIM2.LT.NLIM) GO TO 10
      ELSE
!
!     - THIS SECTION TREATS THE NLTE CASE
!
         if (dbg(16)) then
            print *, 'EMDT::IHIRAC==4 XKTB<=0. :: NLTE:: NOT CHECKED'
            dbg(16) = .false.
         endif
         bb_dif = 0.
         VI = V1P
30       NLIM1 = NLIM2+1
!
         IF (IAFBB.EQ.-1) THEN
            NLIM2 = NLIM
         ELSEIF (IAFBB.EQ.1) THEN
            RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
         ELSEIF (IAFBB.EQ.2) THEN
            EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            VITST = -VIDV
            RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
         ENDIF
!
         DO I = NLIM1, NLIM2
            BB = PLANCK(VI,XKT)
            IF (XKTA.GT.0.) THEN
               bb_dif = PLANCK(VI,XKTA)-BB
            ENDIF            

!              tr(i) contains the layer optical depths at this stage

            ODVI = TR(I)+EXT*RADFN0
!
!              em(i) contains the ratio differences from
!                                              lte of the state populati
            c_nlte = em(i)
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
            if (odvi .lt. -od_lo) odvi = -od_lo
            tr_i   = exp(-odvi)
            if (odvi .le. od_lo) then
               f_i        = rec_6*odvi
               df_dtau_i  = rec_6
               !
               abs_i   = (odvi   - c_nlte) * (1.- 0.5 * odvi)
            else
               f_i        = 1. - 2.*(tr_i  /(tr_i   -1.)    + 1./odvi   )
               df_dtau_i  =     -2.*(tr_i  /(tr_i   -1.)**2 - 1./odvi**2)
               !
               abs_i   = (1. - c_nlte/odvi  ) * (1.-tr_i)
            end if
            TR(i)    = tr_i
            bbeff(i) = bb + bb_dif*f_i
            em(i)    = abs_i * bbeff(i)
!---
! Store quantities for derivative calculation

            fsav(i) = f_i
            dF_dtau(i) = df_dtau_i
            BBSAV(I) = BB
            BBASAV(I) = bb + bb_dif
!---
!
!              Increment interpolation values
!
            EXT = EXT+ADEL
            RADFN0 = RADFN0+RDEL
            VI = VI + DVP
         END DO

!
         IF (NLIM2.LT.NLIM) GO TO 30
!
      ENDIF
! --------------------------------------------------------------
   ELSE
! --------------------------------------------------------------
!
!     - THIS SECTION TREATS THE CASE WHERE THE LAYER CONTRIBUTES
!       TO THE RADIATIVE TRANSFER TWICE:
!
!     - FOR TANGENT PATHS AND FOR THE CASE OF THE REFLECTED ATMOSPHERE
!
      IF (IHIRAC.NE.4) THEN
!
!     - THIS SECTION TREATS THE LTE CASE
!
         bb_dif_a = 0.
         bb_dif_b = 0.
         VI = V1P
50       NLIM1 = NLIM2+1
!
         IF (IAFBB.EQ.-1) THEN
            NLIM2 = NLIM
         ELSEIF (IAFBB.EQ.1) THEN
            RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
         ELSEIF (IAFBB.EQ.2) THEN
            EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            VITST = -VIDV
            RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
         ENDIF
!ccc
!     This calculation  is for specular reflection for the downwelling
!ccc

         if (surf_refl .eq. 's') then
            if (dbg(17)) then
               print *, 'EMDT::IHIRAC!=4 XKTB>0. :: LTE:: specular :: NOT CHECKED'
               dbg(17) = .false.
            endif
!
            DO I = NLIM1, NLIM2
               BB = PLANCK(VI,XKT)
               IF (XKTA.GT.0.) THEN
                  bb_dif_a = PLANCK(VI,XKTA)-BB
               ENDIF            
               IF (XKTB.GT.0.) THEN
                  bb_dif_b = PLANCK(VI,XKTB)-BB
               ENDIF            
!
               ODVI = TR(I)+EXT*RADFN0
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
               if (odvi .lt. -od_lo) odvi = -od_lo
               tr_i   = exp(-odvi)
               if (odvi .le. od_lo) then
                  f_i        = rec_6*odvi
                  df_dtau_i  = rec_6
               else
                  f_i        = 1. - 2.*(tr_i  /(tr_i   -1.)    + 1./odvi   )
                  df_dtau_i  =     -2.*(tr_i  /(tr_i   -1.)**2 - 1./odvi**2)
               end if
               TR(i)    = tr_i
               bbeff(i) = bb + bb_dif_a * f_i
               abs_i    = (1.-tr_i)
               em(i)    = abs_i * bbeff(i)
               emb(i)   = abs_i * (bb + bb_dif_b * f_i)
!---
! Store quantities for derivative calculation

               fsav(i)   = f_i
               dF_dtau(i)= df_dtau_i
               BBSAV(I)  = BB
               BBASAV(I) = bb + bb_dif_a
!---
!
!     Increment interpolation values
!
               EXT = EXT+ADEL
               RADFN0 = RADFN0+RDEL
               VI = VI + DVP
            END DO
   !
         elseif (surf_refl .eq. 'l') then
            if (dbg(18)) then
               print *, 'EMDT::IHIRAC!=4 XKTB>0. :: LTE:: Lambertian :: NOT CHECKED'
               dbg(18) = .false.
            endif
!ccc
!     The following calculation is for an approximation to the
!     downwelling flux for application to Lambertian surfaces. The
!     'diffusivity' approximation is used with the assumption that the
!     dowmwelling flux is isotropic and that the surface scatters
!     isotropically.  with a value obtained from the
!     The value of the diffusivity angle corresponds to a secant of 1.67
!     diffuse_fac is the factor that is used to scale the optical depth
!     on the 'observer side' of the path to that required for the 'back
!ccc
            DO I = NLIM1, NLIM2
               BB = PLANCK(VI,XKT)
               IF (XKTA.GT.0.) THEN
                  bb_dif_a = PLANCK(VI,XKTA)-BB
               ENDIF            
               IF (XKTB.GT.0.) THEN
                  bb_dif_b = PLANCK(VI,XKTB)-BB
               ENDIF            
!
               ODVI = TR(I)+EXT*RADFN0
               odvi_d = diffuse_fac * odvi
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
               if (odvi .lt. -od_lo) odvi = -od_lo
               tr_i   = exp(-odvi)
               tr_d_i = exp(-odvi_d)
               if (odvi .le. od_lo) then
                  f_i       = rec_6*odvi
                  f_d_i     = rec_6*odvi_d
                  df_dtau_i = rec_6
                  !
               else
                  f_i        = 1. - 2.*(tr_i  /(tr_i   -1.)    + 1./odvi   )
                  f_d_i      = 1. - 2.*(tr_d_i/(tr_d_i -1.)    + 1./odvi_d )
                  df_dtau_i  =     -2.*(tr_i  /(tr_i   -1.)**2 - 1./odvi**2)
                  !
               end if
               TR(i)   = tr_i

               bbeff(i) = bb + bb_dif_a * f_i
               em(i)    = (1.-tr_i) * bbeff(i)
               emb(i)   = (1.-tr_d_i) * (bb + bb_dif_b * f_d_i)
!---
! Store quantities for derivative calculation

               fsav(i)    = f_i
               dF_dtau(i) = df_dtau_i
               BBSAV(I)   = BB
               BBASAV(I)  = BB + bb_dif_a
!---
!
!     Increment interpolation values
!
               EXT = EXT+ADEL
               RADFN0 = RADFN0+RDEL
               VI = VI + DVP
            END DO
   !
         else
            WRITE (IPR,906) surf_refl
            STOP 'EMDT::INVALID SURFACE REFLECTIVITY FLAG'
         endif
!
         IF (NLIM2.LT.NLIM) GO TO 50
!
      ELSE
!
!     - THIS SECTION TREATS THE CASE OF NLTE
!
         bb_dif_a = 0.
         bb_dif_b = 0.
         VI = V1P
70       NLIM1 = NLIM2+1
!
         IF (IAFBB.EQ.-1) THEN
            NLIM2 = NLIM
         ELSEIF (IAFBB.EQ.1) THEN
            RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
         ELSEIF (IAFBB.EQ.2) THEN
            EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
            NLIM2 = (VIDV-V1P)/DVP+1.001
            NLIM2 = MIN(NLIM2,NLIM)
            VITST = -VIDV
            RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
         ENDIF
!
!ccc
!     This calculation  is for specular reflection for the downwelling
!ccc

         if (surf_refl .eq. 's') then
            if (dbg(19)) then
               print *, 'EMDT::IHIRAC==4 XKTB>0. :: NLTE:: specular :: NOT FIXED'
               dbg(19) = .false.
            endif
            DO I = NLIM1, NLIM2
               BB = PLANCK(VI,XKT)
               IF (XKTA.GT.0.) THEN
                  bb_dif_a = PLANCK(VI,XKTA)-BB
               ENDIF            
               IF (XKTB.GT.0.) THEN
                  bb_dif_b = PLANCK(VI,XKTB)-BB
               ENDIF            

!     tr(i) contains the layer optical depths at this stage

               ODVI = TR(I)+EXT*RADFN0
!
!     em(i) contains the ratio differences from lte of the state populat
!
               c_nlte = em(i)
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
               if (odvi .lt. -od_lo) odvi = -od_lo
               tr_i   = exp(-odvi)
               if (odvi .le. od_lo) then
                  f_i        = rec_6*odvi
                  df_dtau_i  = rec_6
                  !
                  abs_i   = (odvi   - c_nlte) * (1.- 0.5 * odvi)
               else
                  f_i        = 1. - 2.*(tr_i  /(tr_i   -1.)    + 1./odvi   )
                  df_dtau_i  =     -2.*(tr_i  /(tr_i   -1.)**2 - 1./odvi**2)
                  !
                  abs_i   = (1. - c_nlte/odvi  ) * (1.-tr_i)
               end if
               TR(i)    = tr_i
               bbeff(i) = bb + bb_dif_a * f_i
               em(i)    = abs_i * bbeff(i)
               emb(i)   = abs_i * (bb + bb_dif_b * f_i)

!---
! Other Derivative Terms

               fsav(i)    = f_i
               dF_dtau(i) = dF_dtau_i
               BBSAV(I)   = BB
               BBASAV(I)  = bb + bb_dif_a
!---
!
!     Increment interpolation values
!
               EXT = EXT+ADEL
               RADFN0 = RADFN0+RDEL
               VI = VI + DVP
            END DO
   
         elseif (surf_refl .eq. 'l') then
            if (dbg(20)) then
               print *, 'EMDT::IHIRAC==4 XKTB>0. :: NLTE:: Lambertian :: NOT FIXED'
               dbg(20) = .false.
            endif

!ccc
!     The following calculation is for an approximation to the
!     downwelling flux for application to Lambertian surfaces. The
!     'diffusivity' approximation is used with the assumption that the
!     dowmwelling flux is isotropic and that the surface scatters
!     isotropically.  with a value obtained from the
!     The value of the diffusivity angle corresponds to a secant of 1.67
!     diffuse_fac is the factor that is used to scale the optical depth
!     on the 'observer side' of the path to that required for the 'back
!ccc
            DO I = NLIM1, NLIM2
               BB = PLANCK(VI,XKT)
               IF (XKTA.GT.0.) THEN
                  bb_dif_a = PLANCK(VI,XKTA)-BB
               ENDIF            
               IF (XKTB.GT.0.) THEN
                  bb_dif_b = PLANCK(VI,XKTB)-BB
               ENDIF            

!     tr(i) contains the layer optical depths at this stage

               ODVI = TR(I)+EXT*RADFN0
               odvi_d = diffuse_fac * odvi
!
!     em(i) contains the ratio differences from  lte of the state popula
!
               c_nlte = em(i)
!
!       for odvi ouside the range of the table,  set optical depth to bo
!
               if (odvi .lt. -od_lo) odvi = -od_lo
               tr_i   = exp(-odvi)
               tr_d_i = exp(-odvi_d)
               if (odvi .le. od_lo) then
                  f_i        = rec_6*odvi
                  f_d_i      = rec_6*odvi_d
                  df_dtau_i  = rec_6
                  df_dtau_d_i= rec_6
                  !
                  abs_i   = (odvi   - c_nlte) * (1.- 0.5 * odvi)
                  abs_d_i = (odvi_d - c_nlte) * (1.- 0.5 * odvi_d)
               else
                  f_i        = 1. - 2.*(tr_i  /(tr_i   -1.)    + 1./odvi   )
                  f_d_i      = 1. - 2.*(tr_d_i/(tr_d_i -1.)    + 1./odvi_d )
                  df_dtau_i  =     -2.*(tr_i  /(tr_i   -1.)**2 - 1./odvi**2)
                  df_dtau_d_i=     -2.*(tr_d_i/(tr_d_i -1.)**2 - 1./odvi_d**2)
                  !
                  abs_i   = (1. - c_nlte/odvi  ) * (1.-tr_i)
                  abs_d_i = (1. - c_nlte/odvi_d) * (1.-tr_d_i)
               end if
               TR(i)   = tr_i
               bbeff(i)= bb + bb_dif_a * f_i
               em(i)   = abs_i   * bbeff(i)
               emb(i)  = abs_d_i * (bb + bb_dif_b * f_d_i)
!---
! Other Derivative Terms

               fsav(i)    = f_i
               dF_dtau(i) = dF_dtau_i
               BBSAV(I)   = BB
               BBASAV(I)  = bb + bb_dif_a
!---
!
!     Increment interpolation values
!
               EXT = EXT+ADEL
               RADFN0 = RADFN0+RDEL
               VI = VI + DVP
            END DO
   
         else

            WRITE (IPR,906) surf_refl
            STOP 'EMDT::INVALID SURFACE REFLECTIVITY FLAG'

         endif
!
         IF (NLIM2.LT.NLIM) GO TO 70
!
      ENDIF
   ENDIF
!
   RETURN
!
900 FORMAT ('0EMISSION AND TRANSMISSION  (MOLECULAR) ')
905 FORMAT ('0EMISSION AND TRANSMISSION (AEROSOLS EFFECTS INCLUDED)')
906 FORMAT (' THE SURFACE REFLECTIVITY FLAG OF ', A1, 'IS NOT VALID')
907 FORMAT (' THE SURFACE REFLECTIVITY FLAG OF: ', A1)
!
END SUBROUTINE EMDT
!
!     ---------------------------------------------------------------
!
SUBROUTINE EMADL1 (NPTS,MFILE,JPATHL,TBND)
!
   USE phys_consts, ONLY: radcn2
   USE lblparams, ONLY : dbg, od_lo
   IMPLICIT REAL*8           (V)
!
!     Calculates radiance and radiance derivative for first layer
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!                  IMPLEMENTATION:    P.D. Brown
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   COMMON NEWEM(2410),NEWTR(2410)
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYER,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /EMIHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYDUM,YI1,YID(10),LSTWDF
!
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,pad_3,&
   &    angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!
   COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /RMRG/ XKT,XKTA,XKTB,SECNT

   COMMON /ADRFIL/ KODFIL,kradtot,KTEMP,KFILAD,K_REFTRA,k_rddn_sfc
   COMMON /IADFLG/ NSPCRT,imrgsav
   common /DWNTRMS/ RFTRM(2410),raddwn(2410),tradwn(2410)
   common /pnldum/v1dwn,v2dwn,dvdwn,nlimdwn
   dimension xdwnin(2)
   equivalence (xdwnin(1),v1dwn)
   logical op

   CHARACTER*11 CFORM
   DATA CFORM / 'UNFORMATTED' /

!
   CHARACTER*40 CEXT,CYID
!
   DIMENSION EMLAYB(2410)
   DIMENSION XFILHD(2),OPNLHD(2),XFHDUM(2)
   DIMENSION EMLAYR(2),TRALYR(2)
!      DIMENSION XKMOLC(2),ODACUM(2)
   DIMENSION RPRIME(2410),TRAO(0:5000),OLDEM(0:5000)
!
   EQUIVALENCE (XFILHD(1),XID(1)) , (OPNLHD(1),V1PO)
   EQUIVALENCE (NEWEM(1),EMLAYR(1)) , (NEWTR(1),TRALYR(1)),          &
   &            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),               &
   &            (FSCDID(16),LAYR1)
!
   REAL NEWEM,NEWTR
!
!
!     TBND IS THE BOUNDARY BLACK BODY TEMPERATUE
!
!     IPATHL = -1 IS FOR downwelling radiance REFLECTED ATMOSPHERE
!     IPATHL =  0 IS FOR THE HORIZONTAL PATH CASE (HOMOGENEOUS LAYER)
!     IPATHL =  1 IS for upwelling radiance (TO DENSER LAYERS)
!     IPATHL =  2 IS FOR THE SYMMETRIC TANGENT PATH radiance
!     IPATHL =  3 IS for downwelling radiance (TO LESS DENSE LAYERS)
!     IPATHL = 31 is for downwelling radiance (to more dense layers)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   CALL CPUTIM (TIME)
!
!      ** NOTE ON IPATHL =2
!           THE TANGENT MERGE ROUTINES ARE DIVIDED INTO ANTERIOR (1ST)
!           AND POSTERIOR (2ND) LAYER CROSSINGS.  THIS RECOGNITION IS
!           TRIGGERED BY THE VALUE OF "IANT"
!
!          IF  IANT.EQ. 1  THEN POSTERIOR MERGE
!          IF  IANT.EQ. 0  THEN NORMAL MERGE
!          IF  IANT.EQ.-1  THEN ANTERIOR MERGE
!
   WRITE (IPR,900) TIME
   NPANLS = 0
!
!     Read in file headers for layer absorptance coefficients and
!     layer optical depths and total rad/trans (if there is more than
!     one layer between the present layer and the observer)
!
   ! dummy read to skip header
   CALL BUFIN (KFILE ,KEOF,XFHDUM(1),2)
   CALL BUFIN (KODFIL,KEOF,xfilhd(1),nfhdrf)
   IF (LAYER.LT.NLAYER) CALL BUFIN (kradtot,KEOF,XFHDUM(1),2)

   IF (JPATHL.GE.1) IPATHL = JPATHL
   PLAY = PAVE
   TLAY = TAVE
!
!     FOR AEROSOL RUNS, MOVE EXTID INTO YID
!
   IF (iaersl.ge.1 .and. iaersl.ne.5) THEN
      WRITE (CEXT,'(10A4)') EXTID
      WRITE (CYID,'(5A8)') (YID(I),I=3,7)
      CYID(19:40) = CEXT(19:40)
      READ (CYID,'(5A8)') (YID(I),I=3,7)
   ENDIF
!
   IEMIT = 1
   FACT = 1.
   TIMEM = 0.0
   IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.
   DO 10 MOL = 1, NMOL
      WK(MOL) = WK(MOL)*FACT
10 END DO
   WBROAD = WBROAD*FACT
   LAYR1 = LAYER
   WRITE (IPR,905) LAYR1,LAYER,KFILE,MFILE
   CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)
   CALL BUFOUT (KTEMP,XFILHD(1),NFHDRF)
   DVXM = DV
   XKT = TAVE/RADCN2
   XKTBND = TBND/RADCN2
   IF (IPATHL.EQ.-1) THEN
      XKTA = TZU/RADCN2
      XKTB = TZL/RADCN2
   ENDIF
   IF (IPATHL.EQ.0) THEN
      XKTA = 0.
      XKTB = 0.
   ENDIF
   IF (IPATHL.EQ.1) THEN
      XKTA = TZU/RADCN2
      XKTB = 0.
   ENDIF
   IF (IPATHL.EQ.2) THEN
      XKTA = TZU/RADCN2
      XKTB = TZL/RADCN2
   ENDIF
   IF (IPATHL.EQ.3) THEN
      XKTA = TZL/RADCN2
      XKTB = 0.
   ENDIF
   WRITE (IPR,910) IPATHL,IANT
!
20 CONTINUE
!
!
!     Input emission and transmission, and calculate layer
!     source function
!
!     Call EMDT for temperature derivative
!     Call EMDM for molecular species derivative
!
!
   CALL CPUTIM (TIMEM1)
!
   IF (NSPCRT.LE.0) THEN
      CALL EMDT (V1PO,V2PO,DVPO,NLIMO,KODFIL,EMLAYR,EMLAYB, TRALYR,  &
         KEOF,NPANLS)
   ELSEIF (NSPCRT.GT.0) THEN
      CALL EMDM (V1PO,V2PO,DVPO,NLIMO,KODFIL,EMLAYR,EMLAYB, TRALYR,  &
         KEOF,NPANLS)
   ENDIF
   CALL CPUTIM (TIMEM2)
   TIMEM = TIMEM+TIMEM2-TIMEM1
   IF (KEOF.LE.0) GO TO 80
   VI = V1PO-DVPO
   VIDVBD = VI
   VIDVEM = VI
   VIDVRF = VI
   BBLAST = -1.
   EMLAST = -1.
   RFLAST = -1.

! check to see if this is upwelling case (k_reftra file will be open)
! write panel size to sfc file, initialize emisout,reflout
   inquire(unit=k_reftra,opened=op)
   if (op) then
      call bufin(k_rddn_sfc,keof,xdwnin(1),nphdrf)
      call bufin(k_rddn_sfc,keof,raddwn(1),nlimdwn)
      call bufin(k_rddn_sfc,keof,tradwn(1),nlimdwn)
   endif

   IF (IPATHL.EQ.2.AND.IANT.EQ.0) THEN
      DO 30 J = 1, NLIMO
         TRJ = TRALYR(J)
         NEWEM(J) = EMLAYR(J)+EMLAYB(J)*TRJ
         TRALYR(J) = TRALYR(J)*TRJ
30    CONTINUE
!
   ELSEIF (IPATHL.EQ.3) THEN
!
      NLIM1 = 0
      NLIM2 = 0
      if (dbg(21)) then
         print *, 'EMADL1::(IPATHL.EQ.3): NOT CHECKED'
         dbg(21) = .false.
      endif
      EMDUM = 0.
      BBDUM = 0.
      EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMDUM)
!
      VI = V1PO
40    NLIM1 = NLIM2+1
!
      EMISIV = EMISFN(VI,DVPO,VIDV,EMDEL,EMLAST)
!
      IF (VIDV.GE.9.E+4) THEN
         NLIM2 = NLIMO+1
      ELSE
         NLIM2 = (VIDV-V1PO)/DVPO+1.001
      ENDIF
      NLIM2 = MIN(NLIM2,NLIMO)
!
      DO 50 J = NLIM1, NLIM2
         BB = PLANCK(VI,XKTBND)
         OLDEM(J) = BB*EMISIV
         NEWEM(J) = EMLAYR(J)+TRALYR(J)*OLDEM(J)
!
!           Increment interpolation values
!
         EMISIV = EMISIV+EMDEL
         VI = VI+ DVPO
50    CONTINUE
!
      IF (NLIM2.LT.NLIMO) GO TO 40
!
   ELSEIF (IPATHL.EQ.-1 .or. ipathl.eq.1) THEN
!
      if (dbg(31)) then
         print *, 'EMADL1::(IPATHL.EQ.-1 .or. ipathl.eq.1): CHECKED'
         dbg(31) = .false.
      endif
      NLIM1 = 0
      NLIM2 = 0
      EMDUM = 0.
      RFDUM = 0.
      BBDUM = 0.
      EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMDUM)
      REFLCT = REFLFN(VI,DVPO,VIDVRF,RFDEL,RFDUM)
      IF (VIDVEM.LE.VIDVRF) then
         IEMBB = 1
      ELSE
         IEMBB = 2
      ENDIF   
!
      VI = V1PO
60    NLIM1 = NLIM2+1
!
      IF (IEMBB.EQ.1) THEN
         EMISIV = EMISFN(VI,DVPO,VIDV,EMDEL,EMLAST)
         VIDVRF = -VIDV
         REFLCT = REFLFN(VI,DVPO,VIDVRF,RFDEL,RFLAST)
      ELSE
         REFLCT = REFLFN(VI,DVPO,VIDV,RFDEL,RFLAST)
         VIDVEM = -VIDV
         EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMLAST)
      ENDIF
!
      IF (VIDV.GE.9.E+4) THEN
         NLIM2 = NLIMO+1
      ELSE
         NLIM2 = (VIDV-V1PO)/DVPO+1.001
      ENDIF
      NLIM2 = MIN(NLIM2,NLIMO)
!
      DO 70 J = NLIM1, NLIM2
         BB = PLANCK(VI,XKTBND)
         trao (j) = 1.
         oldem(j) = BB*EMISIV+reflct*raddwn(j)
         NEWEM(J) = EMLAYR(J)+TRALYR(J)*oldem(j)
         rftrm(j) = reflct*tradwn(j)
!
!           Increment interpolation values
!
         BB = BB + BBDEL
         EMISIV = EMISIV+EMDEL
         REFLCT = REFLCT+RFDEL
         VI = VI + DVPO
70    CONTINUE
!
!     write out reflectance * total atmospheric transmittance for AJs
      if (op) then
         write(k_reftra) nlimo,(rftrm(j),j=1,nlimo)
      endif

      IF (NLIM2.LT.NLIMO) GO TO 60
!
   ENDIF
!
!     Output radiance to MFILE
!
   CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,NEWEM,NEWTR,MFILE,NPTS,NPANLS)
!
!     Combine terms of analytic layer radiance derivative
!
!     -------------------------------
!     Analytic Derivative calculation
!     -------------------------------
!
!
!     Call TDERIV for temperature retrieval
!     Call QDERIV for molecular retrieval
!
!     If molecular retrieval:
!     Input layer optical depth and accumulated
!     transmittances and calculate layer derivatives
!       (see lblrtm.f: default is d(log-x) so layer optical depth
!        is needed.  If d(x) is needed, then absorption coeffients
!        rather than optical depth will be used)
!
!     If temperature retrieval:
!     Input Planck function derivative and layer
!     transmittances and calculate layer derivatives
!
   IF     (NSPCRT.EQ.0 .and. ipathl.eq.1) THEN

      CALL TDERIVup (KFILE,kradtot,RPRIME,OLDEM,TRAO,TRALYR, NLIMO,  &
         IPATHL,LAYER,NLAYER,V1PO,DVPO)

   ELSEIF (NSPCRT.GT.0 .and. ipathl.eq.1) THEN

      CALL QDERIVup (KFILE,kradtot,RPRIME,OLDEM,TRAO,TRALYR, NLIMO,  &
         IPATHL,LAYER,NLAYER,V1PO,DVPO)

   ELSEIF (NSPCRT.EQ.0 .and. ipathl.eq.3) THEN

      CALL TDERIVdn (KFILE,kradtot,RPRIME,OLDEM,TRAO,TRALYR, NLIMO,  &
         IPATHL,LAYER,NLAYER,V1PO,DVPO)

   ELSEIF (NSPCRT.GT.0 .and. ipathl.eq.3) THEN

      CALL QDERIVdn (KFILE,kradtot,RPRIME,OLDEM,TRAO,TRALYR, NLIMO,  &
         IPATHL,LAYER,NLAYER,V1PO,DVPO)

   ENDIF
!
!     Output monochromatic radiance derivative to KTEMP
!
   CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,RPRIME,NEWTR,KTEMP,NPTS,NPANLS)
!
   GO TO 20
80 CALL CPUTIM (TIME1)
   TIME = TIME1-TIME
   WRITE (IPR,915) TIME,TIMEM
!
   RETURN
!
900 FORMAT (' TIME AT THE START OF --EMADL1-- ',F10.3)
905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,                &
   &        '0 INPUT FILE =',I5,' OUTPUT FILE =',I5)
910 FORMAT ('0 IPATHL AND IANT',2I5)
915 FORMAT (' TIME REQUIRED FOR --EMADL1-- ',F10.3,                   &
   &        ' --EMDM-- ',F10.3)
!
END SUBROUTINE EMADL1
!
!     ---------------------------------------------------------------
!
SUBROUTINE EMADMG (NPTS,LFILE,MFILE,JPATHL,TBND)
!
   USE phys_consts, ONLY: radcn2
   IMPLICIT REAL*8           (V)
!
!     Merges layer radiances when calculating layer
!     radiance derivatives.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   COMMON RADN(2410),TRAN(2410),RADO(0:5000)
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYER,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),      &
   &               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,  &
   &               EMISIV,FSCDID(17),NMOL,LAYDUM,YI1,YID(10),LSTWDF
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,      &
   &    pad_3,angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!
   COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
   COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /XME/ v1_pad,v2_pad,dv_pad,n_pad,TRAO(0:5000)
   COMMON /RMRG/ XKT,XKTA,XKTB,SECNT
!
   COMMON /EMDXSV/ BBEFF(2410),BBSAV(2410),BBASAV(2410),             &
   &               BBDSAV(2410),BBDLSAV(2410),fsav(2410),dF_dtau(2410)
   COMMON /ADRFIL/ KODFIL,kradtot,KTEMP,KFILAD,K_REFTRA,k_rddn_sfc
   COMMON /IADFLG/ NSPCRT,imrgsav
   common /DWNTRMS/ RFTRM(2410),raddwn(2410),tradwn(2410)
   common /pnldum/v1dwn,v2dwn,dvdwn,nlimdwn
   dimension xdwnin(2)
   equivalence (xdwnin(1),v1dwn)
   logical op
!
   DIMENSION RADLYB(2410)
   DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2),xfhdum(2)
   DIMENSION A1(10),A2(10),A3(10),A4(10)
   DIMENSION RADLYR(2),TRALYR(2),RADOI(2),TRAOI(2)
   DIMENSION WKSAV(35)
   DIMENSION RPRIME(2410)
!
   CHARACTER*40 CYID
!
   EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                 &
   &            (OPNLHD(1),V1PO)
   EQUIVALENCE (RADO(1),RADOI(1)) , (TRAO(1),TRAOI(1)),              &
   &            (RADN(1),RADLYR(1)) , (TRAN(1),TRALYR(1)),            &
   &            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),               &
   &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),               &
   &            (FSCDID(16),LAYR1)
!
   DATA A1 /10*0.0/, A2 /10*0.0/, A3 /10*0.0/, A4 /10*0.0/
!
!
!     IPATHL = -1 IS FOR downwelling radiance REFLECTED ATMOSPHERE
!     IPATHL =  0 IS FOR THE HORIZONTAL PATH CASE (HOMOGENEOUS LAYER)
!     IPATHL =  1 IS for upwelling radiance (TO DENSER LAYERS)
!     IPATHL =  2 IS FOR THE SYMMETRIC TANGENT PATH radiance
!     IPATHL =  3 IS for downwelling radiance (TO LESS DENSE LAYERS)
!     IPATHL = 31 is for downwelling radiance (to more dense layers)
!
!
!
!      ** NOTE ON IPATHL = 2
!            THE TANGENT MERGE ROUTINES ARE DIVIDED INTO ANTERIOR (1ST)
!            AND POSTERIOR (2ND) LAYER CROSSINGS   THIS RECOGNITION IS
!            TRIGGERED BY THE VALUE OF "IANT"
!
!          IF  IANT.EQ. 1  THEN POSTERIOR MERGE
!          IF  IANT.EQ. 0  THEN NORMAL MERGE
!          IF  IANT.EQ.-1  THEN ANTERIOR MERGE
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   CALL CPUTIM (TIME)
   WRITE (IPR,900) TIME
   NPANLS = 0
   TIMEM = 0.0
   TIMRD = 0.0
   TIMOT = 0.0
!
   CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)
   LAY1SV = LAYR1
   DVL = DV
   PL = PAVE
   TL = TAVE
   WTOTL = 0.
!
   DO 10 MOL = 1, NMOL
      WTOTL = WTOTL+WK(MOL)
      WKSAV(MOL) = WK(MOL)
10 END DO
!
   WTOTL = WTOTL+WBROAD
   WN2SAV = WBROAD
!
!     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)
!
   IF (iaersl.ge.1 .and. iaersl.ne.5)                                &
   &                 WRITE (CYID,'(5A8)') (YID(I),I=3,7)
!
!     Read in file headers for layer absorptance coefficients, layer
!     optical depths, and total optical depths (if there is more than
!     one layer between the present layer and the observer)
!
   ! dummy read to skip header
   CALL BUFIN (KFILE ,KEOF,XFHDUM(1),2)
   CALL BUFIN (KODFIL,KEOF,XFILHD(1),NFHDRF)
   IF (LAYER.LT.NLAYER) CALL BUFIN (kradtot,KEOF,XFHDUM(1),2)
!
   IF (iaersl.ge.1 .and. iaersl.ne.5)                                &
   &                 READ (CYID,'(5A8)') (YID(I),I=3,7)
!
   IF (JPATHL.GE.1) IPATHL = JPATHL
   PLAY = PAVE
   TLAY = TAVE
!
   TAVK = TAVE
   DVK = DV
   FACT = 1.
   IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.
!
   IF (DVL.EQ.DVK) THEN
      ITYPE = 0
   ELSEIF (DVL.GT.DVK) THEN
      ITYPE = DVK/(DVL-DVK)+0.5
   ELSE
!
!     DVL.LT.DVK
!
      ITYPE = -INT(DVL/(DVK-DVL)+0.5)
   ENDIF
   IF (ITYPE.LT.0) STOP ' EMADMG; ITYPE LT 0 '
!
   WTOTK = 0.
   LAYR1 = LAY1SV
   WRITE (IPR,905) LAYR1,LAYER,KODFIL,LFILE,MFILE
   IEMIT = 1
   DO 20 MOL = 1, NMOL
      WTOTK = WTOTK+WK(MOL)*FACT
      WK(MOL) = WK(MOL)*FACT+WKSAV(MOL)
20 END DO
   WTOTK = WTOTK+WBROAD*FACT
   WBROAD = WBROAD*FACT+WN2SAV
   PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)
   TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)
   SECANT = 0.
!
!     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES
!
   CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)
   CALL BUFOUT (KTEMP,XFILHD(1),NFHDRF)
   DVXM = DV
   XKT = TAVK/RADCN2
!
   WRITE (IPR,910) IPATHL,IANT
!
   IF (IPATHL.EQ.-1) THEN
      XKTA = TZU/RADCN2
      XKTB = TZL/RADCN2
   ELSEIF (IPATHL.EQ.1) THEN
      XKTA = TZU/RADCN2
      XKTB = 0.
   ELSEIF (IPATHL.EQ.2) THEN
      XKTA = TZU/RADCN2
      XKTB = TZL/RADCN2
   ELSEIF (IPATHL.EQ.3) THEN
      XKTA = TZL/RADCN2
      XKTB = 0.
   ELSE
      STOP ' EMADMG; IPATHL '
   ENDIF
!
   ATYPE = ITYPE
   LL = ITYPE+1
   AP = 1.0/(ATYPE+1.0)
!
!     *** BEGINNING OF LOOP THAT DOES MERGE ***
!
   NPE = 0
   RADO(0) = 0.0
   TRAO(0) = 0.0
   V1PO = 0.0
   V2PO = 0.0
   DVPO = 0.0
!
   inquire(unit=k_reftra,opened=op)
   if (op) rewind(k_reftra)

40 CONTINUE
!
!
!     Input emission and transmission, and calculate layer
!     source function
!
!     Call EMDT for temperature retrieval
!     Call EMDM for molecular retrieval
!
!
   CALL CPUTIM (TIMEM1)
!
   IF (NSPCRT.LE.0) THEN
      CALL EMDT (V1P,V2P,DVP,NLIM,KODFIL,RADLYR,RADLYB,TRALYR,KEOF,  &
         NPANLS)
   ELSEIF (NSPCRT.GT.0) THEN
      CALL EMDM (V1P,V2P,DVP,NLIM,KODFIL,RADLYR,RADLYB,TRALYR,KEOF,  &
         NPANLS)
   ENDIF
   CALL CPUTIM (TIMEM2)
   TIMEM = TIMEM+TIMEM2-TIMEM1
   IF (KEOF.LE.0) GO TO 80
!
   II = 1
!
   IF (V2PO.LE.V2P+DVPO) THEN
50    CALL CPUTIM (TIMEM1)
      CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)
      IF (LEOF.LE.0) GO TO 60
      CALL BUFIN (LFILE,LEOF,RADOI(NPE+1),NLIMO)
      CALL BUFIN (LFILE,LEOF,TRAOI(NPE+1),NLIMO)
      CALL CPUTIM (TIMEM2)
      TIMRD = TIMRD+TIMEM2-TIMEM1
      NPE = NLIMO+NPE
      IF (V2PO.LE.V2P+DVPO) GO TO 50
   ENDIF
!
!     ZERO POINT OF FIRST PANEL
!
60 IF (RADO(0).EQ.0.0.AND.TRAO(0).EQ.0.0) THEN
      RADO(0) = RADO(1)
      TRAO(0) = TRAO(1)
   ENDIF
!
!     END POINT OF LAST PANEL
!
   IF (V2PO+DVPO.GE.V2) THEN
      RADO(NPE+1) = RADO(NPE)
      RADO(NPE+2) = RADO(NPE)
      TRAO(NPE+1) = TRAO(NPE)
      TRAO(NPE+2) = TRAO(NPE)
   ENDIF


   inquire(unit=k_reftra,opened=op)
   if (op) then
      read(k_reftra) nlimsfc,(rftrm(j),j=1,nlimsfc)
      call bufin(k_rddn_sfc,keof,xdwnin(1),nphdrf)
      call bufin(k_rddn_sfc,keof,raddwn(1),nlimdwn)
      call bufin(k_rddn_sfc,keof,tradwn(1),nlimdwn)
   endif
!
!     -------------------------------
!     Analytic Derivative calculation
!     -------------------------------
!
!     Call TDERIV for temperature Jacobian
!     Call QDERIV for molecular Jacobian
!
!     If molecular retrieval:
!     Input layer optical depth and accumulated
!     transmittances and calculate layer derivatives
!       (see lblrtm.f: default is d(log-x) so layer optical depth
!        is needed.  If d(x) is needed, then absorption coeffients
!        rather than optical depth will be used)
!
!     If temperature retrieval:
!     Input Planck function derivative and layer
!     transmittances and calculate layer derivatives
!
   IF     (NSPCRT.EQ.0 .and. ipathl.eq.1) THEN

      CALL TDERIVup (KFILE,kradtot,RPRIME,RADO,TRAO,TRALYR, NLIM,    &
         IPATHL,LAYER,NLAYER,V1P,DVP)

   ELSEIF (NSPCRT.GT.0 .and. ipathl.eq.1) THEN

      CALL QDERIVup (KFILE,kradtot,RPRIME,RADO,TRAO,TRALYR, NLIM,    &
         IPATHL,LAYER,NLAYER,V1P,DVP)

   ELSEIF (NSPCRT.EQ.0 .and. ipathl.eq.3) THEN

      CALL TDERIVdn (KFILE,kradtot,RPRIME,RADO,TRAO,TRALYR, NLIM,    &
         IPATHL,LAYER,NLAYER,V1P,DVP)

   ELSEIF (NSPCRT.GT.0 .and. ipathl.eq.3) THEN

      CALL QDERIVdn (KFILE,kradtot,RPRIME,RADO,TRAO,TRALYR, NLIM,    &
         IPATHL,LAYER,NLAYER,V1P,DVP)

   ENDIF
!
!     Output monochromatic radiance derivatives to KTEMP
!
   CALL EMOUT (V1P,V2P,DVP,NLIM,RPRIME,TRAN,KTEMP,NPTS,NPANLS)
!
!     --------------------------
!     Layer Radiance Calculation
!     --------------------------
!
   NPL = 1
!
!     NPL IS LOCATION OF FIRST ELEMENT ON ARRAYS RADO AND TRAO
!     Combine terms of layer radiative transfer
!
   ipath_flg = ipathl
   if (ipathl .eq. -1) then
      if (surf_refl .eq. 's') ipath_flg = -10
      if (surf_refl .eq. 'l') ipath_flg = -11
   endif
!

   CALL RADNN (RADN,TRAN,RADO,TRAO,RADLYB,NLIM,V1P,DVP,              &
   &           IPATH_flg,A1,A2,A3,A4,LL,NPL)
!
   CALL CPUTIM (TIMEM1)
!
   IF (TBND.GT.0.) CALL EMBND (V1P,V2P,DVP,NLIM,RADN,TRAN,TBND)
!
!     Output radiance to MFILE
!
   CALL EMOUT (V1P,V2P,DVP,NLIM,RADN,TRAN,MFILE,NPTS,NPANLS)
!
   CALL CPUTIM (TIMEM2)
   TIMOT = TIMOT+TIMEM2-TIMEM1
!
!     NPL IS NOW LOCATION OF FIRST ELEMENT TO BE USED FOR NEXT PASS
!
   IPL = -1
   DO 70 NL = NPL, NPE
      IPL = IPL+1
      RADO(IPL) = RADO(NL)
      TRAO(IPL) = TRAO(NL)
70 END DO
!
   NPE = IPL
!
   GO TO 40
!
!     End of loop over panels
!
80 CONTINUE
!
   CALL CPUTIM (TIME1)
   TIM = TIME1-TIME
   WRITE (IPR,915) TIME1,TIM,TIMEM,TIMRD,TIMOT
!
   RETURN
!
!
900 FORMAT ('0 THE TIME AT THE START OF EMADMG IS ',F12.3)
905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,   &
   &        ' MERGED WITH FILE ',I5,' ONTO FILE',I5)
910 FORMAT ('0 IPATHL AND IANT',2I5)
915 FORMAT ('0 THE TIME AT THE END OF EMADMG IS ',F12.3/F12.3,        &
   &        ' SECS WERE REQUIRED FOR THIS MERGE  - EMDM - ',          &
   &        F12.3,' - READ - ',F12.3,' - EMOUT - ',F12.3)
!
END SUBROUTINE EMADMG
!
!     ---------------------------------------------------------------
!
SUBROUTINE QDERIVup (KFILE,kradtot,RPRIME,RADO,TRAO,TRALYR,       &
&                   NLIM,IPATHL,LAYER,NLAYER,V1PO,DVPO)
!
!     This subroutine inputs abosrption coefficient values from
!     KFILE and total transmittance from kradtot (if there is more
!     than one layer between the present layer and the observer),
!     and then calculates the radiance derivatives
!
   USE lblparams, ONLY: NDIM, ND2, IPTS, IPTS2, MXFSC, MXLAY, MXMOL
   IMPLICIT REAL*8           (V)
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
! if this changes, make sure it is changed in subroutines in xmerge.f
!      parameter (ndim=2410, nd2=5000)

   REAL KSUBL(0:ND2)
!
   DIMENSION RADO(0:ND2),TRAO(0:ND2),OPDT(0:ND2)
   DIMENSION RPRIME(NDIM),rdtotdn(ndim),trtotdn(ndim)
   DIMENSION TRALYR(*)
   DIMENSION PNLHDR(2),pnlhdq(2)
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
   COMMON /BUFPNL / V1P,V2P,DVP,NLIMP
   COMMON /BUFPNLq/ V1q,V2q,DVq,NLIMq
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /EMDXSV/ BBEFF(2410),BBSAV(2410),BBASAV(2410),             &
   &               BBDSAV(2410),BBDLSAV(2410),fsav(2410),dF_dtau(2410)
   common /DWNTRMS/ RFTRM(2410),raddwn(2410),tradwn(2410)

   EQUIVALENCE (PNLHDR(1),V1P), (PNLHDq(1),V1q)

   ! used for file check
   logical op


! note: from continuum module
!          ipts  = same dimension as ABSRB
!          ipts2 = same dimension as C
!      parameter (ipts=5050,ipts2=6000)
   common /CDERIV/ icflg,idum,v1absc,v2absc,dvabsc,nptabsc,delT_pert,&
   &    dqh2oC(ipts),dTh2oC(ipts),dUh2o

   real dtaudT(2400)

!********************
   character*20 h_radtot,h_kfile,h_k_od_molec

! for layer2level (if imoldq <> -99)
!      parameter (MXFSC=600, MXLAY=MXFSC+3, MXMOL=39)
   common /dlaydlev/ilevdx,imoldq,iupdwn,                            &
   &    dqdL(mxlay,0:mxmol),dqdU(mxlay,0:mxmol)

!---------------------------------------------------------------------

!       kfile       10   h_kfile     TAPE10                  ksubl
!       kradtot     18   h_radtot    RDDNlayer_00L           T_dn, R_dn
!       kfilad      19               AJ/RDderivUPW_00_001    layer deriv
!       kodfil      17               ODint_001               optical dep
!       ktemp       88               AJ_mono                 mono anal.
!       k_rddn_sfc  90               RDDNlayer_001           downwelling

!---------------------------------------------------------------------
!

! ksubl modified to include continuum term (if present)
! 1: H2O

   CALL BUFIN (kfile,KEOF,PNLHDR(1),NPHDRF)

   IF (KEOF.LE.0) THEN
      WRITE(*,*) 'End of KFILE ',KFILE
      stop 'qderivup'
   ENDIF
!
!     Read in absorptance coefficients
!
!     ksubl includes continuum term (if present)

   CALL BUFIN (kfile,KEOF,KSUBL(1),NLIMP)
!
!     Read in total  downwelling transmittance/radiance if LAYER < NLAYE
!
   IF (LAYER.LT.NLAYER) THEN

      CALL BUFIN (kradtot,KEOF,PNLHDq(1),NPHDRF)

      IF (KEOF.LE.0) THEN
         WRITE(IPR,900) kradtot,KFILE
         STOP 'IN SUBROUTINE QDERIVup: SEE OUTPUT FILE'
      ENDIF
!
!        Read in radiance and transmittance to current layer
!
      CALL BUFIN (kradtot,KEOF,rdtotdn(1),NLIMq)
      CALL BUFIN (kradtot,KEOF,trtotdn(1),NLIMq)
   ENDIF
!
!     Calculate layer derivatives,
!
!           RADO     = upwelling radiance into the layer
!           TRAO     = transmittance from surface to lower layer level
!           radtotdn = downwelling radiance into the layer
!           trtotdn  = (accumulated) total transmittance
!           TRALYR   = layer transmittance
!           KSUBL    = layer optical depth due to selected species
!           FSAV     = linear in tau fn
!           dF_dtau  = dF/dtau change in linear in tau function
!           BBSAV    = BBbar average layer Planck function
!           BBASAV   = BBa level A Planck function
!           BBEFF    = layer Emittance
!
!     When calculating the derivative of the layer nearest the observer,
!     omit the total accumulated transmittance, TRTOTDN(I)

10 CONTINUE
!
   IF (LAYER.lt.NLAYER) THEN

      DO 20 I = 1, NLIM

         betai = BBASAV(I)-BBSAV(I)

         rprime(i) = ksubl(i) * ( trtotdn(i) * ( (BBEFF(I)+betai*    &
            fsav(i)-RADO(I))*TRALYR(I) + (1.0-TRALYR(I))*(betai)*       &
            dF_dtau(I) ) + (BBEFF(I)-rdtotdn(i))*TRALYR(I)*trao(i) *    &
            rftrm(i) )

!            change in layer transmittance:

!            change in linear in tau term (F)

!            contribution from the downwelling radiance derivative refle
!
20    CONTINUE

   ELSE

      DO 30 I = 1, NLIM

         betai = BBASAV(I)-BBSAV(I)

         rprime(i) = ksubl(i) * ( ( (BBEFF(I)+betai*fsav(i)-RADO(I))*&
            TRALYR(I) + (1.0-TRALYR(I))*(betai)*dF_dtau(I) ) + (BBEFF(I)&
            )*TRALYR(I)*trao(i) * rftrm(i) )

!            change in layer transmittance:

!            change in linear in tau term (F)

!            contribution from the downwelling radiance derivative refle
!
30    CONTINUE

   ENDIF
!
   RETURN
!
900 FORMAT ('kradtot, ',I2.2,', reached end prior to end of KFILE, ', &
   &        I2.2)
!
END SUBROUTINE QDERIVup
!
!     ---------------------------------------------------------------
!
SUBROUTINE TDERIVup (KFILE,kradtot,RPRIME,RADO,TRAO,TRALYR,       &
&                   NLIM,IPATHL,LAYER,NLAYER,V1PO,DVPO)
!
!     This subroutine combines the Planck function derivative
!     (calculated in SUBROUTINE EMDT) and the layer transmittance
!     and then calculates the radiance derivatives with respect to
!     temperature
!
   USE lblparams, ONLY: NDIM, ND2, IPTS, IPTS2, MXFSC, MXLAY, MXMOL
   IMPLICIT REAL*8           (V)
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
! if this changes, make sure it is changed in subroutines in xmerge.f
!      parameter (ndim=2410, nd2=5000)

   REAL KSUBL(0:ND2)
!
   DIMENSION RADO(0:ND2),TRAO(0:ND2),OPDT(0:ND2)
   DIMENSION RPRIME(NDIM),rdtotdn(ndim),trtotdn(ndim)
   dimension rtmp(ndim)
   DIMENSION TRALYR(*)
   DIMENSION PNLHDR(2),pnlhdq(2)
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
   COMMON /BUFPNL/ V1P,V2P,DVP,NLIMP
   COMMON /BUFPNLq/ V1q,V2q,DVq,NLIMq
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /EMDXSV/ BBEFF(2410),BBSAV(2410),BBASAV(2410),             &
   &               BBDSAV(2410),BBDLSAV(2410),fsav(2410),dF_dtau(2410)
!
   common /DWNTRMS/ RFTRM(2410),raddwn(2410),tradwn(2410)

   EQUIVALENCE (PNLHDR(1),V1P), (PNLHDq(1),V1q)
!
   ! used for file check
   logical op

! note: from continuum module
!          ipts  = same dimension as ABSRB
!          ipts2 = same dimension as C
!      parameter (ipts=5050,ipts2=6000)
   common /CDERIV/ icflg,idum,v1absc,v2absc,dvabsc,nptabsc,delT_pert,&
   &    dqh2oC(ipts),dTh2oC(ipts),dUh2o

   real dtaudT(2400)

!********************
   character*20 h_radtot,h_kfile

! for layer2level (if imoldq <> -99)
!      parameter (MXFSC=600, MXLAY=MXFSC+3, MXMOL=39)
   common /dlaydlev/ilevdx,imoldq,iupdwn,                            &
   &    dqdL(mxlay,0:mxmol),dqdU(mxlay,0:mxmol)

!---------------------------------------------------------------------

!       kfile       10   h_kfile     TAPE10                  ksubl
!       kradtot     18   h_radtot    RDDNlayer_00L           T_dn, R_dn
!       kfilad      19               AJ/RDderivUPW_00_001    layer deriv
!       kodfil      17               ODint_001               optical dep
!       ktemp       88               AJ_mono                 mono anal.
!       k_rddn_sfc  90               RDDNlayer_001           downwelling

!---------------------------------------------------------------------

!
! ksubl file provides optical depths with temperature perturbation

   CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)

   IF (KEOF.LE.0) THEN
      WRITE(*,*) 'End of KFILE ',KFILE
      stop 'tderiv_up'
   ENDIF
!
   CALL BUFIN (KFILE,KEOF,KSUBL(1),NLIMP)
!
!     Read in total downwelling transmittance/radiance  if LAYER < NLAYE
!
   IF (LAYER.LT.NLAYER) THEN
      CALL BUFIN (kradtot,KEOF,PNLHDq(1),NPHDRF)
      IF (KEOF.LE.0) THEN
         WRITE(IPR,900) kradtot,KFILE
         STOP 'IN SUBROUTINE TDERIVup: SEE OUTPUT FILE'
      ENDIF
!
!        Read in radiance and transmattance to current layer
!
      CALL BUFIN (kradtot,KEOF,rdtotdn(1),NLIMq)
      CALL BUFIN (kradtot,KEOF,trtotdn(1),NLIMq)
   ENDIF
!
!     Calculate layer derivatives,
!
!  drsat/dT = (drdtau)(dtaudT)+(drdb)(dbdT) = srcnon + source
!           RADO     = upwelling radiance into the layer
!           TRAO     = transmittance from surface to lower layer level
!           radtotdn = downwelling radiance into the layer
!           trtotdn  = (accumulated) total transmittance
!           TRALYR   = layer transmittance
!           KSUBL    = layer optical depth with +1K temperature perturba
!           FSAV     = linear in tau fn
!           dF_dtau  = dF/dtau change in linear in tau function
!           BBSAV    = BBbar average layer Planck function
!           BBASAV   = BBa level A Planck function
!           BBEFF    = layer Emittance
!
10 CONTINUE

   IF (LAYER.lt.NLAYER) THEN

      DO 20 I = 1, NLIM

!     it would be better to have stored optical depths at this point!!!!

         if (tralyr(i) .gt. 1.e-06) then
            optdpt = -log(tralyr(i))
         else
            optdpt = ksubl(i)
         endif

!           ksubl is the optical depth with a +1K perturbation

         dtaudT(i) = ksubl(i) - optdpt
         betai = BBASAV(I)-BBSAV(I)

         rprime(i) = trtotdn(i) * ( (bbeff(i)+betai*fsav(i)-rado(i)) &
            * dtaudT(i) * tralyr(i) + bbdsav(i) * (1-tralyr(i)) + (1-   &
            tralyr(i)) * ((bbdlsav(i)-bbdsav(i))*fsav(i) + betai*       &
            dF_dtau(i)*dtaudT(i)) ) + (((bbeff(i)-rdtotdn(i)) * dtaudT( &
            i) * tralyr(i)) + (bbdsav(i) * (1-tralyr(i)))) * trao(i) *  &
            rftrm(i)

!           change in optical depth:

!           change in planck function

!           higher order terms

!            contribution from the downwelling radiance derivative refle

20    CONTINUE

   ELSE

!     this is the case for the top layer (layer .eq. nlayer)
!     When calculating the derivative of the layer nearest the observer,
!     omit the total accumulated transmittance, trtotdn, and term with
!     rdtotdn

      DO 30 I = 1, NLIM

         if (tralyr(i) .gt. 1.e-06) then
            optdpt = -log(tralyr(i))
         else
            optdpt = ksubl(i)
         endif

         dtaudT(i) = ksubl(i) - optdpt

         rprime(i) = ( (bbeff(i)+betai*fsav(i)-rado(i)) * dtaudT(i) *&
            tralyr(i) + bbdsav(i) * (1-tralyr(i)) + (1-tralyr(i)) * ((  &
            bbdlsav(i)-bbdsav(i))*fsav(i) + betai*dF_dtau(i)*dtaudT(i)) &
            ) + (((bbeff(i) ) * dtaudT(i) * tralyr(i)) + (bbdsav(i) * ( &
            1-tralyr(i)))) * trao(i) * rftrm(i)

!           change in optical depth:

!           change in planck function

!           higher order terms

!            ****  rdtotdn(i) is currently undefined for the top layer *
!            contribution from the downwelling radiance derivative refle

30    CONTINUE
   ENDIF
!
   RETURN
!
900 FORMAT ('kradtot, ',I2.2,', reached end prior to end of KFILE, ', &
   &        I2.2)
!
END SUBROUTINE TDERIVup
!
!     ---------------------------------------------------------------
!
SUBROUTINE QDERIVdn (KFILE,kradtot,RPRIME,RADO,TRAO,TRALYR,       &
&                   NLIM,IPATHL,LAYER,NLAYER,V1PO,DVPO)
!
!     This subroutine inputs abosrption coefficient values from
!     KFILE and total transmittance from kradtot (if there is more
!     than one layer between the present layer and the observer),
!     and then calculates the radiance derivatives
!
   USE lblparams, ONLY: NDIM, ND2, IPTS, IPTS2, MXFSC, MXLAY, MXMOL
   IMPLICIT REAL*8           (V)
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
! if this changes, make sure it is changed in subroutines in xmerge.f
!      parameter (ndim=2410, nd2=5000)

   REAL KSUBL(0:ND2)
!
   DIMENSION RADO(0:ND2),TRAO(0:ND2),OPDT(0:ND2)
   DIMENSION RPRIME(NDIM),rdtotdn(ndim),trtotdn(ndim)
   DIMENSION TRALYR(*)
   DIMENSION PNLHDR(2),pnlhdq(2)
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
   COMMON /BUFPNL / V1P,V2P,DVP,NLIMP
   COMMON /BUFPNLq/ V1q,V2q,DVq,NLIMq
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /EMDXSV/ BBEFF(2410),BBSAV(2410),BBASAV(2410),             &
   &               BBDSAV(2410),BBDLSAV(2410),fsav(2410),dF_dtau(2410)
   common /DWNTRMS/ RFTRM(2410),raddwn(2410),tradwn(2410)

   EQUIVALENCE (PNLHDR(1),V1P), (PNLHDq(1),V1q)

   ! used for file check
   logical op


! note: from continuum module
!          ipts  = same dimension as ABSRB
!          ipts2 = same dimension as C
!      parameter (ipts=5050,ipts2=6000)
   common /CDERIV/ icflg,idum,v1absc,v2absc,dvabsc,nptabsc,delT_pert,&
   &    dqh2oC(ipts),dTh2oC(ipts),dUh2o

   real dtaudT(2400)

!********************
   character*20 h_radtot,h_kfile,h_k_od_molec

! for layer2level (if imoldq <> -99)
!      parameter (MXFSC=600, MXLAY=MXFSC+3, MXMOL=39)
   common /dlaydlev/ilevdx,imoldq,iupdwn,                            &
   &    dqdL(mxlay,0:mxmol),dqdU(mxlay,0:mxmol)

!---------------------------------------------------------------------

!       kfile       10   h_kfile     TAPE10                  ksubl
!       kradtot     18   h_radtot    RDDNlayer_00L           T_dn, R_dn
!       kfilad      19               AJ/RDderivUPW_00_001    layer deriv
!       kodfil      17               ODint_001               optical dep
!       ktemp       88               AJ_mono                 mono anal.
!       k_rddn_sfc  90               RDDNlayer_001           downwelling

!---------------------------------------------------------------------
!

! ksubl modified to include continuum term (if present)
! 1: H2O

   CALL BUFIN (kfile,KEOF,PNLHDR(1),NPHDRF)

   IF (KEOF.LE.0) THEN
      WRITE(*,*) 'End of KFILE ',KFILE
      stop 'qderivdn'
   ENDIF
!
!     Read in absorptance coefficients
!
!     ksubl includes continuum term (if present)

   CALL BUFIN (kfile,KEOF,KSUBL(1),NLIMP)
!
!     Read in total  downwelling transmittance/radiance if LAYER < NLAYE
!
   IF (LAYER.LT.NLAYER) THEN

      CALL BUFIN (kradtot,KEOF,PNLHDq(1),NPHDRF)

      IF (KEOF.LE.0) THEN
         WRITE(IPR,900) kradtot,KFILE
         STOP 'IN SUBROUTINE QDERIVup: SEE OUTPUT FILE'
      ENDIF
!
!        Read in radiance and transmittance to current layer
!
      CALL BUFIN (kradtot,KEOF,rdtotdn(1),NLIMq)
      CALL BUFIN (kradtot,KEOF,trtotdn(1),NLIMq)
   ENDIF
!
!     Calculate layer derivatives,
!
!           RADO     = upwelling radiance into the layer
!           TRAO     = transmittance from surface to lower layer level
!           radtotdn = downwelling radiance into the layer
!           trtotdn  = (accumulated) total transmittance
!           TRALYR   = layer transmittance
!           KSUBL    = layer optical depth with +1K temperature perturba
!           FSAV     = linear in tau fn
!           dF_dtau  = dF/dtau change in linear in tau function
!           BBSAV    = BBbar average layer Planck function
!           BBASAV   = BBa level A Planck function
!           BBEFF    = layer Emittance
!
!     When calculating the derivative of the layer nearest the observer,
!     omit the total accumulated transmittance, TRAO(I)

10 CONTINUE

!
   IF (LAYER.eq.1) THEN

      DO 20 I = 1, NLIM

         betai = BBASAV(I)-BBSAV(I)

         rprime(i) = ksubl(i) * ( (BBEFF(I)+betai*fsav(i)-rdtotdn(i))&
            *TRALYR(I) + (1.0-TRALYR(I))*(betai)*dF_dtau(I) )

!            change in layer transmittance:

!            change in linear in tau term (F)
!
20    CONTINUE

   ELSE

      DO 30 I = 1, NLIM

         betai = BBASAV(I)-BBSAV(I)

         rprime(i) = trao(i) * ksubl(i) * ( (BBEFF(I)+betai*fsav(i)- &
            rdtotdn(i))*TRALYR(I) + (1.0-TRALYR(I))*betai*dF_dtau(I)    &
            )

!            change in layer transmittance:

!            change in linear in tau term (F)
!
30    CONTINUE

   ENDIF
!
   RETURN
!
900 FORMAT ('kradtot, ',I2.2,', reached end prior to end of KFILE, ', &
   &        I2.2)
!
END SUBROUTINE QDERIVdn
!
!     ---------------------------------------------------------------
!
SUBROUTINE TDERIVdn (KFILE,kradtot,RPRIME,RADO,TRAO,TRALYR,       &
&                   NLIM,IPATHL,LAYER,NLAYER,V1PO,DVPO)
!
!     This subroutine combines the Planck function derivative
!     (calculated in SUBROUTINE EMDT) and the layer transmittance
!     and then calculates the radiance derivatives with respect to
!     temperature
!
   USE lblparams, ONLY: NDIM, ND2, IPTS, IPTS2, MXFSC, MXLAY, MXMOL
   IMPLICIT REAL*8           (V)
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
! if this changes, make sure it is changed in subroutines in xmerge.f
!      parameter (ndim=2410, nd2=5000)

   REAL KSUBL(0:ND2)
!
   DIMENSION RADO(0:ND2),TRAO(0:ND2),OPDT(0:ND2)
   DIMENSION RPRIME(NDIM),rdtotdn(ndim),trtotdn(ndim)
   dimension rtmp(ndim)
   DIMENSION TRALYR(*)
   DIMENSION PNLHDR(2),pnlhdq(2)
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
   COMMON /BUFPNL/ V1P,V2P,DVP,NLIMP
   COMMON /BUFPNLq/ V1q,V2q,DVq,NLIMq
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /EMDXSV/ BBEFF(2410),BBSAV(2410),BBASAV(2410),             &
   &               BBDSAV(2410),BBDLSAV(2410),fsav(2410),dF_dtau(2410)
!
   common /DWNTRMS/ RFTRM(2410),raddwn(2410),tradwn(2410)

   EQUIVALENCE (PNLHDR(1),V1P), (PNLHDq(1),V1q)
!
   ! used for file check
   logical op

! note: from continuum module
!          ipts  = same dimension as ABSRB
!          ipts2 = same dimension as C
!      parameter (ipts=5050,ipts2=6000)
   common /CDERIV/ icflg,idum,v1absc,v2absc,dvabsc,nptabsc,delT_pert,&
   &    dqh2oC(ipts),dTh2oC(ipts),dUh2o

   real dtaudT(2400)

!********************
   character*20 h_radtot,h_kfile

! for layer2level (if imoldq <> -99)
!      parameter (MXFSC=600, MXLAY=MXFSC+3, MXMOL=39)
   common /dlaydlev/ilevdx,imoldq,iupdwn,                            &
   &    dqdL(mxlay,0:mxmol),dqdU(mxlay,0:mxmol)

!---------------------------------------------------------------------

!       kfile       10   h_kfile     TAPE10                  ksubl
!       kradtot     18   h_radtot    RDDNlayer_00L           T_dn, R_dn
!       kfilad      19               AJ/RDderivUPW_00_001    layer deriv
!       kodfil      17               ODint_001               optical dep
!       ktemp       88               AJ_mono                 mono anal.
!       k_rddn_sfc  90               RDDNlayer_001           downwelling

!---------------------------------------------------------------------

!
! ksubl file provides optical depths with temperature perturbation

   CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)

   IF (KEOF.LE.0) THEN
      WRITE(*,*) 'End of KFILE ',KFILE
      stop 'tderiv_dn'
   ENDIF
!
   CALL BUFIN (KFILE,KEOF,KSUBL(1),NLIMP)
!
!     Read in total downwelling transmittance/radiance  if LAYER < NLAYE
!
   IF (LAYER.LT.NLAYER) THEN
      CALL BUFIN (kradtot,KEOF,PNLHDq(1),NPHDRF)
      IF (KEOF.LE.0) THEN
         WRITE(IPR,900) kradtot,KFILE
         STOP 'IN SUBROUTINE TDERIVup: SEE OUTPUT FILE'
      ENDIF
!
!        Read in radiance and transmattance to current layer
!
      CALL BUFIN (kradtot,KEOF,rdtotdn(1),NLIMq)
      CALL BUFIN (kradtot,KEOF,trtotdn(1),NLIMq)
   ENDIF
!
!     Calculate layer derivatives,
!
!  drsat/dT = (drdtau)(dtaudT)+(drdb)(dbdT) = srcnon + source
!           RADO     = upwelling radiance into the layer
!           TRAO     = transmittance from surface to lower layer level
!           radtotdn = downwelling radiance into the layer
!           trtotdn  = (accumulated) total transmittance
!           TRALYR   = layer transmittance
!           KSUBL    = layer optical depth with +1K temperature perturba
!           FSAV     = linear in tau fn
!           dF_dtau  = dF/dtau change in linear in tau function
!           BBSAV    = BBbar average layer Planck function
!           BBASAV   = BBa level A Planck function
!           BBEFF    = layer Emittance
!
10 CONTINUE

   IF (LAYER.eq.1) THEN

      DO 20 I = 1, NLIM

!     it would be better to have stored optical depths at this point!!!!

         if (tralyr(i) .gt. 1.e-06) then
            optdpt = -log(tralyr(i))
         else
            optdpt = ksubl(i)
         endif

!           ksubl is the optical depth with a +1K perturbation

         dtaudT(i) = ksubl(i) - optdpt
         betai = BBASAV(I)-BBSAV(I)

         rprime(i) =                                                    &

!           change in optical depth:
         &     ( (bbeff(i)+betai*fsav(i)-rdtotdn(i)) * dtaudT(i) * tralyr(i)   &

!           change in planck function
         &     +  bbdsav(i) * (1-tralyr(i))                                    &

!           higher order terms
         &     + (1-tralyr(i)) * ((bbdlsav(i)-bbdsav(i))*fsav(i)               &
         &              + betai*dF_dtau(i)*dtaudT(i))   )

20    CONTINUE

   ELSE

!     this is the case for the top layer (layer .eq. nlayer)
!     When calculating the derivative of the layer nearest the observer,
!     omit the total accumulated transmittance, trtotdn, and term with
!     rdtotdn

      DO 30 I = 1, NLIM

         if (tralyr(i) .gt. 1.e-06) then
            optdpt = -log(tralyr(i))
         else
            optdpt = ksubl(i)
         endif

         dtaudT(i) = ksubl(i) - optdpt
         betai = BBASAV(I)-BBSAV(I)

         rprime(i) = trao(i) *                                         &

!           change in optical depth:
         &     ( (bbeff(i)+betai*fsav(i)-rdtotdn(i)) * dtaudT(i) * tralyr(i)  &

!           change in planck function
         &     +  bbdsav(i) * (1-tralyr(i))                                   &

!           higher order terms
         &     + (1-tralyr(i)) * ((bbdlsav(i)-bbdsav(i))*fsav(i)              &
         &           + betai*dF_dtau(i)*dtaudT(i))   )
!
30    CONTINUE

   ENDIF
!
   RETURN
!
900 FORMAT ('kradtot, ',I2.2,', reached end prior to end of KFILE, ', &
   &        I2.2)
!
END SUBROUTINE TDERIVdn
!
!     ----------------------------------------------------------------
!
SUBROUTINE FLXIN (V1P,V2P,DVP,NLIM,KFILE,EM,TR,KEOF,NPANLS)
!
   USE lblparams, ONLY: NN_TBL, dbg, od_lo
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE FLXIN INPUTS OPTICAL DEPTH VALUES FROM KFILE AND
!     CALCULATES FLUX FOR THE LAYER. THIS VERSION WORKS FOR AEROSOLS.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    14 AUGUST 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     M.J. IACONO
!                                     R.D. WORSHAM
!                                     J.L. MONCET
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
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   &
   &                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         &
   &              NLNGTH,KDUMY,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,        &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF
   COMMON /RMRG/ XKT,XKTA,XKTB,SECNT
!
   DIMENSION PNLHDR(2),EM(*),TR(*)
!
   EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4)
   EQUIVALENCE (PNLHDR(1),V1PBF)
   EQUIVALENCE (FSCDID(4),IAERSL)
!
   data itbl_calc/-99/, aa /0.278/
!
   CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
   IF (KEOF.LE.0) RETURN
   CALL BUFIN (KFILE,KEOF,TR(1),NLIMBF)
!
!     TR CONTAINS THE OPTICAL DEPTHS AT THIS STAGE
!
   IF (IHIRAC.EQ.4) STOP ' IHIRAC=4  FLXIN '
!
!     EM CONTAINS THE OPTICAL DEPTH CORRECTIONS FOR NLTE AT THIS STAGE
!
   IF (NPANLS.LT.1) then
      if (IAERSL.EQ.0 .or. iaersl.eq.5)  then
         WRITE (IPR,900)
      else
         WRITE (IPR,905)
      endif
   ENDIF
!
   EXT = 0.
   ADEL = 0.
   RADFN0 = 0.
   RDEL = 0.
   BB = 0.
   BBDEL = 0.
   BBA = 0.
   BBDLA = 0.
!
   V1P = V1PBF
   V2P = V2PBF
   DVP = DVPBF
   NLIM = NLIMBF
   VI = V1P-DVP
   VIDV = VI
   VIBB = VI
   VAER = VI
   VDUM = VI
   BBLAST = -1.
   BBLXTA = -2.
   RDLAST = -1.
   BBDUM = -4.
   RDDUM = -1.
   NLIM1 = 0
   NLIM2 = 0
!
   rec_6 = 1./6.
!
   IF (iaersl.ge.1 .and. iaersl.ne.5) THEN
      RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDDUM)
      EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
      IF (VITST.LT.VAER) then
         IAFBB = 1
      ELSE
         IAFBB = 2
      ENDIF   
   ELSE
      IAFBB = -1
   ENDIF
!
!     - THIS SECTION TREATS THE CASE WHERE THE LAYER CONTRIBUTES
!       TO THE RADIATIVE TRANSFER ONLY ONCE
!
!     - WITH XKTA=0 THIS ALGORITHM REVERTS TO THE ORIGINAL
!
   IF (XKTB.GT.0.) STOP ' XKTB GT 0.   FLXIN '
   if (dbg(23)) then
      print *, 'FLXIN:: NOT CHECKED'
      dbg(23) = .false.
   endif
!
   bb_dif = 0.
   VI = V1P
10 NLIM1 = NLIM2+1
!
   VI = V1P+ REAL(NLIM1-1)*DVP
   IF (IAFBB.EQ.-1) THEN
      NLIM2 = NLIM
   ELSEIF (IAFBB.EQ.1) THEN
      RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)
      NLIM2 = (VIDV-V1P)/DVP+1.001
      NLIM2 = MIN(NLIM2,NLIM)
      VAER = -VIDV
      EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
   ELSEIF (IAFBB.EQ.2) THEN
      EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
      NLIM2 = (VIDV-V1P)/DVP+1.001
      NLIM2 = MIN(NLIM2,NLIM)
      VITST = -VIDV
      RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)
   ENDIF
!
   if (dbg(26)) then
      print *, 'FLXIN:: NOT CHECKED'
      dbg(26) = .false.
   endif

   DO 20 I = NLIM1, NLIM2
      BB = PLANCK(VI,XKT)
      IF (XKTA.GT.0.) THEN
         bb_dif = PLANCK(VI,XKTA)-BB
      ENDIF            

      ODVI = SECNT*TR(I)+EXT*RADFN0
!
!       for odvi ouside the range of the table,  set optical depth to bound
!
      if (odvi .lt. -od_lo)  odvi = -od_lo
      tr_i   = exp(-odvi)
      if (odvi .le. od_lo) then                                     ! analytic regime
         f_i        = rec_6*odvi
      else                                                          ! use tables
         f_i        = 1. - 2.*(tr_i  /(tr_i   -1.)    + 1./odvi   )
      end if
      em(i) = (1.-tr_i) * (bb + bb_dif*f_i)
      TR(i)   = tr_i
!
!              Increment interpolation values
!
      EXT = EXT+ADEL
      RADFN0 = RADFN0+RDEL
      VI = VI+ DVP
      !
20 CONTINUE
!
   IF (NLIM2.LT.NLIM) GO TO 10
!
   RETURN
!
900 FORMAT ('0EMISSION AND TRANSMISSION  (MOLECULAR) ')
905 FORMAT ('0EMISSION AND TRANSMISSION (AEROSOLS EFFECTS INCLUDED)')
!
END SUBROUTINE FLXIN
!
!     ----------------------------------------------------------------
!
SUBROUTINE FLINIT (NPTS,MFILE,JPATHL,TBND,refl_flg)
!
   USE phys_consts, ONLY: radcn2
   USE lblparams, ONLY : dbg, od_lo
   IMPLICIT REAL*8           (V)
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    14 AUGUST 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     M.J. IACONO
!                                     R.D. WORSHAM
!                                     J.L. MONCET
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
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   COMMON NEWEM(2410),NEWTR(2410)
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   &
   &              AVFIX,LAYER,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,        &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /EMIHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   &
   &                EMISIV,FSCDID(17),NMOL,LAYDUM,YI1,YID(10),LSTWDF
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,        &
   &    pad_3,angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!
   COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,        &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /RMRG/ XKT,XKTA,XKTB,SECNT
!
   CHARACTER*40 CEXT,CYID
!
   integer refl_flg
   real emdown(5000)
   integer ifiledown

   REAL NEWEM,NEWTR
!
   DIMENSION XFILHD(2),OPNLHD(2)
   DIMENSION EMLAYR(2),TRALYR(2)
!
   EQUIVALENCE (XFILHD(1),XID(1)) , (OPNLHD(1),V1PO)
   EQUIVALENCE (NEWEM(1),EMLAYR(1)) , (NEWTR(1),TRALYR(1)),            &
   &            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),                 &
   &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 &
   &            (FSCDID(16),LAYR1)
!
!
!   *******************************************************************
!   ***  THIS SUBROUTINE COMPUTES THE EMISSION FOR THE FIRST LAYER  ***
!   *******************************************************************
!
!     TBND IS THE BOUNDARY BLACK BODY TEMPERATUE
!
!     IPATHL = 1 IS FOR THE DOWNWELLING FLUX CASE
!     IPATHL = 3 IS FOR THE UPWELLING FLUX CASE
!
   CALL CPUTIM (TIME)
!
   WRITE (IPR,900) TIME
   NPANLS = 0
   CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)
   IF (JPATHL.GE.1) IPATHL = JPATHL
   PLAY = PAVE
   TLAY = TAVE
!
!     FOR AEROSOL RUNS, MOVE EXTID INTO YID
!
   IF (iaersl.ge.1 .and. iaersl.ne.5) THEN
      WRITE (CEXT,'(10A4)') EXTID
      WRITE (CYID,'(5A8)') (YID(I),I=3,7)
      CYID(19:40) = CEXT(19:40)
      READ (CYID,'(5A8)') (YID(I),I=3,7)
   ENDIF
!
!     READ IN DIRECTION COSINE
!
   READ (IRD,905) DIRCOS
   SECNT = 1.0/DIRCOS
   SECANT = SECNT
   SECNT0 = SECNT
   WRITE (IPR,910) DIRCOS

   if (dircos .gt. 0.6) then
      ifiledown = 601
   elseif (dircos .gt. 0.3) then
      ifiledown = 602
   else
      ifiledown = 603
   endif
!
!     IF BOUNDARY PROPERTIES ARE SUPPLIED, AND DOWNWARD LOOKING
!     CASE; SET IPATHL TO REFLECTED ATMOSPHERE CASE
!
   IF (IBPROP.EQ.1.AND.IPATHL.EQ.1) IPATHL = -1
   IEMIT = 1
   FACT = 1.
   TIMEM = 0.0
   IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.
   DO 10 MOL = 1, NMOL
      WK(MOL) = WK(MOL)*FACT
10 CONTINUE
   WBROAD = WBROAD*FACT
   LAYR1 = LAYER
   WRITE (IPR,915) LAYR1,LAYER,KFILE,MFILE
   CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)
   DVXM = DV
   XKT = TAVE/RADCN2
   XKTBND = TBND/RADCN2
   IF (IPATHL.EQ.-1) STOP ' IPATH=-1 '
   IF (IPATHL.EQ.0) STOP ' IPATH=0 '
   IF (IPATHL.EQ.1) THEN
      XKTA = TZL/RADCN2
      XKTB = 0.
   ENDIF
   IF (IPATHL.EQ.2) STOP ' IPATH=2 '
   IF (IPATHL.EQ.3) THEN
      XKTA = TZU/RADCN2
      XKTB = 0.
   ENDIF
   WRITE (IPR,920) IPATHL,IANT
!
20 CONTINUE
!
   CALL CPUTIM (TIMEM1)
   CALL FLXIN (V1PO,V2PO,DVPO,NLIMO,KFILE,EMLAYR,TRALYR,KEOF,NPANLS)
   CALL CPUTIM (TIMEM2)
   TIMEM = TIMEM+TIMEM2-TIMEM1
   IF (KEOF.LE.0) GO TO 50
   VI = V1PO-DVPO
   VIDVBD = VI
   VIDVEM = VI
   VIDVRF = VI
   BBLAST = -1.
   EMLAST = -1.
   IF ((IPATHL.EQ.3).AND.(TBND.GT.0.)) THEN
!
      NLIM1 = 0
      NLIM2 = 0
      EMDUM = 0.
      BBDUM = 0.
      EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMDUM)

      if (dbg(24)) then
         print *, 'FLINIT:: ((IPATHL.EQ.3).AND.(TBND.GT.0.)):: NOT CHECKED'
         dbg(24) = .false.
      endif
      if (refl_flg .eq. 1) then
         call bufin (ifiledown, keof, emdown, nlimo)
      endif
!
      VI = V1PO
30    NLIM1 = NLIM2+1
!
      EMISIV = EMISFN(VI,DVPO,VIDV,EMDEL,EMLAST)
!
      IF (VIDV.GE.9.E+4) THEN
         NLIM2 = NLIMO+1
      ELSE
         NLIM2 = (VIDV-V1PO)/DVPO+1.001
      ENDIF
      NLIM2 = MIN(NLIM2,NLIMO)
!
!    Modified (2017) to properly handle downwelling reflected off surface
      DO 40 J = NLIM1, NLIM2
         BB = PLANCK(VI,XKTBND)
         NEWEM(J) = EMLAYR(J) + TRALYR(J) *                            &
         &             (BB*EMISIV + emdown(j)*(1.-emisiv))
!
!           Increment interpolation values
!
         EMISIV = EMISIV+EMDEL
         VI = VI + DVPO
40    CONTINUE
!
      IF (NLIM2.LT.NLIMO) GO TO 30
!
   ENDIF
   CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,NEWEM,NEWTR,MFILE,NPTS,NPANLS)
   GO TO 20
50 CALL CPUTIM (TIME1)
   TIME = TIME1-TIME
   WRITE (IPR,925) TIME,TIMEM
!
   if (refl_flg .eq. 1) then
      close(ifiledown)
   endif

   RETURN
!
900 FORMAT (' TIME AT THE START OF --FLINIT-- ',F10.3)
905 FORMAT (F10.8)
910 FORMAT ('0 ***** DIR. COSINE ***** ',/,7X,F10.8)
915 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,                  &
   &        '0 INPUT FILE =',I5,' OUTPUT FILE =',I5)
920 FORMAT ('0 IPATHL AND IANT',2I5)
925 FORMAT (' TIME REQUIRED FOR --FLINIT-- ',F10.3,' --FLXIN-- ',       &
   &        F10.3)
!
END SUBROUTINE FLINIT
!
!     ----------------------------------------------------------------
!
SUBROUTINE FLUXUP (NPTS,LFILE,MFILE,JPATHL,TBND)
!
   USE phys_consts, ONLY: radcn2
   IMPLICIT REAL*8           (V)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    14 AUGUST 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     M.J. IACONO
!                                     R.D. WORSHAM
!                                     J.L. MONCET
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
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   COMMON RADN(2410),TRAN(2410),RADO(0:5000)
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   &
   &              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),        &
   &               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,    &
   &               EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,        &
   &    pad_3,angle_path, secant_diffuse, secant_path, diffuse_fac
!
   character*1 surf_refl
   character*3 pad_3
!
   COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
   COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,        &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /XME/ v1_pad,v2_pad,dv_pad,n_pad,TRAO(0:5000)
   COMMON /RMRG/ XKT,XKTA,XKTB,SECNT
!
   DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)
   DIMENSION A1(10),A2(10),A3(10),A4(10)
   DIMENSION RADLYR(2),TRALYR(2),RADOI(2),TRAOI(2)
   DIMENSION WKSAV(35)
!
   CHARACTER*40 CYID
!
   EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                   &
   &            (OPNLHD(1),V1PO)
   EQUIVALENCE (RADO(1),RADOI(1)) , (TRAO(1),TRAOI(1))
   EQUIVALENCE (RADN(1),RADLYR(1)) , (TRAN(1),TRALYR(1)),              &
   &            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),                 &
   &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 &
   &            (FSCDID(16),LAYR1)
!
!
!
!
!     ************************************************************
!     ****** THIS SUBROUTINE DOES LAYER MERGE FOR RADIANCE  ******
!     ************************************************************
!
!     IPATHL = 3 IS FOR THE UPWELLING FLUX CASE
!
   CALL CPUTIM (TIME)
   WRITE (IPR,900) TIME
   NPANLS = 0
   TIMEM = 0.0
   TIMRD = 0.0
   TIMOT = 0.0
!
   CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)
   SECNT = SECANT
   LAY1SV = LAYR1
   DVL = DV
   PL = PAVE
   TL = TAVE
   WTOTL = 0.
!
   DO 10 MOL = 1, NMOL
      WTOTL = WTOTL+WK(MOL)
      WKSAV(MOL) = WK(MOL)
10 CONTINUE
!
   WTOTL = WTOTL+WBROAD
   WN2SAV = WBROAD
!
!     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)
!
   IF (iaersl.ge.1 .and. iaersl.ne.5)           &
   &     WRITE (CYID,'(5A8)') (YID(I),I=3,7)
   CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)
   IF (iaersl.ge.1 .and. iaersl.ne.5)                   &
   &                 READ (CYID,'(5A8)') (YID(I),I=3,7)
!
   IF (JPATHL.GE.1) IPATHL = JPATHL
   PLAY = PAVE
   TLAY = TAVE
   TAVK = TAVE
   DVK = DV
   FACT = 1.
!
   IF (DVL.EQ.DVK) THEN
      ITYPE = 0
   ELSEIF (DVL.GT.DVK) THEN
      ITYPE = DVK/(DVL-DVK)+0.5
   ELSE
!
!     DVL.LT.DVK
!
      ITYPE = -INT(DVL/(DVK-DVL)+0.5)
   ENDIF
   IF (ITYPE.LT.0) STOP ' FLUXUP; ITYPE LT 0 '
!
   WTOTK = 0.
   LAYR1 = LAY1SV
   WRITE (IPR,905) LAYR1,LAYER,KFILE,LFILE,MFILE
   IEMIT = 1
   DO 20 MOL = 1, NMOL
      WTOTK = WTOTK+WK(MOL)*FACT
      WK(MOL) = WK(MOL)*FACT+WKSAV(MOL)
20 CONTINUE
   WTOTK = WTOTK+WBROAD*FACT
   WBROAD = WBROAD*FACT+WN2SAV
   PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)
   TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)
   SECANT = SECNT
!
!     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES
!
   CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)
   DVXM = DV
   XKT = TAVK/RADCN2
!
   WRITE (IPR,910) IPATHL,IANT
!
   IF (IPATHL.EQ.3) THEN
      XKTA = TZU/RADCN2
      XKTB = 0.
   ELSE
      STOP ' FLUXUP; IPATHL '
   ENDIF
!
   ATYPE = ITYPE
   LL = ITYPE+1
   AP = 1.0/(ATYPE+1.0)
!
!     A1, A2, A3 AND A4 ARE THE CONSTANTS
!     FOR THE LAGRANGE 4 POINT INTERPOLATION
!
   DO 30 JPG = 1, ITYPE
      APG = JPG
      IPL = JPG+1
      P = 1.0-(AP*APG)
      A1(IPL) = -P*(P-1.0)*(P-2.0)/6.0
      A2(IPL) = (P**2-1.0)*(P-2.0)*0.5
      A3(IPL) = -P*(P+1.0)*(P-2.0)*0.5
      A4(IPL) = P*(P**2-1.0)/6.0
30 CONTINUE
!
!     *** BEGINNING OF LOOP THAT DOES MERGE  ***
!
   NPE = 0
   RADO(0) = 0.0
   TRAO(0) = 0.0
   V1PO = 0.0
   V2PO = 0.0
   DVPO = 0.0
!
40 CONTINUE
!
   CALL CPUTIM (TIMEM1)
   CALL FLXIN (V1P,V2P,DVP,NLIM,KFILE,RADLYR,TRALYR,KEOF,NPANLS)
   CALL CPUTIM (TIMEM2)
   TIMEM = TIMEM+TIMEM2-TIMEM1
   IF (KEOF.LE.0) GO TO 80
   II = 1
!
   IF (V2PO.LE.V2P+DVPO) THEN
50    CALL CPUTIM (TIMEM1)
      CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)
      IF (LEOF.LE.0) GO TO 60
      CALL BUFIN (LFILE,LEOF,RADOI(NPE+1),NLIMO)
      CALL BUFIN (LFILE,LEOF,TRAOI(NPE+1),NLIMO)
      CALL CPUTIM (TIMEM2)
      TIMRD = TIMRD+TIMEM2-TIMEM1
      NPE = NLIMO+NPE
      IF (V2PO.LE.V2P+DVPO) GO TO 50
   ENDIF
!
!     ZERO POINT OF FIRST PANEL
!
60 IF (RADO(0).EQ.0.0.AND.TRAO(0).EQ.0.0) THEN
      RADO(0) = RADO(1)
      TRAO(0) = TRAO(1)
   ENDIF
!
!     END POINT OF LAST PANEL
!
   IF (V2PO+DVPO.GE.V2) THEN
      RADO(NPE+1) = RADO(NPE)
      RADO(NPE+2) = RADO(NPE)
      TRAO(NPE+1) = TRAO(NPE)
      TRAO(NPE+2) = TRAO(NPE)
   ENDIF
!
   NPL = 1
!
!     NPL IS LOCATION OF FIRST ELEMENT ON ARRAYS RADO AND TRAO
!
   CALL FLUXNN (RADN,TRAN,RADO,TRAO,NLIM,V1P,DVP,IPATHL,               &
   &             A1,A2,A3,A4,LL,NPL)
!
   CALL CPUTIM (TIMEM1)
!
   IF (TBND.GT.0.) CALL EMBND (V1P,V2P,DVP,NLIM,RADN,TRAN,TBND)
!
   CALL EMOUT (V1P,V2P,DVP,NLIM,RADN,TRAN,MFILE,NPTS,NPANLS)
   CALL CPUTIM (TIMEM2)
   TIMOT = TIMOT+TIMEM2-TIMEM1
!
!     NPL IS NOW LOCATION OF FIRST ELEMENT TO BE USED FOR NEXT PASS
!
   IPL = -1
   DO 70 NL = NPL, NPE
      IPL = IPL+1
      RADO(IPL) = RADO(NL)
      TRAO(IPL) = TRAO(NL)
70 CONTINUE
!
   NPE = IPL
!
   GO TO 40
80 CONTINUE
!
   CALL CPUTIM (TIME1)
   TIM = TIME1-TIME
   WRITE (IPR,915) TIME1,TIM,TIMEM,TIMRD,TIMOT
!
   RETURN
!
900 FORMAT ('0 THE TIME AT THE START OF FLUXUP IS ',F12.3)
905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,     &
   &        ' MERGED WITH FILE ',I5,' ONTO FILE',I5)
910 FORMAT ('0 IPATHL AND IANT',2I5)
915 FORMAT ('0 THE TIME AT THE END OF FLUXUP IS ',F12.3,/,F12.3,        &
   &        ' SECS WERE REQUIRED FOR THIS MERGE  - FLXIN - ',F12.3,     &
   &        ' - READ - ',F12.3,' - EMOUT - ',F12.3)
!
END SUBROUTINE FLUXUP
!
!     ----------------------------------------------------------------
!
SUBROUTINE FLUXNN (RADLYR,TRALYR,RADO,TRAO,NLIM,V1P,DVP,            &
&                  IPATHL,A1,A2,A3,A4,LL,NPL)
!
   USE lblparams, ONLY: NDIM, ND2
   IMPLICIT REAL*8           (V)
!
!     THIS SUBROUTINE CALCULATES THE NEW RADIANCE AND TRANSMISSION
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    14 AUGUST 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!                       ALGORITHM:    R.D. WORSHAM
!                                     S.A. CLOUGH
!                                     M.J. IACONO
!                                     J.L. MONCET
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

! if this changes, make sure it is changed in subroutines in xmerge.f
!      parameter (ndim=2410, nd2=5000)

   DIMENSION RADLYR(NDIM),TRALYR(NDIM),RADO(0:ND2),TRAO(0:ND2),        &
   &          A1(*),A2(*),A3(*),A4(*)
!
   DATA I_1/1/
!
   LLM1 = LL-1
   LLM1 = MAX(LLM1,I_1)
!
!     LOOPING OVER POINTS WITH SAME WEIGHTS
!
   DO 30 NL = 1, LL
      IPL = (NPL+NL-1)-LLM1
      IF (NL.GT.1) IPL = IPL-1
!
      IF (NL.EQ.1) THEN
!
!     EXACT FREQUENCY - NO INTERPOLATION
!
         DO 10 I = NL, NLIM, LL
            IPL = IPL+LLM1
            RADLYR(I) = TRALYR(I)*RADO(IPL)+RADLYR(I)
            TRALYR(I) = TRALYR(I)*TRAO(IPL)
10       CONTINUE
!
!     NOT EXACT FREQUENCY - INTERPOLATE RESULT
!
      ELSE
!
         A1N = A1(NL)
         A2N = A2(NL)
         A3N = A3(NL)
         A4N = A4(NL)
!
         DO 20 I = NL, NLIM, LL
            IPL = IPL+LLM1
!
!     INTERPOLATE THE OLD RADIANCE
!
            RADLYR(I) = TRALYR(I)*(A1N*RADO(IPL-1)+A2N*RADO(IPL)+      &
            &                     A3N*RADO(IPL+1)+A4N*RADO(IPL+2))+RADLYR(I)
!
!     INTERPOLATE THE OLD TRANSMISSION
!
            TRALYR(I) = TRALYR(I)*(A1N*TRAO(IPL-1)+A2N*TRAO(IPL)+      &
            &                     A3N*TRAO(IPL+1)+A4N*TRAO(IPL+2))
20       CONTINUE
!
!
      ENDIF
!
30 CONTINUE
!
   NPL = IPL
!
   RETURN
!
END SUBROUTINE FLUXNN
!
!     ----------------------------------------------------------------
!
SUBROUTINE FLUXDN (NPTS,LFILE,MFILE,JPATHL,TBND,refl_flg)
!
   USE phys_consts, ONLY: radcn2
   IMPLICIT REAL*8           (V)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    14 AUGUST 1991
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     M.J. IACONO
!                                     R.D. WORSHAM
!                                     J.L. MONCET
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
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   COMMON RADN(2410),TRAN(2410),RADLYR(-1:4818)
   COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   &
   &              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),        &
   &               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,    &
   &               EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
   COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,        &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /XMI/ v1_dum,v2_dum,dv_dum,n_dum,TRALYR(-1:4818)
   COMMON /RMRG/ XKT,XKTA,XKTB,SECNT
!
   DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)
   DIMENSION A1(0:100),A2(0:100),A3(0:100),A4(0:100)
   DIMENSION RADO(2),TRAO(2)
   DIMENSION WKSAV(35)
   dimension raddown(2410)
!
   CHARACTER*40 CYID
   integer refl_flg
   integer ifiledown
!
   EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                   &
   &            (OPNLHD(1),V1PO)
   EQUIVALENCE (RADN(1),RADO(1)) , (TRAN(1),TRAO(1)),                  &
   &            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),                 &
   &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 &
   &            (FSCDID(16),LAYR1)
!
!     ************************************************************
!     ****** THIS SUBROUTINE DOES LAYER MERGE FOR EMISSION  ******
!     ****** USING FOUR POINT GENERAL INTERPOLATION         ******
!     ************************************************************
!
   CALL CPUTIM (TIME)
   WRITE (IPR,900) TIME
   NPANLS = 0
   TIMEM = 0.0
   TIMRD = 0.0
   TIMTB = 0.0
   TIMOT = 0.0
   jdownsum = 0
   iend = 0
   nlimlast = 0
!
   CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)
   SECNT = SECANT
   dircos = 1./secnt
   if (dircos .gt. 0.6) then
      ifiledown = 601
   elseif (dircos .gt. 0.3) then
      ifiledown = 602
   else
      ifiledown = 603
   endif


   DVL = DV
   LAY1SV = LAYR1
   PL = PAVE
   TL = TAVE
   WTOTL = 0.
   DO 10 MOL = 1, NMOL
      WTOTL = WTOTL+WK(MOL)
      WKSAV(MOL) = WK(MOL)
10 CONTINUE
   WTOTL = WTOTL+WBROAD
   WN2SAV = WBROAD
!
!     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)
!
   IF (iaersl.ge.1 .and. iaersl.ne.5)                           &
   &                 WRITE (CYID,'(5A8)') (YID(I),I=3,7)
   CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)
   IF (iaersl.ge.1 .and. iaersl.ne.5)                           &
   &     READ (CYID,'(5A8)') (YID(I),I=3,7)
   IF (JPATHL.GE.1) IPATHL = JPATHL
   PLAY = PAVE
   TLAY = TAVE
!
   IF (IPATHL.NE.1) STOP ' FLUXDN: IPATHL .NE. 1 '
!
   XKT = TAVE/RADCN2
   XKTA = TZL/RADCN2
   XKTB = 0.
   DVK = DV
   LAYR1 = LAY1SV
   FACT = 1.
   IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.
   ATYPE = 9.999E09
   IF (DVK.EQ.DVL) ATYPE = 0.
   IF (DVL.GT.DVK) ATYPE = DVK/(DVL-DVK)+0.5
   IF (DVL.LT.DVK) ATYPE = -DVL/(DVK-DVL)-0.5
!
!     IF (ATYPE .GT. 0) STOP  ' FLUXDN; ATYPE GT 0 '
!
   WTOTK = 0.
   WRITE (IPR,905) LAYR1,LAYER,KFILE,LFILE,MFILE,ATYPE
   IEMIT = 1
   DO 20 MOL = 1, NMOL
      WTOTK = WTOTK+WK(MOL)*FACT
      WK(MOL) = WK(MOL)*FACT+WKSAV(MOL)
20 CONTINUE
   WTOTK = WTOTK+WBROAD*FACT
   WBROAD = WBROAD*FACT+WN2SAV
   PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)
   TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)
   SECANT = SECNT
   DV = DVL
!
!     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES
!
   CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)
   DVXM = DV
!
   IF (ATYPE.EQ.0.) THEN
!
!     1/1 RATIO ONLY
!
!
30    CONTINUE
      CALL CPUTIM (TIMEM1)
      CALL FLXIN (V1P,V2P,DVP,NLIM,KFILE,RADLYR(1),TRALYR(1),KEOF,     &
      &               NPANLS)
      CALL CPUTIM (TIMEM2)
      TIMEM = TIMEM+TIMEM2-TIMEM1
      IF (KEOF.LE.0) GO TO 110
      CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)
      CALL BUFIN (LFILE,LEOF,RADO(1),NLIMO)
      CALL BUFIN (LFILE,LEOF,TRAO(1),NLIMO)
      CALL CPUTIM (TIMEM3)

      TIMRD = TIMRD+TIMEM3-TIMEM2
      DO 40 I = 1, NLIM
         RADN(I) = RADO(I)*TRALYR(I)+RADLYR(I)
         TRAN(I) = TRALYR(I)*TRAO(I)
40    CONTINUE

! Write out downwelling at surface for use by upwelling calculation.
      if (refl_flg .eq. 1) then
         call bufout (ifiledown,radn,nlim)
      endif

      CALL CPUTIM (TIMEM1)
      IF (TBND.GT.0.)                                                  &
      &        CALL EMBND (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,TBND)
!
      CALL CPUTIM (TIMEM2)
      TIMTB = TIMTB+TIMEM2-TIMEM1
      CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,MFILE,NPTS,NPANLS)
      CALL CPUTIM (TIMEM3)
      TIMOT = TIMOT+TIMEM3-TIMEM2
      GO TO 30
!
   ENDIF
!
!     ALL RATIOS EXCEPT 1/1
!
   DO 50 JP = 0, 100
      APG = JP
      P = 0.01*APG
!
!     THE FOLLOW ARE THE CONSTANTS FOR THE LAGRANGE 4 POINT
!     INTERPOLATION
!
      A1(JP) = -P*(P-1.0)*(P-2.0)/6.0
      A2(JP) = (P**2-1.0)*(P-2.0)*0.5
      A3(JP) = -P*(P+1.0)*(P-2.0)*0.5
      A4(JP) = P*(P**2-1.0)/6.0
50 CONTINUE
!
!     *** BEGINNING OF LOOP THAT DOES MERGE  ***
!
   NPE = 0
   RADLYR(0) = 0.0
   TRALYR(0) = 0.0
   V1P = 0.0
   V2P = 0.0
   DVP = 0.0
   V1PO = 0.0
   V2PO = 0.0
   DVPO = 0.0
   KEOF = 1
!
60 CONTINUE
!
   CALL CPUTIM (TIMEM1)
   CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)
   IF (LEOF.LE.0) GO TO 110
   CALL BUFIN (LFILE,LEOF,RADO(1),NLIMO)
   CALL BUFIN (LFILE,LEOF,TRAO(1),NLIMO)
   CALL CPUTIM (TIMEM2)
   TIMRD = TIMRD+TIMEM2-TIMEM1
   II = 1
!
!      IF (V2P.LE.V2PO+DVP .AND. KEOF.GT.0) THEN
   IF (V2P.LE.V2PO+DVP) then
!         if (refl_flg .eq. 1.and.dvp.ne.0.0) then
!            call bufout (ifiledown,raddown,nlim)
!         endif
      if (KEOF.GT.0) THEN
70       CALL CPUTIM (TIMEM2)
         CALL FLXIN (V1P,V2P,DVP,NLIM,KFILE,RADLYR(NPE+1),              &
         &               TRALYR(NPE+1),KEOF,NPANLS)
         CALL CPUTIM (TIMEM3)
         TIMEM = TIMEM+TIMEM3-TIMEM2
         IF (KEOF.LE.0) GO TO 80
         V1P = V1P- REAL(NPE)*DVP

         NPE = NLIM+NPE
      endif
      IF (V2P.LE.V2PO+DVP) GO TO 70
   ENDIF
!
!     ZERO POINT OF FIRST PANEL
!
80 IF (RADLYR(0).EQ.0.0.AND.TRALYR(0).EQ.0.0) THEN
      TRALYR(-1) = TRALYR(1)
      RADLYR(-1) = RADLYR(1)
      TRALYR(0) = TRALYR(1)
      RADLYR(0) = RADLYR(1)
   ENDIF
!
!     END POINT OF LAST PANEL
!
   IF (V2P+DVP.GE.V2) THEN
      TRALYR(NPE+1) = TRALYR(NPE)
      RADLYR(NPE+1) = RADLYR(NPE)
      TRALYR(NPE+2) = TRALYR(NPE)
      RADLYR(NPE+2) = RADLYR(NPE)
   ENDIF
!
!     NPL IS LOCATION OF FIRST ELEMENT ON ARRAYS RADO AND TRAO
!
   NPL = 1
!
   RATDV = DVL/DVK
!
!     FJJ IS OFFSET BY 2. FOR ROUNDING PURPOSES
!
   FJ1DIF = (V1PO-V1P)/DVP+1.+2.
!
!     ***** BEGINNING OF LOOP THAT DOES MERGE  *****
!
   DO 90 II = 1, NLIMO
!
      FJJ = FJ1DIF+RATDV* REAL(II-1)
      JJ =  INT(FJJ)-2
!
      JP = (FJJ- REAL(JJ))*100.-199.5
!
!     INTERPOLATE THE OLD TRANSMISSION
!
      TRNLYR = A1(JP)*TRALYR(JJ-1)+A2(JP)*TRALYR(JJ)+                  &
      &            A3(JP)*TRALYR(JJ+1)+A4(JP)*TRALYR(JJ+2)
!
!     INTERPOLATE THE OLD EMISSION
!
      RADN(II) = RADO(II)*TRNLYR+(A1(JP)*RADLYR(JJ-1)+                 &
      &              A2(JP)*RADLYR(JJ)+A3(JP)*RADLYR(JJ+1)+                &
      &              A4(JP)*RADLYR(JJ+2))
!
      TRAN(II) = TRNLYR*TRAO(II)

! Write out downwelling at surface for use by upwelling calculation.
      if (refl_flg.eq.1 .and.                                          &
      &        jj .ne. int(FJ1DIF+RATDV*REAL(II-2))-2) then
         if (iend.eq.0 .and. nlim.ne.2400) then
            iend = 1
            nlimlast = nlim
         endif

         jdownsum = jdownsum + 1
         raddown(jdownsum)=radn(ii)

         if (jdownsum .eq. 2400) then
! Case where we are at the end of a coarse panel that is not the final panel
            call bufout (ifiledown,raddown,jdownsum)
            jdownsum = 0
         elseif (nlim.ne.2400 .and. jdownsum.eq.nlimlast) then
! Case where the last coarse panel point has been reached
            call bufout (ifiledown,raddown,jdownsum)
            jdownsum = 0
         endif
      elseif (jdownsum .eq. nlimlast-1 .and. ii .eq. nlimo) then
! Use the final fine panel point to match the last coarse point
         jdownsum = jdownsum + 1
         raddown(jdownsum)=radn(ii)
         call bufout (ifiledown,raddown,jdownsum)
      endif
!
90 CONTINUE
!
   NPL = JJ-1
!
   CALL CPUTIM (TIMEM1)
   IF (TBND.GT.0.) CALL EMBND (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,TBND)
!
   CALL CPUTIM (TIMEM2)
   TIMTB = TIMTB+TIMEM2-TIMEM1
   CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,MFILE,NPTS,NPANLS)
   CALL CPUTIM (TIMEM3)
   TIMOT = TIMOT+TIMEM3-TIMEM2
!
!     NPL IS NOW LOCATION OF FIRST ELEMENT TO BE USED FOR NEXT PASS
!
   IPL = -2
   DO 100 NL = NPL, NPE
      IPL = IPL+1
      TRALYR(IPL) = TRALYR(NL)
      RADLYR(IPL) = RADLYR(NL)
100 CONTINUE
!
   V1P = V1P+ REAL(NPL+1)*DVP
   NPE = IPL
!
   GO TO 60
110 CONTINUE
!
   close(ifiledown)
   CALL CPUTIM (TIME1)
   TIM = TIME1-TIME
   WRITE (IPR,910) TIME1,TIM,TIMEM,TIMRD,TIMTB,TIMOT
!
   RETURN
!
900 FORMAT ('0 THE TIME AT THE START OF FLUXDN IS ',F12.3)
905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,     &
   &        ' MERGED WITH FILE ',I5,' ONTO FILE',I5,'  WITH XTYPE=',    &
   &        G15.5)
910 FORMAT ('0 THE TIME AT THE END OF FLUXDN IS ',F12.3/F12.3,          &
   &        ' SECS WERE REQUIRED FOR THIS MERGE  - FLXIN - ',F12.3,     &
   &        ' - READ - ',F12.3,' - EMBND - ',F12.3,' - EMOUT - ',       &
   &        F12.3)
!
END SUBROUTINE FLUXDN
!
!     ----------------------------------------------------------------
!
SUBROUTINE GETEXT (IEXFIL,LYRNOW,IEMITT)
!
   IMPLICIT REAL*8           (V)
!
!     THIS SUBROUTINE HAS BEEN MODIFIED TO ALSO READ IN THE ASYMMETRY
!     FACTORS AND TO SEARCH FOR THE PROPER LAYER IF DESIRED.
!
!                                          JAN 1986 (A.E.R. INC.)
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFLD,        &
   &              NLTEFL,LNFIL4,LNGTH4
!
!     ROUTINE TO BUFFER IN AEROSOL ABSORPTION AND SCATTERRING
!     FROM TAPE 20 INTO COMMON ABSORB SCATTR, AND ASYMT
!
   COMMON /PNLHDR/ V1P,V2P,DVP,NLIM,LDUM
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /FILHDA/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   &
   &                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
   COMMON /ABSORA/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(2025)
   COMMON /SCATTA/ V1SC,V2SC,DVSC,NPTSC,SCTTR(2025)
   COMMON /ASYMMA/ V1AS,V2AS,DVAS,NPTAS,ASYMT(2025)
   DIMENSION APNLHD(2),AFILHD(2)
!
   CHARACTER*40 CEXT
!
   EQUIVALENCE (APNLHD(1),V1P) , (AFILHD(1),XID(1))
!
   IF (IEMITT.EQ.0) THEN
      CALL BUFIN (IEXFIL,IEOF,AFILHD(1),NFHDRF)
!
!     MOVE YID INTO EXTID
!
      WRITE (CEXT,'(5A8)') (YID(I),I=3,7)
      READ (CEXT,'(10A4)') EXTID
   ENDIF
!
   IF (LYRNOW.EQ.1.AND.IEMITT.EQ.1) THEN
      REWIND IEXFIL
      CALL BUFIN (IEXFIL,IEOF,AFILHD(1),NFHDRF)
!
!     MOVE YID INTO EXTID
!
      WRITE (CEXT,'(5A8)') (YID(I),I=3,7)
      READ (CEXT,'(10A4)') EXTID
   ENDIF
!
   LAYER = 0
!
10 DO 20 I = 1, 2025
      ABSRB(I) = 0.
      SCTTR(I) = 0.
      ASYMT(I) = 0.
20 CONTINUE
!
   CALL BUFIN (IEXFIL,IEOF,APNLHD(1),NPHDRF)
!
   LAYER = LAYER+1
   IF (IEOF.LE.0) STOP 'GETEXT; IEXFIL EMPTY'
!
   CALL BUFIN (IEXFIL,IEOF,ABSRB(1),NLIM)
   CALL BUFIN (IEXFIL,IEOF,SCTTR(1),NLIM)
   CALL BUFIN (IEXFIL,IEOF,ASYMT(1),NLIM)
!
   V1ABS = V1P
   V1SC = V1P
   V1AS = V1P
   V2ABS = V2P
   V2SC = V2P
   V2AS = V2P
   DVABS = DVP
   DVSC = DVP
   DVAS = DVP
   NPTABS = NLIM
   NPTSC = NLIM
   NPTAS = NLIM
!
   RETURN
!
END SUBROUTINE GETEXT
!
!     ----------------------------------------------------------------
!
SUBROUTINE ADARSL (NNPTS,IEXFIL,MFILE,IAFIL,IEMIT)
!
   USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2, MXMOL, MXTRAC
   IMPLICIT REAL*8           (V)
!
!     ROUTINE TO ADD ABSORPTION AND SCATTERING TO THE TRANSMITTANCE
!     VALUES AT EACH POINT. THE AEROSOL VALUES ARE STORED IN
!     COMMON ABSORB AND COMMON SCATTR.
!
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                  &
!     &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
!
   COMMON R1(2410)
   COMMON /ABSPNL/ V1P,V2P,DVP,NLIM,NSHFT,NPTS
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      &
   &              EXTID(10)
   CHARACTER*8  EXTID

   COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),            &
   &                MSPNL1(MXLAY),                                      &
   &                MSLAY1,ISFILE,JSFILE,KSFILE,LSFILE,MSFILE,IEFILE,   &
   &                JEFILE,KEFILE
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFLD,IEXFLD,        &
   &              NLTEFL,LNFIL4,LNGTH4
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   &
   &                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
!
   EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHD(1),V1P)
   EQUIVALENCE (FSCDID(4),IAERSL)
!
   CHARACTER*40 CEXT,CYID
!
   DIMENSION PNLHD(4),XFILHD(2)
!
   XKT0 = 0.6951*296.
   CALL GETEXT (IEXFIL,LAYRS,IEMIT)
   CALL BUFIN (MFILE,IEOF,XFILHD(1),NFHDRF)
!
!     FOR AEROSOL RUNS, MOVE EXTID INTO YID
!
   WRITE (CEXT,'(10A4)') EXTID
   WRITE (CYID,'(5A8)') (YID(I),I=3,7)
   CYID(19:40) = CEXT(19:40)
   READ (CYID,'(5A8)') (YID(I),I=3,7)
!
   CALL BUFOUT (IAFIL,XFILHD(1),NFHDRF)
   NPANLS = 0
   IF (NOPR.EQ.0) WRITE (IPR,900) XID,(YID(N),N=1,2)
10 CALL BUFIN (MFILE,IEOF,PNLHD(1),NPHDRF)
   IF (IEOF.LE.0) GO TO 40
   CALL BUFIN (MFILE,IEOF,R1(1),NLIM)
   IF (NPANLS.EQ.0) VIDV = V1P-DVP
   VAER = VIDV
   VITST = VIDV
   NPTS = NNPTS
   RDLAST = -1.
   NLIM1 = 0
   NLIM2 = 0
   RDDUM = 0.
   AF = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
   RDF = RADFNI(VI,DVP,XKT0,VITST,RDEL,RDDUM)
   IAFRD = 0
   IF (VITST.GT.VAER) IAFRD = 1
!
20 NLIM1 = NLIM2+1
!
   VI = V1P+ REAL(NLIM1-1)*DVP
   IF (IAFRD.EQ.0) THEN
      RADFN0 = RADFNI(VI,DVP,XKT0,VIDV,RDEL,RDLAST)
      VAER = -VIDV
      EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)
   ELSE
      EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)
      VITST = -VIDV
      RADFN0 = RADFNI(VI,DVP,XKT0,VITST,RDEL,RDLAST)
   ENDIF
!
   NLIM2 = (VIDV-V1P)/DVP+1.001
   NLIM2 = MIN(NLIM2,NLIM)
!
   DO 30 K = NLIM1, NLIM2
      R1(K) = R1(K)+EXT*RADFN0
!
!        Increment interpolation values
!
      EXT = EXT+ADEL
      RADFN0 = RADFN0+RDEL
30 CONTINUE
!
   IF (NLIM2.LT.NLIM) GO TO 20
!
   CALL ABSOUT (V1P,V2P,DVP,NLIM,1,IAFIL,NPTS,R1,NPANLS)
   GO TO 10
40 CONTINUE
   IAERSL = 2
   RETURN
!
900 FORMAT ('1',5X,' AEROSOLS'//1X,10A8,2X,2(1X,A8,1X)//,5X,            &
   &        'FILE 20 AEROSOL EXTINCTIONS ADDED TO FILE 12 SENT TO ',    &
   &        'FILE 14')
!
END SUBROUTINE ADARSL
!
!     ----------------------------------------------------------------
!
FUNCTION AERF (VI,DVI,VINXT,ADEL,TAUSCT,TDEL,GFACT,GDEL)
!
   IMPLICIT REAL*8           (V)
!
!     THIS FUNCTION CORRELATES THE AEROSOL FREQ. WITH THE LBLRTM
!     FREQ.  AND ADDS THE ABSORPTION TO THE
!     SCATTERING TO PRODUCE THE EXTINCTION
!
!     THIS FUNCTION HAS BEEN MODIFIED TO RETURN THE SCATTERING
!     SEPARATELY, AND TO ALSO RETURN THE ASYMMETRY FACTOR.
!
!                                         JAN 1986 (A.E.R. INC.)
!
   COMMON /ABSORA/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(2025)
   COMMON /SCATTA/ V1SC,V2SC,DVSC,NPTSC,SCTTR(2025)
   COMMON /ASYMMA/ V1AS,V2AS,DVAS,NPTAS,ASYMT(2025)
!
   character*8      XID,       HMOLID,      YID
   real*8               SECANT,       XALTZ
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       &
   &                WK(60),PZL,PZU,TZL,TZU,WN2   ,DV ,V1 ,V2 ,TBOUND,   &
   &                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
!
   EQUIVALENCE (XFILHD(1),XID(1))
   DIMENSION XFILHD(2)
!
   VDIF = VI-V1ABS
   IAER = VDIF/DVABS+1.00001
   VAER = V1ABS+DVABS* REAL(IAER-1)
   DERIVS = (SCTTR(IAER+1)-SCTTR(IAER))/DVABS
   DERIVA = (ASYMT(IAER+1)-ASYMT(IAER))/DVABS
   DERIV = DERIVS+(ABSRB(IAER+1)-ABSRB(IAER))/DVABS
!
!     TAUSCT IS THE SCATTERING TERM
!
   TAUSCT = SCTTR(IAER)+DERIVS*(VI-VAER)
!
!     GFACT IS THE ASYMMETRY FACTOR
!
   GFACT = ASYMT(IAER)+DERIVA*(VI-VAER)
   AERF = SCTTR(IAER)+ABSRB(IAER)+DERIV*(VI-VAER)
!
!     ADEL, TDEL & GDEL ARE THE CHANGE PER DVI
!
   ADEL = DERIV*DVI
   TDEL = DERIVS*DVI
   GDEL = DERIVA*DVI
!
!     Set next point for interpolation to be the point
!     corresponding to the next element of the SCTTR,
!     ASYMT, & ABSRB arrays
!
   VINXT = VAER+DVABS
!
   RETURN
!
END FUNCTION AERF
!
!     ----------------------------------------------------------------
!
SUBROUTINE LINTCO(V1,Z1,V2,Z2,VINT,ZINT,ZDEL)
!
!     Linearly interpolates emission and reflection values which
!     are directly read in from ASCII files
!
   IMPLICIT REAL*8           (V)
!
!     ZDEL is the slope of the line
!
   ZDEL = (Z2-Z1)/(V2-V1)
!
!     ZCEPT is the intercept for V = 0.0
!
   ZCEPT = Z1 - ZDEL*V1
!
!     Calculate ZINT value at VINT
!
   ZINT = ZDEL*VINT + ZCEPT
!
   RETURN
!
END SUBROUTINE LINTCO
!     ----------------------------------------------------------------
