C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
C
C     ----------------------------------------------------------------
C
      SUBROUTINE XMERGE (NPTS,LFILE,MFILE,JPATHL)                         H00010
C                                                                         H00020
      IMPLICIT DOUBLE PRECISION (V)                                     ! H00030
C                                                                         H00040
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H00050
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        H00060
     *              NLTEFL,LNFIL4,LNGTH4                                  H00070
C                                                                         H00080
C     XMERGE CALL ABSMRG,EMINIT,RADMRG                                    H00090
C                                                                         H00100
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H00110
     *              AVFIX,LAYER,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,        H00120
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H00130
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H00140
     *              EXTID(10)                                             H00150
C                                                                         H00160
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H00170
C                                                                         H00180
      COMMON /HVERSN/  HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,
     *                HVROPR,HVRPST,HVRPLT,HVRTST,HVRUTL,HVRXMR
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       H00190
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   H00200
     *                EMISIV,FSCDID(17),NMOL,LAYHDR,YI1,YID(10),LSTWDF    H00210
C                                                                         H00220
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP                   H00230
      COMMON /XME/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(4996)                   H00240
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H00250
C                                                                         H00260
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),                H00270
     *            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),                H00280
     *            (FSCDID(5),IEMIT) , (FSCDID(7),IPLOT),                  H00290
     *            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                  H00300
     *            (FSCDID(11),IMRG) , (FSCDID(16),LAYR1),                 H00310
     *            (FSCDID(17),NLAYHD)                                     H00320
C                                                                         H00330
      CHARACTER*8 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *            HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR
C
C     ASSIGN SCCS VERSION NUMBER TO MODULE 
C
      HVRXMR = '$Revision$'
C
      IOD = 0                                                             H00340
C                                                                         H00350
      IEXFIL = 20                                                         H00360
      IAFIL = 14                                                          H00370
C                                                                         H00380
C     WHEN IAERSL EQUALS 2 CALL ADARSL TO ADD ABSORPTION AND SCATTERING   H00390
C     TO COMMON BLOCKS FOR USE IN A TRANSMITTANCE RUN                     H00400
C                                                                         H00410
      IF (IEMIT.NE.1) THEN                                                H00420
         IF (LAYER.EQ.1) THEN                                             H00430
            CALL ABINIT (NPTS,MFILE,JPATHL)                               H00440
         ELSE                                                             H00450
            WRITE (IPR,900) XID,(YID(M),M=1,2)                            H00460
            CALL ABSMRG (NPTS,LFILE,MFILE,JPATHL)                         H00470
         ENDIF                                                            H00480
      ELSE                                                                H00490
C                                                                         H00500
C     IEMIT = 1 TO REACH THIS STATEMENT                                   H00510
C                                                                         H00520
         WRITE (IPR,900) XID,(YID(M),M=1,2)                               H00530
         IF (LAYER.EQ.1.AND.IAERSL.GE.1) REWIND IEXFIL                    H00540
         IF (IAERSL.GE.1) CALL GETEXT (IEXFIL,LAYER,IEMIT)                H00550
C                                                                         H00560
         TBND = 0.                                                        H00570
         IF (IMRG.NE.36) THEN                                             H00580
            IF (LAYER.EQ.1) THEN                                          H00590
               IF (IPATHL.EQ.1) TBND = TMPBND                             H00600
               CALL EMINIT (NPTS,MFILE,JPATHL,TBND)                       H00610
            ELSE                                                          H00620
               IF (IPATHL.EQ.3.AND.LAYER.EQ.LH2) TBND = TMPBND            H00630
               CALL RADMRG (NPTS,LFILE,MFILE,JPATHL,TBND)                 H00640
            ENDIF                                                         H00650
         ELSE                                                             H00660
            IF (LAYER.EQ.1) THEN                                          H00670
               TBND = TMPBND                                              H00680
               CALL FLINIT (NPTS,MFILE,JPATHL,TBND)                       H00690
            ELSE                                                          H00700
               CALL FLUXUP (NPTS,LFILE,MFILE,JPATHL,TBND)                 H00710
            ENDIF
         ENDIF                                                            H00720
C                                                                         H00730
      ENDIF                                                               H00740
C                                                                         H00750
      RETURN                                                              H00760
C                                                                         H00770
  900 FORMAT (///,1X,10A8,2X,2(1X,A8,1X))                                 H00780
C                                                                         H00790
      END                                                                 H00800
C
C     ----------------------------------------------------------------
C
      SUBROUTINE XMERGI (NPTS,LFILE,MFILE,JPATHL)                         H00810
C                                                                         H00820
      IMPLICIT DOUBLE PRECISION (V)                                     ! H00830
C                                                                         H00840
      PARAMETER (MXFSC=200,MXLAY=MXFSC+3,MXZMD=200,MXPDIM=MXLAY+MXZMD,
     *           IM2=MXPDIM-2,MXMOL=35,MXTRAC=22)
C
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H00850
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        H00860
     *              NLTEFL,LNFIL4,LNGTH4                                  H00870
C                                                                         H00880
C     XMERGI CALL ABINIT,ABSINT,ABSOUT                                    H00890
C                                                                         H00900
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H00910
     *              AVFIX,LAYER,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,        H00920
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H00930
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H00940
     *              EXTID(10)                                             H00950
      COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),
     *                MSPNL1(MXLAY),                                      H00960
     *                MSLAY1,ISFILE,JSFILE,KSFILE,LSFILE,MSFILE,IEFILE,   H00970
     *                JEFILE,KEFILE                                       H00980
C                                                                         H00990
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H01000
C                                                                         H01010
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       H01020
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   H01030
     *                EMISIV,FSCDID(17),NMOL,LAYHDR,YI1,YID(10),LSTWDF    H01040
C                                                                         H01050
      COMMON /XMI/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(4815)                   H01060
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H01070
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP                   H01080
C                                                                         H01090
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),                H01100
     *            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),                H01110
     *            (FSCDID(5),IEMIT) , (FSCDID(7),IPLOT),                  H01120
     *            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                  H01130
     *            (FSCDID(11),IMRG) , (FSCDID(16),LAYR1),                 H01140
     *            (FSCDID(17),NLAYHD)                                     H01150
C                                                                         H01160
      IF (IEMIT.EQ.1) GO TO 10                                            H01170
      IF (LAYER.EQ.LH1.AND.IANT.NE.-1) THEN                               H01180
         CALL ABINIT (NPTS,MFILE,JPATHL)                                  H01190
      ELSE                                                                H01200
C                                                                         H01210
         WRITE (IPR,900) XID,(YID(M),M=1,2)                               H01220
         CALL ABSINT (NPTS,LFILE,MFILE,JPATHL)                            H01230
      ENDIF                                                               H01240
C                                                                         H01250
      GO TO 20                                                            H01260
C                                                                         H01270
C     IEMIT = 1 TO REACH THIS STATEMENT                                   H01280
C                                                                         H01290
   10 CONTINUE                                                            H01300
      WRITE (IPR,900) XID,(YID(M),M=1,2)                                  H01310
      IF (LAYER.EQ.1.AND.IAERSL.GE.1) REWIND IEXFIL                       H01320
      IF (IAERSL.GE.1) CALL GETEXT (IEXFIL,LAYER,IEMIT)                   H01330
C                                                                         H01340
      TBND = 0.                                                           H01350
C                                                                         H01360
      IF (IMRG.NE.35) THEN                                                H01370
         IF (LAYER.EQ.LH1.AND.IANT.NE.-1) THEN                            H01380
            IF (JPATHL.EQ.1.AND.LAYER.EQ.1)   TBND = TMPBND               H01390
            CALL EMINIT (NPTS,MFILE,JPATHL,TBND)                          H01400
         ELSE                                                             H01410
            IF (JPATHL.EQ.3.AND.LAYER.EQ.LH2) TBND = TMPBND               H01420
            CALL RADINT (NPTS,LFILE,MFILE,JPATHL,TBND)                    H01430
         ENDIF                                                            H01440
      ELSE                                                                H01450
         IF (LAYER.EQ.LH1.AND.IANT.NE.-1) THEN                            H01460
            TBND = TMPBND
            CALL FLINIT (NPTS,MFILE,JPATHL,TBND)                          H01470
         ELSE                                                             H01480
            CALL FLUXDN (NPTS,LFILE,MFILE,JPATHL,TBND)                    H01490
         ENDIF                                                            H01500
      ENDIF                                                               H01510
C                                                                         H01520
   20 CONTINUE                                                            H01530
C                                                                         H01540
      RETURN                                                              H01550
C                                                                         H01560
  900 FORMAT (///,1X,10A8,2X,2(1X,A8,1X))                                 H01570
C                                                                         H01580
      END                                                                 H01590
C
C     ----------------------------------------------------------------
C
      SUBROUTINE ABINIT (NPTS,MFILE,JPATHL)                               H01600
C                                                                         H01610
      IMPLICIT DOUBLE PRECISION (V)                                     ! H01620
C                                                                         H01630
      COMMON ODLAY(-2:2407)                                               H01640
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H01650
     *              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       H01660
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H01670
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H01680
     *              EXTID(10)                                             H01690
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H01700
C                                                                         H01710
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H01720
C                                                                         H01730
      COMMON /ABSHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       H01740
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   H01750
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    H01760
      COMMON /ABSPNL/ V1P,V2P,DVP,NLIM,NSHFT,NPNTS                        H01770
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H01780
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        H01790
     *              NLTEFL,LNFIL4,LNGTH4                                  H01800
C                                                                         H01810
      DIMENSION XFILHD(2),PNLHDR(2)                                       H01820
      DIMENSION ODLAYR(2)                                                 H01830
C                                                                         H01840
      CHARACTER*40 CEXT,CYID                                              H01850
C                                                                         H01860
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P)                    H01870
      EQUIVALENCE (ODLAY(1),ODLAYR(1)) , (FSCDID(4),IAERSL),              H01880
     *            (FSCDID(5),IEMIT) , (FSCDID(7),IPLOT),                  H01890
     *            (FSCDID(8),IPATHL) , (FSCDID(16),LAYR1)                 H01900
C                                                                         H01910
C                                                                         H01920
C     ***********************************************************         H01930
C     ****** THIS SUBROUTINE INITALIZES MERGE FOR OPTICAL  ******         H01940
C     ****** DEPTHS                                        ******         H01950
C     ***********************************************************         H01960
C                                                                         H01970
      CALL CPUTIM (TIME)                                                  H01980
      WRITE (IPR,900) TIME                                                H01990
      NPANLS = 0                                                          H02000
C                                                                         H02010
      CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)                            H02020
      IF (JPATHL.GE.1) IPATHL = JPATHL                                    H02030
      PLAY = PAVE                                                         H02040
      TLAY = TAVE                                                         H02050
C                                                                         H02060
C     FOR AEROSOL RUNS, MOVE EXTID INTO YID                               H02070
C                                                                         H02080
      IF (IAERSL.GT.0) THEN                                               H02090
         WRITE (CEXT,'(10A4)') EXTID                                      H02100
         WRITE (CYID,'(5A8)') (YID(I),I=3,7)                              H02110
         CYID(19:40) = CEXT(19:40)                                        H02120
         READ (CYID,'(5A8)') (YID(I),I=3,7)                               H02130
      ENDIF                                                               H02140
C                                                                         H02150
      IEMIT = 0                                                           H02160
      FACT = 1.                                                           H02170
      IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.                            H02180
      DO 10 MOL = 1, NMOL                                                 H02190
         WK(MOL) = WK(MOL)*FACT                                           H02200
   10 CONTINUE                                                            H02210
      WBROAD = WBROAD*FACT                                                H02220
      LAYR1 = LAYER                                                       H02230
      WRITE (IPR,905) LAYR1,KFILE,MFILE                                   H02240
C                                                                         H02250
      CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)                                H02260
      DVXM = DV                                                           H02270
C                                                                         H02280
   20 CONTINUE                                                            H02290
C                                                                         H02300
      CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)                            H02310
      IF (KEOF.LE.0) GO TO 40                                             H02320
      CALL BUFIN (KFILE,KEOF,ODLAYR(1),NLIM)                              H02330
C                                                                         H02340
      IF (FACT.EQ.2.) THEN                                                H02350
         DO 30 I = 1, NLIM                                                H02360
            ODLAYR(I) = ODLAYR(I)+ODLAYR(I)                               H02370
   30    CONTINUE                                                         H02380
      ENDIF                                                               H02390
C                                                                         H02400
      CALL ABSOUT (V1P,V2P,DVP,NLIM,1,MFILE,NPTS,ODLAYR,NPANLS)           H02410
      GO TO 20                                                            H02420
C                                                                         H02430
   40 CONTINUE                                                            H02440
C                                                                         H02450
      CALL CPUTIM (TIME1)                                                 H02460
      TIM = TIME1-TIME                                                    H02470
      WRITE (IPR,910) TIME1,TIM                                           H02480
C                                                                         H02490
      RETURN                                                              H02500
C                                                                         H02510
  900 FORMAT ('0 THE TIME AT THE START OF ABINIT IS ',F12.3)              H02520
  905 FORMAT ('0 INITIAL LAYER',I5,/,'0 FILE ',I5,                        H02530
     *        ' INITIALIZED ONTO FILE',I5)                                H02540
  910 FORMAT ('0 THE TIME AT THE END OF ABINIT IS ',F12.3/F12.3,          H02550
     *        ' SECS WERE REQUIRED FOR THIS MERGE ')                      H02560
C                                                                         H02570
      END                                                                 H02580
C
C     ----------------------------------------------------------------
C
      SUBROUTINE ABSMRG (NPTS,LFILE,MFILE,JPATHL)                         H02590
C                                                                         H02600
      IMPLICIT DOUBLE PRECISION (V)                                     ! H02610
C                                                                         H02620
C     SUBROUTINE ABSMRG MERGES ABSORPTION VALUES FROM KFILE INTO LFILE    H02630
C                                                                         H02640
      COMMON R1(2410),OLDR1(2410)                                         H02650
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H02660
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       H02670
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H02680
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H02690
     *              EXTID(10)                                             H02700
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H02710
C                                                                         H02720
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H02730
C                                                                         H02740
      COMMON /ABSHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       H02750
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   H02760
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    H02770
      COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO                                 H02780
      COMMON /ABSPNL/ V1P,V2P,DVP,NLIM,NSHFT,NPNTS                        H02790
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H02800
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        H02810
     *              NLTEFL,LNFIL4,LNGTH4                                  H02820
      DIMENSION SAVOR1(5),A1(10),A2(10),A3(10),A4(10)                     H02830
      DIMENSION WKSAV(35)                                                 H02840
      DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)                             H02850
C                                                                         H02860
      CHARACTER*40 CYID                                                   H02870
C                                                                         H02880
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                   H02890
     *            (OPNLHD(1),V1PO) , (FSCDID(4),IAERSL),                  H02900
     *            (FSCDID(5),IEMIT) , (FSCDID(8),IPATHL),                 H02910
     *            (FSCDID(16),LAYR1)                                      H02920
C                                                                         H02930
      CALL CPUTIM (TIME)                                                  H02940
      NPANLS = 0                                                          H02950
      IF (NOPR.EQ.0) WRITE (IPR,900) TIME                                 H02960
      CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)                            H02970
      DVL = DV                                                            H02980
      LAY1SV = LAYR1                                                      H02990
      PL = PAVE                                                           H03000
      TL = TAVE                                                           H03010
      WTOTL = 0.                                                          H03020
      DO 10 MOL = 1, NMOL                                                 H03030
         WTOTL = WTOTL+WK(MOL)                                            H03040
         WKSAV(MOL) = WK(MOL)                                             H03050
   10 CONTINUE                                                            H03060
      WTOTL = WTOTL+WBROAD                                                H03070
      WN2SAV = WBROAD                                                     H03080
C                                                                         H03090
C     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)                 H03100
C                                                                         H03110
      IF (IAERSL.GT.0) WRITE (CYID,'(5A8)') (YID(I),I=3,7)                H03120
      CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)                            H03130
      IF (IAERSL.GT.0) READ (CYID,'(5A8)') (YID(I),I=3,7)                 H03140
      IF (JPATHL.GE.1) IPATHL = JPATHL                                    H03150
      PLAY = PAVE                                                         H03160
      TLAY = TAVE                                                         H03170
      DVK = DV                                                            H03180
      LAYR1 = LAY1SV                                                      H03190
      FACT = 1.                                                           H03200
      IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.                            H03210
      IF (DVL.EQ.DVK) ITYPE = 0                                           H03220
      IF (DVL.GT.DVK) ITYPE = DVK/(DVL-DVK)+0.5                           H03230
      IF (DVL.LT.DVK) ITYPE = -INT(DVL/(DVK-DVL)+0.5)                     H03240
      SAVOR1(4) = 0.                                                      H03250
C                                                                         H03260
C     ITYPE .LT. 0  IF DV(K-1) IS LESS THAN DV(K)                         H03270
C                                                                         H03280
      IF (ITYPE.LT.0) STOP ' ABSMRG: ITYPE LT 0 '                         H03290
      ITYPE = IABS(ITYPE)                                                 H03300
      WTOTK = 0.                                                          H03310
      DO 20 MOL = 1, NMOL                                                 H03320
         WTOTK = WTOTK+FACT*WK(MOL)                                       H03330
         WK(MOL) = FACT*WK(MOL)+WKSAV(MOL)                                H03340
   20 CONTINUE                                                            H03350
      WTOTK = WTOTK+FACT*WBROAD                                           H03360
      PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)                          H03370
      TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)                          H03380
      WBROAD = FACT*WBROAD+WN2SAV                                         H03390
      SECANT = 0.                                                         H03400
      IF (NOPR.EQ.0) WRITE (IPR,905) LAYR1,LAYER,KFILE,LFILE,MFILE        H03410
      IEMIT = 0                                                           H03420
C                                                                         H03430
C     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES               H03440
C                                                                         H03450
      CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)                                H03460
      DVXM = DV                                                           H03470
      DO 30 K = 1, 5                                                      H03480
         SAVOR1(K) = 0.                                                   H03490
   30 CONTINUE                                                            H03500
      ATYPE = ITYPE                                                       H03510
      AP = 1.0/(ATYPE+1.0)                                                H03520
      IF (ITYPE.NE.0) GO TO 80                                            H03530
C                                                                         H03540
C     1/1 RATIO ONLY                                                      H03550
C                                                                         H03560
   40 CONTINUE                                                            H03570
      CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)                            H03580
      IF (LEOF.LE.0) GO TO 50                                             H03590
      CALL BUFIN (LFILE,LEOF,OLDR1(1),NLIMO)                              H03600
   50 CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)                            H03610
      IF (KEOF.LE.0) GO TO 250                                            H03620
      CALL BUFIN (KFILE,KEOF,R1(1),NLIM)                                  H03630
C                                                                         H03640
      IF (FACT.EQ.1) THEN                                                 H03650
         DO 60 KOD = 1, NLIM                                              H03660
            R1(KOD) = R1(KOD)+OLDR1(KOD)                                  H03670
   60    CONTINUE                                                         H03680
C                                                                         H03690
      ELSE                                                                H03700
         DO 70 KOD = 1, NLIM                                              H03710
            R1(KOD) = R1(KOD)+R1(KOD)+OLDR1(KOD)                          H03720
   70    CONTINUE                                                         H03730
      ENDIF                                                               H03740
C                                                                         H03750
      CALL ABSOUT (V1P,V2P,DVP,NLIM,1,MFILE,NPTS,R1,NPANLS)               H03760
C                                                                         H03770
      GO TO 40                                                            H03780
C                                                                         H03790
C     ALL RATIOS EXCEPT 1/1                                               H03800
C                                                                         H03810
   80 LL = ITYPE+1                                                        H03820
      DO 90 JPG = 1, ITYPE                                                H03830
         APG = JPG                                                        H03840
         P = 1.0-(AP*APG)                                                 H03850
C                                                                         H03860
C    THE FOLLOWING ARE THE CONSTANTS FOR THE LAGRANGE 4 POINT             H03870
C    INTERPOLATION.                                                       H03880
C                                                                         H03890
         A1(JPG) = -P*(P-1.0)*(P-2.0)/6.0                                 H03900
         A2(JPG) = (P**2-1.0)*(P-2.0)*0.5                                 H03910
         A3(JPG) = -P*(P+1.0)*(P-2.0)*0.5                                 H03920
         A4(JPG) = P*(P**2-1.0)/6.0                                       H03930
   90 CONTINUE                                                            H03940
C                                                                         H03950
C     ********  BEGINNING OF LOOP THAT DOES INTERPOLATION  *********      H03960
C                                                                         H03970
      CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)                            H03980
      CALL BUFIN (LFILE,LEOF,OLDR1(1),NLIMO)                              H03990
      MAXLF = NLIMO                                                       H04000
      NVS = 1                                                             H04010
      CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)                            H04020
      CALL BUFIN (KFILE,KEOF,R1(1),NLIM)                                  H04030
      IF (KEOF.LE.0) GO TO 250                                            H04040
      II = 1                                                              H04050
      DIF = DVP*0.01                                                      H04060
      IF (ABS(V1PO-V1P).LT.DIF) GO TO 120                                 H04070
C                                                                         H04080
C     V1P  <  V1PO  LASER OPTION                                          H04090
C                                                                         H04100
  100 NVS = NVS+1                                                         H04110
      V1PN = V1PO+DVPO*(NVS-1)                                            H04120
      IF (ABS(V1PN-V1P).LT.DIF) GO TO 120                                 H04130
      IF (V1PN.LT.V1P) GO TO 100                                          H04140
C                                                                         H04150
  110 II = II+1                                                           H04160
      V1PP = V1P+DVP*(II-1)                                               H04170
      IF (ABS(V1PN-V1PP).LT.DIF) GO TO 120                                H04180
      IF (V1PP.LT.V1PN) GO TO 110                                         H04190
C                                                                         H04200
      GO TO 100                                                           H04210
  120 R1(II) = FACT*R1(II)+OLDR1(NVS)                                     H04220
      V1PN = V1PO+DVPO*(NVS-1)                                            H04230
      V1PP = V1P+DVP*(II-1)                                               H04240
      II = II+1                                                           H04250
  130 JJ = 1                                                              H04260
C                                                                         H04270
      DO 240 JPG = 1, LL                                                  H04280
         IF (JPG.EQ.LL) GO TO 140                                         H04290
         IF (NVS.EQ.1) GO TO 150                                          H04300
         GO TO 170                                                        H04310
  140    IF (FACT.EQ.1.) THEN                                             H04320
            R1(II) = R1(II)+OLDR1(NVS)                                    H04330
         ELSE                                                             H04340
            R1(II) = R1(II)+R1(II)+OLDR1(NVS)                             H04350
         ENDIF                                                            H04360
         V1PN = V1PO+DVPO*(NVS-1)                                         H04370
         V1PP = V1P+DVP*(II-1)                                            H04380
         GO TO 190                                                        H04390
  150    IF (SAVOR1(4).EQ.0.0) GO TO 160                                  H04400
         OLDR1Y = SAVOR1(4)                                               H04410
         GO TO 180                                                        H04420
  160    OLDR1Y = OLDR1(1)                                                H04430
         GO TO 180                                                        H04440
  170    OLDR1Y = OLDR1(NVS-1)                                            H04450
  180    OLDR1I = A1(JJ)*OLDR1Y+A2(JJ)*OLDR1(NVS)+A3(JJ)*OLDR1(NVS+1)+    H04460
     *            A4(JJ)*OLDR1(NVS+2)                                     H04470
         IF (FACT.EQ.1.) THEN                                             H04480
            R1(II) = R1(II)+OLDR1I                                        H04490
         ELSE                                                             H04500
            R1(II) = R1(II)+R1(II)+OLDR1I                                 H04510
         ENDIF                                                            H04520
  190    NVS = NVS+1                                                      H04530
         IF (NVS.LE.MAXLF-2) GO TO 200                                    H04540
         SAVOR1(1) = OLDR1(NVS-1)                                         H04550
         SAVOR1(2) = OLDR1(NVS)                                           H04560
         SAVOR1(3) = OLDR1(NVS+1)                                         H04570
         SAVOR1(4) = OLDR1(NVS-2)                                         H04580
         CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)                         H04590
         IF (LEOF.LE.0) GO TO 210                                         H04600
         MAXLF = NLIMO+3                                                  H04610
         CALL BUFIN (LFILE,LEOF,OLDR1(4),NLIMO)                           H04620
         OLDR1(1) = SAVOR1(1)                                             H04630
         OLDR1(2) = SAVOR1(2)                                             H04640
         OLDR1(3) = SAVOR1(3)                                             H04650
         NVS = 2                                                          H04660
  200    II = II+1                                                        H04670
         JJ = JJ+1                                                        H04680
         IF (II.GT.NLIM) GO TO 230                                        H04690
         GO TO 240                                                        H04700
  210    II = II+1                                                        H04710
         AVRG = (SAVOR1(3)+SAVOR1(2))*0.5                                 H04720
  220    R1(II) = FACT*R1(II)+AVRG                                        H04730
         II = II+1                                                        H04740
         IF (II.LE.NLIM) GO TO 220                                        H04750
C                                                                         H04760
C     WRITE OUTPUT FILE                                                   H04770
C                                                                         H04780
  230    CALL ABSOUT (V1P,V2P,DVP,NLIM,1,MFILE,NPTS,R1,NPANLS)            H04790
C                                                                         H04800
         CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)                         H04810
         IF (KEOF.LE.0) GO TO 250                                         H04820
         CALL BUFIN (KFILE,KEOF,R1(1),NLIM)                               H04830
         II = 1                                                           H04840
  240 CONTINUE                                                            H04850
      NVS = NVS-1                                                         H04860
      GO TO 130                                                           H04870
  250 CONTINUE                                                            H04880
C                                                                         H04890
      CALL CPUTIM (TIME1)                                                 H04900
      TIM = TIME1-TIME                                                    H04910
      IF (NOPR.EQ.0) WRITE (IPR,910) TIME1,TIM                            H04920
      RETURN                                                              H04930
C                                                                         H04940
  900 FORMAT ('0 THE TIME AT THE START OF ABSMRG IS ',F12.3)              H04950
  905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,     H04960
     *        ' MERGED WITH FILE ',I5,' ONTO FILE',I5)                    H04970
  910 FORMAT (' THE TIME AT THE END OF ABSMRG IS',F12.3/F12.3,            H04980
     *        ' SECS. WERE REQUIRED FOR THIS ADDITION')                   H04990
C                                                                         H05000
      END                                                                 H05010
C
C     ----------------------------------------------------------------
C
      SUBROUTINE ABSINT (NPTS,LFILE,MFILE,JPATHL)                         H05020
C                                                                         H05030
      IMPLICIT DOUBLE PRECISION (V)                                     ! H05040
C                                                                         H05050
      COMMON NEWOD(2410),ODLAY(-2:2407)                                   H05060
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H05070
     *              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       H05080
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H05090
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H05100
     *              EXTID(10)                                             H05110
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H05120
C                                                                         H05130
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H05140
C                                                                         H05150
      COMMON /ABSHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       H05160
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   H05170
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    H05180
      COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO                                 H05190
      COMMON /ABSPNL/ V1P,V2P,DVP,NLIM,NSHFT,NPNTS                        H05200
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H05210
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        H05220
     *              NLTEFL,LNFIL4,LNGTH4                                  H05230
C                                                                         H05240
      DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)                             H05250
      DIMENSION A1(0:100),A2(0:100),A3(0:100),A4(0:100)                   H05260
      DIMENSION OLDOD(2),ODLAYR(2)                                        H05270
      DIMENSION WKSAV(35)                                                 H05280
C                                                                         H05290
      CHARACTER*40 CYID                                                   H05300
      REAL NEWOD                                                          H05310
C                                                                         H05320
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                   H05330
     *            (OPNLHD(1),V1PO)                                        H05340
      EQUIVALENCE (NEWOD(1),OLDOD(1)) , (ODLAY(1),ODLAYR(1)),             H05350
     *            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),                 H05360
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 H05370
     *            (FSCDID(16),LAYR1)                                      H05380
C                                                                         H05390
C     ***********************************************************         H05400
C     ****** THIS SUBROUTINE DOES LAYER MERGE FOR OPTICAL  ******         H05410
C     ****** DEPTHS USING FOUR POINT GENERAL INTERPOLATION ******         H05420
C     ***********************************************************         H05430
C                                                                         H05440
      CALL CPUTIM (TIME)                                                  H05450
      WRITE (IPR,900) TIME                                                H05460
      NPANLS = 0                                                          H05470
C                                                                         H05480
      CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)                            H05490
      DVL = DV                                                            H05500
      LAY1SV = LAYR1                                                      H05510
      PL = PAVE                                                           H05520
      TL = TAVE                                                           H05530
      WTOTL = 0.                                                          H05540
      DO 10 MOL = 1, NMOL                                                 H05550
         WTOTL = WTOTL+WK(MOL)                                            H05560
         WKSAV(MOL) = WK(MOL)                                             H05570
   10 CONTINUE                                                            H05580
      WTOTL = WTOTL+WBROAD                                                H05590
      WN2SAV = WBROAD                                                     H05600
C                                                                         H05610
C     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)                 H05620
C                                                                         H05630
      IF (IAERSL.GT.0) WRITE (CYID,'(5A8)') (YID(I),I=3,7)                H05640
      CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)                            H05650
      IF (IAERSL.GT.0) READ (CYID,'(5A8)') (YID(I),I=3,7)                 H05660
      IF (JPATHL.GE.1) IPATHL = JPATHL                                    H05670
      PLAY = PAVE                                                         H05680
      TLAY = TAVE                                                         H05690
      DVK = DV                                                            H05700
      LAYR1 = LAY1SV                                                      H05710
      FACT = 1.                                                           H05720
C                                                                         H05730
C     IF(IPATHL.EQ.2 .AND. IANT.EQ.0) FACT=2.                             H05740
C                                                                         H05750
      IF (IPATHL.EQ.2.AND.IANT.EQ.0) STOP ' ABSINT: FACT=2.  '            H05760
      ATYPE = 9.999E09                                                    H05770
      IF (DVK.EQ.DVL) ATYPE = 0.                                          H05780
      IF (DVL.GT.DVK) ATYPE = DVK/(DVL-DVK)+0.5                           H05790
      IF (DVL.LT.DVK) ATYPE = -DVL/(DVK-DVL)-0.5                          H05800
      IF (ATYPE.GT.0) STOP ' ABSINT; ATYPE GT 0 '                         H05810
      WTOTK = 0.                                                          H05820
      WRITE (IPR,905) LAYR1,LAYER,KFILE,LFILE,MFILE,ATYPE                 H05830
      IEMIT = 0                                                           H05840
      DO 20 MOL = 1, NMOL                                                 H05850
         WTOTK = WTOTK+WK(MOL)*FACT                                       H05860
         WK(MOL) = WK(MOL)*FACT+WKSAV(MOL)                                H05870
   20 CONTINUE                                                            H05880
      WTOTK = WTOTK+WBROAD*FACT                                           H05890
      WBROAD = WBROAD*FACT+WN2SAV                                         H05900
      PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)                          H05910
      TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)                          H05920
      SECANT = 0.                                                         H05930
      DV = DVL                                                            H05940
C                                                                         H05950
C     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES               H05960
C                                                                         H05970
      CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)                                H05980
      DVXM = DV                                                           H05990
C                                                                         H06000
      IF (ATYPE.EQ.0.) THEN                                               H06010
C                                                                         H06020
C     1/1 RATIO ONLY                                                      H06030
C                                                                         H06040
   30    CONTINUE                                                         H06050
C                                                                         H06060
         CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)                         H06070
         IF (KEOF.LE.0) GO TO 90                                          H06080
         CALL BUFIN (KFILE,KEOF,ODLAYR(1),NLIM)                           H06090
C                                                                         H06100
         CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)                         H06110
         CALL BUFIN (LFILE,LEOF,OLDOD(1),NLIMO)                           H06120
C                                                                         H06130
         DO 40 I = 1, NLIM                                                H06140
            NEWOD(I) = ODLAYR(I)+OLDOD(I)                                 H06150
   40    CONTINUE                                                         H06160
         CALL ABSOUT (V1PO,V2PO,DVPO,NLIMO,1,MFILE,NPTS,NEWOD,NPANLS)     H06170
         GO TO 30                                                         H06180
C                                                                         H06190
      ENDIF                                                               H06200
C                                                                         H06210
C     ALL RATIOS EXCEPT 1/1                                               H06220
C                                                                         H06230
      DO 50 JP = 0, 100                                                   H06240
         APG = JP                                                         H06250
         P = 0.01*APG                                                     H06260
C                                                                         H06270
C     THE FOLLOW ARE THE CONSTANTS FOR THE LAGRANGE 4 POINT               H06280
C     INTERPOLATION                                                       H06290
C                                                                         H06300
         A1(JP) = -P*(P-1.0)*(P-2.0)/6.0                                  H06310
         A2(JP) = (P**2-1.0)*(P-2.0)*0.5                                  H06320
         A3(JP) = -P*(P+1.0)*(P-2.0)*0.5                                  H06330
         A4(JP) = P*(P**2-1.0)/6.0                                        H06340
   50 CONTINUE                                                            H06350
C                                                                         H06360
      CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)                            H06370
      IF (KEOF.LE.0) GO TO 90                                             H06380
      CALL BUFIN (KFILE,KEOF,ODLAYR(1),NLIM)                              H06390
C                                                                         H06400
      ODLAY(-2) = ODLAY(1)                                                H06410
      ODLAY(-1) = ODLAY(1)                                                H06420
      ODLAY(0) = ODLAY(1)                                                 H06430
C                                                                         H06440
      RATDV = DVL/DVK                                                     H06450
C                                                                         H06460
   60 CONTINUE                                                            H06470
C                                                                         H06480
      CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)                            H06490
      IF (LEOF.LE.0) GO TO 90                                             H06500
      CALL BUFIN (LFILE,LEOF,OLDOD(1),NLIMO)                              H06510
C                                                                         H06520
C     FJJ IS OFFSET BY 2. FOR ROUNDING PURPOSES                           H06530
C                                                                         H06540
      FJ1DIF = (V1PO-V1P)/DVP+1.+2.                                       H06550
C                                                                         H06560
C     ***** BEGINNING OF LOOP THAT DOES MERGE  *****                      H06570
C                                                                         H06580
      DO 80 II = 1, NLIMO                                                 H06590
                                                                          H06600
C                                                                         H06610
   70    CONTINUE                                                         H06620
C                                                                         H06630
         FJJ = FJ1DIF+RATDV*FLOAT(II-1)                                   H06640
         JJ = IFIX(FJJ)-2                                                 H06650
C                                                                         H06660
         IF (JJ+2.GT.NLIM) THEN                                           H06670
            ODLAY(-2) = ODLAY(NLIM-2)                                     H06680
            ODLAY(-1) = ODLAY(NLIM-1)                                     H06690
            ODLAY(0) = ODLAY(NLIM)                                        H06700
            V1PST = V1P                                                   H06710
            V2PST = V2P                                                   H06720
            NLIMST = NLIM                                                 H06730
C                                                                         H06740
            CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)                      H06750
C                                                                         H06760
            IF (KEOF.LE.0) THEN                                           H06770
               V1P = V1PST                                                H06780
               DVP = DVK                                                  H06790
               V2P = V2PST+2.*DVP                                         H06800
               NLIM = NLIMST+2                                            H06810
               ODLAY(NLIM-1) = ODLAY(NLIM-2)                              H06820
               ODLAY(NLIM) = ODLAY(NLIM-2)                                H06830
            ELSE                                                          H06840
               CALL BUFIN (KFILE,KEOF,ODLAYR(1),NLIM)                     H06850
            ENDIF                                                         H06860
C                                                                         H06870
            FJ1DIF = (V1PO-V1P)/DVP+1.+2.                                 H06880
            GO TO 70                                                      H06890
         ENDIF                                                            H06900
C                                                                         H06910
C     JP = (FJJ-FLOAT(JJ))*100. + 0.5 - 200.                              H06920
C                                                                         H06930
         JP = (FJJ-FLOAT(JJ))*100.-199.5                                  H06940
         IF (JP.GT.100) THEN                                              H06950
            WRITE (IPR,910) JP,JJ,NLIM                                    H06960
            STOP                                                          H06970
         ENDIF                                                            H06980
C                                                                         H06990
C     INTERPOLATE THE OLD TRANSMISSION                                    H07000
C                                                                         H07010
         ODLAYI = A1(JP)*ODLAY(JJ-1)+A2(JP)*ODLAY(JJ)+                    H07020
     *            A3(JP)*ODLAY(JJ+1)+A4(JP)*ODLAY(JJ+2)                   H07030
         IF (ODLAYI.LT.0.) ODLAYI = 0.                                    H07040
C                                                                         H07050
         NEWOD(II) = ODLAYI+OLDOD(II)                                     H07060
C                                                                         H07070
   80 CONTINUE                                                            H07080
C                                                                         H07090
      CALL ABSOUT (V1PO,V2PO,DVPO,NLIMO,1,MFILE,NPTS,NEWOD,NPANLS)        H07100
C                                                                         H07110
      GO TO 60                                                            H07120
C                                                                         H07130
   90 CONTINUE                                                            H07140
C                                                                         H07150
      CALL CPUTIM (TIME1)                                                 H07160
      TIM = TIME1-TIME                                                    H07170
      WRITE (IPR,915) TIME1,TIM                                           H07180
C                                                                         H07190
      RETURN                                                              H07200
C                                                                         H07210
  900 FORMAT ('0 THE TIME AT THE START OF ABSINT IS ',F12.3)              H07220
  905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,     H07230
     *        ' MERGED WITH FILE ',I5,' ONTO FILE',I5,'  WITH XTYPE=',    H07240
     *        G15.5)                                                      H07250
  910 FORMAT ('0 JP, JJ, NLIM ',3I6)                                      H07260
  915 FORMAT ('0 THE TIME AT THE END OF ABSINT IS ',F12.3/F12.3,          H07270
     *        ' SECS WERE REQUIRED FOR THIS MERGE ')                      H07280
C                                                                         H07290
      END                                                                 H07300
C
C     ----------------------------------------------------------------
C
      SUBROUTINE ABSOUT (V1PO,V2PO,DVPO,NLIMO,JLO,MFILE,NPTS,R1,NPANLS)   H07310
C                                                                         H07320
      IMPLICIT DOUBLE PRECISION (V)                                     ! H07330
C                                                                         H07340
C     SUBROUTINE ABSOUT OUPUTS THE MERGED RESULT (R1) ONTO MFILE          H07350
C                                                                         H07360
      COMMON /ABSPNI/ V1P,V2P,DVP,NLIM                                    H07370
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H07380
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        H07390
     *              NLTEFL,LNFIL4,LNGTH4                                  H07400
      DIMENSION PNLHDR(2),R1(*)                                           H07410
C                                                                         H07420
      EQUIVALENCE (PNLHDR(1),V1P)                                         H07430
C                                                                         H07440
      V1P = V1PO                                                          H07450
      V2P = V2PO                                                          H07460
      DVP = DVPO                                                          H07470
      NLIM = NLIMO                                                        H07480
C                                                                         H07490
      NPANLS = NPANLS+1                                                   H07500
      CALL BUFOUT (MFILE,PNLHDR(1),NPHDRF)                                H07510
      CALL BUFOUT (MFILE,R1(JLO),NLIM)                                    H07520
      IF (NPTS.LE.0) GO TO 20                                             H07530
      IF (NPANLS.EQ.1) WRITE (IPR,900)                                    H07540
      WRITE (IPR,905)                                                     H07550
      NNPTS = NPTS                                                        H07560
      IF (NPTS.GT.(NLIM/2)+1) NNPTS = (NLIM/2)+1                          H07570
      JHILIM = JLO+NLIM-NNPTS                                             H07580
      DO 10 I = 1, NNPTS                                                  H07590
         J = JLO+I-1                                                      H07600
         K = JHILIM+I-1                                                   H07610
         VJ = V1P+FLOAT(J-JLO)*DVP                                        H07620
         VK = V1P+FLOAT(K-JLO)*DVP                                        H07630
         WRITE (IPR,910) J,VJ,R1(J),K,VK,R1(K)                            H07640
   10 CONTINUE                                                            H07650
   20 CONTINUE                                                            H07660
C                                                                         H07670
      RETURN                                                              H07680
C                                                                         H07690
  900 FORMAT ('0 ','LOCATION  WAVENUMBER',2X,'OPT DPTH',27X,              H07700
     *        'LOCATION   WAVENUMBER',2X,'OPT DPTH')                      H07710
  905 FORMAT (' ')                                                        H07720
  910 FORMAT (I8,2X,F12.6,1P,E15.7,0P,20X,I8,2X,F12.6,1P,E15.7)           H07730
C                                                                         H07740
      END                                                                 H07750
C
C     ----------------------------------------------------------------
C
      FUNCTION BBFN (VI,DVI,V2I,XKT,VINEW,BBDEL,BBLAST)                   H07760
C                                                                         H07770
      IMPLICIT DOUBLE PRECISION (V)                                     ! H07780
C                                                                         H07790
C     FUNCTION BBFN CALCULATES BLACK BODY FN FOR WAVENUMBER VALUE VI      H07800
C     AND CALCULATES THE WAVENUMBER VALUE (VINEW) FOR NEXT BBFN CALC.     H07810
C                                                                         H07820
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H07830
C                                                                         H07840
C               LAST MODIFICATION:    23 AUGUST 1991                      H07850
C                                                                         H07860
C                  IMPLEMENTATION:    R.D. WORSHAM                        H07870
C                                                                         H07880
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         H07890
C                                     R.D. WORSHAM                        H07900
C                                     J.L. MONCET                         H07910
C                                                                         H07920
C                                                                         H07930
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H07940
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H07950
C                                                                         H07960
C----------------------------------------------------------------------   H07970
C                                                                         H07980
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H07990
C                                     OFFICE OF ENERGY RESEARCH           H08000
C                                     DEPARTMENT OF ENERGY                H08010
C                                                                         H08020
C                                                                         H08030
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             H08040
C                                                                         H08050
C                                             FASCOD3                     H08060
C                                                                         H08070
C                                                                         H08080
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H08090
C                                                                         H08100
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H08110
C                                                                         H08120
      DATA FACTOR / 0.003 /                                               H08130
C                                                                         H08140
         XVI = VI
         XVIOKT = XVI/XKT
         EXPNEG = EXP(-XVIOKT)
         GNU2 = XVI*XVI
         BG2  = XVIOKT*XVIOKT
C
C     IF FIRST CALL, INITIALIZE BBLAST                                    H08150
C                                                                         H08160
      IF (BBLAST.LT.0.) THEN                                              H08170
         IF (XKT.GT.0.0) THEN                                             H08180
C                                                                         H08190
            IF (XVIOKT.LE.0.01) THEN                                      H08230
               BBLAST = RADCN1*(XVI**2)*XKT/(1.+0.5*XVIOKT)               H08240
            ELSEIF (XVIOKT.LE.80.0) THEN                                  H08250
               BBLAST = RADCN1*(XVI**3)/(EXP(XVIOKT)-1.)                  H08260
            ELSE                                                          H08270
               BBLAST = 0.                                                H08280
            ENDIF                                                         H08290
         ELSE                                                             H08300
            BBLAST = 0.                                                   H08310
         ENDIF                                                            H08320
      ENDIF                                                               H08330
C                                                                         H08340
C     SET BBFN EQUAL TO BLACK BODY FUNCTION AT VI                         H08350
C                                                                         H08360
C     BBLAST IS BBFN(VI) FOR EACH SUBSEQUENT CALL                         H08370
C                                                                         H08380
      BBFN = BBLAST                                                       H08390
C                                                                         H08400
      INTVLS = 1                                                          H08410
      DELTAV2 = V2I - VI
      IF (XKT.GT.0.0) THEN                                                H08420
C                                                                         H08430
         IF (XVIOKT.LE.0.01) THEN                                         H08470
            XDELT = (GNU2 * (4.+4.*XVIOKT + BG2))/
     *        (10.*BG2 - 24.*XVIOKT + 8.)
            DELTAV = SQRT(ABS(FACTOR*XDELT))
            IF (DELTAV .GT. DELTAV2) DELTAV = DELTAV2
            INTVLS = (DELTAV)/DVI
            INTVLS = MAX(INTVLS,1)
            VINEW = VI+DVI*FLOAT(INTVLS)
            XVINEW = VINEW                                                H08580
C                                                                         H08590
            BBNEXT = RADCN1*(XVINEW**2)*XKT/(1.+0.5*XVINEW/XKT)           H08600
         ELSEIF (XVIOKT.LE.80.0) THEN                                     H08610
            FRONT  = XVIOKT/(1.-EXPNEG)
            BOX    = 3.- FRONT
            DELT2C = (1./GNU2)*(2.*BOX-FRONT*(1.+BOX-FRONT*EXPNEG))
            DELTAV = SQRT(ABS(FACTOR/DELT2C))
            IF (DELTAV .GT. DELTAV2) DELTAV = DELTAV2
            INTVLS = (DELTAV)/DVI
            INTVLS = MAX(INTVLS,1)
            VINEW = VI+DVI*FLOAT(INTVLS)
            XVINEW = VINEW                                                H08740
C                                                                         H08750
            BBNEXT = RADCN1*(XVINEW**3)/(EXP(XVINEW/XKT)-1.)              H08760
         ELSE                                                             H08770
            BBNEXT = 0.                                                   H08780
            VINEW = 9.0E+9                                                H08790
         ENDIF                                                            H08800
      ELSE                                                                H08810
         BBNEXT = 0.                                                      H08820
         VINEW = 9.0E+9                                                   H08830
      ENDIF                                                               H08840
C                                                                         H08850
      BBDEL = (BBNEXT-BBFN)/FLOAT(INTVLS)                                 H08860
C                                                                         H08870
      VINEW = VINEW-DVI+0.00001                                           H08880
      BBLAST = BBNEXT                                                     H08890
C                                                                         H08900
      RETURN                                                              H08910
C                                                                         H08920
      END                                                                 H08930
C
C     ----------------------------------------------------------------
C
      FUNCTION EMISFN (VI,DVI,VINEM,EMDEL,EMLAST)                         H08940
C                                                                         H08950
      IMPLICIT DOUBLE PRECISION (V)                                     ! H08960
C                                                                         H08970
C     FUNCTION EMISFN CALCULATES BOUNDARY EMISSIVITY FOR WAVE NUMBER      H08980
C     VALUE CORRESPONDING TO VI AND VINEM, AND THEN CALCULATES THE        H08990
C     LINEAR CHANGE BETWEEN THE EMISSIVITY VALUES AT VI AND VINEM         H09000
C                                                                         H09010
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H09020
C                                                                         H09030
C               LAST MODIFICATION:    23 AUGUST 1991                      H09040
C                                                                         H09050
C                  IMPLEMENTATION:    R.D. WORSHAM                        H09060
C                                                                         H09070
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         H09080
C                                     R.D. WORSHAM                        H09090
C                                     J.L. MONCET                         H09100
C                                                                         H09110
C                                                                         H09120
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H09130
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H09140
C                                                                         H09150
C----------------------------------------------------------------------   H09160
C                                                                         H09170
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H09180
C                                     OFFICE OF ENERGY RESEARCH           H09190
C                                     DEPARTMENT OF ENERGY                H09200
C                                                                         H09210
C                                                                         H09220
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             H09230
C                                                                         H09240
C                                             FASCOD3                     H09250
C                                                                         H09260
C                                                                         H09270
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H09280
C                                                                         H09290
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP                   H09300
      EQUIVALENCE (BNDEMI(1),A) , (BNDEMI(2),B) , (BNDEMI(3),C)           H09310
C                                                                         H09320
      DATA FACTOR / 0.001 /                                               H09330
C                                                                         H09340
C     CHECK FOR CONSTANT E (INDEPENDENT OF VI)                            H09350
C     IF CONSTANT RETURN LARGE VALUE FOR VINEM                            H09360
C                                                                         H09370
      IF (B.EQ.0..AND.C.EQ.0.) THEN                                       H09380
         EMISFN = A                                                       H09390
         VINEM = 9.99E+9                                                  H09400
         EMDEL = 0.0                                                      H09410
         EMLAST = EMISFN                                                  H09420
         RETURN                                                           H09430
      ENDIF                                                               H09440
C                                                                         H09450
      XVI = VI                                                            H09460
      IF (EMLAST.LT.0.) THEN                                              H09470
         EMLAST = A+B*XVI+C*XVI*XVI                                       H09480
      ENDIF                                                               H09490
C                                                                         H09500
C     SET EMISFN EQUAL TO EMISSIVITY AT VI                                H09510
C                                                                         H09520
C     EMLAST IS EMISFN(VI) FOR EACH SUBSEQUENT CALL                       H09530
C                                                                         H09540
      EMISFN = EMLAST                                                     H09550
C                                                                         H09560
      IF (VINEM.GE.0.0) THEN                                              H09570
         XVNEXT = XVI+FACTOR/(B+2.*C*XVI)                                 H09580
         INTVLS = (XVNEXT-XVI)/DVI                                        H09590
         INTVLS = MAX(INTVLS,1)                                           H09600
         XVNEXT = XVI+DVI*FLOAT(INTVLS)                                   H09610
      ELSE                                                                H09620
         XVNEXT = ABS(VINEM)+DVI-0.00001                                  H09630
         INTVLS = (XVNEXT-XVI)/DVI                                        H09640
         INTVLS = MAX(INTVLS,1)                                           H09650
      ENDIF                                                               H09660
C                                                                         H09670
      EMNEXT = A+B*XVNEXT+C*XVNEXT*XVNEXT                                 H09680
C                                                                         H09690
      EMDEL = (EMNEXT-EMISFN)/FLOAT(INTVLS)                               H09700
C                                                                         H09710
      VINEM = XVNEXT-DVI+0.00001                                          H09720
      EMLAST = EMNEXT                                                     H09730
C                                                                         H09740
      RETURN                                                              H09750
C                                                                         H09760
      END                                                                 H09770
C
C     ----------------------------------------------------------------
C
      FUNCTION REFLFN (VI,DVI,VINRF,RFDEL,RFLAST)                         H09780
C                                                                         H09790
      IMPLICIT DOUBLE PRECISION (V)                                     ! H09800
C                                                                         H09810
C     FUNCTION REFLFN CALCULATES BOUNDARY REFLECTIVITY FOR WAVE NUMBER    H09820
C     VALUE CORRESPONDING TO VI AND VINRF, AND THEN CALCULATES THE        H09830
C     LINEAR CHANGE BETWEEN THE REFLECTIVITY VALUES AT VI AND VINRF       H09840
C                                                                         H09850
C                                                                         H09860
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H09870
C                                                                         H09880
C               LAST MODIFICATION:    23 AUGUST 1991                      H09890
C                                                                         H09900
C                  IMPLEMENTATION:    R.D. WORSHAM                        H09910
C                                                                         H09920
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         H09930
C                                     R.D. WORSHAM                        H09940
C                                     J.L. MONCET                         H09950
C                                                                         H09960
C                                                                         H09970
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H09980
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H09990
C                                                                         H10000
C----------------------------------------------------------------------   H10010
C                                                                         H10020
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H10030
C                                     OFFICE OF ENERGY RESEARCH           H10040
C                                     DEPARTMENT OF ENERGY                H10050
C                                                                         H10060
C                                                                         H10070
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             H10080
C                                                                         H10090
C                                             FASCOD3                     H10100
C                                                                         H10110
C                                                                         H10120
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H10130
C                                                                         H10140
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP                   H10150
      EQUIVALENCE (BNDRFL(1),A) , (BNDRFL(2),B) , (BNDRFL(3),C)           H10160
C                                                                         H10170
      DATA FACTOR / 0.001 /                                               H10180
C                                                                         H10190
C     CHECK FOR CONSTANT R (INDEPENDENT OF VI)                            H10200
C     IF CONSTANT RETURN LARGE VALUE FOR VINRF                            H10210
C                                                                         H10220
      IF (B.EQ.0..AND.C.EQ.0.) THEN                                       H10230
         REFLFN = A                                                       H10240
         VINRF = 9.99E+9                                                  H10250
         RFDEL = 0.0                                                      H10260
         RFLAST = REFLFN                                                  H10270
         RETURN                                                           H10280
      ENDIF                                                               H10290
C                                                                         H10300
      XVI = VI                                                            H10310
      IF (RFLAST.LT.0.) THEN                                              H10320
         RFLAST = A+B*XVI+C*XVI*XVI                                       H10330
      ENDIF                                                               H10340
C                                                                         H10350
C     SET REFLFN EQUAL TO REFLECTIVITY AT VI                              H10360
C                                                                         H10370
C     RFLAST IS REFLFN(VI) FOR EACH SUBSEQUENT CALL                       H10380
C                                                                         H10390
      REFLFN = RFLAST                                                     H10400
C                                                                         H10410
      IF (VINRF.GE.0.0) THEN                                              H10420
         XVNEXT = XVI+FACTOR/(B+2.*C*XVI)                                 H10430
         INTVLS = (XVNEXT-XVI)/DVI                                        H10440
         INTVLS = MAX(INTVLS,1)                                           H10450
         XVNEXT = XVI+DVI*FLOAT(INTVLS)                                   H10460
      ELSE                                                                H10470
         XVNEXT = ABS(VINRF)+DVI-0.00001                                  H10480
         INTVLS = (XVNEXT-XVI)/DVI                                        H10490
         INTVLS = MAX(INTVLS,1)                                           H10500
      ENDIF                                                               H10510
C                                                                         H10520
      RFNEXT = A+B*XVNEXT+C*XVNEXT*XVNEXT                                 H10530
C                                                                         H10540
      RFDEL = (RFNEXT-REFLFN)/FLOAT(INTVLS)                               H10550
C                                                                         H10560
      VINRF = XVNEXT-DVI+0.00001                                          H10570
      RFLAST = RFNEXT                                                     H10580
C                                                                         H10590
      RETURN                                                              H10600
C                                                                         H10610
      END                                                                 H10620
C
C     ----------------------------------------------------------------
C
      SUBROUTINE EMIN (V1P,V2P,DVP,NLIM,KFILE,EM,EMB,TR,KEOF,NPANLS)      H10630
C                                                                         H10640
      IMPLICIT DOUBLE PRECISION (V)                                     ! H10650
C                                                                         H10660
C     SUBROUTINE EMIN INPUTS OPTICAL DEPTH VALUES FROM KFILE AND          H10670
C       CALCULATES SOURCE FUNCTION FOR THE LAYER.                         H10680
C       THIS VERSION WORKS FOR AEROSOLS AND NLTE.                         H10690
C                                                                         H10700
C                                                                         H10710
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H10720
C                                                                         H10730
C               LAST MODIFICATION:    14 AUGUST 1991                      H10740
C                                                                         H10750
C                  IMPLEMENTATION:    R.D. WORSHAM                        H10760
C                                                                         H10770
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         H10780
C                                     R.D. WORSHAM                        H10790
C                                     J.L. MONCET                         H10800
C                                                                         H10810
C                                                                         H10820
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H10830
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H10840
C                                                                         H10850
C----------------------------------------------------------------------   H10860
C                                                                         H10870
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H10880
C                                     OFFICE OF ENERGY RESEARCH           H10890
C                                     DEPARTMENT OF ENERGY                H10900
C                                                                         H10910
C                                                                         H10920
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             H10930
C                                                                         H10940
C                                             FASCOD3                     H10950
C                                                                         H10960
C                                                                         H10970
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H10980
C                                                                         H10990
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H11000
C                                                                         H11010
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       H11020
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   H11030
     *                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    H11040
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H11050
     *              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        H11060
     *              NLTEFL,LNFIL4,LNGTH4                                  H11070
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           H11080
      COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF                            H11090
      COMMON /RMRG/ XKT,XKTA,XKTB,SECNT                                   H11100
C                                                                         H11110
      DIMENSION PNLHDR(2),EM(*),EMB(*),TR(*)                              H11120
C                                                                         H11130
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4)                 H11140
      EQUIVALENCE (PNLHDR(1),V1PBF)                                       H11150
      EQUIVALENCE (FSCDID(4),IAERSL)                                      H11160
C                                                                         H11170
      CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)                            H11180
      IF (KEOF.LE.0) RETURN                                               H11190
      CALL BUFIN (KFILE,KEOF,TR(1),NLIMBF)                                H11200
C                                                                         H11210
C     TR CONTAINS THE OPTICAL DEPTHS AT THIS STAGE                        H11220
C                                                                         H11230
      IF (IHIRAC.EQ.4) CALL BUFIN (KFILE,KEOF,EM(1),NLIMBF)               H11240
C                                                                         H11250
C     EM CONTAINS THE OPTICAL DEPTH CORRECTIONS FOR NLTE AT THIS STAGE    H11260
C                                                                         H11270
      IF (NPANLS.LT.1.AND.IAERSL.EQ.0) WRITE (IPR,900)                    H11280
      IF (NPANLS.LT.1.AND.IAERSL.NE.0) WRITE (IPR,905)                    H11290
C                                                                         H11300
      EXT = 0.                                                            H11310
      ADEL = 0.                                                           H11320
      RADFN0 = 0.                                                         H11330
      RDEL = 0.                                                           H11340
      BB = 0.                                                             H11350
      BBDEL = 0.                                                          H11360
      BBA = 0.                                                            H11370
      BBDLA = 0.                                                          H11380
      BBB = 0.                                                            H11390
      BBDLB = 0.                                                          H11400
C                                                                         H11410
      V1P = V1PBF                                                         H11420
      V2P = V2PBF                                                         H11430
      DVP = DVPBF                                                         H11440
      NLIM = NLIMBF                                                       H11450
      VI = V1P-DVP                                                        H11460
      VIDV = VI                                                           H11470
      VIBB = VI                                                           H11480
      VAER = VI                                                           H11490
      VDUM = VI                                                           H11500
      BBLAST = -1.                                                        H11510
      BBLXTA = -2.                                                        H11520
      BBLXTB = -3.                                                        H11530
      RDLAST = -1.                                                        H11540
      BBDUM = -4.                                                         H11550
      RDDUM = -1.                                                         H11560
      NLIM1 = 0                                                           H11570
      NLIM2 = 0                                                           H11580
C                                                                         H11590
      AA = 0.2                                                            H11600
C                                                                         H11610
      IF (IAERSL.NE.0) THEN                                               H11630
         BB = BBFN(VI,DVP,V2P,XKT,VIBB,BBDEL,BBDUM)                       
         RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDDUM)                     H11640
         EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)              H11650
         IAFBB = 0                                                        H11660
         IF (VITST.LT.VAER.AND.VITST.LT.VIBB) IAFBB = 1                   H11670
         IF (VAER.LT.VITST.AND.VAER.LT.VIBB) IAFBB = 2                    H11680
      ELSE                                                                H11690
         IAFBB = -1                                                       H11700
      ENDIF                                                               H11710
C                                                                         H11720
C     - THIS SECTION TREATS THE CASE WHERE THE LAYER CONTRIBUTES          H11730
C       TO THE RADIATIVE TRANSFER ONLY ONCE                               H11740
C                                                                         H11750
C     - WITH XKTA=0 THIS ALGORITHM REVERTS TO THE ORIGINAL                H11760
C                                                                         H11770
      IF (XKTB.LE.0.) THEN                                                H11780
C                                                                         H11790
C     - THIS SECTION TREATS THE LTE CASE                                  H11800
C                                                                         H11810
         IF (IHIRAC.NE.4) THEN                                            H11820
C                                                                         H11830
   10       NLIM1 = NLIM2+1                                               H11840
C                                                                         H11850
            VI = V1P+FLOAT(NLIM1-1)*DVP                                   H11860
            IF (IAFBB.EQ.-1) THEN                                         H11870
               BB = BBFN(VI,DVP,V2P,XKT,VIDV,BBDEL,BBLAST)                H11880
               IF (XKTA.GT.0.) THEN                                       H11890
                  VIBB = -VIDV                                            H11900
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H11910
               ELSE                                                       H11920
                  BBA = BB                                                H11930
                  BBDLA = BBDEL                                           H11940
               ENDIF                                                      H11950
               BB = BB-BBDEL                                              H11960
               BBA = BBA-BBDLA                                            H11970
            ELSEIF (IAFBB.EQ.0) THEN                                      H11980
               BB = BBFN(VI,DVP,V2P,XKT,VIDV,BBDEL,BBLAST)                H11990
               IF (XKTA.GT.0.) THEN                                       H12000
                  VIBB = -VIDV                                            H12010
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H12020
               ELSE                                                       H12030
                  BBA = BB                                                H12040
                  BBDLA = BBDEL                                           H12050
               ENDIF                                                      H12060
               BB = BB-BBDEL                                              H12070
               BBA = BBA-BBDLA                                            H12080
               VITST = -VIDV                                              H12090
               RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)              H12100
               RADFN0 = RADFN0-RDEL                                       H12110
               VAER = -VIDV                                               H12120
               EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)        H12130
               EXT = EXT-ADEL                                             H12140
            ELSEIF (IAFBB.EQ.1) THEN                                      H12150
               RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)               H12160
               RADFN0 = RADFN0-RDEL                                       H12170
               VIBB = -VIDV                                               H12180
               BB = BBFN(VI,DVP,V2P,XKT,VIBB,BBDEL,BBLAST)                H12190
               IF (XKTA.GT.0.) THEN                                       H12200
                  VIBB = -VIDV                                            H12210
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H12220
               ELSE                                                       H12230
                  BBA = BB                                                H12240
                  BBDLA = BBDEL                                           H12250
               ENDIF                                                      H12260
               BB = BB-BBDEL                                              H12270
               BBA = BBA-BBDLA                                            H12280
               VAER = -VIDV                                               H12290
               EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)        H12300
               EXT = EXT-ADEL                                             H12310
            ELSEIF (IAFBB.EQ.2) THEN                                      H12320
               EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)        H12330
               EXT = EXT-ADEL                                             H12340
               VIBB = -VIDV                                               H12350
               BB = BBFN(VI,DVP,V2P,XKT,VIBB,BBDEL,BBLAST)                H12360
               IF (XKTA.GT.0.) THEN                                       H12370
                  VIBB = -VIDV                                            H12380
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H12390
               ELSE                                                       H12400
                  BBA = BB                                                H12410
                  BBDLA = BBDEL                                           H12420
               ENDIF                                                      H12430
               BB = BB-BBDEL                                              H12440
               BBA = BBA-BBDLA                                            H12450
               VITST = -VIDV                                              H12460
               RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)              H12470
               RADFN0 = RADFN0-RDEL                                       H12480
            ENDIF                                                         H12490
C                                                                         H12500
            NLIM2 = (VIDV-V1P)/DVP+1.001                                  H12510
            NLIM2 = MIN(NLIM2,NLIM)                                       H12520
C                                                                         H12530
            DO 20 I = NLIM1, NLIM2                                        H12540
               EXT = EXT+ADEL                                             H12550
               RADFN0 = RADFN0+RDEL                                       H12560
               ODVI = TR(I)+EXT*RADFN0                                    H12570
               BB = BB+BBDEL                                              H12580
               BBA = BBA+BBDLA                                            H12590
C                                                                         H12600
               XX = AA*ODVI                                               H12610
C                                                                         H12620
               TR(I) = EXP(-ODVI)                                         H12630
               EM(I) = (1.-TR(I))*(BB+XX*BBA)/(1.+XX)                     H12640
C                                                                         H12650
   20       CONTINUE                                                      H12660
C                                                                         H12670
            IF (NLIM2.LT.NLIM) GO TO 10                                   H12680
         ELSE                                                             H12690
C                                                                         H12700
C     - THIS SECTION TREATS THE NLTE CASE                                 H12710
C                                                                         H12720
   30       NLIM1 = NLIM2+1                                               H12730
C                                                                         H12740
            VI = V1P+FLOAT(NLIM1-1)*DVP                                   H12750
            IF (IAFBB.EQ.-1) THEN                                         H12760
               BB = BBFN(VI,DVP,V2P,XKT,VIDV,BBDEL,BBLAST)                H12770
               IF (XKTA.GT.0.) THEN                                       H12780
                  VIBB = -VIDV                                            H12790
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H12800
               ELSE                                                       H12810
                  BBA = BB                                                H12820
                  BBDLA = BBDEL                                           H12830
               ENDIF                                                      H12840
               BB = BB-BBDEL                                              H12850
               BBA = BBA-BBDLA                                            H12860
            ELSEIF (IAFBB.EQ.0) THEN                                      H12870
               BB = BBFN(VI,DVP,V2P,XKT,VIDV,BBDEL,BBLAST)                H12880
               IF (XKTA.GT.0.) THEN                                       H12890
                  VIBB = -VIDV                                            H12900
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H12910
               ELSE                                                       H12920
                  BBA = BB                                                H12930
                  BBDLA = BBDEL                                           H12940
               ENDIF                                                      H12950
               BB = BB-BBDEL                                              H12960
               BBA = BBA-BBDLA                                            H12970
               VITST = -VIDV                                              H12980
               RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)              H12990
               RADFN0 = RADFN0-RDEL                                       H13000
               VAER = -VIDV                                               H13010
               EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)        H13020
               EXT = EXT-ADEL                                             H13030
            ELSEIF (IAFBB.EQ.1) THEN                                      H13040
               RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)               H13050
               RADFN0 = RADFN0-RDEL                                       H13060
               VIBB = -VIDV                                               H13070
               BB = BBFN(VI,DVP,V2P,XKT,VIBB,BBDEL,BBLAST)                H13080
               IF (XKTA.GT.0.) THEN                                       H13090
                  VIBB = -VIDV                                            H13100
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H13110
               ELSE                                                       H13120
                  BBA = BB                                                H13130
                  BBDLA = BBDEL                                           H13140
               ENDIF                                                      H13150
               BB = BB-BBDEL                                              H13160
               BBA = BBA-BBDLA                                            H13170
               VAER = -VIDV                                               H13180
               EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)        H13190
               EXT = EXT-ADEL                                             H13200
            ELSEIF (IAFBB.EQ.2) THEN                                      H13210
               EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)        H13220
               EXT = EXT-ADEL                                             H13230
               VIBB = -VIDV                                               H13240
               BB = BBFN(VI,DVP,V2P,XKT,VIBB,BBDEL,BBLAST)                H13250
               IF (XKTA.GT.0.) THEN                                       H13260
                  VIBB = -VIDV                                            H13270
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H13280
               ELSE                                                       H13290
                  BBA = BB                                                H13300
                  BBDLA = BBDEL                                           H13310
               ENDIF                                                      H13320
               BB = BB-BBDEL                                              H13330
               BBA = BBA-BBDLA                                            H13340
               VITST = -VIDV                                              H13350
               RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)              H13360
               RADFN0 = RADFN0-RDEL                                       H13370
            ENDIF                                                         H13380
C                                                                         H13390
            NLIM2 = (VIDV-V1P)/DVP+1.001                                  H13400
            NLIM2 = MIN(NLIM2,NLIM)                                       H13410
C                                                                         H13420
            DO 40 I = NLIM1, NLIM2                                        H13430
               EXT = EXT+ADEL                                             H13440
               RADFN0 = RADFN0+RDEL                                       H13450
               ODVI = TR(I)+EXT*RADFN0                                    H13460
               BB = BB+BBDEL                                              H13470
               BBA = BBA+BBDLA                                            H13480
C                                                                         H13490
               XX = AA*ODVI                                               H13500
C                                                                         H13510
               TR(I) = EXP(-ODVI)                                         H13520
               EM(I) = (1.-TR(I))*(1.0-EM(I)/ODVI)*(BB+XX*BBA)/(1.+XX)    H13530
C                                                                         H13540
   40       CONTINUE                                                      H13550
C                                                                         H13560
            IF (NLIM2.LT.NLIM) GO TO 30                                   H13570
C                                                                         H13580
         ENDIF                                                            H13590
      ELSE                                                                H13600
C                                                                         H13610
C     - THIS SECTION TREATS THE CASE WHERE THE LAYER CONTRIBUTES          H13620
C       TO THE RADIATIVE TRANSFER TWICE:                                  H13630
C                                                                         H13640
C     - FOR TANGENT PATHS AND FOR THE CASE OF THE REFLECTED ATMOSPHERE    H13650
C                                                                         H13660
         IF (IHIRAC.NE.4) THEN                                            H13670
C                                                                         H13680
C     - THIS SECTION TREATS THE LTE CASE                                  H13690
C                                                                         H13700
   50       NLIM1 = NLIM2+1                                               H13710
C                                                                         H13720
            VI = V1P+FLOAT(NLIM1-1)*DVP                                   H13730
            IF (IAFBB.EQ.-1) THEN                                         H13740
               BB = BBFN(VI,DVP,V2P,XKT,VIDV,BBDEL,BBLAST)                H13750
               IF (XKTA.GT.0.) THEN                                       H13760
                  VIBB = -VIDV                                            H13770
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H13780
               ELSE                                                       H13790
                  BBA = BB                                                H13800
                  BBDLA = BBDEL                                           H13810
               ENDIF                                                      H13820
               IF (XKTB.GT.0.) THEN                                       H13830
                  VIBB = -VIDV                                            H13840
                  BBB = BBFN(VI,DVP,V2P,XKTB,VIBB,BBDLB,BBLXTB)           H13850
               ELSE                                                       H13860
                  BBB = BB                                                H13870
                  BBDLB = BBDEL                                           H13880
               ENDIF                                                      H13890
               BB = BB-BBDEL                                              H13900
               BBA = BBA-BBDLA                                            H13910
               BBB = BBB-BBDLB                                            H13920
            ELSEIF (IAFBB.EQ.0) THEN                                      H13930
               BB = BBFN(VI,DVP,V2P,XKT,VIDV,BBDEL,BBLAST)                H13940
               IF (XKTA.GT.0.) THEN                                       H13950
                  VIBB = -VIDV                                            H13960
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H13970
               ELSE                                                       H13980
                  BBA = BB                                                H13990
                  BBDLA = BBDEL                                           H14000
               ENDIF                                                      H14010
               IF (XKTB.GT.0.) THEN                                       H14020
                  VIBB = -VIDV                                            H14030
                  BBB = BBFN(VI,DVP,V2P,XKTB,VIBB,BBDLB,BBLXTB)           H14040
               ELSE                                                       H14050
                  BBB = BB                                                H14060
                  BBDLB = BBDEL                                           H14070
               ENDIF                                                      H14080
               BB = BB-BBDEL                                              H14090
               BBA = BBA-BBDLA                                            H14100
               BBB = BBB-BBDLB                                            H14110
               VITST = -VIDV                                              H14120
               RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)              H14130
               RADFN0 = RADFN0-RDEL                                       H14140
               VAER = -VIDV                                               H14150
               EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)        H14160
               EXT = EXT-ADEL                                             H14170
            ELSEIF (IAFBB.EQ.1) THEN                                      H14180
               RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)               H14190
               RADFN0 = RADFN0-RDEL                                       H14200
               VIBB = -VIDV                                               H14210
               BB = BBFN(VI,DVP,V2P,XKT,VIBB,BBDEL,BBLAST)                H14220
               IF (XKTA.GT.0.) THEN                                       H14230
                  VIBB = -VIDV                                            H14240
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H14250
               ELSE                                                       H14260
                  BBA = BB                                                H14270
                  BBDLA = BBDEL                                           H14280
               ENDIF                                                      H14290
               IF (XKTB.GT.0.) THEN                                       H14300
                  VIBB = -VIDV                                            H14310
                  BBB = BBFN(VI,DVP,V2P,XKTB,VIBB,BBDLB,BBLXTB)           H14320
               ELSE                                                       H14330
                  BBB = BB                                                H14340
                  BBDLB = BBDEL                                           H14350
               ENDIF                                                      H14360
               BB = BB-BBDEL                                              H14370
               BBA = BBA-BBDLA                                            H14380
               BBB = BBB-BBDLB                                            H14390
               VAER = -VIDV                                               H14400
               EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)        H14410
               EXT = EXT-ADEL                                             H14420
            ELSEIF (IAFBB.EQ.2) THEN                                      H14430
               EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)        H14440
               EXT = EXT-ADEL                                             H14450
               VIBB = -VIDV                                               H14460
               BB = BBFN(VI,DVP,V2P,XKT,VIBB,BBDEL,BBLAST)                H14470
               IF (XKTA.GT.0.) THEN                                       H14480
                  VIBB = -VIDV                                            H14490
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H14500
               ELSE                                                       H14510
                  BBA = BB                                                H14520
                  BBDLA = BBDEL                                           H14530
               ENDIF                                                      H14540
               IF (XKTB.GT.0.) THEN                                       H14550
                  VIBB = -VIDV                                            H14560
                  BBB = BBFN(VI,DVP,V2P,XKTB,VIBB,BBDLB,BBLXTB)           H14570
               ELSE                                                       H14580
                  BBB = BB                                                H14590
                  BBDLB = BBDEL                                           H14600
               ENDIF                                                      H14610
               BB = BB-BBDEL                                              H14620
               BBA = BBA-BBDLA                                            H14630
               BBB = BBB-BBDLB                                            H14640
               VITST = -VIDV                                              H14650
               RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)              H14660
               RADFN0 = RADFN0-RDEL                                       H14670
            ENDIF                                                         H14680
C                                                                         H14690
            NLIM2 = (VIDV-V1P)/DVP+1.001                                  H14700
            NLIM2 = MIN(NLIM2,NLIM)                                       H14710
C                                                                         H14720
            DO 60 I = NLIM1, NLIM2                                        H14730
               EXT = EXT+ADEL                                             H14740
               RADFN0 = RADFN0+RDEL                                       H14750
               ODVI = TR(I)+EXT*RADFN0                                    H14760
               BB = BB+BBDEL                                              H14770
               BBA = BBA+BBDLA                                            H14780
               BBB = BBB+BBDLB                                            H14790
C                                                                         H14800
               XX = AA*ODVI                                               H14810
C                                                                         H14820
               TR(I) = EXP(-ODVI)                                         H14830
               EMX = (1.-TR(I))/(1.+XX)                                   H14840
               EM(I) = EMX*(BB+XX*BBA)                                    H14850
               EMB(I) = EMX*(BB+XX*BBB)                                   H14860
C                                                                         H14870
   60       CONTINUE                                                      H14880
C                                                                         H14890
            IF (NLIM2.LT.NLIM) GO TO 50                                   H14900
C                                                                         H14910
         ELSE                                                             H14920
C                                                                         H14930
C     - THIS SECTION TREATS THE CASE OF NLTE                              H14940
C                                                                         H14950
   70       NLIM1 = NLIM2+1                                               H14960
C                                                                         H14970
            VI = V1P+FLOAT(NLIM1-1)*DVP                                   H14980
            IF (IAFBB.EQ.-1) THEN                                         H14990
               BB = BBFN(VI,DVP,V2P,XKT,VIDV,BBDEL,BBLAST)                H15000
               IF (XKTA.GT.0.) THEN                                       H15010
                  VIBB = -VIDV                                            H15020
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H15030
               ELSE                                                       H15040
                  BBA = BB                                                H15050
                  BBDLA = BBDEL                                           H15060
               ENDIF                                                      H15070
               IF (XKTB.GT.0.) THEN                                       H15080
                  VIBB = -VIDV                                            H15090
                  BBB = BBFN(VI,DVP,V2P,XKTB,VIBB,BBDLB,BBLXTB)           H15100
               ELSE                                                       H15110
                  BBB = BB                                                H15120
                  BBDLB = BBDEL                                           H15130
               ENDIF                                                      H15140
               BB = BB-BBDEL                                              H15150
               BBA = BBA-BBDLA                                            H15160
               BBB = BBB-BBDLB                                            H15170
            ELSEIF (IAFBB.EQ.0) THEN                                      H15180
               BB = BBFN(VI,DVP,V2P,XKT,VIDV,BBDEL,BBLAST)                H15190
               IF (XKTA.GT.0.) THEN                                       H15200
                  VIBB = -VIDV                                            H15210
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H15220
               ELSE                                                       H15230
                  BBA = BB                                                H15240
                  BBDLA = BBDEL                                           H15250
               ENDIF                                                      H15260
               IF (XKTB.GT.0.) THEN                                       H15270
                  VIBB = -VIDV                                            H15280
                  BBB = BBFN(VI,DVP,V2P,XKTB,VIBB,BBDLB,BBLXTB)           H15290
               ELSE                                                       H15300
                  BBB = BB                                                H15310
                  BBDLB = BBDEL                                           H15320
               ENDIF                                                      H15330
               BB = BB-BBDEL                                              H15340
               BBA = BBA-BBDLA                                            H15350
               BBB = BBB-BBDLB                                            H15360
               VITST = -VIDV                                              H15370
               RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)              H15380
               RADFN0 = RADFN0-RDEL                                       H15390
               VAER = -VIDV                                               H15400
               EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)        H15410
               EXT = EXT-ADEL                                             H15420
            ELSEIF (IAFBB.EQ.1) THEN                                      H15430
               RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)               H15440
               RADFN0 = RADFN0-RDEL                                       H15450
               VIBB = -VIDV                                               H15460
               BB = BBFN(VI,DVP,V2P,XKT,VIBB,BBDEL,BBLAST)                H15470
               IF (XKTA.GT.0.) THEN                                       H15480
                  VIBB = -VIDV                                            H15490
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H15500
               ELSE                                                       H15510
                  BBA = BB                                                H15520
                  BBDLA = BBDEL                                           H15530
               ENDIF                                                      H15540
               IF (XKTB.GT.0.) THEN                                       H15550
                  VIBB = -VIDV                                            H15560
                  BBB = BBFN(VI,DVP,V2P,XKTB,VIBB,BBDLB,BBLXTB)           H15570
               ELSE                                                       H15580
                  BBB = BB                                                H15590
                  BBDLB = BBDEL                                           H15600
               ENDIF                                                      H15610
               BB = BB-BBDEL                                              H15620
               BBA = BBA-BBDLA                                            H15630
               BBB = BBB-BBDLB                                            H15640
               VAER = -VIDV                                               H15650
               EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)        H15660
               EXT = EXT-ADEL                                             H15670
            ELSEIF (IAFBB.EQ.2) THEN                                      H15680
               EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)        H15690
               EXT = EXT-ADEL                                             H15700
               VIBB = -VIDV                                               H15710
               BB = BBFN(VI,DVP,V2P,XKT,VIBB,BBDEL,BBLAST)                H15720
               IF (XKTA.GT.0.) THEN                                       H15730
                  VIBB = -VIDV                                            H15740
                  BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)           H15750
               ELSE                                                       H15760
                  BBA = BB                                                H15770
                  BBDLA = BBDEL                                           H15780
               ENDIF                                                      H15790
               IF (XKTB.GT.0.) THEN                                       H15800
                  VIBB = -VIDV                                            H15810
                  BBB = BBFN(VI,DVP,V2P,XKTB,VIBB,BBDLB,BBLXTB)           H15820
               ELSE                                                       H15830
                  BBB = BB                                                H15840
                  BBDLB = BBDEL                                           H15850
               ENDIF                                                      H15860
               BB = BB-BBDEL                                              H15870
               BBA = BBA-BBDLA                                            H15880
               BBB = BBB-BBDLB                                            H15890
               VITST = -VIDV                                              H15900
               RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)              H15910
               RADFN0 = RADFN0-RDEL                                       H15920
            ENDIF                                                         H15930
C                                                                         H15940
            NLIM2 = (VIDV-V1P)/DVP+1.001                                  H15950
            NLIM2 = MIN(NLIM2,NLIM)                                       H15960
C                                                                         H15970
            DO 80 I = NLIM1, NLIM2                                        H15980
               EXT = EXT+ADEL                                             H15990
               RADFN0 = RADFN0+RDEL                                       H16000
               ODVI = TR(I)+EXT*RADFN0                                    H16010
               BB = BB+BBDEL                                              H16020
               BBA = BBA+BBDLA                                            H16030
               BBB = BBB+BBDLB                                            H16040
C                                                                         H16050
               XX = AA*ODVI                                               H16060
C                                                                         H16070
               TR(I) = EXP(-ODVI)                                         H16080
               EMX = (1.-TR(I))*(1.0-EM(I)/ODVI)/(1.+XX)                  H16090
               EM(I) = EMX*(BB+XX*BBA)                                    H16100
               EMB(I) = EMX*(BB+XX*BBB)                                   H16110
C                                                                         H16120
   80       CONTINUE                                                      H16130
C                                                                         H16140
            IF (NLIM2.LT.NLIM) GO TO 70                                   H16150
C                                                                         H16160
         ENDIF                                                            H16170
      ENDIF                                                               H16180
C                                                                         H16190
      RETURN                                                              H16200
C                                                                         H16210
  900 FORMAT ('0EMISSION AND TRANSMISSION  (MOLECULAR) ')                 H16220
  905 FORMAT ('0EMISSION AND TRANSMISSION (AEROSOLS EFFECTS INCLUDED)')   H16230
C                                                                         H16240
      END                                                                 H16250
C
C     ----------------------------------------------------------------
C
      SUBROUTINE EMINIT (NPTS,MFILE,JPATHL,TBND)                          H16260
C                                                                         H16270
      IMPLICIT DOUBLE PRECISION (V)                                     ! H16280
C                                                                         H16290
C                                                                         H16300
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H16310
C                                                                         H16320
C               LAST MODIFICATION:    14 AUGUST 1991                      H16330
C                                                                         H16340
C                  IMPLEMENTATION:    R.D. WORSHAM                        H16350
C                                                                         H16360
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         H16370
C                                     R.D. WORSHAM                        H16380
C                                     J.L. MONCET                         H16390
C                                                                         H16400
C                                                                         H16410
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H16420
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H16430
C                                                                         H16440
C----------------------------------------------------------------------   H16450
C                                                                         H16460
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H16470
C                                     OFFICE OF ENERGY RESEARCH           H16480
C                                     DEPARTMENT OF ENERGY                H16490
C                                                                         H16500
C                                                                         H16510
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             H16520
C                                                                         H16530
C                                             FASCOD3                     H16540
C                                                                         H16550
C                                                                         H16560
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H16570
C                                                                         H16580
      COMMON NEWEM(2410),NEWTR(2410)                                      H16590
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H16600
     *              AVFIX,LAYER,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,        H16610
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H16620
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H16630
     *              EXTID(10)                                             H16640
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H16650
C                                                                         H16660
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H16670
C                                                                         H16680
      COMMON /EMIHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       H16690
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   H16700
     *                EMISIV,FSCDID(17),NMOL,LAYDUM,YI1,YID(10),LSTWDF    H16710
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP                   H16720
      COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO                                 H16730
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H16740
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,        H16750
     *              NLTEFL,LNFIL4,LNGTH4                                  H16760
      COMMON /RMRG/ XKT,XKTA,XKTB,SECNT                                   H16770
C                                                                         H16780
      CHARACTER*40 CEXT,CYID                                              H16790
C                                                                         H16800
      DIMENSION EMLAYB(2410)                                              H16810
      DIMENSION XFILHD(2),OPNLHD(2)                                       H16820
      DIMENSION EMLAYR(2),TRLAYR(2)                                       H16830
C                                                                         H16840
      EQUIVALENCE (XFILHD(1),XID(1)) , (OPNLHD(1),V1PO)                   H16850
      EQUIVALENCE (NEWEM(1),EMLAYR(1)) , (NEWTR(1),TRLAYR(1)),            H16860
     *            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),                 H16870
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 H16880
     *            (FSCDID(16),LAYR1)                                      H16890
C                                                                         H16900
      REAL NEWEM,NEWTR                                                    H16910
C                                                                         H16920
C                                                                         H16930
C *********************************************************************   H16940
C ****  THIS SUBROUTINE COMPUTES THE EMISSION FOR THE FIRST LAYER  ****   H16950
C *********************************************************************   H16960
C                                                                         H16970
C     TBND IS THE BOUNDARY BLACK BODY TEMPERATUE                          H16980
C                                                                         H16990
C     IPATHL =-1 IS FOR THE LOOKING DOWN CASE WITH REFLECTED ATMOSPHERE   H17000
C     IPATHL = 0 IS FOR THE HORIZONTAL PATH CASE (HOMOGENEOUS LAYER)      H17010
C     IPATHL = 1 IS FOR THE LOOKING DOWN CASE (TO DENSER LAYERS)          H17020
C     IPATHL = 2 IS FOR THE SYMMETRIC TANGENT PATH CASE                   H17030
C     IPATHL = 3 IS FOR THE LOOKING UP CASE (TO LESS DENSE LAYERS         H17040
C                                                                         H17050
      CALL CPUTIM (TIME)                                                  H17060
C                                                                         H17070
C      ** NOTE ON IPATHL =2                                               H17080
C           THE TANGENT MERGE ROUTINES ARE DIVIDED INTO ANTERIOR (1ST)    H17090
C           AND POSTERIOR (2ND) LAYER CROSSINGS.  THIS RECOGNITION IS     H17100
C           TRIGGERED BY THE VALUE OF "IANT"                              H17110
C                                                                         H17120
C          IF  IANT.EQ. 1  THEN POSTERIOR MERGE                           H17130
C          IF  IANT.EQ. 0  THEN NORMAL MERGE                              H17140
C          IF  IANT.EQ.-1  THEN ANTERIOR MERGE                            H17150
C                                                                         H17160
      WRITE (IPR,900) TIME                                                H17170
      NPANLS = 0                                                          H17180
      CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)                            H17190
      IF (JPATHL.GE.1) IPATHL = JPATHL                                    H17200
      PLAY = PAVE                                                         H17210
      TLAY = TAVE                                                         H17220
C                                                                         H17230
C     FOR AEROSOL RUNS, MOVE EXTID INTO YID                               H17240
C                                                                         H17250
      IF (IAERSL.GT.0) THEN                                               H17260
         WRITE (CEXT,'(10A4)') EXTID                                      H17270
         WRITE (CYID,'(5A8)') (YID(I),I=3,7)                              H17280
         CYID(19:40) = CEXT(19:40)                                        H17290
         READ (CYID,'(5A8)') (YID(I),I=3,7)                               H17300
      ENDIF                                                               H17310
C                                                                         H17320
C     IF BOUNDARY PROPERTIES ARE SUPPLIED, AND DOWNWARD LOOKING           H17330
C     CASE; SET IPATHL TO REFLECTED ATMOSPHERE CASE                       H17340
C                                                                         H17350
      IF (IBPROP.EQ.1.AND.IPATHL.EQ.1) IPATHL = -1                        H17360
      IEMIT = 1                                                           H17370
      FACT = 1.                                                           H17380
      TIMEM = 0.0                                                         H17390
      IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.                            H17400
      DO 10 MOL = 1, NMOL                                                 H17410
         WK(MOL) = WK(MOL)*FACT                                           H17420
   10 CONTINUE                                                            H17430
      WBROAD = WBROAD*FACT                                                H17440
      LAYR1 = LAYER                                                       H17450
      WRITE (IPR,905) LAYR1,LAYER,KFILE,MFILE                             H17460
      CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)                                H17470
      DVXM = DV                                                           H17480
      XKT = TAVE/RADCN2                                                   H17490
      XKTBND = TBND/RADCN2                                                H17500
      IF (IPATHL.EQ.-1) THEN                                              H17510
         XKTA = TZU/RADCN2                                                H17520
         XKTB = TZL/RADCN2                                                H17530
      ENDIF                                                               H17540
      IF (IPATHL.EQ.0) THEN                                               H17550
         XKTA = 0.                                                        H17560
         XKTB = 0.                                                        H17570
      ENDIF                                                               H17580
      IF (IPATHL.EQ.1) THEN                                               H17590
         XKTA = TZU/RADCN2                                                H17600
         XKTB = 0.                                                        H17610
      ENDIF                                                               H17620
      IF (IPATHL.EQ.2) THEN                                               H17630
         XKTA = TZU/RADCN2                                                H17640
         XKTB = TZL/RADCN2                                                H17650
      ENDIF                                                               H17660
      IF (IPATHL.EQ.3) THEN                                               H17670
         XKTA = TZL/RADCN2                                                H17680
         XKTB = 0.                                                        H17690
      ENDIF                                                               H17700
      WRITE (IPR,910) IPATHL,IANT                                         H17710
C                                                                         H17720
   20 CONTINUE                                                            H17730
C                                                                         H17740
      CALL CPUTIM (TIMEM1)                                                H17750
      CALL EMIN (V1PO,V2PO,DVPO,NLIMO,KFILE,EMLAYR,EMLAYB,TRLAYR,KEOF,    H17760
     *           NPANLS)                                                  H17770
      CALL CPUTIM (TIMEM2)                                                H17780
      TIMEM = TIMEM+TIMEM2-TIMEM1                                         H17790
      IF (KEOF.LE.0) GO TO 80                                             H17800
      VI = V1PO-DVPO                                                      H17810
      VIDVBD = VI                                                         H17820
      VIDVEM = VI                                                         H17830
      VIDVRF = VI                                                         H17840
      BBLAST = -1.                                                        H17850
      EMLAST = -1.                                                        H17860
      RFLAST = -1.                                                        H17870
      IF (IPATHL.EQ.2.AND.IANT.EQ.0) THEN                                 H17880
         DO 30 J = 1, NLIMO                                               H17890
            TRJ = TRLAYR(J)                                               H17900
            NEWEM(J) = EMLAYR(J)+EMLAYB(J)*TRJ                            H17910
            TRLAYR(J) = TRLAYR(J)*TRJ                                     H17920
   30    CONTINUE                                                         H17930
      ELSEIF ((IPATHL.EQ.1).AND.(TBND.GT.0.)) THEN                        H17940
C                                                                         H17950
         NLIM1 = 0                                                        H17960
         NLIM2 = 0                                                        H17970
         EMDUM = 0.                                                       H17980
         BBDUM = 0.                                                       H17990
         EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMDUM)                      H18000
         BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDVBD,BBDEL,BBDUM)                H18010
         IEMBB = 0                                                        H18020
         IF (VIDVBD.GT.VIDVEM) IEMBB = 1                                  H18030
C                                                                         H18040
   40    NLIM1 = NLIM2+1                                                  H18050
C                                                                         H18060
         VI = V1PO+FLOAT(NLIM1-1)*DVPO                                    H18070
         IF (IEMBB.EQ.0) THEN                                             H18080
            BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDV,BBDEL,BBLAST)              H18090
            BB = BB-BBDEL                                                 H18100
            VIDVEM = -VIDV                                                H18110
            EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMLAST)                  H18120
            EMISIV = EMISIV-EMDEL                                         H18130
         ELSE                                                             H18140
            EMISIV = EMISFN(VI,DVPO,VIDV,EMDEL,EMLAST)                    H18150
            EMISIV = EMISIV-EMDEL                                         H18160
            VIDVBD = -VIDV                                                H18170
            BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDVBD,BBDEL,BBLAST)            H18180
            BB = BB-BBDEL                                                 H18190
         ENDIF                                                            H18200
C                                                                         H18210
         IF (VIDV.GE.9.E+4) THEN 
            NLIM2 = NLIMO+1
         ELSE
            NLIM2 = (VIDV-V1PO)/DVPO+1.001                                H18220
         ENDIF
         NLIM2 = MIN(NLIM2,NLIMO)                                         H18230
C                                                                         H18240
         DO 50 J = NLIM1, NLIM2                                           H18250
            EMISIV = EMISIV+EMDEL                                         H18260
            BB = BB+BBDEL                                                 H18270
            NEWEM(J) = EMLAYR(J)+TRLAYR(J)*BB*EMISIV                      H18280
   50    CONTINUE                                                         H18290
C                                                                         H18300
         IF (NLIM2.LT.NLIMO) GO TO 40                                     H18310
C                                                                         H18320
      ELSEIF ((IPATHL.EQ.-1).AND.(TBND.GT.0.)) THEN                       H18330
C                                                                         H18340
         NLIM1 = 0                                                        H18350
         NLIM2 = 0                                                        H18360
         EMDUM = 0.                                                       H18370
         RFDUM = 0.                                                       H18380
         BBDUM = 0.                                                       H18390
         EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMDUM)                      H18400
         REFLCT = REFLFN(VI,DVPO,VIDVRF,RFDEL,RFDUM)                      H18410
         BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDVBD,BBDEL,BBDUM)                H18420
         IEMBB = 0                                                        H18430
         IF (VIDVEM.LT.VIDVRF.AND.VIDVEM.LT.VIDVBD) IEMBB = 1             H18440
         IF (VIDVRF.LT.VIDVEM.AND.VIDVRF.LT.VIDVBD) IEMBB = 2             H18450
C                                                                         H18460
   60    NLIM1 = NLIM2+1                                                  H18470
C                                                                         H18480
         VI = V1PO+FLOAT(NLIM1-1)*DVPO                                    H18490
         IF (IEMBB.EQ.0) THEN                                             H18500
            BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDV,BBDEL,BBLAST)              H18510
            BB = BB-BBDEL                                                 H18520
            VIDVEM = -VIDV                                                H18530
            EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMLAST)                  H18540
            EMISIV = EMISIV-EMDEL                                         H18550
            VIDVRF = -VIDV                                                H18560
            REFLCT = REFLFN(VI,DVPO,VIDVRF,RFDEL,RFLAST)                  H18570
            REFLCT = REFLCT-RFDEL                                         H18580
         ELSEIF (IEMBB.EQ.1) THEN                                         H18590
            EMISIV = EMISFN(VI,DVPO,VIDV,EMDEL,EMLAST)                    H18600
            EMISIV = EMISIV-EMDEL                                         H18610
            VIDVBD = -VIDV                                                H18620
            BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDVBD,BBDEL,BBLAST)            H18630
            BB = BB-BBDEL                                                 H18640
            VIDVRF = -VIDV                                                H18650
            REFLCT = REFLFN(VI,DVPO,VIDVRF,RFDEL,RFLAST)                  H18660
            REFLCT = REFLCT-RFDEL                                         H18670
         ELSE                                                             H18680
            REFLCT = REFLFN(VI,DVPO,VIDV,RFDEL,RFLAST)                    H18690
            REFLCT = REFLCT-RFDEL                                         H18700
            VIDVBD = -VIDV                                                H18710
            BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDVBD,BBDEL,BBLAST)            H18720
            BB = BB-BBDEL                                                 H18730
            VIDVEM = -VIDV                                                H18740
            EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMLAST)                  H18750
            EMISIV = EMISIV-EMDEL                                         H18760
         ENDIF                                                            H18770
C                                                                         H18780
         IF (VIDV.GE.9.E+4) THEN 
            NLIM2 = NLIMO+1
         ELSE
            NLIM2 = (VIDV-V1PO)/DVPO+1.001                                H18790
         ENDIF
         NLIM2 = MIN(NLIM2,NLIMO)                                         H18800
C                                                                         H18810
         DO 70 J = NLIM1, NLIM2                                           H18820
            EMISIV = EMISIV+EMDEL                                         H18830
            REFLCT = REFLCT+RFDEL                                         H18840
            BB = BB+BBDEL                                                 H18850
            NEWEM(J) = EMLAYR(J)+EMLAYB(J)*REFLCT*TRLAYR(J)+              H18860
     *                 TRLAYR(J)*BB*EMISIV                                H18870
   70    CONTINUE                                                         H18880
C                                                                         H18890
         IF (NLIM2.LT.NLIMO) GO TO 60                                     H18900
C                                                                         H18910
      ENDIF                                                               H18920
      CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,NEWEM,NEWTR,MFILE,NPTS,NPANLS)     H18930
      GO TO 20                                                            H18940
   80 CALL CPUTIM (TIME1)                                                 H18950
      TIME = TIME1-TIME                                                   H18960
      WRITE (IPR,915) TIME,TIMEM                                          H18970
C                                                                         H18980
      RETURN                                                              H18990
C                                                                         H19000
  900 FORMAT (' TIME AT THE START OF --EMINIT-- ',F10.3)                  H19010
  905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,                  H19020
     *        '0 INPUT FILE =',I5,' OUTPUT FILE =',I5)                    H19030
  910 FORMAT ('0 IPATHL AND IANT',2I5)                                    H19040
  915 FORMAT (' TIME REQUIRED FOR --EMINIT-- ',F10.3,                     H19050
     *        ' --EMIN-- ',F10.3)                                         H19060
C                                                                         H19070
      END                                                                 H19080
C
C     ----------------------------------------------------------------
C
      SUBROUTINE RADMRG (NPTS,LFILE,MFILE,JPATHL,TBND)                    H19090
C                                                                         H19100
      IMPLICIT DOUBLE PRECISION (V)                                     ! H19110
C                                                                         H19120
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H19130
C                                                                         H19140
C               LAST MODIFICATION:    8 APRIL 1991                        H19150
C                                                                         H19160
C                  IMPLEMENTATION:    R.D. WORSHAM                        H19170
C                                                                         H19180
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         H19190
C                                     R.D. WORSHAM                        H19200
C                                     J.L. MONCET                         H19210
C                                                                         H19220
C                                                                         H19230
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H19240
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H19250
C                                                                         H19260
C----------------------------------------------------------------------   H19270
C                                                                         H19280
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H19290
C                                     OFFICE OF ENERGY RESEARCH           H19300
C                                     DEPARTMENT OF ENERGY                H19310
C                                                                         H19320
C                                                                         H19330
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             H19340
C                                                                         H19350
C                                             FASCOD3                     H19360
C                                                                         H19370
C                                                                         H19380
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H19390
C                                                                         H19400
      COMMON RADN(2410),TRAN(2410),RADO(0:5000)                           H19410
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H19420
     *              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       H19430
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H19440
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H19450
     *              EXTID(10)                                             H19460
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H19470
C                                                                         H19480
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H19490
C                                                                         H19500
      COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),        H19510
     *               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,    H19520
     *               EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF     H19530
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP                   H19540
      COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO                                 H19550
      COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS      H19560
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H19570
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        H19580
     *              NLTEFL,LNFIL4,LNGTH4                                  H19590
      COMMON /XME/ TRAO(0:5000)                                           H19600
      COMMON /RMRG/ XKT,XKTA,XKTB,SECNT                                   H19610
C                                                                         H19620
      DIMENSION RADLYB(2410)                                              H19630
      DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)                             H19640
      DIMENSION A1(10),A2(10),A3(10),A4(10)                               H19650
      DIMENSION RADLYR(2),TRALYR(2),RADOI(2),TRAOI(2)                     H19660
      DIMENSION WKSAV(35)                                                 H19670
C                                                                         H19680
      CHARACTER*40 CYID                                                   H19690
C                                                                         H19700
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                   H19710
     *            (OPNLHD(1),V1PO)                                        H19720
      EQUIVALENCE (RADO(1),RADOI(1)) , (TRAO(1),TRAOI(1)),                H19730
     *            (RADN(1),RADLYR(1)) , (TRAN(1),TRALYR(1)),              H19740
     *            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),                 H19750
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 H19760
     *            (FSCDID(16),LAYR1)                                      H19770
C                                                                         H19780
      DATA NDIM / 2410 /,ND2 / 5000 /                                     H19790
C                                                                         H19800
C                                                                         H19810
C                                                                         H19820
C      ************************************************************       H19830
C      ****** THIS SUBROUTINE DOES LAYER MERGE FOR RADIANCE  ******       H19840
C      ************************************************************       H19850
C                                                                         H19860
C     IPATHL =-1 IS FOR THE LOOKING DOWN CASE FOR REFLECTED ATMOSPHERE    H19870
C     IPATHL = 1 IS FOR THE LOOKING DOWN CASE (TO DENSER LAYERS)          H19880
C     IPATHL = 2 IS FOR THE SYMMETRIC TANGENT PATH CASE                   H19890
C     IPATHL = 3 IS FOR THE LOOKING UP CASE (TO LESS DENSE LAYERS)        H19900
C                                                                         H19910
C                                                                         H19920
C      ** NOTE ON IPATHL = 2                                              H19930
C            THE TANGENT MERGE ROUTINES ARE DIVIDED INTO ANTERIOR (1ST)   H19940
C            AND POSTERIOR (2ND) LAYER CROSSINGS   THIS RECOGNITION IS    H19950
C            TRIGGERED BY THE VALUE OF "IANT"                             H19960
C                                                                         H19970
C          IF  IANT.EQ. 1  THEN POSTERIOR MERGE                           H19980
C          IF  IANT.EQ. 0  THEN NORMAL MERGE                              H19990
C          IF  IANT.EQ.-1  THEN ANTERIOR MERGE                            H20000
C                                                                         H20010
      CALL CPUTIM (TIME)                                                  H20020
      WRITE (IPR,900) TIME                                                H20030
      NPANLS = 0                                                          H20040
      TIMEM = 0.0                                                         H20050
      TIMRD = 0.0                                                         H20060
      TIMOT = 0.0                                                         H20070
C                                                                         H20080
      CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)                            H20090
      LAY1SV = LAYR1                                                      H20100
      DVL = DV                                                            H20110
      PL = PAVE                                                           H20120
      TL = TAVE                                                           H20130
      WTOTL = 0.                                                          H20140
C                                                                         H20150
      DO 10 MOL = 1, NMOL                                                 H20160
         WTOTL = WTOTL+WK(MOL)                                            H20170
         WKSAV(MOL) = WK(MOL)                                             H20180
   10 CONTINUE                                                            H20190
C                                                                         H20200
      WTOTL = WTOTL+WBROAD                                                H20210
      WN2SAV = WBROAD                                                     H20220
C                                                                         H20230
C     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)                 H20240
C                                                                         H20250
      IF (IAERSL.GT.0) WRITE (CYID,'(5A8)') (YID(I),I=3,7)                H20260
      CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)                            H20270
      IF (IAERSL.GT.0) READ (CYID,'(5A8)') (YID(I),I=3,7)                 H20280
C                                                                         H20290
      IF (JPATHL.GE.1) IPATHL = JPATHL                                    H20300
      PLAY = PAVE                                                         H20310
      TLAY = TAVE                                                         H20320
C                                                                         H20330
C     IF BOUNDARY PROPERTIES ARE SUPPLIED, AND DOWNWARD LOOKING           H20340
C     CASE; SET IPATHL TO REFLECTED ATMOSPHERE CASE                       H20350
C                                                                         H20360
      IF (IBPROP.EQ.1.AND.IPATHL.EQ.1) IPATHL = -1                        H20370
      TAVK = TAVE                                                         H20380
      DVK = DV                                                            H20390
      FACT = 1.                                                           H20400
      IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.                            H20410
C                                                                         H20420
      IF (DVL.EQ.DVK) THEN                                                H20430
         ITYPE = 0                                                        H20440
      ELSEIF (DVL.GT.DVK) THEN                                            H20450
         ITYPE = DVK/(DVL-DVK)+0.5                                        H20460
      ELSE                                                                H20470
C                                                                         H20480
C     DVL.LT.DVK                                                          H20490
C                                                                         H20500
         ITYPE = -INT(DVL/(DVK-DVL)+0.5)                                  H20510
      ENDIF                                                               H20520
      IF (ITYPE.LT.0) STOP ' RADMRG; ITYPE LT 0 '                         H20530
C                                                                         H20540
      WTOTK = 0.                                                          H20550
      LAYR1 = LAY1SV                                                      H20560
      WRITE (IPR,905) LAYR1,LAYER,KFILE,LFILE,MFILE                       H20570
      IEMIT = 1                                                           H20580
      DO 20 MOL = 1, NMOL                                                 H20590
         WTOTK = WTOTK+WK(MOL)*FACT                                       H20600
         WK(MOL) = WK(MOL)*FACT+WKSAV(MOL)                                H20610
   20 CONTINUE                                                            H20620
      WTOTK = WTOTK+WBROAD*FACT                                           H20630
      WBROAD = WBROAD*FACT+WN2SAV                                         H20640
      PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)                          H20650
      TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)                          H20660
      SECANT = 0.                                                         H20670
C                                                                         H20680
C     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES               H20690
C                                                                         H20700
      CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)                                H20710
      DVXM = DV                                                           H20720
      XKT = TAVK/RADCN2                                                   H20730
C                                                                         H20740
      WRITE (IPR,910) IPATHL,IANT                                         H20750
C                                                                         H20760
      IF (IPATHL.EQ.-1) THEN                                              H20770
         XKTA = TZU/RADCN2                                                H20780
         XKTB = TZL/RADCN2                                                H20790
      ELSEIF (IPATHL.EQ.1) THEN                                           H20800
         XKTA = TZU/RADCN2                                                H20810
         XKTB = 0.                                                        H20820
      ELSEIF (IPATHL.EQ.2) THEN                                           H20830
         XKTA = TZU/RADCN2                                                H20840
         XKTB = TZL/RADCN2                                                H20850
      ELSEIF (IPATHL.EQ.3) THEN                                           H20860
         XKTA = TZL/RADCN2                                                H20870
         XKTB = 0.                                                        H20880
      ELSE                                                                H20890
         STOP ' RADMRG; IPATHL '                                          H20900
      ENDIF                                                               H20910
C                                                                         H20920
      ATYPE = ITYPE                                                       H20930
      LL = ITYPE+1                                                        H20940
      AP = 1.0/(ATYPE+1.0)                                                H20950
C                                                                         H20960
C     A1, A2, A3 AND A4 ARE THE CONSTANTS                                 H20970
C     FOR THE LAGRANGE 4 POINT INTERPOLATION                              H20980
C                                                                         H20990
      DO 30 JPG = 1, ITYPE                                                H21000
         APG = JPG                                                        H21010
         IPL = JPG+1                                                      H21020
         P = 1.0-(AP*APG)                                                 H21030
         A1(IPL) = -P*(P-1.0)*(P-2.0)/6.0                                 H21040
         A2(IPL) = (P**2-1.0)*(P-2.0)*0.5                                 H21050
         A3(IPL) = -P*(P+1.0)*(P-2.0)*0.5                                 H21060
         A4(IPL) = P*(P**2-1.0)/6.0                                       H21070
   30 CONTINUE                                                            H21080
C                                                                         H21090
C     *** BEGINNING OF LOOP THAT DOES MERGE ***                           H21100
C                                                                         H21110
      NPE = 0                                                             H21120
      RADO(0) = 0.0                                                       H21130
      TRAO(0) = 0.0                                                       H21140
      V1PO = 0.0                                                          H21150
      V2PO = 0.0                                                          H21160
      DVPO = 0.0                                                          H21170
C                                                                         H21180
   40 CONTINUE                                                            H21190
C                                                                         H21200
      CALL CPUTIM (TIMEM1)                                                H21210
      CALL EMIN (V1P,V2P,DVP,NLIM,KFILE,RADLYR,RADLYB,TRALYR,KEOF,        H21220
     *           NPANLS)                                                  H21230
      CALL CPUTIM (TIMEM2)                                                H21240
      TIMEM = TIMEM+TIMEM2-TIMEM1                                         H21250
      IF (KEOF.LE.0) GO TO 80                                             H21260
      II = 1                                                              H21270
C                                                                         H21280
      IF (V2PO.LE.V2P+DVPO) THEN                                          H21290
   50    CALL CPUTIM (TIMEM1)                                             H21300
         CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)                         H21310
         IF (LEOF.LE.0) GO TO 60                                          H21320
         CALL BUFIN (LFILE,LEOF,RADOI(NPE+1),NLIMO)                       H21330
         CALL BUFIN (LFILE,LEOF,TRAOI(NPE+1),NLIMO)                       H21340
         CALL CPUTIM (TIMEM2)                                             H21350
         TIMRD = TIMRD+TIMEM2-TIMEM1                                      H21360
         NPE = NLIMO+NPE                                                  H21370
         IF (V2PO.LE.V2P+DVPO) GO TO 50                                   H21380
      ENDIF                                                               H21390
C                                                                         H21400
C     ZERO POINT OF FIRST PANEL                                           H21410
C                                                                         H21420
   60 IF (RADO(0).EQ.0.0.AND.TRAO(0).EQ.0.0) THEN                         H21430
         RADO(0) = RADO(1)                                                H21440
         TRAO(0) = TRAO(1)                                                H21450
      ENDIF                                                               H21460
C                                                                         H21470
C     END POINT OF LAST PANEL                                             H21480
C                                                                         H21490
      IF (V2PO+DVPO.GE.V2) THEN                                           H21500
         RADO(NPE+1) = RADO(NPE)                                          H21510
         RADO(NPE+2) = RADO(NPE)                                          H21520
         TRAO(NPE+1) = TRAO(NPE)                                          H21530
         TRAO(NPE+2) = TRAO(NPE)                                          H21540
      ENDIF                                                               H21550
C                                                                         H21560
      NPL = 1                                                             H21570
C                                                                         H21580
C     NPL IS LOCATION OF FIRST ELEMENT ON ARRAYS RADO AND TRAO            H21590
C                                                                         H21600
      CALL RADNN (RADN,TRAN,RADO,TRAO,RADLYB,NLIM,NDIM,ND2,V1P,DVP,       H21610
     *           IPATHL,A1,A2,A3,A4,LL,NPL)                               H21620
C                                                                         H21630
      CALL CPUTIM (TIMEM1)                                                H21640
C                                                                         H21650
      IF (TBND.GT.0.) CALL EMBND (V1P,V2P,DVP,NLIM,RADN,TRAN,TBND)        H21660
C                                                                         H21670
      CALL EMOUT (V1P,V2P,DVP,NLIM,RADN,TRAN,MFILE,NPTS,NPANLS)           H21680
      CALL CPUTIM (TIMEM2)                                                H21690
      TIMOT = TIMOT+TIMEM2-TIMEM1                                         H21700
C                                                                         H21710
C     NPL IS NOW LOCATION OF FIRST ELEMENT TO BE USED FOR NEXT PASS       H21720
C                                                                         H21730
      IPL = -1                                                            H21740
      DO 70 NL = NPL, NPE                                                 H21750
         IPL = IPL+1                                                      H21760
         RADO(IPL) = RADO(NL)                                             H21770
         TRAO(IPL) = TRAO(NL)                                             H21780
   70 CONTINUE                                                            H21790
C                                                                         H21800
      NPE = IPL                                                           H21810
C                                                                         H21820
      GO TO 40                                                            H21830
   80 CONTINUE                                                            H21840
C                                                                         H21850
      CALL CPUTIM (TIME1)                                                 H21860
      TIM = TIME1-TIME                                                    H21870
      WRITE (IPR,915) TIME1,TIM,TIMEM,TIMRD,TIMOT                         H21880
C                                                                         H21890
      RETURN                                                              H21900
C                                                                         H21910
C                                                                         H21920
  900 FORMAT ('0 THE TIME AT THE START OF RADMRG IS ',F12.3)              H21930
  905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,     H21940
     *        ' MERGED WITH FILE ',I5,' ONTO FILE',I5)                    H21950
  910 FORMAT ('0 IPATHL AND IANT',2I5)                                    H21960
  915 FORMAT ('0 THE TIME AT THE END OF RADMRG IS ',F12.3/F12.3,          H21970
     *        ' SECS WERE REQUIRED FOR THIS MERGE  - EMIN - ',            H21980
     *        F12.3,' - READ - ',F12.3,' - EMOUT - ',F12.3)               H21990
C                                                                         H22000
      END                                                                 H22010
C
C     ----------------------------------------------------------------
C
      SUBROUTINE RADNN (RADLYR,TRALYR,RADO,TRAO,RADLYB,NLIM,NDIM,ND2,     H22020
     *                  V1P,DVP,IPATHL,A1,A2,A3,A4,LL,NPL)                H22030
C                                                                         H22040
      IMPLICIT DOUBLE PRECISION (V)                                     ! H22050
C                                                                         H22060
C     THIS SUBROUTINE CALCULATES THE NEW RADIANCE AND TRANSMISSION        H22070
C                                                                         H22080
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H22090
C                                                                         H22100
C               LAST MODIFICATION:    8 APRIL 1991                        H22110
C                                                                         H22120
C                  IMPLEMENTATION:    R.D. WORSHAM                        H22130
C                                                                         H22140
C                       ALGORITHM:    R.D. WORSHAM                        H22150
C                                     S.A. CLOUGH                         H22160
C                                     J.L. MONCET                         H22170
C                                                                         H22180
C                                                                         H22190
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H22200
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H22210
C                                                                         H22220
C----------------------------------------------------------------------   H22230
C                                                                         H22240
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H22250
C                                     OFFICE OF ENERGY RESEARCH           H22260
C                                     DEPARTMENT OF ENERGY                H22270
C                                                                         H22280
C                                                                         H22290
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H22300
C                                                                         H22310
                                                                          H22320
      DIMENSION RADLYR(NDIM),TRALYR(NDIM),RADO(0:ND2),TRAO(0:ND2),        H22330
     *          RADLYB(NDIM),A1(*),A2(*),A3(*),A4(*)                      H22340
C                                                                         H22350
      LLM1 = LL-1                                                         H22360
      LLM1 = MAX(LLM1,1)                                                  H22370
C                                                                         H22380
C     LOOPING OVER POINTS WITH SAME WEIGHTS                               H22390
C                                                                         H22400
      DO 110 NL = 1, LL                                                   H22410
         IPL = (NPL+NL-1)-LLM1                                            H22420
         IF (NL.GT.1) IPL = IPL-1                                         H22430
C                                                                         H22440
         IF (NL.EQ.1) THEN                                                H22450
C                                                                         H22460
C     EXACT FREQUENCY - NO INTERPOLATION                                  H22470
C                                                                         H22480
            IF (IPATHL.EQ.1) THEN                                         H22490
C                                                                         H22500
               DO 10 I = NL, NLIM, LL                                     H22510
                  IPL = IPL+LLM1                                          H22520
                  RADLYR(I) = TRALYR(I)*RADO(IPL)+RADLYR(I)               H22530
                  TRALYR(I) = TRALYR(I)*TRAO(IPL)                         H22540
   10          CONTINUE                                                   H22550
C                                                                         H22560
            ELSEIF (IPATHL.EQ.2) THEN                                     H22570
C                                                                         H22580
               DO 20 I = NL, NLIM, LL                                     H22590
                  IPL = IPL+LLM1                                          H22600
                  TRTEMP = TRALYR(I)*TRAO(IPL)                            H22610
                  RADLYR(I) = TRALYR(I)*RADO(IPL)+RADLYR(I)+              H22620
     *                        RADLYB(I)*TRTEMP                            H22630
                  TRALYR(I) = TRALYR(I)*TRTEMP                            H22640
   20          CONTINUE                                                   H22650
C                                                                         H22660
            ELSEIF (IPATHL.EQ.3) THEN                                     H22670
C                                                                         H22680
               DO 30 I = NL, NLIM, LL                                     H22690
                  IPL = IPL+LLM1                                          H22700
                  RADLYR(I) = RADO(IPL)+RADLYR(I)*TRAO(IPL)               H22710
                  TRALYR(I) = TRALYR(I)*TRAO(IPL)                         H22720
   30          CONTINUE                                                   H22730
C                                                                         H22740
            ELSEIF (IPATHL.EQ.-1) THEN                                    H22750
C                                                                         H22760
               VI = V1P-DVP                                               H22770
               DVI = DVP*FLOAT(LL)                                        H22780
               VIDVRF = VI                                                H22790
               RFLAST = -1.                                               H22800
               NLIM1 = 0                                                  H22810
               NLIM2 = NL-1                                               H22820
C                                                                         H22830
   40          NLIM1 = NLIM2+1                                            H22840
C                                                                         H22850
               VI = V1P+FLOAT(NLIM1-1)*DVP                                H22860
               REFLCT = REFLFN(VI,DVI,VIDVRF,RFDEL,RFLAST)                H22870
               REFLCT = REFLCT-RFDEL                                      H22880
C                                                                         H22890
               IF (VIDVRF.GE.9.E+4) THEN 
                  NLIM2 = NLIM+1
               ELSE
                  NLIM2 = (VIDVRF+DVI-DVP-V1P)/DVP+1.001                  H22900
               ENDIF
               NLIM2 = MIN(NLIM2,NLIM)                                    H22910
C                                                                         H22920
               DO 50 I = NLIM1, NLIM2, LL                                 H22930
                  IPL = IPL+LLM1                                          H22940
                  REFLCT = REFLCT+RFDEL                                   H22950
                  TRTEMP = TRALYR(I)*TRAO(IPL)                            H22960
                  RADLYR(I) = TRALYR(I)*RADO(IPL)+RADLYR(I)+              H22970
     *                        RADLYB(I)*TRAO(IPL)*TRTEMP*REFLCT           H22980
                  TRALYR(I) = TRTEMP                                      H22990
   50          CONTINUE                                                   H23000
C                                                                         H23010
               IF (NLIM2.LT.NLIM) GO TO 40                                H23020
C                                                                         H23030
            ENDIF                                                         H23040
C                                                                         H23050
C     NOT EXACT FREQUENCY - INTERPOLATE RESULT                            H23060
C                                                                         H23070
         ELSE                                                             H23080
C                                                                         H23090
            A1N = A1(NL)                                                  H23100
            A2N = A2(NL)                                                  H23110
            A3N = A3(NL)                                                  H23120
            A4N = A4(NL)                                                  H23130
C                                                                         H23140
            IF (IPATHL.EQ.1) THEN                                         H23150
               DO 60 I = NL, NLIM, LL                                     H23160
                  IPL = IPL+LLM1                                          H23170
C                                                                         H23180
C     INTERPOLATE THE OLD RADIANCE                                        H23190
C                                                                         H23200
                  RADLYR(I) = TRALYR(I)*(A1N*RADO(IPL-1)+A2N*RADO(IPL)+   H23210
     *                        A3N*RADO(IPL+1)+A4N*RADO(IPL+2))+RADLYR(I)  H23220
C                                                                         H23230
C     INTERPOLATE THE OLD TRANSMISSION                                    H23240
C                                                                         H23250
                  TRALYR(I) = TRALYR(I)*(A1N*TRAO(IPL-1)+A2N*TRAO(IPL)+   H23260
     *                        A3N*TRAO(IPL+1)+A4N*TRAO(IPL+2))            H23270
   60          CONTINUE                                                   H23280
C                                                                         H23290
            ELSEIF (IPATHL.EQ.2) THEN                                     H23300
C                                                                         H23310
               DO 70 I = NL, NLIM, LL                                     H23320
                  IPL = IPL+LLM1                                          H23330
C                                                                         H23340
C     INTERPOLATE THE OLD TRANSMISSION                                    H23350
C                                                                         H23360
                  TRTEMP = TRALYR(I)*(A1N*TRAO(IPL-1)+A2N*TRAO(IPL)+      H23370
     *                     A3N*TRAO(IPL+1)+A4N*TRAO(IPL+2))               H23380
C                                                                         H23390
C     INTERPOLATE THE OLD RADIANCE                                        H23400
C                                                                         H23410
                  RADLYR(I) = TRALYR(I)*(A1N*RADO(IPL-1)+A2N*RADO(IPL)+   H23420
     *                        A3N*RADO(IPL+1)+A4N*RADO(IPL+2))+           H23430
     *                        RADLYR(I)+RADLYB(I)*TRTEMP                  H23440
                  TRALYR(I) = TRALYR(I)*TRTEMP                            H23450
   70          CONTINUE                                                   H23460
C                                                                         H23470
            ELSEIF (IPATHL.EQ.3) THEN                                     H23480
C                                                                         H23490
               DO 80 I = NL, NLIM, LL                                     H23500
                  IPL = IPL+LLM1                                          H23510
C                                                                         H23520
C     INTERPOLATE THE OLD TRANSMISSION                                    H23530
C                                                                         H23540
                  TRAOI = A1N*TRAO(IPL-1)+A2N*TRAO(IPL)+                  H23550
     *                    A3N*TRAO(IPL+1)+A4N*TRAO(IPL+2)                 H23560
C                                                                         H23570
C     INTERPOLATE THE OLD RADIANCE                                        H23580
C                                                                         H23590
                  RADLYR(I) = A1N*RADO(IPL-1)+A2N*RADO(IPL)+              H23600
     *                        A3N*RADO(IPL+1)+A4N*RADO(IPL+2)+            H23610
     *                        RADLYR(I)*TRAOI                             H23620
                  TRALYR(I) = TRALYR(I)*TRAOI                             H23630
   80          CONTINUE                                                   H23640
C                                                                         H23650
            ELSEIF (IPATHL.EQ.-1) THEN                                    H23660
C                                                                         H23670
               VI = V1P-DVP                                               H23680
               DVI = DVP*FLOAT(LL)                                        H23690
               VIDVRF = VI                                                H23700
               RFLAST = -1                                                H23710
               NLIM1 = 0                                                  H23720
               NLIM2 = NL-1                                               H23730
C                                                                         H23740
   90          NLIM1 = NLIM2+1                                            H23750
C                                                                         H23760
               VI = V1P+FLOAT(NLIM1-1)*DVP                                H23770
               REFLCT = REFLFN(VI,DVI,VIDVRF,RFDEL,RFLAST)                H23780
               REFLCT = REFLCT-RFDEL                                      H23790
C                                                                         H23800
               IF (VIDVRF.GE.9.E+4) THEN 
                  NLIM2 = NLIM+1
               ELSE
                  NLIM2 = (VIDVRF+DVI-DVP-V1P)/DVP+1.001                  H23810
               ENDIF
               NLIM2 = MIN(NLIM2,NLIM)                                    H23820
C                                                                         H23830
               DO 100 I = NLIM1, NLIM2, LL                                H23840
                  IPL = IPL+LLM1                                          H23850
                  REFLCT = REFLCT+RFDEL                                   H23860
C                                                                         H23870
C     INTERPOLATE THE OLD TRANSMISSION                                    H23880
C                                                                         H23890
                  TRAOI = A1N*TRAO(IPL-1)+A2N*TRAO(IPL)+                  H23900
     *                    A3N*TRAO(IPL+1)+A4N*TRAO(IPL+2)                 H23910
C                                                                         H23920
                  TRTEMP = TRALYR(I)*TRAOI                                H23930
C                                                                         H23940
C     INTERPOLATE THE OLD RADIANCE                                        H23950
C                                                                         H23960
                  RADLYR(I) = TRALYR(I)*(A1N*RADO(IPL-1)+A2N*RADO(IPL)+   H23970
     *                        A3N*RADO(IPL+1)+A4N*RADO(IPL+2))+           H23980
     *                        RADLYR(I)+RADLYB(I)*TRAOI*TRTEMP*REFLCT     H23990
                  TRALYR(I) = TRTEMP                                      H24000
  100          CONTINUE                                                   H24010
C                                                                         H24020
               IF (NLIM2.LT.NLIM) GO TO 90                                H24030
C                                                                         H24040
            ENDIF                                                         H24050
C                                                                         H24060
         ENDIF                                                            H24070
C                                                                         H24080
  110 CONTINUE                                                            H24090
C                                                                         H24100
      NPL = IPL                                                           H24110
C                                                                         H24120
      RETURN                                                              H24130
C                                                                         H24140
      END                                                                 H24150
C
C     ----------------------------------------------------------------
C
      SUBROUTINE RADINT (NPTS,LFILE,MFILE,JPATHL,TBND)                    H24160
C                                                                         H24170
      IMPLICIT DOUBLE PRECISION (V)                                     ! H24180
C                                                                         H24190
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H24200
C                                                                         H24210
C               LAST MODIFICATION:    5 APRIL 1991                        H24220
C                                                                         H24230
C                  IMPLEMENTATION:    R.D. WORSHAM                        H24240
C                                                                         H24250
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         H24260
C                                     R.D. WORSHAM                        H24270
C                                     J.L. MONCET                         H24280
C                                                                         H24290
C                                                                         H24300
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H24310
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H24320
C                                                                         H24330
C----------------------------------------------------------------------   H24340
C                                                                         H24350
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H24360
C                                     OFFICE OF ENERGY RESEARCH           H24370
C                                     DEPARTMENT OF ENERGY                H24380
C                                                                         H24390
C                                                                         H24400
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             H24410
C                                                                         H24420
C                                             FASCOD3                     H24430
C                                                                         H24440
C                                                                         H24450
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H24460
C                                                                         H24470
      COMMON RADN(2410),TRAN(2410),RADLYR(-1:4818)                        H24480
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H24490
     *              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       H24500
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H24510
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H24520
     *              EXTID(10)                                             H24530
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H24540
C                                                                         H24550
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H24560
C                                                                         H24570
      COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),        H24580
     *               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,    H24590
     *               EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF     H24600
      COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO                                 H24610
      COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS      H24620
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H24630
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        H24640
     *              NLTEFL,LNFIL4,LNGTH4                                  H24650
      COMMON /XMI/ TRALYR(-1:4818)                                        H24660
      COMMON /RMRG/ XKT,XKTA,XKTB,SECNT                                   H24670
C                                                                         H24680
      DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)                             H24690
      DIMENSION A1(0:100),A2(0:100),A3(0:100),A4(0:100)                   H24700
      DIMENSION RADO(2),TRAO(2)                                           H24710
      DIMENSION WKSAV(35)                                                 H24720
C                                                                         H24730
      DIMENSION RADLYB(-1:4818)                                           H24740
C                                                                         H24750
      CHARACTER*40 CYID                                                   H24760
C                                                                         H24770
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                   H24780
     *            (OPNLHD(1),V1PO)                                        H24790
      EQUIVALENCE (RADN(1),RADO(1)) , (TRAN(1),TRAO(1)),                  H24800
     *            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),                 H24810
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 H24820
     *            (FSCDID(16),LAYR1)                                      H24830
C                                                                         H24840
C     ************************************************************        H24850
C     ****** THIS SUBROUTINE DOES LAYER MERGE FOR EMISSION  ******        H24860
C     ****** USING FOUR POINT GENERAL INTERPOLATION         ******        H24870
C     ************************************************************        H24880
C                                                                         H24890
      CALL CPUTIM (TIME)                                                  H24900
      WRITE (IPR,900) TIME                                                H24910
      NPANLS = 0                                                          H24920
      TIMEM = 0.0                                                         H24930
      TIMRD = 0.0                                                         H24940
      TIMTB = 0.0                                                         H24950
      TIMOT = 0.0                                                         H24960
C                                                                         H24970
      CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)                            H24980
      DVL = DV                                                            H24990
      LAY1SV = LAYR1                                                      H25000
      PL = PAVE                                                           H25010
      TL = TAVE                                                           H25020
      WTOTL = 0.                                                          H25030
      DO 10 MOL = 1, NMOL                                                 H25040
         WTOTL = WTOTL+WK(MOL)                                            H25050
         WKSAV(MOL) = WK(MOL)                                             H25060
   10 CONTINUE                                                            H25070
      WTOTL = WTOTL+WBROAD                                                H25080
      WN2SAV = WBROAD                                                     H25090
C                                                                         H25100
C     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)                 H25110
C                                                                         H25120
      IF (IAERSL.GT.0) WRITE (CYID,'(5A8)') (YID(I),I=3,7)                H25130
      CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)                            H25140
      IF (IAERSL.GT.0) READ (CYID,'(5A8)') (YID(I),I=3,7)                 H25150
      IF (JPATHL.GE.1) IPATHL = JPATHL                                    H25160
      PLAY = PAVE                                                         H25170
      TLAY = TAVE                                                         H25180
      XKT = TAVE/RADCN2                                                   H25190
      XKTA = TZU/RADCN2                                                   H25200
      XKTB = 0.                                                           H25210
      DVK = DV                                                            H25220
      LAYR1 = LAY1SV                                                      H25230
      FACT = 1.                                                           H25240
      IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.                            H25250
      ATYPE = 9.999E09                                                    H25260
      IF (DVK.EQ.DVL) ATYPE = 0.                                          H25270
      IF (DVL.GT.DVK) ATYPE = DVK/(DVL-DVK)+0.5                           H25280
      IF (DVL.LT.DVK) ATYPE = -DVL/(DVK-DVL)-0.5                          H25290
C                                                                         H25300
C     IF (ATYPE .GT. 0) STOP  ' RADINT; ATYPE GT 0 '                      H25310
C                                                                         H25320
      WTOTK = 0.                                                          H25330
      WRITE (IPR,905) LAYR1,LAYER,KFILE,LFILE,MFILE,ATYPE                 H25340
      IEMIT = 1                                                           H25350
      DO 20 MOL = 1, NMOL                                                 H25360
         WTOTK = WTOTK+WK(MOL)*FACT                                       H25370
         WK(MOL) = WK(MOL)*FACT+WKSAV(MOL)                                H25380
   20 CONTINUE                                                            H25390
      WTOTK = WTOTK+WBROAD*FACT                                           H25400
      WBROAD = WBROAD*FACT+WN2SAV                                         H25410
      PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)                          H25420
      TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)                          H25430
      SECANT = 0.                                                         H25440
      DV = DVL                                                            H25450
C                                                                         H25460
C     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES               H25470
C                                                                         H25480
      CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)                                H25490
      DVXM = DV                                                           H25500
C                                                                         H25510
      IF (ATYPE.EQ.0.) THEN                                               H25520
C                                                                         H25530
C     1/1 RATIO ONLY                                                      H25540
C                                                                         H25550
   30    CONTINUE                                                         H25560
         CALL CPUTIM (TIMEM1)                                             H25570
         CALL EMIN (V1P,V2P,DVP,NLIM,KFILE,RADLYR(1),RADLYB(1),           H25580
     *              TRALYR(1),KEOF,NPANLS)                                H25590
         CALL CPUTIM (TIMEM2)                                             H25600
         TIMEM = TIMEM+TIMEM2-TIMEM1                                      H25610
         IF (KEOF.LE.0) GO TO 110                                         H25620
         CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)                         H25630
         CALL BUFIN (LFILE,LEOF,RADO(1),NLIMO)                            H25640
         CALL BUFIN (LFILE,LEOF,TRAO(1),NLIMO)                            H25650
         CALL CPUTIM (TIMEM3)                                             H25660
         TIMRD = TIMRD+TIMEM3-TIMEM2                                      H25670
         DO 40 I = 1, NLIM                                                H25680
            RADN(I) = RADO(I)+RADLYR(I)*TRAO(I)                           H25690
            TRAN(I) = TRALYR(I)*TRAO(I)                                   H25700
   40    CONTINUE                                                         H25710
         CALL CPUTIM (TIMEM1)                                             H25720
         IF (TBND.GT.0.)                                                  H25730
     *       CALL EMBND (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,TBND)             H25740
         CALL CPUTIM (TIMEM2)                                             H25750
         TIMTB = TIMTB+TIMEM2-TIMEM1                                      H25760
         CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,MFILE,NPTS,NPANLS)    H25770
         CALL CPUTIM (TIMEM3)                                             H25780
         TIMOT = TIMOT+TIMEM3-TIMEM2                                      H25790
         GO TO 30                                                         H25800
C                                                                         H25810
      ENDIF                                                               H25820
C                                                                         H25830
C     ALL RATIOS EXCEPT 1/1                                               H25840
C                                                                         H25850
      DO 50 JP = 0, 100                                                   H25860
         APG = JP                                                         H25870
         P = 0.01*APG                                                     H25880
C                                                                         H25890
C     THE FOLLOW ARE THE CONSTANTS FOR THE LAGRANGE 4 POINT               H25900
C     INTERPOLATION                                                       H25910
C                                                                         H25920
         A1(JP) = -P*(P-1.0)*(P-2.0)/6.0                                  H25930
         A2(JP) = (P**2-1.0)*(P-2.0)*0.5                                  H25940
         A3(JP) = -P*(P+1.0)*(P-2.0)*0.5                                  H25950
         A4(JP) = P*(P**2-1.0)/6.0                                        H25960
   50 CONTINUE                                                            H25970
C                                                                         H25980
C     *** BEGINNING OF LOOP THAT DOES MERGE  ***                          H25990
C                                                                         H26000
      NPE = 0                                                             H26010
      RADLYR(0) = 0.0                                                     H26020
      TRALYR(0) = 0.0                                                     H26030
      V1P = 0.0                                                           H26040
      V2P = 0.0                                                           H26050
      DVP = 0.0                                                           H26060
      V1PO = 0.0                                                          H26070
      V2PO = 0.0                                                          H26080
      DVPO = 0.0                                                          H26090
      KEOF = 1                                                            H26095
C                                                                         H26100
   60 CONTINUE                                                            H26110
C                                                                         H26120
      CALL CPUTIM (TIMEM1)                                                H26130
      CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)                            H26140
      IF (LEOF.LE.0) GO TO 110                                            H26150
      CALL BUFIN (LFILE,LEOF,RADO(1),NLIMO)                               H26160
      CALL BUFIN (LFILE,LEOF,TRAO(1),NLIMO)                               H26170
      CALL CPUTIM (TIMEM2)                                                H26180
      TIMRD = TIMRD+TIMEM2-TIMEM1                                         H26190
      II = 1                                                              H26200
C                                                                         H26210
      IF (V2P.LE.V2PO+DVP .AND. KEOF.GT.0) THEN                           H26220
   70    CALL CPUTIM (TIMEM2)                                             H26230
         CALL EMIN (V1P,V2P,DVP,NLIM,KFILE,RADLYR(NPE+1),RADLYB(NPE+1),   H26240
     *              TRALYR(NPE+1),KEOF,NPANLS)                            H26250
         CALL CPUTIM (TIMEM3)                                             H26260
         TIMEM = TIMEM+TIMEM3-TIMEM2                                      H26270
         IF (KEOF.LE.0) GO TO 80                                          H26280
         V1P = V1P-FLOAT(NPE)*DVP                                         H26290
         NPE = NLIM+NPE                                                   H26300
         IF (V2P.LE.V2PO+DVP) GO TO 70                                    H26310
      ENDIF                                                               H26320
C                                                                         H26330
C     ZERO POINT OF FIRST PANEL                                           H26340
C                                                                         H26350
   80 IF (RADLYR(0).EQ.0.0.AND.TRALYR(0).EQ.0.0) THEN                     H26360
         TRALYR(-1) = TRALYR(1)                                           H26370
         RADLYR(-1) = RADLYR(1)                                           H26380
         RADLYB(-1) = RADLYB(1)                                           H26390
         TRALYR(0) = TRALYR(1)                                            H26400
         RADLYR(0) = RADLYR(1)                                            H26410
         RADLYB(0) = RADLYB(1)                                            H26420
      ENDIF                                                               H26430
C                                                                         H26440
C     END POINT OF LAST PANEL                                             H26450
C                                                                         H26460
      IF (V2P+DVP.GE.V2) THEN                                             H26470
         TRALYR(NPE+1) = TRALYR(NPE)                                      H26480
         RADLYR(NPE+1) = RADLYR(NPE)                                      H26490
         RADLYB(NPE+1) = RADLYB(NPE)                                      H26500
         TRALYR(NPE+2) = TRALYR(NPE)                                      H26510
         RADLYR(NPE+2) = RADLYR(NPE)                                      H26520
         RADLYB(NPE+2) = RADLYB(NPE)                                      H26530
      ENDIF                                                               H26540
C                                                                         H26550
C     NPL IS LOCATION OF FIRST ELEMENT ON ARRAYS RADO AND TRAO            H26560
C                                                                         H26570
      NPL = 1                                                             H26580
C                                                                         H26590
      RATDV = DVL/DVK                                                     H26600
C                                                                         H26610
C     FJJ IS OFFSET BY 2. FOR ROUNDING PURPOSES                           H26620
C                                                                         H26630
      FJ1DIF = (V1PO-V1P)/DVP+1.+2.                                       H26640
C                                                                         H26650
C     ***** BEGINNING OF LOOP THAT DOES MERGE  *****                      H26660
C                                                                         H26670
      DO 90 II = 1, NLIMO                                                 H26680
C                                                                         H26690
         FJJ = FJ1DIF+RATDV*FLOAT(II-1)                                   H26700
         JJ = IFIX(FJJ)-2                                                 H26710
C                                                                         H26720
         JP = (FJJ-FLOAT(JJ))*100.-199.5                                  H26730
C                                                                         H26740
C     INTERPOLATE THE OLD EMISSION                                        H26750
C                                                                         H26760
         RADN(II) = RADO(II)+(A1(JP)*RADLYR(JJ-1)+A2(JP)*RADLYR(JJ)+      H26770
     *              A3(JP)*RADLYR(JJ+1)+A4(JP)*RADLYR(JJ+2))*TRAO(II)     H26780
C                                                                         H26790
C     INTERPOLATE THE OLD TRANSMISSION                                    H26800
C                                                                         H26810
         TRAN(II) = (A1(JP)*TRALYR(JJ-1)+A2(JP)*TRALYR(JJ)+               H26820
     *              A3(JP)*TRALYR(JJ+1)+A4(JP)*TRALYR(JJ+2))*TRAO(II)     H26830
C                                                                         H26840
   90 CONTINUE                                                            H26850
C                                                                         H26860
      NPL = JJ-1                                                          H26870
C                                                                         H26880
      CALL CPUTIM (TIMEM1)                                                H26890
      IF (TBND.GT.0.) CALL EMBND (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,TBND)    H26900
      CALL CPUTIM (TIMEM2)                                                H26910
      TIMTB = TIMTB+TIMEM2-TIMEM1                                         H26920
      CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,MFILE,NPTS,NPANLS)       H26930
      CALL CPUTIM (TIMEM3)                                                H26940
      TIMOT = TIMOT+TIMEM3-TIMEM2                                         H26950
C                                                                         H26960
C     NPL IS NOW LOCATION OF FIRST ELEMENT TO BE USED FOR NEXT PASS       H26970
C                                                                         H26980
      IPL = -2                                                            H26990
      DO 100 NL = NPL, NPE                                                H27000
         IPL = IPL+1                                                      H27010
         TRALYR(IPL) = TRALYR(NL)                                         H27020
         RADLYR(IPL) = RADLYR(NL)                                         H27030
         RADLYB(IPL) = RADLYB(NL)                                         H27040
  100 CONTINUE                                                            H27050
C                                                                         H27060
      V1P = V1P+FLOAT(NPL+1)*DVP                                          H27070
      NPE = IPL                                                           H27080
C                                                                         H27090
      GO TO 60                                                            H27100
  110 CONTINUE                                                            H27110
C                                                                         H27120
      CALL CPUTIM (TIME1)                                                 H27130
      TIM = TIME1-TIME                                                    H27140
      WRITE (IPR,910) TIME1,TIM,TIMEM,TIMRD,TIMTB,TIMOT                   H27150
C                                                                         H27160
      RETURN                                                              H27170
C                                                                         H27180
  900 FORMAT ('0 THE TIME AT THE START OF RADINT IS ',F12.3)              H27190
  905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,     H27200
     *        ' MERGED WITH FILE ',I5,' ONTO FILE',I5,'  WITH XTYPE=',    H27210
     *        G15.5)                                                      H27220
  910 FORMAT ('0 THE TIME AT THE END OF RADINT IS ',F12.3/F12.3,          H27230
     *        ' SECS WERE REQUIRED FOR THIS MERGE  - EMIN - ',F12.3,      H27240
     *        ' - READ - ',F12.3,' - EMBND - ',F12.3,' - EMOUT - ',       H27250
     *        F12.3)                                                      H27260
C                                                                         H27270
      END                                                                 H27280
C
C     ----------------------------------------------------------------
C
      SUBROUTINE EMBND (V1PO,V2PO,DVPO,NLIMO,NEWEM,NEWTR,TBND)            H27290
C                                                                         H27300
      IMPLICIT DOUBLE PRECISION (V)                                     ! H27310
C                                                                         H27320
C                                                                         H27330
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H27340
C                                                                         H27350
C               LAST MODIFICATION:    9 APRIL 1991                        H27360
C                                                                         H27370
C                  IMPLEMENTATION:    R.D. WORSHAM                        H27380
C                                                                         H27390
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         H27400
C                                     R.D. WORSHAM                        H27410
C                                     J.L. MONCET                         H27420
C                                                                         H27430
C                                                                         H27440
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H27450
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H27460
C                                                                         H27470
C----------------------------------------------------------------------   H27480
C                                                                         H27490
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H27500
C                                     OFFICE OF ENERGY RESEARCH           H27510
C                                     DEPARTMENT OF ENERGY                H27520
C                                                                         H27530
C                                                                         H27540
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             H27550
C                                                                         H27560
C                                             FASCOD3                     H27570
C                                                                         H27580
C                                                                         H27590
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H27600
C                                                                         H27610
      DIMENSION NEWEM(*),NEWTR(*)                                         H27620
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H27630
C                                                                         H27640
      REAL NEWEM,NEWTR                                                    H27650
C                                                                         H27660
      XKTBND = TBND/RADCN2                                                H27670
      VI = V1PO-DVPO                                                      H27680
      VIDVBD = VI                                                         H27690
      VIDVEM = VI                                                         H27700
      BBLAST = -1.                                                        H27710
      EMLAST = -1.                                                        H27720
      NLIM1 = 0                                                           H27730
      NLIM2 = 0                                                           H27740
      EMDUM = 0.                                                          H27750
      BBDUM = 0.                                                          H27760
      EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMDUM)                         H27770
      BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDVBD,BBDEL,BBDUM)                   H27780
      IEMBB = 0                                                           H27790
      IF (VIDVBD.GT.VIDVEM) IEMBB = 1                                     H27800
C                                                                         H27810
   10 NLIM1 = NLIM2+1                                                     H27820
C                                                                         H27830
      VI = V1PO+FLOAT(NLIM1-1)*DVPO                                       H27840
      IF (IEMBB.EQ.0) THEN                                                H27850
         BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDV,BBDEL,BBLAST)                 H27860
         BB = BB-BBDEL                                                    H27870
         VIDVEM = -VIDV                                                   H27880
         EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMLAST)                     H27890
         EMISIV = EMISIV-EMDEL                                            H27900
      ELSE                                                                H27910
         EMISIV = EMISFN(VI,DVPO,VIDV,EMDEL,EMLAST)                       H27920
         EMISIV = EMISIV-EMDEL                                            H27930
         VIDVBD = -VIDV                                                   H27940
         BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDVBD,BBDEL,BBLAST)               H27950
         BB = BB-BBDEL                                                    H27960
      ENDIF                                                               H27970
C                                                                         H27980
      IF (VIDV.GE.9.E+4) THEN 
         NLIM2 = NLIMO+1
      ELSE
         NLIM2 = (VIDV-V1PO)/DVPO+1.001                                   H27990
      ENDIF
      NLIM2 = MIN(NLIM2,NLIMO)                                            H28000
C                                                                         H28010
      DO 20 J = NLIM1, NLIM2                                              H28020
         EMISIV = EMISIV+EMDEL                                            H28030
         BB = BB+BBDEL                                                    H28040
         NEWEM(J) = NEWEM(J)+NEWTR(J)*EMISIV*BB                           H28050
   20 CONTINUE                                                            H28060
C                                                                         H28070
      IF (NLIM2.LT.NLIMO) GO TO 10                                        H28080
C                                                                         H28090
      RETURN                                                              H28100
C                                                                         H28110
      END                                                                 H28120
C
C     ----------------------------------------------------------------
C
      SUBROUTINE EMOUT (V1P,V2P,DVP,NLIM,NEWEM,NEWTR,MFILE,NPTS,NPANLS)   H28130
C                                                                         H28140
      IMPLICIT DOUBLE PRECISION (V)                                     ! H28150
C                                                                         H28160
C     SUBROUTINE EMOUT OUTPUTS MERGED EMISSION AND TRANSMITTANCE RESULT   H28170
C     TO MFILE                                                            H28180
C                                                                         H28190
      COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF                            H28200
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H28210
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        H28220
     *              NLTEFL,LNFIL4,LNGTH4                                  H28230
      DIMENSION PNLHDR(2)                                                 H28240
      DIMENSION NEWEM(*),NEWTR(*)                                         H28250
C                                                                         H28260
      EQUIVALENCE (PNLHDR(1),V1PBF)                                       H28270
C                                                                         H28280
      REAL NEWEM,NEWTR                                                    H28290
C                                                                         H28300
      NPANLS = NPANLS+1                                                   H28310
      V1PBF = V1P                                                         H28320
      V2PBF = V2P                                                         H28330
      DVPBF = DVP                                                         H28340
      NLIMBF = NLIM                                                       H28350
C                                                                         H28360
      CALL BUFOUT (MFILE,PNLHDR(1),NPHDRF)                                H28370
      CALL BUFOUT (MFILE,NEWEM(1),NLIMBF)                                 H28380
      CALL BUFOUT (MFILE,NEWTR(1),NLIMBF)                                 H28390
C                                                                         H28400
      IF (NPTS.GT.0) THEN                                                 H28410
         IF (NPANLS.EQ.1) WRITE (IPR,900)                                 H28420
         WRITE (IPR,905)                                                  H28430
         NNPTS = NPTS                                                     H28440
         IF (NPTS.GT.(NLIMBF/2)+1) NNPTS = (NLIMBF/2)+1                   H28450
         JEND = NLIMBF-NNPTS+1                                            H28460
         DO 10 J = 1, NNPTS                                               H28470
            VJ = V1PBF+FLOAT(J-1)*DVPBF                                   H28480
            K = J+JEND-1                                                  H28490
            VK = V1PBF+FLOAT(K-1)*DVPBF                                   H28500
            WRITE (IPR,910) J,VJ,NEWEM(J),NEWTR(J),                       H28510
     *                      K,VK,NEWEM(K),NEWTR(K)                        H28520
   10    CONTINUE                                                         H28530
      ENDIF                                                               H28540
C                                                                         H28550
      RETURN                                                              H28560
C                                                                         H28570
  900 FORMAT ('0 ','LOCATION  WAVENUMBER',2X,'RADIANCE',7X,               H28580
     *        'TRANSMITTANCE',22X,'LOCATION   WAVENUMBER',2X,             H28590
     *        'RADIANCE',7X,'TRANSMITTANCE')                              H28600
  905 FORMAT (' ')                                                        H28610
  910 FORMAT (I8,2X,F12.6,1P2E15.7,20X,I8,2X,0PF12.6,1P2E15.7)            H28620
C                                                                         H28630
      END                                                                 H28640
C
C     ----------------------------------------------------------------
C
      SUBROUTINE SOLINT(IFILE,LFILE,NPTS,INFLAG,IOTFLG)
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
      IMPLICIT DOUBLE PRECISION (V)
C
C     ------------------------------------------------------------
C     SUBROUTINE SOLINT interpolates solar radiances from the binary
C     file SOLAR.RAD.  The following are input and output options:
C
C       INFLAG = 0   => input transmittance from TAPE12 (default).
C              = 1   => input optical depths from TAPE10 and
C                       convert to transmittance.
C
C       IOTFLG = 0   => attenuate w/transmittance & output (default).
C              = 1   => attenuate and add to radiance from TAPE12
C                       (requires INFLAG = 1).
C
C     Output radiance goes to TAPE11.
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
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID
C
      COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,
     *               EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
      COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
      COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4
C
      DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)
      DIMENSION A1(0:100),A2(0:100),A3(0:100),A4(0:100)
      DIMENSION TRAO(2),TRAN(2410)
      DIMENSION RADO(2),RADN(2410)
      DIMENSION OPTO(2),OPTN(2410)
C
      DIMENSION SOLAR(-1:4818)
      DIMENSION SOLRAD(2410)
C
      CHARACTER*40 CYID
C
      EQUIVALENCE (XFILHD(1),XID(1)),(PNLHDR(1),V1P),
     *            (OPNLHD(1),V1PO)
      EQUIVALENCE (TRAN(1),TRAO(1)),(RADN(1),RADO(1)),
     *            (OPTN(1),OPTO(1)),
     *            (FSCDID(4),IAERSL),(FSCDID(5),IEMIT),
     *            (FSCDID(7),IPLOT),(FSCDID(8),IPATHL),
     *            (FSCDID(16),LAYR1)
C
C     ************************************************************
C     ****** THIS PROGRAM DOES MERGE FOR SOLAR RADIANCE AND ******
C     ****** TRANMITTANCE USING FOUR POINT INTERPOLATION    ******
C     ************************************************************
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
C
      WRITE (IPR,905) ISOLFL,IFILE,LFILE,ATYPE,INFLAG,IOTFLG
      IF (INFLAG.EQ.0) THEN
         WRITE(IPR,920) IFILE
      ELSE
         WRITE(IPR,925) IFILE
      ENDIF
C
C     Test for INFLAG=0, that atmospheric radiance is to be read
C
      IF (IOTFLG.EQ.0) THEN
         WRITE(IPR,930) LFILE
      ELSE
         IF (INFLAG.EQ.1) STOP 'ERROR: INFLAG=1, IOTFLG=1'
         WRITE(IPR,935) LFILE
      ENDIF
      IEMIT = 1
      SECANT = 0.
      DV = DVL
C
C     Output file header
C
      CALL BUFOUT (LFILE,XFILHD(1),NFHDRF)
C
      IF (ATYPE.EQ.0.) THEN
C
C        1/1 ratio only
C
   30    CONTINUE
         CALL CPUTIM (TIMSL1)
         CALL SOLIN (V1P,V2P,DVP,NLIM,ISOLFL,SOLAR(1),LSEOF,NPANLS)
         CALL CPUTIM (TIMSL2)
         TIMRD = TIMRD+TIMSL2-TIMSL1
         IF (LSEOF.LE.0) GO TO 110
C
C        Read file header from transmittance/optical depth file
C
         CALL BUFIN (IFILE,LEOF,OPNLHD(1),NPHDRF)
C
C        If INFLAG = 0, then read radiance and tranmittance
C        If INFLAG = 1, then read optical depth
C
         IF (INFLAG.EQ.0) THEN
            CALL BUFIN (IFILE,LEOF,RADO(1),NLIMO)
            CALL BUFIN (IFILE,LEOF,TRAO(1),NLIMO)
         ELSE
            CALL BUFIN (IFILE,LEOF,OPTO(1),NLIMO)
            DO 35 I = 1,NLIMO
               TRAN(I) = EXP(-OPTN(I))
 35         CONTINUE
         ENDIF
         CALL CPUTIM (TIMSL3)
         TIMRD = TIMRD+TIMSL3-TIMSL2
C
C        If IOTFLG = 0, then calculate attenuated solar radiance
C        If IOTFLG = 1, then calculate attenuated solar radiance
C                       plus atmospheric radiance
C
         IF (IOTFLG.EQ.0) THEN
            DO 40 I = 1, NLIM
               SOLRAD(I) = SOLAR(I)*TRAN(I)
 40         CONTINUE
         ELSE
            DO 41 I = 1, NLIM
               SOLRAD(I) = SOLAR(I)*TRAN(I)+RADN(I)
 41         CONTINUE
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
C     All ratios except 1/1
C
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
   50 CONTINUE
C
C     *** Beginning of loop that does merge  ***
C
      NPE = 0
      SOLAR(0) = 0.0
      V1P = 0.0
      V2P = 0.0
      DVP = 0.0
      V1PO = 0.0
      V2PO = 0.0
      DVPO = 0.0
      LSEOF = 1
C
      ip = 0
   60 CONTINUE
      ip = ip+1
C
C     Read file header from transmittance/optical depth file
C
      CALL CPUTIM(TIMSL1)
      CALL BUFIN(IFILE,LEOF,OPNLHD(1),NPHDRF)
      CALL CPUTIM(TIMSL2)
      TIMRD = TIMRD+TIMSL2-TIMSL1
      IF (LEOF.LE.0) GO TO 110
C
C     If INFLAG = 0, then read radiance and tranmittance
C     If INFLAG = 1, then read optical depth
C
      IF (INFLAG.EQ.0) THEN
         CALL BUFIN (IFILE,LEOF,RADO(1),NLIMO)
         CALL BUFIN (IFILE,LEOF,TRAO(1),NLIMO)
      ELSE
         CALL BUFIN (IFILE,LEOF,OPTO(1),NLIMO)
         DO 65 I = 1,NLIMO
            TRAN(I) = EXP(-OPTN(I))
 65      CONTINUE
      ENDIF
      CALL CPUTIM(TIMSL3)
      TIMRD = TIMRD+TIMSL3-TIMSL2
      II = 1
C
C     Buffer in panels from solar radiance file
C
      IF (V2P.LE.V2PO+DVP .AND.LSEOF.GT.0) THEN
   70    CALL CPUTIM(TIMSL2)
         CALL SOLIN(V1P,V2P,DVP,NLIM,ISOLFL,SOLAR(NPE+1),LSEOF,NPANLS)
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
C     NPL is the location of first element on arrays RADO and TRAO
C
      NPL = 1
C
      RATDV = DVL/DVK
C
C     FJJ is offset by 2. (for rounding purposes)
C
      FJ1DIF = (V1PO-V1P)/DVP+1.+2.
C
C     ***** Beginning of loop that does merge  *****
C
C     If IOTFLG = 0, then calculate attenuated solar radiance
C     If IOTFLG = 1, then calculate attenuated solar radiance
C                    plus atmospheric radiance
C
      IF (IOTFLG.EQ.0) THEN
         DO 90 II = 1, NLIMO
            FJJ = FJ1DIF+RATDV*FLOAT(II-1)
            JJ = IFIX(FJJ)-2
            JP = (FJJ-FLOAT(JJ))*100.-199.5
            SOLRAD(II) = (A1(JP)*SOLAR(JJ-1)+A2(JP)*SOLAR(JJ)+
     *           A3(JP)*SOLAR(JJ+1)+A4(JP)*SOLAR(JJ+2))*TRAN(II)
c
 90      CONTINUE
      ELSE
         DO 91 II = 1, NLIMO
            FJJ = FJ1DIF+RATDV*FLOAT(II-1)
            JJ = IFIX(FJJ)-2
            JP = (FJJ-FLOAT(JJ))*100.-199.5
c            write(*,*) 'TAPE11:'
c            write(*,*) '   ',v1po+dvpo*(ii),rado(ii),trao(ii)
c            write(*,*) 'SOLAR:'
c            write(*,*) '   ',v1p+dvp*(jj-2),solar(jj-1)
c            write(*,*) '   ',v1p+dvp*(jj-1),solar(jj)
c            write(*,*) '   ',v1p+dvp*(jj),solar(jj+1)
c            write(*,*) '   ',v1p+dvp*(jj+1),solar(jj+2)
c
            SOLRAD(II) = (A1(JP)*SOLAR(JJ-1)+A2(JP)*SOLAR(JJ)+
     *           A3(JP)*SOLAR(JJ+1)+A4(JP)*SOLAR(JJ+2))*TRAN(II)+
     *           RADN(II)
c
c            write(*,*) 'SOLRAD:'
c            write(*,*) '   ',v1po+dvpo*(ii),solrad(ii)
c            write(*,*) ' ---'
 91      CONTINUE
      ENDIF
C
      NPL = JJ-1
C
      CALL CPUTIM (TIMSL1)
C
C     Output attenuated radiance
C
      CALL SOLOUT(V1PO,V2PO,DVPO,NLIMO,SOLRAD,LFILE,NPTS,NPANLS)
      CALL CPUTIM (TIMSL2)
      TIMOT = TIMOT+TIMSL2-TIMSL1
C
C     NPL is now location of first element to be used for next pass
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
      GO TO 60
  110 CONTINUE
C
      CALL CPUTIM (TIME1)
      TIM = TIME1-TIME
      WRITE (IPR,910) TIME1,TIM,TIMRD,TIMOT
C
  900 FORMAT ('0 THE TIME AT THE START OF SOLINT IS ',F12.3)
  905 FORMAT ('0 FILE ',I5,' MERGED WITH FILE ',I5,' ONTO FILE',
     *        I5,'  WITH XTYPE=',G15.5,/,'0 INFLAG = ',I5,4X,
     *        'IOTFLG = ',I5)
  910 FORMAT ('0 THE TIME AT THE END OF SOLINT IS ',F12.3/F12.3,
     *        ' SECS WERE REQUIRED FOR THIS SOLAR MERGE',F12.3,
     *        ' - READ - ',F12.3,' - SOLOUT - ',F12.3)
 920  FORMAT ('0 Radiance and Transmittance read in from unit',I5)
 925  FORMAT ('0 Optical Depths read in from unit',I5)
 930  FORMAT ('0 Attenuated solar radiance output to unit',I5,/)
 935  FORMAT ('0 Attenuated solar radiance + atmospheric radiance',
     *        1x,'output to unit',I5,/)
C
      END
C
C     ----------------------------------------------------------------
C
      SUBROUTINE SOLIN (V1P,V2P,DVP,NLIM,KFILE,SOLAR,KEOF)
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
      IMPLICIT DOUBLE PRECISION (V)
C
C     SUBROUTINE SOLIN inputs solar radiation from the file "SOLAR.RAD"
C     for interpolation in SOLINT.
C
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID
C
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,
     *                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4
      COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF
C
      DIMENSION PNLHDR(2),SOLAR(*)
C
      EQUIVALENCE (PNLHDR(1),V1PBF)
C
      CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
      IF (KEOF.LE.0) RETURN
      CALL BUFIN (KFILE,KEOF,SOLAR(1),NLIMBF)
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
      IMPLICIT DOUBLE PRECISION (V)
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
C
C     ----------------------------------------------------------------
C
      SUBROUTINE FLXIN (V1P,V2P,DVP,NLIM,KFILE,EM,TR,KEOF,NPANLS)         H28650
C                                                                         H28660
      IMPLICIT DOUBLE PRECISION (V)                                     ! H28670
C                                                                         H28680
C     SUBROUTINE FLXIN INPUTS OPTICAL DEPTH VALUES FROM KFILE AND         H28690
C     CALCULATES FLUX FOR THE LAYER. THIS VERSION WORKS FOR AEROSOLS.     H28700
C                                                                         H28710
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H28720
C                                                                         H28730
C               LAST MODIFICATION:    14 AUGUST 1991                      H28740
C                                                                         H28750
C                  IMPLEMENTATION:    R.D. WORSHAM                        H28760
C                                                                         H28770
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         H28780
C                                     M.J. IACONO                         H28790
C                                     R.D. WORSHAM                        H28800
C                                     J.L. MONCET                         H28810
C                                                                         H28820
C                                                                         H28830
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H28840
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H28850
C                                                                         H28860
C----------------------------------------------------------------------   H28870
C                                                                         H28880
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H28890
C                                     OFFICE OF ENERGY RESEARCH           H28900
C                                     DEPARTMENT OF ENERGY                H28910
C                                                                         H28920
C                                                                         H28930
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             H28940
C                                                                         H28950
C                                             FASCOD3                     H28960
C                                                                         H28970
C                                                                         H28980
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H28990
C                                                                         H29000
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H29010
C                                                                         H29020
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       H29030
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   H29040
     *                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    H29050
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H29060
     *              NLNGTH,KDUMY,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,        H29070
     *              NLTEFL,LNFIL4,LNGTH4                                  H29080
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           H29090
      COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF                            H29100
      COMMON /RMRG/ XKT,XKTA,XKTB,SECNT                                   H29110
C                                                                         H29120
      DIMENSION PNLHDR(2),EM(*),TR(*)                                     H29130
C                                                                         H29140
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4)                 H29150
      EQUIVALENCE (PNLHDR(1),V1PBF)                                       H29160
      EQUIVALENCE (FSCDID(4),IAERSL)                                      H29170
C                                                                         H29180
      CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)                            H29190
      IF (KEOF.LE.0) RETURN                                               H29200
      CALL BUFIN (KFILE,KEOF,TR(1),NLIMBF)                                H29210
C                                                                         H29220
C     TR CONTAINS THE OPTICAL DEPTHS AT THIS STAGE                        H29230
C                                                                         H29240
      IF (IHIRAC.EQ.4) STOP ' IHIRAC=4  FLXIN '                           H29250
C                                                                         H29260
C     EM CONTAINS THE OPTICAL DEPTH CORRECTIONS FOR NLTE AT THIS STAGE    H29270
C                                                                         H29280
      IF (NPANLS.LT.1.AND.IAERSL.EQ.0) WRITE (IPR,900)                    H29290
      IF (NPANLS.LT.1.AND.IAERSL.NE.0) WRITE (IPR,905)                    H29300
C                                                                         H29310
      EXT = 0.                                                            H29320
      ADEL = 0.                                                           H29330
      RADFN0 = 0.                                                         H29340
      RDEL = 0.                                                           H29350
      BB = 0.                                                             H29360
      BBDEL = 0.                                                          H29370
      BBA = 0.                                                            H29380
      BBDLA = 0.                                                          H29390
C                                                                         H29400
      V1P = V1PBF                                                         H29410
      V2P = V2PBF                                                         H29420
      DVP = DVPBF                                                         H29430
      NLIM = NLIMBF                                                       H29440
      VI = V1P-DVP                                                        H29450
      VIDV = VI                                                           H29460
      VIBB = VI                                                           H29470
      VAER = VI                                                           H29480
      VDUM = VI                                                           H29490
      BBLAST = -1.                                                        H29500
      BBLXTA = -2.                                                        H29510
      RDLAST = -1.                                                        H29520
      BBDUM = -4.                                                         H29530
      RDDUM = -1.                                                         H29540
      NLIM1 = 0                                                           H29550
      NLIM2 = 0                                                           H29560
C                                                                         H29570
      AA = 0.2                                                            H29580
C                                                                         H29590
      BB = BBFN(VI,DVP,V2P,XKT,VIBB,BBDEL,BBDUM)                          H29600
      IF (IAERSL.NE.0) THEN                                               H29610
         RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDDUM)                     H29620
         EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)              H29630
         IAFBB = 0                                                        H29640
         IF (VITST.LT.VAER.AND.VITST.LT.VIBB) IAFBB = 1                   H29650
         IF (VAER.LT.VITST.AND.VAER.LT.VIBB) IAFBB = 2                    H29660
      ELSE                                                                H29670
         IAFBB = -1                                                       H29680
      ENDIF                                                               H29690
C                                                                         H29700
C     - THIS SECTION TREATS THE CASE WHERE THE LAYER CONTRIBUTES          H29710
C       TO THE RADIATIVE TRANSFER ONLY ONCE                               H29720
C                                                                         H29730
C     - WITH XKTA=0 THIS ALGORITHM REVERTS TO THE ORIGINAL                H29740
C                                                                         H29750
      IF (XKTB.GT.0.) STOP ' XKTB GT 0.   FLXIN '                         H29760
C                                                                         H29770
   10 NLIM1 = NLIM2+1                                                     H29780
C                                                                         H29790
      VI = V1P+FLOAT(NLIM1-1)*DVP                                         H29800
      IF (IAFBB.EQ.-1) THEN                                               H29810
         BB = BBFN(VI,DVP,V2P,XKT,VIDV,BBDEL,BBLAST)                      H29820
         IF (XKTA.GT.0.) THEN                                             H29830
            VIBB = -VIDV                                                  H29840
            BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)                 H29850
         ELSE                                                             H29860
            BBA = BB                                                      H29870
            BBDLA = BBDEL                                                 H29880
         ENDIF                                                            H29890
         BB = BB-BBDEL                                                    H29900
         BBA = BBA-BBDLA                                                  H29910
      ELSEIF (IAFBB.EQ.0) THEN                                            H29920
         BB = BBFN(VI,DVP,V2P,XKT,VIDV,BBDEL,BBLAST)                      H29930
         IF (XKTA.GT.0.) THEN                                             H29940
            VIBB = -VIDV                                                  H29950
            BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)                 H29960
         ELSE                                                             H29970
            BBA = BB                                                      H29980
            BBDLA = BBDEL                                                 H29990
         ENDIF                                                            H30000
         BB = BB-BBDEL                                                    H30010
         BBA = BBA-BBDLA                                                  H30020
         VITST = -VIDV                                                    H30030
         RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)                    H30040
         RADFN0 = RADFN0-RDEL                                             H30050
         VAER = -VIDV                                                     H30060
         EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)              H30070
         EXT = EXT-ADEL                                                   H30080
      ELSEIF (IAFBB.EQ.1) THEN                                            H30090
         RADFN0 = RADFNI(VI,DVP,XKT,VIDV,RDEL,RDLAST)                     H30100
         RADFN0 = RADFN0-RDEL                                             H30110
         VIBB = -VIDV                                                     H30120
         BB = BBFN(VI,DVP,V2P,XKT,VIBB,BBDEL,BBLAST)                      H30130
         IF (XKTA.GT.0.) THEN                                             H30140
            VIBB = -VIDV                                                  H30150
            BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)                 H30160
         ELSE                                                             H30170
            BBA = BB                                                      H30180
            BBDLA = BBDEL                                                 H30190
         ENDIF                                                            H30200
         BB = BB-BBDEL                                                    H30210
         BBA = BBA-BBDLA                                                  H30220
         VAER = -VIDV                                                     H30230
         EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)              H30240
         EXT = EXT-ADEL                                                   H30250
      ELSEIF (IAFBB.EQ.2) THEN                                            H30260
         EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)              H30270
         EXT = EXT-ADEL                                                   H30280
         VIBB = -VIDV                                                     H30290
         BB = BBFN(VI,DVP,V2P,XKT,VIBB,BBDEL,BBLAST)                      H30300
         IF (XKTA.GT.0.) THEN                                             H30310
            VIBB = -VIDV                                                  H30320
            BBA = BBFN(VI,DVP,V2P,XKTA,VIBB,BBDLA,BBLXTA)                 H30330
         ELSE                                                             H30340
            BBA = BB                                                      H30350
            BBDLA = BBDEL                                                 H30360
         ENDIF                                                            H30370
         BB = BB-BBDEL                                                    H30380
         BBA = BBA-BBDLA                                                  H30390
         VITST = -VIDV                                                    H30400
         RADFN0 = RADFNI(VI,DVP,XKT,VITST,RDEL,RDLAST)                    H30410
         RADFN0 = RADFN0-RDEL                                             H30420
      ENDIF                                                               H30430
C                                                                         H30440
      NLIM2 = (VIDV-V1P)/DVP+1.001                                        H30450
      NLIM2 = MIN(NLIM2,NLIM)                                             H30460
C                                                                         H30470
      DO 20 I = NLIM1, NLIM2                                              H30480
         EXT = EXT+ADEL                                                   H30490
         RADFN0 = RADFN0+RDEL                                             H30500
         ODVI = SECNT*TR(I)+EXT*RADFN0                                    H30510
         BB = BB+BBDEL                                                    H30520
         BBA = BBA+BBDLA                                                  H30530
C                                                                         H30540
         XX = AA*ODVI                                                     H30550
C                                                                         H30560
         TR(I) = EXP(-ODVI)                                               H30570
         EM(I) = (1.-TR(I))*(BB+XX*BBA)/(1.+XX)                           H30580
C                                                                         H30590
   20 CONTINUE                                                            H30600
C                                                                         H30610
      IF (NLIM2.LT.NLIM) GO TO 10                                         H30620
C                                                                         H30630
      RETURN                                                              H30640
C                                                                         H30650
  900 FORMAT ('0EMISSION AND TRANSMISSION  (MOLECULAR) ')                 H30660
  905 FORMAT ('0EMISSION AND TRANSMISSION (AEROSOLS EFFECTS INCLUDED)')   H30670
C                                                                         H30680
      END                                                                 H30690
C
C     ----------------------------------------------------------------
C
      SUBROUTINE FLINIT (NPTS,MFILE,JPATHL,TBND)                          H30700
C                                                                         H30710
      IMPLICIT DOUBLE PRECISION (V)                                     ! H30720
C                                                                         H30730
C                                                                         H30740
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H30750
C                                                                         H30760
C               LAST MODIFICATION:    14 AUGUST 1991                      H30770
C                                                                         H30780
C                  IMPLEMENTATION:    R.D. WORSHAM                        H30790
C                                                                         H30800
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         H30810
C                                     M.J. IACONO                         H30820
C                                     R.D. WORSHAM                        H30830
C                                     J.L. MONCET                         H30840
C                                                                         H30850
C                                                                         H30860
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H30870
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H30880
C                                                                         H30890
C----------------------------------------------------------------------   H30900
C                                                                         H30910
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H30920
C                                     OFFICE OF ENERGY RESEARCH           H30930
C                                     DEPARTMENT OF ENERGY                H30940
C                                                                         H30950
C                                                                         H30960
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             H30970
C                                                                         H30980
C                                             FASCOD3                     H30990
C                                                                         H31000
C                                                                         H31010
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H31020
C                                                                         H31030
      COMMON NEWEM(2410),NEWTR(2410)                                      H31040
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H31050
     *              AVFIX,LAYER,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,        H31060
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H31070
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H31080
     *              EXTID(10)                                             H31090
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H31100
C                                                                         H31110
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H31120
C                                                                         H31130
      COMMON /EMIHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       H31140
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   H31150
     *                EMISIV,FSCDID(17),NMOL,LAYDUM,YI1,YID(10),LSTWDF    H31160
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP                   H31170
      COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO                                 H31180
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H31190
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,        H31200
     *              NLTEFL,LNFIL4,LNGTH4                                  H31210
      COMMON /RMRG/ XKT,XKTA,XKTB,SECNT                                   H31220
C                                                                         H31230
      CHARACTER*40 CEXT,CYID                                              H31240
C                                                                         H31250
      REAL NEWEM,NEWTR                                                    H31260
C                                                                         H31270
      DIMENSION XFILHD(2),OPNLHD(2)                                       H31280
      DIMENSION EMLAYR(2),TRLAYR(2)                                       H31290
C                                                                         H31300
      EQUIVALENCE (XFILHD(1),XID(1)) , (OPNLHD(1),V1PO)                   H31310
      EQUIVALENCE (NEWEM(1),EMLAYR(1)) , (NEWTR(1),TRLAYR(1)),            H31320
     *            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),                 H31330
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 H31340
     *            (FSCDID(16),LAYR1)                                      H31350
C                                                                         H31360
C                                                                         H31370
C   *******************************************************************   H31380
C   ***  THIS SUBROUTINE COMPUTES THE EMISSION FOR THE FIRST LAYER  ***   H31390
C   *******************************************************************   H31400
C                                                                         H31410
C     TBND IS THE BOUNDARY BLACK BODY TEMPERATUE                          H31420
C                                                                         H31430
C     IPATHL = 1 IS FOR THE DOWNWELLING FLUX CASE                         H31440
C     IPATHL = 3 IS FOR THE UPWELLING FLUX CASE                           H31450
C                                                                         H31460
      CALL CPUTIM (TIME)                                                  H31470
C                                                                         H31480
      WRITE (IPR,900) TIME                                                H31490
      NPANLS = 0                                                          H31500
      CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)                            H31510
      IF (JPATHL.GE.1) IPATHL = JPATHL                                    H31520
      PLAY = PAVE                                                         H31530
      TLAY = TAVE                                                         H31540
C                                                                         H31550
C     FOR AEROSOL RUNS, MOVE EXTID INTO YID                               H31560
C                                                                         H31570
      IF (IAERSL.GT.0) THEN                                               H31580
         WRITE (CEXT,'(10A4)') EXTID                                      H31590
         WRITE (CYID,'(5A8)') (YID(I),I=3,7)                              H31600
         CYID(19:40) = CEXT(19:40)                                        H31610
         READ (CYID,'(5A8)') (YID(I),I=3,7)                               H31620
      ENDIF                                                               H31630
C                                                                         H31640
C     READ IN DIRECTION COSINE                                            H31650
C                                                                         H31660
      READ (IRD,905) DIRCOS                                               H31670
      SECNT = 1.0/DIRCOS                                                  H31680
      SECANT = SECNT                                                      H31690
      SECNT0 = SECNT                                                      H31700
      WRITE (IPR,910) DIRCOS                                              H31710
C                                                                         H31720
C     IF BOUNDARY PROPERTIES ARE SUPPLIED, AND DOWNWARD LOOKING           H31730
C     CASE; SET IPATHL TO REFLECTED ATMOSPHERE CASE                       H31740
C                                                                         H31750
      IF (IBPROP.EQ.1.AND.IPATHL.EQ.1) IPATHL = -1                        H31760
      IEMIT = 1                                                           H31770
      FACT = 1.                                                           H31780
      TIMEM = 0.0                                                         H31790
      IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.                            H31800
      DO 10 MOL = 1, NMOL                                                 H31810
         WK(MOL) = WK(MOL)*FACT                                           H31820
   10 CONTINUE                                                            H31830
      WBROAD = WBROAD*FACT                                                H31840
      LAYR1 = LAYER                                                       H31850
      WRITE (IPR,915) LAYR1,LAYER,KFILE,MFILE                             H31860
      CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)                                H31870
      DVXM = DV                                                           H31880
      XKT = TAVE/RADCN2                                                   H31890
      XKTBND = TBND/RADCN2                                                H31900
      IF (IPATHL.EQ.-1) STOP ' IPATH=-1 '                                 H31910
      IF (IPATHL.EQ.0) STOP ' IPATH=0 '                                   H31920
      IF (IPATHL.EQ.1) THEN                                               H31930
         XKTA = TZL/RADCN2                                                H31940
         XKTB = 0.                                                        H31950
      ENDIF                                                               H31960
      IF (IPATHL.EQ.2) STOP ' IPATH=2 '                                   H31970
      IF (IPATHL.EQ.3) THEN                                               H31980
         XKTA = TZU/RADCN2                                                H31990
         XKTB = 0.                                                        H32000
      ENDIF                                                               H32010
      WRITE (IPR,920) IPATHL,IANT                                         H32020
C                                                                         H32030
   20 CONTINUE                                                            H32040
C                                                                         H32050
      CALL CPUTIM (TIMEM1)                                                H32060
      CALL FLXIN (V1PO,V2PO,DVPO,NLIMO,KFILE,EMLAYR,TRLAYR,KEOF,NPANLS)   H32070
      CALL CPUTIM (TIMEM2)                                                H32080
      TIMEM = TIMEM+TIMEM2-TIMEM1                                         H32090
      IF (KEOF.LE.0) GO TO 50                                             H32100
      VI = V1PO-DVPO                                                      H32110
      VIDVBD = VI                                                         H32120
      VIDVEM = VI                                                         H32130
      VIDVRF = VI                                                         H32140
      BBLAST = -1.                                                        H32150
      EMLAST = -1.                                                        H32160
      IF ((IPATHL.EQ.3).AND.(TBND.GT.0.)) THEN                            H32170
C                                                                         H32180
         NLIM1 = 0                                                        H32190
         NLIM2 = 0                                                        H32200
         EMDUM = 0.                                                       H32210
         BBDUM = 0.                                                       H32220
         EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMDUM)                      H32230
         BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDVBD,BBDEL,BBDUM)                H32240
         IEMBB = 0                                                        H32250
         IF (VIDVBD.GT.VIDVEM) IEMBB = 1                                  H32260
C                                                                         H32270
   30    NLIM1 = NLIM2+1                                                  H32280
C                                                                         H32290
         VI = V1PO+FLOAT(NLIM1-1)*DVPO                                    H32300
         IF (IEMBB.EQ.0) THEN                                             H32310
            BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDV,BBDEL,BBLAST)              H32320
            BB = BB-BBDEL                                                 H32330
            VIDVEM = -VIDV                                                H32340
            EMISIV = EMISFN(VI,DVPO,VIDVEM,EMDEL,EMLAST)                  H32350
            EMISIV = EMISIV-EMDEL                                         H32360
         ELSE                                                             H32370
            EMISIV = EMISFN(VI,DVPO,VIDV,EMDEL,EMLAST)                    H32380
            EMISIV = EMISIV-EMDEL                                         H32390
            VIDVBD = -VIDV                                                H32400
            BB = BBFN(VI,DVPO,V2PO,XKTBND,VIDVBD,BBDEL,BBLAST)            H32410
            BB = BB-BBDEL                                                 H32420
         ENDIF                                                            H32430
C                                                                         H32440
         IF (VIDV.GE.9.E+4) THEN 
            NLIM2 = NLIMO+1
         ELSE
            NLIM2 = (VIDV-V1PO)/DVPO+1.001                                H32450
         ENDIF
         NLIM2 = MIN(NLIM2,NLIMO)                                         H32460
C                                                                         H32470
         DO 40 J = NLIM1, NLIM2                                           H32480
            EMISIV = EMISIV+EMDEL                                         H32490
            BB = BB+BBDEL                                                 H32500
            NEWEM(J) = EMLAYR(J)+TRLAYR(J)*BB*EMISIV                      H32510
   40    CONTINUE                                                         H32520
C                                                                         H32530
         IF (NLIM2.LT.NLIMO) GO TO 30                                     H32540
C                                                                         H32550
      ENDIF                                                               H32560
      CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,NEWEM,NEWTR,MFILE,NPTS,NPANLS)     H32570
      GO TO 20                                                            H32580
   50 CALL CPUTIM (TIME1)                                                 H32590
      TIME = TIME1-TIME                                                   H32600
      WRITE (IPR,925) TIME,TIMEM                                          H32610
C                                                                         H32620
      RETURN                                                              H32630
C                                                                         H32640
  900 FORMAT (' TIME AT THE START OF --FLINIT-- ',F10.3)                  H32650
  905 FORMAT (F10.8)                                                      H32660
  910 FORMAT ('0 ***** DIR. COSINE ***** ',/,7X,F10.8)                    H32670
  915 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,                  H32680
     *        '0 INPUT FILE =',I5,' OUTPUT FILE =',I5)                    H32690
  920 FORMAT ('0 IPATHL AND IANT',2I5)                                    H32700
  925 FORMAT (' TIME REQUIRED FOR --FLINIT-- ',F10.3,' --FLXIN-- ',       H32710
     *        F10.3)                                                      H32720
C                                                                         H32730
      END                                                                 H32740
C
C     ----------------------------------------------------------------
C
      SUBROUTINE FLUXUP (NPTS,LFILE,MFILE,JPATHL,TBND)                    H32750
C                                                                         H32760
      IMPLICIT DOUBLE PRECISION (V)                                     ! H32770
C                                                                         H32780
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H32790
C                                                                         H32800
C               LAST MODIFICATION:    14 AUGUST 1991                      H32810
C                                                                         H32820
C                  IMPLEMENTATION:    R.D. WORSHAM                        H32830
C                                                                         H32840
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         H32850
C                                     M.J. IACONO                         H32860
C                                     R.D. WORSHAM                        H32870
C                                     J.L. MONCET                         H32880
C                                                                         H32890
C                                                                         H32900
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H32910
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H32920
C                                                                         H32930
C----------------------------------------------------------------------   H32940
C                                                                         H32950
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H32960
C                                     OFFICE OF ENERGY RESEARCH           H32970
C                                     DEPARTMENT OF ENERGY                H32980
C                                                                         H32990
C                                                                         H33000
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             H33010
C                                                                         H33020
C                                             FASCOD3                     H33030
C                                                                         H33040
C                                                                         H33050
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H33060
C                                                                         H33070
      COMMON RADN(2410),TRAN(2410),RADO(0:5000)                           H33080
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H33090
     *              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       H33100
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H33110
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H33120
     *              EXTID(10)                                             H33130
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H33140
C                                                                         H33150
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H33160
C                                                                         H33170
      COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),        H33180
     *               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,    H33190
     *               EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF     H33200
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP                   H33210
      COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO                                 H33220
      COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS      H33230
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H33240
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,        H33250
     *              NLTEFL,LNFIL4,LNGTH4                                  H33260
      COMMON /XME/ TRAO(0:5000)                                           H33270
      COMMON /RMRG/ XKT,XKTA,XKTB,SECNT                                   H33280
C                                                                         H33290
      DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)                             H33300
      DIMENSION A1(10),A2(10),A3(10),A4(10)                               H33310
      DIMENSION RADLYR(2),TRALYR(2),RADOI(2),TRAOI(2)                     H33320
      DIMENSION WKSAV(35)                                                 H33330
C                                                                         H33340
      CHARACTER*40 CYID                                                   H33350
C                                                                         H33360
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                   H33370
     *            (OPNLHD(1),V1PO)                                        H33380
      EQUIVALENCE (RADO(1),RADOI(1)) , (TRAO(1),TRAOI(1))                 H33390
      EQUIVALENCE (RADN(1),RADLYR(1)) , (TRAN(1),TRALYR(1)),              H33400
     *            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),                 H33410
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 H33420
     *            (FSCDID(16),LAYR1)                                      H33430
C                                                                         H33440
      DATA NDIM / 2410 /,ND2 / 5000 /                                     H33450
C                                                                         H33460
C                                                                         H33470
C                                                                         H33480
C     ************************************************************        H33490
C     ****** THIS SUBROUTINE DOES LAYER MERGE FOR RADIANCE  ******        H33500
C     ************************************************************        H33510
C                                                                         H33520
C     IPATHL = 3 IS FOR THE UPWELLING FLUX CASE                           H33530
C                                                                         H33540
      CALL CPUTIM (TIME)                                                  H33550
      WRITE (IPR,900) TIME                                                H33560
      NPANLS = 0                                                          H33570
      TIMEM = 0.0                                                         H33580
      TIMRD = 0.0                                                         H33590
      TIMOT = 0.0                                                         H33600
C                                                                         H33610
      CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)                            H33620
      SECNT = SECANT                                                      H33630
      LAY1SV = LAYR1                                                      H33640
      DVL = DV                                                            H33650
      PL = PAVE                                                           H33660
      TL = TAVE                                                           H33670
      WTOTL = 0.                                                          H33680
C                                                                         H33690
      DO 10 MOL = 1, NMOL                                                 H33700
         WTOTL = WTOTL+WK(MOL)                                            H33710
         WKSAV(MOL) = WK(MOL)                                             H33720
   10 CONTINUE                                                            H33730
C                                                                         H33740
      WTOTL = WTOTL+WBROAD                                                H33750
      WN2SAV = WBROAD                                                     H33760
C                                                                         H33770
C     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)                 H33780
C                                                                         H33790
      IF (IAERSL.GT.0) WRITE (CYID,'(5A8)') (YID(I),I=3,7)                H33800
      CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)                            H33810
      IF (IAERSL.GT.0) READ (CYID,'(5A8)') (YID(I),I=3,7)                 H33820
C                                                                         H33830
      IF (JPATHL.GE.1) IPATHL = JPATHL                                    H33840
      PLAY = PAVE                                                         H33850
      TLAY = TAVE                                                         H33860
      TAVK = TAVE                                                         H33870
      DVK = DV                                                            H33880
      FACT = 1.                                                           H33890
C                                                                         H33900
      IF (DVL.EQ.DVK) THEN                                                H33910
         ITYPE = 0                                                        H33920
      ELSEIF (DVL.GT.DVK) THEN                                            H33930
         ITYPE = DVK/(DVL-DVK)+0.5                                        H33940
      ELSE                                                                H33950
C                                                                         H33960
C     DVL.LT.DVK                                                          H33970
C                                                                         H33980
         ITYPE = -INT(DVL/(DVK-DVL)+0.5)                                  H33990
      ENDIF                                                               H34000
      IF (ITYPE.LT.0) STOP ' FLUXUP; ITYPE LT 0 '                         H34010
C                                                                         H34020
      WTOTK = 0.                                                          H34030
      LAYR1 = LAY1SV                                                      H34040
      WRITE (IPR,905) LAYR1,LAYER,KFILE,LFILE,MFILE                       H34050
      IEMIT = 1                                                           H34060
      DO 20 MOL = 1, NMOL                                                 H34070
         WTOTK = WTOTK+WK(MOL)*FACT                                       H34080
         WK(MOL) = WK(MOL)*FACT+WKSAV(MOL)                                H34090
   20 CONTINUE                                                            H34100
      WTOTK = WTOTK+WBROAD*FACT                                           H34110
      WBROAD = WBROAD*FACT+WN2SAV                                         H34120
      PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)                          H34130
      TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)                          H34140
      SECANT = SECNT                                                      H34150
C                                                                         H34160
C     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES               H34170
C                                                                         H34180
      CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)                                H34190
      DVXM = DV                                                           H34200
      XKT = TAVK/RADCN2                                                   H34210
C                                                                         H34220
      WRITE (IPR,910) IPATHL,IANT                                         H34230
C                                                                         H34240
      IF (IPATHL.EQ.3) THEN                                               H34250
         XKTA = TZU/RADCN2                                                H34260
         XKTB = 0.                                                        H34270
      ELSE                                                                H34280
         STOP ' FLUXUP; IPATHL '                                          H34290
      ENDIF                                                               H34300
C                                                                         H34310
      ATYPE = ITYPE                                                       H34320
      LL = ITYPE+1                                                        H34330
      AP = 1.0/(ATYPE+1.0)                                                H34340
C                                                                         H34350
C     A1, A2, A3 AND A4 ARE THE CONSTANTS                                 H34360
C     FOR THE LAGRANGE 4 POINT INTERPOLATION                              H34370
C                                                                         H34380
      DO 30 JPG = 1, ITYPE                                                H34390
         APG = JPG                                                        H34400
         IPL = JPG+1                                                      H34410
         P = 1.0-(AP*APG)                                                 H34420
         A1(IPL) = -P*(P-1.0)*(P-2.0)/6.0                                 H34430
         A2(IPL) = (P**2-1.0)*(P-2.0)*0.5                                 H34440
         A3(IPL) = -P*(P+1.0)*(P-2.0)*0.5                                 H34450
         A4(IPL) = P*(P**2-1.0)/6.0                                       H34460
   30 CONTINUE                                                            H34470
C                                                                         H34480
C     *** BEGINNING OF LOOP THAT DOES MERGE  ***                          H34490
C                                                                         H34500
      NPE = 0                                                             H34510
      RADO(0) = 0.0                                                       H34520
      TRAO(0) = 0.0                                                       H34530
      V1PO = 0.0                                                          H34540
      V2PO = 0.0                                                          H34550
      DVPO = 0.0                                                          H34560
C                                                                         H34570
   40 CONTINUE                                                            H34580
C                                                                         H34590
      CALL CPUTIM (TIMEM1)                                                H34600
      CALL FLXIN (V1P,V2P,DVP,NLIM,KFILE,RADLYR,TRALYR,KEOF,NPANLS)       H34610
      CALL CPUTIM (TIMEM2)                                                H34620
      TIMEM = TIMEM+TIMEM2-TIMEM1                                         H34630
      IF (KEOF.LE.0) GO TO 80                                             H34640
      II = 1                                                              H34650
C                                                                         H34660
      IF (V2PO.LE.V2P+DVPO) THEN                                          H34670
   50    CALL CPUTIM (TIMEM1)                                             H34680
         CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)                         H34690
         IF (LEOF.LE.0) GO TO 60                                          H34700
         CALL BUFIN (LFILE,LEOF,RADOI(NPE+1),NLIMO)                       H34710
         CALL BUFIN (LFILE,LEOF,TRAOI(NPE+1),NLIMO)                       H34720
         CALL CPUTIM (TIMEM2)                                             H34730
         TIMRD = TIMRD+TIMEM2-TIMEM1                                      H34740
         NPE = NLIMO+NPE                                                  H34750
         IF (V2PO.LE.V2P+DVPO) GO TO 50                                   H34760
      ENDIF                                                               H34770
C                                                                         H34780
C     ZERO POINT OF FIRST PANEL                                           H34790
C                                                                         H34800
   60 IF (RADO(0).EQ.0.0.AND.TRAO(0).EQ.0.0) THEN                         H34810
         RADO(0) = RADO(1)                                                H34820
         TRAO(0) = TRAO(1)                                                H34830
      ENDIF                                                               H34840
C                                                                         H34850
C     END POINT OF LAST PANEL                                             H34860
C                                                                         H34870
      IF (V2PO+DVPO.GE.V2) THEN                                           H34880
         RADO(NPE+1) = RADO(NPE)                                          H34890
         RADO(NPE+2) = RADO(NPE)                                          H34900
         TRAO(NPE+1) = TRAO(NPE)                                          H34910
         TRAO(NPE+2) = TRAO(NPE)                                          H34920
      ENDIF                                                               H34930
C                                                                         H34940
      NPL = 1                                                             H34950
C                                                                         H34960
C     NPL IS LOCATION OF FIRST ELEMENT ON ARRAYS RADO AND TRAO            H34970
C                                                                         H34980
      CALL FLUXNN (RADN,TRAN,RADO,TRAO,NLIM,NDIM,ND2,V1P,DVP,IPATHL,      H34990
     *             A1,A2,A3,A4,LL,NPL)                                    H35000
C                                                                         H35010
      CALL CPUTIM (TIMEM1)                                                H35020
C                                                                         H35030
      IF (TBND.GT.0.) CALL EMBND (V1P,V2P,DVP,NLIM,RADN,TRAN,TBND)        H35040
C                                                                         H35050
      CALL EMOUT (V1P,V2P,DVP,NLIM,RADN,TRAN,MFILE,NPTS,NPANLS)           H35060
      CALL CPUTIM (TIMEM2)                                                H35070
      TIMOT = TIMOT+TIMEM2-TIMEM1                                         H35080
C                                                                         H35090
C     NPL IS NOW LOCATION OF FIRST ELEMENT TO BE USED FOR NEXT PASS       H35100
C                                                                         H35110
      IPL = -1                                                            H35120
      DO 70 NL = NPL, NPE                                                 H35130
         IPL = IPL+1                                                      H35140
         RADO(IPL) = RADO(NL)                                             H35150
         TRAO(IPL) = TRAO(NL)                                             H35160
   70 CONTINUE                                                            H35170
C                                                                         H35180
      NPE = IPL                                                           H35190
C                                                                         H35200
      GO TO 40                                                            H35210
   80 CONTINUE                                                            H35220
C                                                                         H35230
      CALL CPUTIM (TIME1)                                                 H35240
      TIM = TIME1-TIME                                                    H35250
      WRITE (IPR,915) TIME1,TIM,TIMEM,TIMRD,TIMOT                         H35260
C                                                                         H35270
      RETURN                                                              H35280
C                                                                         H35290
  900 FORMAT ('0 THE TIME AT THE START OF FLUXUP IS ',F12.3)              H35300
  905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,     H35310
     *        ' MERGED WITH FILE ',I5,' ONTO FILE',I5)                    H35320
  910 FORMAT ('0 IPATHL AND IANT',2I5)                                    H35330
  915 FORMAT ('0 THE TIME AT THE END OF FLUXUP IS ',F12.3,/,F12.3,        H35340
     *        ' SECS WERE REQUIRED FOR THIS MERGE  - FLXIN - ',F12.3,     H35350
     *        ' - READ - ',F12.3,' - EMOUT - ',F12.3)                     H35360
C                                                                         H35370
      END                                                                 H35380
C
C     ----------------------------------------------------------------
C
      SUBROUTINE FLUXNN (RADLYR,TRALYR,RADO,TRAO,NLIM,NDIM,ND2,V1P,DVP,   H35390
     *                  IPATHL,A1,A2,A3,A4,LL,NPL)                        H35400
C                                                                         H35410
      IMPLICIT DOUBLE PRECISION (V)                                     ! H35420
C                                                                         H35430
C     THIS SUBROUTINE CALCULATES THE NEW RADIANCE AND TRANSMISSION        H35440
C                                                                         H35450
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H35460
C                                                                         H35470
C               LAST MODIFICATION:    14 AUGUST 1991                      H35480
C                                                                         H35490
C                  IMPLEMENTATION:    R.D. WORSHAM                        H35500
C                                                                         H35510
C                       ALGORITHM:    R.D. WORSHAM                        H35520
C                                     S.A. CLOUGH                         H35530
C                                     M.J. IACONO                         H35540
C                                     J.L. MONCET                         H35550
C                                                                         H35560
C                                                                         H35570
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H35580
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H35590
C                                                                         H35600
C----------------------------------------------------------------------   H35610
C                                                                         H35620
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H35630
C                                     OFFICE OF ENERGY RESEARCH           H35640
C                                     DEPARTMENT OF ENERGY                H35650
C                                                                         H35660
C                                                                         H35670
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H35680
C                                                                         H35690
                                                                          H35700
      DIMENSION RADLYR(NDIM),TRALYR(NDIM),RADO(0:ND2),TRAO(0:ND2),        H35710
     *          A1(*),A2(*),A3(*),A4(*)                                   H35720
C                                                                         H35730
      LLM1 = LL-1                                                         H35740
      LLM1 = MAX(LLM1,1)                                                  H35750
C                                                                         H35760
C     LOOPING OVER POINTS WITH SAME WEIGHTS                               H35770
C                                                                         H35780
      DO 30 NL = 1, LL                                                    H35790
         IPL = (NPL+NL-1)-LLM1                                            H35800
         IF (NL.GT.1) IPL = IPL-1                                         H35810
C                                                                         H35820
         IF (NL.EQ.1) THEN                                                H35830
C                                                                         H35840
C     EXACT FREQUENCY - NO INTERPOLATION                                  H35850
C                                                                         H35860
            DO 10 I = NL, NLIM, LL                                        H35870
               IPL = IPL+LLM1                                             H35880
               RADLYR(I) = TRALYR(I)*RADO(IPL)+RADLYR(I)                  H35890
               TRALYR(I) = TRALYR(I)*TRAO(IPL)                            H35900
   10       CONTINUE                                                      H35910
C                                                                         H35920
C     NOT EXACT FREQUENCY - INTERPOLATE RESULT                            H35930
C                                                                         H35940
         ELSE                                                             H35950
C                                                                         H35960
            A1N = A1(NL)                                                  H35970
            A2N = A2(NL)                                                  H35980
            A3N = A3(NL)                                                  H35990
            A4N = A4(NL)                                                  H36000
C                                                                         H36010
            DO 20 I = NL, NLIM, LL                                        H36020
               IPL = IPL+LLM1                                             H36030
C                                                                         H36040
C     INTERPOLATE THE OLD RADIANCE                                        H36050
C                                                                         H36060
               RADLYR(I) = TRALYR(I)*(A1N*RADO(IPL-1)+A2N*RADO(IPL)+      H36070
     *                     A3N*RADO(IPL+1)+A4N*RADO(IPL+2))+RADLYR(I)     H36080
C                                                                         H36090
C     INTERPOLATE THE OLD TRANSMISSION                                    H36100
C                                                                         H36110
               TRALYR(I) = TRALYR(I)*(A1N*TRAO(IPL-1)+A2N*TRAO(IPL)+      H36120
     *                     A3N*TRAO(IPL+1)+A4N*TRAO(IPL+2))               H36130
   20       CONTINUE                                                      H36140
C                                                                         H36150
C                                                                         H36160
         ENDIF                                                            H36170
C                                                                         H36180
   30 CONTINUE                                                            H36190
C                                                                         H36200
      NPL = IPL                                                           H36210
C                                                                         H36220
      RETURN                                                              H36230
C                                                                         H36240
      END                                                                 H36250
C
C     ----------------------------------------------------------------
C
      SUBROUTINE FLUXDN (NPTS,LFILE,MFILE,JPATHL,TBND)                    H36260
C                                                                         H36270
      IMPLICIT DOUBLE PRECISION (V)                                     ! H36280
C                                                                         H36290
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H36300
C                                                                         H36310
C               LAST MODIFICATION:    14 AUGUST 1991                      H36320
C                                                                         H36330
C                  IMPLEMENTATION:    R.D. WORSHAM                        H36340
C                                                                         H36350
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         H36360
C                                     M.J. IACONO                         H36370
C                                     R.D. WORSHAM                        H36380
C                                     J.L. MONCET                         H36390
C                                                                         H36400
C                                                                         H36410
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         H36420
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          H36430
C                                                                         H36440
C----------------------------------------------------------------------   H36450
C                                                                         H36460
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     H36470
C                                     OFFICE OF ENERGY RESEARCH           H36480
C                                     DEPARTMENT OF ENERGY                H36490
C                                                                         H36500
C                                                                         H36510
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             H36520
C                                                                         H36530
C                                             FASCOD3                     H36540
C                                                                         H36550
C                                                                         H36560
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   H36570
C                                                                         H36580
      COMMON RADN(2410),TRAN(2410),RADLYR(-1:4818)                        H36590
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H36600
     *              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       H36610
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H36620
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H36630
     *              EXTID(10)                                             H36640
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           H36650
C                                                                         H36660
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H36670
C                                                                         H36680
      COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),        H36690
     *               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,    H36700
     *               EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF     H36710
      COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO                                 H36720
      COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS      H36730
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H36740
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,        H36750
     *              NLTEFL,LNFIL4,LNGTH4                                  H36760
      COMMON /XMI/ TRALYR(-1:4818)                                        H36770
      COMMON /RMRG/ XKT,XKTA,XKTB,SECNT                                   H36780
C                                                                         H36790
      DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)                             H36800
      DIMENSION A1(0:100),A2(0:100),A3(0:100),A4(0:100)                   H36810
      DIMENSION RADO(2),TRAO(2)                                           H36820
      DIMENSION WKSAV(35)                                                 H36830
C                                                                         H36840
      CHARACTER*40 CYID                                                   H36850
C                                                                         H36860
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHDR(1),V1P),                   H36870
     *            (OPNLHD(1),V1PO)                                        H36880
      EQUIVALENCE (RADN(1),RADO(1)) , (TRAN(1),TRAO(1)),                  H36890
     *            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),                 H36900
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 H36910
     *            (FSCDID(16),LAYR1)                                      H36920
C                                                                         H36930
C     ************************************************************        H36940
C     ****** THIS SUBROUTINE DOES LAYER MERGE FOR EMISSION  ******        H36950
C     ****** USING FOUR POINT GENERAL INTERPOLATION         ******        H36960
C     ************************************************************        H36970
C                                                                         H36980
      CALL CPUTIM (TIME)                                                  H36990
      WRITE (IPR,900) TIME                                                H37000
      NPANLS = 0                                                          H37010
      TIMEM = 0.0                                                         H37020
      TIMRD = 0.0                                                         H37030
      TIMTB = 0.0                                                         H37040
      TIMOT = 0.0                                                         H37050
C                                                                         H37060
      CALL BUFIN (LFILE,LEOF,XFILHD(1),NFHDRF)                            H37070
      SECNT = SECANT                                                      H37080
      DVL = DV                                                            H37090
      LAY1SV = LAYR1                                                      H37100
      PL = PAVE                                                           H37110
      TL = TAVE                                                           H37120
      WTOTL = 0.                                                          H37130
      DO 10 MOL = 1, NMOL                                                 H37140
         WTOTL = WTOTL+WK(MOL)                                            H37150
         WKSAV(MOL) = WK(MOL)                                             H37160
   10 CONTINUE                                                            H37170
      WTOTL = WTOTL+WBROAD                                                H37180
      WN2SAV = WBROAD                                                     H37190
C                                                                         H37200
C     FOR AEROSOL RUNS, MOVE YID (LFILE) INTO YID (MFILE)                 H37210
C                                                                         H37220
      IF (IAERSL.GT.0) WRITE (CYID,'(5A8)') (YID(I),I=3,7)                H37230
      CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)                            H37240
      IF (IAERSL.GT.0) READ (CYID,'(5A8)') (YID(I),I=3,7)                 H37250
      IF (JPATHL.GE.1) IPATHL = JPATHL                                    H37260
      PLAY = PAVE                                                         H37270
      TLAY = TAVE                                                         H37280
C                                                                         H37290
      IF (IPATHL.NE.1) STOP ' FLUXDN: IPATHL .NE. 1 '                     H37300
C                                                                         H37310
      XKT = TAVE/RADCN2                                                   H37320
      XKTA = TZL/RADCN2                                                   H37330
      XKTB = 0.                                                           H37340
      DVK = DV                                                            H37350
      LAYR1 = LAY1SV                                                      H37360
      FACT = 1.                                                           H37370
      IF (IPATHL.EQ.2.AND.IANT.EQ.0) FACT = 2.                            H37380
      ATYPE = 9.999E09                                                    H37390
      IF (DVK.EQ.DVL) ATYPE = 0.                                          H37400
      IF (DVL.GT.DVK) ATYPE = DVK/(DVL-DVK)+0.5                           H37410
      IF (DVL.LT.DVK) ATYPE = -DVL/(DVK-DVL)-0.5                          H37420
C                                                                         H37430
C     IF (ATYPE .GT. 0) STOP  ' FLUXDN; ATYPE GT 0 '                      H37440
C                                                                         H37450
      WTOTK = 0.                                                          H37460
      WRITE (IPR,905) LAYR1,LAYER,KFILE,LFILE,MFILE,ATYPE                 H37470
      IEMIT = 1                                                           H37480
      DO 20 MOL = 1, NMOL                                                 H37490
         WTOTK = WTOTK+WK(MOL)*FACT                                       H37500
         WK(MOL) = WK(MOL)*FACT+WKSAV(MOL)                                H37510
   20 CONTINUE                                                            H37520
      WTOTK = WTOTK+WBROAD*FACT                                           H37530
      WBROAD = WBROAD*FACT+WN2SAV                                         H37540
      PAVE = (PL*WTOTL+PAVE*WTOTK)/(WTOTL+WTOTK)                          H37550
      TAVE = (TL*WTOTL+TAVE*WTOTK)/(WTOTL+WTOTK)                          H37560
      SECANT = SECNT                                                      H37570
      DV = DVL                                                            H37580
C                                                                         H37590
C     WK IS NOW THE ACCUMULATED SUM OF THE COLUMN DENSITIES               H37600
C                                                                         H37610
      CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)                                H37620
      DVXM = DV                                                           H37630
C                                                                         H37640
      IF (ATYPE.EQ.0.) THEN                                               H37650
C                                                                         H37660
C     1/1 RATIO ONLY                                                      H37670
C                                                                         H37680
C                                                                         H37690
   30    CONTINUE                                                         H37700
         CALL CPUTIM (TIMEM1)                                             H37710
         CALL FLXIN (V1P,V2P,DVP,NLIM,KFILE,RADLYR(1),TRALYR(1),KEOF,     H37720
     *               NPANLS)                                              H37730
         CALL CPUTIM (TIMEM2)                                             H37740
         TIMEM = TIMEM+TIMEM2-TIMEM1                                      H37750
         IF (KEOF.LE.0) GO TO 110                                         H37760
         CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)                         H37770
         CALL BUFIN (LFILE,LEOF,RADO(1),NLIMO)                            H37780
         CALL BUFIN (LFILE,LEOF,TRAO(1),NLIMO)                            H37790
         CALL CPUTIM (TIMEM3)                                             H37800
         TIMRD = TIMRD+TIMEM3-TIMEM2                                      H37810
         DO 40 I = 1, NLIM                                                H37820
            RADN(I) = RADO(I)*TRALYR(I)+RADLYR(I)                         H37830
            TRAN(I) = TRALYR(I)*TRAO(I)                                   H37840
   40    CONTINUE                                                         H37850
         CALL CPUTIM (TIMEM1)                                             H37860
         IF (TBND.GT.0.)                                                  H37870
     *        CALL EMBND (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,TBND)            H37880
C                                                                         H37890
         CALL CPUTIM (TIMEM2)                                             H37900
         TIMTB = TIMTB+TIMEM2-TIMEM1                                      H37910
         CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,MFILE,NPTS,NPANLS)    H37920
         CALL CPUTIM (TIMEM3)                                             H37930
         TIMOT = TIMOT+TIMEM3-TIMEM2                                      H37940
         GO TO 30                                                         H37950
C                                                                         H37960
      ENDIF                                                               H37970
C                                                                         H37980
C     ALL RATIOS EXCEPT 1/1                                               H37990
C                                                                         H38000
      DO 50 JP = 0, 100                                                   H38010
         APG = JP                                                         H38020
         P = 0.01*APG                                                     H38030
C                                                                         H38040
C     THE FOLLOW ARE THE CONSTANTS FOR THE LAGRANGE 4 POINT               H38050
C     INTERPOLATION                                                       H38060
C                                                                         H38070
         A1(JP) = -P*(P-1.0)*(P-2.0)/6.0                                  H38080
         A2(JP) = (P**2-1.0)*(P-2.0)*0.5                                  H38090
         A3(JP) = -P*(P+1.0)*(P-2.0)*0.5                                  H38100
         A4(JP) = P*(P**2-1.0)/6.0                                        H38110
   50 CONTINUE                                                            H38120
C                                                                         H38130
C     *** BEGINNING OF LOOP THAT DOES MERGE  ***                          H38140
C                                                                         H38150
      NPE = 0                                                             H38160
      RADLYR(0) = 0.0                                                     H38170
      TRALYR(0) = 0.0                                                     H38180
      V1P = 0.0                                                           H38190
      V2P = 0.0                                                           H38200
      DVP = 0.0                                                           H38210
      V1PO = 0.0                                                          H38220
      V2PO = 0.0                                                          H38230
      DVPO = 0.0                                                          H38240
      KEOF = 1                                                            H38245
C                                                                         H38250
   60 CONTINUE                                                            H38260
C                                                                         H38270
      CALL CPUTIM (TIMEM1)                                                H38280
      CALL BUFIN (LFILE,LEOF,OPNLHD(1),NPHDRF)                            H38290
      IF (LEOF.LE.0) GO TO 110                                            H38300
      CALL BUFIN (LFILE,LEOF,RADO(1),NLIMO)                               H38310
      CALL BUFIN (LFILE,LEOF,TRAO(1),NLIMO)                               H38320
      CALL CPUTIM (TIMEM2)                                                H38330
      TIMRD = TIMRD+TIMEM2-TIMEM1                                         H38340
      II = 1                                                              H38350
C                                                                         H38360
      IF (V2P.LE.V2PO+DVP .AND. KEOF.GT.0) THEN                           H38370
   70    CALL CPUTIM (TIMEM2)                                             H38380
         CALL FLXIN (V1P,V2P,DVP,NLIM,KFILE,RADLYR(NPE+1),                H38390
     *               TRALYR(NPE+1),KEOF,NPANLS)                           H38400
         CALL CPUTIM (TIMEM3)                                             H38410
         TIMEM = TIMEM+TIMEM3-TIMEM2                                      H38420
         IF (KEOF.LE.0) GO TO 80                                          H38430
         V1P = V1P-FLOAT(NPE)*DVP                                         H38440
         NPE = NLIM+NPE                                                   H38450
         IF (V2P.LE.V2PO+DVP) GO TO 70                                    H38460
      ENDIF                                                               H38470
C                                                                         H38480
C     ZERO POINT OF FIRST PANEL                                           H38490
C                                                                         H38500
   80 IF (RADLYR(0).EQ.0.0.AND.TRALYR(0).EQ.0.0) THEN                     H38510
         TRALYR(-1) = TRALYR(1)                                           H38520
         RADLYR(-1) = RADLYR(1)                                           H38530
         TRALYR(0) = TRALYR(1)                                            H38540
         RADLYR(0) = RADLYR(1)                                            H38550
      ENDIF                                                               H38560
C                                                                         H38570
C     END POINT OF LAST PANEL                                             H38580
C                                                                         H38590
      IF (V2P+DVP.GE.V2) THEN                                             H38600
         TRALYR(NPE+1) = TRALYR(NPE)                                      H38610
         RADLYR(NPE+1) = RADLYR(NPE)                                      H38620
         TRALYR(NPE+2) = TRALYR(NPE)                                      H38630
         RADLYR(NPE+2) = RADLYR(NPE)                                      H38640
      ENDIF                                                               H38650
C                                                                         H38660
C     NPL IS LOCATION OF FIRST ELEMENT ON ARRAYS RADO AND TRAO            H38670
C                                                                         H38680
      NPL = 1                                                             H38690
C                                                                         H38700
      RATDV = DVL/DVK                                                     H38710
C                                                                         H38720
C     FJJ IS OFFSET BY 2. FOR ROUNDING PURPOSES                           H38730
C                                                                         H38740
      FJ1DIF = (V1PO-V1P)/DVP+1.+2.                                       H38750
C                                                                         H38760
C     ***** BEGINNING OF LOOP THAT DOES MERGE  *****                      H38770
C                                                                         H38780
      DO 90 II = 1, NLIMO                                                 H38790
C                                                                         H38800
         FJJ = FJ1DIF+RATDV*FLOAT(II-1)                                   H38810
         JJ = IFIX(FJJ)-2                                                 H38820
C                                                                         H38830
         JP = (FJJ-FLOAT(JJ))*100.-199.5                                  H38840
C                                                                         H38850
C     INTERPOLATE THE OLD TRANSMISSION                                    H38860
C                                                                         H38870
         TRNLYR = A1(JP)*TRALYR(JJ-1)+A2(JP)*TRALYR(JJ)+                  H38880
     *            A3(JP)*TRALYR(JJ+1)+A4(JP)*TRALYR(JJ+2)                 H38890
C                                                                         H38900
C     INTERPOLATE THE OLD EMISSION                                        H38910
C                                                                         H38920
         RADN(II) = RADO(II)*TRNLYR+(A1(JP)*RADLYR(JJ-1)+                 H38930
     *              A2(JP)*RADLYR(JJ)+A3(JP)*RADLYR(JJ+1)+                H38940
     *              A4(JP)*RADLYR(JJ+2))                                  H38950
C                                                                         H38960
         TRAN(II) = TRNLYR*TRAO(II)                                       H38970
C                                                                         H38980
   90 CONTINUE                                                            H38990
C                                                                         H39000
      NPL = JJ-1                                                          H39010
C                                                                         H39020
      CALL CPUTIM (TIMEM1)                                                H39030
      IF (TBND.GT.0.) CALL EMBND (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,TBND)    H39040
C                                                                         H39050
      CALL CPUTIM (TIMEM2)                                                H39060
      TIMTB = TIMTB+TIMEM2-TIMEM1                                         H39070
      CALL EMOUT (V1PO,V2PO,DVPO,NLIMO,RADN,TRAN,MFILE,NPTS,NPANLS)       H39080
      CALL CPUTIM (TIMEM3)                                                H39090
      TIMOT = TIMOT+TIMEM3-TIMEM2                                         H39100
C                                                                         H39110
C     NPL IS NOW LOCATION OF FIRST ELEMENT TO BE USED FOR NEXT PASS       H39120
C                                                                         H39130
      IPL = -2                                                            H39140
      DO 100 NL = NPL, NPE                                                H39150
         IPL = IPL+1                                                      H39160
         TRALYR(IPL) = TRALYR(NL)                                         H39170
         RADLYR(IPL) = RADLYR(NL)                                         H39180
  100 CONTINUE                                                            H39190
C                                                                         H39200
      V1P = V1P+FLOAT(NPL+1)*DVP                                          H39210
      NPE = IPL                                                           H39220
C                                                                         H39230
      GO TO 60                                                            H39240
  110 CONTINUE                                                            H39250
C                                                                         H39260
      CALL CPUTIM (TIME1)                                                 H39270
      TIM = TIME1-TIME                                                    H39280
      WRITE (IPR,910) TIME1,TIM,TIMEM,TIMRD,TIMTB,TIMOT                   H39290
C                                                                         H39300
      RETURN                                                              H39310
C                                                                         H39320
  900 FORMAT ('0 THE TIME AT THE START OF FLUXDN IS ',F12.3)              H39330
  905 FORMAT ('0 INITIAL LAYER',I5,'  FINAL LAYER',I5,/,'0 FILE ',I5,     H39340
     *        ' MERGED WITH FILE ',I5,' ONTO FILE',I5,'  WITH XTYPE=',    H39350
     *        G15.5)                                                      H39360
  910 FORMAT ('0 THE TIME AT THE END OF FLUXDN IS ',F12.3/F12.3,          H39370
     *        ' SECS WERE REQUIRED FOR THIS MERGE  - FLXIN - ',F12.3,     H39380
     *        ' - READ - ',F12.3,' - EMBND - ',F12.3,' - EMOUT - ',       H39390
     *        F12.3)                                                      H39400
C                                                                         H39410
      END                                                                 H39420
C
C     ----------------------------------------------------------------
C
      SUBROUTINE GETEXT (IEXFIL,LYRNOW,IEMITT)                            H39430
C                                                                         H39440
      IMPLICIT DOUBLE PRECISION (V)                                     ! H39450
C                                                                         H39460
C     THIS SUBROUTINE HAS BEEN MODIFIED TO ALSO READ IN THE ASYMMETRY     H39470
C     FACTORS AND TO SEARCH FOR THE PROPER LAYER IF DESIRED.              H39480
C                                                                         H39490
C                                          JAN 1986 (A.E.R. INC.)         H39500
C                                                                         H39510
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H39520
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFLD,        H39530
     *              NLTEFL,LNFIL4,LNGTH4                                  H39540
C                                                                         H39550
C     ROUTINE TO BUFFER IN AEROSOL ABSORPTION AND SCATTERRING             H39560
C     FROM TAPE 20 INTO COMMON ABSORB SCATTR, AND ASYMT                   H39570
C                                                                         H39580
      COMMON /PNLHDR/ V1P,V2P,DVP,NLIM,LDUM                               H39590
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H39600
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       H39610
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H39620
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H39630
     *              EXTID(10)                                             H39640
C                                                                         H39650
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H39660
C                                                                         H39670
      COMMON /FILHDA/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       H39680
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   H39690
     *                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    H39700
      COMMON /ABSORA/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(2025)                H39710
      COMMON /SCATTA/ V1SC,V2SC,DVSC,NPTSC,SCTTR(2025)                    H39720
      COMMON /ASYMMA/ V1AS,V2AS,DVAS,NPTAS,ASYMT(2025)                    H39730
      DIMENSION APNLHD(2),AFILHD(2)                                       H39740
C                                                                         H39750
      CHARACTER*40 CEXT                                                   H39760
C                                                                         H39770
      EQUIVALENCE (APNLHD(1),V1P) , (AFILHD(1),XID(1))                    H39780
C                                                                         H39790
      IF (IEMITT.EQ.0) THEN                                               H39800
         CALL BUFIN (IEXFIL,IEOF,AFILHD(1),NFHDRF)                        H39810
C                                                                         H39820
C     MOVE YID INTO EXTID                                                 H39830
C                                                                         H39840
         WRITE (CEXT,'(5A8)') (YID(I),I=3,7)                              H39850
         READ (CEXT,'(10A4)') EXTID                                       H39860
      ENDIF                                                               H39870
C                                                                         H39880
      IF (LYRNOW.EQ.1.AND.IEMITT.EQ.1) THEN                               H39890
         REWIND IEXFIL                                                    H39900
         CALL BUFIN (IEXFIL,IEOF,AFILHD(1),NFHDRF)                        H39910
C                                                                         H39920
C     MOVE YID INTO EXTID                                                 H39930
C                                                                         H39940
         WRITE (CEXT,'(5A8)') (YID(I),I=3,7)                              H39950
         READ (CEXT,'(10A4)') EXTID                                       H39960
      ENDIF                                                               H39970
C                                                                         H39980
      LAYER = 0                                                           H39990
C                                                                         H40090
   10 DO 20 I = 1, 2025                                                   H40100
         ABSRB(I) = 0.                                                    H40110
         SCTTR(I) = 0.                                                    H40120
         ASYMT(I) = 0.                                                    H40130
   20 CONTINUE                                                            H40140
C                                                                         H40150
      CALL BUFIN (IEXFIL,IEOF,APNLHD(1),NPHDRF)                           H40160
C                                                                         H40170
      LAYER = LAYER+1                                                     H40180
      IF (IEOF.LE.0) STOP 'GETEXT; IEXFIL EMPTY'                          H40190
C                                                                         H40200
      CALL BUFIN (IEXFIL,IEOF,ABSRB(1),NLIM)                              H40210
      CALL BUFIN (IEXFIL,IEOF,SCTTR(1),NLIM)                              H40220
      CALL BUFIN (IEXFIL,IEOF,ASYMT(1),NLIM)                              H40230
C                                                                         H40240
      V1ABS = V1P                                                         H40260
      V1SC = V1P                                                          H40270
      V1AS = V1P                                                          H40280
      V2ABS = V2P                                                         H40290
      V2SC = V2P                                                          H40300
      V2AS = V2P                                                          H40310
      DVABS = DVP                                                         H40320
      DVSC = DVP                                                          H40330
      DVAS = DVP                                                          H40340
      NPTABS = NLIM                                                       H40350
      NPTSC = NLIM                                                        H40360
      NPTAS = NLIM                                                        H40370
C                                                                         H40380
      RETURN                                                              H40390
C                                                                         H40400
      END                                                                 H40410
C
C     ----------------------------------------------------------------
C
      SUBROUTINE ADARSL (NNPTS,IEXFIL,MFILE,IAFIL,IEMIT)                  H40420
C                                                                         H40430
      IMPLICIT DOUBLE PRECISION (V)                                     ! H40440
C                                                                         H40450
C     ROUTINE TO ADD ABSORPTION AND SCATTERING TO THE TRANSMITTANCE       H40460
C     VALUES AT EACH POINT. THE AEROSOL VALUES ARE STORED IN              H40470
C     COMMON ABSORB AND COMMON SCATTR.                                    H40480
C                                                                         H40490
      PARAMETER (MXFSC=200,MXLAY=MXFSC+3,MXZMD=200,MXPDIM=MXLAY+MXZMD,
     *           IM2=MXPDIM-2,MXMOL=35,MXTRAC=22)
C
      COMMON R1(2410)                                                     H40500
      COMMON /ABSPNL/ V1P,V2P,DVP,NLIM,NSHFT,NPTS                         H40510
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   H40520
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       H40530
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       H40540
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      H40550
     *              EXTID(10)                                             H40560
      COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),
     *                MSPNL1(MXLAY),                                      H40570
     *                MSLAY1,ISFILE,JSFILE,KSFILE,LSFILE,MSFILE,IEFILE,   H40580
     *                JEFILE,KEFILE                                       H40590
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         H40600
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFLD,IEXFLD,        H40610
     *              NLTEFL,LNFIL4,LNGTH4                                  H40620
C                                                                         H40630
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H40640
C                                                                         H40650
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       H40660
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   H40670
     *                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    H40680
C                                                                         H40690
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHD(1),V1P)                     H40700
      EQUIVALENCE (FSCDID(4),IAERSL)                                      H40710
C                                                                         H40720
      CHARACTER*40 CEXT,CYID                                              H40730
C                                                                         H40740
      DIMENSION PNLHD(4),XFILHD(2)                                        H40750
C                                                                         H40760
      XKT0 = 0.6951*296.                                                  H40770
      CALL GETEXT (IEXFIL,LAYRS,IEMIT)                                    H40780
      CALL BUFIN (MFILE,IEOF,XFILHD(1),NFHDRF)                            H40790
C                                                                         H40800
C     FOR AEROSOL RUNS, MOVE EXTID INTO YID                               H40810
C                                                                         H40820
      WRITE (CEXT,'(10A4)') EXTID                                         H40830
      WRITE (CYID,'(5A8)') (YID(I),I=3,7)                                 H40840
      CYID(19:40) = CEXT(19:40)                                           H40850
      READ (CYID,'(5A8)') (YID(I),I=3,7)                                  H40860
C                                                                         H40870
      CALL BUFOUT (IAFIL,XFILHD(1),NFHDRF)                                H40880
      NPANLS = 0                                                          H40890
      IF (NOPR.EQ.0) WRITE (IPR,900) XID,(YID(N),N=1,2)                   H40900
   10 CALL BUFIN (MFILE,IEOF,PNLHD(1),NPHDRF)                             H40910
      IF (IEOF.LE.0) GO TO 40                                             H40920
      CALL BUFIN (MFILE,IEOF,R1(1),NLIM)                                  H40930
      IF (NPANLS.EQ.0) VIDV = V1P-DVP                                     H40940
      VAER = VIDV                                                         H40950
      VITST = VIDV                                                        H40960
      NPTS = NNPTS                                                        H40970
      RDLAST = -1.                                                        H40980
      NLIM1 = 0                                                           H40990
      NLIM2 = 0                                                           H41000
      RDDUM = 0.                                                          H41010
      AF = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)                  H41020
      RDF = RADFNI(VI,DVP,XKT0,VITST,RDEL,RDDUM)                          H41030
      IAFRD = 0                                                           H41040
      IF (VITST.GT.VAER) IAFRD = 1                                        H41050
C                                                                         H41060
   20 NLIM1 = NLIM2+1                                                     H41070
C                                                                         H41080
      VI = V1P+FLOAT(NLIM1-1)*DVP                                         H41090
      IF (IAFRD.EQ.0) THEN                                                H41100
         RADFN0 = RADFNI(VI,DVP,XKT0,VIDV,RDEL,RDLAST)                    H41110
         RADFN0 = RADFN0-RDEL                                             H41120
         VAER = -VIDV                                                     H41130
         EXT = AERF(VI,DVP,VAER,ADEL,TAUSCT,TDEL,GFACT,GDEL)              H41140
         EXT = EXT-ADEL                                                   H41150
      ELSE                                                                H41160
         EXT = AERF(VI,DVP,VIDV,ADEL,TAUSCT,TDEL,GFACT,GDEL)              H41170
         EXT = EXT-ADEL                                                   H41180
         VITST = -VIDV                                                    H41190
         RADFN0 = RADFNI(VI,DVP,XKT0,VITST,RDEL,RDLAST)                   H41200
         RADFN0 = RADFN0-RDEL                                             H41210
      ENDIF                                                               H41220
C                                                                         H41230
      NLIM2 = (VIDV-V1P)/DVP+1.001                                        H41240
      NLIM2 = MIN(NLIM2,NLIM)                                             H41250
C                                                                         H41260
      DO 30 K = NLIM1, NLIM2                                              H41270
         EXT = EXT+ADEL                                                   H41280
         RADFN0 = RADFN0+RDEL                                             H41290
         R1(K) = R1(K)+EXT*RADFN0                                         H41300
   30 CONTINUE                                                            H41310
C                                                                         H41320
      IF (NLIM2.LT.NLIM) GO TO 20                                         H41330
C                                                                         H41340
      CALL ABSOUT (V1P,V2P,DVP,NLIM,1,IAFIL,NPTS,R1,NPANLS)               H41350
      GO TO 10                                                            H41360
   40 CONTINUE                                                            H41370
      IAERSL = 2                                                          H41380
      RETURN                                                              H41390
C                                                                         H41400
  900 FORMAT ('1',5X,' AEROSOLS'//1X,10A8,2X,2(1X,A8,1X)//,5X,            H41410
     *        'FILE 20 AEROSOL EXTINCTIONS ADDED TO FILE 12 SENT TO ',    H41420
     *        'FILE 14')                                                  H41430
C                                                                         H41440
      END                                                                 H41450
C
C     ----------------------------------------------------------------
C
      FUNCTION AERF (VI,DVI,VINXT,ADEL,TAUSCT,TDEL,GFACT,GDEL)            H41460
C                                                                         H41470
      IMPLICIT DOUBLE PRECISION (V)                                     ! H41480
C                                                                         H41490
C     THIS FUNCTION CORRELATES THE AEROSOL FREQ. WITH THE LBLRTM          H41500
C     FREQ.  AND ADDS THE ABSORPTION TO THE                               H41510
C     SCATTERING TO PRODUCE THE EXTINCTION                                H41520
C                                                                         H41530
C     THIS FUNCTION HAS BEEN MODIFIED TO RETURN THE SCATTERING            H41540
C     SEPARATELY, AND TO ALSO RETURN THE ASYMMETRY FACTOR.                H41550
C                                                                         H41560
C                                         JAN 1986 (A.E.R. INC.)          H41570
C                                                                         H41580
      COMMON /ABSORA/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(2025)                H41590
      COMMON /SCATTA/ V1SC,V2SC,DVSC,NPTSC,SCTTR(2025)                    H41600
      COMMON /ASYMMA/ V1AS,V2AS,DVAS,NPTAS,ASYMT(2025)                    H41610
C                                                                         H41620
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & H41630
C                                                                         H41640
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       H41650
     *                WK(60),PZL,PZU,TZL,TZU,WN2   ,DV ,V1 ,V2 ,TBOUND,   H41660
     *                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    H41670
C                                                                         H41680
      EQUIVALENCE (XFILHD(1),XID(1))                                      H41690
      DIMENSION XFILHD(2)                                                 H41700
C                                                                         H41710
      DATA FACTOR / 1.E-03 /                                              H41720
C                                                                         H41730
      VDIF = VI-V1ABS                                                     H41740
      IAER = VDIF/DVABS+1.00001                                           H41750
      VAER = V1ABS+DVABS*FLOAT(IAER-1)                                    H41760
      DERIVS = (SCTTR(IAER+1)-SCTTR(IAER))/DVABS                          H41770
      DERIVA = (ASYMT(IAER+1)-ASYMT(IAER))/DVABS                          H41780
      DERIV = DERIVS+(ABSRB(IAER+1)-ABSRB(IAER))/DVABS                    H41790
C                                                                         H41800
C     TAUSCT IS THE SCATTERING TERM                                       H41810
C                                                                         H41820
      TAUSCT = SCTTR(IAER)+DERIVS*(VI-VAER)                               H41830
C                                                                         H41840
C     GFACT IS THE ASYMMETRY FACTOR                                       H41850
C                                                                         H41860
      GFACT = ASYMT(IAER)+DERIVA*(VI-VAER)                                H41870
      AERF = SCTTR(IAER)+ABSRB(IAER)+DERIV*(VI-VAER)                      H41880
C                                                                         H41890
C     ADEL, TDEL & GDEL ARE THE CHANGE PER DVI                            H41900
C                                                                         H41910
      ADEL = DERIV*DVI                                                    H41920
      TDEL = DERIVS*DVI                                                   H41930
      GDEL = DERIVA*DVI                                                   H41940
C                                                                         H41950
      VINXT = VAER+DVABS                                                  H41960
      VITST = VINXT                                                       H41970
      IF (FACTOR*AERF.LT.ABS(DERIV*1.0E+6))                               H41980
     *   VITST = VI+ABS(FACTOR*AERF/DERIV)                                H41990
      IF (VITST.LT.VINXT) VINXT = VITST                                   H42000
C                                                                         H42010
      RETURN                                                              H42020
C                                                                         H42030
      END                                                                 H42040

