C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      SUBROUTINE PLTLBL (IENDPL)                                          M00010
C                                                                         M00020
      IMPLICIT DOUBLE PRECISION (V)                                     ! M00030
C                                                                         M00040
      COMMON XX(2450),YY(2450)                                            M00050
C                                                                         M00060
      DOUBLE PRECISION XID,SEC   ,HMOL  ,XALTZ,YID,PROGID               & M00070
C                                                                         M00080
      COMMON /HVERSN/  HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,
     *                HVROPR,HVRPST,HVRPLT,HVRTST,HVRUTL,HVRXMR
      COMMON /PLTHDR/ XID(10),SEC,P0,T0,HMOL(60),XALTZ(4),                M00090
     *                W(60),PZL,PZU,TZL,TZU,WBROAD,DVT,V1V,V2V,TBOUND,    M00100
     *                EMISIV,FSCDID(17),NMOL,NLAYER,YID1,YID(10),LSTWDF   M00110
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       M00120
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,     M00130
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      M00140
     *                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES                      M00150
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           M00160
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M00170
      COMMON /POINTS/ XXI,YI                                              M00180
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         M00190
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        M00200
     *              NLTEFL,LNFIL4,LNGTH4                                  M00210
      COMMON /FLFORM/ CFORM                                               M00220
      DIMENSION PROGID(3),PLTHDR(2),PNLHDR(2),DUM(2),JCFN(0:3)            M00230
C                                                                         M00240
      CHARACTER*8 FSTAT,NSTAT,OSTAT                                       M00245
      CHARACTER*8 COPT,CDIF,CRAT                                          M00250
      CHARACTER*8 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *            HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR
      CHARACTER*11 CFORM,BFORM,FFORM                                      M00253
      CHARACTER*25 TAPEJJ,TAPELL,TAPEMM,TAPE29,TAPEST(99),                M00256
     *             JFILEN,LFILEN,MFILEN,JOUTNM,CFILEN(3),BLNKNM           M00260
      CHARACTER CEX*2,CEXST*2,CTAPE*4,HOTHER*6,CPRGID*60                  M00270
C                                                                         M00273
      LOGICAL OP,EX,IFIRST                                                M00276
C                                                                         M00280
      EQUIVALENCE (PLTHDR(1),XID(1)) , (PNLHDR(1),V1P),                   M00290
     *            (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),                M00300
     *            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),                M00310
     *            (FSCDID(5),IEMIT) , (FSCDID(6),ISCAN),                  M00320
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 M00330
     *            (FSCDID(9),JRAD) , (FSCDID(10),ITEST),                  M00340
     *            (FSCDID(11),IMRG) , (FSCDID(12),SCNID),                 M00350
     *            (FSCDID(13),HWHM) , (FSCDID(14),IDABS),                 M00360
     *            (FSCDID(15),IATM) , (FSCDID(16),LAYR1),                 M00370
     *            (FSCDID(17),NLAYFS)                                     M00380
      EQUIVALENCE (CFILEN(1),JFILEN) , (CFILEN(2),LFILEN),                M00383
     *            (CFILEN(3),MFILEN)                                      M00386
C                                                                         M00390
      DATA JCFN/1,1,3,3/                                                  M00393
C                                                                         M00396
      DATA HOTHER / ' OTHER'/                                             M00400
      DATA CDIF / 'DIFRENCE'/,CRAT / ' RATIO  '/                          M00410
      DATA CEXST/'EX'/,CTAPE / 'TAPE'/                                    M00420
      DATA BFORM/'FORMATTED  '/                                           M00421
      DATA NSTAT/'NEW'/,OSTAT/'OLD'/                                      M00422
      DATA TAPEST/99*'                         '/                         M00423
      DATA BLNKNM/'                         '/                            M00424
C                                                                         M00430
C                  PROGRAM TO PLOT LBLRTM RESULTS.                        M00440
C                                                                         M00450
C     CARD 1: NAME, PHONE EXTENSION, ID                                   M00460
C       FORMAT (A30)                                                      M00470
C                                                                         M00480
C     *********************************************                       M00490
C     *** CARDS 2A AND 3A FOR PLOT              ***                       M00500
C     *** CARD 2B FOR FILE DIFFERENCE OR RATIO  ***                       M00510
C     *********************************************                       M00520
C                                                                         M00530
C     CARD 2A: V1,V2,XSIZE,DELV,NUMSBX,NOENDX,LFILE,LSKIPF,SCALE,         M00540
C              IOPT,I4P,IXDEC                                             M00550
C       FORMAT (4F10.4,4I5,F10.3,I2,I3,I5)                                M00560
C                                                                         M00570
C                                                                         M00580
C       V1 IS THE INITIAL WAVENUMBER OF THE PLOT,                         M00590
C       V2 IS THE FINAL WAVENUMBER OF THE PLOT                            M00600
C            (USE INTEGRAL VALUES OF CM-1)                                M00610
C          A NEGATIVE NUMBER FOR V1 WILL TERMINATE PLOTTING SEQUENCE.     M00620
C                                                                         M00630
C       XSIZE IS THE NUMBER OF INCHES OF THE X-AXIS.                      M00640
C                                                                         M00650
C       DELV IS THE NUMBER OF WAVENUMBERS (CM-1) PER MAJOR DIVISION.      M00660
C                                                                         M00670
C       NUMSBX IS THE NUMBER OF SUBDIVISIONS PER MAJOR DIVISION OF        M00680
C          X-AXIS, I.E., 1 GIVES NO SUBDIVISION TIC-MARCS,                M00690
C                        2 GIVES ONE SUBDIVISION MARK                     M00700
C                          (IN OTHER WORDS 2 COMPARTMENTS                 M00710
C                           PER MAJOR DIVISION),                          M00720
C                        ETC.                                             M00730
C                                                                         M00740
C          ** FOR IOPT = 1, NUMSBX BECOMES THE QUANTITY JCAL WHICH  **    M00750
C          ** CONTROLS THE USE OF SYMBOLS ON THE OVERLAYED PLOT:    **    M00760
C          **                                                       **    M00770
C          **        JCAL = 0 WILL PRODUCE A LINE PLOT WITHOUT      **    M00780
C          **                 SYMBOLS                               **    M00790
C          **        JCAL > 0 WILL PRODUCE A LINE PLOT WITH SYMBOLS **    M00800
C          **                 AT EVERY JCALTH POINT                 **    M00810
C          **        JCAL < 0 WILL PRODUCE A POINT PLOT WITH A      **    M00820
C          **                 SYMBOL AT EVERY JCALTH POINT          **    M00830
C                                                                         M00840
C       NOENDX CONTROLS THE NUMBERS AT EITHER END OF THE X-AXIS, I.E.,    M00850
C            1 SUPRESSES THE NUMBERS AT EITHER END OF THE AXIS,           M00860
C            2 SUPRESSES THE BEGINNING NUMBER,                            M00870
C            3 THE ENDING NUMBER, AND                                     M00880
C            0 LEAVES BOTH.                                               M00890
C                                                                         M00900
C          ** FOR IOPT=1, NOENDX BECOMES THE QUANTITY LCAL WHICH    **    M00910
C          ** IS A NUMBER DESCRIBING THE SYMBOL TO BE USED ON THE   **    M00920
C          ** OVERLAYED PLOT.  FOR A LIST OF SYMBOLS SEE THE AFGL   **    M00930
C          ** CYBER PLOTTER HANDOUT, APPENDIX B.                    **    M00940
C                                                                         M00950
C       LFILE IS THE TAPE NUMBER OF FILE TO BE READ FROM LBLRTM.          M00960
C                                                                         M00970
C       LSKIPF IS THE NUMBER OF FILES TO BE SKIPPED IN TAPE LFILE         M00980
C            (NUMBER OF FILE TO BE PLOTTED WILL BE LSKIPF + 1).           M00990
C                                                                         M01000
C       SCALE ENABLES ONE TO ENLARGE OR REDUCE A PLOT.                    M01010
C                                                                         M01020
C       IOPT = 0 FOR PLOT,            = 1 FOR PLOT OVERLAY,               M01030
C            = 2 FOR FILE DIFFERENCE, = 3 FOR FILE RATIO.                 M01040
C                                                                         M01050
C       I4P = 0 FOR LINEAR CONNECTION OF POINTS,                          M01060
C           = 1 FOR FOUR POINT INTERP.                                    M01070
C                                                                         M01080
C       IXDEC IS NUMBER OF FIGURES AFTER DECIMAL POINT ON X-AXIS.         M01090
C                                                                         M01100
C                                                                         M01110
C     *** CARD 3A   (REQUIRED FOR CARD 2A, OTHERWISE OMIT) ***            M01120
C                                                                         M01130
C     CARD 3A: YMIN,YMAX,YSIZE,DELY,NUMSBY,NOENDY,IDEC,JEMIT,JPLOT,       M01140
C                                   LOGPLT,JHDR,JDUMMY,JOUT,JPLTFL        M01150
C         FORMAT (2G10.4,2F10.3,6I5,2(I2,I3))                             M01160
C                                                                         M01170
C                                                                         M01180
C       YMIN IS ORIGIN VALUE OF Y-AXIS,                                   M01190
C       YMAX IS VALUE AT TOP OF Y-AXIS                                    M01200
C            (THESE VALUES SHOULD BE EXPONENTS WHEN LOG PLOT IS           M01210
C             DESIRED, WHERE THE DIFFERENCE YMAX-YMIN WILL COINCIDE       M01220
C             TO THE NUMBER OF DECADES SPANNED).                          M01230
C                                                                         M01240
C          ** FOR IOPT = 1, YMIN DETERMINES THE Y-AXIS OFFSET OF   **     M01250
C          ** THE OVERLAYED PLOT WITH RESPECT TO THE ORIGINAL PLOT **     M01260
C                                                                         M01270
C       YSIZE IS THE NUMBER OF INCHES OF THE Y-AXIS.                      M01280
C                                                                         M01290
C       DELY IS THE NUMBER OF UNITS OF Y VALUES PER MAJOR DIVISION        M01300
C            (= 1. WHEN LOG PLOT IS DESIRED).                             M01310
C                                                                         M01320
C       NUMSBY IS THE NUMBER OF SUBDIVISIONS PER MAJOR DIVISION           M01330
C            OF Y-AXIS (SAME DEFINITION AS FOR NUMSBX).                   M01340
C                                                                         M01350
C       NOENDY CONTROLS THE NUMBERS AT EITHER END OF Y-AXIS               M01360
C            (SAME DEFINITION AS FOR NOENDX).                             M01370
C                                                                         M01380
C       IDEC IS NUMBER OF FIGURES AFTER DECIMAL POINT ON LINEAR Y-AXIS.   M01390
C                                                                         M01400
C       JEMIT = 0 FOR ABSORPTION, = 1 FOR RADIANCE.                       M01410
C                                                                         M01420
C       JPLOT = 0 FOR TRANSMISSION, = 1 FOR OPTICAL DEPTH (JEMIT=0)       M01430
C       JPLOT = 0 FOR WATTS, = 1 FOR TEMPERATURE PLOT (JEMIT=1).          M01440
C       JPLOT = 2 ,JEMIT = 0 DECIBELS SCALE                               M01450
C                                                                         M01460
C       LOGPLT = 0 FOR LINEAR Y-AXIS, = 1 FOR LOG SCALE.                  M01470
C                                                                         M01480
C       JHDR = 0 FOR HEADER PANEL TO BE PLOTTED                           M01490
C       JHDR = 1 FOR HEADER PANEL NOT PLOTTED                             M01500
C                                                                         M01510
C       JDUMMY - NOT CURRENTLY USED                                       M01520
C                                                                         M01530
C       JOUT = 0 FOR PLOT TO SCREEN                                       M01540
C       JOUT = 1 FOR PLOT SAVED TO FILE                                   M01550
C       JOUT = 2 FOR PLOT TO SCREEN AND SAVED TO FILE                     M01560
C       JOUT = 3 FOR PLOT PRINTED (ASCII) TO FILE                         M01570
C       JOUT = 4 FOR PLOT TO SCREEN AND PRINTED (ASCII) TO FILE           M01580
C                                                                         M01590
C       JPLTFL IS THE FILE TO WHICH THE PLOT IS SAVED OR PRINTED          M01600
C          IF JOUT .GE. 1      (DEFAULT = TAPE29)                         M01610
C                                                                         M01620
C                                                                         M01630
C      ** REPEAT CARD 2A & 3A, OR CARD 2B **                              M01640
C                                                                         M01650
C                                                                         M01660
C     CARD 2B: V1,V2,JFILE,JSKIPF,LFILE,LSKIPF,IOPT,MFILE                 M01670
C       FORMAT (2F10.4,20X,4I5,10X,I2,3X,I5)                              M01680
C                                                                         M01690
C       SUBROUTINE FILOPT WILL DIFFERENCE OR RATIO                        M01700
C         TWO LBLRTM OUTPUT FILES.                                        M01710
C                                                                         M01720
C              FILES MUST BE 1) UNFORMATTED LBLRTM FILES                  M01730
C                                                                         M01740
C                            2) SINGLE QUANTITY FILES                     M01750
C                                  I.E. CONTAIN ONLY ONE                  M01760
C                                       OF THE FOLLOWING:                 M01770
C                                                                         M01780
C                                        A) OPTICAL DEPTHS                M01790
C                                        B) TRANSMITTANCE                 M01800
C                                        C) RADIANCE                      M01810
C                                        D) TEMPERATURE                   M01820
C                                                                         M01830
C                  *** NOTE: STANDARD LBLRTM OUTPUT FILES USUALLY         M01840
C                            CONTAIN BOTH TRANSMITTANCE AND RADIANCE.     M01850
C                            THE USER MUST EITHER SCAN OR PLOT THE        M01860
C                            DESIRED QUANTITY TO A FILE FOR USE AS        M01870
C                            INPUT TO THIS ROUTINE.                       M01880
C                                                                         M01890
C       V1 IS THE INITIAL WAVENUMBER OF THE DIFFERENCE OR RATIO           M01900
C       V2 IS THE FINAL WAVENUMBER OF THE DIFFERENCE OR RATIO             M01910
C                                                                         M01920
C            A NEGATIVE NUMBER FOR V1 WILL TERMINATE PLOTTING SEQUENCE.   M01930
C                                                                         M01940
C       JFILE IS THE TAPE NUMBER OF FILE TO BE READ FROM LBLRTM.          M01950
C                    (NO DEFAULT)                                         M01960
C                                                                         M01970
C       JSKIPF IS THE NUMBER OF FILES TO BE SKIPPED IN TAPE JFILE         M01980
C            (NUMBER OF FILE TO BE USED WILL BE JSKIPF + 1).              M01990
C                                                                         M02000
C       LFILE IS THE TAPE NUMBER OF FILE TO BE READ FROM LBLRTM.          M02010
C                                                                         M02020
C       LSKIPF IS THE NUMBER OF FILES TO BE SKIPPED IN TAPE LFILE         M02030
C            (NUMBER OF FILE TO BE USED WILL BE LSKIPF + 1).              M02040
C                                                                         M02050
C       IOPT = 0 FOR PLOT,            = 1 FOR PLOT OVERLAY,               M02060
C            = 2 FOR FILE DIFFERENCE, = 3 FOR FILE RATIO.                 M02070
C                                                                         M02080
C          ***  DIFFERENCE WRITTEN ON MFILE IS (JFILE - LFILE)  ***       M02090
C          ***     RATIO WRITTEN ON MFILE IS (JFILE/LFILE)      ***       M02100
C                                                                         M02110
C       MFILE IS THE TAPE NUMBER OF FILE FOR DIFFERENCE/RATIO OUTPUT.     M02120
C                                                                         M02130
C       ** REPEAT CARDS 2A & 3A, OR CARD 2B **                            M02140
C                                                                         M02150
C**********************************************************************   M02160
C                                                                         M02170
C
C     ASSIGN SCCS VERSION NUMBER TO MODULE 
C
      HVRPLT = '$Revision$' 
C
      IENDPL = 0                                                          M02180
      JHDR = 0                                                            M02200
      JOUT = 0                                                            M02210
      JPLTFL = 0                                                          M02220
      JDUMMY = 0                                                          M02230
      YMIN = 0.                                                           M02240
      YMAX = 0.                                                           M02250
      YSIZ = 0.                                                           M02260
      DELY = 0.                                                           M02270
      NUMSBY = 0.                                                         M02280
      NOENDY = 0.                                                         M02290
      IDEC = 0                                                            M02300
      JEMIT = 0                                                           M02310
      JPLOT = 0                                                           M02320
      LOGPLT = 0                                                          M02330
      NPLT = 0                                                            M02340
      XLOGCN = -ALOG10(EXP(1.))                                           M02350
      X3 = 0.                                                             M02360
      YPL = 0.                                                            M02370
      READ (IRD,900) CPRGID,CEX                                           M02380
      IEXTRC=0                                                            M02383
      IF (CEX.EQ.CEXST) IEXTRC = 1                                        M02386
C                                                                         M02387
   10 WRITE (IPR,902) CPRGID                                              M02390
      IFIRST=.TRUE.                                                       M02391
      JFILEN = BLNKNM                                                     M02392
      LFILEN = BLNKNM                                                     M02393
      MFILEN = BLNKNM                                                     M02394
C                                                                         M02395
      READ (IRD,905) V1,V2,XSIZ,DELV,NUMSBX,NOENDX,LFILE,LSKIPF,          M02400
     *               SCALE,IOPT,I4P,IXDEC                                 M02410
C
      IF (V1.LT.0) GO TO 240                                              M02490
C
      IF (IEXTRC.EQ.1) THEN                                               M02412
         READ (IRD,906) (CFILEN(J),J=1,JCFN(IOPT))                        M02413
      ENDIF                                                               M02414
C                                                                         M02415
      WRITE (TAPELL,930) CTAPE,LFILE                                      M02416
      IF (IOPT.LE.1) LFILEN = JFILEN                                      M02417
      IF (LFILEN.NE.BLNKNM) THEN
         CALL CLJUST(LFILEN,25)                                           M02418
         TAPELL = LFILEN                                                  M02419
      ENDIF
      INQUIRE (UNIT=LFILE,OPENED=OP)                                      M02420
      IF (OP .AND. TAPEST(LFILE).NE.TAPELL) THEN                          M02421
         CLOSE (LFILE)                                                    M02422
         OP=.FALSE.                                                       M02423
      ENDIF                                                               M02424
      IF (.NOT.OP) THEN                                                   M02425
         INQUIRE (FILE=TAPELL,EXIST=EX)                                   M02426
         FSTAT=NSTAT                                                      M02427
         IF (EX) FSTAT=OSTAT                                              M02428
         OPEN (LFILE,FILE=TAPELL,STATUS=FSTAT,FORM=CFORM)                 M02429
         TAPEST(LFILE) = TAPELL                                           M02430
      ENDIF                                                               M02431
C                                                                         M02432
C     FOR IOPT = 0 - PLOT                                                 M02433
C     FOR IOPT = 1 - OVERLAY PLOT                                         M02440
C     FOR IOPT = 2 - DIFFERENCE VALUES                                    M02450
C     FOR IOPT = 3 - RATIO VALUES                                         M02460
C                                                                         M02470
      IF (IOPT.GT.1) THEN                                                 M02480
Csac     IF (V1.LT.0) GO TO 240                                           M02490
         JFILE = NUMSBX                                                   M02500
         JSKIPF = NOENDX                                                  M02510
         MFILE = IXDEC                                                    M02520
         IRDOPT = 0                                                       M02521
         CALL CLJUST(JFILEN,25)                                           M02522
         CALL CLJUST(LFILEN,25)                                           M02523
         CALL CLJUST(MFILEN,25)                                           M02524
         IF (JFILEN.EQ.BLNKNM) IRDOPT = IRDOPT+1                          M02525
         IF (LFILEN.EQ.BLNKNM) IRDOPT = IRDOPT+1                          M02526
         IF (IRDOPT.EQ.1) THEN                                            M02527
            WRITE(IPR,907) IRDOPT,JFILEN,LFILEN                           M02528
            STOP ' IRDOPT '                                               M02529
         ELSE                                                             M02530
            WRITE (TAPEJJ,930) CTAPE,JFILE                                M02531
            IF (JFILEN.NE.BLNKNM) TAPEJJ = JFILEN                         M02532
            INQUIRE (UNIT=JFILE,OPENED=OP)                                M02533
            IF (OP .AND. TAPEST(JFILE).NE.TAPEJJ) THEN                    M02534
               CLOSE (JFILE)                                              M02535
               OP=.FALSE.                                                 M02536
            ENDIF                                                         M02537
            IF (.NOT.OP) THEN                                             M02538
               INQUIRE (FILE=TAPEJJ,EXIST=EX)                             M02539
               FSTAT=NSTAT                                                M02540
               IF (EX) FSTAT=OSTAT                                        M02541
               OPEN (JFILE,FILE=TAPEJJ,STATUS=FSTAT,FORM=CFORM)           M02542
               TAPEST(JFILE) = TAPELL                                     M02543
            ENDIF                                                         M02544
            WRITE (TAPEMM,930) CTAPE,MFILE                                M02545
            IF (MFILEN.NE.BLNKNM) TAPEMM = MFILEN                         M02546
            INQUIRE (UNIT=MFILE,OPENED=OP)                                M02547
            IF (OP .AND. TAPEST(MFILE).NE.TAPEMM) THEN                    M02548
               CLOSE (MFILE)                                              M02549
               OP=.FALSE.                                                 M02550
            ENDIF                                                         M02551
            IF (.NOT.OP) THEN                                             M02552
               INQUIRE (FILE=TAPEMM,EXIST=EX)                             M02553
               FSTAT=NSTAT                                                M02554
               IF (EX) FSTAT=OSTAT                                        M02555
               OPEN (MFILE,FILE=TAPEMM,STATUS=FSTAT,FORM=CFORM)           M02556
               TAPEST(MFILE) = TAPEMM                                     M02557
            ENDIF                                                         M02558
         ENDIF                                                            M02559
         CALL FILOPT (V1,V2,JFILE,JSKIPF,LFILE,LSKIPF,MFILE,LENGTH,       M02560
     *                IOPT)                                               M02561
         GO TO 10                                                         M02562
      ENDIF                                                               M02563
      IF ((V1.LT.0.).AND.(DELV.LT.0.)) THEN                               M02570
         DELV = -DELV                                                     M02580
      ELSE                                                                M02590
         IF (V1.LT.0.) GO TO 240                                          M02600
      ENDIF                                                               M02610
C
      JOUTNM = BLNKNM
      READ (IRD,910) YMINR,YMAX,YSIZ,DELY,NUMSBY,NOENDY,IDEC,JEMIT,       M02620
     *                  JPLOT,LOGPLT,JHDR,JDUMMY,JOUT,JPLTFL              M02630
      IF (IEXTRC.EQ.1) THEN                                               M02615
         READ (IRD,911) JOUTM                                             M02636
      ENDIF                                                               M02638
      IF (JOUT.LT.0.OR.JOUT.GT.4) THEN                                    M02640
         WRITE (IPR,915) JOUT                                             M02650
         JOUT = 0                                                         M02660
      ENDIF                                                               M02670
      IFLAGJ = 1                                                          M02680
      NOAXES = 0                                                          M02690
      IF (IOPT.EQ.1) THEN                                                 M02700
         IF (NPLT.EQ.0) THEN                                              M02710
            WRITE (IPR,920)                                               M02720
            GO TO 10                                                      M02730
         ENDIF                                                            M02740
         NOAXES = 1                                                       M02750
         IFLAGJ = 2                                                       M02760
         JHDR = 1                                                         M02770
         V1 = V1STOR                                                      M02780
         V2 = V2STOR                                                      M02790
         XSIZ = XSIZES                                                    M02800
         SCALE = SCALES                                                   M02810
         YMIN = YMINST-YMINR                                              M02820
         YMAX = YMAXST-YMINR                                              M02830
         WRITE (IPR,927) YMINR,NUMSBX,NOENDX                              M02840
      ELSE                                                                M02850
         V1STOR = V1                                                      M02860
         V2STOR = V2                                                      M02870
         XSIZES = XSIZ                                                    M02880
         SCALES = SCALE                                                   M02890
         YMIN = YMINR                                                     M02900
         YMINST = YMIN                                                    M02910
         YMAXST = YMAX                                                    M02920
      ENDIF                                                               M02930
      IF (JOUT.EQ.1.OR.JOUT.EQ.3) THEN                                    M02940
         IFLAGJ = 0                                                       M02950
      ELSE                                                                M02960
         NPLT = NPLT+1                                                    M02970
      ENDIF                                                               M02980
      IF (NPLT.GT.1.AND.IFLAGJ.EQ.1) CALL PLOT (XSIZE+X3,-YPL,-3)         M02990
      IF (NPLT.GE.1.AND.IFLAGJ.GE.1) IENDPL = 1                           M03000
      IF (SCALE.LT.0.01) SCALE = 1.                                       M03010
      IF (JOUT.GE.1) THEN                                                 M03020
         IF (JPLTFL.EQ.0) JPLTFL = 29                                     M03030
         WRITE (TAPE29,930) CTAPE,JPLTFL                                  M03040
      IF (JOUTNM.NE.BLNKNM) THEN
         CALL CLJUST(JOUTNM,25)                                           M02418
         TAPE29 = JOUTNM                                                  M02419
      ENDIF
c**      CALL CLJUST(JOUTNM,25)                                           M03045
c**      IF (JOUTNM.NE.BLNKNM) TAPE29 = JOUTNM                            M03050
         INQUIRE (UNIT=JPLTFL,OPENED=OP)                                  M03060
         IF (OP .AND. TAPEST(JPLTFL).NE.TAPE29) THEN                      M03070
            CLOSE (JPLTFL)                                                M03071
            OP=.FALSE.                                                    M03072
         ENDIF                                                            M03073
         IF (.NOT.OP) THEN                                                M03074
            INQUIRE (FILE=TAPE29,EXIST=EX)                                M03075
            FSTAT=NSTAT                                                   M03076
            IF (EX) FSTAT=OSTAT                                           M03080
            FFORM=CFORM                                                   M03090
            IF (JOUT.GT.2) FFORM=BFORM                                    M03100
            OPEN (JPLTFL,FILE=TAPE29,STATUS=FSTAT,FORM=FFORM)             M03103
            TAPEST(JPLTFL) = TAPE29                                       M03106
         ENDIF                                                            M03109
      ENDIF                                                               M03110
      WRITE (IPR,932)                                                     M03120
      WRITE (IPR,935) V1,V2,XSIZ,DELV,NUMSBX,NOENDX,LFILE,LSKIPF,SCALE,   M03130
     *                IOPT,I4P,IXDEC                                      M03140
      WRITE (IPR,937)                                                     M03150
      WRITE (IPR,940) YMIN,YMAX,YSIZ,DELY,NUMSBY,NOENDY,IDEC,JEMIT,       M03160
     *                JPLOT,LOGPLT,JHDR,JOUT,JPLTFL                       M03170
      XSIZE = XSIZ*SCALE                                                  M03180
      IF (JOUT.NE.1.AND.JOUT.NE.3.AND.NOAXES.EQ.0) THEN                   M03190
         READ (CPRGID,945) PROGID                                         M03200
         CALL PLTID3 (PROGID,XSIZE+20.,11.0,1.0)                          M03210
      ENDIF                                                               M03220
C                                                                         M03230
      REWIND LFILE                                                        M03240
      CALL SKIPFL (LSKIPF,LFILE,IEOF)                                     M03250
      IEOF = LSKIPF                                                       M03260
      CALL BUFIN (LFILE,LEOF,PLTHDR(1),NFHDRF)                            M03270
C                                                                         M03280
      ICNTNM = MOD(IXSCNT,10)                                             M03290
      IXSECT = IXSCNT/10                                                  M03300
C                                                                         M03310
      WRITE (COPT,950) XID(10)                                            M03320
      IFUNCT = 0                                                          M03330
      IF (COPT.EQ.CDIF) IFUNCT = 1                                        M03340
      IF (COPT.EQ.CRAT) IFUNCT = 2                                        M03350
      IF (ISCAN.GE.1000) THEN                                             M03360
         ISCANT = MOD(ISCAN,100)                                          M03370
         IF (ISCANT.LE.0.OR.SCNID.EQ.-99.) ISCAN = ISCAN-ISCANT           M03380
      ELSE                                                                M03390
         IF (ISCAN.LE.0.OR.SCNID.EQ.-99.) ISCAN = 0                       M03400
      ENDIF                                                               M03410
C                                                                         M03420
      YSIZE = YSIZ*SCALE                                                  M03430
      YPL = (10.0-YSIZE)/2.                                               M03440
      IF (YPL.GT.0.) GO TO 20                                             M03450
      WRITE (IPR,955)                                                     M03460
      YSIZE = 10.                                                         M03470
   20 HGT = 0.140*SCALE                                                   M03480
      IF (NOAXES.EQ.0) DX = (V2-V1)/XSIZE                                 M03490
      DIVLNX = DELV/DX                                                    M03500
      NUMDVX = (V2-V1)/DELV+0.01                                          M03510
      SFY = 1.                                                            M03520
      IF (JPLOT.EQ.0.AND.JEMIT.EQ.1.AND.IFUNCT.NE.2) THEN                 M03530
C                                                                         M03540
C     BBFN BLACK BODY (INPUT TEMPERATURE - PLOT RADIANCE)                 M03550
C                                                                         M03560
         IF (YMAX.GE.2.0) THEN                                            M03570
            CALL BBSCLE                                                   M03580
         ELSE                                                             M03590
            IF (LOGPLT.GT.0) DELY = 1.                                    M03600
            NUMDVY = (YMAX-YMIN)/DELY+0.01                                M03610
         ENDIF                                                            M03620
         IF (LOGPLT.EQ.0) THEN                                            M03630
            NMAX = ALOG10(ABS(YMAX))                                      M03640
            SFY = 10.**(-NMAX+1)                                          M03650
         ENDIF                                                            M03660
      ELSE                                                                M03670
         IF (LOGPLT.GT.0) DELY = 1.                                       M03680
         NUMDVY = (YMAX-YMIN)/DELY+0.01                                   M03690
      ENDIF                                                               M03700
      IF (NOAXES.EQ.0) DY = (YMAX-YMIN)/YSIZE                             M03710
      DIVLNY = DELY/DY                                                    M03720
      WRITE (IPR,937)                                                     M03730
      WRITE (IPR,940) YMIN,YMAX,YSIZ,DELY,NUMSBY,NOENDY,IDEC,JEMIT,       M03740
     *                JPLOT,LOGPLT,JHDR,JOUT,JPLTFL                       M03750
      IF (IXDEC.LT.1) IXDEC = -1                                          M03760
      IF (IDEC.LT.1) IDEC = -1                                            M03770
C                                                                         M03780
      IF (JEMIT.EQ.2.AND.IEMIT.EQ.1.AND.ISCAN.GE.1) THEN                  M03790
         JEMIT = 1                                                        M03800
         IEMIT = 2                                                        M03810
      ENDIF                                                               M03820
      IF (JOUT.GE.1) THEN                                                 M03830
         IF (JOUT.LE.2) THEN                                              M03840
            V1S = V1V                                                     M03850
            V2S = V2V                                                     M03860
            DVS = DVT                                                     M03870
            ISCANS = ISCAN                                                M03880
            IEMITS = IEMIT                                                M03890
            IF (JEMIT.EQ.1.AND.JPLOT.EQ.1.AND.IEMIT.GE.1) IEMIT = 2       M03900
            V1V = V1                                                      M03910
            V2V = V2                                                      M03920
            ISCAN = ISCAN+1000                                            M03930
            IF (I4P.EQ.1) DVT = DVT/4.                                    M03940
            CALL BUFOUT (JPLTFL,PLTHDR(1),NFHDRF)                         M03950
            V1V = V1S                                                     M03960
            V2V = V2S                                                     M03970
            DVT = DVS                                                     M03980
            ISCAN = ISCANS                                                M03990
            IEMIT = IEMITS                                                M04000
         ELSE                                                             M04010
            WRITE (JPLTFL,960) XID,(YID(M),M=1,2)                         M04020
            WRITE (JPLTFL,962) LAYR1,NLAYER                               M04030
            WRITE (JPLTFL,965) SEC,P0,T0,DVT,V1V,V2V                      M04040
            WRITE (JPLTFL,967) HOTHER,WBROAD,(HMOL(I),W(I),I=1,NMOL)      M04050
            WRITE (JPLTFL,970) IHIRAC,ILBLF4,ICNTNM,IXSECT,IAERSL,        M04060
     *                         IEMIT,ISCAN,IPLOT,IPATHL,JRAD,SCNID,       M04070
     *                         HWHM                                       M04080
            IF (JEMIT.EQ.0.AND.JPLOT.EQ.0) WRITE (JPLTFL,972)             M04090
            IF (JEMIT.EQ.0.AND.JPLOT.EQ.1) WRITE (JPLTFL,975)             M04100
            IF (JEMIT.EQ.1.AND.JPLOT.EQ.0) WRITE (JPLTFL,977)             M04110
            IF (JEMIT.EQ.1.AND.JPLOT.EQ.1) WRITE (JPLTFL,980)             M04120
         ENDIF                                                            M04130
      ENDIF                                                               M04140
      WRITE (IPR,982) XID,(YID(M),M=1,2)                                  M04150
      WRITE (IPR,962) LAYR1,NLAYER                                        M04160
      WRITE (IPR,965) SEC,P0,T0,DVT,V1V,V2V                               M04170
      WRITE (IPR,967) HOTHER,WBROAD,(HMOL(I),W(I),I=1,NMOL)               M04180
      WRITE (IPR,970) IHIRAC,ILBLF4,ICNTNM,IXSECT,IAERSL,IEMIT,ISCAN,     M04190
     *                IPLOT,IPATHL,JRAD,SCNID,HWHM                        M04200
      WRITE (IPR,982)                                                     M04210
      ISCANT = MOD(ISCAN,100)                                             M04220
      IF (ISCANT.EQ.0) GO TO 30                                           M04230
      JEMSCN = SCNID/100.                                                 M04240
      IF (JEMIT.EQ.JEMSCN) GO TO 30                                       M04250
      WRITE (IPR,985)                                                     M04260
      GO TO 10                                                            M04270
   30 CONTINUE                                                            M04280
      IF (JOUT.NE.1.AND.JOUT.NE.3.AND.NOAXES.EQ.0) THEN                   M04290
         IF (JHDR.EQ.0) THEN                                              M04300
            CALL HEADER                                                   M04310
         ELSE                                                             M04320
            CALL PLOT (1.0,YPL,-3)                                        M04330
         ENDIF                                                            M04340
      ENDIF                                                               M04350
C                                                                         M04360
      CALL AXES (IGO,IPR,IREJ,SFY)                                        M04370
      IF (IREJ.EQ.1) GO TO 230                                            M04380
C                                                                         M04390
      J = 1                                                               M04400
      IF (I4P.EQ.1) J = 2                                                 M04410
      TIMLOP = 0.                                                         M04420
      TIMLIN = 0.                                                         M04430
      ISCANT = MOD(ISCAN,100)                                             M04440
      IF (ISCANT.EQ.0) GO TO 40                                           M04450
      IF (IGO.EQ.7.OR.IGO.EQ.8) THEN                                      M04460
         WRITE (IPR,985)                                                  M04470
         GO TO 230                                                        M04480
      ENDIF                                                               M04490
      IF (IGO.GE.25.AND.IGO.LE.30) THEN                                   M04500
         WRITE (IPR,985)                                                  M04510
         GO TO 230                                                        M04520
      ENDIF                                                               M04530
      IF (IGO.NE.12.AND.MOD(IGO,3).EQ.3) THEN                             M04540
         WRITE (IPR,985)                                                  M04550
         GO TO 230                                                        M04560
      ENDIF                                                               M04570
C                                                                         M04580
   40 CALL BUFIN (LFILE,LEOF,PNLHDR(1),NPHDRF)                            M04590
      IF (LEOF.LE.0) GO TO 220                                            M04600
      IF (V2P.GE.V1) GO TO 50                                             M04610
      CALL BUFIN (LFILE,LEOF,DUM(1),2)                                    M04620
      IF (ISCAN.GE.1) GO TO 40                                            M04630
      IF (MOD(IGO+1,3).NE.0) GO TO 40                                     M04640
      CALL BUFIN (LFILE,LEOF,DUM(1),2)                                    M04650
      GO TO 40                                                            M04660
   50 NLO = J                                                             M04670
      NHI = NLIM+J-1                                                      M04680
      WRITE (IPR,987) V1P,V2P,DV,NLIM,NLO,NHI                             M04690
      IF (ISCAN.GE.1) GO TO 60                                            M04700
      IF (MOD(IGO+1,3).NE.0) GO TO 60                                     M04710
      IF (MOD(IGO+1,6).EQ.0) GO TO 60                                     M04720
      CALL BUFIN (LFILE,LEOF,DUM(1),2)                                    M04730
   60 CALL BUFIN (LFILE,LEOF,Y(NLO),NLIM)                                 M04740
      IF (ISCAN.GE.1) GO TO 70                                            M04750
      IF (MOD(IGO+1,3).NE.0) GO TO 70                                     M04760
      IF (MOD(IGO+1,6).NE.0) GO TO 70                                     M04770
      CALL BUFIN (LFILE,LEOF,DUM(1),2)                                    M04780
C                                                                         M04790
   70 JMIN = J                                                            M04800
      JMAX = J                                                            M04810
      NEWPTS = 0                                                          M04820
      NST = NLO                                                           M04830
      IF (V1P.LT.V1) NST = (V1-V1P)/DV+0.5+FLOAT(NLO)                     M04840
      NND = NHI                                                           M04850
      IF (V2P.GT.V2) NND = (V2-V1P)/DV+0.5+FLOAT(NLO)                     M04860
      CALL CPUTIM (TIME0)                                                 M04870
      DO 80 I = NST, NND                                                  M04880
         XXI = V1P+DV*FLOAT(I-NLO)                                        M04890
C                                                                         M04900
C     YI=Y(I)                                                             M04910
C   J RANGES FROM NLO TO NLO-1+NEWPTS (=NPTS).                            M04920
C                                                                         M04930
         XX(J) = XXI                                                      M04940
         J = J+1                                                          M04950
   80 CONTINUE                                                            M04960
      J = JMIN                                                            M04970
C                                                                         M04980
      ISCANT = MOD(ISCAN,100)                                             M04990
C                                                                         M05000
      IF (ISCANT.EQ.0)                                                    M05010
     *    GO TO (100,90,210,210,90,210,90,110,210,210,120,90,160,170,     M05020
     *           210,210,210,210,140,130,210,210,130,210,130,150,210,     M05030
     *           210,210,210,180,190,210,210,210,210) IGO                 M05040
      IF (ISCANT.GE.1)                                                    M05050
     *    GO TO (90,90,210,210,90,210,210,210,210,210,120,90,170,170,     M05060
     *           210,210,210,210,130,130,210,210,130,210,210,210,210,     M05070
     *           210,210,210,190,190,210,210,210,210) IGO                 M05080
   90 CALL LINT (NST,NND,J)                                               M05090
      GO TO 200                                                           M05100
  100 CALL EXPT (NST,NND,J)                                               M05110
      GO TO 200                                                           M05120
  110 CALL XNTLOG (NST,NND,J)                                             M05130
      GO TO 200                                                           M05140
  120 CALL TEMPFN (NST,NND,J)                                             M05150
      GO TO 200                                                           M05160
  130 CALL TENLOG (NST,NND,J)                                             M05170
      GO TO 200                                                           M05180
  140 CALL LINSTC (NST,NND,J,XLOGCN)                                      M05190
      GO TO 200                                                           M05200
  150 CALL XLOGLN (NST,NND,J)                                             M05210
      GO TO 200                                                           M05220
  160 CALL DBOD (NST,NND,J)                                               M05230
      GO TO 200                                                           M05240
  170 CALL DBTR (NST,NND,J)                                               M05250
      GO TO 200                                                           M05260
  180 CALL LDBOPD (NST,NND,J)                                             M05270
      GO TO 200                                                           M05280
  190 CALL LDBTR (NST,NND,J)                                              M05290
  200 NEWPTS = J-NLO                                                      M05300
      CALL MNMX (NST,NND,JMIN,JMAX)                                       M05310
      IF (NEWPTS.EQ.0) GO TO 210                                          M05320
      NPTS = NLO-1+NEWPTS                                                 M05330
      WRITE (IPR,990) NPTS,XX(JMIN),YY(JMIN),XX(JMAX),YY(JMAX)            M05340
      IF (IFIRST) THEN                                                    M05341
         YYMIN=YY(JMIN)                                                   M05342
         XXMIN=XX(JMIN)                                                   M05343
         YYMAX=YY(JMAX)                                                   M05344
         XXMAX=XX(JMAX)                                                   M05345
         IFIRST=.FALSE.                                                   M05346
      ELSE                                                                M05347
         IF (YYMIN.GT.YY(JMIN)) THEN                                      M05348
            YYMIN=YY(JMIN)                                                M05349
            XXMIN=XX(JMIN)                                                M05350
         ENDIF                                                            M05351
         IF (YYMAX.LT.YY(JMAX)) THEN                                      M05352
            YYMAX=YY(JMAX)                                                M05353
            XXMAX=XX(JMAX)                                                M05354
         ENDIF                                                            M05355
      ENDIF                                                               M05356
      CALL CPUTIM (TIME1)                                                 M05357
      TIMLOP = TIMLOP+TIME1-TIME0                                         M05360
      IF (I4P.EQ.0) CALL FSCLIN (NLO,NPTS,XX,YY)                          M05370
      IF (I4P.EQ.1) CALL FPLINE (NLO,NPTS)                                M05380
      CALL CPUTIM (TIME2)                                                 M05390
      TIMLIN = TIMLIN+TIME2-TIME1                                         M05400
      J = 2                                                               M05410
      IF (I4P.EQ.1) J = 4                                                 M05420
  210 IF (V2P.LT.V2) GO TO 40                                             M05430
C                                                                         M05440
  220 X3 = 3.                                                             M05450
      WRITE (IPR,992) XXMIN,YYMIN,XXMAX,YYMAX                             M05455
      WRITE (IPR,995) TIMLOP,TIMLIN                                       M05460
  230 CONTINUE                                                            M05470
      IF (JOUT.EQ.1.OR.JOUT.EQ.2) CALL ENDFIL (JPLTFL)                    M05480
      IF (JOUT.GE.1) CLOSE (JPLTFL)                                       M05490
      GO TO 10                                                            M05500
C                                                                         M05510
  240 RETURN                                                              M05520
C                                                                         M05530
  900 FORMAT (A60,18X,A2)                                                 M05540
  902 FORMAT ('1'/30X,A60)                                                M05550
  905 FORMAT (4F10.4,4I5,F10.3,I2,I3,I5)                                  M05560
  906 FORMAT (A25)                                                        M05560
  907 FORMAT (/,' ILLEGAL VALUE FOR IRDOPT =',I2,/,' INDICATES THAT ',    M05562
     *        'ONLY FILE OF THE TWO REQUIRED FILES WAS NAMED. ',/         M05564
     *        ' BOTH FILES MUST EITHER BE NAMED, OR BOTH NOT NAMED.'/,    M05566
     *        ' NAME OF JFILE =',A25,' NAME OF LFILE =',A25,/)            M05568
  910 FORMAT (2G10.4,2G10.3,6I5,2(I2,I3))                                 M05570
  911 FORMAT (A25)                                                        M05570
  915 FORMAT ('0',10X,'INVALID VALUE FOR JOUT = ',I5,', RESET TO ZERO')   M05580
  920 FORMAT (//,'  *** OVERLAY FEATURE REQUIRES A STANDARD PLOT *** ',   M05590
     *        /,'  *** ON WHICH THE LINE CAN BE OVERLAYED ||    *** ')    M05600
  927 FORMAT (//,'  --- OVERLAY PLOT WILL HAVE AN OFFSET OF ',G10.3,/,    M05610
     *        '      WITH JCAL = ',I5,' AND LCAL = ',I5,/)                M05620
  930 FORMAT (A4,I2.2,'                   ')                              M05630
  932 FORMAT ('0'/8X,'V1        V2     XSIZE      DELV  SBX  NDX LFIL',   M05640
     *        ' SKPF     SCALE IOPT  I4P XDEC')                           M05650
  935 FORMAT (4F10.4,4I5,F10.3,3I5)                                       M05660
  937 FORMAT ('0       YMIN       YMAX      YSIZE       DELY  SBY  ',     M05670
     *        'NDY IDEC  JEM JPLT  LOG JHDR JOUT PFIL')                   M05680
  940 FORMAT (1X,1P2G11.4,2G11.3,9I5)                                     M05690
  945 FORMAT (3A10)                                                       M05700
  950 FORMAT (A8)                                                         M05710
  955 FORMAT ('     YSIZE TOO LARGE, HAS BEEN RESET TO 10 INCHES')        M05720
  960 FORMAT ('1'/10X,10A8,2X,2(1X,A8,1X))                                M05730
  962 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)            M05740
  965 FORMAT ('0',10X,'SECANT =',F15.5,/,'0',10X,'PRESS(MB) =',F12.5,/,   M05750
     *        '0',10X,'TEMP =',F8.2,/,'0',10X,'DV =',F12.8,' CM-1',/,     M05760
     *        '0',10X,'V1 =',F12.6,' CM-1',/,'0',10X,'V2 =',F12.6,        M05770
     *        ' CM-1',/,'0',10X)                                          M05780
  967 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//(5X,A6,' = ',         M05790
     *        1PE10.3))                                                   M05800
  970 FORMAT ('0',/4X,'IHIRAC    ILBLF4    ICNTNM    IXSECT    IAERSL',   M05810
     *        '     IEMIT     ISCAN     IPLOT    IPATHL      JRAD    ',   M05820
     *        'JFNVAR      HWHM'/10I10,F10.0,F10.5)                       M05830
  972 FORMAT ('0',5X,'WAVENUMBER',7X,' TRANSMISSION',/)                   M05840
  975 FORMAT ('0',5X,'WAVENUMBER',7X,'OPTICAL DEPTH',/)                   M05850
  977 FORMAT ('0',5X,'WAVENUMBER',7X,'   RADIANCE  ',/)                   M05860
  980 FORMAT ('0',5X,'WAVENUMBER',7X,' TEMPERATURE ',/)                   M05870
  982 FORMAT ('0'/10X,10A8,2X,2(1X,A8,1X))                                M05880
  985 FORMAT ('     RESULT FROM SCANNING FUNCTION INCONSISTENT WITH',     M05890
     *        ' PLOT REQUEST')                                            M05900
  987 FORMAT (10X,' * PANEL *',F10.3,F10.3,F20.6,3I10)                    M05910
  990 FORMAT ('0',9X,'NO. OF PTS =',I5,',    FREQ AT MIN =',F10.3,        M05920
     *        '  MIN VALUE =',1PE11.3,',    FREQ AT MAX =',0PF10.3,       M05930
     *        '  MAX VALUE =',1PE11.3/)                                   M05940
  992 FORMAT ('0',6X,'**** PLOT TOTALS ****    FREQ AT MIN =',F10.3,      M05942
     *        '  MIN VALUE =',1PE11.3,',    FREQ AT MAX =',0PF10.3,       M05944
     *        '  MAX VALUE =',1PE11.3/)                                   M05946
  995 FORMAT (10X,'TIME FOR LOOP =',F8.3,' SECONDS, TIME FOR LINE =',     M05950
     *        F8.3)                                                       M05960
C                                                                         M05970
      END                                                                 M05980
      SUBROUTINE FILOPT (V1,V2,JFILE,JSKIPF,LFILE,LSKIPF,MFILE,LENGTH,    M06040
     *                   IOPT)                                            M06050
C                                                                         M06060
      IMPLICIT DOUBLE PRECISION (V)                                     ! M06070
C                                                                         M06080
C----------------------------------------------------------------------   M06090
C             R.D. WORSHAM                                                M06100
C             ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                M06110
C             JANUARY 1990                                                M06120
C----------------------------------------------------------------------   M06130
C                                                                         M06140
C     THIS SUBROUTINE WILL DIFFERENCE OR RATIO TWO                        M06150
C          LBLRTM OUTPUT OR PLOT FILES                                    M06160
C                                                                         M06170
C          FILES MUST BE    1) UNFORMATTED LBLRTM FILES                   M06180
C                                                                         M06190
C                           2) SINGLE QUANTITY FILES                      M06200
C                              I.E. CONTAIN ONLY ONE OF THE FOLLOWING:    M06210
C                                                                         M06220
C                                    A) OPTICAL DEPTHS                    M06230
C                                    B) TRANSMITTANCE                     M06240
C                                    C) RADIANCE                          M06250
C                                    D) TEMPERATURE                       M06260
C                                                                         M06270
C                           3) THE DIFFERENCE IN DV BETWEEN THE FILES     M06280
C                              MUST BE LESS THAN 1.E-8                    M06290
C                                                                         M06300
C               *** NOTE: STANDARD LBLRTM OUTPUT FILES USUALLY            M06310
C                         CONTAIN BOTH TRANSMITTANCE AND RADIANCE         M06320
C                         THE USER MUST CREATE EITHER SCANNED OR          M06330
C                         PLOT FILES FOR USE AS INPUT TO THIS ROUTINE     M06340
C                                                                         M06350
C----------------------------------------------------------------------   M06360
C                                                                         M06370
      COMMON XX(2450),YY(2450)                                            M06380
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         M06390
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        M06400
     *              NLTEFL,LNFIL4,LNGTH4                                  M06410
C                                                                         M06420
      DOUBLE PRECISION XID,SEC,HMOL,XALTZ,YID,HDATE,HTIME               & M06430
C                                                                         M06440
      COMMON /PLTHDR/ XID(10),SEC,P0  ,T0  ,HMOL(60),XALTZ(4),            M06450
     *                W(60),PZL,PZU,TZL,TZU,WBROAD,DVT,V1V,V2V,TBOUND,    M06460
     *                EMISIV,FSCDID(17),NMOL,NLAYER,YID1,YID(10),LSTWDF   M06470
      COMMON /JYCOM/ V1PJ,V2PJ,DVJ,NLIMJ,YJ(2502)                         M06480
      COMMON /LYCOM/ V1PL,V2PL,DVL,NLIML,YL(2502)                         M06490
      DIMENSION PLTHDR(2),PNLHDJ(2),PNLHDL(2),COPT(2)                     M06500
C                                                                         M06510
      EQUIVALENCE (PLTHDR(1),XID(1)) , (PNLHDJ(1),V1PJ),                  M06520
     *            (PNLHDL(1),V1PL) , (FSCDID(1),IHIRAC),                  M06530
     *            (FSCDID(2),ILBLF4) , (FSCDID(3),IXSCNT),                M06540
     *            (FSCDID(4),IAERSL) , (FSCDID(5),IEMIT),                 M06550
     *            (FSCDID(6),ISCAN) , (FSCDID(7),IPLOT),                  M06560
     *            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                  M06570
     *            (FSCDID(10),ITEST) , (FSCDID(11),IMRG),                 M06580
     *            (FSCDID(12),SCNID) , (FSCDID(13),HWHM),                 M06590
     *            (FSCDID(14),IDABS) , (FSCDID(15),IATM),                 M06600
     *            (FSCDID(16),LAYR1) , (FSCDID(17),NLAYFS),               M06610
     *            (YID(1),HDATE) , (YID(2),HTIME)                         M06620
C                                                                         M06630
      LOGICAL IFIRST                                                      M06640
      CHARACTER COPT*8,HOTHER*6                                           M06650
C                                                                         M06660
      DATA COPT / 'DIFRENCE',' RATIO  '/                                  M06670
      DATA HOTHER / ' OTHER'/                                             M06680
C                                                                         M06690
      REWIND JFILE                                                        M06700
      REWIND LFILE                                                        M06710
      CALL SKIPFL (JSKIPF,JFILE,IJEOF)                                    M06720
      CALL SKIPFL (LSKIPF,LFILE,ILEOF)                                    M06730
      IJEOF = JSKIPF                                                      M06740
      ILEOF = LSKIPF                                                      M06750
      WRITE (IPR,900) JFILE                                               M06760
      CALL BUFIN (JFILE,JEOF,PLTHDR(1),NFHDRF)                            M06770
C                                                                         M06780
      ICNTNM = MOD(IXSCNT,10)                                             M06790
      IXSECT = IXSCNT/10                                                  M06800
C                                                                         M06810
      WRITE (IPR,905) XID,(YID(M),M=1,2)                                  M06820
      WRITE (IPR,910) LAYR1,NLAYER                                        M06830
      WRITE (IPR,915) SEC,P0,T0,DVT,V1V,V2V                               M06840
      WRITE (IPR,920) HOTHER,WBROAD,(HMOL(I),W(I),I=1,NMOL)               M06850
      ISCAN = MAX(ISCAN,0)                                                M06860
      WRITE (IPR,925) IHIRAC,ILBLF4,ICNTNM,IXSECT,IAERSL,IEMIT,ISCAN,     M06870
     *                IPLOT,IPATHL,JRAD,SCNID,HWHM                        M06880
      WRITE (IPR,905)                                                     M06890
C                                                                         M06900
      WRITE (IPR,930) LFILE                                               M06910
      CALL BUFIN (LFILE,LEOF,PLTHDR(1),NFHDRF)                            M06920
C                                                                         M06930
      ICNTNM = MOD(IXSCNT,10)                                             M06940
      IXSECT = IXSCNT/10                                                  M06950
C                                                                         M06960
      READ (COPT(IOPT-1),935) XID(10)                                     M06970
      CALL BUFOUT (MFILE,PLTHDR(1),NFHDRF)                                M06980
      WRITE (IPR,905) XID,(YID(M),M=1,2)                                  M06990
      WRITE (IPR,910) LAYR1,NLAYER                                        M07000
      WRITE (IPR,915) SEC,P0,T0,DVT,V1V,V2V                               M07010
      WRITE (IPR,920) HOTHER,WBROAD,(HMOL(I),W(I),I=1,NMOL)               M07020
      ISCAN = MAX(ISCAN,0)                                                M07030
      WRITE (IPR,925) IHIRAC,ILBLF4,ICNTNM,IXSECT,IAERSL,IEMIT,ISCAN,     M07040
     *                IPLOT,IPATHL,JRAD,SCNID,HWHM                        M07050
      WRITE (IPR,905)                                                     M07060
      IF (IOPT.EQ.2) WRITE (IPR,940)                                      M07070
      IF (IOPT.EQ.3) WRITE (IPR,945)                                      M07080
      J = 0                                                               M07090
      L = 0                                                               M07100
      M = 0                                                               M07110
      LIMOUT = 2400                                                       M07120
      IFIRST = .TRUE.                                                     M07130
      NLOW = 1                                                            M07140
      NLIMJ = 0                                                           M07150
      NLIML = 0                                                           M07160
C                                                                         M07170
   10 IF (J.GE.NLIMJ) THEN                                                M07180
         CALL BUFIN (JFILE,JEOF,PNLHDJ(1),NPHDRF)                         M07190
         WRITE (IPR,950) V1PJ,V2PJ,DVJ,NLIMJ,JFILE,JEOF                   M07200
         IF (JEOF.LE.0.AND.M.NE.0) GO TO 50                               M07210
         IF (JEOF.LE.0) GO TO 60                                          M07220
         CALL BUFIN (JFILE,JEOF,YJ(NLOW),NLIMJ)                           M07230
         J = 0                                                            M07240
      ENDIF                                                               M07250
      IF (L.GE.NLIML) THEN                                                M07260
         CALL BUFIN (LFILE,LEOF,PNLHDL(1),NPHDRF)                         M07270
         WRITE (IPR,950) V1PL,V2PL,DVL,NLIML,LFILE,LEOF                   M07280
         IF (LEOF.LE.0.AND.M.NE.0) GO TO 50                               M07290
         IF (LEOF.LE.0) GO TO 60                                          M07300
         CALL BUFIN (LFILE,LEOF,YL(NLOW),NLIML)                           M07310
         L = 0                                                            M07320
      ENDIF                                                               M07330
C                                                                         M07340
C     CHECK DV TO INSURE SAME POINT SPACING                               M07350
C                                                                         M07360
      IF (IFIRST) THEN                                                    M07370
         V1TST = DVJ*1.E-3                                                M07380
         DVTST = V1TST/2400.                                              M07390
         DVDIF = DVJ-DVL                                                  M07400
         IDVTST = 0                                                       M07410
         IF (ABS(DVDIF).GT.DVTST) THEN                                    M07420
            IDVTST = 1                                                    M07430
            WRITE (IPR,955) DVJ,DVL,DVTST                                 M07440
         ELSE                                                             M07450
C                                                                         M07460
C     IF DV IS OKAY, CHECK STARTING POINTS FOR BEATING                    M07470
C                                                                         M07480
            NTEST = (V1PJ-V1PL)/DVJ                                       M07490
            V1NEW = V1PL+FLOAT(NTEST)*DVJ                                 M07500
            V1DIF = V1NEW-V1PJ                                            M07510
            IF (ABS(V1DIF).GT.V1TST) THEN                                 M07520
               IDVTST = 1                                                 M07530
               WRITE (IPR,960) V1PJ,V1PL,V1TST                            M07540
            ENDIF                                                         M07550
         ENDIF                                                            M07560
         IF (IDVTST.EQ.1) STOP ' FILOPT - PANELS DO NOT MATCH '           M07570
      ENDIF                                                               M07580
      IF (V1PJ.NE.V1PL) THEN                                              M07590
         IF (V1PJ.GT.V1PL) THEN                                           M07600
            VDIF = (V1PJ-V1PL)/DVJ+0.5                                    M07610
            L = VDIF                                                      M07620
            IF (L.GE.NLIML) THEN                                          M07630
               NLIML = 0                                                  M07640
               GO TO 10                                                   M07650
            ENDIF                                                         M07660
         ELSE                                                             M07670
            VDIF = (V1PL-V1PJ)/DVJ+0.5                                    M07680
            J = VDIF                                                      M07690
            IF (J.GE.NLIMJ) THEN                                          M07700
               NLIMJ = 0                                                  M07710
               GO TO 10                                                   M07720
            ENDIF                                                         M07730
         ENDIF                                                            M07740
      ENDIF                                                               M07750
      IF (V2PJ.GE.V1) GO TO 20                                            M07760
      NLIMJ = 0                                                           M07770
   20 IF (V2PL.GE.V1.AND.NLIMJ.NE.0) GO TO 30                             M07780
      NLIML = 0                                                           M07790
      GO TO 10                                                            M07800
C                                                                         M07810
   30 NST = 0                                                             M07820
      IF (V1PJ.LT.V1) NST = (V1-V1PJ)/DVJ+0.5+FLOAT(NLOW)                 M07830
      J = J+NST                                                           M07840
      L = L+NST                                                           M07850
   40 CONTINUE                                                            M07860
      J = J+1                                                             M07870
      IF (J.GT.NLIMJ) GO TO 10                                            M07880
      L = L+1                                                             M07890
      IF (L.GT.NLIML) THEN                                                M07900
         J = J-1                                                          M07910
         GO TO 10                                                         M07920
      ENDIF                                                               M07930
      M = M+1                                                             M07940
C                                                                         M07950
      XX(M) = V1PJ+DVJ*FLOAT(J-1)                                         M07960
C                                                                         M07970
C     FOR IOPT = 2 - DIFFERENCE VALUES                                    M07980
C     FOR IOPT = 3 - RATIO VALUES                                         M07990
C                                                                         M08000
      IF (IOPT.EQ.2) THEN                                                 M08010
         YY(M) = YJ(J)-YL(L)                                              M08020
      ELSE                                                                M08030
         IF (YL(L).NE.0) THEN                                             M08040
            YY(M) = YJ(J)/YL(L)                                           M08050
         ELSE                                                             M08060
            YY(M) = 0.0                                                   M08070
         ENDIF                                                            M08080
      ENDIF                                                               M08090
      IF (IFIRST) THEN                                                    M08100
         YMIN = YY(M)                                                     M08110
         YMAX = YY(M)                                                     M08120
         IFIRST = .FALSE.                                                 M08130
      ELSE                                                                M08140
         YMIN = MIN(YMIN,YY(M))                                           M08150
         YMAX = MAX(YMAX,YY(M))                                           M08160
      ENDIF                                                               M08170
      IF (XX(M).GE.V2) GO TO 50                                           M08180
      IF (J.EQ.NLIMJ) GO TO 50                                            M08190
      IF (M.LE.LIMOUT) GO TO 40                                           M08200
C                                                                         M08210
   50 NPTS = M                                                            M08220
      V1PJS = V1PJ                                                        M08230
      V2PJS = V2PJ                                                        M08240
      NLIMJS = NLIMJ                                                      M08250
      V1PJ = XX(1)                                                        M08260
      V2PJ = XX(NPTS)                                                     M08270
      NLIMJ = NPTS                                                        M08280
      CALL BUFOUT (MFILE,PNLHDJ(1),NPHDRF)                                M08290
      CALL BUFOUT (MFILE,YY(NLOW),NPTS)                                   M08300
      M = 0                                                               M08310
      IF (V2PJ.LT.V2) THEN                                                M08320
         V1PJ = V1PJS                                                     M08330
         V2PJ = V2PJS                                                     M08340
         NLIMJ = NLIMJS                                                   M08350
         GO TO 40                                                         M08360
      ENDIF                                                               M08370
C                                                                         M08380
   60 WRITE (IPR,965) YMIN,YMAX                                           M08390
      CALL ENDFIL (MFILE)                                                 M08400
C                                                                         M08410
      RETURN                                                              M08420
C                                                                         M08430
  900 FORMAT (/,'0',4X,'FILOPT ***** READING FROM FILE',I3,' *****')      M08440
  905 FORMAT ('0',/,10X,10A8,2X,2(1X,A8,1X))                              M08450
  910 FORMAT (//,' INITIAL LAYER = ',I5,'   FINAL LAYER =',I5)            M08460
  915 FORMAT ('0',10X,'SECANT =',F15.5/'0',10X,'PRESS(MB) =',F12.5,/,     M08470
     *        '0',10X,'TEMP =',F8.2/'0',10X,'DV =',F12.8,' CM-1',/,       M08480
     *        '0',10X,'V1 =',F12.6,' CM-1'/'0',10X,'V2 =',F12.6,          M08490
     *        ' CM-1',/,'0',10X)                                          M08500
  920 FORMAT ('0 COLUMN DENSITY (MOLECULES/CM**2)'//(5X,A6,' = ',         M08510
     *        1PE10.3))                                                   M08520
  925 FORMAT ('0',/4X,'IHIRAC    ILBLF4    ICNTNM    IXSECT    ',         M08530
     *        'IAERSL     IEMIT     ISCAN     IPLOT    IPATHL      ',     M08540
     *        'JRAD    JFNVAR      HWHM',/,10(I10),F10.0,F10.5)           M08550
  930 FORMAT (/,'1',4X,'FILOPT ***** READING FROM FILE',I3,' *****')      M08560
  935 FORMAT (A8)                                                         M08570
  940 FORMAT (2(/),5X,'-----------------------------------------',/,      M08580
     *        5X,'OPT = 1 * (FILE 1)-(FILE 2) ** DIFFERENCE',/,           M08590
     *        5X,'-----------------------------------------')             M08600
  945 FORMAT (2(/),5X,'------------------------------------',/,           M08610
     *        5X,'OPT = 2 * (FILE 1)/(FILE 2) ** RATIO',/,                M08620
     *        5X,'------------------------------------')                  M08630
  950 FORMAT (/,5X,'* PANEL *   V1P = ',F10.3,' V2P = ',F10.3,            M08640
     *        ' DVP = ',F13.6,' NLIM = ',I5,' IFILE = ',I3,               M08650
     *        ' IEOF = ',I4)                                              M08660
  955 FORMAT ('  PANELS DO NOT MATCH FOR INPUT FILES ',/,'  DVJ = ',      M08670
     *        F20.10,'  DVL = ',F20.10,' DVTST = ',F14.10)                M08680
  960 FORMAT ('  PANELS DO NOT MATCH FOR INPUT FILES ',/,' V1PJ = ',      M08690
     *        F20.10,' V1PL = ',F20.10,' V1TST = ',F14.10)                M08700
  965 FORMAT (/,6X,'MINIMUM Y VALUE = ',1PE13.6,5X,                       M08710
     *        'MAXIMUM Y VALUE = ',1PE13.6)                               M08720
C                                                                         M08730
      END                                                                 M08740
      SUBROUTINE HEADER                                                   M08750
C                                                                         M08760
      IMPLICIT DOUBLE PRECISION (V)                                     ! M08770
C                                                                         M08780
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       M08790
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,     M08800
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      M08810
     *                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES                      M08820
C                                                                         M08830
      DOUBLE PRECISION XID,SEC   ,HMOL  ,XALTZ,YID                      & M08840
C                                                                         M08850
      COMMON /PLTHDR/ XID(10),SEC,P0  ,T0  ,HMOL(60),XALTZ(4),            M08860
     *                W(60),PZL,PZU,TZL,TZU,WBROAD,DVT,V1V,V2V,TBOUND,    M08870
     *                EMISIV,FSCDID(17),NMOL,NLAYER,YI1,YID(10),LSTWDF    M08880
C                                                                         M08890
      CHARACTER HTEN*2,HSE*5,HPR*11,HTP*10,HDV*11,HV1*11,HV2*11           M08900
      CHARACTER HSC*13,HVA*11,HHW*11,HLR*11,H4C*10,HAMNTS*19,HOTHER*7     M08910
      CHARACTER*10 HSLIT(0:4)                                             M08920
      CHARACTER XID1*72,XID2*30                                           M08930
C                                                                         M08940
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),                M08950
     *            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),                M08960
     *            (FSCDID(5),IEMIT) , (FSCDID(6),ISCAN),                  M08970
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 M08980
     *            (FSCDID(9),JRAD) , (FSCDID(10),ITEST),                  M08990
     *            (FSCDID(11),IMRG) , (FSCDID(12),SCNID),                 M09000
     *            (FSCDID(13),HWHM) , (FSCDID(14),IDABS),                 M09010
     *            (FSCDID(15),IATM) , (FSCDID(16),LAYR1),                 M09020
     *            (FSCDID(17),NLAYFS) , (YID(1),HDATE),                   M09030
     *            (YID(2),HTIME) , (YI1,IMULT)                            M09040
C                                                                         M09050
      DATA HSLIT / ' RECTANGLE',' TRIANGLE ','  GAUSS   ','  SINCSQ  ',   M09060
     *             '   SINC   '/                                          M09070
      DATA XP1 / 0.1 /, XP15 / 0.15 /, X1 / 1.0 /, X2 / 2.0 /,            M09080
     *     X26 / 2.6 /, X27 / 2.7 /, X3 / 3.0 /, X6 / 6.0 /,              M09090
     *     X7 / 7.0 /, X8 / 8.0 /, X8P3 / 8.3 /, X11 / 11. /,             M09100
     *     YP3 / 0.3 /, YSHF / 0.45 /, Y10 / 10.0 /                       M09110
      DATA HTEN / '10'/,                     HSE / 'SEC ='/,              M09120
     *     HPR / 'PRESS (MB)='/,             HTP / 'TEMP (K) ='/,         M09130
     *     HDV / 'DV(CM-1) = '/,             HV1 / 'V1(CM-1) = '/,        M09140
     *     HV2 / 'V2(CM-1) = '/,             HSC / 'SCAN FUNCTION'/,      M09150
     *     HVA / ' (VARIABLE)'/,             HHW / 'HWHM(CM-1)='/,        M09160
     *     HLR / 'NLAYERS  = '/,             H4C / 'H4CXAEISRM'/,         M09170
     *     HAMNTS / 'AMOUNTS (MOL/CM**2)'/,  HOTHER / ' OTHER '/          M09180
C                                                                         M09190
      ICNTNM = MOD(IXSCNT,10)                                             M09200
      IXSECT = IXSCNT/10                                                  M09210
C                                                                         M09220
      WRITE (XID1,'( 9(A8   ))') (XID(I),I=1,9)                           M09230
      WRITE (XID2,'(3(A10  ))') XID(10),(YID(I),I=1,2)                    M09240
      YT = Y10-YSHF                                                       M09250
      CALL SYMBOL (0.0,Y10,XP15,XID1,0.0,72)                              M09260
      CALL SYMBOL (0.0,YT,XP15,XID2,0.0,30)                               M09270
      YT = YT-2.*YSHF                                                     M09280
      CALL SYMBOL (X1,YT,XP15,HSE,0.0,5)                                  M09290
      RSEC = SEC                                                          M09300
      CALL NUMBER (X2,YT,XP15,RSEC,0.0,1)                                 M09310
      YT = YT-YSHF                                                        M09320
      CALL SYMBOL (X1,YT,XP15,HPR,0.0,11)                                 M09330
      CALL NUMBER (X3,YT,XP15,P0,0.0,3)                                   M09340
      YT = YT-YSHF                                                        M09350
      CALL SYMBOL (X1,YT,XP15,HTP,0.0,10)                                 M09360
      CALL NUMBER (X27,YT,XP15,T0,0.0,3)                                  M09370
      YT = YT-YSHF                                                        M09380
      CALL SYMBOL (X1,YT,XP15,HDV,0.0,11)                                 M09390
      CALL NUMBER (X26,YT,XP15,DVT,0.0,5)                                 M09400
      YT = YT-YSHF                                                        M09410
      CALL SYMBOL (X1,YT,XP15,HV1,0.0,11)                                 M09420
      RV1V = V1V                                                          M09430
      CALL NUMBER (X26,YT,XP15,RV1V,0.0,2)                                M09440
      YT = YT-YSHF                                                        M09450
      CALL SYMBOL (X1,YT,XP15,HV2,0.0,11)                                 M09460
      RV2V = V2V                                                          M09470
      CALL NUMBER (X26,YT,XP15,RV2V,0.0,2)                                M09480
      ISCANT = MOD(ISCAN,100)                                             M09490
      IF (ISCANT.EQ.0) GO TO 20                                           M09500
      JFNVAR = SCNID                                                      M09510
      JEMIT = JFNVAR/100                                                  M09520
      JFN = (JFNVAR-100*JEMIT)/10                                         M09530
      JVAR = JFNVAR-100*JEMIT-10*JFN                                      M09540
      YT = YT-YSHF                                                        M09550
      CALL SYMBOL (X1,YT,XP15,HSC,0.0,13)                                 M09560
      IF (JVAR.EQ.0) GO TO 10                                             M09570
      YT = YT-YP3                                                         M09580
      CALL SYMBOL (X1,YT,XP15,HVA,0.0,11)                                 M09590
   10 YT = YT-YP3                                                         M09600
C                                                                         M09610
      CALL SYMBOL (X1,YT,XP15,HSLIT(JFN),0.0,10)                          M09620
      YT = YT-YP3                                                         M09630
      CALL SYMBOL (X1,YT,XP15,HHW,0.0,11)                                 M09640
C                                                                         M09650
      CALL NUMBER (X27,YT,XP15,HWHM,0.0,5)                                M09660
   20 YT = YT-YSHF                                                        M09670
      CALL SYMBOL (X1,YT,XP15,HLR,0.0,11)                                 M09680
      XLAYER = NLAYER                                                     M09690
      CALL NUMBER (X26,YT,XP15,XLAYER,0.0,-1)                             M09700
      YT = YT-YSHF                                                        M09710
      CALL SYMBOL (X1,YT,XP15,H4C,0.0,10)                                 M09720
      YT = YT-YP3                                                         M09730
      JRAD = JRAD+4                                                       M09740
      CONTRL = JRAD+10*(ISCANT+100*(IEMIT+10*(IAERSL+10*(IXSECT+10*       M09750
     *                 (ICNTNM+10*(ILBLF4+10*IHIRAC))))))                 M09760
      CONTRL = CONTRL*10+IMULT                                            M09770
      CALL NUMBER (X1,YT,XP15,CONTRL,0.0,-1)                              M09780
      YT = YT-YP3                                                         M09790
      IF (IATM.NE.0) CALL LBH1H2 (YT,X1,YID)                              M09800
      YT = YT-YP3                                                         M09810
      IF (IAERSL.NE.0) CALL LBAERS (YT,X1,YID)                            M09820
      YT = Y10-2.*YSHF                                                    M09830
      CALL SYMBOL (X6,YT,XP15,HAMNTS,0.0,19)                              M09840
      YT = YT-YSHF                                                        M09850
      YTE = YT+XP1                                                        M09860
      RW = 0.                                                             M09870
      HIP = 0.                                                            M09880
      IF (WBROAD.LT.1.E-23) GO TO 30                                      M09890
      IP = ALOG10(WBROAD)                                                 M09900
      HIP = IP                                                            M09910
      RW = WBROAD/10.**IP                                                 M09920
   30 CALL SYMBOL (X6,YT,XP15,HOTHER,0.0,7)                               M09930
      CALL NUMBER (X7,YT,XP15,RW,0.0,3)                                   M09940
      CALL SYMBOL (X8,YT,XP15,HTEN,0.0,2)                                 M09950
      CALL NUMBER (X8P3,YTE,XP1,HIP,0.0,-1)                               M09960
      DO 40 M = 1, NMOL                                                   M09970
         NM = M                                                           M09980
         IF (W(M).LE.0.) GO TO 40                                         M09990
         YT = YT-YSHF                                                     M10000
         IF (YT.LT.1.0) GO TO 50                                          M10010
         YTE = YT+XP1                                                     M10020
         RW = 0.                                                          M10030
         HIP = 0.                                                         M10040
         IP = ALOG10(W(M))                                                M10050
         HIP = IP                                                         M10060
         RW = W(M)/10.**IP                                                M10070
         CALL SYMBOL (X6,YT,XP15,HMOL(M),0.0,6)                           M10080
         CALL NUMBER (X7,YT,XP15,RW,0.0,3)                                M10090
         CALL SYMBOL (X8,YT,XP15,HTEN,0.0,2)                              M10100
         CALL NUMBER (X8P3,YTE,XP1,HIP,0.0,-1)                            M10110
   40 CONTINUE                                                            M10120
      CALL PLOT (X11,0.0,-3)                                              M10130
      CALL PLOT (1.0,YPL,-3)                                              M10140
C                                                                         M10150
      RETURN                                                              M10160
C                                                                         M10170
   50 YT = Y10-2.*YSHF                                                    M10180
      X12 = X6+X6                                                         M10190
      X13 = X12+X1                                                        M10200
      X14 = X13+X1                                                        M10210
      X14P3 = X14+0.3                                                     M10220
      CALL SYMBOL (X12,YT,XP15,HAMNTS,0.0,19)                             M10230
      DO 60 M = NM, NMOL                                                  M10240
         IF (W(M).LE.0.) GO TO 60                                         M10250
         YT = YT-YSHF                                                     M10260
         YTE = YT+XP1                                                     M10270
         RW = 0.                                                          M10280
         HIP = 0.                                                         M10290
         IP = ALOG10(W(M))                                                M10300
         HIP = IP                                                         M10310
         RW = W(M)/10.**IP                                                M10320
         CALL SYMBOL (X12,YT,XP15,HMOL(M),0.0,6)                          M10330
         CALL NUMBER (X13,YT,XP15,RW,0.0,3)                               M10340
         CALL SYMBOL (X14,YT,XP15,HTEN,0.0,2)                             M10350
         CALL NUMBER (X14P3,YTE,XP1,HIP,0.0,-1)                           M10360
   60 CONTINUE                                                            M10370
      XFN = X14+X3                                                        M10380
      CALL PLOT (XFN,0.0,-3)                                              M10390
      CALL PLOT (1.0,YPL,-3)                                              M10400
C                                                                         M10410
      RETURN                                                              M10420
C                                                                         M10430
      END                                                                 M10440
      SUBROUTINE AXES (IGO,LOUT,IREJ,SFY)                                 M10450
C                                                                         M10460
      IMPLICIT DOUBLE PRECISION (V)                                     ! M10470
C                                                                         M10480
      DOUBLE PRECISION XID,SEC,HMOL,XALTZ,YID,HDATE,HTIME               & M10490
C                                                                         M10500
      COMMON /PLTHDR/ XID(10),SEC,P0  ,T0  ,HMOL(60),XALTZ(4),            M10510
     *                W(60),PZL,PZU,TZL,TZU,WBROAD,DVT,V1V,V2V,TBOUND,    M10520
     *                EMISIV,FSCDID(17),NMOL,NLAYER,YID1,YID(10),LSTWDF   M10530
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       M10540
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,     M10550
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      M10560
     *                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES                      M10570
C                                                                         M10580
      CHARACTER YIDC*43,ASTR*3                                            M10590
C                                                                         M10600
      CHARACTER TITL1*12,TITL2*13,TITL3*50,TITL4*11,XWAVEN*17,XWAVEL*23   M10610
      CHARACTER BLK*2,HSCALE*10,GIGAHZ*15,TITL5*12,TITL6*9,TITL1D*25,     M10620
     *          TITL1R*20,TITL2D*26,TITL2R*21,TITL3D*34,TITL3R*29,        M10630
     *          TITL4D*24,TITL4R*19,TITL5D*23,TITL5R*18,TITL6D*21,        M10640
     *          TITL6R*16                                                 M10650
C                                                                         M10660
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),                M10670
     *            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),                M10680
     *            (FSCDID(5),IEMIT) , (FSCDID(6),ISCAN),                  M10690
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 M10700
     *            (FSCDID(9),JRAD) , (FSCDID(10),ITEST),                  M10710
     *            (FSCDID(11),IMRG) , (FSCDID(12),SCNID),                 M10720
     *            (FSCDID(13),HWHM) , (FSCDID(14),IDABS),                 M10730
     *            (FSCDID(15),IATM) , (FSCDID(16),LAYR1),                 M10740
     *            (FSCDID(17),NLAYFS)                                     M10750
C                                                                         M10760
      DATA NX / -17 /,LX / 23 /,NT1 / 12 /,NT2 / 13 /,NT3 / 37 /,         M10770
     *     NT4 / 11 /,NGH / 15 /,NT6 / 9 /,NT1D / 25 /,NT1R / 20 /,       M10780
     *     NT2D / 26 /,NT2R / 21 /,NT3D / 21 /,NT3R / 16 /,NT4D / 24 /,   M10790
     *     NT4R / 19 /,NT5D / 23 /,NT5R / 18 /,NT6D / 21 /,NT6R / 16 /    M10800
      DATA XWAVEN / 'WAVENUMBER (CM-1)'/,                                 M10810
     *     XWAVEL / 'WAVELENGTH (MICROMETER)'/,                           M10820
     *     GIGAHZ / 'FREQUENCY (GHZ)'/,                                   M10830
     *     BLK / ' '/,                                                    M10840
     *     TITL1 / 'TRANSMISSION'/,                                       M10850
     *     TITL2 / 'OPTICAL DEPTH'/,                                      M10860
     *     TITL3 / ' RADIANCE WATTS /  (CM**2*STER*CM-1)  * '/,           M10870
     *     TITL4 / 'TEMPERATURE'/,                                        M10880
     *     TITL5 / ' ABSORPTION '/,                                       M10890
     *     TITL6 / ' DECIBELS'/                                           M10900
      DATA TITL1D / ' TRANSMISSION DIFFERENCE '/,                         M10910
     *     TITL1R / ' TRANSMISSION RATIO '/,                              M10920
     *     TITL2D / ' OPTICAL DEPTH DIFFERENCE '/,                        M10930
     *     TITL2R / ' OPTICAL DEPTH RATIO '/,                             M10940
     *     TITL3D / ' RADIANCE DIFFERENCE  * '/,                          M10950
     *     TITL3R / ' RADIANCE RATIO  * '/,                               M10960
     *     TITL4D / ' TEMPERATURE DIFFERENCE '/,                          M10970
     *     TITL4R / ' TEMPERATURE RATIO '/,                               M10980
     *     TITL5D / ' ABSORPTION DIFFERENCE '/,                           M10990
     *     TITL5R / ' ABSORPTION RATIO '/,                                M11000
     *     TITL6D / ' DECIBELS DIFFERENCE '/,                             M11010
     *     TITL6R / ' DECIBELS RATIO '/                                   M11020
C                                                                         M11030
      DATA ASTR / '*  '/                                                  M11040
      DATA XP1 / 0.1 /,X11 / 11.0 /,YP5 / 0.5 /                           M11050
      DATA CL / 29.9792458 /                                              M11060
C                                                                         M11070
      IREJ = 0                                                            M11080
      IGO = IEMIT+3*JEMIT+6*JPLOT+18*LOGPLT+1                             M11090
      IF (JOUT.EQ.1.OR.JOUT.EQ.3.OR.NOAXES.EQ.1) THEN                     M11100
         GO TO (10,10,140,140,10,140,10,10,140,140,10,10,10,10,140,       M11110
     *          140,140,140,10,10,140,140,10,140,10,10,140,140,140,       M11120
     *          140,10,10,140,140,140,140) IGO                            M11130
   10    RETURN                                                           M11140
      ENDIF                                                               M11150
C                                                                         M11160
      RV1 = V1                                                            M11170
      CALL AXISL (0.0,0.0,XWAVEN,NX,NUMDVX,DIVLNX,NUMSBX,RV1,DELV,        M11180
     *            IXDEC,0.,HGT,1,0,NOENDX,0,0)                            M11190
      IF (V1.LT.100.) GO TO 40                                            M11200
      CALL AXISL (0.0,YSIZE,BLK,1,NUMDVX,DIVLNX,NUMSBX,RV1,DELV,          M11210
     *            IXDEC,0.,HGT,0,0,NOENDX,1,0)                            M11220
C                                                                         M11230
C   CREATE WAVELENGTH SCALE (MICRONS) OF APPROXIMATELY FIVE VALUES:       M11240
C                                                                         M11250
      IF (YSIZE.GT.8.) GO TO 60                                           M11260
      YSZSH = YSIZE+YP5                                                   M11270
      YINCH = YSZSH+HGT                                                   M11280
      CALL AXISL (0.0,YSZSH,XWAVEL,LX,1,XSIZE,1,RV1,DELV,-1,0.,HGT,0,0,   M11290
     *            1,0,0)                                                  M11300
      TOP = 10000./V1                                                     M11310
      STEP = (TOP-10000./V2)/5.                                           M11320
      NDP = -ALOG10(STEP)+1.                                              M11330
      ISTEP = STEP*(10.**NDP)                                             M11340
      IF (ISTEP.GE.6) ISTEP = 10+10*((ISTEP-5)/10)                        M11350
      STEP = FLOAT(ISTEP)/(10.**NDP)                                      M11360
      IF (NDP.LE.0) NDP = -1                                              M11370
      TOP = TOP-AMOD(TOP,STEP)                                            M11380
   20 BOT = 10000./TOP                                                    M11390
      IF (BOT.GT.V2) GO TO 60                                             M11400
      XINCH = (BOT-V1)/DX                                                 M11410
      IF (XINCH.GT.XSIZE) GO TO 30                                        M11420
      CALL SYMBOL (XINCH,YSZSH,XP1,3,0.0,-1)                              M11430
C                                                                         M11440
C     ('3' AS FOURTH ARGUMENT REPRESENTS PLUS SYMBOL FOR TIC MARK)        M11450
C                                                                         M11460
      DIGITS = 0.                                                         M11470
      IF (TOP.GE.10.) DIGITS = ALOG10(TOP)                                M11480
      DIGITS = AINT(DIGITS+1.0E-12)                                       M11490
      SIZNUM = HGT*(DIGITS+FLOAT(NDP)+1.7)                                M11500
      XPOS = XINCH-0.5*SIZNUM                                             M11510
      IF (XPOS.GT.XSIZE) GO TO 30                                         M11520
      CALL NUMBER (XPOS,YINCH,HGT,TOP,0.0,NDP)                            M11530
   30 TOP = TOP-STEP                                                      M11540
      GO TO 20                                                            M11550
C                                                                         M11560
C     CREATE FREQUENCY SCALE IN GIGAHERTZ                                 M11570
C                                                                         M11580
   40 RV1 = V1                                                            M11590
      CALL AXISL (0.0,YSIZE,GIGAHZ,NGH,1,XSIZE,1,RV1,DELV,-1,0.,HGT,0,    M11600
     *            0,1,0,0)                                                M11610
      IF (DELV.LT.1./CL) GO TO 60                                         M11620
      HLFHGT = 0.5*HGT                                                    M11630
      YPOS = YSIZE+HGT+HLFHGT                                             M11640
      DO 50 I = 1, NUMDVX                                                 M11650
         IG = (V1+FLOAT(I-1)*DELV)*CL                                     M11660
         G = IG                                                           M11670
         IF (IG.NE.0) G = IG+1                                            M11680
         XG = ((G/CL)-V1)/DX                                              M11690
         IF (XG.GT.XSIZE) GO TO 50                                        M11700
         CALL PLOT (XG,YSIZE,3)                                           M11710
         CALL PLOT (XG,YSIZE+HLFHGT,2)                                    M11720
         DIGITS = 0.                                                      M11730
         IF (G.GE.10.) DIGITS = ALOG10(G)                                 M11740
         DIGITS = AINT(DIGITS+1.0E-12)                                    M11750
         SIZNUM = HGT*(DIGITS+0.7)                                        M11760
         XPOS = XG-0.5*SIZNUM                                             M11770
         IF (XPOS.GT.XSIZE) GO TO 50                                      M11780
         CALL NUMBER (XPOS,YPOS,HGT,G,0.0,-1)                             M11790
   50 CONTINUE                                                            M11800
C                                                                         M11810
   60 GO TO (70,70,140,140,100,140,90,90,140,140,130,130,110,110,         M11820
     *       140,140,140,140,70,70,140,140,100,140,90,90,140,140,         M11830
     *       140,140,120,120,140,140,140,140) IGO                         M11840
C                                                                         M11850
   70 IF (IDABS.LT.0) GO TO 80                                            M11860
      IF (IGO.LT.13) THEN                                                 M11870
         IF (IFUNCT.EQ.0) CALL PLTLIN (TITL1,BLK,NT1,1.0)                 M11880
         IF (IFUNCT.EQ.1) CALL PLTLIN (TITL1D,BLK,NT1D,1.0)               M11890
         IF (IFUNCT.EQ.2) CALL PLTLIN (TITL1R,BLK,NT1R,1.0)               M11900
      ELSE                                                                M11910
         IF (IFUNCT.EQ.0) CALL PLTLOG (TITL1,BLK,NT1)                     M11920
         IF (IFUNCT.EQ.1) CALL PLTLOG (TITL1D,BLK,NT1D)                   M11930
         IF (IFUNCT.EQ.2) CALL PLTLOG (TITL1R,BLK,NT1R)                   M11940
      ENDIF                                                               M11950
      GO TO 150                                                           M11960
C                                                                         M11970
   80 IF (IGO.LT.13) THEN                                                 M11980
         IF (IFUNCT.EQ.0) CALL PLTLIN (TITL5,BLK,NT1,1.0)                 M11990
         IF (IFUNCT.EQ.1) CALL PLTLIN (TITL5D,BLK,NT1D,1.0)               M12000
         IF (IFUNCT.EQ.2) CALL PLTLIN (TITL5R,BLK,NT1R,1.0)               M12010
      ELSE                                                                M12020
         IF (IFUNCT.EQ.0) CALL PLTLOG (TITL5,BLK,NT1)                     M12030
         IF (IFUNCT.EQ.1) CALL PLTLOG (TITL5D,BLK,NT1D)                   M12040
         IF (IFUNCT.EQ.2) CALL PLTLOG (TITL5R,BLK,NT1R)                   M12050
      ENDIF                                                               M12060
      GO TO 150                                                           M12070
C                                                                         M12080
   90 IF (IGO.LT.13) THEN                                                 M12090
         IF (IFUNCT.EQ.0) CALL PLTLIN (TITL2,BLK,NT2,1.0)                 M12100
         IF (IFUNCT.EQ.1) CALL PLTLIN (TITL2D,BLK,NT2D,1.0)               M12110
         IF (IFUNCT.EQ.2) CALL PLTLIN (TITL2R,BLK,NT2R,1.0)               M12120
      ELSE                                                                M12130
         IF (IFUNCT.EQ.0) CALL PLTLOG (TITL2,BLK,NT2)                     M12140
         IF (IFUNCT.EQ.1) CALL PLTLOG (TITL2D,BLK,NT2D)                   M12150
         IF (IFUNCT.EQ.2) CALL PLTLOG (TITL2R,BLK,NT2R)                   M12160
      ENDIF                                                               M12170
      GO TO 150                                                           M12180
C                                                                         M12190
  100 IF (IGO.GE.13) THEN                                                 M12200
         IF (IFUNCT.EQ.0) CALL PLTLOG (TITL3,BLK,NT3)                     M12210
         IF (IFUNCT.EQ.1) CALL PLTLOG (TITL3D,BLK,NT3D)                   M12220
         IF (IFUNCT.EQ.2) CALL PLTLOG (TITL3R,BLK,NT3R)                   M12230
         GO TO 150                                                        M12240
      ENDIF                                                               M12250
      WRITE (HSCALE,'(1PE10.0)') SFY                                      M12260
C                                                                         M12270
      NT = 50                                                             M12280
      NTD = 34                                                            M12290
      NTR = 29                                                            M12300
      TITL3(41:50) = HSCALE                                               M12310
      TITL3D(25:34) = HSCALE                                              M12320
      TITL3R(20:29) = HSCALE                                              M12330
      IF (IFUNCT.EQ.0) CALL PLTLIN (TITL3,BLK,NT,SFY)                     M12340
      IF (IFUNCT.EQ.1) CALL PLTLIN (TITL3D,BLK,NTD,SFY)                   M12350
      IF (SFY.EQ.1.) NTR = NT3R                                           M12360
      IF (IFUNCT.EQ.2) CALL PLTLIN (TITL3R,BLK,NTR,SFY)                   M12370
      GO TO 150                                                           M12380
C                                                                         M12390
  110 CALL PLTLIN (TITL6,BLK,NT6,1.0)                                     M12400
      GO TO 150                                                           M12410
  120 CALL PLTLOG (TITL6,BLK,NT6)                                         M12420
      GO TO 150                                                           M12430
  130 IF (IFUNCT.EQ.0) CALL PLTEMP (TITL4,BLK,NT4)                        M12440
      IF (IFUNCT.EQ.1) CALL PLTEMP (TITL4D,BLK,NT4D)                      M12450
      IF (IFUNCT.EQ.2) CALL PLTEMP (TITL4R,BLK,NT4R)                      M12460
      GO TO 150                                                           M12470
  140 IREJ = 1                                                            M12480
      WRITE (LOUT,900) IGO                                                M12490
C                                                                         M12500
      RETURN                                                              M12510
C                                                                         M12520
C     WRITE DATE AND TIME TO BOTTOM OF PLOT                               M12530
C                                                                         M12540
  150 YBIAS = -3.75*1.25*HGT-0.25                                         M12550
      XBIAS = YBIAS-0.5                                                   M12560
      CALL FDATE (HDATE)                                                  M12570
      CALL FTIME (HTIME)                                                  M12580
      WRITE (YIDC,'(2(A10),A3,2(A10))') HDATE,HTIME,ASTR,(YID(I),I=1,2)   M12590
      CALL SYMBOL (XBIAS,YBIAS,0.10,YIDC,0.0,43)                          M12600
C                                                                         M12610
      RETURN                                                              M12620
C                                                                         M12630
  900 FORMAT (10X,'PLOT REQUEST INCOMPATIBLE WITH LBLRMT, IGO =',I2)      M12640
C                                                                         M12650
      END                                                                 M12660
      SUBROUTINE BBSCLE                                                   M12670
C                                                                         M12680
      IMPLICIT DOUBLE PRECISION (V)                                     ! M12690
C                                                                         M12700
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       M12710
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,     M12720
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      M12730
     *                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES                      M12740
C                                                                         M12750
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           M12760
C                                                                         M12770
C     BBFN  BLACK BODY                                                    M12780
C                                                                         M12790
      XKTMX = YMAX/RADCN2                                                 M12800
      XKTMN = YMIN/RADCN2                                                 M12810
      DVDUM = 1.                                                          M12820
      VIDUM = V2                                                          M12830
      BBLST1 = -1.                                                        M12840
      BBLST2 = -1.                                                        M12850
      YMAX1 = BBFN(V2,DVDUM,XKTMX,VIDUM,BBDEL,BBLST2)                     M12860
      YMIN1 = 0.                                                          M12870
      IF (YMIN.GT.0.) THEN                                                M12880
         YMIN1 = BBFN(V1,DVDUM,XKTMN,VIDUM,BBDEL,BBLST1)                  M12890
         RATYF = 0.                                                       M12900
      ELSE                                                                M12910
         RATYF = YMIN/YMAX                                                M12920
      ENDIF                                                               M12930
      NS = ALOG10(YMAX1)                                                  M12940
      IF (YMIN1.GT.0.) THEN                                               M12950
         NB = ALOG10(YMIN1)                                               M12960
      ELSE                                                                M12970
         NB = NS-5                                                        M12980
      ENDIF                                                               M12990
C                                                                         M13000
C     LINEAR RADIANCE                                                     M13010
C                                                                         M13020
      IF (LOGPLT.EQ.0) THEN                                               M13030
         YMAX = 10.**NS                                                   M13040
         YMIN = 10.**NB                                                   M13050
         ZMAX = YMAX1/YMAX                                                M13060
         ZMIN = YMIN1/YMIN                                                M13070
         IF (ZMAX.GT.0.1.AND.ZMAX.LE.0.2) YMAX = .2*YMAX                  M13080
         IF (ZMAX.GT.0.2.AND.ZMAX.LE.0.5) YMAX = .5*YMAX                  M13090
C                                                                         M13100
C        IF(ZMAX.GT.0.5.AND.ZMAX.LE.1.0) YMAX=YMAX                        M13110
C                                                                         M13120
         IF (ZMIN.GT.0.1.AND.ZMIN.LE.0.2) YMIN = .1*YMIN                  M13130
         IF (ZMIN.GT.0.2.AND.ZMIN.LE.0.5) YMIN = .2*YMIN                  M13140
         IF (ZMIN.GT.0.5.AND.ZMIN.LE.1.0) YMIN = .5*YMIN                  M13150
         IF (RATYF.NE.0..AND.IFUNCT.EQ.1) YMIN = RATYF*YMAX               M13160
         DELY = (YMAX-YMIN)/10.                                           M13170
         NUMDVY = 10                                                      M13180
      ELSE                                                                M13190
C                                                                         M13200
C     LOG RADIANCE                                                        M13210
C                                                                         M13220
         YMAX = NS                                                        M13230
         IF (NB.EQ.NS) NB = NS-1                                          M13240
         YMIN = NB                                                        M13250
         NUMDVY = NS-NB                                                   M13260
         DELY = 1.                                                        M13270
      ENDIF                                                               M13280
C                                                                         M13290
      RETURN                                                              M13300
C                                                                         M13310
      END                                                                 M13320
      SUBROUTINE PLTLIN (YTITLE,BLK,NY,SFY)                               M13330
C                                                                         M13340
      IMPLICIT DOUBLE PRECISION (V)                                     ! M13350
C                                                                         M13360
      DIMENSION YTITLE(5)                                                 M13370
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       M13380
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,     M13390
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      M13400
     *                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES                      M13410
C                                                                         M13420
      CALL AXISL (0.0,0.0,YTITLE,NY,NUMDVY,DIVLNY,NUMSBY,YMIN*SFY,        M13430
     *            DELY*SFY,IDEC,90.0,HGT,1,1,NOENDY,0,0)                  M13440
      CALL AXISL (XSIZE,0.0,BLK,-1,NUMDVY,DIVLNY,NUMSBY,YMIN*SFY,         M13450
     *            DELY*SFY,IDEC,90.0,HGT,1,1,NOENDY,1,0)                  M13460
C                                                                         M13470
      RETURN                                                              M13480
C                                                                         M13490
      END                                                                 M13500
      SUBROUTINE PLTEMP (YTITLE,BLK,NY)                                   M13510
C                                                                         M13520
      IMPLICIT DOUBLE PRECISION (V)                                     ! M13530
C                                                                         M13540
      DIMENSION YTITLE(4)                                                 M13550
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       M13560
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,     M13570
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      M13580
     *                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES                      M13590
C                                                                         M13600
      CALL AXISL (0.0,0.0,YTITLE,NY,NUMDVY,DIVLNY,NUMSBY,YMIN,DELY,       M13610
     *            IDEC,90.0,HGT,1,1,NOENDY,0,0)                           M13620
      IF (YMIN.LE.0.0) THEN                                               M13630
         CALL AXISL (XSIZE,0.0,BLK,-1,NUMDVY,DIVLNY,NUMSBY,YMIN,DELY,     M13640
     *               IDEC,90.0,HGT,1,1,NOENDY,1,0)                        M13650
      ELSE                                                                M13660
         CALL AX2                                                         M13670
      ENDIF                                                               M13680
C                                                                         M13690
      RETURN                                                              M13700
C                                                                         M13710
      END                                                                 M13720
      SUBROUTINE PLTLOG (YTITLE,BLK,NY)                                   M13730
C                                                                         M13740
      IMPLICIT DOUBLE PRECISION (V)                                     ! M13750
C                                                                         M13760
      DIMENSION YTITLE(4)                                                 M13770
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       M13780
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,     M13790
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      M13800
     *                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES                      M13810
C                                                                         M13820
      CALL AXLOG (0.0,0.0,YTITLE,NY,NUMDVY,DIVLNY,YMIN,90.0,HGT,1,1,      M13830
     *            NOENDY,0,0)                                             M13840
      CALL AXLOG (XSIZE,0.0,BLK,-1,NUMDVY,DIVLNY,YMIN,90.0,HGT,1,1,       M13850
     *            NOENDY,1,0)                                             M13860
C                                                                         M13870
      RETURN                                                              M13880
C                                                                         M13890
      END                                                                 M13900
      SUBROUTINE AX2                                                      M13910
C                                                                         M13920
      IMPLICIT DOUBLE PRECISION (V)                                     ! M13930
C                                                                         M13940
      CHARACTER TITL3*50,HTEN*2                                           M13950
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       M13960
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,     M13970
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      M13980
     *                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES                      M13990
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           M14000
C                                                                         M14010
      DATA HTEN / '10'/,                                                  M14020
     *     TITL3 / ' RADIANCE WATTS /  (CM**2*STER*CM-1)  * '/,           M14030
     *     NT3 / 37 /                                                     M14040
C                                                                         M14050
      BB(X,TAVE) = (RADCN1*X**3)/(EXP(RADCN2*X/TAVE)-1.0)                 M14060
      TEM(X,BBY) = RADCN2*X/ALOG(RADCN1*X**3/BBY+1.0)                     M14070
C                                                                         M14080
      ISKIP = 0                                                           M14090
      CALL PLOT (XSIZE,0.0,3)                                             M14100
      CALL PLOT (XSIZE,0.0,2)                                             M14110
      CALL PLOT (XSIZE,YSIZE,2)                                           M14120
      CALL PLOT (XSIZE,YSIZE,3)                                           M14130
      X2 = XSIZE+HGT                                                      M14140
      X3 = XSIZE+3.*HGT                                                   M14150
      YOLD = YSIZE+1.001*HGT                                              M14160
      HGTE = HGT*(2./3.)                                                  M14170
      S2 = V2                                                             M14180
      BB1 = BB(S2,YMAX)                                                   M14190
      IB1 = ALOG10(BB1)                                                   M14200
      B1 = 10.**IB1                                                       M14210
   10 DO 30 I = 1, 5                                                      M14220
         FN = 2*(5-I)                                                     M14230
         IF (I.EQ.5) FN = 1.                                              M14240
         BC = B1*FN                                                       M14250
         TMP = TEM(S2,BC)                                                 M14260
         IF (TMP.GT.YMAX) GO TO 30                                        M14270
         YIN = (TMP-YMIN)/DY                                              M14280
         IF (YIN.LT.HGT) GO TO 40                                         M14290
         DELY = YOLD-YIN                                                  M14300
         IF (DELY.LE.HGT) GO TO 20                                        M14310
         IF ((ISKIP.EQ.1).AND.(I.NE.5)) GO TO 30                          M14320
         YOLD = YIN                                                       M14330
         CALL SYMBOL (XSIZE,YIN,HGT,3,0.0,-1)                             M14340
         IF (I.NE.5) GO TO 30                                             M14350
         YIN = YIN-HGT/2.                                                 M14360
         CALL SYMBOL (X2,YIN,HGT,HTEN,0.0,2)                              M14370
         HIP = IB1                                                        M14380
         YINE = YIN+HGTE                                                  M14390
         CALL NUMBER (X3,YINE,HGTE,HIP,0.0,-1)                            M14400
         GO TO 30                                                         M14410
   20    ISKIP = 1                                                        M14420
   30 CONTINUE                                                            M14430
      B1 = B1/10.                                                         M14440
      IB1 = IB1-1                                                         M14450
      GO TO 10                                                            M14460
   40 YLT = NT3*HGT*.333                                                  M14470
      YHV = YSIZE/2.                                                      M14480
      YIN = YHV-YLT                                                       M14490
      X4 = X3+4.5*HGT                                                     M14500
      CALL SYMBOL (X4,YIN,HGT,TITL3,90.0,NT3)                             M14510
C                                                                         M14520
      RETURN                                                              M14530
C                                                                         M14540
      END                                                                 M14550
      SUBROUTINE FPLINE (NLO,NPTS)                                        M14560
C                                                                         M14570
      IMPLICIT DOUBLE PRECISION (V)                                     ! M14580
C                                                                         M14590
      COMMON XX(2450),YY(2450)                                            M14600
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       M14610
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELX,NUMDVY,NUMSBY,     M14620
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      M14630
     *                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES                      M14640
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M14650
      DIMENSION X(2450)                                                   M14660
C                                                                         M14670
      EQUIVALENCE (X(1),XX(1))                                            M14680
C                                                                         M14690
C     FOUR POINT LAGRANGE INTERPOLATION:                                  M14700
C                                                                         M14710
      IF (NLO.EQ.4) GO TO 10                                              M14720
C                                                                         M14730
C     INITIALIZE FIRST PANEL:                                             M14740
C                                                                         M14750
      YY(1) = YY(2)                                                       M14760
      DV1 = DV/4.                                                         M14770
      DV2 = 2.*DV1                                                        M14780
      DV3 = 3.*DV1                                                        M14790
      X00 = -7./128.                                                      M14800
      X01 = 105./128.                                                     M14810
      X02 = 35./128.                                                      M14820
      X03 = -5./128.                                                      M14830
      X10 = -1./16.                                                       M14840
      X11 = 9./16.                                                        M14850
      IF (NPTS.LT.4) GO TO 50                                             M14860
C                                                                         M14870
   10 XXI = XX(NLO)                                                       M14880
      ILO = 2                                                             M14890
   20 IHI = ILO+600                                                       M14900
      IF (IHI.GE.NPTS) IHI = NPTS-1                                       M14910
      IMAX = IHI-1                                                        M14920
      IP0 = 1                                                             M14930
      DO 30 I = ILO, IMAX                                                 M14940
         XI = XXI+DV*FLOAT(I-NLO)                                         M14950
         X(IP0) = XI                                                      M14960
         X(IP0+1) = XI+DV1                                                M14970
         X(IP0+2) = XI+DV2                                                M14980
         X(IP0+3) = XI+DV3                                                M14990
         Y(IP0) = YY(I)                                                   M15000
         Y(IP0+1) = X00*YY(I-1)+X01*YY(I)+X02*YY(I+1)+X03*YY(I+2)         M15010
         Y(IP0+2) = X10*(YY(I-1)+YY(I+2))+X11*(YY(I)+YY(I+1))             M15020
         Y(IP0+3) = X03*YY(I-1)+X02*YY(I)+X01*YY(I+1)+X00*YY(I+2)         M15030
         IP0 = IP0+4                                                      M15040
   30 CONTINUE                                                            M15050
      X(IP0) = XI+DV                                                      M15060
      Y(IP0) = YY(IHI)                                                    M15070
      J1 = 1                                                              M15080
      CALL FSCLIN (J1,IP0,X,Y)                                            M15090
      IF (IHI+1.EQ.NPTS) GO TO 40                                         M15100
      ILO = IHI                                                           M15110
      GO TO 20                                                            M15120
C                                                                         M15130
C     PRESERVE LAST THREE POINTS OF PANEL FOR SUBSEQUENT PANEL:           M15140
C                                                                         M15150
   40 YY(1) = YY(IMAX)                                                    M15160
      YY(2) = YY(IHI)                                                     M15170
      YY(3) = YY(NPTS)                                                    M15180
C                                                                         M15190
      RETURN                                                              M15200
C                                                                         M15210
C     CASE IN WHICH INITIAL PANEL HAS LESS THAN THREE POINTS:             M15220
C                                                                         M15230
   50 YY(3) = YY(NPTS)                                                    M15240
      YY(2) = YY(NPTS-1)                                                  M15250
      YY(1) = YY(2)                                                       M15260
C                                                                         M15270
      RETURN                                                              M15280
C                                                                         M15290
      END                                                                 M15300
      SUBROUTINE LINT (NST,NND,J)                                         M15310
C                                                                         M15320
      IMPLICIT DOUBLE PRECISION (V)                                     ! M15330
C                                                                         M15340
      COMMON XX(2450),YY(2450)                                            M15350
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M15360
C                                                                         M15370
      DO 10 I = NST, NND                                                  M15380
         YY(J) = Y(I)                                                     M15390
         J = J+1                                                          M15400
   10 CONTINUE                                                            M15410
C                                                                         M15420
      RETURN                                                              M15430
C                                                                         M15440
      END                                                                 M15450
      SUBROUTINE EXPT (NST,NND,J)                                         M15460
C                                                                         M15470
      IMPLICIT DOUBLE PRECISION (V)                                     ! M15480
C                                                                         M15490
      COMMON XX(2450),YY(2450)                                            M15500
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M15510
C                                                                         M15520
      DO 10 I = NST, NND                                                  M15530
         YY(J) = 0.                                                       M15540
         IF (Y(I).GT.20.) GO TO 10                                        M15550
         YY(J) = EXP(-Y(I))                                               M15560
   10    J = J+1                                                          M15570
C                                                                         M15580
      RETURN                                                              M15590
C                                                                         M15600
      END                                                                 M15610
      SUBROUTINE XNTLOG (NST,NND,J)                                       M15620
C                                                                         M15630
      IMPLICIT DOUBLE PRECISION (V)                                     ! M15640
C                                                                         M15650
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           M15660
      COMMON XX(2450),YY(2450)                                            M15670
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M15680
C                                                                         M15690
      DO 10 I = NST, NND                                                  M15700
         IF (Y(I).LE.0.) Y(I) = EXPMIN                                    M15710
         YY(J) = -ALOG(Y(I))                                              M15720
         J = J+1                                                          M15730
   10 CONTINUE                                                            M15740
C                                                                         M15750
      RETURN                                                              M15760
C                                                                         M15770
      END                                                                 M15780
      SUBROUTINE TEMPFN (NST,NND,J)                                       M15790
C                                                                         M15800
      IMPLICIT DOUBLE PRECISION (V)                                     ! M15810
C                                                                         M15820
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           M15830
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           M15840
      COMMON XX(2450),YY(2450)                                            M15850
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M15860
C                                                                         M15870
      DO 10 I = NST, NND                                                  M15880
         YY(J) = 0.                                                       M15890
         IF (Y(I).GE.EXPMIN) THEN                                         M15900
            X = RADCN1*(XX(J)**3)/Y(I)                                    M15910
            IF (X.GE.EXPMIN) THEN                                         M15920
               Z = ALOG(X+1.)                                             M15930
               YY(J) = RADCN2*(XX(J)/Z)                                   M15940
            ENDIF                                                         M15950
         ENDIF                                                            M15960
         J = J+1                                                          M15970
   10 CONTINUE                                                            M15980
C                                                                         M15990
      RETURN                                                              M16000
C                                                                         M16010
      END                                                                 M16020
      SUBROUTINE TENLOG (NST,NND,J)                                       M16030
C                                                                         M16040
      IMPLICIT DOUBLE PRECISION (V)                                     ! M16050
C                                                                         M16060
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           M16070
      COMMON XX(2450),YY(2450)                                            M16080
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M16090
C                                                                         M16100
      DO 10 I = NST, NND                                                  M16110
         IF (Y(I).LE.0.) Y(I) = EXPMIN                                    M16120
         YY(J) = ALOG10(Y(I))                                             M16130
         J = J+1                                                          M16140
   10 CONTINUE                                                            M16150
C                                                                         M16160
      RETURN                                                              M16170
C                                                                         M16180
      END                                                                 M16190
      SUBROUTINE DBOD (NST,NND,J)                                         M16200
C                                                                         M16210
      IMPLICIT DOUBLE PRECISION (V)                                     ! M16220
C                                                                         M16230
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           M16240
      COMMON XX(2450),YY(2450)                                            M16250
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M16260
C                                                                         M16270
      DATA CON / 4.3429448 /                                              M16280
C                                                                         M16290
      DO 10 I = NST, NND                                                  M16300
         YY(J) = Y(I)*CON                                                 M16310
         J = J+1                                                          M16320
   10 CONTINUE                                                            M16330
C                                                                         M16340
      RETURN                                                              M16350
C                                                                         M16360
      END                                                                 M16370
      SUBROUTINE DBTR (NST,NND,J)                                         M16380
C                                                                         M16390
      IMPLICIT DOUBLE PRECISION (V)                                     ! M16400
C                                                                         M16410
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           M16420
      COMMON XX(2450),YY(2450)                                            M16430
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M16440
C                                                                         M16450
      DATA CON / 4.3429448 /                                              M16460
C                                                                         M16470
      DO 10 I = NST, NND                                                  M16480
         IF (Y(I).LE.0.) Y(I) = EXPMIN                                    M16490
         YY(J) = -CON*ALOG(Y(I))                                          M16500
         J = J+1                                                          M16510
   10 CONTINUE                                                            M16520
C                                                                         M16530
      RETURN                                                              M16540
C                                                                         M16550
      END                                                                 M16560
      SUBROUTINE LDBOPD (NST,NND,J)                                       M16570
C                                                                         M16580
      IMPLICIT DOUBLE PRECISION (V)                                     ! M16590
C                                                                         M16600
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           M16610
      COMMON XX(2450),YY(2450)                                            M16620
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M16630
C                                                                         M16640
      DATA CON / 4.3429448 /                                              M16650
C                                                                         M16660
      DO 10 I = NST, NND                                                  M16670
         IF (Y(I).LE.0.) Y(I) = EXPMIN                                    M16680
         YY(J) = ALOG10(Y(I)*CON)                                         M16690
         J = J+1                                                          M16700
   10 CONTINUE                                                            M16710
C                                                                         M16720
      RETURN                                                              M16730
C                                                                         M16740
      END                                                                 M16750
      SUBROUTINE LDBTR (NST,NND,J)                                        M16760
C                                                                         M16770
      IMPLICIT DOUBLE PRECISION (V)                                     ! M16780
C                                                                         M16790
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           M16800
      COMMON XX(2450),YY(2450)                                            M16810
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M16820
C                                                                         M16830
      DATA CON / 4.3429448 /                                              M16840
C                                                                         M16850
      DO 10 I = NST, NND                                                  M16860
         IF (Y(I).LE.0.) Y(I) = EXPMIN                                    M16870
         YY(J) = -CON*ALOG(Y(I))                                          M16880
         IF (YY(J).LE.0.) YY(J) = EXPMIN                                  M16890
         YY(J) = ALOG10(YY(J))                                            M16900
         J = J+1                                                          M16910
   10 CONTINUE                                                            M16920
C                                                                         M16930
      RETURN                                                              M16940
C                                                                         M16950
      END                                                                 M16960
      SUBROUTINE LINSTC (NST,NND,J,XLOGCN)                                M16970
C                                                                         M16980
      IMPLICIT DOUBLE PRECISION (V)                                     ! M16990
C                                                                         M17000
      COMMON XX(2450),YY(2450)                                            M17010
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M17020
C                                                                         M17030
      DO 10 I = NST, NND                                                  M17040
         YY(J) = Y(I)*XLOGCN                                              M17050
         J = J+1                                                          M17060
   10 CONTINUE                                                            M17070
C                                                                         M17080
      RETURN                                                              M17090
C                                                                         M17100
      END                                                                 M17110
      SUBROUTINE XLOGLN (NST,NND,J)                                       M17120
C                                                                         M17130
      IMPLICIT DOUBLE PRECISION (V)                                     ! M17140
C                                                                         M17150
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           M17160
      COMMON XX(2450),YY(2450)                                            M17170
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M17180
C                                                                         M17190
      DO 10 I = NST, NND                                                  M17200
         CAY = Y(I)                                                       M17210
         IF (CAY.EQ.0) CAY = EXPMIN                                       M17220
         CAY = -ALOG(CAY)                                                 M17230
         IF (CAY.EQ.0) CAY = EXPMIN                                       M17240
         YY(J) = ALOG10(CAY)                                              M17250
         J = J+1                                                          M17260
   10 CONTINUE                                                            M17270
C                                                                         M17280
      RETURN                                                              M17290
C                                                                         M17300
      END                                                                 M17310
      SUBROUTINE MNMX (NST,NND,JMIN,JMAX)                                 M17320
C                                                                         M17330
      IMPLICIT DOUBLE PRECISION (V)                                     ! M17340
C                                                                         M17350
      COMMON XX(2450),YY(2450)                                            M17360
      COMMON /YCOM/ V1P,V2P,DV,NLIM,Y(2502)                               M17370
      J = JMIN                                                            M17380
      DO 10 I = NST, NND                                                  M17390
         IF (YY(J).LT.YY(JMIN)) JMIN = J                                  M17400
         IF (YY(J).GT.YY(JMAX)) JMAX = J                                  M17410
         J = J+1                                                          M17420
   10 CONTINUE                                                            M17430
      RETURN                                                              M17440
      END                                                                 M17450
      SUBROUTINE FSCLIN (NLO,NPTS,XX,YY)                                  M17460
C                                                                         M17470
      IMPLICIT DOUBLE PRECISION (V)                                     ! M17480
C                                                                         M17490
      DIMENSION XX(*),YY(*),PNLHDR(2)                                     M17500
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         M17510
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        M17520
     *              NLTEFL,LNFIL4,LNGTH4                                  M17530
      COMMON /JPLTFL/ V1P,V2P,DVP,NLIM                                    M17540
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       M17550
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELX,NUMDVY,NUMSBY,     M17560
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      M17570
     *                JOUT,JPLTFL,JHDR,IFUNCT,NOAXES                      M17580
C                                                                         M17590
      EQUIVALENCE (PNLHDR(1),V1P) , (JPLT,NUMSBX) , (LPLT,NOENDX)         M17600
C                                                                         M17610
      YYSTOR = YY(NPTS)                                                   M17620
C                                                                         M17630
      IF (JOUT.GE.1) THEN                                                 M17640
         V1P = XX(NLO)                                                    M17650
         V2P = XX(NPTS)                                                   M17660
         DVP = (V2P-V1P)/FLOAT(NPTS-NLO)                                  M17670
         NLOW = NLO-1                                                     M17680
         NLIM = NPTS-NLOW                                                 M17690
         IF (JOUT.LE.2) THEN                                              M17700
            CALL BUFOUT (JPLTFL,PNLHDR(1),NPHDRF)                         M17710
            CALL BUFOUT (JPLTFL,YY(NLO),NPTS)                             M17720
         ELSE                                                             M17730
            DO 10 II = NLO, NPTS                                          M17740
               VOUT = V1P+FLOAT(II-NLO)*DVP                               M17750
               WRITE (JPLTFL,900) VOUT,YY(II)                             M17760
   10       CONTINUE                                                      M17770
         ENDIF                                                            M17780
         IF (JOUT.EQ.1.OR.JOUT.EQ.3) THEN                                 M17790
            XX(1) = XX(NPTS)                                              M17800
            YY(1) = YYSTOR                                                M17810
            RETURN                                                        M17820
         ENDIF                                                            M17830
      ENDIF                                                               M17840
C                                                                         M17850
C   TRUNCATE POINTS OUTSIDE PLOTTING LIMITS:                              M17860
C                                                                         M17870
      DO 20 J = 1, NPTS                                                   M17880
         YY(J) = MAX(YY(J),YMIN)                                          M17890
         YY(J) = MIN(YY(J),YMAX)                                          M17900
   20 CONTINUE                                                            M17910
C                                                                         M17920
      XX(NPTS+1) = V1                                                     M17930
      YY(NPTS+1) = YMIN                                                   M17940
      XX(NPTS+2) = DX                                                     M17950
      YY(NPTS+2) = DY                                                     M17960
      IF (NOAXES.EQ.1) THEN                                               M17970
         JCAL = JPLT                                                      M17980
         LCAL = LPLT                                                      M17990
      ELSE                                                                M18000
         JCAL = 0                                                         M18010
         LCAL = 0                                                         M18020
      ENDIF                                                               M18030
C                                                                         M18040
      CALL LINE (XX,YY,NPTS,1,JCAL,LCAL)                                  M18050
C                                                                         M18060
C     USE ABOVE CALL TO LINE FOR STANDARD CALCOMP SYSTEM.                 M18070
C     CALL BELOW MAKES USE OF MODIFIED CALCOMP ROUTINE                    M18080
C     AT GEOPHYSICS DIRECTORATE OF THE PHILLIPS LABORATIORY.              M18090
C                                                                         M18100
C>    CALL LINE(XX,YY,NPTS,1,JCAL,LCAL,V1,DX,YMIN,DY,0.08)              > M18110
C                                                                         M18120
      XX(1) = XX(NPTS)                                                    M18130
      YY(1) = YYSTOR                                                      M18140
C                                                                         M18150
      RETURN                                                              M18160
C                                                                         M18170
  900 FORMAT (3X,F15.8,4X,G19.8)                                          M18180
C                                                                         M18190
      END                                                                 M18200
      SUBROUTINE AXISL (X,Y,BCD,N,NUMDIV,DIVLEN,NUMSUB,BEGNUM,DELNUM,     M18210
     *                  NUMDEC,THETA,HEIGHT,NRPT,NTURN,NOEND,LSUPR,       M18220
     *                  LTURN)                                            M18230
C                                                                         M18240
C   MODIFIED VERSION OF AXISL                                             M18250
C                                                                         M18260
C   WRITTEN BY RICHARD L. TAYLOR   RADC/ET   EEC   NOVEMBER 1980          M18270
C                                                                         M18280
C   ROUTINE TO PLOT A LABELLED LINEAR AXIS                                M18290
C                                                                         M18300
C                                                                         M18310
C   X AND Y ARE THE STARTING COORDINATES OF THE AXIS RELATIVE TO THE      M18320
C                                                                         M18330
C        CURRENT ORIGIN                                                   M18340
C                                                                         M18350
C   BCD IS THE LABEL OF THE AXIS EXPRESSED AS A HOLLERITH CONSTANT        M18360
C                                                                         M18370
C   N IS THE NUMBER OF CHARACTERS IN THE LABEL                            M18380
C                                                                         M18390
C        NEGATIVE N PLACES THE LABEL ON THE CLOCKWISE SIDE OF THE AXIS    M18400
C                                                                         M18410
C        POSITIVE N PLACES THE LABEL ON THE COUNTERCLOCKWISE SIDE         M18420
C                                                                         M18430
C   NUMDIV IS THE NUMBER OF MAJOR DIVISIONS                               M18440
C                                                                         M18450
C   DIVLEN IS THE LENGTH IN INCHES OF A MAJOR DIVISION                    M18460
C                                                                         M18470
C   NUMSUB IS THE NUMBER OF MINOR DIVISIONS PER MAJOR DIVISION            M18480
C                                                                         M18490
C        1 GIVES NO SUBDIVISION TICS, 2 GIVES ONE SUBDIVISION TIC, ETC.   M18500
C                                                                         M18510
C   BEGNUM IS THE NUMBER FOR THE BEGINNING OF THE AXIS                    M18520
C                                                                         M18530
C   DELNUM IS THE DELTA NUMBER FOR A MAJOR DIVISION                       M18540
C                                                                         M18550
C   NUMDEC IS THE NUMBER OF DECIMAL PLACES DESIRED                        M18560
C                                                                         M18570
C        NUMDEC EQUAL TO -1 SUPPRESSES THE DECIMAL POINT                  M18580
C                                                                         M18590
C   THETA IS THE ANGLE OF THE AXIS IN DEGREES  (0.0 FOR X, 90.0 FOR Y)    M18600
C                                                                         M18610
C   HEIGHT IS THE HEIGHT OF THE NUMBERS IN INCHES                         M18620
C                                                                         M18630
C   NRPT IS THE REPEAT FACTOR FOR THE SCALE NUMBERS (USUALLY INTEGER 1)   M18640
C                                                                         M18650
C        WHEN NRPT IS ZERO THE SCALE NUMBERS WILL BE SUPPRESSED;          M18660
C                                                                         M18670
C        WHEN NRPT = 2, EVERY 2ND SCALE NUMBER WILL BE PRODUCED; ETC.     M18680
C                                                                         M18690
C   NTURN EQUAL TO 1 TURNS THE AXIS NUMBERS BY 90 DEGREES CLOCKWISE,      M18700
C        -1 TURNS NUMBERS BY 90 DEGREES COUNTERCLOCKWISE, 0 FOR NO TURN   M18710
C                                                                         M18720
C   NOEND EQUAL TO 1 SUPPRESSES THE NUMBERS AT EITHER END OF THE AXIS,    M18730
C                                                                         M18740
C        2 SUPPRESSES THE BEGINNING NUMBER, 3 THE ENDING NUMBER           M18750
C                                                                         M18760
C   LSUPR EQUAL TO 1 SUPPRESSES THE LABEL                                 M18770
C    LTURN NOT USED                                                       M18780
C                                                                         M18790
      COMMON /TITLOC/ XPOS,YPOS                                           M18800
      DIMENSION BCD(*)                                                    M18810
C                                                                         M18820
      THETA1 = THETA-90.*NTURN                                            M18830
      PI = 2.*ASIN(1.)                                                    M18840
      ANGLE = (PI/180.)*THETA                                             M18850
      SINANG = SIN(ANGLE)                                                 M18860
      COSANG = COS(ANGLE)                                                 M18870
      SIGNAX = FLOAT(ISIGN(1,N))                                          M18880
      SIZMAJ = 0.25*HEIGHT+0.05                                           M18890
      OFFST = HEIGHT*1.5                                                  M18900
      DXMAJ = -SIZMAJ*SINANG*SIGNAX                                       M18910
      DYMAJ = SIZMAJ*COSANG*SIGNAX                                        M18920
      DXMIN = 0.5*DXMAJ                                                   M18930
      DYMIN = 0.5*DYMAJ                                                   M18940
      NSUB = NUMSUB                                                       M18950
      IF (NUMSUB.LT.1) NSUB = 1                                           M18960
      SUBDIV = DIVLEN/FLOAT(NSUB)                                         M18970
      BCDSIZ = 1.25*HEIGHT                                                M18980
      YBIAS = (-0.50+SIGN(1.25,SIGNAX))*HEIGHT+DYMAJ                      M18990
      NABS = IABS(N)                                                      M19000
      BCDLEN = (FLOAT(NABS)-0.4)*BCDSIZ                                   M19010
      S = DIVLEN*FLOAT(NUMDIV)                                            M19020
      DIVCOS = DIVLEN*COSANG                                              M19030
      DIVSIN = DIVLEN*SINANG                                              M19040
      SIZMAX = HEIGHT                                                     M19050
C                                                                         M19060
C   DRAW DIVISION NUMBERS                                                 M19070
C                                                                         M19080
      NDIV = 0                                                            M19090
   10 DIGITS = 0.0                                                        M19100
      XTIC = X+DIVCOS*FLOAT(NDIV)                                         M19110
      YTIC = Y+DIVSIN*FLOAT(NDIV)                                         M19120
      IF (NRPT.EQ.0) GO TO 30                                             M19130
      NSUPR = NDIV-(NDIV/NRPT)*NRPT                                       M19140
      IF (NSUPR.NE.0) GO TO 30                                            M19150
      IF ((NOEND.EQ.1.OR.NOEND.EQ.2).AND.NDIV.EQ.0) GO TO 30              M19160
      IF ((NOEND.EQ.1.OR.NOEND.EQ.3).AND.NDIV.EQ.NUMDIV) GO TO 30         M19170
      DIVNUM = BEGNUM+DELNUM*FLOAT(NDIV)                                  M19180
      IF (ABS(DIVNUM).GE.10.0) DIGITS = ALOG10(ABS(DIVNUM))               M19190
      DIGITS = AINT(DIGITS+1.0E-12)                                       M19200
      IF (DIVNUM.LT.0.0) DIGITS = DIGITS+1.0                              M19210
      SIZNUM = (DIGITS+FLOAT(NUMDEC)+1.7)*HEIGHT                          M19220
      XBIAS = -0.5*SIZNUM                                                 M19230
      XBIAS1 = 0.                                                         M19240
      YBIAS1 = 0.                                                         M19250
      IF (NTURN.EQ.0) GO TO 20                                            M19260
      YBIAS1 = YBIAS-SIZNUM-OFFST                                         M19270
      IF (N.LT.0) YBIAS1 = YBIAS+OFFST                                    M19280
      XBIAS1 = XBIAS+HEIGHT*0.5                                           M19290
   20 XPOS = XTIC-YBIAS*SINANG+XBIAS*COSANG+YBIAS1*SINANG-XBIAS1*COSANG   M19300
      YPOS = YTIC+YBIAS*COSANG+XBIAS*SINANG-XBIAS1*SINANG-YBIAS1*COSANG   M19310
      CALL NUMBER (XPOS,YPOS,HEIGHT,DIVNUM,THETA1,NUMDEC)                 M19320
      SIZMAX = MAX(SIZMAX,SIZNUM)                                         M19330
   30 CONTINUE                                                            M19340
C                                                                         M19350
C   DRAW TIC MARKS                                                        M19360
C                                                                         M19370
      CALL PLOT (XTIC,YTIC,3)                                             M19380
      CALL PLOT (XTIC+DXMAJ,YTIC+DYMAJ,2)                                 M19390
      IF (NDIV.EQ.NUMDIV) GO TO 60                                        M19400
      IF (NUMSUB.LE.1) GO TO 50                                           M19410
      DO 40 J = 2, NUMSUB                                                 M19420
         SUBLEN = SUBDIV*FLOAT(J-1)                                       M19430
         XSTIC = XTIC+SUBLEN*COSANG                                       M19440
         YSTIC = YTIC+SUBLEN*SINANG                                       M19450
         CALL PLOT (XSTIC+DXMIN,YSTIC+DYMIN,3)                            M19460
         CALL PLOT (XSTIC,YSTIC,2)                                        M19470
   40 CONTINUE                                                            M19480
   50 NDIV = NDIV+1                                                       M19490
      GO TO 10                                                            M19500
C                                                                         M19510
C   DRAW AXIS                                                             M19520
C                                                                         M19530
   60 CALL PLOT (XTIC,YTIC,3)                                             M19540
      CALL PLOT (X,Y,2)                                                   M19550
C                                                                         M19560
C   DRAW LABEL                                                            M19570
C                                                                         M19580
      IF (LSUPR.EQ.1.OR.NABS.EQ.0) RETURN                                 M19590
      XBIAS1 = 0                                                          M19600
      YBIAS1 = -SIZNUM-OFFST                                              M19610
      IF (N.LT.0) YBIAS1 = -YBIAS1                                        M19620
      IF (NTURN.EQ.0) YBIAS1 = 0                                          M19630
      XBIAS = 0.5*(S-BCDLEN)                                              M19640
      YBIAS = (-0.50+SIGN(3.25,SIGNAX))*BCDSIZ                            M19650
      XPOS = X-YBIAS*SINANG+XBIAS*COSANG+YBIAS1*SINANG-XBIAS1*COSANG      M19660
      YPOS = Y+YBIAS*COSANG+XBIAS*SINANG-XBIAS1*SINANG-YBIAS1*COSANG      M19670
      CALL SYMBOL (XPOS,YPOS,BCDSIZ,BCD,THETA,NABS)                       M19680
C                                                                         M19690
      RETURN                                                              M19700
C                                                                         M19710
      END                                                                 M19720
      SUBROUTINE AXLOG (X,Y,BCD,N,NUMCYC,CYCLEN,BEGEXP,THETA,HEIGHT,      M19730
     *                  NRPT,NTURN,NOEND,LSUPR,LTURN)                     M19740
C                                                                         M19750
C                                                                         M19760
C   WRITTEN BY RICHARD L. TAYLOR   RADC/ET   EEC   NOVEMBER 1980          M19770
C                                                                         M19780
C                                                                         M19790
C                                                                         M19800
C   ROUTINE TO PLOT A LABELLED LOGARITHMIC AXIS                           M19810
C                                                                         M19820
C                                                                         M19830
C                                                                         M19840
C                                                                         M19850
C   X AND Y ARE THE STARTING COORDINATES OF THE AXIS RELATIVE TO THE      M19860
C                                                                         M19870
C        CURRENT ORIGIN                                                   M19880
C                                                                         M19890
C   BCD IS THE LABEL OF THE AXIS EXPRESSED AS A HOLLERITH CONSTANT        M19900
C                                                                         M19910
C   N IS THE NUMBER OF CHARACTERS IN THE LABEL                            M19920
C                                                                         M19930
C        NEGATIVE N PLACES THE LABEL ON THE CLOCKWISE SIDE OF THE AXIS    M19940
C                                                                         M19950
C        POSITIVE N PLACES THE LABEL ON THE COUNTERCLOCKWISE SIDE         M19960
C                                                                         M19970
C   NUMCYC IS THE NUMBER OF CYCLES DESIRED                                M19980
C                                                                         M19990
C   CYCLEN IS THE LENGTH OF ONE CYCLE IN INCHES                           M20000
C                                                                         M20010
C   BEGEXP IS THE EXPONENT FOR THE BEGINNING OF THE AXIS                  M20020
C                                                                         M20030
C   THETA IS THE ANGLE OF THE AXIS IN DEGREES  (0.0 FOR X, 90.0 FOR Y)    M20040
C                                                                         M20050
C   HEIGHT IS THE HEIGHT IN INCHES OF THE TENS                            M20060
C                                                                         M20070
C   NRPT IS THE REPEAT FACTOR FOR THE SCALE NUMBERS (USUALLY INTEGER 1)   M20080
C                                                                         M20090
C        WHEN NRPT IS ZERO THE SCALE NUMBERS WILL BE SUPPRESSED;          M20100
C                                                                         M20110
C        WHEN NRPT = 2, EVERY 2ND SCALE NUMBER WILL BE PRODUCED;          M20120
C                                                                         M20130
C        WHEN NRPT = 3, EVERY 3RD SCALE NUMBER WILL BE PRODUCED; ETC.     M20140
C                                                                         M20150
C    NTURN EQUAL TO 1 TURNS THE AXIS NUMBERS BY 90 DEGREES CLOCKWISE,     M20160
C                  -1 TURNS NUMBERS BY 90 DEGREES COUNTERCLOCKWISE,       M20170
C                   0 FOR NO TURN                                         M20180
C                                                                         M20190
C   NOEND EQUAL TO 1 SUPPRESSES THE NUMBERS AT EITHER END OF THE AXIS,    M20200
C                                                                         M20210
C        NOEND EQUAL TO 2 SUPPRESSES ONLY THE STARTING NUMBER, AND        M20220
C                                                                         M20230
C        NOEND EQUAL TO 3 SUPPRESSES ONLY THE ENDING NUMBER               M20240
C                                                                         M20250
C   LSUPR EQUAL TO 1 SUPPRESSES THE LABEL                                 M20260
C                                                                         M20270
C    LTURN EQUAL TO 1 TURNS THE LABEL BY 90 DEGREES CLOCKWISE,            M20280
C         -1 TURNS LABEL BY 90 DEGREES COUNTERCLOCKWISE, 0 FOR NO TURN    M20290
C                                                                         M20300
C                                                                         M20310
      COMMON /TITLOC/ XPOS,YPOS                                           M20320
      DIMENSION BCD(*),SUBDIV(10),DIVLOG(8)                               M20330
      CHARACTER HTEN*2                                                    M20340
C                                                                         M20350
      DATA DIVLOG / 0.301029995664, 0.477121254720, 0.602059991328,       M20360
     *              0.698970004336, 0.778151250384, 0.845098040014,       M20370
     *              0.903089986992, 0.954242509439 /                      M20380
      DATA HTEN / '10'/                                                   M20390
C                                                                         M20400
      THETA1 = THETA-90.*NTURN                                            M20410
      PI = 2.*ASIN(1.)                                                    M20420
      ANGLE = (PI/180.)*THETA                                             M20430
      SINANG = SIN(ANGLE)                                                 M20440
      COSANG = COS(ANGLE)                                                 M20450
      SIGNAX = FLOAT(ISIGN(1,N))                                          M20460
      SIZMAJ = 0.25*HEIGHT+0.05                                           M20470
      OFFST = HEIGHT*1.5                                                  M20480
      DXMAJ = -SIZMAJ*SINANG*SIGNAX                                       M20490
      DYMAJ = SIZMAJ*COSANG*SIGNAX                                        M20500
      DXMIN = 0.5*DXMAJ                                                   M20510
      DYMIN = 0.5*DYMAJ                                                   M20520
      BCDSIZ = 1.25*HEIGHT                                                M20530
      ENLARG = 1.5                                                        M20540
      EXPSIZ = 0.60*HEIGHT*ENLARG                                         M20550
      NABS = IABS(N)                                                      M20560
      BCDLEN = (FLOAT(NABS)-0.4)*BCDSIZ                                   M20570
      S = CYCLEN*FLOAT(NUMCYC)                                            M20580
      NUMTIC = 2-MIN1(1.0,CYCLEN)                                         M20590
      NUMLOG = 8/NUMTIC                                                   M20600
      XBIAS = 1.85*HEIGHT*ENLARG                                          M20610
      YBIAS = 0.70*HEIGHT*ENLARG                                          M20620
      EXPBX = -YBIAS*SINANG+XBIAS*COSANG                                  M20630
      EXPBY = YBIAS*COSANG+XBIAS*SINANG                                   M20640
      IF (NTURN.EQ.0) GO TO 10                                            M20650
      EXPBX = XBIAS*SINANG+YBIAS*COSANG                                   M20660
      EXPBY = YBIAS*SINANG-XBIAS*COSANG                                   M20670
   10 DO 20 I = 2, 9                                                      M20680
         SUBDIV(I) = DIVLOG(I-1)*CYCLEN                                   M20690
   20 CONTINUE                                                            M20700
      NNUMB = NUMCYC+1                                                    M20710
      SIZMAX = EXPSIZ                                                     M20720
      EXP = BEGEXP                                                        M20730
      DO 30 I = 1, NNUMB                                                  M20740
         DIGITS = 0.0                                                     M20750
         IF (ABS(EXP).GE.10.0) DIGITS = ALOG10(ABS(EXP))                  M20760
         DIGITS = AINT(DIGITS+1.0E-12)+0.7                                M20770
         IF (EXP.LT.0.0) DIGITS = DIGITS+1.0                              M20780
         SIZNUM = DIGITS*EXPSIZ                                           M20790
         SIZMAX = MAX(SIZMAX,SIZNUM)                                      M20800
         EXP = EXP+1.0                                                    M20810
   30 CONTINUE                                                            M20820
      SIZNUM = SIZMAX+2.0*HEIGHT                                          M20830
C                                                                         M20840
C   DRAW CYCLE NUMBERS AND EXPONENTS                                      M20850
C                                                                         M20860
      NCYCLE = 0                                                          M20870
      EXP = BEGEXP                                                        M20880
      XBIAS = -0.7*HEIGHT                                                 M20890
      YBIAS = (-0.75+SIGN(2.00,SIGNAX))*HEIGHT                            M20900
      XBIAS1 = 0.                                                         M20910
      YBIAS1 = 0.                                                         M20920
      IF (NTURN.EQ.0) GO TO 40                                            M20930
      XBIAS1 = XBIAS+HEIGHT*0.5                                           M20940
      YBIAS1 = YBIAS1-OFFST-SIZNUM                                        M20950
      IF (N.LT.0) YBIAS1 = YBIAS+OFFST                                    M20960
   40 TENBX = -YBIAS*SINANG+XBIAS*COSANG+YBIAS1*SINANG-XBIAS1*COSANG      M20970
      TENBY = YBIAS*COSANG+XBIAS*SINANG-YBIAS1*COSANG-XBIAS1*SINANG       M20980
   50 XTIC = X+FLOAT(NCYCLE)*CYCLEN*COSANG                                M20990
      YTIC = Y+FLOAT(NCYCLE)*CYCLEN*SINANG                                M21000
      IF (NRPT.EQ.0) GO TO 60                                             M21010
      NSUPR = NCYCLE-(NCYCLE/NRPT)*NRPT                                   M21020
      IF (NSUPR.NE.0) GO TO 60                                            M21030
      IF ((NOEND.EQ.1.OR.NOEND.EQ.2).AND.NCYCLE.EQ.0) GO TO 60            M21040
      IF ((NOEND.EQ.1.OR.NOEND.EQ.3).AND.NCYCLE.EQ.NUMCYC) GO TO 60       M21050
      XPOS = XTIC+TENBX                                                   M21060
      YPOS = YTIC+TENBY                                                   M21070
      CALL SYMBOL (XPOS,YPOS,HEIGHT,HTEN,THETA1,2)                        M21080
      CALL NUMBER ((XPOS+EXPBX),(YPOS+EXPBY),EXPSIZ,EXP,THETA1,-1)        M21090
   60 CONTINUE                                                            M21100
C                                                                         M21110
C   DRAW TIC MARKS                                                        M21120
C                                                                         M21130
      CALL PLOT (XTIC,YTIC,3)                                             M21140
      CALL PLOT (XTIC+DXMAJ,YTIC+DYMAJ,2)                                 M21150
      IF (NCYCLE.EQ.NUMCYC) GO TO 90                                      M21160
      IF (NRPT.LT.0) GO TO 80                                             M21170
      DO 70 ILOG = 1, NUMLOG                                              M21180
         I = ILOG*NUMTIC+1/NUMTIC                                         M21190
         XLOG = XTIC+SUBDIV(I)*COSANG                                     M21200
         YLOG = YTIC+SUBDIV(I)*SINANG                                     M21210
         CALL PLOT (XLOG+DXMIN,YLOG+DYMIN,3)                              M21220
         CALL PLOT (XLOG,YLOG,2)                                          M21230
   70 CONTINUE                                                            M21240
   80 NCYCLE = NCYCLE+1                                                   M21250
      EXP = EXP+1.0                                                       M21260
      GO TO 50                                                            M21270
C                                                                         M21280
C   DRAW AXIS                                                             M21290
C                                                                         M21300
   90 CALL PLOT (XTIC,YTIC,3)                                             M21310
      CALL PLOT (X,Y,2)                                                   M21320
C                                                                         M21330
C   DRAW LABEL                                                            M21340
C                                                                         M21350
      IF (LSUPR.EQ.1.OR.NABS.EQ.0) RETURN                                 M21360
      XBIAS = 0.5*(S-BCDLEN)                                              M21370
      YBIAS = (-0.50+SIGN(3.25,SIGNAX))*BCDSIZ                            M21380
      THETA2 = THETA-90.*LTURN                                            M21390
      XBIAS2 = 0.                                                         M21400
      OFFST = HEIGHT*2.5                                                  M21410
      YBIAS2 = -SIZNUM-OFFST                                              M21420
      IF (N.LT.0) YBIAS2 = OFFST                                          M21430
      IF (LTURN.EQ.0) GO TO 100                                           M21440
      XBIAS2 = XBIAS-0.5*(S-HEIGHT)                                       M21450
  100 XPOS = X-YBIAS*SINANG+XBIAS*COSANG+YBIAS2*SINANG-XBIAS2*COSANG      M21460
      YPOS = Y+YBIAS*COSANG+XBIAS*SINANG-YBIAS2*COSANG-XBIAS2*SINANG      M21470
      CALL SYMBOL (XPOS,YPOS,BCDSIZ,BCD,THETA2,NABS)                      M21480
C                                                                         M21490
      RETURN                                                              M21500
C                                                                         M21510
      END                                                                 M21520
      SUBROUTINE LBH1H2 (YT,XT,YID)                                       M21530
C                                                                         M21540
C     LABEL H1,H2,ANGLE                                                   M21550
C                                                                         M21560
      DIMENSION YID(10)                                                   M21570
C                                                                         M21580
      DOUBLE PRECISION YID                                              & M21590
C                                                                         M21600
      CHARACTER CARD3*30,ACRD3*30                                         M21610
      CHARACTER GOUT*64                                                   M21620
C                                                                         M21630
      DATA CARD3 / '     H1,       H2,       ANGL'/                       M21640
C                                                                         M21650
      WRITE (GOUT,900) (YID(I),I=3,7)                                     M21660
      READ (GOUT,905) H1,H2,ANGLE                                         M21670
      WRITE (GOUT,910) H1,H2,ANGLE                                        M21680
      READ (GOUT,915) ACRD3                                               M21690
      YT = YT-0.2                                                         M21700
      CALL SYMBOL (XT,YT,0.12,CARD3,0.0,30)                               M21710
      YT = YT-0.2                                                         M21720
      CALL SYMBOL (XT,YT,0.12,ACRD3,0.0,30)                               M21730
C                                                                         M21740
      RETURN                                                              M21750
C                                                                         M21760
  900 FORMAT (8A8)                                                        M21770
  905 FORMAT (3F6.3)                                                      M21780
  910 FORMAT (3F10.3)                                                     M21790
  915 FORMAT (A30)                                                        M21800
C                                                                         M21810
      END                                                                 M21820
      SUBROUTINE LBAERS (YT,XT,YID)                                       M21830
C                                                                         M21840
C     LABEL AERSOL                                                        M21850
C                                                                         M21860
      DIMENSION YID(10)                                                   M21870
C                                                                         M21880
      DOUBLE PRECISION YID                                              & M21890
C                                                                         M21900
      CHARACTER CARD2*85,ACRD2*85                                         M21910
      CHARACTER GOUT*85                                                   M21920
C                                                                         M21930
      DATA CARD2(1:43)/ 'IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,   VIS,'/,   M21940
     *     CARD2(44:85)/ '      WSS,      WHH,      RAINRT,   GNDALT'/    M21950
C                                                                         M21960
      WRITE (GOUT,900) (YID(I),I=3,7)                                     M21970
      READ (GOUT,905) IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,IVIS,IWSS,      M21980
     *                IWHH,IRAINR,IGNDAL                                  M21990
      VIS = IVIS/10                                                       M22000
      WSS = IWSS/10                                                       M22010
      WHH = IWHH/10                                                       M22020
      RAINRT = IRAINR/10                                                  M22030
      GNDALT = IGNDAL/10                                                  M22040
      WRITE (GOUT,910) IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,   M22050
     *                 RAINRT,GNDALT                                      M22060
      READ (GOUT,915) ACRD2                                               M22070
      YT = YT-0.2                                                         M22080
      CALL SYMBOL (XT,YT,0.12,CARD2,0.0,85)                               M22090
      YT = YT-0.2                                                         M22100
      CALL SYMBOL (XT,YT,0.12,ACRD2,0.0,85)                               M22110
C                                                                         M22120
      RETURN                                                              M22130
C                                                                         M22140
  900 FORMAT (8A8)                                                        M22150
  905 FORMAT (18X,4I1,I2,I1,I4,3I3,I2)                                    M22160
  910 FORMAT (2I5,I6,I7,2I5,5F10.3)                                       M22170
  915 FORMAT (A85)                                                        M22180
C                                                                         M22190
      END                                                                 M22200
      SUBROUTINE YDIH1 (H1,H2,ANGLE,YID)                                  M22210
C                                                                         M22220
C     PUT H1 H2 ANGLE INTO YID                                            M22230
C                                                                         M22240
      DIMENSION YID(10)                                                   M22250
C                                                                         M22260
      DOUBLE PRECISION YID                                              & M22270
C                                                                         M22280
      CHARACTER GOUT*64,HOUT*40                                           M22290
      WRITE (HOUT,900) (YID(I),I=3,7)                                     M22300
      IF (H1.LT.100.AND.H2.LT.100.) THEN                                  M22310
         WRITE (GOUT(1:18),905) H1,H2,ANGLE                               M22320
      ELSE                                                                M22330
         WRITE (GOUT(1:18),910) H1,H2,ANGLE                               M22340
      ENDIF                                                               M22350
      GOUT(19:40) = HOUT(19:40)                                           M22360
      READ (GOUT,900) (YID(I),I=3,7)                                      M22370
C                                                                         M22380
      RETURN                                                              M22390
C                                                                         M22400
  900 FORMAT (8A8)                                                        M22410
  905 FORMAT (2F6.3,F6.2)                                                 M22420
  910 FORMAT (3F6.2)                                                      M22430
C                                                                         M22440
      END                                                                 M22450
      SUBROUTINE YDIAR (IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH   M22460
     *                 ,RAINRT,GNDALT,YID)                                M22470
C                                                                         M22480
C     IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT       M22490
C                                                                         M22500
      DIMENSION YID(10)                                                   M22510
C                                                                         M22520
      DOUBLE PRECISION YID                                              & M22530
C                                                                         M22540
      CHARACTER GOUT*64,BLNK*18                                           M22550
      DATA BLNK / '                  '/                                   M22560
      WRITE (GOUT,900) (YID(I),I=3,7)                                     M22570
      IVIS = VIS*10                                                       M22580
      IWSS = WSS*10                                                       M22590
      IWHH = WHH*10                                                       M22600
      IRAINR = RAINRT*10                                                  M22610
      IGNDAL = GNDALT*10                                                  M22620
      WRITE (GOUT(19:40),905) IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,IVIS,   M22630
     *   IWSS,IWHH,IRAINR,IGNDAL                                          M22640
      GOUT(1:18) = BLNK(1:18)                                             M22650
      READ (GOUT,900) (YID(I),I=3,7)                                      M22660
C                                                                         M22670
      RETURN                                                              M22680
C                                                                         M22690
  900 FORMAT (8A8)                                                        M22700
  905 FORMAT (4I1,I2,I1,I4,3I3,I2)                                        M22710
C                                                                         M22720
      END                                                                 M22730
