C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
C
      SUBROUTINE PLTID3 (PROGID,XSIZ,YSIZ,FAC)                            NC0010
C                                                                         NC0020
C----------------------------------------------------------------------   NC0030
C                                                                         NC0040
C     THE PURPOSE OF THIS SUBROUTINE PACKAGE IS TO PROVIDE                NC0050
C     AN INTERFACE FROM THE CALCOMP PEN ROUTINES CURRENTLY                NC0060
C     ON THE AFGL CYBER TO THE NCAR GKS PLOTTING PACKAGE.                 NC0070
C                                                                         NC0080
C     NOTE: THESE ROUTINES WERE WRITTEN TO INTERFACE THE                  NC0090
C           CALCOMP CALLS USED BY LBLRTM, AND DO NOT                      NC0100
C           PROVIDE A COMPLETE CALCOMP TO NCAR INTERFACE.                 NC0110
C                                                                         NC0120
C           (COMPATIBLE WITH VERSION 2.00 OF NCAR GKS)                    NC0130
C                                                                         NC0140
C                                           A.E.R. (DECEMBER 1989)        NC0150
C                                                                         NC0160
C----------------------------------------------------------------------   NC0170
C                                                                         NC0180
C     SUBROUTINE PLTID3(PROGID,XMAX,YMAX,FACTOR)    // AFGL //            NC0190
C                                                                         NC0200
C     SUBROUTINE PLTID3 IS USED IN PEN MODE ONLY.  THIS SUBROUTINE        NC0210
C     MUST BE THE FIRST ROUTINE CALLED AS IT INITIALIZES THE PLOT.        NC0220
C                                                                         NC0230
C     PROGID = AN ARRAY OF ALPHANUMERIC CHARACTERS (HOLLERITH)            NC0240
C              DIMENSIONED BY THREE, TO BE USED AS IDENTIFICATION         NC0250
C              FOR THE PLOT.                                              NC0260
C     XMAX   = MAXIMUM LENGTH OF X (REAL) IN INCHES.                      NC0270
C     YMAX   = MAXIMUM LENGTH OF Y (REAL) IN INCHES.                      NC0280
C     FACTOR = A MULTIPLICATIVE FACTOR TO CHANGE SIZE OT PLOTTING.        NC0290
C                                                                         NC0300
C----------------------------------------------------------------------   NC0310
C                                                                         NC0320
      IMPLICIT REAL*8           (V)                                     [ NC0330
C                                                                         NC0340
      COMMON /CVRNCG/ HVRNCG
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       NC0350
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,     NC0360
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      NC0370
     *                JOUT,JPLTFL,JHDR,IFUNCT,NODUM                       NC0380
C                                                                         NC0390
      COMMON /NCARID/ IFIRST                                              NC0400
C
      CHARACTER*15 HVRNCG
C
      LOGICAL IFIRST                                                      NC0410
      DIMENSION PROGID(*)                                                 NC0420
C                                                                         NC0430
C     ASSIGN CVS VERSION NUMBER TO MODULE 
C
      HVRNCG = '$Revision$' 
C                                                                         NC0450
C     THIS ROUTINE DOES THE INITIALIZATION FOR THE NCAR ROUTINES.         NC0460
C                                                                         NC0470
C     NOTE: ALL FOUR PASSED QUANTITIES ARE IGNORED.                       NC0480
C                                                                         NC0490
      IF (IFIRST) THEN                                                    NC0500
         CALL GOPKS (6,IDUM)                                              NC0510
         CALL GOPWK (1,2,1)                                               NC0520
         CALL GACWK (1)                                                   NC0530
         CALL GSCLIP (0)                                                  NC0540
         IFIRST = .FALSE.                                                 NC0550
      ELSE                                                                NC0560
         CALL FRAME                                                       NC0570
      ENDIF                                                               NC0580
      IF (JHDR.EQ.0) CALL SET (0.1,0.9,0.1,0.65,0.,10.,0.,10.,1)          NC0590
C                                                                         NC0600
      RETURN                                                              NC0610
C                                                                         NC0620
      END                                                                 NC0630
      SUBROUTINE ENDPLT                                                   NC0640
C                                                                         NC0650
C----------------------------------------------------------------------   NC0660
C                                                                         NC0670
C     SUBROUTINE ENDPLT   // AFGL //                                      NC0680
C                                                                         NC0690
C     SUBROUTINE ENDPLT MUST BE THE LAST SUBROUTINE                       NC0700
C     TO BE CALLED IN ALL PLOTTING JOBS.                                  NC0710
C                                                                         NC0720
C----------------------------------------------------------------------   NC0730
C                                                                         NC0740
C     THIS SUBROUTINE CALLS THE NCAR CLOSING ROUTINES                     NC0750
C                                                                         NC0760
C----------------------------------------------------------------------   NC0770
C                                                                         NC0780
      CALL GDAWK (1)                                                      NC0790
      CALL GCLWK (1)                                                      NC0800
      CALL GCLKS                                                          NC0810
C                                                                         NC0820
      RETURN                                                              NC0830
C                                                                         NC0840
      END                                                                 NC0850
      SUBROUTINE PLOT (X,Y,IPEN)                                          NC0860
C                                                                         NC0870
      IMPLICIT REAL*8           (V)                                     [ NC0880
C                                                                         NC0890
C----------------------------------------------------------------------   NC0900
C                                                                         NC0910
C     SUBROUTINE PLOT(X,Y,IC)   // AFGL //                                NC0920
C                                                                         NC0930
C     SUBROUTINE PLOT IS USED TO MOVE THE PEN                             NC0940
C     AND TO REDEFINE A NEW ORIGIN.                                       NC0950
C                                                                         NC0960
C     X   = X COORDINATE, IN INCHES. (REAL)                               NC0970
C     Y   = Y COORDINATE, IN INCHES. (REAL)                               NC0980
C     IC  = IF IC=2, PEN DOWN AS PEN MOVES TO (X,Y).  (INTEGER)           NC0990
C           IF IC=3, PEN UP AS PEN MOVES TO (X,Y).                        NC1000
C           IF IC=-2 OR -3, A NEW ORIGIN IS DEFINED AT (X,Y).             NC1010
C                                                                         NC1020
C----------------------------------------------------------------------   NC1030
C                                                                         NC1040
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       NC1050
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,     NC1060
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      NC1070
     *                JOUT,JPLTFL,JHDR,IFUNCT,NODUM                       NC1080
C                                                                         NC1090
      FAC = 0.141176471                                                   NC1100
      IF (IPEN.LT.0) THEN                                                 NC1110
         IF (X.NE.1.) GO TO 10                                            NC1120
         IF (JHDR.EQ.0) CALL FRAME                                        NC1130
         FL = 0.175                                                       NC1140
         FR = XSIZE*FAC+FL                                                NC1150
         FB = 0.175                                                       NC1160
         FT = YSIZE*FAC+FB                                                NC1170
C                                                                         NC1180
C     IF FR > .85 ADJUST SCALE ON BOTH X AND Y TO RESET                   NC1190
C                                                                         NC1200
         IF (FR.GT.0.85) THEN                                             NC1210
            SX = FR/0.85                                                  NC1220
            FR = 0.85                                                     NC1230
            FT = FT/SX                                                    NC1240
         ENDIF                                                            NC1250
C                                                                         NC1260
C     IF FT > .65 ADJUST SCALE ON BOTH X AND Y TO RESET                   NC1270
C                                                                         NC1280
         IF (FT.GT.0.65) THEN                                             NC1290
            SY = FT/0.65                                                  NC1300
            FT = 0.65                                                     NC1310
            FR = FR/SY                                                    NC1320
         ENDIF                                                            NC1330
C                                                                         NC1340
         UL = 0.                                                          NC1350
         UR = XSIZE                                                       NC1360
         UB = 0.                                                          NC1370
         UT = YSIZE                                                       NC1380
         L = 1                                                            NC1390
         CALL SET (FL,FR,FB,FT,UL,UR,UB,UT,L)                             NC1400
      ELSE                                                                NC1410
         RX = CUFX(X)                                                     NC1420
         RY = CUFY(Y)                                                     NC1430
         IF (IPEN.EQ.2) THEN                                              NC1440
            CALL PLOTIF (RX,RY,1)                                         NC1450
         ELSE                                                             NC1460
            CALL PLOTIF (RX,RY,0)                                         NC1470
         ENDIF                                                            NC1480
      ENDIF                                                               NC1490
C                                                                         NC1500
   10 RETURN                                                              NC1510
C                                                                         NC1520
      END                                                                 NC1530
      SUBROUTINE LINE (X,Y,N,K,J)                                         NC1540
C                                                                         NC1550
      IMPLICIT REAL*8           (V)                                     [ NC1560
C                                                                         NC1570
C----------------------------------------------------------------------   NC1580
C                                                                         NC1590
C     SUBROUTINE LINE(X,Y,N,K,J,L)   // AFGL //                           NC1600
C                                                                         NC1610
C     SUBROUTINE LINE PRODUCES A SINGLE LINE BY CONNECTING THE POINTS     NC1620
C     DEFINED IN THE DIMENSION VARIABLES X AND Y.                         NC1630
C                                                                         NC1640
C     X   = ARRAY OF X COORDINATES. (REAL)                                NC1650
C     Y   = ARRAY OF Y COORDINATES. (REAL)                                NC1660
C     N   = ACTUAL NUMBER OF POINTS TO BE PLOTTED.  (INTEGER)             NC1670
C     K   = REPEAT CYCLE. (INTEGER, USUALLY K=1)                          NC1680
C     J   = CONTROL FOR USING SYMBOLS. (INTEGER)                          NC1690
C           J = 0 WILL PRODUCE A LINE PLOT WITHOUT SYMBOLS.               NC1700
C           J > 0 WILL PRODUCE A LINE PLOT WITH A SYMBOL                  NC1710
C                 AT EVERY JTH POINT.                                     NC1720
C           J < 0 WILL SUPPRESS THE LINE BETWEEN POINTS.                  NC1730
C     L   = A NUMBER DESCRIBING THE SYMBOL TO BE USED.                    NC1740
C                                                                         NC1750
C----------------------------------------------------------------------   NC1760
C                                                                         NC1770
      COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,       NC1780
     *                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,     NC1790
     *                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,      NC1800
     *                JOUT,JPLTFL,JHDR,IFUNCT,NODUM                       NC1810
C                                                                         NC1820
      DIMENSION X(*),Y(*)                                                 NC1830
C                                                                         NC1840
      CALL GETSET (FL,FR,FB,FT,UL,UR,UB,UT,LL)                            NC1850
      UL = V1                                                             NC1860
      UR = V2                                                             NC1870
      UB = YMIN                                                           NC1880
      UT = YMAX                                                           NC1890
      CALL SET (FL,FR,FB,FT,UL,UR,UB,UT,LL)                               NC1900
      CALL CURVE (X,Y,N)                                                  NC1910
C                                                                         NC1920
      RETURN                                                              NC1930
C                                                                         NC1940
      END                                                                 NC1950
      SUBROUTINE SYMBOL (X,Y,H,BCD,T,N)                                   NC1960
C                                                                         NC1970
C----------------------------------------------------------------------   NC1980
C                                                                         NC1990
C     SUBROUTINE SYMBOL (X,Y,H,BCD,T,N)  // AFGL //                       NC2000
C                                                                         NC2010
C     SUBROUTINE SYMBOL WILL DRAW A SERIES OF SYMBOLS.                    NC2020
C                                                                         NC2030
C     X   = X COORDINATE OF THE LOWER LEFT CORNER OF THE FIRST CHARACTE   NC2040
C           IN INCHES, RELATIVE TO THE CURRENT DEFINED ORIGIN.  (REAL)    NC2050
C     Y   = Y COORDINATE OF THE LOWER LEFT CORNER OF THE FIRST CHARACTE   NC2060
C           IN INCHES, RELATIVE TO THE CURRENT DEFINED ORIGIN.  (REAL)    NC2070
C     H   = THE HEIGHT PF THE CHARACTERS, IN INCHES.  (REAL)              NC2080
C     BCD = THIS PARAMETER AND THE PARAMETER N DETERMINE THE TYPE OF      NC2090
C           ANNOTATION THE ROUTINE PRODUCES.  BCD CONTAINS EITHER A       NC2100
C           STRING OF CHARACTERS OR THE INTEGER EQUIVALENT OF A           NC2110
C           DESIRED SYMBOL.  (SEE N BELOW)                                NC2120
C     T   = THE ANGULAR ORIENTATION WITH RESPECT TO THE X AXIS,           NC2130
C           COUNTER-CLOCKWISE IN DEGREES.  (REAL)                         NC2140
C     N   = THIS PARAMETER AND PARAMETER BCD DETERMINE TYPE OF            NC2150
C           LETTERING OR SYMBOLS PRODUCED BY ROUTINE SYMBOL.              NC2160
C           N > 0 - DEFINES CHARACTER COUNT IN BCD, LEFT JUSTIFIED.       NC2170
C           N = 0 - DEFINES SINGLE CHARACTER TO BE PLOTTED,               NC2180
C                                                  RIGHT JUSTIFIED.       NC2190
C           N < 0 - DEFINES BCD TO BE THE INTEGER EQUIVALENT OF A         NC2200
C                   SYMBOL.  FOR N = -1, THE PEN IS UP DURING THE         NC2210
C                   MOVE, FOR N = -2, THE PEN IS DOWN DURING THE          NC2220
C                   MOVE, AFTER WHICH A SYMBOL IS PRODUCED.               NC2230
C                                                                         NC2240
C               (NOTE: BCD = '3' CORRESPONDS TO A '+' SYMBOL.)            NC2250
C                                                                         NC2260
C----------------------------------------------------------------------   NC2270
C                                                                         NC2280
      CHARACTER BCD*(*),PLUS*1,COLON*1,POINT*1,BCDC(43)*1,BCDA*43         NC2290
C                                                                         NC2300
      EQUIVALENCE (BCDC(1),BCDA)                                          NC2310
C                                                                         NC2320
      DATA PLUS / '+'/,COLON / ':'/,POINT / '.'/                          NC2330
C                                                                         NC2340
      CALL PWRITX (1.,1.,'''KRU''',5,1,0,0)                               NC2350
      IS = KUPX(H)/11                                                     NC2360
      IF (H.EQ.0.1) IS = (IS*3)/5                                         NC2370
      IF (N.EQ.43) THEN                                                   NC2380
         READ (BCD,'(A43)') BCDA                                          NC2390
         DO 10 JJ = 1, 43                                                 NC2400
            IF (BCDC(JJ).EQ.COLON) BCDC(JJ) = POINT                       NC2410
   10    CONTINUE                                                         NC2420
         READ (BCDA,'(A43)') BCD                                          NC2430
      ENDIF                                                               NC2440
      IO = T                                                              NC2450
      IF (N.GT.0) THEN                                                    NC2460
         CALL PWRITX (X,Y,BCD,N,IS,IO,-1)                                 NC2470
      ELSE                                                                NC2480
         CALL PWRITX (X,Y,PLUS,1,IS,IO,0)                                 NC2490
      ENDIF                                                               NC2500
C                                                                         NC2510
      RETURN                                                              NC2520
C                                                                         NC2530
      END                                                                 NC2540
C                                                                         NC2550
      SUBROUTINE NUMBER (X,Y,H,F,T,N)                                     NC2560
C                                                                         NC2570
C----------------------------------------------------------------------   NC2580
C                                                                         NC2590
C     SUBROUTINE NUMBER(X,Y,HGHT,FPN,THETA,N)  // AFGL //                 NC2600
C                                                                         NC2610
C     SUBROUTINE NUMBER WILL INTERPRET AND PLOT A REAL ( REALING POINT)   NC2620
C     OR INTEGER NUMBER                                                   NC2630
C                                                                         NC2640
C     X     = X COORDINATE OF LOWER LEFT HAND CORNER OF THE HIGH ORDER    NC2650
C             DIGIT, IN INCHES, RELATIVE TO THE CURRENT ORIGIN.  (REAL)   NC2660
C     Y     = Y COORDINATE OF LOWER LEFT HAND CORNER OF THE HIGH ORDER    NC2670
C             DIGIT, IN INCHES, RELATIVE TO THE CURRENT ORIGIN.  (REAL)   NC2680
C     HGHT  = HEIGHT OF NUMBERS TO BE PLOTTED, IN INCHES.  (REAL)         NC2690
C     FPN   = NUMBER TO BE PLOTTED.  (REAL)                               NC2700
C     THETA = ORIENTATION OF THE NUMBER WITH RESPECT TO THE X AXIS,       NC2710
C             COUNTER-CLOCKWISE IN DEGREES.  (REAL)                       NC2720
C     N     = NUMBER OF DIGITS AFTER THE DECIMAL POINT.  (INTEGER)        NC2730
C             N = -1 WILL SUPPRESS THE DECIMAL POINT.                     NC2740
C                                                                         NC2750
C----------------------------------------------------------------------   NC2760
C                                                                         NC2770
      CHARACTER CHAR(30)*1,CHAR30*30,BLNK*1,POINT*1                       NC2780
C                                                                         NC2790
      EQUIVALENCE (CHAR(1),CHAR30)                                        NC2800
C                                                                         NC2810
      LOGICAL IFIRST                                                      NC2820
C                                                                         NC2830
      DATA BLNK / ' '/,POINT / '.'/                                       NC2840
C                                                                         NC2850
      IFIRST = .TRUE.                                                     NC2860
C                                                                         NC2870
      DO 10 , I = 1, 30                                                   NC2880
         CHAR(I) = BLNK                                                   NC2890
   10 CONTINUE                                                            NC2900
C                                                                         NC2910
      CALL PWRITX (1.,1.,'''KRU''',5,1,0,0)                               NC2920
C                                                                         NC2930
      IS = KUPX(H)/12                                                     NC2940
      IO = T                                                              NC2950
      IF (N.EQ.-1) THEN                                                   NC2960
         INUM = NINT(F)                                                   NC2970
         WRITE (CHAR30,'(I30)') INUM                                      NC2980
      ELSE                                                                NC2990
         INUM = NINT(F*(10.**N))                                          NC3000
         IF (INUM.GE.0) THEN                                              NC3010
            RNUM = ( REAL(INUM)+0.001)/(10.**N)                           NC3020
         ELSE                                                             NC3030
            RNUM = ( REAL(INUM)-0.001)/(10.**N)                           NC3040
         ENDIF                                                            NC3050
         WRITE (CHAR30,'(F30.15)') RNUM                                   NC3060
      ENDIF                                                               NC3070
      IBEG = 0                                                            NC3080
      ICOUNT = 0                                                          NC3090
      DO 20 , I = 1, 30                                                   NC3100
         IF (CHAR(I).EQ.BLNK) THEN                                        NC3110
            GO TO 20                                                      NC3120
         ELSE                                                             NC3130
            ICOUNT = ICOUNT+1                                             NC3140
            IF (IFIRST) THEN                                              NC3150
               IBEG = I                                                   NC3160
               IFIRST = .FALSE.                                           NC3170
            ENDIF                                                         NC3180
            IF (CHAR(I).EQ.POINT) THEN                                    NC3190
               ICOUNT = ICOUNT+N                                          NC3200
               GO TO 30                                                   NC3210
            ENDIF                                                         NC3220
         ENDIF                                                            NC3230
   20 CONTINUE                                                            NC3240
   30 IEND = IBEG+ICOUNT-1                                                NC3250
      CALL PWRITX (X,Y,CHAR30(IBEG:IEND),ICOUNT,IS,IO,-1)                 NC3260
C                                                                         NC3270
      RETURN                                                              NC3280
C                                                                         NC3290
      END                                                                 NC3300
C
C ****************************************************
      BLOCK DATA NCARGKS
C
      LOGICAL IFIRST
      COMMON /NCARID/ IFIRST
C
      DATA IFIRST / .TRUE. /
C
      END
C ****************************************************
