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
SUBROUTINE PLTID3 (PROGID,XSIZ,YSIZ,FAC)
!
!----------------------------------------------------------------------
!
!     THE PURPOSE OF THIS SUBROUTINE PACKAGE IS TO PROVIDE
!     AN INTERFACE FROM THE CALCOMP PEN ROUTINES CURRENTLY
!     ON THE AFGL CYBER TO THE NCAR GKS PLOTTING PACKAGE.
!
!     NOTE: THESE ROUTINES WERE WRITTEN TO INTERFACE THE
!           CALCOMP CALLS USED BY LBLRTM, AND DO NOT
!           PROVIDE A COMPLETE CALCOMP TO NCAR INTERFACE.
!
!           (COMPATIBLE WITH VERSION 2.00 OF NCAR GKS)
!
!                                           A.E.R. (DECEMBER 1989)
!
!----------------------------------------------------------------------
!
!     SUBROUTINE PLTID3(PROGID,XMAX,YMAX,FACTOR)    // AFGL //
!
!     SUBROUTINE PLTID3 IS USED IN PEN MODE ONLY.  THIS SUBROUTINE
!     MUST BE THE FIRST ROUTINE CALLED AS IT INITIALIZES THE PLOT.
!
!     PROGID = AN ARRAY OF ALPHANUMERIC CHARACTERS (HOLLERITH)
!              DIMENSIONED BY THREE, TO BE USED AS IDENTIFICATION
!              FOR THE PLOT.
!     XMAX   = MAXIMUM LENGTH OF X (REAL) IN INCHES.
!     YMAX   = MAXIMUM LENGTH OF Y (REAL) IN INCHES.
!     FACTOR = A MULTIPLICATIVE FACTOR TO CHANGE SIZE OT PLOTTING.
!
!----------------------------------------------------------------------
!
   IMPLICIT REAL*8           (V)
!
   COMMON /CVRNCG/ HNAMNCG,HVRNCG
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NODUM
!
   COMMON /NCARID/ IFIRST
!
   CHARACTER*18 HNAMNCG,HVRNCG
!
   LOGICAL IFIRST
   DIMENSION PROGID(*)
!
!     ASSIGN CVS VERSION NUMBER TO MODULE
!
   HVRNCG = '$Revision$'
!
!     THIS ROUTINE DOES THE INITIALIZATION FOR THE NCAR ROUTINES.
!
!     NOTE: ALL FOUR PASSED QUANTITIES ARE IGNORED.
!
   IF (IFIRST) THEN
      CALL GOPKS (6,IDUM)
      CALL GOPWK (1,2,1)
      CALL GACWK (1)
      CALL GSCLIP (0)
      IFIRST = .FALSE.
   ELSE
      CALL FRAME
   ENDIF
   IF (JHDR.EQ.0) CALL SET (0.1,0.9,0.1,0.65,0.,10.,0.,10.,1)
!
   RETURN
!
end subroutine PLTID3
SUBROUTINE ENDPLT
!
!----------------------------------------------------------------------
!
!     SUBROUTINE ENDPLT   // AFGL //
!
!     SUBROUTINE ENDPLT MUST BE THE LAST SUBROUTINE
!     TO BE CALLED IN ALL PLOTTING JOBS.
!
!----------------------------------------------------------------------
!
!     THIS SUBROUTINE CALLS THE NCAR CLOSING ROUTINES
!
!----------------------------------------------------------------------
!
   CALL GDAWK (1)
   CALL GCLWK (1)
   CALL GCLKS
!
   RETURN
!
end subroutine ENDPLT
SUBROUTINE PLOT (X,Y,IPEN)
!
   IMPLICIT REAL*8           (V)
!
!----------------------------------------------------------------------
!
!     SUBROUTINE PLOT(X,Y,IC)   // AFGL //
!
!     SUBROUTINE PLOT IS USED TO MOVE THE PEN
!     AND TO REDEFINE A NEW ORIGIN.
!
!     X   = X COORDINATE, IN INCHES. (REAL)
!     Y   = Y COORDINATE, IN INCHES. (REAL)
!     IC  = IF IC=2, PEN DOWN AS PEN MOVES TO (X,Y).  (INTEGER)
!           IF IC=3, PEN UP AS PEN MOVES TO (X,Y).
!           IF IC=-2 OR -3, A NEW ORIGIN IS DEFINED AT (X,Y).
!
!----------------------------------------------------------------------
!
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NODUM
!
   FAC = 0.141176471
   IF (IPEN.LT.0) THEN
      IF (X.NE.1.) GO TO 10
      IF (JHDR.EQ.0) CALL FRAME
      FL = 0.175
      FR = XSIZE*FAC+FL
      FB = 0.175
      FT = YSIZE*FAC+FB
!
!     IF FR > .85 ADJUST SCALE ON BOTH X AND Y TO RESET
!
      IF (FR.GT.0.85) THEN
         SX = FR/0.85
         FR = 0.85
         FT = FT/SX
      ENDIF
!
!     IF FT > .65 ADJUST SCALE ON BOTH X AND Y TO RESET
!
      IF (FT.GT.0.65) THEN
         SY = FT/0.65
         FT = 0.65
         FR = FR/SY
      ENDIF
!
      UL = 0.
      UR = XSIZE
      UB = 0.
      UT = YSIZE
      L = 1
      CALL SET (FL,FR,FB,FT,UL,UR,UB,UT,L)
   ELSE
      RX = CUFX(X)
      RY = CUFY(Y)
      IF (IPEN.EQ.2) THEN
         CALL PLOTIF (RX,RY,1)
      ELSE
         CALL PLOTIF (RX,RY,0)
      ENDIF
   ENDIF
!
10 RETURN
!
end subroutine PLOT
SUBROUTINE LINE (X,Y,N,K,J)
!
   IMPLICIT REAL*8           (V)
!
!----------------------------------------------------------------------
!
!     SUBROUTINE LINE(X,Y,N,K,J,L)   // AFGL //
!
!     SUBROUTINE LINE PRODUCES A SINGLE LINE BY CONNECTING THE POINTS
!     DEFINED IN THE DIMENSION VARIABLES X AND Y.
!
!     X   = ARRAY OF X COORDINATES. (REAL)
!     Y   = ARRAY OF Y COORDINATES. (REAL)
!     N   = ACTUAL NUMBER OF POINTS TO BE PLOTTED.  (INTEGER)
!     K   = REPEAT CYCLE. (INTEGER, USUALLY K=1)
!     J   = CONTROL FOR USING SYMBOLS. (INTEGER)
!           J = 0 WILL PRODUCE A LINE PLOT WITHOUT SYMBOLS.
!           J > 0 WILL PRODUCE A LINE PLOT WITH A SYMBOL
!                 AT EVERY JTH POINT.
!           J < 0 WILL SUPPRESS THE LINE BETWEEN POINTS.
!     L   = A NUMBER DESCRIBING THE SYMBOL TO BE USED.
!
!----------------------------------------------------------------------
!
   COMMON /AXISXY/ V1,V2,XSIZE,YMIN,YMAX,YSIZE,IDEC,JEMIT,JPLOT,     &
   &                LOGPLT,NUMDVX,NUMSBX,DIVLNX,DELV,NUMDVY,NUMSBY,   &
   &                DIVLNY,DELY,HGT,YPL,DX,DY,NOENDX,NOENDY,IXDEC,    &
   &                JOUT,JPLTFL,JHDR,IFUNCT,NODUM
!
   DIMENSION X(*),Y(*)
!
   CALL GETSET (FL,FR,FB,FT,UL,UR,UB,UT,LL)
   UL = V1
   UR = V2
   UB = YMIN
   UT = YMAX
   CALL SET (FL,FR,FB,FT,UL,UR,UB,UT,LL)
   CALL CURVE (X,Y,N)
!
   RETURN
!
end subroutine LINE
SUBROUTINE SYMBOL (X,Y,H,BCD,T,N)
!
!----------------------------------------------------------------------
!
!     SUBROUTINE SYMBOL (X,Y,H,BCD,T,N)  // AFGL //
!
!     SUBROUTINE SYMBOL WILL DRAW A SERIES OF SYMBOLS.
!
!     X   = X COORDINATE OF THE LOWER LEFT CORNER OF THE FIRST CHARACTE
!           IN INCHES, RELATIVE TO THE CURRENT DEFINED ORIGIN.  (REAL)
!     Y   = Y COORDINATE OF THE LOWER LEFT CORNER OF THE FIRST CHARACTE
!           IN INCHES, RELATIVE TO THE CURRENT DEFINED ORIGIN.  (REAL)
!     H   = THE HEIGHT PF THE CHARACTERS, IN INCHES.  (REAL)
!     BCD = THIS PARAMETER AND THE PARAMETER N DETERMINE THE TYPE OF
!           ANNOTATION THE ROUTINE PRODUCES.  BCD CONTAINS EITHER A
!           STRING OF CHARACTERS OR THE INTEGER EQUIVALENT OF A
!           DESIRED SYMBOL.  (SEE N BELOW)
!     T   = THE ANGULAR ORIENTATION WITH RESPECT TO THE X AXIS,
!           COUNTER-CLOCKWISE IN DEGREES.  (REAL)
!     N   = THIS PARAMETER AND PARAMETER BCD DETERMINE TYPE OF
!           LETTERING OR SYMBOLS PRODUCED BY ROUTINE SYMBOL.
!           N > 0 - DEFINES CHARACTER COUNT IN BCD, LEFT JUSTIFIED.
!           N = 0 - DEFINES SINGLE CHARACTER TO BE PLOTTED,
!                                                  RIGHT JUSTIFIED.
!           N < 0 - DEFINES BCD TO BE THE INTEGER EQUIVALENT OF A
!                   SYMBOL.  FOR N = -1, THE PEN IS UP DURING THE
!                   MOVE, FOR N = -2, THE PEN IS DOWN DURING THE
!                   MOVE, AFTER WHICH A SYMBOL IS PRODUCED.
!
!               (NOTE: BCD = '3' CORRESPONDS TO A '+' SYMBOL.)
!
!----------------------------------------------------------------------
!
   CHARACTER BCD*(*),PLUS*1,COLON*1,POINT*1,BCDC(43)*1,BCDA*43
!
   EQUIVALENCE (BCDC(1),BCDA)
!
   DATA PLUS / '+'/,COLON / ':'/,POINT / '.'/
!
   CALL PWRITX (1.,1.,'''KRU''',5,1,0,0)
   IS = KUPX(H)/11
   IF (H.EQ.0.1) IS = (IS*3)/5
   IF (N.EQ.43) THEN
      READ (BCD,'(A43)') BCDA
      DO 10 JJ = 1, 43
         IF (BCDC(JJ).EQ.COLON) BCDC(JJ) = POINT
10    CONTINUE
      READ (BCDA,'(A43)') BCD
   ENDIF
   IO = T
   IF (N.GT.0) THEN
      CALL PWRITX (X,Y,BCD,N,IS,IO,-1)
   ELSE
      CALL PWRITX (X,Y,PLUS,1,IS,IO,0)
   ENDIF
!
   RETURN
!
end subroutine SYMBOL
!
SUBROUTINE NUMBER (X,Y,H,F,T,N)
!
!----------------------------------------------------------------------
!
!     SUBROUTINE NUMBER(X,Y,HGHT,FPN,THETA,N)  // AFGL //
!
!     SUBROUTINE NUMBER WILL INTERPRET AND PLOT A REAL ( REALING POINT)
!     OR INTEGER NUMBER
!
!     X     = X COORDINATE OF LOWER LEFT HAND CORNER OF THE HIGH ORDER
!             DIGIT, IN INCHES, RELATIVE TO THE CURRENT ORIGIN.  (REAL)
!     Y     = Y COORDINATE OF LOWER LEFT HAND CORNER OF THE HIGH ORDER
!             DIGIT, IN INCHES, RELATIVE TO THE CURRENT ORIGIN.  (REAL)
!     HGHT  = HEIGHT OF NUMBERS TO BE PLOTTED, IN INCHES.  (REAL)
!     FPN   = NUMBER TO BE PLOTTED.  (REAL)
!     THETA = ORIENTATION OF THE NUMBER WITH RESPECT TO THE X AXIS,
!             COUNTER-CLOCKWISE IN DEGREES.  (REAL)
!     N     = NUMBER OF DIGITS AFTER THE DECIMAL POINT.  (INTEGER)
!             N = -1 WILL SUPPRESS THE DECIMAL POINT.
!
!----------------------------------------------------------------------
!
   CHARACTER CHAR(30)*1,CHAR30*30,BLNK*1,POINT*1
!
   EQUIVALENCE (CHAR(1),CHAR30)
!
   LOGICAL IFIRST
!
   DATA BLNK / ' '/,POINT / '.'/
!
   IFIRST = .TRUE.
!
   DO 10 , I = 1, 30
      CHAR(I) = BLNK
10 END DO
!
   CALL PWRITX (1.,1.,'''KRU''',5,1,0,0)
!
   IS = KUPX(H)/12
   IO = T
   IF (N.EQ.-1) THEN
      INUM = NINT(F)
      WRITE (CHAR30,'(I30)') INUM
   ELSE
      INUM = NINT(F*(10.**N))
      IF (INUM.GE.0) THEN
         RNUM = ( REAL(INUM)+0.001)/(10.**N)
      ELSE
         RNUM = ( REAL(INUM)-0.001)/(10.**N)
      ENDIF
      WRITE (CHAR30,'(F30.15)') RNUM
   ENDIF
   IBEG = 0
   ICOUNT = 0
   DO 20 , I = 1, 30
      IF (CHAR(I).EQ.BLNK) THEN
         GO TO 20
      ELSE
         ICOUNT = ICOUNT+1
         IF (IFIRST) THEN
            IBEG = I
            IFIRST = .FALSE.
         ENDIF
         IF (CHAR(I).EQ.POINT) THEN
            ICOUNT = ICOUNT+N
            GO TO 30
         ENDIF
      ENDIF
20 END DO
30 IEND = IBEG+ICOUNT-1
   CALL PWRITX (X,Y,CHAR30(IBEG:IEND),ICOUNT,IS,IO,-1)
!
   RETURN
!
end subroutine NUMBER
!
! ****************************************************
BLOCK DATA NCARGKS
!
   LOGICAL IFIRST
   COMMON /NCARID/ IFIRST
!
   DATA IFIRST / .TRUE. /
!
end block data NCARGKS
