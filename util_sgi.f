C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      SUBROUTINE BUFIN (IFILE,IEOF,IARRAY,IWORDS)
C
C     THIS SUBROUTINE BUFFERS IN (READS) IWORDS INTO  IARRAY STARTING
C     AT LOCATION IARRAY
C
C     IFILE IS THE FILE DESIGNATION
C                                  
      DIMENSION IARRAY(IWORDS)
C                          
      IEOF = 1             
C                          
C#    BUFFER IN (IFILE,1) (IARRAY(ILO),IARRAY(IHI))
C#    IF (UNIT(IFILE).EQ.0.) GO TO 10              
C                                               
      READ (IFILE,END=10) IARRAY
      ITEST = MIN(IWORDS,4)                 
      IF (IARRAY(ITEST).EQ.-99) IEOF = -99      
C                                               
      RETURN                                    
C                                               
   10 IEOF = 0                                  
C                                               
      RETURN                                    
C                                               
      END                                       
      SUBROUTINE BUFOUT (IFILE,IARRAY,IWORDS)
C                                                 
C     THIS SUBROUTINE BUFFERS OUT (WRITES) IWORDS FROM IARRAY STARTING
C     AT LOCATION IARRAY                                                 
C                                                                     
C     IFILE IS THE FILE DESIGNATION                                   
C                                                                     
      DIMENSION IARRAY(IWORDS)
C                                                   
C#    BUFFER OUT (IFILE,1) (IARRAY(ILO),IARRAY(IHI))
C#    IF (UNIT(IFILE).EQ.0.) STOP ' ERROR IN BUFOUT '
C                                                    
      WRITE (IFILE) IARRAY
C                                                    
      RETURN                                         
C                                                    
      END                                            
      SUBROUTINE FDATE (HDATE)                                           LN05190
C                                                                        LN05200
      DOUBLE PRECISION HDATE                                            &LN05210
C                                                                        LN05220
      CHARACTER GDATE*10                                                 LN05230
C                                                                        LN05240
      INTEGER*4 IARRAY(3)                                               >LN05250
C                                                                        LN05260
      CALL IDATE(IARRAY)                                                >LN05280
c     IARRAY(3)=MOD(IARRAY(3),100)                                      >LN05290
c     WRITE (GDATE,900) IARRAY(3),IARRAY(2),IARRAY(1)                   >LN05300
      WRITE (GDATE,900) IARRAY(3),IARRAY(1),IARRAY(2)                   >LN05300
C                                                                        LN05310
      READ (GDATE,901) HDATE                                             LN05320
C                                                                        LN05330

C     CHANGE THESE TO WORD SIZE AND FORMAT OF ROUTINE DATE               LN05390
C                                                                        LN05400
      RETURN                                                             LN05410
C                                                                        LN05420
  900 FORMAT (   I2,2('/',I2.2))                                         LN05430
  901 FORMAT (A8)                                                        LN05450
C                                                                        LN05460
      END                                                                LN05470
      SUBROUTINE FTIME (HTIME)                                           LN05480
C                                                                        LN05490
      DOUBLE PRECISION HTIME                                            &LN05500
C                                                                        LN05510
c     CHARACTER GTIME*10                                                 LN05520
      CHARACTER GTIME*8                                                  LN05520
C                                                                        LN05530
c     INTEGER*4 IARRAY(3)                                               >LN05540
C                                                                        LN05550
c     CALL ITIME (IARRAY)                                               >LN05580
      CALL  TIME ( clock)                                               >LN05580
c     WRITE (GTIME,900) IARRAY                                          >LN05590
      WRITE (GTIME,900)  clock                                          >LN05590
C                                                                        LN05600
      READ (GTIME,901) HTIME                                             LN05610
C                                                                        LN05620
C     GTIME AND FORMAT ARE FOR CYBER AND CRAY                            LN05630
C                                                                        LN05640
C       -- CYBER REQUIRES FORMAT (1X,A8)                                 LN05650
C       -- CRAY  REQUIRES FORMAT (A8)                                    LN05660
C                                                                        LN05670
C     CHANGE THESE TO WORD SIZE AND FORMAT OF ROUTINE GTIME              LN05680
C                                                                        LN05690
      RETURN                                                             LN05700
C                                                                        LN05710
C>900 FORMAT (1X,I2,2(':',I2.2))                                         LN05720
  900 FORMAT (   I2,2(':',I2.2))                                         LN05720
C>901 FORMAT (1X,A8)                                                    >LN05730
  901 FORMAT (A8)                                                        LN05740
C                                                                        LN05750
      END                                                                LN05760
      SUBROUTINE CPUTIM (TIME)                                           LN05770
C                                                                        LN05780
      COMMON /TIMIN/ A1                                                  LN05790
C                                                                        LN05800
c      REAL*4 ETIME,TARRAY(2)                                            >LN05810
C                                                                        LN05820
C     THIS SUBROUTINE OBTAINS CPU TIME                                   LN05830
C                                                                        LN05840
      IF (A1.LE.0.) THEN                                                 LN05850
C>>      CALL SECOND (TIME)                                              LN05860
C>VAX    A1 = SECNDS(0.0)                                               >LN05870
         A1 = SECNDS(0.0)                                               >LN05870
c        TIME = ETIME(TARRAY)                                           >LN05880
      ELSE                                                               LN05890
C>>      CALL SECOND (TIME)                                              LN05900
C>VAX    TIME = SECNDS(A1)                                              >LN05910
         TIME = SECNDS(A1)                                              >LN05910
c        TIME = ETIME(TARRAY)                                           >LN05920
      ENDIF                                                              LN05930
C                                                                        LN05940
      RETURN                                                             LN05950
C                                                                        LN05960
      END                                                                LN05970
      BLOCK DATA BTIM                                                    LN05980
C                                                                        LN05990
      COMMON /TIMIN/ A1                                                  LN06000
C                                                                        LN06010
      DATA A1 / 0. /                                                     LN06020
C                                                                        LN06030
      END                                                                LN06040
