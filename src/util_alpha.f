C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
C     --------------------------------------------------------------
      SUBROUTINE BUFIN (IFILE,IEOF,IARRAY,IWORDS)
C
C     THIS SUBROUTINE BUFFERS IN (READS) IWORDS INTO  IARRAY STARTING
C     AT LOCATION IARRAY
C
C     IFILE IS THE FILE DESIGNATION
C                                                                         A10830
      COMMON /HVERSN/  HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,
     *                HVROPR,HVRPST,HVRPLT,HVRTST,HVRUTL,HVRXMR,hvnlte
C
      CHARACTER*15 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *            HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR,hvnlte
C                                  
      DIMENSION IARRAY(IWORDS)
C
C     ASSIGN CVS VERSION NUMBER TO MODULE 
C
      HVRUTL = '$Revision$' 
C                          
      IEOF = 1             
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
C
C     --------------------------------------------------------------
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c         note the name change

      SUBROUTINE BUFINln (IFILE,IEOF,IARRAY,IWORDS)
C
C     THIS SUBROUTINE BUFFERS IN (READS) IWORDS INTO  IARRAY STARTING
C     AT LOCATION IARRAY
C
C     IFILE IS THE FILE DESIGNATION
C                                  

      implicit integer*4 (i-n)
      implicit real*4    (a-h,o-z)


      DIMENSION IARRAY(IWORDS)
C                                                                         A10830
C     ASSIGN SCCS VERSION NUMBER TO MODULE 
C
      HVRUTL = '$Revision$' 
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
c______________________________________________________________________________

      SUBROUTINE BUFOUT (IFILE,IARRAY,IWORDS)
C                                                 
C     THIS SUBROUTINE BUFFERS OUT (WRITES) IWORDS FROM IARRAY STARTING
C     AT LOCATION IARRAY                                                 
C                                                                     
C     IFILE IS THE FILE DESIGNATION                                   
C                                                                     
      DIMENSION IARRAY(IWORDS)
C                                                   
      WRITE (IFILE) IARRAY
C                                                    
      RETURN                                         
C                                                    
      END                                            
C
c______________________________________________________________________________

      SUBROUTINE BUFOUTln (IFILE,IARRAY,IWORDS)
C                                                 
C     THIS SUBROUTINE BUFFERS OUT (WRITES) IWORDS FROM IARRAY STARTING
C     AT LOCATION IARRAY                                                 
C                                                                     
C     IFILE IS THE FILE DESIGNATION                                   
C                                                                     
c
      implicit integer*4 (i-n)
      implicit real*4    (a-h,o-z)
c
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
C     --------------------------------------------------------------
      SUBROUTINE LBLDAT(HDATE)                                           LN05190
C                                                                        LN05200
      EXTERNAL IDATE
C
      DOUBLE PRECISION HDATE                                            &LN05210
C                                                                        LN05220
      CHARACTER GDATE*10                                                 LN05230
C                                                                        LN05240
      INTEGER*4 IARRAY(3)                                               >LN05250
C                                                                        LN05260
      CALL IDATE(IARRAY)                                                >LN05280
      IARRAY(3)=MOD(IARRAY(3),100)                                      >LN05290
      WRITE (GDATE,900) IARRAY(3),IARRAY(2),IARRAY(1)                   >LN05300
C                                                                        LN05310
      READ (GDATE,901) HDATE                                             LN05320
C                                                                        LN05400
      RETURN                                                             LN05410
C                                                                        LN05420
  900 FORMAT (I2,2('/',I2.2))                                            LN05430
  901 FORMAT (A8)                                                        LN05450
C                                                                        LN05460
      END                                                                LN05470
C
C     --------------------------------------------------------------
      SUBROUTINE FTIME (HTIME)                                           LN05480
C                                                                        LN05490
      DOUBLE PRECISION HTIME                                            &LN05500
C                                                                        LN05510
      CHARACTER GTIME*10                                                 LN05520
C                                                                        LN05530
      INTEGER*4 IARRAY(5)                                               >LN05540
C                                                                        LN05550
      CALL ITIME (IARRAY)                                               >LN05580
      WRITE (GTIME,900) IARRAY(1),IARRAY(3),IARRAY(5)                   >LN05590
C                                                                        LN05600
      READ (GTIME,901) HTIME                                             LN05610
C                                                                        LN05620
      RETURN                                                             LN05700
C                                                                        LN05710
  900 FORMAT (I2,2(':',I2.2))                                            LN05720
  901 FORMAT (A8)                                                        LN05740
C                                                                        LN05750
      END                                                                LN05760
C
C     --------------------------------------------------------------
      SUBROUTINE CPUTIM (TIME)                                           LN05770
C                                                                        LN05780
      COMMON /TIMIN/ A1                                                  LN05790
C                                                                        LN05800
      REAL*4 ETIME,TARRAY(2)                                            >LN05810
C                                                                        LN05820
C     THIS SUBROUTINE OBTAINS CPU TIME                                   LN05830
C                                                                        LN05840
      IF (A1.LE.0.) THEN                                                 LN05850
         TIME = ETIME(TARRAY)                                           >LN05880
      ELSE                                                               LN05890
         TIME = ETIME(TARRAY)                                           >LN05920
      ENDIF                                                              LN05930
C                                                                        LN05940
      RETURN                                                             LN05950
C                                                                        LN05960
      END                                                                LN05970
C
C     --------------------------------------------------------------
      BLOCK DATA BTIM                                                    LN05980
C                                                                        LN05990
      COMMON /TIMIN/ A1                                                  LN06000
C                                                                        LN06010
      DATA A1 / 0. /                                                     LN06020
C                                                                        LN06030
      END                                                                LN06040

