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
      COMMON /HVERSN/  HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,
     *                HVROPR,HVRPST,HVRPLT,HVRTST,HVRUTL,HVRXMR
C
      CHARACTER*8 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *            HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR
C
      DIMENSION IARRAY(IWORDS)
C
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
      COMMON /HVERSN/  HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,
     *                HVROPR,HVRPST,HVRPLT,HVRTST,HVRUTL,HVRXMR
C
      CHARACTER*8 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *            HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR
C
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
C#    BUFFER OUT (IFILE,1) (IARRAY(ILO),IARRAY(IHI))
C#    IF (UNIT(IFILE).EQ.0.) STOP ' ERROR IN BUFOUT '
C                                                    
      WRITE (IFILE) IARRAY
C                                                    
      RETURN                                         
C                                                    
      END                                            
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
c______________________________________________________________________________

      SUBROUTINE LBLDAT(HDATE)                                           LN05190
C                                                                        LN05200
      DOUBLE PRECISION HDATE                                            &LN05210
C                                                                        LN05220
      CHARACTER GDATE*10                                                 LN05230
C                                                                        LN05240
      INTEGER*4 IARRAY(3)                                               >LN05250
C                                                                        LN05260
C>>   CALL DATE (GDATE)                                                  LN05270
C     CALL IDATE(IARRAY)                                                >LN05280
C     IARRAY(3)=MOD(IARRAY(3),100)                                      >LN05290
C     WRITE (GDATE,900) IARRAY(3),IARRAY(2),IARRAY(1)                   >LN05300
      CALL DATE(GDATE)
C                                                                        LN05310
      READ (GDATE,901) HDATE                                             LN05320
C                                                                        LN05330
C     GDATE AND FORMAT ARE FOR CYBER AND CRAY                            LN05340
C                                                                        LN05350
C       -- CYBER REQUIRES FORMAT (1X,A8)                                 LN05360
C       -- CRAY  REQUIRES FORMAT (A8)                                    LN05370
C                                                                        LN05380
C     CHANGE THESE TO WORD SIZE AND FORMAT OF ROUTINE DATE               LN05390
C                                                                        LN05400
      RETURN                                                             LN05410
C                                                                        LN05420
C>900 FORMAT (1X,I2,2('/',I2.2))                                         LN05430
  900 FORMAT (   I2,2('/',I2.2))                                         LN05430
C>901 FORMAT (1X,A8)                                                    >LN05440
  901 FORMAT (A8)                                                        LN05450
C                                                                        LN05460
      END                                                                LN05470
c******
      SUBROUTINE FTIME (HTIME)                                           LN05480
C                                                                        LN05490
      DOUBLE PRECISION HTIME                                            &LN05500
C                                                                        LN05510
      CHARACTER GTIME*10                                                 LN05520
C                                                                        LN05530
      INTEGER*4 IARRAY(3)                                               >LN05540
C                                                                        LN05550
C>>   CALL CLOCK (GTIME)                                                 LN05560
C>VAX CALL TIME (GTIME)                                                 >LN05570
C     CALL ITIME (IARRAY)                                               >LN05580
C     WRITE (GTIME,900) IARRAY                                          >LN05590
      CALL CLOCKTIME(GTIME)
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
      INTEGER*4 CPU_TIME
C                                                                        LN05800
C     REAL*4 ETIME,TARRAY(2)                                            >LN05810
C                                                                        LN05820
C     THIS SUBROUTINE OBTAINS CPU TIME                                   LN05830
C                                                                        LN05840
      IF (A1.LE.0.) THEN                                                 LN05850
C>>      CALL SECOND (TIME)                                              LN05860
C>VAX    A1 = SECNDS(0.0)                                               >LN05870
C        TIME = ETIME(TARRAY)                                           >LN05880
         CPU_TIME = MCLOCK()
         TIME = FLOAT(CPU_TIME)/100
      ELSE                                                               LN05890
C>>      CALL SECOND (TIME)                                              LN05900
C>VAX    TIME = SECNDS(A1)                                              >LN05910
C        TIME = ETIME(TARRAY)                                           >LN05920
         CPU_TIME = MCLOCK()
         TIME = FLOAT(CPU_TIME)/100
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

      subroutine date(gdate)
c=======================================================================
c Routine to use GETDAY to obtain system date and format it to MM/DD/YY

c CREATED   : 16-MAR-1992 PvD

c ARGUMENTS : GDATE - character variable in which the date is returned
c=======================================================================

      integer*4 jyyddd , jyear
      integer*2 iyear , imonth , iday
      integer*2 month(12)
      character gdate*10

      data month/31,28,31,30,31,30,31,31,30,31,30,31/

      call GETDAY(jyyddd)

c -- Extract year and julian day from packed integer
      iyear = jyyddd/1000
      jyear = iyear
      ijday = mod(jyyddd,1000)
      
c -- Identify a leap year
      if( mod(jyear,4).eq.0 )then
c       LEAP YEAR!
        month(2) = month(2) + 1
      endif

c -- Convert from julian days to calendar days
      imcnt = 0
      julday = 0
10    imcnt = imcnt + 1
      julday = julday + month(imcnt)
      if( julday.ge.ijday )then
        imonth = imcnt
        iday = ijday-julday+month(imcnt)
      else
        goto 10
      endif

c -- Write date in MM/DD/YY format to GDATE
      write(gdate,'(1x,i2.2,2(''/'',i2.2))')imonth,iday,iyear

      return
      end

      subroutine clocktime(gtime)
c=======================================================================
c Routine to use GETTIM to obtain system time and format it to HH:MM:SS 

c CREATED   : 16-MAR-1992 PvD

c ARGUMENTS : GTIME - character variable in which the time is returned
c=======================================================================

      integer*4 jhms
      integer*2 ihour , imin , isec
      character gtime*10

      call GETTIM(jhms)

c -- Extract hours, minutes and seconds from packed integer
      ihour = jhms/10000
      imin  = mod(jhms/100,100)
      isec  = mod(jhms,100)

c -- Write time in HH:MM:SS format to GTIME
      write(gtime,'(1x,i2.2,2('':'',i2.2))')ihour,imin,isec

      return
      end
