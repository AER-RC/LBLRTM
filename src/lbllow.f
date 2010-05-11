C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      SUBROUTINE LOWTRN                                                  FL00000
C
C  --------------------------------------------------------------------------
C |                                                                          |
C |  Copyright 2002 - 2009, Atmospheric & Environmental Research, Inc. (AER).|
C |  This software may be used, copied, or redistributed as long as it is    |
C |  not sold and this copyright notice is reproduced on each copy made.     |
C |  This model is provided as is without any express or implied warranties. |
C |                       (http://www.rtweb.aer.com/)                        |
C |                                                                          |
C  --------------------------------------------------------------------------
C
C                                                                        FL00010
C     CC                                                                 FL00020
C     CC   STRIPPED DOWN VERSION OF LOWTRAN 7 TO RUN AS A SUBROUTINE     FL00030
C     CC   TO SUPPLY LBLRTM WITH AEROSOLS,CLOUDS,FOGS AND RAIN           FL00040
C     CC                                                                 FL00050
c_______________________________________________________________________
c
c     1 May 2003
c
c     The interface between LBLRTM and LOWTRN has been substantially
c     modified by S. A. Clough and M. W. Shephard.
c
c     In general, the overall strategy has been changed to have LOWTRN
c     provide aerosol, cloud and rain properties on the LBLRTM 
c     vertical grid.
c
c     The implications of this are that the any' user suppled' LOWTRN
c     information should be provided on a vertical grid consistent with 
c     that of LBLRTM.
c
c     The results have been checked for a number of options: downwelling 
c     radiance, upwelling radiance, both at arbitrary zenith angles.
c
c     A known error is that the TANGENT case is NOT working properly.
c
c     Continued improvements of the interface will be implemented over 
c     time.
C                                                                        FL00060
C     ****************************************************************** FL00070
C     THIS SUBROUTINE IS ONLY USED FOR AEROSOLS AND CLOUDS               FL00080
C                                                                        FL00090
C     BUILT IN CLOUD AND RAIN MODELS ARE CHOSEN BY ICLD (RECORD 3.1)     FL00100
C                                                                        FL00110
C     USER DEFINED MODEL CAN BE INPUT BY SETTING IAERSL=7 (RECORD 1.2)   FL00120
C                                                                        FL00130
C     FOR A MORE COMPLETE EXPLANATION SEE LBLRTM USER INSTRUCTIONS       FL00140
C                                                                        FL00150
C     ****************************************************************** FL00160
C     PROGRAM ACTIVATED  BY IAERSL = 1 OR 7  (RECORD 1.2)                FL00170
C     RECORD SEQUENCE AS FOLLOWS                                         FL00180
C                                                                        FL00190
C                                                                        FL00200
C     RECORD 3.1   IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINR FL00210
C     GNDALT                                                             FL00220
C     FORMAT(6I5,5F10.3)                                                 FL00230
C                                                                        FL00240
C     RECORD 3.2  CTHIK,CALT,CEXT,ISEED   (ICLD=18,19, OR 20)            FL00250
C     FORMAT(3F10.3,I10)                                                 FL00260
C                                                                        FL00270
C     RECORD 3.3  ZCVSA,ZTVSA,ZINVSA      (IVSA=1)                       FL00280
C     FORMAT(3F10.3)                                                     FL00290
C                                                                        FL00300
C     RECORD 3.4  ML,TITLE                (IAERSL=7)                     FL00310
C     FORMAT(I5,18A4)                                                    FL00320
C                                                                        FL00330
C     RECORD 3.5 IS REPEATED ML TIMES                                    FL00340
C                                                                        FL00350
C     RECORD 3.5   ZMDL,AHAZE,EQLWCZ,RRATZ,IHAZ1,                         FL00360
C     ICLD1,IVUL1,ISEA1,ICHR1                                            FL00370
C     FORMAT (4F10.3,5I5)                                                FL00380
C                                                                        FL00390
C     RECORDS 3.6.1 - 3.6.3 READ IN THE USER DEFINED CLOUD EXTINCTION    FL00400
C     AND ABSORPTION        (IHAZE=7 OR ICLD=11)                         FL00410
C                                                                        FL00420
C     RECORD 3.6.1   (IREG(I),I=1,4)                                     FL00430
C     FORMAT (4I5)                                                       FL00440
C                                                                        FL00450
C     RECORD 3.6.2   AWCCON(N),TITLE(N)                                  FL00460
C     FORMAT (E10.3,18A4)                                                FL00470
C                                                                        FL00480
C     RECORD 3.6.3 (VX(I),EXTC(N,I),ABSC(N,I),ASYM(N,I),I=1,47)          FL00490
C     FORMAT (3(F6.2,2F7.5,F6.4)) ** 16 RECORDS **                       FL00500
C                                                                        FL00510
C                                                                        FL00520
C     ****************************************************************** FL00530
C                                                                        FL00540
C     MODEL IS READ IN LBLATM                                            FL00550
C     M1, M2, AND M3 ARE SET DEPENDING ON MODEL                          FL00560
C                                                                        FL00570
C     ****************************************************************** FL00580
C                                                                        FL00590
C     RECORD 3.1    IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,     FL00600
C     RAINRT,GNDALT                                                      FL00610
C                                                                        FL00620
C     FORMAT(6I5,5F10.3)                                                 FL00630
C                                                                        FL00640
C     IHAZE SELECTS THE TYPE OF EXTINCTION AND A DEFAULT                 FL00650
C     METEROLOGIACL RANGE FOR THE BOUNDARY-LAYER AEROSOL MODEL           FL00660
C     (0 TO 2KM ALTITUDE)                                                FL00670
C     IF VIS IS ALSO SPECIFIED ON RECORD 3.1, IT WILL OVERRIDE           FL00680
C     THE DEFAULT IHAZE VALUE                                            FL00690
C                                                                        FL00700
C     IHAZE=0  NO AEROSOL ATTENUATION INCLUDED IN CALCULATION.           FL00710
C     =1  RURAL EXTINCTION, (23 KM VIS. DEFAULT PROFILE)                 FL00720
C     =2  RURAL EXTINCTION, (5 KM VIS. DEFAULT PROFILE)                  FL00730
C     =3  NAVY MARITIME EXTINCTION,SETS OWN VIS.                         FL00740
C     =4  MARITIME EXTINCTION, 23 KM VIS.    (LOWTRAN 5 MODEL)           FL00750
C     =5  URBAN EXTINCTION, (5 KM VIS DEFAULT PROFILE)                   FL00760
C     =6  TROPOSPHERIC EXTINCTION, (50 KM VIS. DEFAULT PROFILE)          FL00770
C     =7  USER DEFINED  AEROSOL EXTINCTION COEFFICIENTS                  FL00780
C     RECORDS 3.6.1 TO 3.6.3                                             FL00790
C     =8  FOG1 (ADVECTION FOG) EXTINCTION, 0.2 KM VIS.                   FL00800
C     =9  FOG2 (RADIATION FOG) EXTINCTION, 0.5 KM VIS.                   FL00810
C     =10 DESERT EXTINCTION SETS OWN VISIBILITY FROM WIND SPEED          FL00820
C                                                                        FL00830
C                                                                        FL00840
C     ISEASN SELECTS THE SEASONAL DEPENDENCE OF THE PROFILES             FL00850
C     FOR BOTH THE TROPOSPHERIC (2 TO 10 KM) AND                         FL00860
C     STRATOSPHERIC (10 TO 30 KM) AEROSOLS.                              FL00870
C                                                                        FL00880
C     ISEASN=0 DEFAULTS TO SEASON OF MODEL                               FL00890
C     (MODEL 0,1,2,4,6,7) SUMMER                                         FL00900
C     (MODEL 3,5)         WINTER                                         FL00910
C     =1 SPRING-SUMMER                                                   FL00920
C     =2 FALL - WINTER                                                   FL00930
C                                                                        FL00940
C     IVULCN SELECTS BOTH THE PROFILE AND EXTINCTION TYPE                FL00950
C     FOR THE STRATOSPHERIC AEROSOLS AND DETERMINES TRANSITION           FL00960
C     PROFILES ABOVE THE STRATOSPHERE TO 100 KM.                         FL00970
C                                                                        FL00980
C     IVULCN=0 DEFAULT TO STRATOSPHERIC BACKGROUND                       FL00990
C     =1 STRATOSPHERIC BACKGROUND                                        FL01000
C     =2 AGED VOLCANIC EXTINCTION/MODERATE VOLCANIC PROFILE              FL01010
C     =3 FRESH VOLCANIC EXTINCTION/HIGH VOLCANIC PROFILE                 FL01020
C     =4 AGED VOLCANIC EXTINCTION/HIGH VOLCANIC PROFILE                  FL01030
C     =5 FRESH VOLCANIC EXTINCTION/MODERATE VOLCANIC PROFILE             FL01040
C     =6 BACKGROUND STRATOSPHERIC EXTINCTION                             FL01050
C     /MODERATE VOLCANIC PROFILE                                         FL01060
C     =7 BACKGROUND STRATOSPHERIC EXTINCTION                             FL01070
C     /HIGH VOLCANIC PROFILE                                             FL01080
C     =8 FRESH VOLCANIC EXTINCTION/EXTREME VOLCANIC PROFILE              FL01090
C                                                                        FL01100
C     ICSTL IS THE AIR MASS CHARACTER(1 TO 10) ONLY USED WITH            FL01110
C     NAVY MARITIME MODEL (IHAZE=3), DEFAULT VALUE IS 3.                 FL01120
C                                                                        FL01130
C     ICSTL = 1 OPEN OCEAN                                               FL01140
C     .                                                                  FL01150
C     .                                                                  FL01160
C     .                                                                  FL01170
C     10 STRONG CONTINENTAL INFLUENCE                                    FL01180
C                                                                        FL01190
C     ICLD DETERMINES THE INCLUSION OF CIRRUS CLOUD ATTENUATION          FL01200
C     AND GIVES A CHOICE OF FIVE CLOUD MODELS AND 5 RAIN MODELS          FL01210
C                                                                        FL01220
C     ICLD FOR CLOUD AND OR RAIN                                         FL01230
C                                                                        FL01240
C     ICLD = 0  NO CLOUDS OR RAIN                                        FL01250
C     = 1  CUMULUS CLOUD; BASE .66 KM; TOP 3.0 KM                        FL01260
C     = 2  ALTOSTRATUS CLOUD; BASE 2.4 KM; TOP 3.0 KM                    FL01270
C     = 3  STRATUS CLOUD; BASE .33 KM; TOP 1.0 KM                        FL01280
C     = 4  STRATUS/STRATO CU; BASE .66 KM; TOP 2.0 KM                    FL01290
C     = 5  NIMBOSTRATUS CLOUD; BASE .16 KM; TOP .66 KM                   FL01300
C     = 6  2.0 MM/HR DRIZZLE (MODELED WITH CLOUD 3)                      FL01310
C     RAIN 2.0 MM/HR AT 0.0 KM TO 0.22 MM/HR AT 1.5 KM                   FL01320
C     = 7  5.0 MM/HR LIGHT RAIN (MODELED WITH CLOUD 5)                   FL01330
C     RAIN 5.0 MM/HR AT 0.0 KM TO 0.2 MM/HR AT 2.0 KM                    FL01340
C     = 8  12.5 MM/HR MODERATE RAIN (MODELED WITH CLOUD 5)               FL01350
C     RAIN 12.5 MM.HR AT 0.0 KM TO 0.2 MM/HR AT 2.0 KM                   FL01360
C     = 9  25.0 MM/HR HEAVY RAIN (MODELED WITH CLOUD 1)                  FL01370
C     RAIN 25.0 MM/HR AT 0.0 KM TO 0.2 MM/HR AT 3.0 KM                   FL01380
C     =10  75.0 MM/HR EXTREME RAIN (MODELED WITH CLOUD 1)                FL01390
C     RAIN 75.0 MM/HR AT 0.0 KM TO 0.2 MM/HR AT 3.5 KM                   FL01400
C     =11  READ IN USER DEFINED CLOUD EXTINCTION AND ABSORPTION          FL01410
C     =18  STANDARD CIRRUS MODEL                                         FL01420
C     =19  SUB-VISUAL CIRRUS MODEL                                       FL01430
C     =20  NOAA CIRRUS MODEL (LOWTRAN 6 MODEL)                           FL01440
C                                                                        FL01450
C     IVSA DETERMINES THE USE OF THE ARMY VERTICAL STRUCTURE             FL01460
C     ALGORITHM FOR AEROSOLS IN THE BOUNDARY LAYER.                      FL01470
C                                                                        FL01480
C     IVSA=0   NOT USED                                                  FL01490
C     =1   VERTICAL STRUCTURE ALGORITHM                                  FL01500
C                                                                        FL01510
C     VIS =    METEROLOGICAL RANGE (KM) (WHEN SPECIFIED, SUPERSEDES      FL01520
C     DEFAULT VALUE SET BY IHAZE)                                        FL01530
C                                                                        FL01540
C     WSS =    CURRENT WIND SPEED (M/S).                                 FL01550
C     ONLY FOR (IHAZE=3 OR IHAZE=10)                                     FL01560
C     WHH =    24 HOUR AVERAGE WIND SPEED (M/S).  ONLY WITH (IHAZE=3)    FL01570
C                                                                        FL01580
C     RAINRT = RAIN RATE (MM/HR).             DEFAULT VALUE IS ZERO.     FL01590
C     GNDALT = ALTITUDE OF SURFACE RELATIVE TO SEA LEVEL (KM)            FL01600
C     USED TO MODIFY AEROSOL PROFILES BELOW 6 KM ALTITUDE                FL01610
C                                                                        FL01620
C     ****************************************************************** FL01630
C                                                                        FL01640
C     OPTIONAL INPUT RECORDS AFTER RECORD 3.1                            FL01650
C     SELECTED BY PARAMETERS ICLD, IVSA, AND IHAZE ON RECORD 3.1         FL01660
C                                                                        FL01670
C     ****************************************************************** FL01680
C                                                                        FL01690
C     RECORD 3.2     CTHIK,CALT,CEXT,ISEED        (ICLD=18,19,20)        FL01700
C     FORMAT(3F10.3,I10)                                                 FL01710
C     INPUT RECORD FOR CIRRUS ALTITUDE PROFILE                           FL01720
C     SUBROUTINE WHEN ICLD = 18,19,20                                    FL01730
C                                                                        FL01740
C     CHTIK    = CIRRUS THICKNESS (KM)                                   FL01750
C     0  USE THICKNESS STATISTICS                                        FL01760
C     > 0 USER DEFINED THICKNESS                                         FL01770
C                                                                        FL01780
C     CALT     = CIRRUS BASE ALTITUDE (KM)                               FL01790
C     0 USE CALCULATED VALUE                                             FL01800
C     > 0 USER DEFINED BASE ALTITUDE                                     FL01810
C                                                                        FL01820
C     CEXT     = EXTINCTION COEFFIENT (KM-1) AT 0.55 MICRONS             FL01830
C     0 USE 0.14 * CTHIK                                                 FL01840
C     > 0 USER DEFINED EXTINCTION COEFFICIENT                            FL01850
C                                                                        FL01860
C     ISEED    = RANDOM NUMBER INITIALIZATION FLAG.                      FL01870
C     0 USE DEFAULT MEAN VALUES FOR CIRRUS                               FL01880
C     > 0 INITIAL VALUE OF SEED FOR RANDM FUNCTION                        FL01890
C                                                                        FL01900
C                                                                        FL01910
C     ****************************************************************** FL01920
C                                                                        FL01930
C     RECORD 3.3               ZCVSA,ZTVSA,ZINVSA     (IVSA=1)           FL01940
C     FORMAT(3F10.3)                                                     FL01950
C     INPUT RECORD FOR ARMY VERTICAL STRUCTURE                           FL01960
C     ALGORITHM SUBROUTINE WHEN IVSA=1.                                  FL01970
C                                                                        FL01980
C     ZCVSA = CLOUD CEILING HEIGHT (KM) =0 UNKNOWN HEIGHT                FL01990
C                                                                        FL02000
C     ZCVSA > 0  KNOWN CLOUD CEILING                                     FL02010
C     ZCVSA = 0  UNKNOWN CLOUD CEILING HEIGHT                            FL02020
C     PROGRAM CALCULATES CLOUD HEIGHT                                    FL02030
C     ZCVSA < 0  NO CLOUD CEILING                                        FL02040
C                                                                        FL02050
C     ZTVSA = THICKNESS OF CLOUD OR FOG (KM),                            FL02060
C     THICKNESS = 0 DEFAULTS TO 0.2 KM                                   FL02070
C                                                                        FL02080
C     ZINVSA= HEIGHT OF THE INVERSION (KM)                               FL02090
C     = 0   DEFAULTS TO 2 KM (0.2 KM FOR FOG)                            FL02100
C     < 0   NO INVERSION LAYER                                           FL02110
C                                                                        FL02120
C     ****************************************************************** FL02130
C                                                                        FL02140
C     RECORD 3.4     ML,IRD1,IRD2,TITLE       (IAERSL=7)  READ IN LOWTRA FL02150
C     FORMAT(3I5,18A4)                                                   FL02160
C     ADDITIONAL AEROSOL PROFILE                                         FL02170
C                                                                        FL02180
C     ML     = NUMBER OF AEROSOL PROFILES LEVELS TO BE INSERTED          FL02190
C     (MAXIMUM OF 34)                                                    FL02200
C                                                                        FL02210
C     TITLE  = IDENTIFICATION OF NEW MODEL AEROSOL PROFILE               FL02220
C                                                                        FL02230
C                                                                        FL02240
C     RECORD 3.5 IS REPEATED ML TIMES                                    FL02250
C                                                                        FL02260
C     RECORD 3.5                                READ IN AERNSM           FL02270
C     ZMDL,AHAZE,EQLWCZ,RRATZ,IHAZ1,ICLD1,IVUL1,ISEA1,ICHR1               FL02280
C     (IAERSL=7)                                                         FL02290
C     FORMAT(4F10.3,5I5)                                                 FL02300
C                                                                        FL02310
C     ZMDL   = ALTITUDE OF LAYER BOUNDARY (KM)                           FL02320
C                                                                        FL02330
C     AHAZE  = AEROSOL VISIBLE EXTINCTION COFF (KM-1)                    FL02340
C                                                                        FL02350
C     EQLWCZ = LIQUID WATER CONTENT (GM M-3) AT ALT ZMDL                 FL02360
C                                                                        FL02370
C     **** EITHER AHAZE OR EQLWCZ IS ALLOWED ****                        FL02380
C                                                                        FL02390
C     FOR THE AEROSOL, CLOUD OR FOG MODELS                               FL02400
C                                                                        FL02410
C     RRATZ  = RAIN RATE (MM/HR) AT ALT ZMDL                             FL02420
C                                                                        FL02430
C     IHAZ1 AEROSOL MODEL USED FOR SPECTRAL DEPENDENCE OF EXTINCTION      FL02440
C                                                                        FL02450
C     IVUL1 STRATOSPHERIC AERSOL MODEL USED FOR SPECTRAL DEPENDENCE      FL02460
C     OF EXT AT ZMDL                                                     FL02470
C                                                                        FL02480
C     ICLD1 CLOUD MODEL USED FOR SPECTRAL DEPENDENCE OF EXT AT ZMDL      FL02490
C                                                                        FL02500
C     ONLY ONE OF IHAZ1, ICLD1 OR IVUL1 IS ALLOWED                        FL02510
C     IHAZ1 NE 0 OTHERS IGNORED                                           FL02520
C     IHAZ1 EQ 0 AND ICLD1 NE 0 USE ICLD1                                 FL02530
C                                                                        FL02540
C     IF AHAZE AND EQLWCZ ARE BOTH ZERO, DEFAULT PROFILE LOADED          FL02550
C     ACCORDING TO IHAZ1,ICLD1,IVUL1                                     FL02560
C                                                                        FL02570
C     ISEA1 =  AEROSOL SEASON CONTROL FOR THE ALTITUDE ZMDL              FL02580
C                                                                        FL02590
C     ICHR1 =  INDICATES A BOUNDARY CHANGE BETWEEN TWO OR MORE ADJACENT  FL02600
C     USER DEFINED AEROSOL OR CLOUD REGIONS AT ALTITUDE ZMDL             FL02610
C     (REQUIRED FOR IHAZE=7 OR ICLD=11)                                  FL02620
C     NOTE: DEFAULTS TO 0 FOR IHAZE.NE.7 OR ICLD.NE.11                   FL02630
C                                                                        FL02640
C     = 0   NO BOUNDARY CHANGE                                           FL02650
C                                                                        FL02660
C     = 1   SIGNIFIES BOUNDARY CHANGE                                    FL02670
C                                                                        FL02680
C     ****************************************************************** FL02690
C                                                                        FL02700
C     RECORDS 3.6.1 - 3.6.3 READS IN THE USER DEFINED CLOUD EXTINCTION   FL02710
C     AND ABSORPTION        (IHAZE=7 OR ICLD=11)                         FL02720
C                                                                        FL02730
C     RECORD 3.6.1   (IREG(I),I=1,4)                                     FL02740
C     FORMAT (4I5)                                                       FL02750
C                                                                        FL02760
C     IREG   = SPECIFIES WHICH OF THE FOUR ALTITUDE REGIONS A USER       FL02770
C     DEFINED AEROSOL OR CLOUD MODEL WILL USE                            FL02780
C                                                                        FL02790
C     RECORD 3.6.2   AWCCON(N),TITLE(N)                                  FL02800
C     FORMAT (E10.3,18A4)                                                FL02810
C                                                                        FL02820
C     AWCCON(N) = CONVERSION FACTOR FROM EQUIVALENT LIQUID WATER         FL02830
C     CONTENT (GM/M3) TO EXTINCTION COEFFICIENT (KM-1).                  FL02840
C                                                                        FL02850
C     TITLE(N)  = FOR AN AEROSOL OR CLOUD REGION                         FL02860
C                                                                        FL02870
C     RECORD 3.6.3 (VX(I),EXTC(N,I),ABSC(N,I),ASYM(N,I),I=1,47)          FL02880
C     FORMAT (3(F6.2,2F7.5,F6.4)) ** 16 RECORDS **                       FL02890
C                                                                        FL02900
C     VX(I)    = WAVELENGTH OF AEROSOL COEFFICIENT                       FL02910
C     (NOT USED BY PROGRAM BUT CORRESPONDING TO                          FL02920
C     WAVELENGTHS DEFINED IN ARRAY VX2                                   FL02930
C     IN SUBROUTINE EXTDA)                                               FL02940
C                                                                        FL02950
C     EXTC(N,I) = AEROSOL EXTINCTION COEFFICIENT                         FL02960
C     ABSC(N,I) = AEROSOL ABSORPTION COEFFICIENT                         FL02970
C     ASYM(N,I) = AEROSOL ASYMMETRY FACTOR                               FL02980
C                                                                        FL02990
C     *** REPEAT RECORDS 3.6.2 - 3.6.3 N TIMES, WHERE                    FL03000
C     *** N = IREG(1)+IREG(2)+IREG(3)+IREG(4) FROM RECORD 3.6.1          FL03010
C                                                                        FL03020
C     ****************************************************************** FL03030
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUM NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                    FL03040
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),VH(16),TX(16),W(16)       FL03050
      COMMON WPATH(IM2,16),TBBY(IM2)                                     FL03060
      COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),VX2(47),AWCCON(5)
C
      CHARACTER*8      HMOD                                              FL03080
C
      COMMON /CMN/ HMOD(3),ZM(MXZMD),PF(MXZMD),TF(MXZMD),RFNDXM(MXZMD),
     *          ZP(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),
     *          TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),
     *          AMTP(MXMOL,MXPDIM)                                        
C
      COMMON /PATHD/ PBAR(MXLAY),TBAR(MXLAY),AMOUNT(MXMOL,MXLAY),        FA00530
     *               WN2L(MXLAY),DVL(MXLAY),WTOTL(MXLAY),ALBL(MXLAY),    FA00540
     *               ADBL(MXLAY),AVBL(MXLAY),H2OSL(MXLAY),IPATH(MXLAY),  FA00550
     *               ITYL(MXLAY),SECNTA(MXLAY),HT1,HT2,ALTZ(0:MXLAY),    FA00560
     *               PZ(0:MXLAY),TZ(0:MXLAY)                             FA00570
c 
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL03130
     *NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                                 FL03140
     *                     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4       FL03150
      COMMON /LCRD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB   FL03160
      COMMON /ADRIVE/LOWFLG,IREAD,MODELF,ITYPEF,NOZERO,NOPRNF,           FL03170
     * H1F,H2F,ANGLEF,RANGEF,BETAF,LENF,VL1,VL2,RO,IPUNCH,VBAR,          FL03180
     * HMINF,PHIF,IERRF,HSPACE                                           FL03190
      COMMON /LCRD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,    FL03200
     *    RAINRT                                                         FL03210
      COMMON /LCRD2A/ CTHIK,CALT,CEXT                                    FL03220
      COMMON /LCRD2D/ IREG(4),ALTB(4),IREGC(4)                           FL03230
      COMMON /LCRD3/ H1,H2,ANGLE,RANGE,BETA,RE,LEN                       FL03240
      COMMON /LCRD4/ V1,V2,DV                                            FL03250
      REAL*8           V1P,V2P                                           FL03260
      CHARACTER*8       XID,       HMOLID,      YID   
      Real*8                SECANT,       XALTZ
      COMMON /CVRLOW/ HNAMLOW,HVRLOW
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),      FL03280
     *     WK(60),PZL,PZU,TZL,TZU,WN2   ,DVP,V1P,V2P,TBOUNF,EMISIV,      FL03290
     *     FSCDID(17),NMOL,LAYER,YI1,YID(10) ,LSTWDF                     FL03300
      COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                      FL03310
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,N_LVL,JH1              FL03320
      COMMON/MODEL/ ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),                     FL03330
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL03340
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)
      COMMON /MART/ RHH                                                  FL03350
      COMMON /MDLZ/ HMDLZ(10)                                            FL03370
      COMMON /ZVSALY/ ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)             FL03380
      CHARACTER*20 HHAZE,HSEASN,HVULCN,HMET,HMODEL,BLANK                 FL03390
      CHARACTER*24 HTRRAD                                                FL03400
      CHARACTER*18 HNAMLOW,HVRLOW
      COMMON /TITL/ HHAZE(16),HSEASN(2),HVULCN(8),BLANK,                 FL03410
     * HMET(2),HMODEL(8),HTRRAD(4)                                       FL03420
      COMMON /VSBD/ VSB(10)                                              FL03430
C                           
C     ISEED IS INTEGER*4    
C                           
      INTEGER*4 ISEED       
C                                                                        FL03440
C     **   IRD, IPR, AND IPU ARE UNIT NUMBERS FOR INPUT, OUTPUT, AND     FL03500
C     **   TAPE7 RESPECTIVELY                                            FL03510
C                                                                        FL03520
      EQUIVALENCE (FSCDID(5),IEMS),(FSCDID(4),IAERSL)                    FL03530
C                                                                        FL03540
      DATA MAXATM,MAXGEO   /3020, 3014/                                  FL03550
C
      DATA I_1/1/, I_10/10/
C
c%%%%%%%%
c
c     iemsct has been set to 1 here to force the reults from lowtrn
c     to always be available on the lblrtm vertical grid.
c
c     IEMSCT = IEMS                                                      FL03560
      IEMSCT = 1
c
C
C     ASSIGN CVS VERSION NUMBER TO MODULE 
C
      HVRLOW = '$Revision$' 
C                                                                        FL03570
C     ALTITUDE PARAMETERS                                                FL03580
C                                                                        FL03590
C     ZMDL  COMMON/MODEL/  THE ALTITUDES USED IN LOWTRAN                 FL03600
C     ZCVSA,ZTVSA,ZIVSA RECORD 3.3 LOWTRAN FOR VSA INPUT                 FL03610
C     ZM  BLANK COMMON  RETURNS ALTITUDES FOR LBLRTM USE                 FL03620
C     ZP  BLANK COMMON NOT USED BY LOWTRAN                               FL03630
C     ZVSA  NINE ALTITUDES GEN BY VSA ROUTINE                            FL03640
C                                                                        FL03650
      PI = 2.0*ASIN(1.0)                                                 FL03660
      CA = PI/180.                                                       FL03670
      DEG = 1.0/CA                                                       FL03680
C                                                                        FL03690
C     **   GCAIR IS THE GAS CONSTANT FOR AIR IN UNITS OF MB/(GM CM-3 K)  FL03700
C                                                                        FL03710
      GCAIR = 2.87053E+3                                                 FL03720
C                                                                        FL03730
C     **   BIGNUM AND BIGEXP ARE THE LARGEST NUMBER AND THE LARGEST ARGU FL03740
C     **   EXP ALLOWED AND ARE MACHINE DEPENDENT. THE NUMBERS USED HERE  FL03750
C     **   FOR A TYPICAL 32 BIT-WORD COMPUTER.                           FL03760
C                                                                        FL03770
      BIGNUM = 1.0E38                                                    FL03780
      BIGEXP = 87.0                                                      FL03790
      KMAX = 16                                                          FL03800
C                                                                        FL03810
C     **   NL IS THE NUMBER OF BOUNDARIES IN THE STANDARD MODELS 1 TO 6  FL03820
C     **   BOUNDARY    (AT 99999 KM) IS NO LONGER USED                   FL03830
C                                                                        FL03840
      NL = 50                                                            FL03850
      JH1 = 0                                                            FL03860
      IKLO = 1                                                           FL03870
C                                                                        FL03880
C     CC                                                                 FL03890
C     CC    FIX DV TO 5.0 FOR LBLRTM USAGE                               FL03900
C     CC                                                                 FL03910
C                                                                        FL03920
      DV = 5.0                                                           FL03930
C                                                                        FL03940
C     CC                                                                 FL03950
C     CC    OBTAIN PARAMETERS IN COMMON/LCRD3/AND/LCRD4/ FROM COMMON ADR FL03960
C     CC    WHICH PASSED THEM FROM LBLATM                                FL03970
C     CC                                                                 FL03980
C                                                                        FL03990
      DO 10 II = 1, 4                                                    FL04000
         IREG(II) = 0                                                    FL04010
   10 CONTINUE                                                           FL04020
      DO 20 I = 1, 5                                                     FL04030
         DO 18 J = 1, 40                                                 FL04040
            ABSC(I,J) = 0.                                               FL04050
            EXTC(I,J) = 0.                                               FL04060
            ASYM(I,J) = 0.                                               FL04070
 18      CONTINUE 
 20   CONTINUE                                                           FL04080
c
C     CC                                                                 FL04440
C     CC    OBTAIN ITYPE FROM LBLRTM CONTROL AS STORED IN COMMON ADRIVE  FL04450
C     CC                                                                 FL04460
C                                                                        FL04470
      ITYPE = ITYPEF                                                     FL04480
C                                                                        FL04490
      H1 = H1F                                                           FL04090
      NOPRNT = NOPRNF                                                    FL04100
      MODEL = MODELF                                                     FL04110
      MDELS = MODEL                                                      FL04120
      M1 = MODEL                                                         FL04130
      M2 = MODEL                                                         FL04140
      M3 = MODEL                                                         FL04150
c
C                                                                        FL04170
      IF (ITYPE.EQ.1) THEN                                               FL04180
         LENF = 0                                                        FL04190
      ENDIF                                                              FL04200
C                                                                        FL04210
      IF (MODEL.EQ.0) THEN                                               FL04220
         M1 = 0                                                          FL04230
         M2 = 0                                                          FL04240
         M3 = 0                                                          FL04250
      ENDIF                                                              FL04260
      M = MODEL                                                          FL04270
      IM = 0                                                             FL04280
C                                                                        FL04290
      IF (IAERSL.EQ.7) THEN                                              FL04300
         M = 7                                                           FL04310
         IM = 1                                                          FL04320
      ENDIF                                                              FL04330
C                                                                        FL04340
      H2 = H2F                                                           FL04350
      ANGLE = ANGLEF                                                     FL04360
      RANGE = RANGEF                                                     FL04370
      BETA = BETAF                                                       FL04380
      LEN = LENF                                                         FL04390
      V1 = VL1                                                           FL04400
      V2 = VL2                                                           FL04410
      RE = RO                                                            FL04420
C                                                                        FL04430
C     CC                                                                 FL04500
C     CC    SET TBOUND AND SALB TO ZERO NOT UTILIZED HERE                FL04510
C     CC                                                                 FL04520
C                                                                        FL04530
      TBOUND = 0.0                                                       FL04540
      SALB = 0.0                                                         FL04550
C                                                                        FL04560
C     **   START CALCULATION                                             FL04570
C                                                                        FL04580
C                                                                        FL04590
      WRITE (IPR,900)                                                    FL04600
C                                                                        FL04610
C     OBTAIN MODEL PARAMETERS FROM LBLRTM  (     RECORD 2.1)             FL04620
C                                                                        FL04630
      JPRT = 0                                                           FL04640
      WRITE (IPR,905) MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT              FL04650
      NPR = NOPRNT                                                       FL04660
C                                                                        FL04670
C     **   RECORD 3.1 AEROSOL MODEL                                      FL04680
C                                                                        FL04690
      READ (IRD,910) IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,    FL04700
     *   RAINRT,GNDALT                                                   FL04710
      IF (IHAZE.EQ.3) THEN                                               FL04720
         IF (V1.LT.250.0.OR.V2.LT.250.0) IHAZE = 4                       FL04730
         IF (IHAZE.EQ.4) WRITE (IPR,930)                                 FL04740
      ENDIF                                                              FL04750
      WRITE (IPR,920) IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,   FL04760
     *   RAINRT,GNDALT                                                   FL04770
      IF (GNDALT.GT.0.) WRITE (IPR,915) GNDALT                           FL04780
      IF (GNDALT.GE.6.0) THEN                                            FL04790
         WRITE (IPR,925) GNDALT                                          FL04800
         GNDALT = 0.                                                     FL04810
      ENDIF                                                              FL04820
C                                                                        FL04830
      IF (VIS.LE.0.0.AND.IHAZE.GT.0) VIS = VSB(IHAZE)                    FL04840
      RHH = 0.                                                           FL04850
      CALL YDIAR (IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT FL04860
     *   ,GNDALT,YID)                                                    FL04870
      IF (MODEL.EQ.0) GO TO 30                                           FL04880
      IF ((MODEL.EQ.3.OR.MODEL.EQ.5).AND.ISEASN.EQ.0) ISEASN = 2         FL04890
C                                                                        FL04900
C     **WARNING** IF V1 OR V2 LESS THEN 250 CM-1 PROGRAM WILL NOT        FL04910
C     PERMIT USE OF NAVY MARITIME (IHAZE=3) SWITCHES TO IHAZE=4          FL04920
C                                                                        FL04930
      ICH(1) = IHAZE                                                     FL04940
      ICH(2) = 6                                                         FL04950
      ICH(3) = 9+IVULCN                                                  FL04960
   30 IF (RAINRT.EQ.0) GO TO 40                                          FL04970
      WRITE (IPR,935) RAINRT                                             FL04980
   40 ICH(4) = 18                                                        FL04990
      ICH(1) = MAX(ICH(1),I_1)                                             FL05000
      ICH(3) = MAX(ICH(3),I_10)                                            FL05010
      IF (ICLD.GE.1.AND.ICLD.LE.11) THEN                                 FL05020
         ICH(4) = ICH(3)                                                 FL05030
         ICH(3) = ICH(2)                                                 FL05040
         ICH(2) = ICLD                                                   FL05050
      ENDIF                                                              FL05060
C                                                                        FL05070
C     CC   IF(ICH(4).LE.9) ICH(4)=10                                     FL05080
C                                                                        FL05090
      IFLGA = 0                                                          FL05100
      IFLGT = 0                                                          FL05110
      CTHIK = -99.                                                       FL05120
      CALT = -99.                                                        FL05130
      ISEED = -99                                                        FL05140
      IF (ICLD.LT.18) GO TO 50                                           FL05150
C                                                                        FL05160
C     **   RECORD 3.2 CIRRUS CLOUDS                                      FL05170
C                                                                        FL05180
      READ (IRD,940) CTHIK,CALT,CEXT,ISEED                               FL05190
      WRITE (IPR,945) CTHIK,CALT,CEXT,ISEED                              FL05200
   50 CONTINUE                                                           FL05210
C                                                                        FL05220
C     **   RECORD 3.3 VERTICAL STRUCTURE ALGORITHM                       FL05230
C                                                                        FL05240
      ZCVSA = -99.                                                       FL05250
      ZTVSA = -99.                                                       FL05260
      ZINVSA = -99.                                                      FL05270
C                                                                        FL05280
      IF (IVSA.ne.0) then                                                   FL05290
         READ (IRD,950) ZCVSA,ZTVSA,ZINVSA                                  FL05300
         WRITE (IPR,955) ZCVSA,ZTVSA,ZINVSA                                 FL05310
C                                                                        FL05320
         CALL VSA 
     *        (IHAZE,VIS,ZCVSA,ZTVSA,ZINVSA,ZVSA,RHVSA,AHVSA,IHVSA)         FL05330
C                                                                        FL05340
C     END OF VSA MODEL SET-UP                                            FL05350
c
      endif
C                                                                        FL05360
      IF (MODEL.NE.0) ML = NL                                            FL05370
C                                                                        FL05380
      IF (MDELS.NE.0) HMODEL(7) = HMODEL(MDELS)                          FL05390
      IF (MDELS.EQ.0) HMODEL(7) = HMODEL(8)                              FL05400
C                                                                        FL05410
C                                                                        FL05420
      IF (IAERSL.EQ.7) THEN                                              FL05430
C                                                                        FL05440
C        **   RECORD 3.4 USER SUPPLIED AEROSOL AND CLOUD PROFILE         FL05450
C                                                                        FL05460
         READ (IRD,960) ML,HMODEL(7)                                     FL05470
         WRITE (IPR,965) ML,HMODEL(7)                                    FL05480
      ENDIF                                                              FL05490
      M = 7                                                              FL05500
      CALL AERNSM (IAERSL,JPRT,GNDALT)                                   FL05510
      IF (ICLD.LT.20) GO TO 70                                           FL05520
C                                                                        FL05530
C     SET UP CIRRUS MODEL                                                FL05540
C                                                                        FL05550
      IF (CTHIK.NE.0) IFLGT = 1                                          FL05560
      IF (CALT.NE.0) IFLGA = 1                                           FL05570
      IF (ISEED.EQ.0) IFLGT = 2                                          FL05580
      IF (ISEED.EQ.0) IFLGA = 2                                          FL05590
      CALL CIRRUS (CTHIK,CALT,ISEED,CPROB,MDELS)                         FL05600
      WRITE (IPR,970)                                                    FL05610
      IF (IFLGT.EQ.0) WRITE (IPR,975) CTHIK                              FL05620
      IF (IFLGT.EQ.1) WRITE (IPR,980) CTHIK                              FL05630
      IF (IFLGT.EQ.2) WRITE (IPR,985) CTHIK                              FL05640
      IF (IFLGA.EQ.0) WRITE (IPR,990) CALT                               FL05650
      IF (IFLGA.EQ.1) WRITE (IPR,995) CALT                               FL05660
      IF (IFLGA.EQ.2) WRITE (IPR,1000) CALT                              FL05670
      WRITE (IPR,1005) CPROB                                             FL05680
C                                                                        FL05690
C     END OF CIRRUS MODEL SET UP                                         FL05700
C                                                                        FL05710
   70 CONTINUE                                                           FL05720
C                                                                        FL05730
C     **   RECORD 3.6                                                    FL05740
C                                                                        FL05750
      IF ((IHAZE.EQ.7).OR.(ICLD.EQ.11)) THEN                             FL05760
C                                                                        FL05770
C        **   RECORDS 3.6.1 - 3.6.3                                      FL05780
C        **           USER SUPPLIED AEROSOL EXTINCTION AND ABSORPTION    FL05790
C                                                                        FL05800
         CALL RDEXA                                                      FL05810
      ENDIF                                                              FL05820
C                                                                        FL05830
C     WRITE(IPR,1313)H1,H2,ANGLE,RANGE,BETA,RO,LEN                       FL05840
C     1313 FORMAT('0 RECORD 2.2 ****',6F10.3,I5)                         FL05850
C                                                                        FL05860
      GO TO 80                                                           FL05870
C                                                                        FL05880
C     **   RO IS THE RADIUS OF THE EARTH                                 FL05890
C                                                                        FL05900
   80 RE = 6371.23                                                       FL05910
      IF (MODEL.EQ.1) RE = 6378.39                                       FL05920
      IF (MODEL.EQ.4) RE = 6356.91                                       FL05930
      IF (MODEL.EQ.5) RE = 6356.91                                       FL05940
      IF (RO.GT.0.0) RE = RO                                             FL05950
C                                                                        FL05960
C                                                                        FL05970
C     IPH   =-99                                                         FL05980
C     IDAY  =-99                                                         FL05990
C     ISOURC=-99                                                         FL06000
C                                                                        FL06010
C     ANGLEM=-99.                                                        FL06020
C                                                                        FL06030
      WRITE (IPR,1010) HTRRAD(IEMSCT+1)                                  FL06040
      MDEL = MODEL                                                       FL06050
      IF (MDEL.EQ.0) MDEL = 8                                            FL06060
      MM1 = MDEL                                                         FL06070
      MM2 = MDEL                                                         FL06080
      MM3 = MDEL                                                         FL06090
      IF (M1.NE.0) MM1 = M1                                              FL06100
      IF (M2.NE.0) MM2 = M2                                              FL06110
      IF (M3.NE.0) MM3 = M3                                              FL06120
      WRITE (IPR,1015) MM1,HMODEL(MM1),MM2,HMODEL(MM2),MM3,HMODEL(MM3)   FL06130
C                                                                        FL06140
      IF (JPRT.EQ.0) GO TO 90                                            FL06150
      IF (ISEASN.EQ.0) ISEASN = 1                                        FL06160
      IVULCN = MAX(IVULCN,I_1)                                             FL06170
      IHVUL = IVULCN+10                                                  FL06180
      IF (IVULCN.EQ.6) IHVUL = 11                                        FL06190
      IF (IVULCN.EQ.7) IHVUL = 11                                        FL06200
      IF (IVULCN.EQ.8) IHVUL = 13                                        FL06210
      IHMET = 1                                                          FL06220
      IF (IVULCN.GT.1) IHMET = 2                                         FL06230
      IF (IHAZE.EQ.0) GO TO 90                                           FL06240
      WRITE (IPR,1020) HHAZE(IHAZE),VIS,HHAZE(6),HHAZE(6),HSEASN(ISEASN) FL06250
     *   ,HHAZE(IHVUL),HVULCN(IVULCN),HSEASN(ISEASN),HHAZE(15),          FL06260
     *   HMET(IHMET)                                                     FL06270
   90 CONTINUE                                                           FL06280
      IF (ITYPE.EQ.1) WRITE (IPR,1025) H1,RANGE                          FL06290
      IF (ITYPE.EQ.2) WRITE (IPR,1030) H1,H2,ANGLE,RANGE,BETA,LEN        FL06300
      IF (ITYPE.EQ.3) WRITE (IPR,1035) H1,H2,ANGLE                       FL06310
C                                                                        FL06320
C                                                                        FL06330
C                                                                        FL06340
C                                                                        FL06350
      ALAM1 = 1.0E38                                                     FL06360
      IF (V1.GT.0.) ALAM1 = 10000./V1                                    FL06370
      ALAM2 = 10000./V2                                                  FL06380
      DV = MAX(DV,5.)                                                    FL06390
      DV =  REAL(INT(DV/5+0.1))*5.0                                      FL06400
      IF (ALAM1.GT.999999.) ALAM1 = 999999.                              FL06410
      WRITE (IPR,1040) V1,ALAM1,V2,ALAM2,DV                              FL06420
C                                                                        FL06430
C     **   LOAD ATMOSPHERIC PROFILE INTO /MODEL/                         FL06440
C                                                                        FL06450
      CALL STDMDL                                                        FL06460
C                                                                        FL06470
      IF (IEMSCT.EQ.1) CALL NEWMDL (MAXATM)                              FL06480
C                                                                        FL06490
C     **   TRACE PATH THROUGH THE ATMOSPHERE AND CALCULATE ABSORBER AMOU FL06500
C                                                                        FL06510
      ISSGEO = 0                                                         FL06520
      MODEL = MDELS                                                      FL06530
      CALL GEO (IERROR,BENDNG,MAXGEO)                                    FL06540
C                                                                        FL06550
C                                                                        FL06560
C     FINAL SET OF LAYERS                                                FL06570
C                                                                        FL06580
C                                                                        FL06590
C                                                                        FL06600
      IF (IERROR.GT.0) GO TO 100                                         FL06610
C                                                                        FL06620
C                                                                        FL06630
C                                                                        FL06640
C     **   LOAD AEROSOL EXTINCTION AND ABSORPTION COEFFICIENTS           FL06650
C                                                                        FL06660
C     CC                                                                 FL06670
C     CC    LOAD EXTINCTIONS AND ABSORPTIONS FOR 0.2-200.0 UM (1-46)     FL06680
C     CC                                                                 FL06690
C     CC   CALL EXABIN                                                   FL06700
C     CC                                                                 FL06710
C     CC    CALCULATE EQUIVALENT LIQUID WATER CONSTANTS                  FL06720
C     CC                                                                 FL06730
C                                                                        FL06740
      CALL EQULWC                                                        FL06750
C                                                                        FL06760
C                                                                        FL06770
C                                                                        FL06780
      CALL TRANS                                                         FL06790
C                                                                        FL06800
  100 CONTINUE                                                           FL06810
C                                                                        FL06820
      LOWFLG = 0                                                         FL06830
      RETURN                                                             FL06840
C                                                                        FL06850
  900 FORMAT('1',20X,'***** LOWTRAN 7 (MODIFIED)*****')                  FL06860
  905 FORMAT('0 RECORD 2.1 ****',8I5       )                             FL06870
  910 FORMAT(6I5,5F10.3)                                                 FL06880
  915 FORMAT('0','  GNDALT =',F10.2)                                     FL06890
  920 FORMAT('0 RECORD 3.1 ****',6I5,5F10.3)                             FL06900
  925 FORMAT('0 GNDALT GT 6.0 RESET TO ZERO, GNDALT WAS',F10.3)          FL06910
  930 FORMAT('0**WARNING** NAVY MODEL IS NOT USEABLE BELOW 250CM-1'/     FL06920
     * 10X,' PROGRAM WILL SWITCH TO IHAZE=4 LOWTRAN 5 MARITIME'//)       FL06930
  935 FORMAT('0 RAIN MODEL CALLED, RAIN RATE = ',F9.2,' MM/HR')          FL06940
  940 FORMAT(3F10.3,I10)                                                 FL06950
  945 FORMAT('0 RECORD 2A *****',3F10.3,I10)                             FL06960
  950 FORMAT(3F10.3)                                                     FL06970
  955 FORMAT('0 RECORD 3.3 ****',3F10.3)                                 FL06980
  960 FORMAT(I5,18A4)                                                    FL06990
  965 FORMAT('0 RECORD 3.4 ****',I5,18A4)                                FL07000
  970 FORMAT(15X,'CIRRUS ATTENUATION INCLUDED')                          FL07010
  975 FORMAT(15X,'CIRRUS ATTENUTION STATISTICALLY DETERMENED TO BE',     FL07020
     * F10.3,'KM')                                                       FL07030
  980 FORMAT(15X,'CIRRUS THICKNESS USER DETERMINED TO BE',F10.3,'KM')    FL07040
  985 FORMAT(15X,'CIRRUS THICKNESS DEFAULTED TO MEAN VALUE OF    ',      FL07050
     * F10.3,'KM')                                                       FL07060
  990 FORMAT(15X,'CIRRUS BASE ALTITUDE STATISCALLY DETERMINED TO BE',    FL07070
     * F10.3,' KM')                                                      FL07080
  995 FORMAT(15X,'CIRRUS BASE ALTITUDE USER DETERMINED TO BE',           FL07090
     * F10.3,' KM')                                                      FL07100
 1000 FORMAT(15X,'CIRRUS BASE ALTITUDE DEFAULTED TO MEAN VALUE OF',      FL07110
     * F10.3,'KM')                                                       FL07120
 1005 FORMAT(15X,'PROBABILTY OF CLOUD OCCURRING IS',F7.1,                FL07130
     * ' PERCENT')                                                       FL07140
 1010 FORMAT('0 PROGRAM WILL COMPUTE ',A24)                              FL07150
 1015 FORMAT('0 ATMOSPHERIC MODEL',/,                                    FL07160
     * 10X,'TEMPERATURE = ',I4,5X,A20,/,                                 FL07170
     * 10X,'WATER VAPOR = ',I4,5X,A20,/,                                 FL07180
     * 10X,'OZONE       = ',I4,5X,A20)                                   FL07190
 1020 FORMAT('0 AEROSOL MODEL',/,10X,'REGIME',                           FL07200
     * T35,'AEROSOL TYPE',T60,'PROFILE',T85,'SEASON',/,/,                FL07210
     * 10X,'BOUNDARY LAYER (0-2 KM)',T35,A20,T60,F5.1,                   FL07220
     * ' KM VIS AT SEA LEVEL',/,10X,'TROPOSPHERE  (2-10KM)',T35,         FL07230
     * A20,T60,A20,T85,A20,/,10X,'STRATOSPHERE (10-30KM)',               FL07240
     * T35,A20,T60,A20,T85,A20,/,10X,'UPPER ATMOS (30-100KM)',           FL07250
     * T35,A20,T60,A20)                                                  FL07260
 1025 FORMAT('0 HORIZONTAL PATH',/,10X,'ALTITUDE = ',F10.3,' KM',/,      FL07270
     *    10X,'RANGE    = ',F10.3,' KM')                                 FL07280
 1030 FORMAT('0 SLANT PATH, H1 TO H2',/,                                 FL07290
     *    10X,'H1    = ',F10.3,' KM',/,10X,'H2    = ',F10.3,' KM',/,     FL07300
     *    10X,'ANGLE = ',F10.3,' DEG',/,10X,'RANGE = ',F10.3,' KM',/,    FL07310
     *    10X,'BETA  = ',F10.3,' DEG',/,10X,'LEN   = ',I6)               FL07320
 1035 FORMAT('0 SLANT PATH TO SPACE',/,                                  FL07330
     *    10X, 'H1    = ',F10.3,' KM',/,10X,'HMIN  = ',F10.3,' KM',/,    FL07340
     *    10X,'ANGLE = ',F10.3,' DEG')                                   FL07350
 1040 FORMAT('0 FREQUENCY RANGE '/,10X,' V1 = ',F12.1,' CM-1  (',        FL07360
     * F10.2,' MICROMETERS)',/,10X,' V2 = ',F12.1,' CM-1  (',F10.2,      FL07370
     * ' MICROMETERS)',/10X,' DV = ',F12.1,' CM-1')                      FL07380
C                                                                        FL07390
      END                                                                FL07400
      SUBROUTINE AERNSM(IAERSL,JPRT,GNDALT)                              FL07410

      PARAMETER (NCASE=15)

      CHARACTER*1 JCHAR                                                  FL07420
C                                                                        FL07430
C     ****************************************************************** FL07440
C     DEFINES ALTITUDE DEPENDENT VARIABLES Z,P,T,WH,WO AND HAZE          FL07450
C     CLD RAIN  CLDTYPE                                                  FL07460
C     LOADS HAZE INTO APPROPRATE LOCATION                                FL07470
C     ****************************************************************** FL07480
C                                                                        FL07490
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C     
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELFAS(MXZMD),HSTOR(MXZMD),ICH(4),VH(16),TX(16),W(16)       FL07500
      COMMON WPATH(IM2,16),TBBY(IM2)                                     FL07510
      COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),VX2(47),AWCCON(5)
C
      CHARACTER*8      HMOD                                              FL07530
C
      COMMON /CMN/ HMOD(3),ZM(MXZMD),PF(MXZMD),TF(MXZMD),RFNDXM(MXZMD),
     *          ZP(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),
     *          TPSUM(IM2),RHOPSM(IM2),IMMAX,WGM(MXZMD),DENW(MXZMD),
     *          AMTP(MXMOL,MXPDIM)                                        
C
      COMMON /PATHD/ PBAR(MXLAY),TBAR(MXLAY),AMOUNT(MXMOL,MXLAY),        FA00530
     *               WN2L(MXLAY),DVL(MXLAY),WTOTL(MXLAY),ALBL(MXLAY),    FA00540
     *               ADBL(MXLAY),AVBL(MXLAY),H2OSL(MXLAY),IPATH(MXLAY),  FA00550
     *               ITYL(MXLAY),SECNTA(MXLAY),HT1,HT2,ALTZ(0:MXLAY),    FA00560
     *               PZ(0:MXLAY),TZ(0:MXLAY)                             FA00570
c 
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL07580
     *NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                                 FL07590
     *                     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4       FL07600
      COMMON /LCRD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB   FL07610
      COMMON /CARD1B/ JUNIT(NCASE),WMOL(NCASE),WAIR1,JLOW                FL07620
      COMMON /LCRD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,    FL07630
     *    RAINRT                                                         FL07640
      COMMON /LCRD2A/ CTHIK,CALT,CEXT                                    FL07650
      COMMON /LCRD2D/ IREG(4),ALTB(4),IREGC(4)                           FL07660
      COMMON /LCRD3/ H1,H2,ANGLE,RANGE,BETA,RE,LEN                       FL07670
      COMMON /LCRD4/ V1,V2,DV                                            FL07680
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,N_LVL,JH1              FL07690
      COMMON /MART/ RHH                                                  FL07700
c     COMMON /MDATA/ ZDA(MXZMD),P(MXZMD),T(MXZMD),WH(MXZMD),WO(MXZMD),   FL07710
c    *     HMIX(MXZMD),CLD(MXZMD,7),RR(MXZMD,7)                          FL07720
      COMMON /MDATA/                              WH(MXZMD),WO(MXZMD),   FL07710
     *                 CLD(MXZMD,7),RR(MXZMD,7)                          FL07720
      COMMON /MDATA2/ZDA(MXZMD),P(MXZMD),T(MXZMD)

      COMMON /MODEL/ ZMDL(MXZMD),PMM(MXZMD),TMM(MXZMD),                  FL07730
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL07740
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)
      COMMON /ZVSALY/ ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)             FL07750
      COMMON /MDLZ/HMDLZ(10)                                             FL07760
      COMMON /TITL/ HZ(16),SEASN(2),VULCN(8),BLANK,                      FL07790
     *     HMET(2),HMODEL(8),HTRRAD(4)                                   FL07800
      DIMENSION ITY1(MXZMD+1),IH1(MXZMD),IS1(MXZMD),IVL1(MXZMD),
     *     ZGN(MXZMD)                                                    FL07810
      DIMENSION INEW(MXZMD),RELHUM(MXZMD),ZSTF(MXZMD),CLDTOP(10),        FL07820
     *     AHAST(MXZMD)
C                                                                        FL07830
      CHARACTER*20 HZ,SEASN,VULCN,HMET,HMODEL,BLANK
      CHARACTER*24 HTRRAD
      CHARACTER*20 AHOL1,AHOL2,AHOL3,AHLVSA,AHUS                         FL07840
      CHARACTER*20 AHAHOL(NCASE),HHOL                                    FL07850
      DIMENSION  JCHAR(NCASE)                                            FL07860
C
      DATA I_1/1/, I_12/12/, I_32/32/, I_34/34/
C
      DATA AHLVSA/'VSA DEFINED         '/                                FL07870
      DATA  AHUS /'USER DEFINED        '/                                FL07880
      DATA AHAHOL/                                                       FL07890
     * 'CUMULUS             ',                                           FL07900
     * 'ALTOSTRATUS         ',                                           FL07910
     * 'STRATUS             ',                                           FL07920
     * 'STRATUS STRATO CUM  ',                                           FL07930
     * 'NIMBOSTRATUS        ',                                           FL07940
     * 'DRIZZLE 2.0 MM/HR   ',                                           FL07950
     * 'LT RAIN 5.0 MM/HR   ',                                           FL07960
     * 'MOD RAIN 12.5 MM/HR ',                                           FL07970
     * 'HEAVY RAIN 25 MM/HR ',                                           FL07980
     * 'EXTREME RAIN 75MM/HR',                                           FL07990
     * 'USER ATMOSPHERE     ',                                           FL08000
     * 'USER RAIN NO CLOUD  ',                                           FL08010
     * 'CIRRUS CLOUD        ',                                           FL08020
     * 'SUB-VISUAL CIRRUS   ',                                           FL08030
     * 'NOAA CIRRUS MODEL   '/                                           FL08040
      DATA CLDTOP / 3.,3.,1.,2.,.66,1.,.66,.66,3.,3./                    FL08050
C                                                                        FL08060
C     F(A) IS SATURATED WATER WAPOR DENSITY AT TEMP T,A=TZERO/T          FL08070
C                                                                        FL08080
      F(A) = EXP(18.9766-14.9595*A-2.43882*A*A)*A                        FL08090
C                                                                        FL08100
C     ZM ORIGINALLY IS LBLRTM ALT                                        FL08110
C                                                                        FL08120
C     ZGN IS EFFICTIVE ALTITUDE ARRAY                                    FL08130
C     ZDA COMMON   /MDATA/  ALTITUDE OF THE PRESSURES,TEMP IN MDATA      FL08140
C     ZMDL COMMON /MODEL/ FINAL ALTITUDE FOR LOWTRAN                     FL08150
C     ZSTF  STORAGE OF ORIGINAL LBLRTM ALTITUDES                         FL08160
C     ZK  EFFECTIVE ALTITUDE FOR CLOUD                                   FL08170
C     ZSC EFFECTIVE ALTITUDE FOR AEROSOLS                                FL08180
C     ZP  BLANK COMMON  UNUSED                                           FL08190
C     ZM,PM,TM  ARE FOR LBLRTM USE BETWEEN 0 AND 6 KM                    FL08200
C                                                                        FL08210
      IREGC(1) = 0                                                       FL08220
      IREGC(2) = 0                                                       FL08230
      IREGC(3) = 0                                                       FL08240
      IREGC(4) = 0                                                       FL08250
      ICL = 0                                                            FL08260
      MLSV = ML                                                          FL08270
      DO 10 I = 0, n_lvl-1
         ZSTF(I) = ALTZ(I)                                                 FL08290
   10 CONTINUE                                                           FL08300
      ICONV = 1                                                          FL08310
      IRD0 = 1                                                           FL08320
      ICLDL = ICLD                                                       FL08330
      IF ((MODEL.GT.0.).AND.(MODEL.LT.7)) IRD0 = 0                       FL08340
      IF ((IRD0.EQ.1).AND.(IVSA.EQ.1)) THEN                              FL08350
         IRD0 = 0                                                        FL08360
         IRD1 = 0                                                        FL08370
C                                                                        FL08380
C        C         IRD2 = 0                                              FL08390
C        C         IF(IAERSL .EQ. 7) IRD2 = 1                            FL08400
C                                                                        FL08410
         ICONV = 0                                                       FL08420
         ML = ML+10-JLOW                                                 FL08430
         IF (ML.GT.34) WRITE (IPR,905)                                   FL08440
         ML = MIN(ML,I_34)                                                 FL08450
         ZVSA(10) = ZVSA(9)+0.01                                         FL08460
         RHVSA(10) = 0.                                                  FL08470
         AHVSA(10) = 0.                                                  FL08480
         IHVSA(10) = 0                                                   FL08490
         IF (MODEL.EQ.0) WRITE (IPR,900)                                 FL08500
         IF (MODEL.EQ.0) STOP                                            FL08510
      ENDIF                                                              FL08520
      ICL = 0                                                            FL08530
      IDSR = 0                                                           FL08540
      IF (ICLD.EQ.18.OR.ICLD.EQ.19) THEN                                 FL08550
         CALL CIRR18                                                     FL08560

c     cloud is between cl1 and cld2
c     transition regions cloudo>cloud1 and cloud2>cloud3

         cld1 = calt
         cld2 = calt + cthik

         cld0 = cld1 - 0.010
         IF (CLD0.LE.0.) CLD0 = 0.                                       FL08590
         cld3 = cld2 + 0.010
C                                                                        FL08630
      ENDIF                                                              FL08670
c
      CALL FLAYZ 
     *   (ML,MODEL,ICLD,IAERSL,ZMDL,altz,n_lvl,GNDALT,IVSA,IEMSCT)      
c
      DO 20 I = 1, ML                                                    FL08690
         JPRT = 1                                                        FL08700
         IF (MODEL.EQ.0.OR.MODEL.EQ.7) JPRT = 0                          FL08710
         IF (IAERSL.EQ.7) JPRT = 0                                       FL08720
         IF (ICLD.GT.0) JPRT = 0                                         FL08730
         IF (IVSA.GT.0) JPRT = 0                                         FL08740
         HAZEC(I) = 0.0                                                  FL08750
   20 CONTINUE                                                           FL08760
      DO 30 II = 1, 4                                                    FL08770
         ALTB(II) = 0.                                                   FL08780
   30 CONTINUE                                                           FL08790
      T0 = 273.15                                                        FL08800
      IC1 = 1                                                            FL08810
      N = 7                                                              FL08820
      IF (ML.EQ.1) M = 0                                                 FL08830
      IVULCN = MAX(IVULCN,I_1)                                             FL08840
      ISEASN = MAX(ISEASN,I_1)                                             FL08850
      IF (JPRT.EQ.0) THEN                                                FL08860
         WRITE (IPR,950) MODEL,ICLD                                      FL08870
      ENDIF                                                              FL08880
      IF (IAERSL.EQ.7) WRITE (IPR,910)                                   FL08890
C                                                                        FL08900
      KLO = 1                                                            FL08910
C                                                                        FL08920
      IF (IAERSL.NE.7) THEN                                              FL08930
         DO 50 I = 1, ML                                                 FL08940
            INEW(I) = KLO-1                                              FL08950
            IF (ZMDL(I).LT.ALTZ(KLO)) GO TO 50                             FL08960
   40       INEW(I) = KLO                                                FL08970
            KLO = KLO+1                                                  FL08980
            IF (KLO.GT.MLSV) GO TO 50                                    FL08990
            IF (ZMDL(I).GT.ALTZ(KLO)) GO TO 40                             FL09000
   50    CONTINUE                                                        FL09010
      ENDIF                                                              FL09020
C                                                                        FL09030
C                                                                        FL09040
      DO 220 K = 1, ML                                                   FL09050
C                                                                        FL09060
C        LOOP OVER LAYERS                                                FL09070
C                                                                        FL09080
         RH = 0.                                                         FL09090
         WH(K) = 0.                                                      FL09100
         WO(K) = 0.                                                      FL09110
         DP = 0                                                          FL09120
         IHAZ1 = 0                                                        FL09130
         ICLD1 = 0                                                       FL09140
         ISEA1 = 0                                                       FL09150
         IVUL1 = 0                                                       FL09160
         VIS1 = 0.                                                       FL09170
         AHAZE = 0.                                                      FL09180
         EQLWCZ = 0.                                                     FL09190
         RRATZ = 0.                                                      FL09200
         ICHR = 0                                                        FL09210
         DO 60 KM = 1, 15                                                FL09220
            JCHAR(KM) = ' '                                              FL09230
            WMOL(KM) = 0.                                                FL09240
   60    CONTINUE                                                        FL09250
         DO 70 KM = 1, 15                                                FL09260
            JUNIT(KM) = JOU(JCHAR(KM))                                   FL09270
   70    CONTINUE                                                        FL09280
         JUNIT(1) = M1                                                   FL09290
         JUNIT(2) = M1                                                   FL09300
         JUNIT(3) = M2                                                   FL09310
         JUNIT(4) = 6                                                    FL09320
         JUNIT(5) = M3                                                   FL09330
         JUNIT(6) = 6                                                    FL09340
         JUNIT(7) = 6                                                    FL09350
         JUNIT(8) = 6                                                    FL09360
         JUNIT(9) = 6                                                    FL09370
         JUNIT(10) = 6                                                   FL09380
         JUNIT(11) = 6                                                   FL09390
         JUNIT(12) = 6                                                   FL09400
         JUNIT(13) = 6                                                   FL09410
         JUNIT(14) = 6                                                   FL09420
         JUNIT(15) = 6                                                   FL09430
C                                                                        FL09440
C                                                                        FL09450
C        AHAZE =  AEROSOL VISIBLE EXTINCTION COFF (KM-1)                 FL09460
C        AT A WAVELENGTH OF 0.55 MICROMETERS                             FL09470
C                                                                        FL09480
C        EQLWCZ=LIQUID WATER CONTENT (GM M-3) AT ALT Z                   FL09490
C        FOR AEROSOL, CLOUD OR FOG MODELS                                FL09500
C                                                                        FL09510
C        RRATZ=RAIN RATE (MM/HR) AT ALT Z                                FL09520
C                                                                        FL09530
C        IHAZ1 AEROSOL MODEL USED FOR SPECTRAL DEPENDENCE OF EXTINCTION   FL09540
C                                                                        FL09550
C        IVUL1 STRATOSPHERIC AERSOL MODEL USED FOR SPECTRAL DEPENDENCE   FL09560
C        OF EXT AT Z                                                     FL09570
C                                                                        FL09580
C        ICLD1 CLOUD MODEL USED FOR SPECTRAL DEPENDENCE OF EXT AT Z      FL09590
C                                                                        FL09600
C        ONLY ONE OF IHAZ1,ICLD1  OR IVUL1 IS ALLOWED                     FL09610
C        IHAZ1 NE 0 OTHERS IGNORED                                        FL09620
C        IHAZ1 EQ 0 AND ICLD1 NE 0 USE ICLD1                              FL09630
C                                                                        FL09640
C        IF AHAZE AND EQLWCZ ARE BOUTH ZERO                              FL09650
C        DEFAULT PROFILE ARE LOADED FROM IHAZ1,ICLD1,IVUL1                FL09660
C        ISEA1 = AERSOL SEASON CONTROL FOR ALTITUDE Z                    FL09670
C                                                                        FL09680
C        C    IF(IRD2 .EQ. 1) THEN                                       FL09690
C                                                                        FL09700
         IF (IAERSL.EQ.7) THEN                                           FL09710
            READ (IRD,915)  ZMDL(K),AHAZE,EQLWCZ,RRATZ,IHAZ1,ICLD1,
     *         IVUL1,ISEA1,ICHR                                    
            WRITE (IPR,915) ZMDL(K),AHAZE,EQLWCZ,RRATZ,IHAZ1,ICLD1,
     *         IVUL1,ISEA1,ICHR                                    
C                                                                        FL09760
            IF (ICHR.EQ.1) THEN                                          FL09770
               IF (IHAZ1.EQ.0) THEN                                       FL09780
                  IF (ICLD1.NE.11) ICHR = 0                              FL09790
               ELSE                                                      FL09800
                  IF (IHAZ1.NE.7) ICHR = 0                                FL09810
               ENDIF                                                     FL09820
            ENDIF                                                        FL09830
            INEW(K) = KLO-1                                              FL09840
            IF (ZMDL(K).LT.ALTZ(KLO)) GO TO 90                             FL09850
   80       INEW(K) = KLO                                                FL09860
            KLO = KLO+1                                                  FL09870
            IF (KLO.GT.MLSV) GO TO 90                                    FL09880
            IF (ZMDL(K).GT.ALTZ(KLO)) GO TO 80                             FL09890
   90       CONTINUE                                                     FL09900
         ENDIF                                                           FL09910
         IF (IAERSL.NE.7) THEN                                           FL09920
            RRATZ = RAINRT                                               FL09930
            IF (ZMDL(K).GT.6.) RRATZ = 0.                                FL09940
         ENDIF                                                           FL09950
C                                                                        FL09960
C                                                                        FL09970
C        GNDALT NOT ZERO                                                 FL09980
C                                                                        FL09990
         ZSC = ZMDL(K)                                                   FL10000
         IF ((GNDALT.GT.0.).AND.(ZMDL(K).LT.6.0)) THEN                   FL10010
            ASC = 6./(6.-GNDALT)                                         FL10020
            CON = -ASC*GNDALT                                            FL10030
            ZSC = ASC*ZMDL(K)+CON                                        FL10040
            IF (ZSC.LT.0.) ZSC = 0.                                      FL10050
         ENDIF                                                           FL10060
         ZGN(K) = ZSC                                                    FL10070
C                                                                        FL10080
C                                                                        FL10090
         ICLDS = ICLD1                                                   FL10100
         IF (ICLD1.EQ.0) ICLD1 = ICLD                                    FL10110
         ICLDL = ICLD1                                                   FL10120
         IF (ICLD1.GT.11) ICLD1 = 0                                      FL10130
         IF (IHAZ1.NE.0) IVUL1 = 0                                        FL10140
         IF (IHAZ1.NE.0) ICLD1 = 0                                        FL10150
         IF (ICLD1.NE.0) IVUL1 = 0                                       FL10160
         IF ((AHAZE.NE.0.).OR.(EQLWCZ.NE.0.)) GO TO 100                  FL10170
         IF (RRATZ.NE.0.) GO TO 100                                      FL10180
         IF ((IVSA.EQ.1).AND.(ICLD1.EQ.0)) THEN                          FL10190
            CALL LAYVSA (K,RH,AHAZE,IHAZ1,ZSTF)                           FL10200
         ELSE                                                            FL10210
            CALL LAYCLD (K,EQLWCZ,RRATZ,IAERSL,ICLD1,GNDALT,RAINRT)      FL10220
            IF (ICLD1.LT.1) GO TO 100                                    FL10230
            IF (ICLD1.GT.10) GO TO 100                                   FL10240
            IF (ZMDL(K).GT.CLDTOP(ICLD1)+GNDALT) THEN                    FL10250
               RRATZ = 0.                                                FL10260
            ENDIF                                                        FL10270
         ENDIF                                                           FL10280
  100    CONTINUE                                                        FL10290
         ICLDC = ICLD                                                    FL10300
         IF (ICLDS.NE.0) ICLDC = ICLDS                                   FL10310
         IF (ICLDC.EQ.18.OR.ICLDC.EQ.19) THEN                            FL10320
            DENSTY(16,K) = 0.                                            FL10330
            IF (ZMDL(K).GE.CLD1.AND.ZMDL(K).LE.CLD2) DENSTY(16,K) = CEXT FL10340
         ENDIF                                                           FL10350
         CLDAMT(K) = EQLWCZ                                              FL10360
         IF (ICLDS.EQ.0.AND.CLDAMT(K).EQ.0.) ICLD1 = 0                   FL10370
         RRAMT(K) = RRATZ                                                FL10380
         IF (MODEL.NE.0) THEN                                            FL10390
            IF (EQLWCZ.GT.0.0) RH = 100.0                                FL10400
            IF (RRATZ.GT.0.0) RH = 100.0                                 FL10410
         ENDIF                                                           FL10420
         AHAST(K) = AHAZE                                                FL10430
C                                                                        FL10440
C        IHAZ1  IS IHAZE FOR THIS LAYER                                   FL10450
C        ISEA1 IS ISEASN FOR THIS LAYER                                  FL10460
C        IVUL1 IS IVULCN FOR THE LAYER                                   FL10470
C                                                                        FL10480
         IF (ISEA1.EQ.0) ISEA1 = ISEASN                                  FL10490
         ITYAER = IHAZE                                                  FL10500
         IF (IHAZ1.GT.0) ITYAER = IHAZ1                                    FL10510
         IF (IVUL1.GT.0) IVULCN = IVUL1                                  FL10520
         IF (IVUL1.LE.0) IVUL1 = IVULCN                                  FL10530
C                                                                        FL10540
         IF (K.EQ.1) GO TO 130                                           FL10550
         IF (ICHR.EQ.1) GO TO 120                                        FL10560
         IF (ICLD1.NE.IREGC(IC1)) GO TO 110                              FL10570
         IF (IHAZ1.EQ.0.AND.ICLD1.EQ.0) THEN                              FL10580
            IF (ZSC.GT.2.) ITYAER = 6                                    FL10590
            IF (ZSC.GT.10.) ITYAER = IVULCN+10                           FL10600
            IF (ZSC.GT.30.) ITYAER = 19                                  FL10610
            IF (ITYAER.EQ.ICH(IC1)) GO TO 130                            FL10620
         ENDIF                                                           FL10630
         IF (ICLD1.EQ.0.AND.IHAZ1.EQ.0) GO TO 120                         FL10640
         N = 7                                                           FL10650
         IF (IC1.GT.1) N = IC1+10                                        FL10660
         IF (IHAZ1.EQ.0) GO TO 130                                        FL10670
         IF (IHAZ1.NE.ICH(IC1)) GO TO 120                                 FL10680
         GO TO 130                                                       FL10690
  110    IF (ICLD1.NE.0) THEN                                            FL10700
            IF (ICLD1.EQ.IREGC(1)) THEN                                  FL10710
               N = 7                                                     FL10720
               ALTB(1) = ZMDL(K)                                         FL10730
               GO TO 140                                                 FL10740
            ENDIF                                                        FL10750
            IF (IC1.EQ.1) GO TO 120                                      FL10760
            IF (ICLD1.EQ.IREGC(2)) THEN                                  FL10770
               N = 12                                                    FL10780
               ALTB(2) = ZMDL(K)                                         FL10790
               GO TO 140                                                 FL10800
            ENDIF                                                        FL10810
            IF (IC1.EQ.2) GO TO 120                                      FL10820
            IF (ICLD1.EQ.IREGC(3)) THEN                                  FL10830
               N = 13                                                    FL10840
               ALTB(3) = ZMDL(K)                                         FL10850
               GO TO 140                                                 FL10860
            ENDIF                                                        FL10870
         ELSE                                                            FL10880
            IF (IHAZ1.EQ.0.AND.ICLD1.EQ.0) THEN                           FL10890
               IF (ZSC.GT.2.) ITYAER = 6                                 FL10900
               IF (ZSC.GT.10.) ITYAER = IVULCN+10                        FL10910
               IF (ZSC.GT.30.) ITYAER = 19                               FL10920
            ENDIF                                                        FL10930
            IF (ITYAER.EQ.ICH(1)) THEN                                   FL10940
               N = 7                                                     FL10950
               ALTB(1) = ZMDL(K)                                         FL10960
               GO TO 140                                                 FL10970
            ENDIF                                                        FL10980
            IF (IC1.EQ.1) GO TO 120                                      FL10990
            IF (ITYAER.EQ.ICH(2)) THEN                                   FL11000
               N = 12                                                    FL11010
               ALTB(2) = ZMDL(K)                                         FL11020
               GO TO 140                                                 FL11030
            ENDIF                                                        FL11040
            IF (IC1.EQ.2) GO TO 120                                      FL11050
            IF (ITYAER.EQ.ICH(3)) THEN                                   FL11060
               N = 13                                                    FL11070
               ALTB(3) = ZMDL(K)                                         FL11080
               GO TO 140                                                 FL11090
            ENDIF                                                        FL11100
         ENDIF                                                           FL11110
  120    IC1 = IC1+1                                                     FL11120
         ICL = 0                                                         FL11130
C                                                                        FL11140
C                                                                        FL11150
C                                                                        FL11160
         N = IC1+10                                                      FL11170
         IF (RH.GT.0.) RHH = RH                                          FL11180
         IF (IC1.LE.4) GO TO 130                                         FL11190
         IC1 = 4                                                         FL11200
         N = 14                                                          FL11210
         ITYAER = ICH(IC1)                                               FL11220
  130    ICH(IC1) = ITYAER                                               FL11230
         IREGC(IC1) = ICLD1                                              FL11240
         ALTB(IC1) = ZMDL(K)                                             FL11250
C                                                                        FL11260
C        FOR LVSA OR CLD OR RAIN ONLY                                    FL11270
C                                                                        FL11280
  140    IF (IHAZ1.LE.0) IHAZ1 = IHAZE                                     FL11290
C                                                                        FL11300
         DENSTY(7,K) = 0.                                                FL11310
         DENSTY(12,K) = 0.                                               FL11320
         DENSTY(13,K) = 0.                                               FL11330
         DENSTY(14,K) = 0.                                               FL11340
         DENSTY(15,K) = 0.                                               FL11350
C                                                                        FL11360
C        IF((GNDALT.GT.0.).AND.(ZMDL(K).LT.6.0)) THEN                    FL11370
C        J= INT(ZSC+1.0E-6)+1                                            FL11380
C        FAC=ZSC- REAL(J-1)                                              FL11390
C        ELSE                                                            FL11400
C                                                                        FL11410
         J =  INT(ZMDL(K)+1.0E-6)+1                                      FL11420
         IF (ZMDL(K).GE.25.0) J = (ZMDL(K)-25.0)/5.0+26.                 FL11430
         IF (ZMDL(K).GE.50.0) J = (ZMDL(K)-50.0)/20.0+31.                FL11440
         IF (ZMDL(K).GE.70.0) J = (ZMDL(K)-70.0)/30.0+32.                FL11450
         J = MIN(J,I_32)                                                   FL11460
         FAC = ZMDL(K)- REAL(J-1)                                        FL11470
         IF (J.LT.26) GO TO 150                                          FL11480
         FAC = (ZMDL(K)-5.0* REAL(J-26)-25.)/5.                          FL11490
         IF (J.GE.31) FAC = (ZMDL(K)-50.0)/20.                           FL11500
         IF (J.GE.32) FAC = (ZMDL(K)-70.0)/30.                           FL11510
         FAC = MIN(FAC,1.0)                                              FL11520
C                                                                        FL11530
C        ENDIF                                                           FL11540
C                                                                        FL11550
  150    L = J+1                                                         FL11560
         WHN = 0.                                                        FL11570
         IF (MODEL.EQ.0) THEN                                            FL11580
            CALL GETPT (K,ZMDL,P,T,WHN,INEW)                             FL11590
            WH(K) = WHN                                                  FL11600
         ELSE                                                            FL11610
            CALL CHECK (P(K),JUNIT(1),1)                                 FL11620
            CALL CHECK (T(K),JUNIT(2),2)                                 FL11630
            CALL LDEFAL (ZMDL(K),P(K),T(K))                              FL11640
            CALL LCONVR (P(K),T(K))                                      FL11650
            WH(K) = WMOL(1)                                              FL11660
         ENDIF                                                           FL11670
C                                                                        FL11680
         TMP = T(K)-T0                                                   FL11690
C                                                                        FL11700
C        FOR LVSA OR CLD OR RAIN ONLY                                    FL11710
C                                                                        FL11720
         IF (RH.GT.0.0) THEN                                             FL11730
            TA = T0/T(K)                                                 FL11740
            WH(K) = F(TA)*0.01*RH                                        FL11750
            IF (IVSA.EQ.1) GO TO 160                                     FL11760
C                                                                        FL11770
C           WRITE(IPR,800) ZMDL(K),EQLWCZ,ICLD1,RRATZ                    FL11780
C                                                                        FL11790
C                                                                        FL11800
  160    ENDIF                                                           FL11810
C                                                                        FL11820
C        C    IF (M3.GT.0) WO(K,7)=WO(J,M3)*(WO(L,M3)/WO(J,M3))**FAC     FL11830
C                                                                        FL11840
         HSTOR(K) = 0.                                                   FL11850
C                                                                        FL11860
C        IF (HMIX(J).LE.0.) GO TO 40                                     FL11870
C        IF (HMIX(L).LE.0.) GO TO 40                                     FL11880
C        HSTOR(K)=HMIX(J)*(HMIX(L)/HMIX(J))**FAC                         FL11890
C                                                                        FL11900
         DENSTY(7,K) = 0.                                                FL11910
         DENSTY(12,K) = 0.                                               FL11920
         DENSTY(13,K) = 0.                                               FL11930
         DENSTY(14,K) = 0.                                               FL11940
         DENSTY(15,K) = 0.                                               FL11950
C                                                                        FL11960
C        PS=P(K,7)/1013.0                                                FL11970
C                                                                        FL11980
         TS = 273.15/T(K)                                                FL11990
         WTEMP = WH(K)                                                   FL12000
         RELHUM(K) = 0.                                                  FL12010
         IF (WTEMP.LE.0.) GO TO 170                                      FL12020
         RELHUM(K) = 100.0*WTEMP/F(TS)                                   FL12030
         IF (RELHUM(K).GT.100.) WRITE (IPR,920) RELHUM(K),ZMDL(K)        FL12040
         IF (RELHUM(K).GT.100.) RELHUM(K) = 100.                         FL12050
         IF (RELHUM(K).LT.0.) WRITE (IPR,920) RELHUM(K),ZMDL(K)          FL12060
         IF (RELHUM(K).LT.0.) RELHUM(K) = 0.                             FL12070
  170    RHH = RELHUM(K)                                                 FL12080
         RH = RHH                                                        FL12090
         IF (VIS1.LE.0.0) VIS1 = VIS                                     FL12100
         IF (AHAZE.EQ.0.0) GO TO 180                                     FL12110
         DENSTY(N,K) = AHAZE                                             FL12120
         IF (ITYAER.EQ.3) GO TO 180                                      FL12130
C                                                                        FL12140
C        AHAZE IS IN LOWTRAN NUMBER DENSTY UNITS                         FL12150
C                                                                        FL12160
         GO TO 200                                                       FL12170
  180    CONTINUE                                                        FL12180
C                                                                        FL12190
C        AHAZE NOT INPUT OR NAVY MARITIME MODEL IS CALLED                FL12200
C                                                                        FL12210
C        CHECK IF GNDALT NOT ZERO                                        FL12220
C                                                                        FL12230
         IF ((GNDALT.GT.0.).AND.(ZMDL(K).LT.6.0)) THEN                   FL12240
            J =  INT(ZSC+1.0E-6)+1                                       FL12250
            FAC = ZSC- REAL(J-1)                                         FL12260
            L = J+1                                                      FL12270
         ENDIF                                                           FL12280
         IF (ITYAER.EQ.3.AND.ICL.EQ.0) THEN                              FL12290
            CALL MARINE (VIS1,MODEL,WSS,WHH,ICSTL,EXTC,ABSC,IC1)         FL12300
            IREG(IC1) = 1                                                FL12310
            VIS = VIS1                                                   FL12320
            ICL = ICL+1                                                  FL12330
         ENDIF                                                           FL12340
         IF (ITYAER.EQ.10.AND.IDSR.EQ.0) THEN                            FL12350
            CALL DESATT (WSS,VIS1)                                       FL12360
            IREG(IC1) = 1                                                FL12370
            VIS = VIS1                                                   FL12380
            IDSR = IDSR+1                                                FL12390
         ENDIF                                                           FL12400
         IF (AHAZE.GT.0.0) GO TO 200                                     FL12410
         CALL AERPRF (J,K,VIS1,HAZ1,IHAZ1,ICLD1,ISEA1,IVUL1,NN)           FL12420
         CALL AERPRF (L,K,VIS1,HAZ2,IHAZ1,ICLD1,ISEA1,IVUL1,NN)           FL12430
         HAZE = 0.                                                       FL12440
         IF ((HAZ1.LE.0.0).OR.(HAZ2.LE.0.0)) GO TO 190                   FL12450
         HAZE = HAZ1*(HAZ2/HAZ1)**FAC                                    FL12460
         DENSTY(N,K) = HAZE                                              FL12470
  190    CONTINUE                                                        FL12480
         IF (CLDAMT(K).GT.0.0) THEN                                      FL12490
            HAZE = HAZEC(K)                                              FL12500
            IF (HAZE.GT.0.) DENSTY(N,K) = HAZE                           FL12510
         ENDIF                                                           FL12520
  200    CONTINUE                                                        FL12530
         IF (K.EQ.1) GO TO 210                                           FL12540
         IF (CLDAMT(K).LE.0.0.AND.CLDAMT(K-1).GT.0.0) THEN               FL12550
            HAZE = HAZ1*(HAZ2/HAZ1)**FAC                                 FL12560
            IF (HAZE.GT.0.) DENSTY(N,K) = HAZE                           FL12570
         ENDIF                                                           FL12580
  210    CONTINUE                                                        FL12590
         ITY1(K) = ITYAER                                                FL12600
         IH1(K) = IHAZ1                                                   FL12610
         IF (AHAZE.NE.0) IH1(K) = -99                                    FL12620
         IS1(K) = ISEA1                                                  FL12630
         IVL1(K) = IVUL1                                                 FL12640
         WGM(K) = WH(K)                                                  FL12650
  220 CONTINUE                                                           FL12660
C                                                                        FL12670
C     END OF LOOP                                                        FL12680
C                                                                        FL12690
      IF (ML.LT.20) WRITE (IPR,925)                                      FL12700
      IF (ML.GE.20) WRITE (IPR,930)                                      FL12710
      IHH = ICLD                                                         FL12720
      IF (IHH.LE.0) IHH = 12                                             FL12730
      IHH = MIN(IHH,I_12)                                                  FL12740
      IF (ICLD.EQ.18) IHH = 13                                           FL12750
      IF (ICLD.EQ.19) IHH = 14                                           FL12760
      IF (ICLD.EQ.20) IHH = 15                                           FL12770
C                                                                        FL12780
      HHOL = AHAHOL(IHH)                                                 FL12790
      IF (IVSA.NE.0) HHOL = AHLVSA                                       FL12800
C                                                                        FL12810
      IF (ICLD.NE.0) THEN                                                FL12820
         IF (JPRT.EQ.0) WRITE (IPR,935) HHOL,M                           FL12830
      ENDIF                                                              FL12840
      IF (JPRT.EQ.0) WRITE (IPR,940)                                     FL12850
C                                                                        FL12860
      IF (JPRT.EQ.1) GO TO 240                                           FL12870
c
      DO 230 KK = 1, ML                                                  FL12880
c
         K = KK                                                          FL12930
         IF (JPRT.EQ.1) GO TO 230                                        FL12940
C                                                                        FL12950
         AHOL1 = BLANK                                                   FL12960
         AHOL2 = BLANK                                                   FL12970
         AHOL3 = BLANK                                                   FL12980
         ITYAER = ITY1(KK)                                               FL12990
         IF (ITYAER.EQ.0) ITYAER = 1                                     FL13000
         IF (ITYAER.EQ.16) ITYAER = 11                                   FL13010
         IF (ITYAER.EQ.17) ITYAER = 11                                   FL13020
         IF (ITYAER.EQ.18) ITYAER = 13                                   FL13030
         IHAZ1 = IH1(KK)                                                  FL13040
         ISEA1 = IS1(KK)                                                 FL13050
         IVUL1 = IVL1(KK)                                                FL13060
C                                                                        FL13070
         AHOL1 = HZ(ITYAER)                                              FL13080
         IF (IVSA.EQ.1) AHOL1 = HHOL                                     FL13090
         IF (CLDAMT(KK).GT.0.0.OR.RRAMT(KK).GT.0.0) AHOL1 = HHOL         FL13100
         IF (IHAZE.EQ.0) AHOL1 = HHOL                                    FL13110
         AHOL2 = AHUS                                                    FL13120
         IF (AHAST(KK).EQ.0) AHOL2 = AHOL1                               FL13130
         IF (CLDAMT(KK).GT.0.0.OR.RRAMT(KK).GT.0.0) AHOL2 = HHOL         FL13140
         IF (ZGN(KK).GT.2.0) AHOL3 = SEASN(ISEA1)                        FL13150
         WRITE (IPR,945) ZMDL(KK),P(KK),T(KK),RELHUM(KK),WH(KK),         FL13160
     *        CLDAMT(KK),RRAMT(KK),AHOL1,AHOL2,AHOL3                       FL13170
 230  continue
 240  IMMAX = ML                                                         FL13180
      M = 7                                                              FL13190
      IF (ML.EQ.1) WRITE (IPR,925)                                       FL13200
      IF (ML.NE.1) MODEL = M                                             FL13210
      RETURN                                                             FL13220
C                                                                        FL13230
  900 FORMAT('   ERROR MODEL EQ 0 AND ARMY MODEL CANNOT MIX')            FL13240
  905 FORMAT('  ERROR ML GT 24 AND ARMY MODEL TOP LAYER TRUNCATED')      FL13250
  910 FORMAT(/,10X,' MODEL 0 / 7 USER INPUT DATA ',//)                   FL13260
  915      FORMAT    (4F10.3,5I5)                                        FL13270
  920 FORMAT(' ***ERROR RELHUM ' ,E15.4,'  AT ALT  ',F12.3)              FL13280
  925 FORMAT('0 ')                                                       FL13290
  930 FORMAT('1  ')                                                      FL13300
  935 FORMAT(//'0 CLOUD AND OR RAIN TYPE CHOSEN IS   ',A20,              FL13310
     * '  M IS SET TO',I5//)                                             FL13320
  940 FORMAT(T7,'Z',T17,'P',T26,'T',T32,'REL H', T41,'H2O',              FL13330
     * T49,'CLD AMT',T59,'RAIN RATE', T90,'AEROSOL'/,                    FL13340
     * T6,'(KM)',T16,'(MB)',T25,'(K)',T33,'(%)',T39,'(GM M-3)',T49,      FL13350
     * '(GM M-3)',T59,'(MM HR-1)',T69,                                   FL13360
     * 'TYPE', T90,'PROFILE')                                            FL13370
  945 FORMAT(2F10.3,2F8.2,1P3E10.3,1X,3A20)                              FL13380
  950 FORMAT(//,' MODEL ATMOSPHERE NO. ',I5,' ICLD =',I5,//)             FL13390
C                                                                        FL13400
      END                                                                FL13410
C
C     ***********************************************************
C
      SUBROUTINE LAYCLD(K,CLDATZ,RRATZ,IAERSL,ICLD1,GNDALT,RAINRT)       FL13420
C                                                                        FL13430
C     THIS SUBROUTINE RESTRUCTURES THE ATMOSPHERIC PROFILE               FL13440
C     TO PROFIDE FINER LAYERING WITHIN THE FIRST 6 KM.                   FL13450
C                                                                        FL13460
C     ZMDL COMMON /MODEL/ FINAL ALTITUDE FOR LOWTRAN                     FL13470
C     ZK  EFFECTIVE CLOUD ALTITUDES                                      FL13480
C     ZCLD CLOUD ALTITUDE ARRAY                                          FL13490
C     ZDIF  ALT DIFF OF 2 LAYERS                                         FL13500
C     ZDA COMMON /MDATA/ CLD AND RAIN INFO IN THIS COMMON                FL13510
C                                                                        FL13520
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
c     COMMON /MDATA/ ZDA(MXZMD),P(MXZMD),T(MXZMD),WH(MXZMD),WO(MXZMD),   FL13530
c    *     HMIX(MXZMD),CLD(MXZMD,7),RR(MXZMD,7)                          FL13540
      COMMON /MDATA/                              WH(MXZMD),WO(MXZMD),   FL13530
     *                 CLD(MXZMD,7),RR(MXZMD,7)                          FL13540
      COMMON /MDATA2/ZDA(MXZMD),P(MXZMD),T(MXZMD)
      COMMON /MODEL/ ZMDL(MXZMD),PN(MXZMD),TN(MXZMD),                    FL13550
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL13560
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)
      DIMENSION ZCLD(16)                                                 FL13570
      DATA ZCLD/ 0.0,0.16,0.33,0.66,1.0,1.5,2.0,2.4,2.7,                 FL13580
     * 3.0,3.5,4.0,4.5,5.0,5.5,6.0/                                      FL13590
      DATA CLDTP/6.0001/                                                 FL13600
      DATA DELZ /0.02/                                                   FL13610
      ICLD = ICLD1                                                       FL13620
      IF (ICLD.EQ.0) RETURN                                              FL13630
      IF (ICLD.GT.11) RETURN                                             FL13640
      ZK = ZMDL(K)-GNDALT                                                FL13650
      ZK = MAX(ZK,0.)                                                    FL13660
      IF (ZMDL(K).GT.6.) ZK = ZMDL(K)                                    FL13670
      IF (ICLD.GT.5) GO TO 10                                            FL13680
C                                                                        FL13690
C     CC                                                                 FL13700
C     CC    ICLD  IS  1- 5 ONE OF 5 SPECIFIC CLOUD MODELS IS CHOSEN      FL13710
C     CC                                                                 FL13720
C                                                                        FL13730
      MC = ICLD                                                          FL13740
      MR = 6                                                             FL13750
      GO TO 20                                                           FL13760
   10 CONTINUE                                                           FL13770
C                                                                        FL13780
C     CC                                                                 FL13790
C     CC   ICLD  IS  6-10 ONE OF 5 SPECIFIC CLOUD/RAIN MODELS CHOSEN     FL13800
C     CC                                                                 FL13810
C                                                                        FL13820
      IF (ICLD.EQ.6) MC = 3                                              FL13830
      IF (ICLD.EQ.7.OR.ICLD.EQ.8) MC = 5                                 FL13840
      IF (ICLD.GT.8) MC = 1                                              FL13850
      MR = ICLD-5                                                        FL13860
   20 CONTINUE                                                           FL13870
      IF (ZK.GT.CLDTP) GO TO 70                                          FL13880
      CLDATZ = 0.                                                        FL13890
      RRATZ = 0.                                                         FL13900
      IF (ZK.LE.10.) RRATZ = RAINRT                                      FL13910
      IF (MC.LT.1) GO TO 60                                              FL13920
      DO 50 MK = 1, 15                                                   FL13930
         IF (ZK.GE.ZCLD(MK+1)) GO TO 50                                  FL13940
         IF (ZK.LT.ZCLD(MK)) GO TO 50                                    FL13950
         IF (ABS(ZK-ZCLD(MK)).LT.DELZ) GO TO 30                          FL13960
         GO TO 40                                                        FL13970
   30    CLDATZ = CLD(MK,MC)                                             FL13980
         RRATZ = RR(MK,MR)                                               FL13990
         GO TO 60                                                        FL14000
   40    ZDIF = ZCLD(MK+1)-ZCLD(MK)                                      FL14010
         IF (ZDIF.LT.DELZ) GO TO 30                                      FL14020
         FAC = (ZCLD(MK+1)-ZK)/ZDIF                                      FL14030
         CLDATZ = CLD(MK+1,MC)+FAC*(CLD(MK,MC)-CLD(MK+1,MC))             FL14040
         RRATZ = RR(MK+1,MR)+FAC*(RR(MK,MR)-RR(MK+1,MR))                 FL14050
         GO TO 60                                                        FL14060
   50 CONTINUE                                                           FL14070
   60 CLDAMT(K) = CLDATZ                                                 FL14080
      CLD(K,7) = CLDATZ                                                  FL14090
      RR(K,7) = RRATZ                                                    FL14100
      RRAMT(K) = RRATZ                                                   FL14110
      RETURN                                                             FL14120
   70 CONTINUE                                                           FL14130
      CLDAMT(K) = 0.0                                                    FL14140
      RRAMT(K) = 0.0                                                     FL14150
      CLDATZ = 0.0                                                       FL14160
      RRATZ = 0.0                                                        FL14170
      RETURN                                                             FL14180
      END                                                                FL14190
      BLOCK DATA MDTA                                                    FL14200
C                                                                        FL14210
C     >    BLOCK DATA                                                    FL14220
C                                                                        FL14230
C     CLOUD AND RAIN   DATA                                              FL14240
C                                                                        FL14250
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                    FL14260
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)     FL14270
C                                                                        FL14280
c     COMMON /MDATA/ ZDA(MXZMD),P(MXZMD),T(MXZMD),WH(MXZMD),WO(MXZMD),   FL14290
c    *     HMIX(MXZMD),CLD1(MXZMD),CLD2(MXZMD),CLD3(MXZMD),CLD4(MXZMD),  FL14300
c    *     CLD5(MXZMD),CLD6(MXZMD),CLD7(MXZMD),RR1(MXZMD),RR2(MXZMD),    FL14310
c    *     RR3(MXZMD),RR4(MXZMD),RR5(MXZMD),RR6(MXZMD),RR7(MXZMD)        FL14320
      COMMON /MDATA/                              WH(MXZMD),WO(MXZMD),   FL14290
     *                 CLD1(MXZMD),CLD2(MXZMD),CLD3(MXZMD),CLD4(MXZMD),  FL14300
     *     CLD5(MXZMD),CLD6(MXZMD),CLD7(MXZMD),RR1(MXZMD),RR2(MXZMD),    FL14310
     *     RR3(MXZMD),RR4(MXZMD),RR5(MXZMD),RR6(MXZMD),RR7(MXZMD)        FL14320
      COMMON /MDATA2/ZDA(MXZMD),P(MXZMD),T(MXZMD)
C                                                                        FL14330
C     DATA  Z/                                                           FL14340
C     C       0.0,       1.0,       2.0,       3.0,       4.0,           FL14350
C     C       5.0,       6.0,       7.0,       8.0,       9.0,           FL14360
C     C      10.0,      11.0,      12.0,      13.0,      14.0,           FL14370
C     C      15.0,      16.0,      17.0,      18.0,      19.0,           FL14380
C     C      20.0,      21.0,      22.0,      23.0,      24.0,           FL14390
C     C      25.0,      27.5,      30.0,      32.5,      35.0,           FL14400
C     C      37.5,      40.0,      42.5,      45.0,      47.5,           FL14410
C     C      50.0,      55.0,      60.0,      65.0,      70.0,           FL14420
C     C      75.0,      80.0,      85.0,      90.0,      95.0,           FL14430
C     C     100.0,     105.0,     110.0,     115.0,     120.0/           FL14440
C     CC   CLOUD MODELS 1-5                                              FL14450
C                                                                        FL14460
      DATA CLD1/ 0.0,0.0,0.0,0.2,0.35,1.0,1.0,1.0,0.3,0.15,5990*0.0/     FL14470
      DATA CLD2/ 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.3,0.4,0.3,5990*0.0/       FL14480
      DATA CLD3/ 0.0,0.0,0.15,0.30,0.15,5995*0.0/                        FL14490
      DATA CLD4/ 0.0,0.0,0.0,0.10,0.15,0.15,0.10,5993*0.0/               FL14500
      DATA CLD5/ 0.0,0.30,0.65,0.40,5996*0.0/                            FL14510
      DATA CLD6/ 6000*0.0/                                               FL14520
      DATA CLD7/ 6000*0.0/                                               FL14530
C                                                                        FL14540
C     CC   RAIN MODELS 1-5                                               FL14550
C                                                                        FL14560
      DATA RR1/ 2.0,1.78,1.43,1.22,0.86,0.22,5994*0.0/                   FL14570
      DATA RR2/ 5.0,4.0,3.4,2.6,0.8,0.2,5994*0.0/                        FL14580
      DATA RR3/ 12.5,10.5,8.0,6.0,2.5,0.8,0.2,5993*0.0/                  FL14590
      DATA RR4/ 25.0,21.5,17.5,12.0,7.5,4.2,2.5,1.0,0.7,0.2,5990*0.0/    FL14600
      DATA RR5/ 75.0,70.0,65.0,60.0,45.0,20.0,12.5,7.0,3.5,              FL14610
     * 1.0,0.2,5989*0.0/                                                 FL14620
      DATA RR6/ 6000*0.0/                                                FL14630
      DATA RR7/ 6000*0.0/                                                FL14640
C                                                                        FL14650
C     DATA CO2       /                                                   FL14660
C                                                                        FL14670
      END                                                                FL14680
C
C     **************************************************************
C
      SUBROUTINE GETPT(K,ZMDL,P,T,WHN,INEW)                              FL14690
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),VH(16),TX(16),W(16)       
      COMMON WPATH(IM2,16),TBBY(IM2)                                     FL14720
      COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),VX2(47),AWCCON(5)
C
      CHARACTER*8      HMOD                                              FL14730
C
      COMMON /CMN/ HMOD(3),Zc(MXZMD),Pc(MXZMD),Tc(MXZMD),RFNDXM(MXZMD), 
     *          ZP(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),
     *          TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),
     *          AMTP(MXMOL,MXPDIM)                                        
C
C
      COMMON /PATHD/ PBAR(MXLAY),TBAR(MXLAY),AMOUNT(MXMOL,MXLAY),        FA00530
     *               WN2L(MXLAY),DVL(MXLAY),WTOTL(MXLAY),ALBL(MXLAY),    FA00540
     *               ADBL(MXLAY),AVBL(MXLAY),H2OSL(MXLAY),IPATH(MXLAY),  FA00550
     *               ITYL(MXLAY),SECNTA(MXLAY),HT1,HT2,ALTZ(0:MXLAY),    FA00560
     *               PM(0:MXLAY),TM(0:MXLAY)                             FA00570
c 
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,N_LVL,JH1              FL14780
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2,GRAV,CPDAIR,AIRMWT,SECDY 
      DIMENSION INEW( *)                                                 FL14790
      DIMENSION ZMDL( *),P(MXZMD),T(MXZMD)                               FL14800
C                                                                        FL14810
C     ZP BLANK COMMON UNUSED                                             FL14820
C     ALTZ  BLANK COMMON LBLRTM ALTITUDES                                  FL14830
C     ZMDL COMMON /MODEL/ FINAL ALTITUDE FOR LOWTRAN                     FL14840
C                                                                        FL14850
C     THIS ROUTINE INTERPOLATES P,T,AND H2O INTO                         FL14860
C     LOWTRAN LAYERS WHEN MODEL = 7                                      FL14870
C                                                                        FL14880
      WTH2O = 18.015 
      B = WTH2O*1.E06/AVOGAD                                             FL14910
      J = INEW(K)                                                        FL14920
C                                                                        FL14930
      JL = J-1                                                           FL14940
      IF (JL.LT.0) JL = 0                                                FL14950
      JP = JL+1                                                          FL14960
      IF (JP.GT.ML) GO TO 40                                             FL14970
      DIF = ALTZ(JP)-ALTZ(JL)                                                FL14980
      FAC = (ZMDL(K)-ALTZ(JL))/DIF                                         FL14990
      P(K) = PM(JL)                                                      FL15000
      IF (PM(JP).LE.0.0.OR.PM(JL).LE.0.) GO TO 10                        FL15010
      P(K) = PM(JL)*(PM(JP)/PM(JL))**FAC                                 FL15020
   10 T(K) = TM(JL)                                                      FL15030
      IF (TM(JP).LE.0.0.OR.TM(JL).LE.0.) GO TO 20                        FL15040
      T(K) = TM(JL)*(TM(JP)/TM(JL))**FAC                                 FL15050
   20 WHN = DENW(JL+1)                                                     FL15060
      IF (DENW(JP+1).LE.0.0.OR.DENW(JL+1).LE.0.) GO TO 30                    FL15070
      WHN = DENW(JL+1)*(DENW(JP+1)/DENW(JL+1))**FAC                            FL15080
   30 CONTINUE                                                           FL15090
      WHN = WHN*B                                                        FL15100
      RETURN                                                             FL15110
   40 P(K) = PM(JL)                                                      FL15120
      T(K) = TM(JL)                                                      FL15130
      WHN = DENW(JL)*B                                                   FL15140
      RETURN                                                             FL15150
      END                                                                FL15160
      SUBROUTINE CIRR18                                                  FL15170
C                                                                        FL15180
C     ****************************************************************** FL15190
C     *  ROUTINE TO SET CTHIK CALT CEXT  FOR  CIRRUS CLOUDS 18 19        FL15200
C     *  INPUTS!                                                         FL15210
C     *           CHTIK    -  CIRRUS THICKNESS (KM)                      FL15220
C     *                       0 = USE THICKNESS STATISTICS               FL15230
C     *                       .NE. 0 = USER DEFINES THICKNESS            FL15240
C     *                                                                  FL15250
C     *           CALT     -  CIRRUS BASE ALTITUDE (KM)                  FL15260
C     *                       0 = USE CALCULATED VALUE                   FL15270
C     *                       .NE. 0 = USER DEFINES BASE ALTITUDE        FL15280
C     *                                                                  FL15290
C     *           ICLD     -  CIRRUS PRESENCE FLAG                       FL15300
C     *                       0 = NO CIRRUS                              FL15310
C     *                       18  19 = USE CIRRUS PROFILE                FL15320
C     *                                                                  FL15330
C     *           MODEL    -  ATMOSPHERIC MODEL                          FL15340
C     *                       1-5  AS IN MAIN PROGRAM                    FL15350
C     *                       MODEL = 0,6,7 NOT USED SET TO 2            FL15360
C     *                                                                  FL15370
C     *  OUTPUTS!                                                        FL15380
C     *         CTHIK        -  CIRRUS THICKNESS (KM)                    FL15390
C     *         CALT         -  CIRRUS BASE ALTITUDE (KM)                FL15400
C     CEXT IS THE EXTINCTION COEFFIENT(KM-1) AT 0.55                     FL15410
C     DEFAULT VALUE 0.14*CTHIK                                           FL15420
C     *                                                                  FL15430
C     ****************************************************************** FL15440
C                                                                        FL15450
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),VH(16),TX(16),W(16)       FL15460
      COMMON WPATH(IM2,16),TBBY(IM2)                                     FL15470
      COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),VX2(47),AWCCON(5)          FL15480
C
      COMMON /LCRD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB   FL15490
      COMMON /LCRD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,    FL15500
     *     RAINRT                                                        FL15510
      COMMON /LCRD2A/ CTHIK,CALT,CEXT                                    FL15520
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,IMULT,JH1            FL15530
      COMMON /LCRD4/ V1,V2,DV                                            FL15540
      COMMON/MODEL/ ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),                     FL15550
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL15560
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL15570
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL15580
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL15590
      DIMENSION CBASE(5,2),TSTAT(11),PTAB(5),CAMEAN(5)                   FL15600
      DIMENSION CBASE1(5),CBASE2(5)                                      FL15610
      EQUIVALENCE (CBASE1(1),CBASE(1,1)),(CBASE2(1),CBASE(1,2))          FL15620
C                                                                        FL15630
      DATA  CAMEAN           / 11.0, 10.0, 8.0, 7.0, 5.0 /               FL15640
      DATA  PTAB           / 0.8, 0.4, 0.5, 0.45, 0.4/                   FL15650
      DATA  CBASE1            / 7.5, 7.3, 4.5, 4.5, 2.5 /                FL15660
      DATA  CBASE2            /16.5,13.5,14.0, 9.5,10.0 /                FL15670
      DATA  TSTAT             / 0.0,.291,.509,.655,.764,.837,.892,       FL15680
     * 0.928, 0.960, 0.982, 1.00 /                                       FL15690
      MDL = MODEL                                                        FL15700
C                                                                        FL15710
C     CHECK IF USER WANTS TO USE A THICKNESS VALUE HE PROVIDES           FL15720
C     DEFAULTED MEAN CIRRUS THICKNESS IS 1.0KM  OR 0.2 KM.               FL15730
C                                                                        FL15740
      IF (CTHIK.GT.0.0) GO TO 10                                         FL15750
      IF (ICLD.EQ.18) CTHIK = 1.0                                        FL15760
      IF (ICLD.EQ.19) CTHIK = 0.2                                        FL15770
   10 IF (CEXT.EQ.0.) CEXT = 0.14*CTHIK                                  FL15780
C                                                                        FL15790
C     BASE HEIGHT CALCULATIONS                                           FL15800
C                                                                        FL15810
      IF (MODEL.LT.1.OR.MODEL.GT.5) MDL = 2                              FL15820
C                                                                        FL15830
      HMAX = CBASE(MDL,2)-CTHIK                                          FL15840
      BRANGE = HMAX-CBASE(MDL,1)                                         FL15850
      IF (CALT.GT.0.0) GO TO 20                                          FL15860
      CALT = CAMEAN(MDL)                                                 FL15870
C                                                                        FL15880
   20 IF (ICLD.EQ.18) WRITE (IPR,900)                                    FL15890
      IF (ICLD.EQ.19) WRITE (IPR,905)                                    FL15900
      WRITE (IPR,910) CTHIK                                              FL15910
      WRITE (IPR,915) CALT                                               FL15920
      WRITE (IPR,920) CEXT                                               FL15930
C                                                                        FL15940
C     END OF CIRRUS MODEL SET UP                                         FL15950
C                                                                        FL15960
      RETURN                                                             FL15970
C                                                                        FL15980
  900 FORMAT(15X,'CIRRUS ATTENUATION INCLUDED   (STANDARD CIRRUS)')      FL15990
  905 FORMAT(15X,'CIRRUS ATTENUATION INCLUDED   (THIN     CIRRUS)')      FL16000
  910 FORMAT(15X,'CIRRUS THICKNESS ',                                    FL16010
     * F10.3,'KM')                                                       FL16020
  915 FORMAT(15X,'CIRRUS BASE ALTITUDE ',                                FL16030
     * F10.3,' KM')                                                      FL16040
  920   FORMAT(15X,'CIRRUS PROFILE EXTINCT ',F10.3)                      FL16050
C                                                                        FL16060
      END                                                                FL16070
      SUBROUTINE DESATT(WSPD,VIS)                                        FL16080
C                                                                        FL16090
C     ****************************************************************** FL16100
C     *                                                                  FL16110
C     *    THIS SUBROUTINE CALCULATES THE ATTENUATION COEFFICIENTS AND   FL16120
C     *    ASYMMETRY PARAMETER FOR THE DESERT AEROSOL BASED ON THE WIND  FL16130
C     *    SPEED AND METEOROLOGICAL RANGE                                FL16140
C     *                                                                  FL16150
C     *                                                                  FL16160
C     *                                                                  FL16170
C     *    PROGRAMMED BY:  D. R. LONGTIN         OPTIMETRICS, INC.       FL16180
C     *                                          BURLINGTON, MASSACHUSET FL16190
C     *                                          JULY 1987               FL16200
C     *                                                                  FL16210
C     *                                                                  FL16220
C     *    INPUTS:    WSPD    -  WIND SPEED (IN M/S) AT 10 M             FL16230
C     *               VIS     -  METEOROLOGICAL RANGE (KM)               FL16240
C     *                                                                  FL16250
C     *    OUTPUTS:   DESEXT  -  EXTINCTION COEFFICIENT AT 47 WAVELENGTH FL16260
C     *               DESSCA  -  SCATTERING COEFFICIENT AT 47 WAVELENGTH FL16270
C     *    *****      DESABS  -  ABSORPTION COEFFICIENT AT 47 WAVELENGTH FL16280
C     *               DESG    -  ASYMMETRY PARAMETER AT 47 WAVELENGTHS   FL16290
C     *                                                                  FL16300
C     ****************************************************************** FL16310
C                                                                        FL16320
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELHUM(MXZMD),WHNO3(MXZMD),ICH(4),VH(16),TX(16),W(16)       FL16330
      COMMON WPATH(IM2,16),TBBY(IM2)                                     FL16340
C
      COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),VX2(47),AWCCON(5)          FL16350
      COMMON /DESAER/ EXT(47,4),ABS(47,4),G(47,4)                        FL16360
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL16370
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL16380
     *                     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4       FL16390
      DIMENSION DESEXT(47),DESSCA(47),DESABS(47),DESG(47),WIND(4)        FL16400
      REAL      DESEXT    ,DESSCA    ,DESABS    ,DESG    ,WIND           FL16410
      INTEGER WAVEL                                                      FL16420
C
      DATA I_3/3/
C
      DATA WIND/0., 10., 20., 30./                                       FL16430
      DATA RAYSCT / 0.01159 /                                            FL16440
      IF (WSPD.LT.0.) WSPD = 10.                                         FL16450
C                                                                        FL16460
      NWSPD = INT(WSPD/10)+1                                             FL16470
      IF (NWSPD.GE.5) WRITE (IPR,905)                                    FL16480
      NWSPD = MIN(NWSPD,I_3)                                               FL16490
C                                                                        FL16500
C     INTERPOLATE THE RADIATIVE PROPERTIES AT WIND SPEED WSPD            FL16510
C                                                                        FL16520
      DO 10 WAVEL = 1, 47                                                FL16530
C                                                                        FL16540
C        EXTINCTION COEFFICIENT                                          FL16550
C                                                                        FL16560
         SLOPE = LOG(EXT(WAVEL,NWSPD+1)/EXT(WAVEL,NWSPD))/(WIND(NWSPD+1) FL16570
     *      -WIND(NWSPD))                                                FL16580
         B = LOG(EXT(WAVEL,NWSPD+1))-SLOPE*WIND(NWSPD+1)                 FL16590
         DESEXT(WAVEL) = EXP(SLOPE*WSPD+B)                               FL16600
C                                                                        FL16610
C        ABSORPTION COEFFICIENT                                          FL16620
C                                                                        FL16630
         SLOPE = LOG(ABS(WAVEL,NWSPD+1)/ABS(WAVEL,NWSPD))/(WIND(NWSPD+1) FL16640
     *      -WIND(NWSPD))                                                FL16650
         B = LOG(ABS(WAVEL,NWSPD+1))-SLOPE*WIND(NWSPD+1)                 FL16660
         DESABS(WAVEL) = EXP(SLOPE*WSPD+B)                               FL16670
C                                                                        FL16680
C        SCATTERING COEFFICIENT                                          FL16690
C                                                                        FL16700
         DESSCA(WAVEL) = DESEXT(WAVEL)-DESABS(WAVEL)                     FL16710
C                                                                        FL16720
C        ASYMMETRY PARAMETER                                             FL16730
C                                                                        FL16740
         SLOPE = (G(WAVEL,NWSPD+1)-G(WAVEL,NWSPD))/(WIND(NWSPD+1)-       FL16750
     *      WIND(NWSPD))                                                 FL16760
         B = G(WAVEL,NWSPD+1)-SLOPE*(WIND(NWSPD+1))                      FL16770
         DESG(WAVEL) = SLOPE*WSPD+B                                      FL16780
   10 CONTINUE                                                           FL16790
C                                                                        FL16800
      EXT55 = DESEXT(4)                                                  FL16810
C                                                                        FL16820
C     DETERMINE METEROLOGICAL RANGE FROM 0.55 EXTINCTION                 FL16830
C     AND KOSCHMIEDER FORMULA                                            FL16840
C                                                                        FL16850
      IF (VIS.LE.0.) THEN                                                FL16860
         VIS = 3.912/(DESEXT(4)+RAYSCT)                                  FL16870
      ENDIF                                                              FL16880
C                                                                        FL16890
C     RENORMALIZE ATTENUATION COEFFICIENTS TO 1.0 KM-1 AT                FL16900
C     0.55 MICRONS FOR CAPABILTY WITH LOWTRAN                            FL16910
C                                                                        FL16920
      DO 20 WAVEL = 1, 47                                                FL16930
         EXTC(1,WAVEL) = DESEXT(WAVEL)/EXT55                             FL16940
C                                                                        FL16950
C        C          DESSCA(WAVEL) = DESSCA(WAVEL)       /EXT55           FL16960
C                                                                        FL16970
         ABSC(1,WAVEL) = DESABS(WAVEL)/EXT55                             FL16980
         ASYM(1,WAVEL) = DESG(WAVEL)                                     FL16990
   20 CONTINUE                                                           FL17000
      WRITE (IPR,900) VIS,WSPD                                           FL17010
      RETURN                                                             FL17020
C                                                                        FL17030
  900  FORMAT(//,'  VIS = ',F10.3,' WIND = ',F10.3)                      FL17040
  905  FORMAT(' WARNING: WIND SPEED IS BEYOND 30 M/S; RADIATIVE',        FL17050
     *'PROPERTIES',/,'OF THE DESERT AEROSOL HAVE BEEN EXTRAPOLATED')     FL17060
C                                                                        FL17070
      END                                                                FL17080
      BLOCK DATA DSTDTA                                                  FL17090
C                                                                        FL17100
C     >    BLOCK DATA                                                    FL17110
C     ****************************************************************** FL17120
C     *                                                                  FL17130
C     *    DESERT AEROSOL EXTINCTION COEFFICIENTS, ABSORPTION COEFFICIEN FL17140
C     *    AND ASYMMETRY PARAMETERS FOR FOUR WIND SPEEDS: 0 M/S, 10 M/S, FL17150
C     *    20 M/S AND 30 M/S                                             FL17160
C     *                                                                  FL17170
C     *    PROGRAMMED BY:  D. R. LONGTIN         OPTIMETRICS, INC.       FL17180
C     *                                          BURLINGTON, MASSACHUSET FL17190
C     *                                          FEB  1988               FL17200
C     *                                                                  FL17210
C     ****************************************************************** FL17220
C                                                                        FL17230
      COMMON /DESAER/DESEX1(47),DESEX2(47),DESEX3(47),DESEX4(47),        FL17240
     *DESAB1(47),DESAB2(47),DESAB3(47),DESAB4(47),DESG1(47),DESG2(47),   FL17250
     *DESG3(47),DESG4(47)                                                FL17260
C                                                                        FL17270
C     EXTINCTION COEFFICIENTS                                            FL17280
C                                                                        FL17290
      DATA DESEX1 /                                                      FL17300
     * 8.7330E-2, 7.1336E-2, 6.5754E-2, 4.0080E-2, 2.8958E-2, 1.4537E-2, FL17310
     * 7.1554E-3, 4.3472E-3, 3.5465E-3, 2.9225E-3, 2.5676E-3, 4.3573E-3, FL17320
     * 5.7479E-3, 2.9073E-3, 2.0109E-3, 1.8890E-3, 1.8525E-3, 1.8915E-3, FL17330
     * 1.9503E-3, 2.3256E-3, 4.9536E-3, 2.0526E-3, 2.6738E-3, 9.2804E-3, FL17340
     * 1.5352E-2, 6.9396E-3, 2.2455E-3, 1.9840E-3, 1.9452E-3, 1.9019E-3, FL17350
     * 1.8551E-3, 1.9661E-3, 1.9865E-3, 2.4089E-3, 1.7485E-3, 1.4764E-3, FL17360
     * 2.2604E-3, 2.1536E-3, 2.3008E-3, 2.9272E-3, 2.6943E-3, 2.4319E-3, FL17370
     * 1.9199E-3, 1.4887E-3, 8.0630E-4, 4.6950E-4, 2.0792E-4/            FL17380
      DATA DESEX2 /                                                      FL17390
     * 1.0419E-1, 8.8261E-2, 8.2699E-2, 5.7144E-2, 4.6078E-2, 3.1831E-2, FL17400
     * 2.4638E-2, 2.1952E-2, 2.1254E-2, 2.0743E-2, 2.0397E-2, 2.2340E-2, FL17410
     * 2.3848E-2, 2.1104E-2, 2.0422E-2, 2.0462E-2, 2.0591E-2, 2.0843E-2, FL17420
     * 2.1030E-2, 2.1630E-2, 2.2880E-2, 1.9075E-2, 2.0928E-2, 2.9835E-2, FL17430
     * 3.8025E-2, 2.7349E-2, 2.1502E-2, 2.1475E-2, 2.1563E-2, 2.1726E-2, FL17440
     * 2.2265E-2, 2.2580E-2, 2.2708E-2, 2.1705E-2, 2.1230E-2, 2.0523E-2, FL17450
     * 2.6686E-2, 2.5461E-2, 2.3785E-2, 2.6033E-2, 2.6484E-2, 2.6464E-2, FL17460
     * 2.5318E-2, 2.3341E-2, 1.7824E-2, 1.3092E-2, 7.2020E-3/            FL17470
      DATA DESEX3 /                                                      FL17480
     * 2.7337E-1, 2.5795E-1, 2.5252E-1, 2.2773E-1, 2.1710E-1, 2.0402E-1, FL17490
     * 1.9809E-1, 1.9664E-1, 1.9635E-1, 1.9655E-1, 1.9661E-1, 1.9907E-1, FL17500
     * 2.0164E-1, 1.9957E-1, 2.0013E-1, 2.0142E-1, 2.0270E-1, 2.0400E-1, FL17510
     * 2.0501E-1, 2.0665E-1, 2.0573E-1, 1.9165E-1, 2.0121E-1, 2.2402E-1, FL17520
     * 2.4718E-1, 2.2503E-1, 2.0749E-1, 2.0910E-1, 2.0999E-1, 2.1165E-1, FL17530
     * 2.1784E-1, 2.1727E-1, 2.1803E-1, 2.0995E-1, 2.1214E-1, 2.1308E-1, FL17540
     * 2.5226E-1, 2.4234E-1, 2.2638E-1, 2.3991E-1, 2.4680E-1, 2.5176E-1, FL17550
     * 2.5655E-1, 2.5505E-1, 2.3610E-1, 2.1047E-1, 1.5938E-1/            FL17560
      DATA DESEX4 /                                                      FL17570
     * 1.9841E0, 1.9721E0, 1.9676E0, 1.9488E0, 1.9424E0, 1.9377E0,       FL17580
     * 1.9374E0, 1.9484E0, 1.9509E0, 1.9549E0, 1.9570E0, 1.9642E0,       FL17590
     * 1.9737E0, 1.9764E0, 1.9860E0, 1.9944E0, 2.0020E0, 2.0113E0,       FL17600
     * 2.0148E0, 2.0245E0, 2.0283E0, 1.9397E0, 1.9973E0, 2.1039E0,       FL17610
     * 2.2246E0, 2.1587E0, 2.0409E0, 2.0520E0, 2.0613E0, 2.0651E0,       FL17620
     * 2.1194E0, 2.1065E0, 2.1104E0, 2.0651E0, 2.0926E0, 2.1155E0,       FL17630
     * 2.3696E0, 2.2931E0, 2.1828E0, 2.2708E0, 2.3304E0, 2.3762E0,       FL17640
     * 2.4533E0, 2.4915E0, 2.5118E0, 2.4463E0, 2.2122E0/                 FL17650
C                                                                        FL17660
C     ABSORPTION COEFFICIENTS                                            FL17670
C                                                                        FL17680
      DATA DESAB1 /                                                      FL17690
     * 6.4942E-4, 6.1415E-4, 5.8584E-4, 4.4211E-4, 1.3415E-4, 7.8142E-5, FL17700
     * 5.7566E-5, 8.3848E-5, 7.6988E-5, 4.4486E-5, 8.9604E-5, 2.4887E-3, FL17710
     * 3.3444E-3, 6.8781E-4, 1.6387E-4, 3.5236E-4, 3.5340E-4, 4.0930E-4, FL17720
     * 5.0526E-4, 8.2146E-4, 3.7647E-3, 1.0162E-3, 1.3525E-3, 7.7761E-3, FL17730
     * 1.3108E-2, 5.1252E-3, 1.0973E-3, 6.8573E-4, 5.7622E-4, 5.1268E-4, FL17740
     * 7.6834E-4, 5.3793E-4, 5.0611E-4, 1.2828E-3, 6.7827E-4, 4.3826E-4, FL17750
     * 5.1221E-4, 8.8642E-4, 9.5535E-4, 1.0000E-3, 7.5646E-4, 6.1552E-4, FL17760
     * 4.6087E-4, 3.5642E-4, 2.3556E-4, 1.7596E-4, 1.1699E-4/            FL17770
      DATA DESAB2 /                                                      FL17780
     * 4.3569E-3, 4.3413E-3, 4.3277E-3, 4.0649E-3, 3.9091E-4, 8.4594E-5, FL17790
     * 5.8501E-5, 8.4412E-5, 7.7547E-5, 4.6817E-5, 9.2721E-5, 2.5389E-3, FL17800
     * 3.3588E-3, 7.9414E-4, 8.5079E-4, 4.6002E-3, 4.4872E-3, 4.6200E-3, FL17810
     * 5.2973E-3, 4.8910E-3, 8.9899E-3, 5.4745E-3, 3.6375E-3, 1.1862E-2, FL17820
     * 1.5179E-2, 7.0015E-3, 8.4693E-3, 6.9516E-3, 6.3008E-3, 6.3684E-3, FL17830
     * 8.4992E-3, 6.9625E-3, 6.5192E-3, 7.8955E-3, 7.7192E-3, 5.8540E-3, FL17840
     * 5.3263E-3, 9.3004E-3, 7.4848E-3, 3.0952E-3, 1.8219E-3, 1.3078E-3, FL17850
     * 1.0653E-3, 5.5231E-4, 3.2311E-4, 2.2422E-4, 1.3839E-4/            FL17860
      DATA DESAB3 /                                                      FL17870
     * 4.1552E-2, 4.1671E-2, 4.1781E-2, 4.1125E-2, 5.0552E-3, 2.1085E-4, FL17880
     * 7.5703E-5, 9.5531E-5, 8.8354E-5, 9.0588E-5, 1.5058E-4, 3.4972E-3, FL17890
     * 3.6310E-3, 2.6709E-3, 1.2558E-2, 5.9184E-2, 5.8289E-2, 5.9206E-2, FL17900
     * 6.5487E-2, 5.8707E-2, 7.4669E-2, 5.2152E-2, 2.5783E-2, 4.7971E-2, FL17910
     * 3.2378E-2, 2.4739E-2, 8.1225E-2, 7.5085E-2, 7.1232E-2, 7.3042E-2, FL17920
     * 8.0638E-2, 7.8255E-2, 7.4882E-2, 7.8853E-2, 8.1412E-2, 6.5722E-2, FL17930
     * 4.8565E-2, 8.4983E-2, 7.1273E-2, 3.0870E-2, 1.7031E-2, 1.1455E-2, FL17940
     * 1.0554E-2, 4.0418E-3, 2.1509E-3, 1.4115E-3, 7.9698E-4/            FL17950
      DATA DESAB4 /                                                      FL17960
     * 4.1777E-1, 4.1880E-1, 4.2000E-1, 4.1846E-1, 8.6452E-2, 2.6538E-3, FL17970
     * 4.0804E-4, 3.1418E-4, 2.9996E-4, 9.3018E-4, 1.2814E-3, 2.1436E-2, FL17980
     * 8.7553E-3, 3.7670E-2, 2.0849E-1, 7.0914E-1, 7.0420E-1, 7.1379E-1, FL17990
     * 7.6309E-1, 7.1128E-1, 8.2992E-1, 5.3585E-1, 2.4456E-1, 3.8103E-1, FL18000
     * 1.7784E-1, 1.9305E-1, 7.9910E-1, 7.8987E-1, 7.7502E-1, 7.9400E-1, FL18010
     * 7.6332E-1, 8.3629E-1, 8.1581E-1, 8.3122E-1, 8.4901E-1, 7.0150E-1, FL18020
     * 4.4205E-1, 7.7354E-1, 7.1088E-1, 3.9328E-1, 2.3337E-1, 1.6258E-1, FL18030
     * 1.5289E-1, 5.8849E-2, 3.5576E-2, 2.4463E-2, 1.4525E-2/            FL18040
C                                                                        FL18050
C     ASYMMETRY PARAMETER                                                FL18060
C                                                                        FL18070
      DATA DESG1 /                                                       FL18080
     * 0.6603, 0.6581, 0.6547, 0.6383, 0.6276, 0.5997, 0.5829, 0.5873,   FL18090
     * 0.5967, 0.6130, 0.6323, 0.6850, 0.6068, 0.6312, 0.6816, 0.7298,   FL18100
     * 0.7574, 0.7874, 0.8124, 0.8424, 0.8301, 0.8107, 0.6143, 0.6167,   FL18110
     * 0.4892, 0.4917, 0.6662, 0.6334, 0.6298, 0.6498, 0.7470, 0.6711,   FL18120
     * 0.6751, 0.7538, 0.8054, 0.7797, 0.5522, 0.6575, 0.4702, 0.3719,   FL18130
     * 0.3626, 0.3690, 0.3790, 0.3805, 0.3766, 0.3639, 0.3281/           FL18140
      DATA DESG2 /                                                       FL18150
     * 0.6836, 0.6879, 0.6877, 0.6919, 0.6901, 0.7045, 0.7279, 0.7466,   FL18160
     * 0.7522, 0.7568, 0.7629, 0.7700, 0.7567, 0.7617, 0.7781, 0.8289,   FL18170
     * 0.8360, 0.8465, 0.8624, 0.8707, 0.9524, 0.8292, 0.6202, 0.6425,   FL18180
     * 0.5777, 0.5623, 0.7610, 0.7310, 0.7247, 0.7419, 0.7782, 0.7481,   FL18190
     * 0.7446, 0.8090, 0.8415, 0.8110, 0.6120, 0.7106, 0.5739, 0.4421,   FL18200
     * 0.4089, 0.3979, 0.3917, 0.3853, 0.3842, 0.3829, 0.3797/           FL18210
      DATA DESG3 /                                                       FL18220
     * 0.7718, 0.7865, 0.7907, 0.8077, 0.7801, 0.7827, 0.7871, 0.7880,   FL18230
     * 0.7887, 0.7888, 0.7894, 0.7909, 0.7882, 0.7934, 0.8103, 0.8729,   FL18240
     * 0.8766, 0.8844, 0.8979, 0.8997, 0.9698, 0.8318, 0.6197, 0.6420,   FL18250
     * 0.5797, 0.5698, 0.8014, 0.7938, 0.7901, 0.8069, 0.7894, 0.8139,   FL18260
     * 0.8086, 0.8546, 0.8691, 0.8288, 0.6394, 0.7400, 0.6495, 0.5235,   FL18270
     * 0.4793, 0.4583, 0.4376, 0.4169, 0.4006, 0.3941, 0.3875/           FL18280
      DATA DESG4 /                                                       FL18290
     * 0.8290, 0.8407, 0.8443, 0.8500, 0.8087, 0.7994, 0.7988, 0.7987,   FL18300
     * 0.7988, 0.7989, 0.7998, 0.8023, 0.8011, 0.8076, 0.8331, 0.9045,   FL18310
     * 0.9083, 0.9149, 0.9266, 0.9263, 0.9783, 0.8321, 0.6168, 0.6379,   FL18320
     * 0.5706, 0.5673, 0.8196, 0.8324, 0.8347, 0.8549, 0.7940, 0.8621,   FL18330
     * 0.8588, 0.8918, 0.8922, 0.8407, 0.6488, 0.7557, 0.7021, 0.6024,   FL18340
     * 0.5533, 0.5280, 0.5016, 0.4711, 0.4396, 0.4230, 0.4058/           FL18350
      END                                                                FL18360
C
C     *****************************************************************
C
      SUBROUTINE FLAYZ(ML,MODEL,ICLD,IAERSL,ZMDL,ALTZ,n_lvl,
     *     GNDALT,IVSA,IEMSCT) 
C                                                                        FL18380
C     SUBROUTINE TO CREATE FINAL LOWTRAN BOUNDRIES                       FL18390
C                                                                        FL18400
C     ZMDL COMMON /MODEL/ FINAL ALTITUDE FOR LOWTRAN                     FL18410
C     ZCLD CLOUD ALTITUDE                                                FL18420
C     ZK1 USED WITH VSA                                                  FL18430
C     ALTZ BLANK COMMON LBLRTM ALTITUDES                                   FL18440
C     ZNEW ALTITUDES ABOVE THE CLOUD                                     FL18450
C     ZNEWV ALTITUDES ABOVE THE 1ST 9 VSA ALTITUDES                      FL18460
C     ZTST  =ZCLD(J)                                                     FL18470
C     ZVSA  VSA ALTITUDES                                                FL18480
C                                                                        FL18490
      COMMON /LCRD2A/ CTHIK,CALT,CEXT                                    FL18500
      COMMON /ZVSALY/ ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)             FL18510
      DIMENSION ZNEWV(24),ALTZ(0:*),ZMDL( *)                                FL18520
      DIMENSION ZNEW(17),ZCLD(16),ZAER(34),ZST(234)                       FL18530
      DATA ZCLD/ 0.0,0.16,0.33,0.66,1.0,1.5,2.0,2.4,2.7,                 FL18540
     * 3.0,3.5,4.0,4.5,5.0,5.5,6.0/                                      FL18550
      DATA ZNEWV/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,                 FL18560
     * 14.,16.,18.,20.,22.,25.,30.,35.,40.,50.,70.,100./                 FL18570
      DATA ZNEW/ 7.,8.,9.,10.,12.,14.,16.,18.,20.,22.,25.,30.,           FL18580
     * 35.,40.,50.,70.,100./                                             FL18590
      DATA ZAER / 0., 1., 2., 3., 4., 5., 6., 7., 8., 9.,                FL18600
     *           10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,                FL18610
     *           20.,21.,22.,23.,24.,25.,30.,35.,40.,45.,                FL18620
     *           50.,70.,100.,   1000./                                  FL18630
      DATA DELZ /0.02/                                                   FL18640

      parameter (maxlay=205)

      IF (IAERSL.EQ.7) GO TO 300                                         FL18650
C                                                                        FL18660
      IF (MODEL.EQ.0) GO TO 140                                          FL18670
c
      IF (IVSA.EQ.1) THEN                                                FL18680
         DO 10 I = 1, 9                                                  FL18690
            ZMDL(I) = ZVSA(I)                                            FL18700
   10    CONTINUE                                                        FL18710
C                                                                        FL18720
         HMXVSA = ZVSA(9)                                                FL18730
         ZK1 = HMXVSA+0.01                                               FL18740
         IF (HMXVSA.LT.2.) ML = 33                                       FL18750
         IF (HMXVSA.LT.1.) ML = 34                                       FL18760
         IF (HMXVSA.EQ.2.) ML = 32                                       FL18770
         MDEL = 34-ML                                                    FL18780
         DO 20 K = 1, ML                                                 FL18790
            IK = K-10+MDEL                                               FL18800
            IF (IK.GE.1) ZMDL(K) = ZNEWV(IK)                             FL18810
            IF (K.EQ.10) ZMDL(K) = ZK1                                   FL18820
   20    CONTINUE                                                        FL18830
C                                                                        FL18840
         RETURN                                                          FL18850
      ELSE                                                               FL18860
         ML = 46                                                         FL18870
      ENDIF                                                              FL18880
C                                                                        FL18890
      IF (ICLD.GE.1.AND.ICLD.LE.11) GO TO 110                            FL18900
      DO 30 I = 1, n_lvl                                                 FL18910
         IF (ALTZ(I-1).GT.100.) GO TO 40                                     FL18920
         IL = I                                                          FL18930
         ZMDL(I) = ALTZ(I-1)                                                 FL18940
   30 CONTINUE                                                           FL18950
   40 ML = IL                                                            FL18960
C                                                                        FL18970
C     IF(IEMSCT.NE.0) ZMDL(ML)=100.                                      FL18980
C                                                                        FL18990
      IF (GNDALT.LE.0.) GO TO 60                                         FL19000
      DALT = (6.-GNDALT)/6.                                              FL19010
      IF (DALT.LE.0.) GO TO 60                                           FL19020
C                                                                        FL19030
      DO 50 I = 1, 6                                                     FL19040
         ZMDL(I) =  REAL(I-1)*DALT+GNDALT                                FL19050
   50 CONTINUE                                                           FL19060
   60 IF (ICLD.EQ.18.OR.ICLD.EQ.19) THEN                                 FL19070
c******%%%%%%%%
c         CLDD = 0.1*CTHIK                                                FL19080
c         CLD0 = CALT-0.5*CLDD                                            FL19090
c         CLD0 = MAX(CLD0,0.)                                             FL19100
c         CLD1 = CLD0+CLDD                                                FL19110
c         CLD2 = CLD1+CTHIK-CLDD                                          FL19120
c         CLD3 = CLD2+CLDD                                                FL19130

         cld1 = calt
         cld2 =calt + cthik

         cld0 = cld1 - 0.010
         IF (CLD0.LE.0.) CLD0 = 0.                                       FL08590
         cld3 = cld2 + 0.010
c
c     find mdl lvl above cld bottom
         DO 70 I = 1, ML                                                 FL19140
            IJ = I                                                       FL19150
            IF (ZMDL(I).LT.CLD1) GO TO 70                                FL19160
            GO TO 80                                                     FL19170
   70    CONTINUE                                                        FL19180

c     abort if 
         GO TO 300                                                       FL19190
   80    mcld = ij
c     save model levels
         DO 90 I = mcld, ML
            ZST(I) = ZMDL(I)
 90      CONTINUE                                                        FL19230
c     insert cloud- make small layers at cloud bottom and top
c     to trick path intergration in geo
         ZMDL(IJ) = CLD0                                                 FL19240
         ZMDL(IJ+1) = CLD1                                               FL19250
         ZMDL(IJ+2) = CLD2                                               FL19260
         ZMDL(IJ+3) = CLD3                                               FL19270
         IJ = IJ + 3                      
c     restore rest of zmdl above cloud top
         DO 100 I = mcld,ml
            IF (ZST(I).LT.CLD3) GO TO 100                                FL19300
            IJ = IJ+1                                                    FL19310
            ZMDL(IJ) = ZST(I)         
  100    CONTINUE                                                        FL19340
c
      ml = ij
c
      ENDIF                                                              FL19350
c
      GO TO 300                                                          FL19360
c____________________________________________________________
C                                                                        FL19370
C     STAND CLOUD                                                        FL19380
C                                                                        FL19390
  110 continue

c     icld ge 1 and le 11 to reach this point

c     cld model has 16 lvls; upper mdl atmosphere has 17 lvls

      DO 120 I = 1, 16                                                   FL19400
         ZMDL(I) = ZCLD(I)+GNDALT                                        FL19410
  120 CONTINUE                                                           FL19420
C                                                                        FL19440
      ml = 16 + 17

      DO 130 i = 17, ML                                                  FL19450
         ZMDL(I) = ZNEW(i-16)                                               FL19490
  130 CONTINUE                                                           FL19500
c
      GO TO 300                                                          FL19520
c____________________________________________________________
C                                                                        FL19530
C     MODEL 7                                                            FL19540
C                                                                        FL19550
  140 CONTINUE                                                           FL19560
c
      ml = n_lvl
c
      IF (ICLD.EQ.0) GO TO 280                                           FL19570
      IF (ICLD.EQ.20) GO TO 280                                          FL19580
      IF (IVSA.EQ.1) GO TO 280                                           FL19590
      IF (ML.EQ.1) GO TO 280                                             FL19600
      KK = 0                                                             FL19610

      DO 150 I = 0, n_lvl-1                                                  FL19620
         IF (ALTZ(I).GT.6.0) GO TO 160                                     FL19630
         KK = I                                                          FL19640
  150 CONTINUE                                                           FL19650
  160 IF (KK.LT.1) GO TO 200                                             FL19660
C                                                                        FL19670
      I = 1                                                              FL19680
      J = 1                                                              FL19690
      K = 0
  170 ZTST = ZCLD(J)                                                     FL19710
      IF (ZCLD(J).LT.ALTZ(0)) THEN
         J = J+1                                                         FL19730
         GO TO 170                                                       FL19740
      ENDIF                                                              FL19750
      IF (ABS(ZTST-ALTZ(K)).LT.DELZ) GO TO 180                             FL19760
      IF (ZTST.LT.ALTZ(K)) THEN                                            FL19770
         ZMDL(I) = ZTST                                                  FL19780
         I = I+1                                                         FL19790
         J = J+1                                                         FL19800
      ELSE                                                               FL19810
         ZMDL(I) = ALTZ(K)                                                 FL19820
         I = I+1                                                         FL19830
         K = K+1                                                         FL19840
      ENDIF                                                              FL19850
      GO TO 190                                                          FL19860
C                                                                        FL19870
  180 ZMDL(I) = ALTZ(K)                                                    FL19880
      I = I+1                                                            FL19890
      J = J+1                                                            FL19900
      K = K+1                                                            FL19910
C                                                                        FL19920
  190 IF (K.GE.KK) GO TO 200                                             FL19930
      IF (J.GE.17) GO TO 200                                             FL19940
      IF (I.GT.maxlay) GO TO 220                                             FL19950
      GO TO 170                                                          FL19960
C                                                                        FL19970
  200 IF (KK.EQ.0) THEN                                                  FL19980
         I = 1                                                           FL19990
         KK = 1                                                          FL20000
      ENDIF                                                              FL20010
C                                                                        FL20020
      DO 210 M = KK, n_lvl-1                                                  FL20030
         ZMDL(I) = ALTZ(M)                                                 FL20040
         I = I+1                                                         FL20050
         IF (I.GT.maxlay) GO TO 220                                          FL20060
  210 CONTINUE                                                           FL20070
C                                                                        FL20080
  220 continue
      ML = I-1                                                           FL20090

c insert levels just above and below cloud bottom and top

      IF (ICLD.EQ.18.OR.ICLD.EQ.19) THEN                                    FL19070
c******%%%%%%%%
c         CLDD = 0.1*CTHIK                                                FL19080
c         CLD0 = CALT-0.5*CLDD                                            FL19090
c         CLD0 = MAX(CLD0,0.)                                             FL19100
c         CLD1 = CLD0+CLDD                                                FL19110
c         CLD2 = CLD1+CTHIK-CLDD                                          FL19120
c         CLD3 = CLD2+CLDD                                                FL19130

         cld1 = calt
         cld2 =calt + cthik

         cld0 = cld1 - 0.010
         IF (CLD0.LE.0.) CLD0 = 0.                                       FL08590
         cld3 = cld2 + 0.010
c
         DO 225 I = 1, ML                                                 FL19140
            IJ = I                                                       FL19150
            IF (ZMDL(I).LT.CLD1) GO TO 225                                FL19160
            GO TO 230                                                     FL19170
 225     CONTINUE                                                        FL19180
         GO TO 300                                                       FL19190
 230     mcld = ij
c     save model levels
         DO 235 I = mcld, ML
            ZST(I) = ZMDL(I)
 235     CONTINUE                                                        FL19230
c     insert cloud- make small layers at cloud bottom and top
c     to trick path intergration in geo
         ZMDL(IJ) = CLD0                                                 FL19240
         ZMDL(IJ+1) = CLD1                                               FL19250
         ZMDL(IJ+2) = CLD2                                               FL19260
         ZMDL(IJ+3) = CLD3                                               FL19270
         IJ = IJ + 3                      
c     restore rest of zmdl above cloud top
         DO 240 I = mcld,ml
            IF (ZST(I).LT.CLD3) GO TO 240                                FL19300
            IJ = IJ+1                                                    FL19310
            ZMDL(IJ) = ZST(I)         
 240     CONTINUE                                                        FL19340
c
         ml = ij
c     
      ENDIF                                                              FL19350
c
      GO TO 300                                                          FL20100
c____________________________________________________________________
C                                                                        FL20110
 280  CONTINUE
      DO 285 I = 1, ML                                                   FL20120
         ZMDL(I) = ALTZ(I)                                                 FL20130
 285  CONTINUE                                                           FL20140
c____________________________________________________________________
C                                                                        FL20110
 300  CONTINUE


      RETURN  
      END                                                                FL20160
C
C     ******************************************************************
C
      SUBROUTINE TRANS                                                   FL20170
C                                                                        FL20180
C     ****************************************************************** FL20190
C     CALCULATES TRANSMITTANCE VALUES BETWEEN V1 AND V2                  FL20200
C     FOR A GIVEN ATMOSPHERIC SLANT PATH                                 FL20210
C                                                                        FL20220
C     MODIFIED FOR ASYMMETRY CALCULATION - JAN 1986 (A.E.R. INC.)        FL20230
C                                                                        FL20240
C     ****************************************************************** FL20250
C                                                                        FL20260
C     K           WPATH(IK,K)                                            FL20270
C                                                                        FL20280
C     6    MOLECULAR (RAYLIEGH) SCATTERING                               FL20290
C     7    BOUNDRY LAYER AEROSOL (0 TO 2 KM)                             FL20300
C     12    TROPOSPHERIC AEROSOL (2-10 KM)                               FL20310
C     13    STRATOSPHERIC  AEROSOL (10-30)                               FL20320
C     14    UPPER STRATOPHERIC (ABOVE 30KM)                              FL20330
C     15    AEROSOL WEIGHTED RELATIVE HUMITY (0 TO 2 KM)                 FL20340
C     16    CIRRUS CLOUDS                                                FL20350
C     ****************************************************************** FL20360
C                                                                        FL20370
      PARAMETER (MAXDV=2050)                                             FL20380
      INTEGER PHASE,DIST                                                 FL20390
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                    FL20400
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)     FL20410
C
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),VH(16),TX(16),W(16)       FL20420
      COMMON WPATH(IM2,16),TBBY(IM2)                                     FL20430
      COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),VX2(47),AWCCON(5)          FL20470
C
      CHARACTER*8      HMOD                                              FL20440
C
      COMMON /CMN/ HMOD(3),ZM(MXZMD),PF(MXZMD),TF(MXZMD),RFNDXM(MXZMD),
     *          ZP(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),
     *          TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),
     *          AMTP(MXMOL,MXPDIM)                                        
C
      COMMON /PATHD/ PBAR(MXLAY),TBAR(MXLAY),AMOUNT(MXMOL,MXLAY),        FA00530
     *               WN2L(MXLAY),DVL(MXLAY),WTOTL(MXLAY),ALBL(MXLAY),    FA00540
     *               ADBL(MXLAY),AVBL(MXLAY),H2OSL(MXLAY),IPATH(MXLAY),  FA00550
     *               ITYL(MXLAY),SECNTA(MXLAY),HT1,HT2,ALTZ(0:MXLAY),    FA00560
     *               PZ(0:MXLAY),TZ(0:MXLAY)                             FA00570
C
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL20480
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL20490
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL20500
c
      COMMON /MODEL/ ZMDL(MXZMD),PMM(MXZMD),TMM(MXZMD),                  FL07730
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL07740
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)

c
      COMMON /RAIN/ RNPATH(IM2),RRAMTK(IM2)                              FL20510
      COMMON /LCRD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB   FL20520
      COMMON /LCRD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,    FL20530
     *     RAINRT                                                        FL20540
      COMMON /LCRD3/ H1,H2,ANGLE,RANGE,BETA,RE,LEN                       FL20550
      COMMON /LCRD4/ V1,V2,DV                                            FL20560
      COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                      FL20570
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,N_LVL,JH1              FL20580
      COMMON /AER/ XX1,XX2,XX3,XX4,XX5,                                  FL20590
     *     YY1,YY2,YY3,YY4,YY5,ZZ1,ZZ2,ZZ3,ZZ4,ZZ5                       FL20600
      CHARACTER*8      XID,       HMOLID,      YID                     
      Real*8               SECANT,       XALTZ
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),      FL20620
     *     WK(60),PZL,PZU,TZL,TZU,WN2   ,DVP,V1P,V2P,TBOUNF,EMISIV,      FL20630
     *     FSCDID(17),NMOL,LAYER,YI1,YID(10) ,LSTWDF                     FL20640
      REAL*8           VI1,VI2,V1P,V2P,VV,v_mid                                FL20650
      COMMON /LPANEL/ VI1,VI2,DVV,NLIMAP                                 FL20660
      COMMON /ZOUTP/ ZOUT(MXLAY),SOUT(MXLAY),RHOSUM(MXLAY),              FL20670
     *     AMTTOT(MXMOL),AMTCUM(MXMOL),ISKIP(MXMOL)                      FL20680
      EQUIVALENCE (VI1,PNLHDR(1))                                        FL20690
      EQUIVALENCE (FSCDID(17),NLIM)                                      FL20700
      EQUIVALENCE (XID(1),XFILHD(1))                                     FL20710
      DIMENSION XFILHD(2),PNLHDR(2),SRAI(MAXDV)                          FL20720
      DIMENSION ABST(MAXDV),SCTT(MAXDV),ASYT(MAXDV),ASYDM(MAXDV)         FL20730
      DIMENSION ABST_st(MAXDV),SCTT_st(MAXDV),ASYT_st(MAXDV),
     *     ASYDM_st(MAXDV)
      DIMENSION VID(6),VL10(5),SUMEXT(MAXDV)                             FL20740
      DATA VID/0.1,0.2,0.5,1.0,2.0,5.0/                                  FL20750
      DATA VL10/25.,50.,125.,250.,500./                                  FL20760
      DVP = DV                                                           FL20770
      IENT = 0                                                           FL20780
      SUMA = 0.                                                          FL20790
      FACTOR = 0.5                                                       FL20800
C                                                                        FL20810
C     CC                                                                 FL20820
C     CC    FREQUENCY CAN GO BELOW 350 CM-1 FOR LBLRTM                   FL20830
C     CC                                                                 FL20840
C                                                                        FL20850
      V2 = MIN(V2,50000.)                                                FL20860
      DV = MAX(DV,5.)                                                    FL20870
      ICOUNT = 0                                                         FL20880
      IEMISS = 0                                                         FL20890
      IF (IEMSCT.EQ.1.OR.IEMSCT.EQ.2) IEMISS = 1                         FL20900
      TCRRIS = EXP(-W(16)*2.)                                            FL20910
C                                                                        FL20920
C     234  FORMAT(2F8.4)                                                 FL20930
C     CC                                                                 FL20940
C     CC   SET LOWTRAN DV DEPENDING ON FREQUENCY RANGE                   FL20950
C     CC   CAN BE 0.1,0.2,0.5,1.0,2.0 OR 5.0                             FL20960
C     CC                                                                 FL20970
C                                                                        FL20980
      IF (V2.GT.300) THEN                                                FL20990
         V1 =  REAL(INT(V1/5.0+0.1))*5.0                                 FL21000
         V2 =  REAL(INT(V2/5.0+0.1))*5.0                                 FL21010
      ENDIF                                                              FL21020
      V2 = MAX(V2,V1)                                                    FL21030
      VDEL = V2-V1                                                       FL21040
      IF (V1.GT.350.) VIDV = 5.0                                         FL21050
      IF (V1.GT.350.) GO TO 70                                           FL21060
      IF (V1.LE.10.) THEN                                                FL21070
         DO 10 I = 1, 4                                                  FL21080
            IF (VDEL.LE.VL10(I)) GO TO 20                                FL21090
   10    CONTINUE                                                        FL21100
         IC = 5                                                          FL21110
         GO TO 30                                                        FL21120
   20    IC = I                                                          FL21130
   30    CONTINUE                                                        FL21140
         VIDV = VID(IC)                                                  FL21150
      ELSE                                                               FL21160
         DO 40 I = 4, 5                                                  FL21170
            IF (VDEL.LE.VL10(I)) GO TO 50                                FL21180
   40    CONTINUE                                                        FL21190
         IC = 6                                                          FL21200
         GO TO 60                                                        FL21210
   50    IC = I                                                          FL21220
   60    CONTINUE                                                        FL21230
         VIDV = VID(IC)                                                  FL21240
      ENDIF                                                              FL21250
   70 CONTINUE                                                           FL21260
      NLIM = (VDEL/VIDV)+5.                                              FL21270
      DVP = VIDV                                                         FL21280
      V1 = V1-2*DVP                                                      FL21290
      V2 = V2+2*DVP                                                      FL21300
      V1P = V1                                                           FL21310
      V2P = V2                                                           FL21320
      WRITE (IPR,900) V1,V2,DVP                                          FL21330
      RMAXDV =  REAL(MAXDV)                                              FL21340
      IF ((V2-V1)/DV.GT.RMAXDV) STOP 'TRANS; (V2-V1)/DV GT MAXDV '       FL21350
      DO 80 I = 1, MAXDV                                                 FL21360
         ABST(I) = 0.                                                    FL21370
         SCTT(I) = 0.                                                    FL21380
         ASYT(I) = 0.                                                    FL21390
         SRAI(I) = 0.                                                    FL21400
         ASYDM(I) = 0.                                                   FL21410
         ABST_st(I) = 0.                                                    FL21370
         SCTT_st(I) = 0.                                                    FL21380
         ASYT_st(I) = 0.                                                    FL21390
         ASYDM_st(I) = 0.                                                   FL21410
   80 CONTINUE                                                           FL21430
      REWIND IEXFIL                                                      FL21440
      CALL BUFOUT (IEXFIL,XFILHD(1),NFHDRF)                              FL21450
      IF (ICLD.EQ.20.AND.V1.LT.350.) WRITE (IPR,905)                     FL21460
      NLIMAP = NLIM                                                      FL21470
      XKT0 = 0.6951*296.                                                 FL21480
      BETA0 = 1./XKT0                                                    FL21490
C                                                                        FL21500
C     **   BEGINING OF   LAYER   LOOP                                    FL21510
C                                                                        FL21520
      VI1 = V1                                                           FL21530
      VI2 = V2                                                           FL21540
      mid_v = nlim/2
      v_mid = v1 + dvp*real(mid_v-1)

c     initialize flag for buffering out and output layer count

      i_bufout = 1
      laycnt   = 1
C                                                                        FL21560
      do nv = 1,nlim
         sumext(nv) = 0.
      enddo

      DO 130 IK = IKLO, IKMAX                                            FL21570
         W7 = WPATH(IK,7)                                                FL21580
         W12 = WPATH(IK,12)                                              FL21590
         W15 = WPATH(IK,15)                                              FL21600
         IF (W7.GT.0.0.AND.ICH(1).LE.7) W15 = W15/W7                     FL21610
         IF (W12.GT.0.0.AND.ICH(1).GT.7) W15 = W15/W12                   FL21620
C                                                                        FL21630
C        INVERSE OF LOG REL HUM                                          FL21640
C                                                                        FL21650
         W(15) = 100.-EXP(W15)                                           FL21660
         IF (W7.LE.0.0.AND.ICH(1).LE.7) W(15) = 0.                       FL21670
         IF (W12.LE.0.0.AND.ICH(1).GT.7) W(15) = 0.                      FL21680
C                                                                        FL21690
C        **   LOAD AEROSOL EXTINCTION AND ABSORPTION COEFFICIENTS        FL21700
C                                                                        FL21710
C        CC                                                              FL21720
C        CC    LOAD EXTINCTIONS AND ABSORPTIONS FOR 0.2-200.0 UM (1-46)  FL21730
C        CC                                                              FL21740
C                                                                        FL21750
         CALL EXABIN                                                     FL21760
C                                                                        FL21770
C        CC                                                              FL21780
C                                                                        FL21790
         VI = V1                                                         FL21800
         VI = VI-VIDV                                                    FL21810
         NV = 0                                                          FL21820
         XKT = 0.6951*TBBY(IK)                                           FL21830
         BETAR = 1./XKT                                                  FL21840

         rad_mid = RADFN(v_mid,XKT)    

C                                                                        FL21850
C        CC                                                              FL21860
C        CC   BEGINNING OF FREQUENCY LOOP                                FL21870
C        CC                                                              FL21880
C                                                                        FL21890
   90    CONTINUE                                                        FL21900
C                                                                        FL21910
C        CC                                                              FL21920
C                                                                        FL21930
         CSSA = 1.                                                       FL21940
         ASYMR = 1.                                                      FL21950
         NV = NV+1                                                       FL21960
         VI = VI+VIDV                                                    FL21970
         V = ABS(VI)                                                     FL21980
         VV = V                                                          FL21990
C                                                                        FL22000
         SCTMOL = RAYSCT(V)*WPATH(IK,6)                                  FL22010
C                                                                        FL22020
         DVV = VIDV                                                      FL22030
C                                                                        FL22040
         RADFT = RADFN(VV,XKT)                                           FL22050
         RADFT0 = RADFN(VV,XKT0)                                         FL22060
C                                                                        FL22070
C        CC                                                              FL22080
C        CC    AEROSOL ATTENUATIONS                                      FL22090
C        CC                                                              FL22100
C                                                                        FL22110
         TRAIN = 0.0                                                     FL22120
C                                                                        FL22130
         CALL AEREXT (V,IK,RADFT)                                        FL22140
C                                                                        FL22150
         EXT = XX1*WPATH(IK,7)+XX2*WPATH(IK,12)+XX3*WPATH(IK,13)+XX4*    FL22160
     *      WPATH(IK,14)+XX5*WPATH(IK,16)                                FL22170
         ABT = YY1*WPATH(IK,7)+YY2*WPATH(IK,12)+YY3*WPATH(IK,13)+YY4*    FL22180
     *      WPATH(IK,14)+YY5*WPATH(IK,16)                                FL22190
C                                                                        FL22200
C        ASYMMETRY FACTOR IS WEIGHTED AVERAGE                            FL22210
C                                                                        FL22220
C        CC   ASY=(ZZ1*(XX1-YY1)*WPATH(IK,7)+ZZ2*(XX2-YY2)*WPATH(IK,12)+ FL22230
C        CC  + ZZ3*(XX3-YY3)*WPATH(IK,13)+ZZ4*(XX4-YY4)*WPATH(IK,14))/   FL22240
C        CC  + ((XX1-YY1)*WPATH(IK,7)+(XX2-YY2)*WPATH(IK,12)+            FL22250
C        CC  + (XX3-YY3)*WPATH(IK,13)+(XX4-YY4)*WPATH(IK,14)+SCTMOL)     FL22260
C                                                                        FL22270
         ASY = (ZZ1*(XX1-YY1)*WPATH(IK,7)+ZZ2*(XX2-YY2)*WPATH(IK,12)+ZZ3 FL22280
     *      *(XX3-YY3)*WPATH(IK,13)+ZZ4*(XX4-YY4)*WPATH(IK,14)+ZZ5*(XX5- FL22290
     *      YY5)*WPATH(IK,16))                                           FL22300
         SCT = EXT-ABT                                                   FL22310
         IF (VV.GE.350.AND.ICLD.EQ.20) ABT = ABT+(WPATH(IK,16)*2./RADFT) FL22320
C                                                                        FL22330
C        CC                                                              FL22340
C        CC   ADD CONTRIBUTION OF CLOUDS AND RAIN                        FL22350
C        CC                                                              FL22360
C                                                                        FL22370
         IF (RRAMTK(IK).NE.0.0) THEN                                     FL22380
            TRAIN = TNRAIN(RRAMTK(IK),VV,TBBY(IK),RADFT)                 FL22390
            IF (V.LT.250.) THEN                                          FL22400
               IF (ICLD.LE.11) PHASE = 1                                 FL22410
               IF (ICLD.GT.11) PHASE = 2                                 FL22420
               DIST = 1                                                  FL22430
C                                                                        FL22440
C              CALL SCATTERING ROUTINE TO OBTAIN ASYMMTRY FACTOR AND RAT FL22450
C              OF ABSORPTION TO EXTINCTION DUE TO RAIN WITHIN RANGE OF   FL22460
C              19 TO 231 GHZ                                             FL22470
C              EXTRAPOLATE ABOVE AND BELOW THAT FREQ RANGE               FL22480
C                                                                        FL22490
               CALL RNSCAT (V,RRAMTK(IK),TBBY(IK),PHASE,DIST,IK,CSSA,    FL22500
     *            ASYMR,IENT)                                            FL22510
               IENT = IENT+1                                             FL22520
            ELSE                                                         FL22530
               CSSA = 0.5                                                FL22540
               ASYMR = 0.85                                              FL22550
            ENDIF                                                        FL22560
         ENDIF                                                           FL22570
C                                                                        FL22580
C        SET EXT DUE TO RAIN FOR LAYER                                   FL22590
C                                                                        FL22600
         RNEXPT = TRAIN*RNPATH(IK)                                       FL22610
C                                                                        FL22620
C        PUT RADIATION  CLD BACK IN                                      FL22630
C                                                                        FL22640
         SRAI(NV) = SRAI(NV)+RNEXPT*RADFT                                FL22650
C                                                                        FL22660
         ABT = ABT+RNEXPT*CSSA                                           FL22670
         SCT = SCT+RNEXPT*(1.-CSSA)                                      FL22680
         ASY = ASY+ASYMR*(1.-CSSA)*RNEXPT                                FL22690
C                                                                        FL22700
C                                                                        FL22710
         SCT = SCT+SCTMOL                                                FL22720
C                                                                        FL22730
         EXT = SCT+ABT                                                   FL22740
C                                                                        FL22750
         IF (IK.LE.JH1) THEN                                             FL22760
C                                                                        FL22770
C           DOUBLE  TANGENT PATH LAYERS                                  FL22780
C                                                                        FL22790
            SUMEXT(NV) = SUMEXT(NV)+EXT*RADFT*2.0                        FL22800
            IF (IEMISS.EQ.0) THEN                                        FL22810
               EXT = EXT*2.                                              FL22820
               SCT = SCT*2.                                              FL22830
               ABT = ABT*2.                                              FL22840
            ENDIF                                                        FL22850
         ELSE                                                            FL22860
            SUMEXT(NV) = SUMEXT(NV)+EXT*RADFT                            FL22870
         ENDIF                                                           FL22880
C                                                                        FL22890
         IF (VV.GE.1.0) THEN                                             FL22900
            RADRAT = RADFT/RADFT0                                        FL22910
         ELSE                                                            FL22920
            RADRAT = BETAR/BETA0                                         FL22930
         ENDIF                                                           FL22940
C                                                                        FL22950
C        CC                                                              FL22960
C        CC    IF TRANSMISSION STORE THE ACCUMULATED AMOUNTS             FL22970
C        CC    IF EMISSION STORE THE AMOUNTS PER LAYER                   FL22980
C        CC                                                              FL22990
C                                                                        FL23000
         IF (IEMISS.EQ.1) THEN                                           FL23010
            ABST(NV) = ABST(NV)+ABT                                      FL23020
            SCTT(NV) = SCTT(NV)+SCT                                      FL23030
            ASYDM(NV) = ASYDM(NV)+SCT+SCTMOL                             FL23040
            ASYT(NV) = ASYT(NV)+ASY                                      FL23050
         ELSE                                                            FL23060
            ABST(NV) = ABST(NV)+ABT*RADRAT                               FL23070
            SCTT(NV) = SCTT(NV)+SCT*RADRAT                               FL23080
         ENDIF                                                           FL23090

C                                                                        FL23100
C        CC                                                              FL23110
C        CC    CIRRUS CLOUD SHOULD BE ADDED IN LATER                     FL23120
C        CC                                                              FL23130
C                                                                        FL23140

            IF (ASYDM(nv).GT.0.) THEN                                     FL23330
               ASYT(nv) = ASYT(nv)/ASYDM(nv)                              FL23340
            ELSE                                                         FL23350
               ASYT(nv) = 0.                                              FL23360
            ENDIF                                                        FL23370

         IF (VI.LT.V2) GO TO 90                                          FL23150
C                                                                        FL23160
C        CC                                                              FL23170
C        CC    ***END OF FREQUENCY LOOP                                  FL23180
C        CC                                                              FL23190
C        CC   BUFFER OUT ABSORPTION, SCATTERING, AND                     FL23200
C        CC   ASYMMETRY PANELS OF LAYERS BY FREQUENCY                    FL23210
C        CC   TO IEXFIL FOR USE IN LBLRTM                                FL23220
C        CC                                                              FL23230
C                                                                        FL23240
c  _____________________________________________

         IF (((IEMISS.GE.1).OR.(IK.EQ.IKMAX)) ) THEN      

            if (i_bufout .eq. -1) then
c     a thin layer has been identified
c     obtain absorption, scattering and asymmetry for combined layer
               do 110 i=1,nlim

                  abst_st(i) = abst_st(i)+abst(i)
                  sctt_sum_i = sctt_st(i)+sctt(i)

                  if (sctt_sum_i .gt. 0) then
                     f_asym = sctt(i)/sctt_sum_i
c    obtain a weighted average for the asymmetry paramter
c                     asyt_st(i) = (1-f_asym)*asyt_st(i)+f_asym*asyt(i)
                     asyt_st(i) = 
     *                       asyt_st(i) - f_asym*(asyt_st(i)-asyt(i))
                  else
                     asyt_st(i) = 0.
                  endif

                  sctt_st(i) = sctt_sum_i

 110           continue
            else
               do 115 i=1,nlim
                  abst_st(i) = abst(i)
                  sctt_st(i) = sctt(i)
                  asyt_st(i) = asyt(i)
 115           continue
            endif
c  _____________________________________________

            if (ikmax.ne.1 .and. 
     &          abs(zmdl(ik+1)-zout(laycnt+1)).gt.0.001) then
               i_bufout = -1
               go to 122
            else
c     
               CALL BUFOUT (IEXFIL,PNLHDR(1),NPHDRF)
c     
               CALL BUFOUT (IEXFIL,ABST_st(1),NLIM) 
               CALL BUFOUT (IEXFIL,SCTT_st(1),NLIM) 
               CALL BUFOUT (IEXFIL,ASYT_st(1),NLIM) 
c
               if (laycnt.eq.1) write(ipr,935) v_mid
c
               write(ipr,940) laycnt,zout(laycnt),zout(laycnt+1),
     *                 rad_mid*ABST_st(mid_v),rad_mid*sctt_st(mid_v),
     *                 asyt_st(mid_v)
c  
               DO 120 I = 1, nlim
                  ABST_st(I) = 0.                                          FL23460
                  SCTT_st(I) = 0.                                          FL23470
                  ASYT_st(I) = 0.                                          FL23480
                  ASYDM_st(I) = 0.                                         FL23490
 120           CONTINUE                                                    FL23500
               laycnt = laycnt +1
               i_bufout = 1
            endif
         ENDIF                                                           FL23520
c
 122     continue
c
         DO 125 I = 1, MAXDV                                             FL23450
            ABST(I) = 0.                                                 FL23460
            SCTT(I) = 0.                                                 FL23470
            ASYT(I) = 0.                                                 FL23480
            ASYDM(I) = 0.                                                FL23490
 125     CONTINUE                                                        FL23500
c               
C        ***END OF LAYER LOOP***    (IK LOOP)                            FL23540
C                                                                        FL23550
  130 CONTINUE                                                           FL23560
C                                                                        FL23570
C                                                                        FL23580
      REWIND IEXFIL                                                      FL23590
      VI = V1-VIDV                                                       FL23600
      DO 140 NV = 1, NLIM                                                FL23610
         VI = VI+VIDV                                                    FL23620
         IF (ICOUNT.EQ.0.OR.ICOUNT.EQ.50) THEN                           FL23630
            ICOUNT = 0                                                   FL23640
            IF (VI.GT.100.) WRITE (IPR,910)                              FL23650
            IF (VI.LE.100.) WRITE (IPR,915)                              FL23660
         ENDIF                                                           FL23670
         ICOUNT = ICOUNT+1                                               FL23680
         IF (SUMEXT(NV).LE.BIGEXP) THEN                                  FL23690
            TRAN = EXP(-SUMEXT(NV))                                      FL23700
         ELSE                                                            FL23710
            TRAN = 1.0/BIGNUM                                            FL23720
         ENDIF                                                           FL23730
         IF (SRAI(NV).LE.BIGEXP) THEN                                    FL23740
            TR1 = EXP(-SRAI(NV))                                         FL23750
         ELSE                                                            FL23760
            TR1 = 1.0/BIGNUM                                             FL23770
         ENDIF                                                           FL23780
         IF (VI.GT.V1) FACTOR = 1.0                                      FL23790
         IF (VI.GE.V2) FACTOR = 0.5                                      FL23800
         SUMA = SUMA+FACTOR*DVV*(1.0-TRAN)                               FL23810
         IF (VI.GT.100.) ALAM = 1.0E+04/VI                               FL23820
         IF (VI.LE.100.) ALAM = VI*29.979                                FL23830
         WRITE (IPR,920) VI,ALAM,TRAN,TR1                                FL23840
  140 CONTINUE                                                           FL23850
      IF (ICLD.EQ.20) WRITE (IPR,925) TCRRIS                             FL23860
      AB = 1.0-SUMA/(VI-V1)                                              FL23870
      WRITE (IPR,930) V1,VI,SUMA,AB                                      FL23880
c
      RETURN                                                             FL23890
C                                                                        FL23900
C     **   FORMAT STATEMENTS FOR SPECTRAL DATA                           FL23910
C     **   PAGE HEADERS                                                  FL23920
C                                                                        FL23930
  900 FORMAT(/,'1 LOWTRAN WAVENUMBER INTERVAL- V1P, V2P, DVP:',
     *                                                  3F10.3,//)
  905 FORMAT(' CIRRUS NOT DEFINED BELOW 350 CM-1')                       FL23950
  910 FORMAT ('1',/ 1X,'  FREQ WAVELENGTH  TOTAL    RAIN '/)       
  915 FORMAT ('1',/ 1X,'  FREQ FREQUENCY   TOTAL    RAIN  ',        
     *                 /2X,' CM-1    GHZ  ',2(4X,'TRANS')/)         
  920 FORMAT(1X,F7.1,F8.3,10F9.4,F12.3)                                  FL23990
  925 FORMAT('0TRANSMISSION DUE TO CIRRUS = ',F10.4)                     FL24000
  930 FORMAT('0INTEGRATED ABSORPTION FROM ',F9.3,' TO ',F9.3,' CM-1 =',  FL24010
     * F10.2,' CM-1',/,' AVERAGE TRANSMITTANCE =',F6.4,/)                FL24020
 935  FORMAT(//,'    LAYER OPTICAL PROPERTIES AT',f12.3, 
     * ' FOR THE PATH FROM Z(J) TO Z(J+1)',//,
     * T3,'J',T10,'Z(J)',T19,'Z(J+1)',
     * T27,'ABSORB',T37,'SCATTR',T47,'ASYM PAR',/,
     * T10,'(KM)',T20,'(KM)')
 940  format(i4,2f10.3,1p,3e10.2)
C                                                                        FL24030
      END                                                                FL24040
C
C     ******************************************************************
C
      SUBROUTINE RNSCAT(V,R,TT,PHASE,DIST,IK,CSSA,ASYMR,IENT)            FL24050
C                                                                        FL24060
C     ****************************************************************** FL24070
C                                                                        FL24080
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL24090
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL24100
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL24110
      INTEGER PHASE,DIST                                                 FL24120
      DIMENSION SC(3,4)                                                  FL24130
C                                                                        FL24140
C     ARGUMENTS:                                                         FL24150
C                                                                        FL24160
C     F = FREQUENCY (GHZ)                                                FL24170
C     R = RAINFALL RATE (MM/HR)                                          FL24180
C     T = TEMPERATURE (DEGREES CELSIUS)                                  FL24190
C     PHASE = PHASE PARAMETER (1=WATER, 2=ICE)                           FL24200
C     DIST = DROP SIZE DISTRIBUTION PARAMETER                            FL24210
C     (1=MARSHALL-PALMER, 2=BEST)                                        FL24220
C                                                                        FL24230
C     RESULTS:                                                           FL24240
C                                                                        FL24250
C     SC(1) = ABSORPTION COEFFICIENT (1/KM)                              FL24260
C     SC(2) = EXTINCTION COEFFICIENT (1/KM)                              FL24270
C     SC(I),I=3,NSC = LEGENDRE COEFFICIENTS #I-3  (NSC=10)               FL24280
C     2=BAD RAINFALL RATE, 3=BAD TEMPERATURE,                            FL24300
C     4=BAD PHASE PARAMETER, 5=BAD DROP SIZE DISTRIBUTION                FL24310
C                                                                        FL24320
C     THE INTERNAL DATA:                                                 FL24330
C                                                                        FL24340
      DIMENSION FR(9),TEMP(3)                                            FL24350
C                                                                        FL24360
C     FR(I),I=1,NF = TABULATED FREQUENCIES (GHZ)  (NF=9)                 FL24370
C     TEMP(I),I=1,NT = TABULATED TEMPERATURES  (NT=3)                    FL24380
C                                                                        FL24390
C     THE BLOCK-DATA SECTION                                             FL24400
C                                                                        FL24410
      DATA RMIN,RMAX/0.,50./,NF/9/,NT/3/,NSC/4/,MAXI/3/                  FL24420
      DATA TK/273.15/,CMT0/1.0/,C7500/0.5/,G0/0.0/,G7500/0.85/           FL24430
      DATA (TEMP(I),I=1,3)/-10.,0.,10./                                  FL24440
      DATA (FR(I),I=1,9)/19.35,37.,50.3,89.5,100.,118.,130.,183.,231./   FL24450
C                                                                        FL24460
C     THIS SUBROUTINE REQUIRES FREQUENCIES IN GHZ                        FL24470
C                                                                        FL24480
      NOPR = 0                                                           FL24490
      IF (IK.EQ.1) NOPR = 1                                              FL24500
      IF (IENT.GT.1) NOPR = 0                                            FL24510
      F = V*29.97925                                                     FL24520
      FSAV = F                                                           FL24530
      RSAV = R                                                           FL24540
      TSAV = T                                                           FL24550
      INT = 0                                                            FL24560
C                                                                        FL24570
C     CONVERT TEMP TO DEGREES CELSIUS                                    FL24580
C                                                                        FL24590
      T = TT-TK                                                          FL24600
C                                                                        FL24610
C     FREQ RANGE OF DATA 19.35-231 GHZ IF LESS THAN 19.35                FL24620
C     SET UP PARAMETERS FOR INTERPOLATION                                FL24630
C                                                                        FL24640
      IF (F.LT.FR(1)) THEN                                               FL24650
         FL = 0.0                                                        FL24660
         FM = FR(1)                                                      FL24670
         INT = 1                                                         FL24680
         IF (NOPR.GT.0) WRITE (IPR,900)                                  FL24690
      ENDIF                                                              FL24700
C                                                                        FL24710
C     IF MORE THAN 231 GHZ SET UP PARAMETERS FOR EXTRAPOLATION           FL24720
C                                                                        FL24730
      IF (F.GT.FR(NF)) THEN                                              FL24740
         FL = FR(NF)                                                     FL24750
         FM = 7500.                                                      FL24760
         INT = 2                                                         FL24770
         IF (NOPR.GT.0) WRITE (IPR,900)                                  FL24780
      ENDIF                                                              FL24790
C                                                                        FL24800
C     TEMP RANGE OF DATA IS -10 TO +10 DEGREES CELCIUS                   FL24810
C     IF BELOW OR ABOVE EXTREME SET AND DO CALCULATIONS AT EXTREME       FL24820
C                                                                        FL24830
      IF (T.LT.TEMP(1)) THEN                                             FL24840
         T = TEMP(1)                                                     FL24850
         IF (NOPR.GT.0) WRITE (IPR,905)                                  FL24860
      ENDIF                                                              FL24870
C                                                                        FL24880
      IF (T.GT.TEMP(3)) THEN                                             FL24890
         T = TEMP(3)                                                     FL24900
         IF (NOPR.GT.0) WRITE (IPR,905)                                  FL24910
      ENDIF                                                              FL24920
C                                                                        FL24930
C     RAIN RATE OF DATA IS FOR 0-50 MM/HR                                FL24940
C     IF GT 50 TREAT CALCULATIONS AS IF 50 MM/HR WAS INPUT               FL24950
C                                                                        FL24960
      IF (R.GT.50) THEN                                                  FL24970
         R = 50.                                                         FL24980
         IF (NOPR.GT.0) WRITE (IPR,910)                                  FL24990
      ENDIF                                                              FL25000
C                                                                        FL25010
      KI = 1                                                             FL25020
C                                                                        FL25030
C     FIGURE OUT THE SECOND INDEX                                        FL25040
C                                                                        FL25050
   10 J = PHASE+2*DIST                                                   FL25060
C                                                                        FL25070
C                                                                        FL25080
C     GET THE TEMPERATURE INTERPOLATION PARAMETER ST                     FL25090
C     IF NEEDED AND AMEND THE SECOND INDEX                               FL25100
C                                                                        FL25110
      CALL BS (J,T,TEMP,NT,ST)                                           FL25120
C                                                                        FL25130
C     FIGURE OUT THE THIRD INDEX AND THE FREQUENCY INTERPOLATION         FL25140
C     PARAMETER SF                                                       FL25150
C                                                                        FL25160
      CALL BS (K,F,FR,NF,SF)                                             FL25170
C                                                                        FL25180
C     INITIALIZE SC                                                      FL25190
C                                                                        FL25200
      DO 20 I = 1, NSC                                                   FL25210
         SC(KI,I) = 0.                                                   FL25220
   20 CONTINUE                                                           FL25230
      SC(KI,3) = 1.                                                      FL25240
C                                                                        FL25250
C     NOW DO THE CALCULATIONS                                            FL25260
C                                                                        FL25270
C     THE WATER CONTENT IS                                               FL25280
C                                                                        FL25290
      IF (DIST.EQ.1) THEN                                                FL25300
         WC = .0889*R**.84                                               FL25310
      ELSE                                                               FL25320
         WC = .067*R**.846                                               FL25330
      ENDIF                                                              FL25340
C                                                                        FL25350
C     FOR A TEMPERATURE DEPENDENT CASE, I.E.                             FL25360
C                                                                        FL25370
      IF (J.LT.3) THEN                                                   FL25380
         S1 = (1.-SF)*(1.-ST)                                            FL25390
         S2 = (1.-SF)*ST                                                 FL25400
         S3 = SF*(1.-ST)                                                 FL25410
         S4 = SF*ST                                                      FL25420
         DO 30 I = 1, MAXI                                               FL25430
            IF (I.LE.2) THEN                                             FL25440
               ISC = I                                                   FL25450
            ELSE                                                         FL25460
               ISC = I+1                                                 FL25470
            ENDIF                                                        FL25480
            SC(KI,ISC) = S1*TAB(I,J,K,WC)+S2*TAB(I,J+1,K,WC)+S3*TAB(I,J, FL25490
     *         K+1,WC)+S4*TAB(I,J+1,K+1,WC)                              FL25500
   30    CONTINUE                                                        FL25510
C                                                                        FL25520
C        FOR A TEMPERATURE INDEPENDENT CASE                              FL25530
C                                                                        FL25540
      ELSE                                                               FL25550
         S1 = 1.-SF                                                      FL25560
         S2 = SF                                                         FL25570
         DO 40 I = 1, MAXI                                               FL25580
            IF (I.LE.2) THEN                                             FL25590
               ISC = I                                                   FL25600
            ELSE                                                         FL25610
               ISC = I+1                                                 FL25620
            ENDIF                                                        FL25630
            SC(KI,ISC) = S1*TAB(I,J,K,WC)+S2*TAB(I,J,K+1,WC)             FL25640
   40    CONTINUE                                                        FL25650
      ENDIF                                                              FL25660
      F = FSAV                                                           FL25670
      IF (INT.EQ.3) GO TO 50                                             FL25680
      IF (INT.EQ.4) GO TO 60                                             FL25690
      IF (INT.EQ.0) THEN                                                 FL25700
         CSSA = SC(KI,1)/SC(KI,2)                                        FL25710
         CSSA = MIN(CSSA,1.0)                                            FL25720
         ASYMR = SC(KI,4)/3.0                                            FL25730
         F = FSAV                                                        FL25740
         R = RSAV                                                        FL25750
         T = TSAV                                                        FL25760
         RETURN                                                          FL25770
      ENDIF                                                              FL25780
      IF (INT.EQ.1) THEN                                                 FL25790
         INT = 3                                                         FL25800
         F = FM                                                          FL25810
         KI = 2                                                          FL25820
      ENDIF                                                              FL25830
      IF (INT.EQ.2) THEN                                                 FL25840
         INT = 4                                                         FL25850
         F = FL                                                          FL25860
         KI = 3                                                          FL25870
      ENDIF                                                              FL25880
      GO TO 10                                                           FL25890
   50 CONTINUE                                                           FL25900
      FDIF = FM-F                                                        FL25910
      FTOT = FM-FL                                                       FL25920
      CM = SC(KI,1)/SC(KI,2)                                             FL25930
      CM = MIN(CM,1.0)                                                   FL25940
      CL = CMT0                                                          FL25950
      AM = SC(KI,4)/3.0                                                  FL25960
      AL = G0                                                            FL25970
      GO TO 70                                                           FL25980
   60 CONTINUE                                                           FL25990
      FDIF = FM-F                                                        FL26000
      FTOT = FM-FL                                                       FL26010
      CM = C7500                                                         FL26020
      CL = SC(KI,1)/SC(KI,2)                                             FL26030
      CL = MIN(CL,1.0)                                                   FL26040
      AM = G7500                                                         FL26050
      AL = SC(KI,4)/3.0                                                  FL26060
   70 CTOT = CM-CL                                                       FL26070
      CAMT = FDIF*CTOT/FTOT                                              FL26080
      CSSA = CM-CAMT                                                     FL26090
      ATOT = AM-AL                                                       FL26100
      AAMT = FDIF*ATOT/FTOT                                              FL26110
      ASYMR = AM-AAMT                                                    FL26120
      F = FSAV                                                           FL26130
      R = RSAV                                                           FL26140
      T = TSAV                                                           FL26150
      RETURN                                                             FL26160
C                                                                        FL26170
  900 FORMAT(2X,'***  THE ASYMMETRY PARAMETER DUE TO RAIN IS BASED ON',  FL26180
     * 'DATA BETWEEN 19 AND 231 GHZ',                                    FL26190
     * /2X,'***  EXTRAPOLATION IS USED FOR FREQUENCIES LOWER AND',       FL26200
     * 'HIGHER THAN THIS RANGE')                                         FL26210
  905 FORMAT(2X,'***  TEMPERATURE RANGE OF DATA IS -10 TO +10 ',         FL26220
     *'DEGREES CELSIUS',/2X,'***  BEYOND THESE VALUES IT IS ',           FL26230
     *'TREATED AS IF AT THE EXTREMES')                                   FL26240
  910 FORMAT(2X,'***  RAIN RATES BETWEEN 0 AND 50 MM/HR ARE',            FL26250
     *'WITHIN THIS DATA RANGE',/2X,'***  ABOVE THAT THE ASYMMETRY',      FL26260
     *' PARAMETER IS CALCULATED FOR 50 MM/HR')                           FL26270
C                                                                        FL26280
      END                                                                FL26290
C
C     ******************************************************************
C
      SUBROUTINE BS(I,A,B,N,S)                                           FL26300
C                                                                        FL26310
C     ****************************************************************** FL26320
C                                                                        FL26330
      DIMENSION B(9)                                                     FL26340
C                                                                        FL26350
C     THIS SUBROUTINE DOES THE BINARY SEARCH FOR THE INDEX I             FL26360
C     SUCH THAT A IS IN BETWEEN B(I) AND B(I+1)                          FL26370
C     AND CALCULATES THE INTERPOLATION PARAMETER S                       FL26380
C     SUCH THAT A=S*B(I+1)+(1.-S)*B(I)                                   FL26390
C                                                                        FL26400
      I = 1                                                              FL26410
      J = N                                                              FL26420
   10 M = (I+J)/2                                                        FL26430
      IF (A.LE.B(M)) THEN                                                FL26440
         J = M                                                           FL26450
      ELSE                                                               FL26460
         I = M                                                           FL26470
      ENDIF                                                              FL26480
      IF (J.GT.I+1) GO TO 10                                             FL26490
      S = (A-B(I))/(B(I+1)-B(I))                                         FL26500
      RETURN                                                             FL26510
      END                                                                FL26520
      FUNCTION TAB(II,JJ,KK,WC)                                          FL26530
C                                                                        FL26540
C     ****************************************************************** FL26550
C                                                                        FL26560
C     THE INTERNAL DATA:                                                 FL26570
C                                                                        FL26580
      DIMENSION A(9,6,9),ALPHA(9,6,9),A1(5),A2(5),ALPHA1(5),             FL26590
     *    MAXI(6,9)                                                      FL26600
C                                                                        FL26610
C     A(1,J,K),J=1,3 = POWER LAW COEFFICIENT FOR THE ABSORPTION          FL26620
C     COEFICIENT FOR THE MARSHALL-PALMER WATER DROP SIZE                 FL26630
C     DISTRIBUTION FOR TEMPERATURE=10.*(J-2) AND FREQUENCY=FR(K)         FL26640
C     A(2,J,K),J=1,3 = THE SAME FOR THE EXTINCTION COEFFICIENT           FL26650
C     A(I,J,K),J=1,3,I=3,9 = THE SAME FOR THE LEGENDRE                   FL26660
C     COEFFICIENT #I-2                                                   FL26670
C     A(I,4,K),I=1,9 = THE SAME AS A(I,2,K), BUT FOR ICE                 FL26680
C     (NO TEMPERATURE DEPENDENCE)                                        FL26690
C     A(I,5,K),I=1,9 = THE SAME AS A(I,2,K), BUT FOR THE BEST DROP       FL26700
C     SIZE DISTRIBUTION (NO TEMPRATURE DEPENDENCE)                       FL26710
C     A(I,6,K),I=1,9 = THE SAME AS A(I,5,K), BUT FOR ICE                 FL26720
C     ALPHA(I,J,K) = THE POWER EXPONENET CORRESPONDING TO A(I,J,K)       FL26730
C     MAXI(J,K): TAB(I,J,K,WC)=0. IF I.GT.MAXI(J,K)                      FL26740
C     A1, A2 AND ALPHA1 = THE POWER-LINEAR LAW COEFFICIENTS AND          FL26750
C     EXPONENT FOR THE EXCEPTIONAL CASES                                 FL26760
C                                                                        FL26770
C     THE FORMULA:                                                       FL26780
C                                                                        FL26790
C     SC=A*WC**ALPHA IF ABS(A).GT.10.**-8,                               FL26800
C     SC=A1*WC**ALPHA1+A2*WC IF ABS(A).LE.10.**-8,                       FL26810
C     A1, A2 AND ALPHA1 ARE INDEXED BY INT(ALPHA)                        FL26820
C                                                                        FL26830
C     THE BLOCK-DATA SECTION                                             FL26840
C                                                                        FL26850
      DATA ((MAXI(J,K),J=1,6),K=1,9)/4*6,14*7,36*9/                      FL26860
      DATA (A1(I),A2(I),ALPHA1(I),I=1,5)/.611,-.807,1.18,.655,-.772,1.08 FL26870
     * ,.958,-1.,.99,.538,-.696,1.27,1.58,-1.50,1.02/                    FL26880
      DATA ((A(I,J,1),J=1,6),I=1,7)/.284,.285,.294,.001336,.36,.00146,   FL26890
     *.363,.365,.375,.0148,.528,.0317,3*0.,.3147,0.,.438,                FL26900
     *.4908,.487,.482,.528,.478,.538,3*.0350,.0470,.0482,.0647,          FL26910
     *.002,.00205,.00208,.00285,.0037,.0048,4*0.,.00021,.00016/          FL26920
      DATA ((ALPHA(I,J,1),J=1,6),I=1,7)/1.214,1.233,1.25,1.035,1.22,     FL26930
     *1.076,1.291,1.31,1.323,1.63,1.334,1.74,3.1,2.1,1.1,5.005,4.1,.555, FL26940
     *-.009,-.013,-.016,.028,-.019,.031,.398,.399,.4,.473,.461,.525,     FL26950
     *1.06,.97,1.03,1.03,1.18,1.16,4*0.,1.3,1.3/                         FL26960
      DATA ((A(I,J,2),J=1,6),I=1,7)/.8,.77,.73,.00344,.76,.0043,         FL26970
     *1.28,1.27,1.24,.162,1.43,.332,.254,.172,0.,.93,.32,1.29,           FL26980
     *.5,.486,.4706,.69,.481,.8,.0965,.0936,.09,.159,.151,.234,          FL26990
     *.0234,.0228,.0221,.034,.057,.065,2*.0037,.0035,.005,.011,.0106/    FL27000
      DATA ((ALPHA(I,J,2),J=1,6),I=1,7)/2*1.1,1.09,1.13,1.02,1.19,       FL27010
     *2*1.20,1.15,1.66,1.14,1.7,.29,.42,5.1,.39,.66,.44,                 FL27020
     *0.,-.01,-.0199,.12,-.01,.17,.386,.378,.2,.48,.485,.56,             FL27030
     *.92,.91,.90,.97,1.15,1.13,1.32,1.26,1.32,1.41,1.69,1.67/           FL27040
      DATA ((A(I,J,3),J=1,6),I=1,7)/1.11,1.07,1.02,.0059,.92,.00775,     FL27050
     *1.88,1.89,1.87,.43,1.80,.77,.512,.425,.336,1.25,.677,1.55,         FL27060
     *.561,.534,.506,.867,.6,1.07,.175,.165,.156,.300,.292,.49,          FL27070
     *.066,.064,.061,.105,.16,.22,                                       FL27080
     *.0169,.0162,.0156,.023,.055,.056/                                  FL27090
      DATA ((ALPHA(I,J,3),J=1,6),I=1,7)/2*1.01,1.,1.18,.92,1.23,         FL27100
     *3*1.1,1.58,1.,1.57,.264,.320,.445,.27,.416,.27,                    FL27110
     *.048,.033,.018,.168,.09,.224,.429,.417,.402,.501,.528,.62,         FL27120
     *2*.83,.82,.9,1.01,1.11,1.22,1.21,1.2,1.23,1.51,1.53/               FL27130
      DATA ((A(I,J,4),J=1,6),I=1,9)/1.51,1.49,1.44,.0163,1.12,.0194,     FL27140
     *2.73,2.77,2.79,1.61,2.18,1.9,1.14,1.054,.961,1.57,1.36,1.66,       FL27150
     *.99,.93,.87,1.31,1.33,1.63,.594,.557,.516,.77,1.02,1.16,           FL27160
     *.352,.334,.315,.43,.73,.8,.171,.163,.154,.18,.47,.43,              FL27170
     *.084,.081,.077,.106,.29,.32,.037,.036,.034,.029,.16,.11/           FL27180
      DATA ((ALPHA(I,J,4),J=1,6),I=1,9)/.87,.86,.85,1.181,.79,1.16,      FL27190
     *.93,.92,.91,1.3,.84,1.18,.188,.21,.24,.09,.21,.06,                 FL27200
     *2*.2,.19,.175,.275,.2,2*.461,.459,.39,.51,.41,                     FL27210
     *2*.66,.65,.58,.70,.64,2*.94,.93,.84,1.03,1.01,                     FL27220
     *3*1.22,1.09,1.37,1.4,1.58,1.56,1.54,1.5,1.8,1.9/                   FL27230
      DATA ((A(I,J,5),J=1,6),I=1,9)/1.55,1.53,1.49,.0194,1.14,.0225,     FL27240
     *2.82,2.87,2.90,1.91,2.22,2.,1.266,1.184,1.093,1.60,1.48,1.65,      FL27250
     *1.13,1.07,1.,1.4,1.51,1.69,.74,.698,.649,.87,1.24,1.23,            FL27260
     *.465,.444,.418,.52,.94,.91,.248,.238,.225,.24,.65,.53,             FL27270
     *.132,.128,.122,.15,.43,.47,.065,.063,.06,.045,.26,.16/             FL27280
      DATA ((ALPHA(I,J,5),J=1,6),I=1,9)/.85,.84,.83,1.168,.78,1.15,      FL27290
     *.9,.89,.88,1.23,.82,1.11,.172,.191,.216,.071,.181,.04,             FL27300
     *.222,.221,.22,.165,.274,.17,.452,.454,.456,.35,.48,.33,            FL27310
     *.63,.68,.63,.52,.66,.55,3*.89,.76,.94,.86,                         FL27320
     *1.14,1.13,1.12,.96,1.24,1.1,1.44,1.41,1.43,1.31,1.6,1.6/           FL27330
      DATA ((A(I,J,6),J=1,6),I=1,9)/2*1.58,1.54,.0248,1.15,.0279,        FL27340
     *2.94,2.97,3.,2.34,2.25,2.2,1.447,1.374,1.288,1.62,1.64,1.63,       FL27350
     *1.37,1.31,1.234,1.52,1.8,1.77,1.,.96,.898,1.01,1.6,1.3,            FL27360
     *.68,.66,.62,.66,1.31,1.07,.41,.4,.38,.33,.99,.66,                  FL27370
     *.25,.24,.23,.23,.71,.56,.136,.133,.127,.081,.49,.26/               FL27380
      DATA ((ALPHA(I,J,6),J=1,6),I=1,9)/.83,.81,.8,1.145,.762,1.120,     FL27390
     *.87,.86,.85,1.14,.799,1.,.149,.165,.184,.046,.148,.014,            FL27400
     *.232,.236,.238,.146,.255,.13,.428,.433,.438,.28,.44,.23,           FL27410
     *3*.59,.44,.59,.43,3*.81,.64,.83,.66,                               FL27420
     *1.02,2*1.01,.81,1.06,.89,2*1.25,1.24,1.07,1.36,1.3/                FL27430
      DATA ((A(I,J,7),J=1,6),I=1,9)/1.60,1.59,1.56,.0285,1.16,.0314,     FL27440
     *2.98,3.02,3.05,2.6,2.26,2.3,1.546,1.481,1.4,1.63,1.72,1.62,        FL27450
     *1.52,1.464,1.388,1.58,1.97,1.8,1.18,1.13,1.07,1.08,1.82,1.33,      FL27460
     *.84,.82,.78,.75,1.55,1.16,.54,.53,.5,.4,1.22,.74,                  FL27470
     *.34,.33,.32,.3,.93,.67,2*.2,.19,.112,.67,.33/                      FL27480
      DATA ((ALPHA(I,J,7),J=1,6),I=1,9)/.81,.80,.788,1.132,.753,1.105,   FL27490
     *.85,.84,.83,1.09,.788,.95,.136,.153,.167,.033,.131,.004,           FL27500
     *.232,.236,.241,.133,.24,.11,.411,.416,.422,.25,.40,.19,            FL27510
     *3*.56,.4,.55,.38,2*.77,.76,.58,.76,.56,                            FL27520
     *3*.95,.74,.97,.78,1.17,2*1.16,.98,1.23,1.11/                       FL27530
      DATA ((A(I,J,8),J=1,6),I=1,9)/2*1.60,1.58,.045,1.15,.0461,         FL27540
     *3.08,3.09,3.1,3.3,2.27,2.32,1.849,1.81,1.75,1.628,1.98,1.606,      FL27550
     *2.07,2.04,1.98,1.78,2.5,1.946,1.89,1.86,1.81,1.30,2.6,1.508,       FL27560
     *1.58,1.56,1.52,1.11,2.49,1.57,1.22,1.21,1.18,.68,2.2,1.11,         FL27570
     *2*.91,.89,.61,2.,1.18,2*.65,.64,.299,1.6,.73/                      FL27580
      DATA ((ALPHA(I,J,8),J=1,6),I=1,9)/.777,.764,.752,1.092,.729,1.057, FL27590
     *.796,.79,.784,.96,.756,.81,.1,.108,.117,.004,.089,-.006,           FL27600
     *.207,.210,.215,.093,.182,.075,2*.34,.35,.15,.30,.122,              FL27610
     *3*.46,.3,.41,.28,3*.61,.42,.55,.394,                               FL27620
     *3*.75,.56,.7,.55,2*.91,.9,.76,.87,.79/                             FL27630
      DATA ((A(I,J,9),J=1,6),I=1,9)/2*1.58,1.56,.0587,1.13,.0579,        FL27640
     *3.09,2*3.08,3.39,2.26,2.33,2.009,1.99,1.95,1.624,2.11,1.64,        FL27650
     *2.43,2.42,2.38,1.902,2.80,2.078,2*2.42,2.38,1.454,3.09,1.7,        FL27660
     *2*2.2,2.17,1.4,3.1,1.91,1.87,1.88,1.85,.94,3.,1.46,                FL27670
     *2*1.54,1.52,.93,2.8,1.64,2*1.22,1.21,.53,2.5,1.17/                 FL27680
      DATA ((ALPHA(I,J,9),J=1,6),I=1,9)/.757,.746,.736,1.06,.717,1.024,  FL27690
     *.766,.764,.761,.86,.74,.763,.084,.087,.092,-.0018,.069,.007,       FL27700
     *.183,.182,.184,.078,.148,.075,3*.29,.128,.24,.13,                  FL27710
     *.4,2*.39,.264,.33,.256,2*.52,.51,.367,.44,.360,                    FL27720
     *2*.63,.62,.49,.55,.47,.76,2*.75,.67,.67,.66/                       FL27730
C                                                                        FL27740
C                                                                        FL27750
      IF (II.GT.MAXI(JJ,KK)) THEN                                        FL27760
         TAB = 0.                                                        FL27770
         RETURN                                                          FL27780
      ENDIF                                                              FL27790
      IF (ABS(A(II,JJ,KK)).GT.1.E-8) THEN                                FL27800
         TAB = A(II,JJ,KK)*WC**ALPHA(II,JJ,KK)                           FL27810
      ELSE                                                               FL27820
         L = ALPHA(II,JJ,KK)                                             FL27830
         TAB = A1(L)*WC**ALPHA1(L)+A2(L)*WC                              FL27840
      ENDIF                                                              FL27850
      RETURN                                                             FL27860
      END                                                                FL27870
      FUNCTION RAYSCT(V)                                                 FL27880
C                                                                        FL27890
C     RADIATION FLD OUT                                                  FL27900
C     **  MOLECULAR SCATTERING                                           FL27910
C                                                                        FL27920

      IF (V.LE.3000.) then
         RAYSCT = 0.                                                       
      else
c
c     The following statement previosly operative in lbllow.f has  been  ! sac 06/06/02
c     replaced by the expression implemented in SUBROUTINE CONTNM in 
c     module contnm.f.  
c
c***      RAYSCT = V**3/(9.26799E+18-1.07123E+09*V**2)  
c
c     This formulation, adopted from MODTRAN_3.5 (using approximation
c     of Shettle et al., (Appl Opt, 2873-4, 1980) with depolarization
c     = 0.0279, output in km-1 for T=273K & P=1 ATM) has been used.
c
         RAYSCT = V**3/(9.38076E18-1.08426E09*V**2)
C                                                                        FL27960
C     V**4 FOR RADIATION FLD IN                                          FL27970
C                                                                        FL27980
      endif
c
      RETURN                                                             FL27990
      END                                                                FL28000
      FUNCTION TNRAIN(RR,V,TM,RADFLD)                                    FL28010
C                                                                        FL28020
C     CC                                                                 FL28030
C                                                                        FL28040
      COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                      FL28050
      COMMON /LCRD3/ H1,H2,ANGLE,RANGE,BETA,RE,LEN                       FL28060
C                                                                        FL28070
C     CC   CALCULATES TRANSMISSION DUE TO RAIN AS A FUNCTION OF          FL28080
C     CC   RR=RAIN RATE IN MM/HR                                         FL28090
C     CC   OR WITHIN 350CM-1 USES THE MICROWAVE TABLE ROUTINE TO         FL28100
C     CC   OBTAIN THE EXTINCTION DUE TO RAIN                             FL28110
C     CC   RANGE=SLANT RANGE KM                                          FL28120
C     CC                                                                 FL28130
C     CC   ASSUMES A MARSHALL-PALMER RAIN DROP SIZE DISTRIBUTION         FL28140
C     CC   N(D)=NZERO*EXP(-A*D)                                          FL28150
C     CC   NZERO=8.E3 (MM-1)  (M-3)                                      FL28160
C     CC   A=41.*RR**(-0.21)                                             FL28170
C     CC   D=DROP DIAMETER (CM)                                          FL28180
C     CC                                                                 FL28190
C                                                                        FL28200
      REAL*8           V
c
      REAL NZERO                                                         FL28210
      DATA NZERO /8000./                                                 FL28220
C                                                                        FL28230
C     CC                                                                 FL28240
C                                                                        FL28250
      freq = V
c
      A = 41./RR**0.21                                                   FL28260
C                                                                        FL28270
C     CC                                                                 FL28280
C                                                                        FL28290
      IF (RR.LE.0) TNRAIN = 0.                                           FL28300
      IF (RR.LE.0) RETURN                                                FL28310
C                                                                        FL28320
C     CC                                                                 FL28330
C                                                                        FL28340
      IF (FREQ.GE.350.0) THEN                                               FL28350
         TNRAIN = PI*NZERO/A**3                                          FL28360
         TNRAIN = TNRAIN/RADFLD                                          FL28370
      ELSE                                                               FL28380
         TNRAIN = GMRAIN(FREQ,TM,RR)                                        FL28390
      ENDIF                                                              FL28400
      RETURN                                                             FL28410
      END                                                                FL28420
C
C     *****************************************************************
C
      SUBROUTINE LAYVSA(K,RH,AHAZE,IHAZ1,ZSTF)                            FL28430
C                                                                        FL28440
C     RETURNS HAZE FOR VSA OPTION                                        FL28450
C                                                                        FL28460
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,N_LVL,JH1              FL28470
      COMMON/MODEL/ ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),                     FL28480
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)          FL28490
C
      COMMON /LCRD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB   FL28500
      COMMON /LCRD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,    FL28510
     *     RAINRT                                                        FL28520
c     COMMON /MDATA/ ZDA(MXZMD),P(MXZMD),T(MXZMD),WH(MXZMD),WO(MXZMD),   FL28530
c    *     HMIX(MXZMD),CDUM1(MXZMD,7),RDUM2(MXZMD,7)                     FL28540
      COMMON /MDATA/                              WH(MXZMD),WO(MXZMD),   FL28530
     *                 CDUM1(MXZMD,7),RDUM2(MXZMD,7)                     FL28540
      COMMON /MDATA2/ZDA(MXZMD),P(MXZMD),T(MXZMD)
      COMMON /ZVSALY/ ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)             FL28550
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL28560
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL28570
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL28580
C
      DIMENSION ZSTF(MXZMD)                                              FL28590
C
      RH = 0.                                                            FL28600
      AHAZE = 0                                                          FL28610
      IHAZ1 = 0                                                           FL28620
      IF (MODEL.EQ.0) GO TO 10                                           FL28630
      IF (K.GT.9) RETURN                                                 FL28640
      ZMDL(K) = ZVSA(K)                                                  FL28650
      RH = RHVSA(K)                                                      FL28660
      AHAZE = AHVSA(K)                                                   FL28670
      IHAZ1 = IHVSA(K)                                                    FL28680
      RETURN                                                             FL28690
C                                                                        FL28700
C     MODEL 7 CODEING                                                    FL28710
C     OLD LAYERS  AEROSOL RETURNED                                       FL28720
C                                                                        FL28730
   10 CONTINUE                                                           FL28740
      ZVSA(10) = ZVSA(9)+0.01                                            FL28750
      RHVSA(10) = 0.                                                     FL28760
      AHVSA(10) = 0.                                                     FL28770
      IHVSA(10) = 0                                                      FL28780
C                                                                        FL28790
C     JML=ML                                                             FL28800
C                                                                        FL28810
      IF (ML.EQ.1) WRITE (IPR,900)                                       FL28820
      IF (ML.EQ.1) RETURN                                                FL28830
      IF (ZSTF(K).GT.ZVSA(10)) RETURN                                    FL28840
      DO 20 JJ = 1, 9                                                    FL28850
         JL = JJ                                                         FL28860
         IF (ZSTF(K).LT.ZVSA(JJ)) GO TO 20                               FL28870
         JN = JJ+1                                                       FL28880
         IF (ZSTF(K).LT.ZVSA(JN)) GO TO 30                               FL28890
   20 CONTINUE                                                           FL28900
      JN = 10                                                            FL28910
   30 CONTINUE                                                           FL28920
      DIF = ZVSA(JN)-ZVSA(JL)                                            FL28930
      DZ = ZVSA(JN)-ZSTF(K)                                              FL28940
      DLIN = DZ/DIF                                                      FL28950
      IHAZ1 = IHVSA(JL)                                                   FL28960
C                                                                        FL28970
C     FAC=(ZVSA(JL)-ZSTF  ( K))/DIF                                      FL28980
C                                                                        FL28990
      AHAZE = (AHVSA(JN)-AHVSA(JL))*DLIN+AHVSA(JL)                       FL29000
      RETURN                                                             FL29010
C                                                                        FL29020
  900 FORMAT('   ERROR MODEL EQ 0 AND ARMY MODEL CANNOT MIX')            FL29030
C                                                                        FL29040
      END                                                                FL29050
C
C     ******************************************************************
C
      SUBROUTINE STDMDL                                                  FL29060
C                                                                        FL29070
C     ****************************************************************** FL29080
C     LOADS DENSITIES INTO COMMON MODEL AND                              FL29090
C     CALCULATES THE INDEX OF REFRACTION                                 FL29100
C                                                                        FL29110
C     AERSOLS NOW LOADED IN AERNSM                                       FL29120
C                                                                        FL29130
C     ZM COMMON /MODEL/ FINAL ALTITUDE FOR LOWTRAN                       FL29140
C     Z COMMON /MDATA/  ALTITUDE FOR DATA IN MDATA                       FL29150
C     ZN  BLANK COMMON                                                   FL29160
C     ZP  BLANK COMMON                                                   FL29170
C                                                                        FL29180
C     ****************************************************************** FL29190
C                                                                        FL29200
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),VH(16),TX(16),W(16)
      COMMON WPATH(IM2,16),TBBY(IM2)
      COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),VX2(47),AWCCON(5)
C
      CHARACTER*8      HMOD
C
      COMMON /CMN/ HMOD(3),ZN(MXZMD),PN(MXZMD),TN(MXZMD),RFNDXM(MXZMD),
     *         ZP(IM2),PP1(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),
     *          TPSUM(IM2),RHOPSM(IM2),IMMAX,WGM(MXZMD),DEMW(MXZMD),
     *          AMTP(MXMOL,MXPDIM)                                        
C
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL29210
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL29220
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL29230
      COMMON /LCRD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB   FL29240
      COMMON /LCRD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,    FL29250
     *     RAINRT                                                        FL29260
      COMMON /LCRD3/ H1,H2,ANGLE,RANGE,BETA,RE,LEN                       FL29270
      COMMON /LCRD4/ V1,V2,DV                                            FL29280
c     COMMON /MDATA/ ZMDL(MXZMD),P(MXZMD),T(MXZMD),WH(MXZMD),WO(MXZMD),  FL29290
c    *     HMIX(MXZMD),CLD(MXZMD,7),RR(MXZMD,7)                          FL29300
      COMMON /MDATA/                               WH(MXZMD),WO(MXZMD),  FL29290
     *                 CLD(MXZMD,7),RR(MXZMD,7)                          FL29300
      COMMON /MDATA2/ZMDL(MXZMD),P(MXZMD),T(MXZMD)
      COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                      FL29310
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,N_LVL,JH1              FL29320
      COMMON /MODEL/ ZM(MXZMD),PM(MXZMD),TM(MXZMD),
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),
     *     CLDM(MXZMD),RRM(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)
C                                                                        FL29430
C     XLOSCH = LOSCHMIDT'S NUMBER,MOLECULES CM-2,KM-1                    FL29440
C                                                                        FL29450
      DATA PZERO /1013.25/,TZERO/273.15/,XLOSCH/2.6868E24/               FL29460
C                                                                        FL29470
C     RV GAS CONSTANT FOR WATER IN MB/(GM M-3 K)                         FL29480
C     CON CONVERTS WATER VAPOR FROM GM M-3 TO MOLECULES CM-2 KM-1        FL29490
C                                                                        FL29500
      DATA RV/4.6152E-3/,CON/3.3429E21/                                  FL29510
C                                                                        FL29520
C     CONSTANTS FOR INDEX OF REFRACTION, AFTER EDLEN, 1965               FL29530
C                                                                        FL29540
      DATA A0/83.42/,A1/185.08/,A2/4.11/,                                FL29550
     *     B1/1.140E5/,B2/6.24E4/,C0/43.49/,C1/1.70E4/                   FL29560
C                                                                        FL29570
C     F(A) IS SATURATED WATER WAPOR DENSITY AT TEMP T,A=TZERO/T          FL29580
C                                                                        FL29590
      F(A) = EXP(18.9766-14.9595*A-2.43882*A*A)*A                        FL29600
C                                                                        FL29610
C     H20 CONTINUUM IS STORED AT 296 K RHZERO IS AIR DENSITY AT 296 K    FL29620
C     IN UNITS OF LOSCHMIDT'S                                            FL29630
C                                                                        FL29640
C     CALL DRYSTR                                                        FL29650
C                                                                        FL29660
      RHZERO = (273.15/296.0)                                            FL29670
C                                                                        FL29680
      IF (ICLD.GT.20) THEN                                               FL29690
         WRITE (IPR,900) ICLD                                            FL29700
         STOP 'STDMDL: ICLD GT 20'                                       FL29710
      ENDIF                                                              FL29720
C                                                                        FL29730
C     LOAD ATMOSPHERE PROFILE INTO /MODEL/                               FL29740
C                                                                        FL29750
      IF (M.LT.7) ML = NL                                                FL29760
      DO 10 I = 1, ML                                                    FL29770
         IF (M.NE.7) ZM(I) = ZMDL(I)                                     FL29780
         PM(I) = P(I)                                                    FL29790
         TM(I) = T(I)                                                    FL29800
         PP = PM(I)                                                      FL29810
         TT = TM(I)                                                      FL29820
         F1 = (PP/PZERO)/(TT/TZERO)                                      FL29830
         F2 = (PP/PZERO)*SQRT(TZERO/TT)                                  FL29840
         WTEMP = WH(I)                                                   FL29850
         RELHUM(I) = 0.                                                  FL29860
C                                                                        FL29870
C        RELHUM IS CALCULATED ONLY FOR THE BOUNDRY LAYER (0 TO 2 KM)     FL29880
C                                                                        FL29890
C        SCALED H2O DENSITY                                              FL29900
C                                                                        FL29910
         DENSTY(1,I) = 0.1*WTEMP*F2**0.9                                 FL29920
C                                                                        FL29930
C        C    IF (ZM(I).GT.6.0) GO TO 15                                 FL29940
C        C    IF(DENSTY(7,I).LE.0.) GO TO 15                             FL29950
C                                                                        FL29960
         TS = TZERO/TT                                                   FL29970
         RELHUM(I) = 100.0*(WTEMP/F(TS))                                 FL29980
C                                                                        FL29990
C        UNIFORMALY MIXED GASES DENSITYS                                 FL30000
C                                                                        FL30010
         DENSTY(2,I) = F1*F2**0.75                                       FL30020
C                                                                        FL30030
C        UV OZONE                                                        FL30040
C                                                                        FL30050
         DENSTY(8,I) = 46.6667*WO(I)                                     FL30060
C                                                                        FL30070
C        IR OZONE                                                        FL30080
C                                                                        FL30090
         DENSTY(3,I) = DENSTY(8,I)*F2**0.4                               FL30100
C                                                                        FL30110
C        N2 CONTINUUM                                                    FL30120
C                                                                        FL30130
         DENSTY(4,I) = 0.8*F1*F2                                         FL30140
C                                                                        FL30150
C        SELF BROADENED WATER                                            FL30160
C                                                                        FL30170
         RHOAIR = F1                                                     FL30180
         RHOH2O = CON*WTEMP/XLOSCH                                       FL30190
         RHOFRN = RHOAIR-RHOH2O                                          FL30200
         DENSTY(5,I) = XLOSCH*RHOH2O**2/RHZERO                           FL30210
C                                                                        FL30220
C        FOREIGN BROADENED                                               FL30230
C                                                                        FL30240
         DENSTY(10,I) = XLOSCH*RHOH2O*RHOFRN/RHZERO                      FL30250
C                                                                        FL30260
C        MOLECULAR SCATTERING                                            FL30270
C                                                                        FL30280
         DENSTY(6,I) = F1                                                FL30290
C                                                                        FL30300
C        AEROSOL FOR 0 TO 2KM                                            FL30310
C                                                                        FL30320
C                                                                        FL30330
C        RELITIVE HUMIDITY WEIGHTED BY BOUNDRY LAYER AEROSOL (0 TO 2 KM) FL30340
C                                                                        FL30350
         RELH = RELHUM(I)                                                FL30360
         RELH = MIN(RELH,99.)                                            FL30370
         RHLOG =  LOG(100.-RELH)                                         FL30380
C                                                                        FL30390
C        DENSTY(15,I)=RELHUM(I)*DENSTY(7,I)                              FL30400
C                                                                        FL30410
         DENSTY(15,I) = RHLOG*DENSTY(7,I)                                FL30420
C                                                                        FL30430
C        DENSITY (9,I) NO LONGER USED                                    FL30440
C                                                                        FL30450
         DENSTY(9,I) = 0.                                                FL30460
C                                                                        FL30470
C        IF(ICH(1).GT.7) DENSTY(15,I)=RELHUM(I)*DENSTY(12,I)             FL30480
C                                                                        FL30490
         IF (ICH(1).GT.7) DENSTY(15,I) = RHLOG*DENSTY(12,I)              FL30500
C                                                                        FL30510
C        HNO3 IN ATM * CM /KM                                            FL30520
C        DENSTY(11,I)= F1* HMIX(I)*1.0E-4                                FL30530
C                                                                        FL30540
         DENSTY(11,I) = 0.                                               FL30550
C                                                                        FL30560
C        IF(MODEL.EQ.0) DENSTY(11,I)=F1*HSTOR(I)*1.0E-4                  FL30570
C        CIRRUS CLOUD                                                    FL30580
C                                                                        FL30590
         IF (ICLD.LT.18) DENSTY(16,I) = 0.0                              FL30600
C                                                                        FL30610
C        RFNDX = REFRACTIVITY 1-INDEX OF REFRACTION                      FL30620
C        FROM EDLEN, 1966                                                FL30630
C                                                                        FL30640
         PPW = RV*WTEMP*TT                                               FL30650
         AVW = 0.5*(V1+V2)                                               FL30660
         RFNDX(I) = ((A0+A1/(1.-(AVW/B1)**2)+A2/(1.0-(AVW/B2)**2))*(PP/  FL30670
     *      PZERO)*(TZERO+15.0)/TT-(C0-(AVW/C1)**2)*PPW/PZERO)*1.E-6     FL30680
   10 CONTINUE                                                           FL30690
      WRITE (IPR,910)                                                    FL30700
      ZERO = 0.                                                          FL30710
      DO 20 I = 1, ML                                                    FL30720
         WRITE (IPR,905) I,ZM(I),PM(I),TM(I),ZERO,ZERO,DENSTY(7,I),      FL30730
     *      DENSTY(12,I),DENSTY(13,I),DENSTY(14,I),DENSTY(15,I),         FL30740
     *      DENSTY(16,I),RELHUM(I)                                       FL30750
   20 CONTINUE                                                           FL30760
      RETURN                                                             FL30770
C                                                                        FL30780
  900 FORMAT('1',//10X,'ICLD  CANNOT BE GREATER THAN 20 BUT IS',         FL30790
     * I5,//)                                                            FL30800
  905 FORMAT (I4,0PF9.2,F9.3,F7.1,1X,1P9E10.3)                           FL30810
  910 FORMAT('1',/,'  ATMOSPHERIC PROFILES',//,                          FL30820
     * 3X,'I',T10,'Z',T18,'P',T26,'T',T33,'CNTMFRN',T45,'HNO3',          FL30830
     * T53,'AEROSOL 1',T63,'AEROSOL 2', T73,'AEROSOL 3',T83,             FL30840
     * 'AEROSOL 4',T93,'AER1*RH',T103,'CIRRUS',T118,'RH'/,               FL30850
     * T9,'(KM)',T17,'(MB)',T25,'(K)',T31,'MOL/CM2 KM',T42,              FL30860
     * 'ATM CM/KM',T54,'(-)',T64,'(-)',T74,'(-)',T84,'(-)',T94,          FL30870
     * '(-)',T104,'(-)',T113,'(PERCNT)',/)                               FL30880
C                                                                        FL30890
      END                                                                FL30900
C
C     *****************************************************************
C
      SUBROUTINE NEWMDL(MAXATM)                                          FL30910
C                                                                        FL30920
C     CC                                                                 FL30930
C     CC   ROUTINE TO COMBINE LOWTRAN AND LBLRTM LAYERING                FL30940
C                                                                        FL30950
C     ZMTP STORES ZM VALUES                                              FL30960
C     ZOUT COMMON /ZOUTP/ FINAL LBLRTM BOUNDRIES                         FL30970
C     ZMDL COMMON /MODEL/ FINAL ALTITUDE FOR LOWTRAN                     FL30980
C     CC                                                                 FL30990
C                                                                        FL31000
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL31010
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL31020
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL31030
      COMMON /LCRD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB   FL31040
      COMMON/MODEL/ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),                      FL31050
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL31060
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)          FL31070
      COMMON /ZOUTP/ ZOUT(MXLAY),SOUT(MXLAY),RHOSUM(MXLAY),
     *     AMTTOT(MXMOL),AMTCUM(MXMOL),ISKIP(MXMOL)                      FL31080
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISS,N_LVL,JH1                 FL31090
C
      DIMENSION PTMP(MXZMD),TTMP(MXZMD),RTMP(MXZMD),                     FL31100
     *     DENTMP(16,MXZMD),ZMTP(MXZMD),RRAMTJ(MXZMD)                    FL31110
C
      DO 10 I = 1, ML                                                    FL31120
         ZMTP(I) = ZMDL(I)                                               FL31130
         PTMP(I) = PM(I)                                                 FL31140
         TTMP(I) = TM(I)                                                 FL31150
         RTMP(I) = RFNDX(I)                                              FL31160
         RRAMTJ(I) = RRAMT(I)                                            FL31170
         DO  8 K = 1, 16         
            DENTMP(K,I) = DENSTY(K,I)                                    FL31190
 8       CONTINUE      
 10   CONTINUE                                                           FL31200
      IF (ITYPE.EQ.1) GO TO 130                                          FL31210
      IF (ML.LT.2) GO TO 130                                             FL31220
      DO 20 I = 1, N_LVL                                                   FL31230
         DO 18 K = 1, 16 
            DENSTY(K,I) = 0.                                             FL31250
 18      CONTINUE     
 20   CONTINUE                                                           FL31260
      I = 1                                                              FL31270
      L = 1                                                              FL31280
      J1 = 1                                                             FL31290
   30 DO 80 J = J1, ML                                                   FL31300
         IF (ZMDL(J).LT.ZOUT(1)) GO TO 80                                FL31310
         IF (ZMDL(J).LE.ZOUT(I)) GO TO 40                                FL31320
         GO TO 60                                                        FL31330
   40    PM(L) = PTMP(J)                                                 FL31340
         TM(L) = TTMP(J)                                                 FL31350
         RFNDX(L) = RTMP(J)                                              FL31360
         RRAMT(L) = RRAMTJ(J)                                            FL31370
         ZMTP(L) = ZMDL(J)                                               FL31380
         DO 50 K = 1, 16                                                 FL31390
            DENSTY(K,L) = DENTMP(K,J)                                    FL31400
   50    CONTINUE                                                        FL31410
         L = L+1                                                         FL31420
         IF (L.GT.MAXATM) GO TO 100                                      FL31430
         J1 = J+1                                                        FL31440
         IF (ZMDL(J).LT.ZOUT(I)) GO TO 80                                FL31450
         GO TO 90                                                        FL31460
   60    JL = J-1                                                        FL31470
         IF (JL.LT.1) JL = 1                                             FL31480
         JP = JL+1                                                       FL31490
         DIF = ZMDL(JP)-ZMDL(JL)                                         FL31500
         DZ = ZOUT(I)-ZMDL(JL)                                           FL31510
         DLIN = DZ/DIF                                                   FL31520
         PM(L) = (PTMP(JP)-PTMP(JL))*DLIN+PTMP(JL)                       FL31530
         TM(L) = (TTMP(JP)-TTMP(JL))*DLIN+TTMP(JL)                       FL31540
         RFNDX(L) = (RTMP(JP)-RTMP(JL))*DLIN+RTMP(JL)                    FL31550
         RRAMT(L) = (RRAMTJ(JP)-RRAMTJ(JL))*DLIN+RRAMTJ(JL)              FL31560
         ZMTP(L) = ZOUT(I)                                               FL31570
         DO 70 K = 1, 16                                                 FL31580
            DENSTY(K,L) = (DENTMP(K,JP)-DENTMP(K,JL))*DLIN+DENTMP(K,JL)  FL31590
   70    CONTINUE                                                        FL31600
         L = L+1                                                         FL31610
         IF (L.GT.MAXATM) GO TO 100                                      FL31620
         GO TO 90                                                        FL31630
   80 CONTINUE                                                           FL31640
   90 IF (I.EQ.N_LVL) GO TO 110                                            FL31650
      I = I+1                                                            FL31660
      GO TO 30                                                           FL31670
C                                                                        FL31680
C     CC                                                                 FL31690
C     CC    SET LOWTRAN HEIGHTS TO FINAL COMBINED LAYERING OF LBL/LOW    FL31700
C     CC    SET ML TO THE FINAL COUNT OF COMBINED LAYERING               FL31710
C     CC                                                                 FL31720
C                                                                        FL31730
  100 WRITE (IPR,900)                                                    FL31740
      STOP 'NEWMDL; LAYER LIMIT'                                         FL31750
  110 LM = L-1                                                           FL31760
      DO 120 I = 1, LM                                                   FL31770
         ZMDL(I) = ZMTP(I)                                               FL31780
  120 CONTINUE                                                           FL31790
      ML = LM                                                            FL31800
  130 RETURN                                                             FL31810
C                                                                        FL31820
  900 FORMAT(' LAYER LIMIT REACHED  CHANGE AVARAT  2. 10. 20. WORKS' )   FL31830
C                                                                        FL31840
      END                                                                FL31850
C
C     ******************************************************************
C
      SUBROUTINE AERPRF (I,K,VIS,HAZE,IHAZE,ICLD,ISEASN,IVULCN,N)        FL31860
C                                                                        FL31870
C     ****************************************************************** FL31880
C     WILL COMPUTE DENSITY    PROFILES FOR AEROSOLS                      FL31890
C     ****************************************************************** FL31900
C                                                                        FL31910
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
      COMMON/PRFD  / ZHT(34),HZ2K(34,5),FAWI50(34),FAWI23(34),           FL31920
     *     SPSU50(34),SPSU23(34),BASTFW(34),VUMOFW(34),HIVUFW(34),       FL31930
     *     EXVUFW(34),BASTSS(34),VUMOSS(34),HIVUSS(34),EXVUSS(34),       FL31940
     *     UPNATM(34),VUTONO(34),VUTOEX(34),EXUPAT(34)
      COMMON /MODEL/ ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),                    FL31950
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL31960
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)          FL31970
      DIMENSION VS(5)                                                    FL31980
      DATA VS/50.,23.,10.,5.,2./                                         FL31990
      DATA CULWC/7.683E-03/,ASLWC/4.509E-03/,STLWC/5.272E-03/            FL32000
      DATA SCLWC/4.177E-03/,SNLWC/7.518E-03/                             FL32010
      HAZE = 0.                                                          FL32020
      N = 7                                                              FL32030
      IF (IHAZE.EQ.0) THEN                                               FL32040
         IF (ICLD.EQ.0.OR.ICLD.EQ.20) RETURN                             FL32050
      ENDIF                                                              FL32060
      IF (ZHT(I).GT.2.0) GO TO 30                                        FL32070
      DO 10 J = 2, 5                                                     FL32080
         IF (VIS.GE.VS(J)) GO TO 20                                      FL32090
   10 CONTINUE                                                           FL32100
      J = 5                                                              FL32110
   20 CONST = 1./(1./VS(J)-1./VS(J-1))                                   FL32120
      HAZE = CONST*((HZ2K(I,J)-HZ2K(I,J-1))/VIS+HZ2K(I,J-1)/VS(J)-HZ2K(I FL32130
     *   ,J)/VS(J-1))                                                    FL32140
   30 IF (ICLD.GE.1.AND.ICLD.LE.11) GO TO 40                             FL32150
      IF (ZHT(I).GT.2.0) GO TO 100                                       FL32160
      RETURN                                                             FL32170
   40 IF (CLDAMT(K).LE.0.) GO TO 100                                     FL32180
      IH = ICLD                                                          FL32190
      IF (CLDAMT(K).GT.0.0) N = 12                                       FL32200
      GO TO (50,60,70,80,90,70,90,90,50,50,50), IH                       FL32210
   50 HAZEC(K) = CLDAMT(K)/CULWC                                         FL32220
      IF (ZHT(I).GT.2.0) GO TO 100                                       FL32230
      RETURN                                                             FL32240
   60 HAZEC(K) = CLDAMT(K)/ASLWC                                         FL32250
      IF (ZHT(I).GT.2.0) GO TO 100                                       FL32260
      RETURN                                                             FL32270
   70 HAZEC(K) = CLDAMT(K)/STLWC                                         FL32280
      IF (ZHT(I).GT.2.0) GO TO 100                                       FL32290
      RETURN                                                             FL32300
   80 HAZEC(K) = CLDAMT(K)/SCLWC                                         FL32310
      IF (ZHT(I).GT.2.0) GO TO 100                                       FL32320
      RETURN                                                             FL32330
   90 HAZEC(K) = CLDAMT(K)/SNLWC                                         FL32340
      IF (ZHT(I).GT.2.0) GO TO 100                                       FL32350
      RETURN                                                             FL32360
  100 IF (ZHT(I).GT.10.) GO TO 140                                       FL32370
      IF (ICLD.GE.1.AND.ICLD.LE.11) THEN                                 FL32380
         N = 13                                                          FL32390
      ELSE                                                               FL32400
         N = 12                                                          FL32410
      ENDIF                                                              FL32420
      CONST = 1./(1./23.-1./50.)                                         FL32430
      IF (ISEASN.GT.1) GO TO 120                                         FL32440
      IF (VIS.LE.23.) HAZI = SPSU23(I)                                   FL32450
      IF (VIS.LE.23.) GO TO 260                                          FL32460
      IF (ZHT(I).GT.4.0) GO TO 110                                       FL32470
      HAZI = CONST*((SPSU23(I)-SPSU50(I))/VIS+SPSU50(I)/23.-SPSU23(I)/   FL32480
     *   50.)                                                            FL32490
      GO TO 260                                                          FL32500
  110 HAZI = SPSU50(I)                                                   FL32510
      GO TO 260                                                          FL32520
  120 IF (VIS.LE.23.) HAZI = FAWI23(I)                                   FL32530
      IF (VIS.LE.23.) GO TO 260                                          FL32540
      IF (ZHT(I).GT.4.0) GO TO 130                                       FL32550
      HAZI = CONST*((FAWI23(I)-FAWI50(I))/VIS+FAWI50(I)/23.-FAWI23(I)/   FL32560
     *   50.)                                                            FL32570
      GO TO 260                                                          FL32580
  130 HAZI = FAWI50(I)                                                   FL32590
      GO TO 260                                                          FL32600
  140 IF (ZHT(I).GT.30.0) GO TO 240                                      FL32610
      IF (ICLD.GE.1.AND.ICLD.LE.11) THEN                                 FL32620
         N = 14                                                          FL32630
      ELSE                                                               FL32640
         N = 13                                                          FL32650
      ENDIF                                                              FL32660
      HAZI = BASTSS(I)                                                   FL32670
      IF (ISEASN.GT.1) GO TO 190                                         FL32680
      IF (IVULCN.EQ.0) HAZI = BASTSS(I)                                  FL32690
      IF (IVULCN.EQ.0) GO TO 260                                         FL32700
      GO TO (150,160,170,170,160,160,170,180), IVULCN                    FL32710
  150 HAZI = BASTSS(I)                                                   FL32720
      GO TO 260                                                          FL32730
  160 HAZI = VUMOSS(I)                                                   FL32740
      GO TO 260                                                          FL32750
  170 HAZI = HIVUSS(I)                                                   FL32760
      GO TO 260                                                          FL32770
  180 HAZI = EXVUSS(I)                                                   FL32780
      GO TO 260                                                          FL32790
  190 IF (IVULCN.EQ.0) HAZI = BASTFW(I)                                  FL32800
      IF (IVULCN.EQ.0) GO TO 260                                         FL32810
      GO TO (200,210,220,220,210,210,220,230), IVULCN                    FL32820
  200 HAZI = BASTFW(I)                                                   FL32830
      GO TO 260                                                          FL32840
  210 HAZI = VUMOFW(I)                                                   FL32850
      GO TO 260                                                          FL32860
  220 HAZI = HIVUFW(I)                                                   FL32870
      GO TO 260                                                          FL32880
  230 HAZI = EXVUFW(I)                                                   FL32890
      GO TO 260                                                          FL32900
  240 N = 14                                                             FL32910
      IF (IVULCN.GT.1) GO TO 250                                         FL32920
      HAZI = UPNATM(I)                                                   FL32930
      GO TO 260                                                          FL32940
  250 HAZI = VUTONO(I)                                                   FL32950
  260 IF (HAZI.GT.0) HAZE = HAZI                                         FL32960
      END                                                                FL32970
C
C     ******************************************************************
C
      SUBROUTINE GEO(IERROR,BENDNG,MAXGEO)                               FL32980
C                                                                        FL32990
C     ****************************************************************** FL33000
C     THIS SUBROUTINE SERVES AS AN INTERFACE BETWEEN THE MAIN            FL33010
C     LOWTRAN PROGRAM 'LOWTRN' AND THE NEW SET OF SUBROUTINES,           FL33020
C     INCLUDING 'EXPINT', 'FINDSH', 'SCALHT', 'RFPATL', 'FILL',          FL33030
C     AND 'LAYER',  WHICH CALCULATE THE ABSORBER                         FL33040
C     AMOUNTS FOR A REFRACTED PATH THROUGH THE ATMOSPHERE.               FL33050
C     THE INPUT PARAMETERS ITYPE, H1, H2, ANGLE, RANGE, BETA, AND LEN    FL33060
C     ALL FUNCTION IN THE SAME WAY IN THE NEW ROUTINES AS IN THE OLD.    FL33070
C     ****************************************************************** FL33080
C                                                                        FL33090
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
      PARAMETER (MXZ20 = MXZMD+20, MX2Z3 = 2*MXZMD+3)
C
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),VH(16),TX(16),W(16)       FL33100
      COMMON WPATH(IM2,16),TBBY(IM2)                                     FL33110
      COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),VX2(47),AWCCON(5)          FL33120
C
      CHARACTER*8      HMOD                                              FL33130
C
      COMMON /CMN/ HMOD(3),ZN(MXZMD),PN(MXZMD),TN(MXZMD),RFNDXM(MXZMD),
     *          ZP(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),
     *          TPSUM(IM2),RHOPSM(IM2),IMMAX,WGM(MXZMD),DEMW(MXZMD),
     *          AMTP(MXMOL,MXPDIM)                                        
C
C     RFRPTH is dependent upon MXZMD (MXZ20=MXZMD+20;MX2Z3=2*MXZMD+3)
C
      COMMON  /RFRPTH/ ZL(MXZ20),PL(MXZ20),TL(MXZ20),RFNDXL(MXZ20),      FL33180
     *     SL(MXZ20),PPSUML(MXZ20),TPSUML(MXZ20),RHOSML(MXZ20),          FL33190
     *     DENL(16,MXZ20),AMTL(16,MXZ20),LJ(MX2Z3)                       FL33200
      COMMON /RAIN/ RNPATH(IM2),RRAMTK(IM2)                              FL33210
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL33220
     *NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                                 FL33230
     *                     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4       FL33240
      COMMON /LCRD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB   FL33250
      COMMON /LCRD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,    FL33260
     *    RAINRT                                                         FL33270
      COMMON /LCRD3/ H1,H2,ANGLE,RANGE,BETA,REE,LEN                      FL33280
      COMMON /LCRD4/ V1,V2,DV                                            FL33290
      COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                      FL33300
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,N_LVL,JH1              FL33310
      COMMON/MODEL/ ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),                     FL33320
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL33330
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)
      COMMON /PARMLT/ RE,DELTAS,ZMAX,IMAX,IMOD,IBMAX,IPATH               FL33340
      COMMON /ADRIVE/LOWFLG,IREAD,MODELF,ITYPEF,NOZERO,NOPRNF,           FL33350
     * H1F,H2F,ANGLEF,RANGEF,BETAF,LENF,VL1,VL2,RO,IPUNCH,VBAR,          FL33360
     * HMINF,PHIF,IERRF,HSPACE                                           FL33370
      DIMENSION KMOL(16)                                                 FL33380
C                                                                        FL33390
C     **   KMOL(K) IS A POINTER USED TO REORDER THE AMOUNTS WHEN PRINTIN FL33400
C                                                                        FL33410
      DATA KMOL/1,2,3,11,8,5,9,10,4,6,7,12,13,14,16,15/                  FL33420
C                                                                        FL33430
C     **   INITIALIZE CONSTANTS AND CLEAR CUMULATIVE VARIABLES           FL33440
C     **   DELTAS IS THE NOMINAL PATH LENGTH INCRENMENT USED IN THE RAY  FL33450
C                                                                        FL33460
      H1 = H1F                                                           FL33470
      H2 = H2F                                                           FL33480
      ANGLE = ANGLEF                                                     FL33490
      HMIN = HMINF                                                       FL33500
      LEN = LENF                                                         FL33510
      PHI = PHIF                                                         FL33520
      IERROR = IERRF                                                     FL33530
      DELTAS = 5.0                                                       FL33540
      JMAXST = 1                                                         FL33550
      IERROR = 0                                                         FL33560
      RE = REE                                                           FL33570
      IMOD = ML                                                          FL33580
      IMAX = ML                                                          FL33590
C                                                                        FL33600
C     **   ZERO OUT CUMULATIVE VARIABLES                                 FL33610
C                                                                        FL33620
      DO 10 I = 1, 68                                                    FL33630
         LJ(I) = 0                                                       FL33640
         SL(I) = 0.0                                                     FL33650
         PPSUML(I) = 0.0                                                 FL33660
         TPSUML(I) = 0.0                                                 FL33670
         RHOSML(I) = 0.0                                                 FL33680
         DO 8 K = 1, KMAX   
            AMTL(K,I) = 0.0                                              FL33700
 8       CONTINUE     
 10   CONTINUE                                                           FL33710
      ZMAX = ZMDL(IMAX)                                                  FL33720
      IF (ITYPE.GE.2) GO TO 60                                           FL33730
C                                                                        FL33740
C     **   HORIZONTAL PATH, MODEL EQ 1 TO 7:  INTERPOLATE PROFILE TO H1  FL33750
C                                                                        FL33760
      ZL(1) = H1                                                         FL33770
      IF (ML.EQ.1) THEN                                                  FL33780
         TL(1) = TM(1)                                                   FL33790
         LJ(1) = 1                                                       FL33800
         SL(1) = RANGE                                                   FL33810
      ELSE                                                               FL33820
         DO 20 I = 2, ML                                                 FL33830
            I2 = I                                                       FL33840
            IF (H1.LT.ZMDL(I)) GO TO 30                                  FL33850
   20    CONTINUE                                                        FL33860
   30    CONTINUE                                                        FL33870
         I1 = I2-1                                                       FL33880
         FAC = (H1-ZMDL(I1))/(ZMDL(I2)-ZMDL(I1))                         FL33890
         CALL EXPINT (PL(1),PM(I1),PM(I2),FAC)                           FL33900
         TL(1) = TM(I1)+(TM(I2)-TM(I1))*FAC                              FL33910
         II1 = I1                                                        FL33920
         IF (FAC.GT.0.5) II1 = I2                                        FL33930
         LJ(1) = II1                                                     FL33940
         SL(II1) = RANGE                                                 FL33950
         DO 40 K = 1, KMAX                                               FL33960
            CALL EXPINT (DENL(K,1),DENSTY(K,I1),DENSTY(K,I2),FAC)        FL33970
   40    CONTINUE                                                        FL33980
      ENDIF                                                              FL33990
C                                                                        FL34000
C     **   CALCULATE ABSORBER AMOUNTS FOR A HORIZONTAL PATH              FL34010
C                                                                        FL34020
      WRITE (IPR,900) H1,RANGE,MODEL                                     FL34030
      TBBY(1) = TL(1)                                                    FL34040
      IKMAX = 1                                                          FL34050
      DO 50 K = 1, KMAX                                                  FL34060
         IF (ML.EQ.1) DENL(K,1) = DENSTY(K,1)                            FL34070
         W(K) = DENL(K,1)*RANGE                                          FL34080
         WPATH(1,K) = W(K)                                               FL34090
   50 CONTINUE                                                           FL34100
      WTEM = (296.0-TL(1))/(296.0-260.0)                                 FL34110
      WTEM = MAX(WTEM,0.)                                                FL34120
      WTEM = MIN(WTEM,1.)                                                FL34130
      W(9) = W(5)*WTEM                                                   FL34140
      WPATH(1,9) = W(9)                                                  FL34150
      GO TO 170                                                          FL34160
   60 CONTINUE                                                           FL34170
C                                                                        FL34180
C     **   SLANT PATH SELECTED                                           FL34190
C     **   INTERPRET SLANT PATH PARAMETERS                               FL34200
C                                                                        FL34210
      IF (IERROR.EQ.0) GO TO 70                                          FL34220
      IF (ISSGEO.NE.1) WRITE (IPR,905)                                   FL34230
      RETURN                                                             FL34240
   70 CONTINUE                                                           FL34250
C                                                                        FL34260
C     **   CALCULATE THE PATH THROUGH THE ATMOSPHERE                     FL34270
C                                                                        FL34280
      IAMT = 1                                                           FL34290
      CALL RFPATL (H1,H2,ANGLE,PHI,LEN,HMIN,IAMT,RANGE,BETA,BENDNG)      FL34300
C                                                                        FL34310
C     **   UNFOLD LAYER AMOUNTS IN AMTP INTO THE CUMULATIVE              FL34320
C     **   AMOUNTS IN WPATH FROM H1 TO H2                                FL34330
C                                                                        FL34340
      DO 80 I = 1, IPATH                                                 FL34350
         IF (H1.EQ.ZL(I)) IH1 = I                                        FL34360
         IF (H2.EQ.ZL(I)) IH2 = I                                        FL34370
   80 CONTINUE                                                           FL34380
      JMAX = (IPATH-1)+LEN*(MIN0(IH1,IH2)-1)                             FL34390
      IKMAX = JMAX                                                       FL34400
C                                                                        FL34410
C     **   DETERMINE LJ(J), WHICH IS THE NUMBER OF THE LAYER IN AMTP(K,L FL34420
C     **   STARTING FROM HMIN, WHICH CORRESPONDS TO THE LAYER J IN       FL34430
C     **   WPATH(J,K), STARTING FROM H1                                  FL34440
C     **   INITIAL DIRECTION OF PATH IS DOWN                             FL34450
C                                                                        FL34460
      L = IH1                                                            FL34470
      LDEL = -1                                                          FL34480
      IF (LEN.EQ.1.OR.H1.GT.H2) GO TO 90                                 FL34490
C                                                                        FL34500
C     **   INITIAL DIRECTION OF PATH IS UP                               FL34510
C                                                                        FL34520
      L = 0                                                              FL34530
      LDEL = 1                                                           FL34540
   90 CONTINUE                                                           FL34550
      JTURN = 0                                                          FL34560
      JMAXP1 = JMAX+1                                                    FL34570
      DO 110 J = 1, JMAXP1                                               FL34580
C                                                                        FL34590
C        **   TEST FOR REVERSING DIRECTION OF PATH FROM DOWN TO UP       FL34600
C                                                                        FL34610
         IF (L.NE.1.OR.LDEL.NE.-1) GO TO 100                             FL34620
         JTURN = J                                                       FL34630
         L = 0                                                           FL34640
         LDEL = +1                                                       FL34650
  100    CONTINUE                                                        FL34660
         L = L+LDEL                                                      FL34670
         LJ(J) = L                                                       FL34680
  110 CONTINUE                                                           FL34690
C                                                                        FL34700
C     **   LOAD TBBY AND WPATH                                           FL34710
C     **   TBBY IS DENSITY WEIGHTED MEAN TEMPERATURE                     FL34720
C                                                                        FL34730
      AMTTOT = 0.                                                        FL34740
      DO 120 K = 1, KMAX                                                 FL34750
         W(K) = 0.0                                                      FL34760
         WPATH(1,K) = 0.0                                                FL34770
  120 CONTINUE                                                           FL34780
      IMAX = 0                                                           FL34790
      DO 140 J = 1, JMAX                                                 FL34800
         L = LJ(J)                                                       FL34810
         IMAX = MAX(IMAX,L)                                              FL34820
         TBBY(L) = TPSUML(L)/RHOSML(L)                                   FL34830
         AMTTOT = AMTTOT+RHOSML(L)                                       FL34840
         DO 130 K = 1, KMAX                                              FL34850
            IF (K.EQ.9) GO TO 130                                        FL34860
C                                                                        FL34870
C           CC                                                           FL34880
C                                                                        FL34890
            WPATH(L,K) = AMTL(K,L)                                       FL34900
            W(K) = W(K)+WPATH(L,K)                                       FL34910
  130    CONTINUE                                                        FL34920
         WTEM = (296.0-TBBY(L))/(296.0-260.0)                            FL34930
         IF (WTEM.LT.0.0) WTEM = 0.                                      FL34940
         IF (WTEM.GT.1.0) WTEM = 1.0                                     FL34950
         WPATH(L,9) = WTEM*AMTL(5,L)                                     FL34960
         W(9) = W(9)+WPATH(L,9)                                          FL34970
  140 CONTINUE                                                           FL34980
      JMAX = IMAX                                                        FL34990
      JMAXST = IMAX                                                      FL35000
      JMAX = IMAX                                                        FL35010
      IKMAX = IMAX                                                       FL35020
C                                                                        FL35030
C     **   INCLUDE BOUNDARY EMISSION IF:                                 FL35040
C     **       1. TBOUND IS SET TO ZERO IN THIS VERSION OF LOWTRAN       FL35050
C     **       2. SLANT PATH INTERSECTS THE EARTH (TBOUND                FL35060
C     **          SET TO TEMPERATURE OF LOWEST BOUNDARY)                 FL35070
C                                                                        FL35080
      IF (TBOUND.EQ.0.0.AND.H2.EQ.ZMDL(1)) TBOUND = TM(1)                FL35090
C                                                                        FL35100
C     **   PRINT CUMULATIVE ABSORBER AMOUNTS                             FL35110
C                                                                        FL35120
      IF (NPR.EQ.1) GO TO 160                                            FL35130
      WRITE (IPR,910)                                                    FL35140
      DO 150 J = 1, JMAX                                                 FL35150
         LZ = J+1                                                        FL35160
         L1 = LZ-1                                                       FL35170
         IF (NPR.NE.1) WRITE (IPR,915) J,ZL(L1),ZL(LZ),TBBY(J),WPATH(J,  FL35180
     *      KMOL(1)),(WPATH(J,KMOL(K)),K=10,15)                          FL35190
  150 CONTINUE                                                           FL35200
C                                                                        FL35210
C     **   PRINT PATH SUMMARY                                            FL35220
C                                                                        FL35230
  160 WRITE (IPR,920) H1,H2,ANGLE,RANGE,BETA,PHI,HMIN,BENDNG,LEN         FL35240
  170 CONTINUE                                                           FL35250
C                                                                        FL35260
C     **   CALCULATE THE AEROSOL WEIGHTED MEAN RH                        FL35270
C                                                                        FL35280
      IF (W(7).GT.0.0.AND.ICH(1).LE.7) THEN                              FL35290
         W15 = W(15)/W(7)                                                FL35300
C                                                                        FL35310
C        INVERSE OF LOG REL HUM                                          FL35320
C                                                                        FL35330
         W(15) = 100.-EXP(W15)                                           FL35340
         GO TO 180                                                       FL35350
      ENDIF                                                              FL35360
      IF (W(12).GT.0.0.AND.ICH(1).GT.7) THEN                             FL35370
         W15 = W(15)/W(12)                                               FL35380
C                                                                        FL35390
C        INVERSE OF LOG REL HUM                                          FL35400
C                                                                        FL35410
         W(15) = 100.-EXP(W15)                                           FL35420
         GO TO 180                                                       FL35430
      ENDIF                                                              FL35440
      W(15) = 0.                                                         FL35450
  180 CONTINUE                                                           FL35460
C                                                                        FL35470
C     **   PRINT TOTAL PATH AMOUNTS                                      FL35480
C                                                                        FL35490
      WRITE (IPR,925) (W(KMOL(K)),K=10,16)                               FL35500
C                                                                        FL35510
      IF (JMAXST.GT.MAXGEO) THEN                                         FL35520
         WRITE (IPR,930) MAXGEO,JMAXST                                   FL35530
         STOP 'GEO: JMAXST .GT. MAXGEO'                                  FL35540
      ENDIF                                                              FL35550
      DO 190 IK = 1, JMAXST                                              FL35560
         IL = LJ(IK)                                                     FL35570
         RNPATH(IK) = SL(IL)                                             FL35580
         RRAMTK(IK) = RRAMT(IL)                                          FL35590
  190 CONTINUE                                                           FL35600
C                                                                        FL35610
      RETURN                                                             FL35620
C                                                                        FL35630
  900 FORMAT('0HORIZONTAL PATH AT ALTITUDE = ',F10.3,                    FL35640
     *   ' KM WITH RANGE = ',F10.3,' KM, MODEL = ',I3)                   FL35650
  905 FORMAT('0GEO:  IERROR NE 0: END THIS CALCULATION AND SKIP TO'      FL35660
     *    ,' THE NEXT CASE')                                             FL35670
  910 FORMAT(////,'    LAYER   ABSORBER AMOUNTS FOR THE PATH FROM',      FL35680
     *    ' Z(J) TO Z(J+1)',//,T3,'J',T9,'Z(J)',T18,'Z(J+1)',T27,'TBAR', FL35690
     * T37,'H2O',                                                        FL35700
     * T46,'MOL SCAT',T61,'AER 1',T73,'AER 2',T85,'AER 3',T97,'AER 4',   FL35710
     * T109,'CIRRUS',/,                                                  FL35720
     * T8,'(KM)',T17,'(KM)',T28,'(K)',T32,'LOWTRN U.')                   FL35730
  915 FORMAT(I3,2F9.3,F9.2,1P8E12.3)                                     FL35740
  920 FORMAT(//,'0SUMMARY OF THE GEOMETRY CALCULATION',//,               FL35750
     * 10X,'H1      = ',F10.3,' KM',/,10X,'H2      = ',F10.3,' KM',/,    FL35760
     *10X,'ANGLE   = ',F10.3,' DEG',/,10X,'RANGE   = ',F10.3,' KM',/,    FL35770
     *10X,'BETA    = ',F10.3,' DEG',/,10X,'PHI     = ',F10.3,' DEG',/,   FL35780
     * 10X,'HMIN    = ',F10.3,' KM',/,10X,'BENDING = ',F10.3,' DEG',/,   FL35790
     * 10X,'LEN     = ',I10)                                             FL35800
  925 FORMAT(////,' EQUIVALENT SEA LEVEL TOTAL ABSORBER AMOUNTS',//,     FL35810
     *    T15,'       ',T26,'MOL SCAT',T41,'AER 1', T53,'AER 2',         FL35820
     *    T65,'AER 3',T77, 'AER 4',T87,'CIRRUS',T99,'MEAN RH'/,          FL35830
     *    T99,'(PRCNT)',//,22X,1P6E12.3,0PF12.2,//)                      FL35840
  930 FORMAT(//'  CURRENT GEOMETRY DIMENSION ',I5 ,/                     FL35850
     *,' JMAXST = ',I5,' RESET AVTRAT TDIFF1 TDIFF2 TO 2. 10. 20.')      FL35860
C                                                                        FL35870
      END                                                                FL35880
C
C     ******************************************************************
C
      SUBROUTINE FINDSL(H,SH,GAMMA)                                      FL35890
C                                                                        FL35900
C     **   GIVEN AN ALTITUDE H, THIS SUBROUTINE FINDS THE LAYER BOUNDARI FL35910
C     **   ZM(I1) AND ZM(I2) WHICH CONTAIN H,  THEN CALCULATES THE SCALE FL35920
C     **   HEIGHT (SH) AND THE VALUE AT THE GROUND (GAMMA+1) FOR THE     FL35930
C     **   INDEX OF REFRACTION                                           FL35940
C                                                                        FL35950
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
      COMMON /PARMLT/ RE,DELTAS,ZMAX,IMAX,IMOD,IBMAX,IPATH               FL35960
      COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                      FL35970
      COMMON /MODEL/ ZMDL(MXZMD),P(MXZMD),T(MXZMD),                      FL35980
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL35990
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)
C
      DO 10 IM = 2, IMOD                                                 FL36000
         I2 = IM                                                         FL36010
         IF (ZMDL(IM).GE.H) GO TO 20                                     FL36020
   10 CONTINUE                                                           FL36030
      I2 = IMOD                                                          FL36040
   20 CONTINUE                                                           FL36050
      I1 = I2-1                                                          FL36060
      CALL SCALHT (ZMDL(I1),ZMDL(I2),RFNDX(I1),RFNDX(I2),SH,GAMMA)       FL36070
      RETURN                                                             FL36080
      END                                                                FL36090
      SUBROUTINE RFPATL(H1,H2,ANGLE,PHI,LEN,HMIN,IAMT,RANGE,BETA,BENDNG) FL36100
C                                                                        FL36110
C     ****************************************************************** FL36120
C     THIS SUBROUTINE TRACES THE REFRACTED RAY FROM H1 WITH A            FL36130
C     INITIAL ZENITH ANGLE ANGLE TO H2 WHERE THE ZENITH ANGLE IS PHI,    FL36140
C     AND CALCULATES THE ABSORBER AMOUNTS (IF IAMT.EQ.1) ALONG           FL36150
C     THE PATH.  IT STARTS FROM THE LOWEST POINT ALONG THE PATH          FL36160
C     (THE TANGENT HEIGHT HMIN IF LEN = 1 OR HA = MIN(H1,H2) IF LEN = 0) FL36170
C     AND PROCEEDS TO THE HIGHEST POINT.  BETA AND RANGE ARE THE         FL36180
C     EARTH CENTERED ANGLE AND THE TOTAL DISTANCE RESPECTIVELY           FL36190
C     FOR THE REFRACTED PATH FROM H1 TO H2                               FL36200
C     ****************************************************************** FL36210
C                                                                        FL36220
      PARAMETER (MXZMD=6000, MXZ20 = MXZMD+20, MX2Z3 = 2*MXZMD+3)
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL36230
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL36240
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL36250
      COMMON /PARMLT/ RE,DELTAS,ZMAX,IMAX,IMOD,IBMAX,IPATH               FL36260
      COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                      FL36270
C
C     RFRPTH is dependent upon MXZMD (MXZ20=MXZMD+20;MX2Z3=2*MXZMD+3)
C
      COMMON  /RFRPTH/ ZL(MXZ20),PL(MXZ20),TL(MXZ20),RFNDXL(MXZ20),      FL36280
     *     SL(MXZ20),PPSUML(MXZ20),TPSUML(MXZ20),RHOSML(MXZ20),          FL36290
     *     DENL(16,MXZ20),AMTL(16,MXZ20),LJ(MX2Z3)                       FL36300
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,N_LVL,JH1              FL36310
C
      DATA I_2/2/
C
      IF (H1.GT.H2) GO TO 10                                             FL36340
      IORDER = 1                                                         FL36350
      HA = H1                                                            FL36360
      HB = H2                                                            FL36370
      ANGLEA = ANGLE                                                     FL36380
      GO TO 20                                                           FL36390
   10 CONTINUE                                                           FL36400
      IORDER = -1                                                        FL36410
      HA = H2                                                            FL36420
      HB = H1                                                            FL36430
      ANGLEA = PHI                                                       FL36440
   20 CONTINUE                                                           FL36450
      JNEXT = 1                                                          FL36460
C                                                                        FL36470
C     IF(IAMT.EQ.1 .AND. NPR.NE.1)  WRITE(IPR,20)                        FL36480
C                                                                        FL36490
      IF (LEN.EQ.0) GO TO 30                                             FL36500
C                                                                        FL36510
C     **   LONG PATH: FILL IN THE SYMETRIC PART FROM THE TANGENT HEIGHT  FL36520
C     **   TO HA                                                         FL36530
C                                                                        FL36540
      CALL FILL (HMIN,HA,JNEXT)                                          FL36550
      JHA = JNEXT                                                        FL36560
      JH1 = JNEXT-1                                                      FL36570
   30 CONTINUE                                                           FL36580
C                                                                        FL36590
C     **   FILL IN THE REMAINING PATH FROM HA TO HB                      FL36600
C                                                                        FL36610
      IF (HA.EQ.HB) GO TO 40                                             FL36620
      CALL FILL (HA,HB,JNEXT)                                            FL36630
   40 CONTINUE                                                           FL36640
      JMAX = JNEXT                                                       FL36650
      IPATH = JMAX                                                       FL36660
C                                                                        FL36670
C     **   INTEGRATE EACH SEGMENT OF THE PATH                            FL36680
C     **   CALCULATE CPATH SEPERATELY FOR LEN = 0,1                      FL36690
C                                                                        FL36700
      IF (LEN.EQ.1) GO TO 50                                             FL36710
      CALL FINDSL (HA,SH,GAMMA)                                          FL36720
      CPATH = (RE+HA)*ANDEX(HA,SH,GAMMA)*SIN(ANGLEA/DEG)                 FL36730
      GO TO 60                                                           FL36740
   50 CONTINUE                                                           FL36750
      CALL FINDSL (HMIN,SH,GAMMA)                                        FL36760
      CPATH = (RE+HMIN)*ANDEX(HMIN,SH,GAMMA)                             FL36770
   60 CONTINUE                                                           FL36780
      BETA = 0.0                                                         FL36790
      S = 0.0                                                            FL36800
      BENDNG = 0.0                                                       FL36810
      IF (LEN.EQ.0) GO TO 100                                            FL36820
C                                                                        FL36830
C     **   DO SYMETRIC PART, FROM TANGENT HEIGHT(HMIN) TO HA             FL36840
C                                                                        FL36850
      IHLOW = 1                                                          FL36860
      IF (IORDER.EQ.-1) IHLOW = 2                                        FL36870
C                                                                        FL36880
      SINAI = 1.0                                                        FL36910
      COSAI = 0.0                                                        FL36920
      THETA = 90.0                                                       FL36930
      J2 = JHA-1                                                         FL36940
      DO 90 J = 1, J2                                                    FL36950
         CALL SCALHT (ZL(J),ZL(J+1),RFNDXL(J),RFNDXL(J+1),SH,GAMMA)      FL36960
         CALL LOLAYR (J,SINAI,COSAI,CPATH,SH,GAMMA,IAMT,DS,DBEND)        FL36970
         DBEND = DBEND*DEG                                               FL36980
         PHI = ASIN(SINAI)*DEG                                           FL36990
         DBETA = THETA-PHI+DBEND                                         FL37000
         PHI = 180.0-PHI                                                 FL37010
         S = S+DS                                                        FL37020
         BENDNG = BENDNG+DBEND                                           FL37030
         BETA = BETA+DBETA                                               FL37040
         IF (IAMT.NE.1) GO TO 70                                         FL37050
         PBAR = PPSUML(J)/RHOSML(J)                                      FL37060
         TBAR = TPSUML(J)/RHOSML(J)                                      FL37070
         RHOBAR = RHOSML(J)/DS                                           FL37080
C                                                                        FL37090
C        IF(IAMT.EQ.1 .AND. NPR.NE.1) WRITE(IPR,22) J,ZP(J),ZP(J+1),     FL37100
C        1    THETA,DS,S,DBETA,BETA,PHI,DBEND,BENDNG,PBAR,TBAR,RHOBAR    FL37110
C                                                                        FL37120
   70    CONTINUE                                                        FL37130
         IF (ISSGEO.EQ.1) GO TO 80                                       FL37140
C                                                                        FL37150
C        CC   ATHETA(J)=THETA                                            FL37160
C        CC   ADBETA(J)=DBETA                                            FL37170
C                                                                        FL37180
   80    CONTINUE                                                        FL37190
         THETA = 180.0-PHI                                               FL37200
   90 CONTINUE                                                           FL37210
C                                                                        FL37220
C     **   DOUBLE PATH QUANTITIES FOR THE OTHER PART OF THE SYMETRIC PAT FL37230
C                                                                        FL37240
      BENDNG = 2.0*BENDNG                                                FL37250
      BETA = 2.0*BETA                                                    FL37260
      S = 2.0*S                                                          FL37270
C                                                                        FL37280
C     IF(IAMT.EQ.1 .AND. NPR.NE.1) WRITE(IPR,26) S,BETA,BENDNG           FL37290
C                                                                        FL37300
      JNEXT = JHA                                                        FL37310
      GO TO 120                                                          FL37320
  100 CONTINUE                                                           FL37330
C                                                                        FL37340
C     **   SHORT PATH                                                    FL37350
C                                                                        FL37360
      JNEXT = 1                                                          FL37370
C                                                                        FL37380
C     **   ANGLEA IS THE ZENITH ANGLE AT HA IN DEG                       FL37390
C     **   SINAI IS SIN OF THE INCIDENCE ANGLE                           FL37400
C     **   COSAI IS CARRIED SEPERATELY TO AVOID A PRECISION PROBLEM      FL37410
C     **   WHEN SINAI IS CLOSE TO 1.0                                    FL37420
C                                                                        FL37430
      THETA = ANGLEA                                                     FL37440
      IF (ANGLEA.GT.45.0) GO TO 110                                      FL37450
      SINAI = SIN(ANGLEA/DEG)                                            FL37460
      COSAI = -COS(ANGLEA/DEG)                                           FL37470
      GO TO 120                                                          FL37480
  110 CONTINUE                                                           FL37490
      SINAI = COS((90.0-ANGLEA)/DEG)                                     FL37500
      COSAI = -SIN((90.0-ANGLEA)/DEG)                                    FL37510
  120 CONTINUE                                                           FL37520
C                                                                        FL37530
C     **   DO PATH FROM HA TO HB                                         FL37540
C                                                                        FL37550
      IF (HA.EQ.HB) GO TO 160                                            FL37560
      J1 = JNEXT                                                         FL37570
      J2 = JMAX-1                                                        FL37580
      IHLOW = 1                                                          FL37590
      IF (IORDER.EQ.-1) IHLOW = 2                                        FL37600
      IHIGH = MOD(IHLOW,I_2)+1                                             FL37610
C                                                                        FL37620
      DO 150 J = J1, J2                                                  FL37660
         CALL SCALHT (ZL(J),ZL(J+1),RFNDXL(J),RFNDXL(J+1),SH,GAMMA)      FL37670
         CALL LOLAYR (J,SINAI,COSAI,CPATH,SH,GAMMA,IAMT,DS,DBEND)        FL37680
         DBEND = DBEND*DEG                                               FL37690
         PHI = ASIN(SINAI)*DEG                                           FL37700
         DBETA = THETA-PHI+DBEND                                         FL37710
         PHI = 180.0-PHI                                                 FL37720
         S = S+DS                                                        FL37730
         BENDNG = BENDNG+DBEND                                           FL37740
         BETA = BETA+DBETA                                               FL37750
         IF (IAMT.NE.1) GO TO 130                                        FL37760
         PBAR = PPSUML(J)/RHOSML(J)                                      FL37770
         TBAR = TPSUML(J)/RHOSML(J)                                      FL37780
         RHOBAR = RHOSML(J)/DS                                           FL37790
C                                                                        FL37800
C        IF(IAMT.EQ.1 .AND. NPR.NE.1) WRITE(IPR,22) J,ZP(J),ZP(J+1),     FL37810
C        1    THETA,DS,S,DBETA,BETA,PHI,DBEND,BENDNG,PBAR,TBAR,RHOBAR    FL37820
C                                                                        FL37830
  130    CONTINUE                                                        FL37840
         IF (ISSGEO.EQ.1) GO TO 140                                      FL37850
C                                                                        FL37860
C        CC   ADBETA(J)=DBETA                                            FL37870
C        CC   ATHETA(J)=THETA                                            FL37880
C                                                                        FL37890
  140    CONTINUE                                                        FL37900
         THETA = 180.0-PHI                                               FL37910
  150 CONTINUE                                                           FL37920
  160 CONTINUE                                                           FL37930
C                                                                        FL37940
C     CC   IF(ISSGEO.EQ.0) ATHETA(JMAX)=THETA                            FL37950
C                                                                        FL37960
      IF (IORDER.EQ.-1) PHI = ANGLEA                                     FL37970
      RANGE = S                                                          FL37980
      RETURN                                                             FL37990
C                                                                        FL38000
C                                                                        FL38010
      END                                                                FL38020
C
C     ******************************************************************
C
      SUBROUTINE FILL(HA,HB,JNEXT)                                       FL38030
C                                                                        FL38040
C     ****************************************************************** FL38050
C     THIS SUBROUTINE DEFINES THE ATMOSPHERIC BOUNDARIES OF THE PATH     FL38060
C     FROM HA TO HB AND INTERPOLATES (EXTRAPOLATES) THE DENSITIES TO     FL38070
C     THESE BOUNDARIES ASSUMING THE DENSITIES VARY EXPONENTIALLY         FL38080
C     WITH HEIGHT                                                        FL38090
C     ****************************************************************** FL38100
C                                                                        FL38110
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
      PARAMETER (MXZ20 = MXZMD+20, MX2Z3 = 2*MXZMD+3)
C
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL38120
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL38130
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL38140
      COMMON /MODEL/ ZMDL(MXZMD),P(MXZMD),T(MXZMD),                      FL38150
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL38160
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)
      COMMON /PARMLT/ RE,DELTAS,ZMAX,IMAX,IMOD,IBMAX,IPATH               FL38170
      COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                      FL38180
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,N_LVL,JH1              FL38190
C
C     RFRPTH is dependent upon MXZMD (MXZ20=MXZMD+20;MX2Z3=2*MXZMD+3)
C
      COMMON  /RFRPTH/ ZL(MXZ20),PL(MXZ20),TL(MXZ20),RFNDXL(MXZ20),      FL38200
     *     SL(MXZ20),PPSUML(MXZ20),TPSUML(MXZ20),RHOSML(MXZ20),          FL38210
     *     DENL(16,MXZ20),AMTL(16,MXZ20),LJ(MX2Z3)                       FL38220
C
      IF (HA.LT.HB) GO TO 10                                             FL38230
      WRITE (IPR,900) HA,HB,JNEXT                                        FL38240
      STOP                                                               FL38250
   10 CONTINUE                                                           FL38260
C                                                                        FL38270
C     **   FIND ZMDL(IA): THE SMALLEST ZMDL(I).GT.HA                     FL38280
C                                                                        FL38290
      DO 20 I = 1, IMAX                                                  FL38300
         IF (HA.GE.ZMDL(I)) GO TO 20                                     FL38310
         IA = I                                                          FL38320
         GO TO 30                                                        FL38330
   20 CONTINUE                                                           FL38340
      IA = IMAX+1                                                        FL38350
      IB = IA                                                            FL38360
      GO TO 50                                                           FL38370
C                                                                        FL38380
C     **   FIND ZMDL(IB): THE SMALLEST ZMDL(I).GE.HB                     FL38390
C                                                                        FL38400
   30 CONTINUE                                                           FL38410
      DO 40 I = IA, IMAX                                                 FL38420
         IF (HB-ZMDL(I).GT..0001) GO TO 40                               FL38430
         IB = I                                                          FL38440
         GO TO 50                                                        FL38450
   40 CONTINUE                                                           FL38460
      IB = IMAX+1                                                        FL38470
   50 CONTINUE                                                           FL38480
C                                                                        FL38490
C     **   INTERPOLATE DENSITIES TO HA, HB                               FL38500
C                                                                        FL38510
      ZL(JNEXT) = HA                                                     FL38520
      I2 = IA                                                            FL38530
      IF (I2.EQ.1) I2 = 2                                                FL38540
      I2 = MIN(I2,IMAX)                                                  FL38550
      I1 = I2-1                                                          FL38560
      A = (HA-ZMDL(I1))/(ZMDL(I2)-ZMDL(I1))                              FL38570
      CALL EXPINT (PL(JNEXT),P(I1),P(I2),A)                              FL38580
      TL(JNEXT) = T(I1)+(T(I2)-T(I1))*A                                  FL38590
      CALL EXPINT (RFNDXL(JNEXT),RFNDX(I1),RFNDX(I2),A)                  FL38600
      DO 60 K = 1, KMAX                                                  FL38610
         CALL EXPINT (DENL(K,JNEXT),DENSTY(K,I1),DENSTY(K,I2),A)         FL38620
   60 CONTINUE                                                           FL38630
      IF (IA.EQ.IB) GO TO 80                                             FL38640
C                                                                        FL38650
C     **   FILL IN DENSITIES BETWEEN HA AND HB                           FL38660
C                                                                        FL38670
      I1 = IA                                                            FL38680
      I2 = IB-1                                                          FL38690
      DO 70 I = I1, I2                                                   FL38700
         JNEXT = JNEXT+1                                                 FL38710
         ZL(JNEXT) = ZMDL(I)                                             FL38720
         PL(JNEXT) = P(I)                                                FL38730
         TL(JNEXT) = T(I)                                                FL38740
         RFNDXL(JNEXT) = RFNDX(I)                                        FL38750
         DO 68 K = 1, KMAX  
            DENL(K,JNEXT) = DENSTY(K,I)                                  FL38770
 68      CONTINUE
 70   CONTINUE                                                           FL38780
 80   CONTINUE                                                           FL38790
C                                                                        FL38800
C     **   INTERPOLATE THE DENSITIES TO HB                               FL38810
C                                                                        FL38820
      JNEXT = JNEXT+1                                                    FL38830
      ZL(JNEXT) = HB                                                     FL38840
      I2 = IB                                                            FL38850
      IF (I2.EQ.1) I2 = 2                                                FL38860
      I2 = MIN(I2,IMAX)                                                  FL38870
      I1 = I2-1                                                          FL38880
      A = (HB-ZMDL(I1))/(ZMDL(I2)-ZMDL(I1))                              FL38890
      CALL EXPINT (PL(JNEXT),P(I1),P(I2),A)                              FL38900
      TL(JNEXT) = T(I1)+(T(I2)-T(I1))*A                                  FL38910
      CALL EXPINT (RFNDXL(JNEXT),RFNDX(I1),RFNDX(I2),A)                  FL38920
      DO 90 K = 1, KMAX                                                  FL38930
         CALL EXPINT (DENL(K,JNEXT),DENSTY(K,I1),DENSTY(K,I2),A)         FL38940
   90 CONTINUE                                                           FL38950
      RETURN                                                             FL38960
C                                                                        FL38970
  900 FORMAT('0SUBROUTINE FILL- ERROR, HA .GE. HB',//,                   FL38980
     *    10X,'HA, HB, JNEXT = ',2E25.15,I6)                             FL38990
C                                                                        FL39000
      END                                                                FL39010
C
C     *****************************************************************
C
      SUBROUTINE LOLAYR(J,SINAI,COSAI,CPATH,SH,GAMMA,IAMT,S,BEND)        FL39020
C                                                                        FL39030
C     *****************************************************************  FL39040
C     THIS SUBROUTINE CALCULATES THE REFRACTED PATH FROM Z1 TO Z2        FL39050
C     WITH THE SIN OF THE INITIAL INCIDENCE ANGLE SINAI                  FL39060
C     *****************************************************************  FL39070
C                                                                        FL39080
      PARAMETER (MXZMD=6000, MXZ20 = MXZMD+20, MX2Z3 = 2*MXZMD+3)
C
      COMMON /PARMLT/ RE,DELTAS,ZMAX,IMAX,IMOD,IBMAX,IPATH               FL39090
      COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                      FL39100
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,N_LVL,JH1              FL39110
C
C     RFRPTH is dependent upon MXZMD (MXZ20=MXZMD+20;MX2Z3=2*MXZMD+3)
C
      COMMON  /RFRPTH/ ZL(MXZ20),PL(MXZ20),TL(MXZ20),RFNDXL(MXZ20),      FL39120
     *     SL(MXZ20),PPSUML(MXZ20),TPSUML(MXZ20),RHOSML(MXZ20),          FL39130
     *     DENL(16,MXZ20),AMTL(16,MXZ20),LJ(MX2Z3)                       FL39140
C
      DIMENSION HDEN(20),DENA(20),DENB(20)                               FL39150
C
      DATA EPSILN/1.0E-5/                                                FL39160
C                                                                        FL39170
C     **   INITIALIZE LOOP                                               FL39180
C                                                                        FL39190
      N = 0                                                              FL39200
      Z1 = ZL(J)                                                         FL39210
      Z2 = ZL(J+1)                                                       FL39220
      H1 = Z1                                                            FL39230
      R1 = RE+H1                                                         FL39240
      DHMIN = DELTAS**2/(2.0*R1)                                         FL39250
      SINAI1 = SINAI                                                     FL39260
      COSAI1 = COSAI                                                     FL39270
      Y1 = COSAI1**2/2.0+COSAI1**4/8.0+COSAI1**6*3.0/48.0                FL39280
      Y3 = 0.0                                                           FL39290
      X1 = -R1*COSAI1                                                    FL39300
      RATIO1 = R1/RADREF(H1,SH,GAMMA)                                    FL39310
      DSDX1 = 1.0/(1.0-RATIO1*SINAI1**2)                                 FL39320
      DBNDX1 = DSDX1*SINAI1*RATIO1/R1                                    FL39330
      S = 0.0                                                            FL39340
      BEND = 0.0                                                         FL39350
      IF (IAMT.EQ.2) GO TO 50                                            FL39360
C                                                                        FL39370
C     **   INITIALIZE THE VARIABLES FOR THE CALCULATION OF THE           FL39380
C     **   ABSORBER AMOUNTS                                              FL39390
C                                                                        FL39400
      PA = PL(J)                                                         FL39410
      PB = PL(J+1)                                                       FL39420
      TA = TL(J)                                                         FL39430
      TB = TL(J+1)                                                       FL39440
      RHOA = PA/(GCAIR*TA)                                               FL39450
      RHOB = PB/(GCAIR*TB)                                               FL39460
      DZ = ZL(J+1)-ZL(J)                                                 FL39470
      HP = -DZ/ LOG(PB/PA)                                               FL39480
      IF (ABS(RHOB/RHOA-1.0).LT.EPSILN) GO TO 10                         FL39490
      HRHO = -DZ/ LOG(RHOB/RHOA)                                         FL39500
      GO TO 20                                                           FL39510
   10 HRHO = 1.0E30                                                      FL39520
   20 CONTINUE                                                           FL39530
      DO 40 K = 1, KMAX                                                  FL39540
         DENA(K) = DENL(K,J)                                             FL39550
         DENB(K) = DENL(K,J+1)                                           FL39560
         IF (DENA(K).LE.0.0.OR.DENB(K).LE.0.0) GO TO 30                  FL39570
         IF (ABS(1.0-DENA(K)/DENB(K)).LE.EPSILN) GO TO 30                FL39580
C                                                                        FL39590
C        **   USE EXPONENTIAL INTERPOLATION                              FL39600
C                                                                        FL39610
         HDEN(K) = -DZ/ LOG(DENB(K)/DENA(K))                             FL39620
         GO TO 40                                                        FL39630
C                                                                        FL39640
C        **   USE LINEAR INTERPOLATION                                   FL39650
C                                                                        FL39660
   30    HDEN(K) = 0.0                                                   FL39670
   40 CONTINUE                                                           FL39680
   50 CONTINUE                                                           FL39690
C                                                                        FL39700
C     **   LOOP THROUGH PATH                                             FL39710
C     **   INTEGRATE PATH QUANTITIES USING QUADRATIC INTEGRATION WITH    FL39720
C     **   UNEQUALLY SPACED POINTS                                       FL39730
C                                                                        FL39740
   60 CONTINUE                                                           FL39750
      N = N+1                                                            FL39760
      DH = -DELTAS*COSAI1                                                FL39770
      DHMIN = MAX(DH,DHMIN)                                              FL39780
      H3 = H1+DH                                                         FL39790
      H3 = MIN(H3,Z2)                                                    FL39800
      DH = H3-H1                                                         FL39810
      R3 = RE+H3                                                         FL39820
      H2 = H1+DH/2.0                                                     FL39830
      R2 = RE+H2                                                         FL39840
      SINAI2 = CPATH/(ANDEX(H2,SH,GAMMA)*R2)                             FL39850
      SINAI3 = CPATH/(ANDEX(H3,SH,GAMMA)*R3)                             FL39860
      RATIO2 = R2/RADREF(H2,SH,GAMMA)                                    FL39870
      RATIO3 = R3/RADREF(H3,SH,GAMMA)                                    FL39880
      IF ((1.0-SINAI2).GT.EPSILN) GO TO 70                               FL39890
C                                                                        FL39900
C     **   NEAR A TANGENT HEIGHT, COSAI = -SQRT(1-SINAI**2) LOSES        FL39910
C     **   PRECISION. USE THE FOLLOWING ALGORITHM TO GET COSAI.          FL39920
C                                                                        FL39930
      Y3 = Y1+(SINAI1*(1.0-RATIO1)/R1+4.0*SINAI2*(1.0-RATIO2)/R2+SINAI3* FL39940
     *   (1.0-RATIO3)/R3)*DH/6.0                                         FL39950
      COSAI3 = -SQRT(2.0*Y3-Y3**2)                                       FL39960
      X3 = -R3*COSAI3                                                    FL39970
      DX = X3-X1                                                         FL39980
      W1 = 0.5*DX                                                        FL39990
      W2 = 0.0                                                           FL40000
      W3 = 0.5*DX                                                        FL40010
      GO TO 90                                                           FL40020
C                                                                        FL40030
   70 CONTINUE                                                           FL40040
      COSAI2 = -SQRT(1.0-SINAI2**2)                                      FL40050
      COSAI3 = -SQRT(1.0-SINAI3**2)                                      FL40060
      X2 = -R2*COSAI2                                                    FL40070
      X3 = -R3*COSAI3                                                    FL40080
C                                                                        FL40090
C     **   CALCULATE WEIGHTS                                             FL40100
C                                                                        FL40110
      D31 = X3-X1                                                        FL40120
      D32 = X3-X2                                                        FL40130
      D21 = X2-X1                                                        FL40140
      IF (D32.EQ.0.0.OR.D21.EQ.0.0) GO TO 80                             FL40150
      W1 = (2-D32/D21)*D31/6.0                                           FL40160
      W2 = D31**3/(D32*D21*6.0)                                          FL40170
      W3 = (2.0-D21/D32)*D31/6.0                                         FL40180
      GO TO 90                                                           FL40190
   80 CONTINUE                                                           FL40200
      W1 = 0.5*D31                                                       FL40210
      W2 = 0.0                                                           FL40220
      W3 = 0.5*D31                                                       FL40230
C                                                                        FL40240
   90 CONTINUE                                                           FL40250
      DSDX2 = 1.0/(1.0-RATIO2*SINAI2**2)                                 FL40260
      DSDX3 = 1.0/(1.0-RATIO3*SINAI3**2)                                 FL40270
      DBNDX2 = DSDX2*SINAI2*RATIO2/R2                                    FL40280
      DBNDX3 = DSDX3*SINAI3*RATIO3/R3                                    FL40290
C                                                                        FL40300
C     **   INTEGRATE                                                     FL40310
C                                                                        FL40320
      DS = W1*DSDX1+W2*DSDX2+W3*DSDX3                                    FL40330
      S = S+DS                                                           FL40340
      DBEND = W1*DBNDX1+W2*DBNDX2+W3*DBNDX3                              FL40350
      BEND = BEND+DBEND                                                  FL40360
      IF (IAMT.EQ.2) GO TO 150                                           FL40370
C                                                                        FL40380
C     **   CALCULATE AMOUNTS                                             FL40390
C                                                                        FL40400
      DSDZ = DS/DH                                                       FL40410
      PB = PA*EXP(-DH/HP)                                                FL40420
      RHOB = RHOA*EXP(-DH/HRHO)                                          FL40430
      IF ((DH/HRHO).LT.EPSILN) GO TO 100                                 FL40440
      PPSUML(J) = PPSUML(J)+DSDZ*(HP/(1.0+HP/HRHO))*(PA*RHOA-PB*RHOB)    FL40450
      TPSUML(J) = TPSUML(J)+DSDZ*HP*(PA-PB)/GCAIR                        FL40460
      RHOSML(J) = RHOSML(J)+DSDZ*HRHO*(RHOA-RHOB)                        FL40470
      GO TO 110                                                          FL40480
  100 CONTINUE                                                           FL40490
      PPSUML(J) = PPSUML(J)+0.5*DS*(PA*RHOA+PB*RHOB)                     FL40500
      TPSUML(J) = TPSUML(J)+0.5*DS*(PA+PB)/GCAIR                         FL40510
      RHOSML(J) = RHOSML(J)+0.5*DS*(RHOA+RHOB)                           FL40520
  110 CONTINUE                                                           FL40530
      DO 130 K = 1, KMAX                                                 FL40540
         IF (ABS(HDEN(K)).EQ.0.0) GO TO 120                              FL40550
         IF ((DH/HDEN(K)).LT.EPSILN) GO TO 120                           FL40560
C                                                                        FL40570
C        **   EXPONENTIAL INTERPOLATION                                  FL40580
C                                                                        FL40590
         DENB(K) = DENL(K,J)*EXP(-(H3-Z1)/HDEN(K))                       FL40600
         AMTL(K,J) = AMTL(K,J)+DSDZ*HDEN(K)*(DENA(K)-DENB(K))            FL40610
         GO TO 130                                                       FL40620
  120    CONTINUE                                                        FL40630
C                                                                        FL40640
C        **   LINEAR INTERPOLATION                                       FL40650
C                                                                        FL40660
         DENB(K) = DENL(K,J)+(DENL(K,J+1)-DENL(K,J))*(H3-Z1)/DZ          FL40670
         AMTL(K,J) = AMTL(K,J)+0.5*(DENA(K)+DENB(K))*DS                  FL40680
  130 CONTINUE                                                           FL40690
      PA = PB                                                            FL40700
      RHOA = RHOB                                                        FL40710
      DO 140 K = 1, KMAX                                                 FL40720
         DENA(K) = DENB(K)                                               FL40730
  140 CONTINUE                                                           FL40740
  150 CONTINUE                                                           FL40750
      IF (H3.GE.Z2) GO TO 160                                            FL40760
      H1 = H3                                                            FL40770
      R1 = R3                                                            FL40780
      SINAI1 = SINAI3                                                    FL40790
      RATIO1 = RATIO3                                                    FL40800
      Y1 = Y3                                                            FL40810
      COSAI1 = COSAI3                                                    FL40820
      X1 = X3                                                            FL40830
      DSDX1 = DSDX3                                                      FL40840
      DBNDX1 = DBNDX3                                                    FL40850
      GO TO 60                                                           FL40860
  160 CONTINUE                                                           FL40870
      SINAI = SINAI3                                                     FL40880
      COSAI = COSAI3                                                     FL40890
      SL(J) = S                                                          FL40900
      RETURN                                                             FL40910
      END                                                                FL40920
C
C     *****************************************************************
C
      SUBROUTINE EQULWC                                                  FL40930
C                                                                        FL40940
C     CC                                                                 FL40950
C     CC   EQUIVALENT LIQUID  WATER CONSTANTS FOR BEXT (0.55UM)=1.0KM-1  FL40960
C     CC   AWCCON(1-4) IS SET TO ONE OF THE CONSTANTS FOR EACH AEROSOL   FL40970
C     CC   IN SUBROUTINE EXABIN AND MULTIPLIED BY THE BEXT (DENSTY(N,I)) FL40980
C     CC   WHERE N=7,12,13 OR 14 AND I IS THE NUMBER OF LAYERS           FL40990
C     CC                                                                 FL41000
C                                                                        FL41010
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),VH(16),TX(16),W(16)       FL41020
      COMMON WPATH(IM2,16),TBBY(IM2)                                     FL41030
      COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),VX0(47),AWCCON(5)          FL41040
C
      CHARACTER*8      HMOD                                              FL41050
C
      COMMON /CMN/ HMOD(3),ZN(MXZMD),PN(MXZMD),TN(MXZMD),RFNDXM(MXZMD),
     *          ZP(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),
     *          TPSUM(IM2),RHOPSM(IM2),IMMAX,WGM(MXZMD),DEMW(MXZMD),
     *          AMTP(MXMOL,MXPDIM)                                        
C
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISS,N_LVL,JH1
      COMMON /MODEL/ ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),                    FL41100
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL41110
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)          FL41120
C
      DO 10 I = 1, ML                                                    FL41130
         IF (DENSTY(7,I).NE.0.0) EQLWC(I) = DENSTY(7,I)*AWCCON(1)        FL41140
         IF (DENSTY(12,I).NE.0.0) EQLWC(I) = DENSTY(12,I)*AWCCON(2)      FL41150
         IF (DENSTY(13,I).NE.0.0) EQLWC(I) = DENSTY(13,I)*AWCCON(3)      FL41160
         IF (DENSTY(14,I).NE.0.0) EQLWC(I) = DENSTY(14,I)*AWCCON(4)      FL41170
   10 CONTINUE                                                           FL41180
      RETURN                                                             FL41190
      END                                                                FL41200
C
C     *****************************************************************
C
      SUBROUTINE INDX (WAVL,TC,KEY,REIL,AIMAG)                           FL41210
C                                                                        FL41220
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  FL41230
C     * *                                                                FL41240
C     * * WAVELENGTH IS IN CENTIMETERS.  TEMPERATURE IS IN DEG. C.       FL41250
C     * *                                                                FL41260
C     * * KEY IS SET TO 1 IN SUBROUTINE GAMFOG                           FL41270
C     * *                                                                FL41280
C     * * REIL IS THE REAL PART OF THE REFRACTIVE INDEX.                 FL41290
C     * *                                                                FL41300
C     * * AIMAG IS THE IMAGINARY PART OF THE REFRACTIVE INDEX IT IS      FL41310
C     * *                                                                FL41320
C     * * RETURNED NEG. I.E.  M= REAL - I*AIMAG  .                       FL41330
C     * *                                                                FL41340
C     * * A SERIES OF CHECKS ARE MADE AND WARNINGS GIVEN.                FL41350
C     * *                                                                FL41360
C     * * RAY APPLIED OPTICS VOL 11,NO.8,AUG 72, PG. 1836-1844           FL41370
C     * *                                                                FL41380
C     * * CORRECTIONS HAVE BEEN MADE TO RAYS ORIGINAL PAPER              FL41390
C     * *                                                                FL41400
C     * *                                                                FL41410
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  FL41420
C                                                                        FL41430
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL41440
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL41450
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL41460
C
      R1 = 0.0                                                           FL41470
      R2 = 0.0                                                           FL41480
      IF (WAVL.LT..0001) WRITE (IPR,900)                                 FL41490
      IF (TC.LT.-20.) WRITE (IPR,905)                                    FL41500
      CALL DEBYE (WAVL,TC,KEY,REIL,AIMAG)                                FL41510
C                                                                        FL41520
C     * *  TABLE 3 WATER PG. 1840                                        FL41530
C                                                                        FL41540
      IF (WAVL.GT..034) GO TO 10                                         FL41550
      GO TO 20                                                           FL41560
   10 IF (WAVL.GT..1) GO TO 50                                           FL41570
      R2 = DOP(WAVL,1.83899,1639.,52340.4,10399.2,588.24,345005.,259913. FL41580
     *   ,161.29,43319.7,27661.2)                                        FL41590
      R2 = R2+R2*(TC-25.)*.0001*EXP((.000025*WAVL)**.25)                 FL41600
      REIL = REIL*(WAVL-.034)/.066+R2*(.1-WAVL)/.066                     FL41610
      GO TO 50                                                           FL41620
   20 IF (WAVL.GT..0006) GO TO 30                                        FL41630
      GO TO 40                                                           FL41640
   30 REIL = DOP(WAVL,1.83899,1639.,52340.4,10399.2,588.24,345005.,      FL41650
     *   259913.,161.29,43319.7,27661.2)                                 FL41660
      REIL = REIL+REIL*(TC-25.)*.0001*EXP((.000025*WAVL)**.25)           FL41670
      IF (WAVL.GT..0007) GO TO 50                                        FL41680
      R1 = DOP(WAVL,1.79907,3352.27,99.914E+04,15.1963E+04,1639.,50483.5 FL41690
     *   ,9246.27,588.24,84.4697E+04,10.7615E+05)                        FL41700
      R1 = R1+R1*(TC-25.)*.0001*EXP((.000025*WAVL)**.25)                 FL41710
      REIL = R1*(.0007-WAVL)/.0001+REIL*(WAVL-.0006)/.0001               FL41720
      GO TO 50                                                           FL41730
   40 REIL = DOP(WAVL,1.79907,3352.27,99.914E+04,15.1963E+04,1639.,      FL41740
     *   50483.5,9246.27,588.24,84.4697E+04,10.7615E+05)                 FL41750
      REIL = REIL+REIL*(TC-25.)*.0001*EXP((.000025*WAVL)**.25)           FL41760
C                                                                        FL41770
C     * *  TABLE 2 WATER PG. 1840                                        FL41780
C                                                                        FL41790
   50 IF (WAVL.GE..3) GO TO 180                                          FL41800
      IF (WAVL.GE..03) GO TO 60                                          FL41810
      GO TO 70                                                           FL41820
   60 AIMAG = AIMAG+AB(WAVL,.25,300.,.47,3.)+AB(WAVL,.39,17.,.45,1.3)+AB FL41830
     *   (WAVL,.41,62.,.35,1.7)                                          FL41840
      GO TO 180                                                          FL41850
   70 IF (WAVL.GE..0062) GO TO 80                                        FL41860
      GO TO 90                                                           FL41870
   80 AIMAG = AIMAG+AB(WAVL,.41,62.,.35,1.7)+AB(WAVL,.39,17.,.45,1.3)+AB FL41880
     *   (WAVL,.25,300.,.4,2.)                                           FL41890
      GO TO 180                                                          FL41900
   90 IF (WAVL.GE..0017) GO TO 100                                       FL41910
      GO TO 110                                                          FL41920
  100 AIMAG = AIMAG+AB(WAVL,.39,17.,.45,1.3)+AB(WAVL,.41,62.,.22,1.8)+AB FL41930
     *   (WAVL,.25,300.,.4,2.)                                           FL41940
      GO TO 180                                                          FL41950
  110 IF (WAVL.GE..00061) GO TO 120                                      FL41960
      GO TO 130                                                          FL41970
  120 AIMAG = AIMAG+AB(WAVL,.12,6.1,.042,.6)+AB(WAVL,.39,17.,.165,2.4)+  FL41980
     *   AB(WAVL,.41,62.,.22,1.8)                                        FL41990
      GO TO 180                                                          FL42000
  130 IF (WAVL.GE..000495) GO TO 140                                     FL42010
      GO TO 150                                                          FL42020
  140 AIMAG = AIMAG+AB(WAVL,.01,4.95,.05,1.)+AB(WAVL,.12,6.1,.009,2.)    FL42030
      GO TO 180                                                          FL42040
  150 IF (WAVL.GE..000297) GO TO 160                                     FL42050
      GO TO 170                                                          FL42060
  160 AIMAG = AIMAG+AB(WAVL,.27,2.97,.04,2.)+AB(WAVL,.01,4.95,.06,1.)    FL42070
      GO TO 180                                                          FL42080
  170 AIMAG = AIMAG+AB(WAVL,.27,2.97,.025,2.)+AB(WAVL,.01,4.95,.06,1.)   FL42090
  180 CONTINUE                                                           FL42100
      RETURN                                                             FL42110
C                                                                        FL42120
  900 FORMAT(///,30X,'ATTEMPTING TO EVALUATE FOR A WAVELENGTH LESS THAN  FL42130
     *ONE MICRON',//)                                                    FL42140
  905 FORMAT(///,30X,'ATTEMPTING TO EVALUATE FOR A TEMPERATURE LESS THAN FL42150
     * -20. DEGREES CENTIGRADE',//)                                      FL42160
C                                                                        FL42170
      END                                                                FL42180
C
C     *****************************************************************
C
      SUBROUTINE DEBYE(WAVL,TC,KEY,RE,AI)                                FL42190
C                                                                        FL42200
C     CC                                                                 FL42210
C     CC    CALCULATES WAVENUMBER DEPENDENCE OF DIELECTRIC CONSTANT      FL42220
C     CC    OF WATER                                                     FL42230
C     CC                                                                 FL42240
C                                                                        FL42250
      T = TC+273.15                                                      FL42260
      IF (KEY.NE.0) GO TO 10                                             FL42270
      GO TO 20                                                           FL42280
   10 EFIN = 5.27137+.0216474*TC-.00131198*TC*TC                         FL42290
C                                                                        FL42300
C     CC   ALPHA=-16.8129/T+.0609265                                     FL42310
C                                                                        FL42320
      TAU = .00033836*EXP(2513.98/T)                                     FL42330
C                                                                        FL42340
C     CC   SIG=12.5664E+08                                               FL42350
C                                                                        FL42360
      ES = 78.54*(1.-.004579*(TC-25.)+.0000119*(TC-25.)**2-.000000028*   FL42370
     *   (TC-25.)**3)                                                    FL42380
      GO TO 30                                                           FL42390
   20 EFIN = 3.168                                                       FL42400
C                                                                        FL42410
C     CC   ALPHA=.00023*TC*TC+.0052*TC+.288                              FL42420
C     CC   SIG=1.26*EXP(-12500./(T*1.9869))                              FL42430
C                                                                        FL42440
      TAU = 9.990288E-05*EXP(13200./(T*1.9869))                          FL42450
      ES = 3.168+.15*TC*TC+2.5*TC+200.                                   FL42460
   30 C1 = TAU/WAVL                                                      FL42470
C                                                                        FL42480
C     CC                                                                 FL42490
C     CC    TEMPORARY FIX TO CLASSICAL DEBYE EQUATION                    FL42500
C     CC    TO HANDLE ZERO CM-1 PROBLEM                                  FL42510
C     CC                                                                 FL42520
C     CC   ALPHA=0.0                                                     FL42530
C     CC   SIG=0.0                                                       FL42540
C     CC                                                                 FL42550
C     CC   C2=1.5708*ALPHA                                               FL42560
C     CC   DEM=1.+2.*C1**(1.-ALPHA)*SIN(C2)+C1**(2.*(1.-ALPHA))          FL42570
C     CC   E1=EFIN+(ES-EFIN)*(1.+(C1**(1.-ALPHA)*SIN(C2)))/DEM           FL42580
C     CC   IF(KEY.NE.0.AND.WAVL.GE.300.) E1=87.53-0.3956*TC              FL42590
C     CC   IF(KEY.NE.0 .AND. WAVL.GE.300.) E1=ES                         FL42600
C     CC   E2=(ES-EFIN)*C1**(1.-ALPHA)*COS(C2)/DEM+SIG*WAVL/18.8496E+10  FL42610
C     CC                                                                 FL42620
C     CC    PERMANENT FIX TO CLSSICAL DEBYE EQUATION                     FL42630
C     CC    TO HANDLE ZERO CM-1 PROBLEM                                  FL42640
C     CC                                                                 FL42650
C                                                                        FL42660
      E1 = EFIN+(ES-EFIN)/(1.0+C1**2)                                    FL42670
C                                                                        FL42680
C     CC                                                                 FL42690
C                                                                        FL42700
      E2 = ((ES-EFIN)*C1)/(1.0+C1**2)                                    FL42710
C                                                                        FL42720
C     CC                                                                 FL42730
C                                                                        FL42740
      RE = SQRT((E1+SQRT(E1*E1+E2*E2))/2.)                               FL42750
      AI = -E2/(2.*RE)                                                   FL42760
      RETURN                                                             FL42770
      END                                                                FL42780
      FUNCTION DOP(WAVL,A,CEN1,B,C,CEN2,D,E,CEN3,F,G)                    FL42790
C                                                                        FL42800
C     CC                                                                 FL42810
C     CC    DESCRIBES THE REAL PART OF THE DIELECTRIC CONSTANT           FL42820
C     CC                                                                 FL42830
C                                                                        FL42840
      V = 1./WAVL                                                        FL42850
      V2 = V*V                                                           FL42860
      H1 = CEN1**2-V2                                                    FL42870
      H2 = CEN2**2-V2                                                    FL42880
      H3 = CEN3**2-V2                                                    FL42890
      DOP = SQRT(A+B*H1/(H1*H1+C*V2)+D*H2/(H2*H2+E*V2)+F*H3/(H3*H3+G*V2) FL42900
     *   )                                                               FL42910
      RETURN                                                             FL42920
      END                                                                FL42930
      FUNCTION AB(WAVL,A,CEN,B,C)                                        FL42940
C                                                                        FL42950
C     CC                                                                 FL42960
C     CC    DESCRIBES THE IMAGINARY PART OF THE DIELECTRIC CONSTANT      FL42970
C     CC                                                                 FL42980
C                                                                        FL42990
      AB = -A*EXP(-ABS(( LOG10(10000.*WAVL/CEN)/B))**C)                  FL43000
      RETURN                                                             FL43010
      END                                                                FL43020
      FUNCTION GAMFOG(FREQ,T,RHO)                                        FL43030
C                                                                        FL43040
C     COMPUTES ATTENUATION OF EQUIVALENT LIQUID WATER CONTENT            FL43050
C     IN CLOUDS OR FOG IN DB/KM                                          FL43060
C     CONVERTED TO NEPERS BY NEW CONSTANT 1.885                          FL43070
C                                                                        FL43080
C     FREQ = WAVENUMBER (INVERSE CM)                                     FL43090
C     T    = TEMPERATURE (DEGREES KELVIN)                                FL43100
C     RHO  = EQUIVALENT LIQUID CONTENT  (G/CUBIC METER)                  FL43110
C     CINDEX=COMPLEX DIELECTRIC CONSTANT M  FROM INDEX                   FL43120
C     WAVL = WAVELENGTH IN CM                                            FL43130
C                                                                        FL43140
      COMPLEX CINDEX                                                     FL43150
      IF (RHO.GT.0.) GO TO 10                                            FL43160
      GAMFOG = 0.                                                        FL43170
      RETURN                                                             FL43180
   10 CONTINUE                                                           FL43190
      KEY = 1                                                            FL43200
      WAVL = 1.0/FREQ                                                    FL43210
      TC = T-273.2                                                       FL43220
C                                                                        FL43230
C     CC                                                                 FL43240
C     CC    CHANGE TEMP SO THAT MINIMUM IS -20.0 CENT.                   FL43250
C     CC                                                                 FL43260
C                                                                        FL43270
      TC = MAX(TC,-20.0)                                                 FL43280
      CALL INDX (WAVL,TC,KEY,REIL,AIMAK)                                 FL43290
      CINDEX = CMPLX(REIL,AIMAK)                                         FL43300
C                                                                        FL43310
C     CC                                                                 FL43320
C     CC   ATTENUATION = 6.0*PI*FREQ*RHO*IMAG(-K)                        FL43330
C     CC    6.0*PI/10. = 1.885 (THE FACTOR OF 10 IS FOR UNITS CONVERSION FL43340
C     CC                                                                 FL43350
C     GAMFOG=8.1888*FREQ*RHO*AIMAG( -  (CINDEX**2-1)/(CINDEX**2+2))      FL43360
C                                                                        FL43370
      GAMFOG = 1.885*FREQ*RHO*AIMAG(-(CINDEX**2-1)/(CINDEX**2+2))        FL43380
      RETURN                                                             FL43390
      END                                                                FL43400
      FUNCTION AITK(ARG,VAL,X,NDIM)                                      FL43410
C                                                                        FL43420
C     IBM SCIENTIFIC SUBROUTINE                                          FL43430
C     AITKEN INTERPOLATION ROUTINE                                       FL43440
C                                                                        FL43450
      DIMENSION ARG(NDIM),VAL(NDIM)                                      FL43460
c
      IF (NDIM.gt.1) then
C                                                                        FL43480
C     START OF AITKEN-LOOP                                               FL43490
C                                                                        FL43500
   10 DO 30 J = 2, NDIM                                                  FL43510
         IEND = J-1                                                      FL43520
         DO 20 I = 1, IEND                                               FL43530
            H = ARG(I)-ARG(J)                                            FL43540
            IF (H.eq.0) go to 70
            VAL(J) = (VAL(I)*(X-ARG(J))-VAL(J)*(X-ARG(I)))/H             FL43560
 20      continue
 30   CONTINUE                                                           FL43570
C                                                                        FL43580
C     END OF AITKEN-LOOP                                                 FL43590
C                                                                        FL43600
      endif


   40 J = NDIM                                                           FL43610
   50 AITK = VAL(J)                                                      FL43620
   60 RETURN                                                             FL43630
C                                                                        FL43640
C     THERE ARE TWO IDENTICAL ARGUMENT VALUES IN VECTOR ARG              FL43650
C                                                                        FL43660
   70 IER = 3                                                            FL43670
      J = IEND                                                           FL43680
      GO TO 50                                                           FL43690
      END                                                                FL43700
      FUNCTION GMRAIN(FREQ,T,RATE)                                       FL43710
C                                                                        FL43720
C     COMPUTES ATTENUATION OF CONDENSED WATER IN FORM OF RAIN            FL43730
C                                                                        FL43740
C     FREQ = WAVENUMBER (CM-1)                                           FL43750
C     T    = TEMPERATURE (DEGREES KELVIN)                                FL43760
C     RATE = PRECIPITATION RATE (MM/HR)                                  FL43770
C     WVLTH = WAVELENGTH IN CM                                           FL43780
C                                                                        FL43790
C     TABLES ATTAB AND FACTOR CALCULATED FROM FULL MIE THEORY            FL43800
C     UTILIZING MARSHALL-PALMER SIZE DISTRIBUTION WITH RAYS INDEX        FL43810
C     OF REFRACTION                                                      FL43820
C                                                                        FL43830
C     ATTAB IS ATTENUATION DATA TABLE IN NEPERS FOR 20 DEG CELSIUS       FL43840
C     WITH RADIATION FIELD REMOVED                                       FL43850
C                                                                        FL43860
C     WVNTBL IS WAVENUMBER TABLE FOR WAVENUMBERS USED IN TABLE ATTAB     FL43870
C     TMPTAB IS INTERPOLATION DATA TABLE FOR TEMPERATURES IN DEG KELVIN  FL43880
C                                                                        FL43890
C     TLMDA IS INTERPOLATION DATA TABLE FOR WAVELENGTH IN CM             FL43900
C     TFREQ IS INTERPOLATION DATA TABLE FOR WAVENUMBER IN CM-1           FL43910
C                                                                        FL43920
C     RATTAB IS RAIN RATE TABLE IN MM/HR                                 FL43930
C                                                                        FL43940
C     FACTOR IS TABLE OF TEMPERATURE CORRECTION FACTORS FOR              FL43950
C     TABLE ATTAB FOR REPRESENTATIVE RAINS WITHOUT RADIATION FIELD       FL43960
C                                                                        FL43970
C                                                                        FL43980
C     AITKEN INTERPOLATION SCHEME WRITTEN BY                             FL43990
C     E.T. FLORANCE O.N.R. PASADENA CA.                                  FL44000
C                                                                        FL44010
C                                                                        FL44020
      DIMENSION ATTAB1(35),ATTAB2(35),ATTAB3(35),ATTAB4(35),ATTAB5(35)   FL44030
      DIMENSION ATTAB6(35),ATTAB7(35),ATTAB8(35),ATTAB9(35)              FL44040
      DIMENSION ATTAB(35,9),WVLTAB(27),RATTAB(9),FACTOR(5,8,5)           FL44050
      DIMENSION X(4),Y(4),ATTN(4),RATES(4)                               FL44060
C                                                                        FL44070
C     CC   DIMENSION X(3),Y(3),ATTN(3),RATES(3)                          FL44080
C                                                                        FL44090
      DIMENSION TMPTAB(5),TLMDA(6),FACIT(5),TFACT(5)                     FL44100
      DIMENSION TFREQ(8),WVNTBL(35)                                      FL44110
      DIMENSION FACEQ1(5,8),FACEQ2(5,8),FACEQ3(5,8),FACEQ4(5,8)          FL44120
      DIMENSION FACEQ5(5,8)                                              FL44130
      EQUIVALENCE (ATTAB1(1),ATTAB(1,1)),(ATTAB2(1),ATTAB(1,2))          FL44140
      EQUIVALENCE (ATTAB3(1),ATTAB(1,3)),(ATTAB4(1),ATTAB(1,4))          FL44150
      EQUIVALENCE (ATTAB5(1),ATTAB(1,5)),(ATTAB6(1),ATTAB(1,6))          FL44160
      EQUIVALENCE (ATTAB7(1),ATTAB(1,7)),(ATTAB8(1),ATTAB(1,8))          FL44170
      EQUIVALENCE (ATTAB9(1),ATTAB(1,9))                                 FL44180
      EQUIVALENCE (FACEQ1(1,1),FACTOR(1,1,1))                            FL44190
      EQUIVALENCE (FACEQ2(1,1),FACTOR(1,1,2))                            FL44200
      EQUIVALENCE (FACEQ3(1,1),FACTOR(1,1,3))                            FL44210
      EQUIVALENCE (FACEQ4(1,1),FACTOR(1,1,4))                            FL44220
      EQUIVALENCE (FACEQ5(1,1),FACTOR(1,1,5))                            FL44230
c
      DATA WVLTAB/.03,.033,.0375,.043,.05,.06,.075,.1,.15,.2,.25,.3,.5,  FL44240
     *.8,1.,2.,3.,4.,5.,5.5,6.,6.5,7.,8.,9.,10.,15./                     FL44250
      DATA WVNTBL/ 0.0000,                                               FL44260
     *    .0667,.1000,.1111,.1250,.1429,.1538,                           FL44270
     *  .1667,.1818,.2000,.2500,.3333,.5000,1.0000,                      FL44280
     * 1.2500,2.0000,3.3333,4.0000,5.0000,6.6667,10.0000,                FL44290
     * 13.3333,16.6667,20.0000,23.2558,26.6667,30.3030,33.3333,          FL44300
     * 50.0,80.0,120.0,180.0,250.0,300.0,350.0/                          FL44310
      DATA RATTAB /.25,1.25,2.5,5.,12.5,25.,50.,100.,150./               FL44320
      DATA TLMDA/.03,.1,.5,1.25,3.2,10./                                 FL44330
      DATA TFREQ/0.0,0.1,0.3125,0.8,2.0,10.0,33.3333,350.0/              FL44340
      DATA TMPTAB/273.15,283.15,293.15,303.15,313.15/                    FL44350
      DATA ATTAB1/                                                       FL44360
     * 1.272E+00,1.332E+00,1.361E+00,1.368E+00,1.393E+00,1.421E+00,      FL44370
     * 1.439E+00,1.466E+00,1.499E+00,1.541E+00,1.682E+00,1.951E+00,      FL44380
     * 2.571E+00,3.575E+00,3.808E+00,4.199E+00,3.665E+00,3.161E+00,      FL44390
     * 2.462E+00,1.632E+00,8.203E-01,4.747E-01,3.052E-01,2.113E-01,      FL44400
     * 1.551E-01,1.168E-01,8.958E-02,7.338E-02,3.174E-02,1.178E-02,      FL44410
     * 5.016E-03,2.116E-03,1.123E-03,8.113E-04,6.260E-04/                FL44420
      DATA ATTAB2/                                                       FL44430
     * 4.915E+00,5.257E+00,5.518E+00,5.632E+00,5.807E+00,6.069E+00,      FL44440
     * 6.224E+00,6.452E+00,6.756E+00,7.132E+00,8.453E+00,1.132E+01,      FL44450
     * 1.685E+01,2.177E+01,2.246E+01,2.156E+01,1.470E+01,1.167E+01,      FL44460
     * 8.333E+00,5.089E+00,2.356E+00,1.320E+00,8.315E-01,5.705E-01,      FL44470
     * 4.151E-01,3.119E-01,2.385E-01,1.955E-01,8.373E-02,3.138E-02,      FL44480
     * 1.351E-02,5.789E-03,3.090E-03,2.236E-03,1.725E-03/                FL44490
      DATA ATTAB3/                                                       FL44500
     * 8.798E+00,9.586E+00,1.023E+01,1.049E+01,1.093E+01,1.159E+01,      FL44510
     * 1.205E+01,1.263E+01,1.343E+01,1.450E+01,1.832E+01,2.627E+01,      FL44520
     * 3.904E+01,4.664E+01,4.702E+01,4.152E+01,2.542E+01,1.959E+01,      FL44530
     * 1.363E+01,8.087E+00,3.660E+00,2.028E+00,1.274E+00,8.710E-01,      FL44540
     * 6.340E-01,4.757E-01,3.634E-01,2.971E-01,1.275E-01,4.795E-02,      FL44550
     * 2.072E-02,8.936E-03,4.780E-03,3.460E-03,2.670E-03/                FL44560
      DATA ATTAB4/                                                       FL44570
     * 1.575E+01,1.750E+01,1.914E+01,1.991E+01,2.108E+01,2.276E+01,      FL44580
     * 2.399E+01,2.561E+01,2.785E+01,3.097E+01,4.204E+01,6.334E+01,      FL44590
     * 8.971E+01,9.853E+01,9.609E+01,7.718E+01,4.290E+01,3.220E+01,      FL44600
     * 2.188E+01,1.271E+01,5.641E+00,3.110E+00,1.947E+00,1.327E+00,      FL44610
     * 9.657E-01,7.242E-01,5.539E-01,4.528E-01,1.942E-01,7.335E-02,      FL44620
     * 3.181E-02,1.380E-02,7.394E-03,5.354E-03,4.132E-03/                FL44630
      DATA ATTAB5/                                                       FL44640
     * 3.400E+01,3.927E+01,4.523E+01,4.796E+01,5.207E+01,5.886E+01,      FL44650
     * 6.383E+01,7.060E+01,8.005E+01,9.360E+01,1.381E+02,2.069E+02,      FL44660
     * 2.620E+02,2.534E+02,2.366E+02,1.673E+02,8.285E+01,6.059E+01,      FL44670
     * 4.013E+01,2.280E+01,9.939E+00,5.439E+00,3.400E+00,2.315E+00,      FL44680
     * 1.685E+00,1.263E+00,9.664E-01,7.914E-01,3.397E-01,1.288E-01,      FL44690
     * 5.611E-02,2.450E-02,1.316E-02,9.536E-03,7.360E-03/                FL44700
      DATA ATTAB6/                                                       FL44710
     * 6.087E+01,7.347E+01,8.886E+01,9.653E+01,1.081E+02,1.283E+02,      FL44720
     * 1.435E+02,1.649E+02,1.947E+02,2.346E+02,3.543E+02,4.991E+02,      FL44730
     * 5.705E+02,5.048E+02,4.510E+02,2.900E+02,1.335E+02,9.607E+01,      FL44740
     * 6.269E+01,3.520E+01,1.519E+01,8.295E+00,5.182E+00,3.529E+00,      FL44750
     * 2.569E+00,1.927E+00,1.474E+00,1.208E+00,5.191E-01,1.975E-01,      FL44760
     * 8.627E-02,3.784E-02,2.037E-02,1.476E-02,1.139E-02/                FL44770
      DATA ATTAB7/                                                       FL44780
     * 1.090E+02,1.396E+02,1.811E+02,2.029E+02,2.396E+02,3.039E+02,      FL44790
     * 3.536E+02,4.189E+02,5.081E+02,6.217E+02,9.038E+02,1.165E+03,      FL44800
     * 1.212E+03,9.731E+02,8.330E+02,4.901E+02,2.123E+02,1.507E+02,      FL44810
     * 9.718E+01,5.408E+01,2.316E+01,1.264E+01,7.896E+00,5.377E+00,      FL44820
     * 3.915E+00,2.939E+00,2.249E+00,1.844E+00,7.940E-01,3.029E-01,      FL44830
     * 1.327E-01,5.846E-02,3.151E-02,2.284E-02,1.763E-02/                FL44840
      DATA ATTAB8/                                                       FL44850
     * 1.950E+02,2.703E+02,3.904E+02,4.614E+02,5.825E+02,7.909E+02,      FL44860
     * 9.475E+02,1.142E+03,1.380E+03,1.656E+03,2.237E+03,2.610E+03,      FL44870
     * 2.500E+03,1.820E+03,1.491E+03,8.103E+02,3.336E+02,2.344E+02,      FL44880
     * 1.495E+02,8.273E+01,3.524E+01,1.922E+01,1.203E+01,8.182E+00,      FL44890
     * 5.961E+00,4.477E+00,3.429E+00,2.812E+00,1.216E+00,4.651E-01,      FL44900
     * 2.043E-01,9.033E-02,4.874E-02,3.534E-02,2.728E-02/                FL44910
      DATA ATTAB9/                                                       FL44920
     * 2.742E+02,4.012E+02,6.353E+02,7.829E+02,1.027E+03,1.439E+03,      FL44930
     * 1.725E+03,2.071E+03,2.475E+03,2.909E+03,3.738E+03,4.104E+03,      FL44940
     * 3.776E+03,2.589E+03,2.070E+03,1.078E+03,4.326E+02,3.023E+02,      FL44950
     * 1.918E+02,1.059E+02,4.499E+01,2.454E+01,1.539E+01,1.045E+01,      FL44960
     * 7.615E+00,5.722E+00,4.384E+00,3.596E+00,1.561E+00,5.978E-01,      FL44970
     * 2.630E-01,1.165E-01,6.292E-02,4.562E-02,3.522E-02/                FL44980
      DATA FACEQ1/                                                       FL44990
     * 1.606,1.252,1.000, .816, .680,1.603,1.246,1.000, .817, .684,      FL45000
     * 1.444,1.207,1.000, .838, .694,1.016, .985,1.000,1.034,1.058,      FL45010
     *  .950, .976,1.000,1.034,1.068, .922, .956,1.000,1.044,1.090,      FL45020
     *  .932, .966,1.000,1.034,1.068, .957, .978,1.000,1.022,1.044/      FL45030
      DATA FACEQ2/                                                       FL45040
     * 1.606,1.252,1.000, .816, .680,1.612,1.256,1.000, .817, .684,      FL45050
     * 1.193,1.101,1.000, .889, .769, .885, .927,1.000,1.086,1.175,      FL45060
     *  .941, .976,1.000,1.024,1.047, .932, .966,1.000,1.034,1.079,      FL45070
     *  .932, .966,1.000,1.034,1.068, .957, .978,1.000,1.022,1.044/      FL45080
      DATA FACEQ3/                                                       FL45090
     * 1.606,1.252,1.000, .816, .680,1.621,1.256,1.000, .817, .673,      FL45100
     *  .969, .995,1.000, .982, .940, .895, .937,1.000,1.075,1.143,      FL45110
     *  .950, .976,1.000,1.024,1.036, .932, .966,1.000,1.034,1.079,      FL45120
     *  .932, .966,1.000,1.034,1.068, .957, .978,1.000,1.022,1.044/      FL45130
      DATA FACEQ4/                                                       FL45140
     * 1.606,1.252,1.000, .816, .680,1.631,1.265,1.000, .807, .662,      FL45150
     *  .848, .927,1.000,1.044,1.079, .922, .956,1.000,1.055,1.111,      FL45160
     *  .950, .976,1.000,1.013,1.036, .932, .966,1.000,1.034,1.079,      FL45170
     *  .932, .966,1.000,1.034,1.068, .957, .978,1.000,1.022,1.044/      FL45180
      DATA FACEQ5/                                                       FL45190
     * 1.606,1.252,1.000, .816, .680,1.603,1.265,1.000, .807, .662,      FL45200
     *  .820, .918,1.000,1.075,1.132, .941, .966,1.000,1.034,1.079,      FL45210
     *  .960, .976,1.000,1.013,1.036, .932, .966,1.000,1.034,1.079,      FL45220
     *  .932, .966,1.000,1.034,1.068, .957, .978,1.000,1.022,1.044/      FL45230
      DATA RATLIM /.05/                                                  FL45240
C                                                                        FL45250
C     GIVE ZERO ATTN IF RATE FALLS BELOW LIMIT                           FL45260
C                                                                        FL45270
      IF (RATE.GT.RATLIM) GO TO 10                                       FL45280
      GMRAIN = 0.                                                        FL45290
      RETURN                                                             FL45300
   10 WVLTH = 1.0/FREQ                                                   FL45310
C                                                                        FL45320
C     CC   JMAX=3                                                        FL45330
C                                                                        FL45340
      JMAX = 4                                                           FL45350
C                                                                        FL45360
C     CC   IF(WVLTH.GT.WVLTAB(1)) GO TO      14                          FL45370
C     CC   ILOW=0                                                        FL45380
C     CC   JMAX=2                                                        FL45390
C     CC   GO TO 18                                                      FL45400
C     CC   THIS DO LOOP IS 2 LESS THAN NO. OF WVLTAB INPUT               FL45410
C     CC14 DO 15 I=2,25                                                  FL45420
C                                                                        FL45430
      DO 20 I = 3, 33                                                    FL45440
C                                                                        FL45450
C        CC   IF(WVLTH.LT.(.5*(WVLTAB(I)+WVLTAB(I+1)))) GO TO 16         FL45460
C                                                                        FL45470
         IF (FREQ.LT.WVNTBL(I)) GO TO 30                                 FL45480
   20 CONTINUE                                                           FL45490
C                                                                        FL45500
C     CC   SET ILOW EQUAL TO 1 LESS THAN DO MAX                          FL45510
C     CC   ILOW=24                                                       FL45520
C                                                                        FL45530
      I = 34                                                             FL45540
C                                                                        FL45550
C     CC   GO TO 18                                                      FL45560
C     CC16 ILOW = I-2                                                    FL45570
C                                                                        FL45580
   30 ILOW = I-3                                                         FL45590
C                                                                        FL45600
C     CC   DO 190 I=2,7                                                  FL45610
C                                                                        FL45620
      DO 40 K = 3, 7                                                     FL45630
C                                                                        FL45640
C        CC   IF (RATE. LT.(.5*(RATTAB(I)+RATTAB(I+1))))GO TO 195        FL45650
C                                                                        FL45660
         IF (RATE.LT.RATTAB(K)) GO TO 50                                 FL45670
   40 CONTINUE                                                           FL45680
C                                                                        FL45690
C     CC   KMIN=6                                                        FL45700
C                                                                        FL45710
      K = 8                                                              FL45720
C                                                                        FL45730
C     CC   GO TO 198                                                     FL45740
C     C195 KMIN=I-2                                                      FL45750
C                                                                        FL45760
   50 KMIN = K-3                                                         FL45770
      DO 60 J = 1, JMAX                                                  FL45780
         IJ = ILOW+J                                                     FL45790
         X(J) = WVNTBL(IJ)                                               FL45800
   60 CONTINUE                                                           FL45810
C                                                                        FL45820
C     INTERPOLATE                                                        FL45830
C     CC   Z = - LOG(FREQ)                                               FL45840
C     CC   DO 25 K=1,3                                                   FL45850
C                                                                        FL45860
      DO 80 K = 1, 4                                                     FL45870
         KJ = KMIN+K                                                     FL45880
         RATES(K) = RATTAB(KJ)                                           FL45890
         DO 70 J = 1, JMAX                                               FL45900
            IJ = ILOW+J                                                  FL45910
            Y(J) =  LOG(ATTAB(IJ,KJ))                                    FL45920
   70    CONTINUE                                                        FL45930
         ATTN(K) = EXP(AITK(X,Y,FREQ,JMAX))                              FL45940
   80 CONTINUE                                                           FL45950
C                                                                        FL45960
C     APPLY TEMPERATURE CORRECTION                                       FL45970
C                                                                        FL45980
      DO 90 I = 2, 5                                                     FL45990
         IF (T.LT.TMPTAB(I)) GO TO 100                                   FL46000
   90 CONTINUE                                                           FL46010
      ILOW = 4                                                           FL46020
      GO TO 110                                                          FL46030
  100 ILOW = I-1                                                         FL46040
  110 CONTINUE                                                           FL46050
      DO 120 J = 2, 8                                                    FL46060
         IF (FREQ.LT.TFREQ(J)) GO TO 130                                 FL46070
  120 CONTINUE                                                           FL46080
C                                                                        FL46090
C     CC   JLOW IS 2 LESS THAN DO MAX                                    FL46100
C                                                                        FL46110
      JLOW = 6                                                           FL46120
      GO TO 140                                                          FL46130
  130 JLOW = J-2                                                         FL46140
  140 CONTINUE                                                           FL46150
      DO 160 K = 1, 2                                                    FL46160
         DO 150 J = 1, 2                                                 FL46170
C                                                                        FL46180
C           INTERPOLATE IN TEMPERATURE                                   FL46190
C           CC   KJ=(KMIN/2)+K                                           FL46200
C                                                                        FL46210
            KJ = K+(KMIN+1)/2                                            FL46220
            JI = JLOW+J                                                  FL46230
            FAC = ((TMPTAB(ILOW)-T)*FACTOR(ILOW+1,JI,KJ)+(T-TMPTAB(ILOW+ FL46240
     *         1))*FACTOR(ILOW,JI,KJ))/(TMPTAB(ILOW)-TMPTAB(ILOW+1))     FL46250
            JI = JLOW+3-J                                                FL46260
            FACIT(J) = (TFREQ(JI)-FREQ)*FAC                              FL46270
  150    CONTINUE                                                        FL46280
         TFACT(K) = (FACIT(2)-FACIT(1))/(TFREQ(JLOW+1)-TFREQ(JLOW+2))    FL46290
  160 CONTINUE                                                           FL46300
C                                                                        FL46310
C     COMPUTE ATTENUATION (DB/KM)                                        FL46320
C     CC   KJ=2*KMIN/2+1                                                 FL46330
C                                                                        FL46340
      KJ = 2*((KMIN+1)/2)+1                                              FL46350
C                                                                        FL46360
C     CC   GMRAIN=AITK(RATES,ATTN,RATE,3)*                               FL46370
C                                                                        FL46380
      GMRAIN = AITK(RATES,ATTN,RATE,4)*((RATE-RATTAB(KJ))*TFACT(2)+      FL46390
     *   (RATTAB(KJ+2)-RATE)*TFACT(1))/(RATTAB(KJ+2)-RATTAB(KJ))         FL46400
C                                                                        FL46410
C     CC                                                                 FL46420
C     CC    APPLY CONVERSION TO NEPERS                                   FL46430
C     CC                                                                 FL46440
C                                                                        FL46450
      RETURN                                                             FL46460
      END                                                                FL46470
C
C     ******************************************************************
C
      SUBROUTINE CIRRUS(CTHIK,CALT,ISEED,CPROB,MODEL)                    FL46480
C                                                                        FL46490
C     ****************************************************************** FL46500
C     *  ROUTINE TO GENERATE ALTITUDE PROFILES OF CIRRUS DENSITY         FL46510
C     *  PROGRAMMED BY   M.J. POST                                       FL46520
C     *                  R.A. RICHTER        NOAA/WPL                    FL46530
C     *                                      BOULDER, COLORADO           FL46540
C     *                                      01/27/1981                  FL46550
C     *                                                                  FL46560
C     *  INPUTS|                                                         FL46570
C     *           CHTIK    -  CIRRUS THICKNESS (KM)                      FL46580
C     *                       0 = USE THICKNESS STATISTICS               FL46590
C     *                       .NE. 0 = USER DEFINES THICKNESS            FL46600
C     *                                                                  FL46610
C     *           CALT     -  CIRRUS BASE ALTITUDE (KM)                  FL46620
C     *                       0 = USE CALCULATED VALUE                   FL46630
C     *                       .NE. 0 = USER DEFINES BASE ALTITUDE        FL46640
C     *                                                                  FL46650
C     *           ICIR     -  CIRRUS PRESENCE FLAG                       FL46660
C     *                       0 = NO CIRRUS                              FL46670
C     *                       .NE. 0 = USE CIRRUS PROFILE                FL46680
C     *                                                                  FL46690
C     *           MODEL    -  ATMOSPHERIC MODEL                          FL46700
C     *                       1-5  AS IN MAIN PROGRAM                    FL46710
C     *                       MODEL = 0,6,7 NOT USED SET TO 2            FL46720
C     *                                                                  FL46730
C     *           ISEED    -  RANDOM NUMBER INITIALIZATION FLAG.         FL46740
C     *                       0 = USE DEFAULT MEAN VALUES FOR CIRRUS     FL46750
C     *                       .NE. 0 = INITIAL VALUE OF SEED FOR RANF    FL46760
C     *                       FUNCTION. CHANGE SEED VALUE EACH RUN FOR   FL46770
C     *                       DIFFERENT RANDOM NUMBER SEQUENCES. THIS    FL46780
C     *                       PROVIDES FOR STATISTICAL DETERMINATION     FL46790
C     *                       OF CIRRUS BASE ALTITUDE AND THICKNESS.     FL46800
C     *                                                                  FL46810
C     *  OUTPUTS|                                                        FL46820
C     *         CTHIK        -  CIRRUS THICKNESS (KM)                    FL46830
C     *         CALT         -  CIRRUS BASE ALTITUDE (KM)                FL46840
C     *         DENSTY(16,I) -  ARRAY, ALTITUDE PROFILE OF CIRRUS DENSIT FL46850
C     *         CPROB        -  CIRRUS PROBABILITY                       FL46860
C     *                                                                  FL46870
C     ****************************************************************** FL46880
C                                                                        FL46890
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),VH(16),TX(16),W(16)       FL46900
      COMMON WPATH(IM2,16),TBBY(IM2)                                     FL46910
      COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),VX2(47),AWCCON(5)          FL46920
C
      CHARACTER*8      HMOD                                              FL46930
C
C
      COMMON /CMN/ HMOD(3),ZN(MXZMD),PN(MXZMD),TN(MXZMD),RFNDXM(MXZMD),
     *         ZP(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2), PPSUM(IM2),
     *          TPSUM(IM2),RHOPSM(IM2),IMMAX,WGM(MXZMD),DEMW(MXZMD),
     *          AMTP(MXMOL,MXPDIM)                                        
C
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,N_LVL,JH1              FL46980
      COMMON /MODEL/ ZMDL (MXZMD),PM(MXZMD),TM(MXZMD),
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL46990
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)          FL47000
C
      DIMENSION CBASE(5,2),TSTAT(11),PTAB(5),CAMEAN(5)                   FL47010
      DIMENSION CBASE1(5),CBASE2(5)                                      FL47020
C
      EQUIVALENCE (CBASE1(1),CBASE(1,1)),(CBASE2(1),CBASE(1,2))          FL47030
C                                                                        FL47040
C     ISEED IS INTEGER*4                                                 FL47050
C                                                                        FL47060
      INTEGER*4 ISEED,IDUM                                               FL47070
C                                                                        FL47080
      DATA  CAMEAN           / 11.0, 10.0, 8.0, 7.0, 5.0 /               FL47090
      DATA  PTAB           / 0.8, 0.4, 0.5, 0.45, 0.4/                   FL47100
      DATA  CBASE1            / 7.5, 7.3, 4.5, 4.5, 2.5 /                FL47110
      DATA  CBASE2            /16.5,13.5,14.0, 9.5,10.0 /                FL47120
      DATA  TSTAT             / 0.0,.291,.509,.655,.764,.837,.892,       FL47130
     * 0.928, 0.960, 0.982, 1.00 /                                       FL47140
C                                                                        FL47150
C     SET CIRRUS PROBABILITY AND PROFILE TO ALL ZEROES                   FL47160
C                                                                        FL47170
      CPROB = 0.0                                                        FL47180
      MDL = MODEL                                                        FL47190
C                                                                        FL47200
      DO 10 I = 1, 68                                                    FL47210
         DENSTY(16,I) = 0.                                               FL47220
   10 CONTINUE                                                           FL47230
C                                                                        FL47240
C     CHECK IF USER WANTS TO USE A THICKNESS VALUE HE PROVIDES, CALCULAT FL47250
C     A STATISTICAL THICKNESS, OR USE A MEAN THICKNESS (ISEED = 0).      FL47260
C     DEFAULTED MEAN CIRRUS THICKNESS IS 1.0 KM.                         FL47270
C                                                                        FL47280
      IF (CTHIK.GT.0.0) GO TO 40                                         FL47290
      IF (ISEED.NE.0) GO TO 20                                           FL47300
      CTHIK = 1.0                                                        FL47310
      GO TO 40                                                           FL47320
C                                                                        FL47330
C     > CALCULATE CLOUD THICKNESS USING LOWTRAN CIRRUS THICKNESS STATIST FL47340
C     > NOTE - THIS ROUTINE USES A UNIFORM RANDOM NUMBER GENERATOR       FL47350
C     > WHICH RETURNS A NUMBER BETWEEN 0 AND 1.                          FL47360
C     >                                                                  FL47370
C                                                                        FL47380
   20 IDUM = -ISEED                                                      FL47390
C                                                                        FL47400
      URN = RANDM(IDUM)                                                  FL47410
      DO 30 I = 1, 10                                                    FL47420
         IF (URN.GE.TSTAT(I).AND.URN.LT.TSTAT(I+1)) CTHIK = I-1          FL47430
   30 CONTINUE                                                           FL47440
      CTHIK = CTHIK/2.0+RANDM(IDUM)/2.0                                  FL47450
C                                                                        FL47460
C     DENCIR IS CIRRUS DENSITY IN KM-1                                   FL47470
C                                                                        FL47480
   40 DENCIR = 0.07*CTHIK                                                FL47490
C                                                                        FL47500
C     BASE HEIGHT CALCULATIONS                                           FL47510
C                                                                        FL47520
      IF (MODEL.LT.1.OR.MODEL.GT.5) MDL = 2                              FL47530
      CPROB = 100.0*PTAB(MDL)                                            FL47540
C                                                                        FL47550
      HMAX = CBASE(MDL,2)-CTHIK                                          FL47560
      BRANGE = HMAX-CBASE(MDL,1)                                         FL47570
      IF (CALT.GT.0.0) GO TO 60                                          FL47580
      IF (ISEED.NE.0) GO TO 50                                           FL47590
      CALT = CAMEAN(MDL)                                                 FL47600
      GO TO 60                                                           FL47610
   50 CALT = BRANGE*RANDM(IDUM)+CBASE(MDL,1)                             FL47620
C                                                                        FL47630
C     PUT CIRRUS DENSITY IN CORRECT ALTITUDE BINS. IF MODEL = 7,         FL47640
C     INTERPOLATE EH(16,I) FOR NON-STANDARD ALTITUDE BOUNDARIES.         FL47650
C                                                                        FL47660
   60 TOP = CALT+CTHIK                                                   FL47670
      BOTTOM = CALT                                                      FL47680
      IF (TOP.LT.ZMDL(1)) RETURN                                         FL47690
      IF (BOTTOM.GT.ZMDL(ML)) RETURN                                     FL47700
      IML = ML-1                                                         FL47710
      DO 70 I = 1, IML                                                   FL47720
         ZMIN = ZMDL(I)                                                  FL47730
         ZMAX = ZMDL(I+1)                                                FL47740
         DENOM = ZMAX-ZMIN                                               FL47750
         IF (BOTTOM.LE.ZMIN.AND.TOP.GE.ZMAX) DENSTY(16,I) = DENCIR       FL47760
         IF (BOTTOM.GE.ZMIN.AND.TOP.LT.ZMAX) DENSTY(16,I) = DENCIR*CTHIK FL47770
     *      /DENOM                                                       FL47780
         IF (BOTTOM.GE.ZMIN.AND.TOP.GE.ZMAX.AND.BOTTOM.LT.ZMAX)          FL47790
     *       DENSTY(16,I) = DENCIR*(ZMAX-BOTTOM)/DENOM                   FL47800
         IF (BOTTOM.LT.ZMIN.AND.TOP.LE.ZMAX.AND.TOP.GT.ZMIN) DENSTY(16,I FL47810
     *      ) = DENCIR*(TOP-ZMIN)/DENOM                                  FL47820
   70 CONTINUE                                                           FL47830
      RETURN                                                             FL47840
      END                                                                FL47850
C
C     *****************************************************************
C
      SUBROUTINE VSA(IHAZE,VIS,CEILHT,DEPTH,ZINVHT,Z,RH,AHAZE,IH)        FL47860
C                                                                        FL47870
C     VERTICAL STRUCTURE ALGORITHM                                       FL47880
C                                                                        FL47890
C     FROM ATMOSPHERIC SCIENCES LAB (U.S. ARMY)                          FL47900
C     WHITE SANDS N.M.                                                   FL47910
C                                                                        FL47920
C     CREATES A PROFILE OF AEROSOL DENSITY NEAR THE GROUND,INCLUDING     FL47930
C     CLOUDS AND FOG                                                     FL47940
C                                                                        FL47950
C     THESE PROFILES ARE AT 9 HEIGHTS BETWEEN 0 KM AND 2 KM              FL47960
C                                                                        FL47970
C                                                                        FL47980
C     ***VISIBILITY IS ASSUMED TO BE THE SURFACE VISIBILITY***           FL47990
C                                                                        FL48000
C     IHAZE  = THE TYPE OF AEROSOL                                       FL48010
C     VIS    = VISIBILITY IN KM                                          FL48020
C     CEILHT = THE CLOUD CEILING HEIGHT IN KM                            FL48030
C     DEPTH  = THE CLOUD/FOG DEPTH IN KM                                 FL48040
C     ZINVHT = THE HEIGHT OF INVERSION OR BOUNDARY LAYER IN KM           FL48050
C                                                                        FL48060
C     VARIABLES USED IN VSA                                              FL48070
C                                                                        FL48080
C     ZC     = CLOUD CEILING HEIGHT IN M                                 FL48090
C     ZT     = CLOUD DEPTH IN M                                          FL48100
C     ZINV   = INVERSION HEIGHT IN M                                     FL48110
C     SEE BELOW FOR MORE INFORMATION ABOUT ZC, ZT, AND ZINV              FL48120
C     D      = INITIAL EXTINCTION AT THE SURFACE (D=3.912/VIS)           FL48130
C     ZALGO  = THE DEPTH OF THE LAYER FOR THE ALGORITHM                  FL48140
C                                                                        FL48150
C     OUTPUT FROM VSA:                                                   FL48160
C                                                                        FL48170
C     Z      = HEIGHT IN KM                                              FL48180
C     RH     = RELATIVE HUMIDITY AT HEIGHT Z IN PERCENT                  FL48190
C     AHAZE  = EXTINCTION AT HEIGHT Z IN KM**-1                          FL48200
C     IH     = AEROSAL TYPE FOR HEIGHT Z                                 FL48210
C     HMAX   = MAXIMUM HEIGHT IN KM USED IN VSA, NOT NECESSARILY 2.0 KM  FL48220
C                                                                        FL48230
C                                                                        FL48240
C     THE SLANT PATH CALCULATION USES THE FOLLOWING FUNCTION:            FL48250
C                                                                        FL48260
C     EXT55=A*EXP(B*EXP(C*Z))                                            FL48270
C                                                                        FL48280
C     WHERE 'Z' IS THE HEIGHT IN KILOMETERS,                             FL48290
C     'A' IS A FUNCTION OF EXT55 AT Z=0.0 AND IS ALWAYS POSITIVE,        FL48300
C     'B' AND 'C' ARE FUNCTIONS OF CLOUD CONDITIONS AND SURFACE          FL48310
C     VISIBILITY (EITHER A OR B CAN BE POSITIVE OR NEGATIVE),            FL48320
C     'EXT55' IS THE VISIBILE EXTINCTION COEFFICIENT IN KM**-1.          FL48330
C                                                                        FL48340
C     THEREFORE, THERE ARE 4 CASES DEPENDING ON THE SIGNS OF 'B' AND 'C' FL48350
C     CEILHT AND ZINVHT ARE USED AS SWITCHES TO DETERMINE WHICH CASE     FL48360
C     TO USE.  THE SURFACE EXTINCTION 'D' IS CALCULATED FROM THE         FL48370
C     VISIBILITY USING  D=3.912/VIS-0.012 AS FOLLOWS-                    FL48380
C                                                                        FL48390
C     CASE=1  FOG/CLOUD CONDITIONS                                       FL48400
C     'B' LT 0.0, 'C' LT 0.0                                             FL48410
C     'D' GE 7.064 KM**-1                                                FL48420
C     FOR A CLOUD 7.064 KM**-1 IS THE BOUNDARY VALUE AT                  FL48430
C     THE CLOUD BASE AND 'Z' IS THE VERTICAL DISTANCE                    FL48440
C     INTO THE CLOUD.                                                    FL48450
C     VARIABLE USED:   DEPTH                                             FL48460
C     ** DEFAULT:  DEPTH OF FOG/CLOUD IS 0.2 KM WHEN                     FL48470
C     'DEPTH' IS 0.0                                                     FL48480
C                                                                        FL48490
C     =2  CLOUD CEILING PRESENT                                          FL48500
C     'B' GT 0.0, 'C' GT 0.0                                             FL48510
C     'D' GT 0.398 KM**-1 IS CASE 2 FOR HAZY/FOG                         FL48520
C     SURFACE CONDITIONS                                                 FL48530
C     'D' LE 0.398 KM**-1 IS CASE 2' FOR CLEAR/HAZY                      FL48540
C     SURFACE CONDITIONS                                                 FL48550
C     VARIABLE USED:   CEILHT (MUST BE GE 0.0)                           FL48560
C     ** DEFAULTS:  CASE 2 - CEILHT IS CALCULATED FROM                   FL48570
C     SURFACE EXTINCTION OR                                              FL48580
C     CASE 2' - CEILHT IS 1.8 KM WHEN                                    FL48590
C     'CEILHT' IS 0.0                                                    FL48600
C                                                                        FL48610
C     =3  RADIATION FOG OR INVERSION OR BOUNDARY LAYER PRESENT           FL48620
C     'B' LT 0.0, 'C' GT 0.0                                             FL48630
C     VIS LE 2.0 KM DEFAULTS TO A RADIATION FOG AT THE                   FL48640
C     GROUND AND OVERRIDES INPUT BOUNDARY AEROSOL TYPE                   FL48650
C     VIS GT 2.0 KM FOR AN INVERSION OR BOUNDARY LAYER                   FL48660
C     WITH INPUT BOUNDARY AEROSOL TYPE                                   FL48670
C     ** IHAZE=9 (RADIATION FOG) ALWAYS DEFAULTS TO A                    FL48680
C     RADIATION FOG NO MATTER WHAT THE VISIBILITY IS.                    FL48690
C     SWITCH VARIABLE: CEILHT (MUST BE LT 0.0)                           FL48700
C     VARIABLE USED:   ZINVHT (MUST BE GE 0.0)                           FL48710
C     ** CEILHT MUST BE LT 0.0 FOR ZINVHT TO BE USED **                  FL48720
C     HOWEVER, IF DEPTH IS GT 0.0 AND ZINVHT IS EQ 0.0,                  FL48730
C     THE PROGRAM WILL SUBSTITUTE DEPTH FOR ZINVHT.                      FL48740
C     ** DEFAULT:  FOR A RADIATION FOG ZINVHT IS 0.2 KM                  FL48750
C     FOR AN INVERSION LAYER ZINVHT IS 2.0 KM                            FL48760
C                                                                        FL48770
C     =4  NO CLOUD CEILING, INVERSION LAYER, OR BOUNDARY                 FL48780
C     LAYER PRESENT, I.E. CLEAR SKIES                                    FL48790
C     EXTINCTION PROFILE CONSTANT WITH HEIGHT                            FL48800
C                                                                        FL48810
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL48820
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL48830
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL48840
C
      DIMENSION Z(10),RH(10),AHAZE(10),IH(10)                            FL48850
      DIMENSION AA(2),CC(2),EE(4),A(2),B(2),C(2),FAC1(9),FAC2(9)         FL48860
C
      REAL KMTOM                                                         FL48870
C
      DATA AA/135.,0.3981/,CC/-0.030,0.0125/,KMTOM/1000.0/,CR/0.35/      FL48880
C                                                                        FL48890
C     THE LAST 3 VALUES OF EE BELOW ARE EXTINCTIONS FOR VISIBILITIES     FL48900
C     EQUAL TO 5.0, 23.0, AND 50.0 KM, RESPECTIVELY.                     FL48910
C                                                                        FL48920
      DATA EE/7.064,0.7824,0.17009,0.07824/                              FL48930
      DATA FAC1/0.0,0.03,0.05,0.075,0.1,0.18,0.3,0.45,1.0/               FL48940
      DATA FAC2/0.0,0.03,0.1,0.18,0.3,0.45,0.6,0.78,1.0/                 FL48950
C
      WRITE (IPR,900)                                                    FL48960
C                                                                        FL48970
C     UPPER LIMIT ON VERTICAL DISTANCE - 2 KM                            FL48980
C                                                                        FL48990
      ZHIGH = 2000.                                                      FL49000
      HMAX = ZHIGH                                                       FL49010
      IF (VIS.GT.0.0) GO TO 10                                           FL49020
C                                                                        FL49030
C     DEFAULT FOR VISIBILITY DEPENDS ON THE VALUE OF IHAZE.              FL49040
C                                                                        FL49050
      IF (IHAZE.EQ.8) VIS = 0.2                                          FL49060
      IF (IHAZE.EQ.9) VIS = 0.5                                          FL49070
      IF (IHAZE.EQ.2.OR.IHAZE.EQ.5) VIS = 5.0                            FL49080
      IF (IHAZE.EQ.1.OR.IHAZE.EQ.4.OR.IHAZE.EQ.7) VIS = 23.0             FL49090
      IF (IHAZE.EQ.6) VIS = 50.0                                         FL49100
C                                                                        FL49110
C     IF(IHAZE.EQ.3)VIS=??????                                           FL49120
C                                                                        FL49130
   10 D = 3.912/VIS-0.012                                                FL49140
C                                                                        FL49150
      ZC = CEILHT*KMTOM                                                  FL49160
      IF (ZC.GT.CR) THEN                                                 FL49170
         SZ = ZC-CR                                                      FL49180
      ELSE                                                               FL49190
         SZ = 0.                                                         FL49200
      ENDIF                                                              FL49210
      ZT = DEPTH*KMTOM                                                   FL49220
      ZINV = ZINVHT*KMTOM                                                FL49230
C                                                                        FL49240
C     IHAZE=9 (RADIATION FOG) IS ALWAYS CALCULATED AS A RADIATION FOG.   FL49250
C                                                                        FL49260
      IF (IHAZE.EQ.9) ZC = -1.0                                          FL49270
C                                                                        FL49280
C     ALSO, CHECK TO SEE IF THE FOG DEPTH FOR A RADIATION FOG            FL49290
C     WAS INPUT TO DEPTH INSTEAD OF THE CORRECT VARIABLE ZINVHT.         FL49300
C                                                                        FL49310
      IF (IHAZE.EQ.9.AND.ZT.GT.0.0.AND.ZINV.EQ.0.0) ZINV = ZT            FL49320
C                                                                        FL49330
C     'IC' DEFINES WHICH CASE TO USE.                                    FL49340
C                                                                        FL49350
      IC = 2                                                             FL49360
      IF (D.GE.EE(1).AND.ZC.GE.0.0) IC = 1                               FL49370
C                                                                        FL49380
      IF (ZC.LT.0.0.AND.IC.EQ.2) IC = 3                                  FL49390
      IF (ZINV.LT.0.0.AND.IC.EQ.3) IC = 4                                FL49400
C                                                                        FL49410
C     'ICC' IS FOR THE TWO CASES:  2 AND 2'.                             FL49420
C                                                                        FL49430
      ICC = 0                                                            FL49440
      IF (IC.EQ.2) ICC = 1                                               FL49450
      IF (D.LE.AA(2).AND.IC.EQ.2) ICC = 2                                FL49460
      K = 1                                                              FL49470
      IF (ICC.EQ.2) GO TO 40                                             FL49480
      GO TO (20,30,50,60), IC                                            FL49490
C                                                                        FL49500
C     CASE 1:  DEPTH FOG/CLOUD; INCREASING EXTINCTION WITH HEIGHT FROM   FL49510
C     CLOUD/FOG BASE TO CLOUD/FOG TOP.                                   FL49520
C                                                                        FL49530
   20 CONTINUE                                                           FL49540
      IF (ZC.LT.HMAX.AND.IC.EQ.2) K = 2                                  FL49550
C                                                                        FL49560
C     IC=-1 WHEN A CLOUD IS PRESENT AND THE PATH GOES INTO IT.           FL49570
C     USE CASE 2 OR 2' BELOW CLOUD AND CASE 1 INSIDE IT.                 FL49580
C                                                                        FL49590
      IF (K.EQ.2) IC = (-1)                                              FL49600
C                                                                        FL49610
C     THE BASE OF THE CLOUD HAS AN EXTINCTION COEFFICIENT OF 7.064 KM-1. FL49620
C                                                                        FL49630
      IF (K.EQ.2) D = EE(1)                                              FL49640
      A(K) = AA(1)                                                       FL49650
C                                                                        FL49660
C     IF THE SURFACE EXTINCTION IS GREATER THAN THE UPPER LIMIT OF 92.1  FL49670
C     KM**-1, RUN THE ALGORITHM WITH AN UPPER LIMIT OF 'D+10'.           FL49680
C                                                                        FL49690
      IF (D.GE.AA(1)) A(K) = D+10.0                                      FL49700
      C(K) = CC(1)                                                       FL49710
      IF (ZT.LE.0.0) WRITE (IPR,940)                                     FL49720
      IF (ZT.LE.0.0) WRITE (IPR,945)                                     FL49730
      IF (ZT.GT.0.0) WRITE (IPR,955) ZT                                  FL49740
C                                                                        FL49750
C     IF THE DISTANCE FROM THE GROUND TO THE CLOUD/FOG TOP IS LESS       FL49760
C     THAN 2.0 KM, VSA WILL ONLY CALCULATE UP TO THE CLOUD TOP.          FL49770
C                                                                        FL49780
      IF (ZT.LE.0.0) ZT = 200.                                           FL49790
      HMAX =   MIN(ZT+ZC,HMAX)                                           FL49800
      GO TO 60                                                           FL49810
C                                                                        FL49820
C     CASE 2:  HAZY/LIGHTLY FOGGY; INCREASING EXTINCTION WITH HEIGHT     FL49830
C     UP TO THE CLOUD BASE.                                              FL49840
C                                                                        FL49850
   30 A(K) = AA(2)                                                       FL49860
      E = EE(1)                                                          FL49870
      IF (ZC.EQ.0.0) WRITE (IPR,905)                                     FL49880
      IF (ZC.EQ.0.0) CEIL =  LOG( LOG(E/A(K))/( LOG(D/A(K))))/CC(2)      FL49890
      IF (ZC.EQ.0.0) WRITE (IPR,935) CEIL                                FL49900
      IF (ZC.GT.0.0) WRITE (IPR,950) ZC                                  FL49910
      IF (ZC.EQ.0.0) ZC = CEIL                                           FL49920
      GO TO 60                                                           FL49930
C                                                                        FL49940
C     CASE 2':  CLEAR/HAZY; INCREASING EXTINCTION WITH HEIGHT, BUT LESS  FL49950
C     SO THAN CASE 2, UP TO THE CLOUD BASE.                              FL49960
C                                                                        FL49970
   40 A(K) = D*0.9                                                       FL49980
      E = EE(1)                                                          FL49990
      IF (ZC.EQ.0.0) WRITE (IPR,905)                                     FL50000
      IF (ZC.EQ.0.0) WRITE (IPR,920)                                     FL50010
      IF (ZC.GT.0.0) WRITE (IPR,950) ZC                                  FL50020
      IF (ZC.EQ.0.0) ZC = 1800.                                          FL50030
      GO TO 60                                                           FL50040
C                                                                        FL50050
C     CASE 3:  NO CLOUD CEILING BUT A RADIATION FOG OR AN INVERSION      FL50060
C     OR BOUNDARY LAYER PRESENT; DECREASING EXTINCTION WITH              FL50070
C     HEIGHT UP TO THE HEIGHT OF THE FOG OR LAYER.                       FL50080
C                                                                        FL50090
   50 A(K) = D*1.1                                                       FL50100
      E = EE(3)                                                          FL50110
      IF (IHAZE.EQ.2.OR.IHAZE.EQ.5) E = EE(2)                            FL50120
      IF (IHAZE.EQ.6.OR.(VIS.GT.2.0.AND.IHAZE.NE.9)) E = EE(4)           FL50130
      IF (E.GT.D) E = D*0.99999                                          FL50140
      IF (ZT.GT.0.0.AND.ZINV.EQ.0.0.AND.VIS.LE.2.0) ZINV = ZT            FL50150
      IF (ZINV.EQ.0.0.AND.VIS.GT.2.0.AND.IHAZE.NE.9) WRITE (IPR,910)     FL50160
      IF (ZINV.EQ.0.0.AND.(VIS.LE.2.0.OR.IHAZE.EQ.9)) WRITE (IPR,915)    FL50170
      IF (ZINV.EQ.0.0.AND.(VIS.LE.2.0.OR.IHAZE.EQ.9)) WRITE (IPR,945)    FL50180
      IF (ZINV.GT.0.0.AND.VIS.GT.2.0.AND.IHAZE.NE.9) WRITE (IPR,960)     FL50190
     *    ZINV                                                           FL50200
      IF (ZINV.GT.0.0.AND.(VIS.LE.2.0.OR.IHAZE.EQ.9)) WRITE (IPR,965)    FL50210
     *    ZINV                                                           FL50220
      IF (ZINV.EQ.0.0.AND.VIS.GT.2.0.AND.IHAZE.NE.9) ZINV = 2000         FL50230
      IF (ZINV.EQ.0.0.AND.(VIS.LE.2.0.OR.IHAZE.EQ.9)) ZINV = 200         FL50240
      HMAX =   MIN(ZINV,HMAX)                                            FL50250
      ZC = 0.0                                                           FL50260
C                                                                        FL50270
C     CASE 4:  NO CLOUD CEILING OR INVERSION LAYER;                      FL50280
C     CONSTANT EXTINCTION WITH HEIGHT.                                   FL50290
C                                                                        FL50300
   60 IF (IC.NE.4) B(K) =  LOG(D/A(K))                                   FL50310
      IF (IC.EQ.4) WRITE (IPR,970)                                       FL50320
      IF (IC.EQ.2) THEN                                                  FL50330
         C(K) =  LOG( LOG(E/A(K))/B(K))/(ZC-SZ)                          FL50340
      ENDIF                                                              FL50350
      IF (IC.EQ.3) C(K) =  LOG( LOG(E/A(K))/B(K))/ZINV                   FL50360
      IF (ZC.LT.HMAX.AND.K.EQ.1.AND.IC.EQ.2) GO TO 20                    FL50370
      IF (IC.EQ.2) HMAX =   MIN(ZC,HMAX)                                 FL50380
      ZALGO = HMAX                                                       FL50390
      IF (IC.LT.0) ZALGO = ZC                                            FL50400
      WRITE (IPR,925)                                                    FL50410
      IF (IC.LT.0) K = 1                                                 FL50420
C                                                                        FL50430
      DO 70 I = 1, 9                                                     FL50440
         IF (IC.LT.0.AND.I.EQ.5) K = 2                                   FL50450
         IF (IC.LT.0.AND.I.EQ.5) ZALGO = HMAX-ZC                         FL50460
         Z(I) = ZALGO*(1.0-FAC2(10-I))                                   FL50470
         IF (IC.EQ.1) Z(I) = ZALGO*FAC1(I)                               FL50480
         IF (IC.EQ.4) Z(I) = ZALGO* REAL(I-1)/8.0                        FL50490
         IF (IC.LT.0.AND.I.LT.5) Z(I) = ZALGO*(1.0-FAC2(11-2*I))         FL50500
         IF (IC.LT.0.AND.I.GE.5) Z(I) = ZALGO*FAC1(2*I-9)                FL50510
C                                                                        FL50520
C        IF(IC.LT.0.AND.(I.EQ.7.OR.I.EQ.8))Z(I)=ZALGO*FAC1(2*I-10)       FL50530
C        C    IF(IC.NE.4)AHAZE(I)=A(K)*EXP(B(K)*EXP(C(K)*Z(I)))          FL50540
C        C    IF(IC.EQ.4)AHAZE(I)=D                                      FL50550
C                                                                        FL50560
         IF (IC.NE.4) THEN                                               FL50570
            IF (Z(I).GT.SZ) THEN                                         FL50580
               AHAZE(I) = A(K)*EXP(B(K)*EXP(C(K)*Z(I)-SZ))               FL50590
            ELSE                                                         FL50600
               AHAZE(I) = A(K)*EXP(B(K)*EXP(C(K)*Z(I)))                  FL50610
            ENDIF                                                        FL50620
         ELSE                                                            FL50630
            AHAZE(I) = D                                                 FL50640
         ENDIF                                                           FL50650
         IF (IC.LE.0.AND.I.GE.5) Z(I) = Z(I)+ZC                          FL50660
         Z(I) = Z(I)/KMTOM                                               FL50670
         RH(I) = 6.953* LOG(AHAZE(I))+86.407                             FL50680
         IF (AHAZE(I).GE.EE(1)) RH(I) = 100.0                            FL50690
         VISIB = 3.912/(AHAZE(I)+0.012)                                  FL50700
         IH(I) = IHAZE                                                   FL50710
C                                                                        FL50720
C        IF A RADIATION FOG IS PRESENT (I.E. VIS<=2.0 KM AND IC=3),      FL50730
C        IH IS SET TO 9 FOR ALL LEVELS.                                  FL50740
C                                                                        FL50750
         IF (VISIB.LE.2.0.AND.IC.EQ.3) IH(I) = 9                         FL50760
C                                                                        FL50770
C        FOR A DEPTH FOG/CLOUD CASE, IH=8 DENOTING AN ADVECTION FOG.     FL50780
C                                                                        FL50790
         IF (IC.EQ.1.OR.(IC.LT.0.AND.I.GE.5)) IH(I) = 8                  FL50800
         WRITE (IPR,930) Z(I),RH(I),AHAZE(I),VISIB,IH(I)                 FL50810
   70 CONTINUE                                                           FL50820
      HMAX = HMAX/KMTOM                                                  FL50830
      RETURN                                                             FL50840
C                                                                        FL50850
  900 FORMAT('0 VERTICAL STRUCTURE ALGORITHM (VSA) USED')                FL50860
  905 FORMAT(' ',50X,'CLOUD CEILING HEIGHT UNKNOWN')
  910 FORMAT(' ',50X,'INVERSION OR BOUNDARY LAYER HEIGHT UNKNOWN',/,
     *  ' ',50X,'VSA WILL USE A DEFAULT OF 2000.0 METERS',/)
  915 FORMAT(' ',50X,'RADIATION FOG DEPTH UNKNOWN')
  920 FORMAT(' ',50X,'VSA WILL USE A DEFAULT OF 1800.0 METERS',/)     
  925 FORMAT(5X,'HEIGHT(KM)',5X,'R.H.(%)',5X,'EXTINCTION(KM-1)',
     *   5X,'VIS(3.912/EXTN)',5X,'IHAZE',/)  
  930 FORMAT(7X,F7.4,7X,F5.1,8X,E12.4,11X,F7.4,10X,I2)                   FL50940
  935 FORMAT(' ',39X,'VSA WILL USE A CALCULATED VALUE OF ',F7.1,    
     *       ' METERS',/)                                           
  940 FORMAT(' ',50X,'CLOUD DEPTH UNKNOWN')                           
  945 FORMAT(' ',50X,'VSA WILL USE A DEFAULT OF 200.0 METERS',/)      
  950 FORMAT(' ',50X,'CLOUD CEILING HEIGHT IS ',F9.1,' METERS',/)     
  955 FORMAT(' ',50X,'CLOUD DEPTH IS ,F14.1,7H METERS',/)             
  960 FORMAT(' ',50X,'INVERSION OR BOUNDARY LAYER HEIGHT IS ',F7.1,   
     * ' METERS',/)                                                   
  965 FORMAT(' ',50X,'DEPTH OF RADIATION FOG IS ',F7.1,' METERS',/)   
  970 FORMAT(' ',50X,'THERE IS NO INVERSION OR BOUNDARY LAYER OR ',   
     * 'CLOUD PRESENT',/)                                             
C                                                                        FL51060
      END                                                                FL51070
C
C     *****************************************************************
C
      SUBROUTINE EXABIN                                                  FL51080
C                                                                        FL51090
C     LOADS EXTINCTION, ABSORPTION AND ASYMMETRY COEFFICIENTS            FL51100
C     FOR THE FOUR AEROSOL ALTITUDE REGIONS                              FL51110
C                                                                        FL51120
C     MODIFIED FOR ASYMMETRY - JAN 1986 (A.E.R. INC.)                    FL51130
C                                                                        FL51140
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),VH(16),TX(16),W(16)
      COMMON WPATH(IM2,16),TBBY(IM2)
      COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),VX0(47),AWCCON(5)
C
      CHARACTER*8      HMOD
C
      COMMON /CMN/ HMOD(3),ZM(MXZMD),PF(MXZMD),TF(MXZMD),RFNDXM(MXZMD),
     *          ZP(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),
     *          TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),
     *          AMTP(MXMOL,MXPDIM)                                        
C
      COMMON /LCRD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB   FL51150
      COMMON /LCRD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,    FL51160
     *     RAINRT                                                        FL51170
      COMMON /LCRD2D/ IREG(4),ALTB(4),IREGC(4)                           FL51180
      COMMON /LCRD3/ H1,H2,ANGLE,RANGE,BETA,RE,LEN                       FL51190
      COMMON /LCRD4/ V1,V2,DV                                            FL51200
C                                                                        FL51210
      COMMON /EXTD  /  VX2(47),RUREXT(47,4),RURABS(47,4),RURSYM(47,4),   FL51300
     *     URBEXT(47,4),URBABS(47,4),URBSYM(47,4),OCNEXT(47,4),          FL51310
     *     OCNABS(47,4),OCNSYM(47,4),TROEXT(47,4),TROABS(47,4),          FL51320
     *     TROSYM(47,4),FG1EXT(47),FG1ABS(47),FG1SYM(47),                FL51330
     *     FG2EXT(47),FG2ABS(47),FG2SYM(47),BSTEXT(47),BSTABS(47),       FL51340
     *     BSTSYM(47),AVOEXT(47),AVOABS(47),AVOSYM(47),FVOEXT(47),       FL51350
     *     FVOABS(47),FVOSYM(47),DMEEXT(47),DMEABS(47),DMESYM(47),       FL51360
     *     CCUEXT(47),CCUABS(47),CCUSYM(47),CALEXT(47),CALABS(47),       FL51370
     *     CALSYM(47),CSTEXT(47),CSTABS(47),CSTSYM(47),CSCEXT(47),       FL51380
     *     CSCABS(47),CSCSYM(47),CNIEXT(47),CNIABS(47),CNISYM(47)        FL51390
      COMMON/CIRR/ CI64XT(47),CI64AB(47),CI64G(47),                      FL51400
     *     CIR4XT(47),CIR4AB(47),CIR4G(47)                               FL51410
C
      DIMENSION RHZONE(4)                                                FL51420
      DIMENSION ELWCR(4),ELWCU(4),ELWCM(4),ELWCT(4)                      FL51430
C
      DATA RHZONE/0.,70.,80.,99./                                        FL51440
      DATA ELWCR/3.517E-04,3.740E-04,4.439E-04,9.529E-04/                FL51450
      DATA ELWCM/4.675E-04,6.543E-04,1.166E-03,3.154E-03/                FL51460
      DATA ELWCU/3.102E-04,3.802E-04,4.463E-04,9.745E-04/                FL51470
      DATA ELWCT/1.735E-04,1.820E-04,2.020E-04,2.408E-04/                FL51480
      DATA AFLWC/1.295E-02/,RFLWC/1.804E-03/,CULWC/7.683E-03/            FL51490
      DATA ASLWC/4.509E-03/,STLWC/5.272E-03/,SCLWC/4.177E-03/            FL51500
      DATA SNLWC/7.518E-03/,BSLWC/1.567E-04/,FVLWC/5.922E-04/            FL51510
      DATA AVLWC/1.675E-04/,MDLWC/4.775E-04/                             FL51520
C
      DO 10 I = 1, 47                                                    FL51530
         VX0(I) = VX2(I)                                                 FL51540
   10 CONTINUE                                                           FL51550
      I1 = 1                                                             FL51560
      NB = 1                                                             FL51570
      NE = 46                                                            FL51580
C                                                                        FL51590
C     C    IF (IHAZE.EQ.7) I1=2                                          FL51600
C     C    IF(IHAZE.EQ.3) I1 = 2                                         FL51610
C     C    DO 185 M=I1,4                                                 FL51620
C                                                                        FL51630
      DO 260 M = 1, 4                                                    FL51640
C                                                                        FL51650
C        C    IF(ICLD.EQ.11.AND.M.EQ.2) GO TO 185                        FL51660
C                                                                        FL51670
         IF (IREG(M).NE.0) GO TO 260                                     FL51680
         ITA = ICH(M)                                                    FL51690
         ITC = ICH(M)-7                                                  FL51700
         ITAS = ITA                                                      FL51710
 47      IF (IREGC(M).NE.0) GO TO 190                                    FL51720
         WRH = W(15)                                                     FL51730
         IF (ICH(M).EQ.6.AND.M.NE.1) WRH = 70.                           FL51740
C                                                                        FL51750
C        THIS CODING  DOES NOT ALLOW TROP RH DEPENDENT  ABOVE EH(7,I)    FL51760
C        DEFAULTS TO TROPOSPHERIC AT 70. PERCENT                         FL51770
C                                                                        FL51780
         DO 20 I = 2, 4                                                  FL51790
            IF (WRH.LT.RHZONE(I)) GO TO 30                               FL51800
   20    CONTINUE                                                        FL51810
         I = 4                                                           FL51820
   30    II = I-1                                                        FL51830
         IF (WRH.GT.0.0.AND.WRH.LT.99.) X =  LOG(100.0-WRH)              FL51840
         X1 =  LOG(100.0-RHZONE(II))                                     FL51850
         X2 =  LOG(100.0-RHZONE(I))                                      FL51860
         IF (WRH.GE.99.0) X = X2                                         FL51870
         IF (WRH.LE.0.0) X = X1                                          FL51880
         DO 180 N = NB, NE                                               FL51890
            ITA = ITAS                                                   FL51900
            IF (ITA.EQ.3.AND.M.EQ.1) GO TO 40                            FL51910
            ABSC(M,N) = 0.                                               FL51920
            EXTC(M,N) = 0.                                               FL51930
            ASYM(M,N) = 0.0                                              FL51940
            IF (ITA.GT.6) GO TO 110                                      FL51950
            IF (ITA.LE.0) GO TO 180                                      FL51960
   40       IF (N.GE.41.AND.ITA.EQ.3) ITA = 4                            FL51970
C                                                                        FL51980
C           RH DEPENDENT AEROSOLS                                        FL51990
C                                                                        FL52000
            GO TO (50,50,60,70,80,90), ITA                               FL52010
   50       Y2 =  LOG(RUREXT(N,I))                                       FL52020
            Y1 =  LOG(RUREXT(N,II))                                      FL52030
            Z2 =  LOG(RURABS(N,I))                                       FL52040
            Z1 =  LOG(RURABS(N,II))                                      FL52050
            A2 =  LOG(RURSYM(N,I))                                       FL52060
            A1 =  LOG(RURSYM(N,II))                                      FL52070
            E2 =  LOG(ELWCR(I))                                          FL52080
            E1 =  LOG(ELWCR(II))                                         FL52090
            GO TO 100                                                    FL52100
   60       IF (M.GT.1) GO TO 70                                         FL52110
            A2 =  LOG(OCNSYM(N,I))                                       FL52120
            A1 =  LOG(OCNSYM(N,II))                                      FL52130
            A = A1+(A2-A1)*(X-X1)/(X2-X1)                                FL52140
            ASYM(M,N) = EXP(A)                                           FL52150
            E2 =  LOG(ELWCM(I))                                          FL52160
            E1 =  LOG(ELWCM(II))                                         FL52170
C                                                                        FL52180
C           NAVY MARITIME AEROSOL CHANGES TO MARINE IN MICROWAVE         FL52190
C           NO NEED TO DEFINE EQUIVALENT WATER                           FL52200
C                                                                        FL52210
            GO TO 180                                                    FL52220
   70       Y2 =  LOG(OCNEXT(N,I))                                       FL52230
            Y1 =  LOG(OCNEXT(N,II))                                      FL52240
            Z2 =  LOG(OCNABS(N,I))                                       FL52250
            Z1 =  LOG(OCNABS(N,II))                                      FL52260
            A2 =  LOG(OCNSYM(N,I))                                       FL52270
            A1 =  LOG(OCNSYM(N,II))                                      FL52280
            E2 =  LOG(ELWCM(I))                                          FL52290
            E1 =  LOG(ELWCM(II))                                         FL52300
            GO TO 100                                                    FL52310
   80       Y2 =  LOG(URBEXT(N,I))                                       FL52320
            Y1 =  LOG(URBEXT(N,II))                                      FL52330
            Z2 =  LOG(URBABS(N,I))                                       FL52340
            Z1 =  LOG(URBABS(N,II))                                      FL52350
            A2 =  LOG(URBSYM(N,I))                                       FL52360
            A1 =  LOG(URBSYM(N,II))                                      FL52370
            E2 =  LOG(ELWCU(I))                                          FL52380
            E1 =  LOG(ELWCU(II))                                         FL52390
            GO TO 100                                                    FL52400
   90       Y2 =  LOG(TROEXT(N,I))                                       FL52410
            Y1 =  LOG(TROEXT(N,II))                                      FL52420
            Z2 =  LOG(TROABS(N,I))                                       FL52430
            Z1 =  LOG(TROABS(N,II))                                      FL52440
            A2 =  LOG(TROSYM(N,I))                                       FL52450
            A1 =  LOG(TROSYM(N,II))                                      FL52460
            E2 =  LOG(ELWCT(I))                                          FL52470
            E1 =  LOG(ELWCT(II))                                         FL52480
  100       Y = Y1+(Y2-Y1)*(X-X1)/(X2-X1)                                FL52490
            ZK = Z1+(Z2-Z1)*(X-X1)/(X2-X1)                               FL52500
            A = A1+(A2-A1)*(X-X1)/(X2-X1)                                FL52510
            ABSC(M,N) = EXP(ZK)                                          FL52520
            EXTC(M,N) = EXP(Y)                                           FL52530
            ASYM(M,N) = EXP(A)                                           FL52540
            IF (N.EQ.1) EC = E1+(E2-E1)*(X-X1)/(X2-X1)                   FL52550
            IF (N.EQ.1) AWCCON(M) = EXP(EC)                              FL52560
            GO TO 180                                                    FL52570
  110       IF (ITA.GT.19) GO TO 170                                     FL52580
            IF (ITC.LT.1) GO TO 180                                      FL52590
            GO TO (120,130,180,140,150,160,150,160,140,140,160,170), ITC FL52600
  120       ABSC(M,N) = FG1ABS(N)                                        FL52610
            EXTC(M,N) = FG1EXT(N)                                        FL52620
            ASYM(M,N) = FG1SYM(N)                                        FL52630
            IF (N.EQ.1) AWCCON(M) = AFLWC                                FL52640
            GO TO 180                                                    FL52650
  130       ABSC(M,N) = FG2ABS(N)                                        FL52660
            EXTC(M,N) = FG2EXT(N)                                        FL52670
            ASYM(M,N) = FG2SYM(N)                                        FL52680
            IF (N.EQ.1) AWCCON(M) = RFLWC                                FL52690
            GO TO 180                                                    FL52700
  140       ABSC(M,N) = BSTABS(N)                                        FL52710
            EXTC(M,N) = BSTEXT(N)                                        FL52720
            ASYM(M,N) = BSTSYM(N)                                        FL52730
            IF (N.EQ.1) AWCCON(M) = BSLWC                                FL52740
            GO TO 180                                                    FL52750
  150       ABSC(M,N) = AVOABS(N)                                        FL52760
            EXTC(M,N) = AVOEXT(N)                                        FL52770
            ASYM(M,N) = AVOSYM(N)                                        FL52780
            IF (N.EQ.1) AWCCON(M) = AVLWC                                FL52790
            GO TO 180                                                    FL52800
  160       ABSC(M,N) = FVOABS(N)                                        FL52810
            EXTC(M,N) = FVOEXT(N)                                        FL52820
            ASYM(M,N) = FVOSYM(N)                                        FL52830
            ASYM(M,N) = DMESYM(N)                                        FL52840
            IF (N.EQ.1) AWCCON(M) = FVLWC                                FL52850
            GO TO 180                                                    FL52860
  170       ABSC(M,N) = DMEABS(N)                                        FL52870
            EXTC(M,N) = DMEEXT(N)                                        FL52880
            IF (N.EQ.1) AWCCON(M) = MDLWC                                FL52890
  180    CONTINUE                                                        FL52900
         GO TO 260                                                       FL52910
  190    CONTINUE                                                        FL52920
C                                                                        FL52930
C        CC                                                              FL52940
C        CC       SECTION TO LOAD EXTINCTION, ABSORPTION AND ASYMMETRY   FL52950
C        CC       COEFFICIENTS FOR CLOUD AND OR RAIN MODELS              FL52960
C        CC                                                              FL52970
C                                                                        FL52980
         DO 250 N = NB, NE                                               FL52990
            ABSC(M,N) = 0.0                                              FL53000
            EXTC(M,N) = 0.0                                              FL53010
            ASYM(M,N) = 0.0                                              FL53020
            IC = ICLD                                                    FL53030
            GO TO (200,210,220,230,240,220,240,240,200,200,200), IC      FL53040
  200       ABSC(M,N) = CCUABS(N)                                        FL53050
            EXTC(M,N) = CCUEXT(N)                                        FL53060
            ASYM(M,N) = CCUSYM(N)                                        FL53070
            IF (N.EQ.1) AWCCON(M) = CULWC                                FL53080
            GO TO 250                                                    FL53090
  210       ABSC(M,N) = CALABS(N)                                        FL53100
            EXTC(M,N) = CALEXT(N)                                        FL53110
            ASYM(M,N) = CALSYM(N)                                        FL53120
            IF (N.EQ.1) AWCCON(M) = ASLWC                                FL53130
            GO TO 250                                                    FL53140
  220       ABSC(M,N) = CSTABS(N)                                        FL53150
            EXTC(M,N) = CSTEXT(N)                                        FL53160
            ASYM(M,N) = CSTSYM(N)                                        FL53170
            IF (N.EQ.1) AWCCON(M) = STLWC                                FL53180
            GO TO 250                                                    FL53190
  230       ABSC(M,N) = CSCABS(N)                                        FL53200
            EXTC(M,N) = CSCEXT(N)                                        FL53210
            ASYM(M,N) = CSCSYM(N)                                        FL53220
            IF (N.EQ.1) AWCCON(M) = SCLWC                                FL53230
            GO TO 250                                                    FL53240
  240       ABSC(M,N) = CNIABS(N)                                        FL53250
            EXTC(M,N) = CNIEXT(N)                                        FL53260
            ASYM(M,N) = CNISYM(N)                                        FL53270
            IF (N.EQ.1) AWCCON(M) = SNLWC                                FL53280
  250    CONTINUE                                                        FL53290
  260 CONTINUE                                                           FL53300
      DO 270 N = 1, 47                                                   FL53310
         ABSC(5,N) = 0.                                                  FL53320
         EXTC(5,N) = 0.                                                  FL53330
         ASYM(5,N) = 0.                                                  FL53340
         AWCCON(5) = 0.                                                  FL53350
         IF (ICLD.EQ.18) THEN                                            FL53360
            ABSC(5,N) = CI64AB(N)                                        FL53370
            EXTC(5,N) = CI64XT(N)                                        FL53380
            ASYM(5,N) = CI64G(N)                                         FL53390
            AWCCON(5) = ASLWC                                            FL53400
         ENDIF                                                           FL53410
         IF (ICLD.EQ.19) THEN                                            FL53420
            ABSC(5,N) = CIR4AB(N)                                        FL53430
            EXTC(5,N) = CIR4XT(N)                                        FL53440
            ASYM(5,N) = CIR4G(N)                                         FL53450
            AWCCON(5) = ASLWC                                            FL53460
         ENDIF                                                           FL53470
  270 CONTINUE                                                           FL53480
      RETURN                                                             FL53490
C                                                                        FL53500
      END                                                                FL53510
C
C     ******************************************************************
C
      SUBROUTINE AEREXT (V,IK,RADFT)                                     FL53520
C                                                                        FL53530
C     INTERPOLATES AEROSOL EXTINCTION, ABSORPTION, AND ASYMMETRY         FL53540
C     COEFFICIENTS FOR THE WAVENUMBER, V, WITHOUT THE RADIATION FIELD.   FL53550
C                                                                        FL53560
C     MODIFIED FOR ASYMMETRY  - JAN 1986 (A.E.R. INC.)                   FL53570
C                                                                        FL53580
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),VH(16),TX(16),W(16)       FL53590
      COMMON WPATH(IM2,16),TBBY(IM2)                                     FL53600
      COMMON ABSC(5,47),EXTC(5,47),ASYC(5,47),VX2(47),AWCCON(5)          FL53610
C
      CHARACTER*8      HMOD                                              FL53620
C
      COMMON /CMN/ HMOD(3),ZM(MXZMD),PF(MXZMD),TF(MXZMD),RFNDXM(MXZMD),
     *          ZP(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),
     *          TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),
     *          AMTP(MXMOL,MXPDIM)                                        
C
      COMMON /LCRD1/ MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,TBOUND,SALB   FL53670
      COMMON /LCRD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,    FL53680
     *     RAINRT                                                        FL53690
      COMMON /LCRD3/ H1,H2,ANGLE,RANGE,BETA,RE,LEN                       FL53700
      COMMON /LCRD4/ V1,V2,DV                                            FL53710
      COMMON /CNTRL/ KMAX,M,IKMAX,NL,ML,IKLO,ISSGEO,IDUM1,IDUM2          FL53720
      COMMON /MODEL/ ZMDL(MXZMD),PM(MXZMD),TM(MXZMD),                    FL53730
     *     RFNDX(MXZMD),DENSTY(16,MXZMD),                                FL53740
     *     CLDAMT(MXZMD),RRAMT(MXZMD),EQLWC(MXZMD),HAZEC(MXZMD)
      COMMON /AER/ EXTV(5),ABSV(5),ASYV(5)                               FL53750
C                                                                        FL53760
C     CC                                                                 FL53770
C     CC    REDEFINE EXTC(47) AND ABSC(47) IF ALAM GT 200 MICRONS        FL53780
C     CC                                                                 FL53790
C                                                                        FL53800
      IF (V.LE.1.0E-5) GO TO 80                                          FL53810
      IF (RADFT.LE.0.) GO TO 80                                          FL53820
      IF (V.LE.33.333) GO TO 60                                          FL53830
C                                                                        FL53840
C     CC                                                                 FL53850
C     CC   COMPUTE INFRARED ATTENUATION COEFFICIENT                      FL53860
C     CC                                                                 FL53870
C                                                                        FL53880
      IF (V.LE.50.0) THEN                                                FL53890
         DO 10 MR = 1, 5                                                 FL53900
            EXTC(MR,47) = GAMFOG(33.333,TBBY(IK),AWCCON(MR))             FL53910
            ABSC(MR,47) = EXTC(MR,47)                                    FL53920
            ASYC(MR,47) = 0.0                                            FL53930
   10    CONTINUE                                                        FL53940
      ENDIF                                                              FL53950
      DO 20 I = 1, 4                                                     FL53960
         EXTV(I) = 0.                                                    FL53970
         ABSV(I) = 0.                                                    FL53980
         ASYV(I) = 0.                                                    FL53990
   20 CONTINUE                                                           FL54000
C                                                                        FL54010
C     C    IF (IHAZE.EQ.0) RETURN                                        FL54020
C                                                                        FL54030
      ALAM = 1.0E+4/V                                                    FL54040
      DO 30 N = 2, 47                                                    FL54050
         XD = ALAM-VX2(N)                                                FL54060
         IF (XD.lt.0.) go to 40
   30 CONTINUE                                                           FL54080
      N = 47                                                             FL54090
   40 VXD = VX2(N)-VX2(N-1)                                              FL54100
      DO 50 I = 1, 5                                                     FL54110
         ASYV(I) = (ASYC(I,N)-ASYC(I,N-1))*XD/VXD+ASYC(I,N)              FL54120
         EXTV(I) = (EXTC(I,N)-EXTC(I,N-1))*XD/VXD+EXTC(I,N)              FL54130
         ABSV(I) = (ABSC(I,N)-ABSC(I,N-1))*XD/VXD+ABSC(I,N)              FL54140
         EXTV(I) = EXTV(I)/RADFT                                         FL54150
         ABSV(I) = ABSV(I)/RADFT                                         FL54160
   50 CONTINUE                                                           FL54170
      RETURN                                                             FL54180
C                                                                        FL54190
C     CC                                                                 FL54200
C                                                                        FL54210
   60 CONTINUE                                                           FL54220
C                                                                        FL54230
C     CC    COMPUTE MICROWAVE ATTENUATION COEFFICIENTS                   FL54240
C     CC                                                                 FL54250
C                                                                        FL54260
      DO 70 I = 1, 5                                                     FL54270
         EXTV(I) = GAMFOG(V,TBBY(IK),AWCCON(I))                          FL54280
         ABSV(I) = EXTV(I)                                               FL54290
         ASYV(I) = 0.0                                                   FL54300
         EXTV(I) = EXTV(I)/RADFT                                         FL54310
         ABSV(I) = ABSV(I)/RADFT                                         FL54320
   70 CONTINUE                                                           FL54330
      RETURN                                                             FL54340
C                                                                        FL54350
C     CC                                                                 FL54360
C                                                                        FL54370
   80 CONTINUE                                                           FL54380
C                                                                        FL54390
C     CC    CALL FUNCTION TO OBTAIN LIMITING VALUE AS FREQ APPROACHES    FL54400
C     CC    ZERO USING RAY S MODIFIED DEBYE EQUATIONS                    FL54410
C     CC                                                                 FL54420
C     CC   EQL=EQLWC(IL)                                                 FL54430
C                                                                        FL54440
      DO 90 I = 1, 5                                                     FL54450
         EXTV(I) = ABSLIM(TBBY(IK),AWCCON(I))                            FL54460
         ABSV(I) = EXTV(I)                                               FL54470
         ASYV(I) = 0.0                                                   FL54480
C                                                                        FL54490
C        WRITE (IPR,300) I,AWCCON(I)                                     FL54500
C                                                                        FL54510
   90 CONTINUE                                                           FL54520
      RETURN                                                             FL54530
C                                                                        FL54540
C                                                                        FL54550
      END                                                                FL54560
      FUNCTION ABSLIM(TK,AWLWC)                                          FL54570
C                                                                        FL54580
C     CC                                                                 FL54590
C     CC    FOR CLOUD OR AEROSOL ATTENUATION AS FREQ APPROACHES ZERO     FL54600
C     CC    MODIFIED DEBYE EQUATIONS FROM RAY (1972) APPL. OPTICS VOL 11 FL54610
C     CC                                                                 FL54620
C     CC    ANO= 8.0*10**(-2)  (CM-4)                                    FL54630
C     CC    ALM= 41.*RR**(-0.21)  (CM-1)  RR IN (MM/HR)                  FL54640
C     CC                                                                 FL54650
C                                                                        FL54660
      DATA PI/3.14159265/                                                FL54670
C                                                                        FL54680
C     CC   ANO=0.08                                                      FL54690
C     CC   ALM=41./RR**0.21                                              FL54700
C                                                                        FL54710
      TC = TK-273.15                                                     FL54720
C                                                                        FL54730
C     CC                                                                 FL54740
C                                                                        FL54750
      EFIN = 5.27137+.0216474*TC-.00131198*TC*TC                         FL54760
      ES = 78.54*(1.-4.579E-03*(TC-25.)+1.19E-05*(TC-25.)**2-2.8E-08*(TC FL54770
     *   -25.)**3)                                                       FL54780
      SLAMBD = 3.3836E-04*EXP(2513.98/TK)                                FL54790
C                                                                        FL54800
C     CC                                                                 FL54810
C     CC   VOL=PI*ANO*ALM**(-4)                                          FL54820
C                                                                        FL54830
      ESMIE2 = (ES-EFIN)/(ES+2.0)**2                                     FL54840
C                                                                        FL54850
C     CC                                                                 FL54860
C     CC    DIVIDE VOLUME EQUIVALENT LIQUID BY 10 FOR UNITS CONVERSION   FL54870
C     CC                                                                 FL54880
C                                                                        FL54890
      EQLWC = AWLWC/10.0                                                 FL54900
C                                                                        FL54910
C     CC                                                                 FL54920
C                                                                        FL54930
      ABSLIM = 0.6951*TK*36.0*PI*EQLWC*SLAMBD*ESMIE2                     FL54940
C                                                                        FL54950
C     CC                                                                 FL54960
C                                                                        FL54970
      RETURN                                                             FL54980
      END                                                                FL54990
      BLOCK DATA TITLE                                                   FL55000
C                                                                        FL55010
C     >    BLOCK DATA                                                    FL55020
C     TITLE INFORMATION                                                  FL55030
C                                                                        FL55040
      CHARACTER*20 HHAZE,HSEASN,HVULCN,HMET,HMODEL,BLANK                 FL55050
      CHARACTER*24 HTRRAD                                                FL55060
      COMMON /TITL/ HHAZE(16),HSEASN(2),HVULCN(8),BLANK,                 FL55070
     * HMET(2),HMODEL(8),HTRRAD(4)                                       FL55080
      COMMON /VSBD/ VSB(10)                                              FL55090
      DATA VSB /23.,5.,0.,23.,5.,50.,23.,0.2,0.5,0./                     FL55100
      DATA BLANK/'                    '/                                 FL55110
      DATA HHAZE /                                                       FL55120
     * 'RURAL               ',                                           FL55130
     * 'RURAL               ',                                           FL55140
     * 'NAVY MARITIME       ',                                           FL55150
     * 'MARITIME            ',                                           FL55160
     * 'URBAN               ',                                           FL55170
     * 'TROPOSPHERIC        ',                                           FL55180
     * 'USER DEFINED        ',                                           FL55190
     * 'FOG1 (ADVECTTION)   ',                                           FL55200
     * 'FOG2 (RADIATION)    ',                                           FL55210
     * 'DESERT AEROSOL      ',                                           FL55220
     * 'BACKGROUND STRATO   ',                                           FL55230
     * 'AGED VOLCANIC       ',                                           FL55240
     * 'FRESH VOLCANIC      ',                                           FL55250
     * 'AGED VOLCANIC       ',                                           FL55260
     * 'FRESH VOCANIC       ',                                           FL55270
     * 'METEORIC DUST       '/                                           FL55280
      DATA HSEASN /                                                      FL55290
     * 'SPRING-SUMMER       ',                                           FL55300
     * 'FALL-WINTER         '/                                           FL55310
      DATA HVULCN /                                                      FL55320
     * 'BACKGROUND STRATO   ',                                           FL55330
     * 'MODERATE VOLCANIC   ',                                           FL55340
     * 'HIGH VOLCANIC       ',                                           FL55350
     * 'HIGH VOLCANIC       ',                                           FL55360
     * 'MODERATE VOLCANIC   ',                                           FL55370
     * 'MODERATE VOLCANIC   ',                                           FL55380
     * 'HIGH VOLCANIC       ',                                           FL55390
     * 'EXTREME VOLCANIC    '/                                           FL55400
      DATA HMET/                                                         FL55410
     * 'NORMAL              ',                                           FL55420
     * 'TRANSITION          '/                                           FL55430
      DATA HMODEL /                                                      FL55440
     * 'TROPICAL MODEL      ',                                           FL55450
     * 'MIDLATITUDE SUMMER  ',                                           FL55460
     * 'MIDLATITUDE WINTER  ',                                           FL55470
     * 'SUBARCTIC   SUMMER  ',                                           FL55480
     * 'SUBARCTIC   WINTER  ',                                           FL55490
     * '1976 U S STANDARD   ',                                           FL55500
     * '                    ',                                           FL55510
     * 'MODEL=0 HORIZONTAL  '/                                           FL55520
      DATA HTRRAD/                                                       FL55530
     * 'TRANSMITTANCE           ',                                       FL55540
     * 'RADIANCE                ',                                       FL55550
     * 'RADIANCE+SOLAR SCATTERNG',                                       FL55560
     * 'TRANSMITTED SOLAR IRRAD.'/                                       FL55570
      END                                                                FL55580
      BLOCK DATA PRFDTA                                                  FL55590
C                                                                        FL55600
C     >    BLOCK DATA                                                    FL55610
C                                                                        FL55620
C     AEROSOL PROFILE DATA                                               FL55630
C                                                                        FL55640
C     CC         0-2KM                                                   FL55650
C     CC           HZ2K=5 VIS PROFILES- 50KM,23KM,10KM,5KM,2KM           FL55660
C     CC         >2-10KM                                                 FL55670
C     CC           FAWI50=FALL/WINTER   50KM VIS                         FL55680
C     CC           FAWI23=FALL/WINTER    23KM VIS                        FL55690
C     CC           SPSU50=SPRING/SUMMER  50KM VIS                        FL55700
C     CC           SPSU23=SPRING/SUMMER  23KM VIS                        FL55710
C     CC         >10-30KM                                                FL55720
C     CC           BASTFW=BACKGROUND STRATOSPHERIC   FALL/WINTER         FL55730
C     CC           VUMOFW=MODERATE VOLCANIC          FALL/WINTER         FL55740
C     CC           HIVUFW=HIGH VOLCANIC              FALL/WINTER         FL55750
C     CC           EXVUFW=EXTREME VOLCANIC           FALL/WINTER         FL55760
C     CC           BASTSS,VUMOSS,HIVUSS,EXVUSS=      SPRING/SUMMER       FL55770
C     CC         >30-100KM                                               FL55780
C     CC           UPNATM=NORMAL UPPER ATMOSPHERIC                       FL55790
C     CC           VUTONO=TRANSITION FROM VOLCANIC TO NORMAL             FL55800
C     CC           VUTOEX=TRANSITION FROM VOLCANIC TO EXTREME            FL55810
C     CC           EXUPAT=EXTREME UPPER ATMOSPHERIC                      FL55820
C                                                                        FL55830
      COMMON/PRFD  /ZHT(34),HZ2K(34,5),FAWI50(34),FAWI23(34),SPSU50(34), FL55840
     *SPSU23(34),BASTFW(34),VUMOFW(34),HIVUFW(34),EXVUFW(34),BASTSS(34), FL55850
     *VUMOSS(34),HIVUSS(34),EXVUSS(34),UPNATM(34),VUTONO(34),            FL55860
     *VUTOEX(34),EXUPAT(34)                                              FL55870
      DATA ZHT/                                                          FL55880
     *    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,    FL55890
     *    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,    FL55900
     *   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,    FL55910
     *   35.,   40.,   45.,   50.,   70.,  100.,99999./                  FL55920
       DATA HZ2K(1,1),HZ2K(1,2),HZ2K(1,3),HZ2K(1,4),HZ2K(1,5)/           FL55930
     * 6.62E-02, 1.58E-01, 3.79E-01, 7.70E-01, 1.94E+00/                 FL55940
       DATA HZ2K(2,1),HZ2K(2,2),HZ2K(2,3),HZ2K(2,4),HZ2K(2,5)/           FL55950
     * 4.15E-02, 9.91E-02, 3.79E-01, 7.70E-01, 1.94E+00/                 FL55960
       DATA HZ2K(3,1),HZ2K(3,2),HZ2K(3,3),HZ2K(3,4),HZ2K(3,5)/           FL55970
     * 2.60E-02, 6.21E-02, 6.21E-02, 6.21E-02, 6.21E-02/                 FL55980
      DATA FAWI50  /3*0.,                                                FL55990
     * 1.14E-02, 6.43E-03, 4.85E-03, 3.54E-03, 2.31E-03, 1.41E-03,       FL56000
     * 9.80E-04,7.87E-04,23*0./                                          FL56010
      DATA FAWI23              /3*0.,                                    FL56020
     * 2.72E-02, 1.20E-02, 4.85E-03, 3.54E-03, 2.31E-03, 1.41E-03,       FL56030
     * 9.80E-04,7.87E-04, 23*0./                                         FL56040
      DATA  SPSU50              / 3*0.,                                  FL56050
     * 1.46E-02, 1.02E-02, 9.31E-03, 7.71E-03, 6.23E-03, 3.37E-03,       FL56060
     * 1.82E-03  ,1.14E-03,23*0./                                        FL56070
      DATA  SPSU23              / 3*0.,                                  FL56080
     * 3.46E-02, 1.85E-02, 9.31E-03, 7.71E-03, 6.23E-03, 3.37E-03,       FL56090
     * 1.82E-03  ,1.14E-03,23*0./                                        FL56100
      DATA BASTFW       /11*0.,                                          FL56110
     *           7.14E-04, 6.64E-04, 6.23E-04, 6.45E-04, 6.43E-04,       FL56120
     * 6.41E-04, 6.00E-04, 5.62E-04, 4.91E-04, 4.23E-04, 3.52E-04,       FL56130
     * 2.95E-04, 2.42E-04, 1.90E-04, 1.50E-04, 3.32E-05 ,7*0./           FL56140
      DATA    VUMOFW       /11*0.,                                       FL56150
     *           1.79E-03, 2.21E-03, 2.75E-03, 2.89E-03, 2.92E-03,       FL56160
     * 2.73E-03, 2.46E-03, 2.10E-03, 1.71E-03, 1.35E-03, 1.09E-03,       FL56170
     * 8.60E-04, 6.60E-04, 5.15E-04, 4.09E-04, 7.60E-05 ,7*0./           FL56180
      DATA    HIVUFW       /11*0.,                                       FL56190
     *           2.31E-03, 3.25E-03, 4.52E-03, 6.40E-03, 7.81E-03,       FL56200
     * 9.42E-03, 1.07E-02, 1.10E-02, 8.60E-03, 5.10E-03, 2.70E-03,       FL56210
     * 1.46E-03, 8.90E-04, 5.80E-04, 4.09E-04, 7.60E-05 ,7*0./           FL56220
      DATA    EXVUFW       /11*0.,                                       FL56230
     *           2.31E-03, 3.25E-03, 4.52E-03, 6.40E-03, 1.01E-02,       FL56240
     * 2.35E-02, 6.10E-02, 1.00E-01, 4.00E-02, 9.15E-03, 3.13E-03,       FL56250
     * 1.46E-03, 8.90E-04, 5.80E-04, 4.09E-04, 7.60E-05 ,7*0./           FL56260
      DATA    BASTSS       /11*0.,                                       FL56270
     *           7.99E-04, 6.41E-04, 5.17E-04, 4.42E-04, 3.95E-04,       FL56280
     * 3.82E-04, 4.25E-04, 5.20E-04, 5.81E-04, 5.89E-04, 5.02E-04,       FL56290
     * 4.20E-04, 3.00E-04, 1.98E-04, 1.31E-04, 3.32E-05 ,7*0./           FL56300
      DATA    VUMOSS       /11*0.,                                       FL56310
     *           2.12E-03, 2.45E-03, 2.80E-03, 2.89E-03, 2.92E-03,       FL56320
     * 2.73E-03, 2.46E-03, 2.10E-03, 1.71E-03, 1.35E-03, 1.09E-03,       FL56330
     * 8.60E-04, 6.60E-04, 5.15E-04, 4.09E-04, 7.60E-05 ,7*0./           FL56340
      DATA    HIVUSS       /11*0.,                                       FL56350
     *           2.12E-03, 2.45E-03, 2.80E-03, 3.60E-03, 5.23E-03,       FL56360
     * 8.11E-03, 1.20E-02, 1.52E-02, 1.53E-02, 1.17E-02, 7.09E-03,       FL56370
     * 4.50E-03, 2.40E-03, 1.28E-03, 7.76E-04, 7.60E-05 ,7*0./           FL56380
      DATA    EXVUSS       /11*0.,                                       FL56390
     *           2.12E-03, 2.45E-03, 2.80E-03, 3.60E-03, 5.23E-03,       FL56400
     * 8.11E-03, 1.27E-02, 2.32E-02, 4.85E-02, 1.00E-01, 5.50E-02,       FL56410
     * 6.10E-03, 2.40E-03, 1.28E-03, 7.76E-04, 7.60E-05 ,7*0./           FL56420
      DATA UPNATM       /26*0.,                                          FL56430
     * 3.32E-05, 1.64E-05, 7.99E-06, 4.01E-06, 2.10E-06, 1.60E-07,       FL56440
     * 9.31E-10, 0.      /                                               FL56450
      DATA VUTONO       /26*0.,                                          FL56460
     * 7.60E-05, 2.45E-05, 7.99E-06, 4.01E-06, 2.10E-06, 1.60E-07,       FL56470
     * 9.31E-10, 0.      /                                               FL56480
      DATA VUTOEX       /26*0.,                                          FL56490
     * 7.60E-05, 7.20E-05, 6.95E-05, 6.60E-05, 5.04E-05, 1.03E-05,       FL56500
     * 4.50E-07, 0.      /                                               FL56510
      DATA EXUPAT       /26*0.,                                          FL56520
     * 3.32E-05, 4.25E-05, 5.59E-05, 6.60E-05, 5.04E-05, 1.03E-05,       FL56530
     * 4.50E-07, 0.      /                                               FL56540
      END                                                                FL56550
      BLOCK DATA EXTDTA                                                  FL56560
C                                                                        FL56570
C     >    BLOCK DATA                                                    FL56580
C     CC                                                                 FL56590
C     CC   ALTITUDE REGIONS FOR AEROSOL EXTINCTION COEFFICIENTS          FL56600
C     CC                                                                 FL56610
C     CC                                                                 FL56620
C     CC         0-2KM                                                   FL56630
C     CC           RUREXT=RURAL EXTINCTION   RURABS=RURAL ABSORPTION     FL56640
C     CC           RURSYM=RURAL ASYMMETRY FACTORS                        FL56650
C     CC           URBEXT=URBAN EXTINCTION   URBABS=URBAN ABSORPTION     FL56660
C     CC           URBSYM=URBAN ASYMMETRY FACTORS                        FL56670
C     CC           OCNEXT=MARITIME EXTINCTION  OCNABS=MARITIME ABSORPTIO FL56680
C     CC           OCNSYM=MARITIME ASYMMETRY FACTORS                     FL56690
C     CC           TROEXT=TROPSPHER EXTINCTION  TROABS=TROPOSPHER ABSORP FL56700
C     CC           TROSYM=TROPSPHERIC ASYMMETRY FACTORS                  FL56710
C     CC           FG1EXT=FOG1 .2KM VIS EXTINCTION  FG1ABS=FOG1 ABSORPTI FL56720
C     CC           FG1SYM=FOG1 ASYMMETRY FACTORS                         FL56730
C     CC           FG2EXT=FOG2 .5KM VIS EXTINCTION  FG2ABS=FOG2 ABSORPTI FL56740
C     CC           FG2SYM=FOG2 ASYMMETRY FACTORS                         FL56750
C     CC         >2-10KM                                                 FL56760
C     CC           TROEXT=TROPOSPHER EXTINCTION  TROABS=TROPOSPHER ABSOR FL56770
C     CC           TROSYM=TROPOSPHERIC ASYMMETRY FACTORS                 FL56780
C     CC         >10-30KM                                                FL56790
C     CC           BSTEXT=BACKGROUND STRATOSPHERIC EXTINCTION            FL56800
C     CC           BSTABS=BACKGROUND STRATOSPHERIC ABSORPTION            FL56810
C     CC           BSTSYM=BACKGROUND STRATOSPHERIC ASYMMETRY FACTORS     FL56820
C     CC           AVOEXT=AGED VOLCANIC EXTINCTION                       FL56830
C     CC           AVOABS=AGED VOLCANIC ABSORPTION                       FL56840
C     CC           AVOSYM=AGED VOLCANIC ASYMMETRY FACTORS                FL56850
C     CC           FVOEXT=FRESH VOLCANIC EXTINCTION                      FL56860
C     CC           FVOABS=FRESH VOLCANIC ABSORPTION                      FL56870
C     CC           FVOSYM=FRESH VOLCANIC ASYMMETRY FACTORS               FL56880
C     CC         >30-100KM                                               FL56890
C     CC           DMEEXT=METEORIC DUST EXTINCTION                       FL56900
C     CC           DMEABS=METEORIC DUST ABSORPTION                       FL56910
C     CC           DMESYM=METEORIC DUST ASYMMETRY FACTORS                FL56920
C                                                                        FL56930
C     AEROSOL EXTINCTION AND ABSORPTION DATA                             FL56940
C                                                                        FL56950
C     MODIFIED TO INCLUDE ASYMMETRY DATA - JAN 1986 (A.E.R. INC.)        FL56960
C                                                                        FL56970
C     COMMON /EXTD  /VX2(40),RUREXT(40,4),RURABS(40,4),URBEXT(40,4),     FL56980
C     1URBABS(40,4),OCNEXT(40,4),OCNABS(40,4),TROEXT(40,4),TROABS(40,4), FL56990
C     2FG1EXT(40),FG1ABS(40),FG2EXT(40),FG2ABS(40),                      FL57000
C     3   BSTEXT(40),BSTABS(40),AVOEXT(40),AVOABS(40),FVOEXT(40)         FL57010
C     4),FVOABS(40),DMEEXT(40),DMEABS(40)                                FL57020
C                                                                        FL57030
      COMMON /EXTD  / VX2(47),RURE1(47),RURE2(47),RURE3(47),RURE4(47),   FL57040
     * RURA1(47),RURA2(47),RURA3(47),RURA4(47),                          FL57050
     * RURG1(47),RURG2(47),RURG3(47),RURG4(47),                          FL57060
     * URBE1(47),URBE2(47),URBE3(47),URBE4(47),                          FL57070
     * URBA1(47),URBA2(47),URBA3(47),URBA4(47),                          FL57080
     * URBG1(47),URBG2(47),URBG3(47),URBG4(47),                          FL57090
     * OCNE1(47),OCNE2(47),OCNE3(47),OCNE4(47),                          FL57100
     * OCNA1(47),OCNA2(47),OCNA3(47),OCNA4(47),                          FL57110
     * OCNG1(47),OCNG2(47),OCNG3(47),OCNG4(47),                          FL57120
     * TROE1(47),TROE2(47),TROE3(47),TROE4(47),                          FL57130
     * TROA1(47),TROA2(47),TROA3(47),TROA4(47),                          FL57140
     * TROG1(47),TROG2(47),TROG3(47),TROG4(47),                          FL57150
     * FG1EXT(47),FG1ABS(47),FG1SYM(47),FG2EXT(47),FG2ABS(47),           FL57160
     * FG2SYM(47),BSTEXT(47),BSTABS(47),BSTSYM(47),AVOEXT(47),           FL57170
     * AVOABS(47),AVOSYM(47),FVOEXT(47),FVOABS(47),FVOSYM(47),           FL57180
     * DMEEXT(47),DMEABS(47),DMESYM(47),CCUEXT(47),CCUABS(47),           FL57190
     * CCUSYM(47),CALEXT(47),CALABS(47),CALSYM(47),CSTEXT(47),           FL57200
     * CSTABS(47),CSTSYM(47),CSCEXT(47),CSCABS(47),CSCSYM(47),           FL57210
     * CNIEXT(47),CNIABS(47),CNISYM(47)                                  FL57220
C                                                                        FL57230
C     CI64--    STANDARD  CIRRUS  CLOUD  MODEL                           FL57240
C     ICE 64 MICRON MODE RADIUS CIRRUS CLOUD MODEL                       FL57250
C                                                                        FL57260
C     CIR4--    OPTICALLY  THIN  CIRRUS  MODEL                           FL57270
C     ICE  4 MICRON MODE RADIUS CIRRUS CLOUD MODEL                       FL57280
C                                                                        FL57290
       COMMON/CIRR/ CI64XT(47),CI64AB(47),CI64G(47),                     FL57300
     *              CIR4XT(47),CIR4AB(47),CIR4G(47)                      FL57310
      DATA VX2 /                                                         FL57320
     *   .2000,   .3000,   .3371,   .5500,   .6943,  1.0600,  1.5360,    FL57330
     *  2.0000,  2.2500,  2.5000,  2.7000,  3.0000,  3.3923,  3.7500,    FL57340
     *  4.5000,  5.0000,  5.5000,  6.0000,  6.2000,  6.5000,  7.2000,    FL57350
     *  7.9000,  8.2000,  8.7000,  9.0000,  9.2000, 10.0000, 10.5910,    FL57360
     * 11.0000, 11.5000, 12.5000, 14.8000, 15.0000, 16.4000, 17.2000,    FL57370
     * 18.5000, 21.3000, 25.0000, 30.0000, 40.0000, 50.0000, 60.0000,    FL57380
     * 80.0000, 100.000, 150.000, 200.000, 300.000/                      FL57390
      DATA RURE1 /                                                       FL57400
     * 2.09291, 1.74582, 1.60500, 1.00000,  .75203,  .41943,  .24070,    FL57410
     *  .14709,  .13304,  .12234,  .13247,  .11196,  .10437,  .09956,    FL57420
     *  .09190,  .08449,  .07861,  .07025,  .07089,  .07196,  .07791,    FL57430
     *  .04481,  .04399,  .12184,  .12658,  .12829,  .09152,  .08076,    FL57440
     *  .07456,  .06880,  .06032,  .04949,  .05854,  .06000,  .06962,    FL57450
     *  .05722,  .06051,  .05177,  .04589,  .04304,                      FL57460
     *  .03582,  .03155,  .02018,  .01469,  .00798,  .00551, 0./         FL57470
      DATA RURE2 /                                                       FL57480
     * 2.09544, 1.74165, 1.59981, 1.00000,  .75316,  .42171,  .24323,    FL57490
     *  .15108,  .13608,  .12430,  .13222,  .13823,  .11076,  .10323,    FL57500
     *  .09475,  .08728,  .08076,  .07639,  .07797,  .07576,  .07943,    FL57510
     *  .04899,  .04525,  .12165,  .12741,  .12778,  .09032,  .07962,    FL57520
     *  .07380,  .06880,  .06329,  .05791,  .06646,  .06639,  .07443,    FL57530
     *  .06304,  .06443,  .05538,  .04867,  .04519,                      FL57540
     *  .03821,  .03374,  .02173,  .01587,  .00862,  .00594, 0./         FL57550
      DATA RURE3 /                                                       FL57560
     * 2.07082, 1.71456, 1.57962, 1.00000,  .76095,  .43228,  .25348,    FL57570
     *  .16456,  .14677,  .13234,  .13405,  .20316,  .12873,  .11506,    FL57580
     *  .10481,  .09709,  .08918,  .09380,  .09709,  .08791,  .08601,    FL57590
     *  .06247,  .05601,  .11905,  .12595,  .12348,  .08741,  .07703,    FL57600
     *  .07266,  .07044,  .07443,  .08146,  .08810,  .08563,  .08962,    FL57610
     *  .08051,  .07677,  .06658,  .05747,  .05184,                      FL57620
     *  .04572,  .04074,  .02689,  .01981,  .01084,  .00714, 0./         FL57630
      DATA RURE4 /                                                       FL57640
     * 1.66076, 1.47886, 1.40139, 1.00000,  .80652,  .50595,  .32259,    FL57650
     *  .23468,  .20772,  .18532,  .17348,  .35114,  .20006,  .17386,    FL57660
     *  .16139,  .15424,  .14557,  .16215,  .16766,  .14994,  .14032,    FL57670
     *  .12968,  .12601,  .13551,  .13582,  .13228,  .11070,  .09994,    FL57680
     *  .09873,  .10418,  .13241,  .15924,  .16139,  .15949,  .15778,    FL57690
     *  .15184,  .13848,  .12563,  .11076,  .09601,                      FL57700
     *  .09312,  .08720,  .06644,  .05264,  .03181,  .02196, 0.0/        FL57710
      DATA RURA1 /                                                       FL57720
     *  .67196,  .11937,  .08506,  .05930,  .05152,  .05816,  .05006,    FL57730
     *  .01968,  .02070,  .02101,  .05652,  .02785,  .01316,  .00867,    FL57740
     *  .01462,  .01310,  .01627,  .02013,  .02165,  .02367,  .03538,    FL57750
     *  .02823,  .03962,  .06778,  .07285,  .08120,  .04032,  .03177,    FL57760
     *  .02557,  .02342,  .02177,  .02627,  .03943,  .03114,  .03696,    FL57770
     *  .02956,  .03500,  .03241,  .03297,  .03380,                      FL57780
     *  .03170,  .02794,  .01769,  .01305,  .00730,  .00518, 0.0/        FL57790
      DATA RURA2 /                                                       FL57800
     *  .62968,  .10816,  .07671,  .05380,  .04684,  .05335,  .04614,    FL57810
     *  .01829,  .01899,  .01962,  .05525,  .06816,  .01652,  .00867,    FL57820
     *  .01544,  .01373,  .01627,  .02892,  .02829,  .02532,  .03487,    FL57830
     *  .02835,  .03854,  .06684,  .07272,  .08038,  .03987,  .03247,    FL57840
     *  .02816,  .02816,  .03101,  .03741,  .04829,  .04032,  .04399,    FL57850
     *  .03734,  .03956,  .03601,  .03525,  .03563,                      FL57860
     * .03357,  .02965,  .01887,  .01395,  .00782,  .00555, 0.0/         FL57870
      DATA RURA3 /                                                       FL57880
     *  .51899,  .08278,  .05816,  .04082,  .03570,  .04158,  .03620,    FL57890
     *  .01513,  .01481,  .01633,  .05278,  .13690,  .02494,  .00886,    FL57900
     *  .01804,  .01582,  .01677,  .04816,  .04367,  .03013,  .03443,    FL57910
     *  .02930,  .03677,  .06209,  .06911,  .07475,  .03892,  .03494,    FL57920
     *  .03513,  .03968,  .05152,  .06241,  .06937,  .06203,  .06215,    FL57930
     *  .05614,  .05209,  .04608,  .04196,  .04095,                      FL57940
     *  .03916,  .03486,  .02262,  .01686,  .00951,  .00674, 0.0/        FL57950
      DATA RURA4 /                                                       FL57960
     *  .21943,  .02848,  .01943,  .01342,  .01171,  .01437,  .01323,    FL57970
     *  .01152,  .00696,  .01329,  .06108,  .24690,  .05323,  .01430,    FL57980
     *  .03361,  .02949,  .02652,  .09437,  .08506,  .05348,  .04627,    FL57990
     *  .04380,  .04557,  .05380,  .05715,  .05899,  .04861,  .05253,    FL58000
     *  .06171,  .07437,  .10152,  .12019,  .12190,  .11734,  .11411,    FL58010
     *  .10766,  .09487,  .08430,  .07348,  .06861,                      FL58020
     *  .06936,  .06458,  .04735,  .03761,  .02313,  .01668, 0.0/        FL58030
      DATA RURG1 /                                                       FL58040
     *  .7581,   .6785,   .6712,   .6479,   .6342,   .6176,   .6334,     FL58050
     *  .7063,   .7271,   .7463,   .7788,   .7707,   .7424,   .7312,     FL58060
     *  .7442,   .7516,   .7662,   .7940,   .7886,   .7797,   .7664,     FL58070
     *  .8525,   .8700,   .5846,   .5570,   .5992,   .6159,   .6271,     FL58080
     *  .6257,   .6374,   .6546,   .6861,   .6859,   .6120,   .5570,     FL58090
     *  .5813,   .5341,   .5284,   .5137,   .4348,   .4223,   .3775,     FL58100
     *  .3435,   .3182,   .2791,   .2494,   .0000/                       FL58110
      DATA RURG2 /                                                       FL58120
     *  .7632,   .6928,   .6865,   .6638,   .6498,   .6314,   .6440,     FL58130
     *  .7098,   .7303,   .7522,   .7903,   .7804,   .7380,   .7319,     FL58140
     *  .7508,   .7584,   .7738,   .8071,   .7929,   .7843,   .7747,     FL58150
     *  .8507,   .8750,   .6112,   .5851,   .6272,   .6466,   .6616,     FL58160
     *  .6653,   .6798,   .6965,   .7026,   .6960,   .6360,   .5848,     FL58170
     *  .6033,   .5547,   .5445,   .5274,   .4518,   .4318,   .3863,     FL58180
     *  .3516,   .3257,   .2853,   .2548,   .0000/                       FL58190
      DATA RURG3 /                                                       FL58200
     *  .7725,   .7240,   .7197,   .6997,   .6858,   .6650,   .6702,     FL58210
     *  .7181,   .7378,   .7653,   .8168,   .7661,   .7286,   .7336,     FL58220
     *  .7654,   .7735,   .7910,   .8303,   .8025,   .7957,   .7946,     FL58230
     *  .8468,   .8734,   .6831,   .6619,   .6994,   .7250,   .7449,     FL58240
     *  .7547,   .7665,   .7644,   .7265,   .7170,   .6769,   .6409,     FL58250
     *  .6442,   .6031,   .5854,   .5646,   .4977,   .4602,   .4127,     FL58260
     *  .3751,   .3476,   .3048,   .2721,   .0000/                       FL58270
      DATA RURG4 /                                                       FL58280
     *  .7778,   .7793,   .7786,   .7717,   .7628,   .7444,   .7365,     FL58290
     *  .7491,   .7609,   .7921,   .8688,   .7537,   .7294,   .7413,     FL58300
     *  .7928,   .8016,   .8225,   .8761,   .8359,   .8285,   .8385,     FL58310
     *  .8559,   .8654,   .8414,   .8415,   .8527,   .8740,   .8903,     FL58320
     *  .8952,   .8923,   .8611,   .8033,   .7989,   .7758,   .7632,     FL58330
     *  .7508,   .7314,   .7091,   .6867,   .6419,   .5790,   .5259,     FL58340
     *  .4749,   .4415,   .3886,   .3489,   .0000/                       FL58350
      DATA URBE1 /                                                       FL58360
     * 1.88816, 1.63316, 1.51867, 1.00000,  .77785,  .47095,  .30006,    FL58370
     *  .21392,  .19405,  .17886,  .18127,  .16133,  .14785,  .14000,    FL58380
     *  .12715,  .11880,  .11234,  .10601,  .10500,  .10361,  .10342,    FL58390
     *  .08766,  .08652,  .11937,  .12139,  .12297,  .09797,  .09057,    FL58400
     *  .08595,  .08196,  .07563,  .06696,  .07209,  .06842,  .07177,    FL58410
     *  .06354,  .06177,  .05373,  .04728,  .04051,                      FL58420
     *  .03154,  .02771,  .01759,  .01278,  .00693,  .00480, 0.0/        FL58430
      DATA URBE2 /                                                       FL58440
     * 1.95582, 1.64994, 1.53070, 1.00000,  .77614,  .46639,  .29487,    FL58450
     *  .21051,  .18943,  .17285,  .17209,  .21418,  .15354,  .14051,    FL58460
     *  .12728,  .11861,  .11089,  .11329,  .11323,  .10563,  .10247,    FL58470
     *  .08696,  .08361,  .12013,  .12418,  .12304,  .09614,  .08842,    FL58480
     *  .08487,  .08285,  .08361,  .08430,  .08880,  .08449,  .08601,    FL58490
     *  .07835,  .07323,  .06367,  .05500,  .04747,                      FL58500
     *  .03901,  .03454,  .02240,  .01638,  .00891,  .00612, 0.0/        FL58510
      DATA URBE3 /                                                       FL58520
     * 1.96430, 1.64032, 1.52392, 1.00000,  .77709,  .46253,  .28690,    FL58530
     *  .20310,  .17981,  .16101,  .15614,  .26475,  .15456,  .13563,    FL58540
     *  .12215,  .11361,  .10500,  .11715,  .11753,  .10392,  .09766,    FL58550
     *  .08443,  .08057,  .10943,  .11342,  .11063,  .08703,  .08025,    FL58560
     *  .07886,  .08032,  .09101,  .10070,  .10386,  .09943,  .09886,    FL58570
     *  .09152,  .08247,  .07152,  .06089,  .05253,                      FL58580
     *  .04582,  .04091,  .02717,  .02008,  .01103,  .00754, 0.0/        FL58590
      DATA URBE4 /                                                       FL58600
     * 1.41266, 1.33816, 1.29114, 1.00000,  .83646,  .55025,  .35342,    FL58610
     *  .25285,  .21576,  .18310,  .16215,  .37854,  .20494,  .16665,    FL58620
     *  .14778,  .13892,  .12943,  .15525,  .15709,  .13513,  .12481,    FL58630
     *  .11759,  .11494,  .11487,  .11329,  .11108,  .09911,  .09209,    FL58640
     *  .09342,  .10120,  .13177,  .15696,  .15766,  .15513,  .15203,    FL58650
     *  .14532,  .13038,  .11785,  .10411,  .09101,                      FL58660
     *  .08907,  .08399,  .06579,  .05337,  .03372,  .02379, 0.0/        FL58670
      DATA URBA1 /                                                       FL58680
     *  .78437,  .58975,  .54285,  .36184,  .29222,  .20886,  .15658,    FL58690
     *  .12329,  .11462,  .10747,  .11797,  .10025,  .08759,  .08184,    FL58700
     *  .07506,  .07006,  .06741,  .06601,  .06544,  .06449,  .06665,    FL58710
     *  .06278,  .06949,  .07316,  .07462,  .08101,  .05753,  .05272,    FL58720
     *  .04899,  .04734,  .04494,  .04443,  .05133,  .04348,  .04443,    FL58730
     *  .03994,  .03981,  .03633,  .03468,  .03146,                      FL58740
     *  .02809,  .02471,  .01556,  .01145,  .00639,  .00454, 0.0/        FL58750
      DATA URBA2 /                                                       FL58760
     *  .69032,  .49367,  .45165,  .29741,  .24070,  .17399,  .13146,    FL58770
     *  .10354,  .09589,  .09025,  .10411,  .15101,  .07880,  .06949,    FL58780
     *  .06570,  .06095,  .05829,  .07171,  .06797,  .05975,  .06013,    FL58790
     *  .05589,  .06051,  .07139,  .07494,  .07956,  .05525,  .05184,    FL58800
     *  .05089,  .05291,  .05886,  .06380,  .06880,  .06127,  .06019,    FL58810
     *  .05525,  .05070,  .04500,  .04076,  .03741,                      FL58820
     *  .03400,  .03010,  .01926,  .01427,  .00800,  .00567, 0.0/        FL58830
      DATA URBA3 /                                                       FL58840
     *  .54848,  .37101,  .33734,  .21949,  .17785,  .12968,  .09854,    FL58850
     *  .07804,  .07165,  .06791,  .08563,  .19639,  .06722,  .05316,    FL58860
     *  .05316,  .04886,  .04620,  .07570,  .06899,  .05291,  .05101,    FL58870
     *  .04734,  .05025,  .06171,  .06570,  .06854,  .04892,  .04797,    FL58880
     *  .05057,  .05665,  .07127,  .08095,  .08411,  .07728,  .07475,    FL58890
     *  .06886,  .06019,  .05222,  .04538,  .04171,                      FL58900
     *  .03911,  .03486,  .02271,  .01697,  .00961,  .00681, 0.0/        FL58910
      DATA URBA4 /                                                       FL58920
     *  .15975,  .10000,  .09013,  .05785,  .04671,  .03424,  .02633,    FL58930
     *  .02525,  .01975,  .02354,  .06241,  .26690,  .05810,  .02285,    FL58940
     *  .03810,  .03386,  .03044,  .09627,  .08557,  .05405,  .04576,    FL58950
     *  .04392,  .04424,  .04671,  .04791,  .04861,  .04684,  .05177,    FL58960
     *  .06158,  .07475,  .10342,  .12146,  .12177,  .11734,  .11335,    FL58970
     *  .10608,  .09171,  .08063,  .06968,  .06475,                      FL58980
     *  .06559,  .06131,  .04591,  .03714,  .02365,  .01734, 0.0/        FL58990
      DATA URBG1 /                                                       FL59000
     *  .7785,   .7182,   .7067,   .6617,   .6413,   .6166,   .6287,     FL59010
     *  .6883,   .7070,   .7243,   .7370,   .7446,   .7391,   .7371,     FL59020
     *  .7414,   .7435,   .7466,   .7543,   .7498,   .7424,   .7270,     FL59030
     *  .7674,   .7850,   .5880,   .5616,   .5901,   .6159,   .6238,     FL59040
     *  .6240,   .6281,   .6306,   .6298,   .6252,   .5785,   .5378,     FL59050
     *  .5512,   .5072,   .4930,   .4709,   .4009,   .4110,   .3672,     FL59060
     *  .3344,   .3093,   .2717,   .2426,   .0000/                       FL59070
      DATA URBG2 /                                                       FL59080
     *  .7906,   .7476,   .7385,   .6998,   .6803,   .6536,   .6590,     FL59090
     *  .7066,   .7258,   .7484,   .7769,   .7405,   .7351,   .7459,     FL59100
     *  .7625,   .7673,   .7759,   .7910,   .7732,   .7703,   .7644,     FL59110
     *  .7966,   .8142,   .6635,   .6428,   .6700,   .6935,   .7050,     FL59120
     *  .7092,   .7145,   .7094,   .6762,   .6684,   .6316,   .5997,     FL59130
     *  .6013,   .5625,   .5433,   .5198,   .4552,   .4387,   .3928,     FL59140
     *  .3575,   .3310,   .2899,   .2588,   .0000/                       FL59150
      DATA URBG3 /                                                       FL59160
     *  .7949,   .7713,   .7650,   .7342,   .7162,   .6873,   .6820,     FL59170
     *  .7131,   .7312,   .7583,   .8030,   .7171,   .7185,   .7400,     FL59180
     *  .7698,   .7778,   .7923,   .8142,   .7864,   .7867,   .7891,     FL59190
     *  .8147,   .8298,   .7276,   .7136,   .7361,   .7590,   .7729,     FL59200
     *  .7783,   .7808,   .7624,   .7094,   .7022,   .6714,   .6480,     FL59210
     *  .6417,   .6104,   .5887,   .5651,   .5058,   .4692,   .4212,     FL59220
     *  .3825,   .3549,   .3112,   .2778,   .0000/                       FL59230
      DATA URBG4 /                                                       FL59240
     *  .7814,   .7993,   .7995,   .7948,   .7870,   .7682,   .7751,     FL59250
     *  .7501,   .7565,   .7809,   .8516,   .7137,   .7039,   .7241,     FL59260
     *  .7728,   .7846,   .8093,   .8576,   .8125,   .8140,   .8304,     FL59270
     *  .8472,   .8549,   .8525,   .8569,   .8640,   .8853,   .9017,     FL59280
     *  .9061,   .9021,   .8685,   .8126,   .8091,   .7897,   .7802,     FL59290
     *  .7691,   .7550,   .7353,   .7146,   .6754,   .6134,   .5601,     FL59300
     *  .5056,   .4701,   .4134,   .3714,   .0000/                       FL59310
      DATA OCNE1 /                                                       FL59320
     * 1.47576, 1.32614, 1.26171, 1.00000,  .88133,  .70297,  .56487,    FL59330
     *  .46006,  .42044,  .38310,  .35076,  .42266,  .32278,  .28810,    FL59340
     *  .24905,  .21184,  .16734,  .14791,  .21532,  .15076,  .12057,    FL59350
     *  .10038,  .10703,  .15070,  .15665,  .14639,  .10228,  .08367,    FL59360
     *  .07373,  .06829,  .05044,  .04373,  .04962,  .06158,  .07703,    FL59370
     *  .07234,  .06297,  .05481,  .05329,  .08741,                      FL59380
     *  .04608,  .03959,  .02382,  .01712,  .00936,  .00665, 0.0/        FL59390
      DATA OCNE2 /                                                       FL59400
     * 1.36924, 1.25443, 1.20835, 1.00000,  .91367,  .77089,  .64987,    FL59410
     *  .54886,  .50247,  .45038,  .38209,  .50589,  .43766,  .38076,    FL59420
     *  .31658,  .27475,  .22215,  .21019,  .27570,  .21057,  .16949,    FL59430
     *  .14209,  .14215,  .16956,  .17082,  .16025,  .11665,  .09759,    FL59440
     *  .09215,  .09373,  .10532,  .12570,  .13000,  .13633,  .14291,    FL59450
     *  .13506,  .11475,  .09658,  .08291,  .10348,                      FL59460
     *  .06693,  .05786,  .03522,  .02519,  .01358,  .00954, 0.0/        FL59470
      DATA OCNE3 /                                                       FL59480
     * 1.22259, 1.14627, 1.11842, 1.00000,  .94766,  .87538,  .80418,    FL59490
     *  .72930,  .68582,  .62165,  .49962,  .67949,  .66468,  .59253,    FL59500
     *  .49551,  .44671,  .37886,  .35924,  .43367,  .37019,  .30842,    FL59510
     *  .26437,  .25228,  .24905,  .23975,  .22766,  .17804,  .15316,    FL59520
     *  .15373,  .16791,  .22361,  .28348,  .28677,  .29082,  .29038,    FL59530
     *  .27810,  .23867,  .20209,  .16430,  .14943,                      FL59540
     *  .12693,  .11177,  .07095,  .05084,  .02690,  .01838, 0.0/        FL59550
      DATA OCNE4 /                                                       FL59560
     * 1.09133, 1.06601, 1.05620, 1.00000,  .97506,  .94791,  .94203,    FL59570
     *  .93671,  .92867,  .90411,  .80253,  .89222,  .94462,  .92146,    FL59580
     *  .85797,  .82595,  .76747,  .68646,  .78209,  .75266,  .68658,    FL59590
     *  .62722,  .60228,  .56335,  .53728,  .51861,  .43449,  .37196,    FL59600
     *  .35899,  .37316,  .46854,  .58234,  .58690,  .60348,  .60563,    FL59610
     *  .60000,  .55392,  .50367,  .43576,  .35949,                      FL59620
     *  .34729,  .32254,  .23600,  .17953,  .10071,  .06714, 0.0/        FL59630
      DATA OCNA1 /                                                       FL59640
     *  .30987,  .04354,  .02880,  .01797,  .01468,  .01766,  .01582,    FL59650
     *  .00816,  .01146,  .01677,  .03310,  .03380,  .00715,  .00443,    FL59660
     *  .00500,  .00601,  .00753,  .01595,  .02943,  .00994,  .01367,    FL59670
     *  .01671,  .02538,  .03481,  .03405,  .03601,  .01608,  .01310,    FL59680
     *  .01152,  .01082,  .01070,  .01563,  .02063,  .03171,  .03810,    FL59690
     *  .03741,  .03804,  .03759,  .04209,  .07892,                      FL59700
     *  .04347,  .03754,  .02269,  .01649,  .00917,  .00657, 0.0/        FL59710
      DATA OCNA2 /                                                       FL59720
     *  .23367,  .03127,  .02070,  .01297,  .01063,  .01285,  .01190,    FL59730
     *  .00937,  .00911,  .01576,  .05576,  .23487,  .03949,  .00905,    FL59740
     *  .02057,  .01816,  .01665,  .08025,  .08044,  .03677,  .03139,    FL59750
     *  .03190,  .03766,  .04532,  .04544,  .04715,  .03405,  .03614,    FL59760
     *  .04329,  .05424,  .07823,  .09728,  .10057,  .10247,  .10222,    FL59770
     *  .09551,  .08241,  .07158,  .06506,  .09203,                      FL59780
     *  .06133,  .05332,  .03258,  .02366,  .01308,  .00932, 0.0/        FL59790
      DATA OCNA3 /                                                       FL59800
     *  .13025,  .01557,  .01013,  .00646,  .00532,  .00665,  .00722,    FL59810
     *  .01335,  .00728,  .01810,  .09835,  .37329,  .09703,  .01968,    FL59820
     *  .05114,  .04342,  .03709,  .17456,  .16468,  .08785,  .06880,    FL59830
     *  .06589,  .06791,  .07247,  .07329,  .07449,  .07025,  .07962,    FL59840
     *  .09899,  .12481,  .17867,  .22019,  .22228,  .22051,  .21595,    FL59850
     *  .20335,  .17278,  .14677,  .12171,  .12430,                      FL59860
     *  .10890,  .09644,  .06106,  .04465,  .02457,  .01732, 0.0/        FL59870
      DATA OCNA4 /                                                       FL59880
     *  .03506,  .00323,  .00215,  .00139,  .00114,  .00171,  .00532,    FL59890
     *  .03082,  .01101,  .03741,  .20101,  .47608,  .21165,  .05234,    FL59900
     *  .12886,  .11215,  .09684,  .32810,  .31778,  .20513,  .16658,    FL59910
     *  .15956,  .15842,  .15905,  .15968,  .16051,  .16506,  .18323,    FL59920
     *  .21709,  .25652,  .33222,  .39639,  .39854,  .40297,  .40025,    FL59930
     *  .39025,  .35468,  .32006,  .27715,  .25348,                      FL59940
     *  .25632,  .23876,  .17092,  .13198,  .07692,  .05407, 0.0/        FL59950
      DATA OCNG1 /                                                       FL59960
     *  .7516,   .6960,   .6920,   .6756,   .6767,   .6844,   .6936,     FL59970
     *  .7055,   .7110,   .7177,   .7367,   .6287,   .6779,   .6784,     FL59980
     *  .6599,   .6659,   .6859,   .6887,   .6095,   .6558,   .6665,     FL59990
     *  .6697,   .6594,   .5851,   .5644,   .5760,   .5903,   .5991,     FL60000
     *  .6024,   .5979,   .6087,   .5837,   .5763,   .5348,   .4955,     FL60010
     *  .4821,   .4635,   .4373,   .3944,   .2344,   .2754,   .2447,     FL60020
     *  .2266,   .2088,   .1766,   .1481,   .0000/                       FL60030
      DATA OCNG2 /                                                       FL60040
     *  .7708,   .7288,   .7243,   .7214,   .7211,   .7330,   .7445,     FL60050
     *  .7579,   .7649,   .7790,   .8182,   .7673,   .7171,   .7205,     FL60060
     *  .7235,   .7251,   .7397,   .7537,   .6934,   .7137,   .7193,     FL60070
     *  .7206,   .7151,   .6732,   .6620,   .6696,   .6821,   .6895,     FL60080
     *  .6898,   .6819,   .6556,   .5925,   .5869,   .5511,   .5284,     FL60090
     *  .5124,   .4912,   .4646,   .4302,   .3124,   .3101,   .2752,     FL60100
     *  .2529,   .2335,   .2021,   .1738,   .0000/                       FL60110
      DATA OCNG3 /                                                       FL60120
     *  .7954,   .7782,   .7752,   .7717,   .7721,   .7777,   .7872,     FL60130
     *  .8013,   .8089,   .8301,   .8844,   .8332,   .7557,   .7597,     FL60140
     *  .7823,   .7822,   .7944,   .8157,   .7712,   .7738,   .7784,     FL60150
     *  .7807,   .7800,   .7682,   .7659,   .7692,   .7780,   .7828,     FL60160
     *  .7776,   .7621,   .7115,   .6342,   .6294,   .5999,   .5854,     FL60170
     *  .5700,   .5512,   .5265,   .4996,   .4236,   .3765,   .3357,     FL60180
     *  .3066,   .2830,   .2466,   .2184,   .0000/                       FL60190
      DATA OCNG4 /                                                       FL60200
     *  .8208,   .8270,   .8260,   .8196,   .8176,   .8096,   .8096,     FL60210
     *  .8202,   .8255,   .8520,   .9228,   .8950,   .7965,   .7847,     FL60220
     *  .8242,   .8244,   .8376,   .8857,   .8463,   .8332,   .8379,     FL60230
     *  .8441,   .8467,   .8502,   .8534,   .8562,   .8688,   .8789,     FL60240
     *  .8785,   .8683,   .8252,   .7562,   .7519,   .7261,   .7141,     FL60250
     *  .6980,   .6789,   .6540,   .6294,   .5783,   .5100,   .4595,     FL60260
     *  .4164,   .3868,   .3404,   .3042,   .0000/                       FL60270
      DATA TROE1 /                                                       FL60280
     * 2.21222, 1.82753, 1.67032, 1.00000,  .72424,  .35272,  .15234,    FL60290
     *  .05165,  .03861,  .02994,  .04671,  .02462,  .01538,  .01146,    FL60300
     *  .01032,  .00816,  .00861,  .00994,  .01057,  .01139,  .01747,    FL60310
     *  .01494,  .02418,  .03165,  .03386,  .04247,  .01601,  .01215,    FL60320
     *  .00937,  .00861,  .00823,  .01139,  .01924,  .01234,  .01348,    FL60330
     *  .01114,  .01297,  .01266,  .01418,  .01487,                      FL60340
     *  .01543,  .01321,  .00793,  .00582,  .00330,  .00239, 0.0/        FL60350
      DATA TROE2 /                                                       FL60360
     * 2.21519, 1.82266, 1.66557, 1.00000,  .72525,  .35481,  .15449,    FL60370
     *  .05475,  .04044,  .03082,  .04620,  .05272,  .01867,  .01266,    FL60380
     *  .01127,  .00886,  .00886,  .01449,  .01399,  .01228,  .01728,    FL60390
     *  .01475,  .02285,  .03215,  .03494,  .04285,  .01652,  .01304,    FL60400
     *  .01101,  .01120,  .01297,  .01753,  .02468,  .01741,  .01766,    FL60410
     *  .01513,  .01557,  .01456,  .01532,  .01582,                      FL60420
     *  .01619,  .01386,  .00832,  .00610,  .00346,  .00251, 0.0/        FL60430
      DATA TROE3 /                                                       FL60440
     * 2.19082, 1.79462, 1.64456, 1.00000,  .73297,  .36443,  .16278,    FL60450
     *  .06468,  .04658,  .03399,  .04538,  .11892,  .02835,  .01646,    FL60460
     *  .01386,  .01076,  .00968,  .02551,  .02222,  .01468,  .01690,    FL60470
     *  .01437,  .01994,  .03127,  .03513,  .04076,  .01722,  .01513,    FL60480
     *  .01519,  .01791,  .02538,  .03272,  .03816,  .03038,  .02886,    FL60490
     *  .02551,  .02228,  .01937,  .01804,  .01791,                      FL60500
     *  .01798,  .01539,  .00924,  .00678,  .00384,  .00278, 0.0/        FL60510
      DATA TROE4 /                                                       FL60520
     * 1.75696, 1.54829, 1.45962, 1.00000,  .77816,  .43139,  .21778,    FL60530
     *  .11329,  .08101,  .05506,  .04943,  .25291,  .06816,  .03703,    FL60540
     *  .02601,  .01968,  .01468,  .04962,  .04247,  .02234,  .01797,    FL60550
     *  .01532,  .01633,  .02259,  .02487,  .02595,  .01728,  .01892,    FL60560
     *  .02399,  .03247,  .05285,  .06462,  .06608,  .05930,  .05525,    FL60570
     *  .04861,  .03753,  .02968,  .02348,  .02165,                      FL60580
     *  .02152,  .01841,  .01104,  .00809,  .00458,  .00332, 0.0/        FL60590
      DATA TROA1 /                                                       FL60600
     *  .69671,  .09905,  .06563,  .04101,  .03354,  .03627,  .02810,    FL60610
     *  .00873,  .00918,  .00930,  .03215,  .01285,  .00513,  .00316,    FL60620
     *  .00557,  .00494,  .00646,  .00867,  .00937,  .01025,  .01646,    FL60630
     *  .01481,  .02418,  .02886,  .03070,  .04032,  .01494,  .01139,    FL60640
     *  .00873,  .00816,  .00797,  .01133,  .01911,  .01215,  .01329,    FL60650
     *  .01101,  .01291,  .01266,  .01418,  .01487,                      FL60660
     *  .01543,  .01321,  .00793,  .00582,  .00330,  .00239, 0.0/        FL60670
      DATA TROA2 /                                                       FL60680
     *  .65000,  .08791,  .05816,  .03652,  .02994,  .03278,  .02557,    FL60690
     *  .00810,  .00842,  .00867,  .03139,  .03949,  .00646,  .00316,    FL60700
     *  .00595,  .00519,  .00646,  .01304,  .01247,  .01095,  .01620,    FL60710
     *  .01449,  .02278,  .02930,  .03184,  .04063,  .01544,  .01234,    FL60720
     *  .01044,  .01076,  .01272,  .01741,  .02462,  .01722,  .01747,    FL60730
     *  .01506,  .01551,  .01456,  .01532,  .01582,                      FL60740
     *  .01619,  .01386,  .00832,  .00610,  .00346,  .00251, 0.0/        FL60750
      DATA TROA3 /                                                       FL60760
     *  .52804,  .06367,  .04158,  .02633,  .02184,  .02443,  .01937,    FL60770
     *  .00658,  .00646,  .00709,  .02949,  .10013,  .00968,  .00310,    FL60780
     *  .00677,  .00582,  .00646,  .02361,  .01994,  .01266,  .01544,    FL60790
     *  .01386,  .01968,  .02848,  .03203,  .03854,  .01620,  .01449,    FL60800
     *  .01462,  .01747,  .02513,  .03253,  .03797,  .03019,  .02861,    FL60810
     *  .02538,  .02215,  .01930,  .01797,  .01791,                      FL60820
     *  .01797,  .01539,  .00924,  .00677,  .00384,  .00278, 0.0/        FL60830
      DATA TROA4 /                                                       FL60840
     *  .19829,  .01842,  .01215,  .00791,  .00665,  .00778,  .00652,    FL60850
     *  .00361,  .00253,  .00399,  .02570,  .20690,  .01715,  .00316,    FL60860
     *  .00873,  .00728,  .00658,  .04481,  .03525,  .01646,  .01405,    FL60870
     *  .01310,  .01468,  .01956,  .02184,  .02367,  .01608,  .01816,    FL60880
     *  .02342,  .03203,  .05234,  .06399,  .06538,  .05867,  .05456,    FL60890
     *  .04810,  .03715,  .02949,  .02335,  .02158,                      FL60900
     *  .02149,  .01840,  .01104,  .00809,  .00458,  .00332, 0.0/        FL60910
      DATA TROG1 /                                                       FL60920
     *  .7518,   .6710,   .6638,   .6345,   .6152,   .5736,   .5280,     FL60930
     *  .4949,   .4700,   .4467,   .4204,   .4028,   .3777,   .3563,     FL60940
     *  .3150,   .2919,   .2695,   .2465,   .2402,   .2313,   .2101,     FL60950
     *  .1760,   .1532,   .2091,   .2079,   .1843,   .1811,   .1687,     FL60960
     *  .1626,   .1526,   .1356,   .1030,   .0962,   .1024,   .1086,     FL60970
     *  .0928,   .0836,   .0643,   .0451,   .0290,   .0156,   .0118,     FL60980
     *  .0076,   .0050,   .0024,   .0015,   .0000/                       FL60990
      DATA TROG2 /                                                       FL61000
     *  .7571,   .6858,   .6790,   .6510,   .6315,   .5887,   .5418,     FL61010
     *  .5075,   .4829,   .4598,   .4338,   .4043,   .3890,   .3680,     FL61020
     *  .3259,   .3026,   .2800,   .2541,   .2494,   .2414,   .2196,     FL61030
     *  .1873,   .1657,   .2123,   .2110,   .1890,   .1836,   .1709,     FL61040
     *  .1640,   .1534,   .1354,   .1044,   .0984,   .1026,   .1073,     FL61050
     *  .0935,   .0842,   .0661,   .0477,   .0309,   .0171,   .0129,     FL61060
     *  .0084,   .0056,   .0027,   .0017,   .0000/                       FL61070
      DATA TROG3 /                                                       FL61080
     *  .7667,   .7176,   .7128,   .6879,   .6690,   .6255,   .5769,     FL61090
     *  .5403,   .5167,   .4947,   .4703,   .4143,   .4190,   .3993,     FL61100
     *  .3563,   .3325,   .3095,   .2767,   .2751,   .2693,   .2464,     FL61110
     *  .2175,   .1992,   .2247,   .2215,   .2042,   .1952,   .1814,     FL61120
     *  .1726,   .1604,   .1398,   .1111,   .1065,   .1068,   .1086,     FL61130
     *  .0984,   .0888,   .0724,   .0549,   .0358,   .0216,   .0166,     FL61140
     *  .0109,   .0073,   .0036,   .0023,   .0000/                       FL61150
      DATA TROG4 /                                                       FL61160
     *  .7696,   .7719,   .7710,   .7606,   .7478,   .7142,   .6727,     FL61170
     *  .6381,   .6201,   .6050,   .5912,   .4849,   .5137,   .5019,     FL61180
     *  .4625,   .4389,   .4169,   .3696,   .3707,   .3708,   .3473,     FL61190
     *  .3232,   .3112,   .3022,   .2938,   .2850,   .2675,   .2494,     FL61200
     *  .2347,   .2165,   .1857,   .1536,   .1509,   .1441,   .1416,     FL61210
     *  .1354,   .1245,   .1088,   .0905,   .0614,   .0440,   .0354,     FL61220
     *  .0257,   .0179,   .0089,   .0059,   .0000/                       FL61230
      DATA FG1EXT /                                                      FL61240
     *  .98519,  .99155,  .99089, 1.00000, 1.00580, 1.01740, 1.03170,    FL61250
     * 1.04140, 1.04700, 1.05320, 1.05890, 1.04900, 1.06820, 1.07800,    FL61260
     * 1.09270, 1.10370, 1.11680, 1.10430, 1.11370, 1.12900, 1.14990,    FL61270
     * 1.17210, 1.18280, 1.20140, 1.21260, 1.21950, 1.22680, 1.15590,    FL61280
     * 1.05690,  .98291, 1.01120, 1.10910, 1.11460, 1.14670, 1.16250,    FL61290
     * 1.18540, 1.21580, 1.24610, 1.26840, 1.20500, 1.20850, 1.23340,    FL61300
     * 1.19560, 1.06530,  .68949,  .42888, 0.00000/                      FL61310
      DATA FG1ABS /                                                      FL61320
     *  .00012,  .00001,  .00001,  .00000,  .00001,  .00095,  .01515,    FL61330
     *  .10858,  .03890,  .13270,  .47131,  .49695,  .45787,  .17915,    FL61340
     *  .37373,  .34600,  .31866,  .55187,  .55023,  .49984,  .46341,    FL61350
     *  .45944,  .45916,  .46087,  .46240,  .46386,  .47193,  .48902,    FL61360
     *  .51470,  .53099,  .55264,  .58664,  .58897,  .60369,  .61155,    FL61370
     *  .62335,  .64120,  .65627,  .66278,  .66393,  .69344,  .71087,    FL61380
     *  .67625,  .61180,  .42130,  .29086, 0.00000/                      FL61390
      DATA FG1SYM /                                                      FL61400
     *  .8578,   .8726,   .8722,   .8717,   .8703,   .8652,   .8618,     FL61410
     *  .8798,   .8689,   .8918,   .9641,   .9502,   .9297,   .8544,     FL61420
     *  .9007,   .8885,   .8812,   .9604,   .9470,   .9193,   .9039,     FL61430
     *  .9039,   .9057,   .9110,   .9158,   .9194,   .9381,   .9537,     FL61440
     *  .9595,   .9587,   .9418,   .9101,   .9081,   .8957,   .8898,     FL61450
     *  .8812,   .8685,   .8491,   .8246,   .7815,   .7148,   .6480,     FL61460
     *  .5481,   .4725,   .3457,   .2575,   .0000/                       FL61470
      DATA FG2EXT /                                                      FL61480
     *  .94790,  .96213,  .97061, 1.00000, 1.00940, 1.05180, 1.12520,    FL61490
     * 1.29570, 1.39200, 1.41120, 1.04720, 1.10820, 1.43290, 1.45270,    FL61500
     * 1.18710, 1.04370,  .82356,  .71746,  .92406,  .79342,  .60263,    FL61510
     *  .47680,  .43171,  .36732,  .33259,  .31184,  .24137,  .21603,    FL61520
     *  .24005,  .28816,  .42671,  .56861,  .57263,  .58090,  .57164,    FL61530
     *  .54247,  .43983,  .34475,  .24907,  .19291,  .18500,  .15586,    FL61540
     *  .09047,  .06445,  .03533,  .02529, 0.00000/                      FL61550
      DATA FG2ABS /                                                      FL61560
     *  .00002,  .00000,  .00000,  .00000,  .00000,  .00016,  .00245,    FL61570
     *  .01987,  .00619,  .02323,  .17209,  .57930,  .19812,  .03474,    FL61580
     *  .09636,  .07999,  .06585,  .34591,  .32704,  .17023,  .12635,    FL61590
     *  .11817,  .11624,  .11519,  .11538,  .11600,  .12327,  .14468,    FL61600
     *  .18635,  .24056,  .35412,  .44884,  .45092,  .45215,  .44281,    FL61610
     *  .41778,  .34433,  .27826,  .21066,  .17864,  .17626,  .15028,    FL61620
     *  .08844,  .06358,  .03515,  .02523, 0.00000/                      FL61630
      DATA FG2SYM /                                                      FL61640
     *  .8388,   .8459,   .8419,   .8286,   .8224,   .7883,   .7763,     FL61650
     *  .8133,   .8393,   .8767,   .9258,   .8982,   .7887,   .8082,     FL61660
     *  .8319,   .8243,   .8210,   .8282,   .8037,   .7904,   .7728,     FL61670
     *  .7528,   .7436,   .7274,   .7171,   .7100,   .6790,   .6520,     FL61680
     *  .6305,   .6020,   .5475,   .4577,   .4511,   .4084,   .3872,     FL61690
     *  .3566,   .2976,   .2340,   .1711,   .0956,   .0623,   .0454,     FL61700
     *  .0286,   .0190,   .0090,   .0052,   .0000/                       FL61710
      DATA BSTEXT /                                                      FL61720
     * 1.48671, 1.55462, 1.51506, 1.00000,  .70633,  .28867,  .09994,    FL61730
     *  .04184,  .02728,  .01848,  .01335,  .06513,  .08930,  .06532,    FL61740
     *  .04766,  .04278,  .05810,  .05367,  .04392,  .03342,  .04456,    FL61750
     *  .11867,  .14709,  .12734,  .09291,  .08778,  .05019,  .04070,    FL61760
     *  .05734,  .03576,  .01975,  .01892,  .01956,  .03665,  .04152,    FL61770
     *  .01715,  .01620,  .00835,  .00633,  .00589,                      FL61780
     *  .01393,  .01193,  .00716,  .00526,  .00298,  .00216, 0.0/        FL61790
      DATA BSTABS /                                                      FL61800
     * 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,  .00019,    FL61810
     *  .00127,  .00158,  .00291,  .00405,  .05880,  .08297,  .06019,    FL61820
     *  .04519,  .04133,  .05703,  .05266,  .04304,  .03285,  .04437,    FL61830
     *  .11816,  .14633,  .12639,  .09215,  .08722,  .04968,  .04044,    FL61840
     *  .05709,  .03551,  .01962,  .01892,  .01949,  .03665,  .04146,    FL61850
     *  .01709,  .01620,  .00835,  .00633,  .00589,                      FL61860
     *  .01393,  .01193,  .00716,  .00526,  .00298,  .00216, 0.0/        FL61870
      DATA BSTSYM /                                                      FL61880
     *  .6804,   .7134,   .7253,   .7259,   .6943,   .5918,   .4465,     FL61890
     *  .3223,   .2686,   .2233,   .1916,   .1580,   .1299,   .1108,     FL61900
     *  .0780,   .0629,   .0515,   .0454,   .0426,   .0379,   .0287,     FL61910
     *  .0222,   .0204,   .0206,   .0214,   .0202,   .0205,   .0169,     FL61920
     *  .0150,   .0157,   .0124,   .0083,   .0080,   .0063,   .0062,     FL61930
     *  .0062,   .0043,   .0034,   .0024,   .0013,   .0007,   .0005,     FL61940
     *  .0003,   .0002,   .0001,   .0001,   .0000/                       FL61950
      DATA AVOEXT /                                                      FL61960
     * 1.14880, 1.19171, 1.18013, 1.00000,  .84873,  .53019,  .27968,    FL61970
     *  .14551,  .11070,  .08633,  .07184,  .06076,  .04506,  .03399,    FL61980
     *  .02095,  .01538,  .01266,  .01019,  .00994,  .01044,  .01361,    FL61990
     *  .01791,  .02278,  .02918,  .03108,  .03234,  .03456,  .03184,    FL62000
     *  .02772,  .02475,  .01715,  .01563,  .01665,  .01646,  .01734,    FL62010
     *  .01772,  .01076,  .01051,  .01133,  .01329,                      FL62020
     *  .01492,  .01277,  .00766,  .00562,  .00318,  .00231, 0.0/        FL62030
      DATA AVOABS /                                                      FL62040
     *  .44816,  .11259,  .08500,  .05272,  .04082,  .02449,  .01487,    FL62050
     *  .01019,  .00867,  .00842,  .00842,  .00949,  .00741,  .00487,    FL62060
     *  .00316,  .00335,  .00399,  .00449,  .00525,  .00665,  .01114,    FL62070
     *  .01652,  .02177,  .02437,  .02506,  .02658,  .03006,  .02861,    FL62080
     *  .02513,  .02285,  .01620,  .01532,  .01633,  .01620,  .01709,    FL62090
     *  .01741,  .01057,  .01038,  .01127,  .01329,                      FL62100
     *  .01492,  .01277,  .00766,  .00562,  .00318,  .00231, 0.0/        FL62110
      DATA AVOSYM /                                                      FL62120
     *  .8272,   .7148,   .7076,   .6978,   .6886,   .6559,   .6062,     FL62130
     *  .5561,   .5255,   .4958,   .4729,   .4401,   .4015,   .3699,     FL62140
     *  .3125,   .2773,   .2472,   .2173,   .2054,   .1908,   .1623,     FL62150
     *  .1348,   .1233,   .1615,   .1757,   .1712,   .1521,   .1326,     FL62160
     *  .1230,   .1081,   .0801,   .0528,   .0514,   .0461,   .0446,     FL62170
     *  .0449,   .0415,   .0330,   .0198,   .0097,   .0044,   .0032,     FL62180
     *  .0020,   .0013,   .0006,   .0004,   .0000/                       FL62190
      DATA FVOEXT /                                                      FL62200
     *  .88715,  .92532,  .94013, 1.00000, 1.03013, 1.05975, 1.01171,    FL62210
     *  .88677,  .82538,  .76361,  .71563,  .67424,  .60589,  .55057,    FL62220
     *  .45222,  .37646,  .32316,  .25519,  .22728,  .20525,  .17810,    FL62230
     *  .14481,  .14152,  .37639,  .44551,  .44405,  .42222,  .36462,    FL62240
     *  .32551,  .27519,  .16728,  .10627,  .10861,  .10886,  .11665,    FL62250
     *  .13127,  .10108,  .08557,  .06411,  .05741,                      FL62260
     *  .05531,  .04707,  .02792,  .02028,  .01136,  .00820, 0.0/        FL62270
      DATA FVOABS /                                                      FL62280
     *  .41582,  .22892,  .19108,  .14468,  .12475,  .09158,  .06601,    FL62290
     *  .04943,  .04367,  .04342,  .04399,  .05076,  .04133,  .02829,    FL62300
     *  .01924,  .01981,  .02297,  .02475,  .02778,  .03411,  .05335,    FL62310
     *  .07133,  .08816,  .15342,  .18506,  .19354,  .20791,  .18449,    FL62320
     *  .16101,  .13759,  .08456,  .06886,  .07278,  .07367,  .07956,    FL62330
     *  .08785,  .06032,  .05747,  .05133,  .05323,                      FL62340
     *  .05453,  .04657,  .02773,  .02020,  .01135,  .00820, 0.0/        FL62350
      DATA FVOSYM /                                                      FL62360
     *  .9295,   .8115,   .7897,   .7473,   .7314,   .7132,   .7113,     FL62370
     *  .7238,   .7199,   .7165,   .7134,   .6989,   .6840,   .6687,     FL62380
     *  .6409,   .6325,   .6199,   .6148,   .6142,   .6072,   .5853,     FL62390
     *  .5632,   .5486,   .4753,   .4398,   .4329,   .4091,   .4105,     FL62400
     *  .4120,   .4136,   .4140,   .3637,   .3577,   .3344,   .3220,     FL62410
     *  .3052,   .2957,   .2564,   .2055,   .1229,   .0632,   .0483,     FL62420
     *  .0321,   .0216,   .0103,   .0059,   .0000/                       FL62430
      DATA DMEEXT /                                                      FL62440
     * 1.05019, 1.05880, 1.05259, 1.00000,  .94949,  .81456,  .66051,    FL62450
     *  .54380,  .49133,  .44677,  .41671,  .38063,  .34778,  .32804,    FL62460
     *  .29722,  .27506,  .25082,  .22620,  .21652,  .20253,  .17266,    FL62470
     *  .14905,  .14234,  .14082,  .15057,  .16399,  .23608,  .24481,    FL62480
     *  .27791,  .25076,  .15272,  .09601,  .09456,  .14576,  .12373,    FL62490
     *  .18348,  .12190,  .12924,  .08538,  .04108,                      FL62500
     *  .04714,  .04069,  .02480,  .01789,  .00980,  .00693, 0.0/        FL62510
      DATA DMEABS /                                                      FL62520
     *  .00063,  .00152,  .00184,  .00506,  .00791,  .01829,  .03728,    FL62530
     *  .06158,  .07538,  .08943,  .10051,  .11614,  .13310,  .14348,    FL62540
     *  .14633,  .13728,  .12462,  .11184,  .10709,  .10076,  .09006,    FL62550
     *  .08734,  .09000,  .10304,  .11905,  .13437,  .19551,  .20095,    FL62560
     *  .22494,  .18418,  .09285,  .06665,  .06823,  .12329,  .10551,    FL62570
     *  .16184,  .09835,  .10582,  .06759,  .03247,                      FL62580
     *  .04405,  .03816,  .02327,  .01696,  .00946,  .00677, 0.0/        FL62590
      DATA DMESYM /                                                      FL62600
     *  .7173,   .7039,   .7020,   .6908,   .6872,   .6848,   .6891,     FL62610
     *  .6989,   .7046,   .7099,   .7133,   .7159,   .7134,   .7058,     FL62620
     *  .6827,   .6687,   .6583,   .6513,   .6494,   .6475,   .6467,     FL62630
     *  .6496,   .6506,   .6461,   .6334,   .6177,   .5327,   .5065,     FL62640
     *  .4632,   .4518,   .5121,   .5450,   .5467,   .4712,   .4853,     FL62650
     *  .3984,   .4070,   .3319,   .3427,   .3766,   .3288,   .2969,     FL62660
     *  .2808,   .2661,   .2409,   .2098,   .0000/                       FL62670
      DATA CCUEXT /                                                      FL62680
     *  .98081,  .98746,  .98915, 1.00000, 1.00650, 1.02230, 1.04180,    FL62690
     * 1.05830, 1.06780, 1.07870, 1.09780, 1.06440, 1.09750, 1.11300,    FL62700
     * 1.14320, 1.16660, 1.20540, 1.15420, 1.17610, 1.21910, 1.26990,    FL62710
     * 1.30300, 1.31090, 1.31060, 1.29940, 1.28640, 1.16620,  .98693,    FL62720
     *  .88130,  .83429,  .92012, 1.07340, 1.08150, 1.12680, 1.14770,    FL62730
     * 1.17600, 1.19210, 1.19120, 1.14510,  .97814,  .96308,  .94390,    FL62740
     *  .75994,  .56647,  .26801,  .15748, 0.00000/                      FL62750
      DATA CCUABS /                                                      FL62760
     *  .00007,  .00001,  .00000,  .00000,  .00001,  .00059,  .00956,    FL62770
     *  .07224,  .02502,  .08913,  .41512,  .51824,  .41304,  .12614,    FL62780
     *  .29826,  .26739,  .23672,  .55428,  .55642,  .44494,  .38433,    FL62790
     *  .37277,  .37000,  .36872,  .36896,  .36984,  .37868,  .40498,    FL62800
     *  .44993,  .48941,  .54799,  .60964,  .61302,  .63227,  .64074,    FL62810
     *  .65112,  .65367,  .64760,  .61924,  .59000,  .61601,  .61058,    FL62820
     *  .49236,  .38532,  .20641,  .13474, 0.00000/                      FL62830
      DATA CCUSYM /                                                      FL62840
     *  .8557,   .8676,   .8680,   .8658,   .8630,   .8557,   .8496,     FL62850
     *  .8566,   .8464,   .8627,   .9417,   .9458,   .8891,   .8136,     FL62860
     *  .8503,   .8400,   .8453,   .9428,   .9168,   .8759,   .8733,     FL62870
     *  .8841,   .8894,   .8986,   .9044,   .9082,   .9239,   .9342,     FL62880
     *  .9367,   .9331,   .9119,   .8719,   .8692,   .8515,   .8424,     FL62890
     *  .8287,   .8059,   .7742,   .7354,   .6554,   .5557,   .4720,     FL62900
     *  .3713,   .2990,   .1846,   .1156,   .0000/                       FL62910
      DATA CALEXT /                                                      FL62920
     *  .97331,  .98106,  .98472, 1.00000, 1.00850, 1.03090, 1.05770,    FL62930
     * 1.08070, 1.09390, 1.11530, 1.20260, 1.08250, 1.13480, 1.16770,    FL62940
     * 1.26750, 1.33520, 1.41110, 1.18200, 1.28390, 1.38040, 1.38430,    FL62950
     * 1.31200, 1.26540, 1.17160, 1.10410, 1.05640,  .83383,  .66530,    FL62960
     *  .61995,  .62907,  .77190,  .96660,  .97609, 1.02520, 1.04380,    FL62970
     * 1.06270, 1.02550,  .95714,  .82508,  .63464,  .60962,  .54998,    FL62980
     *  .34165,  .22587,  .10647,  .07067, 0.00000/                      FL62990
      DATA CALABS /                                                      FL63000
     *  .00004,  .00000,  .00000,  .00000,  .00000,  .00036,  .00607,    FL63010
     *  .04771,  .01579,  .05734,  .33199,  .54434,  .35157,  .08528,    FL63020
     *  .21785,  .18813,  .15982,  .52068,  .52125,  .35294,  .28359,    FL63030
     *  .26999,  .26668,  .26477,  .26484,  .26565,  .27546,  .30540,    FL63040
     *  .36011,  .41780,  .51479,  .60420,  .60818,  .62781,  .63339,    FL63050
     *  .63544,  .60762,  .56843,  .50067,  .44739,  .45910,  .42486,    FL63060
     *  .27527,  .19352,  .09932,  .06832, 0.00000/                      FL63070
      DATA CALSYM /                                                      FL63080
     *  .8523,   .8632,   .8623,   .8573,   .8532,   .8422,   .8297,     FL63090
     *  .8252,   .8145,   .8317,   .9312,   .9383,   .8291,   .7640,     FL63100
     *  .8202,   .8276,   .8547,   .9224,   .8859,   .8621,   .8706,     FL63110
     *  .8780,   .8804,   .8833,   .8849,   .8858,   .8889,   .8899,     FL63120
     *  .8872,   .8790,   .8513,   .7984,   .7944,   .7683,   .7545,     FL63130
     *  .7333,   .6939,   .6405,   .5727,   .4313,   .3156,   .2437,     FL63140
     *  .1693,   .1185,   .0574,   .0332,   .0000/                       FL63150
      DATA CSTEXT /                                                      FL63160
     *  .97430,  .98324,  .98570, 1.00000, 1.00890, 1.03100, 1.05590,    FL63170
     * 1.08130, 1.09760, 1.12170, 1.16390, 1.07880, 1.13660, 1.16990,    FL63180
     * 1.22930, 1.26720, 1.31080, 1.15290, 1.23270, 1.29770, 1.31180,    FL63190
     * 1.27830, 1.25190, 1.19190, 1.14390, 1.10790,  .91743,  .74497,    FL63200
     *  .68246,  .67604,  .80234,  .98329,  .99219, 1.03880, 1.05710,    FL63210
     * 1.07730, 1.05460, 1.00640,  .90146,  .71967,  .69823,  .65179,    FL63220
     *  .44906,  .30781,  .14114,  .08913, 0.00000/                      FL63230
      DATA CSTABS /                                                      FL63240
     *  .00005,  .00001,  .00000,  .00000,  .00000,  .00042,  .00681,    FL63250
     *  .05317,  .01779,  .06484,  .35033,  .53843,  .36321,  .09457,    FL63260
     *  .23629,  .20663,  .17789,  .52440,  .52484,  .37331,  .30681,    FL63270
     *  .29375,  .29057,  .28880,  .28887,  .28969,  .29913,  .32789,    FL63280
     *  .37961,  .43212,  .51866,  .60025,  .60398,  .62285,  .62874,    FL63290
     *  .63229,  .61185,  .58151,  .52536,  .47993,  .49571,  .47074,    FL63300
     *  .33104,  .24066,  .12346,  .08312, 0.00000/                      FL63310
      DATA CSTSYM /                                                      FL63320
     *  .8519,   .8633,   .8629,   .8590,   .8546,   .8432,   .8328,     FL63330
     *  .8330,   .8251,   .8439,   .9332,   .9388,   .8422,   .7823,     FL63340
     *  .8288,   .8291,   .8482,   .9255,   .8906,   .8613,   .8675,     FL63350
     *  .8772,   .8810,   .8869,   .8905,   .8927,   .9016,   .9069,     FL63360
     *  .9060,   .8989,   .8714,   .8204,   .8168,   .7932,   .7811,     FL63370
     *  .7628,   .7319,   .6905,   .6401,   .5324,   .4233,   .3459,     FL63380
     *  .2636,   .2027,   .1120,   .0663,   .0000/                       FL63390
      DATA CSCEXT /                                                      FL63400
     *  .96965,  .97960,  .98266, 1.00000, 1.01040, 1.03530, 1.06590,    FL63410
     * 1.09980, 1.12280, 1.16020, 1.20330, 1.08630, 1.16840, 1.21860,    FL63420
     * 1.28860, 1.32310, 1.33780, 1.11630, 1.24450, 1.30260, 1.26260,    FL63430
     * 1.17670, 1.12990, 1.04180,  .98070,  .93828,  .74401,  .59962,    FL63440
     *  .56489,  .57976,  .72193,  .90905,  .91772,  .96075,  .97500,    FL63450
     *  .98623,  .93761,  .86388,  .73722,  .56926,  .54699,  .49341,    FL63460
     *  .31131,  .20846,  .09872,  .06531, 0.00000/                      FL63470
      DATA CSCABS /                                                      FL63480
     *  .00004,  .00000,  .00000,  .00000,  .00000,  .00035,  .00553,    FL63490
     *  .04382,  .01430,  .05271,  .30881,  .54982,  .32983,  .07796,    FL63500
     *  .20033,  .17269,  .14662,  .49557,  .49304,  .32632,  .26104,    FL63510
     *  .24829,  .24525,  .24349,  .24358,  .24437,  .25378,  .28239,    FL63520
     *  .33510,  .39227,  .49203,  .58265,  .58638,  .60338,  .60677,    FL63530
     *  .60472,  .56954,  .52556,  .45708,  .40717,  .41646,  .38375,    FL63540
     *  .25009,  .17726,  .09148,  .06291, 0.00000/                      FL63550
      DATA CSCSYM /                                                      FL63560
     *  .8495,   .8597,   .8594,   .8535,   .8479,   .8349,   .8214,     FL63570
     *  .8192,   .8151,   .8395,   .9321,   .9329,   .8156,   .7722,     FL63580
     *  .8270,   .8319,   .8533,   .9138,   .8772,   .8562,   .8628,     FL63590
     *  .8691,   .8713,   .8742,   .8759,   .8768,   .8805,   .8818,     FL63600
     *  .8783,   .8685,   .8362,   .7776,   .7734,   .7458,   .7317,     FL63610
     *  .7106,   .6738,   .6250,   .5655,   .4409,   .3338,   .2655,     FL63620
     *  .1947,   .1427,   .0727,   .0422,   .0000/                       FL63630
      DATA CNIEXT /                                                      FL63640
     *  .97967,  .98623,  .98795, 1.00000, 1.00710, 1.02340, 1.04300,    FL63650
     * 1.06100, 1.07130, 1.08440, 1.10650, 1.06540, 1.10200, 1.12040,    FL63660
     * 1.15490, 1.17990, 1.21730, 1.15000, 1.18140, 1.22610, 1.26770,    FL63670
     * 1.28840, 1.29070, 1.28200, 1.26650, 1.25130, 1.12860,  .95670,    FL63680
     *  .85784,  .81564,  .90486, 1.05950, 1.06760, 1.11240, 1.13250,    FL63690
     * 1.15910, 1.16960, 1.16290, 1.11130,  .94771,  .93251,  .91151,    FL63700
     *  .73279,  .55018,  .26554,  .15656, 0.00000/                      FL63710
      DATA CNIABS /                                                      FL63720
     *  .00007,  .00001,  .00000,  .00000,  .00001,  .00058,  .00948,    FL63730
     *  .07084,  .02436,  .08711,  .40714,  .52024,  .40688,  .12335,    FL63740
     *  .29163,  .26107,  .23098,  .54886,  .55047,  .43579,  .37552,    FL63750
     *  .36411,  .36140,  .36017,  .36043,  .36132,  .37019,  .39640,    FL63760
     *  .44146,  .48184,  .54304,  .60651,  .60988,  .62882,  .63682,    FL63770
     *  .64613,  .64572,  .63682,  .60584,  .57559,  .60014,  .59283,    FL63780
     *  .47587,  .37364,  .20267,  .13269, 0.00000/                      FL63790
      DATA CNISYM /                                                      FL63800
     *  .8550,   .8670,   .8677,   .8645,   .8616,   .8538,   .8474,     FL63810
     *  .8534,   .8439,   .8609,   .9411,   .9449,   .8822,   .8101,     FL63820
     *  .8486,   .8403,   .8475,   .9405,   .9134,   .8749,   .8732,     FL63830
     *  .8833,   .8882,   .8968,   .9025,   .9061,   .9217,   .9322,     FL63840
     *  .9346,   .9308,   .9086,   .8669,   .8641,   .8457,   .8364,     FL63850
     *  .8222,   .7992,   .7677,   .7298,   .6525,   .5558,   .4752,     FL63860
     *  .3796,   .3105,   .1995,   .1287,   .0000/                       FL63870
C                                                                        FL63880
C     EXTINCTION  COEFFICIENTS                                           FL63890
C                                                                        FL63900
      DATA CI64XT/                                                       FL63910
     *   9.947E-01,  9.968E-01,  9.972E-01,  1.000E+00,  1.002E+00,      FL63920
     *   1.005E+00,  1.010E+00,  1.013E+00,  1.016E+00,  1.018E+00,      FL63930
     *   1.019E+00,  1.016E+00,  1.023E+00,  1.026E+00,  1.030E+00,      FL63940
     *   1.033E+00,  1.036E+00,  1.037E+00,  1.038E+00,  1.040E+00,      FL63950
     *   1.043E+00,  1.047E+00,  1.049E+00,  1.051E+00,  1.052E+00,      FL63960
     *   1.053E+00,  1.055E+00,  1.032E+00,  1.034E+00,  1.047E+00,      FL63970
     *   1.060E+00,  1.074E+00,  1.075E+00,  1.081E+00,  1.085E+00,      FL63980
     *   1.090E+00,  1.102E+00,  1.117E+00,  1.131E+00,  1.094E+00,      FL63990
     *   1.168E+00,  1.187E+00,  1.244E+00,  1.297E+00,  1.475E+00,      FL64000
     *   1.695E+00,  1.556E+00 /                                         FL64010
C                                                                        FL64020
C     ABSORPTION  COEFFICIENTS                                           FL64030
C                                                                        FL64040
      DATA CI64AB/                                                       FL64050
     *   7.893E-05,  1.914E-05,  1.450E-05,  5.904E-06,  3.905E-05,      FL64060
     *   1.917E-03,  2.604E-01,  3.732E-01,  8.623E-02,  2.253E-01,      FL64070
     *   4.152E-01,  4.460E-01,  4.660E-01,  4.589E-01,  4.848E-01,      FL64080
     *   4.786E-01,  4.915E-01,  4.944E-01,  4.936E-01,  4.947E-01,      FL64090
     *   4.978E-01,  5.012E-01,  5.028E-01,  5.070E-01,  5.095E-01,      FL64100
     *   5.111E-01,  5.205E-01,  5.126E-01,  4.969E-01,  4.868E-01,      FL64110
     *   4.836E-01,  4.982E-01,  4.999E-01,  5.097E-01,  5.126E-01,      FL64120
     *   5.188E-01,  5.108E-01,  4.915E-01,  5.559E-01,  5.515E-01,      FL64130
     *   5.600E-01,  5.948E-01,  6.225E-01,  6.348E-01,  5.693E-01,      FL64140
     *   3.306E-01,  8.661E-02 /                                         FL64150
C                                                                        FL64160
C     ASYMMETRY  PARAMETER  -  G                                         FL64170
C                                                                        FL64180
      DATA CI64G/                                                        FL64190
     *   .8626,  .8824,  .8851,  .8893,  .8904,  .8913,  .9332,  .9549,  FL64200
     *   .9141,  .9407,  .9763,  .9428,  .9509,  .9580,  .9699,  .9679,  FL64210
     *   .9735,  .9737,  .9717,  .9712,  .9712,  .9715,  .9721,  .9744,  FL64220
     *   .9756,  .9764,  .9822,  .9849,  .9721,  .9530,  .9341,  .9352,  FL64230
     *   .9366,  .9426,  .9425,  .9448,  .9365,  .9256,  .9485,  .9417,  FL64240
     *   .8868,  .8983,  .8589,  .8115,  .6810,  .5923,  .5703 /         FL64250
C                                                                        FL64260
C     EXTINCTION COEFFICIENTS                                            FL64270
C                                                                        FL64280
      DATA CIR4XT/                                                       FL64290
     *   9.685E-01,  9.803E-01,  9.826E-01,  1.000E+00,  1.011E+00,      FL64300
     *   1.038E+00,  1.066E+00,  1.090E+00,  1.118E+00,  1.201E+00,      FL64310
     *   1.374E+00,  1.019E+00,  1.143E+00,  1.198E+00,  1.331E+00,      FL64320
     *   1.434E+00,  1.424E+00,  1.283E+00,  1.298E+00,  1.326E+00,      FL64330
     *   1.287E+00,  1.230E+00,  1.191E+00,  1.048E+00,  9.634E-01,      FL64340
     *   9.093E-01,  6.067E-01,  5.216E-01,  6.953E-01,  8.902E-01,      FL64350
     *   1.083E+00,  1.228E+00,  1.214E+00,  1.076E+00,  1.032E+00,      FL64360
     *   8.881E-01,  6.275E-01,  3.462E-01,  2.118E-01,  3.955E-01,      FL64370
     *   5.089E-01,  3.012E-01,  1.235E-01,  5.377E-02,  2.068E-02,      FL64380
     *   6.996E-03,  1.560E-03 /                                         FL64390
C                                                                        FL64400
C     ABSORPTION  COEFFICIENTS                                           FL64410
C                                                                        FL64420
      DATA CIR4AB/                                                       FL64430
     *   5.316E-06,  1.461E-06,  9.045E-07,  4.431E-07,  2.746E-06,      FL64440
     *   1.413E-04,  2.920E-02,  5.578E-02,  6.844E-03,  2.151E-02,      FL64450
     *   6.322E-02,  5.051E-01,  4.578E-01,  1.360E-01,  3.269E-01,      FL64460
     *   1.572E-01,  2.246E-01,  4.176E-01,  4.282E-01,  3.802E-01,      FL64470
     *   3.517E-01,  3.037E-01,  2.543E-01,  2.410E-01,  2.432E-01,      FL64480
     *   2.438E-01,  2.346E-01,  3.747E-01,  4.839E-01,  5.722E-01,      FL64490
     *   6.368E-01,  5.303E-01,  5.085E-01,  3.920E-01,  3.437E-01,      FL64500
     *   2.481E-01,  1.175E-01,  7.172E-02,  1.108E-01,  3.459E-01,      FL64510
     *   4.044E-01,  2.545E-01,  9.594E-02,  4.410E-02,  1.887E-02,      FL64520
     *   6.433E-03,  1.456E-03 /                                         FL64530
C                                                                        FL64540
C     ASYMMETRY  PARAMETER  -  G                                         FL64550
C                                                                        FL64560
      DATA CIR4G/                                                        FL64570
     *   .8517,  .8654,  .8661,  .8615,  .8574,  .8447,  .8321,  .8248,  FL64580
     *   .8227,  .8612,  .9363,  .9231,  .8419,  .7550,  .8481,  .8358,  FL64590
     *   .8718,  .8953,  .8884,  .8786,  .8731,  .8660,  .8625,  .8652,  FL64600
     *   .8659,  .8658,  .8676,  .8630,  .8434,  .8194,  .7882,  .7366,  FL64610
     *   .7339,  .7161,  .7015,  .6821,  .6383,  .5823,  .4845,  .2977,  FL64620
     *   .2295,  .1716,  .1228,  .0748,  .0329,  .0186,  .0081 /         FL64630
      END                                                                FL64640
C
C     ******************************************************************
C
      SUBROUTINE RDEXA                                                   FL64650
C                                                                        FL64660
C     READ IN USER DEFINED EXTINCTION, ABSORPTION AND                    FL64670
C     ASYMMETRY PARAMETERS                                               FL64680
C                                                                        FL64690
C
C     MXFSC IS THE MAXIMUM NUMBER OF LAYERS FOR OUTPUT TO LBLRTM
C     MXLAY IS THE MAXIMUN NUMBER OF OUTPUT LAYERS
C     MXZMD IS THE MAX NUMBER OF LEVELS IN THE ATMOSPHERIC PROFILE
C         STORED IN ZMDL (INPUT)
C     MXPDIM IS THE MAXIMUM NUMBER OF LEVELS IN THE PROFILE ZPTH
C         OBTAINED BY MERGING ZMDL AND ZOUT
C     MXMOL IS THE MAXIMUM NUMBER OF MOLECULES, KMXNOM IS THE DEFAULT
C
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
C
C     BLANK COMMON FOR ZMDL
C
      COMMON RELHUM(MXZMD),HSTOR(MXZMD),ICH(4),VH(16),TX(16),W(16)       FL64700
      COMMON WPATH(IM2,16),TBBY(IM2)                                     FL64710
      COMMON ABSC(5,47),EXTC(5,47),ASYM(5,47),VX2(47),AWCCON(5)
C
      CHARACTER*8      HMOD                                              FL64730
C
      COMMON /CMN/ HMOD(3),ZM(MXZMD),PF(MXZMD),TF(MXZMD),RFNDXM(MXZMD),
     *          ZP(IM2),PP(IM2),TP(IM2),RFNDXP(IM2),SP(IM2),PPSUM(IM2),
     *          TPSUM(IM2),RHOPSM(IM2),IMLOW,WGM(MXZMD),DENW(MXZMD),
     *          AMTP(MXMOL,MXPDIM)                                        
C
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL64780
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL64790
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL64800
      COMMON /LCRD2D/ IREG(4),ALTB(4),IREGC(4)                           FL64810
C
      DIMENSION TITLE(18),VX(47)                                         FL64820
C                                                                        FL64830
      READ (IRD,900) (IREG(IK),IK=1,4)                                   FL64840
      WRITE (IPR,905) (IREG(IK),IK=1,4)                                  FL64850
C                                                                        FL64860
      DO 10 IHC = 1, 4                                                   FL64870
C                                                                        FL64880
         IF (IREG(IHC).EQ.0) GO TO 10                                    FL64890
         READ (IRD,910) AWCCON(IHC),TITLE                                FL64900
         WRITE (IPR,915) AWCCON(IHC),TITLE                               FL64910
         WRITE (IPR,920)                                                 FL64920
C                                                                        FL64930
         READ (IRD,925) (VX(I),EXTC(IHC,I),ABSC(IHC,I),ASYM(IHC,I),I=1,  FL64940
     *      47)                                                          FL64950
         WRITE (IPR,930) (VX(I),EXTC(IHC,I),ABSC(IHC,I),ASYM(IHC,I),I=1, FL64960
     *      47)                                                          FL64970
   10 CONTINUE                                                           FL64980
      RETURN                                                             FL64990
C                                                                        FL65000
  900 FORMAT(4I5)                                                        FL65010
  905 FORMAT('0 RECORD 3.6.2 *****',4I5)                                 FL65020
  910 FORMAT(E10.3,18A4)                                                 FL65030
  915 FORMAT('0 RECORD 3.6.2 **** EQUIVALENT WATER = ',1PE10.3,18A4)     FL65040
  920 FORMAT('0 RECORD 3.6.3 ****')                                      FL65050
  925 FORMAT(3(F6.2,2F7.5,F6.4))                                         FL65060
  930 FORMAT(2X,F6.2,2F7.5,F6.4,F6.2,2F7.5,F6.4,F6.2,2F7.5,F6.4)         FL65070
C                                                                        FL65080
      END                                                                FL65090
C
C     *****************************************************************
C
      SUBROUTINE MARINE(VIS,MODEL,WS,WH,ICSTL,BEXT,BABS,NL)              FL65100
C                                                                        FL65110
C     THIS SUBROUTINE DETERMINES AEROSOL EXT + ABS COEFFICIENTS          FL65120
C     FOR THE NAVY MARITIME MODEL                                        FL65130
C     CODED BY STU GATHMAN                  -  NRL                       FL65140
C                                                                        FL65150
C     INPUTS-                                                            FL65160
C     WSS = CURRENT WIND SPEED (M/S)                                     FL65170
C     WHH = 24 HOUR AVERAGE WIND SPEED (M/S)                             FL65180
C     RHH = RELATIVE HUMIDITY (PERCENTAGE)                               FL65190
C     VIS = METEOROLOGICAL RANGE (KM)                                    FL65200
C     ICTL = AIR MASS CHARACTER  1 = OPEN OCEAN                          FL65210
C     10 = STRONG CONTINENTAL INFLUENCE                                  FL65220
C     MODEL = MODEL ATMOSPHERE                                           FL65230
C                                                                        FL65240
C     OUTPUTS-                                                           FL65250
C     BEXT = EXTINCTION COEFFICIENT (KM-1)                               FL65260
C     BABS = ABSORPTION COEFFICIENT (KM-1)                               FL65270
C                                                                        FL65280
      COMMON /MART/ RHH                                                  FL65290
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL65300
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL65310
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL65320
      COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP                      FL65330
      COMMON/A/T1QEXT(40,4),T2QEXT(40,4),T3QEXT(40,4),                   FL65340
     *     T1QABS(40,4),T2QABS(40,4),T3QABS(40,4),ALAM(40),AREL(4)       FL65350
C                                                                        FL65360
C     C    COMMON/AER/A1, A2, A3        X(5)                             FL65370
C                                                                        FL65380
      DIMENSION WSPD(8), BEXT(5,47), BABS(5,47)                          FL65390
      DIMENSION RHD(8)                                                   FL65400
C
      DATA WSPD/6.9, 4.1, 4.1, 10.29, 6.69, 12.35, 7.2, 6.9/             FL65410
      DATA RHD/80., 75.63, 76.2, 77.13, 75.24, 80.53, 45.89, 80./        FL65420
C
      PISC = PI/1000.0                                                   FL65430
      WRITE (IPR,900)                                                    FL65440
C                                                                        FL65450
C     CHECK LIMITS OF MODEL VALIDITY                                     FL65460
C                                                                        FL65470
      RH = RHH                                                           FL65480
      IF (RHH.GT.0.) GO TO 10                                            FL65490
      RH = RHD(MODEL+1)                                                  FL65500
   10 WS = MIN(WS,20.)                                                   FL65510
      WH = MIN(WH,20.)                                                   FL65520
      RH = MIN(RH,98.)                                                   FL65530
      IF (RH.LT.50.0.AND.RH.GE.0.0) RH = 50.                             FL65540
      IF (ICSTL.LT.1.OR.ICSTL.GT.10) ICSTL = 3                           FL65550
C                                                                        FL65560
C     FIND SIZE DISTRIBUTION PARAMETERS FROM METEOROLOGY INPUT           FL65570
C                                                                        FL65580
      IF (WH.LE.0.) WRITE (IPR,915)                                      FL65590
      IF (WH.LE.0.0) WH = WSPD(MODEL+1)                                  FL65600
      IF (WS.LE.0.) WRITE (IPR,920)                                      FL65610
      IF (WS.LE.0.0) WS = WH                                             FL65620
      WRITE (IPR,910) WS,WH,RH,ICSTL                                     FL65630
C                                                                        FL65640
C     F IS A RELATIVE HUMIDITY DEPENDENT GROWTH CORRECTION               FL65650
C     TO THE ATTENUATION COEFFICIENT.                                    FL65660
C                                                                        FL65670
      F = ((2.-RH/100.)/(6.*(1.-RH/100.)))**0.33333                      FL65680
      A1 = 2000.0*ICSTL*ICSTL                                            FL65690
      A2 =   MAX(5.866*(WH-2.2),0.5)                                     FL65700
C                                                                        FL65710
C     CC   A3 =   MAX(0.01527*(WS-2.2), 1.14E-5)                         FL65720
C                                                                        FL65730
      A3 = 10**(0.06*WS-2.8)                                             FL65740
C                                                                        FL65750
C     FIND EXTINCTION AT 0.55 MICRONS AND NORMALIZE TO 1.                FL65760
C                                                                        FL65770
C     INTERPOLATE FOR RELATIVE HUMIDITY                                  FL65780
C                                                                        FL65790
      DO 20 J = 2, 4                                                     FL65800
         IF (RH.LE.AREL(J)) GO TO 30                                     FL65810
   20 CONTINUE                                                           FL65820
   30 DELRH = AREL(J)-AREL(J-1)                                          FL65830
      DELRHV = RH-AREL(J-1)                                              FL65840
      RATIO = DELRHV/DELRH                                               FL65850
      QE1 = T1QEXT(4,J-1)+(T1QEXT(4,J)-T1QEXT(4,J-1))*RATIO              FL65860
      QE2 = T2QEXT(4,J-1)+(T2QEXT(4,J)-T2QEXT(4,J-1))*RATIO              FL65870
      QE3 = T3QEXT(4,J-1)+(T3QEXT(4,J)-T3QEXT(4,J-1))*RATIO              FL65880
      TOTAL = A1*10.**QE1+A2*10.**QE2+A3*10.**QE3                        FL65890
      EXT55 = PISC*TOTAL/F                                               FL65900
C                                                                        FL65910
C     IF METEOROLOLICAL RANGE NOT SPECIFIED,FIND FROM METEOR DATA        FL65920
C                                                                        FL65930
      IF (VIS.LE.0.) VIS = 3.912/(EXT55+0.01159)                         FL65940
      C = (1./EXT55)*(PISC/F)                                            FL65950
      A1 = C*A1                                                          FL65960
      A2 = C*A2                                                          FL65970
      A3 = C*A3                                                          FL65980
C                                                                        FL65990
C     CALCULATE NORMALIZED ATTENUATION COEFICIENTS                       FL66000
C                                                                        FL66010
      DO 40 I = 1, 40                                                    FL66020
         T1XV = T1QEXT(I,J-1)+(T1QEXT(I,J)-T1QEXT(I,J-1))*RATIO          FL66030
         T2XV = T2QEXT(I,J-1)+(T2QEXT(I,J)-T2QEXT(I,J-1))*RATIO          FL66040
         T3XV = T3QEXT(I,J-1)+(T3QEXT(I,J)-T3QEXT(I,J-1))*RATIO          FL66050
         T1AV = T1QABS(I,J-1)+(T1QABS(I,J)-T1QABS(I,J-1))*RATIO          FL66060
         T2AV = T2QABS(I,J-1)+(T2QABS(I,J)-T2QABS(I,J-1))*RATIO          FL66070
         T3AV = T3QABS(I,J-1)+(T3QABS(I,J)-T3QABS(I,J-1))*RATIO          FL66080
         BEXT(NL,I) = A1*10**(T1XV)+A2*10**(T2XV)+A3*10**(T3XV)          FL66090
         BABS(NL,I) = A1*10**(T1AV)+A2*10**(T2AV)+A3*10**(T3AV)          FL66100
   40 CONTINUE                                                           FL66110
      WRITE (IPR,905) VIS                                                FL66120
      RETURN                                                             FL66130
C                                                                        FL66140
  900 FORMAT('0MARINE AEROSOL MODEL USED')                               FL66150
  905 FORMAT('0',T10,'VIS = ',F10.2,' KM')                               FL66160
  910 FORMAT(T10,'WIND SPEED = ',F8.2,' M/SEC',/,T10,                    FL66170
     * 'WIND SPEED (24 HR AVERAGE) = ',F8.2,' M/SEC',/,                  FL66180
     * T10,'RELATIVE HUMIDITY = ',F8.2,' PERCENT',/,                     FL66190
     * T10,'AIRMASS CHARACTER =' ,I3)                                    FL66200
  915 FORMAT('0  WS NOT SPECIFIED, A DEFAULT VALUE IS USED')             FL66210
  920 FORMAT('0  WH NOT SPECIFIED, A DEFAULT VALUE IS USED')             FL66220
C                                                                        FL66230
      END                                                                FL66240
      BLOCK DATA MARDTA                                                  FL66250
C                                                                        FL66260
C     >    BLOCK DATA                                                    FL66270
C                                                                        FL66280
C     MARINE AEROSOL EXTINCTION AND ABSORPTION DATA                      FL66290
C     CODED BY STU GATHMAN                  -  NRL                       FL66300
C                                                                        FL66310
      COMMON/A/T1QEXT(40,4),T2QEXT(40,4),T3QEXT(40,4),                   FL66320
     *T1QABS(40,4),T2QABS(40,4),T3QABS(40,4),ALAM(40),AREL(4)            FL66330
      DIMENSION A1(40),A2(40),A3(40),A4(40)                              FL66340
      DIMENSION B1(40),B2(40),B3(40),B4(40)                              FL66350
      DIMENSION C1(40),C2(40),C3(40),C4(40)                              FL66360
      DIMENSION D1(40),D2(40),D3(40),D4(40)                              FL66370
      DIMENSION E1(40),E2(40),E3(40),E4(40)                              FL66380
      DIMENSION F1(40),F2(40),F3(40),F4(40)                              FL66390
      EQUIVALENCE (A1(1), T1QEXT(1,1)), (A2(1), T1QEXT(1,2)),            FL66400
     *            (A3(1), T1QEXT(1,3)), (A4(1), T1QEXT(1,4))             FL66410
      EQUIVALENCE (B1(1), T2QEXT(1,1)), (B2(1), T2QEXT(1,2)),            FL66420
     *            (B3(1), T2QEXT(1,3)), (B4(1), T2QEXT(1,4))             FL66430
      EQUIVALENCE (C1(1), T3QEXT(1,1)), (C2(1), T3QEXT(1,2)),            FL66440
     *            (C3(1), T3QEXT(1,3)), (C4(1), T3QEXT(1,4))             FL66450
      EQUIVALENCE (D1(1), T1QABS(1,1)), (D2(1), T1QABS(1,2)),            FL66460
     *            (D3(1), T1QABS(1,3)), (D4(1), T1QABS(1,4))             FL66470
      EQUIVALENCE (E1(1), T2QABS(1,1)), (E2(1), T2QABS(1,2)),            FL66480
     *            (E3(1), T2QABS(1,3)), (E4(1), T2QABS(1,4))             FL66490
      EQUIVALENCE (F1(1), T3QABS(1,1)), (F2(1), T3QABS(1,2)),            FL66500
     *            (F3(1), T3QABS(1,3)), (F4(1), T3QABS(1,4))             FL66510
      DATA AREL/50.,85.,95.,98./                                         FL66520
      DATA ALAM/                                                         FL66530
     * 0.2000,   0.3000,   0.3371,   0.5500,   0.6943,   1.0600,         FL66540
     * 1.5360,   2.0000,   2.2500,   2.5000,   2.7000,   3.0000,         FL66550
     * 3.3923,   3.7500,   4.5000,   5.0000,   5.5000,   6.0000,         FL66560
     * 6.2000,   6.5000,   7.2000,   7.9000,   8.2000,   8.7000,         FL66570
     * 9.0000,   9.2000,  10.0000,  10.5910,  11.0000,  11.5000,         FL66580
     *12.5000,  14.8000,  15.0000,  16.4000,  17.2000,  18.5000,         FL66590
     *21.3000,  25.0000,  30.0000,  40.0000/                             FL66600
      DATA A1/                                                           FL66610
     *-3.2949,  -3.4662,  -3.5275,  -3.8505,  -4.0388,  -4.4410,         FL66620
     *-4.8584,  -5.1720,  -5.3272,  -5.4342,  -5.2765,  -4.5101,         FL66630
     *-5.3730,  -5.7468,  -5.7579,  -5.8333,  -5.8552,  -5.1780,         FL66640
     *-5.2910,  -5.5959,  -5.6295,  -5.6748,  -5.6051,  -5.5363,         FL66650
     *-5.5330,  -5.5136,  -5.6568,  -5.6040,  -5.5221,  -5.3902,         FL66660
     *-5.1724,  -5.0903,  -5.0901,  -5.1285,  -5.1444,  -5.1963,         FL66670
     *-5.3101,  -5.3994,  -5.4873,  -5.4779/                             FL66680
      DATA A2/                                                           FL66690
     *-2.8302,  -2.9446,  -2.9904,  -3.2510,  -3.4104,  -3.7635,         FL66700
     *-4.1452,  -4.4466,  -4.6160,  -4.7772,  -4.7030,  -3.8461,         FL66710
     *-4.6466,  -5.0105,  -5.0747,  -5.1810,  -5.2705,  -4.5537,         FL66720
     *-4.6594,  -4.9872,  -5.0872,  -5.1229,  -5.0985,  -5.0623,         FL66730
     *-5.0544,  -5.0407,  -5.0793,  -4.9796,  -4.8748,  -4.7298,         FL66740
     *-4.5063,  -4.4260,  -4.4280,  -4.4650,  -4.4912,  -4.5474,         FL66750
     *-4.6672,  -4.7711,  -4.8814,  -4.9073/                             FL66760
      DATA A3/                                                           FL66770
     *-2.3712,  -2.4231,  -2.4512,  -2.6377,  -2.7631,  -3.0569,         FL66780
     *-3.3918,  -3.6682,  -3.8305,  -4.0111,  -4.0467,  -3.2055,         FL66790
     *-3.8717,  -4.1908,  -4.3282,  -4.4495,  -4.5780,  -3.9249,         FL66800
     *-4.0136,  -4.3349,  -4.4674,  -4.5088,  -4.5083,  -4.4973,         FL66810
     *-4.4923,  -4.4845,  -4.4753,  -4.3617,  -4.2509,  -4.1029,         FL66820
     *-3.8779,  -3.7963,  -3.7989,  -3.8345,  -3.8639,  -3.9215,         FL66830
     *-4.0438,  -4.1532,  -4.2719,  -4.3120/                             FL66840
      DATA A4/                                                           FL66850
     *-1.9911,  -1.9989,  -2.0126,  -2.1342,  -2.2283,  -2.4663,         FL66860
     *-2.7552,  -3.0036,  -3.1528,  -3.3328,  -3.4468,  -2.6649,         FL66870
     *-3.1986,  -3.4769,  -3.6571,  -3.7821,  -3.9284,  -3.3776,         FL66880
     *-3.4435,  -3.7436,  -3.8910,  -3.9455,  -3.9573,  -3.9633,         FL66890
     *-3.9639,  -3.9610,  -3.9427,  -3.8304,  -3.7203,  -3.5733,         FL66900
     *-3.3489,  -3.2650,  -3.2675,  -3.3017,  -3.3317,  -3.3893,         FL66910
     *-3.5126,  -3.6243,  -3.7467,  -3.7927/                             FL66920
      DATA B1/                                                           FL66930
     *-0.5781,  -0.5525,  -0.5484,  -0.5147,  -0.5094,  -0.5324,         FL66940
     *-0.6138,  -0.7139,  -0.7776,  -0.8624,  -0.9838,  -0.7720,         FL66950
     *-0.8542,  -0.9535,  -1.0873,  -1.1624,  -1.2647,  -1.2123,         FL66960
     *-1.1811,  -1.2905,  -1.4126,  -1.4643,  -1.5227,  -1.4560,         FL66970
     *-1.4177,  -1.4144,  -1.5746,  -1.6348,  -1.6431,  -1.6023,         FL66980
     *-1.4648,  -1.3910,  -1.3898,  -1.4056,  -1.4196,  -1.4655,         FL66990
     *-1.5795,  -1.6825,  -1.7924,  -1.8224/                             FL67000
      DATA B2/                                                           FL67010
     *-0.1809,  -0.1651,  -0.1566,  -0.1258,  -0.1113,  -0.1046,         FL67020
     *-0.1468,  -0.2157,  -0.2679,  -0.3480,  -0.4988,  -0.2657,         FL67030
     *-0.2991,  -0.3924,  -0.5266,  -0.5983,  -0.7037,  -0.6671,         FL67040
     *-0.6074,  -0.7134,  -0.8352,  -0.9080,  -0.9577,  -0.9579,         FL67050
     *-0.9542,  -0.9629,  -1.0867,  -1.1219,  -1.1032,  -1.0330,         FL67060
     *-0.8663,  -0.7677,  -0.7667,  -0.7768,  -0.7919,  -0.8304,         FL67070
     *-0.9354,  -1.0400,  -1.1640,  -1.2357/                             FL67080
      DATA B3/                                                           FL67090
     * 0.2483,   0.2574,   0.2626,   0.2887,   0.3055,   0.3312,         FL67100
     * 0.3262,   0.2922,   0.2589,   0.1989,   0.0548,   0.2322,         FL67110
     * 0.2487,   0.1816,   0.0685,   0.0090,  -0.0846,  -0.0876,         FL67120
     *-0.0110,  -0.0936,  -0.2013,  -0.2799,  -0.3216,  -0.3575,         FL67130
     *-0.3769,  -0.3944,  -0.5018,  -0.5379,  -0.5179,  -0.4473,         FL67140
     *-0.2822,  -0.1730,  -0.1713,  -0.1737,  -0.1850,  -0.2141,         FL67150
     *-0.3046,  -0.4002,  -0.5221,  -0.6163/                             FL67160
      DATA B4/                                                           FL67170
     * 0.6276,   0.6324,   0.6363,   0.6570,   0.6715,   0.7006,         FL67180
     * 0.7172,   0.7091,   0.6925,   0.6543,   0.5356,   0.6473,         FL67190
     * 0.6924,   0.6516,   0.5661,   0.5206,   0.4440,   0.4091,         FL67200
     * 0.4902,   0.4325,   0.3427,   0.2691,   0.2336,   0.1872,         FL67210
     * 0.1593,   0.1386,   0.0348,  -0.0131,  -0.0031,   0.0566,         FL67220
     * 0.2093,   0.3214,   0.3238,   0.3278,   0.3211,   0.3007,         FL67230
     * 0.2257,   0.1426,   0.0304,  -0.0739/                             FL67240
      DATA C1/                                                           FL67250
     * 2.1434,   2.1454,   2.1469,   2.1539,   2.1577,   2.1673,         FL67260
     * 2.1812,   2.1970,   2.2030,   2.2115,   2.2149,   2.1931,         FL67270
     * 2.2220,   2.2326,   2.2425,   2.2479,   2.2494,   2.2203,         FL67280
     * 2.2382,   2.2473,   2.2380,   2.2373,   2.2179,   2.2310,         FL67290
     * 2.2417,   2.2421,   2.2244,   2.1950,   2.1686,   2.1370,         FL67300
     * 2.1193,   2.1454,   2.1477,   2.1703,   2.1725,   2.1729,         FL67310
     * 2.1580,   2.1324,   2.0878,   2.0131/                             FL67320
      DATA C2/                                                           FL67330
     * 2.5480,   2.5512,   2.5511,   2.5562,   2.5601,   2.5669,         FL67340
     * 2.5792,   2.5874,   2.5950,   2.6022,   2.6081,   2.5875,         FL67350
     * 2.6093,   2.6184,   2.6319,   2.6391,   2.6439,   2.6138,         FL67360
     * 2.6319,   2.6437,   2.6442,   2.6421,   2.6336,   2.6336,         FL67370
     * 2.6353,   2.6325,   2.6075,   2.5680,   2.5340,   2.5025,         FL67380
     * 2.5122,   2.5652,   2.5681,   2.5869,   2.5925,   2.5986,         FL67390
     * 2.5947,   2.5835,   2.5566,   2.4949/                             FL67400
      DATA C3/                                                           FL67410
     * 2.9825,   2.9831,   2.9847,   2.9893,   2.9929,   2.9976,         FL67420
     * 3.0090,   3.0130,   3.0179,   3.0233,   3.0294,   3.0148,         FL67430
     * 3.0293,   3.0357,   3.0481,   3.0563,   3.0627,   3.0410,         FL67440
     * 3.0532,   3.0646,   3.0713,   3.0733,   3.0716,   3.0701,         FL67450
     * 3.0681,   3.0662,   3.0457,   3.0067,   2.9733,   2.9460,         FL67460
     * 2.9643,   3.0156,   3.0182,   3.0337,   3.0399,   3.0477,         FL67470
     * 3.0511,   3.0501,   3.0384,   2.9943/                             FL67480
      DATA C4/                                                           FL67490
     * 3.3635,   3.3621,   3.3652,   3.3699,   3.3729,   3.3768,         FL67500
     * 3.3868,   3.3888,   3.3916,   3.3952,   3.4000,   3.3911,         FL67510
     * 3.4013,   3.4056,   3.4152,   3.4218,   3.4280,   3.4148,         FL67520
     * 3.4222,   3.4312,   3.4393,   3.4442,   3.4452,   3.4463,         FL67530
     * 3.4455,   3.4452,   3.4329,   3.4016,   3.3719,   3.3468,         FL67540
     * 3.3617,   3.4046,   3.4068,   3.4198,   3.4255,   3.4334,         FL67550
     * 3.4402,   3.4447,   3.4428,   3.4144/                             FL67560
      DATA D1/                                                           FL67570
     *-7.7562,  -7.8498,  -7.8630,  -7.8493,  -7.7889,  -7.5044,         FL67580
     *-7.0058,  -6.3955,  -6.3210,  -6.0026,  -5.4176,  -4.5443,         FL67590
     *-5.6380,  -6.2635,  -5.9512,  -5.9860,  -5.9526,  -5.1907,         FL67600
     *-5.3115,  -5.6289,  -5.6502,  -5.6922,  -5.6157,  -5.5462,         FL67610
     *-5.5437,  -5.5234,  -5.6647,  -5.6087,  -5.5250,  -5.3918,         FL67620
     *-5.1733,  -5.0909,  -5.0907,  -5.1291,  -5.1450,  -5.1968,         FL67630
     *-5.3105,  -5.3997,  -5.4875,  -5.4779/                             FL67640
      DATA D2/                                                           FL67650
     *-7.5869,  -7.6977,  -7.7070,  -7.6883,  -7.6227,  -7.2788,         FL67660
     *-6.6637,  -5.9117,  -6.0351,  -5.6292,  -4.8814,  -3.8947,         FL67670
     *-5.0236,  -5.7607,  -5.3390,  -5.4052,  -5.4335,  -4.5711,         FL67680
     *-4.6910,  -5.0400,  -5.1263,  -5.1522,  -5.1200,  -5.0797,         FL67690
     *-5.0708,  -5.0554,  -5.0883,  -4.9842,  -4.8775,  -4.7313,         FL67700
     *-4.5074,  -4.4271,  -4.4290,  -4.4661,  -4.4923,  -4.5484,         FL67710
     *-4.6679,  -4.7716,  -4.8817,  -4.9075/                             FL67720
      DATA D3/                                                           FL67730
     *-7.3806,  -7.5324,  -7.5421,  -7.5190,  -7.4456,  -6.9683,         FL67740
     *-6.1934,  -5.3374,  -5.6261,  -5.1328,  -4.2936,  -3.2785,         FL67750
     *-4.3895,  -5.1770,  -4.7151,  -4.7944,  -4.8513,  -3.9542,         FL67760
     *-4.0698,  -4.4296,  -4.5444,  -4.5647,  -4.5533,  -4.5320,         FL67770
     *-4.5225,  -4.5111,  -4.4899,  -4.3685,  -4.2548,  -4.1053,         FL67780
     *-3.8800,  -3.7987,  -3.8013,  -3.8369,  -3.8663,  -3.9238,         FL67790
     *-4.0456,  -4.1545,  -4.2728,  -4.3123/                             FL67800
      DATA D4/                                                           FL67810
     *-7.1591,  -7.3911,  -7.3998,  -7.3737,  -7.2891,  -6.6133,         FL67820
     *-5.7137,  -4.8091,  -5.1828,  -4.6408,  -3.7712,  -2.7644,         FL67830
     *-3.8361,  -4.6426,  -4.1724,  -4.2573,  -4.3263,  -3.4249,         FL67840
     *-3.5341,  -3.8962,  -4.0222,  -4.0421,  -4.0386,  -4.0258,         FL67850
     *-4.0169,  -4.0077,  -3.9676,  -3.8419,  -3.7270,  -3.5776,         FL67860
     *-3.3529,  -3.2698,  -3.2724,  -3.3066,  -3.3365,  -3.3940,         FL67870
     *-3.5164,  -3.6272,  -3.7486,  -3.7935/                             FL67880
      DATA E1/                                                           FL67890
     *-4.1531,  -4.2017,  -4.0836,  -4.1441,  -4.0515,  -3.7234,         FL67900
     *-3.2022,  -2.5924,  -2.5215,  -2.2244,  -1.7099,  -1.0243,         FL67910
     *-1.8178,  -2.4304,  -2.1483,  -2.1897,  -2.1768,  -1.5025,         FL67920
     *-1.5770,  -1.8688,  -1.9132,  -1.9550,  -1.9023,  -1.8200,         FL67930
     *-1.8019,  -1.7822,  -1.9415,  -1.9082,  -1.8419,  -1.7290,         FL67940
     *-1.5359,  -1.4523,  -1.4511,  -1.4744,  -1.4875,  -1.5339,         FL67950
     *-1.6446,  -1.7377,  -1.8338,  -1.8404/                             FL67960
      DATA E2/                                                           FL67970
     *-4.0237,  -4.0786,  -4.0596,  -4.0117,  -3.9167,  -3.5334,         FL67980
     *-2.8890,  -2.1314,  -2.2533,  -1.8686,  -1.2114,  -0.5112,         FL67990
     *-1.2226,  -1.9313,  -1.5503,  -1.6190,  -1.6646,  -0.9328,         FL68000
     *-0.9892,  -1.2921,  -1.3909,  -1.4236,  -1.4060,  -1.3666,         FL68010
     *-1.3550,  -1.3429,  -1.3966,  -1.3198,  -1.2346,  -1.1147,         FL68020
     *-0.9248,  -0.8332,  -0.8328,  -0.8490,  -0.8658,  -0.9072,         FL68030
     *-1.0110,  -1.1088,  -1.2210,  -1.2642/                             FL68040
      DATA E3/                                                           FL68050
     *-3.8225,  -3.9189,  -3.8934,  -3.8788,  -3.7792,  -3.2584,         FL68060
     *-2.4500,  -1.5859,  -1.8664,  -1.3920,  -0.6602,  -0.0250,         FL68070
     *-0.6305,  -1.3614,  -0.9442,  -1.0200,  -1.0892,  -0.3681,         FL68080
     *-0.4088,  -0.6976,  -0.8140,  -0.8430,  -0.8410,  -0.8268,         FL68090
     *-0.8209,  -0.8142,  -0.8176,  -0.7305,  -0.6447,  -0.5305,         FL68100
     *-0.3534,  -0.2582,  -0.2574,  -0.2661,  -0.2802,  -0.3137,         FL68110
     *-0.4046,  -0.4954,  -0.6063,  -0.6635/                             FL68120
      DATA E4/                                                           FL68130
     *-3.6380,  -3.8218,  -3.8158,  -3.6544,  -3.6442,  -2.9366,         FL68140
     *-1.9981,  -1.0852,  -1.4468,  -0.9222,  -0.1746,   0.3789,         FL68150
     *-0.1326,  -0.8516,  -0.4270,  -0.5021,  -0.5774,   0.1072,         FL68160
     * 0.0779,  -0.1890,  -0.3060,  -0.3330,  -0.3362,  -0.3320,         FL68170
     *-0.3290,  -0.3248,  -0.3123,  -0.2275,  -0.1469,  -0.0421,         FL68180
     * 0.1192,   0.2136,   0.2149,   0.2122,   0.2019,   0.1760,         FL68190
     * 0.0989,   0.0190,  -0.0836,  -0.1437/                             FL68200
      DATA F1/                                                           FL68210
     *-0.5486,  -0.6082,  -0.5956,  -0.5356,  -0.4402,  -0.0871,         FL68220
     * 0.4527,   1.0366,   1.1096,   1.3655,   1.7101,   1.8903,         FL68230
     * 1.6543,   1.2291,   1.4722,   1.4553,   1.4742,   1.8427,         FL68240
     * 1.8260,   1.6925,   1.6714,   1.6561,   1.6818,   1.7408,         FL68250
     * 1.7604,   1.7735,   1.6870,   1.6975,   1.7266,   1.7732,         FL68260
     * 1.8476,   1.8953,   1.8977,   1.9100,   1.9121,   1.9074,         FL68270
     * 1.8820,   1.8553,   1.8167,   1.8034/                             FL68280
      DATA F2/                                                           FL68290
     *-0.4081,  -0.4784,  -0.4660,  -0.4117,  -0.3046,   0.0831,         FL68300
     * 0.7409,   1.4609,   1.3780,   1.7134,   2.1471,   2.2808,         FL68310
     * 2.1315,   1.6742,   1.9804,   1.9449,   1.9238,   2.2748,         FL68320
     * 2.2689,   2.1587,   2.1154,   2.1037,   2.1124,   2.1387,         FL68330
     * 2.1490,   2.1552,   2.1238,   2.1535,   2.1840,   2.2226,         FL68340
     * 2.2790,   2.3247,   2.3268,   2.3387,   2.3422,   2.3439,         FL68350
     * 2.3339,   2.3198,   2.2926,   2.2751/                             FL68360
      DATA F3/                                                           FL68370
     *-0.2242,  -0.3289,  -0.3406,  -0.2786,  -0.1532,   0.3414,         FL68380
     * 1.1618,   1.9783,   1.7412,   2.1629,   2.6182,   2.6999,         FL68390
     * 2.6101,   2.1844,   2.4931,   2.4589,   2.4253,   2.7204,         FL68400
     * 2.7182,   2.6391,   2.5975,   2.5896,   2.5918,   2.6017,         FL68410
     * 2.6055,   2.6086,   2.6049,   2.6334,   2.6574,   2.6843,         FL68420
     * 2.7225,   2.7601,   2.7620,   2.7734,   2.7780,   2.7832,         FL68430
     * 2.7839,   2.7809,   2.7677,   2.7571/                             FL68440
      DATA F4/                                                           FL68450
     *-0.0119,  -0.2110,  -0.2063,  -0.1444,  -0.0667,   0.6542,         FL68460
     * 1.5923,   2.4405,   2.1326,   2.5924,   3.0247,   3.0696,         FL68470
     * 3.0154,   2.6365,   2.9257,   2.8957,   2.8634,   3.1026,         FL68480
     * 3.1009,   3.0465,   3.0135,   3.0090,   3.0097,   3.0147,         FL68490
     * 3.0167,   3.0193,   3.0233,   3.0464,   3.0635,   3.0808,         FL68500
     * 3.1047,   3.1337,   3.1354,   3.1458,   3.1506,   3.1572,         FL68510
     * 3.1639,   3.1680,   3.1651,   3.1624/                             FL68520
      END                                                                FL68530
C
C     *************************************************************
C
      SUBROUTINE LCONVR (P,T)                                            FL68540
C                                                                        FL68550
C     *************************************************************      FL68560
C                                                                        FL68570
C     WRITTEN APR, 1985 TO ACCOMMODATE 'JCHAR' DEFINITIONS FOR           FL68580
C     UNIFORM DATA INPUT -                                               FL68590
C                                                                        FL68600
C     JCHAR    JUNIT                                                     FL68610
C                                                                        FL68620
C     " ",A       10    VOLUME MIXING RATIO (PPMV)                       FL68630
C     B       11    NUMBER DENSITY (CM-3)                                FL68640
C     C       12    MASS MIXING RATIO (GM(K)/KG(AIR))                    FL68650
C     D       13    MASS DENSITY (GM M-3)                                FL68660
C     E       14    PARTIAL PRESSURE (MB)                                FL68670
C     F       15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY               FL68680
C     G       16     "    "     "  (TD IN T(C)) - H2O ONLY               FL68690
C     H       17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY         FL68700
C     I       18    AVAILABLE FOR USER DEFINITION                        FL68710
C     J       19    REQUEST DEFAULT TO SPECIFIED MODEL ATMOSPHERE        FL68720
C                                                                        FL68730
C     ***************************************************************    FL68740
C                                                                        FL68750
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
      PARAMETER (NCASE=15, NCASE2=NCASE-2)
C
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL68760
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL68770
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL68780
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2,GRAV,CPDAIR,AIRMWT,SECDY 
      COMMON /CONSTL/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL) 
      COMMON /CARD1B/ JUNITP,JUNITT,JUNIT1(NCASE2),WMOL1(NCASE),
     *                WAIR1,JLOW                                         FL68810
C
      DATA C1/18.9766/,C2/-14.9595/,C3/-2.43882/                         FL68820
C
C      DENSAT(ATEMP) = ATEMP*B*EXP(C1+C2*ATEMP+C3*ATEMP**2)*1.0E-6       FL68830
C                                                                        FL68840
      RHOAIR = ALOSMT*(P/PZERO)*(TZERO/T)                                FL68850
C                                                                        FL68860
C     NOPRNT = 0                                                         FL68870
C     A = TZERO/T                                                        FL68880
C                                                                        FL68890
      DO 70 K = 1, 12                                                    FL68900
         B = AVOGAD/AMWT(K)                                              FL68910
         R = AIRMWT/AMWT(K)                                              FL68920
         JUNIT = JUNIT1(K)                                               FL68930
         WMOL = WMOL1(K)                                                 FL68940
         IF (K.NE.1) GO TO 10                                            FL68950
         CALL LWATVA (P,T)                                               FL68960
         GO TO 70                                                        FL68970
   10    CONTINUE                                                        FL68980
         IF (JUNIT.GT.10) GO TO 20                                       FL68990
C                                                                        FL69000
C        **   GIVEN VOL. MIXING RATIO                                    FL69010
C                                                                        FL69020
C        C    WMOL1(K)=WMOL*RHOAIR*1.E-6                                 FL69030
C                                                                        FL69040
         GO TO 70                                                        FL69050
   20    IF (JUNIT.NE.11) GO TO 30                                       FL69060
C                                                                        FL69070
C        **   GIVEN NUMBER DENSITY (CM-3)                                FL69080
C                                                                        FL69090
C        C    WMOL1(K) = WMOL                                            FL69100
C                                                                        FL69110
         WMOL1(K) = WMOL/(RHOAIR*1.E-6)                                  FL69120
         GO TO 70                                                        FL69130
   30    CONTINUE                                                        FL69140
         IF (JUNIT.NE.12) GO TO 40                                       FL69150
C                                                                        FL69160
C        **   GIVEN MASS MIXING RATIO (GM KG-1)                          FL69170
C                                                                        FL69180
C        C    WMOL1(K)= R*WMOL*1.0E-3*RHOAIR                             FL69190
C                                                                        FL69200
         WMOL1(K) = R*WMOL*1.0E+3                                        FL69210
         GO TO 70                                                        FL69220
   40    CONTINUE                                                        FL69230
         IF (JUNIT.NE.13) GO TO 50                                       FL69240
C                                                                        FL69250
C        **   GIVEN MASS DENSITY (GM M-3)                                FL69260
C                                                                        FL69270
C        C    WMOL1(K) = B*WMOL*1.0E-6                                   FL69280
C                                                                        FL69290
         WMOL1(K) = B*WMOL/RHOAIR                                        FL69300
         GO TO 70                                                        FL69310
   50    CONTINUE                                                        FL69320
         IF (JUNIT.NE.14) GO TO 60                                       FL69330
C                                                                        FL69340
C        **   GIVEN  PARTIAL PRESSURE (MB)                               FL69350
C                                                                        FL69360
C        C    WMOL1(K)= ALOSMT*(WMOL/PZERO)*(TZERO/T)                    FL69370
C                                                                        FL69380
         WTEM = ALOSMT*(WMOL/PZERO)*(TZERO/T)                            FL69390
         WMOL1(K) = WTEM/(RHOAIR*1.E-6)                                  FL69400
         GO TO 70                                                        FL69410
   60    CONTINUE                                                        FL69420
         IF (JUNIT.GT.14) GO TO 80                                       FL69430
   70 CONTINUE                                                           FL69440
      RETURN                                                             FL69450
   80 CONTINUE                                                           FL69460
      WRITE (IPR,900) JUNIT                                              FL69470
      STOP                                                               FL69480
C                                                                        FL69490
  900 FORMAT(/,'   **** ERROR IN CONVERT ****, JUNIT = ',I5)             FL69500
C                                                                        FL69510
      END                                                                FL69520
C
C     *************************************************************
C
      SUBROUTINE LWATVA(P,T)                                             FL69530
C                                                                        FL69540
C     *************************************************************      FL69550
C                                                                        FL69560
C     WRITTEN APR, 1985 TO ACCOMMODATE 'JCHAR' DEFINITIONS FOR           FL69570
C     UNIFORM DATA INPUT -                                               FL69580
C                                                                        FL69590
C     JCHAR    JUNIT                                                     FL69600
C                                                                        FL69610
C     " ",A       10    VOLUME MIXING RATIO (PPMV)                       FL69620
C     B       11    NUMBER DENSITY (CM-3)                                FL69630
C     C       12    MASS MIXING RATIO (GM(K)/KG(AIR))                    FL69640
C     D       13    MASS DENSITY (GM M-3)                                FL69650
C     E       14    PARTIAL PRESSURE (MB)                                FL69660
C     F       15    DEW POINT TEMP (TD IN T(K)) - H2O ONLY               FL69670
C     G       16     "    "     "  (TD IN T(C)) - H2O ONLY               FL69680
C     H       17    RELATIVE HUMIDITY (RH IN PERCENT) - H2O ONLY         FL69690
C     I       18    AVAILABLE FOR USER DEFINITION                        FL69700
C     J       19    REQUEST DEFAULT TO SPECIFIED MODEL ATMOSPHERE        FL69710
C                                                                        FL69720
C     THIS SUBROUTINE COMPUTES THE WATERVAPOR NUMBER DENSITY (MOL CM-3)  FL69730
C     GIVE HUMIDITY  # TD = DEW POINT TEMP(K,C), RH = RELATIVE           FL69740
C     (PERCENT), PPH2O = WATER VAPOR PARTIAL PRESSURE (MB), DENH2O =     FL69750
C     WATER VAPOR MASS DENSITY (GM M-3),AMSMIX = MASS MIXING RATIO       FL69760
C     (GM/KG).                                                           FL69770
C     THE FUNCTION DENSAT FOR THE SATURATION                             FL69780
C     WATER VAPOR DENSITY OVER WATER IS ACCURATE TO BETTER THAN 1        FL69790
C     PERCENT FROM -50 TO +50 DEG C. (SEE THE LOWTRAN3 OR 5 REPORT)      FL69800
C                                                                        FL69810
C     'JUNIT' GOVERNS CHOICE OF UNITS -                                  FL69820
C                                                                        FL69830
C     ****************************************************************** FL69840
C                                                                        FL69850
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
      PARAMETER (NCASE=15, NCASE2=NCASE-2)
C
      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL69860
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL69870
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL69880
      COMMON /CARD1B/ JUNITP,JUNITT,JUNIT1(NCASE2),WMOL1(NCASE),
     *                WAIR,JLOW                                          FL69890
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2,GRAV,CPDAIR,AIRMWT,SECDY 
      COMMON /CONSTL/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL) 
C
      DATA C1/18.9766/,C2/-14.9595/,C3/-2.43882/                         FL69920
      DATA XLOSCH/2.6868E19/                                             FL69930
C
      DENSAT(ATEMP) = ATEMP*B*EXP(C1+C2*ATEMP+C3*ATEMP**2)*1.0E-6        FL69940
C                                                                        FL69950
      RHOAIR = ALOSMT*(P/PZERO)*(TZERO/T)                                FL69960
      PSS = P/PZERO                                                      FL69970
      A = TZERO/T                                                        FL69980
      WAIR = XLOSCH*PSS*A                                                FL69990
      B = AVOGAD/AMWT(1)                                                 FL70000
      R = AIRMWT/AMWT(1)                                                 FL70010
      JUNIT = JUNIT1(1)                                                  FL70020
      WMOL = WMOL1(1)                                                    FL70030
      IF (JUNIT.NE.10) GO TO 10                                          FL70040
C                                                                        FL70050
C     **   GIVEN VOL. MIXING RATIO                                       FL70060
C                                                                        FL70070
C     C    WMOL1(1)=WMOL*RHOAIR*1.E-6                                    FL70080
C                                                                        FL70090
      GO TO 90                                                           FL70100
   10 IF (JUNIT.NE.11) GO TO 20                                          FL70110
C                                                                        FL70120
C     **   GIVEN NUMBER DENSITY (CM-3)                                   FL70130
C                                                                        FL70140
      WMOL1(1) = WMOL/(RHOAIR*1.E-6)                                     FL70150
      GO TO 90                                                           FL70160
   20 CONTINUE                                                           FL70170
      IF (JUNIT.NE.12) GO TO 30                                          FL70180
C                                                                        FL70190
C     **   GIVEN MASS MIXING RATIO (GM KG-1)                             FL70200
C                                                                        FL70210
C     C    WMOL1(1) = R*WMOL*1.0E-3*RHOAIR                               FL70220
C                                                                        FL70230
      WMOL1(1) = R*WMOL*1.0E+3                                           FL70240
      GO TO 90                                                           FL70250
   30 CONTINUE                                                           FL70260
      IF (JUNIT.NE.13) GO TO 40                                          FL70270
C                                                                        FL70280
C     **   GIVEN MASS DENSITY (GM M-3)                                   FL70290
C                                                                        FL70300
C     C    WMOL1(1) = B*WMOL*1.0E-6                                      FL70310
C                                                                        FL70320
      WMOL1(1) = B*WMOL/RHOAIR                                           FL70330
      GO TO 90                                                           FL70340
   40 CONTINUE                                                           FL70350
      IF (JUNIT.NE.14) GO TO 50                                          FL70360
C                                                                        FL70370
C     **   GIVEN WATER VAPOR PARTIAL PRESSURE (MB)                       FL70380
C                                                                        FL70390
C     C    WMOL1(1) = ALOSMT*(WMOL/PZERO)*(TZERO/T)                      FL70400
C                                                                        FL70410
      WTEM = ALOSMT*(WMOL/PZERO)*(TZERO/T)                               FL70420
      WMOL1(1) = WTEM/(RHOAIR*1.E-6)                                     FL70430
      GO TO 90                                                           FL70440
   50 CONTINUE                                                           FL70450
      IF (JUNIT.NE.15) GO TO 60                                          FL70460
C                                                                        FL70470
C     **   GIVEN DEWPOINT (DEG K)                                        FL70480
C                                                                        FL70490
      ATD = TZERO/(WMOL)                                                 FL70500
C                                                                        FL70510
C     C    WMOL1(1)= DENSAT(ATD)*(WMOL)/T                                FL70520
C                                                                        FL70530
      WTEM = DENSAT(ATD)*(WMOL)/T                                        FL70540
      WMOL1(1) = WTEM/(RHOAIR*1.E-6)                                     FL70550
      GO TO 90                                                           FL70560
   60 CONTINUE                                                           FL70570
      IF (JUNIT.NE.16) GO TO 70                                          FL70580
C                                                                        FL70590
C     **   GIVEN DEWPOINT (DEG C)                                        FL70600
C                                                                        FL70610
      ATD = TZERO/(TZERO+WMOL)                                           FL70620
C                                                                        FL70630
C     C    WMOL1(1) = DENSAT(ATD)*(TZERO+WMOL)/T                         FL70640
C                                                                        FL70650
      WTEM = DENSAT(ATD)*(TZERO+WMOL)/T                                  FL70660
      WMOL1(1) = WTEM/(RHOAIR*1.E-6)                                     FL70670
      GO TO 90                                                           FL70680
   70 CONTINUE                                                           FL70690
      IF (JUNIT.NE.17) GO TO 80                                          FL70700
C                                                                        FL70710
C     **   GIVEN RELATIVE HUMIDITY (PERCENT)                             FL70720
C                                                                        FL70730
C     DENNUM = DENSAT(A)*(WMOL/100.0)/(1.0-(1.0-WMOL/100.0)*DENSAT(A)/   FL70740
C     1    RHOAIR)                                                       FL70750
C     C    WMOL1(1) = DENSAT(A)*(WMOL/100.0)                             FL70760
C                                                                        FL70770
      WTEM = DENSAT(A)*(WMOL/100.0)                                      FL70780
      WMOL1(1) = WTEM/(RHOAIR*1.E-6)                                     FL70790
      GO TO 90                                                           FL70800
   80 WRITE (IPR,900) JUNIT                                              FL70810
      STOP 'JUNIT'                                                       FL70820
   90 CONTINUE                                                           FL70830
      WMOL1(1) = 2.989E-23*WMOL1(1)*WAIR                                 FL70840
      DENST = DENSAT(A)                                                  FL70850
      DENNUM = WMOL1(1)                                                  FL70860
C                                                                        FL70870
C     RHP = 100.0*(DENNUM/DENST)*((RHOAIR-DENST)/(RHOAIR-DENNUM))        FL70880
C                                                                        FL70890
      RHP = 100.0*(DENNUM/DENST)                                         FL70900
      IF (RHP.LE.100.0) GO TO 100                                        FL70910
      WRITE (IPR,905) RHP                                                FL70920
  100 CONTINUE                                                           FL70930
      RETURN                                                             FL70940
C                                                                        FL70950
  900 FORMAT(/,'  **** ERROR IN WATVAP ****, JUNIT = ',I5)               FL70960
  905 FORMAT(/,' ********WARNING (FROM WATVAP) # RELATIVE HUMIDTY = ',   FL70970
     *    G10.3,' IS GREATER THAN 100 PERCENT')                          FL70980
C                                                                        FL70990
      END                                                                FL71000
      BLOCK DATA MLATLB                                                  FL71010
C                                                                        FL71020
C     ****************************************************************** FL71030
C     THIS SUBROUTINE INITIALIZES THE 6 BUILT-IN ATMOSPHERIC PROFILES    FL71040
C     (FROM 'OPTICAL PROPERTIES OF THE ATMOSPHERE, THIRD EDITION'        FL71050
C     AFCRL-72-0497 (AD 753 075), 'U.S. STANDARD ATMOSPHERE 1976' AND    FL71060
C     'SUPPLEMENTS 1966'), PLUS COLLECTED CONSTITUENT PROFILES (REF)     FL71070
C     AND SETS OTHER CONSTANTS RELATED TO THE ATMOSPHERIC PROFILES       FL71080
C     ****************************************************************** FL71090
C                                                                        FL71100
      parameter (mxmol=39)

      CHARACTER*8 CTMNA1,CTMNA2,CTMNA3,CTMNA4,CTMNA5,CTMNA6              FL71110
      COMMON /CLATML/                                                    FL71120
     *CTMNA1(3),CTMNA2(3),CTMNA3(3),CTMNA4(3),CTMNA5(3),CTMNA6(3)        FL71130
      REAL*8           ATMNA1,ATMNA2,ATMNA3,ATMNA4,ATMNA5,ATMNA6
      Character*8                                                HMODS   FL71140
      COMMON /MLATML/ ALT(50),P1(50),P2(50),P3(50),P4(50),P5(50),P6(50), FL71150
     *T1(50),T2(50),T3(50),T4(50),T5(50),T6(50),                         FL71160
     *AMOL11(50),AMOL12(50),AMOL13(50),AMOL14(50),AMOL15(50),AMOL16(50), FL71170
     *AMOL17(50),AMOL18(50),                                             FL71180
     *AMOL21(50),AMOL22(50),AMOL23(50),AMOL24(50),AMOL25(50),AMOL26(50), FL71190
     *AMOL27(50),AMOL28(50),                                             FL71200
     *AMOL31(50),AMOL32(50),AMOL33(50),AMOL34(50),AMOL35(50),AMOL36(50), FL71210
     *AMOL37(50),AMOL38(50),                                             FL71220
     *AMOL41(50),AMOL42(50),AMOL43(50),AMOL44(50),AMOL45(50),AMOL46(50), FL71230
     *AMOL47(50),AMOL48(50),                                             FL71240
     *AMOL51(50),AMOL52(50),AMOL53(50),AMOL54(50),AMOL55(50),AMOL56(50), FL71250
     *AMOL57(50),AMOL58(50),                                             FL71260
     *AMOL61(50),AMOL62(50),AMOL63(50),AMOL64(50),AMOL65(50),AMOL66(50), FL71270
     *AMOL67(50),AMOL68(50),                                             FL71280
     *ATMNA1(3),ATMNA2(3),ATMNA3(3),ATMNA4(3),ATMNA5(3),ATMNA6(3),       FL71290
     *     HMODS(3),ZST(50),PST(50),TST(50),AMOLS(50,mxmol),IDUM            FL71300
C                                                                        FL71320
C     COMMON /TRACL/ TRAC(50,22)                                         FL71330
C                                                                        FL71340
      COMMON /TRACL/ ANO(50),SO2(50),ANO2(50),ANH3(50),HNO3(50),OH(50),  FL71350
     * HF(50),HCL(50),HBR(50),HI(50),CLO(50),OCS(50),H2CO(50),           FL71360
     * HOCL(50),AN2(50),HCN(50),CH3CL(50),H2O2(50),C2H2(50),             FL71370
     * C2H6(50),PH3(50),TDUM(50)                                         FL71380
      DATA CTMNA1  /'TROPICAL','        ','        '/                    FL71390
      DATA CTMNA2  /'MIDLATIT','UDE SUMM','ER      '/                    FL71400
      DATA CTMNA3  /'MIDLATIT','UDE WINT','ER      '/                    FL71410
      DATA CTMNA4  /'SUBARCTI','C SUMMER','        '/                    FL71420
      DATA CTMNA5  /'SUBARCTI','C WINTER','        '/                    FL71430
      DATA CTMNA6  /'U. S. ST','ANDARD, ','1976    '/                    FL71440
C                                                                        FL71450
C     DATA ALT (KM)  /                                                   FL71460
C                                                                        FL71470
      DATA ALT/                                                          FL71480
     *       0.0,       1.0,       2.0,       3.0,       4.0,            FL71490
     *       5.0,       6.0,       7.0,       8.0,       9.0,            FL71500
     *      10.0,      11.0,      12.0,      13.0,      14.0,            FL71510
     *      15.0,      16.0,      17.0,      18.0,      19.0,            FL71520
     *      20.0,      21.0,      22.0,      23.0,      24.0,            FL71530
     *      25.0,      27.5,      30.0,      32.5,      35.0,            FL71540
     *      37.5,      40.0,      42.5,      45.0,      47.5,            FL71550
     *      50.0,      55.0,      60.0,      65.0,      70.0,            FL71560
     *      75.0,      80.0,      85.0,      90.0,      95.0,            FL71570
     *     100.0,     105.0,     110.0,     115.0,     120.0/            FL71580
C                                                                        FL71590
C     DATA PRESSURE  /                                                   FL71600
C                                                                        FL71610
      DATA P1/                                                           FL71620
     * 1.013E+03, 9.040E+02, 8.050E+02, 7.150E+02, 6.330E+02,            FL71630
     * 5.590E+02, 4.920E+02, 4.320E+02, 3.780E+02, 3.290E+02,            FL71640
     * 2.860E+02, 2.470E+02, 2.130E+02, 1.820E+02, 1.560E+02,            FL71650
     * 1.320E+02, 1.110E+02, 9.370E+01, 7.890E+01, 6.660E+01,            FL71660
     * 5.650E+01, 4.800E+01, 4.090E+01, 3.500E+01, 3.000E+01,            FL71670
     * 2.570E+01, 1.763E+01, 1.220E+01, 8.520E+00, 6.000E+00,            FL71680
     * 4.260E+00, 3.050E+00, 2.200E+00, 1.590E+00, 1.160E+00,            FL71690
     * 8.540E-01, 4.560E-01, 2.390E-01, 1.210E-01, 5.800E-02,            FL71700
     * 2.600E-02, 1.100E-02, 4.400E-03, 1.720E-03, 6.880E-04,            FL71710
     * 2.890E-04, 1.300E-04, 6.470E-05, 3.600E-05, 2.250E-05/            FL71720
      DATA P2/                                                           FL71730
     * 1.013E+03, 9.020E+02, 8.020E+02, 7.100E+02, 6.280E+02,            FL71740
     * 5.540E+02, 4.870E+02, 4.260E+02, 3.720E+02, 3.240E+02,            FL71750
     * 2.810E+02, 2.430E+02, 2.090E+02, 1.790E+02, 1.530E+02,            FL71760
     * 1.300E+02, 1.110E+02, 9.500E+01, 8.120E+01, 6.950E+01,            FL71770
     * 5.950E+01, 5.100E+01, 4.370E+01, 3.760E+01, 3.220E+01,            FL71780
     * 2.770E+01, 1.907E+01, 1.320E+01, 9.300E+00, 6.520E+00,            FL71790
     * 4.640E+00, 3.330E+00, 2.410E+00, 1.760E+00, 1.290E+00,            FL71800
     * 9.510E-01, 5.150E-01, 2.720E-01, 1.390E-01, 6.700E-02,            FL71810
     * 3.000E-02, 1.200E-02, 4.480E-03, 1.640E-03, 6.250E-04,            FL71820
     * 2.580E-04, 1.170E-04, 6.110E-05, 3.560E-05, 2.270E-05/            FL71830
      DATA P3/                                                           FL71840
     * 1.018E+03, 8.973E+02, 7.897E+02, 6.938E+02, 6.081E+02,            FL71850
     * 5.313E+02, 4.627E+02, 4.016E+02, 3.473E+02, 2.993E+02,            FL71860
     * 2.568E+02, 2.199E+02, 1.882E+02, 1.611E+02, 1.378E+02,            FL71870
     * 1.178E+02, 1.007E+02, 8.610E+01, 7.360E+01, 6.280E+01,            FL71880
     * 5.370E+01, 4.580E+01, 3.910E+01, 3.340E+01, 2.860E+01,            FL71890
     * 2.440E+01, 1.646E+01, 1.110E+01, 7.560E+00, 5.180E+00,            FL71900
     * 3.600E+00, 2.530E+00, 1.800E+00, 1.290E+00, 9.400E-01,            FL71910
     * 6.830E-01, 3.620E-01, 1.880E-01, 9.500E-02, 4.700E-02,            FL71920
     * 2.220E-02, 1.030E-02, 4.560E-03, 1.980E-03, 8.770E-04,            FL71930
     * 4.074E-04, 2.000E-04, 1.057E-04, 5.980E-05, 3.600E-05/            FL71940
      DATA P4/                                                           FL71950
     * 1.010E+03, 8.960E+02, 7.929E+02, 7.000E+02, 6.160E+02,            FL71960
     * 5.410E+02, 4.740E+02, 4.130E+02, 3.590E+02, 3.108E+02,            FL71970
     * 2.677E+02, 2.300E+02, 1.977E+02, 1.700E+02, 1.460E+02,            FL71980
     * 1.260E+02, 1.080E+02, 9.280E+01, 7.980E+01, 6.860E+01,            FL71990
     * 5.900E+01, 5.070E+01, 4.360E+01, 3.750E+01, 3.228E+01,            FL72000
     * 2.780E+01, 1.923E+01, 1.340E+01, 9.400E+00, 6.610E+00,            FL72010
     * 4.720E+00, 3.400E+00, 2.480E+00, 1.820E+00, 1.340E+00,            FL72020
     * 9.870E-01, 5.370E-01, 2.880E-01, 1.470E-01, 7.100E-02,            FL72030
     * 3.200E-02, 1.250E-02, 4.510E-03, 1.610E-03, 6.060E-04,            FL72040
     * 2.480E-04, 1.130E-04, 6.000E-05, 3.540E-05, 2.260E-05/            FL72050
      DATA P5/                                                           FL72060
     * 1.013E+03, 8.878E+02, 7.775E+02, 6.798E+02, 5.932E+02,            FL72070
     * 5.158E+02, 4.467E+02, 3.853E+02, 3.308E+02, 2.829E+02,            FL72080
     * 2.418E+02, 2.067E+02, 1.766E+02, 1.510E+02, 1.291E+02,            FL72090
     * 1.103E+02, 9.431E+01, 8.058E+01, 6.882E+01, 5.875E+01,            FL72100
     * 5.014E+01, 4.277E+01, 3.647E+01, 3.109E+01, 2.649E+01,            FL72110
     * 2.256E+01, 1.513E+01, 1.020E+01, 6.910E+00, 4.701E+00,            FL72120
     * 3.230E+00, 2.243E+00, 1.570E+00, 1.113E+00, 7.900E-01,            FL72130
     * 5.719E-01, 2.990E-01, 1.550E-01, 7.900E-02, 4.000E-02,            FL72140
     * 2.000E-02, 9.660E-03, 4.500E-03, 2.022E-03, 9.070E-04,            FL72150
     * 4.230E-04, 2.070E-04, 1.080E-04, 6.000E-05, 3.590E-05/            FL72160
      DATA P6/                                                           FL72170
     * 1.013E+03, 8.988E+02, 7.950E+02, 7.012E+02, 6.166E+02,            FL72180
     * 5.405E+02, 4.722E+02, 4.111E+02, 3.565E+02, 3.080E+02,            FL72190
     * 2.650E+02, 2.270E+02, 1.940E+02, 1.658E+02, 1.417E+02,            FL72200
     * 1.211E+02, 1.035E+02, 8.850E+01, 7.565E+01, 6.467E+01,            FL72210
     * 5.529E+01, 4.729E+01, 4.047E+01, 3.467E+01, 2.972E+01,            FL72220
     * 2.549E+01, 1.743E+01, 1.197E+01, 8.010E+00, 5.746E+00,            FL72230
     * 4.150E+00, 2.871E+00, 2.060E+00, 1.491E+00, 1.090E+00,            FL72240
     * 7.978E-01, 4.250E-01, 2.190E-01, 1.090E-01, 5.220E-02,            FL72250
     * 2.400E-02, 1.050E-02, 4.460E-03, 1.840E-03, 7.600E-04,            FL72260
     * 3.200E-04, 1.450E-04, 7.100E-05, 4.010E-05, 2.540E-05/            FL72270
C                                                                        FL72280
C     DATA TEMPERATUR/                                                   FL72290
C                                                                        FL72300
      DATA T1/                                                           FL72310
     *    299.70,    293.70,    287.70,    283.70,    277.00,            FL72320
     *    270.30,    263.60,    257.00,    250.30,    243.60,            FL72330
     *    237.00,    230.10,    223.60,    217.00,    210.30,            FL72340
     *    203.70,    197.00,    194.80,    198.80,    202.70,            FL72350
     *    206.70,    210.70,    214.60,    217.00,    219.20,            FL72360
     *    221.40,    227.00,    232.30,    237.70,    243.10,            FL72370
     *    248.50,    254.00,    259.40,    264.80,    269.60,            FL72380
     *    270.20,    263.40,    253.10,    236.00,    218.90,            FL72390
     *    201.80,    184.80,    177.10,    177.00,    184.30,            FL72400
     *    190.70,    212.00,    241.60,    299.70,    380.00/            FL72410
      DATA T2/                                                           FL72420
     *    294.20,    289.70,    285.20,    279.20,    273.20,            FL72430
     *    267.20,    261.20,    254.70,    248.20,    241.70,            FL72440
     *    235.30,    228.80,    222.30,    215.80,    215.70,            FL72450
     *    215.70,    215.70,    215.70,    216.80,    217.90,            FL72460
     *    219.20,    220.40,    221.60,    222.80,    223.90,            FL72470
     *    225.10,    228.45,    233.70,    239.00,    245.20,            FL72480
     *    251.30,    257.50,    263.70,    269.90,    275.20,            FL72490
     *    275.70,    269.30,    257.10,    240.10,    218.10,            FL72500
     *    196.10,    174.10,    165.10,    165.00,    178.30,            FL72510
     *    190.50,    222.20,    262.40,    316.80,    380.00/            FL72520
      DATA T3/                                                           FL72530
     *    272.20,    268.70,    265.20,    261.70,    255.70,            FL72540
     *    249.70,    243.70,    237.70,    231.70,    225.70,            FL72550
     *    219.70,    219.20,    218.70,    218.20,    217.70,            FL72560
     *    217.20,    216.70,    216.20,    215.70,    215.20,            FL72570
     *    215.20,    215.20,    215.20,    215.20,    215.20,            FL72580
     *    215.20,    215.50,    217.40,    220.40,    227.90,            FL72590
     *    235.50,    243.20,    250.80,    258.50,    265.10,            FL72600
     *    265.70,    260.60,    250.80,    240.90,    230.70,            FL72610
     *    220.40,    210.10,    199.80,    199.50,    208.30,            FL72620
     *    218.60,    237.10,    259.50,    293.00,    333.00/            FL72630
      DATA T4/                                                           FL72640
     *    287.20,    281.70,    276.30,    270.90,    265.50,            FL72650
     *    260.10,    253.10,    246.10,    239.20,    232.20,            FL72660
     *    225.20,    225.20,    225.20,    225.20,    225.20,            FL72670
     *    225.20,    225.20,    225.20,    225.20,    225.20,            FL72680
     *    225.20,    225.20,    225.20,    225.20,    226.60,            FL72690
     *    228.10,    231.00,    235.10,    240.00,    247.20,            FL72700
     *    254.60,    262.10,    269.50,    273.60,    276.20,            FL72710
     *    277.20,    274.00,    262.70,    239.70,    216.60,            FL72720
     *    193.60,    170.60,    161.70,    161.60,    176.80,            FL72730
     *    190.40,    226.00,    270.10,    322.70,    380.00/            FL72740
      DATA T5/                                                           FL72750
     *    257.20,    259.10,    255.90,    252.70,    247.70,            FL72760
     *    240.90,    234.10,    227.30,    220.60,    217.20,            FL72770
     *    217.20,    217.20,    217.20,    217.20,    217.20,            FL72780
     *    217.20,    216.60,    216.00,    215.40,    214.80,            FL72790
     *    214.20,    213.60,    213.00,    212.40,    211.80,            FL72800
     *    211.20,    213.60,    216.00,    218.50,    222.30,            FL72810
     *    228.50,    234.70,    240.80,    247.00,    253.20,            FL72820
     *    259.30,    259.10,    250.90,    248.40,    245.40,            FL72830
     *    234.70,    223.90,    213.10,    202.30,    211.00,            FL72840
     *    218.50,    234.00,    252.60,    288.50,    333.00/            FL72850
      DATA T6/                                                           FL72860
     *    288.20,    281.70,    275.20,    268.70,    262.20,            FL72870
     *    255.70,    249.20,    242.70,    236.20,    229.70,            FL72880
     *    223.30,    216.80,    216.70,    216.70,    216.70,            FL72890
     *    216.70,    216.70,    216.70,    216.70,    216.70,            FL72900
     *    216.70,    217.60,    218.60,    219.60,    220.60,            FL72910
     *    221.60,    224.00,    226.50,    230.00,    236.50,            FL72920
     *    242.90,    250.40,    257.30,    264.20,    270.60,            FL72930
     *    270.70,    260.80,    247.00,    233.30,    219.60,            FL72940
     *    208.40,    198.60,    188.90,    186.90,    188.40,            FL72950
     *    195.10,    208.80,    240.00,    300.00,    360.00/            FL72960
C                                                                        FL72970
C     DATA  H2O      /                                                   FL72980
C                                                                        FL72990
      DATA AMOL11/                                                       FL73000
     * 2.593E+04, 1.949E+04, 1.534E+04, 8.600E+03, 4.441E+03,            FL73010
     * 3.346E+03, 2.101E+03, 1.289E+03, 7.637E+02, 4.098E+02,            FL73020
     * 1.912E+02, 7.306E+01, 2.905E+01, 9.900E+00, 6.220E+00,            FL73030
     * 4.000E+00, 3.000E+00, 2.900E+00, 2.750E+00, 2.600E+00,            FL73040
     * 2.600E+00, 2.650E+00, 2.800E+00, 2.900E+00, 3.200E+00,            FL73050
     * 3.250E+00, 3.600E+00, 4.000E+00, 4.300E+00, 4.600E+00,            FL73060
     * 4.900E+00, 5.200E+00, 5.500E+00, 5.700E+00, 5.900E+00,            FL73070
     * 6.000E+00, 6.000E+00, 6.000E+00, 5.400E+00, 4.500E+00,            FL73080
     * 3.300E+00, 2.100E+00, 1.300E+00, 8.500E-01, 5.400E-01,            FL73090
     * 4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01/            FL73100
      DATA AMOL21/                                                       FL73110
     * 1.876E+04, 1.378E+04, 9.680E+03, 5.984E+03, 3.813E+03,            FL73120
     * 2.225E+03, 1.510E+03, 1.020E+03, 6.464E+02, 4.129E+02,            FL73130
     * 2.472E+02, 9.556E+01, 2.944E+01, 8.000E+00, 5.000E+00,            FL73140
     * 3.400E+00, 3.300E+00, 3.200E+00, 3.150E+00, 3.200E+00,            FL73150
     * 3.300E+00, 3.450E+00, 3.600E+00, 3.850E+00, 4.000E+00,            FL73160
     * 4.200E+00, 4.450E+00, 4.700E+00, 4.850E+00, 4.950E+00,            FL73170
     * 5.000E+00, 5.100E+00, 5.300E+00, 5.450E+00, 5.500E+00,            FL73180
     * 5.500E+00, 5.350E+00, 5.000E+00, 4.400E+00, 3.700E+00,            FL73190
     * 2.950E+00, 2.100E+00, 1.330E+00, 8.500E-01, 5.400E-01,            FL73200
     * 4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01/            FL73210
      DATA AMOL31/                                                       FL73220
     * 4.316E+03, 3.454E+03, 2.788E+03, 2.088E+03, 1.280E+03,            FL73230
     * 8.241E+02, 5.103E+02, 2.321E+02, 1.077E+02, 5.566E+01,            FL73240
     * 2.960E+01, 1.000E+01, 6.000E+00, 5.000E+00, 4.800E+00,            FL73250
     * 4.700E+00, 4.600E+00, 4.500E+00, 4.500E+00, 4.500E+00,            FL73260
     * 4.500E+00, 4.500E+00, 4.530E+00, 4.550E+00, 4.600E+00,            FL73270
     * 4.650E+00, 4.700E+00, 4.750E+00, 4.800E+00, 4.850E+00,            FL73280
     * 4.900E+00, 4.950E+00, 5.000E+00, 5.000E+00, 5.000E+00,            FL73290
     * 4.950E+00, 4.850E+00, 4.500E+00, 4.000E+00, 3.300E+00,            FL73300
     * 2.700E+00, 2.000E+00, 1.330E+00, 8.500E-01, 5.400E-01,            FL73310
     * 4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01/            FL73320
      DATA AMOL41/                                                       FL73330
     * 1.194E+04, 8.701E+03, 6.750E+03, 4.820E+03, 3.380E+03,            FL73340
     * 2.218E+03, 1.330E+03, 7.971E+02, 3.996E+02, 1.300E+02,            FL73350
     * 4.240E+01, 1.330E+01, 6.000E+00, 4.450E+00, 4.000E+00,            FL73360
     * 4.000E+00, 4.000E+00, 4.050E+00, 4.300E+00, 4.500E+00,            FL73370
     * 4.600E+00, 4.700E+00, 4.800E+00, 4.830E+00, 4.850E+00,            FL73380
     * 4.900E+00, 4.950E+00, 5.000E+00, 5.000E+00, 5.000E+00,            FL73390
     * 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00,            FL73400
     * 4.950E+00, 4.850E+00, 4.500E+00, 4.000E+00, 3.300E+00,            FL73410
     * 2.700E+00, 2.000E+00, 1.330E+00, 8.500E-01, 5.400E-01,            FL73420
     * 4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01/            FL73430
      DATA AMOL51/                                                       FL73440
     * 1.405E+03, 1.615E+03, 1.427E+03, 1.166E+03, 7.898E+02,            FL73450
     * 4.309E+02, 2.369E+02, 1.470E+02, 3.384E+01, 2.976E+01,            FL73460
     * 2.000E+01, 1.000E+01, 6.000E+00, 4.450E+00, 4.500E+00,            FL73470
     * 4.550E+00, 4.600E+00, 4.650E+00, 4.700E+00, 4.750E+00,            FL73480
     * 4.800E+00, 4.850E+00, 4.900E+00, 4.950E+00, 5.000E+00,            FL73490
     * 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00,            FL73500
     * 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00,            FL73510
     * 4.950E+00, 4.850E+00, 4.500E+00, 4.000E+00, 3.300E+00,            FL73520
     * 2.700E+00, 2.000E+00, 1.330E+00, 8.500E-01, 5.400E-01,            FL73530
     * 4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01/            FL73540
      DATA AMOL61/                                                       FL73550
     * 7.745E+03, 6.071E+03, 4.631E+03, 3.182E+03, 2.158E+03,            FL73560
     * 1.397E+03, 9.254E+02, 5.720E+02, 3.667E+02, 1.583E+02,            FL73570
     * 6.996E+01, 3.613E+01, 1.906E+01, 1.085E+01, 5.927E+00,            FL73580
     * 5.000E+00, 3.950E+00, 3.850E+00, 3.825E+00, 3.850E+00,            FL73590
     * 3.900E+00, 3.975E+00, 4.065E+00, 4.200E+00, 4.300E+00,            FL73600
     * 4.425E+00, 4.575E+00, 4.725E+00, 4.825E+00, 4.900E+00,            FL73610
     * 4.950E+00, 5.025E+00, 5.150E+00, 5.225E+00, 5.250E+00,            FL73620
     * 5.225E+00, 5.100E+00, 4.750E+00, 4.200E+00, 3.500E+00,            FL73630
     * 2.825E+00, 2.050E+00, 1.330E+00, 8.500E-01, 5.400E-01,            FL73640
     * 4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01/            FL73650
C                                                                        FL73660
C     DATA CO2       /                                                   FL73670
C                                                                        FL73680
      DATA AMOL12/                                                       FL73690
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73700
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73710
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73720
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73730
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73740
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73750
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73760
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73770
     * 3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,            FL73780
     * 1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01/            FL73790
C                                                                        FL73800
C     DATA CO2       /                                                   FL73810
C                                                                        FL73820
      DATA AMOL22/                                                       FL73830
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73840
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73850
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73860
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73870
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73880
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73890
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73900
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73910
     * 3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,            FL73920
     * 1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01/            FL73930
C                                                                        FL73940
C     DATA CO2       /                                                   FL73950
C                                                                        FL73960
      DATA AMOL32/                                                       FL73970
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73980
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL73990
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74000
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74010
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74020
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74030
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74040
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74050
     * 3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,            FL74060
     * 1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01/            FL74070
C                                                                        FL74080
C     DATA CO2       /                                                   FL74090
C                                                                        FL74100
      DATA AMOL42/                                                       FL74110
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74120
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74130
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74140
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74150
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74160
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74170
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74180
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74190
     * 3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,            FL74200
     * 1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01/            FL74210
C                                                                        FL74220
C     DATA CO2       /                                                   FL74230
C                                                                        FL74240
      DATA AMOL52/                                                       FL74250
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74260
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74270
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74280
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74290
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74300
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74310
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74320
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74330
     * 3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,            FL74340
     * 1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01/            FL74350
C                                                                        FL74360
C     DATA CO2       /                                                   FL74370
C                                                                        FL74380
      DATA AMOL62/                                                       FL74390
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74400
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74410
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74420
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74430
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74440
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74450
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74460
     * 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,            FL74470
     * 3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,            FL74480
     * 1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01/            FL74490
C                                                                        FL74500
C     DATA OZONE     /                                                   FL74510
C                                                                        FL74520
      DATA AMOL13/                                                       FL74530
     * 2.869E-02, 3.150E-02, 3.342E-02, 3.504E-02, 3.561E-02,            FL74540
     * 3.767E-02, 3.989E-02, 4.223E-02, 4.471E-02, 5.000E-02,            FL74550
     * 5.595E-02, 6.613E-02, 7.815E-02, 9.289E-02, 1.050E-01,            FL74560
     * 1.256E-01, 1.444E-01, 2.500E-01, 5.000E-01, 9.500E-01,            FL74570
     * 1.400E+00, 1.800E+00, 2.400E+00, 3.400E+00, 4.300E+00,            FL74580
     * 5.400E+00, 7.800E+00, 9.300E+00, 9.850E+00, 9.700E+00,            FL74590
     * 8.800E+00, 7.500E+00, 5.900E+00, 4.500E+00, 3.450E+00,            FL74600
     * 2.800E+00, 1.800E+00, 1.100E+00, 6.500E-01, 3.000E-01,            FL74610
     * 1.800E-01, 3.300E-01, 5.000E-01, 5.200E-01, 5.000E-01,            FL74620
     * 4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04/            FL74630
      DATA AMOL23/                                                       FL74640
     * 3.017E-02, 3.337E-02, 3.694E-02, 4.222E-02, 4.821E-02,            FL74650
     * 5.512E-02, 6.408E-02, 7.764E-02, 9.126E-02, 1.111E-01,            FL74660
     * 1.304E-01, 1.793E-01, 2.230E-01, 3.000E-01, 4.400E-01,            FL74670
     * 5.000E-01, 6.000E-01, 7.000E-01, 1.000E+00, 1.500E+00,            FL74680
     * 2.000E+00, 2.400E+00, 2.900E+00, 3.400E+00, 4.000E+00,            FL74690
     * 4.800E+00, 6.000E+00, 7.000E+00, 8.100E+00, 8.900E+00,            FL74700
     * 8.700E+00, 7.550E+00, 5.900E+00, 4.500E+00, 3.500E+00,            FL74710
     * 2.800E+00, 1.800E+00, 1.300E+00, 8.000E-01, 4.000E-01,            FL74720
     * 1.900E-01, 2.000E-01, 5.700E-01, 7.500E-01, 7.000E-01,            FL74730
     * 4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04/            FL74740
      DATA AMOL33/                                                       FL74750
     * 2.778E-02, 2.800E-02, 2.849E-02, 3.200E-02, 3.567E-02,            FL74760
     * 4.720E-02, 5.837E-02, 7.891E-02, 1.039E-01, 1.567E-01,            FL74770
     * 2.370E-01, 3.624E-01, 5.232E-01, 7.036E-01, 8.000E-01,            FL74780
     * 9.000E-01, 1.100E+00, 1.400E+00, 1.800E+00, 2.300E+00,            FL74790
     * 2.900E+00, 3.500E+00, 3.900E+00, 4.300E+00, 4.700E+00,            FL74800
     * 5.100E+00, 5.600E+00, 6.100E+00, 6.800E+00, 7.100E+00,            FL74810
     * 7.200E+00, 6.900E+00, 5.900E+00, 4.600E+00, 3.700E+00,            FL74820
     * 2.750E+00, 1.700E+00, 1.000E-00, 5.500E-01, 3.200E-01,            FL74830
     * 2.500E-01, 2.300E-01, 5.500E-01, 8.000E-01, 8.000E-01,            FL74840
     * 4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04/            FL74850
      DATA AMOL43/                                                       FL74860
     * 2.412E-02, 2.940E-02, 3.379E-02, 3.887E-02, 4.478E-02,            FL74870
     * 5.328E-02, 6.564E-02, 7.738E-02, 9.114E-02, 1.420E-01,            FL74880
     * 1.890E-01, 3.050E-01, 4.100E-01, 5.000E-01, 6.000E-01,            FL74890
     * 7.000E-01, 8.500E-01, 1.000E+00, 1.300E+00, 1.700E+00,            FL74900
     * 2.100E+00, 2.700E+00, 3.300E+00, 3.700E+00, 4.200E+00,            FL74910
     * 4.500E+00, 5.300E+00, 5.700E+00, 6.900E+00, 7.700E+00,            FL74920
     * 7.800E+00, 7.000E+00, 5.400E+00, 4.200E+00, 3.200E+00,            FL74930
     * 2.500E+00, 1.700E+00, 1.200E+00, 8.000E-01, 4.000E-01,            FL74940
     * 2.000E-01, 1.800E-01, 6.500E-01, 9.000E-01, 8.000E-01,            FL74950
     * 4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04/            FL74960
      DATA AMOL53/                                                       FL74970
     * 1.802E-02, 2.072E-02, 2.336E-02, 2.767E-02, 3.253E-02,            FL74980
     * 3.801E-02, 4.446E-02, 7.252E-02, 1.040E-01, 2.100E-01,            FL74990
     * 3.000E-01, 3.500E-01, 4.000E-01, 6.500E-01, 9.000E-01,            FL75000
     * 1.200E+00, 1.500E+00, 1.900E+00, 2.450E+00, 3.100E+00,            FL75010
     * 3.700E+00, 4.000E+00, 4.200E+00, 4.500E+00, 4.600E+00,            FL75020
     * 4.700E+00, 4.900E+00, 5.400E+00, 5.900E+00, 6.200E+00,            FL75030
     * 6.250E+00, 5.900E+00, 5.100E+00, 4.100E+00, 3.000E+00,            FL75040
     * 2.600E+00, 1.600E+00, 9.500E-01, 6.500E-01, 5.000E-01,            FL75050
     * 3.300E-01, 1.300E-01, 7.500E-01, 8.000E-01, 8.000E-01,            FL75060
     * 4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04/            FL75070
      DATA AMOL63/                                                       FL75080
     * 2.660E-02, 2.931E-02, 3.237E-02, 3.318E-02, 3.387E-02,            FL75090
     * 3.768E-02, 4.112E-02, 5.009E-02, 5.966E-02, 9.168E-02,            FL75100
     * 1.313E-01, 2.149E-01, 3.095E-01, 3.846E-01, 5.030E-01,            FL75110
     * 6.505E-01, 8.701E-01, 1.187E+00, 1.587E+00, 2.030E+00,            FL75120
     * 2.579E+00, 3.028E+00, 3.647E+00, 4.168E+00, 4.627E+00,            FL75130
     * 5.118E+00, 5.803E+00, 6.553E+00, 7.373E+00, 7.837E+00,            FL75140
     * 7.800E+00, 7.300E+00, 6.200E+00, 5.250E+00, 4.100E+00,            FL75150
     * 3.100E+00, 1.800E+00, 1.100E+00, 7.000E-01, 3.000E-01,            FL75160
     * 2.500E-01, 3.000E-01, 5.000E-01, 7.000E-01, 7.000E-01,            FL75170
     * 4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04/            FL75180
C                                                                        FL75190
C     DATA  N2O      /                                                   FL75200
C                                                                        FL75210
      DATA AMOL14/                                                       FL75220
     * 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,            FL75230
     * 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01,            FL75240
     * 3.179E-01, 3.140E-01, 3.095E-01, 3.048E-01, 2.999E-01,            FL75250
     * 2.944E-01, 2.877E-01, 2.783E-01, 2.671E-01, 2.527E-01,            FL75260
     * 2.365E-01, 2.194E-01, 2.051E-01, 1.967E-01, 1.875E-01,            FL75270
     * 1.756E-01, 1.588E-01, 1.416E-01, 1.165E-01, 9.275E-02,            FL75280
     * 6.693E-02, 4.513E-02, 2.751E-02, 1.591E-02, 9.378E-03,            FL75290
     * 4.752E-03, 3.000E-03, 2.065E-03, 1.507E-03, 1.149E-03,            FL75300
     * 8.890E-04, 7.056E-04, 5.716E-04, 4.708E-04, 3.932E-04,            FL75310
     * 3.323E-04, 2.837E-04, 2.443E-04, 2.120E-04, 1.851E-04/            FL75320
C                                                                        FL75330
C     DATA  N2O      /                                                   FL75340
C                                                                        FL75350
      DATA AMOL24/                                                       FL75360
     * 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,            FL75370
     * 3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01, 3.163E-01,            FL75380
     * 3.096E-01, 2.989E-01, 2.936E-01, 2.860E-01, 2.800E-01,            FL75390
     * 2.724E-01, 2.611E-01, 2.421E-01, 2.174E-01, 1.843E-01,            FL75400
     * 1.607E-01, 1.323E-01, 1.146E-01, 1.035E-01, 9.622E-02,            FL75410
     * 8.958E-02, 8.006E-02, 6.698E-02, 4.958E-02, 3.695E-02,            FL75420
     * 2.519E-02, 1.736E-02, 1.158E-02, 7.665E-03, 5.321E-03,            FL75430
     * 3.215E-03, 2.030E-03, 1.397E-03, 1.020E-03, 7.772E-04,            FL75440
     * 6.257E-04, 5.166E-04, 4.352E-04, 3.727E-04, 3.237E-04,            FL75450
     * 2.844E-04, 2.524E-04, 2.260E-04, 2.039E-04, 1.851E-04/            FL75460
C                                                                        FL75470
C     DATA  N2O      /                                                   FL75480
C                                                                        FL75490
      DATA AMOL34/                                                       FL75500
     * 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,            FL75510
     * 3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01, 3.163E-01,            FL75520
     * 3.096E-01, 2.989E-01, 2.936E-01, 2.860E-01, 2.800E-01,            FL75530
     * 2.724E-01, 2.611E-01, 2.421E-01, 2.174E-01, 1.843E-01,            FL75540
     * 1.621E-01, 1.362E-01, 1.230E-01, 1.124E-01, 1.048E-01,            FL75550
     * 9.661E-02, 8.693E-02, 7.524E-02, 6.126E-02, 5.116E-02,            FL75560
     * 3.968E-02, 2.995E-02, 2.080E-02, 1.311E-02, 8.071E-03,            FL75570
     * 4.164E-03, 2.629E-03, 1.809E-03, 1.321E-03, 1.007E-03,            FL75580
     * 7.883E-04, 6.333E-04, 5.194E-04, 4.333E-04, 3.666E-04,            FL75590
     * 3.140E-04, 2.717E-04, 2.373E-04, 2.089E-04, 1.851E-04/            FL75600
C                                                                        FL75610
C     DATA  N2O      /                                                   FL75620
C                                                                        FL75630
      DATA AMOL44/                                                       FL75640
     * 3.100E-01, 3.100E-01, 3.100E-01, 3.100E-01, 3.079E-01,            FL75650
     * 3.024E-01, 2.906E-01, 2.822E-01, 2.759E-01, 2.703E-01,            FL75660
     * 2.651E-01, 2.600E-01, 2.549E-01, 2.494E-01, 2.433E-01,            FL75670
     * 2.355E-01, 2.282E-01, 2.179E-01, 2.035E-01, 1.817E-01,            FL75680
     * 1.567E-01, 1.350E-01, 1.218E-01, 1.102E-01, 9.893E-02,            FL75690
     * 8.775E-02, 7.327E-02, 5.941E-02, 4.154E-02, 3.032E-02,            FL75700
     * 1.949E-02, 1.274E-02, 9.001E-03, 6.286E-03, 4.558E-03,            FL75710
     * 2.795E-03, 1.765E-03, 1.214E-03, 8.866E-04, 6.756E-04,            FL75720
     * 5.538E-04, 4.649E-04, 3.979E-04, 3.459E-04, 3.047E-04,            FL75730
     * 2.713E-04, 2.439E-04, 2.210E-04, 2.017E-04, 1.851E-04/            FL75740
C                                                                        FL75750
C     DATA  N2O      /                                                   FL75760
C                                                                        FL75770
      DATA AMOL54/                                                       FL75780
     * 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,            FL75790
     * 3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01, 3.163E-01,            FL75800
     * 3.096E-01, 2.989E-01, 2.936E-01, 2.860E-01, 2.800E-01,            FL75810
     * 2.724E-01, 2.611E-01, 2.421E-01, 2.174E-01, 1.843E-01,            FL75820
     * 1.621E-01, 1.362E-01, 1.230E-01, 1.122E-01, 1.043E-01,            FL75830
     * 9.570E-02, 8.598E-02, 7.314E-02, 5.710E-02, 4.670E-02,            FL75840
     * 3.439E-02, 2.471E-02, 1.631E-02, 1.066E-02, 7.064E-03,            FL75850
     * 3.972E-03, 2.508E-03, 1.726E-03, 1.260E-03, 9.602E-04,            FL75860
     * 7.554E-04, 6.097E-04, 5.024E-04, 4.210E-04, 3.579E-04,            FL75870
     * 3.080E-04, 2.678E-04, 2.350E-04, 2.079E-04, 1.851E-04/            FL75880
C                                                                        FL75890
C     DATA  N2O      /                                                   FL75900
C                                                                        FL75910
      DATA AMOL64/                                                       FL75920
     * 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,            FL75930
     * 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01,            FL75940
     * 3.179E-01, 3.140E-01, 3.095E-01, 3.048E-01, 2.999E-01,            FL75950
     * 2.944E-01, 2.877E-01, 2.783E-01, 2.671E-01, 2.527E-01,            FL75960
     * 2.365E-01, 2.194E-01, 2.051E-01, 1.967E-01, 1.875E-01,            FL75970
     * 1.756E-01, 1.588E-01, 1.416E-01, 1.165E-01, 9.275E-02,            FL75980
     * 6.693E-02, 4.513E-02, 2.751E-02, 1.591E-02, 9.378E-03,            FL75990
     * 4.752E-03, 3.000E-03, 2.065E-03, 1.507E-03, 1.149E-03,            FL76000
     * 8.890E-04, 7.056E-04, 5.716E-04, 4.708E-04, 3.932E-04,            FL76010
     * 3.323E-04, 2.837E-04, 2.443E-04, 2.120E-04, 1.851E-04/            FL76020
C                                                                        FL76030
C     DATA CO        /                                                   FL76040
C                                                                        FL76050
      DATA AMOL15/                                                       FL76060
     * 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01,            FL76070
     * 1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01,            FL76080
     * 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02,            FL76090
     * 3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02,            FL76100
     * 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02,            FL76110
     * 1.521E-02, 1.722E-02, 1.995E-02, 2.266E-02, 2.487E-02,            FL76120
     * 2.738E-02, 3.098E-02, 3.510E-02, 3.987E-02, 4.482E-02,            FL76130
     * 5.092E-02, 5.985E-02, 6.960E-02, 9.188E-02, 1.938E-01,            FL76140
     * 5.688E-01, 1.549E+00, 3.849E+00, 6.590E+00, 1.044E+01,            FL76150
     * 1.705E+01, 2.471E+01, 3.358E+01, 4.148E+01, 5.000E+01/            FL76160
C                                                                        FL76170
C     DATA CO        /                                                   FL76180
C                                                                        FL76190
      DATA AMOL25/                                                       FL76200
     * 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01,            FL76210
     * 1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01,            FL76220
     * 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02,            FL76230
     * 3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02,            FL76240
     * 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02,            FL76250
     * 1.521E-02, 1.722E-02, 1.995E-02, 2.266E-02, 2.487E-02,            FL76260
     * 2.716E-02, 2.962E-02, 3.138E-02, 3.307E-02, 3.487E-02,            FL76270
     * 3.645E-02, 3.923E-02, 4.673E-02, 6.404E-02, 1.177E-01,            FL76280
     * 2.935E-01, 6.815E-01, 1.465E+00, 2.849E+00, 5.166E+00,            FL76290
     * 1.008E+01, 1.865E+01, 2.863E+01, 3.890E+01, 5.000E+01/            FL76300
C                                                                        FL76310
C     DATA CO        /                                                   FL76320
C                                                                        FL76330
      DATA AMOL35/                                                       FL76340
     * 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01,            FL76350
     * 1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01,            FL76360
     * 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02,            FL76370
     * 3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02,            FL76380
     * 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02,            FL76390
     * 1.498E-02, 1.598E-02, 1.710E-02, 1.850E-02, 1.997E-02,            FL76400
     * 2.147E-02, 2.331E-02, 2.622E-02, 3.057E-02, 3.803E-02,            FL76410
     * 6.245E-02, 1.480E-01, 2.926E-01, 5.586E-01, 1.078E+00,            FL76420
     * 1.897E+00, 2.960E+00, 4.526E+00, 6.862E+00, 1.054E+01,            FL76430
     * 1.709E+01, 2.473E+01, 3.359E+01, 4.149E+01, 5.000E+01/            FL76440
C                                                                        FL76450
C     DATA CO        /                                                   FL76460
C                                                                        FL76470
      DATA AMOL45/                                                       FL76480
     * 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01,            FL76490
     * 1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01,            FL76500
     * 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02,            FL76510
     * 3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02,            FL76520
     * 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02,            FL76530
     * 1.510E-02, 1.649E-02, 1.808E-02, 1.997E-02, 2.183E-02,            FL76540
     * 2.343E-02, 2.496E-02, 2.647E-02, 2.809E-02, 2.999E-02,            FL76550
     * 3.220E-02, 3.650E-02, 4.589E-02, 6.375E-02, 1.176E-01,            FL76560
     * 3.033E-01, 7.894E-01, 1.823E+00, 3.402E+00, 5.916E+00,            FL76570
     * 1.043E+01, 1.881E+01, 2.869E+01, 3.892E+01, 5.000E+01/            FL76580
C                                                                        FL76590
C     DATA CO        /                                                   FL76600
C                                                                        FL76610
      DATA AMOL55/                                                       FL76620
     * 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01,            FL76630
     * 1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01,            FL76640
     * 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02,            FL76650
     * 3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02,            FL76660
     * 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02,            FL76670
     * 1.521E-02, 1.722E-02, 2.037E-02, 2.486E-02, 3.168E-02,            FL76680
     * 4.429E-02, 6.472E-02, 1.041E-01, 1.507E-01, 2.163E-01,            FL76690
     * 3.141E-01, 4.842E-01, 7.147E-01, 1.067E+00, 1.516E+00,            FL76700
     * 2.166E+00, 3.060E+00, 4.564E+00, 6.877E+00, 1.055E+01,            FL76710
     * 1.710E+01, 2.473E+01, 3.359E+01, 4.149E+01, 5.000E+01/            FL76720
C                                                                        FL76730
C     DATA CO        /                                                   FL76740
C                                                                        FL76750
      DATA AMOL65/                                                       FL76760
     * 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01,            FL76770
     * 1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01,            FL76780
     * 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02,            FL76790
     * 3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02,            FL76800
     * 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02,            FL76810
     * 1.498E-02, 1.598E-02, 1.710E-02, 1.850E-02, 2.009E-02,            FL76820
     * 2.220E-02, 2.497E-02, 2.824E-02, 3.241E-02, 3.717E-02,            FL76830
     * 4.597E-02, 6.639E-02, 1.073E-01, 1.862E-01, 3.059E-01,            FL76840
     * 6.375E-01, 1.497E+00, 3.239E+00, 5.843E+00, 1.013E+01,            FL76850
     * 1.692E+01, 2.467E+01, 3.356E+01, 4.148E+01, 5.000E+01/            FL76860
C                                                                        FL76870
C     DATA  CH4      /                                                   FL76880
C                                                                        FL76890
      DATA AMOL16/                                                       FL76900
     * 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00,            FL76910
     * 1.700E+00, 1.700E+00, 1.699E+00, 1.697E+00, 1.693E+00,            FL76920
     * 1.685E+00, 1.675E+00, 1.662E+00, 1.645E+00, 1.626E+00,            FL76930
     * 1.605E+00, 1.582E+00, 1.553E+00, 1.521E+00, 1.480E+00,            FL76940
     * 1.424E+00, 1.355E+00, 1.272E+00, 1.191E+00, 1.118E+00,            FL76950
     * 1.055E+00, 9.870E-01, 9.136E-01, 8.300E-01, 7.460E-01,            FL76960
     * 6.618E-01, 5.638E-01, 4.614E-01, 3.631E-01, 2.773E-01,            FL76970
     * 2.100E-01, 1.651E-01, 1.500E-01, 1.500E-01, 1.500E-01,            FL76980
     * 1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01,            FL76990
     * 1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02/            FL77000
C                                                                        FL77010
C     DATA  CH4      /                                                   FL77020
C                                                                        FL77030
      DATA AMOL26/                                                       FL77040
     * 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.697E+00,            FL77050
     * 1.687E+00, 1.672E+00, 1.649E+00, 1.629E+00, 1.615E+00,            FL77060
     * 1.579E+00, 1.542E+00, 1.508E+00, 1.479E+00, 1.451E+00,            FL77070
     * 1.422E+00, 1.390E+00, 1.356E+00, 1.323E+00, 1.281E+00,            FL77080
     * 1.224E+00, 1.154E+00, 1.066E+00, 9.730E-01, 8.800E-01,            FL77090
     * 7.888E-01, 7.046E-01, 6.315E-01, 5.592E-01, 5.008E-01,            FL77100
     * 4.453E-01, 3.916E-01, 3.389E-01, 2.873E-01, 2.384E-01,            FL77110
     * 1.944E-01, 1.574E-01, 1.500E-01, 1.500E-01, 1.500E-01,            FL77120
     * 1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01,            FL77130
     * 1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02/            FL77140
C                                                                        FL77150
C     DATA  CH4      /                                                   FL77160
C                                                                        FL77170
      DATA AMOL36/                                                       FL77180
     * 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.697E+00,            FL77190
     * 1.687E+00, 1.672E+00, 1.649E+00, 1.629E+00, 1.615E+00,            FL77200
     * 1.579E+00, 1.542E+00, 1.508E+00, 1.479E+00, 1.451E+00,            FL77210
     * 1.422E+00, 1.390E+00, 1.356E+00, 1.323E+00, 1.281E+00,            FL77220
     * 1.224E+00, 1.154E+00, 1.066E+00, 9.730E-01, 8.800E-01,            FL77230
     * 7.931E-01, 7.130E-01, 6.438E-01, 5.746E-01, 5.050E-01,            FL77240
     * 4.481E-01, 3.931E-01, 3.395E-01, 2.876E-01, 2.386E-01,            FL77250
     * 1.944E-01, 1.574E-01, 1.500E-01, 1.500E-01, 1.500E-01,            FL77260
     * 1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01,            FL77270
     * 1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02/            FL77280
C                                                                        FL77290
C     DATA  CH4      /                                                   FL77300
C                                                                        FL77310
      DATA AMOL46/                                                       FL77320
     * 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.697E+00,            FL77330
     * 1.687E+00, 1.672E+00, 1.649E+00, 1.629E+00, 1.615E+00,            FL77340
     * 1.579E+00, 1.542E+00, 1.506E+00, 1.471E+00, 1.434E+00,            FL77350
     * 1.389E+00, 1.342E+00, 1.290E+00, 1.230E+00, 1.157E+00,            FL77360
     * 1.072E+00, 9.903E-01, 9.170E-01, 8.574E-01, 8.013E-01,            FL77370
     * 7.477E-01, 6.956E-01, 6.442E-01, 5.888E-01, 5.240E-01,            FL77380
     * 4.506E-01, 3.708E-01, 2.992E-01, 2.445E-01, 2.000E-01,            FL77390
     * 1.660E-01, 1.500E-01, 1.500E-01, 1.500E-01, 1.500E-01,            FL77400
     * 1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01,            FL77410
     * 1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02/            FL77420
C                                                                        FL77430
C     DATA  CH4      /                                                   FL77440
C                                                                        FL77450
      DATA AMOL56/                                                       FL77460
     * 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.697E+00,            FL77470
     * 1.687E+00, 1.672E+00, 1.649E+00, 1.629E+00, 1.615E+00,            FL77480
     * 1.579E+00, 1.542E+00, 1.506E+00, 1.471E+00, 1.434E+00,            FL77490
     * 1.389E+00, 1.342E+00, 1.290E+00, 1.230E+00, 1.161E+00,            FL77500
     * 1.084E+00, 1.014E+00, 9.561E-01, 9.009E-01, 8.479E-01,            FL77510
     * 7.961E-01, 7.449E-01, 6.941E-01, 6.434E-01, 5.883E-01,            FL77520
     * 5.238E-01, 4.505E-01, 3.708E-01, 3.004E-01, 2.453E-01,            FL77530
     * 1.980E-01, 1.590E-01, 1.500E-01, 1.500E-01, 1.500E-01,            FL77540
     * 1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01,            FL77550
     * 1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02/            FL77560
C                                                                        FL77570
C     DATA  CH4      /                                                   FL77580
C                                                                        FL77590
      DATA AMOL66/                                                       FL77600
     * 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00,            FL77610
     * 1.700E+00, 1.700E+00, 1.699E+00, 1.697E+00, 1.693E+00,            FL77620
     * 1.685E+00, 1.675E+00, 1.662E+00, 1.645E+00, 1.626E+00,            FL77630
     * 1.605E+00, 1.582E+00, 1.553E+00, 1.521E+00, 1.480E+00,            FL77640
     * 1.424E+00, 1.355E+00, 1.272E+00, 1.191E+00, 1.118E+00,            FL77650
     * 1.055E+00, 9.870E-01, 9.136E-01, 8.300E-01, 7.460E-01,            FL77660
     * 6.618E-01, 5.638E-01, 4.614E-01, 3.631E-01, 2.773E-01,            FL77670
     * 2.100E-01, 1.650E-01, 1.500E-01, 1.500E-01, 1.500E-01,            FL77680
     * 1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01,            FL77690
     * 1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02/            FL77700
C                                                                        FL77710
C     DATA O2        /                                                   FL77720
C                                                                        FL77730
      DATA AMOL17/                                                       FL77740
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77750
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77760
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77770
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77780
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77790
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77800
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77810
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77820
     * 2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05,            FL77830
     * 1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04/            FL77840
C                                                                        FL77850
C     DATA O2        /                                                   FL77860
C                                                                        FL77870
      DATA AMOL27/                                                       FL77880
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77890
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77900
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77910
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77920
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77930
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77940
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77950
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL77960
     * 2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05,            FL77970
     * 1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04/            FL77980
C                                                                        FL77990
C     DATA O2        /                                                   FL78000
C                                                                        FL78010
      DATA AMOL37/                                                       FL78020
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78030
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78040
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78050
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78060
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78070
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78080
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78090
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78100
     * 2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05,            FL78110
     * 1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04/            FL78120
C                                                                        FL78130
C     DATA O2        /                                                   FL78140
C                                                                        FL78150
      DATA AMOL47/                                                       FL78160
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78170
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78180
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78190
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78200
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78210
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78220
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78230
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78240
     * 2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05,            FL78250
     * 1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04/            FL78260
C                                                                        FL78270
C     DATA O2        /                                                   FL78280
C                                                                        FL78290
      DATA AMOL57/                                                       FL78300
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78310
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78320
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78330
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78340
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78350
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78360
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78370
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78380
     * 2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05,            FL78390
     * 1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04/            FL78400
C                                                                        FL78410
C     DATA O2        /                                                   FL78420
C                                                                        FL78430
      DATA AMOL67/                                                       FL78440
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78450
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78460
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78470
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78480
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78490
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78500
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78510
     * 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,            FL78520
     * 2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05,            FL78530
     * 1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04/            FL78540
C                                                                        FL78550
C     DATA DENSITY   /                                                   FL78560
C                                                                        FL78570
      DATA AMOL18/                                                       FL78580
     * 2.450E+19, 2.231E+19, 2.028E+19, 1.827E+19, 1.656E+19,            FL78590
     * 1.499E+19, 1.353E+19, 1.218E+19, 1.095E+19, 9.789E+18,            FL78600
     * 8.747E+18, 7.780E+18, 6.904E+18, 6.079E+18, 5.377E+18,            FL78610
     * 4.697E+18, 4.084E+18, 3.486E+18, 2.877E+18, 2.381E+18,            FL78620
     * 1.981E+18, 1.651E+18, 1.381E+18, 1.169E+18, 9.920E+17,            FL78630
     * 8.413E+17, 5.629E+17, 3.807E+17, 2.598E+17, 1.789E+17,            FL78640
     * 1.243E+17, 8.703E+16, 6.147E+16, 4.352E+16, 3.119E+16,            FL78650
     * 2.291E+16, 1.255E+16, 6.844E+15, 3.716E+15, 1.920E+15,            FL78660
     * 9.338E+14, 4.314E+14, 1.801E+14, 7.043E+13, 2.706E+13,            FL78670
     * 1.098E+13, 4.445E+12, 1.941E+12, 8.706E+11, 4.225E+11/            FL78680
      DATA AMOL28/                                                       FL78690
     * 2.496E+19, 2.257E+19, 2.038E+19, 1.843E+19, 1.666E+19,            FL78700
     * 1.503E+19, 1.351E+19, 1.212E+19, 1.086E+19, 9.716E+18,            FL78710
     * 8.656E+18, 7.698E+18, 6.814E+18, 6.012E+18, 5.141E+18,            FL78720
     * 4.368E+18, 3.730E+18, 3.192E+18, 2.715E+18, 2.312E+18,            FL78730
     * 1.967E+18, 1.677E+18, 1.429E+18, 1.223E+18, 1.042E+18,            FL78740
     * 8.919E+17, 6.050E+17, 4.094E+17, 2.820E+17, 1.927E+17,            FL78750
     * 1.338E+17, 9.373E+16, 6.624E+16, 4.726E+16, 3.398E+16,            FL78760
     * 2.500E+16, 1.386E+16, 7.668E+15, 4.196E+15, 2.227E+15,            FL78770
     * 1.109E+15, 4.996E+14, 1.967E+14, 7.204E+13, 2.541E+13,            FL78780
     * 9.816E+12, 3.816E+12, 1.688E+12, 8.145E+11, 4.330E+11/            FL78790
      DATA AMOL38/                                                       FL78800
     * 2.711E+19, 2.420E+19, 2.158E+19, 1.922E+19, 1.724E+19,            FL78810
     * 1.542E+19, 1.376E+19, 1.225E+19, 1.086E+19, 9.612E+18,            FL78820
     * 8.472E+18, 7.271E+18, 6.237E+18, 5.351E+18, 4.588E+18,            FL78830
     * 3.931E+18, 3.368E+18, 2.886E+18, 2.473E+18, 2.115E+18,            FL78840
     * 1.809E+18, 1.543E+18, 1.317E+18, 1.125E+18, 9.633E+17,            FL78850
     * 8.218E+17, 5.536E+17, 3.701E+17, 2.486E+17, 1.647E+17,            FL78860
     * 1.108E+17, 7.540E+16, 5.202E+16, 3.617E+16, 2.570E+16,            FL78870
     * 1.863E+16, 1.007E+16, 5.433E+15, 2.858E+15, 1.477E+15,            FL78880
     * 7.301E+14, 3.553E+14, 1.654E+14, 7.194E+13, 3.052E+13,            FL78890
     * 1.351E+13, 6.114E+12, 2.952E+12, 1.479E+12, 7.836E+11/            FL78900
      DATA AMOL48/                                                       FL78910
     * 2.549E+19, 2.305E+19, 2.080E+19, 1.873E+19, 1.682E+19,            FL78920
     * 1.508E+19, 1.357E+19, 1.216E+19, 1.088E+19, 9.701E+18,            FL78930
     * 8.616E+18, 7.402E+18, 6.363E+18, 5.471E+18, 4.699E+18,            FL78940
     * 4.055E+18, 3.476E+18, 2.987E+18, 2.568E+18, 2.208E+18,            FL78950
     * 1.899E+18, 1.632E+18, 1.403E+18, 1.207E+18, 1.033E+18,            FL78960
     * 8.834E+17, 6.034E+17, 4.131E+17, 2.839E+17, 1.938E+17,            FL78970
     * 1.344E+17, 9.402E+16, 6.670E+16, 4.821E+16, 3.516E+16,            FL78980
     * 2.581E+16, 1.421E+16, 7.946E+15, 4.445E+15, 2.376E+15,            FL78990
     * 1.198E+15, 5.311E+14, 2.022E+14, 7.221E+13, 2.484E+13,            FL79000
     * 9.441E+12, 3.624E+12, 1.610E+12, 7.951E+11, 4.311E+11/            FL79010
      DATA AMOL58/                                                       FL79020
     * 2.855E+19, 2.484E+19, 2.202E+19, 1.950E+19, 1.736E+19,            FL79030
     * 1.552E+19, 1.383E+19, 1.229E+19, 1.087E+19, 9.440E+18,            FL79040
     * 8.069E+18, 6.898E+18, 5.893E+18, 5.039E+18, 4.308E+18,            FL79050
     * 3.681E+18, 3.156E+18, 2.704E+18, 2.316E+18, 1.982E+18,            FL79060
     * 1.697E+18, 1.451E+18, 1.241E+18, 1.061E+18, 9.065E+17,            FL79070
     * 7.742E+17, 5.134E+17, 3.423E+17, 2.292E+17, 1.533E+17,            FL79080
     * 1.025E+17, 6.927E+16, 4.726E+16, 3.266E+16, 2.261E+16,            FL79090
     * 1.599E+16, 8.364E+15, 4.478E+15, 2.305E+15, 1.181E+15,            FL79100
     * 6.176E+14, 3.127E+14, 1.531E+14, 7.244E+13, 3.116E+13,            FL79110
     * 1.403E+13, 6.412E+12, 3.099E+12, 1.507E+12, 7.814E+11/            FL79120
      DATA AMOL68/                                                       FL79130
     * 2.548E+19, 2.313E+19, 2.094E+19, 1.891E+19, 1.704E+19,            FL79140
     * 1.532E+19, 1.373E+19, 1.228E+19, 1.094E+19, 9.719E+18,            FL79150
     * 8.602E+18, 7.589E+18, 6.489E+18, 5.546E+18, 4.739E+18,            FL79160
     * 4.050E+18, 3.462E+18, 2.960E+18, 2.530E+18, 2.163E+18,            FL79170
     * 1.849E+18, 1.575E+18, 1.342E+18, 1.144E+18, 9.765E+17,            FL79180
     * 8.337E+17, 5.640E+17, 3.830E+17, 2.524E+17, 1.761E+17,            FL79190
     * 1.238E+17, 8.310E+16, 5.803E+16, 4.090E+16, 2.920E+16,            FL79200
     * 2.136E+16, 1.181E+16, 6.426E+15, 3.386E+15, 1.723E+15,            FL79210
     * 8.347E+14, 3.832E+14, 1.711E+14, 7.136E+13, 2.924E+13,            FL79220
     * 1.189E+13, 5.033E+12, 2.144E+12, 9.688E+11, 5.114E+11/            FL79230
                                                                         FL79240
      DATA ANO        /                                                  FL79250
     *  3.00E-04,  3.00E-04,  3.00E-04,  3.00E-04,  3.00E-04,            FL79260
     *  3.00E-04,  3.00E-04,  3.00E-04,  3.00E-04,  3.00E-04,            FL79270
     *  3.00E-04,  3.00E-04,  3.00E-04,  2.99E-04,  2.95E-04,            FL79280
     *  2.83E-04,  2.68E-04,  2.52E-04,  2.40E-04,  2.44E-04,            FL79290
     *  2.55E-04,  2.77E-04,  3.07E-04,  3.60E-04,  4.51E-04,            FL79300
     *  6.85E-04,  1.28E-03,  2.45E-03,  4.53E-03,  7.14E-03,            FL79310
     *  9.34E-03,  1.12E-02,  1.19E-02,  1.17E-02,  1.10E-02,            FL79320
     *  1.03E-02,  1.01E-02,  1.01E-02,  1.03E-02,  1.15E-02,            FL79330
     *  1.61E-02,  2.68E-02,  7.01E-02,  2.13E-01,  7.12E-01,            FL79340
     *  2.08E+00,  4.50E+00,  7.98E+00,  1.00E+01,  1.00E+01/            FL79350
      DATA SO2       /                                                   FL79360
     *  3.00E-04,  2.74E-04,  2.36E-04,  1.90E-04,  1.46E-04,            FL79370
     *  1.18E-04,  9.71E-05,  8.30E-05,  7.21E-05,  6.56E-05,            FL79380
     *  6.08E-05,  5.79E-05,  5.60E-05,  5.59E-05,  5.64E-05,            FL79390
     *  5.75E-05,  5.75E-05,  5.37E-05,  4.78E-05,  3.97E-05,            FL79400
     *  3.19E-05,  2.67E-05,  2.28E-05,  2.07E-05,  1.90E-05,            FL79410
     *  1.75E-05,  1.54E-05,  1.34E-05,  1.21E-05,  1.16E-05,            FL79420
     *  1.21E-05,  1.36E-05,  1.65E-05,  2.10E-05,  2.77E-05,            FL79430
     *  3.56E-05,  4.59E-05,  5.15E-05,  5.11E-05,  4.32E-05,            FL79440
     *  2.83E-05,  1.33E-05,  5.56E-06,  2.24E-06,  8.96E-07,            FL79450
     *  3.58E-07,  1.43E-07,  5.73E-08,  2.29E-08,  9.17E-09/            FL79460
      DATA ANO2       /                                                  FL79470
     *  2.30E-05,  2.30E-05,  2.30E-05,  2.30E-05,  2.30E-05,            FL79480
     *  2.30E-05,  2.30E-05,  2.30E-05,  2.30E-05,  2.32E-05,            FL79490
     *  2.38E-05,  2.62E-05,  3.15E-05,  4.45E-05,  7.48E-05,            FL79500
     *  1.71E-04,  3.19E-04,  5.19E-04,  7.71E-04,  1.06E-03,            FL79510
     *  1.39E-03,  1.76E-03,  2.16E-03,  2.58E-03,  3.06E-03,            FL79520
     *  3.74E-03,  4.81E-03,  6.16E-03,  7.21E-03,  7.28E-03,            FL79530
     *  6.26E-03,  4.03E-03,  2.17E-03,  1.15E-03,  6.66E-04,            FL79540
     *  4.43E-04,  3.39E-04,  2.85E-04,  2.53E-04,  2.31E-04,            FL79550
     *  2.15E-04,  2.02E-04,  1.92E-04,  1.83E-04,  1.76E-04,            FL79560
     *  1.70E-04,  1.64E-04,  1.59E-04,  1.55E-04,  1.51E-04/            FL79570
      DATA ANH3       /                                                  FL79580
     *  5.00E-04,  5.00E-04,  4.63E-04,  3.80E-04,  2.88E-04,            FL79590
     *  2.04E-04,  1.46E-04,  9.88E-05,  6.48E-05,  3.77E-05,            FL79600
     *  2.03E-05,  1.09E-05,  6.30E-06,  3.12E-06,  1.11E-06,            FL79610
     *  4.47E-07,  2.11E-07,  1.10E-07,  6.70E-08,  3.97E-08,            FL79620
     *  2.41E-08,  1.92E-08,  1.72E-08,  1.59E-08,  1.44E-08,            FL79630
     *  1.23E-08,  9.37E-09,  6.35E-09,  3.68E-09,  1.82E-09,            FL79640
     *  9.26E-10,  2.94E-10,  8.72E-11,  2.98E-11,  1.30E-11,            FL79650
     *  7.13E-12,  4.80E-12,  3.66E-12,  3.00E-12,  2.57E-12,            FL79660
     *  2.27E-12,  2.04E-12,  1.85E-12,  1.71E-12,  1.59E-12,            FL79670
     *  1.48E-12,  1.40E-12,  1.32E-12,  1.25E-12,  1.19E-12/            FL79680
      DATA HNO3      /                                                   FL79690
     *  5.00E-05,  5.96E-05,  6.93E-05,  7.91E-05,  8.87E-05,            FL79700
     *  9.75E-05,  1.11E-04,  1.26E-04,  1.39E-04,  1.53E-04,            FL79710
     *  1.74E-04,  2.02E-04,  2.41E-04,  2.76E-04,  3.33E-04,            FL79720
     *  4.52E-04,  7.37E-04,  1.31E-03,  2.11E-03,  3.17E-03,            FL79730
     *  4.20E-03,  4.94E-03,  5.46E-03,  5.74E-03,  5.84E-03,            FL79740
     *  5.61E-03,  4.82E-03,  3.74E-03,  2.59E-03,  1.64E-03,            FL79750
     *  9.68E-04,  5.33E-04,  2.52E-04,  1.21E-04,  7.70E-05,            FL79760
     *  5.55E-05,  4.45E-05,  3.84E-05,  3.49E-05,  3.27E-05,            FL79770
     *  3.12E-05,  3.01E-05,  2.92E-05,  2.84E-05,  2.78E-05,            FL79780
     *  2.73E-05,  2.68E-05,  2.64E-05,  2.60E-05,  2.57E-05/            FL79790
      DATA OH        /                                                   FL79800
     *  4.40E-08,  4.40E-08,  4.40E-08,  4.40E-08,  4.40E-08,            FL79810
     *  4.40E-08,  4.40E-08,  4.41E-08,  4.45E-08,  4.56E-08,            FL79820
     *  4.68E-08,  4.80E-08,  4.94E-08,  5.19E-08,  5.65E-08,            FL79830
     *  6.75E-08,  8.25E-08,  1.04E-07,  1.30E-07,  1.64E-07,            FL79840
     *  2.16E-07,  3.40E-07,  5.09E-07,  7.59E-07,  1.16E-06,            FL79850
     *  2.18E-06,  5.00E-06,  1.17E-05,  3.40E-05,  8.35E-05,            FL79860
     *  1.70E-04,  2.85E-04,  4.06E-04,  5.11E-04,  5.79E-04,            FL79870
     *  6.75E-04,  9.53E-04,  1.76E-03,  3.74E-03,  7.19E-03,            FL79880
     *  1.12E-02,  1.13E-02,  6.10E-03,  1.51E-03,  2.42E-04,            FL79890
     *  4.47E-05,  1.77E-05,  1.19E-05,  1.35E-05,  2.20E-05/            FL79900
      DATA HF        /                                                   FL79910
     *  1.00E-08,  1.00E-08,  1.23E-08,  1.97E-08,  3.18E-08,            FL79920
     *  5.63E-08,  9.18E-08,  1.53E-07,  2.41E-07,  4.04E-07,            FL79930
     *  6.57E-07,  1.20E-06,  1.96E-06,  3.12E-06,  4.62E-06,            FL79940
     *  7.09E-06,  1.05E-05,  1.69E-05,  2.57E-05,  4.02E-05,            FL79950
     *  5.77E-05,  7.77E-05,  9.90E-05,  1.23E-04,  1.50E-04,            FL79960
     *  1.82E-04,  2.30E-04,  2.83E-04,  3.20E-04,  3.48E-04,            FL79970
     *  3.72E-04,  3.95E-04,  4.10E-04,  4.21E-04,  4.24E-04,            FL79980
     *  4.25E-04,  4.25E-04,  4.25E-04,  4.25E-04,  4.25E-04,            FL79990
     *  4.25E-04,  4.25E-04,  4.25E-04,  4.25E-04,  4.25E-04,            FL80000
     *  4.25E-04,  4.25E-04,  4.25E-04,  4.25E-04,  4.25E-04/            FL80010
      DATA HCL       /                                                   FL80020
     *  1.00E-03,  7.49E-04,  5.61E-04,  4.22E-04,  3.19E-04,            FL80030
     *  2.39E-04,  1.79E-04,  1.32E-04,  9.96E-05,  7.48E-05,            FL80040
     *  5.68E-05,  4.59E-05,  4.36E-05,  6.51E-05,  1.01E-04,            FL80050
     *  1.63E-04,  2.37E-04,  3.13E-04,  3.85E-04,  4.42E-04,            FL80060
     *  4.89E-04,  5.22E-04,  5.49E-04,  5.75E-04,  6.04E-04,            FL80070
     *  6.51E-04,  7.51E-04,  9.88E-04,  1.28E-03,  1.57E-03,            FL80080
     *  1.69E-03,  1.74E-03,  1.76E-03,  1.79E-03,  1.80E-03,            FL80090
     *  1.80E-03,  1.80E-03,  1.80E-03,  1.80E-03,  1.80E-03,            FL80100
     *  1.80E-03,  1.80E-03,  1.80E-03,  1.80E-03,  1.80E-03,            FL80110
     *  1.80E-03,  1.80E-03,  1.80E-03,  1.80E-03,  1.80E-03/            FL80120
      DATA HBR       /                                                   FL80130
     *  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,            FL80140
     *  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,            FL80150
     *  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,            FL80160
     *  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,            FL80170
     *  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,  1.70E-06,            FL80180
     *  1.71E-06,  1.76E-06,  1.90E-06,  2.26E-06,  2.82E-06,            FL80190
     *  3.69E-06,  4.91E-06,  6.13E-06,  6.85E-06,  7.08E-06,            FL80200
     *  7.14E-06,  7.15E-06,  7.15E-06,  7.15E-06,  7.15E-06,            FL80210
     *  7.15E-06,  7.15E-06,  7.15E-06,  7.15E-06,  7.15E-06,            FL80220
     *  7.15E-06,  7.15E-06,  7.15E-06,  7.15E-06,  7.15E-06/            FL80230
      DATA HI        /                                                   FL80240
     *  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,            FL80250
     *  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,            FL80260
     *  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,            FL80270
     *  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,            FL80280
     *  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,            FL80290
     *  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,            FL80300
     *  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,            FL80310
     *  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,            FL80320
     *  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,            FL80330
     *  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06,  3.00E-06/            FL80340
      DATA CLO       /                                                   FL80350
     *  1.00E-08,  1.00E-08,  1.00E-08,  1.00E-08,  1.00E-08,            FL80360
     *  1.00E-08,  1.00E-08,  1.00E-08,  1.01E-08,  1.05E-08,            FL80370
     *  1.21E-08,  1.87E-08,  3.18E-08,  5.61E-08,  9.99E-08,            FL80380
     *  1.78E-07,  3.16E-07,  5.65E-07,  1.04E-06,  2.04E-06,            FL80390
     *  4.64E-06,  8.15E-06,  1.07E-05,  1.52E-05,  2.24E-05,            FL80400
     *  3.97E-05,  8.48E-05,  1.85E-04,  3.57E-04,  5.08E-04,            FL80410
     *  6.07E-04,  5.95E-04,  4.33E-04,  2.51E-04,  1.56E-04,            FL80420
     *  1.04E-04,  7.69E-05,  6.30E-05,  5.52E-05,  5.04E-05,            FL80430
     *  4.72E-05,  4.49E-05,  4.30E-05,  4.16E-05,  4.03E-05,            FL80440
     *  3.93E-05,  3.83E-05,  3.75E-05,  3.68E-05,  3.61E-05/            FL80450
      DATA OCS       /                                                   FL80460
     *  6.00E-04,  5.90E-04,  5.80E-04,  5.70E-04,  5.62E-04,            FL80470
     *  5.55E-04,  5.48E-04,  5.40E-04,  5.32E-04,  5.25E-04,            FL80480
     *  5.18E-04,  5.09E-04,  4.98E-04,  4.82E-04,  4.60E-04,            FL80490
     *  4.26E-04,  3.88E-04,  3.48E-04,  3.09E-04,  2.74E-04,            FL80500
     *  2.41E-04,  2.14E-04,  1.88E-04,  1.64E-04,  1.37E-04,            FL80510
     *  1.08E-04,  6.70E-05,  2.96E-05,  1.21E-05,  4.31E-06,            FL80520
     *  1.60E-06,  6.71E-07,  4.35E-07,  3.34E-07,  2.80E-07,            FL80530
     *  2.47E-07,  2.28E-07,  2.16E-07,  2.08E-07,  2.03E-07,            FL80540
     *  1.98E-07,  1.95E-07,  1.92E-07,  1.89E-07,  1.87E-07,            FL80550
     *  1.85E-07,  1.83E-07,  1.81E-07,  1.80E-07,  1.78E-07/            FL80560
      DATA H2CO      /                                                   FL80570
     *  2.40E-03,  1.07E-03,  4.04E-04,  2.27E-04,  1.40E-04,            FL80580
     *  1.00E-04,  7.44E-05,  6.04E-05,  5.01E-05,  4.22E-05,            FL80590
     *  3.63E-05,  3.43E-05,  3.39E-05,  3.50E-05,  3.62E-05,            FL80600
     *  3.62E-05,  3.58E-05,  3.50E-05,  3.42E-05,  3.39E-05,            FL80610
     *  3.43E-05,  3.68E-05,  4.03E-05,  4.50E-05,  5.06E-05,            FL80620
     *  5.82E-05,  7.21E-05,  8.73E-05,  1.01E-04,  1.11E-04,            FL80630
     *  1.13E-04,  1.03E-04,  7.95E-05,  4.82E-05,  1.63E-05,            FL80640
     *  5.10E-06,  2.00E-06,  1.05E-06,  6.86E-07,  5.14E-07,            FL80650
     *  4.16E-07,  3.53E-07,  3.09E-07,  2.76E-07,  2.50E-07,            FL80660
     *  2.30E-07,  2.13E-07,  1.98E-07,  1.86E-07,  1.75E-07/            FL80670
      DATA HOCL      /                                                   FL80680
     *  7.70E-06,  1.06E-05,  1.22E-05,  1.14E-05,  9.80E-06,            FL80690
     *  8.01E-06,  6.42E-06,  5.42E-06,  4.70E-06,  4.41E-06,            FL80700
     *  4.34E-06,  4.65E-06,  5.01E-06,  5.22E-06,  5.60E-06,            FL80710
     *  6.86E-06,  8.77E-06,  1.20E-05,  1.63E-05,  2.26E-05,            FL80720
     *  3.07E-05,  4.29E-05,  5.76E-05,  7.65E-05,  9.92E-05,            FL80730
     *  1.31E-04,  1.84E-04,  2.45E-04,  2.96E-04,  3.21E-04,            FL80740
     *  3.04E-04,  2.48E-04,  1.64E-04,  9.74E-05,  4.92E-05,            FL80750
     *  2.53E-05,  1.50E-05,  1.05E-05,  8.34E-06,  7.11E-06,            FL80760
     *  6.33E-06,  5.78E-06,  5.37E-06,  5.05E-06,  4.78E-06,            FL80770
     *  4.56E-06,  4.37E-06,  4.21E-06,  4.06E-06,  3.93E-06/            FL80780
      DATA AN2        /                                                  FL80790
     *  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,            FL80800
     *  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,            FL80810
     *  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,            FL80820
     *  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,            FL80830
     *  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,            FL80840
     *  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,            FL80850
     *  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,            FL80860
     *  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,  7.81E+05,            FL80870
     *  7.81E+05,  7.81E+05,  7.81E+05,  7.80E+05,  7.79E+05,            FL80880
     *  7.77E+05,  7.74E+05,  7.70E+05,  7.65E+05,  7.60E+05/            FL80890
      DATA HCN       /                                                   FL80900
     *  1.70E-04,  1.65E-04,  1.63E-04,  1.61E-04,  1.60E-04,            FL80910
     *  1.60E-04,  1.60E-04,  1.60E-04,  1.60E-04,  1.60E-04,            FL80920
     *  1.60E-04,  1.60E-04,  1.60E-04,  1.59E-04,  1.57E-04,            FL80930
     *  1.55E-04,  1.52E-04,  1.49E-04,  1.45E-04,  1.41E-04,            FL80940
     *  1.37E-04,  1.34E-04,  1.30E-04,  1.25E-04,  1.19E-04,            FL80950
     *  1.13E-04,  1.05E-04,  9.73E-05,  9.04E-05,  8.46E-05,            FL80960
     *  8.02E-05,  7.63E-05,  7.30E-05,  7.00E-05,  6.70E-05,            FL80970
     *  6.43E-05,  6.21E-05,  6.02E-05,  5.88E-05,  5.75E-05,            FL80980
     *  5.62E-05,  5.50E-05,  5.37E-05,  5.25E-05,  5.12E-05,            FL80990
     *  5.00E-05,  4.87E-05,  4.75E-05,  4.62E-05,  4.50E-05/            FL81000
      DATA CH3CL     /                                                   FL81010
     *  7.00E-04,  6.70E-04,  6.43E-04,  6.22E-04,  6.07E-04,            FL81020
     *  6.02E-04,  6.00E-04,  6.00E-04,  5.98E-04,  5.94E-04,            FL81030
     *  5.88E-04,  5.79E-04,  5.66E-04,  5.48E-04,  5.28E-04,            FL81040
     *  5.03E-04,  4.77E-04,  4.49E-04,  4.21E-04,  3.95E-04,            FL81050
     *  3.69E-04,  3.43E-04,  3.17E-04,  2.86E-04,  2.48E-04,            FL81060
     *  1.91E-04,  1.10E-04,  4.72E-05,  1.79E-05,  7.35E-06,            FL81070
     *  3.03E-06,  1.32E-06,  8.69E-07,  6.68E-07,  5.60E-07,            FL81080
     *  4.94E-07,  4.56E-07,  4.32E-07,  4.17E-07,  4.05E-07,            FL81090
     *  3.96E-07,  3.89E-07,  3.83E-07,  3.78E-07,  3.73E-07,            FL81100
     *  3.69E-07,  3.66E-07,  3.62E-07,  3.59E-07,  3.56E-07/            FL81110
      DATA H2O2      /                                                   FL81120
     *  2.00E-04,  1.95E-04,  1.92E-04,  1.89E-04,  1.84E-04,            FL81130
     *  1.77E-04,  1.66E-04,  1.49E-04,  1.23E-04,  9.09E-05,            FL81140
     *  5.79E-05,  3.43E-05,  1.95E-05,  1.08E-05,  6.59E-06,            FL81150
     *  4.20E-06,  2.94E-06,  2.30E-06,  2.24E-06,  2.68E-06,            FL81160
     *  3.68E-06,  5.62E-06,  1.03E-05,  1.97E-05,  3.70E-05,            FL81170
     *  6.20E-05,  1.03E-04,  1.36E-04,  1.36E-04,  1.13E-04,            FL81180
     *  8.51E-05,  6.37E-05,  5.17E-05,  4.44E-05,  3.80E-05,            FL81190
     *  3.48E-05,  3.62E-05,  5.25E-05,  1.26E-04,  3.77E-04,            FL81200
     *  1.12E-03,  2.00E-03,  1.68E-03,  4.31E-04,  4.98E-05,            FL81210
     *  6.76E-06,  8.38E-07,  9.56E-08,  1.00E-08,  1.00E-09/            FL81220
      DATA C2H2      /                                                   FL81230
     *  3.00E-04,  1.72E-04,  9.57E-05,  6.74E-05,  5.07E-05,            FL81240
     *  3.99E-05,  3.19E-05,  2.80E-05,  2.55E-05,  2.40E-05,            FL81250
     *  2.27E-05,  2.08E-05,  1.76E-05,  1.23E-05,  7.32E-06,            FL81260
     *  4.52E-06,  2.59E-06,  1.55E-06,  8.63E-07,  5.30E-07,            FL81270
     *  3.10E-07,  1.89E-07,  1.04E-07,  5.75E-08,  2.23E-08,            FL81280
     *  8.51E-09,  4.09E-09,  2.52E-09,  1.86E-09,  1.52E-09,            FL81290
     *  1.32E-09,  1.18E-09,  1.08E-09,  9.97E-10,  9.34E-10,            FL81300
     *  8.83E-10,  8.43E-10,  8.10E-10,  7.83E-10,  7.60E-10,            FL81310
     *  7.40E-10,  7.23E-10,  7.07E-10,  6.94E-10,  6.81E-10,            FL81320
     *  6.70E-10,  6.59E-10,  6.49E-10,  6.40E-10,  6.32E-10/            FL81330
      DATA C2H6      /                                                   FL81340
     *  2.00E-03,  2.00E-03,  2.00E-03,  2.00E-03,  1.98E-03,            FL81350
     *  1.95E-03,  1.90E-03,  1.85E-03,  1.79E-03,  1.72E-03,            FL81360
     *  1.58E-03,  1.30E-03,  9.86E-04,  7.22E-04,  4.96E-04,            FL81370
     *  3.35E-04,  2.14E-04,  1.49E-04,  1.05E-04,  7.96E-05,            FL81380
     *  6.01E-05,  4.57E-05,  3.40E-05,  2.60E-05,  1.89E-05,            FL81390
     *  1.22E-05,  5.74E-06,  2.14E-06,  8.49E-07,  3.42E-07,            FL81400
     *  1.34E-07,  5.39E-08,  2.25E-08,  1.04E-08,  6.57E-09,            FL81410
     *  4.74E-09,  3.79E-09,  3.28E-09,  2.98E-09,  2.79E-09,            FL81420
     *  2.66E-09,  2.56E-09,  2.49E-09,  2.43E-09,  2.37E-09,            FL81430
     *  2.33E-09,  2.29E-09,  2.25E-09,  2.22E-09,  2.19E-09/            FL81440
      DATA PH3       /                                                   FL81450
     *  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,            FL81460
     *  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,            FL81470
     *  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,            FL81480
     *  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,            FL81490
     *  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,            FL81500
     *  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,            FL81510
     *  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,            FL81520
     *  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,            FL81530
     *  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,            FL81540
     *  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14,  1.00E-14/            FL81550
      END                                                                FL81560
C
C     ******************************************************************
C
      SUBROUTINE LDEFAL  (Z,P,T)                                         FL81570
C                                                                        FL81580
C     ****************************************************************** FL81590
C                                                                        FL81600
C     THIS SUBROUTINE LOADS ONE OF THE 6 BUILT IN ATMOSPHERIC PROFILES   FL81610
C     FROM WHICH IT WILL INTERPOLATE "DEFAULT" VALUES FOR ALTITUDE "Z"   FL81620
C                                                                        FL81630
C                                                                        FL81640
C     ***  THIS SUBROUTINE IS CALLED BY "RDUNIT" WHICH                   FL81650
C     ***  READS USER SUPPLIED INPUT PROFILES OR SINGLE VALUES           FL81660
C     ***  UNDER "MODEL = 0     " SPECIFICATIONS                         FL81670
C                                                                        FL81680
C     *** SEE DOCUMENTATION FOR CLARIFICATION ***                        FL81690
C                                                                        FL81700
C     SUBROUTINE "DEFALT"IS TRIGGERRED WHENEVER ANY ONE OF               FL81710
C     THE INPUT PARAMETERS JCHARP, JCART, (JCHAR(K),K=1,NMOL) IS = 1-6   FL81720
C                                                                        FL81730
C     FOR SIMPLICITY, ALL INTERPOLATIONS ARE DONE AT ONE TIME BECAUSE    FL81740
C     THE LAGRANGE WEIGHTS (4PT), BASED ON (ALT-Z), REMAIN UNCHANGED     FL81750
C                                                                        FL81760
C     JCHAR(K) FOR K<8 ALLOW MODEL-DEPENDENT CHOICES                     FL81770
C                                                                        FL81780
C     JCHAR=JUNIT                                                        FL81790
C                                                                        FL81800
C     1       CHOOSES TROPICAL                                           FL81810
C     2         "     MID-LATITUDE SUMMER                                FL81820
C     3         "     MID-LATITUDE WINTER                                FL81830
C     4         "     HIGH-LAT SUMMER                                    FL81840
C     5         "     HIGH-LAT WINTER                                    FL81850
C     6         "     US STANDARD                                        FL81860
C                                                                        FL81870
C                                                                        FL81880
C     JUNIT(K) FOR K>7 CHOOSES FROM THE SINGLE TRACE CONSTITUENT         FL81890
C     PROFILES, ALL APPRORIATE FOR THE US STD ATMOSPHERE                 FL81900
C                                                                        FL81910
C     ***  NOTE ***  T<0 WILL ALSO PRINT OUT A MESSAGE INDICATING        FL81920
C     ***  A POSSIBLE MISAPPLICATION OF TEMPERATURE UNITS, (K) VS (C)    FL81930
C                                                                        FL81940
C     ****************************************************************** FL81950
C                                                                        FL81960
      PARAMETER (NCASE=15, NCASE2=NCASE-2)
      parameter (mxmol=39)

      COMMON /IFIL/ IRD,IPR,IPU,NPR,NFHDRF,NPHDRF,NFHDRL,                FL81970
     *     NPHDRL,NLNGTH,KFILE,KPANEL,LINFIL,                            FL81980
     *     NFILE,IAFIL,IEXFIL,NLTEFL,LNFIL4,LNGTH4                       FL81990
      COMMON /CARD1B/ JUNITP,JUNITT,JUNIT(NCASE2),WMOL(NCASE),
     *                WAIR,JLOW                                          FL82000
C
      real*8           dum1
      CHARACTER*8      HDUM
C
      COMMON /MLATML/ ALT(50),PMATM(50,6),TMATM(50,6),AMOL(50,8,6),      FL82020
     *                dum1(6,3),HDUM(3),dum2(50,3),DUM3(50,mxmol),IDUM   
      COMMON /TRACL/ TRAC(50,22)                                         FL82040
C                                                                        FL82050
      DATA PZERO /1013.25/,TZERO/273.15/,XLOSCH/2.6868E19/               FL82060
C                                                                        FL82070
C     *** 4PT INTERPOLATION FUNCTION                                     FL82080
C                                                                        FL82090
      VAL(A1,A2,A3,A4,X1,X2,X3,X4) = A1*X1+A2*X2+A3*X3+A4*X4             FL82100
C                                                                        FL82110
C                                                                        FL82120
      NMOL = 1                                                           FL82130
      ILOWER = 0                                                         FL82140
      IUPPER = 0                                                         FL82150
      IM50 = 50                                                          FL82160
      DO 10 IM = 2, IM50                                                 FL82170
         I2 = IM                                                         FL82180
         IF (ALT(IM).GE.Z) GO TO 20                                      FL82190
   10 CONTINUE                                                           FL82200
      I2 = IM50                                                          FL82210
   20 I1 = I2-1                                                          FL82220
      I0 = I2-2                                                          FL82230
      I3 = I2+1                                                          FL82240
      IF (I0.LT.1) GO TO 30                                              FL82250
      IF (I3.GT.IM50) GO TO 40                                           FL82260
C                                                                        FL82270
      GO TO 60                                                           FL82280
C                                                                        FL82290
C     LOWER ENDPOINT CORRECTION                                          FL82300
C                                                                        FL82310
   30 CONTINUE                                                           FL82320
      ILOWER = 1                                                         FL82330
      I0 = I1                                                            FL82340
      I1 = I2                                                            FL82350
      I2 = I3                                                            FL82360
      I3 = I3+1                                                          FL82370
      GO TO 60                                                           FL82380
C                                                                        FL82390
C     UPPER ENDPOINT CORRECTION                                          FL82400
C                                                                        FL82410
   40 CONTINUE                                                           FL82420
      IUPPER = 1                                                         FL82430
      IF (Z.GT.ALT(IM50)) GO TO 50                                       FL82440
      I3 = I2                                                            FL82450
      I2 = I1                                                            FL82460
      I1 = I0                                                            FL82470
      I0 = I1-1                                                          FL82480
      GO TO 60                                                           FL82490
C                                                                        FL82500
C     UPPER ENDPOINT EXTRAPOLATION                                       FL82510
C                                                                        FL82520
   50 CONTINUE                                                           FL82530
      Z0 = ALT(I0)                                                       FL82540
      Z1 = ALT(I1)                                                       FL82550
      Z2 = ALT(I2)                                                       FL82560
      Z3 = Z2+2.*(Z-Z2)                                                  FL82570
      IUPPER = 2                                                         FL82580
      WRITE (IPR,900) Z                                                  FL82590
      STOP 'DEFAULTZ'                                                    FL82600
C                                                                        FL82610
C     I3=I2                                                              FL82620
C     GO TO 31                                                           FL82630
C                                                                        FL82640
C     LAGRANGE CONTINUATION                                              FL82650
C                                                                        FL82660
   60 CONTINUE                                                           FL82670
C                                                                        FL82680
C     LAGRANGE COEF DETERMINATION                                        FL82690
C                                                                        FL82700
      Z1 = ALT(I1)                                                       FL82710
      Z2 = ALT(I2)                                                       FL82720
      Z0 = ALT(I0)                                                       FL82730
      Z3 = ALT(I3)                                                       FL82740
      DEN1 = (Z0-Z1)*(Z0-Z2)*(Z0-Z3)                                     FL82750
      DEN2 = (Z1-Z2)*(Z1-Z3)*(Z1-Z0)                                     FL82760
      DEN3 = (Z2-Z3)*(Z2-Z0)*(Z2-Z1)                                     FL82770
      DEN4 = (Z3-Z0)*(Z3-Z1)*(Z3-Z2)                                     FL82780
      A1 = ((Z-Z1)*(Z-Z2)*(Z-Z3))/DEN1                                   FL82790
      A2 = ((Z-Z2)*(Z-Z3)*(Z-Z0))/DEN2                                   FL82800
      A3 = ((Z-Z3)*(Z-Z0)*(Z-Z1))/DEN3                                   FL82810
      A4 = ((Z-Z0)*(Z-Z1)*(Z-Z2))/DEN4                                   FL82820
C                                                                        FL82830
C                                                                        FL82840
C     TEST INPUT PARAMETERS (JUNIT'S) SEQUENTIALLY FOR TRIGGER           FL82850
C     I.E.  JUNIT(P,T,K) = 1-6                                           FL82860
C                                                                        FL82870
      IF (JUNITP.GT.6) GO TO 70                                          FL82880
      MATM = JUNITP                                                      FL82890
C                                                                        FL82900
C     WRITE (IPR,60) Z,MATM                                              FL82910
C                                                                        FL82920
      X1 =  LOG(PMATM(I0,MATM))                                          FL82930
      X2 =  LOG(PMATM(I1,MATM))                                          FL82940
      X3 =  LOG(PMATM(I2,MATM))                                          FL82950
      X4 =  LOG(PMATM(I3,MATM))                                          FL82960
      IF (IUPPER.EQ.2) X4 = X3+2*(X3-X2)                                 FL82970
      P = VAL(A1,A2,A3,A4,X1,X2,X3,X4)                                   FL82980
      P = EXP(P)                                                         FL82990
   70 IF (JUNITT.GT.6) GO TO 80                                          FL83000
      MATM = JUNITT                                                      FL83010
C                                                                        FL83020
C     WRITE (IPR,65) Z,MATM                                              FL83030
C                                                                        FL83040
      X1 = TMATM(I0,MATM)                                                FL83050
      X2 = TMATM(I1,MATM)                                                FL83060
      X3 = TMATM(I2,MATM)                                                FL83070
      X4 = TMATM(I3,MATM)                                                FL83080
      T = VAL(A1,A2,A3,A4,X1,X2,X3,X4)                                   FL83090
   80 DO 110 K = 1, NMOL                                                 FL83100
         IF (JUNIT(K).GT.6) GO TO 110                                    FL83110
C                                                                        FL83120
         IF (K.GT.7) GO TO 90                                            FL83130
         MATM = JUNIT(K)                                                 FL83140
C                                                                        FL83150
         X1 = AMOL(I0,K,MATM)                                            FL83160
         X2 = AMOL(I1,K,MATM)                                            FL83170
         X3 = AMOL(I2,K,MATM)                                            FL83180
         X4 = AMOL(I3,K,MATM)                                            FL83190
         GO TO 100                                                       FL83200
   90    ITR = K-7                                                       FL83210
         MATM = 6                                                        FL83220
C                                                                        FL83230
         X1 = TRAC(I0,ITR)                                               FL83240
         X2 = TRAC(I1,ITR)                                               FL83250
         X3 = TRAC(I2,ITR)                                               FL83260
         X4 = TRAC(I3,ITR)                                               FL83270
  100    WMOL(K) = VAL(A1,A2,A3,A4,X1,X2,X3,X4)                          FL83280
         JUNIT(K) = 10                                                   FL83290
         GO TO 110                                                       FL83300
C                                                                        FL83310
C        53 JUNIT(K)=10                                                  FL83320
C        WRITE(IPR,54)K                                                  FL83330
C        54 FORMAT('  **** INCONSISTENCY IN THE USER SPECIFICATION',     FL83340
C        A ' , JUNIT = 9 AND WMOL(K) = 0 , K =',I2,/,                    FL83350
C        B '  ****   DENNUM(K) HAS BEEN SET TO 0, NOT DEFAULT VALUE')    FL83360
C                                                                        FL83370
  110 CONTINUE                                                           FL83380
      WMOL(12) = WMOL(12)*1.0E+3                                         FL83390
C                                                                        FL83400
C     THE UNIT FOR NEW PROFILE IS PPMV.                                  FL83410
C                                                                        FL83420
      RETURN                                                             FL83430
C                                                                        FL83440
C     100  CONTINUE                                                      FL83450
C                                                                        FL83460
C                                                                        FL83470
C     STOP'DEFAULT'                                                      FL83480
C                                                                        FL83490
  900 FORMAT(/,'   *** Z IS GREATER THAN 120 KM ***, Z = ',F10.3)        FL83500
C                                                                        FL83510
      END                                                                FL83520
      BLOCK DATA ATMCOL                                                  FL83530
C                                                                        FL83540
C     >    BLOCK DATA                                                    FL83550
C     ****************************************************************** FL83560
C     THIS SUBROUTINE INITIALIZES THE CONSTANTS  USED IN THE             FL83570
C     PROGRAM. CONSTANTS RELATING TO THE ATMOSPHERIC PROFILES ARE STORED FL83580
C     IN BLOCK DATA MLATMB.                                              FL83590
C     ****************************************************************** FL83600
C                                                                        FL83610
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)
C
      COMMON /CONSTL/ PZERO,TZERO,ADCON,ALZERO,AVMWT,AMWT(MXMOL) 
      DATA PZERO/1013.25/,TZERO/273.15/                                  FL83640
C                                                                        FL83680
C     **   ALZERO IS THE MEAN LORENTZ HALFWIDTH AT PZERO AND 296.0 K.    FL83690
C     **   AVMWT IS THE MEAN MOLECULAR WEIGHT USED TO AUTOMATICALLY      FL83700
C     **   GENERATE THE LBLRTM BOUNDARIES IN AUTLAY                      FL83710
C                                                                        FL83720
      DATA ALZERO/0.1/,AVMWT/36.0/                                       FL83730
c      DATA AIRMWT / 28.964 / ,                                          FA12380
      DATA    AMWT /   18.015 ,  44.010 , 47.998 , 44.01 ,               FA12390
     *              28.011 ,  16.043 , 31.999 , 30.01 ,                  FA12400
     *              64.06  ,  46.01  , 17.03  , 63.01 ,                  FA12410
     *              17.00  ,  20.01  , 36.46  , 80.92 ,                  FA12420
     *             127.91  ,  51.45  , 60.08  , 30.03 ,                  FA12430
     *              52.46  ,  28.014 , 27.03  , 50.49 ,                  FA12440
     *              34.01  ,  26.03  , 30.07  , 34.00 ,                  FA12450
     *              66.01  , 146.05  , 34.08  , 46.03 ,                  FA12460
c     approx;
     *              33.00  ,  15.99  , 98.    , 30.00 ,
     *              97.    ,  44.5   , 32.04  /
c
      END                                                                FL83780
      SUBROUTINE YDIAR (IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH   M22460
     *                 ,RAINRT,GNDALT,YID)                                M22470
C                                                                         M22480
C     IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT       M22490
C                                                                         M22500
      DIMENSION YID(10)                                                   M22510
C                                                                         M22520
      CHARACTER*8      YID                                              & M22530
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


