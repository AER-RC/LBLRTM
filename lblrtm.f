C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      PROGRAM LBLRTM                                                      A00010
C                                                                         A00020
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   A00030
C                                                                         A00040
C                                                                         A00060
C                  IMPLEMENTATION:    R.D. WORSHAM                        A00070
C                                                                         A00080
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         A00090
C                                     R.D. WORSHAM                        A00100
C                                     J.L. MONCET                         A00110
C                                                                         A00120
C                                                                         A00130
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         A00140
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          A00150
C                                                                         A00160
C----------------------------------------------------------------------   A00170
C                                                                         A00180
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     A00190
C                                     OFFICE OF ENERGY RESEARCH           A00200
C                                     DEPARTMENT OF ENERGY                A00210
C                                                                         A00220
C                                                                         A00230
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             A00240
C                                                                         A00250
C                                             FASCOD3                     A00260
C                                                                         A00270
C                                                                         A00280
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   A00290
C                                                                         A00300
C                                                                         A00310
C**********************************************************************   A00320
C*                                                                        A00330
C*                           LBLRTM                                       A00340
C*                                                  3 JUNE 1991           A00350
C*                                                                        A00360
C*                FAST ATMOSPHERIC SIGNATURE CODE                         A00370
C*                                                                        A00380
C*                                                                        A00390
C*                AIR FORCE GEOPHYSICS LABORATORY                         A00400
C*                    HANSCOM AFB, MA 01731                               A00410
C*                                                                        A00420
C*                                                                        A00430
C*                    * SHEPARD A. CLOUGH                                 A00440
C*                      FRANCIS X. KNEIZYS                                A00450
C*                      GAIL P. ANDERSON                                  A00460
C*                      JAMES H. CHETWYND JR.                             A00470
C*                    + ROBERT D. WORSHAM                                 A00480
C*                      ERIC P. SHETTLE                                   A00490
C*                      LEONARD W. ABREU                                  A00500
C*                                                                        A00510
C*                                                                        A00520
C*                    * CURRENT ADDRESS:                                  A00530
C*                    + ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.       A00540
C*                      840 MEMORIAL DRIVE                                A00550
C*                      CAMBRIDGE, MA 02139                               A00560
C*                                                                        A00570
C**********************************************************************   A00580
C*                                                                        A00590
C*    DOCUMENTATION AND INFORMATION ABOUT THE PROGRAM  MAY BE OBTAINED    A00600
C*    FROM                                                                A00610
C*                                                                        A00620
C*    S.A. CLOUGH, AER inc., 840 MEMORIAL DRIVE, CAMBRIDGE, MA 02139      A00630
C*                                                                        A00640
C**********************************************************************   A00650
C*                                                                        A00660
C*                  LINE PARAMETER COMPILATIONS                           A00670
C*                                                                        A00680
C*   THE HITRAN DATABASE: 1986 EDITION                                    A00690
C*      L.S. ROTHMAN, R.R. GAMACHE, A. GOLDMAN, L. R. BROWN,              A00700
C*      R. A. TOTH, H. M. PICKETT, R. L. POYNTER, J.-M. FLAUD,            A00710
C*      C. CAMY-PEYRET, A. BARBE, N. HUSSON, C. P. RINSLAND,              A00720
C*      AND M. A. H. SMITH                                                A00730
C*      APPLIED OPTICS, VOL. 26,  P 4058(OCT 1987)                        A00740
C*                                                                        A00750
C*   AFGL ATMOSPHERIC ABSORPTION LINE PARAMETERS COMPILATION:             A00760
C*   1982 EDITION                                                         A00770
C*      L. S. ROTHMAN, R. R. GAMACHE, A. BARBE, A. GOLDMAN,               A00780
C*      J. R. GILLIS, L. R. BROWN, R. A. TOTH, J.-M. FLAUD AND            A00790
C*      C. CAMY-PEYRET                                                    A00800
C*      APPLIED OPTICS, VOL. 22, P 2247(AUG 1983)                         A00810
C*                                                                        A00820
C*   AFGL TRACE GAS COMPILATION: 1982 VERSION                             A00830
C*      L. S. ROTHMAN, A. GOLDMAN, J. R. GILLIS, R. R. GAMACHE,           A00840
C*      H.M. PICKETT, R. L. POYNTER, N. HUSSON AND A. CHEDIN              A00850
C*      APPLIED OPTICS VOL. 22, P 1616(JUN 1983)                          A00860
C*                                                                        A00870
C*   AFCRL ATMOSPHERIC ABSORPTION LINE PARAMETERS COMPILATION             A00880
C*                                                     AFCRL-TR-73-0096   A00890
C*      R. A. MCCLATCHEY, W. S. BENEDICT, S. A. CLOUGH, D. E. BURCH,      A00900
C*      R. F. CALFEE, K. FOX, L.S. ROTHMAN AND J. S. GARING               A00910
C*                                                                        A00920
C**********************************************************************   A00930
C*                                                                        A00940
C*       LBLATM - ATMOSPHERIC OPTICAL PROPERTIES                          A00950
C*                SPHERICAL REFRACTIVE GEOMETRY                           A00960
C*                                                                        A00970
C*    AFGL ATMOSPHERIC CONSTITUENT PROFILES (0-120 KM)  AFGL-TR-86-0110   A00980
C*       G. P. ANDERSON, S. A. CLOUGH, F.X. KNEIZYS, J. H. CHETWYND       A00990
C*       AND E. P. SHETTLE                                                A01000
C*                                                                        A01010
C*    AIR MASS COMPUTER PROGRAM FOR ATMOSPHERIC TRANSMITTANCE/RADIANCE:   A01020
C*    FSCATM                                            AFGL-TR-83-0065   A01030
C*       W. O. GALLERY, F. X. KNEIZYS, AND S. A. CLOUGH                   A01040
C*                                                                        A01050
C**********************************************************************   A01060
C*                                                                        A01070
C*      LOWTRN - AEROSOLS, HYDROMETEORS AND MOLECULAR SCATTERING          A01080
C*                                                                        A01090
C*    ATMOSPHERIC TRANSMITTANCE/RADIANCE:                                 A01100
C*    COMPUTER CODE LOWTRAN 6                           AFGL-TR-83-0187   A01110
C*       F. X. KNEIZYS, E. P. SHETTLE, W. O. GALLERY,                     A01120
C*       J. H. CHETWYND,JR, L. W. ABREU, J. E. A. SELBY,                  A01130
C*       S. A. CLOUGH AND R. W. FENN                                      A01140
C*                                                                        A01150
C*    ATMOSPHERIC TRANSMITTANCE AND RADIANCE: THE LOWTRAN 5 CODE          A01160
C*       F. X. KNEIZYS, E. P. SHETTLE, AND W. O. GALLERY                  A01170
C*       SPIE, V277 (1981) P 116                                          A01180
C*                                                                        A01190
C*    ATMOSPHERIC TRANSMITTANCE/RADIANCE:                                 A01200
C*    COMPUTER CODE LOWTRAN 5                          AFGL-TR-80-00676   A01210
C*       F. X. KNEIZYS, E. P. SHETTLE, W. O. GALLERY,                     A01220
C*       J. H. CHETWYND,JR, L. W. ABREU, J. E. A. SELBY,                  A01230
C*       R. W. FENN AND R. A. MCCLATCHEY                                  A01240
C*                                                                        A01250
C*    ATMOSPHERIC ATTENUATION OF MILLIMETER AND SUBMILLIMETER WAVES:      A01260
C*    MODELS AND COMPUTER CODES                         AFGL-TR-79-0253   A01270
C*       V. J. FALCONE, JR., L. W. ABREU AND E. P. SHETTLE                A01280
C*                                                                        A01290
C*    .    .    .    .    .    .    .    .    .    .    .    .    .       A01300
C*                                                                        A01310
C*    MODELS OF THE ATMOSPHERIC AEROSOLS AND THEIR OPTICAL PROPERTIES,    A01320
C*    IN AGARD PROC. 183, OPTICAL PROPAGATION IN THE ATMOSPHERE (1976)    A01330
C*       E. P. SHETTLE AND R. W. FENN                                     A01340
C*                                                                        A01350
C*    MODELS OF THE AEROSOLS OF THE LOWER ATMOSPHERE AND THE EFFECTS OF   A01360
C*    HUMIDITY VARIATIONS ON THEIR OPTICAL PROPERTIES   AFGL-TR-79-0214   A01370
C*       E. P. SHETTLE AND R. W. FENN                                     A01380
C*                                                                        A01390
C**********************************************************************   A01400
C*                                                                        A01410
C*           NLTE - NON LOCAL THERMODYNAMIC EQUILIBRIUM                   A01420
C*                                                                        A01430
C*                                                                        A01440
C*    ATMOSPHERIC TRANSMITTANCE/RADIANCE: COMPUTER CODE FASCOD2           A01450
C*                                                                        A01460
C*        W. L. RIDGWAY, R. A. MOOSE, AND A. C. COGLEY                    A01470
C*                                                      AFGL-TR-82-0392   A01480
C*                                                      SONICRAFT, INC.   A01490
C*                                                                        A01500
C*    .    .    .   .    .    .    .    .    .    .    .    .    .        A01510
C*                                                                        A01520
C*                        NLTE REFERENCES                                 A01530
C*                                                                        A01540
C*   A USER'S GUIDE TO THE AFGL/VISIDYNE HIGH ALTITUDE INFRARED           A01550
C*   RADIANCE MODEL COMPUTER PROGRAM                    AFGL-TR-85-0015   A01560
C*      T. C. DEGGES AND A. P. D'AGATI                  VISIDYNE/AFGL     A01570
C*                                                                        A01580
C*   A HIGH ALTITUDE INFRARED RADIANCE MODEL            AFGL-TR-77-0271   A01590
C*      T. C. DEGGES AND H. J. P. SMITH                 VISIDYNE,INC      A01600
C*                                                                        A01610
C**********************************************************************   A01620
C*                                                                        A01630
C*                                                                        A01640
C*               GENERAL LBLRTM  REFERENCES -                             A01650
C*                                                                        A01660
C*                                                                        A01670
C*    ATMOSPHERIC RADIANCE AND TRANSMITTANCE: FASCOD2                     A01680
C*        S. A. CLOUGH, F. X. KNEIZYS, E. P. SHETTLE AND G. P. ANDERSON   A01690
C*        PROC. OF THE SIXTH CONFERENCE ON ATMOSPHERIC RADIATION,         A01700
C*        WILLIAMSBURG, VA (1986), P 141                                  A01710
C*                                                                        A01720
C*    LINEAR ABSORPTION AND SCATTERING OF LASER BEAMS   AFGL-TR-84-0265   A01730
C*        F. X. KNEIZYS, S. A. CLOUGH, E. P. SHETTLE, L. S. ROTHMAN       A01740
C*        AND R. W. FENN                                                  A01750
C*                                                                        A01760
C*                                                                        A01770
C*    ATMOSPHERIC ATTENUATION OF LASER RADIATION                          A01780
C*        F. X. KNEIZYS, S. A. CLOUGH, E. P. SHETTLE                      A01790
C*        PROC. OF SPIE VOL. 410, LASER BEAM PROPAGATION(1983) P 13       A01800
C*                                                                        A01810
C*                                                                        A01820
C*    ATMOSPHERIC SPECTRAL TRANSMITTANCE AND RADIANCE-FASCOD1B            A01830
C*        S. A. CLOUGH, F. X. KNEIZYS, L. S. ROTHMAN AND W. O. GALLERY    A01840
C*        PROC. OF SPIE VOL.277 ATMOSPERIC TRANSMISSION(1981) P152        A01850
C*                                                                        A01860
C*                                                                        A01870
C*    CONVOLUTION ALGORITHM FOR THE LORENTZ FUNCTION                      A01880
C*        S. A. CLOUGH AND F. X. KNEIZYS, APPLIED OPTICS 18, 2329(1979)   A01890
C*                                                                        A01900
C*                                                                        A01910
C*    FASCODE - FAST ATMOSPHERIC SIGNATURE CODE                           A01920
C*              (SPECTRAL TRANSMITTANCE AND RADIANCE)   AFGL-TR-78-0081   A01930
C*        H. J. P. SMITH, D. J. DUBE, M. E. GARDNER,                      A01940
C*        S. A. CLOUGH, F. X. KNEIZYS, AND L. S. ROTHMAN                  A01950
C*                                                                        A01960
C*    ALGORITHM FOR THE CALCULATION OF ABSORPTION COEFFICIENT             A01970
C*          - PRESSURE BROADENED MOLECULAR TRANSITIONS  AFGL-TR-77-0164   A01980
C*        S. A. CLOUGH, F. X. KNEIZYS, J. H. CHETWYND, JR.                A01990
C*                                                                        A02000
C**********************************************************************   A02010
C----------------------------------------------------------------------   A02020
C-                                                                        A02030
C-                 FILE ASSIGNMENTS FOR LBLRTM                            A02040
C-                                                                        A02050
C-                                                                        A02060
C-    TAPE3       UNFORMATTED LINE FILE WITH LBLRTM  BLOCKING             A02070
C-                                  EXTERNAL FILE NOT REQUIRED FOR        A02080
C-                                              IHIRAC=0                  A02090
C-                                              IHIRAC=9                  A02100
C-                                              IHIRAC NE 0,9 ; ITEST=1   A02110
C-                                                                        A02120
C-    TAPE4       NLTE VIBRATIONAL TEMPERATURES (POPULATIONS) BY LAYER    A02130
C-                                        ONLY REQUIRED FOR IHIRAC=4      A02140
C-                                                                        A02150
C-    TAPE5       LBLRTM  INPUT FILE                                      A02160
C-                                                                        A02170
C-    TAPE6       LBLRTM  OUTPUT FILE                                     A02180
C-                                                                        A02190
C-    TAPE7       FILE OF MOLECULAR COLUMN AMOUNTS FROM LBLATM            A02200
C-                               ONLY FOR IATM=1; IPUNCH=1 (CARD 2.1)     A02210
C-                                                                        A02220
C-    TAPE9       FILE OF EFFECTIVE LINES FOR LBLF4 CREATED BY LINF4      A02230
C-                                                                        A02240
C-    TAPE10      OPTICAL DEPTH RESULTS FROM LINE BY LINE CALCULATION     A02250
C-                            LAST LAYER     FOR IMRG EQ 0                A02260
C-                            LAYER BY LAYER FOR IMRG EQ 1                A02270
C-                                                                        A02280
C-    TAPE11      SPECTRAL RESULTS FROM SCANFN AND INTRPL                 A02290
C-                            JEMIT=-1: ABSORPTION                        A02300
C-                            JEMIT= 0: TRANSMITTANCE                     A02310
C-                            JEMIT= 1: RADIANCE                          A02320
C-                                                                        A02330
C-    TAPE12      MONOCHROMATIC RESULTS                                   A02340
C-                            IEMIT=0: OPTICAL DEPTH                      A02350
C-                            IEMIT=1: RADIANCE/TRANSMITTANCE             A02360
C-                                 INCLUDES AEROSOL CONTRIBUTION FOR      A02370
C-                                                 IAERSL=1; IEMIT=1      A02380
C-                                                                        A02390
C-    TAPE13      MONOCHROMATIC RESULTS FOR WEIGHTING FUNCTIONS           A02400
C-                            IEMIT=0: OPTICAL DEPTH                      A02410
C-                            IEMIT=1: RADIANCE/TRANSMITTANCE             A02420
C-                                 ONLY CREATED FOR IMRG= 3 TO 18         A02430
C-                                                                        A02440
C-    TAPE14      MONCHROMATIC RESULTS INLUDING AEROSOL CONTRIBUTION      A02450
C-                            IEMIT=0: OPTICAL DEPTH                      A02460
C-                            IEMIT=1: RADIANCE/TRANSMITTANCE             A02470
C-                                 ONLY CREATED FOR IAERSL=1; IEMIT=0     A02480
C-                                                                        A02490
C-    TAPE20      SPECTRAL AEROSOL TRANSMITTANCES                         A02500
C-                            TOTAL AEROSOL  CONTRIBUTION FOR IEMIT=0     A02510
C-                            LAYER BY LAYER CONTRIBUTION FOR IEMIT=1     A02520
C-                                 ONLY CREATED FOR IAERSL=1              A02530
C-                                                                        A02540
C-    TAPE29      FILE CONTAINING VALUES OF Y FOR PLOTTING                A02550
C-                                 ONLY FOR IPLOT EQ 1                    A02560
C-                                                                        A02570
C-                                                                        A02580
C-    TAPE39      AFGL PLOT FILE                                          A02590
C-                                                                        A02600
C----------------------------------------------------------------------   A02610
C-                                                                        A02620
C-                      STATEMENT FLAGS                                   A02630
C-                                                                        A02640
C-    LBLRTM  HAS BEEN STRUCTURED TO HAVE ENHANCED PORTABILITY UNDER      A02650
C-    FORTRAN 77.  FOUR FLAGS (COLUMN73) HAVE BEEN USED TO FACILITATE     A02660
C-    PROGRAM CONVERSION.                                                 A02670
C-                                                                        A02680
C-   &    IDENTIFIES STATEMENTS REQUIRED FOR WORD SIZE LESS THAN 8        A02690
C-               CHAR. ALL STATEMENTS FLAGGED WITH & IN COLUMN 73 HAVE    A02700
C-            C& STARTING IN COLUMN 1. THESE TWO CHARACTERS MUST          A02710
C-               BE CHANGED TO BLANKS FOR COMPUTERS WITH WORD SIZE        A02720
C-               LESS THAN 8 CHARACTERS.                                  A02730
C-                                                                        A02740
C-   !    IDENTIFIES STATEMENTS REQUIRED TO DOUBLE PRECISION THE          A02750
C-               VARIABLES NEEDED FOR CALCULATIONS WHICH NEED MORE        A02760
C-               THAN 32 BITS TO OBTAIN SUFFICIENT ACCURACY (I.E.         A02770
C-               THE FREQUENCIES). STATEMENTS FLAGGED WITH ! HAVE         A02780
C-            C! STARTING IN COLUMN 1. THESE TWO CHARACTERS SHOULD BE     A02790
C-               CHANGED TO BLANKS FOR COMPUTERS HAVING SINGLE            A02800
C-               PRECISION LESS THAN 10 SIGNIFICANT DIGITS.               A02810
C-                                                                        A02820
C-   #    IDENTIFIES STATEMENTS THAT MAY BE USEFUL FOR ACCELERATED        A02830
C-               FILE DATA TRANSFER UNDER CDC AND OTHER OPERATING         A02840
C-               SYSTEMS ALLOWING BUFFERED I/0.                           A02850
C-                                                                        A02860
C-   >    IDENTIFIES STATEMENTS THAT MAY BE USEFUL FOR CONVERSION,        A02870
C-               TYPICALLY SYSTEM SPECIFIC CALLS (I.E. DATE, TIME,        A02880
C-               CPU TIME, RANDOM NUMBER, ETC.).                          A02890
C-                                                                        A02900
C----------------------------------------------------------------------   A02910
C                                                                         A02920
      IMPLICIT DOUBLE PRECISION (V)                                     ! A02930
C
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID,HDATE,HTIME          & A03050
      CHARACTER CXID*80,CFORM*11,XID8*8,IDCNTL*6                          A03430
      CHARACTER*55 CTAPE3
      CHARACTER*9 CT6FIL
      CHARACTER*8 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *                HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR
      CHARACTER*8 HVRSOL
      CHARACTER*1 CONE,CTWO,CTHREE,CFOUR,CA,CB,CC,CDOL,CPRCNT,CBLNK       A03440
      CHARACTER*1 CMRG(2),CXIDA(80)                                       A03450
C                                                                         A02940
      PARAMETER (MXFSC=200,MXLAY=MXFSC+3,MXZMD=200,MXPDIM=MXLAY+MXZMD,
     *                IM2=MXPDIM-2,MXMOL=35,MXTRAC=22)
C
C     -------------------------
C
      DIMENSION IDCNTL(14),IFSDID(17),IWD(2),IWD2(2),IWD3(2),IWD4(2)      A03280
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   A02970
     *                AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     A02980
     *                DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     A02990
     *                ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    A03000
     *                EXTID(10)                                           A03010
      COMMON /HVERSN/  HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,
     *                HVROPR,HVRPST,HVRPLT,HVRTST,HVRUTL,HVRXMR
      COMMON /ARMCM1/ HVRSOL
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       A03070
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   A03080
     *                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    A03090
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP                   A03100
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         A03110
     *                NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      A03120
     *                NLTEFL,LNFIL4,LNGTH4                                A03130
      COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),            A03160
     *                MSPNL1(MXLAY),MSLAY1,ISFILE,JSFILE,KSFILE,          A03170
     *                LSFILE,MSFILE,IEFILE,JEFILE,KEFILE                  A03180
      COMMON /LASIV/ VLAS,ILAS                                            A03190
      COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,NOZERO,NP,H1F,H2F,         A03200
     *                ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH,         A03210
     *                XVBAR, HMINF,PHIF,IERRF,HSPACE                      A03220
      COMMON /MSCONS/ AIRMAS(MXLAY),TGRND,SEMIS(3),HMINMS,HMAXMS,
     *                MSFLAG,                                             A03140
     *                MSWIT,IODFIL,MSTGLE                                 A03150
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           A03020
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           A03030
      COMMON /HDRF/ V1D,V2D,DVD,NLND,IWLD                                 A03240
      COMMON /NGTH/ VD,SD,AD,EPD,MOLD,HWHD,TMPD,PSHD,FLGD,ILS2D           A03250
      COMMON /HDRL/ V1LD,VL2D,NLD,NWDS,ILST3D                             A03260
      COMMON /RCNTRL/ ILNFLG
      COMMON /FLFORM/ CFORM                                               A03270
      COMMON /IODFLG/ DVOUT
      COMMON COMSTR(250,9)                                                A02950
      COMMON R1(3600),R2(900),R3(225)                                     A02960
      EQUIVALENCE (FSCDID(1),IFSDID(1),IHIRAC) , (FSCDID(2),ILBLF4),      A03300
     *            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),                A03310
     *            (FSCDID(5),IEMIT) , (FSCDID(6),ISCAN),                  A03320
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 A03330
     *            (FSCDID(9),JRAD) , (FSCDID(10),ITEST),                  A03340
     *            (FSCDID(11),IMRG) , (FSCDID(12),SCNID),                 A03350
     *            (FSCDID(13),HWHM) , (FSCDID(14),IDABS),                 A03360
     *            (FSCDID(15),IATM) , (FSCDID(16),LAYR1),                 A03370
     *            (FSCDID(17),NLAYFS) , (YID(1),HDATE),                   A03380
     *            (YID(2),HTIME) , (YI1,IMULT)                            A03390
      EQUIVALENCE (IWD(1),XID(1)) , (IWD2(1),V1D) , (IWD3(1),VD) ,        A03400
     *            (IWD4(1),V1LD)                                          A03410
      EQUIVALENCE (CXID,CXIDA(1))                                         A03470
C                                                                         A03480
      DATA IDCNTL / ' HIRAC',' LBLF4',' CNTNM',' AERSL',' EMISS',         A03500
     *              ' SCNFN',' FILTR','  PLOT','  TEST','  IATM',         A03510
     *              '  IMRG','  ILAS',' OPDEP',' XSECT' /                 A03520
      DATA CONE / '1'/,CTWO / '2'/,CTHREE / '3'/,CFOUR / '4'/,            A03530
     *     CA / 'A'/,CB / 'B'/,CC / 'C'/                                  A03540
      DATA CDOL / '$'/,CPRCNT / '%'/,CBLNK / ' '/,CXIDA / 80*' '/         A03550
      DATA XID8 / ' LBLRTM '/                                             A03490
C                                                                         A03560
C     DATA CFORM / 'BUFFERED   '/                                       # A03570
C     DATA CFORM / 'UNFORMATTED'/                                         A03580
C                                                                         A03590
C     Set ILNFLG to default (no line rejection files kept)
C
      ILNFLG = 0
C
C     Set name of output TAPE6, depending upon type calculation
C
      CT6FIL = 'TAPE6    '
      CTAPE3 = 'TAPE3'
C
C     -------------------------
C
C     FILE ASSIGNMENTS                                                    A03600
C                                                                         A03610
      IRD = 55                                                            A03620
      OPEN (IRD,FILE='TAPE5',STATUS='UNKNOWN')                            A03630
      IPR = 66                                                            A03640
      OPEN (IPR,FILE=CT6FIL,STATUS='UNKNOWN')                             A03650
      IPU = 7                                                             A03660
      LNFIL4 = 9                                                          A03670
      OPEN (LNFIL4,FILE='TAPE9',STATUS='UNKNOWN',FORM=CFORM)              A03680
      KFILE = 10                                                          A03690
      OPEN (KFILE,FILE='TAPE10',STATUS='UNKNOWN',FORM=CFORM)              A03700
      LFILE = 11                                                          A03710
      OPEN (LFILE,FILE='TAPE11',STATUS='UNKNOWN',FORM=CFORM)              A03720
      MFILE = 12                                                          A03730
      OPEN (MFILE,FILE='TAPE12',STATUS='UNKNOWN',FORM=CFORM)              A03740
      IODFIL = 19                                                         A03750
      IEXFIL = 20                                                         A03760
      KKSTOR = KFILE                                                      A03770
      LLSTOR = LFILE                                                      A03780
      MMSTOR = MFILE                                                      A03790
C                                                                         A03800
      IENDPL = 0                                                          A03810
      MSFLAG = 0                                                          A03820
      HMINMS = 0.0                                                        A03830
      DVSET = 0.0
      HMAXMS = 15.0                                                       A03840
      MSWIT = 0                                                           A03850
      MSTGLE = 0                                                          A03860
      ONEPL = 1.001                                                       A03870
      ONEMI = 0.999                                                       A03880
      ARGMIN = 34.                                                        A03890
      EXPMIN = EXP(-ARGMIN)                                               A03900
C                                                                         A03910
      REWIND LFILE                                                        A03920
      REWIND MFILE                                                        A03930
      NFHDRF = NWDL(IWD,LSTWDF)                                           A03940
      NPHDRF = NWDL(IWD2,IWLD)                                            A03950
      NLNGTH = NWDL(IWD3,ILS2D)                                           A03960
      NPHDRL = NWDL(IWD4,ILST3D)                                          A03970
      PI = 2.*ASIN(1.)                                                    A03980
      RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07                             A03990
      RADCN2 = PLANCK*CLIGHT/BOLTZ                                        A04000
   10 WRITE (IPR,900)                                                     A04010
C                                                                         A04020
      LOWFLG = 0                                                          A04030
      IREAD = 0                                                           A04040
      NOPR = 0                                                            A04050
C                                                                         A04060
C     XID = 80 CHARACTERS OF USER IDENTIFICATION                          A04070
C                                                                         A04080
   20 READ (IRD,905,END=80) CXID                                          A04090
      IF (CXIDA(1).EQ.CPRCNT) GO TO 90                                    A04100
      IF (CXIDA(1).NE.CDOL) GO TO 20                                      A04110
      CXIDA(1) = CBLNK                                                    A04120
      READ (CXID,910) (XID(I),I=1,10)                                     A04130
      READ (XID8,910) XID(10)                                             A04140
      CALL LBLDAT(HDATE)                                                  A04150
      CALL FTIME (HTIME)                                                  A04160
      WRITE (IPR,915) XID,HDATE,HTIME                                     A04170
      DO 30 I = 1, 17                                                     A04180
         IFSDID(I) = -99                                                  A04190
   30 CONTINUE                                                            A04200
      FSCDID(12) = -99.                                                   A04210
      FSCDID(13) = -99.                                                   A04220
      VLAS = -99.                                                         A04230
      IDABS = 0                                                           A04240
      TIME0 = -99.                                                        A04250
      CALL CPUTIM (TIME0)                                                 A04260
      WRITE (IPR,920) TIME0                                               A04270
C                                                                         A04280
      READ(IRD,925,END=80) IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,             A04290
     *                      ISCAN,IFILTR,IPLOT,ITEST,IATM,CMRG,ILAS,      A04300
     *                      IOD,IXSECT,IRAD,MPTS,NPTS                     A04310
C                                                                         A04320
      IXSCNT = IXSECT*10+ICNTNM                                           A04330
C
C     *********************** SOLAR RADIANCE ***********************
C
C     Check to be sure that no radiative transfer is being done other
C     than solar and calculate attenuated solar radiation.
C
C
      IF (IEMIT.EQ.2) THEN
         IF (IHIRAC+ILBLF4+ICNTNM+IAERSL.GT.0) THEN
            WRITE(IPR,*) 'No radiative transfer calculation allowed',
     *           ' when IEMIT=2 (solar radiance)'
            STOP 'ERROR: IHIRAC+ILBLF4+ICNTNM+IAERSL > 0'
         ENDIF
         INFLAG = 0
         IOTFLG = 0
         READ(IRD,1010) INFLAG,IOTFLG
         IF (INFLAG.EQ.0) THEN
            IFILE = MFILE
         ELSE
            IFILE = KFILE
         ENDIF
         CALL SOLINT(IFILE,LFILE,NPTS,INFLAG,IOTFLG)
         GOTO 60
      ENDIF
C                                                                         A04340
C    OPEN LINFIL DEPENDENT UPON IHIRAC AND ITEST                          A04350
C                                                                         A04360
C     Linefile name specified by CTAPE3
C
      LINFIL = 3                                                          A04370
      IF (IHIRAC.GT.0) THEN                                               A04380
         IF (ITEST.EQ.1) THEN                                             A04390
            OPEN (LINFIL,FILE=CTAPE3,STATUS='NEW',FORM=CFORM)             A04400
         ELSE                                                             A04410
            OPEN (LINFIL,FILE=CTAPE3,STATUS='OLD',FORM=CFORM)             A04420
         ENDIF                                                            A04430
      ENDIF                                                               A04440
C                                                                         A04450
C    CHECK CMRG TO SEE IF QUANTITY IS SINGLE DIGIT, DOUBLE DIGIT          A04460
C    OR CHARACTER                                                         A04470
C                                                                         A04480
      IF (CMRG(2).EQ.CA) THEN                                             A04490
         IMRG = 12                                                        A04500
      ELSEIF (CMRG(2).EQ.CB) THEN                                         A04510
         IMRG = 22                                                        A04520
      ELSEIF (CMRG(2).EQ.CC) THEN                                         A04530
         IMRG = 32                                                        A04540
      ELSE                                                                A04550
         READ (CMRG(2),930) IMRG                                          A04560
         IF (CMRG(1).EQ.CONE) IMRG = IMRG+10                              A04570
         IF (CMRG(1).EQ.CTWO) IMRG = IMRG+20                              A04580
         IF (CMRG(1).EQ.CTHREE) IMRG = IMRG+30                            A04590
         IF (CMRG(1).EQ.CFOUR) IMRG = IMRG+40
      ENDIF                                                               A04600
      IF (IPLOT.GT.0) IENDPL = 1                                          A04610
C                                                                         A04620
C     JRAD= -1  NO RADIATION TERM IN ABSORPTION COEFFICIENTS              A04630
C     JRAD=  0  RADIATION TERM PUT IN BY PANEL                            A04640
C     JRAD=  1  RADIATION TERM INCLUDED IN LINE STRENGTHS                 A04650
C                                                                         A04660
      JRAD = 1                                                            A04670
      IF (IRAD.NE.0) JRAD = -1                                            A04680
      WRITE (IPR,935) (IDCNTL(I),I=1,14)                                  A04690
      WRITE (IPR,940) IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,     A04700
     *                IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT               A04710
C                                                                         A04720
      IF (IHIRAC.EQ.4) THEN                                               A04800
         IF (IEMIT.NE.1) THEN                                             A04810
            WRITE (IPR,950)                                               A04820
            STOP ' IEMIT NE 1 FOR NLTE '                                  A04830
         ENDIF                                                            A04840
      ENDIF                                                               A04850
C
C     Check values of IOD and IMRG and
C     set IMULT equal to IOD, the flag for optical depth DV
C
      IF (IOD.GE.1.AND.IMRG.NE.1) THEN
         WRITE (IPR,955)
         STOP ' IF IOD>=1 THEN IMRG MUST = 1'
      ENDIF
      IMULT = IOD
C
      IF (IAERSL.GE.1) LOWFLG = 1                                         A04870
C                                                                         A05310
      IAFIL = 14                                                          A05400
C                                                                         A05410
C     IEXFIL=20                                                           A05420
C                                                                         A05430
      IF (IAERSL.GE.1) THEN                                               A05440
         OPEN (IAFIL,FILE='TAPE14',STATUS='UNKNOWN',FORM=CFORM)           A05450
         OPEN (IEXFIL,FILE='TAPE20',STATUS='UNKNOWN',FORM=CFORM)          A05460
         REWIND IEXFIL                                                    A05470
         LOWFLG = 1                                                       A05480
      ENDIF                                                               A05490
      NFILE = 13                                                          A05500
      MMRG = MOD(IMRG,10)                                                 A05510
      IF (MMRG.GE.3) THEN                                                 A05520
         OPEN (NFILE,FILE='TAPE13',STATUS='UNKNOWN',FORM=CFORM)           A05530
         REWIND NFILE                                                     A05540
      ENDIF                                                               A05550
      NLTEFL = 4                                                          A05560
      IF (IHIRAC.EQ.4) THEN                                               A05570
         OPEN (NLTEFL,FILE='TAPE4',STATUS='OLD')                          A05580
      ENDIF                                                               A05590
C                                                                         A05600
C     TAPE39  IS AFGL PLOT FILE                                           A05610
C                                                                         A05620
      IPLFL = 39                                                          A05630
C                                                                         A05640
      IF (ITEST.EQ.1) CALL TESTMM (LINFIL)                                A05680
C                                                                         A05690
C     IHIRAC = 1 CALL HIRAC1     VOIGT                                    A05700
C     IHIRAC = 2 CALL HIRACL     LORENTZ    NOT IMPLEMENTED               A05710
C     IHIRAC = 3 CALL HIRACD     DOPPLER    NOT IMPLEMENTED               A05720
C     IHIRAC = 4 CALL HIRACQ     NLTE VOIGT                               A05730
C     IHIRAC = 9 NO LINE BY LINE CALCULATIONS; MECHANICS PURSUED          A05740
C                                                                         A05750
C     IF IEMIT .EQ. 1  PROGRAM WILL COMPUTE EMISSION                      A05760
C                                                                         A05770
C                                                                         A05780
C     PRINT LINE FILE HEADER                                              A05790
C                                                                         A05800
      IF (IHIRAC.EQ.1.OR.IHIRAC.EQ.4.OR.ILBLF4.GE.1) CALL PRLNHD          A05810
C                                                                         A05820
C     PRINT CONTINUUM INFORMATION                                         A05830
C                                                                         A05840
      IF (ICNTNM.NE.0) CALL PRCNTM                                        A05850
C                                                                         A05860
      NOPR = 0                                                            A05870
      IF (MPTS.LT.0) NOPR = 1                                             A05880
      LAYR1 = 1                                                           A05890
      LAYRFX = LAYR1                                                      A05900
C                                                                         A05910
      IF ((IHIRAC+IAERSL+IATM+IMRG).LE.0) GO TO 60                        A05920
      IF ((IHIRAC+IAERSL+IEMIT+IATM+ILAS).GT.0) THEN                      A05930
C                                                                         A05940
         READ (IRD,970,END=80) V1,V2,SAMPLE,DVSET,ALFAL0,AVMASS,DPTMIN,   A05950
     *                         DPTFAC,ILNFLG                              A05960
C                                                                         A05970
C     OPEN LINE REJECTION FILES IF ILNFLG IS ONE OR TWO
C
      IF (ILNFLG.EQ.1) THEN
         OPEN(15,FILE='REJ1',STATUS='NEW',FORM='UNFORMATTED')
         OPEN(16,FILE='REJ4',STATUS='NEW',FORM='UNFORMATTED')
      ENDIF
      IF (ILNFLG.EQ.2) THEN
         OPEN(15,FILE='REJ1',STATUS='OLD',FORM='UNFORMATTED')
         OPEN(16,FILE='REJ4',STATUS='OLD',FORM='UNFORMATTED')
      ENDIF
C
C     IF DPTMIN < 0. SET TO DEFAULT (.0002)                               A05980
C     IF DPTFAC < 0. SET TO DEFAULT (.001)                                A05990
C                                                                         A06000
         IF (DPTMIN.LT.0.) DPTMIN = .0002                                 A06010
         IF (DPTFAC.LT.0.) DPTFAC = .001                                  A06020
         IF (V2.LE.V1.AND.ILAS.EQ.0) ILAS = 1                             A06030
      ENDIF                                                               A06040
      IF (ILAS.GT.0) THEN                                                 A06050
         V2 = V1                                                          A06060
         VLAS = V1                                                        A06070
      ENDIF                                                               A06080
      TBOUND = 0.                                                         A06090
      TMPBND = 0.                                                         A06100
      EMISIV = 0.                                                         A06110
      BNDEMI(1) = 0.                                                      A06120
      BNDEMI(2) = 0.                                                      A06130
      BNDEMI(3) = 0.                                                      A06140
      BNDRFL(1) = 0.                                                      A06150
      BNDRFL(2) = 0.                                                      A06160
      BNDRFL(3) = 0.                                                      A06170
      IBPROP = 0                                                          A06180
      IF (IEMIT.GT.0) THEN                                                A06190
         READ (IRD,970,END=80) TMPBND,(BNDEMI(IBND),IBND=1,3),            A06200
     *                         (BNDRFL(IBND),IBND=1,3)                    A06210
C                                                                         A06220
         BNDTST = ABS(BNDRFL(1))+ABS(BNDRFL(2))+ABS(BNDRFL(3))            A06230
         IF (BNDTST.NE.0.) IBPROP = 1                                     A06240
C                                                                         A06250
C     CHECK TO SEE IF EMISSIVITY AND REFLECTIVITY ARE REASONABLE          A06260
C                                                                         A06270
         XVMID = (V1+V2)/2.                                               A06280
         EMITST = BNDEMI(1)+BNDEMI(2)*XVMID+BNDEMI(3)*XVMID*XVMID         A06290
         REFTST = BNDRFL(1)+BNDRFL(2)*XVMID+BNDRFL(3)*XVMID*XVMID         A06300
         IF (EMITST.LT.0..OR.EMITST.GT.1.) THEN                           A06310
            WRITE (IPR,975) XVMID,EMITST                                  A06320
            STOP 'BNDEMI'                                                 A06330
         ENDIF                                                            A06340
         IF (REFTST.LT.0..OR.REFTST.GT.1.) THEN                           A06350
            WRITE (IPR,980) XVMID,REFTST                                  A06360
            STOP 'BNDRFL'                                                 A06370
         ENDIF                                                            A06380
C                                                                         A06390
C     TBOUND IS THE BOUNDARY TEMPERATURE. TBOUND=0. FOR NO BOUNDARY       A06400
C     EMISIV IS THE BOUNDARY EMISSIVITY                                   A06410
C     SET DEFAULT FOR EMISIV                                              A06420
C                                                                         A06430
         EMITST = ABS(BNDEMI(1))+ABS(BNDEMI(2))+ABS(BNDEMI(3))            A06440
         IF ((TMPBND.GT.0.).AND.(EMITST.EQ.0.)) BNDEMI(1) = 1.            A06450
         EMISIV = BNDEMI(1)                                               A06460
         TBOUND = TMPBND                                                  A06470
         WRITE (IPR,985) V1,V2,TBOUND,(BNDEMI(IBND),IBND=1,3),            A06480
     *                   (BNDRFL(IBND),IBND=1,3)                          A06490
C                                                                         A06500
      ENDIF                                                               A06510
      ILASRD = 0                                                          A06520
   40 CONTINUE                                                            A06530
C                                                                         A06540
      IF (ILASRD.GT.0) THEN                                               A06550
         READ (IRD,990) VLAS                                              A06560
         IF (VLAS.LT.0.) GO TO 70                                         A06570
         V1 = VLAS                                                        A06580
         V2 = VLAS                                                        A06590
      ENDIF                                                               A06600
      IF (ILAS.EQ.2) ILASRD = 1                                           A06610
C                                                                         A06620
C     If IOD = 1, then set DVOUT to DVSET as flag for interpolation
C     in PNLINT, and reset DVSET to zero to avoid ratio error message
C     in SUBROUTINE PATH.
C
      IF (IOD.EQ.1) THEN
         DVOUT = ABS(DVSET)
         DVSET = 0.0
      ENDIF
C
C     If IOD = 3, then set DVOUT to DVSET as flag for interpolation
C     in PNLINT, and preserve value of DVSET for use in
C     SUBROUTINE PATH.
C
      IF (IOD.EQ.3) THEN
         IF (DVSET.EQ.0.) STOP 'DVSET MUST BE NONZERO FOR IOD=3'
         DVOUT = ABS(DVSET)
      ENDIF
C
C     -------------------------------------------------------------
C                              CALL TREE
C     XLAYER                   ---------
C          \_ OPPATH
C                  \_ LBLATM, PATH, & LOWTRAN
C          \_ OPDPTH
C                  \_ CONTNM, LINF4, & HIRAC1
C          \_ SCNMRG
C                  \_ RDSCAN, SHRKSC, CNVRCT,
C                     CONVSC, PNLRCT,  & PANLSC
C     -------------------------------------------------------------
C
      IF (IHIRAC+IATM+IMRG.GT.0)                                          A06630
     *    CALL XLAYER (MPTS,NPTS,LFILE,MFILE,NFILE)                       A06640
C                                                                         A06650
      IF ((IAERSL.EQ.1.OR.IAERSL.EQ.7).AND.IEMIT.EQ.0) THEN               A06760
         REWIND MFILE                                                     A06770
         REWIND IAFIL                                                     A06780
         REWIND IEXFIL                                                    A06790
C                                                                         A06800
         CALL ADARSL (NPTS,IEXFIL,MFILE,IAFIL,IEMIT)                      A06810
C                                                                         A06820
      ENDIF                                                               A06830
C                                                                         A06840
      IF (IMRG.EQ.1) MFILE = KFILE                                        A06850
      IF ((MMRG.GE.3).AND.(MMRG.NE.9)) MFILE = NFILE                      A06860
      IF (ILAS.GT.0) THEN                                                 A06870
         IATM = 9                                                         A06880
         JAERSL = 0                                                       A06890
         IF (IAERSL.GE.1.AND.IEMIT.EQ.1) JAERSL = 1                       A06900
C                                                                         A06910
         CALL LASER (VLAS,MFILE,JAERSL)                                   A06920
C                                                                         A06930
         IF (IAERSL.GE.1.AND.IEMIT.EQ.0) CALL LASER (VLAS,IAFIL,1)        A06940
C                                                                         A06950
         IF (ILAS.GE.2) GO TO 40                                          A06960
      ENDIF                                                               A06970
   60 CONTINUE                                                            A06980
C                                                                         A06990
      REWIND KFILE                                                        A07000
      REWIND LFILE                                                        A07010
      REWIND MFILE                                                        A07020
C                                                                         A07030
      LFILE = LLSTOR                                                      A07040
      IF (ISCAN.EQ.1) CALL SCANFN (MFILE,LFILE)                           A07050
C                                                                         A07060
      IF (ISCAN.EQ.2) CALL INTRPL (MFILE,LFILE)                           A07070
C                                                                         A07080
      IF (ISCAN.EQ.3) CALL FFTSCN (MFILE,LFILE)                           A07081
C                                                                         A07082
      IF (IFILTR.EQ.1) CALL FLTRFN (MFILE)                                A07090
C                                                                         A07100
      IF (IPLOT.NE.0) CALL PLTLBL (IENDPL)                                A07110
C                                                                         A07120
   70 CONTINUE                                                            A07130
      CALL CPUTIM (TIME1)                                                 A07140
      TIME = TIME1-TIME0                                                  A07150
      WRITE (IPR,995) TIME1,TIME                                          A07160
      KFILE = KKSTOR                                                      A07170
      LFILE = LLSTOR                                                      A07180
      MFILE = MMSTOR                                                      A07190
      GO TO 10                                                            A07200
   80 CONTINUE                                                            A07210
      IF (IENDPL.EQ.1) CALL ENDPLT                                        A07220
      STOP ' LBLRTM EXIT; EOF ON TAPE 5 '                                 A07230
C                                                                         A07240
   90 CONTINUE                                                            A07250
      WRITE(IPR,1000) HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *                HVRPST,HVRPLT,HVRTST,HVRXMR,HVRUTL,HVRSOL
      IF (IENDPL.EQ.1) CALL ENDPLT                                        A07260
      STOP ' LBLRTM EXIT '                                                A07270
C                                                                         A07280
  900 FORMAT ('1')                                                        A07290
  905 FORMAT (A80)                                                        A07300
  910 FORMAT (10A8)                                                       A07310
  915 FORMAT ('0',10A8,2X,2(1X,A8,1X))                                    A07320
  920 FORMAT ('0  TIME ENTERING LBLRTM  ',F15.4)                          A07330
  925 FORMAT (10(4X,I1),3X,2A1,3(4X,I1),I1,I4,1X,I4)                      A07340
  930 FORMAT (I1)                                                         A07350
  935 FORMAT (14(A6,3X))                                                  A07360
  940 FORMAT (1X,I4,13I9)                                                 A07370
  945 FORMAT ('0 IOD=1 AND IEMIT=0 IS NOT IMPLEMENTED ',/,                A07380
     *        '  CHANGE IEMIT TO 1 OR IOD TO 0 ')                         A07390
  950 FORMAT ('0 IEMIT=0 IS NOT IMPLEMENTED FOR NLTE ',/,                 A07400
     *        '  CHANGE IEMIT TO 1 OR IHIRAC TO 1 ')                      A07410
  955 FORMAT ('0 IMRG MUST BE ONE WHEN NOT MERGING ',
     *        'OPTICAL DEPTHS (IOD>=1)')
  970 FORMAT (8E10.3,4X,I1)                                               A07460
  975 FORMAT ('0 FOR VNU = ',F10.3,' THE EMISSIVITY = ',E10.3,            A07470
     *        ' AND IS NOT BOUNDED BY (0.,1.) ')                          A07480
  980 FORMAT ('0 FOR VNU = ',F10.3,' THE REFLECTIVITY = ',E10.3,          A07490
     *        ' AND IS NOT BOUNDED BY (0.,1.) ')                          A07500
  985 FORMAT (5(/),'0*********** BOUNDARY PROPERTIES ***********',/,      A07510
     *        '0 V1(CM-1) = ',F12.4,/,'0 V2(CM-1) = ',F12.4,/,            A07520
     *        '0 TBOUND   = ',F12.4,5X,'BOUNDARY EMISSIVITY   = ',        A07530
     *        3(1PE11.3),/,'0',29X,'BOUNDARY REFLECTIVITY = ',            A07540
     *        3(1PE11.3))                                                 A07550
  990 FORMAT (F20.8)                                                      A07560
  995 FORMAT ('0 TIME LEAVING LBLRTM ',F15.4,' TOTAL',F15.4)              A07570
 1000 FORMAT ('0 Modules and versions used in this calculation:',/,/,5X,
     *         'lblrtm.f: ',6X,A8,10X, '  contnm.f: ',6X,A8,/,5X,
     *         'fftscn.f: ',6X,A8,10X, '  lblatm.f: ',6X,A8,/,5X,
     *         'lbllow.f: ',6X,A8,10X, ' ncargks.f: ',6X,A8,/,5X,
     *         ' oprop.f: ',6X,A8,10X, ' postsub.f: ',6X,A8,/,5X,
     *         'pltlbl.f: ',6X,A8,10X, '  testmm.f: ',6X,A8,/,5X,
     *         'xmerge.f: ',6X,A8,10X, 'util_xxx.f: ',6X,A8,/,5X,
     *         ' solar.f: ',6X,A8,10X, '            ',6X,8X,/)
 1010 FORMAT (2I5)
 1015 FORMAT (I5)
C                                                                         A07580
      END                                                                 A07590
      BLOCK DATA                                                          A07600
      COMMON /FLFORM/ CFORM                                               A03270
      COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(203),MSPNL1(203),  A07610
     *                MSLAY1,ISFILE,JSFILE,KSFILE,LSFILE,MSFILE,IEFILE,   A07620
     *                JEFILE,KEFILE                                       A07630
      COMMON /HVERSN/  HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,
     *                HVROPR,HVRPST,HVRPLT,HVRTST,HVRUTL,HVRXMR
      COMMON /ARMCM1/ HVRSOL
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           A07640
C                                                                         A07650
      CHARACTER CFORM*11                                                  A03430
      CHARACTER*8 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *                HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR
      CHARACTER*8 HVRSOL
C
      DATA CFORM / 'UNFORMATTED'/                                         A03580
      DATA PLANCK / 6.626176E-27 /,BOLTZ / 1.380662E-16 /,                A07660
     *     CLIGHT / 2.99792458E10 /,AVOG / 6.022045E23 /                  A07670
      DATA IOD / 0 /,IDIR / 0 /,ITOP / 0 /,ISURF / 0 /,MSPTS / 0 /,       A07680
     *     MSPANL /203*0 /,MSPNL1 /203*0 /,ISFILE / 0 /,JSFILE / 0 /,
     *     KSFILE / 0 /,LSFILE / 0 /,MSFILE / 0 /,IEFILE / 0 /,           A07690
     *     JEFILE / 0 /,KEFILE / 0 /,MSLAY1 / 0 /                         A07700
C
C     ASSIGN SCCS VERSION NUMBER TO MODULES
C
      DATA HVRLBL / '$Revision$' /,      HVRCNT / 'NOT USED' /,
     *     HVRFFT / 'NOT USED' /, HVRATM / 'NOT USED' /,
     *     HVRLOW / 'NOT USED' /, HVRNCG / 'NOT USED' /,
     *     HVROPR / 'NOT USED' /, HVRPST / 'NOT USED' /,
     *     HVRPLT / 'NOT USED' /, HVRTST / 'NOT USED' /,
     *     HVRUTL / 'NOT USED' /, HVRXMR / 'NOT USED' /
      DATA HVRSOL / 'NOT USED' /
C                                                                         A07710
      END                                                                 A07720
      FUNCTION NWDL (IWD,ILAST)                                           A08590
C                                                                         A08600
      DIMENSION IWD(*)                                                    A08610
C                                                                         A08620
      ILAST = -654321                                                     A08630
      DO 10 I = 1, 9000                                                   A08640
         IF (IWD(I).EQ.ILAST) THEN                                        A08650
            NWDL = I-1                                                    A08660
            RETURN                                                        A08670
         ENDIF                                                            A08680
   10 CONTINUE                                                            A08690
C                                                                         A08700
      STOP ' NWDL - IWD,ILAST '                                           A08710
C                                                                         A08720
      END                                                                 A08730
      SUBROUTINE ENDFIL (IFILE)                                           A08740
C                                                                         A08750
      DIMENSION IDUM(6)                                                   A08760
      DATA IDUM / 6*-99 /                                                 A08770
C                                                                         A08780
      CALL BUFOUT (IFILE,IDUM(1),6)                                       A08790
C                                                                         A08800
      RETURN                                                              A08810
C                                                                         A08820
      END                                                                 A08830
      SUBROUTINE SKIPFL (NUMFL,IFILE,IEOF)                                A08840
C                                                                         A08850
      DIMENSION DUM(1)                                                    A08860
C                                                                         A08870
      IF (NUMFL.LE.0) RETURN                                              A08880
      ISKIP = 0                                                           A08890
   10 CALL BUFIN (IFILE,IEOF,DUM(1),1)                                    A08900
      IF (IEOF.EQ.1) GO TO 10                                             A08910
      ISKIP = ISKIP+1                                                     A08920
      IF (ISKIP.LT.NUMFL) GO TO 10                                        A08930
C                                                                         A08940
      RETURN                                                              A08950
C                                                                         A08960
      END                                                                 A08970
      SUBROUTINE COPYFL (NPTS,KFILE,MFILE)                                A08980
C                                                                         A08990
      IMPLICIT DOUBLE PRECISION (V)                                     ! A09000
C                                                                         A09010
      COMMON TR(2410)                                                     A09020
C                                                                         A09030
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & A09040
C                                                                         A09050
      COMMON /EMIHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       A09060
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   A09070
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    A09080
      COMMON /PANL/ V1P,V2P,DVP,NLIMO                                     A09090
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         A09100
     *              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        A09110
     *              NLTEFL,LNFIL4,LNGTH4                                  A09120
C                                                                         A09130
      DIMENSION XFILHD(2),PNLHD(2)                                        A09140
      EQUIVALENCE (FSCDID(5),IEMIT) , (FSCDID(6),ISCAN)                   A09150
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHD(1),V1P)                     A09160
C                                                                         A09170
      CALL CPUTIM (TIME)                                                  A09180
      IF (NOPR.EQ.0) WRITE (IPR,900) TIME                                 A09190
      CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF)                            A09200
      CALL BUFOUT (MFILE,XFILHD(1),NFHDRF)                                A09210
C                                                                         A09220
   10 CALL BUFIN (KFILE,KEOF,PNLHD(1),NPHDRF)                             A09230
      IF (KEOF.LE.0) GO TO 20                                             A09240
      CALL BUFOUT (MFILE,PNLHD(1),NPHDRF)                                 A09250
      CALL BUFIN (KFILE,KEOF,TR(1),NLIMO)                                 A09260
      CALL BUFOUT (MFILE,TR(1),NLIMO)                                     A09270
      IF (IEMIT.EQ.0.OR.ISCAN.GT.0) GO TO 10                              A09280
      CALL BUFIN (KFILE,KEOF,TR(1),NLIMO)                                 A09290
      CALL BUFOUT (MFILE,TR(1),NLIMO)                                     A09300
      GO TO 10                                                            A09310
C                                                                         A09320
   20 CALL CPUTIM (TIME1)                                                 A09330
      TIME = TIME1-TIME                                                   A09340
      IF (NOPR.EQ.0) WRITE (IPR,905) TIME                                 A09350
C                                                                         A09360
      RETURN                                                              A09370
C                                                                         A09380
  900 FORMAT (' TIME AT THE START OF --COPYFL-- ',F10.3)                  A09390
  905 FORMAT (' TIME REQUIRED FOR --COPYFL -- ',F10.3)                    A09400
C                                                                         A09410
      END                                                                 A09420
C
C     -------------------------------------------------------------
C
      SUBROUTINE QNTIFY (CNAME,CFORM)
C
C     This subroutine counts the number of nonblank characters in
C     a string, and creates a format to add two digits to the end
C     of the nonblank string.
C
      CHARACTER*(*) CNAME
      CHARACTER*55 CTEMP
      CHARACTER*10 CFORM
      CHARACTER*1  CTEMP1(55),BLANK
C
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *                NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *                NLTEFL,LNFIL4,LNGTH4
C
      EQUIVALENCE (CTEMP,CTEMP1(1))
C
      DATA BLANK / ' ' /
C
      CTEMP = CNAME
      NCHAR=0
      DO 10 J = 1, 55
         IF (CTEMP1(J).EQ.BLANK) THEN
            NCHAR = J-1
            IF (NCHAR.EQ.1) THEN
               WRITE(IPR,900) CNAME
               STOP 'PATH NOT LEFT JUSTIFIED'
            ENDIF
            GOTO 20
         ENDIF
 10   CONTINUE
C
 20   CONTINUE
      WRITE(CFORM,910) NCHAR
C
      RETURN
C
 900  FORMAT (' QNTIFY: PATH NOT LEFT JUSTIFIED: ',A55)
 910  FORMAT ('(A',I2.2,',I2.2)')
C
      END
C
C     ----------------------------------------------------------------
C
      SUBROUTINE PRLNHD                                                   A09430
C                                                                         A09440
      IMPLICIT DOUBLE PRECISION (V)                                     ! A09450
C                                                                         A09460
C     PRLNHD PRINTS OUT LINE FILE HEADER                                  A09470
C                                                                         A09480
      PARAMETER (MXMOL=35)
C
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         A09490
     *              NLNGTH,KFILE,KPANEL,LINFIL,NDFLE,IAFIL,IEXFIL,        A09500
     *              NLTEFL,LNFIL4,LNGTH4                                  A09510
C                                                                         A09520
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & A09530
C                                                                         A09540
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       A09550
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   A09560
     *                EMISIV,FSCDID(17),NMOL,LAYHDR,YI1,YID(10),LSTWDF    A09570
C                                                                         A09580
      DOUBLE PRECISION HLINID,BMOLID,HID1,HLINHD                        & A09590
C                                                                         A09600
      COMMON /LINHDR/ HLINID(10),BMOLID(64),MOLCNT(64),MCNTLC(64),        A09610
     *                MCNTNL(64),SUMSTR(64),LINMOL,FLINLO,FLINHI,         A09620
     *                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1(2),LSTWDL     A09630
C                                                                         A09640
C     LSTWD (LAST WORD) IS DUMMY, DOES NOT NEED TO BE COUNTED             A09650
C                                                                         A09660
      DIMENSION HLINHD(2),IWD(2)                                          A09670
C                                                                         A09680
      CHARACTER CHID10*8,CHARID*5,CHARDT*2,CHARI*1,CHTST*1                A09690
C                                                                         A09700
      EQUIVALENCE (HLINID(1),HLINHD(1),IWD(1))                            A09710
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),                A09720
     *            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),                A09730
     *            (FSCDID(5),IEMIT) , (FSCDID(7),IPLOT),                  A09740
     *            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                  A09750
     *            (FSCDID(10),ITEST) , (FSCDID(11),IMRG),                 A09760
     *            (FSCDID(16),LAYR1) , (FSCDID(17),NLAYHD)                A09770
C                                                                         A09780
      DATA CHARI / 'I'/                                                   A09790
C                                                                         A09800
      NFHDRL = NWDL(IWD,LSTWDL)                                           A09810
C                                                                         A09820
      REWIND LINFIL                                                       A09830
      CALL BUFIN (LINFIL,LEOF,HLINHD(1),NFHDRL)                           A09840
      IF (LEOF.LE.0) STOP 'LAYER; TAPE3 DOES NOT EXIST'                   A09850
C                                                                         A09860
      DO 10 M = 1, LINMOL                                                 A09870
         HMOLID(M) = BMOLID(M)                                            A09880
   10 CONTINUE                                                            A09890
      WRITE (IPR,900)                                                     A09900
      WRITE (IPR,905) HLINID,HID1                                         A09910
      WRITE (IPR,910)                                                     A09920
      WRITE (IPR,915) (BMOLID(I),MOLCNT(I),MCNTLC(I),MCNTNL(I),           A09930
     *                SUMSTR(I),I=1,LINMOL)                               A09940
C                                                                         A09950
      WRITE (IPR,920) FLINLO,FLINHI,LINCNT                                A09960
C                                                                         A09970
C     CHECK HEADER FOR FLAG INDICATING COMPATIBILITY WITH ISOTOPES        A09980
C                                                                         A09990
 30   WRITE (CHID10,925) HLINID(10)                                       A10000
      READ (CHID10,930) CHARID,CHARDT,CHTST                               A10010
      IF (CHTST.NE.CHARI) THEN                                            A10020
         WRITE (IPR,935) CHARID,CHARDT,CHTST                              A10030
         STOP ' PRLNHD - NO ISOTOPE INFO ON LINFIL '                      A10040
      ENDIF                                                               A10050
C                                                                         A10060
      RETURN                                                              A10070
C                                                                         A10080
  900 FORMAT ('0'/'0',20X,'   LINE FILE INFORMATION ')                    A10090
  905 FORMAT ('0',10A8,2X,2(1X,A8,1X))                                    A10100
  910 FORMAT ('0',/,23X,'COUPLED',4X,'NLTE',3X,'SUM LBLRTM ',/,7X,        A10110
     *        'MOL',5X,'LINES',4X,'LINES',4X,'LINES',4X,'STRENGTHS',/)    A10120
  915 FORMAT (' ',4X,A6,' = ',I6,3X,I6,3X,I6,2X,1PE12.4,0P)               A10130
  920 FORMAT (/,'0 LOWEST LINE = ',F10.3,5X,'  HIGHEST LINE = ',F10.3,    A10140
     *        5X,' TOTAL NUMBER OF LINES =',I8)                           A10150
  925 FORMAT (A8)                                                         A10160
  930 FORMAT (A5,A2,A1)                                                   A10170
  935 FORMAT (3(/),10X,'LINEFILE PROGRAM: ',A5,3X,'VERSION: ',A2,A1,      A10180
     *        3(/),3X,52('*'),/,3X,'* THE LINEFILE (TAPE3) IS NOT ',      A10190
     *        'COMPATIBLE WITH THIS *',/,3X,'* VERSION OF LBLRTM .',      A10200
     *        '  ISOTOPIC INFORMATION (FROM  *',/,3X,'* HITRAN) ',        A10210
     *        'MUST BE PRESERVED ON TAPE3.  USE A TAPE3 *',/,3X,          A10220
     *        '* CREATED WITH THE 91I OR LATER VERSION OF LNFL.   *',     A10230
     *        /,3X,52('*'))                                               A10240
C                                                                         A10250
      END                                                                 A10260
      SUBROUTINE EXPINT (X,X1,X2,A)                                       A10600
C                                                                         A10610
C**********************************************************************   A10620
C     THIS SUBROUTINE EXPONENTIALLY INTERPOLATES X1 AND X2 TO X BY        A10630
C     THE FACTOR A                                                        A10640
C**********************************************************************   A10650
C                                                                         A10660
      IF (X1.EQ.0.0.OR.X2.EQ.0.0) GO TO 10                                A10670
      X = X1*(X2/X1)**A                                                   A10680
C                                                                         A10690
      RETURN                                                              A10700
C                                                                         A10710
   10 X = X1+(X2-X1)*A                                                    A10720
C                                                                         A10730
      RETURN                                                              A10740
C                                                                         A10750
      END                                                                 A10760
      SUBROUTINE XLAYER (MPTS,NPTS,LFILE,MFILE,NFILE)                     A11210
C                                                                         A11220
      IMPLICIT DOUBLE PRECISION (V)                                     ! A11230
C                                                                         A11240
C**********************************************************************
C     XLAYER CONTROLS LAYER BY LAYER CALCULATION                          A11250
C**********************************************************************
C                                                                         A11260
      PARAMETER (MXFSC=200,MXLAY=MXFSC+3,MXZMD=200,MXPDIM=MXLAY+MXZMD,
     *                IM2=MXPDIM-2,MXMOL=35,MXTRAC=22)
C
      CHARACTER*55 PATH1
      CHARACTER*10 HFORM1
C
      COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,NOZERO,NP,H1F,H2F,         A11270
     *                ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH,         A11280
     *                XVBAR, HMINF,PHIF,IERRF,HSPACE                      A11290
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   A11300
     *              AVFIX,LAYER,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,        A11310
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       A11320
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      A11330
     *              EXTID(10)                                             A11340
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           A11350
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         A11360
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,        A11370
     *              NLTEFL,LNFIL4,LNGTH4                                  A11380
      COMMON /MSCONS/ AIRMAS(MXLAY),TGRND,SEMIS(3),HMINMS,HMAXMS,
     *                MSFLAG,                                             A11390
     *                MSWIT,IODFIL,MSTGLE                                 A11400
      COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),            A11410
     *                MSPNL1(MXLAY),MSLAY1,ISFILE,JSFILE,KSFILE,          A11420
     *                LSFILE,MSFILE,IEFILE,JEFILE,KEFILE                  A11430
C                                                                         A11440
C     COMMON /MLTSCT/ TAUGAS(2410),FUPC(2410),RUPC(2410)                  A11450
C                                                                         A11460
      COMMON /LASIV/ VLAS,ILAS                                            A11470
C                                                                         A11480
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & A11490
C                                                                         A11500
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       A11510
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   A11520
     *                EMISIV,FSCDID(17),NMOL,LAYHDR,YI1,YID(10),LSTWDF    A11530
C                                                                         A11540
      DOUBLE PRECISION XI1,SECAN1,HMOLI1,XALT1,Y1D                      & A11550
C                                                                         A11560
      COMMON /FILHD1/ XI1(10),SECAN1,PAV1,TAV1,HMOLI1(60),XALT1(4),       A11570
     *                W1(60),PDL,PDU,TDL,TDU,W12   ,D1 ,VD1,VD2,TBOUN1,   A11580
     *                EMISI1,FSCDI1(17),NMO1,LAYHD1,YD1,Y1D(10),LSTWDD    A11590
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP                   A11600
C                                                                         A11610
      DIMENSION FILDUM(2),FILDU1(2)                                       A11620
      DIMENSION NTAN(160)                                                 A11630
C                                                                         A11640
C     EQUIVALENCE STATEMENTS FOR TANGENT WEIGHTING FNS.                   A11650
C                                                                         A11660
      EQUIVALENCE (Y1D(10),LTNSAV) , (Y1D(9),LH1SAV) , (Y1D(8),LH2SAV)    A11670
C                                                                         A11680
      EQUIVALENCE (XID(1),FILDUM(1)) , (XI1(1),FILDU1(1))                 A11690
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),                A11700
     *            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),                A11710
     *            (FSCDID(5),IEMIT) , (FSCDID(6),ISCNHD),                 A11720
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),                 A11730
     *            (FSCDID(9),JRAD) , (FSCDID(10),ITEST),                  A11740
     *            (FSCDID(11),IMRG) , (FSCDID(12),SCNID),                 A11750
     *            (FSCDID(13),HWHM) , (FSCDID(14),IDABS),                 A11760
     *            (FSCDID(15),IATM) , (FSCDID(16),LAYR1),                 A11770
     *            (FSCDID(17),NLAYHD) , (FSCDI1(8),IPTHD1),               A11780
     *            (FSCDI1(17),NLAYD1)                                     A11790
C                                                                         A11800
      DV = 0.                                                             A11810
C                                                                         A11820
C     IF  IANT.EQ. 1  THEN POSTERIOR MERGE                                A11830
C     IF  IANT.EQ. 0  THEN NORMAL MERGE                                   A11840
C     IF  IANT.EQ.-1  THEN ANTERIOR MERGE                                 A11850
C                                                                         A11860
      IANT = 0                                                            A11870
      JPATHL = 0                                                          A11880
C                                                                         A11890
      IBUF = 1                                                            A11900
      JEMIT = 0                                                           A11910
C                                                                         A11920
      LAYER = 0                                                           A11930
      P0 = 1013.                                                          A11940
      TEMP0 = 296.                                                        A11950
C                                                                         A11960
C**********************************************************************   A11970
C                                                                         A11980
C     IMRG = 0 AND IMRG > 10 ONLY LAST LAYER OPTICAL DEPTH ON KFILE       A11990
C       (EXCEPT 12, 22 AND 32 WHICH ARE SEQUENTIAL FROM PRESTORE)         A12000
C                                                                         A12010
C     1 < IMRG < 10 OPTICAL DEPTHS ARE SEQUENTIAL ON KFILE BY LAYER       A12020
C
C     IMRG = 1 Optical Depths are stored on different files by layer
C
C     IMRG = 40/41 RADIANCE and TRANSMITTANCE are calculated from
C               multiple optical depth files by layer
C                                                                         A12030
C**********************************************************************   A12040
C                                                                         A12050
C      IMRG= 0 STANDARD MERGE                                             A12060
C      IMRG= 1 NO MERGE TAKES PLACE                                       A12070
C      IMRG= 2 SEQUENTIAL OPTICAL DEPTHS FROM KFILE MERGED ONTO MFILE     A12080
C      IMRG= A SEQUENTIAL OPTICAL DEPTHS FROM KFILE USED TO CALCULATE     A12090
C                  MONOCHROMATIC RESULTS WHICH ARE MERGED ONTO MFILE      A12100
C                  FORCES IPATHL = 1 (SAME AS IMRG=12)                    A12110
C      IMRG= B SEQUENTIAL OPTICAL DEPTHS FROM KFILE USED TO CALCULATE     A12120
C                  MONOCHROMATIC RESULTS WHICH ARE MERGED ONTO MFILE      A12130
C                  FORCES IPATHL = 2 (SAME AS IMRG=22)                    A12140
C      IMRG= C SEQUENTIAL OPTICAL DEPTHS FROM KFILE USED TO CALCULATE     A12150
C                  MONOCHROMATIC RESULTS WHICH ARE MERGED ONTO MFILE      A12160
C                  FORCES IPATHL = 3 (SAME AS IMRG=32)                    A12170
C                                                                         A12180
C      IMRG= 9 MULTIPLE RUNS EMISSION FOR AEROSOLS FROM PRESTORE          A12190
C                                                                         A12200
C        ****  WEIGHTING FUNCTIONS - MONOCHROMATIC  ****                  A12210
C                                                                         A12220
C      IMRG= 3 --- SPACE TO GROUND                                        A12230
C      IMRG= 4 --- GROUND TO SPACE                                        A12240
C      IMRG= 5 --- SPACE TO GROUND FROM PRESTORED OPTICAL DEPTHS          A12250
C      IMRG= 6 --- GROUND TO SPACE FROM PRESTORED OPTICAL DEPTHS          A12260
C      IMRG= 7 --- H1 THROUGH H(TAN) TO H2                                A12270
C      IMRG= 8 --- H1 THROUGH H(TAN) TO H2 FROM PRESTORED OPTICAL         A12280
C                  DEPTHS                                                 A12290
C                                                                         A12300
C        ****  WEIGHTING FUNCTIONS - SCANNED ****                         A12310
C                                                                         A12320
C      IMRG=13 --- SPACE TO GROUND                                        A12330
C      IMRG=14 --- GROUND TO SPACE                                        A12340
C      IMRG=15 --- SPACE TO GROUND FROM PRESTORED OPTICAL DEPTHS          A12350
C      IMRG=16 --- GROUND TO SPACE FROM PRESTORED OPTICAL DEPTHS          A12360
C      IMRG=17 --- H1 THROUGH H(TAN) TO H2                                A12370
C      IMRG=18 --- H1 THROUGH H(TAN) TO H2 FROM PRESTORED OPTICAL         A12380
C                  DEPTHS                                                 A12390
C                                                                         A12400
C        ****  WEIGHTING FUNCTIONS - FILTERED ****                        A12410
C                                                                         A12420
C      IMRG=23 --- SPACE TO GROUND                                        A12430
C      IMRG=24 --- GROUND TO SPACE                                        A12440
C      IMRG=25 --- SPACE TO GROUND FROM PRESTORED OPTICAL DEPTHS          A12450
C      IMRG=26 --- GROUND TO SPACE FROM PRESTORED OPTICAL DEPTHS          A12460
C      IMRG=27 --- H1 THROUGH H(TAN) TO H2                                A12470
C      IMRG=28 --- H1 THROUGH H(TAN) TO H2 FROM PRESTORED OPTICAL         A12480
C                  DEPTHS                                                 A12490
C                                                                         A12500
C        ****  FLUX CALCULATIONS -SCANNED ****                            A12510
C                                                                         A12520
C      IMRG=35 --- SPACE TO GROUND FROM PRESTORED OPTICAL DEPTHS          A12530
C      IMRG=36 --- GROUND TO SPACE FROM PRESTORED OPTICAL DEPTHS          A12540
C                                                                         A12550
C        ****  RADIANCE CALCULATIONS  **** 
C                                                                         A12520
C      IMRG=40 --- Downward radiance from prestored optical depths,
C                  monochromatic
C      IMRG=41 --- Upward radiance from prestored optical depths,
C                  monochromatic
C                                                                         A12550
C      WEIGHTING FUNCTION RESULTS ARE ON NFILE SEPARATED BY               A12560
C      INTERNAL 'EOF'                                                     A12570
C                                                                         A12580
C**********************************************************************   A12590
C                                                                         A12600
C     ---------------------
C     Special Merge Options
C     ---------------------
C
C     If IMRG = 1, then calculate optical depths in standard fashion,
C     but output to different files for each layer.
C
C                                            SPECIAL CASE -> IMRG=1
C
C     Start loop over layers
C
      IF (IMRG.EQ.1) THEN
C
C        -----------------------------
C        Initial call to OPPATH, which calls PATH
C        -----------------------------
C
         LAYHDR = LAYER
         CALL OPPATH
         IF (IHIRAC.EQ.0) RETURN
C
C        -----------------------------
C        Begin loop over layers
C        -----------------------------
C
 1       LAYER = LAYER+1
         LAYHDR = LAYER
         CALL OPPATH
         NLAYHD = NLAYER
         CALL OPDPTH (MPTS)
         CALL ENDFIL (KFILE)
         REWIND MFILE
         REWIND LFILE
         IF (LAYER.EQ.NLAYER) RETURN
         GO TO 1
      ENDIF
C
C     ---------------------
C
C     For IMRG = 35,36,40,41 (those options which use
C     precalculated layer optical depths stored on different
C     files for radiative transfer), read in the pathname of
C     the layer optical depths and determine format for the
C     addition of the layer number suffix.
C
      IF (IMRG.GE.35) THEN
         READ (IRD,945) PATH1
         CALL QNTIFY(PATH1,HFORM1)
      ENDIF
C
C     ---------------------
C
C     If IMRG = 40 and IEMIT = 1, then precalculated optical depths
C     on multiple files are combined to produce total downwelling radiance
C     (from space to ground), written to TAPE12.  The results are
C     monochromatic.
C
C     If IMRG = 41 and IEMIT = 1, then precalculated optical depths
C     on multiple files are combined to produce total upwelling radiance
C     (from ground to space), written to TAPE12.  The results are
C     monochromatic.
C
C                                     SPECIAL CASE -> IMRG=40/41, IEMIT=1
C
      IF ((IMRG.GE.40).AND.(IEMIT.EQ.1)) THEN
C
C        -----------------------------
C        Obtain information from KFILE
C        -----------------------------
C
         CALL OPNRAD(1,1,PATH1,HFORM1)
         REWIND KFILE
         CALL BUFIN (KFILE,KEOF,FILDU1(1),NFHDRF)
         LTGNT = LTNSAV
         TBOUND = TMPBND
         EMISIV = BNDEMI(1)
         LH1 = LH1SAV
         LH2 = LH2SAV
         IPATHL = IPTHD1
         NLAYER = NLAYD1
C
C        Check for forced IPATHL, and set layer boundaries as needed
C
         IF (IMRG.EQ.40) THEN
            IF (LH2.NE.1) THEN
               LH1 = MAX(LH1,LH2)
               LH2 = 1
            ENDIF
            JPATHL = 1
         ELSE
            IF (LH1.NE.1) THEN
               LH2 = MAX(LH1,LH2)
               LH1 = 1
            ENDIF
            JPATHL = 3
         ENDIF
C
C        Start of loop over layers
C
         IF (2*(NLAYER/2).NE.NLAYER) GO TO 10
         MSTOR = MFILE
         MFILE = LFILE
         LFILE = MSTOR
 10      LAYER = 0
 11      LAYER = LAYER+1
         REWIND MFILE
         REWIND LFILE
         WRITE (IPR,900)
         WRITE (IPR,905)
         CALL OPNRAD(NLAYER,LAYER,PATH1,HFORM1)
         WRITE (IPR,905)
         CALL XMERGE (NPTS,LFILE,MFILE,JPATHL)
         IF (LAYER.EQ.NLAYER) RETURN
         MSTOR = MFILE
         MFILE = LFILE
         LFILE = MSTOR
C
C        END OF LOOP OVER LAYERS
C
         GO TO 11
      ENDIF
C
C     ---------------------
C
      MMRG = MOD(IMRG,10)                                                 A12610
C                                                                         A12620
C      TESTS FOR PRESTORE                                                 A12630
C                                                                         A12640
      IF (MMRG.EQ.2) GO TO 18                                             A12650
      IF (MMRG.EQ.5) GO TO 18                                             A12660
      IF (MMRG.EQ.6) GO TO 18                                             A12670
      IF (MMRG.EQ.8) GO TO 18                                             A12680
      IF (MMRG.EQ.9) GO TO 18                                             A12690
C                                                                         A12700
C     -------------------------------
C     Call OPPATH, which calls PATH
C     -------------------------------
C
      LAYHDR = LAYER                                                      A12710
      CALL OPPATH                                                         A12720
      IF (IHIRAC.EQ.0) GO TO 170                                          A12730
      GO TO 20                                                            A12740
C
C     -----------------------------
C     Obtain information from KFILE
C     -----------------------------
C                                                                         A12750
 18   REWIND KFILE                                                        A12760
      CALL BUFIN (KFILE,KEOF,FILDU1(1),NFHDRF)                            A12770
      LTGNT = LTNSAV                                                      A12780
      TBOUND = TMPBND                                                     A12790
      EMISIV = BNDEMI(1)                                                  A12800
      LH1 = LH1SAV                                                        A12810
      LH2 = LH2SAV                                                        A12820
      IPATHL = IPTHD1                                                     A12830
      NLAYER = NLAYD1                                                     A12840
C
C     -----------------------------
C     Standard Merge Options Follow
C     -----------------------------
C                                                                         A12850
C     CHECK FOR FORCED IPATHL, AND SET LAYER BOUNDARIES AS NEEDED         A12860
C                                                                         A12870
      IF (IMRG.EQ.12.OR.MMRG.EQ.5) THEN                                   A12880
         IF (LH2.NE.1) THEN                                               A12890
            LH1 = MAX(LH1,LH2)                                            A12900
            LH2 = 1                                                       A12910
         ENDIF                                                            A12920
      ENDIF                                                               A12930
      IF (IMRG.EQ.32.OR.MMRG.EQ.6) THEN                                   A12940
         IF (LH1.NE.1) THEN                                               A12950
            LH2 = MAX(LH1,LH2)                                            A12960
            LH1 = 1                                                       A12970
         ENDIF                                                            A12980
      ENDIF                                                               A12990
      IF (IMRG.EQ.22.OR.MMRG.EQ.8) THEN                                   A13000
         IF (IPATHL.NE.2) THEN                                            A13010
            LH1 = MAX(LH1,LH2)                                            A13020
            LH2 = LH1                                                     A13030
         ENDIF                                                            A13040
      ENDIF                                                               A13050
C                                                                         A13060
C       PRESTORE IMRG TEST, AND TEST FOR SCANNED OR FILTERED              A13070
C                                        WEIGHTING FUNCTIONS              A13080
C                                                                         A13090
   20 IF (MMRG.EQ.2) GO TO 50                                             A13100
      IF (MMRG.EQ.5) GO TO 90                                             A13110
      IF (MMRG.EQ.6) GO TO 50                                             A13120
      IF (MMRG.EQ.8) GO TO 110                                            A13130
      IF (MMRG.EQ.9) GO TO 150                                            A13140
      IF (IMRG.GT.10) GO TO 50                                            A13150
      IF (MMRG.GE.1) GO TO 40                                             A13160
C                                                                         A13170
C                                                                         A13180
C    START OF LOOP OVER LAYERS                    IMRG=0                  A13190
C                                                                         A13200
      IF (2*(NLAYER/2).NE.NLAYER) GO TO 30                                A13250
      MSTOR = MFILE                                                       A13260
      MFILE = LFILE                                                       A13270
      LFILE = MSTOR                                                       A13280
   30 LAYER = LAYER+1                                                     A13290
      REWIND KFILE                                                        A13300
      LAYHDR = LAYER                                                      A13310
      CALL OPPATH                                                         A13320
      NLAYHD = NLAYER                                                     A13330
      CALL OPDPTH (MPTS)                                                  A13340
      REWIND KFILE                                                        A13350
      REWIND MFILE                                                        A13360
      REWIND LFILE                                                        A13370
      CALL XMERGE (NPTS,LFILE,MFILE,JPATHL)                               A13380
      IF (LAYER.EQ.NLAYER) GO TO 170                                      A13390
      MSTOR = MFILE                                                       A13400
      MFILE = LFILE                                                       A13410
      LFILE = MSTOR                                                       A13420
C                                                                         A13430
C    END OF LOOP OVER LAYERS                                              A13440
C                                                                         A13450
      GO TO 30                                                            A13460
C                                                                         A13470
C     START LOOP OVER LAYERS                      IMRG=3,4,7              A13480
C     OPTICAL DEPTH STRUNG OUT ON TAPE10                                  A13490
C                                                                         A13500
   40 LAYER = LAYER+1                                                     A13510
      LAYHDR = LAYER                                                      A13520
      CALL OPPATH                                                         A13530
      NLAYHD = NLAYER                                                     A13540
      CALL OPDPTH (MPTS)                                                  A13550
      CALL ENDFIL (KFILE)                                                 A13560
      REWIND MFILE                                                        A13570
      REWIND LFILE                                                        A13580
      IF (LAYER.EQ.NLAYER) GO TO 50                                       A13590
      GO TO 40                                                            A13600
C                                                                         A13610
C    END OF LOOP OVER LAYERS                                              A13620
C                                                                         A13630
C                                    IMRG=2,4,6,12,14,16,22,24,26,32,36   A13640
C                                                                         A13650
   50 IF (MMRG.EQ.1) GO TO 170                                            A13660
      IF (MMRG.EQ.3) GO TO 90                                             A13670
      IF (MMRG.EQ.7) GO TO 110                                            A13680
C                                                                         A13690
C       READ CARD FOR SCAN OF WEIGHTING FUNCTION                          A13700
C       READ CARD FOR FILTER OF WEIGHTING FUNCTION                        A13710
C                                                                         A13720
      IF (IMRG.EQ.14.OR.IMRG.EQ.16.OR.IMRG.EQ.36)                         A13730
     *   CALL SCANRD (DVINT,JEMIT)                                        A13740
      IF (IMRG.EQ.24.OR.IMRG.EQ.26) CALL FLTRRD                           A13750
C                                                                         A13760
C     START LOOP OVER LAYERS                                              A13770
C     ALL MERGING DONE HERE                                               A13780
C                                                                         A13790
      REWIND KFILE                                                        A13800
      REWIND LFILE                                                        A13810
      REWIND MFILE                                                        A13820
C                                                                         A13830
C    DETERMINE IF FORCED IPATHL, AND SET APPROPRIATELY                    A13840
C                                                                         A13850
      IF (IMRG.NE.2) THEN                                                 A13860
         IF (IMRG.EQ.12) THEN                                             A13870
            JPATHL = 1                                                    A13880
         ELSE                                                             A13890
            JPATHL = 3                                                    A13900
         ENDIF                                                            A13910
      ENDIF                                                               A13920
      NNTAN = 1                                                           A13930
      NTAN(NNTAN) = 1                                                     A13940
      IF (2*(NLAYER/2).NE.NLAYER) GO TO 60                                A13950
      MSTOR = MFILE                                                       A13960
      MFILE = LFILE                                                       A13970
      LFILE = MSTOR                                                       A13980
   60 LAYER = 0                                                           A13990
   70 LAYER = LAYER+1                                                     A14000
      WRITE (IPR,900)                                                     A14010
      WRITE (IPR,905)                                                     A14020
      IF (MMRG.EQ.2) WRITE (IPR,910) NLAYER                               A14030
      IF (MMRG.EQ.4.OR.MMRG.EQ.6) THEN                                    A14040
         WRITE (IPR,915) NLAYER                                           A14050
         WRITE (IPR,920) (NTAN(N),N=1,NNTAN)                              A14060
      ENDIF                                                               A14070
      WRITE (IPR,905)                                                     A14080
      IF (IMRG.EQ.14.OR.IMRG.EQ.24) THEN                                  A14090
         REWIND KFILE                                                     A14100
         LAYHDR = LAYER                                                   A14110
         CALL OPPATH                                                      A14120
         NLAYHD = NLAYER                                                  A14130
         CALL OPDPTH (MPTS)                                               A14140
         REWIND KFILE                                                     A14150
      ENDIF                                                               A14160
      IF (IMRG.EQ.36) CALL OPNRAD(NLAYER,LAYER,PATH1,HFORM1)
      CALL XMERGE (NPTS,LFILE,MFILE,JPATHL)                               A14170
      NNTAN = NNTAN+1                                                     A14180
      NTAN(NNTAN) = NTAN(NNTAN-1)+1                                       A14190
      IF (MMRG.EQ.2) GO TO 80                                             A14200
      REWIND MFILE                                                        A14210
      IF (IMRG.EQ.4.OR.IMRG.EQ.6) CALL COPYFL (NPTS,MFILE,NFILE)          A14220
C                                                                         A14230
C   FOR SCAN CASE, IF DV NOT FINE ENOUGH, FIRST INTERPOLATE               A14240
C                                                                         A14250
      IF (IMRG.EQ.14.OR.IMRG.EQ.16.OR.IMRG.EQ.36) THEN                    A14260
         MMFILE = MFILE                                                   A14270
         IF (DVXM.GT.DVINT) THEN                                          A14280
            CALL SCNINT (MFILE,LFILE,DVINT,JEMIT,NPTS,IBUF)               A14290
            MMFILE = LFILE                                                A14300
         ENDIF                                                            A14310
         CALL SCNMRG (MMFILE,NFILE)                                       A14320
      ENDIF                                                               A14330
      IF (IMRG.EQ.24.OR.IMRG.EQ.26) CALL FLTMRG (MFILE,NFILE)             A14340
      CALL ENDFIL (NFILE)                                                 A14350
   80 IF (LAYER.EQ.NLAYER) GO TO 170                                      A14360
      MSTOR = MFILE                                                       A14370
      MFILE = LFILE                                                       A14380
      LFILE = MSTOR                                                       A14390
      REWIND MFILE                                                        A14400
      REWIND LFILE                                                        A14410
      GO TO 70                                                            A14420
C                                                                         A14430
C                                              IMRG=3,5,13,15,23,25,35    A14440
C                                                                         A14450
C       MODIFIED TO BEGIN WEIGHTING FUNCTION CALC. FROM H1                A14460
C                                                                         A14470
C                                                                         A14480
C     START LOOP OVER LAYERS                                              A14490
C     ALL MERGING DONE HERE                                               A14500
C                                                                         A14510
   90 REWIND KFILE                                                        A14520
      REWIND LFILE                                                        A14530
      REWIND MFILE                                                        A14540
C                                                                         A14550
C       READ CARD FOR SCAN OF WEIGHTING FUNCTION                          A14560
C       READ CARD FOR FILTER OF WEIGHTING FUNCTION                        A14570
C                                                                         A14580
      IF (IMRG.EQ.13.OR.IMRG.EQ.15.OR.IMRG.EQ.35)                         A14590
     *   CALL SCANRD (DVINT,JEMIT)                                        A14600
      IF (IMRG.EQ.23.OR.IMRG.EQ.25) CALL FLTRRD                           A14610
      NNTAN = 1                                                           A14620
      NTAN(NNTAN) = NLAYER                                                A14630
      IF (2*(NLAYER/2).EQ.NLAYER) THEN                                    A14640
         MSTOR = MFILE                                                    A14650
         MFILE = LFILE                                                    A14660
         LFILE = MSTOR                                                    A14670
      ENDIF                                                               A14680
C                                                                         A14690
C    SET VALUE FOR IPATHL                                                 A14700
C                                                                         A14710
      JPATHL = 1                                                          A14720
      LAYER = NLAYER+1                                                    A14730
  100 LAYER = LAYER-1                                                     A14740
      ISKIP = LAYER-1                                                     A14750
      WRITE (IPR,900)                                                     A14760
      WRITE (IPR,905)                                                     A14770
      WRITE (IPR,925) NLAYER,LAYER                                        A14780
      WRITE (IPR,920) (NTAN(N),N=1,NNTAN)                                 A14790
      WRITE (IPR,905)                                                     A14800
      IF (IMRG.EQ.13.OR.IMRG.EQ.23) THEN                                  A14810
         REWIND KFILE                                                     A14820
         LAYHDR = LAYER                                                   A14830
         CALL OPPATH                                                      A14840
         NLAYHD = NLAYER                                                  A14850
         CALL OPDPTH (MPTS)                                               A14860
         REWIND KFILE                                                     A14870
      ELSEIF (IMRG.EQ.35) THEN
         CALL OPNRAD(1,LAYER,PATH1,HFORM1)
      ELSE
         CALL SKIPFL (ISKIP,KFILE,IEOF)                                   A14890
      ENDIF                                                               A14900
      CALL XMERGI (NPTS,LFILE,MFILE,JPATHL)                               A14910
      NNTAN = NNTAN+1                                                     A14920
      NTAN(NNTAN) = NTAN(NNTAN-1)-1                                       A14930
      REWIND MFILE                                                        A14940
      IF (IMRG.EQ.3.OR.IMRG.EQ.5) CALL COPYFL (NPTS,MFILE,NFILE)          A14950
C                                                                         A14960
C     FOR SCAN CASE, IF DV NOT FINE ENOUGH, FIRST INTERPOLATE             A14970
C                                                                         A14980
      IF (IMRG.EQ.13.OR.IMRG.EQ.15.OR.IMRG.EQ.35) THEN                    A14990
         MMFILE = MFILE                                                   A15000
         IF (DVXM.GT.DVINT) THEN                                          A15010
            CALL SCNINT (MFILE,LFILE,DVINT,JEMIT,NPTS,IBUF)               A15020
            MMFILE = LFILE                                                A15030
         ENDIF                                                            A15040
         CALL SCNMRG (MMFILE,NFILE)                                       A15050
      ENDIF                                                               A15060
      IF (IMRG.EQ.23.OR.IMRG.EQ.25) CALL FLTMRG (MFILE,NFILE)             A15070
      CALL ENDFIL (NFILE)                                                 A15080
      IF (LAYER.EQ.1) GO TO 170                                           A15090
      MSTOR = MFILE                                                       A15100
      MFILE = LFILE                                                       A15110
      LFILE = MSTOR                                                       A15120
      REWIND KFILE                                                        A15130
      REWIND LFILE                                                        A15140
      REWIND MFILE                                                        A15150
      GO TO 100                                                           A15160
C                                                                         A15170
C   ** TANGENT **                                 IMRG=7,8,17,18,27,28    A15180
C                                                                         A15190
  110 REWIND KFILE                                                        A15200
      REWIND LFILE                                                        A15210
      REWIND MFILE                                                        A15220
C                                                                         A15230
C     READ CARD FOR SCAN OF WEIGHTING FUNCTION                            A15240
C     READ CARD FOR FILTER OF WEIGHTING FUNCTION                          A15250
C                                                                         A15260
      IF (IMRG.EQ.17.OR.IMRG.EQ.18) CALL SCANRD (DVINT,JEMIT)             A15270
      IF (IMRG.EQ.27.OR.IMRG.EQ.28) CALL FLTRRD                           A15280
C                                                                         A15290
C     **                 DOWN                                             A15300
C                                                                         A15310
      IANT = 1                                                            A15320
      NNTAN = 1                                                           A15330
      NTAN(NNTAN) = LH1                                                   A15340
      NLTOTL = LH1+LH2                                                    A15350
      IF (2*(NLTOTL/2).EQ.NLTOTL) THEN                                    A15360
         MSTOR = MFILE                                                    A15370
         MFILE = LFILE                                                    A15380
         LFILE = MSTOR                                                    A15390
      ENDIF                                                               A15400
C                                                                         A15410
C    SET VALUE FOR IPATHL                                                 A15420
C                                                                         A15430
      JPATHL = 3                                                          A15440
      LAYER = LH1+1                                                       A15450
  120 LAYER = LAYER-1                                                     A15460
      ISKIP = LAYER-1                                                     A15470
      WRITE (IPR,900)                                                     A15480
      WRITE (IPR,905)                                                     A15490
      WRITE (IPR,930) NLAYER,LAYER                                        A15500
      WRITE (IPR,920) (NTAN(N),N=1,NNTAN)                                 A15510
      WRITE (IPR,905)                                                     A15520
      IF (IMRG.EQ.17.OR.IMRG.EQ.27) THEN                                  A15530
         REWIND KFILE                                                     A15540
         LAYHDR = LAYER                                                   A15550
         CALL OPPATH                                                      A15560
         NLAYHD = NLAYER                                                  A15570
         CALL OPDPTH (MPTS)                                               A15580
         REWIND KFILE                                                     A15590
      ELSE                                                                A15600
         CALL SKIPFL (ISKIP,KFILE,IEOF)                                   A15610
      ENDIF                                                               A15620
      CALL XMERGI (NPTS,LFILE,MFILE,JPATHL)                               A15630
      NNTAN = NNTAN+1                                                     A15640
      NTAN(NNTAN) = NTAN(NNTAN-1)-1                                       A15650
      REWIND MFILE                                                        A15660
      IF (IMRG.EQ.7.OR.IMRG.EQ.8) CALL COPYFL (NPTS,MFILE,NFILE)          A15670
C                                                                         A15680
C     FOR SCAN CASE, IF DV NOT FINE ENOUGH, FIRST INTERPOLATE             A15690
C                                                                         A15700
      IF (IMRG.EQ.17.OR.IMRG.EQ.18) THEN                                  A15710
         MMFILE = MFILE                                                   A15720
         IF (DVXM.GT.DVINT) THEN                                          A15730
            CALL SCNINT (MFILE,LFILE,DVINT,JEMIT,NPTS,IBUF)               A15740
            MMFILE = LFILE                                                A15750
         ENDIF                                                            A15760
         CALL SCNMRG (MMFILE,NFILE)                                       A15770
      ENDIF                                                               A15780
      IF (IMRG.EQ.27.OR.IMRG.EQ.28) CALL FLTMRG (MFILE,NFILE)             A15790
      CALL ENDFIL (NFILE)                                                 A15800
      MSTOR = MFILE                                                       A15810
      MFILE = LFILE                                                       A15820
      LFILE = MSTOR                                                       A15830
      REWIND KFILE                                                        A15840
      REWIND LFILE                                                        A15850
      REWIND MFILE                                                        A15860
      IF (LAYER.EQ.1) GO TO 130                                           A15870
      GO TO 120                                                           A15880
C                                                                         A15890
C     **                 UP                                               A15900
C                                                                         A15910
  130 LAYER = 0                                                           A15920
      IANT = -1                                                           A15930
      NTAN(NNTAN) = 1                                                     A15940
  140 LAYER = LAYER+1                                                     A15950
      WRITE (IPR,900)                                                     A15960
      WRITE (IPR,905)                                                     A15970
      WRITE (IPR,930) NLAYER,LAYER                                        A15980
      WRITE (IPR,920) (NTAN(N),N=1,NNTAN)                                 A15990
      WRITE (IPR,905)                                                     A16000
      IF ((IMRG.EQ.17.OR.IMRG.EQ.27).AND.LAYER.NE.1) THEN                 A16010
         REWIND KFILE                                                     A16020
         LAYHDR = LAYER                                                   A16030
         CALL OPPATH                                                      A16040
         NLAYHD = NLAYER                                                  A16050
         CALL OPDPTH (MPTS)                                               A16060
         REWIND KFILE                                                     A16070
      ENDIF                                                               A16080
      IF (LAYER.LE.LH1) CALL XMERGI (NPTS,LFILE,MFILE,JPATHL)             A16090
      IF (LAYER.GT.LH1) CALL XMERGE (NPTS,LFILE,MFILE,JPATHL)             A16100
      NNTAN = NNTAN+1                                                     A16110
      NTAN(NNTAN) = NTAN(NNTAN-1)+1                                       A16120
      REWIND MFILE                                                        A16130
      IF (IMRG.EQ.7.OR.IMRG.EQ.8) CALL COPYFL (NPTS,MFILE,NFILE)          A16140
C                                                                         A16150
C     FOR SCAN CASE, IF DV NOT FINE ENOUGH, FIRST INTERPOLATE             A16160
C                                                                         A16170
      IF (IMRG.EQ.17.OR.IMRG.EQ.18) THEN                                  A16180
         MMFILE = MFILE                                                   A16190
         IF (DVXM.GT.DVINT) THEN                                          A16200
            CALL SCNINT (MFILE,LFILE,DVINT,JEMIT,NPTS,IBUF)               A16210
            MMFILE = LFILE                                                A16220
         ENDIF                                                            A16230
         CALL SCNMRG (MMFILE,NFILE)                                       A16240
      ENDIF                                                               A16250
      IF (IMRG.EQ.27.OR.IMRG.EQ.28) CALL FLTMRG (MFILE,NFILE)             A16260
      CALL ENDFIL (NFILE)                                                 A16270
      IF (LAYER.EQ.LH2) GO TO 170                                         A16280
      MSTOR = MFILE                                                       A16290
      MFILE = LFILE                                                       A16300
      LFILE = MSTOR                                                       A16310
      REWIND LFILE                                                        A16320
      REWIND MFILE                                                        A16330
      GO TO 140                                                           A16340
C                                                                         A16350
C     START LOOP OVER LAYERS                      IMRG=9                  A16360
C     OPTICAL DEPTH STRUNG OUT ON TAPE10                                  A16370
C                                                                         A16380
C                                                                         A16390
C     START LOOP OVER LAYERS                                              A16400
C     ALL MERGING DONE HERE                                               A16410
C                                                                         A16420
  150 REWIND MFILE                                                        A16430
      IF (IAERSL.NE.0) THEN                                               A16440
         REWIND 20                                                        A16450
         IREAD = 0                                                        A16460
         LOWFLG = 3                                                       A16470
         LAYHDR = LAYER                                                   A16480
         CALL OPPATH                                                      A16490
      ENDIF                                                               A16500
      REWIND LFILE                                                        A16510
      REWIND KFILE                                                        A16520
      IF (2*(NLAYER/2).NE.NLAYER) GO TO 160                               A16530
      MSTOR = MFILE                                                       A16540
      MFILE = LFILE                                                       A16550
      LFILE = MSTOR                                                       A16560
      LAYER = 0                                                           A16570
  160 LAYER = LAYER+1                                                     A16580
      WRITE (IPR,900)                                                     A16590
      CALL XMERGE (NPTS,LFILE,MFILE,JPATHL)                               A16600
      IF (LAYER.EQ.NLAYER) GO TO 170                                      A16610
      MSTOR = MFILE                                                       A16620
      MFILE = LFILE                                                       A16630
      LFILE = MSTOR                                                       A16640
      REWIND MFILE                                                        A16650
      REWIND LFILE                                                        A16660
      GO TO 160                                                           A16670
C                                                                         A16680
  170 IF (IMRG.GE.23.AND.IMRG.LE.28) CALL FLTPRT (NFILE)                  A16690
C                                                                         A16700
      RETURN                                                              A16710
C                                                                         A16720
  900 FORMAT ('1')                                                        A16730
  905 FORMAT (/,'*******************************************',/)          A16740
  910 FORMAT ('0  SEQUENTIAL OPTICAL DEPTHS TO KFILE, MERGED TO MFILE',   A16750
     *        I3,' LAYERS')                                               A16760
  915 FORMAT ('0  GROUND TO SPACE WEIGHTING FUNCTION, LAYER 1 TO LAYER'   A16770
     *        ,I3)                                                        A16780
  920 FORMAT (/,'  THE WEIGHTING FUNCTION ACCUMULATION BY LAYER',/,       A16790
     *        (5X,20I3))                                                  A16800
  925 FORMAT (/,'  SPACE TO GROUND WEIGHTING FUNCTION, LAYER',I3,         A16810
     *        ' TO LAYER',I3)                                             A16820
  930 FORMAT (/,'  TANGENT WEIGHTING FUNCTION, LAYER',I3,' TO LAYER',     A16830
     *        I3)                                                         A16840
  940 FORMAT ('TAPE5: IPATHL',I5,', IMRG = ',I5)
  945 FORMAT (A55)
C                                                                         A16850
      END                                                                 A16860
      SUBROUTINE OPPATH                                                   A16870
C                                                                         A16880
      IMPLICIT DOUBLE PRECISION (V)                                     ! A16890
C                                                                         A16900
C     OPPATH CALLS LBLATM AND CALLS PATH FIRST                            A16910
C                                                                         A11260
      PARAMETER (MXFSC=200,MXLAY=MXFSC+3,MXZMD=200,MXPDIM=MXLAY+MXZMD,
     *                IM2=MXPDIM-2,MXMOL=35,MXTRAC=22)
C                                                                         A16920
      COMMON /PATHD/ PAVEL(MXLAY),TAVEL(MXLAY),WKL(35,MXLAY),
     *               WBRODL(MXLAY),DVL(MXLAY),                            A16930
     *               WTOTL(MXLAY),ALBL(MXLAY),ADBL(MXLAY),
     *               AVBL(MXLAY),H2OSL(MXLAY),                            A16940
     *               IPTH(MXLAY),ITYL(MXLAY),SECNTA(MXLAY),
     *               HT1,HT2,ALTZ(0:MXLAY),                               A16950
     *               PZ(0:MXLAY),TZ(0:MXLAY)                              A16960
C                                                                         A16970
C     IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE     A16980
C     MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES       A16990
C     (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR     A17000
C     L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD         A17010
C     MOLECULES.                                                          A17020
C                                                                         A17030
      COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(35),XAMNT(35,MXLAY)              A17040
C                                                                         A17050
C     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES         A17060
C     FOR THE CROSS-SECTION MOLECULES.                                    A17070
C                                                                         A17080
      COMMON /XSECTR/ V1FX(5,35),V2FX(5,35),DVFX(5,35),WXM(35),           A17090
     *                NTEMPF(5,35),NSPECR(35),IXFORM(5,35),               A17100
     *                XSMASS(35),XDOPLR(5,35),NUMXS,IXSBIN                A17105
C                                                                         A17110
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   A17120
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       A17130
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       A17140
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      A17150
     *              EXTID(10)                                             A17160
C                                                                         A17170
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & A17180
C                                                                         A17190
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       A17200
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   A17210
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    A17220
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP                   A17230
C                                                                         A17240
      DOUBLE PRECISION  HLINID,BMOLID,HID1                              & A17250
C                                                                         A17260
      COMMON /LINHDR/ HLINID(10),BMOLID(64),MOLCNT(64),MCNTLC(64),        A17270
     *                MCNTNL(64),SUMSTR(64),LINMOL,FLINLO,FLINHI,         A17280
     *                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1(2),LSTWDL     A17290
C                                                                         A17300
C     LSTWD (LAST WORD) IS DUMMY, DOES NOT NEED TO BE COUNTED             A17310
C                                                                         A17320
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         A17330
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        A17340
     *              NLTEFL,LNFIL4,LNGTH4                                  A17350
      COMMON /MSCONS/ AIRMAS(MXLAY),TGRND,SEMIS(3),HMINMS,HMAXMS,
     *                MSFLAG,                                             A17360
     *                MSWIT,IODFIL,MSTGLE                                 A17370
      COMMON /LASIV/ VLAS,ILAS                                            A17380
      COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,NOZERO,NP,H1F,H2F,         A17390
     *                ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH,         A17400
     *                XVBAR, HMINF,PHIF,IERRF,HSPACE                      A17410
      COMMON /CNTRL/ I1,I2,I3,I4,NBNDL,I6,I7,NBNDF,I9                     A17420
C                                                                         A17430
C     NEW EQUIVALENCE STATEMENTS FOR TANGENT WEIGHTING FNS.               A17440
C                                                                         A17450
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4)                 A17460
      EQUIVALENCE (FSCDID(5),IEMIT) , (FSCDID(4),IAERSL)                  A17470
      EQUIVALENCE (FSCDID(15),IATM) , (FSCDID(8),IPATHL)                  A17480
      EQUIVALENCE (YID(10),LTNSAV) , (YID(9),LH1SAV) , (YID(8),LH2SAV)    A17490
      EQUIVALENCE (XALTZ(1),ALTZL) , (XALTZ(2),ALTZU)                     A17500
C                                                                         A17510
      DATA XDV / 64. /                                                    A17520
      DATA HZ / 6.5 /                                                     A17530
C                                                                         A17540
      NBNDF = 0                                                           A17550
C                                                                         A17560
C     SET AV1 AND AV2 FOR PASS TO LBLATM AND LOWTRN                       A17570
C                                                                         A17580
      AV1 = V1                                                            A17590
      AV2 = V2                                                            A17600
C                                                                         A17610
C     RESET ALFAL0 TO DEFAULT OF 0.04 IN VALUE READ IN .LE. 0.0
C     FOR LINE COUPLING USER SHOULD READ IN VALUE OF 0.05
C
      IF (SAMPLE.LT.0.01) SAMPLE = 4.
      IF (ALFAL0.LE.0.) ALFAL0 = 0.04
      IF (AVMASS.LE.0.) AVMASS = 36.
C
C     Write SAMPLE, ALFAL0, and AVMASS if IHIRAC = 0 and IATM >= 1
C     (printed for informational puroses in case IFXTYP = 1).
C
      IF ((IHIRAC.EQ.0).AND.(IATM.GE.1))
     *     WRITE (IPR,890) SAMPLE,ALFAL0,AVMASS
C
      IF (LAYER.LE.0) THEN                                                A17620
         IF (IAERSL.EQ.0.OR.IAERSL.EQ.9) THEN                             A17640
C                                                                         A17650
C           IF LAYER GT 0 SAVE THE VECTORS                                A17660
C                                                                         A17670
            IF (IATM.GE.1.AND.IATM.LE.5) THEN                             A17680
C                                                                         A17690
C              CALL LBLATM TO COMPLETE GEOMETRY ON FINAL COMBINED LAYERS  A17700
C              IF AEROSOLS PRESENT AND HORIZONTAL PATH SKIP CALL          A17710
C              TO LBLATM                                                  A17720
C                                                                         A17730
               CALL LBLATM                                                A17740
               IF (IHIRAC.EQ.0) RETURN                                    A17750
            ENDIF                                                         A17770
         ELSE                                                             A17780
C                                                                         A17790
            IF (IATM.EQ.0) THEN                                           A17800
               WRITE (IPR,900)                                            A17810
               STOP 'OPPATH; IAER GT 1, IATM = 0 '                        A17820
            ENDIF                                                         A17830
C                                                                         A17840
            IREAD = 0                                                     A17850
C                                                                         A17860
            CALL LBLATM                                                   A17870
C                                                                         A17880
            NBNDF = NLAYRS+1                                              A17890
C                                                                         A17900
C     STORE V1 AND V2 IN V1S AND V2S DURING CALL TO LOWTRN                A17910
C                                                                         A17920
            V1S = V1                                                      A17930
            V2S = V2                                                      A17940
C                                                                         A17950
            CALL LOWTRN                                                   A17960
C                                                                         A17970
C     RESET V1 AND V2 TO V1S AND V2S                                      A17980
C                                                                         A17990
            V1 = V1S                                                      A18000
            V2 = V2S                                                      A18010
C                                                                         A18020
         ENDIF                                                            A18030
C                                                                         A18040
         IF (IATM.LE.5) THEN                                              A18050
C                                                                         A18060
            WRITE (IPR,905)                                               A18070
C                                                                         A18080
            WRITE (IPR,910) V1,V2,SAMPLE,DVSET,ALFAL0,AVMASS,DPTMIN,      A18160
     *                      DPTFAC                                        A18170
C                                                                         A18180
C     BEGINNING AND ENDING WAVENUMBER VALUES                              A18190
C                                                                         A18200
            IF ((V2-V1).GT.2020.) THEN                                    A18210
               WRITE (IPR,915)                                            A18220
               STOP 'OPPATH; V2-V1 GT 2020'                               A18230
            ENDIF                                                         A18240
C                                                                         A18250
            IF ((IHIRAC.EQ.1.OR.IHIRAC.EQ.4.OR.ILBLF4.GE.1).AND.          A18260
     *           (IATM.EQ.1)) THEN                                        A18270
               IF (NMOL.LT.LINMOL) WRITE (IPR,920) (LINMOL-NMOL)          A18280
               IF (NMOL.GT.LINMOL) WRITE (IPR,925) NMOL,LINMOL            A18290
            ENDIF                                                         A18300
C                                                                         A18310
            CALL PATH                                                     A18320
C                                                                         A18330
C  SAVE AIRMASS FACTORS FOR USE WITH MULTIPLE SCATTERING                  A18340
C                                                                         A18350
            DO 10 IAIR = 1, 67                                            A18360
               AIRMAS(IAIR) = SECNTA(IAIR)                                A18370
   10       CONTINUE                                                      A18380
         ENDIF                                                            A18390
         LTNSAV = LTGNT                                                   A18400
         LH2SAV = LH2                                                     A18410
C                                                                         A18420
         NMOLIN = NMOL+1                                                  A18430
         DO 20 MOL = NMOLIN, 35                                           A18440
            DO 20 ILAYR = 1, NLAYRS                                       A18450
               WKL(MOL,ILAYR) = 0.                                        A18460
   20    CONTINUE                                                         A18470
C                                                                         A18480
         RETURN                                                           A18490
      ENDIF                                                               A18500
C                                                                         A18510
      REWIND LINFIL                                                       A18520
C                                                                         A18530
      IF (ILAS.GT.0) THEN                                                 A18540
         DVI = XDV*DVL(1)                                                 A18550
         MM = 64                                                          A18560
         DVC = FLOAT(MM)*DVL(LAYER)                                       A18570
         V2 = VLAS+DVC                                                    A18580
         V1 = VLAS-DVC                                                    A18590
C                                                                         A18600
CPRT     WRITE(IPR,930) V1,V2                                             A18610
C                                                                         A18620
      ENDIF                                                               A18630
C                                                                         A18640
C     SET UP LAYER BOUNDARY PARAMETERS                                    A18650
C                                                                         A18660
      ALTZL = ALTZ(LAYER-1)                                               A18670
      ALTZU = ALTZ(LAYER)                                                 A18680
      PZL = PZ(LAYER-1)                                                   A18690
      PZU = PZ(LAYER)                                                     A18700
      TZL = TZ(LAYER-1)                                                   A18710
      TZU = TZ(LAYER)                                                     A18720
C                                                                         A18730
      PAVE = PAVEL(LAYER)                                                 A18740
      TAVE = TAVEL(LAYER)                                                 A18750
      WBROAD = WBRODL(LAYER)                                              A18760
C                                                                         A18770
      DO 30 M = 1, NMOL                                                   A18780
         WK(M) = WKL(M,LAYER)                                             A18790
         WXM(M) = XAMNT(M,LAYER)                                          A18800
   30 CONTINUE                                                            A18810
C                                                                         A18820
      H2OSLF = H2OSL(LAYER)                                               A18830
      WTOT = WTOTL(LAYER)                                                 A18840
      ALBAR = ALBL(LAYER)                                                 A18850
      ADBAR = ADBL(LAYER)                                                 A18860
      AVBAR = AVBL(LAYER)                                                 A18870
      DV = DVL(LAYER)                                                     A18880
      IPATHL = IPTH(LAYER)                                                A18890
      SECNTL = SECNTA(LAYER)                                              A18900
      SECANT = SECNTA(LAYER)                                              A18910
      IF (ALTZ(LAYER-1).GE.0.) THEN                                       A18920
         ALTAV = ALTZ(LAYER-1)-                                           A18930
     *      HZ*ALOG(.5*(1.+EXP(-(ALTZ(LAYER)-ALTZ(LAYER-1))/HZ)))         A18940
      ELSE                                                                A18950
         ALTAV = ALTZ(LAYER)                                              A18960
      ENDIF                                                               A18970
C                                                                         A18980
      RETURN                                                              A18990
C                                                                         A19000
  890 FORMAT (//,
     *        '0 SAMPLE   =',F13.4,/,
     *        '0 ALFAL0   =',F13.4,/,
     *        '0 AVMASS   =',F13.4,/)
  900 FORMAT ('1 ERROR IATM = 0 IAERSL GT 0 ')                            A19010
  905 FORMAT ('1')                                                        A19020
  910 FORMAT ('0 V1(CM-1) = ',F12.4,/'0 V2(CM-1) = ',F12.4,/,             A19030
     *        '0 SAMPLE   =',F13.4,/'0 DVSET    =',F13.6,/,               A19040
     *        '0 ALFAL0   =',F13.4,/'0 AVMASS   =',F13.4,/,               A19050
     *        '0 DPTMIN   =',1P,E13.4,13X,'  DPTFAC   =',0P,F13.6)        A19060
  915 FORMAT ('0 V2-V1 .GT. 2000. ')                                      A19070
  920 FORMAT ('0',1X,7('*'),' LAST ',I5,' MOLECULES ON LINFIL NOT ',      A19080
     *        'SELECTED')                                                 A19090
  925 FORMAT ('0',1X,53('*'),/,'0',1X,14('*'),I5,' MOLECULES ',           A19100
     *        'REQUESTED',2X,12('*'),/,2X,7('*'),2X,'ONLY ',I5,           A19110
     *        ' MOLECULES ON LINFIL',2X,11('*'),/,'0',1X,53('*'))         A19120
  930 FORMAT (2(/),'  V1 RESET ',F10.3,'  V2 RESET ',F10.3)               A19130
C                                                                         A19140
      END                                                                 A19150
      SUBROUTINE PATH                                                     A19160
C                                                                         A19170
      IMPLICIT DOUBLE PRECISION (V)                                     ! A19180
C                                                                         A19190
      PARAMETER (NTMOL=32,NSPECI=75)                                      A19200
C                                                                         A19210
C     SUBROUTINE PATH INITIALIZES LINFIL AND INPUTS LAYER PARAMETERS      A19220
C     SUBROUTINE PATH INPUTS AND OUTPUTS HEADER FROM LINFIL AND           A19230
C     INPUTS AND OUTPUTS PATH PARAMETERS FOR EACH LAYER                   A19240
C                                                                         A19250
      PARAMETER (MXFSC=200,MXLAY=MXFSC+3,MXZMD=200,MXPDIM=MXLAY+MXZMD,
     *                IM2=MXPDIM-2,MXMOL=35,MXTRAC=22)
C
      COMMON COMSTR(250,9)                                                A19260
      COMMON R1(3600),R2(900),R3(225)                                     A19270
C
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   A19280
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       A19290
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       A19300
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      A19310
     *              EXTID(10)                                             A19320
      COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),
     *                MSPNL1(MXLAY),MSLAY1,ISFILE,JSFILE,KSFILE,
     *                LSFILE,MSFILE,IEFILE,JEFILE,KEFILE
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           A19330
C                                                                         A19340
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & A19350
C                                                                         A19360
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       A19370
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   A19380
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    A19390
      COMMON /PATHD/ PAVEL(MXLAY),TAVEL(MXLAY),WKL(35,MXLAY),
     *               WBRODL(MXLAY),DVL(MXLAY),                            A19400
     *               WTOTL(MXLAY),ALBL(MXLAY),ADBL(MXLAY),
     *               AVBL(MXLAY),H2OSL(MXLAY),                            A19410
     *               IPTH(MXLAY),ITYL(MXLAY),SECNTA(MXLAY),
     *               HT1,HT2,ALTZ(0:MXLAY),                               A19420
     *               PZ(0:MXLAY),TZ(0:MXLAY)                              A19430
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         A19440
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        A19450
     *              NLTEFL,LNFIL4,LNGTH4                                  A19460
C                                                                         A19470
C     IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE     A19480
C     MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES       A19490
C     (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR     A19500
C     L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD         A19510
C     MOLECULES.                                                          A19520
C                                                                         A19530
      COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(35),XAMNT(35,MXLAY)              A19540
C                                                                         A19550
C     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES         A19560
C     FOR THE CROSS-SECTION MOLECULES.                                    A19570
C                                                                         A19580
      CHARACTER*10 XSFILE,XSNAME,ALIAS                                    A19590
      COMMON /XSECTF/ XSFILE(6,5,35),XSNAME(35),ALIAS(4,35)               A19600
      COMMON /XSECTR/ V1FX(5,35),V2FX(5,35),DVFX(5,35),WXM(35),           A19610
     *                NTEMPF(5,35),NSPECR(35),IXFORM(5,35),               A19620
     *                XSMASS(35),XDOPLR(5,35),NUMXS,IXSBIN                A19625
      COMMON /IODFLG/ DVOUT
C                                                                         A19630
      CHARACTER*20 HEAD20
      CHARACTER*6 MOLID                                                   A19640
      COMMON /MOLNAM/ MOLID(0:NTMOL)                                      A19650
      CHARACTER*7 HEAD7
      CHARACTER*6 HOLN2                                                   A19660
      CHARACTER*5 HEAD5
      CHARACTER*4 HEAD4
      CHARACTER*4 HT1HRZ,HT2HRZ,HT1SLT,HT2SLT                             A19670
      CHARACTER*3 CINP,CINPX,CBLNK                                        A19680
      DIMENSION FILHDR(2),AMOUNT(2),AMTSTR(2)                             A19690
      DIMENSION HEDXS(15),WMT(35),SECL(64),WXT(35),WTOTX(67)              A19700
C                                                                         A19710
      EQUIVALENCE (XID(1),FILHDR(1))                                      A19720
      EQUIVALENCE (YID(10),LTNSAV) , (YID(9),LH1SAV) , (YID(8),LH2SAV)    A19730
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),                A19740
     *            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),                A19750
     *            (FSCDID(5),IEMIT) , (FSCDID(7),IPLOT),                  A19760
     *            (FSCDID(9),JRAD) , (FSCDID(11),IMRG),                   A19770
     *            (FSCDID(15),IATM)                                       A19780
      EQUIVALENCE (AMOUNT(1),COMSTR(1,1)) , (AMTSTR(1),COMSTR(21,1))      A19790
C                                                                         A19800
      DATA HOLN2 / ' OTHER'/                                              A19810
      DATA HT1HRZ / ' AT '/,HT2HRZ / ' KM '/,HT1SLT / ' TO '/,            A19820
     *     HT2SLT / ' KM '/                                               A19830
      DATA CBLNK / '   '/                                                 A19840
C                                                                         A19850
      IF (IHIRAC.EQ.0) RETURN                                             A19860
C                                                                         A19870
      ICNTNM = MOD(IXSCNT,10)                                             A19880
      IXSECT = IXSCNT/10                                                  A19890
C                                                                         A19900
      ISET = 0                                                            A19910
      IDVSET = 0                                                          A19920
      IF (DVSET.LT.0.) THEN                                               A19930
         DVSET = -DVSET                                                   A19940
         ISET = 1                                                         A19950
      ENDIF                                                               A19960
C                                                                         A19970
      DO 10 M = 1, 35                                                     A19980
         WMT(M) = 0.                                                      A19990
         WK(M) = 0.                                                       A20000
         WXM(M) = 0.                                                      A20010
         WXT(M) = 0.                                                      A20020
   10 CONTINUE                                                            A20030
C                                                                         A20040
      SUMN2 = 0.                                                          A20050
      ISTOP = 0                                                           A20060
C                                                                         A20070
      SECNT0 = 1.                                                         A20080
      IPATHL = 1                                                          A20090
      LH1 = 1                                                             A20100
      LH2 = 1                                                             A20110
C                                                                         A20120
C     Read in atmospheric definition information
C
      IF (IATM.EQ.0) THEN                                                 A20130
         READ (IRD,901) IFORM,NLAYRS,NMOL,SECNT0,HEAD20,ZH1,HEAD4,ZH2,
     *                  HEAD5,ZANGLE,HEAD7                                A20140
         IF (NMOL.EQ.0) NMOL = 7                                          A20150
         IF (SECNT0.LT.0.) THEN                                           A20160
            IPATHL = 1                                                    A20170
         ELSE                                                             A20180
            IPATHL = 3                                                    A20190
         ENDIF                                                            A20200
         SECNT0 = ABS(SECNT0)                                             A20210
         WRITE (IPR,902) SECNT0,NLAYRS,NMOL,HEAD20,ZH1,HEAD4,ZH2,
     *                   HEAD5,ZANGLE,HEAD7                               A20220
C
C        Put H1, H2, and ANGLE into YID (ANGLE is needed for
C        CHARTS multiple scattering calculation)
C
         CALL YDIH1(ZH1,ZH2,ZANGLE,YID)
         IF (IHIRAC.EQ.9) THEN                                            A20230
            CALL VECISO                                                   A20240
            DO 20 M = 1, NMOL                                             A20250
               READ (MOLID(M),905) HMOLID(M)                              A20260
   20       CONTINUE                                                      A20270
         ENDIF                                                            A20280
      ELSE                                                                A20290
         WRITE (IPR,907) SECNT0                                           A20300
         JTYPE = 0                                                         PDB
         GO TO 50                                                         A20310
      ENDIF                                                               A20320
C                                                                         A20330
      JTYPE = 0                                                           A20340
      DO 30 L = 1, NLAYRS                                                 A20350
C                                                                         A20360
         IF (L.EQ.1) THEN                                                 A20370
            IF (IFORM.EQ.1) THEN
               READ (IRD,910) PAVE,TAVE,SECNTK,CINP,IPTHRK,ALTZ(L-1),     A20380
     *              PZ(L-1),TZ(L-1),ALTZ(L),PZ(L),TZ(L)                   A20390
            ELSE
               READ (IRD,911) PAVE,TAVE,SECNTK,CINP,IPTHRK,ALTZ(L-1),
     *              PZ(L-1),TZ(L-1),ALTZ(L),PZ(L),TZ(L)
            ENDIF
            IF (CINP.NE.CBLNK) WRITE (IPR,912)                            A20400
         ELSE                                                             A20410
            IF (IFORM.EQ.1) THEN
               READ (IRD,915) PAVE,TAVE,SECNTK,CINP,IPTHRK,
     *              ALTZ(L),PZ(L),                                        A20420
     *              TZ(L)                                                 A20430
            ELSE
               READ (IRD,916) PAVE,TAVE,SECNTK,CINP,IPTHRK,
     *              ALTZ(L),PZ(L),
     *              TZ(L)
            ENDIF
            IF ((CINP.EQ.CBLNK.AND.JTYPE.EQ.1).OR.                        A20440
     *          (CINP.NE.CBLNK.AND.JTYPE.EQ.0)) THEN                      A20450
               WRITE (IPR,912)                                            A20460
               STOP ' JTYPE ERROR IN PATH '                               A20470
            ENDIF                                                         A20480
         ENDIF                                                            A20490
         IF (IPTHRK.EQ.0.AND.SECNT0.EQ.1.) THEN                           A20500
            READ (HT1HRZ,917) HT1                                         A20510
            READ (HT2HRZ,917) HT2                                         A20520
         ELSE                                                             A20530
            READ (HT1SLT,917) HT1                                         A20540
            READ (HT2SLT,917) HT2                                         A20550
         ENDIF                                                            A20560
         IF (CINP.NE.CBLNK) THEN                                          A20570
            JTYPE = 1                                                     A20580
            READ (CINP,920) ITYL(L)                                       A20590
         ENDIF                                                            A20600
C                                                                         A20610
         PAVEL(L) = PAVE                                                  A20620
         TAVEL(L) = TAVE                                                  A20630
         SECANT = SECNT0                                                  A20640
         IF (SECNTK.GT.0.) SECANT = SECNTK                                A20650
         SECL(L) = SECANT                                                 A20660
         SECNTA(L) = SECANT                                               A20670
         IF (IPTHRK.NE.0) IPATHL = IPTHRK                                 A20680
         IPTH(L) = IPATHL                                                 A20690
         IF (SECANT.EQ.0.) STOP 'PATH; SECANT = 0'                        A20700
         IF (IFORM.EQ.1) THEN
            READ (IRD,925) (WKL(M,L),M=1,7),WBRODL(L)                     A20720
            IF (NMOL.GT.7) READ (IRD,925) (WKL(M,L),M=8,NMOL)             A20730
         ELSE
            READ (IRD,927) (WKL(M,L),M=1,7),WBRODL(L)
            IF (NMOL.GT.7) READ (IRD,927) (WKL(M,L),M=8,NMOL)
         ENDIF
C
C     --------------------------------------------------------------
C
C                     MIXING RATIO INPUT
C
C
C     First calculate the column amount of dry air ("WDRAIR")
C     Initialize WDRAIR to WBRODL(L) (automatically in column density)
C     Determine if each molecule is in column density.
C        - if so, just add to WDRAIR
C        - if not, convert mixing ratio on the fly in determining WDRAIR
C
C     NOTE that if WKL is greater than one, then column density
C               if WKL is less than or equal to one, then mixing ratio
C
         WDRAIR = WBRODL(L)
         DO 22 M = 2,NMOL
            IF (WKL(M,L).GT.1) THEN
               WDRAIR = WDRAIR + WKL(M,L)
            ELSEIF (WKL(M,L).LT.1) THEN
               WDRAIR = WDRAIR + WKL(M,L)*WBRODL(L)/(1-WKL(M,L))
            ELSE
               WDRAIR = WBRODL(L)
            ENDIF
 22      CONTINUE
C
C     NOW CONVERT ALL OTHER MOLECULES WHICH MAY BE IN MIXING RATIO
C     TO MOLECULAR DENSITY USING WDRAIR
C
         DO 25 M = 1,NMOL
            IF (WKL(M,L).LE.1) WKL(M,L) = WKL(M,L)*WDRAIR
 25      CONTINUE
C
C     --------------------------------------------------------------
C
   30 CONTINUE                                                            A20740
C                                                                         A20750
      IF (IATM.EQ.0.AND.IXSECT.GE.1) THEN                                 A20760
         READ (IRD,930) IXMOLS,IXSBIN                                     A20770
         XV1 = V1                                                         A20780
         XV2 = V2                                                         A20790
         CALL XSREAD (XV1,XV2)                                            A20800
         WRITE (IPR,932) (I,XSNAME(I),I=1,IXMOLS)                         A20810
         READ (IRD,900) IFRMX,NLAYXS,IXMOL,SECNTX,HEDXS                   A20820
         IF (IXMOL.EQ.0) THEN                                             A20830
            WRITE (IPR,935) IXMOL                                         A20840
            STOP ' PATH - IXMOL 0 '                                       A20850
         ENDIF                                                            A20860
         IF (IXMOL.NE.IXMOLS) THEN                                        A20870
            WRITE (IPR,937) IXMOL,IXMOLS                                  A20880
            STOP ' PATH - IXMOL .NE. IXMOLS '                             A20890
         ENDIF                                                            A20900
         IF (NLAYRS.NE.NLAYXS) THEN                                       A20910
            WRITE (IPR,940) NLAYRS,NLAYXS                                 A20920
            STOP ' PATH - NLAYRS .NE. NLAYXS '                            A20930
         ENDIF                                                            A20940
         SECNTX = ABS(SECNTX)                                             A20950
         WRITE (IPR,942) SECNTX,NLAYXS,IXMOLS,HEDXS                       A20960
C                                                                         A20970
         DO 40 L = 1, NLAYXS                                              A20980
C                                                                         A20990
            IF (L.EQ.1) THEN                                              A21000
               IF (IFRMX.EQ.1) THEN
                  READ (IRD,910) PAVX,TAVX,SECKXS,CINPX,IPTHKX,ALTXB,     A21010
     *                 PZXB,TZXB,ALTXT,PZXT,TZXT                          A21020
               ELSE
                  READ (IRD,911) PAVX,TAVX,SECKXS,CINPX,IPTHKX,ALTXB,
     *                 PZXB,TZXB,ALTXT,PZXT,TZXT
               ENDIF
            ELSE                                                          A21030
               IF (IFRMX.EQ.1) THEN
                  READ (IRD,915) PAVX,TAVX,SECKXS,CINPX,IPTHKX,ALTXT,     A21040
     *                 PZXT,TZXT                                          A21050
               ELSE
                  READ (IRD,916) PAVX,TAVX,SECKXS,CINPX,IPTHKX,ALTXT,
     *                 PZXT,TZXT
               ENDIF
            ENDIF                                                         A21060
            IF (IFRMX.EQ.1) THEN
               READ (IRD,925) (XAMNT(M,L),M=1,7),WBRODX                   A21070
               IF (IXMOL.GT.7) READ (IRD,925) (XAMNT(M,L),M=8,IXMOL)      A21080
            ELSE
               READ (IRD,927) (XAMNT(M,L),M=1,7),WBRODX
               IF (IXMOL.GT.7) READ (IRD,927) (XAMNT(M,L),M=8,IXMOL)
            ENDIF
C
C     --------------------------------------------------------------
C
C             MIXING RATIO INPUT FOR CROSS SECTIONAL MOLECULES
C
C
C     The column amount of dry air ("WDRAIR") has already been
C     calculated above, so just convert all cross sectional
C     molecules which may be in mixing ratio to molecular density
C     using WDRAIR
C
C     NOTE that if XAMNT is greater than one, then column density
C               if XAMNT is less than or equal to one, then mixing ratio
C
         DO 35 M = 1,IXMOL
            IF (XAMNT(M,L).LE.1) XAMNT(M,L) = XAMNT(M,L)*WDRAIR
 35      CONTINUE
C
C     --------------------------------------------------------------
C
 40   CONTINUE                                                            A21090
      ENDIF                                                               A21100
C                                                                         A21110
   50 WRITE (IPR,945) XID,(YID(M),M=1,2)                                  A21120
      IF (IFORM.EQ.1) THEN
         WRITE (IPR,950)
      ELSE
         WRITE (IPR,951)
      ENDIF
C                                                                         A21140
C     --------------------------------------------------------------
C
      DO 80 L = 1, NLAYRS                                                 A21150
         IF (IATM.GT.0.) SECL(L) = 1.0                                    A21160
         DO 60 M = 1, NMOL                                                A21180
            WKL(M,L) = WKL(M,L)*SECL(L)                                   A21190
   60    CONTINUE                                                         A21200
         IF (IXSECT.GE.1) THEN                                            A21210
            DO 70 M = 1, IXMOLS                                           A21220
               XAMNT(M,L) = XAMNT(M,L)*SECL(L)                            A21230
   70       CONTINUE                                                      A21240
         ENDIF                                                            A21250
         WBRODL(L) = WBRODL(L)*SECL(L)                                    A21260
         SECL(L) = 1.0                                                    A21270
   80 CONTINUE                                                            A21280
C
C     --------------------------------------------------------------
C                                                                         A21290
C     LTGNT = TOP ANTERIOR LAYER FOR TANGENT VIEWING CASE                 A21300
C     LTGNT = NLAYER FOR SPACE-TO-SPACE PATH OR IPATHL = 1 OR 3           A21310
C     OTHERWISE, LTGNT = LAYER WHICH INCLUDES H1                          A21320
C                              (I.E. OBSERVATION POINT)                   A21330
C                                                                         A21340
      LTGNT = NLAYRS                                                      A21350
      ITCNT = IPTH(1)                                                     A21360
      DO 90 ILAYR = 1, NLAYRS                                             A21370
         IF (IPTH(ILAYR).NE.ITCNT) LTGNT = ILAYR-1                        A21380
         IF (IPTH(ILAYR).EQ.2) THEN                                       A21390
            LH1 = ILAYR                                                   A21400
            LH2 = ILAYR                                                   A21410
         ENDIF                                                            A21420
         IF (IPTH(ILAYR).EQ.3) THEN                                       A21430
            LH2 = ILAYR                                                   A21440
         ELSE                                                             A21450
            LH1 = ILAYR                                                   A21460
         ENDIF                                                            A21470
         IF (IPTH(ILAYR).EQ.1) LH1 = ILAYR                                A21480
         ITCNT = IPTH(ILAYR)                                              A21490
   90 CONTINUE                                                            A21500
C                                                                         A21510
      LTNSAV = LTGNT                                                      A21520
      LH1SAV = LH1                                                        A21530
      LH2SAV = LH2                                                        A21540
C                                                                         A21550
      DV = 0.                                                             A21560
      PWTD = 0.                                                           A21570
      TWTD = 0.                                                           A21580
      WTOT = 0.                                                           A21590
      PWTX = 0.                                                           A21600
      TWTX = 0.                                                           A21610
      WTOX = 0.                                                           A21620
C                                                                         A21630
C     DVTST=DVSET                                                         A21640
C                                                                         A21650
C     Write message if IOD=2 (Optical depth flag) and IMRG = 1
C
      IF (IOD.EQ.2.AND.IMRG.EQ.1) THEN
         WRITE(IPR,953)
         IF (DVSET.NE.0) WRITE(IPR,954) DVSET
      ENDIF
C
C     LOOP OVER LAYERS                                                    A21660
C                                                                         A21670
      DO 130 L = 1, NLAYRS                                                A21680
         IPROB = 0                                                        A21690
         FACTOR = 1.                                                      A21700
         IF ((IPTH(L).EQ.2).AND.(IANT.EQ.0)) FACTOR = 2.                  A21710
         SUMWK = 0.                                                       A21720
         DO 100 M = 1, NMOL                                               A21730
            SUMWK = SUMWK+WKL(M,L)                                        A21740
            WMT(M) = WMT(M)+WKL(M,L)*FACTOR                               A21750
  100    CONTINUE                                                         A21760
         WTOTL(L) = SUMWK+WBRODL(L)                                       A21770
         SUMN2 = SUMN2+WBRODL(L)*FACTOR                                   A21780
         WTOT = WTOT+WTOTL(L)*FACTOR                                      A21790
         PWTD = PWTD+PAVEL(L)*WTOTL(L)*FACTOR                             A21800
         TWTD = TWTD+TAVEL(L)*WTOTL(L)*FACTOR                             A21810
         FRH2O = WKL(1,L)/WTOTL(L)                                        A21820
         ALFCOR = (PAVEL(L)/P0)*SQRT(TEMP0/TAVEL(L))                      A21830
C                                                                         A21840
C     CROSS SECTIONS                                                      A21850
C                                                                         A21860
         IF (IXSECT.GE.1) THEN                                            A21870
            SUMXK = 0.                                                    A21880
            DO 110 M = 1, IXMOLS                                          A21890
               SUMXK = SUMXK+XAMNT(M,L)                                   A21900
               WXT(M) = WXT(M)+XAMNT(M,L)*FACTOR                          A21910
  110       CONTINUE                                                      A21920
            WTOTX(L) = SUMXK+WBRODL(L)                                    A21930
            WTOX = WTOX+WTOTX(L)*FACTOR                                   A21940
            PWTX = PWTX+PAVEL(L)*WTOTX(L)*FACTOR                          A21950
            TWTX = TWTX+TAVEL(L)*WTOTX(L)*FACTOR                          A21960
         ENDIF                                                            A21970
C                                                                         A21980
C     CORRECT FOR WATER SELF BROADENING                                   A21990
C                                                                         A22000
         H2OSLF = (1.-FRH2O+5.*FRH2O)                                     A22010
         H2OSL(L) = H2OSLF                                                A22020
         ALBAR = ALFAL0*ALFCOR*H2OSLF                                     A22030
         ALBL(L) = ALBAR                                                  A22040
C                                                                         A22050
C     3.58115E-07 = SQRT( 2.*ALOG(2.)*AVOG*BOLTZ/(CLIGHT*CLIGHT) )        A22060
C                                                                         A22070
         ADBAR = 3.58115E-07*(0.5*(V1+V2))*SQRT(TAVEL(L)/AVMASS)          A22080
         ADBL(L) = ADBAR                                                  A22090
         AVBAR = 0.5*(ALBAR+SQRT(ALBAR*ALBAR+4.*ADBAR*ADBAR))             A22100
C                                                                         A22110
         AVBL(L) = AVBAR                                                  A22120
C                                                                         A22130
         OLDDV = DV                                                       A22140
         DV = AVBAR/SAMPLE                                                A22150
         IF (IHIRAC.EQ.2) DV = ALBAR/SAMPLE                               A22160
         IF (IHIRAC.EQ.3) DV = ADBAR/SAMPLE                               A22170
C                                                                         A22180
C     SKIP TO NEXT LAYER IF IOD=2 (OPTICAL DEPTH FLAG)
C     AND IMRG = 1
C
         IF (IOD.EQ.2.AND.IMRG.EQ.1) THEN
            DVL(L) = DV
            TYPE = 0.
            ITYPE = -99
            IF (IFORM.EQ.1) THEN
               WRITE (IPR,960) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *              TAVEL(L),ALBL(L),ADBL(L),AVBL(L),ZETA,
     *              DV,H2OSL(L),DV,TYPE,ITYPE,IPTH(L),SECL(L)
            ELSE
               WRITE (IPR,961) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *              TAVEL(L),ALBL(L),ADBL(L),AVBL(L),ZETA,
     *              DV,H2OSL(L),DV,TYPE,ITYPE,IPTH(L),SECL(L)
            ENDIF
            GOTO 130
         ENDIF
C
         IF (DV.LT.DVSET.AND.ISET.EQ.1) THEN                              A22190
            DV = OLDDV                                                    A22200
            IF (L.EQ.1.AND.DV.EQ.0.) DV = DVSET                           A22210
            IDVSET = IDVSET+1                                             A22220
            IF (IDVSET.EQ.1)                                              A22230
     *         WRITE (IPR,955) DV,ALTZ(L-1),HT1,ALTZ(L),HT2               A22240
         ENDIF                                                            A22250
C                                                                         A22260
         DVC = DV                                                         A22270
         TYPE = 0.                                                        A22280
         ITYPE = 99                                                       A22290
         IF (L.EQ.1) THEN                                                 A22300
C                                                                         A22310
C     DV IS ASSUMED TO BE .LT. 1                                          A22320
C     SET DV TO 3 SIGNIFICANT FIGURES                                     A22330
C                                                                         A22340
            ISCAL = ALOG10(DV)-3.                                         A22350
            SCAL = 10.**ISCAL                                             A22360
            IDV = (DV/SCAL)+0.5                                           A22370
C                                                                         A22380
C     SET IDV TO BE EVEN                                                  A22390
C                                                                         A22400
            IF (MOD(IDV,2).GT.0) IDV = IDV+1                              A22410
            DV = SCAL*FLOAT(IDV)                                          A22420
C                                                                         A22430
         ELSE                                                             A22440
C                                                                         A22450
            IF (JTYPE.EQ.1) GO TO 120                                     A22460
            TYPE = OLDDV/DV                                               A22470
            TYPMAX = 2.5                                                  A22480
            IF (TYPE.GT.TYPMAX) THEN                                      A22490
               IPROB = 1                                                  A22500
               ISTOP = 1                                                  A22510
            ELSEIF (TYPE.GE.1.2) THEN                                     A22520
C                                                                         A22530
C     TYPE IS BETWEEN 1.2 AND TYPMAX                                      A22540
C                                                                         A22550
               DV = OLDDV                                                 A22560
               ITYPE = 1./(TYPE-1.)+0.5                                   A22570
               IF (ITYPE.EQ.3) ITYPE = 2                                  A22580
               DV = OLDDV*FLOAT(ITYPE)/FLOAT(ITYPE+1)                     A22590
            ELSEIF (TYPE.GE.0.8) THEN                                     A22600
C                                                                         A22610
C     TYPE IS BETWEEN 0.8 AND 1.2 (SET TO 1.0)                            A22620
C                                                                         A22630
               DV = OLDDV                                                 A22640
               ITYPE = 0                                                  A22650
            ELSE                                                          A22660
C                                                                         A22670
C     TYPE IS LESS THAN 0.8                                               A22680
C                                                                         A22690
               DV = OLDDV                                                 A22700
               ITYPE = 0                                                  A22710
               IF (IEMIT.EQ.0) THEN                                       A22720
                  ITYPE = TYPE/(1.-TYPE)+0.5                              A22730
                  DV = DV*FLOAT(ITYPE+1)/FLOAT(ITYPE)                     A22740
                  ITYPE = -ITYPE                                          A22750
               ENDIF                                                      A22760
            ENDIF                                                         A22770
         ENDIF                                                            A22780
C                                                                         A22790
         DVL(L) = DV                                                      A22800
C                                                                         A22810
         ITYL(L) = ITYPE                                                  A22820
C                                                                         A22830
  120    IF (JTYPE.EQ.1) THEN                                             A22840
            IF (ITYL(L).NE.99) THEN                                       A22850
C                                                                         A22860
               IF (ITYL(L).EQ.0) THEN                                     A22870
                  DV = OLDDV                                              A22880
               ELSE                                                       A22890
                  DV = OLDDV*FLOAT(ITYL(L))/FLOAT(ITYL(L)+1)              A22900
               ENDIF                                                      A22910
               DVL(L) = DV                                                A22920
C                                                                         A22930
            ENDIF                                                         A22940
         ENDIF                                                            A22950
         ZETA = ALBAR/(ALBAR+ADBAR)                                       A22960
C                                                                         A22970
         DV = DVL(L)                                                      A22980
C                                                                         A22990
         IF (IFRMX.EQ.1) THEN
            WRITE (IPR,960) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *           TAVEL(L),                                                A23000
     *           ALBL(L),ADBL(L),AVBL(L),ZETA,DVC,H2OSL(L),               A23010
     *           DVL(L),TYPE,ITYL(L),IPTH(L),SECL(L)                      A23020
         ELSE
            WRITE (IPR,961) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *           TAVEL(L),
     *           ALBL(L),ADBL(L),AVBL(L),ZETA,DVC,H2OSL(L),
     *           DVL(L),TYPE,ITYL(L),IPTH(L),SECL(L)
         ENDIF
         IF (IPROB.GT.0) WRITE (IPR,962) TYPMAX                           A23030
  130 CONTINUE                                                            A23040
C    
C     SKIP TO END WHEN USING EXACT CALCULATED DV FOR OPTICAL
C     DEPTH CALCULATIONS (IOD = 2,IMRG = 1)
C
      IF (IOD.EQ.2.AND.IMRG.EQ.1) GOTO 145
C
      PWTD = PWTD/WTOT                                                    A23050
      TWTD = TWTD/WTOT                                                    A23060
      IF (IXSECT.GE.1) THEN                                               A23070
         PWTX = PWTX/WTOX                                                 A23080
         TWTX = TWTX/WTOX                                                 A23090
      ENDIF                                                               A23100
      IF (ISTOP.EQ.1) WRITE (IPR,965)                                     A23110
      IF (ISTOP.EQ.1) STOP 'PATH; ISTOP EQ 1'                             A23120
C
C     If DVOUT is nonzero (IOD=1,IMRG=1 -> interpolate optical depths to
C     value of DVOUT), then test to be sure that DVOUT is not more than
C     20% larger than the monochromatic DV (and thus ensuring enough
C     monochromatic points are available to reach the V2 endpoint for
C     the interpolated spectrum).
C
      IF (DVOUT.GT.0.) THEN
         RATOUT = 1.
         RATOUT = DVOUT/DV
         IF (RATOUT.GT.1.2) THEN
            WRITE (IPR,968) RATOUT,DVOUT,DV
            STOP 'PATH; DVOUT ERROR, SEE TAPE6'
         ENDIF
      ENDIF
C
C     If DVSET is nonzero (set the final layer DV to the value of
C     DVSET), then test to be sure that DVSET is not more than
C     20% different than the monochromatic DV.
C
      IF (DVSET.GT.0.) THEN                                               A23130
         RATIO = 1.                                                       A23140
         IF (DVSET.GT.0.) RATIO = DVSET/DV                                A23150
         IF (ISET.EQ.0) THEN                                              A23160
            IF (RATIO.GT.1.2.OR.RATIO.LT.0.8) THEN                        A23170
               WRITE (IPR,967) RATIO,DVSET,DV                             A23180
               STOP 'PATH; RATIO ERROR'                                   A23190
            ENDIF                                                         A23200
         ENDIF                                                            A23210
         WRITE (IPR,945) XID,(YID(M),M=1,2)                               A23220
         WRITE (IPR,948)
         IF (IFORM.EQ.1) THEN
            WRITE (IPR,950)
         ELSE
            WRITE (IPR,951)
         ENDIF
         DO 140 L = 1, NLAYRS                                             A23240
            ALBAR = ALBL(L)                                               A23250
            ADBAR = ADBL(L)                                               A23260
            AVBAR = AVBL(L)                                               A23270
C                                                                         A23280
            OLDDV = DV                                                    A23290
            DV = DVL(L)                                                   A23300
            DVC = DV                                                      A23310
            DV = DVL(L)*RATIO                                             A23320
            TYPE = 0.                                                     A23330
            IF (L.GT.1) TYPE = OLDDV/DV                                   A23340
            ITYPE = ITYL(L)                                               A23350
            DVL(L) = DV                                                   A23360
            ZETA = ALBL(L)/(ALBL(L)+ADBL(L))                              A23370
            IF (IFORM.EQ.1) THEN
               WRITE (IPR,960) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),      A23380
     *              TAVEL(L),ALBL(L),ADBL(L),AVBL(L),ZETA,                A23390
     *              DVC,H2OSL(L),DV,TYPE,ITYPE,IPTH(L),SECL(L)            A23400
            ELSE
               WRITE (IPR,961) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *              TAVEL(L),ALBL(L),ADBL(L),AVBL(L),ZETA,
     *              DVC,H2OSL(L),DV,TYPE,ITYPE,IPTH(L),SECL(L)
            ENDIF
  140    CONTINUE                                                         A23410
      ENDIF                                                               A23420
C                                                                         A23430
 145  CONTINUE
C
      IF (NLAYRS.LT.5) THEN                                               A23440
         WRITE (IPR,970)                                                  A23450
      ELSE                                                                A23460
         WRITE (IPR,945) XID,(YID(M),M=1,2)                               A23470
      ENDIF                                                               A23480
C
C     --------------------------------------------------------------
C
C     Write out column densities for molecules to TAPE6
C
C
      IF (IFORM.EQ.1) THEN                                                A23490
         WRITE (IPR,974) (HMOLID(I),I=1,7),HOLN2                          A23500
         DO 150 L = 1, NLAYRS
            WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *           TAVEL(L),                                                A23510
     *           IPTH(L),(WKL(M,L),M=1,7),WBRODL(L)                       A23520
 150     CONTINUE                                                         A23530
         IF (NLAYRS.GT.1) THEN                                            A23540
            WRITE (IPR,985)                                               A23550
            L = NLAYRS                                                    A23560
            WRITE (IPR,990) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTD,TWTD,          A23570
     *                      (WMT(M),M=1,7),SUMN2                          A23580
         ENDIF                                                            A23590
      ELSE
         WRITE (IPR,975) (HMOLID(I),I=1,7),HOLN2
         DO 151 L = 1, NLAYRS
            WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *           TAVEL(L),
     *           IPTH(L),(WKL(M,L),M=1,7),WBRODL(L)
 151     CONTINUE
         IF (NLAYRS.GT.1) THEN
            WRITE (IPR,985)
            L = NLAYRS
            WRITE (IPR,991) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTD,TWTD,
     *                      (WMT(M),M=1,7),SUMN2
         ENDIF
      ENDIF
C
      IF (NMOL.GT.7) THEN                                                 A23600
         DO 170 MLO = 8, NMOL, 8                                          A23610
            MHI = MLO+7                                                   A23620
            MHI = MIN(MHI,NMOL)                                           A23630
            IF (NLAYRS.LT.5) THEN                                         A23640
               WRITE (IPR,970)                                            A23650
            ELSE                                                          A23660
               WRITE (IPR,945) XID,(YID(M),M=1,2)                         A23670
            ENDIF                                                         A23680
            IF (IFORM.EQ.1) THEN
               WRITE (IPR,974) (HMOLID(I),I=MLO,MHI)                      A23690
               DO 160 L = 1, NLAYRS                                       A23700
                  WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),   A23710
     *                 TAVEL(L),IPTH(L),(WKL(M,L),M=MLO,MHI)              A23720
 160           CONTINUE
               IF (NLAYRS.GT.1) THEN
                  WRITE (IPR,985)
                  L = NLAYRS
                  WRITE (IPR,990) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTD,TWTD,
     *                            (WMT(M),M=MLO,MHI)
               ENDIF
            ELSE
               WRITE (IPR,975) (HMOLID(I),I=MLO,MHI)                      A23730
               DO 161 L = 1, NLAYRS                                       A23740
                  WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),   A23750
     *                 TAVEL(L),IPTH(L),(WKL(M,L),M=MLO,MHI)              A23760
 161           CONTINUE                                                   A23770
               IF (NLAYRS.GT.1) THEN                                      A23780
                  WRITE (IPR,985)
                  L = NLAYRS
                  WRITE (IPR,991) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTD,TWTD,
     *                            (WMT(M),M=MLO,MHI)
               ENDIF
            ENDIF                                                         A23790
C     
 170     CONTINUE                                                         A23800
      ENDIF                                                               A23810
C
C     --------------------------------------------------------------
C
C     Write out mixing ratios for molecules to TAPE6 in either
C     15.7 format (IFORM = 1) or 10.4 format (IFORM = 0).
C
C           Reset WDRAIR for each layer (WKL(M,L) now in column density)
C
C
      IF (IFORM.EQ.1) THEN
         WRITE (IPR,976) (HMOLID(I),I=1,7),HOLN2
         DO 172 L = 1, NLAYRS
            WDRAIR = WBRODL(L)
            DO 171 M = 2,NMOL
               WDRAIR = WDRAIR + WKL(M,L)
 171        CONTINUE
            WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *           TAVEL(L),
     *           IPTH(L),(WKL(M,L)/WDRAIR,M=1,7),WBRODL(L)
 172     CONTINUE
      ELSE
         WRITE (IPR,977) (HMOLID(I),I=1,7),HOLN2
         DO 174 L = 1, NLAYRS
            WDRAIR = WBRODL(L)
            DO 173 M = 2,NMOL
               WDRAIR = WDRAIR + WKL(M,L)
 173        CONTINUE
            WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *           TAVEL(L),
     *           IPTH(L),(WKL(M,L)/WDRAIR,M=1,7),WBRODL(L)
 174     CONTINUE
      ENDIF
C
C
      IF (NMOL.GT.7) THEN
         DO 178 MLO = 8, NMOL, 8
            MHI = MLO+7
            MHI = MIN(MHI,NMOL)
            IF (NLAYRS.LT.5) THEN
               WRITE (IPR,970)
            ELSE
               WRITE (IPR,945) XID,(YID(M),M=1,2)
            ENDIF
            IF (IFORM.EQ.1) THEN
               WRITE (IPR,976) (HMOLID(I),I=MLO,MHI)
               DO 176 L = 1, NLAYRS
                  WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *                 TAVEL(L),IPTH(L),(WKL(M,L)/WDRAIR,M=MLO,MHI)
 176           CONTINUE
            ELSE
               WRITE (IPR,977) (HMOLID(I),I=MLO,MHI)
               DO 177 L = 1, NLAYRS
                  WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *                 TAVEL(L),IPTH(L),(WKL(M,L)/WDRAIR,M=MLO,MHI)
 177           CONTINUE
            ENDIF
 178     CONTINUE
      ENDIF
C
C     --------------------------------------------------------------
C
C     Write out column densities for cross sectional molecules to TAPE6
C
C
      IF (IXSECT.GE.1) THEN                                               A23820
         DO 190 MLO = 1, IXMOLS, 8                                        A23830
            MHI = MLO+7                                                   A23840
            MHI = MIN(MHI,IXMOLS)                                         A23850
            IF (NLAYRS.LT.5.AND.MLO.NE.1) THEN                            A23860
               WRITE (IPR,970)                                            A23870
            ELSE                                                          A23880
               WRITE (IPR,995) XID,(YID(M),M=1,2)                         A23890
            ENDIF                                                         A23900
            IF (IFRMX.EQ.1) THEN
               WRITE (IPR,974) (XSNAME(I),I=MLO,MHI)                      A23910
               DO 180 L = 1, NLAYRS                                       A23920
                  WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),   A23930
     *                 TAVEL(L),IPTH(L),(XAMNT(M,L),M=MLO,MHI)            A23940
 180           CONTINUE
               IF (NLAYRS.GT.1) THEN
                  WRITE (IPR,985)
                  L = NLAYRS
                  WRITE (IPR,990) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTX,TWTX,
     *                            (WXT(M),M=MLO,MHI)
               ENDIF
            ELSE
               WRITE (IPR,975) (XSNAME(I),I=MLO,MHI)
               DO 181 L = 1, NLAYRS
                  WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *                 TAVEL(L),IPTH(L),(XAMNT(M,L),M=MLO,MHI)
 181           CONTINUE
               IF (NLAYRS.GT.1) THEN                                      A23950
                  WRITE (IPR,985)                                         A23960
                  L = NLAYRS                                              A23970
                  WRITE (IPR,991) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTX,TWTX,    A23980
     *                            (WXT(M),M=MLO,MHI)                      A23990
               ENDIF                                                      A24000
            ENDIF                                                         A24010
 190     CONTINUE                                                         A24020
C
C     --------------------------------------------------------------
C
C        Write out mixing ratios for cross sectional
C             molecules to TAPE6
C
C
         DO 198 MLO = 1, IXMOLS, 8
            MHI = MLO+7
            MHI = MIN(MHI,IXMOLS)
            IF (IFRMX.EQ.1) THEN
               WRITE (IPR,976) (XSNAME(I),I=MLO,MHI)
               DO 195 L = 1, NLAYRS
                  WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *                 TAVEL(L),IPTH(L),(XAMNT(M,L)/WDRAIR,M=MLO,MHI)
 195           CONTINUE
            ELSE
               WRITE (IPR,977) (XSNAME(I),I=MLO,MHI)
               DO 196 L = 1, NLAYRS
                  WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),
     *                 TAVEL(L),IPTH(L),(XAMNT(M,L)/WDRAIR,M=MLO,MHI)
 196           CONTINUE
            ENDIF
 198     CONTINUE
C
C     --------------------------------------------------------------
C
      ENDIF                                                               A24030
      RETURN                                                              A24040
C                                                                         A24050
  900 FORMAT (1X,I1,I3,I5,F10.2,15A4)                                     A24060
  901 FORMAT (1X,I1,I3,I5,F10.2,A20,F8.2,A4,F8.2,A5,F8.3,A7)
  902 FORMAT ('0 SECANT   =',F13.4,/'0 NLAYRS=',I4,/'0 NMOL=',I4,/'0',    A24070
     *        A20,F8.2,A4,F8.2,A5,F8.3,A7)                                A24080
  905 FORMAT (A6)                                                         A24090
  907 FORMAT ('0 SECANT   =',F13.4)                                       A24100
  910 FORMAT (E15.7,F10.4,F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))
  911 FORMAT (3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))                          A24110
  912 FORMAT ('0   ********* ITYPE(L) IS SET FROM INPUT ******** ')       A24120
  915 FORMAT (E15.7,F10.4,F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))
  916 FORMAT (3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))                          A24130
  917 FORMAT (A4)                                                         A24140
  920 FORMAT (I3)                                                         A24150
  925 FORMAT (8E15.7)
  927 FORMAT (8E10.3)                                                     A24160
  930 FORMAT (I5,5X,I5)                                                   A24170
  932 FORMAT (/,'  THE CROSS-SECTION MOLECULES SELECTED ARE: ',/,/,(5X,   A24180
     *        I5,3X,A))                                                   A24190
  935 FORMAT (/,'***** IXMOL = ',I5,' *****',/)                           A24200
  937 FORMAT (/,'***** IXMOL = ',I5,' .NE. IXMOLS = ',I5,' *****',/)      A24210
  940 FORMAT (/,'***** NLAYRS = ',I5,' .NE. NLAYXS = ',I5,' *****',/)     A24220
  942 FORMAT (/,'0 SECANTX  =',F13.4,/'0 NLAYXS=',I4,/'0 ISMOLS=',I4,/,   A24230
     *        '0',15A4)                                                   A24240
  945 FORMAT ('1'/'0',10A8,2X,2(1X,A8,1X))                                A24250
  948 FORMAT ('0   ****** DVSET is set from input ******')
  950 FORMAT ('0','LAYER',26X,'P(MB)',7X,'T(K)',4X,'ALPHL',4X,'ALPHD',    A24260
     *        4X,'ALPHV',3X,'ZETA',2X,'CALC DV',2X,'H2OSLF',5X,'DV',5X,   A24270
     *        'TYPE',' ITYPE IPATH ',3X,'SECANT'/)                        A24280
  951 FORMAT ('0','LAYER',25X,'P(MB)',3X,'T(K)',4X,'ALPHL',4X,'ALPHD',    A24260
     *        4X,'ALPHV',3X,'ZETA',2X,'CALC DV',2X,'H2OSLF',5X,'DV',5X,   A24270
     *        'TYPE',' ITYPE IPATH ',3X,'SECANT'/)                        A24280
  953 FORMAT ('0 ***** EXACT CALCULATED DV USED IN CALCULATION *****')
  954 FORMAT ('0 *** DVSET =',F10.4,' ***',/,
     *        '  *** OPTICAL DEPTHS FOR EACH LAYER WILL BE ',
     *        'CALCULATED USING THE EXACT DV, SO THE VALUE ',
     *        'OF DVSET HAS BEEN IGNORED ***')
  955 FORMAT (/,'0 **** CALC DV WAS RESET TO PREVIOUS DV',F12.6,/,        A24290
     *        '  AT ALT=  ',2(F7.3,A3),' AND ABOVE')                      A24300
  960 FORMAT ('0',I5,2(F7.3,A3),1P,E15.7,0P,F8.2,3F9.6,F6.3,
     *        F9.6,F8.4,F9.6,
     *        F7.3,2I5,F12.6)                                             A24320
  961 FORMAT ('0',I5,2(F7.3,A3),F10.4,F8.2,3F9.6,F6.3,F9.6,F8.4,F9.6,     A24310
     *        F7.3,2I5,F12.6)                                             A24320
  962 FORMAT (20X,'  DV RATIO  .GT. ',F10.2)                              A24330
  965 FORMAT (/,20X,'  TYPE GT 2.5')                                      A24340
  967 FORMAT ('  RATIO ERROR ',F10.3,'  DVSET = ',F10.4,'  DV=',F10.4)    A24350
  968 FORMAT ('  RATOUT ERROR ',F10.3,'  DVOUT = ',E10.4,'  DV=',E10.4)
  970 FORMAT (////)                                                       A24360
  974 FORMAT ('0',53X,'MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',/,32X,
     *        'P(MB)',6X,'T(K)',3X,'IPATH',5X,5(A10,4X),/,60X,3(A10,4X))
  975 FORMAT ('0',53X,'MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',/,29X,    A24370
     *        'P(MB)',6X,'T(K)',3X,'IPATH',1X,8(1X,A6,3X))
  976 FORMAT (/,'1',54X,'----------------------------------',
     *         /,'0',60X,'MIXING RATIOS BY LAYER ',/,32X,
     *        'P(MB)',6X,'T(K)',3X,'IPATH',5X,8(A10,4X))                  A24380
  977 FORMAT (/,'1',54X,'----------------------------------',
     *         /,'0',60X,'MIXING RATIOS BY LAYER ',/,29X,
     *        'P(MB)',6X,'T(K)',3X,'IPATH',1X,8(1X,A6,3X))
  980 FORMAT ('0',I3,2(F7.3,A3),F15.7,F9.2,I5,2X,1P,5E15.7,0P,/,
     *         54X,1P,3E15.7,0P)
  982 FORMAT ('0',I3,2(F7.3,A3),F12.5,F9.2,I5,2X,1P,8E10.3,0P)            A24390
  985 FORMAT ('0',54X,'ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH')     A24400
  990 FORMAT ('0',I3,2(F7.3,A3),F15.7,F9.2,7X,1P,5E15.7,0P,/,
     *         54X,1P,3E15.7,0P)
  991 FORMAT ('0',I3,2(F7.3,A3),F12.5,F9.2,7X,1P,8E10.3,0P)               A24410
  995 FORMAT ('1'/'0',10A8,2X,2(1X,A8,1X),/,/,'0',53X,                    A24420
     *        '     *****  CROSS SECTIONS  *****      ')                  A24430
C                                                                         A24440
      END                                                                 A24450
      BLOCK DATA BOPDPT                                                   A07600
      COMMON /CONVF/ CHI(251),RDVCHI,RECPI,ZSQBND,A3,B3,JCNVF4            A24720
C                                                                         A24470
      DATA JCNVF4 / 0 /                                                   A24800
C                                                                         A24470
      END                                                                 A24450
      SUBROUTINE OPDPTH (MPTS)                                            A24460
C                                                                         A24470
      IMPLICIT DOUBLE PRECISION (V)                                     ! A24480
C                                                                         A24490
C     OPDPTH CALLS CONTNM,LINF4,HIRAC1,NONLTE                             A24500
C                                                                         A24510
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   A24520
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       A24530
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       A24540
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      A24550
     *              EXTID(10)                                             A24560
C                                                                         A24570
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID                      & A24580
C                                                                         A24590
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       A24600
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   A24610
     *                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    A24620
      COMMON /LASIV/ VLAS,ILAS                                            A24630
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(2030)                A24640
      COMMON /SCATTR/ V1SC,V2SC,DVSC,NPTSC,SCTTR(2025)                    A24650
C                                                                         A24660
      COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)        A24670
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         A24680
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        A24690
     *              NLTEFL,LNFIL4,LNGTH4                                  A24700
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           A24710
      COMMON /CONVF/ CHI(251),RDVCHI,RECPI,ZSQBND,A3,B3,JCNVF4            A24720
C                                                                         A24730
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),                A24740
     *            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),                A24750
     *            (FSCDID(5),IEMIT) , (FSCDID(7),IPLOT),                  A24760
     *            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                  A24770
     *            (FSCDID(11),IMRG)                                       A24780
C                                                                         A24790
C     DATA JCNVF4 / 0 /                                                   A24800
C                                                                         A24810
      CALL CPUTIM (TIME0)                                                 A24820
C                                                                         A24830
      ICNTNM = MOD(IXSCNT,10)                                             A24840
      IXSECT = IXSCNT/10                                                  A24850
C                                                                         A24860
      IEMST = IEMIT                                                       A24870
      IEMIT = 0                                                           A24880
      IPFLAG = 0                                                          A24890
      DPTMST = DPTMIN                                                     A24900
      IF (IEMST.EQ.0.AND.IPATHL.EQ.2) DPTMIN = 2.*DPTMST                  A24910
C                                                                         A24920
C     PRINT LAYER INFORMATION                                             A24930
C                                                                         A24940
      IF (NOPR.EQ.0) THEN                                                 A24950
         IF (IMRG.LE.10) WRITE (IPR,900)                                  A24960
         WRITE (IPR,905) LAYRS                                            A24970
         IF (ILAS.GT.0) WRITE (IPR,910) VLAS,V1,V2                        A24980
         WRITE (IPR,915) XID,(YID(M),M=1,2),TIME0                         A24990
      ENDIF                                                               A25000
C                                                                         A25010
      XKT = TAVE/RADCN2                                                   A25020
C                                                                         A25030
C     JRAD= -1  NO RADIATION TERM IN ABSORPTION COEFFICIENTS              A25040
C     JRAD=  0  RADIATION TERM PUT IN BY PANEL                            A25050
C     JRAD=  1  RADIATION TERM INCLUDED IN LINE STRENGTHS                 A25060
C                                                                         A25070
      IF (((V1/XKT).LT.2.5).AND.(JRAD.NE.-1)) JRAD = 0                    A25080
C                                                                         A25090
C     DVABS IS USED AS A FLAG IN SUBSEQUENT PROGRAMS                      A25100
C                                                                         A25110
      DVABS = 0.                                                          A25120
      IF (ICNTNM.NE.0) THEN                                               A25130
         DVABS = 1.                                                       A25140
         V1ABS = INT(V1)                                                  A25150
         IF (V1.LT.0.) V1ABS = V1ABS-1.                                   A25160
         V1ABS = V1ABS-3.*DVABS                                           A25170
         V2ABS = INT(V2+3.*DVABS+0.5)                                     A25180
         NPTABS = (V2ABS-V1ABS)/DVABS+1.5                                 A25190
         IF (PAVE.LE.0.5) IPFLAG = 1                                      A25200
         DO 10 I = 1, 2030                                                A25210
            ABSRB(I) = 0.                                                 A25220
   10    CONTINUE                                                         A25230
         CALL CONTNM (JRAD)                                               A25240
      ENDIF                                                               A25250
      DVR4 = 0.                                                           A25260
C                                                                         A25270
      IF (ILBLF4.GE.1) THEN                                               A25280
         ALFAV = SAMPLE*DV                                                A25290
         ALFAV4 = 64.*ALFAV                                               A25300
         DVR4 = ALFAV4/SAMPLE                                             A25310
         BOUND4 = 25.                                                     A25320
         IF (ILBLF4.EQ.2.AND.IPFLAG.EQ.1) BOUND4 = 5.                     A25330
         IPTS4 = BOUND4/DVR4                                              A25340
C                                                                         A25350
         IF (NOPR.EQ.0) WRITE (IPR,920) IPTS4,DVR4,BOUND4                 A25360
C                                                                         A25370
         REWIND LINFIL                                                    A25380
         REWIND LNFIL4                                                    A25390
         V1R4 = V1-2.*DVR4                                                A25400
         V2R4 = V2+2.*DVR4                                                A25410
         V1L4 = V1R4-BOUND4-DVR4                                          A25420
         V2L4 = V2R4+BOUND4+2*DVR4                                        A25430
         IF (IHIRAC.EQ.1) CALL LINF4 (V1L4,V2L4)                          A25440
      ENDIF                                                               A25450
C                                                                         A25460
      IF (IHIRAC.EQ.1) CALL HIRAC1 (MPTS)                                 A25470
      IF (IHIRAC.EQ.4) CALL NONLTE (MPTS)                                 A25480
      IF (IHIRAC.EQ.9) CALL HIRAC1 (MPTS)                                 A25490
C                                                                         A25500
      IEMIT = IEMST                                                       A25510
      DPTMIN = DPTMST                                                     A25520
      CALL CPUTIM (TIME1)                                                 A25530
      TIMEO = TIME1-TIME0                                                 A25540
      WRITE (IPR,925) TIME1,TIMEO                                         A25550
C                                                                         A25560
      RETURN                                                              A25570
C                                                                         A25580
  900 FORMAT ('1')                                                        A25590
  905 FORMAT ('0 LAYER = ',I8)                                            A25600
  910 FORMAT ('0 VLAS  ',F20.8,8X,'V1 RESET ',F12.5,8X,'V2 RESET ',       A25610
     *        F12.5)                                                      A25620
  915 FORMAT ('0',10A8,2X,2(1X,A8,1X),/,'0 TIME ENTERING OPDPTH ',        A25630
     *        F15.3)                                                      A25640
  920 FORMAT ('0  IPTS4 FOR LINF4 = ',I5,3X,' DV FOR LINF4 = ',F10.5,     A25650
     *        5X,'BOUND FOR LINF4 =',F10.4)                               A25660
  925 FORMAT ('0 TIME LEAVING OPDPTH ',F15.3,'  TOTAL FOR LAYER ',        A25670
     *        F15.3)                                                      A25680
C                                                                         A25690
      END                                                                 A25700
