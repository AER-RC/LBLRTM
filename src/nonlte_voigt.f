C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

      SUBROUTINE NONLTE(MPTS)                                             600000
      IMPLICIT REAL*8           (V)                                     ! B00030
C
C**********************************************************************
C*                                                                     
C*                                                                     
C*    CALCULATES MONOCHROMATIC ABSORPTION COEFFICIENT FOR SINGLE LAYER 
C*    Under the conditions of NLTE                                     
C*                                                                     
C**********************************************************************
C                                                                      
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      
C                  IMPLEMENTATION:    W.O. GALLERY                     
C                                                                      
C             ALGORITHM REVISIONS:    M.W.SHEPHARD                     
C                                                                      
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.      
C                     131 Hartwell Ave,  Lexington,  MA   02421        
C                                                                      
C----------------------------------------------------------------------
C                                                                      
C               WORK SUPPORTED BY:    JPL, TES                         
C                                                                      
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL FASCOD         
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   B00580
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       B00590
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       B00600
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      B00610
     *              EXTID(10)                                             B00620
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       B00660
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   B00670
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    B00680
      COMMON /VBNLTE/ RATH2O(8),RATCO2(26),RATO3(18),RATCO(3),RATNO(3),   600110
     &               NUMH2O,NUMCO2,NUMO3,NUMCO,NUMNO                      600120
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         B00840
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        B00850
     *              NLTEFL,LNFIL4,LNGTH4                                  B00860
      COMMON /R4SUB/ VLOF4,VHIF4,ILOF4,IST,IHIF4,LIMIN4,LIMOUT,ILAST,     B00820
     *               DPTMN4,DPTFC4,ILIN4,ILIN4T                           B00830
      COMMON /HVERSN/  HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,
     *                HVROPR,HVRPST,HVRPLT,HVRTST,HVRUTL,HVRXMR,hvnlte
c
      CHARACTER*15 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *                HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR,hvnlte
      CHARACTER*15 HVRSOL
      COMMON /ARMCM1/ HVRSOL


      EQUIVALENCE (FSCDID(1),IHIRAC),(FSCDID(2),ILBLF4),
     &  (FSCDID(3),IXSCNT),(FSCDID(4),IAERSL),(FSCDID(5),IEMIT),
     &  (FSCDID(7),IPLOT), (FSCDID(8),IPATHL),(FSCDID(9),JRAD),
     &  (FSCDID(11),IMRG)
C
      CHARACTER HH2O*3,HCO2*5,HO3*3,HCO*1,HNO*1                           600144
      DIMENSION HH2O( 8),EH2O( 8),NDGH2O( 8),                             600150
     &          HCO2(26),ECO2(26),NDGCO2(26),                             600160
     &          HO3 (18),EO3 (18),NDGO3 (18),                             600170
     &          HCO ( 3),ECO ( 3),NDGCO ( 3),                             600180
     &          HNO ( 3),ENO ( 3),NDGNO ( 3)                              600190
      CHARACTER*80 TIT(3)                                                 600200
      CHARACTER*10 TIT1,TIT2,TIT3,TIT4,TIT5                               600210
      DATA TIT1/' H20      '/, TIT2/' CO2      '/,TIT3/' O3       '/,     600220
     &     TIT4/' CO       '/, TIT5/' NO       '/                         600230
C 
      DATA NDGH2O/8*1/,NDGO3/18*1/,NDGCO/3*1/,NDGNO/3*1/                  600250
C                                                                         600270
      DATA  (HH2O(I),EH2O(I),I=1,8)/                                      600280
     1      '000' ,     0.   ,                                            600290
     2      '010' ,  1594.750,                                            600300
     3      '020' ,  3151.630,                                            600310
     4      '100' ,  3657.053,                                            600320
     5      '001' ,  3755.930,                                            600330
     6      '030' ,  4666.793,                                            600340
     7      '110' ,  5234.977,                                            600350
     8      '011' ,  5333.269/                                            600360
C                                                                         600370
      DATA  (HCO2(I),ECO2(I),NDGCO2(I),I=1, 9)/                           600380
     1      '00001' ,    0.  , 1 ,                                        600390
     2      '01101' ,  667.380,2 ,                                        600400
     3      '10002' , 1285.409,1 ,                                        600410
     4      '02201' , 1335.132,2 ,                                        600420
     5      '10001' , 1388.185,1 ,                                        600430
     6      '11102' , 1932.470,2 ,                                        600440
     7      '03301' , 2003.246,2 ,                                        600450
     8      '11101' , 2076.856,2 ,                                        600460
     9      '00011' , 2349.143,1 /                                        600470
      DATA  (HCO2(I),ECO2(I),NDGCO2(I),I=10,26)/                          600380
     Z      '20003' , 2548.366,1 ,                                        600470
     Z      '12202' , 2585.022,2 ,                                        600470
     Z      '20002' , 2671.143,1 ,                                        600470
     Z      '04401' , 2671.717,2 ,                                        600470
     Z      '12201' , 2760.725,2 ,                                        600470
     Z      '20001' , 2797.135,1 ,                                        600470
     A      '01111' , 3004.012,2 ,                                        600480
     1      '10012' , 3612.842,1 ,                                        600490
     2      '02211' , 3659.273,2 ,                                        600500
     3      '10011' , 3714.783,1 ,                                        600510
     Z      '11112' , 4247.706,2 ,                                        600470
     Z      '03311' , 4314.914,2 ,                                        600470
     Z      '11111' , 4390.629,2 ,                                        600470
     4      '20013' , 4853.623,1 ,                                        600520
     Z      '04411' , 4970.931,1 ,                                        600470
     5      '20012' , 4977.834,1 ,                                        600530
     6      '20011' , 5099.660,1 /                                        600540
C                                                                         600550
      DATA (HO3(I),EO3(I),I=1,18)/                                        600560
     1     '000' ,    0.   ,                                              600570
     2     '010' ,  700.931,                                              600580
     3     '001' , 1042.084,                                              600590
     4     '100' , 1103.140,                                              600600
     5     '020' , 1399.275,                                              600610
     6     '011' , 1726.528,                                              600620
     7     '110' , 1796.261,                                              600630
     8     '002' , 2057.892,                                              600640
     9     '101' , 2110.785,                                              600650
     A     '200' , 2201.157,                                              600660
     1     '111' , 2785.245,                                              600670
     2     '003' , 3041.200,                                              600680
     3     '004' , 3988.  ,                                               600690
     4     '005' , 4910.  ,                                               600700
     5     '006' , 5803.  ,                                               600710
     6     '007' , 6665.  ,                                               600720
     7     '008' , 7497.  ,                                               600730
     8     '009' , 8299.  /                                               600740
C                                                                         600750
      DATA (HCO(I),ECO(I),I=1,3)/                                         600760
     1     '0' ,    0.  ,                                                 600770
     2     '1' , 2143.272,                                                600780
     3     '2' , 4260.063/                                                600790
C                                                                         600800
       DATA (HNO(I),ENO(I),I=1,3)/                                        600810
     1      '0' ,    0.  ,                                                600820
     2      '1' , 1878.077,                                               600830
     3      '2' , 3724.067/                                               600840
C
C     ASSIGN SCCS VERSION NUMBER TO MODULE 
C
       hvnlte = '$Revision$'

      NUMH2O =  8
      NUMCO2 = 26
      NUMO3  = 18
      NUMCO  =  3
      NUMNO  =  3
C                                                                    
      REWIND NLTEFL                                                       600860

C NOTE THAT TIT IS A 3-ELEMENT ARRAY OF CHARACTER*80
      READ(NLTEFL,900) TIT                                                600870
  900 FORMAT(A80)                                                         600880
      READ(NLTEFL,910) IVIB,MOLNEQ                                        600890
  910 FORMAT(2I5)                                                         600900
      WRITE(IPR,920) TIT(1),TIT(2),TIT(3)
  920 FORMAT(/2X,'TAPE4 HEADER:',/2X,A80,/2X,A80,/2X,A80)
      WRITE(IPR,921) IVIB,ALTAV,TAVE
  921 FORMAT(/' IVIB =',I5,/'  ALT = ',F10.3,'  TEMP =',F10.3)            600920

      IF(IVIB.EQ.1) THEN
          CALL VIBTMP(XKT,ALTAV,NLTEFL,NUMH2O,HH2O,
     &        NDGH2O(1),EH2O(1),RATH2O(1),TIT1)
          CALL VIBTMP(XKT,ALTAV,NLTEFL,NUMCO2,HCO2,
     &        NDGCO2(1),ECO2(1),RATCO2(1),TIT2)
          CALL VIBTMP(XKT,ALTAV,NLTEFL,NUMO3 ,HO3,
     &        NDGO3(1),EO3(1),RATO3(1),TIT3)
          CALL VIBTMP(XKT,ALTAV,NLTEFL,NUMCO ,HCO,
     &        NDGCO(1),ECO(1),RATCO(1),TIT4)
          CALL VIBTMP(XKT,ALTAV,NLTEFL,NUMNO ,HNO,
     &        NDGNO(1),ENO(1),RATNO(1),TIT5)
      ELSE                                                                600990
          CALL VIBPOP(XKT,ALTAV,NLTEFL,NUMH2O,HH2O,
     &        NDGH2O(1),EH2O(1),RATH2O(1),TIT1)
          CALL VIBPOP(XKT,ALTAV,NLTEFL,NUMCO2,HCO2,
     &        NDGCO2(1),ECO2(1),RATCO2(1),TIT2)
          CALL VIBPOP(XKT,ALTAV,NLTEFL,NUMO3 ,HO3,
     &        NDGO3(1),EO3(1),RATO3(1),TIT3)
          CALL VIBPOP(XKT,ALTAV,NLTEFL,NUMCO ,HCO,
     &        NDGCO(1),ECO(1),RATCO(1),TIT4)
          CALL VIBPOP(XKT,ALTAV,NLTEFL,NUMNO ,HNO,
     &        NDGNO(1),ENO(1),RATNO(1),TIT5)
      ENDIF                                                               601050
      IPFLAG = 0
      IF(PAVE .LE. 0.5) IPFLAG = 1
      ALFAV=SAMPLE*DV
      ALFAV4= 64.*ALFAV
      DVR4=ALFAV4/SAMPLE
      BOUND4=25.
      IF (ILBLF4.EQ.2 .AND. IPFLAG.EQ.1) BOUND4=5.
      V1R4=V1-2.*DVR4
      V2R4=V2+2.*DVR4
      VLOF4=V1R4-BOUND4-DVR4
      VHIF4=V2R4+BOUND4+2*DVR4

C*****Note: this check not needed for lblrtm
C CHECK TO SEE IF VLOF4 OR VHIF4 FALLS IN LINE COUPLING REGION
C      CALL CHKLNC(VLOF4,VHIF4,(SAMPLE*DV),0)

      IF(ILBLF4.GE.1) CALL LINF4Q(VLOF4,VHIF4)
      CALL HIRACQ(MPTS)                                                   601060
      WRITE(IPR,930)                                                      601070
  930 FORMAT(//)                                                          601080
      RETURN                                                              601090
      END                                                                 601100

c ----------------------------------------------------------------
      SUBROUTINE LINF4Q (V1L4,V2L4)                                        D00010
C                                                                         D00020
      IMPLICIT REAL*8           (V)                                     ! D00030
C                                                                         D00040
C     SUBROUTINE LINF4 READS THE LINES AND SHRINKS THE LINES FOR LBLF4    D00050
C                                                                         D00060
      PARAMETER (NTMOL=36,NSPECI=85)                                      D00070
C                                                                         D00080
      COMMON /ISVECT/ ISOVEC(NTMOL),ISO82(NSPECI),ISONM(NTMOL),           D00090
     *                SMASSI(NSPECI)                                      D00100
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           D00110
C                                                                         D00120
      REAL*8            HID,HMOLIL,HID1,HLINHD,VMNCPL,VMXCPL,EXTSPC     & D00130
C                                                                         D00140
      COMMON /BUFID/ HID(10),HMOLIL(64),MOLCNT(64),MCNTLC(64),            D00150
     *               MCNTNL(64),SUMSTR(64),NMOI,FLINLO,FLINHI,            D00160
     *               ILIN,ILINLC,ILINNL,IREC,IRECTL,HID1(2),LSTWDL        D00170
C                                                                         D00180
      COMMON VNU(1250),SP(1250),ALFA0(1250),EPP(1250),MOL(1250),          D00190
     *       SPP(1250),SRAD(1250)                                         D00200
C                                                                         D00210
      COMMON /IOU/ IOUT(250)                                              D00220
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   D00230
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       D00240
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       D00250
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      D00260
     *              EXTID(10)                                             D00270
C                                                                         D00280
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SEC   ,       XALTZ
C                                                                         D00300
      COMMON /FILHDR/ XID(10),SEC   ,PAVE,TAVE,HMOLID(60),XALTZ(4),       D00310
     *                W(60),PZL,PZU,TZL,TZU,WBROAD,DVO,V1 ,V2 ,TBOUND,    D00320
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    D00330
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /R4SUB/ VLO,VHI,ILO,IST,IHI,LIMIN,LIMOUT,ILAST,DPTMN,        D00350
     *               DPTFC,ILIN4,ILIN4T                                   D00360
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         D00370
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        D00380
     *              NLTEFL,LNFIL4,LNGTH4                                  D00390
      COMMON /TPANEL/ VNULO,VNUHI,JLIN,NLNGT4                             D00400
      COMMON /BUFR/ VNUB(250),SB(250),ALB(250),EPPB(250),MOLB(250),       D00410
     *              HWHMB(250),TMPALB(250),PSHIFB(250),IFLG(250)          D00420
      COMMON /NGT4/ VD,SD,AD,EPD,MOLD,SPPD,ILS2D                          D00430
      COMMON /L4TIMG/ L4TIM,L4TMR,L4TMS,L4NLN,L4NLS,LOTHER
      COMMON /VBNLTE/ RATH2O(8),RATCO2(26),RATO3(18),RATCO(3),RATNO(3),   600110
     &               NUMH2O,NUMCO2,NUMO3,NUMCO,NUMNO                      600120
C                                                                         D00440
      REAL L4TIM,L4TMR,L4TMS,LOTHER
      DIMENSION MEFDP(64)                                                 D00450
      DIMENSION SCOR(NSPECI),RHOSLF(NSPECI),ALFD1(NSPECI)                 D00460
      DIMENSION ALFAL(1250),ALFAD(1250),A(4),B(4),TEMPLC(4)
      DIMENSION RCDHDR(2),IWD(2),IWD3(2),HLINHD(2),AMOLB(250)             D00480
C                                                                         D00490
      DIMENSION SABS(2)
      EQUIVALENCE (ALFA0(1),ALFAL(1)) , (EPP(1),ALFAD(1))                 D00500
      EQUIVALENCE (IHIRAC,FSCDID(1)) , (ILBLF4,FSCDID(2))                 D00510
      EQUIVALENCE (VNULO,RCDHDR(1)) , (IWD3(1),VD),                       D00520
     *            (HLINHD(1),HID(1),IWD(1)) , (MOLB(1),AMOLB(1)),          D00530
     &    (SP(1),SABS(1))
C                                                                         D00540
      DATA MEFDP / 64*0 /                                                 D00550
C                                                                         D00560
C     TEMPERATURES FOR LINE COUPLING COEFFICIENTS                         D00570
C                                                                         D00580
      DATA TEMPLC / 200.0,250.0,296.0,340.0 /                             D00590
C                                                                         D00600
C     Initialize timing for the group "OTHER" in the TAPE6 output
C
      LOTHER = 0.0
      TSHRNK = 0.0
      TBUFFR = 0.0
      TMOLN4 = 0.0
C
      CALL CPUTIM (TIMEL0)                                                D00610
C                                                                         D00620
      ILS2D = -654321
      NLNGT4 = NWDL(IWD3,ILS2D)*1250                                      D00630
      LNGTH4 = NLNGT4                                                     D00640
      PAVP0 = PAVE/P0                                                     D00650
      PAVP2 = PAVP0*PAVP0                                                 D00660
      DPTMN = DPTMIN/RADFN(V2,TAVE/RADCN2)                                D00670
      DPTFC = DPTFAC                                                      D00680
      LIMIN = 1000                                                        D00690
      CALL CPUTIM(TPAT0)
      CALL MOLEC (1,SCOR,RHOSLF,ALFD1)                                    D00700
      CALL CPUTIM(TPAT1)
      TMOLN4 = TMOLN4 + TPAT1-TPAT0
C                                                                         D00710
      TIMR = 0.                                                           D00720
      TIMS = 0.                                                           D00730
      SUMS = 0.                                                           D00740
      ILAST = 0                                                           D00750
      ILINLO = 0                                                          D00760
      ILINHI = 0                                                          D00770
      ILO = 1                                                             D00780
      IST = 1                                                             D00790
      NLINS = 0                                                           D00800
      NLIN = 0                                                            D00810
C                                                                         D00820
      VLO = V1L4                                                          D00830
      VHI = V2L4                                                          D00840
C                                                                         D00850
      CALL CPUTIM(TPAT0)
c
      call lnfilhd_4(linfil,lnfil4,v1,v2)
c
      CALL CPUTIM(TPAT1)
      TBUFFR = TBUFFR + TPAT1-TPAT0
C                                                                         D01000
C       TEMPERATURE CORRECTION TO INTENSITY                               D01010
C       TEMPERATURE AND PRESSURE CORRECTION TO HALF-WIDTH                 D01020
C                                                                         D01030
      TRATIO = TAVE/TEMP0                                                 D01040
      RHORAT = (PAVE/P0)*(TEMP0/TAVE)                                     D01050
C                                                                         D01060
      BETA = RADCN2/TAVE                                                  D01070
      BETA0 = RADCN2/TEMP0                                                D01080
      BETACR = BETA-BETA0                                                 D01090
      DELTMP = ABS(TAVE-TEMP0)                                            D01100
      CALL CPUTIM(TPAT0)
      CALL MOLEC (2,SCOR,RHOSLF,ALFD1)                                    D01110
      CALL CPUTIM(TPAT1)
      TMOLN4 = TMOLN4 + TPAT1-TPAT0
C                                                                         D01120
C     FIND CORRECT TEMPERATURE AND INTERPOLATE FOR Y AND G                D01130
C                                                                         D01140
      DO 10 ILC = 1, 4                                                    D01150
         IF (TAVE.LE.TEMPLC(ILC)) GO TO 20                                D01160
   10 CONTINUE                                                            D01170
   20 IF (ILC.EQ.1) ILC = 2                                               D01180
      IF (ILC.EQ.5) ILC = 4                                               D01190
      RECTLC = 1.0/(TEMPLC(ILC)-TEMPLC(ILC-1))                            D01200
      TMPDIF = TAVE-TEMPLC(ILC)                                           D01210
C                                                                         D01220
      IJ = 0                                                              D01230
   30 CALL CPUTIM (TIM0)                                                  D01240
      CALL RDLNFL (IEOF,ILINLO,ILINHI)                                    D01250
      CALL CPUTIM (TIM1)                                                  D01260
      TIMR = TIMR+TIM1-TIM0                                               D01270
C                                                                         D01280
      IF (IEOF.GE.1) GO TO 60                                             D01290
C                                                                         D01300
      DO 50 J = ILINLO, ILINHI                                            D01310
         YI = 0.                                                          D01320
         GI = 0.                                                          D01330
         GAMMA1 = 0.                                                      D01340
         GAMMA2 = 0.                                                      D01350
         I = IOUT(J)                                                      D01360
         IFLAG = IFLG(I)                                                  D01370
         IF (I.LE.0) GO TO 50                                             D01380
C                                                                         D01390
         MFULL = MOLB(I)
         M = MOD(MOLB(I),100)                                             D01400
C                                                                         D01410
C     ISO=(MOD(MOLB(I),1000)-M)/100   IS PROGRAMMED AS:                   D01420
C                                                                         D01430
         ISO = MOD(MOLB(I),1000)/100                                      D01440
         ILOC = ISOVEC(M)+ISO                                             D01450
         IF ((M.GT.NMOL).OR.(M.LT.1)) GO TO 50                            D01460
         SUI = SB(I)*W(M)                                                 D01470
         IF (SUI.EQ.0.) GO TO 50                                          D01480
         IF (VNUB(I).LT.VLO) GO TO 50                                     D01490
         IJ = IJ+1                                                        D01500
C                                                                         D01510
C     Y'S AND G'S ARE STORED IN I+1 POSTION OF VNU,S,ALFA0,EPP...         D01520
C      A(1-4),  B(1-4) CORRESPOND TO TEMPERATURES TEMPLC(1-4) ABOVE       D01530
C                                                                         D01540
         IF (IFLAG.EQ.1.OR.IFLAG.EQ.3) THEN                               D01550
            A(1) = VNUB(I+1)                                              D01560
            B(1) = SB(I+1)                                                D01570
            A(2) = ALB(I+1)                                               D01580
            B(2) = EPPB(I+1)                                              D01590
            A(3) = AMOLB(I+1)                                             D01600
            B(3) = HWHMB(I+1)                                             D01610
            A(4) = TMPALB(I+1)                                            D01620
            B(4) = PSHIFB(I+1)                                            D01630
C                                                                         D01640
C     CALCULATE SLOPE AND EVALUATE                                        D01650
C                                                                         D01660
            SLOPEY = (A(ILC)-A(ILC-1))*RECTLC                             D01670
            SLOPEG = (B(ILC)-B(ILC-1))*RECTLC                             D01680
            IF (IFLAG.EQ.1) THEN                                          D01690
               YI = A(ILC)+SLOPEY*TMPDIF                                  D01700
               GI = B(ILC)+SLOPEG*TMPDIF                                  D01710
            ELSE                                                          D01720
               GAMMA1 = A(ILC)+SLOPEY*TMPDIF                              D01730
               GAMMA2 = B(ILC)+SLOPEG*TMPDIF                              D01740
            ENDIF                                                         D01750
         ENDIF                                                            D01760
C                                                                         D01770
C     IFLAG = 2 IS RESERVED FOR LINE COUPLING COEFFICIENTS ASSOCIATED     D01780
C               WITH AN EXACT TREATMENT (NUMERICAL DIAGONALIZATION)       D01790
C                                                                         D01800
C     IFLAG = 3 TREATS LINE COUPLING IN TERMS OF REDUCED WIDTHS           D01810
C                                                                         D01820
         
         VNU(IJ) = VNUB(I)+RHORAT*PSHIFB(I)                               D01830
         ALFA0(IJ) = ALB(I)                                               D01840
         EPP(IJ) = EPPB(I)                                                D01850
         MOL(IJ) = M                                                      D01860
C                                                                         D01870
         IF (VNU(IJ).EQ.0.) SUI = 2.*SUI                                  D01880
C                                                                         D01890
C     TREAT TRANSITIONS WITH UNKNOWN EPP AS SPECIAL CASE                  D01900
C                                                                         D01910
         IF (EPP(IJ).GE.0.) GO TO 40                                      D01920
         IF (DELTMP.LE.10.) EPP(IJ) = 0.                                  D01930
         IF (DELTMP.GT.10.) MEFDP(M) = MEFDP(M)+1                         D01940
         IF (DELTMP.GT.10.) SUI = 0.                                      D01950
   40    SUI = SUI*SCOR(ILOC)*EXP(-EPP(IJ)*BETACR)*                       D01960
     *         (1.+EXP(-VNU(IJ)*BETA))                                    D01970
C                                                                         D01980
         SUMS = SUMS+SUI                                                  D01990
C                                                                         D02000
C     TEMPERATURE CORRECTION OF THE HALFWIDTH                             D02010
C     SELF TEMP DEPENDENCE TAKEN THE SAME AS FOREIGN                      D02020
C                                                                         D02030
         TMPCOR = TRATIO**TMPALB(I)                                       D02040
         ALFA0I = ALFA0(IJ)*TMPCOR                                        D02050
         HWHMSI = HWHMB(I)*TMPCOR                                         D02060
         ALFAL(IJ) = ALFA0I*(RHORAT-RHOSLF(ILOC))+HWHMSI*RHOSLF(ILOC)     D02070
C                                                                         D02080
         IF (IFLAG.EQ.3)                                                  D02090
     *        ALFAL(IJ) = ALFAL(IJ)*(1.0-GAMMA1*PAVP0-GAMMA2*PAVP2)       D02100
C                                                                         D02110
         ALFAD(IJ) = VNU(IJ)*ALFD1(ILOC)                                  D02120
         NLIN = NLIN+1                                                    D02130
         SP(IJ) = SUI*(1.+GI*PAVP2)                                       D02140
         SPP(IJ) = SUI*YI*PAVP0                                           D02150
         SRAD(IJ) = 0.0

C  For NLTE lines:
         FREQ=VNU(IJ)
         RLOW=1.0
         RUPP=1.0

         IF (MFULL.GE.1000) THEN
             NLOW=MOD(MFULL/1000,100)
             NUPP=MFULL/100000
c             DELTA=EXP(-FREQ/XKT)
c xkt=tave/radcn2=1/beta 
             DELTA=EXP(-FREQ*BETA)
C 
C     PICK OUT MOLECULAR TYPE
C
C     H2O LINE 
C
             IF (M.EQ.1) THEN
                 IF (NLOW.GT.NUMH2O) STOP 'NLOW GT NUMH2O IN LNCO2Q'
                 IF (NLOW.GT.0) RLOW=RATH2O(NLOW)
                 IF (NUPP.GT.NUMH2O) STOP 'NUPP GT NUMH2O IN LNCO2Q'
                 IF (NUPP.GT.0) RUPP=RATH2O(NUPP)
C
C     CO2 LINE
C
             ELSE IF (M.EQ.2) THEN
                 IF (NLOW.GT.NUMCO2) STOP 'NLOW GT NUMCO2 IN LNCO2Q'
                 IF (NLOW.GT.0) RLOW=RATCO2(NLOW)
                 IF (NUPP.GT.NUMCO2) STOP 'NUPP GT NUMCO2 IN LNCO2Q'
                 IF (NUPP.GT.0) RUPP=RATCO2(NUPP)
C
C     O3 LINE
C
             ELSE IF (M.EQ.3) THEN
                 IF (NLOW.GT.NUMO3) STOP 'NLOW GT NUMO3 IN LNCO2Q' 
                 IF (NLOW.GT.0) RLOW=RATO3(NLOW) 
                 IF (NUPP.GT.NUMO3) STOP 'NUPP GT NUMO3 IN LNCO2Q' 
                 IF (NUPP.GT.0) RUPP=RATO3(NUPP)
C
C     CO LINE
C
             ELSE IF (M.EQ.5) THEN 
                 IF (NLOW.GT.NUMCO ) STOP 'NLOW GT NUMCO  IN LNCO2Q'
                 IF (NLOW.GT.0) RLOW=RATCO (NLOW) 
                 IF (NUPP.GT.NUMCO ) STOP 'NUPP GT NUMCO  IN LNCO2Q' 
                 IF (NUPP.GT.0) RUPP=RATCO (NUPP) 
C
C     NO LINE
C
             ELSE IF (M.EQ.8) THEN
                 IF (NLOW.GT.NUMNO) STOP 'NLOW GT NUMNO IN LNCO2Q' 
                 IF (NLOW.GT.0) RLOW=RATNO(NLOW)
                 IF (NUPP.GT.NUMNO) STOP 'NUPP GT NUMNO IN LNCO2Q' 
                 IF (NUPP.GT.0) RUPP=RATNO(NUPP) 
             ELSE
                 RLOW=0.
                 RUPP=0.
             END IF
C
C     RLOW AND RUPP NOW SET
C
             FNLTE=SP(IJ)/(1.0-DELTA)
             SABS(IJ)=FNLTE*(RLOW-RUPP*DELTA)
             SRAD(IJ)=FNLTE*(RLOW-RUPP)
          ENDIF

         IF (VNU(IJ).GT.VHI) THEN                                         D02160
            IEOF = 1                                                      D02170
            GO TO 60                                                      D02180
         ENDIF                                                            D02190

   50 CONTINUE                                                            D02200
      IF (IJ.LT.LIMIN.AND.IEOF.EQ.0) THEN
         CALL CPUTIM (TIM2)
         TIMS = TIMS+TIM2-TIM1
         GO TO 30                                                         D02210
      ENDIF
   60 CALL CPUTIM (TIM2)                                                  D02220
      IHI = IJ                                                            D02230
      TIMS = TIMS+TIM2-TIM1                                               D02240
C                                                                         D02250
      CALL CPUTIM(TPAT0)
      CALL SHRINQ                                                         D02260
      CALL CPUTIM(TPAT1)
      TSHRNK = TSHRNK + TPAT1-TPAT0
      IJ = ILO-1                                                          D02270
      IF (IHI.LT.LIMIN.AND.IEOF.EQ.0) GO TO 30                            D02280
C                                                                         D02290
      VNULO = VNU(1)                                                      D02300
      VNUHI = VNU(IHI)                                                    D02310
      JLIN = IHI                                                          D02320
C                                                                         D02330
      IF (JLIN.GT.0) THEN                                                 D02340
         CALL CPUTIM(TPAT0)
         CALL BUFOUT (LNFIL4,RCDHDR(1),NPHDRL)                            D02350
         CALL BUFOUT (LNFIL4,VNU(1),NLNGT4)                               D02360
         CALL CPUTIM(TPAT1)
         TBUFFR = TBUFFR + TPAT1-TPAT0
      ENDIF                                                               D02370
      NLINS = NLINS+IHI-IST+1                                             D02380
C                                                                         D02390
      IF (IEOF.EQ.1) GO TO 70                                             D02400
      IJ = 0                                                              D02410
      ILO = 1                                                             D02420
      GO TO 30                                                            D02430
   70 CONTINUE                                                            D02440
C                                                                         D02450
      DO 80 M = 1, NMOL                                                   D02460
         IF (MEFDP(M).GT.0) WRITE (IPR,905) MEFDP(M),M                    D02470
   80 CONTINUE                                                            D02480
      CALL CPUTIM (TIMEL1)                                                D02490
      TIME = TIMEL1-TIMEL0                                                D02500
      IF (NOPR.EQ.0) THEN
         WRITE (IPR,910) TIME,TIMR,TIMS,NLIN,NLINS                        D02510
         L4TIM=TIME
         L4TMR=TIMR
         L4TMS=TIMS
         L4NLN=NLIN
         L4NLS=NLINS
         LOTHER = TSHRNK+TBUFFR+TMOLN4
      ENDIF
      RETURN                                                              D02520
C                                                                         D02530
  900 FORMAT ('0  *****  LINF4 - VNU LIMITS DO NOT INTERSECT WITH ',      D02540
     *        'LINFIL - LINF4 NOT USED *****',/,'   VNU = ',F10.3,        D02550
     *        ' - ',F10.3,' CM-1     LINFIL = ',F10.3,' - ',F10.3,        D02560
     *        ' CM-1')                                                    D02570
  905 FORMAT ('0*************************',I5,' STRENGTHS FOR',           D02580
     *        '  TRANSITIONS WITH UNKNOWN EPP FOR MOL =',I5,              D02590
     *        ' SET TO ZERO')                                             D02600
  910 FORMAT ('0',20X,'TIME',11X,'READ',9X,'SHRINQ',6X,'NO. LINES',3X,    D02610
     *        'AFTER SHRINQ',/,2X,'LINF4 ',2X,3F15.3,2I15)                D02620
C                                                                         D02630
      END                                                                 D02640
c ----------------------------------------------------------------
      SUBROUTINE HIRACQ (MPTS)                                            B00010
C                                                                         B00020
      IMPLICIT REAL*8           (V)                                     ! B00030
C                                                                         B00040
C                                                                         B00050
C**********************************************************************   B00060
C*                                                                        B00070
C*                                                                        B00080
C*    CALCULATES MONOCHROMATIC ABSORPTION COEFFICIENT FOR SINGLE LAYER    B00090
C*                                                                        B00100
C*                                                                        B00110
C*            USES APPROXIMATE VOIGT ALGORITHM                            B00120
C*                                                                        B00130
C*                                                                        B00140
C*              VAN VLECK WEISSKOPF LINE SHAPE                            B00150
C*                                                                        B00160
C**********************************************************************   B00170
C                                                                         B00180
C                                                                         B00190
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B00200
C                                                                         B00230
C                  IMPLEMENTATION:    R.D. WORSHAM                        B00240
C                                                                         B00250
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         B00260
C                                     R.D. WORSHAM                        B00270
C                                     J.L. MONCET                         B00280
C                                                                         B00290
C                                                                         B00300
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B00310
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          B00320
C                                                                         B00330
C----------------------------------------------------------------------   B00340
C                                                                         B00350
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     B00360
C                                     OFFICE OF ENERGY RESEARCH           B00370
C                                     DEPARTMENT OF ENERGY                B00380
C                                                                         B00390
C                                                                         B00400
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             B00410
C                                                                         B00420
C                                             FASCOD3                     B00430
C                                                                         B00440
C                                                                         B00450
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B00460
C                                                                         B00470
C                                                                         B00480
C     Common blocks from analytic derivatives
C     -------------------------
      COMMON /ADRPNM/ CDUM1,PTHODI,PTHODT,PTHRDR
C     -------------------------
      COMMON /RCNTRL/ ILNFLG
      COMMON VNU(250),SP(250),ALFA0(250),EPP(250),MOL(250),HWHMS(250),    B00490
     *       TMPALF(250),PSHIFT(250),IFLG(250),SPPSP(250),RECALF(250),    B00500
     *       ZETAI(250),IZETA(250)                                        B00510
C
C     DIMENSION RR1 =  NBOUND   + 1 + DIM(R1)
C     DIMENSION RR2 =  NBOUND/2 + 1 + DIM(R2)
C     DIMENSION RR3 =  NBOUND/4 + 1 + DIM(R3)
C
      COMMON RR1(6099),RR2(2075),RR3(429)                                 B00520
      COMMON /XRNLTE/ RR1s(4050),RR2s(1050),RR3s(300)               
      COMMON /IOU/ IOUT(250)                                              B00530
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(2030)                B00540
      COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,NOZERO,NP,H1F,H2F,         B00550
     *                ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH,         B00560
     *                XVBAR, HMINF,PHIF,IERRF,HSPACE                      B00570
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   B00580
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       B00590
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       B00600
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      B00610
     *              EXTID(10)                                             B00620
C                                                                         B00630
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ
C                                                                         B00650
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       B00660
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   B00670
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    B00680
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /XSUB/ VBOT,VTOP,VFT,LIMIN,ILO,IHI,IEOF,IPANEL,ISTOP,IDATA   B00700
      COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)        B00710
      COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,            B00720
     *                HWF3,DXF3,NX3,N3MAX                                 B00730
      COMMON /SUB1/ MAX1,MAX2,MAX3,NLIM1,NLIM2,NLIM3,NLO,NHI,DVR2,DVR3,   B00740
     *              N1R1,N2R1,N1R2,N2R2,N1R3,N2R3                         B00750
      COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTS       B00760
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,         B00770
     *               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL                      B00780
      COMMON /VOICOM/ AVRAT(102),CGAUSS(102),CF1(102),CF2(102),           B00790
     *                CF3(102),CER(102)                                   B00800
C
      PARAMETER (NFPTS=2001,NFMX=1.3*NFPTS, NUMZ = 101)  
c
      COMMON /FNSHQ/ IFN,F1(NFMX, NUMZ),F2(NFMX, NUMZ),
     $     F3(NFMX, NUMZ), FG(NFMX)

      COMMON /R4SUB/ VLOF4,VHIF4,ILOF4,IST,IHIF4,LIMIN4,LIMOUT,ILAST,     B00820
     *               DPTMN4,DPTFC4,ILIN4,ILIN4T                           B00830
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         B00840
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        B00850
     *              NLTEFL,LNFIL4,LNGTH4                                  B00860
C                                                                         B00870
      PARAMETER (NTMOL=36,NSPECI=85)   
C                                                                         B00890
      COMMON /ISVECT/ ISOVEC(NTMOL),ISO82(NSPECI),ISONM(NTMOL),           B00900
     *                SMASSI(NSPECI)                                      B00910
      COMMON /LNC1/ RHOSLF(NSPECI),ALFD1(NSPECI),SCOR(NSPECI),ALFMAX,     B00920
     *              BETACR,DELTMP,DPTFC,DPTMN,XKT,NMINUS,NPLUS,NLIN,      B00930
     *              LINCNT,NCHNG,SUMALF,SUMZET,TRATIO,RHORAT,PAVP0,       B00940
     *              PAVP2,RECTLC,TMPDIF,ILC                               B00950
      COMMON /FLFORM/ CFORM                                               B00960
      COMMON /L4TIMG/ L4TIM,L4TMR,L4TMS,L4NLN,L4NLS,LOTHER
      COMMON /IODFLG/ DVOUT
c     Total timing array for layer line-by-line calculation
      common /timing_lay_nlte/ time_lay_lbl(20)
C                                                                         B00970
      REAL L4TIM,L4TMR,L4TMS,LOTHER
      CHARACTER*55 CDUM1,PTHODI,PTHODT,PTHRDR
      CHARACTER*10 HFMODL
      CHARACTER CFORM*11,KODLYR*57,PTHODE*55,PTHODD*55                    B00980
      LOGICAL OP                                                          B00990
C                                                                         B01000
      DIMENSION MEFDP(64),FILHDR(2),IWD(2)                                B01010
      DIMENSION R1(4050),R2(1050),R3(300)
      DIMENSION SABS(2),SRAD(2)
C                                                                         B01020
      EQUIVALENCE (IHIRAC,FSCDID(1)) , (ILBLF4,FSCDID(2)),                B01030
     *            (IXSCNT,FSCDID(3)) , (IAERSL,FSCDID(4)),                B01040
     *            (JRAD,FSCDID(9)) , (IMRG,FSCDID(11)),                   B01050
     *            (IATM,FSCDID(15)) , (YI1,IOD) , (XID(1),FILHDR(1)),     B01060
     *            (V1P,IWD(1)) , (NPNLXP,LSTWDX),                         B01070
     &            (SP(1),SABS(1)),(EPP(1),SRAD(1))
      EQUIVALENCE (R1(1), RR1(2049)),(R2(1),RR2(1025)),(R3(1),RR3(129))
C                                                                         B01080
C
C     NOTE that DXFF1 = (HWFF1/(NFPTS-1))
C     and       DXFF2 = (HWFF2/(NFPTS-1))
C     and       DXFF3 = (HWFF3/(NFPTS-1))
C
      DATA HWFF1 /  4. /,DXFF1 / 0.002 /,NXF1 / NFPTS /,NF1MAX / NFMX /   B01090
      DATA HWFF2 / 16. /,DXFF2 / 0.008 /,NXF2 / NFPTS /,NF2MAX / NFMX /   B01100
      DATA HWFF3 / 64. /,DXFF3 / 0.032 /,NXF3 / NFPTS /,NF3MAX / NFMX /   B01110
C                                                                         B01120
      DATA MEFDP / 64*0 /                                                 B01130
C                                                                         B01140
      PTHODI = 'ODint_'
      PTHODE = 'ODexact_'
      PTHODD = 'ODdeflt_'
      DATA KODLYR /
     *     '                                                         '/
      DATA HFMODL /'         '/
C                                                                         B01160
      CALL CPUTIM (TIMEH0)                                                B01170
C                                                                         B01180
C     ASSIGN SCCS VERSION NUMBER TO MODULE 
C
      HVROPR = '????'
C
C     Initialize timing for the group "OTHER" in the TAPE6 output
C
      TLNCOR = 0.0
      TXINT = 0.0
      TSHAPE = 0.0
      TLOOPS = 0.0
      TODFIL = 0.0
      TMOLEC = 0.0
C
      LSTWDX = -654321
      NPNLXP = NWDL(IWD,LSTWDX)                                           B01190
      ICNTNM = MOD(IXSCNT,10)                                             B01200
      IXSECT = IXSCNT/10                                                  B01210
C                                                                         B01220
C     SET INPUT FLAG FOR USE BY X-SECTIONS                                B01230
C                                                                         B01240
      IFST = -99                                                          B01250
      IR4 = 0                                                             B01260
      IENTER = 0                                                          B01270
C                                                                         B01280
C     SET COMMON BLOCK CMSHAP                                             B01290
C                                                                         B01300
      HWF1 = HWFF1                                                        B01310
      DXF1 = DXFF1                                                        B01320
      NX1 = NXF1                                                          B01330
      N1MAX = NF1MAX                                                      B01340
      HWF2 = HWFF2                                                        B01350
      DXF2 = DXFF2                                                        B01360
      NX2 = NXF2                                                          B01370
      N2MAX = NF2MAX                                                      B01380
      HWF3 = HWFF3                                                        B01390
      DXF3 = DXFF3                                                        B01400
      NX3 = NXF3                                                          B01410
      N3MAX = NF3MAX                                                      B01420
C                                                                         B01430
      DPTMN = DPTMIN                                                      B01440
      IF (JRAD.NE.1) DPTMN = DPTMIN/RADFN(V2,TAVE/RADCN2)                 B01450
      DPTFC = DPTFAC                                                      B01460
      ILIN4 = 0                                                           B01470
      ILIN4T = 0                                                          B01480
      NPTS = MPTS                                                         B01490
      LIMIN = 250                                                         B01500
      NSHIFT = 32                                                         B01510
C                                                                         B01520
C     SAMPLE IS AVERAGE ALPHA / DV                                        B01530
C                                                                         B01540
      NBOUND = 4.*(2.*HWF3)*SAMPLE+0.01                                   B01550
      NLIM1 = 2401                                                        B01560
      NLIM2 = (NLIM1/4)+1                                                 B01570
      NLIM3 = (NLIM2/4)+1                                                 B01580
C                                                                         B01590
      IF (IFN.EQ.0) THEN                                                  B01600
         CALL CPUTIM(TPAT0)
         call voigt_init(f1, f2, f3)
         IFN = IFN+1                                                      B01640
         CALL CPUTIM(TPAT1)
         TSHAPE = TSHAPE+TPAT1-TPAT0
      ENDIF                                                               B01650
C                                                                         B01660
      CALL CPUTIM(TPAT0)
      CALL MOLEC (1,SCOR,RHOSLF,ALFD1)                                    B01670
      CALL CPUTIM(TPAT1)
      TMOLEC = TMOLEC+TPAT1-TPAT0
      REWIND LINFIL                                                       B01680
      TIMRDF = 0.                                                         B01690
      TIMCNV = 0.                                                         B01700
      TIMPNL = 0.                                                         B01710
      TF4 = 0.                                                            B01720
      TF4RDF = 0.                                                         B01730
      TF4CNV = 0.                                                         B01740
      TF4PNL = 0.                                                         B01750
      TXS = 0.                                                            B01760
      TXSRDF = 0.                                                         B01770
      TXSCNV = 0.                                                         B01780
      TXSPNL = 0.                                                         B01790
      IEOF = 0                                                            B01800
      ILO = 0                                                             B01810
      IHI = -999                                                          B01820
      NMINUS = 0                                                          B01830
      NPLUS = 0                                                           B01840
C                                                                         B01850
C     NOTE (DXF3/DXF1) IS 16 AND (DXF3/DXF2) IS 4                         B01860
C                                                                         B01870
      DVP = DV                                                            B01880
      DVR2 = (DXF2/DXF1)*DV                                               B01890
      DVR3 = (DXF3/DXF1)*DV                                               B01900
      MAX1 = NSHIFT+NLIM1+(NBOUND/2)                                      B01910
      MAX2 = MAX1/4                                                       B01920
      MAX3 = MAX1/16                                                      B01930
      MAX1 = MAX1+NSHIFT+1+16                                             B01940
      MAX2 = MAX2+NSHIFT+1+4                                              B01950
      MAX3 = MAX3+NSHIFT+1+1                                              B01960
C                                                                         B01970
C     FOR CONSTANTS IN PROGRAM  MAX1=4018  MAX2=1029  MAX3=282            B01980
C                                                                         B01990
      CALL CPUTIM(TPAT0)
      BOUND = FLOAT(NBOUND)*DV/2.                                         B02000
      BOUNF3 = BOUND/2.                                                   B02010
      ALFMAX = BOUND/HWF3                                                 B02020
      NLO = NSHIFT+1                                                      B02030
      NHI = NLIM1+NSHIFT-1                                                B02040
      DO 10 I = 1, MAX1                                                   B02050
         R1(I) = 0.                                                       B02060
         RR1s(I) = 0.                                                     B02060
   10 CONTINUE                                                            B02070
      DO 20 I = 1, MAX2                                                   B02080
         R2(I) = 0.                                                       B02090
         RR2s(I) = 0.                                                     B02090
   20 CONTINUE                                                            B02100
      DO 30 I = 1, MAX3                                                   B02110
         R3(I) = 0.                                                       B02120
         RR3s(I) = 0.                                                     B02120
   30 CONTINUE                                                            B02130
      IF (ILBLF4.EQ.0) THEN                                               B02140
         DO 40 I = 1, 2502                                                B02150
            R4(I) = 0.                                                    B02160
            RR4(I) = 0.                                                   B02160
   40    CONTINUE                                                         B02170
      ENDIF                                                               B02180
C                                                                         B02190
      IF (IATM.GE.1.AND.IATM.LE.5) CALL YDIH1 (H1F,H2F,ANGLEF,YID)        B02200
      CALL CPUTIM(TPAT1)
      TLOOPS = TLOOPS + TPAT1-TPAT0
C                                                                         B02210
C     ---------------------------------------------------------------
C
C     - If IOD = 1 or 4 then calculate optical depths for each
C       layer with DV = DVOUT (using DVSET if IOD=4) and maintain
C       separately. Use PTHODI as the name of the optical depth files.
C       This requires the format HFMODL, which is produced by
C       calling the SUBROUTINE QNTIFY.
C
C     - If IOD = 2 and IMERGE = 1 then calculate optical depths
C       for each layer using the exact DV of each layer
C       Use PTHODE as the name of the optical depth files.
C       This requires the format HFMODL, which is produced by
C       calling the SUBROUTINE QNTIFY.
C
C     - If calculating layer optical depths and cumulative layer
C       optical depths for an analytic derivative calculation
C       (IOD=3, IMRG=10), or when using the same criteria but not
C       calculating the cumulative optical depths (IOD=3),
C       then use PTHODI as the name of the optical depth files.
C       This requires the format HFMODL, which is produced by
C       calling the SUBROUTINE QNTIFY.
C
C     - If calculating layer absorptance coefficients for an
C       analytic derivative calculation (IEMIT=3, IOD=3, and
C       IMRG>40), then use TAPE10 as the name of the layer
C       absorptance coefficient files.
C
C     - If calculating optical depths using the default procedure,
C       sending output to a different file for each layer (IEMIT=0,
C       IOD=0, and IMRG=1), then use PTHODI for the optical depth
C       pathnames.
C
C     - Otherwise, use TAPE10.  For IOD=1, calculate optical depths
C       for each layer with DV = DVOUT (from DVSET in TAPE5, carried
C       in by COMMON BLOCK /IODFLG/ (interpolation in PNLINT).
C
      CALL CPUTIM(TPAT0)
      IF ((IOD.EQ.1).OR.(IOD.EQ.4)) THEN
         CALL QNTIFY(PTHODI,HFMODL)
         WRITE (KODLYR,HFMODL) PTHODI,LAYER                               B02230
         INQUIRE (UNIT=KFILE,OPENED=OP)                                   B02240
         IF (OP) CLOSE (KFILE)                                            B02250
         OPEN (KFILE,FILE=KODLYR,FORM=CFORM,STATUS='UNKNOWN')             B02260
         REWIND KFILE                                                     B02270
         DVSAV = DV
         IF (DVOUT.NE.0.) DV = DVOUT
         CALL BUFOUT (KFILE,FILHDR(1),NFHDRF)                             B02310
         DV = DVSAV
         IF (NOPR.EQ.0) WRITE (IPR,900) KFILE,DV,BOUNF3                   B02320
      ELSEIF (IOD.EQ.2) THEN
         CALL QNTIFY(PTHODE,HFMODL)
         WRITE(KODLYR,HFMODL) PTHODE,LAYER
         INQUIRE (UNIT=KFILE,OPENED=OP)
         IF (OP) CLOSE (KFILE)
         OPEN (KFILE,FILE=KODLYR,FORM=CFORM,STATUS='UNKNOWN')
         REWIND KFILE
         DVOUT = DV
         CALL BUFOUT (KFILE,FILHDR(1),NFHDRF)
         IF (NOPR.EQ.0) WRITE (IPR,900) KFILE,DVOUT,BOUNF3
      ELSEIF (IOD.EQ.3) THEN
         IF ((IMRG.EQ.10).OR.(IMRG.EQ.1)) THEN
            CALL QNTIFY(PTHODI,HFMODL)
            WRITE (KODLYR,HFMODL) PTHODI,LAYER
            INQUIRE (UNIT=KFILE,OPENED=OP)
            IF (OP) CLOSE (KFILE)
            OPEN (KFILE,FILE=KODLYR,FORM=CFORM,STATUS='UNKNOWN')
            REWIND KFILE
            DVSAV = DV
            DV = DVOUT
            CALL BUFOUT (KFILE,FILHDR(1),NFHDRF)
            DV = DVSAV
            IF (NOPR.EQ.0) WRITE (IPR,900) KFILE,DV,BOUNF3
         ELSEIF (IMRG.GE.40) THEN
            DVSAV = DV
            DV = DVOUT
            CALL BUFOUT (KFILE,FILHDR(1),NFHDRF)
            DV = DVSAV
            IF (NOPR.EQ.0) WRITE (IPR,900) KFILE,DV,BOUNF3
         ENDIF
      ELSE
         IF (IMRG.EQ.1) THEN
            CALL QNTIFY(PTHODD,HFMODL)
            WRITE (KODLYR,HFMODL) PTHODD,LAYER
            INQUIRE (UNIT=KFILE,OPENED=OP)
            IF (OP) CLOSE (KFILE)
            OPEN (KFILE,FILE=KODLYR,FORM=CFORM,STATUS='UNKNOWN')
            REWIND KFILE
         ENDIF
         IF (IOD.EQ.1) THEN
            DVSAV = DV
            IF (DVOUT.NE.0.) DV = DVOUT
            CALL BUFOUT (KFILE,FILHDR(1),NFHDRF)
            DV = DVSAV
         ELSE
            DVOUT = 0.0                                                   B02350
            CALL BUFOUT (KFILE,FILHDR(1),NFHDRF)                          B02360
         ENDIF
         IF (NOPR.EQ.0) WRITE (IPR,900) KFILE,DV,BOUNF3                   B02370
      ENDIF                                                               B02380
      CALL CPUTIM(TPAT1)
      TODFIL = TODFIL + TPAT1-TPAT0
C                                                                         B02390
      IF (IHIRAC.EQ.9) THEN                                               B02400
         DO 50 M = 1, NMOL                                                B02410
            WK(M) = 0.                                                    B02420
   50    CONTINUE                                                         B02430
      ENDIF                                                               B02440
C     
C     ---------------------------------------------------------------
C                                                                         B02450
      VFT = V1-FLOAT(NSHIFT)*DV                                           B02460
      VBOT = V1-BOUND                                                     B02470
      VTOP = V2+BOUND                                                     B02480
C CHECK TO SEE IF VBOT OR VTOP FALLS IN LINE COUPLING REGION
C     CALL CHKLNC(VBOT,VTOP,(SAMPLE*DV),1)
C                                                                         B02490
      LINCNT = 0                                                          B02500
      NLIN = 0                                                            B02510
      AVALF = 0.                                                          B02520
      AVZETA = 0.                                                         B02530
      SUMALF = 0.                                                         B02540
      SUMZET = 0.                                                         B02550
      NCHNG = 0                                                           B02560
      NLNCR = 0                                                           B02570
C                                                                         B02580
      V1R4ST = V1R4                                                       B02590
      V2R4ST = V2R4                                                       B02600
      IF (ILBLF4.GE.1) CALL LBLFQ (JRAD,V1R4ST,V2R4ST)                    B02610
C                                                                         B02620
      IFPAN = 1                                                           B02630
C                                                                         B02640
   60 CONTINUE                                                            B02650
C                                                                         B02660
      CALL CPUTIM (TIME0)                                                 B02670
      IF (IEOF.NE.0) GO TO 80                                             B02680
C                                                                         B02690
C     THERE ARE (LIMIN * 9) QUANTITIES READ IN:                           B02700
C     VNU,SP,ALFA0,EPP,MOL,HWHMS,TMPALF,PSHIFT,IFLG                       B02710
C                                                                         B02720
      CALL RDLIN                                                          B02730
C                                                                         B02740
      CALL CPUTIM (TIME)                                                  B02750
      TIMRDF = TIMRDF+TIME-TIME0                                          B02760
C                                                                         B02770
      IF (IEOF.NE.0) GO TO 80                                             B02780
C                                                                         B02790
C     MODIFY LINE DATA FOR TEMPERATURE, PRESSURE, AND COLUMN DENSITY      B02800
C                                                                         B02810
      CALL CPUTIM(TPAT0)
      CALL LNCORQ (NLNCR,IHI,ILO,MEFDP)                                   B02820
      CALL CPUTIM(TPAT1)
      TLNCOR = TLNCOR+TPAT1-TPAT0
C                                                                         B02830
   70 CONTINUE                                                            B02840
C                                                                         B02850
      CALL CNVFNQ (VNU,SABS,SRAD,SPPSP,RECALF,R1,R2,R3,RR1s,RR2s,RR3s,
     &    ZETAI,IZETA)        
C                                                                         B02880
      IF (IPANEL.EQ.0) GO TO 60                                           B02890
C                                                                         B02900
   80 CONTINUE                                                            B02910
C                                                                         B02920
C        FOR FIRST PANEL     N1R1=   1    N1R2=  1    N1R3=  1            B02930
C     FOR SUBSEQUENT PANELS  N1R1=  33   *N1R2= 13   *N1R3=  6            B02940
C         FOR ALL PANELS     N2R1=2432   *N2R2=612   *N2R3=155            B02950
C                                                                         B02960
C            NOTE: THE VALUES FOR N1R2, N1R3, N2R2 AND N2R3 WHICH         B02970
C                  ARE MARKED WITH AN ASTERISK, CONTAIN A 4 POINT         B02980
C                  OFFSET WHICH PROVIDES THE NECESSARY OVERLAP FOR        B02990
C                  THE INTERPOLATION OF R3 INTO R2, AND R2 INTO R1.       B03000
C                                                                         B03010
      IF (IFPAN.EQ.1) THEN                                                B03020
         IFPAN = 0                                                        B03030
         N1R1 = 1                                                         B03040
         N1R2 = 1                                                         B03050
         N1R3 = 1                                                         B03060
      ELSE                                                                B03070
         N1R1 = NSHIFT+1                                                  B03080
         N1R2 = (NSHIFT/4)+1+4                                            B03090
         N1R3 = (NSHIFT/16)+1+3                                           B03100
      ENDIF                                                               B03110
      N2R1 = NLIM1+NSHIFT-1                                               B03120
      N2R2 = NLIM2+(NSHIFT/4)-1+4                                         B03130
      N2R3 = NLIM3+(NSHIFT/16)-1+3                                        B03140
C                                                                         B03150
      IF (VFT.LE.0.) THEN                                                 B03160
         CALL RSYM (R1,DV,VFT)                                            B03170
         CALL RSYM (R2,DVR2,VFT)                                          B03180
         CALL RSYM (R3,DVR3,VFT)                                          B03190
         CALL RSYM (RR1s,DV,VFT)                                          B03170
         CALL RSYM (RR2s,DVR2,VFT)                                        B03180
         CALL RSYM (RR3s,DVR3,VFT)                                        B03190
      ENDIF                                                               B03200
C                                                                         B03210
      IF (IXSECT.GE.1.AND.IR4.EQ.0) THEN                                  B03220
         CALL CPUTIM (TIME0)                                              B03230
         CALL XSECTM (IFST,IR4)                                           B03240
         CALL CPUTIM (TIME)                                               B03250
         TXS = TXS+TIME-TIME0                                             B03260
      ENDIF                                                               B03270
C                                                                         B03280
      CALL CPUTIM(TPAT0)
      IF (ILBLF4.GE.1) THEN                                               B03290
          CALL XINT (V1R4,V2R4,DVR4,R4,1.0,VFT,DVR3,R3,N1R3,N2R3)         B03300
          CALL XINT (V1R4,V2R4,DVR4,RR4,1.0,VFT,DVR3,RR3s,N1R3,N2R3)      B03300
      ENDIF
      IF (ICNTNM.GE.1)                                                    B03310
     *    CALL XINT (V1ABS,V2ABS,DVABS,ABSRB,1.,VFT,DVR3,R3,N1R3,N2R3)    B03320
      CALL CPUTIM(TPAT1)
      TXINT = TXINT + TPAT1-TPAT0
C                                                                         B03330
      CALL PANELQ (R1,R2,R3,RR1s,RR2s,RR3s,KFILE,JRAD,IENTER)             B03340
C                                                                         B03350
      IF (ISTOP.NE.1) THEN                                                B03360
         IF (ILBLF4.GE.1) THEN                                            B03370
            VF1 = VFT-2.*DVR4                                             B03380
            VF2 = VFT+2.*DVR4+FLOAT(N2R3+4)*DVR3                          B03390
            IF (VF2.GT.V2R4.AND.V2R4.NE.V2R4ST) THEN                      B03400
               CALL LBLFQ (JRAD,VF1,V2R4ST)                               B03410
               IF (IXSECT.GE.1.AND.IR4.EQ.1) THEN                         B03420
                  CALL CPUTIM (TIME0)                                     B03430
                  CALL XSECTM (IFST,IR4)                                  B03440
                  CALL CPUTIM (TIME)                                      B03450
                  TXS = TXS+TIME-TIME0                                    B03460
               ENDIF                                                      B03470
            ENDIF                                                         B03480
         ENDIF                                                            B03490
         GO TO 70                                                         B03500
      ENDIF                                                               B03510
C                                                                         B03520
      CALL CPUTIM (TIMEH1)                                                B03530
      TIME = TIMEH1-TIMEH0-TF4-TXS                                        B03540
C                                                                         B03550
      IF (NOPR.NE.1) THEN                                                 B03560
         IF (ILBLF4.GE.1) WRITE (IPR,905) DVR4,BOUND4                     B03570
         IF (NMINUS.GT.0) WRITE (IPR,910) NMINUS                          B03580
         IF (NPLUS.GT.0) WRITE (IPR,915) NPLUS                            B03590
         TOTHHI = TLNCOR+TXINT+TSHAPE+TLOOPS+TODFIL+TMOLEC
         WRITE (IPR,920) L4TIM,L4TMR,L4TMS,LOTHER,L4NLN,L4NLS,
     *                   TXS,TXSRDF,TXSCNV,TXSPNL,                        B03600
     *                   TF4,TF4RDF,TF4CNV,TF4PNL,ILIN4T,ILIN4,           B03610
     *                   TIME,TIMRDF,TIMCNV,TIMPNL,TOTHHI,
     *                   NLIN,LINCNT,NCHNG                                B03620
c         time_lay_lbl(1) = l4tim
c         time_lay_lbl(2) = l4tmr
c         time_lay_lbl(3) = 0.0
c         time_lay_lbl(4) = l4tms
c         time_lay_lbl(5) = lother
c         time_lay_lbl(6) = txs  
c         time_lay_lbl(7) = txsrdf
c         time_lay_lbl(8) = txscnv
c         time_lay_lbl(9) = txspnl
c         time_lay_lbl(10) = 0.0
c         time_lay_lbl(11) = tf4   
c         time_lay_lbl(12) = tf4rdf
c         time_lay_lbl(13) = tf4cnv
c         time_lay_lbl(14) = tf4pnl
c         time_lay_lbl(15) = 0.0
c         time_lay_lbl(16) = time
c         time_lay_lbl(17) = timrdf
c         time_lay_lbl(18) = timcnv
c         time_lay_lbl(19) = timpnl
c         time_lay_lbl(20) = tothhi

c        Accumulate timing array

         DATA  i_time_lay/454545/
         
         if (i_time_lay .eq. 454545) then
            i_time_lay = 676767
            do ltime = 1,20
               time_lay_lbl(ltime) = 0.0
            enddo
         endif

         time_lay_lbl(1) = time_lay_lbl(1) + l4tim
         time_lay_lbl(2) = time_lay_lbl(2) + l4tmr
         time_lay_lbl(3) = 0.0
         time_lay_lbl(4) = time_lay_lbl(4) + l4tms
         time_lay_lbl(5) = time_lay_lbl(5) + lother
         time_lay_lbl(6) = time_lay_lbl(6) + txs  
         time_lay_lbl(7) = time_lay_lbl(7) + txsrdf
         time_lay_lbl(8) = time_lay_lbl(8) + txscnv
         time_lay_lbl(9) =  time_lay_lbl(9) + txspnl
         time_lay_lbl(10) = 0.0
         time_lay_lbl(11) = time_lay_lbl(11) + tf4   
         time_lay_lbl(12) = time_lay_lbl(12) + tf4rdf
         time_lay_lbl(13) = time_lay_lbl(13) + tf4cnv
         time_lay_lbl(14) = time_lay_lbl(14) + tf4pnl
         time_lay_lbl(15) = 0.0
         time_lay_lbl(16) = time_lay_lbl(16) + time
         time_lay_lbl(17) = time_lay_lbl(17) + timrdf
         time_lay_lbl(18) = time_lay_lbl(18) + timcnv
         time_lay_lbl(19) = time_lay_lbl(19) + timpnl
         time_lay_lbl(20) = time_lay_lbl(20)+ tothhi

         IF (LAYER.EQ.NLAYRS) THEN 

             WRITE (IPR,*) 'Total Accumulated Times'
             WRITE (IPR,922) time_lay_lbl(1),time_lay_lbl(2),
     *           time_lay_lbl(4), time_lay_lbl(5),
     *           time_lay_lbl(6),
     *           time_lay_lbl(7),time_lay_lbl(8),time_lay_lbl(9),
     *           time_lay_lbl(11),time_lay_lbl(12),time_lay_lbl(13),                        
     *           time_lay_lbl(14),time_lay_lbl(16),
     *           time_lay_lbl(17),     
     *           time_lay_lbl(18),time_lay_lbl(19),time_lay_lbl(20)
                           
         ENDIF

         WRITE(IPR,935)
         IF (LINCNT.GE.1) THEN                                            B03630
            AVALF = SUMALF/FLOAT(LINCNT)                                  B03640
            AVZETA = SUMZET/FLOAT(LINCNT)                                 B03650
         ENDIF                                                            B03660
         WRITE (IPR,925) AVALF,AVZETA                                     B03670
C                                                                         B03680
         DO 90 M = 1, NMOL                                                B03690
            IF (MEFDP(M).GT.0) WRITE (IPR,930) MEFDP(M),M                 B03700
   90    CONTINUE                                                         B03710
      ENDIF                                                               B03720
C                                                                         B03730
      RETURN                                                              B03740
C                                                                         B03750
  900 FORMAT ('0  * HIRAC1 *  OUTPUT ON FILE ',I5,10X,' DV = ',F12.8,     B03760
     *        10X,' BOUNDF3(CM-1) = ',F8.4)                               B03770
  905 FORMAT ('0 DV FOR LBLF4 = ',F10.5,5X,'BOUND FOR LBLF4 =',F10.4)     B03780
  910 FORMAT ('0 -------------------------',I5,' HALF WIDTH CHANGES')     B03790
  915 FORMAT ('0 +++++++++++++++++++++++++',I5,' HALF WIDTH CHANGES')     B03800
  920 FORMAT ('0',20X,'TIME',11X,'READ',4X,'CONVOLUTION',10X,'PANEL',     B03810
     *        9X,'OTHER+',
     *        6X,'NO. LINES',3X,'AFTER REJECT',5X,'HW CHANGES',/,         B03820
     *        2x,'LINF4',3X,2F15.3,15X,2F15.3,2I15,/,
     *        2X,'XSECT ',2X,4F15.3,/,2X,'LBLF4 ',2X,4F15.3,15X,2I15,/,   B03830
     *        2X,'HIRAC1',2X,5F15.3,3I15)                                 B03840
  922 FORMAT ('0',20X,'TIME',11X,'READ',4X,'CONVOLUTION',10X,'PANEL',     
     *        9X,'OTHER+',/,      
     *        2x,'LINF4',3X,2F15.3,15X,2F15.3,/,
     *        2X,'XSECT ',2X,4F15.3,/,2X,'LBLF4 ',2X,4F15.3,15X,/,   
     *        2X,'HIRAC1',2X,5F15.3)                                 
  925 FORMAT ('0  * HIRAC1 *  AVERAGE WIDTH = ',F8.6,                     B03850
     *        ',  AVERAGE ZETA = ',F8.6)                                  B03860
  930 FORMAT ('0 ********  HIRAC1  ********',I5,' STRENGTHS FOR',         B03870
     *        '  TRANSITIONS WITH UNKNOWN EPP FOR MOL =',I5,              B03880
     *        ' SET TO ZERO')                                             B03890
 935  FORMAT (/,'0     + OTHER timing includes:',/,
     *          '0             In LINF4:  MOLEC, BUFIN, BUFOUT, ',
     *          'NWDL, ENDFIL, and SHRINQ',/,
     *          '0             In HIRAC:  LNCOR, XINT, SHAPEL, ',
     *          'SHAPEG, VERFN, MOLEC, and other loops and ',
     *          'file maintenance within HIRAC',/)
C                                                                         B03900
      END                                                                 B03910
      SUBROUTINE LNCORQ (NLNCR,IHI,ILO,MEFDP)                             B04840
C                                                                         B04850
      IMPLICIT REAL*8           (V)                                     ! B00030
C                                                                         B04870
      CHARACTER*1 FREJ(250),HREJ,HNOREJ
      COMMON /RCNTRL/ ILNFLG
      COMMON VNU(250),S(250),ALFA0(250),EPP(250),MOL(250),HWHMS(250),     B04880
     *       TMPALF(250),PSHIFT(250),IFLG(250),SPPSP(250),RECALF(250),    B04890
     *       ZETAI(250),IZETA(250)                                        B04900
      COMMON /IOU/ IOUT(250)                                              B04920
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   B04930
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       B04940
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       B04950
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      B04960
     *              EXTID(10)                                             B04970
C                                                                         B04980
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ

C                                                                         B05000
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       B05010
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   B05020
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    B05030
      COMMON /XSUB/ VBOT,VTOP,VFT,DUM(7)
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)        B05050
      COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,
     *                HWF3,DXF3,NX3,N3MAX
      COMMON /VOICOM/ AVRAT(102),CGAUSS(102),CF1(102),CF2(102),           B05060
     *                CF3(102),CER(102)                                   B05070

      COMMON /VBNLTE/ RATH2O(8),RATCO2(26),RATO3(18),RATCO(3),RATNO(3),   603420
     &               NUMH2O,NUMCO2,NUMO3,NUMCO,NUMNO                      603430

C                                                                         B05080
      PARAMETER (NTMOL=36,NSPECI=85)                                      B05090
C                                                                         B05100
      COMMON /ISVECT/ ISOVEC(NTMOL),ISO82(NSPECI),ISONM(NTMOL),           B05110
     *                SMASSI(NSPECI)                                      B05120
      COMMON /LNC1/ RHOSLF(NSPECI),ALFD1(NSPECI),SCOR(NSPECI),ALFMAX,     B05130
     *              BETACR,DELTMP,DPTFC,DPTMN,XKT,NMINUS,NPLUS,NLIN,      B05140
     *              LINCNT,NCHNG,SUMALF,SUMZET,TRATIO,RHORAT,PAVP0,       B05150
     *              PAVP2,RECTLC,TMPDIF,ILC                               B05160
      DIMENSION MEFDP(64),FILHDR(2),AMOL(250)                             B05170
      DIMENSION A(4),B(4),TEMPLC(4)                                       B05180
      DIMENSION SABS(250),SRAD(250)
C                                                                         B05190
      EQUIVALENCE (MOL(1),AMOL(1)) , (S(1),SABS(1)), (EPP(1),SRAD(1))     B05200
      EQUIVALENCE (IHIRAC,FSCDID(1)) , (ILBLF4,FSCDID(2)),                B05210
     *            (IXSCNT,FSCDID(3)) , (IAERSL,FSCDID(4)),                B05220
     *            (JRAD,FSCDID(9)) , (XID(1),FILHDR(1))                   B05230
C                                                                         B05240
C     TEMPERATURES FOR LINE COUPLING COEFFICIENTS                         B05250
C                                                                         B05260
      DATA TEMPLC / 200.0,250.0,296.0,340.0 /                             B05270
      DATA HREJ /'0'/,HNOREJ /'1'/
      DATA NWDTH /0/
C                                                                         B05280
      NLNCR = NLNCR+1                                                     B05290
      IF (NLNCR.EQ.1) THEN                                                B05300
C                                                                         B05310
         XKT0 = TEMP0/RADCN2                                              B05320
         XKT = TAVE/RADCN2                                                B05330
         DELTMP = ABS(TAVE-TEMP0)                                         B05340
         BETACR = (1./XKT)-(1./XKT0)                                      B05350
         CALL MOLEC (2,SCOR,RHOSLF,ALFD1)                                 B05360
C                                                                         B05370
         TRATIO = TAVE/TEMP0                                              B05380
         RHORAT = (PAVE/P0)*(TEMP0/TAVE)                                  B05390
C                                                                         B05400
         PAVP0 = PAVE/P0                                                  B05410
         PAVP2 = PAVP0*PAVP0                                              B05420
C                                                                         B05430
C     FIND CORRECT TEMPERATURE AND INTERPOLATE FOR Y AND G                B05440
C                                                                         B05450
         DO 10 IL = 1, 3                                                  B05460
            ILC = IL                                                      B05470
            IF (TAVE.LT.TEMPLC(ILC+1)) GO TO 20                           B05480
   10    CONTINUE                                                         B05490
   20    IF (ILC.EQ.4) ILC = 3                                            B05500
C                                                                         B05510
         RECTLC = 1.0/(TEMPLC(ILC+1)-TEMPLC(ILC))                         B05520
         TMPDIF = TAVE-TEMPLC(ILC)                                        B05530
C                                                                         B05540
      ENDIF                                                               B05550
C
      IF (ILNFLG.EQ.2) READ(15)(FREJ(J),J=ILO,IHI)
C                                                                         B05560
      DO 30 J = ILO, IHI                                                  B05570
         YI = 0.                                                          B05580
         GI = 0.                                                          B05590
         GAMMA1 = 0.                                                      B05600
         GAMMA2 = 0.                                                      B05610
         I = IOUT(J)                                                      B05620
         IFLAG = IFLG(I)                                                  B05630
         MFULL=MOL(I)
         M = MOD(MOL(I),100)                                              B05640
C                                                                         B05650
C     ISO=(MOD(MOL(I),1000)-M)/100   IS PROGRAMMED AS:                    B05660
C                                                                         B05670
         ISO = MOD(MOL(I),1000)/100                                       B05680
         ILOC = ISOVEC(M)+ISO                                             B05690
C
         IF ((M.GT.NMOL).OR.(M.LT.1)) GO TO 25
C
         MOL(I) = M                                                       B05760
         SUI = S(I)*WK(M)                                                 B05770
C
         IF (SUI.EQ.0.) GO TO 25
C
         NLIN = NLIN+1                                                    B05830
C                                                                         B05840
C     Y'S AND G'S ARE STORED IN I+1 POSTION OF VNU,S,ALFA0,EPP...         B05850
C       A(1-4),  B(1-4) CORRESPOND TO TEMPERATURES TEMPLC(1-4) ABOVE      B05860
C                                                                         B05870
         IF (IFLAG.EQ.1.OR.IFLAG.EQ.3) THEN                               B05880
            A(1) = VNU(I+1)                                               B05890
            B(1) = S(I+1)                                                 B05900
            A(2) = ALFA0(I+1)                                             B05910
            B(2) = EPP(I+1)                                               B05920
            A(3) = AMOL(I+1)                                              B05930
            B(3) = HWHMS(I+1)                                             B05940
            A(4) = TMPALF(I+1)                                            B05950
            B(4) = PSHIFT(I+1)                                            B05960
C                                                                         B05970
C     CALCULATE SLOPE AND EVALUATE                                        B05980
C                                                                         B05990
            SLOPEA = (A(ILC+1)-A(ILC))*RECTLC                             B06000
            SLOPEB = (B(ILC+1)-B(ILC))*RECTLC                             B06010
C                                                                         B06020
            IF (IFLAG.EQ.1) THEN                                          B06030
               YI = A(ILC)+SLOPEA*TMPDIF                                  B06040
               GI = B(ILC)+SLOPEB*TMPDIF                                  B06050
            ELSE                                                          B06060
               GAMMA1 = A(ILC)+SLOPEA*TMPDIF                              B06070
               GAMMA2 = B(ILC)+SLOPEB*TMPDIF                              B06080
            ENDIF                                                         B06090
         ENDIF                                                            B06100
C                                                                         B06110
C     IFLAG = 2 IS RESERVED FOR LINE COUPLING COEFFICIENTS ASSOCIATED     B06120
C               WITH AN EXACT TREATMENT (NUMERICAL DIAGONALIZATION)       B06130
C                                                                         B06140
C     IFLAG = 3 TREATS LINE COUPLING IN TERMS OF REDUCED WIDTHS           B06150
C                                                                         B06160
         VNU(I) = VNU(I)+RHORAT*PSHIFT(I)                                 B06170
C                                                                         B06180
C     TEMPERATURE CORRECTION OF THE HALFWIDTH                             B06190
C     SELF TEMP DEPENDENCE TAKEN THE SAME AS FOREIGN                      B06200
C                                                                         B06210
         TMPCOR = TRATIO**TMPALF(I)                                       B06220
         ALFA0I = ALFA0(I)*TMPCOR                                         B06230
         HWHMSI = HWHMS(I)*TMPCOR                                         B06240
         ALFL = ALFA0I*(RHORAT-RHOSLF(ILOC))+HWHMSI*RHOSLF(ILOC)          B06250
C                                                                         B06260
         IF (IFLAG.EQ.3) ALFL = ALFL*(1.0-GAMMA1*PAVP0-GAMMA2*PAVP2)      B06270
C                                                                         B06280
         ALFAD = VNU(I)*ALFD1(ILOC)                                       B06290
         ZETA = ALFL/(ALFL+ALFAD)                                         B06300
         ZETAI(I) = ZETA                                                  B06310
         FZETA = 100.*ZETA
         IZ = FZETA + ONEPL                                               B06320
         IZETA(I) = IZ                                                    B06330
         ZETDIF = FZETA - FLOAT(IZ-1)
         ALFV = (AVRAT(IZ)+ZETDIF*(AVRAT(IZ+1)-AVRAT(IZ)))*(ALFL+ALFAD)   B06340
         IF (ALFV.LT.DV) THEN                                             B06350
            ALFV = DV                                                     B06360
            NMINAD = 1                                                    B06370
         ELSE                                                             B06380
            NMINAD = 0                                                    B06390
         ENDIF                                                            B06400
         IF (ALFV.GT.ALFMAX) THEN                                         B06410
            ALFV = ALFMAX                                                 B06420
            NPLSAD = 1                                                    B06430
         ELSE                                                             B06440
            NPLSAD = 0                                                    B06450
         ENDIF                                                            B06460
C
         IF (HWF3*ALFV+VNU(I) .LT. VFT) GO TO 25
C
         RECALF(I) = 1./ALFV                                              B06470
C                                                                         B06480
C     TREAT TRANSITIONS WITH UNKNOWN EPP AS SPECIAL CASE                  B06490
C                                                                         B06500
         IF (EPP(I).LT.0.) THEN                                           B06510
            IF (DELTMP.LE.10.) THEN                                       B06520
               EPP(I) = 0.                                                B06530
            ELSE                                                          B06540
               MEFDP(M) = MEFDP(M)+1                                      B06560
               GO TO 25
            ENDIF                                                         B06570
         ENDIF                                                            B06580
         IF (JRAD.NE.1) SUI = SUI*SCOR(ILOC)*                             B06590
     *                        EXP(-EPP(I)*BETACR)*(1.+EXP(-VNU(I)/XKT))   B06600
         IF (JRAD.EQ.1) SUI = SUI*SCOR(ILOC)*VNU(I)*                      B06610
     *                        EXP(-EPP(I)*BETACR)*(1.-EXP(-VNU(I)/XKT))   B06620
C                                                                         B06630
C                                                                         B06860
         SPI = SUI*(1.+GI*PAVP2)                                          B06870
         SPPI = SUI*YI*PAVP0                                              B06880
         SPPSP(I) = SPPI/SPI                                              B06890
         SABS(I)=SPI
         SRAD(I)=0.0
C
c ---from nlte:
         IF (MFULL.GE.1000) THEN
             FREQ=VNU(I)    
             NLOW=MOD(MFULL/1000,100)
             NUPP=MFULL/100000
             RLOW=1.0
             RUPP=1.0
             DELTA=EXP(-FREQ/XKT)
C                                 
C     PICK OUT MOLECULAR TYPE 
C                                                                               
C     H2O LINE                
C                                                                               
             IF (M.EQ.1) THEN
                 IF (NLOW.GT.NUMH2O) STOP 'NLOW GT NUMH2O IN LNCORQ'
                 IF (NLOW.GT.0) RLOW=RATH2O(NLOW)
                 IF (NUPP.GT.NUMH2O) STOP 'NUPP GT NUMH2O IN LNCORQ'
                 IF (NUPP.GT.0) RUPP=RATH2O(NUPP)
C                                                                       
C     CO2 LINE  
C                                                      
             ELSE IF (M.EQ.2) THEN
                 IF (NLOW.GT.NUMCO2) STOP 'NLOW GT NUMCO2 IN LNCORQ'
                 IF (NLOW.GT.0) RLOW=RATCO2(NLOW)    
                 IF (NUPP.GT.NUMCO2) STOP 'NUPP GT NUMCO2 IN LNCORQ'
                 IF (NUPP.GT.0) RUPP=RATCO2(NUPP)
C                               
C     O3 LINE 
C                          
             ELSE IF (M.EQ.3) THEN
                 IF (NLOW.GT.NUMO3) STOP 'NLOW GT NUMO3 IN LNCORQ'
                 IF (NLOW.GT.0) RLOW=RATO3(NLOW)
                 IF (NUPP.GT.NUMO3) STOP 'NUPP GT NUMO3 IN LNCORQ' 
                 IF (NUPP.GT.0) RUPP=RATO3(NUPP)
C                                       
C     CO LINE  
C                                               
             ELSE IF (M.EQ.5) THEN
                 IF (NLOW.GT.NUMCO ) STOP 'NLOW GT NUMCO  IN LNCORQ' 
                 IF (NLOW.GT.0) RLOW=RATCO (NLOW)
                 IF (NUPP.GT.NUMCO ) STOP 'NUPP GT NUMCO  IN LNCORQ'
                 IF (NUPP.GT.0) RUPP=RATCO (NUPP) 
C                        
C     NO LINE 
C                
             ELSE IF (M.EQ.8) THEN
                 IF (NLOW.GT.NUMNO) STOP 'NLOW GT NUMNO IN LNCORQ'
                 IF (NLOW.GT.0) RLOW=RATNO(NLOW)
                 IF (NUPP.GT.NUMNO) STOP 'NUPP GT NUMNO IN LNCORQ'
                 IF (NUPP.GT.0) RUPP=RATNO(NUPP)
             ELSE
                 PRINT 900,M  
  900            FORMAT('  MOL  IN TROUBLE',I10) 
                 SABS(I)=0.                                  
                 SRAD(I)=0.  
                 SPPSP(I)=0. 
                 GO TO 30
             END IF 
C                                                                        
C     RLOW AND RUPP NOW SET 
C                                                                         
             FNLTE=SPI/(1.0-DELTA)   
             SABS(I)=FNLTE*(RLOW-RUPP*DELTA) 
             SRAD(I)=FNLTE*(RLOW-RUPP) 
C                                 
             IF (IFLAG .EQ. 0) THEN                                   
                 SPEAK = SABS(I)*ABS(RECALF(I)) 
                 SLFABS = SPEAK         
                 IF(SPEAK.GE.5.)  THEN
                     SLFABS = 1. 
                 ELSE   
                     IF(SPEAK.GT.0.01) SLFABS = 1.-EXP(-SPEAK)
                 ENDIF
                 TEST = SLFABS *(1.-SRAD(I)/SABS(I)) 
                 IF(TEST.LE.DPTMN) GOTO 25                              
             END IF                                           

c --- above from nlte
         ELSE
C                                                               
             IF (IFLAG.EQ.0) THEN
                 IF (ILNFLG.LE.1) THEN
                     FREJ(J) = HNOREJ
                     SPEAK = SUI*RECALF(I)                                B06650
                     IF (DVR4.LE.0.) THEN     
                         IF (SPEAK.LE.DPTMN) THEN 
                             FREJ(J) = HREJ
                             GO TO 25
                         ENDIF
                     ELSE                                                 B06680
                         JJ = (VNU(I)-V1R4)/DVR4+1.                       B06690
                         JJ = MAX(JJ,1)                                   B06700
                         JJ = MIN(JJ,NPTR4)                               B06710
c                         IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ))) THEN
                         IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ)) 
     &                       .and. sppsp(i).eq.0.) THEN
                             FREJ(J) = HREJ
                             GO TO 25
                         ENDIF
                     ENDIF                                                B06730
                 ELSE
C     "ELSE" IS TRUE WHEN "ILNFLG" EQUALS 2
C
                     IF (FREJ(J).EQ.HREJ) GO TO 25
                 ENDIF
             ENDIF                                                        B06790
C                                                                         B06800
         ENDIF    

         NMINUS = NMINUS+NMINAD                                           B06810
         NPLUS = NPLUS+NPLSAD                                             B06820
         SUMALF = SUMALF+ALFV                                             B06830
         SUMZET = SUMZET+ZETA                                             B06840
         LINCNT = LINCNT+1                                                B06850

         GO TO 30
C
  25     SABS(I)=0.0
         SRAD(I)=0.0
         SPPSP(I) = 0.0  
C                                                                         B06900
  30  CONTINUE                                                            B06910
C                                                                         B06920
      NCHNG = NMINUS+NPLUS                                                B06930
      IF (ILNFLG.EQ.1) WRITE(15)(FREJ(J),J=ILO,IHI)
C                                                                         B06940
      RETURN                                                              B06950
C                                                                         B06960
      END                                                                 B06970
      SUBROUTINE CNVFNQ (VNU,SABS,SRAD,SPPSP,RECALF,R1,R2,R3,RR1,
     &    RR2,RR3,ZETAI,IZETA)
C                                                                         B07000
      IMPLICIT REAL*8           (V)                                       B07010
C                                                                         B07020
C     SUBROUTINE CNVFNV PERFORMS THE CONVOLUTION OF THE LINE DATA WITH    B07030
C     THE VOIGT LINE SHAPE (APPROXIMATED)                                 B07040
C                                                                         B07050
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  B07060
C                                                                         B07070
C                                                                         B07110
C     IMPLEMENTATION:    R.D. WORSHAM                                     B07120
C                                                                         B07130
C     ALGORITHM REVISIONS:    S.A. CLOUGH                                 B07140
C     R.D. WORSHAM                                                        B07150
C     J.L. MONCET                                                         B07160
C                                                                         B07170
C                                                                         B07180
C     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.                         B07190
C     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139                          B07200
C                                                                         B07210
C     ------------------------------------------------------------------  B07220
C                                                                         B07230
C     WORK SUPPORTED BY:    THE ARM PROGRAM                               B07240
C     OFFICE OF ENERGY RESEARCH                                           B07250
C     DEPARTMENT OF ENERGY                                                B07260
C                                                                         B07270
C                                                                         B07280
C     SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL              B07290
C                                                                         B07300
C     FASCOD3                                                             B07310
C                                                                         B07320
C                                                                         B07330
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  B07340
C                                                                         B07350
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ
C                                                                         B07370
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       B07380
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   B07390
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    B07400
      COMMON /XSUB/ VBOT,VTOP,VFT,LIMIN,ILO,IHI,IEOF,IPANEL,ISTOP,IDATA   B07410
      COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTS       B07420
      COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,            B07430
     *                HWF3,DXF3,NX3,N3MAX                                 B07440
      COMMON /SUB1/ MAX1,MAX2,MAX3,NLIM1,NLIM2,NLIM3,NLO,NHI,DVR2,DVR3,   B07450
     *              N1R1,N2R1,N1R2,N2R2,N1R3,N2R3                         B07460
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,         B07470
     *               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL                      B07480
      COMMON /VOICOM/ AVRAT(102),CGAUSS(102),CF1(102),CF2(102),           B07490
     *                CF3(102),CER(102)                                   B07500
      COMMON /IOU/ IOUT(250)                                              B07510
C                                                                         B07520
      PARAMETER (NFPTS=2001,NFMX=1.3*NFPTS, NUMZ = 101)  
      COMMON /FNSHQ/ IFN,F1(NFMX, NUMZ),F2(NFMX, NUMZ),F3(NFMX, NUMZ),
     $     FG(NFMX)
c
      DIMENSION VNU(*),SABS(*),SRAD(*),SPPSP(*),RECALF(*)                 B07530
      DIMENSION R1(*),R2(*),R3(*)                                         B07540
      DIMENSION RR1(*),RR2(*),RR3(*)
      DIMENSION IZETA(*),ZETAI(*)                                         B07570
C
      equivalence (zetdif,A)
C                                                                         B07580
      CALL CPUTIM (TIME0)                                                 B07590
C                                                                         B07600
      CLC1 = 4./(FLOAT(NX1-1))                                            B07610
      CLC2 = 16./(FLOAT(NX2-1))                                           B07620
      CLC3 = 64./(FLOAT(NX3-1))                                           B07630
      WAVDXF = DV/DXF1                                                    B07640
      HWDXF = HWF1/DXF1                                                   B07650
      CONF2 = DV/DVR2                                                     B07660
      CONF3 = DV/DVR3                                                     B07670
      ILAST = ILO-1                                                       B07680
C                                                                         B07690
      IF (ILO.LE.IHI) THEN                                                B07700
         DO 30 J = ILO, IHI                                               B07710
            I = IOUT(J)                                                   B07720
            IF (SABS(I).NE.0.) THEN                                       B07730
               DEPTHA = SABS(I)*RECALF(I)                                 B07740
               DEPTHR = SRAD(I)*RECALF(I)
               IZM = IZETA(I)                                             B07750
               ZETDIF = 100.*ZETAI(I)-FLOAT(IZM-1)                        B07840
C                                                                         B07930
               ZSLOPE = RECALF(I)*WAVDXF                                  B07940
               ZINT = (VNU(I)-VFT)/DV                                     B07950
               BHWDXF = HWDXF/ZSLOPE                                      B07960
               JMAX1 = ZINT+BHWDXF+1.5                                    B07970
               IF (JMAX1.GT.MAX1) THEN                                    B07980
                  ILAST = J-1                                             B07990
                  IPANEL = 1                                              B08000
                  GO TO 40                                                B08010
               ENDIF                                                      B08020
               JMIN1 = ZINT-BHWDXF+1.5                                    B08030
               RSHFT = 0.5                                                B08040
               IF (ZINT.LT.0.0) RSHFT = -RSHFT                            B08050
               J2SHFT = ZINT*(1.-CONF2)+RSHFT                             B08060
               J3SHFT = ZINT*(1.-CONF3)+RSHFT                             B08070
               JMIN2 = JMIN1-J2SHFT                                       B08080
               JMIN3 = JMIN1-J3SHFT                                       B08090
               ZF1 = (FLOAT(JMIN1-2)-ZINT)*ZSLOPE                         B08100
               ZF2 = (FLOAT(JMIN2-2)-ZINT*CONF2)*ZSLOPE                   B08110
               ZF3 = (FLOAT(JMIN3-2)-ZINT*CONF3)*ZSLOPE                   B08120
c
               IF (SPPSP(I).eq.0.) then
c
               DO 10 J1 = JMIN1, JMAX1                                    B08160
                  J2 = J1-J2SHFT                                          B08170
                  J3 = J1-J3SHFT                                          B08180
                  ZF3 = ZF3+ZSLOPE                                        B08190
                  ZF2 = ZF2+ZSLOPE                                        B08200
                  ZF1 = ZF1+ZSLOPE                                        B08210
                  IZ3 = ABS(ZF3)+1.5                                      B08220
                  IZ2 = ABS(ZF2)+1.5                                      B08230
                  IZ1 = ABS(ZF1)+1.5                                      B08240
C     ************Using new voigt scheme: 
C     ************Interpolate voigt subfunctions to zeta 

c          A  is equivalenced to ZETDIF

                     x3 = deptha*( (1.-A)*F3(IZ3,IZM)+A*F3(IZ3,IZM+1))
                     x2 = deptha*(( 1.-A)*F2(IZ2,IZM)+A*F2(IZ2,IZM+1))
                     x1 = deptha*( (1.-A)*F1(IZ1,IZM)+A*F1(IZ1,IZM+1))

                     R3(J3) = R3(J3)+ x3
                     R2(J2) = R2(J2)+ x2
                     R1(J1) = R1(J1)+ x1

                  IF (DEPTHR.NE.0) THEN

                     xx3 = depthr*( (1.-A)*F3(IZ3,IZM)+A*F3(IZ3,IZM+1))
                     xx2 = depthr*(( 1.-A)*F2(IZ2,IZM)+A*F2(IZ2,IZM+1))
                     xx1 = depthr*( (1.-A)*F1(IZ1,IZM)+A*F1(IZ1,IZM+1))

                     RR3(J3) = RR3(J3)+ xx3
                     RR2(J2) = RR2(J2)+ xx2
                     RR1(J1) = RR1(J1)+ xx1

                  ENDIF
   10          CONTINUE                                                   B08290
C                                                                         B08300
           else
C                                                                         B08320
C                 THE FOLLOWING DOES LINE COUPLING                        B08330
C                                                                         B08340
C                 SPPSP(I) = SPP(I)/SP(I)                                 B08350
C                                                                         B08360

               dptrat1 = SPPSP(I)*clc1
               dptrat2 = SPPSP(I)*clc2
               dptrat3 = SPPSP(I)*clc3
C                                                                         B08430
C                                                                         B08440
                  DO 20 J1 = JMIN1, JMAX1                                 B08450
                     J2 = J1-J2SHFT                                       B08460
                     J3 = J1-J3SHFT                                       B08470
                     ZF3 = ZF3+ZSLOPE                                    B08480
                     ZF2 = ZF2+ZSLOPE                                    B08490
                     ZF1 = ZF1+ZSLOPE                                    B08500
                     IZ3 = ABS(ZF3)+1.5                                  B08510
                     IZ2 = ABS(ZF2)+1.5                                  B08520
                     IZ1 = ABS(ZF1)+1.5                                  B08530
c                     
C     ************Using new voigt scheme: 
C     ************Interpolate voigt subfunctions to zeta 

                  x3 = DEPTHa*( (1.-A)*F3(IZ3,IZM)+A*F3(IZ3,IZM+1))
                  x2 = DEPTHa*( (1.-A)*F2(IZ2,IZM)+A*F2(IZ2,IZM+1))
                  x1 = DEPTHa*( (1.-A)*F1(IZ1,IZM)+A*F1(IZ1,IZM+1))

                  y3 = dptrat3*x3*ZF3
                  y2 = dptrat2*x2*ZF2
                  y1 = dptrat1*x1*ZF1

                  R3(J3) = R3(J3) + x3 + y3
                  R2(J2) = R2(J2) + x2 + y2
                  R1(J1) = R1(J1) + x1 + y1

                     IF (DEPTHR.NE.0) THEN
                  xx3 = DEPTHr*( (1.-A)*F3(IZ3,IZM)+A*F3(IZ3,IZM+1))
                  xx2 = DEPTHr*( (1.-A)*F2(IZ2,IZM)+A*F2(IZ2,IZM+1))
                  xx1 = DEPTHr*( (1.-A)*F1(IZ1,IZM)+A*F1(IZ1,IZM+1))

                  yy3 = dptrat3*xx3*ZF3
                  yy2 = dptrat2*xx2*ZF2
                  yy1 = dptrat1*xx1*ZF1

                  RR3(J3) = RR3(J3) + xx3 + yy3
                  RR2(J2) = RR2(J2) + xx2 + yy2
                  RR1(J1) = RR1(J1) + xx1 + yy1

                     ENDIF

   20             CONTINUE                                                B08580
C                                                                         B08590
               ENDIF                                                      B08600
            ENDIF                                                         B08610
C                                                                         B08620
   30    CONTINUE                                                         B08630
         ILAST = IHI                                                      B08640
C                                                                         B08650
C        IDATA=0 FOR MORE DATA REQUIRED                                   B08660
C        IDATA=1 IF NO MORE DATA REQUIRED                                 B08670
C                                                                         B08680
         IPANEL = IDATA                                                   B08690
      ELSE                                                                B08700
         IPANEL = 1                                                       B08710
      ENDIF                                                               B08720
C                                                                         B08730
   40 ILO = ILAST+1                                                       B08740
      CALL CPUTIM (TIME)                                                  B08750
      TIMCNV = TIMCNV+TIME-TIME0                                          B08760
      RETURN                                                              B08770
C                                                                         B08780
      END                                                                 B08790
      SUBROUTINE PANELQ (R1,R2,R3,RR1,RR2,RR3,KFILE,JRAD,IENTER)          B09990
C                                                                         B10000
      IMPLICIT REAL*8           (V)                                     ! B10010
C                                                                         B10020
C     SUBROUTINE PANEL COMBINES RESULTS OF R3, R2, AND R1 INTO R1 ARRAY   B10030
C     AND OUTPUTS THE ARRAY R1                                            B10040
C                                                                         B10050
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B10060
C                                                                         B10070
C               LAST MODIFICATION:    28 AUGUST 1992                      B10080
C                                                                         B10090
C                  IMPLEMENTATION:    R.D. WORSHAM                        B10100
C                                                                         B10110
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         B10120
C                                     R.D. WORSHAM                        B10130
C                                     J.L. MONCET                         B10140
C                                                                         B10150
C                                                                         B10160
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B10170
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          B10180
C                                                                         B10190
C----------------------------------------------------------------------   B10200
C                                                                         B10210
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     B10220
C                                     OFFICE OF ENERGY RESEARCH           B10230
C                                     DEPARTMENT OF ENERGY                B10240
C                                                                         B10250
C                                                                         B10260
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             B10270
C                                                                         B10280
C                                             FASCOD3                     B10290
C                                                                         B10300
C                                                                         B10310
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B10320
C                                                                         B10330
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ
C                                                                         B10350
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       B10360
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   B10370
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    B10380
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /XSUB/ VBOT,VTOP,VFT,LIMIN,ILO,IHI,IEOF,IPANEL,ISTOP,IDATA   B10400
      COMMON /SUB1/ MAX1,MAX2,MAX3,NLIM1,NLIM2,NLIM3,NLO,NHI,DVR2,DVR3,   B10410
     *              N1R1,N2R1,N1R2,N2R2,N1R3,N2R3                         B10420
      COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTS       B10430
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,         B10440
     *               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL                      B10450
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         B10460
     *              NLNGTH,KDUMM,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        B10470
     *              NLTEFL,LNFIL4,LNGTH4                                  B10480
      COMMON /IODFLG/ DVOUT
      DIMENSION R1(*),R2(*),R3(*)                                         B10490
      DIMENSION RR1(*),RR2(*),RR3(*)                                      B10490
      DIMENSION PNLHDR(2)                                                 B10500
C                                                                         B10510
      EQUIVALENCE (V1P,PNLHDR(1))                                         B10520
C                                                                         B10530
      CALL CPUTIM (TIME0)                                                 B10540
      X00 = -7./128.                                                      B10550
      X01 = 105./128.                                                     B10560
      X02 = 35./128.                                                      B10570
      X03 = -5./128.                                                      B10580
      X10 = -1./16.                                                       B10590
      X11 = 9./16.                                                        B10600
      ISTOP = 0                                                           B10610
C
C     Test for last panel.  If last, set the last point to one point
C     greater than V1 specified on TAPE5 (to ensure last point for
C     every layer is the same)
C
      IF ((VFT+(NHI-1)*DVP).GT.V2) THEN
         NHI = (V2-VFT)/DVP + 1.
         V2P = VFT+FLOAT(NHI-1)*DVP
         IF (V2P.LT.V2) THEN
            V2P = V2P+DVP
            NHI = NHI+1
         ENDIF
         ISTOP = 1                                                        B10640
      ELSE
         V2P = VFT+FLOAT(NHI-1)*DV
      ENDIF
      NLIM = NHI-NLO+1
      V1P = VFT+FLOAT(NLO-1)*DV
C
      LIMLO = N1R2                                                        B10670
      IF (N1R2.EQ.1) LIMLO = LIMLO+4                                      B10680
      LIMHI = (NHI/4)+1                                                   B10690
C                                                                         B10700
      DO 10 J = LIMLO, LIMHI, 4                                           B10710
         J3 = (J-1)/4+1                                                   B10720
         R2(J) = R2(J)+R3(J3)                                             B10730
         R2(J+1) = R2(J+1)+X00*R3(J3-1)+X01*R3(J3)+X02*R3(J3+1)+          B10740
     *             X03*R3(J3+2)                                           B10750
         R2(J+2) = R2(J+2)+X10*(R3(J3-1)+R3(J3+2))+                       B10760
     *             X11*(R3(J3)+R3(J3+1))                                  B10770
         R2(J+3) = R2(J+3)+X03*R3(J3-1)+X02*R3(J3)+X01*R3(J3+1)+          B10780
     *             X00*R3(J3+2)                                           B10790
         RR2(J) = RR2(J)+RR3(J3)                                          B10730
         RR2(J+1) = RR2(J+1)+X00*RR3(J3-1)+X01*RR3(J3)+X02*RR3(J3+1)+     B10740
     *             X03*RR3(J3+2)                                          B10750
         RR2(J+2) = RR2(J+2)+X10*(RR3(J3-1)+RR3(J3+2))+                   B10760
     *             X11*(RR3(J3)+RR3(J3+1))                                B10770
         RR2(J+3) = RR2(J+3)+X03*RR3(J3-1)+X02*RR3(J3)+X01*RR3(J3+1)+     B10780
     *             X00*RR3(J3+2)                                          B10790
   10 CONTINUE                                                            B10800
      DO 20 J = NLO, NHI, 4                                               B10810
         J2 = (J-1)/4+1                                                   B10820
         R1(J) = R1(J)+R2(J2)                                             B10830
         R1(J+1) = R1(J+1)+X00*R2(J2-1)+X01*R2(J2)+X02*R2(J2+1)+          B10840
     *             X03*R2(J2+2)                                           B10850
         R1(J+2) = R1(J+2)+X10*(R2(J2-1)+R2(J2+2))+                       B10860
     *             X11*(R2(J2)+R2(J2+1))                                  B10870
         R1(J+3) = R1(J+3)+X03*R2(J2-1)+X02*R2(J2)+X01*R2(J2+1)+          B10880
     *             X00*R2(J2+2)                                           B10890
         RR1(J) = RR1(J)+RR2(J2)                                          B10830
         RR1(J+1) = RR1(J+1)+X00*RR2(J2-1)+X01*RR2(J2)+X02*RR2(J2+1)+     B10840
     *             X03*RR2(J2+2)                                          B10850
         RR1(J+2) = RR1(J+2)+X10*(RR2(J2-1)+RR2(J2+2))+                   B10860
     *             X11*(RR2(J2)+RR2(J2+1))                                B10870
         RR1(J+3) = RR1(J+3)+X03*RR2(J2-1)+X02*RR2(J2)+X01*RR2(J2+1)+     B10880
     *             X00*RR2(J2+2)                                          B10890
   20 CONTINUE                                                            B10900
C                                                                         B10910
C     IN THE FOLLOWING SECTION THE ABSORPTION COEFFICIENT IS MULTIPIIED   B10960
C     BY THE RADIATION FIELD                                              B10970
C                                                                         B10980
      IF (JRAD.EQ.0) THEN                                                 B10990
C                                                                         B11000
         XKT = TAVE/RADCN2                                                B11010
         VI = V1P-DV                                                      B11020
         VITST = VI                                                       B11030
         RDLAST = -1.                                                     B11040
         NPTSI1 = NLO-1                                                   B11050
         NPTSI2 = NLO-1                                                   B11060
C                                                                         B11070
   30    NPTSI1 = NPTSI2+1                                                B11080
C                                                                         B11090
         VI = VFT+FLOAT(NPTSI1-1)*DV                                      B11100
         RADVI = RADFNI(VI,DV,XKT,VITST,RDEL,RDLAST)                      B11110
C                                                                         B11130
         NPTSI2 = (VITST-VFT)/DV+1.001                                    B11140
         NPTSI2 = MIN(NPTSI2,NHI)                                         B11150
C                                                                         B11160
         DO 40 I = NPTSI1, NPTSI2                                         B11170
C           VI = VI+DV                                                    B11180
            R1(I) = R1(I)*RADVI                                           B11190
            RADVI = RADVI+RDEL                                            B11200
   40    CONTINUE                                                         B11210
C                                                                         B11220
         IF (NPTSI2.LT.NHI) GO TO 30                                      B11230
C                                                                         B11240
      ENDIF                                                               B11250
C                                                                         B11260
C     V1P IS FIRST FREQ OF PANEL                                          B11270
C     V2P IS LAST FREQ OF PANEL                                           B11280
C                                                                         B11290
C     If DVOUT (carried in from COMMON BLOCK /IODFLG/) is zero,
C     then no interpolation is necessary.  For nozero DVOUT
C     (in case of IOD>0 and IMRG=1), call PNLINT.
C
      IF (DVOUT.EQ.0.) THEN                                               B11300
         CALL BUFOUT (KFILE,PNLHDR(1),NPHDRF)                             B11310
         CALL BUFOUT (KFILE,R1(NLO),NLIM)                                 B11320
         CALL BUFOUT (KFILE,RR1(NLO),NLIM)                                B11320
C                                                                         B11330
         IF (NPTS.GT.0) CALL R1PRNT (V1P,DVP,NLIM,R1,NLO,NPTS,KFILE,
     *                               IENTER)                              B11340
         IF (NPTS.GT.0) CALL R1PRNT (V1P,DVP,NLIM,RR1,NLO,NPTS,KFILE,
     *                               IENTER)                              B11340
      ELSE                                                                B11350
         CALL PNLINT (R1(NLO),IENTER)                                     B11360
         CALL PNLINT (RR1(NLO),IENTER)                                    B11360
      ENDIF                                                               B11370
C                                                                         B11380
      VFT = VFT+FLOAT(NLIM1-1)*DV                                         B11390
      IF (ISTOP.NE.1) THEN                                                B11400
         DO 50 J = NLIM1, MAX1                                            B11420
            R1(J-NLIM1+1) = R1(J)                                         B11430
            RR1(J-NLIM1+1) = RR1(J)                                       B11430
   50    CONTINUE                                                         B11450
         DO 60 J = MAX1-NLIM1+2, MAX1                                     B11460
            R1(J) = 0.                                                    B11470
            RR1(J) = 0.                                                   B11470
   60    CONTINUE                                                         B11480
         DO 70 J = NLIM2, MAX2                                            B11500
            R2(J-NLIM2+1) = R2(J)                                         B11510
            RR2(J-NLIM2+1) = RR2(J)                                       B11510
   70    CONTINUE                                                         B11530
         DO 80 J = MAX2-NLIM2+2, MAX2                                     B11540
            R2(J) = 0.                                                    B11550
            RR2(J) = 0.                                                   B11550
   80    CONTINUE                                                         B11560
         DO 90 J = NLIM3, MAX3                                            B11580
            R3(J-NLIM3+1) = R3(J)                                         B11590
            RR3(J-NLIM3+1) = RR3(J)                                       B11590
   90    CONTINUE                                                         B11610
         DO 100 J = MAX3-NLIM3+2, MAX3                                    B11620
            R3(J) = 0.                                                    B11630
            RR3(J) = 0.                                                   B11630
  100    CONTINUE                                                         B11640
         NLO = NSHIFT+1                                                   B11650
      ENDIF                                                               B11660
      CALL CPUTIM (TIME)                                                  B11670
      TIMPNL = TIMPNL+TIME-TIME0                                          B11680
C                                                                         B11690
      RETURN                                                              B11700
C                                                                         B11710
      END                                                                 B11720
      SUBROUTINE LBLFQ (JRAD,V1,V2)                                       D04440
C                                                                         D04450
      IMPLICIT REAL*8           (V)                                     ! D04460
C                                                                         D04470
C     SUBROUTINE LBLF4 DOES A LINE BY LINE CALCULATION                    D04480
C     USING FUNCTION F4.                                                  D04490
C                                                                         D04500
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           D04510
      COMMON /BUF/ VNU(1250),SABS(1250),ALFAL(1250),ALFAD(1250),
     &    MOL(1250),SPP(1250),SRAD(1250)                                  D04530
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   D04540
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       D04550
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       D04560
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      D04570
     *              EXTID(10)                                             D04580
C                                                                         D04590
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SEC   ,       XALTZ
C                                                                         D04610
      COMMON /FILHDR/ XID(10),SEC   ,PAVE,TAVE,HMOLID(60),XALTZ(4),       D04620
     *                W(60),PZL,PZU,TZL,TZU,WBROAD,DVO,V1H,V2H,TBOUND,    D04630
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    D04640
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,         D04660
     *               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL                      D04670
      COMMON /R4SUB/ VLO,VHI,ILO,IST,IHI,LIMIN,LIMOUT,ILAST,DPTMN,        D04680
     *               DPTFC,ILIN4,ILIN4T                                   D04690
      COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)        D04700
      COMMON /VOICOM/ AVRAT(102),DUMMY(5,102)                             D04710
      COMMON /CONVF/ CHI(251),RDVCHI,RECPI,ZSQBND,A3,B3,JCNVF4            D04720
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         D04730
     *              NLNGTH,KFILE,KPANEL,LINFIO,NFILE,IAFIL,IEXFIL,        D04740
     *              NLTEFL,LNFIL4,LNGTH4                                  D04750
C                                                                         D04760
      EQUIVALENCE (IHIRAC,FSCDID(1)) , (ILBLF4,FSCDID(2))                 D04770
C                                                                         D04780
      DATA JLBLF4 / 0 /                                                   D04790
C                                                                         D04800
      CALL CPUTIM (TIM00)                                                 D04810
C                                                                         D04820
      DPTMN = DPTMIN/RADFN(V2,TAVE/RADCN2)                                D04830
      DPTFC = DPTFAC                                                      D04840
      LIMIN = 1000                                                        D04850
      LIMOUT = 2500                                                       D04860
      JLBLF4 = 1                                                          D04870
C                                                                         D04880
C     SET IEOF EQUAL TO -1 FOR FIRST READ                                 D04890
C                                                                         D04900
      IEOF = -1                                                           D04910
C                                                                         D04920
      V1R4 = V1                                                           D04930
      V2R4 = V2                                                           D04940
      NPTR4 = (V2R4-V1R4)/DVR4+ONEPL                                      D04950
      NPTR4 = MIN(NPTR4,LIMOUT)                                           D04960
      V2R4 = V1R4+DVR4*FLOAT(NPTR4-1)                                     D04970
C                                                                         D04980
      LIMP2 = LIMOUT+2                                                    D04990
      DO 10 I = 1, LIMP2                                                  D05000
         R4(I) = 0.                                                       D05010
         RR4(I) = 0.                                                      D05010
   10 CONTINUE                                                            D05020
      BETA = RADCN2/TAVE                                                  D05030
      VLO = V1R4-BOUND4                                                   D05040
      VHI = V2R4+BOUND4                                                   D05050
   20 CALL CPUTIM (TIM0)                                                  D05060
      CALL RDLN4Q (IEOF)                                                  D05070
      CALL CPUTIM (TIM1)                                                  D05080
C                                                                         D05090
      IF (IEOF.EQ.2) THEN                                                 D05100
         TF4 = TF4+TIM1-TIM00                                             D05110
         RETURN                                                           D05120
      ENDIF                                                               D05130
C                                                                         D05140
      TF4RDF = TF4RDF+TIM1-TIM0                                           D05150
      TIM2 = TIM1                                                         D05160
      IF (IEOF.EQ.1.AND.IHI.EQ.0) GO TO 30                                D05170
C                                                                         D05180
      CALL CNVF4Q (VNU,SABS,ALFAL,ALFAD,MOL,SPP,SRAD)                     D05190
C                                                                         D05200
      CALL CPUTIM (TIM3)                                                  D05210
      TF4CNV = TF4CNV+TIM3-TIM2                                           D05220
C                                                                         D05230
C    IF IHI EQUALS -1 THEN END OF CONVOLUTION                             D05240
C                                                                         D05250
      IF (IHI.EQ.-1) GO TO 30                                             D05260
      GO TO 20                                                            D05270
C                                                                         D05280
   30 CALL CPUTIM (TIM4)                                                  D05290
C                                                                         D05300
      IF (JRAD.EQ.1) THEN                                                 D05310
C                                                                         D05320
C     RADIATION FIELD                                                     D05330
C                                                                         D05340
         XKT = 1./BETA                                                    D05350
         VITST = V1R4-DVR4                                                D05360
         RDLAST = -1.                                                     D05370
         NPTSI1 = 0                                                       D05380
         NPTSI2 = 0                                                       D05390
C                                                                         D05400
   40    NPTSI1 = NPTSI2+1                                                D05410
C                                                                         D05420
         VI = V1R4+DVR4*FLOAT(NPTSI1-1)                                   D05430
         RADVI = RADFNI(VI,DVR4,XKT,VITST,RDEL,RDLAST)                    D05440
C                                                                         D05460
         NPTSI2 = (VITST-V1R4)/DVR4+1.001                                 D05470
         NPTSI2 = MIN(NPTSI2,NPTR4)                                       D05480
C                                                                         D05490
         DO 50 I = NPTSI1, NPTSI2                                         D05500
            VI = VI+DVR4                                                  D05510
            R4(I) = R4(I)*RADVI                                           D05520
            RR4(I) = RR4(I)*RADVI                                         D05530
            RADVI = RADVI+RDEL                                            D05520
   50    CONTINUE                                                         D05540
C                                                                         D05550
         IF (NPTSI2.LT.NPTR4) GO TO 40                                    D05560
      ENDIF                                                               D05570
C                                                                         D05580
      CALL CPUTIM (TIM5)                                                  D05590
      TF4PNL = TF4PNL+TIM5-TIM4                                           D05600
      TF4 = TF4+TIM5-TIM00                                                D05610
C                                                                         D05620
      RETURN                                                              D05630
C                                                                         D05640
      END                                                                 D05650
      SUBROUTINE RDLN4Q (IEOF)                                            D05660
C                                                                         D05670
      IMPLICIT REAL*8           (V)                                     ! D02670
C                                                                         D05690
C     SUBROUTINE RDLIN4Q INPUTS THE LINE DATA FROM LNFIL4                 D05700
C                                                                         D05710
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SEC   ,       XALTZ
C                                                                         D05730
      COMMON /FILHDR/ XID(10),SEC   ,PAV ,TAV ,HMOLID(60),XALTZ(4),       D05740
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   D05750
     *                EMIS ,FSCDID(17),NMOL,LAYRS ,YID1,YID(10),LSTWDF    D05760
      COMMON /R4SUB/ VLO,VHI,ILO,IST,IHI,LIMIN,LIMOUT,ILAST,DPTMN,        D05770
     *               DPTFC,ILIN4,ILIN4T                                   D05780
      COMMON /BUF/ VNU(1250),S(1250),ALFAL(1250),ALFAD(1250),MOL(1250),   D05790
     *             SPP(1250)                                              D05800
      COMMON /BUF2/ VMIN,VMAX,NREC,NWDS                                   D05810
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         D05820
     *              NLNGTH,KFILE,KPANEL,LINDUM,NFILE,IAFIL,IEXFIL,        D05830
     *              NLTEFL,LNFIL4,LNGTH4                                  D05840
      common /eppinfo/ negepp_flag
      integer*4 negepp_flag

      DIMENSION DUM(2),LINPNL(2)                                          D05850
C                                                                         D05860
      EQUIVALENCE (VMIN,LINPNL(1))                                        D05870
C                                                                         D05880
      IF (IEOF.EQ.-1) THEN                                                D05890
C                                                                         D05900
C     BUFFER PAST FILE HEADER                                             D05910
C                                                                         D05920
         REWIND LNFIL4                                                    D05930
         ILIN4T = 0                                                       D05940
         CALL BUFIN (LNFIL4,LEOF,DUM(1),1)                                D05950
         IF (LEOF.EQ.0) STOP 'RDLIN4; TAPE9 EMPTY'                        D05960
         IF (LEOF.EQ.-99) THEN                                            D05970
            IEOF = 2                                                      D05980
C                                                                         D05990
            RETURN                                                        D06000
C                                                                         D06010
         ENDIF                                                            D06020
      ENDIF                                                               D06030
      IEOF = 0                                                            D06040
      ILO = 1                                                             D06050
      IHI = 0                                                             D06060
C                                                                         D06070
   10 CALL BUFIN (LNFIL4,LEOF,LINPNL(1),NPHDRL)                           D06080
      IF (LEOF.EQ.0) GO TO 20                                             D06090
      ILIN4T = ILIN4T+NREC                                                D06100
      IF (VMAX.LT.VLO) THEN                                               D06110
         CALL BUFIN (LNFIL4,LEOF,DUM(1),1)                                D06120
         GO TO 10                                                         D06130
      ELSE                                                                D06140
         CALL BUFIN (LNFIL4,LEOF,VNU(1),NWDS)                             D06150
      ENDIF                                                               D06160
      IHI = NREC                                                          D06170
      IF (VNU(NREC).GT.VHI) GO TO 20                                      D06180
C                                                                         D06190
      RETURN                                                              D06200
C                                                                         D06210
   20 IEOF = 1                                                            D06220
C                                                                         D06230
      RETURN                                                              D06240
C                                                                         D06250
      END                                                                 D06260


      SUBROUTINE CNVF4Q (VNU,SABS,ALFAL,ALFAD,MOL,SPP,SRAD)               D06270
C                                                                         D06280
      IMPLICIT REAL*8           (V)                                     ! D06290
C                                                                         D06300
C     SUBROUTINE CONVF4 CONVOLVES THE LINE DATA WITH FUNCTION F4          D06310
C                                                                         D06320
      CHARACTER*1 FREJ(1250),HREJ,HNOREJ
      COMMON /RCNTRL/ ILNFLG
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           D06330
      COMMON /R4SUB/ VLO,VHI,ILO,IST,IHI,LIMIN,LIMOUT,ILAST,DPTMN,        D06340
     *               DPTFC,ILIN4,ILIN4T                                   D06350
      COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)        D06360
      COMMON /VOICOM/ AVRAT(102),DUMMY(5,102)                             D06370
      COMMON /CONVF/ CHI(251),RDVCHI,RECPI,ZSQBND,A3,B3,JCNVF4            D06380
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         D05820
     *              NLNGTH,KFILE,KPANEL,LINDUM,NFILE,IAFIL,IEXFIL,        D05830
     *              NLTEFL,LNFIL4,LNGTH4                                  D05840
C                                                                         D06390
      parameter (nzeta=101)
      real*8 a_1,a_2,a_3,b_1,b_2,b_3
      common /voigt_cf/
     $     a_1(0:nzeta-1), a_2(0:nzeta-1), a_3(0:nzeta-1),
     $     b_1(0:nzeta-1), b_2(0:nzeta-1), b_3(0:nzeta-1)
C                                                                         D06390
      DIMENSION VNU(*),SABS(*),ALFAL(*),ALFAD(*),MOL(*),SPP(*),SRAD(*)    D06400
C                                                                         D06410
      DATA ZBND / 64. /                                                   D06420
      DATA ASUBL / 0.623 /,BSUBL / 0.410 /                                D06430
      DATA HREJ /'0'/,HNOREJ /'1'/
C                                                                         D06440
      VNULST = V2R4+BOUND4                                                D06450
C                                                                         D06460
      IF (JCNVF4.NE.0) GO TO 20                                           D06470
      JCNVF4 = 1                                                          D06480
C                                                                         D06490
C     SET UP CHI SUB-LORENTZIAN FORM FACTOR FOR CARBON DIOXIDE            D06500
C     POLYNOMIAL MATCHED TO AN EXPONENTIAL AT X0 = 8 CM-1                 D06510
C                                                                         D06520
      X0 = 8.                                                             D06530
      Y0 = EXP(-ASUBL*X0**BSUBL)                                          D06540
      F = ASUBL*BSUBL*X0**(BSUBL-1)                                       D06550
      Y1 = -F*Y0                                                          D06560
      Y2 = Y1*((BSUBL-1)/X0-F)                                            D06570
      Z0 = (Y0-1)/X0**2                                                   D06580
      Z1 = Y1/(2*X0)                                                      D06590
      Z2 = Y2/2.                                                          D06600
      C6 = (Z0-Z1+(Z2-Z1)/4.)/X0**4                                       D06610
      C4 = (Z1-Z0)/X0**2-2.*X0**2*C6                                      D06620
      C2 = Z0-X0**2*C4-X0**4*C6                                           D06630
      DVCHI = 0.1                                                         D06640
      RDVCHI = 1./DVCHI                                                   D06650
C                                                                         D06660
      DO 10 ISUBL = 1, 251                                                D06670
         FI = DVCHI*FLOAT(ISUBL-1)                                        D06680
         IF (FI.LT.X0) THEN                                               D06690
            CHI(ISUBL) = 1.+C2*FI**2+C4*FI**4+C6*FI**6                    D06700
         ELSE                                                             D06710
            FNI = ASUBL*(FI**BSUBL)                                       D06720
            CHI(ISUBL) = EXP(-FNI)                                        D06730
         ENDIF                                                            D06740
   10 CONTINUE                                                            D06750
C                                                                         D06760
C     CONSTANTS FOR FOURTH FUNCTION LINE SHAPE                            D06770
C                                                                         D06780
      RECPI = 1./(2.*ASIN(1.))                                            D06790
      ZSQBND = ZBND*ZBND                                                  D06800
C                                                                         D06830
   20 CONTINUE                                                            D06840
C                                                                         D06850
      BNDSQ = BOUND4*BOUND4                                               D06860
C                                                                         D06870
C     START OF LOOP OVER LINES                                            D06880
C                                                                         D06890
      IF (ILNFLG.EQ.2) READ(16)(FREJ(I),I=ILO,IHI)
C
      DO 60 I = ILO, IHI                                                  D06900
C                                                                         D06910
         IF (SABS(I).EQ.0..AND.SPP(I).EQ.0.) GO TO 60                     D06920
         ALFADI = ALFAD(I)                                                D06930
         ALFALI = ALFAL(I)                                                D06940
         ZETAI = ALFALI/(ALFALI+ALFADI)                                   D06950
         IZ = 100.*ZETAI + ONEPL                                          D06960
         ZETDIF = 100.*ZETAI - FLOAT(IZ-1)
         ALFAVI = ( AVRAT(IZ) + ZETDIF*(AVRAT(IZ+1)-AVRAT(IZ)) ) *        D06970
     x            (ALFALI+ALFADI)
C                                                                         D07010
c        Interpolate coefficients to proper zeta
c        a_3, ... are zero indexed, other variables start at one
         izx=iz-1
         if (izx .eq. 101) then
            izx = 100
            zetdif = 1.0
         endif 
         a3x=a_3(izx)+zetdif*(a_3(izx+1)-a_3(izx))
         b3x=b_3(izx)+zetdif*(b_3(izx+1)-b_3(izx))
         f4_64x= (a3x + b3x*zsqbnd)
c
         RALFVI = 1./ALFAVI                                               D06980
         SIV= SABS(I)*RALFVI
         siv_a = siv*a3x
         siv_b = siv*b3x 
         SRV= SRAD(I)*RALFVI 
         srv_a = srv*a3x
         srv_b = srv*b3x

c nlte line coupling constant
         cupcon=spp(i)/sabs(i)
C                                                                         D07010
         SPEAK = A3x*(ABS(SIV))                                            D07030
C                                                                         D07090
         JJ = (VNU(I)-V1R4)/DVR4+1.                                       D07100
         JJ = MAX(JJ,1)                                                   D07110
         JJ = MIN(JJ,NPTR4)                                               D07120
C
         IF (ILNFLG.LE.1) THEN
            FREJ(I) = HNOREJ
c     No rejection for line-coupled lines (SPP ne. 0)
            IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ)) .and. spp(i).eq.0.) THEN
                FREJ(I) = HREJ
                GO TO 60                                                  D07130
            ENDIF
         ELSE
            IF (FREJ(I).EQ.HREJ) GOTO 60
         ENDIF
C
         ILIN4 = ILIN4+1                                                  D07140
C                                                                         D07150
         VNUI = VNU(I)                                                    D07160
C                                                                         D07170
   30    CONTINUE                                                         D07180
C                                                                         D07190
         XNUI = VNUI-V1R4                                                 D07200
         JMIN = (XNUI-BOUND4)/DVR4+2.                                     D07210
C                                                                         D07220
         IF (VNUI.GE.VNULST) GO TO 70                                     D07230
         IF (JMIN.GT.NPTR4) GO TO 60                                      D07240
         JMIN = MAX(JMIN,1)                                               D07250
         JMAX = (XNUI+BOUND4)/DVR4+1.                                     D07260
         IF (JMAX.LT.JMIN) GO TO 50                                       D07270
         JMAX = MIN(JMAX,NPTR4)                                           D07280
         ALFLI2 = ALFALI*ALFALI                                           D07290
         ALFVI2 = ALFAVI*ALFAVI                                           D07300
         XJJ = FLOAT(JMIN-1)*DVR4                                         D07310
         siv_64 = siv*f4_64x*(alfli2 + alfvi2*zsqbnd)
         srv_64 = srv*f4_64x*(alfli2 + alfvi2*zsqbnd)
         f4bnd  = siv_64/(ALFLI2+bndSQ)
         frbnd  = srv_64/(ALFLI2+bndSQ)
C                                                                         D07340
C                FOURTH FUNCTION CONVOLUTION                              D07350
C                                                                         D07360

         dptrat = spp(i)/(sabs(i)*alfavi)
         dptrat_r =  spp(i)/(srad(i)*alfavi)
         rec_alfvi2 = 1./ALFVI2
C
         DO 40 JJ = JMIN, JMAX                                            D07370
            XM = (XJJ-XNUI)                                               D07380
            XMSQ = XM*XM                                                  D07390
            ZVSQ = XMSQ * rec_ALFVI2 
C                                                                         D07410
            IF (ZVSQ.LE.ZSQBND) THEN                                      D07420
               F4FN = (siv_A + ZVSQ * siv_B) - F4BND
               F4FR = (srv_A + ZVSQ * srv_B) - FRBND
               IF (SPP(I).NE.0.) THEN                                     D07440
                   F4FN = f4fn + xm*dptrat*f4fn
                   F4FR = F4FR + xm*dptrat_r*f4fr
               ENDIF
            ELSE                                                          D07460
                f4fn = siv_64/(alfli2+xmsq) - f4bnd
                f4fr = srv_64/(alfli2+xmsq) - frbnd
               IF (SPP(I).NE.0.) THEN                                     D07480
                   F4FN = F4FN+XM*dptrat*f4fn
                   F4FR = F4FR+XM*dptrat_r*f4fr
               ENDIF
            ENDIF                                                         D07500
C                                                                         D07510
            IF (MOL(I).EQ.2.AND.SPP(I).EQ.0.) THEN                        D07520
C                                                                         D07530
C     ASSIGN ARGUMENT ISUBL OF THE FORM FACTOR FOR CO2 LINES              D07540
C                                                                         D07550
               ISUBL = RDVCHI*ABS(XM)+1.5                                 D07560
               ISUBL = MIN(ISUBL,251)                                     D07570
C                                                                         D07580
               R4(JJ) = R4(JJ)+F4FN*CHI(ISUBL)                            D07590
               RR4(JJ) = RR4(JJ)+F4FR*CHI(ISUBL)                          D07590
            ELSE                                                          D07600
               R4(JJ) = R4(JJ)+F4FN                                       D07610
               RR4(JJ) = RR4(JJ)+F4FR                                     D07610
            ENDIF                                                         D07620
C                                                                         D07630
C                                                                         D07640
            XJJ = XJJ+DVR4                                                D07650
   40    CONTINUE                                                         D07660
C                                                                         D07670
   50    IF (VNUI.GT.0..AND.VNUI.LE.25.) THEN                             D07680
C                                                                         D07690
C     THE CALCULATION FOR NEGATIVE VNU(I) IS FOR VAN VLECK WEISSKOPF      D07700
C                                                                         D07710
            VNUI = -VNU(I)                                                D07720
c
            GO TO 30                                                      D07750
C                                                                         D07760
         ENDIF                                                            D07770
C                                                                         D07780
   60 CONTINUE                                                            D07790
C                                                                         D07800
      IF (ILNFLG.EQ.1) WRITE(16)(FREJ(I),I=ILO,IHI)
      RETURN                                                              D07810
C                                                                         D07820
C     IF END OF CONVOLUTION, SET IHI=-1 AND RETURN                        D07830
C                                                                         D07840
   70 CONTINUE
      IF (ILNFLG.EQ.1) WRITE(16)(FREJ(I),I=ILO,IHI)
      IHI = -1                                                            D07850
C                                                                         D07860
      RETURN                                                              D07870
C                                                                         D07880
      END                                                                 D07890
      SUBROUTINE SHRINQ                                                   D03250
C                                                                         D03260
      IMPLICIT REAL*8           (V)                                     ! D03270
C                                                                         D03280
C     SUBROUTINE SHRINK COMBINES LINES FALLING IN A WAVENUMBER INTERVAL   D03290
C     DVR4/2 INTO A SINGLE EFFECTIVE LINE TO REDUCE COMPUTATION           D03300
C                                                                         D03310
      COMMON VNU(1250),S(1250),ALFAL(1250),ALFAD(1250),MOL(1250),         D03320
     *       SPP(1250)                                                    D03330
      COMMON /R4SUB/ VLO,VHI,ILO,IST,IHI,LIMIN,LIMOUT,ILAST,DPTMN,        D03340
     *               DPTFC,ILIN4,ILIN4T                                   D03350
      COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)        D03360
C                                                                         D03370
      J = ILO-1                                                           D03380
      DV = DVR4/2.                                                        D03390
      VLMT = VNU(ILO)+DV                                                  D03400
C                                                                         D03410
C     INITIALIZE NON-CO2 SUMS                                             D03420
C                                                                         D03430
      SUMAL = 0.                                                          D03440
      SUMAD = 0.                                                          D03450
      SUMS = 0.                                                           D03460
      SUMV = 0.                                                           D03470
      SUMC = 0.                                                           D03480
C                                                                         D03490
C     INITIALIZE CO2 SUMS                                                 D03500
C                                                                         D03510
      SUMAL2 = 0.                                                         D03520
      SUMAD2 = 0.                                                         D03530
      SUMS2 = 0.                                                          D03540
      SUMV2 = 0.                                                          D03550
      SUMC2 = 0.                                                          D03560
C                                                                         D03570
      DO 20 I = ILO, IHI                                                  D03580
C                                                                         D03590
C     IF LINE COUPLING, DON'T SHRINK LINE                                 D03600
C                                                                         D03610
         IF (SPP(I).NE.0.0) THEN                                          D03620
            J = J+1                                                       D03630
            VNU(J) = VNU(I)                                               D03640
            S(J) = S(I)                                                   D03650
            ALFAL(J) = ALFAL(I)                                           D03660
            ALFAD(J) = ALFAD(I)                                           D03670
            SPP(J) = SPP(I)                                               D03680
            MOL(J) = MOL(I)                                               D03690
            IF (MOL(J).NE.2) MOL(J) = 0                                   D03700
            GO TO 10                                                      D03710
         ENDIF                                                            D03720
C                                                                         D03730
C     NON-CO2 LINES OF MOLECULAR INDEX IT.NE.2   ARE LOADED               D03740
C     INTO SUMS IF THE FREQUENCY WITHIN DV GROUP                          D03750
C                                                                         D03760
         IF (MOL(I).NE.2) THEN                                            D03770
            SUMV = SUMV+VNU(I)*S(I)                                       D03780
            SUMS = SUMS+S(I)                                              D03790
            SUMAL = SUMAL+S(I)*ALFAL(I)                                   D03800
            SUMAD = SUMAD+S(I)*ALFAD(I)                                   D03810
            SUMC = SUMC+SPP(I)                                            D03820
         ELSE                                                             D03830
C                                                                         D03840
C     CO2 LINES LOADED     (MOL .EQ. 2)                                   D03850
C                                                                         D03860
            SUMV2 = SUMV2+VNU(I)*S(I)                                     D03870
            SUMS2 = SUMS2+S(I)                                            D03880
            SUMAL2 = SUMAL2+S(I)*ALFAL(I)                                 D03890
            SUMAD2 = SUMAD2+S(I)*ALFAD(I)                                 D03900
            SUMC2 = SUMC2+SPP(I)                                          D03910
         ENDIF                                                            D03920
C                                                                         D03930
C     IF LAST LINE OR VNU GREATER THAN LIMIT THEN STORE SUMS              D03940
C                                                                         D03950
   10    IF (I.LT.IHI) THEN                                               D03960
            IF (VNU(I+1).LE.VLMT) GO TO 20                                D03970
         ENDIF                                                            D03980
C                                                                         D03990
         VLMT = VNU(I)+DV                                                 D04000
C                                                                         D04010
C     ASSIGN NON-CO2 LINE AVERAGES TO 'GROUP' LINE J                      D04020
C                                                                         D04030
         IF (SUMS.GT.0.) THEN                                             D04040
            J = J+1                                                       D04050
            S(J) = SUMS                                                   D04060
            ALFAL(J) = SUMAL/SUMS                                         D04070
            ALFAD(J) = SUMAD/SUMS                                         D04080
            VNU(J) = SUMV/SUMS                                            D04090
            SPP(J) = SUMC                                                 D04100
            MOL(J) = 0                                                    D04110
            SUMAL = 0.                                                    D04120
            SUMAD = 0.                                                    D04130
            SUMS = 0.                                                     D04140
            SUMV = 0.                                                     D04150
            SUMC = 0.                                                     D04160
         ENDIF                                                            D04170
C                                                                         D04180
C     ASSIGN CO2 LINE AVERAGES                                            D04190
C                                                                         D04200
         IF (SUMS2.GT.0.) THEN                                            D04210
            J = J+1                                                       D04220
            S(J) = SUMS2                                                  D04230
            ALFAL(J) = SUMAL2/SUMS2                                       D04240
            ALFAD(J) = SUMAD2/SUMS2                                       D04250
            VNU(J) = SUMV2/SUMS2                                          D04260
            MOL(J) = 2                                                    D04270
            SPP(J) = SUMC2                                                D04280
            SUMAL2 = 0.                                                   D04290
            SUMAD2 = 0.                                                   D04300
            SUMS2 = 0.                                                    D04310
            SUMV2 = 0.                                                    D04320
            SUMC2 = 0.                                                    D04330
         ENDIF                                                            D04340
C                                                                         D04350
   20 CONTINUE                                                            D04360
C                                                                         D04370
      ILO = J+1                                                           D04380
      IHI = J                                                             D04390
C                                                                         D04400
      RETURN                                                              D04410
C                                                                         D04420
      END                                                                 D04430
c ----------------------------------------------------------------

      SUBROUTINE VIBPOP(XKT,HT,NLTEFLAG,NUM,HOL,NDEG,EH,RAT,TITMOL)       604860
C                                                                         604870
C                                                                         604880
C     SUBROUTINE VIBPOP USES THE NON-LTE POPULATION DATA FROM             604890
C     TAPE4 TO CALCULATE THE VIBRATIONAL POPULATION ENHANCEMENT           604900
C     RATIOS FOR SELECTED VIBRATIONAL STATES OF H2O,CO2,NO AND O3.        604910
C                                                                         604920
      IMPLICIT REAL*8           (V)                                     ! B00030

      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         B00840
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        B00850
     *              NLTEFL,LNFIL4,LNGTH4                                  B00860

      CHARACTER*(*) HOL                                                   604960
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       B00660
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   B00670
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    B00680
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      DIMENSION HOL(26),VQNE(26),VQEQ(26),TNE(26),
     *          VPNE1(26),VPNE2(26),VQNEST(26)
      DIMENSION NDEG(*),EH(*),RAT(*)
      CHARACTER*10 TITMOL,HMNLTE                                          605020
C                                                                         605040
C     SKIP TO BEGINNING OF VIBRATIONAL DATA                               605050
C                                                                         605060
       CALL RDSKIP(NLTEFLAG)                                              605070
C                                                                         605080
C   READ NLTE VIB POPULATIONS                                             605090
C                                                                         605100
      XKT= TAVE/RADCN2                                                    605030
      READ (NLTEFLAG,902)  HMNLTE                                         605110
      READ (NLTEFLAG,904)  ALT1,(VPNE1(I),I=1,NUM)                        605120
  10  READ (NLTEFLAG,904)  ALT2,(VPNE2(I),I=1,NUM)                        605140
C*****WOG, 11/06/2000: ALT1 -> AL2:
C     IF( ALT1.LT.HT) THEN
      IF( ALT2.LE.HT) THEN
          ALT1 = ALT2                                                     605170
          DO 20 I=1,NUM                                                   605180
              VPNE1(I) = VPNE2(I)                                         605190
  20      CONTINUE
          GO TO 10                                                        605200
      ENDIF
C     CALL LININT(HT,ALT1,ALT2,NUM,VPNE1,VPNE2,VQNE)                      605220
      A = (HT-ALT1)/(ALT2-ALT1)                                           605230
      DO 30 I=1,NUM                                                       605240
          CALL EXPINT(VQNE(I),VPNE1(I),VPNE2(I),A )                       605250
  30  CONTINUE
C                                                                         605270
      POPEQ=0.0                                                           605280
      POPNE=0.0                                                           605290
C                                                                         605300
      DO 50 LVL=1,NUM                                                     605310
          VQEQ(LVL)=NDEG (LVL)*EXP(-EH(LVL)/XKT)                          605320
          VQNEST(LVL)=VQNE(LVL)                                           605330
          POPEQ=POPEQ + VQEQ(LVL)                                         605340
          POPNE=POPNE + VQNE (LVL)                                        605350
  50  CONTINUE                                                            605360
C                                                                         605370
C    NORMALIZE POPULATIONS AND CALCULATE RATIOS                           605380
C                                                                         605390
      WRITE(IPR,906) TITMOL                                               605400
      WRITE(IPR,935)                                                      605410
      DO 100 LVL=1,NUM                                                    605420
          I=LVL                                                           605430
          VQEQ(LVL)=VQEQ(LVL)/POPEQ                                       605440
          VQNE (LVL)=VQNE(LVL)/POPNE                                      605450
          RAT(I)=VQNE(I)/VQEQ(I)                                          605460
          IF(LVL.EQ.1) THEN                                         
              VST1=VQNE(1)                                                605480
              TNE(I)=TAVE                                                 605490
              WRITE(IPR,920)HOL(I),EH(I),VQEQ(I),VQNE(I),RAT(I),TNE(I),   605500
     &            VQNEST(I)                                               605510
          ELSE                                                            605520
              DEN=(NDEG(1)*VQNE(LVL)/(NDEG(LVL)*VST1))                    605530
              TNE(I)=-RADCN2*EH(LVL)/ LOG(DEN)                            605540
              WRITE(IPR,920)HOL(I),EH(I),VQEQ(I),VQNE(I),RAT(I),TNE(I),   605550
     &            VQNEST(I)                                               605560
          END IF                                                          605570
  100 CONTINUE                                                            605580
C                                                                         605590
      RETURN                                                              605600
C                                                                         605610
  902 FORMAT(A10)                                                         605620
  904 FORMAT (F7.0,1P,7E11.4,     /(18X, 6E11.4))                         605630
  906 FORMAT(//,A10,'  ENERGY LEVELS',10(/,20X,1PE11.4))                  605640
  920 FORMAT(2X,A10,4G15.5,F10.2,G15.5)                                   605660
  935 FORMAT (2X,'VIB',10X,'E(CM-1)',11X,'POP LTE',7X,                    605680
     & 'POP NLTE',6X,'NLTE/LTE',7X,'NLTE TMP',7X,'NLTE POP ORIG')         605690
C                                                                         605700
      END                                                                 605710
c ----------------------------------------------------------------

      SUBROUTINE VIBTMP(XKT,HT,NLTEFLAG,NUM,HOL,NDEG,EH,RAT,TITMOL)       605720
C                                                                         604870
C                                                                         604880
C     SUBROUTINE VIBTMP USES THE NON-LTE TEMPERATURE DATA FROM            604890
C     TAPE4 TO CALCULATE THE VIBRATIONAL POPULATION ENHANCEMENT           604900
C     RATIOS FOR SELECTED VIBRATIONAL STATES OF H2O, CO2, NO AND          604910
C     O3.  THE NLTE VIBRATIONAL TEMPERATURES ARE WITH RESPECT TO          604910
C     THE GROUND VIBRATIONAL STATE.                                       604910
C                                                                         604920
      IMPLICIT REAL*8           (V)                                     ! B00030

      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         B00840
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        B00850
     *              NLTEFL,LNFIL4,LNGTH4                                  B00860

      CHARACTER*(*) HOL                                                   604960
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       B00660
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   B00670
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    B00680
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      DIMENSION HOL(26),VQNE(26),VQEQ(26),TNE(26),TNESAV(26),
     *          TEM1(26),TEM2(26) 
      DIMENSION NDEG(*),EH(*),RAT(*)
      CHARACTER*10 TITMOL,HMNLTE                                          605820
C                                                                         605830
C     SKIP TO BEGINNING OF VIBRATIONAL DATA                               605840
C                                                                         605850
       CALL RDSKIP(NLTEFLAG)                                              605860
C                                                                         605870
C   READ NLTE VIB TEMPERATURE                                             605880
C                                                                         605890
      XKT= TAVE/RADCN2                                                    605900
      READ (NLTEFLAG,902)  HMNLTE                                         605910
      READ (NLTEFLAG,904)  ALT1,(TEM1(I),I=1,NUM)                         605920
  10  READ (NLTEFLAG,904)  ALT2,(TEM2(I),I=1,NUM)                         605930
C*****WOG, 11/06/2000: ALT1 -> AL2:
C     IF( ALT1.LT.HT) THEN
      IF( ALT2.LE.HT) THEN
          ALT1 = ALT2                                                     605960
          DO 20 I=1,NUM                                                   605970
              TEM1(I) = TEM2(I)                                           605980
  20      CONTINUE
          GO TO 10                                                        605990
      ENDIF   
C                                                                         606010
C     SET ZERO TEMP TO AMBIENT TEMP                                       606020
C                                                                         606030
      CALL  STAMB(NUM,TAVE,TEM1)                                          606040
      CALL  STAMB(NUM,TAVE,TEM2)                                          606050
      CALL LININT(HT,ALT1,ALT2,NUM,TEM1,TEM2,TNESAV)                      606060
      WRITE(IPR,906) TITMOL                                               606070
C                                                                               
CC     CORRECT TEMP TO ATMOSPHERIC                                              
C                                                                               
      RATTV = TAVE /TNESAV(1)
c loop 40 corrects the input temperatures when the input
c  level temperature does not match the computed layer
c  temperature
c      write(*,*) ' temperature not corrected'
      
      DO 40 I = 1,NUM 
          TNE(I) = TNESAV(I) * RATTV 
  40  CONTINUE
C                                                                               
      SUMQ=0                                                              606080
      SUMNQ=0                                                             606090
      DO 50 I=1,NUM                                                       606100
          VQNE(I)=1.                                                      606110
          IF (TNE(I).GT.0.0) 
     &        VQNE(I)=NDEG(I)*EXP(-RADCN2*EH(I)/TNE(I))                   606130
          VQEQ(I)=NDEG(I)*EXP(-EH(I)/XKT)                                 606140
          SUMQ=SUMQ+VQEQ(I)                                               606150
          SUMNQ=SUMNQ+VQNE(I)                                             606160
  50  CONTINUE                                                            606170
      WRITE(IPR,935)                                                      606180
      DO 100 I=1,NUM                                                      606190
          VQNE(I)=VQNE(I)/SUMNQ                                           606200
          VQEQ(I)=VQEQ(I)/SUMQ                                            606210
          RAT(I)=VQNE(I)/VQEQ(I)                                          606220
          WRITE(IPR,920)HOL(I),EH(I),VQEQ(I),VQNE(I),RAT(I),
     &        TNESAV(I), TNE(I)                                           606230
  100 CONTINUE                                                            606240
C                                                                         606250
      RETURN                                                              606260
C                                                                         606270
  902 FORMAT(A10)                                                         606280
  904 FORMAT(F7.0,7F11.3/(18X,6F11.3))                                    606290
  906 FORMAT(//,5X,A10,'  ENERGY LEVELS',10(/,20X,1PE11.4))               606300
  920 FORMAT(2X,A10,4G12.5,2F9.2)                                         606320
  935 FORMAT (9X,'VIB E(CM-1)',9X,'POP LTE    POP NLTE NLTE/LTE',         606340
     &   '    NLTE TMP 2-STATE NLTE TMP')                                 606350
C                                                                         606360
      END                                                                 606370

c ----------------------------------------------------------------

      SUBROUTINE RDSKIP(NTAPE)                                            606380
      CHARACTER *1 HRD,HMINUS                                             606390
      DATA HMINUS /'-'/                                                   606400
  10  READ(NTAPE,900) HRD                                                 606410
  900 FORMAT(2X,A1)                                                       606420
      IF(HRD.EQ.HMINUS) RETURN                                            606430
      GO TO 10                                                            606440
      END                                                                 606450

c ----------------------------------------------------------------

      SUBROUTINE LININT(HT,ALT1,ALT2,NUM,T1,T2,TNE)                       606460
      DIMENSION T1(26),T2(26),TNE(26)                                     606470
C*****WOG 11/03/2000
C*****Correct for divide by zero if two altitudes are the same
C*****0.001 = 1 meter, small enough to use average, large enough to
c     prevent numerical error. 
      IF (ABS(ALT2-ALT1) .LE. 0.001) THEN
          AM = 0.5
      ELSE 
          AM = (HT-ALT1)/(ALT2-ALT1) 
      ENDIF
      
      DO 10 I=1,NUM
          TNE(I)=T1(I)+AM*(T2(I)-T1(I))
  10  CONTINUE

C     DO 10 I=1,NUM                                                       606480
C         AM = (T2(I)-T1(I))/(ALT2-ALT1)                                  606490
C         C=T1(I)-AM*ALT1                                                 606500
C         TNE(I)=AM*HT+C                                                  606510
C  10  CONTINUE
      RETURN                                                              606520
      END                                                                 606530

c ----------------------------------------------------------------

      SUBROUTINE  STAMB(NUM,T1AMB,T1)                                     606540
      DIMENSION T1(*)                                                    ?606550

C*****WOG, 11/3/2000
C*****Why start from 2 instead of 1???
C     DO 10 I=2,NUM                                                       606560
      DO 10 I=1,NUM
          IF(T1(I).LE.0.) T1(I)=T1AMB                                     606570
  10  CONTINUE                                                            606580
      RETURN                                                              606590
      END                                                                 606600

