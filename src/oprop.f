C     path: %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      SUBROUTINE HIRAC1 (MPTS)                                            B00010
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
      COMMON /HVERSN/  HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,
     *                HVROPR,HVRPST,HVRPLT,HVRTST,HVRUTL,HVRXMR,hvnlte
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       B00660
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   B00670
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    B00680
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           B00690
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
      PARAMETER (NFPTS=2001,NFMX=1.3*NFPTS)
C
      COMMON /FNSH/ IFN,F1(NFMX),F2(NFMX),F3(NFMX),FG(NFMX),XVER(NFMX)    B00810
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
C                                                                         B00970
      REAL L4TIM,L4TMR,L4TMS,LOTHER
      CHARACTER*55 CDUM1,PTHODI,PTHODT,PTHRDR
      CHARACTER*10 HFMODL
      CHARACTER CFORM*11,KODLYR*57,PTHODE*55,PTHODD*55                    B00980
      CHARACTER*15 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,HVROPR,
     *            HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR,hvnlte
      LOGICAL OP                                                          B00990
C                                                                         B01000
      DIMENSION MEFDP(64),FILHDR(2),IWD(2)                                B01010
      DIMENSION R1(4050),R2(1050),R3(300)
C                                                                         B01020
      EQUIVALENCE (IHIRAC,FSCDID(1)) , (ILBLF4,FSCDID(2)),                B01030
     *            (IXSCNT,FSCDID(3)) , (IAERSL,FSCDID(4)),                B01040
     *            (JRAD,FSCDID(9)) , (IMRG,FSCDID(11)),                   B01050
     *            (IATM,FSCDID(15)) , (YI1,IOD) , (XID(1),FILHDR(1)),     B01060
     *            (V1P,IWD(1)) , (NPNLXP,LSTWDX)                          B01070
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
      PTHODE = 'ODexact_'
      PTHODD = 'ODdeflt_'
      DATA KODLYR /
     *     '                                                         '/
      DATA HFMODL /'         '/
C                                                                         B01160
      CALL CPUTIM (TIMEH0)                                                B01170
C                                                                         B01180
C     ASSIGN CVS VERSION NUMBER TO MODULE 
C
      HVROPR = '$Revision$'
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
         CALL SHAPEL (F1,F2,F3)                                           B01610
         CALL SHAPEG (FG)                                                 B01620
         CALL VERFN (XVER)                                                B01630
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
   10 CONTINUE                                                            B02070
      DO 20 I = 1, MAX2                                                   B02080
         R2(I) = 0.                                                       B02090
   20 CONTINUE                                                            B02100
      DO 30 I = 1, MAX3                                                   B02110
         R3(I) = 0.                                                       B02120
   30 CONTINUE                                                            B02130
      IF (ILBLF4.EQ.0) THEN                                               B02140
         DO 40 I = 1, 2502                                                B02150
            R4(I) = 0.                                                    B02160
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
      IF (ILBLF4.GE.1) CALL LBLF4 (JRAD,V1R4ST,V2R4ST)                    B02610
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
      CALL LNCOR1 (NLNCR,IHI,ILO,MEFDP)                                   B02820
      CALL CPUTIM(TPAT1)
      TLNCOR = TLNCOR+TPAT1-TPAT0
C                                                                         B02830
   70 CONTINUE                                                            B02840
C                                                                         B02850
      CALL CNVFNV (VNU,SP,SPPSP,RECALF,R1,R2,R3,F1,F2,F3,FG,XVER,ZETAI,   B02860
     *             IZETA)                                                 B02870
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
      IF (ILBLF4.GE.1)                                                    B03290
     *    CALL XINT (V1R4,V2R4,DVR4,R4,1.0,VFT,DVR3,R3,N1R3,N2R3)         B03300
      IF (ICNTNM.GE.1)                                                    B03310
     *    CALL XINT (V1ABS,V2ABS,DVABS,ABSRB,1.,VFT,DVR3,R3,N1R3,N2R3)    B03320
      CALL CPUTIM(TPAT1)
      TXINT = TXINT + TPAT1-TPAT0
C                                                                         B03330
      CALL PANEL (R1,R2,R3,KFILE,JRAD,IENTER)                             B03340
C                                                                         B03350
      IF (ISTOP.NE.1) THEN                                                B03360
         IF (ILBLF4.GE.1) THEN                                            B03370
            VF1 = VFT-2.*DVR4                                             B03380
            VF2 = VFT+2.*DVR4+FLOAT(N2R3+4)*DVR3                          B03390
            IF (VF2.GT.V2R4.AND.V2R4.NE.V2R4ST) THEN                      B03400
               CALL LBLF4 (JRAD,VF1,V2R4ST)                               B03410
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
  925 FORMAT ('0  * HIRAC1 *  AVERAGE WIDTH = ',F8.6,                     B03850
     *        ',  AVERAGE ZETA = ',F8.6)                                  B03860
  930 FORMAT ('0 ********  HIRAC1  ********',I5,' STRENGTHS FOR',         B03870
     *        '  TRANSITIONS WITH UNKNOWN EPP FOR MOL =',I5,              B03880
     *        ' SET TO ZERO')                                             B03890
 935  FORMAT (/,'0     + OTHER timing includes:',/,
     *          '0             In LINF4:  MOLEC, BUFIN, BUFOUT, ',
     *          'NWDL, ENDFIL, and SHRINK',/,
     *          '0             In HIRAC:  LNCOR, XINT, SHAPEL, ',
     *          'SHAPEG, VERFN, MOLEC, and other loops and ',
     *          'file maintenance within HIRAC',/)
C                                                                         B03900
      END                                                                 B03910
      BLOCK DATA BHIRAC                                                   B03920
C
      PARAMETER (NFPTS=2001,NFMX=1.3*NFPTS)
C                                                                         B03930
      COMMON /FNSH/ IFN,F1(NFMX),F2(NFMX),F3(NFMX),FG(NFMX),XVER(NFMX)    B03940
C                                                                         B03950
      DATA IFN / 0 /                                                      B03960
C                                                                         B03970
      END                                                                 B03980
      SUBROUTINE RDLIN                                                    B03990
C                                                                         B04000
      IMPLICIT REAL*8           (V)                                     ! B04010
C                                                                         B04020
C     SUBROUTINE RDLIN INPUTS LINE DATA FROM FILE LINFIL                  B04030
C                                                                         B04040
      CHARACTER*8      HLINID,BMOLID,HID1
      CHARACTER*1 CNEGEPP(8)

      integer *4 molcnt,mcntlc,
     *           mcntnl,linmol,
     *           lincnt,ilinlc,ilinnl,irec,irectl
      real *4 sumstr,flinlo,flinhi
c
      COMMON /LINHDR/ HLINID(10),BMOLID(64),MOLCNT(64),MCNTLC(64),
     *                MCNTNL(64),SUMSTR(64),LINMOL,FLINLO,FLINHI,
     *                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1(2),LSTWDL

      COMMON VNU(250),SP(250),ALFA0(250),EPP(250),MOL(250),HWHMS(250),    B04050
     *       TMPALF(250),PSHIFT(250),IFLG(250),SPPSP(250),RECALF(250),    B04060
     *       ZETAI(250),IZETA(250)                                        B04070

      dimension    amol(250)
      equivalence (mol(1),amol(1))

      common /rdlnpnl/ vmin,vmax,nrec,nwds
      integer *4 nrec,nwds,lnfl,leof,npnlhd

      common /rdlnbuf/ vlin(250),str(250),hw_f(250),e_low(250),
     *     mol_id(250),hw_s(250),hw_T(250),shft(250),jflg(250)
      dimension xmol(250)
      equivalence (vmin,rdpnl(1)),(mol_id(1),xmol(1))

      real *4 str,hw_f,e_low,xmol,hw_s,hw_T,shft,rdpnl(2),dum(2)
      integer *4 mol_id,jflg,n_one

      COMMON /XSUB/ VBOT,VTOP,VFT,LIMIN,ILO,IHI,IEOF,IPANEL,ISTOP,IDATA   B04080
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         B04100
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        B04110
     *              NLTEFL,LNFIL4,LNGTH4                                  B04120
      COMMON /IOU/ IOUT(250)                                              B04130
      common /eppinfo/ negepp_flag
      common /bufid2/ n_negepp(64),n_resetepp(64),xspace(4096),lstwdl2
      integer *4 negepp_flag,n_negepp,n_resetepp
      real *4 xspace
C                                                                         B04150

*****************************************************************************
c
      data n_one/ 1/ npnlhd/ 6/


      lnfl = linfil
C                                                                         B04170
C     THERE ARE (LIMIN * 9) QUANTITIES READ IN:                           B04180
C     VNU,SP,ALFA0,EPP,MOL,HWHMS,TMPALF,PSHIFT,IFLG                       B04190
C                                                                         B04200
      ILNGTH = NLNGTH*LIMIN                                               B04210
      IDATA = 0                                                           B04220
C                                                                         B04230
C     BUFFER PAST FILE HEADER if necessary
C                                                                         B04250
      IF (ILO.LE.0) THEN                                                  B04260
         REWIND LNFL
         read (lnfl)    HLINID
         READ (HLINID(7),950) CNEGEPP
         IF (CNEGEPP(8).eq.'^') THEN
            read (lnfl) n_negepp,n_resetepp,xspace
         endif
      ENDIF                                                               B04290
C                                                                         B04300
   10 CALL BUFINln (Lnfl,LEOF,rdpnl(1),npnlhd)
      IF (LEOF.EQ.0) THEN                                                 B04320
         IF (NOPR.EQ.0) WRITE (IPR,900)                                   B04330
         IEOF = 1                                                         B04340
         RETURN                                                           B04350
      ENDIF                                                               B04360
C                                                                         B04370
      IF (NREC.GT.LIMIN) STOP 'RDLIN; NREC GT LIMIN'                      B04380
c
      IF (VMAX.LT.VBOT) THEN                                              B04390
         CALL BUFINln (Lnfl,LEOF,DUM(1),n_one)
         GO TO 10                                                         B04410
      ENDIF                                                               B04420
c
      CALL BUFINln (Lnfl,LEOF,vlin(1),NWDS) 
c
c     precision conversion occurs here:
c     incoming on right: vlin is real*8, others real*4 and integer*4
c
      do 15 i=1,nrec

         IFLG(i)  = jflg(i)     
         VNU(i)   = vlin(i)
         SP(i)    = str(i)   
         ALFA0(i) = hw_f(i)      
         EPP(i)   = e_low(i)    
         if (iflg(i) .ge.  0) then
            MOL(i)   = mol_id(i)    
         else
            amol(i)  = xmol(i)
         endif
         HWHMS(i) = hw_s(i)      
         TMPALF(i)= hw_T(i)       
         PSHIFT(i)= shft(i)       

 15   continue

      IF ((ILO.EQ.0).AND.(VMIN.GT.VBOT)) WRITE (IPR,905)                  B04440
      ILO = 1                                                             B04450
C                                                                         B04460
      IJ = 0                                                              B04470
      DO 20 I = 1, NREC                                                   B04480
         IF (IFLG(I).GE.0) THEN                                           B04490
            IJ = IJ+1                                                     B04500
            IOUT(IJ) = I                                                  B04510
         ENDIF                                                            B04520
   20 CONTINUE                                                            B04530
C                                                                         B04540
      DO 30 I = IJ+1, 250                                                 B04550
         IOUT(I) = NREC                                                   B04560
   30 CONTINUE                                                            B04570
C                                                                         B04580
      IF (VMIN.LT.VBOT) THEN                                              B04590
         DO 40 J = 1, IJ                                                  B04600
            I = IOUT(J)                                                   B04610
            ILO = J                                                       B04620
            IF (VNU(I).GE.VBOT) GO TO 50                                  B04630
   40    CONTINUE                                                         B04640
      ENDIF                                                               B04650
C                                                                         B04660
   50 CONTINUE                                                            B04670
      DO 60 J = ILO, IJ                                                   B04680
         I = IOUT(J)                                                      B04690
         IF (MOL(I).GT.0) THEN                                            B04700
            IHI = J                                                       B04710
            IF (VNU(I).GT.VTOP) GO TO 70                                  B04720
         ENDIF                                                            B04730
   60 CONTINUE                                                            B04740

c     the following test is to see if more data is required
c     idata = 1 means data requirements have been met

   70 IF (IHI.LT.IJ) IDATA = 1                                            B04750
C                                                                         B04760
      RETURN                                                              B04770
C                                                                         B04780
  900 FORMAT ('  EOF ON LINFIL (MORE LINES MAY BE REQUIRED) ')            B04790
  905 FORMAT (                                                            B04800
     *   ' FIRST LINE ON LINFIL USED (MORE LINES MAY BE REQUIRED) ')      B04810
 950  FORMAT (8a1)
C                                                                         B04820
      END                                                                 B04830
      SUBROUTINE LNCOR1 (NLNCR,IHI,ILO,MEFDP)                             B04840
C                                                                         B04850
      IMPLICIT REAL*8           (V)                                     ! B04860
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
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           B05040
      COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)        B05050
      COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,
     *                HWF3,DXF3,NX3,N3MAX
      COMMON /VOICOM/ AVRAT(102),CGAUSS(102),CF1(102),CF2(102),           B05060
     *                CF3(102),CER(102)                                   B05070
C                                                                         B05080
      PARAMETER (NTMOL=36,NSPECI=85) 
C                                                                         B05100
      COMMON /ISVECT/ ISOVEC(NTMOL),ISO82(NSPECI),ISONM(NTMOL),           B05110
     *                SMASSI(NSPECI)                                      B05120
      COMMON /LNC1/ RHOSLF(NSPECI),ALFD1(NSPECI),SCOR(NSPECI),ALFMAX,     B05130
     *              BETACR,DELTMP,DPTFC,DPTMN,XKT,NMINUS,NPLUS,NLIN,      B05140
     *              LINCNT,NCHNG,SUMALF,SUMZET,TRATIO,RHORAT,PAVP0,       B05150
     *              PAVP2,RECTLC,TMPDIF,ILC                               B05160
      DIMENSION MEFDP(64),FILHDR(2),AMOL(250),SP(250)                     B05170
      DIMENSION A(4),B(4),TEMPLC(4)                                       B05180
C                                                                         B05190
      EQUIVALENCE (MOL(1),AMOL(1)) , (S(1),SP(1))                         B05200
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
         IF (IFLAG.EQ.0) THEN                                             B06640
            IF (ILNFLG.LE.1) THEN
               FREJ(J) = HNOREJ
               SPEAK = SUI*RECALF(I)                                      B06650
               IF (DVR4.LE.0.) THEN                                       B06660
                  IF (SPEAK.LE.DPTMN) THEN 
                     FREJ(J) = HREJ
                     GO TO 25
                  ENDIF
               ELSE                                                       B06680
                  JJ = (VNU(I)-V1R4)/DVR4+1.                              B06690
                  JJ = MAX(JJ,1)                                          B06700
                  JJ = MIN(JJ,NPTR4)                                      B06710
                  IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ))) THEN
                     FREJ(J) = HREJ
                     GO TO 25
                  ENDIF
               ENDIF                                                      B06730
            ELSE
C      "ELSE" IS TRUE WHEN "ILNFLG" EQUALS 2
C
               IF (FREJ(J).EQ.HREJ) GO TO 25
            ENDIF
         ENDIF                                                            B06790
C                                                                         B06800
         NMINUS = NMINUS+NMINAD                                           B06810
         NPLUS = NPLUS+NPLSAD                                             B06820
         SUMALF = SUMALF+ALFV                                             B06830
         SUMZET = SUMZET+ZETA                                             B06840
         LINCNT = LINCNT+1                                                B06850
C                                                                         B06860
         SP(I) = SUI*(1.+GI*PAVP2)                                        B06870
         SPPI = SUI*YI*PAVP0                                              B06880
         SPPSP(I) = SPPI/SP(I)                                            B06890
C
         GO TO 30
C
   25    SP(I) = 0.
         SPPSP(I) = 0.  
C                                                                         B06900
   30 CONTINUE                                                            B06910
C                                                                         B06920
      NCHNG = NMINUS+NPLUS                                                B06930
      IF (ILNFLG.EQ.1) WRITE(15)(FREJ(J),J=ILO,IHI)
C                                                                         B06940
      RETURN                                                              B06950
C                                                                         B06960
      END                                                                 B06970
      SUBROUTINE CNVFNV (VNU,SP,SPPSP,RECALF,R1,R2,R3,F1,F2,F3,FG,        B06980
     *                   XVER,ZETAI,IZETA)                                B06990
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
      DIMENSION VNU(*),SP(*),SPPSP(*),RECALF(*)                           B07530
      DIMENSION R1(*),R2(*),R3(*)                                         B07540
      DIMENSION F1(*),F2(*),F3(*)                                         B07550
      DIMENSION FG(*),XVER(*)                                             B07560
      DIMENSION IZETA(*),ZETAI(*)                                         B07570
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
            IF (SP(I).NE.0.) THEN                                         B07730
               DEPTHI = SP(I)*RECALF(I)                                   B07740
               IZM = IZETA(I)                                             B07750
               ZETDIF = 100.*ZETAI(I)-FLOAT(IZM-1)                        B07840
               STRF1 = DEPTHI*(CF1(IZM)+ZETDIF*(CF1(IZM+1)-CF1(IZM)))     B07850
               STRF2 = DEPTHI*(CF2(IZM)+ZETDIF*(CF2(IZM+1)-CF2(IZM)))     B07860
               STRF3 = DEPTHI*(CF3(IZM)+ZETDIF*(CF3(IZM+1)-CF3(IZM)))     B07870
               STRD = DEPTHI*(CGAUSS(IZM)+ZETDIF*(CGAUSS(IZM+1)-          B07880
     *               CGAUSS(IZM)))                                        B07890
               STRVER = DEPTHI*(CER(IZM)+ZETDIF*(CER(IZM+1)-CER(IZM)) )   B07900
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
               ZF1L = (FLOAT(JMIN1-2)-ZINT)*ZSLOPE                        B08100
               ZF2L = (FLOAT(JMIN2-2)-ZINT*CONF2)*ZSLOPE                  B08110
               ZF3L = (FLOAT(JMIN3-2)-ZINT*CONF3)*ZSLOPE                  B08120
               ZF1 = ZF1L                                                 B08130
               ZF2 = ZF2L                                                 B08140
               ZF3 = ZF3L                                                 B08150
               DO 10 J1 = JMIN1, JMAX1                                    B08160
                  J2 = J1-J2SHFT                                          B08170
                  J3 = J1-J3SHFT                                          B08180
                  ZF3 = ZF3+ZSLOPE                                        B08190
                  ZF2 = ZF2+ZSLOPE                                        B08200
                  ZF1 = ZF1+ZSLOPE                                        B08210
                  IZ3 = ABS(ZF3)+1.5                                      B08220
                  IZ2 = ABS(ZF2)+1.5                                      B08230
                  IZ1 = ABS(ZF1)+1.5                                      B08240
                  R3(J3) = R3(J3)+STRF3*F3(IZ3)                           B08250
                  R2(J2) = R2(J2)+STRF2*F2(IZ2)                           B08260
                  R1(J1) = R1(J1)+STRF1*F1(IZ1)+STRD*FG(IZ1)+STRVER*XVER  B08270
     *               (IZ1)                                                B08280
 
c                  if (iz1 .ge. 2601 .or. iz2 .ge. 2601 .or.
c     &                iz3 .ge. 2601 .or. iz1 .le. 0    .or. 
c     &                iz2 .le. 0    .or. iz3 .le. 0) then

c                 if (r3(j3) .lt. 0.0 .or. r2(j2) .lt. 0.0 .or.
                  if (r1(j1) .lt. 0.0) then
                      xxx=0
                  endif

   10          CONTINUE                                                   B08290
C                                                                         B08300
               IF (SPPSP(I).NE.0.) THEN                                   B08310
C                                                                         B08320
C                 THE FOLLOWING DOES LINE COUPLING                        B08330
C                                                                         B08340
C                 SPPSP(I) = SPP(I)/SP(I)                                 B08350
C                                                                         B08360
                  DPTRAT = SPPSP(I)                                       B08370
                  STRF3 = STRF3*CLC3*DPTRAT                               B08380
                  STRF2 = STRF2*CLC2*DPTRAT                               B08390
                  STRF1 = STRF1*CLC1*DPTRAT                               B08400
                  STRD = STRD*CLC1*DPTRAT                                 B08410
                  STRVER = STRVER*CLC1*DPTRAT                             B08420
C                                                                         B08430
C                                                                         B08440
                  DO 20 J1 = JMIN1, JMAX1                                 B08450
                     J2 = J1-J2SHFT                                       B08460
                     J3 = J1-J3SHFT                                       B08470
                     ZF3L = ZF3L+ZSLOPE                                   B08480
                     ZF2L = ZF2L+ZSLOPE                                   B08490
                     ZF1L = ZF1L+ZSLOPE                                   B08500
                     IZ3 = ABS(ZF3L)+1.5                                  B08510
                     IZ2 = ABS(ZF2L)+1.5                                  B08520
                     IZ1 = ABS(ZF1L)+1.5                                  B08530
                     R3(J3) = R3(J3)+STRF3*F3(IZ3)*ZF3L                   B08540
                     R2(J2) = R2(J2)+STRF2*F2(IZ2)*ZF2L                   B08550
                     R1(J1) = R1(J1)+(STRF1*F1(IZ1)+STRD*FG(IZ1)+STRVER*  B08560
     *                  XVER(IZ1))*ZF1L                                   B08570
                     if (r1(j1) .lt. 0.0) then
                         xxx=0
                     endif
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
      do 45 ii=1, 2601
          if (r1(ii) .lt. 0) then
              xxx = 0
          endif
  45  continue

      CALL CPUTIM (TIME)                                                  B08750
      TIMCNV = TIMCNV+TIME-TIME0                                          B08760
      RETURN                                                              B08770
C                                                                         B08780
      END                                                                 B08790
      SUBROUTINE PANEL (R1,R2,R3,KFILE,JRAD,IENTER)                       B09990
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
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           B10390
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
C     (in case of IOD=1,3), call PNLINT.
C
      IF (DVOUT.EQ.0.) THEN                                               B11300
         CALL BUFOUT (KFILE,PNLHDR(1),NPHDRF)                             B11310
         CALL BUFOUT (KFILE,R1(NLO),NLIM)                                 B11320
C                                                                         B11330
         IF (NPTS.GT.0) CALL R1PRNT (V1P,DVP,NLIM,R1,NLO,NPTS,KFILE,
     *                               IENTER)                              B11340
      ELSE                                                                B11350
         CALL PNLINT (R1(NLO),IENTER)                                     B11360
      ENDIF                                                               B11370
C                                                                         B11380
      VFT = VFT+FLOAT(NLIM1-1)*DV                                         B11390
      IF (ISTOP.NE.1) THEN                                                B11400
         DO 50 J = NLIM1, MAX1                                            B11420
            R1(J-NLIM1+1) = R1(J)                                         B11430
   50    CONTINUE                                                         B11450
         DO 60 J = MAX1-NLIM1+2, MAX1                                     B11460
            R1(J) = 0.                                                    B11470
   60    CONTINUE                                                         B11480
         DO 70 J = NLIM2, MAX2                                            B11500
            R2(J-NLIM2+1) = R2(J)                                         B11510
   70    CONTINUE                                                         B11530
         DO 80 J = MAX2-NLIM2+2, MAX2                                     B11540
            R2(J) = 0.                                                    B11550
   80    CONTINUE                                                         B11560
         DO 90 J = NLIM3, MAX3                                            B11580
            R3(J-NLIM3+1) = R3(J)                                         B11590
   90    CONTINUE                                                         B11610
         DO 100 J = MAX3-NLIM3+2, MAX3                                    B11620
            R3(J) = 0.                                                    B11630
  100    CONTINUE                                                         B11640
         NLO = NSHIFT+1                                                   B11650
      ENDIF                                                               B11660
      CALL CPUTIM (TIME)                                                  B11670
      TIMPNL = TIMPNL+TIME-TIME0                                          B11680
C                                                                         B11690
      RETURN                                                              B11700
C                                                                         B11710
      END                                                                 B11720
      SUBROUTINE PNLINT (R1,IENTER)                                       B11730
C                                                                         B11740
      IMPLICIT REAL*8           (V)                                     ! B11750
C                                                                         B11760
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B11770
C                                                                         B11780
C               LAST MODIFICATION:    6 May 1994 pdb                      B11790
C               LAST MODIFICATION:    9 APRIL 1991                        B11790
C                                                                         B11800
C                  IMPLEMENTATION:    R.D. WORSHAM                        B11810
C                                                                         B11820
C                       ALGORITHM:    R.D. WORSHAM                        B11830
C                                     S.A. CLOUGH                         B11840
C                                     J.L. MONCET                         B11850
C                                                                         B11860
C                                                                         B11870
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B11880
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          B11890
C                                                                         B11900
C----------------------------------------------------------------------   B11910
C                                                                         B11920
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     B11930
C                                     OFFICE OF ENERGY RESEARCH           B11940
C                                     DEPARTMENT OF ENERGY                B11950
C                                                                         B11960
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B11970
C                                                                         B11980
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ
C                                                                         B12000
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       B12010
     *                WK(60),PZ1,PZ2,TZ1,TZ2,WBROAD,DV ,V1 ,V2 ,TBOUND,   B12020
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    B12030
      COMMON /XSUB/ VBOT,VTOP,VFT,LIMIN,ILO,IHI,IEOF,IPANEL,ISTOP,IDATA
      COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTS       B12040
      COMMON /XPANO/ V1PO,V2PO,DVPO,NLIM2,RMINO,RMAXO,NPNXPO              B12050
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         B12060
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        B12070
     *              NLTEFL,LNFIL4,LNGTH4                                  B12080
      COMMON /R1SAV/ R1OUT(2410)                                          B12090
      COMMON /IODFLG/ DVOUT
C
C     SAVE statement to preserve value of NLIM1 when returning to
C     subroutine
C
      SAVE NLIM1
C                                                                         B12100
      DIMENSION A1(0:100),A2(0:100),A3(0:100),A4(0:100)                   B12110
      DIMENSION R1(*)                                                     B12120
      DIMENSION PNLHDR(2),PNLHDO(2)                                       B12130
C                                                                         B12140
      EQUIVALENCE (PNLHDR(1),V1P),(PNLHDO(1),V1PO)                        B12150
C                                                                         B12160
      DATA LIMOUT / 2400 /                                                B12180
C                                                                         B12190
C     The data for NM1 and N0 are used instead of directly inserting
C     '-1' and '0' into the subscripts for R1 (lines 1158-9) to avoid     B12200
C     compiler warnings 'CONSTANT SUBSCRIPT IS OUT OF BOUNDS'
C     
      DATA NM1/-1/,N0/0/                                                  B12210
C
C      CALL CPUTIM (TIME)                                                 B12220
C      WRITE (IPR,900) TIME                                               B12230
C                                                                         B12240
C     The value of DVOUT is carried from COMMON BLOCK /IODFLG/
C
      NPNXPO = NPNLXP                                                     B12250
      DVPO = DVOUT                                                        B12260
      NPPANL = 1                                                          B12270
      ATYPE = 9.999E09                                                    B12280
      IF (DVP.EQ.DVOUT) ATYPE = 0.                                        B12290
      IF (DVOUT.GT.DVP) ATYPE = DVP/(DVOUT-DVP)+0.5                       B12300
      IF (DVOUT.LT.DVP) ATYPE = -DVOUT/(DVP-DVOUT)-0.5                    B12310
      IF (ABS(DVOUT-DVP).LT.1.E-8) ATYPE = 0.                             B12320
C                                                                         B12330
C                                                                         B12340
C     1/1 ratio only                                                      B12350
C                                                                         B12360
      IF (ATYPE.EQ.0.) THEN                                               B12370
         CALL PMNMX (R1,NLIM,RMIN,RMAX)                                   B12380
         CALL BUFOUT (KFILE,PNLHDR(1),NPHDRF)                             B12390
         CALL BUFOUT (KFILE,R1(1),NLIM)                                   B12400
C                                                                         B12410
         IF (NPTS.GT.0) CALL R1PRNT (V1P,DVP,NLIM,R1,1,NPTS,KFILE,
     *                               IENTER)                              B12420
C                                                                         B12430
         GO TO 40                                                         B12440
      ENDIF                                                               B12460
C                                                                         B12470
C     All ratios except 1/1                                               B12480
C                                                                         B12490
      DO 10 JP = 0, 100                                                   B12500
         APG = JP                                                         B12510
         P = 0.01*APG                                                     B12520
C                                                                         B12540
C        Constants for the Lagrange 4 point interpolation                 B12560
C                                                                         B12570
         A1(JP) = -P*(P-1.0)*(P-2.0)/6.0                                  B12580
         A2(JP) = (P**2-1.0)*(P-2.0)*0.5                                  B12590
         A3(JP) = -P*(P+1.0)*(P-2.0)*0.5                                  B12600
         A4(JP) = P*(P**2-1.0)/6.0                                        B12610
   10 CONTINUE                                                            B12620
C                                                                         B12630
C     Zero point of first panel                                           B12680
C                                                                         B12690
      IF (V1PO.EQ.0.0) THEN                                               B12700
         R1(NM1) = R1(1)                                                  B12710
         R1(N0) = R1(1)                                                   B12720
         V1PO = V1P                                                       B12730
         NLIM1 = 1                                                        B12740
      ENDIF                                                               B12750
C                                                                         B12760
C     Add points to end of last panel for interpolation                   B12770
C                                                                         B12780
      IF (ISTOP.EQ.1) THEN                                                B12790
         R1(NLIM+1) = R1(NLIM)                                            B12800
         R1(NLIM+2) = R1(NLIM)                                            B12810
         NLIM = NLIM + 2                                                  B12820
         V2P = V2P + 2.*DVP                                               B12830
      ENDIF                                                               B12840
C                                                                         B12850
C     *** BEGINNING OF LOOP THAT DOES INTERPOLATION  ***
C
   20 CONTINUE
C
C     Determine potential last point for the outgoing panel (2400 pts.)
C
      V2PO = V1PO+FLOAT(LIMOUT)*DVOUT                                     B12860
C                                                                         B12870
      IF (V2P.LE.V2PO+DVP.AND.ILAST.EQ.0.AND.NPPANL.LE.0) GO TO 40        B12880
C                                                                         B12890
C     Four possibilities:
C       1a.  Last panel to be done, set the appropriate
C            final output point and total number of points in panel.
C
C       1b.  Would be last panel, but need more incoming points to
C            fill panel.
C
C       2a.  More panels to come, set last point in panel.
C
C       2b.  More panels to come, but need more incoming points to
C            fill panel.
C
      IF ((V1PO+(LIMOUT-1)*DVOUT).GT.V2) THEN                             B12900
         NLIM2 = (V2-V1PO)/DVOUT + 1.                                     B12910
         V2PO = V1PO+FLOAT(NLIM2-1)*DVOUT                                 B12920
         IF (V2PO.LT.V2) THEN
            V2PO = V2PO+DVOUT
            NLIM2 = NLIM2+1
         ENDIF
         ILAST = 1
         IF (V2PO.GT.V2P-DVP) THEN
            NLIM2 = ((V2P-DVP-V1PO)/DVOUT) + 1.
            V2PO = V1PO+FLOAT(NLIM2-1)*DVOUT
            IF (V2PO+DVOUT.LT.V2P-DVP) THEN
               NLIM2 = NLIM2+1
               V2PO = V2PO+DVOUT
            ENDIF
            ILAST = 0
         ENDIF
      ELSE
         NLIM2 = LIMOUT
         V2PO = V1PO+FLOAT(NLIM2-1)*DVOUT                                 B12930
         IF (V2PO.GT.V2P-DVP) THEN                                        B12940
            NLIM2 = ((V2P-DVP-V1PO)/DVOUT) + 1.                           B12950
            V2PO = V1PO+FLOAT(NLIM2-1)*DVOUT                              B12960
            IF (V2PO+DVOUT.LT.V2P-DVP) THEN                               B12970
               NLIM2 = NLIM2+1                                            B12980
               V2PO = V2PO+DVOUT                                          B12990
            ENDIF                                                         B13000
         ENDIF                                                            B13010
         ILAST = 0                                                        B13020
      ENDIF                                                               B13030
C
      RATDV = DVOUT/DVP                                                   B13040
C                                                                         B13050
C     FJJ is offset by 2. for rounding purposes                           B13060
C                                                                         B13070
      FJ1DIF = (V1PO-V1P)/DVP+1.+2.                                       B13080
C                                                                         B13090
C     Interpolate R1 to DVOUT                                             B13100
C                                                                         B13110
      DO 30 II = NLIM1, NLIM2                                             B13120
         FJJ = FJ1DIF+RATDV*FLOAT(II-1)                                   B13130
         JJ  = IFIX(FJJ)-2                                                B13140
         JP  = (FJJ-FLOAT(JJ))*100.-199.5                                 B13150
         R1OUT(II) = A1(JP)*R1(JJ-1)+A2(JP)*R1(JJ)+A3(JP)*R1(JJ+1)+       B13160
     *               A4(JP)*R1(JJ+2)                                      B13170
   30 CONTINUE                                                            B13180
C                                                                         B13190
C     Two possibilities:
C       1.  Buffer out whole panel (NLIM2 = 2400) or the remaining        B13200
C           interpolated points                                           B13210
C
C       2.  Return to PANEL to obtain more incoming points to fill
C           outgoing panel
C
      IF (NLIM2.EQ.LIMOUT.OR.ILAST.EQ.1) THEN                             B13220
         CALL PMNMX (R1OUT,NLIM2,RMINO,RMAXO)                             B13230
         CALL BUFOUT (KFILE,PNLHDO(1),NPHDRF)                             B13240
         CALL BUFOUT (KFILE,R1OUT(1),NLIM2)                               B13250
         IF (NPTS.GT.0) CALL R1PRNT (V1PO,DVOUT,NLIM2,R1OUT,1,NPTS,
     *        KFILE,IENTER)                                               B13260
         NLIM1 = 1                                                        B13270
         NPPANL = 0                                                       B13280
         V1PO = V2PO+DVOUT                                                B13290
         IF ((V1PO+FLOAT(LIMOUT)*DVOUT).GT.(V2P-DVP)) NPPANL = 1          B13300
      ELSE                                                                B13310
         NLIM1 = NLIM2+1                                                  B13320
         NPPANL = -1                                                      B13330
      ENDIF                                                               B13340
C                                                                         B13350
C     If not at last point, continue interpolation                        B13360
C                                                                         B13370
      IF (ILAST.NE.1) GO TO 20                                            B13380
C                                                                         B13390
C     Reset variables                                                     B13400
C                                                                         B13410
      V1PO = 0.0                                                          B13440
      NPPANL = 1                                                          B13450
   40 CONTINUE                                                            B13460
C                                                                         B13470
C      CALL CPUTIM (TIME1)                                                B13480
C      TIM = TIME1-TIME                                                   B13490
C      WRITE (IPR,905) TIME1,TIM                                          B13500
C                                                                         B13510
      RETURN                                                              B13520
C                                                                         B13530
C  900 FORMAT ('0 THE TIME AT THE START OF PNLINT IS ',F12.3)             B13540
C  905 FORMAT ('0 THE TIME AT THE END OF PNLINT IS ',F12.3/F12.3,         B13550
C     *   ' SECS WERE REQUIRED FOR THIS INTERPOLATION ')                  B13560
C                                                                         B13570
      END                                                                 B13580
      SUBROUTINE PMNMX (R1,NLIM,RMIN,RMAX)                                B13590
C                                                                         B13600
      DIMENSION R1(NLIM)                                                  B13610
C                                                                         B13620
      RMIN = R1(1)                                                        B13630
      RMAX = R1(1)                                                        B13640
C                                                                         B13650
      DO 10 I = 2, NLIM                                                   B13660
         RMIN = MIN(RMIN,R1(I))                                           B13670
         RMAX = MAX(RMAX,R1(I))                                           B13680
   10 CONTINUE                                                            B13690
C                                                                         B13700
      RETURN                                                              B13710
C                                                                         B13720
      END                                                                 B13730
      SUBROUTINE SHAPEL (F1,F2,F3)                                        B13740
C                                                                         B13750
C     SUBROUTINE SHAPEL CONSTRUCTS THE SUB-FUNCTIONS FOR THE              B13760
C     LORENTZ LINE SHAPE                                                  B13770
C                                                                         B13780
      COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,            B13790
     *                HWF3,DXF3,NX3,N3MAX                                 B13800
      DIMENSION F1(*),F2(*),F3(*)                                         B13810
C                                                                         B13820
      XLORNZ(XSQ) = 1./(1.+XSQ)                                           B13830
      Q1FN(XSQ) = A1+B1*XSQ                                               B13840
      Q2FN(XSQ) = A2+B2*XSQ                                               B13850
      Q3FN(XSQ) = A3+B3*XSQ                                               B13860
C                                                                         B13870
      A(Z0) = (1.+2.*Z0*Z0)/(1.+Z0*Z0)**2                                 B13880
      B(Z0) = -1./(1.+Z0*Z0)**2                                           B13890
      RECPI = 1./(2.*ASIN(1.))                                            B13900
      TOTAL = 0.                                                          B13910
      A1 = A(HWF1)                                                        B13920
      B1 = B(HWF1)                                                        B13930
C                                                                         B13940
      A2 = A(HWF2)                                                        B13950
      B2 = B(HWF2)                                                        B13960
C                                                                         B13970
      A3 = A(HWF3)                                                        B13980
      B3 = B(HWF3)                                                        B13990
C                                                                         B14000
      DO 10 I = 1, N1MAX                                                  B14010
         F1(I) = 0.                                                       B14020
   10 CONTINUE                                                            B14030
      F1(1) = RECPI*(XLORNZ(0.)-Q1FN(0.))                                 B14040
      SUM = F1(1)                                                         B14050
      DO 20 JJ = 2, NX1                                                   B14060
         X = FLOAT(JJ-1)*DXF1                                             B14070
         XSQ = X*X                                                        B14080
         F1(JJ) = RECPI*(XLORNZ(XSQ)-Q1FN(XSQ))                           B14090
         SUM = SUM+F1(JJ)*2.                                              B14100
   20 CONTINUE                                                            B14110
      F1(NX1) = 0.                                                        B14120
      SUM = SUM*DXF1                                                      B14130
      TOTAL = TOTAL+SUM                                                   B14140
C                                                                         B14150
      DO 30 I = 1, N2MAX                                                  B14160
         F2(I) = 0.                                                       B14170
   30 CONTINUE                                                            B14180
      F2(1) = RECPI*(Q1FN(0.)-Q2FN(0.))                                   B14190
      SUM = F2(1)                                                         B14200
      J1LIM = HWF1/DXF2+1.001                                             B14210
      DO 40 JJ = 2, J1LIM                                                 B14220
         X = FLOAT(JJ-1)*DXF2                                             B14230
         XSQ = X*X                                                        B14240
         F2(JJ) = RECPI*(Q1FN(XSQ)-Q2FN(XSQ))                             B14250
         SUM = SUM+F2(JJ)*2.                                              B14260
   40 CONTINUE                                                            B14270
      J1LIMP = J1LIM+1                                                    B14280
      DO 50 JJ = J1LIMP, NX2                                              B14290
         X = FLOAT(JJ-1)*DXF2                                             B14300
         XSQ = X*X                                                        B14310
         F2(JJ) = RECPI*(XLORNZ(XSQ)-Q2FN(XSQ))                           B14320
         SUM = SUM+F2(JJ)*2.                                              B14330
   50 CONTINUE                                                            B14340
      F2(NX2) = 0.                                                        B14350
      SUM = SUM*DXF2                                                      B14360
      TOTAL = TOTAL+SUM                                                   B14370
C                                                                         B14380
      DO 60 I = 1, N3MAX                                                  B14390
         F3(I) = 0.                                                       B14400
   60 CONTINUE                                                            B14410
      F3(1) = RECPI*(Q2FN(0.)-Q3FN(0.))                                   B14420
      SUM = F3(1)                                                         B14430
      J2LIM = HWF2/DXF3+1.001                                             B14440
      DO 70 JJ = 2, J2LIM                                                 B14450
         X = FLOAT(JJ-1)*DXF3                                             B14460
         XSQ = X*X                                                        B14470
         F3(JJ) = RECPI*(Q2FN(XSQ)-Q3FN(XSQ))                             B14480
         SUM = SUM+F3(JJ)*2.                                              B14490
   70 CONTINUE                                                            B14500
      J2LIMP = J2LIM+1                                                    B14510
      DO 80 JJ = J2LIMP, NX3                                              B14520
         X = FLOAT(JJ-1)*DXF3                                             B14530
         XSQ = X*X                                                        B14540
         F3(JJ) = RECPI*(XLORNZ(XSQ)-Q3FN(XSQ))                           B14550
         SUM = SUM+F3(JJ)*2.                                              B14560
   80 CONTINUE                                                            B14570
      SUM = SUM*DXF3                                                      B14580
      TOTAL = TOTAL+SUM                                                   B14590
C                                                                         B14600
      RETURN                                                              B14610
C                                                                         B14620
      END                                                                 B14630
      SUBROUTINE SHAPEG (FG)                                              B14640
C                                                                         B14650
C     SUBROUTINE SHAPEG CONSTRUCTS THE FUNCTION FOR THE DOPPLER PROFILE   B14660
C                                                                         B14670
      COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,            B14680
     *                HWF3,DXF3,NX3,N3MAX                                 B14690
      DIMENSION FG(*)                                                     B14700
C                                                                         B14710
      FGAUSS(XSQ) = EXP(-FLN2*XSQ)                                        B14720
      FLN2 = ALOG(2.)                                                     B14730
      RECPI = 1./(2.*ASIN(1.))                                            B14740
      FGNORM = SQRT(FLN2*RECPI)                                           B14750
      TOTAL = 0.                                                          B14760
      DO 10 I = 1, N1MAX                                                  B14770
         FG(I) = 0.                                                       B14780
   10 CONTINUE                                                            B14790
      FG(1) = FGNORM*FGAUSS(0.)                                           B14800
      SUM = FG(1)                                                         B14810
      DO 20 JJ = 2, NX1                                                   B14820
         X = FLOAT(JJ-1)*DXF1                                             B14830
         XSQ = X*X                                                        B14840
         FG(JJ) = FGNORM*FGAUSS(XSQ)                                      B14850
         SUM = SUM+FG(JJ)*2.                                              B14860
   20 CONTINUE                                                            B14870
      FG(NX1) = 0.                                                        B14880
      SUM = SUM*DXF1                                                      B14890
      TOTAL = TOTAL+SUM                                                   B14900
C                                                                         B14910
      RETURN                                                              B14920
C                                                                         B14930
      END                                                                 B14940
      SUBROUTINE VERFN (XVER)                                             B14950
C                                                                         B14960
C     VERFN IS A FUNCTION USED TO IMPROVE THE ACCURACY OF THE             B14970
C     VOIGT APPROXIMATION IN THE DOMAIN 0 - 4 HALFWIDTHS.                 B14980
C                                                                         B14990
      COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,            B15000
     *                HWF3,DXF3,NX3,N3MAX                                 B15010
      DIMENSION XVER(*)                                                   B15020
C                                                                         B15030
C     FOR ZETA = 0.3                                                      B15040
C                                                                         B15050
      DATA CEXP,CE0,CE2,CE4 / 0.45,1.,-.20737285249,-.00872684335747 /    B15060
      DATA SUM0,SUM2,SUM4,SUMER / 4*0. /                                  B15070
C                                                                         B15080
      ERFN(Z2) = (1./(CE0+CE2+CE4))*(CE0+CE2*AE2*Z2+CE4*AE4*Z2*Z2)*XE0    B15090
C                                                                         B15100
      IF (SUMER.NE.0.) RETURN                                             B15110
      PI = 2.*ASIN(1.)                                                    B15120
      SE0 = SQRT(CEXP/PI)                                                 B15130
      AE0 = 1.                                                            B15140
      AE2 = 2.*CEXP                                                       B15150
      AE4 = AE2*AE2/3.                                                    B15160
      FACTOR = 1.                                                         B15170
C                                                                         B15180
      DO 10 I = 1, N1MAX                                                  B15190
         XVER(I) = 0.                                                     B15200
   10 CONTINUE                                                            B15210
C                                                                         B15220
      DO 20 I = 1, N1MAX                                                  B15230
         Z = DXF1*FLOAT(I-1)                                              B15240
         Z2 = Z*Z                                                         B15250
         XE0 = SE0*EXP(-CEXP*Z2)                                          B15260
         XE2 = AE2*Z2*XE0                                                 B15270
         XE4 = AE4*Z2*Z2*XE0                                              B15280
         XVER(I) = ERFN(Z2)                                               B15290
         SUM0 = SUM0+FACTOR*DXF1*XE0                                      B15300
         SUM2 = SUM2+FACTOR*DXF1*XE2                                      B15310
         SUM4 = SUM4+FACTOR*DXF1*XE4                                      B15320
         SUMER = SUMER+FACTOR*DXF1*XVER(I)                                B15330
         FACTOR = 2.                                                      B15340
   20 CONTINUE                                                            B15350
C                                                                         B15360
CPRT  WRITE (IPR,900) Z,SUM0,SUM2,SUM4,SUMER                              B15370
C                                                                         B15380
      RETURN                                                              B15390
C                                                                         B15400
  900 FORMAT (F10.3,6F15.10)                                              B15410
C                                                                         B15420
      END                                                                 B15430
      BLOCK DATA VOICON                                                   B15440
C                                                                         B15450
C     AVRAT CONTAINS THE PARAMTERS AS A FUNCTION OF ZETA USED TO          B15460
C     OBTAIN THE VOIGTS' WIDTH FROM THE LORENTZ AND DOPPLER WIDTHS.       B15470
C                                                                         B15480
C     COMMON /VOICOM/ AVRAT(102),                                         B15490
C    C                CGAUSS(102),CF1(102),CF2(102),CF3(102),CER(102)     B15500
C                                                                         B15510
      COMMON /VOICOM/ AV01(50),AV51(52),CG01(50),CG51(52),CFA01(50),      B15520
     *                CFA51(52),CFB01(50),CFB51(52),CFC01(50),            B15530
     *                CFC51(52),CER01(50),CER51(52)                       B15540
C                                                                         B15550
       DATA AV01/                                                         B15560
     *  .10000E+01,  .99535E+00,  .99073E+00,  .98613E+00,  .98155E+00,   B15570
     *  .97700E+00,  .97247E+00,  .96797E+00,  .96350E+00,  .95905E+00,   B15580
     *  .95464E+00,  .95025E+00,  .94589E+00,  .94156E+00,  .93727E+00,   B15590
     *  .93301E+00,  .92879E+00,  .92460E+00,  .92045E+00,  .91634E+00,   B15600
     *  .91227E+00,  .90824E+00,  .90425E+00,  .90031E+00,  .89641E+00,   B15610
     *  .89256E+00,  .88876E+00,  .88501E+00,  .88132E+00,  .87768E+00,   B15620
     *  .87410E+00,  .87058E+00,  .86712E+00,  .86372E+00,  .86039E+00,   B15630
     *  .85713E+00,  .85395E+00,  .85083E+00,  .84780E+00,  .84484E+00,   B15640
     *  .84197E+00,  .83919E+00,  .83650E+00,  .83390E+00,  .83141E+00,   B15650
     *  .82901E+00,  .82672E+00,  .82454E+00,  .82248E+00,  .82053E+00/   B15660
       DATA AV51/                                                         B15670
     *  .81871E+00,  .81702E+00,  .81547E+00,  .81405E+00,  .81278E+00,   B15680
     *  .81166E+00,  .81069E+00,  .80989E+00,  .80925E+00,  .80879E+00,   B15690
     *  .80851E+00,  .80842E+00,  .80852E+00,  .80882E+00,  .80932E+00,   B15700
     *  .81004E+00,  .81098E+00,  .81214E+00,  .81353E+00,  .81516E+00,   B15710
     *  .81704E+00,  .81916E+00,  .82154E+00,  .82418E+00,  .82708E+00,   B15720
     *  .83025E+00,  .83370E+00,  .83742E+00,  .84143E+00,  .84572E+00,   B15730
     *  .85029E+00,  .85515E+00,  .86030E+00,  .86573E+00,  .87146E+00,   B15740
     *  .87747E+00,  .88376E+00,  .89035E+00,  .89721E+00,  .90435E+00,   B15750
     *  .91176E+00,  .91945E+00,  .92741E+00,  .93562E+00,  .94409E+00,   B15760
     *  .95282E+00,  .96179E+00,  .97100E+00,  .98044E+00,  .99011E+00,   B15770
     *  .10000E+01,  .10000E+01/                                          B15780
      DATA CG01 /                                                         B15790
     *  1.00000E+00, 1.01822E+00, 1.03376E+00, 1.04777E+00, 1.06057E+00,  B15800
     *  1.07231E+00, 1.08310E+00, 1.09300E+00, 1.10204E+00, 1.11025E+00,  B15810
     *  1.11766E+00, 1.12429E+00, 1.13014E+00, 1.13523E+00, 1.13955E+00,  B15820
     *  1.14313E+00, 1.14595E+00, 1.14803E+00, 1.14936E+00, 1.14994E+00,  B15830
     *  1.14978E+00, 1.14888E+00, 1.14723E+00, 1.14484E+00, 1.14170E+00,  B15840
     *  1.13782E+00, 1.13319E+00, 1.12782E+00, 1.12171E+00, 1.11485E+00,  B15850
     *  1.10726E+00, 1.09893E+00, 1.08986E+00, 1.08006E+00, 1.06953E+00,  B15860
     *  1.05828E+00, 1.04631E+00, 1.03363E+00, 1.02024E+00, 1.00617E+00,  B15870
     *  9.91403E-01, 9.75966E-01, 9.59868E-01, 9.43123E-01, 9.25745E-01,  B15880
     *  9.07752E-01, 8.89162E-01, 8.69994E-01, 8.50272E-01, 8.30017E-01/  B15890
      DATA CG51 /                                                         B15900
     *  8.09256E-01, 7.88017E-01, 7.66327E-01, 7.44219E-01, 7.21726E-01,  B15910
     *  6.98886E-01, 6.75729E-01, 6.52299E-01, 6.28637E-01, 6.04787E-01,  B15920
     *  5.80794E-01, 5.56704E-01, 5.32565E-01, 5.08428E-01, 4.84343E-01,  B15930
     *  4.60364E-01, 4.36543E-01, 4.12933E-01, 3.89589E-01, 3.66564E-01,  B15940
     *  3.43913E-01, 3.21688E-01, 2.99940E-01, 2.78720E-01, 2.58077E-01,  B15950
     *  2.38056E-01, 2.18701E-01, 2.00053E-01, 1.82148E-01, 1.65021E-01,  B15960
     *  1.48701E-01, 1.33213E-01, 1.18579E-01, 1.04815E-01, 9.19338E-02,  B15970
     *  7.99428E-02, 6.88453E-02, 5.86399E-02, 4.93211E-02, 4.08796E-02,  B15980
     *  3.33018E-02, 2.65710E-02, 2.06669E-02, 1.55667E-02, 1.12449E-02,  B15990
     *  7.67360E-03, 4.82345E-03, 2.66344E-03, 1.16151E-03, 2.84798E-04,  B16000
     *  0.         , 0.         /                                         B16010
      DATA CFA01 /                                                        B16020
     *  0.         ,-2.56288E-03,-3.05202E-03,-2.50689E-03,-1.18504E-03,  B16030
     *  7.84668E-04, 3.32528E-03, 6.38605E-03, 9.93124E-03, 1.39345E-02,  B16040
     *  1.83758E-02, 2.32392E-02, 2.85120E-02, 3.41837E-02, 4.02454E-02,  B16050
     *  4.66897E-02, 5.35099E-02, 6.07003E-02, 6.82556E-02, 7.61711E-02,  B16060
     *  8.44422E-02, 9.30647E-02, 1.02034E-01, 1.11348E-01, 1.21000E-01,  B16070
     *  1.30988E-01, 1.41307E-01, 1.51952E-01, 1.62921E-01, 1.74208E-01,  B16080
     *  1.85808E-01, 1.97716E-01, 2.09927E-01, 2.22436E-01, 2.35236E-01,  B16090
     *  2.48321E-01, 2.61684E-01, 2.75318E-01, 2.89215E-01, 3.03367E-01,  B16100
     *  3.17764E-01, 3.32399E-01, 3.47260E-01, 3.62338E-01, 3.77620E-01,  B16110
     *  3.93096E-01, 4.08752E-01, 4.24575E-01, 4.40550E-01, 4.56665E-01/  B16120
      DATA CFA51 /                                                        B16130
     *  4.72901E-01, 4.89244E-01, 5.05675E-01, 5.22177E-01, 5.38731E-01,  B16140
     *  5.55315E-01, 5.71913E-01, 5.88502E-01, 6.05059E-01, 6.21561E-01,  B16150
     *  6.37986E-01, 6.54308E-01, 6.70504E-01, 6.86549E-01, 7.02417E-01,  B16160
     *  7.18083E-01, 7.33520E-01, 7.48703E-01, 7.63607E-01, 7.78204E-01,  B16170
     *  7.92472E-01, 8.06384E-01, 8.19918E-01, 8.33050E-01, 8.45759E-01,  B16180
     *  8.58025E-01, 8.69828E-01, 8.81151E-01, 8.91979E-01, 9.02298E-01,  B16190
     *  9.12097E-01, 9.21366E-01, 9.30098E-01, 9.38289E-01, 9.45935E-01,  B16200
     *  9.53036E-01, 9.59594E-01, 9.65613E-01, 9.71101E-01, 9.76064E-01,  B16210
     *  9.80513E-01, 9.84460E-01, 9.87919E-01, 9.90904E-01, 9.93432E-01,  B16220
     *  9.95519E-01, 9.97184E-01, 9.98445E-01, 9.99322E-01, 9.99834E-01,  B16230
     *  1.00000E+00, 1.00000E+00/                                         B16240
      DATA CFB01 /                                                        B16250
     *  0.         , 1.15907E-02, 2.32978E-02, 3.51022E-02, 4.69967E-02,  B16260
     *  5.89773E-02, 7.10411E-02, 8.31858E-02, 9.54097E-02, 1.07711E-01,  B16270
     *  1.20089E-01, 1.32541E-01, 1.45066E-01, 1.57663E-01, 1.70330E-01,  B16280
     *  1.83065E-01, 1.95868E-01, 2.08737E-01, 2.21669E-01, 2.34664E-01,  B16290
     *  2.47718E-01, 2.60830E-01, 2.73998E-01, 2.87219E-01, 3.00491E-01,  B16300
     *  3.13812E-01, 3.27178E-01, 3.40587E-01, 3.54035E-01, 3.67520E-01,  B16310
     *  3.81037E-01, 3.94583E-01, 4.08155E-01, 4.21747E-01, 4.35356E-01,  B16320
     *  4.48978E-01, 4.62606E-01, 4.76237E-01, 4.89864E-01, 5.03482E-01,  B16330
     *  5.17086E-01, 5.30669E-01, 5.44225E-01, 5.57746E-01, 5.71226E-01,  B16340
     *  5.84657E-01, 5.98032E-01, 6.11342E-01, 6.24580E-01, 6.37736E-01/  B16350
      DATA CFB51 /                                                        B16360
     *  6.50802E-01, 6.63769E-01, 6.76626E-01, 6.89365E-01, 7.01974E-01,  B16370
     *  7.14444E-01, 7.26764E-01, 7.38924E-01, 7.50912E-01, 7.62717E-01,  B16380
     *  7.74328E-01, 7.85735E-01, 7.96925E-01, 8.07888E-01, 8.18612E-01,  B16390
     *  8.29087E-01, 8.39302E-01, 8.49246E-01, 8.58910E-01, 8.68284E-01,  B16400
     *  8.77358E-01, 8.86125E-01, 8.94577E-01, 9.02706E-01, 9.10506E-01,  B16410
     *  9.17972E-01, 9.25100E-01, 9.31885E-01, 9.38325E-01, 9.44419E-01,  B16420
     *  9.50166E-01, 9.55568E-01, 9.60625E-01, 9.65340E-01, 9.69718E-01,  B16430
     *  9.73763E-01, 9.77481E-01, 9.80878E-01, 9.83962E-01, 9.86741E-01,  B16440
     *  9.89223E-01, 9.91419E-01, 9.93337E-01, 9.94989E-01, 9.96385E-01,  B16450
     *  9.97536E-01, 9.98452E-01, 9.99146E-01, 9.99628E-01, 9.99909E-01,  B16460
     *  1.00000E+00, 1.00000E+00/                                         B16470
      DATA CFC01 /                                                        B16480
     *  0.         , 9.88700E-03, 1.98515E-02, 2.99036E-02, 4.00474E-02,  B16490
     *  5.02856E-02, 6.06200E-02, 7.10521E-02, 8.15830E-02, 9.22137E-02,  B16500
     *  1.02945E-01, 1.13778E-01, 1.24712E-01, 1.35749E-01, 1.46889E-01,  B16510
     *  1.58132E-01, 1.69478E-01, 1.80928E-01, 1.92480E-01, 2.04136E-01,  B16520
     *  2.15894E-01, 2.27754E-01, 2.39716E-01, 2.51780E-01, 2.63943E-01,  B16530
     *  2.76205E-01, 2.88564E-01, 3.01020E-01, 3.13571E-01, 3.26214E-01,  B16540
     *  3.38948E-01, 3.51771E-01, 3.64679E-01, 3.77670E-01, 3.90741E-01,  B16550
     *  4.03888E-01, 4.17108E-01, 4.30397E-01, 4.43750E-01, 4.57162E-01,  B16560
     *  4.70628E-01, 4.84142E-01, 4.97700E-01, 5.11293E-01, 5.24915E-01,  B16570
     *  5.38560E-01, 5.52218E-01, 5.65882E-01, 5.79542E-01, 5.93190E-01/  B16580
      DATA CFC51 /                                                        B16590
     *  6.06816E-01, 6.20408E-01, 6.33957E-01, 6.47451E-01, 6.60877E-01,  B16600
     *  6.74223E-01, 6.87477E-01, 7.00624E-01, 7.13651E-01, 7.26544E-01,  B16610
     *  7.39288E-01, 7.51868E-01, 7.64268E-01, 7.76474E-01, 7.88470E-01,  B16620
     *  8.00240E-01, 8.11768E-01, 8.23041E-01, 8.34042E-01, 8.44756E-01,  B16630
     *  8.55171E-01, 8.65271E-01, 8.75044E-01, 8.84478E-01, 8.93562E-01,  B16640
     *  9.02285E-01, 9.10639E-01, 9.18616E-01, 9.26210E-01, 9.33414E-01,  B16650
     *  9.40227E-01, 9.46644E-01, 9.52666E-01, 9.58293E-01, 9.63528E-01,  B16660
     *  9.68373E-01, 9.72833E-01, 9.76915E-01, 9.80625E-01, 9.83973E-01,  B16670
     *  9.86967E-01, 9.89617E-01, 9.91935E-01, 9.93933E-01, 9.95622E-01,  B16680
     *  9.97015E-01, 9.98125E-01, 9.98965E-01, 9.99549E-01, 9.99889E-01,  B16690
     *  1.00000E+00, 1.00000E+00/                                         B16700
      DATA CER01 /                                                        B16710
     *  0.         ,-2.11394E-02,-4.08818E-02,-5.97585E-02,-7.79266E-02,  B16720
     * -9.54663E-02,-1.12425E-01,-1.28834E-01,-1.44713E-01,-1.60076E-01,  B16730
     * -1.74933E-01,-1.89289E-01,-2.03149E-01,-2.16515E-01,-2.29388E-01,  B16740
     * -2.41768E-01,-2.53653E-01,-2.65043E-01,-2.75936E-01,-2.86328E-01,  B16750
     * -2.96217E-01,-3.05601E-01,-3.14476E-01,-3.22839E-01,-3.30686E-01,  B16760
     * -3.38015E-01,-3.44822E-01,-3.51105E-01,-3.56859E-01,-3.62083E-01,  B16770
     * -3.66773E-01,-3.70928E-01,-3.74546E-01,-3.77625E-01,-3.80164E-01,  B16780
     * -3.82161E-01,-3.83618E-01,-3.84534E-01,-3.84911E-01,-3.84749E-01,  B16790
     * -3.84051E-01,-3.82821E-01,-3.81062E-01,-3.78778E-01,-3.75976E-01,  B16800
     * -3.72663E-01,-3.68845E-01,-3.64532E-01,-3.59733E-01,-3.54461E-01/  B16810
      DATA CER51 /                                                        B16820
     * -3.48726E-01,-3.42543E-01,-3.35927E-01,-3.28893E-01,-3.21461E-01,  B16830
     * -3.13650E-01,-3.05477E-01,-2.96967E-01,-2.88142E-01,-2.79029E-01,  B16840
     * -2.69652E-01,-2.60040E-01,-2.50221E-01,-2.40225E-01,-2.30084E-01,  B16850
     * -2.19829E-01,-2.09493E-01,-1.99109E-01,-1.88712E-01,-1.78335E-01,  B16860
     * -1.68014E-01,-1.57782E-01,-1.47673E-01,-1.37721E-01,-1.27957E-01,  B16870
     * -1.18414E-01,-1.09120E-01,-1.00105E-01,-9.13939E-02,-8.30122E-02,  B16880
     * -7.49818E-02,-6.73226E-02,-6.00518E-02,-5.31840E-02,-4.67313E-02,  B16890
     * -4.07029E-02,-3.51053E-02,-2.99424E-02,-2.52153E-02,-2.09229E-02,  B16900
     * -1.70614E-02,-1.36249E-02,-1.06056E-02,-7.99360E-03,-5.77750E-03,  B16910
     * -3.94443E-03,-2.48028E-03,-1.36995E-03,-5.97540E-04,-1.46532E-04,  B16920
     *  0.         , 0.         /                                         B16930
C                                                                         B16940
      END                                                                 B16950
      SUBROUTINE RSYM (R,DV,VFT)                                          B16960
C                                                                         B16970
      IMPLICIT REAL*8           (V)                                     ! B16980
C                                                                         B16990
      DIMENSION R(*)                                                      B17000
C                                                                         B17010
      IP = (-VFT/DV)+1.-.000001                                           B17020
      IP = IP+1                                                           B17030
      P = (FLOAT(IP-1)+VFT/DV)*2.                                         B17040
      PST = P                                                             B17050
      IF (P.GT.1.) P = P-1.                                               B17060
C                                                                         B17070
C     VFT/DV- INT(VFT/DV)= 0. TO 0.5                                      B17080
C                                                                         B17090
      WN1 = -P*(P-1.)*(P-2.)/6.                                           B17100
      W0 = (P*P-1.)*(P-2.)/2.                                             B17110
      W1 = -P*(P+1.)*(P-2.)/2.                                            B17120
      W2 = P*(P*P-1.)/6.                                                  B17130
      K = IP                                                              B17140
      IPMAX = IP+IP-1                                                     B17150
      IF (PST.LE.1.) GO TO 20                                             B17160
      B1 = R(IP-2)                                                        B17170
      B2 = R(IP-1)                                                        B17180
      C1 = R(K)                                                           B17190
      DO 10 I = IP, IPMAX                                                 B17200
         K = K-1                                                          B17210
         C2 = C1                                                          B17220
         IF (K.LT.1) GO TO 40                                             B17230
         C1 = R(K)                                                        B17240
         R(K) = R(K)+WN1*R(I+1)+W0*R(I)+W1*B2+W2*B1                       B17250
         B1 = B2                                                          B17260
         B2 = R(I)                                                        B17270
         IF (K.LE.2) GO TO 10                                             B17280
         R(I) = R(I)+WN1*C2+W0*C1+W1*R(K-1)+W2*R(K-2)                     B17290
   10 CONTINUE                                                            B17300
      GO TO 40                                                            B17310
C                                                                         B17320
C    VFT/DV- INT(VFT/DV) = 0.5 TO 1.0                                     B17330
C                                                                         B17340
   20 C1 = R(IP)                                                          B17350
      C2 = R(IP+1)                                                        B17360
      B2 = R(IP-1)                                                        B17370
      DO 30 I = IP, IPMAX                                                 B17380
         K = K-1                                                          B17390
         B1 = B2                                                          B17400
         B2 = R(I)                                                        B17410
         IF (K.LE.1) GO TO 40                                             B17420
         R(I) = R(I)+WN1*C2+W0*C1+W1*R(K)+W2*R(K-1)                       B17430
         C2 = C1                                                          B17440
         C1 = R(K)                                                        B17450
         R(K) = R(K)+WN1*R(I+2)+W0*R(I+1)+W1*B2+W2*B1                     B17460
   30 CONTINUE                                                            B17470
C                                                                         B17480
   40 RETURN                                                              B17490
C                                                                         B17500
      END                                                                 B17510
      SUBROUTINE XINT (V1A,V2A,DVA,A,AFACT,VFT,DVR3,R3,N1R3,N2R3)         B17520
C                                                                         B17530
      IMPLICIT REAL*8           (V)                                     ! B17540
C                                                                         B17550
C     THIS SUBROUTINE INTERPOLATES THE A ARRAY STORED                     B17560
C     FROM V1A TO V2A IN INCREMENTS OF DVA USING A MULTIPLICATIVE         B17570
C     FACTOR AFACT, INTO THE R3 ARRAY FROM LOCATION N1R3 TO N2R3 IN       B17580
C     INCREMENTS OF DVR3                                                  B17590
C                                                                         B17600
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           B17610
      DIMENSION A(*),R3(*)                                                B17620
C                                                                         B17630
      RECDVA = 1./DVA                                                     B17640
      ILO = (V1A+DVA-VFT)/DVR3+1.+ONEMI                                   B17650
      ILO = MAX(ILO,N1R3)                                                 B17660
      IHI = (V2A-DVA-VFT)/DVR3+ONEMI                                      B17670
      IHI = MIN(IHI,N2R3)                                                 B17680
C                                                                         B17690
      DO 10 I = ILO, IHI                                                  B17700
         VI = VFT+DVR3*FLOAT(I-1)                                         B17710
         J = (VI-V1A)*RECDVA+ONEPL                                        B17720
         VJ = V1A+DVA*FLOAT(J-1)                                          B17730
         P = RECDVA*(VI-VJ)                                               B17740
         C = (3.-2.*P)*P*P                                                B17750
         B = 0.5*P*(1.-P)                                                 B17760
         B1 = B*(1.-P)                                                    B17770
         B2 = B*P                                                         B17780
         CONTI = -A(J-1)*B1+A(J)*(1.-C+B2)+A(J+1)*(C+B1)-A(J+2)*B2        B17790
         R3(I) = R3(I)+CONTI*AFACT                                        B17800
   10 CONTINUE                                                            B17810
C                                                                         B17820
      RETURN                                                              B17830
C                                                                         B17840
      END                                                                 B17850
      FUNCTION RADFN (VI,XKT)                                             B17860
C                                                                         B17870
      IMPLICIT REAL*8           (V)                                     ! B17880
C                                                                         B17890
C     FUNCTION RADFN CALCULATES THE RADIATION TERM FOR THE LINE SHAPE     B17900
C                                                                         B17910
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B17920
C                                                                         B17930
C               LAST MODIFICATION:    12 AUGUST 1991                      B17940
C                                                                         B17950
C                  IMPLEMENTATION:    R.D. WORSHAM                        B17960
C                                                                         B17970
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         B17980
C                                     R.D. WORSHAM                        B17990
C                                     J.L. MONCET                         B18000
C                                                                         B18010
C                                                                         B18020
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B18030
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          B18040
C                                                                         B18050
C----------------------------------------------------------------------   B18060
C                                                                         B18070
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     B18080
C                                     OFFICE OF ENERGY RESEARCH           B18090
C                                     DEPARTMENT OF ENERGY                B18100
C                                                                         B18110
C                                                                         B18120
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             B18130
C                                                                         B18140
C                                             FASCOD3                     B18150
C                                                                         B18160
C                                                                         B18170
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B18180
C                                                                         B18190
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           B18200
C                                                                         B18210
C      IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED                         B18220
C                                                                         B18230
      XVI = VI                                                            B18240
C                                                                         B18250
      IF (XKT.GT.0.0) THEN                                                B18260
C                                                                         B18270
         XVIOKT = XVI/XKT                                                 B18280
C                                                                         B18290
         IF (XVIOKT.LE.0.01) THEN                                         B18300
            RADFN = 0.5*XVIOKT*XVI                                        B18310
C                                                                         B18320
         ELSEIF (XVIOKT.LE.10.0) THEN                                     B18330
            EXPVKT = EXP(-XVIOKT)                                         B18340
            RADFN = XVI*(1.-EXPVKT)/(1.+EXPVKT)                           B18350
C                                                                         B18360
         ELSE                                                             B18370
            RADFN = XVI                                                   B18380
         ENDIF                                                            B18390
C                                                                         B18400
      ELSE                                                                B18410
         RADFN = XVI                                                      B18420
      ENDIF                                                               B18430
C                                                                         B18440
      RETURN                                                              B18450
C                                                                         B18460
      END                                                                 B18470
      FUNCTION RADFNI (VI,DVI,XKT,VINEW,RDEL,RDLAST)                      B18480
C                                                                         B18490
      IMPLICIT REAL*8           (V)                                     ! B18500
C                                                                         B18510
C     FUNCTION RADFNI CALCULATES THE RADIATION TERM FOR THE LINE SHAPE    B18520
C     OVER INTERVAL VI TO VINEW                                           B18530
C                                                                         B18540
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B18550
C                                                                         B18560
C               LAST MODIFICATION:    23 AUGUST 1991                      B18570
C                                                                         B18580
C                  IMPLEMENTATION:    R.D. WORSHAM                        B18590
C                                                                         B18600
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         B18610
C                                     R.D. WORSHAM                        B18620
C                                     J.L. MONCET                         B18630
C                                                                         B18640
C                                                                         B18650
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B18660
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          B18670
C                                                                         B18680
C----------------------------------------------------------------------   B18690
C                                                                         B18700
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     B18710
C                                     OFFICE OF ENERGY RESEARCH           B18720
C                                     DEPARTMENT OF ENERGY                B18730
C                                                                         B18740
C                                                                         B18750
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             B18760
C                                                                         B18770
C                                             FASCOD3                     B18780
C                                                                         B18790
C                                                                         B18800
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B18810
C                                                                         B18820
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           B18830
      DATA FACT1 / 3.0E-03 /                                              B18840
C                                                                         B18850
C     RADFNI IS COMPUTED AT VI AND AND CALCULATES THE                     B18860
C     WAVENUMBER VALUE (VINEW) FOR NEXT RADFNI CALC.                      B18870
C                                                                         B18880
C     IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED                          B18890
C                                                                         B18900
      XVI = VI                                                            B18910
C                                                                         B18920
C     IF FIRST CALL, INITIALIZE RDLAST                                    B18930
C                                                                         B18940
      IF (RDLAST.LT.0.) THEN                                              B18950
         IF (XKT.GT.0.0) THEN                                             B18960
            XVIOKT = XVI/XKT                                              B18970
C                                                                         B18980
            IF (XVIOKT.LE.0.01) THEN                                      B18990
               RDLAST = 0.5*XVIOKT*XVI                                    B19000
C                                                                         B19010
            ELSEIF (XVIOKT.LE.10.0) THEN                                  B19020
               EXPVKT = EXP(-XVIOKT)                                      B19030
               RDLAST = XVI*(1.-EXPVKT)/(1.+EXPVKT)                       B19040
C                                                                         B19050
            ELSE                                                          B19060
               RDLAST = XVI                                               B19070
            ENDIF                                                         B19080
         ELSE                                                             B19090
            RDLAST = XVI                                                  B19100
         ENDIF                                                            B19110
      ENDIF                                                               B19120
C                                                                         B19130
C     SET RADFNI EQUAL TO RADIATION FUNCTION AT VI                        B19140
C                                                                         B19150
C     RDLAST IS RADFNI(VI) FOR EACH SUBSEQUENT CALL                       B19160
C                                                                         B19170
      RADFNI = RDLAST                                                     B19180
C                                                                         B19190
      INTVLS = 1                                                          B19200
      IF (XKT.GT.0.0) THEN                                                B19210
C                                                                         B19220
         XVIOKT = XVI/XKT                                                 B19230
C                                                                         B19240
         IF (XVIOKT.LE.0.01) THEN                                         B19250
            IF (VINEW.GE.0.0) THEN                                        B19260
               VINEW = VI+FACT1*0.5*XVI                                   B19270
               INTVLS = (VINEW-VI)/DVI                                    B19280
               INTVLS = MAX(INTVLS,1)                                     B19290
               VINEW = VI+DVI*FLOAT(INTVLS)                               B19300
            ELSE                                                          B19310
               VINEW = ABS(VINEW)                                         B19320
               INTVLS = (VINEW-VI)/DVI                                    B19330
               INTVLS = MAX(INTVLS,1)                                     B19340
            ENDIF                                                         B19350
            XVINEW = VINEW                                                B19360
C                                                                         B19370
            RDNEXT = 0.5*XVIOKT*XVINEW                                    B19380
C                                                                         B19390
         ELSEIF (XVIOKT.LE.10.0) THEN                                     B19400
            EXPVKT = EXP(-XVIOKT)                                         B19410
            XMINUS = 1.-EXPVKT                                            B19420
            XPLUS = 1.+EXPVKT                                             B19430
            IF (VINEW.GE.0.0) THEN                                        B19440
               CVIKT = XVIOKT*EXPVKT                                      B19450
               VINEW = VI+FACT1*XVI/(1.+(CVIKT/XMINUS+CVIKT/XPLUS))       B19460
               INTVLS = (VINEW-VI)/DVI                                    B19470
               INTVLS = MAX(INTVLS,1)                                     B19480
               VINEW = VI+DVI*FLOAT(INTVLS)                               B19490
            ELSE                                                          B19500
               VINEW = ABS(VINEW)                                         B19510
               INTVLS = (VINEW-VI)/DVI                                    B19520
               INTVLS = MAX(INTVLS,1)                                     B19530
            ENDIF                                                         B19540
            XVINEW = VINEW                                                B19550
C                                                                         B19560
            RDNEXT = XVINEW*XMINUS/XPLUS                                  B19570
C                                                                         B19580
         ELSE                                                             B19590
            IF (VINEW.GE.0.0) THEN                                        B19600
               VINEW = VI+(FACT1*XVI)                                     B19610
               INTVLS = (VINEW-VI)/DVI                                    B19620
               INTVLS = MAX(INTVLS,1)                                     B19630
               VINEW = VI+DVI*FLOAT(INTVLS)                               B19640
            ELSE                                                          B19650
               VINEW = ABS(VINEW)                                         B19660
               INTVLS = (VINEW-VI)/DVI                                    B19670
               INTVLS = MAX(INTVLS,1)                                     B19680
            ENDIF                                                         B19690
            XVINEW = VINEW                                                B19700
C                                                                         B19710
            RDNEXT = XVINEW                                               B19720
         ENDIF                                                            B19730
      ELSE                                                                B19740
         IF (VINEW.GE.0.0) THEN                                           B19750
            VINEW = VI+(FACT1*XVI)                                        B19760
            INTVLS = (VINEW-VI)/DVI                                       B19770
            INTVLS = MAX(INTVLS,1)                                        B19780
            VINEW = VI+DVI*FLOAT(INTVLS)                                  B19790
         ELSE                                                             B19800
            VINEW = ABS(VINEW)                                            B19810
            INTVLS = (VINEW-VI)/DVI                                       B19820
            INTVLS = MAX(INTVLS,1)                                        B19830
         ENDIF                                                            B19840
         XVINEW = VINEW                                                   B19850
C                                                                         B19860
         RDNEXT = XVI                                                     B19870
      ENDIF                                                               B19880
C                                                                         B19890
      RDEL = (RDNEXT-RADFNI)/FLOAT(INTVLS)                                B19900
C                                                                         B19910
      RDLAST = RDNEXT                                                     B19930
C                                                                         B19940
      RETURN                                                              B19950
C                                                                         B19960
      END                                                                 B19970
      SUBROUTINE MOLEC (IND,SCOR,RHOSLF,ALFD1)                            C00010
C                                                                         C00020
      IMPLICIT REAL*8           (V)                                     ! C00030
C                                                                         C00040
      PARAMETER (NTMOL=36,Nspeci=85)
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4
      COMMON /ISVECT/ ISOVEC(NTMOL),ISO82(NSPECI),ISONM(NTMOL),           C00060
     *                SMASSI(NSPECI)                                      C00070
      COMMON /QTOT/ QCOEF(NSPECI,3,5),Q296(NSPECI),AQ(NSPECI),            C00080
     *              BQ(NSPECI),GJ(NSPECI)                                 C00090
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   C00100
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       C00110
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       C00120
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      C00130
     *              EXTID(10)                                             C00140
C                                                                         C00150
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ
C                                                                         C00170
      COMMON /FILHDR/ XID(10),SECANT,P   ,TEMP,HMOLID(60),XALTZ(4),       C00180
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   C00190
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    C00200
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           C00210
      DIMENSION SCOR(*),RHOSLF(*),ALFD1(*)                                C00220
      COMMON /SMOLEC/ W(42,9),ND(42,9),FAD                                C00230
      COMMON /XMOLEC/ NV(42),IVIB(42,2,9),XR(42),ROTFAC(42),QV0(42)
      COMMON /MOLNAM/ MOLID(0:NTMOL)                                      C00260
      CHARACTER*6 MOLID                                                   C00270
C                                                                         C00280
C     IS EQUIV. TO THE FOLLOWING DIMENSION AND EQUIVALENT STATEMENTS      C00290
C                                                                         C00300
C*****wog, 2/1/2001. 
C     DIMENSION IV(2)
      DIMENSION IV(42*2*9)
      EQUIVALENCE (IV(1),IVIB(1,1,1))                                     C00320
C                                                                         C00330
      DATA MDIM / 42 /,NVDIM / 9 /
C                                                                         C00350
c
c
c    Program TIPS written by R.R. Gamache
c
c     This is an updated version of TIPS: TIPS_97
c     obtained from R.R. Gamache on 28 April 1998
c
c***************
c
c     Modifications have been made to the partition sums for
c     hno3, c2h6, sf6, o, and clono2 by tony clough on 30 april 98,
c     based on data provided by R. Gamache.
c    
c***************
c
c    This program enables the user to calculate the Total Internal
c    Partition Sum (TIPS) for a given molecule, isotopic variant, and
c    temperature.  Current limitations are the molecular species on the
c    HITRAN molecular database and the temperature range 70 - 3000 K.
c
c...This program calculates the TIPS by the formula
c...     Q(T) = A + B*T + C*T*T + D*T*T*T + E*T*T*T*T
c...Reference: R.R. Gamache, R.L. Hawkins, and L.S. Rothman,
c    J.Mol.Spectrosc. 142, 205-219 (1990)
c
c     Program modified by Tony Clough
c      
c     Converted to single precision
c
c     Substantive changes have retained original record as comment with c%%
c
C     THIS PROGRAM ENABLES THE USER TO CALCULATE THE TOTAL INTERNAL       C00450
C     PARTITION SUM (TIPS) FOR A GIVEN MOLECULE, ISOTOPIC SPECIES,        C00460
C     AND TEMPERATURE.  CURRENT LIMITATIONS ARE THE MOLECULAR SPECIES     C00470
C     ON THE HITRAN MOLECULAR DATABASE AND THE TEMPERATURE RANGE          C00480
C     70 - 3000 K.                                                        C00490
C                                                                         C00500
C     MOLEC MAKES THE MOLECULAR IDENTIFICATIONS                           C00510
C                                                                         C00520
C     SCOR IS THE FACTOR BY WHICH THE LINE INTENSITY IS CHANGED DUE TO    C00530
C        TEMPERATURE DEPENDENCE OF THE VIB AND ROT PARTITION SUMS         C00540
C                                                                         C00550
C     RHOSLF IS A QUANTITY ('PARTIAL DENSITY') FOR CORRECTING THE         C00560
C        COLLISION WIDTH FOR SELF BROADENING                              C00570
C                                                                         C00580
C     ALFD1 CONTAINS THE DOPPLER WIDTHS AT 1 CM-1                         C00590
C                                                                         C00600
      IF (IND.EQ.1) THEN 
c
c...Set up molecule isotope vectors:
c
         CALL vecIso
c
c 5       WRITE (IPR,'(8x,A,/,8x,A////////////)')
c     +     'This program calculates the Total Internal Partition Sum'
c     +,'       for the molecular species on the HITRAN database.'
c%%   
C     
         DO 10 M = 1, NMOL                                                C00670
            READ (MOLID(M),900) HMOLID(M)                                 C00680
 10      CONTINUE                                                         C00690
C     
         FLN2 = ALOG(2.)                                                  C00710
         FAD = FLN2*2.*AVOG*BOLTZ/(CLIGHT*CLIGHT)                         C00720
         XKT0 = TEMP0/RADCN2                                              C00730
C     
         DO 30 M = 1, MDIM                                                C00760
            DO 20 K = 1, NVDIM                                            C00770
               LOC = 2*MDIM*(K-1)+2*(M-1)                                 C00780
               W(M,K) = IV(LOC+1)                                         C00790
               ND(M,K) = IV(LOC+2)                                        C00800
 20         CONTINUE                                                      C00810
            NVM = NV(M)                                                   C00820
 30      CONTINUE
         RETURN                                                           C00860
      ELSE                                                                C00870
C     

         RHORAT = (P/P0)*(TEMP0/TEMP)                                     C00890
         XKT = TEMP/RADCN2                                                C00900
C     
         DO 50 M = 1, NMOL
C     
            DO 40 ISO = 1, ISONM(M)
               ILOC = ISOVEC(M)+ISO                                       C00980
c     
               CALL QOFT (M,ISO,296.,QT_296) 
               CALL QOFT (M,ISO,TEMP,QT) 

               SCOR(iloc) = QT_296/QT


c              Stop program if covering a molecule and isotope
c              not valid for T > 500K.

               if ((iloc.ge.29).and.(wk(m).gt.0.).and.
     *              (TEMP.GT.500.)) then
                  write(ipr,*) 'TIPS calculation of Isotope ',
     *                 iso82(iloc),' for molecule ',m,
     *                 ' not valid for T = ',TEMP
                  write(*,*) 'TIPS calculation of Isotope ',
     *                 iso82(iloc),' for molecule ',m,
     *                 ' not valid for T = ',TEMP
                  stop 'SUBROUTINE MOLEC'
               endif

c              Stop program if covering a molecule of nonzero amount
c              and which has an isotope with no TIPs coefficients.

               if ((scor(iloc).lt.0.).and.(wk(m).gt.0.)) then
                  write(ipr,*) 'Isotope ',iso82(iloc),' for molecule ',
     *                 m,' not included in TIPs'
                  write(*,*) 'Isotope ',iso82(iloc),' for molecule ',
     *                 m,' not included in TIPs'
                  stop 'SUBROUTINE MOLEC'
               endif

               RHOSLF(ILOC) = RHORAT*WK(M)/WTOT                    
               ALFD1(ILOC) = SQRT(FAD*TEMP/SMASSI(ILOC))           
 40         CONTINUE                                               
 50      CONTINUE                                                
C     RETURN                                                              C01250
C     
      ENDIF                                                               C01270
C     
 900  FORMAT (A6)                                                         C01290
C                                                                         C01300
      END                                                                 C01310
      BLOCK DATA BMOLEC                                                   C01320
C                                                                         C01330
      COMMON /XMOLEC/                                                     C01340
     2  NV1(7),NV2(7),NV3(7),NV4(7),NV5(7),NV6(7),
     3  IV11(14),IV12(14),IV13(14),IV14(14),IV15(14),IV16(14),  
     4  IV21(14),IV22(14),IV23(14),IV24(14),IV25(14),IV26(14),  
     5  IV31(14),IV32(14),IV33(14),IV34(14),IV35(14),IV36(14),  
     6  IV41(14),IV42(14),IV43(14),IV44(14),IV45(14),IV46(14),  
     7  IV51(14),IV52(14),IV53(14),IV54(14),IV55(14),IV56(14),  
     8  IV61(14),IV62(14),IV63(14),IV64(14),IV65(14),IV66(14),  
     9  IV71(14),IV72(14),IV73(14),IV74(14),IV75(14),IV76(14),  
     *  IV81(14),IV82(14),IV83(14),IV84(14),IV85(14),IV86(14),  
     1  IV91(14),IV92(14),IV93(14),IV94(14),IV95(14),IV96(14),  
     2  XR1(7),XR2(7),XR3(7),XR4(7),XR5(7),XR6(7),                     
     3  ROTFC1(7),ROTFC2(7),ROTFC3(7),ROTFC4(7),ROTFC5(7),ROTFC6(7),      
     4  QV0(42)
C                                                                         C01480
      DATA NV1,IV11,IV21,IV31,IV41,IV51,IV61,IV71,IV81,IV91,XR1,ROTFC1/   C01490
C                                                                         C01500
C          H2O      CO2       O3      N2O       CO      CH4       O2      C01510
     C       3 ,      3 ,      3 ,      3 ,      1 ,      4 ,      1 ,    C01520
     1  3657,1 , 1388,1 , 1103,1 , 1285,1 , 2143,1 , 2917,1 , 1556,1 ,    C01530
     2  1595,1 ,  667,2 ,  701,1 ,  589,2 ,    0,0 , 1533,2 ,    0,0 ,    C01540
     3  3756,1 , 2349,1 , 1042,1 , 2224,1 ,    0,0 , 3019,3 ,    0,0 ,    C01550
     4     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 , 1311,3 ,    0,0 ,    C01560
     5     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C01570
     6     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C01580
     7     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C01590
     8     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C01600
     9     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C01610
     X     0.5 ,    0.25,    0.5 ,    0.5 ,    0.5 ,    0.5 ,    0.5 ,    C01620
     Y     1.5 ,    1.0 ,    1.5 ,    1.0 ,    1.0 ,    1.5 ,    1.0 /    C01630
C                                                                         C01640
      DATA NV2,IV12,IV22,IV32,IV42,IV52,IV62,IV72,IV82,IV92,XR2,ROTFC2/   C01650
C                                                                         C01660
C           NO      SO2      NO2      NH3     HNO3       OH       HF      C01670
     C       1 ,      3 ,      3 ,      4 ,      9 ,      1 ,      1 ,    C01680
     1  1876,1 , 1152,1 , 1318,1 , 3337,1 , 3550,1 , 3569,1 , 3961,1 ,    C01690
     2     0,0 ,  518,1 ,  750,1 ,  950,1 , 1710,1 ,    0,0 ,    0,0 ,    C01700
     3     0,0 , 1362,1 , 1617,1 , 3444,2 , 1331,1 ,    0,0 ,    0,0 ,    C01710
     4     0,0 ,    0,0 ,    0,0 , 1627,2 , 1325,1 ,    0,0 ,    0,0 ,    C01720
     5     0,0 ,    0,0 ,    0,0 ,    0,0 ,  879,1 ,    0,0 ,    0,0 ,    C01730
     6     0,0 ,    0,0 ,    0,0 ,    0,0 ,  647,1 ,    0,0 ,    0,0 ,    C01740
     7     0,0 ,    0,0 ,    0,0 ,    0,0 ,  579,1 ,    0,0 ,    0,0 ,    C01750
     8     0,0 ,    0,0 ,    0,0 ,    0,0 ,  762,1 ,    0,0 ,    0,0 ,    C01760
     9     0,0 ,    0,0 ,    0,0 ,    0,0 ,  456,1 ,    0,0 ,    0,0 ,    C01770
     X     0.5 ,    0.5 ,    0.5 ,    0.5 ,    0.5 ,    0.5 ,    0.5 ,    C01780
     Y     1.0 ,    1.5 ,    1.5 ,    1.5 ,    1.5 ,    1.0 ,    1.0 /    C01790
C                                                                         C01800
      DATA NV3,IV13,IV23,IV33,IV43,IV53,IV63,IV73,IV83,IV93,XR3,ROTFC3/   C01810
C                                                                         C01820
C          HCL      HBR       HI      CLO      OCS     H2CO     HOCL      C01830
     C       1 ,      1 ,      1 ,      1 ,      3 ,      6 ,      3 ,    C01840
     1  2885,1 , 2558,1 , 2229,1 ,  842,1 ,  859,1 , 2782,1 , 3609,1 ,    C01850
     2     0,0 ,    0,0 ,    0,0 ,    0,0 ,  520,2 , 1746,1 , 1238,1 ,    C01860
     3     0,0 ,    0,0 ,    0,0 ,    0,0 , 2062,1 , 1500,1 ,  740,1 ,    C01870
     4     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 , 1167,1 ,    0,0 ,    C01880
     5     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 , 2843,1 ,    0,0 ,    C01890
     6     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 , 1249,1 ,    0,0 ,    C01900
     7     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C01910
     8     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C01920
     9     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C01930
     X     0.5 ,    0.5 ,    0.5 ,    0.25,    0.5 ,    0.5 ,    0.5 ,    C01940
     Y     1.0 ,    1.0 ,    1.0 ,    1.0 ,    1.0 ,    1.5 ,    1.5 /    C01950
C                                                                         C01960
      DATA NV4,IV14,IV24,IV34,IV44,IV54,IV64,IV74,IV84,IV94,XR4,ROTFC4/   C01970
C                                                                         C01980
C           N2      HCN    CH3CL     H2O2     C2H2     C2H6      PH3      C01990
     C       1 ,      3 ,      6 ,      6 ,      5 ,      9 ,      4 ,    C02000
     1  2330,1 , 2089,1 , 2968,1 , 3607,1 , 3374,1 , 2899,1 , 2327,1 ,    C02010
     2     0,0 ,  713,2 , 1355,1 , 1394,1 , 1974,1 , 1375,1 ,  992,1 ,    C02020
     3     0,0 , 3311,1 ,  732,1 ,  864,1 , 3295,1 ,  993,1 , 1118,2 ,    C02030
     4     0,0 ,    0,0 , 3039,2 ,  317,1 ,  612,2 ,  275,1 , 2421,2 ,    C02040
     5     0,0 ,    0,0 , 1455,2 , 3608,1 ,  730,2 , 2954,1 ,    0,0 ,    C02050
     6     0,0 ,    0,0 , 1015,2 , 1269,1 ,    0,0 , 1379,1 ,    0,0 ,    C02060
     7     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 , 2994,2 ,    0,0 ,    C02070
     8     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 , 1486,1 ,    0,0 ,    C02080
     9     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,  822,2 ,    0,0 ,    C02090
     X     0.5 ,    0.5 ,    0.5 ,    0.5 ,    0.5 ,    0.5 ,    0.5 ,    C02100
     Y     1.0 ,    1.0 ,    1.5 ,    1.5 ,    1.0 ,    1.5 ,    1.5 /    C02110
C                                                                         C02120
      DATA NV5,IV15,IV25,IV35,IV45,IV55,IV65,IV75,IV85,IV95,XR5,ROTFC5/   C02130
C                                                                         C02140
C         COF2      SF6      H2S    HCOOH      HO2        O   CLONO2      C02150
     C       0 ,      0 ,      0 ,      0 ,      0 ,      0 ,      0 ,    C02160
     1  0000,1 , 0000,0 , 0000,0 , 0000,0 , 0000,0 , 0000,0 , 0000,0 ,    C02170
     2     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02180
     3     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02190
     4     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02200
     5     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02210
     6     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02220
     7     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02230
     8     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02240
     9     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02250
     X     0.0 ,    0.0 ,    0.0 ,    0.0 ,    0.0 ,    0.0 ,    0.0 ,    C02260
     Y     0.0 ,    0.0 ,    0.0 ,    0.0 ,    0.0 ,    0.0 ,    0.0 /    C02270
C                                                                         C02280
      DATA NV6,IV16,IV26,IV36,IV46,IV56,IV66,IV76,IV86,IV96,XR6,ROTFC6/ 
C                                                                       
C          NO+      ???      ???      ???      ???      ???      ??? 
     C       0 ,      0 ,      0 ,      0 ,      0 ,      0 ,      0 ,    C02160
     1  0000,1 , 0000,0 , 0000,0 , 0000,0 , 0000,0 , 0000,0 , 0000,0 ,    C02170
     2     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02180
     3     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02190
     4     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02200
     5     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02210
     6     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02220
     7     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02230
     8     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02240
     9     0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    0,0 ,    C02250
     X     0.0 ,    0.0 ,    0.0 ,    0.0 ,    0.0 ,    0.0 ,    0.0 ,    C02260
     Y     0.0 ,    0.0 ,    0.0 ,    0.0 ,    0.0 ,    0.0 ,    0.0 /    C02270
C                                                                         C02280
      END                                                                 C02290
      FUNCTION QV (M,XKT,W,ND,NV,MDIM,NVDIM)                              C02300
C                                                                         C02310
C     FUNCTION QV CALCULATES THE VIBRATIONAL PARTITION SUM                C02320
C                                                                         C02330
      DIMENSION W(MDIM,NVDIM),ND(MDIM,NVDIM)                              C02340
      QV = 1.                                                             C02350
      DO 10 I = 1, NV                                                     C02360
         SV = 1.-EXP(-W(M,I)/XKT)                                         C02370
         IF (ND(M,I).GT.1) SV = SV**ND(M,I)                               C02380
         QV = QV/SV                                                       C02390
   10 CONTINUE                                                            C02400
C                                                                         C02410
      RETURN                                                              C02420
C                                                                         C02430
      END                                                                 C02440
c***********************************************
      SUBROUTINE vecIso
c***********************************************
c
      PARAMETER (NMOL=36,Nspeci=85)
      COMMON /ISVECT/ ISOVEC(NMOL),ISO82(Nspeci),ISONM(NMOL),
     *     sdum(Nspeci)
c
c...Isotope vector information
c     Set up ISOVEC:
         ISOVEC(1) = 0
         DO 20 I = 2,NMOL
          ISOVEC(I) = 0
          DO 10 J = 1,I-1
           ISOVEC(I) = ISOVEC(I)+ISONM(J)
   10     CONTINUE
   20    CONTINUE
c
      RETURN
      END
c  ****************************************
      BLOCK DATA Isotop
c  ****************************************
C$$   IMPLICIT DOUBLE PRECISION (a-h,o-z)
c
      PARAMETER (NMOL=36,Nspeci=85)
      COMMON /ISVECT/ ISOVEC(NMOL),ISO82(Nspeci),ISONM(NMOL),
     *     smassi(Nspeci)
c
c    The number of isotopes for a particular molecule:
      DATA (ISONM(I),I=1,NMOL)/
c     H2O, CO2, O3, N2O, CO, CH4, O2,
     +  4,   8,  5,   5,  6,   3,  3,
c      NO, SO2, NO2, NH3, HNO3, OH, HF, HCl, HBr, HI,
     +  3,   2,   1,   2,    1,  3,  1,   2,   2,  1,
c     ClO, OCS, H2CO, HOCl, N2, HCN, CH3Cl, H2O2, C2H2, C2H6, PH3
     +  2,   4,    3,    2,  1,   3,     2,    1,    2,    1,   1,
c     COF2, SF6, H2S, HCOOH, HO2, O, ClONO2,  NO+
     +  1,   1,   3,     1,   1,  1,     2,    1/
c
      DATA ISO82/
c       H2O
     +  161,181,171,162,                                           
c       CO2
     +  626,636,628,627,638,637,828,728,
c       O3
     +  666,668,686,667,676,
c       N2O
     +  446,456,546,448,447,
c       CO,                 CH4
     +  26,36,28,27,38,37,  211,311,212,
c       O2,        NO,        SO2
     +  66,68,67,  46,56,48  ,626,646,
c      NO2,   NH3,        HNO3
     + 646,   4111,5111,  146,
c       OH,        HF,  HCl,    HBr,    HI
     +  61,81,62,  19,  15,17,  19,11,  17,
c       ClO,    OCS,              H2CO
     +  56,76,  622,624,632,822,  126,136,128,
c       HOCl,     N2,  HCN
     +  165,167,  44,  124,134,125
c      CH3Cl,    H2O2,  C2H2,       C2H6,  PH3
     +,215,217,  1661,  1221,1231,  1221,  1111,
c     COF2, SF6, H2S,            HCOOH,  HO2, O,   ClONO2      NO+
     + 269,  29,  121,141,131,   126,    166, 6,   5646,7646,  46/
c
C                                                                         C03620
C     MOLECULAR MASSES FOR EACH ISOTOPE                                   C03630
C                                                                         C03640
      DATA SMASSI/                                                        C03650
C     H2O:   161,   181,   171,   162                                     C03660
     *       18.01, 20.01, 19.01, 19.01,                                  C03670
C     CO2:   626,   636,   628,   627,   638,   637,   828,   728         C03680
     *       43.98, 44.98, 45.98, 44.98, 46.98, 45.98, 47.98, 46.98,      C03690
C     O3:    666,   668,   686    667    676                              ******
     *       47.97, 49.97, 49.97, 48.99, 48.99,                           ******
C     N2O:   446,   456,   546,   448,   447                              C03720
     *       43.99, 44.99, 44.99, 45.99, 44.99,                           C03730
C     CO:    26,    36,    28,    27,    38     37                        ******
     *       27.99, 28.99, 29.99, 28.99, 30.99, 30.00,                    ******
C     CH4:   211,   311,   212;   O2:    66,    68,    67                 C03760
     *       16.04, 17.04, 17.04,        31.98, 33.98, 32.98,             C03770
C     NO:    46,    56,    48;    SO2:   626,   646                       C03780
     *       29.99, 30.99, 31.99,        63.95, 65.95,                    C03790
C     NO2:   646;   NH3:   4111,  5111;  HNO3:  146                       C03800
     *       45.98,        17.03, 18.03,        62.98,                    C03810
C     OH:    61,    81,    62;    HF:    19                               C03820
     *       17.00, 19.00, 18.00,        20.01,                           C03830
C     HCL:   15,    17;    HBR:   19,    11;    HI:    17                 C03840
     *       35.98, 37.98,        79.92, 81.92,        127.91,            C03850
C     CLO:   56,    76;    OCS:   622,   624,   632,   822                C03860
     *       50.96, 52.96,        59.96, 61.96, 60.96, 61.96,             C03870
C     H2CO:  126,   136,   128;   HOCL:  165,   167                       C03880
     *       30.01, 31.01, 32.01,        51.97, 53.97,                    C03890
C     N2:    44;    HCN:   124,   134,   125                              C03900
     *       28.00,        27.01, 28.01, 28.01,                           C03910
C     CH3CL: 215,   217;   H2O2:  1661;  C2H2:  1221,  1231               C03920
     *       50.00, 52.00,        34.00,        26.02, 27.02,             C03930
C     C2H6:  1221;  PH3:   1111;  COF2:  269;   SF6:   29                 C03940
     *       30.06,        34.00,        65.99,        145.97,            C03950
C     H2S:   121    141    131;   HCOOH: 126;   HO2:   166                ******
     *       33.99, 35.98, 34.99,        46.00,        33.00,             ******
C     O:     6      ClONO2:5646   7646;  NO+:   46                        ******
     *       15.99,        96.96, 98.95,        30.00 /                   ******
C
      END

c
C*******************************************
      SUBROUTINE QofT(Mol, Iso, Tout, QT)
C*******************************************
c...date last changed:  September 19, 1997
c
c...input - Mol, Iso, and a temperature Tout
c...output - Total internal partition sum, QT, at T=Tout
c
c$$$      IMPLICIT DOUBLE PRECISION (a-h,o-z)
c
c++
      PARAMETER (NMOL=36,Nspeci=85)
c++
      COMMON /ISVECT/ ISOVEC(NMOL),ISO82(Nspeci),ISONM(NMOL),
     *     sdum(Nspeci)
c++:  bd-QT
      common/Qtot/ Qcoef(Nspeci,3,5), Q296(Nspeci), aQ(Nspeci), 
     + bQ(Nspeci), gj(Nspeci)
c
      ivec = isovec(Mol) + iso
c
c...value depends on temperature range
      IF (Tout.lt.70. .OR. Tout.gt.3005.) THEN
        QT = -1.
        WRITE (*,*) '  OUT OF TEMPERATURE RANGE'
        GOTO 99
      ENDIF
      IF (Tout .ge.  70. .and. Tout .le.  500.) irange = 1
      IF (Tout .gt. 500. .and. Tout .le. 1500.) irange = 2
      IF (Tout .gt. 1500.) irange = 3
c
        QT = Qcoef(ivec,irange,1)
     +       + Qcoef(ivec,irange,2)*Tout
     +       + Qcoef(ivec,irange,3)*Tout*Tout
     +       + Qcoef(ivec,irange,4)*Tout*Tout*Tout
     +       + Qcoef(ivec,irange,5)*Tout*Tout*Tout*Tout
c
   99 RETURN
      END

c**************************************
      BLOCK DATA QTdata
c**************************************
c$$$      IMPLICIT DOUBLE PRECISION (a-h,o-z)
c
c
c++
      PARAMETER (NMOL=36,Nspeci=85)
c++:  bd-QT
      common/Qtot/ Qcoef(Nspeci,3,5), Q296(Nspeci), aQ(Nspeci), 
     + bQ(Nspeci), gj(Nspeci)
c
c...State independent degeneracy factors: gj
c...(includes general nuclear factor (P(2I+1)), (2S+1), and (2-dl0)
c... when applicable)
c              H2O              CO2                     O3
      DATA gj/1.,1.,6.,6.,  1.,2.,1.,6.,2.,12.,1.,6.,  1.,1.,1.,6.,6.,
c              N2O                CO             CH4        O2
     +  9.,6.,6.,9.,54.,  1.,2.,1.,6.,2.,12.,  1.,2.,3.,  1.,1.,6.,
c          NO        SO2    NO2    NH3   HNO3    OH        HF    HCl
     +  1.,1.,1.,  1.,1.,  3.,  3.,2.,  6.,  1.,1.,1.,  4.,  8.,8.,
c        HBr    HI     ClO      OCS          H2CO      HOCl         N2
     +  8.,8.,  12.,  1.,1.,  1.,1.,2.,1.,  1.,2.,1.,  8.,8.,  1,
c        HCN        CH3Cl  H2O2  C2H2   C2H6  PH3
     +  6.,12.,4.,  4.,4.,  1.,  1.,8.,  64.,  2.,
c      COF2,  SF6,  H2S,  HCOOH
     +  1.,  1.,  1.,1.,4.,  4.,
c        HO2, O,   ClONO2    NO+
     +   2.,  1.,   12.,12.,   3./
c
c...Total internal partition sums for T<=70 to <=500 K range:
c...   H2O  --   161
      DATA (Qcoef( 1,1,j),j=1,5)/-.44405E+01, .27678E+00,
     +                .12536E-02,-.48938E-06, 0.   /
c...   H2O  --   181
      DATA (Qcoef( 2,1,j),j=1,5)/-.43624E+01, .27647E+00,
     +                .12802E-02,-.52046E-06, 0.   /
c...   H2O  --   171
      DATA (Qcoef( 3,1,j),j=1,5)/-.25767E+02, .16458E+01,
     +                .76905E-02,-.31668E-05, 0.   /
c...   H2O  --   162
      DATA (Qcoef( 4,1,j),j=1,5)/-.23916E+02, .13793E+01,
     +                .61246E-02,-.21530E-05, 0.   /
c...   CO2  --   626
      DATA (Qcoef( 5,1,j),j=1,5)/-.13617E+01, .94899E+00,
     +               -.69259E-03, .25974E-05, 0.   /
c...   CO2  --   636
      DATA (Qcoef( 6,1,j),j=1,5)/-.20631E+01, .18873E+01,
     +               -.13669E-02, .54032E-05, 0.   /
c...   CO2  --   628
      DATA (Qcoef( 7,1,j),j=1,5)/-.29175E+01, .20114E+01,
     +               -.14786E-02, .55941E-05, 0.   /
c...   CO2  --   627
      DATA (Qcoef( 8,1,j),j=1,5)/-.16558E+02, .11733E+02,
     +               -.85844E-02, .32379E-04, 0.   /
c...   CO2  --   638
      DATA (Qcoef( 9,1,j),j=1,5)/-.44685E+01, .40330E+01,
     +               -.29590E-02, .11770E-04, 0.   /
c...   CO2  --   637
      DATA (Qcoef(10,1,j),j=1,5)/-.26263E+02, .23350E+02,
     +               -.17032E-01, .67532E-04, 0.   /
c...   CO2  --   828
      DATA (Qcoef(11,1,j),j=1,5)/-.14811E+01, .10667E+01,
     +               -.78758E-03, .30133E-05, 0.   /
c...   CO2  --   728
      DATA (Qcoef(12,1,j),j=1,5)/-.17600E+02, .12445E+02,
     +               -.91837E-02, .34915E-04, 0.   /
c...    O3  --   666
      DATA (Qcoef(13,1,j),j=1,5)/-.16443E+03, .69047E+01,
     +                .10396E-01, .26669E-04, 0.   /
c...    O3  --   668
      DATA (Qcoef(14,1,j),j=1,5)/-.35222E+03, .14796E+02,
     +                .21475E-01, .59891E-04, 0.   /
c...    O3  --   686
      DATA (Qcoef(15,1,j),j=1,5)/-.17466E+03, .72912E+01,
     +                .10093E-01, .29991E-04, 0.   /
c...    O3  --   667
      DATA (Qcoef(16,1,j),j=1,5)/-.20540E+04, .85998E+02,
     +                .12667E+00, .33026E-03, 0.   /
c...    O3  --   676
      DATA (Qcoef(17,1,j),j=1,5)/-.10148E+04, .42494E+02,
     +                .62586E-01, .16319E-03, 0.   /
c...   N2O  --   446
      DATA (Qcoef(18,1,j),j=1,5)/ .24892E+02, .14979E+02,
     +               -.76213E-02, .46310E-04, 0.   /
c...   N2O  --   456
      DATA (Qcoef(19,1,j),j=1,5)/ .36318E+02, .95497E+01,
     +               -.23943E-02, .26842E-04, 0.   /
c...   N2O  --   546
      DATA (Qcoef(20,1,j),j=1,5)/ .24241E+02, .10179E+02,
     +               -.43002E-02, .30425E-04, 0.   /
c...   N2O  --   448
      DATA (Qcoef(21,1,j),j=1,5)/ .67708E+02, .14878E+02,
     +               -.10730E-02, .34254E-04, 0.   /
c...   N2O  --   447
      DATA (Qcoef(22,1,j),j=1,5)/ .50069E+03, .84526E+02,
     +                .83494E-02, .17154E-03, 0.   /
c...    CO  --    26
      DATA (Qcoef(23,1,j),j=1,5)/ .27758E+00, .36290E+00,
     +               -.74669E-05, .14896E-07, 0.   /
c...    CO  --    36
      DATA (Qcoef(24,1,j),j=1,5)/ .53142E+00, .75953E+00,
     +               -.17810E-04, .35160E-07, 0.   /
c...    CO  --    28
      DATA (Qcoef(25,1,j),j=1,5)/ .26593E+00, .38126E+00,
     +               -.92083E-05, .18086E-07, 0.   /
c...    CO  --    27                              
      DATA (Qcoef(26,1,j),j=1,5)/ .16376E+01, .22343E+01,
     +               -.49025E-04, .97389E-07, 0.   /
c...    CO  --    38                              
      DATA (Qcoef(27,1,j),j=1,5)/ .51216E+00, .79978E+00,
     +               -.21784E-04, .42749E-07, 0.   /
c...    CO  --    37                              
      DATA (Qcoef(28,1,j),j=1,5)/ .32731E+01, .46577E+01,
     +               -.69833E-04, .18853E-06, 0.   /
c...   CH4  --   211                              
      DATA (Qcoef(29,1,j),j=1,5)/-.26479E+02, .11557E+01,
     +                .26831E-02, .15117E-05, 0.   /
c...   CH4  --   311                              
      DATA (Qcoef(30,1,j),j=1,5)/-.52956E+02, .23113E+01,
     +                .53659E-02, .30232E-05, 0.   /
c...   CH4  --   212                              
      DATA (Qcoef(31,1,j),j=1,5)/-.21577E+03, .93318E+01,
     +                .21779E-01, .12183E-04, 0.   /
c...    O2  --    66                              
      DATA (Qcoef(32,1,j),j=1,5)/ .35923E+00, .73534E+00,
     +               -.64870E-04, .13073E-06, 0.   /
c...    O2  --    68                              
      DATA (Qcoef(33,1,j),j=1,5)/-.40039E+01, .15595E+01,
     +               -.15357E-03, .30969E-06, 0.   /
c...    O2  --    67                              
      DATA (Qcoef(34,1,j),j=1,5)/-.23325E+02, .90981E+01,
     +               -.84435E-03, .17062E-05, 0.   /
c...    NO  --    46                              
      DATA (Qcoef(35,1,j),j=1,5)/-.25296E+02, .26349E+01,
     +                .58517E-02,-.52020E-05, 0.   /
c...    NO  --    56                              
      DATA (Qcoef(36,1,j),j=1,5)/-.14990E+02, .18240E+01,
     +                .40261E-02,-.35648E-05, 0.   /
c...    NO  --    48                              
      DATA (Qcoef(37,1,j),j=1,5)/-.26853E+02, .27816E+01,
     +                .61493E-02,-.54410E-05, 0.   /
c...   SO2  --   626                              
      DATA (Qcoef(38,1,j),j=1,5)/-.24056E+03, .11101E+02,
     +                .22164E-01, .52334E-04, 0.   /
c...   SO2  --   646                              
      DATA (Qcoef(39,1,j),j=1,5)/-.24167E+03, .11151E+02,
     +                .22270E-01, .52550E-04, 0.   /
c...   NO2  --   646                              
      DATA (Qcoef(40,1,j),j=1,5)/-.53042E+03, .24216E+02,
     +                .66856E-01, .43823E-04, 0.   /
c...   NH3  --  4111                              
      DATA (Qcoef(41,1,j),j=1,5)/-.42037E+02, .25976E+01,
     +                .13073E-01,-.62230E-05, 0.   /
c...   NH3  --  5111                              
      DATA (Qcoef(42,1,j),j=1,5)/-.28609E+02, .17272E+01,
     +                .87529E-02,-.41714E-05, 0.   /

c**********
c...  HNO3  --   146                              
c      DATA (Qcoef(43,1,j),j=1,5)/-.10000E+01, .00000E+00,
c     +                .00000E+00, .00000E+00, 0.   /

      DATA (Qcoef(43,1,j),j=1,5)/ -6.672718E+4, 1.462506E+3,
     +               -5.981021E+0, 1.414328E-2, 0.  /

c 1.414328E-2*x^3 + -5.981021E+0*x^2 + 1.462506E+3*x + -6.672718E+4
 
c...    OH  --    61                              
      DATA (Qcoef(44,1,j),j=1,5)/ .87390E+01, .15977E+00,
     +                .38291E-03,-.35669E-06, 0.   /
c...    OH  --    81                              
      DATA (Qcoef(45,1,j),j=1,5)/ .86770E+01, .16175E+00,
     +                .38223E-03,-.35466E-06, 0.   /
c...    OH  --    62                              
      DATA (Qcoef(46,1,j),j=1,5)/ .10239E+02, .43783E+00,
     +                .10477E-02,-.94570E-06, 0.   /
c...    HF  --    19                              
      DATA (Qcoef(47,1,j),j=1,5)/ .15486E+01, .13350E+00,
     +                .59154E-05,-.46889E-08, 0.   /
c...   HCl  --    15                              
      DATA (Qcoef(48,1,j),j=1,5)/ .28627E+01, .53122E+00,
     +                .67464E-05,-.16730E-08, 0.   /
c...   HCl  --    17                              
      DATA (Qcoef(49,1,j),j=1,5)/ .28617E+01, .53203E+00,
     +                .66553E-05,-.15168E-08, 0.   /
c...   HBr  --    19                              
      DATA (Qcoef(50,1,j),j=1,5)/ .27963E+01, .66532E+00,
     +                .34255E-05, .52274E-08, 0.   /
c...   HBr  --    11                              
      DATA (Qcoef(51,1,j),j=1,5)/ .27953E+01, .66554E+00,
     +                .32931E-05, .54823E-08, 0.   /
c...    HI  --    17                              
      DATA (Qcoef(52,1,j),j=1,5)/ .40170E+01, .13003E+01,
     +               -.11409E-04, .40026E-07, 0.   /
c...   ClO  --    56                              
      DATA (Qcoef(53,1,j),j=1,5)/ .90968E+02, .70918E+01,
     +                .11639E-01, .30145E-05, 0.   /
c...   ClO  --    76                              
      DATA (Qcoef(54,1,j),j=1,5)/ .92598E+02, .72085E+01,
     +                .11848E-01, .31305E-05, 0.   /
c...   OCS  --   622                              
      DATA (Qcoef(55,1,j),j=1,5)/-.93697E+00, .36090E+01,
     +               -.34552E-02, .17462E-04, 0.   /
c...   OCS  --   624                              
      DATA (Qcoef(56,1,j),j=1,5)/-.11536E+01, .37028E+01,
     +               -.35582E-02, .17922E-04, 0.   /
c...   OCS  --   632                              
      DATA (Qcoef(57,1,j),j=1,5)/-.61015E+00, .72200E+01,
     +               -.70044E-02, .36708E-04, 0.   /
c...   OCS  --   822                              
      DATA (Qcoef(58,1,j),j=1,5)/-.21569E+00, .38332E+01,
     +               -.36783E-02, .19177E-04, 0.   /
c...  H2CO  --   126                              
      DATA (Qcoef(59,1,j),j=1,5)/-.11760E+03, .46885E+01,
     +                .15088E-01, .35367E-05, 0.   /
c...  H2CO  --   136                              
      DATA (Qcoef(60,1,j),j=1,5)/-.24126E+03, .96134E+01,
     +                .30938E-01, .72579E-05, 0.   /
c...  H2CO  --   128                              
      DATA (Qcoef(61,1,j),j=1,5)/-.11999E+03, .52912E+01,
     +                .14686E-01, .43505E-05, 0.   /
c...  HOCl  --   165                              
      DATA (Qcoef(62,1,j),j=1,5)/-.73640E+03, .34149E+02,
     +                .93554E-01, .67409E-04, 0.   /
c...  HOCl  --   167                              
      DATA (Qcoef(63,1,j),j=1,5)/-.74923E+03, .34747E+02,
     +                .95251E-01, .68523E-04, 0.   /
c...    N2  --    44                              
      DATA (Qcoef(64,1,j),j=1,5)/ .13684E+01, .15756E+01,
     +               -.18511E-04, .38960E-07, 0.   /
c...   HCN  --   124                              
      DATA (Qcoef(65,1,j),j=1,5)/-.13992E+01, .29619E+01,
     +               -.17464E-02, .65937E-05, 0.   /
c...   HCN  --   134                              
      DATA (Qcoef(66,1,j),j=1,5)/-.25869E+01, .60744E+01,
     +               -.35719E-02, .13654E-04, 0.   /
c...   HCN  --   125                              
      DATA (Qcoef(67,1,j),j=1,5)/-.11408E+01, .20353E+01,
     +               -.12159E-02, .46375E-05, 0.   /
c... CH3Cl  --   215                              
      DATA (Qcoef(68,1,j),j=1,5)/-.91416E+03, .34081E+02,
     +                .75461E-02, .17933E-03, 0.   /
c... CH3Cl  --   217                              
      DATA (Qcoef(69,1,j),j=1,5)/-.92868E+03, .34621E+02,
     +                .76674E-02, .18217E-03, 0.   /
c...  H2O2  --  1661                              
      DATA (Qcoef(70,1,j),j=1,5)/-.36499E+03, .13712E+02,
     +                .38658E-01, .23052E-04, 0.   /
c...  C2H2  --  1221                              
      DATA (Qcoef(71,1,j),j=1,5)/-.83088E+01, .14484E+01,
     +               -.25946E-02, .84612E-05, 0.   /
c...  C2H2  --  1231                              
      DATA (Qcoef(72,1,j),j=1,5)/-.66736E+02, .11592E+02,
     +               -.20779E-01, .67719E-04, 0.   /


c**********
c...  C2H6  --  1221                              
c      DATA (Qcoef(73,1,j),j=1,5)/-.10000E+01, .00000E+00,
c     +                .00000E+00, .00000E+00, 0.   /
      DATA (Qcoef(73,1,j),j=1,5)/-1.198174E+4, 2.799115E+2,
     +              -7.416173E-1, 1.846289E-3, 0.   /

c 1.846289E-3*x^3 + -7.416173E-1*x^2 + 2.799115E+2*x + -1.198174E+4


c...   PH3  --  1111                              
      DATA (Qcoef(74,1,j),j=1,5)/-.15068E+03, .64718E+01,
     +                .12588E-01, .14759E-04, 0.   /
c...  COF2  --   269                              
      DATA (Qcoef(75,1,j),j=1,5)/-.54180E+04, .18868E+03,
     +               -.33139E+00, .18650E-02, 0.   /

c**********
c...   SF6  --    29                              
c      DATA (Qcoef(76,1,j),j=1,5)/-.10000E+01, .00000E+00,
c     +                .00000E+00, .00000E+00, 0.   /
      DATA (Qcoef(76,1,j),j=1,5)/7.892502E+6, -1.839926E+5,
     +             1.433742E+3, -4.510500E+0,  5.163613E-3/

c 5.163613E-3*x^4 + -4.510500E+0*x^3 + 1.433742E+3*x^2 + -1.839926E+5*x + 7.892502E+6


c...   H2S  --   121                              
      DATA (Qcoef(77,1,j),j=1,5)/-.15521E+02, .83130E+00,
     +                .33656E-02,-.85691E-06, 0.   /
c...   H2S  --   141                              
      DATA (Qcoef(78,1,j),j=1,5)/-.15561E+02, .83337E+00,
     +                .33744E-02,-.85937E-06, 0.   /
c...   H2S  --   131                              
      DATA (Qcoef(79,1,j),j=1,5)/-.62170E+02, .33295E+01,
     +                .13480E-01,-.34323E-05, 0.   /
c... HCOOH  --   126                              
      DATA (Qcoef(80,1,j),j=1,5)/-.29550E+04, .10349E+03,
     +               -.13146E+00, .87787E-03, 0.   /
c...   HO2  --   166                              
      DATA (Qcoef(81,1,j),j=1,5)/-.15684E+03, .74450E+01,
     +                .26011E-01,-.92704E-06, 0.   /

c**********
c...     O  --     6                              
c      DATA (Qcoef(82,1,j),j=1,5)/-.10000E+01, .00000E+00,
c     +                .00000E+00, .00000E+00, 0.   /

      DATA (Qcoef(82,1,j),j=1,5)/  1.0       , .00000E+00,
     +                 .00000E+00, .00000E+00, 0.   /


c**********
c...ClONO2  --  5646                              
c      DATA (Qcoef(83,1,j),j=1,5)/-.10000E+01, .00000E+00,
c     +                .00000E+00, .00000E+00, 0.   /
      DATA (Qcoef(83,1,j),j=1,5)/2.516742E+6, -5.996989E+4,
     +             5.111204E+2, -1.665663E+0,  2.399248E-3/

c 2.399248E-3*x^4 + -1.665663E+0*x^3 + 5.111204E+2*x^2 + -5.996989E+4*x + 2.516742E+6


c**********
c...ClONO2  --  7646                              
c      DATA (Qcoef(84,1,j),j=1,5)/-.10000E+01, .00000E+00,
c     +                .00000E+00, .00000E+00, 0.   /
      DATA (Qcoef(84,1,j),j=1,5)/2.601684E+6, -6.187640E+4,
     +            5.260888E+2,  -1.711186E+0,  2.461028E-3/

c 2.461028E-3*x^4 + -1.711186E+0*x^3 + 5.260888E+2*x^2 + -6.187640E+4*x + 2.601684E+6

c...   NO+  --    46                              
      DATA (Qcoef(85,1,j),j=1,5)/ .91798E+00, .10416E+01,
     +               -.11614E-04, .24499E-07, 0.   /
c
c...Total internal partition sums for T>500 to <=1500 K range:
c...   H2O  --   161                              
      DATA (Qcoef( 1,2,j),j=1,5)/-.94327E+02, .81903E+00,
     +                .74005E-04, .42437E-06, 0.   /
c...   H2O  --   181                              
      DATA (Qcoef( 2,2,j),j=1,5)/-.95686E+02, .82839E+00,
     +                .68311E-04, .42985E-06, 0.   /
c...   H2O  --   171                              
      DATA (Qcoef( 3,2,j),j=1,5)/-.57133E+03, .49480E+01,
     +                .41517E-03, .25599E-05, 0.   /
c...   H2O  --   162                              
      DATA (Qcoef( 4,2,j),j=1,5)/-.53366E+03, .44246E+01,
     +               -.46935E-03, .29548E-05, 0.   /
c...   CO2  --   626                              
      DATA (Qcoef( 5,2,j),j=1,5)/-.50925E+03, .32766E+01,
     +               -.40601E-02, .40907E-05, 0.   /
c...   CO2  --   636                              
      DATA (Qcoef( 6,2,j),j=1,5)/-.11171E+04, .70346E+01,
     +               -.89063E-02, .88249E-05, 0.   /
c...   CO2  --   628                              
      DATA (Qcoef( 7,2,j),j=1,5)/-.11169E+04, .71299E+01,
     +               -.89194E-02, .89268E-05, 0.   /
c...   CO2  --   627                              
      DATA (Qcoef( 8,2,j),j=1,5)/-.66816E+04, .42402E+02,
     +               -.53269E-01, .52774E-04, 0.   /
c...   CO2  --   638                              
      DATA (Qcoef( 9,2,j),j=1,5)/-.25597E+04, .15855E+02,
     +               -.20440E-01, .19855E-04, 0.   /
c...   CO2  --   637                              
      DATA (Qcoef(10,2,j),j=1,5)/-.14671E+05, .91204E+02,
     +               -.11703E+00, .11406E-03, 0.   /
c...   CO2  --   828                              
      DATA (Qcoef(11,2,j),j=1,5)/-.63775E+03, .40047E+01,
     +               -.50950E-02, .50023E-05, 0.   /
c...   CO2  --   728                              
      DATA (Qcoef(12,2,j),j=1,5)/-.73235E+04, .46140E+02,
     +               -.58473E-01, .57573E-04, 0.   /
c...    O3  --   666                              
      DATA (Qcoef(13,2,j),j=1,5)/-.11725E+05, .66515E+02,
     +               -.96010E-01, .94001E-04, 0.   /
c...    O3  --   668                              
      DATA (Qcoef(14,2,j),j=1,5)/-.25409E+05, .14393E+03,
     +               -.20850E+00, .20357E-03, 0.   /
c...    O3  --   686                              
      DATA (Qcoef(15,2,j),j=1,5)/-.12624E+05, .71391E+02,
     +               -.10383E+00, .10106E-03, 0.   /
c...    O3  --   667                              
      DATA (Qcoef(16,2,j),j=1,5)/-.14000E+06, .79825E+03,
     +               -.11465E+01, .11372E-02, 0.   /
c...    O3  --   676                              
      DATA (Qcoef(17,2,j),j=1,5)/-.69175E+05, .39442E+03,
     +               -.56650E+00, .56189E-03, 0.   /
c...   N2O  --   446                              
      DATA (Qcoef(18,2,j),j=1,5)/-.12673E+05, .75128E+02,
     +               -.10092E+00, .95557E-04, 0.   /
c...   N2O  --   456                              
      DATA (Qcoef(19,2,j),j=1,5)/-.90045E+04, .52833E+02,
     +               -.71771E-01, .67297E-04, 0.   /
c...   N2O  --   546                              
      DATA (Qcoef(20,2,j),j=1,5)/-.89960E+04, .53096E+02,
     +               -.71784E-01, .67592E-04, 0.   /
c...   N2O  --   448                              
      DATA (Qcoef(21,2,j),j=1,5)/-.13978E+05, .82338E+02,
     +               -.11167E+00, .10507E-03, 0.   /
c...   N2O  --   447                              
      DATA (Qcoef(22,2,j),j=1,5)/-.79993E+05, .47265E+03,
     +               -.63804E+00, .60218E-03, 0.   /
c...    CO  --    26                              
      DATA (Qcoef(23,2,j),j=1,5)/ .90723E+01, .33263E+00,
     +                .11806E-04, .27035E-07, 0.   /
c...    CO  --    36                              
      DATA (Qcoef(24,2,j),j=1,5)/ .20651E+02, .68810E+00,
     +                .34217E-04, .55823E-07, 0.   /
c...    CO  --    28                              
      DATA (Qcoef(25,2,j),j=1,5)/ .98497E+01, .34713E+00,
     +                .15290E-04, .28766E-07, 0.   /
c...    CO  --    27                              
      DATA (Qcoef(26,2,j),j=1,5)/ .58498E+02, .20351E+01,
     +                .87684E-04, .16554E-06, 0.   /
c...    CO  --    38                              
      DATA (Qcoef(27,2,j),j=1,5)/ .23511E+02, .71565E+00,
     +                .46681E-04, .58223E-07, 0.   /
c...    CO  --    37                              
      DATA (Qcoef(28,2,j),j=1,5)/ .11506E+03, .42727E+01,
     +                .17494E-03, .34413E-06, 0.   /
c...   CH4  --   211                              
      DATA (Qcoef(29,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   CH4  --   311                              
      DATA (Qcoef(30,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   CH4  --   212                              
      DATA (Qcoef(31,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...    O2  --    66                              
      DATA (Qcoef(32,2,j),j=1,5)/ .36539E+02, .57015E+00,
     +                .16332E-03, .45568E-07, 0.   /
c...    O2  --    68                              
      DATA (Qcoef(33,2,j),j=1,5)/ .77306E+02, .11818E+01,
     +                .38661E-03, .89415E-07, 0.   /
c...    O2  --    67                              
      DATA (Qcoef(34,2,j),j=1,5)/ .44281E+03, .69531E+01,
     +                .21669E-02, .53053E-06, 0.   /
c...    NO  --    46                              
      DATA (Qcoef(35,2,j),j=1,5)/-.78837E+02, .39173E+01,
     +                .80657E-03, .22042E-06, 0.   /
c...    NO  --    56                              
      DATA (Qcoef(36,2,j),j=1,5)/-.67000E+02, .27874E+01,
     +                .45181E-03, .21161E-06, 0.   /
c...    NO  --    48                              
      DATA (Qcoef(37,2,j),j=1,5)/-.98460E+02, .42347E+01,
     +                .71550E-03, .32213E-06, 0.   /
c...   SO2  --   626                              
      DATA (Qcoef(38,2,j),j=1,5)/-.21162E+05, .11846E+03,
     +               -.16648E+00, .16825E-03, 0.   /
c...   SO2  --   646                              
      DATA (Qcoef(39,2,j),j=1,5)/-.21251E+05, .11896E+03,
     +               -.16717E+00, .16895E-03, 0.   /
c...   NO2  --   646                              
      DATA (Qcoef(40,2,j),j=1,5)/-.27185E+05, .16489E+03,
     +               -.19540E+00, .22024E-03, 0.   /
c...   NH3  --  4111                              
      DATA (Qcoef(41,2,j),j=1,5)/-.47139E+03, .54035E+01,
     +                .64491E-02,-.72674E-06, 0.   /
c...   NH3  --  5111                              
      DATA (Qcoef(42,2,j),j=1,5)/-.31638E+03, .36086E+01,
     +                .43087E-02,-.48207E-06, 0.   /


c**********
c...  HNO3  --   146                              
c      DATA (Qcoef(43,2,j),j=1,5)/-.10000E+01, .00000E+00,
c     +                .00000E+00, .00000E+00, 0.   /

      DATA (Qcoef(43,2,j),j=1,5)/-2.155375E+8, 8.804353E+5,
     +              -1.165682E+3, 5.170532E-1, 0.   /

c 5.170532E-1*x^3 + -1.165682E+3*x^2 + 8.804353E+5*x + -2.155375E+8

c...    OH  --    61                              
      DATA (Qcoef(44,2,j),j=1,5)/-.88840E+01, .30202E+00,
     +               -.15565E-04, .14330E-07, 0.   /
c...    OH  --    81                              
      DATA (Qcoef(45,2,j),j=1,5)/-.51535E+01, .29076E+00,
     +               -.72340E-05, .20702E-07, 0.   /
c...    OH  --    62                              
      DATA (Qcoef(46,2,j),j=1,5)/-.41683E+02, .83890E+00,
     +               -.36063E-04, .38083E-07, 0.   /
c...    HF  --    19                              
      DATA (Qcoef(47,2,j),j=1,5)/-.36045E-01, .14220E+00,
     +               -.10755E-04, .65523E-08, 0.   /
c...   HCl  --    15                              
      DATA (Qcoef(48,2,j),j=1,5)/ .25039E+01, .54430E+00,
     +               -.38656E-04, .39793E-07, 0.   /
c...   HCl  --    17                              
      DATA (Qcoef(49,2,j),j=1,5)/ .14998E+01, .54847E+00,
     +               -.42209E-04, .41029E-07, 0.   /
c...   HBr  --    19                              
      DATA (Qcoef(50,2,j),j=1,5)/ .67229E+01, .66356E+00,
     +               -.33749E-04, .54818E-07, 0.   /
c...   HBr  --    11                              
      DATA (Qcoef(51,2,j),j=1,5)/ .67752E+01, .66363E+00,
     +               -.33655E-04, .54823E-07, 0.   /
c...    HI  --    17                              
      DATA (Qcoef(52,2,j),j=1,5)/ .29353E+02, .12220E+01,
     +                .10209E-04, .10719E-06, 0.   /
c...   ClO  --    56                              
      DATA (Qcoef(53,2,j),j=1,5)/ .22662E+03, .61093E+01,
     +                .14454E-01, .16928E-06, 0.   /
c...   ClO  --    76                              
      DATA (Qcoef(54,2,j),j=1,5)/ .23304E+03, .61805E+01,
     +                .14797E-01, .16629E-06, 0.   /
c...   OCS  --   622                              
      DATA (Qcoef(55,2,j),j=1,5)/-.54125E+04, .29749E+02,
     +               -.44698E-01, .38878E-04, 0.   /
c...   OCS  --   624                              
      DATA (Qcoef(56,2,j),j=1,5)/-.55472E+04, .30489E+02,
     +               -.45809E-01, .39847E-04, 0.   /
c...   OCS  --   632                              
      DATA (Qcoef(57,2,j),j=1,5)/-.11863E+05, .64745E+02,
     +               -.98318E-01, .84563E-04, 0.   /
c...   OCS  --   822                              
      DATA (Qcoef(58,2,j),j=1,5)/-.61288E+04, .33520E+02,
     +               -.50734E-01, .43792E-04, 0.   /
c...  H2CO  --   126                              
      DATA (Qcoef(59,2,j),j=1,5)/-.17628E+05, .91794E+02,
     +               -.13055E+00, .89336E-04, 0.   /
c...  H2CO  --   136                              
      DATA (Qcoef(60,2,j),j=1,5)/-.36151E+05, .18825E+03,
     +               -.26772E+00, .18321E-03, 0.   /
c...  H2CO  --   128                              
      DATA (Qcoef(61,2,j),j=1,5)/-.17628E+05, .91794E+02,
     +               -.13055E+00, .89336E-04, 0.   /
c...  HOCl  --   165                              
      DATA (Qcoef(62,2,j),j=1,5)/-.24164E+05, .15618E+03,
     +               -.13206E+00, .21900E-03, 0.   /
c...  HOCl  --   167                              
      DATA (Qcoef(63,2,j),j=1,5)/-.24592E+05, .15895E+03,
     +               -.13440E+00, .22289E-03, 0.   /
c...    N2  --    44                              
      DATA (Qcoef(64,2,j),j=1,5)/ .27907E+02, .14972E+01,
     +               -.70424E-05, .11734E-06, 0.   /
c...   HCN  --   124                              
      DATA (Qcoef(65,2,j),j=1,5)/-.78078E+03, .61725E+01,
     +               -.53816E-02, .73379E-05, 0.   /
c...   HCN  --   134                              
      DATA (Qcoef(66,2,j),j=1,5)/-.16309E+04, .12801E+02,
     +               -.11242E-01, .15268E-04, 0.   /
c...   HCN  --   125                              
      DATA (Qcoef(67,2,j),j=1,5)/-.56301E+03, .43794E+01,
     +               -.38928E-02, .52467E-05, 0.   /
c... CH3Cl  --   215                              
      DATA (Qcoef(68,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c... CH3Cl  --   217                              
      DATA (Qcoef(69,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...  H2O2  --  1661                              
      DATA (Qcoef(70,2,j),j=1,5)/-.27583E+05, .15064E+03,
     +               -.19917E+00, .16977E-03, 0.   /
c...  C2H2  --  1221                              
      DATA (Qcoef(71,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...  C2H2  --  1231                              
      DATA (Qcoef(72,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /

c**********
c...  C2H6  --  1221                              
      DATA (Qcoef(73,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   PH3  --  1111                              
      DATA (Qcoef(74,2,j),j=1,5)/-.28390E+05, .14463E+03,
     +               -.21473E+00, .14346E-03, 0.   /
c...  COF2  --   269                              
      DATA (Qcoef(75,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   SF6  --    29                              
      DATA (Qcoef(76,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   H2S  --   121                              
      DATA (Qcoef(77,2,j),j=1,5)/-.37572E+03, .29157E+01,
     +               -.98642E-03, .24113E-05, 0.   /
c...   H2S  --   141                              
      DATA (Qcoef(78,2,j),j=1,5)/-.37668E+03, .29231E+01,
     +               -.98894E-03, .24174E-05, 0.   /
c...   H2S  --   131                              
      DATA (Qcoef(79,2,j),j=1,5)/-.15049E+04, .11678E+02,
     +               -.39510E-02, .96579E-05, 0.   /
c... HCOOH  --   126                              
      DATA (Qcoef(80,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   HO2  --   166                              
      DATA (Qcoef(81,2,j),j=1,5)/-.32576E+04, .25539E+02,
     +               -.12803E-01, .29358E-04, 0.   /

c**********
c...     O  --     6                              
c      DATA (Qcoef(82,2,j),j=1,5)/-.10000E+01, .00000E+00,
c     +                 .00000E+00, .00000E+00, 0.   /

      DATA (Qcoef(82,2,j),j=1,5)/  1.0       , .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /

c...ClONO2  --  5646                              
      DATA (Qcoef(83,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...ClONO2  --  7646                              
      DATA (Qcoef(84,2,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   NO+  --    46                              
      DATA (Qcoef(85,2,j),j=1,5)/ .17755E+02, .99262E+00,
     +               -.70814E-05, .76699E-07, 0.   /
c
c...total internal partition sums for T>1500 <=3005 K range:
c...   H2O  --   161                              
      DATA (Qcoef( 1,3,j),j=1,5)/-.11727E+04, .29261E+01,
     +               -.13299E-02, .74356E-06, 0.   /
c...   H2O  --   181                              
      DATA (Qcoef( 2,3,j),j=1,5)/-.17914E+04, .39835E+01,
     +               -.19288E-02, .86144E-06, 0.   /
c...   H2O  --   171                              
      DATA (Qcoef( 3,3,j),j=1,5)/-.10665E+05, .23729E+02,
     +               -.11474E-01, .51294E-05, 0.   /
c...   H2O  --   162                              
      DATA (Qcoef( 4,3,j),j=1,5)/-.12585E+05, .26707E+02,
     +               -.14454E-01, .59457E-05, 0.   /
c...   CO2  --   626                              
      DATA (Qcoef( 5,3,j),j=1,5)/-.34938E+05, .66965E+02,
     +               -.44010E-01, .12662E-04, 0.   /
c...   CO2  --   636                              
      DATA (Qcoef( 6,3,j),j=1,5)/-.76420E+05, .14638E+03,
     +               -.96343E-01, .27589E-04, 0.   /
c...   CO2  --   628                              
      DATA (Qcoef( 7,3,j),j=1,5)/-.76677E+05, .14693E+03,
     +               -.96622E-01, .27746E-04, 0.   /
c...   CO2  --   627                              
      DATA (Qcoef( 8,3,j),j=1,5)/-.44040E+06, .84397E+03,
     +               -.55484E+00, .15946E-03, 0.   /
c...   CO2  --   638                              
      DATA (Qcoef( 9,3,j),j=1,5)/-.16856E+06, .32278E+03,
     +               -.21259E+00, .60747E-04, 0.   /
c...   CO2  --   637                              
      DATA (Qcoef(10,3,j),j=1,5)/-.96531E+06, .18487E+04,
     +               -.12172E+01, .34817E-03, 0.   /
c...   CO2  --   828                              
      DATA (Qcoef(11,3,j),j=1,5)/-.42074E+05, .80599E+02,
     +               -.53035E-01, .15202E-04, 0.   /
c...   CO2  --   728                              
      DATA (Qcoef(12,3,j),j=1,5)/-.48298E+06, .92535E+03,
     +               -.60873E+00, .17463E-03, 0.   /
c...    O3  --   666                              
      DATA (Qcoef(13,3,j),j=1,5)/-.61205E+06, .11896E+04,
     +               -.80924E+00, .24833E-03, 0.   /
c...    O3  --   668                              
      DATA (Qcoef(14,3,j),j=1,5)/-.13289E+07, .25826E+04,
     +               -.17574E+01, .53877E-03, 0.   /
c...    O3  --   686                              
      DATA (Qcoef(15,3,j),j=1,5)/-.66163E+06, .12857E+04,
     +               -.87521E+00, .26802E-03, 0.   /
c...    O3  --   667                              
      DATA (Qcoef(16,3,j),j=1,5)/-.70636E+07, .13772E+05,
     +               -.94024E+01, .29276E-02, 0.   /
c...    O3  --   676                              
      DATA (Qcoef(17,3,j),j=1,5)/-.34902E+07, .68051E+04,
     +               -.46459E+01, .14466E-02, 0.   /
c...   N2O  --   446                              
      DATA (Qcoef(18,3,j),j=1,5)/-.83406E+06, .15951E+04,
     +               -.10534E+01, .29849E-03, 0.   /
c...   N2O  --   456                              
      DATA (Qcoef(19,3,j),j=1,5)/-.59281E+06, .11334E+04,
     +               -.74907E+00, .21164E-03, 0.   /
c...   N2O  --   546                              
      DATA (Qcoef(20,3,j),j=1,5)/-.59301E+06, .11339E+04,
     +               -.74918E+00, .21193E-03, 0.   /
c...   N2O  --   448                              
      DATA (Qcoef(21,3,j),j=1,5)/-.92317E+06, .17651E+04,
     +               -.11664E+01, .32984E-03, 0.   /
c...   N2O  --   447                              
      DATA (Qcoef(22,3,j),j=1,5)/-.52739E+07, .10085E+05,
     +               -.66623E+01, .18858E-02, 0.   /
c...    CO  --    26                              
      DATA (Qcoef(23,3,j),j=1,5)/ .63418E+02, .20760E+00,
     +                .10895E-03, .19844E-08, 0.   /
c...    CO  --    36                              
      DATA (Qcoef(24,3,j),j=1,5)/ .13265E+03, .43434E+00,
     +                .22794E-03, .41523E-08, 0.   /
c...    CO  --    28                              
      DATA (Qcoef(25,3,j),j=1,5)/ .66581E+02, .21800E+00,
     +                .11441E-03, .20839E-08, 0.   /
c...    CO  --    27                              
      DATA (Qcoef(26,3,j),j=1,5)/ .39033E+03, .12780E+01,
     +                .67066E-03, .12218E-07, 0.   /
c...    CO  --    38                              
      DATA (Qcoef(27,3,j),j=1,5)/ .13959E+03, .45717E+00,
     +                .23991E-03, .43712E-08, 0.   /
c...    CO  --    37                              
      DATA (Qcoef(28,3,j),j=1,5)/ .81756E+03, .26767E+01,
     +                .14046E-02, .25378E-07, 0.   /
c...   CH4  --   211                              
      DATA (Qcoef(29,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   CH4  --   311                              
      DATA (Qcoef(30,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   CH4  --   212                              
      DATA (Qcoef(31,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...    O2  --    66                              
      DATA (Qcoef(32,3,j),j=1,5)/ .76324E+01, .58006E+00,
     +                .18941E-03, .27822E-07, 0.   /
c...    O2  --    68                              
      DATA (Qcoef(33,3,j),j=1,5)/ .16157E+02, .12282E+01,
     +                .40112E-03, .58919E-07, 0.   /
c...    O2  --    67                              
      DATA (Qcoef(34,3,j),j=1,5)/ .94397E+02, .71717E+01,
     +                .23423E-02, .34425E-06, 0.   /
c...    NO  --    46                              
      DATA (Qcoef(35,3,j),j=1,5)/ .52033E+03, .26381E+01,
     +                .17177E-02, .17933E-07, 0.   /
c...    NO  --    56                              
      DATA (Qcoef(36,3,j),j=1,5)/ .35655E+03, .18125E+01,
     +                .12126E-02, .12110E-07, 0.   /
c...    NO  --    48                              
      DATA (Qcoef(37,3,j),j=1,5)/ .54190E+03, .27581E+01,
     +                .18694E-02, .18266E-07, 0.   /
c...   SO2  --   626                              
      DATA (Qcoef(38,3,j),j=1,5)/-.10718E+07, .20831E+04,
     +               -.14138E+01, .43806E-03, 0.   /
c...   SO2  --   646                              
      DATA (Qcoef(39,3,j),j=1,5)/-.10762E+07, .20918E+04,
     +               -.14196E+01, .43988E-03, 0.   /
c...   NO2  --   646                              
      DATA (Qcoef(40,3,j),j=1,5)/-.12837E+07, .25067E+04,
     +               -.16761E+01, .53904E-03, 0.   /
c...   NH3  --  4111                              
      DATA (Qcoef(41,3,j),j=1,5)/-.17334E+04, .80988E+01,
     +                .44771E-02,-.24084E-06, 0.   /
c...   NH3  --  5111                              
      DATA (Qcoef(42,3,j),j=1,5)/-.11656E+04, .54254E+01,
     +                .29809E-02,-.15750E-06, 0.   /

c**********
c...  HNO3  --   146                              
c      DATA (Qcoef(43,3,j),j=1,5)/-.10000E+01, .00000E+00,
c     +                .00000E+00, .00000E+00, 0.   /

      DATA (Qcoef(43,3,j),j=1,5)/-1.425802E+11, 2.311441E+8,
     +             -1.248681E+5,  2.268254E+1 , 0.  /

c 2.268254E+1*x^3 + -1.248681E+5*x^2 + 2.311441E+8*x + -1.425802E+11

c...    OH  --    61                              
      DATA (Qcoef(44,3,j),j=1,5)/ .33750E+02, .22130E+00,
     +                .35953E-04, .32366E-08, 0.   /
c...    OH  --    81                              
      DATA (Qcoef(45,3,j),j=1,5)/ .42716E+02, .18843E+00,
     +                .67175E-04, .23903E-08, 0.   /
c...    OH  --    62                              
      DATA (Qcoef(46,3,j),j=1,5)/ .72913E+02, .62430E+00,
     +                .99073E-04, .94503E-08, 0.   /
c...    HF  --    19                              
      DATA (Qcoef(47,3,j),j=1,5)/ .18423E+02, .10799E+00,
     +                .10568E-04, .20752E-08, 0.   /
c...   HCl  --    15                              
      DATA (Qcoef(48,3,j),j=1,5)/ .92445E+02, .35539E+00,
     +                .96272E-04, .71602E-08, 0.   /
c...   HCl  --    17                              
      DATA (Qcoef(49,3,j),j=1,5)/ .92519E+02, .35592E+00,
     +                .96492E-04, .71775E-08, 0.   /
c...   HBr  --    19                              
      DATA (Qcoef(50,3,j),j=1,5)/ .11692E+03, .42161E+00,
     +                .14690E-03, .92595E-08, 0.   /
c...   HBr  --    11                              
      DATA (Qcoef(51,3,j),j=1,5)/ .11700E+03, .42161E+00,
     +                .14703E-03, .92525E-08, 0.   /
c...    HI  --    17                              
      DATA (Qcoef(52,3,j),j=1,5)/ .22138E+03, .78595E+00,
     +                .34579E-03, .20348E-07, 0.   /
c...   ClO  --    56                              
      DATA (Qcoef(53,3,j),j=1,5)/ .37348E+03, .56800E+01,
     +                .14805E-01, .23168E-06, 0.   /
c...   ClO  --    76                              
      DATA (Qcoef(54,3,j),j=1,5)/ .38100E+03, .57530E+01,
     +                .15142E-01, .23652E-06, 0.   /
c...   OCS  --   622                              
      DATA (Qcoef(55,3,j),j=1,5)/-.37301E+06, .71169E+03,
     +               -.47328E+00, .13049E-03, 0.   /
c...   OCS  --   624                              
      DATA (Qcoef(56,3,j),j=1,5)/-.38232E+06, .72945E+03,
     +               -.48509E+00, .13375E-03, 0.   /
c...   OCS  --   632                              
      DATA (Qcoef(57,3,j),j=1,5)/-.82204E+06, .15682E+04,
     +               -.10435E+01, .28668E-03, 0.   /
c...   OCS  --   822                              
      DATA (Qcoef(58,3,j),j=1,5)/-.42390E+06, .80869E+03,
     +               -.53803E+00, .14798E-03, 0.   /
c...  H2CO  --   126                              
      DATA (Qcoef(59,3,j),j=1,5)/-.24906E+07, .45519E+04,
     +               -.28336E+01, .64198E-03, 0.   /
c...  H2CO  --   136                              
      DATA (Qcoef(60,3,j),j=1,5)/-.51075E+07, .93349E+04,
     +               -.58110E+01, .13165E-02, 0.   /
c...  H2CO  --   128                              
      DATA (Qcoef(61,3,j),j=1,5)/-.24906E+07, .45519E+04,
     +               -.28336E+01, .64198E-03, 0.   /
c...  HOCl  --   165                              
      DATA (Qcoef(62,3,j),j=1,5)/-.11326E+07, .22197E+04,
     +               -.14357E+01, .49952E-03, 0.   /
c...  HOCl  --   167                              
      DATA (Qcoef(63,3,j),j=1,5)/-.11527E+07, .22590E+04,
     +               -.14612E+01, .50838E-03, 0.   /
c...    N2  --    44                              
      DATA (Qcoef(64,3,j),j=1,5)/ .27986E+03, .93070E+00,
     +                .42409E-03, .95573E-08, 0.   /
c...   HCN  --   124                              
      DATA (Qcoef(65,3,j),j=1,5)/-.51989E+05, .10057E+03,
     +               -.64310E-01, .19844E-04, 0.   /
c...   HCN  --   134                              
      DATA (Qcoef(66,3,j),j=1,5)/-.10838E+06, .20960E+03,
     +               -.13411E+00, .41345E-04, 0.   /
c...   HCN  --   125                              
      DATA (Qcoef(67,3,j),j=1,5)/-.37363E+05, .72225E+02,
     +               -.46251E-01, .14237E-04, 0.   /
c... CH3Cl  --   215                              
      DATA (Qcoef(68,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c... CH3Cl  --   217                              
      DATA (Qcoef(69,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...  H2O2  --  1661                              
      DATA (Qcoef(70,3,j),j=1,5)/-.26863E+07, .49815E+04,
     +               -.31584E+01, .78351E-03, 0.   /
c...  C2H2  --  1221                              
      DATA (Qcoef(71,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...  C2H2  --  1231                              
      DATA (Qcoef(72,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /

c**********
c...  C2H6  --  1221                              
      DATA (Qcoef(73,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /

c...   PH3  --  1111                              
      DATA (Qcoef(74,3,j),j=1,5)/-.44074E+07, .80563E+04,
     +               -.50179E+01, .11272E-02, 0.   /
c...  COF2  --   269                              
      DATA (Qcoef(75,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   SF6  --    29                              
      DATA (Qcoef(76,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   H2S  --   121                              
      DATA (Qcoef(77,3,j),j=1,5)/-.10043E+05, .20827E+02,
     +               -.12249E-01, .48236E-05, 0.   /
c...   H2S  --   141                              
      DATA (Qcoef(78,3,j),j=1,5)/-.10069E+05, .20881E+02,
     +               -.12280E-01, .48359E-05, 0.   /
c...   H2S  --   131                              
      DATA (Qcoef(79,3,j),j=1,5)/-.40225E+05, .83420E+02,
     +               -.49061E-01, .19320E-04, 0.   /
c... HCOOH  --   126                              
      DATA (Qcoef(80,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   HO2  --   166                              
      DATA (Qcoef(81,3,j),j=1,5)/-.13056E+06, .26188E+03,
     +               -.16161E+00, .61250E-04, 0.   /

c**********
c...     O  --     6                              
c      DATA (Qcoef(82,3,j),j=1,5)/-.10000E+01, .00000E+00,
c     +                .00000E+00, .00000E+00, 0.   /

      DATA (Qcoef(82,3,j),j=1,5)/  1.0       , .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /

c...ClONO2  --  5646                              
      DATA (Qcoef(83,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...ClONO2  --  7646                              
      DATA (Qcoef(84,3,j),j=1,5)/-.10000E+01, .00000E+00,
     +                .00000E+00, .00000E+00, 0.   /
c...   NO+  --    46                              
      DATA (Qcoef(85,3,j),j=1,5)/ .18634E+03, .61771E+00,
     +                .27607E-03, .45828E-08, 0.   /
c     
c                    H2O 161,         181,         171,                         
      DATA q296/ .174626E+03, .176141E+03, .105306E+04,
c               162;     CO2 626,         636,         628,         627,        
     +  .865122E+03, .286219E+03, .576928E+03, .607978E+03, .354389E+04,
c               638,         637,         828,         728;      O3 666,        
     +  .123528E+04, .714432E+04, .323407E+03, .376700E+04, .348186E+04,
c               668,         686,         667,         676;     N2O 446,        
     +  .746207E+04, .364563E+04, .430647E+05, .212791E+05, .499183E+04,
c               456,         546,         448,         447;       CO 26,        
     +  .334938E+04, .344940E+04, .526595E+04, .307008E+05, .107428E+03,
c                36,          28,          27,          38,          37;        
     +  .224704E+03, .112781E+03, .661209E+03, .236447E+03, .138071E+04,
c           CH4 211,         311,         212;       O2 66,          68,        
     +  .589908E+03, .117974E+04, .477061E+04, .215726E+03, .452188E+03,
c                67;      NO  46,          56,          48;     SO2 626,        
     +  .263998E+04, .113243E+04, .785200E+03, .119417E+04, .634449E+04,
c               646;     NO2 646;    NH3 4111,        5111;    HNO3 146;        
     +  .637321E+04, .136318E+05, .171089E+04, .114134E+04, .213822E+06,
c             OH 61,         81,           62;       HF 19;      HCl 15,        
     +  .803295E+02, .808460E+02, .207108E+03, .414625E+02, .160650E+03,
c                17;      HBr 19,          11;       HI 17;      ClO 56,        
     +  .160887E+03, .200165E+03, .200227E+03, .388948E+03, .328810E+04,
c                76;     OCS 622,         624,         632,         822;        
     +  .334560E+04, .121746E+04, .124793E+04, .247482E+04, .130948E+04,
c          H2CO 126,         136,         128;    HOCl 165,         167;        
     +  .268388E+04, .550322E+04, .284573E+04, .193166E+05, .196584E+05,
c             N2 44;     HCN 124,         134,         125;   CH3Cl 215,        
     +  .467136E+03, .893323E+03, .183657E+04, .615046E+03, .144858E+05,
c               217;   H2O2 1661;   C2H2 1221,        1231,   C2H6 1221;        
     +  .147153E+05, .767871E+04, .412519E+03, .330014E+04, .546265E+05,
c          PH3 1111;    COF2 269;      SF6 29;     H2S 121,         141,        
     +  .325067E+04, .697632E+05, .162242E+07, .503204E+03, .504486E+03,
c               131;   HCOOH 126;     HO2 166;        O 6;  ClONO2 5646,        
     +  .201546E+04, .389257E+05, .430184E+04,-.100000E+01, .212829E+07,
c              7646;      NO+ 46;                                               
     +  .218246E+07, .308855E+03/
      END

      BLOCK DATA BDMOL                                                    C04000
C                                                                         C04010
C     LAST MODIFIED JANUARY 17, 1991                                      C04020
C                                                                         C04030
      PARAMETER (NTMOL=36,NSPECI=85)                                      C04040
      CHARACTER*6 MOLID                                                   C04050
      COMMON /MOLNAM/ MOLID(0:NTMOL)                                      C04060
C                                                                         C04070
      DATA (MOLID(I),I=0,NTMOL)/   '  ALL ',                              C04080
     *  '  H2O ','  CO2 ','   O3 ','  N2O ','   CO ','  CH4 ','   O2 ',   C04090
     *  '   NO ','  SO2 ','  NO2 ','  NH3 ',' HNO3 ','   OH ','   HF ',   C04100
     *  '  HCL ','  HBR ','   HI ','  CLO ','  OCS ',' H2CO ',' HOCL ',   C04110
     *  '   N2 ','  HCN ','CH3CL ',' H2O2 ',' C2H2 ',' C2H6 ','  PH3 ',   C04120
     *  ' COF2 ','  SF6 ','  H2S ','HCOOH ','  HO2 ','    O ','CLONO2',
     *  '  NO+ ' /
C                                                                         C04140
      END                                                                 C04150
      SUBROUTINE R1PRNT (V1P,DVP,NLIM,R1,JLO,NPTS,MFILE,IENTER)           C12550
C                                                                         C12560
      IMPLICIT REAL*8           (V)                                     ! C12570
C                                                                         C12580
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         C12600
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        C12610
     *              NLTEFL,LNFIL4,LNGTH4                                  C12620
      DIMENSION R1(*)                                                     C12630
C                                                                         C12640
C     THIS SUBROUTINE PRINTS THE FIRST NPTS VALUES STARTING AT JLO        C12650
C     AND THE LAST NPTS VALUES ENDING AT NLIM OF THE R1 ARRAY             C12660
C                                                                         C12670
      IF (NPTS.LE.0) RETURN                                               C12680
      IF (IENTER.LT.1) WRITE (IPR,900)                                    C12690
      WRITE (IPR,905)                                                     C12700
      IENTER = IENTER+1                                                   C12710
      JHI = JLO+NLIM-1                                                    C12720
      NNPTS = NPTS                                                        C12730
      IF (NPTS.GT.(NLIM/2)+1) NNPTS = (NLIM/2)+1                          C12740
      JLOLIM = JLO+NNPTS-1                                                C12750
      JHILIM = JHI-NNPTS+1                                                C12760
      DO 10 KK = 1, NNPTS                                                 C12770
         J = JLO+KK-1                                                     C12780
         VJ = V1P+FLOAT(J-JLO)*DVP                                        C12790
         IK = JHILIM+KK-1                                                 C12800
         VK = V1P+FLOAT(IK-JLO)*DVP                                       C12810
         WRITE (IPR,910) J,VJ,R1(J),IK,VK,R1(IK)                          C12820
   10 CONTINUE                                                            C12830
C                                                                         C12840
      RETURN                                                              C12850
C                                                                         C12860
  900 FORMAT ('0 ','LOCATION  WAVENUMBER',2X,'OPT. DEPTH',27X,            C12870
     *        'LOCATION   WAVENUMBER',2X,'OPT. DEPTH')                    C12880
  905 FORMAT (' ')                                                        C12890
  910 FORMAT (I8,2X,F12.6,1P,E15.7,0P,20X,I8,2X,F12.6,1P,E15.7)           C12900
C                                                                         C12910
      END                                                                 C12920
      SUBROUTINE LINF4 (V1L4,V2L4)                                        D00010
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
      REAL*8             HID,HMOLIL,HID1,HLINHD                         & D00130
C                                                                         D00140
      COMMON /BUFID/ HID(10),HMOLIL(64),MOLCNT(64),MCNTLC(64),            D00150
     *               MCNTNL(64),SUMSTR(64),NMOI,FLINLO,FLINHI,            D00160
     *               ILIN,ILINLC,ILINNL,IREC,IRECTL,HID1(2),LSTWDL        D00170
C                                                                         D00180
      COMMON VNU(1250),SP(1250),ALFA0(1250),EPP(1250),MOL(1250),          D00190
     *       SPP(1250)                                                    D00200
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
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           D00340
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
C                                                                         D00440
      REAL L4TIM,L4TMR,L4TMS,LOTHER
      DIMENSION MEFDP(64)                                                 D00450
      DIMENSION SCOR(NSPECI),RHOSLF(NSPECI),ALFD1(NSPECI)                 D00460
      DIMENSION ALFAL(1250),ALFAD(1250),A(4),B(4),TEMPLC(4)               D00470
      DIMENSION RCDHDR(2),IWD(2),IWD3(2),HLINHD(2),AMOLB(250)             D00480
C                                                                         D00490
      EQUIVALENCE (ALFA0(1),ALFAL(1)) , (EPP(1),ALFAD(1))                 D00500
      EQUIVALENCE (IHIRAC,FSCDID(1)) , (ILBLF4,FSCDID(2))                 D00510
      EQUIVALENCE (VNULO,RCDHDR(1)) , (IWD3(1),VD),                       D00520
     *            (HLINHD(1),HID(1),IWD(1)) , (MOLB(1),AMOLB(1))          D00530
                                                                          D00540
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
      CALL SHRINK                                                         D02260
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
  905 FORMAT ('0*************************',I5,' STRENGTHS FOR',           D02580
     *        '  TRANSITIONS WITH UNKNOWN EPP FOR MOL =',I5,              D02590
     *        ' SET TO ZERO')                                             D02600
  910 FORMAT ('0',20X,'TIME',11X,'READ',9X,'SHRINK',6X,'NO. LINES',3X,    D02610
     *        'AFTER SHRINK',/,2X,'LINF4 ',2X,3F15.3,2I15)                D02620
C                                                                         D02630
      END                                                                 D02640
c***********************************************************************
      subroutine lnfilhd_4(linfil,linfil4,v1,v2)
c
c     this subroutine buffers past the file header for LINF4
c
      IMPLICIT REAL*8           (V) 
c
      CHARACTER*8      HLINID,BMOLID,HID1
      CHARACTER*1 CNEGEPP(8)
C                                                                         A09600
      integer *4 molcnt,mcntlc,
     *           mcntnl,linmol,
     *           lincnt,ilinlc,ilinnl,irec,irectl
c
      COMMON /LINHDR/ HLINID(10),BMOLID(64),MOLCNT(64),MCNTLC(64),        A09610
     *                MCNTNL(64),SUMSTR(64),LINMOL,FLINLO,FLINHI,         A09620
     *                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1(2),LSTWDL     A09630
      common /bufid2/ n_negepp(64),n_resetepp(64),xspace(4096),lstwdl2

      COMMON /IFIL/ Idum1,IPR,Idum2,Ndum1,Ndum2F,Ndum3,Ndum4,Ndum5,
     *              Ndum6,Kdum1,Kdum2,Ldum1,Ndum7,Idum3,Idum4,
     *              Ndum8,Ldum2,Ldum3
      common /eppinfo/ negepp_flag

      real *4 sumstr,flinlo,flinhi

      integer *4 negepp_flag,n_negepp,n_resetepp
      real *4 xspace
C                                                                         A09640
      lnfil = linfil
      lnfil4= linfil4

c
      REWIND lnfil
c
      read (lnfil,end=777)    HLINID,BMOLID,MOLCNT,MCNTLC,       
     *                MCNTNL,SUMSTR,LINMOL,FLINLO,FLINHI,  
     *                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1
c
      READ (HLINID(7),950) CNEGEPP
      IF (CNEGEPP(8).eq.'^') THEN
         read (lnfil) n_negepp,n_resetepp,xspace
      endif

      go to 5
c
 777  STOP 'Linf4: LINFIL DOES NOT EXIST'                   
c
 5    continue

C     
      IF (V1.GT.FLINHI.OR.V2.LT.FLINLO) THEN                              D00940
         CALL ENDFIL_4 (LNFIL4)
         WRITE (IPR,900) V1,V2,FLINLO,FLINHI                              D00960
         RETURN                                                           D00970
      ENDIF                                                               D00980
c
      write (lnfil4)    HLINID,BMOLID,MOLCNT,MCNTLC,       
     *                MCNTNL,SUMSTR,LINMOL,FLINLO,FLINHI,  
     *                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1

      IF (CNEGEPP(8).eq.'^') THEN
         write (lnfil4) n_negepp,n_resetepp,xspace
         negepp_flag = 1
      endif
C
      return
c
  900 FORMAT ('0  *****  LINF4 - VNU LIMITS DO NOT INTERSECT WITH ',      D02540
     *        'LINFIL - LINF4 NOT USED *****',/,'   VNU = ',F10.3,        D02550
     *        ' - ',F10.3,' CM-1     LINFIL = ',F10.3,' - ',F10.3,        D02560
     *        ' CM-1')                                                    D02570
 950  FORMAT (8a1)
c
      end
c*******************************************************************
      SUBROUTINE RDLNFL (IEOF,ILO,IHI)                                    D02650
C                                                                         D02660
      IMPLICIT REAL*8           (V)                                     ! D02670
C                                                                         D02680
C     SUBROUTINE RDLNFL INPUTS THE LINE DATA FROM LINFIL                  D02690
C                                                                         D02700
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SEC   ,       XALTZ
C                                                                         D02720
      COMMON /FILHDR/ XID(10),SEC   ,PAV ,TAV ,HMOLID(60),XALTZ(4),       D02730
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   D02740
     *                EMIS ,FSCDID(17),NMOL,LAYRS ,YID1,YID(10),LSTWDF    D02750
      COMMON /R4SUB/ VLO,VHI,ILD,IST,IHD,LIMIN,LIMOUT,ILAST,DPTMN,        D02760
     *               DPTFC,ILIN4,ILIN4T                                   D02770
      COMMON /BUFR/ VNUB(250),SB(250),ALB(250),EPPB(250),MOLB(250),       D02790
     *              HWHMB(250),TMPALB(250),PSHIFB(250),IFLG(250)          D02800
c
      dimension amolb(250)
      equivalence (molb(1),amolb(1))
c
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         D02810
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        D02820
     *              NLTEFL,LNFIL4,LNGTH4                                  D02830
      COMMON /IOU/ IOUT(250)                                              D02840
c
      common /rdlnpnl/ vmin,vmax,nrec,nwds
      integer *4 nrec,nwds,lnfl,leof,npnlhd
c
      common /rdlnbuf/ vlin(250),str(250),hw_f(250),e_low(250),
     *     mol_id(250),hw_s(250),hw_T(250),shft(250),jflg(250)

      real *4 str,hw_f,e_low,hw_s,hw_T,shft,rdpnl(2),dum(2),xmol(2)
      integer *4 mol_id,jflg

      equivalence (vmin,rdpnl(1)), (mol_id(1),xmol(1))
c
      IPASS = 1                                                           D02890
      IF (ILO.GT.0) IPASS = 2                                             D02900
C                                                                         D02910
      ILNGTH = NLNGTH*250                                                 D02920
C                                                                         D02930
      IEOF = 0                                                            D02940
      ILO = 1                                                             D02950
      IHI = 0                                                             D02960
c
      lnfl = linfil	
      npnlhd = 6
c
   10 CALL BUFINln (Lnfl,LEOF,rdpnl(1),npnlhd) 
c
      IF (LEOF.EQ.0) GO TO 30                                             D02980
      IF (VMAX.LT.VLO) THEN                                               D02990
         CALL BUFINln (lnfl,LEOF,dum(1),1) 
         GO TO 10                                                         D03010
      ELSE                                                                D03020
         CALL BUFINln (Lnfl,LEOF,vlin(1),NWDS) 
      ENDIF                                                               D03040
c
      IF ((IPASS.EQ.1).AND.(Vlin(1).GT.VLO)) WRITE (IPR,900)
C                                                                         D03060
      IJ = 0                                                              D03070
C                                                                         D03080
c     precision conversion occurs here:
c     incoming on right: vlin is real*8;  others are real*4 and integer*4
c
      do 15 i=1,nrec
         IFLG(i)  = jflg(i)     
         VNUB(i)   = vlin(i)
         SB(i)    = str(i)   
         ALB(i)   = hw_f(i)      
         EPPB(i)   = e_low(i)    
         if (iflg(i) .ge.  0) then
            MOLB(i)  = mol_id(i)    
         else
            amolb(i)  = xmol(i)
         endif
         HWHMB(i) = hw_s(i)      
         TMPALB(i)= hw_T(i)       
         PSHIFB(i)= shft(i)       
 15   continue
c
      DO 20 J = 1, NREC                                                   D03090
         IF (IFLG(J).GE.0) THEN                                           D03100
            IJ = IJ+1                                                     D03110
            IOUT(IJ) = J                                                  D03120
         ENDIF                                                            D03130
   20 CONTINUE                                                            D03140
      IHI = IJ                                                            D03150
      RETURN                                                              D03160
   30 IF (NOPR.EQ.0) WRITE (IPR,905)                                      D03170
      IEOF = 1                                                            D03180
      RETURN                                                              D03190
C                                                                         D03200
  900 FORMAT ('0 FIRST LINE USED IN RDLNFL--- CHECK THE LINEFIL  ')       D03210
  905 FORMAT ('0 EOF ON LINFIL IN RDLNFL -- CHECK THE LINFIL ')           D03220
C                                                                         D03230
      END                                                                 D03240
      SUBROUTINE SHRINK                                                   D03250
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
c
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
      SUBROUTINE LBLF4 (JRAD,V1,V2)                                       D04440
C                                                                         D04450
      IMPLICIT REAL*8           (V)                                     ! D04460
C                                                                         D04470
C     SUBROUTINE LBLF4 DOES A LINE BY LINE CALCULATION                    D04480
C     USING FUNCTION F4.                                                  D04490
C                                                                         D04500
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           D04510
      COMMON /BUF/ VNU(1250),S(1250),ALFAL(1250),ALFAD(1250),MOL(1250),   D04520
     *             SPP(1250)                                              D04530
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
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           D04650
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
   10 CONTINUE                                                            D05020
      BETA = RADCN2/TAVE                                                  D05030
      VLO = V1R4-BOUND4                                                   D05040
      VHI = V2R4+BOUND4                                                   D05050
   20 CALL CPUTIM (TIM0)                                                  D05060
      CALL RDLIN4 (IEOF)                                                  D05070
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
      CALL CONVF4 (VNU,S,ALFAL,ALFAD,MOL,SPP)                             D05190
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
            RADVI = RADVI+RDEL                                            D05530
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
      SUBROUTINE RDLIN4 (IEOF)                                            D05660
C                                                                         D05670
      IMPLICIT REAL*8           (V)                                     ! D05680
C                                                                         D05690
C     SUBROUTINE RDLIN4 INPUTS THE LINE DATA FROM LNFIL4                  D05700
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
         if (negepp_flag.eq.1) CALL BUFIN (LNFIL4,LEOF,DUM(1),1)
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
 950  FORMAT (8a1)

      RETURN                                                              D06240
C                                                                         D06250
      END                                                                 D06260
      SUBROUTINE CONVF4 (VNU,S,ALFAL,ALFAD,MOL,SPP)                       D06270
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
C                                                                         D06390
      DIMENSION VNU(*),S(*),ALFAL(*),ALFAD(*),MOL(*),SPP(*)               D06400
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
      A3 = (1.+2.*ZSQBND)/(1.+ZSQBND)**2                                  D06810
      B3 = -1./(1.+ZSQBND)**2                                             D06820
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
         IF (S(I).EQ.0..AND.SPP(I).EQ.0.) GO TO 60                        D06920
         ALFADI = ALFAD(I)                                                D06930
         ALFALI = ALFAL(I)                                                D06940
         ZETAI = ALFALI/(ALFALI+ALFADI)                                   D06950
         IZ = 100.*ZETAI + ONEPL                                          D06960
         ZETDIF = 100.*ZETAI - FLOAT(IZ-1)
         ALFAVI = ( AVRAT(IZ) + ZETDIF*(AVRAT(IZ+1)-AVRAT(IZ)) ) *        D06970
     x            (ALFALI+ALFADI)
         RALFVI = 1./ALFAVI                                               D06980
         SIL = S(I)*RECPI*ALFALI                                          D06990
         SIV = (ALFALI*RALFVI)*S(I)*RECPI*RALFVI                          D07000
C                                                                         D07010
         IF (SPP(I).EQ.0.) THEN                                           D07020
            SPEAK = A3*(ABS(SIV))                                         D07030
         ELSE                                                             D07040
            SILX = SPP(I)*RECPI                                           D07050
            SIVX = ((ALFALI*RALFVI)*SPP(I)*RECPI*RALFVI)*RALFVI           D07060
            SPEAK = A3*(ABS(SIV)+ABS(SIVX))                               D07070
         ENDIF                                                            D07080
C                                                                         D07090
         JJ = (VNU(I)-V1R4)/DVR4+1.                                       D07100
         JJ = MAX(JJ,1)                                                   D07110
         JJ = MIN(JJ,NPTR4)                                               D07120
C
         IF (ILNFLG.LE.1) THEN
            FREJ(I) = HNOREJ
            IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ))) THEN
               FREJ(I) = HREJ
               GO TO 60                                                   D07130
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
         F4BND = SIL/(ALFLI2+BNDSQ)                                       D07320
         IF (SPP(I).NE.0.) F4BNDX = SILX/(ALFLI2+BNDSQ)                   D07330
C                                                                         D07340
C                FOURTH FUNCTION CONVOLUTION                              D07350
C                                                                         D07360
         DO 40 JJ = JMIN, JMAX                                            D07370
            XM = (XJJ-XNUI)                                               D07380
            XMSQ = XM*XM                                                  D07390
            ZVSQ = XMSQ/ALFVI2                                            D07400
C                                                                         D07410
            IF (ZVSQ.LE.ZSQBND) THEN                                      D07420
               F4FN = SIV*(A3+ZVSQ*B3)-F4BND                              D07430
               IF (SPP(I).NE.0.)                                          D07440
     *             F4FN = F4FN+XM*(SIVX*(A3+ZVSQ*B3)-F4BNDX)              D07450
            ELSE                                                          D07460
               F4FN = SIL/(ALFLI2+XMSQ)-F4BND                             D07470
               IF (SPP(I).NE.0.)                                          D07480
     *             F4FN = F4FN+XM*(SILX/(ALFLI2+XMSQ)-F4BNDX)             D07490
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
            ELSE                                                          D07600
               R4(JJ) = R4(JJ)+F4FN                                       D07610
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
            SIVX = -SIVX                                                  D07730
            SILX = -SILX                                                  D07740
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
      SUBROUTINE XSREAD (XV1,XV2)                                         E00010
C                                                                         E00020
      IMPLICIT REAL*8           (V)                                     ! E00030
C                                                                         E00040
C**********************************************************************   E00050
C     THIS SUBROUTINE READS IN THE DESIRED "CROSS-SECTION"                E00060
C     MOLECULES WHICH ARE THEN MATCHED TO THE DATA CONTAINED              E00070
C     ON INPUT FILE FSCDXS.                                               E00080
C**********************************************************************   E00090
C                                                                         E00100
C     IFIL CARRIES FILE INFORMATION                                       E00110
C                                                                         E00120
      PARAMETER (MXFSC=200, MXLAY=MXFSC+3,MXZMD=3400,
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=35,MXTRAC=22)
C
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         E00130
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        E00140
     *              NLTEFL,LNFIL4,LNGTH4                                  E00150
C                                                                         E00160
C     IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE     E00170
C     MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES       E00180
C     (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR     E00190
C     L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD         E00200
C     MOLECULES.                                                          E00210
C                                                                         E00220
      COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(35),XAMNT(35,MXLAY)              E00230
C                                                                         E00240
C     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES         E00250
C     FOR THE CROSS-SECTION MOLECULES.                                    E00260
C     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES          E00270
C                                                                         E00280
      CHARACTER*10 XSFILE,XSNAME,ALIAS,XNAME,XFILS(6),BLANK               E00290
      COMMON /XSECTF/ XSFILE(6,5,35),XSNAME(35),ALIAS(4,35)               E00300
      COMMON /XSECTR/ V1FX(5,35),V2FX(5,35),DVFX(5,35),WXM(35),           E00310
     *                NTEMPF(5,35),NSPECR(35),IXFORM(5,35),               E00320
     *                XSMASS(35),XDOPLR(5,35),NUMXS,IXSBIN                E00325
C                                                                         E00330
      DIMENSION IXFLG(35)                                                 E00340
C                                                                         E00350
      CHARACTER*120 XSREC                                                 E00360
      CHARACTER*1 CFLG,CASTSK,CPRCNT,CFRM,CN,CF                           E00370
      EQUIVALENCE (CFLG,XSREC)                                            E00380
C                                                                         E00390
      DATA CASTSK / '*'/,CPRCNT / '%'/,CN / 'N'/,CF / 'F'/                E00400
      DATA BLANK / '          '/                                          E00410
C                                                                         E00411
C     T296 IS TEMPERATURE FOR INITAL CALCULATIN OF DOPPLER WIDTHS         E00412
C                                                                         E00413
      DATA T296 / 296.0 /                                                 E00414
C                                                                         E00420
      IXMAX = 35                                                          E00430
      DO 10 I = 1, IXMAX                                                  E00440
         XSNAME(I) = BLANK                                                E00450
   10 CONTINUE                                                            E00460
C                                                                         E00470
C     READ IN THE NAMES OF THE MOLECULES                                  E00480
C                                                                         E00490
      IF (IXMOLS.GT.7) THEN                                               E00500
         READ (IRD,'(7A10)') (XSNAME(I),I=1,7)                            E00510
         READ (IRD,'(8A10)') (XSNAME(I),I=8,IXMOLS)                       E00520
      ELSE                                                                E00530
         READ (IRD,'(7A10)') (XSNAME(I),I=1,IXMOLS)                       E00540
      ENDIF                                                               E00550
C                                                                         E00560
C     Left-justify all inputed names                                      E00570
C                                                                         E00580
      DO 15 I=1,IXMOLS                                                    E00582
         CALL CLJUST (XSNAME(I),10)                                       E00590
 15   CONTINUE
C                                                                         E00600
CPRT  WRITE(IPR,'(/,''  THE FOLLOWING MOLECULES ARE REQUESTED:'',//,      E00610
CPRT 1    (5X,I5,2X,A))') (I,XSNAME(I),I=1,IXMOLS)                        E00620
C                                                                         E00630
C     MATCH THE NAMES READ IN AGAINST THE NAMES STORED IN ALIAS           E00640
C     AND DETERMINE THE INDEX VALUE.  STOP IF NO MATCH IS FOUND.          E00650
C     NAME MUST BE ALL IN CAPS.                                           E00660
C                                                                         E00670
      DO 40 I = 1, IXMOLS                                                 E00680
         DO 20 J = 1, IXMAX                                               E00690
            IF ((XSNAME(I).EQ.ALIAS(1,J)) .OR.                            E00700
     *          (XSNAME(I).EQ.ALIAS(2,J)) .OR.                            E00710
     *          (XSNAME(I).EQ.ALIAS(3,J)) .OR.                            E00720
     *          (XSNAME(I).EQ.ALIAS(4,J))) THEN                           E00730
               IXINDX(I) = J                                              E00740
               GO TO 30                                                   E00750
            ENDIF                                                         E00760
   20    CONTINUE                                                         E00770
C                                                                         E00780
C         NO MATCH FOUND                                                  E00790
C                                                                         E00800
         WRITE (IPR,900) XSNAME(I)                                        E00810
         STOP 'STOPPED IN XSREAD'                                         E00820
C                                                                         E00830
   30    CONTINUE                                                         E00840
         IXFLG(I) = 0                                                     E00850
   40 CONTINUE                                                            E00860
C                                                                         E00870
C     READ IN "CROSS SECTION" MASTER FILE FSCDXS                          E00880
C                                                                         E00890
      IXFIL = 8                                                           E00900
      OPEN (IXFIL,FILE='FSCDXS',STATUS='OLD',FORM='FORMATTED')            E00910
      REWIND IXFIL                                                        E00920
      READ (IXFIL,905)                                                    E00930
C                                                                         E00940
   50 READ (IXFIL,910,END=80) XSREC                                       E00950
C                                                                         E00960
      IF (CFLG.EQ.CASTSK) GO TO 50                                        E00970
      IF (CFLG.EQ.CPRCNT) GO TO 80                                        E00980
C                                                                         E00990
      READ (XSREC,915) XNAME,V1X,V2X,DVX,NTEMP,IFRM,CFRM,                 E01000
     *                 (XFILS(I),I=1,NTEMP)                               E01010
C                                                                         E01020
C     LEFT-JUSTIFY INPUTED NAME                                           E01030
C                                                                         E01040
      CALL CLJUST (XNAME,10)                                              E01050
C                                                                         E01060
C     CHECK MASTER FILE FOR CROSS SECTION AND STORE DATA                  E01070
C                                                                         E01080
      NUMXS = IXMOLS                                                      E01090
      DO 70 I = 1, IXMOLS                                                 E01100
         IF ((XNAME.EQ.ALIAS(1,IXINDX(I))) .OR.                           E01110
     *       (XNAME.EQ.ALIAS(2,IXINDX(I))) .OR.                           E01120
     *       (XNAME.EQ.ALIAS(3,IXINDX(I))) .OR.                           E01130
     *       (XNAME.EQ.ALIAS(4,IXINDX(I)))) THEN                          E01140
            IXFLG(I) = 1                                                  E01150
            IF (V2X.GT.XV1.AND.V1X.LT.XV2) THEN                           E01160
               NSPECR(I) = NSPECR(I)+1                                    E01170
               IF (NSPECR(I).GT.6) THEN                                   E01180
                  WRITE (IPR,920) I,XSNAME(I),NSPECR(I)                   E01190
                  STOP ' XSREAD - NSPECR .GT. 6'                          E01200
               ENDIF                                                      E01210
               IXFORM(NSPECR(I),I) = 91                                   E01220
               IF (IFRM.EQ.86) IXFORM(NSPECR(I),I) = IFRM                 E01230
               IF (CFRM.NE.CN)                                            E01240
     *             IXFORM(NSPECR(I),I) = IXFORM(NSPECR(I),I)+100          E01250
               IF (CFRM.EQ.CF)                                            E01260
     *             IXFORM(NSPECR(I),I) = -IXFORM(NSPECR(I),I)             E01270
               NTEMPF(NSPECR(I),I) = NTEMP                                E01280
               V1FX(NSPECR(I),I) = V1X                                    E01290
               V2FX(NSPECR(I),I) = V2X                                    E01300
C                                                                         E01301
C     3.58115E-07 = SQRT( 2.*ALOG(2.)*AVOG*BOLTZ/(CLIGHT*CLIGHT) )        E01302
C                                                                         E01303
               XDOPLR(NSPECR(I),I)=3.58115E-07*(0.5*(V1X+V2X))*           E01304
     *                             SQRT(T296/XSMASS(IXINDX(I)))           E01305
C                                                                         E01306
               DO 60 J = 1, NTEMP                                         E01310
                  XSFILE(J,NSPECR(I),I) = XFILS(J)                        E01320
   60          CONTINUE                                                   E01330
            ENDIF                                                         E01340
         ENDIF                                                            E01350
   70 CONTINUE                                                            E01360
C                                                                         E01370
      GO TO 50                                                            E01380
C                                                                         E01390
   80 IXFLAG = 0                                                          E01400
      DO 90 I = 1, IXMOLS                                                 E01410
         IF (IXFLG(I).EQ.0) THEN                                          E01420
            WRITE (IPR,925) XSNAME(I)                                     E01430
            IXFLAG = 1                                                    E01440
         ENDIF                                                            E01450
   90 CONTINUE                                                            E01460
      IF (IXFLAG.EQ.1) STOP ' IXFLAG - XSREAD '                           E01470
C                                                                         E01480
      RETURN                                                              E01490
C                                                                         E01500
  900 FORMAT (/,'  THE NAME: ',A10, ' IS NOT ONE OF THE ',                E01510
     *        'CROSS SECTION MOLECULES. CHECK THE SPELLING.')             E01520
  905 FORMAT (/)                                                          E01530
  910 FORMAT (A120)                                                       E01540
  915 FORMAT (A10,2F10.4,F10.8,I5,5X,I5,A1,4X,6A10)                       E01550
  920 FORMAT (/,'******* ERROR IN XSREAD ** MOLECULE SECLECTED -',A10,    E01560
     *        '- HAS ',I2,' SPECTRAL REGIONS ON FILE FSCDXS, BUT THE',    E01570
     *        ' MAXIMUM ALLOWED IS 6 *******',/)                          E01580
  925 FORMAT (/,'******* MOLECULE SELECTED -',A10,'- IS NOT FOUND ON',    E01590
     *        ' FILE FSCDXS *******',/)                                   E01600
C                                                                         E01610
      END                                                                 E01620
      BLOCK DATA BXSECT                                                   E01630
C                                                                         E01640
      IMPLICIT REAL*8           (V)                                     ! E01650
C                                                                         E01660
C**   XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES          E01670
C**            (NOTE: ALL NAMES ARE LEFT-JUSTIFIED)                       E01680
C                                                                         E01690
      CHARACTER*10 XSFILE,XSNAME,ALIAS                                    E01700
      COMMON /XSECTI/ XSMAX(6,5,35),XSTEMP(6,5,35),NPTSFX(5,35),          E02850
     *                NFILEX(5,35),NLIMX                                  E02860
      COMMON /XSECTF/ XSFILE(6,5,35),XSNAME(35),ALIAS(4,35)               E01710
      COMMON /XSECTR/ V1FX(5,35),V2FX(5,35),DVFX(5,35),WXM(35),           E01720
     *                NTEMPF(5,35),NSPECR(35),IXFORM(5,35),               E01730
     *                XSMASS(35),XDOPLR(5,35),NUMXS,IXSBIN                E01740
      COMMON /XSECTS/ JINPUT,NMODES,NPANEL,NDUM,V1XS,V2XS,DVXS,NPTSXS     E02870
C                                                                         E01750
      DATA NMODES / 1 /,NPANEL / 0 /,V1XS / 0.0 /,V2XS / 0.0 /,           E02990
     *     DVXS / 0.0 /,NPTSXS / 0 /                                      E03000
      DATA XSMAX / 1050*0.0 /                                             E03010
      DATA (ALIAS(1,I),I=1,35)/                                           E01760
     *    'CLONO2    ', 'HNO4      ', 'CHCL2F    ', 'CCL4      ',         E01770
     *    'CCL3F     ', 'CCL2F2    ', 'C2CL2F4   ', 'C2CL3F3   ',         E01780
     *    'N2O5      ', 'HNO3      ', 'CF4       ', 'CHCLF2    ',         E01790
     *    'CCLF3     ', 'C2CLF5    ', 21*' ZZZZZZZZ ' /                   E01800
      DATA (ALIAS(2,I),I=1,35)/                                           E01810
     *    'CLNO3     ', ' ZZZZZZZZ ', 'CFC21     ', ' ZZZZZZZZ ',         E01820
     *    'CFCL3     ', 'CF2CL2    ', 'C2F4CL2   ', 'C2F3CL3   ',         E01830
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CHF2CL    ',         E01840
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 21*' ZZZZZZZZ ' /                   E01850
      DATA (ALIAS(3,I),I=1,35)/                                           E01860
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CFC21     ', ' ZZZZZZZZ ',         E01870
     *    'CFC11     ', 'CFC12     ', 'CFC114    ', 'CFC113    ',         E01880
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CFC14     ', 'CFC22     ',         E01890
     *    'CFC13     ', 'CFC115    ', 21*' ZZZZZZZZ ' /                   E01900
      DATA (ALIAS(4,I),I=1,35)/                                           E01910
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'F21       ', ' ZZZZZZZZ ',         E01920
     *    'F11       ', 'F12       ', 'F114      ', 'F113      ',         E01930
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'F14       ', 'F22       ',         E01940
     *    'F13       ', 'F115      ', 21*' ZZZZZZZZ ' /                   E01950
C                                                                         E01960
C     XSMASS IS MASS OF EACH CROSS-SECTION                                E01961
C                                                                         E01962
      DATA XSMASS/                                                        E01963
     1      97.46     ,   79.01     ,  102.92     ,  153.82     ,         E01964
     2     137.37     ,  120.91     ,  170.92     ,  187.38     ,         E01965
     3     108.01     ,   63.01     ,   88.00     ,   86.47     ,         E01966
     4     104.46     ,  154.47     ,  21*0.00 /                          E01967
C                                                                         E01969
      DATA V1FX / 175*0.0 /,V2FX / 175*0.0 /,DVFX / 175*0.0 /,            E01970
     *     WXM / 35*0.0 /                                                 E01980
      DATA NTEMPF / 175*0 /,NSPECR / 35*0 /,IXFORM / 175*0 /,             E01990
     *     NUMXS / 0 /                                                    E02000
C                                                                         E02010
      END                                                                 E02020
      SUBROUTINE CLJUST (CNAME,NCHAR)                                     E02030
C                                                                         E02040
C     THIS SUBROUTINE LEFT-JUSTIFIES THE CHARACTER CNAME                  E02050
C                                                                         E02060
      CHARACTER*(*) CNAME                                                 E02070
      CHARACTER*25 CTEMP                                                  E02070
      CHARACTER*1  CTEMP1(25),BLANK                                       E02080
      EQUIVALENCE (CTEMP,CTEMP1(1))                                       E02090
C                                                                         E02100
      DATA BLANK / ' '/                                                   E02110
C                                                                         E02120
         CTEMP = CNAME                                                    E02140
         JJ=0                                                             E02145
         DO 10 J = 1, NCHAR                                               E02150
            IF (CTEMP1(J).NE.BLANK) THEN                                  E02160
               JJ = J                                                     E02170
               IF (JJ.EQ.1) GO TO 50                                      E02180
               GO TO 20                                                   E02190
            ENDIF                                                         E02200
   10    CONTINUE                                                         E02210
         IF (JJ .EQ. 0) GO TO 50                                          E02215
C                                                                         E02220
   20    KCNT = 0                                                         E02230
         DO 30 K = JJ, NCHAR                                              E02240
            KCNT = KCNT+1                                                 E02250
            CTEMP1(KCNT) = CTEMP1(K)                                      E02260
   30    CONTINUE                                                         E02270
C                                                                         E02280
         KK = NCHAR-JJ+2                                                  E02290
         DO 40 L = KK,NCHAR                                               E02300
            CTEMP1(L) = BLANK                                             E02310
   40    CONTINUE                                                         E02320
         CNAME = CTEMP                                                    E02330
   50 CONTINUE                                                            E02340
C                                                                         E02350
      RETURN                                                              E02360
C                                                                         E02370
      END                                                                 E02380
      SUBROUTINE XSECTM (IFST,IR4)                                        E02390
C                                                                         E02400
      IMPLICIT REAL*8           (V)                                     ! E02410
C                                                                         E02420
C     THIS SUBROUTINE MOVES THE CROSS SECTIONS INTO                       E02430
C     THE APPROPRIATE ARRAY R1, R2, R3, R4, OR ABSRB                      E02440
C                                                                         E02450
C                       A.E.R. INC.    (AUGUST 1990)                      E02460
C                                                                         E02470
      COMMON VNU(250),SP(250),ALFA0(250),EPP(250),MOL(250),HWHMS(250),    E02480
     *       TMPALF(250),PSHIFT(250),IFLG(250),SPPSP(250),RECALF(250),    E02490
     *       ZETAI(250),IZETA(250)                                        E02500
      COMMON RR1(6099),RR2(2075),RR3(429)
      COMMON /IOU/ IOUT(250)                                              E02520
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(2030)                E02530
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   E02540
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       E02550
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       E02560
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      E02570
     *              EXTID(10)                                             E02580
C                                                                         E02590
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ
C                                                                         E02610
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       E02620
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   E02630
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    E02640
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           E02650
      COMMON /XSUB/ VBOT,VTOP,VFT,LIMIN,ILO,IHI,IEOF,IPANEL,ISTOP,IDATA   E02660
      COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)        E02670
      COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,            E02680
     *                HWF3,DXF3,NX3,N3MAX                                 E02690
      COMMON /SUB1/ MAX1,MAX2,MAX3,NLIM1,NLIM2,NLIM3,NLO,NHI,DVR2,DVR3,   E02700
     *              N1R1,N2R1,N1R2,N2R2,N1R3,N2R3                         E02710
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,         E02720
     *               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL                      E02730
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         E02740
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        E02750
     *              NLTEFL,LNFIL4,LNGTH4                                  E02760
      COMMON /XSHEAD/ HEADT1(35)                                          E02770
      COMMON /XSTMPR/ PF,TF,PDX(6,5,35),DVXPR(5,35),IXBIN(5,35),          E02780
     *                IXSBN(5,35)                                         E02790
      COMMON /XSECTP/ V1X,V2X,DVX,NPTSX,RX(13000)                         E02800
      COMMON /XSECTD/ V1DX,V2DX,DVDX,NPTSDX,RDX1(520),RDX2(520)           E02810
      COMMON /XSECTF/ XSFILE(6,5,35),XSNAME(35),ALIAS(4,35)               E02820
      COMMON /XSECTR/ V1FX(5,35),V2FX(5,35),DVFX(5,35),WXM(35),           E02830
     *                NTEMPF(5,35),NSPECR(35),IXFORM(5,35),               E02840
     *                XSMASS(35),XDOPLR(5,35),NUMXS,IXSBIN                E02845
      COMMON /XSECTI/ XSMAX(6,5,35),XSTEMP(6,5,35),NPTSFX(5,35),          E02850
     *                NFILEX(5,35),NLIMX                                  E02860
      COMMON /XSECTS/ JINPUT,NMODES,NPANEL,NDUM,V1XS,V2XS,DVXS,NPTSXS     E02870
      CHARACTER*10 XSFILE,XSNAME,ALIAS                                    E02880
      CHARACTER HEADT1*100                                                E02890
C                                                                         E02900
      DIMENSION R1(4050),R2(1050),R3(300)
      DIMENSION FILHDR(2)                                                 E02910
      LOGICAL OPCL                                                        E02920
C                                                                         E02930
      EQUIVALENCE (R1(1),RR1(2049)),(R2(1),RR2(1025)),(R3(1),RR3(129))
      EQUIVALENCE (IHIRAC,FSCDID(1)) , (ILBLF4,FSCDID(2)),                E02940
     *            (IXSCNT,FSCDID(3)) , (IAERSL,FSCDID(4)),                E02950
     *            (JRAD,FSCDID(9)) , (IATM,FSCDID(15)),                   E02960
     *            (XID(1),FILHDR(1))                                      E02970
C                                                                         E02980
      DATA IFILE,JFILE / 91,92 /                                          E03020
C                                                                         E03030
C**********************************************************************   E03040
C     NUMXS IS THE NUMBER OF 'CROSS SECTION' MOLECULES TO BE USED         E03050
C                                                                         E03060
C     XSFILE(ITEMP,ISPEC,NXS) IS THE NAME OF THE FILE CONTAINING THE      E03070
C                             'CROSS SECTION' DATA.  THE THREE INDICES    E03080
C                             ARE DEFINED AS FOLLOWS:                     E03090
C                                                                         E03100
C                             ITEMP - DENOTES THE TEMPERATURE FOR WHICH   E03110
C                                     THE 'CROSS SECTION' IS VALID        E03120
C                                     (IMPLEMENTED FOR HITRAN 91 TAPE)    E03130
C                             ISPEC - DENOTES THE SECTRAL REGION FOR      E03140
C                                     WHICH THE FILE PERTAINS             E03150
C                             NXS   - IS THE INCREMENT FOR THE 'CROSS     E03160
C                                     SECTION' INDEX                      E03170
C                                                                         E03180
C     NTEMPF(ISPEC,NXS) IS THE NUMBER OF TEMPERATURE FILES TO BE USED     E03190
C                       FOR EACH SPECTRAL REGION OF EACH MOLECULE         E03200
C                                                                         E03210
C     NSPECR(NXS) IS THE NUMBER OF SPECTRAL REGIONS FOR THE MOLECULE NX   E03220
C**********************************************************************   E03230
C                                                                         E03240
      NLIMX = 510                                                         E03250
      LIMOUT = 13000                                                      E03260
      IRPEAT = 0                                                          E03270
C                                                                         E03280
      PF = PAVE                                                           E03290
      TF = TAVE                                                           E03300
C                                                                         E03310
C     IF FIRST ENTRANCE, RESET QUANTITES                                  E03320
C                                                                         E03330
      IF (IFST.EQ.-99) THEN                                               E03340
         IFST = 0                                                         E03350
         JINPUT = -1                                                      E03360
         NMODES = 1                                                       E03370
C                                                                         E03380
         V1X = V1XS                                                       E03390
         V2X = V2XS                                                       E03400
         DVX = DVXS                                                       E03410
         NPTSX = NPTSXS                                                   E03420
C                                                                         E03430
         DO 10 NI = 1, NUMXS                                              E03440
            DO  9 NS = 1, NSPECR(NI)                                      E03450
               IXBIN(NS,NI) = 1                                           E03460
               IXSBN(NS,NI) = 0                                           E03470
               NFILEX(NS,NI) = ABS(NFILEX(NS,NI))                         E03480
 9          continue
 10      CONTINUE                                                         E03490
      ENDIF                                                               E03500
C                                                                         E03510
C     CHECK V1X FOR INPUT                                                 E03520
C                                                                         E03530
   20 VFX2 = VFT+2.*DVX+FLOAT(NHI)*DV                                     E03540
      IF (IR4.EQ.1) VFX2 = V2R4+2.*DVX                                    E03550
      IF (V1X.GT.VFX2) GO TO 140                                          E03560
C                                                                         E03570
      IF (JINPUT.EQ.-1) THEN                                              E03580
         JINPUT = 1                                                       E03590
      ELSE                                                                E03600
         VFX2 = MIN(VFX2,V2+2.*DVX)                                       E03610
         IF (VFX2.GT.V2X) THEN                                            E03620
            JINPUT = 1                                                    E03630
            IF (IRPEAT.EQ.1) THEN                                         E03640
               V1X = V2X-2.*DVX                                           E03650
            ELSE                                                          E03660
               V1X = VFT-2.*DVX                                           E03670
               IF (IR4.EQ.1) V1X = V1R4-2.*DVX                            E03680
            ENDIF                                                         E03690
            V2X = V1X+FLOAT(LIMOUT-1)*DVX                                 E03700
            IF (V2X.GT.V2) V2X = V1X+FLOAT(INT((V2-V1X)/DVX)+3)*DVX       E03710
            NPTSX = (V2X-V1X)/DVX+1                                       E03720
         ENDIF                                                            E03730
         IFL = 0                                                          E03740
         V1XT = V1X+2.*DVX                                                E03750
         V2XT = V2X-2.*DVX                                                E03760
         IF (V1XT.GT.V2XT) GO TO 140                                      E03770
         DO 30 NI = 1, NUMXS                                              E03780
            DO 29 NS = 1, NSPECR(NI)                                      E03790
               IF (NFILEX(NS,NI).EQ.0) GO TO 30                           E03800
               IF (V1FX(NS,NI).LE.V1XT.AND.V2FX(NS,NI).GE.V1XT) IFL = 1   E03810
               IF (V1FX(NS,NI).LE.V2XT.AND.V2FX(NS,NI).GE.V2XT) IFL = 1   E03820
               IF (V1FX(NS,NI).GE.V1XT.AND.V2FX(NS,NI).LE.V2XT) IFL = 1   E03830
 29         CONTINUE                                                         E03840
 30      CONTINUE                                                         E03840
         IF (IFL.EQ.0) GO TO 140                                          E03850
      ENDIF                                                               E03860
C                                                                         E03870
C     READ IN CROSS SECTION                                               E03880
C                                                                         E03890
      IF (JINPUT.EQ.1) THEN                                               E03900
         JINPUT = 0                                                       E03910
         DO 40 I = 1, LIMOUT                                              E03920
            RX(I) = 0.0                                                   E03930
 40      CONTINUE                                                         E03940
         DO 50 I = 1, NLIMX+10                                            E03950
            RDX1(I) = 0.0                                                 E03960
            RDX2(I) = 0.0                                                 E03970
 50      CONTINUE                                                         E03980
C                                                                         E03990
C     FOR NPANEL = 0, READ IN FILE HEADERS                                E04000
C                                                                         E04010
         IF (NPANEL.EQ.0) THEN                                            E04020
            DVXMIN = V2-V1                                                E04030
            V1XMIN = V2                                                   E04040
            IMAX = 0                                                      E04050
            NT2 = 0                                                       E04060
            NMODE = 0                                                     E04070
            DO 60 NI = 1, NUMXS                                           E04080
               DO 59 NS = 1, NSPECR(NI)                                   E04090
                  DO 58 NT1 = 1, NTEMPF(NS,NI)                            E04100
                     NFILEX(NS,NI) = 1                                    E04110
                     CALL CPUTIM (TIME0)                                  E04120
                     CALL XSECIN (NPANEL,NI,NS,NT1,NT2,NMODE,NSKIP,       E04130
     *                            IMAX,NEOF)                              E04140
                     CALL CPUTIM (TIME)                                   E04150
                     TXSRDF = TXSRDF+TIME-TIME0                           E04160
C                                                                         E04170
C     CHECK FOR WAVENUMBER BOUNDS AND SMALLEST DV                         E04180
C                                                                         E04190
                     IF (V1DX.GT.V2.OR.V2DX.LT.V1.OR.NEOF.EQ.1) THEN      E04200
                        NFILEX(NS,NI) = 0                                 E04210
                     ELSE                                                 E04220
                        DVXMIN = MIN(DVXMIN,DVDX)                         E04230
                        V1XMIN = MIN(V1XMIN,V1DX)                         E04240
                        V1FX(NS,NI) = V1DX                                E04250
                        V2FX(NS,NI) = V2DX                                E04260
                        DVFX(NS,NI) = DVDX                                E04270
                        NPTSFX(NS,NI) = NPTSDX                            E04280
                     ENDIF                                                E04290
C                                                                         E04291
Cmji  CHECK FOR TEMPERATURES; MUST BE IN ASCENDING ORDER                  E04292
C                                                                         E04293
                     IF (NT1.GT.1.AND.XSTEMP(NT1,NS,NI).LT.               E04294
     *                                XSTEMP(NT1-1,NS,NI)) THEN           E04295
                        WRITE(IPR,900)                                    E04296
                        STOP 'XSTEMP - XSECTM'                            E04297
                     ENDIF                                                E04298
 58               CONTINUE                                                
 59            CONTINUE                                                      E04300
 60         CONTINUE                                                      E04300
            DVX = DVXMIN                                                  E04310
            V1X = MAX(VFT,V1XMIN)                                         E04320
            V1X = V1X-2.*DVX                                              E04330
            V2X = V1X+FLOAT(LIMOUT-1)*DVX                                 E04340
            IF (V2X.GT.V2) V2X = V1X+FLOAT(INT((V2-V1X)/DVX)+2)*DVX       E04350
            NPTSX = (V2X-V1X)/DVX+1                                       E04360
            V1XS = V1X                                                    E04370
            V2XS = V2X                                                    E04380
            DVXS = DVX                                                    E04390
            NPTSXS = NPTSX                                                E04400
            IF (V1X.GT.VFX2) THEN                                         E04410
               JINPUT = 1                                                 E04420
               NPANEL = -1                                                E04430
               GO TO 140                                                  E04440
            ENDIF                                                         E04450
         ENDIF                                                            E04460
C                                                                         E04470
         NFILET = 0                                                       E04480
         NMODES = 0                                                       E04490
C                                                                         E04500
         DO 110 NI = 1, NUMXS                                             E04510
            DO 109 NS = 1, NSPECR(NI)                                     E04520
               NPANEL = -1                                                E04530
               IF (NFILEX(NS,NI).LE.0) GO TO 105                          E04540
               IF (V1FX(NS,NI).GT.V2X) GO TO 105                          E04550
               IF (V2FX(NS,NI).LT.V1X) THEN                               E04560
                  NFILEX(NS,NI) = -NFILEX(NS,NI)                          E04570
                  GO TO 105                                               E04580
               ENDIF                                                      E04590
C                                                                         E04600
C     DETERMINE TEMPERATURE FILES AND TEST ON DPTMIN                      E04610
C                                                                         E04620
               CALL XSNTMP (NI,NS,NT1,NT2,NMODE)                          E04630
C                                                                         E04640
C     DPTMIN TEST - IF NMODE = 0, SKIP CROSS SECTION                      E04650
C                                                                         E04660
               IF (NMODE.EQ.0) GO TO 105                                  E04670
               NMODES = NMODES+NMODE                                      E04680
C                                                                         E04690
C     FOR PRESSURE BROADENED CROSS-SECTION                                E04700
C     CREATE TEMPERATURE AVERAGED BINARY FILE                             E04710
C                                                                         E04720
               IF (IXSBIN.EQ.0.AND.IXBIN(NS,NI).EQ.1) THEN                E04730
                  CALL CPUTIM (TIME0)                                     E04740
                  CALL XSBINF (NI,NS,NT1,NT2,NMODE)                       E04750
                  CALL CPUTIM (TIME)                                      E04760
                  TXSCNV = TXSCNV+TIME-TIME0                              E04770
                  IXBIN(NS,NI) = 0                                        E04780
               ENDIF                                                      E04790
               IF (IXSBN(NS,NI).EQ.1) THEN                                E04800
                  DVFXX = DVXPR(NS,NI)                                    E04810
               ELSE                                                       E04820
                  DVFXX = DVFX(NS,NI)                                     E04830
               ENDIF                                                      E04840
               NFILET = NFILET+NFILEX(NS,NI)                              E04850
C                                                                         E04860
               NNSKIP = (V1X-V1FX(NS,NI))/DVFXX                           E04870
               NSKIP = (NNSKIP-3)/10                                      E04880
               NSKIP = MAX(NSKIP,0)                                       E04890
               NRSKIP = NSKIP*10                                          E04900
               NBSKIP = NSKIP                                             E04910
C                                                                         E04920
C     FOR BLOCKED DATA, V1FP MUST REFLECT SHORT RECORD                    E04930
C                                                                         E04940
               IAFORM = ABS(IXFORM(NS,NI))                                E04950
               IF (IXSBN(NS,NI).EQ.1) IAFORM = IAFORM+100                 E04960
               IF (IAFORM.GT.100) THEN                                    E04970
                  NBSKIP = NSKIP/51                                       E04980
                  NRSKIP = (NBSKIP-1)*510+500                             E04990
                  NRSKIP = MAX(NRSKIP,0)                                  E05000
               ENDIF                                                      E05010
               V1FP = V1FX(NS,NI)+FLOAT(NRSKIP)*DVFXX                     E05020
               V2FP = V2X+2.0*DVFXX                                       E05030
               V2FP = MIN(V2FP,V2FX(NS,NI))                               E05040
               NMAX = (V2FP-V1FP)/DVFXX+1.                                E05050
               NPAN = (NMAX+NLIMX-1)/NLIMX                                E05060
               IF (IAFORM.GT.100.AND.NPANEL.LE.0.AND.NBSKIP.EQ.0) THEN    E05070
                  NPTST = NMAX-500-(NPAN-1)*NLIMX                         E05080
                  IF (NPTST.GT.0) NPAN = NPAN+1                           E05090
                  IF (NMAX.GT.500) NMAX = NMAX+10                         E05100
               ENDIF                                                      E05110
               N2RX = ((V1FP-4.*DVFXX-V1X)/DVX+0.999)-1.                  E05120
               N2RX = MAX(N2RX,0)                                         E05130
C                                                                         E05140
C     IMAX = -4 TO PLACE THE FIRST PANEL V1 AT ARRAY LOCATION 1           E05150
C                                                                         E05160
               IMAX = -4                                                  E05170
               DO 100 NP = 1, NPAN                                        E05180
                  V1FP = V1FP+FLOAT(IMAX)*DVFXX                           E05190
                  IMAX = NMAX-(NP-1)*NLIMX                                E05200
                  IF (IAFORM.GT.100.AND.NPANEL.LE.0.AND.                  E05210
     *                NBSKIP.EQ.0.AND.IMAX.GT.500) IMAX = 500             E05220
                  IMAX = MIN(IMAX,NLIMX)                                  E05230
C                                                                         E05240
C     FOR V2FP IMAX + 3 GIVES US ARRAY LOCATION 514                       E05250
C             (504 FOR FIRST PANEL OF BLOCKED DATA)                       E05260
C                                                                         E05270
                  V2FP = V1FP+FLOAT(IMAX+3)*DVFXX                         E05280
C                                                                         E05290
                  IF (NP.GT.1) THEN                                       E05300
                     DO 70 JI = 1, 4                                      E05310
                        RDX1(JI) = RDX1(JI+IMAXSV)                        E05320
                        RDX2(JI) = RDX2(JI+IMAXSV)                        E05330
   70                CONTINUE                                             E05340
                  ENDIF                                                   E05350
                  IMAXSV = IMAX                                           E05360
C                                                                         E05370
                  CALL CPUTIM (TIME0)                                     E05380
                  CALL XSECIN (NPANEL,NI,NS,NT1,NT2,NMODE,NSKIP,IMAX,     E05390
     *               NEOF)                                                E05400
                  CALL CPUTIM (TIME)                                      E05410
                  TXSRDF = TXSRDF+TIME-TIME0                              E05420
C                                                                         E05430
                  IF (NP.EQ.1) THEN                                       E05440
                     DO 80 JI = 1, 4                                      E05450
                        RDX1(JI) = RDX1(5)                                E05460
                        RDX2(JI) = RDX2(5)                                E05470
   80                CONTINUE                                             E05480
                     NPANEL = ABS(NPANEL)                                 E05490
                     NSKIP = 0                                            E05500
                  ENDIF                                                   E05510
C                                                                         E05520
C     IF LAST PANEL OF FILE, FILL IN ADDITIONAL POINTS                    E05530
C     TO ENSURE ENDPOINT CAPTURE                                          E05540
C                                                                         E05550
                  IF (NP.EQ.NPAN) THEN                                    E05560
                     JMAX = IMAX+4                                        E05570
                     DO 90 JI = 1, 4                                      E05580
                        KI = JMAX+JI                                      E05590
                        RDX1(KI) = RDX1(JMAX)                             E05600
                        RDX2(KI) = RDX2(JMAX)                             E05610
   90                CONTINUE                                             E05620
C                                                                         E05630
                     V2FP = V2FP+3.*DVFXX                                 E05640
                  ENDIF                                                   E05650
C                                                                         E05660
                  N1RX = MAX(1,N2RX+1)                                    E05670
                  N2RX = (V2FP-DVFXX-V1X)/DVX+.999                        E05680
                  N2RX = MIN(N2RX,LIMOUT)                                 E05690
C                                                                         E05700
                  WXM1 = WXM(NI)                                          E05710
C                                                                         E05720
C     FOR TWO TEMPERATURES LINEARLY INTERPOLATE FACTOR                    E05730
C                                                                         E05740
                  CALL CPUTIM (TIME0)                                     E05750
                  IF (NMODE.EQ.2.AND.IXBIN(NS,NI).EQ.1) THEN              E05760
                     TFACT2 = (TAVE-XSTEMP(NT1,NS,NI))/                   E05770
     *                        (XSTEMP(NT2,NS,NI)-XSTEMP(NT1,NS,NI))       E05780
                     TFACT1 = 1.-TFACT2                                   E05790
                     WXM1 = WXM(NI)*TFACT1                                E05800
                     WXM2 = WXM(NI)*TFACT2                                E05810
                     CALL XINT (V1FP,V2FP,DVFXX,RDX2,WXM2,V1X,DVX,RX,     E05820
     *                  N1RX,N2RX)                                        E05830
                  ENDIF                                                   E05840
                  CALL XINT (V1FP,V2FP,DVFXX,RDX1,WXM1,V1X,DVX,RX,N1RX,   E05850
     *               N2RX)                                                E05860
                  CALL CPUTIM (TIME)                                      E05870
                  TXSPNL = TXSPNL+TIME-TIME0                              E05880
C                                                                         E05890
 100           CONTINUE                                                   E05900
C
C              Continue for GOTO statements at E04540, E04550, E04580, &
C              E04670.
C
 105           CONTINUE
C                                                                         E05910
 109        CONTINUE                                                         E05920
 110     CONTINUE                                                         E05920
         IF (NFILET.EQ.0) GO TO 140                                       E05930
C                                                                         E05940
C     FACTOR OUT RADIATION FIELD IF REQUIRED                              E05950
C                                                                         E05960
         IF (JRAD.EQ.0) THEN                                              E05970
            CALL CPUTIM (TIME0)                                           E05980
            XKT = TAVE/RADCN2                                             E05990
            VI = V1X-DVX                                                  E06000
            VITST = VI                                                    E06010
            RDLAST = -1.                                                  E06020
            NPTSI1 = 0                                                    E06030
            NPTSI2 = 0                                                    E06040
C                                                                         E06050
  120       NPTSI1 = NPTSI2+1                                             E06060
            NPTSX = (V2X-V1X)/DVX+1                                       E06062
C                                                                         E06070
            VI = V1X+FLOAT(NPTSI1-1)*DVX                                  E06080
            RADVI = RADFNI(VI,DVX,XKT,VITST,RDEL,RDLAST)                  E06090
C                                                                         E06110
            NPTSI2 = (VITST-V1X)/DVX+1.001                                E06120
            NPTSI2 = MIN(NPTSI2,NPTSX)                                    E06130
C                                                                         E06140
            DO 130 I = NPTSI1, NPTSI2                                     E06150
               VI = VI+DVX                                                E06160
               RX(I) = RX(I)/RADVI                                        E06170
               RADVI = RADVI+RDEL                                         E06180
  130       CONTINUE                                                      E06190
C                                                                         E06200
            IF (NPTSI2.LT.NPTSX) GO TO 120                                E06210
            CALL CPUTIM (TIME)                                            E06220
            TXSPNL = TXSPNL+TIME-TIME0                                    E06230
C                                                                         E06240
         ENDIF                                                            E06250
      ENDIF                                                               E06260
      IF (NMODES.EQ.0) GO TO 140                                          E06270
C                                                                         E06280
C     DETERMINE TARGET ARRAY                                              E06290
C                                                                         E06300
C      ===> R1                                                            E06310
C                                                                         E06320
      CALL CPUTIM (TIME0)                                                 E06330
      IF (DVX.LT.DVR2) THEN                                               E06340
         CALL XINT (V1X,V2X,DVX,RX,1.0,VFT,DV,R1,N1R1,N2R1)               E06350
         IR4 = 0                                                          E06360
C                                                                         E06370
C      ===> R2                                                            E06380
C                                                                         E06390
      ELSEIF (DVX.LT.DVR3) THEN                                           E06400
         CALL XINT (V1X,V2X,DVX,RX,1.0,VFT,DVR2,R2,N1R2,N2R2)             E06410
         IR4 = 0                                                          E06420
C                                                                         E06430
C      ===> R3                                                            E06440
C                                                                         E06450
      ELSEIF (DVX.LT.DVR4.OR.ILBLF4.EQ.0) THEN                            E06460
         CALL XINT (V1X,V2X,DVX,RX,1.0,VFT,DVR3,R3,N1R3,N2R3)             E06470
         IR4 = 0                                                          E06480
C                                                                         E06490
C      ===> R4                                                            E06500
C                                                                         E06510
      ELSE                                                                E06520
         CALL XINT (V1X,V2X,DVX,RX,1.0,V1R4,DVR4,R4,1,NPTR4)              E06530
         IF (IR4.EQ.0) VFX2 = V2R4+2.*DVX                                 E06540
         IR4 = 1                                                          E06550
      ENDIF                                                               E06560
      CALL CPUTIM (TIME)                                                  E06570
      TXSPNL = TXSPNL+TIME-TIME0                                          E06580
C                                                                         E06590
      IRPEAT = 1                                                          E06600
      IF (VFX2.GT.V2X) GO TO 20                                           E06610
C                                                                         E06620
  140 INQUIRE (UNIT=IFILE,OPENED=OPCL)                                    E06630
      IF (OPCL) CLOSE (IFILE)                                             E06640
      INQUIRE (UNIT=JFILE,OPENED=OPCL)                                    E06650
      IF (OPCL) CLOSE (JFILE)                                             E06660
      IF (NMODES.EQ.0) IR4 = 1                                            E06670
C                                                                         E06680
 900  FORMAT(/,'******* ERROR IN XSECTM *******',/                        E06682
     *         'CROSS-SECTION FILES MUST BE IN ASCENDING ORDER ',         E06683
     *         'BY TEMPERATURE IN FSCDXS.')                               E06684
      RETURN                                                              E06690
C                                                                         E06700
      END                                                                 E06710
      SUBROUTINE XSECIN (NPANEL,NI,NS,NT1,NT2,NMODE,NSKIP,NMAX,IEOF)      E06720
C                                                                         E06730
      IMPLICIT REAL*8           (V)                                     ! E06740
C                                                                         E06750
C     THIS SUBROUTINE READS IN THE DESIRED CROSS SECTIONS                 E06760
C                                                                         E06770
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         E06780
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        E06790
     *              NLTEFL,LNFIL4,LNGTH4                                  E06800
C                                                                         E06810
      CHARACTER*8      XI1,       HMOLI1,      YID1   
      Real*8               SECAN1,       XALT1
C                                                                         E06830
      COMMON /FXSHDR/ XI1(10),SECAN1,PAV1,TAV1,HMOLI1(60),XALT1(4),       E06840
     *                W1(60),P1L,P1U,T1L,T1U,WBROA1,DVB,V1B,V2B,TBOUN1,   E06850
     *                EMISI1,FSCDI1(17),NMO1,LAYER1,Y11,YID1(10),LSTWD1   E06860
      COMMON /PXSHDR/ V1PX,V2PX,DVPX,NLIMPX,RBX(2050)                     E06870
      COMMON /XSECTP/ V1X,V2X,DVX,NPTSX,RX(13000)                         E06880
      COMMON /XSECTD/ V1DX,V2DX,DVDX,NPTSDX,RDX1(520),RDX2(520)           E06890
      COMMON /XSECTF/ XSFILE(6,5,35),XSNAME(35),ALIAS(4,35)               E06900
      COMMON /XSECTR/ V1FX(5,35),V2FX(5,35),DVFX(5,35),WXM(35),           E06910
     *                NTEMPF(5,35),NSPECR(35),IXFORM(5,35),               E06920
     *                XSMASS(35),XDOPLR(5,35),NUMXS,IXSBIN                E06925
      COMMON /XSECTI/ XSMAX(6,5,35),XSTEMP(6,5,35),NPTSFX(5,35),          E06930
     *                NFILEX(5,35),NLIMX                                  E06940
      COMMON /XSHEAD/ HEADT1(35)                                          E06950
      COMMON /XSTMPR/ PF,TF,PDX(6,5,35),DVXPR(5,35),IXBIN(5,35),          E06960
     *                IXSBN(5,35)                                         E06970
      COMMON /FLFORM/ CFORM                                               E06980
C                                                                         E06990
      CHARACTER*10 XSFILE,XSNAME,ALIAS,SOURCE(3),CTORR                    E07000
      CHARACTER AMOL*8,BMOL*6,HEADER*100,HEADT1*100                       E07010
      CHARACTER XSFIL1*10,XSFIL2*10,XSTMP*4,XSNUM*3,CI*1                  E07020
      CHARACTER CFORM*11,BFRM*10,UNBFRM*10,BLKFRM*10,BFORM*9              E07030
      LOGICAL OP,OPCL                                                     E07040
      DIMENSION RDXX1(516),RDXX2(516),RDXA1(510),RDXA2(510),RDXH1(500),   E07050
     *          RDXH2(500),FILHDR(2),PNLHDR(2),DUM(2)                     E07060
C                                                                         E07070
      EQUIVALENCE (RDX1(5),RDXX1(1),RDXA1(1),RDXH1(1)),                   E07080
     *            (RDX2(5),RDXX2(1),RDXA2(1),RDXH2(1)),                   E07090
     *            (XI1(1),FILHDR(1)) , (V1PX,PNLHDR(1))                   E07100
C                                                                         E07110
      DATA XSTMP / 'TMPX'/,LIMXX / 516 /,BFORM / 'FORMATTED'/             E07120
      DATA IFILE,JFILE / 91,92 /                                          E07130
      DATA UNBFRM / '(10E10.3)'/,BLKFRM / '(510E10.3)'/                   E07140
      DATA CTORR / '      TORR'/
C                                                                         E07150
C     DEFINE PRESSURE CONVERSIONS                                         E07160
C                                                                         E07170
C        PTORMB = 1013. MB / 760. TORR  (TORR TO MILLIBARS)               E07180
C        PATMMB = 1013. MB / 1.0  ATM   (ATMOPHERES TO MILLIBARS)         E07190
C                                                                         E07200
      PTORMB = 1013./760.                                                 E07210
      PATMMB = 1013.                                                      E07220
C                                                                         E07230
      IEOF = 0                                                            E07240
      ISFORM = IXFORM(NS,NI)                                              E07250
      NXMODE = NMODE                                                      E07260
      IF (IXSBN(NS,NI).EQ.1) THEN                                         E07270
         IF (ABS(ISFORM).LT.100) ISFORM = ABS(ISFORM)+100                 E07280
         ISFORM = -ISFORM                                                 E07290
         NXMODE = 1                                                       E07300
      ENDIF                                                               E07310
      IAFORM = ABS(ISFORM)                                                E07320
      IMFORM = MOD(IAFORM,100)                                            E07330
C                                                                         E07340
C     IF NPANEL <= 0, OPEN FILE AND READ HEADER                           E07350
C                                                                         E07360
      IF (NPANEL.LE.0) THEN                                               E07370
         IF (IXSBN(NS,NI).EQ.0) THEN                                      E07380
            XSFIL1 = XSFILE(NT1,NS,NI)                                    E07390
         ELSE                                                             E07400
            WRITE (XSNUM,'(I1,I2.2)') NS,NI                               E07410
            XSFIL1 = XSTMP//XSNUM                                         E07420
         ENDIF                                                            E07430
C                                                                         E07440
         INQUIRE (FILE=XSFIL1,OPENED=OP)                                  E07450
         INQUIRE (UNIT=IFILE,OPENED=OPCL)                                 E07460
         IF (.NOT.OP.AND.OPCL) CLOSE (IFILE)                              E07470
         IF (.NOT.OP) THEN                                                E07480
            IF (ISFORM.GT.0) THEN                                         E07490
               OPEN (IFILE,FILE=XSFIL1,STATUS='OLD',FORM=BFORM)           E07500
            ELSE                                                          E07510
               OPEN (IFILE,FILE=XSFIL1,STATUS='OLD',FORM=CFORM)           E07520
            ENDIF                                                         E07530
         ENDIF                                                            E07540
         REWIND IFILE                                                     E07550
         IF (NXMODE.EQ.2) THEN                                            E07560
            XSFIL2 = XSFILE(NT2,NS,NI)                                    E07570
C                                                                         E07580
            INQUIRE (FILE=XSFIL2,OPENED=OP)                               E07590
            INQUIRE (UNIT=JFILE,OPENED=OPCL)                              E07600
            IF (.NOT.OP.AND.OPCL) CLOSE (JFILE)                           E07610
            IF (.NOT.OP) THEN                                             E07620
               IF (ISFORM.GT.0) THEN                                      E07630
                  OPEN (JFILE,FILE=XSFIL2,STATUS='OLD',FORM=BFORM)        E07640
               ELSE                                                       E07650
                  OPEN (JFILE,FILE=XSFIL2,STATUS='OLD',FORM=CFORM)        E07660
               ENDIF                                                      E07670
            ENDIF                                                         E07680
            REWIND JFILE                                                  E07690
         ENDIF                                                            E07700
C                                                                         E07710
C     HEADER: 86 FORMAT                                                   E07720
C                                                                         E07730
C             AMOL,V1,V2,NPTS,BMOL,PRES,ICM,ITEMP,SOURCE                  E07740
C                                                                         E07750
C     HEADER: 91 FORMAT                                                   E07760
C                                                                         E07770
C             AMOL,V1,V2,NPTS,TEMP,PRES,SMAX,SOURCE                       E07780
C                                                                         E07790
C                                                                         E07800
C     IAFORM < 100, UNBLOCKED DATA (100 CHARACTERS/RECORD)                E07810
C                                                                         E07820
         IF (IAFORM.LT.100) THEN                                          E07830
            READ (IFILE,900,END=30) HEADER                                E07840
            HEADT1(NI) = HEADER                                           E07850
            IF (IMFORM.EQ.86) THEN                                        E07860
               READ (HEADER,905) AMOL,V1DX,V2DX,NPTSDX,BMOL,PRES,ICM,     E07870
     *                           ITEMP,SOURCE                             E07880
               IF (NPANEL.EQ.0) THEN                                      E07890
                  XSTEMP(NT1,NS,NI) = FLOAT(ITEMP)+273.15                 E07900
                  XSMAX(NT1,NS,NI) = 0.0                                  E07910
                  PDX(NT1,NS,NI) = PRES*PTORMB                            E07920
               ENDIF                                                      E07930
            ELSE                                                          E07940
               READ (HEADER,910) AMOL,V1DX,V2DX,NPTSDX,TEMP,PRES,SMAX,    E07950
     *                           SOURCE                                   E07960
               IF (NPANEL.EQ.0) THEN                                      E07970
                  XSTEMP(NT1,NS,NI) = TEMP                                E07980
                  XSMAX(NT1,NS,NI) = SMAX                                 E07990
                  IF (SOURCE(3).EQ.CTORR) THEN
                     PDX(NT1,NS,NI) = PRES*PTORMB                         E08000
                  ELSE
                     PDX(NT1,NS,NI) = PRES
                  ENDIF
               ENDIF                                                      E08010
            ENDIF                                                         E08020
            IF (NXMODE.EQ.2) THEN                                         E08030
               READ (JFILE,900,END=30) HEADER                             E08040
               IF (IMFORM.EQ.86) THEN                                     E08050
                  READ (HEADER,905) AMOL,V1DX,V2DX,NPTSDX,BMOL,PRES,      E08060
     *                              ICM,ITEMP,SOURCE                      E08070
                  IF (NPANEL.EQ.0) THEN                                   E08080
                     XSTEMP(NT2,NS,NI) = FLOAT(ITEMP)+273.15              E08090
                     XSMAX(NT2,NS,NI) = 0.0                               E08100
                     PDX(NT2,NS,NI) = PRES*PTORMB                         E08110
                  ENDIF                                                   E08120
               ELSE                                                       E08130
                  READ (HEADER,910) AMOL,V1DX,V2DX,NPTSDX,TEMP,PRES,      E08140
     *                              SMAX,SOURCE                           E08150
                  IF (NPANEL.EQ.0) THEN                                   E08160
                     XSTEMP(NT2,NS,NI) = TEMP                             E08170
                     XSMAX(NT2,NS,NI) = SMAX                              E08180
                     IF (SOURCE(3).EQ.CTORR) THEN
                        PDX(NT2,NS,NI) = PRES*PTORMB                      E08190
                     ELSE
                        PDX(NT2,NS,NI) = PRES
                     ENDIF
                  ENDIF                                                   E08200
               ENDIF                                                      E08210
            ENDIF                                                         E08220
         ELSE                                                             E08230
C                                                                         E08240
C     IAFORM > 100, BLOCKED DATA (51*100 CHARACTERS/RECORD)               E08250
C                                                                         E08260
            IF (ISFORM.GT.0) THEN                                         E08270
               READ (IFILE,915,END=30) HEADER,(RDXX1(J),J=1,500)          E08280
            ELSE                                                          E08290
               CALL BUFIN (IFILE,IEOF,FILHDR(1),NFHDRF)                   E08300
               IF (IEOF.LE.0) GO TO 30                                    E08310
               WRITE (HEADER,'(10A8)') XI1                                E08320
               CALL BUFIN (IFILE,IEOF,PNLHDR(1),NPHDRF)                   E08330
               IF (IEOF.LE.0) GO TO 30                                    E08340
               CALL BUFIN (IFILE,IEOF,RDXH1(1),NLIMPX)                    E08350
            ENDIF                                                         E08360
            HEADT1(NI) = HEADER                                           E08370
            IF (IMFORM.EQ.86) THEN                                        E08380
               READ (HEADER,905) AMOL,V1DX,V2DX,NPTSDX,BMOL,PRES,ICM,     E08390
     *                           ITEMP,SOURCE                             E08400
               IF (NPANEL.EQ.0) THEN                                      E08410
                  XSTEMP(NT1,NS,NI) = FLOAT(ITEMP)+273.15                 E08420
                  XSMAX(NT1,NS,NI) = 0.0                                  E08430
                  PDX(NT1,NS,NI) = PRES*PTORMB                            E08440
               ENDIF                                                      E08450
            ELSE                                                          E08460
               READ (HEADER,910) AMOL,V1DX,V2DX,NPTSDX,TEMP,PRES,SMAX,    E08470
     *                           SOURCE                                   E08480
               IF (NPANEL.EQ.0) THEN                                      E08490
                  XSTEMP(NT1,NS,NI) = TEMP                                E08500
                  XSMAX(NT1,NS,NI) = SMAX                                 E08510
                  IF (SOURCE(3).EQ.CTORR) THEN
                     PDX(NT1,NS,NI) = PRES*PTORMB                         E08520
                  ELSE
                     PDX(NT1,NS,NI) = PRES
                  ENDIF
               ENDIF                                                      E08530
            ENDIF                                                         E08540
            IF (NXMODE.EQ.2) THEN                                         E08550
               IF (ISFORM.GT.0) THEN                                      E08560
                  READ (JFILE,915,END=30) HEADER,(RDXX2(J),J=1,500)       E08570
               ELSE                                                       E08580
                  CALL BUFIN (JFILE,JEOF,FILHDR(1),NFHDRF)                E08590
                  IF (JEOF.LE.0) GO TO 30                                 E08600
                  WRITE (HEADER,'(10A8)') XI1                             E08610
                  CALL BUFIN (JFILE,JEOF,PNLHDR(1),NPHDRF)                E08620
                  IF (JEOF.LE.0) GO TO 30                                 E08630
                  CALL BUFIN (JFILE,JEOF,RDXH2(1),NLIMPX)                 E08640
               ENDIF                                                      E08650
               IF (IMFORM.EQ.86) THEN                                     E08660
                  READ (HEADER,905) AMOL,V1DX,V2DX,NPTSDX,BMOL,PRES,      E08670
     *                              ICM,ITEMP,SOURCE                      E08680
                  IF (NPANEL.EQ.0) THEN                                   E08690
                     XSTEMP(NT2,NS,NI) = FLOAT(ITEMP)+273.15              E08700
                     XSMAX(NT2,NS,NI) = 0.0                               E08710
                     PDX(NT2,NS,NI) = PRES*PTORMB                         E08720
                  ENDIF                                                   E08730
               ELSE                                                       E08740
                  READ (HEADER,910) AMOL,V1DX,V2DX,NPTSDX,TEMP,PRES,      E08750
     *                              SMAX,SOURCE                           E08760
                  IF (NPANEL.EQ.0) THEN                                   E08770
                     XSTEMP(NT2,NS,NI) = TEMP                             E08780
                     XSMAX(NT2,NS,NI) = SMAX                              E08790
                     IF (SOURCE(3).EQ.CTORR) THEN
                        PDX(NT2,NS,NI) = PRES*PTORMB                      E08800
                     ELSE
                        PDX(NT2,NS,NI) = PRES
                     ENDIF
                  ENDIF                                                   E08810
               ENDIF                                                      E08820
            ENDIF                                                         E08830
         ENDIF                                                            E08840
         DVDX = (V2DX-V1DX)/FLOAT(NPTSDX-1)                               E08850
C                                                                         E08860
C     FOR NPANEL = -1, SKIP REQUIRED NUMBER OF RECORDS                    E08870
C                                                                         E08880
         IF (NPANEL.EQ.-1) THEN                                           E08890
            NSTRT = 1                                                     E08900
            IF (IAFORM.GT.100) THEN                                       E08910
               NSKIP = NSKIP/51                                           E08920
               NSTRT = 2                                                  E08930
               IF (NSKIP.EQ.0) RETURN                                     E08940
            ENDIF                                                         E08950
C                                                                         E08960
            DO 10 I = NSTRT, NSKIP                                        E08970
               IF (ISFORM.GT.0) THEN                                      E08980
                  READ (IFILE,920,END=30) CI                              E08990
                  IF (NXMODE.EQ.2) READ (JFILE,920,END=30) CI             E09000
               ELSE                                                       E09010
                  CALL BUFIN (IFILE,IEOF,PNLHDR(1),NPHDRF)                E09020
                  IF (IEOF.LE.0) GO TO 30                                 E09030
                  CALL BUFIN (IFILE,IEOF,DUM(1),1)                        E09040
                  IF (NXMODE.EQ.2) THEN                                   E09050
                     CALL BUFIN (JFILE,JEOF,PNLHDR(1),NPHDRF)             E09060
                     IF (JEOF.LE.0) GO TO 30                              E09070
                     CALL BUFIN (JFILE,JEOF,DUM(1),1)                     E09080
                  ENDIF                                                   E09090
               ENDIF                                                      E09100
   10       CONTINUE                                                      E09110
         ENDIF                                                            E09120
      ENDIF                                                               E09130
C                                                                         E09140
C     FOR ABS(NPANEL) > 0, READ IN MORE DATA                              E09150
C                                                                         E09160
      IF (ABS(NPANEL).GT.0) THEN                                          E09170
         BFRM = UNBFRM                                                    E09180
         IF (IAFORM.GT.100) BFRM = BLKFRM                                 E09190
         IF (ISFORM.GT.0) THEN                                            E09200
            READ (IFILE,BFRM,END=30) (RDXX1(J),J=1,NMAX)                  E09210
            IF (NXMODE.EQ.2)                                              E09220
     *         READ (JFILE,BFRM,END=30) (RDXX2(J),J=1,NMAX)               E09230
         ELSE                                                             E09240
            CALL BUFIN (IFILE,IEOF,PNLHDR(1),NPHDRF)                      E09250
            IF (IEOF.LE.0) GO TO 30                                       E09260
            CALL BUFIN (IFILE,IEOF,RDXA1(1),NLIMPX)                       E09270
            IF (NXMODE.EQ.2) THEN                                         E09280
               CALL BUFIN (JFILE,JEOF,PNLHDR(1),NPHDRF)                   E09290
               IF (JEOF.LE.0) GO TO 30                                    E09300
               CALL BUFIN (JFILE,JEOF,RDXA2(1),NLIMPX)                    E09310
            ENDIF                                                         E09320
         ENDIF                                                            E09330
      ENDIF                                                               E09340
C                                                                         E09350
      DO 20 I = NMAX+1, LIMXX                                             E09360
         RDXX1(I) = 0.0                                                   E09370
         IF (NXMODE.EQ.2) RDXX2(I) = 0.0                                  E09380
   20 CONTINUE                                                            E09390
C                                                                         E09400
      RETURN                                                              E09410
C                                                                         E09420
   30 IEOF = 1                                                            E09430
C                                                                         E09440
      DO 40 I = NMAX+1, LIMXX                                             E09450
         RDXX1(I) = 0.0                                                   E09460
         IF (NXMODE.EQ.2) RDXX2(I) = 0.0                                  E09470
   40 CONTINUE                                                            E09480
C                                                                         E09490
      RETURN                                                              E09500
C                                                                         E09510
  900 FORMAT (A100)                                                       E09520
  905 FORMAT (A8,2F10.4,I10,1X,A6,F4.2,5X,I4,3X,I5,3A10)                  E09530
  910 FORMAT (A10,2F10.4,I10,3G10.3,3A10)                                 E09540
  915 FORMAT (A100,50(10E10.3))                                           E09550
  920 FORMAT (A1)                                                         E09560
C                                                                         E09570
      END                                                                 E09580
      SUBROUTINE XSNTMP (NI,NS,NT1,NT2,NMODE)                             E09590
C                                                                         E09600
      IMPLICIT REAL*8           (V)                                     ! E09610
C                                                                         E09620
C     THIS SUBROUTINE DETERMINES THE CORRECT MODE                         E09630
C     AND BRACKETS THE LAYER TEMPERATURE                                  E09640
C                                                                         E09650
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ
C                                                                         E09670
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       E09680
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   E09690
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    E09700
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,   E09710
     *              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,       E09720
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,       E09730
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,      E09740
     *              EXTID(10)                                             E09750
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           E09760
      COMMON /XSECTP/ V1X,V2X,DVX,NPTSX,RX(13000)                         E09770
      COMMON /XSECTD/ V1DX,V2DX,DVDX,NPTSDX,RDX1(520),RDX2(520)           E09780
      COMMON /XSECTF/ XSFILE(6,5,35),XSNAME(35),ALIAS(4,35)               E09790
      COMMON /XSECTR/ V1FX(5,35),V2FX(5,35),DVFX(5,35),WXM(35),           E09800
     *                NTEMPF(5,35),NSPECR(35),IXFORM(5,35),               E09810
     *                XSMASS(35),XDOPLR(5,35),NUMXS,IXSBIN                E09815
      COMMON /XSECTI/ XSMAX(6,5,35),XSTEMP(6,5,35),NPTSFX(5,35),          E09820
     *                NFILEX(5,35),NLIMX                                  E09830
      CHARACTER*10 XSFILE,XSNAME,ALIAS                                    E09840
C                                                                         E09850
      EQUIVALENCE (JRAD,FSCDID(9))                                        E09860
C                                                                         E09870
      NT1 = 0                                                             E09880
      NT2 = 0                                                             E09890
C                                                                         E09900
      IF (NTEMPF(NS,NI).LE.1) THEN                                        E09910
         NMODE = 1                                                        E09920
         NT1 = 1                                                          E09930
         NT2 = 1                                                          E09935
      ELSE                                                                E09940
         NMODE = 2                                                        E09950
         DO 10 I = 2, NTEMPF(NS,NI)                                       E09960
            IF (TAVE.LT.XSTEMP(I,NS,NI)) THEN                             E09970
               NT1 = I-1                                                  E09980
               NT2 = I                                                    E09990
               GO TO 20                                                   E10000
            ENDIF                                                         E10010
   10    CONTINUE                                                         E10020
      ENDIF                                                               E10030
C                                                                         E10040
   20 IF (NT1.EQ.0) THEN                                                  E10050
         NT2 = NTEMPF(NS,NI)                                              E10060
         NT1 = NT2-1                                                      E10070
      ENDIF                                                               E10080
C                                                                         E10090
C     CHECK VERSUS DPTMIN                                                 E10100
C                                                                         E10110
      IF (XSMAX(NT1,NS,NI).NE.0.0) THEN                                   E10120
         WXM1 = WXM(NI)*XSMAX(NT1,NS,NI)                                  E10130
         IF (WXM1.LT.DPTMIN) IDPTMN = IDPTMN+1                            E10132
         WXM2 = 0.0                                                       E10140
         IF (NMODE.EQ.2) THEN                                             E10142
            WXM2 = WXM(NI)*XSMAX(NT2,NS,NI)                               E10150
            IF (WXM2.LT.DPTMIN) IDPTMN = IDPTMN+1                         E10152
         ENDIF                                                            E10154
         IF (JRAD.EQ.0) THEN                                              E10160
            XKT1 = XSTEMP(NT1,NS,NI)/RADCN2                               E10170
            IF (NMODE.EQ.2) XKT2 = XSTEMP(NT2,NS,NI)/RADCN2               E10180
            VI = V1FX(NS,NI)                                              E10190
C                                                                         E10200
            RADVI1 = RADFN(VI,XKT1)                                       E10210
            IF (NMODE.EQ.2) RADVI2 = RADFN(VI,XKT2)                       E10220
            WXM1 = WXM(NI)*XSMAX(NT1,NS,NI)/RADVI1                        E10230
            IF (NMODE.EQ.2) WXM2 = WXM(NI)*XSMAX(NT2,NS,NI)/RADVI2        E10240
         ENDIF                                                            E10250
C                                                                         E10260
C     DETERMINE IDPTMN --- IF IDPTMN = NMODE  ==> BELOW THRESHOLD         E10270
C                                                 SKIP CROSS SECTION      E10280
C                                                                         E10290
         IDPTMN = 0                                                       E10300
         IF (IDPTMN.EQ.NMODE) NMODE = 0                                   E10330
      ENDIF                                                               E10340
C                                                                         E10350
      RETURN                                                              E10360
C                                                                         E10370
      END                                                                 E10380
      SUBROUTINE XSBINF (NI,NS,NT1,NT2,NMODE)                             E10390
C                                                                         E10400
      IMPLICIT REAL*8           (V)                                     ! E10410
C                                                                         E10420
C     THIS SUBROUTINE PERFORMS A TEMPERATURE DEPENDENT CONVOLUTION        E10430
C     ON THE CROSS-SECTIONS PRODUCING A BINARY INTERMEDIATE FILE          E10440
C                                                                         E10450
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ
C                                                                         E10470
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       E10480
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV0,V10,V20,TBOUND,   E10490
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    E10500
C                                                                         E10510
      CHARACTER*8      XI1,       HMOLI1,      YID1   
      Real*8               SECAN1,       XALT1
C                                                                         E10530
      COMMON /FXSHDR/ XI1(10),SECAN1,PAV1,TAV1,HMOLI1(60),XALT1(4),       E10540
     *                W1(60),P1L,P1U,T1L,T1U,WBROA1,DVB,V1B,V2B,TBOUN1,   E10550
     *                EMISI1,FSCDI1(17),NMO1,LAYER1,Y11,YID1(10),LSTWD1   E10560
      COMMON /PXSHDR/ V1PX,V2PX,DVPX,NLIMPX,RBX(2050)                     E10570
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         E10580
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        E10590
     *              NLTEFL,LNFIL4,LNGTH4                                  E10600
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   E10610
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       E10620
     *               DVSC,XDUM,V1SHFT                                     E10630
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,         E10640
     *               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL                      E10650
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 E10660
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
      COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,               E10670
     *                HWF3,DXF3,NX3,N3MAX                                 E10680
      COMMON /XSCINF/ HWHM,JEMIT,JFN,SAMPLE,SCANID,NPTS,XF(851)           E10690
C                                                                         E10700
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2           E10710
      COMMON /XSECTP/ V1X,V2X,DVX,NPTSX,RX(13000)                         E10720
      COMMON /XSECTD/ V1DX,V2DX,DVDX,NPTSDX,RDX1(520),RDX2(520)           E10730
      COMMON /XSECTF/ XSFILE(6,5,35),XSNAME(35),ALIAS(4,35)               E10740
      COMMON /XSECTR/ V1FX(5,35),V2FX(5,35),DVFX(5,35),WXM(35),           E10750
     *                NTEMPF(5,35),NSPECR(35),IXFORM(5,35),               E10760
     *                XSMASS(35),XDOPLR(5,35),NUMXS,IXSBIN                E10765
      COMMON /XSECTI/ XSMAX(6,5,35),XSTEMP(6,5,35),NPTSFX(5,35),          E10770
     *                NFILEX(5,35),NLIMX                                  E10780
      COMMON /XSTMPR/ PF,TF,PDX(6,5,35),DVXPR(5,35),IXBIN(5,35),          E10790
     *                IXSBN(5,35)                                         E10800
      COMMON /XSHEAD/ HEADT1(35)                                          E10810
      COMMON /FLFORM/ CFORM                                               E10820
C                                                                         E10830
      CHARACTER*10 XSFILE,XSNAME,ALIAS                                    E10840
      CHARACTER HEADT1*100,CFORM*11,XSFIL*10,XSTMP*4,XSNUM*3              E10850
      LOGICAL OP,OPCL                                                     E10860
      DIMENSION FILHDR(2),PNLHDR(2),FILHDS(2)                             E10870
C                                                                         E10880
      EQUIVALENCE (IHIRAC,FSCDI1(1)) , (ILBLF4,FSCDI1(2)),                E10890
     *            (IXSCNT,FSCDI1(3)) , (IAERSL,FSCDI1(4)),                E10900
     *            (IEMIT,FSCDI1(5)) , (ISCHDR,FSCDI1(6)),                 E10910
     *            (JRAD,FSCDI1(9)) , (XSCID,FSCDI1(12)),                  E10920
     *            (XHWHM,FSCDI1(13)) , (IDABS,FSCDI1(14)),                E10930
     *            (IATM,FSCDI1(15)) , (LAYR1,FSCDI1(16)),                 E10940
     *            (XI1(1),FILHDR(1)) , (V1PX,PNLHDR(1)),                  E10950
     *            (XID(1),FILHDS(1))                                      E10960
C                                                                         E10970
      DATA HWJ / 16. /,DXJ / 0.02 /,NJ / 801 /,NJMX / 851 /,              E10980
     *     SMPLJ / 4. /,XSCAL / 0. /                                      E10990
C                                                                         E11000
      DATA XSTMP / 'TMPX'/                                                E11010
      DATA IFILEO,JFILEO / 93,91 /                                        E11020
C                                                                         E11030
C     STANDARD PRESSURE AND TEMPERATURE                                   E11040
C                                                                         E11050
      DATA P0 / 1013. /,T0 / 273.15 /                                     E11060
C                                                                         E11070
C     ASSUMED MEAN HALFWIDTH AT P0 IS 0.10                                E11080
C     DOPPLER VALUES ARE INITIALIZED AT T296.                             E11085
C                                                                         E11090
      DATA HWHM0 / 0.10 /, T296 / 296.0 /                                 E11100
C                                                                         E11110
C     INITIALIZE IXSBN AND TEMPERATURE RATIO                              E11120
C                                                                         E11130
      IXSBN(NS,NI) = 0                                                    E11140
C                                                                         E11150
      FAC1 = 1.0                                                          E11160
      FAC2 = 0.0                                                          E11170
      IF (NMODE.EQ.2) THEN                                                E11180
         FAC2 = (TAVE-XSTEMP(NT1,NS,NI))/                                 E11190
     *               (XSTEMP(NT2,NS,NI)-XSTEMP(NT1,NS,NI))                E11200
         FAC1 = 1.0-FAC2                                                  E11210
      ENDIF                                                               E11220
      PD = PDX(NT1,NS,NI)*FAC1+PDX(NT2,NS,NI)*FAC2                        E11230
C                                                                         E11280
C     NOTE THAT AT THIS POINT, THE CROSS-SECTIONS HAVE BEEN               E11290
C     LINEARLY INTERPOLATED IN TEMPERATURE (HENCE TD=TF),                 E11300
C     AND THE CONVOLUTION WILL BE DONE ONLY FOR PRESSURE                  E11310
C                                                                         E11320
      HWHMF = HWHM0*(PF/P0)*(T0/TF)                                       E11330
      HWHMD = HWHM0*(PD/P0)*(T0/TF)                                       E11340
C                                                                         E11341
C     SET MINIMUM HALF-WIDTH TO DOPPLER                                   E11342
C                                                                         E11343
      HDOPLR=XDOPLR(NS,NI)*SQRT(TF/T296)                                  E11344
      HWHMD=MAX(HDOPLR,HWHMD)                                             E11345
      HWHMSC = HWHMF-HWHMD                                                E11350
      IF (HWHMSC/HWHMD.LT.0.1) GO TO 30                                   E11360
      HWHM = HWHMSC                                                       E11370
C                                                                         E11380
C     BOUND AT THIS POINT IS THE WAVENUMBER VALUE                         E11390
C     OF HALF THE SCANNING FUNCTION                                       E11400
C                                                                         E11410
      DVO = HWHMF/2.0                                                     E11420
      IF (HWHMSC/2.0.LT.DVFX(NS,NI)) GO TO 30                             E11430
      SAMPLE = HWHM/DVO                                                   E11440
      XHWHM = HWHM                                                        E11450
C                                                                         E11460
C     OPEN FILE AND SET FILHDR                                            E11470
C                                                                         E11480
      INQUIRE (FILE='TMPXBIN',OPENED=OP)                                  E11490
      INQUIRE (UNIT=IFILEO,OPENED=OPCL)                                   E11500
      IF (.NOT.OP) THEN                                                   E11510
         IF (OPCL) CLOSE (IFILEO)                                         E11520
         OPEN (IFILEO,FILE='TMPXBIN',STATUS='UNKNOWN',FORM=CFORM)         E11530
         REWIND IFILEO                                                    E11540
         CALL BUFOUT (IFILEO,FILHDS(1),NFHDRF)                            E11550
         REWIND IFILEO                                                    E11560
         CALL BUFIN (IFILEO,IEOF,FILHDR(1),NFHDRF)                        E11570
         READ (HEADT1(NI),'(10A8)') XI1                                   E11580
      ENDIF                                                               E11590
      REWIND IFILEO                                                       E11600
C                                                                         E11610
      NLIMX = 510                                                         E11630
      IOTPAN = 1                                                          E11640
      LPMAX = 0                                                           E11650
C                                                                         E11660
      NPAN = (NPTSFX(NS,NI)+9)/NLIMX+1                                    E11670
      V1PX = V1FX(NS,NI)                                                  E11680
      NPANEL = -1                                                         E11690
      NSKIP = 0                                                           E11700
C                                                                         E11710
      DO 20 NP = 1, NPAN                                                  E11720
         IF (NP.NE.1) NPANEL = 1                                          E11730
         NMAX = 510                                                       E11740
         IF (NP.EQ.1) NMAX = 500                                          E11750
         V2PX = V1PX+FLOAT(LPMAX+NMAX-1)*DVFX(NS,NI)                      E11760
         V2PX = MIN(V2PX,V2FX(NS,NI))                                     E11770
         NMAX = ((V2PX-V1PX)/DVFX(NS,NI)+ONEPL)-LPMAX                     E11780
C                                                                         E11790
         IMAX = MIN(NMAX,NLIMX)                                           E11800
C                                                                         E11810
         CALL CPUTIM (TIME0)                                              E11820
         CALL XSECIN (NPANEL,NI,NS,NT1,NT2,NMODE,NSKIP,IMAX,NEOF)         E11830
         CALL CPUTIM (TIME)                                               E11840
         TXSRDF = TXSRDF+TIME-TIME0                                       E11850
         TXSCNV = TXSCNV-TIME+TIME0                                       E11860
C                                                                         E11870
         DO 10 JI = 1, NMAX                                               E11880
            JJ = JI+4                                                     E11890
            RBX(JI+LPMAX) = FAC1*RDX1(JJ)                                 E11900
            IF (NMODE.EQ.2) RBX(JI+LPMAX) = RBX(JI+LPMAX)+FAC2*RDX2(JJ)   E11910
   10    CONTINUE                                                         E11920
C                                                                         E11930
         LPMAX = LPMAX+NMAX                                               E11940
         IOTPAN = IOTPAN+1                                                E11950
C                                                                         E11960
         IF (NP.EQ.1) THEN                                                E11970
            V1B = V1FX(NS,NI)                                             E11980
            V2B = V2FX(NS,NI)                                             E11990
            DVB = DVFX(NS,NI)                                             E12000
            XSCID = -99                                                   E12010
            ISCHDR = 0                                                    E12020
            IEMIT = 0                                                     E12030
            CALL BUFOUT (IFILEO,FILHDR(1),NFHDRF)                         E12040
         ENDIF                                                            E12050
C                                                                         E12060
         IF (IOTPAN.EQ.5.OR.NP.EQ.NPAN) THEN                              E12070
            IOTPAN = 1                                                    E12080
            DVPX = DVFX(NS,NI)                                            E12090
            NLIMPX = LPMAX                                                E12100
            CALL BUFOUT (IFILEO,PNLHDR(1),NPHDRF)                         E12110
            CALL BUFOUT (IFILEO,RBX(1),NLIMPX)                            E12120
            LPMAX = 0                                                     E12130
            V1PX = V2PX+DVFX(NS,NI)                                       E12140
         ENDIF                                                            E12150
C                                                                         E12160
   20 CONTINUE                                                            E12170
C                                                                         E12180
      NSHIFT = 0                                                          E12190
C                                                                         E12200
      V1 = V1FX(NS,NI)                                                    E12210
      V2 = V2FX(NS,NI)                                                    E12220
C                                                                         E12230
      JEMIT = 0                                                           E12240
      JABS = 0                                                            E12250
C                                                                         E12260
      HWFSV = HWF                                                         E12270
      DXFSV = DXF                                                         E12280
      NFSV = NF                                                           E12290
      NFMXSV = NFMAX                                                      E12300
C                                                                         E12310
      HWF = HWJ                                                           E12320
      DXF = DXJ                                                           E12330
      NF = NJ                                                             E12340
      NFMAX = NJMX                                                        E12350
CCP   XSCALE = XSCAL                                                      E12360
      CALL SLRENZ (XF)                                                    E12370
C                                                                         E12380
      WRITE (XSNUM,'(I1,I2.2)') NS,NI                                     E12390
      XSFIL = XSTMP//XSNUM                                                E12400
C                                                                         E12410
      INQUIRE (FILE=XSFIL,OPENED=OP)                                      E12420
      INQUIRE (UNIT=JFILEO,OPENED=OPCL)                                   E12430
      IF (.NOT.OP.AND.OPCL) CLOSE (JFILEO)                                E12440
      IF (.NOT.OP) OPEN (JFILEO,FILE=XSFIL,STATUS='UNKNOWN',FORM=CFORM)   E12450
      REWIND JFILEO                                                       E12460
C                                                                         E12470
      CALL XSCNVN (IFILEO,JFILEO,NS,NI)                                   E12480
C                                                                         E12490
      CLOSE (IFILEO)                                                      E12500
      CLOSE (JFILEO)                                                      E12510
C                                                                         E12520
      HWF = HWFSV                                                         E12530
      DXF = DXFSV                                                         E12540
      NF = NFSV                                                           E12550
      NFMAX = NFMXSV                                                      E12560
C                                                                         E12570
      IXSBN(NS,NI) = 1                                                    E12580
C                                                                         E12590
   30 RETURN                                                              E12600
C                                                                         E12610
      END                                                                 E12650
      SUBROUTINE XSCNVN (IFILE,JFILE,NS,NI)                               E12660
C                                                                         E12670
      IMPLICIT REAL*8           (V)                                     ! E12680
C                                                                         E12690
C     DRIVER FOR CONVOLVING SPECTRUM WITH INSTRUMENTAL SCANNING FUNCTIO   E12700
C                                                                         E12710
      COMMON /XSCONV/ S(2050),R1(1025)                                    E12720
C                                                                         E12730
      CHARACTER*8      XID,       HMOLID,      YID   
      Real*8               SECANT,       XALTZ
C                                                                         E12750
      COMMON /SCNHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       E12760
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1C,V2C,TBOUND,   E12770
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    E12780
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   E12790
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       E12800
     *               DVSC,XDUM,V1SHFT                                     E12810
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 E12820
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,tdum(8)
      COMMON /RSCAN/ V1I,V2I,DVI,NNI                                      E12840
      COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,               E12850
     *                HWF3,DXF3,NX3,N3MAX                                 E12860
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         E12870
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        E12880
     *              NLTEFL,LNFIL4,LNGTH4                                  E12890
      COMMON /XSCINF/ HWHM,JEMIT,JFN,SAMPLE,SCANID,NPTS,XF(851)           E12900
C                                                                         E12910
      DIMENSION FILHDR(2)                                                 E12920
      DIMENSION SUMR(4)                                                   E12930
C                                                                         E12940
      EQUIVALENCE (FILHDR(1),XID(1)) , (FSCDID(5),IEMIT),                 E12950
     *            (FSCDID(6),ISCHDR) , (FSCDID(12),XSCID),                E12960
     *            (FSCDID(13),XHWHM) , (FSCDID(14),IDABS),                E12970
     *            (FSCDID(16),LAYR1)                                      E12980
C                                                                         E12990
C     IUNIT INPUT FILE                                                    E13000
C     JUNIT OUTPUT FILE                                                   E13010
C                                                                         E13020
      IUNIT = IFILE                                                       E13030
      JUNIT = JFILE                                                       E13040
      NREN = 0                                                            E13050
      IDABS = 0                                                           E13060
      JVAR = 0                                                            E13070
      IPRT = 0                                                            E13080
C                                                                         E13090
      IF (JEMIT.LT.0) THEN                                                E13100
         JABS = 1                                                         E13110
         JEMIT = 0                                                        E13120
         IDABS = -1                                                       E13130
      ENDIF                                                               E13140
      IDABST = IDABS                                                      E13150
      IFILST = 1                                                          E13160
      NIFILS = 9999                                                       E13170
C                                                                         E13180
      SUMOUT = 0.                                                         E13190
      SMIN = 999999.                                                      E13200
      SMAX = -99999.                                                      E13210
      DVOSAV = 0.                                                         E13220
      SUMR(1) = SUMOUT                                                    E13230
      SUMR(2) = SMIN                                                      E13240
      SUMR(3) = SMAX                                                      E13250
      SUMR(4) = DVOSAV                                                    E13260
C                                                                         E13270
      REWIND IUNIT                                                        E13280
      CALL BUFIN (IUNIT,IEOF,FILHDR(1),NFHDRF)                            E13290
      IF (IEOF.EQ.0) GO TO 50                                             E13300
C                                                                         E13310
      DVSAV = DV                                                          E13320
      IDABS = IDABST                                                      E13330
C                                                                         E13340
      ISCAN = ISCHDR                                                      E13350
      JTREM = 3                                                           E13360
C                                                                         E13370
C     JTREM=3   SCANFN CONVOLVED WITH OPTICAL DEPTH                       E13380
C                                                                         E13390
      DVI = DV                                                            E13400
C                                                                         E13410
C     BOUND AT THIS POINT IS THE WAVENUMBER VALUE                         E13420
C     OF HALF THE SCANNING FUNCTION                                       E13430
C                                                                         E13440
      BOUND = HWF*HWHM                                                    E13450
      DV = DVO                                                            E13460
      V1C = V1                                                            E13470
      V2C = V2                                                            E13480
      XHWHM = HWHM                                                        E13490
      IEMIT = 0                                                           E13500
      CALL BUFOUT (JUNIT,FILHDR(1),NFHDRF)                                E13510
      NBOUND = (2.*HWF)*SAMPLE+0.01                                       E13520
C                                                                         E13530
C     BOUND AT THIS POINT IS THE WAVENUMBER VALUE                         E13540
C     OF THE FULL SCANNING FUNCTION                                       E13550
C                                                                         E13560
      BOUND = FLOAT(NBOUND)*DVO/2.                                        E13570
C                                                                         E13580
      NXPAN = 500                                                         E13590
      NLO = NBOUND+1                                                      E13600
      NLIMF = NLO+NXPAN-NSHIFT                                            E13610
      NHI = NLIMF+NSHIFT-1                                                E13620
      MAXF = NLIMF+2*NBOUND                                               E13630
C                                                                         E13640
      TIMRDF = 0.                                                         E13650
      TIMCNV = 0.                                                         E13660
      TIMPNL = 0.                                                         E13670
      IEOFSC = 1                                                          E13680
      SUMIN = 0.                                                          E13690
      DO 10 I = 1, MAXF                                                   E13700
         R1(I) = 0.                                                       E13710
   10 CONTINUE                                                            E13720
      INIT = 0                                                            E13730
      IDATA = -1                                                          E13740
      VFT = V1-2.*BOUND                                                   E13750
      VBOT = V1-BOUND                                                     E13760
      VTOP = V2+BOUND                                                     E13770
C                                                                         E13780
   20 CALL CPUTIM (TIME0)                                                 E13790
      IF (IEOFSC.LE.0) GO TO 40                                           E13800
      CALL RDSCAN (S,JTREM,IUNIT,ISCAN,IPRT)                              E13810
C                                                                         E13820
      CALL CPUTIM (TIME)                                                  E13830
      TIMRDF = TIMRDF+TIME-TIME0                                          E13840
C                                                                         E13850
      IF (IEOFSC.LE.0) GO TO 40                                           E13860
      CALL SHRKSC (INIT,HWHM)                                             E13870
C                                                                         E13880
C     SHRKSC MAY SHRINK (COMPRESS) THE DATA;                              E13890
C     DVI IS MODIFIED ACCORDINGLY                                         E13900
C                                                                         E13910
   30 CONTINUE                                                            E13920
      CALL CONVSC (S,HWHM,R1,XF)                                          E13930
C                                                                         E13940
      IF (IPANEL.EQ.0) GO TO 20                                           E13950
C                                                                         E13960
   40 CALL PNLCNV (R1,JUNIT,SUMR,NPTS,NS,NI)                              E13970
      IF ((ISTOP.NE.1).AND.(IEOFSC.GT.0)) GO TO 30                        E13980
      IF (ISTOP.NE.1) GO TO 40                                            E13990
      CALL CPUTIM (TIME)                                                  E14000
C                                                                         E14010
      SUMIN = SUMIN*DVSAV                                                 E14020
C                                                                         E14030
      IF (IEOFSC.EQ.1) CALL SKIPFL (1,IUNIT,IEOFSC)                       E14040
C                                                                         E14050
      IEOFT = IEOFT+1                                                     E14060
C                                                                         E14070
      SUMOUT = SUMR(1)                                                    E14080
      SMIN = SUMR(2)                                                      E14090
      SMAX = SUMR(3)                                                      E14100
      DVOSAV = SUMR(4)                                                    E14110
C                                                                         E14120
      SUMOUT = SUMOUT*DVOSAV                                              E14130
C                                                                         E14140
   50 RETURN                                                              E14150
C                                                                         E14160
      END                                                                 E14170
      SUBROUTINE PNLCNV (R1,JFILE,SUMR,NPTS,NS,NI)                        E14180
C                                                                         E14190
      IMPLICIT REAL*8           (V)                                     ! E14200
C                                                                         E14210
C     SUBROUTINE PNLCNV OUTPUTS THE RESULTS OF THE CONVOLUTION            E14220
C     TO FILE JFILE                                                       E14230
C                                                                         E14240
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         E14250
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        E14260
     *              NLTEFL,LNFIL4,LNGTH4                                  E14270
      COMMON /SSUBS/ VFT,VBOT,VTOP,V1,V2,DVO,NLIMF,NSHIFT,MAXF,ILO,IHI,   E14280
     *               NLO,NHI,RATIO,SUMIN,IRATSH,SRATIO,IRATM1,NREN,       E14290
     *               DVSC,XDUM,V1SHFT                                     E14300
      COMMON /CONTRL/ IEOFSC,IPANEL,ISTOP,IDATA,JVAR,JABS                 E14310
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,tdum(8)
      COMMON /SPANEL/ V1P,V2P,DV,NLIM                                     E14330
      COMMON /XSTMPR/ PF,TF,PDX(6,5,35),DVXPR(5,35),IXBIN(5,35),          E14340
     *                IXSBN(5,35)                                         E14350
      COMMON /XSHEAD/ HEADT1(35)                                          E14360
      CHARACTER HEADT1*100                                                E14370
      DIMENSION PNLHDR(2)                                                 E14380
      DIMENSION R1(*),SUMR(*)                                             E14390
C                                                                         E14400
      EQUIVALENCE (PNLHDR(1),V1P)                                         E14410
C                                                                         E14420
      CALL CPUTIM (TIME0)                                                 E14430
C                                                                         E14440
      SUMOUT = SUMR(1)                                                    E14450
      SMIN = SUMR(2)                                                      E14460
      SMAX = SUMR(3)                                                      E14470
      DV = DVO                                                            E14480
      ISTOP = 0                                                           E14490
      NNHI = (V2-VFT)/DV+1.5                                              E14500
      IF (NHI.GE.NNHI) THEN                                               E14510
         ISTOP = 1                                                        E14520
         NHI = NNHI                                                       E14530
      ENDIF                                                               E14540
      NLIM = NHI-NLO+1                                                    E14550
      V1P = VFT+FLOAT(NLO-1)*DV                                           E14560
      V2P = VFT+FLOAT(NHI-1)*DV                                           E14570
C                                                                         E14580
C     V1P IS FIRST FREQ OF PANEL                                          E14590
C     V2P IS LAST  FREQ OF PANEL                                          E14600
C                                                                         E14610
      CALL BUFOUT (JFILE,PNLHDR(1),NPHDRF)                                E14620
      CALL BUFOUT (JFILE,R1(NLO),NLIM)                                    E14630
C                                                                         E14640
      VFT = VFT+FLOAT(NLIMF-1)*DV                                         E14650
      DVXPR(NS,NI) = DV                                                   E14660
      NLIMHI = NLIM+NLO-1                                                 E14670
      DO 10 I = NLO, NLIMHI                                               E14680
         SMIN = MIN(SMIN,R1(I))                                           E14690
         SMAX = MAX(SMAX,R1(I))                                           E14700
         SUMOUT = SUMOUT+R1(I)                                            E14710
   10 CONTINUE                                                            E14720
      IF (ISTOP.EQ.1) GO TO 40                                            E14730
      JF = 1                                                              E14740
      DO 20 J = NLIMF, MAXF                                               E14750
         R1(JF) = R1(J)                                                   E14760
         JF = JF+1                                                        E14770
   20 CONTINUE                                                            E14780
      DO 30 J = JF, MAXF                                                  E14790
         R1(J) = 0.                                                       E14800
   30 CONTINUE                                                            E14810
      NLIMF = 511                                                         E14820
      NLO = NSHIFT+1                                                      E14830
      NHI = NLIMF+NSHIFT-1                                                E14840
   40 SUMR(1) = SUMOUT                                                    E14850
      SUMR(2) = SMIN                                                      E14860
      SUMR(3) = SMAX                                                      E14870
      SUMR(4) = DVO                                                       E14880
      CALL CPUTIM (TIME)                                                  E14890
      TIMPNL = TIMPNL+TIME-TIME0                                          E14900
C                                                                         E14910
      RETURN                                                              E14920
C                                                                         E14930
      END                                                                 E14940
      SUBROUTINE SLRENZ (XF)                                              E14950
CCP   SUBROUTINE SLRENZ (XF,XSCALE)                                       E14951
C                                                                         E14960
C     SUBROUTINE SLRENZ SETS UP THE LORENZ SCANNING FUNCTION              E14970
C                                                                         E14980
      COMMON /CMSHAP/ HWF,DXF,NF,NFMAX,HWF2,DXF2,NX2,N2MAX,               E14990
     *                HWF3,DXF3,NX3,N3MAX                                 E15000
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         E15010
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        E15020
     *              NLTEFL,LNFIL4,LNGTH4                                  E15030
      DIMENSION XF(*)                                                     E15040
C                                                                         E15050
      PI = 2.*ASIN(1.)                                                    E15060
      XNORM = 1.0/PI                                                      E15070
      DO 10 I = 1, NFMAX                                                  E15080
         XF(I) = 0.                                                       E15090
   10 CONTINUE                                                            E15100
      XF(1) = XNORM                                                       E15110
      SUM = XF(1)                                                         E15120
      DO 20 I = 2, NF                                                     E15130
         X = FLOAT(I-1)*DXF                                               E15140
         XF(I) = XNORM*(1./(1.+X**2))                                     E15150
         SUM = SUM+2.*XF(I)                                               E15160
   20 CONTINUE                                                            E15170
      SUM = SUM*DXF                                                       E15180
C                                                                         E15190
C     RENORMALIZE                                                         E15200
C                                                                         E15210
      XNORM = 1.0/SUM                                                     E15220
      DO 30 I = 1, NF                                                     E15230
         XF(I) = XNORM*XF(I)                                              E15240
   30 CONTINUE                                                            E15250
C                                                                         E15260
C     WRITE(IPR,900) NF,DXF,SUM                                           E15270
C                                                                         E15280
      RETURN                                                              E15290
C                                                                         E15300
      END                                                                 E15310
