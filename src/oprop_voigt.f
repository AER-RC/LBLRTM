C     path: %Source%
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
C     *            USES TABULATED VOIGT ALGORITHM           
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
C                                     M.W.SHEPHARD, W.O.GALLERY, 
C                                     &  S.A. CLOUGH 
C                                                                         B00300
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B00310
C                     131 Hartwell Ave,  Lexington,  MA   02421           B00320
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
      COMMON /CVROPR/ HNAMOPR,HVROPR
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
      COMMON /FNSH/ IFN,F1(NFMX, NUMZ),F2(NFMX, NUMZ),
     $     F3(NFMX, NUMZ), FG(NFMX)

      COMMON /R4SUB/ VLOF4,VHIF4,ILOF4,IST,IHIF4,LIMIN4,LIMOUT,ILAST,     B00820
     *               DPTMN4,DPTFC4,ILIN4,ILIN4T                           B00830
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         B00840
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        B00850
     *              NLTEFL,LNFIL4,LNGTH4                                  B00860
C                                                                         B00870
      PARAMETER (NTMOL=38)   
C                                                                         B00890
      COMMON /ISVECT/ ISO_MAX(NTMOL),SMASSI(ntmol,9)
      COMMON /LNC1/ RHOSLF(ntmol),ALFD1(42,9),SCOR(42,9),ALFMAX,  
     *              BETACR,DELTMP,DPTFC,DPTMN,XKT,NMINUS,NPLUS,NLIN,      B00930
     *              LINCNT,NCHNG,SUMALF,SUMZET,TRATIO,RHORAT,PAVP0,       B00940
     *              PAVP2,RECTLC,TMPDIF,ILC                               B00950
      COMMON /FLFORM/ CFORM                                               B00960
      COMMON /L4TIMG/ L4TIM,L4TMR,L4TMS,L4NLN,L4NLS,LOTHER
      COMMON /IODFLG/ DVOUT

c     Total timing array for layer line-by-line calculation
      common /timing_lay/ time_lay_lbl(20)
C                                                                         B00970
      REAL L4TIM,L4TMR,L4TMS,LOTHER
      CHARACTER*55 CDUM1,PTHODI,PTHODT,PTHRDR
      CHARACTER*10 HFMODL
      CHARACTER CFORM*11,KODLYR*57,PTHODE*55,PTHODD*55                    B00980
      CHARACTER*18 HNAMOPR,HVROPR
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
      DATA I_10/10/
C
      PTHODE = 'ODexact_'
      PTHODD = 'ODdeflt_'
      DATA KODLYR /
     *     '                                                         '/
      DATA HFMODL /'         '/
C                                                                         B01160
      CALL CPUTIM (TIMEH0)                                                B01170
C                                                                         B01180
C     ASSIGN NAME and CVS VERSION NUMBER TO MODULE 
C
      HNAMOPR = '    oprop_voigt.f:'
      HVROPR  = '$Revision$'
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
      ICNTNM = MOD(IXSCNT,I_10)                                             B01200
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
      BOUND =  REAL(NBOUND)*DV/2.                                         B02000
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
      VFT = V1- REAL(NSHIFT)*DV                                           B02460
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
      CALL CNVFNV (VNU,SP,SPPSP,RECALF,R1,R2,R3,ZETAI,IZETA)
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
            VF2 = VFT+2.*DVR4+ REAL(N2R3+4)*DVR3                          B03390
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

c        Fill timing array
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

         DATA i_time_lay/454545/
         
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
            AVALF = SUMALF/ REAL(LINCNT)                                  B03640
            AVZETA = SUMZET/ REAL(LINCNT)                                 B03650
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
 921  FORMAT (2x,'LINF4',3X,2F15.3,15X,2F15.3,
     *        2X,'XSECT ',2X,4F15.3,
     *        2X,'LBLF4 ',2X,4F15.3,
     *        2X,'HIRAC1',2X,5F15.3)
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
     *          'NWDL, ENDFIL, and SHRINK',/,
     *          '0             In HIRAC:  LNCOR, XINT, SHAPEL, ',
     *          'SHAPEG, VERFN, MOLEC, and other loops and ',
     *          'file maintenance within HIRAC',/)
C                                                                         B03900
      END                                                                 B03910
      BLOCK DATA BHIRAC                                                   B03920
C
      PARAMETER (NFPTS=2001,NFMX=1.3*NFPTS, NUMZ = 101)  
c
      COMMON /FNSH/ IFN,F1(NFMX, NUMZ),F2(NFMX, NUMZ),
     $     F3(NFMX, NUMZ), FG(NFMX)
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
   10 CALL BUFIN_sgl(Lnfl,LEOF,rdpnl(1),npnlhd)
      IF (LEOF.EQ.0) THEN                                                 B04320
         IF (NOPR.EQ.0) WRITE (IPR,900)                                   B04330
         IEOF = 1                                                         B04340
         RETURN                                                           B04350
      ENDIF                                                               B04360
C                                                                         B04370
      IF (NREC.GT.LIMIN) STOP 'RDLIN; NREC GT LIMIN'                      B04380
c
      IF (VMAX.LT.VBOT) THEN                                              B04390
         CALL BUFIN_sgl(Lnfl,LEOF,DUM(1),n_one)
         GO TO 10                                                         B04410
      ENDIF                                                               B04420
c
      CALL BUFIN_sgl(Lnfl,LEOF,vlin(1),NWDS) 
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
      COMMON /IFIL/   IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,         B04100
     *                NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        B04110
     *                NLTEFL,LNFIL4,LNGTH4                                  B04120
      COMMON /XSUB/   VBOT,VTOP,VFT,DUM(7)
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)        B05050
      COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,
     *                HWF3,DXF3,NX3,N3MAX
      COMMON /VOICOM/ AVRAT(102),CGAUSS(102),CF1(102),CF2(102),           B05060
     *                CF3(102),CER(102)                                   B05070
C                                                                         B05080
      PARAMETER (NTMOL=38) 
C                                                                         B05100
      COMMON /ISVECT/ ISO_MAX(NTMOL),SMASSI(ntmol,9)
      COMMON /LNC1/ RHOSLF(ntmol),ALFD1(42,9),SCOR(42,9),ALFMAX, 
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
c                                                                          D00540
      character*8 h_lncor1
c
      data h_lncor1/' lncor1 '/
C                                                                         B05240
C     TEMPERATURES FOR LINE COUPLING COEFFICIENTS                         B05250
C                                                                         B05260
      DATA TEMPLC / 200.0,250.0,296.0,340.0 /                             B05270
      DATA HREJ /'0'/,HNOREJ /'1'/
      DATA NWDTH /0/
C                                                                         B05280
      DATA I_1/1/, I_100/100/, I_1000/1000/
C
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
         M = MOD(MOL(I),I_100)                                              B05640
C                                                                         B05650
C     ISO=(MOD(MOL(I),1000)-M)/100   IS PROGRAMMED AS:                    B05660
C                                                                         B05670
         ISO = MOD(MOL(I),I_1000)/100                                       B05680
C
c     check if lines are within allowed molecular and isotopic limits
c
         if (m.gt.ntmol .or. m.lt. 1) then
            call line_exception (1,ipr,h_lncor1,m,nmol,iso,iso_max)
            go to 25
         else if (iso .gt. iso_max(m)) then
            call line_exception (2,ipr,h_lncor1,m,nmol,iso,iso_max)
            go to 25
         endif
C
         MOL(I) = M                                                       B05760
C
         SUI = S(I)*WK(M)                                                 B05770
         IF (JRAD.EQ.1) SUI = SUI*VNU(I)
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
         ALFL = ALFA0I*(RHORAT-RHOSLF(m))+HWHMSI*RHOSLF(m)  
C                                                                         B06260
         IF (IFLAG.EQ.3) ALFL = ALFL*(1.0-GAMMA1*PAVP0-GAMMA2*PAVP2)      B06270
C                                                                         B06280
         ALFAD = VNU(I)*ALFD1(m,iso)                                       B06290
         ZETA = ALFL/(ALFL+ALFAD)                                         B06300
         ZETAI(I) = ZETA                                                  B06310
         FZETA = 100.*ZETA
         IZ = FZETA + ONEPL                                               B06320
         IZETA(I) = IZ                                                    B06330
         ZETDIF = FZETA -  REAL(IZ-1)
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
C     TREAT TRANSITIONS WITH negative EPP AS SPECIAL CASE 
C                                                                         B06500
c>>   an epp value between -0.9999 and 0.  cm-1 is taken as valid 
c
c>>   an epp value of -1. is assumed set by hitran indicating an unknown 
c     value: no temperature correction is performed
c
c>>   for an epp value of less than -1., it is assumed that value has
c     been provided as a reasonable value to be used for purposes of 
c     temperature correction.  epp is set positive
c
         if (epp(i).le.-1.001)  epp(i) = abs(epp(i))

         if (epp(i).le.-0.999)  MEFDP(M) = MEFDP(M)+1 

c     temperature correction:

         if (epp(i) .gt. -0.999) then
            SUI = SUI*SCOR(m,iso)*  
     *              EXP(-EPP(I)*BETACR)*(1.+EXP(-VNU(I)/XKT))
         endif
C                                                                         B06860
         SP(I) = SUI*(1.+GI*PAVP2)                                        B06870
         SPPI = SUI*YI*PAVP0                                              B06880
         SPPSP(I) = SPPI/SP(I)                                            B06890
C
c     SPEAK is used for line rejection (no LCPL lines (sppsp ne 1) are rejected)

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
                  JJ = MAX(JJ,I_1)                                          B06700
                  JJ = MIN(JJ,NPTR4)                                      B06710
                  IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ)) 
     &                .and. sppsp(i).eq.0.) THEN
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
c-----------------------------------------------------------------------
c
      SUBROUTINE CNVFNV (VNU,SP,SPPSP,RECALF,R1,R2,R3,ZETAI,IZETA) 
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
      COMMON /FNSH/ IFN,F1(NFMX, NUMZ),F2(NFMX, NUMZ),F3(NFMX, NUMZ),
     $     FG(NFMX)
c
      DIMENSION VNU(*),SP(*),SPPSP(*),RECALF(*)                           B07530
      DIMENSION R1(*),R2(*),R3(*)                                         B07540
      DIMENSION IZETA(*),ZETAI(*)                                         B07570
C                                                                         B07580
      equivalence (zetdif,A)
c
      CALL CPUTIM (TIME0)                                                 B07590
C                                                                         B07600
      CLC1 = 4./( REAL(NX1-1))                                            B07610
      CLC2 = 16./( REAL(NX2-1))                                           B07620
      CLC3 = 64./( REAL(NX3-1))                                           B07630
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
               ZETDIF = 100.*ZETAI(I)- REAL(IZM-1)                        B07840
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
               ZF1 = ( REAL(JMIN1-2)-ZINT)*ZSLOPE                         B08100
               ZF2 = ( REAL(JMIN2-2)-ZINT*CONF2)*ZSLOPE                   B08110
               ZF3 = ( REAL(JMIN3-2)-ZINT*CONF3)*ZSLOPE                   B08120
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
c                     
C     ************Using new voigt scheme: 
C     ************Interpolate voigt subfunctions to zeta 

c          A  is equivalenced to ZETDIF

                     x3 = depthi*( (1.-A)*F3(IZ3,IZM)+A*F3(IZ3,IZM+1))
                     x2 = depthi*(( 1.-A)*F2(IZ2,IZM)+A*F2(IZ2,IZM+1))
                     x1 = depthi*( (1.-A)*F1(IZ1,IZM)+A*F1(IZ1,IZM+1))

                     R3(J3) = R3(J3)+ x3
                     R2(J2) = R2(J2)+ x2
                     R1(J1) = R1(J1)+ x1

 10               CONTINUE                                                B08290
C     
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
                  DO 20 J1 = JMIN1, JMAX1                                 B08450
c
                     J2 = J1-J2SHFT                                       B08460
                     J3 = J1-J3SHFT                                       B08470
                     ZF3 = ZF3+ZSLOPE                                     B08480
                     ZF2 = ZF2+ZSLOPE                                     B08490
                     ZF1 = ZF1+ZSLOPE                                     B08500
                     IZ3 = ABS(ZF3)+1.5                                   B08510
                     IZ2 = ABS(ZF2)+1.5                                   B08520
                     IZ1 = ABS(ZF1)+1.5                                   B08530
c                     
C     ************Using new voigt scheme: 
C     ************Interpolate voigt subfunctions to zeta 

                  x3 = DEPTHI*( (1.-A)*F3(IZ3,IZM)+A*F3(IZ3,IZM+1))
                  x2 = DEPTHI*( (1.-A)*F2(IZ2,IZM)+A*F2(IZ2,IZM+1))
                  x1 = DEPTHI*( (1.-A)*F1(IZ1,IZM)+A*F1(IZ1,IZM+1))

                  y3 = dptrat3*x3*ZF3
                  y2 = dptrat2*x2*ZF2
                  y1 = dptrat1*x1*ZF1

                  R3(J3) = R3(J3) + x3 + y3
                  R2(J2) = R2(J2) + x2 + y2
                  R1(J1) = R1(J1) + x1 + y1

 20               CONTINUE                                                B08580
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
         V2P = VFT+ REAL(NHI-1)*DVP
         IF (V2P.LT.V2) THEN
            V2P = V2P+DVP
            NHI = NHI+1
         ENDIF
         ISTOP = 1                                                        B10640
      ELSE
         V2P = VFT+ REAL(NHI-1)*DV
      ENDIF
      NLIM = NHI-NLO+1
      V1P = VFT+ REAL(NLO-1)*DV
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
         VI = VFT+ REAL(NPTSI1-1)*DV                                      B11100
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
      VFT = VFT+ REAL(NLIM1-1)*DV                                         B11390
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
      V2PO = V1PO+ REAL(LIMOUT)*DVOUT                                     B12860
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
         V2PO = V1PO+ REAL(NLIM2-1)*DVOUT                                 B12920
         IF (V2PO.LT.V2) THEN
            V2PO = V2PO+DVOUT
            NLIM2 = NLIM2+1
         ENDIF
         ILAST = 1
         IF (V2PO.GT.V2P-DVP) THEN
            NLIM2 = ((V2P-DVP-V1PO)/DVOUT) + 1.
            V2PO = V1PO+ REAL(NLIM2-1)*DVOUT
            IF (V2PO+DVOUT.LT.V2P-DVP) THEN
               NLIM2 = NLIM2+1
               V2PO = V2PO+DVOUT
            ENDIF
            ILAST = 0
         ENDIF
      ELSE
         NLIM2 = LIMOUT
         V2PO = V1PO+ REAL(NLIM2-1)*DVOUT                                 B12930
         IF (V2PO.GT.V2P-DVP) THEN                                        B12940
            NLIM2 = ((V2P-DVP-V1PO)/DVOUT) + 1.                           B12950
            V2PO = V1PO+ REAL(NLIM2-1)*DVOUT                              B12960
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
         FJJ = FJ1DIF+RATDV* REAL(II-1)                                   B13130
         JJ  =  INT(FJJ)-2                                                B13140
         JP  = (FJJ- REAL(JJ))*100.-199.5                                 B13150
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
         IF ((V1PO+ REAL(LIMOUT)*DVOUT).GT.(V2P-DVP)) NPPANL = 1          B13300
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
         X =  REAL(JJ-1)*DXF1                                             B14070
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
         X =  REAL(JJ-1)*DXF2                                             B14230
         XSQ = X*X                                                        B14240
         F2(JJ) = RECPI*(Q1FN(XSQ)-Q2FN(XSQ))                             B14250
         SUM = SUM+F2(JJ)*2.                                              B14260
   40 CONTINUE                                                            B14270
      J1LIMP = J1LIM+1                                                    B14280
      DO 50 JJ = J1LIMP, NX2                                              B14290
         X =  REAL(JJ-1)*DXF2                                             B14300
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
         X =  REAL(JJ-1)*DXF3                                             B14460
         XSQ = X*X                                                        B14470
         F3(JJ) = RECPI*(Q2FN(XSQ)-Q3FN(XSQ))                             B14480
         SUM = SUM+F3(JJ)*2.                                              B14490
   70 CONTINUE                                                            B14500
      J2LIMP = J2LIM+1                                                    B14510
      DO 80 JJ = J2LIMP, NX3                                              B14520
         X =  REAL(JJ-1)*DXF3                                             B14530
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
      FLN2 =  LOG(2.)                                                     B14730
      RECPI = 1./(2.*ASIN(1.))                                            B14740
      FGNORM = SQRT(FLN2*RECPI)                                           B14750
      TOTAL = 0.                                                          B14760
      DO 10 I = 1, N1MAX                                                  B14770
         FG(I) = 0.                                                       B14780
   10 CONTINUE                                                            B14790
      FG(1) = FGNORM*FGAUSS(0.)                                           B14800
      SUM = FG(1)                                                         B14810
      DO 20 JJ = 2, NX1                                                   B14820
         X =  REAL(JJ-1)*DXF1                                             B14830
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
c
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
      P = ( REAL(IP-1)+VFT/DV)*2.                                         B17040
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
         VI = VFT+DVR3* REAL(I-1)                                         B17710
         J = (VI-V1A)*RECDVA+ONEPL                                        B17720
         VJ = V1A+DVA* REAL(J-1)                                          B17730
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
      DATA I_1/1/
C
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
               INTVLS = MAX(INTVLS,I_1)                                     B19290
               VINEW = VI+DVI* REAL(INTVLS)                               B19300
            ELSE                                                          B19310
               VINEW = ABS(VINEW)                                         B19320
               INTVLS = (VINEW-VI)/DVI                                    B19330
               INTVLS = MAX(INTVLS,I_1)                                     B19340
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
               INTVLS = MAX(INTVLS,I_1)                                     B19480
               VINEW = VI+DVI* REAL(INTVLS)                               B19490
            ELSE                                                          B19500
               VINEW = ABS(VINEW)                                         B19510
               INTVLS = (VINEW-VI)/DVI                                    B19520
               INTVLS = MAX(INTVLS,I_1)                                     B19530
            ENDIF                                                         B19540
            XVINEW = VINEW                                                B19550
C                                                                         B19560
            RDNEXT = XVINEW*XMINUS/XPLUS                                  B19570
C                                                                         B19580
         ELSE                                                             B19590
            IF (VINEW.GE.0.0) THEN                                        B19600
               VINEW = VI+(FACT1*XVI)                                     B19610
               INTVLS = (VINEW-VI)/DVI                                    B19620
               INTVLS = MAX(INTVLS,I_1)                                     B19630
               VINEW = VI+DVI* REAL(INTVLS)                               B19640
            ELSE                                                          B19650
               VINEW = ABS(VINEW)                                         B19660
               INTVLS = (VINEW-VI)/DVI                                    B19670
               INTVLS = MAX(INTVLS,I_1)                                     B19680
            ENDIF                                                         B19690
            XVINEW = VINEW                                                B19700
C                                                                         B19710
            RDNEXT = XVINEW                                               B19720
         ENDIF                                                            B19730
      ELSE                                                                B19740
         IF (VINEW.GE.0.0) THEN                                           B19750
            VINEW = VI+(FACT1*XVI)                                        B19760
            INTVLS = (VINEW-VI)/DVI                                       B19770
            INTVLS = MAX(INTVLS,I_1)                                        B19780
            VINEW = VI+DVI* REAL(INTVLS)                                  B19790
         ELSE                                                             B19800
            VINEW = ABS(VINEW)                                            B19810
            INTVLS = (VINEW-VI)/DVI                                       B19820
            INTVLS = MAX(INTVLS,I_1)                                        B19830
         ENDIF                                                            B19840
         XVINEW = VINEW                                                   B19850
C                                                                         B19860
         RDNEXT = XVI                                                     B19870
      ENDIF                                                               B19880
C                                                                         B19890
      RDEL = (RDNEXT-RADFNI)/ REAL(INTVLS)                                B19900
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
      PARAMETER (NTMOL=38)
c
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4
      COMMON /ISVECT/ ISO_MAX(NTMOL),SMASSI(ntmol,9)
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
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      DIMENSION SCOR(42,9),RHOSLF(*),ALFD1(42,9)    
      COMMON /SMOLEC/ W(42,9),ND(42,9),FAD                                C00230
      COMMON /XMOLEC/ NV(42),IVIB(42,2,9),XR(42),ROTFAC(42),QV0(42)
      COMMON /MOLNAM/ MOLID(0:NTMOL)                                      C00260
      CHARACTER*6 MOLID                                                   C00270
C                                                                         C00280
C     IS EQUIV. TO THE FOLLOWING DIMENSION AND EQUIVALENT STATEMENTS      C00290
C                                                                         C00300
      DIMENSION IV(2)
      EQUIVALENCE (IV(1),IVIB(1,1,1))                                     C00320
C                                                                         C00330
      DATA MDIM / 42 /,NVDIM / 9 /
C                                                                         C00350
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
c
      IF (IND.EQ.1) THEN 
c
         DO 10 M = 1, NMOL                                                C00670
            READ (MOLID(M),900) HMOLID(M)                                 C00680
 10      CONTINUE                                                         C00690
C     
         FLN2 =  LOG(2.)                                                  C00710
         FAD = FLN2*2.*AVOGAD*BOLTZ/(CLIGHT*CLIGHT)
         XKT0 = TEMP0/RADCN2                                              C00730
C     
         RETURN                                                           C00860
      ELSE                                                                C00870
C     
         RHORAT = (P/P0)*(TEMP0/TEMP)                                     C00890
         XKT = TEMP/RADCN2                                                C00900
C     
c     call tips to get partition sum correction to intensities
c
         call tips_2003(nmol,iso_max,temp,scor)
c
         DO 50 M = 1, NMOL
            DO 40 ISO = 1, ISO_MAX(M)
               RHOSLF(m) = RHORAT*WK(M)/WTOT                    
               ALFD1(m,iso) = SQRT(FAD*TEMP/SMASSI(m,iso))           
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
c            1        2        3        4        5        6        7
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
c            8        9       10       11       12       13       14
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
c           15       16       17       18       19       20       21
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
c           22       23       24       25       26       27       28
C           N2      HCN    CH3Cl     H2O2     C2H2     C2H6      PH3      C01990
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
c           29       30       31       32       33       34       35  
C         COF2      SF6      H2S    HCOOH      HO2        O   ClONO2      C02150
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
c           36       37       38
C          NO+     HOBr     C2H4      ???      ???      ???      ??? 
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
c  ****************************************
      BLOCK DATA Isotop
c  ****************************************
c
      PARAMETER (NMOL=38)
      COMMON /ISVECT/ ISO_MAX(NMOL),SMASS(nmol,9)
      common /iso_id/ iso_82(97)
c
c    The number of isotopes for a particular molecule:
      DATA (ISO_MAX(I),I=1,NMOL)/
c     H2O, CO2, O3, N2O, CO, CH4, O2,
     +  6,   9,  9,   5,  6,   3,  3,
c      NO, SO2, NO2, NH3, HNO3, OH, HF, HCl, HBr, HI,
     +  3,   2,   1,   2,    1,  3,  1,   2,   2,  1,
c     ClO, OCS, H2CO, HOCl, N2, HCN, CH3Cl, H2O2, C2H2, C2H6, PH3
     +  2,   5,    3,    2,  1,   3,     2,    1,    2,    1,   1,
c     COF2, SF6, H2S, HCOOH, HO2, O, ClONO2, NO+, HOBr, C2H4
     +  1,   1,   3,     1,   1,  1,     2,    1,    2,    2/
c
      DATA ISO_82/
c       H2O
     +   161,181,171,162,182,172,
c       CO2
     +  626,636,628,627,638,637,828,728,727,
c       O3
     +  666,668,686,667,676,886,868,678,768,
c       N2O
     +  446,456,546,448,447,
c       CO,                 CH4
     +  26,36,28,27,38,37,  211,311,212,
c       O2,        NO,        SO2
     +  66,68,67,  46,56,48  ,626,646,
c       NO2,   NH3,        HNO3
     +  646,   4111,5111,  146,
c       OH,        HF,  HCl,    HBr,    HI
     +  61,81,62,  19,  15,17,  19,11,  17,
c       ClO,    OCS,                 H2CO
     +  56,76,  622,624,632,623,822,  126,136,128,
c       HOCl,     N2,  HCN
     +  165,167,  44,  124,134,125
c      CH3Cl,    H2O2,  C2H2,       C2H6,  PH3
     +, 215,217,  1661,  1221,1231,  1221,  1111,
c       COF2, SF6, H2S,           HCOOH,  HO2, O,   ClONO2      NO+
     +  269,  29,  121,141,131,   126,    166, 6,   5646,7646,  46,
c       HOBr,      C2H4
     +  169,161,   221,231/  
c
C                                                                         C03620
C     MOLECULAR MASSES FOR EACH ISOTOPE                                   C03630
C                                                                         C03640
      data (smass(1,i),i=1,6)
C  H2O:   161,   181,   171,   162,   182,   172                      
     * /  18.01, 20.01, 19.01, 19.01, 21.02, 20.02/                   
      data (smass(2,i),i=1,9)
C  CO2:   626,   636,   628,   627,   638,   637,   828,   728,   727    ok 
     * /  43.99, 44.99, 45.99, 44.99, 47.00, 46.00, 48.00, 47.00, 46.00/ ok 
      data (smass(3,i),i=1,9)
C   O3:   666,   668,   686    667    676    886    868    678    768    ok
     * /  47.98, 49.99, 49.99, 48.99, 48.99, 51.99, 51.99, 50.99, 50.99/ ok
      data (smass(4,i),i=1,5)
C  N2O:   446,   456,   546,   448,   447                                ok
     * /  44.00, 45.00, 45.00, 46.00, 45.00/                             ok
      data (smass(5,i),i=1,6)
C   CO:   26,    36,    28,    27,    38     37                          ok
     * /  27.99, 28.99, 29.99, 29.00, 31.00, 30.00/                      ok
      data (smass(6,i),i=1,3)
C  CH4:   211,   311,   212
     * /  16.03, 17.03, 17.03/
      data (smass(7,i),i=1,3)
C   O2:    66,    68,    67     
     * /  31.99, 33.99, 32.99/   
      data (smass(8,i),i=1,3)
C   NO:    46,    56,    48,
     * /  30.00, 31.00, 32.00/
      data (smass(9,i),i=1,2)
C  SO2:   626,   646  
     * /  63.96, 65.96/    
      data (smass(10,i),i=1,1)
C  NO2:   646
     * /  45.99/
      data (smass(11,i),i=1,2)
C  NH3:  4111,  5111;
     * /  17.03, 18.02/
      data (smass(12,i),i=1,1)
C HNO3:
     * /  62.99/
      data (smass(13,i),i=1,3)
C   OH:    61,    81,    62
     * /   17.00, 19.01, 18.01/
      data (smass(14,i),i=1,1)
C   HF:    19           
     * /   20.01/          
      data (smass(15,i),i=1,2)
C  HCL:   15,    17;
     * /  35.98, 37.97/
      data (smass(16,i),i=1,2)
C  HBr:   19,    11;
     * /  79.92, 81.92/
      data (smass(17,i),i=1,1)
C   HI:   17  
     * /  127.91/
      data (smass(18,i),i=1,2)
C  ClO:   56,    76;
     * /  50.96, 52.96/
      data (smass(19,i),i=1,5)       
C  OCS:   622,   624,   632,   623,   822  
     * /  59.97, 61.96, 60.97, 60.97, 61.97/
      data (smass(20,i),i=1,3)
C H2CO:  126,   136,   128;
     * /  30.01, 31.01, 32.01/
      data (smass(21,i),i=1,2)
C HOCl:  165,   167    
     * /  51.97, 53.97/
      data (smass(22,i),i=1,1)
C   N2:   44;
     * /  28.01/
      data (smass(23,i),i=1,3)
C  HCN:   124,   134,   125,
     * /  27.01, 28.01, 28.01/ 
      data (smass(24,i),i=1,2)
C CH3CL:  215,   217; 
     * /  49.99, 51.99/
      data (smass(25,i),i=1,1)
C  H2O2:  1661;
     * /  34.01/
      data (smass(26,i),i=1,2)
C  C2H2: 1221,  1231
     * /  26.01, 27.02/
      data (smass(27,i),i=1,1)
C  C2H6: 1221;
     * /  30.05/
      data (smass(28,i),i=1,1)
C   PH3:   1111;
     * /  34.00/
      data (smass(29,i),i=1,1)
C  COF2:  269;
     * /  65.99/
      data (smass(30,i),i=1,1)
C   SF6:   29 
     * /  145.96/
      data (smass(31,i),i=1,3)
C   H2S:  121    141    131;
     * /  33.99, 35.98, 34.99/
      data (smass(32,i),i=1,1)
C HCOOH: 126;
     * /  46.01/ 
      data (smass(33,i),i=1,1)
C   HO2:   166 
     * / 33.00/
      data (smass(34,i),i=1,1)
C     O:   6
     * /  15.99/
      data (smass(35,i),i=1,2)
C ClONO2: 5646   7646;
     * /  96.96, 98.95/
      data (smass(36,i),i=1,1)
C   NO+:   46      
     * /  30.00 /
      data (smass(37,i),i=1,2)
c  HOBr:  169,    161
     * /  95.92,  97.92/
      data (smass(38,i),i=1,2)
c   C2H4: 221,   231;
     * /  44.03, 45.03/
C
      END
c**************************************
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
         VJ = V1P+ REAL(J-JLO)*DVP                                        C12790
         IK = JHILIM+KK-1                                                 C12800
         VK = V1P+ REAL(IK-JLO)*DVP                                       C12810
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
      PARAMETER (NTMOL=38) 
C                                                                         D00080
      COMMON /ISVECT/ ISO_MAX(NTMOL),SMASSI(ntmol,9)
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
C                                                                         D00440
      REAL L4TIM,L4TMR,L4TMS,LOTHER
      DIMENSION MEFDP(64)                                                 D00450
      DIMENSION SCOR(42,9),RHOSLF(ntmol),ALFD1(42,9)
      DIMENSION ALFAL(1250),ALFAD(1250),A(4),B(4),TEMPLC(4)               D00470
      DIMENSION RCDHDR(2),IWD(2),IWD3(2),HLINHD(2),AMOLB(250)             D00480
C                                                                         D00490
      EQUIVALENCE (ALFA0(1),ALFAL(1)) , (EPP(1),ALFAD(1))                 D00500
      EQUIVALENCE (IHIRAC,FSCDID(1)) , (ILBLF4,FSCDID(2))                 D00510
      EQUIVALENCE (VNULO,RCDHDR(1)) , (IWD3(1),VD),                       D00520
     *            (HLINHD(1),HID(1),IWD(1)) , (MOLB(1),AMOLB(1))          D00530
c                                                                          D00540
      character*8 h_linf4
c
      data h_linf4/' linf4  '/
      DATA MEFDP / 64*0 /                                                 D00550
C                                                                         D00560
C     TEMPERATURES FOR LINE COUPLING COEFFICIENTS                         D00570
C                                                                         D00580
      DATA TEMPLC / 200.0,250.0,296.0,340.0 /                             D00590
C                                                                         D00600
      DATA I_100/100/, I_1000/1000/
C
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
         M = MOD(MOLB(I),I_100)                                             D01400
C                                                                         D01410
C     ISO=(MOD(MOLB(I),I_1000)-M)/100   IS PROGRAMMED AS:                   D01420
C                                                                         D01430
         ISO = MOD(MOLB(I),I_1000)/100                                      D01440
C
c     check if lines are within allowed molecular and isotopic limits
c
         if (m.gt.ntmol .or. m.lt. 1) then
            call line_exception (1,ipr,h_linf4,m,nmol,iso,iso_max)
            go to 50
         else if (iso .gt. iso_max(m)) then
            call line_exception (2,ipr,h_linf4,m,nmol,iso,iso_max)
            go to 50
         endif
C
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
c
         IF (JRAD.EQ.1) SUI = SUI*VNU(ij)
C                                                                         D01870
         IF (VNU(IJ).EQ.0.) SUI = 2.*SUI                                  D01880
C                                                                         D01890
C     TREAT TRANSITIONS WITH UNKNOWN EPP AS SPECIAL CASE                  D01900
C                                                                         D01910
c>>   an epp value between -0.9999 and 0.  cm-1 is taken as valid 
c
c>>   an epp value of -1. is assumed set by hitran indicating an unknown 
c     value: no temperature correction is performed
c
c>>   for an epp value of less than -1., it is assumed that value has
c     been provided as a reasonable value to be used for purposes of 
c     temperature correction.  epp is set positive
c
         if (epp(ij).le.-1.001)  epp(ij) = abs(epp(ij))

         if (epp(ij).le.-0.999)  MEFDP(M) = MEFDP(M)+1 

c     temperature correction:

         if (epp(ij) .gt. -0.999) then
            SUI = SUI*SCOR(m,iso)*  
     *              EXP(-EPP(ij)*BETACR)*(1.+EXP(-VNU(ij)*BETA))
         endif
C                                                                         D01980
         SUMS = SUMS+SUI                                                  D01990
C                                                                         D02000
C     TEMPERATURE CORRECTION OF THE HALFWIDTH                             D02010
C     SELF TEMP DEPENDENCE TAKEN THE SAME AS FOREIGN                      D02020
C                                                                         D02030
         TMPCOR = TRATIO**TMPALB(I)                                       D02040
         ALFA0I = ALFA0(IJ)*TMPCOR                                        D02050
         HWHMSI = HWHMB(I)*TMPCOR                                         D02060
         ALFAL(IJ) = ALFA0I*(RHORAT-RHOSLF(m))+HWHMSI*RHOSLF(m)
C                                                                         D02080
         IF (IFLAG.EQ.3)                                                  D02090
     *        ALFAL(IJ) = ALFAL(IJ)*(1.0-GAMMA1*PAVP0-GAMMA2*PAVP2)       D02100
C                                                                         D02110
         ALFAD(IJ) = VNU(IJ)*ALFD1(m,iso)                                  D02120
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

      integer *4 lnfil,lnfil4

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
   10 CALL BUFIN_sgl(Lnfl,LEOF,rdpnl(1),npnlhd) 
c
      IF (LEOF.EQ.0) GO TO 30                                             D02980
      IF (VMAX.LT.VLO) THEN                                               D02990
         CALL BUFIN_sgl(lnfl,LEOF,dum(1),1) 
         GO TO 10                                                         D03010
      ELSE                                                                D03020
         CALL BUFIN_sgl(Lnfl,LEOF,vlin(1),NWDS) 
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
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,         D04660
     *               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL                      D04670
      COMMON /R4SUB/ VLO,VHI,ILO,IST,IHI,LIMIN,LIMOUT,ILAST,DPTMN,        D04680
     *               DPTFC,ILIN4,ILIN4T                                   D04690
      COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)        D04700
      COMMON /VOICOM/ AVRAT(102),DUMMY(5,102)                             D04710
      COMMON /CONVF/ CHI(0:250),RDVCHI,RECPI,ZSQBND,A3,B3,JCNVF4            D04720
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
      V2R4 = V1R4+DVR4* REAL(NPTR4-1)                                     D04970
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
         VI = V1R4+DVR4* REAL(NPTSI1-1)                                   D05430
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
      SUBROUTINE CONVF4 (VNU,SP,ALFAL,ALFAD,MOL,SPP)
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
      COMMON /CONVF/ CHI(0:250),RDVCHI,RECPI,ZSQBND,A3,B3,JCNVF4            D06380
C                                                                         D06390
      parameter (nzeta=101)
      real*8 a_1,a_2,a_3,b_1,b_2,b_3
      common /voigt_cf/
     $     a_1(0:nzeta-1), a_2(0:nzeta-1), a_3(0:nzeta-1),
     $     b_1(0:nzeta-1), b_2(0:nzeta-1), b_3(0:nzeta-1)
C 
      DIMENSION VNU(*),SP(*),ALFAL(*),ALFAD(*),MOL(*),SPP(*)              D06400
C                                                                         D06410
      DATA ZBND / 64. /                                                   D06420
      DATA HREJ /'0'/,HNOREJ /'1'/
C                                                                         D06440
      DATA I_1/1/, I_250/250/
C
      VNULST = V2R4+BOUND4                                                D06450
C                                                                         D06460
      IF (JCNVF4.NE.1234) then
c
         JCNVF4 = 1234       
c
C        Obtain CHI SUB-LORENTZIAN FORM FACTOR FOR CARBON DIOXIDE: 
c
         dvchi  = 0.1
         rdvchi = 1./dvchi

         call chi_fn (chi,dvchi)
C                                                                            D06760
C        CONSTANTS FOR FOURTH FUNCTION LINE SHAPE:                           D06770
C                                                                            D06780
         RECPI = 1./(2.*ASIN(1.))                                            D06790
         ZSQBND = ZBND*ZBND                                                  D06800
         A3 = (1.+2.*ZSQBND)/(1.+ZSQBND)**2                                  D06810
         B3 = -1./(1.+ZSQBND)**2                                             D06820
C                                                                            D06830
      endif
C                                                                         D06850
      BNDSQ = BOUND4*BOUND4                                               D06860
C                                                                         D06870
C     START OF LOOP OVER LINES                                            D06880
C                                                                         D06890
      IF (ILNFLG.EQ.2) READ(16)(FREJ(I),I=ILO,IHI)
C
      DO 60 I = ILO, IHI                                                  D06900
C                                                                         D06910
         IF (SP(I).EQ.0..AND.SPP(I).EQ.0.) GO TO 60                        D06920
         ALFADI = ALFAD(I)                                                D06930
         ALFALI = ALFAL(I)                                                D06940
         ZETAI = ALFALI/(ALFALI+ALFADI)                                   D06950
         IZ = 100.*ZETAI + ONEPL                                          D06960
         ZETDIF = 100.*ZETAI -  REAL(IZ-1)
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
         SIV = SP(I)*RALFVI
         siv_a = siv*a3x
         siv_b = siv*b3x
c
         SPEAK = A3x*(ABS(SIV))
C                                                                         D07090
         JJ = (VNU(I)-V1R4)/DVR4+1.                                       D07100
         JJ = MAX(JJ,I_1)                                                   D07110
         JJ = MIN(JJ,NPTR4)                                               D07120
C
C     SPEAK is used for line rejection
         IF (ILNFLG.LE.1) THEN
            FREJ(I) = HNOREJ
c     No rejection for line-coupled lines (SPP ne. 0)
            IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ)) .and. spp(i).eq.0.) THEN
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
         JMIN = MAX(JMIN,I_1)                                               D07250
         JMAX = (XNUI+BOUND4)/DVR4+1.                                     D07260
         IF (JMAX.LT.JMIN) GO TO 50                                       D07270
         JMAX = MIN(JMAX,NPTR4)                                           D07280
         ALFLI2 = ALFALI*ALFALI                                           D07290
         ALFVI2 = ALFAVI*ALFAVI                                           D07300
         XJJ =  REAL(JMIN-1)*DVR4                                         D07310
c
         siv_64 = siv*f4_64x*(alfli2 + alfvi2*zsqbnd)
         f4bnd  = siv_64/(ALFLI2+bndSQ)
C                                                                         D07340
C                FOURTH FUNCTION CONVOLUTION                              D07350
C                                                                         D07360

         dptrat = spp(i)/(sp(i)*alfavi)
         rec_alfvi2 = 1./ALFVI2
         siv_a3 = SIV*A3
         siv_b3 = SIV*B3
c
         DO 40 JJ = JMIN, JMAX                                            D07370
            XM = (XJJ-XNUI)                                               D07380
            XMSQ = XM*XM                                                  D07390
            ZVSQ = XMSQ * rec_alfvi2  
C                                                                         D07410
            IF (ZVSQ.LE.ZSQBND) THEN                                      D07420
               F4FN = (siv_A + ZVSQ * siv_B) - F4BND
               IF (SPP(I).NE.0.) then
                   f4fn = f4fn + xm*dptrat*f4fn
               endif 
            ELSE                                                          D07460
                f4fn = siv_64/(alfli2+xmsq) - f4bnd
                if (SPP(I).NE.0.) then
                    F4FN = F4FN + XM*dptrat*f4fn
                endif
            ENDIF                                                         D07500
C                                                                         D07510
            IF (MOL(I).EQ.2.AND.SPP(I).EQ.0.) THEN                        D07520
C                                                                         D07530
C     ASSIGN ARGUMENT ISUBL OF THE FORM FACTOR FOR CO2 LINES              D07540
C                                                                         D07550
               ISUBL = RDVCHI*ABS(XM)+0.5                                 D07560
               ISUBL = MIN(ISUBL,i_250)                                     D07570
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
C
            SPP(I) = -SPP(I)
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
c____________________________________________________________________
c
      subroutine chi_fn (chi,dvchi)
c
      dimension chi(0:250)
c
      DATA ASUBL / 0.800 /,BSUBL / 10.0 /
C                                                                         D06490
C     SET UP CHI SUB-LORENTZIAN FORM FACTOR FOR CARBON DIOXIDE            D06500
C     POLYNOMIAL MATCHED TO AN EXPONENTIAL AT X0 = 10 CM-1                 D06510
C                                                                         D06520
c     0 - 25 cm-1 on 0.1 cm-1 grid

      X0 = 10.                                                             D06530
      Y0 = asubl*EXP(-x0/bsubl) 
      F  = 1./bsubl
      Y1 = -F*Y0                                                          D06560
      Y2 = Y1*((BSUBL-1)/X0-F)                                            D06570
      Z0 = (Y0-1)/X0**2                                                   D06580
      Z1 = Y1/(2*X0)                                                      D06590
      Z2 = Y2/2.                                                          D06600
      C6 = (Z0-Z1+(Z2-Z1)/4.)/X0**4                                       D06610
      C4 = (Z1-Z0)/X0**2-2.*X0**2*C6                                      D06620
      C2 = Z0-X0**2*C4-X0**4*C6                                           D06630
C
      DO 10 ISUBL = 0, 250                                                D06670
         FI = DVCHI* REAL(ISUBL)                                          D06680
         IF (FI.LT.X0) THEN                                               D06690
            CHI(ISUBL) = 1.+C2*FI**2+C4*FI**4+C6*FI**6                    D06700
         ELSE                                                             D06710
            CHI(ISUBL) = asubl*EXP(-fi/bsubl) 
         ENDIF                                                            D06740
c
   10 CONTINUE                                                            D06750
c
      return
c
      end
c____________________________________________________________________
c
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
     *           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=38,MXTRAC=22)
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
      COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(38),XAMNT(38,MXLAY)              E00230
C                                                                         E00240
C     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES         E00250
C     FOR THE CROSS-SECTION MOLECULES.                                    E00260
C     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES          E00270
C                                                                         E00280
c%%%%%LINUX_PGI90 (-i8)%%%%%      integer*4 iostat
      CHARACTER*10 XSFILE,XSNAME,ALIAS,XNAME,XFILS(6),BLANK               E00290
      COMMON /XSECTF/ XSFILE(6,5,38),XSNAME(38),ALIAS(4,38)               E00300
      COMMON /XSECTR/ V1FX(5,38),V2FX(5,38),DVFX(5,38),WXM(38),           E00310
     *                NTEMPF(5,38),NSPECR(38),IXFORM(5,38),               E00320
     *                XSMASS(38),XDOPLR(5,38),NUMXS,IXSBIN                E00325
C                                                                         E00330
      DIMENSION IXFLG(38)                                                 E00340
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
      IXMAX = 38                                                          E00430
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
CPRT 1    (5X,I5,2X,A1))') (I,XSNAME(I),I=1,IXMOLS)                        E00620
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
      OPEN (IXFIL,FILE='FSCDXS',STATUS='OLD',FORM='FORMATTED',
     *       IOSTAT=iostat)
        if (IOSTAT.gt.0) stop 'FSCDXS does not exist - oprop.f'
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
C     3.58115E-07 = SQRT( 2.* LOG(2.)*AVOGAD*BOLTZ/(CLIGHT*CLIGHT) ) 
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
      COMMON /XSECTI/ XSMAX(6,5,38),XSTEMP(6,5,38),NPTSFX(5,38),          E02850
     *                NFILEX(5,38),NLIMX                                  E02860
      COMMON /XSECTF/ XSFILE(6,5,38),XSNAME(38),ALIAS(4,38)               E01710
      COMMON /XSECTR/ V1FX(5,38),V2FX(5,38),DVFX(5,38),WXM(38),           E01720
     *                NTEMPF(5,38),NSPECR(38),IXFORM(5,38),               E01730
     *                XSMASS(38),XDOPLR(5,38),NUMXS,IXSBIN                E01740
      COMMON /XSECTS/ JINPUT,NMODES,NPANEL,NDUM,V1XS,V2XS,DVXS,NPTSXS     E02870
C                                                                         E01750
      DATA NMODES / 1 /,NPANEL / 0 /,V1XS / 0.0 /,V2XS / 0.0 /,           E02990
     *     DVXS / 0.0 /,NPTSXS / 0 /                                      E03000
      DATA XSMAX / 1140*0.0 /                                             E03010
      DATA (ALIAS(1,I),I=1,38)/                                           E01760
     *    'CLONO2    ', 'HNO4      ', 'CHCL2F    ', 'CCL4      ',         E01770
     *    'CCL3F     ', 'CCL2F2    ', 'C2CL2F4   ', 'C2CL3F3   ',         E01780
     *    'N2O5      ', 'HNO3      ', 'CF4       ', 'CHCLF2    ',         E01790
     *    'CCLF3     ', 'C2CLF5    ', 24*' ZZZZZZZZ ' /                   E01800
      DATA (ALIAS(2,I),I=1,38)/                                           E01810
     *    'CLNO3     ', ' ZZZZZZZZ ', 'CFC21     ', ' ZZZZZZZZ ',         E01820
     *    'CFCL3     ', 'CF2CL2    ', 'C2F4CL2   ', 'C2F3CL3   ',         E01830
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CHF2CL    ',         E01840
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 24*' ZZZZZZZZ ' /                   E01850
      DATA (ALIAS(3,I),I=1,38)/                                           E01860
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CFC21     ', ' ZZZZZZZZ ',         E01870
     *    'CFC11     ', 'CFC12     ', 'CFC114    ', 'CFC113    ',         E01880
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'CFC14     ', 'CFC22     ',         E01890
     *    'CFC13     ', 'CFC115    ', 24*' ZZZZZZZZ ' /                   E01900
      DATA (ALIAS(4,I),I=1,38)/                                           E01910
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'F21       ', ' ZZZZZZZZ ',         E01920
     *    'F11       ', 'F12       ', 'F114      ', 'F113      ',         E01930
     *    ' ZZZZZZZZ ', ' ZZZZZZZZ ', 'F14       ', 'F22       ',         E01940
     *    'F13       ', 'F115      ', 24*' ZZZZZZZZ ' /                   E01950
C                                                                         E01960
C     XSMASS IS MASS OF EACH CROSS-SECTION                                E01961
C                                                                         E01962
      DATA XSMASS/                                                        E01963
     1      97.46     ,   79.01     ,  102.92     ,  153.82     ,         E01964
     2     137.37     ,  120.91     ,  170.92     ,  187.38     ,         E01965
     3     108.01     ,   63.01     ,   88.00     ,   86.47     ,         E01966
     4     104.46     ,  154.47     ,  24*0.00 /                          E01967
C                                                                         E01969
      DATA V1FX / 190*0.0 /,V2FX / 190*0.0 /,DVFX / 190*0.0 /,            E01970
     *     WXM / 38*0.0 /                                                 E01980
      DATA NTEMPF / 190*0 /,NSPECR / 38*0 /,IXFORM / 190*0 /,             E01990
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
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
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
      COMMON /XSHEAD/ HEADT1(38)                                          E02770
      COMMON /XSTMPR/ PF,TF,PDX(6,5,38),DVXPR(5,38),IXBIN(5,38),          E02780
     *                IXSBN(5,38)                                         E02790
      COMMON /XSECTP/ V1X,V2X,DVX,NPTSX,RX(13000)                         E02800
      COMMON /XSECTD/ V1DX,V2DX,DVDX,NPTSDX,RDX1(520),RDX2(520)           E02810
      COMMON /XSECTF/ XSFILE(6,5,38),XSNAME(38),ALIAS(4,38)               E02820
      COMMON /XSECTR/ V1FX(5,38),V2FX(5,38),DVFX(5,38),WXM(38),           E02830
     *                NTEMPF(5,38),NSPECR(38),IXFORM(5,38),               E02840
     *                XSMASS(38),XDOPLR(5,38),NUMXS,IXSBIN                E02845
      COMMON /XSECTI/ XSMAX(6,5,38),XSTEMP(6,5,38),NPTSFX(5,38),          E02850
     *                NFILEX(5,38),NLIMX                                  E02860
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
      DATA I_0/0/
C
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
   20 VFX2 = VFT+2.*DVX+ REAL(NHI)*DV                                     E03540
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
            V2X = V1X+ REAL(LIMOUT-1)*DVX                                 E03700
            IF (V2X.GT.V2) V2X = V1X+ REAL(INT((V2-V1X)/DVX)+3)*DVX       E03710
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
            V2X = V1X+ REAL(LIMOUT-1)*DVX                                 E04340
            IF (V2X.GT.V2) V2X = V1X+ REAL(INT((V2-V1X)/DVX)+2)*DVX       E04350
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
               NSKIP = MAX(NSKIP,I_0)                                       E04890
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
                  NRSKIP = MAX(NRSKIP,I_0)                                  E05000
               ENDIF                                                      E05010
               V1FP = V1FX(NS,NI)+ REAL(NRSKIP)*DVFXX                     E05020
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
               N2RX = MAX(N2RX,I_0)                                         E05130
C                                                                         E05140
C     IMAX = -4 TO PLACE THE FIRST PANEL V1 AT ARRAY LOCATION 1           E05150
C                                                                         E05160
               IMAX = -4                                                  E05170
               DO 100 NP = 1, NPAN                                        E05180
                  V1FP = V1FP+ REAL(IMAX)*DVFXX                           E05190
                  IMAX = NMAX-(NP-1)*NLIMX                                E05200
                  IF (IAFORM.GT.100.AND.NPANEL.LE.0.AND.                  E05210
     *                NBSKIP.EQ.0.AND.IMAX.GT.500) IMAX = 500             E05220
                  IMAX = MIN(IMAX,NLIMX)                                  E05230
C                                                                         E05240
C     FOR V2FP IMAX + 3 GIVES US ARRAY LOCATION 514                       E05250
C             (504 FOR FIRST PANEL OF BLOCKED DATA)                       E05260
C                                                                         E05270
                  V2FP = V1FP+ REAL(IMAX+3)*DVFXX                         E05280
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
            VI = V1X+ REAL(NPTSI1-1)*DVX                                  E06080
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
      COMMON /XSECTF/ XSFILE(6,5,38),XSNAME(38),ALIAS(4,38)               E06900
      COMMON /XSECTR/ V1FX(5,38),V2FX(5,38),DVFX(5,38),WXM(38),           E06910
     *                NTEMPF(5,38),NSPECR(38),IXFORM(5,38),               E06920
     *                XSMASS(38),XDOPLR(5,38),NUMXS,IXSBIN                E06925
      COMMON /XSECTI/ XSMAX(6,5,38),XSTEMP(6,5,38),NPTSFX(5,38),          E06930
     *                NFILEX(5,38),NLIMX                                  E06940
      COMMON /XSHEAD/ HEADT1(38)                                          E06950
      COMMON /XSTMPR/ PF,TF,PDX(6,5,38),DVXPR(5,38),IXBIN(5,38),          E06960
     *                IXSBN(5,38)                                         E06970
      COMMON /FLFORM/ CFORM                                               E06980
C                                                                         E06990
      CHARACTER*10 XSFILE,XSNAME,ALIAS,SOURCE(3),CTORR                    E07000
      CHARACTER AMOL*8,BMOL*6,HEADER*100,HEADT1*100                       E07010
      CHARACTER XSFIL1*10,XSFIL2*10,XSTMP*4,XSNUM*3,CI*1                  E07020
      CHARACTER CFORM*11,BFRM*10,UNBFRM*10,BLKFRM*10,BFORM*9              E07030
      LOGICAL OP,OPCL                                                     E07040
c%%%%%LINUX_PGI90 (-i8)%%%%%      integer*4 iostat
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
C
      DATA I_100/100/
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
      IMFORM = MOD(IAFORM,I_100)                                            E07330
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
               OPEN (IFILE,FILE=XSFIL1,STATUS='OLD',FORM=BFORM,
     *          IOSTAT=iostat)
               if (IOSTAT.gt.0) stop 'in oprop - No file XSFIL1' 
            ELSE                                                          E07510
               OPEN (IFILE,FILE=XSFIL1,STATUS='OLD',FORM=CFORM,
     *          IOSTAT=iostat)
               if (IOSTAT.gt.0) stop 'in oprop - No file XSFIL1 ' 
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
                  OPEN (JFILE,FILE=XSFIL2,STATUS='OLD',FORM=BFORM,
     *          IOSTAT=iostat)
                  if (IOSTAT.gt.0) stop 'in oprop - No file XSFIL2 ' 
               ELSE                                                       E07650
                  OPEN (JFILE,FILE=XSFIL2,STATUS='OLD',FORM=CFORM,
     *          IOSTAT=iostat)
                  if (IOSTAT.gt.0) stop 'in oprop - No file XSFIL2 ' 
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
                  XSTEMP(NT1,NS,NI) =  REAL(ITEMP)+273.15                 E07900
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
                     XSTEMP(NT2,NS,NI) =  REAL(ITEMP)+273.15              E08090
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
                  XSTEMP(NT1,NS,NI) =  REAL(ITEMP)+273.15                 E08420
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
                     XSTEMP(NT2,NS,NI) =  REAL(ITEMP)+273.15              E08700
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
         DVDX = (V2DX-V1DX)/ REAL(NPTSDX-1)                               E08850
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
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /XSECTP/ V1X,V2X,DVX,NPTSX,RX(13000)                         E09770
      COMMON /XSECTD/ V1DX,V2DX,DVDX,NPTSDX,RDX1(520),RDX2(520)           E09780
      COMMON /XSECTF/ XSFILE(6,5,38),XSNAME(38),ALIAS(4,38)               E09790
      COMMON /XSECTR/ V1FX(5,38),V2FX(5,38),DVFX(5,38),WXM(38),           E09800
     *                NTEMPF(5,38),NSPECR(38),IXFORM(5,38),               E09810
     *                XSMASS(38),XDOPLR(5,38),NUMXS,IXSBIN                E09815
      COMMON /XSECTI/ XSMAX(6,5,38),XSTEMP(6,5,38),NPTSFX(5,38),          E09820
     *                NFILEX(5,38),NLIMX                                  E09830
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
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /XSECTP/ V1X,V2X,DVX,NPTSX,RX(13000)                         E10720
      COMMON /XSECTD/ V1DX,V2DX,DVDX,NPTSDX,RDX1(520),RDX2(520)           E10730
      COMMON /XSECTF/ XSFILE(6,5,38),XSNAME(38),ALIAS(4,38)               E10740
      COMMON /XSECTR/ V1FX(5,38),V2FX(5,38),DVFX(5,38),WXM(38),           E10750
     *                NTEMPF(5,38),NSPECR(38),IXFORM(5,38),               E10760
     *                XSMASS(38),XDOPLR(5,38),NUMXS,IXSBIN                E10765
      COMMON /XSECTI/ XSMAX(6,5,38),XSTEMP(6,5,38),NPTSFX(5,38),          E10770
     *                NFILEX(5,38),NLIMX                                  E10780
      COMMON /XSTMPR/ PF,TF,PDX(6,5,38),DVXPR(5,38),IXBIN(5,38),          E10790
     *                IXSBN(5,38)                                         E10800
      COMMON /XSHEAD/ HEADT1(38)                                          E10810
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
         V2PX = V1PX+ REAL(LPMAX+NMAX-1)*DVFX(NS,NI)                      E11760
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
      BOUND =  REAL(NBOUND)*DVO/2.                                        E13570
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
      COMMON /XSTMPR/ PF,TF,PDX(6,5,38),DVXPR(5,38),IXBIN(5,38),          E14340
     *                IXSBN(5,38)                                         E14350
      COMMON /XSHEAD/ HEADT1(38)                                          E14360
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
      V1P = VFT+ REAL(NLO-1)*DV                                           E14560
      V2P = VFT+ REAL(NHI-1)*DV                                           E14570
C                                                                         E14580
C     V1P IS FIRST FREQ OF PANEL                                          E14590
C     V2P IS LAST  FREQ OF PANEL                                          E14600
C                                                                         E14610
      CALL BUFOUT (JFILE,PNLHDR(1),NPHDRF)                                E14620
      CALL BUFOUT (JFILE,R1(NLO),NLIM)                                    E14630
C                                                                         E14640
      VFT = VFT+ REAL(NLIMF-1)*DV                                         E14650
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
         X =  REAL(I-1)*DXF                                               E15140
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
C***********************
      Subroutine TIPS_2003 (mol_max,iso_max,temp_lbl,scor)
C*********************** 
c
c    This program calculates the total internal
c    partition sum (TIPS) for a given molecule, isotopic species, and
c    temperature.  Current limitations are the molecular species on the
c    HITRAN molecular database and the temperature range 70 - 3000 K.
c
c...date last changed 18 March, 2003
c
C    Program TIPS_2003 written by R.R. Gamache
C 
c..  JQSRT - HITRAN Special Issue, 2003.
c..     J. Fischer(a) R.R. Gamache(a), A. Goldman(b), L.S. Rothman(c), 
c..     and A. Perrin(d)
c..
c..     (a) Department of Environmental, Earth, and Atmospheric Sciences,
c..         University of Massachusetts Lowell, Lowell, MA 01854, U.S.A.
c..
c..     (b) Department of Physics, 
c..         University of Denver, Denver, CO 80208, U.S.A.
c..
c..     (c) Harvard-Smithsonian Center for Astrophysics, 
c..         60 Garden St, Cambridge, MA 02138 USA
c..
c..     (d) Laboratoire de Photophysique Moleculaire, 
c..         Universite Paris Sud, 91405 Orsay, FRANCE
c..
c..  Corresponding author. Email address: Robert_Gamache@uml.edu
c..
c.. Abstract: Total internal partition sums (TIPS) are calculated for all 
c.. molecular species in the 2000 HITRAN database.  In addition, the TIPS for 13 
c.. other isotopomers/isotopologues of ozone and carbon dioxide are presented.  
c.. The calculations address the corrections suggested by Goldman et al. 
c.. (JQSRT 2000;66:455-86).  The calculations consider the temperature range 
c.. 70-3000 K to be applicable to a variety of remote sensing needs.  The method 
c.. of calculation for each molecular species is stated and comparisons with data
c.. from the literature are discussed.  A new method of recall for the partition 
c.. sums, Lagrange 4-point interpolation, is developed.  This method, unlike 
c.. previous versions of the TIPS code, allows all molecular species to be 
c.. considered.  
c_________________________________________________________________________________
c
c...This program calculates the TIPS by 4-point LaGrange interpolation
c
c++
c     Max_ISO here is the number that TIPS treats: LBLRTM is limited to 9
c
      PARAMETER (NMOL=38,Max_ISO=20)
c++
      PARAMETER (NT=119)
c++
      COMMON /ISO_data/ISO82(NMOL,Max_ISO),ISONM(NMOL)
c++:  bd-MOL
      CHARACTER*6 MOLID
c++:  bd-MOL
      COMMON/MOLNAM/MOLID(0:NMOL)
c++:  bd-QT
      COMMON/Temperatures/tdat(NT)

      dimension iso_max(nmol),scor(42,9)
 
      character*30 stopNgo
c
      data Tdat/  60.,  85., 110., 135., 160., 185., 210., 235.,
     + 260., 285., 310., 335., 360., 385., 410., 435., 460., 485.,
     + 510., 535., 560., 585., 610., 635., 660., 685., 710., 735.,
     + 760., 785., 810., 835., 860., 885., 910., 935., 960., 985.,
     +1010.,1035.,1060.,1085.,1110.,1135.,1160.,1185.,1210.,1235.,
     +1260.,1285.,1310.,1335.,1360.,1385.,1410.,1435.,1460.,1485.,
     +1510.,1535.,1560.,1585.,1610.,1635.,1660.,1685.,1710.,1735.,
     +1760.,1785.,1810.,1835.,1860.,1885.,1910.,1935.,1960.,1985.,
     +2010.,2035.,2060.,2085.,2110.,2135.,2160.,2185.,2210.,2235.,
     +2260.,2285.,2310.,2335.,2360.,2385.,2410.,2435.,2460.,2485.,
     +2510.,2535.,2560.,2585.,2610.,2635.,2660.,2685.,2710.,2735.,
     +2760.,2785.,2810.,2835.,2860.,2885.,2910.,2935.,2960.,2985.,
     +3010./
C
C
c       CALL ISO_82_TO_85 (MOL,NSO82,ISO,iso_max)
C

      do 110 mol= 1,mol_max
         do 110 iso= 1,iso_max(mol)

            do 105 itemp=1,2
               if (itemp .eq. 1) temp = 296.
               if (itemp .eq. 2) temp = temp_lbl

C     GO TO PARTICULAR MOLECULE:  
 
      IF(MOL.EQ.1) THEN 
      CALL QT_H2O(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.2) THEN 
      CALL QT_CO2(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.3) THEN 
      CALL QT_O3(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.4) THEN 
      CALL QT_N2O(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.5) THEN 
      CALL QT_CO(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.6) THEN 
      CALL QT_CH4(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.7) THEN 
      CALL QT_O2(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.8) THEN 
      CALL QT_NO(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.9) THEN 
      CALL QT_SO2(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.10) THEN
      CALL QT_NO2(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.11) THEN
      CALL QT_NH3(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.12) THEN
      CALL QT_HNO3(Temp,ISO,gi,QT)
      go to 100 
      ENDIF 
C 
      IF(MOL.EQ.13) THEN
      CALL QT_OH(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.14) THEN
      CALL QT_HF(Temp,ISO,gi,QT) 
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.15) THEN
      CALL QT_HCL(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.16) THEN
      CALL QT_HBR(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.17) THEN
      CALL QT_HI(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.18) THEN
      CALL QT_CLO(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.19) THEN
      CALL QT_OCS(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.20) THEN
      CALL QT_H2CO(Temp,ISO,gi,QT) 
      go to 100 
      ENDIF 
C 
      IF(MOL.EQ.21) THEN
      CALL QT_HOCL(Temp,ISO,gi,QT) 
      go to 100 
      ENDIF 
C 
      IF(MOL.EQ.22) THEN
      CALL QT_N2(Temp,ISO,gi,QT) 
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.23) THEN
      CALL QT_HCN(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
C 
      IF(MOL.EQ.24) THEN
      CALL QT_CH3CL(Temp,ISO,gi,QT)
      go to 100 
      ENDIF 
C 
      IF(MOL.EQ.25) THEN
      CALL QT_H2O2(Temp,ISO,gi,QT) 
      go to 100 
      ENDIF 
C 
      IF(MOL.EQ.26) THEN
      CALL QT_C2H2(Temp,ISO,gi,QT)
      go to 100 
      ENDIF 
C 
      IF(MOL.EQ.27) THEN
      CALL QT_C2H6(Temp,ISO,gi,QT)
      go to 100 
      ENDIF 
C 
      IF(MOL.EQ.28) THEN
      CALL QT_PH3(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
c
      IF(MOL.EQ.29) THEN 
      CALL QT_COF2(Temp,ISO,gi,QT)
      go to 100 
      ENDIF 
c
      IF(MOL.EQ.30) THEN 
      CALL QT_SF6(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
c
      IF(MOL.EQ.31) THEN
      CALL QT_H2S(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
c
      IF(MOL.EQ.32) THEN
      CALL QT_HCOOH(Temp,ISO,gi,QT)
      go to 100 
      ENDIF 
c
      IF(MOL.EQ.33) THEN 
      CALL QT_HO2(Temp,ISO,gi,QT)
      go to 100
      ENDIF 
c
      IF(MOL.EQ.34) THEN 
c...not applicable to O
      gi=0
      QT = 0.
      go to 100
      ENDIF 
c
      IF(MOL.EQ.35) THEN 
      CALL QT_ClONO2(Temp,ISO,gi,QT)
      go to 100 
      ENDIF 
c
      IF(MOL.EQ.36) THEN 
      CALL QT_NOp(Temp,ISO,gi,QT)
      go to 100
      ENDIF     
C 
      IF(MOL.EQ.37) THEN 
      CALL QT_HOBr(Temp,ISO,gi,QT)
      go to 100
      ENDIF     
C 
      IF(MOL.EQ.38) THEN 
      CALL QT_C2H4(Temp,ISO,gi,QT)
      go to 100
      ENDIF     
C 
 100  continue
c
      if (QT .le. 0.) stop ' partition sum less than 0.'
c
c      if(QT .gt. 0.) then
c        Write(*,900) MOLID(Mol), mol, iso, iso82(mol,iso),Temp, QT, gi
c  900 format(5x,A6, 2i4,i8,
c     * '   Q(',f7.1,' K) = ',1pE12.6,5x,0pf7.1)
c      endif

      if (itemp .eq.1) qt_296 = qt
      if (itemp .eq.2) qt_temp = qt

 105  continue

      scor(mol,iso) = qt_296/qt_temp

 110  continue
c
      END 
  
  
C**************************************
      SUBROUTINE CLEAR
C**************************************
C
      DO 10 I = 1,24
       WRITE (*,'(1X)')
   10 CONTINUE
      RETURN
      END
  
c**************************************
      BLOCK DATA BDMol
c**************************************
c
c      March 20, 1990
c...date last changed 19 February, 2002
c
c++
      PARAMETER (NMOL=38)
c++:  bd-MOL
      CHARACTER*6 MOLID
c++:  bd-MOL
      COMMON/MOLNAM/MOLID(0:NMOL)
c
      DATA (MOLID(I),I=0,NMOL)/ '   All',
     1'   H2O','   CO2','    O3','   N2O','    CO','   CH4','    O2',
     2'    NO','   SO2','   NO2','   NH3','  HNO3','    OH','    HF',
     3'   HCl','   HBr','    HI','   ClO','   OCS','  H2CO','  HOCl',
     4'    N2','   HCN',' CH3Cl','  H2O2','  C2H2','  C2H6','   PH3',
     5'  COF2','   SF6','   H2S',' HCOOH','   HO2','     O','ClONO2',
     6'   NO+','  HOBr','  C2H4'/
c
      END

  
c  ****************************************
      Block Data ISO_2002
c  ****************************************
c
c++
      PARAMETER (NMOL=38,Max_ISO=20)
c++
      COMMON /ISO_data/ ISO82(NMOL,Max_ISO),ISONM(NMOL)
c
c    The number of isotopes for a particular molecule:
      DATA (ISONM(I),I=1,NMOL)/
c     H2O, CO2, O3, N2O, CO, CH4, O2,
     +  6,   9,  18,   5,  6,   3,  3,
c      NO, SO2, NO2, NH3, HNO3, OH, HF, HCl, HBr, HI,
     +  3,   2,   1,   2,    1,  3,  1,   2,   2,  1,
c     ClO, OCS, H2CO, HOCl, N2, HCN, CH3Cl, H2O2, C2H2, C2H6, PH3
     +  2,   5,    3,    2,  1,   3,     2,    1,    2,    1,   1,       4/24/97
c     COF2, SF6, H2S, HCOOH, HO2, O, ClONO2,  NO+, HOBr, C2H4
     +   1,   1,   3,     1,   1, 1,      2,    1,    2,    2/
c

c       H2O
      DATA (ISO82(1,J),J=1,6)/
     +  161,181,171,162,182,172/

      DATA (ISO82(2,J),J=1,9)/
c       CO2
     +  626,636,628,627,638,637,828,728,727/
c       O3
      DATA (ISO82(3,J),J=1,18)/
     +  666,668,686,667,676,886,868,678,768,
     +  786,776,767,888,887,878,778,787,777/

      DATA (ISO82(4,J),J=1,5)/
c       N2O
     +  446,456,546,448,447/

      DATA (ISO82(5,J),J=1,6)/
c       CO,                 
     +  26,36,28,27,38,37/  

      DATA (ISO82(6,J),J=1,3)/
c       CH4
     +  211,311,212/

      DATA (ISO82(7,J),J=1,3)/
c       O2        
     +  66,68,67/  

      DATA (ISO82(8,J),J=1,3)/
c       NO
     +  46,56,48/

      DATA (ISO82(9,J),J=1,2)/
c       SO2
     +  626,646/

      DATA (ISO82(10,J),J=1,1)/
c      NO2 
     + 646/  

      DATA (ISO82(11,J),J=1,2)/
c       NH3
     +  4111,5111/

      DATA (ISO82(12,J),J=1,1)/
c       HNO3
     +  146/

      DATA (ISO82(13,J),J=1,3)/
c       OH
     +  61,81,62/

      DATA (ISO82(14,J),J=1,1)/
c       HF
     +  19/

      DATA (ISO82(15,J),J=1,2)/
c       HCl
     +  15,17/

      DATA (ISO82(16,J),J=1,2)/
c       HBr
     +  19,11/

      DATA (ISO82(17,J),J=1,1)/
c       HI
     +  17/

      DATA (ISO82(18,J),J=1,2)/
c       ClO,  
     +  56,76/

      DATA (ISO82(19,J),J=1,5)/
c       OCS              
     +  622,624,632,623,822/

      DATA (ISO82(20,J),J=1,3)/
c       H2CO
     +  126,136,128/

      DATA (ISO82(21,J),J=1,2)/
c       HOCl,    
     +  165,167/

      DATA (ISO82(22,J),J=1,1)/
c       N2
     +  44/

      DATA (ISO82(23,J),J=1,3)/
c       HCN
     +  124,134,125/

      DATA (ISO82(24,J),J=1,2)/
c       CH3Cl
     +  215,217/

      DATA (ISO82(25,J),J=1,1)/
c       H2O2
     +  1661/

      DATA (ISO82(26,J),J=1,2)/
c       C2H2       
     +  1221,1231/

      DATA (ISO82(27,J),J=1,1)/
c       C2H6
     +  1221/ 

      DATA (ISO82(28,J),J=1,1)/
c       PH3
     =  1111/

      DATA (ISO82(29,J),J=1,1)/
c     COF2
     + 269/ 

      DATA (ISO82(30,J),J=1,1)/
c       SF6
     +  29/

      DATA (ISO82(31,J),J=1,3)/
c       H2S         
     +  121,141,131/

      DATA (ISO82(32,J),J=1,1)/
c       HCOOH 
     +  126/

      DATA (ISO82(33,J),J=1,1)/
c       HO2
     +  166/

      DATA (ISO82(34,J),J=1,1)/
c       O  
     +  6/

      DATA (ISO82(35,J),J=1,2)/
c       ClONO2     
     +  5646,7646/

      DATA (ISO82(36,J),J=1,1)/
c       NO+
     +  46/

      DATA (ISO82(37,J),J=1,2)/
c      HOBr
     + 169,161/

      DATA (ISO82(38,J),J=1,2)/
c       C2H4
     +  221,231/

      end

c
c...date last changed 25 Nov 2001
c  ****************************************
      SUBROUTINE ISO_82_to_85 (
     I  MOLEC, 
     O NSO82, ISO85, iso_max)
c  ****************************************
c
c
c++
      PARAMETER (NMOL=38,Max_ISO=20)
c++
      COMMON /ISO_data/ ISO82(NMOL,Max_ISO),ISONM(NMOL)
c
c
         iso_max = ISONM(MOLEC)
c
      RETURN
      END  
c
c     *****************
      Subroutine QT_H2O   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 6), QofT( 6,119),Q(NT)
      data xgj/ 1.,1.,6.,6.,6.,36/
c...      H2O
c...        --       161
      data (QofT( 1,J),J=1,119)/ 0.16824E+02, 0.27771E+02, 0.40408E+02,
     + 0.54549E+02, 0.70054E+02, 0.86817E+02, 0.10475E+03, 0.12380E+03,
     + 0.14391E+03, 0.16503E+03, 0.18714E+03, 0.21021E+03, 0.23425E+03,
     + 0.25924E+03, 0.28518E+03, 0.31209E+03, 0.33997E+03, 0.36883E+03,
     + 0.39870E+03, 0.42959E+03, 0.46152E+03, 0.49452E+03, 0.52860E+03,
     + 0.56380E+03, 0.60015E+03, 0.63766E+03, 0.67637E+03, 0.71631E+03,
     + 0.75750E+03, 0.79999E+03, 0.84380E+03, 0.88897E+03, 0.93553E+03,
     + 0.98353E+03, 0.10330E+04, 0.10840E+04, 0.11365E+04, 0.11906E+04,
     + 0.12463E+04, 0.13037E+04, 0.13628E+04, 0.14237E+04, 0.14863E+04,
     + 0.15509E+04, 0.16173E+04, 0.16856E+04, 0.17559E+04, 0.18283E+04,
     + 0.19028E+04, 0.19793E+04, 0.20581E+04, 0.21391E+04, 0.22224E+04,
     + 0.23080E+04, 0.24067E+04, 0.24975E+04, 0.25908E+04, 0.26867E+04,
     + 0.27853E+04, 0.28865E+04, 0.29904E+04, 0.30972E+04, 0.32068E+04,
     + 0.33194E+04, 0.34349E+04, 0.35535E+04, 0.36752E+04, 0.38001E+04,
     + 0.39282E+04, 0.40597E+04, 0.41945E+04, 0.43327E+04, 0.44745E+04,
     + 0.46199E+04, 0.47688E+04, 0.49215E+04, 0.50780E+04, 0.52384E+04,
     + 0.54027E+04, 0.55710E+04, 0.57434E+04, 0.59200E+04, 0.61008E+04,
     + 0.62859E+04, 0.64754E+04, 0.66693E+04, 0.68679E+04, 0.70710E+04,
     + 0.72788E+04, 0.74915E+04, 0.77090E+04, 0.79315E+04, 0.81590E+04,
     + 0.83917E+04, 0.86296E+04, 0.88728E+04, 0.91214E+04, 0.93755E+04,
     + 0.96351E+04, 0.99005E+04, 0.10171E+05, 0.10448E+05, 0.10731E+05,
     + 0.11020E+05, 0.11315E+05, 0.11617E+05, 0.11924E+05, 0.12238E+05,
     + 0.12559E+05, 0.12886E+05, 0.13220E+05, 0.13561E+05, 0.13909E+05,
     + 0.14263E+05, 0.14625E+05, 0.14995E+05, 0.15371E+05, 0.15755E+05,
     + 0.16147E+05/
c...        --       181
      data (QofT( 2,J),J=1,119)/ 0.15960E+02, 0.26999E+02, 0.39743E+02,
     + 0.54003E+02, 0.69639E+02, 0.86543E+02, 0.10463E+03, 0.12384E+03,
     + 0.14412E+03, 0.16542E+03, 0.18773E+03, 0.21103E+03, 0.23531E+03,
     + 0.26057E+03, 0.28681E+03, 0.31406E+03, 0.34226E+03, 0.37130E+03,
     + 0.40135E+03, 0.43243E+03, 0.46456E+03, 0.49777E+03, 0.53206E+03,
     + 0.56748E+03, 0.60405E+03, 0.64179E+03, 0.68074E+03, 0.72093E+03,
     + 0.76238E+03, 0.80513E+03, 0.84922E+03, 0.89467E+03, 0.94152E+03,
     + 0.98982E+03, 0.10396E+04, 0.10909E+04, 0.11437E+04, 0.11982E+04,
     + 0.12543E+04, 0.13120E+04, 0.13715E+04, 0.14328E+04, 0.14959E+04,
     + 0.15608E+04, 0.16276E+04, 0.16964E+04, 0.17672E+04, 0.18401E+04,
     + 0.19151E+04, 0.19922E+04, 0.20715E+04, 0.21531E+04, 0.22370E+04,
     + 0.23232E+04, 0.24118E+04, 0.25030E+04, 0.25967E+04, 0.26929E+04,
     + 0.27918E+04, 0.28934E+04, 0.29978E+04, 0.31050E+04, 0.32151E+04,
     + 0.33281E+04, 0.34441E+04, 0.35632E+04, 0.36854E+04, 0.38108E+04,
     + 0.39395E+04, 0.40715E+04, 0.42070E+04, 0.43459E+04, 0.44883E+04,
     + 0.46343E+04, 0.47840E+04, 0.49374E+04, 0.50946E+04, 0.52558E+04,
     + 0.54209E+04, 0.55900E+04, 0.57632E+04, 0.59407E+04, 0.61224E+04,
     + 0.63084E+04, 0.64988E+04, 0.66938E+04, 0.68933E+04, 0.70975E+04,
     + 0.73064E+04, 0.75202E+04, 0.77389E+04, 0.79625E+04, 0.81913E+04,
     + 0.84252E+04, 0.86644E+04, 0.89089E+04, 0.91588E+04, 0.94143E+04,
     + 0.96754E+04, 0.99422E+04, 0.10215E+05, 0.10493E+05, 0.10778E+05,
     + 0.11068E+05, 0.11365E+05, 0.11668E+05, 0.11977E+05, 0.12293E+05,
     + 0.12616E+05, 0.12945E+05, 0.13281E+05, 0.13624E+05, 0.13973E+05,
     + 0.14330E+05, 0.14694E+05, 0.15066E+05, 0.15445E+05, 0.15831E+05,
     + 0.16225E+05/
c...        --       171
      data (QofT( 3,J),J=1,119)/ 0.95371E+02, 0.16134E+03, 0.23750E+03,
     + 0.32273E+03, 0.41617E+03, 0.51722E+03, 0.62540E+03, 0.74036E+03,
     + 0.86185E+03, 0.98970E+03, 0.11238E+04, 0.12642E+04, 0.14097E+04,
     + 0.15599E+04, 0.17159E+04, 0.18777E+04, 0.20453E+04, 0.22188E+04,
     + 0.23983E+04, 0.25840E+04, 0.27760E+04, 0.29743E+04, 0.31792E+04,
     + 0.33907E+04, 0.36091E+04, 0.38346E+04, 0.40672E+04, 0.43072E+04,
     + 0.45547E+04, 0.48100E+04, 0.50732E+04, 0.53446E+04, 0.56244E+04,
     + 0.59128E+04, 0.62100E+04, 0.65162E+04, 0.68317E+04, 0.71567E+04,
     + 0.74915E+04, 0.78363E+04, 0.81914E+04, 0.85571E+04, 0.89335E+04,
     + 0.93211E+04, 0.97200E+04, 0.10131E+05, 0.10553E+05, 0.10988E+05,
     + 0.11435E+05, 0.11895E+05, 0.12368E+05, 0.12855E+05, 0.13356E+05,
     + 0.13870E+05, 0.14399E+05, 0.14943E+05, 0.15502E+05, 0.16076E+05,
     + 0.16666E+05, 0.17272E+05, 0.17895E+05, 0.18534E+05, 0.19191E+05,
     + 0.19865E+05, 0.20557E+05, 0.21267E+05, 0.21996E+05, 0.22744E+05,
     + 0.23512E+05, 0.24299E+05, 0.25106E+05, 0.25935E+05, 0.26784E+05,
     + 0.27655E+05, 0.28547E+05, 0.29462E+05, 0.30400E+05, 0.31361E+05,
     + 0.32345E+05, 0.33353E+05, 0.34386E+05, 0.35444E+05, 0.36527E+05,
     + 0.37637E+05, 0.38772E+05, 0.39934E+05, 0.41124E+05, 0.42341E+05,
     + 0.43587E+05, 0.44861E+05, 0.46165E+05, 0.47498E+05, 0.48862E+05,
     + 0.50256E+05, 0.51682E+05, 0.53139E+05, 0.54629E+05, 0.56152E+05,
     + 0.57708E+05, 0.59299E+05, 0.60923E+05, 0.62583E+05, 0.64279E+05,
     + 0.66011E+05, 0.67779E+05, 0.69585E+05, 0.71429E+05, 0.73312E+05,
     + 0.75234E+05, 0.77195E+05, 0.79197E+05, 0.81240E+05, 0.83325E+05,
     + 0.85452E+05, 0.87622E+05, 0.89835E+05, 0.92093E+05, 0.94395E+05,
     + 0.96743E+05/
c...        --       162
      data (QofT( 4,J),J=1,119)/ 0.75792E+02, 0.12986E+03, 0.19244E+03,
     + 0.26253E+03, 0.33942E+03, 0.42259E+03, 0.51161E+03, 0.60619E+03,
     + 0.70609E+03, 0.81117E+03, 0.92132E+03, 0.10365E+04, 0.11567E+04,
     + 0.12820E+04, 0.14124E+04, 0.15481E+04, 0.16891E+04, 0.18355E+04,
     + 0.19876E+04, 0.21455E+04, 0.23092E+04, 0.24791E+04, 0.26551E+04,
     + 0.28376E+04, 0.30268E+04, 0.32258E+04, 0.34288E+04, 0.36392E+04,
     + 0.38571E+04, 0.40828E+04, 0.43165E+04, 0.45584E+04, 0.48089E+04,
     + 0.50681E+04, 0.53363E+04, 0.56139E+04, 0.59009E+04, 0.61979E+04,
     + 0.65049E+04, 0.68224E+04, 0.71506E+04, 0.74898E+04, 0.78403E+04,
     + 0.82024E+04, 0.85765E+04, 0.89628E+04, 0.93618E+04, 0.97736E+04,
     + 0.10199E+05, 0.10637E+05, 0.11090E+05, 0.11557E+05, 0.12039E+05,
     + 0.12535E+05, 0.13047E+05, 0.13575E+05, 0.14119E+05, 0.14679E+05,
     + 0.15257E+05, 0.15851E+05, 0.16464E+05, 0.17094E+05, 0.17743E+05,
     + 0.18411E+05, 0.19098E+05, 0.19805E+05, 0.20532E+05, 0.21280E+05,
     + 0.22049E+05, 0.22840E+05, 0.23652E+05, 0.24487E+05, 0.25345E+05,
     + 0.26227E+05, 0.27132E+05, 0.28062E+05, 0.29016E+05, 0.29997E+05,
     + 0.31002E+05, 0.32035E+05, 0.33094E+05, 0.34180E+05, 0.35295E+05,
     + 0.36438E+05, 0.37610E+05, 0.38812E+05, 0.40044E+05, 0.41306E+05,
     + 0.42600E+05, 0.43926E+05, 0.45284E+05, 0.46675E+05, 0.48100E+05,
     + 0.49559E+05, 0.51053E+05, 0.52583E+05, 0.54148E+05, 0.55750E+05,
     + 0.57390E+05, 0.59067E+05, 0.60783E+05, 0.62539E+05, 0.64334E+05,
     + 0.66170E+05, 0.68047E+05, 0.69967E+05, 0.71929E+05, 0.73934E+05,
     + 0.75983E+05, 0.78078E+05, 0.80217E+05, 0.82403E+05, 0.84636E+05,
     + 0.86917E+05, 0.89246E+05, 0.91625E+05, 0.94053E+05, 0.96533E+05,
     + 0.99064E+05/
c...        --       182
      data (QofT( 5,J),J=1,119)/ 0.82770E+02, 0.13749E+03, 0.20083E+03,
     + 0.27176E+03, 0.34955E+03, 0.43370E+03, 0.52376E+03, 0.61944E+03,
     + 0.72050E+03, 0.82679E+03, 0.93821E+03, 0.10547E+04, 0.11763E+04,
     + 0.13031E+04, 0.14350E+04, 0.15723E+04, 0.17150E+04, 0.18633E+04,
     + 0.20172E+04, 0.21770E+04, 0.23429E+04, 0.25149E+04, 0.26934E+04,
     + 0.28784E+04, 0.30702E+04, 0.32690E+04, 0.34750E+04, 0.36885E+04,
     + 0.39096E+04, 0.41386E+04, 0.43758E+04, 0.46213E+04, 0.48755E+04,
     + 0.51386E+04, 0.54109E+04, 0.56927E+04, 0.59841E+04, 0.62856E+04,
     + 0.65973E+04, 0.69197E+04, 0.72529E+04, 0.75973E+04, 0.79533E+04,
     + 0.83210E+04, 0.87009E+04, 0.90933E+04, 0.94985E+04, 0.99168E+04,
     + 0.10348E+05, 0.10794E+05, 0.11254E+05, 0.11728E+05, 0.12217E+05,
     + 0.12722E+05, 0.13242E+05, 0.13778E+05, 0.14331E+05, 0.14900E+05,
     + 0.15486E+05, 0.16091E+05, 0.16713E+05, 0.17353E+05, 0.18012E+05,
     + 0.18691E+05, 0.19389E+05, 0.20108E+05, 0.20847E+05, 0.21607E+05,
     + 0.22388E+05, 0.23191E+05, 0.24017E+05, 0.24866E+05, 0.25738E+05,
     + 0.26633E+05, 0.27553E+05, 0.28498E+05, 0.29468E+05, 0.30464E+05,
     + 0.31486E+05, 0.32536E+05, 0.33612E+05, 0.34716E+05, 0.35849E+05,
     + 0.37011E+05, 0.38202E+05, 0.39424E+05, 0.40676E+05, 0.41959E+05,
     + 0.43274E+05, 0.44622E+05, 0.46002E+05, 0.47416E+05, 0.48864E+05,
     + 0.50348E+05, 0.51866E+05, 0.53421E+05, 0.55012E+05, 0.56640E+05,
     + 0.58307E+05, 0.60012E+05, 0.61757E+05, 0.63541E+05, 0.65366E+05,
     + 0.67233E+05, 0.69141E+05, 0.71092E+05, 0.73087E+05, 0.75125E+05,
     + 0.77209E+05, 0.79338E+05, 0.81513E+05, 0.83736E+05, 0.86006E+05,
     + 0.88324E+05, 0.90693E+05, 0.93111E+05, 0.95580E+05, 0.98100E+05,
     + 0.10067E+06/
c...        --       172
      data (QofT( 6,J),J=1,119)/ 0.49379E+03, 0.82021E+03, 0.11980E+04,
     + 0.16211E+04, 0.20851E+04, 0.25870E+04, 0.31242E+04, 0.36949E+04,
     + 0.42977E+04, 0.49317E+04, 0.55963E+04, 0.62911E+04, 0.70164E+04,
     + 0.77722E+04, 0.85591E+04, 0.93777E+04, 0.10228E+05, 0.11112E+05,
     + 0.12030E+05, 0.12983E+05, 0.13971E+05, 0.14997E+05, 0.16061E+05,
     + 0.17163E+05, 0.18306E+05, 0.19491E+05, 0.20719E+05, 0.21991E+05,
     + 0.23309E+05, 0.24673E+05, 0.26086E+05, 0.27549E+05, 0.29064E+05,
     + 0.30631E+05, 0.32254E+05, 0.33932E+05, 0.35669E+05, 0.37464E+05,
     + 0.39321E+05, 0.41242E+05, 0.43227E+05, 0.45279E+05, 0.47399E+05,
     + 0.49589E+05, 0.51852E+05, 0.54189E+05, 0.56602E+05, 0.59094E+05,
     + 0.61666E+05, 0.64320E+05, 0.67058E+05, 0.69883E+05, 0.72796E+05,
     + 0.75801E+05, 0.78899E+05, 0.82092E+05, 0.85382E+05, 0.88773E+05,
     + 0.92266E+05, 0.95863E+05, 0.99568E+05, 0.10338E+06, 0.10731E+06,
     + 0.11135E+06, 0.11551E+06, 0.11979E+06, 0.12419E+06, 0.12871E+06,
     + 0.13337E+06, 0.13815E+06, 0.14307E+06, 0.14812E+06, 0.15331E+06,
     + 0.15865E+06, 0.16412E+06, 0.16975E+06, 0.17553E+06, 0.18146E+06,
     + 0.18754E+06, 0.19379E+06, 0.20020E+06, 0.20678E+06, 0.21352E+06,
     + 0.22044E+06, 0.22753E+06, 0.23480E+06, 0.24226E+06, 0.24990E+06,
     + 0.25773E+06, 0.26575E+06, 0.27397E+06, 0.28239E+06, 0.29102E+06,
     + 0.29985E+06, 0.30889E+06, 0.31814E+06, 0.32762E+06, 0.33731E+06,
     + 0.34724E+06, 0.35739E+06, 0.36777E+06, 0.37840E+06, 0.38926E+06,
     + 0.40038E+06, 0.41174E+06, 0.42335E+06, 0.43523E+06, 0.44737E+06,
     + 0.45977E+06, 0.47245E+06, 0.48540E+06, 0.49863E+06, 0.51214E+06,
     + 0.52595E+06, 0.54005E+06, 0.55444E+06, 0.56914E+06, 0.58415E+06,
     + 0.59947E+06/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_CO2   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 9), QofT( 9,119),Q(NT)
      data xgj/ 1.,2.,1.,6.,2.,12.,1.,6.,1/
c...      CO2
c...        --       626
      data (QofT( 1,J),J=1,119)/ 0.53642E+02, 0.75947E+02, 0.98292E+02,
     + 0.12078E+03, 0.14364E+03, 0.16714E+03, 0.19160E+03, 0.21731E+03,
     + 0.24454E+03, 0.27355E+03, 0.30456E+03, 0.33778E+03, 0.37343E+03,
     + 0.41170E+03, 0.45280E+03, 0.49692E+03, 0.54427E+03, 0.59505E+03,
     + 0.64948E+03, 0.70779E+03, 0.77019E+03, 0.83693E+03, 0.90825E+03,
     + 0.98440E+03, 0.10656E+04, 0.11522E+04, 0.12445E+04, 0.13427E+04,
     + 0.14471E+04, 0.15580E+04, 0.16759E+04, 0.18009E+04, 0.19334E+04,
     + 0.20739E+04, 0.22225E+04, 0.23798E+04, 0.25462E+04, 0.27219E+04,
     + 0.29074E+04, 0.31032E+04, 0.33097E+04, 0.35272E+04, 0.37564E+04,
     + 0.39976E+04, 0.42514E+04, 0.45181E+04, 0.47985E+04, 0.50929E+04,
     + 0.54019E+04, 0.57260E+04, 0.60659E+04, 0.64221E+04, 0.67952E+04,
     + 0.71859E+04, 0.75946E+04, 0.80222E+04, 0.84691E+04, 0.89362E+04,
     + 0.94241E+04, 0.99335E+04, 0.10465E+05, 0.11020E+05, 0.11598E+05,
     + 0.12201E+05, 0.12828E+05, 0.13482E+05, 0.14163E+05, 0.14872E+05,
     + 0.15609E+05, 0.16376E+05, 0.17173E+05, 0.18001E+05, 0.18861E+05,
     + 0.19754E+05, 0.20682E+05, 0.21644E+05, 0.22643E+05, 0.23678E+05,
     + 0.24752E+05, 0.25865E+05, 0.27018E+05, 0.28212E+05, 0.29449E+05,
     + 0.30730E+05, 0.32055E+05, 0.33426E+05, 0.34845E+05, 0.36312E+05,
     + 0.37828E+05, 0.39395E+05, 0.41015E+05, 0.42688E+05, 0.44416E+05,
     + 0.46199E+05, 0.48041E+05, 0.49942E+05, 0.51902E+05, 0.53925E+05,
     + 0.56011E+05, 0.58162E+05, 0.60379E+05, 0.62664E+05, 0.65019E+05,
     + 0.67444E+05, 0.69942E+05, 0.72515E+05, 0.75163E+05, 0.77890E+05,
     + 0.80695E+05, 0.83582E+05, 0.86551E+05, 0.89605E+05, 0.92746E+05,
     + 0.95975E+05, 0.99294E+05, 0.10271E+06, 0.10621E+06, 0.10981E+06,
     + 0.11351E+06/
c...        --       636
      data (QofT( 2,J),J=1,119)/ 0.10728E+03, 0.15189E+03, 0.19659E+03,
     + 0.24164E+03, 0.28753E+03, 0.33486E+03, 0.38429E+03, 0.43643E+03,
     + 0.49184E+03, 0.55104E+03, 0.61449E+03, 0.68263E+03, 0.75589E+03,
     + 0.83468E+03, 0.91943E+03, 0.10106E+04, 0.11085E+04, 0.12137E+04,
     + 0.13266E+04, 0.14477E+04, 0.15774E+04, 0.17163E+04, 0.18649E+04,
     + 0.20237E+04, 0.21933E+04, 0.23743E+04, 0.25673E+04, 0.27729E+04,
     + 0.29917E+04, 0.32245E+04, 0.34718E+04, 0.37345E+04, 0.40132E+04,
     + 0.43087E+04, 0.46218E+04, 0.49533E+04, 0.53041E+04, 0.56749E+04,
     + 0.60668E+04, 0.64805E+04, 0.69171E+04, 0.73774E+04, 0.78626E+04,
     + 0.83736E+04, 0.89114E+04, 0.94772E+04, 0.10072E+05, 0.10697E+05,
     + 0.11353E+05, 0.12042E+05, 0.12765E+05, 0.13523E+05, 0.14317E+05,
     + 0.15148E+05, 0.16019E+05, 0.16930E+05, 0.17883E+05, 0.18879E+05,
     + 0.19920E+05, 0.21008E+05, 0.22143E+05, 0.23328E+05, 0.24563E+05,
     + 0.25852E+05, 0.27195E+05, 0.28594E+05, 0.30051E+05, 0.31568E+05,
     + 0.33146E+05, 0.34788E+05, 0.36496E+05, 0.38271E+05, 0.40115E+05,
     + 0.42031E+05, 0.44021E+05, 0.46086E+05, 0.48230E+05, 0.50453E+05,
     + 0.52759E+05, 0.55150E+05, 0.57628E+05, 0.60195E+05, 0.62854E+05,
     + 0.65608E+05, 0.68459E+05, 0.71409E+05, 0.74461E+05, 0.77618E+05,
     + 0.80883E+05, 0.84258E+05, 0.87746E+05, 0.91350E+05, 0.95073E+05,
     + 0.98918E+05, 0.10289E+06, 0.10698E+06, 0.11121E+06, 0.11558E+06,
     + 0.12008E+06, 0.12472E+06, 0.12950E+06, 0.13443E+06, 0.13952E+06,
     + 0.14475E+06, 0.15015E+06, 0.15571E+06, 0.16143E+06, 0.16732E+06,
     + 0.17338E+06, 0.17962E+06, 0.18604E+06, 0.19264E+06, 0.19943E+06,
     + 0.20642E+06, 0.21360E+06, 0.22098E+06, 0.22856E+06, 0.23636E+06,
     + 0.24436E+06/
c...        --       628
      data (QofT( 3,J),J=1,119)/ 0.11368E+03, 0.16096E+03, 0.20833E+03,
     + 0.25603E+03, 0.30452E+03, 0.35442E+03, 0.40640E+03, 0.46110E+03,
     + 0.51910E+03, 0.58093E+03, 0.64709E+03, 0.71804E+03, 0.79422E+03,
     + 0.87607E+03, 0.96402E+03, 0.10585E+04, 0.11600E+04, 0.12689E+04,
     + 0.13857E+04, 0.15108E+04, 0.16449E+04, 0.17883E+04, 0.19416E+04,
     + 0.21054E+04, 0.22803E+04, 0.24668E+04, 0.26655E+04, 0.28770E+04,
     + 0.31021E+04, 0.33414E+04, 0.35956E+04, 0.38654E+04, 0.41516E+04,
     + 0.44549E+04, 0.47761E+04, 0.51160E+04, 0.54755E+04, 0.58555E+04,
     + 0.62568E+04, 0.66804E+04, 0.71273E+04, 0.75982E+04, 0.80944E+04,
     + 0.86169E+04, 0.91666E+04, 0.97446E+04, 0.10352E+05, 0.10990E+05,
     + 0.11660E+05, 0.12363E+05, 0.13101E+05, 0.13874E+05, 0.14683E+05,
     + 0.15531E+05, 0.16418E+05, 0.17347E+05, 0.18317E+05, 0.19332E+05,
     + 0.20392E+05, 0.21499E+05, 0.22654E+05, 0.23859E+05, 0.25116E+05,
     + 0.26426E+05, 0.27792E+05, 0.29214E+05, 0.30695E+05, 0.32236E+05,
     + 0.33840E+05, 0.35508E+05, 0.37242E+05, 0.39045E+05, 0.40917E+05,
     + 0.42862E+05, 0.44881E+05, 0.46977E+05, 0.49152E+05, 0.51407E+05,
     + 0.53746E+05, 0.56171E+05, 0.58683E+05, 0.61286E+05, 0.63981E+05,
     + 0.66772E+05, 0.69661E+05, 0.72650E+05, 0.75742E+05, 0.78940E+05,
     + 0.82246E+05, 0.85664E+05, 0.89196E+05, 0.92845E+05, 0.96613E+05,
     + 0.10050E+06, 0.10452E+06, 0.10867E+06, 0.11295E+06, 0.11736E+06,
     + 0.12191E+06, 0.12661E+06, 0.13145E+06, 0.13643E+06, 0.14157E+06,
     + 0.14687E+06, 0.15232E+06, 0.15794E+06, 0.16372E+06, 0.16968E+06,
     + 0.17580E+06, 0.18211E+06, 0.18859E+06, 0.19526E+06, 0.20213E+06,
     + 0.20918E+06, 0.21643E+06, 0.22388E+06, 0.23154E+06, 0.23941E+06,
     + 0.24750E+06/
c...        --       627
      data (QofT( 4,J),J=1,119)/ 0.66338E+03, 0.93923E+03, 0.12156E+04,
     + 0.14938E+04, 0.17766E+04, 0.20676E+04, 0.23705E+04, 0.26891E+04,
     + 0.30267E+04, 0.33866E+04, 0.37714E+04, 0.41839E+04, 0.46267E+04,
     + 0.51023E+04, 0.56132E+04, 0.61618E+04, 0.67508E+04, 0.73827E+04,
     + 0.80603E+04, 0.87863E+04, 0.95636E+04, 0.10395E+05, 0.11284E+05,
     + 0.12233E+05, 0.13246E+05, 0.14326E+05, 0.15477E+05, 0.16702E+05,
     + 0.18005E+05, 0.19390E+05, 0.20861E+05, 0.22422E+05, 0.24077E+05,
     + 0.25832E+05, 0.27689E+05, 0.29655E+05, 0.31734E+05, 0.33931E+05,
     + 0.36250E+05, 0.38698E+05, 0.41280E+05, 0.44002E+05, 0.46869E+05,
     + 0.49886E+05, 0.53062E+05, 0.56400E+05, 0.59909E+05, 0.63594E+05,
     + 0.67462E+05, 0.71521E+05, 0.75777E+05, 0.80238E+05, 0.84911E+05,
     + 0.89804E+05, 0.94925E+05, 0.10028E+06, 0.10588E+06, 0.11173E+06,
     + 0.11785E+06, 0.12423E+06, 0.13090E+06, 0.13785E+06, 0.14510E+06,
     + 0.15265E+06, 0.16053E+06, 0.16873E+06, 0.17727E+06, 0.18615E+06,
     + 0.19540E+06, 0.20501E+06, 0.21501E+06, 0.22540E+06, 0.23619E+06,
     + 0.24740E+06, 0.25904E+06, 0.27112E+06, 0.28365E+06, 0.29664E+06,
     + 0.31012E+06, 0.32409E+06, 0.33856E+06, 0.35356E+06, 0.36908E+06,
     + 0.38516E+06, 0.40180E+06, 0.41902E+06, 0.43683E+06, 0.45525E+06,
     + 0.47429E+06, 0.49397E+06, 0.51431E+06, 0.53532E+06, 0.55702E+06,
     + 0.57943E+06, 0.60256E+06, 0.62644E+06, 0.65107E+06, 0.67648E+06,
     + 0.70269E+06, 0.72972E+06, 0.75758E+06, 0.78629E+06, 0.81588E+06,
     + 0.84636E+06, 0.87775E+06, 0.91008E+06, 0.94337E+06, 0.97763E+06,
     + 0.10129E+07, 0.10492E+07, 0.10865E+07, 0.11249E+07, 0.11644E+07,
     + 0.12050E+07, 0.12467E+07, 0.12896E+07, 0.13337E+07, 0.13789E+07,
     + 0.14255E+07/
c...        --       638
      data (QofT( 5,J),J=1,119)/ 0.22737E+03, 0.32194E+03, 0.41671E+03,
     + 0.51226E+03, 0.60963E+03, 0.71017E+03, 0.81528E+03, 0.92628E+03,
     + 0.10444E+04, 0.11707E+04, 0.13061E+04, 0.14518E+04, 0.16085E+04,
     + 0.17772E+04, 0.19588E+04, 0.21542E+04, 0.23644E+04, 0.25903E+04,
     + 0.28330E+04, 0.30934E+04, 0.33726E+04, 0.36717E+04, 0.39918E+04,
     + 0.43342E+04, 0.47001E+04, 0.50907E+04, 0.55074E+04, 0.59515E+04,
     + 0.64244E+04, 0.69276E+04, 0.74626E+04, 0.80310E+04, 0.86344E+04,
     + 0.92744E+04, 0.99528E+04, 0.10671E+05, 0.11432E+05, 0.12236E+05,
     + 0.13086E+05, 0.13984E+05, 0.14932E+05, 0.15932E+05, 0.16985E+05,
     + 0.18096E+05, 0.19265E+05, 0.20495E+05, 0.21788E+05, 0.23148E+05,
     + 0.24576E+05, 0.26075E+05, 0.27648E+05, 0.29298E+05, 0.31027E+05,
     + 0.32839E+05, 0.34736E+05, 0.36721E+05, 0.38798E+05, 0.40970E+05,
     + 0.43240E+05, 0.45611E+05, 0.48087E+05, 0.50671E+05, 0.53368E+05,
     + 0.56180E+05, 0.59111E+05, 0.62165E+05, 0.65347E+05, 0.68659E+05,
     + 0.72107E+05, 0.75694E+05, 0.79425E+05, 0.83303E+05, 0.87334E+05,
     + 0.91522E+05, 0.95872E+05, 0.10039E+06, 0.10507E+06, 0.10994E+06,
     + 0.11498E+06, 0.12021E+06, 0.12563E+06, 0.13125E+06, 0.13707E+06,
     + 0.14309E+06, 0.14933E+06, 0.15579E+06, 0.16247E+06, 0.16938E+06,
     + 0.17653E+06, 0.18392E+06, 0.19156E+06, 0.19946E+06, 0.20761E+06,
     + 0.21604E+06, 0.22473E+06, 0.23371E+06, 0.24298E+06, 0.25254E+06,
     + 0.26240E+06, 0.27258E+06, 0.28307E+06, 0.29388E+06, 0.30502E+06,
     + 0.31651E+06, 0.32834E+06, 0.34052E+06, 0.35307E+06, 0.36599E+06,
     + 0.37929E+06, 0.39298E+06, 0.40706E+06, 0.42155E+06, 0.43645E+06,
     + 0.45178E+06, 0.46753E+06, 0.48373E+06, 0.50038E+06, 0.51748E+06,
     + 0.53506E+06/
c...        --       637
      data (QofT( 6,J),J=1,119)/ 0.13267E+04, 0.18785E+04, 0.24314E+04,
     + 0.29888E+04, 0.35566E+04, 0.41426E+04, 0.47550E+04, 0.54013E+04,
     + 0.60886E+04, 0.68232E+04, 0.76109E+04, 0.84574E+04, 0.93678E+04,
     + 0.10348E+05, 0.11402E+05, 0.12536E+05, 0.13755E+05, 0.15065E+05,
     + 0.16471E+05, 0.17980E+05, 0.19598E+05, 0.21330E+05, 0.23184E+05,
     + 0.25166E+05, 0.27283E+05, 0.29543E+05, 0.31953E+05, 0.34521E+05,
     + 0.37256E+05, 0.40164E+05, 0.43256E+05, 0.46541E+05, 0.50026E+05,
     + 0.53723E+05, 0.57641E+05, 0.61790E+05, 0.66180E+05, 0.70823E+05,
     + 0.75729E+05, 0.80910E+05, 0.86378E+05, 0.92145E+05, 0.98224E+05,
     + 0.10463E+06, 0.11137E+06, 0.11846E+06, 0.12592E+06, 0.13375E+06,
     + 0.14198E+06, 0.15062E+06, 0.15969E+06, 0.16920E+06, 0.17916E+06,
     + 0.18959E+06, 0.20052E+06, 0.21196E+06, 0.22392E+06, 0.23642E+06,
     + 0.24949E+06, 0.26314E+06, 0.27740E+06, 0.29227E+06, 0.30779E+06,
     + 0.32398E+06, 0.34085E+06, 0.35842E+06, 0.37673E+06, 0.39579E+06,
     + 0.41563E+06, 0.43626E+06, 0.45772E+06, 0.48003E+06, 0.50322E+06,
     + 0.52730E+06, 0.55232E+06, 0.57829E+06, 0.60524E+06, 0.63320E+06,
     + 0.66219E+06, 0.69226E+06, 0.72342E+06, 0.75571E+06, 0.78916E+06,
     + 0.82380E+06, 0.85966E+06, 0.89678E+06, 0.93518E+06, 0.97490E+06,
     + 0.10160E+07, 0.10585E+07, 0.11023E+07, 0.11477E+07, 0.11946E+07,
     + 0.12430E+07, 0.12929E+07, 0.13445E+07, 0.13977E+07, 0.14526E+07,
     + 0.15093E+07, 0.15677E+07, 0.16280E+07, 0.16901E+07, 0.17541E+07,
     + 0.18200E+07, 0.18880E+07, 0.19579E+07, 0.20300E+07, 0.21042E+07,
     + 0.21805E+07, 0.22591E+07, 0.23400E+07, 0.24232E+07, 0.25087E+07,
     + 0.25967E+07, 0.26871E+07, 0.27801E+07, 0.28757E+07, 0.29739E+07,
     + 0.30747E+07/
c...        --       828
      data (QofT( 7,J),J=1,119)/ 0.60334E+02, 0.85430E+02, 0.11058E+03,
     + 0.13590E+03, 0.16167E+03, 0.18821E+03, 0.21588E+03, 0.24502E+03,
     + 0.27595E+03, 0.30896E+03, 0.34431E+03, 0.38225E+03, 0.42301E+03,
     + 0.46684E+03, 0.51397E+03, 0.56464E+03, 0.61907E+03, 0.67753E+03,
     + 0.74027E+03, 0.80753E+03, 0.87961E+03, 0.95676E+03, 0.10393E+04,
     + 0.11275E+04, 0.12217E+04, 0.13222E+04, 0.14293E+04, 0.15434E+04,
     + 0.16648E+04, 0.17940E+04, 0.19312E+04, 0.20769E+04, 0.22315E+04,
     + 0.23954E+04, 0.25691E+04, 0.27529E+04, 0.29474E+04, 0.31530E+04,
     + 0.33702E+04, 0.35995E+04, 0.38414E+04, 0.40965E+04, 0.43654E+04,
     + 0.46484E+04, 0.49464E+04, 0.52598E+04, 0.55892E+04, 0.59353E+04,
     + 0.62988E+04, 0.66803E+04, 0.70804E+04, 0.74998E+04, 0.79394E+04,
     + 0.83998E+04, 0.88817E+04, 0.93859E+04, 0.99132E+04, 0.10464E+05,
     + 0.11040E+05, 0.11642E+05, 0.12270E+05, 0.12925E+05, 0.13609E+05,
     + 0.14321E+05, 0.15064E+05, 0.15838E+05, 0.16643E+05, 0.17482E+05,
     + 0.18355E+05, 0.19263E+05, 0.20207E+05, 0.21188E+05, 0.22208E+05,
     + 0.23267E+05, 0.24366E+05, 0.25508E+05, 0.26692E+05, 0.27921E+05,
     + 0.29195E+05, 0.30516E+05, 0.31886E+05, 0.33304E+05, 0.34773E+05,
     + 0.36294E+05, 0.37869E+05, 0.39499E+05, 0.41185E+05, 0.42929E+05,
     + 0.44732E+05, 0.46596E+05, 0.48522E+05, 0.50513E+05, 0.52569E+05,
     + 0.54692E+05, 0.56884E+05, 0.59146E+05, 0.61481E+05, 0.63890E+05,
     + 0.66375E+05, 0.68937E+05, 0.71578E+05, 0.74301E+05, 0.77107E+05,
     + 0.79998E+05, 0.82976E+05, 0.86043E+05, 0.89201E+05, 0.92452E+05,
     + 0.95799E+05, 0.99242E+05, 0.10278E+06, 0.10643E+06, 0.11018E+06,
     + 0.11403E+06, 0.11799E+06, 0.12206E+06, 0.12625E+06, 0.13055E+06,
     + 0.13497E+06/
c...        --       728
      data (QofT( 8,J),J=1,119)/ 0.70354E+03, 0.99615E+03, 0.12893E+04,
     + 0.15846E+04, 0.18848E+04, 0.21940E+04, 0.25162E+04, 0.28554E+04,
     + 0.32152E+04, 0.35991E+04, 0.40099E+04, 0.44507E+04, 0.49242E+04,
     + 0.54332E+04, 0.59802E+04, 0.65681E+04, 0.71996E+04, 0.78776E+04,
     + 0.86050E+04, 0.93847E+04, 0.10220E+05, 0.11114E+05, 0.12070E+05,
     + 0.13091E+05, 0.14182E+05, 0.15345E+05, 0.16585E+05, 0.17906E+05,
     + 0.19311E+05, 0.20805E+05, 0.22393E+05, 0.24078E+05, 0.25865E+05,
     + 0.27760E+05, 0.29768E+05, 0.31893E+05, 0.34140E+05, 0.36516E+05,
     + 0.39025E+05, 0.41674E+05, 0.44469E+05, 0.47416E+05, 0.50520E+05,
     + 0.53789E+05, 0.57229E+05, 0.60847E+05, 0.64650E+05, 0.68645E+05,
     + 0.72840E+05, 0.77242E+05, 0.81859E+05, 0.86699E+05, 0.91770E+05,
     + 0.97081E+05, 0.10264E+06, 0.10846E+06, 0.11454E+06, 0.12090E+06,
     + 0.12754E+06, 0.13447E+06, 0.14171E+06, 0.14927E+06, 0.15715E+06,
     + 0.16536E+06, 0.17392E+06, 0.18284E+06, 0.19213E+06, 0.20179E+06,
     + 0.21185E+06, 0.22231E+06, 0.23319E+06, 0.24450E+06, 0.25625E+06,
     + 0.26845E+06, 0.28112E+06, 0.29427E+06, 0.30791E+06, 0.32206E+06,
     + 0.33674E+06, 0.35196E+06, 0.36772E+06, 0.38406E+06, 0.40098E+06,
     + 0.41850E+06, 0.43663E+06, 0.45539E+06, 0.47480E+06, 0.49488E+06,
     + 0.51564E+06, 0.53710E+06, 0.55928E+06, 0.58219E+06, 0.60586E+06,
     + 0.63029E+06, 0.65553E+06, 0.68157E+06, 0.70844E+06, 0.73616E+06,
     + 0.76476E+06, 0.79424E+06, 0.82464E+06, 0.85597E+06, 0.88826E+06,
     + 0.92153E+06, 0.95580E+06, 0.99108E+06, 0.10274E+07, 0.10648E+07,
     + 0.11033E+07, 0.11429E+07, 0.11837E+07, 0.12256E+07, 0.12687E+07,
     + 0.13131E+07, 0.13586E+07, 0.14055E+07, 0.14536E+07, 0.15031E+07,
     + 0.15539E+07/
c...        --       727
      data (QofT( 9,J),J=1,119)/ 0.20518E+04, 0.29051E+04, 0.37601E+04,
     + 0.46209E+04, 0.54961E+04, 0.63969E+04, 0.73353E+04, 0.83227E+04,
     + 0.93698E+04, 0.10486E+05, 0.11681E+05, 0.12962E+05, 0.14337E+05,
     + 0.15815E+05, 0.17403E+05, 0.19110E+05, 0.20942E+05, 0.22909E+05,
     + 0.25018E+05, 0.27278E+05, 0.29699E+05, 0.32290E+05, 0.35060E+05,
     + 0.38019E+05, 0.41177E+05, 0.44545E+05, 0.48135E+05, 0.51957E+05,
     + 0.56023E+05, 0.60346E+05, 0.64938E+05, 0.69812E+05, 0.74981E+05,
     + 0.80461E+05, 0.86264E+05, 0.92406E+05, 0.98902E+05, 0.10577E+06,
     + 0.11302E+06, 0.12067E+06, 0.12875E+06, 0.13726E+06, 0.14622E+06,
     + 0.15566E+06, 0.16559E+06, 0.17604E+06, 0.18702E+06, 0.19855E+06,
     + 0.21066E+06, 0.22336E+06, 0.23669E+06, 0.25065E+06, 0.26528E+06,
     + 0.28061E+06, 0.29664E+06, 0.31342E+06, 0.33096E+06, 0.34930E+06,
     + 0.36845E+06, 0.38845E+06, 0.40933E+06, 0.43111E+06, 0.45383E+06,
     + 0.47751E+06, 0.50219E+06, 0.52790E+06, 0.55466E+06, 0.58252E+06,
     + 0.61151E+06, 0.64166E+06, 0.67300E+06, 0.70558E+06, 0.73943E+06,
     + 0.77458E+06, 0.81108E+06, 0.84896E+06, 0.88827E+06, 0.92904E+06,
     + 0.97131E+06, 0.10151E+07, 0.10605E+07, 0.11076E+07, 0.11563E+07,
     + 0.12068E+07, 0.12590E+07, 0.13130E+07, 0.13689E+07, 0.14267E+07,
     + 0.14865E+07, 0.15483E+07, 0.16121E+07, 0.16781E+07, 0.17462E+07,
     + 0.18165E+07, 0.18892E+07, 0.19641E+07, 0.20415E+07, 0.21213E+07,
     + 0.22036E+07, 0.22884E+07, 0.23759E+07, 0.24661E+07, 0.25590E+07,
     + 0.26547E+07, 0.27533E+07, 0.28549E+07, 0.29594E+07, 0.30670E+07,
     + 0.31778E+07, 0.32918E+07, 0.34090E+07, 0.35296E+07, 0.36536E+07,
     + 0.37812E+07, 0.39123E+07, 0.40470E+07, 0.41855E+07, 0.43278E+07,
     + 0.44739E+07/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_O3    (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
c 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj(18), QofT(18,119),Q(NT)
      data xgj/ 1.,1.,1.,6.,6.,1.,1.,6.,6.,6.,36.,
     + 1.,1.,6.,6.,36.,1.,6./
c...       O3
c...        --       666
      data (QofT( 1,J),J=1,119)/ 0.30333E+03, 0.51126E+03, 0.75274E+03,
     + 0.10241E+04, 0.13236E+04, 0.16508E+04, 0.20068E+04, 0.23935E+04,
     + 0.28136E+04, 0.32703E+04, 0.37672E+04, 0.43082E+04, 0.48975E+04,
     + 0.55395E+04, 0.62386E+04, 0.69996E+04, 0.78272E+04, 0.87264E+04,
     + 0.97026E+04, 0.10761E+05, 0.11907E+05, 0.13146E+05, 0.14485E+05,
     + 0.15929E+05, 0.17484E+05, 0.19158E+05, 0.20957E+05, 0.22887E+05,
     + 0.24956E+05, 0.27172E+05, 0.29541E+05, 0.32072E+05, 0.34773E+05,
     + 0.37652E+05, 0.40718E+05, 0.43979E+05, 0.47444E+05, 0.51123E+05,
     + 0.55026E+05, 0.59161E+05, 0.63540E+05, 0.68172E+05, 0.73069E+05,
     + 0.78240E+05, 0.83698E+05, 0.89453E+05, 0.95517E+05, 0.10190E+06,
     + 0.10862E+06, 0.11569E+06, 0.12311E+06, 0.13091E+06, 0.13909E+06,
     + 0.14767E+06, 0.15666E+06, 0.16608E+06, 0.17594E+06, 0.18626E+06,
     + 0.19706E+06, 0.20834E+06, 0.22012E+06, 0.23242E+06, 0.24526E+06,
     + 0.25866E+06, 0.27262E+06, 0.28717E+06, 0.30233E+06, 0.31811E+06,
     + 0.33453E+06, 0.35161E+06, 0.36937E+06, 0.38784E+06, 0.40702E+06,
     + 0.42694E+06, 0.44762E+06, 0.46909E+06, 0.49135E+06, 0.51444E+06,
     + 0.53838E+06, 0.56318E+06, 0.58887E+06, 0.61548E+06, 0.64303E+06,
     + 0.67153E+06, 0.70102E+06, 0.73153E+06, 0.76306E+06, 0.79566E+06,
     + 0.82934E+06, 0.86413E+06, 0.90006E+06, 0.93716E+06, 0.97545E+06,
     + 0.10150E+07, 0.10557E+07, 0.10977E+07, 0.11411E+07, 0.11858E+07,
     + 0.12318E+07, 0.12792E+07, 0.13281E+07, 0.13784E+07, 0.14302E+07,
     + 0.14835E+07, 0.15384E+07, 0.15948E+07, 0.16529E+07, 0.17126E+07,
     + 0.17740E+07, 0.18371E+07, 0.19020E+07, 0.19686E+07, 0.20371E+07,
     + 0.21074E+07, 0.21797E+07, 0.22538E+07, 0.23300E+07, 0.24081E+07,
     + 0.24883E+07/
c...        --       668
      data (QofT( 2,J),J=1,119)/ 0.64763E+03, 0.10916E+04, 0.16073E+04,
     + 0.21870E+04, 0.28271E+04, 0.35272E+04, 0.42900E+04, 0.51197E+04,
     + 0.60225E+04, 0.70057E+04, 0.80771E+04, 0.92455E+04, 0.10520E+05,
     + 0.11911E+05, 0.13427E+05, 0.15079E+05, 0.16878E+05, 0.18834E+05,
     + 0.20960E+05, 0.23267E+05, 0.25767E+05, 0.28472E+05, 0.31397E+05,
     + 0.34553E+05, 0.37957E+05, 0.41620E+05, 0.45559E+05, 0.49790E+05,
     + 0.54327E+05, 0.59187E+05, 0.64387E+05, 0.69944E+05, 0.75877E+05,
     + 0.82203E+05, 0.88943E+05, 0.96114E+05, 0.10374E+06, 0.11184E+06,
     + 0.12043E+06, 0.12954E+06, 0.13918E+06, 0.14939E+06, 0.16018E+06,
     + 0.17159E+06, 0.18362E+06, 0.19632E+06, 0.20970E+06, 0.22380E+06,
     + 0.23863E+06, 0.25423E+06, 0.27063E+06, 0.28786E+06, 0.30594E+06,
     + 0.32490E+06, 0.34478E+06, 0.36561E+06, 0.38743E+06, 0.41026E+06,
     + 0.43413E+06, 0.45909E+06, 0.48517E+06, 0.51241E+06, 0.54084E+06,
     + 0.57049E+06, 0.60141E+06, 0.63365E+06, 0.66722E+06, 0.70219E+06,
     + 0.73858E+06, 0.77644E+06, 0.81581E+06, 0.85674E+06, 0.89927E+06,
     + 0.94345E+06, 0.98932E+06, 0.10369E+07, 0.10863E+07, 0.11375E+07,
     + 0.11906E+07, 0.12457E+07, 0.13027E+07, 0.13618E+07, 0.14229E+07,
     + 0.14862E+07, 0.15517E+07, 0.16194E+07, 0.16894E+07, 0.17618E+07,
     + 0.18366E+07, 0.19139E+07, 0.19937E+07, 0.20761E+07, 0.21612E+07,
     + 0.22490E+07, 0.23395E+07, 0.24330E+07, 0.25293E+07, 0.26286E+07,
     + 0.27309E+07, 0.28363E+07, 0.29449E+07, 0.30568E+07, 0.31720E+07,
     + 0.32905E+07, 0.34125E+07, 0.35381E+07, 0.36672E+07, 0.38000E+07,
     + 0.39366E+07, 0.40770E+07, 0.42213E+07, 0.43696E+07, 0.45220E+07,
     + 0.46785E+07, 0.48392E+07, 0.50043E+07, 0.51737E+07, 0.53476E+07,
     + 0.55261E+07/
c...        --       686
      data (QofT( 3,J),J=1,119)/ 0.31656E+03, 0.53355E+03, 0.78557E+03,
     + 0.10688E+04, 0.13815E+04, 0.17235E+04, 0.20960E+04, 0.25011E+04,
     + 0.29420E+04, 0.34223E+04, 0.39459E+04, 0.45172E+04, 0.51408E+04,
     + 0.58213E+04, 0.65639E+04, 0.73735E+04, 0.82555E+04, 0.92152E+04,
     + 0.10259E+05, 0.11391E+05, 0.12619E+05, 0.13949E+05, 0.15387E+05,
     + 0.16940E+05, 0.18614E+05, 0.20417E+05, 0.22357E+05, 0.24440E+05,
     + 0.26675E+05, 0.29070E+05, 0.31633E+05, 0.34374E+05, 0.37299E+05,
     + 0.40420E+05, 0.43746E+05, 0.47285E+05, 0.51049E+05, 0.55047E+05,
     + 0.59289E+05, 0.63788E+05, 0.68554E+05, 0.73598E+05, 0.78932E+05,
     + 0.84568E+05, 0.90519E+05, 0.96796E+05, 0.10341E+06, 0.11039E+06,
     + 0.11772E+06, 0.12544E+06, 0.13356E+06, 0.14208E+06, 0.15103E+06,
     + 0.16041E+06, 0.17026E+06, 0.18057E+06, 0.19137E+06, 0.20268E+06,
     + 0.21450E+06, 0.22687E+06, 0.23979E+06, 0.25328E+06, 0.26736E+06,
     + 0.28206E+06, 0.29738E+06, 0.31336E+06, 0.33000E+06, 0.34733E+06,
     + 0.36537E+06, 0.38414E+06, 0.40366E+06, 0.42396E+06, 0.44505E+06,
     + 0.46696E+06, 0.48971E+06, 0.51332E+06, 0.53782E+06, 0.56323E+06,
     + 0.58958E+06, 0.61689E+06, 0.64518E+06, 0.67448E+06, 0.70482E+06,
     + 0.73623E+06, 0.76872E+06, 0.80234E+06, 0.83710E+06, 0.87303E+06,
     + 0.91017E+06, 0.94853E+06, 0.98816E+06, 0.10291E+07, 0.10713E+07,
     + 0.11149E+07, 0.11599E+07, 0.12063E+07, 0.12541E+07, 0.13034E+07,
     + 0.13542E+07, 0.14066E+07, 0.14606E+07, 0.15161E+07, 0.15733E+07,
     + 0.16322E+07, 0.16928E+07, 0.17552E+07, 0.18194E+07, 0.18854E+07,
     + 0.19532E+07, 0.20230E+07, 0.20947E+07, 0.21684E+07, 0.22441E+07,
     + 0.23219E+07, 0.24018E+07, 0.24838E+07, 0.25680E+07, 0.26545E+07,
     + 0.27432E+07/
c...        --       667
      data (QofT( 4,J),J=1,119)/ 0.37657E+04, 0.63472E+04, 0.93454E+04,
     + 0.12715E+05, 0.16435E+05, 0.20502E+05, 0.24929E+05, 0.29742E+05,
     + 0.34975E+05, 0.40668E+05, 0.46868E+05, 0.53624E+05, 0.60990E+05,
     + 0.69018E+05, 0.77768E+05, 0.87296E+05, 0.97666E+05, 0.10894E+06,
     + 0.12118E+06, 0.13446E+06, 0.14885E+06, 0.16441E+06, 0.18123E+06,
     + 0.19938E+06, 0.21894E+06, 0.23998E+06, 0.26261E+06, 0.28690E+06,
     + 0.31295E+06, 0.34084E+06, 0.37068E+06, 0.40256E+06, 0.43659E+06,
     + 0.47287E+06, 0.51151E+06, 0.55262E+06, 0.59632E+06, 0.64272E+06,
     + 0.69194E+06, 0.74412E+06, 0.79937E+06, 0.85783E+06, 0.91963E+06,
     + 0.98492E+06, 0.10538E+07, 0.11265E+07, 0.12031E+07, 0.12837E+07,
     + 0.13686E+07, 0.14579E+07, 0.15517E+07, 0.16502E+07, 0.17536E+07,
     + 0.18621E+07, 0.19758E+07, 0.20949E+07, 0.22196E+07, 0.23501E+07,
     + 0.24866E+07, 0.26292E+07, 0.27783E+07, 0.29339E+07, 0.30963E+07,
     + 0.32658E+07, 0.34425E+07, 0.36266E+07, 0.38184E+07, 0.40181E+07,
     + 0.42260E+07, 0.44422E+07, 0.46671E+07, 0.49008E+07, 0.51437E+07,
     + 0.53959E+07, 0.56578E+07, 0.59296E+07, 0.62116E+07, 0.65040E+07,
     + 0.68071E+07, 0.71213E+07, 0.74468E+07, 0.77838E+07, 0.81328E+07,
     + 0.84939E+07, 0.88676E+07, 0.92541E+07, 0.96536E+07, 0.10067E+08,
     + 0.10493E+08, 0.10934E+08, 0.11390E+08, 0.11860E+08, 0.12345E+08,
     + 0.12846E+08, 0.13363E+08, 0.13895E+08, 0.14445E+08, 0.15011E+08,
     + 0.15595E+08, 0.16196E+08, 0.16815E+08, 0.17453E+08, 0.18110E+08,
     + 0.18786E+08, 0.19482E+08, 0.20198E+08, 0.20934E+08, 0.21691E+08,
     + 0.22470E+08, 0.23270E+08, 0.24093E+08, 0.24939E+08, 0.25807E+08,
     + 0.26699E+08, 0.27616E+08, 0.28556E+08, 0.29522E+08, 0.30514E+08,
     + 0.31531E+08/
c...        --       676
      data (QofT( 5,J),J=1,119)/ 0.18608E+04, 0.31363E+04, 0.46177E+04,
     + 0.62826E+04, 0.81202E+04, 0.10129E+05, 0.12316E+05, 0.14693E+05,
     + 0.17277E+05, 0.20089E+05, 0.23153E+05, 0.26492E+05, 0.30133E+05,
     + 0.34103E+05, 0.38430E+05, 0.43145E+05, 0.48277E+05, 0.53858E+05,
     + 0.59920E+05, 0.66497E+05, 0.73624E+05, 0.81336E+05, 0.89671E+05,
     + 0.98668E+05, 0.10836E+06, 0.11880E+06, 0.13002E+06, 0.14207E+06,
     + 0.15500E+06, 0.16884E+06, 0.18365E+06, 0.19947E+06, 0.21636E+06,
     + 0.23438E+06, 0.25356E+06, 0.27398E+06, 0.29568E+06, 0.31873E+06,
     + 0.34318E+06, 0.36911E+06, 0.39656E+06, 0.42561E+06, 0.45632E+06,
     + 0.48877E+06, 0.52302E+06, 0.55914E+06, 0.59722E+06, 0.63732E+06,
     + 0.67952E+06, 0.72390E+06, 0.77055E+06, 0.81954E+06, 0.87097E+06,
     + 0.92491E+06, 0.98146E+06, 0.10407E+07, 0.11027E+07, 0.11677E+07,
     + 0.12356E+07, 0.13066E+07, 0.13807E+07, 0.14582E+07, 0.15390E+07,
     + 0.16233E+07, 0.17113E+07, 0.18029E+07, 0.18984E+07, 0.19978E+07,
     + 0.21012E+07, 0.22089E+07, 0.23208E+07, 0.24372E+07, 0.25581E+07,
     + 0.26837E+07, 0.28141E+07, 0.29494E+07, 0.30898E+07, 0.32354E+07,
     + 0.33864E+07, 0.35428E+07, 0.37049E+07, 0.38728E+07, 0.40466E+07,
     + 0.42264E+07, 0.44125E+07, 0.46050E+07, 0.48040E+07, 0.50098E+07,
     + 0.52224E+07, 0.54420E+07, 0.56689E+07, 0.59031E+07, 0.61449E+07,
     + 0.63943E+07, 0.66517E+07, 0.69172E+07, 0.71909E+07, 0.74731E+07,
     + 0.77639E+07, 0.80635E+07, 0.83721E+07, 0.86900E+07, 0.90172E+07,
     + 0.93541E+07, 0.97008E+07, 0.10058E+08, 0.10424E+08, 0.10802E+08,
     + 0.11190E+08, 0.11589E+08, 0.11999E+08, 0.12420E+08, 0.12853E+08,
     + 0.13298E+08, 0.13755E+08, 0.14223E+08, 0.14705E+08, 0.15199E+08,
     + 0.15706E+08/
c...        --       886
      data (QofT( 6,J),J=1,119)/ 0.67639E+03, 0.11401E+04, 0.16787E+04,
     + 0.22843E+04, 0.29532E+04, 0.36856E+04, 0.44842E+04, 0.53545E+04,
     + 0.63030E+04, 0.73381E+04, 0.84686E+04, 0.97040E+04, 0.11054E+05,
     + 0.12530E+05, 0.14143E+05, 0.15903E+05, 0.17823E+05, 0.19915E+05,
     + 0.22190E+05, 0.24663E+05, 0.27346E+05, 0.30254E+05, 0.33400E+05,
     + 0.36800E+05, 0.40469E+05, 0.44423E+05, 0.48678E+05, 0.53251E+05,
     + 0.58160E+05, 0.63423E+05, 0.69058E+05, 0.75085E+05, 0.81524E+05,
     + 0.88395E+05, 0.95719E+05, 0.10352E+06, 0.11181E+06, 0.12063E+06,
     + 0.12999E+06, 0.13991E+06, 0.15043E+06, 0.16157E+06, 0.17335E+06,
     + 0.18580E+06, 0.19895E+06, 0.21283E+06, 0.22746E+06, 0.24288E+06,
     + 0.25911E+06, 0.27619E+06, 0.29415E+06, 0.31301E+06, 0.33283E+06,
     + 0.35362E+06, 0.37542E+06, 0.39827E+06, 0.42221E+06, 0.44726E+06,
     + 0.47348E+06, 0.50089E+06, 0.52954E+06, 0.55947E+06, 0.59072E+06,
     + 0.62332E+06, 0.65733E+06, 0.69279E+06, 0.72973E+06, 0.76821E+06,
     + 0.80827E+06, 0.84996E+06, 0.89332E+06, 0.93840E+06, 0.98526E+06,
     + 0.10339E+07, 0.10845E+07, 0.11370E+07, 0.11914E+07, 0.12479E+07,
     + 0.13065E+07, 0.13672E+07, 0.14302E+07, 0.14953E+07, 0.15628E+07,
     + 0.16327E+07, 0.17050E+07, 0.17798E+07, 0.18571E+07, 0.19371E+07,
     + 0.20197E+07, 0.21051E+07, 0.21933E+07, 0.22844E+07, 0.23785E+07,
     + 0.24755E+07, 0.25757E+07, 0.26790E+07, 0.27855E+07, 0.28954E+07,
     + 0.30086E+07, 0.31253E+07, 0.32455E+07, 0.33693E+07, 0.34967E+07,
     + 0.36280E+07, 0.37631E+07, 0.39021E+07, 0.40451E+07, 0.41922E+07,
     + 0.43435E+07, 0.44990E+07, 0.46589E+07, 0.48232E+07, 0.49920E+07,
     + 0.51654E+07, 0.53436E+07, 0.55265E+07, 0.57143E+07, 0.59071E+07,
     + 0.61050E+07/
c...        --       868
      data (QofT( 7,J),J=1,119)/ 0.34615E+03, 0.58348E+03, 0.85915E+03,
     + 0.11692E+04, 0.15117E+04, 0.18868E+04, 0.22960E+04, 0.27419E+04,
     + 0.32278E+04, 0.37579E+04, 0.43366E+04, 0.49686E+04, 0.56591E+04,
     + 0.64134E+04, 0.72369E+04, 0.81354E+04, 0.91148E+04, 0.10181E+05,
     + 0.11341E+05, 0.12600E+05, 0.13966E+05, 0.15446E+05, 0.17046E+05,
     + 0.18775E+05, 0.20640E+05, 0.22649E+05, 0.24810E+05, 0.27132E+05,
     + 0.29624E+05, 0.32295E+05, 0.35154E+05, 0.38211E+05, 0.41475E+05,
     + 0.44958E+05, 0.48670E+05, 0.52621E+05, 0.56823E+05, 0.61288E+05,
     + 0.66026E+05, 0.71052E+05, 0.76376E+05, 0.82011E+05, 0.87972E+05,
     + 0.94271E+05, 0.10092E+06, 0.10794E+06, 0.11534E+06, 0.12313E+06,
     + 0.13134E+06, 0.13997E+06, 0.14905E+06, 0.15858E+06, 0.16859E+06,
     + 0.17909E+06, 0.19010E+06, 0.20164E+06, 0.21373E+06, 0.22638E+06,
     + 0.23962E+06, 0.25346E+06, 0.26792E+06, 0.28302E+06, 0.29879E+06,
     + 0.31524E+06, 0.33240E+06, 0.35029E+06, 0.36892E+06, 0.38833E+06,
     + 0.40853E+06, 0.42956E+06, 0.45142E+06, 0.47416E+06, 0.49778E+06,
     + 0.52233E+06, 0.54781E+06, 0.57427E+06, 0.60172E+06, 0.63019E+06,
     + 0.65971E+06, 0.69031E+06, 0.72201E+06, 0.75485E+06, 0.78886E+06,
     + 0.82405E+06, 0.86048E+06, 0.89815E+06, 0.93711E+06, 0.97739E+06,
     + 0.10190E+07, 0.10620E+07, 0.11065E+07, 0.11523E+07, 0.11997E+07,
     + 0.12485E+07, 0.12990E+07, 0.13510E+07, 0.14046E+07, 0.14599E+07,
     + 0.15169E+07, 0.15756E+07, 0.16361E+07, 0.16984E+07, 0.17626E+07,
     + 0.18287E+07, 0.18966E+07, 0.19666E+07, 0.20386E+07, 0.21126E+07,
     + 0.21887E+07, 0.22669E+07, 0.23474E+07, 0.24300E+07, 0.25150E+07,
     + 0.26022E+07, 0.26919E+07, 0.27839E+07, 0.28784E+07, 0.29753E+07,
     + 0.30749E+07/
c...        --       678
      data (QofT( 8,J),J=1,119)/ 0.39745E+04, 0.66993E+04, 0.98642E+04,
     + 0.13422E+05, 0.17352E+05, 0.21652E+05, 0.26339E+05, 0.31442E+05,
     + 0.37000E+05, 0.43058E+05, 0.49669E+05, 0.56885E+05, 0.64766E+05,
     + 0.73372E+05, 0.82765E+05, 0.93011E+05, 0.10418E+06, 0.11633E+06,
     + 0.12955E+06, 0.14390E+06, 0.15946E+06, 0.17632E+06, 0.19455E+06,
     + 0.21424E+06, 0.23547E+06, 0.25835E+06, 0.28296E+06, 0.30939E+06,
     + 0.33776E+06, 0.36816E+06, 0.40070E+06, 0.43549E+06, 0.47264E+06,
     + 0.51228E+06, 0.55451E+06, 0.59947E+06, 0.64728E+06, 0.69807E+06,
     + 0.75198E+06, 0.80915E+06, 0.86971E+06, 0.93381E+06, 0.10016E+07,
     + 0.10733E+07, 0.11489E+07, 0.12287E+07, 0.13128E+07, 0.14015E+07,
     + 0.14948E+07, 0.15930E+07, 0.16961E+07, 0.18045E+07, 0.19183E+07,
     + 0.20378E+07, 0.21629E+07, 0.22942E+07, 0.24316E+07, 0.25754E+07,
     + 0.27258E+07, 0.28831E+07, 0.30475E+07, 0.32192E+07, 0.33984E+07,
     + 0.35855E+07, 0.37805E+07, 0.39838E+07, 0.41956E+07, 0.44162E+07,
     + 0.46458E+07, 0.48847E+07, 0.51332E+07, 0.53916E+07, 0.56601E+07,
     + 0.59390E+07, 0.62286E+07, 0.65292E+07, 0.68412E+07, 0.71647E+07,
     + 0.75002E+07, 0.78479E+07, 0.82081E+07, 0.85813E+07, 0.89676E+07,
     + 0.93676E+07, 0.97814E+07, 0.10209E+08, 0.10652E+08, 0.11110E+08,
     + 0.11583E+08, 0.12071E+08, 0.12576E+08, 0.13097E+08, 0.13635E+08,
     + 0.14190E+08, 0.14763E+08, 0.15354E+08, 0.15963E+08, 0.16592E+08,
     + 0.17239E+08, 0.17906E+08, 0.18593E+08, 0.19301E+08, 0.20030E+08,
     + 0.20780E+08, 0.21553E+08, 0.22347E+08, 0.23165E+08, 0.24006E+08,
     + 0.24870E+08, 0.25759E+08, 0.26673E+08, 0.27612E+08, 0.28577E+08,
     + 0.29568E+08, 0.30585E+08, 0.31631E+08, 0.32704E+08, 0.33805E+08,
     + 0.34936E+08/
c...        --       768
      data (QofT( 9,J),J=1,119)/ 0.40228E+04, 0.67808E+04, 0.99842E+04,
     + 0.13586E+05, 0.17564E+05, 0.21919E+05, 0.26665E+05, 0.31833E+05,
     + 0.37461E+05, 0.43596E+05, 0.50286E+05, 0.57589E+05, 0.65562E+05,
     + 0.74264E+05, 0.83761E+05, 0.94115E+05, 0.10540E+06, 0.11767E+06,
     + 0.13102E+06, 0.14550E+06, 0.16121E+06, 0.17822E+06, 0.19661E+06,
     + 0.21646E+06, 0.23788E+06, 0.26094E+06, 0.28574E+06, 0.31239E+06,
     + 0.34097E+06, 0.37160E+06, 0.40437E+06, 0.43941E+06, 0.47683E+06,
     + 0.51673E+06, 0.55925E+06, 0.60451E+06, 0.65262E+06, 0.70374E+06,
     + 0.75799E+06, 0.81550E+06, 0.87643E+06, 0.94092E+06, 0.10091E+07,
     + 0.10812E+07, 0.11572E+07, 0.12375E+07, 0.13221E+07, 0.14112E+07,
     + 0.15050E+07, 0.16037E+07, 0.17074E+07, 0.18164E+07, 0.19307E+07,
     + 0.20507E+07, 0.21765E+07, 0.23084E+07, 0.24464E+07, 0.25909E+07,
     + 0.27421E+07, 0.29001E+07, 0.30652E+07, 0.32377E+07, 0.34177E+07,
     + 0.36055E+07, 0.38014E+07, 0.40055E+07, 0.42182E+07, 0.44397E+07,
     + 0.46703E+07, 0.49102E+07, 0.51597E+07, 0.54191E+07, 0.56886E+07,
     + 0.59686E+07, 0.62593E+07, 0.65611E+07, 0.68742E+07, 0.71989E+07,
     + 0.75356E+07, 0.78846E+07, 0.82461E+07, 0.86206E+07, 0.90083E+07,
     + 0.94097E+07, 0.98249E+07, 0.10254E+08, 0.10699E+08, 0.11158E+08,
     + 0.11632E+08, 0.12123E+08, 0.12629E+08, 0.13152E+08, 0.13691E+08,
     + 0.14248E+08, 0.14823E+08, 0.15416E+08, 0.16027E+08, 0.16657E+08,
     + 0.17307E+08, 0.17976E+08, 0.18665E+08, 0.19375E+08, 0.20106E+08,
     + 0.20858E+08, 0.21633E+08, 0.22430E+08, 0.23250E+08, 0.24093E+08,
     + 0.24960E+08, 0.25851E+08, 0.26767E+08, 0.27709E+08, 0.28676E+08,
     + 0.29670E+08, 0.30691E+08, 0.31739E+08, 0.32815E+08, 0.33919E+08,
     + 0.35053E+08/
c...        --       786
      data (QofT(10,J),J=1,119)/ 0.39315E+04, 0.66267E+04, 0.97569E+04,
     + 0.13276E+05, 0.17162E+05, 0.21414E+05, 0.26048E+05, 0.31094E+05,
     + 0.36590E+05, 0.42581E+05, 0.49120E+05, 0.56260E+05, 0.64061E+05,
     + 0.72580E+05, 0.81882E+05, 0.92031E+05, 0.10309E+06, 0.11514E+06,
     + 0.12824E+06, 0.14247E+06, 0.15791E+06, 0.17463E+06, 0.19272E+06,
     + 0.21226E+06, 0.23333E+06, 0.25604E+06, 0.28047E+06, 0.30673E+06,
     + 0.33490E+06, 0.36510E+06, 0.39743E+06, 0.43200E+06, 0.46892E+06,
     + 0.50831E+06, 0.55029E+06, 0.59498E+06, 0.64251E+06, 0.69301E+06,
     + 0.74662E+06, 0.80347E+06, 0.86370E+06, 0.92747E+06, 0.99491E+06,
     + 0.10662E+07, 0.11414E+07, 0.12208E+07, 0.13046E+07, 0.13928E+07,
     + 0.14856E+07, 0.15833E+07, 0.16860E+07, 0.17939E+07, 0.19072E+07,
     + 0.20261E+07, 0.21508E+07, 0.22814E+07, 0.24182E+07, 0.25614E+07,
     + 0.27112E+07, 0.28679E+07, 0.30316E+07, 0.32026E+07, 0.33811E+07,
     + 0.35674E+07, 0.37617E+07, 0.39642E+07, 0.41752E+07, 0.43950E+07,
     + 0.46237E+07, 0.48618E+07, 0.51094E+07, 0.53668E+07, 0.56343E+07,
     + 0.59123E+07, 0.62009E+07, 0.65005E+07, 0.68113E+07, 0.71338E+07,
     + 0.74681E+07, 0.78147E+07, 0.81737E+07, 0.85457E+07, 0.89308E+07,
     + 0.93295E+07, 0.97420E+07, 0.10169E+08, 0.10610E+08, 0.11066E+08,
     + 0.11538E+08, 0.12025E+08, 0.12528E+08, 0.13048E+08, 0.13584E+08,
     + 0.14138E+08, 0.14709E+08, 0.15298E+08, 0.15906E+08, 0.16532E+08,
     + 0.17178E+08, 0.17843E+08, 0.18528E+08, 0.19234E+08, 0.19961E+08,
     + 0.20710E+08, 0.21480E+08, 0.22272E+08, 0.23088E+08, 0.23926E+08,
     + 0.24789E+08, 0.25675E+08, 0.26587E+08, 0.27523E+08, 0.28485E+08,
     + 0.29474E+08, 0.30489E+08, 0.31532E+08, 0.32603E+08, 0.33701E+08,
     + 0.34829E+08/
c...        --       776
      data (QofT(11,J),J=1,119)/ 0.23106E+05, 0.38945E+05, 0.57342E+05,
     + 0.78021E+05, 0.10085E+06, 0.12582E+06, 0.15302E+06, 0.18262E+06,
     + 0.21482E+06, 0.24989E+06, 0.28812E+06, 0.32983E+06, 0.37535E+06,
     + 0.42501E+06, 0.47919E+06, 0.53825E+06, 0.60258E+06, 0.67256E+06,
     + 0.74862E+06, 0.83118E+06, 0.92069E+06, 0.10176E+07, 0.11223E+07,
     + 0.12354E+07, 0.13574E+07, 0.14887E+07, 0.16299E+07, 0.17816E+07,
     + 0.19443E+07, 0.21187E+07, 0.23052E+07, 0.25047E+07, 0.27176E+07,
     + 0.29447E+07, 0.31866E+07, 0.34441E+07, 0.37179E+07, 0.40087E+07,
     + 0.43173E+07, 0.46444E+07, 0.49910E+07, 0.53578E+07, 0.57456E+07,
     + 0.61554E+07, 0.65880E+07, 0.70444E+07, 0.75255E+07, 0.80322E+07,
     + 0.85656E+07, 0.91266E+07, 0.97163E+07, 0.10336E+08, 0.10986E+08,
     + 0.11668E+08, 0.12383E+08, 0.13133E+08, 0.13918E+08, 0.14739E+08,
     + 0.15598E+08, 0.16496E+08, 0.17435E+08, 0.18415E+08, 0.19438E+08,
     + 0.20505E+08, 0.21619E+08, 0.22779E+08, 0.23987E+08, 0.25246E+08,
     + 0.26556E+08, 0.27920E+08, 0.29337E+08, 0.30811E+08, 0.32343E+08,
     + 0.33934E+08, 0.35585E+08, 0.37300E+08, 0.39079E+08, 0.40924E+08,
     + 0.42837E+08, 0.44819E+08, 0.46873E+08, 0.49001E+08, 0.51203E+08,
     + 0.53483E+08, 0.55842E+08, 0.58282E+08, 0.60805E+08, 0.63414E+08,
     + 0.66109E+08, 0.68894E+08, 0.71770E+08, 0.74740E+08, 0.77806E+08,
     + 0.80970E+08, 0.84234E+08, 0.87600E+08, 0.91072E+08, 0.94651E+08,
     + 0.98339E+08, 0.10214E+09, 0.10605E+09, 0.11009E+09, 0.11424E+09,
     + 0.11851E+09, 0.12291E+09, 0.12744E+09, 0.13209E+09, 0.13688E+09,
     + 0.14180E+09, 0.14687E+09, 0.15207E+09, 0.15742E+09, 0.16291E+09,
     + 0.16855E+09, 0.17435E+09, 0.18030E+09, 0.18641E+09, 0.19268E+09,
     + 0.19912E+09/
c...        --       767
      data (QofT(12,J),J=1,119)/ 0.11692E+05, 0.19707E+05, 0.29017E+05,
     + 0.39482E+05, 0.51038E+05, 0.63680E+05, 0.77450E+05, 0.92432E+05,
     + 0.10873E+06, 0.12649E+06, 0.14584E+06, 0.16694E+06, 0.18996E+06,
     + 0.21507E+06, 0.24245E+06, 0.27229E+06, 0.30478E+06, 0.34013E+06,
     + 0.37853E+06, 0.42020E+06, 0.46536E+06, 0.51424E+06, 0.56708E+06,
     + 0.62411E+06, 0.68559E+06, 0.75178E+06, 0.82296E+06, 0.89939E+06,
     + 0.98137E+06, 0.10692E+07, 0.11631E+07, 0.12636E+07, 0.13708E+07,
     + 0.14851E+07, 0.16069E+07, 0.17365E+07, 0.18742E+07, 0.20206E+07,
     + 0.21758E+07, 0.23404E+07, 0.25148E+07, 0.26992E+07, 0.28943E+07,
     + 0.31004E+07, 0.33179E+07, 0.35474E+07, 0.37892E+07, 0.40440E+07,
     + 0.43121E+07, 0.45940E+07, 0.48904E+07, 0.52017E+07, 0.55285E+07,
     + 0.58713E+07, 0.62306E+07, 0.66071E+07, 0.70014E+07, 0.74140E+07,
     + 0.78456E+07, 0.82967E+07, 0.87681E+07, 0.92604E+07, 0.97742E+07,
     + 0.10310E+08, 0.10869E+08, 0.11452E+08, 0.12059E+08, 0.12691E+08,
     + 0.13348E+08, 0.14033E+08, 0.14745E+08, 0.15484E+08, 0.16253E+08,
     + 0.17052E+08, 0.17881E+08, 0.18741E+08, 0.19634E+08, 0.20560E+08,
     + 0.21520E+08, 0.22515E+08, 0.23546E+08, 0.24613E+08, 0.25718E+08,
     + 0.26862E+08, 0.28046E+08, 0.29270E+08, 0.30536E+08, 0.31845E+08,
     + 0.33197E+08, 0.34594E+08, 0.36037E+08, 0.37527E+08, 0.39065E+08,
     + 0.40652E+08, 0.42289E+08, 0.43977E+08, 0.45719E+08, 0.47514E+08,
     + 0.49363E+08, 0.51270E+08, 0.53233E+08, 0.55255E+08, 0.57337E+08,
     + 0.59480E+08, 0.61686E+08, 0.63956E+08, 0.66290E+08, 0.68691E+08,
     + 0.71160E+08, 0.73699E+08, 0.76307E+08, 0.78988E+08, 0.81743E+08,
     + 0.84572E+08, 0.87478E+08, 0.90462E+08, 0.93525E+08, 0.96669E+08,
     + 0.99896E+08/
c...        --       888
      data (QofT(13,J),J=1,119)/ 0.36175E+03, 0.60978E+03, 0.89790E+03,
     + 0.12219E+04, 0.15802E+04, 0.19728E+04, 0.24016E+04, 0.28696E+04,
     + 0.33807E+04, 0.39394E+04, 0.45506E+04, 0.52196E+04, 0.59521E+04,
     + 0.67538E+04, 0.76308E+04, 0.85894E+04, 0.96361E+04, 0.10777E+05,
     + 0.12021E+05, 0.13373E+05, 0.14841E+05, 0.16434E+05, 0.18158E+05,
     + 0.20023E+05, 0.22037E+05, 0.24208E+05, 0.26547E+05, 0.29061E+05,
     + 0.31762E+05, 0.34659E+05, 0.37762E+05, 0.41083E+05, 0.44632E+05,
     + 0.48421E+05, 0.52462E+05, 0.56766E+05, 0.61346E+05, 0.66215E+05,
     + 0.71386E+05, 0.76873E+05, 0.82688E+05, 0.88848E+05, 0.95365E+05,
     + 0.10226E+06, 0.10954E+06, 0.11722E+06, 0.12532E+06, 0.13387E+06,
     + 0.14286E+06, 0.15233E+06, 0.16229E+06, 0.17275E+06, 0.18374E+06,
     + 0.19528E+06, 0.20737E+06, 0.22006E+06, 0.23335E+06, 0.24726E+06,
     + 0.26182E+06, 0.27705E+06, 0.29297E+06, 0.30960E+06, 0.32696E+06,
     + 0.34509E+06, 0.36399E+06, 0.38371E+06, 0.40425E+06, 0.42566E+06,
     + 0.44794E+06, 0.47114E+06, 0.49527E+06, 0.52036E+06, 0.54644E+06,
     + 0.57354E+06, 0.60169E+06, 0.63091E+06, 0.66124E+06, 0.69270E+06,
     + 0.72533E+06, 0.75916E+06, 0.79421E+06, 0.83053E+06, 0.86814E+06,
     + 0.90708E+06, 0.94737E+06, 0.98907E+06, 0.10322E+07, 0.10768E+07,
     + 0.11229E+07, 0.11705E+07, 0.12197E+07, 0.12705E+07, 0.13230E+07,
     + 0.13771E+07, 0.14330E+07, 0.14906E+07, 0.15501E+07, 0.16114E+07,
     + 0.16745E+07, 0.17397E+07, 0.18067E+07, 0.18759E+07, 0.19470E+07,
     + 0.20203E+07, 0.20957E+07, 0.21733E+07, 0.22532E+07, 0.23353E+07,
     + 0.24198E+07, 0.25067E+07, 0.25960E+07, 0.26878E+07, 0.27821E+07,
     + 0.28790E+07, 0.29785E+07, 0.30807E+07, 0.31857E+07, 0.32934E+07,
     + 0.34040E+07/
c...        --       887
      data (QofT(14,J),J=1,119)/ 0.42000E+04, 0.70796E+04, 0.10424E+05,
     + 0.14186E+05, 0.18342E+05, 0.22896E+05, 0.27866E+05, 0.33285E+05,
     + 0.39199E+05, 0.45659E+05, 0.52720E+05, 0.60444E+05, 0.68895E+05,
     + 0.78139E+05, 0.88246E+05, 0.99288E+05, 0.11134E+06, 0.12447E+06,
     + 0.13877E+06, 0.15431E+06, 0.17119E+06, 0.18949E+06, 0.20930E+06,
     + 0.23071E+06, 0.25383E+06, 0.27875E+06, 0.30558E+06, 0.33442E+06,
     + 0.36539E+06, 0.39861E+06, 0.43418E+06, 0.47224E+06, 0.51291E+06,
     + 0.55632E+06, 0.60260E+06, 0.65189E+06, 0.70434E+06, 0.76008E+06,
     + 0.81927E+06, 0.88206E+06, 0.94862E+06, 0.10191E+07, 0.10937E+07,
     + 0.11725E+07, 0.12558E+07, 0.13436E+07, 0.14363E+07, 0.15340E+07,
     + 0.16368E+07, 0.17450E+07, 0.18588E+07, 0.19784E+07, 0.21040E+07,
     + 0.22358E+07, 0.23741E+07, 0.25190E+07, 0.26708E+07, 0.28297E+07,
     + 0.29961E+07, 0.31700E+07, 0.33518E+07, 0.35417E+07, 0.37400E+07,
     + 0.39469E+07, 0.41628E+07, 0.43878E+07, 0.46224E+07, 0.48667E+07,
     + 0.51210E+07, 0.53858E+07, 0.56611E+07, 0.59475E+07, 0.62451E+07,
     + 0.65544E+07, 0.68755E+07, 0.72089E+07, 0.75550E+07, 0.79139E+07,
     + 0.82861E+07, 0.86720E+07, 0.90719E+07, 0.94861E+07, 0.99151E+07,
     + 0.10359E+08, 0.10819E+08, 0.11294E+08, 0.11786E+08, 0.12294E+08,
     + 0.12820E+08, 0.13363E+08, 0.13924E+08, 0.14503E+08, 0.15101E+08,
     + 0.15719E+08, 0.16356E+08, 0.17013E+08, 0.17690E+08, 0.18389E+08,
     + 0.19109E+08, 0.19851E+08, 0.20616E+08, 0.21404E+08, 0.22215E+08,
     + 0.23050E+08, 0.23910E+08, 0.24794E+08, 0.25704E+08, 0.26640E+08,
     + 0.27603E+08, 0.28593E+08, 0.29610E+08, 0.30656E+08, 0.31731E+08,
     + 0.32835E+08, 0.33969E+08, 0.35133E+08, 0.36329E+08, 0.37556E+08,
     + 0.38816E+08/
c...        --       878
      data (QofT(15,J),J=1,119)/ 0.21250E+04, 0.35820E+04, 0.52744E+04,
     + 0.71778E+04, 0.92814E+04, 0.11586E+05, 0.14102E+05, 0.16845E+05,
     + 0.19839E+05, 0.23108E+05, 0.26680E+05, 0.30588E+05, 0.34861E+05,
     + 0.39534E+05, 0.44642E+05, 0.50219E+05, 0.56305E+05, 0.62937E+05,
     + 0.70155E+05, 0.78001E+05, 0.86516E+05, 0.95747E+05, 0.10574E+06,
     + 0.11653E+06, 0.12819E+06, 0.14075E+06, 0.15427E+06, 0.16881E+06,
     + 0.18441E+06, 0.20114E+06, 0.21906E+06, 0.23823E+06, 0.25871E+06,
     + 0.28056E+06, 0.30386E+06, 0.32867E+06, 0.35507E+06, 0.38312E+06,
     + 0.41291E+06, 0.44450E+06, 0.47799E+06, 0.51344E+06, 0.55095E+06,
     + 0.59060E+06, 0.63248E+06, 0.67667E+06, 0.72327E+06, 0.77238E+06,
     + 0.82409E+06, 0.87850E+06, 0.93571E+06, 0.99583E+06, 0.10590E+07,
     + 0.11252E+07, 0.11947E+07, 0.12675E+07, 0.13438E+07, 0.14237E+07,
     + 0.15072E+07, 0.15946E+07, 0.16859E+07, 0.17814E+07, 0.18810E+07,
     + 0.19849E+07, 0.20934E+07, 0.22064E+07, 0.23242E+07, 0.24469E+07,
     + 0.25747E+07, 0.27076E+07, 0.28459E+07, 0.29897E+07, 0.31391E+07,
     + 0.32944E+07, 0.34557E+07, 0.36231E+07, 0.37968E+07, 0.39770E+07,
     + 0.41639E+07, 0.43576E+07, 0.45583E+07, 0.47663E+07, 0.49816E+07,
     + 0.52045E+07, 0.54352E+07, 0.56739E+07, 0.59207E+07, 0.61759E+07,
     + 0.64396E+07, 0.67121E+07, 0.69936E+07, 0.72844E+07, 0.75845E+07,
     + 0.78943E+07, 0.82139E+07, 0.85436E+07, 0.88837E+07, 0.92342E+07,
     + 0.95956E+07, 0.99680E+07, 0.10352E+08, 0.10747E+08, 0.11154E+08,
     + 0.11573E+08, 0.12004E+08, 0.12448E+08, 0.12904E+08, 0.13374E+08,
     + 0.13857E+08, 0.14353E+08, 0.14864E+08, 0.15388E+08, 0.15927E+08,
     + 0.16481E+08, 0.17050E+08, 0.17634E+08, 0.18234E+08, 0.18849E+08,
     + 0.19481E+08/
c...        --       778
      data (QofT(16,J),J=1,119)/ 0.24692E+05, 0.41621E+05, 0.61284E+05,
     + 0.83394E+05, 0.10782E+06, 0.13457E+06, 0.16375E+06, 0.19554E+06,
     + 0.23020E+06, 0.26801E+06, 0.30930E+06, 0.35443E+06, 0.40375E+06,
     + 0.45763E+06, 0.51650E+06, 0.58075E+06, 0.65080E+06, 0.72711E+06,
     + 0.81012E+06, 0.90030E+06, 0.99815E+06, 0.11042E+07, 0.12189E+07,
     + 0.13428E+07, 0.14765E+07, 0.16206E+07, 0.17757E+07, 0.19423E+07,
     + 0.21212E+07, 0.23129E+07, 0.25181E+07, 0.27377E+07, 0.29721E+07,
     + 0.32223E+07, 0.34890E+07, 0.37729E+07, 0.40750E+07, 0.43959E+07,
     + 0.47365E+07, 0.50978E+07, 0.54807E+07, 0.58860E+07, 0.63147E+07,
     + 0.67678E+07, 0.72463E+07, 0.77512E+07, 0.82836E+07, 0.88445E+07,
     + 0.94351E+07, 0.10056E+08, 0.10710E+08, 0.11396E+08, 0.12117E+08,
     + 0.12873E+08, 0.13666E+08, 0.14497E+08, 0.15367E+08, 0.16279E+08,
     + 0.17232E+08, 0.18229E+08, 0.19271E+08, 0.20359E+08, 0.21495E+08,
     + 0.22681E+08, 0.23917E+08, 0.25206E+08, 0.26549E+08, 0.27948E+08,
     + 0.29404E+08, 0.30920E+08, 0.32496E+08, 0.34135E+08, 0.35838E+08,
     + 0.37608E+08, 0.39445E+08, 0.41353E+08, 0.43332E+08, 0.45385E+08,
     + 0.47514E+08, 0.49721E+08, 0.52007E+08, 0.54376E+08, 0.56829E+08,
     + 0.59367E+08, 0.61995E+08, 0.64712E+08, 0.67523E+08, 0.70429E+08,
     + 0.73432E+08, 0.76535E+08, 0.79740E+08, 0.83050E+08, 0.86467E+08,
     + 0.89993E+08, 0.93632E+08, 0.97385E+08, 0.10126E+09, 0.10525E+09,
     + 0.10936E+09, 0.11360E+09, 0.11796E+09, 0.12246E+09, 0.12709E+09,
     + 0.13186E+09, 0.13677E+09, 0.14182E+09, 0.14701E+09, 0.15236E+09,
     + 0.15785E+09, 0.16350E+09, 0.16931E+09, 0.17528E+09, 0.18141E+09,
     + 0.18771E+09, 0.19418E+09, 0.20082E+09, 0.20764E+09, 0.21465E+09,
     + 0.22183E+09/
c...        --       787
      data (QofT(17,J),J=1,119)/ 0.12211E+05, 0.20582E+05, 0.30305E+05,
     + 0.41237E+05, 0.53314E+05, 0.66536E+05, 0.80957E+05, 0.96672E+05,
     + 0.11380E+06, 0.13250E+06, 0.15292E+06, 0.17524E+06, 0.19965E+06,
     + 0.22632E+06, 0.25546E+06, 0.28728E+06, 0.32199E+06, 0.35980E+06,
     + 0.40094E+06, 0.44565E+06, 0.49417E+06, 0.54676E+06, 0.60366E+06,
     + 0.66516E+06, 0.73152E+06, 0.80305E+06, 0.88002E+06, 0.96276E+06,
     + 0.10516E+07, 0.11468E+07, 0.12488E+07, 0.13578E+07, 0.14743E+07,
     + 0.15987E+07, 0.17312E+07, 0.18723E+07, 0.20225E+07, 0.21820E+07,
     + 0.23514E+07, 0.25310E+07, 0.27214E+07, 0.29230E+07, 0.31362E+07,
     + 0.33616E+07, 0.35997E+07, 0.38509E+07, 0.41158E+07, 0.43949E+07,
     + 0.46887E+07, 0.49980E+07, 0.53231E+07, 0.56647E+07, 0.60234E+07,
     + 0.63998E+07, 0.67946E+07, 0.72084E+07, 0.76418E+07, 0.80955E+07,
     + 0.85702E+07, 0.90666E+07, 0.95854E+07, 0.10127E+08, 0.10693E+08,
     + 0.11284E+08, 0.11900E+08, 0.12542E+08, 0.13211E+08, 0.13907E+08,
     + 0.14633E+08, 0.15388E+08, 0.16173E+08, 0.16990E+08, 0.17838E+08,
     + 0.18720E+08, 0.19636E+08, 0.20586E+08, 0.21573E+08, 0.22596E+08,
     + 0.23657E+08, 0.24757E+08, 0.25896E+08, 0.27077E+08, 0.28299E+08,
     + 0.29565E+08, 0.30874E+08, 0.32229E+08, 0.33630E+08, 0.35079E+08,
     + 0.36576E+08, 0.38123E+08, 0.39721E+08, 0.41371E+08, 0.43075E+08,
     + 0.44833E+08, 0.46647E+08, 0.48518E+08, 0.50448E+08, 0.52438E+08,
     + 0.54489E+08, 0.56603E+08, 0.58780E+08, 0.61023E+08, 0.63332E+08,
     + 0.65710E+08, 0.68157E+08, 0.70676E+08, 0.73266E+08, 0.75931E+08,
     + 0.78672E+08, 0.81490E+08, 0.84386E+08, 0.87363E+08, 0.90422E+08,
     + 0.93564E+08, 0.96791E+08, 0.10011E+09, 0.10351E+09, 0.10700E+09,
     + 0.11059E+09/
c...        --       777
      data (QofT(18,J),J=1,119)/ 0.71750E+05, 0.12094E+06, 0.17807E+06,
     + 0.24230E+06, 0.31324E+06, 0.39088E+06, 0.47550E+06, 0.56764E+06,
     + 0.66800E+06, 0.77740E+06, 0.89677E+06, 0.10271E+07, 0.11694E+07,
     + 0.13249E+07, 0.14945E+07, 0.16796E+07, 0.18813E+07, 0.21009E+07,
     + 0.23396E+07, 0.25989E+07, 0.28801E+07, 0.31847E+07, 0.35140E+07,
     + 0.38698E+07, 0.42535E+07, 0.46669E+07, 0.51115E+07, 0.55893E+07,
     + 0.61019E+07, 0.66513E+07, 0.72393E+07, 0.78680E+07, 0.85395E+07,
     + 0.92558E+07, 0.10019E+08, 0.10832E+08, 0.11696E+08, 0.12614E+08,
     + 0.13588E+08, 0.14621E+08, 0.15716E+08, 0.16875E+08, 0.18100E+08,
     + 0.19395E+08, 0.20762E+08, 0.22205E+08, 0.23726E+08, 0.25328E+08,
     + 0.27015E+08, 0.28789E+08, 0.30654E+08, 0.32614E+08, 0.34671E+08,
     + 0.36830E+08, 0.39093E+08, 0.41465E+08, 0.43949E+08, 0.46549E+08,
     + 0.49269E+08, 0.52112E+08, 0.55084E+08, 0.58188E+08, 0.61428E+08,
     + 0.64809E+08, 0.68335E+08, 0.72010E+08, 0.75840E+08, 0.79828E+08,
     + 0.83979E+08, 0.88299E+08, 0.92792E+08, 0.97463E+08, 0.10232E+09,
     + 0.10736E+09, 0.11260E+09, 0.11803E+09, 0.12367E+09, 0.12952E+09,
     + 0.13559E+09, 0.14187E+09, 0.14839E+09, 0.15513E+09, 0.16212E+09,
     + 0.16935E+09, 0.17683E+09, 0.18457E+09, 0.19257E+09, 0.20085E+09,
     + 0.20940E+09, 0.21824E+09, 0.22736E+09, 0.23678E+09, 0.24651E+09,
     + 0.25655E+09, 0.26691E+09, 0.27759E+09, 0.28861E+09, 0.29997E+09,
     + 0.31167E+09, 0.32374E+09, 0.33616E+09, 0.34896E+09, 0.36214E+09,
     + 0.37571E+09, 0.38967E+09, 0.40404E+09, 0.41882E+09, 0.43403E+09,
     + 0.44966E+09, 0.46573E+09, 0.48226E+09, 0.49923E+09, 0.51668E+09,
     + 0.53460E+09, 0.55301E+09, 0.57191E+09, 0.59131E+09, 0.61123E+09,
     + 0.63167E+09/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_N2O   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 5), QofT( 5,119),Q(NT)
      data xgj/ 9.,6.,6.,9.,54. /
c...      N2O
c...        --       446
      data (QofT( 1,J),J=1,119)/ 0.89943E+03, 0.12734E+04, 0.16489E+04,
     + 0.20293E+04, 0.24205E+04, 0.28289E+04, 0.32609E+04, 0.37222E+04,
     + 0.42180E+04, 0.47529E+04, 0.53312E+04, 0.59572E+04, 0.66348E+04,
     + 0.73683E+04, 0.81616E+04, 0.90190E+04, 0.99450E+04, 0.10944E+05,
     + 0.12021E+05, 0.13180E+05, 0.14426E+05, 0.15766E+05, 0.17203E+05,
     + 0.18745E+05, 0.20396E+05, 0.22162E+05, 0.24051E+05, 0.26069E+05,
     + 0.28222E+05, 0.30517E+05, 0.32962E+05, 0.35564E+05, 0.38331E+05,
     + 0.41271E+05, 0.44393E+05, 0.47704E+05, 0.51214E+05, 0.54932E+05,
     + 0.58868E+05, 0.63030E+05, 0.67429E+05, 0.72075E+05, 0.76979E+05,
     + 0.82151E+05, 0.87604E+05, 0.93348E+05, 0.99395E+05, 0.10576E+06,
     + 0.11245E+06, 0.11948E+06, 0.12686E+06, 0.13461E+06, 0.14275E+06,
     + 0.15128E+06, 0.16021E+06, 0.16958E+06, 0.17938E+06, 0.18964E+06,
     + 0.20037E+06, 0.21159E+06, 0.22331E+06, 0.23556E+06, 0.24834E+06,
     + 0.26169E+06, 0.27561E+06, 0.29012E+06, 0.30525E+06, 0.32101E+06,
     + 0.33743E+06, 0.35452E+06, 0.37230E+06, 0.39080E+06, 0.41004E+06,
     + 0.43004E+06, 0.45082E+06, 0.47241E+06, 0.49483E+06, 0.51810E+06,
     + 0.54225E+06, 0.56730E+06, 0.59329E+06, 0.62022E+06, 0.64814E+06,
     + 0.67707E+06, 0.70703E+06, 0.73806E+06, 0.77018E+06, 0.80342E+06,
     + 0.83781E+06, 0.87338E+06, 0.91016E+06, 0.94818E+06, 0.98748E+06,
     + 0.10281E+07, 0.10700E+07, 0.11133E+07, 0.11581E+07, 0.12042E+07,
     + 0.12519E+07, 0.13010E+07, 0.13517E+07, 0.14040E+07, 0.14579E+07,
     + 0.15134E+07, 0.15707E+07, 0.16297E+07, 0.16905E+07, 0.17530E+07,
     + 0.18175E+07, 0.18838E+07, 0.19521E+07, 0.20224E+07, 0.20947E+07,
     + 0.21690E+07, 0.22455E+07, 0.23242E+07, 0.24050E+07, 0.24881E+07,
     + 0.25735E+07/
c...        --       456
      data (QofT( 2,J),J=1,119)/ 0.59966E+03, 0.84903E+03, 0.10995E+04,
     + 0.13538E+04, 0.16158E+04, 0.18903E+04, 0.21815E+04, 0.24934E+04,
     + 0.28295E+04, 0.31927E+04, 0.35862E+04, 0.40128E+04, 0.44752E+04,
     + 0.49763E+04, 0.55189E+04, 0.61059E+04, 0.67404E+04, 0.74256E+04,
     + 0.81646E+04, 0.89609E+04, 0.98180E+04, 0.10740E+05, 0.11729E+05,
     + 0.12791E+05, 0.13930E+05, 0.15149E+05, 0.16453E+05, 0.17847E+05,
     + 0.19335E+05, 0.20922E+05, 0.22614E+05, 0.24416E+05, 0.26333E+05,
     + 0.28371E+05, 0.30535E+05, 0.32833E+05, 0.35269E+05, 0.37851E+05,
     + 0.40585E+05, 0.43478E+05, 0.46537E+05, 0.49769E+05, 0.53182E+05,
     + 0.56783E+05, 0.60580E+05, 0.64582E+05, 0.68796E+05, 0.73232E+05,
     + 0.77898E+05, 0.82803E+05, 0.87957E+05, 0.93369E+05, 0.99048E+05,
     + 0.10501E+06, 0.11125E+06, 0.11780E+06, 0.12465E+06, 0.13182E+06,
     + 0.13933E+06, 0.14718E+06, 0.15539E+06, 0.16396E+06, 0.17291E+06,
     + 0.18226E+06, 0.19201E+06, 0.20218E+06, 0.21278E+06, 0.22383E+06,
     + 0.23534E+06, 0.24733E+06, 0.25980E+06, 0.27278E+06, 0.28628E+06,
     + 0.30032E+06, 0.31491E+06, 0.33007E+06, 0.34581E+06, 0.36216E+06,
     + 0.37912E+06, 0.39673E+06, 0.41499E+06, 0.43392E+06, 0.45355E+06,
     + 0.47389E+06, 0.49496E+06, 0.51678E+06, 0.53937E+06, 0.56276E+06,
     + 0.58695E+06, 0.61199E+06, 0.63788E+06, 0.66464E+06, 0.69231E+06,
     + 0.72090E+06, 0.75044E+06, 0.78094E+06, 0.81244E+06, 0.84496E+06,
     + 0.87853E+06, 0.91316E+06, 0.94889E+06, 0.98573E+06, 0.10237E+07,
     + 0.10629E+07, 0.11033E+07, 0.11449E+07, 0.11877E+07, 0.12319E+07,
     + 0.12773E+07, 0.13241E+07, 0.13723E+07, 0.14219E+07, 0.14729E+07,
     + 0.15254E+07, 0.15793E+07, 0.16349E+07, 0.16919E+07, 0.17506E+07,
     + 0.18109E+07/
c...        --       546
      data (QofT( 3,J),J=1,119)/ 0.62051E+03, 0.87856E+03, 0.11377E+04,
     + 0.14003E+04, 0.16705E+04, 0.19529E+04, 0.22518E+04, 0.25713E+04,
     + 0.29149E+04, 0.32859E+04, 0.36873E+04, 0.41220E+04, 0.45929E+04,
     + 0.51028E+04, 0.56547E+04, 0.62515E+04, 0.68963E+04, 0.75923E+04,
     + 0.83428E+04, 0.91511E+04, 0.10021E+05, 0.10956E+05, 0.11960E+05,
     + 0.13036E+05, 0.14190E+05, 0.15425E+05, 0.16746E+05, 0.18158E+05,
     + 0.19664E+05, 0.21271E+05, 0.22984E+05, 0.24806E+05, 0.26745E+05,
     + 0.28806E+05, 0.30995E+05, 0.33317E+05, 0.35780E+05, 0.38389E+05,
     + 0.41151E+05, 0.44073E+05, 0.47162E+05, 0.50425E+05, 0.53871E+05,
     + 0.57505E+05, 0.61338E+05, 0.65375E+05, 0.69628E+05, 0.74102E+05,
     + 0.78808E+05, 0.83755E+05, 0.88951E+05, 0.94407E+05, 0.10013E+06,
     + 0.10614E+06, 0.11243E+06, 0.11902E+06, 0.12593E+06, 0.13316E+06,
     + 0.14072E+06, 0.14862E+06, 0.15689E+06, 0.16552E+06, 0.17453E+06,
     + 0.18394E+06, 0.19376E+06, 0.20399E+06, 0.21466E+06, 0.22578E+06,
     + 0.23737E+06, 0.24942E+06, 0.26198E+06, 0.27503E+06, 0.28861E+06,
     + 0.30273E+06, 0.31741E+06, 0.33265E+06, 0.34848E+06, 0.36492E+06,
     + 0.38197E+06, 0.39967E+06, 0.41803E+06, 0.43706E+06, 0.45679E+06,
     + 0.47723E+06, 0.49840E+06, 0.52033E+06, 0.54303E+06, 0.56653E+06,
     + 0.59084E+06, 0.61599E+06, 0.64200E+06, 0.66888E+06, 0.69667E+06,
     + 0.72539E+06, 0.75506E+06, 0.78569E+06, 0.81733E+06, 0.84998E+06,
     + 0.88369E+06, 0.91846E+06, 0.95433E+06, 0.99132E+06, 0.10295E+07,
     + 0.10688E+07, 0.11093E+07, 0.11511E+07, 0.11941E+07, 0.12384E+07,
     + 0.12840E+07, 0.13310E+07, 0.13793E+07, 0.14291E+07, 0.14803E+07,
     + 0.15329E+07, 0.15871E+07, 0.16428E+07, 0.17000E+07, 0.17589E+07,
     + 0.18194E+07/
c...        --       448
      data (QofT( 4,J),J=1,119)/ 0.95253E+03, 0.13487E+04, 0.17465E+04,
     + 0.21498E+04, 0.25648E+04, 0.29986E+04, 0.34580E+04, 0.39493E+04,
     + 0.44779E+04, 0.50488E+04, 0.56669E+04, 0.63366E+04, 0.70625E+04,
     + 0.78488E+04, 0.87003E+04, 0.96216E+04, 0.10617E+05, 0.11692E+05,
     + 0.12852E+05, 0.14102E+05, 0.15447E+05, 0.16893E+05, 0.18446E+05,
     + 0.20112E+05, 0.21898E+05, 0.23811E+05, 0.25856E+05, 0.28042E+05,
     + 0.30377E+05, 0.32866E+05, 0.35520E+05, 0.38345E+05, 0.41351E+05,
     + 0.44545E+05, 0.47939E+05, 0.51540E+05, 0.55359E+05, 0.59405E+05,
     + 0.63689E+05, 0.68222E+05, 0.73015E+05, 0.78078E+05, 0.83424E+05,
     + 0.89064E+05, 0.95012E+05, 0.10128E+06, 0.10788E+06, 0.11482E+06,
     + 0.12213E+06, 0.12981E+06, 0.13788E+06, 0.14635E+06, 0.15524E+06,
     + 0.16456E+06, 0.17433E+06, 0.18457E+06, 0.19530E+06, 0.20652E+06,
     + 0.21827E+06, 0.23055E+06, 0.24338E+06, 0.25679E+06, 0.27079E+06,
     + 0.28541E+06, 0.30066E+06, 0.31656E+06, 0.33314E+06, 0.35042E+06,
     + 0.36841E+06, 0.38715E+06, 0.40666E+06, 0.42695E+06, 0.44805E+06,
     + 0.46999E+06, 0.49279E+06, 0.51649E+06, 0.54109E+06, 0.56664E+06,
     + 0.59315E+06, 0.62066E+06, 0.64919E+06, 0.67877E+06, 0.70943E+06,
     + 0.74121E+06, 0.77413E+06, 0.80822E+06, 0.84351E+06, 0.88004E+06,
     + 0.91783E+06, 0.95693E+06, 0.99737E+06, 0.10392E+07, 0.10824E+07,
     + 0.11270E+07, 0.11732E+07, 0.12208E+07, 0.12700E+07, 0.13208E+07,
     + 0.13732E+07, 0.14272E+07, 0.14830E+07, 0.15405E+07, 0.15999E+07,
     + 0.16610E+07, 0.17240E+07, 0.17890E+07, 0.18559E+07, 0.19248E+07,
     + 0.19957E+07, 0.20687E+07, 0.21439E+07, 0.22213E+07, 0.23009E+07,
     + 0.23828E+07, 0.24671E+07, 0.25537E+07, 0.26428E+07, 0.27343E+07,
     + 0.28284E+07/
c...        --       447
      data (QofT( 5,J),J=1,119)/ 0.55598E+04, 0.78718E+04, 0.10193E+05,
     + 0.12546E+05, 0.14966E+05, 0.17495E+05, 0.20171E+05, 0.23031E+05,
     + 0.26106E+05, 0.29426E+05, 0.33018E+05, 0.36908E+05, 0.41121E+05,
     + 0.45684E+05, 0.50622E+05, 0.55962E+05, 0.61731E+05, 0.67958E+05,
     + 0.74671E+05, 0.81902E+05, 0.89681E+05, 0.98043E+05, 0.10702E+06,
     + 0.11665E+06, 0.12697E+06, 0.13801E+06, 0.14983E+06, 0.16244E+06,
     + 0.17591E+06, 0.19028E+06, 0.20558E+06, 0.22188E+06, 0.23920E+06,
     + 0.25762E+06, 0.27718E+06, 0.29793E+06, 0.31993E+06, 0.34323E+06,
     + 0.36791E+06, 0.39401E+06, 0.42160E+06, 0.45074E+06, 0.48151E+06,
     + 0.51397E+06, 0.54819E+06, 0.58424E+06, 0.62221E+06, 0.66215E+06,
     + 0.70416E+06, 0.74832E+06, 0.79470E+06, 0.84340E+06, 0.89450E+06,
     + 0.94808E+06, 0.10042E+07, 0.10631E+07, 0.11247E+07, 0.11892E+07,
     + 0.12567E+07, 0.13272E+07, 0.14009E+07, 0.14779E+07, 0.15583E+07,
     + 0.16422E+07, 0.17298E+07, 0.18211E+07, 0.19163E+07, 0.20154E+07,
     + 0.21187E+07, 0.22263E+07, 0.23382E+07, 0.24546E+07, 0.25757E+07,
     + 0.27016E+07, 0.28324E+07, 0.29683E+07, 0.31095E+07, 0.32560E+07,
     + 0.34081E+07, 0.35659E+07, 0.37295E+07, 0.38991E+07, 0.40750E+07,
     + 0.42572E+07, 0.44459E+07, 0.46414E+07, 0.48437E+07, 0.50531E+07,
     + 0.52698E+07, 0.54939E+07, 0.57257E+07, 0.59653E+07, 0.62129E+07,
     + 0.64688E+07, 0.67331E+07, 0.70061E+07, 0.72880E+07, 0.75790E+07,
     + 0.78792E+07, 0.81891E+07, 0.85086E+07, 0.88382E+07, 0.91780E+07,
     + 0.95283E+07, 0.98893E+07, 0.10261E+08, 0.10644E+08, 0.11039E+08,
     + 0.11445E+08, 0.11864E+08, 0.12294E+08, 0.12738E+08, 0.13194E+08,
     + 0.13663E+08, 0.14145E+08, 0.14641E+08, 0.15151E+08, 0.15675E+08,
     + 0.16214E+08/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_CO    (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 6), QofT( 6,119),Q(NT)
      data xgj/ 1.,2.,1.,6.,2.,12./
c...       CO
c...        --        26
      data (QofT( 1,J),J=1,119)/ 0.21942E+02, 0.30949E+02, 0.39960E+02,
     + 0.48975E+02, 0.57993E+02, 0.67015E+02, 0.76040E+02, 0.85069E+02,
     + 0.94102E+02, 0.10314E+03, 0.11218E+03, 0.12123E+03, 0.13029E+03,
     + 0.13936E+03, 0.14845E+03, 0.15756E+03, 0.16669E+03, 0.17585E+03,
     + 0.18505E+03, 0.19429E+03, 0.20359E+03, 0.21293E+03, 0.22233E+03,
     + 0.23181E+03, 0.24135E+03, 0.25096E+03, 0.26066E+03, 0.27045E+03,
     + 0.28032E+03, 0.29030E+03, 0.30037E+03, 0.31053E+03, 0.32081E+03,
     + 0.33120E+03, 0.34170E+03, 0.35231E+03, 0.36304E+03, 0.37388E+03,
     + 0.38486E+03, 0.39595E+03, 0.40717E+03, 0.41852E+03, 0.42999E+03,
     + 0.44160E+03, 0.45334E+03, 0.46522E+03, 0.47723E+03, 0.48937E+03,
     + 0.50165E+03, 0.51406E+03, 0.52662E+03, 0.53932E+03, 0.55216E+03,
     + 0.56513E+03, 0.57825E+03, 0.59152E+03, 0.60492E+03, 0.61847E+03,
     + 0.63217E+03, 0.64601E+03, 0.65999E+03, 0.67412E+03, 0.68840E+03,
     + 0.70282E+03, 0.71739E+03, 0.73211E+03, 0.74698E+03, 0.76200E+03,
     + 0.77716E+03, 0.79247E+03, 0.80793E+03, 0.82355E+03, 0.83931E+03,
     + 0.85522E+03, 0.87128E+03, 0.88749E+03, 0.90386E+03, 0.92037E+03,
     + 0.93703E+03, 0.95385E+03, 0.97082E+03, 0.98794E+03, 0.10052E+04,
     + 0.10226E+04, 0.10402E+04, 0.10580E+04, 0.10758E+04, 0.10939E+04,
     + 0.11121E+04, 0.11304E+04, 0.11489E+04, 0.11675E+04, 0.11864E+04,
     + 0.12053E+04, 0.12244E+04, 0.12437E+04, 0.12631E+04, 0.12827E+04,
     + 0.13024E+04, 0.13223E+04, 0.13423E+04, 0.13625E+04, 0.13829E+04,
     + 0.14034E+04, 0.14240E+04, 0.14448E+04, 0.14658E+04, 0.14870E+04,
     + 0.15082E+04, 0.15297E+04, 0.15513E+04, 0.15730E+04, 0.15949E+04,
     + 0.16170E+04, 0.16392E+04, 0.16616E+04, 0.16841E+04, 0.17068E+04,
     + 0.17296E+04/
c...        --        36
      data (QofT( 2,J),J=1,119)/ 0.45875E+02, 0.64721E+02, 0.83574E+02,
     + 0.10243E+03, 0.12130E+03, 0.14018E+03, 0.15906E+03, 0.17795E+03,
     + 0.19685E+03, 0.21576E+03, 0.23468E+03, 0.25362E+03, 0.27257E+03,
     + 0.29156E+03, 0.31059E+03, 0.32966E+03, 0.34879E+03, 0.36799E+03,
     + 0.38727E+03, 0.40665E+03, 0.42614E+03, 0.44575E+03, 0.46549E+03,
     + 0.48539E+03, 0.50544E+03, 0.52566E+03, 0.54606E+03, 0.56665E+03,
     + 0.58744E+03, 0.60843E+03, 0.62965E+03, 0.65108E+03, 0.67275E+03,
     + 0.69466E+03, 0.71681E+03, 0.73921E+03, 0.76187E+03, 0.78478E+03,
     + 0.80796E+03, 0.83141E+03, 0.85512E+03, 0.87912E+03, 0.90339E+03,
     + 0.92795E+03, 0.95279E+03, 0.97792E+03, 0.10033E+04, 0.10291E+04,
     + 0.10551E+04, 0.10814E+04, 0.11080E+04, 0.11349E+04, 0.11621E+04,
     + 0.11896E+04, 0.12174E+04, 0.12455E+04, 0.12739E+04, 0.13027E+04,
     + 0.13317E+04, 0.13611E+04, 0.13908E+04, 0.14208E+04, 0.14510E+04,
     + 0.14817E+04, 0.15126E+04, 0.15438E+04, 0.15754E+04, 0.16073E+04,
     + 0.16395E+04, 0.16720E+04, 0.17049E+04, 0.17380E+04, 0.17715E+04,
     + 0.18053E+04, 0.18394E+04, 0.18739E+04, 0.19086E+04, 0.19437E+04,
     + 0.19792E+04, 0.20149E+04, 0.20510E+04, 0.20874E+04, 0.21241E+04,
     + 0.21611E+04, 0.21985E+04, 0.22362E+04, 0.22742E+04, 0.23126E+04,
     + 0.23513E+04, 0.23903E+04, 0.24296E+04, 0.24693E+04, 0.25093E+04,
     + 0.25496E+04, 0.25903E+04, 0.26312E+04, 0.26725E+04, 0.27142E+04,
     + 0.27562E+04, 0.27985E+04, 0.28411E+04, 0.28841E+04, 0.29274E+04,
     + 0.29710E+04, 0.30150E+04, 0.30593E+04, 0.31039E+04, 0.31489E+04,
     + 0.31942E+04, 0.32398E+04, 0.32858E+04, 0.33321E+04, 0.33787E+04,
     + 0.34257E+04, 0.34730E+04, 0.35207E+04, 0.35686E+04, 0.36170E+04,
     + 0.36656E+04/
c...        --        28
      data (QofT( 3,J),J=1,119)/ 0.23024E+02, 0.32483E+02, 0.41946E+02,
     + 0.51412E+02, 0.60882E+02, 0.70356E+02, 0.79834E+02, 0.89315E+02,
     + 0.98801E+02, 0.10829E+03, 0.11779E+03, 0.12729E+03, 0.13681E+03,
     + 0.14634E+03, 0.15589E+03, 0.16546E+03, 0.17506E+03, 0.18470E+03,
     + 0.19438E+03, 0.20411E+03, 0.21389E+03, 0.22374E+03, 0.23365E+03,
     + 0.24364E+03, 0.25371E+03, 0.26386E+03, 0.27411E+03, 0.28444E+03,
     + 0.29489E+03, 0.30543E+03, 0.31608E+03, 0.32685E+03, 0.33773E+03,
     + 0.34873E+03, 0.35986E+03, 0.37111E+03, 0.38249E+03, 0.39400E+03,
     + 0.40565E+03, 0.41742E+03, 0.42934E+03, 0.44139E+03, 0.45359E+03,
     + 0.46592E+03, 0.47841E+03, 0.49103E+03, 0.50380E+03, 0.51672E+03,
     + 0.52979E+03, 0.54300E+03, 0.55637E+03, 0.56989E+03, 0.58356E+03,
     + 0.59738E+03, 0.61136E+03, 0.62549E+03, 0.63977E+03, 0.65421E+03,
     + 0.66881E+03, 0.68357E+03, 0.69847E+03, 0.71354E+03, 0.72877E+03,
     + 0.74415E+03, 0.75969E+03, 0.77540E+03, 0.79126E+03, 0.80728E+03,
     + 0.82346E+03, 0.83981E+03, 0.85631E+03, 0.87297E+03, 0.88980E+03,
     + 0.90679E+03, 0.92394E+03, 0.94125E+03, 0.95873E+03, 0.97636E+03,
     + 0.99417E+03, 0.10121E+04, 0.10303E+04, 0.10485E+04, 0.10670E+04,
     + 0.10856E+04, 0.11044E+04, 0.11234E+04, 0.11425E+04, 0.11617E+04,
     + 0.11812E+04, 0.12008E+04, 0.12206E+04, 0.12405E+04, 0.12606E+04,
     + 0.12809E+04, 0.13013E+04, 0.13219E+04, 0.13427E+04, 0.13636E+04,
     + 0.13847E+04, 0.14060E+04, 0.14274E+04, 0.14490E+04, 0.14708E+04,
     + 0.14927E+04, 0.15148E+04, 0.15371E+04, 0.15595E+04, 0.15821E+04,
     + 0.16049E+04, 0.16278E+04, 0.16509E+04, 0.16742E+04, 0.16976E+04,
     + 0.17212E+04, 0.17450E+04, 0.17690E+04, 0.17931E+04, 0.18174E+04,
     + 0.18418E+04/
c...        --        27
      data (QofT( 4,J),J=1,119)/ 0.13502E+03, 0.19046E+03, 0.24593E+03,
     + 0.30143E+03, 0.35694E+03, 0.41248E+03, 0.46804E+03, 0.52362E+03,
     + 0.57922E+03, 0.63485E+03, 0.69052E+03, 0.74623E+03, 0.80201E+03,
     + 0.85786E+03, 0.91382E+03, 0.96991E+03, 0.10262E+04, 0.10826E+04,
     + 0.11393E+04, 0.11963E+04, 0.12536E+04, 0.13112E+04, 0.13692E+04,
     + 0.14276E+04, 0.14865E+04, 0.15459E+04, 0.16057E+04, 0.16662E+04,
     + 0.17272E+04, 0.17888E+04, 0.18510E+04, 0.19139E+04, 0.19774E+04,
     + 0.20416E+04, 0.21066E+04, 0.21722E+04, 0.22386E+04, 0.23058E+04,
     + 0.23737E+04, 0.24424E+04, 0.25118E+04, 0.25821E+04, 0.26532E+04,
     + 0.27251E+04, 0.27978E+04, 0.28714E+04, 0.29458E+04, 0.30211E+04,
     + 0.30972E+04, 0.31742E+04, 0.32520E+04, 0.33307E+04, 0.34104E+04,
     + 0.34908E+04, 0.35722E+04, 0.36545E+04, 0.37376E+04, 0.38217E+04,
     + 0.39066E+04, 0.39925E+04, 0.40793E+04, 0.41670E+04, 0.42556E+04,
     + 0.43451E+04, 0.44355E+04, 0.45269E+04, 0.46191E+04, 0.47124E+04,
     + 0.48065E+04, 0.49016E+04, 0.49976E+04, 0.50945E+04, 0.51923E+04,
     + 0.52912E+04, 0.53909E+04, 0.54916E+04, 0.55932E+04, 0.56957E+04,
     + 0.57993E+04, 0.59037E+04, 0.60091E+04, 0.61155E+04, 0.62228E+04,
     + 0.63310E+04, 0.64402E+04, 0.65504E+04, 0.66615E+04, 0.67735E+04,
     + 0.68866E+04, 0.70005E+04, 0.71154E+04, 0.72313E+04, 0.73481E+04,
     + 0.74660E+04, 0.75847E+04, 0.77045E+04, 0.78251E+04, 0.79468E+04,
     + 0.80694E+04, 0.81930E+04, 0.83175E+04, 0.84431E+04, 0.85695E+04,
     + 0.86970E+04, 0.88254E+04, 0.89548E+04, 0.90852E+04, 0.92164E+04,
     + 0.93487E+04, 0.94820E+04, 0.96163E+04, 0.97515E+04, 0.98877E+04,
     + 0.10025E+05, 0.10163E+05, 0.10302E+05, 0.10442E+05, 0.10583E+05,
     + 0.10725E+05/
c...        --        38
      data (QofT( 5,J),J=1,119)/ 0.48251E+02, 0.68086E+02, 0.87930E+02,
     + 0.10778E+03, 0.12764E+03, 0.14751E+03, 0.16738E+03, 0.18727E+03,
     + 0.20716E+03, 0.22706E+03, 0.24698E+03, 0.26692E+03, 0.28688E+03,
     + 0.30687E+03, 0.32691E+03, 0.34701E+03, 0.36717E+03, 0.38742E+03,
     + 0.40776E+03, 0.42821E+03, 0.44880E+03, 0.46951E+03, 0.49039E+03,
     + 0.51142E+03, 0.53264E+03, 0.55404E+03, 0.57565E+03, 0.59747E+03,
     + 0.61951E+03, 0.64179E+03, 0.66430E+03, 0.68707E+03, 0.71008E+03,
     + 0.73336E+03, 0.75691E+03, 0.78073E+03, 0.80483E+03, 0.82922E+03,
     + 0.85390E+03, 0.87887E+03, 0.90413E+03, 0.92970E+03, 0.95558E+03,
     + 0.98176E+03, 0.10082E+04, 0.10351E+04, 0.10622E+04, 0.10896E+04,
     + 0.11174E+04, 0.11455E+04, 0.11739E+04, 0.12026E+04, 0.12317E+04,
     + 0.12611E+04, 0.12908E+04, 0.13209E+04, 0.13512E+04, 0.13820E+04,
     + 0.14130E+04, 0.14444E+04, 0.14762E+04, 0.15082E+04, 0.15407E+04,
     + 0.15734E+04, 0.16065E+04, 0.16400E+04, 0.16737E+04, 0.17079E+04,
     + 0.17424E+04, 0.17772E+04, 0.18123E+04, 0.18479E+04, 0.18837E+04,
     + 0.19199E+04, 0.19565E+04, 0.19934E+04, 0.20306E+04, 0.20683E+04,
     + 0.21062E+04, 0.21445E+04, 0.21832E+04, 0.22222E+04, 0.22615E+04,
     + 0.23013E+04, 0.23413E+04, 0.23817E+04, 0.24225E+04, 0.24636E+04,
     + 0.25051E+04, 0.25470E+04, 0.25892E+04, 0.26317E+04, 0.26746E+04,
     + 0.27179E+04, 0.27615E+04, 0.28054E+04, 0.28498E+04, 0.28945E+04,
     + 0.29395E+04, 0.29849E+04, 0.30306E+04, 0.30768E+04, 0.31232E+04,
     + 0.31701E+04, 0.32173E+04, 0.32648E+04, 0.33127E+04, 0.33610E+04,
     + 0.34096E+04, 0.34586E+04, 0.35080E+04, 0.35577E+04, 0.36077E+04,
     + 0.36582E+04, 0.37090E+04, 0.37601E+04, 0.38116E+04, 0.38635E+04,
     + 0.39158E+04/
c...        --        37
      data (QofT( 6,J),J=1,119)/ 0.28263E+03, 0.39877E+03, 0.51497E+03,
     + 0.63121E+03, 0.74749E+03, 0.86382E+03, 0.98020E+03, 0.10966E+04,
     + 0.12131E+04, 0.13296E+04, 0.14462E+04, 0.15629E+04, 0.16797E+04,
     + 0.17966E+04, 0.19138E+04, 0.20313E+04, 0.21490E+04, 0.22672E+04,
     + 0.23858E+04, 0.25049E+04, 0.26247E+04, 0.27452E+04, 0.28665E+04,
     + 0.29886E+04, 0.31116E+04, 0.32356E+04, 0.33606E+04, 0.34868E+04,
     + 0.36141E+04, 0.37427E+04, 0.38725E+04, 0.40036E+04, 0.41361E+04,
     + 0.42700E+04, 0.44054E+04, 0.45422E+04, 0.46805E+04, 0.48204E+04,
     + 0.49618E+04, 0.51049E+04, 0.52495E+04, 0.53958E+04, 0.55438E+04,
     + 0.56934E+04, 0.58448E+04, 0.59979E+04, 0.61527E+04, 0.63092E+04,
     + 0.64675E+04, 0.66276E+04, 0.67895E+04, 0.69531E+04, 0.71187E+04,
     + 0.72859E+04, 0.74551E+04, 0.76261E+04, 0.77989E+04, 0.79736E+04,
     + 0.81501E+04, 0.83285E+04, 0.85088E+04, 0.86910E+04, 0.88751E+04,
     + 0.90610E+04, 0.92488E+04, 0.94385E+04, 0.96302E+04, 0.98238E+04,
     + 0.10019E+05, 0.10217E+05, 0.10416E+05, 0.10617E+05, 0.10820E+05,
     + 0.11025E+05, 0.11233E+05, 0.11441E+05, 0.11652E+05, 0.11865E+05,
     + 0.12080E+05, 0.12297E+05, 0.12516E+05, 0.12736E+05, 0.12959E+05,
     + 0.13184E+05, 0.13410E+05, 0.13639E+05, 0.13869E+05, 0.14102E+05,
     + 0.14336E+05, 0.14573E+05, 0.14811E+05, 0.15051E+05, 0.15294E+05,
     + 0.15538E+05, 0.15784E+05, 0.16033E+05, 0.16283E+05, 0.16535E+05,
     + 0.16790E+05, 0.17046E+05, 0.17304E+05, 0.17564E+05, 0.17827E+05,
     + 0.18091E+05, 0.18357E+05, 0.18626E+05, 0.18896E+05, 0.19168E+05,
     + 0.19443E+05, 0.19719E+05, 0.19997E+05, 0.20277E+05, 0.20560E+05,
     + 0.20844E+05, 0.21130E+05, 0.21419E+05, 0.21709E+05, 0.22002E+05,
     + 0.22296E+05/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_CH4   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 3), QofT( 3,119),Q(NT)
      data xgj/ 1.,2.,3./
c...      CH4
c...        --       211
      data (QofT( 1,J),J=1,119)/ 0.54791E+02, 0.91526E+02, 0.13410E+03,
     + 0.18179E+03, 0.23412E+03, 0.29072E+03, 0.35137E+03, 0.41596E+03,
     + 0.48450E+03, 0.55713E+03, 0.63413E+03, 0.71585E+03, 0.80281E+03,
     + 0.89552E+03, 0.99464E+03, 0.11009E+04, 0.12150E+04, 0.13377E+04,
     + 0.14700E+04, 0.16129E+04, 0.17672E+04, 0.19341E+04, 0.21149E+04,
     + 0.23106E+04, 0.25227E+04, 0.27525E+04, 0.30017E+04, 0.32719E+04,
     + 0.35649E+04, 0.38826E+04, 0.42271E+04, 0.46005E+04, 0.50054E+04,
     + 0.54442E+04, 0.59197E+04, 0.64348E+04, 0.69927E+04, 0.75969E+04,
     + 0.82510E+04, 0.89588E+04, 0.97245E+04, 0.10553E+05, 0.11448E+05,
     + 0.12416E+05, 0.13462E+05, 0.14591E+05, 0.15810E+05, 0.17127E+05,
     + 0.18546E+05, 0.20078E+05, 0.21729E+05, 0.23508E+05, 0.25425E+05,
     + 0.27489E+05, 0.29710E+05, 0.32101E+05, 0.34673E+05, 0.37438E+05,
     + 0.40410E+05, 0.43603E+05, 0.47032E+05, 0.50713E+05, 0.54663E+05,
     + 0.58901E+05, 0.63446E+05, 0.68317E+05, 0.73536E+05, 0.79127E+05,
     + 0.85113E+05, 0.91519E+05, 0.98374E+05, 0.10570E+06, 0.11354E+06,
     + 0.12192E+06, 0.13086E+06, 0.14042E+06, 0.15062E+06, 0.16151E+06,
     + 0.17312E+06, 0.18550E+06, 0.19871E+06, 0.21277E+06, 0.22776E+06,
     + 0.24372E+06, 0.26071E+06, 0.27878E+06, 0.29802E+06, 0.31847E+06,
     + 0.34021E+06, 0.36332E+06, 0.38787E+06, 0.41394E+06, 0.44162E+06,
     + 0.47100E+06, 0.50218E+06, 0.53524E+06, 0.57030E+06, 0.60747E+06,
     + 0.64685E+06, 0.68857E+06, 0.73276E+06, 0.77954E+06, 0.82906E+06,
     + 0.88145E+06, 0.93687E+06, 0.99548E+06, 0.10574E+07, 0.11229E+07,
     + 0.11921E+07, 0.12652E+07, 0.13424E+07, 0.14238E+07, 0.15098E+07,
     + 0.16005E+07, 0.16962E+07, 0.17972E+07, 0.19036E+07, 0.20157E+07,
     + 0.21339E+07/
c...        --       311
      data (QofT( 2,J),J=1,119)/ 0.10958E+03, 0.18304E+03, 0.26818E+03,
     + 0.36356E+03, 0.46820E+03, 0.58141E+03, 0.70270E+03, 0.83186E+03,
     + 0.96893E+03, 0.11142E+04, 0.12682E+04, 0.14316E+04, 0.16055E+04,
     + 0.17909E+04, 0.19891E+04, 0.22016E+04, 0.24297E+04, 0.26752E+04,
     + 0.29399E+04, 0.32255E+04, 0.35342E+04, 0.38680E+04, 0.42294E+04,
     + 0.46208E+04, 0.50449E+04, 0.55046E+04, 0.60030E+04, 0.65434E+04,
     + 0.71293E+04, 0.77646E+04, 0.84535E+04, 0.92004E+04, 0.10010E+05,
     + 0.10888E+05, 0.11838E+05, 0.12869E+05, 0.13984E+05, 0.15193E+05,
     + 0.16501E+05, 0.17916E+05, 0.19448E+05, 0.21104E+05, 0.22895E+05,
     + 0.24830E+05, 0.26921E+05, 0.29180E+05, 0.31618E+05, 0.34250E+05,
     + 0.37090E+05, 0.40152E+05, 0.43454E+05, 0.47012E+05, 0.50845E+05,
     + 0.54973E+05, 0.59416E+05, 0.64197E+05, 0.69340E+05, 0.74870E+05,
     + 0.80813E+05, 0.87198E+05, 0.94055E+05, 0.10142E+06, 0.10932E+06,
     + 0.11779E+06, 0.12688E+06, 0.13662E+06, 0.14706E+06, 0.15824E+06,
     + 0.17021E+06, 0.18302E+06, 0.19673E+06, 0.21139E+06, 0.22706E+06,
     + 0.24381E+06, 0.26171E+06, 0.28082E+06, 0.30122E+06, 0.32299E+06,
     + 0.34621E+06, 0.37097E+06, 0.39737E+06, 0.42551E+06, 0.45548E+06,
     + 0.48739E+06, 0.52136E+06, 0.55752E+06, 0.59598E+06, 0.63688E+06,
     + 0.68036E+06, 0.72657E+06, 0.77566E+06, 0.82780E+06, 0.88316E+06,
     + 0.94191E+06, 0.10043E+07, 0.10704E+07, 0.11405E+07, 0.12148E+07,
     + 0.12936E+07, 0.13770E+07, 0.14654E+07, 0.15589E+07, 0.16579E+07,
     + 0.17627E+07, 0.18736E+07, 0.19908E+07, 0.21147E+07, 0.22456E+07,
     + 0.23840E+07, 0.25301E+07, 0.26844E+07, 0.28474E+07, 0.30193E+07,
     + 0.32007E+07, 0.33921E+07, 0.35939E+07, 0.38067E+07, 0.40310E+07,
     + 0.42673E+07/
c...        --       212
      data (QofT( 3,J),J=1,119)/ 0.44079E+03, 0.73786E+03, 0.10822E+04,
     + 0.14679E+04, 0.18912E+04, 0.23493E+04, 0.28400E+04, 0.33627E+04,
     + 0.39174E+04, 0.45053E+04, 0.51285E+04, 0.57900E+04, 0.64938E+04,
     + 0.72443E+04, 0.80465E+04, 0.89064E+04, 0.98299E+04, 0.10824E+05,
     + 0.11895E+05, 0.13051E+05, 0.14300E+05, 0.15652E+05, 0.17115E+05,
     + 0.18699E+05, 0.20416E+05, 0.22277E+05, 0.24294E+05, 0.26481E+05,
     + 0.28853E+05, 0.31425E+05, 0.34214E+05, 0.37237E+05, 0.40515E+05,
     + 0.44067E+05, 0.47916E+05, 0.52087E+05, 0.56604E+05, 0.61495E+05,
     + 0.66790E+05, 0.72521E+05, 0.78720E+05, 0.85426E+05, 0.92675E+05,
     + 0.10051E+06, 0.10898E+06, 0.11812E+06, 0.12799E+06, 0.13865E+06,
     + 0.15014E+06, 0.16254E+06, 0.17591E+06, 0.19031E+06, 0.20583E+06,
     + 0.22254E+06, 0.24053E+06, 0.25989E+06, 0.28071E+06, 0.30310E+06,
     + 0.32716E+06, 0.35301E+06, 0.38077E+06, 0.41058E+06, 0.44257E+06,
     + 0.47688E+06, 0.51367E+06, 0.55311E+06, 0.59537E+06, 0.64064E+06,
     + 0.68910E+06, 0.74098E+06, 0.79647E+06, 0.85583E+06, 0.91928E+06,
     + 0.98710E+06, 0.10595E+07, 0.11369E+07, 0.12195E+07, 0.13076E+07,
     + 0.14017E+07, 0.15019E+07, 0.16088E+07, 0.17227E+07, 0.18441E+07,
     + 0.19733E+07, 0.21108E+07, 0.22572E+07, 0.24129E+07, 0.25785E+07,
     + 0.27545E+07, 0.29416E+07, 0.31404E+07, 0.33515E+07, 0.35756E+07,
     + 0.38135E+07, 0.40659E+07, 0.43336E+07, 0.46174E+07, 0.49183E+07,
     + 0.52372E+07, 0.55750E+07, 0.59327E+07, 0.63114E+07, 0.67123E+07,
     + 0.71365E+07, 0.75852E+07, 0.80596E+07, 0.85612E+07, 0.90914E+07,
     + 0.96515E+07, 0.10243E+08, 0.10868E+08, 0.11527E+08, 0.12224E+08,
     + 0.12958E+08, 0.13733E+08, 0.14550E+08, 0.15411E+08, 0.16319E+08,
     + 0.17275E+08/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_O2    (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 3), QofT( 3,119),Q(NT)
      data xgj/ 1.,1.,6./
c...       O2
c...        --        66
      data (QofT( 1,J),J=1,119)/ 0.44334E+02, 0.62460E+02, 0.80596E+02,
     + 0.98738E+02, 0.11688E+03, 0.13503E+03, 0.15319E+03, 0.17136E+03,
     + 0.18954E+03, 0.20775E+03, 0.22600E+03, 0.24431E+03, 0.26270E+03,
     + 0.28119E+03, 0.29981E+03, 0.31857E+03, 0.33750E+03, 0.35662E+03,
     + 0.37594E+03, 0.39550E+03, 0.41529E+03, 0.43535E+03, 0.45568E+03,
     + 0.47630E+03, 0.49722E+03, 0.51844E+03, 0.53998E+03, 0.56185E+03,
     + 0.58406E+03, 0.60660E+03, 0.62949E+03, 0.65274E+03, 0.67635E+03,
     + 0.70031E+03, 0.72465E+03, 0.74936E+03, 0.77444E+03, 0.79990E+03,
     + 0.82574E+03, 0.85197E+03, 0.87858E+03, 0.90558E+03, 0.93297E+03,
     + 0.96076E+03, 0.98895E+03, 0.10175E+04, 0.10465E+04, 0.10759E+04,
     + 0.11057E+04, 0.11359E+04, 0.11665E+04, 0.11976E+04, 0.12290E+04,
     + 0.12609E+04, 0.12931E+04, 0.13258E+04, 0.13590E+04, 0.13925E+04,
     + 0.14265E+04, 0.14609E+04, 0.14958E+04, 0.15311E+04, 0.15669E+04,
     + 0.16031E+04, 0.16397E+04, 0.16768E+04, 0.17144E+04, 0.17524E+04,
     + 0.17909E+04, 0.18298E+04, 0.18692E+04, 0.19091E+04, 0.19495E+04,
     + 0.19904E+04, 0.20318E+04, 0.20736E+04, 0.21160E+04, 0.21588E+04,
     + 0.22022E+04, 0.22461E+04, 0.22905E+04, 0.23354E+04, 0.23809E+04,
     + 0.24268E+04, 0.24734E+04, 0.25204E+04, 0.25680E+04, 0.26162E+04,
     + 0.26649E+04, 0.27142E+04, 0.27641E+04, 0.28145E+04, 0.28655E+04,
     + 0.29171E+04, 0.29693E+04, 0.30221E+04, 0.30755E+04, 0.31295E+04,
     + 0.31841E+04, 0.32393E+04, 0.32951E+04, 0.33516E+04, 0.34087E+04,
     + 0.34665E+04, 0.35249E+04, 0.35839E+04, 0.36436E+04, 0.37040E+04,
     + 0.37650E+04, 0.38267E+04, 0.38891E+04, 0.39522E+04, 0.40159E+04,
     + 0.40804E+04, 0.41455E+04, 0.42114E+04, 0.42780E+04, 0.43452E+04,
     + 0.44132E+04/
c...        --        68
      data (QofT( 2,J),J=1,119)/ 0.89206E+02, 0.12759E+03, 0.16600E+03,
     + 0.20442E+03, 0.24285E+03, 0.28128E+03, 0.31973E+03, 0.35821E+03,
     + 0.39672E+03, 0.43530E+03, 0.47398E+03, 0.51281E+03, 0.55183E+03,
     + 0.59108E+03, 0.63062E+03, 0.67051E+03, 0.71078E+03, 0.75148E+03,
     + 0.79265E+03, 0.83435E+03, 0.87659E+03, 0.91941E+03, 0.96285E+03,
     + 0.10069E+04, 0.10517E+04, 0.10971E+04, 0.11432E+04, 0.11901E+04,
     + 0.12377E+04, 0.12861E+04, 0.13352E+04, 0.13851E+04, 0.14358E+04,
     + 0.14872E+04, 0.15395E+04, 0.15926E+04, 0.16466E+04, 0.17013E+04,
     + 0.17569E+04, 0.18134E+04, 0.18706E+04, 0.19288E+04, 0.19877E+04,
     + 0.20476E+04, 0.21083E+04, 0.21698E+04, 0.22323E+04, 0.22956E+04,
     + 0.23598E+04, 0.24248E+04, 0.24908E+04, 0.25576E+04, 0.26253E+04,
     + 0.26940E+04, 0.27635E+04, 0.28339E+04, 0.29052E+04, 0.29775E+04,
     + 0.30506E+04, 0.31247E+04, 0.31997E+04, 0.32756E+04, 0.33524E+04,
     + 0.34302E+04, 0.35089E+04, 0.35885E+04, 0.36691E+04, 0.37506E+04,
     + 0.38331E+04, 0.39166E+04, 0.40010E+04, 0.40864E+04, 0.41727E+04,
     + 0.42601E+04, 0.43484E+04, 0.44377E+04, 0.45280E+04, 0.46193E+04,
     + 0.47116E+04, 0.48049E+04, 0.48992E+04, 0.49946E+04, 0.50909E+04,
     + 0.51883E+04, 0.52868E+04, 0.53863E+04, 0.54868E+04, 0.55884E+04,
     + 0.56911E+04, 0.57949E+04, 0.58997E+04, 0.60056E+04, 0.61126E+04,
     + 0.62207E+04, 0.63298E+04, 0.64401E+04, 0.65516E+04, 0.66641E+04,
     + 0.67778E+04, 0.68926E+04, 0.70085E+04, 0.71256E+04, 0.72439E+04,
     + 0.73633E+04, 0.74839E+04, 0.76056E+04, 0.77286E+04, 0.78527E+04,
     + 0.79781E+04, 0.81046E+04, 0.82324E+04, 0.83613E+04, 0.84915E+04,
     + 0.86229E+04, 0.87556E+04, 0.88895E+04, 0.90247E+04, 0.91611E+04,
     + 0.92988E+04/
c...        --        67
      data (QofT( 3,J),J=1,119)/ 0.52071E+03, 0.74484E+03, 0.96908E+03,
     + 0.11934E+04, 0.14177E+04, 0.16422E+04, 0.18667E+04, 0.20913E+04,
     + 0.23161E+04, 0.25413E+04, 0.27671E+04, 0.29936E+04, 0.32212E+04,
     + 0.34501E+04, 0.36806E+04, 0.39130E+04, 0.41476E+04, 0.43846E+04,
     + 0.46242E+04, 0.48668E+04, 0.51125E+04, 0.53615E+04, 0.56140E+04,
     + 0.58701E+04, 0.61300E+04, 0.63938E+04, 0.66617E+04, 0.69337E+04,
     + 0.72099E+04, 0.74904E+04, 0.77754E+04, 0.80647E+04, 0.83586E+04,
     + 0.86571E+04, 0.89602E+04, 0.92680E+04, 0.95805E+04, 0.98977E+04,
     + 0.10220E+05, 0.10547E+05, 0.10878E+05, 0.11215E+05, 0.11556E+05,
     + 0.11903E+05, 0.12254E+05, 0.12611E+05, 0.12972E+05, 0.13338E+05,
     + 0.13710E+05, 0.14086E+05, 0.14468E+05, 0.14855E+05, 0.15247E+05,
     + 0.15644E+05, 0.16046E+05, 0.16453E+05, 0.16866E+05, 0.17283E+05,
     + 0.17706E+05, 0.18135E+05, 0.18568E+05, 0.19007E+05, 0.19452E+05,
     + 0.19901E+05, 0.20356E+05, 0.20817E+05, 0.21283E+05, 0.21754E+05,
     + 0.22231E+05, 0.22713E+05, 0.23201E+05, 0.23695E+05, 0.24194E+05,
     + 0.24699E+05, 0.25209E+05, 0.25725E+05, 0.26247E+05, 0.26775E+05,
     + 0.27308E+05, 0.27847E+05, 0.28393E+05, 0.28944E+05, 0.29500E+05,
     + 0.30063E+05, 0.30632E+05, 0.31207E+05, 0.31788E+05, 0.32375E+05,
     + 0.32968E+05, 0.33568E+05, 0.34173E+05, 0.34785E+05, 0.35403E+05,
     + 0.36028E+05, 0.36659E+05, 0.37296E+05, 0.37939E+05, 0.38590E+05,
     + 0.39246E+05, 0.39909E+05, 0.40579E+05, 0.41256E+05, 0.41939E+05,
     + 0.42629E+05, 0.43325E+05, 0.44029E+05, 0.44739E+05, 0.45456E+05,
     + 0.46180E+05, 0.46911E+05, 0.47649E+05, 0.48394E+05, 0.49146E+05,
     + 0.49905E+05, 0.50671E+05, 0.51445E+05, 0.52226E+05, 0.53014E+05,
     + 0.53809E+05/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_NO    (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 3), QofT( 3,119),Q(NT)
      data xgj/ 3.,2.,3./
c...       NO
c...        --        46
      data (QofT( 1,J),J=1,119)/ 0.15840E+03, 0.23971E+03, 0.33080E+03,
     + 0.42907E+03, 0.53251E+03, 0.63972E+03, 0.74975E+03, 0.86195E+03,
     + 0.97582E+03, 0.10911E+04, 0.12074E+04, 0.13248E+04, 0.14430E+04,
     + 0.15621E+04, 0.16820E+04, 0.18027E+04, 0.19243E+04, 0.20468E+04,
     + 0.21703E+04, 0.22948E+04, 0.24204E+04, 0.25472E+04, 0.26753E+04,
     + 0.28046E+04, 0.29354E+04, 0.30676E+04, 0.32013E+04, 0.33365E+04,
     + 0.34734E+04, 0.36120E+04, 0.37522E+04, 0.38942E+04, 0.40379E+04,
     + 0.41835E+04, 0.43310E+04, 0.44803E+04, 0.46316E+04, 0.47849E+04,
     + 0.49400E+04, 0.50972E+04, 0.52564E+04, 0.54176E+04, 0.55809E+04,
     + 0.57462E+04, 0.59137E+04, 0.60832E+04, 0.62548E+04, 0.64286E+04,
     + 0.66045E+04, 0.67825E+04, 0.69628E+04, 0.71451E+04, 0.73297E+04,
     + 0.75164E+04, 0.77053E+04, 0.78964E+04, 0.80897E+04, 0.82853E+04,
     + 0.84830E+04, 0.86830E+04, 0.88852E+04, 0.90896E+04, 0.92963E+04,
     + 0.95052E+04, 0.97164E+04, 0.99297E+04, 0.10145E+05, 0.10363E+05,
     + 0.10583E+05, 0.10806E+05, 0.11031E+05, 0.11258E+05, 0.11487E+05,
     + 0.11718E+05, 0.11952E+05, 0.12188E+05, 0.12426E+05, 0.12667E+05,
     + 0.12910E+05, 0.13155E+05, 0.13403E+05, 0.13652E+05, 0.13905E+05,
     + 0.14159E+05, 0.14416E+05, 0.14675E+05, 0.14936E+05, 0.15199E+05,
     + 0.15465E+05, 0.15733E+05, 0.16004E+05, 0.16277E+05, 0.16552E+05,
     + 0.16829E+05, 0.17109E+05, 0.17391E+05, 0.17675E+05, 0.17962E+05,
     + 0.18251E+05, 0.18542E+05, 0.18836E+05, 0.19131E+05, 0.19430E+05,
     + 0.19730E+05, 0.20033E+05, 0.20338E+05, 0.20646E+05, 0.20955E+05,
     + 0.21268E+05, 0.21582E+05, 0.21899E+05, 0.22218E+05, 0.22539E+05,
     + 0.22863E+05, 0.23189E+05, 0.23518E+05, 0.23848E+05, 0.24181E+05,
     + 0.24517E+05/
c...        --        56
      data (QofT( 2,J),J=1,119)/ 0.10942E+03, 0.16560E+03, 0.22856E+03,
     + 0.29647E+03, 0.36795E+03, 0.44204E+03, 0.51808E+03, 0.59561E+03,
     + 0.67432E+03, 0.75396E+03, 0.83439E+03, 0.91551E+03, 0.99725E+03,
     + 0.10796E+04, 0.11625E+04, 0.12460E+04, 0.13302E+04, 0.14150E+04,
     + 0.15005E+04, 0.15868E+04, 0.16739E+04, 0.17618E+04, 0.18506E+04,
     + 0.19404E+04, 0.20311E+04, 0.21229E+04, 0.22158E+04, 0.23098E+04,
     + 0.24050E+04, 0.25013E+04, 0.25989E+04, 0.26976E+04, 0.27977E+04,
     + 0.28991E+04, 0.30018E+04, 0.31058E+04, 0.32112E+04, 0.33180E+04,
     + 0.34262E+04, 0.35358E+04, 0.36468E+04, 0.37593E+04, 0.38732E+04,
     + 0.39885E+04, 0.41054E+04, 0.42237E+04, 0.43436E+04, 0.44649E+04,
     + 0.45877E+04, 0.47121E+04, 0.48379E+04, 0.49654E+04, 0.50943E+04,
     + 0.52248E+04, 0.53568E+04, 0.54904E+04, 0.56255E+04, 0.57622E+04,
     + 0.59004E+04, 0.60403E+04, 0.61816E+04, 0.63246E+04, 0.64692E+04,
     + 0.66152E+04, 0.67630E+04, 0.69123E+04, 0.70631E+04, 0.72156E+04,
     + 0.73696E+04, 0.75253E+04, 0.76825E+04, 0.78414E+04, 0.80018E+04,
     + 0.81638E+04, 0.83275E+04, 0.84927E+04, 0.86596E+04, 0.88280E+04,
     + 0.89981E+04, 0.91698E+04, 0.93430E+04, 0.95180E+04, 0.96945E+04,
     + 0.98726E+04, 0.10052E+05, 0.10234E+05, 0.10417E+05, 0.10601E+05,
     + 0.10788E+05, 0.10975E+05, 0.11165E+05, 0.11356E+05, 0.11549E+05,
     + 0.11743E+05, 0.11939E+05, 0.12137E+05, 0.12336E+05, 0.12537E+05,
     + 0.12739E+05, 0.12943E+05, 0.13149E+05, 0.13356E+05, 0.13565E+05,
     + 0.13776E+05, 0.13988E+05, 0.14202E+05, 0.14418E+05, 0.14635E+05,
     + 0.14853E+05, 0.15074E+05, 0.15296E+05, 0.15520E+05, 0.15745E+05,
     + 0.15972E+05, 0.16200E+05, 0.16431E+05, 0.16663E+05, 0.16896E+05,
     + 0.17131E+05/
c...        --        48
      data (QofT( 3,J),J=1,119)/ 0.16695E+03, 0.25269E+03, 0.34876E+03,
     + 0.45239E+03, 0.56148E+03, 0.67455E+03, 0.79059E+03, 0.90891E+03,
     + 0.10290E+04, 0.11506E+04, 0.12733E+04, 0.13971E+04, 0.15219E+04,
     + 0.16476E+04, 0.17742E+04, 0.19017E+04, 0.20302E+04, 0.21598E+04,
     + 0.22904E+04, 0.24223E+04, 0.25553E+04, 0.26897E+04, 0.28255E+04,
     + 0.29628E+04, 0.31016E+04, 0.32420E+04, 0.33842E+04, 0.35280E+04,
     + 0.36736E+04, 0.38211E+04, 0.39704E+04, 0.41217E+04, 0.42750E+04,
     + 0.44302E+04, 0.45876E+04, 0.47469E+04, 0.49084E+04, 0.50720E+04,
     + 0.52378E+04, 0.54058E+04, 0.55759E+04, 0.57483E+04, 0.59230E+04,
     + 0.60999E+04, 0.62791E+04, 0.64605E+04, 0.66443E+04, 0.68304E+04,
     + 0.70187E+04, 0.72095E+04, 0.74026E+04, 0.75980E+04, 0.77958E+04,
     + 0.79960E+04, 0.81986E+04, 0.84036E+04, 0.86109E+04, 0.88207E+04,
     + 0.90328E+04, 0.92474E+04, 0.94644E+04, 0.96839E+04, 0.99057E+04,
     + 0.10130E+05, 0.10357E+05, 0.10586E+05, 0.10817E+05, 0.11052E+05,
     + 0.11288E+05, 0.11527E+05, 0.11768E+05, 0.12012E+05, 0.12259E+05,
     + 0.12507E+05, 0.12759E+05, 0.13012E+05, 0.13269E+05, 0.13527E+05,
     + 0.13788E+05, 0.14052E+05, 0.14318E+05, 0.14587E+05, 0.14858E+05,
     + 0.15131E+05, 0.15408E+05, 0.15686E+05, 0.15967E+05, 0.16251E+05,
     + 0.16537E+05, 0.16825E+05, 0.17116E+05, 0.17410E+05, 0.17706E+05,
     + 0.18004E+05, 0.18305E+05, 0.18609E+05, 0.18915E+05, 0.19224E+05,
     + 0.19535E+05, 0.19848E+05, 0.20164E+05, 0.20483E+05, 0.20804E+05,
     + 0.21127E+05, 0.21453E+05, 0.21782E+05, 0.22113E+05, 0.22447E+05,
     + 0.22783E+05, 0.23122E+05, 0.23463E+05, 0.23807E+05, 0.24153E+05,
     + 0.24502E+05, 0.24853E+05, 0.25207E+05, 0.25563E+05, 0.25922E+05,
     + 0.26283E+05/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_SO2   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 2), QofT( 2,119),Q(NT)
      data xgj/ 1.,1./
c...      SO2
c...        --       626
      data (QofT( 1,J),J=1,119)/ 0.52899E+03, 0.89171E+03, 0.13139E+04,
     + 0.17915E+04, 0.23246E+04, 0.29155E+04, 0.35675E+04, 0.42848E+04,
     + 0.50723E+04, 0.59352E+04, 0.68794E+04, 0.79109E+04, 0.90366E+04,
     + 0.10264E+05, 0.11599E+05, 0.13052E+05, 0.14629E+05, 0.16340E+05,
     + 0.18193E+05, 0.20199E+05, 0.22366E+05, 0.24704E+05, 0.27225E+05,
     + 0.29938E+05, 0.32855E+05, 0.35987E+05, 0.39346E+05, 0.42944E+05,
     + 0.46794E+05, 0.50909E+05, 0.55302E+05, 0.59986E+05, 0.64977E+05,
     + 0.70288E+05, 0.75934E+05, 0.81931E+05, 0.88294E+05, 0.95040E+05,
     + 0.10219E+06, 0.10975E+06, 0.11774E+06, 0.12619E+06, 0.13511E+06,
     + 0.14452E+06, 0.15443E+06, 0.16487E+06, 0.17586E+06, 0.18742E+06,
     + 0.19957E+06, 0.21234E+06, 0.22573E+06, 0.23978E+06, 0.25451E+06,
     + 0.26995E+06, 0.28611E+06, 0.30302E+06, 0.32071E+06, 0.33920E+06,
     + 0.35852E+06, 0.37869E+06, 0.39974E+06, 0.42171E+06, 0.44461E+06,
     + 0.46848E+06, 0.49334E+06, 0.51922E+06, 0.54617E+06, 0.57419E+06,
     + 0.60334E+06, 0.63363E+06, 0.66511E+06, 0.69780E+06, 0.73174E+06,
     + 0.76696E+06, 0.80349E+06, 0.84138E+06, 0.88066E+06, 0.92136E+06,
     + 0.96352E+06, 0.10072E+07, 0.10524E+07, 0.10992E+07, 0.11475E+07,
     + 0.11976E+07, 0.12493E+07, 0.13028E+07, 0.13580E+07, 0.14151E+07,
     + 0.14741E+07, 0.15349E+07, 0.15977E+07, 0.16625E+07, 0.17293E+07,
     + 0.17982E+07, 0.18693E+07, 0.19425E+07, 0.20180E+07, 0.20958E+07,
     + 0.21758E+07, 0.22583E+07, 0.23432E+07, 0.24305E+07, 0.25204E+07,
     + 0.26129E+07, 0.27080E+07, 0.28058E+07, 0.29064E+07, 0.30097E+07,
     + 0.31159E+07, 0.32250E+07, 0.33371E+07, 0.34522E+07, 0.35705E+07,
     + 0.36918E+07, 0.38164E+07, 0.39442E+07, 0.40754E+07, 0.42099E+07,
     + 0.43479E+07/
c...        --       646
      data (QofT( 2,J),J=1,119)/ 0.53140E+03, 0.89578E+03, 0.13199E+04,
     + 0.17997E+04, 0.23353E+04, 0.29288E+04, 0.35837E+04, 0.43043E+04,
     + 0.50953E+04, 0.59621E+04, 0.69104E+04, 0.79465E+04, 0.90772E+04,
     + 0.10310E+05, 0.11651E+05, 0.13110E+05, 0.14694E+05, 0.16413E+05,
     + 0.18274E+05, 0.20289E+05, 0.22465E+05, 0.24814E+05, 0.27345E+05,
     + 0.30070E+05, 0.33000E+05, 0.36145E+05, 0.39519E+05, 0.43133E+05,
     + 0.46999E+05, 0.51132E+05, 0.55544E+05, 0.60248E+05, 0.65260E+05,
     + 0.70594E+05, 0.76264E+05, 0.82287E+05, 0.88678E+05, 0.95453E+05,
     + 0.10263E+06, 0.11022E+06, 0.11825E+06, 0.12674E+06, 0.13569E+06,
     + 0.14514E+06, 0.15510E+06, 0.16558E+06, 0.17662E+06, 0.18823E+06,
     + 0.20043E+06, 0.21325E+06, 0.22670E+06, 0.24081E+06, 0.25561E+06,
     + 0.27111E+06, 0.28733E+06, 0.30432E+06, 0.32208E+06, 0.34065E+06,
     + 0.36005E+06, 0.38031E+06, 0.40145E+06, 0.42351E+06, 0.44651E+06,
     + 0.47047E+06, 0.49544E+06, 0.52144E+06, 0.54849E+06, 0.57664E+06,
     + 0.60591E+06, 0.63633E+06, 0.66794E+06, 0.70077E+06, 0.73485E+06,
     + 0.77022E+06, 0.80691E+06, 0.84496E+06, 0.88440E+06, 0.92527E+06,
     + 0.96761E+06, 0.10115E+07, 0.10568E+07, 0.11038E+07, 0.11524E+07,
     + 0.12027E+07, 0.12546E+07, 0.13083E+07, 0.13638E+07, 0.14211E+07,
     + 0.14803E+07, 0.15414E+07, 0.16045E+07, 0.16695E+07, 0.17366E+07,
     + 0.18059E+07, 0.18772E+07, 0.19507E+07, 0.20265E+07, 0.21046E+07,
     + 0.21850E+07, 0.22678E+07, 0.23531E+07, 0.24408E+07, 0.25310E+07,
     + 0.26239E+07, 0.27194E+07, 0.28176E+07, 0.29186E+07, 0.30224E+07,
     + 0.31290E+07, 0.32386E+07, 0.33512E+07, 0.34668E+07, 0.35855E+07,
     + 0.37074E+07, 0.38324E+07, 0.39608E+07, 0.40925E+07, 0.42276E+07,
     + 0.43662E+07/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_NO2   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 3./
c...      NO2
c...        --       646
      data (QofT( 1,J),J=1,119)/ 0.12046E+04, 0.20297E+04, 0.29875E+04,
     + 0.40626E+04, 0.52463E+04, 0.65350E+04, 0.79286E+04, 0.94298E+04,
     + 0.11043E+05, 0.12776E+05, 0.14634E+05, 0.16627E+05, 0.18765E+05,
     + 0.21056E+05, 0.23511E+05, 0.26143E+05, 0.28961E+05, 0.31979E+05,
     + 0.35209E+05, 0.38663E+05, 0.42355E+05, 0.46300E+05, 0.50510E+05,
     + 0.55001E+05, 0.59787E+05, 0.64884E+05, 0.70308E+05, 0.76075E+05,
     + 0.82201E+05, 0.88704E+05, 0.95602E+05, 0.10291E+06, 0.11065E+06,
     + 0.11884E+06, 0.12750E+06, 0.13665E+06, 0.14631E+06, 0.15650E+06,
     + 0.16724E+06, 0.17856E+06, 0.19047E+06, 0.20301E+06, 0.21618E+06,
     + 0.23002E+06, 0.24456E+06, 0.25981E+06, 0.27580E+06, 0.29256E+06,
     + 0.31012E+06, 0.32850E+06, 0.34773E+06, 0.36784E+06, 0.38886E+06,
     + 0.41082E+06, 0.43374E+06, 0.45766E+06, 0.48262E+06, 0.50863E+06,
     + 0.53574E+06, 0.56398E+06, 0.59339E+06, 0.62398E+06, 0.65581E+06,
     + 0.68891E+06, 0.72331E+06, 0.75905E+06, 0.79617E+06, 0.83470E+06,
     + 0.87469E+06, 0.91617E+06, 0.95919E+06, 0.10038E+07, 0.10500E+07,
     + 0.10979E+07, 0.11474E+07, 0.11988E+07, 0.12519E+07, 0.13068E+07,
     + 0.13636E+07, 0.14224E+07, 0.14831E+07, 0.15459E+07, 0.16107E+07,
     + 0.16776E+07, 0.17467E+07, 0.18180E+07, 0.18916E+07, 0.19675E+07,
     + 0.20458E+07, 0.21265E+07, 0.22097E+07, 0.22954E+07, 0.23837E+07,
     + 0.24747E+07, 0.25684E+07, 0.26648E+07, 0.27641E+07, 0.28662E+07,
     + 0.29713E+07, 0.30794E+07, 0.31905E+07, 0.33048E+07, 0.34223E+07,
     + 0.35430E+07, 0.36670E+07, 0.37944E+07, 0.39253E+07, 0.40597E+07,
     + 0.41976E+07, 0.43393E+07, 0.44846E+07, 0.46337E+07, 0.47867E+07,
     + 0.49437E+07, 0.51046E+07, 0.52696E+07, 0.54388E+07, 0.56122E+07,
     + 0.57900E+07/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_NH3   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 2), QofT( 2,119),Q(NT)
      data xgj/ 3.,2./
c...      NH3
c...        --      4111
      data (QofT( 1,J),J=1,119)/ 0.16013E+03, 0.26692E+03, 0.39067E+03,
     + 0.52933E+03, 0.68153E+03, 0.84641E+03, 0.10234E+04, 0.12125E+04,
     + 0.14136E+04, 0.16272E+04, 0.18537E+04, 0.20937E+04, 0.23481E+04,
     + 0.26177E+04, 0.29035E+04, 0.32065E+04, 0.35279E+04, 0.38688E+04,
     + 0.42304E+04, 0.46141E+04, 0.50212E+04, 0.54531E+04, 0.59114E+04,
     + 0.63976E+04, 0.69133E+04, 0.74602E+04, 0.80401E+04, 0.86549E+04,
     + 0.93066E+04, 0.99971E+04, 0.10729E+05, 0.11504E+05, 0.12324E+05,
     + 0.13193E+05, 0.14112E+05, 0.15085E+05, 0.16114E+05, 0.17201E+05,
     + 0.18352E+05, 0.19567E+05, 0.20851E+05, 0.22208E+05, 0.23640E+05,
     + 0.25152E+05, 0.26747E+05, 0.28430E+05, 0.30205E+05, 0.32077E+05,
     + 0.34050E+05, 0.36128E+05, 0.38317E+05, 0.40623E+05, 0.43050E+05,
     + 0.45605E+05, 0.48292E+05, 0.51119E+05, 0.54091E+05, 0.57215E+05,
     + 0.60498E+05, 0.63947E+05, 0.67569E+05, 0.71372E+05, 0.75364E+05,
     + 0.79552E+05, 0.83946E+05, 0.88553E+05, 0.93384E+05, 0.98447E+05,
     + 0.10375E+06, 0.10931E+06, 0.11513E+06, 0.12122E+06, 0.12760E+06,
     + 0.13427E+06, 0.14125E+06, 0.14855E+06, 0.15619E+06, 0.16417E+06,
     + 0.17250E+06, 0.18121E+06, 0.19031E+06, 0.19981E+06, 0.20973E+06,
     + 0.22008E+06, 0.23088E+06, 0.24215E+06, 0.25390E+06, 0.26615E+06,
     + 0.27892E+06, 0.29223E+06, 0.30610E+06, 0.32055E+06, 0.33559E+06,
     + 0.35125E+06, 0.36756E+06, 0.38453E+06, 0.40219E+06, 0.42056E+06,
     + 0.43967E+06, 0.45953E+06, 0.48019E+06, 0.50165E+06, 0.52396E+06,
     + 0.54714E+06, 0.57122E+06, 0.59622E+06, 0.62218E+06, 0.64913E+06,
     + 0.67710E+06, 0.70613E+06, 0.73624E+06, 0.76748E+06, 0.79988E+06,
     + 0.83347E+06, 0.86829E+06, 0.90439E+06, 0.94180E+06, 0.98056E+06,
     + 0.10207E+07/
c...        --      5111
      data (QofT( 2,J),J=1,119)/ 0.10697E+03, 0.17832E+03, 0.26100E+03,
     + 0.35364E+03, 0.45533E+03, 0.56549E+03, 0.68377E+03, 0.81007E+03,
     + 0.94447E+03, 0.10872E+04, 0.12385E+04, 0.13988E+04, 0.15688E+04,
     + 0.17490E+04, 0.19399E+04, 0.21424E+04, 0.23571E+04, 0.25848E+04,
     + 0.28264E+04, 0.30828E+04, 0.33548E+04, 0.36434E+04, 0.39496E+04,
     + 0.42745E+04, 0.46190E+04, 0.49845E+04, 0.53720E+04, 0.57828E+04,
     + 0.62182E+04, 0.66796E+04, 0.71684E+04, 0.76862E+04, 0.82344E+04,
     + 0.88149E+04, 0.94292E+04, 0.10079E+05, 0.10767E+05, 0.11494E+05,
     + 0.12262E+05, 0.13074E+05, 0.13932E+05, 0.14839E+05, 0.15796E+05,
     + 0.16806E+05, 0.17872E+05, 0.18997E+05, 0.20183E+05, 0.21434E+05,
     + 0.22752E+05, 0.24141E+05, 0.25604E+05, 0.27145E+05, 0.28767E+05,
     + 0.30475E+05, 0.32271E+05, 0.34160E+05, 0.36146E+05, 0.38234E+05,
     + 0.40428E+05, 0.42733E+05, 0.45154E+05, 0.47696E+05, 0.50364E+05,
     + 0.53163E+05, 0.56100E+05, 0.59180E+05, 0.62408E+05, 0.65792E+05,
     + 0.69339E+05, 0.73053E+05, 0.76943E+05, 0.81016E+05, 0.85279E+05,
     + 0.89740E+05, 0.94406E+05, 0.99287E+05, 0.10439E+06, 0.10972E+06,
     + 0.11530E+06, 0.12112E+06, 0.12720E+06, 0.13355E+06, 0.14018E+06,
     + 0.14711E+06, 0.15433E+06, 0.16186E+06, 0.16971E+06, 0.17791E+06,
     + 0.18645E+06, 0.19534E+06, 0.20462E+06, 0.21428E+06, 0.22434E+06,
     + 0.23481E+06, 0.24572E+06, 0.25706E+06, 0.26887E+06, 0.28116E+06,
     + 0.29393E+06, 0.30722E+06, 0.32103E+06, 0.33539E+06, 0.35031E+06,
     + 0.36581E+06, 0.38191E+06, 0.39864E+06, 0.41600E+06, 0.43403E+06,
     + 0.45274E+06, 0.47215E+06, 0.49230E+06, 0.51319E+06, 0.53487E+06,
     + 0.55734E+06, 0.58064E+06, 0.60478E+06, 0.62981E+06, 0.65574E+06,
     + 0.68260E+06/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_HNO3  (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 6./
c...     HNO3
c...        --       146
      data (QofT( 1,J),J=1,119)/ 0.15010E+05, 0.25316E+05, 0.37374E+05,
     + 0.51216E+05, 0.67105E+05, 0.85473E+05, 0.10688E+06, 0.13201E+06,
     + 0.16165E+06, 0.19671E+06, 0.23825E+06, 0.28749E+06, 0.34583E+06,
     + 0.41490E+06, 0.49657E+06, 0.59302E+06, 0.70673E+06, 0.84054E+06,
     + 0.99775E+06, 0.11821E+07, 0.13978E+07, 0.16498E+07, 0.19436E+07,
     + 0.22855E+07, 0.26825E+07, 0.31428E+07, 0.36753E+07, 0.42903E+07,
     + 0.49993E+07, 0.58151E+07, 0.67523E+07, 0.78269E+07, 0.90572E+07,
     + 0.10463E+08, 0.12067E+08, 0.13895E+08, 0.15973E+08, 0.18333E+08,
     + 0.21009E+08, 0.24039E+08, 0.27464E+08, 0.31331E+08, 0.35690E+08,
     + 0.40597E+08, 0.46115E+08, 0.52310E+08, 0.59257E+08, 0.67037E+08,
     + 0.75739E+08, 0.85461E+08, 0.96310E+08, 0.10840E+09, 0.12186E+09,
     + 0.13683E+09, 0.15346E+09, 0.17191E+09, 0.19236E+09, 0.21501E+09,
     + 0.24006E+09, 0.26774E+09, 0.29830E+09, 0.33200E+09, 0.36914E+09,
     + 0.41002E+09, 0.45498E+09, 0.50438E+09, 0.55862E+09, 0.61812E+09,
     + 0.68332E+09, 0.75473E+09, 0.83286E+09, 0.91828E+09, 0.10116E+10,
     + 0.11134E+10, 0.12245E+10, 0.13456E+10, 0.14775E+10, 0.16210E+10,
     + 0.17771E+10, 0.19467E+10, 0.21309E+10, 0.23309E+10, 0.25477E+10,
     + 0.27827E+10, 0.30372E+10, 0.33127E+10, 0.36107E+10, 0.39329E+10,
     + 0.42809E+10, 0.46567E+10, 0.50623E+10, 0.54997E+10, 0.59711E+10,
     + 0.64789E+10, 0.70257E+10, 0.76140E+10, 0.82468E+10, 0.89269E+10,
     + 0.96575E+10, 0.10442E+11, 0.11284E+11, 0.12187E+11, 0.13155E+11,
     + 0.14193E+11, 0.15304E+11, 0.16494E+11, 0.17767E+11, 0.19129E+11,
     + 0.20585E+11, 0.22140E+11, 0.23802E+11, 0.25576E+11, 0.27469E+11,
     + 0.29489E+11, 0.31642E+11, 0.33937E+11, 0.36382E+11, 0.38985E+11,
     + 0.41757E+11/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_OH    (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 3), QofT( 3,119),Q(NT)
      data xgj/ 2.,2.,3./
c...       OH
c...        --        61
      data (QofT( 1,J),J=1,119)/ 0.20066E+02, 0.24774E+02, 0.30309E+02,
     + 0.36357E+02, 0.42745E+02, 0.49371E+02, 0.56168E+02, 0.63093E+02,
     + 0.70116E+02, 0.77217E+02, 0.84380E+02, 0.91594E+02, 0.98850E+02,
     + 0.10614E+03, 0.11346E+03, 0.12081E+03, 0.12818E+03, 0.13557E+03,
     + 0.14298E+03, 0.15041E+03, 0.15785E+03, 0.16531E+03, 0.17278E+03,
     + 0.18027E+03, 0.18778E+03, 0.19530E+03, 0.20284E+03, 0.21040E+03,
     + 0.21797E+03, 0.22556E+03, 0.23318E+03, 0.24082E+03, 0.24848E+03,
     + 0.25617E+03, 0.26389E+03, 0.27163E+03, 0.27941E+03, 0.28721E+03,
     + 0.29505E+03, 0.30292E+03, 0.31084E+03, 0.31878E+03, 0.32677E+03,
     + 0.33480E+03, 0.34287E+03, 0.35099E+03, 0.35915E+03, 0.36736E+03,
     + 0.37561E+03, 0.38391E+03, 0.39227E+03, 0.40067E+03, 0.40913E+03,
     + 0.41764E+03, 0.42620E+03, 0.43482E+03, 0.44350E+03, 0.45223E+03,
     + 0.46102E+03, 0.46987E+03, 0.47878E+03, 0.48775E+03, 0.49679E+03,
     + 0.50588E+03, 0.51503E+03, 0.52425E+03, 0.53354E+03, 0.54288E+03,
     + 0.55229E+03, 0.56177E+03, 0.57132E+03, 0.58092E+03, 0.59060E+03,
     + 0.60035E+03, 0.61016E+03, 0.62004E+03, 0.62999E+03, 0.64001E+03,
     + 0.65010E+03, 0.66025E+03, 0.67049E+03, 0.68078E+03, 0.69115E+03,
     + 0.70160E+03, 0.71211E+03, 0.72269E+03, 0.73335E+03, 0.74408E+03,
     + 0.75488E+03, 0.76576E+03, 0.77671E+03, 0.78773E+03, 0.79883E+03,
     + 0.81000E+03, 0.82124E+03, 0.83256E+03, 0.84396E+03, 0.85542E+03,
     + 0.86696E+03, 0.87858E+03, 0.89027E+03, 0.90204E+03, 0.91389E+03,
     + 0.92580E+03, 0.93781E+03, 0.94988E+03, 0.96203E+03, 0.97425E+03,
     + 0.98656E+03, 0.99893E+03, 0.10114E+04, 0.10239E+04, 0.10365E+04,
     + 0.10492E+04, 0.10620E+04, 0.10748E+04, 0.10878E+04, 0.11007E+04,
     + 0.11138E+04/
c...        --        81
      data (QofT( 2,J),J=1,119)/ 0.20124E+02, 0.24876E+02, 0.30457E+02,
     + 0.36553E+02, 0.42991E+02, 0.49666E+02, 0.56513E+02, 0.63489E+02,
     + 0.70563E+02, 0.77715E+02, 0.84929E+02, 0.92195E+02, 0.99504E+02,
     + 0.10685E+03, 0.11423E+03, 0.12164E+03, 0.12907E+03, 0.13654E+03,
     + 0.14403E+03, 0.15154E+03, 0.15909E+03, 0.16666E+03, 0.17427E+03,
     + 0.18191E+03, 0.18959E+03, 0.19731E+03, 0.20507E+03, 0.21287E+03,
     + 0.22073E+03, 0.22863E+03, 0.23658E+03, 0.24459E+03, 0.25266E+03,
     + 0.26078E+03, 0.26897E+03, 0.27722E+03, 0.28554E+03, 0.29393E+03,
     + 0.30238E+03, 0.31091E+03, 0.31952E+03, 0.32820E+03, 0.33696E+03,
     + 0.34579E+03, 0.35471E+03, 0.36371E+03, 0.37279E+03, 0.38196E+03,
     + 0.39121E+03, 0.40055E+03, 0.40998E+03, 0.41949E+03, 0.42910E+03,
     + 0.43879E+03, 0.44858E+03, 0.45845E+03, 0.46843E+03, 0.47849E+03,
     + 0.48865E+03, 0.49890E+03, 0.50924E+03, 0.51969E+03, 0.53022E+03,
     + 0.54086E+03, 0.55159E+03, 0.56242E+03, 0.57335E+03, 0.58437E+03,
     + 0.59550E+03, 0.60673E+03, 0.61805E+03, 0.62947E+03, 0.64100E+03,
     + 0.65263E+03, 0.66435E+03, 0.67618E+03, 0.68811E+03, 0.70014E+03,
     + 0.71228E+03, 0.72451E+03, 0.73685E+03, 0.74929E+03, 0.76184E+03,
     + 0.77449E+03, 0.78724E+03, 0.80009E+03, 0.81306E+03, 0.82612E+03,
     + 0.83929E+03, 0.85256E+03, 0.86594E+03, 0.87942E+03, 0.89301E+03,
     + 0.90670E+03, 0.92050E+03, 0.93440E+03, 0.94841E+03, 0.96253E+03,
     + 0.97675E+03, 0.99108E+03, 0.10055E+04, 0.10201E+04, 0.10347E+04,
     + 0.10495E+04, 0.10643E+04, 0.10793E+04, 0.10944E+04, 0.11096E+04,
     + 0.11248E+04, 0.11402E+04, 0.11558E+04, 0.11714E+04, 0.11871E+04,
     + 0.12029E+04, 0.12189E+04, 0.12349E+04, 0.12511E+04, 0.12673E+04,
     + 0.12837E+04/
c...        --        62
      data (QofT( 3,J),J=1,119)/ 0.41032E+02, 0.54704E+02, 0.70201E+02,
     + 0.86985E+02, 0.10469E+03, 0.12306E+03, 0.14194E+03, 0.16119E+03,
     + 0.18075E+03, 0.20054E+03, 0.22053E+03, 0.24068E+03, 0.26096E+03,
     + 0.28135E+03, 0.30183E+03, 0.32241E+03, 0.34305E+03, 0.36376E+03,
     + 0.38453E+03, 0.40535E+03, 0.42622E+03, 0.44714E+03, 0.46811E+03,
     + 0.48913E+03, 0.51019E+03, 0.53131E+03, 0.55246E+03, 0.57368E+03,
     + 0.59495E+03, 0.61627E+03, 0.63766E+03, 0.65912E+03, 0.68064E+03,
     + 0.70223E+03, 0.72390E+03, 0.74565E+03, 0.76749E+03, 0.78941E+03,
     + 0.81143E+03, 0.83355E+03, 0.85578E+03, 0.87810E+03, 0.90054E+03,
     + 0.92310E+03, 0.94577E+03, 0.96857E+03, 0.99149E+03, 0.10145E+04,
     + 0.10377E+04, 0.10611E+04, 0.10845E+04, 0.11081E+04, 0.11319E+04,
     + 0.11558E+04, 0.11798E+04, 0.12040E+04, 0.12284E+04, 0.12529E+04,
     + 0.12776E+04, 0.13025E+04, 0.13275E+04, 0.13527E+04, 0.13781E+04,
     + 0.14036E+04, 0.14293E+04, 0.14552E+04, 0.14813E+04, 0.15076E+04,
     + 0.15340E+04, 0.15606E+04, 0.15874E+04, 0.16144E+04, 0.16416E+04,
     + 0.16690E+04, 0.16965E+04, 0.17243E+04, 0.17522E+04, 0.17804E+04,
     + 0.18087E+04, 0.18373E+04, 0.18660E+04, 0.18949E+04, 0.19241E+04,
     + 0.19534E+04, 0.19829E+04, 0.20127E+04, 0.20426E+04, 0.20727E+04,
     + 0.21031E+04, 0.21336E+04, 0.21644E+04, 0.21954E+04, 0.22266E+04,
     + 0.22579E+04, 0.22895E+04, 0.23213E+04, 0.23534E+04, 0.23856E+04,
     + 0.24180E+04, 0.24506E+04, 0.24835E+04, 0.25166E+04, 0.25499E+04,
     + 0.25834E+04, 0.26171E+04, 0.26510E+04, 0.26852E+04, 0.27195E+04,
     + 0.27541E+04, 0.27889E+04, 0.28239E+04, 0.28592E+04, 0.28946E+04,
     + 0.29303E+04, 0.29661E+04, 0.30023E+04, 0.30386E+04, 0.30751E+04,
     + 0.31119E+04/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_HF    (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 4./
c...       HF
c...        --        19
      data (QofT( 1,J),J=1,119)/ 0.95958E+01, 0.12933E+02, 0.16295E+02,
     + 0.19666E+02, 0.23043E+02, 0.26425E+02, 0.29809E+02, 0.33195E+02,
     + 0.36584E+02, 0.39974E+02, 0.43366E+02, 0.46759E+02, 0.50154E+02,
     + 0.53550E+02, 0.56947E+02, 0.60346E+02, 0.63746E+02, 0.67148E+02,
     + 0.70550E+02, 0.73955E+02, 0.77361E+02, 0.80769E+02, 0.84179E+02,
     + 0.87591E+02, 0.91006E+02, 0.94424E+02, 0.97846E+02, 0.10127E+03,
     + 0.10470E+03, 0.10813E+03, 0.11157E+03, 0.11502E+03, 0.11847E+03,
     + 0.12193E+03, 0.12540E+03, 0.12888E+03, 0.13236E+03, 0.13586E+03,
     + 0.13936E+03, 0.14288E+03, 0.14641E+03, 0.14995E+03, 0.15351E+03,
     + 0.15708E+03, 0.16066E+03, 0.16426E+03, 0.16788E+03, 0.17151E+03,
     + 0.17516E+03, 0.17882E+03, 0.18251E+03, 0.18621E+03, 0.18994E+03,
     + 0.19368E+03, 0.19745E+03, 0.20123E+03, 0.20504E+03, 0.20887E+03,
     + 0.21272E+03, 0.21659E+03, 0.22049E+03, 0.22441E+03, 0.22836E+03,
     + 0.23233E+03, 0.23632E+03, 0.24034E+03, 0.24439E+03, 0.24846E+03,
     + 0.25255E+03, 0.25668E+03, 0.26083E+03, 0.26501E+03, 0.26921E+03,
     + 0.27344E+03, 0.27770E+03, 0.28199E+03, 0.28631E+03, 0.29066E+03,
     + 0.29503E+03, 0.29944E+03, 0.30387E+03, 0.30833E+03, 0.31282E+03,
     + 0.31735E+03, 0.32190E+03, 0.32648E+03, 0.33110E+03, 0.33574E+03,
     + 0.34042E+03, 0.34512E+03, 0.34986E+03, 0.35463E+03, 0.35943E+03,
     + 0.36426E+03, 0.36913E+03, 0.37402E+03, 0.37895E+03, 0.38391E+03,
     + 0.38891E+03, 0.39393E+03, 0.39899E+03, 0.40408E+03, 0.40921E+03,
     + 0.41436E+03, 0.41955E+03, 0.42478E+03, 0.43004E+03, 0.43533E+03,
     + 0.44065E+03, 0.44601E+03, 0.45140E+03, 0.45683E+03, 0.46229E+03,
     + 0.46779E+03, 0.47332E+03, 0.47888E+03, 0.48448E+03, 0.49011E+03,
     + 0.49578E+03/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_HCl   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 2), QofT( 2,119),Q(NT)
      data xgj/ 8.,8./
c...      HCl
c...        --        15
      data (QofT( 1,J),J=1,119)/ 0.34775E+02, 0.48060E+02, 0.61370E+02,
     + 0.74692E+02, 0.88024E+02, 0.10136E+03, 0.11471E+03, 0.12806E+03,
     + 0.14141E+03, 0.15478E+03, 0.16814E+03, 0.18151E+03, 0.19489E+03,
     + 0.20827E+03, 0.22166E+03, 0.23506E+03, 0.24847E+03, 0.26189E+03,
     + 0.27533E+03, 0.28878E+03, 0.30225E+03, 0.31575E+03, 0.32928E+03,
     + 0.34284E+03, 0.35645E+03, 0.37009E+03, 0.38378E+03, 0.39753E+03,
     + 0.41134E+03, 0.42521E+03, 0.43914E+03, 0.45316E+03, 0.46725E+03,
     + 0.48142E+03, 0.49568E+03, 0.51003E+03, 0.52448E+03, 0.53902E+03,
     + 0.55368E+03, 0.56843E+03, 0.58330E+03, 0.59829E+03, 0.61339E+03,
     + 0.62862E+03, 0.64396E+03, 0.65944E+03, 0.67504E+03, 0.69078E+03,
     + 0.70665E+03, 0.72265E+03, 0.73880E+03, 0.75508E+03, 0.77151E+03,
     + 0.78809E+03, 0.80481E+03, 0.82168E+03, 0.83870E+03, 0.85587E+03,
     + 0.87320E+03, 0.89068E+03, 0.90832E+03, 0.92611E+03, 0.94407E+03,
     + 0.96218E+03, 0.98046E+03, 0.99889E+03, 0.10175E+04, 0.10363E+04,
     + 0.10552E+04, 0.10743E+04, 0.10936E+04, 0.11130E+04, 0.11326E+04,
     + 0.11524E+04, 0.11723E+04, 0.11924E+04, 0.12127E+04, 0.12332E+04,
     + 0.12538E+04, 0.12746E+04, 0.12956E+04, 0.13168E+04, 0.13381E+04,
     + 0.13597E+04, 0.13814E+04, 0.14032E+04, 0.14253E+04, 0.14475E+04,
     + 0.14700E+04, 0.14926E+04, 0.15153E+04, 0.15383E+04, 0.15615E+04,
     + 0.15848E+04, 0.16083E+04, 0.16320E+04, 0.16559E+04, 0.16800E+04,
     + 0.17043E+04, 0.17287E+04, 0.17533E+04, 0.17782E+04, 0.18032E+04,
     + 0.18284E+04, 0.18538E+04, 0.18794E+04, 0.19051E+04, 0.19311E+04,
     + 0.19573E+04, 0.19836E+04, 0.20102E+04, 0.20369E+04, 0.20638E+04,
     + 0.20910E+04, 0.21183E+04, 0.21458E+04, 0.21735E+04, 0.22014E+04,
     + 0.22295E+04/
c...        --        17
      data (QofT( 2,J),J=1,119)/ 0.34823E+02, 0.48128E+02, 0.61458E+02,
     + 0.74801E+02, 0.88152E+02, 0.10151E+03, 0.11488E+03, 0.12825E+03,
     + 0.14162E+03, 0.15500E+03, 0.16839E+03, 0.18178E+03, 0.19518E+03,
     + 0.20858E+03, 0.22199E+03, 0.23541E+03, 0.24884E+03, 0.26228E+03,
     + 0.27574E+03, 0.28921E+03, 0.30270E+03, 0.31622E+03, 0.32977E+03,
     + 0.34336E+03, 0.35698E+03, 0.37065E+03, 0.38436E+03, 0.39813E+03,
     + 0.41196E+03, 0.42585E+03, 0.43981E+03, 0.45384E+03, 0.46796E+03,
     + 0.48215E+03, 0.49644E+03, 0.51081E+03, 0.52528E+03, 0.53986E+03,
     + 0.55453E+03, 0.56932E+03, 0.58421E+03, 0.59922E+03, 0.61435E+03,
     + 0.62960E+03, 0.64498E+03, 0.66048E+03, 0.67611E+03, 0.69187E+03,
     + 0.70777E+03, 0.72381E+03, 0.73998E+03, 0.75630E+03, 0.77276E+03,
     + 0.78936E+03, 0.80612E+03, 0.82302E+03, 0.84007E+03, 0.85727E+03,
     + 0.87463E+03, 0.89215E+03, 0.90982E+03, 0.92765E+03, 0.94563E+03,
     + 0.96378E+03, 0.98209E+03, 0.10006E+04, 0.10192E+04, 0.10380E+04,
     + 0.10570E+04, 0.10761E+04, 0.10954E+04, 0.11149E+04, 0.11345E+04,
     + 0.11543E+04, 0.11743E+04, 0.11945E+04, 0.12148E+04, 0.12353E+04,
     + 0.12560E+04, 0.12768E+04, 0.12979E+04, 0.13191E+04, 0.13405E+04,
     + 0.13620E+04, 0.13838E+04, 0.14057E+04, 0.14278E+04, 0.14501E+04,
     + 0.14726E+04, 0.14952E+04, 0.15180E+04, 0.15410E+04, 0.15642E+04,
     + 0.15876E+04, 0.16112E+04, 0.16349E+04, 0.16589E+04, 0.16830E+04,
     + 0.17073E+04, 0.17318E+04, 0.17565E+04, 0.17814E+04, 0.18064E+04,
     + 0.18317E+04, 0.18572E+04, 0.18828E+04, 0.19086E+04, 0.19346E+04,
     + 0.19609E+04, 0.19873E+04, 0.20139E+04, 0.20406E+04, 0.20676E+04,
     + 0.20948E+04, 0.21222E+04, 0.21498E+04, 0.21775E+04, 0.22055E+04,
     + 0.22337E+04/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_HBr   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 2), QofT( 2,119),Q(NT)
      data xgj/ 8.,8./
c...      HBr
c...        --        19
      data (QofT( 1,J),J=1,119)/ 0.42744E+02, 0.59373E+02, 0.76023E+02,
     + 0.92685E+02, 0.10936E+03, 0.12604E+03, 0.14272E+03, 0.15942E+03,
     + 0.17612E+03, 0.19282E+03, 0.20954E+03, 0.22626E+03, 0.24299E+03,
     + 0.25973E+03, 0.27648E+03, 0.29325E+03, 0.31004E+03, 0.32686E+03,
     + 0.34371E+03, 0.36060E+03, 0.37753E+03, 0.39451E+03, 0.41156E+03,
     + 0.42868E+03, 0.44587E+03, 0.46314E+03, 0.48051E+03, 0.49798E+03,
     + 0.51556E+03, 0.53325E+03, 0.55106E+03, 0.56900E+03, 0.58708E+03,
     + 0.60530E+03, 0.62367E+03, 0.64219E+03, 0.66088E+03, 0.67972E+03,
     + 0.69874E+03, 0.71793E+03, 0.73730E+03, 0.75685E+03, 0.77659E+03,
     + 0.79652E+03, 0.81664E+03, 0.83696E+03, 0.85748E+03, 0.87820E+03,
     + 0.89914E+03, 0.92028E+03, 0.94163E+03, 0.96319E+03, 0.98498E+03,
     + 0.10070E+04, 0.10292E+04, 0.10516E+04, 0.10743E+04, 0.10972E+04,
     + 0.11203E+04, 0.11437E+04, 0.11673E+04, 0.11911E+04, 0.12151E+04,
     + 0.12394E+04, 0.12640E+04, 0.12887E+04, 0.13137E+04, 0.13390E+04,
     + 0.13645E+04, 0.13902E+04, 0.14162E+04, 0.14424E+04, 0.14689E+04,
     + 0.14956E+04, 0.15226E+04, 0.15498E+04, 0.15773E+04, 0.16050E+04,
     + 0.16330E+04, 0.16612E+04, 0.16897E+04, 0.17185E+04, 0.17475E+04,
     + 0.17767E+04, 0.18062E+04, 0.18360E+04, 0.18660E+04, 0.18963E+04,
     + 0.19269E+04, 0.19577E+04, 0.19888E+04, 0.20202E+04, 0.20518E+04,
     + 0.20837E+04, 0.21158E+04, 0.21482E+04, 0.21809E+04, 0.22139E+04,
     + 0.22471E+04, 0.22806E+04, 0.23143E+04, 0.23484E+04, 0.23827E+04,
     + 0.24173E+04, 0.24521E+04, 0.24873E+04, 0.25227E+04, 0.25584E+04,
     + 0.25943E+04, 0.26306E+04, 0.26671E+04, 0.27039E+04, 0.27409E+04,
     + 0.27783E+04, 0.28159E+04, 0.28538E+04, 0.28920E+04, 0.29305E+04,
     + 0.29693E+04/
c...        --        11
      data (QofT( 2,J),J=1,119)/ 0.42756E+02, 0.59390E+02, 0.76045E+02,
     + 0.92713E+02, 0.10939E+03, 0.12607E+03, 0.14277E+03, 0.15947E+03,
     + 0.17617E+03, 0.19288E+03, 0.20960E+03, 0.22633E+03, 0.24306E+03,
     + 0.25981E+03, 0.27656E+03, 0.29334E+03, 0.31014E+03, 0.32696E+03,
     + 0.34381E+03, 0.36071E+03, 0.37764E+03, 0.39464E+03, 0.41169E+03,
     + 0.42881E+03, 0.44601E+03, 0.46329E+03, 0.48066E+03, 0.49813E+03,
     + 0.51572E+03, 0.53341E+03, 0.55123E+03, 0.56918E+03, 0.58727E+03,
     + 0.60549E+03, 0.62387E+03, 0.64240E+03, 0.66109E+03, 0.67994E+03,
     + 0.69896E+03, 0.71816E+03, 0.73754E+03, 0.75710E+03, 0.77684E+03,
     + 0.79678E+03, 0.81691E+03, 0.83724E+03, 0.85776E+03, 0.87850E+03,
     + 0.89943E+03, 0.92058E+03, 0.94194E+03, 0.96352E+03, 0.98531E+03,
     + 0.10073E+04, 0.10295E+04, 0.10520E+04, 0.10747E+04, 0.10976E+04,
     + 0.11207E+04, 0.11441E+04, 0.11677E+04, 0.11915E+04, 0.12156E+04,
     + 0.12399E+04, 0.12644E+04, 0.12892E+04, 0.13142E+04, 0.13395E+04,
     + 0.13650E+04, 0.13907E+04, 0.14167E+04, 0.14429E+04, 0.14694E+04,
     + 0.14961E+04, 0.15231E+04, 0.15504E+04, 0.15778E+04, 0.16056E+04,
     + 0.16336E+04, 0.16618E+04, 0.16903E+04, 0.17191E+04, 0.17481E+04,
     + 0.17773E+04, 0.18069E+04, 0.18367E+04, 0.18667E+04, 0.18970E+04,
     + 0.19276E+04, 0.19584E+04, 0.19895E+04, 0.20209E+04, 0.20525E+04,
     + 0.20844E+04, 0.21166E+04, 0.21490E+04, 0.21817E+04, 0.22147E+04,
     + 0.22479E+04, 0.22814E+04, 0.23152E+04, 0.23492E+04, 0.23835E+04,
     + 0.24181E+04, 0.24530E+04, 0.24882E+04, 0.25236E+04, 0.25593E+04,
     + 0.25952E+04, 0.26315E+04, 0.26680E+04, 0.27048E+04, 0.27419E+04,
     + 0.27793E+04, 0.28169E+04, 0.28549E+04, 0.28931E+04, 0.29316E+04,
     + 0.29703E+04/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_HI    (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 12./
c...       HI
c...        --        17
      data (QofT( 1,J),J=1,119)/ 0.82031E+02, 0.11447E+03, 0.14694E+03,
     + 0.17943E+03, 0.21194E+03, 0.24445E+03, 0.27699E+03, 0.30953E+03,
     + 0.34209E+03, 0.37466E+03, 0.40725E+03, 0.43986E+03, 0.47249E+03,
     + 0.50517E+03, 0.53789E+03, 0.57068E+03, 0.60354E+03, 0.63650E+03,
     + 0.66957E+03, 0.70278E+03, 0.73614E+03, 0.76967E+03, 0.80340E+03,
     + 0.83735E+03, 0.87153E+03, 0.90596E+03, 0.94067E+03, 0.97566E+03,
     + 0.10110E+04, 0.10466E+04, 0.10826E+04, 0.11189E+04, 0.11555E+04,
     + 0.11926E+04, 0.12300E+04, 0.12679E+04, 0.13061E+04, 0.13448E+04,
     + 0.13839E+04, 0.14235E+04, 0.14635E+04, 0.15039E+04, 0.15448E+04,
     + 0.15862E+04, 0.16280E+04, 0.16704E+04, 0.17132E+04, 0.17565E+04,
     + 0.18003E+04, 0.18446E+04, 0.18894E+04, 0.19347E+04, 0.19806E+04,
     + 0.20269E+04, 0.20738E+04, 0.21212E+04, 0.21691E+04, 0.22176E+04,
     + 0.22666E+04, 0.23162E+04, 0.23662E+04, 0.24169E+04, 0.24680E+04,
     + 0.25198E+04, 0.25720E+04, 0.26249E+04, 0.26783E+04, 0.27322E+04,
     + 0.27867E+04, 0.28418E+04, 0.28975E+04, 0.29537E+04, 0.30105E+04,
     + 0.30678E+04, 0.31258E+04, 0.31843E+04, 0.32434E+04, 0.33031E+04,
     + 0.33633E+04, 0.34242E+04, 0.34856E+04, 0.35477E+04, 0.36103E+04,
     + 0.36735E+04, 0.37373E+04, 0.38018E+04, 0.38668E+04, 0.39324E+04,
     + 0.39986E+04, 0.40654E+04, 0.41329E+04, 0.42009E+04, 0.42696E+04,
     + 0.43388E+04, 0.44087E+04, 0.44792E+04, 0.45503E+04, 0.46221E+04,
     + 0.46944E+04, 0.47674E+04, 0.48410E+04, 0.49152E+04, 0.49901E+04,
     + 0.50656E+04, 0.51417E+04, 0.52185E+04, 0.52959E+04, 0.53739E+04,
     + 0.54526E+04, 0.55319E+04, 0.56118E+04, 0.56924E+04, 0.57736E+04,
     + 0.58555E+04, 0.59380E+04, 0.60212E+04, 0.61050E+04, 0.61895E+04,
     + 0.62746E+04/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_ClO   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 2), QofT( 2,119),Q(NT)
      data xgj/ 4.,4./
c...      ClO
c...        --        56
      data (QofT( 1,J),J=1,119)/ 0.53847E+03, 0.76580E+03, 0.10017E+04,
     + 0.12511E+04, 0.15168E+04, 0.18001E+04, 0.21014E+04, 0.24206E+04,
     + 0.27577E+04, 0.31127E+04, 0.34857E+04, 0.38765E+04, 0.42854E+04,
     + 0.47124E+04, 0.51575E+04, 0.56208E+04, 0.61025E+04, 0.66026E+04,
     + 0.71211E+04, 0.76582E+04, 0.82138E+04, 0.87882E+04, 0.93813E+04,
     + 0.99932E+04, 0.10624E+05, 0.11273E+05, 0.11942E+05, 0.12629E+05,
     + 0.13336E+05, 0.14061E+05, 0.14806E+05, 0.15570E+05, 0.16353E+05,
     + 0.17155E+05, 0.17976E+05, 0.18816E+05, 0.19676E+05, 0.20555E+05,
     + 0.21453E+05, 0.22371E+05, 0.23308E+05, 0.24264E+05, 0.25240E+05,
     + 0.26236E+05, 0.27250E+05, 0.28284E+05, 0.29338E+05, 0.30412E+05,
     + 0.31505E+05, 0.32617E+05, 0.33749E+05, 0.34901E+05, 0.36072E+05,
     + 0.37263E+05, 0.38474E+05, 0.39705E+05, 0.40955E+05, 0.42225E+05,
     + 0.43515E+05, 0.44825E+05, 0.46154E+05, 0.47504E+05, 0.48873E+05,
     + 0.50262E+05, 0.51672E+05, 0.53101E+05, 0.54549E+05, 0.56019E+05,
     + 0.57508E+05, 0.59017E+05, 0.60546E+05, 0.62095E+05, 0.63665E+05,
     + 0.65254E+05, 0.66864E+05, 0.68494E+05, 0.70144E+05, 0.71814E+05,
     + 0.73504E+05, 0.75215E+05, 0.76946E+05, 0.78698E+05, 0.80470E+05,
     + 0.82261E+05, 0.84074E+05, 0.85907E+05, 0.87760E+05, 0.89633E+05,
     + 0.91527E+05, 0.93442E+05, 0.95377E+05, 0.97333E+05, 0.99309E+05,
     + 0.10131E+06, 0.10332E+06, 0.10536E+06, 0.10742E+06, 0.10950E+06,
     + 0.11160E+06, 0.11372E+06, 0.11586E+06, 0.11802E+06, 0.12020E+06,
     + 0.12241E+06, 0.12463E+06, 0.12688E+06, 0.12914E+06, 0.13143E+06,
     + 0.13374E+06, 0.13607E+06, 0.13842E+06, 0.14079E+06, 0.14318E+06,
     + 0.14559E+06, 0.14802E+06, 0.15048E+06, 0.15295E+06, 0.15545E+06,
     + 0.15797E+06/
c...        --        76
      data (QofT( 2,J),J=1,119)/ 0.54775E+03, 0.77899E+03, 0.10189E+04,
     + 0.12726E+04, 0.15430E+04, 0.18313E+04, 0.21378E+04, 0.24627E+04,
     + 0.28059E+04, 0.31674E+04, 0.35472E+04, 0.39454E+04, 0.43621E+04,
     + 0.47972E+04, 0.52508E+04, 0.57232E+04, 0.62143E+04, 0.67242E+04,
     + 0.72531E+04, 0.78010E+04, 0.83678E+04, 0.89537E+04, 0.95589E+04,
     + 0.10183E+05, 0.10827E+05, 0.11490E+05, 0.12172E+05, 0.12874E+05,
     + 0.13595E+05, 0.14335E+05, 0.15095E+05, 0.15875E+05, 0.16674E+05,
     + 0.17493E+05, 0.18332E+05, 0.19190E+05, 0.20068E+05, 0.20965E+05,
     + 0.21882E+05, 0.22820E+05, 0.23776E+05, 0.24753E+05, 0.25750E+05,
     + 0.26766E+05, 0.27803E+05, 0.28859E+05, 0.29935E+05, 0.31032E+05,
     + 0.32148E+05, 0.33284E+05, 0.34441E+05, 0.35617E+05, 0.36814E+05,
     + 0.38031E+05, 0.39267E+05, 0.40524E+05, 0.41802E+05, 0.43099E+05,
     + 0.44417E+05, 0.45755E+05, 0.47113E+05, 0.48492E+05, 0.49891E+05,
     + 0.51310E+05, 0.52750E+05, 0.54210E+05, 0.55690E+05, 0.57191E+05,
     + 0.58713E+05, 0.60255E+05, 0.61817E+05, 0.63400E+05, 0.65004E+05,
     + 0.66628E+05, 0.68272E+05, 0.69938E+05, 0.71624E+05, 0.73331E+05,
     + 0.75058E+05, 0.76806E+05, 0.78575E+05, 0.80364E+05, 0.82175E+05,
     + 0.84006E+05, 0.85858E+05, 0.87731E+05, 0.89625E+05, 0.91539E+05,
     + 0.93475E+05, 0.95431E+05, 0.97409E+05, 0.99407E+05, 0.10143E+06,
     + 0.10347E+06, 0.10553E+06, 0.10761E+06, 0.10972E+06, 0.11184E+06,
     + 0.11399E+06, 0.11615E+06, 0.11834E+06, 0.12055E+06, 0.12278E+06,
     + 0.12503E+06, 0.12731E+06, 0.12960E+06, 0.13192E+06, 0.13425E+06,
     + 0.13661E+06, 0.13899E+06, 0.14139E+06, 0.14382E+06, 0.14626E+06,
     + 0.14873E+06, 0.15121E+06, 0.15372E+06, 0.15625E+06, 0.15880E+06,
     + 0.16138E+06/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_OCS   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 5), QofT( 5,119),Q(NT)
      data xgj/ 1.,1.,2.,4.,1./
c...      OCS
c...        --       622
      data (QofT( 1,J),J=1,119)/ 0.20609E+03, 0.29199E+03, 0.37861E+03,
     + 0.46737E+03, 0.56024E+03, 0.65929E+03, 0.76649E+03, 0.88361E+03,
     + 0.10123E+04, 0.11541E+04, 0.13105E+04, 0.14829E+04, 0.16728E+04,
     + 0.18818E+04, 0.21113E+04, 0.23629E+04, 0.26383E+04, 0.29391E+04,
     + 0.32672E+04, 0.36245E+04, 0.40128E+04, 0.44343E+04, 0.48911E+04,
     + 0.53853E+04, 0.59193E+04, 0.64956E+04, 0.71166E+04, 0.77849E+04,
     + 0.85033E+04, 0.92746E+04, 0.10102E+05, 0.10988E+05, 0.11936E+05,
     + 0.12949E+05, 0.14032E+05, 0.15186E+05, 0.16416E+05, 0.17726E+05,
     + 0.19120E+05, 0.20601E+05, 0.22173E+05, 0.23842E+05, 0.25611E+05,
     + 0.27484E+05, 0.29468E+05, 0.31566E+05, 0.33783E+05, 0.36124E+05,
     + 0.38595E+05, 0.41202E+05, 0.43949E+05, 0.46842E+05, 0.49888E+05,
     + 0.53092E+05, 0.56460E+05, 0.59999E+05, 0.63716E+05, 0.67616E+05,
     + 0.71708E+05, 0.75997E+05, 0.80491E+05, 0.85197E+05, 0.90124E+05,
     + 0.95278E+05, 0.10067E+06, 0.10630E+06, 0.11219E+06, 0.11833E+06,
     + 0.12475E+06, 0.13144E+06, 0.13842E+06, 0.14570E+06, 0.15328E+06,
     + 0.16117E+06, 0.16940E+06, 0.17795E+06, 0.18686E+06, 0.19611E+06,
     + 0.20574E+06, 0.21574E+06, 0.22613E+06, 0.23692E+06, 0.24813E+06,
     + 0.25975E+06, 0.27182E+06, 0.28433E+06, 0.29730E+06, 0.31074E+06,
     + 0.32467E+06, 0.33909E+06, 0.35403E+06, 0.36950E+06, 0.38551E+06,
     + 0.40207E+06, 0.41920E+06, 0.43691E+06, 0.45522E+06, 0.47415E+06,
     + 0.49370E+06, 0.51390E+06, 0.53476E+06, 0.55629E+06, 0.57852E+06,
     + 0.60146E+06, 0.62513E+06, 0.64954E+06, 0.67471E+06, 0.70067E+06,
     + 0.72742E+06, 0.75499E+06, 0.78339E+06, 0.81265E+06, 0.84279E+06,
     + 0.87381E+06, 0.90576E+06, 0.93863E+06, 0.97246E+06, 0.10073E+07,
     + 0.10431E+07/
c...        --       624
      data (QofT( 2,J),J=1,119)/ 0.21125E+03, 0.29930E+03, 0.38809E+03,
     + 0.47911E+03, 0.57437E+03, 0.67603E+03, 0.78610E+03, 0.90643E+03,
     + 0.10387E+04, 0.11846E+04, 0.13456E+04, 0.15231E+04, 0.17188E+04,
     + 0.19342E+04, 0.21709E+04, 0.24304E+04, 0.27145E+04, 0.30250E+04,
     + 0.33638E+04, 0.37328E+04, 0.41339E+04, 0.45694E+04, 0.50415E+04,
     + 0.55524E+04, 0.61045E+04, 0.67004E+04, 0.73427E+04, 0.80340E+04,
     + 0.87773E+04, 0.95755E+04, 0.10432E+05, 0.11349E+05, 0.12330E+05,
     + 0.13380E+05, 0.14500E+05, 0.15696E+05, 0.16970E+05, 0.18327E+05,
     + 0.19770E+05, 0.21305E+05, 0.22934E+05, 0.24663E+05, 0.26497E+05,
     + 0.28439E+05, 0.30495E+05, 0.32669E+05, 0.34968E+05, 0.37396E+05,
     + 0.39958E+05, 0.42661E+05, 0.45510E+05, 0.48511E+05, 0.51669E+05,
     + 0.54993E+05, 0.58487E+05, 0.62159E+05, 0.66014E+05, 0.70061E+05,
     + 0.74306E+05, 0.78757E+05, 0.83421E+05, 0.88305E+05, 0.93418E+05,
     + 0.98767E+05, 0.10436E+06, 0.11021E+06, 0.11632E+06, 0.12270E+06,
     + 0.12936E+06, 0.13631E+06, 0.14355E+06, 0.15111E+06, 0.15898E+06,
     + 0.16718E+06, 0.17572E+06, 0.18460E+06, 0.19385E+06, 0.20346E+06,
     + 0.21346E+06, 0.22385E+06, 0.23464E+06, 0.24585E+06, 0.25748E+06,
     + 0.26956E+06, 0.28209E+06, 0.29509E+06, 0.30856E+06, 0.32252E+06,
     + 0.33699E+06, 0.35198E+06, 0.36750E+06, 0.38357E+06, 0.40020E+06,
     + 0.41741E+06, 0.43521E+06, 0.45362E+06, 0.47264E+06, 0.49231E+06,
     + 0.51263E+06, 0.53362E+06, 0.55529E+06, 0.57768E+06, 0.60078E+06,
     + 0.62462E+06, 0.64922E+06, 0.67459E+06, 0.70075E+06, 0.72773E+06,
     + 0.75554E+06, 0.78419E+06, 0.81372E+06, 0.84413E+06, 0.87546E+06,
     + 0.90771E+06, 0.94092E+06, 0.97509E+06, 0.10103E+07, 0.10464E+07,
     + 0.10837E+07/
c...        --       632
      data (QofT( 3,J),J=1,119)/ 0.41351E+03, 0.58591E+03, 0.76004E+03,
     + 0.93907E+03, 0.11273E+04, 0.13289E+04, 0.15481E+04, 0.17884E+04,
     + 0.20533E+04, 0.23459E+04, 0.26692E+04, 0.30264E+04, 0.34205E+04,
     + 0.38547E+04, 0.43323E+04, 0.48565E+04, 0.54309E+04, 0.60592E+04,
     + 0.67451E+04, 0.74928E+04, 0.83064E+04, 0.91903E+04, 0.10149E+05,
     + 0.11187E+05, 0.12310E+05, 0.13523E+05, 0.14831E+05, 0.16240E+05,
     + 0.17756E+05, 0.19384E+05, 0.21132E+05, 0.23005E+05, 0.25011E+05,
     + 0.27157E+05, 0.29449E+05, 0.31896E+05, 0.34506E+05, 0.37286E+05,
     + 0.40245E+05, 0.43392E+05, 0.46735E+05, 0.50284E+05, 0.54048E+05,
     + 0.58038E+05, 0.62263E+05, 0.66733E+05, 0.71460E+05, 0.76455E+05,
     + 0.81728E+05, 0.87292E+05, 0.93159E+05, 0.99341E+05, 0.10585E+06,
     + 0.11270E+06, 0.11991E+06, 0.12748E+06, 0.13543E+06, 0.14378E+06,
     + 0.15255E+06, 0.16174E+06, 0.17137E+06, 0.18146E+06, 0.19202E+06,
     + 0.20308E+06, 0.21465E+06, 0.22674E+06, 0.23937E+06, 0.25257E+06,
     + 0.26635E+06, 0.28073E+06, 0.29573E+06, 0.31137E+06, 0.32767E+06,
     + 0.34466E+06, 0.36235E+06, 0.38076E+06, 0.39992E+06, 0.41985E+06,
     + 0.44057E+06, 0.46211E+06, 0.48450E+06, 0.50775E+06, 0.53189E+06,
     + 0.55695E+06, 0.58295E+06, 0.60992E+06, 0.63789E+06, 0.66688E+06,
     + 0.69693E+06, 0.72806E+06, 0.76030E+06, 0.79368E+06, 0.82823E+06,
     + 0.86399E+06, 0.90097E+06, 0.93923E+06, 0.97878E+06, 0.10197E+07,
     + 0.10619E+07, 0.11056E+07, 0.11506E+07, 0.11972E+07, 0.12453E+07,
     + 0.12949E+07, 0.13460E+07, 0.13988E+07, 0.14533E+07, 0.15094E+07,
     + 0.15673E+07, 0.16270E+07, 0.16884E+07, 0.17518E+07, 0.18170E+07,
     + 0.18842E+07, 0.19533E+07, 0.20245E+07, 0.20978E+07, 0.21732E+07,
     + 0.22507E+07/
c...        --       623
      data (QofT( 4,J),J=1,119)/ 0.83485E+03, 0.11828E+04, 0.15337E+04,
     + 0.18934E+04, 0.22697E+04, 0.26712E+04, 0.31059E+04, 0.35809E+04,
     + 0.41030E+04, 0.46785E+04, 0.53133E+04, 0.60135E+04, 0.67850E+04,
     + 0.76338E+04, 0.85663E+04, 0.95888E+04, 0.10708E+05, 0.11931E+05,
     + 0.13265E+05, 0.14718E+05, 0.16298E+05, 0.18012E+05, 0.19870E+05,
     + 0.21881E+05, 0.24054E+05, 0.26399E+05, 0.28926E+05, 0.31646E+05,
     + 0.34570E+05, 0.37710E+05, 0.41077E+05, 0.44685E+05, 0.48545E+05,
     + 0.52672E+05, 0.57078E+05, 0.61780E+05, 0.66790E+05, 0.72125E+05,
     + 0.77801E+05, 0.83833E+05, 0.90239E+05, 0.97036E+05, 0.10424E+06,
     + 0.11188E+06, 0.11996E+06, 0.12850E+06, 0.13754E+06, 0.14708E+06,
     + 0.15715E+06, 0.16777E+06, 0.17896E+06, 0.19076E+06, 0.20317E+06,
     + 0.21623E+06, 0.22996E+06, 0.24438E+06, 0.25953E+06, 0.27543E+06,
     + 0.29211E+06, 0.30959E+06, 0.32791E+06, 0.34710E+06, 0.36718E+06,
     + 0.38820E+06, 0.41017E+06, 0.43314E+06, 0.45713E+06, 0.48219E+06,
     + 0.50835E+06, 0.53564E+06, 0.56409E+06, 0.59376E+06, 0.62468E+06,
     + 0.65688E+06, 0.69041E+06, 0.72530E+06, 0.76161E+06, 0.79937E+06,
     + 0.83862E+06, 0.87941E+06, 0.92179E+06, 0.96581E+06, 0.10115E+07,
     + 0.10589E+07, 0.11081E+07, 0.11591E+07, 0.12120E+07, 0.12669E+07,
     + 0.13237E+07, 0.13825E+07, 0.14435E+07, 0.15066E+07, 0.15718E+07,
     + 0.16394E+07, 0.17093E+07, 0.17815E+07, 0.18562E+07, 0.19334E+07,
     + 0.20132E+07, 0.20956E+07, 0.21807E+07, 0.22685E+07, 0.23592E+07,
     + 0.24528E+07, 0.25494E+07, 0.26490E+07, 0.27517E+07, 0.28576E+07,
     + 0.29667E+07, 0.30792E+07, 0.31951E+07, 0.33145E+07, 0.34374E+07,
     + 0.35640E+07, 0.36943E+07, 0.38285E+07, 0.39665E+07, 0.41085E+07,
     + 0.42546E+07/
c...        --       822
      data (QofT( 5,J),J=1,119)/ 0.21967E+03, 0.31126E+03, 0.40370E+03,
     + 0.49862E+03, 0.59823E+03, 0.70481E+03, 0.82050E+03, 0.94724E+03,
     + 0.10868E+04, 0.12409E+04, 0.14112E+04, 0.15993E+04, 0.18067E+04,
     + 0.20353E+04, 0.22866E+04, 0.25624E+04, 0.28645E+04, 0.31950E+04,
     + 0.35558E+04, 0.39490E+04, 0.43767E+04, 0.48413E+04, 0.53452E+04,
     + 0.58909E+04, 0.64810E+04, 0.71182E+04, 0.78053E+04, 0.85454E+04,
     + 0.93413E+04, 0.10196E+05, 0.11114E+05, 0.12098E+05, 0.13151E+05,
     + 0.14277E+05, 0.15480E+05, 0.16764E+05, 0.18133E+05, 0.19592E+05,
     + 0.21144E+05, 0.22794E+05, 0.24548E+05, 0.26409E+05, 0.28383E+05,
     + 0.30475E+05, 0.32689E+05, 0.35033E+05, 0.37511E+05, 0.40128E+05,
     + 0.42892E+05, 0.45808E+05, 0.48882E+05, 0.52121E+05, 0.55532E+05,
     + 0.59121E+05, 0.62895E+05, 0.66861E+05, 0.71028E+05, 0.75402E+05,
     + 0.79991E+05, 0.84803E+05, 0.89847E+05, 0.95130E+05, 0.10066E+06,
     + 0.10645E+06, 0.11251E+06, 0.11883E+06, 0.12545E+06, 0.13236E+06,
     + 0.13957E+06, 0.14710E+06, 0.15495E+06, 0.16313E+06, 0.17166E+06,
     + 0.18055E+06, 0.18980E+06, 0.19944E+06, 0.20946E+06, 0.21989E+06,
     + 0.23073E+06, 0.24200E+06, 0.25371E+06, 0.26587E+06, 0.27850E+06,
     + 0.29161E+06, 0.30521E+06, 0.31931E+06, 0.33394E+06, 0.34910E+06,
     + 0.36482E+06, 0.38109E+06, 0.39795E+06, 0.41541E+06, 0.43348E+06,
     + 0.45217E+06, 0.47151E+06, 0.49151E+06, 0.51219E+06, 0.53356E+06,
     + 0.55565E+06, 0.57847E+06, 0.60204E+06, 0.62637E+06, 0.65149E+06,
     + 0.67742E+06, 0.70417E+06, 0.73176E+06, 0.76023E+06, 0.78957E+06,
     + 0.81982E+06, 0.85100E+06, 0.88313E+06, 0.91622E+06, 0.95031E+06,
     + 0.98541E+06, 0.10216E+07, 0.10587E+07, 0.10970E+07, 0.11364E+07,
     + 0.11769E+07/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_H2CO  (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 3), QofT( 3,119),Q(NT)
      data xgj/ 1.,2.,1./
c...     H2CO
c...        --       126
      data (QofT( 1,J),J=1,119)/ 0.25934E+03, 0.43623E+03, 0.64143E+03,
     + 0.87152E+03, 0.11241E+04, 0.13975E+04, 0.16906E+04, 0.20029E+04,
     + 0.23344E+04, 0.26857E+04, 0.30577E+04, 0.34518E+04, 0.38698E+04,
     + 0.43138E+04, 0.47860E+04, 0.52890E+04, 0.58256E+04, 0.63985E+04,
     + 0.70109E+04, 0.76660E+04, 0.83673E+04, 0.91184E+04, 0.99230E+04,
     + 0.10785E+05, 0.11710E+05, 0.12700E+05, 0.13762E+05, 0.14900E+05,
     + 0.16119E+05, 0.17425E+05, 0.18823E+05, 0.20320E+05, 0.21923E+05,
     + 0.23637E+05, 0.25471E+05, 0.27432E+05, 0.29527E+05, 0.31765E+05,
     + 0.34155E+05, 0.36706E+05, 0.39428E+05, 0.42330E+05, 0.45424E+05,
     + 0.48720E+05, 0.52231E+05, 0.55968E+05, 0.59945E+05, 0.64175E+05,
     + 0.68672E+05, 0.73450E+05, 0.78526E+05, 0.83915E+05, 0.89634E+05,
     + 0.95701E+05, 0.10213E+06, 0.10895E+06, 0.11618E+06, 0.12383E+06,
     + 0.13193E+06, 0.14049E+06, 0.14956E+06, 0.15914E+06, 0.16927E+06,
     + 0.17997E+06, 0.19127E+06, 0.20320E+06, 0.21578E+06, 0.22906E+06,
     + 0.24306E+06, 0.25782E+06, 0.27336E+06, 0.28974E+06, 0.30698E+06,
     + 0.32513E+06, 0.34422E+06, 0.36430E+06, 0.38542E+06, 0.40761E+06,
     + 0.43093E+06, 0.45542E+06, 0.48114E+06, 0.50813E+06, 0.53646E+06,
     + 0.56617E+06, 0.59733E+06, 0.63000E+06, 0.66423E+06, 0.70010E+06,
     + 0.73767E+06, 0.77701E+06, 0.81818E+06, 0.86127E+06, 0.90635E+06,
     + 0.95349E+06, 0.10028E+07, 0.10543E+07, 0.11082E+07, 0.11644E+07,
     + 0.12232E+07, 0.12845E+07, 0.13485E+07, 0.14154E+07, 0.14851E+07,
     + 0.15578E+07, 0.16337E+07, 0.17127E+07, 0.17952E+07, 0.18810E+07,
     + 0.19705E+07, 0.20637E+07, 0.21607E+07, 0.22617E+07, 0.23669E+07,
     + 0.24763E+07, 0.25901E+07, 0.27085E+07, 0.28316E+07, 0.29596E+07,
     + 0.30926E+07/
c...        --       136
      data (QofT( 2,J),J=1,119)/ 0.53173E+03, 0.89447E+03, 0.13153E+04,
     + 0.17871E+04, 0.23051E+04, 0.28658E+04, 0.34669E+04, 0.41073E+04,
     + 0.47872E+04, 0.55074E+04, 0.62702E+04, 0.70785E+04, 0.79357E+04,
     + 0.88462E+04, 0.98147E+04, 0.10846E+05, 0.11946E+05, 0.13121E+05,
     + 0.14377E+05, 0.15721E+05, 0.17159E+05, 0.18699E+05, 0.20349E+05,
     + 0.22118E+05, 0.24013E+05, 0.26045E+05, 0.28222E+05, 0.30555E+05,
     + 0.33055E+05, 0.35733E+05, 0.38601E+05, 0.41671E+05, 0.44958E+05,
     + 0.48474E+05, 0.52235E+05, 0.56255E+05, 0.60552E+05, 0.65142E+05,
     + 0.70043E+05, 0.75275E+05, 0.80856E+05, 0.86808E+05, 0.93152E+05,
     + 0.99913E+05, 0.10711E+06, 0.11478E+06, 0.12293E+06, 0.13161E+06,
     + 0.14083E+06, 0.15063E+06, 0.16104E+06, 0.17209E+06, 0.18382E+06,
     + 0.19626E+06, 0.20945E+06, 0.22343E+06, 0.23825E+06, 0.25394E+06,
     + 0.27054E+06, 0.28812E+06, 0.30671E+06, 0.32636E+06, 0.34713E+06,
     + 0.36907E+06, 0.39224E+06, 0.41671E+06, 0.44252E+06, 0.46975E+06,
     + 0.49845E+06, 0.52872E+06, 0.56060E+06, 0.59418E+06, 0.62954E+06,
     + 0.66676E+06, 0.70591E+06, 0.74710E+06, 0.79040E+06, 0.83591E+06,
     + 0.88373E+06, 0.93395E+06, 0.98669E+06, 0.10421E+07, 0.11001E+07,
     + 0.11611E+07, 0.12250E+07, 0.12920E+07, 0.13622E+07, 0.14357E+07,
     + 0.15128E+07, 0.15934E+07, 0.16779E+07, 0.17662E+07, 0.18587E+07,
     + 0.19554E+07, 0.20565E+07, 0.21621E+07, 0.22725E+07, 0.23879E+07,
     + 0.25084E+07, 0.26342E+07, 0.27655E+07, 0.29026E+07, 0.30456E+07,
     + 0.31947E+07, 0.33502E+07, 0.35124E+07, 0.36814E+07, 0.38575E+07,
     + 0.40410E+07, 0.42321E+07, 0.44311E+07, 0.46382E+07, 0.48538E+07,
     + 0.50782E+07, 0.53116E+07, 0.55544E+07, 0.58068E+07, 0.60693E+07,
     + 0.63421E+07/
c...        --       128
      data (QofT( 3,J),J=1,119)/ 0.27198E+03, 0.45755E+03, 0.67282E+03,
     + 0.91421E+03, 0.11792E+04, 0.14660E+04, 0.17735E+04, 0.21012E+04,
     + 0.24490E+04, 0.28175E+04, 0.32077E+04, 0.36212E+04, 0.40598E+04,
     + 0.45256E+04, 0.50211E+04, 0.55488E+04, 0.61116E+04, 0.67127E+04,
     + 0.73552E+04, 0.80426E+04, 0.87783E+04, 0.95663E+04, 0.10410E+05,
     + 0.11315E+05, 0.12285E+05, 0.13324E+05, 0.14438E+05, 0.15632E+05,
     + 0.16911E+05, 0.18281E+05, 0.19748E+05, 0.21319E+05, 0.23000E+05,
     + 0.24799E+05, 0.26723E+05, 0.28780E+05, 0.30978E+05, 0.33326E+05,
     + 0.35834E+05, 0.38510E+05, 0.41365E+05, 0.44410E+05, 0.47656E+05,
     + 0.51115E+05, 0.54798E+05, 0.58719E+05, 0.62891E+05, 0.67329E+05,
     + 0.72047E+05, 0.77060E+05, 0.82385E+05, 0.88039E+05, 0.94039E+05,
     + 0.10040E+06, 0.10715E+06, 0.11431E+06, 0.12189E+06, 0.12991E+06,
     + 0.13841E+06, 0.14740E+06, 0.15691E+06, 0.16696E+06, 0.17759E+06,
     + 0.18882E+06, 0.20067E+06, 0.21318E+06, 0.22639E+06, 0.24032E+06,
     + 0.25501E+06, 0.27049E+06, 0.28680E+06, 0.30398E+06, 0.32207E+06,
     + 0.34111E+06, 0.36114E+06, 0.38221E+06, 0.40436E+06, 0.42765E+06,
     + 0.45211E+06, 0.47781E+06, 0.50479E+06, 0.53311E+06, 0.56283E+06,
     + 0.59400E+06, 0.62669E+06, 0.66097E+06, 0.69688E+06, 0.73451E+06,
     + 0.77393E+06, 0.81520E+06, 0.85840E+06, 0.90360E+06, 0.95090E+06,
     + 0.10004E+07, 0.10521E+07, 0.11061E+07, 0.11626E+07, 0.12216E+07,
     + 0.12833E+07, 0.13476E+07, 0.14148E+07, 0.14849E+07, 0.15581E+07,
     + 0.16344E+07, 0.17140E+07, 0.17969E+07, 0.18834E+07, 0.19735E+07,
     + 0.20674E+07, 0.21651E+07, 0.22669E+07, 0.23729E+07, 0.24832E+07,
     + 0.25980E+07, 0.27174E+07, 0.28416E+07, 0.29708E+07, 0.31050E+07,
     + 0.32446E+07/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_HOCl  (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 2), QofT( 2,119),Q(NT)
      data xgj/ 8.,8./
c...     HOCl
c...        --       165
      data (QofT( 1,J),J=1,119)/ 0.17041E+04, 0.28708E+04, 0.42250E+04,
     + 0.57456E+04, 0.74211E+04, 0.92470E+04, 0.11225E+05, 0.13359E+05,
     + 0.15657E+05, 0.18129E+05, 0.20785E+05, 0.23637E+05, 0.26696E+05,
     + 0.29974E+05, 0.33484E+05, 0.37239E+05, 0.41252E+05, 0.45536E+05,
     + 0.50105E+05, 0.54973E+05, 0.60152E+05, 0.65659E+05, 0.71507E+05,
     + 0.77711E+05, 0.84286E+05, 0.91249E+05, 0.98614E+05, 0.10640E+06,
     + 0.11462E+06, 0.12330E+06, 0.13244E+06, 0.14208E+06, 0.15222E+06,
     + 0.16289E+06, 0.17411E+06, 0.18589E+06, 0.19825E+06, 0.21123E+06,
     + 0.22483E+06, 0.23908E+06, 0.25400E+06, 0.26962E+06, 0.28596E+06,
     + 0.30303E+06, 0.32087E+06, 0.33950E+06, 0.35895E+06, 0.37923E+06,
     + 0.40038E+06, 0.42243E+06, 0.44539E+06, 0.46930E+06, 0.49419E+06,
     + 0.52008E+06, 0.54700E+06, 0.57498E+06, 0.60406E+06, 0.63426E+06,
     + 0.66562E+06, 0.69816E+06, 0.73192E+06, 0.76692E+06, 0.80322E+06,
     + 0.84083E+06, 0.87979E+06, 0.92014E+06, 0.96192E+06, 0.10052E+07,
     + 0.10499E+07, 0.10961E+07, 0.11440E+07, 0.11934E+07, 0.12445E+07,
     + 0.12973E+07, 0.13518E+07, 0.14081E+07, 0.14661E+07, 0.15261E+07,
     + 0.15879E+07, 0.16516E+07, 0.17174E+07, 0.17851E+07, 0.18550E+07,
     + 0.19269E+07, 0.20010E+07, 0.20773E+07, 0.21559E+07, 0.22367E+07,
     + 0.23200E+07, 0.24056E+07, 0.24936E+07, 0.25842E+07, 0.26773E+07,
     + 0.27730E+07, 0.28714E+07, 0.29724E+07, 0.30763E+07, 0.31829E+07,
     + 0.32924E+07, 0.34049E+07, 0.35203E+07, 0.36387E+07, 0.37603E+07,
     + 0.38850E+07, 0.40129E+07, 0.41441E+07, 0.42786E+07, 0.44165E+07,
     + 0.45579E+07, 0.47028E+07, 0.48512E+07, 0.50033E+07, 0.51592E+07,
     + 0.53187E+07, 0.54822E+07, 0.56495E+07, 0.58208E+07, 0.59961E+07,
     + 0.61755E+07/
c...        --       167
      data (QofT( 2,J),J=1,119)/ 0.17342E+04, 0.29215E+04, 0.42998E+04,
     + 0.58473E+04, 0.75524E+04, 0.94107E+04, 0.11423E+05, 0.13595E+05,
     + 0.15935E+05, 0.18450E+05, 0.21154E+05, 0.24056E+05, 0.27168E+05,
     + 0.30505E+05, 0.34077E+05, 0.37899E+05, 0.41983E+05, 0.46343E+05,
     + 0.50993E+05, 0.55947E+05, 0.61218E+05, 0.66822E+05, 0.72774E+05,
     + 0.79088E+05, 0.85780E+05, 0.92866E+05, 0.10036E+06, 0.10829E+06,
     + 0.11665E+06, 0.12548E+06, 0.13479E+06, 0.14460E+06, 0.15492E+06,
     + 0.16578E+06, 0.17719E+06, 0.18918E+06, 0.20177E+06, 0.21497E+06,
     + 0.22881E+06, 0.24332E+06, 0.25851E+06, 0.27440E+06, 0.29102E+06,
     + 0.30840E+06, 0.32656E+06, 0.34552E+06, 0.36531E+06, 0.38595E+06,
     + 0.40748E+06, 0.42991E+06, 0.45328E+06, 0.47762E+06, 0.50295E+06,
     + 0.52929E+06, 0.55669E+06, 0.58517E+06, 0.61477E+06, 0.64550E+06,
     + 0.67741E+06, 0.71053E+06, 0.74489E+06, 0.78052E+06, 0.81745E+06,
     + 0.85573E+06, 0.89539E+06, 0.93645E+06, 0.97897E+06, 0.10230E+07,
     + 0.10685E+07, 0.11156E+07, 0.11643E+07, 0.12146E+07, 0.12666E+07,
     + 0.13203E+07, 0.13757E+07, 0.14330E+07, 0.14921E+07, 0.15531E+07,
     + 0.16160E+07, 0.16809E+07, 0.17478E+07, 0.18168E+07, 0.18878E+07,
     + 0.19611E+07, 0.20365E+07, 0.21141E+07, 0.21941E+07, 0.22764E+07,
     + 0.23611E+07, 0.24482E+07, 0.25378E+07, 0.26300E+07, 0.27248E+07,
     + 0.28222E+07, 0.29223E+07, 0.30251E+07, 0.31308E+07, 0.32393E+07,
     + 0.33508E+07, 0.34652E+07, 0.35827E+07, 0.37032E+07, 0.38269E+07,
     + 0.39539E+07, 0.40840E+07, 0.42176E+07, 0.43545E+07, 0.44948E+07,
     + 0.46387E+07, 0.47861E+07, 0.49372E+07, 0.50920E+07, 0.52506E+07,
     + 0.54130E+07, 0.55793E+07, 0.57496E+07, 0.59239E+07, 0.61024E+07,
     + 0.62850E+07/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_N2    (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 1./
c...       N2
c...        --        44
      data (QofT( 1,J),J=1,119)/ 0.95487E+02, 0.13466E+03, 0.17386E+03,
     + 0.21307E+03, 0.25230E+03, 0.29154E+03, 0.33080E+03, 0.37008E+03,
     + 0.40937E+03, 0.44868E+03, 0.48800E+03, 0.52736E+03, 0.56674E+03,
     + 0.60616E+03, 0.64562E+03, 0.68515E+03, 0.72475E+03, 0.76445E+03,
     + 0.80426E+03, 0.84420E+03, 0.88430E+03, 0.92457E+03, 0.96505E+03,
     + 0.10057E+04, 0.10467E+04, 0.10879E+04, 0.11293E+04, 0.11711E+04,
     + 0.12132E+04, 0.12556E+04, 0.12984E+04, 0.13416E+04, 0.13851E+04,
     + 0.14291E+04, 0.14734E+04, 0.15182E+04, 0.15635E+04, 0.16091E+04,
     + 0.16553E+04, 0.17019E+04, 0.17490E+04, 0.17965E+04, 0.18446E+04,
     + 0.18932E+04, 0.19422E+04, 0.19918E+04, 0.20419E+04, 0.20926E+04,
     + 0.21437E+04, 0.21954E+04, 0.22477E+04, 0.23004E+04, 0.23538E+04,
     + 0.24077E+04, 0.24621E+04, 0.25171E+04, 0.25727E+04, 0.26288E+04,
     + 0.26856E+04, 0.27428E+04, 0.28007E+04, 0.28591E+04, 0.29181E+04,
     + 0.29777E+04, 0.30379E+04, 0.30986E+04, 0.31600E+04, 0.32219E+04,
     + 0.32844E+04, 0.33475E+04, 0.34112E+04, 0.34755E+04, 0.35404E+04,
     + 0.36059E+04, 0.36720E+04, 0.37387E+04, 0.38060E+04, 0.38739E+04,
     + 0.39424E+04, 0.40115E+04, 0.40812E+04, 0.41515E+04, 0.42224E+04,
     + 0.42939E+04, 0.43661E+04, 0.44388E+04, 0.45122E+04, 0.45861E+04,
     + 0.46607E+04, 0.47359E+04, 0.48117E+04, 0.48882E+04, 0.49652E+04,
     + 0.50428E+04, 0.51211E+04, 0.52000E+04, 0.52795E+04, 0.53596E+04,
     + 0.54404E+04, 0.55217E+04, 0.56037E+04, 0.56863E+04, 0.57695E+04,
     + 0.58533E+04, 0.59378E+04, 0.60229E+04, 0.61086E+04, 0.61950E+04,
     + 0.62819E+04, 0.63695E+04, 0.64577E+04, 0.65465E+04, 0.66360E+04,
     + 0.67261E+04, 0.68168E+04, 0.69081E+04, 0.70001E+04, 0.70927E+04,
     + 0.71859E+04/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_HCN   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 3), QofT( 3,119),Q(NT)
      data xgj/ 6.,12.,4./
c...      HCN
c...        --       124
      data (QofT( 1,J),J=1,119)/ 0.17143E+03, 0.24209E+03, 0.31285E+03,
     + 0.38392E+03, 0.45582E+03, 0.52929E+03, 0.60515E+03, 0.68424E+03,
     + 0.76731E+03, 0.85505E+03, 0.94805E+03, 0.10468E+04, 0.11519E+04,
     + 0.12637E+04, 0.13826E+04, 0.15090E+04, 0.16435E+04, 0.17863E+04,
     + 0.19378E+04, 0.20985E+04, 0.22689E+04, 0.24492E+04, 0.26401E+04,
     + 0.28418E+04, 0.30550E+04, 0.32801E+04, 0.35176E+04, 0.37680E+04,
     + 0.40318E+04, 0.43097E+04, 0.46021E+04, 0.49097E+04, 0.52330E+04,
     + 0.55727E+04, 0.59294E+04, 0.63038E+04, 0.66964E+04, 0.71081E+04,
     + 0.75396E+04, 0.79915E+04, 0.84646E+04, 0.89596E+04, 0.94774E+04,
     + 0.10019E+05, 0.10585E+05, 0.11176E+05, 0.11793E+05, 0.12437E+05,
     + 0.13108E+05, 0.13809E+05, 0.14540E+05, 0.15301E+05, 0.16094E+05,
     + 0.16919E+05, 0.17779E+05, 0.18673E+05, 0.19603E+05, 0.20570E+05,
     + 0.21575E+05, 0.22619E+05, 0.23704E+05, 0.24831E+05, 0.26000E+05,
     + 0.27213E+05, 0.28472E+05, 0.29778E+05, 0.31131E+05, 0.32534E+05,
     + 0.33987E+05, 0.35493E+05, 0.37052E+05, 0.38666E+05, 0.40336E+05,
     + 0.42064E+05, 0.43852E+05, 0.45701E+05, 0.47612E+05, 0.49587E+05,
     + 0.51629E+05, 0.53738E+05, 0.55916E+05, 0.58165E+05, 0.60486E+05,
     + 0.62883E+05, 0.65355E+05, 0.67905E+05, 0.70536E+05, 0.73249E+05,
     + 0.76045E+05, 0.78927E+05, 0.81897E+05, 0.84957E+05, 0.88108E+05,
     + 0.91354E+05, 0.94696E+05, 0.98136E+05, 0.10168E+06, 0.10532E+06,
     + 0.10907E+06, 0.11292E+06, 0.11689E+06, 0.12096E+06, 0.12516E+06,
     + 0.12946E+06, 0.13389E+06, 0.13844E+06, 0.14311E+06, 0.14791E+06,
     + 0.15284E+06, 0.15790E+06, 0.16310E+06, 0.16843E+06, 0.17391E+06,
     + 0.17953E+06, 0.18529E+06, 0.19120E+06, 0.19726E+06, 0.20348E+06,
     + 0.20986E+06/
c...        --       134
      data (QofT( 2,J),J=1,119)/ 0.35186E+03, 0.49693E+03, 0.64221E+03,
     + 0.78815E+03, 0.93585E+03, 0.10868E+04, 0.12428E+04, 0.14056E+04,
     + 0.15766E+04, 0.17574E+04, 0.19491E+04, 0.21528E+04, 0.23695E+04,
     + 0.26002E+04, 0.28457E+04, 0.31068E+04, 0.33845E+04, 0.36795E+04,
     + 0.39926E+04, 0.43249E+04, 0.46770E+04, 0.50500E+04, 0.54447E+04,
     + 0.58621E+04, 0.63032E+04, 0.67690E+04, 0.72606E+04, 0.77789E+04,
     + 0.83252E+04, 0.89005E+04, 0.95062E+04, 0.10143E+05, 0.10813E+05,
     + 0.11517E+05, 0.12256E+05, 0.13032E+05, 0.13846E+05, 0.14699E+05,
     + 0.15593E+05, 0.16530E+05, 0.17511E+05, 0.18538E+05, 0.19612E+05,
     + 0.20734E+05, 0.21908E+05, 0.23134E+05, 0.24414E+05, 0.25750E+05,
     + 0.27145E+05, 0.28599E+05, 0.30115E+05, 0.31694E+05, 0.33340E+05,
     + 0.35054E+05, 0.36838E+05, 0.38694E+05, 0.40625E+05, 0.42633E+05,
     + 0.44720E+05, 0.46889E+05, 0.49142E+05, 0.51481E+05, 0.53910E+05,
     + 0.56430E+05, 0.59045E+05, 0.61757E+05, 0.64568E+05, 0.67482E+05,
     + 0.70502E+05, 0.73630E+05, 0.76869E+05, 0.80223E+05, 0.83694E+05,
     + 0.87285E+05, 0.91000E+05, 0.94843E+05, 0.98815E+05, 0.10292E+06,
     + 0.10716E+06, 0.11155E+06, 0.11608E+06, 0.12075E+06, 0.12558E+06,
     + 0.13056E+06, 0.13570E+06, 0.14100E+06, 0.14647E+06, 0.15211E+06,
     + 0.15793E+06, 0.16392E+06, 0.17009E+06, 0.17646E+06, 0.18301E+06,
     + 0.18976E+06, 0.19671E+06, 0.20387E+06, 0.21123E+06, 0.21881E+06,
     + 0.22660E+06, 0.23462E+06, 0.24287E+06, 0.25135E+06, 0.26007E+06,
     + 0.26903E+06, 0.27824E+06, 0.28771E+06, 0.29743E+06, 0.30742E+06,
     + 0.31767E+06, 0.32820E+06, 0.33901E+06, 0.35011E+06, 0.36150E+06,
     + 0.37319E+06, 0.38518E+06, 0.39749E+06, 0.41010E+06, 0.42304E+06,
     + 0.43631E+06/
c...        --       125
      data (QofT( 3,J),J=1,119)/ 0.11863E+03, 0.16755E+03, 0.21653E+03,
     + 0.26576E+03, 0.31559E+03, 0.36656E+03, 0.41926E+03, 0.47428E+03,
     + 0.53214E+03, 0.59333E+03, 0.65824E+03, 0.72727E+03, 0.80074E+03,
     + 0.87898E+03, 0.96227E+03, 0.10509E+04, 0.11452E+04, 0.12454E+04,
     + 0.13518E+04, 0.14647E+04, 0.15844E+04, 0.17112E+04, 0.18455E+04,
     + 0.19875E+04, 0.21377E+04, 0.22962E+04, 0.24636E+04, 0.26402E+04,
     + 0.28263E+04, 0.30224E+04, 0.32289E+04, 0.34461E+04, 0.36745E+04,
     + 0.39145E+04, 0.41667E+04, 0.44314E+04, 0.47092E+04, 0.50005E+04,
     + 0.53059E+04, 0.56259E+04, 0.59609E+04, 0.63116E+04, 0.66785E+04,
     + 0.70622E+04, 0.74633E+04, 0.78823E+04, 0.83200E+04, 0.87769E+04,
     + 0.92536E+04, 0.97509E+04, 0.10269E+05, 0.10810E+05, 0.11373E+05,
     + 0.11959E+05, 0.12570E+05, 0.13205E+05, 0.13866E+05, 0.14554E+05,
     + 0.15268E+05, 0.16011E+05, 0.16782E+05, 0.17583E+05, 0.18415E+05,
     + 0.19279E+05, 0.20174E+05, 0.21103E+05, 0.22067E+05, 0.23065E+05,
     + 0.24100E+05, 0.25172E+05, 0.26282E+05, 0.27432E+05, 0.28622E+05,
     + 0.29853E+05, 0.31127E+05, 0.32445E+05, 0.33807E+05, 0.35215E+05,
     + 0.36670E+05, 0.38174E+05, 0.39727E+05, 0.41330E+05, 0.42986E+05,
     + 0.44695E+05, 0.46459E+05, 0.48278E+05, 0.50155E+05, 0.52091E+05,
     + 0.54086E+05, 0.56143E+05, 0.58263E+05, 0.60447E+05, 0.62696E+05,
     + 0.65013E+05, 0.67399E+05, 0.69856E+05, 0.72384E+05, 0.74986E+05,
     + 0.77663E+05, 0.80416E+05, 0.83249E+05, 0.86161E+05, 0.89156E+05,
     + 0.92233E+05, 0.95397E+05, 0.98648E+05, 0.10199E+06, 0.10542E+06,
     + 0.10894E+06, 0.11256E+06, 0.11627E+06, 0.12009E+06, 0.12400E+06,
     + 0.12802E+06, 0.13214E+06, 0.13636E+06, 0.14070E+06, 0.14515E+06,
     + 0.14971E+06/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_CH3Cl (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 2), QofT( 2,119),Q(NT)
      data xgj/ 4.,4./
c...    CH3Cl
c...        --       215
      data (QofT( 1,J),J=1,119)/ 0.10106E+05, 0.17025E+05, 0.25055E+05,
     + 0.34073E+05, 0.44011E+05, 0.54858E+05, 0.66651E+05, 0.79468E+05,
     + 0.93426E+05, 0.10867E+06, 0.12538E+06, 0.14375E+06, 0.16401E+06,
     + 0.18641E+06, 0.21121E+06, 0.23871E+06, 0.26925E+06, 0.30317E+06,
     + 0.34086E+06, 0.38274E+06, 0.42928E+06, 0.48098E+06, 0.53840E+06,
     + 0.60213E+06, 0.67284E+06, 0.75125E+06, 0.83814E+06, 0.93438E+06,
     + 0.10409E+07, 0.11587E+07, 0.12890E+07, 0.14328E+07, 0.15916E+07,
     + 0.17668E+07, 0.19599E+07, 0.21727E+07, 0.24069E+07, 0.26645E+07,
     + 0.29478E+07, 0.32590E+07, 0.36006E+07, 0.39754E+07, 0.43864E+07,
     + 0.48367E+07, 0.53297E+07, 0.58692E+07, 0.64591E+07, 0.71038E+07,
     + 0.78078E+07, 0.85762E+07, 0.94143E+07, 0.10328E+08, 0.11323E+08,
     + 0.12406E+08, 0.13585E+08, 0.14867E+08, 0.16260E+08, 0.17772E+08,
     + 0.19414E+08, 0.21195E+08, 0.23126E+08, 0.25218E+08, 0.27484E+08,
     + 0.29936E+08, 0.32589E+08, 0.35456E+08, 0.38555E+08, 0.41901E+08,
     + 0.45513E+08, 0.49409E+08, 0.53610E+08, 0.58137E+08, 0.63014E+08,
     + 0.68264E+08, 0.73913E+08, 0.79989E+08, 0.86521E+08, 0.93539E+08,
     + 0.10108E+09, 0.10917E+09, 0.11785E+09, 0.12716E+09, 0.13714E+09,
     + 0.14783E+09, 0.15928E+09, 0.17154E+09, 0.18466E+09, 0.19869E+09,
     + 0.21369E+09, 0.22972E+09, 0.24685E+09, 0.26513E+09, 0.28465E+09,
     + 0.30547E+09, 0.32768E+09, 0.35136E+09, 0.37659E+09, 0.40346E+09,
     + 0.43208E+09, 0.46254E+09, 0.49495E+09, 0.52943E+09, 0.56608E+09,
     + 0.60504E+09, 0.64643E+09, 0.69039E+09, 0.73707E+09, 0.78660E+09,
     + 0.83916E+09, 0.89491E+09, 0.95401E+09, 0.10167E+10, 0.10830E+10,
     + 0.11533E+10, 0.12278E+10, 0.13066E+10, 0.13900E+10, 0.14782E+10,
     + 0.15715E+10/
c...        --       217
      data (QofT( 2,J),J=1,119)/ 0.10265E+05, 0.17294E+05, 0.25452E+05,
     + 0.34612E+05, 0.44707E+05, 0.55726E+05, 0.67706E+05, 0.80727E+05,
     + 0.94906E+05, 0.11039E+06, 0.12737E+06, 0.14603E+06, 0.16661E+06,
     + 0.18936E+06, 0.21456E+06, 0.24250E+06, 0.27352E+06, 0.30798E+06,
     + 0.34626E+06, 0.38881E+06, 0.43609E+06, 0.48861E+06, 0.54694E+06,
     + 0.61168E+06, 0.68351E+06, 0.76316E+06, 0.85144E+06, 0.94920E+06,
     + 0.10574E+07, 0.11771E+07, 0.13094E+07, 0.14556E+07, 0.16169E+07,
     + 0.17949E+07, 0.19910E+07, 0.22072E+07, 0.24451E+07, 0.27068E+07,
     + 0.29946E+07, 0.33107E+07, 0.36578E+07, 0.40385E+07, 0.44560E+07,
     + 0.49134E+07, 0.54143E+07, 0.59624E+07, 0.65617E+07, 0.72166E+07,
     + 0.79318E+07, 0.87124E+07, 0.95638E+07, 0.10492E+08, 0.11503E+08,
     + 0.12603E+08, 0.13801E+08, 0.15103E+08, 0.16518E+08, 0.18055E+08,
     + 0.19723E+08, 0.21532E+08, 0.23494E+08, 0.25619E+08, 0.27921E+08,
     + 0.30412E+08, 0.33106E+08, 0.36020E+08, 0.39167E+08, 0.42567E+08,
     + 0.46236E+08, 0.50194E+08, 0.54462E+08, 0.59061E+08, 0.64015E+08,
     + 0.69349E+08, 0.75088E+08, 0.81261E+08, 0.87896E+08, 0.95026E+08,
     + 0.10268E+09, 0.11090E+09, 0.11972E+09, 0.12918E+09, 0.13932E+09,
     + 0.15018E+09, 0.16181E+09, 0.17427E+09, 0.18759E+09, 0.20185E+09,
     + 0.21709E+09, 0.23337E+09, 0.25077E+09, 0.26935E+09, 0.28918E+09,
     + 0.31033E+09, 0.33289E+09, 0.35695E+09, 0.38258E+09, 0.40988E+09,
     + 0.43895E+09, 0.46990E+09, 0.50283E+09, 0.53785E+09, 0.57509E+09,
     + 0.61467E+09, 0.65672E+09, 0.70138E+09, 0.74880E+09, 0.79912E+09,
     + 0.85252E+09, 0.90915E+09, 0.96919E+09, 0.10328E+10, 0.11003E+10,
     + 0.11717E+10, 0.12473E+10, 0.13274E+10, 0.14121E+10, 0.15017E+10,
     + 0.15965E+10/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_H2O2  (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 1./
c...     H2O2
c...        --      1661
      data (QofT( 1,J),J=1,119)/ 0.62392E+03, 0.10958E+04, 0.16692E+04,
     + 0.23492E+04, 0.31427E+04, 0.40574E+04, 0.51014E+04, 0.62840E+04,
     + 0.76157E+04, 0.91085E+04, 0.10776E+05, 0.12633E+05, 0.14696E+05,
     + 0.16983E+05, 0.19515E+05, 0.22312E+05, 0.25396E+05, 0.28792E+05,
     + 0.32526E+05, 0.36625E+05, 0.41118E+05, 0.46036E+05, 0.51410E+05,
     + 0.57275E+05, 0.63667E+05, 0.70623E+05, 0.78185E+05, 0.86394E+05,
     + 0.95295E+05, 0.10493E+06, 0.11536E+06, 0.12662E+06, 0.13878E+06,
     + 0.15188E+06, 0.16600E+06, 0.18118E+06, 0.19750E+06, 0.21503E+06,
     + 0.23383E+06, 0.25398E+06, 0.27556E+06, 0.29864E+06, 0.32333E+06,
     + 0.34970E+06, 0.37784E+06, 0.40786E+06, 0.43985E+06, 0.47392E+06,
     + 0.51018E+06, 0.54874E+06, 0.58972E+06, 0.63324E+06, 0.67943E+06,
     + 0.72843E+06, 0.78037E+06, 0.83540E+06, 0.89366E+06, 0.95530E+06,
     + 0.10205E+07, 0.10894E+07, 0.11622E+07, 0.12391E+07, 0.13202E+07,
     + 0.14057E+07, 0.14959E+07, 0.15909E+07, 0.16910E+07, 0.17963E+07,
     + 0.19072E+07, 0.20237E+07, 0.21463E+07, 0.22750E+07, 0.24102E+07,
     + 0.25522E+07, 0.27012E+07, 0.28575E+07, 0.30213E+07, 0.31931E+07,
     + 0.33730E+07, 0.35615E+07, 0.37588E+07, 0.39653E+07, 0.41813E+07,
     + 0.44072E+07, 0.46433E+07, 0.48901E+07, 0.51479E+07, 0.54171E+07,
     + 0.56982E+07, 0.59915E+07, 0.62976E+07, 0.66167E+07, 0.69495E+07,
     + 0.72963E+07, 0.76577E+07, 0.80342E+07, 0.84262E+07, 0.88343E+07,
     + 0.92591E+07, 0.97011E+07, 0.10161E+08, 0.10639E+08, 0.11136E+08,
     + 0.11652E+08, 0.12189E+08, 0.12746E+08, 0.13325E+08, 0.13926E+08,
     + 0.14550E+08, 0.15198E+08, 0.15870E+08, 0.16566E+08, 0.17289E+08,
     + 0.18038E+08, 0.18814E+08, 0.19619E+08, 0.20452E+08, 0.21315E+08,
     + 0.22209E+08/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_C2H2  (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 2), QofT( 2,119),Q(NT)
      data xgj/ 1.,8./
c...     C2H2
c...        --      1221
      data (QofT( 1,J),J=1,119)/ 0.71617E+02, 0.10121E+03, 0.13092E+03,
     + 0.16104E+03, 0.19218E+03, 0.22509E+03, 0.26062E+03, 0.29959E+03,
     + 0.34281E+03, 0.39103E+03, 0.44503E+03, 0.50558E+03, 0.57346E+03,
     + 0.64950E+03, 0.73457E+03, 0.82960E+03, 0.93557E+03, 0.10535E+04,
     + 0.11846E+04, 0.13301E+04, 0.14911E+04, 0.16692E+04, 0.18658E+04,
     + 0.20825E+04, 0.23211E+04, 0.25833E+04, 0.28711E+04, 0.31867E+04,
     + 0.35323E+04, 0.39102E+04, 0.43230E+04, 0.47735E+04, 0.52645E+04,
     + 0.57991E+04, 0.63807E+04, 0.70127E+04, 0.76988E+04, 0.84430E+04,
     + 0.92495E+04, 0.10123E+05, 0.11067E+05, 0.12088E+05, 0.13191E+05,
     + 0.14381E+05, 0.15664E+05, 0.17047E+05, 0.18536E+05, 0.20137E+05,
     + 0.21859E+05, 0.23710E+05, 0.25696E+05, 0.27827E+05, 0.30112E+05,
     + 0.32561E+05, 0.35183E+05, 0.37990E+05, 0.40991E+05, 0.44199E+05,
     + 0.47626E+05, 0.51285E+05, 0.55189E+05, 0.59353E+05, 0.63791E+05,
     + 0.68518E+05, 0.73551E+05, 0.78908E+05, 0.84604E+05, 0.90661E+05,
     + 0.97095E+05, 0.10393E+06, 0.11118E+06, 0.11888E+06, 0.12704E+06,
     + 0.13569E+06, 0.14486E+06, 0.15457E+06, 0.16485E+06, 0.17572E+06,
     + 0.18722E+06, 0.19938E+06, 0.21223E+06, 0.22581E+06, 0.24014E+06,
     + 0.25527E+06, 0.27123E+06, 0.28807E+06, 0.30582E+06, 0.32452E+06,
     + 0.34423E+06, 0.36498E+06, 0.38683E+06, 0.40982E+06, 0.43401E+06,
     + 0.45944E+06, 0.48618E+06, 0.51428E+06, 0.54380E+06, 0.57480E+06,
     + 0.60735E+06, 0.64151E+06, 0.67735E+06, 0.71495E+06, 0.75436E+06,
     + 0.79568E+06, 0.83898E+06, 0.88434E+06, 0.93184E+06, 0.98158E+06,
     + 0.10336E+07, 0.10881E+07, 0.11451E+07, 0.12047E+07, 0.12670E+07,
     + 0.13321E+07, 0.14002E+07, 0.14713E+07, 0.15455E+07, 0.16231E+07,
     + 0.17040E+07/
c...        --      1231
      data (QofT( 2,J),J=1,119)/ 0.28647E+03, 0.40486E+03, 0.52369E+03,
     + 0.64419E+03, 0.76874E+03, 0.90040E+03, 0.10425E+04, 0.11984E+04,
     + 0.13713E+04, 0.15642E+04, 0.17802E+04, 0.20223E+04, 0.22939E+04,
     + 0.25981E+04, 0.29384E+04, 0.33185E+04, 0.37424E+04, 0.42142E+04,
     + 0.47386E+04, 0.53203E+04, 0.59646E+04, 0.66769E+04, 0.74634E+04,
     + 0.83302E+04, 0.92845E+04, 0.10333E+05, 0.11485E+05, 0.12747E+05,
     + 0.14129E+05, 0.15641E+05, 0.17292E+05, 0.19094E+05, 0.21058E+05,
     + 0.23197E+05, 0.25523E+05, 0.28051E+05, 0.30796E+05, 0.33773E+05,
     + 0.36999E+05, 0.40492E+05, 0.44270E+05, 0.48354E+05, 0.52765E+05,
     + 0.57525E+05, 0.62658E+05, 0.68189E+05, 0.74144E+05, 0.80551E+05,
     + 0.87439E+05, 0.94840E+05, 0.10279E+06, 0.11131E+06, 0.12045E+06,
     + 0.13025E+06, 0.14074E+06, 0.15196E+06, 0.16397E+06, 0.17680E+06,
     + 0.19051E+06, 0.20514E+06, 0.22076E+06, 0.23742E+06, 0.25517E+06,
     + 0.27408E+06, 0.29421E+06, 0.31564E+06, 0.33842E+06, 0.36265E+06,
     + 0.38839E+06, 0.41572E+06, 0.44474E+06, 0.47553E+06, 0.50818E+06,
     + 0.54278E+06, 0.57945E+06, 0.61829E+06, 0.65940E+06, 0.70289E+06,
     + 0.74890E+06, 0.79754E+06, 0.84894E+06, 0.90324E+06, 0.96057E+06,
     + 0.10211E+07, 0.10849E+07, 0.11523E+07, 0.12233E+07, 0.12981E+07,
     + 0.13769E+07, 0.14599E+07, 0.15473E+07, 0.16393E+07, 0.17361E+07,
     + 0.18378E+07, 0.19447E+07, 0.20571E+07, 0.21752E+07, 0.22992E+07,
     + 0.24294E+07, 0.25661E+07, 0.27094E+07, 0.28598E+07, 0.30175E+07,
     + 0.31828E+07, 0.33560E+07, 0.35374E+07, 0.37274E+07, 0.39264E+07,
     + 0.41346E+07, 0.43525E+07, 0.45805E+07, 0.48188E+07, 0.50681E+07,
     + 0.53286E+07, 0.56008E+07, 0.58852E+07, 0.61823E+07, 0.64924E+07,
     + 0.68162E+07/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_C2H6  (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 1./
c...     C2H6
c...        --      1221
      data (QofT( 1,J),J=1,119)/ 0.47267E+04, 0.80011E+04, 0.11928E+05,
     + 0.16564E+05, 0.21985E+05, 0.28287E+05, 0.35590E+05, 0.44049E+05,
     + 0.53862E+05, 0.65277E+05, 0.78597E+05, 0.94191E+05, 0.11250E+06,
     + 0.13407E+06, 0.15952E+06, 0.18962E+06, 0.22526E+06, 0.26751E+06,
     + 0.31763E+06, 0.37714E+06, 0.44780E+06, 0.53174E+06, 0.63145E+06,
     + 0.74989E+06, 0.89056E+06, 0.10576E+07, 0.12559E+07, 0.14912E+07,
     + 0.17704E+07, 0.21013E+07, 0.24936E+07, 0.29582E+07, 0.35083E+07,
     + 0.41591E+07, 0.49286E+07, 0.58379E+07, 0.69116E+07, 0.81787E+07,
     + 0.96728E+07, 0.11433E+08, 0.13506E+08, 0.15945E+08, 0.18812E+08,
     + 0.22180E+08, 0.26134E+08, 0.30770E+08, 0.36204E+08, 0.42565E+08,
     + 0.50008E+08, 0.58708E+08, 0.68868E+08, 0.80725E+08, 0.94548E+08,
     + 0.11065E+09, 0.12940E+09, 0.15119E+09, 0.17652E+09, 0.20593E+09,
     + 0.24003E+09, 0.27956E+09, 0.32533E+09, 0.37829E+09, 0.43951E+09,
     + 0.51021E+09, 0.59180E+09, 0.68588E+09, 0.79427E+09, 0.91904E+09,
     + 0.10625E+10, 0.12275E+10, 0.14168E+10, 0.16341E+10, 0.18831E+10,
     + 0.21684E+10, 0.24949E+10, 0.28684E+10, 0.32951E+10, 0.37823E+10,
     + 0.43382E+10, 0.49719E+10, 0.56938E+10, 0.65156E+10, 0.74502E+10,
     + 0.85125E+10, 0.97190E+10, 0.11088E+11, 0.12641E+11, 0.14401E+11,
     + 0.16393E+11, 0.18648E+11, 0.21198E+11, 0.24079E+11, 0.27332E+11,
     + 0.31003E+11, 0.35142E+11, 0.39807E+11, 0.45060E+11, 0.50972E+11,
     + 0.57620E+11, 0.65091E+11, 0.73483E+11, 0.82902E+11, 0.93467E+11,
     + 0.10531E+12, 0.11858E+12, 0.13343E+12, 0.15005E+12, 0.16864E+12,
     + 0.18941E+12, 0.21260E+12, 0.23849E+12, 0.26737E+12, 0.29957E+12,
     + 0.33545E+12, 0.37541E+12, 0.41987E+12, 0.46934E+12, 0.52432E+12,
     + 0.58542E+12/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_PH3   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 2./
c...      PH3
c...        --      1111
      data (QofT( 1,J),J=1,119)/ 0.29652E+03, 0.49643E+03, 0.72810E+03,
     + 0.98777E+03, 0.12729E+04, 0.15820E+04, 0.19145E+04, 0.22708E+04,
     + 0.26520E+04, 0.30600E+04, 0.34971E+04, 0.39662E+04, 0.44701E+04,
     + 0.50123E+04, 0.55964E+04, 0.62261E+04, 0.69055E+04, 0.76389E+04,
     + 0.84308E+04, 0.92861E+04, 0.10210E+05, 0.11207E+05, 0.12284E+05,
     + 0.13447E+05, 0.14703E+05, 0.16057E+05, 0.17518E+05, 0.19093E+05,
     + 0.20790E+05, 0.22619E+05, 0.24588E+05, 0.26707E+05, 0.28987E+05,
     + 0.31439E+05, 0.34073E+05, 0.36903E+05, 0.39940E+05, 0.43200E+05,
     + 0.46695E+05, 0.50442E+05, 0.54455E+05, 0.58753E+05, 0.63351E+05,
     + 0.68269E+05, 0.73527E+05, 0.79144E+05, 0.85143E+05, 0.91545E+05,
     + 0.98374E+05, 0.10566E+06, 0.11342E+06, 0.12168E+06, 0.13048E+06,
     + 0.13984E+06, 0.14980E+06, 0.16038E+06, 0.17163E+06, 0.18357E+06,
     + 0.19625E+06, 0.20970E+06, 0.22396E+06, 0.23908E+06, 0.25510E+06,
     + 0.27207E+06, 0.29003E+06, 0.30904E+06, 0.32914E+06, 0.35040E+06,
     + 0.37286E+06, 0.39659E+06, 0.42164E+06, 0.44809E+06, 0.47599E+06,
     + 0.50542E+06, 0.53645E+06, 0.56914E+06, 0.60359E+06, 0.63986E+06,
     + 0.67804E+06, 0.71822E+06, 0.76049E+06, 0.80493E+06, 0.85164E+06,
     + 0.90073E+06, 0.95229E+06, 0.10064E+07, 0.10633E+07, 0.11229E+07,
     + 0.11855E+07, 0.12511E+07, 0.13199E+07, 0.13920E+07, 0.14675E+07,
     + 0.15466E+07, 0.16294E+07, 0.17161E+07, 0.18068E+07, 0.19017E+07,
     + 0.20010E+07, 0.21047E+07, 0.22132E+07, 0.23265E+07, 0.24449E+07,
     + 0.25685E+07, 0.26976E+07, 0.28323E+07, 0.29728E+07, 0.31195E+07,
     + 0.32724E+07, 0.34319E+07, 0.35981E+07, 0.37714E+07, 0.39519E+07,
     + 0.41399E+07, 0.43357E+07, 0.45395E+07, 0.47517E+07, 0.49724E+07,
     + 0.52022E+07/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_COF2  (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 1./
c...     COF2
c...        --       269
      data (QofT( 1,J),J=1,119)/ 0.54999E+04, 0.92749E+04, 0.13668E+05,
     + 0.18643E+05, 0.24224E+05, 0.30487E+05, 0.37547E+05, 0.45543E+05,
     + 0.54639E+05, 0.65019E+05, 0.76886E+05, 0.90462E+05, 0.10600E+06,
     + 0.12377E+06, 0.14407E+06, 0.16723E+06, 0.19363E+06, 0.22367E+06,
     + 0.25780E+06, 0.29650E+06, 0.34031E+06, 0.38982E+06, 0.44568E+06,
     + 0.50859E+06, 0.57932E+06, 0.65872E+06, 0.74770E+06, 0.84724E+06,
     + 0.95844E+06, 0.10825E+07, 0.12205E+07, 0.13741E+07, 0.15446E+07,
     + 0.17336E+07, 0.19428E+07, 0.21742E+07, 0.24296E+07, 0.27113E+07,
     + 0.30214E+07, 0.33626E+07, 0.37373E+07, 0.41484E+07, 0.45989E+07,
     + 0.50921E+07, 0.56313E+07, 0.62202E+07, 0.68626E+07, 0.75628E+07,
     + 0.83251E+07, 0.91542E+07, 0.10055E+08, 0.11033E+08, 0.12093E+08,
     + 0.13242E+08, 0.14486E+08, 0.15831E+08, 0.17284E+08, 0.18853E+08,
     + 0.20546E+08, 0.22371E+08, 0.24335E+08, 0.26450E+08, 0.28724E+08,
     + 0.31167E+08, 0.33790E+08, 0.36605E+08, 0.39623E+08, 0.42856E+08,
     + 0.46318E+08, 0.50022E+08, 0.53983E+08, 0.58215E+08, 0.62735E+08,
     + 0.67558E+08, 0.72702E+08, 0.78186E+08, 0.84028E+08, 0.90247E+08,
     + 0.96865E+08, 0.10390E+09, 0.11138E+09, 0.11933E+09, 0.12777E+09,
     + 0.13672E+09, 0.14622E+09, 0.15629E+09, 0.16695E+09, 0.17825E+09,
     + 0.19021E+09, 0.20287E+09, 0.21625E+09, 0.23039E+09, 0.24534E+09,
     + 0.26113E+09, 0.27779E+09, 0.29538E+09, 0.31392E+09, 0.33348E+09,
     + 0.35409E+09, 0.37580E+09, 0.39867E+09, 0.42274E+09, 0.44806E+09,
     + 0.47470E+09, 0.50271E+09, 0.53215E+09, 0.56308E+09, 0.59557E+09,
     + 0.62968E+09, 0.66548E+09, 0.70304E+09, 0.74243E+09, 0.78374E+09,
     + 0.82703E+09, 0.87240E+09, 0.91992E+09, 0.96967E+09, 0.10218E+10,
     + 0.10763E+10/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_SF6   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 1./
c...      SF6
c...        --        29
      data (QofT( 1,J),J=1,119)/ 0.46373E+05, 0.78844E+05, 0.11939E+06,
     + 0.17183E+06, 0.24247E+06, 0.34059E+06, 0.47963E+06, 0.67906E+06,
     + 0.96713E+06, 0.13848E+07, 0.19911E+07, 0.28714E+07, 0.41481E+07,
     + 0.59956E+07, 0.86617E+07, 0.12496E+08, 0.17991E+08, 0.25832E+08,
     + 0.36971E+08, 0.52724E+08, 0.74895E+08, 0.10595E+09, 0.14923E+09,
     + 0.20925E+09, 0.29208E+09, 0.40582E+09, 0.56124E+09, 0.77259E+09,
     + 0.10586E+10, 0.14439E+10, 0.19605E+10, 0.26500E+10, 0.35662E+10,
     + 0.47781E+10, 0.63747E+10, 0.84689E+10, 0.11205E+11, 0.14765E+11,
     + 0.19378E+11, 0.25336E+11, 0.32998E+11, 0.42819E+11, 0.55361E+11,
     + 0.71323E+11, 0.91569E+11, 0.11716E+12, 0.14941E+12, 0.18992E+12,
     + 0.24065E+12, 0.30398E+12, 0.38283E+12, 0.48069E+12, 0.60182E+12,
     + 0.75136E+12, 0.93546E+12, 0.11615E+13, 0.14384E+13, 0.17767E+13,
     + 0.21890E+13, 0.26903E+13, 0.32984E+13, 0.40344E+13, 0.49232E+13,
     + 0.59942E+13, 0.72819E+13, 0.88272E+13, 0.10678E+14, 0.12889E+14,
     + 0.15527E+14, 0.18666E+14, 0.22397E+14, 0.26823E+14, 0.32062E+14,
     + 0.38253E+14, 0.45558E+14, 0.54161E+14, 0.64277E+14, 0.76153E+14,
     + 0.90072E+14, 0.10636E+15, 0.12539E+15, 0.14759E+15, 0.17345E+15,
     + 0.20354E+15, 0.23848E+15, 0.27902E+15, 0.32597E+15, 0.38028E+15,
     + 0.44303E+15, 0.51542E+15, 0.59883E+15, 0.69482E+15, 0.80516E+15,
     + 0.93182E+15, 0.10770E+16, 0.12434E+16, 0.14336E+16, 0.16511E+16,
     + 0.18992E+16, 0.21821E+16, 0.25043E+16, 0.28709E+16, 0.32875E+16,
     + 0.37604E+16, 0.42968E+16, 0.49046E+16, 0.55925E+16, 0.63704E+16,
     + 0.72492E+16, 0.82411E+16, 0.93596E+16, 0.10620E+17, 0.12038E+17,
     + 0.13633E+17, 0.15425E+17, 0.17438E+17, 0.19694E+17, 0.22224E+17,
     + 0.25057E+17/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_H2S   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 3), QofT( 3,119),Q(NT)
      data xgj/ 1.,1.,4./
c...      H2S
c...        --       121
      data (QofT( 1,J),J=1,119)/ 0.47192E+02, 0.78671E+02, 0.11510E+03,
     + 0.15589E+03, 0.20061E+03, 0.24896E+03, 0.30070E+03, 0.35571E+03,
     + 0.41386E+03, 0.47513E+03, 0.53951E+03, 0.60703E+03, 0.67772E+03,
     + 0.75167E+03, 0.82896E+03, 0.90969E+03, 0.99396E+03, 0.10819E+04,
     + 0.11736E+04, 0.12692E+04, 0.13689E+04, 0.14727E+04, 0.15809E+04,
     + 0.16937E+04, 0.18111E+04, 0.19333E+04, 0.20606E+04, 0.21931E+04,
     + 0.23309E+04, 0.24744E+04, 0.26236E+04, 0.27788E+04, 0.29403E+04,
     + 0.31081E+04, 0.32825E+04, 0.34638E+04, 0.36522E+04, 0.38478E+04,
     + 0.40510E+04, 0.42619E+04, 0.44808E+04, 0.47080E+04, 0.49437E+04,
     + 0.51881E+04, 0.54415E+04, 0.57042E+04, 0.59764E+04, 0.62584E+04,
     + 0.65505E+04, 0.68529E+04, 0.71660E+04, 0.74899E+04, 0.78251E+04,
     + 0.81718E+04, 0.85303E+04, 0.89008E+04, 0.92838E+04, 0.96795E+04,
     + 0.10088E+05, 0.10510E+05, 0.10946E+05, 0.11396E+05, 0.11860E+05,
     + 0.12339E+05, 0.12833E+05, 0.13342E+05, 0.13867E+05, 0.14408E+05,
     + 0.14966E+05, 0.15540E+05, 0.16132E+05, 0.16741E+05, 0.17368E+05,
     + 0.18013E+05, 0.18677E+05, 0.19361E+05, 0.20064E+05, 0.20786E+05,
     + 0.21529E+05, 0.22293E+05, 0.23078E+05, 0.23885E+05, 0.24714E+05,
     + 0.25565E+05, 0.26439E+05, 0.27337E+05, 0.28258E+05, 0.29204E+05,
     + 0.30174E+05, 0.31170E+05, 0.32191E+05, 0.33239E+05, 0.34313E+05,
     + 0.35414E+05, 0.36543E+05, 0.37700E+05, 0.38886E+05, 0.40101E+05,
     + 0.41346E+05, 0.42621E+05, 0.43926E+05, 0.45263E+05, 0.46631E+05,
     + 0.48033E+05, 0.49466E+05, 0.50934E+05, 0.52435E+05, 0.53971E+05,
     + 0.55542E+05, 0.57149E+05, 0.58792E+05, 0.60472E+05, 0.62190E+05,
     + 0.63946E+05, 0.65740E+05, 0.67574E+05, 0.69448E+05, 0.71362E+05,
     + 0.73318E+05/
c...        --       141
      data (QofT( 2,J),J=1,119)/ 0.47310E+02, 0.78869E+02, 0.11539E+03,
     + 0.15628E+03, 0.20112E+03, 0.24959E+03, 0.30147E+03, 0.35661E+03,
     + 0.41491E+03, 0.47634E+03, 0.54088E+03, 0.60857E+03, 0.67945E+03,
     + 0.75359E+03, 0.83107E+03, 0.91201E+03, 0.99649E+03, 0.10846E+04,
     + 0.11766E+04, 0.12724E+04, 0.13724E+04, 0.14765E+04, 0.15850E+04,
     + 0.16980E+04, 0.18157E+04, 0.19382E+04, 0.20658E+04, 0.21987E+04,
     + 0.23369E+04, 0.24807E+04, 0.26303E+04, 0.27859E+04, 0.29478E+04,
     + 0.31160E+04, 0.32909E+04, 0.34727E+04, 0.36615E+04, 0.38576E+04,
     + 0.40613E+04, 0.42728E+04, 0.44923E+04, 0.47200E+04, 0.49563E+04,
     + 0.52013E+04, 0.54554E+04, 0.57188E+04, 0.59917E+04, 0.62744E+04,
     + 0.65672E+04, 0.68704E+04, 0.71843E+04, 0.75090E+04, 0.78451E+04,
     + 0.81926E+04, 0.85520E+04, 0.89236E+04, 0.93075E+04, 0.97042E+04,
     + 0.10114E+05, 0.10537E+05, 0.10974E+05, 0.11425E+05, 0.11890E+05,
     + 0.12370E+05, 0.12866E+05, 0.13376E+05, 0.13903E+05, 0.14445E+05,
     + 0.15004E+05, 0.15580E+05, 0.16173E+05, 0.16784E+05, 0.17412E+05,
     + 0.18059E+05, 0.18725E+05, 0.19410E+05, 0.20115E+05, 0.20839E+05,
     + 0.21584E+05, 0.22350E+05, 0.23137E+05, 0.23946E+05, 0.24777E+05,
     + 0.25630E+05, 0.26507E+05, 0.27407E+05, 0.28330E+05, 0.29278E+05,
     + 0.30251E+05, 0.31249E+05, 0.32273E+05, 0.33324E+05, 0.34401E+05,
     + 0.35505E+05, 0.36637E+05, 0.37797E+05, 0.38985E+05, 0.40204E+05,
     + 0.41451E+05, 0.42729E+05, 0.44038E+05, 0.45379E+05, 0.46751E+05,
     + 0.48155E+05, 0.49593E+05, 0.51064E+05, 0.52569E+05, 0.54109E+05,
     + 0.55684E+05, 0.57295E+05, 0.58943E+05, 0.60627E+05, 0.62349E+05,
     + 0.64109E+05, 0.65908E+05, 0.67747E+05, 0.69625E+05, 0.71544E+05,
     + 0.73505E+05/
c...        --       131
      data (QofT( 3,J),J=1,119)/ 0.18901E+03, 0.31509E+03, 0.46102E+03,
     + 0.62437E+03, 0.80349E+03, 0.99713E+03, 0.12044E+04, 0.14247E+04,
     + 0.16576E+04, 0.19030E+04, 0.21609E+04, 0.24313E+04, 0.27145E+04,
     + 0.30106E+04, 0.33202E+04, 0.36436E+04, 0.39811E+04, 0.43332E+04,
     + 0.47005E+04, 0.50835E+04, 0.54827E+04, 0.58987E+04, 0.63321E+04,
     + 0.67836E+04, 0.72538E+04, 0.77434E+04, 0.82532E+04, 0.87838E+04,
     + 0.93360E+04, 0.99106E+04, 0.10508E+05, 0.11130E+05, 0.11777E+05,
     + 0.12449E+05, 0.13147E+05, 0.13874E+05, 0.14628E+05, 0.15412E+05,
     + 0.16225E+05, 0.17070E+05, 0.17947E+05, 0.18857E+05, 0.19801E+05,
     + 0.20780E+05, 0.21795E+05, 0.22847E+05, 0.23937E+05, 0.25067E+05,
     + 0.26236E+05, 0.27448E+05, 0.28702E+05, 0.29999E+05, 0.31342E+05,
     + 0.32730E+05, 0.34166E+05, 0.35650E+05, 0.37184E+05, 0.38769E+05,
     + 0.40406E+05, 0.42097E+05, 0.43842E+05, 0.45644E+05, 0.47503E+05,
     + 0.49421E+05, 0.51399E+05, 0.53439E+05, 0.55542E+05, 0.57709E+05,
     + 0.59942E+05, 0.62242E+05, 0.64611E+05, 0.67051E+05, 0.69563E+05,
     + 0.72148E+05, 0.74808E+05, 0.77545E+05, 0.80360E+05, 0.83255E+05,
     + 0.86232E+05, 0.89291E+05, 0.92435E+05, 0.95667E+05, 0.98986E+05,
     + 0.10240E+06, 0.10590E+06, 0.10949E+06, 0.11318E+06, 0.11697E+06,
     + 0.12086E+06, 0.12484E+06, 0.12893E+06, 0.13313E+06, 0.13743E+06,
     + 0.14184E+06, 0.14637E+06, 0.15100E+06, 0.15575E+06, 0.16062E+06,
     + 0.16560E+06, 0.17071E+06, 0.17594E+06, 0.18129E+06, 0.18677E+06,
     + 0.19238E+06, 0.19813E+06, 0.20400E+06, 0.21002E+06, 0.21617E+06,
     + 0.22246E+06, 0.22890E+06, 0.23548E+06, 0.24221E+06, 0.24909E+06,
     + 0.25612E+06, 0.26331E+06, 0.27065E+06, 0.27816E+06, 0.28583E+06,
     + 0.29366E+06/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_HCOOH (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 4./
c...    HCOOH
c...        --       126
      data (QofT( 1,J),J=1,119)/ 0.31899E+04, 0.53773E+04, 0.79205E+04,
     + 0.10792E+05, 0.13993E+05, 0.17550E+05, 0.21509E+05, 0.25930E+05,
     + 0.30885E+05, 0.36460E+05, 0.42750E+05, 0.49864E+05, 0.57926E+05,
     + 0.67071E+05, 0.77453E+05, 0.89243E+05, 0.10263E+06, 0.11783E+06,
     + 0.13507E+06, 0.15462E+06, 0.17676E+06, 0.20183E+06, 0.23018E+06,
     + 0.26221E+06, 0.29836E+06, 0.33911E+06, 0.38501E+06, 0.43664E+06,
     + 0.49467E+06, 0.55981E+06, 0.63286E+06, 0.71470E+06, 0.80628E+06,
     + 0.90865E+06, 0.10230E+07, 0.11505E+07, 0.12927E+07, 0.14509E+07,
     + 0.16269E+07, 0.18225E+07, 0.20396E+07, 0.22804E+07, 0.25472E+07,
     + 0.28425E+07, 0.31692E+07, 0.35301E+07, 0.39285E+07, 0.43681E+07,
     + 0.48525E+07, 0.53858E+07, 0.59727E+07, 0.66178E+07, 0.73265E+07,
     + 0.81042E+07, 0.89571E+07, 0.98918E+07, 0.10915E+08, 0.12035E+08,
     + 0.13259E+08, 0.14597E+08, 0.16057E+08, 0.17650E+08, 0.19387E+08,
     + 0.21279E+08, 0.23339E+08, 0.25579E+08, 0.28016E+08, 0.30663E+08,
     + 0.33536E+08, 0.36655E+08, 0.40037E+08, 0.43701E+08, 0.47671E+08,
     + 0.51967E+08, 0.56614E+08, 0.61639E+08, 0.67068E+08, 0.72930E+08,
     + 0.79257E+08, 0.86082E+08, 0.93439E+08, 0.10137E+09, 0.10990E+09,
     + 0.11909E+09, 0.12898E+09, 0.13960E+09, 0.15102E+09, 0.16329E+09,
     + 0.17646E+09, 0.19059E+09, 0.20575E+09, 0.22200E+09, 0.23941E+09,
     + 0.25806E+09, 0.27802E+09, 0.29938E+09, 0.32223E+09, 0.34666E+09,
     + 0.37276E+09, 0.40064E+09, 0.43041E+09, 0.46218E+09, 0.49607E+09,
     + 0.53221E+09, 0.57074E+09, 0.61179E+09, 0.65551E+09, 0.70206E+09,
     + 0.75159E+09, 0.80430E+09, 0.86034E+09, 0.91992E+09, 0.98324E+09,
     + 0.10505E+10, 0.11219E+10, 0.11977E+10, 0.12782E+10, 0.13635E+10,
     + 0.14540E+10/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_HO2   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 2./
c...      HO2
c...        --       166
      data (QofT( 1,J),J=1,119)/ 0.39277E+03, 0.66062E+03, 0.97123E+03,
     + 0.13194E+04, 0.17014E+04, 0.21148E+04, 0.25578E+04, 0.30296E+04,
     + 0.35297E+04, 0.40585E+04, 0.46167E+04, 0.52055E+04, 0.58264E+04,
     + 0.64809E+04, 0.71707E+04, 0.78978E+04, 0.86641E+04, 0.94715E+04,
     + 0.10322E+05, 0.11218E+05, 0.12161E+05, 0.13154E+05, 0.14198E+05,
     + 0.15296E+05, 0.16449E+05, 0.17661E+05, 0.18933E+05, 0.20267E+05,
     + 0.21666E+05, 0.23133E+05, 0.24669E+05, 0.26277E+05, 0.27960E+05,
     + 0.29720E+05, 0.31560E+05, 0.33482E+05, 0.35489E+05, 0.37584E+05,
     + 0.39769E+05, 0.42048E+05, 0.44423E+05, 0.46898E+05, 0.49475E+05,
     + 0.52157E+05, 0.54948E+05, 0.57850E+05, 0.60868E+05, 0.64003E+05,
     + 0.67261E+05, 0.70643E+05, 0.74154E+05, 0.77797E+05, 0.81575E+05,
     + 0.85492E+05, 0.89553E+05, 0.93760E+05, 0.98118E+05, 0.10263E+06,
     + 0.10730E+06, 0.11213E+06, 0.11713E+06, 0.12230E+06, 0.12765E+06,
     + 0.13317E+06, 0.13888E+06, 0.14478E+06, 0.15086E+06, 0.15715E+06,
     + 0.16363E+06, 0.17032E+06, 0.17723E+06, 0.18434E+06, 0.19168E+06,
     + 0.19924E+06, 0.20704E+06, 0.21506E+06, 0.22333E+06, 0.23185E+06,
     + 0.24061E+06, 0.24963E+06, 0.25891E+06, 0.26846E+06, 0.27828E+06,
     + 0.28838E+06, 0.29876E+06, 0.30943E+06, 0.32039E+06, 0.33166E+06,
     + 0.34323E+06, 0.35512E+06, 0.36732E+06, 0.37985E+06, 0.39271E+06,
     + 0.40590E+06, 0.41944E+06, 0.43333E+06, 0.44758E+06, 0.46219E+06,
     + 0.47717E+06, 0.49252E+06, 0.50826E+06, 0.52439E+06, 0.54091E+06,
     + 0.55784E+06, 0.57518E+06, 0.59293E+06, 0.61112E+06, 0.62973E+06,
     + 0.64878E+06, 0.66828E+06, 0.68824E+06, 0.70866E+06, 0.72955E+06,
     + 0.75091E+06, 0.77276E+06, 0.79511E+06, 0.81795E+06, 0.84131E+06,
     + 0.86518E+06/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_O     (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 0./
c...        O
c...        --         6
      data (QofT( 1,J),J=1,119)/ 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     + 0.00000E+00/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_ClONO2(                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 2), QofT( 2,119),Q(NT)
      data xgj/ 12.,12./
c...   ClONO2
c...        --      5646
      data (QofT( 1,J),J=1,119)/ 0.11444E+06, 0.21121E+06, 0.34858E+06,
     + 0.53934E+06, 0.80041E+06, 0.11539E+07, 0.16286E+07, 0.22614E+07,
     + 0.30992E+07, 0.42015E+07, 0.56426E+07, 0.75152E+07, 0.99344E+07,
     + 0.13042E+08, 0.17012E+08, 0.22058E+08, 0.28437E+08, 0.36463E+08,
     + 0.46514E+08, 0.59042E+08, 0.74589E+08, 0.93801E+08, 0.11744E+09,
     + 0.14643E+09, 0.18181E+09, 0.22486E+09, 0.27705E+09, 0.34009E+09,
     + 0.41598E+09, 0.50705E+09, 0.61599E+09, 0.74590E+09, 0.90037E+09,
     + 0.10835E+10, 0.13001E+10, 0.15554E+10, 0.18556E+10, 0.22079E+10,
     + 0.26200E+10, 0.31012E+10, 0.36615E+10, 0.43126E+10, 0.50675E+10,
     + 0.59409E+10, 0.69492E+10, 0.81110E+10, 0.94469E+10, 0.10980E+11,
     + 0.12736E+11, 0.14745E+11, 0.17037E+11, 0.19649E+11, 0.22620E+11,
     + 0.25994E+11, 0.29819E+11, 0.34150E+11, 0.39044E+11, 0.44568E+11,
     + 0.50794E+11, 0.57799E+11, 0.65672E+11, 0.74506E+11, 0.84408E+11,
     + 0.95490E+11, 0.10788E+12, 0.12171E+12, 0.13713E+12, 0.15431E+12,
     + 0.17342E+12, 0.19465E+12, 0.21822E+12, 0.24435E+12, 0.27329E+12,
     + 0.30530E+12, 0.34069E+12, 0.37976E+12, 0.42286E+12, 0.47034E+12,
     + 0.52262E+12, 0.58012E+12, 0.64330E+12, 0.71267E+12, 0.78875E+12,
     + 0.87214E+12, 0.96344E+12, 0.10633E+13, 0.11725E+13, 0.12918E+13,
     + 0.14220E+13, 0.15640E+13, 0.17188E+13, 0.18873E+13, 0.20706E+13,
     + 0.22700E+13, 0.24866E+13, 0.27218E+13, 0.29771E+13, 0.32538E+13,
     + 0.35537E+13, 0.38784E+13, 0.42299E+13, 0.46100E+13, 0.50208E+13,
     + 0.54645E+13, 0.59435E+13, 0.64603E+13, 0.70175E+13, 0.76180E+13,
     + 0.82647E+13, 0.89608E+13, 0.97097E+13, 0.10515E+14, 0.11380E+14,
     + 0.12310E+14, 0.13307E+14, 0.14378E+14, 0.15526E+14, 0.16756E+14,
     + 0.18075E+14/
c...        --      7646
      data (QofT( 2,J),J=1,119)/ 0.11735E+06, 0.21659E+06, 0.35745E+06,
     + 0.55307E+06, 0.82078E+06, 0.11833E+07, 0.16700E+07, 0.23189E+07,
     + 0.31781E+07, 0.43084E+07, 0.57862E+07, 0.77065E+07, 0.10187E+08,
     + 0.13374E+08, 0.17445E+08, 0.22619E+08, 0.29161E+08, 0.37391E+08,
     + 0.47698E+08, 0.60545E+08, 0.76487E+08, 0.96188E+08, 0.12043E+09,
     + 0.15015E+09, 0.18644E+09, 0.23059E+09, 0.28410E+09, 0.34874E+09,
     + 0.42657E+09, 0.51995E+09, 0.63167E+09, 0.76489E+09, 0.92329E+09,
     + 0.11111E+10, 0.13331E+10, 0.15950E+10, 0.19029E+10, 0.22641E+10,
     + 0.26867E+10, 0.31801E+10, 0.37547E+10, 0.44224E+10, 0.51965E+10,
     + 0.60921E+10, 0.71261E+10, 0.83174E+10, 0.96873E+10, 0.11260E+11,
     + 0.13061E+11, 0.15120E+11, 0.17471E+11, 0.20149E+11, 0.23196E+11,
     + 0.26656E+11, 0.30578E+11, 0.35019E+11, 0.40038E+11, 0.45703E+11,
     + 0.52087E+11, 0.59270E+11, 0.67343E+11, 0.76403E+11, 0.86556E+11,
     + 0.97921E+11, 0.11062E+12, 0.12481E+12, 0.14062E+12, 0.15824E+12,
     + 0.17783E+12, 0.19961E+12, 0.22377E+12, 0.25057E+12, 0.28024E+12,
     + 0.31308E+12, 0.34936E+12, 0.38943E+12, 0.43362E+12, 0.48232E+12,
     + 0.53593E+12, 0.59489E+12, 0.65968E+12, 0.73081E+12, 0.80883E+12,
     + 0.89434E+12, 0.98797E+12, 0.10904E+13, 0.12024E+13, 0.13247E+13,
     + 0.14582E+13, 0.16038E+13, 0.17625E+13, 0.19353E+13, 0.21233E+13,
     + 0.23278E+13, 0.25499E+13, 0.27911E+13, 0.30528E+13, 0.33366E+13,
     + 0.36442E+13, 0.39772E+13, 0.43376E+13, 0.47273E+13, 0.51486E+13,
     + 0.56036E+13, 0.60948E+13, 0.66248E+13, 0.71962E+13, 0.78119E+13,
     + 0.84751E+13, 0.91889E+13, 0.99569E+13, 0.10783E+14, 0.11670E+14,
     + 0.12623E+14, 0.13646E+14, 0.14744E+14, 0.15921E+14, 0.17183E+14,
     + 0.18535E+14/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_NOp   (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 1), QofT( 1,119),Q(NT)
      data xgj/ 3./
c...      NO+
c...        --        46
      data (QofT( 1,J),J=1,119)/ 0.63956E+02, 0.90185E+02, 0.11642E+03,
     + 0.14265E+03, 0.16889E+03, 0.19513E+03, 0.22138E+03, 0.24763E+03,
     + 0.27388E+03, 0.30013E+03, 0.32639E+03, 0.35266E+03, 0.37894E+03,
     + 0.40523E+03, 0.43155E+03, 0.45790E+03, 0.48429E+03, 0.51074E+03,
     + 0.53725E+03, 0.56383E+03, 0.59052E+03, 0.61731E+03, 0.64422E+03,
     + 0.67127E+03, 0.69846E+03, 0.72582E+03, 0.75335E+03, 0.78108E+03,
     + 0.80901E+03, 0.83715E+03, 0.86552E+03, 0.89413E+03, 0.92298E+03,
     + 0.95208E+03, 0.98144E+03, 0.10111E+04, 0.10410E+04, 0.10712E+04,
     + 0.11017E+04, 0.11325E+04, 0.11636E+04, 0.11950E+04, 0.12268E+04,
     + 0.12588E+04, 0.12912E+04, 0.13239E+04, 0.13570E+04, 0.13903E+04,
     + 0.14241E+04, 0.14581E+04, 0.14926E+04, 0.15273E+04, 0.15624E+04,
     + 0.15979E+04, 0.16337E+04, 0.16699E+04, 0.17065E+04, 0.17434E+04,
     + 0.17806E+04, 0.18183E+04, 0.18563E+04, 0.18947E+04, 0.19334E+04,
     + 0.19725E+04, 0.20120E+04, 0.20519E+04, 0.20921E+04, 0.21327E+04,
     + 0.21737E+04, 0.22151E+04, 0.22568E+04, 0.22990E+04, 0.23415E+04,
     + 0.23844E+04, 0.24276E+04, 0.24713E+04, 0.25153E+04, 0.25598E+04,
     + 0.26046E+04, 0.26497E+04, 0.26953E+04, 0.27413E+04, 0.27876E+04,
     + 0.28343E+04, 0.28815E+04, 0.29290E+04, 0.29769E+04, 0.30251E+04,
     + 0.30738E+04, 0.31229E+04, 0.31723E+04, 0.32222E+04, 0.32724E+04,
     + 0.33230E+04, 0.33740E+04, 0.34254E+04, 0.34772E+04, 0.35294E+04,
     + 0.35819E+04, 0.36349E+04, 0.36883E+04, 0.37420E+04, 0.37961E+04,
     + 0.38507E+04, 0.39056E+04, 0.39609E+04, 0.40166E+04, 0.40727E+04,
     + 0.41292E+04, 0.41861E+04, 0.42434E+04, 0.43010E+04, 0.43591E+04,
     + 0.44176E+04, 0.44764E+04, 0.45357E+04, 0.45953E+04, 0.46554E+04,
     + 0.47158E+04/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_HOBr  (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 2), QofT( 2,119),Q(NT)
      data xgj/ 8.,8./
c...     HOBr
c...        --       169
      data (QofT( 1,J),J=1,119)/ 0.24445E+04, 0.41206E+04, 0.60683E+04,
     + 0.82610E+04, 0.10689E+05, 0.13352E+05, 0.16261E+05, 0.19427E+05,
     + 0.22867E+05, 0.26600E+05, 0.30643E+05, 0.35018E+05, 0.39745E+05,
     + 0.44844E+05, 0.50338E+05, 0.56249E+05, 0.62599E+05, 0.69410E+05,
     + 0.76706E+05, 0.84509E+05, 0.92845E+05, 0.10174E+06, 0.11121E+06,
     + 0.12128E+06, 0.13199E+06, 0.14335E+06, 0.15540E+06, 0.16815E+06,
     + 0.18165E+06, 0.19591E+06, 0.21096E+06, 0.22684E+06, 0.24358E+06,
     + 0.26120E+06, 0.27974E+06, 0.29922E+06, 0.31969E+06, 0.34118E+06,
     + 0.36372E+06, 0.38735E+06, 0.41210E+06, 0.43800E+06, 0.46511E+06,
     + 0.49345E+06, 0.52307E+06, 0.55400E+06, 0.58628E+06, 0.61997E+06,
     + 0.65509E+06, 0.69170E+06, 0.72984E+06, 0.76954E+06, 0.81087E+06,
     + 0.85386E+06, 0.89856E+06, 0.94502E+06, 0.99329E+06, 0.10434E+07,
     + 0.10955E+07, 0.11495E+07, 0.12055E+07, 0.12636E+07, 0.13238E+07,
     + 0.13862E+07, 0.14508E+07, 0.15177E+07, 0.15870E+07, 0.16587E+07,
     + 0.17328E+07, 0.18095E+07, 0.18888E+07, 0.19707E+07, 0.20554E+07,
     + 0.21428E+07, 0.22331E+07, 0.23263E+07, 0.24225E+07, 0.25217E+07,
     + 0.26241E+07, 0.27296E+07, 0.28385E+07, 0.29506E+07, 0.30662E+07,
     + 0.31853E+07, 0.33079E+07, 0.34341E+07, 0.35641E+07, 0.36979E+07,
     + 0.38355E+07, 0.39771E+07, 0.41228E+07, 0.42725E+07, 0.44265E+07,
     + 0.45848E+07, 0.47474E+07, 0.49145E+07, 0.50862E+07, 0.52624E+07,
     + 0.54435E+07, 0.56293E+07, 0.58201E+07, 0.60159E+07, 0.62168E+07,
     + 0.64229E+07, 0.66343E+07, 0.68511E+07, 0.70734E+07, 0.73013E+07,
     + 0.75349E+07, 0.77742E+07, 0.80196E+07, 0.82709E+07, 0.85283E+07,
     + 0.87920E+07, 0.90620E+07, 0.93385E+07, 0.96215E+07, 0.99112E+07,
     + 0.10208E+08/
c...        --       161
      data (QofT( 2,J),J=1,119)/ 0.24350E+04, 0.41047E+04, 0.60448E+04,
     + 0.82291E+04, 0.10648E+05, 0.13301E+05, 0.16200E+05, 0.19355E+05,
     + 0.22784E+05, 0.26504E+05, 0.30534E+05, 0.34895E+05, 0.39607E+05,
     + 0.44691E+05, 0.50169E+05, 0.56063E+05, 0.62394E+05, 0.69186E+05,
     + 0.76461E+05, 0.84243E+05, 0.92555E+05, 0.10142E+06, 0.11087E+06,
     + 0.12091E+06, 0.13159E+06, 0.14292E+06, 0.15494E+06, 0.16766E+06,
     + 0.18112E+06, 0.19534E+06, 0.21036E+06, 0.22620E+06, 0.24289E+06,
     + 0.26047E+06, 0.27896E+06, 0.29840E+06, 0.31882E+06, 0.34025E+06,
     + 0.36274E+06, 0.38630E+06, 0.41099E+06, 0.43683E+06, 0.46387E+06,
     + 0.49215E+06, 0.52169E+06, 0.55255E+06, 0.58475E+06, 0.61836E+06,
     + 0.65340E+06, 0.68992E+06, 0.72796E+06, 0.76757E+06, 0.80880E+06,
     + 0.85169E+06, 0.89628E+06, 0.94263E+06, 0.99079E+06, 0.10408E+07,
     + 0.10927E+07, 0.11466E+07, 0.12025E+07, 0.12605E+07, 0.13205E+07,
     + 0.13828E+07, 0.14472E+07, 0.15140E+07, 0.15831E+07, 0.16546E+07,
     + 0.17286E+07, 0.18051E+07, 0.18842E+07, 0.19660E+07, 0.20504E+07,
     + 0.21377E+07, 0.22277E+07, 0.23207E+07, 0.24167E+07, 0.25157E+07,
     + 0.26178E+07, 0.27231E+07, 0.28317E+07, 0.29436E+07, 0.30589E+07,
     + 0.31777E+07, 0.33001E+07, 0.34260E+07, 0.35557E+07, 0.36892E+07,
     + 0.38265E+07, 0.39678E+07, 0.41131E+07, 0.42626E+07, 0.44162E+07,
     + 0.45741E+07, 0.47364E+07, 0.49031E+07, 0.50744E+07, 0.52503E+07,
     + 0.54309E+07, 0.56164E+07, 0.58067E+07, 0.60021E+07, 0.62025E+07,
     + 0.64081E+07, 0.66191E+07, 0.68354E+07, 0.70572E+07, 0.72846E+07,
     + 0.75177E+07, 0.77565E+07, 0.80013E+07, 0.82521E+07, 0.85090E+07,
     + 0.87721E+07, 0.90415E+07, 0.93173E+07, 0.95997E+07, 0.98888E+07,
     + 0.10185E+08/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_C2H4  (                       
     I T, 	! temperature in K 
     I iso, 	! isotope code (HITRAN INDEX)
     O gsi, 	! state independent nuclear degeneracyfactor
     O QT) 	! Total Internal Partition Function
 
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 2), QofT( 2,119),Q(NT)
      data xgj/ 1.,2./
c...     C2H4
c...        --       221
      data (QofT( 1,J),J=1,119)/ 0.95843E+03, 0.16137E+04, 0.23744E+04,
     + 0.32285E+04, 0.41694E+04, 0.51963E+04, 0.63143E+04, 0.75337E+04,
     + 0.88702E+04, 0.10344E+05, 0.11978E+05, 0.13802E+05, 0.15846E+05,
     + 0.18145E+05, 0.20740E+05, 0.23675E+05, 0.27000E+05, 0.30770E+05,
     + 0.35048E+05, 0.39905E+05, 0.45420E+05, 0.51680E+05, 0.58786E+05,
     + 0.66850E+05, 0.75997E+05, 0.86369E+05, 0.98123E+05, 0.11144E+06,
     + 0.12651E+06, 0.14356E+06, 0.16284E+06, 0.18463E+06, 0.20923E+06,
     + 0.23699E+06, 0.26831E+06, 0.30360E+06, 0.34334E+06, 0.38808E+06,
     + 0.43840E+06, 0.49495E+06, 0.55847E+06, 0.62976E+06, 0.70973E+06,
     + 0.79935E+06, 0.89973E+06, 0.10121E+07, 0.11378E+07, 0.12782E+07,
     + 0.14351E+07, 0.16102E+07, 0.18055E+07, 0.20231E+07, 0.22656E+07,
     + 0.25354E+07, 0.28356E+07, 0.31692E+07, 0.35398E+07, 0.39511E+07,
     + 0.44074E+07, 0.49132E+07, 0.54736E+07, 0.60940E+07, 0.67803E+07,
     + 0.75392E+07, 0.83776E+07, 0.93035E+07, 0.10325E+08, 0.11452E+08,
     + 0.12694E+08, 0.14062E+08, 0.15567E+08, 0.17224E+08, 0.19045E+08,
     + 0.21046E+08, 0.23243E+08, 0.25655E+08, 0.28300E+08, 0.31200E+08,
     + 0.34377E+08, 0.37856E+08, 0.41662E+08, 0.45826E+08, 0.50378E+08,
     + 0.55351E+08, 0.60781E+08, 0.66707E+08, 0.73172E+08, 0.80219E+08,
     + 0.87899E+08, 0.96262E+08, 0.10537E+09, 0.11527E+09, 0.12604E+09,
     + 0.13775E+09, 0.15047E+09, 0.16428E+09, 0.17927E+09, 0.19553E+09,
     + 0.21316E+09, 0.23226E+09, 0.25296E+09, 0.27537E+09, 0.29963E+09,
     + 0.32587E+09, 0.35425E+09, 0.38492E+09, 0.41805E+09, 0.45383E+09,
     + 0.49246E+09, 0.53413E+09, 0.57908E+09, 0.62754E+09, 0.67977E+09,
     + 0.73602E+09, 0.79660E+09, 0.86179E+09, 0.93194E+09, 0.10074E+10,
     + 0.10885E+10/
c...        --       231
      data (QofT( 2,J),J=1,119)/ 0.39228E+04, 0.66051E+04, 0.97190E+04,
     + 0.13215E+05, 0.17066E+05, 0.21270E+05, 0.25846E+05, 0.30838E+05,
     + 0.36309E+05, 0.42341E+05, 0.49032E+05, 0.56496E+05, 0.64862E+05,
     + 0.74275E+05, 0.84897E+05, 0.96912E+05, 0.11052E+06, 0.12595E+06,
     + 0.14347E+06, 0.16335E+06, 0.18592E+06, 0.21155E+06, 0.24064E+06,
     + 0.27365E+06, 0.31109E+06, 0.35354E+06, 0.40166E+06, 0.45615E+06,
     + 0.51785E+06, 0.58765E+06, 0.66657E+06, 0.75575E+06, 0.85646E+06,
     + 0.97011E+06, 0.10983E+07, 0.12428E+07, 0.14055E+07, 0.15886E+07,
     + 0.17945E+07, 0.20260E+07, 0.22861E+07, 0.25779E+07, 0.29052E+07,
     + 0.32721E+07, 0.36830E+07, 0.41429E+07, 0.46573E+07, 0.52323E+07,
     + 0.58744E+07, 0.65912E+07, 0.73906E+07, 0.82816E+07, 0.92740E+07,
     + 0.10379E+08, 0.11607E+08, 0.12973E+08, 0.14490E+08, 0.16174E+08,
     + 0.18042E+08, 0.20112E+08, 0.22406E+08, 0.24945E+08, 0.27755E+08,
     + 0.30861E+08, 0.34293E+08, 0.38083E+08, 0.42266E+08, 0.46878E+08,
     + 0.51961E+08, 0.57560E+08, 0.63724E+08, 0.70504E+08, 0.77959E+08,
     + 0.86150E+08, 0.95145E+08, 0.10502E+09, 0.11585E+09, 0.12772E+09,
     + 0.14072E+09, 0.15496E+09, 0.17054E+09, 0.18759E+09, 0.20622E+09,
     + 0.22658E+09, 0.24880E+09, 0.27306E+09, 0.29952E+09, 0.32837E+09,
     + 0.35981E+09, 0.39404E+09, 0.43131E+09, 0.47186E+09, 0.51595E+09,
     + 0.56387E+09, 0.61594E+09, 0.67247E+09, 0.73382E+09, 0.80038E+09,
     + 0.87255E+09, 0.95076E+09, 0.10355E+10, 0.11272E+10, 0.12265E+10,
     + 0.13339E+10, 0.14501E+10, 0.15756E+10, 0.17113E+10, 0.18577E+10,
     + 0.20159E+10, 0.21865E+10, 0.23705E+10, 0.25688E+10, 0.27826E+10,
     + 0.30129E+10, 0.32608E+10, 0.35277E+10, 0.38149E+10, 0.41237E+10,
     + 0.44557E+10/
  
      eps=0.01
c
      gsi = xgj(iso)
	do I=1,NT
	  Q(I)=QofT(iso,I)
	enddo
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
        Qt = -1.
        write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c
c***************************
      SUBROUTINE AtoB(aa,bb,A,B,npt)
c***************************
c...LaGrange 3- and 4-point interpolation
c...arrays A and B are the npt data points,  given aa, a value of the 
c...A variable, the routine will find the corresponding bb value
c
c...input:  aa
c...output: bb 
      Parameter (Nmax=600)
      dimension A(Nmax),B(Nmax)
c
C 
c
      DO 50 I=2,npt
      IF(A(I).GE.aa)THEN 
      IF(I.LT.3 .OR. I.EQ.npt) THEN
C     LaGrange three point interpolation 
      J = I
      IF(I.LT.3) J = 3
      IF(I.EQ.npT) J = npt
c.....do not devide by zero
            A0D1=A(J-2)-A(J-1)
            IF(A0D1.EQ.0.) A0D1=0.0001
            A0D2=A(J-2)-A(J)
            IF(A0D2.EQ.0.) A0D2=0.0001
            A1D1=A(J-1)-A(J-2)
            IF(A1D1.EQ.0.) A1D1=0.0001
            A1D2=A(J-1)-A(J)
            IF(A1D2.EQ.0.) A1D2=0.0001
            A2D1=A(J)-A(J-2)
            IF(A2D1.EQ.0.) A2D1=0.0001
            A2D2=A(J)-A(J-1)
            IF(A2D2.EQ.0.) A2D2=0.0001
c
      A0=(aa-A(J-1))*(aa-A(J))/(A0D1*A0D2)
      A1=(aa-A(J-2))*(aa-A(J))/(A1D1*A1D2)
      A2=(aa-A(J-2))*(aa-A(J-1))/(A2D1*A2D2)
c
      bb = A0*B(J-2) + A1*B(J-1) + A2*B(J)
c
      ELSE
C     LaGrange four point interpolation 
      J = I
c.....do not devide by zero
            A0D1=A(J-2)-A(J-1)
            IF(A0D1.EQ.0.) A0D1=0.0001
            A0D2=A(J-2)-A(J)
            IF(A0D2.EQ.0.) A0D2=0.0001
            A0D3 = (A(J-2)-A(J+1))
            IF(A0D3.EQ.0.) A0D3=0.0001
c
            A1D1=A(J-1)-A(J-2)
            IF(A1D1.EQ.0.) A1D1=0.0001
            A1D2=A(J-1)-A(J)
            IF(A1D2.EQ.0.) A1D2=0.0001
            A1D3 = A(J-1)-A(J+1)
            IF(A1D3.EQ.0.) A1D3=0.0001
c
            A2D1=A(J)-A(J-2)
            IF(A2D1.EQ.0.) A2D1=0.0001
            A2D2=A(J)-A(J-1)
            IF(A2D2.EQ.0.) A2D2=0.0001
            A2D3 = A(J)-A(J+1)
            IF(A2D3.EQ.0.) A2D3=0.0001
c
            A3D1 = A(J+1)-A(J-2)
            IF(A3D1.EQ.0.) A3D1=0.0001
            A3D2 = A(J+1)-A(J-1)
            IF(A3D2.EQ.0.) A3D2=0.0001
            A3D3 = A(J+1)-A(J)
            IF(A3D3.EQ.0.) A3D3=0.0001
c
      A0=(aa-A(J-1))*(aa-A(J))*(aa-A(J+1))
      A0=A0/(A0D1*A0D2*A0D3)
      A1=(aa-A(J-2))*(aa-A(J))*(aa-A(J+1))
      A1=A1/(A1D1*A1D2*A1D3)
      A2=(aa-A(J-2))*(aa-A(J-1))*(aa-A(J+1))
      A2=A2/(A2D1*A2D2*A2D3)
      A3=(aa-A(J-2))*(aa-A(J-1))*(aa-A(J))
      A3=A3/(A3D1*A3D2*A3D3)
c
      bb = A0*B(J-2) + A1*B(J-1) + A2*B(J) + A3*B(J+1)
      ENDIF 
c
      GO TO 100
      ENDIF 
   50 CONTINUE
  100 CONTINUE
ccc      write(2,*) 'F1, F2, F3, H1, H2, H3 =',B(J-2),B(J-1),B(J),
ccc     + A(J-2), A(J-1), A(J)
ccc      write(2,*) 'A0, A1, A2, bb =',A0,A1,A2,bb
c
      RETURN
      END 
c
c------------------------------------------------------------------------------
c
      Subroutine voigt_init (f1, f2, f3)
c*******************************************************************************
c       This subroutine initializes the functions f1, f2, f3 used in the 
c       decomposition of the voigt function.  
c       These functions are defined respectively over the regions (domains)
c       from 0 to 4, 16 and 64 voigt halfwidths from the line center.  At 
c       the end of the domain boundary, the value and first derivatives
c       are all zero.  
c       Each function is defined on a grid of nx points so each function has a
c       different frequency resolution.  
c       The calculations of the coefficients need to be performed in double 
c       precision, otherwise instabilities occur at large values of zeta.
c       The subfunctions f1,f2,f3,f4 can be calculated in single precision

c*******************************************************************************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
C                                                                        
C                  IMPLEMENTATION:    W.O. GALLERY                        
C                                                                        
C             ALGORITHM REVISIONS:    M.W. SHEPHARD & S.A. CLOUGH                        
C                                                                        
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.        
C                     131 Hartwell Ave,  Lexington,  MA   02421          
C                                                                         
C----------------------------------------------------------------------   

      parameter (nx=2001, nxmx=1.3*nx, nzeta=101, nq=3, ndom=4)

      implicit real*8 (a-e,g-h, o-z)

	dimension  x(0:nx), z(0:nx)
	dimension xf_voigt_armstr(0:nx)

c       Need fv and dfv at point nx in order to compute the second derivative 
c       as was previously done.  Store dfv at last 3 points for each of 3 separate domains
	dimension xfv(0:nx, 0:nzeta-1)
	dimension xfv_tmp(0:4, 0:ndom-1, 0:nzeta-1)
	dimension dfv(0:ndom-1, 0:nzeta-1)
c	dimension d2fv(0:ndom-1, 0:nzeta-1)
c       Store fv(nx-1) for domain 1: need it to calculate the a coefficients
c       For the other domains, it can be obtained from the fv for domain 4.  
	dimension domain(0:ndom-1)

c       f1, f2, f3 dimensioned nxmx=1.3*nx, extra points are to allow for 
c       subscript overrun without overwriting other arrays.

	dimension f1(0:nxmx-1, 0:nzeta-1), f2(0:nxmx-1, 0:nzeta-1),
	1    f3(0:nxmx-1, 0:nzeta-1)
c
	dimension q1(0:nx-1), q2(0:nx-1), q3(0:nx-1)
	dimension             AVC(0:nzeta)
	common /voigt_cf/
	1    a_1(0:nzeta-1), a_2(0:nzeta-1), a_3(0:nzeta-1),
	2    b_1(0:nzeta-1), b_2(0:nzeta-1), b_3(0:nzeta-1)
C
C       AVC transforms lorentz and doppler halfwidths to voigt
C       halfwidths: AV = AVC*(AL + AD)

        DATA AVC/                                                         
	1    .10000E+01,  .99535E+00,  .99073E+00,  .98613E+00,  .98155E+00, 
	2    .97700E+00,  .97247E+00,  .96797E+00,  .96350E+00,  .95905E+00, 
	3    .95464E+00,  .95025E+00,  .94589E+00,  .94156E+00,  .93727E+00,
	4    .93301E+00,  .92879E+00,  .92460E+00,  .92045E+00,  .91634E+00,
	5    .91227E+00,  .90824E+00,  .90425E+00,  .90031E+00,  .89641E+00,
	6    .89256E+00,  .88876E+00,  .88501E+00,  .88132E+00,  .87768E+00,
	7    .87410E+00,  .87058E+00,  .86712E+00,  .86372E+00,  .86039E+00,
	8    .85713E+00,  .85395E+00,  .85083E+00,  .84780E+00,  .84484E+00,
	9    .84197E+00,  .83919E+00,  .83650E+00,  .83390E+00,  .83141E+00,
	1    .82901E+00,  .82672E+00,  .82454E+00,  .82248E+00,  .82053E+00,
	2    .81871E+00,  .81702E+00,  .81547E+00,  .81405E+00,  .81278E+00,
	3    .81166E+00,  .81069E+00,  .80989E+00,  .80925E+00,  .80879E+00,
	4    .80851E+00,  .80842E+00,  .80852E+00,  .80882E+00,  .80932E+00,
	5    .81004E+00,  .81098E+00,  .81214E+00,  .81353E+00,  .81516E+00,
	6    .81704E+00,  .81916E+00,  .82154E+00,  .82418E+00,  .82708E+00,
	7    .83025E+00,  .83370E+00,  .83742E+00,  .84143E+00,  .84572E+00,
	8    .85029E+00,  .85515E+00,  .86030E+00,  .86573E+00,  .87146E+00,
	9    .87747E+00,  .88376E+00,  .89035E+00,  .89721E+00,  .90435E+00,
	1    .91176E+00,  .91945E+00,  .92741E+00,  .93562E+00,  .94409E+00,
	2    .95282E+00,  .96179E+00,  .97100E+00,  .98044E+00,  .99011E+00,
	3    .10000E+01,  .10000E+01/                                       

C       Define 3 domains: 4, 16, and 64 voigt halfwidths from the line center 
	DATA domain / 4.0, 16.0, 64.0, 256.0/

	pi = 2.0d0*asin(1.0d0)

c       This constant is related to Doppler halfwidth
	const = sqrt(log(2.0D0))

C************************************************************************
c       The scaled frequency variable x = const*v/AD 
c       x is required by the humlicek and Armstrong routine where
c       v is the frequency in wavenumbers and AD is the
c       Doppler halfwidth at half height. 
c       .  
c       For LBLRTM we want the function converted to the variable z = v/AV.
c       Thus, x=const*z*(AV/AD)
C************************************************************************
	
c       Loop over zeta from 0 to 1 in steps of 1/(nzeta-1) ;0.01
c       dble is just a function to convert from integer to double precision

	do 100 izet = 0,nzeta-1
	   zeta = dble(izet)/dble(nzeta-1)
	   
c       Force Doppler Halfwidth to 1.

	   if (zeta .lt. 1.00) then 
	      AD = 1.
	      AL = zeta/(1-zeta)
	   end if
	   
c       zeta = 1: pure lorenz, AD=0
	   if (zeta .eq. 1.00) then 
	      AD = 0.
	      AL = 1.00
	   end if
	   
c       Calculate the voigt halwidth from AVC, AD, AL
c       The following works only for nzeta=101
	   AV = AVC(izet)*(AD + AL)

c       AD1 is doppler width at the 1/e point
	   AD1 = AD/const
	   if (AD1 .ne. 0) then
	      ep = AL/AD1
	   else
	      ep = 0.0
	   end if
c       Normalization factors: 
c       Note: Fucntion returned is the generalized voigt function u(zeta,x)
c       Voigt line shape function = sqrt(log(2.0)/pi)/AD *u(zeta,x)

	   cnorm1 = const/sqrt(pi)*AV

c       Loop over the 3 domains 4, 16, and 64 voigt widths calculating f1, f2,
c       f3 respectively
c
	   do 70 id=0, ndom-2
c
C       SETUP PARAMETERS FOR CALL TO VOIGT GENERATOR 
	      
c       Define the range of x: nx frequency points equals domain(id) voigt
c       halfwidths. vmax = nvm*AV, xmax = const*domain(id)*AV/AD = dx*(nx-1)
c       dx = (const*AV)/AD*domain(id)/ REAL(nx-1)
c       where dx = (const*AV)/AD*dz;   therefore dz = domain(id)/ REAL(nx-1)

	      if (AD .ne. 0.0) then 
		 dx = (const*AV)/AD*(domain(id)/ REAL(nx-1))
	      else 
		 dx = (const*AV)*(domain(id)/ REAL(nx-1))
	      endif 

	      dz=domain(id)/ REAL(nx-1)

	      do 20 i=0,nx
		 x(i) = dble(i)*dx
		 z(i) = dble(i)*domain(id)/dble(nx-1)
 20	      continue

C*****************************************************************************
C       Compute the Voigt Function
C*****************************************************************************

C       FOR ARMSTRONG:
C       see paper by B. H. Armstrong, J. Quant. Spectrosc. Radiat. Transfer. 
c       Vol. 7, pp 61-88, 1967 for the Armstrong routine
c       returns back a real part of the voigt function "v_arm"
	      
	      call armstrong(nx+1,x,zeta,xf_voigt_armstr)
	      
c       Normalize the Voigt Profile and get the first derivative
c       The scaled frequency variable x = const*v/AD is required by the humlicek 
c       and Armstrong routine.  We want it scaled by z by using the cnorm1
c       (first derivative actually only needed at 3 points: i=nx-1, nx+1, nx)
	      if (zeta .lt. 1.00) then
		  do 30 i = 0,nx
		      xfv(i,izet) = xf_voigt_armstr(i)*cnorm1
  30		  continue

c       First derivative computed using ( symmetric finite difference = f(y+x)-f(y-x)/2x )
		  do 36 i=0,2
		      i2=i+nx-2
		      xfv_tmp(i,id,izet) = xfv(i2,izet)
   36		  continue

		  dfv(id,izet) =.5*(xfv_tmp(2,id,izet) 
	1	      - xfv_tmp(0,id,izet))/(dz)	
	      endif
c.................................................................................
c       FOR ZETA=1
c       Just use the straight Lorenz function at zeta=1 (AD=0, AL=1) 
c       xfv(i,izet) = (1./pi)*(1./(1.+z(i)**2))
c.................................................................................
	      
	      if (zeta .eq. 1.00) then 
		  do 40 i = 0,nx
			  xfv(i,izet) = (1./pi)*(1./(1.+z(i)**2))
  40		  continue

c       Pick the last few points as this is where you want the functions 
c       to merge (4,16,64,256)
c
		  do 46 i=0,2
		      i2=i+nx-2
		      xfv_tmp(i,id,izet) = xfv(i2,izet)
   46		  continue

		  dfv(id,izet) =.5*(xfv_tmp(2,id,izet) 
	1	      - xfv_tmp(0,id,izet))/(dz)	
	      endif

c       For now, set xfvbnd = 0
c	      xfvbnd(izet) = 0.0

c       Now decompose voigt profile into subfunctions

	      if (id .eq.0) then
		  b_1(izet) = dfv(0,izet)/(2*domain(0))
		  a_1(izet) = xfv(nx-1,izet) - b_1(izet)
	1	      *domain(0)**2 
	      else if (id .eq. 1) then
		  b_2(izet) = dfv(1,izet)/(2*domain(1))
		  a_2(izet) = xfv(nx-1,izet) - b_2(izet)
	1	      *domain(1)**2 
	      else if (id .eq. 2) then
		  b_3(izet) = dfv(2,izet)/(2*domain(2))
		  a_3(izet) = xfv(nx-1,izet) - b_3(izet)
	1	      *domain(2)**2
	      endif

c       Compute the subfunctions f1, f2, f3, f4 (f4 not used but carried for debugging)
c       Note: the subfunctions should equal zero at the domain boundaries.
c       In practice, they may have a small finite value, possibly negative,
c       due to roundoff or truncation error.  Here they will be set
c       identically to zero by restricting the loop to nx-2. The assumption
c       is that the arrays are initialized to zero.
 
	      do 60 iz = 0,nx-2
c
		  q1(iz) = a_1(izet) + b_1(izet)*z(iz)**2 
		  q2(iz) = a_2(izet) + b_2(izet)*z(iz)**2 
		  q3(iz) = a_3(izet) + b_3(izet)*z(iz)**2 
c		  
		  if (id .eq. 0) then 
		      f1(iz, izet) = xfv(iz, izet) - q1(iz)
		  else if (id .eq. 1) then
		      if (z(iz) .le. domain(0)) then
			  f2(iz, izet) = q1(iz) - q2(iz)
		      else
			  f2(iz, izet) = xfv(iz, izet) - q2(iz)
		      endif
		  else if (id .eq. 2) then
		      if (z(iz) .le. domain(1)) then 
			  f3(iz, izet) = q2(iz) - q3(iz)
		      else 
			  f3(iz, izet) = xfv(iz, izet) - q3(iz)
		      endif 
		  endif
  60	      continue
c             end of iz loop
  70	  continue
c         end of id loop
100	continue
c       end of zeta loop

	end

C*************************************************************************
C*************************************************************************
	Subroutine armstrong(nx,xx,zeta,f_voigt)
C
C
c       Calculates the voigt function using the Armstrong algorithm
C
c       Parameters:  
c       arg1 = nx: number of points                                    in
c       arg2 = x : grid points: sqrt(ln(2.)/pi)/alpha_d*freq(in cm-1)  in
c       arg3 = zeta: Voigt parameter = alpha_l/(alpha_l+alpha_d)       in
c       arg4 = voigt function                                         out
C
c       Implemented and adapted by MWS Oct,2001
c       Obtained from WOG who adapted the code provided 
c       from B. H. Armstrong, J. Quant. Spectrosc. Radiat. Transfer. Vol. 7
c       pp 61-88, 1967
C
C************************************************************************
C************************************************************************

	implicit real*8 (a-h, o-z)
	REAL*8 K1, K2, K3
	
	COMMON /armstrng/ W(10),T(10)
c
	dimension xx(nx), f_voigt(nx)
	Dimension C(34)

c     y = Lorenz width / Doppler width (1/e point) 
	if (zeta .ne. 1.0) then
            y = sqrt(log(2.0d0))*zeta/(1.d0-zeta)
	else 
            y = 100.0
	endif 
c     
	do 100 iz=1, nx
	    x = xx(iz)
c     
	    IF(Y .LT. 1.0 .AND. X .LT. 4.0 .OR. Y .LT. 1.8/(X+1.0)) then 
c     
c     F3(T)=EXP(T**2-X**2)
		Y2=Y**2
		IF((X**2-Y2) .GT. 70.0) GO TO 1002
		U1=EXP(-X**2+Y2)*COS(2.0*X*Y)
                GO TO 1005
 1002           U1=0.0
 1005           IF(X.GT.5.0) GO TO 1100
C     FROM HERE TO STATEMENT 30 WE CALCULATE DAWSONS FUNCTION
C     ENTER HUMMERS CHEBYSHEV COEFFICIENTS C(I)
	data C/0.1999999999972224, -0.1840000000029998, 
	1    0.1558399999965025, -0.1216640000043988,
	2    0.0877081599940391, -0.0585141248086907, 
	3    0.0362157301623914, -0.0208497654398036, 
	4    0.0111960116346270, -0.56231896167109D-2,
	5    0.26487634172265D-2,-0.11732670757704D-2,
	6    0.4899519978088D-3, -0.1933630801528D-3, 
	7    0.722877446788D-4,  -0.256555124979D-4, 
	8    0.86620736841D-5,   -0.27876379719D-5, 
	9    0.8566873627D-6,    -0.2518433784D-6, 
	1    0.709360221D-7,     -0.191732257D-7,
	2    0.49801256D-9,      -0.12447734D-8,
	3    0.2997777D-9,       -0.696450D-10, 
	4    0.156262D-10,       -0.33897D-11, 
	5    0.7116D-12,         -0.1447D-12, 
	6    0.285D-13,          -0.55D-14,
	7    0.10D-14,           -0.2D-15/
C       CLENSHAWS ALGORITHM AS GIVEN BY HUMMER
	BN01=0.0D0
	BN02=0.0D0
	X1=X/5.0D0
	COEF=4.0D0*X1**2-2.0D0
	DO 20 I=1,34
	   II=35-I
	   BN=COEF*BN01-BN02+C(II)
	   BN02=BN01
	   BN01=BN
 20	CONTINUE
 	F=X1*(BN-BN02)
 	DN01=1.0D0-2.0D0*X*F
 	DN02=F
	GO TO 1200
 1100	DN01=-(.5/X**2+.75/X**4+1.875/X**6+6.5625/X**8+29.53125/X**10+
	1    162.4218/X**12+1055.7421/X**14)
	DN02 = (1.0D0-DN01)/(2.0D0*X)
 1200	FUNCT=Y*DN01
	IF(Y .LE. 1.0D-08) GO TO 1800
	Q=1.0
	YN=Y
	DO 1600 I=2,50
	   DN=(X*DN01+DN02)*(-2.)/ REAL(I)
	   DN02=DN01
	   DN01=DN
	   IF(MOD(I,2)) 1600, 1600, 1500
 1500	   Q=-Q
	   YN=YN*Y2
	   G=DN*YN
	   FUNCT=FUNCT+Q*G
	   IF(ABS(G/FUNCT) .LE. 1.0E-08) GO TO 1800
 1600	CONTINUE
 1800	K1=U1-1.12837917D0*FUNCT
c
	f_voigt(iz) = K1
c
	    else IF(Y .LT. 2.5 .AND. X .LT. 4.0) then 

	Y2=Y**2
	G=0.0
	DO 2100 I=1,10
	   R=T(I)-X
	   S=T(I)+X
	   G=G+(4.0*T(I)**2-2.0)*(R*ATAN(R/Y)+S*ATAN(S/Y)-
 	1	0.5*Y*(LOG(Y2+R**2)+
	2	LOG(Y2+S**2)))*W(I) 
 2100	CONTINUE
	K2=0.318309886D0*G 	
c
	f_voigt(iz) = K2
c
	    else 
	Y2=Y**2
	G=0.0
	DO 3100 I=1,10
	   G=G+(1.0D0/((X-T(I))**2+Y2)+1.0/((X+T(I))**2+Y2))*W(I)
 3100	CONTINUE
	K3=0.318309886D0*Y*G
c
	f_voigt(iz) = K3
c
	    endif 
c
 100	continue
	return 
	end 

c==============================================================================
	BLOCK DATA ARMCOEFF

	implicit real*8 (a-h, o-z)

c*******************************************************************************

	COMMON /armstrng/ W(10),T(10)
	DATA W/4.62243670E-1, 2.86675505E-1, 1.09017206E-1, 
	1      2.48105209E-2, 3.24377334E-3, 2.28338636E-4,  
	2      7.80255648E-6, 1.08606937E-7, 4.39934099E-10, 
	3      2.22939365E-13/
        DATA T/0.245340708, 0.737473729, 1.23407622, 1.73853771, 
	1      2.25497400,  2.78880606,  3.34785457, 3.94476404, 
	2      4.60368245,  5.38748089/
	END 

