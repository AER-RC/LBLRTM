C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      SUBROUTINE READAL(N)
C
C     Reads layer optical depth files in for three temperatures
C     and puts them into one array with dimension (30000,3)
C
      PARAMETER (NLAYR=33,NGNU=1000000,NTEMP=3)
C
      IMPLICIT DOUBLE PRECISION (V)
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID,HDATE,HTIME          & A03050
      DOUBLE PRECISION XID1,SECNT1,HMOLD1,XALTZ1,YID1
      DOUBLE PRECISION XID2,SECNT2,HMOLD2,XALTZ2,YID2
      DOUBLE PRECISION XID3,SECNT3,HMOLD3,XALTZ3,YID3
C
      DIMENSION XFILHD(2),PNLHD(2)                                        A09140
      DIMENSION IWD(2),IWD2(2)
      DIMENSION XOD1(2410),XOD2(2410),XOD3(2410)
C
      CHARACTER*55 PATH1,PATH2,PATH3
      CHARACTER*10 HFORM1,HFORM2,HFORM3
      CHARACTER CXID*80,CFORM*11,XID8*8,IDCNTL*6                          A03430
C
      COMMON /HDRF/ V1D,V2D,DVD,NLND,IWLD                                 A03240
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       A03070
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   A03080
     *                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF    A03090
      COMMON /EMIHD1/ XID1(10),SECNT1,PAVE1,TAVE1,HMOLD1(60),XALTZ1(4),
     *                WK1(60),PZL1,PZU1,TZL1,TZU1,WBROD1,DV1,V11,V21,
     *                TBOND1,
     *                EMISV1,FSCDD1(17),NMOL1,LAYDM1,YI11,YID1(10),
     *                LSTWD1
      COMMON /EMIHD2/ XID2(10),SECNT2,PAVE2,TAVE2,HMOLD2(60),XALTZ2(4),
     *                WK2(60),PZL2,PZU2,TZL2,TZU2,WBROD2,DV2,V12,V22,
     *                TBOND2,
     *                EMISV2,FSCDD2(17),NMOL2,LAYDM2,YI12,YID2(10),
     *                LSTWD2
      COMMON /EMIHD3/ XID3(10),SECNT3,PAVE3,TAVE3,HMOLD3(60),XALTZ3(4),
     *                WK3(60),PZL3,PZU3,TZL3,TZU3,WBROD3,DV3,V13,V23,
     *                TBOND3,
     *                EMISV3,FSCDD3(17),NMOL3,LAYDM3,YI13,YID3(10),
     *                LSTWD3
      COMMON /PNLHDR/ V1P,V2P,DVP,NLIMO                                     A09090
      COMMON /ABSCOF/ V1POD(5000),V2POD(5000),DVPOD(5000),NLIMOD(5000),
     *     NDONE,XOPTDP(NGNU,NTEMP)
      DIMENSION XFILH3(2)
      DIMENSION XFILH1(2)
      DIMENSION XFILH2(2)
C
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHD(1),V1P)                     A09160
      EQUIVALENCE (IWD(1),XID(1)) , (IWD2(1),V1D)
      EQUIVALENCE (FSCDID(5),IEMIT)
      EQUIVALENCE (XFILH3(1),XID3(1))
      EQUIVALENCE (XFILH1(1),XID1(1))
      EQUIVALENCE (XFILH2(1),XID2(1))
C
C     DATA statements for names of optical depth files
C
      DATA PATH1 / "ODint-15_" /
      DATA PATH2 / "ODint+00_" /
      DATA PATH3 / "ODint+15_" /
C
C     Determine word length for hardware
C
      NFHDRF = NWDL(IWD,LSTWDF)
      NPHDRF = NWDL(IWD2,IWLD)
C
C     Save incoming value of iemit
C
      IEMSV = IEMIT
C
C     Get format for layer number suffix addition to
C     BASE PATHNAME
C
      CALL QNTIFY(PATH1,HFORM1)
      CALL QNTIFY(PATH2,HFORM2)
      CALL QNTIFY(PATH3,HFORM3)
C
      CALL OPNOD(NLAYR,N,PATH1,HFORM1,16)
      CALL OPNOD(NLAYR,N,PATH2,HFORM2,17)
      CALL OPNOD(NLAYR,N,PATH3,HFORM3,18)
      CALL BUFIN(16,IEOF,XFILH1(1),NFHDRF)
      CALL BUFIN(17,IEOF,XFILH2(1),NFHDRF)
      CALL BUFIN(18,IEOF,XFILH3(1),NFHDRF)
C
C     Set DV = DV1
C
      DV = DV1
C
      NDONE = 0
      NPAN = 1
 110  CONTINUE
      CALL BUFIN(16,IEOF,PNLHD(1),NPHDRF)
      V1POD(NPAN) = V1P
      V2POD(NPAN) = V2P
      DVPOD(NPAN) = DVP
      NLIMOD(NPAN) = NLIMO
      IF (IEOF.LE.0) GOTO 200
      CALL BUFIN(17,IEOF,PNLHD(1),NPHDRF)
      CALL BUFIN(18,IEOF,PNLHD(1),NPHDRF)
      CALL BUFIN(16,IEOF,XOPTDP(1+NDONE,1),NLIMO)
      CALL BUFIN(17,IEOF,XOPTDP(1+NDONE,2),NLIMO)
      CALL BUFIN(18,IEOF,XOPTDP(1+NDONE,3),NLIMO)
      NDONE = NDONE+NLIMO
      NPAN = NPAN+1
      GOTO 110
 200  CONTINUE
C
C     Restore value of iemit
C
      IEMIT = IEMSV
      RETURN
C
      END
C
C
C     ---------------------------------------------------------------
C
      SUBROUTINE OPNOD(NLAYER,LAYER,PTHODL,HFMODL,IFILE)
C
C     This subroutine opens file for calculating the radiance using
C     PRECALCULATED ABSORPTANCE COEFFICICIENTS
C
      LOGICAL OP
      CHARACTER*57 FILE1
      CHARACTER*55 PTHODL
      CHARACTER*11 CFORM
      CHARACTER*10 HFMODL
C
C     -------------------------
C     Common block for analytic derivative
C     -------------------------
      COMMON /ADRFIL/ KODFIL,KODTOT,KTEMP,KFILAD
C     -------------------------
C
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4
C
C           123456789-123456789-123456789-123456789-123456789-1234567
      DATA FILE1 /
     *     '                                                         '/
      DATA CFORM / 'UNFORMATTED' /
C
C
C     For calculation of analytic derivative, open previously
C     calculated optical depth files as unit KODFIL.  Otherwise,
C     open as KFILE.
C
C      WRITE(*,910) LAYER,NLAYER
      INQUIRE (UNIT=IFILE,OPENED=OP)
      IF (OP) CLOSE (IFILE)
      WRITE(FILE1,HFMODL) PTHODL,LAYER
      OPEN(UNIT=IFILE,FILE=FILE1,FORM=CFORM,STATUS='OLD')
C
C     Write procedure
C
C      WRITE(*,900) FILE1
C
      RETURN
C
 900  FORMAT ('          OPENED LAYER ABSORPTANCE COEF. FILE:  ',A57)
 910  FORMAT ('LAYER ',I5,' OF ',I5,':')
C
      END
C
C     ----------------------------------------------------------------

