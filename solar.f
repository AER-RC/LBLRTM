C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
C
C     ----------------------------------------------------------------
C
      SUBROUTINE SOLINT(IFILE,LFILE,NPTS,INFLAG,IOTFLG)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      
C               LAST MODIFICATION:    3 April 1994                   
C                                                                      
C                  IMPLEMENTATION:    P.D. Brown
C                                                                      
C             ALGORITHM REVISIONS:    S.A. Clough
C                                     P.D. Brown
C                                                                      
C                                                                      
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.      
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139       
C                                                                      
C----------------------------------------------------------------------
C                                                                      
C               WORK SUPPORTED BY:    THE ARM PROGRAM                  
C                                     OFFICE OF ENERGY RESEARCH        
C                                     DEPARTMENT OF ENERGY             
C                                                                      
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT DOUBLE PRECISION (V)
C
C     ------------------------------------------------------------
C     SUBROUTINE SOLINT interpolates solar radiances from the binary
C     file SOLAR.RAD.  The following are input and output options:
C
C       INFLAG = 0   => input transmittance from TAPE12 (default).
C              = 1   => input optical depths from TAPE10 and
C                       convert to transmittance.
C
C       IOTFLG = 0   => attenuate w/transmittance & output (default).
C              = 1   => attenuate and add to radiance from TAPE12
C                       (requires INFLAG = 1).
C
C     Output radiance goes to TAPE11.
C
C     ------------------------------------------------------------
C
      COMMON /MANE/ P0,TEMP0,NLAYER,DDUM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR,
     *              AVFIX,LAYDUM,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,
     *              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,
     *              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,
     *              EXTID(10)
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2
C
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID
C
      COMMON /EMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *               WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,
     *               EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
      COMMON /OPANL/ V1PO,V2PO,DVPO,NLIMO
      COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTSS
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4
      COMMON /ARMCM1/ HVRSOL
C
      DIMENSION XFILHD(2),PNLHDR(2),OPNLHD(2)
      DIMENSION A1(0:100),A2(0:100),A3(0:100),A4(0:100)
      DIMENSION TRAO(2),TRAN(2410)
      DIMENSION RADO(2),RADN(2410)
      DIMENSION OPTO(2),OPTN(2410)
C
      DIMENSION SOLAR(-1:4818)
      DIMENSION SOLRAD(2410)
C
      CHARACTER*40 CYID
      CHARACTER*8 HVRSOL
C
      EQUIVALENCE (XFILHD(1),XID(1)),(PNLHDR(1),V1P),
     *            (OPNLHD(1),V1PO)
      EQUIVALENCE (TRAN(1),TRAO(1)),(RADN(1),RADO(1)),
     *            (OPTN(1),OPTO(1)),
     *            (FSCDID(4),IAERSL),(FSCDID(5),IEMIT),
     *            (FSCDID(7),IPLOT),(FSCDID(8),IPATHL),
     *            (FSCDID(16),LAYR1)
C
C     ************************************************************
C     ****** THIS PROGRAM DOES MERGE FOR SOLAR RADIANCE AND ******
C     ****** TRANMITTANCE USING FOUR POINT INTERPOLATION    ******
C     ************************************************************
C
C
C     ASSIGN SCCS VERSION NUMBER TO MODULE 
C
      HVRSOL = '$Revision$'
C
C     -------------------
C
C     Open file SOLAR.RAD
C
      ISOLFL = 19
      OPEN(UNIT=ISOLFL,FILE='SOLAR.RAD',FORM='UNFORMATTED',
     *     STATUS='OLD')
C
      CALL CPUTIM (TIME)
      WRITE (IPR,900) TIME
      NPANLS = 0
      TIMRD = 0.0
      TIMOT = 0.0
C
C     FOR AEROSOL RUNS, MOVE YID (IFILE) INTO YID (LFILE)
C
C     Read file header of solar radiance file and determine dv ratio
C
      IF (IAERSL.GT.0) WRITE (CYID,'(5A8)') (YID(I),I=3,7)
      CALL BUFIN (ISOLFL,KEOF,XFILHD(1),NFHDRF)
cpdb      IF (IAERSL.GT.0) READ (CYID,'(5A8)') (YID(I),I=3,7)
cpdb      IF (JPATHL.GE.1) IPATHL = JPATHL
      DVK = DV
C
C     Read in file header of transmittance/optical depth file
C
      CALL BUFIN (IFILE,LEOF,XFILHD(1),NFHDRF)
      DVL = DV
C
      ATYPE = 9.999E09
      IF (DVK.EQ.DVL) ATYPE = 0.
      IF (DVL.GT.DVK) ATYPE = DVK/(DVL-DVK)+0.5
      IF (DVL.LT.DVK) ATYPE = -DVL/(DVK-DVL)-0.5
C
C     IF (ATYPE .GT. 0) STOP  ' SOLINT; ATYPE GT 0 '
C
C
C     Write file information out to TAPE6
C
      WRITE (IPR,905) ISOLFL,IFILE,LFILE,ATYPE,INFLAG,IOTFLG
      IF (INFLAG.EQ.0) THEN
         WRITE(IPR,920) IFILE
      ELSE
         WRITE(IPR,925) IFILE
      ENDIF
C
C     Test for INFLAG=0, that atmospheric radiance is to be read
C
      IF (IOTFLG.EQ.0) THEN
         WRITE(IPR,930) LFILE
      ELSE
         IF (INFLAG.EQ.1) STOP 'ERROR: INFLAG=1, IOTFLG=1'
         WRITE(IPR,935) LFILE
      ENDIF
      IEMIT = 1
      SECANT = 0.
      DV = DVL
C
C     Output file header
C
      CALL BUFOUT (LFILE,XFILHD(1),NFHDRF)
C
      IF (ATYPE.EQ.0.) THEN
C
C        1/1 ratio only
C
   30    CONTINUE
         CALL CPUTIM (TIMSL1)
         CALL SOLIN (V1P,V2P,DVP,NLIM,ISOLFL,SOLAR(1),LSEOF,NPANLS)
         CALL CPUTIM (TIMSL2)
         TIMRD = TIMRD+TIMSL2-TIMSL1
         IF (LSEOF.LE.0) GO TO 110
C
C        Read file header from transmittance/optical depth file
C
         CALL BUFIN (IFILE,LEOF,OPNLHD(1),NPHDRF)
C
C        If INFLAG = 0, then read radiance and tranmittance
C        If INFLAG = 1, then read optical depth
C
         IF (INFLAG.EQ.0) THEN
            CALL BUFIN (IFILE,LEOF,RADO(1),NLIMO)
            CALL BUFIN (IFILE,LEOF,TRAO(1),NLIMO)
         ELSE
            CALL BUFIN (IFILE,LEOF,OPTO(1),NLIMO)
            DO 35 I = 1,NLIMO
               TRAN(I) = EXP(-OPTN(I))
 35         CONTINUE
         ENDIF
         CALL CPUTIM (TIMSL3)
         TIMRD = TIMRD+TIMSL3-TIMSL2
C
C        If IOTFLG = 0, then calculate attenuated solar radiance
C        If IOTFLG = 1, then calculate attenuated solar radiance
C                       plus atmospheric radiance
C
         IF (IOTFLG.EQ.0) THEN
            DO 40 I = 1, NLIM
               SOLRAD(I) = SOLAR(I)*TRAN(I)
 40         CONTINUE
         ELSE
            DO 41 I = 1, NLIM
               SOLRAD(I) = SOLAR(I)*TRAN(I)+RADN(I)
 41         CONTINUE
         ENDIF
C
         CALL CPUTIM (TIMSL2)
         CALL SOLOUT(V1PO,V2PO,DVPO,NLIMO,SOLRAD,LFILE,NPTS,NPANLS)
         CALL CPUTIM (TIMSL3)
         TIMOT = TIMOT+TIMSL3-TIMSL2
         GO TO 30
C
      ENDIF
C
C     All ratios except 1/1
C
      DO 50 JP = 0,100
         APG = JP
         P = 0.01*APG
C
C        The following are the constants for the Lagrange
C        4 point interpolation
C
         A1(JP) = -P*(P-1.0)*(P-2.0)/6.0
         A2(JP) = (P**2-1.0)*(P-2.0)*0.5
         A3(JP) = -P*(P+1.0)*(P-2.0)*0.5
         A4(JP) = P*(P**2-1.0)/6.0
   50 CONTINUE
C
C     *** Beginning of loop that does merge  ***
C
      NPE = 0
      SOLAR(0) = 0.0
      V1P = 0.0
      V2P = 0.0
      DVP = 0.0
      V1PO = 0.0
      V2PO = 0.0
      DVPO = 0.0
      LSEOF = 1
C
      ip = 0
   60 CONTINUE
      ip = ip+1
C
C     Read file header from transmittance/optical depth file
C
      CALL CPUTIM(TIMSL1)
      CALL BUFIN(IFILE,LEOF,OPNLHD(1),NPHDRF)
      CALL CPUTIM(TIMSL2)
      TIMRD = TIMRD+TIMSL2-TIMSL1
      IF (LEOF.LE.0) GO TO 110
C
C     If INFLAG = 0, then read radiance and tranmittance
C     If INFLAG = 1, then read optical depth
C
      IF (INFLAG.EQ.0) THEN
         CALL BUFIN (IFILE,LEOF,RADO(1),NLIMO)
         CALL BUFIN (IFILE,LEOF,TRAO(1),NLIMO)
      ELSE
         CALL BUFIN (IFILE,LEOF,OPTO(1),NLIMO)
         DO 65 I = 1,NLIMO
            TRAN(I) = EXP(-OPTN(I))
 65      CONTINUE
      ENDIF
      CALL CPUTIM(TIMSL3)
      TIMRD = TIMRD+TIMSL3-TIMSL2
      II = 1
C
C     Buffer in panels from solar radiance file
C
      IF (V2P.LE.V2PO+DVP .AND.LSEOF.GT.0) THEN
   70    CALL CPUTIM(TIMSL2)
         CALL SOLIN(V1P,V2P,DVP,NLIM,ISOLFL,SOLAR(NPE+1),LSEOF,NPANLS)
         CALL CPUTIM(TIMSL3)
         TIMRD = TIMRD+TIMSL3-TIMSL2
         IF (LSEOF.LE.0) GO TO 80
         V1P = V1P-FLOAT(NPE)*DVP
         NPE = NLIM+NPE
         IF (V2P.LE.V2PO+DVP) GO TO 70
      ENDIF
C
C     Zero point of first panel
C
   80 IF (SOLAR(0).EQ.0.0) THEN
         SOLAR(-1) = SOLAR(1)
         SOLAR(0) = SOLAR(1)
      ENDIF
C
C     End point of last panel
C
      IF (V2P+DVP.GE.V2) THEN
         SOLAR(NPE+1) = SOLAR(NPE)
         SOLAR(NPE+2) = SOLAR(NPE)
      ENDIF
C
C     NPL is the location of first element on arrays RADO and TRAO
C
      NPL = 1
C
      RATDV = DVL/DVK
C
C     FJJ is offset by 2. (for rounding purposes)
C
      FJ1DIF = (V1PO-V1P)/DVP+1.+2.
C
C     ***** Beginning of loop that does merge  *****
C
C     If IOTFLG = 0, then calculate attenuated solar radiance
C     If IOTFLG = 1, then calculate attenuated solar radiance
C                    plus atmospheric radiance
C
      IF (IOTFLG.EQ.0) THEN
         DO 90 II = 1, NLIMO
            FJJ = FJ1DIF+RATDV*FLOAT(II-1)
            JJ = IFIX(FJJ)-2
            JP = (FJJ-FLOAT(JJ))*100.-199.5
            SOLRAD(II) = (A1(JP)*SOLAR(JJ-1)+A2(JP)*SOLAR(JJ)+
     *           A3(JP)*SOLAR(JJ+1)+A4(JP)*SOLAR(JJ+2))*TRAN(II)
c
 90      CONTINUE
      ELSE
         DO 91 II = 1, NLIMO
            FJJ = FJ1DIF+RATDV*FLOAT(II-1)
            JJ = IFIX(FJJ)-2
            JP = (FJJ-FLOAT(JJ))*100.-199.5
c            write(*,*) 'TAPE11:'
c            write(*,*) '   ',v1po+dvpo*(ii),rado(ii),trao(ii)
c            write(*,*) 'SOLAR:'
c            write(*,*) '   ',v1p+dvp*(jj-2),solar(jj-1)
c            write(*,*) '   ',v1p+dvp*(jj-1),solar(jj)
c            write(*,*) '   ',v1p+dvp*(jj),solar(jj+1)
c            write(*,*) '   ',v1p+dvp*(jj+1),solar(jj+2)
c
            SOLRAD(II) = (A1(JP)*SOLAR(JJ-1)+A2(JP)*SOLAR(JJ)+
     *           A3(JP)*SOLAR(JJ+1)+A4(JP)*SOLAR(JJ+2))*TRAN(II)+
     *           RADN(II)
c
c            write(*,*) 'SOLRAD:'
c            write(*,*) '   ',v1po+dvpo*(ii),solrad(ii)
c            write(*,*) ' ---'
 91      CONTINUE
      ENDIF
C
      NPL = JJ-1
C
      CALL CPUTIM (TIMSL1)
C
C     Output attenuated radiance
C
      CALL SOLOUT(V1PO,V2PO,DVPO,NLIMO,SOLRAD,LFILE,NPTS,NPANLS)
      CALL CPUTIM (TIMSL2)
      TIMOT = TIMOT+TIMSL2-TIMSL1
C
C     NPL is now location of first element to be used for next pass
C
      IPL = -2
      DO 100 NL = NPL, NPE
         IPL = IPL+1
         SOLAR(IPL) = SOLAR(NL)
  100 CONTINUE
C
      V1P = V1P+FLOAT(NPL+1)*DVP
      NPE = IPL
C
      GO TO 60
  110 CONTINUE
C
      CALL CPUTIM (TIME1)
      TIM = TIME1-TIME
      WRITE (IPR,910) TIME1,TIM,TIMRD,TIMOT
C
  900 FORMAT ('0 THE TIME AT THE START OF SOLINT IS ',F12.3)
  905 FORMAT ('0 FILE ',I5,' MERGED WITH FILE ',I5,' ONTO FILE',
     *        I5,'  WITH XTYPE=',G15.5,/,'0 INFLAG = ',I5,4X,
     *        'IOTFLG = ',I5)
  910 FORMAT ('0 THE TIME AT THE END OF SOLINT IS ',F12.3/F12.3,
     *        ' SECS WERE REQUIRED FOR THIS SOLAR MERGE',F12.3,
     *        ' - READ - ',F12.3,' - SOLOUT - ',F12.3)
 920  FORMAT ('0 Radiance and Transmittance read in from unit',I5)
 925  FORMAT ('0 Optical Depths read in from unit',I5)
 930  FORMAT ('0 Attenuated solar radiance output to unit',I5,/)
 935  FORMAT ('0 Attenuated solar radiance + atmospheric radiance',
     *        1x,'output to unit',I5,/)
C
      END
C
C     ----------------------------------------------------------------
C
      SUBROUTINE SOLIN (V1P,V2P,DVP,NLIM,KFILE,SOLAR,KEOF)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      
C               LAST MODIFICATION:    3 April 1994                   
C                                                                      
C                  IMPLEMENTATION:    P.D. Brown
C                                                                      
C             ALGORITHM REVISIONS:    S.A. Clough
C                                     P.D. Brown
C                                                                      
C                                                                      
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.      
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139       
C                                                                      
C----------------------------------------------------------------------
C                                                                      
C               WORK SUPPORTED BY:    THE ARM PROGRAM                  
C                                     OFFICE OF ENERGY RESEARCH        
C                                     DEPARTMENT OF ENERGY             
C                                                                      
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT DOUBLE PRECISION (V)
C
C     SUBROUTINE SOLIN inputs solar radiation from the file "SOLAR.RAD"
C     for interpolation in SOLINT.
C
      DOUBLE PRECISION XID,SECANT,HMOLID,XALTZ,YID
C
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,
     *                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4
      COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF
C
      DIMENSION PNLHDR(2),SOLAR(*)
C
      EQUIVALENCE (PNLHDR(1),V1PBF)
C
      CALL BUFIN (KFILE,KEOF,PNLHDR(1),NPHDRF)
      IF (KEOF.LE.0) RETURN
      CALL BUFIN (KFILE,KEOF,SOLAR(1),NLIMBF)
C
      V1P = V1PBF
      V2P = V2PBF
      DVP = DVPBF
      NLIM = NLIMBF
C
      RETURN
C
      END
C
C     ----------------------------------------------------------------
C
      SUBROUTINE SOLOUT (V1P,V2P,DVP,NLIM,SOLRAD,LFILE,NPTS,NPANLS)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      
C               LAST MODIFICATION:    3 April 1994                   
C                                                                      
C                  IMPLEMENTATION:    P.D. Brown
C                                                                      
C             ALGORITHM REVISIONS:    S.A. Clough
C                                     P.D. Brown
C                                                                      
C                                                                      
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.      
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139       
C                                                                      
C----------------------------------------------------------------------
C                                                                      
C               WORK SUPPORTED BY:    THE ARM PROGRAM                  
C                                     OFFICE OF ENERGY RESEARCH        
C                                     DEPARTMENT OF ENERGY             
C                                                                      
C                                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT DOUBLE PRECISION (V)
C
C     SUBROUTINE SOLOUT OUTPUTS ATTENUATED SOLAR RADIANCE (INTERPOLATED)
C     TO LFILE
C
      COMMON /BUFPNL/ V1PBF,V2PBF,DVPBF,NLIMBF
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,
     *              NLTEFL,LNFIL4,LNGTH4
      DIMENSION PNLHDR(2)
      DIMENSION SOLRAD(*)
C
      EQUIVALENCE (PNLHDR(1),V1PBF)
C
      REAL SOLRAD
C
      NPANLS = NPANLS+1
      V1PBF = V1P
      V2PBF = V2P
      DVPBF = DVP
      NLIMBF = NLIM
C
      CALL BUFOUT (LFILE,PNLHDR(1),NPHDRF)
      CALL BUFOUT (LFILE,SOLRAD(1),NLIMBF)
C
      IF (NPTS.GT.0) THEN
         IF (NPANLS.EQ.1) WRITE (IPR,900)
         WRITE (IPR,905)
         NNPTS = NPTS
         IF (NPTS.GT.(NLIMBF/2)+1) NNPTS = (NLIMBF/2)+1
         JEND = NLIMBF-NNPTS+1
         DO 10 J = 1, NNPTS
            VJ = V1PBF+FLOAT(J-1)*DVPBF
            K = J+JEND-1
            VK = V1PBF+FLOAT(K-1)*DVPBF
            WRITE (IPR,910) J,VJ,SOLRAD(J),K,VK,SOLRAD(K)
   10    CONTINUE
      ENDIF
C
      RETURN
C
  900 FORMAT ('0 ','LOCATION  WAVENUMBER',4X,'RADIANCE',13X,
     *        'LOCATION   WAVENUMBER',4X,
     *        'RADIANCE')
  905 FORMAT (' ')
  910 FORMAT (I8,2X,F12.6,1P,E15.7,0P,9X,I8,2X,F12.6,1P,E15.7,0P)
C
      END

