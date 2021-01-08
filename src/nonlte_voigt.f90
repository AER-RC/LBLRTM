!     path:      $HeadURL$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!
!  --------------------------------------------------------------------------
! |  Copyright �, Atmospheric and Environmental Research, Inc., 2015         |
! |                                                                          |
! |  All rights reserved. This source code is part of the LBLRTM software    |
! |  and is designed for scientific and research purposes. Atmospheric and   |
! |  Environmental Research, Inc. (AER) grants USER the right to download,   |
! |  install, use and copy this software for scientific and research         |
! |  purposes only. This software may be redistributed as long as this       |
! |  copyright notice is reproduced on any copy made and appropriate         |
! |  acknowledgment is given to AER. This software or any modified version   |
! |  of this software may not be incorporated into proprietary software or   |
! |  commercial software offered for sale without the express written        |
! |  consent of AER.                                                         |
! |                                                                          |
! |  This software is provided as is without any express or implied          |
! |  warranties.                                                             |
! |                       (http://www.rtweb.aer.com/)                        |
!  --------------------------------------------------------------------------
!
SUBROUTINE NONLTE(MPTS)
!
   USE lblparams
!      include 'lblparams.inc'
   IMPLICIT REAL*8           (V)
!
!**********************************************************************
!*
!*
!*    CALCULATES MONOCHROMATIC ABSORPTION COEFFICIENT FOR SINGLE LAYER
!*    Under the conditions of NLTE
!*
!**********************************************************************
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!                  IMPLEMENTATION:    W.O. GALLERY
!
!             ALGORITHM REVISIONS:    M.W.SHEPHARD
!                                     D. WEISENSTEIN
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     131 Hartwell Ave,  Lexington,  MA   02421
!
!----------------------------------------------------------------------
!
!               WORK SUPPORTED BY:    JPL, TES
!                                     JCSDA
!
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL FASCOD
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID

   CHARACTER*8      XID,       HMOLID,      YID
   Real*8               SECANT,       XALTZ
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /VBNLTE/ RATSTATE(MAXSTATE*Max_ISO,MXMOL),NUMSTATE(MXMOL)
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /R4SUB/ VLOF4,VHIF4,ILOF4,IST,IHIF4,LIMIN4,LIMOUT,ILAST,   &
   &               DPTMN4,DPTFC4,ILIN4,ILIN4T
   COMMON /ISVECT/ ISO_MAX(MXMOL),SMASSI(mxmol,10)
   common /cmol_nam/ cmol(mxmol),cspc(mxspc)
   CHARACTER*6  CMOL,CSPC
!
   CHARACTER*18 hnmnlte,hvnlte
   COMMON /CVNLTE/ HNMNLTE,HVNLTE

   EQUIVALENCE (FSCDID(1),IHIRAC),(FSCDID(2),ILBLF4),                &
   &  (FSCDID(3),IXSCNT),(FSCDID(4),IAERSL),(FSCDID(5),IEMIT),        &
   &  (FSCDID(7),IPLOT), (FSCDID(8),IPATHL),(FSCDID(9),JRAD),         &
   &  (FSCDID(11),IMRG)
!
   CHARACTER*5 IDSTATE
   INTEGER INDEX(MXMOL)
   INTEGER ISORDER(MAXSTATE*Max_ISO)
   DIMENSION IDSTATE(MAXSTATE,MXMOL),                                &
   &     EESTATE(MAXSTATE*Max_ISO,MXMOL),                             &
   &     NDGSTATE(MAXSTATE*Max_ISO,MXMOL)
   CHARACTER*80 TIT(3),TEXTLINE
   LOGICAL ISODATA
   CHARACTER*6  TXTISO,BLANKS
   DATA BLANKS/'      '/
!
   CALL DEFNLTEDAT(NUMSTATE,IDSTATE,EESTATE,NDGSTATE,RATSTATE)
!
   ISODATA=.FALSE.
!
   REWIND NLTEFL

! READ UP TO 20 LINES OF TEXT AT BEGINNING OF TAPE4
   WRITE(IPR,890)
890 FORMAT(/2X,'TAPE4 HEADER:')
   ! 20 MAX TEXT LINES AT BEGINNING OF TAPE4
   DO I=1,20
      READ(NLTEFL,900) TEXTLINE
900   FORMAT(A80)
      IF(TEXTLINE(1:1).NE.'!') GO TO 915
      WRITE(IPR,910) TEXTLINE
910   FORMAT(2X,A80)
   END DO
915 READ(TEXTLINE,920) IVIB,MOLNEQ
920 FORMAT(2I5)
   WRITE(IPR,921) IVIB,ALTAV,TAVE
921 FORMAT(/' IVIB =',I5,/'  ALT = ',F10.3,'  TEMP =',F10.3)

   READ(NLTEFL,900) TEXTLINE
   write(ipr,940) textline
940 FORMAT(A80)
   DO I=1,74
      IF(TEXTLINE(I:I+6).EQ.'VIBRATI') ISODATA=.TRUE.
   END DO
   IF(ISODATA) THEN
10    CALL GETINDEX(TEXTLINE,CMOL,MXMOL,ID,TXTISO)
!           END OF DATA ENCOUNTERED
      IF(ID.EQ.0) GO TO 30
      CALL RDNLTE(NLTEFL,TEXTLINE,TXTISO,NUMSTATE(ID), IDSTATE(1,ID),&
         EESTATE(1,ID),NDGSTATE(1,ID), ISORDER)
      IF(IVIB.EQ.1) THEN
         CALL VIBTMP(XKT,ALTAV,NLTEFL,NUMSTATE(ID), IDSTATE(1,ID),   &
            NDGSTATE(1,ID),EESTATE(1,ID), RATSTATE(1,ID),TXTISO,        &
            TEXTLINE,ISORDER)
      ELSE
         CALL VIBPOP(XKT,ALTAV,NLTEFL,NUMSTATE(ID), IDSTATE(1,ID),   &
            NDGSTATE(1,ID),EESTATE(1,ID), RATSTATE(1,ID),TXTISO,        &
            TEXTLINE,ISORDER)
      END IF
!                 IF END OF DATA ENCOUNTERED
      DO I=1,70
         IF(TEXTLINE(I:I+10).EQ.'END OF DATA') GO TO 30
      END DO
!              READ DATA FOR NEXT SPECIE
      GO TO 10
   ELSE
      DO ID=1,MXMOL
         IF(NUMSTATE(ID).GT.0) THEN
            CALL DROPSPACE(CMOL(ID),TXTISO)
            IF(IVIB.EQ.1) THEN
               CALL VIBTMP(XKT,ALTAV,NLTEFL,NUMSTATE(ID), IDSTATE(1,ID),&
                  NDGSTATE(1,ID),EESTATE(1,ID), RATSTATE(1,ID),TXTISO,     &
                  TEXTLINE,ISORDER)
            ELSE
               CALL VIBPOP(XKT,ALTAV,NLTEFL,NUMSTATE(ID), IDSTATE(1,ID),&
                  NDGSTATE(1,ID),EESTATE(1,ID), RATSTATE(1,ID),TXTISO,     &
                  TEXTLINE,ISORDER)
            END IF
         END IF
      END DO
   END IF
!
30 IPFLAG = 0
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

!*****Note: this check not needed for lblrtm
! CHECK TO SEE IF VLOF4 OR VHIF4 FALLS IN LINE COUPLING REGION
!      CALL CHKLNC(VLOF4,VHIF4,(SAMPLE*DV),0)

   IF(ILBLF4.GE.1) CALL LINF4Q(VLOF4,VHIF4)
   CALL HIRACQ(MPTS)
   WRITE(IPR,930)
930 FORMAT(//)
   RETURN
end subroutine NONLTE
!
! --------------------------
SUBROUTINE GETINDEX(TEXTLINE,CMOL,MXMOL,ID,TXTISO)
! ------- GET MOLECULE INDEX FROM CMOL ARRAY FOR VIBRATIONAL DATA
   CHARACTER*80 TEXTLINE
   CHARACTER*6 CMOL(MXMOL),TXTISO,TXTMOL,BLANKS
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFLE,LNFIL4,LNGTH4
   DATA BLANKS/'      '/

   I1=0
   I2=0
   DO I=1,80
      IF(TEXTLINE(I:I).NE.'-' .AND. TEXTLINE(I:I).NE.' ') THEN
         IF(I1.EQ.0) I1=I
      END IF
      IF(I1.GT.0 .AND. TEXTLINE(I:I).EQ.' ') THEN
         IF(I2.EQ.0) I2=I-1
      END IF
   END DO
   TXTISO=TEXTLINE(I1:I2)//BLANKS
   ID=0
   DO ISP=1,MXMOL
      CALL DROPSPACE(CMOL(ISP),TXTMOL)
      IF(TXTISO.EQ.TXTMOL) THEN
         ID=ISP
         WRITE(IPR,*) 'GETINDEX: ',TXTISO,'  INDEX=',ID
         RETURN
      END IF
   END DO
   WRITE(IPR,*) 'READING NONLTE DATA HEADER FROM TAPE4'
   WRITE(IPR,*) 'BUT THIS SPECIE ',TXTISO,                           &
   &     ' NOT FOUND IN MOLECULE LIST'
!##     STOP 'ERROR READING NLTE DATA FROM TAPE4'
   RETURN
end subroutine GETINDEX
! ----------------------------------------------------------------
SUBROUTINE RDNLTE(NLTEFL,TEXTLINE,TXTISO,IMAX,                    &
&     IDX,EEX,NDG,ISORDER)
!
   USE lblparams
!      include 'lblparams.inc'
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFLE,LNFIL4,LNGTH4
   DIMENSION IDX(MAXSTATE),ISORDER(MAXSTATE*Max_ISO)
   DIMENSION EEX(MAXSTATE*Max_ISO),NDG(MAXSTATE*Max_ISO)
   CHARACTER*80 TEXTLINE
   CHARACTER*6 TXTISO
   CHARACTER*5 IDX,IDENT,BLANKS
   CHARACTER*1 QUOTE
   DATA BLANKS/' '/
   ! QUOTE = '
   QUOTE=CHAR(39)
!       INITIALIZE ARRAY TO HOLD ISOTOPE ORDER INFO
   DO I=1,MAXSTATE*Max_ISO
      ISORDER(I)=0
   END DO

   INUM=0
50 READ(NLTEFL,940) TEXTLINE
940 FORMAT(A80)
   IF(TEXTLINE(1:2).EQ.'--') RETURN
   INUM=INUM+1
   I1=0
   I2=0
   DO I=1,80
      IF(TEXTLINE(I:I).EQ.QUOTE) THEN
         IF(I1.EQ.0) THEN
            I1=I
         ELSE IF(I2.EQ.0) THEN
            I2=I
         END IF
      END IF
   END DO
   READ(TEXTLINE(1:I1-1),*) NN,ISOTOPE
   IF(NN.NE.INUM .AND. ISOTOPE.EQ.1) THEN
      WRITE(IPR,*) 'WARNING IN TAPE 4: ISOTOPE DATA FOR ',TXTISO
      WRITE(IPR,*) 'LINE ',INUM,' NOT IN ORDER'
!         STOP 'NN.NE.INUM IN NONLTE DATA FROM TAPE4'
   END IF
   IF(ISOTOPE.EQ.1 .AND. INUM.GT.IMAX) THEN
      WRITE(IPR,*) 'READING NONLTE DATA FROM TAPE4'
      WRITE(IPR,*) 'ISOTOPE DATA FOR ',TXTISO,                       &
      &        ' CONTAINS MORE STATES THAN DEFAULT ALLOWANCE'
      STOP 'NUMBER OF NONLTE STATES IN TAPE4 EXCEEDS DEFAULT'
   END IF
   IF(ISOTOPE.GT.Max_ISO) THEN
      WRITE(IPR,*) 'ERROR IN TAPE4:  ISOTOPE NUMBER ',ISOTOPE,       &
      &        ' GREATER THAN Max_ISO=',Max_ISO
      STOP 'TAPE4 ERROR: ISOTOPE NUMBER TOO LARGE'
   END IF
   IDENT=TEXTLINE(I1+1:I2-1)//BLANKS
   INDEX=0
   DO I=1,MAXSTATE
      IF(IDENT.EQ.IDX(I)) INDEX=I
   END DO
   IF(INDEX.EQ.0) THEN
      WRITE(IPR,*) 'ERROR IN TAPE 4 ISOTOPE DATA FOR ',TXTISO,       &
         ' ON LINE ',INUM
      WRITE(IPR,*) '   ISOTOPE NUMBER ',ISOTOPE,' IDENTIFIER ',IDENT,&
      &       ' NOT CONSISTENT WITH EXPECTED INPUT'
      STOP 'ERROR IN ISOTOPE INPUT IN NONLTE DATA FROM TAPE4'
   END IF
   INDEX2=(ISOTOPE-1)*MAXSTATE + INDEX
   ISORDER(INUM)=INDEX2
   READ(TEXTLINE(I2+1:80),*) EEX(INDEX2),NDG(INDEX2)
   IF(INUM.EQ.1 .AND. ABS(EEX(INDEX2)).GE.1.E-25) THEN
      WRITE(IPR,*) 'RDNLTE: GROUND STATE NOT SPECIFIED ',            &
      &        'FOR MOLECULE ',TXTISO
      STOP 'TAPE4 GROUND STATE MISSING'
   END IF
!      WRITE(IPR,*) 'RDNLTE, LINE ',inum,' state=',index,'  index=',
!     $     index2,'  isorder=',isorder(inum)
   GO TO 50
end subroutine RDNLTE
!
! ----------------------------------------------------------------
!
SUBROUTINE VIBPOP(XKT,HT,NLTEFLAG,NUM,IDX,NDEG,EH,RAT,            &
&   TITMOL,TEXTLINE,ISORDER)
!
!
!     SUBROUTINE VIBPOP USES THE NON-LTE POPULATION DATA FROM
!     TAPE4 TO CALCULATE THE VIBRATIONAL POPULATION ENHANCEMENT
!     RATIOS FOR SELECTED VIBRATIONAL STATES OF H2O,CO2,NO AND O3.
!
   USE phys_consts, ONLY: radcn2
   USE lblparams
!      include 'lblparams.inc'
   IMPLICIT REAL*8           (V)

   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4

   CHARACTER*5 IDX
   CHARACTER*80 TEXTLINE
   CHARACTER*8      XID,       HMOLID,      YID
   Real*8               SECANT,       XALTZ
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   DIMENSION IDX(MAXSTATE),TNE(MAXSTATE*Max_ISO)
   DIMENSION VQNE(MAXSTATE*Max_ISO),VQEQ(MAXSTATE*Max_ISO)
   DIMENSION VPNE1(MAXSTATE*Max_ISO),VPNE2(MAXSTATE*Max_ISO),        &
   &     VQNEST(MAXSTATE*Max_ISO),VQNEIN(MAXSTATE*Max_ISO)
   DIMENSION NDEG(MAXSTATE*Max_ISO),EH(MAXSTATE*Max_ISO),            &
   &     RAT(MAXSTATE*Max_ISO),ISORDER(MAXSTATE*Max_ISO)
   CHARACTER*6 TITMOL,HMNLTE,BLANKS
   DATA BLANKS/' '/
!
!   READ NLTE VIB POPULATIONS
!
   XKT= TAVE/RADCN2
   DO I=1,MAXSTATE*Max_ISO
      IF(ISORDER(I).GT.0) NUMIN=I
   end do
!      write(ipr,*) 'called vibpop with txtmol=',titmol,'  alt=',ht
!      write(ipr,*) 'num=',num,'  numin=',numin
!      write(ipr,*) (isorder(i),i=1,numin)
!
   I1=0
   I2=0
   DO I=1,80
      IF(TEXTLINE(I:I).NE.'-' .AND. TEXTLINE(I:I).NE.' ') THEN
         IF(I1.EQ.0) I1=I
      END IF
      IF(I1.GT.0 .AND. TEXTLINE(I:I).EQ.' ') THEN
         IF(I2.EQ.0) I2=I-1
      END IF
   END DO
   HMNLTE=TEXTLINE(I1:I2)//BLANKS
   IF(HMNLTE.NE.TITMOL) THEN
      READ (NLTEFLAG,902) HMNLTE
      i1=0
      DO I=1,6
         IF(HMNLTE(I:I).EQ.' ') I1=I
      END DO
      HMNLTE=HMNLTE(I1+1:6)//BLANKS
   END IF
   IF(HMNLTE.NE.TITMOL) THEN
      WRITE(IPR,*) 'READING VIBPOP DATA FROM TAPE4'
      WRITE(IPR,*) 'EXPECTED PROFILE DATA FOR ',TOTMOL
      WRITE(IPR,*) 'BUT READ SPECIE ',HMNLTE
      STOP 'ERROR READING NLTE DATA FROM TAPE4'
   END IF
   READ (NLTEFLAG,*)  ALT1,(VPNE1(I),I=1,NUMIN)
10 READ (NLTEFLAG,*)  ALT2,(VPNE2(I),I=1,NUMIN)
!*****WOG, 11/06/2000: ALT1 -> AL2:
!     IF( ALT1.LT.HT) THEN
   IF( ALT2.LE.HT) THEN
      ALT1 = ALT2
      DO 20 I=1,NUMIN
         VPNE1(I) = VPNE2(I)
20    CONTINUE
      GO TO 10
   ENDIF
!     CALL LININT(HT,ALT1,ALT2,NUMIN,VPNE1,VPNE2,VQNEIN)
   A = (HT-ALT1)/(ALT2-ALT1)
   DO 30 I=1,NUMIN
      CALL EXPINT(VQNEIN(I),VPNE1(I),VPNE2(I),A )
30 END DO
   DO I=1,NUMIN
      INDEX=ISORDER(I)
      VQNE(INDEX)=VQNEIN(I)
   END DO
!
   POPEQ=0.0
   POPNE=0.0
!
   DO 50 LVL=1,MAXSTATE*Max_ISO
      VQEQ(LVL)=NDEG (LVL)*EXP(-EH(LVL)/XKT)
      VQNEST(LVL)=VQNE(LVL)
      POPEQ=POPEQ + VQEQ(LVL)
      POPNE=POPNE + VQNE (LVL)
50 END DO
!
!    NORMALIZE POPULATIONS AND CALCULATE RATIOS
!
   DO 100 LVL=1,MAXSTATE*Max_ISO
      I=LVL
      VQEQ(LVL)=VQEQ(LVL)/POPEQ
      VQNE (LVL)=VQNE(LVL)/POPNE
      RAT(I)=VQNE(I)/VQEQ(I)
      IF(LVL.EQ.1) THEN
         VST1=VQNE(1)
         TNE(I)=TAVE
      ELSE
         DEN=(NDEG(1)*VQNE(LVL)/(NDEG(LVL)*VST1))
         TNE(I)=-RADCN2*EH(LVL)/ LOG(DEN)
      END IF
100 END DO
!
   WRITE(IPR,906) TITMOL
   WRITE(IPR,935)
   DO J=1,NUMIN
      I=ISORDER(J)
      ISOTOPE=(I-1)/MAXSTATE + 1
      K=I-(ISOTOPE-1)*MAXSTATE
      WRITE(IPR,920) ISOTOPE,IDX(K),EH(I),VQEQ(I),VQNE(I),RAT(I),    &
      &        TNE(I),VQNEST(I)
   END DO
!
!     READ TO THE END OF THE VIBRATIONAL DATA
!
   CALL RDSKIP(NLTEFLAG,TEXTLINE)
   RETURN
!
902 FORMAT(4x,A6)
904 FORMAT (F7.0,1P,7E11.4,     /(18X, 6E11.4))
906 FORMAT(//,A10,'  ENERGY LEVELS',10(/,20X,1PE11.4))
920 FORMAT(I7,2X,A10,4G15.5,F10.2,G15.5)
935 FORMAT ('ISOTOPE',2X,'VIB',10X,'E(CM-1)',11X,'POP LTE',7X,        &
   & 'POP NLTE',6X,'NLTE/LTE',7X,'NLTE TMP',7X,'NLTE POP ORIG')
!
end subroutine VIBPOP
!
! ----------------------------------------------------------------
!
SUBROUTINE VIBTMP(XKT,HT,NLTEFLAG,NUM,IDX,NDEG,EH,RAT,            &
&  TITMOL,TEXTLINE,ISORDER)
!
!
!     SUBROUTINE VIBTMP USES THE NON-LTE TEMPERATURE DATA FROM
!     TAPE4 TO CALCULATE THE VIBRATIONAL POPULATION ENHANCEMENT
!     RATIOS FOR SELECTED VIBRATIONAL STATES OF H2O, CO2, NO AND
!     O3.  THE NLTE VIBRATIONAL TEMPERATURES ARE WITH RESPECT TO
!     THE GROUND VIBRATIONAL STATE.
!
   USE phys_consts, ONLY: radcn2
   USE lblparams
   IMPLICIT REAL*8           (V)

   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4

   CHARACTER*5 IDX
   CHARACTER*80 TEXTLINE
   CHARACTER*8      XID,       HMOLID,      YID
   Real*8               SECANT,       XALTZ
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID

   DIMENSION IDX(MAXSTATE),VQNE(MAXSTATE*Max_ISO),                   &
   &     VQEQ(MAXSTATE*Max_ISO),RAT(MAXSTATE*Max_ISO)
   DIMENSION TNE(MAXSTATE*Max_ISO),TNESAV(MAXSTATE*Max_ISO),         &
   &     TEM1(MAXSTATE*Max_ISO),TEM2(MAXSTATE*Max_ISO)
   DIMENSION NDEG(MAXSTATE*Max_ISO),EH(MAXSTATE*Max_ISO),            &
   &     ISORDER(MAXSTATE*Max_ISO),TNEIN(MAXSTATE*Max_ISO)
   CHARACTER*6 TITMOL,HMNLTE,BLANKS
   DATA BLANKS/' '/
!
!   READ NLTE VIB TEMPERATURE
!       If this routine is called repeatedly for different altitudes,
!       this code is very inefficient, reading the input files every time
!       to search for the proper altitude
!
   XKT= TAVE/RADCN2
   DO I=1,MAXSTATE*Max_ISO
      IF(ISORDER(I).GT.0) NUMIN=I
   end do
   write(ipr,*) 'called vibtmp with txtmol=',titmol,'  alt=',ht
   write(ipr,*) 'num=',num,'  numin=',numin
   write(ipr,*) 'isorder=',(isorder(i),i=1,numin)

   I1=0
   I2=0
   DO I=1,80
      IF(TEXTLINE(I:I).NE.'-' .AND. TEXTLINE(I:I).NE.' ') THEN
         IF(I1.EQ.0) I1=I
      END IF
      IF(I1.GT.0 .AND. TEXTLINE(I:I).EQ.' ') THEN
         IF(I2.EQ.0) I2=I-1
      END IF
   END DO
   HMNLTE=TEXTLINE(I1:I2)//BLANKS
   IF(HMNLTE.NE.TITMOL) THEN
      READ (NLTEFLAG,902) HMNLTE
      i1=0
      DO I=1,6
         IF(HMNLTE(I:I).EQ.' ') I1=I
      END DO
      HMNLTE=HMNLTE(I1+1:6)//BLANKS
   END IF
   IF(HMNLTE.NE.TITMOL) THEN
      WRITE(IPR,*) 'READING VIBTMP DATA FROM TAPE4'
      WRITE(IPR,*) 'EXPECTED PROFILE DATA FOR ',TITMOL
      WRITE(IPR,*) 'BUT READ SPECIE ',HMNLTE
      STOP 'ERROR READING NLTE DATA FROM TAPE4'
   END IF
   READ (NLTEFLAG,*)  ALT1,(TEM1(I),I=1,NUMIN)
10 READ (NLTEFLAG,*)  ALT2,(TEM2(I),I=1,NUMIN)
   IF(TEM1(1).LE.1.E-12 .OR. TEM2(1).LE.1.E-12) THEN
      WRITE(IPR,*) 'VIBRATIONAL TEMPRATURE BASELINE FOR ',TITMOL,    &
         ' MISSSING'
      STOP 'KINETIC BASELINE TEMPERATURE MISSING'
   END IF
!*****WOG, 11/06/2000: ALT1 -> AL2:
!     IF( ALT1.LT.HT) THEN
   IF( ALT2.LE.HT) THEN
      ALT1 = ALT2
      DO 20 I=1,NUMIN
         TEM1(I) = TEM2(I)
20    CONTINUE
      GO TO 10
   ENDIF
!
!     SET ZERO TEMP TO AMBIENT TEMP (NOW DONE AFTER INTERPOLATION)
!
!      write(ipr,941) 'tem1',alt1,(tem1(i),i=1,numin)
!      write(ipr,941) 'tem2',alt2,(tem2(i),i=1,numin)
   CALL LININT(HT,ALT1,ALT2,NUMIN,TEM1,TEM2,TNEIN)
!      write(ipr,941) 'TNEIN',ht,(tnein(i),i=1,numin)
   DO I=1,MAXSTATE
      ! this rescales to ambient
      TNESAV(I)=TNEIN(1)
   END DO
   DO I=1,NUMIN
      INDEX=ISORDER(I)
      ! FOR ISOTOPE #1
      IF(INDEX.LE.MAXSTATE) THEN
         IF(TNEIN(I).GT.1.E-25) TNESAV(INDEX)=TNEIN(I)
      END IF
   END DO
!       FOR ISOTOPES>1, DEFAULT TO ISOTOPE #1 TEMPERATURE DATA
   DO ISOTOPE=2,Max_ISO
      DO K=1,MAXSTATE
         I= K + (ISOTOPE-1)*MAXSTATE
         TNESAV(I)=TNESAV(K)
      END DO
   END DO
   NUMISO=0
   DO I=1,NUMIN
      INDEX=ISORDER(I)
      ISOTOPE=(INDEX-1)/MAXSTATE +1
      IF(ISOTOPE.GT.NUMISO) NUMISO=ISOTOPE
      IF(ISOTOPE.GT.1) THEN
         IF(TNEIN(I).GT.1.E-25) THEN
            TNESAV(INDEX)=TNEIN(I)
         END IF
      END IF
   END DO
   WRITE(IPR,*) 'NUMBER OF ISOTOPES=',NUMISO
!      write(ipr,941) 'TNESAV',ht,(tnesav(i),i=1,MAXSTATE*NUMISO)
941 format(a6,2x,27f8.3/26F8.3/26F8.3)
!
!C     CORRECT TEMP TO ATMOSPHERIC
!
   RATTV = TAVE /TNESAV(1)
! loop 40 corrects the input temperatures when the input
!  level temperature does not match the computed layer
!  temperature
!      write(*,*) ' temperature not corrected'

   DO 40 I = 1,MAXSTATE*Max_ISO
      TNE(I) = TNESAV(I) * RATTV
40 END DO
!      write(ipr,941) 'TNE',ht,(tne(i),i=1,MAXSTATE*NUMISO)
!
   SUMQ=0
   SUMNQ=0
   DO 50 I=1,MAXSTATE*Max_ISO
      VQNE(I)=1.
      IF (TNE(I).GT.0.0) VQNE(I)=NDEG(I)*EXP(-RADCN2*EH(I)/TNE(I))
      VQEQ(I)=NDEG(I)*EXP(-EH(I)/XKT)
      SUMQ=SUMQ+VQEQ(I)
      SUMNQ=SUMNQ+VQNE(I)
50 END DO
   DO 100 I=1,MAXSTATE*Max_ISO
      VQNE(I)=VQNE(I)/SUMNQ
      VQEQ(I)=VQEQ(I)/SUMQ
      IF(VQNE(I).GT.0.) THEN
         RAT(I)=VQNE(I)/VQEQ(I)
      ELSE
         RAT(I)=1.0
      END IF
100 END DO

   WRITE(IPR,906) TITMOL
   WRITE(IPR,935)
   DO J=1,NUMIN
      I=ISORDER(J)
      ISOTOPE=(I-1)/MAXSTATE + 1
      K=I-(ISOTOPE-1)*MAXSTATE
      WRITE(IPR,920)ISOTOPE,IDX(K),EH(I),VQEQ(I),VQNE(I),RAT(I),     &
      &        TNESAV(I), TNE(I)
   END DO
!
!     READ TO THE END OF THE VIBRATIONAL DATA
!
   CALL RDSKIP(NLTEFLAG,TEXTLINE)
   RETURN
!
902 FORMAT(4x,A6)
904 FORMAT(F7.0,7F11.3/(18X,6F11.3))
906 FORMAT(//,5X,A10,'  ENERGY LEVELS',10(/,20X,1PE11.4))
920 FORMAT(I7,2X,A10,4G12.5,2F9.3)
935 FORMAT ('ISOTOPE',9X,'VIB E(CM-1)',9X,'POP LTE    POP NLTE ',     &
   &     'NLTE/LTE    NLTE TMP 2-STATE NLTE TMP')
!
end subroutine VIBTMP

! ----------------------------------------------------------------

SUBROUTINE RDSKIP(NTAPE,TEXTLINE)
   CHARACTER *1 HRD,HMINUS
   CHARACTER*80 TEXTLINE
   DATA HMINUS /'-'/
10 READ(NTAPE,900,END=50) TEXTLINE
900 FORMAT(A80)
   IF(TEXTLINE(1:2).EQ.'--') RETURN
   GO TO 10
50 RETURN
end subroutine RDSKIP

! ----------------------------------------------------------------

SUBROUTINE LININT(HT,ALT1,ALT2,NUM,T1,T2,TNE)
   USE lblparams
!      include 'lblparams.inc'
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   DIMENSION T1(*),T2(*),TNE(*)
!*****WOG 11/03/2000
!*****Correct for divide by zero if two altitudes are the same
!*****0.001 = 1 meter, small enough to use average, large enough to
!     prevent numerical error.
   IF (ABS(ALT2-ALT1) .LE. 0.001) THEN
      AM = 0.5
   ELSE
      AM = (HT-ALT1)/(ALT2-ALT1)
   ENDIF

   DO 10 I=1,NUM
      IF(ABS(T1(I)).LE.1.E-25 .AND. ABS(T2(I)).GT.1.E-25) GOTO 50
      IF(ABS(T1(I)).GT.1.E-25 .AND. ABS(T2(I)).LE.1.E-25) GOTO 50
      TNE(I)=T1(I)+AM*(T2(I)-T1(I))
10 END DO

!     DO 10 I=1,NUM
!         AM = (T2(I)-T1(I))/(ALT2-ALT1)
!         C=T1(I)-AM*ALT1
!         TNE(I)=AM*HT+C
!  10  CONTINUE
   RETURN
50 WRITE(IPR,*) 'ERROR IN LININT FOR Tvib INPUT'
   WRITE(IPR,*) 'Tvib=0 AT SOME ALTITUDES AND NOT OTHERS'
   write(ipr,*) 'ht=',ht,'  alt1,alt2=',alt1,alt2
   write(ipr,*) 'T1=',(t1(i),i=1,num)
   write(ipr,*) 'T2=',(t2(i),i=1,num)
   STOP 'ERROR IN LININT FOR Tvib INPUT'
end subroutine LININT

! ----------------------------------------------------------------

SUBROUTINE STAMB(NUM,T1AMB,T1)
   DIMENSION T1(*)

!*****WOG, 11/3/2000
!*****Why start from 2 instead of 1???
!     DO 10 I=2,NUM
   DO 10 I=1,NUM
      IF(T1(I).LE.0.) T1(I)=T1AMB
10 END DO
   RETURN
end subroutine STAMB
!
! ----------------------------------------------------------------
!
SUBROUTINE HIRACQ (MPTS)
!
   USE phys_consts, ONLY: radcn2
   USE lblparams
   IMPLICIT REAL*8           (V)
!
!
!**********************************************************************
!*
!*
!*    CALCULATES MONOCHROMATIC ABSORPTION COEFFICIENT FOR SINGLE LAYER
!*
!*
!*            USES APPROXIMATE VOIGT ALGORITHM
!*
!*
!*              VAN VLECK WEISSKOPF LINE SHAPE
!*
!**********************************************************************
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
!                                     M.W.SHEPHARD & W.O.GALLERY
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     131 Hartwell Ave,  Lexington,  MA   02421
!
!----------------------------------------------------------------------
!
!               WORK SUPPORTED BY:    THE ARM PROGRAM
!                                     OFFICE OF ENERGY RESEARCH
!                                     DEPARTMENT OF ENERGY
!
!
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!      include 'lblparams.inc'
!
!     Common blocks from analytic derivatives
!     -------------------------
   COMMON /ADRPNM/ CDUM1,PTHODI,PTHODTU,PTHODTD
   COMMON /ADRPTH/ PTHDIR,PTHRDRU,PTHRDRD,AJID
!     -------------------------
   COMMON /RCNTRL/ ILNFLG
   COMMON VNU(250),SP(250),ALFA0(250),EPP(250),MOL(250),HWHMS(250),  &
   &       TMPALF(250),PSHIFT(250),IFLG(250),SPPSP(250),RECALF(250),  &
   &       ZETAI(250),IZETA(250)
!
!     DIMENSION RR1 =  NBOUND   + 1 + DIM(R1)
!     DIMENSION RR2 =  NBOUND/2 + 1 + DIM(R2)
!     DIMENSION RR3 =  NBOUND/4 + 1 + DIM(R3)
!
   COMMON RR1(-8704:11169),RR2(-2560:3177),RR3(-1024:1179)
   COMMON /XRNLTE/ RR1s(-8704:11169),RR2s(-2560:3177),RR3s(-1024:1179)
   COMMON /IOU/ IOUT(250)
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,n_zero,NP,H1F,H2F,       &
   &                ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH,       &
   &                XVBAR, HMINF,PHIF,IERRF,HSPACE
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   CHARACTER*8      XID,       HMOLID,      YID
   Real*8               SECANT,       XALTZ
!
   COMMON /CVNLTE/ HNMNLTE,HVNLTE
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /XSUB/ VBOT,VTOP,VFT,LIMIN,ILO,IHI,IEOF,IPANEL,ISTOP,IDATA
   COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)
   COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,          &
   &                HWF3,DXF3,NX3,N3MAX
   COMMON /SUB1/ MAX1,MAX2,MAX3,NLIM1,NLIM2,NLIM3,NLO,NHI,DVR2,DVR3, &
   &              N1R1,N2R1,N1R2,N2R2,N1R3,N2R3
   COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /VOICOM/ AVRAT(102),CGAUSS(102),CF1(102),CF2(102),         &
   &                CF3(102),CER(102)
!
!      PARAMETER (NUMZ = 101)
!
   COMMON /FNSHQ/ IFN,F1(NFMX, NUMZ),F2(NFMX, NUMZ),                 &
   &     F3(NFMX, NUMZ), FG(NFMX)

   COMMON /R4SUB/ VLOF4,VHIF4,ILOF4,IST,IHIF4,LIMIN4,LIMOUT,ILAST,   &
   &               DPTMN4,DPTFC4,ILIN4,ILIN4T
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
!
   COMMON /ISVECT/ ISO_MAX(MXMOL),SMASSI(mxmol,10)
   COMMON /LNC1/ RHOSLF(mxmol),ALFD1(42,10),SCOR(42,10),ALFMAX,      &
   &              BETACR,DELTMP,DPTFC,DPTMN,XKT,NMINUS,NPLUS,NLIN,    &
   &              LINCNT,NCHNG,SUMALF,SUMZET,TRATIO,RHORAT,PAVP0,     &
   &              PAVP2,RECTLC,TMPDIF,ILC
   COMMON /FLFORM/ CFORM
   COMMON /L4TIMG/ L4TIM,L4TMR,L4TMS,L4NLN,L4NLS,LOTHER
   COMMON /IODFLG/ DVOUT

!     Total timing array for layer line-by-line calculation
   common /timing_lay_nlte/ time_lay_lbl(20)
!
   REAL L4TIM,L4TMR,L4TMS,LOTHER
   CHARACTER*55 CDUM1,PTHODI,PTHODTU,PTHODTD
   CHARACTER*11 PTHRDRU,PTHRDRD
   CHARACTER*3  PTHDIR,AJID
   CHARACTER*10 HFMODL
   CHARACTER CFORM*11,KODLYR*57,PTHODE*55,PTHODD*55
   CHARACTER*18 HNMNLTE,HVNLTE
   LOGICAL OP
!
   DIMENSION MEFDP(64),FILHDR(2),IWD(2)
   DIMENSION R1(4050),R2(1050),R3(300)
   DIMENSION SRAD(2)
!
   EQUIVALENCE (IHIRAC,FSCDID(1)) , (ILBLF4,FSCDID(2)),              &
   &            (IXSCNT,FSCDID(3)) , (IAERSL,FSCDID(4)),              &
   &            (JRAD,FSCDID(9)) , (IMRG,FSCDID(11)),                 &
   &            (IATM,FSCDID(15)) , (YI1,IOD) , (XID(1),FILHDR(1)),   &
   &            (V1P,IWD(1)) , (NPNLXP,LSTWDX),                       &
   &            (EPP(1),SRAD(1))
   EQUIVALENCE (R1(1), RR1(1)),(R2(1),RR2(1)),(R3(1),RR3(1))
!
!
!     NOTE that DXFF1 = (HWFF1/(NFPTS-1))
!     and       DXFF2 = (HWFF2/(NFPTS-1))
!     and       DXFF3 = (HWFF3/(NFPTS-1))
!
   DATA HWFF1 /  4. /,DXFF1 / 0.002 /,NXF1 / NFPTS /,NF1MAX / NFMX /
   DATA HWFF2 / 16. /,DXFF2 / 0.008 /,NXF2 / NFPTS /,NF2MAX / NFMX /
   DATA HWFF3 / 64. /,DXFF3 / 0.032 /,NXF3 / NFPTS /,NF3MAX / NFMX /
!
   DATA MEFDP / 64*0 /
!
   DATA I_10/10/
!
   PTHODE = 'ODexact_'
   PTHODD = 'ODdeflt_'
   DATA KODLYR /                                                     &
   &     '                                                         '/
   DATA HFMODL /'         '/
!
   CALL CPUTIM (TIMEH0)
!
!     ASSIGN NAME and CVS VERSION NUMBER TO MODULE
!
   HNMNLTE = '         nonlte.f:'
   HVNLTE  = '$Revision$'
!
!     Initialize timing for the group "OTHER" in the TAPE6 output
!
   TLNCOR = 0.0
   TXINT = 0.0
   TSHAPE = 0.0
   TLOOPS = 0.0
   TODFIL = 0.0
   TMOLEC = 0.0
!
   LSTWDX = -654321
   NPNLXP = NWDL(IWD,LSTWDX)
   ICNTNM = MOD(IXSCNT,I_10)
   IXSECT = IXSCNT/10
!
!     SET INPUT FLAG FOR USE BY X-SECTIONS
!
   IFST = -99
   IR4 = 0
   IENTER = 0
!
!     SET COMMON BLOCK CMSHAP
!
   HWF1 = HWFF1
   DXF1 = DXFF1
   NX1 = NXF1
   N1MAX = NF1MAX
   HWF2 = HWFF2
   DXF2 = DXFF2
   NX2 = NXF2
   N2MAX = NF2MAX
   HWF3 = HWFF3
   DXF3 = DXFF3
   NX3 = NXF3
   N3MAX = NF3MAX
!
   DPTMN = DPTMIN
   IF (JRAD.NE.1) DPTMN = DPTMIN/RADFN(V2,TAVE/RADCN2)
   DPTFC = DPTFAC
   ILIN4 = 0
   ILIN4T = 0
   NPTS = MPTS
   LIMIN = 250
   NSHIFT = 32
!
!     SAMPLE IS AVERAGE ALPHA / DV
! "*0.04/ALFAL0" is a temporary solution to keep ALFMAX to be close to the
! original default value. It may enlarge the ALFMAX in regions dominated
! Doppler broadening where ALFAL0 is supposed to have no effect if user
! sets ALFAL0 to a value less than 0.04.
!    previous expression:  NBOUND = 4.*(2.*HWF3)*SAMPLE+0.01
!
   ALFMAX = 4*SAMPLE*DV * 0.04/ALFAL0
   nALFMAX = ALFMAX/DV
   NBOUND = (2.*HWF3)*nALFMAX+0.01
!
   NLIM1 = 2401
   NLIM2 = (NLIM1/4)+1
   NLIM3 = (NLIM2/4)+1
!
   IF (IFN.EQ.0) THEN
      CALL CPUTIM(TPAT0)
      call voigt_init(f1, f2, f3)
      IFN = IFN+1
      CALL CPUTIM(TPAT1)
      TSHAPE = TSHAPE+TPAT1-TPAT0
   ENDIF
!
   CALL CPUTIM(TPAT0)
   CALL MOLEC (1,SCOR,RHOSLF,ALFD1)
   CALL CPUTIM(TPAT1)
   TMOLEC = TMOLEC+TPAT1-TPAT0
   REWIND LINFIL
   TIMRDF = 0.
   TIMCNV = 0.
   TIMPNL = 0.
   TF4 = 0.
   TF4RDF = 0.
   TF4CNV = 0.
   TF4PNL = 0.
   TXS = 0.
   TXSRDF = 0.
   TXSCNV = 0.
   TXSPNL = 0.
   IEOF = 0
   ILO = 0
   IHI = -999
   NMINUS = 0
   NPLUS = 0
!
!     NOTE (DXF3/DXF1) IS 16 AND (DXF3/DXF2) IS 4
!
   DVP = DV
   DVR2 = (DXF2/DXF1)*DV
   DVR3 = (DXF3/DXF1)*DV
! previous code
!      MAX1 = NSHIFT+NLIM1+(NBOUND/2)
!      MAX2 = MAX1/4
!      MAX3 = MAX1/16
!      MAX1 = MAX1+NSHIFT+1+16
!      MAX2 = MAX2+NSHIFT+1+4
!      MAX3 = MAX3+NSHIFT+1+1
   MAX1 = NSHIFT+NLIM1+NSHIFT+NBOUND/2
   MAX2 = MAX1/4
   MAX3 = MAX1/16
   MAX1 = MAX1 +  4*nALFMAX
   MAX2 = MAX2 + 16*nALFMAX/4 + 1  !+1 to account for rounding error
   MAX3 = MAX3 + 64*nALFMAX/16 + 1 !+1 to account for rounding error
!
!     FOR CONSTANTS IN PROGRAM  MAX1=4018  MAX2=1029  MAX3=282
!
   CALL CPUTIM(TPAT0)
   BOUND =  REAL(NBOUND)*DV/2.
   BOUNF3 = BOUND/2.
!      ALFMAX = BOUND/HWF3
   NLO = NSHIFT+1
   NHI = NLIM1+NSHIFT-1
   DO 10 I = 1, MAX1
      R1(I) = 0.
      RR1s(I) = 0.
10 END DO
   DO 20 I = 1, MAX2
      R2(I) = 0.
      RR2s(I) = 0.
20 END DO
   DO 30 I = 1, MAX3
      R3(I) = 0.
      RR3s(I) = 0.
30 END DO
   IF (ILBLF4.EQ.0) THEN
      DO 40 I = 1, 2502
         R4(I) = 0.
         RR4(I) = 0.
40    CONTINUE
   ENDIF
!
   IF (IATM.GE.1.AND.IATM.LE.5) CALL YDIH1 (H1F,H2F,ANGLEF,YID)
   CALL CPUTIM(TPAT1)
   TLOOPS = TLOOPS + TPAT1-TPAT0
!
!     ---------------------------------------------------------------
!
!     - If IOD = 1 or 4 then calculate optical depths for each
!       layer with DV = DVOUT (using DVSET if IOD=4) and maintain
!       separately. Use PTHODI as the name of the optical depth files.
!       This requires the format HFMODL, which is produced by
!       calling the SUBROUTINE QNTIFY.
!
!     - If IOD = 2 and IMERGE = 1 then calculate optical depths
!       for each layer using the exact DV of each layer
!       Use PTHODE as the name of the optical depth files.
!       This requires the format HFMODL, which is produced by
!       calling the SUBROUTINE QNTIFY.
!
!     - If IOD=3 and IMRG=1 then calculate layer optical depths and
!       and interpolate all layers to the dv of the final layer
!       (used for analytic derivative calculation)
!
!     - If calculating optical depths using the default procedure,
!       sending output to a different file for each layer (IEMIT=0,
!       IOD=0, and IMRG=1), then use PTHODI for the optical depth
!       pathnames.
!
!     - Otherwise, use TAPE10.  For IOD=1, calculate optical depths
!       for each layer with DV = DVOUT (from DVSET in TAPE5, carried
!       in by COMMON BLOCK /IODFLG/ (interpolation in PNLINT).
!
   CALL CPUTIM(TPAT0)
   IF ((IOD.EQ.1).OR.(IOD.EQ.4)) THEN
      CALL QNTIFY(PTHODI,HFMODL)
      WRITE (KODLYR,HFMODL) PTHODI,LAYER
      INQUIRE (UNIT=KFILE,OPENED=OP)
      IF (OP) CLOSE (KFILE)
      OPEN (KFILE,FILE=KODLYR,FORM=CFORM,STATUS='UNKNOWN')
      REWIND KFILE
      DVSAV = DV
      IF (DVOUT.NE.0.) DV = DVOUT
      CALL BUFOUT (KFILE,FILHDR(1),NFHDRF)
      DV = DVSAV
      IF (NOPR.EQ.0) WRITE (IPR,900) KFILE,DV,BOUNF3
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
      IF (IMRG.EQ.1) THEN
         CALL QNTIFY(PTHODI,HFMODL)
         WRITE (KODLYR,HFMODL) PTHODI,LAYER
         INQUIRE (UNIT=KFILE,OPENED=OP)
         IF (OP) CLOSE (KFILE)
         OPEN (KFILE,FILE=KODLYR,FORM=CFORM,STATUS='UNKNOWN')
         REWIND KFILE
         dv_lbl = DV
         DV = DVOUT
         CALL BUFOUT (KFILE,FILHDR(1),NFHDRF)
         DV = dv_lbl
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
         DVOUT = 0.0
         CALL BUFOUT (KFILE,FILHDR(1),NFHDRF)
      ENDIF
      IF (NOPR.EQ.0) WRITE (IPR,900) KFILE,DV,BOUNF3
   ENDIF
   CALL CPUTIM(TPAT1)
   TODFIL = TODFIL + TPAT1-TPAT0
!
   IF (IHIRAC.EQ.9) THEN
      DO 50 M = 1, NMOL
         WK(M) = 0.
50    CONTINUE
   ENDIF
!
!     ---------------------------------------------------------------
!
   VFT = V1- REAL(NSHIFT)*DV
   VBOT = V1-BOUND
   VTOP = V2+BOUND
!
   LINCNT = 0
   NLIN = 0
   AVALF = 0.
   AVZETA = 0.
   SUMALF = 0.
   SUMZET = 0.
   NCHNG = 0
   NLNCR = 0
!
   V1R4ST = V1R4
   V2R4ST = V2R4
   IF (ILBLF4.GE.1) CALL LBLF4Q (JRAD,V1R4ST,V2R4ST)
!
   IFPAN = 1
!
60 CONTINUE
!
   CALL CPUTIM (TIME0)
   IF (IEOF.NE.0) GO TO 80
!
!     THERE ARE (LIMIN * 9) QUANTITIES READ IN:
!     VNU,SP,ALFA0,EPP,MOL,HWHMS,TMPALF,PSHIFT,IFLG
!
   CALL RDLIN
!
   CALL CPUTIM (TIME)
   TIMRDF = TIMRDF+TIME-TIME0
!
   IF (IEOF.NE.0) GO TO 80
!
!     MODIFY LINE DATA FOR TEMPERATURE, PRESSURE, AND COLUMN DENSITY
!
   CALL CPUTIM(TPAT0)
   CALL LNCORQ (NLNCR,IHI,ILO,MEFDP)
   CALL CPUTIM(TPAT1)
   TLNCOR = TLNCOR+TPAT1-TPAT0
!
70 CONTINUE
!
   CALL CNVFNQ (VNU,SP,SRAD,SPPSP,RECALF,R1,R2,R3,RR1s,RR2s,RR3s,    &
   &    ZETAI,IZETA)
!
   IF (IPANEL.EQ.0) GO TO 60
!
80 CONTINUE
!
!        FOR FIRST PANEL     N1R1=   1    N1R2=  1    N1R3=  1
!     FOR SUBSEQUENT PANELS  N1R1=  33   *N1R2= 13   *N1R3=  6
!         FOR ALL PANELS     N2R1=2432   *N2R2=612   *N2R3=155
!
!            NOTE: THE VALUES FOR N1R2, N1R3, N2R2 AND N2R3 WHICH
!                  ARE MARKED WITH AN ASTERISK, CONTAIN A 4 POINT
!                  OFFSET WHICH PROVIDES THE NECESSARY OVERLAP FOR
!                  THE INTERPOLATION OF R3 INTO R2, AND R2 INTO R1.
!
   IF (IFPAN.EQ.1) THEN
      IFPAN = 0
      N1R1 = 1
      N1R2 = 1
      N1R3 = 1
   ELSE
      N1R1 = NSHIFT+1
      N1R2 = (NSHIFT/4)+1+4
      N1R3 = (NSHIFT/16)+1+3
   ENDIF
   N2R1 = NLIM1+NSHIFT-1
   N2R2 = NLIM2+(NSHIFT/4)-1+4
   N2R3 = NLIM3+(NSHIFT/16)-1+3
!
   IF (VFT.LE.0.) THEN
      CALL RSYM (R1,DV,VFT)
      CALL RSYM (R2,DVR2,VFT)
      CALL RSYM (R3,DVR3,VFT)
      CALL RSYM (RR1s,DV,VFT)
      CALL RSYM (RR2s,DVR2,VFT)
      CALL RSYM (RR3s,DVR3,VFT)
   ENDIF
!
   IF (IXSECT.GE.1.AND.IR4.EQ.0) THEN
      CALL CPUTIM (TIME0)
      CALL XSECTM (IFST,IR4)
      CALL CPUTIM (TIME)
      TXS = TXS+TIME-TIME0
   ENDIF
!
   CALL CPUTIM(TPAT0)
   IF (ILBLF4.GE.1) THEN
      CALL XINT (V1R4,V2R4,DVR4,R4,1.0,VFT,DVR3,R3,N1R3,N2R3)
      CALL XINT (V1R4,V2R4,DVR4,RR4,1.0,VFT,DVR3,RR3s,N1R3,N2R3)
   ENDIF
   IF (ICNTNM.GE.1)                                                  &
   &    CALL XINT (V1ABS,V2ABS,DVABS,ABSRB,1.,VFT,DVR3,R3,N1R3,N2R3)
   CALL CPUTIM(TPAT1)
   TXINT = TXINT + TPAT1-TPAT0
!
   CALL PANELQ (R1,R2,R3,RR1s,RR2s,RR3s,KFILE,JRAD,IENTER)
!
   IF (ISTOP.NE.1) THEN
      IF (ILBLF4.GE.1) THEN
         VF1 = VFT-2.*DVR4
! Matt Alvarado 20150819 Put in Yingtao fix for VF1
         VF1 = V1+floor((VF1-V1)/DVR4)*DVR4
         VF2 = VFT+2.*DVR4+ REAL(N2R3+4)*DVR3
         IF (VF2.GT.V2R4.AND.V2R4.NE.V2R4ST) THEN
            CALL LBLF4Q (JRAD,VF1,V2R4ST)
            IF (IXSECT.GE.1.AND.IR4.EQ.1) THEN
               CALL CPUTIM (TIME0)
               CALL XSECTM (IFST,IR4)
               CALL CPUTIM (TIME)
               TXS = TXS+TIME-TIME0
            ENDIF
         ENDIF
      ENDIF
      GO TO 70
   ENDIF
!
   CALL CPUTIM (TIMEH1)
   TIME = TIMEH1-TIMEH0-TF4-TXS
!
   IF (NOPR.NE.1) THEN
      IF (ILBLF4.GE.1) WRITE (IPR,905) DVR4,BOUND4
      IF (NMINUS.GT.0) WRITE (IPR,910) NMINUS
      IF (NPLUS.GT.0) WRITE (IPR,915) NPLUS
      TOTHHI = TLNCOR+TXINT+TSHAPE+TLOOPS+TODFIL+TMOLEC
      WRITE (IPR,920) L4TIM,L4TMR,L4TMS,LOTHER,L4NLN,L4NLS, TXS,     &
         TXSRDF,TXSCNV,TXSPNL, TF4,TF4RDF,TF4CNV,TF4PNL,ILIN4T,ILIN4,   &
         TIME,TIMRDF,TIMCNV,TIMPNL,TOTHHI, NLIN,LINCNT,NCHNG

!        Fill timing array
!         time_lay_lbl(1) = l4tim
!         time_lay_lbl(2) = l4tmr
!         time_lay_lbl(3) = 0.0
!         time_lay_lbl(4) = l4tms
!         time_lay_lbl(5) = lother
!         time_lay_lbl(6) = txs
!         time_lay_lbl(7) = txsrdf
!         time_lay_lbl(8) = txscnv
!         time_lay_lbl(9) = txspnl
!         time_lay_lbl(10) = 0.0
!         time_lay_lbl(11) = tf4
!         time_lay_lbl(12) = tf4rdf
!         time_lay_lbl(13) = tf4cnv
!         time_lay_lbl(14) = tf4pnl
!         time_lay_lbl(15) = 0.0
!         time_lay_lbl(16) = time
!         time_lay_lbl(17) = timrdf
!         time_lay_lbl(18) = timcnv
!         time_lay_lbl(19) = timpnl
!         time_lay_lbl(20) = tothhi

!        Accumulate timing array

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
      time_lay_lbl(9) = time_lay_lbl(9) + txspnl
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
         WRITE (IPR,922) time_lay_lbl(1),time_lay_lbl(2), time_lay_lbl(4), &
            time_lay_lbl(5), time_lay_lbl(6), time_lay_lbl(7),          &
            time_lay_lbl(8),time_lay_lbl(9), time_lay_lbl(11),          &
            time_lay_lbl(12),time_lay_lbl(13), time_lay_lbl(14),        &
            time_lay_lbl(16), time_lay_lbl(17), time_lay_lbl(18),       &
            time_lay_lbl(19),time_lay_lbl(20)

      ENDIF

      WRITE(IPR,935)
      IF (LINCNT.GE.1) THEN
         AVALF = SUMALF/ REAL(LINCNT)
         AVZETA = SUMZET/ REAL(LINCNT)
      ENDIF
      WRITE (IPR,925) AVALF,AVZETA
!
      DO 90 M = 1, NMOL
         IF (MEFDP(M).GT.0) WRITE (IPR,930) MEFDP(M),M
90    CONTINUE
   ENDIF
!
   RETURN
!
900 FORMAT ('0  * HIRAC1 *  OUTPUT ON FILE ',I5,10X,' DV = ',F12.8,   &
   &        10X,' BOUNDF3(CM-1) = ',F8.4)
905 FORMAT ('0 DV FOR LBLF4 = ',F10.5,5X,'BOUND FOR LBLF4 =',F10.4)
910 FORMAT ('0 -------------------------',I5,' HALF WIDTH CHANGES')
915 FORMAT ('0 +++++++++++++++++++++++++',I5,' HALF WIDTH CHANGES')
920 FORMAT ('0',20X,'TIME',11X,'READ',4X,'CONVOLUTION',10X,'PANEL',   &
   &        9X,'OTHER+',                                              &
   &        6X,'NO. LINES',3X,'AFTER REJECT',5X,'HW CHANGES',/,       &
   &        2x,'LINF4',3X,2F15.3,15X,2F15.3,2I15,/,                   &
   &        2X,'XSECT ',2X,4F15.3,/,2X,'LBLF4 ',2X,4F15.3,15X,2I15,/, &
   &        2X,'HIRAC1',2X,5F15.3,3I15)
921 FORMAT (2x,'LINF4',3X,2F15.3,15X,2F15.3,                          &
   &        2X,'XSECT ',2X,4F15.3,                                    &
   &        2X,'LBLF4 ',2X,4F15.3,                                    &
   &        2X,'HIRAC1',2X,5F15.3)
922 FORMAT ('0',20X,'TIME',11X,'READ',4X,'CONVOLUTION',10X,'PANEL',   &
   &        9X,'OTHER+',/,                                            &
   &        2x,'LINF4',3X,2F15.3,15X,2F15.3,/,                        &
   &        2X,'XSECT ',2X,4F15.3,/,2X,'LBLF4 ',2X,4F15.3,15X,/,      &
   &        2X,'HIRAC1',2X,5F15.3)
925 FORMAT ('0  * HIRAC1 *  AVERAGE WIDTH = ',F8.6,                   &
   &        ',  AVERAGE ZETA = ',F8.6)
930 FORMAT ('0 ********  HIRAC1  ********',I5,' STRENGTHS FOR',       &
   &        '  TRANSITIONS WITH UNKNOWN EPP FOR MOL =',I5,            &
   &        ' SET TO ZERO')
935 FORMAT (/,'0     + OTHER timing includes:',/,                     &
   &          '0             In LINF4:  MOLEC, BUFIN, BUFOUT, ',      &
   &          'NWDL, ENDFIL, and SHRINQ',/,                           &
   &          '0             In HIRAC:  LNCOR, XINT, SHAPEL, ',       &
   &          'SHAPEG, VERFN, MOLEC, and other loops and ',           &
   &          'file maintenance within HIRAC',/)
!
end subroutine HIRACQ
SUBROUTINE LNCORQ (NLNCR,IHI,ILO,MEFDP)
!
   USE phys_consts, ONLY: radcn2
   USE lblparams
!      include 'lblparams.inc'
   USE struct_types, ONLY: mxbrdmol, nlinerec
   IMPLICIT REAL*8           (V)
!
   CHARACTER*1 FREJ(nlinerec),HREJ,HNOREJ
   COMMON /RCNTRL/ ILNFLG
   COMMON VNU(nlinerec),S(nlinerec),ALFA0(nlinerec),EPP(nlinerec),MOL(nlinerec),HWHMS(nlinerec),   &
   &       TMPALF(nlinerec),PSHIFT(nlinerec),IFLG(nlinerec),SPPSP(nlinerec),RECALF(nlinerec),  &
   &       ZETAI(nlinerec),IZETA(nlinerec)
   common /brdmoldat/ brd_mol_flg(mxbrdmol,nlinerec),       &
   &     brd_mol_hw(mxbrdmol,nlinerec),brd_mol_tmp(MXBRDMOL,nlinerec),    &
   &     brd_mol_shft(mxbrdmol,nlinerec),sdep(nlinerec)
   integer*4 brd_mol_flg

   DIMENSION TMPCOR_ARR(MXBRDMOL),ALFA_TMP(MXBRDMOL)
   Real*8 ALFSUM

   COMMON /IOU/ IOUT(nlinerec)
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   CHARACTER*8      XID,       HMOLID,      YID
   Real*8               SECANT,       XALTZ
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /IFIL/   IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,     &
   &                NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,    &
   &                NLTEFL,LNFIL4,LNGTH4,IBRD
   COMMON /XSUB/   VBOT,VTOP,VFT,DUM(7)
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)
   COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,          &
   &                HWF3,DXF3,NX3,N3MAX
   COMMON /VOICOM/ AVRAT(102),CGAUSS(102),CF1(102),CF2(102),         &
   &                CF3(102),CER(102)

   COMMON /VBNLTE/ RATSTATE(MAXSTATE*Max_ISO,MXMOL),NUMSTATE(MXMOL)
!
   COMMON /ISVECT/ ISO_MAX(MXMOL),SMASSI(mxmol,10)
!
   COMMON /LNC1/ RHOSLF(mxmol),ALFD1(42,10),SCOR(42,10),ALFMAX,      &
   &              BETACR,DELTMP,DPTFC,DPTMN,XKT,NMINUS,NPLUS,NLIN,    &
   &              LINCNT,NCHNG,SUMALF,SUMZET,TRATIO,RHORAT,PAVP0,     &
   &              PAVP2,RECTLC,TMPDIF,ILC
   DIMENSION MEFDP(64),FILHDR(2),AMOL(250),SP(250)
   DIMENSION A(4),B(4),TEMPLC(4)
   DIMENSION SRAD(250)
!
   COMMON /PATH_ISOTPL/ ISOTPL,NISOTPL,                              &
   &                     ISOTPL_FLAG(MXMOL,MXISOTPL),                 &
   &                     ISOTPL_MAIN_FLAG(MXMOL),                     &
   &                     MOLNUM(MXMOL*MXISOTPL),                      &
   &                     ISOTPLNUM(MXMOL*MXISOTPL),                   &
   &                     WKI(MXMOL,MXISOTPL)
!
   EQUIVALENCE (MOL(1),AMOL(1)) , (S(1),SP(1)) , (EPP(1),SRAD(1))
   EQUIVALENCE (IHIRAC,FSCDID(1)) , (ILBLF4,FSCDID(2)),              &
   &            (IXSCNT,FSCDID(3)) , (IAERSL,FSCDID(4)),              &
   &            (JRAD,FSCDID(9)) , (XID(1),FILHDR(1))
!
   character*8 h_lncor1
!
   data h_lncor1/' lncor1 '/
!
!     TEMPERATURES FOR LINE COUPLING COEFFICIENTS
!
   DATA TEMPLC / 200.0,250.0,296.0,340.0 /
   DATA HREJ /'0'/,HNOREJ /'1'/
   DATA NWDTH /0/
!
   DATA I_1/1/, I_100/100/, I_1000/1000/
!
   NLNCR = NLNCR+1
   IF (NLNCR.EQ.1) THEN
!
      XKT0 = TEMP0/RADCN2
      XKT = TAVE/RADCN2
      DELTMP = ABS(TAVE-TEMP0)
      BETACR = (1./XKT)-(1./XKT0)
      CALL MOLEC (2,SCOR,RHOSLF,ALFD1)
!
      TRATIO = TAVE/TEMP0
      RHORAT = (PAVE/P0)*(TEMP0/TAVE)
!
      PAVP0 = PAVE/P0
      PAVP2 = PAVP0*PAVP0
!
!     FIND CORRECT TEMPERATURE AND INTERPOLATE FOR Y AND G
!
      DO 10 IL = 1, 3
         ILC = IL
         IF (TAVE.LT.TEMPLC(ILC+1)) GO TO 20
10    CONTINUE
20    IF (ILC.EQ.4) ILC = 3
!
      RECTLC = 1.0/(TEMPLC(ILC+1)-TEMPLC(ILC))
      TMPDIF = TAVE-TEMPLC(ILC)
!
   ENDIF
!
   IF (ILNFLG.EQ.2) READ(15)(FREJ(J),J=ILO,IHI)
!
   DO 30 J = ILO, IHI
      YI = 0.
      GI = 0.
      GAMMA1 = 0.
      GAMMA2 = 0.
      I = IOUT(J)
      IFLAG = IFLG(I)
      MFULL=MOL(I)
!           Molecule number for this line, last 2 digits of MFULL
      M = MOD(MOL(I),I_100)
!
!     ISO=(MOD(MOL(I),1000)-M)/100   IS PROGRAMMED AS:
!
!           Isotope number for this line, 3rd digit from right of MFULL
      ISO = MOD(MOL(I),I_1000)/100
!
!     check if lines are within allowed molecular and isotopic limits
!
      if (m.gt.mxmol .or. m.lt. 1) then
         call line_exception (1,ipr,h_lncor1,m,nmol,iso,iso_max)
         go to 25
      else if (iso .gt. iso_max(m)) then
         call line_exception (2,ipr,h_lncor1,m,nmol,iso,iso_max)
         go to 25
      endif
!
      MOL(I) = M
!
      IF (ISOTPL_FLAG(M,ISO).EQ.0) THEN
         SUI = S(I)*WK(M)
      ELSE
         SUI = S(I)*WKI(M,ISO)
      ENDIF
!
!MJA, 20150821 - using the VNU(I) approximation for the radiation term
!                can cause issues for wavenumbers around 1000 cm-1,
!                so use the full rad term instead
!         IF (JRAD.EQ.1) SUI = SUI*VNU(I)
      IF (JRAD.EQ.1) SUI = SUI*RADFN(VNU(I),XKT)
!
      IF (SUI.EQ.0.) GO TO 25
!
      NLIN = NLIN+1
!
!     Y'S AND G'S ARE STORED IN I+1 POSTION OF VNU,S,ALFA0,EPP...
!       A(1-4),  B(1-4) CORRESPOND TO TEMPERATURES TEMPLC(1-4) ABOVE
!
      IF (IFLAG.EQ.1.OR.IFLAG.EQ.3) THEN
         A(1) = VNU(I+1)
         B(1) = S(I+1)
         A(2) = ALFA0(I+1)
         B(2) = EPP(I+1)
         A(3) = AMOL(I+1)
         B(3) = HWHMS(I+1)
         A(4) = TMPALF(I+1)
         B(4) = PSHIFT(I+1)
!
!     CALCULATE SLOPE AND EVALUATE
!
         SLOPEA = (A(ILC+1)-A(ILC))*RECTLC
         SLOPEB = (B(ILC+1)-B(ILC))*RECTLC
!
         IF (IFLAG.EQ.1) THEN
            YI = A(ILC)+SLOPEA*TMPDIF
            GI = B(ILC)+SLOPEB*TMPDIF
         ELSE
            GAMMA1 = A(ILC)+SLOPEA*TMPDIF
            GAMMA2 = B(ILC)+SLOPEB*TMPDIF
         ENDIF
      ENDIF
!
!     IFLAG = 2 IS RESERVED FOR LINE COUPLING COEFFICIENTS ASSOCIATED
!               WITH AN EXACT TREATMENT (NUMERICAL DIAGONALIZATION)
!
!     IFLAG = 3 TREATS LINE COUPLING IN TERMS OF REDUCED WIDTHS
!
      VNU(I) = VNU(I)+RHORAT*PSHIFT(I)
      if(sum(brd_mol_flg(:,i)).gt.0.AND.ibrd.gt.0) then
         vnu(i) = vnu(i)+sum(rhoslf(1:mxbrdmol)*brd_mol_flg(:,i)* &
         &           (brd_mol_shft(:,i)-pshift(i)))
      endif
!
!     TEMPERATURE CORRECTION OF THE HALFWIDTH
!     SELF TEMP DEPENDENCE TAKEN THE SAME AS FOREIGN
!
      TMPCOR = TRATIO**TMPALF(I)
      ALFA0I = ALFA0(I)*TMPCOR
      HWHMSI = HWHMS(I)*TMPCOR
      ALFL = ALFA0I*(RHORAT-RHOSLF(m))+HWHMSI*RHOSLF(m)

      if(sum(brd_mol_flg(:,i)).gt.0.AND.ibrd.gt.0) then
         tmpcor_arr = tratio**brd_mol_tmp(:,i)
         alfa_tmp = brd_mol_hw(:,i)*tmpcor_arr
         alfsum = sum(rhoslf(1:mxbrdmol)*brd_mol_flg(:,i)*alfa_tmp)
         alfl = (rhorat-sum(rhoslf(1:mxbrdmol)*brd_mol_flg(:,i))) &
         &           *alfa0i + alfsum
         if(brd_mol_flg(m,i).eq.0)   &
         &           alfl = alfl + rhoslf(m)*(hwhmsi-alfa0i)
      end if

!
      IF (IFLAG.EQ.3) ALFL = ALFL*(1.0-GAMMA1*PAVP0-GAMMA2*PAVP2)
!
      ALFAD = VNU(I)*ALFD1(m,iso)
      ZETA = ALFL/(ALFL+ALFAD)
      ZETAI(I) = ZETA
      FZETA = 100.*ZETA
      IZ = FZETA + ONEPL
      IZETA(I) = IZ
      ZETDIF = FZETA - REAL(IZ-1)
      ALFV = (AVRAT(IZ)+ZETDIF*(AVRAT(IZ+1)-AVRAT(IZ)))*(ALFL+ALFAD)
      IF (ALFV.LT.DV) THEN
         ALFV = DV
         NMINAD = 1
      ELSE
         NMINAD = 0
      ENDIF
      IF (ALFV.GT.ALFMAX) THEN
         ALFV = ALFMAX
         NPLSAD = 1
      ELSE
         NPLSAD = 0
      ENDIF
!
      IF (HWF3*ALFV+VNU(I) .LT. VFT) GO TO 25
!
      RECALF(I) = 1./ALFV
!
!     TREAT TRANSITIONS WITH negative EPP AS SPECIAL CASE
!
!>>   an epp value between -0.9999 and 0.  cm-1 is taken as valid
!
!>>   an epp value of -1. is assumed set by hitran indicating an unknown
!     value: no temperature correction is performed
!
!>>   for an epp value of less than -1., it is assumed that value has
!     been provided as a reasonable value to be used for purposes of
!     temperature correction.  epp is set positive
!
      if (epp(i).le.-1.001) epp(i) = abs(epp(i))

      if (epp(i).le.-0.999) MEFDP(M) = MEFDP(M)+1

!     temperature correction:

      if (epp(i) .gt. -0.999) then
         SUI = SUI*SCOR(m,iso)* EXP(-EPP(I)*BETACR)*(1.+EXP(-VNU(I)/XKT)&
            )
      endif
!
      SP(I) = SUI*(1.+GI*PAVP2)
      SPPI = SUI*YI*PAVP0
      SPPSP(I) = SPPI/SP(I)
      SRAD(I)=0.0
!
! ---from nlte:
      IF (MFULL.GE.1000) THEN
         FREQ=VNU(I)
!              NLOW is 4th and 5th digit from right of MFULL
         NLOW=MOD(MFULL/1000,100)
!              NUPP is 6th and 7th digit from right of MFULL
         NUPP=MFULL/100000
         RLOW=1.0
         RUPP=1.0
         DELTA=EXP(-FREQ/XKT)
!
!     PICK OUT MOLECULAR TYPES WITH VIBRATIONAL STATES
!
         IF(NUMSTATE(M).GT.0) THEN
            IF (NLOW.GT.NUMSTATE(M)) STOP 'NLOW GT NUMSTATE IN LNCORQ'
            INDLOW=NLOW + (ISO-1)*MAXSTATE
            IF (NLOW.GT.0) RLOW=RATSTATE(INDLOW,M)
            IF (NUPP.GT.NUMSTATE(M)) STOP 'NUPP GT NUMSTATE IN LNCORQ'
            INDUPP=NUPP + (ISO-1)*MAXSTATE
            IF (NUPP.GT.0) RUPP=RATSTATE(INDUPP,M)
         ELSE
            PRINT 900,M
900         FORMAT('LNCORQ: MOL IN TROUBLE',I10)
            SP(I)=0.
            SRAD(I)=0.
            SPPSP(I)=0.
            GO TO 30
         END IF
!
!     RLOW AND RUPP NOW SET
!
         FNLTE=SP(I)/(1.0-DELTA)
         SP(I)=FNLTE*(RLOW-RUPP*DELTA)
         SRAD(I)=FNLTE*(RLOW-RUPP)
!
         IF (IFLAG .EQ. 0) THEN
            SPEAK = SP(I)*ABS(RECALF(I))
            SLFABS = SPEAK
            IF(SPEAK.GE.5.) THEN
               SLFABS = 1.
            ELSE
               IF(SPEAK.GT.0.01) SLFABS = 1.-EXP(-SPEAK)
            ENDIF
            TEST = SLFABS *(1.-SRAD(I)/SP(I))
            IF(TEST.LE.DPTMN) GOTO 25
         END IF

! --- above from nlte
      ELSE
!
         IF (IFLAG.EQ.0) THEN
            IF (ILNFLG.LE.1) THEN
               FREJ(J) = HNOREJ
               SPEAK = SUI*RECALF(I)
               IF (DVR4.LE.0.) THEN
                  IF (SPEAK.LE.DPTMN) THEN
                     FREJ(J) = HREJ
                     GO TO 25
                  ENDIF
               ELSE
                  JJ = (VNU(I)-V1R4)/DVR4+1.
                  JJ = MAX(JJ,1)
                  JJ = MIN(JJ,NPTR4)
!                         IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ))) THEN
                  IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ)) .and. sppsp(i)   &
                     .eq.0.) THEN
                     FREJ(J) = HREJ
                     GO TO 25
                  ENDIF
               ENDIF
            ELSE
!     "ELSE" IS TRUE WHEN "ILNFLG" EQUALS 2
!
               IF (FREJ(J).EQ.HREJ) GO TO 25
            ENDIF
         ENDIF
      ENDIF

      NMINUS = NMINUS+NMINAD
      NPLUS = NPLUS+NPLSAD
      SUMALF = SUMALF+ALFV
      SUMZET = SUMZET+ZETA
      LINCNT = LINCNT+1
!
      GO TO 30
!
25    SP(I)=0.0
      SRAD(I)=0.0
      SPPSP(I) = 0.0
!
30 END DO
!
   NCHNG = NMINUS+NPLUS
   IF (ILNFLG.EQ.1) WRITE(15)(FREJ(J),J=ILO,IHI)
!
   RETURN
!
end subroutine LNCORQ
!-----------------------------------------------------------------------
SUBROUTINE CNVFNQ (VNU,SP,SRAD,SPPSP,RECALF,R1,R2,R3,RR1,         &
&    RR2,RR3,ZETAI,IZETA)
!
   USE lblparams
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE CNVFNV PERFORMS THE CONVOLUTION OF THE LINE DATA WITH
!     THE VOIGT LINE SHAPE (APPROXIMATED)
!
!     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!
!     IMPLEMENTATION:    R.D. WORSHAM
!
!     ALGORITHM REVISIONS:    S.A. CLOUGH
!     R.D. WORSHAM
!     J.L. MONCET
!
!
!     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!     ------------------------------------------------------------------
!
!     WORK SUPPORTED BY:    THE ARM PROGRAM
!     OFFICE OF ENERGY RESEARCH
!     DEPARTMENT OF ENERGY
!
!
!     SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!     FASCOD3
!
!
!     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   CHARACTER*8      XID,       HMOLID,      YID
   Real*8               SECANT,       XALTZ
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /XSUB/ VBOT,VTOP,VFT,LIMIN,ILO,IHI,IEOF,IPANEL,ISTOP,IDATA
   COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTS
   COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,          &
   &                HWF3,DXF3,NX3,N3MAX
   COMMON /SUB1/ MAX1,MAX2,MAX3,NLIM1,NLIM2,NLIM3,NLO,NHI,DVR2,DVR3, &
   &              N1R1,N2R1,N1R2,N2R2,N1R3,N2R3
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /VOICOM/ AVRAT(102),CGAUSS(102),CF1(102),CF2(102),         &
   &                CF3(102),CER(102)
   COMMON /IOU/ IOUT(250)
!
!      PARAMETER (NUMZ = 101)
   COMMON /FNSHQ/ IFN,F1(NFMX, NUMZ),F2(NFMX, NUMZ),F3(NFMX, NUMZ),  &
   &     FG(NFMX)
!
   DIMENSION VNU(*),SP(*),SRAD(*),SPPSP(*),RECALF(*)
   DIMENSION R1(*),R2(*),R3(*)
   DIMENSION RR1(-8704:11169),RR2(-2560:3177),RR3(-1024:1179)
!      DIMENSION RR1(*),RR2(*),RR3(*)
   DIMENSION IZETA(*),ZETAI(*)
!
   equivalence (zetdif,A)
!
   CALL CPUTIM (TIME0)
!
   CLC1 = 4./( REAL(NX1-1))
   CLC2 = 16./( REAL(NX2-1))
   CLC3 = 64./( REAL(NX3-1))
   WAVDXF = DV/DXF1
   HWDXF = HWF1/DXF1
   CONF2 = DV/DVR2
   CONF3 = DV/DVR3
   ILAST = ILO-1
!
   IF (ILO.LE.IHI) THEN
      DO 30 J = ILO, IHI
         I = IOUT(J)
         IF (SP(I).NE.0.) THEN
            DEPTHA = SP(I)*RECALF(I)
            DEPTHR = SRAD(I)*RECALF(I)
            IZM = IZETA(I)
            ZETDIF = 100.*ZETAI(I)- REAL(IZM-1)
!
            ZSLOPE = RECALF(I)*WAVDXF
            ZINT = (VNU(I)-VFT)/DV
            BHWDXF = HWDXF/ZSLOPE
            JMAX1 = ZINT+BHWDXF+1.5
            IF (JMAX1.GT.MAX1) THEN
               ILAST = J-1
               IPANEL = 1
               GO TO 40
            ENDIF
            JMIN1 = ZINT-BHWDXF+1.5
            RSHFT = 0.5
            IF (ZINT.LT.0.0) RSHFT = -RSHFT
            J2SHFT = ZINT*(1.-CONF2)+RSHFT
            J3SHFT = ZINT*(1.-CONF3)+RSHFT
            JMIN2 = JMIN1-J2SHFT
            JMIN3 = JMIN1-J3SHFT
            ZF1 = ( REAL(JMIN1-2)-ZINT)*ZSLOPE
            ZF2 = ( REAL(JMIN2-2)-ZINT*CONF2)*ZSLOPE
            ZF3 = ( REAL(JMIN3-2)-ZINT*CONF3)*ZSLOPE
!
            IF (SPPSP(I).eq.0.) then
!
               DO 10 J1 = JMIN1, JMAX1
                  J2 = J1-J2SHFT
                  J3 = J1-J3SHFT
                  ZF3 = ZF3+ZSLOPE
                  ZF2 = ZF2+ZSLOPE
                  ZF1 = ZF1+ZSLOPE
                  IZ3 = ABS(ZF3)+1.5
                  IZ2 = ABS(ZF2)+1.5
                  IZ1 = ABS(ZF1)+1.5
!     ************Using new voigt scheme:
!     ************Interpolate voigt subfunctions to zeta

!          A  is equivalenced to ZETDIF

                  x3 = deptha*( (1.-A)*F3(IZ3,IZM)+A*F3(IZ3,IZM+1))
                  x2 = deptha*(( 1.-A)*F2(IZ2,IZM)+A*F2(IZ2,IZM+1))
                  x1 = deptha*( (1.-A)*F1(IZ1,IZM)+A*F1(IZ1,IZM+1))

                  R3(J3) = R3(J3)+ x3
                  R2(J2) = R2(J2)+ x2
                  R1(J1) = R1(J1)+ x1

                  IF (DEPTHR.NE.0) THEN

                     xx3 = depthr*( (1.-A)*F3(IZ3,IZM)+A*F3(IZ3,IZM+ &
                        1))
                     xx2 = depthr*(( 1.-A)*F2(IZ2,IZM)+A*F2(IZ2,IZM+ &
                        1))
                     xx1 = depthr*( (1.-A)*F1(IZ1,IZM)+A*F1(IZ1,IZM+ &
                        1))

                     RR3(J3) = RR3(J3)+ xx3
                     RR2(J2) = RR2(J2)+ xx2
                     RR1(J1) = RR1(J1)+ xx1
                  ENDIF
!
10             CONTINUE
!
            else
!
!                 THE FOLLOWING DOES LINE COUPLING
!
!                 SPPSP(I) = SPP(I)/SP(I)
!
               dptrat1 = SPPSP(I)*clc1
               dptrat2 = SPPSP(I)*clc2
               dptrat3 = SPPSP(I)*clc3
!
               DO 20 J1 = JMIN1, JMAX1
!
                  J2 = J1-J2SHFT
                  J3 = J1-J3SHFT
                  ZF3 = ZF3+ZSLOPE
                  ZF2 = ZF2+ZSLOPE
                  ZF1 = ZF1+ZSLOPE
                  IZ3 = ABS(ZF3)+1.5
                  IZ2 = ABS(ZF2)+1.5
                  IZ1 = ABS(ZF1)+1.5
!
!     ************Using new voigt scheme:
!     ************Interpolate voigt subfunctions to zeta

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
                     xx3 = DEPTHr*( (1.-A)*F3(IZ3,IZM)+A*F3(IZ3,IZM+ &
                        1))
                     xx2 = DEPTHr*( (1.-A)*F2(IZ2,IZM)+A*F2(IZ2,IZM+ &
                        1))
                     xx1 = DEPTHr*( (1.-A)*F1(IZ1,IZM)+A*F1(IZ1,IZM+ &
                        1))

                     yy3 = dptrat3*xx3*ZF3
                     yy2 = dptrat2*xx2*ZF2
                     yy1 = dptrat1*xx1*ZF1

                     RR3(J3) = RR3(J3) + xx3 + yy3
                     RR2(J2) = RR2(J2) + xx2 + yy2
                     RR1(J1) = RR1(J1) + xx1 + yy1

                  ENDIF

20             CONTINUE
!
            ENDIF
         ENDIF
!
30    CONTINUE
      ILAST = IHI
!
!        IDATA=0 FOR MORE DATA REQUIRED
!        IDATA=1 IF NO MORE DATA REQUIRED
!
      IPANEL = IDATA
   ELSE
      IPANEL = 1
   ENDIF
!
40 ILO = ILAST+1
   CALL CPUTIM (TIME)
   TIMCNV = TIMCNV+TIME-TIME0
   RETURN
!
end subroutine CNVFNQ
SUBROUTINE PANELQ (R1,R2,R3,RR1,RR2,RR3,KFILE,JRAD,IENTER)
!
   USE phys_consts, ONLY: radcn2
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE PANEL COMBINES RESULTS OF R3, R2, AND R1 INTO R1 ARRAY
!     AND OUTPUTS THE ARRAY R1
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!               LAST MODIFICATION:    28 AUGUST 1992
!
!                  IMPLEMENTATION:    R.D. WORSHAM
!
!             ALGORITHM REVISIONS:    S.A. CLOUGH
!                                     R.D. WORSHAM
!                                     J.L. MONCET
!
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!----------------------------------------------------------------------
!
!               WORK SUPPORTED BY:    THE ARM PROGRAM
!                                     OFFICE OF ENERGY RESEARCH
!                                     DEPARTMENT OF ENERGY
!
!
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
   CHARACTER*8      XID,       HMOLID,      YID
   Real*8               SECANT,       XALTZ
!
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /XSUB/ VBOT,VTOP,VFT,LIMIN,ILO,IHI,IEOF,IPANEL,ISTOP,IDATA
   COMMON /SUB1/ MAX1,MAX2,MAX3,NLIM1,NLIM2,NLIM3,NLO,NHI,DVR2,DVR3, &
   &              N1R1,N2R1,N1R2,N2R2,N1R3,N2R3
   COMMON /XPANEL/ V1P,V2P,DVP,NLIM,RMIN,RMAX,NPNLXP,NSHIFT,NPTS
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KDUMM,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   COMMON /IODFLG/ DVOUT
   DIMENSION R1(*),R2(*),R3(*)
   DIMENSION RR1(-8704:11169),RR2(-2560:3177),RR3(-1024:1179)
!     DIMENSION RR1(*),RR2(*),RR3(*)
   DIMENSION PNLHDR(2)
!
   EQUIVALENCE (V1P,PNLHDR(1))
!
   CALL CPUTIM (TIME0)
   X00 = -7./128.
   X01 = 105./128.
   X02 = 35./128.
   X03 = -5./128.
   X10 = -1./16.
   X11 = 9./16.
   ISTOP = 0
!
!     Test for last panel.  If last, set the last point to one point
!     greater than V1 specified on TAPE5 (to ensure last point for
!     every layer is the same)
!
   IF ((VFT+(NHI-1)*DVP).GT.V2) THEN
      NHI = (V2-VFT)/DVP + 1.
      V2P = VFT+ REAL(NHI-1)*DVP
      IF (V2P.LT.V2) THEN
         V2P = V2P+DVP
         NHI = NHI+1
      ENDIF
      ISTOP = 1
   ELSE
      V2P = VFT+ REAL(NHI-1)*DV
   ENDIF
   NLIM = NHI-NLO+1
   V1P = VFT+ REAL(NLO-1)*DV
!
   LIMLO = N1R2
   IF (N1R2.EQ.1) LIMLO = LIMLO+4
   LIMHI = (NHI/4)+1
!
   DO 10 J = LIMLO, LIMHI, 4
      J3 = (J-1)/4+1
      R2(J) = R2(J)+R3(J3)
      R2(J+1) = R2(J+1)+X00*R3(J3-1)+X01*R3(J3)+X02*R3(J3+1)+        &
         X03*R3(J3+2)
      R2(J+2) = R2(J+2)+X10*(R3(J3-1)+R3(J3+2))+ X11*(R3(J3)+R3(J3+1)&
         )
      R2(J+3) = R2(J+3)+X03*R3(J3-1)+X02*R3(J3)+X01*R3(J3+1)+        &
         X00*R3(J3+2)
      RR2(J) = RR2(J)+RR3(J3)
      RR2(J+1) = RR2(J+1)+X00*RR3(J3-1)+X01*RR3(J3)+X02*RR3(J3+1)+   &
         X03*RR3(J3+2)
      RR2(J+2) = RR2(J+2)+X10*(RR3(J3-1)+RR3(J3+2))+ X11*(RR3(J3)+   &
         RR3(J3+1))
      RR2(J+3) = RR2(J+3)+X03*RR3(J3-1)+X02*RR3(J3)+X01*RR3(J3+1)+   &
         X00*RR3(J3+2)
10 END DO
   !    !--- If the last panel, interpolate the first point of the next 4-DV2 segment
   !    ! The first point is exactly aligned, so interpolation is simply taking the 
   !    ! corresponding R3 value.
   if (LIMLO<=LIMHI .and. ISTOP==1) then
      J3 = J3 + 1
      R2(J) = R2(J)+R3(J3)
      R2(J+1) = R2(J+1)+X00*R3(J3-1)+X01*R3(J3)+X02*R3(J3+1)+X03*R3(J3+2)

      RR2(J) = RR2(J)+RR3(J3)
      RR2(J+1) = RR2(J+1)+X00*RR3(J3-1)+X01*RR3(J3)+X02*RR3(J3+1)+X03*RR3(J3+2)
   endif
   DO 20 J = NLO, NHI, 4
      J2 = (J-1)/4+1
      R1(J) = R1(J)+R2(J2)
      R1(J+1) = R1(J+1)+X00*R2(J2-1)+X01*R2(J2)+X02*R2(J2+1)+        &
         X03*R2(J2+2)
      R1(J+2) = R1(J+2)+X10*(R2(J2-1)+R2(J2+2))+ X11*(R2(J2)+R2(J2+1)&
         )
      R1(J+3) = R1(J+3)+X03*R2(J2-1)+X02*R2(J2)+X01*R2(J2+1)+        &
         X00*R2(J2+2)
      RR1(J) = RR1(J)+RR2(J2)
      RR1(J+1) = RR1(J+1)+X00*RR2(J2-1)+X01*RR2(J2)+X02*RR2(J2+1)+   &
         X03*RR2(J2+2)
      RR1(J+2) = RR1(J+2)+X10*(RR2(J2-1)+RR2(J2+2))+ X11*(RR2(J2)+   &
         RR2(J2+1))
      RR1(J+3) = RR1(J+3)+X03*RR2(J2-1)+X02*RR2(J2)+X01*RR2(J2+1)+   &
         X00*RR2(J2+2)
20 END DO
!
!     IN THE FOLLOWING SECTION THE ABSORPTION COEFFICIENT IS MULTIPIIED
!     BY THE RADIATION FIELD
!
   IF (JRAD.EQ.0) THEN
!
      XKT = TAVE/RADCN2
      VI = V1P-DV
      VITST = VI
      RDLAST = -1.
      NPTSI1 = NLO-1
      NPTSI2 = NLO-1
!
30    NPTSI1 = NPTSI2+1
!
      VI = VFT+ REAL(NPTSI1-1)*DV
      RADVI = RADFNI(VI,DV,XKT,VITST,RDEL,RDLAST)
!
!         NPTSI2 = (VITST-VFT)/DV+1.001
!MJA 20150819 Implementing Yingtao Fix
      NPTSI2 = (VITST-VFT)/DV+0.001
      NPTSI2 = MIN(NPTSI2,NHI)
!
      DO 40 I = NPTSI1, NPTSI2
!           VI = VI+DV
         R1(I) = R1(I)*RADVI
         rr1(i) = rr1(i)*radvi
         RADVI = RADVI+RDEL
40    CONTINUE
!
      IF (NPTSI2.LT.NHI) GO TO 30
!
   ENDIF
!
!     V1P IS FIRST FREQ OF PANEL
!     V2P IS LAST FREQ OF PANEL
!
   IF (DVOUT.EQ.0.) THEN
      CALL BUFOUT (KFILE,PNLHDR(1),NPHDRF)
      CALL BUFOUT (KFILE,R1(NLO),NLIM)
      CALL BUFOUT (KFILE,RR1(NLO),NLIM)
!
      IF (NPTS.GT.0) CALL R1PRNT (V1P,DVP,NLIM,R1,NLO,NPTS,KFILE,    &
         IENTER)
      IF (NPTS.GT.0) CALL R1PRNT (V1P,DVP,NLIM,RR1,NLO,NPTS,KFILE,   &
         IENTER)
   ELSE
      CALL PNLINT (R1(NLO),IENTER)
      CALL PNLINT (RR1(NLO),IENTER)
   ENDIF
!
   VFT = VFT+ REAL(NLIM1-1)*DV
   IF (ISTOP.NE.1) THEN
      DO 50 J = NLIM1, MAX1
         R1(J-NLIM1+1) = R1(J)
         RR1(J-NLIM1+1) = RR1(J)
50    CONTINUE
      DO 60 J = MAX1-NLIM1+2, MAX1
         R1(J) = 0.
         RR1(J) = 0.
60    CONTINUE
      DO 70 J = NLIM2, MAX2
         R2(J-NLIM2+1) = R2(J)
         RR2(J-NLIM2+1) = RR2(J)
70    CONTINUE
      DO 80 J = MAX2-NLIM2+2, MAX2
         R2(J) = 0.
         RR2(J) = 0.
80    CONTINUE
      DO 90 J = NLIM3, MAX3
         R3(J-NLIM3+1) = R3(J)
         RR3(J-NLIM3+1) = RR3(J)
90    CONTINUE
      DO 100 J = MAX3-NLIM3+2, MAX3
         R3(J) = 0.
         RR3(J) = 0.
100   CONTINUE
      NLO = NSHIFT+1
   ENDIF
   CALL CPUTIM (TIME)
   TIMPNL = TIMPNL+TIME-TIME0
!
   RETURN
!
end subroutine PANELQ

! ----------------------------------------------------------------
SUBROUTINE LINF4Q (V1L4,V2L4)
!
   USE phys_consts, ONLY: radcn2
   USE lblparams
   USE struct_types
!
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE LINF4 READS THE LINES AND SHRINKS THE LINES FOR LBLF4
!
   COMMON /ISVECT/ ISO_MAX(MXMOL),SMASSI(mxmol,10)
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
!
   REAL*8            HID,HMOLIL,HID1,HLINHD,VMNCPL,VMXCPL,EXTSPC
!
   COMMON /BUFID/ HID(10),HMOLIL(64),MOLCNT(64),MCNTLC(64),          &
   &               MCNTNL(64),SUMSTR(64),NMOI,FLINLO,FLINHI,          &
   &               ILIN,ILINLC,ILINNL,IREC,IRECTL,HID1(2),LSTWDL
!
   TYPE(LINE_SHRINK)  :: SHRUNK
!
!      COMMON VNU(1250),SP(1250),ALFA0(1250),EPP(1250),MOL(1250),        &
!     &       SPP(1250),SRAD(1250)
!
   COMMON /IOU/ IOUT(250)

   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   CHARACTER*8      XID,       HMOLID,      YID
   Real*8               SEC   ,       XALTZ
!
   COMMON /FILHDR/ XID(10),SEC   ,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                W(60),PZL,PZU,TZL,TZU,WBROAD,DVO,V1 ,V2 ,TBOUND,  &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /R4SUB/ VLO,VHI,ILO,IST,IHI,LIMIN,LIMOUT,ILAST,DPTMN,      &
   &               DPTFC,ILIN4,ILIN4T
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4,IBRD
   COMMON /TPANEL/ VNULO,VNUHI,JLIN,NLNGT4,lstdum
!
   TYPE(LINE_DATA)  :: BUFR
!      COMMON /BUFR/ VNUB(250),SB(250),ALB(250),EPPB(250),MOLB(250),     &
!     &              HWHMB(250),TMPALB(250),PSHIFB(250),IFLG(250)
   COMMON /NGT4/ VD,SD,AD,EPD,MOLD,SPPD,ILS2D
   COMMON /L4TIMG/ L4TIM,L4TMR,L4TMS,L4NLN,L4NLS,LOTHER
   COMMON /VBNLTE/ RATSTATE(MAXSTATE*Max_ISO,MXMOL),NUMSTATE(MXMOL)
!
   DIMENSION TMPCOR_ARR(MXBRDMOL),ALFA_TMP(MXBRDMOL)
   Real*8 ALFSUM
!
   REAL L4TIM,L4TMR,L4TMS,LOTHER
   DIMENSION MEFDP(64)
   DIMENSION SCOR(42,10),RHOSLF(mxmol),ALFD1(42,10)
   DIMENSION A(4),B(4),TEMPLC(4)
!      DIMENSION ALFAL(1250),ALFAD(1250),A(4),B(4),TEMPLC(4)
   DIMENSION RCDHDR(2),IWD(2),IWD3(2),HLINHD(2) !,AMOLB(250)
!
   COMMON /PATH_ISOTPL/ ISOTPL,NISOTPL,                              &
   &                     ISOTPL_FLAG(MXMOL,MXISOTPL),                 &
   &                     ISOTPL_MAIN_FLAG(MXMOL),                     &
   &                     MOLNUM(MXMOL*MXISOTPL),                      &
   &                     ISOTPLNUM(MXMOL*MXISOTPL),                   &
   &                     WKI(MXMOL,MXISOTPL)
!
   EQUIVALENCE (IHIRAC,FSCDID(1)) , (ILBLF4,FSCDID(2))
   EQUIVALENCE (VNULO,RCDHDR(1)) , (IWD3(1),VD),                     &
!    &            (HLINHD(1),HID(1),IWD(1)) , (MOLB(1),AMOLB(1)),       &
!    &    (SP(1),SABS(1))
   &            (HLINHD(1),HID(1),IWD(1))
!
   character*8 h_linf4
!
   data h_linf4/' linf4  '/
   DATA MEFDP / 64*0 /
!
   data jrad4 /0/
!     the fourth function is always computed without the radiation field
!
!     TEMPERATURES FOR LINE COUPLING COEFFICIENTS
!
   DATA TEMPLC / 200.0,250.0,296.0,340.0 /
!
   DATA I_100/100/, I_1000/1000/
!
!     Initialize timing for the group "OTHER" in the TAPE6 output
!
   LOTHER = 0.0
   TSHRNK = 0.0
   TBUFFR = 0.0
   TMOLN4 = 0.0
!
   CALL CPUTIM (TIMEL0)
!
   ILS2D = -654321
   NLNGT4 = NWDL(IWD3,ILS2D)*1250
   LNGTH4 = NLNGT4
   PAVP0 = PAVE/P0
   PAVP2 = PAVP0*PAVP0
   DPTMN = DPTMIN/RADFN(V2,TAVE/RADCN2)
   DPTFC = DPTFAC
   LIMIN = 1000
!
   CALL CPUTIM(TPAT0)
   CALL MOLEC (1,SCOR,RHOSLF,ALFD1)
   CALL CPUTIM(TPAT1)
   TMOLN4 = TMOLN4 + TPAT1-TPAT0
!
   TIMR = 0.
   TIMS = 0.
   SUMS = 0.
   ILAST = 0
   ILINLO = 0
   ILINHI = 0
   ILO = 1
   IST = 1
   NLINS = 0
   NLIN = 0
!
   VLO = V1L4
   VHI = V2L4
!
   CALL CPUTIM(TPAT0)
!
   call lnfilhd_4(linfil,lnfil4,v1,v2)
!
   CALL CPUTIM(TPAT1)
   TBUFFR = TBUFFR + TPAT1-TPAT0
!
!       TEMPERATURE CORRECTION TO INTENSITY
!       TEMPERATURE AND PRESSURE CORRECTION TO HALF-WIDTH
!
   TRATIO = TAVE/TEMP0
   RHORAT = (PAVE/P0)*(TEMP0/TAVE)
!
   BETA = RADCN2/TAVE
   BETA0 = RADCN2/TEMP0
   BETACR = BETA-BETA0
   DELTMP = ABS(TAVE-TEMP0)
   CALL CPUTIM(TPAT0)
   CALL MOLEC (2,SCOR,RHOSLF,ALFD1)
   CALL CPUTIM(TPAT1)
   TMOLN4 = TMOLN4 + TPAT1-TPAT0
!
!     FIND CORRECT TEMPERATURE AND INTERPOLATE FOR Y AND G
!
   DO 10 ILC = 1, 4
      IF (TAVE.LE.TEMPLC(ILC)) GO TO 20
10 END DO
20 IF (ILC.EQ.1) ILC = 2
   IF (ILC.EQ.5) ILC = 4
   RECTLC = 1.0/(TEMPLC(ILC)-TEMPLC(ILC-1))
   TMPDIF = TAVE-TEMPLC(ILC)
!
   IJ = 0
30 CALL CPUTIM (TIM0)

   CALL RDLNFL (IEOF,ILINLO,ILINHI, BUFR)
   CALL CPUTIM (TIM1)
   TIMR = TIMR+TIM1-TIM0
!
   IF (IEOF.GE.1) GO TO 60
!
   DO 50 J = ILINLO, ILINHI
      YI = 0.
      GI = 0.
      GAMMA1 = 0.
      GAMMA2 = 0.
      I = IOUT(J)
      IFLAG = BUFR%IFLG(I)
      IF (I.LE.0) GO TO 50
!
      MFULL = BUFR%MOL(I)
      M = MOD(BUFR%MOL(I),I_100)
!
!     ISO=(MOD(BUFR%MOL(I),I_1000)-M)/100   IS PROGRAMMED AS:
!
      ISO = MOD(BUFR%MOL(I),I_1000)/100
!
!     check if lines are within allowed molecular and isotopic limits
!
      if (m.gt.mxmol .or. m.lt. 1) then
         call line_exception (1,ipr,h_linf4,m,nmol,iso,iso_max)
         go to 50
      else if (iso .gt. iso_max(m)) then
         call line_exception (2,ipr,h_linf4,m,nmol,iso,iso_max)
         go to 50
      endif
!
      IF (ISOTPL_FLAG(M,ISO).EQ.0) THEN
         SUI = BUFR%SP(I)*W(M)
      ELSE
         SUI = BUFR%SP(I)*WKI(M,ISO)
      ENDIF
!
      IF (SUI.EQ.0.) GO TO 50
      IF (BUFR%VNU(I).LT.VLO) GO TO 50
      IJ = IJ+1
!
!     Y'S AND G'S ARE STORED IN I+1 POSTION OF VNU,S,ALFA0,EPP...
!      A(1-4),  B(1-4) CORRESPOND TO TEMPERATURES TEMPLC(1-4) ABOVE
!
      IF (IFLAG.EQ.1.OR.IFLAG.EQ.3) THEN
         A(1) = BUFR%VNU(I+1)
         B(1) = BUFR%SP(I+1)
         A(2) = BUFR%ALFA(I+1)
         B(2) = BUFR%EPP(I+1)
         A(3) = TRANSFER(BUFR%MOL(I+1),A(3) )  ! real representation of mol
         B(3) = BUFR%HWHM(I+1)
         A(4) = BUFR%TMPALF(I+1)
         B(4) = BUFR%PSHIFT(I+1)
!
!     CALCULATE SLOPE AND EVALUATE
!
         SLOPEY = (A(ILC)-A(ILC-1))*RECTLC
         SLOPEG = (B(ILC)-B(ILC-1))*RECTLC
         IF (IFLAG.EQ.1) THEN
            YI = A(ILC)+SLOPEY*TMPDIF
            GI = B(ILC)+SLOPEG*TMPDIF
         ELSE
            GAMMA1 = A(ILC)+SLOPEY*TMPDIF
            GAMMA2 = B(ILC)+SLOPEG*TMPDIF
         ENDIF
      ENDIF
!
!     IFLAG = 2 IS RESERVED FOR LINE COUPLING COEFFICIENTS ASSOCIATED
!               WITH AN EXACT TREATMENT (NUMERICAL DIAGONALIZATION)
!
!     IFLAG = 3 TREATS LINE COUPLING IN TERMS OF REDUCED WIDTHS
!
      SHRUNK%VNU(IJ) = BUFR%VNU(I) + RHORAT * BUFR%PSHIFT(I)
      SHRUNK%ALFA(IJ) = BUFR%ALFA(I)
      SHRUNK%EPP(IJ) = BUFR%EPP(I)
      SHRUNK%MOL(IJ) = M

      if(sum(bufr%brd_mol_flg(:,i)).gt.0.AND.ibrd.gt.0) then
         shrunk%vnu(ij) = shrunk%vnu(ij)+sum(rhoslf(1:mxbrdmol)*bufr%brd_mol_flg(:,i) &
         &           *(bufr%brd_mol_shft(:,i)-bufr%pshift(i)))
      endif

!
      IF (jrad4.EQ.1) SUI = SUI*SHRUNK%VNU(IJ)
!
      IF (SHRUNK%VNU(IJ).EQ.0.) SUI = 2.*SUI
!
!     TREAT TRANSITIONS WITH UNKNOWN EPP AS SPECIAL CASE
!
!>>   an epp value between -0.9999 and 0.  cm-1 is taken as valid
!
!>>   an epp value of -1. is assumed set by hitran indicating an unknown
!     value: no temperature correction is performed
!
!>>   for an epp value of less than -1., it is assumed that value has
!     been provided as a reasonable value to be used for purposes of
!     temperature correction.  epp is set positive
!
      if (shrunk%epp(ij).le.-1.001) shrunk%epp(ij) = abs(shrunk%epp(ij))

      if (shrunk%epp(ij).le.-0.999) MEFDP(M) = MEFDP(M)+1

!     temperature correction:

      if (shrunk%epp(ij) .gt. -0.999) then
         SUI = SUI*SCOR(m,iso)* EXP(-SHRUNK%EPP(ij)*BETACR)*  &
            (1.+EXP(-SHRUNK%VNU(ij)* BETA))
      endif
!
      SUMS = SUMS+SUI
!
!     TEMPERATURE CORRECTION OF THE HALFWIDTH
!     SELF TEMP DEPENDENCE TAKEN THE SAME AS FOREIGN
!
      TMPCOR = TRATIO**BUFR%TMPALF(I)
      ALFA0I = SHRUNK%ALFA(IJ)*TMPCOR
      HWHMSI = BUFR%HWHM(I)*TMPCOR
      SHRUNK%ALFA(IJ) = ALFA0I*(RHORAT-RHOSLF(m))+HWHMSI*RHOSLF(m)

      if(sum(bufr%brd_mol_flg(:,i)).gt.0.AND.ibrd.gt.0) then
         tmpcor_arr = tratio**bufr%brd_mol_tmp(:,i)
         alfa_tmp = bufr%brd_mol_hw(:,i)*tmpcor_arr
         alfsum = sum(rhoslf(1:mxbrdmol)*bufr%brd_mol_flg(:,i)*alfa_tmp)
         shrunk%alfa(ij) = (rhorat-sum(rhoslf(1:mxbrdmol)* &
         &           bufr%brd_mol_flg(:,i)))*alfa0i + alfsum
         if(bufr%brd_mol_flg(m,i).eq.0)   &
         &           shrunk%alfa(ij) = shrunk%alfa(ij) + rhoslf(m)*(hwhmsi-alfa0i)
      end if
!
      IF (IFLAG.EQ.3) SHRUNK%ALFA(IJ) = SHRUNK%ALFA(IJ)* &
      &      (1.0-GAMMA1*PAVP0-GAMMA2*PAVP2)
!
      SHRUNK%EPP(IJ) = SHRUNK%VNU(IJ)*ALFD1(m,iso)
      NLIN = NLIN+1
      SHRUNK%SP(IJ) = SUI*(1.+GI*PAVP2)
      SHRUNK%SPP(IJ) = SUI*YI*PAVP0
      SHRUNK%SRAD(IJ) = 0.0

!  For NLTE lines:
      FREQ=SHRUNK%VNU(IJ)
      RLOW=1.0
      RUPP=1.0

!
!     PICK OUT MOLECULAR TYPES WITH VIBRATIONAL STATES
!
      IF(NUMSTATE(M).GT.0) THEN
         NLOW=MOD(MFULL/1000,100)
         NUPP=MFULL/100000
!             DELTA=EXP(-FREQ/XKT)
!     xkt=tave/radcn2=1/beta
         DELTA=EXP(-FREQ*BETA)
         IF (NLOW.GT.NUMSTATE(M)) STOP 'NLOW GT NUMSTATE IN LINF4Q'
         INDLOW=NLOW + (ISO-1)*MAXSTATE
         IF (NLOW.GT.0) RLOW=RATSTATE(INDLOW,M)
         IF (NUPP.GT.NUMSTATE(M)) STOP 'NUPP GT NUMSTATE IN LINF4Q'
         INDUPP=NUPP + (ISO-1)*MAXSTATE
         IF (NUPP.GT.0) RUPP=RATSTATE(INDUPP,M)
         FNLTE=SHRUNK%SP(IJ)/(1.0-DELTA)
         SHRUNK%SP(IJ)=FNLTE*(RLOW-RUPP*DELTA)
         SHRUNK%SRAD(IJ)=FNLTE*(RLOW-RUPP)
      ELSE
         RLOW=0.
         RUPP=0.
      END IF
!
!     RLOW AND RUPP NOW SET
!
      IF (SHRUNK%VNU(IJ).GT.VHI) THEN
         IEOF = 1
         GO TO 60
      ENDIF

50 END DO
   IF (IJ.LT.LIMIN.AND.IEOF.EQ.0) THEN
      CALL CPUTIM (TIM2)
      TIMS = TIMS+TIM2-TIM1
      GO TO 30
   ENDIF
60 CALL CPUTIM (TIM2)
   IHI = IJ
   TIMS = TIMS+TIM2-TIM1
!
   CALL CPUTIM(TPAT0)

   CALL SHRINQ(SHRUNK)

   CALL CPUTIM(TPAT1)
   TSHRNK = TSHRNK + TPAT1-TPAT0
   IJ = ILO-1
   IF (IHI.LT.LIMIN.AND.IEOF.EQ.0) GO TO 30
!
   VNULO = SHRUNK%VNU(1)
   VNUHI = SHRUNK%VNU(IHI)
   JLIN = IHI
!
   IF (JLIN.GT.0) THEN
      CALL CPUTIM(TPAT0)
      CALL BUFOUT (LNFIL4,RCDHDR(1),NPHDRL)
      CALL BUFOUT (LNFIL4,SHRUNK,NLNGT4)
      CALL CPUTIM(TPAT1)
      TBUFFR = TBUFFR + TPAT1-TPAT0
   ENDIF
   NLINS = NLINS+IHI-IST+1
!
   IF (IEOF.EQ.1) GO TO 70
   IJ = 0
   ILO = 1
   GO TO 30
70 CONTINUE
!
   DO 80 M = 1, NMOL
      IF (MEFDP(M).GT.0) WRITE (IPR,905) MEFDP(M),M
80 END DO
   CALL CPUTIM (TIMEL1)
   TIME = TIMEL1-TIMEL0
   IF (NOPR.EQ.0) THEN
      WRITE (IPR,910) TIME,TIMR,TIMS,NLIN,NLINS
      L4TIM=TIME
      L4TMR=TIMR
      L4TMS=TIMS
      L4NLN=NLIN
      L4NLS=NLINS
      LOTHER = TSHRNK+TBUFFR+TMOLN4
   ENDIF
   RETURN
!
900 FORMAT ('0  *****  LINF4 - VNU LIMITS DO NOT INTERSECT WITH ',    &
   &        'LINFIL - LINF4 NOT USED *****',/,'   VNU = ',F10.3,      &
   &        ' - ',F10.3,' CM-1     LINFIL = ',F10.3,' - ',F10.3,      &
   &        ' CM-1')
905 FORMAT ('0*************************',I5,' STRENGTHS FOR',         &
   &        '  TRANSITIONS WITH UNKNOWN EPP FOR MOL =',I5,            &
   &        ' SET TO ZERO')
910 FORMAT ('0',20X,'TIME',11X,'READ',9X,'SHRINQ',6X,'NO. LINES',3X,  &
   &        'AFTER SHRINQ',/,2X,'LINF4 ',2X,3F15.3,2I15)
!
end subroutine LINF4Q
!___________________________________________________________________
!
SUBROUTINE SHRINQ(SHRUNK)
!
   USE struct_types
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE SHRINK COMBINES LINES FALLING IN A WAVENUMBER INTERVAL
!     DVR4/2 INTO A SINGLE EFFECTIVE LINE TO REDUCE COMPUTATION
!
   TYPE(LINE_SHRINK)  :: SHRUNK
!      COMMON VNU(1250),S(1250),ALFAL(1250),ALFAD(1250),MOL(1250),       &
!     &       SPP(1250),SRAD(1250)
   COMMON /R4SUB/ VLO,VHI,ILO,IST,IHI,LIMIN,LIMOUT,ILAST,DPTMN,      &
   &               DPTFC,ILIN4,ILIN4T
   COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)
!
   J = ILO-1
   DV = DVR4/2.
   VLMT = SHRUNK%VNU(ILO)+DV
!
!     INITIALIZE NON-CO2 SUMS
!
   SUMAL = 0.
   SUMAD = 0.
   SUMS = 0.
   SUMV = 0.
   SUMC = 0.
!
!     INITIALIZE CO2 SUMS
!
   SUMAL2 = 0.
   SUMAD2 = 0.
   SUMS2 = 0.
   SUMV2 = 0.
   SUMC2 = 0.
!
   DO 20 I = ILO, IHI

!     To prevent underflow issues in CONVF4 we set S < 1.0e-35 to zero
      IF (SHRUNK%SP(I).lt.1.0e-35) THEN
         SHRUNK%SP(I)= 0.0
         SHRUNK%SPP(I) =0.0
      ENDIF

!      skip shrink if desired
      IF (.NOT. USESHRINK) THEN
         J = J+1
         GO TO 10
      ENDIF
!
!     IF LINE COUPLING, DON'T SHRINK LINE
!
      IF (SHRUNK%SPP(I).NE.0.0) THEN
         J = J+1
         SHRUNK%VNU(J) = SHRUNK%VNU(I)
         SHRUNK%SP(J) = SHRUNK%SP(I)
         SHRUNK%ALFA(J) = SHRUNK%ALFA(I)
         SHRUNK%EPP(J) = SHRUNK%EPP(I)
         SHRUNK%SPP(J) = SHRUNK%SPP(I)
         SHRUNK%MOL(J) = SHRUNK%MOL(I)
         IF (SHRUNK%MOL(J).NE.2) SHRUNK%MOL(J) = 0
!
         GO TO 10
      ENDIF
!
!     NON-CO2 LINES OF MOLECULAR INDEX IT.NE.2   ARE LOADED
!     INTO SUMS IF THE FREQUENCY WITHIN DV GROUP
!
      IF (SHRUNK%MOL(I).NE.2) THEN
         SUMV = SUMV+SHRUNK%VNU(I)*SHRUNK%SP(I)
         SUMS = SUMS+SHRUNK%SP(I)
         SUMAL = SUMAL+SHRUNK%SP(I)*SHRUNK%ALFA(I)
         SUMAD = SUMAD+SHRUNK%SP(I)*SHRUNK%EPP(I)
         SUMC = SUMC+SHRUNK%SPP(I)
      ELSE
!
!     CO2 LINES LOADED     (MOL .EQ. 2)
!
         SUMV2 = SUMV2+SHRUNK%VNU(I)*SHRUNK%SP(I)
         SUMS2 = SUMS2+SHRUNK%SP(I)
         SUMAL2 = SUMAL2+SHRUNK%SP(I)*SHRUNK%ALFA(I)
         SUMAD2 = SUMAD2+SHRUNK%SP(I)*SHRUNK%EPP(I)
         SUMC2 = SUMC2+SHRUNK%SPP(I)
      ENDIF
!
!     IF LAST LINE OR VNU GREATER THAN LIMIT THEN STORE SUMS
!
10    IF (I.LT.IHI) THEN
         IF (SHRUNK%VNU(I+1).LE.VLMT) GO TO 20
      ENDIF
!
      VLMT = SHRUNK%VNU(I)+DV
!
!     ASSIGN NON-CO2 LINE AVERAGES TO 'GROUP' LINE J
!
      IF (SUMS.GT.0.) THEN
         J = J+1
         SHRUNK%SP(J) = SUMS
         SHRUNK%ALFA(J) = SUMAL/SUMS
         SHRUNK%EPP(J) = SUMAD/SUMS
         SHRUNK%VNU(J) = SUMV/SUMS
         SHRUNK%SPP(J) = SUMC
         SHRUNK%MOL(J) = 0
         SUMAL = 0.
         SUMAD = 0.
         SUMS = 0.
         SUMV = 0.
         SUMC = 0.
      ENDIF
!
!     ASSIGN CO2 LINE AVERAGES
!
      IF (SUMS2.GT.0.) THEN
         J = J+1
         SHRUNK%SP(J) = SUMS2
         SHRUNK%ALFA(J) = SUMAL2/SUMS2
         SHRUNK%EPP(J) = SUMAD2/SUMS2
         SHRUNK%VNU(J) = SUMV2/SUMS2
         SHRUNK%MOL(J) = 2
         SHRUNK%SPP(J) = SUMC2
         SUMAL2 = 0.
         SUMAD2 = 0.
         SUMS2 = 0.
         SUMV2 = 0.
         SUMC2 = 0.
      ENDIF
!
20 END DO
!
   ILO = J+1
   IHI = J
!
   RETURN
!
end subroutine SHRINQ
!_________________________________________________________________
!
SUBROUTINE LBLF4Q (JRAD,V1,V2)
!
   USE phys_consts, ONLY: radcn2
   USE struct_types
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE LBLF4 DOES A LINE BY LINE CALCULATION
!     USING FUNCTION F4.
!
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
!
   type(LINE_SHRINK)  ::  LINE
   COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
   &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
   &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
   &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
   &              EXTID(10)
   CHARACTER*8  EXTID
!
   CHARACTER*8      XID,       HMOLID,      YID
   Real*8               SEC   ,       XALTZ
!
   COMMON /FILHDR/ XID(10),SEC   ,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                W(60),PZL,PZU,TZL,TZU,WBROAD,DVO,V1H,V2H,TBOUND,  &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /XTIME/ TIME,TIMRDF,TIMCNV,TIMPNL,TF4,TF4RDF,TF4CNV,       &
   &               TF4PNL,TXS,TXSRDF,TXSCNV,TXSPNL
   COMMON /R4SUB/ VLO,VHI,ILO,IST,IHI,LIMIN,LIMOUT,ILAST,DPTMN,      &
   &               DPTFC,ILIN4,ILIN4T
   COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)
   COMMON /VOICOM/ AVRAT(102),DUMMY(5,102)
   COMMON /CONVF/ CHI(0:250),RDVCHI,RECPI,ZSQBND,A3,B3,JCNVF4
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIO,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
   EQUIVALENCE (IHIRAC,FSCDID(1)) , (ILBLF4,FSCDID(2))
!
   DATA JLBLF4 / 0 /
!
   CALL CPUTIM (TIM00)
!
   DPTMN = DPTMIN/RADFN(V2,TAVE/RADCN2)
   DPTFC = DPTFAC
   LIMIN = 1000
   LIMOUT = 2500
   JLBLF4 = 1
!
!     SET IEOF EQUAL TO -1 FOR FIRST READ
!
   IEOF = -1
!
   V1R4 = V1
   V2R4 = V2
   NPTR4 = (V2R4-V1R4)/DVR4+ONEPL
   NPTR4 = MIN(NPTR4,LIMOUT)
   V2R4 = V1R4+DVR4* REAL(NPTR4-1)
!
   LIMP2 = LIMOUT+2
   DO 10 I = 1, LIMP2
      R4(I) = 0.
      RR4(I) = 0.
10 END DO
   BETA = RADCN2/TAVE
   VLO = V1R4-BOUND4
   VHI = V2R4+BOUND4
20 CALL CPUTIM (TIM0)
   CALL RDLN4Q (LINE, IEOF)
   CALL CPUTIM (TIM1)
!
   IF (IEOF.EQ.2) THEN
      TF4 = TF4+TIM1-TIM00
      RETURN
   ENDIF
!
   TF4RDF = TF4RDF+TIM1-TIM0
   TIM2 = TIM1
   IF (IEOF.EQ.1.AND.IHI.EQ.0) GO TO 30
!
   CALL CNVF4Q (LINE)
!
   CALL CPUTIM (TIM3)
   TF4CNV = TF4CNV+TIM3-TIM2
!
!    IF IHI EQUALS -1 THEN END OF CONVOLUTION
!
   IF (IHI.EQ.-1) GO TO 30
   GO TO 20
!
30 CALL CPUTIM (TIM4)
!
   IF (JRAD.EQ.1) THEN
!
!     RADIATION FIELD
!
      XKT = 1./BETA
      VITST = V1R4-DVR4
      RDLAST = -1.
      NPTSI1 = 0
      NPTSI2 = 0
!
40    NPTSI1 = NPTSI2+1
!
      VI = V1R4+DVR4* REAL(NPTSI1-1)
      RADVI = RADFNI(VI,DVR4,XKT,VITST,RDEL,RDLAST)
!
      NPTSI2 = (VITST-V1R4)/DVR4+0.001
      NPTSI2 = MIN(NPTSI2,NPTR4)
!
      DO 50 I = NPTSI1, NPTSI2
         VI = VI+DVR4
         R4(I) = R4(I)*RADVI
         RR4(I) = RR4(I)*RADVI
         RADVI = RADVI+RDEL
50    CONTINUE
!
      IF (NPTSI2.LT.NPTR4) GO TO 40
   ENDIF
!
   CALL CPUTIM (TIM5)
   TF4PNL = TF4PNL+TIM5-TIM4
   TF4 = TF4+TIM5-TIM00
!
   RETURN
!
end subroutine LBLF4Q
SUBROUTINE RDLN4Q (LINE,IEOF)
!
   USE struct_types
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE RDLIN4Q INPUTS THE LINE DATA FROM LNFIL4
!
   CHARACTER*8      XID,       HMOLID,      YID
   Real*8               SEC   ,       XALTZ
!
   COMMON /FILHDR/ XID(10),SEC   ,PAV ,TAV ,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMIS ,FSCDID(17),NMOL,LAYRS ,YID1,YID(10),LSTWDF
   COMMON /R4SUB/ VLO,VHI,ILO,IST,IHI,LIMIN,LIMOUT,ILAST,DPTMN,      &
   &               DPTFC,ILIN4,ILIN4T
!
   type(LINE_SHRINK)  ::  LINE
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINDUM,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   common /eppinfo/ negepp_flag
   integer*4 negepp_flag

   Real DUM(2),LINPNL(2)
!
!      EQUIVALENCE (VMIN,LINPNL(1))
!
   IF (IEOF.EQ.-1) THEN
!
!     BUFFER PAST FILE HEADER
!
      REWIND LNFIL4
      ILIN4T = 0
      CALL BUFIN (LNFIL4,LEOF,DUM(1),1)
      IF (LEOF.EQ.0) STOP 'RDLIN4; TAPE9 EMPTY'
      IF (LEOF.EQ.-99) THEN
         IEOF = 2
!
         RETURN
!
      ENDIF
      if (negepp_flag.eq.1) CALL BUFIN (LNFIL4,LEOF,DUM(1),1)
   ENDIF
   IEOF = 0
   ILO = 1
   IHI = 0
!
10 READ(LNFIL4,END=20) VMIN,VMAX,NREC,NWDS
!   10 CALL BUFIN (LNFIL4,LEOF,LINPNL(1),NPHDRL)
!      IF (LEOF.EQ.0) GO TO 20
   ILIN4T = ILIN4T+NREC
   IF (VMAX.LT.VLO) THEN
      CALL BUFIN (LNFIL4,LEOF,DUM(1),1)
      GO TO 10
   ELSE
      CALL BUFIN (LNFIL4,LEOF,LINE,NWDS)
   ENDIF
   IHI = NREC
   IF (LINE%VNU(NREC).GT.VHI) GO TO 20
!
   RETURN
!
20 IEOF = 1
!
950 FORMAT (8a1)

   RETURN
!
end subroutine RDLN4Q

SUBROUTINE CNVF4Q (SHRUNK)  ! (VNU,SABS,ALFAL,ALFAD,MOL,SPP,SRAD)
!
   USE struct_types
   IMPLICIT REAL*8           (V)
   type(LINE_SHRINK)  :: SHRUNK
!
!     SUBROUTINE CNVF4Q CONVOLVES THE LINE DATA WITH FUNCTION F4
!
   CHARACTER*1 FREJ(1250),HREJ,HNOREJ
   COMMON /RCNTRL/ ILNFLG
   COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN
   COMMON /R4SUB/ VLO,VHI,ILO,IST,IHI,LIMIN,LIMOUT,ILAST,DPTMN,      &
   &               DPTFC,ILIN4,ILIN4T
   COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502)
   COMMON /VOICOM/ AVRAT(102),DUMMY(5,102)
   COMMON /CONVF/ CHI(0:250),RDVCHI,RECPI,ZSQBND,A3,B3,JCNVF4
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINDUM,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4,IBRD
!
   parameter (nzeta=101)
   real*8 a_1,a_2,a_3,b_1,b_2,b_3
   common /voigt_cf/                                                 &
   &     a_1(0:nzeta-1), a_2(0:nzeta-1), a_3(0:nzeta-1),              &
   &     b_1(0:nzeta-1), b_2(0:nzeta-1), b_3(0:nzeta-1)
!
   DATA ZBND / 64. /
   DATA HREJ /'0'/,HNOREJ /'1'/
!
   DATA I_1/1/, I_250/250/
!
   VNULST = V2R4+BOUND4
!
   IF (JCNVF4.NE.1234) then
!
      JCNVF4 = 1234
!
!        Obtain CHI SUB-LORENTZIAN FORM FACTOR FOR CARBON DIOXIDE:
!
      dvchi = 0.1
      rdvchi = 1./dvchi

      call chi_fn (chi,dvchi)
!
!        CONSTANTS FOR FOURTH FUNCTION LINE SHAPE:
!
      RECPI = 1./(2.*ASIN(1.))
      ZSQBND = ZBND*ZBND
      A3 = (1.+2.*ZSQBND)/(1.+ZSQBND)**2
      B3 = -1./(1.+ZSQBND)**2
!
   endif
!
   BNDSQ = BOUND4*BOUND4
   r_bndsq = 1./bndsq
!
!     START OF LOOP OVER LINES
!
   IF (ILNFLG.EQ.2) READ(16)(FREJ(I),I=ILO,IHI)
!
   DO 60 I = ILO, IHI
!
      IF (SHRUNK%SP(I).EQ.0..AND.SHRUNK%SPP(I).EQ.0.) GO TO 60
      ALFADI = SHRUNK%EPP(I)
      ALFALI = SHRUNK%ALFA(I)
      ZETAI = ALFALI/(ALFALI+ALFADI)
      IZ = 100.*ZETAI + ONEPL
      ZETDIF = 100.*ZETAI - REAL(IZ-1)
      ALFAVI = ( AVRAT(IZ) + ZETDIF*(AVRAT(IZ+1)-AVRAT(IZ)) ) *   &
         (ALFALI+ALFADI)
!
!        Interpolate coefficients to proper zeta
!        a_3, ... are zero indexed, other variables start at one
      izx=iz-1
      if (izx .eq. 101) then
         izx = 100
         zetdif = 1.0
      endif
      a3x=a_3(izx)+zetdif*(a_3(izx+1)-a_3(izx))
      b3x=b_3(izx)+zetdif*(b_3(izx+1)-b_3(izx))
      f4_64x= (a3x + b3x*zsqbnd)
!
      RALFVI = 1./ALFAVI
      SIV= SHRUNK%SP(I)*RALFVI
      siv_a = siv*a3x
      siv_b = siv*b3x
      SRV= SHRUNK%SRAD(I)*RALFVI
      srv_a = srv*a3x
      srv_b = srv*b3x

! nlte line coupling constant
      cupcon=shrunk%spp(i)/shrunk%sp(i)
!
      SPEAK = A3x*(ABS(SIV))
!
      JJ = (SHRUNK%VNU(I)-V1R4)/DVR4+1.
      JJ = MAX(JJ,I_1)
      JJ = MIN(JJ,NPTR4)
!
!     SPEAK is used for line rejection
      IF (ILNFLG.LE.1) THEN
         FREJ(I) = HNOREJ
!     No rejection for line-coupled lines (SHRUNK%SPP ne. 0)
         IF (SPEAK.LE.(DPTMN+DPTFC*R4(JJ)) .and. shrunk%spp(i).eq.0.)    &
            THEN
            FREJ(I) = HREJ
            GO TO 60
         ENDIF
      ELSE
         IF (FREJ(I).EQ.HREJ) GOTO 60
      ENDIF
!
      ILIN4 = ILIN4+1
!
      VNUI = SHRUNK%VNU(I)
!
30    CONTINUE
!
      XNUI = VNUI-V1R4
      JMIN = (XNUI-BOUND4)/DVR4+2.
!
      IF (VNUI.GE.VNULST) GO TO 70
      IF (JMIN.GT.NPTR4) GO TO 60
      JMIN = MAX(JMIN,I_1)
      JMAX = (XNUI+BOUND4)/DVR4+1.
      IF (JMAX.LT.JMIN) GO TO 50
      JMAX = MIN(JMAX,NPTR4)
      ALFLI2 = ALFALI*ALFALI
      ALFVI2 = ALFAVI*ALFAVI
      XJJ = REAL(JMIN-1)*DVR4
!
      siv_64 = siv*f4_64x*(alfli2 + alfvi2*zsqbnd)
      srv_64 = srv*f4_64x*(alfli2 + alfvi2*zsqbnd)
      f4bnd = siv_64/(ALFLI2+bndSQ)
      frbnd = srv_64/(ALFLI2+bndSQ)
!
!                FOURTH FUNCTION CONVOLUTION
!

      dptrat = shrunk%spp(i)/(shrunk%sp(i)*alfavi)
      dptrat_r = shrunk%spp(i)/(shrunk%srad(i)*alfavi)

      rec_alfvi2 = 1./ALFVI2
      siv_a3 = SIV*A3
      siv_b3 = SIV*B3
!
      IF (SHRUNK%MOL(I).EQ.2.) THEN

         DO 40 JJ = JMIN, JMAX
            XM = (XJJ-XNUI)
            XMSQ = XM*XM
            ZVSQ = XMSQ * rec_alfvi2
            fcnt_fn = (2.-(xmsq*r_bndsq))*f4bnd
            fcntr_fn= (2.-(xmsq*r_bndsq))*frbnd
!
            IF (ZVSQ.LE.ZSQBND) THEN
               F4FN = (siv_A3 + ZVSQ * siv_B3) - fcnt_fn
               F4FR = (srv_A3 + ZVSQ * srv_B3) - fcntr_fn
               IF (SHRUNK%SPP(I).NE.0.) THEN
                  F4FN = f4fn + xm*dptrat*f4fn
                  F4FR = F4FR + xm*dptrat_r*f4fr
               ENDIF
            ELSE
               F4FN = SIL/(ALFLI2+XMSQ) - fcnt_fn
               F4FR = SRL/(ALFLI2+XMSQ) - fcntr_fn
               IF (SHRUNK%SPP(I).NE.0.) THEN
                  F4FN = F4FN+XM*dptrat*f4fn
                  F4FR = F4FR+XM*dptrat_r*f4fr
               ENDIF
            ENDIF
!
            IF (SHRUNK%SPP(I).EQ.0.) THEN
!
!         ASSIGN ARGUMENT ISUBL OF THE FORM FACTOR FOR CO2 LINES
!
               ISUBL = RDVCHI*ABS(XM)+0.5
               ISUBL = MIN(ISUBL,i_250)
!
               R4(JJ) = R4(JJ)+F4FN*CHI(ISUBL)
               RR4(JJ) = RR4(JJ)+F4FR*CHI(ISUBL)
            ELSE
               R4(JJ) = R4(JJ)+F4FN
            ENDIF
!
!
            XJJ = XJJ+DVR4
40       CONTINUE
!
      ELSE
!
!        all molecules other than co2:

         DO 45 JJ = JMIN, JMAX
            XM = (XJJ-XNUI)
            XMSQ = XM*XM
            ZVSQ = XMSQ * rec_alfvi2
!
            IF (ZVSQ.LE.ZSQBND) THEN
               F4FN = (siv_A3 + ZVSQ * siv_B3) - F4BND
               F4FR = (srv_A3 + ZVSQ * srv_B3) - FRBND
               IF (SHRUNK%SPP(I).NE.0.) THEN
                  F4FN = f4fn + xm*dptrat*f4fn
                  F4FR = F4FR + xm*dptrat_r*f4fr
               ENDIF
            ELSE
               F4FN = SIL/(ALFLI2+XMSQ)-F4BND
               F4FR = SRL/(ALFLI2+XMSQ)-FRBND
               IF (SHRUNK%SPP(I).NE.0.) THEN
                  F4FN = F4FN+XM*dptrat*f4fn
                  F4FR = F4FR+XM*dptrat_r*f4fr
               ENDIF
            ENDIF
!
            IF (SHRUNK%SPP(I).EQ.0.) THEN
!
!         ASSIGN ARGUMENT ISUBL OF THE FORM FACTOR FOR CO2 LINES
!
               ISUBL = RDVCHI*ABS(XM)+0.5
               ISUBL = MIN(ISUBL,i_250)
!
               R4(JJ) = R4(JJ)+F4FN*CHI(ISUBL)
               RR4(JJ) = RR4(JJ)+F4FR*CHI(ISUBL)
            ELSE
               R4(JJ) = R4(JJ)+F4FN
            ENDIF
!
!
            XJJ = XJJ+DVR4
45       CONTINUE

      ENDIF
!
50    IF (VNUI.GT.0..AND.VNUI.LE.25.) THEN
!
!     THE CALCULATION FOR NEGATIVE VNU(I) IS FOR VAN VLECK WEISSKOPF
!
         VNUI = -SHRUNK%VNU(I)
!
         SHRUNK%SPP(I) = -SHRUNK%SPP(I)
!
         GO TO 30
!
      ENDIF
!
60 END DO
!
   IF (ILNFLG.EQ.1) WRITE(16)(FREJ(I),I=ILO,IHI)
   RETURN
!
!     IF END OF CONVOLUTION, SET IHI=-1 AND RETURN
!
70 CONTINUE
   IF (ILNFLG.EQ.1) WRITE(16)(FREJ(I),I=ILO,IHI)
   IHI = -1
!
   RETURN
!
end subroutine CNVF4Q
!---------------------
SUBROUTINE DEFNLTEDAT(NUMSTATE,IDSTATE,EESTATE,NDGSTATE,RATSTATE)

   USE lblparams
!      include 'lblparams.inc'
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
   CHARACTER*5 IDSTATE
   DIMENSION IDSTATE(MAXSTATE,MXMOL),NUMISO(MXMOL),                  &
   &     EESTATE(MAXSTATE*Max_ISO,MXMOL),                             &
   &     NDGSTATE(MAXSTATE*Max_ISO,MXMOL),NUMSTATE(MXMOL),            &
   &     RATSTATE(MAXSTATE*Max_ISO,MXMOL)
   PARAMETER (MAXH2O=8,MAXCO2=26,MAXO3=18,MAXCO=3,MAXNO=3)
   CHARACTER*5 AH2O,ACO2,AO3,ACO,ANO
   common /cmol_nam/ cmol(mxmol),cspc(mxspc)
   CHARACTER*6  CMOL,CSPC,TXTISO
   CHARACTER*80 ISOMOL(5)
   DATA ISOMOL/'H2O','CO2','O3','CO','NO'/
   DIMENSION JSOH2O(MAXH2O),AH2O(MAXH2O),FH2O(MAXH2O),MDGH2O(MAXH2O),&
   &          JSOCO2(MAXCO2),ACO2(MAXCO2),FCO2(MAXCO2),MDGCO2(MAXCO2),&
   &          JSOO3(MAXO3)  ,AO3 (MAXO3) ,FO3 (MAXO3) ,MDGO3 (MAXO3), &
   &          JSOCO(MAXCO)  ,ACO (MAXCO) ,FCO (MAXCO) ,MDGCO (MAXCO), &
   &          JSONO(MAXNO)  ,ANO (MAXNO) ,FNO (MAXNO) ,MDGNO (MAXNO)
!
   DATA  (JSOH2O(I),AH2O(I),FH2O(I),MDGH2O(I),I=1,8)/                &
   &     1, '000' ,     0.   , 1,                                     &
   &     1, '010' ,  1594.750, 1,                                     &
   &     1, '020' ,  3151.630, 1,                                     &
   &     1, '100' ,  3657.053, 1,                                     &
   &     1, '001' ,  3755.930, 1,                                     &
   &     1, '030' ,  4666.793, 1,                                     &
   &     1, '110' ,  5234.977, 1,                                     &
   &     1, '011' ,  5333.269, 1/
!
   DATA  (JSOCO2(I),ACO2(I),FCO2(I),MDGCO2(I),I=1, 9)/               &
   &     1, '00001' ,    0.   , 1 ,                                   &
   &     1, '01101' ,  667.380, 2 ,                                   &
   &     1, '10002' , 1285.409, 1 ,                                   &
   &     1, '02201' , 1335.132, 2 ,                                   &
   &     1, '10001' , 1388.185, 1 ,                                   &
   &     1, '11102' , 1932.470, 2 ,                                   &
   &     1, '03301' , 2003.246, 2 ,                                   &
   &     1, '11101' , 2076.856, 2 ,                                   &
   &     1, '00011' , 2349.143, 1 /
   DATA  (JSOCO2(I),ACO2(I),FCO2(I),MDGCO2(I),I=10,26)/              &
   &     1, '20003' , 2548.366, 1 ,                                   &
   &     1, '12202' , 2585.022, 2 ,                                   &
   &     1, '20002' , 2671.143, 1 ,                                   &
   &     1, '04401' , 2671.717, 2 ,                                   &
   &     1, '12201' , 2760.725, 2 ,                                   &
   &     1, '20001' , 2797.135, 1 ,                                   &
   &     1, '01111' , 3004.012, 2 ,                                   &
   &     1, '10012' , 3612.842, 1 ,                                   &
   &     1, '02211' , 3659.273, 2 ,                                   &
   &     1, '10011' , 3714.783, 1 ,                                   &
   &     1, '11112' , 4247.706, 2 ,                                   &
   &     1, '03311' , 4314.914, 2 ,                                   &
   &     1, '11111' , 4390.629, 2 ,                                   &
   &     1, '20013' , 4853.623, 1 ,                                   &
   &     1, '04411' , 4970.931, 1 ,                                   &
   &     1, '20012' , 4977.834, 1 ,                                   &
   &     1, '20011' , 5099.660, 1 /
!
   DATA (JSOO3(I),AO3(I),FO3(I),MDGO3(I),I=1,18)/                    &
   &     1, '000' ,    0.   , 1,                                      &
   &     1, '010' ,  700.931, 1,                                      &
   &     1, '001' , 1042.084, 1,                                      &
   &     1, '100' , 1103.140, 1,                                      &
   &     1, '020' , 1399.275, 1,                                      &
   &     1, '011' , 1726.528, 1,                                      &
   &     1, '110' , 1796.261, 1,                                      &
   &     1, '002' , 2057.892, 1,                                      &
   &     1, '101' , 2110.785, 1,                                      &
   &     1, '200' , 2201.157, 1,                                      &
   &     1, '111' , 2785.245, 1,                                      &
   &     1, '003' , 3041.200, 1,                                      &
   &     1, '004' , 3988.   , 1,                                      &
   &     1, '005' , 4910.   , 1,                                      &
   &     1, '006' , 5803.   , 1,                                      &
   &     1, '007' , 6665.   , 1,                                      &
   &     1, '008' , 7497.   , 1,                                      &
   &     1, '009' , 8299.   , 1/
!
   DATA (JSOCO(I),ACO(I),FCO(I),MDGCO(I),I=1,3)/                     &
   &     1, '0' ,    0.   , 1,                                        &
   &     1, '1' , 2143.272, 1,                                        &
   &     1, '2' , 4260.063, 1/
!
   DATA (JSONO(I),ANO(I),FNO(I),MDGNO(I),I=1,3)/                     &
   &     1, '0' ,    0.   , 1,                                        &
   &     1, '1' , 1878.077, 1,                                        &
   &     1, '2' , 3724.067, 1/
!
   DO J=1,MXMOL
      NUMSTATE(J)=0
      NUMISO(J)=0
      DO I=1,MAXSTATE
         IDSTATE(I,J)='     '
         DO ISO=1,Max_ISO
            INDEX=(ISO-1)*MAXSTATE + I
            EESTATE(INDEX,J)=0.
            NDGSTATE(INDEX,J)=0
            RATSTATE(INDEX,J)=1.
         END DO
      END DO
   END DO
!
!---H2O
   CALL GETINDEX(ISOMOL(1),CMOL,MXMOL,ID,TXTISO)
   DO I=1,MAXH2O
      ISOTOPE=JSOH2O(I)
      IF(ISOTOPE.LT.1 .OR. ISOTOPE.GT.Max_ISO)                       &
      &        STOP 'ERROR IN DEFNLTEDAT FOR ISOTOPE NUMBER'
      IF(ISOTOPE.GT.NUMISO(ID)) NUMISO(ID)=ISOTOPE
      IF(ISOTOPE.EQ.1) THEN
         NUMSTATE(ID) = NUMSTATE(ID)+1
         INDEX=NUMSTATE(ID)
         IDSTATE(INDEX,ID)=AH2O(I)
         ISO1=1
         ISO2=Max_ISO
      ELSE
         INDEX=0
         DO J=1,NUMSTATE(ID)
            IF(AH2O(I).EQ.IDSTATE(J,ID)) INDEX=J
         END DO
         IF(INDEX.EQ.0) STOP 'ERROR IN DEFNLTEDAT FOR H2O ISOTOPE>1'
         ISO1=ISOTOPE
         ISO2=ISOTOPE
      END IF
      DO J=ISO1,ISO2
         K=(J-1)*MAXSTATE + INDEX
         NDGSTATE(K,ID)=MDGH2O(I)
         EESTATE(K,ID)=FH2O(I)
      END DO
   END DO
!
!---CO2
   CALL GETINDEX(ISOMOL(2),CMOL,MXMOL,ID,TXTISO)
   DO I=1,MAXCO2
      ISOTOPE=JSOCO2(I)
      IF(ISOTOPE.LT.1 .OR. ISOTOPE.GT.Max_ISO)                       &
      &        STOP 'ERROR IN DEFNLTEDAT FOR ISOTOPE NUMBER'
      IF(ISOTOPE.GT.NUMISO(ID)) NUMISO(ID)=ISOTOPE
      IF(ISOTOPE.EQ.1) THEN
         NUMSTATE(ID) = NUMSTATE(ID)+1
         INDEX=NUMSTATE(ID)
         IDSTATE(INDEX,ID)=ACO2(I)
         ISO1=1
         ISO2=Max_ISO
      ELSE
         INDEX=0
         DO J=1,NUMSTATE(ID)
            IF(ACO2(I).EQ.IDSTATE(J,ID)) INDEX=J
         END DO
         IF(INDEX.EQ.0) STOP 'ERROR IN DEFNLTEDAT FOR CO2 ISOTOPE>1'
         ISO1=ISOTOPE
         ISO2=ISOTOPE
      END IF
      DO J=ISO1,ISO2
         K=(J-1)*MAXSTATE + INDEX
         NDGSTATE(K,ID)=MDGCO2(I)
         EESTATE(K,ID)=FCO2(I)
      END DO
   END DO
!
!---O3
   CALL GETINDEX(ISOMOL(3),CMOL,MXMOL,ID,TXTISO)
   DO I=1,MAXO3
      ISOTOPE=JSOO3(I)
      IF(ISOTOPE.LT.1 .OR. ISOTOPE.GT.Max_ISO)                       &
      &        STOP 'ERROR IN DEFNLTEDAT FOR ISOTOPE NUMBER'
      IF(ISOTOPE.GT.NUMISO(ID)) NUMISO(ID)=ISOTOPE
      IF(ISOTOPE.EQ.1) THEN
         NUMSTATE(ID) = NUMSTATE(ID)+1
         INDEX=NUMSTATE(ID)
         IDSTATE(INDEX,ID)=AO3(I)
         ISO1=1
         ISO2=Max_ISO
      ELSE
         INDEX=0
         DO J=1,NUMSTATE(ID)
            IF(AO3(I).EQ.IDSTATE(J,ID)) INDEX=J
         END DO
         IF(INDEX.EQ.0) STOP 'ERROR IN DEFNLTEDAT FOR O3 ISOTOPE>1'
         ISO1=ISOTOPE
         ISO2=ISOTOPE
      END IF
      DO J=ISO1,ISO2
         K=(J-1)*MAXSTATE + INDEX
         NDGSTATE(K,ID)=MDGO3(I)
         EESTATE(K,ID)=FO3(I)
      END DO
   END DO
!
!---CO
   CALL GETINDEX(ISOMOL(4),CMOL,MXMOL,ID,TXTISO)
   DO I=1,MAXCO
      ISOTOPE=JSOCO(I)
      IF(ISOTOPE.LT.1 .OR. ISOTOPE.GT.Max_ISO)                       &
      &        STOP 'ERROR IN DEFNLTEDAT FOR ISOTOPE NUMBER'
      IF(ISOTOPE.GT.NUMISO(ID)) NUMISO(ID)=ISOTOPE
      IF(ISOTOPE.EQ.1) THEN
         NUMSTATE(ID) = NUMSTATE(ID)+1
         INDEX=NUMSTATE(ID)
         IDSTATE(INDEX,ID)=ACO(I)
         ISO1=1
         ISO2=Max_ISO
      ELSE
         INDEX=0
         DO J=1,NUMSTATE(ID)
            IF(ACO(I).EQ.IDSTATE(J,ID)) INDEX=J
         END DO
         IF(INDEX.EQ.0) STOP 'ERROR IN DEFNLTEDAT FOR CO ISOTOPE>1'
         ISO1=ISOTOPE
         ISO2=ISOTOPE
      END IF
      DO J=ISO1,ISO2
         K=(J-1)*MAXSTATE + INDEX
         NDGSTATE(K,ID)=MDGCO(I)
         EESTATE(K,ID)=FCO(I)
      END DO
   END DO
!
!---NO
   CALL GETINDEX(ISOMOL(5),CMOL,MXMOL,ID,TXTISO)
   DO I=1,MAXNO
      ISOTOPE=JSONO(I)
      IF(ISOTOPE.LT.1 .OR. ISOTOPE.GT.Max_ISO)                       &
      &        STOP 'ERROR IN DEFNLTEDAT FOR ISOTOPE NUMBER'
      IF(ISOTOPE.GT.NUMISO(ID)) NUMISO(ID)=ISOTOPE
      IF(ISOTOPE.EQ.1) THEN
         NUMSTATE(ID) = NUMSTATE(ID)+1
         INDEX=NUMSTATE(ID)
         IDSTATE(INDEX,ID)=ANO(I)
         ISO1=1
         ISO2=Max_ISO
      ELSE
         INDEX=0
         DO J=1,NUMSTATE(ID)
            IF(ANO(I).EQ.IDSTATE(J,ID)) INDEX=J
         END DO
         IF(INDEX.EQ.0) STOP 'ERROR IN DEFNLTEDAT FOR NO ISOTOPE>1'
         ISO1=ISOTOPE
         ISO2=ISOTOPE
      END IF
      DO J=ISO1,ISO2
         K=(J-1)*MAXSTATE + INDEX
         NDGSTATE(K,ID)=MDGNO(I)
         EESTATE(K,ID)=FNO(I)
      END DO
   END DO
!
   DO I=1,MXMOL
      IF(NUMSTATE(I).GT.MAXSTATE) THEN
         WRITE(IPR,*) 'ERROR IN DEFNLTEDAT: MAXSTATE NEEDS TO BE '
         WRITE(IPR,*) 'INCREASED TO ',I,' TO ACCOMODATE DEFAULT ',      &
            'ISOTOPE STATE DATA'
         STOP 'ERROR IN DEFNLTEDAT'
      END IF
   END DO
!
   DO ID=1,MXMOL
      IF(NUMSTATE(ID).GT.0) THEN
         write(ipr,*) 'Defaults for Molecule',cmol(Id),'  id=',id,   &
         &           '  NUMSTATE= ',NUMSTATE(ID),'  NUMISO=',NUMISO(ID)
         write(ipr,901) (idstate(j,id),j=1,numstate(id))
901      format('IDSTATE= ',26a6)
         do iso=1,NUMISO(ID)
            j0=(iso-1)*maxstate
            write(ipr,902) (ndgstate(j+j0,id),j=1,numstate(id))
902         format('NDGSTATE= ',26I3)
         end do
         do iso=1,NUMISO(ID)
            j0=(iso-1)*maxstate
            write(ipr,903) (eestate(j+j0,id),j=1,numstate(id))
903         FORMAT('EESTATE=',26f9.3)
         end do
      END IF
   end do
!
   RETURN
end subroutine DEFNLTEDAT
!--------------------------------------
SUBROUTINE DROPSPACE(TXTIN,TXTOUT)
   CHARACTER*6 TXTIN,TXTOUT,BLANKS
   DATA BLANKS/'      '/
!
   DO I=1,6
      IF(TXTIN(I:I).NE.' ') THEN
         TXTOUT=TXTIN(I:6)//BLANKS
         RETURN
      END IF
   END DO
end subroutine DROPSPACE
