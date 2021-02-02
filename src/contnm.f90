!     path:      $HeadURL: https://svn.aer.com/svn/aer/project/RD/LBLRTM/trunk/src/contnm.f90 $
!     author:    $Author: jmascio $
!     revision:  $Revision: 32919 $
!     created:   $Date: 2019-08-23 09:54:07 -0400 (Fri, 23 Aug 2019) $
!
!  --------------------------------------------------------------------------
! |  Copyright �, Atmospheric and Environmental Research, Inc., 2017         |
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
SUBROUTINE CONTNM(JRAD)
!
   Use lblparams, ONLY: n_absrb, ipts, ipts2
   USE phys_consts, ONLY: radcn2
   IMPLICIT REAL*8           (V)
!
!     SUBROUTINE CONTNM CONTAINS THE CONTINUUM DATA
!     WHICH IS INTERPOLATED INTO THE ARRAY ABSRB
!
!********************************************
   COMMON /cnth2o/ V1h,V2h,DVh,NPTh,                                 &
   &                Ch(n_absrb),csh2o(n_absrb),cfh2o(n_absrb)
!********************************************
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   COMMON /XCONT/ V1C,V2C,DVC,NPTC,C(6000)
!
   CHARACTER*8      XID,       HMOLID,      YID
   REAL*8               SECANT,       XALTZ
!
   COMMON /CVRCNT/ HNAMCNT,HVRCNT
   COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
   &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
   &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4

   common /cntscl/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL
!
!------------------------------------
! for analytic derivative calculation
! note: ipts  = same dimension as ABSRB
!       ipts2 = same dimension as C
   common /CDERIV/ icflg,iuf,v1absc,v2absc,dvabsc,nptabsc,delT_pert, &
   &    dqh2oC(ipts),dTh2oC(ipts),dUh2o

   real cself(ipts),cfrgn_aj(ipts)
!------------------------------------
!
! for cloud calculation
! note: ipts  = same dimension as ABSRB
!
   DIMENSION c_cld(ipts)
!
!------------------------------------
!
   DIMENSION C0(n_absrb),C1(n_absrb),C2(n_absrb)
   DIMENSION SH2OT0(n_absrb),SH2OT1(n_absrb),FH2O(n_absrb),          &
             CN2T0(n_absrb),FCO2(n_absrb),CT1(n_absrb),CT2(n_absrb)
   DIMENSION CCH0(5150),CCH1(5150),CCH2(5150)
   DIMENSION CN0(5150),CN1(5150),CN2(5150)
!
   REAL ABSBSV(n_absrb)
!
   CHARACTER*18 HNAMCNT,HVRCNT
!
   equivalence (fscdid(4), iaersl)
!
   EQUIVALENCE (C0,SH2OT0,CN2T0,FCO2) , (C1,SH2OT1,CT1),             &
   &            (C2,FH2O,CT2)
!
   DATA P0 / 1013. /,T0 / 296. /
   DATA XLOSMT / 2.68675E+19 /
!
   DIMENSION XFACCO2(500)
!     Correction factors for CO2 from 2000-3000 cm-1 (mt_ckd_2.5)
!     (stored every 2 cm-1 - same as CO2 continuum).
   DATA XFACCO2/                                                     &
   &    1.0000,0.9998,0.9997,0.9996,0.9995,0.9994,0.9992,0.9991,      &
   &    0.9991,0.9990,0.9990,0.9989,0.9988,0.9988,0.9987,0.9986,      &
   &    0.9985,0.9984,0.9983,0.9982,0.9981,0.9980,0.9979,0.9978,      &
   &    0.9976,                                                       &
   &    0.9975,0.9973,0.9972,0.9970,0.9969,0.9967,0.9965,0.9963,      &
   &    0.9961,0.9958,0.9956,0.9954,0.9951,0.9948,0.9946,0.9943,      &
   &    0.9940,0.9936,0.9933,0.9929,0.9926,0.9922,0.9918,0.9913,      &
   &    0.9909,                                                       &
   &    0.9904,0.9899,0.9894,0.9889,0.9884,0.9878,0.9872,0.9866,      &
   &    0.9859,0.9853,0.9846,0.9838,0.9831,0.9823,0.9815,0.9806,      &
   &    0.9798,0.9789,0.9779,0.9770,0.9759,0.9749,0.9738,0.9727,      &
   &    0.9716,                                                       &
   &    0.9704,0.9691,0.9679,0.9666,0.9652,0.9638,0.9624,0.9609,      &
   &    0.9594,0.9578,0.9562,0.9546,0.9529,0.9511,0.9493,0.9475,      &
   &    0.9456,0.9436,0.9417,0.9396,0.9375,0.9354,0.9332,0.9310,      &
   &    0.9287,                                                       &
   &    0.9264,0.9240,0.9216,0.9191,0.9166,0.9140,0.9114,0.9087,      &
   &    0.9060,0.9032,0.9004,0.8976,0.8947,0.8917,0.8887,0.8857,      &
   &    0.8827,0.8796,0.8764,0.8732,0.8700,0.8668,0.8635,0.8602,      &
   &    0.8568,                                                       &
   &    0.8534,0.8500,0.8466,0.8432,0.8397,0.8362,0.8327,0.8292,      &
   &    0.8257,0.8221,0.8186,0.8151,0.8115,0.8080,0.8044,0.8009,      &
   &    0.7973,0.7938,0.7903,0.7868,0.7833,0.7799,0.7764,0.7730,      &
   &    0.7697,                                                       &
   &    0.7663,0.7630,0.7597,0.7565,0.7533,0.7502,0.7471,0.7441,      &
   &    0.7411,0.7382,0.7354,0.7326,0.7298,0.7272,0.7246,0.7221,      &
   &    0.7197,0.7173,0.7150,0.7129,0.7108,0.7088,0.7068,0.7050,      &
   &    0.7033,                                                       &
   &    0.7016,0.7001,0.6986,0.6973,0.6961,0.6949,0.6939,0.6930,      &
   &    0.6921,0.6914,0.6908,0.6903,0.6899,0.6897,0.6895,0.6895,      &
   &    0.6895,0.6895,0.6895,0.6895,0.6908,0.7014,0.7121,0.7227,      &
   &    0.7552,                                                       &
   &    0.8071,0.8400,0.9012,0.9542,1.0044,1.0330,1.0554,1.0766,      &
   &    1.0967,1.1160,1.1346,1.1525,1.1700,1.1869,1.2035,1.2196,      &
   &    1.2354,1.2509,1.2662,1.2811,1.2958,1.3103,1.3245,1.3386,      &
   &    1.3525,                                                       &
   &    1.3661,1.3796,1.3930,1.4062,1.4193,1.4322,1.4449,1.4576,      &
   &    1.4701,1.4825,1.4949,1.5070,1.5191,1.5311,1.5430,1.5548,      &
   &    1.5550,1.5550,1.5550,1.5550,1.5550,1.5550,1.5550,1.5550,      &
   &    1.5550,                                                       &
   &    1.5550,1.5550,1.5550,1.5550,1.5550,1.5550,1.5550,1.5550,      &
   &    1.5550,1.5550,1.5550,1.5549,1.5547,1.5543,1.5539,1.5532,      &
   &    1.5525,1.5516,1.5506,1.5494,1.5481,1.5467,1.5452,1.5435,      &
   &    1.5417,                                                       &
   &    1.5397,1.5377,1.5355,1.5332,1.5308,1.5282,1.5255,1.5228,      &
   &    1.5199,1.5169,1.5137,1.5105,1.5072,1.5037,1.5002,1.4966,      &
   &    1.4929,1.4890,1.4851,1.4811,1.4771,1.4729,1.4686,1.4643,      &
   &    1.4599,                                                       &
   &    1.4555,1.4509,1.4463,1.4417,1.4370,1.4322,1.4274,1.4225,      &
   &    1.4176,1.4126,1.4076,1.4025,1.3974,1.3923,1.3872,1.3820,      &
   &    1.3768,1.3716,1.3663,1.3611,1.3558,1.3505,1.3452,1.3400,      &
   &    1.3347,                                                       &
   &    1.3294,1.3241,1.3188,1.3135,1.3083,1.3030,1.2978,1.2926,      &
   &    1.2874,1.2822,1.2771,1.2720,1.2669,1.2618,1.2568,1.2518,      &
   &    1.2468,1.2419,1.2370,1.2322,1.2274,1.2227,1.2180,1.2133,      &
   &    1.2087,                                                       &
   &    1.2041,1.1996,1.1952,1.1907,1.1864,1.1821,1.1778,1.1737,      &
   &    1.1695,1.1654,1.1614,1.1575,1.1536,1.1497,1.1460,1.1422,      &
   &    1.1386,1.1350,1.1314,1.1280,1.1246,1.1212,1.1179,1.1147,      &
   &    1.1115,                                                       &
   &    1.1084,1.1053,1.1024,1.0994,1.0966,1.0938,1.0910,1.0883,      &
   &    1.0857,1.0831,1.0806,1.0781,1.0757,1.0734,1.0711,1.0688,      &
   &    1.0667,1.0645,1.0624,1.0604,1.0584,1.0565,1.0546,1.0528,      &
   &    1.0510,                                                       &
   &    1.0493,1.0476,1.0460,1.0444,1.0429,1.0414,1.0399,1.0385,      &
   &    1.0371,1.0358,1.0345,1.0332,1.0320,1.0308,1.0296,1.0285,      &
   &    1.0275,1.0264,1.0254,1.0244,1.0235,1.0226,1.0217,1.0208,      &
   &    1.0200,                                                       &
   &    1.0192,1.0184,1.0177,1.0170,1.0163,1.0156,1.0150,1.0143,      &
   &    1.0137,1.0132,1.0126,1.0121,1.0116,1.0111,1.0106,1.0101,      &
   &    1.0097,1.0092,1.0088,1.0084,1.0081,1.0077,1.0074,1.0070,      &
   &    1.0067,                                                       &
   &    1.0064,1.0061,1.0058,1.0055,1.0053,1.0050,1.0048,1.0046,      &
   &    1.0043,1.0041,1.0039,1.0037,1.0036,1.0034,1.0032,1.0030,      &
   &    1.0029,1.0027,1.0026,1.0025,1.0023,1.0022,1.0021,1.0020,      &
   &    1.0019,                                                       &
   &    1.0018,1.0017,1.0016,1.0015,1.0014,1.0014,1.0013,1.0012,      &
   &    1.0011,1.0011,1.0010,1.0010,1.0009,1.0009,1.0008,1.0007,      &
   &    1.0006,1.0005,1.0004,1.0003,1.0002,1.0001,1.0000,1.0000,      &
   &    1.0000/

   DIMENSION XFACREV(0:14),XFAC_RHU(-1:61)

!     Self correction factors for 820-960 cm-1.
   DATA (XFACREV(I),I=0,14)/                                         &
   &     1.003,1.009,1.015,1.023,1.029,1.033,                         &
   &     1.037,1.039,1.040,1.046,1.036,1.027,                         &
   &     1.01,1.002,1.00/
!
!     Foreign correction factors from joint RHUBC-II/RHUBC-I
!     analysis (mt_ckd_3.0).
!     Modified from new MWR analysis (Payne et al., 2021).
!     (mt_ckd_3.5)
   DATA (XFAC_RHU(I),I=-1,61)/                                       &
   &     0.7620,0.7840,                                               &
   &     0.7820,0.7840,0.7620,0.7410,0.7970,                          &
!   &     0.7810,0.8330,                                               &
!   &     0.8500,0.8330,0.7810,0.7540,0.8180,                          &
   &     0.9140,0.9980,0.9830,0.9330,0.8850,                          &
   &     0.8420,0.8070,0.8000,0.8010,0.8100,                          &
   &     0.8090,0.8320,0.8180,0.7970,0.8240,                          &
   &     0.8640,0.8830,0.8830,0.8470,0.8380,                          &
   &     0.8660,0.9410,1.0400,1.0680,1.1410,                          &
   &     1.0800,1.0340,1.1550,1.0990,1.0270,                          &
   &     0.9500,0.8950,0.8150,0.7830,0.7700,                          &
   &     0.7000,0.7650,0.7750,0.8500,0.9000,                          &
   &     0.9050,0.9540,1.0200,1.0200,1.0250,                          &
   &     1.0200,1.1000,1.1250,1.1200,1.1110,                          &
   &     1.1370,1.1600,1.1490,1.1070,1.0640,                          &
   &     1.0450/
!
!     ASSIGN SCCS VERSION NUMBER TO MODULE
!
!     Continuum calculation flags:
!     ---------------------------
!     ICNTNM Value      Self     Foreign    Rayleigh     Others
!           0            no        no          no          no
!           1            yes       yes         yes         yes
!           2            no        yes         yes         yes
!           3            yes       no          yes         yes
!           4            no        no          yes         yes
!           5            yes       yes         no          yes
!           6   READ IN XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN,
!               and XRAYL in Record 1.2a
!
!     ASSIGN CVS VERSION NUMBER TO MODULE
!
   HVRCNT = '$Revision: 32919 $'
!
   RHOAVE = (PAVE/P0)*(T0/TAVE)
   XKT = TAVE/RADCN2

!     the amagat value is used for the broadenening component for a
!     number of the collision induced continua
!
   amagat = (Pave/P0)*(273./Tave)
!
   WTOT = WBROAD
   DO 10 M = 1, NMOL
      WTOT = WTOT+WK(M)
10 CONTINUE

   x_vmr_h2o = wk(1)/wtot
   x_vmr_o2  = wk(7)/wtot
   x_vmr_n2  = 1. - x_vmr_h2o - x_vmr_o2

   wn2 = x_vmr_n2 * wtot

!     H2O continuum derivatives are computed w.r.t. ln(q)
!        dqh2o must be returned with the radiation field included

   if (icflg.gt.0) then

!     amounts except for species of interest and n2 have been set to
!     zero in lblrtm.
!     wn2 must be set to zero here:

      wn2 = 0.

!     zero derivative arrays and initialize panel information

      do j=1,ipts
         cself(j) = 0.0
         cfrgn_aj(j) = 0.0
      enddo

      v1absc=v1abs
      v2absc=v2abs
      dvabsc=dvabs
      nptabsc=nptabs
   endif
!
!=======================================================================
!
!**** CLOUD EFFECTIVE OPTICAL DEPTH  FROM "in_lblrtm_cld" file  ********
!=======================================================================
!
   if (iaersl.eq.5) then

      call cld_od (V1C,V2C,DVC,NPTC,c_cld,layer,xkt)

!        ---------------------------------------------------------
!        Radiation field
!
      if (jrad.eq.1) then

         do j = 1, nptc
            vj = v1c +real(j-1)*dvc
            c_cld(j) = c_cld(j)*RADFN(VJ,XKT)
         enddo

      endif
!        ---------------------------------------------------------

!        Interpolate to total optical depth grid

      CALL XINT (V1C,V2C,DVC,c_cld,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)

   endif
!
!=======================================================================

!               ********    WATER VAPOR   ********
!
!=======================================================================
!
   h2o_fac  = WK(1)/Wtot
   Rself    =     h2o_fac  * RHOave * 1.e-20 * xself
   Rfrgn    = (1.-h2o_fac) * RHOave * 1.e-20 * xfrgn
   Rfrgn_aj =     h2o_fac  * RHOave * 1.e-20 * xfrgn
!
!=======================================================================
!
!     CORRECTION TO THE WATER VAPOR CONTINUUM    mt_ckd_2.4   Nov 2008
!
!     The following modifications to the water vapor continuum arise
!     from new analyses of ARM measurements in the microwave and far-IR
!     regions. Analyses of measurements in the microwave are based
!     primarily on the two-channel MWR (23.8 and 31.4 GHz) at SGP,
!     with supporting evidence from 150 GHz MWRHF measurements during
!     the COPS campaign and from 170 GHz GVRP measurements at SGP (V. H.
!     Payne, E. J. Mlawer and S. A. Clough). Measurements in the far-IR
!     were from the AERO_ext at the NSA site, in the time surrounding
!     and including the RHUBC-I campaign (J. Delamere and S. A. Clough).
!
!=======================================================================
!
!                             SELF

!     Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1
!
   if ((V2.gt.-20.0).and.(V1.lt.20000.) .and. xself.gt.0.) then
      sh2ot0 = 0.
      sh2ot1 = 0.
!
      CALL SL296 (V1C,V2C,DVC,NPTC,SH2OT0,v1ss,v2ss)
      CALL SL260 (V1C,V2C,DVC,NPTC,SH2OT1,v1ss,v2ss)
!
!           Loop calculating self continuum optical depth
!
      TFAC = (TAVE-T0)/(260.-T0)

!  MT_CKD_3.5  All previous IR corrections now included in stored coefficients
!  rather than correction functions.

      DO 20 J = 1, NPTC
         VJ = V1C+DVC* REAL(J-1)
         SH2O = 0.
         IF (SH2OT0(J).GT.0.) THEN
            SH2O = SH2OT0(J)*(SH2OT1(J)/SH2OT0(J))**TFAC
         ENDIF
!              ---------------------------------------------------------
!
         cself(j) = WK(1)*(SH2O*Rself)
!
!********************************************
         v1h=V1C
         dvh=DVC
         npth=NPTC
!
         csh2o(j)=1.e-20 * sh2o * xself
!********************************************
!
!              ---------------------------------------------------------
!              Radiation field
!
         IF (JRAD.EQ.1) cself(j) = cself(j)*RADFN(VJ,XKT)
!              ---------------------------------------------------------

20    CONTINUE
!
!           Interpolate to total optical depth grid

      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,cself,1.0,V1ABS,DVABS,ABSRB,ist,last)

   endif
!
!=======================================================================
!                             FOREIGN
!
!--------------------------------------------------------------------
!
!        Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1
!
   if ((V2.gt.-20.0).and.(V1.lt.20000.) .and. xfrgn.gt.0.) then
      fh2o = 0.
!-----------------------------------------------------------------------
!           CORRECTION TO FOREIGN CONTINUUM   mt_ckd_2.4  Nov 2008   sac
!-----------------------------------------------------------------------

      f0     = 0.06
      V0F1   = 255.67
      HWSQ1  = 240.**2
      BETA1  = 57.83
      C_1    = -0.42
      N_1    = 8

      C_2    = 0.3
      BETA2  = 630.
      N_2    = 8
!-----------------------------------------------------------------------
!           mt_ckd_2.8    March 2016     Mlawer and Alvarado
!-----------------------------------------------------------------------
!     Extensive changes to the foreign continuum were made for
!     mt_ckd_2.8.  Based on measurements by Baranov and Lafferty
!     (2012) from 850-1160 cm-1 and Mondelain et al.(2014) at
!     4255 cm-1, a revised MT_CKD foreign continuum formulation
!     was derived and has been implemented in window regions >
!     4000 cm-1 (blended with previous coefficients in transition
!     regions between windows and bands). For 1800-3000 cm-1,
!     this formulation guided the spectral shape of the foreign
!     coefficients, but the values were reduced to obtain
!     agreement with measurements in this window by Baranov and
!     Lafferty (2012) and IASI measurements from 1900-2150 cm-1.
!     Coefficients in this region were derived simultaneously
!     with N2-H2O CIA coefficients and water vapor self continuum
!     coefficents.

!
      CALL FRN296 (V1C,V2C,DVC,NPTC,FH2O,v1ss,v2ss)
!
      DO 24 J = 1, NPTC
         VJ = V1C+DVC* REAL(J-1)
         IF (VJ .LE. 600.) THEN
            JFAC = (VJ +10.)/10. + 0.00001
            FSCAL = XFAC_RHU(JFAC)
         ELSE
!
            vdelsq1  = (VJ-V0F1)**2
            vdelmsq1 = (VJ+V0F1)**2
            VF1  = ((VJ-V0F1)/beta1)**N_1
            VmF1 = ((VJ+V0F1)/beta1)**N_1
            VF2  = ((VJ     )/beta2)**N_2

            FSCAL = 1. +                                          &
            &                (f0 + C_1*( (HWSQ1/(VDELSQ1 +HWSQ1+VF1))  +       &
            &                            (HWSQ1/(VDELmSQ1+HWSQ1+VmF1)) ) ) /   &
            &                                                 (1.+C_2*VF2)
         ENDIF

         FH2O(J)=FH2O(J)*FSCAL
!
         c_f = WK(1) * FH2O(J)
!
!********************************************
         cfh2o(j)=1.e-20 * fh2o(j) * xfrgn
!********************************************
!              ---------------------------------------------------------
!              Radiation field
!
         IF (JRAD.EQ.1) c_f = c_f * RADFN(VJ,XKT)
!              ---------------------------------------------------------

         C(J)        = c_f * RFRGN
         cfrgn_aj(j) = c_f * rfrgn_aj
!
24    CONTINUE
!
      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)
!
!           ------------------------------------------------------------
!
      if  (icflg.eq.1) then

         do j=1,nptc
            c(j) =  cself(j)-cfrgn_aj(j)
         enddo

         call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

         CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)

      endif

!           ------------------------------------------------------------
!
   endif

!=======================================================================


!     ********    CARBON DIOXIDE   ********
!
!
!        Only calculate if V2 > -20. cm-1 and V1 <  10000. cm-1
!
   if ((V2.gt.-20.0).and.(V1.lt.10000.) .and. xco2c.gt.0) then
      fco2 = 0.
!
      WCO2 = WK(2) * RHOAVE * 1.0E-20 * xco2c
!
      CALL FRNCO2 (V1C,V2C,DVC,NPTC,FCO2,tave,v1ss,v2ss)
!
      DO 30 J = 1, NPTC
         VJ = V1C+DVC* REAL(J-1)

!**mt_ckd_2.0      11 July 2007    sac
!             This continuum differs from mt_ck_1.3 in that an entirely
!             new co2 continuum has been developed based on the line
!             coupling parameters from Hartmann's group as distributed
!             with hitran.  This continuum must be used with lblrtm_v11
!             and spectral line parameters including Hartmann's line
!             parameters for co2.
!             Based on recent validation studies, a scaling of the
!             continuum for v3 is required to achieve an acceptable
!             result at 2385 cm-1, the 'bandhead' of v3.
!             Clough et al., presentation at EGU 2007
!   *** mt_ckd_2.5  Adjustment to the original scaling made.
!                   (temperature dependence of continuum also added)

         CFAC = 1.
         IF (VJ .GE. 2000. .AND. VJ .LE. 2998.) THEN
            JFAC = (VJ - 1998.)/2. + 0.00001
            CFAC = XFACCO2(JFAC)
         ENDIF
         fco2(j) = cfac*fco2(j)
!**********
!
         C(J) = FCO2(J)*WCO2

!
!              Radiation field
!
         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)

30    CONTINUE
      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)

   endif
!
!     ********    DIFFUSE OZONE  ********
!
!     Smoothing coefficients from 8920.0-9165.0 cm-1 and from
!     24570.0-24665.0 cm-1.  Data covers 9170.0-24565.0 cm-1
!     region.
!
   IF (V2.GT.8920.0.AND.V1.LE.24665.0.and.xo3cn.gt.0.) THEN
      cch0(:) = 0.
      cch1(:) = 0.
      cch2(:) = 0.

      WO3 = WK(3) * 1.0E-20 * xo3cn
      CALL XO3CHP (V1C,V2C,DVC,NPTO3,CCH0,CCH1,CCH2,v1ss,v2ss)
!
      DT=TAVE-273.15
!
      DO 50 J = 1, NPTO3
         CCH0(J)=(CCH0(J)+(CCH1(J)+CCH2(J)*DT)*DT)*WO3
         VJ = V1C+DVC* REAL(J-1)
         IF (JRAD.EQ.1) CCH0(J) = CCH0(J)*RADFN(VJ,XKT)
50    CONTINUE
      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,CCH0,1.0,V1ABS,DVABS,ABSRB,ist,last)
   ENDIF
!
   IF (V2.GT.27370..AND.V1.LT.40800. .and.xo3cn.gt.0.) THEN
      c0(:) = 0.
      ct1(:) = 0.
      ct2(:) = 0.

      WO3 = WK(3) * 1.E-20 * xo3cn
      TC = TAVE-273.15
      CALL O3HHT0 (V1C,V2C,DVC,NPTO3,C0,v1ss,v2ss)
      CALL O3HHT1 (V1T1,V2T1,DVT1,NPT1,CT1)
      CALL O3HHT2 (V1T2,V2T2,DVT2,NPT2,CT2)
!
      DO 60 J = 1, NPTO3
         C(J) = C0(J)*WO3
!
         VJ = V1C+DVC* REAL(J-1)
         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)
         C(J) = C(J)*(1.+CT1(J)*TC+CT2(J)*TC*TC)
60    CONTINUE
!
!           Save non-Hartley Huggins optical depth contribution to
!           prevent double counting for wavenumber region beyond
!           40800 cm-1.
!
      IF ((VJ.GT.40815.).AND.(V2.GT.40800) .and.xo3cn.gt.0.) THEN
         I_FIX = (40800.-V1ABS)/DVABS+1.001
         DO 62 I=I_FIX,NPTABS
            ABSBSV(I) = ABSRB(I)
62       CONTINUE
      ENDIF
!
!           Combine Hartley Huggins with previous optical depths
!
      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)
!
!           If V2 > 40800 cm-1, replace points with previously
!           saved values (non-Hartley Huggins contribution)
!
      IF ((VJ.GT.40815.).AND.(V2.GT.40800).and.xo3cn.gt.0.)THEN
         DO 64 I=I_FIX,NPTABS
            ABSRB(I) = ABSBSV(I)
64       CONTINUE
      ENDIF
   ENDIF
!
!        If V2 > 40800 cm-1, add UV Hartley Huggins contribution
!
   IF (V2.GT.40800..AND.V1.LT.54000. .and.xo3cn.gt.0.) THEN
      c0(:) = 0.

      WO3 = WK(3) * xo3cn
      CALL O3HHUV (V1C,V2C,DVC,NPTO3,C0,v1ss,v2ss)
!
      DO 70 J = 1, NPTO3
         C(J) = C0(J)*WO3
         VJ = V1C+DVC* REAL(J-1)
         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)
70    CONTINUE
!
!           Save non-Hartley Huggins UV optical depth contribution to
!           prevent double counting for wavenumber region before
!           40800 cm-1.
!
      IF (V1.LT.40800) THEN
         I_FIX = (40800.-V1ABS)/DVABS+1.001
         DO 72 I=1,I_FIX-1
            ABSBSV(I) = ABSRB(I)
72       CONTINUE
      ENDIF
!
!           Combine UV Hartley Huggins with previous optical depths
!
      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)
!
!           If V1 < 40800 cm-1, replace points with previously
!           saved values (non-Hartley Huggins UV contribution)
!
      IF (V1.LT.40800) THEN
         DO 74 I=1,I_FIX-1
            ABSRB(I) = ABSBSV(I)
74       CONTINUE
      ENDIF
!
   ENDIF
!
!     ********    O2 OXYGEN COLLISION INDUCED FUNDAMENTAL  ***********
!
!     version_1 of the Oxygen Collision Induced Fundamental
!
!     F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, J.-M. Hartmann,
!        and Ch. Boulet,
!        Infrared collision-induced absorption by O2 near 6.4 microns
!        for atmospheric applications: measurements and emprirical
!        modeling, Appl. Optics, 35, 5911-5917, (1996).

!
!        Only calculate if V2 > 1340. cm-1 and V1 <  1850. cm-1

   if ((V2.gt.1340.0).and.(V1.lt.1850.).and. xo2cn.gt.0.) then
      c0(:) = 0.
!
      tau_fac = xo2cn *  Wk(7) * 1.e-20 * amagat
!
!           Wk(7) is the oxygen column amount in units of molec/cm2
!           amagat is in units of amagats (air)
!
!           The temperature correction is done in the subroutine o2_ver_
!
      call o2_ver_1 (v1c,v2c,dvc,nptc,c0,tave,v1ss,v2ss)
!
!           c0 are the oxygen absorption coefficients at temperature
!           tave
!              - these absorption coefficients are in units of
!                   [(cm^2/molec) 10^20)]/(cm-1  amagat)
!              - cm-1 in the denominator arises through the removal
!                   of the radiation field
!              - for this case, an amagat is interpreted as one
!                   loschmidt of air (273K)
!
      DO 80 J = 1, NPTC
         VJ = V1C+DVC* REAL(J-1)
         C(J) = tau_fac * c0(J)
!
!              Radiation field
!
         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)
!
80    CONTINUE

      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)
   endif

!        ********    O2 Collision Induced   ********
!
!        O2 continuum formulated by Mate et al. over the spectral region
!        7550-8486 cm-1:  "Absolute Intensities for the O2 1.27 micron
!        continuum absorption", B. Mate, C. Lugez, G.T. Fraser, and
!        W.J. Lafferty, J. Geophys. Res., 104, 30,585-30,590, 1999.
!
!        Units of these coefficients are 1 / (amagat_O2*amagat_air)
!
!        Also, refer to the paper "Observed  Atmospheric
!        Collision Induced Absorption in Near Infrared Oxygen Bands",
!        Mlawer, Clough, Brown, Stephen, Landry, Goldman, & Murcray,
!        Journal of Geophysical Research (1998).
!
!        Only calculate if V2 > 7536. cm-1 and V1 <  8500. cm-1
!
   if ((V2.gt.7536.0).and.(V1.lt.8500.).and. xo2cn.gt.0.) then
      c0(:) = 0.
!
      a_o2  = 1./0.446
      a_n2  = 0.3/0.446
      a_h2o = 1.

      tau_fac = xo2cn * (Wk(7)/xlosmt) * amagat *                 &
      &           (a_o2*x_vmr_o2+a_n2*x_vmr_n2+a_h2o*x_vmr_h2o)

!
      CALL O2INF1 (V1C,V2C,DVC,NPTC,C0,v1ss,v2ss)
!
      DO 92 J = 1, NPTC
         C(J) = tau_fac * C0(J)
         VJ = V1C+DVC* REAL(J-1)
!
!              Radiation field
!
         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)

92    CONTINUE
!
      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)
!
   endif
!
!        O2 continuum formulated by Mlawer et al. over the spectral
!        region 9100-11000 cm-1. Refer to the paper "Observed
!        Atmospheric Collision Induced Absorption in Near Infrared
!        Oxygen Bands", Mlawer, Clough, Brown, Stephen, Landry, Goldman,
!        & Murcray, Journal of Geophysical Research (1998).
!
!        Only calculate if V2 > 9100. cm-1 and V1 <  11000. cm-1
!
   if ((V2.gt.9100.0).and.(V1.lt.11000.).and. xo2cn.gt.0.) then
      c0(:) = 0.
!
      CALL O2INF2 (V1C,V2C,DVC,NPTC,C0,v1ss,v2ss)
      WO2 = xo2cn * (WK(7)*1.e-20) * RHOAVE
      ADJWO2 = (WK(7)/WTOT) * (1./0.209) * WO2
!
      DO 93 J = 1, NPTC
         C(J) = C0(J)*ADJWO2
         VJ = V1C+DVC* REAL(J-1)
!
!              Radiation field
!
         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)
!
93    CONTINUE
!
      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)
!
   endif

! *******
!        O2 A-band continuum formulated by Mlawer based on solar FTS measurements.
!
!        Only calculate if V2 > 12961.5 cm-1 and V1 < 13221.5 cm-1
!
   if ((V2.gt.12961.5).and.(V1.lt.13221.5).and. xo2cn.gt.0.) then
      c0(:) = 0.
!
      tau_fac = xo2cn * (Wk(7)/xlosmt) * amagat
!
!
      CALL O2INF3 (V1C,V2C,DVC,NPTC,C0,v1ss,v2ss)
!
      DO 94 J = 1, NPTC
         C(J) = tau_fac * C0(J)
         VJ = V1C+DVC* REAL(J-1)
!
!              Radiation field
!
         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)
94    CONTINUE
!
      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)
!
   endif

! *******
!        O2 continuum formulated by Greenblatt et al. over the spectral
!        region 8797-29870 cm-1:  "Absorption Coefficients of Oxygen
!        Between 330 and 1140 nm, G.D. Greenblatt, J.J. Orlando, J.B.
!        Burkholder, and A.R. Ravishabkara,  J. Geophys. Res., 95,
!        18577-18582, 1990.
!
!        The units conversion to (cm^2/molec)/atm(o2)  has been done in
!        subroutine o2_vis.
!
!        Only calculate if V2 > 15000. cm-1 and V1 <  29870. cm-1
!
   if ((V2.gt.15000.0).and.(V1.lt.29870.).and. xo2cn.gt.0.) then
      c0(:) = 0.
!
      WO2 = WK(7) * 1.e-20 * ((pave/1013.)*(273./tave)) * xo2cn
      CHIO2 =  WK(7)/WTOT
      ADJWO2 = chio2 * WO2
!
      CALL O2_vis (V1C,V2C,DVC,NPTC,C0,v1ss,v2ss)
!
      DO 96 J = 1, NPTC
         C(J) = C0(J)*ADJWO2
         VJ = V1C+DVC* REAL(J-1)
!
!              Radiation field
!
         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)
!
96    CONTINUE
!
      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)
!
   endif
!
!        Only calculate if V2 > 36000. cm-1

   if (V2.gt.36000.0 .and. xo2cn.gt.0.) then
      c0(:) = 0.
      WO2 = WK(7) * 1.e-20 * xo2cn
!
      CALL O2HERZ (V1C,V2C,DVC,NPTC,C0,TAVE,PAVE,v1ss,v2ss)
      DO 90 J = 1, NPTC
         C(J) = C0(J)*WO2
         VJ = V1C+DVC* REAL(J-1)
!
!              Radiation field
!
         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)
90    CONTINUE
      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)

   endif
!
!  ******
!        O2 in far-UV (large portion is referred to as Schumann-Runge continuum)
!        Only calculate if V2 > 56740. cm-1.

   if (V2.gt.56740.0 .and. xo2cn.gt.0.) then
      c0(:) = 0.
      WO2 = WK(7) * 1.e-20 * xo2cn
!
      CALL O2FUV (V1C,V2C,DVC,NPTC,C0,v1ss,v2ss)
      DO 97 J = 1, NPTC
         VJ = V1C+DVC* REAL(J-1)
         C(J) = C0(J)*WO2
!
!              Radiation field
!
         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)
97    CONTINUE

      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)

   endif
!
!     *********************  NITROGEN CONTINUA  ********************
!
!
!     ******** NITROGEN COLLISION INDUCED PURE ROTATION BAND  ********
!
!        Model used:
!         Borysow, A, and L. Frommhold, "Collision-induced
!            rototranslational absorption spectra of N2-N2
!            pairs for temperatures from 50 to 300 K", The
!            Astrophysical Journal, 311, 1043-1057, 1986.
!
!     Uodated 2004/09/22 based on:
!
!      Boissoles, J., C. Boulet, R.H. Tipping, A. Brown and Q. Ma,
!         Theoretical Calculations of the Translation-Rotation
!         Collision-Induced Absorption in N2-N2, O2-O2 and N2-O2 Pairs,
!         J.Quant. Spec. Rad. Transfer, 82,505 (2003).
!
!         The temperature dependence between the two reference
!         temperatures has been assumed the same as that for the
!         original continuum.
!
!        THIS NITROGEN CONTINUUM IS IN UNITS OF 1./(CM AMAGAT^2)
!
!        Only calculate if V2 > -10. cm-1 and V1 <  350. cm-1
!
   if ((V2.gt.-10.0).and.(V1.lt.350.).and. xn2cn.gt.0.) then
      c0(:) = 0.
      c1(:) = 0.
!
!           The following puts WXN2 units in 1./(CM AMAGAT)
!
!           c1(j) represents the relative broadening efficiency of o2
!     a_h2o represents the relative broadening efficiency of h2o

      a_h2o = 1.

!     correct formulation for consistency with LBLRTM (per molec/cm^2)
!
      tau_fac =  xn2cn * (Wn2/xlosmt) * amagat
!
      CALL xn2_r (V1C,V2C,DVC,NPTC,c0,c1,Tave,v1ss,v2ss)
!
!           c1 is  ~ the ratio of alpha(n2-o2)/alpha(n2-n2)
!           Eq's 7 and 8 in the Boissoles paper.
!
      DO 40 J = 1, NPTC
         VJ = V1C+DVC* REAL(J-1)
!
         C(J) = tau_fac * c0(J) *                                 &
         &             (x_vmr_n2 + c1(j)*x_vmr_o2 + a_h2o*x_vmr_h2o)
!
!              Radiation field
!

         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)

40    CONTINUE

      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)

   endif
!
!
!        ********    NITROGEN COLLISION INDUCED FUNDAMENTAL ********
!
!        version_1 of the Nitrogen Collision Induced Fundamental
!
!        Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and
!        J._M. Hartmann, Infrared collision-induced absorption by
!        N2 near 4.3 microns for atmospheric applications:
!        Measurements and emprirical modeling, Appl. Optics, 35,
!        5911-5917, (1996).
!
!        Only calculate if V2 > 2001.77 cm-1 and V1 < 2897.59 cm-1.
!        Bounds for original implmentation were from Lafferty (2085-
!        2670 cm-1).  N2-N2 coefficients were analytically extended
!        to 2000-2900 cm-1 to allow implementation of N2-H2O, which
!        has been assumed to have a wider spectral shape to allow
!        agreement with  Baranov and Lafferty (2012).
!
   if ((V2.gt.2001.77).and.(V1.lt.2897.59).and.xn2cn.gt.0.) then
      cn0 = 0.
      cn1 = 0.
      cn2 = 0.
!
!           The absorption coefficients from the Lafferty et al.
!           reference are for pure nitrogen (absorber and broadener).
!
!     correct formulation for consistency with LBLRTM (per molec/cm^2)
!
      tau_fac =  xn2cn* (Wn2/xlosmt) * amagat
!
!           Wn2 is in units of molec/cm2
!           amagat is in units of amagats (air)
!
!           The temperature correction of the absorption coefficient and the
!           adjustments for relative broadening efficiency are done in
!           subroutine n2_ver_1:
!
      call n2_ver_1 (v1c,v2c,dvc,nptc,cn0,cn1,cn2,tave,v1ss,v2ss)
!
!           cn0 are the nitrogen absorption coefficients at
!           temperature tave
!              - these absorption coefficients are in units of
!                   [(cm^2/molec) 10^20)]/(cm-1  amagat)
!              - cm-1 in the denominator arises through the removal
!                   of the radiation field
!              - for this case, an amagat is interpreted as one
!                   loschmidt of air (273K)
!           cn1 are N2 absorption coefficents with collision partner O2
!           cn2 are N2 absorption coefficents with collision partner H2O
!
      DO 45 J = 1, NPTC
         VJ = V1C+DVC* REAL(J-1)
         C(J) = tau_fac * (x_vmr_n2*cn0(J) + x_vmr_o2*cn1(j) +    &
         &             x_vmr_h2o*cn2(j))
!              Radiation field
!
         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)

45    CONTINUE
!
      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)

   endif


!        ********  NITROGEN COLLISION INDUCED FIRST OVERTONE
!
!        version_1 of the Nitrogen Collision Induced First Overtone
!
!        Shapiro and Gush (1966) modified by Mlawer and Gombos (2015)
!        based on comparisons with measurements from SGP solar FTS
!        (TCCON network).
!
!        Only calculate if V2 > 4340. cm-1 and V1 <  4910. cm-1.
!
   if ((V2.gt.4340.0).and.(V1.lt.4910.).and. xn2cn.gt.0.) then
      c0(:) = 0.
!
!        All species are assumed to have the same broadening efficiency.
!           a_o2  represents the relative broadening efficiency of o2
!           a_h2o represents the relative broadening efficiency of h2o

      a_o2  = 1.
      a_h2o = 1.

!     correct formulation for consistency with LBLRTM (per molec/cm^2)
!
      tau_fac =  xn2cn* (Wn2/xlosmt) *                            &
      &           amagat * (x_vmr_n2+a_o2*x_vmr_o2+a_h2o*x_vmr_h2o)

!
!           Wn2 is in units of molec/cm2
!           amagat is in units of amagats (air)
!
!           The absorption coefficients are assumed to have no
!           temperature dependence.
!
      call n2_overtone1 (v1c,v2c,dvc,nptc,c0,v1ss,v2ss)
!
!           c0 are the nitrogen absorption coefficients
!              - these absorption coefficients are in units of
!                   [(cm^2/molec) 10^20)]/(cm-1  amagat)
!              - cm-1 in the denominator arises through the removal
!                   of the radiation field
!              - for this case, an amagat is interpreted as one
!                   loschmidt of air (273K)
!
      DO 48 J = 1, NPTC
         VJ = V1C+DVC* REAL(J-1)
         C(J) = tau_fac * c0(J)
!              Radiation field
!
         IF (JRAD.EQ.1) C(J) = C(J)*RADFN(VJ,XKT)

48    CONTINUE
!
      call pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

      CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,ist,last)
!            CALL XINT (V1C,V2C,DVC,C,1.0,V1ABS,DVABS,ABSRB,1,NPTABS)

   endif
!
!     ********** Rayleigh Scattering calculation **********
!     (sac)
!     The effects of Rayleigh scattering are also included in module
!     lbllow.f with the aerosol/cloud properties.  In the case that
!     that lbllow.f is selected (iaersl .ne. 0), the decision has been
!     made to include this effect with the scattering processes,
!     even though it is molecular in nature.  Otherwise the effects
!     of Rayleigh scattering are included here.
!
!     The formulation, adopted from MODTRAN_3.5 (using approximation
!     of Shettle et al., (Appl Opt, 2873-4, 1980) with depolarization
!     = 0.0279, output in km-1 for T=273K & P=1 ATM) is as follows:
!
!     The rayleigh extinction coefficient (scattering out of the direct
!     beam), ray_ext, can be defined as
!
!         ray_ext = (vrayleigh**4/(9.38076E18-1.08426E09*vrayleigh**2))
!     *        *wmol_tot*conv_cm2mol
!
!     where vrayleigh is the wavenumber value, wmol_tot is the total
!     molecular amount in the layer, and conv_cm2mol is the conversion
!     factor derived by multiplying air density (2.68675E19 mol/cm3)
!     at 273 K with the number of km per cm (1.e-5 km/cm).
!
!     The total layer amount of all molecules is calculated above as
!     WTOT. For numerical purposes a factor of 1.e-20  has been
!     included in conv_cm2mol and the same factor has been included
!     in the air density in the denominator as well.
!     In addition, a factor of 10,000 (cm-1) has been
!     divided out of vrayleigh. Finally, the radiation field is
!     excluded, so xvrayleigh**4 is replaced by xvrayleigh**3. When
!     JRAD=1, the radiation field is put in by multiplying the
!     absorption by xvrayleigh.
!
!     Rayleigh scattering in the direct beam is only calculated for
!     model runs > 820 cm-1.
!
   If ((iaersl.eq.0 .or. iaersl.eq.5).and. v2.ge.820.            &
   &       .and. xrayl.gt.0.) then
!
!        Thus the current formulation is
      conv_cm2mol = xrayl*1.E-20/(2.68675e-1*1.e5)
!
      do 95 i=1,nptabs
         vrayleigh = v1abs+(i-1)*dvabs
         xvrayleigh = vrayleigh/1.e4
         ray_ext = (xvrayleigh**3/(9.38076E2-10.8426*xvrayleigh**2))  &
         &                 *(wtot*conv_cm2mol)

!           Radiation field
         IF (JRAD.EQ.1) then
            ray_ext = ray_ext*xvrayleigh
         elseif (JRAD.eq. 0) then
!              ray_ext = ray_ext/1.e4
            ray_ext = ray_ext*xvrayleigh/radfn(vrayleigh,xkt)
         endif

         absrb(i) = absrb(i)+ray_ext

95    continue
!
   endif
!
100 continue

   RETURN
!
900 FORMAT (/,'0    *********************************************',/, &
   &          '     *      BYPASS O2 CONTINUUM TO HERZBERG      *',/, &
   &          '     *       AS A RESULT of TAVE > 350. K        *',/, &
   &          '     *********************************************',/)
!
end subroutine CONTNM
!
!     --------------------------------------------------------------

subroutine pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,last)

   IMPLICIT REAL*8           (V)

!   Set up needed variables for call to XINT
!   Output variables
!     v1abs_loc - wavenumber of first value to be processed in XINT
!     nptabs_loc - number of values to be processed in XINT
!     ist - index of first value to be processed in XINT

   nbnd_v1c =  2 +  (v1ss-v1abs)/dvabs + 1.e-5
   ist = max(1,nbnd_v1c)
   v1abs_loc = v1abs + dvabs * float(ist-1)

   nbnd_v2c = 1 + (v2ss-v1abs)/dvabs + 1.e-5
   last = min(nptabs,nbnd_v2c)

   return
end subroutine pre_xint

!     --------------------------------------------------------------

!
!
SUBROUTINE PRCNTM
!
!     THIS SUBROUTINE PRINTS THE CONTINUUM INFORMATION TO FILE IPR
!
   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NDFLE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
   COMMON /CNTPR/ CINFO1,CINFO2,cnam3,CINFO3,cnam4,CINFO4,CHEADING
!
   CHARACTER*18 cnam3(9),cnam4(42)
   CHARACTER*51 CINFO1(2,12),CINFO2(2,14),CINFO3(2,9),CINFO4(2,42)
   CHARACTER*40 CHEADING(3,2)
!
   WRITE (IPR,910) ((CINFO1(I,J),I=1,2),J=1,12)
   WRITE (IPR,910) ((CINFO2(I,J),I=1,2),J=1,14)
   WRITE (IPR,918) ((CHEADING(I,J),I=1,3),J=1,2)
   WRITE (IPR,915) (cnam3(j),(CINFO3(I,J),I=1,2),J=1,9)
   WRITE (IPR,915) (cnam4(j),(CINFO4(I,J),I=1,2),J=1,42)
!
   RETURN
!
910 FORMAT (18x,2A51)
915 FORMAT (a18,2A51)
918 FORMAT (3A40)
!

end subroutine PRCNTM
!
!     --------------------------------------------------------------
BLOCK DATA CNTINF
!
!     Continuum information for output to TAPE6 in SUBROUTINE PRCNTM
!
   COMMON /CNTPR/ CINFO1,CINFO2,CNAM3,CINFO3,CNAM4,CINFO4,CHEADING
   CHARACTER*18 cnam3(9),cnam4(42)
   CHARACTER*51 CINFO1(2,12),CINFO2(2,14),CINFO3(2,9),CINFO4(2,42)
   CHARACTER*40 CHEADING(3,2)
!
   DATA cnam3/                                                       &
   &     '                  ',                                        &
   &     '                  ',                                        &
   &     '                  ',                                        &
   &     ' ckd_1.0      2.2 ',                                        &
   &     '     "            ',                                        &
   &     ' ckd_2.1      3.3 ',                                        &
   &     '     "            ',                                        &
   &     ' ckd_2.2      3.7 ',                                        &
   &     '     "            '/
!           123456789-123456789-123456789-123456789-123456789-1
!
   DATA cnam4/                                                       &
   &     '     "            ',                                        &
   &     ' ckd_2.2.2    3.12',                                        &
   &     ' ckd_2.4.1    5.12',                                        &
   &     '     "            ',                                        &
   &     '     "            ',                                        &
   &     '     "            ',                                        &
   &     '     "            ',                                        &
   &     '     "            ',                                        &
   &     '     "            ',                                        &
   &     ' ckd_2.4.2    5.17',                                        &
   &     '     "            ',                                        &
   &     ' mt_ckd_1.00  7.01',                                        &
   &     '     "            ',                                        &
   &     ' mt_ckd_1.1   9.1 ',                                        &
   &     ' mt_ckd_1.2   9.2 ',                                        &
   &     ' mt_ckd_1.3  10.0 ',                                        &
   &     ' mt_ckd_2.0  11.1 ',                                        &
   &     ' mt_ckd_2.01 11.2 ',                                        &
   &     ' mt_ckd_2.1  11.3 ',                                        &
   &     '     "            ',                                        &
   &     ' mt_ckd_2.2  11.4 ',                                        &
   &     ' mt_ckd_2.3  11.5 ',                                        &
   &     ' mt_ckd_2.4  11.6 ',                                        &
   &     ' mt_ckd_2.5  11.7 ',                                        &
   &     '     "            ',                                        &
   &     '     "            ',                                        &
   &     ' mt_ckd_2.5.1 11.8',                                        &
   &     ' mt_ckd_2.5.2 12.0',                                        &
   &     ' mt_ckd_2.5.3  -  ',                                        &
   &     ' mt_ckd_2.7  12.4 ',                                        &
   &     '     "            ',                                        &
   &     ' mt_ckd_2.8  12.5 ',                                        &
   &     ' mt_ckd_3.0  12.6 ',                                        &
   &     ' mt_ckd_3.1  12.7 ',                                        &
   &     ' mt_ckd_3.2  12.8 ',                                        &
   &     ' mt_ckd_3.3  12.9 ',                                        &
   &     ' mt_ckd_3.4  12.10',                                        &
   &     ' mt_ckd_3.5  12.11',                                        &
   &     '     "            ',                                        &
   &     '     "            ',                                        &
   &     '                  ',                                        &
   &     '                  '/
!           123456789-123456789-123456789-123456789-123456789-1
!
   DATA CINFO1/                                                      &
   &     '                                                   ',       &
   &     '                                                   ',       &
   &     '                                                   ',       &
   &     '                                                   ',       &
   &     '*** CONTINUA mt_ckd_3.5                            ',       &
   &     '                                                   ',       &
   &     '                                                   ',       &
   &     '            Most recent significant change         ',       &
   &     '       H2O  SELF  (T)     0 - 20000 CM-1           ',       &
   &     '   mt_ckd_3.5 - modify Cs,Tdep < 800 cm-1(Oct 2020)',       &
   &     '            AIR           0 - 20000 CM-1           ',       &
   &     '   mt_ckd_3.5 - modify 0-40 cm-1         (Oct 2020)',       &
   &     '       CO2  AIR           0 - 10000 CM-1           ',       &
   &     '   mt_ckd_2.5 - modify 2000-3000 cm-1    (Jan 2010)',       &
   &     '            AIR   (T)  2386 -  2434 CM-1           ',       &
   &     '   mt_ckd_2.5 - temperature dep. added   (Jan 2010)',       &
   &     '       N2   SELF  (T)     0 -   350 CM-1           ',       &
   &     '   ckd_2.2    - Borysow Fromhold         (Feb 1996)',       &
   &     '            AIR   (T)     0 -   350 CM-1           ',       &
   &     '   mt_ckd_2.1 - incl. O2 rel. efficiency (Sep 2004)',       &
   &     '            AIR   (T)  1993 -  2905 CM-1           ',       &
   &     '   mt_ckd_2.8 - N2-H2O Mlawer/Alvarado   (Mar 2016)',       &
   &     '            AIR        4340 -  4910 CM-1           ',       &
   &     '   mt_ckd_2.7 - Mlawer and Gombos        (May 2015)' /
!           123456789-123456789-123456789-123456789-123456789-1
!
   DATA CINFO2/                                                      &
   &     '       O2   AIR   (T)  1340 -  1850 CM-1 ',                 &
   &     '   ckd_2.4.1  - Thibault et al.          (Mar 1998)',       &
   &     '            O2/N2      7550 -  8486 CM-1',                  &
   &     '   ckd_2.4.1  - Mate et al.              (Feb 2000)',       &
   &     '            AIR        9100 - 11000 CM-1',                  &
   &     '   ckd_2.4.1  - Mlawer et al.            (Aug 1999)',       &
   &     '            AIR       12950 - 13220 CM-1 A-band    ',       &
   &     '   mt_ckd_3.4 - Mlawer                   (May 2020)',       &
   &     '            O2        15000 - 29870 CM-1           ',       &
   &     '   ckd_2.4.2  - Greenblatt et al.        (May 2000)',       &
   &     '            O2/N2     36000 -  >>>> CM-1 HERZBERG  ',       &
   &     '   ckd_0                                           ',       &
   &     '            AIR       55000 - 87000 CM-1 FAR-UV    ',       &
   &     '   mt_ckd_3.3 - Lu et al.                (Aug 2018)',       &
   &     '       O3   AIR   (T)  9170 - 24565 CM-1 CHAPP/WULF',       &
   &     '   ckd_2.2                               (Feb 1996)',       &
   &     '                  (T) 27370 - 40800 CM-1 HARTL/HUGG',       &
   &     '   ckd_0                     ',                             &
   &     '                      40800 - 54000 CM-1 HARTL/HUGG',       &
   &     '   ckd_0                        ',                          &
   &     8*'                                                 '/
!
   DATA CHEADING/                                                    &
   &     'Continuum    LBL_n                      ',                  &
   &     'Description of Modification             ',                  &
   &     '                               Date     ',                  &
   &     '---------    -----                    --',                  &
   &     '-----------------------------           ',                  &
   &     '                               ----     '/
   DATA CINFO3/                                                      &
   &     '  H2O SELF HAS BEEN REDUCED IN THE 800-1200 CM-1 RE',       &
   &     'GION                                  (01 SEP 1985)',       &
   &     '  03       TEMPERATURE DEPENDENCE HAS BEEN CORRECTE',       &
   &     'D                                     (01 MAY 1987)',       &
   &     '  02       (1390-1760) HAS BEEN REDUCED (FACTOR = 0',       &
   &     '.78)                                  (07 MAR 1990)',       &
   &     '  H2O SELF HAS BEEN REDUCED IN THE 1100-1500 CM-1 R',       &
   &     'EGION                                 (01 APR 1993)',       &
   &     '  H2O FOREIGN HAS BEEN REDUCED AT ~1300 CM-1 AND IN',       &
   &     ' ALL THE WINDOW REGIONS               (01 APR 1993)',       &
   &     '  H2O SELF HAS BEEN MODIFIED IN THE 700-1500 CM-1 R',       &
   &     'EGION                                 (01 MAY 1994)',       &
   &     '  H2O FOREIGN HAS BEEN MODIFIED IN THE ENTIRE 1200-',       &
   &     '2200 CM-1 SPECTRAL RANGE              (01 MAY 1994)',       &
   &     '  H2O SELF HAS BEEN INCREASED 30% IN THE MICROWAVE ',       &
   &     'REGION                                (09 FEB 1996)',       &
   &     '  N2 COLLISION INDUCED PURE ROTATION BAND ADDED    ',       &
   &     '                                      (09 FEB 1996)'/

   DATA CINFO4/                                                      &
   &     '  O3 CHAPPUIS CHANGED TO VALUES FROM MODTRAN3      ',       &
   &     '                                      (09 FEB 1996)',       &
   &     '  INTERPOLATION EFFECTS BETWEEN HARTLEY HUGGINS DAT',       &
   &     'A AROUND 40800 CM-1                   (18 SEP 1996)',       &
   &     '  O2 COLLISION INDUCED FUNDAMENTAL BAND (1340-1850 ',       &
   &     'CM-1) CHANGED TO THIBAULT ET AL., 1996   (MAR 1998)',       &
   &     '  N2 COLLISION INDUCED FUNDAMENTAL BAND (2085-2670 ',       &
   &     'CM-1) CHANGED TO LAFFERTY ET AL., 1996   (MAR 1998)',       &
   &     '  H2O FOREIGN HAS BEEN MODIFIED IN THE 0-800 CM-1 A',       &
   &     'ND 1200-2200 CM-1 RANGE BASED ON      (04 JUN 1999)',       &
   &     '           AERI-ER FROM SHEBA, TOBIN ET AL., 1998  ',       &
   &     '                                                   ',       &
   &     '  H2O SELF HAS BEEN INCREASED IN THE 0-200 CM-1 REG',       &
   &     'ION                                   (04 JUN 1999)',       &
   &     '  O2 COLLISION INDUCED BAND HAS BEEN ADDED  (9100 -',       &
   &     ' 11000 CM-1); MLAWER ET AL. 1998      (16 AUG 1999)',       &
   &     '  O2 COLLISION INDUCED BAND HAS BEEN CHANGED (7555 ',       &
   &     '- 8486 CM-1); MATE ET AL. 1999        (10 FEB 2000)',       &
   &     '  O2 COLLISION INDUCED BANDS HAVE BEEN ADDED (15000',       &
   &     ' - 29870 CM-1); GREENBLATT ET AL. 1990 (8 MAY 2000)',       &
   &     '  The nu2 CO2 increased by a factor of 7: based on ',       &
   &     'U  Wisc. AFWEX and AERI_xr at ARM NSA    (Jul 2002)',       &
   &     '  *********   The entire water vapor continuum has ',       &
   &     'been revised   *********                 (Dec 2002)',       &
   &     '  *********   This continuum is based on a new line',       &
   &     ' shape/collision induced formulation, mt_ckd ******',       &
   &     '  H2O foreign modified in the 250-550 cm-1 region; ',       &
   &     'results now consistent with ckd_2.4.1    (Aug 2004)',       &
   &     '  Collision induced nitrogen 0-350 cm-1 increased (',       &
   &     '~1.35):  Boissoles at al., 2003          (Sep 2004)',       &
   &     '  Nu2 CO2: with P-R line mixing included; factor of',       &
   &     ' 7 increase has been reduceed to 4.5     (May 2006)',       &
   &     '  CO2: Based on Hartmann P-Q-R line mixing.  Modifi',       &
   &     'cation to v3 band based on AIRS data.    (Jul 2007)',       &
   &     '  H2O foreign modified in 250-550 cm-1 region based',       &
   &     ' on analyses of nsa aeri_er data.        (Sep 2007)',       &
   &     '  CO2: Fundamental change in the lblrtm fourth func',       &
   &     'tion with consequent changes in continuum(Nov 2007)',       &
   &     '  Bug fix impacting the nitrogen continuum in the 0',       &
   &     '-350 cm-1 region.                        (Nov 2007)',       &
   &     '  Analytic Derivative (species retrievals): n2 cont',       &
   &     'inuum removed from Jacobian calculation  (Mar 2008)',       &
   &     '  Bug fix: corrects error in which Analytic Derivat',       &
   &     'ive result depends on starting wavenumber(Aug 2008)',       &
   &     '  H2O: modification to self and foreign continuum (',       &
   &     'microwave and IR ARM data 0-600 cm-1)    (Nov 2008)',       &
   &     '  CO2: modification from 2000-3200 cm-1 (AERI(ARM), ',      &
   &     'IASI AIRS measurements); Temp. dep. added(Jan 2010)',       &
   &     '  H2O: modification to self cont. 2000-3000 cm-1   ',       &
   &     '(IASI data, fit to near-IR results of              ',       &
   &     '  Bicknell et al., 2006 and Fulghum and Tilleman,  ',       &
   &     '1991                                               ',       &
   &     '  Updated common block for constants               ',       &
   &     '                                         (May 2010)',       &
   &     '  Minor bug fixes, updated continuum summary info, ',       &
   &     ' made consistent with monortm contnm.f   (Jan 2011)',       &
   &     '  Cosmetic changes                                 ',       &
   &     '                                         (Mar 2015)',       &
   &     '  N2 1st overtone - Shapiro & Gush (1966) modified ',       &
   &     'by Mlawer and Gombos based on solar FTS  (Apr 2015)',       &
   &     '  Bug fix for analytic jacobians                   ',       &
   &     '                                         (Feb 2016)',       &
   &     '  N2-H2O from 2000-3000 cm-1,H2O self for 2000-3190',       &
   &     ' cm-1,H2O foreign in windows >1880 cm-1  (Jul 2016)',       &
   &     '  H2O: foreign/self 0-700/0-100 cm-1; RHUBC-II/-I; ',       &
   &     'REFIR,AERI,SAO-FTS; Mlawer,Turner,Paine  (Nov 2016)',       &
   &     '  O2 A-band: Mlawer and Gombos based on analysis of',       &
   &     'TCCON FTS measurements for OCO-2         (Jul 2017)',       &
   &     '  H2O self > 2000 cm-1 and tdep 1800-3500 cm-1: Cam',       &
   &     'pargue group and IASI obs; Mlawer et al. (Aug 2017)',       &
   &     '  O2 Far-UV: 55000-87000 cm-1 from Lu et al. (2010)',       &
   &     '                                         (Aug 2018)',       &
   &     '  O2 A-band: Mlawer re-analysis of TCCON FTS measur',       &
   &     'ements consistent with OCO-2 ABSCO 5.1   (May 2020)',       &
   &     '  H2O self: MW - Payne et al. 2021 MWR study, FIR -',       &
   &     ' Odintsova et al. 2020; self T: MW -     (Oct 2020)',       &
   &     '  Payne et al., Katkov et al. 1995, and Tretyakov e',       &
   &     't al. 2016, FIR - Odintsova et al. and             ',       &
   &     '  Burch/Grynvak 1979, IR - Burch/Alt 1984 and Burch',       &
   &     '/Grynvak; H2O foreign: MW - Payne et al.           ',       &
   &     '  -------------------------------------------------',       &
   &     '---------------------------------------------------',       &
   &     '  -------------------------------------------------',       &
   &     '---------------------------------------------------'/
!
end block data CNTINF
!
!
SUBROUTINE SL296 (V1C,V2C,DVC,NPTC,C,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   COMMON /SH2O/ V1S,V2S,DVS,NPTS,S(2003)
   DIMENSION C(*)
!
   DVC = DVS
   v1ss = v1s
   v2ss = v2s
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      C(J) = S(I)

10 END DO
!
   RETURN
!
end subroutine SL296
!
!     --------------------------------------------------------------
!
BLOCK DATA BS296
!
   IMPLICIT REAL*8           (V)
!
!               UNITS OF (CM**3/MOL)*1.E-20
!
   COMMON /SH2O/ V1,V2,DV,NPT,                                       &
   &              S0000( 2),S0001(50),S0051(50),S0101(50),S0151(50),  &
   &              S0201(50),S0251(50),S0301(50),S0351(50),S0401(50),  &
   &              S0451(50),S0501(50),S0551(50),S0601(50),S0651(50),  &
   &              S0701(50),S0751(50),S0801(50),S0851(50),S0901(50),  &
   &              S0951(50),S1001(50),S1051(50),S1101(50),S1151(50),  &
   &              S1201(50),S1251(50),S1301(50),S1351(50),S1401(50),  &
   &              S1451(50),S1501(50),S1551(50),S1601(50),S1651(50),  &
   &              S1701(50),S1751(50),S1801(50),S1851(50),S1901(50),  &
   &              S1951(50), S2001(1)
!
   DATA V1,V2,DV,NPT / -20.0, 20000.0, 10.0, 2003/
   DATA S0000/                                                       &
   &  2.731e-01, 2.855e-01/  
   DATA S0001/                                                       &
   &  2.877E-01,  2.855E-01,  2.731E-01,  2.490E-01,  2.178E-01,      &
   &  1.863E-01,  1.595E-01,  1.381E-01,  1.221E-01,  1.095E-01,      &
   &  9.882E-02,  8.926E-02,  8.041E-02,  7.217E-02,  6.450E-02,      &
   &  5.740E-02,  5.088E-02,  4.492E-02,  3.952E-02,  3.466E-02,      &
   &  3.031E-02,  2.644E-02,  2.302E-02,  2.001E-02,  1.738E-02,      &
   &  1.509E-02,  1.311E-02,  1.139E-02,  9.913E-03,  8.639E-03,      &
   &  7.544E-03,  6.603E-03,  5.794E-03,  5.096E-03,  4.495E-03,      &
   &  3.974E-03,  3.523E-03,  3.130E-03,  2.786E-03,  2.485E-03,      &
   &  2.221E-03,  1.987E-03,  1.780E-03,  1.596E-03,  1.433E-03,      &
   &  1.287E-03,  1.156E-03,  1.040E-03,  9.349E-04,  8.410E-04/   
   DATA S0051/                                                       &
   &  7.568E-04,  6.810E-04,  6.130E-04,  5.517E-04,  4.967E-04,      &
   &  4.471E-04,  4.025E-04,  3.624E-04,  3.262E-04,  2.937E-04,      &
   &  2.644E-04,  2.381E-04,  2.143E-04,  1.930E-04,  1.737E-04,      &
   &  1.564E-04,  1.413E-04,  1.284E-04,  1.176E-04,  1.087E-04,      &
   &  1.013E-04,  9.489E-05,  8.961E-05,  8.451E-05,  7.943E-05,      &
   &  7.430E-05,  6.912E-05,  6.447E-05,  6.025E-05,  5.642E-05,      &
   &  5.292E-05,  4.971E-05,  4.690E-05,  4.444E-05,  4.216E-05,      &
   &  4.012E-05,  3.814E-05,  3.621E-05,  3.442E-05,  3.268E-05,      &
   &  3.100E-05,  2.957E-05,  2.780E-05,  2.617E-05,  2.445E-05,      &
   &  2.305E-05,  2.187E-05,  2.080E-05,  1.980E-05,  1.885E-05/    
   DATA S0101/                                                       &
   &  1.796e-05,  1.711e-05,  1.633e-05,  1.559e-05,  1.490e-05,      &
   &  1.426e-05,  1.367e-05,  1.312e-05,  1.263e-05,  1.218e-05,      &
   &  1.178e-05,  1.143e-05,  1.112e-05,  1.088e-05,  1.070e-05,      &
   &  1.057e-05,  1.050e-05,  1.051e-05,  1.059e-05,  1.076e-05,      &
   &  1.100e-05,  1.133e-05,  1.180e-05,  1.237e-05,  1.308e-05,      &
   &  1.393e-05,  1.483e-05,  1.614e-05,  1.758e-05,  1.930e-05,      &
   &  2.123e-05,  2.346e-05,  2.647e-05,  2.930e-05,  3.279e-05,      &
   &  3.745e-05,  4.152e-05,  4.813e-05,  5.477e-05,  6.203e-05,      &
   &  7.331e-05,  8.056e-05,  9.882e-05,  1.050e-04,  1.210e-04,      &
   &  1.341e-04,  1.572e-04,  1.698e-04,  1.968e-04,  2.175e-04/
   DATA S0151/                                                       &
   &  2.431e-04,  2.735e-04,  2.867e-04,  3.190e-04,  3.371e-04,      &
   &  3.554e-04,  3.726e-04,  3.837e-04,  3.878e-04,  3.864e-04,      &
   &  3.858e-04,  3.841e-04,  3.852e-04,  3.815e-04,  3.762e-04,      &
   &  3.618e-04,  3.579e-04,  3.450e-04,  3.202e-04,  3.018e-04,      &
   &  2.785e-04,  2.602e-04,  2.416e-04,  2.097e-04,  1.939e-04,      &
   &  1.689e-04,  1.498e-04,  1.308e-04,  1.170e-04,  1.011e-04,      &
   &  9.237e-05,  7.909e-05,  7.006e-05,  6.112e-05,  5.401e-05,      &
   &  4.914e-05,  4.266e-05,  3.963e-05,  3.316e-05,  3.037e-05,      &
   &  2.598e-05,  2.294e-05,  2.066e-05,  1.813e-05,  1.583e-05,      &
   &  1.423e-05,  1.247e-05,  1.116e-05,  9.760e-06,  8.596e-06/
   DATA S0201/                                                       &
   &  7.720E-06,  7.091E-06,  6.499E-06,  5.801E-06,  5.150E-06,      &
   &  4.500E-06,  4.000E-06,  3.509E-06,  3.100E-06,  2.714E-06,      &
   &  2.406E-06,  2.108E-06,  1.875E-06,  1.669E-06,  1.555E-06,      &
   &  1.444E-06,  1.345E-06,  1.257E-06,  1.175E-06,  1.102E-06,      &
   &  1.035E-06,  9.745E-07,  9.194E-07,  8.685E-07,  8.222E-07,      &
   &  7.798E-07,  7.412E-07,  7.052E-07,  6.724E-07,  6.420E-07,      &
   &  6.143E-07,  5.886E-07,  5.649E-07,  5.430E-07,  5.228E-07,      &
   &  5.041E-07,  4.868E-07,  4.708E-07,  4.559E-07,  4.421E-07,      &
   &  4.293E-07,  4.174E-07,  4.064E-07,  3.961E-07,  3.865E-07,      &
   &  3.775E-07,  3.692E-07,  3.614E-07,  3.542E-07,  3.474E-07/
   DATA S0251/                                                       &
   &  3.411E-07,  3.353E-07,  3.298E-07,  3.248E-07,  3.201E-07,      &
   &  3.158E-07,  3.119E-07,  3.084E-07,  3.053E-07,  3.025E-07,      &
   &  3.003E-07,  2.988E-07,  2.976E-07,  2.977E-07,  2.989E-07,      &
   &  3.015E-07,  3.057E-07,  3.112E-07,  3.161E-07,  3.186E-07,      &
   &  3.176E-07,  3.142E-07,  3.099E-07,  3.065E-07,  3.011E-07,      &
   &  2.955E-07,  2.914E-07,  2.884E-07,  2.864E-07,  2.858E-07,      &
   &  2.851E-07,  2.857E-07,  2.870E-07,  2.888E-07,  2.914E-07,      &
   &  2.948E-07,  2.990E-07,  3.042E-07,  3.105E-07,  3.185E-07,      &
   &  3.276E-07,  3.385E-07,  3.511E-07,  3.646E-07,  3.858E-07,      &
   &  4.081E-07,  4.335E-07,  4.643E-07,  5.011E-07,  5.494E-07/
   DATA S0301/                                                       &
   &  5.953E-07,  6.354E-07,  7.116E-07,  7.498E-07,  8.613E-07,      &
   &  9.350E-07,  1.002E-06,  1.090E-06,  1.174E-06,  1.230E-06,      &
   &  1.288E-06,  1.334E-06,  1.383E-06,  1.380E-06,  1.390E-06,      &
   &  1.373E-06,  1.354E-06,  1.335E-06,  1.343E-06,  1.323E-06,      &
   &  1.349E-06,  1.353E-06,  1.362E-06,  1.344E-06,  1.329E-06,      &
   &  1.336E-06,  1.327E-06,  1.325E-06,  1.359E-06,  1.374E-06,      &
   &  1.415E-06,  1.462E-06,  1.526E-06,  1.619E-06,  1.735E-06,      &
   &  1.863E-06,  2.034E-06,  2.265E-06,  2.482E-06,  2.756E-06,      &
   &  3.103E-06,  3.466E-06,  3.832E-06,  4.378E-06,  4.913E-06,      &
   &  5.651E-06,  6.311E-06,  7.169E-06,  8.057E-06,  9.253E-06/
   DATA S0351/                                                       &
   &  1.047E-05,  1.212E-05,  1.360E-05,  1.569E-05,  1.776E-05,      &
   &  2.020E-05,  2.281E-05,  2.683E-05,  2.994E-05,  3.488E-05,      &
   &  3.896E-05,  4.499E-05,  5.175E-05,  6.035E-05,  6.340E-05,      &
   &  7.281E-05,  7.923E-05,  8.348E-05,  9.631E-05,  1.044E-04,      &
   &  1.102E-04,  1.176E-04,  1.244E-04,  1.283E-04,  1.326E-04,      &
   &  1.400E-04,  1.395E-04,  1.387E-04,  1.363E-04,  1.314E-04,      &
   &  1.241E-04,  1.228E-04,  1.148E-04,  1.086E-04,  1.018E-04,      &
   &  8.890E-05,  8.316E-05,  7.292E-05,  6.452E-05,  5.625E-05,      &
   &  5.045E-05,  4.380E-05,  3.762E-05,  3.290E-05,  2.836E-05,      &
   &  2.485E-05,  2.168E-05,  1.895E-05,  1.659E-05,  1.453E-05/
   DATA S0401/                                                       &
   &  1.282E-05,  1.132E-05,  1.001E-05,  8.836E-06,  7.804E-06,      &
   &  6.922E-06,  6.116E-06,  5.429E-06,  4.824E-06,  4.278E-06,      &
   &  3.788E-06,  3.371E-06,  2.985E-06,  2.649E-06,  2.357E-06,      &
   &  2.090E-06,  1.858E-06,  1.647E-06,  1.462E-06,  1.299E-06,      &
   &  1.155E-06,  1.028E-06,  9.142E-07,  8.132E-07,  7.246E-07,      &
   &  6.418E-07,  5.726E-07,  5.221E-07,  4.890E-07,  4.666E-07,      &
   &  4.443E-07,  4.203E-07,  3.978E-07,  3.770E-07,  3.575E-07,      &
   &  3.395E-07,  3.227E-07,  3.070E-07,  2.924E-07,  2.788E-07,      &
   &  2.661E-07,  2.543E-07,  2.433E-07,  2.329E-07,  2.233E-07,      &
   &  2.143E-07,  2.058E-07,  1.980E-07,  1.906E-07,  1.837E-07/
   DATA S0451/                                                       &
   &  1.773E-07,  1.713E-07,  1.657E-07,  1.605E-07,  1.557E-07,      &
   &  1.514E-07,  1.472E-07,  1.435E-07,  1.404E-07,  1.377E-07,      &
   &  1.354E-07,  1.332E-07,  1.308E-07,  1.282E-07,  1.256E-07,      &
   &  1.227E-07,  1.192E-07,  1.159E-07,  1.134E-07,  1.114E-07,      &
   &  1.098E-07,  1.081E-07,  1.066E-07,  1.055E-07,  1.045E-07,      &
   &  1.036E-07,  1.030E-07,  1.026E-07,  1.025E-07,  1.025E-07,      &
   &  1.028E-07,  1.033E-07,  1.041E-07,  1.049E-07,  1.060E-07,      &
   &  1.075E-07,  1.093E-07,  1.114E-07,  1.139E-07,  1.166E-07,      &
   &  1.199E-07,  1.239E-07,  1.282E-07,  1.316E-07,  1.355E-07,      &
   &  1.436E-07,  1.588E-07,  1.807E-07,  2.060E-07,  2.337E-07/
   DATA S0501/                                                       &
   &  2.645E-07,  2.996E-07,  3.393E-07,  3.843E-07,  4.363E-07,      &
   &  4.935E-07,  5.607E-07,  6.363E-07,  7.242E-07,  8.230E-07,      &
   &  9.411E-07,  1.071E-06,  1.232E-06,  1.402E-06,  1.600E-06,      &
   &  1.820E-06,  2.128E-06,  2.386E-06,  2.781E-06,  3.242E-06,      &
   &  3.653E-06,  4.323E-06,  4.747E-06,  5.321E-06,  5.919E-06,      &
   &  6.681E-06,  7.101E-06,  7.983E-06,  8.342E-06,  8.741E-06,      &
   &  9.431E-06,  9.952E-06,  1.026E-05,  1.055E-05,  1.095E-05,      &
   &  1.095E-05,  1.087E-05,  1.056E-05,  1.026E-05,  9.715E-06,      &
   &  9.252E-06,  8.452E-06,  7.958E-06,  7.268E-06,  6.295E-06,      &
   &  6.003E-06,  5.000E-06,  4.591E-06,  3.983E-06,  3.479E-06/
   DATA S0551/                                                       &
   &  3.058E-06,  2.667E-06,  2.293E-06,  1.995E-06,  1.747E-06,      &
   &  1.517E-06,  1.335E-06,  1.165E-06,  1.028E-06,  9.007E-07,      &
   &  7.956E-07,  7.015E-07,  6.192E-07,  5.491E-07,  4.859E-07,      &
   &  4.297E-07,  3.799E-07,  3.380E-07,  3.002E-07,  2.659E-07,      &
   &  2.366E-07,  2.103E-07,  1.861E-07,  1.655E-07,  1.469E-07,      &
   &  1.309E-07,  1.162E-07,  1.032E-07,  9.198E-08,  8.181E-08,      &
   &  7.294E-08,  6.516E-08,  5.787E-08,  5.163E-08,  4.612E-08,      &
   &  4.096E-08,  3.654E-08,  3.299E-08,  3.052E-08,  2.856E-08,      &
   &  2.683E-08,  2.507E-08,  2.346E-08,  2.194E-08,  2.053E-08,      &
   &  1.921E-08,  1.799E-08,  1.686E-08,  1.580E-08,  1.481E-08/
   DATA S0601/                                                       &
   &  1.389E-08,  1.304E-08,  1.226E-08,  1.151E-08,  1.083E-08,      &
   &  1.020E-08,  9.618E-09,  9.080E-09,  8.575E-09,  8.080E-09,      &
   &  7.624E-09,  7.176E-09,  6.761E-09,  6.370E-09,  6.004E-09,      &
   &  5.672E-09,  5.371E-09,  5.097E-09,  4.845E-09,  4.608E-09,      &
   &  4.390E-09,  4.192E-09,  4.011E-09,  3.846E-09,  3.698E-09,      &
   &  3.568E-09,  3.450E-09,  3.351E-09,  3.268E-09,  3.195E-09,      &
   &  3.131E-09,  3.091E-09,  3.060E-09,  3.039E-09,  3.034E-09,      &
   &  3.045E-09,  3.072E-09,  3.116E-09,  3.176E-09,  3.253E-09,      &
   &  3.347E-09,  3.467E-09,  3.605E-09,  3.763E-09,  3.946E-09,      &
   &  4.162E-09,  4.399E-09,  4.677E-09,  4.983E-09,  5.336E-09/
   DATA S0651/                                                       &
   &  5.728E-09,  6.168E-09,  6.684E-09,  7.214E-09,  7.861E-09,      &
   &  8.620E-09,  9.481E-09,  1.032E-08,  1.135E-08,  1.264E-08,      &
   &  1.396E-08,  1.567E-08,  1.740E-08,  1.977E-08,  2.200E-08,      &
   &  2.563E-08,  2.935E-08,  3.386E-08,  3.993E-08,  4.828E-08,      &
   &  5.714E-08,  6.877E-08,  8.135E-08,  9.779E-08,  1.132E-07,      &
   &  1.344E-07,  1.544E-07,  1.820E-07,  2.207E-07,  2.509E-07,      &
   &  2.854E-07,  3.026E-07,  3.278E-07,  3.474E-07,  3.693E-07,      &
   &  3.930E-07,  4.104E-07,  4.220E-07,  4.439E-07,  4.545E-07,      &
   &  4.778E-07,  4.812E-07,  5.018E-07,  4.899E-07,  5.075E-07,      &
   &  5.073E-07,  5.171E-07,  5.131E-07,  5.250E-07,  5.617E-07/
   DATA S0701/                                                       &
   &  5.846E-07,  6.239E-07,  6.696E-07,  7.398E-07,  8.073E-07,      &
   &  9.150E-07,  1.009E-06,  1.116E-06,  1.264E-06,  1.439E-06,      &
   &  1.644E-06,  1.856E-06,  2.147E-06,  2.317E-06,  2.713E-06,      &
   &  2.882E-06,  2.990E-06,  3.489E-06,  3.581E-06,  4.033E-06,      &
   &  4.260E-06,  4.543E-06,  4.840E-06,  4.826E-06,  5.013E-06,      &
   &  5.252E-06,  5.277E-06,  5.306E-06,  5.234E-06,  5.111E-06,      &
   &  5.134E-06,  4.748E-06,  4.434E-06,  4.109E-06,  3.605E-06,      &
   &  3.252E-06,  2.731E-06,  2.369E-06,  1.922E-06,  1.558E-06,      &
   &  1.321E-06,  1.058E-06,  8.515E-07,  6.811E-07,  5.431E-07,      &
   &  4.348E-07,  3.510E-07,  2.853E-07,  2.344E-07,  1.949E-07/
   DATA S0751/                                                       &
   &  1.649E-07,  1.409E-07,  1.232E-07,  1.098E-07,  9.866E-08,      &
   &  8.899E-08,  8.135E-08,  7.409E-08,  6.809E-08,  6.237E-08,      &
   &  5.728E-08,  5.267E-08,  4.865E-08,  4.482E-08,  4.169E-08,      &
   &  3.841E-08,  3.573E-08,  3.304E-08,  3.065E-08,  2.851E-08,      &
   &  2.634E-08,  2.448E-08,  2.275E-08,  2.114E-08,  1.965E-08,      &
   &  1.830E-08,  1.707E-08,  1.591E-08,  1.484E-08,  1.386E-08,      &
   &  1.294E-08,  1.209E-08,  1.131E-08,  1.056E-08,  9.901E-09,      &
   &  9.263E-09,  8.680E-09,  8.137E-09,  7.627E-09,  7.158E-09,      &
   &  6.715E-09,  6.305E-09,  5.927E-09,  5.568E-09,  5.242E-09,      &
   &  4.935E-09,  4.648E-09,  4.383E-09,  4.137E-09,  3.909E-09/
   DATA S0801/                                                       &
   &  3.697E-09,  3.500E-09,  3.321E-09,  3.156E-09,  3.006E-09,      &
   &  2.868E-09,  2.745E-09,  2.638E-09,  2.545E-09,  2.469E-09,      &
   &  2.405E-09,  2.369E-09,  2.344E-09,  2.350E-09,  2.371E-09,      &
   &  2.491E-09,  2.535E-09,  2.730E-09,  3.067E-09,  3.322E-09,      &
   &  3.630E-09,  4.045E-09,  4.524E-09,  4.981E-09,  5.613E-09,      &
   &  6.051E-09,  6.538E-09,  6.892E-09,  7.601E-09,  8.194E-09,      &
   &  8.806E-09,  9.450E-09,  9.745E-09,  1.009E-08,  1.017E-08,      &
   &  1.039E-08,  1.061E-08,  1.076E-08,  1.086E-08,  1.111E-08,      &
   &  1.142E-08,  1.160E-08,  1.174E-08,  1.180E-08,  1.187E-08,      &
   &  1.194E-08,  1.192E-08,  1.224E-08,  1.245E-08,  1.246E-08/
   DATA S0851/                                                       &
   &  1.318E-08,  1.377E-08,  1.471E-08,  1.582E-08,  1.713E-08,      &
   &  1.853E-08,  2.063E-08,  2.270E-08,  2.567E-08,  2.891E-08,      &
   &  3.264E-08,  3.744E-08,  4.286E-08,  4.915E-08,  5.623E-08,      &
   &  6.336E-08,  7.293E-08,  8.309E-08,  9.319E-08,  1.091E-07,      &
   &  1.243E-07,  1.348E-07,  1.449E-07,  1.620E-07,  1.846E-07,      &
   &  1.937E-07,  2.040E-07,  2.179E-07,  2.298E-07,  2.433E-07,      &
   &  2.439E-07,  2.464E-07,  2.611E-07,  2.617E-07,  2.582E-07,      &
   &  2.453E-07,  2.401E-07,  2.349E-07,  2.203E-07,  2.066E-07,      &
   &  1.939E-07,  1.780E-07,  1.558E-07,  1.391E-07,  1.203E-07,      &
   &  1.048E-07,  9.464E-08,  8.306E-08,  7.239E-08,  6.317E-08/
   DATA S0901/                                                       &
   &  5.520E-08,  4.847E-08,  4.282E-08,  3.796E-08,  3.377E-08,      &
   &  2.996E-08,  2.678E-08,  2.400E-08,  2.134E-08,  1.904E-08,      &
   &  1.705E-08,  1.523E-08,  1.350E-08,  1.204E-08,  1.070E-08,      &
   &  9.408E-09,  8.476E-09,  7.470E-09,  6.679E-09,  5.929E-09,      &
   &  5.267E-09,  4.711E-09,  4.172E-09,  3.761E-09,  3.288E-09,      &
   &  2.929E-09,  2.609E-09,  2.315E-09,  2.042E-09,  1.844E-09,      &
   &  1.640E-09,  1.470E-09,  1.310E-09,  1.176E-09,  1.049E-09,      &
   &  9.377E-10,  8.462E-10,  7.616E-10,  6.818E-10,  6.119E-10,      &
   &  5.558E-10,  5.152E-10,  4.832E-10,  4.546E-10,  4.267E-10,      &
   &  4.009E-10,  3.770E-10,  3.546E-10,  3.340E-10,  3.149E-10/
   DATA S0951/                                                       &
   &  2.971E-10,  2.805E-10,  2.652E-10,  2.509E-10,  2.377E-10,      &
   &  2.255E-10,  2.142E-10,  2.038E-10,  1.943E-10,  1.855E-10,      &
   &  1.776E-10,  1.703E-10,  1.640E-10,  1.586E-10,  1.536E-10,      &
   &  1.494E-10,  1.466E-10,  1.446E-10,  1.440E-10,  1.439E-10,      &
   &  1.455E-10,  1.461E-10,  1.479E-10,  1.514E-10,  1.538E-10,      &
   &  1.590E-10,  1.652E-10,  1.750E-10,  1.811E-10,  1.880E-10,      &
   &  1.933E-10,  2.027E-10,  2.078E-10,  2.170E-10,  2.249E-10,      &
   &  2.323E-10,  2.408E-10,  2.461E-10,  2.536E-10,  2.581E-10,      &
   &  2.615E-10,  2.721E-10,  2.782E-10,  2.866E-10,  2.938E-10,      &
   &  2.958E-10,  3.011E-10,  3.155E-10,  3.237E-10,  3.589E-10/
   DATA S1001/                                                       &
   &  3.838E-10,  4.083E-10,  4.354E-10,  4.745E-10,  5.148E-10,      &
   &  5.648E-10,  6.194E-10,  6.748E-10,  7.566E-10,  8.421E-10,      &
   &  9.427E-10,  1.042E-09,  1.165E-09,  1.341E-09,  1.518E-09,      &
   &  1.686E-09,  2.036E-09,  2.362E-09,  2.811E-09,  3.296E-09,      &
   &  4.070E-09,  4.695E-09,  6.152E-09,  7.572E-09,  9.557E-09,      &
   &  1.135E-08,  1.279E-08,  1.364E-08,  1.436E-08,  1.540E-08,      &
   &  1.672E-08,  1.793E-08,  1.906E-08,  2.036E-08,  2.144E-08,      &
   &  2.292E-08,  2.371E-08,  2.493E-08,  2.606E-08,  2.706E-08,      &
   &  2.866E-08,  3.036E-08,  3.136E-08,  3.405E-08,  3.665E-08,      &
   &  3.837E-08,  4.229E-08,  4.748E-08,  5.320E-08,  5.763E-08/
   DATA S1051/                                                       &
   &  6.677E-08,  7.216E-08,  7.716E-08,  8.958E-08,  9.419E-08,      &
   &  1.036E-07,  1.108E-07,  1.189E-07,  1.246E-07,  1.348E-07,      &
   &  1.310E-07,  1.361E-07,  1.364E-07,  1.363E-07,  1.343E-07,      &
   &  1.293E-07,  1.254E-07,  1.235E-07,  1.158E-07,  1.107E-07,      &
   &  9.961E-08,  9.011E-08,  7.910E-08,  6.916E-08,  6.338E-08,      &
   &  5.564E-08,  4.827E-08,  4.198E-08,  3.695E-08,  3.276E-08,      &
   &  2.929E-08,  2.633E-08,  2.391E-08,  2.192E-08,  2.021E-08,      &
   &  1.890E-08,  1.772E-08,  1.667E-08,  1.603E-08,  1.547E-08,      &
   &  1.537E-08,  1.492E-08,  1.515E-08,  1.479E-08,  1.450E-08,      &
   &  1.513E-08,  1.495E-08,  1.529E-08,  1.565E-08,  1.564E-08/
   DATA S1101/                                                       &
   &  1.553E-08,  1.569E-08,  1.584E-08,  1.570E-08,  1.538E-08,      &
   &  1.531E-08,  1.346E-08,  1.199E-08,  1.044E-08,  8.892E-09,      &
   &  7.076E-09,  5.884E-09,  5.039E-09,  4.129E-09,  3.663E-09,      &
   &  2.993E-09,  2.622E-09,  2.232E-09,  1.884E-09,  1.569E-09,      &
   &  1.348E-09,  1.144E-09,  1.009E-09,  8.833E-10,  7.907E-10,      &
   &  7.079E-10,  6.514E-10,  5.823E-10,  5.312E-10,  4.883E-10,      &
   &  4.458E-10,  4.094E-10,  3.734E-10,  3.458E-10,  3.205E-10,      &
   &  2.962E-10,  2.733E-10,  2.532E-10,  2.351E-10,  2.187E-10,      &
   &  2.034E-10,  1.889E-10,  1.764E-10,  1.648E-10,  1.541E-10,      &
   &  1.443E-10,  1.355E-10,  1.274E-10,  1.198E-10,  1.133E-10/
   DATA S1151/                                                       &
   &  1.071E-10,  1.015E-10,  9.648E-11,  9.195E-11,  8.791E-11,      &
   &  8.435E-11,  8.125E-11,  7.859E-11,  7.621E-11,  7.481E-11,      &
   &  7.351E-11,  7.228E-11,  7.295E-11,  7.339E-11,  7.588E-11,      &
   &  7.701E-11,  7.886E-11,  8.511E-11,  9.061E-11,  9.926E-11,      &
   &  1.083E-10,  1.241E-10,  1.469E-10,  1.736E-10,  1.966E-10,      &
   &  2.132E-10,  2.280E-10,  2.473E-10,  2.718E-10,  2.922E-10,      &
   &  3.128E-10,  3.361E-10,  3.641E-10,  3.910E-10,  4.196E-10,      &
   &  4.501E-10,  4.932E-10,  5.258E-10,  5.755E-10,  6.253E-10,      &
   &  6.664E-10,  7.344E-10,  7.985E-10,  8.877E-10,  1.005E-09,      &
   &  1.118E-09,  1.251E-09,  1.428E-09,  1.610E-09,  1.888E-09/
   DATA S1201/                                                       &
   &  2.077E-09,  2.331E-09,  2.751E-09,  3.061E-09,  3.522E-09,      &
   &  3.805E-09,  4.181E-09,  4.575E-09,  5.167E-09,  5.634E-09,      &
   &  6.007E-09,  6.501E-09,  6.829E-09,  7.211E-09,  7.262E-09,      &
   &  7.696E-09,  7.832E-09,  7.799E-09,  7.651E-09,  7.304E-09,      &
   &  7.150E-09,  6.977E-09,  6.603E-09,  6.209E-09,  5.690E-09,      &
   &  5.432E-09,  4.764E-09,  4.189E-09,  3.640E-09,  3.203E-09,      &
   &  2.848E-09,  2.510E-09,  2.194E-09,  1.946E-09,  1.750E-09,      &
   &  1.567E-09,  1.426E-09,  1.302E-09,  1.197E-09,  1.109E-09,      &
   &  1.035E-09,  9.719E-10,  9.207E-10,  8.957E-10,  8.578E-10,      &
   &  8.262E-10,  8.117E-10,  7.987E-10,  7.875E-10,  7.741E-10/
   DATA S1251/                                                       &
   &  7.762E-10,  7.537E-10,  7.424E-10,  7.474E-10,  7.294E-10,      &
   &  7.216E-10,  7.233E-10,  7.140E-10,  7.085E-10,  6.434E-10,      &
   &  5.585E-10,  5.015E-10,  3.956E-10,  3.337E-10,  2.545E-10,      &
   &  2.141E-10,  1.761E-10,  1.685E-10,  1.507E-10,  1.241E-10,      &
   &  1.081E-10,  9.675E-11,  7.964E-11,  6.930E-11,  5.790E-11,      &
   &  5.069E-11,  4.452E-11,  3.970E-11,  3.551E-11,  3.228E-11,      &
   &  2.939E-11,  2.682E-11,  2.471E-11,  2.278E-11,  2.106E-11,      &
   &  1.950E-11,  1.811E-11,  1.682E-11,  1.563E-11,  1.458E-11,      &
   &  1.363E-11,  1.278E-11,  1.200E-11,  1.129E-11,  1.066E-11,      &
   &  1.007E-11,  9.545E-12,  9.075E-12,  8.656E-12,  8.285E-12/
   DATA S1301/                                                       &
   &  7.958E-12,  7.674E-12,  7.431E-12,  7.227E-12,  7.061E-12,      &
   &  6.934E-12,  6.843E-12,  6.790E-12,  6.775E-12,  6.798E-12,      &
   &  6.859E-12,  6.960E-12,  7.102E-12,  7.287E-12,  7.517E-12,      &
   &  7.795E-12,  8.123E-12,  8.504E-12,  8.945E-12,  9.450E-12,      &
   &  1.003E-11,  1.070E-11,  1.148E-11,  1.232E-11,  1.327E-11,      &
   &  1.441E-11,  1.535E-11,  1.655E-11,  1.814E-11,  2.020E-11,      &
   &  2.307E-11,  2.622E-11,  2.962E-11,  3.369E-11,  3.819E-11,      &
   &  4.329E-11,  4.932E-11,  5.589E-11,  6.364E-11,  7.284E-11,      &
   &  8.236E-11,  9.447E-11,  1.078E-10,  1.229E-10,  1.417E-10,      &
   &  1.614E-10,  1.843E-10,  2.107E-10,  2.406E-10,  2.728E-10/
   DATA S1351/                                                       &
   &  3.195E-10,  3.595E-10,  4.153E-10,  4.736E-10,  5.410E-10,      &
   &  6.088E-10,  6.769E-10,  7.691E-10,  8.545E-10,  9.621E-10,      &
   &  1.047E-09,  1.161E-09,  1.296E-09,  1.424E-09,  1.576E-09,      &
   &  1.739E-09,  1.893E-09,  2.080E-09,  2.336E-09,  2.604E-09,      &
   &  2.760E-09,  3.001E-09,  3.365E-09,  3.550E-09,  3.895E-09,      &
   &  4.183E-09,  4.614E-09,  4.846E-09,  5.068E-09,  5.427E-09,      &
   &  5.541E-09,  5.864E-09,  5.997E-09,  5.997E-09,  6.061E-09,      &
   &  5.944E-09,  5.855E-09,  5.661E-09,  5.523E-09,  5.374E-09,      &
   &  4.940E-09,  4.688E-09,  4.170E-09,  3.913E-09,  3.423E-09,      &
   &  2.997E-09,  2.598E-09,  2.253E-09,  1.946E-09,  1.710E-09/
   DATA S1401/                                                       &
   &  1.507E-09,  1.336E-09,  1.190E-09,  1.068E-09,  9.623E-10,      &
   &  8.772E-10,  8.007E-10,  7.420E-10,  6.884E-10,  6.483E-10,      &
   &  6.162E-10,  5.922E-10,  5.688E-10,  5.654E-10,  5.637E-10,      &
   &  5.701E-10,  5.781E-10,  5.874E-10,  6.268E-10,  6.357E-10,      &
   &  6.525E-10,  7.137E-10,  7.441E-10,  8.024E-10,  8.485E-10,      &
   &  9.143E-10,  9.536E-10,  9.717E-10,  1.018E-09,  1.042E-09,      &
   &  1.054E-09,  1.092E-09,  1.093E-09,  1.072E-09,  1.046E-09,      &
   &  9.283E-10,  8.053E-10,  6.680E-10,  5.630E-10,  4.933E-10,      &
   &  3.488E-10,  3.324E-10,  2.826E-10,  2.429E-10,  2.051E-10,      &
   &  1.679E-10,  1.411E-10,  1.240E-10,  1.098E-10,  9.668E-11/
   DATA S1451/                                                       &
   &  8.529E-11,  7.848E-11,  7.098E-11,  6.117E-11,  5.246E-11,      &
   &  4.521E-11,  4.059E-11,  3.783E-11,  3.502E-11,  3.271E-11,      &
   &  2.956E-11,  2.636E-11,  2.464E-11,  2.190E-11,  1.986E-11,      &
   &  1.791E-11,  1.638E-11,  1.491E-11,  1.383E-11,  1.271E-11,      &
   &  1.166E-11,  1.082E-11,  1.002E-11,  9.317E-12,  8.697E-12,      &
   &  8.157E-12,  7.678E-12,  7.252E-12,  6.874E-12,  6.539E-12,      &
   &  6.243E-12,  5.985E-12,  5.763E-12,  5.575E-12,  5.421E-12,      &
   &  5.300E-12,  5.212E-12,  5.157E-12,  5.138E-12,  5.156E-12,      &
   &  5.214E-12,  5.313E-12,  5.457E-12,  5.686E-12,  5.953E-12,      &
   &  6.284E-12,  6.856E-12,  7.351E-12,  8.006E-12,  8.957E-12/
   DATA S1501/                                                       &
   &  1.031E-11,  1.170E-11,  1.345E-11,  1.623E-11,  1.999E-11,      &
   &  2.507E-11,  3.004E-11,  3.378E-11,  3.688E-11,  4.118E-11,      &
   &  4.569E-11,  5.025E-11,  5.660E-11,  6.231E-11,  6.881E-11,      &
   &  7.996E-11,  8.526E-11,  9.694E-11,  1.106E-10,  1.222E-10,      &
   &  1.355E-10,  1.525E-10,  1.775E-10,  1.924E-10,  2.181E-10,      &
   &  2.379E-10,  2.662E-10,  2.907E-10,  3.154E-10,  3.366E-10,      &
   &  3.579E-10,  3.858E-10,  4.046E-10,  4.196E-10,  4.166E-10,      &
   &  4.457E-10,  4.466E-10,  4.404E-10,  4.337E-10,  4.150E-10,      &
   &  4.083E-10,  3.910E-10,  3.723E-10,  3.514E-10,  3.303E-10,      &
   &  2.847E-10,  2.546E-10,  2.230E-10,  1.994E-10,  1.733E-10/
   DATA S1551/                                                       &
   &  1.488E-10,  1.297E-10,  1.144E-10,  1.004E-10,  8.741E-11,      &
   &  7.928E-11,  7.034E-11,  6.323E-11,  5.754E-11,  5.250E-11,      &
   &  4.850E-11,  4.502E-11,  4.286E-11,  4.028E-11,  3.899E-11,      &
   &  3.824E-11,  3.761E-11,  3.804E-11,  3.839E-11,  3.845E-11,      &
   &  4.244E-11,  4.382E-11,  4.582E-11,  4.847E-11,  5.209E-11,      &
   &  5.384E-11,  5.887E-11,  6.371E-11,  6.737E-11,  7.168E-11,      &
   &  7.415E-11,  7.827E-11,  8.037E-11,  8.120E-11,  8.154E-11,      &
   &  7.981E-11,  7.783E-11,  6.900E-11,  6.030E-11,  4.900E-11,      &
   &  3.878E-11,  3.068E-11,  2.533E-11,  2.255E-11,  2.056E-11,      &
   &  1.569E-11,  1.398E-11,  1.149E-11,  9.249E-12,  7.706E-12/
   DATA S1601/                                                       &
   &  6.436E-12,  5.514E-12,  4.782E-12,  4.202E-12,  3.721E-12,      &
   &  3.336E-12,  3.015E-12,  2.744E-12,  2.512E-12,  2.312E-12,      &
   &  2.138E-12,  1.987E-12,  1.855E-12,  1.741E-12,  1.642E-12,      &
   &  1.557E-12,  1.483E-12,  1.422E-12,  1.370E-12,  1.329E-12,      &
   &  1.297E-12,  1.274E-12,  1.259E-12,  1.253E-12,  1.255E-12,      &
   &  1.266E-12,  1.284E-12,  1.312E-12,  1.348E-12,  1.393E-12,      &
   &  1.448E-12,  1.512E-12,  1.587E-12,  1.672E-12,  1.770E-12,      &
   &  1.880E-12,  1.990E-12,  2.110E-12,  2.274E-12,  2.516E-12,      &
   &  2.844E-12,  3.231E-12,  3.661E-12,  4.153E-12,  4.717E-12,      &
   &  5.360E-12,  6.094E-12,  6.930E-12,  7.882E-12,  8.966E-12/
   DATA S1651/                                                       &
   &  1.020E-11,  1.162E-11,  1.324E-11,  1.510E-11,  1.720E-11,      &
   &  1.965E-11,  2.237E-11,  2.560E-11,  2.927E-11,  3.371E-11,      &
   &  3.842E-11,  4.429E-11,  5.139E-11,  5.798E-11,  6.697E-11,      &
   &  7.626E-11,  8.647E-11,  1.022E-10,  1.136E-10,  1.300E-10,      &
   &  1.481E-10,  1.672E-10,  1.871E-10,  2.126E-10,  2.357E-10,      &
   &  2.583E-10,  2.997E-10,  3.289E-10,  3.702E-10,  4.012E-10,      &
   &  4.319E-10,  4.527E-10,  5.001E-10,  5.448E-10,  5.611E-10,      &
   &  5.760E-10,  5.965E-10,  6.079E-10,  6.207E-10,  6.276E-10,      &
   &  6.222E-10,  6.137E-10,  6.000E-10,  5.814E-10,  5.393E-10,      &
   &  5.350E-10,  4.947E-10,  4.629E-10,  4.117E-10,  3.712E-10/
   DATA S1701/                                                       &
   &  3.372E-10,  2.923E-10,  2.550E-10,  2.232E-10,  1.929E-10,      &
   &  1.679E-10,  1.460E-10,  1.289E-10,  1.130E-10,  9.953E-11,      &
   &  8.763E-11,  7.760E-11,  6.900E-11,  6.160E-11,  5.525E-11,      &
   &  4.958E-11,  4.489E-11,  4.072E-11,  3.728E-11,  3.438E-11,      &
   &  3.205E-11,  3.006E-11,  2.848E-11,  2.766E-11,  2.688E-11,      &
   &  2.664E-11,  2.670E-11,  2.696E-11,  2.786E-11,  2.861E-11,      &
   &  3.009E-11,  3.178E-11,  3.389E-11,  3.587E-11,  3.819E-11,      &
   &  4.054E-11,  4.417E-11,  4.703E-11,  5.137E-11,  5.460E-11,      &
   &  6.055E-11,  6.333E-11,  6.773E-11,  7.219E-11,  7.717E-11,      &
   &  8.131E-11,  8.491E-11,  8.574E-11,  9.010E-11,  9.017E-11/
   DATA S1751/                                                       &
   &  8.999E-11,  8.959E-11,  8.838E-11,  8.579E-11,  8.162E-11,      &
   &  8.098E-11,  7.472E-11,  7.108E-11,  6.559E-11,  5.994E-11,      &
   &  5.172E-11,  4.424E-11,  3.951E-11,  3.340E-11,  2.902E-11,      &
   &  2.541E-11,  2.215E-11,  1.945E-11,  1.716E-11,  1.503E-11,      &
   &  1.339E-11,  1.185E-11,  1.050E-11,  9.336E-12,  8.307E-12,      &
   &  7.312E-12,  6.550E-12,  5.836E-12,  5.178E-12,  4.600E-12,      &
   &  4.086E-12,  3.639E-12,  3.247E-12,  2.904E-12,  2.604E-12,      &
   &  2.341E-12,  2.112E-12,  1.914E-12,  1.744E-12,  1.598E-12,      &
   &  1.476E-12,  1.374E-12,  1.293E-12,  1.230E-12,  1.185E-12,      &
   &  1.158E-12,  1.147E-12,  1.154E-12,  1.177E-12,  1.219E-12/
   DATA S1801/                                                       &
   &  1.280E-12,  1.360E-12,  1.463E-12,  1.591E-12,  1.750E-12,      &
   &  1.940E-12,  2.156E-12,  2.430E-12,  2.748E-12,  3.052E-12,      &
   &  3.533E-12,  3.967E-12,  4.471E-12,  5.041E-12,  5.860E-12,      &
   &  6.664E-12,  7.522E-12,  8.342E-12,  9.412E-12,  1.072E-11,      &
   &  1.213E-11,  1.343E-11,  1.496E-11,  1.664E-11,  1.822E-11,      &
   &  2.029E-11,  2.233E-11,  2.457E-11,  2.709E-11,  2.928E-11,      &
   &  3.115E-11,  3.356E-11,  3.592E-11,  3.818E-11,  3.936E-11,      &
   &  4.061E-11,  4.149E-11,  4.299E-11,  4.223E-11,  4.251E-11,      &
   &  4.287E-11,  4.177E-11,  4.094E-11,  3.942E-11,  3.772E-11,      &
   &  3.614E-11,  3.394E-11,  3.222E-11,  2.791E-11,  2.665E-11/
   DATA S1851/                                                       &
   &  2.309E-11,  2.032E-11,  1.740E-11,  1.535E-11,  1.323E-11,      &
   &  1.151E-11,  9.803E-12,  8.650E-12,  7.540E-12,  6.619E-12,      &
   &  5.832E-12,  5.113E-12,  4.503E-12,  3.975E-12,  3.520E-12,      &
   &  3.112E-12,  2.797E-12,  2.500E-12,  2.240E-12,  2.013E-12,      &
   &  1.819E-12,  1.653E-12,  1.513E-12,  1.395E-12,  1.299E-12,      &
   &  1.225E-12,  1.168E-12,  1.124E-12,  1.148E-12,  1.107E-12,      &
   &  1.128E-12,  1.169E-12,  1.233E-12,  1.307E-12,  1.359E-12,      &
   &  1.543E-12,  1.686E-12,  1.794E-12,  2.028E-12,  2.210E-12,      &
   &  2.441E-12,  2.653E-12,  2.828E-12,  3.093E-12,  3.280E-12,      &
   &  3.551E-12,  3.677E-12,  3.803E-12,  3.844E-12,  4.068E-12/
   DATA S1901/                                                       &
   &  4.093E-12,  4.002E-12,  3.904E-12,  3.624E-12,  3.633E-12,      &
   &  3.622E-12,  3.443E-12,  3.184E-12,  2.934E-12,  2.476E-12,      &
   &  2.212E-12,  1.867E-12,  1.594E-12,  1.370E-12,  1.192E-12,      &
   &  1.045E-12,  9.211E-13,  8.170E-13,  7.290E-13,  6.550E-13,      &
   &  5.929E-13,  5.415E-13,  4.995E-13,  4.661E-13,  4.406E-13,      &
   &  4.225E-13,  4.116E-13,  4.075E-13,  4.102E-13,  4.198E-13,      &
   &  4.365E-13,  4.606E-13,  4.925E-13,  5.326E-13,  5.818E-13,      &
   &  6.407E-13,  7.104E-13,  7.920E-13,  8.868E-13,  9.964E-13,      &
   &  1.123E-12,  1.268E-12,  1.434E-12,  1.626E-12,  1.848E-12,      &
   &  2.107E-12,  2.422E-12,  2.772E-12,  3.145E-12,  3.704E-12/
   DATA S1951/                                                       &
   &  4.270E-12,  4.721E-12,  5.361E-12,  6.083E-12,  7.095E-12,      &
   &  7.968E-12,  9.228E-12,  1.048E-11,  1.187E-11,  1.336E-11,      &
   &  1.577E-11,  1.772E-11,  2.017E-11,  2.250E-11,  2.630E-11,      &
   &  2.911E-11,  3.356E-11,  3.820E-11,  4.173E-11,  4.811E-11,      &
   &  5.254E-11,  5.839E-11,  6.187E-11,  6.805E-11,  7.118E-11,      &
   &  7.369E-11,  7.664E-11,  7.794E-11,  7.947E-11,  8.036E-11,      &
   &  7.954E-11,  7.849E-11,  7.518E-11,  7.462E-11,  6.926E-11,      &
   &  6.531E-11,  6.197E-11,  5.421E-11,  4.777E-11,  4.111E-11,      &
   &  3.679E-11,  3.166E-11,  2.786E-11,  2.436E-11,  2.144E-11,      &
   &  1.859E-11,  1.628E-11,  1.414E-11,  1.237E-11,  1.093E-11/
   DATA S2001/                                                       &
   &  9.558e-12/
!
end block data BS296
!
!     --------------------------------------------------------------
!
SUBROUTINE SL260 (V1C,V2C,DVC,NPTC,C,v1ss,v2ss)

   Use lblparams, ONLY: n_absrb
!
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   COMMON /S260/ V1S,V2S,DVS,NPTS,S(2003)
   DIMENSION C(*)
!
   DVC = DVS
   v1ss = v1s
   v2ss = v2s
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      C(J) = S(I)
10 END DO
!
   RETURN
!
end subroutine SL260
!
!     --------------------------------------------------------------
!
BLOCK DATA BS260
!
   IMPLICIT REAL*8           (V)
!
!               UNITS OF (CM**3/MOL)*1.E-20
!
   COMMON /S260/ V1,V2,DV,NPT,                                       &
   &              S0000( 2),S0001(50),S0051(50),S0101(50),S0151(50),  &
   &              S0201(50),S0251(50),S0301(50),S0351(50),S0401(50),  &
   &              S0451(50),S0501(50),S0551(50),S0601(50),S0651(50),  &
   &              S0701(50),S0751(50),S0801(50),S0851(50),S0901(50),  &
   &              S0951(50),S1001(50),S1051(50),S1101(50),S1151(50),  &
   &              S1201(50),S1251(50),S1301(50),S1351(50),S1401(50),  &
   &              S1451(50),S1501(50),S1551(50),S1601(50),S1651(50),  &
   &              S1701(50),S1751(50),S1801(50),S1851(50),S1901(50),  &
   &              S1951(50), S2001(1)
!
   DATA V1,V2,DV,NPT / -20.0, 20000.0, 10.0, 2003/
   DATA S0000/                                                       &
   &      5.998e-01, 6.382e-01/
   DATA S0001/                                                       &
   &  6.431E-01,  6.382E-01,  5.998E-01,  5.314E-01,  4.497E-01,      &
   &  3.735E-01,  3.107E-01,  2.615E-01,  2.250E-01,  1.965E-01,      &
   &  1.728E-01,  1.522E-01,  1.338E-01,  1.173E-01,  1.024E-01,      &
   &  8.916E-02,  7.735E-02,  6.688E-02,  5.766E-02,  4.960E-02,      &
   &  4.257E-02,  3.647E-02,  3.120E-02,  2.667E-02,  2.280E-02,      &
   &  1.950E-02,  1.669E-02,  1.430E-02,  1.228E-02,  1.057E-02,      &
   &  9.119E-03,  7.892E-03,  6.894E-03,  6.049E-03,  5.334E-03,      &
   &  4.717E-03,  4.180E-03,  3.728E-03,  3.343E-03,  3.002E-03,      &
   &  2.703E-03,  2.437E-03,  2.201E-03,  1.990E-03,  1.802E-03,      &
   &  1.637E-03,  1.481E-03,  1.345E-03,  1.222E-03,  1.109E-03/     
   DATA S0051/                                                       &
   &  1.007E-03,  9.138E-04,  8.298E-04,  7.536E-04,  6.851E-04,      &
   &  6.215E-04,  5.646E-04,  5.131E-04,  4.662E-04,  4.239E-04,      &
   &  3.859E-04,  3.524E-04,  3.233E-04,  2.963E-04,  2.716E-04,      &
   &  2.493E-04,  2.295E-04,  2.123E-04,  1.980E-04,  1.864E-04,      &
   &  1.769E-04,  1.685E-04,  1.622E-04,  1.558E-04,  1.491E-04,      &
   &  1.426E-04,  1.360E-04,  1.301E-04,  1.243E-04,  1.184E-04,      &
   &  1.120E-04,  1.054E-04,  9.967E-05,  9.444E-05,  8.984E-05,      &
   &  8.537E-05,  8.108E-05,  7.700E-05,  7.306E-05,  6.936E-05,      &
   &  6.565E-05,  6.254E-05,  5.863E-05,  5.513E-05,  5.119E-05,      &
   &  4.833E-05,  4.568E-05,  4.324E-05,  4.104E-05,  3.884E-05/     
   DATA S0101/                                                       &
   &  3.702e-05,  3.511e-05,  3.339e-05,  3.177e-05,  3.026e-05,      &
   &  2.886e-05,  2.756e-05,  2.636e-05,  2.527e-05,  2.427e-05,      &
   &  2.337e-05,  2.257e-05,  2.185e-05,  2.127e-05,  2.080e-05,      &
   &  2.041e-05,  2.013e-05,  2.000e-05,  1.997e-05,  2.009e-05,      &
   &  2.031e-05,  2.068e-05,  2.124e-05,  2.189e-05,  2.267e-05,      &
   &  2.364e-05,  2.463e-05,  2.618e-05,  2.774e-05,  2.937e-05,      &
   &  3.144e-05,  3.359e-05,  3.695e-05,  4.002e-05,  4.374e-05,      &
   &  4.947e-05,  5.431e-05,  6.281e-05,  7.169e-05,  8.157e-05,      &
   &  9.728e-05,  1.079e-04,  1.337e-04,  1.442e-04,  1.683e-04,      &
   &  1.879e-04,  2.223e-04,  2.425e-04,  2.838e-04,  3.143e-04/
   DATA S0151/                                                       &
   &  3.527e-04,  4.012e-04,  4.237e-04,  4.747e-04,  5.057e-04,      &
   &  5.409e-04,  5.734e-04,  5.944e-04,  6.077e-04,  6.175e-04,      &
   &  6.238e-04,  6.226e-04,  6.248e-04,  6.192e-04,  6.098e-04,      &
   &  5.818e-04,  5.709e-04,  5.465e-04,  5.043e-04,  4.699e-04,      &
   &  4.294e-04,  3.984e-04,  3.672e-04,  3.152e-04,  2.883e-04,      &
   &  2.503e-04,  2.211e-04,  1.920e-04,  1.714e-04,  1.485e-04,      &
   &  1.358e-04,  1.156e-04,  1.021e-04,  8.887e-05,  7.842e-05,      &
   &  7.120e-05,  6.186e-05,  5.730e-05,  4.792e-05,  4.364e-05,      &
   &  3.720e-05,  3.280e-05,  2.946e-05,  2.591e-05,  2.261e-05,      &
   &  2.048e-05,  1.813e-05,  1.630e-05,  1.447e-05,  1.282e-05/
   DATA S0201/                                                       &
   &  1.156E-05,  1.072E-05,  9.922E-06,  8.919E-06,  7.974E-06,      &
   &  7.050E-06,  6.301E-06,  5.559E-06,  4.978E-06,  4.371E-06,      &
   &  3.941E-06,  3.479E-06,  3.119E-06,  2.795E-06,  2.624E-06,      &
   &  2.458E-06,  2.321E-06,  2.192E-06,  2.077E-06,  1.966E-06,      &
   &  1.863E-06,  1.775E-06,  1.691E-06,  1.617E-06,  1.545E-06,      &
   &  1.484E-06,  1.422E-06,  1.371E-06,  1.320E-06,  1.275E-06,      &
   &  1.233E-06,  1.187E-06,  1.145E-06,  1.118E-06,  1.082E-06,      &
   &  1.048E-06,  1.020E-06,  9.889E-07,  9.649E-07,  9.368E-07,      &
   &  9.174E-07,  8.960E-07,  8.763E-07,  8.566E-07,  8.381E-07,      &
   &  8.192E-07,  8.033E-07,  7.864E-07,  7.715E-07,  7.563E-07/
   DATA S0251/                                                       &
   &  7.429E-07,  7.282E-07,  7.149E-07,  7.011E-07,  6.901E-07,      &
   &  6.771E-07,  6.671E-07,  6.577E-07,  6.476E-07,  6.376E-07,      &
   &  6.313E-07,  6.233E-07,  6.171E-07,  6.163E-07,  6.162E-07,      &
   &  6.186E-07,  6.241E-07,  6.337E-07,  6.420E-07,  6.421E-07,      &
   &  6.368E-07,  6.269E-07,  6.141E-07,  6.023E-07,  5.886E-07,      &
   &  5.730E-07,  5.601E-07,  5.496E-07,  5.414E-07,  5.354E-07,      &
   &  5.285E-07,  5.233E-07,  5.199E-07,  5.168E-07,  5.145E-07,      &
   &  5.125E-07,  5.112E-07,  5.125E-07,  5.135E-07,  5.209E-07,      &
   &  5.256E-07,  5.357E-07,  5.503E-07,  5.636E-07,  5.891E-07,      &
   &  6.194E-07,  6.511E-07,  6.915E-07,  7.401E-07,  8.004E-07/
   DATA S0301/                                                       &
   &  8.606E-07,  9.087E-07,  1.003E-06,  1.043E-06,  1.190E-06,      &
   &  1.293E-06,  1.398E-06,  1.508E-06,  1.649E-06,  1.735E-06,      &
   &  1.852E-06,  1.955E-06,  2.056E-06,  2.079E-06,  2.127E-06,      &
   &  2.117E-06,  2.109E-06,  2.101E-06,  2.122E-06,  2.097E-06,      &
   &  2.142E-06,  2.139E-06,  2.138E-06,  2.097E-06,  2.069E-06,      &
   &  2.072E-06,  2.038E-06,  2.020E-06,  2.052E-06,  2.046E-06,      &
   &  2.077E-06,  2.103E-06,  2.163E-06,  2.244E-06,  2.388E-06,      &
   &  2.552E-06,  2.773E-06,  3.067E-06,  3.360E-06,  3.714E-06,      &
   &  4.172E-06,  4.643E-06,  5.116E-06,  5.831E-06,  6.522E-06,      &
   &  7.478E-06,  8.312E-06,  9.407E-06,  1.051E-05,  1.207E-05/
   DATA S0351/                                                       &
   &  1.359E-05,  1.568E-05,  1.759E-05,  2.029E-05,  2.284E-05,      &
   &  2.602E-05,  2.940E-05,  3.483E-05,  3.928E-05,  4.618E-05,      &
   &  5.240E-05,  6.132E-05,  7.183E-05,  8.521E-05,  9.111E-05,      &
   &  1.070E-04,  1.184E-04,  1.264E-04,  1.475E-04,  1.612E-04,      &
   &  1.704E-04,  1.818E-04,  1.924E-04,  1.994E-04,  2.061E-04,      &
   &  2.180E-04,  2.187E-04,  2.200E-04,  2.196E-04,  2.131E-04,      &
   &  2.015E-04,  1.988E-04,  1.847E-04,  1.729E-04,  1.597E-04,      &
   &  1.373E-04,  1.262E-04,  1.087E-04,  9.439E-05,  8.061E-05,      &
   &  7.093E-05,  6.049E-05,  5.120E-05,  4.435E-05,  3.817E-05,      &
   &  3.340E-05,  2.927E-05,  2.573E-05,  2.291E-05,  2.040E-05/
   DATA S0401/                                                       &
   &  1.827E-05,  1.636E-05,  1.463E-05,  1.309E-05,  1.170E-05,      &
   &  1.047E-05,  9.315E-06,  8.328E-06,  7.458E-06,  6.665E-06,      &
   &  5.940E-06,  5.316E-06,  4.752E-06,  4.252E-06,  3.825E-06,      &
   &  3.421E-06,  3.064E-06,  2.746E-06,  2.465E-06,  2.216E-06,      &
   &  1.990E-06,  1.790E-06,  1.609E-06,  1.449E-06,  1.306E-06,      &
   &  1.171E-06,  1.056E-06,  9.738E-07,  9.212E-07,  8.894E-07,      &
   &  8.561E-07,  8.183E-07,  7.824E-07,  7.484E-07,  7.168E-07,      &
   &  6.872E-07,  6.590E-07,  6.325E-07,  6.070E-07,  5.833E-07,      &
   &  5.610E-07,  5.400E-07,  5.197E-07,  5.005E-07,  4.823E-07,      &
   &  4.652E-07,  4.483E-07,  4.329E-07,  4.177E-07,  4.029E-07/
   DATA S0451/                                                        &
   &  3.888E-07,  3.749E-07,  3.618E-07,  3.492E-07,  3.377E-07,      &
   &  3.269E-07,  3.166E-07,  3.081E-07,  3.010E-07,  2.947E-07,      &
   &  2.890E-07,  2.844E-07,  2.789E-07,  2.731E-07,  2.670E-07,      &
   &  2.604E-07,  2.525E-07,  2.444E-07,  2.384E-07,  2.328E-07,      &
   &  2.280E-07,  2.230E-07,  2.183E-07,  2.144E-07,  2.105E-07,      &
   &  2.073E-07,  2.046E-07,  2.022E-07,  2.002E-07,  1.985E-07,      &
   &  1.975E-07,  1.965E-07,  1.960E-07,  1.954E-07,  1.954E-07,      &
   &  1.960E-07,  1.971E-07,  1.984E-07,  2.003E-07,  2.026E-07,      &
   &  2.058E-07,  2.099E-07,  2.142E-07,  2.171E-07,  2.206E-07,      &
   &  2.312E-07,  2.523E-07,  2.834E-07,  3.195E-07,  3.587E-07/
   DATA S0501/                                                       &
   &  4.015E-07,  4.497E-07,  5.049E-07,  5.665E-07,  6.366E-07,      &
   &  7.121E-07,  7.996E-07,  8.946E-07,  1.002E-06,  1.117E-06,      &
   &  1.262E-06,  1.416E-06,  1.611E-06,  1.807E-06,  2.056E-06,      &
   &  2.351E-06,  2.769E-06,  3.138E-06,  3.699E-06,  4.386E-06,      &
   &  5.041E-06,  6.074E-06,  6.812E-06,  7.790E-06,  8.855E-06,      &
   &  1.014E-05,  1.095E-05,  1.245E-05,  1.316E-05,  1.390E-05,      &
   &  1.504E-05,  1.583E-05,  1.617E-05,  1.652E-05,  1.713E-05,      &
   &  1.724E-05,  1.715E-05,  1.668E-05,  1.629E-05,  1.552E-05,      &
   &  1.478E-05,  1.340E-05,  1.245E-05,  1.121E-05,  9.575E-06,      &
   &  8.956E-06,  7.345E-06,  6.597E-06,  5.612E-06,  4.818E-06/
   DATA S0551/                                                       &
   &  4.165E-06,  3.579E-06,  3.041E-06,  2.623E-06,  2.290E-06,      &
   &  1.984E-06,  1.748E-06,  1.534E-06,  1.369E-06,  1.219E-06,      &
   &  1.092E-06,  9.800E-07,  8.762E-07,  7.896E-07,  7.104E-07,      &
   &  6.364E-07,  5.691E-07,  5.107E-07,  4.575E-07,  4.090E-07,      &
   &  3.667E-07,  3.287E-07,  2.931E-07,  2.633E-07,  2.356E-07,      &
   &  2.111E-07,  1.895E-07,  1.697E-07,  1.525E-07,  1.369E-07,      &
   &  1.233E-07,  1.114E-07,  9.988E-08,  9.004E-08,  8.149E-08,      &
   &  7.311E-08,  6.587E-08,  6.014E-08,  5.618E-08,  5.320E-08,      &
   &  5.052E-08,  4.768E-08,  4.514E-08,  4.267E-08,  4.032E-08,      &
   &  3.815E-08,  3.605E-08,  3.409E-08,  3.219E-08,  3.038E-08/
   DATA S0601/                                                       &
   &  2.869E-08,  2.710E-08,  2.563E-08,  2.418E-08,  2.288E-08,      &
   &  2.166E-08,  2.051E-08,  1.945E-08,  1.844E-08,  1.745E-08,      &
   &  1.650E-08,  1.558E-08,  1.471E-08,  1.387E-08,  1.308E-08,      &
   &  1.235E-08,  1.169E-08,  1.107E-08,  1.050E-08,  9.963E-09,      &
   &  9.461E-09,  9.004E-09,  8.580E-09,  8.191E-09,  7.835E-09,      &
   &  7.525E-09,  7.242E-09,  7.004E-09,  6.798E-09,  6.617E-09,      &
   &  6.463E-09,  6.355E-09,  6.267E-09,  6.194E-09,  6.156E-09,      &
   &  6.144E-09,  6.159E-09,  6.210E-09,  6.288E-09,  6.398E-09,      &
   &  6.536E-09,  6.719E-09,  6.929E-09,  7.176E-09,  7.459E-09,      &
   &  7.787E-09,  8.147E-09,  8.573E-09,  9.024E-09,  9.547E-09/
   DATA S0651/                                                       &
   &  1.011E-08,  1.073E-08,  1.146E-08,  1.217E-08,  1.306E-08,      &
   &  1.411E-08,  1.525E-08,  1.637E-08,  1.776E-08,  1.958E-08,      &
   &  2.136E-08,  2.372E-08,  2.612E-08,  2.944E-08,  3.249E-08,      &
   &  3.740E-08,  4.244E-08,  4.845E-08,  5.630E-08,  6.728E-08,      &
   &  7.908E-08,  9.511E-08,  1.128E-07,  1.358E-07,  1.587E-07,      &
   &  1.904E-07,  2.224E-07,  2.660E-07,  3.278E-07,  3.775E-07,      &
   &  4.355E-07,  4.672E-07,  5.110E-07,  5.461E-07,  5.828E-07,      &
   &  6.233E-07,  6.509E-07,  6.672E-07,  6.969E-07,  7.104E-07,      &
   &  7.439E-07,  7.463E-07,  7.708E-07,  7.466E-07,  7.668E-07,      &
   &  7.549E-07,  7.586E-07,  7.384E-07,  7.439E-07,  7.785E-07/
   DATA S0701/                                                       &
   &  7.915E-07,  8.310E-07,  8.745E-07,  9.558E-07,  1.038E-06,      &
   &  1.173E-06,  1.304E-06,  1.452E-06,  1.671E-06,  1.931E-06,      &
   &  2.239E-06,  2.578E-06,  3.032E-06,  3.334E-06,  3.980E-06,      &
   &  4.300E-06,  4.518E-06,  5.321E-06,  5.508E-06,  6.211E-06,      &
   &  6.590E-06,  7.046E-06,  7.555E-06,  7.558E-06,  7.875E-06,      &
   &  8.319E-06,  8.433E-06,  8.590E-06,  8.499E-06,  8.284E-06,      &
   &  8.276E-06,  7.587E-06,  7.014E-06,  6.397E-06,  5.520E-06,      &
   &  4.888E-06,  4.024E-06,  3.433E-06,  2.726E-06,  2.169E-06,      &
   &  1.811E-06,  1.438E-06,  1.156E-06,  9.272E-07,  7.455E-07,      &
   &  6.056E-07,  4.980E-07,  4.116E-07,  3.425E-07,  2.878E-07/
   DATA S0751/                                                       &
   &  2.455E-07,  2.114E-07,  1.862E-07,  1.667E-07,  1.507E-07,      &
   &  1.370E-07,  1.260E-07,  1.153E-07,  1.065E-07,  9.812E-08,      &
   &  9.045E-08,  8.334E-08,  7.719E-08,  7.133E-08,  6.653E-08,      &
   &  6.160E-08,  5.753E-08,  5.349E-08,  5.005E-08,  4.693E-08,      &
   &  4.388E-08,  4.127E-08,  3.883E-08,  3.653E-08,  3.443E-08,      &
   &  3.256E-08,  3.075E-08,  2.899E-08,  2.735E-08,  2.586E-08,      &
   &  2.443E-08,  2.305E-08,  2.177E-08,  2.049E-08,  1.938E-08,      &
   &  1.833E-08,  1.734E-08,  1.638E-08,  1.548E-08,  1.464E-08,      &
   &  1.386E-08,  1.310E-08,  1.239E-08,  1.173E-08,  1.110E-08,      &
   &  1.051E-08,  9.947E-09,  9.414E-09,  8.911E-09,  8.435E-09/
   DATA S0801/                                                       &
   &  7.989E-09,  7.567E-09,  7.174E-09,  6.804E-09,  6.460E-09,      &
   &  6.132E-09,  5.828E-09,  5.550E-09,  5.286E-09,  5.049E-09,      &
   &  4.831E-09,  4.658E-09,  4.505E-09,  4.413E-09,  4.362E-09,      &
   &  4.489E-09,  4.487E-09,  4.761E-09,  5.278E-09,  5.630E-09,      &
   &  6.056E-09,  6.633E-09,  7.337E-09,  7.954E-09,  8.869E-09,      &
   &  9.500E-09,  1.027E-08,  1.087E-08,  1.204E-08,  1.304E-08,      &
   &  1.413E-08,  1.529E-08,  1.586E-08,  1.649E-08,  1.661E-08,      &
   &  1.703E-08,  1.741E-08,  1.766E-08,  1.779E-08,  1.816E-08,      &
   &  1.866E-08,  1.889E-08,  1.904E-08,  1.897E-08,  1.893E-08,      &
   &  1.888E-08,  1.868E-08,  1.895E-08,  1.899E-08,  1.876E-08/
   DATA S0851/                                                       &
   &  1.960E-08,  2.020E-08,  2.121E-08,  2.239E-08,  2.379E-08,      &
   &  2.526E-08,  2.766E-08,  2.994E-08,  3.332E-08,  3.703E-08,      &
   &  4.158E-08,  4.774E-08,  5.499E-08,  6.355E-08,  7.349E-08,      &
   &  8.414E-08,  9.846E-08,  1.143E-07,  1.307E-07,  1.562E-07,      &
   &  1.817E-07,  2.011E-07,  2.192E-07,  2.485E-07,  2.867E-07,      &
   &  3.035E-07,  3.223E-07,  3.443E-07,  3.617E-07,  3.793E-07,      &
   &  3.793E-07,  3.839E-07,  4.081E-07,  4.117E-07,  4.085E-07,      &
   &  3.920E-07,  3.851E-07,  3.754E-07,  3.490E-07,  3.229E-07,      &
   &  2.978E-07,  2.691E-07,  2.312E-07,  2.029E-07,  1.721E-07,      &
   &  1.472E-07,  1.308E-07,  1.132E-07,  9.736E-08,  8.458E-08/
   DATA S0901/                                                       &
   &  7.402E-08,  6.534E-08,  5.811E-08,  5.235E-08,  4.762E-08,      &
   &  4.293E-08,  3.896E-08,  3.526E-08,  3.165E-08,  2.833E-08,      &
   &  2.551E-08,  2.288E-08,  2.036E-08,  1.820E-08,  1.626E-08,      &
   &  1.438E-08,  1.299E-08,  1.149E-08,  1.030E-08,  9.148E-09,      &
   &  8.122E-09,  7.264E-09,  6.425E-09,  5.777E-09,  5.060E-09,      &
   &  4.502E-09,  4.013E-09,  3.567E-09,  3.145E-09,  2.864E-09,      &
   &  2.553E-09,  2.311E-09,  2.087E-09,  1.886E-09,  1.716E-09,      &
   &  1.556E-09,  1.432E-09,  1.311E-09,  1.196E-09,  1.091E-09,      &
   &  1.006E-09,  9.428E-10,  8.900E-10,  8.438E-10,  8.005E-10,      &
   &  7.597E-10,  7.201E-10,  6.818E-10,  6.550E-10,  6.286E-10/
   DATA S0951/                                                       &
   &  6.004E-10,  5.730E-10,  5.466E-10,  5.203E-10,  4.968E-10,      &
   &  4.741E-10,  4.520E-10,  4.314E-10,  4.124E-10,  3.944E-10,      &
   &  3.777E-10,  3.619E-10,  3.477E-10,  3.344E-10,  3.217E-10,      &
   &  3.096E-10,  2.999E-10,  2.898E-10,  2.826E-10,  2.750E-10,      &
   &  2.701E-10,  2.651E-10,  2.641E-10,  2.671E-10,  2.670E-10,      &
   &  2.748E-10,  2.863E-10,  3.051E-10,  3.163E-10,  3.284E-10,      &
   &  3.408E-10,  3.601E-10,  3.710E-10,  3.878E-10,  4.028E-10,      &
   &  4.159E-10,  4.296E-10,  4.371E-10,  4.463E-10,  4.520E-10,      &
   &  4.554E-10,  4.716E-10,  4.791E-10,  4.900E-10,  5.001E-10,      &
   &  4.999E-10,  5.061E-10,  5.291E-10,  5.416E-10,  5.994E-10/
   DATA S1001/                                                       &
   &  6.394E-10,  6.807E-10,  7.250E-10,  7.847E-10,  8.430E-10,      &
   &  9.146E-10,  9.907E-10,  1.062E-09,  1.170E-09,  1.277E-09,      &
   &  1.401E-09,  1.524E-09,  1.674E-09,  1.888E-09,  2.112E-09,      &
   &  2.322E-09,  2.808E-09,  3.269E-09,  3.902E-09,  4.614E-09,      &
   &  5.763E-09,  6.747E-09,  8.972E-09,  1.118E-08,  1.428E-08,      &
   &  1.713E-08,  1.948E-08,  2.090E-08,  2.211E-08,  2.362E-08,      &
   &  2.556E-08,  2.729E-08,  2.880E-08,  3.046E-08,  3.167E-08,      &
   &  3.367E-08,  3.457E-08,  3.590E-08,  3.711E-08,  3.826E-08,      &
   &  4.001E-08,  4.211E-08,  4.315E-08,  4.661E-08,  5.010E-08,      &
   &  5.249E-08,  5.840E-08,  6.628E-08,  7.512E-08,  8.253E-08/
   DATA S1051/                                                       &
   &  9.722E-08,  1.067E-07,  1.153E-07,  1.347E-07,  1.428E-07,      &
   &  1.577E-07,  1.694E-07,  1.833E-07,  1.938E-07,  2.108E-07,      &
   &  2.059E-07,  2.157E-07,  2.185E-07,  2.208E-07,  2.182E-07,      &
   &  2.093E-07,  2.014E-07,  1.962E-07,  1.819E-07,  1.713E-07,      &
   &  1.510E-07,  1.340E-07,  1.154E-07,  9.890E-08,  8.880E-08,      &
   &  7.673E-08,  6.599E-08,  5.730E-08,  5.081E-08,  4.567E-08,      &
   &  4.147E-08,  3.773E-08,  3.460E-08,  3.194E-08,  2.953E-08,      &
   &  2.759E-08,  2.594E-08,  2.442E-08,  2.355E-08,  2.283E-08,      &
   &  2.279E-08,  2.231E-08,  2.279E-08,  2.239E-08,  2.210E-08,      &
   &  2.309E-08,  2.293E-08,  2.352E-08,  2.415E-08,  2.430E-08/
   DATA S1101/                                                       &
   &  2.426E-08,  2.465E-08,  2.500E-08,  2.496E-08,  2.465E-08,      &
   &  2.475E-08,  2.179E-08,  1.934E-08,  1.676E-08,  1.415E-08,      &
   &  1.115E-08,  9.148E-09,  7.692E-09,  6.193E-09,  5.401E-09,      &
   &  4.341E-09,  3.734E-09,  3.127E-09,  2.619E-09,  2.181E-09,      &
   &  1.887E-09,  1.625E-09,  1.456E-09,  1.304E-09,  1.190E-09,      &
   &  1.083E-09,  1.012E-09,  9.153E-10,  8.434E-10,  7.828E-10,      &
   &  7.210E-10,  6.665E-10,  6.136E-10,  5.728E-10,  5.351E-10,      &
   &  4.985E-10,  4.633E-10,  4.337E-10,  4.069E-10,  3.825E-10,      &
   &  3.596E-10,  3.374E-10,  3.186E-10,  3.006E-10,  2.837E-10,      &
   &  2.681E-10,  2.538E-10,  2.405E-10,  2.277E-10,  2.167E-10/
   DATA S1151/                                                       &
   &  2.060E-10,  1.961E-10,  1.871E-10,  1.788E-10,  1.711E-10,      &
   &  1.639E-10,  1.574E-10,  1.515E-10,  1.456E-10,  1.411E-10,      &
   &  1.365E-10,  1.315E-10,  1.296E-10,  1.262E-10,  1.264E-10,      &
   &  1.257E-10,  1.252E-10,  1.327E-10,  1.392E-10,  1.519E-10,      &
   &  1.664E-10,  1.911E-10,  2.276E-10,  2.702E-10,  3.086E-10,      &
   &  3.378E-10,  3.632E-10,  3.957E-10,  4.360E-10,  4.701E-10,      &
   &  5.030E-10,  5.381E-10,  5.793E-10,  6.190E-10,  6.596E-10,      &
   &  7.004E-10,  7.561E-10,  7.934E-10,  8.552E-10,  9.142E-10,      &
   &  9.570E-10,  1.027E-09,  1.097E-09,  1.193E-09,  1.334E-09,      &
   &  1.470E-09,  1.636E-09,  1.871E-09,  2.122E-09,  2.519E-09/
   DATA S1201/                                                       &
   &  2.806E-09,  3.203E-09,  3.846E-09,  4.362E-09,  5.114E-09,      &
   &  5.643E-09,  6.305E-09,  6.981E-09,  7.983E-09,  8.783E-09,      &
   &  9.419E-09,  1.017E-08,  1.063E-08,  1.121E-08,  1.130E-08,      &
   &  1.201E-08,  1.225E-08,  1.232E-08,  1.223E-08,  1.177E-08,      &
   &  1.151E-08,  1.116E-08,  1.047E-08,  9.698E-09,  8.734E-09,      &
   &  8.202E-09,  7.041E-09,  6.074E-09,  5.172E-09,  4.468E-09,      &
   &  3.913E-09,  3.414E-09,  2.975E-09,  2.650E-09,  2.406E-09,      &
   &  2.173E-09,  2.009E-09,  1.861E-09,  1.727E-09,  1.612E-09,      &
   &  1.514E-09,  1.430E-09,  1.362E-09,  1.333E-09,  1.288E-09,      &
   &  1.249E-09,  1.238E-09,  1.228E-09,  1.217E-09,  1.202E-09/
   DATA S1251/                                                       &
   &  1.209E-09,  1.177E-09,  1.157E-09,  1.165E-09,  1.142E-09,      &
   &  1.131E-09,  1.138E-09,  1.127E-09,  1.131E-09,  1.039E-09,      &
   &  9.049E-10,  8.118E-10,  6.369E-10,  5.329E-10,  4.016E-10,      &
   &  3.329E-10,  2.698E-10,  2.539E-10,  2.236E-10,  1.814E-10,      &
   &  1.560E-10,  1.385E-10,  1.135E-10,  9.918E-11,  8.391E-11,      &
   &  7.496E-11,  6.715E-11,  6.105E-11,  5.585E-11,  5.160E-11,      &
   &  4.764E-11,  4.401E-11,  4.099E-11,  3.816E-11,  3.565E-11,      &
   &  3.342E-11,  3.142E-11,  2.954E-11,  2.776E-11,  2.624E-11,      &
   &  2.485E-11,  2.355E-11,  2.237E-11,  2.128E-11,  2.031E-11,      &
   &  1.942E-11,  1.859E-11,  1.785E-11,  1.718E-11,  1.658E-11/
   DATA S1301/                                                       &
   &  1.603E-11,  1.555E-11,  1.514E-11,  1.477E-11,  1.447E-11,      &
   &  1.423E-11,  1.405E-11,  1.393E-11,  1.388E-11,  1.388E-11,      &
   &  1.396E-11,  1.409E-11,  1.430E-11,  1.459E-11,  1.494E-11,      &
   &  1.538E-11,  1.590E-11,  1.650E-11,  1.719E-11,  1.800E-11,      &
   &  1.890E-11,  1.996E-11,  2.118E-11,  2.247E-11,  2.391E-11,      &
   &  2.563E-11,  2.695E-11,  2.866E-11,  3.096E-11,  3.392E-11,      &
   &  3.806E-11,  4.261E-11,  4.748E-11,  5.323E-11,  5.935E-11,      &
   &  6.619E-11,  7.418E-11,  8.294E-11,  9.260E-11,  1.039E-10,      &
   &  1.156E-10,  1.297E-10,  1.460E-10,  1.641E-10,  1.858E-10,      &
   &  2.100E-10,  2.383E-10,  2.724E-10,  3.116E-10,  3.538E-10/
   DATA S1351/                                                       &
   &  4.173E-10,  4.727E-10,  5.503E-10,  6.337E-10,  7.320E-10,      &
   &  8.298E-10,  9.328E-10,  1.059E-09,  1.176E-09,  1.328E-09,      &
   &  1.445E-09,  1.593E-09,  1.770E-09,  1.954E-09,  2.175E-09,      &
   &  2.405E-09,  2.622E-09,  2.906E-09,  3.294E-09,  3.713E-09,      &
   &  3.980E-09,  4.384E-09,  4.987E-09,  5.311E-09,  5.874E-09,      &
   &  6.337E-09,  7.027E-09,  7.390E-09,  7.769E-09,  8.374E-09,      &
   &  8.605E-09,  9.165E-09,  9.415E-09,  9.511E-09,  9.704E-09,      &
   &  9.588E-09,  9.450E-09,  9.086E-09,  8.798E-09,  8.469E-09,      &
   &  7.697E-09,  7.168E-09,  6.255E-09,  5.772E-09,  4.970E-09,      &
   &  4.271E-09,  3.653E-09,  3.154E-09,  2.742E-09,  2.435E-09/
   DATA S1401/                                                       &
   &  2.166E-09,  1.936E-09,  1.731E-09,  1.556E-09,  1.399E-09,      &
   &  1.272E-09,  1.157E-09,  1.066E-09,  9.844E-10,  9.258E-10,      &
   &  8.787E-10,  8.421E-10,  8.083E-10,  8.046E-10,  8.067E-10,      &
   &  8.181E-10,  8.325E-10,  8.517E-10,  9.151E-10,  9.351E-10,      &
   &  9.677E-10,  1.071E-09,  1.126E-09,  1.219E-09,  1.297E-09,      &
   &  1.408E-09,  1.476E-09,  1.517E-09,  1.600E-09,  1.649E-09,      &
   &  1.678E-09,  1.746E-09,  1.764E-09,  1.741E-09,  1.703E-09,      &
   &  1.506E-09,  1.298E-09,  1.066E-09,  8.872E-10,  7.676E-10,      &
   &  5.313E-10,  4.978E-10,  4.163E-10,  3.527E-10,  2.929E-10,      &
   &  2.372E-10,  1.992E-10,  1.767E-10,  1.588E-10,  1.419E-10/
   DATA S1451/                                                       &
   &  1.275E-10,  1.186E-10,  1.082E-10,  9.390E-11,  8.110E-11,      &
   &  7.039E-11,  6.366E-11,  5.968E-11,  5.547E-11,  5.202E-11,      &
   &  4.730E-11,  4.240E-11,  3.969E-11,  3.535E-11,  3.219E-11,      &
   &  2.916E-11,  2.677E-11,  2.445E-11,  2.279E-11,  2.107E-11,      &
   &  1.950E-11,  1.830E-11,  1.716E-11,  1.617E-11,  1.536E-11,      &
   &  1.463E-11,  1.394E-11,  1.331E-11,  1.273E-11,  1.220E-11,      &
   &  1.171E-11,  1.129E-11,  1.091E-11,  1.057E-11,  1.028E-11,      &
   &  1.002E-11,  9.801E-12,  9.618E-12,  9.471E-12,  9.371E-12,      &
   &  9.313E-12,  9.311E-12,  9.387E-12,  9.547E-12,  9.799E-12,      &
   &  1.013E-11,  1.082E-11,  1.141E-11,  1.220E-11,  1.352E-11/
   DATA S1501/                                                       &
   &  1.538E-11,  1.732E-11,  1.978E-11,  2.366E-11,  2.896E-11,      &
   &  3.601E-11,  4.259E-11,  4.737E-11,  5.089E-11,  5.592E-11,      &
   &  6.109E-11,  6.628E-11,  7.381E-11,  8.088E-11,  8.966E-11,      &
   &  1.045E-10,  1.120E-10,  1.287E-10,  1.486E-10,  1.662E-10,      &
   &  1.866E-10,  2.133E-10,  2.524E-10,  2.776E-10,  3.204E-10,      &
   &  3.559E-10,  4.028E-10,  4.448E-10,  4.882E-10,  5.244E-10,      &
   &  5.605E-10,  6.018E-10,  6.328E-10,  6.579E-10,  6.541E-10,      &
   &  7.024E-10,  7.074E-10,  7.068E-10,  7.009E-10,  6.698E-10,      &
   &  6.545E-10,  6.209E-10,  5.834E-10,  5.412E-10,  5.001E-10,      &
   &  4.231E-10,  3.727E-10,  3.211E-10,  2.833E-10,  2.447E-10/
   DATA S1551/                                                       &
   &  2.097E-10,  1.843E-10,  1.639E-10,  1.449E-10,  1.270E-10,      &
   &  1.161E-10,  1.033E-10,  9.282E-11,  8.407E-11,  7.639E-11,      &
   &  7.023E-11,  6.474E-11,  6.142E-11,  5.760E-11,  5.568E-11,      &
   &  5.472E-11,  5.390E-11,  5.455E-11,  5.540E-11,  5.587E-11,      &
   &  6.230E-11,  6.490E-11,  6.868E-11,  7.382E-11,  8.022E-11,      &
   &  8.372E-11,  9.243E-11,  1.004E-10,  1.062E-10,  1.130E-10,      &
   &  1.176E-10,  1.244E-10,  1.279E-10,  1.298E-10,  1.315E-10,      &
   &  1.308E-10,  1.284E-10,  1.138E-10,  9.898E-11,  7.972E-11,      &
   &  6.261E-11,  4.860E-11,  3.916E-11,  3.416E-11,  3.043E-11,      &
   &  2.277E-11,  1.999E-11,  1.632E-11,  1.326E-11,  1.127E-11/
   DATA S1601/                                                       &
   &  9.680E-12,  8.474E-12,  7.490E-12,  6.688E-12,  6.012E-12,      &
   &  5.468E-12,  5.010E-12,  4.623E-12,  4.285E-12,  3.994E-12,      &
   &  3.741E-12,  3.519E-12,  3.326E-12,  3.157E-12,  3.011E-12,      &
   &  2.885E-12,  2.773E-12,  2.683E-12,  2.604E-12,  2.542E-12,      &
   &  2.493E-12,  2.458E-12,  2.434E-12,  2.423E-12,  2.424E-12,      &
   &  2.440E-12,  2.465E-12,  2.506E-12,  2.560E-12,  2.626E-12,      &
   &  2.707E-12,  2.803E-12,  2.914E-12,  3.039E-12,  3.184E-12,      &
   &  3.345E-12,  3.501E-12,  3.670E-12,  3.910E-12,  4.273E-12,      &
   &  4.770E-12,  5.347E-12,  5.978E-12,  6.682E-12,  7.467E-12,      &
   &  8.340E-12,  9.293E-12,  1.035E-11,  1.152E-11,  1.285E-11/
   DATA S1651/                                                       &
   &  1.428E-11,  1.586E-11,  1.764E-11,  1.972E-11,  2.214E-11,      &
   &  2.478E-11,  2.776E-11,  3.151E-11,  3.591E-11,  4.103E-11,      &
   &  4.660E-11,  5.395E-11,  6.306E-11,  7.172E-11,  8.358E-11,      &
   &  9.670E-11,  1.110E-10,  1.325E-10,  1.494E-10,  1.736E-10,      &
   &  2.007E-10,  2.296E-10,  2.608E-10,  3.004E-10,  3.361E-10,      &
   &  3.727E-10,  4.373E-10,  4.838E-10,  5.483E-10,  6.006E-10,      &
   &  6.535E-10,  6.899E-10,  7.687E-10,  8.444E-10,  8.798E-10,      &
   &  9.135E-10,  9.532E-10,  9.757E-10,  9.968E-10,  1.006E-09,      &
   &  9.949E-10,  9.789E-10,  9.564E-10,  9.215E-10,  8.510E-10,      &
   &  8.394E-10,  7.707E-10,  7.152E-10,  6.274E-10,  5.598E-10/
   DATA S1701/                                                       &
   &  5.028E-10,  4.300E-10,  3.710E-10,  3.245E-10,  2.809E-10,      &
   &  2.461E-10,  2.154E-10,  1.910E-10,  1.685E-10,  1.487E-10,      &
   &  1.313E-10,  1.163E-10,  1.031E-10,  9.172E-11,  8.221E-11,      &
   &  7.382E-11,  6.693E-11,  6.079E-11,  5.581E-11,  5.167E-11,      &
   &  4.811E-11,  4.506E-11,  4.255E-11,  4.083E-11,  3.949E-11,      &
   &  3.881E-11,  3.861E-11,  3.858E-11,  3.951E-11,  4.045E-11,      &
   &  4.240E-11,  4.487E-11,  4.806E-11,  5.133E-11,  5.518E-11,      &
   &  5.919E-11,  6.533E-11,  7.031E-11,  7.762E-11,  8.305E-11,      &
   &  9.252E-11,  9.727E-11,  1.045E-10,  1.117E-10,  1.200E-10,      &
   &  1.275E-10,  1.341E-10,  1.362E-10,  1.438E-10,  1.450E-10/
   DATA S1751/                                                       &
   &  1.455E-10,  1.455E-10,  1.434E-10,  1.381E-10,  1.301E-10,      &
   &  1.276E-10,  1.163E-10,  1.089E-10,  9.911E-11,  8.943E-11,      &
   &  7.618E-11,  6.424E-11,  5.717E-11,  4.866E-11,  4.257E-11,      &
   &  3.773E-11,  3.331E-11,  2.958E-11,  2.629E-11,  2.316E-11,      &
   &  2.073E-11,  1.841E-11,  1.635E-11,  1.464E-11,  1.310E-11,      &
   &  1.160E-11,  1.047E-11,  9.408E-12,  8.414E-12,  7.521E-12,      &
   &  6.705E-12,  5.993E-12,  5.371E-12,  4.815E-12,  4.338E-12,      &
   &  3.921E-12,  3.567E-12,  3.265E-12,  3.010E-12,  2.795E-12,      &
   &  2.613E-12,  2.464E-12,  2.346E-12,  2.256E-12,  2.195E-12,      &
   &  2.165E-12,  2.166E-12,  2.198E-12,  2.262E-12,  2.364E-12/
   DATA S1801/                                                       &
   &  2.502E-12,  2.682E-12,  2.908E-12,  3.187E-12,  3.533E-12,      &
   &  3.946E-12,  4.418E-12,  5.013E-12,  5.708E-12,  6.379E-12,      &
   &  7.430E-12,  8.390E-12,  9.510E-12,  1.078E-11,  1.259E-11,      &
   &  1.438E-11,  1.630E-11,  1.814E-11,  2.055E-11,  2.348E-11,      &
   &  2.664E-11,  2.956E-11,  3.300E-11,  3.677E-11,  4.032E-11,      &
   &  4.494E-11,  4.951E-11,  5.452E-11,  6.014E-11,  6.500E-11,      &
   &  6.915E-11,  7.450E-11,  7.971E-11,  8.468E-11,  8.726E-11,      &
   &  8.995E-11,  9.182E-11,  9.509E-11,  9.333E-11,  9.386E-11,      &
   &  9.457E-11,  9.210E-11,  9.019E-11,  8.680E-11,  8.298E-11,      &
   &  7.947E-11,  7.460E-11,  7.082E-11,  6.132E-11,  5.855E-11/
   DATA S1851/                                                       &
   &  5.073E-11,  4.464E-11,  3.825E-11,  3.375E-11,  2.911E-11,      &
   &  2.535E-11,  2.160E-11,  1.907E-11,  1.665E-11,  1.463E-11,      &
   &  1.291E-11,  1.133E-11,  9.997E-12,  8.836E-12,  7.839E-12,      &
   &  6.943E-12,  6.254E-12,  5.600E-12,  5.029E-12,  4.529E-12,      &
   &  4.102E-12,  3.737E-12,  3.428E-12,  3.169E-12,  2.959E-12,      &
   &  2.798E-12,  2.675E-12,  2.582E-12,  2.644E-12,  2.557E-12,      &
   &  2.614E-12,  2.717E-12,  2.874E-12,  3.056E-12,  3.187E-12,      &
   &  3.631E-12,  3.979E-12,  4.248E-12,  4.817E-12,  5.266E-12,      &
   &  5.836E-12,  6.365E-12,  6.807E-12,  7.470E-12,  7.951E-12,      &
   &  8.636E-12,  8.972E-12,  9.314E-12,  9.445E-12,  1.003E-11/
   DATA S1901/                                                       &
   &  1.013E-11,  9.937E-12,  9.729E-12,  9.064E-12,  9.119E-12,      &
   &  9.124E-12,  8.704E-12,  8.078E-12,  7.470E-12,  6.329E-12,      &
   &  5.674E-12,  4.808E-12,  4.119E-12,  3.554E-12,  3.103E-12,      &
   &  2.731E-12,  2.415E-12,  2.150E-12,  1.926E-12,  1.737E-12,      &
   &  1.578E-12,  1.447E-12,  1.340E-12,  1.255E-12,  1.191E-12,      &
   &  1.146E-12,  1.121E-12,  1.114E-12,  1.126E-12,  1.156E-12,      &
   &  1.207E-12,  1.278E-12,  1.372E-12,  1.490E-12,  1.633E-12,      &
   &  1.805E-12,  2.010E-12,  2.249E-12,  2.528E-12,  2.852E-12,      &
   &  3.228E-12,  3.658E-12,  4.153E-12,  4.728E-12,  5.394E-12,      &
   &  6.176E-12,  7.126E-12,  8.188E-12,  9.328E-12,  1.103E-11/
   DATA S1951/                                                       &
   &  1.276E-11,  1.417E-11,  1.615E-11,  1.840E-11,  2.155E-11,      &
   &  2.429E-11,  2.826E-11,  3.222E-11,  3.664E-11,  4.140E-11,      &
   &  4.906E-11,  5.536E-11,  6.327E-11,  7.088E-11,  8.316E-11,      &
   &  9.242E-11,  1.070E-10,  1.223E-10,  1.341E-10,  1.553E-10,      &
   &  1.703E-10,  1.900E-10,  2.022E-10,  2.233E-10,  2.345E-10,      &
   &  2.438E-10,  2.546E-10,  2.599E-10,  2.661E-10,  2.703E-10,      &
   &  2.686E-10,  2.662E-10,  2.560E-10,  2.552E-10,  2.378E-10,      &
   &  2.252E-10,  2.146E-10,  1.885E-10,  1.668E-10,  1.441E-10,      &
   &  1.295E-10,  1.119E-10,  9.893E-11,  8.687E-11,  7.678E-11,      &
   &  6.685E-11,  5.879E-11,  5.127E-11,  4.505E-11,  3.997E-11/
   DATA S2001/                                                       &
   &  3.511e-11/
!
end block data BS260
!
!     --------------------------------------------------------------
!
SUBROUTINE FRN296 (V1C,V2C,DVC,NPTC,C,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   COMMON /FH2O/ V1S,V2S,DVS,NPTS,S(2003)
   DIMENSION C(*)
!
   DVC = DVS
   v1ss = v1s
   v2ss = v2s
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.GE.1).AND.(I.LE.NPTS)) THEN
         C(J) = S(I)
      ENDIF
10 END DO
!
   RETURN
!
end subroutine FRN296
!
!     --------------------------------------------------------------
!
BLOCK DATA BFH2O
!
   IMPLICIT REAL*8           (V)
!
!               11/18/02
!               UNITS OF (CM**3/MOL)*1.E-20
!
   COMMON /FH2O/ V1,V2,DV,NPT,                                       &
   &              F0000( 2),F0001(50),F0051(50),F0101(50),F0151(50),  &
   &              F0201(50),F0251(50),F0301(50),F0351(50),F0401(50),  &
   &              F0451(50),F0501(50),F0551(50),F0601(50),F0651(50),  &
   &              F0701(50),F0751(50),F0801(50),F0851(50),F0901(50),  &
   &              F0951(50),F1001(50),F1051(50),F1101(50),F1151(50),  &
   &              F1201(50),F1251(50),F1301(50),F1351(50),F1401(50),  &
   &              F1451(50),F1501(50),F1551(50),F1601(50),F1651(50),  &
   &              F1701(50),F1751(50),F1801(50),F1851(50),F1901(50),  &
   &              F1951(50), F2001(1)
!
   DATA V1,V2,DV,NPT / -20.0, 20000.0, 10.0, 2003/
!
   DATA F0000/                                                       &
   &      1.205e-02, 1.126e-02/
   DATA F0001/                                                       &
   &  1.095e-02,  1.126e-02,  1.205e-02,  1.322e-02,  1.430e-02,      &
   &  1.506e-02,  1.548e-02,  1.534e-02,  1.486e-02,  1.373e-02,      &
   &  1.262e-02,  1.134e-02,  1.001e-02,  8.702e-03,  7.475e-03,      &
   &  6.481e-03,  5.480e-03,  4.600e-03,  3.833e-03,  3.110e-03,      &
   &  2.543e-03,  2.049e-03,  1.680e-03,  1.374e-03,  1.046e-03,      &
   &  8.193e-04,  6.267e-04,  4.968e-04,  3.924e-04,  2.983e-04,      &
   &  2.477e-04,  1.997e-04,  1.596e-04,  1.331e-04,  1.061e-04,      &
   &  8.942e-05,  7.168e-05,  5.887e-05,  4.848e-05,  3.817e-05,      &
   &  3.170e-05,  2.579e-05,  2.162e-05,  1.768e-05,  1.490e-05,      &
   &  1.231e-05,  1.013e-05,  8.555e-06,  7.328e-06,  6.148e-06/
   DATA F0051/                                                       &
   &  5.207e-06,  4.387e-06,  3.741e-06,  3.220e-06,  2.753e-06,      &
   &  2.346e-06,  1.985e-06,  1.716e-06,  1.475e-06,  1.286e-06,      &
   &  1.122e-06,  9.661e-07,  8.284e-07,  7.057e-07,  6.119e-07,      &
   &  5.290e-07,  4.571e-07,  3.948e-07,  3.432e-07,  2.983e-07,      &
   &  2.589e-07,  2.265e-07,  1.976e-07,  1.704e-07,  1.456e-07,      &
   &  1.260e-07,  1.101e-07,  9.648e-08,  8.415e-08,  7.340e-08,      &
   &  6.441e-08,  5.643e-08,  4.940e-08,  4.276e-08,  3.703e-08,      &
   &  3.227e-08,  2.825e-08,  2.478e-08,  2.174e-08,  1.898e-08,      &
   &  1.664e-08,  1.458e-08,  1.278e-08,  1.126e-08,  9.891e-09,      &
   &  8.709e-09,  7.652e-09,  6.759e-09,  5.975e-09,  5.310e-09/
   DATA F0101/                                                       &
   &  4.728e-09,  4.214e-09,  3.792e-09,  3.463e-09,  3.226e-09,      &
   &  2.992e-09,  2.813e-09,  2.749e-09,  2.809e-09,  2.913e-09,      &
   &  3.037e-09,  3.413e-09,  3.738e-09,  4.189e-09,  4.808e-09,      &
   &  5.978e-09,  7.088e-09,  8.071e-09,  9.610e-09,  1.210e-08,      &
   &  1.500e-08,  1.764e-08,  2.221e-08,  2.898e-08,  3.948e-08,      &
   &  5.068e-08,  6.227e-08,  7.898e-08,  1.033e-07,  1.437e-07,      &
   &  1.889e-07,  2.589e-07,  3.590e-07,  4.971e-07,  7.156e-07,      &
   &  9.983e-07,  1.381e-06,  1.929e-06,  2.591e-06,  3.453e-06,      &
   &  4.570e-06,  5.930e-06,  7.552e-06,  9.556e-06,  1.183e-05,      &
   &  1.425e-05,  1.681e-05,  1.978e-05,  2.335e-05,  2.668e-05/
   DATA F0151/                                                          &
   &  3.022e-05,  3.371e-05,  3.715e-05,  3.967e-05,  4.060e-05,      &
   &  4.010e-05,  3.809e-05,  3.491e-05,  3.155e-05,  2.848e-05,      &
   &  2.678e-05,  2.660e-05,  2.811e-05,  3.071e-05,  3.294e-05,      &
   &  3.459e-05,  3.569e-05,  3.560e-05,  3.434e-05,  3.186e-05,      &
   &  2.916e-05,  2.622e-05,  2.275e-05,  1.918e-05,  1.620e-05,      &
   &  1.373e-05,  1.182e-05,  1.006e-05,  8.556e-06,  7.260e-06,      &
   &  6.110E-06,  5.040E-06,  4.230E-06,  3.456E-06,  2.903E-06,      &
   &  2.486E-06,  2.039E-06,  1.672E-06,  1.285E-06,  1.054E-06,      &
   &  8.302E-07,  6.667E-07,  5.503E-07,  4.562E-07,  3.948E-07,      &
   &  3.198E-07,  2.586E-07,  2.225E-07,  1.807E-07,  1.530E-07/
   DATA F0201/                                                       &
   &  1.294E-07,  1.126E-07,  9.604E-08,  7.850E-08,  6.813E-08,      &
   &  5.583E-08,  4.690E-08,  3.996E-08,  3.373E-08,  2.930E-08,      &
   &  2.417E-08,  2.061E-08,  1.743E-08,  1.475E-08,  1.273E-08,      &
   &  1.084E-08,  9.368E-09,  7.985E-09,  6.785E-09,  5.804E-09,      &
   &  4.975E-09,  4.311E-09,  3.738E-09,  3.275E-09,  2.847E-09,      &
   &  2.469E-09,  2.149E-09,  1.884E-09,  1.631E-09,  1.393E-09,      &
   &  1.201E-09,  1.027E-09,  8.807E-10,  7.521E-10,  6.436E-10,      &
   &  5.509E-10,  4.729E-10,  4.055E-10,  3.490E-10,  3.006E-10,      &
   &  2.599E-10,  2.246E-10,  1.942E-10,  1.685E-10,  1.464E-10,      &
   &  1.273E-10,  1.115E-10,  9.794E-11,  8.729E-11,  7.893E-11/
   DATA F0251/                                                       &
   &  7.313E-11,  7.069E-11,  7.190E-11,  7.828E-11,  9.295E-11,      &
   &  1.174E-10,  1.578E-10,  2.184E-10,  3.053E-10,  4.212E-10,      &
   &  5.733E-10,  7.497E-10,  9.487E-10,  1.153E-09,  1.343E-09,      &
   &  1.503E-09,  1.623E-09,  1.696E-09,  1.746E-09,  1.789E-09,      &
   &  1.852E-09,  1.862E-09,  1.886E-09,  1.912E-09,  1.942E-09,      &
   &  1.938E-09,  1.914E-09,  1.898E-09,  1.864E-09,  1.770E-09,      &
   &  1.635E-09,  1.464E-09,  1.279E-09,  1.099E-09,  9.528E-10,      &
   &  8.433E-10,  7.794E-10,  7.793E-10,  8.241E-10,  9.417E-10,      &
   &  1.179E-09,  1.616E-09,  2.144E-09,  2.992E-09,  4.344E-09,      &
   &  6.415E-09,  9.242E-09,  1.310E-08,  1.847E-08,  2.567E-08/
   DATA F0301/                                                       &
   &  3.390E-08,  4.357E-08,  5.301E-08,  6.364E-08,  7.438E-08,      &
   &  8.381E-08,  9.294E-08,  1.013E-07,  1.103E-07,  1.170E-07,      &
   &  1.200E-07,  1.190E-07,  1.133E-07,  1.040E-07,  9.475E-08,      &
   &  8.601E-08,  8.074E-08,  8.023E-08,  8.473E-08,  9.256E-08,      &
   &  9.989E-08,  1.056E-07,  1.079E-07,  1.079E-07,  1.049E-07,      &
   &  1.000E-07,  9.350E-08,  8.662E-08,  7.880E-08,  7.137E-08,      &
   &  6.438E-08,  5.706E-08,  5.132E-08,  4.815E-08,  4.736E-08,      &
   &  4.709E-08,  4.918E-08,  5.227E-08,  5.603E-08,  6.023E-08,      &
   &  6.567E-08,  7.403E-08,  8.695E-08,  9.929E-08,  1.163E-07,      &
   &  1.379E-07,  1.658E-07,  2.010E-07,  2.425E-07,  2.920E-07/
   DATA F0351/                                                       &
   &  3.562E-07,  4.274E-07,  5.173E-07,  6.285E-07,  7.787E-07,      &
   &  9.563E-07,  1.194E-06,  1.517E-06,  1.934E-06,  2.511E-06,      &
   &  3.197E-06,  4.024E-06,  4.981E-06,  6.090E-06,  7.229E-06,      &
   &  8.439E-06,  9.367E-06,  1.016E-05,  1.057E-05,  1.088E-05,      &
   &  1.108E-05,  1.129E-05,  1.159E-05,  1.190E-05,  1.210E-05,      &
   &  1.210E-05,  1.200E-05,  1.190E-05,  1.180E-05,  1.190E-05,      &
   &  1.230E-05,  1.280E-05,  1.300E-05,  1.290E-05,  1.210E-05,      &
   &  1.090E-05,  9.440E-06,  7.790E-06,  6.220E-06,  4.640E-06,      &
   &  3.400E-06,  2.470E-06,  1.780E-06,  1.290E-06,  9.260E-07,      &
   &  7.000E-07,  5.360E-07,  4.260E-07,  3.410E-07,  2.790E-07/
   DATA F0401/                                                       &
   &  2.360E-07,  2.010E-07,  1.750E-07,  1.520E-07,  1.330E-07,      &
   &  1.160E-07,  1.010E-07,  8.830E-08,  7.710E-08,  6.770E-08,      &
   &  6.010E-08,  5.340E-08,  4.720E-08,  4.200E-08,  3.690E-08,      &
   &  3.270E-08,  2.920E-08,  2.570E-08,  2.270E-08,  2.000E-08,      &
   &  1.770E-08,  1.560E-08,  1.370E-08,  1.200E-08,  1.050E-08,      &
   &  9.180E-09,  8.000E-09,  6.980E-09,  6.080E-09,  5.270E-09,      &
   &  4.560E-09,  3.940E-09,  3.400E-09,  2.930E-09,  2.510E-09,      &
   &  2.150E-09,  1.840E-09,  1.570E-09,  1.340E-09,  1.140E-09,      &
   &  9.660E-10,  8.190E-10,  6.930E-10,  5.860E-10,  4.950E-10,      &
   &  4.190E-10,  3.560E-10,  3.060E-10,  2.670E-10,  2.410E-10/
   DATA F0451/                                                       &
   &  2.260E-10,  2.250E-10,  2.360E-10,  2.570E-10,  2.840E-10,      &
   &  3.180E-10,  3.490E-10,  3.710E-10,  3.810E-10,  3.840E-10,      &
   &  3.800E-10,  3.590E-10,  3.370E-10,  3.080E-10,  2.770E-10,      &
   &  2.470E-10,  2.200E-10,  2.050E-10,  1.990E-10,  2.030E-10,      &
   &  2.130E-10,  2.250E-10,  2.340E-10,  2.410E-10,  2.470E-10,      &
   &  2.480E-10,  2.500E-10,  2.520E-10,  2.580E-10,  2.640E-10,      &
   &  2.690E-10,  2.790E-10,  2.950E-10,  3.130E-10,  3.340E-10,      &
   &  3.570E-10,  3.840E-10,  4.200E-10,  4.660E-10,  5.210E-10,      &
   &  5.960E-10,  6.890E-10,  8.010E-10,  9.320E-10,  1.080E-09,      &
   &  1.240E-09,  1.450E-09,  1.710E-09,  2.010E-09,  2.360E-09/
   DATA F0501/                                                       &
   &  2.820E-09,  3.380E-09,  4.010E-09,  4.790E-09,  5.740E-09,      &
   &  6.920E-09,  8.338E-09,  1.015E-08,  1.246E-08,  1.512E-08,      &
   &  1.851E-08,  2.320E-08,  2.988E-08,  3.989E-08,  5.351E-08,      &
   &  7.256E-08,  1.011E-07,  1.388E-07,  1.914E-07,  2.512E-07,      &
   &  3.257E-07,  4.075E-07,  4.991E-07,  5.914E-07,  6.638E-07,      &
   &  7.160E-07,  7.450E-07,  7.590E-07,  7.540E-07,  7.470E-07,      &
   &  7.530E-07,  7.940E-07,  8.460E-07,  8.950E-07,  9.360E-07,      &
   &  9.840E-07,  1.020E-06,  1.050E-06,  1.060E-06,  1.070E-06,      &
   &  1.070E-06,  1.040E-06,  9.790E-07,  8.800E-07,  7.660E-07,      &
   &  6.380E-07,  5.142E-07,  3.982E-07,  2.940E-07,  2.177E-07/
   DATA F0551/                                                       &
   &  1.612E-07,  1.185E-07,  8.966E-08,  6.735E-08,  5.219E-08,      &
   &  4.140E-08,  3.403E-08,  2.811E-08,  2.358E-08,  1.969E-08,      &
   &  1.684E-08,  1.462E-08,  1.285E-08,  1.134E-08,  9.916E-09,      &
   &  8.650E-09,  7.550E-09,  6.620E-09,  5.830E-09,  5.140E-09,      &
   &  4.510E-09,  3.980E-09,  3.530E-09,  3.120E-09,  2.780E-09,      &
   &  2.480E-09,  2.180E-09,  1.930E-09,  1.700E-09,  1.500E-09,      &
   &  1.330E-09,  1.160E-09,  1.030E-09,  9.010E-10,  7.860E-10,      &
   &  6.880E-10,  6.020E-10,  5.250E-10,  4.570E-10,  3.960E-10,      &
   &  3.440E-10,  2.980E-10,  2.570E-10,  2.210E-10,  1.900E-10,      &
   &  1.630E-10,  1.400E-10,  1.210E-10,  1.040E-10,  9.090E-11/
   DATA F0601/                                                       &
   &  7.960E-11,  7.050E-11,  6.300E-11,  5.670E-11,  5.090E-11,      &
   &  4.550E-11,  4.040E-11,  3.550E-11,  3.090E-11,  2.640E-11,      &
   &  2.260E-11,  1.900E-11,  1.600E-11,  1.340E-11,  1.140E-11,      &
   &  1.010E-11,  9.170E-12,  8.690E-12,  8.470E-12,  8.330E-12,      &
   &  8.200E-12,  8.040E-12,  7.950E-12,  7.980E-12,  8.250E-12,      &
   &  8.470E-12,  8.740E-12,  8.870E-12,  8.920E-12,  8.710E-12,      &
   &  8.300E-12,  7.890E-12,  7.590E-12,  7.500E-12,  7.630E-12,      &
   &  8.060E-12,  8.890E-12,  1.000E-11,  1.140E-11,  1.280E-11,      &
   &  1.450E-11,  1.640E-11,  1.870E-11,  2.120E-11,  2.440E-11,      &
   &  2.850E-11,  3.310E-11,  3.830E-11,  4.460E-11,  5.230E-11/
   DATA F0651/                                                       &
   &  6.150E-11,  7.350E-11,  8.840E-11,  1.070E-10,  1.290E-10,      &
   &  1.590E-10,  1.990E-10,  2.490E-10,  3.120E-10,  3.910E-10,      &
   &  4.980E-10,  6.280E-10,  7.860E-10,  9.660E-10,  1.160E-09,      &
   &  1.410E-09,  1.700E-09,  2.060E-09,  2.540E-09,  3.170E-09,      &
   &  4.080E-09,  5.240E-09,  6.800E-09,  8.920E-09,  1.130E-08,      &
   &  1.400E-08,  1.700E-08,  2.050E-08,  2.370E-08,  2.610E-08,      &
   &  2.770E-08,  2.850E-08,  2.860E-08,  2.780E-08,  2.720E-08,      &
   &  2.740E-08,  2.830E-08,  3.010E-08,  3.220E-08,  3.520E-08,      &
   &  3.780E-08,  4.030E-08,  4.210E-08,  4.330E-08,  4.390E-08,      &
   &  4.340E-08,  4.250E-08,  4.070E-08,  3.880E-08,  3.670E-08/
   DATA F0701/                                                       &
   &  3.560E-08,  3.520E-08,  3.670E-08,  4.020E-08,  4.640E-08,      &
   &  5.630E-08,  6.890E-08,  8.820E-08,  1.110E-07,  1.410E-07,      &
   &  1.730E-07,  2.120E-07,  2.570E-07,  2.990E-07,  3.450E-07,      &
   &  3.790E-07,  4.090E-07,  4.230E-07,  4.300E-07,  4.330E-07,      &
   &  4.350E-07,  4.440E-07,  4.540E-07,  4.620E-07,  4.670E-07,      &
   &  4.710E-07,  4.770E-07,  4.930E-07,  5.190E-07,  5.540E-07,      &
   &  5.900E-07,  6.120E-07,  6.150E-07,  5.810E-07,  5.210E-07,      &
   &  4.510E-07,  3.690E-07,  2.890E-07,  2.140E-07,  1.550E-07,      &
   &  1.110E-07,  7.780E-08,  5.460E-08,  3.930E-08,  2.930E-08,      &
   &  2.260E-08,  1.820E-08,  1.520E-08,  1.310E-08,  1.150E-08/
   DATA F0751/                                                       &
   &  1.020E-08,  9.010E-09,  7.950E-09,  7.040E-09,  6.200E-09,      &
   &  5.490E-09,  4.860E-09,  4.340E-09,  3.880E-09,  3.460E-09,      &
   &  3.070E-09,  2.710E-09,  2.380E-09,  2.090E-09,  1.820E-09,      &
   &  1.590E-09,  1.380E-09,  1.200E-09,  1.040E-09,  9.080E-10,      &
   &  7.880E-10,  6.850E-10,  5.960E-10,  5.200E-10,  4.530E-10,      &
   &  3.940E-10,  3.430E-10,  2.970E-10,  2.580E-10,  2.230E-10,      &
   &  1.920E-10,  1.660E-10,  1.430E-10,  1.240E-10,  1.060E-10,      &
   &  9.100E-11,  7.790E-11,  6.650E-11,  5.680E-11,  4.840E-11,      &
   &  4.110E-11,  3.510E-11,  2.990E-11,  2.540E-11,  2.170E-11,      &
   &  1.860E-11,  1.590E-11,  1.370E-11,  1.180E-11,  1.040E-11/
   DATA F0801/                                                       &
   &  9.240E-12,  8.390E-12,  7.860E-12,  7.700E-12,  7.850E-12,      &
   &  8.410E-12,  9.620E-12,  1.160E-11,  1.470E-11,  1.940E-11,      &
   &  2.580E-11,  3.510E-11,  4.730E-11,  6.260E-11,  7.980E-11,      &
   &  9.950E-11,  1.210E-10,  1.400E-10,  1.580E-10,  1.770E-10,      &
   &  1.980E-10,  2.220E-10,  2.470E-10,  2.790E-10,  3.140E-10,      &
   &  3.550E-10,  4.000E-10,  4.510E-10,  5.150E-10,  5.730E-10,      &
   &  6.250E-10,  6.650E-10,  6.910E-10,  6.980E-10,  6.830E-10,      &
   &  6.670E-10,  6.570E-10,  6.600E-10,  6.840E-10,  7.390E-10,      &
   &  8.040E-10,  8.730E-10,  9.440E-10,  9.950E-10,  1.030E-09,      &
   &  1.010E-09,  9.800E-10,  9.320E-10,  8.670E-10,  7.990E-10/
   DATA F0851/                                                       &
   &  7.290E-10,  6.670E-10,  6.140E-10,  5.780E-10,  5.650E-10,      &
   &  5.650E-10,  5.880E-10,  6.430E-10,  7.500E-10,  9.160E-10,      &
   &  1.160E-09,  1.510E-09,  1.990E-09,  2.700E-09,  3.660E-09,      &
   &  4.950E-09,  6.410E-09,  8.160E-09,  1.020E-08,  1.240E-08,      &
   &  1.460E-08,  1.640E-08,  1.790E-08,  1.870E-08,  1.930E-08,      &
   &  1.930E-08,  1.930E-08,  1.950E-08,  2.070E-08,  2.230E-08,      &
   &  2.340E-08,  2.440E-08,  2.550E-08,  2.670E-08,  2.760E-08,      &
   &  2.810E-08,  2.900E-08,  2.940E-08,  2.910E-08,  2.770E-08,      &
   &  2.520E-08,  2.210E-08,  1.840E-08,  1.490E-08,  1.150E-08,      &
   &  8.640E-09,  6.440E-09,  4.790E-09,  3.640E-09,  2.750E-09/
   DATA F0901/                                                       &
   &  2.140E-09,  1.750E-09,  1.480E-09,  1.300E-09,  1.160E-09,      &
   &  1.050E-09,  9.460E-10,  8.550E-10,  7.660E-10,  6.790E-10,      &
   &  6.010E-10,  5.300E-10,  4.670E-10,  4.130E-10,  3.670E-10,      &
   &  3.290E-10,  2.920E-10,  2.570E-10,  2.270E-10,  1.960E-10,      &
   &  1.700E-10,  1.470E-10,  1.270E-10,  1.100E-10,  9.190E-11,      &
   &  7.790E-11,  6.590E-11,  5.610E-11,  4.810E-11,  4.070E-11,      &
   &  3.540E-11,  3.040E-11,  2.620E-11,  2.300E-11,  1.970E-11,      &
   &  1.710E-11,  1.480E-11,  1.280E-11,  1.100E-11,  9.440E-12,      &
   &  8.080E-12,  6.940E-12,  5.930E-12,  5.060E-12,  4.310E-12,      &
   &  3.660E-12,  3.120E-12,  2.670E-12,  2.280E-12,  1.950E-12/
   DATA F0951/                                                       &
   &  1.670E-12,  1.430E-12,  1.220E-12,  1.050E-12,  9.060E-13,      &
   &  7.860E-13,  6.870E-13,  6.090E-13,  5.500E-13,  5.070E-13,      &
   &  4.780E-13,  4.660E-13,  4.770E-13,  5.210E-13,  5.940E-13,      &
   &  7.440E-13,  9.710E-13,  1.350E-12,  1.870E-12,  2.630E-12,      &
   &  3.670E-12,  4.950E-12,  6.390E-12,  8.000E-12,  9.880E-12,      &
   &  1.150E-11,  1.280E-11,  1.370E-11,  1.430E-11,  1.470E-11,      &
   &  1.450E-11,  1.420E-11,  1.400E-11,  1.410E-11,  1.430E-11,      &
   &  1.510E-11,  1.650E-11,  1.830E-11,  2.070E-11,  2.250E-11,      &
   &  2.380E-11,  2.420E-11,  2.380E-11,  2.270E-11,  2.100E-11,      &
   &  1.980E-11,  1.890E-11,  1.790E-11,  1.720E-11,  1.620E-11/
   DATA F1001/                                                       &
   &  1.540E-11,  1.430E-11,  1.350E-11,  1.350E-11,  1.420E-11,      &
   &  1.570E-11,  1.770E-11,  2.060E-11,  2.420E-11,  2.850E-11,      &
   &  3.430E-11,  4.160E-11,  5.240E-11,  6.820E-11,  8.920E-11,      &
   &  1.200E-10,  1.600E-10,  2.150E-10,  2.890E-10,  3.770E-10,      &
   &  4.800E-10,  5.930E-10,  7.250E-10,  8.530E-10,  9.600E-10,      &
   &  1.040E-09,  1.080E-09,  1.100E-09,  1.090E-09,  1.070E-09,      &
   &  1.070E-09,  1.120E-09,  1.210E-09,  1.330E-09,  1.480E-09,      &
   &  1.640E-09,  1.790E-09,  1.940E-09,  2.050E-09,  2.170E-09,      &
   &  2.280E-09,  2.400E-09,  2.530E-09,  2.740E-09,  3.030E-09,      &
   &  3.420E-09,  3.870E-09,  4.470E-09,  5.190E-09,  5.910E-09/
   DATA F1051/                                                       &
   &  6.770E-09,  7.500E-09,  8.280E-09,  8.830E-09,  9.360E-09,      &
   &  9.910E-09,  1.030E-08,  1.090E-08,  1.140E-08,  1.170E-08,      &
   &  1.180E-08,  1.180E-08,  1.180E-08,  1.210E-08,  1.260E-08,      &
   &  1.350E-08,  1.430E-08,  1.490E-08,  1.490E-08,  1.400E-08,      &
   &  1.250E-08,  1.070E-08,  8.610E-09,  6.570E-09,  4.800E-09,      &
   &  3.470E-09,  2.470E-09,  1.780E-09,  1.330E-09,  1.050E-09,      &
   &  8.720E-10,  7.510E-10,  6.760E-10,  6.270E-10,  6.000E-10,      &
   &  5.990E-10,  6.260E-10,  6.860E-10,  7.580E-10,  8.490E-10,      &
   &  9.530E-10,  1.050E-09,  1.140E-09,  1.200E-09,  1.250E-09,      &
   &  1.280E-09,  1.300E-09,  1.310E-09,  1.320E-09,  1.330E-09/
   DATA F1101/                                                       &
   &  1.330E-09,  1.310E-09,  1.270E-09,  1.230E-09,  1.190E-09,      &
   &  1.200E-09,  1.250E-09,  1.330E-09,  1.420E-09,  1.480E-09,      &
   &  1.490E-09,  1.410E-09,  1.260E-09,  1.080E-09,  8.800E-10,      &
   &  6.780E-10,  4.950E-10,  3.580E-10,  2.540E-10,  1.780E-10,      &
   &  1.260E-10,  9.280E-11,  7.130E-11,  5.560E-11,  4.500E-11,      &
   &  3.750E-11,  3.210E-11,  2.760E-11,  2.410E-11,  2.090E-11,      &
   &  1.820E-11,  1.590E-11,  1.360E-11,  1.180E-11,  1.030E-11,      &
   &  9.010E-12,  7.940E-12,  6.900E-12,  6.020E-12,  5.260E-12,      &
   &  4.620E-12,  4.060E-12,  3.550E-12,  3.130E-12,  2.760E-12,      &
   &  2.430E-12,  2.150E-12,  1.900E-12,  1.680E-12,  1.490E-12/
   DATA F1151/                                                       &
   &  1.320E-12,  1.170E-12,  1.040E-12,  9.350E-13,  8.470E-13,      &
   &  7.790E-13,  7.300E-13,  6.960E-13,  6.880E-13,  7.110E-13,      &
   &  7.600E-13,  8.670E-13,  1.040E-12,  1.310E-12,  1.710E-12,      &
   &  2.200E-12,  2.930E-12,  3.800E-12,  4.970E-12,  6.220E-12,      &
   &  7.630E-12,  9.290E-12,  1.100E-11,  1.240E-11,  1.340E-11,      &
   &  1.410E-11,  1.440E-11,  1.440E-11,  1.390E-11,  1.380E-11,      &
   &  1.410E-11,  1.490E-11,  1.640E-11,  1.860E-11,  2.110E-11,      &
   &  2.380E-11,  2.650E-11,  2.890E-11,  3.040E-11,  3.140E-11,      &
   &  3.290E-11,  3.430E-11,  3.600E-11,  3.910E-11,  4.440E-11,      &
   &  5.210E-11,  6.430E-11,  8.180E-11,  1.060E-10,  1.400E-10/
   DATA F1201/                                                       &
   &  1.780E-10,  2.240E-10,  2.770E-10,  3.350E-10,  3.950E-10,      &
   &  4.470E-10,  4.920E-10,  5.200E-10,  5.390E-10,  5.460E-10,      &
   &  5.530E-10,  5.660E-10,  6.020E-10,  6.430E-10,  6.680E-10,      &
   &  6.910E-10,  7.120E-10,  7.420E-10,  7.550E-10,  7.750E-10,      &
   &  8.090E-10,  8.300E-10,  8.350E-10,  7.950E-10,  7.220E-10,      &
   &  6.300E-10,  5.230E-10,  4.200E-10,  3.190E-10,  2.390E-10,      &
   &  1.790E-10,  1.360E-10,  1.050E-10,  8.320E-11,  6.890E-11,      &
   &  5.890E-11,  5.140E-11,  4.600E-11,  4.190E-11,  3.940E-11,      &
   &  3.900E-11,  4.070E-11,  4.390E-11,  4.820E-11,  5.290E-11,      &
   &  5.750E-11,  6.090E-11,  6.280E-11,  6.310E-11,  6.190E-11/
   DATA F1251/                                                       &
   &  6.040E-11,  5.850E-11,  5.710E-11,  5.680E-11,  5.680E-11,      &
   &  5.640E-11,  5.580E-11,  5.510E-11,  5.480E-11,  5.450E-11,      &
   &  5.570E-11,  5.770E-11,  5.960E-11,  6.050E-11,  5.890E-11,      &
   &  5.430E-11,  4.770E-11,  4.020E-11,  3.250E-11,  2.500E-11,      &
   &  1.820E-11,  1.310E-11,  9.320E-12,  6.580E-12,  4.640E-12,      &
   &  3.350E-12,  2.530E-12,  1.960E-12,  1.580E-12,  1.310E-12,      &
   &  1.110E-12,  9.580E-13,  8.320E-13,  7.310E-13,  6.430E-13,      &
   &  5.620E-13,  4.890E-13,  4.270E-13,  3.730E-13,  3.250E-13,      &
   &  2.840E-13,  2.500E-13,  2.210E-13,  1.950E-13,  1.720E-13,      &
   &  1.520E-13,  1.340E-13,  1.180E-13,  1.050E-13,  9.260E-14/
   DATA F1301/                                                       &
   &  8.210E-14,  7.290E-14,  6.480E-14,  5.770E-14,  5.150E-14,      &
   &  4.620E-14,  4.170E-14,  3.790E-14,  3.490E-14,  3.260E-14,      &
   &  3.090E-14,  3.000E-14,  2.980E-14,  3.040E-14,  3.200E-14,      &
   &  3.490E-14,  3.960E-14,  4.680E-14,  5.780E-14,  7.420E-14,      &
   &  9.800E-14,  1.320E-13,  1.790E-13,  2.230E-13,  2.650E-13,      &
   &  3.020E-13,  3.380E-13,  3.740E-13,  4.050E-13,  4.420E-13,      &
   &  4.940E-13,  5.480E-13,  6.110E-13,  6.770E-13,  7.560E-13,      &
   &  8.490E-13,  9.590E-13,  1.090E-12,  1.270E-12,  1.480E-12,      &
   &  1.730E-12,  2.080E-12,  2.480E-12,  3.030E-12,  3.880E-12,      &
   &  4.860E-12,  6.250E-12,  8.050E-12,  1.040E-11,  1.370E-11/
   DATA F1351/                                                       &
   &  1.750E-11,  2.200E-11,  2.750E-11,  3.350E-11,  4.080E-11,      &
   &  4.790E-11,  5.440E-11,  5.990E-11,  6.460E-11,  6.820E-11,      &
   &  7.130E-11,  7.620E-11,  8.360E-11,  9.590E-11,  1.150E-10,      &
   &  1.380E-10,  1.670E-10,  1.990E-10,  2.350E-10,  2.740E-10,      &
   &  3.110E-10,  3.510E-10,  3.890E-10,  4.210E-10,  4.490E-10,      &
   &  4.660E-10,  4.860E-10,  4.960E-10,  5.080E-10,  5.200E-10,      &
   &  5.290E-10,  5.300E-10,  5.220E-10,  5.160E-10,  5.180E-10,      &
   &  5.350E-10,  5.680E-10,  6.070E-10,  6.420E-10,  6.600E-10,      &
   &  6.390E-10,  5.830E-10,  5.030E-10,  4.110E-10,  3.200E-10,      &
   &  2.350E-10,  1.700E-10,  1.200E-10,  8.540E-11,  6.310E-11/
   DATA F1401/                                                       &
   &  4.920E-11,  4.010E-11,  3.380E-11,  2.970E-11,  2.710E-11,      &
   &  2.550E-11,  2.450E-11,  2.420E-11,  2.430E-11,  2.490E-11,      &
   &  2.600E-11,  2.730E-11,  2.910E-11,  3.160E-11,  3.450E-11,      &
   &  3.810E-11,  4.140E-11,  4.540E-11,  5.020E-11,  5.480E-11,      &
   &  6.110E-11,  6.720E-11,  7.370E-11,  7.990E-11,  8.500E-11,      &
   &  9.010E-11,  9.300E-11,  9.520E-11,  9.640E-11,  9.590E-11,      &
   &  9.390E-11,  9.140E-11,  8.950E-11,  8.970E-11,  9.300E-11,      &
   &  9.910E-11,  1.060E-10,  1.120E-10,  1.140E-10,  1.100E-10,      &
   &  1.000E-10,  8.680E-11,  7.140E-11,  5.610E-11,  4.170E-11,      &
   &  3.050E-11,  2.180E-11,  1.550E-11,  1.110E-11,  8.290E-12/
   DATA F1451/                                                       &
   &  6.380E-12,  5.070E-12,  4.200E-12,  3.640E-12,  3.280E-12,      &
   &  3.030E-12,  2.900E-12,  2.780E-12,  2.670E-12,  2.520E-12,      &
   &  2.340E-12,  2.150E-12,  1.950E-12,  1.750E-12,  1.540E-12,      &
   &  1.320E-12,  1.130E-12,  9.550E-13,  7.900E-13,  6.440E-13,      &
   &  5.180E-13,  4.140E-13,  3.320E-13,  2.670E-13,  2.180E-13,      &
   &  1.830E-13,  1.560E-13,  1.350E-13,  1.180E-13,  1.050E-13,      &
   &  9.310E-14,  8.350E-14,  7.550E-14,  6.890E-14,  6.360E-14,      &
   &  5.960E-14,  5.690E-14,  5.560E-14,  5.600E-14,  5.860E-14,      &
   &  6.420E-14,  7.330E-14,  8.850E-14,  1.140E-13,  1.490E-13,      &
   &  1.970E-13,  2.690E-13,  3.610E-13,  4.790E-13,  6.040E-13/
   DATA F1501/                                                       &
   &  7.510E-13,  9.200E-13,  1.090E-12,  1.260E-12,  1.390E-12,      &
   &  1.530E-12,  1.640E-12,  1.710E-12,  1.770E-12,  1.880E-12,      &
   &  2.080E-12,  2.370E-12,  2.820E-12,  3.420E-12,  4.180E-12,      &
   &  5.180E-12,  6.470E-12,  7.980E-12,  9.970E-12,  1.220E-11,      &
   &  1.470E-11,  1.750E-11,  2.050E-11,  2.380E-11,  2.660E-11,      &
   &  2.910E-11,  3.090E-11,  3.210E-11,  3.290E-11,  3.310E-11,      &
   &  3.370E-11,  3.490E-11,  3.690E-11,  3.820E-11,  3.920E-11,      &
   &  4.050E-11,  4.170E-11,  4.300E-11,  4.420E-11,  4.610E-11,      &
   &  4.780E-11,  4.830E-11,  4.690E-11,  4.310E-11,  3.810E-11,      &
   &  3.210E-11,  2.610E-11,  2.040E-11,  1.560E-11,  1.200E-11/
   DATA F1551/                                                       &
   &  9.200E-12,  7.100E-12,  5.650E-12,  4.630E-12,  3.860E-12,      &
   &  3.270E-12,  2.830E-12,  2.480E-12,  2.220E-12,  2.030E-12,      &
   &  1.910E-12,  1.840E-12,  1.820E-12,  1.820E-12,  1.860E-12,      &
   &  1.940E-12,  2.110E-12,  2.360E-12,  2.600E-12,  2.940E-12,      &
   &  3.310E-12,  3.730E-12,  4.170E-12,  4.570E-12,  4.980E-12,      &
   &  5.260E-12,  5.540E-12,  5.790E-12,  6.050E-12,  6.380E-12,      &
   &  6.710E-12,  6.930E-12,  7.030E-12,  7.050E-12,  7.050E-12,      &
   &  7.160E-12,  7.370E-12,  7.790E-12,  8.300E-12,  8.650E-12,      &
   &  8.760E-12,  8.320E-12,  7.510E-12,  6.480E-12,  5.310E-12,      &
   &  4.180E-12,  3.070E-12,  2.220E-12,  1.570E-12,  1.090E-12/
   DATA F1601/                                                       &
   &  7.390E-13,  5.110E-13,  3.660E-13,  2.680E-13,  2.000E-13,      &
   &  1.560E-13,  1.270E-13,  1.060E-13,  9.050E-14,  7.850E-14,      &
   &  6.900E-14,  6.100E-14,  5.420E-14,  4.840E-14,  4.320E-14,      &
   &  3.870E-14,  3.470E-14,  3.120E-14,  2.800E-14,  2.520E-14,      &
   &  2.270E-14,  2.050E-14,  1.860E-14,  1.690E-14,  1.550E-14,      &
   &  1.430E-14,  1.330E-14,  1.250E-14,  1.200E-14,  1.160E-14,      &
   &  1.150E-14,  1.160E-14,  1.190E-14,  1.250E-14,  1.340E-14,      &
   &  1.450E-14,  1.600E-14,  1.780E-14,  2.010E-14,  2.270E-14,      &
   &  2.590E-14,  2.960E-14,  3.400E-14,  3.920E-14,  4.530E-14,      &
   &  5.240E-14,  6.080E-14,  7.090E-14,  8.310E-14,  9.810E-14/
   DATA F1651/                                                       &
   &  1.170E-13,  1.420E-13,  1.780E-13,  2.250E-13,  2.790E-13,      &
   &  3.500E-13,  4.400E-13,  5.530E-13,  6.920E-13,  8.920E-13,      &
   &  1.150E-12,  1.480E-12,  1.920E-12,  2.540E-12,  3.370E-12,      &
   &  4.370E-12,  5.710E-12,  7.260E-12,  8.980E-12,  1.100E-11,      &
   &  1.340E-11,  1.590E-11,  1.860E-11,  2.140E-11,  2.400E-11,      &
   &  2.690E-11,  2.980E-11,  3.330E-11,  3.770E-11,  4.250E-11,      &
   &  4.750E-11,  5.110E-11,  5.370E-11,  5.560E-11,  5.610E-11,      &
   &  5.670E-11,  5.680E-11,  5.780E-11,  5.860E-11,  5.910E-11,      &
   &  5.930E-11,  5.910E-11,  5.920E-11,  5.930E-11,  5.940E-11,      &
   &  5.870E-11,  5.650E-11,  5.210E-11,  4.540E-11,  3.820E-11/
   DATA F1701/                                                       &
   &  3.070E-11,  2.370E-11,  1.780E-11,  1.320E-11,  9.820E-12,      &
   &  7.340E-12,  5.590E-12,  4.330E-12,  3.420E-12,  2.740E-12,      &
   &  2.230E-12,  1.860E-12,  1.560E-12,  1.350E-12,  1.170E-12,      &
   &  1.050E-12,  9.540E-13,  8.830E-13,  8.340E-13,  8.100E-13,      &
   &  8.110E-13,  8.200E-13,  8.460E-13,  8.830E-13,  9.330E-13,      &
   &  9.920E-13,  1.070E-12,  1.190E-12,  1.350E-12,  1.580E-12,      &
   &  1.860E-12,  2.170E-12,  2.520E-12,  2.920E-12,  3.360E-12,      &
   &  3.830E-12,  4.300E-12,  4.830E-12,  5.370E-12,  5.810E-12,      &
   &  6.230E-12,  6.480E-12,  6.770E-12,  6.990E-12,  7.170E-12,      &
   &  7.470E-12,  7.660E-12,  7.770E-12,  7.770E-12,  7.720E-12/
   DATA F1751/                                                       &
   &  7.770E-12,  7.970E-12,  8.400E-12,  8.900E-12,  9.420E-12,      &
   &  9.710E-12,  9.550E-12,  8.850E-12,  7.680E-12,  6.390E-12,      &
   &  5.020E-12,  3.740E-12,  2.680E-12,  1.870E-12,  1.330E-12,      &
   &  9.480E-13,  7.000E-13,  5.370E-13,  4.170E-13,  3.330E-13,      &
   &  2.760E-13,  2.350E-13,  2.040E-13,  1.780E-13,  1.590E-13,      &
   &  1.420E-13,  1.230E-13,  1.060E-13,  8.950E-14,  7.450E-14,      &
   &  6.060E-14,  5.020E-14,  4.220E-14,  3.610E-14,  3.130E-14,      &
   &  2.750E-14,  2.450E-14,  2.190E-14,  1.970E-14,  1.790E-14,      &
   &  1.630E-14,  1.490E-14,  1.380E-14,  1.290E-14,  1.220E-14,      &
   &  1.170E-14,  1.140E-14,  1.140E-14,  1.180E-14,  1.260E-14/
   DATA F1801/                                                       &
   &  1.410E-14,  1.650E-14,  2.040E-14,  2.650E-14,  3.590E-14,      &
   &  4.870E-14,  6.410E-14,  8.470E-14,  1.070E-13,  1.320E-13,      &
   &  1.660E-13,  2.070E-13,  2.620E-13,  3.210E-13,  3.940E-13,      &
   &  4.840E-13,  5.810E-13,  6.840E-13,  7.730E-13,  8.670E-13,      &
   &  9.620E-13,  1.060E-12,  1.190E-12,  1.350E-12,  1.560E-12,      &
   &  1.850E-12,  2.210E-12,  2.660E-12,  3.140E-12,  3.600E-12,      &
   &  3.990E-12,  4.220E-12,  4.360E-12,  4.350E-12,  4.250E-12,      &
   &  4.100E-12,  3.990E-12,  3.950E-12,  3.870E-12,  3.830E-12,      &
   &  3.820E-12,  3.890E-12,  3.980E-12,  4.060E-12,  4.190E-12,      &
   &  4.220E-12,  4.170E-12,  3.900E-12,  3.500E-12,  3.030E-12/
   DATA F1851/                                                       &
   &  2.490E-12,  1.990E-12,  1.510E-12,  1.150E-12,  8.680E-13,      &
   &  6.530E-13,  4.890E-13,  3.780E-13,  2.920E-13,  2.250E-13,      &
   &  1.740E-13,  1.350E-13,  1.070E-13,  8.540E-14,  7.020E-14,      &
   &  5.950E-14,  5.130E-14,  4.480E-14,  3.890E-14,  3.340E-14,      &
   &  2.940E-14,  2.630E-14,  2.420E-14,  2.300E-14,  2.270E-14,      &
   &  2.370E-14,  2.560E-14,  2.780E-14,  3.160E-14,  3.800E-14,      &
   &  4.860E-14,  6.250E-14,  8.150E-14,  1.010E-13,  1.240E-13,      &
   &  1.490E-13,  1.730E-13,  1.990E-13,  2.220E-13,  2.440E-13,      &
   &  2.610E-13,  2.720E-13,  2.850E-13,  2.940E-13,  3.100E-13,      &
   &  3.290E-13,  3.500E-13,  3.620E-13,  3.660E-13,  3.700E-13/
   DATA F1901/                                                       &
   &  3.720E-13,  3.830E-13,  4.010E-13,  4.260E-13,  4.510E-13,      &
   &  4.640E-13,  4.490E-13,  4.050E-13,  3.450E-13,  2.790E-13,      &
   &  2.120E-13,  1.460E-13,  9.670E-14,  6.180E-14,  3.980E-14,      &
   &  2.630E-14,  1.810E-14,  1.320E-14,  1.030E-14,  8.420E-15,      &
   &  7.220E-15,  6.390E-15,  5.780E-15,  5.330E-15,  5.000E-15,      &
   &  4.760E-15,  4.610E-15,  4.540E-15,  4.550E-15,  4.640E-15,      &
   &  4.810E-15,  5.060E-15,  5.410E-15,  5.870E-15,  6.440E-15,      &
   &  7.150E-15,  8.020E-15,  9.100E-15,  1.050E-14,  1.220E-14,      &
   &  1.460E-14,  1.800E-14,  2.280E-14,  2.980E-14,  4.020E-14,      &
   &  5.530E-14,  7.810E-14,  1.040E-13,  1.310E-13,  1.580E-13/
   DATA F1951/                                                       &
   &  1.890E-13,  2.200E-13,  2.540E-13,  3.010E-13,  3.670E-13,      &
   &  4.430E-13,  5.420E-13,  6.850E-13,  8.440E-13,  1.040E-12,      &
   &  1.310E-12,  1.610E-12,  1.950E-12,  2.380E-12,  2.880E-12,      &
   &  3.420E-12,  4.030E-12,  4.650E-12,  5.240E-12,  5.700E-12,      &
   &  6.090E-12,  6.390E-12,  6.480E-12,  6.580E-12,  6.660E-12,      &
   &  6.700E-12,  6.660E-12,  6.650E-12,  6.800E-12,  7.090E-12,      &
   &  7.540E-12,  8.060E-12,  8.620E-12,  8.890E-12,  8.680E-12,      &
   &  8.040E-12,  7.000E-12,  5.830E-12,  4.580E-12,  3.490E-12,      &
   &  2.710E-12,  2.100E-12,  1.650E-12,  1.310E-12,  1.060E-12,      &
   &  8.350E-13,  6.450E-13,  5.020E-13,  3.940E-13,  3.090E-13/
   DATA F2001/                                                       &
   &  2.460E-13/
!
end block data BFH2O
!
!     --------------------------------------------------------------
!
SUBROUTINE FRNCO2 (V1C,V2C,DVC,NPTC,C,tave,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   dimension tdep_bandhead(1196:1220)

!     Temparature dependence coefficients for wavenumbers between 2386
!     and 2434. Computed based on (line-coupled) continuum coefficients
!     at 250K and 296K, set to unity at T_eff (determined by invariance
!     of calculations in this region for IASI low PWV cases).
   data (tdep_bandhead(i),i=1196,1220)/                              &
   &    1.44E-01,3.61E-01,5.71E-01,7.63E-01,8.95E-01,                 &
   &    9.33E-01,8.75E-01,7.30E-01,5.47E-01,3.79E-01,                 &
   &    2.55E-01,1.78E-01,1.34E-01,1.07E-01,9.06E-02,                 &
   &    7.83E-02,6.83E-02,6.00E-02,5.30E-02,4.72E-02,                 &
   &    4.24E-02,3.83E-02,3.50E-02,3.23E-02,3.01E-02/
   data t_eff/246./

   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   COMMON /FCO2/ V1S,V2S,DVS,NPTS,S(5003)
   DIMENSION C(*)
!
   trat = tave/t_eff

   DVC = DVS
   v1ss = v1s
   v2ss = v2s
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.GE.1).AND.(I.LE.NPTS)) THEN
         tcor = 1.
         if (i .ge. 1196 .and. i .le. 1220) tcor = (trat)**       &
            tdep_bandhead(i)
         C(J) = tcor * S(I)
      ENDIF
10 END DO
!
   RETURN
!
end subroutine FRNCO2
!
!     --------------------------------------------------------------
!
BLOCK DATA BFCO2
!
   IMPLICIT REAL*8           (V)
!
!               09/14/07     Hartmann line coupling: isotopes 1 and 2
!
!               UNITS OF (CM**3/MOL)*1.E-20
!
   COMMON /FCO2/ V1,V2,DV,NPT,                                       &
   &              F0000( 2),F0001(50),F0051(50),F0101(50),F0151(50),  &
   &              F0201(50),F0251(50),F0301(50),F0351(50),F0401(50),  &
   &              F0451(50),F0501(50),F0551(50),F0601(50),F0651(50),  &
   &              F0701(50),F0751(50),F0801(50),F0851(50),F0901(50),  &
   &              F0951(50),F1001(50),F1051(50),F1101(50),F1151(50),  &
   &              F1201(50),F1251(50),F1301(50),F1351(50),F1401(50),  &
   &              F1451(50),F1501(50),F1551(50),F1601(50),F1651(50),  &
   &              F1701(50),F1751(50),F1801(50),F1851(50),F1901(50),  &
   &              F1951(50),F2001(50),F2051(50),F2101(50),F2151(50),  &
   &              F2201(50),F2251(50),F2301(50),F2351(50),F2401(50),  &
   &              F2451(50),F2501(50),F2551(50),F2601(50),F2651(50),  &
   &              F2701(50),F2751(50),F2801(50),F2851(50),F2901(50),  &
   &              F2951(50),F3001(50),F3051(50),F3101(50),F3151(50),  &
   &              F3201(50),F3251(50),F3301(50),F3351(50),F3401(50),  &
   &              F3451(50),F3501(50),F3551(50),F3601(50),F3651(50),  &
   &              F3701(50),F3751(50),F3801(50),F3851(50),F3901(50),  &
   &              F3951(50),F4001(50),F4051(50),F4101(50),F4151(50),  &
   &              F4201(50),F4251(50),F4301(50),F4351(50),F4401(50),  &
   &              F4451(50),F4501(50),F4551(50),F4601(50),F4651(50),  &
   &              F4701(50),F4751(50),F4801(50),F4851(50),F4901(50),  &
   &              F4951(50), F5001(1)
!
   DATA V1,V2,DV,NPT / -4.0, 10000.0, 2.0, 5003/
!
   DATA F0000/                                                       &
   &  8.391e-13,  8.359e-13/
   DATA F0001/                                                       &
   &  8.345e-13,  8.359e-13,  8.391e-13,  8.439e-13,  8.500e-13,      &
   &  8.573e-13,  8.655e-13,  8.750e-13,  8.859e-13,  8.981e-13,      &
   &  9.119e-13,  9.273e-13,  9.442e-13,  9.626e-13,  9.825e-13,      &
   &  1.004e-12,  1.027e-12,  1.052e-12,  1.078e-12,  1.106e-12,      &
   &  1.135e-12,  1.166e-12,  1.199e-12,  1.233e-12,  1.270e-12,      &
   &  1.308e-12,  1.347e-12,  1.389e-12,  1.432e-12,  1.478e-12,      &
   &  1.525e-12,  1.574e-12,  1.625e-12,  1.678e-12,  1.734e-12,      &
   &  1.791e-12,  1.850e-12,  1.912e-12,  1.976e-12,  2.042e-12,      &
   &  2.110e-12,  2.181e-12,  2.254e-12,  2.329e-12,  2.408e-12,      &
   &  2.488e-12,  2.572e-12,  2.658e-12,  2.746e-12,  2.838e-12/
   DATA F0051/                                                       &
   &  2.933e-12,  3.030e-12,  3.130e-12,  3.234e-12,  3.341e-12,      &
   &  3.451e-12,  3.564e-12,  3.681e-12,  3.801e-12,  3.925e-12,      &
   &  4.053e-12,  4.184e-12,  4.320e-12,  4.459e-12,  4.603e-12,      &
   &  4.750e-12,  4.902e-12,  5.059e-12,  5.220e-12,  5.386e-12,      &
   &  5.556e-12,  5.732e-12,  5.912e-12,  6.098e-12,  6.290e-12,      &
   &  6.487e-12,  6.689e-12,  6.897e-12,  7.112e-12,  7.333e-12,      &
   &  7.560e-12,  7.793e-12,  8.034e-12,  8.281e-12,  8.535e-12,      &
   &  8.797e-12,  9.067e-12,  9.344e-12,  9.629e-12,  9.923e-12,      &
   &  1.023e-11,  1.054e-11,  1.086e-11,  1.119e-11,  1.152e-11,      &
   &  1.187e-11,  1.223e-11,  1.260e-11,  1.298e-11,  1.338e-11/
   DATA F0101/                                                       &
   &  1.378e-11,  1.419e-11,  1.462e-11,  1.506e-11,  1.552e-11,      &
   &  1.598e-11,  1.646e-11,  1.696e-11,  1.747e-11,  1.799e-11,      &
   &  1.854e-11,  1.909e-11,  1.967e-11,  2.026e-11,  2.087e-11,      &
   &  2.150e-11,  2.214e-11,  2.281e-11,  2.350e-11,  2.421e-11,      &
   &  2.494e-11,  2.569e-11,  2.647e-11,  2.727e-11,  2.810e-11,      &
   &  2.895e-11,  2.983e-11,  3.073e-11,  3.167e-11,  3.263e-11,      &
   &  3.363e-11,  3.466e-11,  3.572e-11,  3.681e-11,  3.794e-11,      &
   &  3.910e-11,  4.031e-11,  4.155e-11,  4.283e-11,  4.416e-11,      &
   &  4.552e-11,  4.694e-11,  4.840e-11,  4.991e-11,  5.147e-11,      &
   &  5.308e-11,  5.474e-11,  5.646e-11,  5.824e-11,  6.008e-11/
   DATA F0151/                                                       &
   &  6.199e-11,  6.395e-11,  6.599e-11,  6.809e-11,  7.027e-11,      &
   &  7.253e-11,  7.486e-11,  7.728e-11,  7.978e-11,  8.237e-11,      &
   &  8.505e-11,  8.782e-11,  9.070e-11,  9.368e-11,  9.677e-11,      &
   &  9.997e-11,  1.033e-10,  1.067e-10,  1.103e-10,  1.140e-10,      &
   &  1.178e-10,  1.218e-10,  1.259e-10,  1.302e-10,  1.347e-10,      &
   &  1.393e-10,  1.441e-10,  1.491e-10,  1.542e-10,  1.596e-10,      &
   &  1.652e-10,  1.710e-10,  1.770e-10,  1.833e-10,  1.899e-10,      &
   &  1.966e-10,  2.037e-10,  2.111e-10,  2.187e-10,  2.267e-10,      &
   &  2.350e-10,  2.437e-10,  2.527e-10,  2.621e-10,  2.719e-10,      &
   &  2.821e-10,  2.928e-10,  3.039e-10,  3.155e-10,  3.276e-10/
   DATA F0201/                                                       &
   &  3.403e-10,  3.535e-10,  3.674e-10,  3.818e-10,  3.969e-10,      &
   &  4.127e-10,  4.293e-10,  4.466e-10,  4.647e-10,  4.837e-10,      &
   &  5.036e-10,  5.245e-10,  5.464e-10,  5.695e-10,  5.937e-10,      &
   &  6.191e-10,  6.459e-10,  6.742e-10,  7.039e-10,  7.353e-10,      &
   &  7.683e-10,  8.032e-10,  8.400e-10,  8.785e-10,  9.213e-10,      &
   &  9.661e-10,  1.013e-09,  1.062e-09,  1.114e-09,  1.169e-09,      &
   &  1.226e-09,  1.286e-09,  1.350e-09,  1.417e-09,  1.488e-09,      &
   &  1.564e-09,  1.643e-09,  1.728e-09,  1.818e-09,  1.914e-09,      &
   &  2.016e-09,  2.126e-09,  2.244e-09,  2.378e-09,  2.522e-09,      &
   &  2.677e-09,  2.843e-09,  3.024e-09,  3.220e-09,  3.437e-09/
   DATA F0251/                                                       &
   &  3.676e-09,  3.941e-09,  4.236e-09,  4.565e-09,  4.931e-09,      &
   &  5.336e-09,  5.784e-09,  6.276e-09,  6.817e-09,  7.410e-09,      &
   &  8.108e-09,  9.238e-09,  1.044e-08,  1.167e-08,  1.292e-08,      &
   &  1.419e-08,  1.548e-08,  1.679e-08,  1.816e-08,  1.956e-08,      &
   &  2.103e-08,  2.258e-08,  2.425e-08,  2.609e-08,  2.815e-08,      &
   &  3.051e-08,  3.326e-08,  3.653e-08,  4.040e-08,  4.546e-08,      &
   &  5.160e-08,  5.909e-08,  6.836e-08,  7.995e-08,  9.457e-08,      &
   &  1.130e-07,  1.366e-07,  1.731e-07,  2.156e-07,  2.655e-07,      &
   &  3.241e-07,  3.927e-07,  4.723e-07,  5.639e-07,  6.685e-07,      &
   &  7.872e-07,  9.510e-07,  1.199e-06,  1.528e-06,  1.847e-06/
   DATA F0301/                                                       &
   &  2.155e-06,  2.453e-06,  2.742e-06,  3.026e-06,  3.308e-06,      &
   &  3.594e-06,  3.895e-06,  4.227e-06,  4.609e-06,  5.073e-06,      &
   &  5.660e-06,  6.406e-06,  7.410e-06,  8.801e-06,  1.052e-05,      &
   &  1.261e-05,  1.512e-05,  1.808e-05,  2.149e-05,  2.535e-05,      &
   &  2.970e-05,  3.413e-05,  4.573e-05,  6.122e-05,  7.642e-05,      &
   &  9.099e-05,  1.048e-04,  1.176e-04,  1.294e-04,  1.398e-04,      &
   &  1.488e-04,  1.561e-04,  1.616e-04,  1.653e-04,  1.670e-04,      &
   &  1.666e-04,  1.642e-04,  1.599e-04,  1.538e-04,  1.460e-04,      &
   &  1.365e-04,  1.256e-04,  1.132e-04,  9.969e-05,  8.499e-05,      &
   &  6.925e-05,  5.260e-05,  4.198e-05,  3.475e-05,  2.933e-05/
   DATA F0351/                                                       &
   &  2.457e-05,  2.043e-05,  1.690e-05,  1.394e-05,  1.151e-05,      &
   &  9.559e-06,  8.021e-06,  6.826e-06,  5.920e-06,  5.232e-06,      &
   &  4.703e-06,  4.285e-06,  3.941e-06,  3.643e-06,  3.368e-06,      &
   &  3.106e-06,  2.847e-06,  2.585e-06,  2.316e-06,  2.040e-06,      &
   &  1.755e-06,  1.463e-06,  1.165e-06,  9.093e-07,  8.070e-07,      &
   &  6.991e-07,  6.009e-07,  5.118e-07,  4.319e-07,  3.612e-07,      &
   &  2.996e-07,  2.465e-07,  2.012e-07,  1.633e-07,  1.370e-07,      &
   &  1.172e-07,  1.015e-07,  8.892e-08,  7.876e-08,  7.053e-08,      &
   &  6.378e-08,  5.817e-08,  5.367e-08,  4.988e-08,  4.660e-08,      &
   &  4.372e-08,  4.113e-08,  3.874e-08,  3.651e-08,  3.439e-08/
   DATA F0401/                                                       &
   &  3.234e-08,  3.033e-08,  2.836e-08,  2.641e-08,  2.448e-08,      &
   &  2.255e-08,  2.062e-08,  1.869e-08,  1.677e-08,  1.532e-08,      &
   &  1.417e-08,  1.317e-08,  1.226e-08,  1.143e-08,  1.068e-08,      &
   &  1.001e-08,  9.396e-09,  8.846e-09,  8.350e-09,  7.902e-09,      &
   &  7.496e-09,  7.125e-09,  6.784e-09,  6.467e-09,  6.171e-09,      &
   &  5.892e-09,  5.629e-09,  5.379e-09,  5.161e-09,  4.952e-09,      &
   &  4.756e-09,  4.571e-09,  4.396e-09,  4.232e-09,  4.077e-09,      &
   &  3.930e-09,  3.792e-09,  3.662e-09,  3.541e-09,  3.426e-09,      &
   &  3.320e-09,  3.221e-09,  3.129e-09,  3.045e-09,  2.967e-09,      &
   &  2.898e-09,  2.839e-09,  2.788e-09,  2.743e-09,  2.706e-09/
   DATA F0451/                                                       &
   &  2.676e-09,  2.653e-09,  2.638e-09,  2.633e-09,  2.640e-09,      &
   &  2.661e-09,  2.700e-09,  2.760e-09,  2.844e-09,  2.953e-09,      &
   &  3.092e-09,  3.260e-09,  3.459e-09,  3.685e-09,  3.936e-09,      &
   &  4.205e-09,  4.482e-09,  4.756e-09,  5.015e-09,  5.255e-09,      &
   &  5.479e-09,  5.696e-09,  5.915e-09,  6.139e-09,  6.369e-09,      &
   &  6.596e-09,  6.811e-09,  6.999e-09,  7.148e-09,  7.244e-09,      &
   &  7.280e-09,  7.252e-09,  7.159e-09,  7.007e-09,  6.799e-09,      &
   &  6.547e-09,  6.260e-09,  5.948e-09,  5.618e-09,  5.279e-09,      &
   &  4.935e-09,  4.584e-09,  4.224e-09,  3.849e-09,  3.457e-09,      &
   &  3.055e-09,  2.658e-09,  2.285e-09,  1.948e-09,  1.661e-09/
   DATA F0501/                                                       &
   &  1.429e-09,  1.254e-09,  1.133e-09,  1.060e-09,  1.031e-09,      &
   &  1.038e-09,  1.079e-09,  1.151e-09,  1.254e-09,  1.391e-09,      &
   &  1.564e-09,  1.777e-09,  2.033e-09,  2.335e-09,  2.682e-09,      &
   &  3.072e-09,  3.497e-09,  3.947e-09,  4.407e-09,  4.859e-09,      &
   &  5.293e-09,  5.707e-09,  6.108e-09,  6.510e-09,  6.919e-09,      &
   &  7.337e-09,  7.762e-09,  8.178e-09,  8.573e-09,  8.923e-09,      &
   &  9.211e-09,  9.419e-09,  9.536e-09,  9.554e-09,  9.473e-09,      &
   &  9.300e-09,  9.043e-09,  8.715e-09,  8.330e-09,  7.904e-09,      &
   &  7.447e-09,  6.968e-09,  6.471e-09,  5.955e-09,  5.414e-09,      &
   &  4.843e-09,  4.246e-09,  3.644e-09,  3.054e-09,  2.504e-09/
   DATA F0551/                                                       &
   &  2.010e-09,  1.584e-09,  1.234e-09,  9.578e-10,  7.490e-10,      &
   &  5.962e-10,  4.887e-10,  4.149e-10,  3.652e-10,  3.309e-10,      &
   &  3.067e-10,  2.892e-10,  2.758e-10,  2.650e-10,  2.560e-10,      &
   &  2.482e-10,  2.413e-10,  2.351e-10,  2.295e-10,  2.244e-10,      &
   &  2.197e-10,  2.154e-10,  2.113e-10,  2.075e-10,  2.039e-10,      &
   &  2.005e-10,  1.973e-10,  1.942e-10,  1.912e-10,  1.884e-10,      &
   &  1.857e-10,  1.831e-10,  1.805e-10,  1.781e-10,  1.757e-10,      &
   &  1.734e-10,  1.712e-10,  1.691e-10,  1.670e-10,  1.649e-10,      &
   &  1.629e-10,  1.610e-10,  1.591e-10,  1.573e-10,  1.555e-10,      &
   &  1.537e-10,  1.520e-10,  1.503e-10,  1.487e-10,  1.471e-10/
   DATA F0601/                                                       &
   &  1.455e-10,  1.440e-10,  1.425e-10,  1.410e-10,  1.396e-10,      &
   &  1.381e-10,  1.368e-10,  1.354e-10,  1.341e-10,  1.328e-10,      &
   &  1.315e-10,  1.303e-10,  1.290e-10,  1.278e-10,  1.267e-10,      &
   &  1.255e-10,  1.244e-10,  1.233e-10,  1.222e-10,  1.211e-10,      &
   &  1.201e-10,  1.190e-10,  1.180e-10,  1.170e-10,  1.161e-10,      &
   &  1.151e-10,  1.142e-10,  1.133e-10,  1.124e-10,  1.115e-10,      &
   &  1.106e-10,  1.098e-10,  1.089e-10,  1.081e-10,  1.073e-10,      &
   &  1.065e-10,  1.058e-10,  1.050e-10,  1.043e-10,  1.036e-10,      &
   &  1.028e-10,  1.022e-10,  1.015e-10,  1.008e-10,  1.001e-10,      &
   &  9.950e-11,  9.887e-11,  9.825e-11,  9.765e-11,  9.705e-11/
   DATA F0651/                                                       &
   &  9.647e-11,  9.591e-11,  9.535e-11,  9.480e-11,  9.427e-11,      &
   &  9.375e-11,  9.324e-11,  9.274e-11,  9.225e-11,  9.177e-11,      &
   &  9.131e-11,  9.085e-11,  9.040e-11,  8.997e-11,  8.954e-11,      &
   &  8.913e-11,  8.872e-11,  8.833e-11,  8.794e-11,  8.757e-11,      &
   &  8.720e-11,  8.684e-11,  8.650e-11,  8.616e-11,  8.583e-11,      &
   &  8.551e-11,  8.520e-11,  8.490e-11,  8.461e-11,  8.432e-11,      &
   &  8.405e-11,  8.378e-11,  8.352e-11,  8.327e-11,  8.303e-11,      &
   &  8.280e-11,  8.257e-11,  8.236e-11,  8.215e-11,  8.195e-11,      &
   &  8.176e-11,  8.158e-11,  8.140e-11,  8.123e-11,  8.107e-11,      &
   &  8.092e-11,  8.078e-11,  8.064e-11,  8.051e-11,  8.039e-11/
   DATA F0701/                                                       &
   &  8.028e-11,  8.018e-11,  8.008e-11,  7.999e-11,  7.991e-11,      &
   &  7.983e-11,  7.976e-11,  7.971e-11,  7.965e-11,  7.961e-11,      &
   &  7.957e-11,  7.954e-11,  7.952e-11,  7.951e-11,  7.950e-11,      &
   &  7.950e-11,  7.951e-11,  7.952e-11,  7.955e-11,  7.958e-11,      &
   &  7.961e-11,  7.966e-11,  7.971e-11,  7.977e-11,  7.984e-11,      &
   &  7.991e-11,  8.000e-11,  8.009e-11,  8.018e-11,  8.029e-11,      &
   &  8.040e-11,  8.052e-11,  8.065e-11,  8.079e-11,  8.093e-11,      &
   &  8.108e-11,  8.124e-11,  8.141e-11,  8.158e-11,  8.176e-11,      &
   &  8.195e-11,  8.215e-11,  8.235e-11,  8.257e-11,  8.279e-11,      &
   &  8.302e-11,  8.326e-11,  8.350e-11,  8.376e-11,  8.402e-11/
   DATA F0751/                                                       &
   &  8.429e-11,  8.457e-11,  8.486e-11,  8.515e-11,  8.546e-11,      &
   &  8.577e-11,  8.609e-11,  8.642e-11,  8.676e-11,  8.711e-11,      &
   &  8.747e-11,  8.783e-11,  8.821e-11,  8.859e-11,  8.898e-11,      &
   &  8.939e-11,  8.980e-11,  9.022e-11,  9.065e-11,  9.109e-11,      &
   &  9.155e-11,  9.201e-11,  9.248e-11,  9.296e-11,  9.345e-11,      &
   &  9.395e-11,  9.446e-11,  9.499e-11,  9.552e-11,  9.606e-11,      &
   &  9.662e-11,  9.719e-11,  9.776e-11,  9.835e-11,  9.895e-11,      &
   &  9.956e-11,  1.002e-10,  1.008e-10,  1.015e-10,  1.021e-10,      &
   &  1.028e-10,  1.035e-10,  1.042e-10,  1.049e-10,  1.056e-10,      &
   &  1.063e-10,  1.071e-10,  1.079e-10,  1.086e-10,  1.094e-10/
   DATA F0801/                                                       &
   &  1.102e-10,  1.111e-10,  1.119e-10,  1.127e-10,  1.136e-10,      &
   &  1.145e-10,  1.154e-10,  1.163e-10,  1.172e-10,  1.181e-10,      &
   &  1.191e-10,  1.201e-10,  1.211e-10,  1.221e-10,  1.231e-10,      &
   &  1.241e-10,  1.252e-10,  1.263e-10,  1.274e-10,  1.285e-10,      &
   &  1.296e-10,  1.308e-10,  1.319e-10,  1.331e-10,  1.344e-10,      &
   &  1.356e-10,  1.368e-10,  1.381e-10,  1.394e-10,  1.407e-10,      &
   &  1.421e-10,  1.435e-10,  1.448e-10,  1.463e-10,  1.477e-10,      &
   &  1.492e-10,  1.506e-10,  1.522e-10,  1.537e-10,  1.553e-10,      &
   &  1.568e-10,  1.585e-10,  1.601e-10,  1.618e-10,  1.635e-10,      &
   &  1.652e-10,  1.670e-10,  1.688e-10,  1.706e-10,  1.725e-10/
   DATA F0851/                                                       &
   &  1.744e-10,  1.763e-10,  1.782e-10,  1.802e-10,  1.823e-10,      &
   &  1.843e-10,  1.864e-10,  1.886e-10,  1.907e-10,  1.929e-10,      &
   &  1.952e-10,  1.975e-10,  1.998e-10,  2.022e-10,  2.046e-10,      &
   &  2.071e-10,  2.096e-10,  2.122e-10,  2.148e-10,  2.174e-10,      &
   &  2.202e-10,  2.229e-10,  2.257e-10,  2.286e-10,  2.315e-10,      &
   &  2.345e-10,  2.375e-10,  2.406e-10,  2.438e-10,  2.470e-10,      &
   &  2.503e-10,  2.536e-10,  2.570e-10,  2.605e-10,  2.641e-10,      &
   &  2.677e-10,  2.714e-10,  2.752e-10,  2.791e-10,  2.831e-10,      &
   &  2.871e-10,  2.913e-10,  2.955e-10,  2.999e-10,  3.043e-10,      &
   &  3.089e-10,  3.136e-10,  3.184e-10,  3.234e-10,  3.285e-10/
   DATA F0901/                                                       &
   &  3.337e-10,  3.391e-10,  3.447e-10,  3.505e-10,  3.564e-10,      &
   &  3.626e-10,  3.689e-10,  3.754e-10,  3.822e-10,  3.892e-10,      &
   &  3.965e-10,  4.041e-10,  4.120e-10,  4.203e-10,  4.290e-10,      &
   &  4.381e-10,  4.478e-10,  4.581e-10,  4.691e-10,  4.807e-10,      &
   &  4.931e-10,  5.062e-10,  5.200e-10,  5.346e-10,  5.499e-10,      &
   &  5.659e-10,  5.826e-10,  6.002e-10,  6.189e-10,  6.392e-10,      &
   &  6.606e-10,  6.834e-10,  7.080e-10,  7.348e-10,  7.643e-10,      &
   &  7.974e-10,  8.351e-10,  8.791e-10,  9.301e-10,  9.894e-10,      &
   &  1.059e-09,  1.140e-09,  1.235e-09,  1.346e-09,  1.472e-09,      &
   &  1.616e-09,  1.777e-09,  1.954e-09,  2.147e-09,  2.347e-09/
   DATA F0951/                                                       &
   &  2.549e-09,  2.747e-09,  2.941e-09,  3.126e-09,  3.309e-09,      &
   &  3.589e-09,  3.848e-09,  4.066e-09,  4.238e-09,  4.365e-09,      &
   &  4.445e-09,  4.481e-09,  4.474e-09,  4.426e-09,  4.341e-09,      &
   &  4.223e-09,  4.076e-09,  3.905e-09,  3.717e-09,  3.516e-09,      &
   &  3.308e-09,  3.100e-09,  2.897e-09,  2.701e-09,  2.514e-09,      &
   &  2.332e-09,  2.156e-09,  1.985e-09,  1.819e-09,  1.659e-09,      &
   &  1.635e-09,  1.630e-09,  1.634e-09,  1.642e-09,  1.654e-09,      &
   &  1.668e-09,  1.685e-09,  1.705e-09,  1.728e-09,  1.754e-09,      &
   &  1.785e-09,  1.819e-09,  1.859e-09,  1.903e-09,  1.952e-09,      &
   &  2.007e-09,  2.068e-09,  2.134e-09,  2.207e-09,  2.287e-09/
   DATA F1001/                                                       &
   &  2.374e-09,  2.469e-09,  2.573e-09,  2.686e-09,  2.810e-09,      &
   &  2.947e-09,  3.098e-09,  3.285e-09,  3.494e-09,  3.729e-09,      &
   &  3.999e-09,  4.311e-09,  4.677e-09,  5.110e-09,  5.628e-09,      &
   &  6.239e-09,  6.963e-09,  7.806e-09,  8.772e-09,  9.860e-09,      &
   &  1.106e-08,  1.236e-08,  1.375e-08,  1.519e-08,  1.666e-08,      &
   &  1.817e-08,  1.962e-08,  2.240e-08,  2.516e-08,  2.762e-08,      &
   &  2.977e-08,  3.157e-08,  3.303e-08,  3.411e-08,  3.479e-08,      &
   &  3.524e-08,  3.533e-08,  3.506e-08,  3.445e-08,  3.353e-08,      &
   &  3.236e-08,  3.096e-08,  2.939e-08,  2.769e-08,  2.589e-08,      &
   &  2.404e-08,  2.215e-08,  2.024e-08,  1.832e-08,  1.640e-08/
   DATA F1051/                                                       &
   &  1.449e-08,  1.261e-08,  1.203e-08,  1.188e-08,  1.185e-08,      &
   &  1.188e-08,  1.195e-08,  1.205e-08,  1.216e-08,  1.228e-08,      &
   &  1.250e-08,  1.278e-08,  1.307e-08,  1.337e-08,  1.369e-08,      &
   &  1.402e-08,  1.436e-08,  1.474e-08,  1.513e-08,  1.556e-08,      &
   &  1.601e-08,  1.651e-08,  1.703e-08,  1.760e-08,  1.821e-08,      &
   &  1.886e-08,  1.957e-08,  2.032e-08,  2.119e-08,  2.213e-08,      &
   &  2.313e-08,  2.420e-08,  2.534e-08,  2.656e-08,  2.785e-08,      &
   &  2.923e-08,  3.070e-08,  3.228e-08,  3.396e-08,  3.577e-08,      &
   &  3.772e-08,  3.981e-08,  4.207e-08,  4.451e-08,  4.715e-08,      &
   &  5.001e-08,  5.313e-08,  5.654e-08,  6.027e-08,  6.438e-08/
   DATA F1101/                                                       &
   &  6.891e-08,  7.395e-08,  7.959e-08,  8.594e-08,  9.315e-08,      &
   &  1.014e-07,  1.109e-07,  1.220e-07,  1.351e-07,  1.506e-07,      &
   &  1.691e-07,  1.914e-07,  2.184e-07,  2.511e-07,  2.910e-07,      &
   &  3.396e-07,  3.986e-07,  4.703e-07,  5.565e-07,  6.598e-07,      &
   &  7.824e-07,  9.260e-07,  1.092e-06,  1.282e-06,  1.495e-06,      &
   &  1.730e-06,  1.984e-06,  2.254e-06,  2.533e-06,  2.816e-06,      &
   &  3.100e-06,  3.386e-06,  3.683e-06,  3.999e-06,  4.343e-06,      &
   &  4.721e-06,  5.140e-06,  5.603e-06,  6.119e-06,  6.699e-06,      &
   &  7.361e-06,  8.132e-06,  9.049e-06,  1.017e-05,  1.155e-05,      &
   &  1.329e-05,  1.549e-05,  1.827e-05,  2.178e-05,  2.619e-05/
   DATA F1151/                                                       &
   &  3.169e-05,  3.845e-05,  4.668e-05,  5.656e-05,  6.822e-05,      &
   &  8.175e-05,  9.719e-05,  1.145e-04,  1.334e-04,  1.537e-04,      &
   &  1.748e-04,  1.960e-04,  2.167e-04,  2.363e-04,  2.549e-04,      &
   &  2.731e-04,  2.911e-04,  3.094e-04,  3.277e-04,  3.458e-04,      &
   &  3.629e-04,  3.782e-04,  3.910e-04,  4.003e-04,  4.058e-04,      &
   &  4.069e-04,  4.038e-04,  3.964e-04,  3.853e-04,  3.709e-04,      &
   &  3.538e-04,  3.348e-04,  3.143e-04,  2.928e-04,  2.707e-04,      &
   &  2.481e-04,  2.250e-04,  2.010e-04,  1.760e-04,  1.504e-04,      &
   &  1.250e-04,  1.010e-04,  7.907e-05,  6.005e-05,  4.430e-05,      &
   &  3.180e-05,  2.231e-05,  1.544e-05,  1.069e-05,  7.481e-06/
   DATA F1201/                                                       &
   &  5.397e-06,  4.055e-06,  3.178e-06,  2.581e-06,  2.153e-06,      &
   &  1.827e-06,  1.568e-06,  1.357e-06,  1.182e-06,  1.035e-06,      &
   &  9.114e-07,  8.059e-07,  7.156e-07,  6.378e-07,  5.706e-07,      &
   &  5.122e-07,  4.612e-07,  4.165e-07,  3.771e-07,  3.424e-07,      &
   &  3.116e-07,  2.843e-07,  2.599e-07,  2.381e-07,  2.186e-07,      &
   &  2.010e-07,  1.852e-07,  1.709e-07,  1.580e-07,  1.463e-07,      &
   &  1.356e-07,  1.259e-07,  1.171e-07,  1.090e-07,  1.017e-07,      &
   &  9.489e-08,  8.868e-08,  8.297e-08,  7.772e-08,  7.288e-08,      &
   &  6.841e-08,  6.428e-08,  6.045e-08,  5.691e-08,  5.362e-08,      &
   &  5.057e-08,  4.773e-08,  4.508e-08,  4.262e-08,  4.032e-08/
   DATA F1251/                                                       &
   &  3.818e-08,  3.617e-08,  3.430e-08,  3.255e-08,  3.090e-08,      &
   &  2.936e-08,  2.791e-08,  2.655e-08,  2.527e-08,  2.407e-08,      &
   &  2.294e-08,  2.187e-08,  2.087e-08,  1.992e-08,  1.902e-08,      &
   &  1.818e-08,  1.738e-08,  1.662e-08,  1.590e-08,  1.522e-08,      &
   &  1.458e-08,  1.397e-08,  1.339e-08,  1.284e-08,  1.232e-08,      &
   &  1.183e-08,  1.136e-08,  1.091e-08,  1.048e-08,  1.008e-08,      &
   &  9.691e-09,  9.322e-09,  8.971e-09,  8.636e-09,  8.316e-09,      &
   &  8.011e-09,  7.720e-09,  7.441e-09,  7.175e-09,  6.921e-09,      &
   &  6.677e-09,  6.444e-09,  6.221e-09,  6.008e-09,  5.803e-09,      &
   &  5.607e-09,  5.419e-09,  5.239e-09,  5.066e-09,  4.900e-09/
   DATA F1301/                                                       &
   &  4.740e-09,  4.587e-09,  4.440e-09,  4.299e-09,  4.163e-09,      &
   &  4.033e-09,  3.907e-09,  3.787e-09,  3.671e-09,  3.559e-09,      &
   &  3.451e-09,  3.347e-09,  3.247e-09,  3.151e-09,  3.058e-09,      &
   &  2.969e-09,  2.883e-09,  2.799e-09,  2.719e-09,  2.642e-09,      &
   &  2.567e-09,  2.495e-09,  2.425e-09,  2.357e-09,  2.292e-09,      &
   &  2.229e-09,  2.169e-09,  2.110e-09,  2.053e-09,  1.998e-09,      &
   &  1.945e-09,  1.893e-09,  1.843e-09,  1.795e-09,  1.748e-09,      &
   &  1.703e-09,  1.659e-09,  1.617e-09,  1.576e-09,  1.536e-09,      &
   &  1.497e-09,  1.460e-09,  1.424e-09,  1.388e-09,  1.354e-09,      &
   &  1.321e-09,  1.289e-09,  1.258e-09,  1.227e-09,  1.198e-09/
   DATA F1351/                                                       &
   &  1.169e-09,  1.142e-09,  1.115e-09,  1.088e-09,  1.063e-09,      &
   &  1.038e-09,  1.014e-09,  9.908e-10,  9.681e-10,  9.460e-10,      &
   &  9.246e-10,  9.037e-10,  8.834e-10,  8.636e-10,  8.444e-10,      &
   &  8.257e-10,  8.074e-10,  7.897e-10,  7.724e-10,  7.556e-10,      &
   &  7.393e-10,  7.233e-10,  7.078e-10,  6.927e-10,  6.779e-10,      &
   &  6.636e-10,  6.496e-10,  6.359e-10,  6.226e-10,  6.097e-10,      &
   &  5.970e-10,  5.847e-10,  5.727e-10,  5.610e-10,  5.495e-10,      &
   &  5.384e-10,  5.275e-10,  5.169e-10,  5.065e-10,  4.964e-10,      &
   &  4.865e-10,  4.769e-10,  4.675e-10,  4.583e-10,  4.493e-10,      &
   &  4.406e-10,  4.320e-10,  4.236e-10,  4.155e-10,  4.075e-10/
   DATA F1401/                                                       &
   &  3.997e-10,  3.921e-10,  3.846e-10,  3.774e-10,  3.703e-10,      &
   &  3.633e-10,  3.565e-10,  3.499e-10,  3.434e-10,  3.370e-10,      &
   &  3.308e-10,  3.248e-10,  3.188e-10,  3.130e-10,  3.073e-10,      &
   &  3.018e-10,  2.963e-10,  2.910e-10,  2.858e-10,  2.807e-10,      &
   &  2.757e-10,  2.708e-10,  2.660e-10,  2.613e-10,  2.567e-10,      &
   &  2.523e-10,  2.479e-10,  2.436e-10,  2.393e-10,  2.352e-10,      &
   &  2.312e-10,  2.272e-10,  2.233e-10,  2.195e-10,  2.158e-10,      &
   &  2.121e-10,  2.086e-10,  2.051e-10,  2.016e-10,  1.983e-10,      &
   &  1.950e-10,  1.917e-10,  1.886e-10,  1.855e-10,  1.824e-10,      &
   &  1.794e-10,  1.765e-10,  1.736e-10,  1.708e-10,  1.680e-10/
   DATA F1451/                                                       &
   &  1.653e-10,  1.627e-10,  1.601e-10,  1.575e-10,  1.550e-10,      &
   &  1.525e-10,  1.501e-10,  1.478e-10,  1.454e-10,  1.431e-10,      &
   &  1.409e-10,  1.387e-10,  1.366e-10,  1.344e-10,  1.324e-10,      &
   &  1.303e-10,  1.283e-10,  1.264e-10,  1.244e-10,  1.225e-10,      &
   &  1.207e-10,  1.189e-10,  1.171e-10,  1.153e-10,  1.136e-10,      &
   &  1.119e-10,  1.102e-10,  1.086e-10,  1.070e-10,  1.054e-10,      &
   &  1.039e-10,  1.023e-10,  1.008e-10,  9.937e-11,  9.793e-11,      &
   &  9.651e-11,  9.512e-11,  9.375e-11,  9.241e-11,  9.109e-11,      &
   &  8.979e-11,  8.852e-11,  8.726e-11,  8.603e-11,  8.483e-11,      &
   &  8.364e-11,  8.247e-11,  8.132e-11,  8.020e-11,  7.909e-11/
   DATA F1501/                                                       &
   &  7.800e-11,  7.693e-11,  7.588e-11,  7.485e-11,  7.384e-11,      &
   &  7.284e-11,  7.186e-11,  7.090e-11,  6.995e-11,  6.902e-11,      &
   &  6.811e-11,  6.721e-11,  6.633e-11,  6.546e-11,  6.461e-11,      &
   &  6.377e-11,  6.295e-11,  6.214e-11,  6.135e-11,  6.057e-11,      &
   &  5.980e-11,  5.905e-11,  5.831e-11,  5.758e-11,  5.687e-11,      &
   &  5.617e-11,  5.548e-11,  5.481e-11,  5.415e-11,  5.350e-11,      &
   &  5.286e-11,  5.224e-11,  5.163e-11,  5.103e-11,  5.045e-11,      &
   &  4.989e-11,  4.934e-11,  4.881e-11,  4.830e-11,  4.781e-11,      &
   &  4.734e-11,  4.689e-11,  4.646e-11,  4.605e-11,  4.565e-11,      &
   &  4.527e-11,  4.489e-11,  4.453e-11,  4.417e-11,  4.381e-11/
   DATA F1551/                                                       &
   &  4.346e-11,  4.312e-11,  4.279e-11,  4.247e-11,  4.216e-11,      &
   &  4.187e-11,  4.160e-11,  4.136e-11,  4.114e-11,  4.097e-11,      &
   &  4.084e-11,  4.075e-11,  4.073e-11,  4.077e-11,  4.089e-11,      &
   &  4.110e-11,  4.142e-11,  4.185e-11,  4.238e-11,  4.301e-11,      &
   &  4.372e-11,  4.448e-11,  4.527e-11,  4.604e-11,  4.677e-11,      &
   &  4.742e-11,  4.798e-11,  4.844e-11,  4.879e-11,  4.905e-11,      &
   &  4.924e-11,  4.936e-11,  4.940e-11,  4.939e-11,  4.937e-11,      &
   &  4.935e-11,  4.938e-11,  4.946e-11,  4.963e-11,  4.987e-11,      &
   &  5.018e-11,  5.057e-11,  5.100e-11,  5.147e-11,  5.195e-11,      &
   &  5.242e-11,  5.285e-11,  5.321e-11,  5.346e-11,  5.358e-11/
   DATA F1601/                                                       &
   &  5.352e-11,  5.327e-11,  5.283e-11,  5.217e-11,  5.133e-11,      &
   &  5.034e-11,  4.918e-11,  4.786e-11,  4.643e-11,  4.492e-11,      &
   &  4.338e-11,  4.186e-11,  4.042e-11,  3.908e-11,  3.789e-11,      &
   &  3.685e-11,  3.598e-11,  3.527e-11,  3.472e-11,  3.430e-11,      &
   &  3.402e-11,  3.385e-11,  3.379e-11,  3.382e-11,  3.395e-11,      &
   &  3.420e-11,  3.466e-11,  3.519e-11,  3.579e-11,  3.645e-11,      &
   &  3.718e-11,  3.798e-11,  3.886e-11,  3.982e-11,  4.088e-11,      &
   &  4.204e-11,  4.331e-11,  4.471e-11,  4.626e-11,  4.800e-11,      &
   &  4.996e-11,  5.219e-11,  5.475e-11,  5.769e-11,  6.108e-11,      &
   &  6.494e-11,  6.930e-11,  7.414e-11,  7.941e-11,  8.501e-11/
   DATA F1651/                                                       &
   &  9.086e-11,  9.693e-11,  1.029e-10,  1.087e-10,  1.142e-10,      &
   &  1.193e-10,  1.240e-10,  1.282e-10,  1.342e-10,  1.412e-10,      &
   &  1.485e-10,  1.559e-10,  1.639e-10,  1.725e-10,  1.821e-10,      &
   &  1.928e-10,  2.047e-10,  2.175e-10,  2.312e-10,  2.454e-10,      &
   &  2.596e-10,  2.735e-10,  2.867e-10,  2.988e-10,  3.094e-10,      &
   &  3.182e-10,  3.248e-10,  3.291e-10,  3.308e-10,  3.298e-10,      &
   &  3.260e-10,  3.195e-10,  3.101e-10,  3.000e-10,  2.891e-10,      &
   &  2.769e-10,  2.635e-10,  2.493e-10,  2.347e-10,  2.204e-10,      &
   &  2.068e-10,  1.944e-10,  1.836e-10,  1.746e-10,  1.675e-10,      &
   &  1.622e-10,  1.586e-10,  1.565e-10,  1.559e-10,  1.564e-10/
   DATA F1701/                                                       &
   &  1.579e-10,  1.604e-10,  1.636e-10,  1.675e-10,  1.721e-10,      &
   &  1.774e-10,  1.833e-10,  1.899e-10,  1.973e-10,  2.055e-10,      &
   &  2.146e-10,  2.247e-10,  2.360e-10,  2.486e-10,  2.626e-10,      &
   &  2.784e-10,  2.960e-10,  3.158e-10,  3.380e-10,  3.630e-10,      &
   &  3.913e-10,  4.234e-10,  4.599e-10,  5.020e-10,  5.507e-10,      &
   &  6.075e-10,  6.737e-10,  7.509e-10,  8.409e-10,  9.452e-10,      &
   &  1.066e-09,  1.203e-09,  1.360e-09,  1.535e-09,  1.730e-09,      &
   &  1.945e-09,  2.179e-09,  2.436e-09,  2.721e-09,  3.044e-09,      &
   &  3.417e-09,  3.853e-09,  4.362e-09,  4.955e-09,  5.639e-09,      &
   &  6.417e-09,  7.287e-09,  8.241e-09,  9.269e-09,  1.035e-08/
   DATA F1751/                                                       &
   &  1.145e-08,  1.256e-08,  1.365e-08,  1.475e-08,  1.589e-08,      &
   &  1.710e-08,  1.842e-08,  1.988e-08,  2.149e-08,  2.326e-08,      &
   &  2.520e-08,  2.731e-08,  2.962e-08,  3.215e-08,  3.498e-08,      &
   &  3.823e-08,  4.199e-08,  4.646e-08,  5.180e-08,  5.819e-08,      &
   &  6.582e-08,  7.484e-08,  8.539e-08,  9.756e-08,  1.114e-07,      &
   &  1.270e-07,  1.443e-07,  1.635e-07,  1.848e-07,  2.090e-07,      &
   &  2.371e-07,  2.706e-07,  3.111e-07,  3.601e-07,  4.190e-07,      &
   &  4.892e-07,  5.717e-07,  6.668e-07,  7.748e-07,  8.947e-07,      &
   &  1.025e-06,  1.162e-06,  1.301e-06,  1.439e-06,  1.570e-06,      &
   &  1.692e-06,  1.808e-06,  1.920e-06,  2.034e-06,  2.150e-06/
   DATA F1801/                                                       &
   &  2.269e-06,  2.386e-06,  2.499e-06,  2.600e-06,  2.684e-06,      &
   &  2.746e-06,  2.781e-06,  2.788e-06,  2.768e-06,  2.720e-06,      &
   &  2.648e-06,  2.556e-06,  2.448e-06,  2.328e-06,  2.198e-06,      &
   &  2.063e-06,  1.923e-06,  1.779e-06,  1.628e-06,  1.469e-06,      &
   &  1.300e-06,  1.126e-06,  9.526e-07,  7.885e-07,  6.394e-07,      &
   &  5.112e-07,  4.064e-07,  3.262e-07,  2.698e-07,  2.352e-07,      &
   &  2.206e-07,  2.235e-07,  2.429e-07,  2.784e-07,  3.302e-07,      &
   &  3.995e-07,  4.880e-07,  5.968e-07,  7.277e-07,  8.811e-07,      &
   &  1.056e-06,  1.251e-06,  1.461e-06,  1.680e-06,  1.901e-06,      &
   &  2.116e-06,  2.321e-06,  2.519e-06,  2.715e-06,  2.913e-06/
   DATA F1851/                                                       &
   &  3.115e-06,  3.321e-06,  3.526e-06,  3.724e-06,  3.905e-06,      &
   &  4.059e-06,  4.181e-06,  4.261e-06,  4.295e-06,  4.283e-06,      &
   &  4.225e-06,  4.126e-06,  3.991e-06,  3.826e-06,  3.637e-06,      &
   &  3.432e-06,  3.215e-06,  2.989e-06,  2.755e-06,  2.511e-06,      &
   &  2.254e-06,  1.985e-06,  1.706e-06,  1.429e-06,  1.166e-06,      &
   &  9.250e-07,  7.135e-07,  5.356e-07,  3.921e-07,  2.814e-07,      &
   &  1.988e-07,  1.398e-07,  9.915e-08,  7.137e-08,  5.280e-08,      &
   &  4.038e-08,  3.192e-08,  2.601e-08,  2.164e-08,  1.836e-08,      &
   &  1.578e-08,  1.371e-08,  1.199e-08,  1.056e-08,  9.353e-09,      &
   &  8.321e-09,  7.436e-09,  6.672e-09,  6.010e-09,  5.433e-09/
   DATA F1901/                                                       &
   &  4.929e-09,  4.487e-09,  4.098e-09,  3.754e-09,  3.448e-09,      &
   &  3.174e-09,  2.929e-09,  2.708e-09,  2.507e-09,  2.325e-09,      &
   &  2.159e-09,  2.006e-09,  1.866e-09,  1.738e-09,  1.619e-09,      &
   &  1.510e-09,  1.409e-09,  1.316e-09,  1.230e-09,  1.149e-09,      &
   &  1.074e-09,  1.003e-09,  9.372e-10,  8.757e-10,  8.186e-10,      &
   &  7.661e-10,  7.179e-10,  6.740e-10,  6.339e-10,  5.974e-10,      &
   &  5.641e-10,  5.335e-10,  5.053e-10,  4.791e-10,  4.548e-10,      &
   &  4.321e-10,  4.109e-10,  3.909e-10,  3.722e-10,  3.546e-10,      &
   &  3.380e-10,  3.224e-10,  3.076e-10,  2.936e-10,  2.804e-10,      &
   &  2.680e-10,  2.562e-10,  2.451e-10,  2.346e-10,  2.246e-10/
   DATA F1951/                                                       &
   &  2.153e-10,  2.064e-10,  1.981e-10,  1.901e-10,  1.826e-10,      &
   &  1.755e-10,  1.687e-10,  1.623e-10,  1.562e-10,  1.503e-10,      &
   &  1.448e-10,  1.395e-10,  1.344e-10,  1.296e-10,  1.250e-10,      &
   &  1.206e-10,  1.165e-10,  1.125e-10,  1.087e-10,  1.050e-10,      &
   &  1.015e-10,  9.822e-11,  9.506e-11,  9.203e-11,  8.915e-11,      &
   &  8.640e-11,  8.377e-11,  8.127e-11,  7.888e-11,  7.660e-11,      &
   &  7.442e-11,  7.234e-11,  7.035e-11,  6.844e-11,  6.661e-11,      &
   &  6.486e-11,  6.316e-11,  6.154e-11,  5.998e-11,  5.850e-11,      &
   &  5.710e-11,  5.582e-11,  5.460e-11,  5.344e-11,  5.234e-11,      &
   &  5.132e-11,  5.040e-11,  4.961e-11,  4.896e-11,  4.845e-11/
   DATA F2001/                                                       &
   &  4.806e-11,  4.774e-11,  4.742e-11,  4.702e-11,  4.652e-11,      &
   &  4.594e-11,  4.527e-11,  4.451e-11,  4.366e-11,  4.274e-11,      &
   &  4.172e-11,  4.063e-11,  3.946e-11,  3.822e-11,  3.693e-11,      &
   &  3.560e-11,  3.428e-11,  3.293e-11,  3.155e-11,  3.014e-11,      &
   &  2.872e-11,  2.734e-11,  2.603e-11,  2.482e-11,  2.376e-11,      &
   &  2.286e-11,  2.211e-11,  2.146e-11,  2.088e-11,  2.032e-11,      &
   &  1.980e-11,  1.930e-11,  1.882e-11,  1.836e-11,  1.792e-11,      &
   &  1.749e-11,  1.708e-11,  1.669e-11,  1.630e-11,  1.593e-11,      &
   &  1.557e-11,  1.522e-11,  1.488e-11,  1.456e-11,  1.424e-11,      &
   &  1.392e-11,  1.362e-11,  1.333e-11,  1.304e-11,  1.276e-11/
   DATA F2051/                                                       &
   &  1.249e-11,  1.223e-11,  1.197e-11,  1.172e-11,  1.147e-11,      &
   &  1.124e-11,  1.100e-11,  1.078e-11,  1.056e-11,  1.034e-11,      &
   &  1.013e-11,  9.924e-12,  9.723e-12,  9.528e-12,  9.336e-12,      &
   &  9.149e-12,  8.967e-12,  8.789e-12,  8.615e-12,  8.445e-12,      &
   &  8.278e-12,  8.116e-12,  7.957e-12,  7.802e-12,  7.650e-12,      &
   &  7.502e-12,  7.357e-12,  7.215e-12,  7.076e-12,  6.940e-12,      &
   &  6.807e-12,  6.677e-12,  6.550e-12,  6.426e-12,  6.304e-12,      &
   &  6.185e-12,  6.068e-12,  5.954e-12,  5.842e-12,  5.732e-12,      &
   &  5.625e-12,  5.520e-12,  5.417e-12,  5.316e-12,  5.217e-12,      &
   &  5.121e-12,  5.026e-12,  4.933e-12,  4.842e-12,  4.753e-12/
   DATA F2101/                                                       &
   &  4.665e-12,  4.580e-12,  4.496e-12,  4.414e-12,  4.333e-12,      &
   &  4.254e-12,  4.176e-12,  4.100e-12,  4.025e-12,  3.952e-12,      &
   &  3.880e-12,  3.810e-12,  3.741e-12,  3.673e-12,  3.607e-12,      &
   &  3.542e-12,  3.478e-12,  3.415e-12,  3.353e-12,  3.293e-12,      &
   &  3.233e-12,  3.175e-12,  3.118e-12,  3.062e-12,  3.007e-12,      &
   &  2.953e-12,  2.899e-12,  2.847e-12,  2.796e-12,  2.746e-12,      &
   &  2.696e-12,  2.648e-12,  2.600e-12,  2.554e-12,  2.508e-12,      &
   &  2.463e-12,  2.418e-12,  2.375e-12,  2.332e-12,  2.290e-12,      &
   &  2.249e-12,  2.208e-12,  2.168e-12,  2.129e-12,  2.090e-12,      &
   &  2.053e-12,  2.015e-12,  1.979e-12,  1.943e-12,  1.908e-12/
   DATA F2151/                                                       &
   &  1.873e-12,  1.839e-12,  1.805e-12,  1.772e-12,  1.740e-12,      &
   &  1.708e-12,  1.677e-12,  1.646e-12,  1.616e-12,  1.586e-12,      &
   &  1.557e-12,  1.528e-12,  1.500e-12,  1.472e-12,  1.445e-12,      &
   &  1.418e-12,  1.391e-12,  1.365e-12,  1.340e-12,  1.314e-12,      &
   &  1.290e-12,  1.265e-12,  1.241e-12,  1.218e-12,  1.195e-12,      &
   &  1.172e-12,  1.150e-12,  1.128e-12,  1.106e-12,  1.085e-12,      &
   &  1.064e-12,  1.043e-12,  1.023e-12,  1.003e-12,  9.835e-13,      &
   &  9.643e-13,  9.453e-13,  9.267e-13,  9.084e-13,  8.904e-13,      &
   &  8.727e-13,  8.553e-13,  8.382e-13,  8.213e-13,  8.048e-13,      &
   &  7.885e-13,  7.725e-13,  7.568e-13,  7.414e-13,  7.262e-13/
   DATA F2201/                                                       &
   &  7.113e-13,  6.966e-13,  6.822e-13,  6.680e-13,  6.541e-13,      &
   &  6.404e-13,  6.270e-13,  6.138e-13,  6.009e-13,  5.882e-13,      &
   &  5.757e-13,  5.634e-13,  5.514e-13,  5.396e-13,  5.280e-13,      &
   &  5.167e-13,  5.055e-13,  4.946e-13,  4.839e-13,  4.734e-13,      &
   &  4.631e-13,  4.531e-13,  4.432e-13,  4.336e-13,  4.241e-13,      &
   &  4.149e-13,  4.059e-13,  3.971e-13,  3.885e-13,  3.801e-13,      &
   &  3.720e-13,  3.640e-13,  3.562e-13,  3.487e-13,  3.414e-13,      &
   &  3.342e-13,  3.274e-13,  3.207e-13,  3.142e-13,  3.080e-13,      &
   &  3.020e-13,  2.963e-13,  2.908e-13,  2.856e-13,  2.806e-13,      &
   &  2.759e-13,  2.715e-13,  2.673e-13,  2.635e-13,  2.600e-13/
   DATA F2251/                                                       &
   &  2.569e-13,  2.542e-13,  2.519e-13,  2.500e-13,  2.487e-13,      &
   &  2.479e-13,  2.477e-13,  2.483e-13,  2.498e-13,  2.524e-13,      &
   &  2.562e-13,  2.614e-13,  2.684e-13,  2.776e-13,  2.894e-13,      &
   &  3.045e-13,  3.235e-13,  3.472e-13,  3.765e-13,  4.123e-13,      &
   &  4.554e-13,  5.064e-13,  5.654e-13,  6.321e-13,  7.055e-13,      &
   &  7.844e-13,  8.666e-13,  9.498e-13,  1.032e-12,  1.111e-12,      &
   &  1.186e-12,  1.257e-12,  1.322e-12,  1.377e-12,  1.473e-12,      &
   &  1.580e-12,  1.690e-12,  1.806e-12,  1.930e-12,  2.066e-12,      &
   &  2.217e-12,  2.386e-12,  2.572e-12,  2.771e-12,  2.983e-12,      &
   &  3.199e-12,  3.414e-12,  3.619e-12,  3.810e-12,  3.980e-12/
   DATA F2301/                                                       &
   &  4.123e-12,  4.235e-12,  4.311e-12,  4.348e-12,  4.343e-12,      &
   &  4.296e-12,  4.204e-12,  4.069e-12,  3.891e-12,  3.711e-12,      &
   &  3.520e-12,  3.310e-12,  3.085e-12,  2.854e-12,  2.622e-12,      &
   &  2.400e-12,  2.198e-12,  2.022e-12,  1.879e-12,  1.772e-12,      &
   &  1.703e-12,  1.673e-12,  1.680e-12,  1.723e-12,  1.799e-12,      &
   &  1.911e-12,  2.055e-12,  2.231e-12,  2.440e-12,  2.683e-12,      &
   &  2.960e-12,  3.271e-12,  3.620e-12,  4.012e-12,  4.456e-12,      &
   &  4.958e-12,  5.520e-12,  6.142e-12,  6.818e-12,  7.538e-12,      &
   &  8.289e-12,  9.052e-12,  9.821e-12,  1.060e-11,  1.140e-11,      &
   &  1.226e-11,  1.320e-11,  1.426e-11,  1.547e-11,  1.686e-11/
   DATA F2351/                                                       &
   &  1.846e-11,  2.030e-11,  2.243e-11,  2.487e-11,  2.767e-11,      &
   &  3.089e-11,  3.452e-11,  3.858e-11,  4.305e-11,  4.791e-11,      &
   &  5.309e-11,  5.851e-11,  6.412e-11,  6.995e-11,  7.609e-11,      &
   &  8.271e-11,  8.987e-11,  9.758e-11,  1.058e-10,  1.146e-10,      &
   &  1.237e-10,  1.331e-10,  1.428e-10,  1.530e-10,  1.638e-10,      &
   &  1.755e-10,  1.886e-10,  2.036e-10,  2.212e-10,  2.418e-10,      &
   &  2.662e-10,  2.949e-10,  3.287e-10,  3.677e-10,  4.130e-10,      &
   &  4.648e-10,  5.232e-10,  5.877e-10,  6.574e-10,  7.307e-10,      &
   &  8.058e-10,  8.802e-10,  9.540e-10,  1.027e-09,  1.102e-09,      &
   &  1.181e-09,  1.268e-09,  1.365e-09,  1.476e-09,  1.601e-09/
   DATA F2401/                                                       &
   &  1.744e-09,  1.908e-09,  2.099e-09,  2.323e-09,  2.589e-09,      &
   &  2.904e-09,  3.278e-09,  3.721e-09,  4.237e-09,  4.828e-09,      &
   &  5.491e-09,  6.214e-09,  6.977e-09,  7.752e-09,  8.508e-09,      &
   &  9.225e-09,  9.899e-09,  1.055e-08,  1.120e-08,  1.187e-08,      &
   &  1.256e-08,  1.325e-08,  1.395e-08,  1.461e-08,  1.520e-08,      &
   &  1.569e-08,  1.605e-08,  1.626e-08,  1.631e-08,  1.620e-08,      &
   &  1.594e-08,  1.553e-08,  1.501e-08,  1.439e-08,  1.370e-08,      &
   &  1.296e-08,  1.219e-08,  1.139e-08,  1.055e-08,  9.668e-09,      &
   &  8.719e-09,  7.720e-09,  6.699e-09,  5.699e-09,  4.763e-09,      &
   &  3.920e-09,  3.193e-09,  2.592e-09,  2.116e-09,  1.754e-09/
   DATA F2451/                                                       &
   &  1.496e-09,  1.322e-09,  1.222e-09,  1.179e-09,  1.188e-09,      &
   &  1.242e-09,  1.341e-09,  1.488e-09,  1.691e-09,  1.963e-09,      &
   &  2.320e-09,  2.783e-09,  3.375e-09,  4.124e-09,  5.059e-09,      &
   &  6.212e-09,  7.615e-09,  9.297e-09,  1.128e-08,  1.359e-08,      &
   &  1.623e-08,  1.920e-08,  2.246e-08,  2.598e-08,  2.967e-08,      &
   &  3.341e-08,  3.709e-08,  4.060e-08,  4.393e-08,  4.715e-08,      &
   &  5.035e-08,  5.359e-08,  5.689e-08,  6.020e-08,  6.340e-08,      &
   &  6.637e-08,  6.892e-08,  7.091e-08,  7.222e-08,  7.275e-08,      &
   &  7.249e-08,  7.146e-08,  6.970e-08,  6.732e-08,  6.443e-08,      &
   &  6.112e-08,  5.751e-08,  5.370e-08,  4.974e-08,  4.570e-08/
   DATA F2501/                                                       &
   &  4.155e-08,  3.726e-08,  3.279e-08,  2.817e-08,  2.355e-08,      &
   &  1.910e-08,  1.501e-08,  1.141e-08,  8.412e-09,  6.023e-09,      &
   &  4.215e-09,  2.923e-09,  2.033e-09,  1.455e-09,  1.093e-09,      &
   &  8.740e-10,  7.431e-10,  6.659e-10,  6.264e-10,  6.148e-10,      &
   &  6.290e-10,  6.704e-10,  7.425e-10,  8.507e-10,  1.003e-09,      &
   &  1.207e-09,  1.473e-09,  1.814e-09,  2.240e-09,  2.764e-09,      &
   &  3.393e-09,  4.135e-09,  4.986e-09,  5.942e-09,  6.983e-09,      &
   &  8.082e-09,  9.200e-09,  1.030e-08,  1.135e-08,  1.236e-08,      &
   &  1.334e-08,  1.434e-08,  1.536e-08,  1.642e-08,  1.749e-08,      &
   &  1.854e-08,  1.954e-08,  2.042e-08,  2.114e-08,  2.167e-08/
   DATA F2551/                                                       &
   &  2.198e-08,  2.205e-08,  2.189e-08,  2.153e-08,  2.099e-08,      &
   &  2.030e-08,  1.949e-08,  1.859e-08,  1.763e-08,  1.662e-08,      &
   &  1.556e-08,  1.445e-08,  1.325e-08,  1.197e-08,  1.060e-08,      &
   &  9.205e-09,  7.832e-09,  6.544e-09,  5.383e-09,  4.378e-09,      &
   &  3.539e-09,  2.863e-09,  2.326e-09,  1.907e-09,  1.578e-09,      &
   &  1.311e-09,  1.088e-09,  8.975e-10,  7.318e-10,  5.895e-10,      &
   &  4.686e-10,  3.686e-10,  2.875e-10,  2.242e-10,  1.750e-10,      &
   &  1.380e-10,  1.104e-10,  8.948e-11,  7.378e-11,  6.180e-11,      &
   &  5.251e-11,  4.520e-11,  3.940e-11,  3.475e-11,  3.099e-11,      &
   &  2.792e-11,  2.537e-11,  2.323e-11,  2.141e-11,  1.985e-11/
   DATA F2601/                                                       &
   &  1.850e-11,  1.734e-11,  1.633e-11,  1.545e-11,  1.468e-11,      &
   &  1.399e-11,  1.338e-11,  1.281e-11,  1.228e-11,  1.179e-11,      &
   &  1.133e-11,  1.091e-11,  1.055e-11,  1.019e-11,  9.838e-12,      &
   &  9.494e-12,  9.161e-12,  8.845e-12,  8.544e-12,  8.258e-12,      &
   &  7.982e-12,  7.708e-12,  7.428e-12,  7.141e-12,  6.853e-12,      &
   &  6.574e-12,  6.319e-12,  6.101e-12,  5.933e-12,  5.825e-12,      &
   &  5.790e-12,  5.845e-12,  6.018e-12,  6.355e-12,  6.862e-12,      &
   &  7.521e-12,  8.372e-12,  9.455e-12,  1.076e-11,  1.231e-11,      &
   &  1.415e-11,  1.635e-11,  1.900e-11,  2.233e-11,  2.688e-11,      &
   &  3.335e-11,  4.316e-11,  5.277e-11,  6.203e-11,  7.088e-11/
   DATA F2651/                                                       &
   &  7.925e-11,  8.706e-11,  9.420e-11,  1.005e-10,  1.059e-10,      &
   &  1.101e-10,  1.130e-10,  1.145e-10,  1.146e-10,  1.132e-10,      &
   &  1.104e-10,  1.063e-10,  1.008e-10,  9.418e-11,  8.642e-11,      &
   &  7.766e-11,  6.805e-11,  5.779e-11,  4.724e-11,  3.706e-11,      &
   &  2.875e-11,  2.362e-11,  1.918e-11,  1.520e-11,  1.168e-11,      &
   &  8.671e-12,  6.195e-12,  4.263e-12,  2.844e-12,  1.873e-12,      &
   &  1.259e-12,  8.916e-13,  6.808e-13,  5.494e-13,  4.535e-13,      &
   &  3.758e-13,  3.113e-13,  2.571e-13,  2.110e-13,  1.714e-13,      &
   &  1.372e-13,  1.074e-13,  8.140e-14,  5.888e-14,  3.992e-14,      &
   &  2.574e-14,  2.098e-14,  1.985e-14,  1.894e-14,  1.804e-14/
   DATA F2701/                                                       &
   &  1.721e-14,  1.645e-14,  1.570e-14,  1.509e-14,  1.449e-14,      &
   &  1.389e-14,  1.336e-14,  1.291e-14,  1.238e-14,  1.200e-14,      &
   &  1.155e-14,  1.117e-14,  1.079e-14,  1.042e-14,  1.004e-14,      &
   &  9.736e-15,  9.434e-15,  9.132e-15,  8.830e-15,  8.604e-15,      &
   &  8.302e-15,  8.075e-15,  7.849e-15,  7.623e-15,  7.381e-15,      &
   &  7.170e-15,  6.974e-15,  6.777e-15,  6.589e-15,  6.408e-15,      &
   &  6.234e-15,  6.068e-15,  5.909e-15,  5.751e-15,  5.600e-15,      &
   &  5.457e-15,  5.321e-15,  5.185e-15,  5.057e-15,  4.928e-15,      &
   &  4.808e-15,  4.694e-15,  4.581e-15,  4.475e-15,  4.370e-15,      &
   &  4.272e-15,  4.174e-15,  4.083e-15,  3.992e-15,  3.909e-15/
   DATA F2751/                                                       &
   &  3.834e-15,  3.758e-15,  3.691e-15,  3.623e-15,  3.562e-15,      &
   &  3.509e-15,  3.464e-15,  3.426e-15,  3.396e-15,  3.374e-15,      &
   &  3.351e-15,  3.343e-15,  3.336e-15,  3.343e-15,  3.358e-15,      &
   &  3.374e-15,  3.396e-15,  3.434e-15,  3.472e-15,  3.517e-15,      &
   &  3.562e-15,  3.608e-15,  3.653e-15,  3.706e-15,  3.751e-15,      &
   &  3.437e-14,  7.989e-14,  1.282e-13,  1.776e-13,  2.263e-13,      &
   &  2.729e-13,  3.172e-13,  3.605e-13,  4.054e-13,  4.547e-13,      &
   &  5.102e-13,  5.724e-13,  6.399e-13,  7.079e-13,  7.639e-13,      &
   &  8.028e-13,  8.255e-13,  8.329e-13,  8.264e-13,  8.071e-13,      &
   &  7.764e-13,  7.358e-13,  6.868e-13,  6.310e-13,  5.700e-13/
   DATA F2801/                                                       &
   &  5.051e-13,  4.376e-13,  3.681e-13,  2.968e-13,  2.234e-13,      &
   &  1.473e-13,  6.888e-14,  3.200e-15,  3.147e-15,  3.087e-15,      &
   &  3.011e-15,  2.898e-15,  2.755e-15,  2.611e-15,  2.483e-15,      &
   &  2.392e-15,  2.340e-15,  2.317e-15,  2.317e-15,  2.332e-15,      &
   &  2.347e-15,  2.385e-15,  2.430e-15,  2.491e-15,  2.558e-15,      &
   &  2.619e-15,  3.956e-14,  8.813e-14,  1.397e-13,  1.926e-13,      &
   &  2.457e-13,  2.970e-13,  3.453e-13,  3.918e-13,  4.390e-13,      &
   &  4.899e-13,  5.472e-13,  6.123e-13,  6.845e-13,  7.615e-13,      &
   &  8.299e-13,  8.795e-13,  9.111e-13,  9.260e-13,  9.252e-13,      &
   &  9.099e-13,  8.817e-13,  8.421e-13,  7.925e-13,  7.350e-13/
   DATA F2851/                                                       &
   &  6.711e-13,  6.022e-13,  5.297e-13,  4.544e-13,  3.767e-13,      &
   &  2.965e-13,  2.132e-13,  1.265e-13,  3.808e-14,  2.408e-15,      &
   &  2.302e-15,  2.158e-15,  1.985e-15,  1.796e-15,  1.615e-15,      &
   &  1.449e-15,  1.313e-15,  1.177e-15,  1.057e-15,  9.509e-16,      &
   &  8.604e-16,  7.849e-16,  7.336e-16,  6.928e-16,  6.589e-16,      &
   &  6.257e-16,  5.947e-16,  5.698e-16,  5.494e-16,  5.321e-16,      &
   &  5.155e-16,  5.004e-16,  4.875e-16,  4.770e-16,  4.679e-16,      &
   &  4.604e-16,  4.528e-16,  4.475e-16,  4.430e-16,  4.408e-16,      &
   &  4.385e-16,  4.377e-16,  4.385e-16,  4.400e-16,  4.438e-16,      &
   &  4.498e-16,  4.574e-16,  4.664e-16,  4.777e-16,  4.936e-16/
   DATA F2901/                                                       &
   &  5.155e-16,  5.442e-16,  5.781e-16,  6.166e-16,  6.581e-16,      &
   &  7.011e-16,  7.479e-16,  8.000e-16,  8.528e-16,  9.057e-16,      &
   &  9.509e-16,  9.887e-16,  1.019e-15,  1.049e-15,  1.079e-15,      &
   &  1.117e-15,  1.170e-15,  1.230e-15,  1.291e-15,  1.351e-15,      &
   &  1.404e-15,  1.457e-15,  1.509e-15,  1.555e-15,  1.592e-15,      &
   &  1.623e-15,  1.645e-15,  1.660e-15,  1.660e-15,  1.660e-15,      &
   &  1.660e-15,  1.660e-15,  1.653e-15,  1.638e-15,  1.623e-15,      &
   &  1.608e-15,  1.585e-15,  1.555e-15,  1.532e-15,  1.502e-15,      &
   &  1.487e-15,  1.487e-15,  1.502e-15,  1.509e-15,  1.525e-15,      &
   &  1.540e-15,  1.547e-15,  1.547e-15,  1.555e-15,  1.570e-15/
   DATA F2951/                                                       &
   &  1.592e-15,  1.638e-15,  1.691e-15,  1.751e-15,  1.826e-15,      &
   &  1.917e-15,  2.015e-15,  2.136e-15,  2.264e-15,  2.408e-15,      &
   &  2.543e-15,  2.687e-15,  2.830e-15,  2.981e-15,  3.140e-15,      &
   &  3.321e-15,  3.517e-15,  1.344e-14,  5.071e-14,  9.298e-14,      &
   &  1.386e-13,  1.864e-13,  2.353e-13,  2.838e-13,  3.317e-13,      &
   &  3.798e-13,  4.276e-13,  4.770e-13,  5.299e-13,  5.879e-13,      &
   &  6.531e-13,  7.277e-13,  8.135e-13,  9.128e-13,  1.028e-12,      &
   &  1.163e-12,  1.320e-12,  1.506e-12,  1.723e-12,  1.983e-12,      &
   &  2.292e-12,  2.657e-12,  3.080e-12,  3.560e-12,  4.093e-12,      &
   &  4.664e-12,  5.252e-12,  5.837e-12,  6.395e-12,  6.924e-12/
   DATA F3001/                                                       &
   &  7.424e-12,  7.909e-12,  8.394e-12,  8.890e-12,  9.398e-12,      &
   &  9.911e-12,  1.042e-11,  1.092e-11,  1.140e-11,  1.184e-11,      &
   &  1.226e-11,  1.268e-11,  1.311e-11,  1.360e-11,  1.420e-11,      &
   &  1.499e-11,  1.603e-11,  1.740e-11,  1.916e-11,  2.136e-11,      &
   &  2.401e-11,  2.709e-11,  3.051e-11,  3.419e-11,  3.796e-11,      &
   &  4.166e-11,  4.519e-11,  4.858e-11,  5.189e-11,  5.528e-11,      &
   &  5.883e-11,  6.257e-11,  6.644e-11,  7.032e-11,  7.407e-11,      &
   &  7.751e-11,  8.044e-11,  8.273e-11,  8.424e-11,  8.491e-11,      &
   &  8.472e-11,  8.371e-11,  8.193e-11,  7.954e-11,  7.660e-11,      &
   &  7.327e-11,  6.964e-11,  6.579e-11,  6.175e-11,  5.748e-11/
   DATA F3051/                                                       &
   &  5.295e-11,  4.807e-11,  4.289e-11,  3.757e-11,  3.230e-11,      &
   &  2.731e-11,  2.278e-11,  1.881e-11,  1.546e-11,  1.273e-11,      &
   &  1.058e-11,  8.955e-12,  7.746e-12,  6.902e-12,  6.331e-12,      &
   &  5.989e-12,  5.820e-12,  5.808e-12,  5.927e-12,  6.175e-12,      &
   &  6.552e-12,  7.059e-12,  7.708e-12,  8.512e-12,  9.506e-12,      &
   &  1.072e-11,  1.220e-11,  1.402e-11,  1.623e-11,  1.889e-11,      &
   &  2.204e-11,  2.573e-11,  2.996e-11,  3.477e-11,  4.016e-11,      &
   &  4.617e-11,  5.289e-11,  6.051e-11,  6.933e-11,  7.972e-11,      &
   &  9.209e-11,  1.068e-10,  1.243e-10,  1.447e-10,  1.684e-10,      &
   &  1.953e-10,  2.254e-10,  2.582e-10,  2.933e-10,  3.297e-10/
   DATA F3101/                                                       &
   &  3.662e-10,  4.014e-10,  4.345e-10,  4.655e-10,  4.951e-10,      &
   &  5.248e-10,  5.554e-10,  5.872e-10,  6.200e-10,  6.526e-10,      &
   &  6.838e-10,  7.115e-10,  7.342e-10,  7.503e-10,  7.590e-10,      &
   &  7.597e-10,  7.526e-10,  7.380e-10,  7.170e-10,  6.905e-10,      &
   &  6.597e-10,  6.255e-10,  5.889e-10,  5.506e-10,  5.107e-10,      &
   &  4.691e-10,  4.252e-10,  3.786e-10,  3.296e-10,  2.797e-10,      &
   &  2.307e-10,  1.848e-10,  1.440e-10,  1.091e-10,  8.078e-11,      &
   &  5.881e-11,  4.268e-11,  3.141e-11,  2.402e-11,  1.953e-11,      &
   &  1.711e-11,  1.615e-11,  1.625e-11,  1.722e-11,  1.901e-11,      &
   &  2.170e-11,  2.541e-11,  3.033e-11,  3.670e-11,  4.480e-11/
   DATA F3151/                                                       &
   &  5.496e-11,  6.753e-11,  8.287e-11,  1.013e-10,  1.232e-10,      &
   &  1.486e-10,  1.776e-10,  2.100e-10,  2.453e-10,  2.826e-10,      &
   &  3.207e-10,  3.585e-10,  3.947e-10,  4.295e-10,  4.634e-10,      &
   &  4.975e-10,  5.327e-10,  5.691e-10,  6.066e-10,  6.440e-10,      &
   &  6.801e-10,  7.130e-10,  7.410e-10,  7.624e-10,  7.761e-10,      &
   &  7.814e-10,  7.782e-10,  7.670e-10,  7.483e-10,  7.233e-10,      &
   &  6.931e-10,  6.588e-10,  6.215e-10,  5.819e-10,  5.402e-10,      &
   &  4.966e-10,  4.506e-10,  4.017e-10,  3.503e-10,  2.981e-10,      &
   &  2.468e-10,  1.987e-10,  1.552e-10,  1.179e-10,  8.705e-11,      &
   &  6.278e-11,  4.457e-11,  3.136e-11,  2.219e-11,  1.601e-11/
   DATA F3201/                                                       &
   &  1.191e-11,  9.181e-12,  7.312e-12,  6.003e-12,  5.040e-12,      &
   &  4.304e-12,  3.722e-12,  3.250e-12,  2.861e-12,  2.533e-12,      &
   &  2.254e-12,  2.016e-12,  1.812e-12,  1.637e-12,  1.491e-12,      &
   &  1.370e-12,  1.275e-12,  1.205e-12,  1.164e-12,  1.153e-12,      &
   &  1.179e-12,  1.248e-12,  1.370e-12,  1.561e-12,  1.839e-12,      &
   &  2.232e-12,  2.761e-12,  3.467e-12,  4.391e-12,  5.580e-12,      &
   &  7.078e-12,  8.932e-12,  1.117e-11,  1.382e-11,  1.688e-11,      &
   &  2.031e-11,  2.403e-11,  2.792e-11,  3.188e-11,  3.573e-11,      &
   &  3.941e-11,  4.294e-11,  4.645e-11,  5.000e-11,  5.367e-11,      &
   &  5.743e-11,  6.123e-11,  6.491e-11,  6.837e-11,  7.141e-11/
   DATA F3251/                                                       &
   &  7.392e-11,  7.576e-11,  7.687e-11,  7.720e-11,  7.678e-11,      &
   &  7.567e-11,  7.396e-11,  7.187e-11,  6.938e-11,  6.665e-11,      &
   &  6.373e-11,  6.070e-11,  5.757e-11,  5.426e-11,  5.070e-11,      &
   &  4.683e-11,  4.269e-11,  3.841e-11,  3.420e-11,  3.020e-11,      &
   &  2.654e-11,  2.330e-11,  2.050e-11,  1.812e-11,  1.611e-11,      &
   &  1.439e-11,  1.291e-11,  1.158e-11,  1.037e-11,  9.206e-12,      &
   &  8.084e-12,  7.010e-12,  6.045e-12,  5.105e-12,  4.203e-12,      &
   &  3.363e-12,  2.610e-12,  1.964e-12,  1.429e-12,  1.009e-12,      &
   &  6.869e-13,  4.567e-13,  2.909e-13,  1.790e-13,  1.033e-13,      &
   &  5.183e-14,  1.809e-14,  8.906e-15,  8.151e-15,  7.472e-15/
   DATA F3301/                                                       &
   &  6.860e-15,  6.332e-15,  5.902e-15,  5.555e-15,  5.268e-15,      &
   &  5.019e-15,  4.808e-15,  4.657e-15,  4.551e-15,  4.483e-15,      &
   &  4.453e-15,  4.460e-15,  4.498e-15,  4.581e-15,  4.702e-15,      &
   &  4.860e-15,  5.049e-15,  5.275e-15,  5.540e-15,  5.864e-15,      &
   &  6.272e-15,  3.285e-14,  8.652e-14,  1.564e-13,  2.476e-13,      &
   &  3.706e-13,  5.518e-13,  8.262e-13,  1.265e-12,  1.708e-12,      &
   &  2.148e-12,  2.584e-12,  3.012e-12,  3.432e-12,  3.836e-12,      &
   &  4.219e-12,  4.571e-12,  4.879e-12,  5.133e-12,  5.322e-12,      &
   &  5.439e-12,  5.482e-12,  5.452e-12,  5.352e-12,  5.185e-12,      &
   &  4.955e-12,  4.666e-12,  4.323e-12,  3.935e-12,  3.511e-12/
   DATA F3351/                                                       &
   &  3.069e-12,  2.641e-12,  2.302e-12,  2.118e-12,  1.978e-12,      &
   &  1.874e-12,  1.809e-12,  1.787e-12,  1.812e-12,  1.888e-12,      &
   &  2.015e-12,  2.195e-12,  2.429e-12,  2.721e-12,  3.074e-12,      &
   &  3.495e-12,  3.990e-12,  4.564e-12,  5.216e-12,  5.943e-12,      &
   &  6.726e-12,  7.566e-12,  8.462e-12,  9.412e-12,  1.040e-11,      &
   &  1.142e-11,  1.244e-11,  1.342e-11,  1.433e-11,  1.519e-11,      &
   &  1.602e-11,  1.688e-11,  1.782e-11,  1.890e-11,  2.010e-11,      &
   &  2.139e-11,  2.268e-11,  2.375e-11,  2.445e-11,  2.481e-11,      &
   &  2.486e-11,  2.464e-11,  2.417e-11,  2.349e-11,  2.262e-11,      &
   &  2.157e-11,  2.038e-11,  1.908e-11,  1.770e-11,  1.626e-11/
   DATA F3401/                                                       &
   &  1.477e-11,  1.323e-11,  1.165e-11,  9.999e-12,  8.290e-12,      &
   &  6.576e-12,  4.935e-12,  3.475e-12,  2.280e-12,  1.419e-12,      &
   &  9.163e-13,  7.470e-13,  7.556e-13,  7.834e-13,  8.307e-13,      &
   &  8.969e-13,  9.822e-13,  1.088e-12,  1.215e-12,  1.367e-12,      &
   &  1.546e-12,  1.755e-12,  1.999e-12,  2.283e-12,  2.614e-12,      &
   &  2.998e-12,  3.446e-12,  3.968e-12,  4.576e-12,  5.280e-12,      &
   &  6.095e-12,  7.031e-12,  8.103e-12,  9.320e-12,  1.070e-11,      &
   &  1.225e-11,  1.400e-11,  1.600e-11,  1.825e-11,  2.082e-11,      &
   &  2.377e-11,  2.719e-11,  3.115e-11,  3.576e-11,  4.110e-11,      &
   &  4.725e-11,  5.425e-11,  6.219e-11,  7.118e-11,  8.133e-11/
   DATA F3451/                                                       &
   &  9.264e-11,  1.052e-10,  1.190e-10,  1.341e-10,  1.506e-10,      &
   &  1.684e-10,  1.879e-10,  2.094e-10,  2.339e-10,  2.624e-10,      &
   &  2.956e-10,  3.345e-10,  3.795e-10,  4.309e-10,  4.883e-10,      &
   &  5.499e-10,  6.166e-10,  6.887e-10,  7.664e-10,  8.495e-10,      &
   &  9.372e-10,  1.028e-09,  1.119e-09,  1.209e-09,  1.293e-09,      &
   &  1.372e-09,  1.448e-09,  1.526e-09,  1.612e-09,  1.711e-09,      &
   &  1.822e-09,  1.944e-09,  2.070e-09,  2.175e-09,  2.246e-09,      &
   &  2.287e-09,  2.298e-09,  2.285e-09,  2.248e-09,  2.190e-09,      &
   &  2.111e-09,  2.016e-09,  1.906e-09,  1.785e-09,  1.657e-09,      &
   &  1.521e-09,  1.382e-09,  1.238e-09,  1.089e-09,  9.341e-10/
   DATA F3501/                                                       &
   &  7.733e-10,  6.100e-10,  4.528e-10,  3.101e-10,  1.906e-10,      &
   &  1.028e-10,  4.896e-11,  2.904e-11,  2.392e-11,  1.984e-11,      &
   &  1.661e-11,  1.403e-11,  1.193e-11,  1.023e-11,  8.815e-12,      &
   &  7.641e-12,  6.658e-12,  5.828e-12,  5.123e-12,  4.521e-12,      &
   &  4.004e-12,  3.558e-12,  3.171e-12,  2.834e-12,  2.539e-12,      &
   &  2.280e-12,  2.052e-12,  1.850e-12,  1.671e-12,  1.512e-12,      &
   &  1.370e-12,  1.243e-12,  1.129e-12,  1.026e-12,  9.339e-13,      &
   &  8.504e-13,  7.748e-13,  7.063e-13,  6.440e-13,  5.872e-13,      &
   &  5.354e-13,  4.881e-13,  4.448e-13,  4.051e-13,  3.686e-13,      &
   &  3.351e-13,  3.042e-13,  2.757e-13,  2.494e-13,  2.252e-13/
   DATA F3551/                                                       &
   &  2.027e-13,  1.818e-13,  1.625e-13,  1.446e-13,  1.280e-13,      &
   &  1.125e-13,  9.815e-14,  8.481e-14,  7.242e-14,  6.095e-14,      &
   &  5.034e-14,  4.060e-14,  3.174e-14,  2.390e-14,  1.717e-14,      &
   &  1.206e-14,  1.034e-14,  9.811e-15,  9.358e-15,  8.906e-15,      &
   &  8.453e-15,  8.075e-15,  7.698e-15,  7.328e-15,  6.989e-15,      &
   &  6.664e-15,  6.355e-15,  6.068e-15,  5.796e-15,  5.540e-15,      &
   &  5.298e-15,  5.064e-15,  4.845e-15,  4.642e-15,  4.445e-15,      &
   &  4.257e-15,  4.075e-15,  3.909e-15,  3.743e-15,  3.592e-15,      &
   &  3.449e-15,  3.306e-15,  3.177e-15,  3.057e-15,  2.936e-15,      &
   &  2.823e-15,  2.709e-15,  2.611e-15,  2.513e-15,  2.415e-15/
   DATA F3601/                                                       &
   &  2.325e-15,  2.242e-15,  2.158e-15,  2.083e-15,  2.008e-15,      &
   &  1.940e-15,  1.872e-15,  1.804e-15,  1.743e-15,  1.691e-15,      &
   &  1.638e-15,  1.585e-15,  1.540e-15,  1.494e-15,  1.457e-15,      &
   &  1.419e-15,  1.389e-15,  1.358e-15,  1.343e-15,  1.321e-15,      &
   &  1.298e-15,  1.283e-15,  1.268e-15,  1.253e-15,  1.245e-15,      &
   &  5.862e-15,  3.288e-14,  6.265e-14,  9.326e-14,  1.234e-13,      &
   &  1.522e-13,  1.796e-13,  2.058e-13,  2.317e-13,  2.575e-13,      &
   &  2.833e-13,  3.087e-13,  3.328e-13,  3.550e-13,  3.742e-13,      &
   &  3.895e-13,  4.003e-13,  4.058e-13,  4.058e-13,  4.003e-13,      &
   &  3.895e-13,  3.742e-13,  3.548e-13,  3.324e-13,  3.077e-13/
   DATA F3651/                                                       &
   &  2.818e-13,  2.549e-13,  2.274e-13,  1.990e-13,  1.691e-13,      &
   &  1.374e-13,  1.040e-13,  6.996e-14,  3.648e-14,  5.560e-15,      &
   &  6.709e-16,  6.468e-16,  6.196e-16,  5.917e-16,  5.645e-16,      &
   &  5.396e-16,  5.162e-16,  4.928e-16,  4.717e-16,  4.528e-16,      &
   &  4.385e-16,  4.279e-16,  4.211e-16,  4.181e-16,  4.166e-16,      &
   &  4.166e-16,  4.189e-16,  4.226e-16,  4.279e-16,  4.370e-16,      &
   &  4.475e-16,  4.611e-16,  4.755e-16,  4.936e-16,  5.155e-16,      &
   &  5.426e-16,  5.736e-16,  6.060e-16,  6.453e-16,  6.928e-16,      &
   &  2.296e-14,  5.327e-14,  8.601e-14,  1.202e-13,  1.548e-13,      &
   &  1.892e-13,  2.235e-13,  2.581e-13,  2.945e-13,  3.339e-13/
   DATA F3701/                                                       &
   &  3.774e-13,  4.261e-13,  4.803e-13,  5.411e-13,  6.089e-13,      &
   &  6.851e-13,  7.711e-13,  8.691e-13,  9.818e-13,  1.112e-12,      &
   &  1.263e-12,  1.436e-12,  1.634e-12,  1.856e-12,  2.099e-12,      &
   &  2.358e-12,  2.625e-12,  2.889e-12,  3.139e-12,  3.373e-12,      &
   &  3.595e-12,  3.811e-12,  4.032e-12,  4.260e-12,  4.497e-12,      &
   &  4.736e-12,  4.970e-12,  5.185e-12,  5.369e-12,  5.511e-12,      &
   &  5.600e-12,  5.630e-12,  5.600e-12,  5.512e-12,  5.369e-12,      &
   &  5.181e-12,  4.955e-12,  4.699e-12,  4.424e-12,  4.134e-12,      &
   &  3.833e-12,  3.520e-12,  3.193e-12,  2.847e-12,  2.481e-12,      &
   &  2.103e-12,  1.730e-12,  1.375e-12,  1.053e-12,  7.751e-13/
   DATA F3751/                                                       &
   &  5.442e-13,  3.613e-13,  2.236e-13,  1.239e-13,  5.619e-14,      &
   &  1.470e-14,  4.747e-15,  4.483e-15,  4.279e-15,  4.158e-15,      &
   &  4.136e-15,  4.204e-15,  4.355e-15,  1.840e-14,  4.496e-14,      &
   &  8.112e-14,  1.284e-13,  1.892e-13,  2.667e-13,  3.638e-13,      &
   &  4.862e-13,  6.386e-13,  8.266e-13,  1.056e-12,  1.334e-12,      &
   &  1.665e-12,  2.055e-12,  2.506e-12,  3.021e-12,  3.595e-12,      &
   &  4.225e-12,  4.900e-12,  5.605e-12,  6.318e-12,  7.018e-12,      &
   &  7.688e-12,  8.326e-12,  8.949e-12,  9.578e-12,  1.023e-11,      &
   &  1.089e-11,  1.157e-11,  1.224e-11,  1.286e-11,  1.339e-11,      &
   &  1.380e-11,  1.407e-11,  1.418e-11,  1.412e-11,  1.390e-11/
   DATA F3801/                                                       &
   &  1.354e-11,  1.306e-11,  1.247e-11,  1.181e-11,  1.108e-11,      &
   &  1.032e-11,  9.520e-12,  8.700e-12,  7.857e-12,  6.985e-12,      &
   &  6.077e-12,  5.144e-12,  4.213e-12,  3.324e-12,  2.518e-12,      &
   &  1.820e-12,  1.251e-12,  8.129e-13,  4.967e-13,  2.819e-13,      &
   &  1.460e-13,  6.531e-14,  1.977e-14,  6.506e-15,  5.547e-15,      &
   &  4.740e-15,  4.113e-15,  3.645e-15,  3.291e-15,  3.004e-15,      &
   &  2.755e-15,  2.536e-15,  2.392e-15,  2.302e-15,  2.249e-15,      &
   &  2.219e-15,  2.226e-15,  2.264e-15,  2.340e-15,  2.445e-15,      &
   &  2.581e-15,  8.149e-15,  4.005e-14,  8.285e-14,  1.376e-13,      &
   &  2.056e-13,  2.892e-13,  3.902e-13,  5.099e-13,  6.488e-13/
   DATA F3851/                                                       &
   &  8.062e-13,  9.798e-13,  1.166e-12,  1.358e-12,  1.551e-12,      &
   &  1.737e-12,  1.915e-12,  2.088e-12,  2.261e-12,  2.440e-12,      &
   &  2.626e-12,  2.818e-12,  3.010e-12,  3.195e-12,  3.365e-12,      &
   &  3.508e-12,  3.618e-12,  3.689e-12,  3.718e-12,  3.705e-12,      &
   &  3.652e-12,  3.564e-12,  3.447e-12,  3.305e-12,  3.146e-12,      &
   &  2.973e-12,  2.790e-12,  2.600e-12,  2.400e-12,  2.189e-12,      &
   &  1.963e-12,  1.722e-12,  1.473e-12,  1.226e-12,  9.909e-13,      &
   &  7.778e-13,  5.938e-13,  4.407e-13,  3.182e-13,  2.235e-13,      &
   &  1.522e-13,  9.812e-14,  5.676e-14,  2.415e-14,  1.683e-15,      &
   &  1.434e-15,  1.238e-15,  1.079e-15,  9.434e-16,  8.302e-16/
   DATA F3901/                                                       &
   &  7.306e-16,  6.400e-16,  5.691e-16,  5.117e-16,  4.626e-16,      &
   &  4.166e-16,  3.743e-16,  3.411e-16,  3.132e-16,  2.891e-16,      &
   &  2.664e-16,  2.453e-16,  2.279e-16,  2.136e-16,  2.008e-16,      &
   &  1.887e-16,  1.781e-16,  1.698e-16,  1.623e-16,  1.562e-16,      &
   &  1.509e-16,  1.472e-16,  1.434e-16,  1.419e-16,  1.411e-16,      &
   &  1.419e-16,  1.434e-16,  1.457e-16,  1.502e-16,  1.570e-16,      &
   &  1.668e-16,  1.811e-16,  1.985e-16,  2.181e-16,  2.415e-16,      &
   &  2.664e-16,  2.966e-16,  3.313e-16,  3.675e-16,  4.030e-16,      &
   &  4.347e-16,  4.619e-16,  4.875e-16,  5.117e-16,  5.336e-16,      &
   &  5.540e-16,  5.706e-16,  5.849e-16,  5.977e-16,  6.113e-16/
   DATA F3951/                                                       &
   &  1.061e-14,  2.170e-14,  3.323e-14,  4.525e-14,  5.762e-14,      &
   &  7.011e-14,  8.221e-14,  9.352e-14,  1.034e-13,  1.115e-13,      &
   &  1.176e-13,  1.216e-13,  1.234e-13,  1.229e-13,  1.206e-13,      &
   &  1.170e-13,  1.127e-13,  1.081e-13,  1.038e-13,  1.002e-13,      &
   &  9.738e-14,  9.552e-14,  9.443e-14,  9.389e-14,  9.366e-14,      &
   &  9.377e-14,  9.433e-14,  9.567e-14,  9.768e-14,  1.005e-13,      &
   &  1.046e-13,  1.106e-13,  1.191e-13,  1.304e-13,  1.446e-13,      &
   &  1.611e-13,  1.789e-13,  1.960e-13,  2.106e-13,  2.203e-13,      &
   &  2.255e-13,  2.263e-13,  2.231e-13,  2.165e-13,  2.068e-13,      &
   &  1.945e-13,  1.800e-13,  1.637e-13,  1.462e-13,  1.278e-13/
   DATA F4001/                                                       &
   &  1.088e-13,  8.916e-14,  6.893e-14,  4.796e-14,  2.634e-14,      &
   &  4.974e-15,  8.151e-16,  8.302e-16,  8.453e-16,  8.679e-16,      &
   &  8.906e-16,  9.132e-16,  9.358e-16,  9.660e-16,  1.004e-15,      &
   &  1.057e-15,  1.125e-15,  1.200e-15,  1.291e-15,  1.389e-15,      &
   &  6.595e-15,  2.605e-14,  4.954e-14,  7.688e-14,  1.083e-13,      &
   &  1.441e-13,  1.841e-13,  2.286e-13,  2.772e-13,  3.294e-13,      &
   &  3.845e-13,  4.415e-13,  4.993e-13,  5.575e-13,  6.175e-13,      &
   &  6.813e-13,  7.524e-13,  8.335e-13,  9.256e-13,  1.029e-12,      &
   &  1.141e-12,  1.253e-12,  1.358e-12,  1.460e-12,  1.561e-12,      &
   &  1.663e-12,  1.770e-12,  1.884e-12,  2.005e-12,  2.136e-12/
   DATA F4051/                                                       &
   &  2.276e-12,  2.426e-12,  2.585e-12,  2.750e-12,  2.918e-12,      &
   &  3.082e-12,  3.243e-12,  3.398e-12,  3.558e-12,  3.735e-12,      &
   &  3.947e-12,  4.208e-12,  4.525e-12,  4.897e-12,  5.316e-12,      &
   &  5.752e-12,  6.174e-12,  6.589e-12,  7.013e-12,  7.462e-12,      &
   &  7.957e-12,  8.518e-12,  9.167e-12,  9.929e-12,  1.082e-11,      &
   &  1.187e-11,  1.309e-11,  1.449e-11,  1.607e-11,  1.782e-11,      &
   &  1.972e-11,  2.172e-11,  2.376e-11,  2.578e-11,  2.771e-11,      &
   &  2.956e-11,  3.140e-11,  3.336e-11,  3.556e-11,  3.809e-11,      &
   &  4.093e-11,  4.398e-11,  4.703e-11,  4.986e-11,  5.188e-11,      &
   &  5.313e-11,  5.366e-11,  5.352e-11,  5.277e-11,  5.147e-11/
   DATA F4101/                                                       &
   &  4.970e-11,  4.753e-11,  4.503e-11,  4.228e-11,  3.936e-11,      &
   &  3.631e-11,  3.316e-11,  2.993e-11,  2.661e-11,  2.317e-11,      &
   &  1.961e-11,  1.604e-11,  1.264e-11,  9.610e-12,  7.158e-12,      &
   &  5.460e-12,  4.584e-12,  4.502e-12,  5.024e-12,  5.734e-12,      &
   &  6.605e-12,  7.649e-12,  8.877e-12,  1.031e-11,  1.195e-11,      &
   &  1.382e-11,  1.593e-11,  1.828e-11,  2.090e-11,  2.377e-11,      &
   &  2.691e-11,  3.027e-11,  3.380e-11,  3.740e-11,  4.095e-11,      &
   &  4.430e-11,  4.743e-11,  5.043e-11,  5.352e-11,  5.687e-11,      &
   &  6.059e-11,  6.473e-11,  6.903e-11,  7.342e-11,  7.661e-11,      &
   &  7.866e-11,  7.963e-11,  7.959e-11,  7.862e-11,  7.683e-11/
   DATA F4151/                                                       &
   &  7.428e-11,  7.108e-11,  6.733e-11,  6.313e-11,  5.862e-11,      &
   &  5.389e-11,  4.904e-11,  4.412e-11,  3.915e-11,  3.410e-11,      &
   &  2.889e-11,  2.350e-11,  1.814e-11,  1.303e-11,  8.549e-12,      &
   &  4.950e-12,  2.450e-12,  1.118e-12,  8.048e-13,  6.484e-13,      &
   &  5.257e-13,  4.281e-13,  3.496e-13,  2.858e-13,  2.334e-13,      &
   &  1.900e-13,  1.538e-13,  1.235e-13,  9.783e-14,  7.609e-14,      &
   &  5.762e-14,  4.193e-14,  2.872e-14,  1.785e-14,  9.478e-15,      &
   &  6.702e-15,  6.045e-15,  5.502e-15,  5.034e-15,  4.611e-15,      &
   &  4.219e-15,  3.857e-15,  3.547e-15,  3.268e-15,  3.019e-15,      &
   &  2.792e-15,  2.574e-15,  2.385e-15,  2.219e-15,  2.060e-15/
   DATA F4201/                                                       &
   &  1.917e-15,  1.781e-15,  1.660e-15,  1.555e-15,  1.457e-15,      &
   &  1.358e-15,  1.268e-15,  1.192e-15,  1.117e-15,  1.049e-15,      &
   &  9.887e-16,  9.283e-16,  8.755e-16,  8.226e-16,  7.774e-16,      &
   &  7.343e-16,  6.928e-16,  6.551e-16,  6.204e-16,  5.872e-16,      &
   &  5.562e-16,  5.268e-16,  4.996e-16,  4.740e-16,  4.506e-16,      &
   &  4.279e-16,  4.060e-16,  3.864e-16,  3.675e-16,  3.502e-16,      &
   &  3.336e-16,  3.177e-16,  3.026e-16,  2.891e-16,  2.762e-16,      &
   &  2.634e-16,  2.513e-16,  2.400e-16,  2.302e-16,  2.196e-16,      &
   &  2.106e-16,  2.015e-16,  1.925e-16,  1.849e-16,  1.766e-16,      &
   &  1.698e-16,  1.623e-16,  1.562e-16,  1.494e-16,  1.434e-16/
   DATA F4251/                                                       &
   &  1.381e-16,  1.321e-16,  1.275e-16,  1.223e-16,  1.177e-16,      &
   &  1.132e-16,  1.087e-16,  1.049e-16,  1.004e-16,  9.660e-17,      &
   &  9.358e-17,  8.981e-17,  8.679e-17,  8.377e-17,  8.075e-17,      &
   &  7.774e-17,  7.479e-17,  7.215e-17,  6.966e-17,  6.725e-17,      &
   &  6.491e-17,  6.264e-17,  6.053e-17,  5.849e-17,  5.653e-17,      &
   &  5.457e-17,  5.275e-17,  5.102e-17,  4.936e-17,  4.777e-17,      &
   &  4.619e-17,  4.468e-17,  4.325e-17,  4.189e-17,  4.053e-17,      &
   &  3.925e-17,  3.804e-17,  3.683e-17,  3.570e-17,  3.457e-17,      &
   &  3.351e-17,  3.245e-17,  3.147e-17,  3.057e-17,  2.958e-17,      &
   &  2.875e-17,  2.785e-17,  2.702e-17,  2.626e-17,  2.551e-17/
   DATA F4301/                                                       &
   &  2.475e-17,  2.400e-17,  2.332e-17,  2.264e-17,  2.196e-17,      &
   &  2.136e-17,  2.075e-17,  2.015e-17,  1.962e-17,  1.902e-17,      &
   &  1.849e-17,  1.804e-17,  1.751e-17,  1.706e-17,  1.660e-17,      &
   &  1.615e-17,  1.570e-17,  1.525e-17,  1.487e-17,  1.442e-17,      &
   &  1.404e-17,  1.366e-17,  1.336e-17,  1.298e-17,  1.260e-17,      &
   &  1.230e-17,  1.200e-17,  1.170e-17,  1.140e-17,  1.109e-17,      &
   &  1.079e-17,  1.057e-17,  1.026e-17,  1.004e-17,  9.736e-18,      &
   &  9.509e-18,  9.283e-18,  9.057e-18,  8.830e-18,  8.604e-18,      &
   &  8.377e-18,  8.226e-18,  8.000e-18,  7.774e-18,  7.623e-18,      &
   &  7.419e-18,  7.245e-18,  7.072e-18,  6.906e-18,  6.740e-18/
   DATA F4351/                                                       &
   &  6.581e-18,  6.423e-18,  6.279e-18,  6.128e-18,  5.985e-18,      &
   &  5.849e-18,  5.713e-18,  5.585e-18,  5.457e-18,  5.328e-18,      &
   &  5.208e-18,  5.087e-18,  4.974e-18,  4.860e-18,  4.755e-18,      &
   &  4.649e-18,  4.543e-18,  4.445e-18,  4.347e-18,  4.249e-18,      &
   &  4.158e-18,  4.068e-18,  3.977e-18,  3.887e-18,  3.804e-18,      &
   &  3.721e-18,  3.645e-18,  3.562e-18,  3.487e-18,  3.411e-18,      &
   &  3.343e-18,  3.268e-18,  3.200e-18,  3.132e-18,  3.072e-18,      &
   &  3.004e-18,  2.943e-18,  2.883e-18,  2.823e-18,  2.762e-18,      &
   &  2.709e-18,  2.649e-18,  2.596e-18,  2.543e-18,  2.491e-18,      &
   &  2.445e-18,  2.392e-18,  2.347e-18,  2.302e-18,  2.257e-18/
   DATA F4401/                                                       &
   &  2.211e-18,  2.166e-18,  2.121e-18,  2.083e-18,  2.038e-18,      &
   &  2.000e-18,  1.962e-18,  1.925e-18,  1.887e-18,  1.849e-18,      &
   &  1.811e-18,  1.781e-18,  1.743e-18,  1.713e-18,  1.683e-18,      &
   &  1.645e-18,  1.615e-18,  1.585e-18,  1.555e-18,  1.525e-18,      &
   &  1.502e-18,  1.472e-18,  1.442e-18,  1.419e-18,  1.389e-18,      &
   &  1.366e-18,  1.343e-18,  1.321e-18,  1.291e-18,  1.268e-18,      &
   &  1.245e-18,  1.223e-18,  1.208e-18,  1.185e-18,  1.162e-18,      &
   &  1.140e-18,  1.125e-18,  1.102e-18,  1.087e-18,  1.064e-18,      &
   &  1.049e-18,  1.026e-18,  1.011e-18,  9.962e-19,  9.811e-19,      &
   &  9.585e-19,  9.434e-19,  9.283e-19,  9.132e-19,  8.981e-19/
   DATA F4451/                                                       &
   &  8.830e-19,  8.679e-19,  8.604e-19,  8.453e-19,  8.302e-19,      &
   &  8.151e-19,  8.075e-19,  7.925e-19,  7.774e-19,  7.698e-19,      &
   &  7.547e-19,  7.442e-19,  7.328e-19,  7.215e-19,  7.109e-19,      &
   &  7.011e-19,  6.906e-19,  6.808e-19,  6.709e-19,  6.611e-19,      &
   &  6.521e-19,  6.430e-19,  6.340e-19,  6.257e-19,  6.174e-19,      &
   &  6.091e-19,  6.008e-19,  5.932e-19,  5.857e-19,  5.781e-19,      &
   &  5.706e-19,  5.638e-19,  5.570e-19,  5.502e-19,  5.442e-19,      &
   &  5.374e-19,  5.313e-19,  5.253e-19,  5.200e-19,  5.140e-19,      &
   &  5.087e-19,  5.034e-19,  4.981e-19,  4.936e-19,  4.891e-19,      &
   &  4.845e-19,  4.800e-19,  4.755e-19,  4.717e-19,  4.679e-19/
   DATA F4501/                                                       &
   &  4.642e-19,  4.604e-19,  4.566e-19,  4.536e-19,  4.506e-19,      &
   &  4.475e-19,  4.445e-19,  4.423e-19,  4.400e-19,  4.370e-19,      &
   &  4.355e-19,  4.332e-19,  4.317e-19,  4.294e-19,  4.287e-19,      &
   &  4.272e-19,  4.257e-19,  4.249e-19,  4.242e-19,  4.234e-19,      &
   &  4.234e-19,  4.226e-19,  4.226e-19,  4.226e-19,  4.234e-19,      &
   &  4.234e-19,  4.242e-19,  4.249e-19,  4.264e-19,  4.272e-19,      &
   &  4.287e-19,  4.302e-19,  4.325e-19,  4.347e-19,  4.370e-19,      &
   &  4.392e-19,  4.423e-19,  4.453e-19,  4.483e-19,  4.521e-19,      &
   &  4.558e-19,  4.596e-19,  4.642e-19,  4.687e-19,  4.732e-19,      &
   &  4.785e-19,  4.838e-19,  4.898e-19,  4.958e-19,  5.026e-19/
   DATA F4551/                                                       &
   &  5.094e-19,  5.162e-19,  5.238e-19,  5.313e-19,  5.396e-19,      &
   &  5.487e-19,  5.577e-19,  5.668e-19,  5.774e-19,  5.872e-19,      &
   &  5.985e-19,  6.098e-19,  6.219e-19,  6.347e-19,  6.475e-19,      &
   &  6.611e-19,  6.755e-19,  6.906e-19,  7.064e-19,  7.223e-19,      &
   &  7.396e-19,  7.547e-19,  7.774e-19,  7.925e-19,  8.151e-19,      &
   &  8.377e-19,  8.604e-19,  8.830e-19,  9.057e-19,  9.283e-19,      &
   &  9.585e-19,  9.887e-19,  1.011e-18,  1.042e-18,  1.072e-18,      &
   &  1.109e-18,  1.140e-18,  1.177e-18,  1.215e-18,  1.253e-18,      &
   &  1.298e-18,  1.343e-18,  1.389e-18,  1.434e-18,  1.479e-18,      &
   &  1.532e-18,  1.592e-18,  1.645e-18,  1.706e-18,  1.774e-18/
   DATA F4601/                                                       &
   &  1.842e-18,  1.909e-18,  1.985e-18,  2.060e-18,  2.143e-18,      &
   &  2.226e-18,  2.325e-18,  2.423e-18,  2.521e-18,  2.634e-18,      &
   &  2.747e-18,  2.868e-18,  3.004e-18,  3.132e-18,  3.283e-18,      &
   &  3.434e-18,  3.600e-18,  3.774e-18,  3.955e-18,  4.151e-18,      &
   &  4.362e-18,  4.589e-18,  4.823e-18,  5.072e-18,  5.343e-18,      &
   &  5.638e-18,  5.947e-18,  6.279e-18,  6.634e-18,  7.011e-18,      &
   &  7.434e-18,  7.925e-18,  8.377e-18,  8.906e-18,  9.434e-18,      &
   &  1.004e-17,  1.072e-17,  1.147e-17,  1.223e-17,  1.306e-17,      &
   &  1.404e-17,  1.509e-17,  1.615e-17,  1.736e-17,  1.872e-17,      &
   &  2.030e-17,  2.204e-17,  2.377e-17,  2.581e-17,  2.815e-17/
   DATA F4651/                                                       &
   &  3.094e-17,  3.389e-17,  3.706e-17,  4.068e-17,  4.498e-17,      &
   &  5.042e-17,  5.638e-17,  6.279e-17,  7.019e-17,  7.925e-17,      &
   &  9.132e-17,  1.049e-16,  1.208e-16,  1.389e-16,  1.592e-16,      &
   &  1.811e-16,  2.060e-16,  2.340e-16,  2.634e-16,  2.943e-16,      &
   &  3.260e-16,  3.600e-16,  3.955e-16,  4.309e-16,  1.103e-14,      &
   &  2.506e-14,  4.104e-14,  5.874e-14,  7.793e-14,  9.821e-14,      &
   &  1.190e-13,  1.397e-13,  1.596e-13,  1.784e-13,  1.965e-13,      &
   &  2.147e-13,  2.340e-13,  2.551e-13,  2.783e-13,  3.032e-13,      &
   &  3.282e-13,  3.517e-13,  3.712e-13,  3.839e-13,  3.902e-13,      &
   &  3.907e-13,  3.857e-13,  3.758e-13,  3.617e-13,  3.440e-13/
   DATA F4701/                                                       &
   &  3.234e-13,  3.005e-13,  2.758e-13,  2.500e-13,  2.232e-13,      &
   &  1.957e-13,  1.672e-13,  1.376e-13,  1.068e-13,  7.534e-14,      &
   &  4.448e-14,  1.646e-14,  7.117e-16,  6.853e-16,  6.611e-16,      &
   &  6.423e-16,  6.242e-16,  6.060e-16,  5.925e-16,  5.887e-16,      &
   &  6.008e-16,  6.279e-16,  6.679e-16,  7.192e-16,  7.849e-16,      &
   &  8.604e-16,  9.434e-16,  1.237e-14,  2.704e-14,  4.486e-14,      &
   &  6.621e-14,  9.162e-14,  1.215e-13,  1.567e-13,  1.979e-13,      &
   &  2.458e-13,  3.012e-13,  3.647e-13,  4.369e-13,  5.184e-13,      &
   &  6.090e-13,  7.084e-13,  8.155e-13,  9.281e-13,  1.044e-12,      &
   &  1.160e-12,  1.272e-12,  1.377e-12,  1.480e-12,  1.584e-12/
   DATA F4751/                                                       &
   &  1.698e-12,  1.828e-12,  1.975e-12,  2.140e-12,  2.312e-12,      &
   &  2.455e-12,  2.557e-12,  2.620e-12,  2.647e-12,  2.640e-12,      &
   &  2.602e-12,  2.536e-12,  2.446e-12,  2.334e-12,  2.206e-12,      &
   &  2.063e-12,  1.911e-12,  1.750e-12,  1.584e-12,  1.413e-12,      &
   &  1.237e-12,  1.054e-12,  8.642e-13,  6.704e-13,  4.821e-13,      &
   &  3.093e-13,  1.651e-13,  5.853e-14,  2.657e-15,  2.340e-15,      &
   &  2.068e-15,  1.849e-15,  1.645e-15,  1.472e-15,  1.321e-15,      &
   &  1.208e-15,  1.132e-15,  1.094e-15,  1.079e-15,  4.198e-15,      &
   &  1.473e-14,  2.792e-14,  4.378e-14,  6.256e-14,  8.443e-14,      &
   &  1.096e-13,  1.382e-13,  1.702e-13,  2.055e-13,  2.437e-13/
   DATA F4801/                                                       &
   &  2.843e-13,  3.262e-13,  3.682e-13,  4.090e-13,  4.477e-13,      &
   &  4.847e-13,  5.220e-13,  5.617e-13,  6.059e-13,  6.555e-13,      &
   &  7.100e-13,  7.681e-13,  8.206e-13,  8.615e-13,  8.884e-13,      &
   &  9.020e-13,  9.031e-13,  8.928e-13,  8.720e-13,  8.425e-13,      &
   &  8.051e-13,  7.613e-13,  7.123e-13,  6.594e-13,  6.037e-13,      &
   &  5.458e-13,  4.862e-13,  4.248e-13,  3.610e-13,  2.945e-13,      &
   &  2.263e-13,  1.588e-13,  9.623e-14,  4.206e-14,  1.162e-15,      &
   &  1.064e-15,  9.509e-16,  8.226e-16,  6.913e-16,  5.683e-16,      &
   &  4.687e-16,  3.955e-16,  3.404e-16,  2.974e-16,  2.619e-16,      &
   &  2.287e-16,  1.985e-16,  1.758e-16,  1.570e-16,  1.419e-16/
   DATA F4851/                                                       &
   &  1.268e-16,  1.132e-16,  1.019e-16,  9.283e-17,  8.453e-17,      &
   &  7.698e-17,  6.966e-17,  6.362e-17,  5.842e-17,  5.381e-17,      &
   &  4.943e-17,  4.536e-17,  4.189e-17,  3.879e-17,  3.600e-17,      &
   &  3.336e-17,  3.087e-17,  2.875e-17,  2.679e-17,  2.506e-17,      &
   &  2.332e-17,  2.181e-17,  2.038e-17,  1.909e-17,  1.796e-17,      &
   &  1.683e-17,  1.585e-17,  1.487e-17,  1.404e-17,  1.321e-17,      &
   &  1.245e-17,  1.177e-17,  1.109e-17,  1.049e-17,  9.962e-18,      &
   &  9.434e-18,  8.906e-18,  8.453e-18,  8.000e-18,  7.623e-18,      &
   &  7.238e-18,  6.875e-18,  6.543e-18,  6.226e-18,  5.932e-18,      &
   &  5.653e-18,  5.381e-18,  5.132e-18,  4.898e-18,  4.679e-18/
   DATA F4901/                                                       &
   &  4.468e-18,  4.272e-18,  4.083e-18,  3.902e-18,  3.736e-18,      &
   &  3.577e-18,  3.426e-18,  3.283e-18,  3.147e-18,  3.019e-18,      &
   &  2.891e-18,  2.777e-18,  2.664e-18,  2.558e-18,  2.460e-18,      &
   &  2.362e-18,  2.272e-18,  2.181e-18,  2.098e-18,  2.023e-18,      &
   &  1.947e-18,  1.872e-18,  1.804e-18,  1.736e-18,  1.675e-18,      &
   &  1.615e-18,  1.555e-18,  1.502e-18,  1.449e-18,  1.396e-18,      &
   &  1.351e-18,  1.298e-18,  1.253e-18,  1.215e-18,  1.170e-18,      &
   &  1.132e-18,  1.094e-18,  1.057e-18,  1.026e-18,  9.887e-19,      &
   &  9.585e-19,  9.283e-19,  8.981e-19,  8.679e-19,  8.377e-19,      &
   &  8.151e-19,  7.849e-19,  7.623e-19,  7.396e-19,  7.162e-19/
   DATA F4951/                                                       &
   &  6.943e-19,  6.732e-19,  6.528e-19,  6.332e-19,  6.143e-19,      &
   &  5.955e-19,  5.781e-19,  5.608e-19,  5.442e-19,  5.283e-19,      &
   &  5.132e-19,  4.981e-19,  4.838e-19,  4.702e-19,  4.566e-19,      &
   &  4.438e-19,  4.309e-19,  4.189e-19,  4.075e-19,  3.962e-19,      &
   &  3.849e-19,  3.743e-19,  3.645e-19,  3.547e-19,  3.449e-19,      &
   &  3.358e-19,  3.268e-19,  3.177e-19,  3.094e-19,  3.011e-19,      &
   &  2.936e-19,  2.853e-19,  2.785e-19,  2.709e-19,  2.642e-19,      &
   &  2.574e-19,  2.506e-19,  2.438e-19,  2.377e-19,  2.317e-19,      &
   &  2.257e-19,  2.204e-19,  2.151e-19,  2.098e-19,  2.045e-19,      &
   &  1.992e-19,  1.970e-19,  1.970e-19,  1.962e-19,  1.902e-19/
   DATA F5001/                                                       &
   &  1.887e-19/
!
end block data BFCO2
!
!     --------------------------------------------------------------
!
SUBROUTINE xn2_r (V1C,V2C,DVC,NPTC,C,fo2,Tave,v1ss,v2ss)
!
!     Model used:
!      Borysow, A, and L. Frommhold, "Collision-induced
!         rototranslational absorption spectra of N2-N2
!         pairs for temperatures from 50 to 300 K", The
!         Astrophysical Journal, 311, 1043-1057, 1986.
!
!     Updated 2004/09/22 based on:
!
!      Boissoles, J., C. Boulet, R.H. Tipping, A. Brown and Q. Ma,
!         Theoretical CAlculations of the Translation-Rotation
!         Collision-Induced Absorption in N2-N2, O2-O2 and N2-O2 Pairs,
!         J.Quant. Spec. Rad. Transfer, 82,505 (2003).
!
!     The scale factors are reported to account for the efect of o2-o2
!     and n2-o2 collision induced effects.
!     The values for scale factor values (sf296) for 296K are based on
!     linear interpolation of Boissoles at al. values at 250K and 300K
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/  V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   COMMON /N2RT296/ V1S,V2S,DVS,NPTS,C_296(73),sf_296(73)
   COMMON /N2RT220/ V1b,V2b,DVb,NPTb,C_220(73),sf_220(73)
   DIMENSION C(*),fo2(*)
!
   data xo2 / 0.21/, xn2 / 0.79/, T_296 / 296./, T_220 / 220./

   tfac = (TAVE-T_296)/(T_220-T_296)

   DVC = DVS
   v1ss = v1s
   v2ss = v2s
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
!*******  ABSORPTION COEFFICIENT IN UNITS OF CM-1 AMAGAT-2
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      C(J) = C_296(i)*(( C_220(i)/ C_296(i))**tfac)
      sf_T = sf_296(i)*((sf_220(i)/sf_296(i))**tfac)

!        correct for incorporation of air mixing ratios in sf
!        fo2 is now ~ the ratio of alpha(n2-o2)/alpha(n2-n2)
!        Eq's 7 and 8 in the Boissoles paper.

!        fo2(i) = (sf_T - 1.)*(xn2**2)/(xn2*xo2)
      fo2(j) = (sf_T - 1.)*(xn2)/(xo2)

10 END DO
!
   RETURN
!
end subroutine xn2_r
!
BLOCK DATA BN2T296
!
   IMPLICIT REAL*8           (V)
!
!*******  ABSORPTION COEFFICIENT IN UNITS OF CM-1 AMAGAT-2
!
!           THESE DATA ARE FOR 296K
!
   COMMON /N2RT296/ V1N2CR,V2N2CR,DVN2CR,NPTN2C,CT296(73),sf_296(73)
!
   DATA V1N2CR,V2N2CR,DVN2CR,NPTN2C /                                &
   &      -10.,  350.,  5.0,   73 /
!
   DATA CT296/                                                       &
   &     0.4303E-06, 0.4850E-06, 0.4979E-06, 0.4850E-06, 0.4303E-06,  &
   &     0.3715E-06, 0.3292E-06, 0.3086E-06, 0.2920E-06, 0.2813E-06,  &
   &     0.2804E-06, 0.2738E-06, 0.2726E-06, 0.2724E-06, 0.2635E-06,  &
   &     0.2621E-06, 0.2547E-06, 0.2428E-06, 0.2371E-06, 0.2228E-06,  &
   &     0.2100E-06, 0.1991E-06, 0.1822E-06, 0.1697E-06, 0.1555E-06,  &
   &     0.1398E-06, 0.1281E-06, 0.1138E-06, 0.1012E-06, 0.9078E-07,  &
   &     0.7879E-07, 0.6944E-07, 0.6084E-07, 0.5207E-07, 0.4540E-07,  &
   &     0.3897E-07, 0.3313E-07, 0.2852E-07, 0.2413E-07, 0.2045E-07,  &
   &     0.1737E-07, 0.1458E-07, 0.1231E-07, 0.1031E-07, 0.8586E-08,  &
   &     0.7162E-08, 0.5963E-08, 0.4999E-08, 0.4226E-08, 0.3607E-08,  &
   &     0.3090E-08, 0.2669E-08, 0.2325E-08, 0.2024E-08, 0.1783E-08,  &
   &     0.1574E-08, 0.1387E-08, 0.1236E-08, 0.1098E-08, 0.9777E-09,  &
   &     0.8765E-09, 0.7833E-09, 0.7022E-09, 0.6317E-09, 0.5650E-09,  &
   &     0.5100E-09, 0.4572E-09, 0.4115E-09, 0.3721E-09, 0.3339E-09,  &
   &     0.3005E-09, 0.2715E-09, 0.2428E-09/
!
   DATA sf_296/                                                      &
   &         1.3534,     1.3517,     1.3508,     1.3517,     1.3534,  &
   &         1.3558,     1.3584,     1.3607,     1.3623,     1.3632,  &
   &         1.3634,     1.3632,     1.3627,     1.3620,     1.3612,  &
   &         1.3605,     1.3597,     1.3590,     1.3585,     1.3582,  &
   &         1.3579,     1.3577,     1.3577,     1.3580,     1.3586,  &
   &         1.3594,     1.3604,     1.3617,     1.3633,     1.3653,  &
   &         1.3677,     1.3706,     1.3742,     1.3780,     1.3822,  &
   &         1.3868,     1.3923,     1.3989,     1.4062,     1.4138,  &
   &         1.4216,     1.4298,     1.4388,     1.4491,     1.4604,  &
   &         1.4718,     1.4829,     1.4930,     1.5028,     1.5138,  &
   &         1.5265,     1.5392,     1.5499,     1.5577,     1.5639,  &
   &         1.5714,     1.5816,     1.5920,     1.6003,     1.6051,  &
   &         1.6072,     1.6097,     1.6157,     1.6157,     1.6157,  &
   &         1.6157,     1.6157,     1.6157,     1.6157,     1.6157,  &
   &         1.6157,     1.6157,     1.6157/
!
end block data BN2T296
!
BLOCK DATA BN2T220
!
   IMPLICIT REAL*8           (V)
!
!*******  ABSORPTION COEFFICIENT IN UNITS OF CM-1 AMAGAT-2
!
!         THESE DATA ARE FOR 220K
!
   COMMON /N2RT220/ V1N2CR,V2N2CR,DVN2CR,NPTN2C,CT220(73),sf_220(73)
!
   DATA V1N2CR,V2N2CR,DVN2CR,NPTN2C / -10., 350., 5.0, 73 /
!
   DATA CT220/                                                       &
   &     0.4946E-06, 0.5756E-06, 0.5964E-06, 0.5756E-06, 0.4946E-06,  &
   &     0.4145E-06, 0.3641E-06, 0.3482E-06, 0.3340E-06, 0.3252E-06,  &
   &     0.3299E-06, 0.3206E-06, 0.3184E-06, 0.3167E-06, 0.2994E-06,  &
   &     0.2943E-06, 0.2794E-06, 0.2582E-06, 0.2468E-06, 0.2237E-06,  &
   &     0.2038E-06, 0.1873E-06, 0.1641E-06, 0.1474E-06, 0.1297E-06,  &
   &     0.1114E-06, 0.9813E-07, 0.8309E-07, 0.7059E-07, 0.6068E-07,  &
   &     0.5008E-07, 0.4221E-07, 0.3537E-07, 0.2885E-07, 0.2407E-07,  &
   &     0.1977E-07, 0.1605E-07, 0.1313E-07, 0.1057E-07, 0.8482E-08,  &
   &     0.6844E-08, 0.5595E-08, 0.4616E-08, 0.3854E-08, 0.3257E-08,  &
   &     0.2757E-08, 0.2372E-08, 0.2039E-08, 0.1767E-08, 0.1548E-08,  &
   &     0.1346E-08, 0.1181E-08, 0.1043E-08, 0.9110E-09, 0.8103E-09,  &
   &     0.7189E-09, 0.6314E-09, 0.5635E-09, 0.4976E-09, 0.4401E-09,  &
   &     0.3926E-09, 0.3477E-09, 0.3085E-09, 0.2745E-09, 0.2416E-09,  &
   &     0.2155E-09, 0.1895E-09, 0.1678E-09, 0.1493E-09, 0.1310E-09,  &
   &     0.1154E-09, 0.1019E-09, 0.8855E-10/
!
   DATA sf_220/                                                      &
   &         1.3536,     1.3515,     1.3502,     1.3515,     1.3536,  &
   &         1.3565,     1.3592,     1.3612,     1.3623,     1.3626,  &
   &         1.3623,     1.3616,     1.3609,     1.3600,     1.3591,  &
   &         1.3583,     1.3576,     1.3571,     1.3571,     1.3572,  &
   &         1.3574,     1.3578,     1.3585,     1.3597,     1.3616,  &
   &         1.3640,     1.3666,     1.3698,     1.3734,     1.3776,  &
   &         1.3828,     1.3894,     1.3969,     1.4049,     1.4127,  &
   &         1.4204,     1.4302,     1.4427,     1.4562,     1.4687,  &
   &         1.4798,     1.4894,     1.5000,     1.5142,     1.5299,  &
   &         1.5441,     1.5555,     1.5615,     1.5645,     1.5730,  &
   &         1.5880,     1.6028,     1.6121,     1.6133,     1.6094,  &
   &         1.6117,     1.6244,     1.6389,     1.6485,     1.6513,  &
   &         1.6468,     1.6438,     1.6523,     1.6523,     1.6523,  &
   &         1.6523,     1.6523,     1.6523,     1.6523,     1.6523,  &
   &         1.6523,     1.6523,     1.6523/
!
end block data BN2T220
!
!     --------------------------------------------------------------
!
subroutine n2_ver_1 (v1c,v2c,dvc,nptc,c,c1,c2,T,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8 (v)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
!
   COMMON /n2_f/V1S,V2S,DVS,NPTS,xn2_272(228),xn2_228(228),a_h2o(228)
!
   dimension c(*),c1(*),c2(*)
!
!     Nitrogen Collision Induced Fundamental

!     Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and J._M.
!        Hartmann, Infrared collision-induced absorption by N2 near 4.3
!        microns for atmospheric applications: measurements and
!         emprirical modeling, Appl. Optics, 35, 5911-5917, (1996).
!
!     mt_ckd_2.8: Coefficients for N2-H2O relative efficiency determined
!     by Mlawer and Alvarado based primarily on measurements from
!     Baranov and Lafferty (2012). These coefficients were derived
!     simultaneously with water vapor foreign and self continuum
!     coefficients from 1800-2600 cm-1 using IASI measurements as in
!     Alvarado et al. (2012).
   DATA  T_272/ 272./, T_228/ 228./
!
   xtfac  = ((1./T)-(1./T_272))/((1./T_228)-(1./T_272))
   xt_lin = (T-T_272)/(T_228-T_272)
!
!     a_o2  represents the relative broadening efficiency of o2
   a_o2  = 1.294 - 0.4545*T/296.

!     a_h2o represents the relative broadening efficiency of h2o.  It
!     has spectral dependence and is stored on same grid as xn2.

!     The absorption coefficients from the Lafferty et al. reference
!     are for pure nitrogen (absorber and broadener)
!
   DVC = DVS
   v1ss = v1s
   v2ss = v2s
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   do 10 j=1,nptc
      i = i1+(j-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      VJ = V1C+DVC* REAL(J-1)
!
      if ((xn2_272(i).gt.0.) .and. (xn2_228(i) .gt. 0.)) then
!           logarithmic interpolation in reciprical of temperature

         c(j) = xn2_272(i) *(xn2_228(i)/xn2_272(i))**xtfac

      else
!           linear interpolation  (xn2_272 or xn2_228 = 0 to get here)

         c(j) = xn2_272(i) + (xn2_228(i)-xn2_272(i))*xt_lin

      endif

!     the radiation field is removed with 1/vj

      c(j) = c(j)/vj
      c1(j) = a_o2 * c(j)
!     the factor 9/7 is a modification to a preliminary formulation
!     of n2-h2o absorption.
      c2(j) = (9./7.) * a_h2o(i) * c(j)
!
10 end do

   return

end subroutine n2_ver_1

BLOCK DATA bn2f

   IMPLICIT REAL*8 (v)

   COMMON /n2_f/ V1n2f,V2n2f,DVn2f,NPTn2f,                           &
   &          xn2_272(228),xn2_228(228),a_h2o(228)

   DATA V1n2f,V2n2f,DVn2f,NPTn2f                                     &
   &     /1997.784896, 2901.576661, 3.981461525, 228/
   DATA xn2_272/                                                     &
   &      0.000E+00,                                                  &
   &      4.691E-11,  5.960E-11,  7.230E-11,  9.435E-11,  1.171E-10,  &
   &      1.472E-10,  1.874E-10,  2.276E-10,  2.960E-10,  3.671E-10,  &
   &      4.605E-10,  5.874E-10,  7.144E-10,  9.293E-10,  1.155E-09,  &
   &      1.447E-09,  1.847E-09,  2.247E-09,  2.919E-09,  3.635E-09,  &
   &      4.594E-09,  6.003E-09,  7.340E-09,  9.190E-09,  1.130E-08,  &
   &      1.370E-08,  1.650E-08,  1.960E-08,  2.310E-08,  2.710E-08,  &
   &      3.160E-08,  3.660E-08,  4.230E-08,  4.860E-08,  5.570E-08,  &
   &      6.350E-08,  7.230E-08,  8.200E-08,  9.270E-08,  1.050E-07,  &
   &      1.180E-07,  1.320E-07,  1.480E-07,  1.650E-07,  1.840E-07,  &
   &      2.040E-07,  2.270E-07,  2.510E-07,  2.770E-07,  3.060E-07,  &
   &      3.360E-07,  3.670E-07,  4.010E-07,  4.330E-07,  4.710E-07,  &
   &      5.050E-07,  5.450E-07,  5.790E-07,  6.200E-07,  6.540E-07,  &
   &      6.940E-07,  7.240E-07,  7.610E-07,  7.880E-07,  8.220E-07,  &
   &      8.440E-07,  8.720E-07,  8.930E-07,  9.190E-07,  9.370E-07,  &
   &      9.620E-07,  9.870E-07,  1.020E-06,  1.060E-06,  1.110E-06,  &
   &      1.180E-06,  1.280E-06,  1.400E-06,  1.570E-06,  1.750E-06,  &
   &      1.880E-06,  2.020E-06,  2.080E-06,  2.060E-06,  1.960E-06,  &
   &      1.860E-06,  1.710E-06,  1.570E-06,  1.490E-06,  1.440E-06,  &
   &      1.410E-06,  1.390E-06,  1.380E-06,  1.380E-06,  1.390E-06,  &
   &      1.390E-06,  1.410E-06,  1.420E-06,  1.430E-06,  1.420E-06,  &
   &      1.430E-06,  1.410E-06,  1.400E-06,  1.370E-06,  1.350E-06,  &
   &      1.310E-06,  1.270E-06,  1.220E-06,  1.170E-06,  1.120E-06,  &
   &      1.060E-06,  1.010E-06,  9.470E-07,  8.910E-07,  8.290E-07,  &
   &      7.740E-07,  7.160E-07,  6.620E-07,  6.090E-07,  5.600E-07,  &
   &      5.130E-07,  4.680E-07,  4.290E-07,  3.900E-07,  3.560E-07,  &
   &      3.240E-07,  2.950E-07,  2.680E-07,  2.440E-07,  2.230E-07,  &
   &      2.030E-07,  1.850E-07,  1.690E-07,  1.540E-07,  1.410E-07,  &
   &      1.290E-07,  1.180E-07,  1.080E-07,  9.950E-08,  9.100E-08,  &
   &      8.380E-08,  7.700E-08,  7.100E-08,  6.510E-08,  6.010E-08,  &
   &      5.550E-08,  5.110E-08,  4.710E-08,  4.340E-08,  3.980E-08,  &
   &      3.660E-08,  3.380E-08,  3.110E-08,  2.840E-08,  2.610E-08,  &
   &      2.390E-08,  2.210E-08,  2.010E-08,  1.830E-08,  1.710E-08,  &
   &      1.550E-08,  1.450E-08,  1.320E-08,  1.208E-08,  1.112E-08,  &
   &      1.015E-08,  9.339E-09,  8.597E-09,  7.873E-09,  7.247E-09,  &
   &      6.620E-09,  6.074E-09,  5.570E-09,  5.081E-09,  4.676E-09,  &
   &      4.272E-09,  3.919E-09,  3.595E-09,  3.279E-09,  3.019E-09,  &
   &      2.758E-09,  2.529E-09,  2.320E-09,  2.115E-09,  1.948E-09,  &
   &      1.780E-09,  1.632E-09,  1.497E-09,  1.365E-09,  1.257E-09,  &
   &      1.149E-09,  1.053E-09,  9.663E-10,  8.806E-10,  8.111E-10,  &
   &      7.416E-10,  6.795E-10,  6.237E-10,  5.682E-10,  5.233E-10,  &
   &      4.785E-10,  4.383E-10,  4.024E-10,  3.666E-10,  3.378E-10,  &
   &      3.090E-10,  2.829E-10,  2.598E-10,  2.366E-10,  2.180E-10,  &
   &      1.994E-10,  1.825E-10,  1.676E-10,  1.527E-10,  1.406E-10,  &
   &      1.287E-10,  1.178E-10,  1.082E-10,  9.859E-11,  9.076E-11,  &
   &      8.305E-11,  7.599E-11,  6.981E-11,  6.363E-11,  5.857E-11,  &
   &      5.362E-11,  0.000E+00/

   DATA xn2_228/                                                     &
   &      0.000E+00,                                                  &
   &      5.736E-11,  7.296E-11,  8.856E-11,  1.154E-10,  1.431E-10,  &
   &      1.799E-10,  2.291E-10,  2.783E-10,  3.623E-10,  4.497E-10,  &
   &      5.642E-10,  7.195E-10,  8.749E-10,  1.137E-09,  1.413E-09,  &
   &      1.769E-09,  2.259E-09,  2.749E-09,  3.568E-09,  4.440E-09,  &
   &      5.549E-09,  7.097E-09,  8.645E-09,  1.120E-08,  1.395E-08,  &
   &      1.650E-08,  1.880E-08,  2.130E-08,  2.400E-08,  2.690E-08,  &
   &      3.010E-08,  3.360E-08,  3.750E-08,  4.180E-08,  4.670E-08,  &
   &      5.210E-08,  5.830E-08,  6.520E-08,  7.290E-08,  8.170E-08,  &
   &      9.150E-08,  1.030E-07,  1.150E-07,  1.290E-07,  1.440E-07,  &
   &      1.610E-07,  1.800E-07,  2.020E-07,  2.250E-07,  2.510E-07,  &
   &      2.790E-07,  3.090E-07,  3.430E-07,  3.770E-07,  4.160E-07,  &
   &      4.540E-07,  4.990E-07,  5.370E-07,  5.850E-07,  6.250E-07,  &
   &      6.750E-07,  7.130E-07,  7.610E-07,  7.970E-07,  8.410E-07,  &
   &      8.720E-07,  9.100E-07,  9.380E-07,  9.720E-07,  9.940E-07,  &
   &      1.020E-06,  1.050E-06,  1.080E-06,  1.120E-06,  1.170E-06,  &
   &      1.240E-06,  1.340E-06,  1.470E-06,  1.660E-06,  1.870E-06,  &
   &      2.040E-06,  2.220E-06,  2.300E-06,  2.290E-06,  2.160E-06,  &
   &      2.050E-06,  1.870E-06,  1.710E-06,  1.620E-06,  1.580E-06,  &
   &      1.550E-06,  1.540E-06,  1.540E-06,  1.550E-06,  1.560E-06,  &
   &      1.570E-06,  1.590E-06,  1.590E-06,  1.600E-06,  1.580E-06,  &
   &      1.570E-06,  1.540E-06,  1.510E-06,  1.470E-06,  1.430E-06,  &
   &      1.370E-06,  1.310E-06,  1.250E-06,  1.180E-06,  1.110E-06,  &
   &      1.040E-06,  9.740E-07,  9.020E-07,  8.360E-07,  7.650E-07,  &
   &      7.050E-07,  6.430E-07,  5.860E-07,  5.320E-07,  4.820E-07,  &
   &      4.370E-07,  3.950E-07,  3.570E-07,  3.220E-07,  2.910E-07,  &
   &      2.630E-07,  2.390E-07,  2.160E-07,  1.960E-07,  1.780E-07,  &
   &      1.620E-07,  1.480E-07,  1.330E-07,  1.220E-07,  1.120E-07,  &
   &      1.020E-07,  9.280E-08,  8.420E-08,  7.700E-08,  6.990E-08,  &
   &      6.390E-08,  5.880E-08,  5.380E-08,  4.840E-08,  4.380E-08,  &
   &      4.020E-08,  3.690E-08,  3.290E-08,  3.050E-08,  2.720E-08,  &
   &      2.490E-08,  2.260E-08,  2.020E-08,  1.810E-08,  1.620E-08,  &
   &      1.500E-08,  1.359E-08,  1.232E-08,  1.111E-08,  1.011E-08,  &
   &      9.115E-09,  8.273E-09,  7.497E-09,  6.753E-09,  6.148E-09,  &
   &      5.543E-09,  5.029E-09,  4.558E-09,  4.105E-09,  3.738E-09,  &
   &      3.371E-09,  3.057E-09,  2.771E-09,  2.494E-09,  2.272E-09,  &
   &      2.049E-09,  1.858E-09,  1.685E-09,  1.516E-09,  1.381E-09,  &
   &      1.246E-09,  1.129E-09,  1.024E-09,  9.215E-10,  8.396E-10,  &
   &      7.578E-10,  6.865E-10,  6.228E-10,  5.601E-10,  5.105E-10,  &
   &      4.609E-10,  4.173E-10,  3.786E-10,  3.403E-10,  3.102E-10,  &
   &      2.802E-10,  2.536E-10,  2.302E-10,  2.069E-10,  1.886E-10,  &
   &      1.704E-10,  1.542E-10,  1.400E-10,  1.257E-10,  1.147E-10,  &
   &      1.036E-10,  9.371E-11,  8.509E-11,  7.647E-11,  6.970E-11,  &
   &      6.298E-11,  5.695E-11,  5.172E-11,  4.650E-11,  4.237E-11,  &
   &      3.829E-11,  3.462E-11,  3.145E-11,  2.828E-11,  2.576E-11,  &
   &      2.329E-11,  2.104E-11,  1.912E-11,  1.720E-11,  1.566E-11,  &
   &      1.416E-11,  0.000E+00/

   DATA a_h2o/                                                       &
   &      200.00,                                                     &
   &      199.95,  179.14,  158.32,  143.10,  128.29,                 &
   &      115.21,  104.51,  93.816,  85.716,  77.876,                 &
   &      70.838,  65.017,  59.196,  54.652,  50.272,                 &
   &      46.284,  42.951,  39.619,  36.945,  34.377,                 &
   &      32.009,  30.010,  28.012,  26.370,  24.797,                 &
   &      23.332,  22.085,  20.839,  19.794,  18.795,                 &
   &      17.855,  17.051,  16.247,  15.562,  14.909,                 &
   &      14.290,  13.756,  13.223,  12.763,  12.325,                 &
   &      11.906,  11.544,  11.182,  10.870,  10.575,                 &
   &      10.292,  10.045,  9.7985,  9.5833,  9.3802,                 &
   &      9.1847,  9.0175,  8.8503,  8.7025,  8.5632,                 &
   &      8.4290,  8.3176,  8.2061,  8.1081,  8.0165,                 &
   &      7.9282,  7.8565,  7.7848,  7.7236,  7.6678,                 &
   &      7.6138,  7.5700,  7.5262,  7.4899,  7.4581,                 &
   &      7.4271,  7.4032,  7.3794,  7.3603,  7.3444,                 &
   &      7.3292,  7.3212,  7.3133,  7.3076,  7.3036,                 &
   &      7.2997,  7.2957,  7.2917,  7.2900,  7.2900,                 &
   &      7.2900,  7.2900,  7.2900,  7.2921,  7.2961,                 &
   &      7.3000,  7.3000,  7.3000,  7.3040,  7.3120,                 &
   &      7.3200,  7.3359,  7.3518,  7.3736,  7.4015,                 &
   &      7.4293,  7.4728,  7.5166,  7.5697,  7.6334,                 &
   &      7.6971,  7.7874,  7.8790,  7.9846,  8.1081,                 &
   &      8.2315,  8.3920,  8.5552,  8.7419,  8.9609,                 &
   &      9.1798,  9.4567,  9.7394,  10.059,  10.438,                 &
   &      10.816,  11.300,  11.798,  12.358,  13.023,                 &
   &      13.687,  14.541,  15.425,  16.402,  17.553,                 &
   &      18.704,  20.098,  21.539,  23.061,  24.749,                 &
   &      26.437,  28.236,  30.059,  31.710,  32.964,                 &
   &      34.218,  35.570,  36.948,  38.366,  39.886,                 &
   &      41.407,  43.045,  44.717,  46.432,  48.271,                 &
   &      50.111,  52.090,  54.116,  56.190,  58.420,                 &
   &      60.649,  63.040,  65.493,  67.997,  70.697,                 &
   &      73.396,  76.287,  79.262,  82.291,  85.560,                 &
   &      88.828,  92.320,  95.920,  99.578,  103.54,                 &
   &      107.50,  111.72,  116.08,  120.50,  125.30,                 &
   &      130.09,  135.20,  140.48,  145.83,  151.64,                 &
   &      157.45,  163.62,  170.01,  176.47,  183.51,                 &
   &      190.55,  198.01,  205.75,  213.56,  222.09,                 &
   &      230.61,  239.62,  249.01,  258.43,  268.76,                 &
   &      279.09,  289.98,  301.35,  312.75,  325.26,                 &
   &      337.77,  350.92,  364.69,  378.46,  393.61,                 &
   &      408.77,  424.67,  441.35,  458.03,  476.35,                 &
   &      494.70,  513.92,  534.12,  554.33,  576.47,                 &
   &      598.70,  621.92,  646.40,  670.87,  697.63,                 &
   &      724.56,  752.63,  782.27,  811.91,  844.26,                 &
   &      876.88,  876.88/

end block data bn2f
!
!     --------------------------------------------------------------
!
subroutine n2_overtone1 (v1c,v2c,dvc,nptc,c,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8 (v)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
!
   COMMON /n2_f1/ V1S,V2S,DVS,NPTS,xn2(191)
!
   dimension c(*)
!
!     Nitrogen Collision Induced First Overtone

!     Shapiro and Gush (1966) modified by Mlawer and Gombos (2015).

!     The absorption coefficients are for pure nitrogen (absorber and
!     broadener.
!
   DVC = DVS
   v1ss = v1s
   v2ss = v2s
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   do 10 j=1,nptc
      i = i1+(j-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      VJ = V1C+DVC* REAL(J-1)
!
      c(j) = xn2(i)

!     the radiation field is removed with 1/vj

      c(j) = c(j)/vj
!
10 end do

   return

end subroutine n2_overtone1

BLOCK DATA bn2f1

   IMPLICIT REAL*8 (v)

   COMMON /n2_f1/ V1n2f,V2n2f,DVn2f,NPTn2f,xn2(191)

   DATA V1n2f,V2n2f,DVn2f,NPTn2f                                     &
   &     /4340.0, 4910.0, 3.0, 191/
   DATA xn2/                                                    &
   &      0.000E+00, 3.709E-11, 7.418E-11, 1.113E-10, 1.484E-10, &
   &      1.843E-10, 2.163E-10, 2.482E-10, 2.802E-10, 3.122E-10, &
   &      3.442E-10, 3.640E-10, 3.776E-10, 3.912E-10, 4.048E-10, &
   &      4.183E-10, 4.334E-10, 4.626E-10, 4.918E-10, 5.210E-10, &
   &      5.357E-10, 5.411E-10, 5.465E-10, 5.520E-10, 5.593E-10, &
   &      5.853E-10, 6.114E-10, 6.375E-10, 6.635E-10, 6.855E-10, &
   &      6.926E-10, 6.997E-10, 7.069E-10, 7.140E-10, 7.211E-10, &
   &      7.283E-10, 7.380E-10, 7.551E-10, 7.722E-10, 7.893E-10, &
   &      8.064E-10, 8.235E-10, 8.419E-10, 8.627E-10, 8.835E-10, &
   &      9.043E-10, 9.251E-10, 9.779E-10, 1.066E-09, 1.154E-09, &
   &      1.242E-09, 1.344E-09, 1.446E-09, 1.549E-09, 1.653E-09, &
   &      1.759E-09, 1.865E-09, 1.977E-09, 2.103E-09, 2.228E-09, &
   &      2.348E-09, 2.467E-09, 2.586E-09, 2.705E-09, 2.824E-09, &
   &      2.944E-09, 3.066E-09, 3.188E-09, 3.309E-09, 3.426E-09, &
   &      3.543E-09, 3.660E-09, 3.813E-09, 3.976E-09, 4.135E-09, &
   &      4.309E-09, 4.499E-09, 4.700E-09, 4.905E-09, 5.105E-09, &
   &      5.332E-09, 5.575E-09, 5.856E-09, 6.175E-09, 6.421E-09, &
   &      6.640E-09, 7.086E-09, 7.508E-09, 7.906E-09, 8.304E-09, &
   &      8.930E-09, 9.480E-09, 9.921E-09, 1.051E-08, 1.070E-08, &
   &      1.090E-08, 1.090E-08, 1.090E-08, 1.090E-08, 1.070E-08, &
   &      1.050E-08, 1.030E-08, 1.010E-08, 9.940E-09, 9.760E-09, &
   &      9.580E-09, 9.400E-09, 9.220E-09, 9.040E-09, 8.860E-09, &
   &      8.680E-09, 8.510E-09, 8.370E-09, 8.250E-09, 8.150E-09, &
   &      8.070E-09, 8.010E-09, 7.950E-09, 7.890E-09, 7.830E-09, &
   &      7.759E-09, 7.553E-09, 7.347E-09, 7.141E-09, 6.935E-09, &
   &      6.729E-09, 6.523E-09, 6.317E-09, 6.111E-09, 5.905E-09, &
   &      5.699E-09, 5.493E-09, 5.287E-09, 5.081E-09, 4.876E-09, &
   &      4.670E-09, 4.464E-09, 4.258E-09, 4.052E-09, 3.846E-09, &
   &      3.640E-09, 3.450E-09, 3.268E-09, 3.092E-09, 2.917E-09, &
   &      2.741E-09, 2.566E-09, 2.392E-09, 2.219E-09, 2.076E-09, &
   &      1.959E-09, 1.841E-09, 1.723E-09, 1.617E-09, 1.527E-09, &
   &      1.437E-09, 1.346E-09, 1.256E-09, 1.172E-09, 1.093E-09, &
   &      1.013E-09, 9.335E-10, 8.539E-10, 7.979E-10, 7.514E-10, &
   &      7.050E-10, 6.586E-10, 6.121E-10, 5.687E-10, 5.413E-10, &
   &      5.138E-10, 4.864E-10, 4.589E-10, 4.315E-10, 4.040E-10, &
   &      3.770E-10, 3.504E-10, 3.238E-10, 2.972E-10, 2.706E-10, &
   &      2.409E-10, 2.099E-10, 1.788E-10, 1.478E-10, 1.225E-10, &
   &      1.021E-10, 8.165E-11, 6.123E-11, 4.082E-11, 2.041E-11, &
   &      0.000E+00/
end block data bn2f1
!
!     --------------------------------------------------------------
!
SUBROUTINE XO3CHP (V1C,V2C,DVC,NPTC,C0,C1,C2,v1ss,v2ss)

   Use lblparams, ONLY: n_absrb
!
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   COMMON /O3CHAP/ V1S,V2S,DVS,NPTS,X(3150),Y(3150),Z(3150)
   DIMENSION C0(*),C1(*),C2(*)
!
   DVC = DVS
   v1ss = v1s
   v2ss = v2s
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      IF ((I.LT.1).OR.(I.GT.NPTS)) THEN
         C0(J) = 0.
         C1(J)=0.
         C2(J)=0.
      ELSE
!
!            Remove radiation field from diffuse ozone
!
         VJ = V1C+DVC* REAL(J-1)
         C0(J)=X(I)/VJ
         C1(J)=Y(I)/VJ
         C2(J)=Z(I)/VJ
      ENDIF
10 END DO
!
   RETURN
!
end subroutine XO3CHP
!
!     --------------------------------------------------------------
BLOCK DATA O3CH
!
!     CHAPPUIS AND WULF BAND
!
!     BEGINNING AND ENDING FREQUENCIES FROM DATA (CM-1):
!
!                        9170.0 24565.0
!
!     Added points at beginning and end (X,Y,Z(1:50) and
!     X,Y,Z(3130:3150)).  Zeroed values of Y,Z(1:789) to eliminate
!     ringing from interpolations done in MODTRAN.  Changed coefficients
!     X(32:50,3130:3150) Y(821:841,3130:3150), & Z(821:841,3130:3150) to
!     smooth coefficients to zero.
!     Smoothing coefficient frequencies (cm-1):
!
!             9075.0 -  9165.0  and 24570.0 - 24665.0 for X
!            13020.0 - 13120.0  and 24570.0 - 24665.0 for Y
!            13020.0 - 13120.0  and 24570.0 - 24665.0 for Z
!
!
!
!     CROSS-SECTIONS IN CM^2 TIMES 1.0E20
!     FORMULA FOR CROSS SECTION:  X+Y*DT+Z*DT*DT, DT=T-273.15
!     THE OUTPUT OF THIS ROUTINE IS C0=X, CT1=Y AND CT2=Z.
!
   IMPLICIT REAL*8           (V)
   COMMON /O3CHAP/VBEG,VEND,DVINCR,NMAX,X(3150),Y(3150),Z(3150)
!
   DATA VBEG, VEND, DVINCR, NMAX /8920.0, 24665.0, 5.0, 3150/
   DATA (X(I),I=    1,  50)/                                         &
   &      0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  &
   &      0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  &
   &      0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  &
   &      0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  &
   &      0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  &
   &      0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  0.00000  ,  &
   &      0.00000  ,  0.075E-05,  0.150E-05,  0.225E-05,  0.300E-05,  &
   &      0.400E-05,  0.500E-05,  0.600E-05,  0.700E-05,  0.850E-05,  &
   &      1.000E-05,  1.200E-05,  1.430E-05,  1.680E-05,  1.980E-05,  &
   &      2.280E-05,  2.630E-05,  2.980E-05,  3.376E-05,  3.826E-05/
   DATA (X(I),I= 51, 100)/                                           &
   &      4.276E-05,  4.775E-05,  5.825E-05,  6.908E-05,  7.299E-05,  &
   &      7.116E-05,  7.388E-05,  7.965E-05,  7.689E-05,  6.900E-05,  &
   &      7.008E-05,  6.945E-05,  7.083E-05,  7.053E-05,  6.908E-05,  &
   &      6.923E-05,  6.770E-05,  7.146E-05,  7.749E-05,  8.464E-05,  &
   &      8.441E-05,  8.754E-05,  8.795E-05,  9.971E-05,  9.632E-05,  &
   &      9.539E-05,  1.037E-04,  1.085E-04,  1.058E-04,  1.077E-04,  &
   &      1.121E-04,  1.193E-04,  1.292E-04,  1.364E-04,  1.526E-04,  &
   &      1.658E-04,  1.808E-04,  1.861E-04,  1.786E-04,  1.804E-04,  &
   &      1.885E-04,  1.972E-04,  2.218E-04,  2.408E-04,  2.317E-04,  &
   &      2.098E-04,  1.938E-04,  1.851E-04,  1.896E-04,  1.875E-04/
   DATA (X(I),I=  101,  150)/                                        &
   &      1.708E-04,  1.710E-04,  1.796E-04,  1.865E-04,  1.943E-04,  &
   &      1.881E-04,  1.885E-04,  2.136E-04,  2.255E-04,  2.267E-04,  &
   &      2.234E-04,  2.418E-04,  2.695E-04,  2.710E-04,  2.738E-04,  &
   &      3.066E-04,  3.269E-04,  3.465E-04,  3.986E-04,  4.410E-04,  &
   &      4.719E-04,  5.051E-04,  5.211E-04,  5.132E-04,  5.125E-04,  &
   &      5.159E-04,  5.549E-04,  6.562E-04,  7.168E-04,  6.502E-04,  &
   &      5.140E-04,  4.161E-04,  3.620E-04,  3.264E-04,  3.004E-04,  &
   &      2.815E-04,  2.650E-04,  2.527E-04,  2.424E-04,  2.292E-04,  &
   &      2.155E-04,  2.072E-04,  1.992E-04,  1.943E-04,  1.914E-04,  &
   &      1.855E-04,  1.813E-04,  1.724E-04,  1.687E-04,  1.676E-04/
   DATA (X(I),I= 151, 200)/                                          &
   &      1.601E-04,  1.503E-04,  1.518E-04,  1.436E-04,  1.455E-04,  &
   &      1.448E-04,  1.410E-04,  1.406E-04,  1.425E-04,  1.407E-04,  &
   &      1.405E-04,  1.436E-04,  1.369E-04,  1.355E-04,  1.331E-04,  &
   &      1.328E-04,  1.350E-04,  1.394E-04,  1.372E-04,  1.444E-04,  &
   &      1.490E-04,  1.455E-04,  1.460E-04,  1.523E-04,  1.559E-04,  &
   &      1.654E-04,  1.766E-04,  1.843E-04,  1.911E-04,  1.881E-04,  &
   &      1.894E-04,  1.927E-04,  2.043E-04,  2.106E-04,  2.215E-04,  &
   &      2.268E-04,  2.249E-04,  2.230E-04,  2.302E-04,  2.408E-04,  &
   &      2.518E-04,  2.625E-04,  2.753E-04,  2.788E-04,  2.701E-04,  &
   &      2.746E-04,  2.935E-04,  3.173E-04,  3.457E-04,  3.452E-04/
   DATA (X(I),I=  201,  250)/                                        &
   &      3.329E-04,  3.443E-04,  3.706E-04,  4.079E-04,  4.403E-04,  &
   &      4.343E-04,  4.172E-04,  4.448E-04,  5.132E-04,  5.635E-04,  &
   &      5.590E-04,  5.419E-04,  6.007E-04,  6.912E-04,  7.258E-04,  &
   &      7.146E-04,  7.529E-04,  8.706E-04,  9.465E-04,  9.923E-04,  &
   &      1.134E-03,  1.286E-03,  1.351E-03,  1.485E-03,  1.709E-03,  &
   &      1.897E-03,  2.086E-03,  2.186E-03,  2.195E-03,  2.185E-03,  &
   &      2.199E-03,  2.336E-03,  2.666E-03,  3.076E-03,  3.075E-03,  &
   &      2.543E-03,  1.920E-03,  1.498E-03,  1.283E-03,  1.165E-03,  &
   &      1.070E-03,  9.833E-04,  9.018E-04,  8.207E-04,  7.451E-04,  &
   &      6.811E-04,  6.178E-04,  5.661E-04,  5.199E-04,  4.868E-04/
   DATA (X(I),I= 251, 300)/                                          &
   &      4.541E-04,  4.291E-04,  4.135E-04,  3.990E-04,  3.878E-04,  &
   &      3.815E-04,  3.722E-04,  3.691E-04,  3.726E-04,  3.711E-04,  &
   &      3.744E-04,  3.778E-04,  3.808E-04,  3.826E-04,  3.852E-04,  &
   &      3.919E-04,  3.975E-04,  3.990E-04,  4.053E-04,  4.176E-04,  &
   &      4.232E-04,  4.291E-04,  4.414E-04,  4.541E-04,  4.723E-04,  &
   &      4.887E-04,  5.058E-04,  5.226E-04,  5.501E-04,  5.836E-04,  &
   &      6.059E-04,  6.238E-04,  6.469E-04,  6.711E-04,  7.046E-04,  &
   &      7.448E-04,  7.794E-04,  8.054E-04,  8.222E-04,  8.371E-04,  &
   &      8.538E-04,  8.612E-04,  8.698E-04,  8.914E-04,  9.122E-04,  &
   &      9.305E-04,  9.562E-04,  9.844E-04,  1.018E-03,  1.053E-03/
   DATA (X(I),I=  301,  350)/                                        &
   &      1.091E-03,  1.136E-03,  1.187E-03,  1.233E-03,  1.289E-03,  &
   &      1.336E-03,  1.372E-03,  1.405E-03,  1.435E-03,  1.470E-03,  &
   &      1.504E-03,  1.517E-03,  1.511E-03,  1.541E-03,  1.619E-03,  &
   &      1.728E-03,  1.848E-03,  1.955E-03,  2.044E-03,  2.128E-03,  &
   &      2.254E-03,  2.396E-03,  2.527E-03,  2.660E-03,  2.832E-03,  &
   &      3.010E-03,  3.182E-03,  3.340E-03,  3.504E-03,  3.673E-03,  &
   &      3.822E-03,  3.923E-03,  3.997E-03,  4.042E-03,  4.061E-03,  &
   &      4.035E-03,  3.979E-03,  3.901E-03,  3.785E-03,  3.642E-03,  &
   &      3.494E-03,  3.339E-03,  3.173E-03,  3.004E-03,  2.849E-03,  &
   &      2.703E-03,  2.556E-03,  2.432E-03,  2.310E-03,  2.191E-03/
   DATA (X(I),I= 351, 400)/                                          &
   &      2.076E-03,  1.969E-03,  1.883E-03,  1.818E-03,  1.753E-03,  &
   &      1.705E-03,  1.672E-03,  1.643E-03,  1.617E-03,  1.616E-03,  &
   &      1.629E-03,  1.648E-03,  1.662E-03,  1.667E-03,  1.669E-03,  &
   &      1.664E-03,  1.655E-03,  1.645E-03,  1.643E-03,  1.642E-03,  &
   &      1.632E-03,  1.629E-03,  1.632E-03,  1.638E-03,  1.644E-03,  &
   &      1.647E-03,  1.646E-03,  1.642E-03,  1.638E-03,  1.632E-03,  &
   &      1.628E-03,  1.626E-03,  1.628E-03,  1.635E-03,  1.642E-03,  &
   &      1.649E-03,  1.653E-03,  1.656E-03,  1.660E-03,  1.669E-03,  &
   &      1.685E-03,  1.705E-03,  1.730E-03,  1.755E-03,  1.779E-03,  &
   &      1.804E-03,  1.830E-03,  1.861E-03,  1.896E-03,  1.931E-03/
   DATA (X(I),I=  401,  450)/                                        &
   &      1.962E-03,  1.991E-03,  2.024E-03,  2.068E-03,  2.131E-03,  &
   &      2.207E-03,  2.285E-03,  2.357E-03,  2.423E-03,  2.490E-03,  &
   &      2.564E-03,  2.649E-03,  2.743E-03,  2.842E-03,  2.943E-03,  &
   &      3.044E-03,  3.146E-03,  3.248E-03,  3.350E-03,  3.452E-03,  &
   &      3.555E-03,  3.664E-03,  3.785E-03,  3.927E-03,  4.083E-03,  &
   &      4.250E-03,  4.418E-03,  4.570E-03,  4.708E-03,  4.835E-03,  &
   &      4.961E-03,  5.088E-03,  5.218E-03,  5.348E-03,  5.471E-03,  &
   &      5.594E-03,  5.713E-03,  5.828E-03,  5.933E-03,  6.026E-03,  &
   &      6.100E-03,  6.152E-03,  6.186E-03,  6.193E-03,  6.182E-03,  &
   &      6.149E-03,  6.093E-03,  6.011E-03,  5.914E-03,  5.799E-03/
   DATA (X(I),I= 451, 500)/                                          &
   &      5.676E-03,  5.553E-03,  5.438E-03,  5.330E-03,  5.233E-03,  &
   &      5.151E-03,  5.080E-03,  5.025E-03,  4.987E-03,  4.972E-03,  &
   &      4.976E-03,  4.991E-03,  5.013E-03,  5.032E-03,  5.043E-03,  &
   &      5.043E-03,  5.032E-03,  5.010E-03,  4.980E-03,  4.950E-03,  &
   &      4.913E-03,  4.879E-03,  4.838E-03,  4.786E-03,  4.723E-03,  &
   &      4.652E-03,  4.578E-03,  4.503E-03,  4.433E-03,  4.366E-03,  &
   &      4.306E-03,  4.247E-03,  4.191E-03,  4.135E-03,  4.083E-03,  &
   &      4.035E-03,  3.997E-03,  3.968E-03,  3.945E-03,  3.923E-03,  &
   &      3.904E-03,  3.886E-03,  3.867E-03,  3.856E-03,  3.848E-03,  &
   &      3.845E-03,  3.848E-03,  3.860E-03,  3.878E-03,  3.897E-03/
   DATA (X(I),I=  501,  550)/                                        &
   &      3.915E-03,  3.941E-03,  3.971E-03,  4.008E-03,  4.057E-03,  &
   &      4.113E-03,  4.176E-03,  4.243E-03,  4.325E-03,  4.418E-03,  &
   &      4.518E-03,  4.626E-03,  4.723E-03,  4.812E-03,  4.891E-03,  &
   &      4.972E-03,  5.058E-03,  5.155E-03,  5.263E-03,  5.382E-03,  &
   &      5.516E-03,  5.657E-03,  5.810E-03,  5.974E-03,  6.145E-03,  &
   &      6.331E-03,  6.532E-03,  6.752E-03,  6.993E-03,  7.247E-03,  &
   &      7.507E-03,  7.768E-03,  8.036E-03,  8.304E-03,  8.579E-03,  &
   &      8.862E-03,  9.148E-03,  9.442E-03,  9.744E-03,  1.006E-02,  &
   &      1.038E-02,  1.071E-02,  1.104E-02,  1.137E-02,  1.168E-02,  &
   &      1.195E-02,  1.220E-02,  1.242E-02,  1.264E-02,  1.283E-02/
   DATA (X(I),I= 551, 600)/                                          &
   &      1.303E-02,  1.322E-02,  1.339E-02,  1.356E-02,  1.371E-02,  &
   &      1.385E-02,  1.398E-02,  1.408E-02,  1.415E-02,  1.417E-02,  &
   &      1.415E-02,  1.408E-02,  1.395E-02,  1.376E-02,  1.353E-02,  &
   &      1.326E-02,  1.295E-02,  1.262E-02,  1.228E-02,  1.194E-02,  &
   &      1.161E-02,  1.128E-02,  1.097E-02,  1.067E-02,  1.038E-02,  &
   &      1.011E-02,  9.859E-03,  9.625E-03,  9.409E-03,  9.208E-03,  &
   &      9.022E-03,  8.843E-03,  8.668E-03,  8.505E-03,  8.348E-03,  &
   &      8.207E-03,  8.088E-03,  7.987E-03,  7.909E-03,  7.842E-03,  &
   &      7.782E-03,  7.727E-03,  7.675E-03,  7.619E-03,  7.570E-03,  &
   &      7.526E-03,  7.488E-03,  7.459E-03,  7.440E-03,  7.429E-03/
   DATA (X(I),I=  601,  650)/                                        &
   &      7.429E-03,  7.429E-03,  7.440E-03,  7.455E-03,  7.474E-03,  &
   &      7.500E-03,  7.529E-03,  7.563E-03,  7.593E-03,  7.622E-03,  &
   &      7.649E-03,  7.675E-03,  7.715E-03,  7.771E-03,  7.846E-03,  &
   &      7.939E-03,  8.039E-03,  8.147E-03,  8.255E-03,  8.367E-03,  &
   &      8.482E-03,  8.605E-03,  8.746E-03,  8.903E-03,  9.078E-03,  &
   &      9.271E-03,  9.472E-03,  9.677E-03,  9.889E-03,  1.011E-02,  &
   &      1.034E-02,  1.059E-02,  1.085E-02,  1.113E-02,  1.143E-02,  &
   &      1.174E-02,  1.207E-02,  1.242E-02,  1.277E-02,  1.313E-02,  &
   &      1.350E-02,  1.388E-02,  1.425E-02,  1.464E-02,  1.503E-02,  &
   &      1.544E-02,  1.586E-02,  1.628E-02,  1.670E-02,  1.713E-02/
   DATA (X(I),I= 651, 700)/                                          &
   &      1.755E-02,  1.796E-02,  1.837E-02,  1.875E-02,  1.911E-02,  &
   &      1.945E-02,  1.975E-02,  2.002E-02,  2.028E-02,  2.050E-02,  &
   &      2.070E-02,  2.089E-02,  2.104E-02,  2.117E-02,  2.126E-02,  &
   &      2.132E-02,  2.135E-02,  2.135E-02,  2.130E-02,  2.123E-02,  &
   &      2.114E-02,  2.101E-02,  2.087E-02,  2.072E-02,  2.053E-02,  &
   &      2.032E-02,  2.010E-02,  1.986E-02,  1.963E-02,  1.939E-02,  &
   &      1.915E-02,  1.891E-02,  1.868E-02,  1.845E-02,  1.821E-02,  &
   &      1.798E-02,  1.773E-02,  1.746E-02,  1.719E-02,  1.692E-02,  &
   &      1.666E-02,  1.643E-02,  1.621E-02,  1.598E-02,  1.576E-02,  &
   &      1.558E-02,  1.542E-02,  1.529E-02,  1.519E-02,  1.509E-02/
   DATA (X(I),I=  701,  750)/                                        &
   &      1.501E-02,  1.493E-02,  1.484E-02,  1.477E-02,  1.473E-02,  &
   &      1.471E-02,  1.469E-02,  1.468E-02,  1.468E-02,  1.470E-02,  &
   &      1.473E-02,  1.475E-02,  1.476E-02,  1.477E-02,  1.480E-02,  &
   &      1.484E-02,  1.489E-02,  1.497E-02,  1.507E-02,  1.520E-02,  &
   &      1.533E-02,  1.543E-02,  1.551E-02,  1.557E-02,  1.563E-02,  &
   &      1.569E-02,  1.575E-02,  1.581E-02,  1.591E-02,  1.602E-02,  &
   &      1.614E-02,  1.625E-02,  1.637E-02,  1.654E-02,  1.676E-02,  &
   &      1.698E-02,  1.719E-02,  1.740E-02,  1.762E-02,  1.784E-02,  &
   &      1.805E-02,  1.827E-02,  1.852E-02,  1.878E-02,  1.906E-02,  &
   &      1.935E-02,  1.962E-02,  1.985E-02,  2.005E-02,  2.024E-02/
   DATA (X(I),I= 751, 800)/                                          &
   &      2.047E-02,  2.073E-02,  2.106E-02,  2.137E-02,  2.167E-02,  &
   &      2.194E-02,  2.223E-02,  2.256E-02,  2.287E-02,  2.316E-02,  &
   &      2.344E-02,  2.373E-02,  2.405E-02,  2.445E-02,  2.490E-02,  &
   &      2.542E-02,  2.592E-02,  2.640E-02,  2.684E-02,  2.729E-02,  &
   &      2.778E-02,  2.828E-02,  2.876E-02,  2.914E-02,  2.948E-02,  &
   &      2.979E-02,  3.010E-02,  3.039E-02,  3.066E-02,  3.091E-02,  &
   &      3.116E-02,  3.138E-02,  3.153E-02,  3.155E-02,  3.152E-02,  &
   &      3.146E-02,  3.144E-02,  3.138E-02,  3.126E-02,  3.110E-02,  &
   &      3.092E-02,  3.073E-02,  3.054E-02,  3.033E-02,  3.008E-02,  &
   &      2.980E-02,  2.947E-02,  2.910E-02,  2.870E-02,  2.832E-02/
   DATA (X(I),I=  801,  850)/                                        &
   &      2.795E-02,  2.765E-02,  2.735E-02,  2.706E-02,  2.680E-02,  &
   &      2.656E-02,  2.637E-02,  2.620E-02,  2.604E-02,  2.587E-02,  &
   &      2.570E-02,  2.551E-02,  2.533E-02,  2.520E-02,  2.514E-02,  &
   &      2.513E-02,  2.513E-02,  2.512E-02,  2.509E-02,  2.506E-02,  &
   &      2.504E-02,  2.501E-02,  2.498E-02,  2.494E-02,  2.493E-02,  &
   &      2.497E-02,  2.508E-02,  2.522E-02,  2.535E-02,  2.545E-02,  &
   &      2.550E-02,  2.558E-02,  2.568E-02,  2.578E-02,  2.587E-02,  &
   &      2.592E-02,  2.598E-02,  2.605E-02,  2.619E-02,  2.631E-02,  &
   &      2.621E-02,  2.617E-02,  2.629E-02,  2.642E-02,  2.654E-02,  &
   &      2.669E-02,  2.685E-02,  2.700E-02,  2.716E-02,  2.734E-02/
   DATA (X(I),I= 851, 900)/                                          &
   &      2.752E-02,  2.772E-02,  2.792E-02,  2.813E-02,  2.834E-02,  &
   &      2.858E-02,  2.885E-02,  2.913E-02,  2.941E-02,  2.973E-02,  &
   &      3.005E-02,  3.038E-02,  3.075E-02,  3.117E-02,  3.159E-02,  &
   &      3.202E-02,  3.246E-02,  3.290E-02,  3.335E-02,  3.384E-02,  &
   &      3.438E-02,  3.493E-02,  3.547E-02,  3.603E-02,  3.660E-02,  &
   &      3.718E-02,  3.772E-02,  3.826E-02,  3.879E-02,  3.931E-02,  &
   &      3.987E-02,  4.042E-02,  4.098E-02,  4.151E-02,  4.199E-02,  &
   &      4.243E-02,  4.287E-02,  4.316E-02,  4.344E-02,  4.369E-02,  &
   &      4.392E-02,  4.405E-02,  4.417E-02,  4.429E-02,  4.436E-02,  &
   &      4.436E-02,  4.438E-02,  4.437E-02,  4.427E-02,  4.416E-02/
   DATA (X(I),I=  901,  950)/                                        &
   &      4.405E-02,  4.394E-02,  4.383E-02,  4.372E-02,  4.359E-02,  &
   &      4.344E-02,  4.329E-02,  4.312E-02,  4.299E-02,  4.289E-02,  &
   &      4.278E-02,  4.269E-02,  4.258E-02,  4.242E-02,  4.227E-02,  &
   &      4.213E-02,  4.202E-02,  4.194E-02,  4.183E-02,  4.178E-02,  &
   &      4.179E-02,  4.177E-02,  4.175E-02,  4.174E-02,  4.174E-02,  &
   &      4.175E-02,  4.177E-02,  4.183E-02,  4.191E-02,  4.199E-02,  &
   &      4.207E-02,  4.214E-02,  4.219E-02,  4.226E-02,  4.232E-02,  &
   &      4.239E-02,  4.247E-02,  4.254E-02,  4.261E-02,  4.270E-02,  &
   &      4.279E-02,  4.287E-02,  4.298E-02,  4.311E-02,  4.323E-02,  &
   &      4.338E-02,  4.357E-02,  4.375E-02,  4.393E-02,  4.413E-02/
   DATA (X(I),I= 951, 1000)/                                         &
   &      4.433E-02,  4.454E-02,  4.475E-02,  4.499E-02,  4.525E-02,  &
   &      4.551E-02,  4.579E-02,  4.613E-02,  4.645E-02,  4.678E-02,  &
   &      4.712E-02,  4.749E-02,  4.785E-02,  4.820E-02,  4.856E-02,  &
   &      4.894E-02,  4.928E-02,  4.963E-02,  4.991E-02,  5.016E-02,  &
   &      5.042E-02,  5.072E-02,  5.109E-02,  5.144E-02,  5.183E-02,  &
   &      5.218E-02,  5.251E-02,  5.282E-02,  5.315E-02,  5.351E-02,  &
   &      5.391E-02,  5.430E-02,  5.471E-02,  5.510E-02,  5.548E-02,  &
   &      5.588E-02,  5.628E-02,  5.673E-02,  5.722E-02,  5.771E-02,  &
   &      5.821E-02,  5.874E-02,  5.927E-02,  5.980E-02,  6.034E-02,  &
   &      6.088E-02,  6.144E-02,  6.197E-02,  6.250E-02,  6.303E-02/
   DATA (X(I),I= 1001, 1050)/                                        &
   &      6.352E-02,  6.404E-02,  6.452E-02,  6.493E-02,  6.537E-02,  &
   &      6.578E-02,  6.617E-02,  6.653E-02,  6.688E-02,  6.722E-02,  &
   &      6.747E-02,  6.768E-02,  6.788E-02,  6.808E-02,  6.827E-02,  &
   &      6.842E-02,  6.859E-02,  6.875E-02,  6.884E-02,  6.889E-02,  &
   &      6.896E-02,  6.900E-02,  6.915E-02,  6.927E-02,  6.942E-02,  &
   &      6.956E-02,  6.972E-02,  6.990E-02,  7.009E-02,  7.024E-02,  &
   &      7.043E-02,  7.062E-02,  7.081E-02,  7.102E-02,  7.127E-02,  &
   &      7.151E-02,  7.175E-02,  7.199E-02,  7.225E-02,  7.248E-02,  &
   &      7.274E-02,  7.300E-02,  7.325E-02,  7.351E-02,  7.375E-02,  &
   &      7.402E-02,  7.429E-02,  7.458E-02,  7.488E-02,  7.514E-02/
   DATA (X(I),I= 1051, 1100)/                                        &
   &      7.546E-02,  7.575E-02,  7.605E-02,  7.634E-02,  7.667E-02,  &
   &      7.698E-02,  7.732E-02,  7.767E-02,  7.803E-02,  7.841E-02,  &
   &      7.879E-02,  7.916E-02,  7.959E-02,  8.001E-02,  8.042E-02,  &
   &      8.082E-02,  8.118E-02,  8.152E-02,  8.188E-02,  8.224E-02,  &
   &      8.270E-02,  8.318E-02,  8.367E-02,  8.415E-02,  8.467E-02,  &
   &      8.518E-02,  8.569E-02,  8.621E-02,  8.673E-02,  8.725E-02,  &
   &      8.779E-02,  8.831E-02,  8.887E-02,  8.945E-02,  9.003E-02,  &
   &      9.060E-02,  9.123E-02,  9.187E-02,  9.254E-02,  9.317E-02,  &
   &      9.382E-02,  9.444E-02,  9.506E-02,  9.570E-02,  9.634E-02,  &
   &      9.702E-02,  9.769E-02,  9.838E-02,  9.904E-02,  9.968E-02/
   DATA (X(I),I= 1101, 1150)/                                        &
   &      1.003E-01,  1.010E-01,  1.016E-01,  1.022E-01,  1.028E-01,  &
   &      1.033E-01,  1.039E-01,  1.046E-01,  1.053E-01,  1.060E-01,  &
   &      1.067E-01,  1.075E-01,  1.082E-01,  1.089E-01,  1.096E-01,  &
   &      1.103E-01,  1.110E-01,  1.117E-01,  1.125E-01,  1.132E-01,  &
   &      1.139E-01,  1.147E-01,  1.154E-01,  1.162E-01,  1.169E-01,  &
   &      1.177E-01,  1.184E-01,  1.191E-01,  1.197E-01,  1.203E-01,  &
   &      1.209E-01,  1.215E-01,  1.221E-01,  1.227E-01,  1.232E-01,  &
   &      1.238E-01,  1.244E-01,  1.249E-01,  1.254E-01,  1.259E-01,  &
   &      1.263E-01,  1.268E-01,  1.273E-01,  1.278E-01,  1.283E-01,  &
   &      1.288E-01,  1.292E-01,  1.296E-01,  1.300E-01,  1.305E-01/
   DATA (X(I),I= 1151, 1200)/                                        &
   &      1.309E-01,  1.314E-01,  1.319E-01,  1.324E-01,  1.329E-01,  &
   &      1.335E-01,  1.340E-01,  1.345E-01,  1.351E-01,  1.356E-01,  &
   &      1.362E-01,  1.368E-01,  1.374E-01,  1.380E-01,  1.386E-01,  &
   &      1.392E-01,  1.399E-01,  1.406E-01,  1.413E-01,  1.420E-01,  &
   &      1.427E-01,  1.434E-01,  1.442E-01,  1.449E-01,  1.457E-01,  &
   &      1.465E-01,  1.472E-01,  1.479E-01,  1.487E-01,  1.495E-01,  &
   &      1.502E-01,  1.509E-01,  1.516E-01,  1.523E-01,  1.530E-01,  &
   &      1.539E-01,  1.547E-01,  1.555E-01,  1.563E-01,  1.571E-01,  &
   &      1.580E-01,  1.588E-01,  1.596E-01,  1.605E-01,  1.614E-01,  &
   &      1.623E-01,  1.632E-01,  1.641E-01,  1.649E-01,  1.658E-01/
   DATA (X(I),I= 1201, 1250)/                                        &
   &      1.666E-01,  1.675E-01,  1.684E-01,  1.692E-01,  1.701E-01,  &
   &      1.710E-01,  1.719E-01,  1.728E-01,  1.737E-01,  1.746E-01,  &
   &      1.756E-01,  1.764E-01,  1.774E-01,  1.783E-01,  1.792E-01,  &
   &      1.801E-01,  1.810E-01,  1.820E-01,  1.829E-01,  1.838E-01,  &
   &      1.848E-01,  1.857E-01,  1.866E-01,  1.876E-01,  1.885E-01,  &
   &      1.893E-01,  1.902E-01,  1.911E-01,  1.920E-01,  1.928E-01,  &
   &      1.936E-01,  1.945E-01,  1.953E-01,  1.961E-01,  1.969E-01,  &
   &      1.978E-01,  1.986E-01,  1.994E-01,  2.002E-01,  2.010E-01,  &
   &      2.018E-01,  2.026E-01,  2.034E-01,  2.041E-01,  2.049E-01,  &
   &      2.057E-01,  2.065E-01,  2.073E-01,  2.081E-01,  2.089E-01/
   DATA (X(I),I= 1251, 1300)/                                        &
   &      2.097E-01,  2.105E-01,  2.113E-01,  2.121E-01,  2.129E-01,  &
   &      2.137E-01,  2.146E-01,  2.154E-01,  2.163E-01,  2.172E-01,  &
   &      2.180E-01,  2.190E-01,  2.198E-01,  2.207E-01,  2.216E-01,  &
   &      2.225E-01,  2.234E-01,  2.243E-01,  2.251E-01,  2.260E-01,  &
   &      2.269E-01,  2.277E-01,  2.285E-01,  2.294E-01,  2.302E-01,  &
   &      2.311E-01,  2.320E-01,  2.328E-01,  2.337E-01,  2.346E-01,  &
   &      2.355E-01,  2.364E-01,  2.372E-01,  2.381E-01,  2.390E-01,  &
   &      2.398E-01,  2.407E-01,  2.416E-01,  2.424E-01,  2.432E-01,  &
   &      2.440E-01,  2.448E-01,  2.456E-01,  2.464E-01,  2.473E-01,  &
   &      2.482E-01,  2.491E-01,  2.500E-01,  2.509E-01,  2.517E-01/
   DATA (X(I),I= 1301, 1350)/                                        &
   &      2.525E-01,  2.533E-01,  2.541E-01,  2.550E-01,  2.559E-01,  &
   &      2.568E-01,  2.577E-01,  2.587E-01,  2.597E-01,  2.607E-01,  &
   &      2.617E-01,  2.626E-01,  2.636E-01,  2.645E-01,  2.654E-01,  &
   &      2.663E-01,  2.672E-01,  2.682E-01,  2.692E-01,  2.703E-01,  &
   &      2.713E-01,  2.724E-01,  2.734E-01,  2.744E-01,  2.754E-01,  &
   &      2.764E-01,  2.774E-01,  2.784E-01,  2.795E-01,  2.806E-01,  &
   &      2.816E-01,  2.827E-01,  2.838E-01,  2.850E-01,  2.861E-01,  &
   &      2.872E-01,  2.884E-01,  2.895E-01,  2.907E-01,  2.918E-01,  &
   &      2.930E-01,  2.942E-01,  2.954E-01,  2.967E-01,  2.980E-01,  &
   &      2.993E-01,  3.005E-01,  3.017E-01,  3.029E-01,  3.041E-01/
   DATA (X(I),I= 1351, 1400)/                                        &
   &      3.052E-01,  3.064E-01,  3.076E-01,  3.088E-01,  3.100E-01,  &
   &      3.112E-01,  3.124E-01,  3.136E-01,  3.149E-01,  3.161E-01,  &
   &      3.173E-01,  3.185E-01,  3.196E-01,  3.208E-01,  3.219E-01,  &
   &      3.230E-01,  3.242E-01,  3.253E-01,  3.265E-01,  3.277E-01,  &
   &      3.289E-01,  3.300E-01,  3.312E-01,  3.323E-01,  3.334E-01,  &
   &      3.345E-01,  3.356E-01,  3.367E-01,  3.378E-01,  3.389E-01,  &
   &      3.400E-01,  3.410E-01,  3.421E-01,  3.431E-01,  3.441E-01,  &
   &      3.452E-01,  3.462E-01,  3.472E-01,  3.482E-01,  3.493E-01,  &
   &      3.503E-01,  3.513E-01,  3.524E-01,  3.534E-01,  3.545E-01,  &
   &      3.555E-01,  3.566E-01,  3.577E-01,  3.588E-01,  3.599E-01/
   DATA (X(I),I= 1401, 1450)/                                        &
   &      3.611E-01,  3.622E-01,  3.632E-01,  3.643E-01,  3.653E-01,  &
   &      3.664E-01,  3.674E-01,  3.683E-01,  3.692E-01,  3.702E-01,  &
   &      3.711E-01,  3.721E-01,  3.732E-01,  3.740E-01,  3.751E-01,  &
   &      3.762E-01,  3.774E-01,  3.781E-01,  3.792E-01,  3.800E-01,  &
   &      3.811E-01,  3.818E-01,  3.826E-01,  3.837E-01,  3.844E-01,  &
   &      3.853E-01,  3.863E-01,  3.870E-01,  3.881E-01,  3.889E-01,  &
   &      3.898E-01,  3.908E-01,  3.916E-01,  3.928E-01,  3.937E-01,  &
   &      3.946E-01,  3.957E-01,  3.967E-01,  3.976E-01,  3.983E-01,  &
   &      3.995E-01,  4.002E-01,  4.013E-01,  4.023E-01,  4.032E-01,  &
   &      4.042E-01,  4.050E-01,  4.062E-01,  4.073E-01,  4.084E-01/
   DATA (X(I),I= 1451, 1500)/                                        &
   &      4.095E-01,  4.106E-01,  4.117E-01,  4.130E-01,  4.144E-01,  &
   &      4.155E-01,  4.165E-01,  4.178E-01,  4.189E-01,  4.200E-01,  &
   &      4.215E-01,  4.226E-01,  4.241E-01,  4.252E-01,  4.266E-01,  &
   &      4.280E-01,  4.293E-01,  4.306E-01,  4.321E-01,  4.335E-01,  &
   &      4.347E-01,  4.362E-01,  4.376E-01,  4.388E-01,  4.403E-01,  &
   &      4.418E-01,  4.433E-01,  4.450E-01,  4.465E-01,  4.480E-01,  &
   &      4.495E-01,  4.510E-01,  4.524E-01,  4.540E-01,  4.555E-01,  &
   &      4.571E-01,  4.585E-01,  4.603E-01,  4.619E-01,  4.634E-01,  &
   &      4.652E-01,  4.668E-01,  4.683E-01,  4.700E-01,  4.716E-01,  &
   &      4.734E-01,  4.750E-01,  4.764E-01,  4.782E-01,  4.797E-01/
   DATA (X(I),I= 1501, 1550)/                                        &
   &      4.812E-01,  4.830E-01,  4.845E-01,  4.861E-01,  4.875E-01,  &
   &      4.893E-01,  4.908E-01,  4.920E-01,  4.935E-01,  4.950E-01,  &
   &      4.964E-01,  4.976E-01,  4.990E-01,  5.003E-01,  5.016E-01,  &
   &      5.028E-01,  5.043E-01,  5.054E-01,  5.064E-01,  5.073E-01,  &
   &      5.082E-01,  5.094E-01,  5.104E-01,  5.111E-01,  5.119E-01,  &
   &      5.126E-01,  5.133E-01,  5.141E-01,  5.146E-01,  5.149E-01,  &
   &      5.154E-01,  5.159E-01,  5.162E-01,  5.166E-01,  5.167E-01,  &
   &      5.168E-01,  5.170E-01,  5.171E-01,  5.171E-01,  5.169E-01,  &
   &      5.166E-01,  5.162E-01,  5.159E-01,  5.155E-01,  5.151E-01,  &
   &      5.146E-01,  5.139E-01,  5.133E-01,  5.128E-01,  5.120E-01/
   DATA (X(I),I= 1551, 1600)/                                        &
   &      5.113E-01,  5.103E-01,  5.092E-01,  5.080E-01,  5.070E-01,  &
   &      5.059E-01,  5.048E-01,  5.036E-01,  5.022E-01,  5.011E-01,  &
   &      4.997E-01,  4.984E-01,  4.969E-01,  4.954E-01,  4.939E-01,  &
   &      4.924E-01,  4.909E-01,  4.893E-01,  4.877E-01,  4.863E-01,  &
   &      4.845E-01,  4.831E-01,  4.814E-01,  4.798E-01,  4.781E-01,  &
   &      4.766E-01,  4.748E-01,  4.732E-01,  4.718E-01,  4.703E-01,  &
   &      4.688E-01,  4.673E-01,  4.658E-01,  4.643E-01,  4.630E-01,  &
   &      4.615E-01,  4.600E-01,  4.586E-01,  4.572E-01,  4.559E-01,  &
   &      4.548E-01,  4.536E-01,  4.524E-01,  4.512E-01,  4.501E-01,  &
   &      4.491E-01,  4.483E-01,  4.475E-01,  4.468E-01,  4.459E-01/
   DATA (X(I),I= 1601, 1650)/                                        &
   &      4.450E-01,  4.444E-01,  4.438E-01,  4.431E-01,  4.424E-01,  &
   &      4.416E-01,  4.412E-01,  4.409E-01,  4.405E-01,  4.401E-01,  &
   &      4.397E-01,  4.394E-01,  4.392E-01,  4.390E-01,  4.389E-01,  &
   &      4.386E-01,  4.386E-01,  4.384E-01,  4.384E-01,  4.385E-01,  &
   &      4.385E-01,  4.385E-01,  4.385E-01,  4.387E-01,  4.387E-01,  &
   &      4.387E-01,  4.387E-01,  4.387E-01,  4.387E-01,  4.390E-01,  &
   &      4.391E-01,  4.394E-01,  4.398E-01,  4.398E-01,  4.402E-01,  &
   &      4.406E-01,  4.410E-01,  4.413E-01,  4.417E-01,  4.421E-01,  &
   &      4.425E-01,  4.428E-01,  4.432E-01,  4.440E-01,  4.443E-01,  &
   &      4.448E-01,  4.452E-01,  4.459E-01,  4.467E-01,  4.471E-01/
   DATA (X(I),I= 1651, 1700)/                                        &
   &      4.479E-01,  4.486E-01,  4.491E-01,  4.498E-01,  4.505E-01,  &
   &      4.512E-01,  4.519E-01,  4.525E-01,  4.532E-01,  4.539E-01,  &
   &      4.547E-01,  4.554E-01,  4.562E-01,  4.569E-01,  4.577E-01,  &
   &      4.584E-01,  4.592E-01,  4.599E-01,  4.606E-01,  4.614E-01,  &
   &      4.621E-01,  4.629E-01,  4.636E-01,  4.640E-01,  4.648E-01,  &
   &      4.655E-01,  4.662E-01,  4.670E-01,  4.675E-01,  4.682E-01,  &
   &      4.689E-01,  4.697E-01,  4.701E-01,  4.708E-01,  4.712E-01,  &
   &      4.718E-01,  4.724E-01,  4.729E-01,  4.735E-01,  4.739E-01,  &
   &      4.742E-01,  4.745E-01,  4.748E-01,  4.751E-01,  4.753E-01,  &
   &      4.755E-01,  4.757E-01,  4.757E-01,  4.757E-01,  4.756E-01/
   DATA (X(I),I= 1701, 1750)/                                        &
   &      4.756E-01,  4.756E-01,  4.753E-01,  4.752E-01,  4.749E-01,  &
   &      4.747E-01,  4.744E-01,  4.741E-01,  4.737E-01,  4.734E-01,  &
   &      4.730E-01,  4.725E-01,  4.721E-01,  4.715E-01,  4.708E-01,  &
   &      4.701E-01,  4.693E-01,  4.686E-01,  4.681E-01,  4.673E-01,  &
   &      4.663E-01,  4.657E-01,  4.649E-01,  4.641E-01,  4.632E-01,  &
   &      4.623E-01,  4.615E-01,  4.606E-01,  4.596E-01,  4.588E-01,  &
   &      4.579E-01,  4.569E-01,  4.561E-01,  4.551E-01,  4.542E-01,  &
   &      4.532E-01,  4.524E-01,  4.513E-01,  4.506E-01,  4.498E-01,  &
   &      4.487E-01,  4.479E-01,  4.472E-01,  4.461E-01,  4.454E-01,  &
   &      4.443E-01,  4.435E-01,  4.428E-01,  4.418E-01,  4.411E-01/
   DATA (X(I),I= 1751, 1800)/                                        &
   &      4.400E-01,  4.388E-01,  4.380E-01,  4.368E-01,  4.357E-01,  &
   &      4.347E-01,  4.338E-01,  4.328E-01,  4.316E-01,  4.305E-01,  &
   &      4.294E-01,  4.283E-01,  4.272E-01,  4.261E-01,  4.249E-01,  &
   &      4.235E-01,  4.222E-01,  4.212E-01,  4.201E-01,  4.186E-01,  &
   &      4.171E-01,  4.159E-01,  4.145E-01,  4.130E-01,  4.115E-01,  &
   &      4.100E-01,  4.085E-01,  4.070E-01,  4.057E-01,  4.042E-01,  &
   &      4.028E-01,  4.014E-01,  3.998E-01,  3.982E-01,  3.967E-01,  &
   &      3.950E-01,  3.935E-01,  3.919E-01,  3.904E-01,  3.892E-01,  &
   &      3.878E-01,  3.863E-01,  3.848E-01,  3.833E-01,  3.818E-01,  &
   &      3.803E-01,  3.789E-01,  3.775E-01,  3.761E-01,  3.746E-01/
   DATA (X(I),I= 1801, 1850)/                                        &
   &      3.731E-01,  3.718E-01,  3.706E-01,  3.694E-01,  3.681E-01,  &
   &      3.669E-01,  3.657E-01,  3.646E-01,  3.635E-01,  3.624E-01,  &
   &      3.613E-01,  3.603E-01,  3.592E-01,  3.581E-01,  3.571E-01,  &
   &      3.561E-01,  3.550E-01,  3.540E-01,  3.530E-01,  3.520E-01,  &
   &      3.510E-01,  3.503E-01,  3.495E-01,  3.487E-01,  3.479E-01,  &
   &      3.472E-01,  3.464E-01,  3.457E-01,  3.450E-01,  3.443E-01,  &
   &      3.436E-01,  3.429E-01,  3.422E-01,  3.415E-01,  3.409E-01,  &
   &      3.403E-01,  3.398E-01,  3.392E-01,  3.386E-01,  3.380E-01,  &
   &      3.375E-01,  3.369E-01,  3.364E-01,  3.358E-01,  3.353E-01,  &
   &      3.347E-01,  3.342E-01,  3.337E-01,  3.332E-01,  3.327E-01/
   DATA (X(I),I= 1851, 1900)/                                        &
   &      3.323E-01,  3.318E-01,  3.313E-01,  3.308E-01,  3.304E-01,  &
   &      3.300E-01,  3.295E-01,  3.291E-01,  3.286E-01,  3.282E-01,  &
   &      3.278E-01,  3.275E-01,  3.271E-01,  3.267E-01,  3.264E-01,  &
   &      3.260E-01,  3.256E-01,  3.251E-01,  3.245E-01,  3.240E-01,  &
   &      3.235E-01,  3.230E-01,  3.224E-01,  3.219E-01,  3.213E-01,  &
   &      3.207E-01,  3.202E-01,  3.196E-01,  3.190E-01,  3.185E-01,  &
   &      3.179E-01,  3.174E-01,  3.169E-01,  3.163E-01,  3.158E-01,  &
   &      3.152E-01,  3.146E-01,  3.139E-01,  3.132E-01,  3.125E-01,  &
   &      3.117E-01,  3.110E-01,  3.102E-01,  3.095E-01,  3.087E-01,  &
   &      3.079E-01,  3.071E-01,  3.063E-01,  3.055E-01,  3.048E-01/
   DATA (X(I),I= 1901, 1950)/                                        &
   &      3.039E-01,  3.031E-01,  3.022E-01,  3.014E-01,  3.005E-01,  &
   &      2.996E-01,  2.988E-01,  2.979E-01,  2.970E-01,  2.961E-01,  &
   &      2.952E-01,  2.944E-01,  2.935E-01,  2.927E-01,  2.920E-01,  &
   &      2.913E-01,  2.906E-01,  2.900E-01,  2.893E-01,  2.886E-01,  &
   &      2.880E-01,  2.874E-01,  2.869E-01,  2.863E-01,  2.858E-01,  &
   &      2.852E-01,  2.847E-01,  2.842E-01,  2.838E-01,  2.834E-01,  &
   &      2.830E-01,  2.826E-01,  2.822E-01,  2.818E-01,  2.815E-01,  &
   &      2.813E-01,  2.811E-01,  2.809E-01,  2.807E-01,  2.805E-01,  &
   &      2.803E-01,  2.802E-01,  2.803E-01,  2.803E-01,  2.803E-01,  &
   &      2.803E-01,  2.803E-01,  2.804E-01,  2.804E-01,  2.805E-01/
   DATA (X(I),I= 1951, 2000)/                                        &
   &      2.806E-01,  2.807E-01,  2.808E-01,  2.809E-01,  2.810E-01,  &
   &      2.810E-01,  2.809E-01,  2.808E-01,  2.808E-01,  2.807E-01,  &
   &      2.806E-01,  2.805E-01,  2.804E-01,  2.801E-01,  2.799E-01,  &
   &      2.796E-01,  2.794E-01,  2.791E-01,  2.789E-01,  2.785E-01,  &
   &      2.780E-01,  2.775E-01,  2.770E-01,  2.765E-01,  2.760E-01,  &
   &      2.755E-01,  2.749E-01,  2.741E-01,  2.734E-01,  2.726E-01,  &
   &      2.718E-01,  2.710E-01,  2.703E-01,  2.694E-01,  2.684E-01,  &
   &      2.674E-01,  2.664E-01,  2.654E-01,  2.644E-01,  2.634E-01,  &
   &      2.624E-01,  2.612E-01,  2.601E-01,  2.589E-01,  2.578E-01,  &
   &      2.566E-01,  2.555E-01,  2.543E-01,  2.530E-01,  2.517E-01/
   DATA (X(I),I= 2001, 2050)/                                        &
   &      2.503E-01,  2.490E-01,  2.477E-01,  2.463E-01,  2.450E-01,  &
   &      2.436E-01,  2.423E-01,  2.410E-01,  2.396E-01,  2.383E-01,  &
   &      2.370E-01,  2.356E-01,  2.342E-01,  2.328E-01,  2.314E-01,  &
   &      2.299E-01,  2.285E-01,  2.271E-01,  2.257E-01,  2.243E-01,  &
   &      2.230E-01,  2.218E-01,  2.205E-01,  2.192E-01,  2.179E-01,  &
   &      2.166E-01,  2.153E-01,  2.140E-01,  2.126E-01,  2.113E-01,  &
   &      2.100E-01,  2.086E-01,  2.073E-01,  2.060E-01,  2.048E-01,  &
   &      2.036E-01,  2.025E-01,  2.013E-01,  2.001E-01,  1.989E-01,  &
   &      1.977E-01,  1.967E-01,  1.957E-01,  1.947E-01,  1.937E-01,  &
   &      1.927E-01,  1.917E-01,  1.907E-01,  1.898E-01,  1.890E-01/
   DATA (X(I),I= 2051, 2100)/                                        &
   &      1.883E-01,  1.875E-01,  1.867E-01,  1.859E-01,  1.851E-01,  &
   &      1.844E-01,  1.836E-01,  1.829E-01,  1.821E-01,  1.813E-01,  &
   &      1.806E-01,  1.798E-01,  1.791E-01,  1.786E-01,  1.781E-01,  &
   &      1.776E-01,  1.771E-01,  1.766E-01,  1.761E-01,  1.756E-01,  &
   &      1.751E-01,  1.747E-01,  1.743E-01,  1.739E-01,  1.735E-01,  &
   &      1.731E-01,  1.726E-01,  1.722E-01,  1.718E-01,  1.714E-01,  &
   &      1.710E-01,  1.706E-01,  1.703E-01,  1.698E-01,  1.695E-01,  &
   &      1.691E-01,  1.686E-01,  1.682E-01,  1.677E-01,  1.673E-01,  &
   &      1.668E-01,  1.664E-01,  1.660E-01,  1.655E-01,  1.650E-01,  &
   &      1.646E-01,  1.642E-01,  1.637E-01,  1.633E-01,  1.628E-01/
   DATA (X(I),I= 2101, 2150)/                                        &
   &      1.624E-01,  1.619E-01,  1.615E-01,  1.611E-01,  1.607E-01,  &
   &      1.602E-01,  1.598E-01,  1.593E-01,  1.588E-01,  1.583E-01,  &
   &      1.578E-01,  1.574E-01,  1.569E-01,  1.564E-01,  1.560E-01,  &
   &      1.555E-01,  1.551E-01,  1.548E-01,  1.544E-01,  1.540E-01,  &
   &      1.536E-01,  1.532E-01,  1.528E-01,  1.525E-01,  1.523E-01,  &
   &      1.520E-01,  1.518E-01,  1.516E-01,  1.514E-01,  1.511E-01,  &
   &      1.510E-01,  1.511E-01,  1.513E-01,  1.514E-01,  1.516E-01,  &
   &      1.518E-01,  1.519E-01,  1.521E-01,  1.523E-01,  1.526E-01,  &
   &      1.528E-01,  1.531E-01,  1.534E-01,  1.537E-01,  1.540E-01,  &
   &      1.543E-01,  1.547E-01,  1.551E-01,  1.555E-01,  1.560E-01/
   DATA (X(I),I= 2151, 2200)/                                        &
   &      1.564E-01,  1.568E-01,  1.572E-01,  1.575E-01,  1.579E-01,  &
   &      1.581E-01,  1.584E-01,  1.586E-01,  1.589E-01,  1.592E-01,  &
   &      1.594E-01,  1.596E-01,  1.598E-01,  1.599E-01,  1.600E-01,  &
   &      1.601E-01,  1.602E-01,  1.603E-01,  1.604E-01,  1.604E-01,  &
   &      1.602E-01,  1.600E-01,  1.598E-01,  1.596E-01,  1.594E-01,  &
   &      1.592E-01,  1.590E-01,  1.585E-01,  1.580E-01,  1.574E-01,  &
   &      1.568E-01,  1.562E-01,  1.555E-01,  1.549E-01,  1.543E-01,  &
   &      1.535E-01,  1.527E-01,  1.518E-01,  1.510E-01,  1.501E-01,  &
   &      1.493E-01,  1.484E-01,  1.475E-01,  1.464E-01,  1.453E-01,  &
   &      1.442E-01,  1.431E-01,  1.420E-01,  1.409E-01,  1.398E-01/
   DATA (X(I),I= 2201, 2250)/                                        &
   &      1.386E-01,  1.374E-01,  1.362E-01,  1.350E-01,  1.338E-01,  &
   &      1.326E-01,  1.314E-01,  1.302E-01,  1.290E-01,  1.278E-01,  &
   &      1.265E-01,  1.253E-01,  1.241E-01,  1.229E-01,  1.217E-01,  &
   &      1.204E-01,  1.192E-01,  1.180E-01,  1.169E-01,  1.157E-01,  &
   &      1.145E-01,  1.133E-01,  1.122E-01,  1.110E-01,  1.099E-01,  &
   &      1.089E-01,  1.079E-01,  1.069E-01,  1.060E-01,  1.050E-01,  &
   &      1.040E-01,  1.031E-01,  1.021E-01,  1.013E-01,  1.005E-01,  &
   &      9.973E-02,  9.897E-02,  9.820E-02,  9.743E-02,  9.664E-02,  &
   &      9.588E-02,  9.524E-02,  9.462E-02,  9.400E-02,  9.339E-02,  &
   &      9.279E-02,  9.217E-02,  9.158E-02,  9.098E-02,  9.049E-02/
   DATA (X(I),I= 2251, 2300)/                                        &
   &      9.002E-02,  8.958E-02,  8.913E-02,  8.869E-02,  8.827E-02,  &
   &      8.783E-02,  8.742E-02,  8.712E-02,  8.690E-02,  8.670E-02,  &
   &      8.648E-02,  8.629E-02,  8.607E-02,  8.588E-02,  8.568E-02,  &
   &      8.547E-02,  8.525E-02,  8.503E-02,  8.482E-02,  8.462E-02,  &
   &      8.440E-02,  8.418E-02,  8.397E-02,  8.379E-02,  8.369E-02,  &
   &      8.359E-02,  8.349E-02,  8.341E-02,  8.332E-02,  8.322E-02,  &
   &      8.316E-02,  8.305E-02,  8.288E-02,  8.269E-02,  8.251E-02,  &
   &      8.232E-02,  8.214E-02,  8.195E-02,  8.178E-02,  8.158E-02,  &
   &      8.133E-02,  8.108E-02,  8.083E-02,  8.057E-02,  8.031E-02,  &
   &      8.003E-02,  7.976E-02,  7.949E-02,  7.917E-02,  7.874E-02/
   DATA (X(I),I= 2301, 2350)/                                        &
   &      7.830E-02,  7.789E-02,  7.744E-02,  7.704E-02,  7.662E-02,  &
   &      7.620E-02,  7.579E-02,  7.549E-02,  7.519E-02,  7.490E-02,  &
   &      7.460E-02,  7.432E-02,  7.404E-02,  7.377E-02,  7.347E-02,  &
   &      7.333E-02,  7.329E-02,  7.323E-02,  7.318E-02,  7.315E-02,  &
   &      7.310E-02,  7.307E-02,  7.303E-02,  7.304E-02,  7.320E-02,  &
   &      7.336E-02,  7.351E-02,  7.369E-02,  7.386E-02,  7.402E-02,  &
   &      7.417E-02,  7.435E-02,  7.458E-02,  7.481E-02,  7.505E-02,  &
   &      7.529E-02,  7.556E-02,  7.582E-02,  7.607E-02,  7.634E-02,  &
   &      7.661E-02,  7.695E-02,  7.728E-02,  7.763E-02,  7.797E-02,  &
   &      7.830E-02,  7.864E-02,  7.895E-02,  7.929E-02,  7.952E-02/
   DATA (X(I),I= 2351, 2400)/                                        &
   &      7.975E-02,  7.996E-02,  8.018E-02,  8.039E-02,  8.061E-02,  &
   &      8.081E-02,  8.102E-02,  8.115E-02,  8.105E-02,  8.096E-02,  &
   &      8.086E-02,  8.074E-02,  8.062E-02,  8.048E-02,  8.033E-02,  &
   &      8.019E-02,  7.988E-02,  7.946E-02,  7.907E-02,  7.866E-02,  &
   &      7.824E-02,  7.784E-02,  7.743E-02,  7.702E-02,  7.658E-02,  &
   &      7.607E-02,  7.552E-02,  7.497E-02,  7.442E-02,  7.386E-02,  &
   &      7.328E-02,  7.271E-02,  7.214E-02,  7.149E-02,  7.072E-02,  &
   &      6.995E-02,  6.917E-02,  6.840E-02,  6.762E-02,  6.684E-02,  &
   &      6.605E-02,  6.527E-02,  6.453E-02,  6.382E-02,  6.311E-02,  &
   &      6.240E-02,  6.169E-02,  6.099E-02,  6.027E-02,  5.956E-02/
   DATA (X(I),I= 2401, 2450)/                                        &
   &      5.886E-02,  5.817E-02,  5.747E-02,  5.676E-02,  5.607E-02,  &
   &      5.537E-02,  5.469E-02,  5.400E-02,  5.331E-02,  5.265E-02,  &
   &      5.205E-02,  5.146E-02,  5.088E-02,  5.028E-02,  4.970E-02,  &
   &      4.912E-02,  4.853E-02,  4.794E-02,  4.738E-02,  4.685E-02,  &
   &      4.634E-02,  4.583E-02,  4.532E-02,  4.480E-02,  4.431E-02,  &
   &      4.380E-02,  4.331E-02,  4.291E-02,  4.260E-02,  4.230E-02,  &
   &      4.199E-02,  4.168E-02,  4.137E-02,  4.107E-02,  4.078E-02,  &
   &      4.047E-02,  4.029E-02,  4.016E-02,  4.004E-02,  3.991E-02,  &
   &      3.980E-02,  3.970E-02,  3.959E-02,  3.949E-02,  3.937E-02,  &
   &      3.934E-02,  3.933E-02,  3.933E-02,  3.934E-02,  3.935E-02/
   DATA (X(I),I= 2451, 2500)/                                        &
   &      3.935E-02,  3.936E-02,  3.936E-02,  3.936E-02,  3.931E-02,  &
   &      3.925E-02,  3.918E-02,  3.912E-02,  3.905E-02,  3.897E-02,  &
   &      3.889E-02,  3.881E-02,  3.874E-02,  3.866E-02,  3.855E-02,  &
   &      3.846E-02,  3.837E-02,  3.826E-02,  3.818E-02,  3.807E-02,  &
   &      3.795E-02,  3.786E-02,  3.769E-02,  3.748E-02,  3.727E-02,  &
   &      3.706E-02,  3.686E-02,  3.664E-02,  3.643E-02,  3.622E-02,  &
   &      3.601E-02,  3.581E-02,  3.561E-02,  3.542E-02,  3.522E-02,  &
   &      3.503E-02,  3.484E-02,  3.465E-02,  3.446E-02,  3.427E-02,  &
   &      3.407E-02,  3.386E-02,  3.364E-02,  3.343E-02,  3.322E-02,  &
   &      3.301E-02,  3.280E-02,  3.259E-02,  3.238E-02,  3.221E-02/
   DATA (X(I),I= 2501, 2550)/                                        &
   &      3.209E-02,  3.198E-02,  3.186E-02,  3.175E-02,  3.164E-02,  &
   &      3.153E-02,  3.143E-02,  3.132E-02,  3.126E-02,  3.136E-02,  &
   &      3.148E-02,  3.159E-02,  3.170E-02,  3.182E-02,  3.194E-02,  &
   &      3.206E-02,  3.219E-02,  3.232E-02,  3.253E-02,  3.275E-02,  &
   &      3.298E-02,  3.320E-02,  3.343E-02,  3.366E-02,  3.389E-02,  &
   &      3.412E-02,  3.435E-02,  3.455E-02,  3.475E-02,  3.495E-02,  &
   &      3.515E-02,  3.535E-02,  3.555E-02,  3.573E-02,  3.593E-02,  &
   &      3.612E-02,  3.625E-02,  3.628E-02,  3.631E-02,  3.634E-02,  &
   &      3.637E-02,  3.639E-02,  3.641E-02,  3.643E-02,  3.646E-02,  &
   &      3.646E-02,  3.636E-02,  3.623E-02,  3.611E-02,  3.599E-02/
   DATA (X(I),I= 2551, 2600)/                                        &
   &      3.586E-02,  3.572E-02,  3.559E-02,  3.545E-02,  3.531E-02,  &
   &      3.509E-02,  3.481E-02,  3.453E-02,  3.425E-02,  3.398E-02,  &
   &      3.369E-02,  3.341E-02,  3.312E-02,  3.284E-02,  3.254E-02,  &
   &      3.217E-02,  3.180E-02,  3.143E-02,  3.106E-02,  3.069E-02,  &
   &      3.031E-02,  2.994E-02,  2.956E-02,  2.919E-02,  2.882E-02,  &
   &      2.845E-02,  2.808E-02,  2.771E-02,  2.734E-02,  2.696E-02,  &
   &      2.660E-02,  2.622E-02,  2.585E-02,  2.549E-02,  2.512E-02,  &
   &      2.476E-02,  2.440E-02,  2.404E-02,  2.368E-02,  2.332E-02,  &
   &      2.297E-02,  2.261E-02,  2.225E-02,  2.193E-02,  2.164E-02,  &
   &      2.135E-02,  2.106E-02,  2.076E-02,  2.048E-02,  2.019E-02/
   DATA (X(I),I= 2601, 2650)/                                        &
   &      1.990E-02,  1.962E-02,  1.935E-02,  1.916E-02,  1.898E-02,  &
   &      1.881E-02,  1.863E-02,  1.846E-02,  1.828E-02,  1.811E-02,  &
   &      1.794E-02,  1.777E-02,  1.764E-02,  1.757E-02,  1.749E-02,  &
   &      1.742E-02,  1.735E-02,  1.728E-02,  1.721E-02,  1.714E-02,  &
   &      1.708E-02,  1.701E-02,  1.699E-02,  1.699E-02,  1.699E-02,  &
   &      1.699E-02,  1.699E-02,  1.699E-02,  1.699E-02,  1.699E-02,  &
   &      1.699E-02,  1.699E-02,  1.700E-02,  1.702E-02,  1.703E-02,  &
   &      1.704E-02,  1.706E-02,  1.707E-02,  1.708E-02,  1.709E-02,  &
   &      1.710E-02,  1.710E-02,  1.706E-02,  1.701E-02,  1.696E-02,  &
   &      1.692E-02,  1.687E-02,  1.683E-02,  1.678E-02,  1.673E-02/
   DATA (X(I),I= 2651, 2700)/                                        &
   &      1.668E-02,  1.661E-02,  1.651E-02,  1.642E-02,  1.632E-02,  &
   &      1.622E-02,  1.612E-02,  1.602E-02,  1.592E-02,  1.582E-02,  &
   &      1.572E-02,  1.560E-02,  1.545E-02,  1.531E-02,  1.517E-02,  &
   &      1.503E-02,  1.489E-02,  1.474E-02,  1.460E-02,  1.446E-02,  &
   &      1.432E-02,  1.420E-02,  1.408E-02,  1.397E-02,  1.386E-02,  &
   &      1.375E-02,  1.363E-02,  1.352E-02,  1.341E-02,  1.329E-02,  &
   &      1.318E-02,  1.313E-02,  1.310E-02,  1.308E-02,  1.305E-02,  &
   &      1.303E-02,  1.300E-02,  1.298E-02,  1.295E-02,  1.292E-02,  &
   &      1.290E-02,  1.293E-02,  1.297E-02,  1.302E-02,  1.307E-02,  &
   &      1.311E-02,  1.316E-02,  1.320E-02,  1.325E-02,  1.330E-02/
   DATA (X(I),I= 2701, 2750)/                                        &
   &      1.334E-02,  1.341E-02,  1.349E-02,  1.357E-02,  1.366E-02,  &
   &      1.374E-02,  1.382E-02,  1.390E-02,  1.398E-02,  1.406E-02,  &
   &      1.414E-02,  1.421E-02,  1.427E-02,  1.433E-02,  1.438E-02,  &
   &      1.444E-02,  1.450E-02,  1.456E-02,  1.462E-02,  1.467E-02,  &
   &      1.473E-02,  1.475E-02,  1.474E-02,  1.473E-02,  1.472E-02,  &
   &      1.471E-02,  1.470E-02,  1.468E-02,  1.467E-02,  1.466E-02,  &
   &      1.465E-02,  1.461E-02,  1.452E-02,  1.442E-02,  1.433E-02,  &
   &      1.423E-02,  1.414E-02,  1.405E-02,  1.395E-02,  1.386E-02,  &
   &      1.377E-02,  1.367E-02,  1.356E-02,  1.344E-02,  1.333E-02,  &
   &      1.321E-02,  1.310E-02,  1.298E-02,  1.287E-02,  1.275E-02/
   DATA (X(I),I= 2751, 2800)/                                        &
   &      1.264E-02,  1.252E-02,  1.236E-02,  1.220E-02,  1.204E-02,  &
   &      1.187E-02,  1.171E-02,  1.155E-02,  1.138E-02,  1.122E-02,  &
   &      1.106E-02,  1.089E-02,  1.072E-02,  1.055E-02,  1.038E-02,  &
   &      1.021E-02,  1.004E-02,  9.872E-03,  9.701E-03,  9.531E-03,  &
   &      9.361E-03,  9.190E-03,  9.029E-03,  8.896E-03,  8.763E-03,  &
   &      8.634E-03,  8.503E-03,  8.370E-03,  8.240E-03,  8.108E-03,  &
   &      7.977E-03,  7.847E-03,  7.717E-03,  7.622E-03,  7.540E-03,  &
   &      7.457E-03,  7.373E-03,  7.291E-03,  7.209E-03,  7.126E-03,  &
   &      7.041E-03,  6.961E-03,  6.878E-03,  6.813E-03,  6.782E-03,  &
   &      6.751E-03,  6.721E-03,  6.690E-03,  6.658E-03,  6.628E-03/
   DATA (X(I),I= 2801, 2850)/                                        &
   &      6.599E-03,  6.567E-03,  6.536E-03,  6.508E-03,  6.499E-03,  &
   &      6.495E-03,  6.493E-03,  6.490E-03,  6.488E-03,  6.484E-03,  &
   &      6.480E-03,  6.478E-03,  6.474E-03,  6.470E-03,  6.471E-03,  &
   &      6.473E-03,  6.476E-03,  6.480E-03,  6.482E-03,  6.486E-03,  &
   &      6.489E-03,  6.493E-03,  6.496E-03,  6.499E-03,  6.501E-03,  &
   &      6.488E-03,  6.465E-03,  6.444E-03,  6.422E-03,  6.401E-03,  &
   &      6.381E-03,  6.359E-03,  6.337E-03,  6.316E-03,  6.294E-03,  &
   &      6.266E-03,  6.203E-03,  6.135E-03,  6.067E-03,  6.002E-03,  &
   &      5.932E-03,  5.867E-03,  5.797E-03,  5.732E-03,  5.664E-03,  &
   &      5.596E-03,  5.528E-03,  5.457E-03,  5.388E-03,  5.318E-03/
   DATA (X(I),I= 2851, 2900)/                                        &
   &      5.248E-03,  5.177E-03,  5.107E-03,  5.036E-03,  4.966E-03,  &
   &      4.895E-03,  4.825E-03,  4.781E-03,  4.755E-03,  4.729E-03,  &
   &      4.703E-03,  4.675E-03,  4.648E-03,  4.622E-03,  4.595E-03,  &
   &      4.569E-03,  4.542E-03,  4.516E-03,  4.514E-03,  4.517E-03,  &
   &      4.520E-03,  4.523E-03,  4.527E-03,  4.530E-03,  4.533E-03,  &
   &      4.539E-03,  4.541E-03,  4.545E-03,  4.550E-03,  4.568E-03,  &
   &      4.588E-03,  4.609E-03,  4.628E-03,  4.649E-03,  4.669E-03,  &
   &      4.689E-03,  4.710E-03,  4.731E-03,  4.750E-03,  4.771E-03,  &
   &      4.796E-03,  4.822E-03,  4.847E-03,  4.869E-03,  4.896E-03,  &
   &      4.921E-03,  4.945E-03,  4.970E-03,  4.995E-03,  5.020E-03/
   DATA (X(I),I= 2901, 2950)/                                        &
   &      5.038E-03,  5.040E-03,  5.038E-03,  5.037E-03,  5.036E-03,  &
   &      5.036E-03,  5.034E-03,  5.034E-03,  5.031E-03,  5.031E-03,  &
   &      5.031E-03,  5.019E-03,  4.979E-03,  4.934E-03,  4.892E-03,  &
   &      4.848E-03,  4.805E-03,  4.763E-03,  4.718E-03,  4.676E-03,  &
   &      4.632E-03,  4.590E-03,  4.541E-03,  4.475E-03,  4.405E-03,  &
   &      4.336E-03,  4.268E-03,  4.198E-03,  4.130E-03,  4.060E-03,  &
   &      3.990E-03,  3.922E-03,  3.852E-03,  3.782E-03,  3.715E-03,  &
   &      3.646E-03,  3.577E-03,  3.508E-03,  3.439E-03,  3.370E-03,  &
   &      3.301E-03,  3.232E-03,  3.163E-03,  3.094E-03,  3.026E-03,  &
   &      2.971E-03,  2.919E-03,  2.868E-03,  2.816E-03,  2.764E-03/
   DATA (X(I),I= 2951, 3000)/                                        &
   &      2.712E-03,  2.661E-03,  2.609E-03,  2.557E-03,  2.505E-03,  &
   &      2.454E-03,  2.416E-03,  2.386E-03,  2.356E-03,  2.326E-03,  &
   &      2.297E-03,  2.267E-03,  2.237E-03,  2.207E-03,  2.177E-03,  &
   &      2.148E-03,  2.118E-03,  2.096E-03,  2.087E-03,  2.078E-03,  &
   &      2.070E-03,  2.061E-03,  2.052E-03,  2.043E-03,  2.034E-03,  &
   &      2.025E-03,  2.016E-03,  2.007E-03,  2.000E-03,  2.000E-03,  &
   &      2.001E-03,  2.002E-03,  2.003E-03,  2.004E-03,  2.005E-03,  &
   &      2.006E-03,  2.007E-03,  2.007E-03,  2.008E-03,  2.009E-03,  &
   &      2.008E-03,  2.006E-03,  2.003E-03,  2.001E-03,  1.999E-03,  &
   &      1.997E-03,  1.994E-03,  1.992E-03,  1.990E-03,  1.988E-03/
   DATA (X(I),I= 3001, 3050)/                                        &
   &      1.985E-03,  1.980E-03,  1.968E-03,  1.956E-03,  1.944E-03,  &
   &      1.932E-03,  1.919E-03,  1.907E-03,  1.895E-03,  1.883E-03,  &
   &      1.871E-03,  1.859E-03,  1.846E-03,  1.827E-03,  1.805E-03,  &
   &      1.783E-03,  1.761E-03,  1.740E-03,  1.718E-03,  1.696E-03,  &
   &      1.674E-03,  1.652E-03,  1.631E-03,  1.609E-03,  1.585E-03,  &
   &      1.557E-03,  1.529E-03,  1.500E-03,  1.472E-03,  1.444E-03,  &
   &      1.416E-03,  1.388E-03,  1.359E-03,  1.331E-03,  1.303E-03,  &
   &      1.275E-03,  1.251E-03,  1.228E-03,  1.206E-03,  1.183E-03,  &
   &      1.161E-03,  1.139E-03,  1.116E-03,  1.094E-03,  1.072E-03,  &
   &      1.049E-03,  1.027E-03,  1.007E-03,  1.003E-03,  9.991E-04/
   DATA (X(I),I= 3051, 3100)/                                        &
   &      9.959E-04,  9.926E-04,  9.891E-04,  9.857E-04,  9.826E-04,  &
   &      9.790E-04,  9.758E-04,  9.725E-04,  9.692E-04,  9.720E-04,  &
   &      9.842E-04,  9.967E-04,  1.009E-03,  1.022E-03,  1.034E-03,  &
   &      1.047E-03,  1.059E-03,  1.071E-03,  1.084E-03,  1.096E-03,  &
   &      1.109E-03,  1.118E-03,  1.126E-03,  1.133E-03,  1.141E-03,  &
   &      1.148E-03,  1.156E-03,  1.163E-03,  1.171E-03,  1.178E-03,  &
   &      1.186E-03,  1.193E-03,  1.200E-03,  1.194E-03,  1.185E-03,  &
   &      1.176E-03,  1.166E-03,  1.157E-03,  1.148E-03,  1.138E-03,  &
   &      1.129E-03,  1.120E-03,  1.110E-03,  1.101E-03,  1.091E-03,  &
   &      1.076E-03,  1.060E-03,  1.044E-03,  1.028E-03,  1.012E-03/
   DATA (X(I),I= 3101, 3150)/                                        &
   &      9.955E-04,  9.794E-04,  9.632E-04,  9.473E-04,  9.310E-04,  &
   &      9.151E-04,  8.982E-04,  8.770E-04,  8.556E-04,  8.339E-04,  &
   &      8.121E-04,  7.907E-04,  7.689E-04,  7.471E-04,  7.255E-04,  &
   &      7.038E-04,  6.822E-04,  6.606E-04,  6.376E-04,  6.083E-04,  &
   &      5.783E-04,  5.478E-04,  5.177E-04,  4.875E-04,  4.574E-04,  &
   &      4.272E-04,  3.969E-04,  3.668E-04,  3.365E-04,  3.063E-04,  &
   &      2.750E-04,  2.500E-04,  2.250E-04,  2.000E-04,  1.850E-04,  &
   &      1.700E-04,  1.550E-04,  1.400E-04,  1.250E-04,  1.100E-04,  &
   &      0.950E-04,  0.825E-04,  0.700E-04,  0.575E-04,  0.400E-04,  &
   &      0.275E-04,  0.175E-04,  0.100E-04,  0.040E-04,  0.00000  /
!
   DATA (Y(I),I=    1,   50)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I= 51, 100)/                                           &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I=  101,  150)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I= 151, 200)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I=  201,  250)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I= 251, 300)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I=  301,  350)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I= 351, 400)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I=  401,  450)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I= 451, 500)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I=  501,  550)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I= 551, 600)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I=  601,  650)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I= 651, 700)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I=  701,  750)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I= 751, 800)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Y(I),I=  801,  850)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.025e-05,  0.050e-05,  0.075e-05,  0.100e-05,  0.125e-05,  &
   &      0.150e-05,  0.175e-05,  0.225e-05,  0.250e-05,  0.300e-05,  &
   &      0.350e-05,  0.400e-05,  0.500e-05,  0.650e-05,  0.800e-05,  &
   &      1.100e-05,  1.400e-05,  1.800e-05,  2.250e-05,  2.650e-05,  &
   &      3.000e-05,  3.104E-05,  3.136E-05,  3.152E-05,  3.186E-05,  &
   &      3.213E-05,  3.229E-05,  3.206E-05,  3.156E-05,  3.063E-05/
   DATA (Y(I),I= 851, 900)/                                          &
   &      3.098E-05,  3.197E-05,  3.271E-05,  3.315E-05,  3.262E-05,  &
   &      3.201E-05,  3.129E-05,  3.148E-05,  3.206E-05,  3.175E-05,  &
   &      3.148E-05,  3.167E-05,  3.159E-05,  3.120E-05,  3.117E-05,  &
   &      3.109E-05,  3.041E-05,  2.995E-05,  3.011E-05,  3.004E-05,  &
   &      2.972E-05,  2.933E-05,  2.887E-05,  2.801E-05,  2.731E-05,  &
   &      2.735E-05,  2.656E-05,  2.712E-05,  2.576E-05,  2.565E-05,  &
   &      2.449E-05,  2.450E-05,  2.454E-05,  2.414E-05,  2.488E-05,  &
   &      2.413E-05,  2.297E-05,  2.298E-05,  2.335E-05,  2.259E-05,  &
   &      2.142E-05,  2.257E-05,  2.259E-05,  2.220E-05,  2.259E-05,  &
   &      2.336E-05,  2.299E-05,  2.262E-05,  2.298E-05,  2.298E-05/
   DATA (Y(I),I=  901,  950)/                                        &
   &      2.298E-05,  2.337E-05,  2.414E-05,  2.491E-05,  2.570E-05,  &
   &      2.492E-05,  2.531E-05,  2.567E-05,  2.647E-05,  2.839E-05,  &
   &      2.916E-05,  2.917E-05,  2.916E-05,  2.956E-05,  3.033E-05,  &
   &      3.110E-05,  3.108E-05,  3.071E-05,  3.226E-05,  3.303E-05,  &
   &      3.302E-05,  3.417E-05,  3.419E-05,  3.458E-05,  3.458E-05,  &
   &      3.457E-05,  3.613E-05,  3.769E-05,  3.768E-05,  3.845E-05,  &
   &      3.843E-05,  3.882E-05,  3.920E-05,  3.959E-05,  3.998E-05,  &
   &      4.037E-05,  4.036E-05,  4.074E-05,  4.074E-05,  4.036E-05,  &
   &      4.115E-05,  4.192E-05,  4.192E-05,  4.269E-05,  4.270E-05,  &
   &      4.346E-05,  4.346E-05,  4.346E-05,  4.463E-05,  4.423E-05/
   DATA (Y(I),I= 951, 1000)/                                         &
   &      4.502E-05,  4.539E-05,  4.501E-05,  4.465E-05,  4.388E-05,  &
   &      4.503E-05,  4.578E-05,  4.617E-05,  4.657E-05,  4.657E-05,  &
   &      4.617E-05,  4.657E-05,  4.655E-05,  4.695E-05,  4.808E-05,  &
   &      4.807E-05,  4.768E-05,  4.771E-05,  4.807E-05,  4.808E-05,  &
   &      4.808E-05,  4.731E-05,  4.731E-05,  4.772E-05,  4.808E-05,  &
   &      4.846E-05,  4.769E-05,  4.807E-05,  4.846E-05,  4.847E-05,  &
   &      4.845E-05,  4.885E-05,  4.886E-05,  4.848E-05,  4.810E-05,  &
   &      4.848E-05,  4.769E-05,  4.773E-05,  4.656E-05,  4.690E-05,  &
   &      4.770E-05,  4.730E-05,  4.690E-05,  4.728E-05,  4.687E-05,  &
   &      4.728E-05,  4.649E-05,  4.730E-05,  4.690E-05,  4.766E-05/
   DATA (Y(I),I= 1001, 1050)/                                        &
   &      4.804E-05,  4.727E-05,  4.649E-05,  4.727E-05,  4.647E-05,  &
   &      4.723E-05,  4.760E-05,  4.800E-05,  4.836E-05,  4.762E-05,  &
   &      4.760E-05,  4.834E-05,  4.837E-05,  4.796E-05,  4.796E-05,  &
   &      4.796E-05,  4.876E-05,  4.955E-05,  5.030E-05,  5.067E-05,  &
   &      5.106E-05,  5.183E-05,  5.104E-05,  5.185E-05,  5.224E-05,  &
   &      5.340E-05,  5.301E-05,  5.340E-05,  5.456E-05,  5.532E-05,  &
   &      5.532E-05,  5.531E-05,  5.608E-05,  5.762E-05,  5.799E-05,  &
   &      5.877E-05,  5.837E-05,  5.839E-05,  5.878E-05,  5.916E-05,  &
   &      6.032E-05,  6.109E-05,  6.226E-05,  6.342E-05,  6.305E-05,  &
   &      6.305E-05,  6.265E-05,  6.226E-05,  6.298E-05,  6.452E-05/
   DATA (Y(I),I= 1051, 1100)/                                        &
   &      6.531E-05,  6.572E-05,  6.572E-05,  6.614E-05,  6.611E-05,  &
   &      6.574E-05,  6.690E-05,  6.648E-05,  6.688E-05,  6.688E-05,  &
   &      6.764E-05,  6.841E-05,  6.843E-05,  6.918E-05,  6.919E-05,  &
   &      6.959E-05,  6.955E-05,  6.993E-05,  6.990E-05,  7.031E-05,  &
   &      6.914E-05,  6.835E-05,  6.912E-05,  7.105E-05,  7.222E-05,  &
   &      7.261E-05,  7.258E-05,  7.258E-05,  7.258E-05,  7.260E-05,  &
   &      7.179E-05,  7.257E-05,  7.180E-05,  7.142E-05,  7.181E-05,  &
   &      7.217E-05,  7.256E-05,  7.217E-05,  7.217E-05,  7.219E-05,  &
   &      7.256E-05,  7.180E-05,  7.178E-05,  7.255E-05,  7.216E-05,  &
   &      7.175E-05,  7.175E-05,  7.137E-05,  7.177E-05,  7.214E-05/
   DATA (Y(I),I= 1101, 1150)/                                        &
   &      7.175E-05,  7.134E-05,  7.094E-05,  7.014E-05,  7.014E-05,  &
   &      7.016E-05,  7.013E-05,  6.937E-05,  6.938E-05,  6.900E-05,  &
   &      6.863E-05,  6.860E-05,  6.859E-05,  6.741E-05,  6.700E-05,  &
   &      6.700E-05,  6.740E-05,  6.742E-05,  6.704E-05,  6.665E-05,  &
   &      6.549E-05,  6.432E-05,  6.427E-05,  6.463E-05,  6.581E-05,  &
   &      6.662E-05,  6.625E-05,  6.552E-05,  6.473E-05,  6.473E-05,  &
   &      6.471E-05,  6.469E-05,  6.391E-05,  6.433E-05,  6.475E-05,  &
   &      6.513E-05,  6.554E-05,  6.552E-05,  6.592E-05,  6.552E-05,  &
   &      6.553E-05,  6.476E-05,  6.440E-05,  6.483E-05,  6.443E-05,  &
   &      6.443E-05,  6.522E-05,  6.522E-05,  6.599E-05,  6.599E-05/
   DATA (Y(I),I= 1151, 1200)/                                        &
   &      6.560E-05,  6.520E-05,  6.521E-05,  6.558E-05,  6.556E-05,  &
   &      6.557E-05,  6.596E-05,  6.596E-05,  6.596E-05,  6.634E-05,  &
   &      6.595E-05,  6.555E-05,  6.555E-05,  6.594E-05,  6.515E-05,  &
   &      6.553E-05,  6.590E-05,  6.551E-05,  6.474E-05,  6.549E-05,  &
   &      6.589E-05,  6.628E-05,  6.709E-05,  6.672E-05,  6.634E-05,  &
   &      6.555E-05,  6.478E-05,  6.475E-05,  6.395E-05,  6.472E-05,  &
   &      6.472E-05,  6.510E-05,  6.509E-05,  6.508E-05,  6.511E-05,  &
   &      6.434E-05,  6.357E-05,  6.356E-05,  6.356E-05,  6.278E-05,  &
   &      6.239E-05,  6.161E-05,  6.199E-05,  6.163E-05,  6.161E-05,  &
   &      6.161E-05,  6.085E-05,  5.931E-05,  5.888E-05,  5.888E-05/
   DATA (Y(I),I= 1201, 1250)/                                        &
   &      5.846E-05,  5.805E-05,  5.809E-05,  5.890E-05,  5.892E-05,  &
   &      5.896E-05,  5.823E-05,  5.746E-05,  5.709E-05,  5.708E-05,  &
   &      5.707E-05,  5.668E-05,  5.709E-05,  5.666E-05,  5.627E-05,  &
   &      5.628E-05,  5.552E-05,  5.550E-05,  5.552E-05,  5.588E-05,  &
   &      5.627E-05,  5.588E-05,  5.627E-05,  5.591E-05,  5.477E-05,  &
   &      5.515E-05,  5.517E-05,  5.558E-05,  5.599E-05,  5.591E-05,  &
   &      5.668E-05,  5.507E-05,  5.429E-05,  5.311E-05,  5.279E-05,  &
   &      5.281E-05,  5.402E-05,  5.406E-05,  5.485E-05,  5.562E-05,  &
   &      5.562E-05,  5.481E-05,  5.484E-05,  5.483E-05,  5.443E-05,  &
   &      5.404E-05,  5.365E-05,  5.442E-05,  5.444E-05,  5.521E-05/
   DATA (Y(I),I= 1251, 1300)/                                        &
   &      5.481E-05,  5.443E-05,  5.363E-05,  5.286E-05,  5.283E-05,  &
   &      5.358E-05,  5.356E-05,  5.317E-05,  5.240E-05,  5.279E-05,  &
   &      5.240E-05,  5.241E-05,  5.321E-05,  5.322E-05,  5.327E-05,  &
   &      5.329E-05,  5.331E-05,  5.329E-05,  5.290E-05,  5.251E-05,  &
   &      5.213E-05,  5.135E-05,  5.057E-05,  5.056E-05,  5.056E-05,  &
   &      5.094E-05,  5.133E-05,  5.133E-05,  5.132E-05,  5.133E-05,  &
   &      5.054E-05,  4.934E-05,  4.815E-05,  4.733E-05,  4.735E-05,  &
   &      4.777E-05,  4.819E-05,  4.784E-05,  4.825E-05,  4.790E-05,  &
   &      4.830E-05,  4.793E-05,  4.873E-05,  4.874E-05,  4.790E-05,  &
   &      4.709E-05,  4.627E-05,  4.467E-05,  4.427E-05,  4.469E-05/
   DATA (Y(I),I= 1301, 1350)/                                        &
   &      4.434E-05,  4.431E-05,  4.473E-05,  4.477E-05,  4.365E-05,  &
   &      4.369E-05,  4.527E-05,  4.528E-05,  4.603E-05,  4.641E-05,  &
   &      4.520E-05,  4.480E-05,  4.363E-05,  4.404E-05,  4.444E-05,  &
   &      4.408E-05,  4.448E-05,  4.452E-05,  4.336E-05,  4.184E-05,  &
   &      4.223E-05,  4.263E-05,  4.301E-05,  4.262E-05,  4.299E-05,  &
   &      4.143E-05,  4.105E-05,  3.989E-05,  3.995E-05,  3.919E-05,  &
   &      3.999E-05,  4.118E-05,  4.039E-05,  3.997E-05,  3.880E-05,  &
   &      3.762E-05,  3.645E-05,  3.493E-05,  3.337E-05,  3.222E-05,  &
   &      3.296E-05,  3.412E-05,  3.255E-05,  3.257E-05,  3.103E-05,  &
   &      3.029E-05,  2.876E-05,  2.878E-05,  2.800E-05,  2.883E-05/
   DATA (Y(I),I= 1351, 1400)/                                        &
   &      2.881E-05,  2.806E-05,  2.729E-05,  2.650E-05,  2.574E-05,  &
   &      2.380E-05,  2.262E-05,  2.108E-05,  2.031E-05,  1.842E-05,  &
   &      1.765E-05,  1.648E-05,  1.646E-05,  1.685E-05,  1.529E-05,  &
   &      1.451E-05,  1.258E-05,  1.104E-05,  9.506E-06,  9.546E-06,  &
   &      8.010E-06,  6.431E-06,  4.851E-06,  4.067E-06,  2.472E-06,  &
   &      8.919E-07, -2.698E-07, -2.356E-07, -6.024E-07, -1.335E-06,  &
   &     -2.450E-06, -3.996E-06, -5.582E-06, -6.779E-06, -7.956E-06,  &
   &     -9.542E-06, -1.072E-05, -1.150E-05, -1.227E-05, -1.305E-05,  &
   &     -1.382E-05, -1.500E-05, -1.580E-05, -1.738E-05, -1.779E-05,  &
   &     -1.823E-05, -1.900E-05, -2.050E-05, -2.086E-05, -2.159E-05/
   DATA (Y(I),I= 1401, 1450)/                                        &
   &     -2.191E-05, -2.268E-05, -2.308E-05, -2.276E-05, -2.393E-05,  &
   &     -2.551E-05, -2.669E-05, -2.709E-05, -2.632E-05, -2.709E-05,  &
   &     -2.709E-05, -2.749E-05, -3.016E-05, -2.709E-05, -2.709E-05,  &
   &     -2.709E-05, -2.708E-05, -2.709E-05, -2.709E-05, -2.709E-05,  &
   &     -2.709E-05, -2.709E-05, -2.709E-05, -2.709E-05, -2.709E-05,  &
   &     -3.113E-05, -2.709E-05, -2.709E-05, -3.477E-05, -2.709E-05,  &
   &     -2.685E-05, -3.453E-05, -2.684E-05, -2.685E-05, -3.065E-05,  &
   &     -2.684E-05, -2.684E-05, -3.065E-05, -2.685E-05, -2.685E-05,  &
   &     -2.685E-05, -2.684E-05, -2.684E-05, -3.065E-05, -2.685E-05,  &
   &     -3.453E-05, -2.685E-05, -2.685E-05, -2.684E-05, -2.685E-05/
   DATA (Y(I),I= 1451, 1500)/                                        &
   &     -2.684E-05, -2.684E-05, -2.685E-05, -3.086E-05, -3.466E-05,  &
   &     -3.467E-05, -3.453E-05, -3.854E-05, -3.854E-05, -3.072E-05,  &
   &     -3.854E-05, -3.840E-05, -3.854E-05, -3.840E-05, -3.452E-05,  &
   &     -4.221E-05, -3.840E-05, -4.221E-05, -4.221E-05, -3.833E-05,  &
   &     -4.220E-05, -4.221E-05, -3.833E-05, -4.221E-05, -4.220E-05,  &
   &     -4.221E-05, -4.221E-05, -3.833E-05, -3.833E-05, -3.833E-05,  &
   &     -3.834E-05, -4.601E-05, -4.601E-05, -3.833E-05, -3.833E-05,  &
   &     -4.221E-05, -4.221E-05, -4.601E-05, -4.221E-05, -4.221E-05,  &
   &     -4.221E-05, -4.626E-05, -4.625E-05, -4.650E-05, -5.054E-05,  &
   &     -4.650E-05, -5.054E-05, -5.040E-05, -4.650E-05, -5.041E-05/
   DATA (Y(I),I= 1501, 1550)/                                        &
   &     -5.027E-05, -4.637E-05, -4.636E-05, -5.041E-05, -5.026E-05,  &
   &     -5.405E-05, -5.404E-05, -5.027E-05, -5.041E-05, -5.041E-05,  &
   &     -4.637E-05, -5.027E-05, -5.027E-05, -5.432E-05, -5.028E-05,  &
   &     -5.419E-05, -6.200E-05, -5.418E-05, -6.224E-05, -5.418E-05,  &
   &     -5.443E-05, -5.442E-05, -6.248E-05, -6.248E-05, -5.480E-05,  &
   &     -5.480E-05, -6.249E-05, -6.248E-05, -5.443E-05, -5.442E-05,  &
   &     -5.834E-05, -6.235E-05, -6.222E-05, -6.221E-05, -5.429E-05,  &
   &     -5.819E-05, -5.452E-05, -5.429E-05, -5.429E-05, -6.220E-05,  &
   &     -5.453E-05, -5.453E-05, -5.452E-05, -5.453E-05, -5.453E-05,  &
   &     -6.246E-05, -5.477E-05, -5.845E-05, -5.477E-05, -5.477E-05/
   DATA (Y(I),I= 1551, 1600)/                                        &
   &     -5.478E-05, -5.453E-05, -4.671E-05, -4.672E-05, -4.633E-05,  &
   &     -4.633E-05, -4.647E-05, -4.671E-05, -4.633E-05, -4.647E-05,  &
   &     -4.269E-05, -4.672E-05, -4.671E-05, -4.671E-05, -4.671E-05,  &
   &     -4.671E-05, -4.671E-05, -4.695E-05, -4.291E-05, -4.304E-05,  &
   &     -4.266E-05, -4.280E-05, -3.899E-05, -3.512E-05, -3.875E-05,  &
   &     -3.488E-05, -3.107E-05, -2.720E-05, -2.720E-05, -2.719E-05,  &
   &     -2.720E-05, -2.719E-05, -2.720E-05, -2.719E-05, -3.121E-05,  &
   &     -3.120E-05, -3.121E-05, -2.353E-05, -2.755E-05, -1.961E-05,  &
   &     -2.730E-05, -2.730E-05, -2.328E-05, -1.156E-05, -1.169E-05,  &
   &     -1.561E-05, -1.170E-05, -1.938E-05, -1.951E-05, -1.975E-05/
   DATA (Y(I),I= 1601, 1650)/                                        &
   &     -1.574E-05, -1.975E-05, -1.951E-05, -1.951E-05, -1.183E-05,  &
   &     -1.183E-05, -7.924E-06, -7.916E-06, -7.921E-06, -7.921E-06,  &
   &     -7.924E-06, -7.925E-06, -1.194E-05, -7.919E-06, -1.169E-05,  &
   &     -7.920E-06, -4.014E-06, -1.193E-05, -8.024E-06, -7.778E-06,  &
   &     -7.778E-06, -7.778E-06, -7.778E-06, -7.540E-06, -1.145E-05,  &
   &     -1.145E-05, -1.145E-05, -3.770E-06, -3.770E-06, -7.538E-06,  &
   &     -3.765E-06, -1.145E-05, -7.543E-06, -3.772E-06, -3.772E-06,  &
   &     -3.768E-06, -3.769E-06, -3.769E-06, -3.769E-06, -3.767E-06,  &
   &     -3.766E-06, -3.769E-06, -3.772E-06, -1.159E-05, -1.145E-05,  &
   &     -7.817E-06, -7.815E-06, -7.816E-06, -1.563E-05, -7.577E-06/
   DATA (Y(I),I= 1651, 1700)/                                        &
   &     -7.570E-06, -1.525E-05, -7.328E-06, -7.572E-06, -7.574E-06,  &
   &     -1.525E-05, -1.120E-05, -1.145E-05, -1.145E-05, -1.145E-05,  &
   &     -1.145E-05, -1.913E-05, -1.145E-05, -1.145E-05, -1.145E-05,  &
   &     -1.145E-05, -1.145E-05, -1.145E-05, -1.145E-05, -2.304E-05,  &
   &     -2.304E-05, -2.304E-05, -2.304E-05, -1.536E-05, -1.536E-05,  &
   &     -1.536E-05, -2.304E-05, -2.304E-05, -1.924E-05, -1.923E-05,  &
   &     -2.692E-05, -2.692E-05, -1.924E-05, -2.691E-05, -1.923E-05,  &
   &     -2.716E-05, -2.328E-05, -2.716E-05, -3.096E-05, -3.096E-05,  &
   &     -3.096E-05, -3.120E-05, -3.120E-05, -3.107E-05, -2.301E-05,  &
   &     -3.107E-05, -2.315E-05, -2.315E-05, -2.315E-05, -2.301E-05/
   DATA (Y(I),I= 1701, 1750)/                                        &
   &     -2.301E-05, -3.069E-05, -2.301E-05, -3.069E-05, -2.301E-05,  &
   &     -3.093E-05, -3.055E-05, -3.056E-05, -2.287E-05, -2.301E-05,  &
   &     -3.069E-05, -3.093E-05, -3.093E-05, -3.069E-05, -2.301E-05,  &
   &     -2.301E-05, -2.301E-05, -2.301E-05, -3.107E-05, -3.106E-05,  &
   &     -2.301E-05, -2.314E-05, -3.083E-05, -3.083E-05, -3.107E-05,  &
   &     -2.706E-05, -3.083E-05, -2.681E-05, -2.301E-05, -2.682E-05,  &
   &     -2.706E-05, -1.913E-05, -2.682E-05, -1.914E-05, -2.277E-05,  &
   &     -1.913E-05, -1.522E-05, -1.509E-05, -1.522E-05, -2.290E-05,  &
   &     -1.509E-05, -2.290E-05, -2.290E-05, -1.522E-05, -1.522E-05,  &
   &     -1.910E-05, -1.522E-05, -1.522E-05, -1.155E-05, -1.155E-05/
   DATA (Y(I),I= 1751, 1800)/                                        &
   &     -1.156E-05, -1.142E-05, -7.541E-06, -7.536E-06, -7.544E-06,  &
   &     -1.142E-05, -1.522E-05, -7.542E-06, -7.538E-06, -7.539E-06,  &
   &     -7.541E-06, -7.541E-06, -7.545E-06, -7.538E-06, -3.636E-06,  &
   &      1.426E-07, -3.494E-06, -3.633E-06, -3.630E-06, -3.631E-06,  &
   &     -3.634E-06,  4.149E-07, -3.629E-06, -3.629E-06, -3.633E-06,  &
   &     -3.632E-06, -3.634E-06, -3.628E-06, -3.385E-06,  4.292E-06,  &
   &      2.437E-07,  2.431E-07,  4.155E-06, -3.766E-06, -3.771E-06,  &
   &      2.430E-07,  2.438E-07,  3.915E-06,  3.909E-06,  7.824E-06,  &
   &      3.913E-06,  3.911E-06,  3.913E-06,  3.909E-06,  3.910E-06,  &
   &      3.910E-06,  3.914E-06,  4.150E-06,  1.184E-05,  1.184E-05/
   DATA (Y(I),I= 1801, 1850)/                                        &
   &      1.183E-05,  1.108E-05,  1.031E-05,  9.561E-06,  9.950E-06,  &
   &      1.227E-05,  1.420E-05,  1.534E-05,  1.572E-05,  1.571E-05,  &
   &      1.492E-05,  1.490E-05,  1.445E-05,  1.442E-05,  1.525E-05,  &
   &      1.487E-05,  1.452E-05,  1.604E-05,  1.840E-05,  2.070E-05,  &
   &      2.146E-05,  2.144E-05,  2.105E-05,  2.104E-05,  2.101E-05,  &
   &      2.100E-05,  2.175E-05,  2.095E-05,  2.097E-05,  1.981E-05,  &
   &      1.982E-05,  2.058E-05,  2.095E-05,  2.098E-05,  1.979E-05,  &
   &      1.860E-05,  1.782E-05,  1.934E-05,  1.972E-05,  2.083E-05,  &
   &      2.201E-05,  2.281E-05,  2.476E-05,  2.480E-05,  2.480E-05,  &
   &      2.407E-05,  2.371E-05,  2.295E-05,  2.218E-05,  2.141E-05/
   DATA (Y(I),I= 1851, 1900)/                                        &
   &      2.218E-05,  2.294E-05,  2.371E-05,  2.369E-05,  2.328E-05,  &
   &      2.287E-05,  2.169E-05,  2.128E-05,  2.126E-05,  2.004E-05,  &
   &      2.004E-05,  2.085E-05,  2.045E-05,  2.045E-05,  2.006E-05,  &
   &      2.085E-05,  2.082E-05,  2.121E-05,  2.198E-05,  2.275E-05,  &
   &      2.238E-05,  2.083E-05,  2.046E-05,  1.928E-05,  1.893E-05,  &
   &      1.815E-05,  1.815E-05,  1.812E-05,  1.929E-05,  2.046E-05,  &
   &      2.087E-05,  2.244E-05,  2.283E-05,  2.322E-05,  2.327E-05,  &
   &      2.291E-05,  2.329E-05,  2.482E-05,  2.558E-05,  2.593E-05,  &
   &      2.514E-05,  2.436E-05,  2.279E-05,  2.282E-05,  2.363E-05,  &
   &      2.441E-05,  2.524E-05,  2.599E-05,  2.759E-05,  2.874E-05/
   DATA (Y(I),I= 1901, 1950)/                                        &
   &      2.954E-05,  3.028E-05,  3.064E-05,  3.063E-05,  3.021E-05,  &
   &      2.866E-05,  2.788E-05,  2.788E-05,  2.749E-05,  2.711E-05,  &
   &      2.866E-05,  3.016E-05,  3.250E-05,  3.323E-05,  3.366E-05,  &
   &      3.328E-05,  3.367E-05,  3.366E-05,  3.483E-05,  3.521E-05,  &
   &      3.599E-05,  3.603E-05,  3.606E-05,  3.685E-05,  3.764E-05,  &
   &      3.960E-05,  4.076E-05,  4.079E-05,  4.041E-05,  3.926E-05,  &
   &      3.848E-05,  4.004E-05,  4.161E-05,  4.350E-05,  4.426E-05,  &
   &      4.386E-05,  4.305E-05,  4.344E-05,  4.457E-05,  4.570E-05,  &
   &      4.685E-05,  4.567E-05,  4.450E-05,  4.337E-05,  4.297E-05,  &
   &      4.414E-05,  4.527E-05,  4.644E-05,  4.721E-05,  4.798E-05/
   DATA (Y(I),I= 1951, 2000)/                                        &
   &      4.873E-05,  4.952E-05,  4.797E-05,  4.798E-05,  4.681E-05,  &
   &      4.721E-05,  4.760E-05,  4.724E-05,  4.802E-05,  4.805E-05,  &
   &      4.766E-05,  4.690E-05,  4.727E-05,  4.800E-05,  4.911E-05,  &
   &      4.947E-05,  4.636E-05,  4.401E-05,  4.164E-05,  4.089E-05,  &
   &      4.167E-05,  4.245E-05,  4.248E-05,  4.136E-05,  3.944E-05,  &
   &      3.832E-05,  3.753E-05,  3.867E-05,  3.982E-05,  4.022E-05,  &
   &      3.829E-05,  3.714E-05,  3.481E-05,  3.365E-05,  3.438E-05,  &
   &      3.478E-05,  3.474E-05,  3.395E-05,  3.277E-05,  3.085E-05,  &
   &      2.966E-05,  3.045E-05,  3.044E-05,  3.123E-05,  3.042E-05,  &
   &      2.965E-05,  2.969E-05,  2.966E-05,  3.082E-05,  3.273E-05/
   DATA (Y(I),I= 2001, 2050)/                                        &
   &      3.388E-05,  3.310E-05,  3.195E-05,  3.196E-05,  3.080E-05,  &
   &      3.004E-05,  3.004E-05,  2.927E-05,  2.927E-05,  3.044E-05,  &
   &      3.081E-05,  3.198E-05,  3.197E-05,  3.237E-05,  3.200E-05,  &
   &      3.236E-05,  3.277E-05,  3.239E-05,  3.238E-05,  3.276E-05,  &
   &      3.199E-05,  3.161E-05,  3.200E-05,  3.199E-05,  3.236E-05,  &
   &      3.352E-05,  3.351E-05,  3.429E-05,  3.388E-05,  3.350E-05,  &
   &      3.426E-05,  3.541E-05,  3.618E-05,  3.697E-05,  3.697E-05,  &
   &      3.696E-05,  3.737E-05,  3.737E-05,  3.736E-05,  3.854E-05,  &
   &      3.893E-05,  3.738E-05,  3.624E-05,  3.395E-05,  3.318E-05,  &
   &      3.550E-05,  3.819E-05,  4.092E-05,  4.204E-05,  4.245E-05/
   DATA (Y(I),I= 2051, 2100)/                                        &
   &      4.284E-05,  4.325E-05,  4.324E-05,  4.363E-05,  4.403E-05,  &
   &      4.442E-05,  4.365E-05,  4.284E-05,  4.246E-05,  4.244E-05,  &
   &      4.360E-05,  4.552E-05,  4.629E-05,  4.548E-05,  4.431E-05,  &
   &      4.316E-05,  4.238E-05,  4.353E-05,  4.428E-05,  4.504E-05,  &
   &      4.581E-05,  4.621E-05,  4.543E-05,  4.583E-05,  4.620E-05,  &
   &      4.659E-05,  4.622E-05,  4.661E-05,  4.621E-05,  4.583E-05,  &
   &      4.583E-05,  4.503E-05,  4.582E-05,  4.541E-05,  4.578E-05,  &
   &      4.577E-05,  4.655E-05,  4.691E-05,  4.769E-05,  4.693E-05,  &
   &      4.693E-05,  4.539E-05,  4.462E-05,  4.462E-05,  4.425E-05,  &
   &      4.351E-05,  4.391E-05,  4.468E-05,  4.588E-05,  4.741E-05/
   DATA (Y(I),I= 2101, 2150)/                                        &
   &      4.818E-05,  5.011E-05,  5.046E-05,  5.198E-05,  5.124E-05,  &
   &      5.086E-05,  4.931E-05,  4.776E-05,  4.700E-05,  4.623E-05,  &
   &      4.469E-05,  4.472E-05,  4.625E-05,  4.856E-05,  5.087E-05,  &
   &      5.239E-05,  5.318E-05,  5.474E-05,  5.547E-05,  5.626E-05,  &
   &      5.706E-05,  5.668E-05,  5.704E-05,  5.629E-05,  5.361E-05,  &
   &      5.208E-05,  5.017E-05,  5.209E-05,  5.554E-05,  5.864E-05,  &
   &      6.132E-05,  6.251E-05,  6.253E-05,  6.293E-05,  6.371E-05,  &
   &      6.489E-05,  6.606E-05,  6.608E-05,  6.644E-05,  6.453E-05,  &
   &      6.336E-05,  6.178E-05,  6.177E-05,  6.327E-05,  6.557E-05,  &
   &      6.711E-05,  6.864E-05,  6.901E-05,  6.937E-05,  6.975E-05/
   DATA (Y(I),I= 2151, 2200)/                                        &
   &      6.821E-05,  6.704E-05,  6.554E-05,  6.438E-05,  6.398E-05,  &
   &      6.318E-05,  6.319E-05,  6.280E-05,  6.048E-05,  5.933E-05,  &
   &      5.817E-05,  5.660E-05,  5.620E-05,  5.619E-05,  5.616E-05,  &
   &      5.501E-05,  5.226E-05,  4.917E-05,  4.647E-05,  4.528E-05,  &
   &      4.565E-05,  4.602E-05,  4.641E-05,  4.524E-05,  4.213E-05,  &
   &      3.943E-05,  3.713E-05,  3.711E-05,  3.866E-05,  3.982E-05,  &
   &      4.137E-05,  3.944E-05,  3.674E-05,  3.368E-05,  3.139E-05,  &
   &      3.100E-05,  3.139E-05,  3.212E-05,  3.251E-05,  3.135E-05,  &
   &      2.904E-05,  2.633E-05,  2.441E-05,  2.478E-05,  2.480E-05,  &
   &      2.557E-05,  2.594E-05,  2.477E-05,  2.404E-05,  2.404E-05/
   DATA (Y(I),I= 2201, 2250)/                                        &
   &      2.325E-05,  2.365E-05,  2.286E-05,  2.365E-05,  2.442E-05,  &
   &      2.404E-05,  2.406E-05,  2.401E-05,  2.403E-05,  2.480E-05,  &
   &      2.480E-05,  2.598E-05,  2.675E-05,  2.675E-05,  2.793E-05,  &
   &      2.793E-05,  2.793E-05,  2.796E-05,  2.758E-05,  2.761E-05,  &
   &      2.721E-05,  2.875E-05,  2.991E-05,  3.147E-05,  3.301E-05,  &
   &      3.266E-05,  3.187E-05,  3.148E-05,  3.110E-05,  3.261E-05,  &
   &      3.417E-05,  3.532E-05,  3.686E-05,  3.684E-05,  3.606E-05,  &
   &      3.568E-05,  3.493E-05,  3.645E-05,  3.797E-05,  4.030E-05,  &
   &      4.144E-05,  4.221E-05,  4.181E-05,  4.184E-05,  4.142E-05,  &
   &      4.181E-05,  4.257E-05,  4.372E-05,  4.412E-05,  4.334E-05/
   DATA (Y(I),I= 2251, 2300)/                                        &
   &      4.220E-05,  4.146E-05,  3.991E-05,  4.067E-05,  4.223E-05,  &
   &      4.377E-05,  4.608E-05,  4.609E-05,  4.608E-05,  4.571E-05,  &
   &      4.569E-05,  4.532E-05,  4.608E-05,  4.648E-05,  4.684E-05,  &
   &      4.764E-05,  4.725E-05,  4.803E-05,  4.805E-05,  4.726E-05,  &
   &      4.806E-05,  4.806E-05,  4.810E-05,  4.770E-05,  4.731E-05,  &
   &      4.617E-05,  4.577E-05,  4.424E-05,  4.464E-05,  4.502E-05,  &
   &      4.577E-05,  4.538E-05,  4.619E-05,  4.696E-05,  4.696E-05,  &
   &      4.696E-05,  4.655E-05,  4.657E-05,  4.540E-05,  4.537E-05,  &
   &      4.618E-05,  4.770E-05,  4.847E-05,  5.041E-05,  5.041E-05,  &
   &      4.928E-05,  4.773E-05,  4.619E-05,  4.502E-05,  4.582E-05/
   DATA (Y(I),I= 2301, 2350)/                                        &
   &      4.658E-05,  4.581E-05,  4.658E-05,  4.734E-05,  4.813E-05,  &
   &      4.809E-05,  4.927E-05,  4.848E-05,  4.814E-05,  4.813E-05,  &
   &      4.736E-05,  4.931E-05,  4.968E-05,  5.124E-05,  5.317E-05,  &
   &      5.242E-05,  5.087E-05,  5.013E-05,  4.899E-05,  4.819E-05,  &
   &      5.129E-05,  5.321E-05,  5.514E-05,  5.667E-05,  5.666E-05,  &
   &      5.669E-05,  5.745E-05,  5.708E-05,  5.744E-05,  5.746E-05,  &
   &      5.823E-05,  5.862E-05,  5.745E-05,  5.706E-05,  5.593E-05,  &
   &      5.513E-05,  5.512E-05,  5.512E-05,  5.586E-05,  5.663E-05,  &
   &      5.701E-05,  5.662E-05,  5.619E-05,  5.657E-05,  5.618E-05,  &
   &      5.502E-05,  5.385E-05,  5.307E-05,  5.191E-05,  5.113E-05/
   DATA (Y(I),I= 2351, 2400)/                                        &
   &      5.074E-05,  5.033E-05,  5.036E-05,  4.994E-05,  4.801E-05,  &
   &      4.723E-05,  4.609E-05,  4.531E-05,  4.724E-05,  4.992E-05,  &
   &      5.185E-05,  5.298E-05,  5.028E-05,  4.607E-05,  4.225E-05,  &
   &      3.800E-05,  3.608E-05,  3.610E-05,  3.454E-05,  3.455E-05,  &
   &      3.339E-05,  3.302E-05,  3.263E-05,  3.264E-05,  3.186E-05,  &
   &      3.262E-05,  3.417E-05,  3.493E-05,  3.647E-05,  3.530E-05,  &
   &      3.375E-05,  3.221E-05,  3.068E-05,  3.031E-05,  2.992E-05,  &
   &      3.029E-05,  3.029E-05,  3.108E-05,  3.108E-05,  3.031E-05,  &
   &      3.031E-05,  3.031E-05,  3.031E-05,  3.031E-05,  3.070E-05,  &
   &      2.994E-05,  3.111E-05,  3.109E-05,  3.184E-05,  3.300E-05/
   DATA (Y(I),I= 2401, 2450)/                                        &
   &      3.377E-05,  3.379E-05,  3.378E-05,  3.378E-05,  3.380E-05,  &
   &      3.456E-05,  3.535E-05,  3.572E-05,  3.651E-05,  3.650E-05,  &
   &      3.766E-05,  3.766E-05,  3.843E-05,  3.959E-05,  3.961E-05,  &
   &      4.037E-05,  4.039E-05,  4.038E-05,  4.077E-05,  3.923E-05,  &
   &      3.806E-05,  3.654E-05,  3.537E-05,  3.690E-05,  3.804E-05,  &
   &      3.956E-05,  4.189E-05,  4.264E-05,  4.228E-05,  4.150E-05,  &
   &      4.113E-05,  4.152E-05,  4.189E-05,  4.189E-05,  4.265E-05,  &
   &      4.304E-05,  4.227E-05,  4.189E-05,  4.035E-05,  3.921E-05,  &
   &      3.844E-05,  3.920E-05,  3.996E-05,  4.072E-05,  4.265E-05,  &
   &      4.265E-05,  4.265E-05,  4.189E-05,  4.148E-05,  4.151E-05/
   DATA (Y(I),I= 2451, 2500)/                                        &
   &      4.151E-05,  4.188E-05,  4.188E-05,  4.188E-05,  4.306E-05,  &
   &      4.226E-05,  4.304E-05,  4.385E-05,  4.357E-05,  4.342E-05,  &
   &      4.250E-05,  4.239E-05,  4.224E-05,  4.160E-05,  4.217E-05,  &
   &      4.281E-05,  4.298E-05,  4.359E-05,  4.302E-05,  4.287E-05,  &
   &      4.195E-05,  4.139E-05,  4.195E-05,  4.233E-05,  4.270E-05,  &
   &      4.229E-05,  4.268E-05,  4.234E-05,  4.191E-05,  4.138E-05,  &
   &      4.084E-05,  4.065E-05,  4.062E-05,  4.052E-05,  4.045E-05,  &
   &      4.053E-05,  4.100E-05,  4.151E-05,  4.206E-05,  4.260E-05,  &
   &      4.292E-05,  4.304E-05,  4.320E-05,  4.333E-05,  4.345E-05,  &
   &      4.377E-05,  4.404E-05,  4.432E-05,  4.460E-05,  4.453E-05/
   DATA (Y(I),I= 2501, 2550)/                                        &
   &      4.380E-05,  4.316E-05,  4.239E-05,  4.170E-05,  4.214E-05,  &
   &      4.302E-05,  4.387E-05,  4.472E-05,  4.523E-05,  4.454E-05,  &
   &      4.369E-05,  4.300E-05,  4.215E-05,  4.239E-05,  4.346E-05,  &
   &      4.450E-05,  4.554E-05,  4.650E-05,  4.638E-05,  4.595E-05,  &
   &      4.568E-05,  4.529E-05,  4.501E-05,  4.520E-05,  4.527E-05,  &
   &      4.534E-05,  4.541E-05,  4.559E-05,  4.584E-05,  4.613E-05,  &
   &      4.626E-05,  4.628E-05,  4.530E-05,  4.413E-05,  4.292E-05,  &
   &      4.171E-05,  4.101E-05,  4.133E-05,  4.153E-05,  4.180E-05,  &
   &      4.208E-05,  4.158E-05,  4.078E-05,  4.003E-05,  3.923E-05,  &
   &      3.854E-05,  3.888E-05,  3.936E-05,  3.993E-05,  4.038E-05/
   DATA (Y(I),I= 2551, 2600)/                                        &
   &      4.018E-05,  3.852E-05,  3.686E-05,  3.515E-05,  3.349E-05,  &
   &      3.280E-05,  3.284E-05,  3.285E-05,  3.277E-05,  3.273E-05,  &
   &      3.205E-05,  3.128E-05,  3.044E-05,  2.968E-05,  2.902E-05,  &
   &      2.895E-05,  2.891E-05,  2.895E-05,  2.895E-05,  2.876E-05,  &
   &      2.837E-05,  2.803E-05,  2.760E-05,  2.726E-05,  2.707E-05,  &
   &      2.707E-05,  2.700E-05,  2.700E-05,  2.689E-05,  2.685E-05,  &
   &      2.681E-05,  2.674E-05,  2.667E-05,  2.651E-05,  2.628E-05,  &
   &      2.602E-05,  2.575E-05,  2.552E-05,  2.544E-05,  2.567E-05,  &
   &      2.591E-05,  2.614E-05,  2.633E-05,  2.610E-05,  2.557E-05,  &
   &      2.511E-05,  2.457E-05,  2.416E-05,  2.450E-05,  2.512E-05/
   DATA (Y(I),I= 2601, 2650)/                                        &
   &      2.578E-05,  2.635E-05,  2.678E-05,  2.655E-05,  2.628E-05,  &
   &      2.594E-05,  2.563E-05,  2.548E-05,  2.606E-05,  2.667E-05,  &
   &      2.725E-05,  2.779E-05,  2.802E-05,  2.768E-05,  2.737E-05,  &
   &      2.699E-05,  2.661E-05,  2.681E-05,  2.746E-05,  2.816E-05,  &
   &      2.877E-05,  2.939E-05,  2.951E-05,  2.935E-05,  2.924E-05,  &
   &      2.913E-05,  2.897E-05,  2.916E-05,  2.936E-05,  2.951E-05,  &
   &      2.978E-05,  2.990E-05,  2.997E-05,  3.001E-05,  3.004E-05,  &
   &      3.004E-05,  3.000E-05,  2.981E-05,  2.957E-05,  2.938E-05,  &
   &      2.918E-05,  2.903E-05,  2.918E-05,  2.930E-05,  2.937E-05,  &
   &      2.948E-05,  2.944E-05,  2.898E-05,  2.848E-05,  2.806E-05/
   DATA (Y(I),I= 2651, 2700)/                                        &
   &      2.756E-05,  2.733E-05,  2.744E-05,  2.751E-05,  2.763E-05,  &
   &      2.767E-05,  2.751E-05,  2.716E-05,  2.682E-05,  2.635E-05,  &
   &      2.605E-05,  2.585E-05,  2.581E-05,  2.578E-05,  2.570E-05,  &
   &      2.574E-05,  2.566E-05,  2.554E-05,  2.551E-05,  2.543E-05,  &
   &      2.535E-05,  2.535E-05,  2.544E-05,  2.548E-05,  2.560E-05,  &
   &      2.564E-05,  2.569E-05,  2.585E-05,  2.589E-05,  2.597E-05,  &
   &      2.601E-05,  2.617E-05,  2.641E-05,  2.656E-05,  2.671E-05,  &
   &      2.695E-05,  2.710E-05,  2.730E-05,  2.746E-05,  2.769E-05,  &
   &      2.785E-05,  2.800E-05,  2.812E-05,  2.831E-05,  2.842E-05,  &
   &      2.858E-05,  2.878E-05,  2.889E-05,  2.908E-05,  2.920E-05/
   DATA (Y(I),I= 2701, 2750)/                                        &
   &      2.939E-05,  2.935E-05,  2.935E-05,  2.931E-05,  2.923E-05,  &
   &      2.927E-05,  2.919E-05,  2.915E-05,  2.915E-05,  2.911E-05,  &
   &      2.903E-05,  2.895E-05,  2.864E-05,  2.837E-05,  2.814E-05,  &
   &      2.790E-05,  2.763E-05,  2.740E-05,  2.716E-05,  2.693E-05,  &
   &      2.666E-05,  2.639E-05,  2.612E-05,  2.584E-05,  2.550E-05,  &
   &      2.526E-05,  2.499E-05,  2.464E-05,  2.437E-05,  2.414E-05,  &
   &      2.391E-05,  2.360E-05,  2.336E-05,  2.313E-05,  2.293E-05,  &
   &      2.266E-05,  2.247E-05,  2.227E-05,  2.196E-05,  2.177E-05,  &
   &      2.157E-05,  2.130E-05,  2.118E-05,  2.107E-05,  2.095E-05,  &
   &      2.075E-05,  2.063E-05,  2.044E-05,  2.036E-05,  2.017E-05/
   DATA (Y(I),I= 2751, 2800)/                                        &
   &      2.005E-05,  1.994E-05,  1.982E-05,  1.971E-05,  1.959E-05,  &
   &      1.944E-05,  1.933E-05,  1.925E-05,  1.906E-05,  1.898E-05,  &
   &      1.887E-05,  1.876E-05,  1.864E-05,  1.856E-05,  1.849E-05,  &
   &      1.830E-05,  1.827E-05,  1.815E-05,  1.807E-05,  1.800E-05,  &
   &      1.788E-05,  1.778E-05,  1.770E-05,  1.766E-05,  1.762E-05,  &
   &      1.755E-05,  1.755E-05,  1.751E-05,  1.740E-05,  1.736E-05,  &
   &      1.736E-05,  1.728E-05,  1.725E-05,  1.721E-05,  1.717E-05,  &
   &      1.712E-05,  1.712E-05,  1.709E-05,  1.704E-05,  1.693E-05,  &
   &      1.693E-05,  1.685E-05,  1.685E-05,  1.680E-05,  1.680E-05,  &
   &      1.672E-05,  1.672E-05,  1.672E-05,  1.672E-05,  1.668E-05/
   DATA (Y(I),I= 2801, 2850)/                                        &
   &      1.660E-05,  1.659E-05,  1.659E-05,  1.651E-05,  1.655E-05,  &
   &      1.655E-05,  1.655E-05,  1.651E-05,  1.647E-05,  1.647E-05,  &
   &      1.647E-05,  1.643E-05,  1.643E-05,  1.636E-05,  1.636E-05,  &
   &      1.640E-05,  1.640E-05,  1.636E-05,  1.640E-05,  1.640E-05,  &
   &      1.640E-05,  1.637E-05,  1.637E-05,  1.641E-05,  1.637E-05,  &
   &      1.641E-05,  1.638E-05,  1.626E-05,  1.622E-05,  1.619E-05,  &
   &      1.611E-05,  1.611E-05,  1.607E-05,  1.604E-05,  1.600E-05,  &
   &      1.596E-05,  1.589E-05,  1.573E-05,  1.562E-05,  1.550E-05,  &
   &      1.538E-05,  1.531E-05,  1.519E-05,  1.499E-05,  1.496E-05,  &
   &      1.480E-05,  1.469E-05,  1.461E-05,  1.445E-05,  1.437E-05/
   DATA (Y(I),I= 2851, 2900)/                                        &
   &      1.422E-05,  1.418E-05,  1.403E-05,  1.399E-05,  1.383E-05,  &
   &      1.375E-05,  1.360E-05,  1.356E-05,  1.356E-05,  1.360E-05,  &
   &      1.360E-05,  1.364E-05,  1.360E-05,  1.360E-05,  1.360E-05,  &
   &      1.360E-05,  1.356E-05,  1.364E-05,  1.364E-05,  1.368E-05,  &
   &      1.380E-05,  1.384E-05,  1.388E-05,  1.392E-05,  1.403E-05,  &
   &      1.403E-05,  1.411E-05,  1.411E-05,  1.419E-05,  1.423E-05,  &
   &      1.431E-05,  1.427E-05,  1.435E-05,  1.439E-05,  1.439E-05,  &
   &      1.443E-05,  1.451E-05,  1.454E-05,  1.454E-05,  1.458E-05,  &
   &      1.458E-05,  1.450E-05,  1.450E-05,  1.446E-05,  1.446E-05,  &
   &      1.446E-05,  1.438E-05,  1.442E-05,  1.441E-05,  1.433E-05/
   DATA (Y(I),I= 2901, 2950)/                                        &
   &      1.429E-05,  1.422E-05,  1.418E-05,  1.406E-05,  1.398E-05,  &
   &      1.387E-05,  1.379E-05,  1.367E-05,  1.360E-05,  1.352E-05,  &
   &      1.340E-05,  1.329E-05,  1.313E-05,  1.310E-05,  1.290E-05,  &
   &      1.275E-05,  1.263E-05,  1.248E-05,  1.233E-05,  1.225E-05,  &
   &      1.210E-05,  1.198E-05,  1.187E-05,  1.171E-05,  1.156E-05,  &
   &      1.144E-05,  1.133E-05,  1.118E-05,  1.106E-05,  1.098E-05,  &
   &      1.081E-05,  1.073E-05,  1.055E-05,  1.046E-05,  1.035E-05,  &
   &      1.025E-05,  1.013E-05,  1.003E-05,  9.912E-06,  9.800E-06,  &
   &      9.691E-06,  9.582E-06,  9.470E-06,  9.362E-06,  9.253E-06,  &
   &      9.175E-06,  9.097E-06,  9.027E-06,  8.950E-06,  8.883E-06/
   DATA (Y(I),I= 2951, 3000)/                                        &
   &      8.802E-06,  8.732E-06,  8.658E-06,  8.588E-06,  8.510E-06,  &
   &      8.440E-06,  8.385E-06,  8.342E-06,  8.307E-06,  8.267E-06,  &
   &      8.228E-06,  8.189E-06,  8.149E-06,  8.106E-06,  8.063E-06,  &
   &      8.024E-06,  7.984E-06,  7.957E-06,  7.949E-06,  7.944E-06,  &
   &      7.928E-06,  7.920E-06,  7.908E-06,  7.907E-06,  7.895E-06,  &
   &      7.883E-06,  7.874E-06,  7.870E-06,  7.854E-06,  7.858E-06,  &
   &      7.853E-06,  7.850E-06,  7.853E-06,  7.849E-06,  7.853E-06,  &
   &      7.849E-06,  7.845E-06,  7.852E-06,  7.848E-06,  7.844E-06,  &
   &      7.836E-06,  7.825E-06,  7.810E-06,  7.798E-06,  7.787E-06,  &
   &      7.776E-06,  7.753E-06,  7.741E-06,  7.726E-06,  7.715E-06/
   DATA (Y(I),I= 3001, 3050)/                                        &
   &      7.703E-06,  7.684E-06,  7.642E-06,  7.604E-06,  7.565E-06,  &
   &      7.527E-06,  7.492E-06,  7.454E-06,  7.416E-06,  7.370E-06,  &
   &      7.332E-06,  7.297E-06,  7.251E-06,  7.197E-06,  7.132E-06,  &
   &      7.070E-06,  7.005E-06,  6.935E-06,  6.870E-06,  6.805E-06,  &
   &      6.743E-06,  6.674E-06,  6.608E-06,  6.543E-06,  6.481E-06,  &
   &      6.393E-06,  6.312E-06,  6.223E-06,  6.143E-06,  6.058E-06,  &
   &      5.977E-06,  5.897E-06,  5.816E-06,  5.723E-06,  5.643E-06,  &
   &      5.558E-06,  5.493E-06,  5.420E-06,  5.358E-06,  5.285E-06,  &
   &      5.212E-06,  5.150E-06,  5.073E-06,  5.004E-06,  4.939E-06,  &
   &      4.866E-06,  4.804E-06,  4.743E-06,  4.731E-06,  4.727E-06/
   DATA (Y(I),I= 3051, 3100)/                                        &
   &      4.723E-06,  4.719E-06,  4.715E-06,  4.707E-06,  4.699E-06,  &
   &      4.699E-06,  4.695E-06,  4.691E-06,  4.680E-06,  4.695E-06,  &
   &      4.745E-06,  4.798E-06,  4.844E-06,  4.890E-06,  4.944E-06,  &
   &      4.986E-06,  5.040E-06,  5.086E-06,  5.136E-06,  5.189E-06,  &
   &      5.236E-06,  5.262E-06,  5.270E-06,  5.289E-06,  5.297E-06,  &
   &      5.305E-06,  5.324E-06,  5.332E-06,  5.347E-06,  5.355E-06,  &
   &      5.370E-06,  5.381E-06,  5.389E-06,  5.347E-06,  5.297E-06,  &
   &      5.231E-06,  5.181E-06,  5.127E-06,  5.065E-06,  5.011E-06,  &
   &      4.950E-06,  4.896E-06,  4.838E-06,  4.780E-06,  4.726E-06,  &
   &      4.661E-06,  4.595E-06,  4.525E-06,  4.456E-06,  4.390E-06/
   DATA (Y(I),I= 3101, 3150)/                                        &
   &      4.328E-06,  4.255E-06,  4.193E-06,  4.127E-06,  4.058E-06,  &
   &      3.988E-06,  3.926E-06,  3.849E-06,  3.760E-06,  3.679E-06,  &
   &      3.598E-06,  3.520E-06,  3.439E-06,  3.357E-06,  3.272E-06,  &
   &      3.191E-06,  3.118E-06,  3.033E-06,  2.944E-06,  2.824E-06,  &
   &      2.693E-06,  2.561E-06,  2.434E-06,  2.306E-06,  2.179E-06,  &
   &      2.043E-06,  1.917E-06,  1.790E-06,  1.663E-06,  1.533E-06,  &
   &      1.300E-06,  1.100E-06,  0.910E-06,  0.725E-06,  0.575E-06,  &
   &      0.450E-06,  0.370E-06,  0.295E-06,  0.225E-06,  0.160E-06,  &
   &      0.100E-06,  0.055E-06,  0.035E-06,  0.025E-06,  0.020E-06,  &
   &      0.015E-06,  0.010E-06,  0.006E-06,  0.003E-06,  0.00000  /
!
   DATA (Z(I),I=    1,   50)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I= 51, 100)/                                           &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I=  101,  150)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I= 151, 200)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I=  201,  250)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I= 251, 300)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I=  301,  350)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I= 351, 400)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I=  401,  450)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I= 451, 500)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I=  501,  550)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I= 551, 600)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I=  601,  650)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I= 651, 700)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I=  701,  750)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I= 751, 800)/                                          &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    /
   DATA (Z(I),I=  801,  850)/                                        &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &      0.000    ,  0.000    ,  0.000    ,  0.000    ,  0.000    ,  &
   &     -5.000e-10, -1.000e-09, -1.500e-09, -2.500e-09, -5.000e-09,  &
   &     -7.500e-09, -1.000e-08, -1.375e-08, -1.750e-08, -2.325e-08,  &
   &     -3.200e-08, -4.100e-08, -5.000e-08, -6.000e-08, -7.000e-08,  &
   &     -8.000e-08, -9.500e-08, -1.200e-07, -2.400e-07, -3.300e-07,  &
   &     -3.400e-07, -3.476E-07, -3.507E-07, -3.498E-07, -3.441E-07,  &
   &     -3.405E-07, -3.384E-07, -3.439E-07, -3.538E-07, -3.699E-07/
   DATA (Z(I),I= 851, 900)/                                          &
   &     -3.651E-07, -3.462E-07, -3.292E-07, -3.177E-07, -3.239E-07,  &
   &     -3.321E-07, -3.447E-07, -3.404E-07, -3.284E-07, -3.301E-07,  &
   &     -3.295E-07, -3.214E-07, -3.194E-07, -3.230E-07, -3.202E-07,  &
   &     -3.182E-07, -3.245E-07, -3.272E-07, -3.181E-07, -3.133E-07,  &
   &     -3.115E-07, -3.134E-07, -3.144E-07, -3.226E-07, -3.292E-07,  &
   &     -3.239E-07, -3.248E-07, -3.088E-07, -3.269E-07, -3.166E-07,  &
   &     -3.283E-07, -3.143E-07, -3.064E-07, -3.140E-07, -2.938E-07,  &
   &     -2.938E-07, -3.014E-07, -2.873E-07, -2.809E-07, -2.809E-07,  &
   &     -3.092E-07, -2.748E-07, -2.809E-07, -2.744E-07, -2.809E-07,  &
   &     -2.669E-07, -2.666E-07, -2.870E-07, -2.806E-07, -2.806E-07/
   DATA (Z(I),I=  901,  950)/                                        &
   &     -2.806E-07, -2.870E-07, -2.730E-07, -2.590E-07, -2.512E-07,  &
   &     -2.792E-07, -2.856E-07, -2.860E-07, -2.781E-07, -2.565E-07,  &
   &     -2.425E-07, -2.626E-07, -2.767E-07, -2.624E-07, -2.484E-07,  &
   &     -2.343E-07, -2.484E-07, -2.689E-07, -2.268E-07, -2.329E-07,  &
   &     -2.470E-07, -2.193E-07, -2.254E-07, -2.318E-07, -2.318E-07,  &
   &     -2.459E-07, -2.240E-07, -2.021E-07, -2.161E-07, -2.021E-07,  &
   &     -2.162E-07, -2.226E-07, -2.022E-07, -2.086E-07, -2.083E-07,  &
   &     -2.148E-07, -2.288E-07, -2.285E-07, -2.286E-07, -2.288E-07,  &
   &     -2.209E-07, -2.069E-07, -2.069E-07, -2.131E-07, -2.131E-07,  &
   &     -1.991E-07, -1.991E-07, -1.991E-07, -1.707E-07, -1.851E-07/
   DATA (Z(I),I= 951, 1000)/                                         &
   &     -1.772E-07, -1.775E-07, -1.912E-07, -1.976E-07, -2.117E-07,  &
   &     -1.974E-07, -1.772E-07, -1.837E-07, -1.693E-07, -1.693E-07,  &
   &     -1.769E-07, -1.693E-07, -1.632E-07, -1.489E-07, -1.352E-07,  &
   &     -1.492E-07, -1.495E-07, -1.556E-07, -1.492E-07, -1.352E-07,  &
   &     -1.352E-07, -1.492E-07, -1.492E-07, -1.349E-07, -1.285E-07,  &
   &     -1.147E-07, -1.287E-07, -1.150E-07, -1.147E-07, -1.007E-07,  &
   &     -9.455E-08, -8.022E-08, -6.618E-08, -7.990E-08, -8.016E-08,  &
   &     -7.992E-08, -8.779E-08, -7.986E-08, -1.015E-07, -9.567E-08,  &
   &     -7.373E-08, -8.135E-08, -9.567E-08, -7.520E-08, -8.955E-08,  &
   &     -7.521E-08, -1.033E-07, -8.133E-08, -8.895E-08, -7.493E-08/
   DATA (Z(I),I= 1001, 1050)/                                        &
   &     -5.446E-08, -6.850E-08, -9.655E-08, -6.846E-08, -9.040E-08,  &
   &     -7.640E-08, -7.669E-08, -6.239E-08, -5.600E-08, -7.612E-08,  &
   &     -6.998E-08, -4.985E-08, -5.597E-08, -7.032E-08, -7.032E-08,  &
   &     -7.031E-08, -6.243E-08, -5.453E-08, -5.458E-08, -4.817E-08,  &
   &     -5.461E-08, -4.059E-08, -6.867E-08, -4.671E-08, -5.316E-08,  &
   &     -3.889E-08, -3.915E-08, -3.888E-08, -3.131E-08, -1.732E-08,  &
   &     -1.729E-08, -3.134E-08, -1.735E-08, -9.496E-09, -9.835E-09,  &
   &     -1.596E-08, -3.031E-08, -3.644E-08, -4.288E-08, -2.243E-08,  &
   &     -1.488E-08, -8.478E-10,  2.750E-08,  3.508E-08,  3.537E-08,  &
   &      3.540E-08,  2.780E-08, -6.711E-09, -6.391E-10,  2.739E-08/
   DATA (Z(I),I= 1051, 1100)/                                        &
   &      3.530E-08,  4.290E-08,  4.289E-08,  5.107E-08,  4.316E-08,  &
   &      4.345E-08,  5.102E-08,  2.937E-08,  3.696E-08,  3.699E-08,  &
   &      3.695E-08,  5.099E-08,  4.486E-08,  4.480E-08,  4.479E-08,  &
   &      5.240E-08,  4.448E-08,  6.496E-08,  7.107E-08,  8.540E-08,  &
   &      5.706E-08,  2.900E-08,  4.300E-08,  7.128E-08,  9.295E-08,  &
   &      9.322E-08,  9.936E-08,  9.935E-08,  9.934E-08,  1.134E-07,  &
   &      9.146E-08,  1.195E-07,  1.055E-07,  9.177E-08,  9.207E-08,  &
   &      9.848E-08,  9.876E-08,  9.846E-08,  9.846E-08,  1.125E-07,  &
   &      1.122E-07,  1.122E-07,  1.183E-07,  1.324E-07,  1.321E-07,  &
   &      1.245E-07,  1.245E-07,  1.108E-07,  1.251E-07,  1.315E-07/
   DATA (Z(I),I= 1101, 1150)/                                        &
   &      1.380E-07,  1.303E-07,  1.227E-07,  1.148E-07,  1.148E-07,  &
   &      1.289E-07,  1.350E-07,  1.210E-07,  1.351E-07,  1.556E-07,  &
   &      1.559E-07,  1.479E-07,  1.681E-07,  1.465E-07,  1.591E-07,  &
   &      1.591E-07,  1.667E-07,  1.808E-07,  1.670E-07,  1.667E-07,  &
   &      1.592E-07,  1.376E-07,  1.498E-07,  1.562E-07,  1.986E-07,  &
   &      2.206E-07,  2.209E-07,  2.147E-07,  2.069E-07,  2.069E-07,  &
   &      2.130E-07,  1.990E-07,  2.051E-07,  2.268E-07,  2.552E-07,  &
   &      2.689E-07,  2.973E-07,  3.034E-07,  3.177E-07,  3.101E-07,  &
   &      3.242E-07,  3.102E-07,  3.038E-07,  3.053E-07,  2.909E-07,  &
   &      2.909E-07,  2.988E-07,  2.988E-07,  3.128E-07,  3.128E-07/
   DATA (Z(I),I= 1151, 1200)/                                        &
   &      3.193E-07,  3.257E-07,  3.257E-07,  3.461E-07,  3.321E-07,  &
   &      3.321E-07,  3.257E-07,  3.256E-07,  3.257E-07,  3.394E-07,  &
   &      3.391E-07,  3.248E-07,  3.248E-07,  3.251E-07,  3.172E-07,  &
   &      3.175E-07,  3.239E-07,  3.236E-07,  3.028E-07,  3.229E-07,  &
   &      3.306E-07,  3.309E-07,  3.528E-07,  3.531E-07,  3.394E-07,  &
   &      3.315E-07,  3.175E-07,  3.236E-07,  3.017E-07,  3.157E-07,  &
   &      3.157E-07,  3.295E-07,  3.154E-07,  3.154E-07,  3.093E-07,  &
   &      2.952E-07,  2.812E-07,  2.671E-07,  2.672E-07,  2.734E-07,  &
   &      2.731E-07,  2.792E-07,  2.997E-07,  2.933E-07,  2.792E-07,  &
   &      2.792E-07,  2.653E-07,  2.372E-07,  2.289E-07,  2.290E-07/
   DATA (Z(I),I= 1201, 1250)/                                        &
   &      2.276E-07,  2.132E-07,  2.211E-07,  2.431E-07,  2.369E-07,  &
   &      2.448E-07,  2.388E-07,  2.247E-07,  2.043E-07,  2.043E-07,  &
   &      1.903E-07,  1.967E-07,  2.251E-07,  2.235E-07,  2.300E-07,  &
   &      2.301E-07,  2.234E-07,  2.295E-07,  2.368E-07,  2.365E-07,  &
   &      2.368E-07,  2.089E-07,  2.093E-07,  1.894E-07,  1.757E-07,  &
   &      1.895E-07,  2.035E-07,  2.319E-07,  2.394E-07,  2.236E-07,  &
   &      2.376E-07,  1.937E-07,  1.658E-07,  1.374E-07,  1.456E-07,  &
   &      1.395E-07,  1.690E-07,  1.769E-07,  1.849E-07,  1.988E-07,  &
   &      1.988E-07,  1.769E-07,  1.708E-07,  1.976E-07,  2.042E-07,  &
   &      2.105E-07,  2.171E-07,  2.310E-07,  2.249E-07,  2.389E-07/
   DATA (Z(I),I= 1251, 1300)/                                        &
   &      2.246E-07,  2.109E-07,  1.822E-07,  1.682E-07,  1.536E-07,  &
   &      1.467E-07,  1.529E-07,  1.593E-07,  1.521E-07,  1.866E-07,  &
   &      1.930E-07,  1.729E-07,  1.947E-07,  1.745E-07,  1.623E-07,  &
   &      1.562E-07,  1.636E-07,  1.629E-07,  1.626E-07,  1.623E-07,  &
   &      1.688E-07,  1.750E-07,  1.676E-07,  1.878E-07,  1.946E-07,  &
   &      2.151E-07,  2.086E-07,  2.086E-07,  1.946E-07,  2.087E-07,  &
   &      2.007E-07,  1.784E-07,  1.630E-07,  1.613E-07,  1.551E-07,  &
   &      1.425E-07,  1.507E-07,  1.242E-07,  1.184E-07,  1.259E-07,  &
   &      1.605E-07,  1.540E-07,  1.962E-07,  2.102E-07,  1.804E-07,  &
   &      1.585E-07,  1.225E-07,  9.271E-08,  7.834E-08,  1.067E-07/
   DATA (Z(I),I= 1301, 1350)/                                        &
   &      1.144E-07,  1.205E-07,  1.489E-07,  1.568E-07,  1.162E-07,  &
   &      1.039E-07,  1.399E-07,  1.196E-07,  1.263E-07,  1.267E-07,  &
   &      9.704E-08,  8.956E-08,  6.793E-08,  9.629E-08,  1.241E-07,  &
   &      1.383E-07,  1.662E-07,  1.742E-07,  1.256E-07,  7.070E-08,  &
   &      6.420E-08,  3.083E-08,  5.131E-08,  5.781E-08,  7.816E-08,  &
   &      6.300E-08,  8.356E-08,  7.602E-08,  9.796E-08,  9.801E-08,  &
   &      1.402E-07,  1.758E-07,  1.478E-07,  1.261E-07,  6.339E-08,  &
   &      2.789E-08, -2.071E-08, -2.791E-08, -4.312E-08, -3.657E-08,  &
   &     -9.756E-09,  3.890E-08,  2.843E-09, -3.339E-09, -3.133E-08,  &
   &     -5.142E-08, -7.951E-08, -5.144E-08, -4.518E-08, -9.156E-09/
   DATA (Z(I),I= 1351, 1400)/                                        &
   &      1.096E-08,  1.085E-08,  3.110E-08,  5.755E-08,  7.769E-08,  &
   &      6.960E-08,  6.819E-08,  4.017E-08,  2.609E-08, -8.452E-09,  &
   &     -2.233E-08, -4.403E-08, -1.706E-08,  1.748E-08,  2.296E-09,  &
   &      1.513E-08, -6.365E-09, -3.449E-08, -6.250E-08, -5.452E-08,  &
   &     -8.259E-08, -1.186E-07, -1.342E-07, -1.281E-07, -1.579E-07,  &
   &     -1.736E-07, -1.812E-07, -1.666E-07, -1.662E-07, -1.655E-07,  &
   &     -1.787E-07, -1.864E-07, -2.023E-07, -2.179E-07, -2.260E-07,  &
   &     -2.417E-07, -2.499E-07, -2.504E-07, -2.510E-07, -2.516E-07,  &
   &     -2.657E-07, -2.738E-07, -2.889E-07, -3.182E-07, -3.259E-07,  &
   &     -3.346E-07, -3.486E-07, -3.754E-07, -3.819E-07, -3.947E-07/
   DATA (Z(I),I= 1401, 1450)/                                        &
   &     -3.998E-07, -4.139E-07, -4.214E-07, -4.162E-07, -4.378E-07,  &
   &     -4.672E-07, -4.888E-07, -4.964E-07, -4.825E-07, -4.964E-07,  &
   &     -4.964E-07, -5.040E-07, -5.524E-07, -4.964E-07, -4.963E-07,  &
   &     -4.963E-07, -4.963E-07, -4.963E-07, -4.964E-07, -4.963E-07,  &
   &     -4.964E-07, -4.963E-07, -4.963E-07, -4.963E-07, -4.963E-07,  &
   &     -5.724E-07, -4.964E-07, -4.963E-07, -6.365E-07, -4.963E-07,  &
   &     -5.577E-07, -6.978E-07, -5.575E-07, -5.576E-07, -6.950E-07,  &
   &     -5.577E-07, -5.576E-07, -6.951E-07, -5.577E-07, -5.577E-07,  &
   &     -5.577E-07, -5.576E-07, -5.576E-07, -6.950E-07, -5.577E-07,  &
   &     -6.978E-07, -5.576E-07, -5.577E-07, -5.576E-07, -5.577E-07/
   DATA (Z(I),I= 1451, 1500)/                                        &
   &     -5.576E-07, -5.577E-07, -5.576E-07, -7.010E-07, -8.381E-07,  &
   &     -8.383E-07, -6.977E-07, -8.411E-07, -8.411E-07, -5.604E-07,  &
   &     -8.411E-07, -7.005E-07, -8.410E-07, -7.006E-07, -6.978E-07,  &
   &     -8.379E-07, -7.006E-07, -8.381E-07, -8.378E-07, -8.349E-07,  &
   &     -8.380E-07, -8.380E-07, -8.351E-07, -8.379E-07, -8.379E-07,  &
   &     -8.380E-07, -8.379E-07, -8.350E-07, -8.350E-07, -8.350E-07,  &
   &     -8.351E-07, -9.751E-07, -9.751E-07, -8.350E-07, -8.351E-07,  &
   &     -8.380E-07, -8.379E-07, -9.751E-07, -8.377E-07, -8.379E-07,  &
   &     -8.379E-07, -9.140E-07, -9.139E-07, -8.527E-07, -9.286E-07,  &
   &     -8.526E-07, -9.286E-07, -7.881E-07, -8.526E-07, -7.882E-07/
   DATA (Z(I),I= 1501, 1550)/                                        &
   &     -6.477E-07, -7.122E-07, -7.121E-07, -7.883E-07, -6.475E-07,  &
   &     -8.524E-07, -8.522E-07, -6.478E-07, -7.881E-07, -7.883E-07,  &
   &     -7.124E-07, -6.475E-07, -6.476E-07, -7.238E-07, -6.478E-07,  &
   &     -5.834E-07, -8.637E-07, -5.832E-07, -8.025E-07, -5.832E-07,  &
   &     -5.219E-07, -5.218E-07, -7.412E-07, -7.413E-07, -6.010E-07,  &
   &     -6.011E-07, -7.413E-07, -7.410E-07, -5.221E-07, -5.219E-07,  &
   &     -4.576E-07, -6.007E-07, -4.604E-07, -4.603E-07, -3.814E-07,  &
   &     -3.169E-07, -3.199E-07, -3.813E-07, -3.813E-07, -4.601E-07,  &
   &     -3.201E-07, -3.203E-07, -3.201E-07, -3.200E-07, -3.202E-07,  &
   &     -3.991E-07, -2.587E-07, -2.560E-07, -2.589E-07, -2.587E-07/
   DATA (Z(I),I= 1551, 1600)/                                        &
   &     -2.591E-07, -3.202E-07, -3.950E-08, -3.960E-08,  3.970E-08,  &
   &      3.965E-08, -1.009E-07, -3.950E-08,  3.965E-08, -1.009E-07,  &
   &      1.038E-07, -3.955E-08, -3.950E-08, -3.955E-08, -3.960E-08,  &
   &     -3.950E-08, -3.955E-08,  2.176E-08,  9.783E-08, -4.263E-08,  &
   &      3.668E-08, -1.039E-07,  3.329E-08,  3.605E-08, -2.791E-08,  &
   &     -2.515E-08,  1.121E-07,  1.149E-07,  1.150E-07,  1.150E-07,  &
   &      1.148E-07,  1.150E-07,  1.149E-07,  1.150E-07, -2.838E-08,  &
   &     -2.823E-08, -2.833E-08,  1.118E-07, -3.157E-08,  4.727E-08,  &
   &     -9.266E-08, -9.276E-08,  5.061E-08,  2.666E-07,  1.263E-07,  &
   &      1.906E-07,  1.262E-07, -1.403E-08, -1.544E-07, -9.308E-08/
   DATA (Z(I),I= 1601, 1650)/                                        &
   &      5.024E-08, -9.318E-08, -1.543E-07, -1.545E-07, -1.435E-08,  &
   &     -1.440E-08, -7.873E-08, -7.863E-08, -7.878E-08, -7.883E-08,  &
   &     -7.878E-08, -7.878E-08, -2.221E-07, -7.873E-08, -2.833E-07,  &
   &     -7.878E-08, -1.432E-07, -2.220E-07, -2.865E-07, -3.478E-07,  &
   &     -3.478E-07, -3.478E-07, -3.478E-07, -4.091E-07, -3.448E-07,  &
   &     -3.448E-07, -3.448E-07, -2.046E-07, -2.046E-07, -4.090E-07,  &
   &     -2.044E-07, -3.446E-07, -4.092E-07, -2.045E-07, -2.046E-07,  &
   &     -2.045E-07, -2.046E-07, -2.046E-07, -2.044E-07, -2.047E-07,  &
   &     -2.045E-07, -2.046E-07, -2.046E-07, -4.852E-07, -3.446E-07,  &
   &     -2.806E-07, -2.806E-07, -2.805E-07, -5.611E-07, -3.418E-07/
   DATA (Z(I),I= 1651, 1700)/                                        &
   &     -3.416E-07, -4.820E-07, -4.031E-07, -3.418E-07, -3.418E-07,  &
   &     -4.819E-07, -4.059E-07, -3.447E-07, -3.447E-07, -3.447E-07,  &
   &     -3.447E-07, -4.848E-07, -3.447E-07, -3.446E-07, -3.448E-07,  &
   &     -3.447E-07, -3.446E-07, -3.447E-07, -3.447E-07, -4.202E-07,  &
   &     -4.203E-07, -4.205E-07, -4.204E-07, -2.801E-07, -2.802E-07,  &
   &     -2.801E-07, -4.203E-07, -4.203E-07, -2.830E-07, -2.830E-07,  &
   &     -4.230E-07, -4.233E-07, -2.831E-07, -4.230E-07, -2.829E-07,  &
   &     -3.619E-07, -3.590E-07, -3.618E-07, -4.992E-07, -4.991E-07,  &
   &     -4.992E-07, -4.378E-07, -4.378E-07, -2.973E-07, -7.816E-08,  &
   &     -2.974E-07, -2.186E-07, -2.186E-07, -2.186E-07, -7.800E-08/
   DATA (Z(I),I= 1701, 1750)/                                        &
   &     -7.800E-08, -2.181E-07, -7.816E-08, -2.182E-07, -7.810E-08,  &
   &     -1.568E-07, -7.758E-08, -7.774E-08,  6.250E-08, -7.810E-08,  &
   &     -2.183E-07, -1.568E-07, -1.569E-07, -2.179E-07, -7.810E-08,  &
   &     -7.800E-08, -7.795E-08, -7.800E-08, -2.973E-07, -2.973E-07,  &
   &     -7.805E-08, -2.185E-07, -3.586E-07, -3.588E-07, -2.973E-07,  &
   &     -1.541E-07, -3.586E-07, -2.152E-07, -7.805E-08, -2.154E-07,  &
   &     -1.541E-07, -7.518E-08, -2.154E-07, -7.523E-08, -1.393E-07,  &
   &     -7.518E-08, -1.396E-07,  7.826E-10, -1.396E-07, -2.798E-07,  &
   &      9.391E-10, -2.798E-07, -2.799E-07, -1.396E-07, -1.396E-07,  &
   &     -1.425E-07, -1.397E-07, -1.397E-07, -1.430E-07, -1.427E-07/
   DATA (Z(I),I= 1751, 1800)/                                        &
   &     -1.430E-07, -2.609E-09,  4.174E-10,  5.217E-10,  3.652E-10,  &
   &     -2.348E-09, -1.397E-07,  3.652E-10,  5.739E-10,  4.696E-10,  &
   &      4.696E-10,  4.696E-10,  2.609E-10,  4.696E-10, -6.407E-08,  &
   &      1.406E-07,  7.654E-08, -6.402E-08, -6.402E-08, -6.412E-08,  &
   &     -6.402E-08,  1.205E-08, -6.396E-08, -6.407E-08, -6.402E-08,  &
   &     -6.407E-08, -6.402E-08, -6.402E-08, -1.252E-07,  1.487E-08,  &
   &     -6.130E-08, -6.130E-08, -1.256E-07, -2.046E-07, -2.045E-07,  &
   &     -6.130E-08, -6.130E-08, -6.438E-08, -6.449E-08, -1.288E-07,  &
   &     -6.438E-08, -6.449E-08, -6.449E-08, -6.443E-08, -6.443E-08,  &
   &     -6.438E-08, -6.438E-08, -1.258E-07,  1.435E-08,  1.440E-08/
   DATA (Z(I),I= 1801, 1850)/                                        &
   &      1.440E-08, -2.645E-08, -4.043E-08, -8.129E-08, -8.791E-08,  &
   &     -7.268E-08, -5.113E-08, -5.770E-08, -7.821E-08, -9.219E-08,  &
   &     -1.203E-07, -1.483E-07, -1.639E-07, -1.921E-07, -1.559E-07,  &
   &     -1.422E-07, -1.278E-07, -7.962E-08, -1.560E-08,  4.664E-08,  &
   &      6.057E-08,  4.649E-08,  1.883E-08,  4.696E-09, -2.332E-08,  &
   &     -3.751E-08, -3.751E-08, -5.937E-08, -6.553E-08, -7.983E-08,  &
   &     -7.983E-08, -6.584E-08, -5.932E-08, -6.548E-08, -1.079E-07,  &
   &     -1.436E-07, -1.717E-07, -1.576E-07, -1.782E-07, -1.719E-07,  &
   &     -1.502E-07, -1.489E-07, -9.939E-08, -9.814E-08, -9.819E-08,  &
   &     -1.110E-07, -1.174E-07, -1.313E-07, -1.454E-07, -1.594E-07/
   DATA (Z(I),I= 1851, 1900)/                                        &
   &     -1.454E-07, -1.314E-07, -1.173E-07, -1.044E-07, -1.121E-07,  &
   &     -1.198E-07, -1.554E-07, -1.630E-07, -1.502E-07, -1.798E-07,  &
   &     -1.797E-07, -1.578E-07, -1.653E-07, -1.654E-07, -1.657E-07,  &
   &     -1.578E-07, -1.517E-07, -1.514E-07, -1.374E-07, -1.233E-07,  &
   &     -1.297E-07, -1.516E-07, -1.513E-07, -1.798E-07, -1.793E-07,  &
   &     -1.933E-07, -1.934E-07, -2.013E-07, -1.796E-07, -1.513E-07,  &
   &     -1.437E-07, -1.145E-07, -1.007E-07, -9.318E-08, -9.198E-08,  &
   &     -9.835E-08, -1.049E-07, -9.078E-08, -9.083E-08, -9.845E-08,  &
   &     -1.266E-07, -1.545E-07, -2.108E-07, -2.169E-07, -1.949E-07,  &
   &     -2.011E-07, -1.853E-07, -1.853E-07, -1.555E-07, -1.413E-07/
   DATA (Z(I),I= 1901, 1950)/                                        &
   &     -1.333E-07, -1.474E-07, -1.478E-07, -1.619E-07, -1.834E-07,  &
   &     -2.255E-07, -2.536E-07, -2.536E-07, -2.472E-07, -2.677E-07,  &
   &     -2.255E-07, -2.055E-07, -1.555E-07, -1.495E-07, -1.480E-07,  &
   &     -1.684E-07, -1.748E-07, -1.955E-07, -1.880E-07, -1.945E-07,  &
   &     -2.007E-07, -1.926E-07, -1.989E-07, -1.910E-07, -1.831E-07,  &
   &     -1.469E-07, -1.393E-07, -1.454E-07, -1.658E-07, -2.003E-07,  &
   &     -2.352E-07, -2.133E-07, -1.981E-07, -1.844E-07, -1.705E-07,  &
   &     -1.780E-07, -1.932E-07, -1.930E-07, -1.724E-07, -1.521E-07,  &
   &     -1.176E-07, -1.393E-07, -1.676E-07, -1.812E-07, -1.955E-07,  &
   &     -1.672E-07, -1.536E-07, -1.252E-07, -1.113E-07, -9.715E-08/
   DATA (Z(I),I= 1951, 2000)/                                        &
   &     -9.725E-08, -6.918E-08, -1.112E-07, -9.720E-08, -1.256E-07,  &
   &     -1.111E-07, -8.343E-08, -8.309E-08, -5.502E-08, -2.687E-08,  &
   &     -2.721E-08, -2.723E-08, -2.079E-08, -1.463E-08,  5.113E-09,  &
   &      1.158E-08, -5.241E-08, -1.165E-07, -1.602E-07, -1.805E-07,  &
   &     -1.456E-07, -1.176E-07, -1.169E-07, -1.166E-07, -1.583E-07,  &
   &     -1.580E-07, -1.658E-07, -1.382E-07, -1.103E-07, -8.259E-08,  &
   &     -1.109E-07, -1.252E-07, -1.611E-07, -1.754E-07, -1.626E-07,  &
   &     -1.483E-07, -1.561E-07, -1.574E-07, -1.790E-07, -2.005E-07,  &
   &     -2.363E-07, -2.081E-07, -2.082E-07, -1.801E-07, -2.020E-07,  &
   &     -2.160E-07, -2.082E-07, -2.020E-07, -1.741E-07, -1.325E-07/
   DATA (Z(I),I= 2001, 2050)/                                        &
   &     -1.047E-07, -9.861E-08, -1.127E-07, -9.884E-08, -1.131E-07,  &
   &     -1.270E-07, -1.406E-07, -1.546E-07, -1.478E-07, -1.194E-07,  &
   &     -1.333E-07, -1.050E-07, -1.190E-07, -1.114E-07, -1.319E-07,  &
   &     -1.321E-07, -1.246E-07, -1.450E-07, -1.590E-07, -1.656E-07,  &
   &     -1.935E-07, -2.140E-07, -2.205E-07, -2.345E-07, -2.550E-07,  &
   &     -2.475E-07, -2.615E-07, -2.335E-07, -2.411E-07, -2.414E-07,  &
   &     -2.273E-07, -1.928E-07, -1.788E-07, -1.709E-07, -1.709E-07,  &
   &     -1.850E-07, -1.909E-07, -1.908E-07, -2.049E-07, -1.967E-07,  &
   &     -2.032E-07, -2.452E-07, -2.589E-07, -3.211E-07, -3.352E-07,  &
   &     -3.133E-07, -2.710E-07, -2.275E-07, -2.278E-07, -2.203E-07/
   DATA (Z(I),I= 2051, 2100)/                                        &
   &     -2.267E-07, -2.392E-07, -2.600E-07, -2.665E-07, -2.729E-07,  &
   &     -2.793E-07, -2.934E-07, -3.153E-07, -3.291E-07, -3.229E-07,  &
   &     -3.153E-07, -2.871E-07, -2.730E-07, -2.950E-07, -3.166E-07,  &
   &     -3.510E-07, -3.791E-07, -3.446E-07, -3.446E-07, -3.447E-07,  &
   &     -3.307E-07, -3.231E-07, -3.512E-07, -3.436E-07, -3.439E-07,  &
   &     -3.503E-07, -3.500E-07, -3.565E-07, -3.708E-07, -3.845E-07,  &
   &     -3.845E-07, -4.065E-07, -3.986E-07, -4.061E-07, -4.065E-07,  &
   &     -4.205E-07, -3.925E-07, -3.861E-07, -3.580E-07, -3.720E-07,  &
   &     -3.720E-07, -4.001E-07, -4.141E-07, -4.140E-07, -4.138E-07,  &
   &     -4.338E-07, -4.263E-07, -4.123E-07, -3.900E-07, -3.620E-07/
   DATA (Z(I),I= 2101, 2150)/                                        &
   &     -3.480E-07, -3.265E-07, -3.408E-07, -3.268E-07, -3.470E-07,  &
   &     -3.674E-07, -4.095E-07, -4.515E-07, -4.723E-07, -4.863E-07,  &
   &     -5.143E-07, -5.205E-07, -4.925E-07, -4.571E-07, -4.151E-07,  &
   &     -4.011E-07, -3.932E-07, -3.713E-07, -3.652E-07, -3.574E-07,  &
   &     -3.495E-07, -3.632E-07, -3.635E-07, -3.977E-07, -4.602E-07,  &
   &     -5.084E-07, -5.569E-07, -5.286E-07, -4.930E-07, -4.431E-07,  &
   &     -4.216E-07, -4.061E-07, -4.262E-07, -4.187E-07, -4.248E-07,  &
   &     -4.234E-07, -4.018E-07, -4.219E-07, -4.223E-07, -4.640E-07,  &
   &     -4.924E-07, -5.283E-07, -5.424E-07, -5.223E-07, -4.802E-07,  &
   &     -4.522E-07, -4.241E-07, -4.178E-07, -4.181E-07, -3.976E-07/
   DATA (Z(I),I= 2151, 2200)/                                        &
   &     -4.256E-07, -4.472E-07, -4.674E-07, -4.750E-07, -4.825E-07,  &
   &     -4.702E-07, -4.562E-07, -4.498E-07, -4.716E-07, -4.792E-07,  &
   &     -4.867E-07, -4.884E-07, -4.820E-07, -4.618E-07, -4.355E-07,  &
   &     -4.290E-07, -4.664E-07, -5.023E-07, -5.177E-07, -5.191E-07,  &
   &     -4.852E-07, -4.446E-07, -4.101E-07, -3.975E-07, -4.272E-07,  &
   &     -4.426E-07, -4.505E-07, -4.235E-07, -3.815E-07, -3.330E-07,  &
   &     -2.910E-07, -3.124E-07, -3.279E-07, -3.699E-07, -3.978E-07,  &
   &     -3.914E-07, -3.569E-07, -3.166E-07, -2.821E-07, -2.896E-07,  &
   &     -3.115E-07, -3.269E-07, -3.485E-07, -3.280E-07, -3.140E-07,  &
   &     -2.999E-07, -2.863E-07, -3.079E-07, -3.139E-07, -3.139E-07/
   DATA (Z(I),I= 2201, 2250)/                                        &
   &     -3.218E-07, -3.075E-07, -3.496E-07, -3.417E-07, -3.277E-07,  &
   &     -3.415E-07, -3.476E-07, -3.695E-07, -3.757E-07, -3.617E-07,  &
   &     -3.617E-07, -3.603E-07, -3.463E-07, -3.463E-07, -3.448E-07,  &
   &     -3.448E-07, -3.448E-07, -3.509E-07, -3.714E-07, -3.775E-07,  &
   &     -3.918E-07, -3.840E-07, -3.764E-07, -3.545E-07, -3.265E-07,  &
   &     -3.531E-07, -3.812E-07, -4.156E-07, -4.361E-07, -4.362E-07,  &
   &     -4.143E-07, -4.208E-07, -3.928E-07, -4.068E-07, -4.349E-07,  &
   &     -4.553E-07, -4.895E-07, -4.755E-07, -4.616E-07, -4.256E-07,  &
   &     -4.119E-07, -4.181E-07, -4.325E-07, -4.386E-07, -4.670E-07,  &
   &     -4.734E-07, -4.734E-07, -4.659E-07, -4.723E-07, -5.004E-07/
   DATA (Z(I),I= 2251, 2300)/                                        &
   &     -5.141E-07, -5.342E-07, -5.763E-07, -5.623E-07, -5.404E-07,  &
   &     -5.124E-07, -4.905E-07, -4.905E-07, -4.905E-07, -5.109E-07,  &
   &     -5.250E-07, -5.455E-07, -5.315E-07, -5.379E-07, -5.382E-07,  &
   &     -5.303E-07, -5.239E-07, -4.958E-07, -5.020E-07, -5.098E-07,  &
   &     -4.879E-07, -4.879E-07, -4.800E-07, -4.876E-07, -4.879E-07,  &
   &     -5.016E-07, -5.159E-07, -5.439E-07, -5.296E-07, -5.158E-07,  &
   &     -5.159E-07, -5.094E-07, -4.875E-07, -4.735E-07, -4.735E-07,  &
   &     -4.735E-07, -4.811E-07, -4.671E-07, -4.954E-07, -4.893E-07,  &
   &     -4.673E-07, -4.534E-07, -4.394E-07, -4.037E-07, -4.037E-07,  &
   &     -4.241E-07, -4.662E-07, -4.942E-07, -5.226E-07, -5.147E-07/
   DATA (Z(I),I= 2301, 2350)/                                        &
   &     -5.007E-07, -5.147E-07, -5.007E-07, -5.007E-07, -4.928E-07,  &
   &     -5.008E-07, -4.791E-07, -5.072E-07, -5.198E-07, -5.338E-07,  &
   &     -5.478E-07, -5.324E-07, -5.327E-07, -5.108E-07, -4.892E-07,  &
   &     -5.093E-07, -5.514E-07, -5.715E-07, -6.060E-07, -6.481E-07,  &
   &     -5.982E-07, -5.766E-07, -5.483E-07, -5.203E-07, -5.344E-07,  &
   &     -5.405E-07, -5.265E-07, -5.469E-07, -5.405E-07, -5.466E-07,  &
   &     -5.326E-07, -5.324E-07, -5.607E-07, -5.542E-07, -5.747E-07,  &
   &     -5.826E-07, -5.966E-07, -5.966E-07, -5.765E-07, -5.625E-07,  &
   &     -5.487E-07, -5.423E-07, -5.438E-07, -5.300E-07, -5.236E-07,  &
   &     -5.311E-07, -5.527E-07, -5.466E-07, -5.541E-07, -5.480E-07/
   DATA (Z(I),I= 2351, 2400)/                                        &
   &     -5.415E-07, -5.289E-07, -5.009E-07, -4.883E-07, -5.099E-07,  &
   &     -5.037E-07, -4.973E-07, -4.911E-07, -4.286E-07, -3.660E-07,  &
   &     -3.035E-07, -2.489E-07, -2.844E-07, -3.340E-07, -3.760E-07,  &
   &     -4.335E-07, -4.410E-07, -4.269E-07, -4.488E-07, -4.348E-07,  &
   &     -4.423E-07, -4.286E-07, -4.221E-07, -4.081E-07, -4.019E-07,  &
   &     -4.019E-07, -3.599E-07, -3.459E-07, -3.179E-07, -3.395E-07,  &
   &     -3.473E-07, -3.754E-07, -4.034E-07, -4.031E-07, -4.033E-07,  &
   &     -3.969E-07, -3.969E-07, -3.890E-07, -3.891E-07, -4.031E-07,  &
   &     -4.031E-07, -4.031E-07, -4.031E-07, -4.031E-07, -4.095E-07,  &
   &     -4.235E-07, -4.019E-07, -4.160E-07, -3.958E-07, -3.882E-07/
   DATA (Z(I),I= 2401, 2450)/                                        &
   &     -3.742E-07, -3.804E-07, -3.944E-07, -4.011E-07, -4.073E-07,  &
   &     -4.073E-07, -3.994E-07, -3.997E-07, -3.919E-07, -4.059E-07,  &
   &     -3.983E-07, -3.983E-07, -4.045E-07, -3.902E-07, -3.963E-07,  &
   &     -3.963E-07, -4.025E-07, -4.165E-07, -4.230E-07, -4.510E-07,  &
   &     -4.793E-07, -5.276E-07, -5.492E-07, -5.212E-07, -5.074E-07,  &
   &     -4.935E-07, -4.576E-07, -4.576E-07, -4.640E-07, -4.921E-07,  &
   &     -4.985E-07, -4.982E-07, -4.985E-07, -4.985E-07, -4.944E-07,  &
   &     -4.952E-07, -5.176E-07, -5.243E-07, -5.594E-07, -5.868E-07,  &
   &     -6.078E-07, -6.009E-07, -5.939E-07, -5.869E-07, -5.502E-07,  &
   &     -5.515E-07, -5.487E-07, -5.599E-07, -5.661E-07, -5.695E-07/
   DATA (Z(I),I= 2451, 2500)/                                        &
   &     -5.667E-07, -5.574E-07, -5.560E-07, -5.532E-07, -5.221E-07,  &
   &     -5.314E-07, -5.061E-07, -4.870E-07, -4.811E-07, -4.715E-07,  &
   &     -4.773E-07, -4.683E-07, -4.601E-07, -4.733E-07, -4.546E-07,  &
   &     -4.414E-07, -4.385E-07, -4.191E-07, -4.363E-07, -4.309E-07,  &
   &     -4.395E-07, -4.499E-07, -4.353E-07, -4.183E-07, -4.019E-07,  &
   &     -4.052E-07, -3.913E-07, -3.900E-07, -3.902E-07, -3.931E-07,  &
   &     -3.953E-07, -3.960E-07, -3.944E-07, -3.943E-07, -3.948E-07,  &
   &     -3.911E-07, -3.817E-07, -3.704E-07, -3.583E-07, -3.475E-07,  &
   &     -3.403E-07, -3.394E-07, -3.356E-07, -3.333E-07, -3.309E-07,  &
   &     -3.257E-07, -3.206E-07, -3.148E-07, -3.096E-07, -3.108E-07/
   DATA (Z(I),I= 2501, 2550)/                                        &
   &     -3.247E-07, -3.364E-07, -3.517E-07, -3.635E-07, -3.555E-07,  &
   &     -3.400E-07, -3.238E-07, -3.082E-07, -3.003E-07, -3.170E-07,  &
   &     -3.359E-07, -3.518E-07, -3.707E-07, -3.698E-07, -3.543E-07,  &
   &     -3.388E-07, -3.232E-07, -3.092E-07, -3.154E-07, -3.275E-07,  &
   &     -3.359E-07, -3.486E-07, -3.578E-07, -3.578E-07, -3.607E-07,  &
   &     -3.643E-07, -3.672E-07, -3.668E-07, -3.645E-07, -3.628E-07,  &
   &     -3.626E-07, -3.652E-07, -3.859E-07, -4.088E-07, -4.338E-07,  &
   &     -4.595E-07, -4.744E-07, -4.748E-07, -4.758E-07, -4.761E-07,  &
   &     -4.772E-07, -4.909E-07, -5.102E-07, -5.308E-07, -5.501E-07,  &
   &     -5.653E-07, -5.538E-07, -5.375E-07, -5.190E-07, -5.028E-07/
   DATA (Z(I),I= 2551, 2600)/                                        &
   &     -4.997E-07, -5.220E-07, -5.441E-07, -5.685E-07, -5.907E-07,  &
   &     -5.939E-07, -5.815E-07, -5.711E-07, -5.615E-07, -5.499E-07,  &
   &     -5.515E-07, -5.545E-07, -5.589E-07, -5.625E-07, -5.621E-07,  &
   &     -5.540E-07, -5.438E-07, -5.328E-07, -5.212E-07, -5.144E-07,  &
   &     -5.105E-07, -5.065E-07, -5.034E-07, -4.994E-07, -4.926E-07,  &
   &     -4.837E-07, -4.761E-07, -4.672E-07, -4.611E-07, -4.521E-07,  &
   &     -4.445E-07, -4.370E-07, -4.294E-07, -4.233E-07, -4.186E-07,  &
   &     -4.159E-07, -4.126E-07, -4.093E-07, -4.018E-07, -3.901E-07,  &
   &     -3.776E-07, -3.659E-07, -3.542E-07, -3.515E-07, -3.558E-07,  &
   &     -3.600E-07, -3.643E-07, -3.651E-07, -3.546E-07, -3.379E-07/
   DATA (Z(I),I= 2601, 2650)/                                        &
   &     -3.204E-07, -3.058E-07, -2.933E-07, -2.947E-07, -2.975E-07,  &
   &     -3.010E-07, -3.052E-07, -3.066E-07, -2.933E-07, -2.793E-07,  &
   &     -2.667E-07, -2.563E-07, -2.499E-07, -2.561E-07, -2.603E-07,  &
   &     -2.686E-07, -2.762E-07, -2.712E-07, -2.599E-07, -2.459E-07,  &
   &     -2.360E-07, -2.240E-07, -2.218E-07, -2.260E-07, -2.294E-07,  &
   &     -2.328E-07, -2.370E-07, -2.342E-07, -2.306E-07, -2.292E-07,  &
   &     -2.250E-07, -2.243E-07, -2.236E-07, -2.216E-07, -2.216E-07,  &
   &     -2.210E-07, -2.225E-07, -2.246E-07, -2.290E-07, -2.318E-07,  &
   &     -2.354E-07, -2.376E-07, -2.313E-07, -2.245E-07, -2.198E-07,  &
   &     -2.129E-07, -2.102E-07, -2.160E-07, -2.210E-07, -2.253E-07/
   DATA (Z(I),I= 2651, 2700)/                                        &
   &     -2.297E-07, -2.306E-07, -2.237E-07, -2.176E-07, -2.121E-07,  &
   &     -2.066E-07, -2.062E-07, -2.077E-07, -2.092E-07, -2.143E-07,  &
   &     -2.144E-07, -2.140E-07, -2.119E-07, -2.085E-07, -2.058E-07,  &
   &     -2.010E-07, -1.997E-07, -1.970E-07, -1.936E-07, -1.915E-07,  &
   &     -1.896E-07, -1.874E-07, -1.838E-07, -1.816E-07, -1.794E-07,  &
   &     -1.778E-07, -1.750E-07, -1.706E-07, -1.698E-07, -1.662E-07,  &
   &     -1.648E-07, -1.632E-07, -1.602E-07, -1.580E-07, -1.586E-07,  &
   &     -1.557E-07, -1.549E-07, -1.527E-07, -1.525E-07, -1.489E-07,  &
   &     -1.487E-07, -1.479E-07, -1.486E-07, -1.484E-07, -1.491E-07,  &
   &     -1.489E-07, -1.488E-07, -1.494E-07, -1.493E-07, -1.485E-07/
   DATA (Z(I),I= 2701, 2750)/                                        &
   &     -1.484E-07, -1.518E-07, -1.526E-07, -1.547E-07, -1.569E-07,  &
   &     -1.575E-07, -1.611E-07, -1.632E-07, -1.640E-07, -1.654E-07,  &
   &     -1.682E-07, -1.696E-07, -1.740E-07, -1.783E-07, -1.811E-07,  &
   &     -1.855E-07, -1.885E-07, -1.927E-07, -1.956E-07, -1.992E-07,  &
   &     -2.028E-07, -2.051E-07, -2.067E-07, -2.069E-07, -2.099E-07,  &
   &     -2.094E-07, -2.109E-07, -2.126E-07, -2.141E-07, -2.136E-07,  &
   &     -2.131E-07, -2.140E-07, -2.135E-07, -2.124E-07, -2.105E-07,  &
   &     -2.093E-07, -2.081E-07, -2.055E-07, -2.058E-07, -2.039E-07,  &
   &     -2.028E-07, -2.023E-07, -1.996E-07, -1.984E-07, -1.951E-07,  &
   &     -1.959E-07, -1.940E-07, -1.928E-07, -1.901E-07, -1.889E-07/
   DATA (Z(I),I= 2751, 2800)/                                        &
   &     -1.870E-07, -1.843E-07, -1.823E-07, -1.782E-07, -1.755E-07,  &
   &     -1.728E-07, -1.702E-07, -1.647E-07, -1.634E-07, -1.593E-07,  &
   &     -1.566E-07, -1.518E-07, -1.491E-07, -1.464E-07, -1.436E-07,  &
   &     -1.423E-07, -1.381E-07, -1.368E-07, -1.334E-07, -1.299E-07,  &
   &     -1.287E-07, -1.258E-07, -1.232E-07, -1.203E-07, -1.163E-07,  &
   &     -1.157E-07, -1.129E-07, -1.088E-07, -1.081E-07, -1.047E-07,  &
   &     -1.019E-07, -1.006E-07, -9.780E-08, -9.514E-08, -9.455E-08,  &
   &     -9.330E-08, -9.128E-08, -8.923E-08, -8.864E-08, -8.879E-08,  &
   &     -8.536E-08, -8.615E-08, -8.414E-08, -8.288E-08, -8.226E-08,  &
   &     -8.305E-08, -8.305E-08, -8.244E-08, -8.183E-08, -8.198E-08/
   DATA (Z(I),I= 2801, 2850)/                                        &
   &     -8.478E-08, -8.417E-08, -8.356E-08, -8.435E-08, -8.370E-08,  &
   &     -8.370E-08, -8.432E-08, -8.575E-08, -8.651E-08, -8.651E-08,  &
   &     -8.510E-08, -8.648E-08, -8.648E-08, -8.788E-08, -8.928E-08,  &
   &     -8.785E-08, -8.645E-08, -8.770E-08, -8.627E-08, -8.627E-08,  &
   &     -8.627E-08, -8.685E-08, -8.685E-08, -8.609E-08, -8.606E-08,  &
   &     -8.462E-08, -8.257E-08, -8.394E-08, -8.330E-08, -8.186E-08,  &
   &     -8.265E-08, -8.186E-08, -8.122E-08, -7.978E-08, -7.913E-08,  &
   &     -7.910E-08, -7.641E-08, -7.517E-08, -7.453E-08, -7.326E-08,  &
   &     -7.060E-08, -6.998E-08, -6.732E-08, -6.886E-08, -6.541E-08,  &
   &     -6.417E-08, -6.353E-08, -6.230E-08, -6.308E-08, -6.246E-08/
   DATA (Z(I),I= 2851, 2900)/                                        &
   &     -6.325E-08, -6.126E-08, -6.205E-08, -6.000E-08, -6.078E-08,  &
   &     -5.815E-08, -5.893E-08, -5.896E-08, -5.896E-08, -5.960E-08,  &
   &     -5.960E-08, -5.823E-08, -6.027E-08, -6.013E-08, -6.048E-08,  &
   &     -6.090E-08, -6.196E-08, -6.042E-08, -6.145E-08, -6.167E-08,  &
   &     -5.982E-08, -6.005E-08, -6.021E-08, -6.029E-08, -5.844E-08,  &
   &     -6.004E-08, -5.883E-08, -5.981E-08, -5.916E-08, -5.857E-08,  &
   &     -5.779E-08, -5.930E-08, -5.790E-08, -5.801E-08, -5.863E-08,  &
   &     -5.798E-08, -5.719E-08, -5.717E-08, -5.717E-08, -5.714E-08,  &
   &     -5.653E-08, -5.933E-08, -5.872E-08, -5.808E-08, -5.948E-08,  &
   &     -5.887E-08, -5.965E-08, -5.963E-08, -5.901E-08, -6.121E-08/
   DATA (Z(I),I= 2901, 2950)/                                        &
   &     -6.056E-08, -6.135E-08, -5.991E-08, -6.146E-08, -6.146E-08,  &
   &     -6.222E-08, -6.300E-08, -6.376E-08, -6.315E-08, -6.455E-08,  &
   &     -6.530E-08, -6.465E-08, -6.544E-08, -6.339E-08, -6.414E-08,  &
   &     -6.493E-08, -6.367E-08, -6.506E-08, -6.377E-08, -6.377E-08,  &
   &     -6.315E-08, -6.390E-08, -6.283E-08, -6.306E-08, -6.258E-08,  &
   &     -6.160E-08, -6.109E-08, -6.061E-08, -6.010E-08, -5.830E-08,  &
   &     -5.820E-08, -5.667E-08, -5.671E-08, -5.521E-08, -5.468E-08,  &
   &     -5.347E-08, -5.241E-08, -5.126E-08, -5.032E-08, -4.923E-08,  &
   &     -4.822E-08, -4.715E-08, -4.619E-08, -4.512E-08, -4.411E-08,  &
   &     -4.294E-08, -4.198E-08, -4.082E-08, -3.972E-08, -3.848E-08/
   DATA (Z(I),I= 2951, 3000)/                                        &
   &     -3.765E-08, -3.656E-08, -3.531E-08, -3.421E-08, -3.319E-08,  &
   &     -3.209E-08, -3.131E-08, -3.082E-08, -3.010E-08, -2.947E-08,  &
   &     -2.890E-08, -2.833E-08, -2.776E-08, -2.718E-08, -2.669E-08,  &
   &     -2.626E-08, -2.569E-08, -2.517E-08, -2.505E-08, -2.494E-08,  &
   &     -2.504E-08, -2.485E-08, -2.481E-08, -2.454E-08, -2.444E-08,  &
   &     -2.453E-08, -2.435E-08, -2.416E-08, -2.425E-08, -2.419E-08,  &
   &     -2.441E-08, -2.454E-08, -2.448E-08, -2.470E-08, -2.463E-08,  &
   &     -2.491E-08, -2.499E-08, -2.485E-08, -2.506E-08, -2.514E-08,  &
   &     -2.521E-08, -2.494E-08, -2.501E-08, -2.488E-08, -2.488E-08,  &
   &     -2.460E-08, -2.475E-08, -2.468E-08, -2.470E-08, -2.462E-08/
   DATA (Z(I),I= 3001, 3050)/                                        &
   &     -2.442E-08, -2.429E-08, -2.417E-08, -2.398E-08, -2.373E-08,  &
   &     -2.353E-08, -2.314E-08, -2.295E-08, -2.276E-08, -2.264E-08,  &
   &     -2.245E-08, -2.212E-08, -2.200E-08, -2.163E-08, -2.146E-08,  &
   &     -2.114E-08, -2.091E-08, -2.081E-08, -2.064E-08, -2.047E-08,  &
   &     -2.009E-08, -1.999E-08, -1.982E-08, -1.959E-08, -1.907E-08,  &
   &     -1.885E-08, -1.841E-08, -1.804E-08, -1.760E-08, -1.724E-08,  &
   &     -1.674E-08, -1.630E-08, -1.579E-08, -1.564E-08, -1.513E-08,  &
   &     -1.483E-08, -1.431E-08, -1.401E-08, -1.350E-08, -1.320E-08,  &
   &     -1.289E-08, -1.238E-08, -1.228E-08, -1.190E-08, -1.146E-08,  &
   &     -1.116E-08, -1.064E-08, -1.033E-08, -1.040E-08, -1.028E-08/
   DATA (Z(I),I= 3051, 3100)/                                        &
   &     -1.042E-08, -1.036E-08, -1.023E-08, -1.037E-08, -1.045E-08,  &
   &     -1.025E-08, -1.033E-08, -1.026E-08, -1.040E-08, -1.054E-08,  &
   &     -1.080E-08, -1.111E-08, -1.150E-08, -1.190E-08, -1.207E-08,  &
   &     -1.260E-08, -1.272E-08, -1.317E-08, -1.342E-08, -1.368E-08,  &
   &     -1.399E-08, -1.426E-08, -1.466E-08, -1.493E-08, -1.533E-08,  &
   &     -1.568E-08, -1.586E-08, -1.621E-08, -1.654E-08, -1.688E-08,  &
   &     -1.714E-08, -1.741E-08, -1.775E-08, -1.805E-08, -1.794E-08,  &
   &     -1.819E-08, -1.814E-08, -1.810E-08, -1.848E-08, -1.844E-08,  &
   &     -1.868E-08, -1.864E-08, -1.874E-08, -1.885E-08, -1.887E-08,  &
   &     -1.864E-08, -1.835E-08, -1.819E-08, -1.789E-08, -1.760E-08/
   DATA (Z(I),I= 3101, 3150)/                                        &
   &     -1.722E-08, -1.706E-08, -1.663E-08, -1.639E-08, -1.610E-08,  &
   &     -1.580E-08, -1.536E-08, -1.480E-08, -1.459E-08, -1.403E-08,  &
   &     -1.346E-08, -1.304E-08, -1.248E-08, -1.199E-08, -1.157E-08,  &
   &     -1.115E-08, -1.045E-08, -1.002E-08, -9.538E-09, -8.929E-09,  &
   &     -8.535E-09, -8.001E-09, -7.532E-09, -6.920E-09, -6.422E-09,  &
   &     -6.094E-09, -5.472E-09, -4.954E-09, -4.367E-09, -3.874E-09,  &
   &     -3.500E-09, -3.300E-09, -3.000E-09, -2.700E-09, -2.500E-09,  &
   &     -2.300E-09, -2.100E-09, -1.900E-09, -1.700E-09, -1.500E-09,  &
   &     -1.300E-09, -1.125E-09, -0.950E-09, -0.775E-09, -0.600E-09,  &
   &     -0.425E-09, -0.250E-09, -0.100E-09, -0.005E-09,  0.00000  /
!
end block data O3CH
!
!     --------------------------------------------------------------
!
SUBROUTINE O3HHT0 (V1C,V2C,DVC,NPTC,C,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   COMMON /O3HH0/ V1S,V2S,DVS,NPTS,S(2687)
   DIMENSION C(*)
!
   DVC = DVS
   v1ss = v1s
   v2ss = v2s
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      VJ = V1C+DVC* REAL(J-1)
      C(J) = S(I)/VJ
!
!     RADIATION FLD REMOVED FROM DIFFUSE OZONE
!
10 END DO
!
   RETURN
!
end subroutine O3HHT0
!
!     --------------------------------------------------------------
!
BLOCK DATA BO3HH0
!
   IMPLICIT REAL*8           (V)
!
!     O3HH0 CONTAINS O3 HARTLEY HUGGINS CROSS SECTIONS FOR 273K
!               UNITS OF (CM**2/MOL)*1.E-20
!
!     NOW INCLUDES MOLINA & MOLINA AT 273K WITH THE TEMPERATURE
!     DEPENDENCE DETERMINED FROM THE 195K HARVARD MEASUREMENTS,
!     EMPLOYING THE BASS ALGORITHM
!
!              (CO(1+C1*(T-273.15)+C2*(T-273.15)**2);
!
!     THIS IS ONLY FOR THE WAVELENGTH RANGE FROM .34 TO .35 MICRONS;
!     OTHERWISE, THE BASS DATA ALONE HAVE BEEN EMPLOYED BETWEEN
!     .34 AND .245 MICRONS.
!
!     NEW T-DEPENDENT X-SECTIONS BETWEEN .345 AND .36 MICRONS
!     HAVE NOW BEEN ADDED, BASED ON WORK BY CACCIANI, DISARRA
!     AND FIOCCO, UNIVERSITY OF ROME, 1987.  QUADRATIC TEMP
!     HAS BEEN DERIVED, AS ABOVE.
!
!     MOLINA & MOLINA HAVE AGAIN BEEN USED BETWEEN .245 AND .185
!     MICRONS (NO TEMPERATURE DEPENDENCE)
!
!     AGREEMENT AMONGST THE FOUR DATA SETS IS REASONABLE (<10%)
!     AND OFTEN EXCELLENT (0-3%)
!
!
   COMMON /O3HH0/ V1C,V2C,DVC,NC,                                    &
   &               O30001(80),O30081(80),O30161(80),O30241(80),       &
   &               O30321(80),O30401( 7),                             &
   &               C00001(80),C00081(80),C00161(80),C00241(80),       &
   &               C00321(80),C00401(80),C00481(80),C00561(80),       &
   &               C00641(80),C00721(80),C00801(80),C00881(80),       &
   &               C00961(80),C01041(80),C01121(80),C01201(80),       &
   &               C01281(80),C01361(80),C01441(80),C01521(80),       &
   &               C01601(80),C01681(80),C01761(80),C01841(80),       &
   &               C01921(80),C02001(80),C02081(80),C02161(80),       &
   &               C02241(40)
!
!     DATA V1C  /27370./,V2C  /29400./,DVC  /5./,NC  /407/ INN & TANAKA
!         DATA FROM INN & TANAKA, HANDBOOK OF GEOPHYSICS, 1957, P 16-24
!                LINEARLY INTERPOLATED BY SAC, JUNE 1985
!                CONVERSION: (I&T)/(LOSCHMIDT 1 1987*1.2)
!
!     DATA V1C /29405./, V2C /40800./ ,DVC /5./, NC /2280/  BASS
!         DATA FROM BASS, JUNE 1985
!
   DATA V1C /27370./, V2C /40800./ ,DVC /5./, NC /2687/
!
!
!    X 2.08858E-03, 1.98947E-03, 1.89037E-03, 1.79126E-03, 1.69215E-03,
!     THIS LINE OF DATA HAS BEEN REPLACED BY MONOTONICALLY INCREASING
!     VALUES
!
   DATA O30001/                                                      &
   & 1.00000E-03, 1.15000E-03, 1.25000E-03, 1.40000E-03, 1.50000E-03, &
   & 1.59304E-03, 1.62396E-03, 1.76216E-03, 1.90036E-03, 2.03856E-03, &
   & 2.16538E-03, 2.02324E-03, 1.88110E-03, 1.73896E-03, 1.59682E-03, &
   & 1.45468E-03, 1.31253E-03, 1.17039E-03, 1.02825E-03, 8.86108E-04, &
   & 7.43963E-04, 6.01821E-04, 4.59679E-04, 5.14820E-04, 5.73044E-04, &
   & 6.31269E-04, 6.89493E-04, 7.47718E-04, 8.05942E-04, 8.64167E-04, &
   & 9.22392E-04, 9.80617E-04, 1.03884E-03, 1.09707E-03, 1.15528E-03, &
   & 1.21351E-03, 1.27173E-03, 1.32996E-03, 1.38818E-03, 1.44641E-03, &
   & 1.50463E-03, 1.56286E-03, 1.62108E-03, 1.67931E-03, 1.73753E-03, &
   & 1.79575E-03, 1.85398E-03, 1.91220E-03, 1.97043E-03, 2.02865E-03, &
   & 2.08688E-03, 2.14510E-03, 2.20333E-03, 2.26155E-03, 2.31978E-03, &
   & 2.37800E-03, 2.43623E-03, 2.49444E-03, 2.55267E-03, 2.61089E-03, &
   & 2.66912E-03, 2.72734E-03, 2.78557E-03, 2.84379E-03, 2.90202E-03, &
   & 2.96024E-03, 3.01847E-03, 3.07669E-03, 3.13491E-03, 3.19313E-03, &
   & 3.25136E-03, 3.30958E-03, 3.36781E-03, 3.31660E-03, 3.21583E-03, &
   & 3.11505E-03, 3.22165E-03, 3.46058E-03, 3.69953E-03, 3.93846E-03/
   DATA O30081/                                                      &
   & 4.17739E-03, 4.41633E-03, 4.42256E-03, 4.13791E-03, 4.17894E-03, &
   & 4.25583E-03, 4.33273E-03, 4.40963E-03, 4.49259E-03, 4.44532E-03, &
   & 4.17540E-03, 3.84814E-03, 3.41823E-03, 3.11003E-03, 2.86548E-03, &
   & 2.73912E-03, 2.70800E-03, 2.70882E-03, 2.70866E-03, 2.70816E-03, &
   & 2.71228E-03, 2.78044E-03, 2.86135E-03, 3.00163E-03, 3.15222E-03, &
   & 3.33394E-03, 3.48231E-03, 3.64966E-03, 3.83242E-03, 3.97733E-03, &
   & 4.10299E-03, 4.26332E-03, 4.41165E-03, 4.54040E-03, 4.65544E-03, &
   & 4.91897E-03, 5.23429E-03, 5.45390E-03, 5.74420E-03, 5.96314E-03, &
   & 6.07198E-03, 6.07338E-03, 5.99162E-03, 5.95079E-03, 6.04655E-03, &
   & 6.18239E-03, 6.56998E-03, 6.93885E-03, 7.38561E-03, 7.73029E-03, &
   & 7.90493E-03, 7.72072E-03, 7.40226E-03, 6.53860E-03, 5.30328E-03, &
   & 4.23000E-03, 3.45735E-03, 3.21167E-03, 3.16694E-03, 3.30966E-03, &
   & 3.47431E-03, 3.68089E-03, 3.92006E-03, 4.05246E-03, 4.16408E-03, &
   & 4.08710E-03, 3.98224E-03, 4.07316E-03, 4.19498E-03, 4.44990E-03, &
   & 4.77881E-03, 5.08270E-03, 5.37384E-03, 5.70240E-03, 5.91906E-03, &
   & 5.96745E-03, 5.92363E-03, 5.80363E-03, 5.60812E-03, 5.37450E-03/
   DATA O30161/                                                      &
   & 5.16202E-03, 4.98389E-03, 4.95294E-03, 5.04930E-03, 5.17576E-03, &
   & 5.26042E-03, 5.22957E-03, 5.32404E-03, 5.39630E-03, 5.53353E-03, &
   & 5.68057E-03, 5.78679E-03, 5.83795E-03, 5.93810E-03, 6.09330E-03, &
   & 6.40001E-03, 6.69056E-03, 7.04863E-03, 7.41339E-03, 7.87421E-03, &
   & 8.35570E-03, 8.97672E-03, 9.58486E-03, 1.01972E-02, 1.08463E-02, &
   & 1.14105E-02, 1.18935E-02, 1.22404E-02, 1.25053E-02, 1.28759E-02, &
   & 1.32169E-02, 1.37796E-02, 1.46488E-02, 1.57324E-02, 1.68897E-02, &
   & 1.78560E-02, 1.87101E-02, 1.92197E-02, 1.94106E-02, 1.90711E-02, &
   & 1.86585E-02, 1.82149E-02, 1.82219E-02, 1.85639E-02, 1.91924E-02, &
   & 2.01342E-02, 2.12312E-02, 2.26362E-02, 2.39610E-02, 2.55156E-02, &
   & 2.71338E-02, 2.87904E-02, 3.04268E-02, 3.17055E-02, 3.28248E-02, &
   & 3.36026E-02, 3.36867E-02, 3.26393E-02, 2.99356E-02, 2.56607E-02, &
   & 2.11545E-02, 1.79508E-02, 1.59757E-02, 1.49569E-02, 1.46214E-02, &
   & 1.46214E-02, 1.48217E-02, 1.51379E-02, 1.53816E-02, 1.58087E-02, &
   & 1.62186E-02, 1.66627E-02, 1.70961E-02, 1.76101E-02, 1.81759E-02, &
   & 1.86154E-02, 1.88889E-02, 1.89577E-02, 1.89316E-02, 1.88826E-02/
   DATA O30241/                                                      &
   & 1.90915E-02, 1.95550E-02, 2.02707E-02, 2.11620E-02, 2.21844E-02, &
   & 2.30920E-02, 2.37270E-02, 2.37422E-02, 2.33578E-02, 2.20358E-02, &
   & 1.96239E-02, 1.73329E-02, 1.57013E-02, 1.50566E-02, 1.49248E-02, &
   & 1.52044E-02, 1.57658E-02, 1.63436E-02, 1.68986E-02, 1.74180E-02, &
   & 1.78192E-02, 1.80677E-02, 1.79927E-02, 1.77900E-02, 1.75599E-02, &
   & 1.74982E-02, 1.76674E-02, 1.81633E-02, 1.87826E-02, 1.96898E-02, &
   & 2.06898E-02, 2.17167E-02, 2.28231E-02, 2.40702E-02, 2.55084E-02, &
   & 2.69701E-02, 2.86915E-02, 3.05796E-02, 3.22328E-02, 3.42637E-02, &
   & 3.61708E-02, 3.79118E-02, 3.94418E-02, 4.07333E-02, 4.17158E-02, &
   & 4.17081E-02, 4.01127E-02, 3.65411E-02, 3.25123E-02, 2.98737E-02, &
   & 2.83616E-02, 2.79907E-02, 2.80571E-02, 2.84778E-02, 2.91698E-02, &
   & 2.99500E-02, 3.07468E-02, 3.13903E-02, 3.19811E-02, 3.24616E-02, &
   & 3.26503E-02, 3.26829E-02, 3.27688E-02, 3.36446E-02, 3.55133E-02, &
   & 3.88447E-02, 4.28854E-02, 4.55381E-02, 4.77161E-02, 4.93567E-02, &
   & 4.95127E-02, 5.00492E-02, 5.06233E-02, 5.12739E-02, 5.20327E-02, &
   & 5.29001E-02, 5.38677E-02, 5.49272E-02, 5.60703E-02, 5.72886E-02/
   DATA O30321/                                                      &
   & 5.85739E-02, 5.99178E-02, 6.13170E-02, 6.28474E-02, 6.46499E-02, &
   & 6.68672E-02, 6.96421E-02, 7.31174E-02, 7.74361E-02, 8.27413E-02, &
   & 8.91756E-02, 9.67018E-02, 1.04844E-01, 1.13063E-01, 1.20818E-01, &
   & 1.27567E-01, 1.32771E-01, 1.35888E-01, 1.36377E-01, 1.33780E-01, &
   & 1.28385E-01, 1.20887E-01, 1.11978E-01, 1.02354E-01, 9.27108E-02, &
   & 8.37418E-02, 7.61423E-02, 7.06032E-02, 6.74255E-02, 6.62092E-02, &
   & 6.64813E-02, 6.77689E-02, 6.95995E-02, 7.15004E-02, 7.29991E-02, &
   & 7.36229E-02, 7.29641E-02, 7.11015E-02, 6.83345E-02, 6.49638E-02, &
   & 6.12897E-02, 5.76125E-02, 5.42326E-02, 5.14504E-02, 4.95645E-02, &
   & 4.87078E-02, 4.87234E-02, 4.94254E-02, 5.06280E-02, 5.21454E-02, &
   & 5.37919E-02, 5.53818E-02, 5.67293E-02, 5.76709E-02, 5.82319E-02, &
   & 5.85334E-02, 5.86968E-02, 5.88439E-02, 5.90963E-02, 5.95756E-02, &
   & 6.04035E-02, 6.17016E-02, 6.35548E-02, 6.59664E-02, 6.89282E-02, &
   & 7.24326E-02, 7.64718E-02, 8.10380E-02, 8.61236E-02, 9.17211E-02, &
   & 9.78192E-02, 1.04353E-01, 1.11218E-01, 1.18308E-01, 1.25519E-01, &
   & 1.32745E-01, 1.39881E-01, 1.46821E-01, 1.53461E-01, 1.59687E-01/
!
!    X 1.64187E-01, 1.69368E-01, 1.74549E-01, 1.79731E-01, 1.84912E-01,
!      1.90094E-01, 1.95275E-01/
!   THE VALUE AT 29400. HAS BEEN CHANGED TO PROVIDE A SMOOTH TRANSITION
!    X 1.90094E-01, 1.85275E-01/
!
   DATA O30401/                                                      &
   & 1.65365E-01, 1.70353E-01, 1.74507E-01, 1.77686E-01, 1.79748E-01, &
   & 1.80549E-01, 1.79948E-01/
!
!
!    FOLLOWING DATA ARE FROM BASS JUNE 1985
!
   DATA C00001 /                                                     &
   & 1.81094E-01, 1.57760E-01, 1.37336E-01, 1.19475E-01, 1.17191E-01, &
   & 1.14331E-01, 1.15984E-01, 1.10412E-01, 1.12660E-01, 1.16014E-01, &
   & 1.15060E-01, 1.12041E-01, 1.11611E-01, 1.00378E-01, 9.54850E-02, &
   & 9.87528E-02, 9.46153E-02, 9.53093E-02, 9.72653E-02, 9.66468E-02, &
   & 9.39750E-02, 1.03552E-01, 1.01361E-01, 1.04315E-01, 1.12842E-01, &
   & 1.02800E-01, 1.09576E-01, 1.05577E-01, 1.17334E-01, 1.25763E-01, &
   & 1.27597E-01, 1.34267E-01, 1.44799E-01, 1.57366E-01, 1.67369E-01, &
   & 1.81778E-01, 1.89207E-01, 2.01376E-01, 2.10310E-01, 2.21721E-01, &
   & 2.43162E-01, 2.55542E-01, 2.75312E-01, 2.88576E-01, 3.02505E-01, &
   & 3.15141E-01, 3.28908E-01, 3.49000E-01, 3.56620E-01, 3.59852E-01, &
   & 3.57517E-01, 3.12924E-01, 2.63610E-01, 2.50854E-01, 2.25642E-01, &
   & 2.15954E-01, 2.12099E-01, 2.13039E-01, 2.12286E-01, 2.17214E-01, &
   & 2.28784E-01, 2.28276E-01, 2.34677E-01, 2.30730E-01, 2.16107E-01, &
   & 1.99471E-01, 1.85629E-01, 1.72730E-01, 1.56229E-01, 1.38156E-01, &
   & 1.37641E-01, 1.33169E-01, 1.32759E-01, 1.30102E-01, 1.35396E-01, &
   & 1.37976E-01, 1.41571E-01, 1.46448E-01, 1.44508E-01, 1.47612E-01/
   DATA C00081 /                                                     &
   & 1.47424E-01, 1.48173E-01, 1.52936E-01, 1.58908E-01, 1.58808E-01, &
   & 1.59860E-01, 1.73936E-01, 1.84109E-01, 1.95143E-01, 2.08267E-01, &
   & 2.19256E-01, 2.31653E-01, 2.46400E-01, 2.60437E-01, 2.70792E-01, &
   & 2.79749E-01, 2.91068E-01, 2.98080E-01, 3.10421E-01, 3.24540E-01, &
   & 3.39003E-01, 3.58322E-01, 3.81520E-01, 4.02798E-01, 4.35972E-01, &
   & 4.56220E-01, 4.79037E-01, 5.02597E-01, 5.24648E-01, 5.33964E-01, &
   & 5.39211E-01, 5.43613E-01, 5.28793E-01, 4.94103E-01, 4.34481E-01, &
   & 3.76792E-01, 3.37161E-01, 3.15750E-01, 3.11042E-01, 3.08745E-01, &
   & 3.09195E-01, 3.05859E-01, 3.01443E-01, 2.88111E-01, 2.81303E-01, &
   & 2.75329E-01, 2.60812E-01, 2.59337E-01, 2.45576E-01, 2.40470E-01, &
   & 2.39705E-01, 2.45389E-01, 2.49801E-01, 2.53235E-01, 2.54387E-01, &
   & 2.64311E-01, 2.74146E-01, 2.89737E-01, 2.96673E-01, 3.07337E-01, &
   & 3.24380E-01, 3.42266E-01, 3.59522E-01, 3.78005E-01, 3.97178E-01, &
   & 4.23351E-01, 4.45925E-01, 4.63029E-01, 4.94843E-01, 5.19418E-01, &
   & 5.49928E-01, 5.69115E-01, 6.02396E-01, 6.43471E-01, 6.76401E-01, &
   & 7.14024E-01, 7.42425E-01, 7.60916E-01, 7.83319E-01, 7.98299E-01/
   DATA C00161 /                                                     &
   & 7.76672E-01, 7.22769E-01, 6.45967E-01, 5.80850E-01, 5.76514E-01, &
   & 5.79380E-01, 5.90359E-01, 6.21721E-01, 6.37540E-01, 6.52572E-01, &
   & 6.63442E-01, 6.69026E-01, 6.69038E-01, 6.53319E-01, 6.21950E-01, &
   & 5.47619E-01, 4.58994E-01, 4.14888E-01, 3.97736E-01, 3.88775E-01, &
   & 3.87424E-01, 3.93567E-01, 4.03442E-01, 4.05217E-01, 4.12848E-01, &
   & 4.12246E-01, 4.16620E-01, 4.13195E-01, 4.08467E-01, 4.13104E-01, &
   & 4.24498E-01, 4.32002E-01, 4.46361E-01, 4.61131E-01, 4.77228E-01, &
   & 4.96519E-01, 5.16764E-01, 5.38966E-01, 5.54187E-01, 5.73748E-01, &
   & 6.07260E-01, 6.34358E-01, 6.60286E-01, 6.95533E-01, 7.37090E-01, &
   & 7.83894E-01, 8.19557E-01, 8.49244E-01, 8.91832E-01, 9.44885E-01, &
   & 9.86271E-01, 1.02262E+00, 1.07242E+00, 1.12162E+00, 1.18287E+00, &
   & 1.22402E+00, 1.24978E+00, 1.24392E+00, 1.19668E+00, 1.11562E+00, &
   & 1.03983E+00, 9.31884E-01, 8.35307E-01, 7.92620E-01, 7.81980E-01, &
   & 7.89623E-01, 8.05987E-01, 8.27344E-01, 8.57514E-01, 8.66302E-01, &
   & 8.72092E-01, 8.66840E-01, 8.40536E-01, 7.87360E-01, 7.35743E-01, &
   & 6.92039E-01, 6.64032E-01, 6.48360E-01, 6.46288E-01, 6.49505E-01/
   DATA C00241 /                                                     &
   & 6.69937E-01, 6.81006E-01, 7.00969E-01, 7.19834E-01, 7.26964E-01, &
   & 7.50591E-01, 7.73600E-01, 8.00673E-01, 8.20347E-01, 8.37855E-01, &
   & 8.66780E-01, 9.04297E-01, 9.46300E-01, 9.69134E-01, 9.97928E-01, &
   & 1.06388E+00, 1.11032E+00, 1.15221E+00, 1.21324E+00, 1.24462E+00, &
   & 1.31978E+00, 1.35617E+00, 1.38792E+00, 1.39196E+00, 1.35161E+00, &
   & 1.29381E+00, 1.30295E+00, 1.32965E+00, 1.37024E+00, 1.44064E+00, &
   & 1.50484E+00, 1.57200E+00, 1.62097E+00, 1.67874E+00, 1.72676E+00, &
   & 1.73383E+00, 1.66091E+00, 1.54936E+00, 1.35454E+00, 1.20070E+00, &
   & 1.14609E+00, 1.13642E+00, 1.13784E+00, 1.14609E+00, 1.14531E+00, &
   & 1.16024E+00, 1.16891E+00, 1.16111E+00, 1.14192E+00, 1.09903E+00, &
   & 1.05745E+00, 1.02341E+00, 1.00121E+00, 1.00036E+00, 1.00576E+00, &
   & 1.02405E+00, 1.04379E+00, 1.07623E+00, 1.11347E+00, 1.17305E+00, &
   & 1.20016E+00, 1.22697E+00, 1.27479E+00, 1.32572E+00, 1.38690E+00, &
   & 1.43768E+00, 1.48379E+00, 1.55317E+00, 1.64020E+00, 1.71268E+00, &
   & 1.77183E+00, 1.85824E+00, 1.95131E+00, 2.04609E+00, 2.13151E+00, &
   & 2.17777E+00, 2.22832E+00, 2.26886E+00, 2.19775E+00, 2.05087E+00/
   DATA C00321 /                                                     &
   & 1.96103E+00, 1.95554E+00, 1.98037E+00, 2.05440E+00, 2.11629E+00, &
   & 2.17893E+00, 2.24384E+00, 2.30464E+00, 2.32525E+00, 2.29945E+00, &
   & 2.21712E+00, 2.03430E+00, 1.82139E+00, 1.70354E+00, 1.64631E+00, &
   & 1.62164E+00, 1.61356E+00, 1.63900E+00, 1.66313E+00, 1.67409E+00, &
   & 1.69143E+00, 1.70181E+00, 1.69165E+00, 1.67699E+00, 1.67879E+00, &
   & 1.67312E+00, 1.68133E+00, 1.70002E+00, 1.72500E+00, 1.76308E+00, &
   & 1.80634E+00, 1.87548E+00, 1.94924E+00, 1.99812E+00, 2.05333E+00, &
   & 2.14035E+00, 2.21847E+00, 2.27412E+00, 2.29752E+00, 2.30750E+00, &
   & 2.36165E+00, 2.44394E+00, 2.52782E+00, 2.61343E+00, 2.71640E+00, &
   & 2.81613E+00, 2.93679E+00, 3.01577E+00, 3.15995E+00, 3.15931E+00, &
   & 2.96658E+00, 2.73295E+00, 2.67480E+00, 2.66652E+00, 2.69393E+00, &
   & 2.75102E+00, 2.86503E+00, 2.99163E+00, 2.99576E+00, 3.02603E+00, &
   & 2.98415E+00, 2.79309E+00, 2.65337E+00, 2.50962E+00, 2.43207E+00, &
   & 2.34812E+00, 2.34872E+00, 2.35186E+00, 2.39477E+00, 2.42629E+00, &
   & 2.48068E+00, 2.55087E+00, 2.55952E+00, 2.56497E+00, 2.64323E+00, &
   & 2.67961E+00, 2.66263E+00, 2.70243E+00, 2.74911E+00, 2.81786E+00/
   DATA C00401 /                                                     &
   & 2.88684E+00, 2.97790E+00, 3.04305E+00, 3.13053E+00, 3.23857E+00, &
   & 3.35582E+00, 3.40654E+00, 3.38117E+00, 3.36296E+00, 3.39480E+00, &
   & 3.49066E+00, 3.60832E+00, 3.71817E+00, 3.83924E+00, 3.96355E+00, &
   & 4.03656E+00, 4.00518E+00, 3.90389E+00, 3.74790E+00, 3.61385E+00, &
   & 3.57066E+00, 3.59438E+00, 3.66182E+00, 3.71176E+00, 3.75255E+00, &
   & 3.79101E+00, 3.85278E+00, 3.85027E+00, 3.81112E+00, 3.72553E+00, &
   & 3.61017E+00, 3.54384E+00, 3.52406E+00, 3.54097E+00, 3.59375E+00, &
   & 3.66312E+00, 3.72632E+00, 3.76825E+00, 3.86798E+00, 3.92916E+00, &
   & 3.95610E+00, 4.00120E+00, 4.05865E+00, 4.11981E+00, 4.14634E+00, &
   & 4.19109E+00, 4.20317E+00, 4.25754E+00, 4.35131E+00, 4.48573E+00, &
   & 4.58716E+00, 4.67462E+00, 4.78228E+00, 4.91196E+00, 5.01871E+00, &
   & 5.10663E+00, 5.17780E+00, 5.21393E+00, 5.18144E+00, 5.04379E+00, &
   & 4.86504E+00, 4.78569E+00, 4.72717E+00, 4.69132E+00, 4.65797E+00, &
   & 4.60305E+00, 4.59798E+00, 4.65300E+00, 4.69707E+00, 4.74790E+00, &
   & 4.82581E+00, 4.80953E+00, 4.80517E+00, 4.82685E+00, 4.82321E+00, &
   & 4.84806E+00, 4.88591E+00, 4.91759E+00, 4.98074E+00, 5.07071E+00/
   DATA C00481 /                                                     &
   & 5.18733E+00, 5.30567E+00, 5.38670E+00, 5.43942E+00, 5.51797E+00, &
   & 5.62652E+00, 5.71228E+00, 5.82347E+00, 5.91434E+00, 6.00171E+00, &
   & 6.06977E+00, 6.13040E+00, 6.21990E+00, 6.29980E+00, 6.37206E+00, &
   & 6.48233E+00, 6.53068E+00, 6.53275E+00, 6.56858E+00, 6.54577E+00, &
   & 6.50472E+00, 6.41504E+00, 6.33853E+00, 6.31184E+00, 6.21253E+00, &
   & 6.22034E+00, 6.26918E+00, 6.28982E+00, 6.29461E+00, 6.35418E+00, &
   & 6.40956E+00, 6.38020E+00, 6.39784E+00, 6.45383E+00, 6.50134E+00, &
   & 6.56808E+00, 6.58850E+00, 6.58882E+00, 6.65097E+00, 6.75259E+00, &
   & 6.83256E+00, 6.92593E+00, 6.98083E+00, 7.03632E+00, 7.11147E+00, &
   & 7.15622E+00, 7.21106E+00, 7.27319E+00, 7.33382E+00, 7.38601E+00, &
   & 7.48971E+00, 7.61459E+00, 7.70134E+00, 7.76194E+00, 7.85534E+00, &
   & 7.99519E+00, 8.12227E+00, 8.25461E+00, 8.34670E+00, 8.42733E+00, &
   & 8.51806E+00, 8.57638E+00, 8.56481E+00, 8.55461E+00, 8.55593E+00, &
   & 8.58756E+00, 8.50070E+00, 8.54400E+00, 8.57575E+00, 8.62083E+00, &
   & 8.60684E+00, 8.67824E+00, 8.72069E+00, 8.79127E+00, 8.85479E+00, &
   & 8.86770E+00, 8.90574E+00, 8.91531E+00, 8.94800E+00, 9.00167E+00/
   DATA C00561 /                                                     &
   & 9.14051E+00, 9.25421E+00, 9.39694E+00, 9.50896E+00, 9.53190E+00, &
   & 9.55977E+00, 9.53482E+00, 9.49662E+00, 9.53359E+00, 9.54007E+00, &
   & 9.49809E+00, 9.49373E+00, 9.53282E+00, 9.63757E+00, 9.67855E+00, &
   & 9.67633E+00, 9.67045E+00, 9.79481E+00, 9.93420E+00, 1.00234E+01, &
   & 1.01372E+01, 1.02577E+01, 1.05056E+01, 1.07873E+01, 1.09967E+01, &
   & 1.10873E+01, 1.11624E+01, 1.13006E+01, 1.14875E+01, 1.16106E+01, &
   & 1.16744E+01, 1.17582E+01, 1.17709E+01, 1.18537E+01, 1.19623E+01, &
   & 1.19763E+01, 1.19879E+01, 1.20384E+01, 1.20763E+01, 1.20826E+01, &
   & 1.20449E+01, 1.19747E+01, 1.20227E+01, 1.21805E+01, 1.23134E+01, &
   & 1.24042E+01, 1.25614E+01, 1.26828E+01, 1.26645E+01, 1.26963E+01, &
   & 1.28226E+01, 1.28720E+01, 1.28981E+01, 1.29462E+01, 1.29363E+01, &
   & 1.29199E+01, 1.29797E+01, 1.28860E+01, 1.29126E+01, 1.30205E+01, &
   & 1.31327E+01, 1.31722E+01, 1.31901E+01, 1.33189E+01, 1.34833E+01, &
   & 1.36228E+01, 1.37474E+01, 1.38548E+01, 1.39450E+01, 1.40926E+01, &
   & 1.43099E+01, 1.44836E+01, 1.46257E+01, 1.47755E+01, 1.49163E+01, &
   & 1.51038E+01, 1.53308E+01, 1.54194E+01, 1.54852E+01, 1.55968E+01/
   DATA C00641 /                                                     &
   & 1.57025E+01, 1.58667E+01, 1.60365E+01, 1.61427E+01, 1.62967E+01, &
   & 1.64735E+01, 1.66123E+01, 1.67268E+01, 1.67673E+01, 1.67825E+01, &
   & 1.68898E+01, 1.68178E+01, 1.68216E+01, 1.68574E+01, 1.68799E+01, &
   & 1.70317E+01, 1.70767E+01, 1.71508E+01, 1.72965E+01, 1.73421E+01, &
   & 1.73937E+01, 1.74420E+01, 1.74535E+01, 1.75110E+01, 1.75497E+01, &
   & 1.75149E+01, 1.75955E+01, 1.78260E+01, 1.78271E+01, 1.79750E+01, &
   & 1.80600E+01, 1.81597E+01, 1.83454E+01, 1.85243E+01, 1.87382E+01, &
   & 1.88904E+01, 1.90395E+01, 1.92759E+01, 1.95398E+01, 1.97712E+01, &
   & 1.98487E+01, 1.99522E+01, 2.02363E+01, 2.03271E+01, 2.07090E+01, &
   & 2.09195E+01, 2.10974E+01, 2.11702E+01, 2.12964E+01, 2.14339E+01, &
   & 2.15764E+01, 2.17351E+01, 2.18486E+01, 2.19700E+01, 2.21663E+01, &
   & 2.24244E+01, 2.24813E+01, 2.25248E+01, 2.26357E+01, 2.26457E+01, &
   & 2.27249E+01, 2.27172E+01, 2.27123E+01, 2.26859E+01, 2.27216E+01, &
   & 2.29306E+01, 2.30711E+01, 2.31374E+01, 2.31815E+01, 2.33423E+01, &
   & 2.33810E+01, 2.36430E+01, 2.36807E+01, 2.36676E+01, 2.38607E+01, &
   & 2.41559E+01, 2.43413E+01, 2.44401E+01, 2.45968E+01, 2.47927E+01/
   DATA C00721 /                                                     &
   & 2.50743E+01, 2.53667E+01, 2.55749E+01, 2.57357E+01, 2.58927E+01, &
   & 2.61523E+01, 2.64110E+01, 2.66650E+01, 2.68829E+01, 2.70635E+01, &
   & 2.72797E+01, 2.75064E+01, 2.77229E+01, 2.80341E+01, 2.82003E+01, &
   & 2.83346E+01, 2.83909E+01, 2.86212E+01, 2.88006E+01, 2.89577E+01, &
   & 2.90965E+01, 2.91834E+01, 2.93224E+01, 2.94094E+01, 2.94848E+01, &
   & 2.96584E+01, 2.96749E+01, 2.97760E+01, 2.99163E+01, 3.00238E+01, &
   & 3.01290E+01, 3.02307E+01, 3.03663E+01, 3.05897E+01, 3.07937E+01, &
   & 3.10403E+01, 3.11778E+01, 3.13271E+01, 3.15799E+01, 3.18435E+01, &
   & 3.21614E+01, 3.25097E+01, 3.27701E+01, 3.29600E+01, 3.32583E+01, &
   & 3.36348E+01, 3.40282E+01, 3.41751E+01, 3.44128E+01, 3.46199E+01, &
   & 3.49363E+01, 3.52087E+01, 3.54056E+01, 3.55596E+01, 3.56694E+01, &
   & 3.58104E+01, 3.60276E+01, 3.62818E+01, 3.63505E+01, 3.66069E+01, &
   & 3.67544E+01, 3.70664E+01, 3.72525E+01, 3.73491E+01, 3.76006E+01, &
   & 3.77102E+01, 3.78970E+01, 3.81254E+01, 3.82728E+01, 3.81720E+01, &
   & 3.82781E+01, 3.84982E+01, 3.87202E+01, 3.89958E+01, 3.94148E+01, &
   & 3.98434E+01, 3.98952E+01, 4.01573E+01, 4.06014E+01, 4.09651E+01/
   DATA C00801 /                                                     &
   & 4.12821E+01, 4.16849E+01, 4.19899E+01, 4.22719E+01, 4.27736E+01, &
   & 4.32254E+01, 4.33883E+01, 4.39831E+01, 4.39414E+01, 4.42613E+01, &
   & 4.46503E+01, 4.49027E+01, 4.50384E+01, 4.52929E+01, 4.57269E+01, &
   & 4.56433E+01, 4.57350E+01, 4.60128E+01, 4.60487E+01, 4.61183E+01, &
   & 4.64397E+01, 4.68211E+01, 4.70706E+01, 4.72821E+01, 4.74972E+01, &
   & 4.78253E+01, 4.81615E+01, 4.84480E+01, 4.85703E+01, 4.87397E+01, &
   & 4.90015E+01, 4.93673E+01, 4.97291E+01, 4.99836E+01, 5.02975E+01, &
   & 5.05572E+01, 5.08226E+01, 5.13433E+01, 5.17112E+01, 5.19703E+01, &
   & 5.23128E+01, 5.27305E+01, 5.30599E+01, 5.34555E+01, 5.39625E+01, &
   & 5.43627E+01, 5.45446E+01, 5.49263E+01, 5.53511E+01, 5.57270E+01, &
   & 5.60904E+01, 5.63875E+01, 5.68475E+01, 5.73172E+01, 5.81134E+01, &
   & 5.86399E+01, 5.90384E+01, 5.91417E+01, 5.90883E+01, 5.93610E+01, &
   & 5.95794E+01, 5.99600E+01, 5.98493E+01, 5.99441E+01, 6.02748E+01, &
   & 6.04778E+01, 6.05233E+01, 6.07194E+01, 6.11589E+01, 6.13324E+01, &
   & 6.17685E+01, 6.23166E+01, 6.31055E+01, 6.38211E+01, 6.42320E+01, &
   & 6.45195E+01, 6.51125E+01, 6.56765E+01, 6.59286E+01, 6.62716E+01/
   DATA C00881 /                                                     &
   & 6.65693E+01, 6.68906E+01, 6.72246E+01, 6.75177E+01, 6.78476E+01, &
   & 6.82599E+01, 6.84400E+01, 6.89072E+01, 6.95720E+01, 7.01410E+01, &
   & 7.05519E+01, 7.09367E+01, 7.13975E+01, 7.22128E+01, 7.28222E+01, &
   & 7.33808E+01, 7.38828E+01, 7.44496E+01, 7.49983E+01, 7.54178E+01, &
   & 7.60554E+01, 7.62484E+01, 7.67892E+01, 7.71262E+01, 7.76235E+01, &
   & 7.81413E+01, 7.85694E+01, 7.91248E+01, 7.94715E+01, 7.96200E+01, &
   & 8.00270E+01, 8.03783E+01, 8.07100E+01, 8.11929E+01, 8.17375E+01, &
   & 8.18410E+01, 8.23341E+01, 8.26754E+01, 8.30893E+01, 8.34232E+01, &
   & 8.35533E+01, 8.36017E+01, 8.38589E+01, 8.43366E+01, 8.47593E+01, &
   & 8.51614E+01, 8.55271E+01, 8.58979E+01, 8.64892E+01, 8.74367E+01, &
   & 8.82440E+01, 8.89105E+01, 8.90980E+01, 8.97266E+01, 9.04886E+01, &
   & 9.12709E+01, 9.21243E+01, 9.26673E+01, 9.31331E+01, 9.38190E+01, &
   & 9.44877E+01, 9.50636E+01, 9.57445E+01, 9.65211E+01, 9.68623E+01, &
   & 9.75356E+01, 9.81991E+01, 9.88881E+01, 9.94554E+01, 9.99292E+01, &
   & 1.00357E+02, 1.00670E+02, 1.01227E+02, 1.01529E+02, 1.01889E+02, &
   & 1.02033E+02, 1.02254E+02, 1.02731E+02, 1.02914E+02, 1.03120E+02/
   DATA C00961 /                                                     &
   & 1.03674E+02, 1.03768E+02, 1.04146E+02, 1.04850E+02, 1.05525E+02, &
   & 1.06263E+02, 1.06653E+02, 1.07084E+02, 1.07461E+02, 1.08052E+02, &
   & 1.08793E+02, 1.09395E+02, 1.09811E+02, 1.10079E+02, 1.10656E+02, &
   & 1.11575E+02, 1.12544E+02, 1.13453E+02, 1.14440E+02, 1.15292E+02, &
   & 1.15869E+02, 1.16925E+02, 1.17854E+02, 1.18723E+02, 1.19574E+02, &
   & 1.19940E+02, 1.21108E+02, 1.21807E+02, 1.22490E+02, 1.23278E+02, &
   & 1.24094E+02, 1.24816E+02, 1.25469E+02, 1.26217E+02, 1.26878E+02, &
   & 1.27536E+02, 1.28168E+02, 1.28682E+02, 1.29076E+02, 1.30171E+02, &
   & 1.30667E+02, 1.31242E+02, 1.31665E+02, 1.31961E+02, 1.32347E+02, &
   & 1.32805E+02, 1.33152E+02, 1.33869E+02, 1.34261E+02, 1.34498E+02, &
   & 1.35028E+02, 1.36049E+02, 1.36577E+02, 1.37491E+02, 1.38078E+02, &
   & 1.38389E+02, 1.38819E+02, 1.39653E+02, 1.39770E+02, 1.40812E+02, &
   & 1.40926E+02, 1.41267E+02, 1.41872E+02, 1.42233E+02, 1.43447E+02, &
   & 1.44641E+02, 1.45500E+02, 1.45996E+02, 1.47040E+02, 1.48767E+02, &
   & 1.48785E+02, 1.49525E+02, 1.50266E+02, 1.50814E+02, 1.51443E+02, &
   & 1.52272E+02, 1.52846E+02, 1.54000E+02, 1.54629E+02, 1.54907E+02/
   DATA C01041 /                                                     &
   & 1.55527E+02, 1.56642E+02, 1.57436E+02, 1.59036E+02, 1.59336E+02, &
   & 1.59661E+02, 1.60287E+02, 1.61202E+02, 1.62410E+02, 1.63040E+02, &
   & 1.62872E+02, 1.63248E+02, 1.63776E+02, 1.64313E+02, 1.65782E+02, &
   & 1.65692E+02, 1.66049E+02, 1.66701E+02, 1.67786E+02, 1.69150E+02, &
   & 1.69996E+02, 1.71634E+02, 1.71137E+02, 1.71372E+02, 1.72525E+02, &
   & 1.73816E+02, 1.75219E+02, 1.76091E+02, 1.78260E+02, 1.79299E+02, &
   & 1.79904E+02, 1.81718E+02, 1.83807E+02, 1.85488E+02, 1.85929E+02, &
   & 1.86787E+02, 1.88282E+02, 1.89546E+02, 1.91489E+02, 1.92646E+02, &
   & 1.93399E+02, 1.93838E+02, 1.94406E+02, 1.95829E+02, 1.96745E+02, &
   & 1.96978E+02, 1.97243E+02, 1.97636E+02, 1.98025E+02, 1.98227E+02, &
   & 1.99552E+02, 2.00304E+02, 2.01031E+02, 2.01788E+02, 2.02432E+02, &
   & 2.03817E+02, 2.04866E+02, 2.05561E+02, 2.06180E+02, 2.07024E+02, &
   & 2.08303E+02, 2.09426E+02, 2.10575E+02, 2.11637E+02, 2.12559E+02, &
   & 2.13361E+02, 2.14191E+02, 2.15264E+02, 2.16366E+02, 2.17316E+02, &
   & 2.17717E+02, 2.17154E+02, 2.19172E+02, 2.20346E+02, 2.20849E+02, &
   & 2.21539E+02, 2.22810E+02, 2.22740E+02, 2.22824E+02, 2.23285E+02/
   DATA C01121 /                                                     &
   & 2.23696E+02, 2.23864E+02, 2.23968E+02, 2.23544E+02, 2.24804E+02, &
   & 2.25953E+02, 2.26753E+02, 2.27732E+02, 2.29505E+02, 2.30108E+02, &
   & 2.31232E+02, 2.32552E+02, 2.33979E+02, 2.36677E+02, 2.38481E+02, &
   & 2.41797E+02, 2.44025E+02, 2.45113E+02, 2.47373E+02, 2.47258E+02, &
   & 2.48617E+02, 2.49790E+02, 2.50562E+02, 2.51198E+02, 2.51289E+02, &
   & 2.52509E+02, 2.54136E+02, 2.55335E+02, 2.55808E+02, 2.56567E+02, &
   & 2.57977E+02, 2.58987E+02, 2.59622E+02, 2.60170E+02, 2.61127E+02, &
   & 2.60655E+02, 2.62129E+02, 2.64020E+02, 2.65659E+02, 2.67086E+02, &
   & 2.67615E+02, 2.69800E+02, 2.71452E+02, 2.73314E+02, 2.76972E+02, &
   & 2.78005E+02, 2.79815E+02, 2.81709E+02, 2.84043E+02, 2.87070E+02, &
   & 2.88842E+02, 2.90555E+02, 2.92401E+02, 2.94314E+02, 2.96074E+02, &
   & 2.97103E+02, 2.98037E+02, 2.98113E+02, 2.97705E+02, 2.97350E+02, &
   & 2.97329E+02, 2.97016E+02, 2.96752E+02, 2.96599E+02, 2.96637E+02, &
   & 2.97057E+02, 2.97585E+02, 2.98179E+02, 2.98997E+02, 3.00012E+02, &
   & 3.00806E+02, 3.00908E+02, 3.02369E+02, 3.04063E+02, 3.05325E+02, &
   & 3.06737E+02, 3.08066E+02, 3.09694E+02, 3.11530E+02, 3.13132E+02/
   DATA C01201 /                                                     &
   & 3.13296E+02, 3.15513E+02, 3.16887E+02, 3.17682E+02, 3.18296E+02, &
   & 3.18654E+02, 3.18912E+02, 3.19236E+02, 3.19626E+02, 3.20020E+02, &
   & 3.20186E+02, 3.20709E+02, 3.21628E+02, 3.22625E+02, 3.23504E+02, &
   & 3.25479E+02, 3.26825E+02, 3.28146E+02, 3.29404E+02, 3.30512E+02, &
   & 3.32634E+02, 3.34422E+02, 3.35602E+02, 3.36833E+02, 3.39372E+02, &
   & 3.43446E+02, 3.46374E+02, 3.48719E+02, 3.50881E+02, 3.53160E+02, &
   & 3.54890E+02, 3.57162E+02, 3.59284E+02, 3.60876E+02, 3.62295E+02, &
   & 3.63987E+02, 3.64835E+02, 3.65257E+02, 3.65738E+02, 3.65904E+02, &
   & 3.65976E+02, 3.66460E+02, 3.67087E+02, 3.67377E+02, 3.69079E+02, &
   & 3.70694E+02, 3.70940E+02, 3.70557E+02, 3.72693E+02, 3.73852E+02, &
   & 3.75679E+02, 3.77863E+02, 3.79964E+02, 3.81368E+02, 3.82716E+02, &
   & 3.85556E+02, 3.89072E+02, 3.91796E+02, 3.92766E+02, 3.96551E+02, &
   & 3.97833E+02, 3.97285E+02, 4.01929E+02, 4.02158E+02, 4.04553E+02, &
   & 4.06451E+02, 4.06236E+02, 4.08135E+02, 4.07797E+02, 4.08415E+02, &
   & 4.10111E+02, 4.11781E+02, 4.12735E+02, 4.11547E+02, 4.11606E+02, &
   & 4.13548E+02, 4.12557E+02, 4.12923E+02, 4.12866E+02, 4.13009E+02/
   DATA C01281 /                                                     &
   & 4.14447E+02, 4.16032E+02, 4.17032E+02, 4.19064E+02, 4.22458E+02, &
   & 4.26021E+02, 4.25192E+02, 4.25684E+02, 4.27536E+02, 4.29972E+02, &
   & 4.31994E+02, 4.36037E+02, 4.39132E+02, 4.40363E+02, 4.40716E+02, &
   & 4.40342E+02, 4.42063E+02, 4.44408E+02, 4.45454E+02, 4.47835E+02, &
   & 4.48256E+02, 4.48831E+02, 4.50257E+02, 4.51427E+02, 4.52373E+02, &
   & 4.53899E+02, 4.55496E+02, 4.56311E+02, 4.57314E+02, 4.59922E+02, &
   & 4.61048E+02, 4.59840E+02, 4.62144E+02, 4.63152E+02, 4.64565E+02, &
   & 4.66715E+02, 4.69380E+02, 4.70751E+02, 4.72012E+02, 4.73482E+02, &
   & 4.75524E+02, 4.79307E+02, 4.82035E+02, 4.84423E+02, 4.86712E+02, &
   & 4.88754E+02, 4.90102E+02, 4.92047E+02, 4.94150E+02, 4.95375E+02, &
   & 4.95828E+02, 4.97555E+02, 4.98559E+02, 4.97618E+02, 4.99265E+02, &
   & 4.99979E+02, 5.00681E+02, 5.01386E+02, 5.00868E+02, 5.01935E+02, &
   & 5.03151E+02, 5.04329E+02, 5.05546E+02, 5.08259E+02, 5.09222E+02, &
   & 5.09818E+02, 5.11397E+02, 5.12391E+02, 5.13326E+02, 5.14329E+02, &
   & 5.15443E+02, 5.16533E+02, 5.21417E+02, 5.25071E+02, 5.26581E+02, &
   & 5.27762E+02, 5.29274E+02, 5.31704E+02, 5.34310E+02, 5.35727E+02/
   DATA C01361 /                                                     &
   & 5.36838E+02, 5.37082E+02, 5.36733E+02, 5.36170E+02, 5.36063E+02, &
   & 5.36451E+02, 5.37870E+02, 5.40475E+02, 5.42268E+02, 5.41972E+02, &
   & 5.42532E+02, 5.44764E+02, 5.46844E+02, 5.47525E+02, 5.49150E+02, &
   & 5.52049E+02, 5.55423E+02, 5.56259E+02, 5.57424E+02, 5.59189E+02, &
   & 5.61167E+02, 5.64512E+02, 5.66753E+02, 5.68183E+02, 5.69628E+02, &
   & 5.73474E+02, 5.76192E+02, 5.78058E+02, 5.79588E+02, 5.81619E+02, &
   & 5.83530E+02, 5.84852E+02, 5.85326E+02, 5.88130E+02, 5.90570E+02, &
   & 5.91785E+02, 5.91371E+02, 5.90931E+02, 5.90942E+02, 5.91168E+02, &
   & 5.91291E+02, 5.89791E+02, 5.91146E+02, 5.90804E+02, 5.87847E+02, &
   & 5.89067E+02, 5.91027E+02, 5.90951E+02, 5.89227E+02, 5.93389E+02, &
   & 5.92921E+02, 5.92739E+02, 5.94544E+02, 5.98941E+02, 6.02302E+02, &
   & 6.03908E+02, 6.04265E+02, 6.06737E+02, 6.08560E+02, 6.11272E+02, &
   & 6.14992E+02, 6.18595E+02, 6.20930E+02, 6.22107E+02, 6.22957E+02, &
   & 6.26710E+02, 6.28657E+02, 6.30132E+02, 6.31543E+02, 6.33043E+02, &
   & 6.36932E+02, 6.38248E+02, 6.37126E+02, 6.41648E+02, 6.48274E+02, &
   & 6.52638E+02, 6.53922E+02, 6.56647E+02, 6.59351E+02, 6.60525E+02/
   DATA C01441 /                                                     &
   & 6.60130E+02, 6.61375E+02, 6.62660E+02, 6.63976E+02, 6.65181E+02, &
   & 6.64820E+02, 6.64458E+02, 6.64927E+02, 6.66555E+02, 6.66759E+02, &
   & 6.68218E+02, 6.70323E+02, 6.72703E+02, 6.76085E+02, 6.79180E+02, &
   & 6.80850E+02, 6.80017E+02, 6.79928E+02, 6.80886E+02, 6.82038E+02, &
   & 6.82271E+02, 6.84057E+02, 6.85309E+02, 6.86816E+02, 6.90180E+02, &
   & 6.93205E+02, 6.95870E+02, 6.98794E+02, 7.03776E+02, 7.04010E+02, &
   & 7.05041E+02, 7.07254E+02, 7.07432E+02, 7.10736E+02, 7.13791E+02, &
   & 7.15542E+02, 7.16468E+02, 7.17412E+02, 7.17783E+02, 7.17340E+02, &
   & 7.18184E+02, 7.18716E+02, 7.18809E+02, 7.18282E+02, 7.20317E+02, &
   & 7.18568E+02, 7.16274E+02, 7.19119E+02, 7.20852E+02, 7.21727E+02, &
   & 7.22607E+02, 7.26369E+02, 7.26412E+02, 7.27101E+02, 7.29404E+02, &
   & 7.30786E+02, 7.30910E+02, 7.30656E+02, 7.30566E+02, 7.33408E+02, &
   & 7.37064E+02, 7.39178E+02, 7.36713E+02, 7.37365E+02, 7.40861E+02, &
   & 7.45281E+02, 7.46178E+02, 7.46991E+02, 7.48035E+02, 7.49777E+02, &
   & 7.54665E+02, 7.56585E+02, 7.57408E+02, 7.58131E+02, 7.58155E+02, &
   & 7.60838E+02, 7.64792E+02, 7.68161E+02, 7.69263E+02, 7.73166E+02/
   DATA C01521 /                                                     &
   & 7.79006E+02, 7.82037E+02, 7.83109E+02, 7.84674E+02, 7.87444E+02, &
   & 7.89510E+02, 7.90130E+02, 7.91364E+02, 7.95225E+02, 8.03599E+02, &
   & 8.06340E+02, 8.05105E+02, 8.05120E+02, 8.08515E+02, 8.10907E+02, &
   & 8.11388E+02, 8.13432E+02, 8.12579E+02, 8.10564E+02, 8.08719E+02, &
   & 8.07682E+02, 8.05009E+02, 8.01754E+02, 8.01013E+02, 7.99926E+02, &
   & 7.99067E+02, 7.98369E+02, 7.94090E+02, 7.92883E+02, 7.94244E+02, &
   & 7.98220E+02, 7.98201E+02, 7.98332E+02, 7.99289E+02, 8.02355E+02, &
   & 8.03621E+02, 8.05302E+02, 8.08368E+02, 8.09983E+02, 8.11529E+02, &
   & 8.13068E+02, 8.14717E+02, 8.16441E+02, 8.19241E+02, 8.22944E+02, &
   & 8.23768E+02, 8.25030E+02, 8.26103E+02, 8.26374E+02, 8.28331E+02, &
   & 8.32620E+02, 8.38618E+02, 8.43666E+02, 8.45212E+02, 8.46324E+02, &
   & 8.48536E+02, 8.50192E+02, 8.53083E+02, 8.56653E+02, 8.59614E+02, &
   & 8.62000E+02, 8.64593E+02, 8.67678E+02, 8.70908E+02, 8.73408E+02, &
   & 8.74779E+02, 8.74005E+02, 8.76718E+02, 8.80445E+02, 8.84365E+02, &
   & 8.83806E+02, 8.84292E+02, 8.85539E+02, 8.87474E+02, 8.84905E+02, &
   & 8.84039E+02, 8.85105E+02, 8.83733E+02, 8.82224E+02, 8.79865E+02/
   DATA C01601 /                                                     &
   & 8.75663E+02, 8.75575E+02, 8.73144E+02, 8.68602E+02, 8.70278E+02, &
   & 8.69659E+02, 8.68701E+02, 8.69250E+02, 8.71057E+02, 8.72860E+02, &
   & 8.74361E+02, 8.74458E+02, 8.77576E+02, 8.81613E+02, 8.84358E+02, &
   & 8.87440E+02, 8.91549E+02, 8.96568E+02, 8.99836E+02, 9.02880E+02, &
   & 9.05428E+02, 9.06891E+02, 9.07349E+02, 9.10151E+02, 9.15917E+02, &
   & 9.16197E+02, 9.18571E+02, 9.21219E+02, 9.20292E+02, 9.21949E+02, &
   & 9.24509E+02, 9.27454E+02, 9.29474E+02, 9.31348E+02, 9.32818E+02, &
   & 9.32658E+02, 9.36280E+02, 9.39512E+02, 9.39667E+02, 9.44078E+02, &
   & 9.47196E+02, 9.48291E+02, 9.46150E+02, 9.46918E+02, 9.49093E+02, &
   & 9.51372E+02, 9.53109E+02, 9.56308E+02, 9.61335E+02, 9.58214E+02, &
   & 9.56188E+02, 9.55660E+02, 9.58633E+02, 9.57541E+02, 9.54879E+02, &
   & 9.51663E+02, 9.52839E+02, 9.52055E+02, 9.49253E+02, 9.50187E+02, &
   & 9.50323E+02, 9.50937E+02, 9.54362E+02, 9.55855E+02, 9.56350E+02, &
   & 9.55908E+02, 9.57963E+02, 9.61866E+02, 9.66948E+02, 9.69786E+02, &
   & 9.74302E+02, 9.79061E+02, 9.82465E+02, 9.86019E+02, 9.89930E+02, &
   & 9.94294E+02, 9.97011E+02, 9.98207E+02, 9.98607E+02, 1.00175E+03/
   DATA C01681 /                                                     &
   & 1.00275E+03, 1.00284E+03, 1.00294E+03, 1.00485E+03, 1.00593E+03, &
   & 1.00524E+03, 1.00415E+03, 1.00335E+03, 1.00278E+03, 1.00185E+03, &
   & 9.99982E+02, 9.98177E+02, 9.97959E+02, 9.99161E+02, 9.98810E+02, &
   & 9.95415E+02, 9.94342E+02, 9.92998E+02, 9.91340E+02, 9.90900E+02, &
   & 9.90407E+02, 9.89232E+02, 9.85447E+02, 9.86312E+02, 9.87461E+02, &
   & 9.86090E+02, 9.86670E+02, 9.85534E+02, 9.81877E+02, 9.84946E+02, &
   & 9.86392E+02, 9.86709E+02, 9.88086E+02, 9.90269E+02, 9.92566E+02, &
   & 9.94029E+02, 9.95795E+02, 9.97788E+02, 1.00005E+03, 1.00287E+03, &
   & 1.00566E+03, 1.00833E+03, 1.00982E+03, 1.01348E+03, 1.01862E+03, &
   & 1.02322E+03, 1.02786E+03, 1.03179E+03, 1.03339E+03, 1.03833E+03, &
   & 1.04317E+03, 1.04598E+03, 1.04753E+03, 1.04981E+03, 1.05321E+03, &
   & 1.05492E+03, 1.05721E+03, 1.05978E+03, 1.06033E+03, 1.06107E+03, &
   & 1.06155E+03, 1.06035E+03, 1.05838E+03, 1.05649E+03, 1.05553E+03, &
   & 1.05498E+03, 1.05387E+03, 1.05171E+03, 1.04877E+03, 1.04725E+03, &
   & 1.04748E+03, 1.04733E+03, 1.04704E+03, 1.04643E+03, 1.04411E+03, &
   & 1.04435E+03, 1.04520E+03, 1.04233E+03, 1.04047E+03, 1.03992E+03/
   DATA C01761 /                                                     &
   & 1.04192E+03, 1.04171E+03, 1.04140E+03, 1.04197E+03, 1.04415E+03, &
   & 1.04548E+03, 1.04533E+03, 1.04616E+03, 1.04705E+03, 1.04800E+03, &
   & 1.05025E+03, 1.05219E+03, 1.05412E+03, 1.05808E+03, 1.06062E+03, &
   & 1.06292E+03, 1.06780E+03, 1.07219E+03, 1.07610E+03, 1.07913E+03, &
   & 1.08405E+03, 1.08798E+03, 1.08835E+03, 1.09140E+03, 1.09447E+03, &
   & 1.09676E+03, 1.10015E+03, 1.10272E+03, 1.10410E+03, 1.10749E+03, &
   & 1.10991E+03, 1.11121E+03, 1.10981E+03, 1.10981E+03, 1.11063E+03, &
   & 1.10714E+03, 1.10500E+03, 1.10357E+03, 1.10093E+03, 1.09898E+03, &
   & 1.09679E+03, 1.09188E+03, 1.09088E+03, 1.09040E+03, 1.08586E+03, &
   & 1.08178E+03, 1.07752E+03, 1.07243E+03, 1.07178E+03, 1.07084E+03, &
   & 1.06693E+03, 1.06527E+03, 1.06405E+03, 1.06285E+03, 1.06287E+03, &
   & 1.06276E+03, 1.06221E+03, 1.06464E+03, 1.06579E+03, 1.06498E+03, &
   & 1.06596E+03, 1.06812E+03, 1.07159E+03, 1.07361E+03, 1.07556E+03, &
   & 1.07751E+03, 1.08128E+03, 1.08523E+03, 1.08927E+03, 1.09193E+03, &
   & 1.09612E+03, 1.10133E+03, 1.10435E+03, 1.10781E+03, 1.11168E+03, &
   & 1.11641E+03, 1.12217E+03, 1.12839E+03, 1.13298E+03, 1.13575E+03/
   DATA C01841 /                                                     &
   & 1.13742E+03, 1.13929E+03, 1.14132E+03, 1.14340E+03, 1.14518E+03, &
   & 1.14742E+03, 1.14943E+03, 1.14935E+03, 1.14975E+03, 1.15086E+03, &
   & 1.15420E+03, 1.15267E+03, 1.15007E+03, 1.15155E+03, 1.14982E+03, &
   & 1.14663E+03, 1.14301E+03, 1.13986E+03, 1.13676E+03, 1.13307E+03, &
   & 1.12898E+03, 1.12516E+03, 1.12284E+03, 1.12068E+03, 1.11855E+03, &
   & 1.11632E+03, 1.11464E+03, 1.11318E+03, 1.11180E+03, 1.11163E+03, &
   & 1.11160E+03, 1.11035E+03, 1.11178E+03, 1.11395E+03, 1.11447E+03, &
   & 1.11439E+03, 1.11440E+03, 1.11582E+03, 1.11560E+03, 1.11478E+03, &
   & 1.11448E+03, 1.11454E+03, 1.11494E+03, 1.11607E+03, 1.11736E+03, &
   & 1.11854E+03, 1.11875E+03, 1.11989E+03, 1.12165E+03, 1.12427E+03, &
   & 1.12620E+03, 1.12758E+03, 1.12774E+03, 1.12870E+03, 1.13001E+03, &
   & 1.13006E+03, 1.13078E+03, 1.13172E+03, 1.12971E+03, 1.12857E+03, &
   & 1.12810E+03, 1.12740E+03, 1.12659E+03, 1.12564E+03, 1.12338E+03, &
   & 1.12117E+03, 1.11902E+03, 1.11878E+03, 1.11855E+03, 1.11828E+03, &
   & 1.11791E+03, 1.11784E+03, 1.11815E+03, 1.11957E+03, 1.12046E+03, &
   & 1.12042E+03, 1.11929E+03, 1.12074E+03, 1.12708E+03, 1.12600E+03/
   DATA C01921 /                                                     &
   & 1.12538E+03, 1.12871E+03, 1.13167E+03, 1.13388E+03, 1.13444E+03, &
   & 1.13595E+03, 1.13801E+03, 1.14096E+03, 1.14230E+03, 1.14304E+03, &
   & 1.14421E+03, 1.14580E+03, 1.14767E+03, 1.15000E+03, 1.15126E+03, &
   & 1.15181E+03, 1.15197E+03, 1.15364E+03, 1.15626E+03, 1.15538E+03, &
   & 1.15636E+03, 1.15908E+03, 1.16024E+03, 1.16188E+03, 1.16411E+03, &
   & 1.16310E+03, 1.16430E+03, 1.16927E+03, 1.17035E+03, 1.17052E+03, &
   & 1.17013E+03, 1.16968E+03, 1.16969E+03, 1.17106E+03, 1.17123E+03, &
   & 1.17006E+03, 1.16536E+03, 1.16087E+03, 1.15691E+03, 1.15608E+03, &
   & 1.15388E+03, 1.15077E+03, 1.14967E+03, 1.14793E+03, 1.14554E+03, &
   & 1.14212E+03, 1.13908E+03, 1.13654E+03, 1.13499E+03, 1.13308E+03, &
   & 1.13033E+03, 1.13051E+03, 1.13073E+03, 1.12898E+03, 1.12941E+03, &
   & 1.13051E+03, 1.13086E+03, 1.13189E+03, 1.13304E+03, 1.13192E+03, &
   & 1.13131E+03, 1.13110E+03, 1.13499E+03, 1.13914E+03, 1.14359E+03, &
   & 1.14383E+03, 1.14390E+03, 1.14435E+03, 1.14540E+03, 1.14646E+03, &
   & 1.14716E+03, 1.14880E+03, 1.15062E+03, 1.15170E+03, 1.15093E+03, &
   & 1.14926E+03, 1.15133E+03, 1.15167E+03, 1.15043E+03, 1.15134E+03/
   DATA C02001 /                                                     &
   & 1.15135E+03, 1.15000E+03, 1.15087E+03, 1.15118E+03, 1.14935E+03, &
   & 1.14780E+03, 1.14647E+03, 1.14560E+03, 1.14404E+03, 1.14238E+03, &
   & 1.14406E+03, 1.14245E+03, 1.13781E+03, 1.13664E+03, 1.13653E+03, &
   & 1.13778E+03, 1.13813E+03, 1.13794E+03, 1.13681E+03, 1.13515E+03, &
   & 1.13328E+03, 1.13132E+03, 1.13080E+03, 1.13130E+03, 1.13400E+03, &
   & 1.13526E+03, 1.13494E+03, 1.13193E+03, 1.12898E+03, 1.12654E+03, &
   & 1.12739E+03, 1.12849E+03, 1.12774E+03, 1.12733E+03, 1.12733E+03, &
   & 1.12943E+03, 1.13014E+03, 1.12967E+03, 1.12731E+03, 1.12671E+03, &
   & 1.12885E+03, 1.13050E+03, 1.13201E+03, 1.13345E+03, 1.13488E+03, &
   & 1.13605E+03, 1.13530E+03, 1.13737E+03, 1.14186E+03, 1.14250E+03, &
   & 1.14305E+03, 1.14383E+03, 1.14510E+03, 1.14659E+03, 1.14848E+03, &
   & 1.14949E+03, 1.14995E+03, 1.14934E+03, 1.15058E+03, 1.15368E+03, &
   & 1.15435E+03, 1.15422E+03, 1.15296E+03, 1.15228E+03, 1.15189E+03, &
   & 1.15198E+03, 1.15081E+03, 1.14881E+03, 1.14562E+03, 1.14276E+03, &
   & 1.14030E+03, 1.13637E+03, 1.13254E+03, 1.12942E+03, 1.12653E+03, &
   & 1.12362E+03, 1.11987E+03, 1.11712E+03, 1.11522E+03, 1.11403E+03/
   DATA C02081 /                                                     &
   & 1.11226E+03, 1.10947E+03, 1.10956E+03, 1.10976E+03, 1.10748E+03, &
   & 1.10673E+03, 1.10688E+03, 1.10675E+03, 1.10533E+03, 1.10230E+03, &
   & 1.10384E+03, 1.10496E+03, 1.10274E+03, 1.10197E+03, 1.10196E+03, &
   & 1.10278E+03, 1.10257E+03, 1.10147E+03, 1.10205E+03, 1.10308E+03, &
   & 1.10478E+03, 1.10358E+03, 1.10197E+03, 1.10305E+03, 1.10390E+03, &
   & 1.10456E+03, 1.10526E+03, 1.10588E+03, 1.10640E+03, 1.10747E+03, &
   & 1.10904E+03, 1.11214E+03, 1.11350E+03, 1.11359E+03, 1.11604E+03, &
   & 1.11706E+03, 1.11594E+03, 1.11600E+03, 1.11616E+03, 1.11561E+03, &
   & 1.11556E+03, 1.11547E+03, 1.11370E+03, 1.11289E+03, 1.11276E+03, &
   & 1.11338E+03, 1.11437E+03, 1.11595E+03, 1.11309E+03, 1.10958E+03, &
   & 1.10887E+03, 1.10573E+03, 1.10068E+03, 1.10194E+03, 1.10165E+03, &
   & 1.09813E+03, 1.09973E+03, 1.10233E+03, 1.10121E+03, 1.10097E+03, &
   & 1.10149E+03, 1.10162E+03, 1.10222E+03, 1.10389E+03, 1.10315E+03, &
   & 1.10158E+03, 1.10193E+03, 1.10186E+03, 1.10135E+03, 1.10336E+03, &
   & 1.10500E+03, 1.10459E+03, 1.10592E+03, 1.10784E+03, 1.10076E+03, &
   & 1.09615E+03, 1.09496E+03, 1.09422E+03, 1.09350E+03, 1.09244E+03/
   DATA C02161 /                                                     &
   & 1.08955E+03, 1.08535E+03, 1.08379E+03, 1.08184E+03, 1.07889E+03, &
   & 1.07563E+03, 1.07238E+03, 1.07042E+03, 1.06882E+03, 1.06761E+03, &
   & 1.06816E+03, 1.06772E+03, 1.06327E+03, 1.06313E+03, 1.06563E+03, &
   & 1.06254E+03, 1.06072E+03, 1.06095E+03, 1.06173E+03, 1.06269E+03, &
   & 1.06361E+03, 1.06438E+03, 1.06501E+03, 1.06465E+03, 1.06481E+03, &
   & 1.06685E+03, 1.06642E+03, 1.06447E+03, 1.06701E+03, 1.06791E+03, &
   & 1.06612E+03, 1.06471E+03, 1.06403E+03, 1.06774E+03, 1.06823E+03, &
   & 1.06524E+03, 1.06479E+03, 1.06453E+03, 1.06346E+03, 1.06175E+03, &
   & 1.05958E+03, 1.05941E+03, 1.05936E+03, 1.05938E+03, 1.05736E+03, &
   & 1.05449E+03, 1.05307E+03, 1.05180E+03, 1.05074E+03, 1.04810E+03, &
   & 1.04536E+03, 1.04477E+03, 1.04389E+03, 1.04272E+03, 1.04006E+03, &
   & 1.03739E+03, 1.03533E+03, 1.03476E+03, 1.03516E+03, 1.03275E+03, &
   & 1.03093E+03, 1.03062E+03, 1.02997E+03, 1.02919E+03, 1.02993E+03, &
   & 1.02983E+03, 1.02837E+03, 1.02611E+03, 1.02386E+03, 1.02426E+03, &
   & 1.02542E+03, 1.02750E+03, 1.02638E+03, 1.02496E+03, 1.02608E+03, &
   & 1.02568E+03, 1.02388E+03, 1.02522E+03, 1.02692E+03, 1.02834E+03/
   DATA C02241 /                                                     &
   & 1.02828E+03, 1.02716E+03, 1.02667E+03, 1.02607E+03, 1.02503E+03, &
   & 1.02723E+03, 1.03143E+03, 1.02881E+03, 1.02646E+03, 1.02500E+03, &
   & 1.02569E+03, 1.02743E+03, 1.02608E+03, 1.02548E+03, 1.02620E+03, &
   & 1.02733E+03, 1.02839E+03, 1.02575E+03, 1.02432E+03, 1.02471E+03, &
   & 1.02392E+03, 1.02267E+03, 1.02077E+03, 1.01964E+03, 1.01957E+03, &
   & 1.01848E+03, 1.01704E+03, 1.01524E+03, 1.01352E+03, 1.01191E+03, &
   & 1.01066E+03, 1.00952E+03, 1.00849E+03, 1.00660E+03, 1.00368E+03, &
   & 9.99713E+02, 9.95921E+02, 9.94845E+02, 9.93286E+02, 9.91204E+02/
!
end block data BO3HH0
!
!     --------------------------------------------------------------
!
SUBROUTINE O3HHT1 (V1C,V2C,DVC,NPTC,C)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   COMMON /O3HH1/ V1S,V2S,DVS,NPTS,S(2687)
   DIMENSION C(*)
!
   DVC = DVS
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      C(J) = S(I)
10 END DO
!
   RETURN
!
end subroutine O3HHT1
!
!     --------------------------------------------------------------
!
BLOCK DATA BO3HH1
!
   IMPLICIT REAL*8           (V)
!
!     RATIO (C1/C0)
!     DATA FROM BASS 1985
!
!     NOW INCLUDES MOLINA & MOLINA AT 273K WITH THE TEMPERATURE
!     DEPENDENCE DETERMINED FROM THE 195K HARVARD MEASUREMENTS,
!     EMPLOYING THE BASS ALGORITHM
!
!              (CO(1+C1*(T-273.15)+C2*(T-273.15)**2);
!
!     THIS IS ONLY FOR THE WAVELENGTH RANGE FROM .34 TO .35 MICRONS;
!     OTHERWISE, THE BASS DATA ALONE HAVE BEEN EMPLOYED BETWEEN
!     .34 AND .245 MICRONS.
!
!     NEW T-DEPENDENT X-SECTIONS BETWEEN .345 AND .36 MICRONS
!     HAVE NOW BEEN ADDED, BASED ON WORK BY CACCIANI, DISARRA
!     AND FIOCCO, UNIVERSITY OF ROME, 1987.  QUADRATIC TEMP
!     HAS BEEN DERIVED, AS ABOVE.
!
!     AGREEMENT AMONGST THE FOUR DATA SETS IS REASONABLE (<10%)
!     AND OFTEN EXCELLENT (0-3%)
!
!
   COMMON /O3HH1/ V1C,V2C,DVC,NC,                                    &
   &               O31001(85),C10086(80),C10166(80),C10246(65),       &
   &               C10311(16),C10327(80),C10407( 1),                  &
   &               C10001(80),C10081(80),C10161(80),C10241(80),       &
   &               C10321(80),C10401(80),C10481(80),C10561(80),       &
   &               C10641(80),C10721(80),C10801(80),C10881(80),       &
   &               C10961(80),C11041(80),C11121(80),C11201(80),       &
   &               C11281(80),C11361(80),C11441(80),C11521(80),       &
   &               C11601(80),C11681(80),C11761(80),C11841(80),       &
   &               C11921(80),C12001(80),C12081(80),C12161(80),       &
   &               C12241(40)
!
!     DATA V1C /29405./, V2C /40800./ ,DVC /5./, NC /2280/   BASS
!
   DATA V1C /27370./, V2C /40800./ ,DVC /5./, NC /2687/
!
   DATA O31001/85*1.3E-3/
!
   DATA C10086/                                                      &
   & 1.37330E-03, 1.62821E-03, 2.01703E-03, 2.54574E-03, 3.20275E-03, &
   & 3.89777E-03, 4.62165E-03, 5.26292E-03, 5.86986E-03, 6.41494E-03, &
   & 6.96761E-03, 7.48539E-03, 7.89600E-03, 7.87305E-03, 7.81981E-03, &
   & 7.63864E-03, 7.67455E-03, 7.72586E-03, 7.69784E-03, 7.57367E-03, &
   & 7.27336E-03, 7.14064E-03, 7.24207E-03, 7.09851E-03, 6.93654E-03, &
   & 6.89385E-03, 7.05768E-03, 6.85578E-03, 6.58301E-03, 6.50848E-03, &
   & 6.52083E-03, 6.46590E-03, 6.70692E-03, 6.92053E-03, 7.17734E-03, &
   & 7.05364E-03, 6.63440E-03, 6.54702E-03, 6.27173E-03, 5.98150E-03, &
   & 5.66579E-03, 5.51549E-03, 5.50291E-03, 5.93271E-03, 6.36950E-03, &
   & 7.18562E-03, 7.51767E-03, 6.53815E-03, 7.22341E-03, 8.63056E-03, &
   & 9.11740E-03, 8.80903E-03, 8.59902E-03, 7.74287E-03, 7.33509E-03, &
   & 7.50180E-03, 7.81686E-03, 7.85635E-03, 8.08554E-03, 7.21968E-03, &
   & 7.99028E-03, 9.90724E-03, 1.29121E-02, 1.54686E-02, 1.60876E-02, &
   & 1.59530E-02, 1.57040E-02, 1.59499E-02, 1.63961E-02, 1.72670E-02, &
   & 1.81634E-02, 1.95519E-02, 2.14181E-02, 2.28670E-02, 2.33506E-02, &
   & 2.22736E-02, 2.14296E-02, 2.15271E-02, 2.30730E-02, 2.36220E-02/
   DATA C10166/                                                      &
   & 2.44466E-02, 2.44476E-02, 2.39223E-02, 2.41386E-02, 2.53687E-02, &
   & 2.67491E-02, 2.80425E-02, 2.77558E-02, 2.82626E-02, 2.86776E-02, &
   & 2.88781E-02, 2.89248E-02, 2.89983E-02, 2.85534E-02, 2.87102E-02, &
   & 2.83695E-02, 2.76719E-02, 2.76091E-02, 2.90733E-02, 2.80388E-02, &
   & 2.73706E-02, 2.65055E-02, 2.61268E-02, 2.45892E-02, 2.37213E-02, &
   & 2.22542E-02, 2.10116E-02, 2.02852E-02, 1.97635E-02, 1.94079E-02, &
   & 1.90997E-02, 1.85598E-02, 1.79221E-02, 1.77887E-02, 1.73709E-02, &
   & 1.67263E-02, 1.60932E-02, 1.50775E-02, 1.39563E-02, 1.23691E-02, &
   & 1.07402E-02, 9.35859E-03, 8.43786E-03, 7.92075E-03, 7.33239E-03, &
   & 6.73638E-03, 6.28740E-03, 5.85640E-03, 5.85384E-03, 6.10577E-03, &
   & 7.26050E-03, 9.66384E-03, 1.29629E-02, 1.69596E-02, 2.03465E-02, &
   & 2.26429E-02, 2.39653E-02, 2.47970E-02, 2.51993E-02, 2.51383E-02, &
   & 2.52014E-02, 2.47766E-02, 2.47171E-02, 2.47478E-02, 2.43986E-02, &
   & 2.43498E-02, 2.40537E-02, 2.40574E-02, 2.40446E-02, 2.40847E-02, &
   & 2.39400E-02, 2.42127E-02, 2.47123E-02, 2.52914E-02, 2.52103E-02, &
   & 2.51421E-02, 2.43229E-02, 2.37902E-02, 2.30865E-02, 2.28174E-02/
   DATA C10246/                                                      &
   & 2.28830E-02, 2.33671E-02, 2.38274E-02, 2.46699E-02, 2.56739E-02, &
   & 2.61408E-02, 2.62898E-02, 2.64228E-02, 2.55561E-02, 2.47095E-02, &
   & 2.39071E-02, 2.34319E-02, 2.28738E-02, 2.23434E-02, 2.18888E-02, &
   & 2.13639E-02, 2.11937E-02, 2.10110E-02, 2.07672E-02, 2.00697E-02, &
   & 1.97605E-02, 1.91208E-02, 1.82056E-02, 1.73945E-02, 1.64542E-02, &
   & 1.53969E-02, 1.41816E-02, 1.35665E-02, 1.27109E-02, 1.18254E-02, &
   & 1.11489E-02, 1.03984E-02, 1.00760E-02, 9.86649E-03, 9.76766E-03, &
   & 9.41662E-03, 9.19082E-03, 9.44272E-03, 1.04547E-02, 1.24713E-02, &
   & 1.49310E-02, 1.70272E-02, 1.86057E-02, 1.93555E-02, 1.98350E-02, &
   & 2.00041E-02, 2.01233E-02, 2.01917E-02, 1.98918E-02, 1.96649E-02, &
   & 1.95162E-02, 2.01044E-02, 2.06711E-02, 2.08881E-02, 2.04812E-02, &
   & 1.92249E-02, 1.80188E-02, 1.69496E-02, 1.60488E-02, 1.52865E-02, &
   & 1.46940E-02, 1.41067E-02, 1.35675E-02, 1.31094E-02, 1.27542E-02/
   DATA C10311/                                                      &
   &                                                     1.3073E-02,  &
   & 1.2795E-02,  1.2753E-02,  1.2868E-02,  1.2885E-02,  1.2554E-02,  &
   & 1.2106E-02,  1.1616E-02,  1.1394E-02,  1.1092E-02,  1.0682E-02,  &
   & 1.0519E-02,  9.7219E-03,  9.3434E-03,  8.5260E-03,  8.3333E-03/
   DATA C10327/                                                      &
   & 7.8582E-03,  6.8295E-03,  6.7963E-03,  6.7516E-03,  6.2930E-03,  &
   & 6.1615E-03,  6.1250E-03,  5.9011E-03,  5.7823E-03,  5.4688E-03,  &
   & 5.0978E-03,  4.4526E-03,  3.8090E-03,  3.2310E-03,  3.0128E-03,  &
   & 3.9063E-03,  6.7911E-03,  9.3161E-03,  1.0256E-02,  1.0183E-02,  &
   & 9.8289E-03,  9.5683E-03,  9.0406E-03,  8.7148E-03,  8.5284E-03,  &
   & 8.6149E-03,  8.7238E-03,  9.3679E-03,  1.0683E-02,  1.2016E-02,  &
   & 1.3097E-02,  1.3610E-02,  1.3588E-02,  1.3805E-02,  1.3928E-02,  &
   & 1.3903E-02,  1.3446E-02,  1.3258E-02,  1.3194E-02,  1.2703E-02,  &
   & 1.2393E-02,  1.2487E-02,  1.2341E-02,  1.2388E-02,  1.2061E-02,  &
   & 1.2122E-02,  1.1850E-02,  1.2032E-02,  1.1806E-02,  1.1810E-02,  &
   & 1.1572E-02,  1.1397E-02,  1.0980E-02,  1.1012E-02,  1.0524E-02,  &
   & 1.0518E-02,  1.0227E-02,  9.6837E-03,  9.6425E-03,  8.9938E-03,  &
   & 9.1488E-03,  8.8595E-03,  8.5976E-03,  8.4447E-03,  8.0731E-03,  &
   & 8.0283E-03,  7.7827E-03,  7.7638E-03,  7.2438E-03,  6.8246E-03,  &
   & 6.3457E-03,  5.6632E-03,  5.2500E-03,  4.3593E-03,  3.9431E-03,  &
   & 3.1580E-03,  2.2298E-03,  1.7818E-03,  1.4513E-03,  1.3188E-03/
   DATA C10407/                                                      &
   & 2.1034E-03/
   DATA C10001 /                                                     &
   & 6.45621E-03, 7.11308E-03, 1.06130E-02, 1.36338E-02, 1.27746E-02, &
   & 1.42152E-02, 1.41144E-02, 1.64830E-02, 1.67110E-02, 1.57368E-02, &
   & 1.54644E-02, 1.45248E-02, 1.43206E-02, 1.56946E-02, 1.54268E-02, &
   & 1.37500E-02, 1.50224E-02, 1.60919E-02, 1.49099E-02, 1.53960E-02, &
   & 1.61871E-02, 1.46539E-02, 1.38258E-02, 1.32571E-02, 1.21580E-02, &
   & 1.39596E-02, 1.16029E-02, 1.47042E-02, 1.07441E-02, 1.08999E-02, &
   & 1.05562E-02, 1.00589E-02, 9.60711E-03, 9.36950E-03, 7.65303E-03, &
   & 6.86216E-03, 7.05344E-03, 6.90728E-03, 6.78627E-03, 6.97435E-03, &
   & 5.75456E-03, 5.81685E-03, 5.00915E-03, 4.90259E-03, 4.42545E-03, &
   & 4.14633E-03, 3.61657E-03, 3.08178E-03, 2.91680E-03, 2.94554E-03, &
   & 3.35794E-03, 5.49025E-03, 7.09867E-03, 6.82592E-03, 8.84835E-03, &
   & 9.15718E-03, 9.17935E-03, 8.31848E-03, 7.79481E-03, 7.75125E-03, &
   & 6.95844E-03, 7.34506E-03, 7.53823E-03, 7.03272E-03, 7.57051E-03, &
   & 9.20239E-03, 1.10864E-02, 1.16188E-02, 1.30029E-02, 1.44364E-02, &
   & 1.29292E-02, 1.36031E-02, 1.35967E-02, 1.30412E-02, 1.29874E-02, &
   & 1.14829E-02, 1.18009E-02, 1.20829E-02, 1.17831E-02, 1.21489E-02/
   DATA C10081 /                                                     &
   & 1.27019E-02, 1.25557E-02, 1.23812E-02, 1.20158E-02, 1.26749E-02, &
   & 1.17139E-02, 1.14552E-02, 1.11268E-02, 9.79143E-03, 8.79741E-03, &
   & 8.85709E-03, 8.57653E-03, 8.93908E-03, 8.46205E-03, 8.56506E-03, &
   & 8.14319E-03, 8.14415E-03, 7.74205E-03, 7.80727E-03, 7.49886E-03, &
   & 7.71114E-03, 6.55963E-03, 6.87550E-03, 6.39162E-03, 5.55359E-03, &
   & 5.43275E-03, 4.90649E-03, 4.41165E-03, 4.21875E-03, 3.62592E-03, &
   & 3.40700E-03, 2.40267E-03, 2.61479E-03, 2.75677E-03, 4.10842E-03, &
   & 5.79601E-03, 7.10708E-03, 8.07826E-03, 8.16166E-03, 8.72620E-03, &
   & 8.85878E-03, 8.72755E-03, 8.25811E-03, 8.12100E-03, 7.78534E-03, &
   & 7.39762E-03, 8.43880E-03, 8.53789E-03, 9.90072E-03, 1.01668E-02, &
   & 1.00827E-02, 9.73556E-03, 9.57462E-03, 1.01289E-02, 1.10670E-02, &
   & 1.03508E-02, 1.00929E-02, 9.10236E-03, 9.39459E-03, 8.79601E-03, &
   & 8.67936E-03, 8.53862E-03, 7.95459E-03, 8.04037E-03, 7.95361E-03, &
   & 7.87432E-03, 6.99165E-03, 7.37107E-03, 6.09187E-03, 6.21030E-03, &
   & 5.33277E-03, 5.04633E-03, 4.45811E-03, 4.34153E-03, 3.98596E-03, &
   & 3.84225E-03, 3.41943E-03, 3.60535E-03, 2.81691E-03, 2.49771E-03/
   DATA C10161 /                                                     &
   & 2.35046E-03, 2.50947E-03, 3.75462E-03, 4.92349E-03, 5.09294E-03, &
   & 4.98312E-03, 5.19325E-03, 4.41827E-03, 4.25192E-03, 4.46745E-03, &
   & 4.08731E-03, 3.84776E-03, 3.67507E-03, 3.76845E-03, 3.69210E-03, &
   & 4.59864E-03, 6.42677E-03, 7.83255E-03, 7.89247E-03, 8.10883E-03, &
   & 8.00825E-03, 8.40322E-03, 7.97108E-03, 8.24714E-03, 8.39006E-03, &
   & 8.68787E-03, 8.61108E-03, 8.81552E-03, 9.36996E-03, 9.08243E-03, &
   & 9.69116E-03, 9.66185E-03, 9.22856E-03, 9.65086E-03, 9.35398E-03, &
   & 9.06358E-03, 8.76851E-03, 8.43072E-03, 7.85659E-03, 7.93936E-03, &
   & 7.49712E-03, 7.20199E-03, 6.94581E-03, 6.64086E-03, 6.12627E-03, &
   & 6.13967E-03, 5.67310E-03, 5.09928E-03, 4.59112E-03, 3.95257E-03, &
   & 3.67652E-03, 3.28781E-03, 2.77471E-03, 2.74494E-03, 2.15529E-03, &
   & 1.95283E-03, 1.75043E-03, 1.60419E-03, 1.82688E-03, 2.34667E-03, &
   & 2.92502E-03, 3.88322E-03, 4.39984E-03, 4.67814E-03, 4.80395E-03, &
   & 4.69130E-03, 4.54564E-03, 4.46773E-03, 4.59178E-03, 4.37498E-03, &
   & 4.12706E-03, 4.18299E-03, 4.57267E-03, 5.60127E-03, 6.51936E-03, &
   & 7.10498E-03, 7.49870E-03, 7.89554E-03, 7.97428E-03, 8.21044E-03/
   DATA C10241 /                                                     &
   & 8.06324E-03, 7.76648E-03, 7.62238E-03, 7.77675E-03, 7.46905E-03, &
   & 7.61115E-03, 7.42715E-03, 7.28461E-03, 7.51514E-03, 7.38782E-03, &
   & 6.97206E-03, 6.52738E-03, 6.10147E-03, 5.87553E-03, 5.49218E-03, &
   & 4.94873E-03, 4.47920E-03, 4.25005E-03, 3.98094E-03, 3.92084E-03, &
   & 3.41707E-03, 3.30501E-03, 3.09208E-03, 3.19686E-03, 3.55283E-03, &
   & 4.20775E-03, 4.11155E-03, 3.72193E-03, 3.52000E-03, 3.13572E-03, &
   & 2.87629E-03, 2.64251E-03, 2.33451E-03, 2.22426E-03, 2.05800E-03, &
   & 1.75214E-03, 2.32530E-03, 2.68651E-03, 3.66315E-03, 4.93904E-03, &
   & 5.32850E-03, 5.43978E-03, 5.32656E-03, 5.15649E-03, 5.42096E-03, &
   & 5.37193E-03, 5.23454E-03, 5.34557E-03, 5.50533E-03, 6.13216E-03, &
   & 6.65129E-03, 7.09357E-03, 7.46042E-03, 7.68605E-03, 7.91866E-03, &
   & 7.52953E-03, 7.48272E-03, 7.17800E-03, 6.80060E-03, 6.60427E-03, &
   & 6.43049E-03, 6.45975E-03, 6.20534E-03, 5.93094E-03, 5.67360E-03, &
   & 5.38584E-03, 5.19364E-03, 4.92599E-03, 4.60655E-03, 4.24669E-03, &
   & 3.94253E-03, 3.55894E-03, 3.24256E-03, 2.92974E-03, 2.62760E-03, &
   & 2.52238E-03, 2.24714E-03, 2.26350E-03, 2.44380E-03, 3.03798E-03/
   DATA C10321 /                                                     &
   & 3.50000E-03, 3.55416E-03, 3.43661E-03, 3.19814E-03, 3.02155E-03, &
   & 2.73890E-03, 2.50078E-03, 2.34595E-03, 2.18282E-03, 2.19285E-03, &
   & 2.49482E-03, 3.13434E-03, 4.18947E-03, 4.72069E-03, 5.29712E-03, &
   & 5.39004E-03, 5.44846E-03, 5.37952E-03, 5.09935E-03, 5.08741E-03, &
   & 5.05257E-03, 5.10339E-03, 5.17968E-03, 5.31841E-03, 5.58106E-03, &
   & 5.65031E-03, 5.65680E-03, 5.76184E-03, 5.71213E-03, 5.48515E-03, &
   & 5.32168E-03, 5.18505E-03, 4.99640E-03, 4.78746E-03, 4.57244E-03, &
   & 4.32728E-03, 4.14464E-03, 3.97659E-03, 4.01874E-03, 4.10588E-03, &
   & 3.99644E-03, 3.84584E-03, 3.64222E-03, 3.39590E-03, 3.00386E-03, &
   & 2.73790E-03, 2.45095E-03, 2.29068E-03, 1.64530E-03, 1.68602E-03, &
   & 2.32934E-03, 3.14851E-03, 3.65706E-03, 3.70878E-03, 3.75103E-03, &
   & 3.79183E-03, 3.32032E-03, 2.42604E-03, 2.48775E-03, 2.34603E-03, &
   & 2.36349E-03, 3.33744E-03, 3.44617E-03, 4.27280E-03, 4.61076E-03, &
   & 5.20165E-03, 5.14851E-03, 5.22909E-03, 5.08278E-03, 5.16125E-03, &
   & 5.01572E-03, 4.51685E-03, 4.67541E-03, 4.83421E-03, 4.57546E-03, &
   & 4.55111E-03, 5.03093E-03, 4.67838E-03, 4.44282E-03, 4.40774E-03/
   DATA C10401 /                                                     &
   & 4.48123E-03, 4.24410E-03, 4.03559E-03, 3.73969E-03, 3.45458E-03, &
   & 3.18217E-03, 3.16115E-03, 3.36877E-03, 3.62026E-03, 3.69898E-03, &
   & 3.49845E-03, 3.13839E-03, 2.77731E-03, 2.40106E-03, 2.03935E-03, &
   & 1.84377E-03, 2.07757E-03, 2.39550E-03, 2.86272E-03, 3.27900E-03, &
   & 3.42304E-03, 3.50211E-03, 3.29197E-03, 3.24784E-03, 3.20864E-03, &
   & 3.28063E-03, 3.01328E-03, 3.00379E-03, 3.19562E-03, 3.45113E-03, &
   & 3.75149E-03, 3.98520E-03, 4.19181E-03, 4.15773E-03, 4.02490E-03, &
   & 3.95936E-03, 3.79001E-03, 3.77647E-03, 3.48528E-03, 3.55768E-03, &
   & 3.62812E-03, 3.48650E-03, 3.35434E-03, 3.20088E-03, 3.25316E-03, &
   & 3.04467E-03, 3.12633E-03, 3.23602E-03, 3.07723E-03, 2.80070E-03, &
   & 2.72498E-03, 2.74752E-03, 2.58943E-03, 2.32482E-03, 2.20218E-03, &
   & 2.10846E-03, 2.05991E-03, 2.01844E-03, 2.16224E-03, 2.48456E-03, &
   & 2.88022E-03, 2.93939E-03, 3.01176E-03, 2.98886E-03, 2.96947E-03, &
   & 3.38082E-03, 3.61657E-03, 3.42654E-03, 3.41274E-03, 3.22475E-03, &
   & 2.97658E-03, 3.21944E-03, 3.32032E-03, 3.33273E-03, 3.58854E-03, &
   & 3.67023E-03, 3.64069E-03, 3.74557E-03, 3.77703E-03, 3.64042E-03/
   DATA C10481 /                                                     &
   & 3.39468E-03, 3.22657E-03, 3.16466E-03, 3.24224E-03, 3.24801E-03, &
   & 3.19487E-03, 3.40155E-03, 3.16940E-03, 2.92293E-03, 3.00998E-03, &
   & 2.82851E-03, 2.60381E-03, 2.59242E-03, 2.48530E-03, 2.76677E-03, &
   & 2.45506E-03, 2.21845E-03, 2.30407E-03, 2.28136E-03, 2.37278E-03, &
   & 2.25313E-03, 2.47836E-03, 2.77858E-03, 2.89803E-03, 2.86131E-03, &
   & 3.14118E-03, 3.14119E-03, 2.88881E-03, 3.19502E-03, 2.99538E-03, &
   & 2.91212E-03, 3.22739E-03, 3.05960E-03, 3.18901E-03, 3.05805E-03, &
   & 3.12205E-03, 2.95636E-03, 3.24111E-03, 3.29433E-03, 3.09206E-03, &
   & 3.06696E-03, 2.97735E-03, 2.90897E-03, 2.88979E-03, 2.75105E-03, &
   & 2.92156E-03, 3.03445E-03, 2.91664E-03, 2.85559E-03, 2.98405E-03, &
   & 2.95376E-03, 2.80234E-03, 2.78349E-03, 2.73421E-03, 2.70035E-03, &
   & 2.60074E-03, 2.34840E-03, 2.37626E-03, 2.32927E-03, 2.20842E-03, &
   & 2.31080E-03, 2.42771E-03, 2.43339E-03, 2.53280E-03, 2.37093E-03, &
   & 2.37377E-03, 2.73453E-03, 2.60836E-03, 2.55568E-03, 2.44062E-03, &
   & 2.71093E-03, 2.64421E-03, 2.66969E-03, 2.55560E-03, 2.71800E-03, &
   & 2.79534E-03, 2.59070E-03, 2.55373E-03, 2.45272E-03, 2.55571E-03/
   DATA C10561 /                                                     &
   & 2.54606E-03, 2.57349E-03, 2.46807E-03, 2.35634E-03, 2.44470E-03, &
   & 2.47050E-03, 2.57131E-03, 2.71649E-03, 2.58800E-03, 2.54524E-03, &
   & 2.69505E-03, 2.89122E-03, 2.77399E-03, 2.63306E-03, 2.82269E-03, &
   & 2.95684E-03, 3.07415E-03, 2.70594E-03, 2.65650E-03, 2.90613E-03, &
   & 2.96666E-03, 2.94767E-03, 2.81765E-03, 2.64829E-03, 2.43062E-03, &
   & 2.33816E-03, 2.38210E-03, 2.45701E-03, 2.38508E-03, 2.40746E-03, &
   & 2.49779E-03, 2.28209E-03, 2.26185E-03, 2.26604E-03, 2.19232E-03, &
   & 2.19160E-03, 2.32246E-03, 2.11108E-03, 2.26220E-03, 2.26849E-03, &
   & 2.34787E-03, 2.49323E-03, 2.46872E-03, 2.52974E-03, 2.35858E-03, &
   & 2.36865E-03, 2.33533E-03, 2.21338E-03, 2.24610E-03, 2.24776E-03, &
   & 2.24423E-03, 2.29276E-03, 2.18487E-03, 2.27621E-03, 2.31141E-03, &
   & 2.44095E-03, 2.45198E-03, 2.56919E-03, 2.56823E-03, 2.41982E-03, &
   & 2.39968E-03, 2.62447E-03, 2.55339E-03, 2.51556E-03, 2.47477E-03, &
   & 2.50276E-03, 2.48381E-03, 2.48484E-03, 2.48316E-03, 2.38541E-03, &
   & 2.41183E-03, 2.55888E-03, 2.42810E-03, 2.43356E-03, 2.25996E-03, &
   & 2.34736E-03, 2.10305E-03, 2.13870E-03, 2.17472E-03, 2.05354E-03/
   DATA C10641 /                                                     &
   & 2.11572E-03, 2.19557E-03, 2.09545E-03, 2.07831E-03, 1.94425E-03, &
   & 1.89333E-03, 1.98025E-03, 1.98328E-03, 2.01702E-03, 1.98333E-03, &
   & 2.01150E-03, 2.02484E-03, 2.10759E-03, 2.11892E-03, 2.10175E-03, &
   & 2.05314E-03, 2.13338E-03, 2.25764E-03, 2.19055E-03, 2.10818E-03, &
   & 2.05100E-03, 2.05685E-03, 2.10843E-03, 2.10228E-03, 2.10646E-03, &
   & 2.22640E-03, 2.31253E-03, 2.31230E-03, 2.21885E-03, 2.19568E-03, &
   & 2.23583E-03, 2.34754E-03, 2.28622E-03, 2.21876E-03, 2.26679E-03, &
   & 2.30828E-03, 2.24944E-03, 2.13851E-03, 2.02938E-03, 1.96770E-03, &
   & 2.05953E-03, 2.13814E-03, 2.03158E-03, 2.24655E-03, 1.95119E-03, &
   & 2.12979E-03, 2.08581E-03, 2.02434E-03, 1.98926E-03, 1.98792E-03, &
   & 1.97237E-03, 1.93397E-03, 1.92360E-03, 1.90805E-03, 1.89300E-03, &
   & 1.83548E-03, 1.87215E-03, 1.85589E-03, 1.85718E-03, 1.79361E-03, &
   & 1.77984E-03, 1.91506E-03, 2.04256E-03, 2.04095E-03, 1.94031E-03, &
   & 1.90447E-03, 2.02049E-03, 1.98360E-03, 2.04364E-03, 2.02519E-03, &
   & 2.20802E-03, 1.96964E-03, 1.94559E-03, 2.09922E-03, 2.11184E-03, &
   & 2.05706E-03, 2.02257E-03, 2.01781E-03, 2.01055E-03, 1.86538E-03/
   DATA C10721 /                                                     &
   & 1.86899E-03, 1.76798E-03, 1.85871E-03, 1.95363E-03, 1.96404E-03, &
   & 1.84169E-03, 1.82851E-03, 1.84582E-03, 1.81997E-03, 1.76461E-03, &
   & 1.68384E-03, 1.65530E-03, 1.73550E-03, 1.62463E-03, 1.68793E-03, &
   & 1.60472E-03, 1.67560E-03, 1.67431E-03, 1.61779E-03, 1.66446E-03, &
   & 1.66403E-03, 1.55724E-03, 1.62351E-03, 1.71545E-03, 1.69645E-03, &
   & 1.59540E-03, 1.62948E-03, 1.66784E-03, 1.66416E-03, 1.66131E-03, &
   & 1.71502E-03, 1.76555E-03, 1.75182E-03, 1.72327E-03, 1.72338E-03, &
   & 1.69993E-03, 1.78819E-03, 1.73517E-03, 1.74802E-03, 1.81751E-03, &
   & 1.70973E-03, 1.65075E-03, 1.70784E-03, 1.73655E-03, 1.71670E-03, &
   & 1.67367E-03, 1.69338E-03, 1.61772E-03, 1.54914E-03, 1.56009E-03, &
   & 1.59467E-03, 1.60761E-03, 1.57117E-03, 1.54045E-03, 1.53102E-03, &
   & 1.44516E-03, 1.49898E-03, 1.56048E-03, 1.60087E-03, 1.62636E-03, &
   & 1.62472E-03, 1.53931E-03, 1.55536E-03, 1.61649E-03, 1.66493E-03, &
   & 1.86915E-03, 1.59984E-03, 1.60483E-03, 1.66549E-03, 1.73449E-03, &
   & 1.73673E-03, 1.68393E-03, 1.67434E-03, 1.77880E-03, 1.76154E-03, &
   & 1.43028E-03, 1.69651E-03, 1.60934E-03, 1.69413E-03, 1.70514E-03/
   DATA C10801 /                                                     &
   & 1.62471E-03, 1.74854E-03, 1.76480E-03, 1.63495E-03, 1.59364E-03, &
   & 1.39603E-03, 1.47897E-03, 1.49509E-03, 1.70002E-03, 1.63048E-03, &
   & 1.44807E-03, 1.45071E-03, 1.53998E-03, 1.45276E-03, 1.29129E-03, &
   & 1.52900E-03, 1.64444E-03, 1.37450E-03, 1.42574E-03, 1.47355E-03, &
   & 1.51202E-03, 1.54376E-03, 1.51421E-03, 1.43989E-03, 1.45732E-03, &
   & 1.42912E-03, 1.59906E-03, 1.56748E-03, 1.52383E-03, 1.47665E-03, &
   & 1.51465E-03, 1.55582E-03, 1.54521E-03, 1.55189E-03, 1.56772E-03, &
   & 1.45401E-03, 1.55775E-03, 1.43120E-03, 1.39659E-03, 1.41451E-03, &
   & 1.45157E-03, 1.48303E-03, 1.42540E-03, 1.26387E-03, 1.37479E-03, &
   & 1.46381E-03, 1.38134E-03, 1.32733E-03, 1.38030E-03, 1.44619E-03, &
   & 1.41344E-03, 1.31982E-03, 1.24944E-03, 1.20096E-03, 1.21107E-03, &
   & 1.27999E-03, 1.22523E-03, 1.22193E-03, 1.35957E-03, 1.41427E-03, &
   & 1.35679E-03, 1.15438E-03, 1.41184E-03, 1.49093E-03, 1.32193E-03, &
   & 1.25009E-03, 1.37625E-03, 1.49022E-03, 1.44180E-03, 1.27628E-03, &
   & 1.29670E-03, 1.31636E-03, 1.28874E-03, 1.31177E-03, 1.35732E-03, &
   & 1.33854E-03, 1.30253E-03, 1.31374E-03, 1.27379E-03, 1.18339E-03/
   DATA C10881 /                                                     &
   & 1.22016E-03, 1.26551E-03, 1.26371E-03, 1.28180E-03, 1.36024E-03, &
   & 1.45759E-03, 1.29413E-03, 1.35858E-03, 1.26528E-03, 1.18623E-03, &
   & 1.21812E-03, 1.28799E-03, 1.37028E-03, 1.29268E-03, 1.27639E-03, &
   & 1.19487E-03, 1.23542E-03, 1.25010E-03, 1.17418E-03, 1.13914E-03, &
   & 1.21951E-03, 1.13780E-03, 1.16443E-03, 1.17883E-03, 1.11982E-03, &
   & 1.05708E-03, 1.04865E-03, 1.05884E-03, 1.06599E-03, 1.13828E-03, &
   & 1.10373E-03, 1.07739E-03, 1.04632E-03, 1.06118E-03, 1.15445E-03, &
   & 1.17300E-03, 1.00675E-03, 1.04235E-03, 1.08398E-03, 1.06587E-03, &
   & 1.05536E-03, 1.08614E-03, 1.09026E-03, 1.09141E-03, 1.13051E-03, &
   & 1.08667E-03, 1.04016E-03, 1.04897E-03, 1.08894E-03, 1.09682E-03, &
   & 1.09638E-03, 9.79254E-04, 1.00668E-03, 1.02569E-03, 1.00581E-03, &
   & 9.74433E-04, 9.66321E-04, 9.78440E-04, 9.01587E-04, 1.02149E-03, &
   & 9.87464E-04, 9.41872E-04, 9.05021E-04, 8.59547E-04, 9.03963E-04, &
   & 8.66415E-04, 8.84726E-04, 8.77087E-04, 8.70584E-04, 8.81338E-04, &
   & 8.97658E-04, 8.97586E-04, 9.19028E-04, 8.82438E-04, 9.00710E-04, &
   & 9.54329E-04, 9.54490E-04, 9.10940E-04, 9.95472E-04, 9.50134E-04/
   DATA C10961 /                                                     &
   & 9.17127E-04, 9.70916E-04, 9.87575E-04, 9.65026E-04, 9.71779E-04, &
   & 1.00967E-03, 1.00053E-03, 9.26063E-04, 9.34721E-04, 9.76354E-04, &
   & 9.78436E-04, 9.36012E-04, 9.64448E-04, 9.95903E-04, 9.89960E-04, &
   & 9.41143E-04, 9.04393E-04, 8.84719E-04, 8.41396E-04, 8.67234E-04, &
   & 8.55864E-04, 8.63314E-04, 8.72317E-04, 8.40899E-04, 7.79593E-04, &
   & 7.88481E-04, 8.21075E-04, 7.38342E-04, 7.56537E-04, 7.57278E-04, &
   & 7.35854E-04, 7.32765E-04, 6.67398E-04, 7.45338E-04, 7.33094E-04, &
   & 7.01840E-04, 6.85595E-04, 6.95740E-04, 7.24015E-04, 7.00907E-04, &
   & 7.28498E-04, 6.89410E-04, 6.91728E-04, 7.40601E-04, 7.62775E-04, &
   & 7.40912E-04, 7.35021E-04, 7.07799E-04, 7.54113E-04, 8.44845E-04, &
   & 8.53956E-04, 6.42186E-04, 7.40557E-04, 7.54340E-04, 7.55544E-04, &
   & 7.88986E-04, 7.97902E-04, 6.98460E-04, 7.74873E-04, 6.81178E-04, &
   & 7.15567E-04, 7.56723E-04, 7.98438E-04, 8.83150E-04, 8.45671E-04, &
   & 7.40924E-04, 7.35498E-04, 7.77829E-04, 6.93566E-04, 5.10188E-04, &
   & 7.52717E-04, 6.94185E-04, 6.71928E-04, 6.73286E-04, 6.89415E-04, &
   & 7.22917E-04, 7.89448E-04, 8.53812E-04, 7.45132E-04, 7.68732E-04/
   DATA C11041 /                                                     &
   & 8.10104E-04, 7.55615E-04, 7.09145E-04, 6.80676E-04, 7.54594E-04, &
   & 7.89416E-04, 7.88579E-04, 7.49805E-04, 6.13534E-04, 7.22491E-04, &
   & 7.95410E-04, 7.80604E-04, 7.74283E-04, 7.93224E-04, 6.86522E-04, &
   & 8.06038E-04, 8.30285E-04, 8.37763E-04, 8.03863E-04, 7.33526E-04, &
   & 7.42588E-04, 6.31046E-04, 8.16153E-04, 8.95391E-04, 8.61330E-04, &
   & 8.38726E-04, 8.16761E-04, 8.16118E-04, 6.37058E-04, 6.30868E-04, &
   & 7.26410E-04, 7.03464E-04, 5.93454E-04, 6.01985E-04, 6.51157E-04, &
   & 6.68569E-04, 6.56297E-04, 6.58732E-04, 5.99721E-04, 5.34301E-04, &
   & 5.33271E-04, 5.57992E-04, 5.70096E-04, 5.59932E-04, 5.32110E-04, &
   & 5.64713E-04, 6.25026E-04, 6.38973E-04, 6.05323E-04, 7.17460E-04, &
   & 6.19407E-04, 5.90228E-04, 5.43682E-04, 5.38446E-04, 6.56146E-04, &
   & 6.09081E-04, 6.04737E-04, 6.45526E-04, 6.46978E-04, 5.89738E-04, &
   & 5.63852E-04, 6.18018E-04, 5.71768E-04, 5.75433E-04, 6.05766E-04, &
   & 5.93065E-04, 5.31708E-04, 5.41187E-04, 5.76985E-04, 5.78176E-04, &
   & 5.75339E-04, 6.85426E-04, 5.51038E-04, 6.02049E-04, 6.20406E-04, &
   & 5.80169E-04, 5.36399E-04, 5.59608E-04, 5.46575E-04, 5.66979E-04/
   DATA C11121 /                                                     &
   & 5.94982E-04, 6.18469E-04, 6.56281E-04, 8.22124E-04, 7.81716E-04, &
   & 7.29616E-04, 7.14460E-04, 7.08969E-04, 6.53794E-04, 7.33138E-04, &
   & 8.29513E-04, 8.99395E-04, 9.05526E-04, 7.98257E-04, 7.86935E-04, &
   & 6.10797E-04, 4.63912E-04, 4.05675E-04, 3.66230E-04, 4.86472E-04, &
   & 5.31818E-04, 5.15865E-04, 4.87344E-04, 4.99857E-04, 5.35479E-04, &
   & 5.27561E-04, 4.99000E-04, 4.77056E-04, 4.74242E-04, 4.66595E-04, &
   & 4.66325E-04, 4.94704E-04, 5.12842E-04, 5.01795E-04, 4.80789E-04, &
   & 5.73709E-04, 5.65214E-04, 5.11321E-04, 4.55242E-04, 4.29330E-04, &
   & 5.09792E-04, 4.70489E-04, 4.82859E-04, 4.99195E-04, 4.07724E-04, &
   & 4.99951E-04, 4.55755E-04, 4.42528E-04, 4.19433E-04, 3.31325E-04, &
   & 3.70517E-04, 3.77708E-04, 2.97923E-04, 2.27470E-04, 2.47389E-04, &
   & 2.38324E-04, 2.56706E-04, 2.45046E-04, 2.62539E-04, 3.37054E-04, &
   & 3.33930E-04, 3.01390E-04, 3.08028E-04, 3.41464E-04, 3.70574E-04, &
   & 3.47893E-04, 3.28433E-04, 3.46976E-04, 3.60351E-04, 3.50559E-04, &
   & 3.56070E-04, 3.62782E-04, 3.37330E-04, 3.33763E-04, 3.57046E-04, &
   & 3.08784E-04, 2.93898E-04, 2.80842E-04, 2.54114E-04, 2.38198E-04/
   DATA C11201 /                                                     &
   & 3.48753E-04, 2.97334E-04, 2.82929E-04, 2.94150E-04, 3.07875E-04, &
   & 3.21129E-04, 3.38335E-04, 3.49826E-04, 3.47647E-04, 3.35438E-04, &
   & 3.58145E-04, 3.72391E-04, 3.59372E-04, 3.64755E-04, 4.16867E-04, &
   & 3.43614E-04, 3.34932E-04, 3.12782E-04, 3.28220E-04, 4.32595E-04, &
   & 3.49513E-04, 3.51861E-04, 3.81166E-04, 3.91194E-04, 3.38944E-04, &
   & 2.63445E-04, 2.49520E-04, 2.46184E-04, 2.33203E-04, 2.16315E-04, &
   & 1.89536E-04, 1.95730E-04, 1.99664E-04, 1.77139E-04, 1.27969E-04, &
   & 5.17216E-05, 7.60445E-05, 1.24418E-04, 1.30989E-04, 2.31539E-04, &
   & 2.21334E-04, 2.08757E-04, 2.18351E-04, 2.46202E-04, 2.29824E-04, &
   & 2.28909E-04, 2.88826E-04, 3.58039E-04, 2.60800E-04, 2.33025E-04, &
   & 2.52667E-04, 2.61394E-04, 2.31384E-04, 2.29388E-04, 2.54701E-04, &
   & 2.21158E-04, 1.61506E-04, 1.36752E-04, 1.69481E-04, 8.64539E-05, &
   & 1.64407E-04, 3.65674E-04, 3.18233E-04, 4.00755E-04, 3.33375E-04, &
   & 2.62930E-04, 2.87052E-04, 2.51395E-04, 2.85274E-04, 2.66915E-04, &
   & 2.10866E-04, 1.89517E-04, 1.67378E-04, 2.79951E-04, 2.97224E-04, &
   & 1.89222E-04, 3.33825E-04, 3.56386E-04, 3.89727E-04, 4.30407E-04/
   DATA C11281 /                                                     &
   & 4.45922E-04, 4.23446E-04, 4.41347E-04, 4.06723E-04, 3.00181E-04, &
   & 1.85243E-04, 3.13176E-04, 4.08991E-04, 4.24776E-04, 3.56412E-04, &
   & 3.84760E-04, 2.30602E-04, 1.77702E-04, 2.62329E-04, 2.49442E-04, &
   & 3.76212E-04, 3.69176E-04, 2.97681E-04, 2.71662E-04, 2.05694E-04, &
   & 2.11418E-04, 2.25439E-04, 2.27013E-04, 2.47845E-04, 3.14603E-04, &
   & 2.68802E-04, 2.04334E-04, 2.77399E-04, 2.68273E-04, 2.04991E-04, &
   & 2.24441E-04, 3.55074E-04, 2.90135E-04, 3.35680E-04, 3.59358E-04, &
   & 3.44716E-04, 3.24496E-04, 3.48146E-04, 3.49042E-04, 3.54848E-04, &
   & 3.86418E-04, 3.59198E-04, 3.47608E-04, 3.20522E-04, 2.78401E-04, &
   & 2.64579E-04, 2.23694E-04, 2.34370E-04, 2.52559E-04, 1.88475E-04, &
   & 2.01258E-04, 1.63979E-04, 1.45384E-04, 1.91215E-04, 1.76958E-04, &
   & 1.69167E-04, 1.71767E-04, 1.86595E-04, 2.14969E-04, 2.48345E-04, &
   & 2.46691E-04, 2.25234E-04, 2.26755E-04, 1.64112E-04, 1.87750E-04, &
   & 2.22984E-04, 2.00443E-04, 2.38863E-04, 2.77590E-04, 2.91953E-04, &
   & 2.80611E-04, 3.08215E-04, 1.79095E-04, 1.46920E-04, 2.29177E-04, &
   & 2.54685E-04, 2.68866E-04, 2.13346E-04, 1.20122E-04, 5.55240E-05/
   DATA C11361 /                                                     &
   & 5.99017E-05, 1.07768E-04, 1.67810E-04, 2.06886E-04, 2.36232E-04, &
   & 2.24598E-04, 2.30792E-04, 2.71274E-04, 1.29062E-04, 1.92624E-04, &
   & 2.38438E-04, 1.98994E-04, 1.81687E-04, 2.55733E-04, 2.84379E-04, &
   & 2.54459E-04, 2.30884E-04, 2.68873E-04, 3.07231E-04, 3.15063E-04, &
   & 2.46725E-04, 2.60370E-04, 2.66391E-04, 2.50708E-04, 2.04296E-04, &
   & 1.66011E-04, 1.19164E-04, 1.06700E-04, 1.77576E-04, 1.91741E-04, &
   & 1.66618E-04, 1.49824E-04, 1.80699E-04, 2.20905E-04, 1.38754E-04, &
   & 6.27971E-05, 7.52567E-05, 1.89995E-04, 1.72489E-04, 1.40424E-04, &
   & 1.52384E-04, 1.63942E-04, 1.19901E-04, 1.49234E-04, 2.68313E-04, &
   & 2.08815E-04, 1.17218E-04, 1.42235E-04, 2.71237E-04, 1.38192E-04, &
   & 2.15643E-04, 2.84476E-04, 2.78117E-04, 2.19234E-04, 1.59128E-04, &
   & 1.78819E-04, 2.67785E-04, 2.66786E-04, 2.58545E-04, 2.68476E-04, &
   & 2.88542E-04, 2.59726E-04, 3.00936E-04, 3.11237E-04, 2.61275E-04, &
   & 1.37136E-04, 2.76566E-04, 3.82888E-04, 3.97564E-04, 4.43655E-04, &
   & 3.15415E-04, 2.60869E-04, 3.19171E-04, 3.34205E-04, 2.02914E-04, &
   & 1.16223E-04, 1.14737E-04, 6.10978E-05,-8.03695E-06,-1.07062E-05/
   DATA C11441 /                                                     &
   & 6.50664E-05, 1.12586E-04, 1.56727E-04, 1.57927E-04, 1.05762E-04, &
   & 1.03646E-04, 1.72520E-04, 2.23668E-04, 2.12775E-04, 2.33525E-04, &
   & 2.75558E-04, 2.34256E-04, 5.10062E-05, 1.76007E-04, 1.70850E-04, &
   & 1.43266E-04, 1.89626E-04, 2.97283E-04, 3.02773E-04, 2.74401E-04, &
   & 3.00754E-04, 3.66813E-04, 3.54383E-04, 2.90580E-04, 2.32206E-04, &
   & 1.58405E-04, 1.54663E-04, 1.84598E-04, 1.26408E-04, 2.14481E-04, &
   & 2.00791E-04, 1.05796E-04, 2.39794E-04, 1.66105E-04, 7.88615E-05, &
   & 4.30615E-05, 7.37518E-05, 1.24926E-04, 1.38295E-04, 8.54356E-05, &
   & 6.12641E-05, 6.54466E-05, 6.17727E-05, 1.30688E-05, 6.00462E-05, &
   & 1.52612E-04, 2.11656E-04, 9.67692E-05, 8.67858E-05, 1.34888E-04, &
   & 1.90899E-04, 1.03234E-04, 1.03837E-04, 1.49767E-04, 2.19058E-04, &
   & 2.26549E-04, 2.11506E-04, 1.85238E-04, 1.53774E-04, 1.32313E-04, &
   & 6.10658E-05, 2.37782E-05, 1.24450E-04, 1.87610E-04, 1.44775E-04, &
   & 5.60937E-05, 6.64032E-05, 1.28073E-04, 1.77512E-04, 1.84684E-04, &
   & 5.73677E-05, 5.29679E-05, 9.95510E-05, 1.61423E-04, 3.19036E-04, &
   & 3.17383E-04, 2.36505E-04, 1.80844E-04, 1.63722E-04, 1.21478E-04/
   DATA C11521 /                                                     &
   & 6.85823E-05, 7.42058E-05, 1.14838E-04, 1.21131E-04, 8.01009E-05, &
   & 1.52058E-04, 2.18368E-04, 2.53416E-04, 2.27116E-04, 1.25336E-04, &
   & 6.26421E-05, 5.32471E-05, 1.34705E-04, 2.07005E-05,-5.18630E-05, &
   &-3.25696E-05,-8.06171E-05,-1.09430E-04,-1.05637E-04,-4.96066E-05, &
   &-7.76138E-05,-4.85930E-05, 3.65111E-06,-2.86933E-05,-4.61366E-05, &
   &-4.88820E-05,-3.08816E-05, 8.43778E-05, 1.40484E-04, 1.31125E-04, &
   & 3.55198E-05, 8.47412E-05, 1.23408E-04, 1.36799E-04, 1.21147E-04, &
   & 1.25585E-04, 1.32337E-04, 1.34092E-04, 1.26652E-04, 1.12131E-04, &
   & 1.00927E-04, 1.13828E-04, 1.06053E-04, 9.43643E-05, 8.33628E-05, &
   & 8.65842E-05, 7.59315E-05, 8.28623E-05, 1.39681E-04, 1.80492E-04, &
   & 1.65779E-04, 1.03843E-04, 3.10284E-05, 1.94408E-05, 4.57525E-05, &
   & 1.02436E-04, 1.39750E-04, 1.43342E-04, 1.11999E-04, 2.94197E-05, &
   & 2.76980E-05, 5.51685E-05, 9.39909E-05, 1.16108E-04, 7.72703E-05, &
   & 4.37409E-05, 1.13925E-04, 8.18872E-05, 2.87657E-05,-2.41413E-05, &
   & 1.24699E-05, 2.19589E-05,-5.88247E-06,-9.66151E-05,-2.06255E-05, &
   &-1.83148E-06,-5.63625E-05,-8.65590E-05,-8.26020E-05,-5.06239E-05/
   DATA C11601 /                                                     &
   & 1.28065E-05,-1.34669E-05, 1.59701E-05, 9.44755E-05, 1.63032E-05, &
   & 2.51304E-05, 7.38226E-05, 1.28405E-04, 1.17413E-04, 9.92387E-05, &
   & 9.51533E-05, 2.17008E-04, 2.25854E-04, 1.90448E-04, 1.77207E-04, &
   & 1.80844E-04, 1.53501E-04, 9.80430E-05, 1.27404E-04, 1.16465E-04, &
   & 9.98611E-05, 1.25556E-04, 1.73627E-04, 1.12347E-04,-7.73523E-05, &
   & 5.66599E-05, 5.36347E-05, 1.20227E-06, 6.96325E-05, 4.79010E-05, &
   &-1.09886E-05,-9.16457E-05,-7.09170E-05,-5.31410E-05,-2.68376E-05, &
   & 6.32641E-05, 8.06052E-06,-4.99262E-05,-2.56644E-05,-8.76854E-05, &
   &-8.21360E-05,-5.02403E-06, 4.66629E-05, 6.93127E-05, 5.53828E-05, &
   &-2.32399E-05,-2.07514E-05,-7.33240E-05,-2.10483E-04,-1.53757E-04, &
   &-7.13861E-05,-1.07356E-05,-1.26578E-04,-7.48854E-05, 3.25418E-06, &
   & 2.97068E-05, 3.35685E-05, 3.15022E-05, 2.68904E-05, 3.87401E-05, &
   & 5.12522E-05, 5.12172E-05, 1.05053E-05, 1.65321E-05, 3.47537E-05, &
   & 5.62503E-05, 4.18666E-05, 3.13970E-05, 3.11750E-05, 7.21547E-05, &
   & 2.55262E-05,-2.76061E-05, 5.43449E-06,-5.20575E-05,-1.08627E-04, &
   &-1.40475E-04,-1.59926E-04,-1.32237E-04,-8.15458E-05,-1.31738E-04/
   DATA C11681 /                                                     &
   &-1.64036E-04,-1.69351E-04,-1.24797E-04,-1.61950E-04,-2.01904E-04, &
   &-2.22995E-04,-1.87647E-04,-1.70817E-04,-1.64583E-04,-1.12811E-04, &
   &-8.38306E-05,-8.62707E-05,-1.54362E-04,-1.98090E-04,-2.12920E-04, &
   &-1.89358E-04,-2.02988E-04,-1.72791E-04,-1.02863E-04,-1.09877E-04, &
   &-1.04257E-04,-8.20734E-05,-2.18346E-05,-2.94593E-05,-4.18226E-05, &
   &-1.86891E-05,-6.14620E-05,-3.21912E-05, 1.00844E-04, 6.92419E-05, &
   & 3.16713E-05, 5.62042E-07, 5.18900E-05, 7.48835E-05, 8.03381E-05, &
   & 7.24685E-05, 9.55588E-05, 9.22801E-05, 2.87159E-05, 2.26234E-05, &
   & 2.62790E-05, 3.58332E-05, 6.23297E-05, 5.01998E-05, 1.81446E-05, &
   & 3.33564E-05, 3.97765E-06,-2.60624E-05, 7.01802E-06,-4.16797E-05, &
   &-8.70108E-05,-8.22182E-05,-6.64886E-05,-7.88704E-05,-1.28305E-04, &
   &-1.29990E-04,-1.12646E-04,-8.68394E-05,-1.29584E-04,-1.44352E-04, &
   &-1.42082E-04,-1.33790E-04,-1.27963E-04,-1.21233E-04,-1.09965E-04, &
   &-1.02233E-04,-1.03804E-04,-1.19503E-04,-7.74707E-05,-4.66805E-05, &
   &-3.52201E-05,-4.07406E-05,-4.66887E-05,-5.05962E-05,-3.30333E-05, &
   &-3.47981E-05,-3.60962E-05, 1.44242E-05, 4.10478E-05, 3.68984E-05/
   DATA C11761 /                                                     &
   &-2.81300E-05, 2.83171E-05, 7.48062E-05, 4.29333E-05, 8.50076E-06, &
   & 4.98135E-06, 4.44854E-05, 2.51860E-05, 3.12189E-05, 6.39424E-05, &
   & 7.20715E-05, 9.89688E-05, 1.33768E-04, 1.07781E-04, 9.76731E-05, &
   & 9.21479E-05, 6.72624E-05, 5.41295E-05, 4.89022E-05, 5.28039E-05, &
   &-4.48737E-06,-5.15409E-05,-3.57396E-05,-1.94752E-05,-2.09521E-05, &
   &-5.13096E-05,-2.62781E-05,-2.75451E-05,-6.98423E-05,-1.25462E-04, &
   &-1.68362E-04,-1.97456E-04,-1.90669E-04,-2.06890E-04,-2.36699E-04, &
   &-1.97732E-04,-1.76504E-04,-1.67505E-04,-1.60694E-04,-1.85851E-04, &
   &-2.01567E-04,-9.82507E-05,-1.33338E-04,-1.95199E-04,-1.40781E-04, &
   &-8.90988E-05,-3.63239E-05, 2.16510E-05,-1.56807E-05,-4.21285E-05, &
   & 5.50505E-06, 6.78937E-07, 3.12346E-06, 3.64202E-05, 3.50651E-05, &
   & 6.20423E-05, 1.38667E-04, 7.74738E-05, 6.77036E-05, 1.38367E-04, &
   & 1.17359E-04, 1.06637E-04, 1.12404E-04, 9.78586E-05, 1.03178E-04, &
   & 1.28717E-04, 1.56642E-04, 1.62544E-04, 1.50109E-04, 1.43214E-04, &
   & 1.33651E-04, 1.24352E-04, 1.41420E-04, 1.36340E-04, 1.18769E-04, &
   & 1.31656E-04, 8.81533E-05, 1.55214E-05,-3.68736E-07,-1.76213E-05/
   DATA C11841 /                                                     &
   &-2.85341E-05, 4.65155E-06, 5.41350E-06,-7.01247E-06, 6.57918E-06, &
   &-2.45784E-05,-6.89104E-05,-6.90953E-05,-6.23937E-05,-6.72978E-05, &
   &-1.39547E-04,-1.44228E-04,-1.42543E-04,-2.31080E-04,-2.12756E-04, &
   &-1.62089E-04,-1.66063E-04,-1.61872E-04,-1.59764E-04,-1.80217E-04, &
   &-1.38355E-04,-8.45661E-05,-7.58308E-05,-4.65144E-05,-2.76855E-05, &
   &-7.48714E-05,-8.28561E-05,-6.45277E-05,-7.08509E-06,-1.05566E-05, &
   &-1.96352E-05, 3.55561E-05, 2.24676E-05,-1.25648E-05,-1.87661E-05, &
   & 6.99061E-06, 2.33676E-05,-5.25111E-05,-3.86758E-05, 1.03585E-06, &
   &-1.65901E-05,-1.04855E-05, 5.03694E-06, 1.25937E-05,-8.31340E-06, &
   &-4.37906E-05,-7.91444E-05,-4.62167E-05, 5.14238E-06,-4.52863E-05, &
   &-5.86455E-05,-4.98093E-05,-3.03495E-05,-5.09377E-05,-8.88116E-05, &
   &-6.21360E-05,-7.38148E-05,-1.07502E-04,-7.55276E-05,-6.39257E-05, &
   &-6.86921E-05,-8.05504E-05,-9.24178E-05,-1.03991E-04,-1.00468E-04, &
   &-6.71447E-05,-3.84897E-06,-5.99067E-06,-2.21894E-05,-5.21766E-05, &
   &-3.93796E-05,-4.06712E-05,-6.21649E-05,-1.13073E-04,-1.20560E-04, &
   &-5.92397E-05, 5.24432E-05, 9.41628E-05,-3.47458E-07, 5.33267E-05/
   DATA C11921 /                                                     &
   & 8.92961E-05, 2.75694E-05,-7.48460E-06,-2.15504E-05, 1.05501E-06, &
   & 6.30910E-06, 5.94620E-07,-2.45194E-05,-1.59657E-05, 7.93610E-07, &
   &-1.05319E-05,-2.36584E-05,-3.95700E-05,-6.57225E-05,-5.23797E-05, &
   &-1.82588E-05,-1.43240E-05,-3.29989E-05,-6.48909E-05,-2.41326E-05, &
   &-1.89195E-05,-4.64607E-05,-1.00739E-05,-1.35033E-05,-6.49945E-05, &
   &-5.19986E-05,-6.68505E-05,-1.31530E-04,-1.45464E-04,-1.46815E-04, &
   &-1.39684E-04,-1.23205E-04,-1.26738E-04,-1.93822E-04,-2.37508E-04, &
   &-2.52917E-04,-1.91110E-04,-1.36217E-04,-9.41093E-05,-1.20601E-04, &
   &-1.17295E-04,-9.57420E-05,-1.57227E-04,-1.62795E-04,-1.12201E-04, &
   &-1.20419E-04,-1.10597E-04,-7.61223E-05,-6.27167E-05,-5.54733E-05, &
   &-5.50437E-05,-5.14148E-05,-3.59591E-05, 1.09906E-05, 5.94396E-06, &
   &-1.38597E-05,-8.80857E-06,-3.13101E-05,-6.31715E-05,-4.04264E-05, &
   &-1.66405E-05, 7.94396E-06,-3.41772E-06,-4.03175E-05,-1.06888E-04, &
   &-9.50526E-05,-7.46111E-05,-5.09617E-05,-6.70981E-05,-7.93529E-05, &
   &-5.58423E-05,-1.01523E-04,-1.62269E-04,-1.69958E-04,-1.37786E-04, &
   &-8.79862E-05,-1.46838E-04,-1.66938E-04,-1.51380E-04,-1.62184E-04/
   DATA C12001 /                                                     &
   &-1.61105E-04,-1.42088E-04,-1.57033E-04,-1.65294E-04,-1.45079E-04, &
   &-9.76982E-05,-6.09891E-05,-1.01719E-04,-1.03049E-04,-8.85546E-05, &
   &-1.47754E-04,-1.44542E-04,-8.34620E-05,-8.99440E-05,-7.11901E-05, &
   &-1.57480E-05,-8.81797E-05,-1.56314E-04,-1.65952E-04,-1.80986E-04, &
   &-2.04610E-04,-2.58669E-04,-2.16016E-04,-1.21582E-04,-1.44929E-04, &
   &-1.72886E-04,-2.05950E-04,-1.93829E-04,-1.67518E-04,-1.22969E-04, &
   &-1.13060E-04,-1.14854E-04,-1.26198E-04,-1.24288E-04,-1.19519E-04, &
   &-1.50456E-04,-1.53286E-04,-1.32231E-04,-7.42672E-05,-2.23129E-05, &
   & 1.79115E-05, 1.42073E-05,-1.21676E-05,-7.56567E-05,-1.03423E-04, &
   &-1.10373E-04,-8.77244E-05,-6.43485E-05,-4.05156E-05,-6.24405E-05, &
   &-5.70375E-05,-2.36695E-06,-3.75929E-05,-7.97119E-05,-6.70419E-05, &
   &-6.99475E-05,-8.19748E-05,-1.06895E-04,-1.31422E-04,-1.55438E-04, &
   &-1.61937E-04,-1.62626E-04,-1.54977E-04,-1.77814E-04,-2.00386E-04, &
   &-1.87407E-04,-2.07243E-04,-2.44672E-04,-2.19014E-04,-2.13695E-04, &
   &-2.32440E-04,-1.85194E-04,-1.51172E-04,-1.69834E-04,-1.73780E-04, &
   &-1.75232E-04,-2.00698E-04,-1.82826E-04,-1.27786E-04,-1.33633E-04/
   DATA C12081 /                                                     &
   &-1.21317E-04,-7.50390E-05,-1.06743E-04,-1.40805E-04,-1.06336E-04, &
   &-9.46654E-05,-9.78182E-05,-1.19906E-04,-1.14160E-04,-7.28186E-05, &
   &-1.07652E-04,-1.20978E-04,-3.79658E-05,-3.16113E-05,-6.02417E-05, &
   &-7.51148E-05,-5.56145E-05,-6.77421E-06,-1.74321E-05,-4.67952E-05, &
   &-1.05000E-04,-6.29932E-05,-4.74356E-06,-2.83397E-05,-4.65192E-05, &
   &-6.04574E-05,-4.33970E-05,-3.18311E-05,-3.02321E-05,-4.49667E-05, &
   &-6.85347E-05,-1.11375E-04,-1.16293E-04,-9.38757E-05,-1.38594E-04, &
   &-1.60483E-04,-1.48344E-04,-1.33436E-04,-1.27387E-04,-1.59508E-04, &
   &-1.74026E-04,-1.72170E-04,-1.49196E-04,-1.33233E-04,-1.22382E-04, &
   &-1.78156E-04,-2.21349E-04,-2.41846E-04,-2.06549E-04,-1.68283E-04, &
   &-1.89512E-04,-1.44523E-04,-4.67953E-05,-1.00334E-04,-1.23478E-04, &
   &-8.14024E-05,-9.18016E-05,-1.17536E-04,-1.36160E-04,-1.38780E-04, &
   &-1.27749E-04,-1.45598E-04,-1.55964E-04,-1.45120E-04,-1.25544E-04, &
   &-1.05692E-04,-1.17639E-04,-1.24142E-04,-1.24749E-04,-1.63878E-04, &
   &-1.97021E-04,-1.98617E-04,-2.69136E-04,-3.68357E-04,-2.33702E-04, &
   &-1.61830E-04,-1.78578E-04,-2.01839E-04,-2.28731E-04,-2.63606E-04/
   DATA C12161 /                                                     &
   &-2.44698E-04,-1.86451E-04,-2.20546E-04,-2.22752E-04,-1.55169E-04, &
   &-1.25100E-04,-1.09794E-04,-9.59016E-05,-1.03857E-04,-1.35573E-04, &
   &-1.73780E-04,-1.82457E-04,-9.39821E-05,-1.18245E-04,-2.11563E-04, &
   &-1.37392E-04,-9.28173E-05,-9.71073E-05,-9.72535E-05,-9.39557E-05, &
   &-7.50117E-05,-6.70754E-05,-7.01186E-05,-5.76151E-05,-5.18785E-05, &
   &-7.14209E-05,-7.01682E-05,-5.61614E-05,-8.92769E-05,-1.06238E-04, &
   &-9.70294E-05,-6.70229E-05,-4.69214E-05,-1.53105E-04,-2.02326E-04, &
   &-1.90395E-04,-2.04367E-04,-2.16787E-04,-2.08725E-04,-1.78119E-04, &
   &-1.31043E-04,-1.32204E-04,-1.51522E-04,-2.05143E-04,-1.77144E-04, &
   &-1.16130E-04,-1.44440E-04,-1.66010E-04,-1.78206E-04,-1.61163E-04, &
   &-1.46351E-04,-1.96722E-04,-2.27027E-04,-2.37243E-04,-2.25235E-04, &
   &-1.99552E-04,-1.40238E-04,-1.26311E-04,-1.42746E-04,-1.19028E-04, &
   &-1.18750E-04,-1.72076E-04,-1.72120E-04,-1.48285E-04,-1.85116E-04, &
   &-1.98602E-04,-1.74016E-04,-1.37913E-04,-1.01221E-04,-9.69581E-05, &
   &-1.08794E-04,-1.39433E-04,-1.38575E-04,-1.32088E-04,-1.37431E-04, &
   &-1.30033E-04,-1.10829E-04,-1.35604E-04,-1.66515E-04,-1.98167E-04/
   DATA C12241 /                                                     &
   &-1.97716E-04,-1.74019E-04,-1.64719E-04,-1.64779E-04,-1.85725E-04, &
   &-2.28526E-04,-2.84329E-04,-1.82449E-04,-1.30747E-04,-1.93620E-04, &
   &-2.28529E-04,-2.47361E-04,-1.90001E-04,-1.66278E-04,-2.02540E-04, &
   &-2.31811E-04,-2.53772E-04,-2.08629E-04,-1.85021E-04,-1.93989E-04, &
   &-2.16568E-04,-2.38288E-04,-1.94453E-04,-1.87154E-04,-2.30493E-04, &
   &-2.34696E-04,-2.30351E-04,-2.60562E-04,-2.86427E-04,-3.06699E-04, &
   &-2.79131E-04,-2.49392E-04,-3.03389E-04,-3.10346E-04,-2.61782E-04, &
   &-2.30905E-04,-2.11669E-04,-2.37680E-04,-2.38194E-04,-2.10955E-04/
!
end block data BO3HH1
!
!     --------------------------------------------------------------
!
SUBROUTINE O3HHT2 (V1C,V2C,DVC,NPTC,C)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   COMMON /O3HH2/ V1S,V2S,DVS,NPTS,S(2687)
   DIMENSION C(*)
!
   DVC = DVS
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      C(J) = S(I)
10 END DO
!
   RETURN
!
end subroutine O3HHT2
!
!     --------------------------------------------------------------
!
BLOCK DATA BO3HH2
!
   IMPLICIT REAL*8           (V)
!
!     RATIO (C2/C0)
!     DATA FROM BASS 1985
!
!     NOW INCLUDES MOLINA & MOLINA AT 273K WITH THE TEMPERATURE
!     DEPENDENCE DETERMINED FROM THE 195K HARVARD MEASUREMENTS,
!     EMPLOYING THE BASS ALGORITHM
!
!              (CO(1+C1*(T-273.15)+C2*(T-273.15)**2);
!
!     THIS IS ONLY FOR THE WAVELENGTH RANGE FROM .34 TO .35 MICRONS;
!     OTHERWISE, THE BASS DATA ALONE HAVE BEEN EMPLOYED BETWEEN
!     .34 AND .245 MICRONS.
!
!     NEW T-DEPENDENT X-SECTIONS BETWEEN .345 AND .36 MICRONS
!     HAVE NOW BEEN ADDED, BASED ON WORK BY CACCIANI, DISARRA
!     AND FIOCCO, UNIVERSITY OF ROME, 1987.  QUADRATIC TEMP
!     HAS BEEN DERIVED, AS ABOVE.
!
!     AGREEMENT AMONGST THE FOUR DATA SETS IS REASONABLE (<10%)
!     AND OFTEN EXCELLENT (0-3%)
!
!
   COMMON /O3HH2/ V1C,V2C,DVC,NC,                                    &
   &               O32001(85),C20086(80),C20166(80),C20246(65),       &
   &               C20311(16),C20327(80),C20407( 1),                  &
   &               C20001(80),C20081(80),C20161(80),C20241(80),       &
   &               C20321(80),C20401(80),C20481(80),C20561(80),       &
   &               C20641(80),C20721(80),C20801(80),C20881(80),       &
   &               C20961(80),C21041(80),C21121(80),C21201(80),       &
   &               C21281(80),C21361(80),C21441(80),C21521(80),       &
   &               C21601(80),C21681(80),C21761(80),C21841(80),       &
   &               C21921(80),C22001(80),C22081(80),C22161(80),       &
   &               C22241(40)
!
!     DATA V1C /29405./, V2C /40800./ ,DVC /5./, NC /2280/   BASS
!
   DATA V1C /27370./, V2C /40800./ ,DVC /5./, NC /2687/
!
   DATA O32001/85*1.0E-5/
!
   DATA C20086/                                                      &
   & 1.29359E-05, 1.55806E-05, 2.00719E-05, 2.64912E-05, 3.48207E-05, &
   & 4.36986E-05, 5.31318E-05, 6.13173E-05, 6.89465E-05, 7.56793E-05, &
   & 8.26345E-05, 8.90916E-05, 9.38759E-05, 9.22998E-05, 9.03184E-05, &
   & 8.65369E-05, 8.58531E-05, 8.55635E-05, 8.40418E-05, 8.11983E-05, &
   & 7.58246E-05, 7.29282E-05, 7.32629E-05, 7.04060E-05, 6.71451E-05, &
   & 6.56515E-05, 6.68943E-05, 6.32785E-05, 5.88386E-05, 5.70860E-05, &
   & 5.64435E-05, 5.49441E-05, 5.70845E-05, 5.89357E-05, 6.14433E-05, &
   & 5.91790E-05, 5.31727E-05, 5.14007E-05, 4.74318E-05, 4.35356E-05, &
   & 3.93903E-05, 3.70963E-05, 3.63867E-05, 4.05296E-05, 4.48891E-05, &
   & 5.37190E-05, 5.70440E-05, 4.60408E-05, 5.25778E-05, 6.81728E-05, &
   & 7.27275E-05, 6.81353E-05, 6.48386E-05, 5.46521E-05, 4.93098E-05, &
   & 5.04438E-05, 5.30309E-05, 5.28788E-05, 5.47387E-05, 4.52523E-05, &
   & 5.29451E-05, 7.42215E-05, 1.08971E-04, 1.40085E-04, 1.46553E-04, &
   & 1.43526E-04, 1.39051E-04, 1.40983E-04, 1.45564E-04, 1.55589E-04, &
   & 1.66142E-04, 1.82840E-04, 2.06486E-04, 2.24339E-04, 2.29268E-04, &
   & 2.13109E-04, 2.00305E-04, 1.99955E-04, 2.18566E-04, 2.24182E-04/
   DATA C20166/                                                      &
   & 2.33505E-04, 2.31824E-04, 2.22666E-04, 2.23905E-04, 2.38131E-04, &
   & 2.54322E-04, 2.69548E-04, 2.62953E-04, 2.67609E-04, 2.70567E-04, &
   & 2.70689E-04, 2.68251E-04, 2.66029E-04, 2.60053E-04, 2.61689E-04, &
   & 2.56582E-04, 2.43655E-04, 2.38792E-04, 2.45309E-04, 2.31061E-04, &
   & 2.22837E-04, 2.16440E-04, 2.19032E-04, 1.85634E-04, 1.74638E-04, &
   & 1.51767E-04, 1.38480E-04, 1.32506E-04, 1.28317E-04, 1.26855E-04, &
   & 1.27123E-04, 1.24040E-04, 1.19202E-04, 1.28649E-04, 1.36271E-04, &
   & 1.42080E-04, 1.47804E-04, 1.39534E-04, 1.27284E-04, 1.09554E-04, &
   & 8.69470E-05, 6.72096E-05, 5.23407E-05, 5.12433E-05, 5.15794E-05, &
   & 4.94683E-05, 4.95809E-05, 4.07499E-05, 3.14984E-05, 1.46457E-05, &
   & 6.98660E-06, 1.85313E-05, 5.48879E-05, 1.09447E-04, 1.52536E-04, &
   & 1.78778E-04, 1.91128E-04, 1.99161E-04, 2.02937E-04, 1.95527E-04, &
   & 1.92214E-04, 1.83376E-04, 1.81710E-04, 1.82283E-04, 1.75182E-04, &
   & 1.72406E-04, 1.68170E-04, 1.67400E-04, 1.69469E-04, 1.69092E-04, &
   & 1.65985E-04, 1.66912E-04, 1.74226E-04, 1.85036E-04, 1.85517E-04, &
   & 1.85805E-04, 1.73809E-04, 1.67628E-04, 1.57690E-04, 1.54952E-04/
   DATA C20246/                                                      &
   & 1.53707E-04, 1.57710E-04, 1.58175E-04, 1.67253E-04, 1.82079E-04, &
   & 1.91285E-04, 1.96564E-04, 2.03822E-04, 1.93736E-04, 1.82924E-04, &
   & 1.73610E-04, 1.69904E-04, 1.66725E-04, 1.63747E-04, 1.63129E-04, &
   & 1.62435E-04, 1.67218E-04, 1.69507E-04, 1.70744E-04, 1.65839E-04, &
   & 1.72077E-04, 1.67734E-04, 1.51487E-04, 1.43770E-04, 1.37435E-04, &
   & 1.25172E-04, 1.12395E-04, 1.07991E-04, 1.00345E-04, 9.36611E-05, &
   & 9.59763E-05, 9.26600E-05, 1.00120E-04, 1.04746E-04, 1.10222E-04, &
   & 1.03308E-04, 8.97457E-05, 7.91634E-05, 7.50275E-05, 8.30832E-05, &
   & 1.01191E-04, 1.21560E-04, 1.34840E-04, 1.38712E-04, 1.41746E-04, &
   & 1.39578E-04, 1.37052E-04, 1.33850E-04, 1.26641E-04, 1.21342E-04, &
   & 1.17669E-04, 1.25973E-04, 1.33623E-04, 1.33839E-04, 1.24427E-04, &
   & 1.02462E-04, 8.76101E-05, 8.27912E-05, 8.29040E-05, 7.78590E-05, &
   & 7.39042E-05, 6.45765E-05, 5.70151E-05, 5.11846E-05, 4.83163E-05/
   DATA C20311/                                                      &
   &                                                     5.4470E-05,  &
   & 5.3312E-05,  5.3135E-05,  5.3619E-05,  5.3686E-05,  5.2308E-05,  &
   & 5.0441E-05,  4.8402E-05,  4.7476E-05,  4.6215E-05,  4.4507E-05,  &
   & 4.3830E-05,  4.0508E-05,  3.8931E-05,  3.5525E-05,  3.4722E-05/
   DATA C20327/                                                      &
   & 3.2743E-05,  2.8456E-05,  2.8318E-05,  2.8132E-05,  2.6221E-05,  &
   & 2.5673E-05,  2.5521E-05,  2.4588E-05,  2.4093E-05,  2.2787E-05,  &
   & 2.1241E-05,  1.8553E-05,  1.5871E-05,  1.3462E-05,  1.2553E-05,  &
   & 1.6276E-05,  2.8296E-05,  3.8817E-05,  4.2733E-05,  4.2429E-05,  &
   & 4.0954E-05,  3.9868E-05,  3.7669E-05,  3.6312E-05,  3.5535E-05,  &
   & 3.5895E-05,  3.6349E-05,  3.9033E-05,  4.4512E-05,  5.0066E-05,  &
   & 5.4572E-05,  5.6710E-05,  5.6615E-05,  5.7520E-05,  5.8034E-05,  &
   & 5.7927E-05,  5.6027E-05,  5.5242E-05,  5.4974E-05,  5.2927E-05,  &
   & 5.1638E-05,  5.2027E-05,  5.1420E-05,  5.1618E-05,  5.0253E-05,  &
   & 5.0509E-05,  4.9376E-05,  5.0135E-05,  4.9191E-05,  4.9210E-05,  &
   & 4.8216E-05,  4.7487E-05,  4.5749E-05,  4.5884E-05,  4.3852E-05,  &
   & 4.3824E-05,  4.2612E-05,  4.0349E-05,  4.0177E-05,  3.7474E-05,  &
   & 3.8120E-05,  3.6915E-05,  3.5823E-05,  3.5186E-05,  3.3638E-05,  &
   & 3.3451E-05,  3.2428E-05,  3.2349E-05,  3.0183E-05,  2.8436E-05,  &
   & 2.6440E-05,  2.3597E-05,  2.1875E-05,  1.8164E-05,  1.6430E-05,  &
   & 1.3159E-05,  9.2907E-06,  7.4243E-06,  6.0469E-06,  5.4951E-06/
   DATA C20407/                                                      &
   & 8.7642E-06/
   DATA C20001 /                                                     &
   & 2.16295E-05, 1.69111E-05, 5.39633E-05, 1.01866E-04, 8.28657E-05, &
   & 9.16593E-05, 8.88666E-05, 1.37764E-04, 1.44322E-04, 1.20659E-04, &
   & 1.10332E-04, 1.01317E-04, 9.09964E-05, 1.17148E-04, 1.18000E-04, &
   & 7.21801E-05, 1.10550E-04, 1.32672E-04, 1.02474E-04, 1.10434E-04, &
   & 1.38759E-04, 8.92135E-05, 9.18239E-05, 9.08256E-05, 7.02969E-05, &
   & 1.12827E-04, 8.25561E-05, 1.39555E-04, 6.72239E-05, 7.82804E-05, &
   & 8.56258E-05, 8.61068E-05, 7.16732E-05, 6.25720E-05, 5.23957E-05, &
   & 3.78801E-05, 4.37281E-05, 4.99821E-05, 5.96976E-05, 7.19070E-05, &
   & 3.89579E-05, 5.30171E-05, 3.92507E-05, 4.93901E-05, 4.53047E-05, &
   & 4.89955E-05, 4.61649E-05, 3.75742E-05, 3.14124E-05, 2.37893E-05, &
   & 3.34899E-06, 3.08080E-05, 5.35883E-05, 3.39838E-05, 7.02334E-05, &
   & 7.24784E-05, 7.46533E-05, 6.22257E-05, 6.38945E-05, 6.73423E-05, &
   & 4.51321E-05, 5.91854E-05, 5.51601E-05, 4.41923E-05, 3.59217E-05, &
   & 4.08520E-05, 6.15981E-05, 6.66549E-05, 8.26031E-05, 1.13556E-04, &
   & 8.72988E-05, 9.71052E-05, 9.31839E-05, 8.73745E-05, 8.61717E-05, &
   & 6.05645E-05, 6.51131E-05, 6.93393E-05, 7.01096E-05, 6.43565E-05/
   DATA C20081 /                                                     &
   & 7.36929E-05, 7.66881E-05, 7.60815E-05, 7.13570E-05, 8.40487E-05, &
   & 8.51489E-05, 7.54168E-05, 6.72694E-05, 4.75508E-05, 3.59379E-05, &
   & 4.24698E-05, 4.17850E-05, 4.56047E-05, 4.12779E-05, 4.55933E-05, &
   & 4.27941E-05, 4.42230E-05, 3.68525E-05, 3.83392E-05, 3.83722E-05, &
   & 4.64904E-05, 3.33878E-05, 3.53027E-05, 3.54694E-05, 2.36233E-05, &
   & 2.99641E-05, 2.56097E-05, 2.14134E-05, 2.74403E-05, 2.83896E-05, &
   & 3.17082E-05, 1.75526E-05, 2.80382E-05, 3.18009E-05, 4.08715E-05, &
   & 4.77807E-05, 5.00609E-05, 5.12459E-05, 4.44062E-05, 4.74942E-05, &
   & 4.99882E-05, 5.18837E-05, 5.03246E-05, 5.55168E-05, 5.35853E-05, &
   & 4.81834E-05, 6.66231E-05, 5.26670E-05, 6.84700E-05, 6.53412E-05, &
   & 5.71740E-05, 4.61076E-05, 3.90239E-05, 4.72924E-05, 6.32194E-05, &
   & 5.20868E-05, 4.81039E-05, 3.71748E-05, 4.37492E-05, 3.63959E-05, &
   & 3.79823E-05, 3.72225E-05, 3.02360E-05, 3.22961E-05, 3.43398E-05, &
   & 3.57176E-05, 2.65446E-05, 3.29388E-05, 1.65455E-05, 2.66173E-05, &
   & 1.74277E-05, 1.74324E-05, 1.27879E-05, 1.46247E-05, 1.92378E-05, &
   & 2.20049E-05, 1.44790E-05, 2.49244E-05, 2.29209E-05, 1.76192E-05/
   DATA C20161 /                                                     &
   & 1.84528E-05, 2.54350E-05, 3.33972E-05, 3.69190E-05, 2.92139E-05, &
   & 2.47666E-05, 2.86764E-05, 1.48163E-05, 1.80461E-05, 2.84545E-05, &
   & 2.41064E-05, 2.85721E-05, 3.31996E-05, 3.75973E-05, 3.73874E-05, &
   & 4.69293E-05, 5.12665E-05, 5.35607E-05, 4.64577E-05, 3.59887E-05, &
   & 3.39168E-05, 3.89746E-05, 3.12196E-05, 3.70907E-05, 3.95172E-05, &
   & 4.61642E-05, 4.26029E-05, 4.17856E-05, 4.51437E-05, 4.04189E-05, &
   & 4.19251E-05, 4.53977E-05, 3.69860E-05, 4.20904E-05, 3.69735E-05, &
   & 3.57898E-05, 3.47729E-05, 3.14280E-05, 2.71197E-05, 3.34380E-05, &
   & 2.69843E-05, 2.88036E-05, 2.51912E-05, 2.45699E-05, 2.23184E-05, &
   & 2.50563E-05, 2.24493E-05, 1.77101E-05, 1.64763E-05, 1.34978E-05, &
   & 1.57081E-05, 1.45966E-05, 1.02722E-05, 2.07177E-05, 1.47662E-05, &
   & 1.50721E-05, 1.24431E-05, 1.51572E-05, 1.92210E-05, 2.06047E-05, &
   & 2.02921E-05, 3.22062E-05, 2.37112E-05, 1.94803E-05, 2.40726E-05, &
   & 2.11531E-05, 1.89158E-05, 2.46957E-05, 2.63175E-05, 2.57747E-05, &
   & 2.22047E-05, 2.52755E-05, 2.80848E-05, 3.75157E-05, 4.09915E-05, &
   & 4.04853E-05, 3.21661E-05, 3.15652E-05, 3.21576E-05, 3.67060E-05/
   DATA C20241 /                                                     &
   & 3.13071E-05, 2.84939E-05, 2.71169E-05, 2.99559E-05, 2.94631E-05, &
   & 3.26716E-05, 2.99028E-05, 2.60045E-05, 3.15375E-05, 3.12895E-05, &
   & 2.77767E-05, 2.43976E-05, 2.10764E-05, 2.22725E-05, 2.04581E-05, &
   & 1.63509E-05, 1.60028E-05, 1.60294E-05, 1.62366E-05, 1.89293E-05, &
   & 1.79675E-05, 1.89259E-05, 1.68300E-05, 1.99460E-05, 2.42370E-05, &
   & 2.64738E-05, 1.93137E-05, 1.39460E-05, 1.32222E-05, 1.38752E-05, &
   & 1.62071E-05, 1.79652E-05, 1.63772E-05, 1.56251E-05, 1.81918E-05, &
   & 1.46111E-05, 2.92174E-05, 2.94263E-05, 2.46180E-05, 2.93333E-05, &
   & 3.13657E-05, 2.97686E-05, 2.78387E-05, 2.40924E-05, 2.93369E-05, &
   & 2.93747E-05, 2.77665E-05, 3.00814E-05, 3.01068E-05, 3.62275E-05, &
   & 3.56613E-05, 3.66913E-05, 3.56280E-05, 3.52856E-05, 3.63928E-05, &
   & 2.96738E-05, 2.90314E-05, 2.62972E-05, 2.15250E-05, 1.97910E-05, &
   & 2.02314E-05, 2.20209E-05, 2.05131E-05, 2.12017E-05, 1.96689E-05, &
   & 1.61907E-05, 1.57662E-05, 1.58239E-05, 1.54650E-05, 1.46376E-05, &
   & 1.32891E-05, 1.30511E-05, 1.17635E-05, 1.28585E-05, 1.12887E-05, &
   & 1.32627E-05, 1.31833E-05, 1.68679E-05, 1.98092E-05, 2.70744E-05/
   DATA C20321 /                                                     &
   & 2.22033E-05, 1.63430E-05, 1.61104E-05, 1.50865E-05, 1.54382E-05, &
   & 1.55654E-05, 1.67924E-05, 1.89185E-05, 1.96791E-05, 2.14894E-05, &
   & 2.76137E-05, 2.67339E-05, 2.79423E-05, 2.54664E-05, 3.10707E-05, &
   & 2.72745E-05, 2.60940E-05, 2.47736E-05, 2.21105E-05, 2.20357E-05, &
   & 2.26499E-05, 2.34137E-05, 2.29537E-05, 2.36157E-05, 2.48244E-05, &
   & 2.26667E-05, 2.07781E-05, 2.11702E-05, 1.91214E-05, 1.62172E-05, &
   & 1.61285E-05, 1.63952E-05, 1.68156E-05, 1.61236E-05, 1.56611E-05, &
   & 1.47697E-05, 1.50856E-05, 1.44169E-05, 1.63816E-05, 1.74283E-05, &
   & 1.49853E-05, 1.62444E-05, 1.70007E-05, 1.60371E-05, 1.22713E-05, &
   & 1.45518E-05, 1.35051E-05, 1.40787E-05,-1.54925E-05,-2.15204E-05, &
   &-4.04516E-06, 2.22439E-05, 3.21262E-05, 3.83792E-05, 4.44462E-05, &
   & 4.44192E-05, 2.77328E-05, 4.10549E-06, 4.48758E-06,-1.27771E-05, &
   &-2.17204E-05,-1.23979E-05,-1.04928E-05, 7.43085E-06, 1.55350E-05, &
   & 3.15204E-05, 3.17601E-05, 2.93677E-05, 3.42485E-05, 3.87087E-05, &
   & 3.61242E-05, 2.62406E-05, 3.31686E-05, 3.54314E-05, 2.50625E-05, &
   & 2.60444E-05, 4.10729E-05, 3.47247E-05, 3.31716E-05, 3.34778E-05/
   DATA C20401 /                                                     &
   & 4.03029E-05, 4.09241E-05, 3.96717E-05, 3.53410E-05, 2.81048E-05, &
   & 1.98891E-05, 1.92314E-05, 2.82525E-05, 3.76641E-05, 4.34135E-05, &
   & 4.24570E-05, 3.98429E-05, 3.29417E-05, 2.16679E-05, 8.88085E-06, &
   &-5.05319E-06,-8.14815E-06,-5.01930E-06, 7.13565E-06, 2.00949E-05, &
   & 2.65988E-05, 2.77656E-05, 2.09299E-05, 1.98968E-05, 2.04835E-05, &
   & 1.75254E-05, 6.48674E-06, 3.14323E-06, 1.93242E-06, 3.86745E-06, &
   & 1.39727E-05, 2.10731E-05, 2.66432E-05, 2.69551E-05, 2.57453E-05, &
   & 2.72834E-05, 2.58860E-05, 2.51266E-05, 1.76048E-05, 2.03072E-05, &
   & 2.61960E-05, 2.36230E-05, 1.81172E-05, 1.33972E-05, 1.60959E-05, &
   & 1.61081E-05, 2.34099E-05, 2.64979E-05, 2.36894E-05, 2.13665E-05, &
   & 2.16774E-05, 2.52566E-05, 1.99785E-05, 1.40414E-05, 1.39948E-05, &
   & 1.32637E-05, 7.24742E-06, 1.11395E-06,-1.27323E-06, 4.56637E-07, &
   & 6.93250E-06, 5.07198E-06, 7.90632E-06, 9.08149E-06, 1.03602E-05, &
   & 2.17425E-05, 2.71741E-05, 2.16875E-05, 1.95088E-05, 1.56568E-05, &
   & 8.41152E-06, 1.26749E-05, 1.17673E-05, 9.96037E-06, 1.21982E-05, &
   & 1.31854E-05, 1.50216E-05, 1.72214E-05, 2.02773E-05, 2.09625E-05/
   DATA C20481 /                                                     &
   & 1.66656E-05, 1.45666E-05, 1.66608E-05, 2.04989E-05, 2.21395E-05, &
   & 2.35993E-05, 2.69390E-05, 2.13921E-05, 1.72643E-05, 1.70995E-05, &
   & 1.78241E-05, 1.85308E-05, 1.80360E-05, 1.48619E-05, 1.90369E-05, &
   & 1.51089E-05, 1.22705E-05, 1.62608E-05, 1.41637E-05, 1.23786E-05, &
   & 7.02677E-06, 8.89811E-06, 1.07379E-05, 1.23677E-05, 1.48196E-05, &
   & 2.05770E-05, 1.70994E-05, 1.00072E-05, 1.76119E-05, 1.41779E-05, &
   & 1.34358E-05, 1.58674E-05, 1.65837E-05, 1.69569E-05, 1.40381E-05, &
   & 1.46118E-05, 1.30556E-05, 1.97204E-05, 1.97488E-05, 1.64524E-05, &
   & 1.73764E-05, 1.66355E-05, 1.64419E-05, 1.65486E-05, 1.21523E-05, &
   & 1.51513E-05, 1.60354E-05, 1.38528E-05, 1.45538E-05, 1.71702E-05, &
   & 1.56336E-05, 1.31279E-05, 1.47346E-05, 1.70719E-05, 1.75588E-05, &
   & 1.55187E-05, 1.29598E-05, 1.38463E-05, 1.35382E-05, 1.16062E-05, &
   & 1.37014E-05, 1.34487E-05, 1.15536E-05, 1.33597E-05, 9.24478E-06, &
   & 7.28477E-06, 1.40321E-05, 1.31518E-05, 1.03118E-05, 8.59764E-06, &
   & 1.57138E-05, 1.20792E-05, 1.49440E-05, 1.34375E-05, 1.54686E-05, &
   & 1.65346E-05, 1.33823E-05, 1.37238E-05, 1.36128E-05, 1.46206E-05/
   DATA C20561 /                                                     &
   & 1.40777E-05, 1.59980E-05, 1.30180E-05, 1.01390E-05, 1.12366E-05, &
   & 9.86099E-06, 1.10702E-05, 1.26783E-05, 9.51072E-06, 8.07299E-06, &
   & 1.22955E-05, 1.53506E-05, 1.29711E-05, 9.78759E-06, 1.28800E-05, &
   & 1.39702E-05, 1.64832E-05, 1.06473E-05, 1.15419E-05, 1.63795E-05, &
   & 1.69837E-05, 1.72726E-05, 1.77231E-05, 1.62337E-05, 1.20881E-05, &
   & 1.13210E-05, 1.20531E-05, 1.31374E-05, 1.22259E-05, 1.27802E-05, &
   & 1.38962E-05, 8.87355E-06, 9.42264E-06, 1.02075E-05, 7.91816E-06, &
   & 9.66835E-06, 1.24921E-05, 8.43227E-06, 1.10637E-05, 1.03958E-05, &
   & 9.40996E-06, 1.22922E-05, 1.21088E-05, 1.30116E-05, 1.18776E-05, &
   & 1.42245E-05, 1.34745E-05, 1.11165E-05, 1.29914E-05, 1.29801E-05, &
   & 1.10895E-05, 1.12331E-05, 9.03490E-06, 9.33726E-06, 9.63923E-06, &
   & 1.11299E-05, 9.53481E-06, 1.21708E-05, 1.11951E-05, 7.22558E-06, &
   & 6.66928E-06, 1.08926E-05, 1.07870E-05, 9.23485E-06, 8.50452E-06, &
   & 9.41914E-06, 8.74027E-06, 8.93322E-06, 9.79061E-06, 8.26490E-06, &
   & 8.37630E-06, 1.17064E-05, 1.10176E-05, 1.11587E-05, 9.45563E-06, &
   & 1.18352E-05, 7.79327E-06, 9.22766E-06, 1.01868E-05, 8.23925E-06/
   DATA C20641 /                                                     &
   & 9.23706E-06, 1.04428E-05, 8.80392E-06, 9.37098E-06, 7.43126E-06, &
   & 7.01424E-06, 9.29360E-06, 8.97171E-06, 9.31718E-06, 9.87118E-06, &
   & 8.11419E-06, 8.77416E-06, 9.96927E-06, 8.87533E-06, 9.33163E-06, &
   & 7.41505E-06, 9.39988E-06, 1.17932E-05, 1.03287E-05, 9.17415E-06, &
   & 8.43035E-06, 8.00040E-06, 8.33346E-06, 7.66787E-06, 7.18411E-06, &
   & 1.06236E-05, 1.05559E-05, 8.49187E-06, 9.22472E-06, 8.16512E-06, &
   & 8.35687E-06, 1.06325E-05, 9.80273E-06, 9.01599E-06, 9.20499E-06, &
   & 9.98417E-06, 9.23191E-06, 6.98769E-06, 5.17748E-06, 4.57130E-06, &
   & 8.18492E-06, 9.98095E-06, 7.52148E-06, 1.33038E-05, 8.17630E-06, &
   & 1.02454E-05, 9.62706E-06, 9.44304E-06, 8.86704E-06, 8.88116E-06, &
   & 8.79062E-06, 8.20042E-06, 8.55789E-06, 9.26249E-06, 1.00467E-05, &
   & 7.96012E-06, 9.08773E-06, 1.01481E-05, 8.84360E-06, 7.94928E-06, &
   & 6.68425E-06, 8.56576E-06, 1.05282E-05, 1.10647E-05, 9.91625E-06, &
   & 7.95356E-06, 8.66443E-06, 9.13551E-06, 1.04870E-05, 9.79244E-06, &
   & 1.26214E-05, 8.42148E-06, 8.13468E-06, 1.11338E-05, 1.06780E-05, &
   & 8.54578E-06, 7.82119E-06, 8.33258E-06, 8.23644E-06, 5.95583E-06/
   DATA C20721 /                                                     &
   & 5.85592E-06, 4.05898E-06, 6.39260E-06, 8.43280E-06, 8.76251E-06, &
   & 6.70423E-06, 6.81368E-06, 7.43506E-06, 7.14376E-06, 6.51065E-06, &
   & 5.65633E-06, 5.42394E-06, 7.10817E-06, 4.78831E-06, 6.29380E-06, &
   & 4.87344E-06, 6.81764E-06, 6.51611E-06, 5.70526E-06, 6.50590E-06, &
   & 6.61568E-06, 5.39248E-06, 6.32002E-06, 7.98976E-06, 7.73795E-06, &
   & 4.85788E-06, 5.83443E-06, 6.11694E-06, 5.40408E-06, 5.00946E-06, &
   & 5.62153E-06, 6.30263E-06, 6.05764E-06, 5.53274E-06, 5.80664E-06, &
   & 5.18684E-06, 6.85555E-06, 6.22889E-06, 6.06959E-06, 6.49228E-06, &
   & 5.64064E-06, 4.92690E-06, 5.77661E-06, 7.18450E-06, 7.38658E-06, &
   & 6.77379E-06, 5.74668E-06, 6.68355E-06, 6.13655E-06, 6.43266E-06, &
   & 7.08896E-06, 7.71187E-06, 7.37273E-06, 6.75882E-06, 6.39307E-06, &
   & 4.59520E-06, 5.10323E-06, 5.80178E-06, 6.88172E-06, 6.68825E-06, &
   & 7.50416E-06, 6.14975E-06, 6.51422E-06, 7.74942E-06, 8.11492E-06, &
   & 1.19607E-05, 7.92722E-06, 4.47848E-06, 6.02524E-06, 9.74067E-06, &
   & 1.02429E-05, 8.60819E-06, 8.57044E-06, 1.09196E-05, 1.02048E-05, &
   & 3.86222E-06, 9.26104E-06, 7.33341E-06, 9.08181E-06, 1.05569E-05/
   DATA C20801 /                                                     &
   & 1.06776E-05, 1.10247E-05, 1.04520E-05, 8.78328E-06, 7.60679E-06, &
   & 7.27896E-06, 9.72776E-06, 5.16039E-06, 1.03134E-05, 1.09088E-05, &
   & 8.12575E-06, 7.61685E-06, 8.16346E-06, 5.91269E-06, 3.61448E-06, &
   & 8.74336E-06, 1.03990E-05, 6.25691E-06, 7.04541E-06, 7.94348E-06, &
   & 8.39807E-06, 8.67342E-06, 8.32173E-06, 7.56015E-06, 8.31782E-06, &
   & 6.36556E-06, 6.99328E-06, 6.24490E-06, 6.73080E-06, 6.95852E-06, &
   & 7.55508E-06, 7.74168E-06, 7.90414E-06, 8.94934E-06, 7.99809E-06, &
   & 6.12528E-06, 9.04115E-06, 7.14535E-06, 5.88625E-06, 6.43941E-06, &
   & 7.11566E-06, 7.47425E-06, 8.23805E-06, 6.19919E-06, 7.31614E-06, &
   & 8.24852E-06, 6.82172E-06, 5.45362E-06, 6.66115E-06, 8.44300E-06, &
   & 8.07530E-06, 7.22735E-06, 5.85614E-06, 5.13900E-06, 6.03215E-06, &
   & 6.59491E-06, 4.81592E-06, 4.48587E-06, 7.11525E-06, 8.36201E-06, &
   & 7.11669E-06, 2.80033E-06, 6.50756E-06, 9.43974E-06, 5.22402E-06, &
   & 3.82334E-06, 7.29963E-06, 8.62313E-06, 7.42018E-06, 4.56506E-06, &
   & 5.29972E-06, 5.62787E-06, 4.63852E-06, 5.18329E-06, 7.01884E-06, &
   & 7.24888E-06, 5.18157E-06, 5.40219E-06, 5.92412E-06, 4.97977E-06/
   DATA C20881 /                                                     &
   & 5.29040E-06, 5.33812E-06, 4.76620E-06, 4.65759E-06, 5.10546E-06, &
   & 6.49525E-06, 4.43416E-06, 5.30223E-06, 3.27044E-06, 2.55324E-06, &
   & 4.85017E-06, 7.46556E-06, 8.04448E-06, 5.14009E-06, 6.09755E-06, &
   & 5.38381E-06, 6.41959E-06, 6.59233E-06, 4.83160E-06, 3.81289E-06, &
   & 5.37013E-06, 5.69212E-06, 5.54983E-06, 5.73495E-06, 4.00639E-06, &
   & 2.33817E-06, 2.55751E-06, 3.29627E-06, 3.59845E-06, 6.20623E-06, &
   & 4.47088E-06, 3.49267E-06, 3.09273E-06, 3.32506E-06, 4.83353E-06, &
   & 6.39001E-06, 3.78074E-06, 4.07848E-06, 4.01811E-06, 3.19767E-06, &
   & 3.34053E-06, 4.34246E-06, 3.68003E-06, 3.01090E-06, 3.98545E-06, &
   & 2.72338E-06, 1.90024E-06, 2.77553E-06, 3.73381E-06, 2.58685E-06, &
   & 1.70987E-06,-5.48480E-07, 1.64591E-06, 2.43481E-06, 2.52116E-06, &
   & 2.19316E-06, 1.32392E-06, 1.75370E-06, 2.65409E-07, 2.22278E-06, &
   & 2.53079E-06, 2.87260E-06, 1.87600E-06,-3.84453E-07, 1.80836E-06, &
   & 9.28123E-07, 1.94986E-06, 2.40483E-06, 2.79865E-06, 2.86361E-06, &
   & 2.63868E-06, 3.34704E-06, 3.32132E-06, 2.58463E-06, 2.45684E-06, &
   & 3.35043E-06, 3.19848E-06, 1.73037E-06, 2.98206E-06, 2.77491E-06/
   DATA C20961 /                                                     &
   & 6.51674E-07, 2.52219E-06, 2.97136E-06, 1.96700E-06, 2.29350E-06, &
   & 3.01956E-06, 3.20811E-06, 1.30467E-06, 1.68172E-06, 2.56264E-06, &
   & 2.46167E-06, 1.78221E-06, 2.31647E-06, 2.69480E-06, 2.63619E-06, &
   & 1.81319E-06, 1.83448E-06, 2.23432E-06, 8.14045E-07, 8.75863E-07, &
   & 1.61350E-06, 1.59796E-06, 2.08419E-06, 1.89665E-06, 6.93584E-07, &
   & 1.09880E-06, 3.79031E-07,-3.36470E-07, 1.04326E-06, 1.06497E-06, &
   & 2.15108E-07, 3.28774E-07,-5.17613E-07, 1.27762E-06, 8.22924E-07, &
   & 4.92835E-07, 2.24698E-08,-1.99111E-07, 1.30262E-06,-3.81299E-07, &
   & 9.55084E-07, 2.17641E-07,-6.03874E-08, 8.44121E-07, 1.72391E-06, &
   & 1.66921E-06, 2.19855E-06, 1.17655E-06, 1.79637E-06, 3.31670E-06, &
   & 3.40206E-06, 6.05670E-07, 2.08299E-06, 2.10121E-06, 1.68598E-06, &
   & 2.21155E-06, 2.43221E-06, 5.81282E-08, 1.62613E-06,-5.49850E-07, &
   & 2.14143E-07, 1.21751E-06, 2.30470E-06, 4.27911E-06, 2.96622E-06, &
   & 8.67534E-07, 9.12041E-07, 2.48797E-06, 9.43519E-07,-3.60949E-06, &
   & 2.01928E-06, 1.88873E-06, 8.06749E-07, 7.33519E-07, 1.17440E-06, &
   & 1.69744E-06, 3.64492E-06, 3.11556E-06, 8.89471E-07, 1.93064E-06/
   DATA C21041 /                                                     &
   & 3.02787E-06, 1.92575E-06, 1.73720E-06,-1.32700E-07, 1.41743E-06, &
   & 2.24632E-06, 2.47945E-06, 2.05151E-06,-9.56031E-07, 2.57317E-07, &
   & 3.00980E-06, 3.07981E-06, 2.78202E-06, 3.02555E-06, 5.48784E-09, &
   & 2.37693E-06, 2.90011E-06, 2.93608E-06, 2.14837E-06, 6.55832E-07, &
   & 3.41155E-07,-2.13884E-06, 2.52553E-06, 4.27109E-06, 3.33766E-06, &
   & 3.07708E-06, 2.66405E-06, 3.22850E-06,-5.78879E-07,-6.06194E-07, &
   & 1.72864E-06, 1.57072E-06,-3.39701E-07, 7.21540E-08, 1.67012E-06, &
   & 2.48568E-06, 2.70214E-06, 3.62383E-06, 2.20408E-06, 1.19395E-06, &
   & 1.53825E-06, 2.37511E-06, 2.66754E-06, 1.77020E-06, 5.40420E-07, &
   & 2.01156E-06, 3.27498E-06, 3.04720E-06, 1.96213E-06, 3.71633E-06, &
   & 2.07886E-06, 1.60069E-06, 5.33370E-07, 1.33966E-07, 2.16073E-06, &
   & 8.81457E-07, 1.12880E-06, 2.40509E-06, 2.94252E-06, 2.22899E-06, &
   & 1.80941E-06, 2.68577E-06, 2.44584E-06, 2.51720E-06, 2.64857E-06, &
   & 2.24182E-06, 1.62007E-06, 2.60421E-06, 3.09782E-06, 3.11099E-06, &
   & 3.81513E-06, 6.91606E-06, 3.28767E-06, 3.44175E-06, 4.16771E-06, &
   & 3.75452E-06, 2.21050E-06, 2.99939E-06, 2.86993E-06, 2.47080E-06/
   DATA C21121 /                                                     &
   & 2.33607E-06, 2.68568E-06, 3.39344E-06, 6.09518E-06, 5.10422E-06, &
   & 4.04027E-06, 4.01363E-06, 4.53142E-06, 2.94424E-06, 4.76694E-06, &
   & 6.44206E-06, 7.86435E-06, 8.55564E-06, 6.00857E-06, 5.48073E-06, &
   & 1.56287E-06,-1.16619E-06,-1.85215E-06,-3.04762E-06,-3.45420E-07, &
   & 2.48111E-07,-1.39302E-07,-6.27593E-07,-5.26792E-07, 4.81454E-08, &
   &-3.08631E-08,-1.02976E-06,-1.54919E-06,-9.34044E-07,-1.02507E-06, &
   &-1.39794E-06,-1.15709E-06,-1.04875E-06,-1.64379E-06,-2.97514E-06, &
   &-3.22236E-07,-1.18978E-06,-2.85325E-06,-3.93143E-06,-4.15349E-06, &
   &-2.33228E-06,-3.27125E-06,-2.44987E-06,-1.44460E-06,-3.59727E-06, &
   &-7.18516E-07,-1.53237E-06,-1.53526E-06,-1.56450E-06,-2.91088E-06, &
   &-8.52134E-07,-1.44575E-07,-1.50350E-06,-2.92806E-06,-2.47710E-06, &
   &-9.71202E-07,-9.82754E-07,-1.09924E-06,-6.08199E-07, 3.62885E-07, &
   &-6.67372E-07,-1.00033E-06,-1.12001E-06,-1.06624E-06,-9.23789E-07, &
   &-9.83788E-07,-2.11656E-06,-2.45001E-06,-2.75874E-06,-3.36003E-06, &
   &-3.38364E-06,-2.63747E-06,-3.11047E-06,-3.75258E-06,-3.83211E-06, &
   &-3.52833E-06,-3.48464E-06,-3.77021E-06,-4.26887E-06,-4.23917E-06/
   DATA C21201 /                                                     &
   &-1.42438E-06,-2.48477E-06,-2.84719E-06,-2.70247E-06,-2.50588E-06, &
   &-2.22900E-06,-1.78393E-06,-1.76826E-06,-2.16396E-06,-2.67543E-06, &
   &-2.23706E-06,-2.31793E-06,-2.87590E-06,-3.07803E-06,-2.50493E-06, &
   &-4.54223E-06,-5.15511E-06,-5.39690E-06,-4.89633E-06,-3.33710E-06, &
   &-4.56583E-06,-4.78877E-06,-3.93508E-06,-3.29027E-06,-4.95668E-06, &
   &-6.01801E-06,-5.76016E-06,-5.34657E-06,-5.29080E-06,-5.57133E-06, &
   &-5.73135E-06,-5.39374E-06,-5.09808E-06,-5.12874E-06,-5.20269E-06, &
   &-7.30702E-06,-7.04220E-06,-5.96514E-06,-5.74802E-06,-4.53961E-06, &
   &-4.42127E-06,-4.63922E-06,-4.80622E-06,-4.69659E-06,-5.96786E-06, &
   &-6.29800E-06,-4.75452E-06,-2.85907E-06,-5.33662E-06,-5.31681E-06, &
   &-5.04646E-06,-5.21729E-06,-5.93409E-06,-5.73462E-06,-5.44926E-06, &
   &-6.43325E-06,-7.74451E-06,-7.83147E-06,-5.51568E-06,-7.37048E-06, &
   &-4.25726E-06, 2.32917E-06,-5.61131E-07, 2.05234E-06, 3.74631E-07, &
   &-7.66493E-07, 1.42689E-06,-7.79683E-07, 9.06809E-07, 5.13642E-07, &
   &-1.52504E-06,-2.12058E-06,-2.50316E-06, 1.03637E-08, 5.60002E-07, &
   &-1.48075E-06, 1.94155E-06, 1.91846E-06, 2.78507E-06, 3.90146E-06/
   DATA C21281 /                                                     &
   & 3.61409E-06, 3.23677E-06, 4.00022E-06, 3.19157E-06, 4.03034E-07, &
   &-2.03929E-06, 1.23366E-06, 3.28589E-06, 3.94168E-06, 3.94672E-06, &
   & 3.84619E-06, 2.30400E-07,-2.07799E-06,-1.75115E-06,-5.71958E-07, &
   & 2.33425E-06, 2.01664E-06, 6.05673E-07, 9.57363E-07,-8.89924E-07, &
   &-4.71331E-07, 2.82826E-07, 5.10859E-07, 3.63512E-07, 9.86288E-07, &
   &-4.86309E-07,-2.23163E-06,-1.23370E-06,-2.43131E-07,-2.11498E-06, &
   &-1.56756E-06, 2.70905E-06, 1.87606E-08, 7.83721E-08, 1.58444E-06, &
   & 2.88574E-06, 1.40306E-06, 2.40883E-06, 2.84063E-06, 3.13820E-06, &
   & 3.71016E-06, 3.12975E-06, 3.21981E-06, 2.56191E-06, 1.04624E-06, &
   & 1.87464E-07, 7.25329E-07, 1.03650E-06, 7.23663E-07,-4.18739E-07, &
   & 9.95744E-07,-1.80878E-07,-1.04044E-06, 3.86965E-07,-9.36186E-07, &
   &-4.02271E-07,-2.00231E-07,-5.94965E-07, 4.94038E-07, 3.34585E-07, &
   & 4.82255E-07, 1.12599E-06, 2.11763E-06, 2.66807E-07, 2.29324E-07, &
   & 7.07005E-07, 3.41907E-07,-1.17115E-07, 9.03089E-07, 1.76844E-06, &
   & 1.87134E-06, 2.64057E-06, 4.00395E-07,-4.19679E-07, 6.30769E-07, &
   & 1.02725E-06, 1.05876E-06,-4.08660E-07,-2.32668E-06,-2.73468E-06/
   DATA C21361 /                                                     &
   &-2.40600E-06,-1.81203E-06,-7.96431E-07, 7.40789E-07, 2.73188E-07, &
   & 1.68367E-07,-1.27227E-07,-1.05041E-06,-3.51726E-06,-1.64956E-06, &
   &-5.63840E-07,-1.61242E-06,-1.33264E-06, 1.56604E-06, 2.35083E-06, &
   & 9.26708E-07, 5.41983E-07, 3.54277E-07, 8.53743E-07, 1.54196E-06, &
   & 1.19902E-06, 1.10552E-06, 1.63179E-06, 1.96366E-06, 7.82848E-07, &
   &-3.34741E-08,-7.90842E-07,-6.45131E-07, 1.36158E-06, 1.62453E-06, &
   & 6.68965E-07,-4.86203E-08, 6.83561E-07, 1.89652E-06,-2.80988E-07, &
   &-2.30536E-06,-1.90777E-06, 1.31617E-06, 1.27309E-06, 5.90825E-07, &
   & 5.65686E-07, 1.23631E-07,-1.70279E-06,-1.60768E-06, 9.69543E-07, &
   & 1.01108E-07,-2.02473E-06,-1.75146E-06, 6.33201E-07,-3.59110E-06, &
   &-9.71706E-07, 9.16822E-07, 1.40681E-07,-7.16745E-07,-2.11376E-06, &
   &-1.00951E-06, 2.12465E-06, 1.06982E-06, 1.44032E-06, 1.49692E-06, &
   & 1.07277E-06, 1.37006E-06, 1.66932E-06, 1.75820E-06, 1.41859E-06, &
   &-5.84947E-08, 2.17349E-06, 4.27053E-06, 5.27286E-06, 5.87085E-06, &
   & 2.42692E-06, 2.39305E-06, 6.19573E-06, 5.12518E-06, 1.27171E-06, &
   &-6.81963E-07, 4.16199E-08,-1.36608E-06,-2.53272E-06,-2.37700E-06/
   DATA C21441 /                                                     &
   &-7.96719E-07, 3.85367E-07,-1.08393E-07,-9.04587E-07,-1.54917E-06, &
   &-3.11945E-06,-5.58484E-07, 1.61347E-06, 1.11736E-06, 2.11889E-06, &
   & 2.43534E-06, 1.46709E-06,-1.05429E-06, 1.09978E-06, 7.22493E-07, &
   & 8.53307E-08, 1.22733E-06, 2.99380E-06, 3.62416E-06, 3.81404E-06, &
   & 4.46735E-06, 4.70753E-06, 4.54494E-06, 3.83002E-06, 2.28067E-06, &
   & 2.03102E-06, 2.43844E-06, 2.93132E-06, 2.17555E-06, 3.92919E-06, &
   & 3.53089E-06, 1.61388E-06, 5.09498E-06, 3.40067E-06, 1.58876E-06, &
   & 1.17367E-06, 1.13344E-06, 1.17798E-06, 1.10976E-06, 7.90635E-07, &
   &-4.15989E-07,-1.00581E-06,-9.60236E-07,-1.79111E-07,-5.70733E-07, &
   & 1.49766E-06, 3.44374E-06, 6.45914E-07, 1.00532E-06, 2.01068E-06, &
   & 2.59092E-06, 9.35770E-08, 6.00121E-07, 1.54409E-06, 2.03537E-06, &
   & 8.10358E-07, 1.34126E-06, 1.88873E-06, 1.43283E-06,-2.05029E-07, &
   &-1.09782E-06,-6.56149E-07, 2.01650E-06, 1.84770E-06, 4.39586E-08, &
   &-2.03588E-06,-1.46366E-06,-3.45189E-07, 4.02577E-07, 3.10362E-07, &
   &-2.16073E-06,-1.91861E-06,-2.90520E-07, 2.03692E-06, 3.47996E-06, &
   & 4.21761E-06, 3.89000E-06, 1.86138E-06, 1.56143E-06, 4.88964E-07/
   DATA C21521 /                                                     &
   &-9.28184E-07,-4.34315E-07, 8.74954E-07, 1.58417E-06, 1.36880E-06, &
   & 2.65016E-06, 4.62623E-06, 5.81990E-06, 4.72139E-06, 1.95905E-06, &
   & 1.54151E-06, 2.95768E-06, 4.71536E-06, 2.62359E-06, 9.11513E-07, &
   & 4.75677E-07,-1.53801E-06,-2.32382E-06,-2.25220E-06,-1.46641E-06, &
   &-2.23014E-06,-2.12604E-06,-1.66259E-06,-2.48856E-06,-2.38895E-06, &
   &-2.18158E-06,-1.95841E-06, 4.43899E-07, 1.08517E-06, 1.66370E-07, &
   &-2.42342E-06,-7.19331E-07, 3.19532E-07, 3.58690E-07,-2.01979E-07, &
   & 5.07242E-07, 1.10562E-06, 1.00419E-06, 1.22379E-06, 7.05180E-07, &
   & 1.42283E-07, 8.61092E-07, 8.95236E-07, 1.18043E-07,-1.23589E-06, &
   &-6.16316E-07,-1.18947E-06,-1.45838E-06,-1.47522E-09, 1.33867E-06, &
   & 9.18310E-07,-8.98949E-07,-2.27314E-06,-1.71510E-06,-7.16704E-07, &
   & 8.60666E-09, 5.68015E-07, 1.31219E-06, 1.75478E-06, 5.11790E-07, &
   & 3.35270E-07, 5.39243E-07, 9.08467E-07, 1.39382E-06, 1.08806E-06, &
   & 1.18589E-06, 3.58461E-06, 2.78668E-06, 1.25964E-06,-2.72255E-07, &
   & 1.72305E-06, 1.82937E-06, 7.46252E-07,-1.10555E-06, 2.24967E-07, &
   & 6.45674E-07,-1.87591E-07,-8.84068E-07,-1.75433E-06,-2.17670E-06/
   DATA C21601 /                                                     &
   &-1.37112E-06,-2.31722E-06,-2.23860E-06,-1.16796E-06,-2.23765E-06, &
   &-1.86406E-06,-1.03517E-06,-5.90824E-07,-6.57710E-07,-7.00941E-07, &
   &-4.46064E-07, 1.77205E-06, 2.45066E-06, 2.39371E-06, 2.30736E-06, &
   & 2.35355E-06, 1.85070E-06, 9.62711E-07, 2.59644E-06, 2.05304E-06, &
   & 9.70090E-07, 1.50942E-06, 3.79439E-06, 2.94597E-06,-1.91789E-06, &
   & 6.44324E-08,-3.92094E-07,-1.55398E-06, 4.46701E-08,-4.78760E-07, &
   &-1.70061E-06,-3.17252E-06,-2.93173E-06,-2.01455E-06,-7.76298E-07, &
   &-2.74577E-07,-1.39907E-06,-2.16470E-06,-1.26010E-06,-2.76845E-06, &
   &-2.38226E-06,-5.49068E-08, 9.65258E-07, 1.08650E-06, 5.64738E-07, &
   &-5.78379E-07,-5.68918E-07,-1.90177E-06,-5.08874E-06,-3.03648E-06, &
   &-1.30527E-06,-4.87669E-07,-2.83326E-06,-1.97823E-06,-5.94313E-07, &
   &-1.50961E-07,-1.15908E-06,-1.43260E-06,-9.29331E-07,-1.39459E-06, &
   &-1.27237E-06,-1.50189E-06,-3.79292E-06,-3.92038E-06,-3.58490E-06, &
   &-3.26439E-06,-2.42138E-06,-2.70516E-06,-3.58080E-06,-1.71822E-06, &
   &-2.41567E-06,-3.50193E-06,-2.62394E-06,-3.08424E-06,-3.89604E-06, &
   &-4.84127E-06,-4.41385E-06,-3.22673E-06,-1.80987E-06,-2.93027E-06/
   DATA C21681 /                                                     &
   &-3.17366E-06,-2.79721E-06,-1.78848E-06,-2.80254E-06,-3.55572E-06, &
   &-3.34632E-06,-2.83979E-06,-2.48022E-06,-2.15090E-06,-1.08311E-06, &
   &-6.15216E-07,-7.13008E-07,-1.70841E-06,-2.96098E-06,-3.57134E-06, &
   &-3.04405E-06,-3.35280E-06,-2.97780E-06,-1.97966E-06,-2.33197E-06, &
   &-2.76708E-06,-2.70409E-06,-4.51094E-07,-1.43068E-06,-2.83719E-06, &
   &-2.98921E-06,-4.14949E-06,-3.63780E-06,-8.10138E-07,-1.61597E-06, &
   &-2.25394E-06,-2.58110E-06,-1.57781E-06,-1.71520E-06,-2.30016E-06, &
   &-2.61268E-06,-1.96696E-06,-1.86744E-06,-3.15645E-06,-3.59354E-06, &
   &-3.61015E-06,-3.21793E-06,-2.57436E-06,-2.74347E-06,-3.33319E-06, &
   &-2.93400E-06,-3.25986E-06,-3.46384E-06,-2.22114E-06,-2.92650E-06, &
   &-3.73666E-06,-3.70485E-06,-2.75963E-06,-2.40652E-06,-2.93107E-06, &
   &-1.77517E-06,-1.57096E-06,-2.17533E-06,-2.80190E-06,-2.27942E-06, &
   &-1.37371E-06,-1.65974E-06,-1.26079E-06,-8.08050E-07,-8.41278E-07, &
   &-1.53860E-06,-1.66687E-06,-6.56592E-07,-3.05110E-08, 1.08623E-07, &
   &-2.87222E-07,-2.63555E-07,-7.89575E-07,-1.56059E-06,-6.42174E-07, &
   &-9.43333E-07,-1.38671E-06, 6.50443E-07, 1.35301E-06, 9.27981E-07/
   DATA C21761 /                                                     &
   &-1.21705E-06,-9.63848E-08, 8.73593E-07,-3.47278E-08,-1.79042E-06, &
   &-2.15544E-06,-4.48668E-07,-1.17414E-06,-1.35437E-06,-8.90688E-07, &
   &-4.54757E-07, 2.41484E-09, 3.88010E-07,-1.85349E-08, 1.58011E-07, &
   & 3.70566E-07,-7.30268E-07,-8.42354E-07,-4.13738E-07, 3.96796E-07, &
   &-5.55763E-07,-1.26877E-06,-2.89854E-07, 5.78676E-07, 9.51356E-07, &
   & 5.56912E-07, 1.05014E-06, 9.75896E-07, 5.91573E-08,-6.15073E-07, &
   &-1.48803E-06,-2.53397E-06,-1.77027E-06,-2.08546E-06,-3.10452E-06, &
   &-1.65227E-06,-1.15981E-06,-1.25849E-06,-9.65711E-07,-1.90319E-06, &
   &-2.71831E-06,-5.71559E-08,-1.20368E-06,-3.16820E-06,-2.22766E-06, &
   &-1.19828E-06,-2.82573E-07, 2.53850E-07,-9.10547E-07,-1.65529E-06, &
   &-6.00138E-07,-4.98898E-07,-3.45799E-07, 2.25160E-07, 1.14332E-07, &
   & 3.16082E-07, 1.12681E-06,-6.04876E-07,-7.24616E-07, 1.48177E-06, &
   & 1.05680E-06, 5.91076E-07, 2.07187E-07, 3.82385E-07, 5.91560E-07, &
   & 8.26519E-07, 1.22139E-06, 1.63501E-06, 2.06423E-06, 2.50038E-06, &
   & 2.38037E-06, 1.91688E-06, 2.46702E-06, 2.45066E-06, 2.16732E-06, &
   & 3.13517E-06, 2.68221E-06, 1.39877E-06, 8.58945E-07, 6.83181E-07/
   DATA C21841 /                                                     &
   & 8.46816E-07, 1.73491E-06, 1.98732E-06, 1.94059E-06, 2.19284E-06, &
   & 1.73215E-06, 1.06865E-06, 1.14117E-06, 1.43213E-06, 1.42275E-06, &
   &-4.15449E-07,-2.39911E-07, 3.46498E-08,-2.75022E-06,-2.43736E-06, &
   &-1.06489E-06,-7.81941E-07,-8.04801E-07,-1.04984E-06,-1.65734E-06, &
   &-1.03167E-06,-3.18255E-08, 5.70283E-07, 6.19050E-07, 2.92257E-07, &
   &-6.01436E-07,-7.04005E-07,-3.70875E-07, 4.12830E-07, 1.31319E-07, &
   &-1.61570E-07, 9.76170E-07, 7.99907E-07, 1.41860E-07,-1.98022E-07, &
   & 3.13766E-07, 7.43982E-07,-6.11287E-07,-5.21146E-07, 1.11156E-07, &
   & 3.91719E-07, 5.45566E-07, 6.39059E-07, 7.29515E-07, 4.59167E-07, &
   & 6.13179E-08,-3.48146E-08, 5.32046E-07, 1.19736E-06, 3.83982E-07, &
   & 1.73267E-07, 3.54304E-07, 9.34657E-07, 5.53819E-07,-2.86678E-07, &
   & 2.01853E-08,-1.56159E-07,-6.08130E-07,-2.14929E-07, 1.66317E-08, &
   & 9.32462E-08,-4.83623E-07,-9.16323E-07,-1.22772E-06,-1.61586E-06, &
   &-1.27409E-06,-1.98119E-07,-3.69182E-08,-1.41061E-07,-5.12562E-07, &
   &-4.55495E-07,-8.12132E-07,-1.71772E-06,-2.70741E-06,-2.98751E-06, &
   &-2.19520E-06, 3.01900E-07, 1.17806E-06,-1.23067E-06, 4.17086E-07/
   DATA C21921 /                                                     &
   & 1.68113E-06, 4.81677E-07,-1.55187E-07,-3.35287E-07, 2.94916E-07, &
   & 4.57124E-07, 3.38692E-07,-2.49203E-07,-3.62585E-07,-2.39653E-07, &
   & 3.72675E-08,-7.79964E-09,-2.83285E-07,-9.74713E-07,-6.91171E-07, &
   & 1.21925E-07, 3.39940E-07, 3.68441E-08,-5.82188E-07, 2.12605E-07, &
   & 4.65144E-07, 2.17190E-07, 7.50119E-07, 8.62008E-07, 4.63016E-07, &
   & 1.25620E-06, 1.04567E-06,-8.17037E-07,-1.20023E-06,-1.06224E-06, &
   &-3.77100E-07,-1.28057E-07,-2.76183E-07,-1.24304E-06,-2.56776E-06, &
   &-3.36699E-06,-1.49408E-06,-1.01189E-07, 7.41870E-07,-6.45425E-07, &
   &-7.47111E-07, 4.79055E-10,-1.32339E-06,-1.86135E-06,-1.61074E-06, &
   &-1.82039E-06,-1.68040E-06,-1.08025E-06,-8.61965E-07,-7.00131E-07, &
   &-5.63105E-07,-8.09843E-07,-8.09221E-07, 1.69474E-07,-1.33941E-07, &
   &-7.49558E-07,-5.19013E-07,-8.53534E-07,-1.33703E-06,-3.11161E-07, &
   & 8.99037E-07, 2.25330E-06, 1.44822E-06, 3.07437E-07,-1.22366E-06, &
   &-7.64217E-07, 2.13156E-08, 1.07909E-06, 6.10755E-07, 1.81483E-07, &
   & 8.12405E-07,-9.13283E-08,-1.35885E-06,-1.58366E-06,-7.88594E-07, &
   & 4.48283E-07,-1.23754E-06,-1.65105E-06,-8.93014E-07,-1.48622E-06/
   DATA C22001 /                                                     &
   &-1.67948E-06,-1.24310E-06,-1.54411E-06,-1.65677E-06,-1.04998E-06, &
   &-1.46985E-07, 4.61778E-07,-4.87832E-07,-4.89452E-07,-1.24840E-07, &
   &-1.70101E-06,-1.66976E-06,-1.48528E-07,-1.12621E-07,-2.30607E-08, &
   & 1.82301E-07,-8.58152E-07,-1.89794E-06,-2.46464E-06,-2.32745E-06, &
   &-2.02112E-06,-2.07656E-06,-1.43824E-06,-5.16583E-07,-1.80702E-06, &
   &-2.93490E-06,-3.89216E-06,-3.36211E-06,-2.41393E-06,-9.53406E-07, &
   &-1.16269E-06,-1.66431E-06,-1.77150E-06,-1.82496E-06,-1.93095E-06, &
   &-2.75759E-06,-2.83618E-06,-2.27908E-06,-6.33348E-07, 5.61257E-07, &
   & 1.00142E-06, 7.73337E-07, 3.17721E-07,-3.69804E-07,-8.82058E-07, &
   &-1.17364E-06,-4.53480E-07,-2.47824E-07,-4.79624E-07,-5.17032E-07, &
   &-3.46498E-07, 1.42669E-07,-1.59168E-07,-5.06580E-07,-3.18573E-07, &
   &-2.74092E-07,-2.68860E-07, 1.32811E-07,-2.35567E-09,-6.71971E-07, &
   &-9.75302E-07,-8.70978E-07,-3.59071E-08,-3.01726E-07,-8.27641E-07, &
   &-1.14899E-06,-1.50160E-06,-1.83660E-06,-1.26290E-06,-1.07659E-06, &
   &-1.34878E-06,-5.24626E-07,-7.85094E-08,-8.79473E-07,-1.19291E-06, &
   &-1.33298E-06,-1.59750E-06,-1.31836E-06,-5.73079E-07,-1.10349E-06/
   DATA C22081 /                                                     &
   &-1.11807E-06,-1.99530E-07,-8.10496E-07,-1.42679E-06,-5.34617E-07, &
   &-2.05001E-07,-2.51690E-07,-1.01740E-06,-1.02841E-06,-7.48750E-08, &
   &-1.01770E-06,-1.50413E-06, 1.80898E-07, 3.63788E-07,-1.97900E-07, &
   &-1.16721E-06,-1.05497E-06,-2.07218E-08,-1.90590E-07,-8.25501E-07, &
   &-2.21142E-06,-1.19905E-06, 2.16271E-07,-2.52574E-07,-4.35837E-07, &
   &-3.95272E-07, 5.97065E-08, 2.76639E-07, 9.22569E-08, 1.20142E-07, &
   &-2.95030E-09,-1.08216E-06,-1.32386E-06,-9.62248E-07,-1.99430E-06, &
   &-2.13890E-06,-9.56082E-07,-6.94022E-07,-7.75721E-07,-1.31048E-06, &
   &-1.50080E-06,-1.35873E-06,-7.48378E-07,-4.83436E-07,-4.69624E-07, &
   &-1.51156E-06,-2.48221E-06,-3.30134E-06,-2.79114E-06,-2.08976E-06, &
   &-2.24768E-06,-1.06947E-06, 1.17462E-06,-2.51423E-07,-7.85729E-07, &
   & 5.37467E-07,-9.39876E-08,-1.11303E-06,-7.46860E-07,-9.36220E-07, &
   &-1.59880E-06,-1.61420E-06,-1.54368E-06,-1.41036E-06,-7.20350E-07, &
   & 1.35544E-07, 3.14481E-07, 6.29265E-07, 1.09161E-06,-1.36044E-07, &
   &-1.22932E-06,-1.29847E-06,-3.26429E-06,-6.01062E-06,-2.09945E-06, &
   & 1.26878E-07,-2.88050E-08,-6.82802E-07,-1.39340E-06,-1.82986E-06/
   DATA C22161 /                                                     &
   &-1.67208E-06,-1.07994E-06,-1.89195E-06,-2.10782E-06,-1.04519E-06, &
   &-3.27672E-07, 1.95516E-07, 1.63838E-07,-2.29575E-07,-1.01609E-06, &
   &-2.19286E-06,-2.71850E-06,-9.77485E-07,-1.48830E-06,-3.37826E-06, &
   &-1.59130E-06,-5.74498E-07,-8.27962E-07,-9.92211E-07,-1.14422E-06, &
   &-1.41420E-06,-1.11629E-06,-2.51575E-07, 1.60805E-07, 1.82934E-07, &
   &-7.28868E-07,-2.57062E-07, 1.06520E-06, 4.16488E-07, 2.97049E-08, &
   & 6.62797E-08, 8.29435E-07, 1.29657E-06,-2.27961E-06,-3.40386E-06, &
   &-1.88594E-06,-2.29732E-06,-2.72594E-06,-2.09847E-06,-1.31771E-06, &
   &-4.23693E-07,-4.96348E-07,-9.40209E-07,-2.08707E-06,-1.21368E-06, &
   & 4.79409E-07,-1.12548E-08,-4.57316E-07,-8.40885E-07,-5.03210E-07, &
   &-1.61036E-07,-1.05835E-06,-1.66417E-06,-1.97827E-06,-1.63737E-06, &
   &-1.11711E-06,-3.16081E-07,-6.81746E-07,-1.82599E-06,-1.12895E-06, &
   &-9.19712E-07,-1.91707E-06,-2.14767E-06,-2.03629E-06,-2.86441E-06, &
   &-3.07735E-06,-2.28656E-06,-1.40256E-06,-5.50649E-07,-3.11627E-07, &
   &-7.90261E-07,-2.10728E-06,-1.89739E-06,-1.53762E-06,-2.39947E-06, &
   &-2.28765E-06,-1.27564E-06,-2.15154E-06,-3.17932E-06,-3.84234E-06/
   DATA C22241 /                                                     &
   &-3.65102E-06,-2.84055E-06,-2.48744E-06,-2.27683E-06,-2.33087E-06, &
   &-3.44460E-06,-5.19613E-06,-2.85882E-06,-1.39921E-06,-2.00579E-06, &
   &-2.80593E-06,-3.65940E-06,-2.39526E-06,-1.70389E-06,-2.03532E-06, &
   &-2.71522E-06,-3.42227E-06,-2.23606E-06,-1.77845E-06,-2.42071E-06, &
   &-2.61515E-06,-2.56413E-06,-1.49601E-06,-1.23245E-06,-2.08440E-06, &
   &-2.11121E-06,-1.93424E-06,-2.27439E-06,-2.58183E-06,-2.84705E-06, &
   &-2.32183E-06,-1.80966E-06,-3.04089E-06,-3.14334E-06,-1.91331E-06, &
   &-1.51037E-06,-1.43610E-06,-2.11316E-06,-2.45184E-06,-2.42262E-06/
!
end block data BO3HH2
!
!     --------------------------------------------------------------
!
SUBROUTINE O3HHUV (V1C,V2C,DVC,NPTC,C,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   COMMON /O3HUV/ V1S,V2S,DVS,NPTS,S(133)
   DIMENSION C(*)
!
   DVC = DVS
   v1ss = v1s
   v2ss = v2s
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      VJ = V1C+DVC* REAL(J-1)
      C(J) = S(I)/VJ
!
!     RADIATION FLD REMOVED FROM U.V.    OZONE
!
10 END DO
!
   RETURN
!
end subroutine O3HHUV
!
!     --------------------------------------------------------------
!
BLOCK DATA BO3HUV
!
   IMPLICIT REAL*8           (V)
!
!     DATA DERIVED FROM MOLINA & MOLINA, JGR,91,14501-14508,1986.
!     VALUES BETWEEN 245 AND 185NM (40800 AND 54054CM-1) USED AS
!     DIRECT AVERAGE WITH NO TEMPERATURE DEPENDENCE.
!
   COMMON /O3HUV/ V1C,V2C,DVC,NC,                                    &
   &               C02281(80),C02361(53)
!
   DATA V1C /40800./, V2C /54000./ ,DVC /100./, NC /133/
!
   DATA C02281/                                                      &
   & 9.91204E-18, 9.76325E-18, 9.72050E-18, 9.51049E-18, 9.23530E-18, &
   & 9.02306E-18, 8.90510E-18, 8.60115E-18, 8.39094E-18, 8.27926E-18, &
   & 7.95525E-18, 7.73583E-18, 7.55018E-18, 7.31076E-18, 7.10415E-18, &
   & 6.87747E-18, 6.66639E-18, 6.39484E-18, 6.27101E-18, 6.01019E-18, &
   & 5.77594E-18, 5.60403E-18, 5.40837E-18, 5.21289E-18, 4.99329E-18, &
   & 4.81742E-18, 4.61608E-18, 4.45707E-18, 4.28261E-18, 4.09672E-18, &
   & 3.93701E-18, 3.77835E-18, 3.61440E-18, 3.45194E-18, 3.30219E-18, &
   & 3.15347E-18, 3.01164E-18, 2.87788E-18, 2.74224E-18, 2.61339E-18, &
   & 2.48868E-18, 2.36872E-18, 2.25747E-18, 2.14782E-18, 2.03997E-18, &
   & 1.94281E-18, 1.84525E-18, 1.75275E-18, 1.67151E-18, 1.58813E-18, &
   & 1.50725E-18, 1.43019E-18, 1.35825E-18, 1.28878E-18, 1.22084E-18, &
   & 1.15515E-18, 1.09465E-18, 1.03841E-18, 9.83780E-19, 9.31932E-19, &
   & 8.83466E-19, 8.38631E-19, 7.96631E-19, 7.54331E-19, 7.13805E-19, &
   & 6.78474E-19, 6.44340E-19, 6.13104E-19, 5.81777E-19, 5.53766E-19, &
   & 5.27036E-19, 5.03555E-19, 4.82633E-19, 4.61483E-19, 4.42014E-19, &
   & 4.23517E-19, 4.07774E-19, 3.93060E-19, 3.80135E-19, 3.66348E-19/
   DATA C02361/                                                      &
   & 3.53665E-19, 3.47884E-19, 3.39690E-19, 3.34288E-19, 3.29135E-19, &
   & 3.23104E-19, 3.18875E-19, 3.16800E-19, 3.15925E-19, 3.12932E-19, &
   & 3.12956E-19, 3.15522E-19, 3.14950E-19, 3.15924E-19, 3.19059E-19, &
   & 3.23109E-19, 3.27873E-19, 3.33788E-19, 3.39804E-19, 3.44925E-19, &
   & 3.50502E-19, 3.55853E-19, 3.59416E-19, 3.68933E-19, 3.78284E-19, &
   & 3.86413E-19, 3.98049E-19, 4.04700E-19, 4.12958E-19, 4.23482E-19, &
   & 4.31203E-19, 4.41885E-19, 4.52651E-19, 4.61492E-19, 4.70493E-19, &
   & 4.80497E-19, 4.90242E-19, 4.99652E-19, 5.10316E-19, 5.21510E-19, &
   & 5.32130E-19, 5.43073E-19, 5.56207E-19, 5.61756E-19, 5.66799E-19, &
   & 5.85545E-19, 5.92409E-19, 5.96168E-19, 6.12497E-19, 6.20231E-19, &
   & 6.24621E-19, 6.34160E-19, 6.43622E-19/
!
end block data BO3HUV
!
!     --------------------------------------------------------------

subroutine o2_ver_1 (v1c,v2c,dvc,nptc,c,T,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8 (v)

   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)

   COMMON /o2_f  / V1S,V2S,DVS,NPTS,xo2(103),xo2t(103)

   dimension c(*)

!
!     Oxygen Collision Induced Fundamental

!     F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, J.-M. Hartmann,
!     and Ch. Boulet
!     Infrared collision-induced absorption by O2 near 6.4 microns for
!     atmospheric applications: measurements and emprirical modeling,
!     Appl. Optics, 35, 5911-5917, (1996).

   DATA T_0/ 296./, xlosmt/ 2.68675e+19/
!
   xktfac = (1./T_0)-(1./T)
!
!     correct formulation for consistency with LBLRTM:
!
   factor = (1.e+20 /xlosmt)
!
!     A factor of 0.21, the mixing ratio of oxygen, in the Thibault et
!     al. formulation is not included here.  This factor is in the
!     column amount.
!
   DVC = DVS
   v1ss = v1s
   v2ss = v2s
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   do 10 j=1,nptc
      i = i1+(j-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      VJ = V1C+DVC* REAL(J-1)
!     the radiation field is removed with 1/vj
!
      c(j) = factor * xo2(i)* exp(xo2t(i)*xktfac) / vj
!
10 end do
!
920 format (f10.2,1p,e12.2,0p,f10.2,1p2e12.2)
!
   return
end subroutine o2_ver_1

BLOCK DATA bo2f

   IMPLICIT REAL*8 (V)

   COMMON /o2_f  / V1S,V2S,DVS,NPTS,                                 &
   &          o0001(50),o0051(50),o0101(03),                          &
   &          ot0001(50),ot0051(50),ot0101(03)

   DATA V1S,V2S,DVS,NPTS /1340.000,1850.000,   5.000,  103/

   DATA o0001/                                                       &
   &      0.000E+00,  9.744E-09,  2.256E-08,  3.538E-08,  4.820E-08,  &
   &      6.100E-08,  7.400E-08,  8.400E-08,  9.600E-08,  1.200E-07,  &
   &      1.620E-07,  2.080E-07,  2.460E-07,  2.850E-07,  3.140E-07,  &
   &      3.800E-07,  4.440E-07,  5.000E-07,  5.710E-07,  6.730E-07,  &
   &      7.680E-07,  8.530E-07,  9.660E-07,  1.100E-06,  1.210E-06,  &
   &      1.330E-06,  1.470E-06,  1.590E-06,  1.690E-06,  1.800E-06,  &
   &      1.920E-06,  2.040E-06,  2.150E-06,  2.260E-06,  2.370E-06,  &
   &      2.510E-06,  2.670E-06,  2.850E-06,  3.070E-06,  3.420E-06,  &
   &      3.830E-06,  4.200E-06,  4.450E-06,  4.600E-06,  4.530E-06,  &
   &      4.280E-06,  3.960E-06,  3.680E-06,  3.480E-06,  3.350E-06/
   DATA o0051/                                                       &
   &      3.290E-06,  3.250E-06,  3.230E-06,  3.230E-06,  3.210E-06,  &
   &      3.190E-06,  3.110E-06,  3.030E-06,  2.910E-06,  2.800E-06,  &
   &      2.650E-06,  2.510E-06,  2.320E-06,  2.130E-06,  1.930E-06,  &
   &      1.760E-06,  1.590E-06,  1.420E-06,  1.250E-06,  1.110E-06,  &
   &      9.900E-07,  8.880E-07,  7.910E-07,  6.780E-07,  5.870E-07,  &
   &      5.240E-07,  4.640E-07,  4.030E-07,  3.570E-07,  3.200E-07,  &
   &      2.900E-07,  2.670E-07,  2.420E-07,  2.150E-07,  1.820E-07,  &
   &      1.600E-07,  1.460E-07,  1.280E-07,  1.030E-07,  8.700E-08,  &
   &      8.100E-08,  7.100E-08,  6.400E-08,  5.807E-08,  5.139E-08,  &
   &      4.496E-08,  3.854E-08,  3.212E-08,  2.569E-08,  1.927E-08/
   DATA o0101/                                                       &
   &      1.285E-08,  6.423E-09,  0.000E+00/

   DATA ot0001/                                                      &
   &      4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  &
   &      4.670E+02,  4.000E+02,  3.150E+02,  3.790E+02,  3.680E+02,  &
   &      4.750E+02,  5.210E+02,  5.310E+02,  5.120E+02,  4.420E+02,  &
   &      4.440E+02,  4.300E+02,  3.810E+02,  3.350E+02,  3.240E+02,  &
   &      2.960E+02,  2.480E+02,  2.150E+02,  1.930E+02,  1.580E+02,  &
   &      1.270E+02,  1.010E+02,  7.100E+01,  3.100E+01, -6.000E+00,  &
   &     -2.600E+01, -4.700E+01, -6.300E+01, -7.900E+01, -8.800E+01,  &
   &     -8.800E+01, -8.700E+01, -9.000E+01, -9.800E+01, -9.900E+01,  &
   &     -1.090E+02, -1.340E+02, -1.600E+02, -1.670E+02, -1.640E+02,  &
   &     -1.580E+02, -1.530E+02, -1.510E+02, -1.560E+02, -1.660E+02/
   DATA ot0051/                                                      &
   &     -1.680E+02, -1.730E+02, -1.700E+02, -1.610E+02, -1.450E+02,  &
   &     -1.260E+02, -1.080E+02, -8.400E+01, -5.900E+01, -2.900E+01,  &
   &      4.000E+00,  4.100E+01,  7.300E+01,  9.700E+01,  1.230E+02,  &
   &      1.590E+02,  1.980E+02,  2.200E+02,  2.420E+02,  2.560E+02,  &
   &      2.810E+02,  3.110E+02,  3.340E+02,  3.190E+02,  3.130E+02,  &
   &      3.210E+02,  3.230E+02,  3.100E+02,  3.150E+02,  3.200E+02,  &
   &      3.350E+02,  3.610E+02,  3.780E+02,  3.730E+02,  3.380E+02,  &
   &      3.190E+02,  3.460E+02,  3.220E+02,  2.910E+02,  2.900E+02,  &
   &      3.500E+02,  3.710E+02,  5.040E+02,  4.000E+02,  4.000E+02,  &
   &      4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02/
   DATA ot0101/                                                      &
   &      4.000E+02,  4.000E+02,  4.000E+02/

end block data bo2f

!     --------------------------------------------------------------
!
SUBROUTINE O2INF1 (V1C,V2C,DVC,NPTC,C,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   DIMENSION C(*)

   COMMON /o2inf1_mate/ V1S,V2S,DVS,NPTS,xo2inf1(483)

!
!        O2 continuum formulated by Mate et al. over the spectral region
!        7550-8486 cm-1:  "Absolute Intensities for the O2 1.27 micron
!        continuum absorption", B. Mate, C. Lugez, G.T. Fraser, and
!        W.J. Lafferty, J. Geophys. Res., 104, 30,585-30,590, 1999.
!
!        The units of these continua coefficients are
!         1 / (amagat_O2*amagat_air).
!        Also, refer to the paper "Observed  Atmospheric
!        Collision Induced Absorption in Near Infrared Oxygen Bands",
!        Mlawer, Clough, Brown, Stephen, Landry, Goldman, & Murcray,
!        Journal of Geophysical Research (1997).
!   ***********

   DVC = DVS
   v1ss = v1s
   v2ss = v2s
!
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      vj = v1c + dvc* REAL(j-1)
      C(J) = xo2inf1(I)/vj
10 END DO
!
   RETURN
!
end subroutine O2INF1

!     --------------------------------------------------------------
!
BLOCK DATA bo2inf1

   IMPLICIT REAL*8 (V)

   COMMON /o2inf1_mate/ V1,V2,DV,NPT,                                &
   &          o0001(50),o0051(50),o0101(50),o0151(50),o0201(50),      &
   &          o0251(50),o0301(50),o0351(50),o0401(50),o0451(33)

   DATA V1,V2,DV,NPT /7536.000,8500.000,   2.000,  483/

   DATA o0001/                                                       &
   &      0.000E+00,  4.355E-11,  8.709E-11,  1.742E-10,  3.484E-10,  &
   &      6.968E-10,  1.394E-09,  2.787E-09,  3.561E-09,  3.314E-09,  &
   &      3.368E-09,  3.435E-09,  2.855E-09,  3.244E-09,  3.447E-09,  &
   &      3.891E-09,  4.355E-09,  3.709E-09,  4.265E-09,  4.772E-09,  &
   &      4.541E-09,  4.557E-09,  4.915E-09,  4.688E-09,  5.282E-09,  &
   &      5.755E-09,  5.096E-09,  5.027E-09,  4.860E-09,  4.724E-09,  &
   &      5.048E-09,  5.248E-09,  5.473E-09,  4.852E-09,  5.362E-09,  &
   &      6.157E-09,  6.150E-09,  6.347E-09,  6.388E-09,  6.213E-09,  &
   &      6.521E-09,  8.470E-09,  8.236E-09,  8.269E-09,  8.776E-09,  &
   &      9.122E-09,  9.189E-09,  9.778E-09,  8.433E-09,  9.964E-09/
   DATA o0051/                                                       &
   &      9.827E-09,  1.064E-08,  1.063E-08,  1.031E-08,  1.098E-08,  &
   &      1.156E-08,  1.295E-08,  1.326E-08,  1.467E-08,  1.427E-08,  &
   &      1.452E-08,  1.456E-08,  1.554E-08,  1.605E-08,  1.659E-08,  &
   &      1.754E-08,  1.757E-08,  1.876E-08,  1.903E-08,  1.876E-08,  &
   &      1.869E-08,  2.036E-08,  2.203E-08,  2.221E-08,  2.284E-08,  &
   &      2.288E-08,  2.394E-08,  2.509E-08,  2.663E-08,  2.720E-08,  &
   &      2.839E-08,  2.923E-08,  2.893E-08,  2.949E-08,  2.962E-08,  &
   &      3.057E-08,  3.056E-08,  3.364E-08,  3.563E-08,  3.743E-08,  &
   &      3.813E-08,  3.946E-08,  4.082E-08,  4.201E-08,  4.297E-08,  &
   &      4.528E-08,  4.587E-08,  4.704E-08,  4.962E-08,  5.115E-08/
   DATA o0101/                                                       &
   &      5.341E-08,  5.365E-08,  5.557E-08,  5.891E-08,  6.084E-08,  &
   &      6.270E-08,  6.448E-08,  6.622E-08,  6.939E-08,  7.233E-08,  &
   &      7.498E-08,  7.749E-08,  8.027E-08,  8.387E-08,  8.605E-08,  &
   &      8.888E-08,  9.277E-08,  9.523E-08,  9.880E-08,  1.037E-07,  &
   &      1.076E-07,  1.114E-07,  1.151E-07,  1.203E-07,  1.246E-07,  &
   &      1.285E-07,  1.345E-07,  1.408E-07,  1.465E-07,  1.519E-07,  &
   &      1.578E-07,  1.628E-07,  1.685E-07,  1.760E-07,  1.847E-07,  &
   &      1.929E-07,  2.002E-07,  2.070E-07,  2.177E-07,  2.262E-07,  &
   &      2.365E-07,  2.482E-07,  2.587E-07,  2.655E-07,  2.789E-07,  &
   &      2.925E-07,  3.023E-07,  3.153E-07,  3.296E-07,  3.409E-07/
   DATA o0151/                                                       &
   &      3.532E-07,  3.680E-07,  3.859E-07,  3.951E-07,  4.074E-07,  &
   &      4.210E-07,  4.381E-07,  4.588E-07,  4.792E-07,  4.958E-07,  &
   &      5.104E-07,  5.271E-07,  5.501E-07,  5.674E-07,  5.913E-07,  &
   &      6.243E-07,  6.471E-07,  6.622E-07,  6.831E-07,  6.987E-07,  &
   &      7.159E-07,  7.412E-07,  7.698E-07,  7.599E-07,  7.600E-07,  &
   &      7.918E-07,  8.026E-07,  8.051E-07,  8.049E-07,  7.914E-07,  &
   &      7.968E-07,  7.945E-07,  7.861E-07,  7.864E-07,  7.741E-07,  &
   &      7.675E-07,  7.592E-07,  7.400E-07,  7.362E-07,  7.285E-07,  &
   &      7.173E-07,  6.966E-07,  6.744E-07,  6.597E-07,  6.413E-07,  &
   &      6.265E-07,  6.110E-07,  5.929E-07,  5.717E-07,  5.592E-07/
   DATA o0201/                                                       &
   &      5.411E-07,  5.235E-07,  5.061E-07,  4.845E-07,  4.732E-07,  &
   &      4.593E-07,  4.467E-07,  4.328E-07,  4.161E-07,  4.035E-07,  &
   &      3.922E-07,  3.820E-07,  3.707E-07,  3.585E-07,  3.475E-07,  &
   &      3.407E-07,  3.317E-07,  3.226E-07,  3.134E-07,  3.016E-07,  &
   &      2.969E-07,  2.894E-07,  2.814E-07,  2.749E-07,  2.657E-07,  &
   &      2.610E-07,  2.536E-07,  2.467E-07,  2.394E-07,  2.337E-07,  &
   &      2.302E-07,  2.241E-07,  2.191E-07,  2.140E-07,  2.093E-07,  &
   &      2.052E-07,  1.998E-07,  1.963E-07,  1.920E-07,  1.862E-07,  &
   &      1.834E-07,  1.795E-07,  1.745E-07,  1.723E-07,  1.686E-07,  &
   &      1.658E-07,  1.629E-07,  1.595E-07,  1.558E-07,  1.523E-07/
   DATA o0251/                                                       &
   &      1.498E-07,  1.466E-07,  1.452E-07,  1.431E-07,  1.408E-07,  &
   &      1.381E-07,  1.362E-07,  1.320E-07,  1.298E-07,  1.262E-07,  &
   &      1.247E-07,  1.234E-07,  1.221E-07,  1.197E-07,  1.176E-07,  &
   &      1.142E-07,  1.121E-07,  1.099E-07,  1.081E-07,  1.073E-07,  &
   &      1.061E-07,  1.041E-07,  1.019E-07,  9.969E-08,  9.727E-08,  &
   &      9.642E-08,  9.487E-08,  9.318E-08,  9.116E-08,  9.046E-08,  &
   &      8.827E-08,  8.689E-08,  8.433E-08,  8.324E-08,  8.204E-08,  &
   &      8.036E-08,  7.951E-08,  7.804E-08,  7.524E-08,  7.392E-08,  &
   &      7.227E-08,  7.176E-08,  6.975E-08,  6.914E-08,  6.859E-08,  &
   &      6.664E-08,  6.506E-08,  6.368E-08,  6.262E-08,  6.026E-08/
   DATA o0301/                                                       &
   &      6.002E-08,  5.866E-08,  5.867E-08,  5.641E-08,  5.589E-08,  &
   &      5.499E-08,  5.309E-08,  5.188E-08,  5.139E-08,  4.991E-08,  &
   &      4.951E-08,  4.833E-08,  4.640E-08,  4.524E-08,  4.479E-08,  &
   &      4.304E-08,  4.228E-08,  4.251E-08,  4.130E-08,  3.984E-08,  &
   &      3.894E-08,  3.815E-08,  3.732E-08,  3.664E-08,  3.512E-08,  &
   &      3.463E-08,  3.503E-08,  3.218E-08,  3.253E-08,  3.107E-08,  &
   &      2.964E-08,  2.920E-08,  2.888E-08,  2.981E-08,  2.830E-08,  &
   &      2.750E-08,  2.580E-08,  2.528E-08,  2.444E-08,  2.378E-08,  &
   &      2.413E-08,  2.234E-08,  2.316E-08,  2.199E-08,  2.088E-08,  &
   &      1.998E-08,  1.920E-08,  1.942E-08,  1.859E-08,  1.954E-08/
   DATA o0351/                                                       &
   &      1.955E-08,  1.749E-08,  1.720E-08,  1.702E-08,  1.521E-08,  &
   &      1.589E-08,  1.469E-08,  1.471E-08,  1.543E-08,  1.433E-08,  &
   &      1.298E-08,  1.274E-08,  1.226E-08,  1.204E-08,  1.201E-08,  &
   &      1.298E-08,  1.220E-08,  1.220E-08,  1.096E-08,  1.080E-08,  &
   &      9.868E-09,  9.701E-09,  1.130E-08,  9.874E-09,  9.754E-09,  &
   &      9.651E-09,  9.725E-09,  8.413E-09,  7.705E-09,  7.846E-09,  &
   &      8.037E-09,  9.163E-09,  8.098E-09,  8.160E-09,  7.511E-09,  &
   &      7.011E-09,  6.281E-09,  6.502E-09,  7.323E-09,  7.569E-09,  &
   &      5.941E-09,  5.867E-09,  5.676E-09,  4.840E-09,  5.063E-09,  &
   &      5.207E-09,  4.917E-09,  5.033E-09,  5.356E-09,  3.795E-09/
   DATA o0401/                                                       &
   &      4.983E-09,  4.600E-09,  3.635E-09,  3.099E-09,  2.502E-09,  &
   &      3.823E-09,  3.464E-09,  4.332E-09,  3.612E-09,  3.682E-09,  &
   &      3.709E-09,  3.043E-09,  3.593E-09,  3.995E-09,  4.460E-09,  &
   &      3.583E-09,  3.290E-09,  3.132E-09,  2.812E-09,  3.109E-09,  &
   &      3.874E-09,  3.802E-09,  4.024E-09,  3.901E-09,  2.370E-09,  &
   &      1.821E-09,  2.519E-09,  4.701E-09,  3.855E-09,  4.685E-09,  &
   &      5.170E-09,  4.387E-09,  4.148E-09,  4.043E-09,  3.545E-09,  &
   &      3.392E-09,  3.609E-09,  4.635E-09,  3.467E-09,  2.558E-09,  &
   &      3.389E-09,  2.672E-09,  2.468E-09,  1.989E-09,  2.816E-09,  &
   &      4.023E-09,  2.664E-09,  2.219E-09,  3.169E-09,  1.654E-09/
   DATA o0451/                                                       &
   &      3.189E-09,  2.535E-09,  2.618E-09,  3.265E-09,  2.138E-09,  &
   &      1.822E-09,  2.920E-09,  2.002E-09,  1.300E-09,  3.764E-09,  &
   &      3.212E-09,  3.222E-09,  2.961E-09,  2.108E-09,  1.708E-09,  &
   &      2.636E-09,  2.937E-09,  2.939E-09,  2.732E-09,  2.218E-09,  &
   &      1.046E-09,  6.419E-10,  1.842E-09,  1.112E-09,  1.265E-09,  &
   &      4.087E-09,  2.044E-09,  1.022E-09,  5.109E-10,  2.554E-10,  &
   &      1.277E-10,  6.386E-11,  0.000E+00/

end block data bo2inf1


!     --------------------------------------------------------------
!
SUBROUTINE O2INF2 (V1C,V2C,DVC,NPTC,C,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   DIMENSION C(*)
!
   DATA V1_osc /9375./, HW1 /58.96/, V2_osc /9439./, HW2 /45.04/
   DATA S1 /1.166E-04/, S2 /3.086E-05/
!
   V1S = 9100.
   v2s = 11000.
   DVS = 2.
   DVC = DVS
   v1ss = v1s
   v2ss = v2s
!
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
! The following lines prevent a possible problem that can only occur if
! v2abs-v1abs >> 2000 cm-1 (i.e. in the standalone continuum code).
   if (v1c .lt. v1s) v1c = v1s - 2. * dvs
   if (v2c .gt. v2s) v2c = v2s + 2. * dvs
   NPTC = (v2c-v1c)/dvc + 3.01
   V2C = V1C+DVc* REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      C(J) = 0.
      VJ = V1C+DVC* REAL(J-1)
      IF ((Vj.gt.v1s) .and. (Vj.lt.v2s)) then
         DV1 = Vj - V1_osc
         DV2 = Vj - V2_osc
         IF (DV1 .LT. 0.0) THEN
            DAMP1 = EXP (DV1 / 176.1)
         ELSE
            DAMP1 = 1.0
         ENDIF
         IF (DV2 .LT. 0.0) THEN
            DAMP2 = EXP (DV2 / 176.1)
         ELSE
            DAMP2 = 1.0
         ENDIF
         O2INF = 0.31831 * (((S1 * DAMP1 / HW1)/(1. + (DV1/HW1)**2)) &
            + ((S2 * DAMP2 / HW2)/(1. + (DV2/HW2)**2))) * 1.054
         C(J) = O2INF/VJ
      endif
10 END DO
!
   RETURN
!
end subroutine O2INF2
!     --------------------------------------------------------------
!
SUBROUTINE O2INF3 (V1C,V2C,DVC,NPTC,C,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   DIMENSION C(*)

   COMMON /o2inf3_aband/ V1S,V2S,DVS,NPTS,xo2inf3(261)
!
!        O2 A-band continuum formulated by Mlawer based on solar FTS measurements.
!        See Payne et al. (2020) "Absorption Coefficient (ABSCO) Tables for the
!        Orbiting Carbon Observatories: Version 5.1", JQSRT
!
!        Units of these coefficients are 1 / (amagat_O2*amagat_air)
!        Spectral range of coefficients is 12961.5 - 13221.5 cm-1.
!
!   ***********

   DVC = DVS
   v1ss = v1s
   v2ss = v2s
!
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      vj = v1c + dvc* REAL(j-1)
      C(J) = xo2inf3(I)/vj
10 END DO
!
   RETURN
!
end subroutine O2INF3
!     --------------------------------------------------------------
!
BLOCK DATA bo2inf3

   IMPLICIT REAL*8 (V)

   COMMON /o2inf3_aband/ V1,V2,DV,NPT,x02inf3(261)

   DATA V1,V2,DV,NPT /12961.5, 13221.5, 1.0, 261/

   DATA x02inf3/                                                     &
   &      0.000e+00,  1.253e-10,  2.785e-10,  4.316e-10,  5.848e-10,  &
   &      7.379e-10,  8.911e-10,  1.044e-09,  1.197e-09,  1.351e-09,  &
   &      1.504e-09,  1.657e-09,  1.810e-09,  1.963e-09,  2.116e-09,  &
   &      2.269e-09,  2.423e-09,  2.576e-09,  2.729e-09,  2.882e-09,  &
   &      3.035e-09,  3.188e-09,  3.342e-09,  3.495e-09,  3.648e-09,  &
   &      3.801e-09,  3.954e-09,  4.107e-09,  4.260e-09,  4.414e-09,  &
   &      4.567e-09,  4.720e-09,  4.873e-09,  5.026e-09,  5.179e-09,  &
   &      5.333e-09,  5.486e-09,  5.639e-09,  5.792e-09,  5.945e-09,  &
   &      6.098e-09,  6.251e-09,  6.405e-09,  6.552e-09,  6.738e-09,  &
   &      7.008e-09,  7.216e-09,  7.353e-09,  7.427e-09,  7.687e-09,  &
   &      7.917e-09,  8.114e-09,  8.273e-09,  8.403e-09,  8.496e-09,  &
   &      8.591e-09,  8.946e-09,  9.307e-09,  9.641e-09,  9.988e-09,  &
   &      1.033e-08,  1.065e-08,  1.105e-08,  1.156e-08,  1.221e-08,  &
   &      1.301e-08,  1.395e-08,  1.490e-08,  1.588e-08,  1.688e-08,  &
   &      1.783e-08,  1.888e-08,  2.000e-08,  2.114e-08,  2.246e-08,  &
   &      2.374e-08,  2.498e-08,  2.631e-08,  2.760e-08,  2.899e-08,  &
   &      3.043e-08,  3.188e-08,  3.349e-08,  3.543e-08,  3.732e-08,  &
   &      3.912e-08,  4.081e-08,  4.252e-08,  4.454e-08,  4.628e-08,  &
   &      4.784e-08,  4.935e-08,  5.103e-08,  5.292e-08,  5.486e-08,  &
   &      5.702e-08,  5.920e-08,  6.045e-08,  6.317e-08,  6.630e-08,  &
   &      6.945e-08,  7.282e-08,  7.599e-08,  8.020e-08,  8.432e-08,  &
   &      8.814e-08,  9.196e-08,  9.582e-08,  9.978e-08,  1.038e-07,  &
   &      1.078e-07,  1.115e-07,  1.143e-07,  1.182e-07,  1.222e-07,  &
   &      1.264e-07,  1.307e-07,  1.349e-07,  1.389e-07,  1.432e-07,  &
   &      1.488e-07,  1.543e-07,  1.596e-07,  1.650e-07,  1.702e-07,  &
   &      1.753e-07,  1.795e-07,  1.835e-07,  1.885e-07,  1.945e-07,  &
   &      2.005e-07,  2.056e-07,  2.103e-07,  2.149e-07,  2.195e-07,  &
   &      2.234e-07,  2.267e-07,  2.287e-07,  2.310e-07,  2.322e-07,  &
   &      2.335e-07,  2.346e-07,  2.349e-07,  2.345e-07,  2.345e-07,  &
   &      2.338e-07,  2.328e-07,  2.313e-07,  2.290e-07,  2.262e-07,  &
   &      2.221e-07,  2.164e-07,  2.094e-07,  1.993e-07,  1.866e-07,  &
   &      1.680e-07,  1.478e-07,  1.301e-07,  1.210e-07,  1.176e-07,  &
   &      1.186e-07,  1.231e-07,  1.292e-07,  1.386e-07,  1.562e-07,  &
   &      1.758e-07,  1.934e-07,  2.154e-07,  2.415e-07,  2.582e-07,  &
   &      2.794e-07,  2.968e-07,  3.101e-07,  3.216e-07,  3.337e-07,  &
   &      3.480e-07,  3.600e-07,  3.730e-07,  3.840e-07,  3.930e-07,  &
   &      4.000e-07,  4.040e-07,  4.080e-07,  4.111e-07,  4.128e-07,  &
   &      4.127e-07,  4.116e-07,  4.088e-07,  4.021e-07,  3.931e-07,  &
   &      3.782e-07,  3.619e-07,  3.457e-07,  3.263e-07,  3.031e-07,  &
   &      2.806e-07,  2.633e-07,  2.461e-07,  2.288e-07,  2.147e-07,  &
   &      1.980e-07,  1.806e-07,  1.535e-07,  1.331e-07,  1.039e-07,  &
   &      6.677e-08,  4.973e-08,  4.163e-08,  3.696e-08,  3.430e-08,  &
   &      3.243e-08,  2.983e-08,  2.845e-08,  2.709e-08,  2.574e-08,  &
   &      2.438e-08,  2.342e-08,  2.273e-08,  2.228e-08,  2.197e-08,  &
   &      2.167e-08,  2.134e-08,  2.079e-08,  2.011e-08,  1.942e-08,  &
   &      1.938e-08,  1.933e-08,  1.927e-08,  1.922e-08,  1.917e-08,  &
   &      1.912e-08,  1.907e-08,  1.902e-08,  1.897e-08,  1.892e-08,  &
   &      1.887e-08,  1.882e-08,  1.877e-08,  1.872e-08,  1.867e-08,  &
   &      1.862e-08,  1.857e-08,  1.852e-08,  1.847e-08,  1.842e-08,  &
   &      1.837e-08,  1.832e-08,  1.827e-08,  1.822e-08,  1.817e-08,  &
   &      1.812e-08,  1.807e-08,  1.802e-08,  1.797e-08,  1.792e-08,  &
   &      1.787e-08,  1.782e-08,  1.777e-08,  1.772e-08,  1.767e-08,  &
   &      0.000e+00/

end block data bo2inf3

!     --------------------------------------------------------------
!
SUBROUTINE O2_vis (V1C,V2C,DVC,NPTC,C,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   DIMENSION C(*)
   COMMON /o2_o2_vis/ V1s,V2s,DVs,NPTs, s(1474)

   DATA XLOSMT / 2.68675E+19 /
!
!        O2 continuum formulated by Greenblatt et al. over the spectral
!        region 8797-29870 cm-1:  "Absorption Coefficients of Oxygen
!        Between 330 and 1140 nm, G.D. Greenblatt, J.J. Orlando, J.B.
!        Burkholder and A.R. Ravishabkara,  J. Geophys. Res., 95,
!        18577-18582, 1990.
!
!        The units conversion  is to (cm^2/molec)/atm(o2)
!
!      These are the conditions reported in the paper by Greenblatt et
!      al. for the spectrum of Fig. 1.
!
!     conditions:  55 atm.; 296 K; 89.5 cm path
!
   factor = 1./((xlosmt*1.e-20*(55.*273./296.)**2)*89.5)
!
   DVC = DVS
   v1ss = v1s
   v2ss = v2s
!
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      vj = v1c + dvc* REAL(j-1)
!
      C(J) = factor*S(I)/vj

10 END DO
!
   RETURN
!
end subroutine O2_vis
!
!     --------------------------------------------------------------
!
BLOCK DATA bo2in_vis

   IMPLICIT REAL*8 (V)

   COMMON /o2_o2_vis/ V1,V2,DV,NPT,                                  &
   &  o2vis0001(36),o2vis0051(50),o2vis0101(50),o2vis0151(50),        &
   &  o2vis0201(50),o2vis0251(50),o2vis0301(50),o2vis0351(50),        &
   &  o2vis0401(50),o2vis0451(50),o2vis0501(50),o2vis0551(50),        &
   &  o2vis0601(50),o2vis0651(50),o2vis0701(50),o2vis0751(50),        &
   &  o2vis0801(50),o2vis0851(50),o2vis0901(50),o2vis0951(50),        &
   &  o2vis1001(50),o2vis1051(50),o2vis1101(50),o2vis1151(50),        &
   &  o2vis1201(50),o2vis1251(50),o2vis1301(50),o2vis1351(50),        &
   &  o2vis1401(50),o2vis1451(38)

   DATA V1,V2,DV,NPT /15140.0, 29870.0, 10.0,  1474/

   DATA o2vis0001/                                                   &
   &      0.00E+00,                                                   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   6.06E-04,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.49E-03,   3.00E-03,   3.00E-03,   3.00E-03,   4.00E-03,   &
   &      4.00E-03,   5.00E-03,   5.00E-03,   6.00E-03,   7.00E-03/
   DATA o2vis0051/                                                   &
   &      8.00E-03,   9.00E-03,   1.00E-02,   1.10E-02,   1.25E-02,   &
   &      1.46E-02,   1.60E-02,   1.80E-02,   2.00E-02,   2.23E-02,   &
   &      2.50E-02,   2.69E-02,   3.00E-02,   3.30E-02,   3.63E-02,   &
   &      4.01E-02,   4.42E-02,   4.67E-02,   5.14E-02,   5.55E-02,   &
   &      5.96E-02,   6.43E-02,   6.94E-02,   7.37E-02,   7.88E-02,   &
   &      8.38E-02,   8.86E-02,   9.37E-02,   9.89E-02,   1.03E-01,   &
   &      1.07E-01,   1.10E-01,   1.14E-01,   1.16E-01,   1.18E-01,   &
   &      1.19E-01,   1.20E-01,   1.21E-01,   1.20E-01,   1.20E-01,   &
   &      1.19E-01,   1.17E-01,   1.16E-01,   1.13E-01,   1.10E-01,   &
   &      1.07E-01,   1.03E-01,   9.97E-02,   9.58E-02,   9.15E-02/
   DATA o2vis0101/                                                   &
   &      8.80E-02,   8.41E-02,   7.94E-02,   7.53E-02,   7.17E-02,   &
   &      6.83E-02,   6.43E-02,   6.08E-02,   5.69E-02,   5.31E-02,   &
   &      5.02E-02,   4.77E-02,   4.40E-02,   4.23E-02,   3.94E-02,   &
   &      3.70E-02,   3.51E-02,   3.30E-02,   3.10E-02,   2.90E-02,   &
   &      2.79E-02,   2.60E-02,   2.50E-02,   2.32E-02,   2.20E-02,   &
   &      2.10E-02,   2.00E-02,   1.90E-02,   1.80E-02,   1.70E-02,   &
   &      1.65E-02,   1.50E-02,   1.40E-02,   1.30E-02,   1.30E-02,   &
   &      1.20E-02,   1.10E-02,   1.10E-02,   1.00E-02,   1.00E-02,   &
   &      9.00E-03,   9.00E-03,   9.00E-03,   8.00E-03,   8.00E-03,   &
   &      7.01E-03,   7.00E-03,   7.00E-03,   6.98E-03,   6.00E-03/
   DATA o2vis0151/                                                   &
   &      5.80E-03,   5.00E-03,   5.00E-03,   5.00E-03,   4.00E-03,   &
   &      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   &
   &      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   3.00E-03,   &
   &      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   &
   &      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   &
   &      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   &
   &      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   &
   &      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   5.00E-03,   &
   &      5.00E-03,   6.00E-03,   6.00E-03,   7.00E-03,   7.41E-03,   &
   &      8.15E-03,   9.00E-03,   1.01E-02,   1.10E-02,   1.20E-02/
   DATA o2vis0201/                                                   &
   &      1.40E-02,   1.50E-02,   1.70E-02,   1.85E-02,   1.97E-02,   &
   &      2.24E-02,   2.47E-02,   2.74E-02,   3.06E-02,   3.36E-02,   &
   &      3.70E-02,   4.05E-02,   4.49E-02,   4.93E-02,   5.47E-02,   &
   &      6.01E-02,   6.52E-02,   7.23E-02,   7.89E-02,   8.80E-02,   &
   &      9.61E-02,   1.05E-01,   1.17E-01,   1.26E-01,   1.39E-01,   &
   &      1.49E-01,   1.60E-01,   1.68E-01,   1.74E-01,   1.79E-01,   &
   &      1.82E-01,   1.84E-01,   1.85E-01,   1.84E-01,   1.83E-01,   &
   &      1.81E-01,   1.80E-01,   1.77E-01,   1.74E-01,   1.71E-01,   &
   &      1.68E-01,   1.64E-01,   1.60E-01,   1.55E-01,   1.51E-01,   &
   &      1.46E-01,   1.40E-01,   1.36E-01,   1.30E-01,   1.25E-01/
   DATA o2vis0251/                                                   &
   &      1.20E-01,   1.14E-01,   1.09E-01,   1.05E-01,   9.93E-02,   &
   &      9.30E-02,   8.88E-02,   8.38E-02,   7.94E-02,   7.51E-02,   &
   &      7.08E-02,   6.66E-02,   6.32E-02,   6.01E-02,   5.55E-02,   &
   &      5.24E-02,   4.93E-02,   4.63E-02,   4.41E-02,   4.15E-02,   &
   &      3.90E-02,   3.63E-02,   3.50E-02,   3.26E-02,   3.05E-02,   &
   &      2.94E-02,   2.73E-02,   2.62E-02,   2.46E-02,   2.36E-02,   &
   &      2.25E-02,   2.10E-02,   2.00E-02,   1.90E-02,   1.80E-02,   &
   &      1.76E-02,   1.70E-02,   1.60E-02,   1.50E-02,   1.49E-02,   &
   &      1.40E-02,   1.30E-02,   1.30E-02,   1.22E-02,   1.20E-02,   &
   &      1.20E-02,   1.10E-02,   1.10E-02,   1.10E-02,   1.00E-02/
   DATA o2vis0301/                                                   &
   &      1.00E-02,   1.00E-02,   1.00E-02,   9.16E-03,   9.00E-03,   &
   &      9.00E-03,   9.00E-03,   9.00E-03,   8.49E-03,   8.00E-03,   &
   &      8.00E-03,   8.00E-03,   8.00E-03,   8.00E-03,   8.00E-03,   &
   &      8.00E-03,   7.00E-03,   8.00E-03,   7.00E-03,   7.00E-03,   &
   &      7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   &
   &      7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   &
   &      7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   &
   &      7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   &
   &      7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   &
   &      7.00E-03,   7.00E-03,   8.00E-03,   8.00E-03,   8.00E-03/
   DATA o2vis0351/                                                   &
   &      8.00E-03,   8.00E-03,   8.00E-03,   9.00E-03,   9.00E-03,   &
   &      9.00E-03,   9.07E-03,   1.00E-02,   1.00E-02,   1.00E-02,   &
   &      1.10E-02,   1.10E-02,   1.20E-02,   1.22E-02,   1.30E-02,   &
   &      1.31E-02,   1.40E-02,   1.50E-02,   1.60E-02,   1.70E-02,   &
   &      1.82E-02,   2.00E-02,   2.01E-02,   2.10E-02,   2.20E-02,   &
   &      2.28E-02,   2.30E-02,   2.30E-02,   2.30E-02,   2.30E-02,   &
   &      2.30E-02,   2.30E-02,   2.30E-02,   2.20E-02,   2.20E-02,   &
   &      2.20E-02,   2.10E-02,   2.10E-02,   2.00E-02,   2.00E-02,   &
   &      1.90E-02,   1.90E-02,   1.82E-02,   1.80E-02,   1.74E-02,   &
   &      1.70E-02,   1.63E-02,   1.60E-02,   1.50E-02,   1.49E-02/
   DATA o2vis0401/                                                   &
   &      1.40E-02,   1.37E-02,   1.30E-02,   1.30E-02,   1.21E-02,   &
   &      1.20E-02,   1.13E-02,   1.09E-02,   1.00E-02,   9.34E-03,   &
   &      9.00E-03,   8.43E-03,   8.00E-03,   7.39E-03,   7.00E-03,   &
   &      6.00E-03,   6.00E-03,   5.74E-03,   5.00E-03,   5.00E-03,   &
   &      5.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   &
   &      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   &
   &      3.17E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   &
   &      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   &
   &      3.00E-03,   3.00E-03,   3.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03/
   DATA o2vis0451/                                                   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      1.04E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03/
   DATA o2vis0501/                                                   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.41E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   1.98E-03,   1.46E-03,   &
   &      1.05E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03/
   DATA o2vis0551/                                                   &
   &      1.00E-03,   1.00E-03,   1.71E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   3.00E-03,   3.00E-03,   &
   &      3.82E-03,   4.00E-03,   4.17E-03,   5.00E-03,   6.00E-03,   &
   &      7.00E-03,   7.73E-03,   8.07E-03,   9.70E-03,   1.17E-02,   &
   &      1.31E-02,   1.47E-02,   1.64E-02,   1.81E-02,   2.07E-02,   &
   &      2.37E-02,   2.70E-02,   2.97E-02,   3.27E-02,   3.70E-02,   &
   &      4.13E-02,   4.49E-02,   4.89E-02,   5.38E-02,   5.98E-02,   &
   &      6.45E-02,   6.94E-02,   7.41E-02,   8.01E-02,   8.51E-02,   &
   &      9.00E-02,   9.49E-02,   9.88E-02,   1.01E-01,   1.04E-01,   &
   &      1.07E-01,   1.07E-01,   1.06E-01,   1.03E-01,   1.00E-01/
   DATA o2vis0601/                                                   &
   &      9.66E-02,   8.93E-02,   8.35E-02,   7.92E-02,   7.33E-02,   &
   &      6.84E-02,   6.40E-02,   5.91E-02,   5.57E-02,   5.26E-02,   &
   &      5.03E-02,   4.75E-02,   4.48E-02,   4.26E-02,   4.07E-02,   &
   &      3.83E-02,   3.69E-02,   3.47E-02,   3.24E-02,   3.11E-02,   &
   &      2.85E-02,   2.69E-02,   2.55E-02,   2.42E-02,   2.21E-02,   &
   &      2.09E-02,   1.93E-02,   1.77E-02,   1.62E-02,   1.60E-02,   &
   &      1.44E-02,   1.36E-02,   1.30E-02,   1.16E-02,   1.10E-02,   &
   &      1.00E-02,   1.00E-02,   9.00E-03,   8.27E-03,   8.00E-03,   &
   &      7.45E-03,   7.00E-03,   7.00E-03,   6.18E-03,   6.00E-03,   &
   &      6.00E-03,   5.00E-03,   5.00E-03,   5.00E-03,   5.00E-03/
   DATA o2vis0651/                                                   &
   &      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   &
   &      4.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   &
   &      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   2.07E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   1.28E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03/
   DATA o2vis0701/                                                   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   &
   &      4.00E-03,   4.00E-03,   4.00E-03,   4.57E-03,   5.00E-03,   &
   &      5.00E-03,   5.64E-03,   6.00E-03,   6.67E-03,   7.00E-03,   &
   &      7.35E-03,   8.00E-03,   8.36E-03,   9.00E-03,   9.00E-03,   &
   &      1.00E-02,   1.00E-02,   1.00E-02,   1.00E-02,   1.00E-02,   &
   &      1.00E-02,   1.00E-02,   9.65E-03,   9.00E-03,   9.00E-03,   &
   &      8.00E-03,   8.00E-03,   7.69E-03,   7.00E-03,   7.00E-03/
   DATA o2vis0751/                                                   &
   &      6.44E-03,   6.00E-03,   6.00E-03,   6.00E-03,   5.00E-03,   &
   &      5.00E-03,   5.00E-03,   5.00E-03,   5.00E-03,   4.00E-03,   &
   &      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   3.98E-03,   &
   &      3.01E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   &
   &      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   2.54E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03/
   DATA o2vis0801/                                                   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      1.33E-03,   1.89E-03,   1.07E-03,   1.06E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03/
   DATA o2vis0851/                                                   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03/
   DATA o2vis0901/                                                   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   5.50E-04,   &
   &      0.00E+00,   0.00E+00,   1.00E-03,   1.00E-03,   7.51E-04,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00/
   DATA o2vis0951/                                                   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   1.34E-04,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   7.65E-05,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03/
   DATA o2vis1001/                                                   &
   &      1.00E-03,   1.20E-04,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00/
   DATA o2vis1051/                                                   &
   &      0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   0.00E+00,   &
   &      0.00E+00,   6.09E-04,   3.47E-04,   6.97E-04,   2.60E-04,   &
   &      7.81E-04,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.68E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.76E-03,   3.00E-03,   3.00E-03,   3.00E-03/
   DATA o2vis1101/                                                   &
   &      3.80E-03,   4.00E-03,   4.82E-03,   5.00E-03,   5.84E-03,   &
   &      6.00E-03,   6.85E-03,   7.85E-03,   8.86E-03,   9.86E-03,   &
   &      1.09E-02,   1.19E-02,   1.29E-02,   1.47E-02,   1.59E-02,   &
   &      1.77E-02,   1.97E-02,   2.09E-02,   2.27E-02,   2.47E-02,   &
   &      2.67E-02,   2.87E-02,   3.07E-02,   3.26E-02,   3.38E-02,   &
   &      3.56E-02,   3.68E-02,   3.86E-02,   3.90E-02,   3.98E-02,   &
   &      4.07E-02,   4.10E-02,   4.10E-02,   4.03E-02,   3.93E-02,   &
   &      3.83E-02,   3.73E-02,   3.64E-02,   3.48E-02,   3.34E-02,   &
   &      3.18E-02,   2.99E-02,   2.85E-02,   2.70E-02,   2.50E-02,   &
   &      2.31E-02,   2.11E-02,   1.92E-02,   1.76E-02,   1.63E-02/
   DATA o2vis1151/                                                   &
   &      1.47E-02,   1.34E-02,   1.17E-02,   1.07E-02,   9.78E-03,   &
   &      8.81E-03,   7.84E-03,   6.88E-03,   6.00E-03,   5.94E-03,   &
   &      5.00E-03,   5.00E-03,   4.05E-03,   4.00E-03,   3.13E-03,   &
   &      3.00E-03,   3.00E-03,   2.24E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.00E-03,   1.54E-03,   &
   &      1.41E-03,   1.64E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03/
   DATA o2vis1201/                                                   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   1.00E-03,   &
   &      1.00E-03,   1.15E-03,   2.00E-03,   2.00E-03,   2.00E-03,   &
   &      2.00E-03,   2.00E-03,   2.00E-03,   2.56E-03,   3.00E-03,   &
   &      3.00E-03,   3.30E-03,   4.00E-03,   4.00E-03,   4.04E-03,   &
   &      4.95E-03,   5.85E-03,   6.00E-03,   6.67E-03,   7.58E-03,   &
   &      8.48E-03,   9.39E-03,   1.03E-02,   1.14E-02,   1.31E-02/
   DATA o2vis1251/                                                   &
   &      1.40E-02,   1.58E-02,   1.76E-02,   1.94E-02,   2.12E-02,   &
   &      2.30E-02,   2.56E-02,   2.89E-02,   3.16E-02,   3.44E-02,   &
   &      3.80E-02,   4.16E-02,   4.52E-02,   4.87E-02,   5.23E-02,   &
   &      5.59E-02,   5.91E-02,   6.20E-02,   6.53E-02,   6.71E-02,   &
   &      6.89E-02,   6.98E-02,   7.07E-02,   7.10E-02,   7.10E-02,   &
   &      7.06E-02,   6.97E-02,   6.89E-02,   6.80E-02,   6.71E-02,   &
   &      6.54E-02,   6.43E-02,   6.29E-02,   6.11E-02,   5.94E-02,   &
   &      5.74E-02,   5.48E-02,   5.31E-02,   5.05E-02,   4.86E-02,   &
   &      4.62E-02,   4.41E-02,   4.23E-02,   4.03E-02,   3.78E-02,   &
   &      3.61E-02,   3.43E-02,   3.26E-02,   3.08E-02,   2.91E-02/
   DATA o2vis1301/                                                   &
   &      2.73E-02,   2.58E-02,   2.49E-02,   2.31E-02,   2.22E-02,   &
   &      2.07E-02,   1.95E-02,   1.86E-02,   1.77E-02,   1.69E-02,   &
   &      1.60E-02,   1.51E-02,   1.43E-02,   1.40E-02,   1.35E-02,   &
   &      1.27E-02,   1.18E-02,   1.10E-02,   1.10E-02,   1.02E-02,   &
   &      1.00E-02,   1.00E-02,   9.67E-03,   8.81E-03,   8.05E-03,   &
   &      8.90E-03,   8.24E-03,   8.00E-03,   7.53E-03,   7.00E-03,   &
   &      7.00E-03,   7.00E-03,   7.00E-03,   7.00E-03,   6.42E-03,   &
   &      6.00E-03,   6.00E-03,   6.00E-03,   6.00E-03,   5.18E-03,   &
   &      5.00E-03,   5.00E-03,   5.00E-03,   4.80E-03,   4.04E-03,   &
   &      4.89E-03,   4.27E-03,   4.00E-03,   4.00E-03,   4.00E-03/
   DATA o2vis1351/                                                   &
   &      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   &
   &      4.00E-03,   4.00E-03,   4.00E-03,   3.20E-03,   3.00E-03,   &
   &      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   &
   &      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   &
   &      3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   3.00E-03,   &
   &      3.00E-03,   3.75E-03,   4.00E-03,   4.00E-03,   4.00E-03,   &
   &      4.00E-03,   4.00E-03,   4.69E-03,   5.00E-03,   5.00E-03,   &
   &      5.15E-03,   5.97E-03,   6.00E-03,   6.61E-03,   7.43E-03,   &
   &      8.00E-03,   8.06E-03,   8.88E-03,   9.70E-03,   1.05E-02,   &
   &      1.13E-02,   1.21E-02,   1.30E-02,   1.38E-02,   1.52E-02/
   DATA o2vis1401/                                                   &
   &      1.64E-02,   1.72E-02,   1.80E-02,   1.88E-02,   1.96E-02,   &
   &      2.04E-02,   2.10E-02,   2.10E-02,   2.10E-02,   2.10E-02,   &
   &      2.10E-02,   2.10E-02,   2.10E-02,   2.10E-02,   2.10E-02,   &
   &      2.05E-02,   2.00E-02,   1.99E-02,   1.91E-02,   1.90E-02,   &
   &      1.85E-02,   1.80E-02,   1.79E-02,   1.71E-02,   1.63E-02,   &
   &      1.55E-02,   1.47E-02,   1.40E-02,   1.40E-02,   1.33E-02,   &
   &      1.25E-02,   1.20E-02,   1.19E-02,   1.11E-02,   1.03E-02,   &
   &      1.00E-02,   9.75E-03,   9.00E-03,   9.00E-03,   8.37E-03,   &
   &      8.00E-03,   8.00E-03,   8.00E-03,   7.22E-03,   7.00E-03,   &
   &      7.00E-03,   6.86E-03,   6.07E-03,   6.00E-03,   6.00E-03/
   DATA o2vis1451/                                                   &
   &      6.00E-03,   5.93E-03,   5.15E-03,   5.00E-03,   5.00E-03,   &
   &      5.00E-03,   5.00E-03,   5.00E-03,   5.00E-03,   5.00E-03,   &
   &      5.00E-03,   5.00E-03,   5.00E-03,   5.00E-03,   5.00E-03,   &
   &      5.00E-03,   5.00E-03,   5.00E-03,   4.68E-03,   4.00E-03,   &
   &      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   &
   &      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   &
   &      4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   4.00E-03,   &
   &      1.00E-03,   2.00E-04,   0./
!
end block data bo2in_vis

!     --------------------------------------------------------------
!
SUBROUTINE O2HERZ (V1C,V2C,DVC,NPTC,C,T,P,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   DIMENSION C(*)
!
   V1S = 36000.
   v1ss = V1S
   v2ss = 99999.

   DVS = 10.
   DVC = DVS
!
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 0.01
   end if
!
   V1C = V1S + DVS*REAL(I1-1)
   I2 = (V2C-V1S)/DVS + 0.01
   NPTC = I2-I1+3
!         IF (NPTC.GT.NPTS) NPTC=NPTS+4
!        mja, 10-27-2011 - this seems to be redundant as the
!        Herzberg O2 continuum is a function, not block data
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF (I.LT.1) GO TO 10
      VJ = V1C+DVC* REAL(J-1)
      CALL HERTDA (HERZ,VJ)
      CALL HERPRS (HERZ,T,P)
      C(J) = HERZ/VJ
10 END DO
!
   RETURN
!
end subroutine O2HERZ
!
!     --------------------------------------------------------------
!
SUBROUTINE HERTDA (HERZ,V)
!
   IMPLICIT REAL*8           (V)
!
!     HERZBERG O2 ABSORPTION
!     HALL,1987 PRIVATE COMMUNICATION, BASED ON:
!
!     REF. JOHNSTON, ET AL., JGR,89,11661-11665,1984
!          NICOLET, 1987 (RECENT STUDIES IN ATOMIC
!                         & MOLECULAR PROCESSES,
!                         PLENUM PUBLISHING CORP, NY 1987)
!
!     AND YOSHINO, ET AL., 1988 (PREPRINT OF "IMPROVED ABSORPTION
!          CROSS SECTIONS OF OXYGEN IN THE WAVELENGTH REGION 205-240NM
!          OF THE HERZBERG CONTINUUM")
!
!     **** NOTE:  CROSS SECTION AT 0 PRESSURE  ***
!     THE PRESSURE DEPENDENT TERM IS IN SUBROUTINE HERPRS
!
!C    COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP
!
   HERZ = 0.0
   IF (V.LE.36000.00) RETURN
!
!     EXTRAPOLATE SMOOTHLY THROUGH THE HERZBERG BAND REGION
!     NOTE: HERZBERG BANDS ARE NOT CORRECTLY INCLUDED
!
   CORR = 0.
!      IF (V.LE.40000.) CORR = ((40000.-V)/4000.)*7.917E-27
!     factor of 1.e-20 removed; put in front factor
   IF (V.LE.40000.) CORR = ((40000.-V)/4000.)*7.917E-07
!
!     UNITS ARE (CM2)
!
!     HALL'S NEW HERZBERG  (LEAST SQRS FIT, LN(P))
!
!     YRATIO=2048.7/WL(I)  ****IN ANGSTOMS****
!           =.20487/WN(I)     IN MICRONS
!           =WCM(I)/48811.0   IN CM-1
!
   YRATIO = V/48811.0
!     HERZ = 6.884E-24*(YRATIO)*EXP(-69.738*( LOG(YRATIO))**2)-CORR
!     factor of 1.e-20 removed; put in front factor
   HERZ = 6.884E-04*(YRATIO)*EXP(-69.738*( LOG(YRATIO))**2)-CORR
!
   RETURN
!
end subroutine HERTDA
!
!     --------------------------------------------------------------
!
SUBROUTINE HERPRS (HERZ,T,P)
!
!     CORRECT THE HERZBERG CONTINUUM CROSS SECTION FOR PRESSURE
!     DEPENDENCE; BASED ON SHARDANAND, JQRST, 18, 525-530, 1977.
!                 FOR UN2| BROADENING
!                 AND YOSHINO ET AL 1988 FOR UO2| BROADENING
!
!     PO2= PARTIAL PRESSURE OF O2
!     PN2= PARTIAL PRESSURE OF N2; BN2=.45*BO2
!
!     DATA BO2 / 1.72E-3 /
!
!     Changed in Herzberg continuum pressure,
!     Reference:
!     "Atmospheric Propagation in the UV, Visible, IR and MM-wave
!     Region and Related Systems Aspects".
!     G.P. Anderson,F.X. Kneizys, E.P. Shettle, L.W. Abreu,
!     J.H. Chetwynd, R.E. Huffman, and L.A. Hall; Conference
!     Proceedings No. 454 of the Advisory Group for Aerospace
!     Research & Development; 1990.
!
   DATA BO2 / 1.81E-3 /
   DATA PO / 1013. /,TO / 273.16 /
!
!     NOTE:  THE HERZBERG CONTINUUM OBEYS BEER'S LAW
!            OPTICAL DEPTH(TOTAL)=SUM OVER LAYER O.D.(I)
!
!     BO2= RATIO OF SIGMA(O2-O2)/(SIGMA(O2)) * 760(TORR)*.2095
!     BN2=.45*BO2= RATIO OF SIGMA(O2-N2)/(SIGMA(O2)) * 760(TORR)*.78
!
!     BO2*760*(.2095+.45*.78) = .73 , AS BELOW
!
!     Changed Herzberg continuum pressure (see above reference)
!
!     BO2*760*(.2095+.45*.78) = .83 , AS BELOW
!
!
   HERZ = HERZ*(1.+.83*(P/PO)*(TO/T))
!
   RETURN
!
end subroutine HERPRS
!
!     --------------------------------------------------------------
!
SUBROUTINE O2FUV (V1C,V2C,DVC,NPTC,C,v1ss,v2ss)
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)
   DIMENSION C(*)
   COMMON /o2_fuv/ V1s,V2s,DVs,NPTs, s(1512)
!
!     O2 continuum in the far-UV is from Lu et al. (2010)
!
   DVC = DVS
   v1ss = v1s
   v2ss = v2s
!
   V1C = V1ABS-DVC
   V2C = V2ABS+DVC
!
   IF (V1C.LT.V1S) then
      I1 = -1
   else
      I1 = (V1C-V1S)/DVS + 1.e-5
   end if
   V1C = V1S + DVS*REAL(I1-1)

   I2 = (V2C-V1S)/DVS + 1.e-5

   NPTC = I2-I1+3
   IF (NPTC.GT.NPTS) NPTC=NPTS+4
   V2C = V1C + DVS*REAL(NPTC-1)
!
   DO 10 J = 1, NPTC
      I = I1+(J-1)
      C(J) = 0.
      IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10
      vj = v1c + dvc* REAL(j-1)
      !
      C(J) = S(I)/vj
10 END DO
!
   RETURN
!
end subroutine O2FUV
!
!     --------------------------------------------------------------
!
BLOCK DATA bo2in_fuv

   IMPLICIT REAL*8 (V)
!     Only valid up to 87000 cm-1.

   COMMON /o2_fuv/ V1,V2,DV,NPT,                                     &
   &  o2fuv0001(13),o2fuv0051(50),o2fuv0101(50),o2fuv0151(50),        &
   &  o2fuv0201(50),o2fuv0251(50),o2fuv0301(50),o2fuv0351(50),        &
   &  o2fuv0401(50),o2fuv0451(50),o2fuv0501(50),o2fuv0551(50),        &
   &  o2fuv0601(50),o2fuv0651(50),o2fuv0701(50),o2fuv0751(50),        &
   &  o2fuv0801(50),o2fuv0851(50),o2fuv0901(50),o2fuv0951(50),        &
   &  o2fuv1001(50),o2fuv1051(50),o2fuv1101(50),o2fuv1151(50),        &
   &  o2fuv1201(50),o2fuv1251(50),o2fuv1301(50),o2fuv1351(50),        &
   &  o2fuv1401(50),o2fuv1451(50),o2fuv1501(49)

   DATA V1,V2,DV,NPT /56740.0, 86960.0, 20.0,  1512/

   DATA o2fuv0001/                                                   &
   &      0.00E+00,   1.50E-01,   6.00E-01,                           &
   &      2.81E+00,   8.44E+00,   1.05E+01,   1.89E+01,   1.93E+01,   &
   &      1.57E+01,   2.49E+01,   2.97E+01,   2.29E+01,   2.33E+01/
   DATA o2fuv0051/                                                   &
   &      3.18E+01,   3.42E+01,   3.86E+01,   3.90E+01,   4.62E+01,   &
   &      4.71E+01,   4.68E+01,   4.82E+01,   4.89E+01,   4.97E+01,   &
   &      5.05E+01,   5.18E+01,   5.15E+01,   5.27E+01,   5.40E+01,   &
   &      5.47E+01,   5.48E+01,   5.60E+01,   5.81E+01,   5.86E+01,   &
   &      5.95E+01,   5.94E+01,   6.18E+01,   6.21E+01,   6.28E+01,   &
   &      6.38E+01,   6.53E+01,   6.62E+01,   6.68E+01,   6.74E+01,   &
   &      6.88E+01,   6.90E+01,   7.13E+01,   7.25E+01,   7.25E+01,   &
   &      7.41E+01,   7.55E+01,   7.69E+01,   7.76E+01,   7.75E+01,   &
   &      7.96E+01,   7.64E+01,   7.78E+01,   7.98E+01,   8.05E+01,   &
   &      7.95E+01,   8.91E+01,   7.85E+01,   8.26E+01,   8.71E+01/
   DATA o2fuv0101/                                                   &
   &      8.64E+01,   9.15E+01,   8.38E+01,   9.18E+01,   8.85E+01,   &
   &      9.26E+01,   8.80E+01,   8.57E+01,   9.37E+01,   8.91E+01,   &
   &      9.19E+01,   1.01E+02,   9.63E+01,   1.01E+02,   1.00E+02,   &
   &      9.85E+01,   1.01E+02,   1.01E+02,   1.05E+02,   1.06E+02,   &
   &      1.08E+02,   1.09E+02,   1.06E+02,   1.12E+02,   1.07E+02,   &
   &      1.11E+02,   1.13E+02,   1.10E+02,   1.06E+02,   1.09E+02,   &
   &      1.12E+02,   1.17E+02,   1.19E+02,   1.14E+02,   1.13E+02,   &
   &      1.15E+02,   1.18E+02,   1.18E+02,   1.24E+02,   1.25E+02,   &
   &      1.18E+02,   1.28E+02,   1.27E+02,   1.25E+02,   1.26E+02,   &
   &      1.33E+02,   1.35E+02,   1.33E+02,   1.28E+02,   1.27E+02/
   DATA o2fuv0151/                                                   &
   &      1.37E+02,   1.39E+02,   1.41E+02,   1.39E+02,   1.42E+02,   &
   &      1.44E+02,   1.40E+02,   1.49E+02,   1.49E+02,   1.56E+02,   &
   &      1.58E+02,   1.54E+02,   1.51E+02,   1.55E+02,   1.53E+02,   &
   &      1.58E+02,   1.59E+02,   1.59E+02,   1.64E+02,   1.58E+02,   &
   &      1.59E+02,   1.67E+02,   1.68E+02,   1.63E+02,   1.68E+02,   &
   &      1.76E+02,   1.67E+02,   1.79E+02,   1.82E+02,   1.78E+02,   &
   &      1.75E+02,   1.75E+02,   1.79E+02,   1.87E+02,   1.85E+02,   &
   &      1.79E+02,   1.91E+02,   1.86E+02,   1.94E+02,   1.95E+02,   &
   &      1.94E+02,   1.92E+02,   1.98E+02,   1.99E+02,   2.10E+02,   &
   &      2.08E+02,   2.01E+02,   2.07E+02,   2.09E+02,   2.14E+02/
   DATA o2fuv0201/                                                   &
   &      2.10E+02,   2.10E+02,   2.18E+02,   2.21E+02,   2.20E+02,   &
   &      2.16E+02,   2.23E+02,   2.22E+02,   2.29E+02,   2.25E+02,   &
   &      2.33E+02,   2.36E+02,   2.29E+02,   2.29E+02,   2.37E+02,   &
   &      2.41E+02,   2.44E+02,   2.44E+02,   2.47E+02,   2.48E+02,   &
   &      2.47E+02,   2.50E+02,   2.57E+02,   2.50E+02,   2.60E+02,   &
   &      2.56E+02,   2.58E+02,   2.56E+02,   2.65E+02,   2.69E+02,   &
   &      2.68E+02,   2.73E+02,   2.72E+02,   2.76E+02,   2.74E+02,   &
   &      2.78E+02,   2.82E+02,   2.84E+02,   2.89E+02,   2.93E+02,   &
   &      2.90E+02,   2.92E+02,   2.96E+02,   2.96E+02,   3.00E+02,   &
   &      2.96E+02,   3.00E+02,   3.01E+02,   3.15E+02,   3.09E+02/
   DATA o2fuv0251/                                                   &
   &      3.15E+02,   3.13E+02,   3.16E+02,   3.25E+02,   3.22E+02,   &
   &      3.18E+02,   3.21E+02,   3.27E+02,   3.29E+02,   3.33E+02,   &
   &      3.27E+02,   3.30E+02,   3.44E+02,   3.47E+02,   3.45E+02,   &
   &      3.45E+02,   3.44E+02,   3.45E+02,   3.52E+02,   3.50E+02,   &
   &      3.58E+02,   3.60E+02,   3.61E+02,   3.68E+02,   3.70E+02,   &
   &      3.69E+02,   3.70E+02,   3.73E+02,   3.78E+02,   3.79E+02,   &
   &      3.84E+02,   3.86E+02,   3.87E+02,   3.90E+02,   3.91E+02,   &
   &      3.91E+02,   3.95E+02,   4.05E+02,   4.03E+02,   4.01E+02,   &
   &      4.07E+02,   4.15E+02,   4.11E+02,   4.24E+02,   4.24E+02,   &
   &      4.23E+02,   4.32E+02,   4.24E+02,   4.31E+02,   4.35E+02/
   DATA o2fuv0301/                                                   &
   &      4.34E+02,   4.41E+02,   4.41E+02,   4.48E+02,   4.48E+02,   &
   &      4.47E+02,   4.48E+02,   4.55E+02,   4.58E+02,   4.66E+02,   &
   &      4.68E+02,   4.66E+02,   4.73E+02,   4.80E+02,   4.76E+02,   &
   &      4.77E+02,   4.82E+02,   4.92E+02,   4.89E+02,   5.00E+02,   &
   &      4.85E+02,   5.02E+02,   4.89E+02,   5.09E+02,   5.13E+02,   &
   &      5.15E+02,   5.25E+02,   5.24E+02,   5.23E+02,   5.23E+02,   &
   &      5.22E+02,   5.28E+02,   5.41E+02,   5.38E+02,   5.38E+02,   &
   &      5.39E+02,   5.45E+02,   5.51E+02,   5.55E+02,   5.58E+02,   &
   &      5.57E+02,   5.67E+02,   5.68E+02,   5.64E+02,   5.75E+02,   &
   &      5.74E+02,   5.76E+02,   5.81E+02,   5.84E+02,   5.91E+02/
   DATA o2fuv0351/                                                   &
   &      5.97E+02,   5.95E+02,   6.00E+02,   5.99E+02,   6.08E+02,   &
   &      6.11E+02,   6.15E+02,   6.15E+02,   6.13E+02,   6.21E+02,   &
   &      6.24E+02,   6.23E+02,   6.33E+02,   6.35E+02,   6.39E+02,   &
   &      6.46E+02,   6.50E+02,   6.47E+02,   6.45E+02,   6.51E+02,   &
   &      6.56E+02,   6.64E+02,   6.77E+02,   6.71E+02,   6.74E+02,   &
   &      6.82E+02,   6.85E+02,   6.81E+02,   6.81E+02,   6.88E+02,   &
   &      6.93E+02,   7.00E+02,   7.01E+02,   7.00E+02,   7.07E+02,   &
   &      7.17E+02,   7.17E+02,   7.19E+02,   7.23E+02,   7.21E+02,   &
   &      7.30E+02,   7.40E+02,   7.47E+02,   7.43E+02,   7.48E+02,   &
   &      7.50E+02,   7.47E+02,   7.53E+02,   7.70E+02,   7.56E+02/
   DATA o2fuv0401/                                                   &
   &      7.68E+02,   7.68E+02,   7.74E+02,   7.71E+02,   7.86E+02,   &
   &      7.92E+02,   7.84E+02,   7.93E+02,   7.99E+02,   8.05E+02,   &
   &      8.09E+02,   8.12E+02,   8.04E+02,   8.12E+02,   8.14E+02,   &
   &      8.14E+02,   8.19E+02,   8.26E+02,   8.31E+02,   8.34E+02,   &
   &      8.41E+02,   8.45E+02,   8.48E+02,   8.48E+02,   8.53E+02,   &
   &      8.57E+02,   8.54E+02,   8.58E+02,   8.62E+02,   8.71E+02,   &
   &      8.76E+02,   8.90E+02,   8.94E+02,   8.90E+02,   8.92E+02,   &
   &      8.98E+02,   9.06E+02,   9.07E+02,   9.10E+02,   9.15E+02,   &
   &      9.20E+02,   9.21E+02,   9.19E+02,   9.19E+02,   9.28E+02,   &
   &      9.37E+02,   9.40E+02,   9.38E+02,   9.48E+02,   9.51E+02/
   DATA o2fuv0451/                                                   &
   &      9.52E+02,   9.58E+02,   9.49E+02,   9.67E+02,   9.66E+02,   &
   &      9.72E+02,   9.78E+02,   9.85E+02,   9.81E+02,   9.78E+02,   &
   &      9.86E+02,   9.90E+02,   9.96E+02,   1.01E+03,   1.00E+03,   &
   &      1.00E+03,   1.01E+03,   1.01E+03,   1.02E+03,   1.01E+03,   &
   &      1.02E+03,   1.03E+03,   1.04E+03,   1.04E+03,   1.04E+03,   &
   &      1.04E+03,   1.04E+03,   1.05E+03,   1.06E+03,   1.05E+03,   &
   &      1.06E+03,   1.06E+03,   1.06E+03,   1.07E+03,   1.07E+03,   &
   &      1.08E+03,   1.08E+03,   1.09E+03,   1.09E+03,   1.09E+03,   &
   &      1.09E+03,   1.10E+03,   1.11E+03,   1.11E+03,   1.10E+03,   &
   &      1.11E+03,   1.11E+03,   1.13E+03,   1.12E+03,   1.12E+03/
   DATA o2fuv0501/                                                   &
   &      1.13E+03,   1.13E+03,   1.13E+03,   1.13E+03,   1.14E+03,   &
   &      1.14E+03,   1.14E+03,   1.15E+03,   1.16E+03,   1.16E+03,   &
   &      1.18E+03,   1.17E+03,   1.16E+03,   1.17E+03,   1.18E+03,   &
   &      1.18E+03,   1.18E+03,   1.18E+03,   1.20E+03,   1.19E+03,   &
   &      1.19E+03,   1.21E+03,   1.20E+03,   1.21E+03,   1.21E+03,   &
   &      1.21E+03,   1.22E+03,   1.22E+03,   1.22E+03,   1.24E+03,   &
   &      1.23E+03,   1.23E+03,   1.23E+03,   1.24E+03,   1.25E+03,   &
   &      1.26E+03,   1.26E+03,   1.25E+03,   1.25E+03,   1.26E+03,   &
   &      1.26E+03,   1.27E+03,   1.27E+03,   1.28E+03,   1.28E+03,   &
   &      1.28E+03,   1.27E+03,   1.28E+03,   1.29E+03,   1.29E+03/
   DATA o2fuv0551/                                                   &
   &      1.30E+03,   1.30E+03,   1.29E+03,   1.30E+03,   1.30E+03,   &
   &      1.31E+03,   1.31E+03,   1.30E+03,   1.32E+03,   1.31E+03,   &
   &      1.32E+03,   1.32E+03,   1.33E+03,   1.33E+03,   1.33E+03,   &
   &      1.33E+03,   1.33E+03,   1.34E+03,   1.34E+03,   1.36E+03,   &
   &      1.35E+03,   1.34E+03,   1.34E+03,   1.35E+03,   1.36E+03,   &
   &      1.36E+03,   1.36E+03,   1.35E+03,   1.36E+03,   1.35E+03,   &
   &      1.35E+03,   1.38E+03,   1.40E+03,   1.39E+03,   1.39E+03,   &
   &      1.38E+03,   1.39E+03,   1.39E+03,   1.39E+03,   1.39E+03,   &
   &      1.40E+03,   1.40E+03,   1.41E+03,   1.41E+03,   1.41E+03,   &
   &      1.41E+03,   1.40E+03,   1.41E+03,   1.40E+03,   1.42E+03/
   DATA o2fuv0601/                                                   &
   &      1.42E+03,   1.43E+03,   1.43E+03,   1.43E+03,   1.42E+03,   &
   &      1.43E+03,   1.44E+03,   1.42E+03,   1.43E+03,   1.43E+03,   &
   &      1.44E+03,   1.44E+03,   1.43E+03,   1.45E+03,   1.45E+03,   &
   &      1.44E+03,   1.45E+03,   1.45E+03,   1.45E+03,   1.46E+03,   &
   &      1.47E+03,   1.49E+03,   1.46E+03,   1.47E+03,   1.48E+03,   &
   &      1.47E+03,   1.46E+03,   1.47E+03,   1.48E+03,   1.48E+03,   &
   &      1.47E+03,   1.47E+03,   1.47E+03,   1.48E+03,   1.49E+03,   &
   &      1.48E+03,   1.49E+03,   1.48E+03,   1.47E+03,   1.48E+03,   &
   &      1.50E+03,   1.50E+03,   1.50E+03,   1.50E+03,   1.50E+03,   &
   &      1.51E+03,   1.50E+03,   1.50E+03,   1.50E+03,   1.50E+03/
   DATA o2fuv0651/                                                   &
   &      1.50E+03,   1.50E+03,   1.49E+03,   1.51E+03,   1.52E+03,   &
   &      1.51E+03,   1.51E+03,   1.52E+03,   1.51E+03,   1.51E+03,   &
   &      1.52E+03,   1.53E+03,   1.53E+03,   1.52E+03,   1.53E+03,   &
   &      1.52E+03,   1.52E+03,   1.51E+03,   1.51E+03,   1.52E+03,   &
   &      1.53E+03,   1.53E+03,   1.53E+03,   1.52E+03,   1.52E+03,   &
   &      1.53E+03,   1.54E+03,   1.54E+03,   1.52E+03,   1.52E+03,   &
   &      1.52E+03,   1.52E+03,   1.54E+03,   1.54E+03,   1.53E+03,   &
   &      1.52E+03,   1.52E+03,   1.54E+03,   1.56E+03,   1.55E+03,   &
   &      1.55E+03,   1.55E+03,   1.54E+03,   1.55E+03,   1.55E+03,   &
   &      1.54E+03,   1.52E+03,   1.53E+03,   1.55E+03,   1.54E+03/
   DATA o2fuv0701/                                                   &
   &      1.55E+03,   1.54E+03,   1.54E+03,   1.53E+03,   1.52E+03,   &
   &      1.55E+03,   1.54E+03,   1.56E+03,   1.56E+03,   1.55E+03,   &
   &      1.57E+03,   1.55E+03,   1.55E+03,   1.55E+03,   1.55E+03,   &
   &      1.57E+03,   1.54E+03,   1.55E+03,   1.56E+03,   1.54E+03,   &
   &      1.53E+03,   1.54E+03,   1.55E+03,   1.55E+03,   1.55E+03,   &
   &      1.56E+03,   1.56E+03,   1.56E+03,   1.55E+03,   1.55E+03,   &
   &      1.55E+03,   1.55E+03,   1.54E+03,   1.54E+03,   1.54E+03,   &
   &      1.53E+03,   1.54E+03,   1.53E+03,   1.54E+03,   1.55E+03,   &
   &      1.55E+03,   1.55E+03,   1.54E+03,   1.54E+03,   1.54E+03,   &
   &      1.54E+03,   1.54E+03,   1.54E+03,   1.54E+03,   1.54E+03/
   DATA o2fuv0751/                                                   &
   &      1.53E+03,   1.54E+03,   1.53E+03,   1.53E+03,   1.53E+03,   &
   &      1.52E+03,   1.54E+03,   1.53E+03,   1.53E+03,   1.53E+03,   &
   &      1.52E+03,   1.54E+03,   1.53E+03,   1.53E+03,   1.54E+03,   &
   &      1.55E+03,   1.54E+03,   1.54E+03,   1.54E+03,   1.53E+03,   &
   &      1.52E+03,   1.53E+03,   1.54E+03,   1.53E+03,   1.51E+03,   &
   &      1.51E+03,   1.52E+03,   1.52E+03,   1.52E+03,   1.53E+03,   &
   &      1.51E+03,   1.50E+03,   1.50E+03,   1.52E+03,   1.51E+03,   &
   &      1.50E+03,   1.49E+03,   1.50E+03,   1.52E+03,   1.51E+03,   &
   &      1.52E+03,   1.52E+03,   1.48E+03,   1.50E+03,   1.50E+03,   &
   &      1.50E+03,   1.50E+03,   1.48E+03,   1.49E+03,   1.48E+03/
   DATA o2fuv0801/                                                   &
   &      1.47E+03,   1.48E+03,   1.48E+03,   1.50E+03,   1.50E+03,   &
   &      1.48E+03,   1.45E+03,   1.47E+03,   1.47E+03,   1.46E+03,   &
   &      1.47E+03,   1.47E+03,   1.46E+03,   1.45E+03,   1.45E+03,   &
   &      1.45E+03,   1.46E+03,   1.45E+03,   1.44E+03,   1.45E+03,   &
   &      1.44E+03,   1.43E+03,   1.42E+03,   1.42E+03,   1.43E+03,   &
   &      1.43E+03,   1.43E+03,   1.40E+03,   1.41E+03,   1.40E+03,   &
   &      1.38E+03,   1.39E+03,   1.38E+03,   1.38E+03,   1.37E+03,   &
   &      1.36E+03,   1.37E+03,   1.35E+03,   1.33E+03,   1.34E+03,   &
   &      1.35E+03,   1.33E+03,   1.32E+03,   1.32E+03,   1.29E+03,   &
   &      1.28E+03,   1.27E+03,   1.25E+03,   1.25E+03,   1.24E+03/
   DATA o2fuv0851/                                                   &
   &      1.23E+03,   1.21E+03,   1.19E+03,   1.18E+03,   1.17E+03,   &
   &      1.15E+03,   1.13E+03,   1.12E+03,   1.11E+03,   1.07E+03,   &
   &      1.06E+03,   1.04E+03,   1.02E+03,   1.01E+03,   9.90E+02,   &
   &      9.68E+02,   9.44E+02,   9.14E+02,   9.00E+02,   8.85E+02,   &
   &      8.67E+02,   8.46E+02,   8.28E+02,   8.16E+02,   8.07E+02,   &
   &      8.00E+02,   7.83E+02,   7.75E+02,   7.75E+02,   7.70E+02,   &
   &      7.60E+02,   7.55E+02,   7.54E+02,   7.28E+02,   7.44E+02,   &
   &      7.59E+02,   7.53E+02,   7.41E+02,   7.51E+02,   7.64E+02,   &
   &      7.63E+02,   7.65E+02,   7.67E+02,   7.67E+02,   7.58E+02,   &
   &      7.38E+02,   7.20E+02,   7.05E+02,   7.01E+02,   6.71E+02/
   DATA o2fuv0901/                                                   &
   &      6.49E+02,   6.48E+02,   6.25E+02,   5.98E+02,   5.71E+02,   &
   &      5.45E+02,   5.06E+02,   4.84E+02,   4.73E+02,   4.45E+02,   &
   &      4.24E+02,   4.00E+02,   3.69E+02,   3.53E+02,   3.39E+02,   &
   &      3.27E+02,   3.17E+02,   2.98E+02,   2.82E+02,   2.72E+02,   &
   &      2.69E+02,   2.58E+02,   2.51E+02,   2.52E+02,   2.48E+02,   &
   &      2.42E+02,   2.35E+02,   2.31E+02,   2.38E+02,   2.41E+02,   &
   &      2.35E+02,   2.28E+02,   2.28E+02,   2.30E+02,   2.31E+02,   &
   &      2.37E+02,   2.43E+02,   2.47E+02,   2.50E+02,   2.48E+02,   &
   &      2.44E+02,   2.41E+02,   2.37E+02,   2.33E+02,   2.32E+02,   &
   &      2.30E+02,   2.30E+02,   2.34E+02,   2.47E+02,   2.48E+02/
   DATA o2fuv0951/                                                   &
   &      2.45E+02,   2.39E+02,   2.48E+02,   2.53E+02,   2.50E+02,   &
   &      2.28E+02,   2.31E+02,   2.39E+02,   2.45E+02,   2.40E+02,   &
   &      2.33E+02,   2.25E+02,   2.19E+02,   2.18E+02,   2.16E+02,   &
   &      2.09E+02,   2.06E+02,   2.07E+02,   2.10E+02,   2.04E+02,   &
   &      2.03E+02,   2.03E+02,   1.98E+02,   1.91E+02,   1.84E+02,   &
   &      1.78E+02,   1.68E+02,   1.65E+02,   1.67E+02,   1.62E+02,   &
   &      1.57E+02,   1.51E+02,   1.47E+02,   1.41E+02,   1.36E+02,   &
   &      1.30E+02,   1.26E+02,   1.23E+02,   1.21E+02,   1.13E+02,   &
   &      1.09E+02,   1.06E+02,   1.01E+02,   9.54E+01,   9.08E+01,   &
   &      8.81E+01,   8.50E+01,   8.12E+01,   7.68E+01,   7.57E+01/
   DATA o2fuv1001/                                                   &
   &      7.21E+01,   6.72E+01,   6.38E+01,   6.10E+01,   5.76E+01,   &
   &      5.19E+01,   5.02E+01,   4.84E+01,   4.53E+01,   4.06E+01,   &
   &      4.77E+01,   5.24E+01,   3.21E+01,   2.91E+01,   3.15E+01,   &
   &      4.32E+01,   3.42E+01,   2.19E+01,   3.81E+01,   3.69E+01,   &
   &      2.94E+01,   3.39E+01,   4.05E+01,   2.71E+01,   3.47E+01,   &
   &      3.13E+01,   3.00E+01,   2.37E+01,   3.80E+01,   2.79E+01,   &
   &      3.18E+01,   5.54E+01,   4.15E+01,   4.95E+01,   6.05E+01,   &
   &      2.71E+01,   4.48E+01,   3.24E+01,   4.10E+01,   3.09E+01,   &
   &      5.75E+01,   4.17E+01,   4.27E+01,   3.34E+01,   4.28E+01,   &
   &      4.43E+01,   6.25E+01,   5.29E+01,   4.58E+01,   5.00E+01/
   DATA o2fuv1051/                                                   &
   &      5.35E+01,   5.64E+01,   7.03E+01,   5.49E+01,   5.02E+01,   &
   &      4.99E+01,   5.43E+01,   5.11E+01,   6.67E+01,   4.96E+01,   &
   &      3.67E+01,   7.64E+01,   5.68E+01,   5.43E+01,   6.72E+01,   &
   &      5.98E+01,   7.08E+01,   7.38E+01,   8.53E+01,   5.81E+01,   &
   &      5.27E+01,   6.31E+01,   5.17E+01,   4.52E+01,   5.86E+01,   &
   &      6.72E+01,   5.70E+01,   3.92E+01,   6.05E+01,   4.40E+01,   &
   &      7.16E+01,   5.53E+01,   6.79E+01,   5.38E+01,   4.91E+01,   &
   &      4.50E+01,   4.47E+01,   5.20E+01,   6.93E+01,   4.93E+01,   &
   &      4.83E+01,   5.13E+01,   3.73E+01,   4.88E+01,   4.62E+01,   &
   &      4.17E+01,   3.50E+01,   3.12E+01,   3.65E+01,   3.89E+01/
   DATA o2fuv1101/                                                   &
   &      5.29E+01,   4.43E+01,   5.46E+01,   4.33E+01,   3.36E+01,   &
   &      3.41E+01,   3.88E+01,   2.80E+01,   2.87E+01,   2.81E+01,   &
   &      2.77E+01,   3.40E+01,   2.41E+01,   3.97E+01,   2.70E+01,   &
   &      1.94E+01,   1.94E+01,   1.07E+01,   4.82E+01,   3.65E+01,   &
   &      1.86E+01,   2.14E+01,   4.18E+01,   1.57E+01,   3.23E+01,   &
   &      1.25E+01,   1.45E+01,   2.26E+01,   1.22E+01,   1.83E+01,   &
   &      1.84E+01,   1.60E+01,   6.62E+00,   9.56E+00,   2.12E+00,   &
   &      2.24E+00,   3.89E+00,   1.04E+01,   1.00E+01,   1.58E+01,   &
   &      2.70E+01,   2.25E+01,   1.05E+01,   2.22E+01,   1.93E+01,   &
   &      1.39E+00,   2.08E+01,   3.38E+01,   1.24E+01,   3.37E+01/
   DATA o2fuv1151/                                                   &
   &      7.77E+00,   8.82E+00,   2.21E+01,   7.28E+00,   2.12E+01,   &
   &      3.22E+01,   1.18E+01,   2.04E+01,   4.08E+01,   4.33E+01,   &
   &      4.11E+01,   4.78E+01,   5.46E+01,   3.50E+01,   5.61E+01,   &
   &      3.45E+01,   2.66E+01,   4.45E+01,   4.21E+01,   4.36E+01,   &
   &      3.78E+01,   5.85E+01,   3.82E+01,   4.76E+01,   5.26E+01,   &
   &      6.05E+01,   6.20E+01,   7.09E+01,   6.54E+01,   6.20E+01,   &
   &      8.90E+01,   7.33E+01,   6.40E+01,   4.73E+01,   8.59E+01,   &
   &      7.83E+01,   8.50E+01,   8.20E+01,   9.01E+01,   9.17E+01,   &
   &      8.37E+01,   9.41E+01,   1.01E+02,   9.00E+01,   1.11E+02,   &
   &      1.05E+02,   1.00E+02,   1.32E+02,   1.51E+02,   1.56E+02/
   DATA o2fuv1201/                                                   &
   &      1.76E+02,   1.65E+02,   1.35E+02,   1.56E+02,   1.80E+02,   &
   &      1.57E+02,   1.95E+02,   1.66E+02,   1.84E+02,   1.90E+02,   &
   &      2.31E+02,   2.52E+02,   2.98E+02,   3.92E+02,   5.95E+02,   &
   &      9.08E+02,   1.42E+03,   3.61E+03,   5.76E+03,   3.74E+03,   &
   &      4.55E+03,   4.36E+03,   2.24E+03,   1.01E+03,   7.27E+02,   &
   &      2.36E+02,   8.03E+01,   5.40E+01,   5.39E+01,   3.16E+01,   &
   &      2.22E+01,   3.46E+01,   2.98E+01,   5.38E+01,   2.48E+01,   &
   &      4.09E+01,   3.63E+01,   3.31E+01,   6.51E+01,   5.23E+01,   &
   &      4.38E+01,   5.59E+01,   7.41E+01,   3.86E+01,   5.56E+01,   &
   &      5.50E+01,   3.96E+01,   6.23E+01,   5.22E+01,   7.83E+01/
   DATA o2fuv1251/                                                   &
   &      6.92E+01,   4.53E+01,   4.24E+01,   4.84E+01,   1.02E+02,   &
   &      5.88E+01,   7.17E+01,   5.04E+01,   4.49E+01,   1.41E+01,   &
   &      5.84E+01,   4.82E+01,   6.24E+01,   2.86E+01,   5.23E+01,   &
   &      3.14E+01,   3.55E+01,   3.80E+01,   1.53E+01,   5.39E+01,   &
   &      4.29E+01,   2.64E+01,   2.14E+01,   5.24E+01,   5.56E+01,   &
   &      3.31E+01,   5.20E+01,   3.90E+01,   3.50E+01,   4.66E-01,   &
   &      2.82E+01,   9.34E+00,   3.27E+01,   1.87E+01,   4.05E+01,   &
   &      1.37E+01,   1.91E+01,   3.45E+01,   1.73E+01,   2.99E+01,   &
   &      3.11E+01,   4.81E+01,   2.16E+01,   3.75E+01,   3.03E+01,   &
   &      2.94E+01,   2.07E+01,   3.12E+01,   1.89E+01,   8.00E+00/
   DATA o2fuv1301/                                                   &
   &      4.80E-01,   2.03E+01,   1.21E+01,   4.53E+01,   1.46E+01,   &
   &      3.63E+01,   6.59E+00,   1.58E+01,   1.63E+01,   1.30E+01,   &
   &      1.14E+00,   1.76E+01,   2.24E+01,   2.29E+01,   2.01E+01,   &
   &      4.47E+00,   3.82E+01,   6.67E+00,   1.33E+01,   2.35E+01,   &
   &      2.82E+01,   2.68E+01,   3.97E+01,   4.01E+01,   4.64E+01,   &
   &      6.94E+01,   9.20E+01,   9.07E+01,   1.56E+02,   1.17E+02,   &
   &      1.34E+02,   2.02E+02,   2.20E+02,   3.23E+02,   4.28E+02,   &
   &      3.73E+02,   5.95E+02,   7.50E+02,   9.79E+02,   1.21E+03,   &
   &      1.26E+03,   1.32E+03,   1.30E+03,   1.44E+03,   1.41E+03,   &
   &      1.53E+03,   1.73E+03,   1.81E+03,   1.84E+03,   1.70E+03/
   DATA o2fuv1351/                                                   &
   &      1.53E+03,   1.30E+03,   1.15E+03,   9.53E+02,   8.98E+02,   &
   &      6.83E+02,   7.02E+02,   5.10E+02,   4.85E+02,   4.98E+02,   &
   &      3.62E+02,   4.05E+02,   2.97E+02,   2.60E+02,   2.95E+02,   &
   &      2.54E+02,   2.49E+02,   2.19E+02,   1.89E+02,   1.33E+02,   &
   &      1.92E+02,   1.60E+02,   1.47E+02,   1.33E+02,   1.12E+02,   &
   &      8.22E+01,   1.25E+02,   7.80E+01,   1.09E+02,   1.00E+02,   &
   &      1.13E+02,   8.59E+01,   8.27E+01,   4.05E+01,   7.16E+01,   &
   &      5.32E+01,   7.80E+01,   8.18E+01,   4.63E+01,   3.39E+01,   &
   &      3.88E+01,   1.78E+01,   4.35E+01,   3.87E+01,   2.20E+00,   &
   &      1.15E+01,   6.25E+00,   7.12E+01,   2.73E+01,   2.58E+01/
   DATA o2fuv1401/                                                   &
   &      9.90E+00,   7.16E+00,   1.57E+01,   1.28E+01,   7.76E+00,   &
   &      2.79E+00,   5.46E+00,   7.73E+00,   2.65E+01,   7.39E+00,   &
   &      4.60E+00,   1.15E+00,   6.40E+00,   7.02E+00,   5.32E+00,   &
   &      9.56E+00,   3.76E+00,   1.45E+01,   8.59E+00,   1.30E+01,   &
   &      3.54E+00,   9.19E-01,   7.70E+00,   2.95E+01,   1.08E+01,   &
   &      2.27E+00,   1.06E+01,   6.45E+00,   1.74E+01,   1.28E+01,   &
   &      1.82E+01,   4.31E+01,   2.13E+01,   2.55E+01,   2.26E+01,   &
   &      1.75E+01,   1.83E+01,   3.11E+01,   1.62E+01,   3.31E+01,   &
   &      3.33E+01,   2.44E+01,   2.23E+01,   2.86E+01,   2.70E+01,   &
   &      2.83E+01,   4.04E+01,   5.29E+01,   6.31E+01,   6.27E+01/
   DATA o2fuv1451/                                                   &
   &      6.03E+01,   7.34E+01,   8.07E+01,   8.87E+01,   8.87E+01,   &
   &      1.02E+02,   1.36E+02,   1.38E+02,   1.47E+02,   1.50E+02,   &
   &      1.63E+02,   1.65E+02,   1.48E+02,   1.40E+02,   1.12E+02,   &
   &      1.20E+02,   1.13E+02,   8.94E+01,   5.64E+01,   5.32E+01,   &
   &      3.46E+01,   3.41E+01,   3.26E+01,   1.12E+01,   1.28E+01,   &
   &      1.04E+00,   1.64E+01,   2.49E+00,   1.33E+01,   8.67E+00,   &
   &      9.16E+00,   1.24E-01,   5.54E+00,   4.67E-01,   4.39E+00,   &
   &      4.54E+00,   1.69E+00,   7.87E+00,   1.42E+01,   9.46E+00,   &
   &      1.41E+01,   3.61E+01,   9.26E+01,   1.78E+02,   1.47E+02,   &
   &      2.36E+02,   1.03E+02,   1.44E+02,   1.86E+02,   1.65E+02/
   DATA o2fuv1501/                                                   &
   &      1.89E+02,   1.47E+02,   1.18E+02,   3.07E+02,   1.11E+02,   &
   &      7.62E+01,   2.89E+01,   1.57E+01,   1.69E-01,   2.56E+00,   &
   &      5.14E+00,   2.00E+00,   3.91E-01,   8.80E+00,   2.83E+00,   &
   &      3.39E+00,   2.22E+00,   6.03E+00,   5.87E+00,   1.14E+01,   &
   &      4.12E-01,   4.53E+00,   1.42E+01,   7.51E+00,   1.00E-01,   &
   &      2.96E-01,   6.34E+00,   8.74E+00,   1.33E+01,   2.80E+01,   &
   &      1.40E+02,   3.14E+02,   4.33E+02,   3.61E+02,   3.65E+02,   &
   &      3.48E+02,   4.25E+02,   2.92E+02,   2.33E+02,   2.07E+02,   &
   &      3.27E+02,   4.09E+02,   3.90E+02,   2.57E+02,   1.47E+02,   &
   &      1.09E+02,   9.49E+01,   5.31E+01,   0.00E+00/

end block data bo2in_fuv



!     --------------------------------------------------------------
subroutine cld_od(V1C,V2C,DVC,NPTC,C,layer,xkt)
!     --------------------------------------------------------------
!
   Use lblparams, ONLY: n_absrb
   IMPLICIT REAL*8           (V)
!
   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)

   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
   &              NLTEFL,LNFIL4,LNGTH4
!
   parameter (n_lyr=200,n_cld=500)
!
   COMMON /cld_rd/ n_freq, n_align, v_cloud_freq(n_cld),             &
   &                                      cloudodlayer(n_lyr,n_cld)
   DIMENSION C(*)
!
   logical EX
   character*55 in_cld_file
   dimension i_layer(n_lyr), pres_layer_dum(n_lyr), v_cntnm(n_lyr)
!
   data in_cld_file /'in_lblrtm_cld'/
   data dvs /5./
!
!     ----------------------------------------------------------
!     Read in TES cloud effective optical depth file
!     ----------------------------------------------------------
!
   if (layer .eq. 1) then

      open (35,FILE=in_cld_file,STATUS='OLD')
      read (35,*) n_freq
      read (35,*) (v_cloud_freq(j),j=1,n_freq)

      write (ipr,*)
      write (ipr,*)                                                  &
      &           '** iaersl=5; Cloud Information from "in_cld_file" **'
      write (ipr,'(" n_freq = ",i5)') n_freq
      write (ipr,'(5x,10f10.4)') (v_cloud_freq(j),j=1,n_freq)

      read (35,*) n_layer

      write (ipr,'(" n_layer = ",i5)') n_layer

      do l =1,n_layer
         read (35,*) i_layer(l), pres_layer_dum(l)
         read (35,*) (cloudodlayer(l,j),j=1,n_freq)
         write (ipr,'(i5,f12.5)') i_layer(l), pres_layer_dum(l)
         write (ipr,'(5x,10f10.4)') (cloudodlayer(l,j),j=1,n_freq)
      enddo
      close (35)

   endif

!
!        ----------------------------------------------------------
!        Generated output continuum grid
!        ----------------------------------------------------------
!
   DVC = DVS
   V1C = V1ABS-10.
   V2C = V2ABS+10.
   NPTC = ((V2C - V1C)/DVC) + 1

   do J = 1, NPTC
      v_cntnm(j) = V1C+DVC* REAL(J-1)
   enddo
!
!        ----------------------------------------------------------
!        Linearly interpolate TES cloud effective OD onto continuum grid
!        ----------------------------------------------------------

   ilo = 1
   do j=1, NPTC
      IF (v_cntnm(j).LE.v_cloud_freq(1))      THEN
         C(j) = cloudodlayer(layer,1)
         GO TO 10
      ELSE IF (v_cntnm(j).GT.v_cloud_freq(n_freq)) THEN
         C(j) = cloudodlayer(layer,n_freq)
         GO TO 10
      END IF

      do i=ilo, n_freq
         IF (v_cntnm(j).LE.v_cloud_freq(i))  THEN
            v_m=(cloudodlayer(layer,i)-cloudodlayer(layer,i-1))/ (         &
               v_cloud_freq(i)-v_cloud_freq(i-1))
            C(j) = cloudodlayer(layer,i-1)+ (v_cntnm(j)-v_cloud_freq(i-1))*&
               v_m
            ilo = i-1
            GO TO 10
         END IF
      enddo

10    continue
!
      if (v_cntnm(j).eq.0.) then
         C(j) = 0.
      else
         C(j) = C(j)/RADFN(v_cntnm(j),XKT)
      endif

      if (ilo.lt.1) ilo = 1
!
   enddo

!
   RETURN
end subroutine cld_od
