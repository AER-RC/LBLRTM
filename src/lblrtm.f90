!     path:      $HeadURL$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!
!  --------------------------------------------------------------------------
! |  Copyright �, Atmospheric and Environmental Research, Inc., 2015        |
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
      PROGRAM LBLRTM 
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!                                                                       
!                                                                       
!                  IMPLEMENTATION:    R.D. WORSHAM                      
!                                                                       
!             ALGORITHM REVISIONS:    S.A. CLOUGH                       
!                                     R.D. WORSHAM                      
!                                     J.L. MONCET                       
!                                     M.W. SHEPHARD  
!                                     D. WEISENSTEIN                   
!                                                                       
!                                                                       
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.       
!                     131 Hartwell Ave, Lexington, MA, 02421            
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!               WORK SUPPORTED BY:    THE ARM PROGRAM                   
!                                     OFFICE OF ENERGY RESEARCH         
!                                     DEPARTMENT OF ENERGY  
!                                     
!                                     AND
!    
!                                     THE JOINT CENTER FOR
!                                      SATELLITE DATA ASSIMILATION            
!                                                                       
!                                                                       
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL           
!                                                                       
!                                             FASCOD3                   
!                                                                       
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!                                                                       
!                                                                       
!********************************************************************** 
!*                                                                      
!*                           LBLRTM                                     
!*                                                                      
!*                FAST ATMOSPHERIC SIGNATURE CODE                       
!*                                                                      
!*                                                                      
!*                AIR FORCE GEOPHYSICS LABORATORY                       
!*                    HANSCOM AFB, MA 01731                             
!*                                                                      
!*                                                                      
!*                      SHEPARD A. CLOUGH                               
!*                      FRANCIS X. KNEIZYS                              
!*                      GAIL P. ANDERSON                                
!*                      JAMES H. CHETWYND JR.                           
!*                      ROBERT D. WORSHAM                               
!*                      ERIC P. SHETTLE                                 
!*                      LEONARD W. ABREU                                
!*                                                                      
!********************************************************************** 
!*                                                                      
!*    DOCUMENTATION AND INFORMATION ABOUT THE PROGRAM  MAY BE OBTAINED  
!*    FROM                                                              
!*                                                                      
!*    AER Inc., 131 Hartwell Ave., Lexington, MA 02421 USA 
!*    <aer_lblrtm@aer.com>
!*                                                                      
!********************************************************************** 
!*                                                                      
!*                  LINE PARAMETER COMPILATIONS                         
!*                        
!*    THE HITRAN 2008 MOLECULAR SPECTROSCOPIC DATABASE
!*       L. S. ROTHMAN ET AL.
!*       J. QUANT. SPECTROSC. RADIAT. TRANSFER, VOL 110, P533-572 (2009) 
!*
!*    THE HITRAN 2004 MOLECULAR SPECTROSCOPIC DATABASE
!*       L. S. ROTHMAN ET AL.
!*       J. QUANT. SPECTROSC. RADIAT. TRANSFER, VOL 96, P139-204 (2005) 
!*
!*    THE HITRAN MOLECULAR DATABASE: EDITION OF 2000 INCLUDING UPDATES 
!*    THROUGH 2001
!*       L. S. ROTHMAN ET AL.
!*       J. QUANT. SPECTROSC. RADIAT. TRANSFER (2003) 
!*       
!*    THE HITRAN MOLECULAR DATABASE: 1996 EDITION
!*       L. S. ROTHMAN ET AL.
!*       J. QUANT. SPECTROSC. RADIAT. TRANSFER, VOL 60 (1998)        
!*
!*   THE HITRAN DATABASE: 1986 EDITION                                  
!*      L.S. ROTHMAN, R.R. GAMACHE, A. GOLDMAN, L. R. BROWN,            
!*      R. A. TOTH, H. M. PICKETT, R. L. POYNTER, J.-M. FLAUD,          
!*      C. CAMY-PEYRET, A. BARBE, N. HUSSON, C. P. RINSLAND,            
!*      AND M. A. H. SMITH                                              
!*      APPLIED OPTICS, VOL. 26,  P 4058(OCT 1987)                      
!*                                                                      
!*   AFGL ATMOSPHERIC ABSORPTION LINE PARAMETERS COMPILATION:           
!*   1982 EDITION                                                       
!*      L. S. ROTHMAN, R. R. GAMACHE, A. BARBE, A. GOLDMAN,             
!*      J. R. GILLIS, L. R. BROWN, R. A. TOTH, J.-M. FLAUD AND          
!*      C. CAMY-PEYRET                                                  
!*      APPLIED OPTICS, VOL. 22, P 2247(AUG 1983)                       
!*                                                                      
!*   AFGL TRACE GAS COMPILATION: 1982 VERSION                           
!*      L. S. ROTHMAN, A. GOLDMAN, J. R. GILLIS, R. R. GAMACHE,         
!*      H.M. PICKETT, R. L. POYNTER, N. HUSSON AND A. CHEDIN            
!*      APPLIED OPTICS VOL. 22, P 1616(JUN 1983)                        
!*                                                                      
!*   AFCRL ATMOSPHERIC ABSORPTION LINE PARAMETERS COMPILATION           
!*                                                     AFCRL-TR-73-0096 
!*      R. A. MCCLATCHEY, W. S. BENEDICT, S. A. CLOUGH, D. E. BURCH,    
!*      R. F. CALFEE, K. FOX, L.S. ROTHMAN AND J. S. GARING             
!*                                                                      
!********************************************************************** 
!*                                                                      
!*       LBLATM - ATMOSPHERIC OPTICAL PROPERTIES                        
!*                SPHERICAL REFRACTIVE GEOMETRY                         
!*                                                                      
!*    AFGL ATMOSPHERIC CONSTITUENT PROFILES (0-120 KM)  AFGL-TR-86-0110 
!*       G. P. ANDERSON, S. A. CLOUGH, F.X. KNEIZYS, J. H. CHETWYND     
!*       AND E. P. SHETTLE                                              
!*                                                                      
!*    AIR MASS COMPUTER PROGRAM FOR ATMOSPHERIC TRANSMITTANCE/RADIANCE: 
!*    FSCATM                                            AFGL-TR-83-0065 
!*       W. O. GALLERY, F. X. KNEIZYS, AND S. A. CLOUGH                 
!*                                                                      
!********************************************************************** 
!*                                                                      
!*      LOWTRN - AEROSOLS, HYDROMETEORS AND MOLECULAR SCATTERING        
!*                                                                      
!*    ATMOSPHERIC TRANSMITTANCE/RADIANCE:                               
!*    COMPUTER CODE LOWTRAN 6                           AFGL-TR-83-0187 
!*       F. X. KNEIZYS, E. P. SHETTLE, W. O. GALLERY,                   
!*       J. H. CHETWYND,JR, L. W. ABREU, J. E. A. SELBY,                
!*       S. A. CLOUGH AND R. W. FENN                                    
!*                                                                      
!*    ATMOSPHERIC TRANSMITTANCE AND RADIANCE: THE LOWTRAN 5 CODE        
!*       F. X. KNEIZYS, E. P. SHETTLE, AND W. O. GALLERY                
!*       SPIE, V277 (1981) P 116                                        
!*                                                                      
!*    ATMOSPHERIC TRANSMITTANCE/RADIANCE:                               
!*    COMPUTER CODE LOWTRAN 5                          AFGL-TR-80-00676 
!*       F. X. KNEIZYS, E. P. SHETTLE, W. O. GALLERY,                   
!*       J. H. CHETWYND,JR, L. W. ABREU, J. E. A. SELBY,                
!*       R. W. FENN AND R. A. MCCLATCHEY                                
!*                                                                      
!*    ATMOSPHERIC ATTENUATION OF MILLIMETER AND SUBMILLIMETER WAVES:    
!*    MODELS AND COMPUTER CODES                         AFGL-TR-79-0253 
!*       V. J. FALCONE, JR., L. W. ABREU AND E. P. SHETTLE              
!*                                                                      
!*    .    .    .    .    .    .    .    .    .    .    .    .    .     
!*                                                                      
!*    MODELS OF THE ATMOSPHERIC AEROSOLS AND THEIR OPTICAL PROPERTIES,  
!*    IN AGARD PROC. 183, OPTICAL PROPAGATION IN THE ATMOSPHERE (1976)  
!*       E. P. SHETTLE AND R. W. FENN                                   
!*                                                                      
!*    MODELS OF THE AEROSOLS OF THE LOWER ATMOSPHERE AND THE EFFECTS OF 
!*    HUMIDITY VARIATIONS ON THEIR OPTICAL PROPERTIES   AFGL-TR-79-0214 
!*       E. P. SHETTLE AND R. W. FENN                                   
!*                                                                      
!********************************************************************** 
!*                                                                      
!*           NLTE - NON LOCAL THERMODYNAMIC EQUILIBRIUM                 
!*                                                                      
!*                                                                      
!*    ATMOSPHERIC TRANSMITTANCE/RADIANCE: COMPUTER CODE FASCOD2         
!*                                                                      
!*        W. L. RIDGWAY, R. A. MOOSE, AND A. C. COGLEY                  
!*                                                      AFGL-TR-82-0392 
!*                                                      SONICRAFT, INC. 
!*                                                                      
!*    .    .    .   .    .    .    .    .    .    .    .    .    .      
!*                                                                      
!*                        NLTE REFERENCES                               
!*                                                                      
!*   A USER'S GUIDE TO THE AFGL/VISIDYNE HIGH ALTITUDE INFRARED         
!*   RADIANCE MODEL COMPUTER PROGRAM                    AFGL-TR-85-0015 
!*      T. C. DEGGES AND A. P. D'AGATI                  VISIDYNE/AFGL   
!*                                                                      
!*   A HIGH ALTITUDE INFRARED RADIANCE MODEL            AFGL-TR-77-0271 
!*      T. C. DEGGES AND H. J. P. SMITH                 VISIDYNE,INC    
!*                                                                      
!********************************************************************** 
!*                                                                      
!*         Analytic Jacobians (IEMIT=3 and IMRG=40-43)                  
!*                                                                      
!*  May 2004:  This work was funded by Eumetsat (Stephen Tjemkes)       
!*             and is based upon work done for the NASA-EOS-TES project 
!*                                                                      
!********************************************************************** 
!*                                                                      
!*                                                                      
!*               GENERAL LBLRTM  REFERENCES -                           
!*
!*    Shephard, M.W., S.A. Clough, V.H. Payne, W. L. Smith, S. Kireev, 
!*     and K. E. Cady-Pereira, Performance of the line-by-line radiative
!*     transfer model (LBLRTM) for temperature and species retrievals: 
!*     IASI case studies from JAIVEx, Atmos. Chem. Phys. Discuss., 9, 
!*     9313-9366, 2009.
!*
!*    Delamere, J. S., S.A. Clough, V. H. Payne, E. J. Mlawer, D. D.Turner 
!*     and R. R. Gamache, A far-infrared radiative closure study in the 
!*     Arctic: Application to water vapor,  J. Geophys. Res., 
!*     doi:10.1029/2009JD012968, 2010.
!*    
!*    Clough, S. A., M. W. Shephard, E. J. Mlawer, J. S. Delamere, 
!*     M. J. Iacono, K. Cady-Pereira, S. Boukabara, and P. D. Brown, 
!*     Atmospheric radiative transfer modeling: a summary of the AER codes, 
!*     Short Communication, J. Quant. Spectrosc. Radiat. Transfer, 91, 
!*     233-244, 2005.
!*                                                                  
!*    Clough, S.A., and M.J. Iacono, Line-by-line calculations of       
!*      atmospheric fluxes and cooling rates II: Application to carbon  
!*      dioxide, ozone, methane, nitrous oxide, and the halocarbons. J. 
!*      Geophys. Res., 100, 16,519-16,535, 1995.                        
!*                                                                      
!*    Clough, S.A., M.J. Iacono, and J.-L. Moncet, Line-by-line         
!*      calculation of atmospheric fluxes and cooling rates:  Applicatio
!*      to water vapor. J. Geophys. Res., 97, 15761-15785, 1992.        
!*                                                                      
!*    ATMOSPHERIC RADIANCE AND TRANSMITTANCE: FASCOD2                   
!*        S. A. CLOUGH, F. X. KNEIZYS, E. P. SHETTLE AND G. P. ANDERSON 
!*        PROC. OF THE SIXTH CONFERENCE ON ATMOSPHERIC RADIATION,       
!*        WILLIAMSBURG, VA (1986), P 141                                
!*                                                                      
!*    LINEAR ABSORPTION AND SCATTERING OF LASER BEAMS   AFGL-TR-84-0265 
!*        F. X. KNEIZYS, S. A. CLOUGH, E. P. SHETTLE, L. S. ROTHMAN     
!*        AND R. W. FENN                                                
!*                                                                      
!*                                                                      
!*    ATMOSPHERIC ATTENUATION OF LASER RADIATION                        
!*        F. X. KNEIZYS, S. A. CLOUGH, E. P. SHETTLE                    
!*        PROC. OF SPIE VOL. 410, LASER BEAM PROPAGATION(1983) P 13     
!*                                                                      
!*                                                                      
!*    ATMOSPHERIC SPECTRAL TRANSMITTANCE AND RADIANCE-FASCOD1B          
!*        S. A. CLOUGH, F. X. KNEIZYS, L. S. ROTHMAN AND W. O. GALLERY  
!*        PROC. OF SPIE VOL.277 ATMOSPERIC TRANSMISSION(1981) P152      
!*                                                                      
!*    Clough, S.A., F.X. Kneizys, R. Davis, R. Gamache and R. Tipping   
!*      (1980): Theoretical line shape for H2O vapor:  Application to th
!*      continuum.  Atmospheric Water Vapor, edited by A. Deepak, T.D.  
!*      Wilkerson and L.H. Ruhnke, 52,  Academic Press, New York.       
!*                                                                      
!*    CONVOLUTION ALGORITHM FOR THE LORENTZ FUNCTION                    
!*        S. A. CLOUGH AND F. X. KNEIZYS, APPLIED OPTICS 18, 2329(1979) 
!*                                                                      
!*                                                                      
!*    FASCODE - FAST ATMOSPHERIC SIGNATURE CODE                         
!*              (SPECTRAL TRANSMITTANCE AND RADIANCE)   AFGL-TR-78-0081 
!*        H. J. P. SMITH, D. J. DUBE, M. E. GARDNER,                    
!*        S. A. CLOUGH, F. X. KNEIZYS, AND L. S. ROTHMAN                
!*                                                                      
!*    ALGORITHM FOR THE CALCULATION OF ABSORPTION COEFFICIENT           
!*          - PRESSURE BROADENED MOLECULAR TRANSITIONS  AFGL-TR-77-0164 
!*        S. A. CLOUGH, F. X. KNEIZYS, J. H. CHETWYND, JR.              
!*                                                                      
!********************************************************************** 
!---------------------------------------------------------------------- 
!-                                                                      
!-                 FILE ASSIGNMENTS FOR LBLRTM                          
!-                                                                      
!-                                                                      
!-    TAPE3       UNFORMATTED LINE FILE WITH LBLRTM  BLOCKING           
!-                                  EXTERNAL FILE NOT REQUIRED FOR      
!-                                              IHIRAC=0                
!-                                              IHIRAC=9                
!-                                              IHIRAC NE 0,9 ; ITEST=1 
!-                                                                      
!-    TAPE4       NLTE VIBRATIONAL TEMPERATURES (POPULATIONS) BY LAYER  
!-                                        ONLY REQUIRED FOR IHIRAC=4    
!-                                                                      
!-    TAPE5       LBLRTM  INPUT FILE                                    
!-                                                                      
!-    TAPE6       LBLRTM  OUTPUT FILE                                   
!-                                                                      
!-    TAPE7       FILE OF MOLECULAR COLUMN AMOUNTS FROM LBLATM          
!-                               ONLY FOR IATM=1; IPUNCH=1 (CARD 2.1)   
!-                                                                      
!-    TAPE9       FILE OF EFFECTIVE LINES FOR LBLF4 CREATED BY LINF4    
!-                                                                      
!-    TAPE10      OPTICAL DEPTH RESULTS FROM LINE BY LINE CALCULATION   
!-                            LAST LAYER     FOR IMRG EQ 0              
!-                            LAYER BY LAYER FOR IMRG EQ 1              
!-                                                                      
!-    TAPE11      SPECTRAL RESULTS FROM SCANFN AND INTRPL               
!-                            JEMIT=-1: ABSORPTION                      
!-                            JEMIT= 0: TRANSMITTANCE                   
!-                            JEMIT= 1: RADIANCE                        
!-                                                                      
!-    TAPE12      MONOCHROMATIC RESULTS                                 
!-                            IEMIT=0: OPTICAL DEPTH                    
!-                            IEMIT=1: RADIANCE/TRANSMITTANCE           
!-                                 INCLUDES AEROSOL CONTRIBUTION FOR    
!-                                                 IAERSL=1; IEMIT=1    
!-                                                                      
!-    TAPE13      MONOCHROMATIC RESULTS FOR WEIGHTING FUNCTIONS         
!-                            IEMIT=0: OPTICAL DEPTH                    
!-                            IEMIT=1: RADIANCE/TRANSMITTANCE           
!-                                 ONLY CREATED FOR IMRG= 3 TO 18       
!-                                                                      
!-    TAPE14      MONCHROMATIC RESULTS INLUDING AEROSOL CONTRIBUTION    
!-                            IEMIT=0: OPTICAL DEPTH                    
!-                            IEMIT=1: RADIANCE/TRANSMITTANCE           
!-                                 ONLY CREATED FOR IAERSL=1; IEMIT=0   
!-                                                                      
!-    TAPE20      SPECTRAL AEROSOL TRANSMITTANCES                       
!-                            TOTAL AEROSOL  CONTRIBUTION FOR IEMIT=0   
!-                            LAYER BY LAYER CONTRIBUTION FOR IEMIT=1   
!-                                 ONLY CREATED FOR IAERSL=1            
!-                                                                      
!-    TAPE29      FILE CONTAINING VALUES OF Y FOR PLOTTING              
!-                                 ONLY FOR IPLOT EQ 1                  
!-                                                                      
!-                                                                      
!-    TAPE39      AFGL PLOT FILE                                        
!-                                                                      
!---------------------------------------------------------------------- 
!-                                                                      
!-                      STATEMENT FLAGS                                 
!-                                                                      
!-    LBLRTM  HAS BEEN STRUCTURED TO HAVE ENHANCED PORTABILITY UNDER    
!-    FORTRAN 77.  FOUR FLAGS (COLUMN73) HAVE BEEN USED TO FACILITATE   
!-    PROGRAM CONVERSION.                                               
!-                                                                      
!-   &    IDENTIFIES STATEMENTS REQUIRED FOR WORD SIZE LESS THAN 8      
!-               CHAR. ALL STATEMENTS FLAGGED WITH & IN COLUMN 73 HAVE  
!-            C& STARTING IN COLUMN 1. THESE TWO CHARACTERS MUST        
!-               BE CHANGED TO BLANKS FOR COMPUTERS WITH WORD SIZE      
!-               LESS THAN 8 CHARACTERS.                                
!-                                                                      
!-   !    IDENTIFIES STATEMENTS REQUIRED TO DOUBLE PRECISION THE        
!-               VARIABLES NEEDED FOR CALCULATIONS WHICH NEED MORE      
!-               THAN 32 BITS TO OBTAIN SUFFICIENT ACCURACY (I.E.       
!-               THE FREQUENCIES). STATEMENTS FLAGGED WITH ! HAVE       
!-            C! STARTING IN COLUMN 1. THESE TWO CHARACTERS SHOULD BE   
!-               CHANGED TO BLANKS FOR COMPUTERS HAVING SINGLE          
!-               PRECISION LESS THAN 10 SIGNIFICANT DIGITS.             
!-                                                                      
!-   #    IDENTIFIES STATEMENTS THAT MAY BE USEFUL FOR ACCELERATED      
!-               FILE DATA TRANSFER UNDER CDC AND OTHER OPERATING       
!-               SYSTEMS ALLOWING BUFFERED I/0.                         
!-                                                                      
!-   >    IDENTIFIES STATEMENTS THAT MAY BE USEFUL FOR CONVERSION,      
!-               TYPICALLY SYSTEM SPECIFIC CALLS (I.E. DATE, TIME,      
!-               CPU TIME, RANDOM NUMBER, ETC.).                        
!-                                                                      
!---------------------------------------------------------------------- 
!                                                                       
      USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
     &                     MXISOTPL,                                    &
     &                     MXMOL, MX_XS, MXTRAC, MXSPC, NMAXCO,         &
     &                     IPTS, IPTS2
!                                                                       
      IMPLICIT REAL*8           (V) 
!                                                                       
      character*8      XID,       HMOLID,      YID,HDATE,HTIME 
      real*8               SECANT,       XALTZ 
      real             solvar(2)
!                                                                       
      LOGICAL OP 
!%%%%%LINUX_PGI90 (-i8)%%%%%      integer*4 iostat                      
      CHARACTER CXID*80,CFORM*11,XID8*8,IDCNTL*6 
      CHARACTER*55 PTHT3M,PTHODI,PTHODTU,PTHODTD,CTAPE3 
      CHARACTER*11 PTHRDRU,PTHRDRD 
      CHARACTER*3 PTHDIR,AJID 
                                             ! change if PTHDIR//PTHRDRD
      CHARACTER*17 FULLPTH 
      CHARACTER*10 HFMODI,HFMODTU,HFMODTD,HFMRDR 
      CHARACTER*9 CT6FIL 
      CHARACTER*18 HNAMLBL,HNAMCNT,HNAMFFT,HNAMATM,HNAMLOW,HNAMNCG,     &
     &             HNAMOPR,HNAMPLT,HNAMPST,HNAMTST,HNAMUTL,HNAMXMR,     &
     &             hnmnlte                                              
      CHARACTER*18 HNAMSOL 
      CHARACTER*18 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,           &
     &             HVROPR,HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR,           &
     &             hvnlte                                               
      CHARACTER*18 HVRSOL 
!                                                                       
      CHARACTER*1 CONE,CTWO,CTHREE,CFOUR,CA,CB,CC,CDOL,CPRCNT,CBLNK 
      CHARACTER*1 CMRG(2),CXIDA(80) 
!                                                                       

!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &                MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,mx_xs=38,&
!     &                MXTRAC=22,MXSPC=5)                                
!                                                                       
!     ----------------------------------------------------------------  
                                                                        
!     unlabeled common                                                  
                                                                        
      COMMON COMSTR(250,9) 
      COMMON R1(3600),R2(900),R3(225) 
                                                                        
!     ----------------------------------------------------------------  
                                                                        
!     Parameter and common blocks for direct input of emissivity and    
!     reflectivity function values                                      
!                                                                       
!      PARAMETER (NMAXCO=4040) 
      COMMON /EMSFIN/ V1EMIS,V2EMIS,DVEMIS,NLIMEM,ZEMIS(NMAXCO) 
      COMMON /RFLTIN/ V1RFLT,V2RFLT,DVRFLT,NLIMRF,ZRFLT(NMAXCO) 
!     ----------------------------------------------------------------  
!                                                                       
!     -------------------------                                         
      CHARACTER*6  CMOL,CSPC 
      common /cmol_nam/ cmol(mxmol),cspc(mxspc) 
!                                                                       
!     -------------------------                                         
!     Common blocks for analytic derivative                             
! note: comments may not be consistent - if doing a search, use both    
!  "derivative" and "jacobian"                                          
!     -------------------------                                         
      COMMON /ADRPNM/ PTHT3M,PTHODI,PTHODTU,PTHODTD 
      COMMON /ADRPTH/ PTHDIR,PTHRDRU,PTHRDRD,AJID 
      COMMON /ADRFRM/ HFMODI,HFMODTU,HFMODTD,HFMRDR 
      COMMON /IADFLG/ NSPCRT,IMRGSAV 
      COMMON /ADRFIL/ KODFIL,kradtot,KTEMP,KFILAD,K_REFTRA,k_rddn_sfc 
                                                                        
      common /dlaydlev/ilevdx,imoldx,iup_dn,                            &
     &    dxdL(mxlay,0:mxmol),dxdU(mxlay,0:mxmol)                       
                                                                        
! note: from continuum module                                           
!          ipts  = same dimension as ABSRB                              
!          ipts2 = same dimension as C                                  
!      parameter (ipts=5050,ipts2=6000) 
      common /CDERIV/ icflg,idum,v1absc,v2absc,dvabsc,nptabsc,delT_pert,&
     &    dqh2oC(ipts),dTh2oC(ipts),dUh2o                               
                                                                        
!     -------------------------                                         
!                                                                       
      DIMENSION IDCNTL(16),IFSDID(17),IWD(2),IWD2(2),IWD3(2),IWD4(2) 
!                                                                       
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
     &                AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,   &
     &                DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,   &
     &                ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,  &
     &                EXTID(10)                                         
      CHARACTER*8  EXTID 
                                                                        
      COMMON /CVRLBL/ HNAMLBL,HVRLBL 
      COMMON /CVRCNT/ HNAMCNT,HVRCNT 
      COMMON /CVRFFT/ HNAMFFT,HVRFFT 
      COMMON /CVRATM/ HNAMATM,HVRATM 
      COMMON /CVRLOW/ HNAMLOW,HVRLOW 
      COMMON /CVRNCG/ HNAMNCG,HVRNCG 
      COMMON /CVROPR/ HNAMOPR,HVROPR 
      COMMON /CVRPST/ HNAMPST,HVRPST 
      COMMON /CVRPLT/ HNAMPLT,HVRPLT 
      COMMON /CVRTST/ HNAMTST,HVRTST 
      COMMON /CVRUTL/ HNAMUTL,HVRUTL 
      COMMON /CVRXMR/ HNAMXMR,HVRXMR 
      COMMON /CVNLTE/ HNMNLTE,HVNLTE 
      COMMON /CVRSOL/ HNAMSOL,HVRSOL 
                                                                        
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,pad_3,&
     &                angle_path,secant_diffuse,secant_path,diffuse_fac 
!                                                                       
      character*1 hmol_scal,h_xs_scal 
      character*1 surf_refl 
      character*3 pad_3 
!                                                                       
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &                NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,    &
     &                NLTEFL,LNFIL4,LNGTH4,IBRD                              
      COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),          &
     &                MSPNL1(MXLAY),MSLAY1,ISFILE,JSFILE,KSFILE,        &
     &                LSFILE,MSFILE,IEFILE,JEFILE,KEFILE                
      COMMON /LASIV/ VLAS,ILAS 
      COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,NOZERO,NP,H1F,H2F,       &
     &                ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH,       &
     &                XVBAR, HMINF,PHIF,IERRF,HSPACE                    
      COMMON /MSCONS/ AIRMAS(MXLAY),TGRND,SEMIS(3),HMINMS,HMAXMS,       &
     &                MSFLAG,                                           &
     &                MSWIT,IODFIL,MSTGLE                               
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN 
      COMMON /HDRF/ V1D,V2D,DVD,NLND,IWLD 
      COMMON /NGTH/ VD,SD,AD,EPD,MOLD,HWHD,TMPD,PSHD,FLGD,ILS2D 
      COMMON /HDRL/ V1LD,VL2D,NLD,NWDS,ILST3D 
      COMMON /RCNTRL/ ILNFLG 
      COMMON /FLFORM/ CFORM 
      COMMON /IODFLG/ DVOUT 
      COMMON /CNTSCL/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL 
      common /profil_scal/ nmol_scal,hmol_scal(64),xmol_scal(64),       &
     &                     n_xs_scal,h_xs_scal(64),x_xs_scal(64)        
      COMMON /PATH_ISOTPL/ ISOTPL,NISOTPL,                              &
     &                     ISOTPL_FLAG(MXMOL,MXISOTPL),                 &
     &                     ISOTPL_MAIN_FLAG(MXMOL),                     &
     &                     MOLNUM(MXMOL*MXISOTPL),                      &
     &                     ISOTPLNUM(MXMOL*MXISOTPL),                   &
     &                     WKI(MXMOL,MXISOTPL)             
!                                                                       
                                                                        
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
     &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
     &                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF  
                                                                        
      EQUIVALENCE (FSCDID( 1),IFSDID(1),IHIRAC), (FSCDID( 2),ILBLF4),   &
     &            (FSCDID( 3),IXSCNT),           (FSCDID( 4),IAERSL),   &
     &            (FSCDID( 5),IEMIT),            (FSCDID( 6),ISCAN),    &
     &            (FSCDID( 7),IPLOT),            (FSCDID( 8),IPATHL),   &
     &            (FSCDID( 9),JRAD),             (FSCDID(10),ITEST),    &
     &            (FSCDID(11),IMRG),             (FSCDID(12),SCNID),    &
     &            (FSCDID(13),HWHM),             (FSCDID(14),IDABS),    &
     &            (FSCDID(15),IATM),             (FSCDID(16),LAYR1),    &
     &            (FSCDID(17),NLAYFS),           (YI1,   IMULT),        &
     &            (YID(1),HDATE),                (YID(2),HTIME),        &
     &            (YID(8),LH2SAV),               (YID(9),LH1SAV),       &
     &            (YID(10),LTNSAV,dv_lbl)                               
                                                                        
!                  yid(3) through yid(7) cotains information from LBLLOW
                                                                        
                                                                        
      EQUIVALENCE (IWD(1),XID(1)) , (IWD2(1),V1D) , (IWD3(1),VD) ,      &
     &            (IWD4(1),V1LD)                                        
      EQUIVALENCE (CXID,CXIDA(1)) 
!                                                                       
      DATA IDCNTL / ' HIRAC',' LBLF4',' CNTNM',' AERSL',' EMISS',       &
     &              ' SCNFN',' FILTR','  PLOT','  TEST','  IATM',       &
     &              '  IMRG','  ILAS',' OPDEP',' XSECT','ISOTPL' ,      &
     &              '  IBRD'/
!                                                                       
      DATA CONE / '1'/,CTWO / '2'/,CTHREE / '3'/,CFOUR / '4'/,          &
     &     CA / 'A'/,CB / 'B'/,CC / 'C'/                                
      DATA CDOL / '$'/,CPRCNT / '%'/,CBLNK / ' '/,CXIDA / 80*' '/ 
      DATA XID8 / ' LBLRTM '/ 
!                                                                       
!     DATA CFORM / 'BUFFERED   '/                                       
!     DATA CFORM / 'UNFORMATTED'/                                       
!                                                                       
      DATA I_10/10/ 
!                                                                       
!     set the svn version number                                        
!                                                                       
      HVRLBL  = '$Revision$'
!                                                                       
!     Set ILNFLG to default (no line rejection files kept)              
!                                                                       
      ILNFLG = 0 
!                                                                       
!     Set name of output TAPE6, depending upon type calculation         
!                                                                       
      CT6FIL = 'TAPE6    ' 
!                                                                       
!     -------------------------                                         
!                                                                       
      KODFIL = 17 
      kradtot = 18 
      KFILAD = 19 
      KTEMP  = 88 
      K_REFTRA = 89 
      k_rddn_sfc = 90 
      CALL QNTIFY(PTHODI,HFMODI) 
      CALL QNTIFY(PTHODTU,HFMODTU) 
      CALL QNTIFY(PTHODTD,HFMODTD) 
      fullpth=pthdir//pthrdru//ajid 
      CALL QNTIFY(fullpth,HFMRDR) 
      CTAPE3 = PTHT3M 
                                                                        
      IMRGSAV = 0 
                                                                        
      DO I = 1,MXISOTPL
         DO M = 1,MXMOL
            ISOTPL_FLAG(M,I) = 0
         END DO
      END DO
!                                                                       
! analytic jacobians:                                                   
!   set flag for layer2level conversion                                 
!   this will be reset if it is appropriate to use layer2level conversio
      imoldx=-999 
                                                                        
!     open default files: read, print, punch                            
      IRD = 55 
      OPEN (IRD,FILE='TAPE5',STATUS='UNKNOWN') 
      IPR = 66 
      OPEN (IPR,FILE=CT6FIL,STATUS='UNKNOWN') 
      IPU = 7 
                                                                        
   10 WRITE (IPR,900) 
!                                                                       
!     -------------------------                                         
!                                                                       
!     OTHER FILE ASSIGNMENTS                                            
!                                                                       
      LNFIL4 = 9 
      OPEN (LNFIL4,FILE='TAPE9',STATUS='UNKNOWN',FORM=CFORM) 
      KFILE = 10 
      OPEN (KFILE,FILE='TAPE10',STATUS='UNKNOWN',FORM=CFORM) 
      LFILE = 11 
      OPEN (LFILE,FILE='TAPE11',STATUS='UNKNOWN',FORM=CFORM) 
      MFILE = 12 
      OPEN (MFILE,FILE='TAPE12',STATUS='UNKNOWN',FORM=CFORM) 
      IODFIL = 19 
      IEXFIL = 20 
      KKSTOR = KFILE 
      LLSTOR = LFILE 
      MMSTOR = MFILE 
!                                                                       
      IENDPL = 0 
      MSFLAG = 0 
      HMINMS = 0.0 
      DVSET = 0.0 
      HMAXMS = 15.0 
      MSWIT = 0 
      MSTGLE = 0 
      ONEPL = 1.001 
      ONEMI = 0.999 
      ARGMIN = 34. 
      EXPMIN = EXP(-ARGMIN) 
!                                                                       
      REWIND LFILE 
      REWIND MFILE 
      LSTWDF = -654321 
      NFHDRF = NWDL(IWD,LSTWDF) 
      IWLD = -654321 
      NPHDRF = NWDL(IWD2,IWLD) 
      ILS2D = -654321 
      NLNGTH = NWDL(IWD3,ILS2D) 
      ILST3D = -654321 
      NPHDRL = NWDL(IWD4,ILST3D) 
                                                                        
!                                                                       
      LOWFLG = 0 
      IREAD = 0 
      NOPR = 0 
!                                                                       
!     XID = 80 CHARACTERS OF USER IDENTIFICATION                        
!                                                                       
   20 READ (IRD,905,END=80) CXID 
      IF (CXIDA(1).EQ.CPRCNT) GO TO 90 
      IF (CXIDA(1).NE.CDOL) GO TO 20 
      CXIDA(1) = CBLNK 
      READ (CXID,910) (XID(I),I=1,10) 
      READ (XID8,910) XID(10) 
      CALL LBLDAT(HDATE) 
      CALL FTIME (HTIME) 
      WRITE (IPR,915) XID,HDATE,HTIME 
      DO 30 I = 1, 17 
         IFSDID(I) = -99 
   30 END DO 
      FSCDID(12) = -99. 
      FSCDID(13) = -99. 
      VLAS = -99. 
      IDABS = 0 
      TIME0 = -99. 
      CALL CPUTIM (TIME0) 
      WRITE (IPR,920) TIME0 
!                                                                       
      READ(IRD,925,END=80) IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,           &
     &                      ISCAN,IFILTR,IPLOT,ITEST,IATM,CMRG,ILAS,    &
     &                      IOD,IXSECT,IRAD,MPTS,NPTS,ISOTPL,IBRD            
!                                                                       
                                                                        
      ICNTNM_sav = ICNTNM 
                                                                        
!     Set continuum flags as needed                                     
                                                                        
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
!                                                                       
      IF (ICNTNM.EQ.0) THEN 
         XSELF = 0.0 
         XFRGN = 0.0 
         XCO2C = 0.0 
         XO3CN = 0.0 
         XO2CN = 0.0 
         XN2CN = 0.0 
         XRAYL = 0.0 
      ELSEIF (ICNTNM.EQ.1) THEN 
         XSELF = 1.0 
         XFRGN = 1.0 
         XCO2C = 1.0 
         XO3CN = 1.0 
         XO2CN = 1.0 
         XN2CN = 1.0 
         XRAYL = 1.0 
      ELSEIF (ICNTNM.EQ.2) THEN 
         XSELF = 0.0 
         XFRGN = 1.0 
         XCO2C = 1.0 
         XO3CN = 1.0 
         XO2CN = 1.0 
         XN2CN = 1.0 
         XRAYL = 1.0 
         ICNTNM = 1 
      ELSEIF (ICNTNM.EQ.3) THEN 
         XSELF = 1.0 
         XFRGN = 0.0 
         XCO2C = 1.0 
         XO3CN = 1.0 
         XO2CN = 1.0 
         XN2CN = 1.0 
         XRAYL = 1.0 
         ICNTNM = 1 
      ELSEIF (ICNTNM.EQ.4) THEN 
         XSELF = 0.0 
         XFRGN = 0.0 
         XCO2C = 1.0 
         XO3CN = 1.0 
         XO2CN = 1.0 
         XN2CN = 1.0 
         XRAYL = 1.0 
         ICNTNM = 1 
      ELSEIF (ICNTNM.EQ.5) THEN 
         XSELF = 1.0 
         XFRGN = 1.0 
         XCO2C = 1.0 
         XO3CN = 1.0 
         XO2CN = 1.0 
         XN2CN = 1.0 
         XRAYL = 0.0 
         ICNTNM = 1 
      ELSEIF (ICNTNM.EQ.6) THEN 
         READ(IRD,*) XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN, XRAYL 
         ICNTNM = 1 
      ENDIF 
                                                                        
      IXSCNT = IXSECT*10+ICNTNM 
!                                                                       
!     *********************** SOLAR RADIANCE ***********************    
!                                                                       
!     Check to be sure that no radiative transfer is being done other   
!     than solar and calculate attenuated solar radiation.              
!                                                                       
!                                                                       
      IF (IEMIT.EQ.2) THEN 
         IF (IHIRAC+ILBLF4+ICNTNM+IAERSL.GT.0) THEN 
            WRITE(IPR,*) 'No radiative transfer calculation allowed',   &
            ' when IEMIT=2 (solar radiance)'                            
            STOP 'ERROR: IHIRAC+ILBLF4+ICNTNM+IAERSL > 0' 
         ENDIF 
         INFLAG = 0 
         IOTFLG = 0 
         READ(IRD,1010) INFLAG,IOTFLG,JULDAT,ISOLVAR,SCON,SOLCYCFRAC,SOLVAR
         IF (INFLAG.EQ.1) THEN 
!            IFILE = KFILE                                              
            IFILE = MFILE 
         ELSE 
            IFILE = MFILE 
         ENDIF 
         INQUIRE (UNIT=13,OPENED=OP) 
         IF (OP) CLOSE(13) 
         OPEN(UNIT=13,FILE='TAPE13',FORM='UNFORMATTED') 
         CALL SOLINT(IFILE,13,NPTS,INFLAG,IOTFLG,JULDAT,ISOLVAR,   &
                     SCON,SOLCYCFRAC,SOLVAR) 
         REWIND 13 
         GOTO 60 
      ENDIF 
!                                                                       
!     ***************************************************************   
!                                                                       
!                                                                       
!    OPEN LINFIL DEPENDENT UPON IHIRAC AND ITEST                        
!                                                                       
!     Linefile name specified by CTAPE3                                 
!                                                                       
      LINFIL = 3 
      IF (IHIRAC.GT.0) THEN 
         IF (ITEST.EQ.1) THEN 
            OPEN (LINFIL,FILE=CTAPE3,STATUS='NEW',FORM=CFORM) 
         ELSE 
            OPEN (LINFIL,FILE=CTAPE3,STATUS='OLD',FORM=CFORM) 
         ENDIF 
      ENDIF 
!                                                                       
!    CHECK CMRG TO SEE IF QUANTITY IS SINGLE DIGIT, DOUBLE DIGIT        
!    OR CHARACTER                                                       
!                                                                       
      IF (CMRG(2).EQ.CA) THEN 
         IMRG = 12 
      ELSEIF (CMRG(2).EQ.CB) THEN 
         IMRG = 22 
      ELSEIF (CMRG(2).EQ.CC) THEN 
         IMRG = 32 
      ELSE 
         READ (CMRG(2),930) IMRG 
         IF (CMRG(1).EQ.CONE) IMRG = IMRG+10 
         IF (CMRG(1).EQ.CTWO) IMRG = IMRG+20 
         IF (CMRG(1).EQ.CTHREE) IMRG = IMRG+30 
         IF (CMRG(1).EQ.CFOUR) IMRG = IMRG+40 
      ENDIF 
      IF (IPLOT.GT.0) IENDPL = 1 
!                                                                       
!     JRAD= -1  NO RADIATION TERM IN ABSORPTION COEFFICIENTS            
!     JRAD=  0  RADIATION TERM PUT IN BY PANEL                          
!     JRAD=  1  RADIATION TERM INCLUDED IN LINE STRENGTHS               
!                                                                       
      JRAD = 1 
      IF (IRAD.NE.0) JRAD = -1 
      WRITE (IPR,935) (IDCNTL(I),I=1,16) 
      WRITE (IPR,940) IHIRAC,ILBLF4,ICNTNM_sav,IAERSL,IEMIT,ISCAN,      &
     &    IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,ISOTPL,IBRD           
!                                                                       
      IF (IHIRAC.EQ.4) THEN 
         IF (IEMIT.NE.1) THEN 
            WRITE (IPR,950) 
            STOP ' IEMIT NE 1 FOR NLTE ' 
         ENDIF 
      ENDIF 
!                                                                       
!     Set IMULT equal to IOD, the flag for optical depth DV             
!                                                                       
      IMULT = IOD 
!                                                                       
      IF (IAERSL.GE.1 .and. iaersl.ne.5) LOWFLG = 1 
!                                                                       
      IAFIL = 14 
!                                                                       
!     IEXFIL=20                                                         
!                                                                       
      IF (IAERSL.GE.1 .and. iaersl.ne.5) THEN 
         OPEN (IAFIL,FILE='TAPE14',STATUS='UNKNOWN',FORM=CFORM) 
         OPEN (IEXFIL,FILE='TAPE20',STATUS='UNKNOWN',FORM=CFORM) 
         REWIND IEXFIL 
         LOWFLG = 1 
      ENDIF 
                                                                        
      NFILE = 13 
      MMRG = MOD(IMRG,I_10) 
      IF (MMRG.GE.3) THEN 
         OPEN (NFILE,FILE='TAPE13',STATUS='UNKNOWN',FORM=CFORM) 
         REWIND NFILE 
      ENDIF 
                                                                        
      NLTEFL = 4 
      IF (IHIRAC.EQ.4) THEN 
         OPEN (NLTEFL,FILE='TAPE4',STATUS='OLD') 
      ENDIF 
!                                                                       
!     TAPE39  IS AFGL PLOT FILE                                         
!                                                                       
      IPLFL = 39 
!                                                                       
      IF (ITEST.EQ.1) CALL TESTMM (LINFIL) 
!                                                                       
!     IHIRAC = 1 CALL HIRAC1     VOIGT                                  
!     IHIRAC = 2 CALL HIRACL     LORENTZ    NOT IMPLEMENTED             
!     IHIRAC = 3 CALL HIRACD     DOPPLER    NOT IMPLEMENTED             
!     IHIRAC = 4 CALL HIRACQ     NLTE VOIGT                             
!     IHIRAC = 9 NO LINE BY LINE CALCULATIONS; MECHANICS PURSUED        
!                                                                       
!     IF IEMIT .EQ. 1  PROGRAM WILL COMPUTE EMISSION                    
!                                                                       
!                                                                       
!     PRINT LINE FILE HEADER                                            
!                                                                       
      IF (IHIRAC.EQ.1.OR.IHIRAC.EQ.4.OR.ILBLF4.GE.1) CALL PRLNHD 
!                                                                       
!     PRINT CONTINUUM INFORMATION                                       
!                                                                       
      IF (ICNTNM.NE.0) THEN 
         CALL PRCNTM 
         WRITE(IPR,1020) XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL 
      ENDIF 
                                                                        
!                                                                       
      NOPR = 0 
      IF (MPTS.LT.0) NOPR = 1 
      LAYR1 = 1 
      LAYRFX = LAYR1 
!                                                                       
      IF ((IHIRAC+IAERSL+IATM+IMRG).LE.0) GO TO 60 
      IF ((IHIRAC+IAERSL+IEMIT+IATM+ILAS).GT.0) THEN 
!                                                                       
         READ (IRD,970,END=80,err=82) V1,V2,SAMPLE,DVSET,ALFAL0,AVMASS, &
         DPTMIN, DPTFAC,ILNFLG,DVOUT, nmol_scal,n_xs_scal               
         IF (ALFAL0.LE.0.) ALFAL0 = 0.04  
         IF (SAMPLE.LE.0) SAMPLE = 4
         IF (AVMASS.LE.0.) AVMASS = 36.0                                                                     
!_______________________________________________________________________
!     read information to scale the entire profile for indicated species
!                                                                       
         if (nmol_scal .gt. 0 ) then 
         if (nmol_scal .gt. mxmol) stop ' nmol_scal .gt. mxmol ' 
         read (ird,972) (hmol_scal(m),m=1,nmol_scal) 
         read (ird,973) (xmol_scal(m),m=1,nmol_scal) 
         endif 
!                                                                       
         if (n_xs_scal .gt. 0 ) then 
         if (n_xs_scal .gt. mx_xs) stop ' nmol_scal .gt. mx_xs ' 
         read (ird,972) (h_xs_scal(m),m=1,n_xs_scal) 
         read (ird,973) (x_xs_scal(m),m=1,n_xs_scal) 
         endif 
!_______________________________________________________________________
!                                                                       
!     OPEN LINE REJECTION FILES IF ILNFLG IS ONE OR TWO                 
!                                                                       
         IF (ILNFLG.EQ.1) THEN 
            OPEN(15,FILE='REJ1',STATUS='NEW',FORM='UNFORMATTED',        &
            IOSTAT=iostat)                                              
            if (IOSTAT.gt.0) stop 'REJ1 file is already created' 
            OPEN(16,FILE='REJ4',STATUS='NEW',FORM='UNFORMATTED') 
         ENDIF 
         IF (ILNFLG.EQ.2) THEN 
            OPEN(15,FILE='REJ1',STATUS='OLD',FORM='UNFORMATTED',        &
            IOSTAT=iostat)                                              
            if (IOSTAT.gt.0) stop 'REJ1 file does not exist' 
            OPEN(16,FILE='REJ4',STATUS='OLD',FORM='UNFORMATTED',        &
            IOSTAT=iostat)                                              
            if (IOSTAT.gt.0) stop 'REJ4 file does not exist' 
         ENDIF 
!                                                                       
!     IF DPTMIN < 0. SET TO DEFAULT (.0002)                             
!     IF DPTFAC < 0. SET TO DEFAULT (.001)                              
!                                                                       
         IF (DPTMIN.LT.0.) DPTMIN = .0002 
         IF (DPTFAC.LT.0.) DPTFAC = .001 
         IF (V2.LE.V1.AND.ILAS.EQ.0) ILAS = 1 
      ENDIF 
                                                                        
      IF (ILAS.GT.0) THEN 
         V2 = V1 
         VLAS = V1 
      ENDIF 
                                                                        
      TBOUND = 0. 
      TMPBND = 0. 
      EMISIV = 0. 
      BNDEMI(1) = 0. 
      BNDEMI(2) = 0. 
      BNDEMI(3) = 0. 
      BNDRFL(1) = 0. 
      BNDRFL(2) = 0. 
      BNDRFL(3) = 0. 
      IBPROP = 0 
                                                                        
!     Default surface type set to specular for all calculations,        
!     unless changed by following read.                                 
      surf_refl = 's' 
                                                                        
      IF (IEMIT.GT.0) THEN 
         READ (IRD,971,END=80) TMPBND,(BNDEMI(IBND),IBND=1,3), (BNDRFL( &
         IBND),IBND=1,3), surf_refl                                     
!                                                                       
         BNDTST = ABS(BNDRFL(1))+ABS(BNDRFL(2))+ABS(BNDRFL(3)) 
         IF (BNDTST.NE.0.) IBPROP = 1 
!                                                                       
!     **************************************************************    
!        If BNDEMI(1) < 0, read in coefficients from file 'EMISSIVITY'  
!        If BNDEMI(1) > 0, check to see if emissivity is reasonable     
!                                                                       
!        UNIT ICOEF used for input files                                
!                                                                       
         ICOEF = 13 
!                                                                       
         IF (BNDEMI(1).LT.0) THEN 
            OPEN (UNIT=ICOEF,FILE='EMISSIVITY', STATUS='OLD',IOSTAT=    &
            iostat)                                                     
                                                                        
            if ( iostat .gt. 0) &
                stop "FILE 'EMISSIVITY' FOR PATH BOUNDARY NOT FOUND"
!                                                                       
            CALL READEM(ICOEF) 
            CLOSE (ICOEF) 
         ELSE 
            XVMID = (V1+V2)/2. 
            EMITST = BNDEMI(1)+BNDEMI(2)*XVMID+BNDEMI(3)*XVMID*XVMID 
            IF (EMITST.LT.0..OR.EMITST.GT.1.) THEN 
               WRITE (IPR,975) XVMID,EMITST 
               STOP 'BNDEMI' 
            ENDIF 
         ENDIF 
!                                                                       
!        If BNDRFL(1) < 0, read in coefficients from file 'REFLECTIVITY'
!        If BNDRFL(1) > 0, check to see if reflectivity is reasonable   
!                                                                       
         IF (BNDRFL(1).LT.0) THEN 
            OPEN (UNIT=ICOEF,FILE='REFLECTIVITY', STATUS='OLD',IOSTAT=  &
            iostat)                                                     
                                                                        
            if ( iostat .gt. 0) &
                 stop "FILE 'REFLECTIVITY' FOR PATH BOUNDARY NOT FOUND"                                         
!                                                                       
            CALL READRF(ICOEF) 
            CLOSE (ICOEF) 
         ELSE 
            REFTST = BNDRFL(1)+BNDRFL(2)*XVMID+BNDRFL(3)*XVMID*XVMID 
            IF (REFTST.LT.0..OR.REFTST.GT.1.) THEN 
               WRITE (IPR,980) XVMID,REFTST 
               STOP 'BNDRFL' 
            ENDIF 
                                                                        
         ENDIF 
                                                                        
!     **************************************************************    
!                                                                       
!     TBOUND IS THE BOUNDARY TEMPERATURE. TBOUND=0. FOR NO BOUNDARY     
!     EMISIV IS THE BOUNDARY EMISSIVITY                                 
!     SET DEFAULT FOR EMISIV and Surface Reflection                     
!                                                                       
         EMITST = ABS(BNDEMI(1))+ABS(BNDEMI(2))+ABS(BNDEMI(3)) 
                                                                        
         EMISIV = BNDEMI(1) 
         TBOUND = TMPBND 
!                                                                       
!        set default surface reflection settings                        
         IF (surf_refl .eq. ' ') surf_refl = 's' 
!                                                                       
         WRITE (IPR,985) V1,V2,TBOUND,(BNDEMI(IBND),IBND=1,3), (BNDRFL( &
         IBND),IBND=1,3), surf_refl                                     
!                                                                       
      ENDIF 
!                                                                       
      ILASRD = 0 
   40 CONTINUE 
!                                                                       
      IF (ILASRD.GT.0) THEN 
         READ (IRD,990) VLAS 
         IF (VLAS.LT.0.) GO TO 70 
         V1 = VLAS 
         V2 = VLAS 
      ENDIF 
      IF (ILAS.EQ.2) ILASRD = 1 
!                                                                       
!     check that dvset and dvout are set properly                       
!                                                                       
      IF (IOD.EQ.0) THEN 
         IF (DVOUT.ne.0.) STOP 'DVOUT MUST BE ZERO FOR IOD=0' 
      ENDIF 
      IF (IOD.EQ.1) THEN 
         IF (DVSET.NE.0) STOP 'DVSET MUST BE ZERO FOR IOD=1' 
         IF (DVOUT.EQ.0.) STOP 'DVOUT MUST BE NONZERO FOR IOD=1' 
      ENDIF 
      IF (IOD.EQ.2) THEN 
         IF (DVSET.NE.0) STOP 'DVSET MUST BE ZERO FOR IOD=2' 
         IF (DVOUT.ne.0.) STOP 'DVOUT MUST BE ZERO FOR IOD=2' 
      ENDIF 
      IF (IOD.EQ.3) THEN 
         IF (DVOUT.ne.0.) STOP 'DVOUT MUST BE ZERO FOR IOD=3 ' 
      ENDIF 
      IF (IOD.EQ.4) THEN 
         IF (DVOUT.EQ.0.) STOP 'DVOUT MUST BE NONZERO FOR IOD=4' 
      ENDIF 
!                                                                       
!     -------------------------------------------------------------     
!                              CALL TREE                                
!     XLAYER                   ---------                                
!          \_ OPPATH                                                    
!                  \_ LBLATM, PATH, & LOWTRAN                           
!          \_ OPDPTH                                                    
!                  \_ CONTNM, LINF4, & HIRAC1                           
!          \_ SCNMRG                                                    
!                  \_ RDSCAN, SHRKSC, CNVRCT,                           
!                     CONVSC, PNLRCT,  & PANLSC                         
!     -------------------------------------------------------------     
!                                                                       
      IF (IHIRAC+IATM+IMRG.GT.0)                                        &
     &    CALL XLAYER (MPTS,NPTS,LFILE,MFILE,NFILE)                     
!                                                                       
      IF ((IAERSL.EQ.1.OR.IAERSL.EQ.7) .AND. (IEMIT.EQ.0)               &
     &                        .and. (ihirac+imrg.gt.0)) THEN            
         REWIND MFILE 
         REWIND IAFIL 
         REWIND IEXFIL 
!                                                                       
         CALL ADARSL (NPTS,IEXFIL,MFILE,IAFIL,IEMIT) 
!                                                                       
      ENDIF 
!                                                                       
      IF (IMRG.EQ.1) MFILE = KFILE 
      IF ((MMRG.GE.3).AND.(MMRG.NE.9)) MFILE = NFILE 
      IF (ILAS.GT.0) THEN 
         IATM = 9 
         JAERSL = 0 
         IF (IAERSL.GE.1 .and. iaersl.ne.5 .AND. IEMIT.EQ.1) JAERSL = 1 
!                                                                       
         CALL LASER (VLAS,MFILE,JAERSL) 
!                                                                       
         IF (IAERSL.GE.1 .and. iaersl.ne.5 .AND.IEMIT.EQ.0) CALL LASER (&
         VLAS,IAFIL,1)                                                  
!                                                                       
         IF (ILAS.GE.2) GO TO 40 
      ENDIF 
   60 CONTINUE 
!                                                                       
      REWIND KFILE 
      REWIND LFILE 
      REWIND MFILE 
!                                                                       
      LFILE = LLSTOR 
      IF (ISCAN.EQ.1) CALL SCANFN (MFILE,LFILE) 
!                                                                       
      IF (ISCAN.EQ.2) CALL INTRPL (MFILE,LFILE) 
!                                                                       
      IF (ISCAN.EQ.3) CALL FFTSCN (MFILE,LFILE) 
!                                                                       
      IF (IFILTR.EQ.1) CALL FLTRFN (MFILE) 
!                                                                       
      IF (IPLOT.NE.0) CALL PLTLBL (IENDPL) 
!                                                                       
   70 CONTINUE 
      CALL CPUTIM (TIME1) 
      TIME = TIME1-TIME0 
      WRITE (IPR,995) TIME1,TIME 
      KFILE = KKSTOR 
      LFILE = LLSTOR 
      MFILE = MMSTOR 
                                                                        
      close(kfile) 
      close(lfile) 
      close(mfile) 
      close(nfile) 
                                                                        
      GO TO 10 
                                                                        
   80 CONTINUE 
      IF (IENDPL.EQ.1) CALL ENDPLT 
      STOP ' LBLRTM EXIT; EOF ON TAPE 5 ' 
                                                                        
   82 continue 
      stop ' Error Reading Record with V1,V2, ... ' 
                                                                        
!                                                                       
   90 CONTINUE 
      WRITE (IPR,915) XID,HDATE,HTIME 
      WRITE(IPR,1000) HNAMLBL,HVRLBL,HNAMCNT,HVRCNT,                    &
     &                HNAMFFT,HVRFFT,HNAMATM,HVRATM,                    &
     &                HNAMLOW,HVRLOW,HNAMNCG,HVRNCG,                    &
     &                HNAMOPR,HVROPR,HNAMPST,HVRPST,                    &
     &                HNAMPLT,HVRPLT,HNAMTST,HVRTST,                    &
     &                HNAMXMR,HVRXMR,HNAMUTL,HVRUTL,                    &
     &                HNAMSOL,HVRSOL,hnmnlte,hvnlte                     
      IF (IENDPL.EQ.1) CALL ENDPLT 
      STOP ' LBLRTM EXIT ' 
!                                                                       
  900 FORMAT ('1') 
  905 FORMAT (A80) 
  910 FORMAT (10A8) 
  915 FORMAT ('0',10A8,2X,2(1X,A8,1X)) 
  920 FORMAT ('0  TIME ENTERING LBLRTM  ',F15.4) 
  925 FORMAT (10(4X,I1),3X,2A1,3(4X,I1),I1,I4,1X,I4,4X,I1,4x,I1,4x,I1) 
!                                                                       
  930 FORMAT (I1) 
  935 FORMAT (16(A6,3X)) 
  940 FORMAT (1X,I4,15I9) 
!                                                                       
  950 FORMAT ('0 IEMIT=0 IS NOT IMPLEMENTED FOR NLTE ',/,               &
     &        '  CHANGE IEMIT TO 1 OR IHIRAC TO 1 ')                    
  970 FORMAT (8E10.3,4X,I1,5x,e10.3,3X,i2,3x,i2) 
  971 FORMAT (7E10.3,4X,A1) 
  972 FORMAT (64a1) 
  973 FORMAT (8e15.7) 
  975 FORMAT ('0 FOR VNU = ',F10.3,' THE EMISSIVITY = ',E10.3,          &
     &        ' AND IS NOT BOUNDED BY (0.,1.) ')                        
  980 FORMAT ('0 FOR VNU = ',F10.3,' THE REFLECTIVITY = ',E10.3,        &
     &        ' AND IS NOT BOUNDED BY (0.,1.) ')                        
  985 FORMAT (5(/),'0*********** BOUNDARY PROPERTIES ***********',/,    &
     &        '0 V1(CM-1) = ',F12.4,/,'0 V2(CM-1) = ',F12.4,/,          &
     &        '0 TBOUND   = ',F12.4,5X,'BOUNDARY EMISSIVITY   = ',      &
     &        3(1PE11.3),/,'0',29X,'BOUNDARY REFLECTIVITY = ',          &
     &        3(1PE11.3),/,'0',29X,' SURFACE REFLECTIVITY = ', A1)      
  990 FORMAT (F20.8) 
  995 FORMAT ('0 TIME  LEAVING LBLRTM ',F15.4,' TOTAL',F15.4) 
 1000 FORMAT ('0 Modules and versions used in this calculation:',/,/,   &
     &         7(5X,a18,2X,A18,10X, a18,2X,A18,/))                      
 1010 FORMAT (2I5,2X,I3,I5,f10.4,5x,f10.4,5x,2f10.5) 
 1020 FORMAT (/,'  The continuum scale factors are as follows: ',       &
     &         /,5x,'H2O Self:    ',f10.3,                              &
     &         /,5x,'H2O Foreign: ',f10.3,                              &
     &         /,5x,'CO2:         ',f10.3,                              &
     &         /,5x,'O3:          ',f10.3,                              &
     &         /,5x,'O2:          ',f10.3,                              &
     &         /,5x,'N2:          ',f10.3,                              &
     &         /,5x,'Rayleigh:    ',f10.3,/)                            
                                                                        
!                                                                       
end program LBLRTM
!********************************************************************** 
!                                                                       
      BLOCK DATA 
      USE lblparams, ONLY: MXFSC, MXLAY, MXMOL, MXSPC, IPTS, IPTS2
      IMPLICIT REAL*8           (V) 
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3) 
      COMMON /FLFORM/ CFORM 
      COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),          &
     &                MSPNL1(MXLAY),MSLAY1,ISFILE,JSFILE,KSFILE,        &
     &                LSFILE,MSFILE,IEFILE,JEFILE,KEFILE                
!                                                                       
      COMMON /CVRLBL/ HNAMLBL,HVRLBL 
      COMMON /CVRCNT/ HNAMCNT,HVRCNT 
      COMMON /CVRFFT/ HNAMFFT,HVRFFT 
      COMMON /CVRATM/ HNAMATM,HVRATM 
      COMMON /CVRLOW/ HNAMLOW,HVRLOW 
      COMMON /CVRNCG/ HNAMNCG,HVRNCG 
      COMMON /CVROPR/ HNAMOPR,HVROPR 
      COMMON /CVRPST/ HNAMPST,HVRPST 
      COMMON /CVRPLT/ HNAMPLT,HVRPLT 
      COMMON /CVRTST/ HNAMTST,HVRTST 
      COMMON /CVRUTL/ HNAMUTL,HVRUTL 
      COMMON /CVRXMR/ HNAMXMR,HVRXMR 
      COMMON /CVNLTE/ HNMNLTE,HVNLTE 
      COMMON /CVRSOL/ HNAMSOL,HVRSOL 
                                                                        
      COMMON /ADRPNM/ PTHT3M,PTHODI,PTHODTU,PTHODTD 
      COMMON /ADRPTH/ PTHDIR,PTHRDRU,PTHRDRD,AJID 
      COMMON /CNTSCL/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL 
                                                                        
!     -------------------------                                         
!      PARAMETER (MXMOL=39,MXSPC=5) 
      common /cmol_nam/ cmol(mxmol),cspc(mxspc) 
      CHARACTER*6  CMOL,CSPC 
!     -------------------------                                         
!                                                                       
      CHARACTER CFORM*11 
                                                                        
!                                                                       
      CHARACTER*18 HNAMLBL,HNAMCNT,HNAMFFT,HNAMATM,HNAMLOW,HNAMNCG,     &
     &             HNAMOPR,HNAMPLT,HNAMPST,HNAMTST,HNAMUTL,HNAMXMR,     &
     &             hnmnlte                                              
      CHARACTER*18 HNAMSOL 
                                                                        
      CHARACTER*18 HVRLBL,HVRCNT,HVRFFT,HVRATM,HVRLOW,HVRNCG,           &
     &             HVROPR,HVRPLT,HVRPST,HVRTST,HVRUTL,HVRXMR,           &
     &             hvnlte                                               
      CHARACTER*18 HVRSOL 
!                                                                       
      CHARACTER*55 PTHT3M,PTHODI,PTHODTU,PTHODTD 
      CHARACTER*11 PTHRDRU,PTHRDRD 
      CHARACTER*3  PTHDIR,AJID 
!                                                                       
      DATA PTHT3M /'TAPE3'/ 
                                                                        
! *note that PTHRDRU and PTHRDRD must have the same length string       
! *also note that you must change the size declaration for PTHDIR if    
! you change the name of the path (in all routines that carry the common
      DATA PTHODI/'ODint_'/,                                            &
     &     PTHRDRD/'RDderivDNW_'/,PTHRDRU/'RDderivUPW_'/,               &
     &     PTHDIR/'AJ/'/,AJID/'xx_'/                                    
                                                                        
                                                                        
      COMMON /EMDTSV/ BBDSAV(2410),FSAV(2410),BBDLSAV(2410) 
      data bbdsav,fsav,bbdlsav/2410*0.0,2410*0.0,2410*0.0/ 
                                                                        
      DATA CFORM / 'UNFORMATTED'/ 
!                                                                       
      DATA IOD / 0 /,IDIR / 0 /,ITOP / 0 /,ISURF / 0 /,MSPTS / 0 /,     &
     &     MSPANL /MXLAY*0/,MSPNL1 /MXLAY*0/,ISFILE / 0 /,JSFILE / 0 /, &
     &     KSFILE / 0 /,LSFILE / 0 /,MSFILE / 0 /,IEFILE / 0 /,         &
     &     JEFILE / 0 /,KEFILE / 0 /,MSLAY1 / 0 /                       
!                                                                       
      DATA XSELF / 1 /,XFRGN / 1 /, XCO2C / 1 /, XO3CN / 1 /,           &
     &     XO2CN / 1 /,XN2CN / 1 /, XRAYL / 1 /                         
!                                                                       
!     ASSIGN DEFAULT MODULE NAMES                                       
!                                                                       
      DATA HNAMLBL / '         lblrtm.f:' /,                            &
     &     HNAMCNT / '         contnm.f:' /,                            &
     &     HNAMFFT / '         fftscn.f:' /,                            &
     &     HNAMATM / '         lblatm.f:' /,                            &
     &     HNAMLOW / '         lbllow.f:' /,                            &
     &     HNAMNCG / '        ncargks.f:' /,                            &
     &     HNAMOPR / '          oprop.f:' /,                            &
     &     HNAMPST / '        postsub.f:' /,                            &
     &     HNAMPLT / '         pltlbl.f:' /,                            &
     &     HNAMTST / '         testmm.f:' /,                            &
     &     HNAMUTL / '       util_xxx.f:' /,                            &
     &     HNAMXMR / '         xmerge.f:' /,                            &
     &     hnmnlte / '         nonlte.f:' /                             
      DATA HNAMSOL / '          solar.f:' / 
!                                                                       
!     ASSIGN CVS VERSION NUMBER TO MODULES                              
!                                                                       
      DATA HVRLBL / '   NOT USED       ' /,                             &
     &     HVRCNT / '   NOT USED       ' /,                             &
     &     HVRFFT / '   NOT USED       ' /,                             &
     &     HVRATM / '   NOT USED       ' /,                             &
     &     HVRLOW / '   NOT USED       ' /,                             &
     &     HVRNCG / '   NOT USED       ' /,                             &
     &     HVROPR / '   NOT USED       ' /,                             &
     &     HVRPST / '   NOT USED       ' /,                             &
     &     HVRPLT / '   NOT USED       ' /,                             &
     &     HVRTST / '   NOT USED       ' /,                             &
     &     HVRUTL / '   NOT USED       ' /,                             &
     &     HVRXMR / '   NOT USED       ' /,                             &
     &     hvnlte / '   NOT USED       ' /                              
      DATA HVRSOL / '   NOT USED       ' / 
!                                                                       
!     -------------------------                                         
!     Variables for analytic derivative                                 
!                                                                       
      DATA CMOL   /                                                     &
     &     '  H2O ','  CO2 ','   O3 ','  N2O ','   CO ','  CH4 ',       &
     &     '   O2 ','   NO ','  SO2 ','  NO2 ','  NH3 ',' HNO3 ',       &
     &     '   OH ','   HF ','  HCL ','  HBR ','   HI ','  CLO ',       &
     &     '  OCS ',' H2CO ',' HOCL ','   N2 ','  HCN ','CH3CL ',       &
     &     ' H2O2 ',' C2H2 ',' C2H6 ','  PH3 ',' COF2 ','  SF6 ',       &
     &     '  H2S ','HCOOH ','  HO2 ','    O ','ClONO2','  NO+ ',       &
     &     ' HOBr ',' C2H4 ','CH3OH ',' CH3Br',' CH3CN','  CF4 ',       &
     &     ' C4H2 ',' HC3N ','   H2 ','   CS ','  SO3 '/               
                                                                        
      DATA CSPC   / 'T LAYR','T SURF','SFC EM','SFC RF','LOW PR' / 
!                                                                       
! --- for analytic derivative ---                                       
! note: from continuum module                                           
!          ipts  = same dimension as ABSRB                              
!          ipts2 = same dimension as C                                  
!      parameter (ipts=5050,ipts2=6000) 
      common /CDERIV/ icflg,idum,v1absc,v2absc,dvabsc,nptabsc,delT_pert,&
     &    dqh2oC(ipts),dTh2oC(ipts),dUh2o                               
      data icflg /-999/ 
! --- for analytic derivative ---                                       
                                                                        
end block data
      FUNCTION NWDL (IWD,ILAST) 
!                                                                       
      DIMENSION IWD(*) 
!                                                                       
      DO 10 I = 1, 9000 
         IF (IWD(I).EQ.ILAST) THEN 
            NWDL = I-1 
            RETURN 
         ENDIF 
   10 END DO 
!                                                                       
      STOP ' NWDL - IWD,ILAST ' 
!                                                                       
end function NWDL
!                                                                       
!     -------------------------------------------------------------     
!                                                                       
      SUBROUTINE ENDFIL (IFILE) 
!                                                                       
      DIMENSION IDUM(6) 
      DATA IDUM / 6*-99 / 
!                                                                       
      CALL BUFOUT (IFILE,IDUM(1),6) 
!                                                                       
      RETURN 
!                                                                       
end subroutine ENDFIL
!                                                                       
!     -------------------------------------------------------------     
!                                                                       
      subroutine endfil_4 (ifile) 
!                                                                       
      integer*4 ifile 
      integer*4 idum(6) 
      data idum / 6*-99 / 
                                                                        
!                                                                       
      write(ifile) (idum(i),i=1,6) 
!                                                                       
      return 
!                                                                       
end subroutine endfil_4
!     -------------------------------------------------------------     
                                                                        
      SUBROUTINE SKIPFL (NUMFL,IFILE,IEOF) 
!                                                                       
      DIMENSION DUM(1) 
!                                                                       
      IF (NUMFL.LE.0) RETURN 
      ISKIP = 0 
   10 CALL BUFIN (IFILE,IEOF,DUM(1),1) 
      IF (IEOF.EQ.1) GO TO 10 
      ISKIP = ISKIP+1 
      IF (ISKIP.LT.NUMFL) GO TO 10 
!                                                                       
      RETURN 
!                                                                       
end subroutine SKIPFL
!                                                                       
!     -------------------------------------------------------------     
!                                                                       
      SUBROUTINE COPYFL (NPTS,KFILE,MFILE) 
!                                                                       
      IMPLICIT REAL*8           (V) 
!                                                                       
      COMMON TR(2410) 
!                                                                       
      character*8      XID,       HMOLID,      YID 
      real*8               SECANT,       XALTZ 
!                                                                       
      COMMON /EMIHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
     &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
     &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF  
      COMMON /PANL/ V1P,V2P,DVP,NLIMO 
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
     &              NLTEFL,LNFIL4,LNGTH4                                
!                                                                       
      DIMENSION XFILHD(2),PNLHD(2) 
      EQUIVALENCE (FSCDID(5),IEMIT) , (FSCDID(6),ISCAN) 
      equivalence (fscdid(12), xscid) 
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHD(1),V1P) 
!                                                                       
      CALL CPUTIM (TIME) 
      IF (NOPR.EQ.0) WRITE (IPR,900) TIME 
      CALL BUFIN (KFILE,KEOF,XFILHD(1),NFHDRF) 
      CALL BUFOUT (MFILE,XFILHD(1),NFHDRF) 
!                                                                       
   10 CALL BUFIN (KFILE,KEOF,PNLHD(1),NPHDRF) 
      IF (KEOF.LE.0) GO TO 20 
      CALL BUFOUT (MFILE,PNLHD(1),NPHDRF) 
      CALL BUFIN (KFILE,KEOF,TR(1),NLIMO) 
      CALL BUFOUT (MFILE,TR(1),NLIMO) 
!     Testing to see if there is just 1 data blocks (i.e.trans) or      
!     if there are 2 (i.e. trans. and rad) to read in                   
!     IF (IEMIT.EQ.0.OR.ISCAN.GT.0) GO TO 10                            
      if (iemit.eq.0.or.xscid.gt.0) go to 10 
      CALL BUFIN (KFILE,KEOF,TR(1),NLIMO) 
      CALL BUFOUT (MFILE,TR(1),NLIMO) 
      GO TO 10 
!                                                                       
   20 CALL CPUTIM (TIME1) 
      TIME = TIME1-TIME 
      IF (NOPR.EQ.0) WRITE (IPR,905) TIME 
!                                                                       
      RETURN 
!                                                                       
  900 FORMAT (' TIME AT THE START OF --COPYFL-- ',F10.3) 
  905 FORMAT (' TIME REQUIRED FOR --COPYFL -- ',F10.3) 
!                                                                       
end subroutine COPYFL
!                                                                       
!     -------------------------------------------------------------     
!                                                                       
      SUBROUTINE QNTIFY (CNAME,CFORM) 
!                                                                       
!     This subroutine counts the number of nonblank characters in       
!     a string, and creates a format to add two digits to the end       
!     of the nonblank string.                                           
!                                                                       
      CHARACTER*(*) CNAME 
      CHARACTER*55 CTEMP 
      CHARACTER*10 CFORM 
      CHARACTER*1  CTEMP1(55),BLANK 
!                                                                       
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &                NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,    &
     &                NLTEFL,LNFIL4,LNGTH4                              
!                                                                       
      EQUIVALENCE (CTEMP,CTEMP1(1)) 
!                                                                       
      DATA BLANK / ' ' / 
!                                                                       
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
   10 END DO 
!                                                                       
   20 CONTINUE 
      WRITE(CFORM,910) NCHAR 
!                                                                       
      RETURN 
!                                                                       
  900 FORMAT (' QNTIFY: PATH NOT LEFT JUSTIFIED: ',A55) 
  910 FORMAT ('(A',I2.2,',I3.3)') 
!                                                                       
end subroutine QNTIFY
!                                                                       
!     ---------------------------------------------------------------   
!                                                                       
      SUBROUTINE OPNODF(NLAYER,LAYER,PTHODL,HFMODL,kflg) 
!                                                                       
!     This subroutine opens file for calculating the radiance using     
!     precalculated optical depths                                      
!     (IEMIT = 1,IMRG=A/12,B/22,C/32,40,41)                             
!                                                                       
      LOGICAL OP 
      CHARACTER*57 FILE1 
      CHARACTER*55 PTHODL 
      CHARACTER*11 CFORM 
      CHARACTER*10 HFMODL 
!                                                                       
!     -------------------------                                         
!     Common block for analytic derivative                              
!     -------------------------                                         
      COMMON /ADRFIL/ KODFIL,kradtot,KTEMP,KFILAD,K_REFTRA,k_rddn_sfc 
!     -------------------------                                         
!                                                                       
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
     &              NLTEFL,LNFIL4,LNGTH4                                
!                                                                       
!           123456789-123456789-123456789-123456789-123456789-1234567   
      DATA FILE1 /                                                      &
     &     '                                                         '/ 
      DATA CFORM / 'UNFORMATTED' / 
!                                                                       
      if (kflg.lt.0) then 
         kopen = kodfil 
      else 
         kopen = kfile 
      endif 
                                                                        
      WRITE(IPR,910) LAYER,NLAYER 
      INQUIRE (UNIT=KOPEN,OPENED=OP) 
      IF (OP) CLOSE (KOPEN) 
                                                                        
      WRITE(FILE1,HFMODL) PTHODL,LAYER 
      OPEN(UNIT=KOPEN,FILE=FILE1,FORM=CFORM,STATUS='OLD') 
!                                                                       
!     Write procedure                                                   
!                                                                       
      WRITE(IPR,900) kopen, FILE1 
!                                                                       
      RETURN 
!                                                                       
  900 FORMAT ('       Opened layer optical depth file:  ',i5,' = ',A57) 
  910 FORMAT ('LAYER ',I5,' OF ',I5,':') 
!                                                                       
end subroutine OPNODF
!                                                                       
!     ---------------------------------------------------------------   
!                                                                       
      SUBROUTINE OPNRAD(nfile,NLAYER,LAYER) 
!                                                                       
!     This subroutine opens file for calculating the layer radiances    
!     (IEMIT = 1; IMRG=45,46)                                           
!                                                                       
      LOGICAL OP 
      CHARACTER*57 FILE1 
      CHARACTER*11 CFORM 
      CHARACTER*55 PTHRAD 
      CHARACTER*10 HFMRAD 
!                                                                       
!     Common block for layer radiances                                  
!     -------------------------                                         
      COMMON /RADLAY/ PTHRAD,HFMRAD 
!     -------------------------                                         
!                                                                       
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &              NLNGTH,KFILE,KPANEL,LINFIL,NFIL ,IAFIL,IEXFIL,      &
     &              NLTEFL,LNFIL4,LNGTH4                                
!                                                                       
!           123456789-123456789-123456789-123456789-123456789-1234567   
      DATA FILE1 /                                                      &
     &     '                                                         '/ 
      DATA CFORM / 'UNFORMATTED' / 
!                                                                       
      INQUIRE (UNIT=NFILE,OPENED=OP) 
      IF (OP) CLOSE (NFILE) 
      WRITE(FILE1,HFMRAD) PTHRAD,LAYER 
      OPEN(UNIT=NFILE,FILE=FILE1,FORM=CFORM,STATUS='UNKNOWN') 
!                                                                       
!     Write procedure                                                   
!                                                                       
      WRITE(IPR,900) FILE1 
!                                                                       
      RETURN 
!                                                                       
  900 FORMAT ('          Opened layer radiance file:  ',A57,/) 
  910 FORMAT ('LAYER ',I5,' OF ',I5,':') 
!                                                                       
end subroutine OPNRAD
!                                                                       
!     ---------------------------------------------------------------   
!                                                                       
      SUBROUTINE OPNDRV(NLAYER,LAYER,LAYTOT,ipathl) 
!                                                                       
!     This subroutine opens file for calculating the layer derivatives  
!     (IEMIT = 3)                                                       
!                                                                       
      IMPLICIT REAL*8 (V) 
                                                                        
      LOGICAL OP 
      CHARACTER*57 FILE1,FILE2,FILE3,FILE4 
      CHARACTER*11 CFORM 
      CHARACTER*55 CDUM1,PTHODI,PTHODTU,PTHODTD,pthrad 
      CHARACTER*11 PTHRDRU,PTHRDRD 
      CHARACTER*3  PTHDIR,AJID 
                            ! change if PTHDIR//PTHRDRD//AJID changes si
      CHARACTER*17 FULLPTH 
                                                                        
      CHARACTER*10 HFMODI,HFMODTU,HFMODTD,HFMRDR,hfmrad 
!                                                                       
!     Common block for analytic derivatives                             
!     -------------------------                                         
      COMMON /ADRPNM/ CDUM1,PTHODI,PTHODTU,PTHODTD 
      COMMON /ADRPTH/ PTHDIR,PTHRDRU,PTHRDRD,AJID 
      COMMON /ADRFRM/ HFMODI,HFMODTU,HFMODTD,HFMRDR 
      COMMON /ADRFIL/ KODFIL,kradtot,KTEMP,KFILAD,K_REFTRA,k_rddn_sfc 
      COMMON /IADFLG/ NSPCRT,imrgsav 
      COMMON /RADLAY/ PTHRAD,HFMRAD 
!     -------------------------                                         
!                                                                       
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
     &              NLTEFL,LNFIL4,LNGTH4                                
                                                                        
!---                                                                    
! this common must be changed if the FILHDR common is changed           
! it is only here for a dummy read to make sure nothing important       
! gets changed by accident                                              
                                                                        
      character*8     XIDj,HMOLIDj,YIDj 
      real*8          SECANTj,XALTZj 
      COMMON /DUMHDR/ XIDj(10),SECANTj,PAVEj,TAVEj,HMOLIDj(60),         &
     &    XALTZj(4),WKj(60),PZLj,PZUj,TZLj,TZUj,WBROADj,DVj,V1j,        &
     &    V2j,TBOUNDj,EMISIVj,FSCDIDj(17),NMOLj,LAYRSj,YI1j,            &
     &    YIDj(10),LSTWDF                                               
      dimension FILHDRj(2) 
      equivalence (filhdrj(1),xidj(1)) 
!---                                                                    
                                                                        
!                                                                       
!           123456789-123456789-123456789-123456789-123456789-1234567   
      DATA FILE1 /                                                      &
     &     '                                                         '/,&
     &     FILE2 /                                                      &
     &     '                                                         '/,&
     &     FILE3 /                                                      &
     &     '                                                         '/ &
     &     FILE4 /                                                      &
     &     '                                                         '/ 
      DATA CFORM / 'UNFORMATTED' / 
!                                                                       
                                                                        
      WRITE(IPR,910) LAYER,NLAYER 
      INQUIRE (UNIT=KODFIL,OPENED=OP) 
      IF (OP) CLOSE (KODFIL) 
      WRITE(FILE1,HFMODI) PTHODI,LAYER 
      OPEN(UNIT=KODFIL,FILE=FILE1,FORM=CFORM,STATUS='OLD') 
                                                                        
!      write (*,*) 'kodfil', kodfil, file1                              
!                                                                       
                                                                        
! open total rad files                                                  
!----                                                                   
                                                                        
      INQUIRE (UNIT=kradtot,OPENED=OP) 
      IF (OP) CLOSE (kradtot) 
                                                                        
! upwelling                                                             
                                                                        
      if (ipathl.eq.1 .and. layer.lt.nlayer) then 
         WRITE(FILE2,hfmrad) pthrad,LAYTOT 
         OPEN(UNIT=kradtot,FILE=FILE2,FORM=CFORM,STATUS='OLD') 
      endif 
                                                                        
! downwelling                                                           
                                                                        
      if (ipathl.eq.3 .and. layer.lt.nlayer) then 
         WRITE(FILE2,hfmrad) pthrad,LAYTOT 
         OPEN(UNIT=kradtot,FILE=FILE2,FORM=CFORM,STATUS='OLD') 
      endif 
!                                                                       
      INQUIRE (UNIT=KTEMP,OPENED=OP) 
      IF (OP) CLOSE (KTEMP) 
      OPEN(UNIT=KTEMP,FILE='AJ_mono',FORM=CFORM,STATUS='unknown') 
                                                                        
! analytic derivative file                                              
!  (not necessary for surface terms)                                    
                                                                        
      if (nspcrt.ge.0) then 
         INQUIRE (UNIT=KFILAD,OPENED=OP) 
         IF (OP) CLOSE (KFILAD) 
                                                                        
! downlooking/upwelling: ipathl = 1;  uplooking/downwelling: ipathl = 3 
         if (ipathl.eq.1) then 
            FULLPTH=PTHDIR//PTHRDRU//AJID 
            WRITE(FILE3,HFMRDR) FULLPTH,LAYER 
         else 
            FULLPTH=PTHDIR//PTHRDRD//AJID 
            WRITE(FILE3,HFMRDR) FULLPTH,LAYER 
         endif 
         OPEN(UNIT=KFILAD,FILE=FILE3,FORM=CFORM,                        &
     &        STATUS='UNKNOWN',ERR=881)                                 
                                                                        
! if upwelling, open appropriate downwelling file for merge             
! also open temporary file for surface terms                            
                                                                        
         IF (IPATHL.EQ.1) THEN 
                                                                        
         INQUIRE (UNIT=K_REFTRA,OPENED=OP) 
                                                                        
         IF (OP) close (K_REFTRA) 
                                                                        
         OPEN(UNIT=K_REFTRA,FILE='AJ_sfc_rfl.atm_tr', FORM=CFORM,STATUS=&
         'UNKNOWN')                                                     
                                                                        
         ENDIF 
                                                                        
! write file info to TAPE6                                              
         WRITE(IPR,900) FILE1,FILE2,FILE3 
      else 
         WRITE(IPR,902) FILE1,FILE2 
      endif 
!                                                                       
      RETURN 
                                                                        
! reach this if file error                                              
  881 continue 
      ierrmsg=881 
      goto 888 
                                                                        
  882 continue 
      ierrmsg=882 
      goto 888 
                                                                        
  888 continue 
                                                                        
      write(ipr,*) '**************************************' 
      write(ipr,*) '  ' 
      write(ipr,*) 'error opening file - check to see that' 
      write(ipr,*) ' AJ/ directory exists' 
      write(ipr,*) '  ' 
      write(ipr,*) '    see also xmerge.f error = ',ierrmsg 
      write(ipr,*) '**************************************' 
                                                                        
      write(*,*) '**************************************' 
      write(*,*) '  ' 
      write(*,*) 'error opening file - check to see that' 
      write(*,*) ' AJ/ directory exists' 
      write(*,*) '  ' 
      write(*,*) '    see also xmerge.f error = ',ierrmsg 
      write(*,*) '**************************************' 
                                                                        
!                                                                       
  900 FORMAT ('          Opened layer optical depth file:  ',A57,/,     &
     &        '    Opened accumulated optical depth file:  ',A57,/,     &
     &        '    Opened layer analytic derivative file:  ',A57)       
                                                                        
  902 FORMAT ('          Opened layer optical depth file:  ',A57,/,     &
     &        '    Opened accumulated optical depth file:  ',A57)       
  910 FORMAT ('LAYER ',I5,' OF ',I5,':') 
!                                                                       
end subroutine OPNDRV
!                                                                       
!     ----------------------------------------------------------------  
!                                                                       
      SUBROUTINE PRLNHD 
!                                                                       
      USE lblparams, ONLY: MXMOL
      IMPLICIT REAL*8           (V) 
!                                                                       
!     PRLNHD PRINTS OUT LINE FILE HEADER                                
!                                                                       
!      PARAMETER (MXMOL=39) 
!     -------------------------                                         
!                                                                       
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &              NLNGTH,KFILE,KPANEL,LINFIL,NDFLE,IAFIL,IEXFIL,      &
     &              NLTEFL,LNFIL4,LNGTH4                                
!                                                                       
      character*8      XID,       HMOLID,      YID 
      real*8               SECANT,       XALTZ 
!                                                                       
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
     &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
     &                EMISIV,FSCDID(17),NMOL,LAYHDR,YI1,YID(10),LSTWDF  
!                                                                       
      CHARACTER*8      HLINID,BMOLID,HID1,HLINHD 
!                                                                       
      integer *4 molcnt,mcntlc,                                         &
     &           mcntnl,linmol,                                         &
     &           lincnt,ilinlc,ilinnl,irec,irectl                       
!                                                                       
      COMMON /LINHDR/ HLINID(10),BMOLID(64),MOLCNT(64),MCNTLC(64),      &
     &                MCNTNL(64),SUMSTR(64),LINMOL,FLINLO,FLINHI,       &
     &                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1(2),LSTWDL   
      common /bufid2/ n_negepp(64),n_resetepp(64),xspace(4096),lstwdl2 
      common /eppinfo/ negepp_flag 
                                                                        
      real *4 sumstr,flinlo,flinhi 
      integer *4 lnfil 
      integer *4 negepp_flag,n_negepp,n_resetepp 
      real *4 xspace 
!                                                                       
!     LSTWD (LAST WORD) IS DUMMY, DOES NOT NEED TO BE COUNTED           
!                                                                       
      DIMENSION HLINHD(2),IWD(2) 
!                                                                       
      CHARACTER CHID10*8,CHARID*5,CHARDT*2,CHARI*1,CHTST*1 
      CHARACTER*1 CNEGEPP(8) 
      CHARACTER*6 CDUM,SPCRT 
!                                                                       
      EQUIVALENCE (HLINID(1),HLINHD(1),IWD(1)) 
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),              &
     &            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),              &
     &            (FSCDID(5),IEMIT) , (FSCDID(7),IPLOT),                &
     &            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                &
     &            (FSCDID(10),ITEST) , (FSCDID(11),IMRG),               &
     &            (FSCDID(16),LAYR1) , (FSCDID(17),NLAYHD)              
!                                                                       
      DATA CHARI / 'I'/ 
!                                                                       
      REWIND LINFIL 
                                                                        
      lnfil = linfil 
      negepp_flag = 0 
                                                                        
      read (lnfil,end=777)    HLINID,BMOLID,MOLCNT,MCNTLC,              &
     &                MCNTNL,SUMSTR,LINMOL,FLINLO,FLINHI,               &
     &                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1             
                                                                        
                                                                        
!     Test for negative values of ENERGY identified in lnfl             
!     and read in second header for line information, if needed         
                                                                        
      READ (HLINID(7),950) CNEGEPP 
      IF (CNEGEPP(8).eq.'^') THEN 
         negepp_flag = 1 
         read (lnfil) n_negepp,n_resetepp,xspace 
         endif 
!                                                                       
         go to 5 
!                                                                       
  777    STOP 'LAYER; TAPE3 DOES NOT EXIST' 
!                                                                       
    5    continue 
!                                                                       
         DO 10 M = 1, LINMOL 
            HMOLID(M) = BMOLID(M) 
   10    END DO 
         WRITE (IPR,900) 
         WRITE (IPR,905) HLINID,HID1 
                                                                        
!     Output header information regarding lines; if negative values of  
!     ENERGY were identified in lnfl, output extra header information   
                                                                        
         if (CNEGEPP(8).eq.'^') THEN 
         WRITE (IPR,960) 
         WRITE (IPR,965) (BMOLID(I),MOLCNT(I),MCNTLC(I),MCNTNL(I),      &
         N_NEGEPP(I),N_RESETEPP(I), SUMSTR(I),I=1,LINMOL)               
         else 
         WRITE (IPR,910) 
         WRITE (IPR,915) (BMOLID(I),MOLCNT(I),MCNTLC(I),MCNTNL(I),      &
         SUMSTR(I),I=1,LINMOL)                                          
         endif 
!                                                                       
         WRITE (IPR,920) FLINLO,FLINHI,LINCNT 
!                                                                       
!     When calculating derivative, check make sure the                  
!     appropriate molecule is included in the linefile.                 
!     If not, then stop and issue message.                              
!                                                                       
!     CHECK HEADER FOR FLAG INDICATING COMPATIBILITY WITH ISOTOPES      
!                                                                       
   30    WRITE (CHID10,925) HLINID(10) 
         READ (CHID10,930) CHARID,CHARDT,CHTST 
         IF (CHTST.NE.CHARI) THEN 
            WRITE (IPR,935) CHARID,CHARDT,CHTST 
            STOP ' PRLNHD - NO ISOTOPE INFO ON LINFIL ' 
         ENDIF 
!                                                                       
         RETURN 
!                                                                       
  900 FORMAT ('0'/'0',20X,'   LINE FILE INFORMATION ') 
  905 FORMAT ('0',10A8,2X,2(1X,A8,1X)) 
  910 FORMAT ('0',/,23X,'COUPLED',4X,'NLTE',3X,'SUM LBLRTM ',/,7X,      &
     &        'MOL',5X,'LINES',4X,'LINES',4X,'LINES',4X,'STRENGTHS',/)  
  915 FORMAT (' ',4X,A6,' = ',I6,3X,I6,3X,I6,2X,1PE12.4,0P) 
  920 FORMAT (/,'0 LOWEST LINE = ',F10.3,5X,'  HIGHEST LINE = ',F10.3,  &
     &        5X,' TOTAL NUMBER OF LINES =',I8)                         
  925 FORMAT (A8) 
  930 FORMAT (A5,A2,A1) 
  935 FORMAT (3(/),10X,'LINEFILE PROGRAM: ',A5,3X,'VERSION: ',A2,A1,    &
     &        3(/),3X,52('*'),/,3X,'* THE LINEFILE (TAPE3) IS NOT ',    &
     &        'COMPATIBLE WITH THIS *',/,3X,'* VERSION OF LBLRTM .',    &
     &        '  ISOTOPIC INFORMATION (FROM  *',/,3X,'* HITRAN) ',      &
     &        'MUST BE PRESERVED ON TAPE3.  USE A TAPE3 *',/,3X,        &
     &        '* CREATED WITH THE 91I OR LATER VERSION OF LNFL.   *',   &
     &        /,3X,52('*'))                                             
  940 FORMAT (' Molecule to be retrieved: ',A6,' not in linefile.',/,   &
     &        ' Molecules in linefile: ')                               
  945 FORMAT (24X,A6) 
  950 FORMAT (8a1) 
  960 FORMAT ('0',/,23X,'COUPLED',4X,'NLTE',3X,'NEGATIVE',3X,           &
     &        'RESET',4X,'SUM LBLRTM',/,7X,'MOL',5X,'LINES',4X,         &
     &        'LINES',4X,'LINES',6X,'EPP',6X,'EPP',                     &
     &        6X,'STRENGTHS',/)                                         
  965 FORMAT (' ',4X,A6,' = ',I6,                                       &
     &        3X,I6,3X,I6,3X,I6,3X,i6,3X,1PE12.4)                       
!                                                                       
end subroutine PRLNHD
!                                                                       
!     -------------------------------------------------------------     
!                                                                       
      SUBROUTINE EXPINT (X,X1,X2,A) 
!                                                                       
!********************************************************************** 
!     THIS SUBROUTINE EXPONENTIALLY INTERPOLATES X1 AND X2 TO X BY      
!     THE FACTOR A                                                      
!********************************************************************** 
!                                                                       
      IF (X1.EQ.0.0.OR.X2.EQ.0.0) GO TO 10 
      X = X1*(X2/X1)**A 
!                                                                       
      RETURN 
!                                                                       
   10 X = X1+(X2-X1)*A 
!                                                                       
      RETURN 
!                                                                       
end subroutine EXPINT
!                                                                       
!     -------------------------------------------------------------     
!                                                                       
      SUBROUTINE XLAYER (MPTS,NPTS,LFILE,MFILE,NFILE) 
!                                                                       
      USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
                           MXMOL, MX_XS, MXTRAC, MXSPC, IPTS,           &
                           IPTS2
      IMPLICIT REAL*8           (V) 
!                                                                       
!********************************************************************** 
!     XLAYER CONTROLS LAYER BY LAYER CALCULATION                        
!********************************************************************** 
!                                                                       
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &     MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,mx_xs=38,MXTRAC=22, &
!     &     mxspc=5)                                                     
!                                                                       
      CHARACTER*55 CDUM1,PTHODI,PTHODTU,PTHODTD 
      CHARACTER*55 PTHRAD,PATH1,PATH2 
      CHARACTER*11 PTHRDRU,PTHRDRD 
      CHARACTER*3  PTHDIR,AJID 
      CHARACTER*10 HFMODI,HFMODTU,HFMODTD,HFMRDR,HFMRAD,HFORM1,HFORM2 
!                                                                       
!     -------------------------                                         
      common /cmol_nam/ cmol(mxmol),cspc(mxspc) 
      CHARACTER*6  CMOL,CSPC,SPCRT 
!                                                                       
!     -------------------------                                         
!                                                                       
!     Common block for analytic derivative                              
!     -------------------------                                         
      COMMON /IADFLG/ NSPCRT,IMRGSAV 
      COMMON /ADRPNM/ CDUM1,PTHODI,PTHODTU,PTHODTD 
      COMMON /ADRPTH/ PTHDIR,PTHRDRU,PTHRDRD,AJID 
      COMMON /ADRFRM/ HFMODI,HFMODTU,HFMODTD,HFMRDR 
      COMMON /ADRFIL/ KODFIL,kradtot,KTEMP,KFILAD,K_REFTRA,k_rddn_sfc 
! note: from continuum module                                           
!          ipts  = same dimension as ABSRB                              
!          ipts2 = same dimension as C                                  
!      parameter (ipts=5050,ipts2=6000) 
      common /CDERIV/ icflg,idum,v1absc,v2absc,dvabsc,nptabsc,delT_pert,&
     &    dqh2oC(ipts),dTh2oC(ipts),dUh2o                               
      logical op 
      common /dlaydlev/ilevdx,imoldx,iup_dn,                            &
     &    dxdL(mxlay,0:mxmol),dxdU(mxlay,0:mxmol)                       
!---                                                                    
! this common must be changed if the FILHDR common is changed           
! it is only here for a dummy read to make sure nothing important       
! gets changed by accident                                              
      character*8     XIDj,HMOLIDj,YIDj 
      real*8          SECANTj,XALTZj 
      COMMON /DUMHDR/ XIDj(10),SECANTj,PAVEj,TAVEj,HMOLIDj(60),         &
     &    XALTZj(4),WKj(60),PZLj,PZUj,TZLj,TZUj,WBROADj,DVj,V1j,        &
     &    V2j,TBOUNDj,EMISIVj,FSCDIDj(17),NMOLj,LAYRSj,YI1j,            &
     &    YIDj(10),LSTWDFj                                              
      dimension FILHDRj(2) 
      equivalence (filhdrj(1),xidj(1)) 
!---                                                                    
                                                                        
                                                                        
!                                                                       
!     -------------------------                                         
!     Common blocks for layer radiances                                 
!     -------------------------                                         
!                                                                       
      COMMON /RADLAY/ PTHRAD,HFMRAD 
!     -------------------------                                         
!                                                                       
      COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,NOZERO,NP,H1F,H2F,       &
     &                ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH,       &
     &                XVBAR, HMINF,PHIF,IERRF,HSPACE                    
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
     &              AVFIX,LAYER,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,      &
     &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
     &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
     &              EXTID(10)                                           
      CHARACTER*8  EXTID 

      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &              NLNGTH,KFILE,KPANEL,LINFIL,NFILA,IAFIL,IEXFIL,      &
     &              NLTEFL,LNFIL4,LNGTH4                                
      COMMON /MSCONS/ AIRMAS(MXLAY),TGRND,SEMIS(3),HMINMS,HMAXMS,       &
     &                MSFLAG,                                           &
     &                MSWIT,IODFIL,MSTGLE                               
      COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),          &
     &                MSPNL1(MXLAY),MSLAY1,ISFILE,JSFILE,KSFILE,        &
     &                LSFILE,MSFILE,IEFILE,JEFILE,KEFILE                
      COMMON /SCSHAP/ HWFS,DXFS,NFS,NFMAXS 
      COMMON /CMSHAP/ HWF1,DXF1,NX1,N1MAX,HWF2,DXF2,NX2,N2MAX,          &
     &                HWF3,DXF3,NX3,N3MAX                               
!                                                                       
      COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(mx_xs),XAMNT(mx_xs,MXLAY) 
!                                                                       
!     COMMON /MLTSCT/ TAUGAS(2410),FUPC(2410),RUPC(2410)                
!                                                                       
      COMMON /LASIV/ VLAS,ILAS 
!                                                                       
      character*8      XID,       HMOLID,      YID 
      real*8               SECANT,       XALTZ 
!                                                                       
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
     &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
     &                EMISIV,FSCDID(17),NMOL,LAYHDR,YI1,YID(10),LSTWDF  
!                                                                       
      character*8      XI1,       HMOLI1,      Y1D 
      real*8               SECAN1,       XALT1 
!                                                                       
      COMMON /FILHD1/ XI1(10),SECAN1,PAV1,TAV1,HMOLI1(60),XALT1(4),     &
     &                W1(60),PDL,PDU,TDL,TDU,Wbrd1 ,D1 ,VD1,VD2,TBOUN1, &
     &                EMISI1,FSCDI1(17),NMO1,LAYHD1,YD1,Y1D(10),LSTWDD  
      COMMON /XSECTR/ V1FX(5,MX_XS),V2FX(5,MX_XS),DVFX(5,MX_XS),        &
     &                WXM(MX_XS),NTEMPF(5,MX_XS),NSPECR(MX_XS),         &
     &                IXFORM(5,MX_XS),XSMASS(MX_XS),XDOPLR(5,MX_XS),    &
     &                NUMXS,IXSBIN                                      
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,pad_3,&
     &                angle_path,secant_diffuse,secant_path,diffuse_fac 
!                                                                       
      COMMON /IODFLG/ DVOUT 
!                                                                       
      character*1 surf_refl,surf_refl_sav,h_blank 
      character*3 pad_3 
!                                                                       
      integer *4 molcnt,mcntlc,                                         &
     &           mcntnl,linmol,                                         &
     &           lincnt,ilinlc,ilinnl,irec,irectl                       
!                                                                       
      real *4    sumstr,flinlo,flinhi 
!                                                                       
      CHARACTER*8     HLINID,BMOLID,HID1 
!                                                                       
      COMMON /LINHDR/ HLINID(10),BMOLID(64),MOLCNT(64),MCNTLC(64),      &
     &                MCNTNL(64),SUMSTR(64),LINMOL,FLINLO,FLINHI,       &
     &                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1(2),LSTWDL   
!                                                                       
      DIMENSION FILDUM(2),FILDU1(2) 
      DIMENSION NTAN(160) 
!                                                                       
      EQUIVALENCE (XID(1),FILDUM(1)) , (XI1(1),FILDU1(1)) 
                                                                        
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),              &
     &            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),              &
     &            (FSCDID(5),IEMIT)  , (FSCDID(6),ISCNHD),              &
     &            (FSCDID(7),IPLOT)  , (FSCDID(8),IPATHL),              &
     &            (FSCDID(9),JRAD)   , (FSCDID(10),ITEST),              &
     &            (FSCDID(11),IMRG)  , (FSCDID(12),SCNID),              &
     &            (FSCDID(13),HWHM)  , (FSCDID(14),IDABS),              &
     &            (FSCDID(15),IATM)  , (FSCDID(16),LAYR1),              &
     &            (FSCDID(17),NLAYHD),                                  &
     &            (Y1D(8),LH2SAV)    , (Y1D(9),LH1SAV)!,                &
!    &            (Y1D(10),LTNSAV, dv_lbl)                              
                                                                        
      equivalence (FSCDI1(8),IPTHD1),(FSCDI1(17),NLAYD1),              &
     &            (Y1D(10),LTNSAV, dv_lbl1)! ,                           &
!    &            (Y1D(8),LH2SAV),  (Y1D(9),LH1SAV)                     
!                                                                       
      DATA I_10/10/ 
      data h_blank /' '/ 
!                                                                       
      DV = 0. 
!                                                                       
!     IF  IANT.EQ. 1  THEN POSTERIOR MERGE                              
!     IF  IANT.EQ. 0  THEN NORMAL MERGE                                 
!     IF  IANT.EQ.-1  THEN ANTERIOR MERGE                               
!                                                                       
      IANT = 0 
      JPATHL = 0 
!                                                                       
      IBUF = 1 
      JEMIT = 0 
!                                                                       
      LAYER = 0 
      P0 = 1013.25 
      TEMP0 = 296. 
!                                                                       
!********************************************************************** 
!                                                                       
!     IMRG = 0 AND IMRG > 10 ONLY LAST LAYER OPTICAL DEPTH ON KFILE     
!       (EXCEPT 12, 22 AND 32 WHICH ARE SEQUENTIAL FROM PRESTORE)       
!                                                                       
!     1 < IMRG < 10 OPTICAL DEPTHS ARE SEQUENTIAL ON KFILE BY LAYER     
!                                                                       
!     IMRG = 1 Optical Depths are stored on different files by layer    
!                                                                       
!     IMRG = 10 OPTICAL DEPTHS are stored on multiple KFILEs by layer   
!               and accumulated Optical Depths are calculated by layer  
!     IMRG = 40-46 RADIANCE and TRANSMITTANCE are calculated from       
!               multiple optical depth files by layer                   
!                                                                       
!********************************************************************** 
!                                                                       
!      IMRG= 0 STANDARD MERGE                                           
!      IMRG= 1 NO MERGE TAKES PLACE                                     
!      IMRG= 2 SEQUENTIAL OPTICAL DEPTHS FROM KFILE MERGED ONTO MFILE   
!      IMRG= A SEQUENTIAL OPTICAL DEPTHS FROM KFILE USED TO CALCULATE   
!                  MONOCHROMATIC RESULTS WHICH ARE MERGED ONTO MFILE    
!                  FORCES IPATHL = 1 (SAME AS IMRG=12)                  
!      IMRG= B SEQUENTIAL OPTICAL DEPTHS FROM KFILE USED TO CALCULATE   
!                  MONOCHROMATIC RESULTS WHICH ARE MERGED ONTO MFILE    
!                  FORCES IPATHL = 2 (SAME AS IMRG=22)                  
!      IMRG= C SEQUENTIAL OPTICAL DEPTHS FROM KFILE USED TO CALCULATE   
!                  MONOCHROMATIC RESULTS WHICH ARE MERGED ONTO MFILE    
!                  FORCES IPATHL = 3 (SAME AS IMRG=32)                  
!                                                                       
!      IMRG= 9 MULTIPLE RUNS EMISSION FOR AEROSOLS FROM PRESTORE        
!                                                                       
!        ****  WEIGHTING FUNCTIONS - MONOCHROMATIC  ****                
!                                                                       
!      IMRG= 3 --- SPACE TO GROUND                                      
!      IMRG= 4 --- GROUND TO SPACE                                      
!      IMRG= 5 --- SPACE TO GROUND FROM PRESTORED OPTICAL DEPTHS        
!      IMRG= 6 --- GROUND TO SPACE FROM PRESTORED OPTICAL DEPTHS        
!      IMRG= 7 --- H1 THROUGH H(TAN) TO H2                              
!      IMRG= 8 --- H1 THROUGH H(TAN) TO H2 FROM PRESTORED OPTICAL       
!                  DEPTHS                                               
!                                                                       
!        ****  WEIGHTING FUNCTIONS - SCANNED ****                       
!                                                                       
!      IMRG=13 --- SPACE TO GROUND                                      
!      IMRG=14 --- GROUND TO SPACE                                      
!      IMRG=15 --- SPACE TO GROUND FROM PRESTORED OPTICAL DEPTHS        
!      IMRG=16 --- GROUND TO SPACE FROM PRESTORED OPTICAL DEPTHS        
!      IMRG=17 --- H1 THROUGH H(TAN) TO H2                              
!      IMRG=18 --- H1 THROUGH H(TAN) TO H2 FROM PRESTORED OPTICAL       
!                  DEPTHS                                               
!                                                                       
!        ****  WEIGHTING FUNCTIONS - FILTERED ****                      
!                                                                       
!      IMRG=23 --- SPACE TO GROUND                                      
!      IMRG=24 --- GROUND TO SPACE                                      
!      IMRG=25 --- SPACE TO GROUND FROM PRESTORED OPTICAL DEPTHS        
!      IMRG=26 --- GROUND TO SPACE FROM PRESTORED OPTICAL DEPTHS        
!      IMRG=27 --- H1 THROUGH H(TAN) TO H2                              
!      IMRG=28 --- H1 THROUGH H(TAN) TO H2 FROM PRESTORED OPTICAL       
!                  DEPTHS                                               
!                                                                       
!        ****  FLUX CALCULATIONS -SCANNED ****                          
!                                                                       
!      IMRG=35 --- SPACE TO GROUND MERGE FROM PRESTORED OPTICAL DEPTHS  
!      IMRG=36 --- GROUND TO SPACE MERGE FROM PRESTORED OPTICAL DEPTHS  
!                                                                       
!        ****  RADIANCE/DERIVATIVE CALCULATIONS  ****                   
!                                                                       
!      IMRG=40 --- Upwelling radiance from prestored optical depths,    
!                  monochromatic                                        
!      IMRG=41 --- Downwelling radiance from prestored optical depths,  
!                  monochromatic                                        
!                                                                       
!      IMRG=42 --- Upwelling radiance from prestored optical depths,    
!                  scanned                                              
!      IMRG=43 --- Downwelling radiance from prestored optical depths,  
!                  scanned                                              
!                                                                       
!  Note for IMRG=40,41:  Monochromatic radiance & derivative calculation
!  Note for IMRG=42,43:  Only the derivative calculations are scanned   
!                                                                       
!        ****  FLUX CALCULATIONS (Layer Radiance) -MONOCHROMATIC ****   
!                                                                       
!      IMRG=45 --- Space to ground merge from prestored optical depths  
!      IMRG=46 --- Ground to space merge from prestored optical depths  
!                                                                       
!      WEIGHTING FUNCTION RESULTS ARE ON NFILE SEPARATED BY             
!      INTERNAL 'EOF'                                                   
!                                                                       
!********************************************************************** 
!                                                                       
!     ---------------------                                             
!     Special Merge Options                                             
!     ---------------------                                             
!                                                                       
!     If IMRG = 1, then calculate optical depths in standard fashion,   
!     but output to different files for each layer.                     
!                                                                       
!                                            SPECIAL CASE -> IMRG=1     
!                                                                       
!     Start loop over layers                                            
!                                                                       
      IF (IMRG.EQ.1) THEN 
!                                                                       
!        -----------------------------                                  
!        Initial call to OPPATH, which calls PATH                       
!        -----------------------------                                  
!                                                                       
         LAYHDR = LAYER 
         CALL OPPATH 
         IF (IHIRAC.EQ.0) RETURN 
!                                                                       
!        -----------------------------                                  
!        Begin loop over layers                                         
!        -----------------------------                                  
!                                                                       
    1    LAYER = LAYER+1 
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
!                                                                       
!     ---------------------                                             
!                                                                       
!     If IMRG = 10, then calculate optical depths on multiple           
!     files, and calculate accumulated optical depths to                
!     multiple files. Default names for input files KFILE "ODDV##"      
!     and output files MFILE "ODMG##" are used (layer number ##).       
!     Spacing of optical depths for all layers should be equal.         
!                                                                       
!                                            SPECIAL CASE -> IMRG=10    
!                                                                       
!     First portion mimicks IMRG = 1 procedure, and then when           
!     LAYER = NLAYER, merge takes place.                                
!                                                                       
!                                                                       
      IF (IMRG.EQ.10) THEN 
!                                                                       
!        -----------------------------                                  
!        Initial call to OPPATH, which calls PATH                       
!        -----------------------------                                  
!                                                                       
         LAYHDR = LAYER 
         CALL OPPATH 
         IF (IHIRAC.EQ.0) RETURN 
!                                                                       
!        -----------------------------                                  
!        Begin loop over layers                                         
!        -----------------------------                                  
!                                                                       
    4    LAYER = LAYER+1 
         LAYHDR = LAYER 
         CALL OPPATH 
         NLAYHD = NLAYER 
         CALL OPDPTH (MPTS) 
         CALL ENDFIL (KFILE) 
         CLOSE (KFILE) 
         REWIND MFILE 
         REWIND LFILE 
         IF (LAYER.EQ.NLAYER) THEN 
                                                                        
            if ((ipathl.ne.1).and.(ipathl.ne.3)) then 
            STOP 'XLAYER: IPATHL NOT VALID' 
            endif 
                                                                        
!           First, copy the farthest layer optical depths to the        
!           pathname for the total optical depths up to the             
!           first layer. Then, add the sum of the previous              
!           L layers to the (L+1)'st layer.  Output procedure           
!           (title written here, rest written in OPNMRG).               
                                                                        
!           Note: do downlooking first and then uplooking               
!                 do both up and down regardless of ipathl              
                                                                        
! downlooking/upwelling (ipathl = 1)                                    
            PATH1 = PTHODI 
            PATH2 = PTHODTD 
            HFORM1 = HFMODI 
            HFORM2 = HFMODTD 
            WRITE(IPR,935) 
            CALL OPNMRG(LFILE,PATH1,NLAYER,HFORM1,PATH1,NLAYER, HFORM1, &
            MFILE,PATH2,HFORM2)                                         
            CALL COPYFL(NPTS,LFILE,MFILE) 
            CALL ENDFIL (MFILE) 
            CLOSE(MFILE) 
            CLOSE(LFILE) 
            CLOSE(KFILE) 
            DO 5 L = 2,NLAYER 
               CALL OPNMRG(LFILE,PATH2,NLAYER-L+2,HFORM2,PATH1, NLAYER- &
               L+1,HFORM1,MFILE,PATH2,HFORM2)                           
               CALL XMERGE (NPTS,LFILE,MFILE,JPATHL) 
               CLOSE(MFILE) 
               CLOSE(LFILE) 
               CLOSE(KFILE) 
    5       CONTINUE 
                                                                        
! uplooking/downwelling (ipathl = 3)                                    
            PATH1 = PTHODI 
            PATH2 = PTHODTU 
            HFORM1 = HFMODI 
            HFORM2 = HFMODTU 
            WRITE(IPR,936) 
            CALL OPNMRG(LFILE,PATH1,1,HFORM1,PATH1,1,HFORM1, MFILE,     &
            PATH2,HFORM2)                                               
            CALL COPYFL(NPTS,LFILE,MFILE) 
            CALL ENDFIL (MFILE) 
            CLOSE(MFILE) 
            CLOSE(LFILE) 
            CLOSE(KFILE) 
            DO 7 L = 2,NLAYER 
               CALL OPNMRG(LFILE,PATH2,L-1,HFORM2,PATH1, L,HFORM1,MFILE,&
               PATH2,HFORM2)                                            
               CALL XMERGE (NPTS,LFILE,MFILE,JPATHL) 
               CLOSE(MFILE) 
               CLOSE(LFILE) 
               CLOSE(KFILE) 
    7       CONTINUE 
                                                                        
! done with total optical depth calculation                             
            RETURN 
         ENDIF 
         GO TO 4 
!                                                                       
!        -----------------------------                                  
!        End loop over layers                                           
!        -----------------------------                                  
!                                                                       
      ENDIF 
!                                                                       
!     ---------------------                                             
!                                                                       
!     For IMRG = 35,36,40,41,42,43,45,46 (those options which           
!     use precalculated layer optical depths stored on different        
!     files for radiative transfer), read in the pathname of            
!     the layer optical depths and determine format for the             
!     addition of the layer number suffix.                              
!                                                                       
      IF (IMRG.GE.35 .and. iemit.ne.3) THEN 
                                                                        
!     Get the pathname for the optical depth files                      
         READ (IRD,945) PATH1, laytot 
         CALL QNTIFY(PATH1,HFORM1) 
                                                                        
!                    nlayr, layr,                                       
         CALL OPNODF( 1, 1,PATH1,HFORM1,-iemit) 
                                                                        
!     Get the maximum number of layers from the optical depth file for l
         REWIND kodfil 
         CALL BUFIN (kodfil,KEOF,FILDU1(1),NFHDRF) 
                                                                        
         close (kodfil) 
                                                                        
         nlayer = nlayd1 
                                                                        
      ENDIF 
!                                                                       
!***********************************************************************
!***********************************************************************
                                                                        
! radiance merge from space to surface using precalculated ODint files  
                                                                        
      IF ((IMRG.eq.40.or.IMRG.eq.42) .AND. (IEMIT.EQ.1)) THEN 
                                                                        
         iemit_sav = iemit 
         iemit = 1 
!                                                                       
         lh1 = nlayer 
         lh2 = 1 
                                                                        
         tmpbndsav = tmpbnd 
         tmpbnd = 0. 
                                                                        
         ipathl_sav = ipathl 
         jpathl_sav = jpathl 
                                                                        
         surf_refl_sav = surf_refl 
         surf_refl = h_blank 
                                                                        
!     ipathl is set to 31 for algorithm in which the RT (downwelling) is
!            direction of the loop overlayers,  nlayer > 1.             
                                                                        
         ipathl = 31 
         jpathl = ipathl 
                                                                        
!     set up format for file name                                       
         CALL QNTIFY(PATH1,HFORM1) 
                                                                        
!     set path name for rad down files                                  
         PTHRAD = 'RDDNlayer_' 
         CALL QNTIFY(PTHRAD,HFMRAD) 
                                                                        
         kfile = kodfil 
                                                                        
         IF (2*(NLAYER/2).eq.NLAYER) then 
            MSTOR = MFILE 
            MFILE = LFILE 
            LFILE = MSTOR 
            endif 
                                                                        
            LAYER = nlayer+1 
                                                                        
!***************Loop over Layers ************                           
  705       LAYER = LAYER-1 
                                                                        
            REWIND MFILE 
            REWIND LFILE 
            WRITE (IPR,900) 
            WRITE (IPR,905) 
                                                                        
!     open appropriate optical depth file                               
            CALL OPNODF(NLAYER,LAYER,PATH1,HFORM1,-iemit) 
            WRITE (IPR,905) 
                                                                        
            CALL XMERGE (NPTS,LFILE,MFILE,JPATHL) 
                                                                        
            CALL OPNRAD(nfile,NLAYER,LAYER) 
                                                                        
            rewind mfile 
            CALL COPYFL (NPTS,MFILE,NFILE) 
            close (nfile) 
                                                                        
            IF (LAYER.EQ.lh2) go to 708 
                                                                        
            MSTOR = MFILE 
            MFILE = LFILE 
            LFILE = MSTOR 
!                                                                       
!        END OF LOOP OVER LAYERS                                        
!                                                                       
            GO TO 705 
                                                                        
  708       continue 
                                                                        
!        RDDN files have been created                                   
                                                                        
            tmpbnd = tmpbndsav 
            surf_refl = surf_refl_sav 
                                                                        
         ENDIF 
                                                                        
!***********************************************************************
!***********************************************************************
                                                                        
                                                                        
! radiance merge from surface to space using precalculated ODint files  
                                                                        
! downwelling radiance at the surface must first be calculated for refle
                                                                        
!     ---------------------                                             
!                                     SPECIAL CASE -> IMRG=41-43, IEMIT=
!                                                                       
         IF (((IMRG.eq.41.or.IMRG.eq.43)) .AND. (IEMIT.EQ.1)) THEN 
                                                                        
!     First compute the downwelling radiance at the lower boundary, e.g.
                                                                        
            iemit_sav = iemit 
            iemit = 1 
!                                                                       
            lh1 = 1 
            lh2 = nlayer 
                                                                        
            tmpbndsav = tmpbnd 
                                                                        
            ipathl_sav = ipathl 
            jpathl_sav = jpathl 
                                                                        
            surf_refl_sav = surf_refl 
                                                                        
            ipathl = 3 
            jpathl = ipathl 
                                                                        
!     set up format for file name                                       
            CALL QNTIFY(PATH1,HFORM1) 
                                                                        
!     set path name for rad up files                                    
            PTHRAD = 'RDUPlayer_' 
            CALL QNTIFY(PTHRAD,HFMRAD) 
                                                                        
            kfile = kodfil 
                                                                        
            IF (2*(NLAYER/2).eq.NLAYER) then 
               MSTOR = MFILE 
               MFILE = LFILE 
               LFILE = MSTOR 
               endif 
                                                                        
               LAYER = 0 
                                                                        
!***************Loop over Layers ************                           
  725          LAYER = LAYER+1 
                                                                        
               REWIND MFILE 
               REWIND LFILE 
               WRITE (IPR,900) 
               WRITE (IPR,905) 
                                                                        
!     open appropriate optical depth file                               
               CALL OPNODF(NLAYER,LAYER,PATH1,HFORM1,-iemit) 
               WRITE (IPR,905) 
                                                                        
               CALL XMERGE (NPTS,LFILE,MFILE,JPATHL) 
                                                                        
               IF (LAYER.EQ.lh2) go to 728 
                                                                        
               MSTOR = MFILE 
               MFILE = LFILE 
               LFILE = MSTOR 
!                                                                       
!        END OF LOOP OVER LAYERS                                        
!                                                                       
               GO TO 725 
                                                                        
  728          continue 
                                                                        
               k_rddn = 27 
                                                                        
               OPEN(UNIT=k_rddn,FILE='RDDN_sfc',FORM='unformatted',     &
               STATUS='UNKNOWN')                                        
                                                                        
               rewind mfile 
               CALL COPYFL (NPTS,MFILE,k_rddn) 
               close (k_rddn) 
!_______________________________________________________________        
                                                                        
!     Now create upwelling radiance files                               
                                                                        
               iemit = 1 
               ipathl = 1 
!                                                                       
               lh1 = 1 
               lh2 = nlayer 
                                                                        
               tmpbndsav = tmpbnd 
                                                                        
               ipathl_sav = ipathl1 
               jpathl_sav = jpathl 
                                                                        
               surf_refl_sav = surf_refl 
                                                                        
               jpathl = ipathl 
                                                                        
!     set up format for file name                                       
               CALL QNTIFY(PATH1,HFORM1) 
                                                                        
!     set path name for rad up files                                    
               PTHRAD = 'RDUPlayer_' 
               CALL QNTIFY(PTHRAD,HFMRAD) 
                                                                        
               kfile = kodfil 
                                                                        
               IF (2*(NLAYER/2).eq.NLAYER) then 
                  MSTOR = MFILE 
                  MFILE = LFILE 
                  LFILE = MSTOR 
                  endif 
                                                                        
                  LAYER = 0 
                                                                        
!***************Loop over Layers ************                           
  735             LAYER = LAYER+1 
                                                                        
                  REWIND MFILE 
                  REWIND LFILE 
                  WRITE (IPR,900) 
                  WRITE (IPR,905) 
                                                                        
!     open appropriate optical depth file                               
                  CALL OPNODF(NLAYER,LAYER,PATH1,HFORM1,-iemit) 
                  WRITE (IPR,905) 
                                                                        
                  CALL XMERGE (NPTS,LFILE,MFILE,JPATHL) 
                                                                        
                  CALL OPNRAD(nfile,NLAYER,LAYER) 
                                                                        
                  rewind mfile 
                  CALL COPYFL (NPTS,MFILE,NFILE) 
                  close (nfile) 
                                                                        
                  IF (LAYER.EQ.lh2) go to 738 
                                                                        
                  MSTOR = MFILE 
                  MFILE = LFILE 
                  LFILE = MSTOR 
!                                                                       
!        END OF LOOP OVER LAYERS                                        
!                                                                       
                  GO TO 735 
                                                                        
  738             continue 
                                                                        
!        RDUP files have been created                                   
                                                                        
                  tmpbnd = tmpbndsav 
                  surf_refl = surf_refl_sav 
                                                                        
               ENDIF 
                                                                        
!***********************************************************************
!                                                                       
!                                                                       
!     ****************** ANALYTIC DERIVATIVE ***********************    
!                                                                       
!     Assign name to number of the species selected.                    
!     If the species is a molecule, then assign the appropriate         
!     molecule name from CMOL.  If the species is something other       
!     than a molecule (e.g., layer temperature, surface temperature,    
!     etc.), then assign the appropriate species name from CSPC.        
!                                                                       
               IF (IEMIT.EQ.3) THEN 
                                                                        
                  INQUIRE (UNIT=k_rddn_sfc,OPENED=OP) 
                  IF (OP) CLOSE (k_rddn_sfc) 
                  OPEN(UNIT=k_rddn_sfc,FILE='RDDNlayer_001', FORM=      &
                  'unformatted',STATUS='unknown')                       
                                                                        
                  READ (IRD,1015) NSPCRT 
                  IF ((NSPCRT.GT.0).AND.(NSPCRT.LE.MXMOL)) THEN 
                     SPCRT = CMOL(NSPCRT) 
                  ELSEIF ((NSPCRT.LE.0).AND.(ABS(NSPCRT).LE.MXSPC))     &
                  THEN                                                  
                     SPCRT = CSPC(abs(NSPCRT)+1) 
                  ELSE 
                     write(ipr,*) 'FATAL ERROR ON NSPCRT SPECIFICATION' 
            write(ipr,*) 'value must be between ',                      &
     &           (-mxspc),' and ',mxmol                                 
                     STOP 
                  ENDIF 
                                                                        
                                                                        
!     write information on variable for analytic derivative             
                                                                        
         write (ipr,*) '                                              ' 
                  write (ipr,*)                                         &
                  '**********************************************'      
         write (ipr,*) '****** Analytic Derivative *******************' 
         write (ipr,*) '******                           *************' 
         write (ipr,*) '****** variable number =    ',nspcrt,           &
     &                                                '  *************' 
         write (ipr,*) '******   variable name = ',spcrt,               &
     &                                                '  *************' 
                  write (ipr,*)                                         &
                  '***********************************************'     
         write (ipr,*) '                                              ' 
!                                                                       
!     check that molec ule is in the line file:                         
                                                                        
                  IF ((NSPCRT.GT.0).and.(NSPCRT.LE.MXMOL)) THEN 
                     DO M = 1,LINMOL 
                     if (m.eq.nspcrt .and. molcnt(m).gt.0) go to 739 
                     ENDDO 
                                                                        
            write (ipr,*)  'Molecule to be retrieved not in line file' 
            WRITE(IPR,*)           'nspcrt,  spcrt,  linmol' 
                     WRITE(IPR,'(i5,a8,i5)') nspcrt, spcrt, linmol 
                     STOP 'Molecule to be retrieved not in line file' 
                  ENDIF 
!                                                                       
  739             continue 
!                                                                       
                  if ((nspcrt.ge.0).and.(nspcrt.le.9)) then 
                                                    ! need to change if 
                  write(ajid,'(a1,i1,a1)') '0',nspcrt,'_' 
                  else 
                                             ! need to change if ajid si
                  write(ajid,'(i2,a1)') nspcrt,'_' 
                  endif 
                                                                        
                  icflg=nspcrt 
                                                                        
               ENDIF 
                                                                        
!     *********     Do the surface derivatives   ****************       
                                                                        
!     dL/demis, dL/drefl and dL/dTsfc                                   
                                                                        
               IF (((IMRG.eq.41.or.IMRG.eq.43)) .AND. (IEMIT.EQ.3)      &
               .and. icflg.eq.-1) then                                  
                                                                        
                  INQUIRE (UNIT=k_rddn_sfc,OPENED=OP) 
                  IF (OP) CLOSE (k_rddn_sfc) 
                  OPEN(UNIT=k_rddn_sfc,FILE='RDDNlayer_001', FORM=      &
                  'unformatted',STATUS='unknown')                       
                                                                        
                  call sfcderiv(k_rddn_sfc,tbound) 
                                                                        
                  close (k_rddn_sfc) 
                                                                        
                  return 
                                                                        
                  endif 
                                                                        
!     *************  Surface derivatives completed  ****************    
                                                                        
!                                                                       
!     ---------------------                                             
!                                     SPECIAL CASE -> IMRG=40-43, IEMIT=
! analytic derivatives/jacobians                                        
!                                                                       
!     If IMRG = 41/43 and IEMIT = 3, then precalculated optical depths  
!     on multiple files, precalculated cumulative optical depths on     
!     multiple files, and just-calculated layer optical depths          
!     are combined to produce analytic layer radiance                   
!     derivatives from both space-to-ground and ground-to-space         
!     (written to PTHRDRd and PTHRDRu) as well as total upwelling       
!     radiance (written to TAPE12).  The results are monochromatic for  
!     IMRG = 41, scanned for IMRG = 43.                                 
!                                                                       
!     If IMRG = 40/42 and IEMIT = 3, then precalculated optical depths  
!     on multiple files, precalculated cumulative optical depths on     
!     multiple files, and just-calculated layer optical depths          
!     are combined to produce analytic layer radiance                   
!     derivatives from space-to-ground (written to PTHRDRd) as well     
!     as total downwelling radiance (written to TAPE12).                
!     The results are monochromatic for IMRG = 40, scanned for IMRG = 42
!                                                                       
                                                                        
!                                                                       
!---------------------------------------------------------------------  
                                                                        
!       kfile       10   h_kfile     TAPE10                  ksubl      
!       kradtot     18   h_radtot    RDDNlayer_00L           T_dn, R_dn 
!       kfilad      19               AJ/RDderivUPW_00_001    layer deriv
!       kodfil      17               ODint_001               optical dep
!       ktemp       88               AJ_mono                 mono anal. 
!       k_rddn_sfc  90               RDDNlayer_001           downwelling
                                                                        
!---------------------------------------------------------------------  
                                                                        
!****************************************************                   
! arrive at this point if doing upwelling Jacobians                     
                                                                        
                  IF (((IMRG.eq.41.or.IMRG.eq.43)) .AND. (IEMIT.EQ.3)   &
                  .and. icflg.ge.0) then                                
                                                                        
                     iemit = 1 
                                                                        
! iup_dn = 1 selects upwelling radiance                                 
                                                                        
                                  ! use for layer to level conversion (i
                     iup_dn = 1 
!                                                                       
!        Read card for scan for IMRG = 42,43                            
!                                                                       
                     IF (IMRG.GE.42) then 
                        CALL SCANRD_aj 
                        IMRGSAV=IMRG 
                        endif 
                                                                        
!     Get the pathname for the optical depth files                      
                        READ (IRD,945) PATH1 
                                                                        
                        CALL QNTIFY(PATH1,HFORM1) 
                                                                        
!                    nlayr, layr,                                       
                        CALL OPNODF( 1, 1,PATH1,HFORM1,-iemit) 
                                                                        
!     Get the maximum number of layers from the optical depth file for l
                                                                        
                        REWIND kodfil 
                        CALL BUFIN (kodfil,KEOF,FILDU1(1),NFHDRF) 
                                                                        
                        nlayer = nlayd1 
                        dv_lbl = dv_lbl1 
                                                                        
                        lh1 = nlayer 
                        lh2 = 1 
                                                                        
                        ipathl_sav = ipathl 
                        jpathl_sav = jpathl 
                                                                        
                        surf_refl_sav = surf_refl 
                        surf_refl = h_blank 
                                                                        
!     ipathl is set to 31 for algorithm in which the RT (downwelling) is
!            direction of the loop overlayers,  nlayer > 1.             
                                                                        
                        ipathl = 31 
                        jpathl = ipathl 
                                                                        
!     set up format for file name                                       
                        CALL QNTIFY(PATH1,HFORM1) 
                                                                        
!     set path name for rad down files                                  
                        PTHRAD = 'RDDNlayer_' 
                        CALL QNTIFY(PTHRAD,HFMRAD) 
                                                                        
!     now check for cross sections                                      
                                                                        
                        ixsect_sav = ixscnt/10 
                        ixsect = 0 
                                                                        
                        if (ixsect_sav .eq. 1 .and. nspcrt.eq.0) ixsect &
                        = 1                                             
                                                                        
                        IF (IXSECT.GE.1) THEN 
                                                                        
                           open (20,file='AJ_xs_amnts', form =          &
                           'unformatted',status='old')                  
                           read (20) IXMAX,IXMOLS, ( IXINDX(mol),(XAMNT(&
                           mol,l),l=1,nlayer),mol=1,ixmols )            
                           close (20) 
                           numxs = ixmols 
                           call xs_set(v1,v2) 
                                                                        
                        ENDIF 
!                                                                       
!****************     Now do the derivatives   ********************     
!                                                                       
!        Call OPPATH with layer = 0  which sets atmospheric path in subr
!                                                                       
                        layer = 0 
                        LAYHDR = LAYER 
                                                                        
                        iemit = 3 
                        ipathl = 1 
                        jpathl = 1 
!                                                                       
                        IF (2*(NLAYER/2).eq.NLAYER) then 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                           endif 
                                                                        
                           open(kfile,file='TAPE10',status='unknown',   &
                           form='unformatted')                          
                                                                        
                           lh1 = 1 
                           lh2 = nlayer 
                                                                        
                           LAYER = 0 
                                                                        
!***********   Loop over Layers   ************                          
                                                                        
    9                      LAYER = LAYER+1 
                           REWIND KFILE 
                                                                        
                           LAYHDR = LAYER 
!                                                                       
                           CALL OPNODF(nlayer,layer,PATH1,HFORM1,-iemit) 
                                                                        
!     Get the maximum number of layers from the optical depth file for l
                           REWIND kodfil 
                           CALL BUFIN (kodfil,KEOF,FILDU1(1),NFHDRF) 
                                                                        
!        d1 comes from FILHD1 and specifies the grid for the interpolate
!        dv_lbl1 comes from FILHD1 and specifies the grid for the lbl ca
                                                                        
                           dvout = d1 
                                                                        
                           dv_lbl= dv_lbl1 
                           dv = dv_lbl1 
                                                                        
                           V1 = VD1 
                           V2 = VD2 
                           SECANT = SECAN1 
                           secnt0 = secan1 
                                                                        
! ignore imrg=42/43 options for now (they are taken care of with IMRGSAV
                                                                        
                           imrg=41 
                                                                        
! make sure downwelling TAPE12 file is at beginning and                 
! then read the file header (data not used, so put into dummy common)   
                                                                        
                           rewind(k_rddn_sfc) 
                           call bufin(k_rddn_sfc,keof,filhdrj(1),nfhdrf) 
                                                                        
                           NLAYHD = NLAYER 
!                                                                       
                                                                        
!     SET UP LAYER BOUNDARY PARAMETERS                                  
!                                                                       
                           ALTZL = xALT1(1) 
                           ALTZU = xALT1(2) 
                           PZL = PDL 
                           PZU = PDU 
                           TZL = TDL 
                           TZU = TDU 
!                                                                       
                           PAVE = PAV1 
                           TAVE = TAV1 
                           WBROAD = wbrd1 
                                                                        
                           nmol = nmo1 
                                                                        
                           do m=1,nmol 
                           wk(m) = w1(m) 
                           enddo 
                                                                        
                           sample = 4. 
                                                                        
                           wtot = wbroad 
                           do m=1,nmol 
                           wtot = wtot + wk(m) 
                           enddo 
!                                                                       
                           IF (NSPCRT.EQ.0) THEN 
                                                                        
!    *******   set parameters for OD calculation for temperature  ****  
!                                                                       
!           set layer temperature for forward finite difference calculat
!           to be used for analytic derivative results (includes continu
!                                                                       
                           DELT_PERT = 1. 
                           TAVE = TAVE + DELT_PERT 
                                                                        
!           set cross section amounts for layer optical depth calculatio
                                                                        
                           IF (IXSECT.GE.1) THEN 
                           DO M = 1, IXMOLS 
                           WXM(M) = XAMNT(M,LAYER) 
                           ENDDO 
                           ENDIF 
                                                                        
                           ENDIF 
                                                                        
                           IF ((NSPCRT.GT.0).and.(NSPCRT.LE.MXMOL))     &
                           THEN                                         
                                                                        
!    ********   set parameters for OD calculation for selected molecule 
!                                                                       
!          set column amount for SPECIES analytic derivative            
!                                                                       
!          save column amount, zero all amounts, and reset old amount fo
!           molecule to be retrieved                                    
!                                                                       
!          w.r.t. amount (adjust dry air column for change in water vapo
!          (mass of layer must be held constant to maintain pressure lev
!                dUh2o=-1.0/(1.0+(1.609/wq))                            
!                                                                       
!**%%$$ ?????                                                           
                           frh2o =  wk(1)/wtot !mja, 10-27-2011
                           if (nspcrt.eq.1) then 
                           dUh2o = -frh2o/(frh2o+1.609) 
                           endif 
                                                                        
!           molecules (save one of interest)                            
!           also, wbroad is principally composed of o2 and n2;          
!           if o2 (mol=7) is set to zero, wbroad must include o2 to reta
                           wklsav = wk(nspcrt) 
                           wkl_7 = wk(7) 
                                                                        
                           DO M = 1,NMOL 
                           WK(M) = 0.0 
                           ENDDO 
                                                                        
                                                                        
! if change these two lines for abs coef or o.d., change write 1020/1021
!           WK(NSPCRT) = 1.0E20  ! use for absorption coef scaling      
                                                                        
                                                             ! use if o.
                           wk(nspcrt) = wklsav 
                           wbroad = wbrd1 + wkl_7 
                                                                        
                           ENDIF 
                                                                        
                           CALL OPDPTH (MPTS) 
!                                                                       
                           REWIND KFILE 
                           REWIND MFILE 
                           REWIND LFILE 
!                                                                       
!        Check to ensure derivative and radiance calculations           
!        are going in right direction.                                  
!                                                                       
! 40/42 = uplooking/downwelling  -> need ipathl=3                       
! 41/43 = downlooking/upwelling  -> need ipathl=1                       
                                                                        
                           IF (((IMRG.EQ.41).OR.(IMRG.EQ.43)).AND.(     &
                           IPATHL.EQ.31)) THEN                          
                           WRITE(IPR,*) 'XLAYER ERROR: IPATH PROBLEM' 
                           WRITE(IPR,940) IPATHL,IMRG 
                           WRITE(*,*) 'IPATHL,IMRG: ',IPATHL,IMRG 
                           STOP 'XLAYER ERROR: IPATH PROBLEM' 
                           ENDIF 
!                                                                       
!        Open files appropriate to derivative calculation:              
!                                                                       
!        KODFIL = optical depth file for layer (all molecules)          
!        kradtot = total radiance/transmittance file to layer (all molec
!        KTEMP  = outgoing monochromatic layer analytic derivatives     
!        KFILAD = outgoing scanned layer analytic derivatives           
!                                                                       
!        Other files appropriate to derivative calculation:             
!                                                                       
!        KFILE  = Absorptance coefficient file for layer & molecule (TAP
!        LFILE  = Incoming acculmulated radiance and transmittance      
!        MFILE  = Outgoing acculmulated radiance and transmittance      
!                                                                       
                           IF (IPATHL.EQ.1) THEN 
                           CALL OPNDRV(NLAYER,LAYER,LAYER+1,ipathl) 
                        ELSEIF (IPATHL.EQ.3) THEN 
                           CALL OPNDRV(1,NLAYER-LAYER+1,NLAYER-LAYER,   &
                           ipathl)                                      
                        ELSE 
                           STOP 'XLAYER: IPATHL NOT VALID' 
                           ENDIF 
!                                                                       
                                ! use to pass correct ipathl to emadl1 a
                           jpathl=ipathl 
                                                                        
                           CALL XMERGE (NPTS,LFILE,MFILE,JPATHL) 
                                                                        
                           REWIND MFILE 
                           MMFILE = KTEMP 
                           REWIND MMFILE 
                                          ! MOVE FROM KTEMP TO KFILAD   
                           CALL COPYFL (NPTS,MMFILE,KFILAD) 
                                ! PUTS -99 IN LAST LINE OF FILE         
                           CALL ENDFIL (KFILAD) 
                                                                        
                                                                        
! IMRGSAV=42,43 for scanning of jacobian files                          
!   If scanning, reset values of HWF1,DXF1,NX1,N1MAX which may          
!   have been been changed in HIRAC1 after having been read in          
!   in SCANRD, but before being used in SCNMRG.                         
!                                                                       
!   Scanning is done after all files are complete, since the            
!   downwelling term is needed in the calculation of the upwelling      
                                                                        
                           IF (IMRGSAV.GE.42) THEN 
                                                                        
                           HWF1 = HWFS 
                           DXF1 = DXFS 
                           NX1 = NFS 
                           N1MAX = NFMAXS 
                                                                        
                                ! MAKE SURE LAST FILE IS CLOSED         
                           CLOSE(KFILAD) 
                           CLOSE(KRADTOT) 
                                                                        
                           CALL SCANRD (DVINT,JEMIT,1) 
                           CALL SCNMRG_AJ(NLAYER,IUP_DN) 
                                                                        
                           ENDIF 
                                                                        
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                                                                        
                                       ! DO NEXT LAYER                  
                           IF (LAYER.LT. NLAYER) GO TO 9 
                                                                        
                                                                        
!*************  END LOOP OVER LAYERS   *************                    
                                                                        
                           CLOSE (KRADTOT) 
                           close (k_rddn_sfc) 
!                                                                       
! ------------------------------------------------                      
! IF DERIVATIVE CALCULATION (40 <= IMRG <= 43)                          
!  -AND - LAYER INPUT (IMOLDX > -99)                                    
! [ THESE TWO IMPLY THAT IMOLDX > -99 ]                                 
!  -AND - NOT A SURFACE TERM (NSPCRT < 0) THEN                          
! DO THE LAYER TO LEVEL CONVERSION OF THE DERIVATIVE FILES              
                                                                        
                           IF (NSPCRT.GE.0) THEN 
                                                                        
                           ilevdx = nlayer 
                                                                        
                           CALL LAYER2LEVEL 
                           ENDIF 
! ------------------------------------------------                      
                                                                        
                           return 
                                                                        
                        ENDIF 
                                                                        
!****************************************************                   
! arrive at this point if doing downwelling Jacobians                   
                                                                        
                        IF (((IMRG.eq.40.or.IMRG.eq.42)) .AND. (        &
                        IEMIT.EQ.3) .and. icflg.ge.0) then              
                                                                        
                           iemit = 1 
                                                                        
! iup_dn = -1 selects downwelling radiance                              
                                                                        
                                   ! use for layer to level conversion (
                           iup_dn = -1 
!                                                                       
!        Read card for scan for IMRG = 42,43                            
!                                                                       
                           IF (IMRG.GE.42) then 
                           CALL SCANRD_aj 
                           IMRGSAV=IMRG 
                           endif 
                                                                        
!     Get the pathname for the optical depth files                      
                           READ (IRD,945) PATH1 
                                                                        
                           CALL QNTIFY(PATH1,HFORM1) 
                                                                        
!                    nlayr, layr,                                       
                           CALL OPNODF( 1, 1,PATH1,HFORM1,-iemit) 
                                                                        
!     Get the maximum number of layers from the optical depth file for l
                                                                        
                           REWIND kodfil 
                           CALL BUFIN (kodfil,KEOF,FILDU1(1),NFHDRF) 
                                                                        
                           nlayer = nlayd1 
                           dv_lbl = dv_lbl1 
                                                                        
                           lh1 = nlayer 
                           lh2 = 1 
                                                                        
                           ipathl_sav = ipathl 
                           jpathl_sav = jpathl 
                                                                        
!     set up format for file name                                       
                           CALL QNTIFY(PATH1,HFORM1) 
                                                                        
!     set path name for rad down files                                  
                           PTHRAD = 'RDDNlayer_' 
                           CALL QNTIFY(PTHRAD,HFMRAD) 
                                                                        
!     now check for cross sections                                      
                                                                        
                           ixsect_sav = ixscnt/10 
                           ixsect = 0 
                                                                        
                           if (ixsect_sav .eq. 1 .and. nspcrt.eq.0)     &
                           ixsect = 1                                   
                                                                        
                           IF (IXSECT.GE.1) THEN 
                                                                        
                           open (20,file='AJ_xs_amnts', form =          &
                           'unformatted',status='old')                  
                           read (20) IXMAX,IXMOLS, ( IXINDX(mol),(XAMNT(&
                           mol,l),l=1,nlayer),mol=1,ixmols )            
                           close (20) 
                           numxs = ixmols 
                           call xs_set(v1,v2) 
                                                                        
                           ENDIF 
!                                                                       
!****************     Now do the derivatives   ********************     
!                                                                       
!        Call OPPATH with layer = 0  which sets atmospheric path in subr
!                                                                       
                           layer = 0 
                           LAYHDR = LAYER 
                                                                        
                           iemit = 3 
                           ipathl = 3 
                           jpathl = 3 
!                                                                       
                           IF (2*(NLAYER/2).eq.NLAYER) then 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                           endif 
                                                                        
                           open(kfile,file='TAPE10',status='unknown',   &
                           form='unformatted')                          
                                                                        
                           lh1 = 1 
                           lh2 = nlayer 
                                                                        
                           LAYER = 0 
                                                                        
!***********   Loop over Layers   ************                          
                                                                        
  719                      LAYER = LAYER+1 
                           REWIND KFILE 
                                                                        
                           LAYHDR = LAYER 
!                                                                       
                           CALL OPNODF(nlayer,layer,PATH1,HFORM1,-iemit) 
                                                                        
!     Get the maximum number of layers from the optical depth file for l
                           REWIND kodfil 
                           CALL BUFIN (kodfil,KEOF,FILDU1(1),NFHDRF) 
                                                                        
!        d1 comes from FILHD1 and specifies the grid for the interpolate
!        dv_lbl1 comes from FILHD1 and specifies the grid for the lbl ca
                                                                        
                           dvout = d1 
                                                                        
                           dv_lbl= dv_lbl1 
                           dv = dv_lbl1 
                                                                        
                           V1 = VD1 
                           V2 = VD2 
                           SECANT = SECAN1 
                           secnt0 = secan1 
                                                                        
! ignore imrg=42/43 options for now (they are taken care of with IMRGSAV
                                                                        
                           imrg=40 
                                                                        
! make sure downwelling TAPE12 file is at beginning and                 
! then read the file header (data not used, so put into dummy common)   
                                                                        
                           rewind(k_rddn_sfc) 
                           call bufin(k_rddn_sfc,keof,filhdrj(1),nfhdrf) 
                                                                        
                           NLAYHD = NLAYER 
!                                                                       
                                                                        
!     SET UP LAYER BOUNDARY PARAMETERS                                  
!                                                                       
                           ALTZL = xALT1(1) 
                           ALTZU = xALT1(2) 
                           PZL = PDL 
                           PZU = PDU 
                           TZL = TDL 
                           TZU = TDU 
!                                                                       
                           PAVE = PAV1 
                           TAVE = TAV1 
                           WBROAD = wbrd1 
                                                                        
                           nmol = nmo1 
                                                                        
                           do m=1,nmol 
                           wk(m) = w1(m) 
                           enddo 
                                                                        
                           sample = 4. 
                                                                        
                           wtot = wbroad 
                           do m=1,nmol 
                           wtot = wtot + wk(m) 
                           enddo 
!                                                                       
                           IF (NSPCRT.EQ.0) THEN 
                                                                        
!    *******   set parameters for OD calculation for temperature  ****  
!                                                                       
!           set layer temperature for forward finite difference calculat
!           to be used for analytic derivative results (includes continu
!                                                                       
                           DELT_PERT = 1. 
                           TAVE = TAVE + DELT_PERT 
                                                                        
!           set cross section amounts for layer optical depth calculatio
                                                                        
                           IF (IXSECT.GE.1) THEN 
                           DO M = 1, IXMOLS 
                           WXM(M) = XAMNT(M,LAYER) 
                           ENDDO 
                           ENDIF 
                                                                        
                           ENDIF 
                                                                        
                           IF ((NSPCRT.GT.0).and.(NSPCRT.LE.MXMOL))     &
                           THEN                                         
                                                                        
!    ********   set parameters for OD calculation for selected molecule 
!                                                                       
!          set column amount for SPECIES analytic derivative            
!                                                                       
!          save column amount, zero all amounts, and reset old amount fo
!           molecule to be retrieved                                    
!                                                                       
!          w.r.t. amount (adjust dry air column for change in water vapo
!          (mass of layer must be held constant to maintain pressure lev
!                dUh2o=-1.0/(1.0+(1.609/wq))                            
!                                                                       
!**%%$$ ?????                                                           
                           frh2o =  wk(1)/wtot !mja, 10-27-2011
                           if (nspcrt.eq.1) then 
                           dUh2o = -frh2o/(frh2o+1.609) 
                           endif 
                                                                        
!           molecules (save one of interest)                            
!           also, wbroad is principally composed of o2 and n2;          
!           if o2 (mol=7) is set to zero, wbroad must include o2 to reta
                           wklsav = wk(nspcrt) 
                           wkl_7 = wk(7) 
                                                                        
                           DO M = 1,NMOL 
                           WK(M) = 0.0 
                           ENDDO 
                                                                        
                                                                        
! if change these two lines for abs coef or o.d., change write 1020/1021
!           WK(NSPCRT) = 1.0E20  ! use for absorption coef scaling      
                                                                        
                                ! use if o.d. desired (dR/dlnx)         
                           wk(nspcrt) = wklsav 
                           wbroad = wbrd1 + wkl_7 
                                                                        
                           ENDIF 
                                                                        
                           CALL OPDPTH (MPTS) 
!                                                                       
                           REWIND KFILE 
                           REWIND MFILE 
                           REWIND LFILE 
!                                                                       
!        Check to ensure derivative and radiance calculations           
!        are going in right direction.                                  
!                                                                       
! 40/42 = uplooking/downwelling  -> need ipathl=3                       
                                                                        
                           IF (((IMRG.EQ.40).OR.(IMRG.EQ.42)).AND.(     &
                           IPATHL.NE.3)) then                           
                           WRITE(IPR,*) 'XLAYER ERROR: IPATH PROBLEM' 
                           WRITE(IPR,940) IPATHL,IMRG 
                           WRITE(*,*) 'IPATHL,IMRG: ',IPATHL,IMRG 
                           STOP 'XLAYER ERROR: IPATH PROBLEM' 
                           ENDIF 
!                                                                       
!        Open files appropriate to derivative calculation:              
!                                                                       
!        KODFIL = optical depth file for layer (all molecules)          
!        kradtot = total radiance/transmittance file to layer (all molec
!        KTEMP  = outgoing monochromatic layer analytic derivatives     
!        KFILAD = outgoing scanned layer analytic derivatives           
!                                                                       
!        Other files appropriate to derivative calculation:             
!                                                                       
!        KFILE  = Absorptance coefficient file for layer & molecule (TAP
!        LFILE  = Incoming acculmulated radiance and transmittance      
!        MFILE  = Outgoing acculmulated radiance and transmittance      
!                                                                       
                           CALL OPNDRV(nlayer,layer,layer+1,ipathl) 
!                                                                       
                                ! use to pass correct ipathl to emadl1 a
                           jpathl=ipathl 
                                                                        
                           CALL XMERGE (NPTS,LFILE,MFILE,JPATHL) 
                                                                        
                           REWIND MFILE 
                           MMFILE = KTEMP 
                           REWIND MMFILE 
                                          ! MOVE FROM KTEMP TO KFILAD   
                           CALL COPYFL (NPTS,MMFILE,KFILAD) 
                                ! PUTS -99 IN LAST LINE OF FILE         
                           CALL ENDFIL (KFILAD) 
                                                                        
                                                                        
! IMRGSAV=42,43 for scanning of jacobian files                          
!   If scanning, reset values of HWF1,DXF1,NX1,N1MAX which may          
!   have been been changed in HIRAC1 after having been read in          
!   in SCANRD, but before being used in SCNMRG.                         
!                                                                       
!   Scanning is done after all files are complete, since the            
!   downwelling term is needed in the calculation of the upwelling      
                                                                        
                           IF (IMRGSAV.GE.42) THEN 
                                                                        
                           HWF1 = HWFS 
                           DXF1 = DXFS 
                           NX1 = NFS 
                           N1MAX = NFMAXS 
                                                                        
                                ! MAKE SURE LAST FILE IS CLOSED         
                           CLOSE(KFILAD) 
                           CLOSE(KRADTOT) 
                                                                        
                           CALL SCANRD (DVINT,JEMIT,1) 
                           CALL SCNMRG_AJ(NLAYER,IUP_DN) 
                                                                        
                           ENDIF 
                                                                        
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                                                                        
                                         ! DO NEXT LAYER                
                           IF (LAYER.LT. NLAYER) GO TO 719 
                                                                        
!*************  END LOOP OVER LAYERS   *************                    
                                                                        
                           CLOSE (KRADTOT) 
                           close (k_rddn_sfc) 
!                                                                       
! ------------------------------------------------                      
! IF DERIVATIVE CALCULATION (40 <= IMRG <= 43)                          
!  -AND - LAYER INPUT (IMOLDX > -99)                                    
! [ THESE TWO IMPLY THAT IMOLDX > -99 ]                                 
!  -AND - NOT A SURFACE TERM (NSPCRT < 0) THEN                          
! DO THE LAYER TO LEVEL CONVERSION OF THE DERIVATIVE FILES              
                                                                        
                           IF (NSPCRT.GE.0) THEN 
                                                                        
                           ilevdx = nlayer 
                                                                        
                           CALL LAYER2LEVEL 
                           ENDIF 
! ------------------------------------------------                      
                                                                        
                           return 
                                                                        
                           ENDIF 
                                                                        
!****************************************************                   
                                                                        
                                                                        
!**************  END OF DERIVATIVE SECTION  *********                   
                                                                        
                                                                        
!*****************************************************                  
!                                                                       
!     ---------------------                                             
!                                                                       
!     If IMRG = 40 and IEMIT = 1, then precalculated optical depths     
!     on multiple files are combined to produce total downwelling radian
!     (from space to ground), written to TAPE12.  The results are       
!     monochromatic.                                                    
!                                                                       
!     If IMRG = 41 and IEMIT = 1, then precalculated optical depths     
!     on multiple files are combined to produce total upwelling radiance
!     (from ground to space), written to TAPE12.  The results are       
!     monochromatic.                                                    
!                                                                       
!                                     SPECIAL CASE -> IMRG=40/41, IEMIT=
!                                                                       
!**%%$$                                                                 
!       force program to skip this section as it conflicts with the     
!          current analytic derivative implementation                   
                                                                        
!      IF ((IMRG.EQ.40.OR.IMRG.EQ.41).AND.(IEMIT.EQ.1)) THEN            
                                                                        
                           iskip = 1 
                           IF ((IMRG.EQ.40.OR.IMRG.EQ.41).AND.(         &
                           iskip.EQ.0)) THEN                            
!                                                                       
!        -----------------------------                                  
!        Obtain information from KFILE                                  
!        -----------------------------                                  
!                                                                       
                           REWIND KFILE 
                           CALL BUFIN (KFILE,KEOF,FILDU1(1),NFHDRF) 
                           LTGNT = LTNSAV 
                           TBOUND = TMPBND 
                           EMISIV = BNDEMI(1) 
                           LH1 = LH1SAV 
                           LH2 = LH2SAV 
                                                                        
                           IPATHL = IPTHD1 
!                                                                       
!        Set number of layers upon which to perform radiative transfer  
!        to LAYTOT, read in from TAPE5 RECORD1.6a.                      
!                                                                       
                           NLAYER = LAYTOT 
!                                                                       
!        Test number of layer read in from TAPE5 RECORD1.6a to total    
!        number of layers information extracted from the fileheader     
!        from the first layer optical depth file.  If they do not       
!        agree, then issue a warning to TAPE6, and set LH2 to LAYTOT    
!        as well (for use with boundary temperature in XMERGE).         
!                                                                       
                           IF (LAYTOT.NE.NLAYD1) THEN 
                           WRITE(IPR,950) NLAYER,NLAYD1 
                           LH2 = LAYTOT 
                           ENDIF 
!                                                                       
!        Check for forced IPATHL, and set layer boundaries as needed    
!                                                                       
                           IF (IMRG.EQ.41) THEN 
                           IF (LH2.NE.1) THEN 
                           LH1 = MAX(LH1,LH2) 
                           LH2 = 1 
                           ENDIF 
                           IPATHL = 1 
                           JPATHL = 1 
                        ELSE 
                           IF (LH1.NE.1) THEN 
                           LH2 = MAX(LH1,LH2) 
                           LH1 = 1 
                           ENDIF 
                           IPATHL = 3 
                           JPATHL = 3 
                           ENDIF 
!                                                                       
!                           IF (LH1RD.GT.0) LH1=LH1RD !mja, these 
!                           variables aren't defined anywhere
!                           IF (LH2RD.GT.0) LH2=LH2RD 
                                                                        
                           IF (LH1.LT.LH2) IPATHL = 3 
                           JPATHL = IPATHL 
                                                                        
!                                                                       
!        Start of loop over layers                                      
!                                                                       
                           IF (2*(NLAYER/2).EQ.NLAYER) THEN 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                           ENDIF 
                                                                        
   10                      LAYER = 0 
   11                      LAYER = LAYER+1 
                                                                        
                           REWIND MFILE 
                           REWIND LFILE 
                           WRITE (IPR,900) 
                           WRITE (IPR,905) 
                           CALL OPNODF(NLAYER,LAYER,PATH1,HFORM1,IEMIT) 
                           WRITE (IPR,905) 
                           CALL XMERGE (NPTS,LFILE,MFILE,JPATHL) 
                           IF (LAYER.EQ.NLAYER) RETURN 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
!                                                                       
!        END OF LOOP OVER LAYERS                                        
!                                                                       
                           GO TO 11 
                           ENDIF 
!                                                                       
!     ---------------------                                             
!                                                                       
                           MMRG = MOD(IMRG,I_10) 
!                                                                       
!      TESTS FOR PRESTORE                                               
!                                                                       
                           IF (MMRG.EQ.2) GO TO 18 
                           IF (MMRG.EQ.5) GO TO 18 
                           IF (MMRG.EQ.6) GO TO 18 
                           IF (MMRG.EQ.8) GO TO 18 
                           IF (MMRG.EQ.9) GO TO 18 
!                                                                       
!     -------------------------------                                   
!     Call OPPATH, which calls PATH                                     
!     -------------------------------                                   
!                                                                       
                           LAYHDR = LAYER 
                           CALL OPPATH 
                           IF (IHIRAC.EQ.0) GO TO 170 
!                                                                       
!     NORMAL CALCULATIONS:                                              
                                                                        
!     INITIALIZE LH1, LH2 AND JPATHL                                    
                                                                        
                           if (ipathl .eq.3) then 
!        upwelling                                                      
                           LH1 = NLAYER 
                           LH2 = 1 
                           else 
!        downwelling                                                    
                           LH1 = 1 
                           LH2 = NLAYER 
                           endif 
                                                                        
                           JPATHL = IPATHL 
                                                                        
                           GO TO 20 
!                                                                       
!     -----------------------------                                     
!     Obtain information from KFILE                                     
!     -----------------------------                                     
!                                                                       
   18                      REWIND KFILE 
                           CALL BUFIN (KFILE,KEOF,FILDU1(1),NFHDRF) 
                           LTGNT = LTNSAV 
                           TBOUND = TMPBND 
                           EMISIV = BNDEMI(1) 
                           LH1 = LH1SAV 
                           LH2 = LH2SAV 
                           IPATHL = IPTHD1 
                           NLAYER = NLAYD1 
!                                                                       
!     For IMRG = 35,36,45,46, reset NLAYER to LAYTOT, read in from      
!     TAPE5 RECORD1.6a.  Test number of layer read in from TAPE5        
!     RECORD1.6a to total number of layers information extracted        
!     from the fileheader from the first layer optical depth file.      
!     If they do not agree, then issue a warning to TAPE6.              
!                                                                       
                           IF (IMRG.GE.35) THEN 
                           NLAYER = LAYTOT 
                           LH1 = NLAYER 
                           IF (LAYTOT.NE.NLAYD1) WRITE(IPR,950) NLAYER, &
                           NLAYD1                                       
                           ENDIF 
!                                                                       
!     -----------------------------                                     
!     Standard Merge Options Follow                                     
!     -----------------------------                                     
!                                                                       
!     CHECK FOR FORCED IPATHL, AND SET LAYER BOUNDARIES AS NEEDED       
!                                                                       
                           IF (IMRG.EQ.12.OR.MMRG.EQ.5) THEN 
                           IF (LH2.NE.1) THEN 
                           LH1 = MAX(LH1,LH2) 
                           LH2 = 1 
                           ENDIF 
                           ENDIF 
                           IF (IMRG.EQ.32.OR.MMRG.EQ.6) THEN 
                           IF (LH1.NE.1) THEN 
                           LH2 = MAX(LH1,LH2) 
                           LH1 = 1 
                           ENDIF 
                           ENDIF 
                           IF (IMRG.EQ.22.OR.MMRG.EQ.8) THEN 
                           IF (IPATHL.NE.2) THEN 
                           LH1 = MAX(LH1,LH2) 
                           LH2 = LH1 
                           ENDIF 
                           ENDIF 
!                                                                       
!       PRESTORE IMRG TEST, AND TEST FOR SCANNED OR FILTERED            
!                                        WEIGHTING FUNCTIONS            
!                                                                       
   20                      IF (MMRG.EQ.2) GO TO 50 
                           IF (MMRG.EQ.5) GO TO 90 
                           IF (MMRG.EQ.6) GO TO 50 
                           IF (MMRG.EQ.8) GO TO 110 
                           IF (MMRG.EQ.9) GO TO 150 
                           IF (IMRG.GT.10) GO TO 50 
                           IF (MMRG.GE.1) GO TO 40 
!                                                                       
!                                                                       
!    START OF LOOP OVER LAYERS                    IMRG=0                
!                                                                       
                           IF (2*(NLAYER/2).NE.NLAYER) GO TO 30 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
   30                      LAYER = LAYER+1 
                           REWIND KFILE 
                           LAYHDR = LAYER 
                           CALL OPPATH 
                           NLAYHD = NLAYER 
                           CALL OPDPTH (MPTS) 
                           REWIND KFILE 
                           REWIND MFILE 
                           REWIND LFILE 
                           CALL XMERGE (NPTS,LFILE,MFILE,JPATHL) 
                           IF (LAYER.EQ.NLAYER) GO TO 170 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
!                                                                       
!    END OF LOOP OVER LAYERS                                            
!                                                                       
                           GO TO 30 
!                                                                       
!     START LOOP OVER LAYERS                      IMRG=3,4,7            
!     OPTICAL DEPTH STRUNG OUT ON TAPE10                                
!                                                                       
   40                      LAYER = LAYER+1 
                           LAYHDR = LAYER 
                           CALL OPPATH 
                           NLAYHD = NLAYER 
                           CALL OPDPTH (MPTS) 
                           CALL ENDFIL (KFILE) 
                           REWIND MFILE 
                           REWIND LFILE 
                           IF (LAYER.EQ.NLAYER) GO TO 50 
                           GO TO 40 
!                                                                       
!    END OF LOOP OVER LAYERS                                            
!                                                                       
!                                 IMRG=2,4,6,12,14,16,22,24,26,32,36,46 
!                                                                       
   50                      IF (MMRG.EQ.1) GO TO 170 
                           IF (MMRG.EQ.3) GO TO 90 
                           IF (MMRG.EQ.7) GO TO 110 
!                                                                       
!       READ CARD FOR SCAN OF WEIGHTING FUNCTION                        
!       READ CARD FOR FILTER OF WEIGHTING FUNCTION                      
!                                                                       
                           IF (IMRG.EQ.14.OR.IMRG.EQ.16.OR.IMRG.EQ.36)  &
                           CALL SCANRD (DVINT,JEMIT,0)                  
                           IF (IMRG.EQ.24.OR.IMRG.EQ.26) CALL FLTRRD 
!                                                                       
!     START LOOP OVER LAYERS                                            
!     ALL MERGING DONE HERE                                             
!                                                                       
                           REWIND KFILE 
                           REWIND LFILE 
                           REWIND MFILE 
!                                                                       
!    DETERMINE IF FORCED IPATHL, AND SET APPROPRIATELY                  
!     This occurs when you are reading in precalculated KFILE           
!     with imrg options for merged (not sequential) output              
!     imrg=12 upwelling, imrg=22 tangent, imrg=32 downwelling           
!                                                                       
                           IF (IMRG.NE.2) THEN 
                           IF (IMRG.EQ.12) THEN 
                           JPATHL = 1 
                           else if (imrg.eq.22) then 
                           JPATHL = 2 
                        ELSE 
                           JPATHL = 3 
                           ENDIF 
                           ENDIF 
!                                                                       
                           NNTAN = 1 
                           NTAN(NNTAN) = 1 
                           IF (2*(NLAYER/2).NE.NLAYER) GO TO 60 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
   60                      LAYER = 0 
   70                      LAYER = LAYER+1 
                           WRITE (IPR,900) 
                           WRITE (IPR,905) 
                           IF (MMRG.EQ.2) WRITE (IPR,910) NLAYER 
                           IF (MMRG.EQ.4.OR.MMRG.EQ.6) THEN 
                           WRITE (IPR,915) NLAYER 
                           WRITE (IPR,920) (NTAN(N),N=1,NNTAN) 
                           ENDIF 
                           WRITE (IPR,905) 
                           IF (IMRG.EQ.14.OR.IMRG.EQ.24) THEN 
                           REWIND KFILE 
                           LAYHDR = LAYER 
                           CALL OPPATH 
                           NLAYHD = NLAYER 
                           CALL OPDPTH (MPTS) 
                           REWIND KFILE 
                           ENDIF 
!                                                                       
!     Open layer optical depth file for IMRG=36,46.  If IMRG=46, then   
!     open output layer radiance file.                                  
!                                                                       
                           IF ((IMRG.EQ.36).OR.(IMRG.EQ.46)) THEN 
                           CALL OPNODF(NLAYER,LAYER,PATH1,HFORM1,IEMIT) 
                           IF (IMRG.EQ.46) THEN 
                           PTHRAD = 'RDUPlayer_' 
                           CALL QNTIFY(PTHRAD,HFMRAD) 
                           CALL OPNRAD(nfile,NLAYER,LAYER) 
                           ENDIF 
                           ENDIF 
                           CALL XMERGE (NPTS,LFILE,MFILE,JPATHL) 
                           NNTAN = NNTAN+1 
                           NTAN(NNTAN) = NTAN(NNTAN-1)+1 
                           IF (MMRG.EQ.2) GO TO 80 
                           REWIND MFILE 
                           IF (IMRG.EQ.4.OR.IMRG.EQ.6.OR.IMRG.EQ.46)    &
                           CALL COPYFL (NPTS,MFILE,NFILE)               
!                                                                       
!   FOR SCAN CASE, IF DV NOT FINE ENOUGH, FIRST INTERPOLATE             
!                                                                       
                           IF (IMRG.EQ.14.OR.IMRG.EQ.16.OR.IMRG.EQ.36)  &
                           THEN                                         
                           MMFILE = MFILE 
                           IF (DVXM.GT.DVINT) THEN 
                           CALL SCNINT (MFILE,LFILE,DVINT,JEMIT,NPTS,   &
                           IBUF)                                        
                           MMFILE = LFILE 
                           ENDIF 

!****************************************************************** 
!        MJA, 10-18-2012: Added this to the SCNMRG call for 
!        IMRG = 14, 16, or 36 to match the call below 
!        for IMRG = 13, 15, or 35   
!                                       
!        If scanning, reset values of HWF1,DXF1,NX1,N1MAX which may     
!        have been been changed in HIRAC1 after having been read in     
!        in SCANRD, but before being used in SCNMRG.                    
!                                                                       
                           HWF1 = HWFS 
                           DXF1 = DXFS 
                           NX1 = NFS 
                           N1MAX = NFMAXS
                           CALL SCNMRG (MMFILE,NFILE) 
!******************************************************************
                           ENDIF 
                           IF (IMRG.EQ.24.OR.IMRG.EQ.26) CALL FLTMRG (  &
                           MFILE,NFILE)                                 
                           CALL ENDFIL (NFILE) 
   80                      IF (LAYER.EQ.NLAYER) GO TO 170 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                           REWIND MFILE 
                           REWIND LFILE 
                           GO TO 70 
!                                                                       
!                                            IMRG=3,5,13,15,23,25,35,45 
!                                                                       
!       MODIFIED TO BEGIN WEIGHTING FUNCTION CALC. FROM H1              
!                                                                       
!                                                                       
!     START LOOP OVER LAYERS                                            
!     ALL MERGING DONE HERE                                             
!                                                                       
   90                      REWIND KFILE 
                           REWIND LFILE 
                           REWIND MFILE 
!                                                                       
!       READ CARD FOR SCAN OF WEIGHTING FUNCTION                        
!       READ CARD FOR FILTER OF WEIGHTING FUNCTION                      
!                                                                       
                           IF (IMRG.EQ.13.OR.IMRG.EQ.15.OR.IMRG.EQ.35)  &
                           CALL SCANRD (DVINT,JEMIT,0)                  
                           IF (IMRG.EQ.23.OR.IMRG.EQ.25) CALL FLTRRD 
                           NNTAN = 1 
                           NTAN(NNTAN) = NLAYER 
                           IF (2*(NLAYER/2).EQ.NLAYER) THEN 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                           ENDIF 
!                                                                       
!    SET VALUE FOR IPATHL                                               
!                                                                       
                           JPATHL = 1 
                           LAYER = NLAYER+1 
  100                      LAYER = LAYER-1 
                           ISKIP = LAYER-1 
                           WRITE (IPR,900) 
                           WRITE (IPR,905) 
                           WRITE (IPR,925) NLAYER,LAYER 
                           WRITE (IPR,920) (NTAN(N),N=1,NNTAN) 
                           WRITE (IPR,905) 
!                                                                       
!     Open layer optical depth file for IMRG=35,45.  If IMRG=45, then   
!     open output layer radiance file.                                  
!                                                                       
                           IF (IMRG.EQ.13.OR.IMRG.EQ.23) THEN 
                           REWIND KFILE 
                           LAYHDR = LAYER 
                           CALL OPPATH 
                           NLAYHD = NLAYER 
                           CALL OPDPTH (MPTS) 
                           REWIND KFILE 
                        ELSEIF (IMRG.EQ.35.OR.IMRG.EQ.45) THEN 
                           CALL OPNODF(1,LAYER,PATH1,HFORM1,IEMIT) 
                           IF (IMRG.EQ.45) THEN 
                           PTHRAD = 'RDDNlayer_' 
                           CALL QNTIFY(PTHRAD,HFMRAD) 
                           CALL OPNRAD(nfile,NLAYER,LAYER) 
                           ENDIF 
                        ELSE 
                           CALL SKIPFL (ISKIP,KFILE,IEOF) 
                           ENDIF 
                           CALL XMERGI (NPTS,LFILE,MFILE,JPATHL) 
                           NNTAN = NNTAN+1 
                           NTAN(NNTAN) = NTAN(NNTAN-1)-1 
                           REWIND MFILE 
                           IF (IMRG.EQ.3.OR.IMRG.EQ.5.OR.IMRG.EQ.45)    &
                           CALL COPYFL (NPTS,MFILE,NFILE)               
!                                                                       
!     FOR SCAN CASE, IF DV NOT FINE ENOUGH, FIRST INTERPOLATE           
!                                                                       
                           IF (IMRG.EQ.13.OR.IMRG.EQ.15.OR.IMRG.EQ.35)  &
                           THEN                                         
                           MMFILE = MFILE 
                           IF (DVXM.GT.DVINT) THEN 
                           CALL SCNINT (MFILE,LFILE,DVINT,JEMIT,NPTS,   &
                           IBUF)                                        
                           MMFILE = LFILE 
                           ENDIF 
!                                                                       
!        If scanning, reset values of HWF1,DXF1,NX1,N1MAX which may     
!        have been been changed in HIRAC1 after having been read in     
!        in SCANRD, but before being used in SCNMRG.                    
!                                                                       
                           HWF1 = HWFS 
                           DXF1 = DXFS 
                           NX1 = NFS 
                           N1MAX = NFMAXS 
                           CALL SCNMRG (MMFILE,NFILE) 
                           ENDIF 
                           IF (IMRG.EQ.23.OR.IMRG.EQ.25) CALL FLTMRG (  &
                           MFILE,NFILE)                                 
                           CALL ENDFIL (NFILE) 
                           IF (LAYER.EQ.1) GO TO 170 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                           REWIND KFILE 
                           REWIND LFILE 
                           REWIND MFILE 
                           GO TO 100 
!                                                                       
!   ** TANGENT **                                 IMRG=7,8,17,18,27,28  
!                                                                       
  110                      REWIND KFILE 
                           REWIND LFILE 
                           REWIND MFILE 
!                                                                       
!     READ CARD FOR SCAN OF WEIGHTING FUNCTION                          
!     READ CARD FOR FILTER OF WEIGHTING FUNCTION                        
!                                                                       
                           IF (IMRG.EQ.17.OR.IMRG.EQ.18) CALL SCANRD (  &
                           DVINT,JEMIT,0)                               
                           IF (IMRG.EQ.27.OR.IMRG.EQ.28) CALL FLTRRD 
!                                                                       
!     **                 DOWN                                           
!                                                                       
                           IANT = 1 
                           NNTAN = 1 
                           NTAN(NNTAN) = LH1 
                           NLTOTL = LH1+LH2 
                           IF (2*(NLTOTL/2).EQ.NLTOTL) THEN 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                           ENDIF 
!                                                                       
!    SET VALUE FOR IPATHL                                               
!                                                                       
                           JPATHL = 3 
                           LAYER = LH1+1 
  120                      LAYER = LAYER-1 
                           ISKIP = LAYER-1 
                           WRITE (IPR,900) 
                           WRITE (IPR,905) 
                           WRITE (IPR,930) NLAYER,LAYER 
                           WRITE (IPR,920) (NTAN(N),N=1,NNTAN) 
                           WRITE (IPR,905) 
                           IF (IMRG.EQ.17.OR.IMRG.EQ.27) THEN 
                           REWIND KFILE 
                           LAYHDR = LAYER 
                           CALL OPPATH 
                           NLAYHD = NLAYER 
                           CALL OPDPTH (MPTS) 
                           REWIND KFILE 
                        ELSE 
                           CALL SKIPFL (ISKIP,KFILE,IEOF) 
                           ENDIF 
                           CALL XMERGI (NPTS,LFILE,MFILE,JPATHL) 
                           NNTAN = NNTAN+1 
                           NTAN(NNTAN) = NTAN(NNTAN-1)-1 
                           REWIND MFILE 
                           IF (IMRG.EQ.7.OR.IMRG.EQ.8) CALL COPYFL (    &
                           NPTS,MFILE,NFILE)                            
!                                                                       
!     FOR SCAN CASE, IF DV NOT FINE ENOUGH, FIRST INTERPOLATE           
!                                                                       
                           IF (IMRG.EQ.17.OR.IMRG.EQ.18) THEN 
                           MMFILE = MFILE 
                           IF (DVXM.GT.DVINT) THEN 
                           CALL SCNINT (MFILE,LFILE,DVINT,JEMIT,NPTS,   &
                           IBUF)                                        
                           MMFILE = LFILE 
                           ENDIF 
                           CALL SCNMRG (MMFILE,NFILE) 
                           ENDIF 
                           IF (IMRG.EQ.27.OR.IMRG.EQ.28) CALL FLTMRG (  &
                           MFILE,NFILE)                                 
                           CALL ENDFIL (NFILE) 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                           REWIND KFILE 
                           REWIND LFILE 
                           REWIND MFILE 
                           IF (LAYER.EQ.1) GO TO 130 
                           GO TO 120 
!                                                                       
!     **                 UP                                             
!                                                                       
  130                      LAYER = 0 
                           IANT = -1 
                           NTAN(NNTAN) = 1 
  140                      LAYER = LAYER+1 
                           WRITE (IPR,900) 
                           WRITE (IPR,905) 
                           WRITE (IPR,930) NLAYER,LAYER 
                           WRITE (IPR,920) (NTAN(N),N=1,NNTAN) 
                           WRITE (IPR,905) 
                           IF ((IMRG.EQ.17.OR.IMRG.EQ.27)               &
                           .AND.LAYER.NE.1) THEN                        
                           REWIND KFILE 
                           LAYHDR = LAYER 
                           CALL OPPATH 
                           NLAYHD = NLAYER 
                           CALL OPDPTH (MPTS) 
                           REWIND KFILE 
                           ENDIF 
                           IF (LAYER.LE.LH1) CALL XMERGI (NPTS,LFILE,   &
                           MFILE,JPATHL)                                
                           IF (LAYER.GT.LH1) CALL XMERGE (NPTS,LFILE,   &
                           MFILE,JPATHL)                                
                           NNTAN = NNTAN+1 
                           NTAN(NNTAN) = NTAN(NNTAN-1)+1 
                           REWIND MFILE 
                           IF (IMRG.EQ.7.OR.IMRG.EQ.8) CALL COPYFL (    &
                           NPTS,MFILE,NFILE)                            
!                                                                       
!     FOR SCAN CASE, IF DV NOT FINE ENOUGH, FIRST INTERPOLATE           
!                                                                       
                           IF (IMRG.EQ.17.OR.IMRG.EQ.18) THEN 
                           MMFILE = MFILE 
                           IF (DVXM.GT.DVINT) THEN 
                           CALL SCNINT (MFILE,LFILE,DVINT,JEMIT,NPTS,   &
                           IBUF)                                        
                           MMFILE = LFILE 
                           ENDIF 
                           CALL SCNMRG (MMFILE,NFILE) 
                           ENDIF 
                           IF (IMRG.EQ.27.OR.IMRG.EQ.28) CALL FLTMRG (  &
                           MFILE,NFILE)                                 
                           CALL ENDFIL (NFILE) 
                           IF (LAYER.EQ.LH2) GO TO 170 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                           REWIND LFILE 
                           REWIND MFILE 
                           GO TO 140 
!                                                                       
!     START LOOP OVER LAYERS                      IMRG=9                
!     OPTICAL DEPTH STRUNG OUT ON TAPE10                                
!                                                                       
!                                                                       
!     START LOOP OVER LAYERS                                            
!     ALL MERGING DONE HERE                                             
!                                                                       
  150                      REWIND MFILE 
                           IF (IAERSL.GE.1 .and. iaersl.ne.5 ) THEN 
                           REWIND 20 
                           IREAD = 0 
                           LOWFLG = 3 
                           LAYHDR = LAYER 
                           CALL OPPATH 
                           ENDIF 
                           REWIND LFILE 
                           REWIND KFILE 
                           IF (2*(NLAYER/2).NE.NLAYER) GO TO 160 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                           LAYER = 0 
  160                      LAYER = LAYER+1 
                           WRITE (IPR,900) 
                           CALL XMERGE (NPTS,LFILE,MFILE,JPATHL) 
                           IF (LAYER.EQ.NLAYER) GO TO 170 
                           MSTOR = MFILE 
                           MFILE = LFILE 
                           LFILE = MSTOR 
                           REWIND MFILE 
                           REWIND LFILE 
                           GO TO 160 
!                                                                       
  170                      IF (IMRG.GE.23.AND.IMRG.LE.28) CALL FLTPRT ( &
                           NFILE)                                       
!                                                                       
                           RETURN 
!                                                                       
  900 FORMAT ('1') 
  905 FORMAT (/,'*******************************************',/) 
  910 FORMAT ('0  SEQUENTIAL OPTICAL DEPTHS TO KFILE, MERGED TO MFILE', &
     &        I3,' LAYERS')                                             
  915 FORMAT ('0  GROUND TO SPACE WEIGHTING FUNCTION, LAYER 1 TO LAYER' &
     &        ,I3)                                                      
  920 FORMAT (/,'  THE WEIGHTING FUNCTION ACCUMULATION BY LAYER',/,     &
     &        (5X,20I3))                                                
  925 FORMAT (/,'  SPACE TO GROUND WEIGHTING FUNCTION, LAYER',I3,       &
     &        ' TO LAYER',I3)                                           
  930 FORMAT (/,'  TANGENT WEIGHTING FUNCTION, LAYER',I3,' TO LAYER',   &
     &        I3)                                                       
  935 FORMAT (/,' ---------------------------------------------- ',/,   &
     &          '  Results of Optical Depth Merging - Downlooking',/,   &
     &          ' ----------------------------------------------',//)   
  936 FORMAT (/,' ---------------------------------------------- ',/,   &
     &          '  Results of Optical Depth Merging - Uplooking',/,     &
     &          ' ----------------------------------------------',//)   
  940 FORMAT ('TAPE5: IPATHL',I5,', IMRG = ',I5) 
  945 FORMAT (A55,1X,I4,2I5) 
  946 FORMAT (A55) 
  950 FORMAT ('***********************************************',/,      &
     &        '*                   WARNING                   *',/,      &
     &        '*                  =========                  *',/,      &
     &        '*   TOTAL NUMER OF LAYERS LAYTOT FROM TAPE5   *',/,      &
     &        '*       NOT EQUAL TO LAYER TOTAL NLAYD1       *',/,      &
     &        '*   EXTRACTED FROM OPTICAL DEPTH FILEHEADER   *',/,      &
     &        '*                                             *',/,      &
     &        '***********************************************',//,     &
     &        'LAYTOT = ',I4,'NLAYD1 = ',I4)                            
 1015 FORMAT (3x,I2) 
!                                                                       
end subroutine XLAYER
!                                                                       
!     -------------------------------------------------------------     
!                                                                       
      SUBROUTINE OPPATH 
!                                                                       
      USE phys_consts, ONLY: pi
      USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,             &
     &                     MXISOTPL, ISOTPL_ABD, ISOTPL_NUM,             &
     &                     MXMOL, MX_XS, MXTRAC, IPTS
      IMPLICIT REAL*8           (V) 
!                                                                       
!     OPPATH CALLS LBLATM AND CALLS PATH FIRST                          
!                                                                       
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &     MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,mx_xs=38,MXTRAC=22) 
!                                                                       
      COMMON /PATHD/ PAVEL(MXLAY),TAVEL(MXLAY),WKL(MXMOL,MXLAY),        &
     &               WBRODL(MXLAY),DVL(MXLAY),                          &
     &               WTOTL(MXLAY),ALBL(MXLAY),ADBL(MXLAY),              &
     &               AVBL(MXLAY),H2OSL(MXLAY),                          &
     &               IPTH(MXLAY),ITYL(MXLAY),SECNTA(MXLAY),             &
     &               HT1,HT2,ALTZ(0:MXLAY),                             &
     &               PZ(0:MXLAY),TZ(0:MXLAY)                            
!                                                                       
!     IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE   
!     MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES     
!     (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR   
!     L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD       
!     MOLECULES.                                                        
!                                                                       
      COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(mx_xs),XAMNT(mx_xs,MXLAY) 
!                                                                       
!     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES       
!     FOR THE CROSS-SECTION MOLECULES.                                  
!                                                                       
      COMMON /XSECTR/ V1FX(5,MX_XS),V2FX(5,MX_XS),DVFX(5,MX_XS),        &
     &                WXM(MX_XS),NTEMPF(5,MX_XS),NSPECR(MX_XS),         &
     &                IXFORM(5,MX_XS),XSMASS(MX_XS),XDOPLR(5,MX_XS),    &
     &                NUMXS,IXSBIN                                      
!                                                                       
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
     &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
     &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
     &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
     &              EXTID(10)                                           
      CHARACTER*8  EXTID 
!                                                                       
      COMMON /PATH_ISOTPL/ ISOTPL,NISOTPL,                              &
     &                     ISOTPL_FLAG(MXMOL,MXISOTPL),                 &
     &                     ISOTPL_MAIN_FLAG(MXMOL),                     &
     &                     MOLNUM(MXMOL*MXISOTPL),                      &
     &                     ISOTPLNUM(MXMOL*MXISOTPL),                   &
     &                     WKI(MXMOL,MXISOTPL)             
      COMMON /OPPATH_ISOTPL/ WKL_ISOTPL(MXMOL,MXISOTPL,MXLAY)
!                                                                       
      common /lbl_geo/ zh1,zh2,zangle 
!                                                                       
      character*8      XID,       HMOLID,      YID 
      real*8               SECANT,       XALTZ 
!                                                                       
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
     &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
     &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF  
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,pad_3,&
     &                angle_path,secant_diffuse,secant_path,diffuse_fac 
!                                                                       
      character*1 surf_refl 
      character*3 pad_3 
!                                                                       
      CHARACTER*8       HLINID,BMOLID,HID1 
!                                                                       
      integer *4 molcnt,mcntlc,                                         &
     &           mcntnl,linmol,                                         &
     &           lincnt,ilinlc,ilinnl,irec,irectl                       
!                                                                       
      real *4 sumstr,flinlo,flinhi 
!                                                                       
      COMMON /LINHDR/ HLINID(10),BMOLID(64),MOLCNT(64),MCNTLC(64),      &
     &                MCNTNL(64),SUMSTR(64),LINMOL,FLINLO,FLINHI,       &
     &                LINCNT,ILINLC,ILINNL,IREC,IRECTL,HID1(2),LSTWDL   
!                                                                       
!     LSTWD (LAST WORD) IS DUMMY, DOES NOT NEED TO BE COUNTED           
!                                                                       
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
     &              NLTEFL,LNFIL4,LNGTH4                                
      COMMON /MSCONS/ AIRMAS(MXLAY),TGRND,SEMIS(3),HMINMS,HMAXMS,       &
     &                MSFLAG,                                           &
     &                MSWIT,IODFIL,MSTGLE                               
      COMMON /LASIV/ VLAS,ILAS 
      COMMON /ADRIVE/ LOWFLG,IREAD,MODEL,ITYPE,NOZERO,NP,H1F,H2F,       &
     &                ANGLEF,RANGEF,BETAF,LENF,AV1,AV2,RO,IPUNCH,       &
     &                XVBAR, HMINF,PHIF,IERRF,HSPACE                    
      COMMON /CNTRL/ I1,I2,I3,I4,NBNDL,I6,I7,NBNDF,I9 
!                                                                       
!     Common blocks for analytic derivative                             
!     -------------------------                                         
      COMMON /IADFLG/ NSPCRT,imrgsav 
!                                                                       
!      parameter (ipts=5050) 
      common /CDERIV/ icflg,idum,v1absc,v2absc,dvabsc,nptabsc,delT_pert,&
     &    dqh2oC(ipts),dTh2oC(ipts),dUh2o                               
!     -------------------------                                         
!                                                                       
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),              &
     &            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),              &
     &            (FSCDID(5),IEMIT) , (FSCDID(6),ISCNHD),               &
     &            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL),               &
     &            (FSCDID(9),JRAD) , (FSCDID(10),ITEST),                &
     &            (FSCDID(11),IMRG) , (FSCDID(12),SCNID),               &
     &            (FSCDID(13),HWHM) , (FSCDID(14),IDABS),               &
     &            (FSCDID(15),IATM) , (FSCDID(16),LAYR1),               &
     &            (FSCDID(17),NLAYHD)                                   
!                                                                       
!     NEW EQUIVALENCE STATEMENTS FOR TANGENT WEIGHTING FNS.             
!                                                                       
      EQUIVALENCE (YID(10),LTNSAV) , (YID(9),LH1SAV) , (YID(8),LH2SAV) 
      EQUIVALENCE (XALTZ(1),ALTZL) , (XALTZ(2),ALTZU) 
!                                                                       
      character*4 ht1,ht2 
!                                                                       
      DATA XDV / 64. / 
      DATA HZ / 6.5 / 
!                                                                       
      NBNDF = 0 
!                                                                       
!     SET AV1 AND AV2 FOR PASS TO LBLATM AND LOWTRN                     
!                                                                       
      AV1 = V1 
      AV2 = V2 
!                                                                       
!     RESET ALFAL0 TO DEFAULT OF 0.04 IN VALUE READ IN .LE. 0.0         
!     FOR LINE COUPLING USER SHOULD READ IN VALUE OF 0.05               
!                                                                       
      IF (SAMPLE.LT.0.01) SAMPLE = 4. 
      IF (ALFAL0.LE.0.) ALFAL0 = 0.04 
      IF (AVMASS.LE.0.) AVMASS = 36. 
!                                                                       
!     Write SAMPLE, ALFAL0, and AVMASS if IHIRAC = 0 and IATM >= 1      
!     (printed for informational puroses in case IFXTYP = 1).           
!                                                                       
      IF ((IHIRAC.EQ.0).AND.(IATM.GE.1))                                &
     &     WRITE (IPR,890) SAMPLE,ALFAL0,AVMASS                         
!                                                                       
      IF (LAYER.LE.0) THEN 
         IF (IAERSL.EQ.0 .or. iaersl.eq.5 .OR.IAERSL.EQ.9) THEN 
!                                                                       
!           IF LAYER GT 0 SAVE THE VECTORS                              
!                                                                       
            IF (IATM.GE.1.AND.IATM.LE.5) THEN 
!                                                                       
!              CALL LBLATM TO COMPLETE GEOMETRY ON FINAL COMBINED LAYERS
!              IF AEROSOLS PRESENT AND HORIZONTAL PATH SKIP CALL        
!              TO LBLATM                                                
!                                                                       
               CALL LBLATM 
               IF (IHIRAC.EQ.0) RETURN 
            ENDIF 
         ELSE 
!                                                                       
            IF (IATM.EQ.0) THEN 
               WRITE (IPR,900) 
               STOP 'OPPATH; IAER GT 1, IATM = 0 ' 
            ENDIF 
!                                                                       
            IREAD = 0 
!                                                                       
            CALL LBLATM 
!                                                                       
            NBNDF = NLAYRS+1 
!                                                                       
!     STORE V1 AND V2 IN V1S AND V2S DURING CALL TO LOWTRN              
!                                                                       
            V1S = V1 
            V2S = V2 
!                                                                       
            CALL LOWTRN 
!                                                                       
!     RESET V1 AND V2 TO V1S AND V2S                                    
!                                                                       
            V1 = V1S 
            V2 = V2S 
!                                                                       
         ENDIF 
!                                                                       
         IF (IATM.LE.5) THEN 
!                                                                       
            WRITE (IPR,905) 
!                                                                       
            WRITE (IPR,910) V1,V2,SAMPLE,DVSET,ALFAL0,AVMASS,DPTMIN,    &
            DPTFAC                                                      
!                                                                       
!     BEGINNING AND ENDING WAVENUMBER VALUES                            
!                                                                       
            IF ((V2-V1).GT.2020.) THEN 
               WRITE (IPR,915) 
               STOP 'OPPATH; V2-V1 GT 2020' 
            ENDIF 
!                                                                       
            IF ((IHIRAC.EQ.1.OR.IHIRAC.EQ.4.OR.ILBLF4.GE.1).AND.        &
            (IATM.EQ.1)) THEN                                           
               IF (NMOL.LT.LINMOL) WRITE (IPR,920) (LINMOL-NMOL) 
               IF (NMOL.GT.LINMOL) WRITE (IPR,925) NMOL,LINMOL 
            ENDIF 
!                                                                       
            CALL PATH 
!                                                                       
!     *************************************************************     
!     Compute the diffuse_fac for lambertian surface reflection         
!                                                                       
!     angle_path is the effective angle for the calculation from H2 to H
!            For IATM.eq.0, this angle is read in on card 2.1 as zangle 
!            For IATM.ne.0, this angle is obtained from the lblatm ray t
!                                                                       
!     For a lambertian surface, the flux is obtained from a radiance cal
!     at the diffusivity angle (secant=1.67) using the information avail
!     from the calculation from H2 to H1 (at the effective angle, angle_
!                                                                       
!     diffuse_fac is the ratio of this secant value (1.67) to the secant
!     associated with angle_path.  This factor is used to obtain the opt
!     (subroutine EMIN) and the total ray transmittances (see module XME
!     for the downwelling contribution at the diffusivity angle.        
!                                                                       
            if (surf_refl .ne. 's') then 
            IF (IATM.EQ.0) THEN 
               angle_path = zangle 
                                                                        
               if ( ((angle_path-180.) .gt. 0.001) .or. ((angle_path-   &
               90.) .lt. -0.001) ) then                                 
                        write(*,*) '     For lambertian, surf_refl = l' 
               stop 'zangle must be between (90<zangle<=180deg)' 
               endif 
            ELSE 
               angle_path = angle 
                                                                        
               if ( ((angle_path-180.) .gt. 0.001) .or. ((angle_path-   &
               90.) .lt. -0.001) ) then                                 
                       write(*,*) '     For lambertian, surf_refl = l' 
               stop 'angle must be between (90<angle<=180 deg)' 
               endif 
            ENDIF 
!                                                                       
            secant_diffuse = 1.67 
            secant_path = 1. / cos(abs(angle_path-180.)*PI/180.) 
            diffuse_fac = secant_diffuse / secant_path 
            endif 
!                                                                       
!     *************************************************************     
!                                                                       
!  SAVE AIRMASS FACTORS FOR USE WITH MULTIPLE SCATTERING                
!                                                                       
            DO 10 IAIR = 1, MXLAY 
               AIRMAS(IAIR) = SECNTA(IAIR) 
   10       CONTINUE 
         ENDIF 
         LTNSAV = LTGNT 
         LH2SAV = LH2 
!                                                                       
         if (nmol.lt.mxmol) then 
         NMOL_1 = NMOL+1 
         DO 20 MOL = NMOL_1, MXMOL 
            DO 19 ILAYR = 1, NLAYRS 
               WKL(MOL,ILAYR) = 0. 
   19       continue 
   20    CONTINUE 
         endif 
!                                                                       
         RETURN 
      ENDIF 
!                                                                       
      REWIND LINFIL 
!                                                                       
      IF (ILAS.GT.0) THEN 
         DVI = XDV*DVL(1) 
         MM = 64 
         DVC = REAL(MM)*DVL(LAYER) 
         V2 = VLAS+DVC 
         V1 = VLAS-DVC 
!                                                                       
!PRT     WRITE(IPR,930) V1,V2                                           
!                                                                       
      ENDIF 
!                                                                       
!     SET UP LAYER BOUNDARY PARAMETERS                                  
!                                                                       
      ALTZL = ALTZ(LAYER-1) 
      ALTZU = ALTZ(LAYER) 
      PZL = PZ(LAYER-1) 
      PZU = PZ(LAYER) 
      TZL = TZ(LAYER-1) 
      TZU = TZ(LAYER) 
!                                                                       
      PAVE = PAVEL(LAYER) 
      TAVE = TAVEL(LAYER) 
      WBROAD = WBRODL(LAYER) 
!                                                                       
      DO 30 M = 1, NMOL 
         WK(M) = WKL(M,LAYER) 
   30 END DO 
      do m=1,mx_xs 
         WXM(M) = XAMNT(M,LAYER) 
      enddo 

! Isotopologue column density scaled by Hitran ratio
      DO I = 1,NISOTPL
         M = MOLNUM(I)
!         ISOTPL = ISOTPLNUM(I)
         DO J = 1, ISOTPL_NUM(M)
            WKI(M,J) = WKL_ISOTPL(M,J,LAYER) &
     &                 / ISOTPL_ABD(M,J)
         ENDDO
      ENDDO
!                                                                       
      H2OSLF = H2OSL(LAYER) 
      WTOT = WTOTL(LAYER) 
      ALBAR = ALBL(LAYER) 
      ADBAR = ADBL(LAYER) 
      AVBAR = AVBL(LAYER) 
      DV = DVL(LAYER) 
      IPATHL = IPTH(LAYER) 
      SECNTL = SECNTA(LAYER) 
      SECANT = SECNTA(LAYER) 
      IF (ALTZ(LAYER-1).GE.0.) THEN 
         ALTAV = ALTZ(LAYER-1)- HZ* LOG(.5*(1.+EXP(-(ALTZ(LAYER)-ALTZ(  &
         LAYER-1))/HZ)))                                                
      ELSE 
         ALTAV = ALTZ(LAYER) 
      ENDIF 
!                                                                       
      RETURN 
!                                                                       
  890 FORMAT (//,                                                       &
     &        '0 SAMPLE   =',F13.4,/,                                   &
     &        '0 ALFAL0   =',F13.4,/,                                   &
     &        '0 AVMASS   =',F13.4,/)                                   
  900 FORMAT ('1 ERROR IATM = 0 IAERSL GT 0 ') 
  905 FORMAT ('1') 
  910 FORMAT ('0 V1(CM-1) = ',F12.4,/'0 V2(CM-1) = ',F12.4,/,           &
     &        '0 SAMPLE   =',F13.4,/'0 DVSET    =',F13.6,/,             &
     &        '0 ALFAL0   =',F13.4,/'0 AVMASS   =',F13.4,/,             &
     &        '0 DPTMIN   =',1P,E13.4,13X,'  DPTFAC   =',0P,F13.6)      
  915 FORMAT ('0 V2-V1 .GT. 2020. ') 
  920 FORMAT ('0',1X,7('*'),' LAST ',I5,' MOLECULES ON LINFIL NOT ',    &
     &        'SELECTED')                                               
  925 FORMAT ('0',1X,53('*'),/,'0',1X,14('*'),I5,' MOLECULES ',         &
     &        'REQUESTED',2X,12('*'),/,2X,7('*'),2X,'ONLY ',I5,         &
     &        ' MOLECULES ON LINFIL',2X,11('*'),/,'0',1X,53('*'))       
  930 FORMAT (2(/),'  V1 RESET ',F10.3,'  V2 RESET ',F10.3) 
 1020 FORMAT (55X,'************************************************',/, &
     &        55X,'**                                            **',/, &
     &        55X,'**            DERIVATIVE CALCULATION          **',/, &
     &        55X,'**            ----------------------          **',/, &
     &        55X,'**   All molecular amounts were retained for  **',/, &
     &        55X,'**   molecular broadening purposes.           **',/, &
     &        55X,'**   Molecular amounts except those for       **',/, &
     &        55X,'**   which derivatives are to be calculated   **',/, &
     &        55X,'**   will now be zeroed.  Absorptance         **',/, &
     &        55X,'**   coefficients will be calculated in       **',/, &
     &        55X,'**   place of optical depths on TAPE10.       **',/, &
     &        55X,'**                                            **',/, &
     &        55X,'**   (See subroutine oppath in lblrtm.f)      **',/, &
     &        55X,'**                                            **',/, &
     &        55X,'************************************************',/) 
                                                                        
 1021 FORMAT (55X,'************************************************',/, &
     &        55X,'**                                            **',/, &
     &        55X,'**            DERIVATIVE CALCULATION          **',/, &
     &        55X,'**            ----------------------          **',/, &
     &        55X,'**   All molecular amounts were retained for  **',/, &
     &        55X,'**   molecular broadening purposes.           **',/, &
     &        55X,'**   Molecular amounts except those for       **',/, &
     &        55X,'**   which derivatives are to be calculated   **',/, &
     &        55X,'**   will now be zeroed.  Optical depth for   **',/, &
     &        55X,'**   the derivative species will be stored    **',/, &
     &        55X,'**   on TAPE10.                               **',/, &
     &        55X,'**                                            **',/, &
     &        55X,'**   (See subroutine oppath in lblrtm.f)      **',/, &
     &        55X,'**                                            **',/, &
     &        55X,'************************************************',/) 
 1022 FORMAT (55X,'************************************************',/, &
     &        55X,'**                                            **',/, &
     &        55X,'** Note that XS species were also zeroed...   **',/, &
     &        55X,'**                                            **',/, &
     &        55X,'************************************************',/) 
                                                                        
!                                                                       
end subroutine OPPATH
!                                                                       
!     -------------------------------------------------------------     
!                                                                       
      SUBROUTINE PATH 
!                                                                       
      USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,             &
     &                     MXISOTPL, ISOTPL_ABD, ISOTPL_NUM,             &
     &                     MXMOL, MX_XS, MXTRAC
      IMPLICIT REAL*8           (V) 
!                                                                       
!                                                                       
!     SUBROUTINE PATH INITIALIZES LINFIL AND INPUTS LAYER PARAMETERS    
!     SUBROUTINE PATH INPUTS AND OUTPUTS HEADER FROM LINFIL AND         
!     INPUTS AND OUTPUTS PATH PARAMETERS FOR EACH LAYER                 
!                                                                       
!      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MXZMD=6000,                   &
!     &     MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,mx_xs=38,MXTRAC=22) 
!                                                                       
      COMMON COMSTR(250,9) 
      COMMON R1(3600),R2(900),R3(225) 
!                                                                       
      COMMON /MANE/ P0,TEMP0,NLAYRS,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
     &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
     &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
     &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
     &              EXTID(10)                                           
      CHARACTER*8  EXTID 
                                                                        
      COMMON /BNDPRP/ TMPBND,BNDEMI(3),BNDRFL(3),IBPROP,surf_refl,pad_3,&
     &                angle_path,secant_diffuse,secant_path,diffuse_fac 
      common /profil_scal/ nmol_scal,hmol_scal(64),xmol_scal(64),       &
     &                     n_xs_scal,h_xs_scal(64),x_xs_scal(64)        
!                                                                       
      character*1 hmol_scal,h_xs_scal 
      character*1 surf_refl 
      character*3 pad_3 
!                                                                       
      common /lbl_geo/ zh1,zh2,zangle 
      character*8 hol_angle,blank_angle 
!                                                                       
      COMMON /MSACCT/ IOD,IDIR,ITOP,ISURF,MSPTS,MSPANL(MXLAY),          &
     &                MSPNL1(MXLAY),MSLAY1,ISFILE,JSFILE,KSFILE,        &
     &                LSFILE,MSFILE,IEFILE,JEFILE,KEFILE                
!                                                                       
      character*8      XID,       HMOLID,      YID 
      real*8               SECANT,       XALTZ 
!                                                                       
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
     &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
     &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF  
      COMMON /PATHD/ PAVEL(MXLAY),TAVEL(MXLAY),WKL(MXMOL,MXLAY),        &
     &               WBRODL(MXLAY),DVL(MXLAY),                          &
     &               WTOTL(MXLAY),ALBL(MXLAY),ADBL(MXLAY),              &
     &               AVBL(MXLAY),H2OSL(MXLAY),                          &
     &               IPTH(MXLAY),ITYL(MXLAY),SECNTA(MXLAY),             &
     &               HT1,HT2,ALTZ(0:MXLAY),                             &
     &               PZ(0:MXLAY),TZ(0:MXLAY)                            
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
     &              NLTEFL,LNFIL4,LNGTH4                                
!                                                                       
!                                                                       
!     IXMAX=MAX NUMBER OF X-SECTION MOLECULES, IXMOLS=NUMBER OF THESE   
!     MOLECULES SELECTED, IXINDX=INDEX VALUES OF SELECTED MOLECULES     
!     (E.G. 1=CLONO2), XAMNT(I,L)=LAYER AMOUNTS FOR I'TH MOLECULE FOR   
!     L'TH LAYER, ANALOGOUS TO AMOUNT IN /PATHD/ FOR THE STANDARD       
!     MOLECULES.                                                        
!                                                                       
      COMMON /PATHX/ IXMAX,IXMOLS,IXINDX(mx_xs),XAMNT(mx_xs,MXLAY) 
!                                                                       
!     COMMON BLOCKS AND PARAMETERS FOR THE PROFILES AND DENSITIES       
!     FOR THE CROSS-SECTION MOLECULES.                                  
!                                                                       
      CHARACTER*10 XSFILE,XSNAME,ALIAS 
      COMMON /XSECTF/ XSFILE(6,5,mx_xs),XSNAME(mx_xs),ALIAS(4,mx_xs) 
      COMMON /XSECTR/ V1FX(5,MX_XS),V2FX(5,MX_XS),DVFX(5,MX_XS),        &
     &                WXM(MX_XS),NTEMPF(5,MX_XS),NSPECR(MX_XS),         &
     &                IXFORM(5,MX_XS),XSMASS(MX_XS),XDOPLR(5,MX_XS),    &
     &                NUMXS,IXSBIN                                      
      COMMON /IODFLG/ DVOUT 
!                                                                       
      COMMON /PATH_ISOTPL/ ISOTPL,NISOTPL,                              &
     &                     ISOTPL_FLAG(MXMOL,MXISOTPL),                 &
     &                     ISOTPL_MAIN_FLAG(MXMOL),                     &
     &                     MOLNUM(MXMOL*MXISOTPL),                      &
     &                     ISOTPLNUM(MXMOL*MXISOTPL),                   &
     &                     WKI(MXMOL,MXISOTPL)             
      COMMON /OPPATH_ISOTPL/ WKL_ISOTPL(MXMOL,MXISOTPL,MXLAY)
!                                                                       
      INTEGER :: ISOTPL_HCODE(MXMOL*MXISOTPL)
      INTEGER :: INPTYP(MXMOL*MXISOTPL),MOLLST(MXMOL)
      INTEGER :: IDXM(MXMOL*MXISOTPL),IDXJ(MXMOL*MXISOTPL)
      REAL :: ISOTPL_ABD_SUM(MXMOL),ISOTPL_ABD_SUBSUM(MXMOL)
      REAL :: ISOTPL_AMNT(MXMOL,MXISOTPL,MXLAY)
      REAL :: ISOTPL_AMNT_SUM(MXMOL)
!                                                                       
      CHARACTER*20 HEAD20 
      CHARACTER*6 MOLID 
      COMMON /MOLNAM/ MOLID(0:MXMOL) 
      CHARACTER*2 MOLSTR 
      CHARACTER*7 HEAD7 
      CHARACTER*6 HOLN2 
      CHARACTER*5 HEAD5 
      CHARACTER*4 HEAD4,hedxs 
      CHARACTER*4 HT1HRZ,HT2HRZ,HT1SLT,HT2SLT,  ht1,ht2 
      CHARACTER*3 CINP,CINPX,CBLNK 
      DIMENSION FILHDR(2),AMOUNT(2),AMTSTR(2) 
      DIMENSION HEDXS(15),WMT(mxmol),SECL(MXFSC),WXT(mx_xs),WTOTX(MXLAY) 
      DIMENSION WIT(MXISOTPL),WTOTI(MXLAY)
!
      DIMENSION WDRAIR(MXLAY) 
!                                                                       
      EQUIVALENCE (XID(1),FILHDR(1)) 
      EQUIVALENCE (YID(10),LTNSAV) , (YID(9),LH1SAV) , (YID(8),LH2SAV) 
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),              &
     &            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),              &
     &            (FSCDID(5),IEMIT) , (FSCDID(7),IPLOT),                &
     &            (FSCDID(9),JRAD) , (FSCDID(11),IMRG),                 &
     &            (FSCDID(15),IATM)                                     
      EQUIVALENCE (AMOUNT(1),COMSTR(1,1)) , (AMTSTR(1),COMSTR(21,1)) 
!                                                                       
      DATA HOLN2 / ' OTHER'/ 
      DATA HT1HRZ / ' AT '/,HT2HRZ / ' KM '/,HT1SLT / ' TO '/,          &
     &     HT2SLT / ' KM '/                                             
      DATA CBLNK / '   '/ 
!                                                                       
      DATA I_2/2/, I_10/10/ 
!                                                                       
      IF (IHIRAC.EQ.0) RETURN 
!                                                                       
      iform = 1 
      ifrmx = 1 
!                                                                       
      ICNTNM = MOD(IXSCNT,I_10) 
      IXSECT = IXSCNT/10 
!                                                                       
      ISET = 0 
      IDVSET = 0 
      IF (DVSET.LT.0.) THEN 
         DVSET = -DVSET 
         ISET = 1 
      ENDIF 
!                                                                       
      DO M = 1, MXMOL 
         WMT(M) = 0. 
         WK(M) = 0. 
      ENDDO 
                                                                        
      DO M = 1, mx_xs 
         WXM(M) = 0. 
         WXT(M) = 0. 
      ENDDO 

      DO M = 1, MXISOTPL
         WIT(M) = 0. 
      ENDDO 
!                                                                        
      SUMN2 = 0. 
      ISTOP = 0 
!                                                                       
      SECNT0 = 1. 
      IPATHL = 1 
      LH1 = 1 
      LH2 = 1 
!                                                                       
!     Obtain atmospheric definition information                         
!                                                                       
      IF (IATM.EQ.0) THEN 
         READ (IRD,901) IFORM,NLAYRS,NMOL,SECNT0,HEAD20,ZH1,HEAD4,ZH2,  &
         HEAD5,hol_angle,HEAD7                                          
                                                                        
!       test if path angle was read in for lambertain surface reflection
         data blank_angle /'        '/ 
         if ((hol_angle.eq.blank_angle).and.(surf_refl.ne.'s')) then 
         stop 'must input value for zangle (Record 2.1)' 
         endif 
!                                                                       
         read (hol_angle,903) zangle 
!                                                                       
         IF (NMOL.EQ.0) NMOL = 7 
         IF (SECNT0.LT.0.) THEN 
            IPATHL = 1 
         ELSE 
            IPATHL = 3 
         ENDIF 
         SECNT0 = ABS(SECNT0) 
         WRITE (IPR,902) SECNT0,NLAYRS,NMOL,HEAD20,ZH1,HEAD4,ZH2,       &
         HEAD5,ZANGLE,HEAD7                                             
!                                                                       
!        Put H1, H2, and ANGLE into YID (ANGLE is needed for            
!        CHARTS multiple scattering calculation)                        
!                                                                       
         CALL YDIH1(ZH1,ZH2,ZANGLE,YID) 
!                                                                       
         IF (IHIRAC.EQ.9) THEN 
!                                                                       
            DO 20 M = 1, NMOL 
               READ (MOLID(M),905) HMOLID(M) 
   20       CONTINUE 
         ENDIF 
      ELSE 
         WRITE (IPR,907) SECNT0 
         JTYPE = 0 
         GO TO 50 
      ENDIF 
!                                                                       
      JTYPE = 0 
                                                                        
      DO 30 L = 1, NLAYRS 
!                                                                       
         IF (L.EQ.1) THEN 
            IF (IFORM.EQ.1) THEN 
               READ (IRD,910) PAVE,TAVE,SECNTK,CINP,IPTHRK,ALTZ(L-1),   &
               PZ(L-1),TZ(L-1),ALTZ(L),PZ(L),TZ(L)                      
            ELSE 
               READ (IRD,911) PAVE,TAVE,SECNTK,CINP,IPTHRK,ALTZ(L-1),   &
               PZ(L-1),TZ(L-1),ALTZ(L),PZ(L),TZ(L)                      
            ENDIF 
            IF (CINP.NE.CBLNK) WRITE (IPR,912) 
         ELSE 
            IF (IFORM.EQ.1) THEN 
               READ (IRD,915) PAVE,TAVE,SECNTK,CINP,IPTHRK, ALTZ(L),PZ( &
               L), TZ(L)                                                
            ELSE 
               READ (IRD,916) PAVE,TAVE,SECNTK,CINP,IPTHRK, ALTZ(L),PZ( &
               L), TZ(L)                                                
            ENDIF 
            IF ((CINP.EQ.CBLNK.AND.JTYPE.EQ.1).OR. (                    &
            CINP.NE.CBLNK.AND.JTYPE.EQ.0)) THEN                         
               WRITE (IPR,912) 
               STOP ' JTYPE ERROR IN PATH ' 
            ENDIF 
         ENDIF 
         IF (IPTHRK.EQ.0.AND.SECNT0.EQ.1.) THEN 
            READ (HT1HRZ,917) HT1 
            READ (HT2HRZ,917) HT2 
         ELSE 
            READ (HT1SLT,917) HT1 
            READ (HT2SLT,917) HT2 
         ENDIF 
         IF (CINP.NE.CBLNK) THEN 
            JTYPE = 1 
            READ (CINP,920) ITYL(L) 
         ENDIF 
!                                                                       
!        If TZ(L) = 0 then reset TZ(L) = TAVE to avoid errors when      
!        calculating EMB in SUBROUTINE EMIN                             
!                                                                       
         IF (TZ(L).EQ.0.) TZ(L) = TAVE 
!                                                                       
         PAVEL(L) = PAVE 
         TAVEL(L) = TAVE 
         SECANT = SECNT0 
         IF (SECNTK.GT.0.) SECANT = SECNTK 
         SECL(L) = SECANT 
         SECNTA(L) = SECANT 
         IF (IPTHRK.NE.0) IPATHL = IPTHRK 
         IPTH(L) = IPATHL 
         IF (SECANT.EQ.0.) STOP 'PATH; SECANT = 0' 
         IF (IFORM.EQ.1) THEN 
            READ (IRD,925) (WKL(M,L),M=1,7),WBRODL(L) 
            IF (NMOL.GT.7) READ (IRD,925) (WKL(M,L),M=8,NMOL) 
         ELSE 
            READ (IRD,927) (WKL(M,L),M=1,7),WBRODL(L) 
            IF (NMOL.GT.7) READ (IRD,927) (WKL(M,L),M=8,NMOL) 
         ENDIF 
!                                                                       
!     --------------------------------------------------------------    
!                                                                       
!                     MIXING RATIO INPUT                                
!                                                                       
!                                                                       
!     First calculate the column amount of dry air ("WDRAIR")           
!     Initialize WDNSTY to WBRODL(L) (always in column density)         
!     Determine if each molecule is in column density.                  
!        - if so, just add to WDNSTY                                    
!        - if not, add to WMXRAT                                        
!                                                                       
!     NOTE that if WKL is greater than one, then column density         
!               if WKL is less than one, then mixing ratio              
!                                                                       
         WDNSTY = WBRODL(L) 
         WMXRAT = 0.0 
         WDRAIR(L) = 0.0 
                                                                        
         DO 22 M = 2,NMOL 
            IF (WKL(M,L).GT.1) THEN 
               WDNSTY = WDNSTY + WKL(M,L) 
            ELSE 
               WMXRAT = WMXRAT + WKL(M,L) 
            ENDIF 
   22    CONTINUE 
                                                                        
!                                                                       
!        EXECUTE TESTS TO ENSURE ALL COMBINATION OF COLUMN DENSITIES    
!        AND MIXING RATIOS FOR EACH LAYER HAVE BEEN PROPERLY SPECIFIED. 
                                                                        
!        IF THE LAYER SUM OF MIXING RATIOS IS LESS THAN ONE (WHICH      
!        IT SHOULD BE, GIVEN THAT WBROAD CONTRIBUTES TO THE DRY AIR     
!        MIXING RATIO), THEN COMPUTE DRY AIR BY DIVIDING THE TOTAL      
!        MOLECULAR AMOUNTS GIVEN IN DENSITY BY THE FRACTION OF DRY      
!        AIR (MIXING RATIO) THOSE MOLECULES COMPRISE.                   
!                                                                       
!        IF THE LAYER SUM OF MIXING RATIOS IS GREATER THAN OR EQUAL     
!        TO ONE, THAN AN ERROR HAS OCCURRED, SO STOP THE PROGRAM.       
!        WBROAD IS ALWAYS LISTED IN COLUMN DENSITY, SO THE SUM OF       
!        THE GIVEN MIXING RATIOS MUST ALWAYS BE LESS THAN ONE.          
!                                                                       
                                                                        
         IF (WBRODL(L).LT.1.0 .AND. WBRODL(L).NE.0.0) THEN 
            WRITE(IPR,918) L 
            WRITE(*,918) L 
            STOP 
         ENDIF 
                                                                        
         IF (WDNSTY.EQ.0.0 .AND. WMXRAT.NE.0.0) THEN 
            WRITE(IPR,921) L,WDNSTY,WMXRAT 
            WRITE(*,921) L,WDNSTY,WMXRAT 
            STOP 'WMXRAT AND/OR WDNSTY NOT PROPERLY SPECIFIED IN PATH' 
         ENDIF 
                                                                        
         IF (WMXRAT.LT.1.0) THEN 
            WDRAIR(L) = WDNSTY/(1.0-WMXRAT) 
         ELSE 
            WRITE(IPR,921) L,WMXRAT, WDNSTY 
            WRITE(*,921) L,WMXRAT, WDNSTY 
            STOP 'WMXRAT EXCEEDS 1.0' 
         ENDIF 
                                                                        
         IF (WKL(1,L).LE.1.0 .AND. WKL(1,L) .NE. 0.0 .AND. WDRAIR(L)    &
         .EQ.0.0) THEN                                                  
            WRITE(IPR,921) L,WKL(1,L),WDRAIR(L) 
            WRITE(*,921) L,WKL(1,L),WDRAIR(L) 
            STOP 'WMXRAT NOT PROPERLY SPECIFIED IN PATH' 
         ENDIF 
                                                                        
!                                                                       
!     NOW CONVERT ALL OTHER MOLECULES WHICH MAY BE IN MIXING RATIO      
!     TO MOLECULAR DENSITY USING WDRAIR(L)                              
!                                                                       
         DO 25 M = 1,NMOL 
            IF (WKL(M,L).LT.1.) WKL(M,L) = WKL(M,L)*WDRAIR(L) 
   25    CONTINUE 
!                                                                       
!     --------------------------------------------------------------    
!                                                                       
   30 END DO 
!                                                                       
      IF (IATM.EQ.0.AND.IXSECT.GE.1) THEN 
         READ (IRD,930) IXMOLS,IXSBIN 
         XV1 = V1 
         XV2 = V2 
         CALL XSREAD (XV1,XV2) 
         WRITE (IPR,932) (I,XSNAME(I),I=1,IXMOLS) 
         READ (IRD,900) IFRMX,NLAYXS,IXMOL,SECNTX,HEDXS 
         IF (IXMOL.EQ.0) THEN 
            WRITE (IPR,935) IXMOL 
            STOP ' PATH - IXMOL 0 ' 
         ENDIF 
         IF (IXMOL.NE.IXMOLS) THEN 
            WRITE (IPR,937) IXMOL,IXMOLS 
            STOP ' PATH - IXMOL .NE. IXMOLS ' 
         ENDIF 
         IF (NLAYRS.NE.NLAYXS) THEN 
            WRITE (IPR,940) NLAYRS,NLAYXS 
            STOP ' PATH - NLAYRS .NE. NLAYXS ' 
         ENDIF 
         SECNTX = ABS(SECNTX) 
         WRITE (IPR,942) SECNTX,NLAYXS,IXMOLS,HEDXS 
!                                                                       
         DO 40 L = 1, NLAYXS 
!                                                                       
            IF (L.EQ.1) THEN 
               IF (IFRMX.EQ.1) THEN 
                  READ (IRD,910) PAVX,TAVX,SECKXS,CINPX,IPTHKX,ALTXB,   &
                  PZXB,TZXB,ALTXT,PZXT,TZXT                             
               ELSE 
                  READ (IRD,911) PAVX,TAVX,SECKXS,CINPX,IPTHKX,ALTXB,   &
                  PZXB,TZXB,ALTXT,PZXT,TZXT                             
               ENDIF 
            ELSE 
               IF (IFRMX.EQ.1) THEN 
                  READ (IRD,915) PAVX,TAVX,SECKXS,CINPX,IPTHKX,ALTXT,   &
                  PZXT,TZXT                                             
               ELSE 
                  READ (IRD,916) PAVX,TAVX,SECKXS,CINPX,IPTHKX,ALTXT,   &
                  PZXT,TZXT                                             
               ENDIF 
            ENDIF 
            IF (IFRMX.EQ.1) THEN 
               READ (IRD,925) (XAMNT(M,L),M=1,7),WBRODX 
               IF (IXMOL.GT.7) READ (IRD,925) (XAMNT(M,L),M=8,IXMOL) 
            ELSE 
               READ (IRD,927) (XAMNT(M,L),M=1,7),WBRODX 
               IF (IXMOL.GT.7) READ (IRD,927) (XAMNT(M,L),M=8,IXMOL) 
            ENDIF 
!                                                                       
!     --------------------------------------------------------------    
!                                                                       
!             MIXING RATIO INPUT FOR CROSS SECTIONAL MOLECULES          
!                                                                       
!                                                                       
!     The column amount of dry air ("WDRAIR") has already been          
!     calculated above, so just convert all cross sectional             
!     molecules which may be in mixing ratio to molecular density       
!     using WDRAIR(L)                                                   
!                                                                       
!     NOTE that if XAMNT is greater than one, then column density       
!               if XAMNT is less than one, then mixing ratio            
!                                                                       
            DO 35 M = 1,IXMOL 
               IF (WDRAIR(L).EQ.0.0 .AND. XAMNT(M,L).LT.1 .AND. XAMNT(M,&
               L).NE.0.0) THEN                                          
                  WRITE(IPR,921) L,XAMNT(M,L),WDRAIR(L) 
                  WRITE(*,921) L,XAMNT(M,L),WDRAIR(L) 
                  STOP 'XAMNT NOT PROPERLY SPECIFIED IN PATH' 
               ENDIF 
               IF (XAMNT(M,L).LT.1) XAMNT(M,L) = XAMNT(M,L)*WDRAIR(L) 
   35       CONTINUE 
!                                                                       
!     --------------------------------------------------------------    
!                                                                       
   40    CONTINUE 
      ENDIF 
!                                                                       
   50 continue 
                                                                        
      WRITE (IPR,945) XID,(YID(M),M=1,2) 
                                                                        
!_______________________________________________________________________
                                                                        
!     at this point scale profile if option selected                    
                                                                        
      if (nmol_scal.gt.0 .or. n_xs_scal.gt.0) then 
                                                                        
!  *** It should be noted that no attempt has been made to keep the     
!      mass in a given layer constant, i.e. level pressure nor retained 
                                                                        
!      obtain accumulated amounts by molecule                           
                                                                        
         do m = 1, mxmol 
            wmt(m) = 0. 
         enddo 
                                                                        
         do m = 1, nmol 
            do l = 1, nlayrs 
               wmt(m) = wmt(m) + wkl(m,l) 
            enddo 
         enddo 
                                                                        
         wsum_brod = 0. 
         do l = 1, nlayrs 
            wsum_brod = wsum_brod + wbrodl(l) 
         enddo 
                                                                        
!        obtain dry air sum                                             
!             check to see if nitrogen is included in the selected molec
                                                                        
         if (wmt(22).ge.1.) then 
            wsum_drair = 0.0 
         else 
            wsum_drair = wsum_brod 
         endif 
                                                                        
         do m = 2, nmol 
            wsum_drair = wsum_drair + wmt(m) 
         enddo 
                                                                        
         write (ipr,*) 
         write (ipr,*) '   ',                                           &
     &         '******************************************************' 
         write (ipr,*) 
         write (ipr,*) '               Profile Scaling          ' 
                                                                        
         write (ipr,956) 
  956    format (/,4x,' molecule',                                      &
     &        2x,'hmol_scale',3x, 'xmol_scal_in',3x, 'scale factor',/)  
                                                                        
         do m = 1, nmol_scal 
            xmol_scal_m = xmol_scal(m) 
            xmol_scal(m) = -999. 
            if (hmol_scal(m).eq.' ') xmol_scal(m) = 1. 
            if (hmol_scal(m).eq.'0') xmol_scal(m) = 0. 
            if (hmol_scal(m).eq.'1') xmol_scal(m) = xmol_scal_m 
                                                                        
            if (hmol_scal(m).eq.'C' .or. hmol_scal(m).eq.'c')           &
     &           xmol_scal(m) = xmol_scal_m/wmt(m)                      
                                                                        
            if (hmol_scal(m).eq.'M' .or. hmol_scal(m).eq.'m')           &
     &           xmol_scal(m) = xmol_scal_m/(wmt(m)/wsum_drair)         
                                                                        
            if ((hmol_scal(m).eq.'P' .or. hmol_scal(m).eq.'p')          &
     &                                                .and. m.eq.1)     &
     &           xmol_scal(m) = (xmol_scal_m/2.99150e-23)/wmt(m)        
!                value from vpayne 2006/07/24                           
                                                                        
            if (hmol_scal(m).eq.'D' .or. hmol_scal(m).eq.'d')           &
     &           xmol_scal(m) =  (xmol_scal_m*2.68678e16)/wmt(m)       
                                                                        
            if ((hmol_scal(m).eq.'P' .or. hmol_scal(m).eq.'p')          &
     &                                                .and. m.ne.1) then
                                                                        
               write (ipr,*) 'm = ', m 
               stop ' (hmol_scal(m).eq."P" .and. m.ne.1) ' 
            endif 
                                                                        
            write (ipr,957) m, hmol_scal(m), xmol_scal_m, xmol_scal(m) 
  957       format (5x,i5,9x,a1,5x,1p, 4e15.7) 
                                                                        
           if (xmol_scal(m) .lt. -998.) then 
               write (ipr,*) 'm = ', m,' h_mol_scal(m) not valid ' 
               write (  *,*) 'm = ', m,' h_mol_scal(m) not valid ' 
               stop 
            endif 
                                                                        
            do l = 1, nlayrs 
               wkl(m,l) = wkl(m,l) * xmol_scal(m) 
            enddo 
                                                                        
            wmt(m) = wmt(m)*xmol_scal(m) 
                                                                        
         enddo 
                                                                        
         write (ipr,*) 
         write (ipr,*) '   ',                                           &
     &           '*****************************************************'
         write (ipr,*) 
                                                                        
      endif 
!                                                                       
!     at this point scale cross section profiles if option selected     
                                                                        
      if (n_xs_scal.gt.0 .and. ixsect.ge.1) then 
                                                                        
!  *** It should be noted that no attempt has been made to keep the     
!      mass in a given layer constant, i.e. level pressure nor retained 
                                                                        
!      obtain accumulated amounts by cross section species              
                                                                        
         do m = 1, n_xs_scal 
            wmt(m) = 0. 
            do l = 1, nlayrs 
               wmt(m) = wmt(m) + xamnt(m,l) 
            enddo 
         enddo 
                                                                        
         write (ipr,958) 
  958    format (/,4x,'  species',                                      &
     &        2x,'h_xs_scale',3x, 'x_xs_scal_in',3x, 'scale factor',/)  
                                                                        
         do m = 1, n_xs_scal 
            x_xs_scal_m = x_xs_scal(m) 
            x_xs_scal(m) = -999. 
            if (h_xs_scal(m).eq.' ') x_xs_scal(m) = 1. 
            if (h_xs_scal(m).eq.'0') x_xs_scal(m) = 0. 
            if (h_xs_scal(m).eq.'1') x_xs_scal(m) = x_xs_scal_m 
                                                                        
            if (h_xs_scal(m).eq.'C' .or. h_xs_scal(m).eq.'c')           &
     &           x_xs_scal(m) = x_xs_scal_m/wmt(m)                      
                                                                        
            if (h_xs_scal(m).eq.'M' .or. h_xs_scal(m).eq.'m')           &
     &           x_xs_scal(m) = x_xs_scal_m/(wmt(m)/wsum_drair)         
                                                                        
            if (h_xs_scal(m).eq.'D' .or. h_xs_scal(m).eq.'d')           &
     &           x_xs_scal(m) =  (x_xs_scal_m*2.68678e16)/wmt(m)       
                                                                        
            if ((h_xs_scal(m).eq.'P' .or. h_xs_scal(m).eq.'p')          &
     &                                           .and. m.ne.1) then     
               write (ipr,*) 'm = ', m 
               stop ' (h_xs_scal(m).eq."P" .and. m.ne.1) ' 
            endif 
                                                                        
           if (x_xs_scal(m) .lt. -998.) then 
               write (ipr,*) 'm = ', m,' h_xs_scal(m) not valid ' 
               write (  *,*) 'm = ', m,' h_xs_scal(m) not valid ' 
               stop 
            endif 
                                                                        
            write (ipr,959) m, h_xs_scal(m), x_xs_scal(m), x_xs_scal_m 
  959       format (5x,i5,9x,a1,5x,1p, 4e15.7) 
                                                                        
            x_xs_scal(m) = x_xs_scal_m 
                                                                        
            do l = 1, nlayrs 
               xamnt(m,l) = xamnt(m,l) * x_xs_scal(m) 
            enddo 
                                                                        
         enddo 
                                                                        
         wmt(m) = wmt(m) * x_xs_scal(m) 
                                                                        
         write (ipr,*) 
         write (ipr,*) '   ',                                           &
     &         '******************************************************' 
         write (ipr,*) 
                                                                        
!        reset the accumulated array to zero                            
                                                                        
         do m = 1,mxmol 
            wmt(m) = 0. 
         enddo 
                                                                        
      endif 
!_______________________________________________________________________
!     end of profile scaling                                            

!                                                                       
!     Obtain isotopologue information                         
!                                                                       
      IF (ISOTPL.EQ.1) THEN 
         READ (IRD,928) NISOTPL
         READ (IRD,928) (ISOTPL_HCODE(I),I=1,NISOTPL)
         MOLCOUNT = 1
         DO I = 1,NISOTPL
            MOLNUM(I) = FLOOR(0.1*ISOTPL_HCODE(I))
            ISOTPLNUM(I) = INT((0.1*ISOTPL_HCODE(I)-MOLNUM(I))*10.)
            IF (ISOTPLNUM(I).EQ.0) ISOTPLNUM(I)=10
            ISOTPL_FLAG(MOLNUM(I),ISOTPLNUM(I)) = 1
            ISOTPL_MAIN_FLAG(MOLNUM(I)) = 1
            IF (I.EQ.1) THEN
               MOLLST(MOLCOUNT) = MOLNUM(I)
            ELSE 
               IF (MOLNUM(I).NE.MOLNUM(I-1)) THEN
                  MOLCOUNT = MOLCOUNT + 1
                  IF (MOLCOUNT.GT.MXMOL) THEN
                     STOP 'MOLCOUNT GREATER THAN MXMOL IN PATH'
                  ENDIF
                  MOLLST(MOLCOUNT) = MOLNUM(I)
               ENDIF
            ENDIF
         ENDDO
         READ (IRD,900) IFRMI, NLAYIS
         DO L = 1,NLAYIS
            IF (IFRMI.EQ.1) THEN
               READ (IRD,925) (ISOTPL_AMNT(MOLNUM(I),ISOTPLNUM(I),L),&
     &                      I=1,NISOTPL)
            ELSE
               READ (IRD,927) (ISOTPL_AMNT(MOLNUM(I),ISOTPLNUM(I),L),&
     &                      I=1,NISOTPL)
            ENDIF
!     Check for negative input, which is permitted to identify when amounts
!     for isotopologues not directly specified as column amount are to be
!     set to their original default values (i.e. WKL * HITRAN ratio).
            DO I = 1,NISOTPL
               IF (ISOTPL_AMNT(MOLNUM(I),ISOTPLNUM(I),L).LT.0.) THEN 
                  ISOTPL_AMNT(MOLNUM(I),ISOTPLNUM(I),L) = &
     &              WKL(MOLNUM(I),L) * ISOTPL_ABD(MOLNUM(I),ISOTPLNUM(I))
               ENDIF
            ENDDO

!     Check input for consistent type for each specified molecular species;
            DO I = 1,NISOTPL
               INPTYP(I) = 0
               IF (ISOTPL_AMNT(MOLNUM(I),ISOTPLNUM(I),L).LE.1.) THEN 
                  INPTYP(I) = 1
               ENDIF
            ENDDO
            DO I = 1,NISOTPL-1
               DO J = I+1,NISOTPL
                 IF (MOLNUM(I).EQ.MOLNUM(J).AND.INPTYP(I).NE.INPTYP(J)) THEN 
                   WRITE(IPR,921) L,ISOTPL_AMNT(MOLNUM(I),ISOTPLNUM(I),L), &
     &                              ISOTPL_AMNT(MOLNUM(J),ISOTPLNUM(J),L)
                   WRITE(*,921) L,ISOTPL_AMNT(MOLNUM(I),ISOTPLNUM(I),L), &
     &                              ISOTPL_AMNT(MOLNUM(J),ISOTPLNUM(J),L)
                   WRITE(MOLSTR,'(I2)') MOLNUM(I)
                  STOP 'ISOTPL_AMNT INPUT TYPES DIFFERENT'
                ENDIF
               ENDDO
            ENDDO
         ENDDO
!     Get sum of isotopologue Hitran ratios for each molecule
! h
         DO I=1,MXMOL
            ISOTPL_ABD_SUM(I) = 0.0
            DO J=1,ISOTPL_NUM(I)
               ISOTPL_ABD_SUM(I) = ISOTPL_ABD_SUM(I) + ISOTPL_ABD(I,J)
            ENDDO
         ENDDO
!                                                                       
!     --------------------------------------------------------------    
!                                                                       
!             COLUMN AMOUNT INPUT FOR ISOTOPOLOGUE MOLECULES  
!                                                                       
!             OR CONVERSION FROM FRACTION TO COLUMN AMOUNT
!                                                                       
!                                                                       
!     Isotopologue amounts are specified either as column amounts
!     or as fractions (in terms of moleulces, not mass) for that
!     species.  For column amount input, the total column amount
!     for that species is allowed to increase.  For fractional 
!     input, fractions of all isotopologues are adjusted to force 
!     the total column amount for that species to remain fixed, 
!     and all isotopologues are activated for that species.
!     Adjusted fractions are used to define the column amount for
!     each isotopologue from the column amount of the primary 
!     molecule amount input earlier. 
!                                                                       
!     NOTE that if ISOTPL_AMNT greater than one, then column density,
!                  ISOTPL_AMNT is less than or equal to one, then fraction
!                                                                       
         DO L = 1,NLAYIS

            DO I = 1,MOLCOUNT
               M = MOLLST(I)
!               ISOTPL = ISOTPLNUM(I)
! Generate sums needed for fractional input
               ISOTPL_ABD_SUBSUM(M) = ISOTPL_ABD_SUM(M)
               ISOTPL_AMNT_SUM(M) = 0.0
               DO J = 1, ISOTPL_NUM(M)
                  IF (ISOTPL_FLAG(M,J).EQ.1 .AND. &
     &                ISOTPL_AMNT(M,J,L).LE.1.) THEN
! p
                      ISOTPL_ABD_SUBSUM(M) = ISOTPL_ABD_SUBSUM(M) - &
     &                                       ISOTPL_ABD(M,J)
! s
                      ISOTPL_AMNT_SUM(M) = ISOTPL_AMNT_SUM(M) + &
     &                                     ISOTPL_AMNT(M,J,L)  
                  ENDIF
               ENDDO

! Test to check sum of specified fractions less than sum of Hitran ratios
               IF (ISOTPL_AMNT_SUM(M) .GT. ISOTPL_ABD_SUM(M)) THEN 
                  WRITE(IPR,'(2I3,2E15.7)') L, I, ISOTPL_AMNT_SUM(M), &
     &                           ISOTPL_ABD_SUM(M) 
                  WRITE(*,'(2I3,2E15.7)') L, I, ISOTPL_AMNT_SUM(M), &
     &                         ISOTPL_ABD_SUM(M)
                  STOP 'ISOTPL_AMNT NOT PROPERLY SPECIFIED IN PATH; ' &
     &                //'SUM EXCEEDS SUM OF HITRAN RATIOS'
               ENDIF 

               DO J = 1, ISOTPL_NUM(M)
! Generate isotopologue column amounts
! Isotopologue column density input for selected molecules/isotopologues
! Total column amount for parent species increases
                  IF (ISOTPL_AMNT(M,J,L).GT.1. .AND. &
     &                ISOTPL_FLAG(M,J).EQ.1) THEN
                    WKL_ISOTPL(M,J,L) = ISOTPL_AMNT(M,J,L)
                  ENDIF

! Isotopologue fraction input
! Total column amount for parent species remains constant
! Isotopologue fraction scaled by Hitran ratio and converted to column density
                  IF (ISOTPL_AMNT(M,J,L).LE.1.) THEN
! Selected molecules/isotopologues; unmodified fractions
                    IF (ISOTPL_FLAG(M,J).EQ.1) THEN
! g
                      ADJ_FRAC = ISOTPL_AMNT(M,J,L)
                    ENDIF

! Unselected molecules/isotopologues; modified fractions
                    IF (ISOTPL_FLAG(M,J).EQ.0) THEN
! Test to prevent divide by zero, though code is not likely to reach this 
! point in that condition
                      IF (ISOTPL_ABD_SUBSUM(M) .EQ. 0.0) THEN 
                        WRITE(IPR,'(2I3,E15.7)') L, I, &
     &                           ISOTPL_ABD_SUBSUM(M) 
                        WRITE(*,'(2I3,E15.7)') L, I, &
     &                           ISOTPL_ABD_SUBSUM(M) 
                        STOP 'ISOTPL_ABD_SUBSUM IS ZERO IN PATH'
                      ENDIF 
! g
                      ADJ_FRAC = &
     &                    (ISOTPL_ABD_SUM(M) - ISOTPL_AMNT_SUM(M)) &
     &                  * (ISOTPL_ABD(M,J) / ISOTPL_ABD_SUBSUM(M))
                    ENDIF

! Generate isotopologue column density with fractions, g, from
! column density for primary molecule input earlier
                    WKL_ISOTPL(M,J,L) = WKL(M,L) * ADJ_FRAC

                  ENDIF

               ENDDO

            ENDDO

         ENDDO

! Activate all isotopologues for molecules with selected isotopologues
! provided as fractional input
         DO I = 1,NISOTPL
            M = MOLNUM(I)
            IF (MAXVAL(ISOTPL_FLAG(M,1:MXISOTPL)).GT.0.AND.&
     &          INPTYP(I).EQ.1) THEN
               DO J = 1, ISOTPL_NUM(M)
                  ISOTPL_FLAG(M,J) = 1 
               ENDDO
            ENDIF
         ENDDO
! Get total number of activated isotoplogues
         ISOCOUNT = SUM(ISOTPL_FLAG)
! Update list of Hitran codes for all activated isotopologues
         IDXCNT = 0
         DO I = 1,MOLCOUNT
            M = MOLLST(I)
            DO J = 1, ISOTPL_NUM(M)
               IF (ISOTPL_FLAG(M,J).EQ.1) THEN
                  IDXCNT = IDXCNT + 1
                  ISOTPL_HCODE(IDXCNT) = M*10+J
                  IF (J.EQ.10) ISOTPL_HCODE(IDXCNT) = M*10+J-10
                  IDXM(IDXCNT) = M
                  IDXJ(IDXCNT) = J
               ENDIF
            ENDDO
         ENDDO

      ENDIF                                                                        
!                                                                       
!     --------------------------------------------------------------    
!                                                                       
      IF (IFORM.EQ.1) THEN 
         WRITE (IPR,950) 
      ELSE 
         WRITE (IPR,951) 
      ENDIF 
!                                                                       
      DO 80 L = 1, NLAYRS 
         IF (IATM.GT.0.) SECL(L) = 1.0 
!                                                                       
         DO 60 M = 1, NMOL 
            WKL(M,L) = WKL(M,L)*SECL(L) 
   60    CONTINUE 
         IF (IXSECT.GE.1) THEN 
            DO 70 M = 1, IXMOLS 
               XAMNT(M,L) = XAMNT(M,L)*SECL(L) 
   70       CONTINUE 
         ENDIF 
         WBRODL(L) = WBRODL(L)*SECL(L) 
         SECL(L) = 1.0 
   80 END DO 
!                                                                       
!     LTGNT = TOP ANTERIOR LAYER FOR TANGENT VIEWING CASE               
!     LTGNT = NLAYER FOR SPACE-TO-SPACE PATH OR IPATHL = 1 OR 3         
!     OTHERWISE, LTGNT = LAYER WHICH INCLUDES H1                        
!                              (I.E. OBSERVATION POINT)                 
!                                                                       
      LTGNT = NLAYRS 
      ITCNT = IPTH(1) 
      DO 90 ILAYR = 1, NLAYRS 
         IF (IPTH(ILAYR).NE.ITCNT) LTGNT = ILAYR-1 
         IF (IPTH(ILAYR).EQ.2) THEN 
            LH1 = ILAYR 
            LH2 = ILAYR 
         ENDIF 
         IF (IPTH(ILAYR).EQ.3) THEN 
            LH2 = ILAYR 
         ELSE 
            LH1 = ILAYR 
         ENDIF 
         IF (IPTH(ILAYR).EQ.1) LH1 = ILAYR 
         ITCNT = IPTH(ILAYR) 
   90 END DO 
!                                                                       
      LTNSAV = LTGNT 
      LH1SAV = LH1 
      LH2SAV = LH2 
!                                                                       
      DV = 0. 
      PWTD = 0. 
      TWTD = 0. 
      WTOT = 0. 
      PWTX = 0. 
      TWTX = 0. 
      WTOX = 0. 
      PWTI = 0. 
      TWTI = 0. 
      WTOI = 0. 
!                                                                       
!     Write message if IOD=2 (Optical depth flag) and IMRG = 1          
!                                                                       
      IF (IOD.EQ.2.AND.IMRG.EQ.1) THEN 
         WRITE(IPR,953) 
      ENDIF 
!                                                                       
!     LOOP OVER LAYERS                                                  
!                                                                       
      do m = 1,nmol 
         wmt(m) = 0. 
      enddo 
                                                                        
      DO 130 L = 1, NLAYRS 
         IPROB = 0 
         FACTOR = 1. 
         IF ((IPTH(L).EQ.2).AND.(IANT.EQ.0)) FACTOR = 2. 
         SUMWK = 0. 
         DO 100 M = 1, NMOL 
            SUMWK = SUMWK+WKL(M,L) 
            WMT(M) = WMT(M)+WKL(M,L)*FACTOR 
  100    CONTINUE 
         WTOTL(L) = SUMWK+WBRODL(L) 
         SUMN2 = SUMN2+WBRODL(L)*FACTOR 
         WTOT = WTOT+WTOTL(L)*FACTOR 
         PWTD = PWTD+PAVEL(L)*WTOTL(L)*FACTOR 
         TWTD = TWTD+TAVEL(L)*WTOTL(L)*FACTOR 
         FRH2O = WKL(1,L)/WTOTL(L) 
         ALFCOR = (PAVEL(L)/P0)*SQRT(TEMP0/TAVEL(L)) 
!                                                                       
!     CROSS SECTIONS                                                    
!                                                                       
         IF (IXSECT.GE.1) THEN 
            SUMXK = 0. 
            DO 110 M = 1, IXMOLS 
               SUMXK = SUMXK+XAMNT(M,L) 
               WXT(M) = WXT(M)+XAMNT(M,L)*FACTOR 
  110       CONTINUE 
            WTOTX(L) = SUMXK+WBRODL(L) 
            WTOX = WTOX+WTOTX(L)*FACTOR 
            PWTX = PWTX+PAVEL(L)*WTOTX(L)*FACTOR 
            TWTX = TWTX+TAVEL(L)*WTOTX(L)*FACTOR 
         ENDIF 
!                                                                       
!     ISOTOPOLOGUES
!                                                                       
         IF (ISOTPL.GE.1) THEN 
            SUMIK = 0. 
!            DO I = 1, NISOTPL
!               SUMIK = SUMIK+WKL_ISOTPL(MOLNUM(I),ISOTPLNUM(I),L) 
!               WIT(I) = WIT(I)+WKL_ISOTPL(MOLNUM(I),ISOTPLNUM(I),L)*&
!     &                  FACTOR 
!            ENDDO
            DO I = 1, ISOCOUNT
               SUMIK = SUMIK+WKL_ISOTPL(IDXM(I),IDXJ(I),L) 
               WIT(I) = WIT(I)+WKL_ISOTPL(IDXM(I),IDXJ(I),L)*&
     &                  FACTOR 
            ENDDO
            WTOTI(L) = SUMIK+WBRODL(L) 
            WTOI = WTOI+WTOTI(L)*FACTOR 
            PWTI = PWTI+PAVEL(L)*WTOTI(L)*FACTOR 
            TWTI = TWTI+TAVEL(L)*WTOTI(L)*FACTOR 
         ENDIF 
!                                                                       
!     CORRECT FOR WATER SELF BROADENING                                 
!                                                                       
         H2OSLF = (1.-FRH2O+5.*FRH2O) 
         H2OSL(L) = H2OSLF 
         ALBAR = ALFAL0*ALFCOR*H2OSLF 
         ALBL(L) = ALBAR 
!                                                                       
!     3.58115E-07 = SQRT( 2.* LOG(2.)*AVOGAD*BOLTZ/(CLIGHT*CLIGHT) )    
!                                                                       
         ADBAR = 3.58115E-07*(0.5*(V1+V2))*SQRT(TAVEL(L)/AVMASS) 
         ADBL(L) = ADBAR 
         AVBAR = 0.5*(ALBAR+SQRT(ALBAR*ALBAR+4.*ADBAR*ADBAR)) 
!                                                                       
         AVBL(L) = AVBAR 
!                                                                       
         OLDDV = DV 
         DV = AVBAR/SAMPLE 
         IF (IHIRAC.EQ.2) DV = ALBAR/SAMPLE 
         IF (IHIRAC.EQ.3) DV = ADBAR/SAMPLE 
!                                                                       
!     Skip to next layer if IOD=2 (OPTICAL DEPTH FLAG FOR EXACT DV)     
!     and IMRG = 1, or if IOD=1 and DVOUT nonzero (OPTICAL DEPTH FLAG   
!     FOR INTERPOLATED DV).  This bypasses ITYPE assignment from one    
!     layer to the next.                                                
!                                                                       
         IF ( (IOD.EQ.2.AND.IMRG.EQ.1) .OR. (IOD.EQ.1.) ) THEN 
            DVL(L) = DV 
            TYPE = 0. 
            ITYPE = -99 
            ZETA = ALBL(L)/(ALBL(L)+ADBL(L)) 
            IF (IFORM.EQ.1) THEN 
               WRITE (IPR,960) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),    &
               TAVEL(L),ALBL(L),ADBL(L),AVBL(L),ZETA, DV,H2OSL(L),DV,   &
               TYPE,ITYPE,IPTH(L),SECL(L)                               
            ELSE 
               WRITE (IPR,961) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),    &
               TAVEL(L),ALBL(L),ADBL(L),AVBL(L),ZETA, DV,H2OSL(L),DV,   &
               TYPE,ITYPE,IPTH(L),SECL(L)                               
            ENDIF 
            GOTO 130 
         ENDIF 
!                                                                       
         IF (DV.LT.DVSET.AND.ISET.EQ.1) THEN 
            DV = OLDDV 
            IF (L.EQ.1.AND.DV.EQ.0.) DV = DVSET 
            IDVSET = IDVSET+1 
            IF (IDVSET.EQ.1) WRITE (IPR,955) DV,ALTZ(L-1),HT1,ALTZ(L),  &
            HT2                                                         
         ENDIF 
!                                                                       
         DVC = DV 
         TYPE = 0. 
         ITYPE = 99 
         IF (L.EQ.1) THEN 
!                                                                       
!     DV IS ASSUMED TO BE .LT. 1                                        
!     SET DV TO 3 SIGNIFICANT FIGURES                                   
!                                                                       
            ISCAL = LOG10(DV)-3. 
            SCAL = 10.**ISCAL 
            IDV = (DV/SCAL)+0.5 
!                                                                       
!     SET IDV TO BE EVEN                                                
!                                                                       
            IF (MOD(IDV,I_2).GT.0) IDV = IDV+1 
            DV = SCAL* REAL(IDV) 
!                                                                       
         ELSE 
!                                                                       
            IF (JTYPE.EQ.1) GO TO 120 
            TYPE = OLDDV/DV 
            TYPMAX = 2.5 
            IF (TYPE.GT.TYPMAX) THEN 
               IPROB = 1 
               ISTOP = 1 
            ELSEIF (TYPE.GE.1.2) THEN 
!                                                                       
!     TYPE IS BETWEEN 1.2 AND TYPMAX                                    
!                                                                       
               DV = OLDDV 
               ITYPE = 1./(TYPE-1.)+0.5 
               IF (ITYPE.EQ.3) ITYPE = 2 
               DV = OLDDV* REAL(ITYPE)/ REAL(ITYPE+1) 
            ELSEIF (TYPE.GE.0.8) THEN 
!                                                                       
!     TYPE IS BETWEEN 0.8 AND 1.2 (SET TO 1.0)                          
!                                                                       
               DV = OLDDV 
               ITYPE = 0 
            ELSE 
!                                                                       
!     TYPE IS LESS THAN 0.8                                             
!                                                                       
               DV = OLDDV 
               ITYPE = 0 
               IF (IEMIT.EQ.0) THEN 
                  ITYPE = TYPE/(1.-TYPE)+0.5 
                  DV = DV* REAL(ITYPE+1)/ REAL(ITYPE) 
                  ITYPE = -ITYPE 
               ENDIF 
            ENDIF 
         ENDIF 
!                                                                       
         DVL(L) = DV 
!                                                                       
         ITYL(L) = ITYPE 
!                                                                       
  120    IF (JTYPE.EQ.1) THEN 
            IF (ITYL(L).NE.99) THEN 
!                                                                       
               IF (ITYL(L).EQ.0) THEN 
                  DV = OLDDV 
               ELSE 
                  DV = OLDDV* REAL(ITYL(L))/ REAL(ITYL(L)+1) 
               ENDIF 
               DVL(L) = DV 
!                                                                       
            ENDIF 
         ENDIF 
         ZETA = ALBAR/(ALBAR+ADBAR) 
!                                                                       
         DV = DVL(L) 
!                                                                       
         IF (iform.EQ.1) THEN 
            WRITE (IPR,960) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),       &
            TAVEL(L), ALBL(L),ADBL(L),AVBL(L),ZETA,DVC,H2OSL(L),        &
            DVL(L),TYPE,ITYL(L),IPTH(L),SECL(L)                         
         ELSE 
            WRITE (IPR,961) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),       &
            TAVEL(L), ALBL(L),ADBL(L),AVBL(L),ZETA,DVC,H2OSL(L),        &
            DVL(L),TYPE,ITYL(L),IPTH(L),SECL(L)                         
         ENDIF 
         IF (IPROB.GT.0) WRITE (IPR,962) TYPMAX 
  130 END DO 
!                                                                       
!     SKIP TO END WHEN USING EXACT CALCULATED DV FOR OPTICAL            
!     DEPTH CALCULATIONS (IOD = 2,IMRG = 1)                             
!                                                                       
      IF (IOD.EQ.2.AND.IMRG.EQ.1) GOTO 142 
!                                                                       
      PWTD = PWTD/WTOT 
      TWTD = TWTD/WTOT 
      IF (IXSECT.GE.1) THEN 
         PWTX = PWTX/WTOX 
         TWTX = TWTX/WTOX 
      ENDIF 
      IF (ISOTPL.GE.1) THEN 
         PWTI = PWTI/WTOI 
         TWTI = TWTI/WTOI
      ENDIF 
!
      IF (ISTOP.EQ.1) WRITE (IPR,965) 
      IF (ISTOP.EQ.1) STOP 'PATH; ISTOP EQ 1' 
!                                                                       
!     If DVOUT is nonzero (IOD=1 or 4:  interpolate optical depths to   
!     value of DVOUT), then test to be sure that DVOUT is finer than    
!     the monochromatic DV (and thus ensuring enough monochromatic point
!     are available to reach the V2 endpoint for the interpolated       
!     spectrum).                                                        
                                                                        
      IF (IOD.EQ.1.AND.DV.LT.DVOUT) THEN 
         WRITE (IPR,968) DVOUT,DV 
         endif 
!                                                                       
! Special case:  if imrg=40-43, the user is doing analytic jacobians    
!    ---> non zero dvset implies the user wants a specific value,       
!         (this is tested for above)                                    
!         so use that one rather than re-setting to dv                  
!         (for this case at this point in the code, dvset>0 and         
!          dvout=dvset, so nothing needs to be done)                    
!                                                                       
         IF (IOD.ge.3) THEN 
            if ((imrg.eq.1).or.((imrg.ge.40).and.(imrg.le.43))) then 
            if (iset.eq.0) then 
            dvout = dv 
            write(ipr,*) 'new dvout, dvset, and dv: ', dvout,dvset,dv 
            endif 
         ENDIF 
      ENDIF 
!                                                                       
!     If DVSET is nonzero (set the final layer DV to the value of       
!     DVSET), then test to be sure that DVSET is not more than          
!     20% different than the monochromatic DV.                          
!                                                                       
      IF (DVSET.GT.0.) THEN 
         RATIO = 1. 
         IF (DVSET.GT.0.) RATIO = DVSET/DV 
         IF (ISET.EQ.0) THEN 
            IF (RATIO.GT.1.2.OR.RATIO.LT.0.8) THEN 
               WRITE (IPR,967) RATIO,DVSET,DV 
               STOP 'PATH; RATIO ERROR' 
            ENDIF 
         ENDIF 
         WRITE (IPR,945) XID,(YID(M),M=1,2) 
         WRITE (IPR,948) 
         IF (IFORM.EQ.1) THEN 
            WRITE (IPR,950) 
         ELSE 
            WRITE (IPR,951) 
         ENDIF 
         DO 140 L = 1, NLAYRS 
            ALBAR = ALBL(L) 
            ADBAR = ADBL(L) 
            AVBAR = AVBL(L) 
!                                                                       
            OLDDV = DV 
            DV = DVL(L) 
            DVC = DV 
            DV = DVL(L)*RATIO 
            TYPE = 0. 
            IF (L.GT.1) TYPE = OLDDV/DV 
            ITYPE = ITYL(L) 
            DVL(L) = DV 
            ZETA = ALBL(L)/(ALBL(L)+ADBL(L)) 
            IF (IFORM.EQ.1) THEN 
               WRITE (IPR,960) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),    &
               TAVEL(L),ALBL(L),ADBL(L),AVBL(L),ZETA, DVC,H2OSL(L),DV,  &
               TYPE,ITYPE,IPTH(L),SECL(L)                               
            ELSE 
               WRITE (IPR,961) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),    &
               TAVEL(L),ALBL(L),ADBL(L),AVBL(L),ZETA, DVC,H2OSL(L),DV,  &
               TYPE,ITYPE,IPTH(L),SECL(L)                               
            ENDIF 
  140    CONTINUE 
      ENDIF 
!                                                                       
  142 CONTINUE 
!                                                                       
      IF (NLAYRS.LT.5) THEN 
         WRITE (IPR,970) 
      ELSE 
         WRITE (IPR,945) XID,(YID(M),M=1,2) 
      ENDIF 
!                                                                       
!     --------------------------------------------------------------    
!                                                                       
!     Write out column densities for molecules to TAPE6                 
!                                                                       
!                                                                       
      IF (IFORM.EQ.1) THEN 
         WRITE (IPR,974) (HMOLID(I),I=1,7),HOLN2 
         DO 150 L = 1, NLAYRS 
            WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),       &
            TAVEL(L), IPTH(L),(WKL(M,L),M=1,7),WBRODL(L)                
  150    CONTINUE 
         IF (NLAYRS.GT.1) THEN 
            WRITE (IPR,985) 
            L = NLAYRS 
            WRITE (IPR,990) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTD,TWTD,        &
            (WMT(M),M=1,7),SUMN2                                        
         ENDIF 
      ELSE 
         WRITE (IPR,975) (HMOLID(I),I=1,7),HOLN2 
         DO 151 L = 1, NLAYRS 
            WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),       &
            TAVEL(L), IPTH(L),(WKL(M,L),M=1,7),WBRODL(L)                
  151    CONTINUE 
         IF (NLAYRS.GT.1) THEN 
            WRITE (IPR,985) 
            L = NLAYRS 
            WRITE (IPR,991) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTD,TWTD,        &
            (WMT(M),M=1,7),SUMN2                                        
         ENDIF 
      ENDIF 
!                                                                       
      IF (NMOL.GT.7) THEN 
         DO 170 MLO = 8, NMOL, 8 
            MHI = MLO+7 
            MHI = MIN(MHI,NMOL) 
            IF (NLAYRS.LT.5) THEN 
               WRITE (IPR,970) 
            ELSE 
               WRITE (IPR,945) XID,(YID(M),M=1,2) 
            ENDIF 
            IF (IFORM.EQ.1) THEN 
               WRITE (IPR,974) (HMOLID(I),I=MLO,MHI) 
               DO 160 L = 1, NLAYRS 
                  WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L), &
                  TAVEL(L),IPTH(L),(WKL(M,L),M=MLO,MHI)                 
  160          CONTINUE 
               IF (NLAYRS.GT.1) THEN 
                  WRITE (IPR,985) 
                  L = NLAYRS 
                  WRITE (IPR,990) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTD,TWTD,  &
                  (WMT(M),M=MLO,MHI)                                    
               ENDIF 
            ELSE 
               WRITE (IPR,975) (HMOLID(I),I=MLO,MHI) 
               DO 161 L = 1, NLAYRS 
                  WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L), &
                  TAVEL(L),IPTH(L),(WKL(M,L),M=MLO,MHI)                 
  161          CONTINUE 
               IF (NLAYRS.GT.1) THEN 
                  WRITE (IPR,985) 
                  L = NLAYRS 
                  WRITE (IPR,991) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTD,TWTD,  &
                  (WMT(M),M=MLO,MHI)                                    
               ENDIF 
            ENDIF 
!                                                                       
  170    CONTINUE 
      ENDIF 
!                                                                       
!     --------------------------------------------------------------    
!                                                                       
!     Write out mixing ratios for molecules to TAPE6 in either          
!     15.7 format (IFORM = 1) or 10.4 format (IFORM = 0).               
!                                                                       
!           Reset WDRAIR(L) for each layer                              
!           (WKL(M,L) now in column density)                            
!                                                                       
!                                                                       
      IF (IFORM.EQ.1) THEN 
         WRITE (IPR,976) (HMOLID(I),I=1,7),HOLN2 
         DO 172 L = 1, NLAYRS 
            WDRAIR(L) = WBRODL(L) 
            DO 171 M = 2,NMOL 
               WDRAIR(L) = WDRAIR(L) + WKL(M,L) 
  171       CONTINUE 
            IF (WDRAIR(L).EQ.0.0) THEN 
               WRITE(IPR,979) 
            ELSE 
               WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),    &
               TAVEL(L), IPTH(L),(WKL(M,L)/WDRAIR(L),M=1,7),WBRODL(L)   
            ENDIF 
  172    CONTINUE 
      ELSE 
         WRITE (IPR,977) (HMOLID(I),I=1,7),HOLN2 
         DO 174 L = 1, NLAYRS 
            WDRAIR(L) = WBRODL(L) 
            DO 173 M = 2,NMOL 
               WDRAIR(L) = WDRAIR(L) + WKL(M,L) 
  173       CONTINUE 
            IF (WDRAIR(L).EQ.0.0) THEN 
               WRITE(IPR,979) 
            ELSE 
               WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),    &
               TAVEL(L), IPTH(L),(WKL(M,L)/WDRAIR(L),M=1,7),WBRODL(L)   
            ENDIF 
  174    CONTINUE 
      ENDIF 
!                                                                       
!                                                                       
                                                                        
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
                  IF (WDRAIR(L).EQ.0.0) THEN 
                     WRITE(IPR,979) 
                  ELSE 
                     WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,       &
                     PAVEL(L),TAVEL(L),IPTH(L), (WKL(M,L)/WDRAIR(L),M=  &
                     MLO,MHI)                                           
                  ENDIF 
  176          CONTINUE 
            ELSE 
               WRITE (IPR,977) (HMOLID(I),I=MLO,MHI) 
               DO 177 L = 1, NLAYRS 
                  IF (WDRAIR(L).EQ.0.0) THEN 
                     WRITE(IPR,979) 
                  ELSE 
                     WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,       &
                     PAVEL(L),TAVEL(L),IPTH(L), (WKL(M,L)/WDRAIR(L),M=  &
                     MLO,MHI)                                           
                  ENDIF 
  177          CONTINUE 
            ENDIF 
  178    CONTINUE 
      ENDIF 
!                                                                       
!     --------------------------------------------------------------    
!                                                                       
!     Write out column densities for cross sectional molecules to TAPE6 
!                                                                       
!                                                                       
      IF (IXSECT.GE.1) THEN 
         DO 190 MLO = 1, IXMOLS, 8 
            MHI = MLO+7 
            MHI = MIN(MHI,IXMOLS) 
            IF (NLAYRS.LT.5.AND.MLO.NE.1) THEN 
               WRITE (IPR,970) 
            ELSE 
               WRITE (IPR,995) XID,(YID(M),M=1,2) 
            ENDIF 
            IF (IFRMX.EQ.1) THEN 
               WRITE (IPR,974) (XSNAME(I),I=MLO,MHI) 
               DO 180 L = 1, NLAYRS 
                  WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L), &
                  TAVEL(L),IPTH(L),(XAMNT(M,L),M=MLO,MHI)               
  180          CONTINUE 
               IF (NLAYRS.GT.1) THEN 
                  WRITE (IPR,985) 
                  L = NLAYRS 
                  WRITE (IPR,990) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTX,TWTX,  &
                  (WXT(M),M=MLO,MHI)                                    
               ENDIF 
            ELSE 
               WRITE (IPR,975) (XSNAME(I),I=MLO,MHI) 
               DO 181 L = 1, NLAYRS 
                  WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L), &
                  TAVEL(L),IPTH(L),(XAMNT(M,L),M=MLO,MHI)               
  181          CONTINUE 
               IF (NLAYRS.GT.1) THEN 
                  WRITE (IPR,985) 
                  L = NLAYRS 
                  WRITE (IPR,991) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTX,TWTX,  &
                  (WXT(M),M=MLO,MHI)                                    
               ENDIF 
            ENDIF 
  190    CONTINUE 
!                                                                       
!     --------------------------------------------------------------    
!                                                                       
!        Write out mixing ratios for cross sectional                    
!             molecules to TAPE6                                        
!                                                                       
!                                                                       
         DO 198 MLO = 1, IXMOLS, 8 
            MHI = MLO+7 
            MHI = MIN(MHI,IXMOLS) 
            IF (IFRMX.EQ.1) THEN 
               WRITE (IPR,976) (XSNAME(I),I=MLO,MHI) 
               DO 195 L = 1, NLAYRS 
                  WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L), &
                  TAVEL(L),IPTH(L),(XAMNT(M,L)/WDRAIR(L),M=MLO,MHI)     
  195          CONTINUE 
            ELSE 
               WRITE (IPR,977) (XSNAME(I),I=MLO,MHI) 
               DO 196 L = 1, NLAYRS 
                  WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L), &
                  TAVEL(L),IPTH(L),(XAMNT(M,L)/WDRAIR(L),M=MLO,MHI)     
  196          CONTINUE 
            ENDIF 
  198    CONTINUE 

      ENDIF 
!                                                                       
!     --------------------------------------------------------------    
!                                                                       
!                                                                       
!     Write out column densities for isotolologue molecules to TAPE6                 
!                                                                       
!                                                                       
      IF (ISOTPL.GE.1) THEN

      IF (ISOCOUNT.LE.8) THEN 
      IF (NLAYIS.LT.5) THEN 
          WRITE (IPR,970) 
      ELSE 
          WRITE (IPR,945) XID,(YID(M),M=1,2) 
      ENDIF 
      IF (IFORM.EQ.1) THEN 
         WRITE (IPR,996) (ISOTPL_HCODE(I),I=1,ISOCOUNT)
         DO L = 1, NLAYIS 
            WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),       &
            TAVEL(L),IPTH(L),(WKL_ISOTPL(IDXM(I),IDXJ(I),L),     &
     &      I=1,ISOCOUNT)
         ENDDO
         IF (NLAYIS.GT.1) THEN 
            WRITE (IPR,985) 
            L = NLAYIS 
            WRITE (IPR,990) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTI,TWTI,        &
            (WIT(I),I=1,ISOCOUNT)
         ENDIF 
      ELSE 
         WRITE (IPR,997) (ISOTPL_HCODE(I),I=1,ISOCOUNT)
         DO L = 1, NLAYIS 
            WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),       &
            TAVEL(L),IPTH(L),(WKL_ISOTPL(IDXM(I),IDXJ(I),L),     &
     &      I=1,ISOCOUNT)
         ENDDO
         IF (NLAYIS.GT.1) THEN 
            WRITE (IPR,985) 
            L = NLAYIS 
            WRITE (IPR,991) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTI,TWTI,        &
            (WIT(I),I=1,ISOCOUNT)
         ENDIF 
      ENDIF 
      ENDIF 
!                                                                       
      IF (ISOCOUNT.GT.8) THEN 
         DO MLO = 1, ISOCOUNT, 8 
            MHI = MLO+7 
            MHI = MIN(MHI,ISOCOUNT) 
            IF (NLAYIS.LT.5) THEN 
               WRITE (IPR,970) 
            ELSE 
               WRITE (IPR,945) XID,(YID(M),M=1,2) 
            ENDIF 
            IF (IFORM.EQ.1) THEN 
               WRITE (IPR,996) (ISOTPL_HCODE(I),I=MLO,MHI)
               DO L = 1, NLAYIS 
                  WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L), &
                  TAVEL(L),IPTH(L),(WKL_ISOTPL(IDXM(I),IDXJ(I),L)&
     &            ,I=MLO,MHI)                 
               ENDDO
               IF (NLAYIS.GT.1) THEN 
                  WRITE (IPR,985) 
                  L = NLAYIS 
                  WRITE (IPR,990) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTI,TWTI,  &
                  (WIT(I),I=MLO,MHI)                                    
               ENDIF 
            ELSE 
               WRITE (IPR,997) (ISOTPL_HCODE(I),I=MLO,MHI)
               DO L = 1, NLAYIS 
                  WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L), &
                  TAVEL(L),IPTH(L),(WKL_ISOTPL(IDXM(I),IDXJ(I),L)&
     &            ,I=MLO,MHI)                 
               ENDDO
               IF (NLAYIS.GT.1) THEN 
                  WRITE (IPR,985) 
                  L = NLAYIS 
                  WRITE (IPR,991) L,ALTZ(0),HT1,ALTZ(L),HT2,PWTI,TWTI,  &
                  (WIT(I),I=MLO,MHI)                                    
               ENDIF 
            ENDIF 
!                                                                       
         ENDDO
      ENDIF 
!                                                                       
!     --------------------------------------------------------------    
!                                                                       
!     Write out mixing ratios for isotopologue molecules to TAPE6 in either          
!     15.7 format (IFORM = 1) or 10.4 format (IFORM = 0).               
!                                                                       
!           Reset WDRAIR(L) for each layer                              
!           Convert WKL_ISOTPL(M,ISO,L) to mixing ratio
!                                                                       
!                                                                       
      IF (ISOCOUNT.LE.8) THEN
      IF (IFORM.EQ.1) THEN 
         WRITE (IPR,998) (ISOTPL_HCODE(I),I=1,ISOCOUNT)
         DO L = 1, NLAYIS 
            WDRAIR(L) = WBRODL(L) 
            DO M = 2,NMOL
               WDRAIR(L) = WDRAIR(L) + WKL(M,L) 
            ENDDO
            IF (WDRAIR(L).EQ.0.0) THEN 
               WRITE(IPR,979) 
            ELSE 
               WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),    &
               TAVEL(L),IPTH(L),(WKL_ISOTPL(IDXM(I),IDXJ(I),L)          &
     &         /WDRAIR(L),I=1,ISOCOUNT)                                           
            ENDIF 
         ENDDO
      ELSE 
         WRITE (IPR,999) (ISOTPL_HCODE(I),I=1,ISOCOUNT)
         DO L = 1, NLAYIS 
            WDRAIR(L) = WBRODL(L) 
            DO M = 2,NMOL 
               WDRAIR(L) = WDRAIR(L) + WKL(M,L) 
            ENDDO
            IF (WDRAIR(L).EQ.0.0) THEN 
               WRITE(IPR,979) 
            ELSE 
               WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,PAVEL(L),    &
               TAVEL(L),IPTH(L),(WKL_ISOTPL(IDXM(I),IDXJ(I),L)          &
     &         /WDRAIR(L),I=1,ISOCOUNT)                                           
            ENDIF 
         ENDDO
      ENDIF 
      ENDIF 
!                                                                       
!                                                                       
      IF (ISOCOUNT.GT.8) THEN 
         DO MLO = 1, ISOCOUNT, 8 
            MHI = MLO+7 
            MHI = MIN(MHI,ISOCOUNT) 
            IF (NLAYIS.LT.5) THEN 
               WRITE (IPR,970) 
            ELSE 
               WRITE (IPR,945) XID,(YID(M),M=1,2) 
            ENDIF 
            IF (IFORM.EQ.1) THEN 
               WRITE (IPR,998) (ISOTPL_HCODE(I),I=MLO,MHI)
               DO L = 1, NLAYIS 
                  WDRAIR(L) = WBRODL(L) 
                  DO M = 2,NMOL
                     WDRAIR(L) = WDRAIR(L) + WKL(M,L) 
                  ENDDO
                  IF (WDRAIR(L).EQ.0.0) THEN 
                     WRITE(IPR,979) 
                  ELSE 
                     WRITE (IPR,980) L,ALTZ(L-1),HT1,ALTZ(L),HT2,       &
                     PAVEL(L),TAVEL(L),IPTH(L),(WKL_ISOTPL(IDXM(I),     &
     &               IDXJ(I),L)/WDRAIR(L),I=MLO,MHI)                    
                  ENDIF 
               ENDDO
            ELSE 
               WRITE (IPR,999) (ISOTPL_HCODE(I),I=MLO,MHI)
               DO L = 1, NLAYIS 
                  WDRAIR(L) = WBRODL(L) 
                  DO M = 2,NMOL
                     WDRAIR(L) = WDRAIR(L) + WKL(M,L) 
                  ENDDO
                  IF (WDRAIR(L).EQ.0.0) THEN 
                     WRITE(IPR,979) 
                  ELSE 
                     WRITE (IPR,982) L,ALTZ(L-1),HT1,ALTZ(L),HT2,       &
                     PAVEL(L),TAVEL(L),IPTH(L),(WKL_ISOTPL(IDXM(I),     &
     &               IDXJ(I),L)/WDRAIR(L),I=MLO,MHI)                    
                  ENDIF 
               ENDDO
            ENDIF 
         ENDDO
      ENDIF 

      ENDIF 

      RETURN 
!                                                                       
  900 FORMAT (1X,I1,I3,I5,F10.2,15A4) 
  901 FORMAT (1X,I1,I3,I5,F10.2,A20,F8.2,A4,F8.2,A5,A8,A7) 
  902 FORMAT ('0 SECANT   =',F13.4,/'0 NLAYRS=',I4,/'0 NMOL=',I4,/'0',  &
     &        A20,F8.2,A4,F8.2,A5,F8.3,A7)                              
  903 format (f8.3) 
  905 FORMAT (A6) 
  907 FORMAT ('0 SECANT   =',F13.4) 
  910 FORMAT (E15.7,F10.4,F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2)) 
  911 FORMAT (3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2)) 
  912 FORMAT ('0   ********* ITYPE(L) IS SET FROM INPUT ******** ') 
  915 FORMAT (E15.7,F10.4,F10.4,A3,I2,23X,(F7.2,F8.3,F7.2)) 
  916 FORMAT (3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2)) 
  917 FORMAT (A4) 
  918 FORMAT ('0 WBROD FOR LAYER=',I3,' MUST BE SPECIFIED               &
     &       IN COLUMN DENSITY')                                        
  920 FORMAT (I3) 
  921 FORMAT (I3,2E15.7) 
  922 FORMAT (A16,2I5) 
!                                                                       
  925 FORMAT (8E15.7) 
  927 FORMAT (8E10.3) 
  928 FORMAT (10I5) 
!                                                                       
  930 FORMAT (I5,5X,I5) 
  932 FORMAT (/,'  THE CROSS-SECTION MOLECULES SELECTED ARE: ',/,/,(5X, &
     &        I5,3X,A))                                                 
  935 FORMAT (/,'***** IXMOL = ',I5,' *****',/) 
  937 FORMAT (/,'***** IXMOL = ',I5,' .NE. IXMOLS = ',I5,' *****',/) 
  940 FORMAT (/,'***** NLAYRS = ',I5,' .NE. NLAYXS = ',I5,' *****',/) 
  942 FORMAT (/,'0 SECANTX  =',F13.4,/'0 NLAYXS=',I4,/'0 ISMOLS=',I4,/, &
     &        '0',15A4)                                                 
  945 FORMAT ('1'/'0',10A8,2X,2(1X,A8,1X)) 
  948 FORMAT ('0   ****** DVSET is set from input ******') 
  950 FORMAT ('0','LAYER',26X,'P(MB)',7X,'T(K)',4X,'ALPHL',4X,'ALPHD',  &
     &        4X,'ALPHV',3X,'ZETA',2X,'CALC DV',2X,'H2OSLF',5X,'DV',5X, &
     &        'TYPE',' ITYPE IPATH ',3X,'SECANT'/)                      
  951 FORMAT ('0','LAYER',25X,'P(MB)',3X,'T(K)',4X,'ALPHL',4X,'ALPHD',  &
     &        4X,'ALPHV',3X,'ZETA',2X,'CALC DV',2X,'H2OSLF',5X,'DV',5X, &
     &        'TYPE',' ITYPE IPATH ',3X,'SECANT'/)                      
  953 FORMAT ('0 ***** EXACT CALCULATED DV USED IN CALCULATION *****') 
  955 FORMAT (/,'0 **** CALC DV WAS RESET TO PREVIOUS DV',F12.6,/,      &
     &        '  AT ALT=  ',2(F7.3,A3),' AND ABOVE')                    
  960 FORMAT ('0',I5,2(F7.3,A3),1P,E15.7,0P,F8.2,3F9.6,F6.3,            &
     &        F9.6,F8.4,F9.6,F7.3,2I5,F12.6)                            
  961 FORMAT ('0',I5,2(F7.3,A3),   F10.4,   F8.2,3F9.6,F6.3,            &
     &        F9.6,F8.4,F9.6,F7.3,2I5,F12.6)                            
  962 FORMAT (20X,'  DV RATIO  .GT. ',F10.2) 
  965 FORMAT (/,20X,'  TYPE GT 2.5') 
  967 FORMAT ('  RATIO ERROR ',F10.3,'  DVSET = ',F10.4,'  DV=',F10.4) 
  968 FORMAT ('  DVOUT MUST BE < DV ','  DVOUT = ',E10.4,'  DV=',E10.4) 
  970 FORMAT (////) 
  974 FORMAT ('0',53X,'MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',/,32X,  &
     &        'P(MB)',6X,'T(K)',3X,'IPATH',5X,8(A10,4X))                
!     format (  60X,3(A10,4X))                                          
  975 FORMAT ('0',53X,'MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',/,29X,  &
     &        'P(MB)',6X,'T(K)',3X,'IPATH',1X,8(1X,A6,3X))              
  976 FORMAT (/,'1',54X,'----------------------------------',           &
     &         /,'0',60X,'MIXING RATIOS BY LAYER ',/,32X,               &
     &        'P(MB)',6X,'T(K)',3X,'IPATH',5X,8(A10,4X))                
  977 FORMAT (/,'1',54X,'----------------------------------',           &
     &         /,'0',60X,'MIXING RATIOS BY LAYER ',/,29X,               &
     &        'P(MB)',6X,'T(K)',3X,'IPATH',1X,8(1X,A6,3X))              
  979 FORMAT (/,'0','  MIXING RATIO IS UNDEFINED. DRYAIR DENSITY=0.0') 
  980 FORMAT ('0',I3,2(F7.3,A3),F15.7,F9.2,I5,2X,1P,8E15.7) 
!     format (    (55X,1P,8E15.7,0P))                                   
  982 FORMAT ('0',I3,2(F7.3,A3),F12.5,F9.2,I5,2X,1P,8E10.3,0P) 
  985 FORMAT ('0',54X,'ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH') 
  990 FORMAT ('0',I3,2(F7.3,A3),F15.7,F9.2,   7X,1P,8E15.7) 
!     format (    (55X,1P,8E15.7,0P))                                   
  991 FORMAT ('0',I3,2(F7.3,A3),F12.5,F9.2,7X,1P,8E10.3,0P) 
  995 FORMAT ('1'/'0',10A8,2X,2(1X,A8,1X),/,/,'0',53X,                  &
     &        '     *****  CROSS SECTIONS  *****      ')                
  996 FORMAT (/,'0',40X,'ISOTOPOLOGUE MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',/  &
     &        ,32X,'P(MB)',6X,'T(K)',3X,'IPATH',5X,8(I10,4X))                
  997 FORMAT (/,'0',40X,'ISOTOPOLOGUE MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',/  &
     &        ,29X,'P(MB)',6X,'T(K)',3X,'IPATH',1X,8(1X,I6,3X))              
  998 FORMAT (/,'1',54X,'----------------------------------',           &
     &        /,'0',47X,'ISOTOPOLOGUE MIXING RATIOS BY LAYER ',/,32X,   &
     &        'P(MB)',6X,'T(K)',3X,'IPATH',5X,8(I10,4X))                
  999 FORMAT (/,'1',54X,'----------------------------------',           &
     &        /,'0',47X,'ISOTOPOLOGUE MIXING RATIOS BY LAYER ',/,29X,   &
     &        'P(MB)',6X,'T(K)',3X,'IPATH',1X,8(1X,I6,3X))              
!
 1000 FORMAT ('Layer',I2,': Changing molecule ',I2,' from ',E10.3,      &
     &          ' to 1.000E+20.')                                       
 1001 FORMAT (' ************************************************',/     &
     &        '  ERROR in SUBROUTINE PATH: Sum of mixing ratios ',/     &
     &        '         greater than or equal to 1.0:           ',/     &
     &        '       Layer #:',I3,/                                    &
     &        '       Total mixing ratio: ',E10.3,/                     &
     &        '       Total column density: ',E10.3,/                   &
     &        ' ************************************************')      
 1010 FORMAT (2F7.3,E15.7,F10.4) 
 1011 FORMAT (' ************************************************',/     &
     &        '    ERROR: Sum of molecular densities in layer  ',/      &
     &        '           equal to zero; Layer #: ',I3,/                &
     &        ' ************************************************')      
!                                                                       
end subroutine PATH
      BLOCK DATA BOPDPT 
      COMMON /CONVF/ CHI(251),RDVCHI,RECPI,ZSQBND,A3,B3,JCNVF4 
!                                                                       
      DATA JCNVF4 / 0 / 
!                                                                       
end block data BOPDPT
!                                                                       
!     -------------------------------------------------------------     
!                                                                       
      SUBROUTINE OPDPTH (MPTS) 
!                                                                       
      USE phys_consts, ONLY: radcn2
      USE lblparams,   ONLY: N_ABSRB
      IMPLICIT REAL*8           (V) 
!                                                                       
!      parameter (n_absrb=5050) 
                                                                        
!     OPDPTH CALLS CONTNM,LINF4,HIRAC1,NONLTE                           
!                                                                       
      COMMON /MANE/ P0,TEMP0,NLAYER,DVXM,H2OSLF,WTOT,ALBAR,ADBAR,AVBAR, &
     &              AVFIX,LAYRFX,SECNT0,SAMPLE,DVSET,ALFAL0,AVMASS,     &
     &              DPTMIN,DPTFAC,ALTAV,AVTRAT,TDIFF1,TDIFF2,ALTD1,     &
     &              ALTD2,ANGLE,IANT,LTGNT,LH1,LH2,IPFLAG,PLAY,TLAY,    &
     &              EXTID(10)                                           
      CHARACTER*8  EXTID 
!                                                                       
      character*8      XID,       HMOLID,      YID 
      real*8               SECANT,       XALTZ 
!                                                                       
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
     &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
     &                EMISIV,FSCDID(17),NMOL,LAYRS ,YI1,YID(10),LSTWDF  
      COMMON /LASIV/ VLAS,ILAS 
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb) 
      COMMON /SCATTR/ V1SC,V2SC,DVSC,NPTSC,SCTTR(n_absrb) 
!                                                                       
      COMMON /IODFLG/ DVOUT 
      COMMON /RCNTRL/ ILNFLG 
      COMMON /LBLF/ V1R4,V2R4,DVR4,NPTR4,BOUND4,R4(2502),RR4(2502) 
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
     &              NLTEFL,LNFIL4,LNGTH4                                
      COMMON /CONVF/ CHI(251),RDVCHI,RECPI,ZSQBND,A3,B3,JCNVF4 
!                                                                       
      COMMON /CNTSCL/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL 
!                                                                       
      EQUIVALENCE (FSCDID(1),IHIRAC) , (FSCDID(2),ILBLF4),              &
     &            (FSCDID(3),IXSCNT) , (FSCDID(4),IAERSL),              &
     &            (FSCDID(5),IEMIT) , (FSCDID(7),IPLOT),                &
     &            (FSCDID(8),IPATHL) , (FSCDID(9),JRAD),                &
     &            (FSCDID(11),IMRG)                                     
!                                                                       
!     DATA JCNVF4 / 0 /                                                 
!                                                                       
      DATA I_10/10/ 
!                                                                       
      CALL CPUTIM (TIME0) 
!                                                                       
      ICNTNM = MOD(IXSCNT,I_10) 
      IXSECT = IXSCNT/10 
!                                                                       
      IEMST = IEMIT 
      IEMIT = 0 
      IPFLAG = 0 
      DPTMST = DPTMIN 
      IF (IEMST.EQ.0.AND.IPATHL.EQ.2) DPTMIN = 2.*DPTMST 
!                                                                       
!     PRINT LAYER INFORMATION                                           
!                                                                       
      IF (NOPR.EQ.0) THEN 
         IF (IMRG.LE.10) WRITE (IPR,900) 
         WRITE (IPR,905) LAYRS 
         IF (ILAS.GT.0) WRITE (IPR,910) VLAS,V1,V2 
         WRITE (IPR,915) XID,(YID(M),M=1,2),TIME0 
      ENDIF 
!                                                                       
      XKT = TAVE/RADCN2 
!                                                                       
!     JRAD= -1  NO RADIATION TERM IN ABSORPTION COEFFICIENTS            
!     JRAD=  0  RADIATION TERM PUT IN BY PANEL                          
!     JRAD=  1  RADIATION TERM INCLUDED IN LINE STRENGTHS               
!                                                                       
      IF (((V1/XKT).LT. 5.).AND.(JRAD.NE.-1)) JRAD = 0 
!                                                                       
!     DVABS IS USED AS A FLAG IN SUBSEQUENT PROGRAMS                    
!                                                                       
      DVABS = 0. 
      IF (ICNTNM.NE.0) THEN 
         DVABS = 1. 
         V1ABS = INT(V1) 
         IF (V1.LT.0.) V1ABS = V1ABS-1. 
         V1ABS = V1ABS-3.*DVABS 
         V2ABS = INT(V2+3.*DVABS+0.5) 
         NPTABS = (V2ABS-V1ABS)/DVABS+1.5 
         IF (PAVE.LE.0.5) IPFLAG = 1 
         DO 10 I = 1, n_absrb 
            ABSRB(I) = 0. 
   10    CONTINUE 
         CALL CONTNM (JRAD) 
      ENDIF 
      DVR4 = 0. 
!                                                                       
      IF (ILBLF4.GE.1) THEN 
         ALFAV = SAMPLE*DV 
         ALFAV4 = 64.*ALFAV 
!     Read in DVR4 from REJ file                                        
         IF (ILNFLG.EQ.2) THEN 
            READ(16) LAYRS, DVR4 
         ELSE 
!     Compute DVR4                                                      
            DVR4 = ALFAV4/SAMPLE 
            IF (ILNFLG.EQ.1) WRITE(16) LAYRS, DVR4 
         ENDIF 
         BOUND4 = 25. 
         IF (ILBLF4.EQ.2.AND.IPFLAG.EQ.1) BOUND4 = 5. 
         IPTS4 = BOUND4/DVR4 
!                                                                       
         IF (NOPR.EQ.0) WRITE (IPR,920) IPTS4,DVR4,BOUND4 
!                                                                       
         REWIND LINFIL 
         REWIND LNFIL4 
         V1R4 = V1-2.*DVR4 
         V2R4 = V2+2.*DVR4 
         V1L4 = V1R4-BOUND4-DVR4 
         V2L4 = V2R4+BOUND4+2*DVR4 
         IF ((IHIRAC.EQ.1).OR.(IHIRAC.EQ.9)) CALL LINF4 (V1L4,V2L4) 
      ENDIF 
!                                                                       
!    Write out DV to REJ1 file                                          
      IF (ILNFLG.EQ.1) WRITE(15) LAYRS, DV 
!    Read in DV from  REJ1 file                                         
      IF (ILNFLG.EQ.2) READ(15) LAYRS,DV 
!                                                                       
      IF (IHIRAC.EQ.1) CALL HIRAC1 (MPTS) 
      IF (IHIRAC.EQ.4) CALL NONLTE (MPTS) 
      IF (IHIRAC.EQ.9) CALL HIRAC1 (MPTS) 
!                                                                       
      IEMIT = IEMST 
      DPTMIN = DPTMST 
      CALL CPUTIM (TIME1) 
      TIMEO = TIME1-TIME0 
      WRITE (IPR,925) TIME1,TIMEO 
!                                                                       
      RETURN 
!                                                                       
  900 FORMAT ('1') 
  905 FORMAT ('0 LAYER = ',I8) 
  910 FORMAT ('0 VLAS  ',F20.8,8X,'V1 RESET ',F12.5,8X,'V2 RESET ',     &
     &        F12.5)                                                    
  915 FORMAT ('0',10A8,2X,2(1X,A8,1X),/,'0 TIME ENTERING OPDPTH ',      &
     &        F15.3)                                                    
  920 FORMAT ('0  IPTS4 FOR LINF4 = ',I5,3X,' DV FOR LINF4 = ',F10.5,   &
     &        5X,'BOUND FOR LINF4 =',F10.4)                             
  925 FORMAT ('0 TIME LEAVING OPDPTH ',F15.3,'  TOTAL FOR LAYER ',      &
     &        F15.3)                                                    
!                                                                       
end subroutine OPDPTH
!                                                                       
!     -------------------------------------------------------------     
!                                                                       
      SUBROUTINE READEM(ICOEF) 
!                                                                       
!     Reads in emission function values directly from file "EMISSIVITY" 
!                                                                       
      USE lblparams, ONLY: NMAXCO 
      IMPLICIT REAL*8           (V) 
!                                                                       
!     ----------------------------------------------------------------  
!     Parameter and common blocks for direct input of emission          
!     function values                                                   
!                                                                       
!      PARAMETER (NMAXCO=4040) 
      COMMON /EMSFIN/ V1EMIS,V2EMIS,DVEMIS,NLIMEM,ZEMIS(NMAXCO) 
!     ----------------------------------------------------------------  
!                                                                       
!     Read header information                                           
!                                                                       
      READ (ICOEF,900) V1EMIS,V2EMIS,DVEMIS,NLIMEM 
                                                                        
      if (nlimem.gt.nmaxco) then 
         print *, '*********************************************' 
         print *, ' Number of points on EMISSIVITY file > nmaxco' 
         print *, ' Also, check the REFLECTIVITY file' 
         stop 
      endif 
                                                                        
                                                                        
!                                                                       
!     Read in emissivity values                                         
!                                                                       
      DO 100 NGNU = 1,NLIMEM 
         READ (ICOEF,*) ZEMIS(NGNU) 
  100 END DO 
!                                                                       
      RETURN 
!                                                                       
!     FORMAT statements                                                 
!                                                                       
  900 FORMAT (3E10.3,5X,I5) 
  910 FORMAT (E15.7) 
!                                                                       
end subroutine READEM
!     -------------------------------------------------------------     
!                                                                       
      SUBROUTINE READRF(ICOEF) 
!                                                                       
!     Reads in reflection function values directly from file "REFLECTIVI
!                                                                       
      USE lblparams, ONLY: NMAXCO 
      IMPLICIT REAL*8           (V) 
!                                                                       
!     ----------------------------------------------------------------  
!     Parameter and common blocks for direct input of reflection        
!     function values                                                   
!                                                                       
!      PARAMETER (NMAXCO=4040) 
      COMMON /RFLTIN/ V1RFLT,V2RFLT,DVRFLT,NLIMRF,ZRFLT(NMAXCO) 
!     ----------------------------------------------------------------  
!                                                                       
!     Read header information                                           
!                                                                       
      READ (ICOEF,900) V1RFLT,V2RFLT,DVRFLT,NLIMRF 
                                                                        
      if (nlimrf.gt.nmaxco) then 
         print *, '*********************************************' 
         print *, ' Number of points on REFLECTIVITY file > nmaxco' 
         stop 
      endif 
!                                                                       
!     Read in reflectivity values                                       
!                                                                       
      DO 100 NGNU = 1,NLIMRF 
         READ (ICOEF,*) ZRFLT(NGNU) 
  100 END DO 
!                                                                       
      RETURN 
!                                                                       
!     FORMAT statements                                                 
!                                                                       
  900 FORMAT (3E10.3,5X,I5) 
  910 FORMAT (E15.7) 
!                                                                       
end subroutine READRF
!-----------------------------------------------------------------------
!                                                                       
      subroutine line_exception(ind,ipr,h_sub,mol,nmol,iso,iso_max) 

      character*8 h_sub 
      dimension iso_max(*) 
                                                                        
      data  mol_max_pr_1/-99/, iso_max_pr_1/-99/ 
                                                                        
      if ((ind.eq.1 .and. mol_max_pr_1.lt.0) .or.                       &
     &    (ind.eq.2 .and. iso_max_pr_1.lt.0)) then                      
         write (*,*) 
         write (*,*) 'Line file exception encountered in', h_sub 
         write (*,*) 'This message only written for first exception',   &
     &               ' for molecule and isotope cases'                  
         write (*,*) 'Other exceptions may exist' 
                                                                        
         write (ipr,*) '****************************************' 
         write (ipr,*) 'Line file exception encountered' 
         write (ipr,*) 'This message only written for first exception' 
         write (ipr,*) 'Other exceptions may exist' 
      endif 
!                                                                       
      if (ind .eq. 1) then 
          if (mol_max_pr_1 .lt. 0) then 
             mol_max_pr_1 = 11 
             write (*,*) 
             write (*,*)   ' tape3: molecule number ', mol,             &
     &             ' greater than ', nmol,' encountered and skipped'    
             write (ipr,*) ' tape3: molecule number ', mol,             &
     &             ' greater than ', nmol,' encountered and skipped'    
               write (*,*) 
            endif 
            go to 25 
!                                                                       
         else if (ind .eq. 2) then 
            if (iso_max_pr_1 .lt. 0) then 
               iso_max_pr_1 = 11 
               write (*,*) 
               write (*,*)   ' tape3: molecule number ', mol 
               write (ipr,*) ' tape3: molecule number ', mol 
                                                                        
               write (*,*)   ' tape3: isotope number ', iso,            &
     &                       ' greater than ', iso_max(mol),            &
     &                       ' encountered and skipped'                 
               write (ipr,*) ' tape3: isotope number ', iso,            &
     &                       ' greater than ', iso_max(mol),            &
     &                       ' encountered and skipped'                 
               write (*,*) 
            endif 
            go to 25 
         endif 
!                                                                       
   25    continue 
                                                                        
         return 
end subroutine line_exception
!-----------------------------------------------------------------------
!                                                                       
      subroutine scnmrg_aj(nlayer,iup_dn) 
!                                                                       
! subroutine to call scnmrg for the layer jacobian files                
!                                                                       
! nlayer is number of layers                                            
! iup_dn is used to determine what files are needed                     
!  1 = upwelling                                                        
! -1 = downwelling                                                      
!                                                                       
      USE lblparams, ONLY: MXFSC, MXLAY, MXMOL 
      IMPLICIT REAL*8 (V) 
                                                                        
!      PARAMETER (MXFSC=600,MXLAY=MXFSC+3,MXMOL=39) 
                                                                        
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
     &              NLTEFL,LNFIL4,LNGTH4                                
                                                                        
      COMMON /IADFLG/ NSPCRT,imrgsav 
      COMMON /ADRFIL/ KODFIL,kradtot,KTEMP,KFILAD,K_REFTRA,k_rddn_sfc 
                                                                        
      CHARACTER*61 FILE1,FILE2 
      DATA FILE1 /                                                      &
     &     '                                                         '/,&
     &     FILE2 /'scnfile.TMP'/                                        
                                                                        
                                  ! i/o file unit numbers               
      data iu1,iu2,iu3/82,83,84/ 
                                                                        
      CHARACTER*55 PTHT3M,PTHODI,PTHODTU,PTHODTD 
      CHARACTER*11 PTHRDRU,PTHRDRD 
      CHARACTER*3  PTHDIR,AJID 
                            ! change if PTHDIR//PTHRDRD//AJID changes si
      CHARACTER*17 FULLPTH 
                                                                        
      COMMON /ADRPNM/ PTHT3M,PTHODI,PTHODTU,PTHODTD 
      COMMON /ADRPTH/ PTHDIR,PTHRDRU,PTHRDRD,AJID 
                                                                        
      CHARACTER*10 HFMODI,HFMODTU,HFMODTD,HFMRDR 
      COMMON /ADRFRM/ HFMODI,HFMODTU,HFMODTD,HFMRDR 
                                                                        
!------                                                                 
                                                                        
! loop over layers                                                      
      do 100 ilay=1,nlayer 
                                                                        
! loop over type                                                        
!  iup_dn = 1 => upwelling                                              
!  iup_dn =-1 => downwelling                                            
                                                                        
! construct filenames                                                   
         if (iup_dn.eq.-1) then 
            FULLPTH=PTHDIR//PTHRDRd//AJID 
            WRITE(FILE1,HFMRDR) FULLPTH,ILAY 
         else 
            FULLPTH=PTHDIR//PTHRDRu//AJID 
            WRITE(FILE1,HFMRDR) FULLPTH,ILAY 
         endif 
         write(ipr,*) ' ' 
         write(ipr,*) 'scnmrg_aj Files ->' 
         write(ipr,'(a61)') file1 
         write(ipr,'(a61)') file2 
         write(ipr,*) '<- scnmrg_aj Files' 
                                                                        
         open(iu1,file=file1,form='unformatted') 
         open(iu2,form='unformatted',status='SCRATCH') 
                                                                        
                                ! scan iu1 into iu2                     
         call scnmrg(iu1,iu2) 
                                                                        
                                    ! file deleted on close             
         close(iu1,status='delete') 
                                ! don't want this deleted yet           
         rewind(iu2) 
                                                                        
         open(iu1,file=file1,form='unformatted',status='NEW') 
                                ! dummy number not needed               
         npts=0 
                                   ! copy iu2 into iu1                  
         call copyfl(npts,iu2,iu1) 
                                ! puts -99 in last line of file         
         call endfil (iu1) 
                                                                        
                                ! file contains scanned derivative      
         close(iu1) 
                                ! file will be deleted on close         
         close(iu2) 
                                                                        
                                                                        
! done with down or upwelling                                           
                                                                        
                ! ***** ilev loop over levels *****                     
  100 continue 
                                                                        
      return 
end subroutine scnmrg_aj
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                                                       
      subroutine layer2level 
!                                                                       
! subroutine to convert layer derivatives to level derivatives          
!                                                                       
      USE lblparams, ONLY: MXFSC, MXLAY, MXZMD, MXPDIM, IM2,            &
                           MXMOL, MXTRAC, IPTS, IPTS2
      IMPLICIT REAL*8 (V) 
                                                                        
 !     PARAMETER (MXFSC=600,MXLAY=MXFSC+3,MXZMD=6000,                    &
 !    &           MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXMOL=39,MXTRAC=22)    
                                                                        
! iup_dn is used to determine what to map                               
!  1 = upwelling                                                        
! -1 = downwelling                                                      
      common /dlaydlev/ilevdx,imoldx,iup_dn,                            &
     &    dxdL(mxlay,0:mxmol),dxdU(mxlay,0:mxmol)                       
                                                                        
      COMMON /IADFLG/ NSPCRT,imrgsav 
      COMMON /ADRFIL/ KODFIL,kradtot,KTEMP,KFILAD,K_REFTRA,k_rddn_sfc 
                                                                        
      character*8      XID,       HMOLID,      YID 
      real*8               SECANT,       XALTZ 
!                                                                       
      COMMON /EMIHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
     &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &
     &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF  
      COMMON /PANL/ V1P,V2P,DVP,NLIMO 
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
     &              NLTEFL,LNFIL4,LNGTH4                                
!                                                                       
      COMMON /PATHD/ PBAR(MXLAY),TBAR(MXLAY),AMOUNT(MXMOL,MXLAY),       &
     &               WN2L(MXLAY),DVL(MXLAY),WTOTL(MXLAY),ALBL(MXLAY),   &
     &               ADBL(MXLAY),AVBL(MXLAY),H2OSL(MXLAY),IPATH(MXLAY), &
     &               ITYL(MXLAY),SECNTA(MXLAY),HT1,HT2,ALTZ(0:MXLAY),   &
     &               PZ(0:MXLAY),TZ(0:MXLAY)                            
      COMMON /DEAMT/ DENM(MXMOL,MXZMD),DENP(MXMOL,MXPDIM),DRYAIR(MXZMD) 
!                                                                       
      COMMON /BNDRY/ ZBND(MXFSC),PBND(MXFSC),TBND(MXFSC),ALORNZ(MXFSC), &
     &               ADOPP(MXFSC),AVOIGT(MXFSC)                         
!                                                                       
      DIMENSION XFILHD(2),PNLHD(2) 
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHD(1),V1P) 
                                                                        
      character*4 ht1,ht2 
!                                                                       
      LOGICAL op 
      CHARACTER*61 FILE1,FILE2,FILE3 
                                                                        
      DATA FILE1 /                                                      &
     &     '                                                         '/,&
     &     FILE2 /                                                      &
     &     '                                                         '/,&
     &     FILE3 /                                                      &
     &     '                                                         '/ 
                                                                        
                                  ! i/o file unit numbers               
      data iu1,iu2,iu3/82,83,84/ 
                                                                        
      CHARACTER*55 CDUM1,PTHODI,PTHODTU,PTHODTD 
      CHARACTER*11 PTHRDRU,PTHRDRD 
      CHARACTER*3  PTHDIR,AJID 
                            ! change if PTHDIR//PTHRDRD//AJID changes si
      CHARACTER*17 FULLPTH 
                            ! also need to change levrdru,levrdrd (below
                                                                        
      COMMON /ADRPNM/ CDUM1,PTHODI,PTHODTU,PTHODTD 
      COMMON /ADRPTH/ PTHDIR,PTHRDRU,PTHRDRD,AJID 
                                                                        
      CHARACTER*10 HFMODI,HFMODTU,HFMODTD,HFMRDR 
      COMMON /ADRFRM/ HFMODI,HFMODTU,HFMODTD,HFMRDR 
                                                                        
      character*8 hmod 
      character*4 txtlev 
      data txtlev/"LEV_"/ 
! change the following if PTHDIR//txtlev//PTHRDRD//AJID changes size    
      character*21 LEVRDRu,LEVRDRd 
      character*10 hfmrdru,hfmrdrd 
                                                                        
      dimension hmod(2) 
      dimension RDL(2410),RDU(2410),RDLEV(2410) 
      dimension DUMRD(2410) 
                                                                        
! use this block only for icflg                                         
! note: from continuum module                                           
!          ipts  = same dimension as ABSRB                              
!          ipts2 = same dimension as C                                  
!      parameter (ipts=5050,ipts2=6000) 
      common /CDERIV/ icflg,idum,v1absc,v2absc,dvabsc,nptabsc,delT_pert,&
     &    dqh2oC(ipts),dTh2oC(ipts),dUh2o                               
                                                                        
!------                                                                 
! bring in the profile information to form the layer weights            
                                                                        
      inquire (97,opened=op) 
      if (op .eqv. .false.)                                             &
     &     open(97,file='AJ_atmosphere',status='old',form='unformatted')
                                                                        
      rewind (97) 
      read   (97)  xid 
      read   (97)  LMAX,NMOL,SECNT0,(HMOD(I),I=1,2),H1,H2,ANGLE,LEN 
      read   (97)  ibmax,(pbar(l),tbar(l),l=1,ibmax-1) 
      read   (97)  (pbnd(l),tbnd(l),(denm(k,l),k=1,nmol),l=1,ibmax) 

!***********************************
!MJA 10-17-2012 Uncomment these lines to get text version of AJ_atmosphere
!      open(10,file='AJ_atmosphere.txt',status='replace',action='write')
!      write   (10,*)  xid 
!      write   (10,*)  LMAX,NMOL,SECNT0,(HMOD(I),I=1,2),H1,H2,ANGLE,LEN 
!      write  (10,*)  ibmax
!      write  (10,1111) (pbar(l),tbar(l),l=1,ibmax-1) 
!1111  format (2F10.4) 
!      write   (10,1112)  (pbnd(l),tbnd(l),(denm(k,l),k=1,nmol),l=1,ibmax) 
!1112  format (2F10.4,6E15.5) 
!***********************************  

                                  
!-----------------------------------------------------------            
! compute layer-to-level conversion for analytical jacobians            
! pbar,tbar                                                             
! only go into this if imoldx was set in lblrtm                         
!                                                                       
! note that the dxdL and dxdU arrays are indexed by mol-id              
! number with the "0" index reserved for temperature                    
!                                                                       
      ilevdx=ibmax-1 
      imoldx=nmol 
                                                                        
      do 500 l=1,ilevdx 
                                                                        
         rhoU=pbnd(l+1)/(tbnd(l+1)*1.3806503E-19) 
         rhoL=pbnd(l)/(tbnd(l)*1.3806503E-19) 
         alpha=rhoU/rhoL 
         alphaT=-(tbnd(l+1)-tbnd(l))/alog(alpha) 
                                                                        
! temperature                                                           
         dxdL(l,0)=((tbar(l)-alphaT)/tbnd(l))                           &
     &        *(rhoL/(rhoL-rhoU))                                       &
     &        +(1.0-alphaT/tbnd(l))/alog(alpha)                         
                                                                        
         dxdU(l,0)=((tbar(l)-alphaT)/tbnd(l+1))                         &
     &        *(-rhoU/(rhoL-rhoU))                                      &
     &        -(1.0-alphaT/tbnd(l+1))/alog(alpha)                       
                                                                        
! molecules                                                             
         do  k=1,nmol 
                                                                        
            if (denm(k,l).ne.0.0) then 
                                                                        
               ratU=denm(k,l+1)/rhoU 
               ratL=denm(k,l)/rhoL 
                                                                        
               dxdL(l,k)=(ratL/(ratL-alpha*ratU))                       &
     &              +1.0/alog(alpha*ratU/ratL)                          
                                                                        
               dxdU(l,k)=((-alpha*ratU)/(ratL-alpha*ratU))              &
     &              -1.0/alog(alpha*ratU/ratL)                          
            else 
               dxdL(l,k)=0.0 
               dxdU(l,k)=0.0 
                                                                        
                                                                        
! check to be sure molecular amount non-zero for molecular jacobian     
               if (k.eq.nspcrt) then 
                  write(*,*) ' --- FATAL ERROR ---' 
                  write(*,*) 'molecular amount for species ',k 
                  write(*,*) '     must be non-zero ' 
                  write(*,*) 'for analytic jacobian #',nspcrt 
                  write(*,*) ' -------------------' 
                  STOP 
               endif 
                                                                        
            endif 
                                                                        
         enddo 
                                                                        
  500 continue 
                                                                        
!------                                                                 
! now obtain the AJs for levels using layer weights just calculated     
                                                                        
! initialize arrays                                                     
                                                                        
      do i=1,2410 
         rdl(i)  = -999. 
         rdu(i)  = -999. 
         dumrd(i)= -999. 
      enddo 
                                                                        
! need to close last layer derivative files that were used              
      close(KFILAD) 
                                                                        
! form level derivative file prefix and get format statement            
                                             ! upwelling                
      levrdru=pthdir//txtlev//pthrdrU//ajid 
      call qntify(levrdru,hfmrdru) 
                                             ! downwelling              
      levrdrd=pthdir//txtlev//pthrdrD//ajid 
      call qntify(levrdrd,hfmrdrd) 
                                                                        
! set derivative flag based on nspcrt                                   
      ideriv=nspcrt 
                                ! temperature derivatives               
      if (nspcrt.eq.0) ideriv=0 
                                                                        
! loop over levels                                                      
      do 100 ilev=1,ilevdx+1 
                                                                        
                                ! upper layer relative to level         
         ilay_p =ilev 
                                ! lower layer relative to level         
         ilay_m =ilev-1 
                                                                        
! loop over type                                                        
!  iup_dn = 1 => upwelling                                              
!  iup_dn =-1 => downwelling                                            
                                                                        
         wt_p = dxdL(ilay_p,ideriv) 
         wt_m = dxdU(ilay_m,ideriv) 
                                                                        
! construct filenames                                                   
         if (iup_dn.eq.-1) then 
            fullpth=pthdir//pthrdrd//ajid 
            WRITE(FILE1,HFMRDR) fullpth,ilay_m 
            fullpth=pthdir//pthrdrd//ajid 
            WRITE(FILE2,HFMRDR) fullpth,ilay_p 
            WRITE(FILE3,HFMRDRd) LEVRDRd,ilev 
         else 
            fullpth=pthdir//pthrdru//ajid 
            WRITE(FILE1,HFMRDR) fullpth,ilay_m 
            fullpth=pthdir//pthrdru//ajid 
            WRITE(FILE2,HFMRDR) fullpth,ilay_p 
            WRITE(FILE3,HFMRDRu) LEVRDRu,ilev 
         endif 
                                                                        
!         write(  *,*) ' '                                              
!         write(  *,*) 'Layer to Level Conversion Files ->'             
!         write(  *,910) file1,wt_m                                     
!         write(  *,910) file2,wt_p                                     
!         write(  *,'(a61)') file3                                      
!         write(  *,*) '<- Layer to Level Conversion Files'             
                                                                        
!         write(ipr,*) ' '                                              
!         write(ipr,*) 'Layer to Level Conversion Files ->'             
!         write(ipr,910) file1,wt_m                                     
!         write(ipr,910) file2,wt_p                                     
!         write(ipr,'(a61)') file3                                      
!         write(ipr,*) '<- Layer to Level Conversion Files'             
                                                                        
  910    format(a61,'  weight:',f10.6) 
                                                                        
! open files and read/write header                                      
!   valid layers:  1 -> ilevdx                                          
!   upper/lower boundary layers: (ilay_p=1) and (ilay_p=ilevdx)         
!   each level: layer below = (ilay_p-1), layer above = (ilay_p)        
                                                                        
         if (ilay_m.ge.1) then 
            open(unit=iu1,file=file1,                                   &
     &           form='unformatted',status='old')                       
            call bufin (iu1,keof,xfilhd(1),nfhdrf) 
         endif 
                                                                        
         if (ilay_p.le.ilevdx) then 
            open(unit=iu2,file=file2,                                   &
     &           form='unformatted',status='old')                       
            call bufin (iu2,keof,xfilhd(1),nfhdrf) 
         endif 
                                                                        
         open(unit=iu3,file=file3,                                      &
     &        form='unformatted',status='unknown')                      
         call bufout (iu3,xfilhd(1),nfhdrf) 
                                                                        
! begin loop over panels                                                
   10    continue 
                                                                        
         if (ilay_m.ge.1) then 
            CALL BUFIN (iu1,KEOF,PNLHD(1),NPHDRF) 
            IF (KEOF.LE.0) GO TO 20 
            CALL BUFIN (iu1,KEOF,RDL(1),NLIMO) 
            IF (IMRGSAV.LT.42) CALL BUFIN(iu1,KEOF,DUMRD(1),NLIMO) 
         endif 
                                                                        
         if (ilay_p.le.ilevdx) then 
            CALL BUFIN (iu2,KEOF,PNLHD(1),NPHDRF) 
            IF (KEOF.LE.0) GO TO 20 
            CALL BUFIN (iu2,KEOF,RDU(1),NLIMO) 
            IF (IMRGSAV.LT.42) CALL BUFIN(iu2,KEOF,DUMRD(1),NLIMO) 
         endif 
                                                                        
! if not upper or lower level, combine rdl and rdu                      
! if lower level use rdu                                                
! if upper level use rdl                                                
! dxdL is conversion of layer to lower level of layer                   
! dxdU is conversion of layer to upper level of layer                   
                                                                        
         if ((ilay_m.ge.1).and.(ilay_p.le.ilevdx)) then 
            do j=1,nlimO 
               rdlev(j) = rdl(j)*wt_m                                   &
     &              + rdu(j)*wt_p                                       
            enddo 
         else 
            if (ilay_m.eq.0) then 
               do j=1,nlimO 
                  rdlev(j)=rdu(j)*wt_p 
               enddo 
            else 
               do j=1,nlimO 
                  rdlev(j)=rdl(j)*wt_m 
               enddo 
            endif 
         endif 
                                                                        
! write panel and level derivative                                      
         call bufout (iu3,pnlhd(1),nphdrf) 
         call bufout (iu3,rdlev(1),nlimo) 
         IF (IMRGSAV.LT.42) call bufout(iu3,dumrd(1),nlimo) 
                                                                        
! end loop over panels                                                  
         goto 10 
                                                                        
! done with this set of panels                                          
   20    continue 
                                                                        
! close files                                                           
                                ! puts -99 in last line of file         
         call endfil (iu3) 
         close(iu1) 
         close(iu2) 
         close(iu3) 
                                                                        
! done with down or upwelling case                                      
                                                                        
                                ! ***** ilev loop over levels *****     
  100 continue 
                                                                        
      return 
end subroutine layer2level
!-----------------------------------------------------------------------
      subroutine sfcderiv(k_rddn_sfc,tbound) 
!                                                                       
! subroutine to compute surface derivatives                             
!                                                                       
      USE phys_consts, ONLY: radcn2
   USE lblparams, ONLY : dbg
      IMPLICIT REAL*8 (V) 
      character*8      XID,       HMOLID,      YID 
      real*8               SECANT,       XALTZ 
                                                                        
      COMMON /DUMHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),     &
     &                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2,TBOUNDx, &
     &                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF  
      COMMON /DUMPAN/ V1P,V2P,DVP,NLIM 
      COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
     &              NLNGTH,KDUMY,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
     &              NLTEFL,LNFIL4,LNGTH4                                
!                                                                       
      DIMENSION XFILHD(2),PNLHD(2) 
      EQUIVALENCE (XFILHD(1),XID(1)) , (PNLHD(1),V1P) 
                                                                        
      dimension RADDWN(2410),TRADWN(2410),                              &
     &     DERVOUTt(2410),DERVOUTe(2410),DERVOUTr(2410),DERVOUTe_r(2410)
                                                                        
      data iut,iue,iur,iue_r/85,86,87,88/ 
                                                                        
      character*24 tsffil,emifil,rflfil,e_rfil,                         &
     &     filoutt,filoute,filoutr,filoute_r                            
      data tsffil,emifil,rflfil,e_rfil                                  &
     &    /'AJ/LEV_RDderivTSF_-1_000',                                  &
     &     'AJ/LEV_RDderivEMI_-1_000',                                  &
     &     'AJ/LEV_RDderivRFL_-1_000',                                  &
     &     'AJ/LEV_RDderivE-R_-1_000'/                                  
                                                                        
                                                                        
!--------------------------------------------------------------------   
      filoutt   = tsffil 
      filoute   = emifil 
      filoutr   = rflfil 
      filoute_r = e_rfil 
                                                                        
      write(ipr,*) ' ' 
      write(ipr,*) 'Surface Property Derivative Output:' 
      write(ipr,'(a24)') filoutt 
      write(ipr,'(a24)') filoute 
      write(ipr,'(a24)') filoutr 
      write(ipr,'(a24)') filoute_r 
      write(ipr,*) ' ' 
                                                                        
! set up files                                                          
      rewind(k_rddn_sfc) 
                                                                        
      open(unit=iut,file=filoutt,form='unformatted',status='unknown') 
                                                                        
      open(unit=iue,file=filoute,form='unformatted',                    &
     &     status='unknown')                                            
                                                                        
      open(unit=iur,file=filoutr,form='unformatted',                    &
     &     status='unknown')                                            
                                                                        
      open(unit=iue_r,file=filoute_r,form='unformatted',                &
     &     status='unknown')                                            
                                                                        
! read header info and write to output file                             
!  reset tboundx to tbound for output file                              
      call bufin (k_rddn_sfc,keof,xfilhd(1),nfhdrf) 
      tboundx=tbound 
      call bufout (iut,  xfilhd(1),nfhdrf) 
      call bufout (iue,  xfilhd(1),nfhdrf) 
      call bufout (iur,  xfilhd(1),nfhdrf) 
      call bufout (iue_r,xfilhd(1),nfhdrf) 
                                                                        
! begin loop over panels                                                
   10 continue 
                                                                        
      CALL BUFIN (k_rddn_sfc,KEOF,PNLHD(1),NPHDRF) 
      IF (KEOF.LE.0) GO TO 20 
      CALL BUFIN (k_rddn_sfc,KEOF,RADDWN(1),NLIM) 
      CALL BUFIN (k_rddn_sfc,KEOF,TRADWN(1),NLIM) 
                                                                        
!     surface temperature, emissivity and reflectivity derivatives      
                                                                        
      XKTBND = TBOUND/RADCN2 
      VI = V1P-DVP 
   EMLAST = -1.
      VIDVEM   = VI 
      VIDVRF   = VI 
      VIDVBD   = VI 
      BBdum    = 0. 
      BBlast   = -1. 
      VIDD     = VI 
      BBdTdum  = 0. 
      BBdTlast = -1. 
                                                                        
      vidd = vi 
      VDdel  = VI 
                                                                        
      NLIM1 = 0 
      NLIM2 = 0 
      EMDUM = 0. 
   if (dbg(25)) then
      print *, 'sfcderiv:: :: CHECKED'
      dbg(25) = .false.
   endif
      EMISIV = EMISFN   (VI,DVP,VIDVEM,EMDEL,EMDUM) 
!                                                                       
   VI = V1P
   40 NLIM1 = NLIM2+1 
!                                                                       
      EMISIV = EMISFN (VI,DVP,VIDV, EMDEL,EMLAST) 
!                                                                       
      IF (VIDV.GE.9.E+4) THEN 
         NLIM2 = NLIM+1 
      ELSE 
         NLIM2 = (VIDV-V1P)/DVP+1.001 
      ENDIF 
      NLIM2 = MIN(NLIM2,NLIM) 
!                                                                       
      DO J = NLIM1, NLIM2
         BB   = PLANCK   (VI,XKTBND)
         BBdT = PLANCK_DT(VI,XKTBND, BB)
         DERVOUTt(J) = EMISIV*TRADWN(J)*BBdT 
         DERVOUTe(J) = BB*TRADWN(J) 
         dervoutr(j) = raddwn(j)*tradwn(j) 
         dervoute_r(j) = (bb-raddwn(j))*tradwn(j) 
!                                                                       
!     Increment interpolation values                                    
!                                                                       
         EMISIV = EMISIV+EMDEL 
         VI = VI+ DVP
      END DO
!                                                                       
      IF (NLIM2.LT.NLIM) GO TO 40 
                                                                        
!     write panel, derivative and transmission                          
!     (same format as RDderiv* files)                                   
      CALL BUFOUT (iut,PNLHD(1),NPHDRF) 
      CALL BUFOUT (iut,DERVOUTt(1),NLIM) 
      CALL BUFOUT (iut,TRADWN(1),NLIM) 
                                                                        
                                                                        
      CALL BUFOUT (iue,PNLHD(1),NPHDRF) 
      CALL BUFOUT (iue,DERVOUTe(1),NLIM) 
      CALL BUFOUT (iue,TRADWN(1),NLIM) 
                                                                        
      CALL BUFOUT (iur,PNLHD(1),NPHDRF) 
      CALL BUFOUT (iur,DERVOUTr(1),NLIM) 
      CALL BUFOUT (iur,TRADWN(1),NLIM) 
                                                                        
      CALL BUFOUT (iue_r,PNLHD(1),NPHDRF) 
      CALL BUFOUT (iue_r,DERVOUTe_r(1),NLIM) 
      CALL BUFOUT (iue_r,TRADWN(1),NLIM) 
                                                                        
! end loop over panels                                                  
      goto 10 
                                                                        
! done with all panels                                                  
   20 continue 
                                                                        
                                ! puts -99 in last line of file         
      call endfil (iut) 
                                ! puts -99 in last line of file         
      call endfil (iue) 
                                ! puts -99 in last line of file         
      call endfil (iur) 
                                ! puts -99 in last line of file         
      call endfil (iue_r) 
                                                                        
! close/rewind files                                                    
      close(iut) 
      close(iue) 
      close(iur) 
      close(iue_r) 
      rewind(k_rddn_sfc) 
                                                                        
      return 
                                                                        
end subroutine sfcderiv
