!     path:      $HeadURL$
!     author:    $Author$
!     revision:  $Revision: 10658 $
!     created:   $Date$
!
!  --------------------------------------------------------------------------
! |  Copyright �, Atmospheric and Environmental Research, Inc., 2012         |
! |                                                                          |
! |  All rights reserved. This source code is part of the LBLRTM software    |
! |  and is designed for scientific and research purposes. Atmospheric and   |
! |  Environmental Research, Inc. (AER) grants USER the right to download,   |
! |  install, use and copy this software for scientific and research         |
! |  purposes only. This software may be redistributed as long as this       |
! |  copyright notice is reproduced on any copy made and appropriate         |
! |  acknowledgment is given to AER. This software or any modified version   |
! |  of this software may not be incorporated into proprietary software or   |
! |  commercial software offered for sale.                                   |
! |                                                                          |
! |  This software is provided as is without any express or implied          |
! |  warranties.                                                             |
! |                       (http://www.rtweb.aer.com/)                        |
!  --------------------------------------------------------------------------
!
MODULE lblparams   ! Parameters for array dimensions in lblrtm

   implicit none

   integer, parameter :: MXMOL=47, MXSPC=5, Max_ISO=20, MXISOTPL=10
!
   integer, parameter :: MXFSC=600, MXLAY=MXFSC+3, MX_XS=40
   integer, parameter :: MXZMD=7000, MXPDIM=MXLAY+MXZMD
   integer, parameter :: IM2=MXPDIM-2, MXTRAC=41
   integer, parameter :: NFPTS=2001, NFMX=1.3*NFPTS
   integer, parameter :: NMAXCO=4040, NUMZ = 101
   integer, parameter :: IPTS=5050, IPTS2=6000
   integer, parameter :: N_ABSRB=5050, nzeta=101
   integer, parameter :: NT=119, Nmax=600
   integer, parameter :: NN_TBL=10000, NDIM=2410, ND2=5000
   integer, parameter :: MAXSTATE=26
   integer, parameter :: NFLTPT=3001
   
!    parameter that control usage output 
   logical ::  dbg(100) = .FALSE.

!    transfered from  create_fn_tbls of xmerge.f90
!    data od_switch /0.06/
  real, parameter      :: od_lo = 0.06


! HITRAN 2012 isotopologue information
   real :: isotpl_num(mxmol)   ! Number of isotopologues for each species
   real :: isotpl_abd(mxmol,mxisotpl) ! Isotopologue abundance by species
   integer :: isotpl_code(mxmol,mxisotpl)  ! Isotopologue code by species

! HITRAN 2012 isotopologue abundance values are based on:
! De Bievre, P., N.E. Holden, and I.L. Barnes, Isotopic Abundances and
! Atomic Weights of the Elements, J. Phys. Chem. Ref. Data, 13,
! 809-891, 1984.

! Molecule names
!      DATA CMOL   /                                                     &
!     &     '  H2O ','  CO2 ','   O3 ','  N2O ','   CO ','  CH4 ',       &
!     &     '   O2 ','   NO ','  SO2 ','  NO2 ','  NH3 ',' HNO3 ',       &
!     &     '   OH ','   HF ','  HCL ','  HBR ','   HI ','  CLO ',       &
!     &     '  OCS ',' H2CO ',' HOCL ','   N2 ','  HCN ','CH3CL ',       &
!     &     ' H2O2 ',' C2H2 ',' C2H6 ','  PH3 ',' COF2 ','  SF6 ',       &
!     &     '  H2S ','HCOOH ','  HO2 ','    O ','ClONO2','  NO+ ',       &
!     &     ' HOBr ',' C2H4 ','CH3OH ',' CH3Br',' CH3CN','  CF4 ',       &
!     &     ' C4H2 ',' HC3N ','   H2 ','   CS ','  SO3 '/

! Number of isotopologues for each species
   data isotpl_num(:) /                                              &
   &           6,      10,       5,       5,       6,       4,        &
   &           3,       3,       2,       1,       2,       2,        &
   &           3,       2,       4,       4,       2,       2,        &
   &           5,       3,       2,       2,       3,       2,        &
   &           1,       3,       2,       1,       2,       1,        &
   &           3,       1,       1,       1,       2,       1,        &
   &           2,       2,       1,       2,       1,       1,        &
   &           1,       1,       2,       4,       1/

! HITRAN 2012/AFGL isotopologue codes for each species
! H2O
   data isotpl_code(1,:) /                                           &
   &         161,       181,      171,      162,      182,            &
   &         172,         0,        0,        0,        0/
! CO2
   data isotpl_code(2,:) /                                           &
   &         626,       636,      628,      627,      638,            &
   &         637,       828,      827,      727,      838/
! O3
   data isotpl_code(3,:) /                                           &
   &         666,       668,      686,      667,      676,            &
   &           0,         0,        0,        0,        0/
! N2O
   data isotpl_code(4,:) /                                           &
   &         446,       456,      546,      448,      447,            &
   &           0,         0,        0,        0,        0/
! CO
   data isotpl_code(5,:) /                                           &
   &          26,        36,       28,       27,       38,            &
   &          37,         0,        0,        0,        0/
! CH4
   data isotpl_code(6,:) /                                           &
   &         211,       311,      212,      312,        0,            &
   &           0,         0,        0,        0,        0/
! O2
   data isotpl_code(7,:) /                                           &
   &          66,        68,       67,        0,        0,            &
   &           0,         0,        0,        0,        0/
! NO
   data isotpl_code(8,:) /                                           &
   &          46,        56,       48,        0,        0,            &
   &           0,         0,        0,        0,        0/
! SO2
   data isotpl_code(9,:) /                                           &
   &         626,       646,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! NO2
   data isotpl_code(10,:) /                                          &
   &         646,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! NH3
   data isotpl_code(11,:) /                                          &
   &        4111,      5111,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! HNO3
   data isotpl_code(12,:) /                                          &
   &         146,       156,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! OH
   data isotpl_code(13,:) /                                          &
   &          61,        81,       62,        0,        0,            &
   &           0,         0,        0,        0,        0/
! HF
   data isotpl_code(14,:) /                                          &
   &          19,        29,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! HCl
   data isotpl_code(15,:) /                                          &
   &          15,        17,       25,       27,        0,            &
   &           0,         0,        0,        0,        0/
! HBr
   data isotpl_code(16,:) /                                          &
   &          19,        11,       29,       21,        0,            &
   &           0,         0,        0,        0,        0/
! HI
   data isotpl_code(17,:) /                                          &
   &          17,        27,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! ClO
   data isotpl_code(18,:) /                                          &
   &          56,        76,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! OCS
   data isotpl_code(19,:) /                                          &
   &         622,       624,      632,      623,      822,            &
   &           0,         0,        0,        0,        0/
! H2CO
   data isotpl_code(20,:) /                                          &
   &         126,       136,      128,        0,        0,            &
   &           0,         0,        0,        0,        0/
! HOCl
   data isotpl_code(21,:) /                                          &
   &         165,       167,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! N2
   data isotpl_code(22,:) /                                          &
   &          44,        45,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! HCN
   data isotpl_code(23,:) /                                          &
   &         124,       134,      125,        0,        0,            &
   &           0,         0,        0,        0,        0/
! CH3Cl
   data isotpl_code(24,:) /                                          &
   &         215,       217,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! H2O2
   data isotpl_code(25,:) /                                          &
   &        1661,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! C2H2
   data isotpl_code(26,:) /                                          &
   &        1221,      1231,     1222,        0,        0,            &
   &           0,         0,        0,        0,        0/
! C2H6
   data isotpl_code(27,:) /                                          &
   &        1221,      1231,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! PH3
   data isotpl_code(28,:) /                                          &
   &        1111,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! COF2
   data isotpl_code(29,:) /                                          &
   &         269,       369,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! SF6
   data isotpl_code(30,:) /                                          &
   &          29,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! H2S
   data isotpl_code(31,:) /                                          &
   &         121,       141,      131,        0,        0,            &
   &           0,         0,        0,        0,        0/
! HCOOH
   data isotpl_code(32,:) /                                          &
   &         126,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! HO2
   data isotpl_code(33,:) /                                          &
   &         166,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! O
   data isotpl_code(34,:) /                                          &
   &           6,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! ClONO2
   data isotpl_code(35,:) /                                          &
   &        5646,      7646,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! NO+
   data isotpl_code(36,:) /                                          &
   &          46,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! HOBr
   data isotpl_code(37,:) /                                          &
   &         169,       161,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! C2H4
   data isotpl_code(38,:) /                                          &
   &         221,       231,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! CH3OH
   data isotpl_code(39,:) /                                          &
   &        2161,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! CH3Br
   data isotpl_code(40,:) /                                          &
   &         219,       211,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! CH3CN
   data isotpl_code(41,:) /                                          &
   &        2124,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! CF4
   data isotpl_code(42,:) /                                          &
   &          29,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! C4H2
   data isotpl_code(43,:) /                                          &
   &        2211,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! HC3N
   data isotpl_code(44,:) /                                          &
   &        1224,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! H2
   data isotpl_code(45,:) /                                          &
   &          11,        12,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/
! CS
   data isotpl_code(46,:) /                                          &
   &          22,        24,       32,       23,        0,            &
   &           0,         0,        0,        0,        0/
! SO3
   data isotpl_code(47,:) /                                          &
   &          26,         0,        0,        0,        0,            &
   &           0,         0,        0,        0,        0/

! HITRAN 2012 isotopologue abundance for each species
! H2O
   data isotpl_abd(1,:) /                                            &
   &     9.973e-1, 1.999e-3, 3.719e-4, 3.107e-4, 6.230e-7,            &
   &     1.158e-7, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CO2
   data isotpl_abd(2,:) /                                            &
   &     9.842e-1, 1.106e-2, 3.947e-3, 7.340e-4, 4.434e-5,            &
   &     8.246e-6, 3.957e-6, 1.472e-6, 1.368e-7, 4.446e-8/
! O3
   data isotpl_abd(3,:) /                                            &
   &     9.929e-1, 3.982e-3, 1.991e-3, 7.405e-4, 3.702e-4,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! N2O
   data isotpl_abd(4,:) /                                            &
   &     9.903e-1, 3.641e-3, 3.641e-3, 1.986e-3, 3.693e-4,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CO
   data isotpl_abd(5,:) /                                            &
   &     9.865e-1, 1.108e-2, 1.978e-3, 3.679e-4, 2.223e-5,            &
   &     4.133e-6, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CH4
   data isotpl_abd(6,:) /                                            &
   &     9.883e-1, 1.110e-2, 6.158e-4, 6.918e-6, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! O2
   data isotpl_abd(7,:) /                                            &
   &     9.953e-1, 3.991e-3, 7.422e-4, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! NO
   data isotpl_abd(8,:) /                                            &
   &     9.940e-1, 3.654e-3, 1.993e-3, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! SO2
   data isotpl_abd(9,:) /                                            &
   &     9.457e-1, 4.195e-2, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! NO2
   data isotpl_abd(10,:) /                                           &
   &     9.916e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! NH3
   data isotpl_abd(11,:) /                                           &
   &     9.959e-1, 3.661e-3, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HNO3
   data isotpl_abd(12,:) /                                           &
   &     9.891e-1, 3.636e-3, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! OH
   data isotpl_abd(13,:) /                                           &
   &     9.975e-1, 2.000e-3, 1.554e-4, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HF
   data isotpl_abd(14,:) /                                           &
   &     9.998e-1, 1.557e-4, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HCl
   data isotpl_abd(15,:) /                                           &
   &     7.576e-1, 2.423e-1, 1.180e-4, 3.774e-5, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HBr
   data isotpl_abd(16,:) /                                           &
   &     5.068e-1, 4.931e-1, 7.894e-5, 7.680e-5, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HI
   data isotpl_abd(17,:) /                                           &
   &     9.998e-1, 1.557e-4, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! ClO
   data isotpl_abd(18,:) /                                           &
   &     7.559e-1, 2.417e-1, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! OCS
   data isotpl_abd(19,:) /                                           &
   &     9.374e-1, 4.158e-2, 1.053e-2, 7.399e-3, 1.880e-3,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! H2CO
   data isotpl_abd(20,:) /                                           &
   &     9.862e-1, 1.108e-2, 1.978e-3, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HOCl
   data isotpl_abd(21,:) /                                           &
   &     7.558e-1, 2.417e-1, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! N2
   data isotpl_abd(22,:) /                                           &
   &     9.927e-1, 7.478e-3, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HCN
   data isotpl_abd(23,:) /                                           &
   &     9.851e-1, 1.107e-2, 3.622e-3, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CH3Cl
   data isotpl_abd(24,:) /                                           &
   &     7.489e-1, 2.395e-1, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! H2O2
   data isotpl_abd(25,:) /                                           &
   &     9.950e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! C2H2
   data isotpl_abd(26,:) /                                           &
   &     9.776e-1, 2.197e-2, 3.046e-4, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! C2H6
   data isotpl_abd(27,:) /                                           &
   &     9.770e-1, 2.195e-2, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! PH3
   data isotpl_abd(28,:) /                                           &
   &     9.995e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! COF2
   data isotpl_abd(29,:) /                                           &
   &     9.865e-1, 1.108e-2, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! SF6
   data isotpl_abd(30,:) /                                           &
   &     9.502e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! H2S
   data isotpl_abd(31,:) /                                           &
   &     9.499e-1, 4.214e-2, 7.498e-3, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HCOOH
   data isotpl_abd(32,:) /                                           &
   &     9.839e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HO2
   data isotpl_abd(33,:) /                                           &
   &     9.951e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! O
   data isotpl_abd(34,:) /                                           &
   &     9.976e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! ClONO2
   data isotpl_abd(35,:) /                                           &
   &     7.496e-1, 2.397e-1, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! NO+
   data isotpl_abd(36,:) /                                           &
   &     9.940e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HOBr
   data isotpl_abd(37,:) /                                           &
   &     5.056e-1, 4.919e-1, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! C2H4
   data isotpl_abd(38,:) /                                           &
   &     9.773e-1, 2.196e-2, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CH3OH
   data isotpl_abd(39,:) /                                           &
   &     9.859e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CH3Br
   data isotpl_abd(40,:) /                                           &
   &     5.010e-1, 4.874e-1, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CH3CN
   data isotpl_abd(41,:) /                                           &
   &     9.739e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CF4
   data isotpl_abd(42,:) /                                           &
   &     9.889e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! C4H2
   data isotpl_abd(43,:) /                                           &
   &     9.560e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! HC3N
   data isotpl_abd(44,:) /                                           &
   &     9.633e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! H2
   data isotpl_abd(45,:) /                                           &
   &     9.997e-1, 3.114e-4, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! CS
   data isotpl_abd(46,:) /                                           &
   &     9.396e-1, 4.168e-2, 1.056e-2, 7.417e-3, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
! SO3
   data isotpl_abd(47,:) /                                           &
   &     9.434e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00,            &
   &     0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00/
!

end module lblparams
