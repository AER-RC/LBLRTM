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

   public :: decode_mol_id
   ! 
   integer, parameter :: MXMOL=47     ! max number of molecules
   ! MXISOTPL	sets the isotope dimension of ISOTPL_ABD, ISOTPL_CODE, ISOTPL_FLAG, WKI, and WKL_ISOTPL. 
   ! The line calculation checks ISOTPL_FLAG(M,ISO) to choose between the parent-molecule amount WK(M) 
   ! and the isotope-specific amount WKI(M,ISO). It checks that flag even when the user has not specified 
   ! isotope-specific profiles.
   !
   ! Max_ISO  Sets the isotope dimension of the ISO82 identification table and reserves 
   ! MAXSTATE × Max_ISO vibrational-state slots per molecule. The non-LTE TAPE4 reader checks 
   ! the input isotope against this limit. It does not set the sizes of the mass or partition-sum 
   ! correction arrays.   
   integer, parameter :: Max_ISO=20   ! max number of isotopes (used in to create arrays)
   integer, parameter :: MXISOTPL=12  ! max number of isotopes (used in to create arrays)
   integer, parameter :: MXSPC=5
!
   integer, parameter :: MXFSC=600, MXLAY=MXFSC+3, MX_XS=50
   integer, parameter :: XS_SIZE = 6*5*MX_XS

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
   real :: isotpl_abd(mxmol,mxisotpl) ! Isotopologue abundance by species
   integer :: isotpl_code(mxmol,mxisotpl)  ! Isotopologue code by species
   integer, private :: i

! HITRAN 2012 isotopologue abundance values are based on:
! De Bievre, P., N.E. Holden, and I.L. Barnes, Isotopic Abundances and
! Atomic Weights of the Elements, J. Phys. Chem. Ref. Data, 13,
! 809-891, 1984.

! Molecule names
!      DATA CMOL   /                                                     &
!     &     '  H2O ','  CO2 ','   O3 ','  N2O ','   CO ','  CH4 ',       &
!     &     '   H2O','   CO2','    O3','   N2O','    CO','   CH4',       &
!     &     '    O2','    NO','   SO2','   NO2','   NH3','  HNO3',       &
!     &     '    OH','    HF','   HCl','   HBr','    HI','   ClO',       &
!     &     '   OCS','  H2CO','  HOCl','    N2','   HCN',' CH3Cl',       &
!     &     '  H2O2','  C2H2','  C2H6','   PH3','  COF2','   SF6',       &
!     &     '   H2S',' HCOOH','   HO2','     O','ClONO2','   NO+',       &
!     &     '  HOBr','  C2H4',' CH3OH',' CH3Br',' CH3CN','   CF4',       &
!     &     '  C4H2','  HC3N','    H2','    CS','   SO3/
! Number of isotopologues for each species
   integer, parameter :: isotpl_num(mxmol) =  (/                                               &
     &            7,      12,       5,       5,       6,       4,            &
     &            3,       3,       4,       3,       2,       2,            &
     &            3,       2,       4,       4,       2,       2,            &
     &            6,       3,       2,       2,       3,       2,            &
     &            1,       3,       2,       1,       2,       1,            &
     &            3,       2,       1,       1,       2,       1,            &
     &            2,       2,       1,       2,       1,       1,            &
     &            1,       1,       2,       4,       1/)
    ! HITRAN 2024/AFGL isotopologue codes for each species

! H2O
data (isotpl_code(1,i),i=1,isotpl_num(1)) /                        &
     &          161,     181,     171,     162,     182,     172,            &
     &          262/
! CO2
data (isotpl_code(2,i),i=1,isotpl_num(2)) /                        &
     &          626,     636,     628,     627,     638,     637,            &
     &          828,     827,     727,     838,     837,     737/
! O3
data (isotpl_code(3,i),i=1,isotpl_num(3)) /                        &
     &          666,     668,     686,     667,     676/
! N2O
data (isotpl_code(4,i),i=1,isotpl_num(4)) /                        &
     &          446,     456,     546,     448,     447/
! CO
data (isotpl_code(5,i),i=1,isotpl_num(5)) /                        &
     &           26,      36,      28,      27,      38,      37/
! CH4
data (isotpl_code(6,i),i=1,isotpl_num(6)) /                        &
     &          211,     311,     212,     312/
! O2
data (isotpl_code(7,i),i=1,isotpl_num(7)) /                        &
     &           66,      68,      67/
! NO
data (isotpl_code(8,i),i=1,isotpl_num(8)) /                        &
     &           46,      56,      48/
! SO2
data (isotpl_code(9,i),i=1,isotpl_num(9)) /                        &
     &          626,     646,     636,     628/
! NO2
data (isotpl_code(10,i),i=1,isotpl_num(10)) /                        &
     &          646,     656,     648/
! NH3
data (isotpl_code(11,i),i=1,isotpl_num(11)) /                        &
     &         4111,    5111/
! HNO3
data (isotpl_code(12,i),i=1,isotpl_num(12)) /                        &
     &          146,     156/
! OH
data (isotpl_code(13,i),i=1,isotpl_num(13)) /                        &
     &           61,      81,      62/
! HF
data (isotpl_code(14,i),i=1,isotpl_num(14)) /                        &
     &           19,      29/
! HCl
data (isotpl_code(15,i),i=1,isotpl_num(15)) /                        &
     &           15,      17,      25,      27/
! HBr
data (isotpl_code(16,i),i=1,isotpl_num(16)) /                        &
     &           19,      11,      29,      21/
! HI
data (isotpl_code(17,i),i=1,isotpl_num(17)) /                        &
     &           17,      27/
! ClO
data (isotpl_code(18,i),i=1,isotpl_num(18)) /                        &
     &           56,      76/
! OCS
data (isotpl_code(19,i),i=1,isotpl_num(19)) /                        &
     &          622,     624,     632,     623,     822,     634/
! H2CO
data (isotpl_code(20,i),i=1,isotpl_num(20)) /                        &
     &          126,     136,     128/
! HOCl
data (isotpl_code(21,i),i=1,isotpl_num(21)) /                        &
     &          165,     167/
! N2
data (isotpl_code(22,i),i=1,isotpl_num(22)) /                        &
     &           44,      45/
! HCN
data (isotpl_code(23,i),i=1,isotpl_num(23)) /                        &
     &          124,     134,     125/
! CH3Cl
data (isotpl_code(24,i),i=1,isotpl_num(24)) /                        &
     &          215,     217/
! H2O2
data (isotpl_code(25,i),i=1,isotpl_num(25)) /                        &
     &         1661/
! C2H2
data (isotpl_code(26,i),i=1,isotpl_num(26)) /                        &
     &         1221,    1231,    1222/
! C2H6
data (isotpl_code(27,i),i=1,isotpl_num(27)) /                        &
     &         1221,    1231/
! PH3
data (isotpl_code(28,i),i=1,isotpl_num(28)) /                        &
     &         1111/
! COF2
data (isotpl_code(29,i),i=1,isotpl_num(29)) /                        &
     &          269,     369/
! SF6
data (isotpl_code(30,i),i=1,isotpl_num(30)) /                        &
     &           29/
! H2S
data (isotpl_code(31,i),i=1,isotpl_num(31)) /                        &
     &          121,     141,     131/
! HCOOH
data (isotpl_code(32,i),i=1,isotpl_num(32)) /                        &
     &          126,     136/
! HO2
data (isotpl_code(33,i),i=1,isotpl_num(33)) /                        &
     &          166/
! O
data (isotpl_code(34,i),i=1,isotpl_num(34)) /                        &
     &            6/
! ClONO2
data (isotpl_code(35,i),i=1,isotpl_num(35)) /                        &
     &         5646,    7646/
! NO+
data (isotpl_code(36,i),i=1,isotpl_num(36)) /                        &
     &           46/
! HOBr
data (isotpl_code(37,i),i=1,isotpl_num(37)) /                        &
     &          169,     161/
! C2H4
data (isotpl_code(38,i),i=1,isotpl_num(38)) /                        &
     &          221,     231/
! CH3OH
data (isotpl_code(39,i),i=1,isotpl_num(39)) /                        &
     &         2161/
! CH3Br
data (isotpl_code(40,i),i=1,isotpl_num(40)) /                        &
     &          219,     211/
! CH3CN
data (isotpl_code(41,i),i=1,isotpl_num(41)) /                        &
     &         2124/
! CF4
data (isotpl_code(42,i),i=1,isotpl_num(42)) /                        &
     &           29/
! C4H2
data (isotpl_code(43,i),i=1,isotpl_num(43)) /                        &
     &         2211/
! HC3N
data (isotpl_code(44,i),i=1,isotpl_num(44)) /                        &
     &         1224/
! H2
data (isotpl_code(45,i),i=1,isotpl_num(45)) /                        &
     &           11,      12/
! CS
data (isotpl_code(46,i),i=1,isotpl_num(46)) /                        &
     &           22,      24,      32,      23/
! SO3
data (isotpl_code(47,i),i=1,isotpl_num(47)) /                        &
     &           26/
! HITRAN 2024 isotopologue abundance for each species
! H2O
data (isotpl_abd(1,i),i=1,isotpl_num(1)) /                        &
     &     9.9732e-1,1.9998e-3,3.7188e-4,3.1069e-4,6.2300e-7,1.1585e-7,            &
     &     2.4197e-8/
! CO2
data (isotpl_abd(2,i),i=1,isotpl_num(2)) /                        &
     &     9.8420e-1,1.1057e-2,3.9471e-3,7.3399e-4,4.4345e-5,8.2462e-6,            &
     &     3.9573e-6,1.4718e-6,1.3685e-7,4.4460e-8,1.6535e-8,1.5374e-9/
! O3
data (isotpl_abd(3,i),i=1,isotpl_num(3)) /                        &
     &     9.9290e-1,3.9819e-3,1.9910e-3,7.4047e-4,3.7024e-4/
! N2O
data (isotpl_abd(4,i),i=1,isotpl_num(4)) /                        &
     &     9.9033e-1,3.6409e-3,3.6409e-3,1.9858e-3,3.6928e-4/
! CO
data (isotpl_abd(5,i),i=1,isotpl_num(5)) /                        &
     &     9.8654e-1,1.1084e-2,1.9782e-3,3.6787e-4,2.2225e-5,4.1329e-6/
! CH4
data (isotpl_abd(6,i),i=1,isotpl_num(6)) /                        &
     &     9.8827e-1,1.1103e-2,6.1575e-4,6.9179e-6/
! O2
data (isotpl_abd(7,i),i=1,isotpl_num(7)) /                        &
     &     9.9526e-1,3.9914e-3,7.4224e-4/
! NO
data (isotpl_abd(8,i),i=1,isotpl_num(8)) /                        &
     &     9.9397e-1,3.6543e-3,1.9931e-3/
! SO2
data (isotpl_abd(9,i),i=1,isotpl_num(9)) /                        &
     &     9.4568e-1,4.1950e-2,7.4645e-3,3.7926e-3/
! NO2
data (isotpl_abd(10,i),i=1,isotpl_num(10)) /                        &
     &     9.9162e-1,3.6456e-3,3.9768e-3/
! NH3
data (isotpl_abd(11,i),i=1,isotpl_num(11)) /                        &
     &     9.9587e-1,3.6613e-3/
! HNO3
data (isotpl_abd(12,i),i=1,isotpl_num(12)) /                        &
     &     9.8911e-1,3.6364e-3/
! OH
data (isotpl_abd(13,i),i=1,isotpl_num(13)) /                        &
     &     9.9747e-1,2.0001e-3,1.5537e-4/
! HF
data (isotpl_abd(14,i),i=1,isotpl_num(14)) /                        &
     &     9.9984e-1,1.5574e-4/
! HCl
data (isotpl_abd(15,i),i=1,isotpl_num(15)) /                        &
     &     7.5759e-1,2.4226e-1,1.1801e-4,3.7735e-5/
! HBr
data (isotpl_abd(16,i),i=1,isotpl_num(16)) /                        &
     &     5.0678e-1,4.9306e-1,7.8938e-5,7.6802e-5/
! HI
data (isotpl_abd(17,i),i=1,isotpl_num(17)) /                        &
     &     9.9984e-1,1.5574e-4/
! ClO
data (isotpl_abd(18,i),i=1,isotpl_num(18)) /                        &
     &     7.5591e-1,2.4172e-1/
! OCS
data (isotpl_abd(19,i),i=1,isotpl_num(19)) /                        &
     &     9.3739e-1,4.1583e-2,1.0531e-2,7.3991e-3,1.8797e-3,4.6718e-4/
! H2CO
data (isotpl_abd(20,i),i=1,isotpl_num(20)) /                        &
     &     9.8624e-1,1.1080e-2,1.9776e-3/
! HOCl
data (isotpl_abd(21,i),i=1,isotpl_num(21)) /                        &
     &     7.5579e-1,2.4168e-1/
! N2
data (isotpl_abd(22,i),i=1,isotpl_num(22)) /                        &
     &     9.9269e-1,7.2992e-3/
! HCN
data (isotpl_abd(23,i),i=1,isotpl_num(23)) /                        &
     &     9.8511e-1,1.1068e-2,3.6217e-3/
! CH3Cl
data (isotpl_abd(24,i),i=1,isotpl_num(24)) /                        &
     &     7.4894e-1,2.3949e-1/
! H2O2
data (isotpl_abd(25,i),i=1,isotpl_num(25)) /                        &
     &     9.9495e-1/
! C2H2
data (isotpl_abd(26,i),i=1,isotpl_num(26)) /                        &
     &     9.7760e-1,2.1966e-2,3.0455e-4/
! C2H6
data (isotpl_abd(27,i),i=1,isotpl_num(27)) /                        &
     &     9.7699e-1,2.1953e-2/
! PH3
data (isotpl_abd(28,i),i=1,isotpl_num(28)) /                        &
     &     9.9953e-1/
! COF2
data (isotpl_abd(29,i),i=1,isotpl_num(29)) /                        &
     &     9.8654e-1,1.1084e-2/
! SF6
data (isotpl_abd(30,i),i=1,isotpl_num(30)) /                        &
     &     9.5018e-1/
! H2S
data (isotpl_abd(31,i),i=1,isotpl_num(31)) /                        &
     &     9.4988e-1,4.2137e-2,7.4977e-3/
! HCOOH
data (isotpl_abd(32,i),i=1,isotpl_num(32)) /                        &
     &     9.8390e-1,1.1054e-2/
! HO2
data (isotpl_abd(33,i),i=1,isotpl_num(33)) /                        &
     &     9.9511e-1/
! O
data (isotpl_abd(34,i),i=1,isotpl_num(34)) /                        &
     &     9.9763e-1/
! ClONO2
data (isotpl_abd(35,i),i=1,isotpl_num(35)) /                        &
     &     7.4957e-1,2.3969e-1/
! NO+
data (isotpl_abd(36,i),i=1,isotpl_num(36)) /                        &
     &     9.9397e-1/
! HOBr
data (isotpl_abd(37,i),i=1,isotpl_num(37)) /                        &
     &     5.0558e-1,4.9189e-1/
! C2H4
data (isotpl_abd(38,i),i=1,isotpl_num(38)) /                        &
     &     9.7729e-1,2.1959e-2/
! CH3OH
data (isotpl_abd(39,i),i=1,isotpl_num(39)) /                        &
     &     9.8593e-1/
! CH3Br
data (isotpl_abd(40,i),i=1,isotpl_num(40)) /                        &
     &     5.0099e-1,4.8743e-1/
! CH3CN
data (isotpl_abd(41,i),i=1,isotpl_num(41)) /                        &
     &     9.7387e-1/
! CF4
data (isotpl_abd(42,i),i=1,isotpl_num(42)) /                        &
     &     9.8889e-1/
! C4H2
data (isotpl_abd(43,i),i=1,isotpl_num(43)) /                        &
     &     9.5600e-1/
! HC3N
data (isotpl_abd(44,i),i=1,isotpl_num(44)) /                        &
     &     9.6335e-1/
! H2
data (isotpl_abd(45,i),i=1,isotpl_num(45)) /                        &
     &     9.9969e-1,3.1143e-4/
! CS
data (isotpl_abd(46,i),i=1,isotpl_num(46)) /                        &
     &     9.3962e-1,4.1682e-2,1.0556e-2,7.4167e-3/
! SO3
data (isotpl_abd(47,i),i=1,isotpl_num(47)) /                        &
     &     9.4343e-1/
     !
contains
    subroutine decode_mol_id(molid, mol, iso)
        implicit none
        integer, intent(in) :: molid
        integer, intent(out):: mol
        integer, intent(out):: iso
        ! LNFL encodes molid as 
        ! MOLID = 100*isotope + molecule
        mol = MOD(molid, 100)
        iso = molid/100
    end subroutine decode_mol_id
end module lblparams
