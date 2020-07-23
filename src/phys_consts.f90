!     path:      $HeadURL$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!
!  --------------------------------------------------------------------------
! |  Copyright ©, Atmospheric and Environmental Research, Inc., 2015         |
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
MODULE phys_consts   ! Physical constants

   implicit none
!
!    Constants from NIST May 2010
!
!     units are generally cgs
!
   real, parameter  :: PI = 3.1415926535898
   real, parameter  :: PLANCK = 6.62606876E-27
   real, parameter  :: BOLTZ = 1.3806503E-16
   real, parameter  :: CLIGHT = 2.99792458E+10
   real, parameter  :: AVOGAD = 6.02214199E+23
   real, parameter  :: ALOSMT = 2.6867775E+19
   real, parameter  :: GASCON = 8.314472E+07
   real, parameter  :: RADCN1 = 1.191042722E-12
   real, parameter  :: RADCN2 = 1.4387752
!
!     The first and second radiation constants are taken from NIST.
!     They were previously obtained from the relations:
!                            RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07
!                            RADCN2 = PLANCK*CLIGHT/BOLTZ

end module phys_consts
