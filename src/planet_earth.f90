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
MODULE planet_consts   ! Physical constants for Earth

   implicit none

   real, parameter :: AIRMWT = 28.964  ! air molecular weight (grams/mole)
   real, parameter :: XMASS_DRY = AIRMWT*1.E-3   ! previously was 0.0289654

CONTAINS

   FUNCTION GRAV_CONST(LATITUDE)

      USE phys_consts, ONLY: pi
      REAL, INTENT(IN), OPTIONAL  :: LATITUDE      ! in degrees
      REAL                        :: GRAV_CONST    ! in meters/s^2
      REAL                        :: REF_LAT
      REAL, PARAMETER             :: DEFAULT_LAT= 45.   ! in degrees

!         Latitude for which gravitational constant desired
      IF (PRESENT (LATITUDE) ) THEN
         REF_LAT = LATITUDE
      ELSE
         REF_LAT = DEFAULT_LAT
      END IF
!         Gravitational constant for Earth in meters/s^2
      GRAV_CONST = 9.80665 - 0.02586*COS(2.0*PI*REF_LAT/180.0)

   end function GRAV_CONST

end module planet_consts
