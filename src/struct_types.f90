!     path:      $HeadURL$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!
!  --------------------------------------------------------------------------
! |  Copyright ©, Atmospheric and Environmental Research, Inc., 2012         |
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
MODULE STRUCT_TYPES

  TYPE :: INPUT_HEADER
     SEQUENCE
     REAL(8)    :: VMIN, VMAX
     INTEGER(4) :: NREC, NWDS
  END TYPE INPUT_HEADER

  TYPE :: INPUT_BLOCK
     SEQUENCE
     REAL(8), DIMENSION(250)     :: VNU
     REAL(4), DIMENSION(250)     :: SP, ALFA, EPP
     INTEGER(4), DIMENSION(250)  :: MOL
     REAL(4), DIMENSION(250)     :: HWHM, TMPALF, PSHIFT
     INTEGER(4), DIMENSION(250)  :: IFLG
  END TYPE INPUT_BLOCK

  TYPE :: LINE_DATA
     REAL(8), DIMENSION(250)   :: VNU
     REAL, DIMENSION(250)     :: SP, ALFA, EPP
     INTEGER, DIMENSION(250)  :: MOL
     REAL, DIMENSION(250)     :: HWHM, TMPALF, PSHIFT
     INTEGER, DIMENSION(250)  :: IFLG
     REAL, DIMENSION(250)     :: SPPSP, RECALF, ZETAI
     INTEGER, DIMENSION(250)  :: IZETA
  END TYPE LINE_DATA

  TYPE :: LINE_SHRINK
     REAL(8), DIMENSION(1250)  :: VNU
     REAL, DIMENSION(1250)    :: SP, ALFA, EPP
     INTEGER, DIMENSION(1250)  :: MOL
     REAL, DIMENSION(1250)    :: SPP, SRAD
  END TYPE LINE_SHRINK

END MODULE STRUCT_TYPES
