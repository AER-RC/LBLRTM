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
MODULE STRUCT_TYPES

   PARAMETER (MXBRDMOL=7,MXBRDMOL3=MXBRDMOL*3,NLINEREC=250)
   TYPE :: INPUT_HEADER
      SEQUENCE
      REAL(8)    :: VMIN, VMAX
      INTEGER(4) :: NREC, NWDS
   END TYPE INPUT_HEADER

   TYPE :: INPUT_BLOCK
      SEQUENCE
      REAL(8), DIMENSION(NLINEREC)     :: VNU
      REAL(4), DIMENSION(NLINEREC)     :: SP, ALFA, EPP
      INTEGER(4), DIMENSION(NLINEREC)  :: MOL
      REAL(4), DIMENSION(NLINEREC)     :: HWHM, TMPALF, PSHIFT
      INTEGER(4), DIMENSION(NLINEREC)  :: IFLG
      INTEGER(4), DIMENSION(MXBRDMOL,NLINEREC):: BRD_MOL_FLG_IN
      REAL(4), DIMENSION(MXBRDMOL3,NLINEREC) :: BRD_MOL_DAT
      REAL(4), DIMENSION(NLINEREC)     :: SPEED_DEP
   END TYPE INPUT_BLOCK

   TYPE :: LINE_DATA
      REAL(8), DIMENSION(NLINEREC)   :: VNU
      REAL, DIMENSION(NLINEREC)     :: SP, ALFA, EPP
      INTEGER, DIMENSION(NLINEREC)  :: MOL
      REAL, DIMENSION(NLINEREC)     :: HWHM, TMPALF, PSHIFT
      INTEGER, DIMENSION(NLINEREC)  :: IFLG
      REAL, DIMENSION(NLINEREC)     :: SPPSP, RECALF, ZETAI
      INTEGER, DIMENSION(NLINEREC)  :: IZETA
      INTEGER, DIMENSION(MXBRDMOL,NLINEREC)  :: BRD_MOL_FLG
      REAL, DIMENSION(MXBRDMOL,NLINEREC)     :: BRD_MOL_HW, BRD_MOL_TMP, BRD_MOL_SHFT
      REAL, DIMENSION(NLINEREC)     :: SPEED_DEP

   END TYPE LINE_DATA

   TYPE :: LINE_SHRINK
      REAL(8), DIMENSION(1250)  :: VNU
      REAL, DIMENSION(1250)    :: SP, ALFA, EPP
      INTEGER, DIMENSION(1250)  :: MOL
      REAL, DIMENSION(1250)    :: SPP, SRAD
   END TYPE LINE_SHRINK
   LOGICAL, PARAMETER :: USESHRINK = .TRUE.
end module STRUCT_TYPES
