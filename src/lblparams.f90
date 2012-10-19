!     path:      $HeadURL$
!     author:    $Author$
!     revision:  $Revision: 10658 $
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
MODULE lblparams   ! Parameters for array dimensions in lblrtm

  implicit none

      integer, parameter :: MXMOL=39, MXSPC=5, Max_ISO=20
      integer, parameter :: MXFSC=600, MXLAY=MXFSC+3, MX_XS=38
      integer, parameter :: MXZMD=6000, MXPDIM=MXLAY+MXZMD
      integer, parameter :: IM2=MXPDIM-2, MXTRAC=22
      integer, parameter :: NFPTS=2001, NFMX=1.3*NFPTS
      integer, parameter :: NMAXCO=4040, NUMZ = 101
      integer, parameter :: IPTS=5050, IPTS2=6000
      integer, parameter :: N_ABSRB=5050, nzeta=101
      integer, parameter :: NT=119, Nmax=600
      integer, parameter :: NN_TBL=10000, NDIM=2410, ND2=5000
      integer, parameter :: MAXSTATE=26
      integer, parameter :: NFLTPT=3001

END MODULE lblparams
