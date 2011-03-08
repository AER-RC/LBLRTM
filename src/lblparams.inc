C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$
C
C  --------------------------------------------------------------------------
C |  Copyright ©, Atmospheric and Environmental Research, Inc., 2011         |
C |                                                                          |
C |  All rights reserved. This source code is part of the LBLRTM software    |
C |  and is designed for scientific and research purposes. Atmospheric and   |
C |  Environmental Research, Inc. (AER) grants USER the right to download,   |
C |  install, use and copy this software for scientific and research         |
C |  purposes only. This software may be redistributed as long as this       |
C |  copyright notice is reproduced on any copy made and appropriate         |
C |  acknowledgment is given to AER. This software or any modified version   |
C |  of this software may not be incorporated into proprietary software or   |
C |  commercial software offered for sale.                                   |
C |                                                                          |
C |  This software is provided as is without any express or implied          |
C |  warranties.                                                             |
C |                       (http://www.rtweb.aer.com/)                        |
C  --------------------------------------------------------------------------
C
      PARAMETER (MXMOL=39,MXSPC=5,Max_ISO=20)
      PARAMETER (MXFSC=600, MXLAY=MXFSC+3,MX_XS=38)
      PARAMETER (MXZMD=6000,MXPDIM=MXLAY+MXZMD,IM2=MXPDIM-2,MXTRAC=22)
      PARAMETER (NFPTS=2001,NFMX=1.3*NFPTS)
      PARAMETER (NMAXCO=4040)
      PARAMETER (IPTS=5050,IPTS2=6000)
      PARAMETER (N_ABSRB=5050)
      PARAMETER (NT=119)
      PARAMETER (NN_TBL=10000, NDIM=2410, ND2=5000)
      PARAMETER (MAXSTATE=26)
