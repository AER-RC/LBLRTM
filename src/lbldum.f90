!     path:      $HeadURL$
!     revision:  $Revision$
!     created:   $Date$
!     presently: %H%  %T%
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
!     SUBROUTINE NONLTE (MPTS)
!     PRINT * , ' NONLTE NOT IMPLEMENTED'
!     STOP
!     END

      SUBROUTINE LASER(VLAS,MFILE,JAERSL)
         PRINT * , ' LASER NOT IMPLEMENTED'
         STOP
      end subroutine LASER
      FUNCTION RANDM(IRAND)
         PRINT * , ' RANDM NOT IMPLEMENTED PROPERLY'
         RANDM=0.5
         IRAND=IABS(IRAND)
         RETURN
      end function RANDM
      SUBROUTINE PLTID3(PROGID,XSIZ,YSIZ,SCAL)
         DIMENSION PROGID(3)
         PRINT * , ' PLTID3 NOT IMPLEMENTED'
         RETURN
      end subroutine PLTID3
      SUBROUTINE ENDPLT
         PRINT * , ' ENDPLT NOT IMPLEMENTED'
         RETURN
      end subroutine ENDPLT
      SUBROUTINE NUMBER(X,Y,A,B,C,II)
         PRINT * , ' NUMBER NOT IMPLEMENTED'
         RETURN
      end subroutine NUMBER
      SUBROUTINE SYMBOL(X,Y,A,B,C,II)
         PRINT * , ' SYMBOL NOT IMPLEMENTED'
         RETURN
      end subroutine SYMBOL
      SUBROUTINE PLOT(X,Y,II)
         PRINT * , ' PLOT NOT IMPLEMENTED'
         RETURN
      end subroutine PLOT
      SUBROUTINE LINE(X,Y,NPTS,I,J,K)
         PRINT * , ' LINE NOT IMPLEMENTED'
         RETURN
      end subroutine LINE
