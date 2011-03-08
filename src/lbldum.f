C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
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
C     SUBROUTINE NONLTE (MPTS)
C     PRINT * , ' NONLTE NOT IMPLEMENTED'
C     STOP
C     END

      SUBROUTINE LASER(VLAS,MFILE,JAERSL)
      PRINT * , ' LASER NOT IMPLEMENTED'
      STOP
      END
      FUNCTION RANDM(IRAND)
      PRINT * , ' RANDM NOT IMPLEMENTED PROPERLY'
      RANDM=0.5
      IRAND=IABS(IRAND)
      RETURN
      END
      SUBROUTINE PLTID3(PROGID,XSIZ,YSIZ,SCAL)
      DIMENSION PROGID(3)
      PRINT * , ' PLTID3 NOT IMPLEMENTED'
      RETURN
      END
      SUBROUTINE ENDPLT
      PRINT * , ' ENDPLT NOT IMPLEMENTED'
      RETURN
      END
      SUBROUTINE NUMBER(X,Y,A,B,C,II)
      PRINT * , ' NUMBER NOT IMPLEMENTED'
      RETURN
      END
      SUBROUTINE SYMBOL(X,Y,A,B,C,II)
      PRINT * , ' SYMBOL NOT IMPLEMENTED'
      RETURN
      END
      SUBROUTINE PLOT(X,Y,II)
      PRINT * , ' PLOT NOT IMPLEMENTED'
      RETURN
      END
      SUBROUTINE LINE(X,Y,NPTS,I,J,K)
      PRINT * , ' LINE NOT IMPLEMENTED'
      RETURN
      END
