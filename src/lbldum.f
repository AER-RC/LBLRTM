C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
C
C  --------------------------------------------------------------------------
C |                                                                          |
C |  Copyright 2002, 2003, Atmospheric & Environmental Research, Inc. (AER). |
C |  This software may be used, copied, or redistributed as long as it is    |
C |  not sold and this copyright notice is reproduced on each copy made.     |
C |  This model is provided as is without any express or implied warranties. |
C |                       (http://www.rtweb.aer.com/)                        |
C |                                                                          |
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
