C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      SUBROUTINE NONLTE (MPTS)                                                  
      PRINT * , ' NONLTE NOT IMPLEMENTED'                                       
      STOP                                                                      
      END                                                                       
      SUBROUTINE XLAYMS(MPS,NPTS,LFILE,MFILE,NFILE)                             
      PRINT * , ' XLAYMS NOT IMPLEMENTED'                                       
      STOP                                                                      
      END                                                                       
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
