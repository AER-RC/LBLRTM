module voigt_module
! module contains functionality to create a voigt shape based on the
! lorenzian and dopplers shapes
! LBLRTM approach
! approximate vight as a sum of Gaussian, Lorenzian and correction term
! Humliceck - Voight mostly Pade approaximation 
implicit none
private

public LSF_VOIGT  ! compute Voigt using the Humlicek approximation
public LSF_VOIGTQ ! compute Voigt using the Humlicek approximation for Non-LTE case
public VOIGT_FLAG ! flag to use either Voigt shape (true) or LBL approximation (false)

logical :: dbg = .false.
logical :: VOIGT_FLAG = .true.

REAL, PARAMETER  :: PI = ACOS(-1.)
REAL, PARAMETER  :: SQRTLOG2   = sqrt(log(2.))
REAL, PARAMETER  :: SQRTLOG2PI = SQRTLOG2/sqrt(PI)

contains
!
! LBLRTM verion of CNVFN
!

Subroutine LSF_VOIGT(DV, VNU, SP, SPPSP, RECALF, ALPHAD, ALPHAL, VFT, &
                     R1, R2, R3, ILO, IHI, IOUT, MAX1, &
                     IPANEL, IDATA, ILAST, SDEP)
   implicit none
   integer, intent(inout)     :: ILO
   integer, intent(in)        :: IHI
   integer, intent(in)        :: MAX1
   integer, intent(in)        :: IOUT(250)
   integer, intent(in)        :: IDATA

   integer, intent(inout)     :: IPANEL
   integer, intent(inout)     :: ILAST

   real*8,  intent(in)        :: VNU(250)
   real*8,  intent(in)        :: VFT
   real,    intent(in)        :: DV
   real,    intent(in)        :: SP(250)
   real,    intent(in)        :: SPPSP(250)
   real,    intent(in)        :: RECALF(250)
   real,    intent(in)        :: ALPHAD(250)
   real,    intent(in)        :: ALPHAL(250)
   real, optional, intent(in) :: SDEP(250)

   ! real,    intent(inout)     :: R1(-4050:4050)
   ! real,    intent(inout)     :: R2(-1050:1050)
   ! real,    intent(inout)     :: R3(-300:300)
   real,    intent(inout)     :: R1(4050)
   real,    intent(inout)     :: R2(1050)
   real,    intent(inout)     :: R3(300)

   real*8  :: dVNU1 
   real*8  :: dVNU2 
   real*8  :: dVNU3 
   real*8  :: VNU_SHIFT

   real*8  :: d04
   real*8  :: d16
   real*8  :: d64
   
   real    :: RSHFT 
   real    :: BHWDXF 
   real    :: ZSLOPE 
   real    :: ZINT 
   real    :: DEPTHI 
   real    :: ALPHAV 
   real    :: DPTRAT
   real    :: STRF1
   real    :: STRF2
   real    :: STRF3

   integer :: J, I
   integer :: JMAX1
   integer :: JMIN1
   integer :: J1, J2, J3
   integer :: J2SHFT, J3SHFT


   real, parameter    :: HWF1 =  4.
   real, parameter    :: HWF2 = 16.
   real, parameter    :: HWF3 = 64.

   real, parameter    :: DXF1 = 0.002
   real               :: ZETA
   real               :: vg04
   real               :: vg16
   real               :: vg64
   real               :: ZF1L
   real               :: ZF2L
   real               :: ZF3L
   integer, parameter :: NFPTS = 2001
   real,    parameter :: CLC1  = HWF1/REAL(NFPTS-1)
   real,    parameter :: CLC2  = HWF2/REAL(NFPTS-1)
   real,    parameter :: CLC3  = HWF3/REAL(NFPTS-1)

   real,    parameter :: CONF2  = 1./ 4.
   real,    parameter :: CONF3  = 1./16.

   real               :: SDEP_LOC
   real               :: voigtShape, Y1

   ! real, parameter    :: DVR2 = 4.*DV
   ! real, parameter    :: DVR3 = 16.*DV

   real*8  :: X
   real    :: X0
   logical :: lineCoupling
   real :: QFN
   QFN(X, X0, DV) = (DV**2 + 2.*X0**2 - x**2)/(DV**2 + X0**2)

   IF (ILO.LE.IHI) THEN
      DO J = ILO, IHI
         I = IOUT(J)
         IF (SP(I).NE.0.) THEN
            ALPHAV = 1./RECALF(I)
            VNU_SHIFT = VNU(I)-VFT

            ZINT   = VNU_SHIFT/DV
            BHWDXF = HWF1*ALPHAV/DV
            JMAX1  = ZINT+BHWDXF+1.5

            IF (JMAX1.GT.MAX1) THEN
               ILAST  = J-1
               IPANEL = 1
               GO TO 40
            ENDIF

            DEPTHI = SP(I)
            DPTRAT = SP(I)*SPPSP(I)/ALPHAV
            lineCoupling = SPPSP(I).NE.0.

            JMIN1 = ZINT-BHWDXF+1.5

            RSHFT = SIGN(0.5, ZINT)

            J2SHFT =  3./ 4. * ZINT + RSHFT
            J3SHFT = 15./16. * ZINT + RSHFT
!
            d04 = HWF1*ALPHAV;  vg04 = SDVOIGT(d04, ALPHAL(I), ALPHAD(I), SDEP(I))
            d16 = HWF2*ALPHAV;  vg16 = SDVOIGT(d16, ALPHAL(I), ALPHAD(I), SDEP(I))
            d64 = HWF3*ALPHAV;  vg64 = SDVOIGT(d64, ALPHAL(I), ALPHAD(I), SDEP(I))

            DO J1    = JMIN1, JMAX1
               J2    =  J1 - J2SHFT
               J3    =  J1 - J3SHFT
               dVNU1 =     (J1-1)*DV - VNU_SHIFT
               dVNU2 =  4.*(J2-1)*DV - VNU_SHIFT
               dVNU3 = 16.*(J3-1)*DV - VNU_SHIFT
               !
               ! F1 defined -d04 +d04
               !
               if (ABS(dVNU1) <= d04) then
                  voigtShape = SDVOIGT(dVNU1, ALPHAL(I), ALPHAD(I), SDEP(I)) - &
                                               vg04*QFN(dVNU1, d04, DV)
                  if (lineCoupling) then
                     Y1 = DEPTHI + dVNU1 * DPTRAT
                  else
                     Y1 = DEPTHI
                  ENDIF
                  R1(J1) = R1(J1) + Y1*voigtShape
               end if
               !
               ! F2 defined -d16 +d16
               !
               if (ABS(dVNU2) <= d16) then
                  if (ABS(dVNU2) <= d04) then
                     ! this region  - approximation
                     voigtShape = vg04*QFN(dVNU2, d04, DV)
                  else 
                     ! this region  - accurate SD Voigt 
                     voigtShape = SDVOIGT(dVNU2, ALPHAL(I), ALPHAD(I), SDEP(I))
                  endif
                  ! remove next region approximation
                  voigtShape = voigtShape - vg16*QFN(dVNU2, d16, DV)
                  if (lineCoupling) then
                     Y1 = DEPTHI + dVNU2 * DPTRAT
                  else
                     Y1 = DEPTHI
                  ENDIF
                  R2(J2) = R2(J2) + Y1*voigtShape
               ENDIF
               !
               ! F3 defined -d64 +d64
               !
               if (ABS(dVNU3) <= d64) then
                  if (ABS(dVNU3) <= d16) then
                     voigtShape = vg16*QFN(dVNU3, d16, DV)
                  else 
                     voigtShape = SDVOIGT(dVNU3, ALPHAL(I), ALPHAD(I), SDEP(I))
                  endif
                  voigtShape = voigtShape - vg64*QFN(dVNU3, d64, DV)
                  if (lineCoupling) then
                     Y1 = DEPTHI + dVNU3 * DPTRAT
                  else
                     Y1 = DEPTHI
                  ENDIF
                  R3(J3) = R3(J3) + Y1 * voigtShape
               endif   
            END DO
         ENDIF
      ENDDO
      ILAST = IHI
!
!        IDATA=0 FOR MORE DATA REQUIRED
!        IDATA=1 IF NO MORE DATA REQUIRED
!
      IPANEL = IDATA
   ELSE
      IPANEL = 1
   ENDIF
!
40 ILO = ILAST+1
   RETURN
END Subroutine LSF_VOIGT
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
FUNCTION VOIGT(DELTNU,ALPHAL,ALPHAD)
   REAL*8, intent(in) :: DELTNU
   REAL,   intent(in) :: ALPHAL
   REAL,   intent(in) :: ALPHAD

   REAL               :: VOIGT
   ! local
   REAL X, Y
   !---case of pure lorent fast return
   if ( ALPHAD < epsilon(1.0) ) then
      VOIGT = ALPHAL/PI/(ALPHAL**2 + DELTNU**2)
   else
     !---SETUP PARAMETERS FOR CALL TO VOIGT GENERATOR (HUMLICEK)
     x = SQRTLOG2*DELTNU/ALPHAD
     y = SQRTLOG2*ALPHAL/ALPHAD

     !---CALL the Humlicek subroutine
     VOIGT  = SQRTLOG2PI/ALPHAD*W4(x,y)
   end if
   RETURN
end FUNCTION VOIGT
!------------------------------------------------------------------------
! 
!  speed dependent voigt
!
!
FUNCTION SDVOIGT(DELTNU,ALPHAL,ALPHAD,SDEP)
      REAL*8,           intent(in) :: DELTNU
      REAL,             intent(in) :: ALPHAL
      REAL,             intent(in) :: ALPHAD
      REAL,   optional, intent(in) :: SDEP

      REAL SDVOIGT

      REAL v
      REAL ZETA,ALPHAV,DNU
      REAL AL,AD,ANORM1,DZETA
      !MJA 20130517 Speed Dependence
      REAL alfa, beta, delta, temp, x1, x2, y1, y2
      REAL delta2, gamma2
      REAL alfadelta
      REAL, parameter :: TINY = 1.0e-4
      LOGICAL  :: isSD       ! speed dependent Voigt flag
      LOGICAL  :: isLorentz  ! Lorentz flag
      real, parameter :: PI = acos(-1.0)

      isSD = present(SDEP)
      isLorentz = ((ABS(deltnu) > (100.*alphad)).or.( alphal > (99. * alphad)))
      if (isSD) then
         isSD =  (ABS(SDEP) > TINY) .and. (.not. isLorentz)
      endif

      !---SETUP PARAMETERS FOR CALL TO VOIGT GENERATOR (HUMLICEK)
      IF (isSD) THEN !SD Voigt

         !---Speed Dependence follows Boone et al., 2011, An efficient analytical
         !---approach for calculating line mixing in atmospheric remote sensing
         !---applications, JQSRT, 112, 980-989.

         !SDEP is Benner's speed dependence definition
         !Therefore Boone's gamma2 = alphal*SDEP and his alfa = (1/SDEP) - 1.5
         gamma2    = 2.0*alphal*SDEP
         alfa      = 1./SDEP - 1.5                       !Boone et al., 2011 Eq 13
         beta      = deltnu/gamma2                       !Boone et al., 2011 Eq 13
         delta     = alphad/gamma2/SQRTLOG2              !Boone et al., 2011 Eq 13
         delta2    = delta**2                            !Boone et al., 2011 Eq 13
         alfadelta = (alfa + delta**2)/2.0

         !Boone et al., 2011, Eq. 12
         temp = sqrt(alfadelta**2 + beta**2) !"temp" means temporary, used below
         !Real part of z1, z2, Boone et al., 2011, Eq. 12
         x2 = sqrt(temp + alfadelta)
         x1 = x2 - delta
         x2 = x2 + delta
         !Imag part of z1, z2, Boone et al., 2011, Eq. 12
         if (beta > 0.0) then
            y1 =  sqrt(temp-alfadelta)
         else if (beta < 0.0) then
            y1 = -sqrt(temp-alfadelta)
         else
            y1 = 0.0
         endif

         !---CALL the modified Humlicek subroutine to calc speed dependent Voigt
         v=SD_Humlicek(y1,x1,y1,x2)
         IF (v < 0.0) STOP "SD_Humlicek :: negative v detected"

      ELSE !Normal Voigt
         !---case of pure lorentz
         if (isLorentz) then
            SDVOIGT = ALPHAL/PI/(ALPHAL**2 + DELTNU**2)
            RETURN
         end if

         !---SETUP PARAMETERS FOR CALL TO VOIGT GENERATOR (HUMLICEK)
         !-----GENERATE VOIGT PROFILE SUCH THAT THE VOIGT HALFWIDTH = 1.
         x1 = SQRTLOG2*DELTNU/ALPHAD
         y1 = SQRTLOG2*ALPHAL/ALPHAD
         !---CALL the Humlicek subroutine
         v = W4(x1, y1)
      ENDIF
      SDVOIGT= v*SQRTLOG2PI/ALPHAD !Boone et al., 2011 Eq 10
      if (isNAN(SDVOIGT)) then
         print *, "SDVOIGT::error:: ", SDVOIGT, DELTNU, ALPHAL, ALPHAD
         STOP "NAN detected"
      endif
      RETURN
   end FUNCTION SDVOIGT

!*************************************************************************
!     FORTRAN function for the complex probability function w(z).
!     COMPUTES THE COMPLEX PROBABILITY FUNCTION W(Z)=EXP(-Z**2)*ERFC(-I*Z)
!     IN THE UPPER HALF-PLANE Z=X+I*Y (I.E. FOR Y>=0)
!     MAXIMUM RELATIVE ERROR OF BOTH REAL AND IMAGINARY PARTS IS <1*10**(-4)
!     (current version return just the REAl part)
!     Reference:
!     Humlicek, J., 1982; Optimized Computation of the Voigt and Complex
!     Probability Functions, J. Quant. Spectrosc. Radiat. Transfer, 27, 437-444.
!     S-A Boukabara AER Inc, 2000
!*************************************************************************
function w4_region1(X,Y)
   real, intent(in) :: x
   real, intent(in) :: y
   real             :: w4_region1
   ! local
   REAL             :: XY

   ! T=CMPLX(Y,-X)
   ! CW4=T*.5641896/(.5+T*T)

   XY = 0.5 + Y**2 + X**2
   w4_region1 = .5641896*Y*XY/(XY**2 - 2.*X**2)

end function w4_region1
!-----------------------------------
function w4_region2(X,Y)
   real, intent(in) :: x
   real, intent(in) :: y
   real             :: w4_region2
   ! local
   COMPLEX :: T,U
   COMPLEX :: CW4

   T   = CMPLX(Y,-X)
   U   = T*T
   CW4 = T*(1.410474+U*.5641896)/(.75+U*(3.+U))
   w4_region2 = real(CW4)

end function w4_region2
!-----------------------------------
function w4_region3(X,Y)
   real, intent(in) :: x
   real, intent(in) :: y
   real             :: w4_region3
   ! local
   COMPLEX :: T
   COMPLEX :: CW4

   T=CMPLX(Y,-X)
   CW4=(16.4955+T*(20.20933+T*(11.96482+ &
     T*(3.778987+T*.5642236))))/ &
     (16.4955+T*(38.82363+T*(39.27121+ &
     T*(21.69274+T*(6.699398+T)))))

   w4_region3 = real(CW4)

end function w4_region3
!-----------------------------------
function w4_region4(X,Y)
   real, intent(in) :: x
   real, intent(in) :: y
   real             :: w4_region4
   ! local
   COMPLEX :: T, U
   COMPLEX :: CW4

   T=CMPLX(Y,-X)
   U=T*T
   CW4=CEXP(U)-T*(36183.31-U*(3321.9905- &
     U*(1540.787-U*(219.0313-U* &
     (35.76683-U*(1.320522-U*.56419))))))/ &
     (32066.6-U*(24322.84-U* &
     (9022.228-U*(2186.181-U*(364.2191- &
     U*(61.57037-U*(1.841439-U)))))))

   w4_region4 = real(CW4)

end function w4_region4
!-----------------------------------

FUNCTION W4(x,y)
   real, intent(in) :: x
   real, intent(in) :: y
   real             :: W4

   REAL    :: S

   S=ABS(X)+Y
   IF (S >= 15.) THEN 
     !     ***   REGION I
     if (dbg) write(900, *) 1 
     W4 = w4_region1(X,Y)
     RETURN

   ELSE IF(S >= 5.5) THEN
     if (dbg) write(900, *) 2 
     !     ***   REGION II
     W4 = w4_region2(X,Y)
   
   ELSE IF (Y >= 0.195*ABS(X)-0.176) THEN
     if (dbg) write(900, *) 3 
     !     ***   REGION III
     W4 = w4_region3(X,Y)
   ELSE
     !     ***   REGION IV
     if (dbg) write(900, *) 4 
     W4 = w4_region4(X,Y)
   ENDIF
   RETURN
END FUNCTION W4
!*************************************************************************
function get_region(y, x) result(r)
   real, intent(in) :: y
   real, intent(in) :: x

   integer          :: r
   ! local
   real :: xx1
   real :: s

   s = abs(x) + y
   xx1 = 0.195*ABS(x)-0.176

   IF(s >= 15.0) THEN
      r = 1
   ELSE IF (s >= 6.0) THEN
      r = 2
   ELSE 
      IF (y < x) THEN
         r = 4
      ELSE
         r = 3
      ENDIF
   ENDIF
end function get_region

!*************************************************************************
!     FORTRAN function for the difference of two complex probability
!     functions w(z1)-w(z2)..
!     COMPUTES THE COMPLEX PROBABILITY FUNCTION W(Z)=EXP(-Z**2)*ERFC(-I*Z)
!     IN THE UPPER HALF-PLANE Z=X+I*Y (I.E. FOR Y>=0)
!     MAXIMUM RELATIVE ERROR OF BOTH REAL AND IMAGINARY PARTS IS <1*10**(-4)
!     Modified for speed dependece as recommended in Boone et al., 2011.
!     References:
!     Humlicek, J., 1982; Optimized Computation of the Voigt and Complex
!     Probability Functions, J. Quant. Spectrosc. Radiat. Transfer, 27, 437-444.
!     Boone et al., 2011, An efficient analytical
!     approach for calculating line mixing in atmospheric remote sensing
!     applications, JQSRT, 112, 980-989.
!     Original Humlicek Implementation:
!     S-A Boukabara AER Inc, 2000
!     Modifications for speed dependence:
!     M. J. Alvarado, AER, 2013
!*************************************************************************
   FUNCTION SD_Humlicek(x1,y1,x2,y2)
      REAL, intent(in) :: X1
      REAL, intent(in) :: Y1
      REAL, intent(in) :: X2
      REAL, intent(in) :: Y2
      REAL             :: SD_Humlicek
    
      REAL    :: W1, W2
      REAL    :: XX1, XX2
      INTEGER :: Region1, Region2, Region

!     Find correct Humlicek region for each line
      Region1 = get_region(Y1, X1)
      Region2 = get_region(Y2, X2)

      !Use Largest of two regions
      Region = MAX(Region1, Region2)
      IF(Region == 1) THEN
        !     ***   REGION I
        W1=w4_region1(X1,Y1)
        W2=w4_region1(X2,Y2)
      ELSE IF(Region == 2) THEN
        !     Change in Region 2 and 3 boundary recommended on page 985 of Boone et al. 2011
        !     (see 4th paragraph of second column)
        !     ***   REGION II
        W1=w4_region2(X1,Y1)
        W2=w4_region2(X2,Y2)
      ELSE IF (Region == 3) THEN
        !     ***   REGION III
        W1=w4_region3(X1,Y1)
        W2=w4_region3(X2,Y2)
      ELSE
  !     ***   REGION IV
  !     You get errors if you use Region IV approximations for lines outside
  !     Region IV, so use Region III instead (MJA, 08062013)
        IF(Region1 == 4) THEN
           W1=w4_region4(X1,Y1)
        ELSE
           W1=w4_region3(X1,Y1)
        ENDIF
        IF(Region2 == 4) THEN
           W2=w4_region4(X2,Y2)
        ELSE
           W2=w4_region3(X2,Y2)
        ENDIF
      ENDIF
      SD_Humlicek = W1 - W2
      RETURN
   END FUNCTION SD_Humlicek

!-----------------------------------------------------------------------
!     SUBROUTINE CNVFNV PERFORMS THE CONVOLUTION OF THE LINE DATA WITH
!     THE VOIGT LINE SHAPE (APPROXIMATED)
!
!     IMPLEMENTATION:    R.D. WORSHAM
!
!     ALGORITHM REVISIONS:    S.A. CLOUGH
!     R.D. WORSHAM
!     J.L. MONCET
!
!
!     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139
!
!     WORK SUPPORTED BY:    THE ARM PROGRAM
!     OFFICE OF ENERGY RESEARCH
!     DEPARTMENT OF ENERGY
!
!
!     SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!     FASCOD3
!
!-----------------------------------------------------------------------
      SUBROUTINE LSF_VOIGTQ( VNU,SABS,SRAD,SPPSP,RECALF,&
                         R1,R2,R3,RR1,RR2,RR3,lbR1,lbR2,lbR3, &
                         DV,DVR2,DVR3,HWF1,HWF2,HWF3,DXF1,NX1,NX2,NX3,&
                         VFT,ILO,IHI,MAX1,IOUT,IPANEL,IDATA, &
                         ALPHAL, ALPHAD, SDEP )

!-----------------------------------------------------------------------
      IMPLICIT NONE !REAL*8           (V)

      real*8    ,intent(in)    :: VNU(:)
      real,      intent(in)    :: SABS(:)
      real,      intent(in)    :: SRAD(:)
      real,      intent(in)    :: SPPSP(:)
      real,      intent(in)    :: RECALF(:)
      real      ,intent(inout) :: R1(lbR1:), R2(lbR2:), R3(lbR3:)
      real      ,intent(inout) :: RR1(lbR1:), RR2(lbR2:), RR3(lbR3:)
      integer   ,intent(in)    :: lbR1,lbR2,lbR3
      !
      real      ,intent(in)    :: DV, DVR2, DVR3
      real      ,intent(in)    :: HWF1
      real      ,intent(in)    :: HWF2
      real      ,intent(in)    :: HWF3
      real      ,intent(in)    :: DXF1
      integer   ,intent(in)    :: NX1, NX2, NX3
      real*8    ,intent(in)    :: VFT
      integer   ,intent(inout) :: ILO,IHI
      integer   ,intent(in)    :: MAX1
      integer   ,intent(in)    :: IOUT(250)
      integer   ,intent(out)   :: IPANEL
      integer   ,intent(in)    :: IDATA
      real      ,intent(in)    :: ALPHAL(250)    ! Lorentz half-width
      real      ,intent(in)    :: ALPHAD (250)   ! Doppler half-width
      real      ,intent(in)    :: SDEP(250)      ! speed dependence


      !--- Local variables
      INTEGER  :: I,      ILAST,  IZ1,    IZ2
      INTEGER  :: IZ3,    IZM,    J,      J1
      INTEGER  :: J2,     J2SHFT, J3,     J3SHFT
      INTEGER  :: JMAX1,  JMIN1,  JMIN2,  JMIN3
      REAL     :: BHWDXF, CLC1,   CLC2,   CLC3
      REAL     :: CONF2,  CONF3,  DEPTHA, DEPTHR
      REAL     :: DPTRAT, HWDXF,  RSHFT,  STRDA
      REAL     :: STRDR,  STRFA1, STRFA2, STRFA3
      REAL     :: STRFR1, STRFR2, STRFR3, STRVRA
      REAL     :: STRVRR, WAVDXF, ZETDIF
      REAL     :: ZF1,    ZF1L,   ZF2,    ZF2L
      REAL     :: ZF3,    ZF3L,   ZINT,   ZSLOPE

      real     :: voigtShape, Y1
      logical  :: lineCoupling
      logical  :: nlteFlag
      real     :: ALPHAV 

      real    :: vg04
      real    :: vg16
      real    :: vg64
      real*8  :: d04
      real*8  :: d16
      real*8  :: d64
      real*8  :: dVNU1 
      real*8  :: dVNU2 
      real*8  :: dVNU3 
      real*8  :: VNU_SHIFT

      real :: QFN, X, X0 
      QFN(X, X0, DV) = (DV**2 +2.*X0**2 - x**2)/(DV**2+X0**2)

      CLC1 =  4./( REAL(NX1-1))
      CLC2 = 16./( REAL(NX2-1))
      CLC3 = 64./( REAL(NX3-1))
      WAVDXF = DV/DXF1
      HWDXF = HWF1/DXF1
      CONF2 = DV/DVR2
      CONF3 = DV/DVR3
      ILAST  = ILO-1

      IF (ILO.LE.IHI) THEN
         DO J = ILO, IHI
            I = IOUT(J)
            IF (SABS(I).NE.0.) THEN

               ALPHAV    = 1./RECALF(I)
               VNU_SHIFT = VNU(I)-VFT

               ZINT   = VNU_SHIFT/DV
               BHWDXF = HWF1*ALPHAV/DV
               JMAX1  = ZINT+BHWDXF+1.5
               ! 
               ! quick exit
               !
               IF (JMAX1.GT.MAX1) THEN
                  ILAST  = J-1
                  IPANEL = 1
                  GO TO 40
               ENDIF

               DEPTHA = SABS(I)
               DEPTHR = SRAD(I)
               DPTRAT = SPPSP(I)/ALPHAV
               
               nlteFlag     = DEPTHR   /= 0.
               lineCoupling = SPPSP(I) /= 0.

               JMIN1 = ZINT-BHWDXF+1.5

               RSHFT = SIGN(0.5, ZINT)

               J2SHFT =  3./ 4. * ZINT + RSHFT
               J3SHFT = 15./16. * ZINT + RSHFT
               !
               d04 = HWF1*ALPHAV;  vg04 = SDVOIGT(d04, ALPHAL(I), ALPHAD(I), SDEP(I))
               d16 = HWF2*ALPHAV;  vg16 = SDVOIGT(d16, ALPHAL(I), ALPHAD(I), SDEP(I))
               d64 = HWF3*ALPHAV;  vg64 = SDVOIGT(d64, ALPHAL(I), ALPHAD(I), SDEP(I))
               !
               DO J1 = JMIN1, JMAX1
                  dVNU1 = (J1-1)*DV - VNU_SHIFT
                  !
                  ! F1 defined -d04 +d04
                  !
                  if (ABS(dVNU1) <= d04) then
                     voigtShape = (SDVOIGT(dVNU1, ALPHAL(I), ALPHAD(I), SDEP(I)) - &
                                                  vg04*QFN(dVNU1, d04, DV))
                     if (lineCoupling) then
                        Y1 = (1. + dVNU1 * DPTRAT)*voigtShape
                     else
                        Y1 = voigtShape
                     ENDIF
                     R1(J1) = R1(J1) + DEPTHA*Y1
                     if (nlteFlag) RR1(J1) = RR1(J1)+DEPTHR*Y1
                  end if
                  !
                  ! F2 defined -d16 +d16
                  !
                  J2    = J1 - J2SHFT
                  dVNU2 =  4.*(J2-1)*DV - VNU_SHIFT
                  if (ABS(dVNU2) <= d16) then
                     if (ABS(dVNU2) <= d04) then
                        voigtShape = vg04*QFN(dVNU2, d04, DV)
                     else 
                        voigtShape = SDVOIGT(dVNU2, ALPHAL(I), ALPHAD(I), SDEP(I))
                     endif
                     voigtShape = voigtShape - vg16*QFN(dVNU2, d16, DV)
                     if (lineCoupling) then
                        Y1 = (1. + dVNU2 * DPTRAT)*voigtShape
                     else
                        Y1 = voigtShape
                     ENDIF
                     R2(J2) = R2(J2) + DEPTHA*Y1
                     if (nlteFlag) RR2(J2) = RR2(J2)+DEPTHR*Y1
                  ENDIF
                  !
                  ! F3 defined -d64 +d64
                  !
                  J3    = J1 - J3SHFT
                  dVNU3 = 16.*(J3-1)*DV - VNU_SHIFT
                  if (ABS(dVNU3) <= d64) then
                     if (ABS(dVNU3) <= d16) then
                        voigtShape = vg16*QFN(dVNU3, d16, DV)
                     else 
                        voigtShape = SDVOIGT(dVNU3, ALPHAL(I), ALPHAD(I), SDEP(I))
                     endif
                     voigtShape = voigtShape - vg64*QFN(dVNU3, d64, DV)
                     if (lineCoupling) then
                        Y1 = (1. + dVNU3 * DPTRAT)*voigtShape
                     else
                        Y1 = voigtShape
                     ENDIF
                     R3(J3) = R3(J3) + DEPTHA*Y1
                     if (nlteFlag) RR3(J3) = RR3(J3)+DEPTHR*Y1
                  endif   
               END DO
            ENDIF

         END DO
         ILAST = IHI

!        IDATA=0 FOR MORE DATA REQUIRED
!        IDATA=1 IF NO MORE DATA REQUIRED

         IPANEL = IDATA
      ELSE
         IPANEL = 1
      ENDIF

   40 ILO = ILAST+1
      !CALL CPUTIM (TIME)
      !TIMCNV = TIMCNV+TIME-TIME0
      RETURN

      END  SUBROUTINE LSF_VOIGTQ

end module voigt_module
