module pCqSDHC_module
implicit none
private

public:: pCqSDHC

complex, parameter :: zone = complex(1.d0,0.D0)
complex, parameter :: zi   = complex(0.d0,1.D0)
integer, parameter      :: tt_len = 15
real, parameter :: tt(tt_len)= (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5/)
real, parameter :: pipwoeronehalf = 0.564189583547756d0

contains    
    subroutine  pCqSDHC(&
                   sg0,            &!<--: Unperturbed line position in cm-1 (Input).
                   GamD,           &!<--: Doppler HWHM in cm-1 (Input)
                   Gam0,           &!<--: Speed-averaged line-width in cm-1 (Input).    
                   Gam2,           &!<--: Speed dependence of the line-width in cm-1 (Input).
                   anuVC,          &!<--: Velocity-changing frequency in cm-1 (Input).
                   eta,            &!<--: Correlation parameter, No unit (Input).
                   Shift0,         &!<--: Speed-averaged line-shift in cm-1 (Input).
                   Shift2,         &!<--: Speed dependence of the line-shift in cm-1 (Input)     
                   sg,             &!<--: Current WaveNumber of the Computation in cm-1 (Input).
                   LS_pCqSDHC_R,   &!-->: Real part of the normalized spectral shape (cm) (Output)
                   LS_pCqSDHC_I    &!-->: Imaginary part of the normalized spectral shape (cm) (Output)
                   )

!-------------------------------------------------
!   "pCqSDHC": partially-Correlated quadratic-Speed-Dependent Hard-Collision
!   Subroutine to Compute the complex normalized spectral shape of an 
!   isolated line by the pCqSDHC model following the two references:
!   [1] Ngo NH, Lisak D, Tran H, Hartmann J-M. An isolated line-shape model
!   to go beyond the Voigt profile in spectroscopic databases and radiative 
!   transfer codes. J Quant Radiat Transfer 2013;129:89-100.        
!   [2] Tran H, Ngo NH, Hartmann J-M. Efficient computation of some speed-dependent 
!   isolated line profiles. J Quant Radiat Transfer 2013;129:199-203.
!
!   Input/Output Parameters of Routine (Arguments or Common)
!   ---------------------------------
!   T       : Temperature in Kelvin (Input).
!   amM1    : Molar mass of the absorber in g/mol(Input).
!   sg0     : Unperturbed line position in cm-1 (Input).
!   GamD    : Doppler HWHM in cm-1 (Input)
!   Gam0    : Speed-averaged line-width in cm-1 (Input).    
!   Gam2    : Speed dependence of the line-width in cm-1 (Input).
!   anuVC   : Velocity-changing frequency in cm-1 (Input).
!   eta     : Correlation parameter, No unit (Input).
!   Shift0  : Speed-averaged line-shift in cm-1 (Input).
!   Shift2  : Speed dependence of the line-shift in cm-1 (Input)     
!   sg      : Current WaveNumber of the Computation in cm-1 (Input).
!
!   Output Quantities (through Common Statements)
!   -----------------
!   LS_pCqSDHC_R: Real part of the normalized spectral shape (cm)
!   LS_pCqSDHC_I: Imaginary part of the normalized spectral shape (cm)
!
!   Called Routines: 'CPF'  (Complex Probability Function)
!   ---------------  'CPF3' (Complex Probability Function for the region 3)
!
!-------------------------------------------------
    real, intent(in)  :: sg0            !<--: Unperturbed line position in cm-1 (Input).
    real, intent(in)  :: GamD           !<--: Doppler HWHM in cm-1 (Input)
    real, intent(in)  :: Gam0           !<--: Speed-averaged line-width in cm-1 (Input).    
    real, intent(in)  :: Gam2           !<--: Speed dependence of the line-width in cm-1 (Input).
    real, intent(in)  :: anuVC          !<--: Velocity-changing frequency in cm-1 (Input).
    real, intent(in)  :: eta            !<--: Correlation parameter, No unit (Input).
    real, intent(in)  :: Shift0         !<--: Speed-averaged line-shift in cm-1 (Input).
    real, intent(in)  :: Shift2         !<--: Speed dependence of the line-shift in cm-1 (Input)     
    real, intent(in)  :: sg             !<--: Current WaveNumber of the Computation in cm-1 (Input).
    real, intent(out) :: LS_pCqSDHC_R   !-->: Real part of the normalized spectral shape (cm) (Output)
    real, intent(out) :: LS_pCqSDHC_I   !-->: Imaginary part of the normalized spectral shape (cm) (Output)

    ! local variables
    real    :: cte
    real    :: xz1,xz2,yz1,yz2,xXb,yXb
    real    :: wr1,wi1,wr2,wi2,wrb,wib
    real    :: SZ1,SZ2,DSZ,SZmx,SZmn
    complex :: c0,c2,c0t,c2t
    complex :: X,Y,Z1,Z2,csqrtY
    complex :: Aterm,Bterm,LS_pCqSDHC

    real, parameter ::  pi = acos(-1.)
    real, parameter :: rpi = sqrt(pi)
!
!-------------------------------------------------
!
    cte = sqrt(log(2.))/GamD
! Calculating the different parameters 
    c0  = cmplx(Gam0,Shift0)
    c2  = cmplx(Gam2,Shift2)
    c0t = (1.-eta)*(c0-1.5*c2)+anuVC
    c2t = (1.-eta)*c2
!
    if (abs(C2t) == 0.) then
        ! when C2t=0
        Z1 = (zi*(sg0-sg)+C0t)*cte
        xZ1=-imag(Z1)
        yZ1= real(Z1)
        Call CPF ( xZ1, yZ1, WR1, WI1 )
        Aterm=rpi*cte*cmplx(WR1,WI1)
        if (abs(Z1) <= 4e3) then
            Bterm=rpi*cte*((1.-Z1**2)*cmplx(WR1,WI1)+Z1/rpi)
        else
            Bterm=cte*(rpi*cmplx(WR1,WI1)+0.5/Z1-0.75/(Z1**3))
        endif
    else   
        X = (zi*(sg0-sg)+c0t)/c2t
        Y = 1./(2.*cte*C2t)**2          
        csqrtY=(Gam2-zi*Shift2)/(2.d0*cte*(1.-eta)*(Gam2**2+Shift2**2))
        
        if (abs(X) <=  3e-8 *abs(Y)) then
        ! when abs(Y) is much larger than abs(X)
            Z1 = (zi*(sg0-sg)+C0t)*cte
            Z2 = sqrt(X+Y)+csqrtY
            xZ1=-imag(z1)
            yZ1= real(z1)
            xZ2=-imag(z2)
            yZ2= real(z2)
            Call CPF ( xZ1, yZ1, WR1, WI1 )
            Call CPF ( xZ2, yZ2, WR2, WI2 ) 
            Aterm=rpi*cte*(cmplx(wr1,wi1) - cmplx(wr2,wi2))
            Bterm=( -1.0 + &
                   rpi/(2.*csqrtY)*(1.-Z1**2)*cmplx(wr1,wi1)- &
                   rpi/(2.*csqrtY)*(1.-Z2**2)*cmplx(wr2,wi2) )/C2t
        else    
            if (abs(Y) >  1e-15*abs(X)) then
            ! calculating Z1 and Z2
                Z1=sqrt(X+Y)-csqrtY
                Z2=Z1+2.*csqrtY
            ! calculating the real and imaginary parts of Z1 and Z2
                xZ1 = -imag(Z1)
                yZ1 =  real(Z1)
                xZ2 = -imag(Z2)
                yZ2 =  real(Z2)
            ! check if Z1 and Z2 are close to each other
                SZ1 = sqrt(xZ1**2+yZ1**2)
                SZ2 = sqrt(xZ2**2+yZ2**2)
                DSZ = abs(SZ1-SZ2)
                SZmx= max(SZ1, SZ2)
                SZmn= min(SZ1, SZ2)
            ! when Z1 and Z2 are close to each other, ensure that they are in 
            ! the same interval of CPF 
                if ( (DSZ <= 1.).and.(SZmx > 8.).and.(SZmn <= 8.0) )  then
                    Call CPF3 ( xZ1, yZ1, WR1, WI1 ) 
                    Call CPF3 ( xZ2, yZ2, WR2, WI2 ) 
                else    
                    Call CPF ( xZ1, yZ1, WR1, WI1 ) 
                    Call CPF ( xZ2, yZ2, WR2, WI2 ) 
                endif
            ! calculating the A and B terms of the profile
                Aterm=rpi*cte*(cmplx(wr1,wi1) - cmplx(wr2,wi2))
                Bterm=( -1. + &
                       rpi/(2.*csqrtY)*(1.-Z1**2)*cmplx(wr1,wi1)- &
                       rpi/(2.*csqrtY)*(1.-Z2**2)*cmplx(wr2,wi2) )/C2t
            else
            ! when abs(X) is much larger than abs(Y)
                xZ1 =-imag(sqrt(X+Y))
                yZ1 = real(sqrt(X+Y))
                Call CPF ( xZ1, yZ1, WR1, WI1 ) 
                if (abs(sqrt(X)) <= 4.e3) then
                  xXb =-imag(sqrt(X))
                  yXb = real(sqrt(X))
                  Call CPF ( xXb, yXb, WRb, WIb ) 
                  Aterm=(2.*rpi/C2t)*(1./rpi-sqrt(X)*cmplx(WRb,WIb))
                  Bterm=(1./C2t)*(-1.+ &
                         2.*rpi*(1.-X-2.*Y)*(1./rpi-sqrt(X)*cmplx(wrb,wib))+ &
                         2.*rpi*sqrt(X+Y)*cmplx(wr1,wi1))
            ! and when abs(X) is much larger than 1
                else
                  Aterm=(1.0/C2t)*( 1./X-1.5/(X**2))
                  Bterm=(1.0/C2t)*(-1.+(1.-X-2.*Y)*(1./X-1.5/(X**2))+ &
                                   2.*rpi*sqrt(X+Y)*cmplx(wr1,wi1))
                endif
            endif
        endif
    endif
!
10    continue
!
    LS_pCqSDHC =(1./pi)*(Aterm/(1.-(anuVC-eta*(C0-1.5*C2))*Aterm+eta*C2*Bterm))

    LS_pCqSDHC_R = real(LS_pCqSDHC)
    LS_pCqSDHC_I = imag(LS_pCqSDHC)
   
    Return
    End Subroutine pCqSDHC

    Subroutine CPF(X,Y,WR,WI)
!-------------------------------------------------
! "CPF": Complex Probability Function
!  translated cpf function of hapi.py 
! .........................................................
!         .       Subroutine to Compute the Complex       .
!         .        Probability Function W(z=X+iY)         .
!         .     W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0      .
!         .    Which Appears when Convoluting a Complex   .
!         .     Lorentzian Profile by a Gaussian Shape    .
!         .................................................
!
!             WR : Real Part of W(z)
!             WI : Imaginary Part of W(z)
!
! This Routine was Taken from the Paper by J. Humlicek, which 
! is Available in Page 309 of Volume 21 of the 1979 Issue of
! the Journal of Quantitative Spectroscopy and Radiative Transfer
! Please Refer to this Paper for More Information
!
!-------------------------------------------------
!      
    Implicit None
    real, intent(in)  :: X, Y
    real, intent(out) :: WR,WI

    Integer :: I
    complex :: zm1,zm2,zterm,zsum
    real    :: Y1,Y2,Y3,R,R2,D,D1,D2,D3,D4
!   
    integer, parameter :: sz = 6   
    real, parameter :: T(sz) = (/0.314240376, 0.947788391,  1.59768264,   2.27950708,    3.02063703,    3.8897249/)
    real, parameter :: U(sz) = (/1.01172805 ,-0.75197147,   1.2557727e-2, 1.00220082e-2,-2.42068135e-4, 5.00848061e-7/)
    real, parameter :: S(sz) = (/1.393237   , 0.231152406, -0.155351466,  6.21836624e-3, 9.19082986e-5,-6.27525958e-7/)

! new Region 3
    if(sqrt(x**2 + y**2) > 8.) then
        call CPF3(X,Y,WR,WI)
        return
    end if
!
      WR = 0.
      WI = 0.
      Y1 = Y+1.5
      Y2 = Y1**2
      If( (Y > 0.85) .OR. (ABS(X) < (18.1*Y+1.65)) ) then
    !
    !       Region 1
    !
          DO I=1,sz
              R=X-T(I)
              D=1./(R**2+Y2)
              D1=Y1*D
              D2=R*D
              R=X+T(I)
              D=1./(R**2+Y2)
              D3=Y1*D
              D4=R*D
              WR=WR+U(I)*(D1+D3)-S(I)*(D2-D4)
              WI=WI+U(I)*(D2+D4)+S(I)*(D1-D3)
          END DO
      ELSE  
!
!       Region 2
!
          If (ABS(X) < 12.) WR = EXP(-X**2)
          Y3=Y+3.
          DO I=1,sz
              R =X-T(I)
              R2=R**2
              D =1./(R2+Y2)
              D1=Y1*D
              D2=R*D
              WR=WR+Y*(U(I)*(R*D2-1.5*D1)+S(I)*Y3*D2)/(R2+2.25)
              R =X+T(I)
              R2=R**2
              D =1./(R2+Y2)
              D3=Y1*D
              D4=R*D
              WR=WR+Y*(U(I)*(R*D4-1.5*D3)-S(I)*Y3*D4)/(R2+2.25)
              WI=WI+U(I)*(D2+D4)+S(I)*(D1-D3)
          end do  
      END IF
  END subroutine CPF  

  Subroutine CPF3(X,Y,WR,WI)
!-------------------------------------------------
! "CPF": Complex Probability Function
! .........................................................
!         .       Subroutine to Compute the Complex       .
!         .        Probability Function W(z=X+iY)         .
!         .     W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0      .
!         .    Which Appears when Convoluting a Complex   .
!         .     Lorentzian Profile by a Gaussian Shape    .
!         .................................................
!
!             WR : Real Part of W(z)
!             WI : Imaginary Part of W(z)
!
! This Routine takes into account the region 3 only, i.e. when sqrt(x**2+y**2)>8. 
!-------------------------------------------------
!   
    real, intent(in)  :: X, Y
    real, intent(out) :: WR,WI

    Integer :: I
    complex :: zm1,zm2,zterm,zsum

! Region 3
    zm1  = zone/cmplx(x,y)
    zm2  = zm1**2
    ! zsum = zone
    ! zterm= zone
    ! do i=1,tt_len
    !    zterm = zterm*zm2*tt(i)
    !    zsum  = zsum+zterm
    ! end do
    ! zsum = zsum*zi*zm1*pipwoeronehalf

    zsum = zi*zm1*pipwoeronehalf
    zterm= zsum
    do i=1,tt_len
       zterm = zterm*zm2*tt(i)
       zsum  = zsum+zterm
    end do
    wr   = real(zsum)
    wi   = imag(zsum)
    return
    End Subroutine CPF3

end module pCqSDHC_module


! ! code to compare pCqSDHC with
! ! speed dependent Voigt

! program test
!     use pCqSDHC_module
!     use voigt_module
!     implicit none
    
!     real :: VU0 = 1000.
!     real :: GammaD = 0.005
!     real :: Gamma0 = 0.52
!     real :: Gamma2
!     real :: Delta0 = 0.002
!     ! real :: Delta2 = 0.001 * Delta0
!     ! real :: nuVC = 0.2
!     ! real :: eta = 0.5
!     real :: DVU = 1.
!     real :: VU
!     real :: nuVC = 0.
!     real :: eta = 0.
!     real :: Delta2 = 0.
    
!     real :: ALPHAL
!     real :: ALPHAD
!     real :: SDEP
!     real :: LS_pCqSDHC_R
!     real :: LS_pCqSDHC_I
!     real :: SDV
    
    
!     integer :: k
    
!     Gamma2 = 0.01 * Gamma0
    
!     ALPHAL = Gamma0
!     ALPHAD = GammaD
!     SDEP   = Gamma2/ALPHAL
!     Delta0 = 0.
    
!     VU = VU0-DVU
!     do while(VU < VU0+DVU)
!         call pCqSDHC(&
!                    VU0,            &!<--: Unperturbed line position in cm-1 (Input).
!                    GammaD,         &!<--: Doppler HWHM in cm-1 (Input)
!                    Gamma0,         &!<--: Speed-averaged line-width in cm-1 (Input).    
!                    Gamma2,         &!<--: Speed dependence of the line-width in cm-1 (Input).
!                    nuVC,           &!<--: Velocity-changing frequency in cm-1 (Input).
!                    eta,            &!<--: Correlation parameter, No unit (Input).
!                    Delta0,         &!<--: Speed-averaged line-shift in cm-1 (Input).
!                    Delta2,         &!<--: Speed dependence of the line-shift in cm-1 (Input)     
!                    VU,             &!<--: Current WaveNumber of the Computation in cm-1 (Input).
!                    LS_pCqSDHC_R,   &!-->: Real part of the normalized spectral shape (cm) (Output)
!                    LS_pCqSDHC_I    &!-->: Imaginary part of the normalized spectral shape (cm) (Output)
!                    )
!         SDV=SDVOIGT(real(VU-VU0, kind=8),ALPHAL,ALPHAD,SDEP)
!         print *, VU, LS_pCqSDHC_R, SDV, (LS_pCqSDHC_R/SDV-1)*100.0
!         VU=VU+0.01
!     end do
    
! end program test