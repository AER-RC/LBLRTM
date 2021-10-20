!<f90File>**************************************************************
!
! CONTACT:
!
!   Atmospheric & Environmental Research, Inc
!   131 Hartwell Ave
!   Lexington ,MA 02421-3126 USA
!   Phone: 781.761.2288
!   E-mail: ipolonsk@aer.com
!
! COPYRIGHT NOTICE:
!
!   Copyright AER, Inc 2001-2019, All Right Reserved
!   See the file README-DATARIGHTS.txt included with this release
!   for additional details.
!
!*************************************************************</f90File>


! Calculates wv continuum coefficients: reads from netCDF file: wv-mt-ckd.nc

Module h2oCont
   
     USE read_file
!                                                                       
      IMPLICIT none
      private 
      public h2oAbsrb
!              
      type(data2read)  :: dat
      character(len=*), parameter :: fDataname = "wv-mt-ckd.nc"

      real,parameter :: p0=1013.,t0 =296.,xlosmt=2.68675E+19
      real ::onepl = 1.001,onemi = 0.999 

   contains
   subroutine h2oAbsrb(pave,tave,xself,xfrgn,v1abs,v2abs,dvabs,wv_self_abs,wv_for_abs, &
                    radflag,radcn2) 

! Inputs
   real,dimension(:),intent(inout)  :: wv_self_abs,wv_for_abs 
   real, intent(in) :: pave,tave,dvabs,xself,xfrgn
   real, intent(in):: v1abs,v2abs
   logical,optional :: radflag
   real, optional :: radcn2

! Local variables
   integer :: ncoeff,nptabs,i1,i2,ist,lst,i
   real :: tfac,xkt,dvc

   integer,save :: ncoeffin
   real,dimension(:), allocatable,save :: wvn
   real,dimension(:,:), allocatable,save :: coeff

   real,dimension(:), allocatable :: sh2o,fh2o,rad
   integer :: iret
   logical,save :: lread=.False.
   
!#ifdef USENETCDF
   if (.not. lread) then
      lread = .True.
      if (getData(fDataname,dat)) STOP   
      ncoeffin = size(dat%wavenumber)
      if (allocated(wvn)) deallocate(wvn)
      if (allocated(coeff)) deallocate(coeff)
      allocate (wvn(ncoeffin))
      allocate (coeff(ncoeffin,3))
      wvn = dat%wavenumber(:)
      coeff(:,1) = dat%wvself_260(:)
      coeff(:,2) = dat%wvself_296(:)
      coeff(:,3) = dat%wvforeign(:)
   endif
!#endif

   if (.not. present(radflag)) then
      radflag = .TRUE.
   endif
   if(.not. present(radcn2)) then
     radcn2 = 1.4387752
   endif

! Find coeff wavenumber range that brackets [v1abs,v2abs]
   dvc = dat%wavenumber(2)-dat%wavenumber(1)
   i=1
   do while (wvn(i) < (v1abs-2*dvc))
      i = i+1
   enddo
   i1=i-1 
   do while (wvn(i) < (v2abs+2*dvc))
      i = i+1
   enddo
   i2=i
   ncoeff = i2-i1+1

! Set up arrays 
   
   if (allocated (sh2o)) deallocate (sh2o)
   if (allocated (fh2o)) deallocate (fh2o)
   if (allocated (rad)) deallocate (rad)
   allocate (sh2o(ncoeff))
   allocate (fh2o(ncoeff))
   allocate (rad(ncoeff))
   
! Define some atmospheric parameters
   xkt = tave/radcn2 

! Apply temperature dependence for water vapor self and scale
    tfac = (tave-t0)/(260.-t0)
    sh2o = coeff(i1:i2,2)*(coeff(i1:i2,1)/coeff(i1:i2,2))**tfac
    sh2o = sh2o*xself

! Multiply by radiation term if requested
    if (radflag) then 
       iret = myradfn(wvn(i1:i2),xkt,ncoeff,rad)
       sh2o = sh2o*rad
    endif

! Interpolate water vapor self to output grid
   nptabs = (v2abs-v1abs)/dvabs+1
   call pre_xint(wvn(1),wvn(ncoeffin),v1abs,dvabs,nptabs,ist,lst)
   call myxint(wvn(i1),wvn(i2),dvc,sh2o,1.0,v1abs,dvabs,wv_self_abs,ist,lst)

! Generate  water vapor foreign
  fh2o(:) = coeff(i1:i2,3)*xfrgn
   
! Multiply by radiation term if requested
    if (radflag) then 
       fh2o = fh2o*rad
    endif

! Interpolate water vapor foreign to output grid
   call myxint(wvn(i1),wvn(i2),dvc,fh2o,1.0,v1abs,dvabs,wv_for_abs,ist,lst)


   end subroutine h2oAbsrb


   integer function myradfn(vi,xkt,nvi,rad)

! Input variables

   real, dimension(:),intent(in) :: vi
   real, intent(in) :: xkt
   integer,intent(in) :: nvi

! Output variable
  real, dimension(:),intent(out) :: rad


! Local variables
   real :: xvi(nvi), xviokt(nvi),expvkt(nvi)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     FUNCTION RADFN CALCULATES THE RADIATION TERM FOR THE LINE SHAPE
!
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.
!
!
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL
!
!                                             FASCOD3
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!      IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED
!
   xvi = vi
   rad = xvi
   xviokt = xvi/xkt
!
   where (xviokt.le.0.01) 
      rad = 0.5*xviokt*xvi
   elsewhere (xviokt.le.10.)
      expvkt = exp(-xviokt)
      rad = xvi*(1-expvkt)/(1.+expvkt)
   endwhere

   myradfn=1
!
   RETURN
!
   end function myradfn


   subroutine myxint (v1a,v2a,dva,a,afact,vft,dvr3,r3,n1r3,n2r3)
!
!
!     THIS SUBROUTINE INTERPOLATES THE A ARRAY STORED
!     FROM V1A TO V2A IN INCREMENTS OF DVA USING A MULTIPLICATIVE
!     FACTOR AFACT, INTO THE R3 ARRAY FROM LOCATION N1R3 TO N2R3 IN
!     INCREMENTS OF DVR3
!
! Input variables
   real,dimension(:),intent(in) :: a
    real,intent(in) :: v1a,v2a,vft
   real, intent(in) :: dva,afact,dvr3
   integer, intent(in) :: n1r3,n2r3
 
! Input/output variables
   real,dimension(:),intent(inout) :: r3

! Local variables
   real :: recdva,p,c,b,b1,b2,conti
   integer :: i,j
   integer :: ilo,ihi
   real :: vi,vj

!
   RECDVA = 1./DVA
   ILO = (V1A+DVA-VFT)/DVR3+1.+ONEMI
   ILO = MAX(ILO,N1R3)
   IHI = (V2A-DVA-VFT)/DVR3+ONEMI
   IHI = MIN(IHI,N2R3)
!
   DO 10 I = ILO, IHI
      VI = VFT+DVR3* REAL(I-1)
      J = (VI-V1A)*RECDVA+ONEPL
      VJ = V1A+DVA* REAL(J-1)
      P = RECDVA*(VI-VJ)
      C = (3.-2.*P)*P*P
      B = 0.5*P*(1.-P)
      B1 = B*(1.-P)
      B2 = B*P 
      CONTI = -A(J-1)*B1+A(J)*(1.-C+B2)+A(J+1)*(C+B1)-A(J+2)*B2
      R3(I) = R3(I)+CONTI*AFACT
10 END DO
!
   RETURN
!
   end subroutine myxint

   subroutine pre_xint(v1ss,v2ss,v1abs,dvabs,nptabs,ist,lst)


! Input variables
    real,intent(in) :: v1ss,v2ss,v1abs
   integer, intent(in) :: nptabs
   real, intent(in) :: dvabs

! Output variables
   integer, intent(out) :: ist,lst

! Local variables
   integer :: nbnd_v1c,nbnd_v2c
   real :: v1abs_loc

!   Set up needed variables for call to XINT
!   Output variables
!     ist - index of first value to be processed in XINT
!     lst - index of last value to be processed in XINT

   nbnd_v1c =  2 +  (v1ss-v1abs)/dvabs + 1.e-5
   ist = max(1,nbnd_v1c)
   v1abs_loc = v1abs + dvabs * float(ist-1)

   nbnd_v2c = 1 + (v2ss-v1abs)/dvabs + 1.e-5
   lst = min(nptabs,nbnd_v2c)

   return
  
   end subroutine pre_xint

end module h2oCont
