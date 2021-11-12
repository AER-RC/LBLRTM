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

MODULE read_file
  implicit none
  include "netcdf.inc"
  !
  ! Module that provides FORTRAN 90 style netCDF reading support
  private
  public getData
  public data2read
  type data2read
     real(kind=8), allocatable, dimension(:) :: wavenumber
     real(kind=8), allocatable, dimension(:) :: wvforeign
     real(kind=8), allocatable, dimension(:) :: wvself_260
     real(kind=8), allocatable, dimension(:) :: wvself_296
  end type data2read

  logical, parameter                      :: dbg = .FALSE.
  interface readVarNC
     module procedure readReal1D
     module procedure readDouble1D
  end interface

  contains
  function getData(fname, dat) result(isError)
                      
    character(len=*), intent(in)     :: fname
    type(data2read),  intent(inout)  :: dat
    logical                          :: isError

    integer(kind=4)   :: ncid
    integer(kind=4)   :: nWavenumbers
    integer(kind=4)   :: stat

    ! check on the file
    inquire(file = fname, EXIST=isError) 
    isError = .NOT. isError
    if (isError) then
       print '("ERROR::read_file:: file not found ",A)', trim(fname)
       return
    endif

    if (dbg) print *, 'reading: ', trim(fname)
    call check( nf_open(fname, nf_nowrite, ncid) )

    if (.not. inqDim(ncid, "wavenumbers",  dimLen=nWavenumbers)) then
      call check( nf_close(ncid) )
      isError = .false.
      return
    end if
    ! allocate structure
    if (allocated(dat%wavenumber))   deallocate(dat%wavenumber)
    if (allocated(dat%wvforeign))    deallocate(dat%wvforeign)
    if (allocated(dat%wvself_260))   deallocate(dat%wvself_260)
    if (allocated(dat%wvself_296))   deallocate(dat%wvself_296)

    allocate(dat%wavenumber(nWavenumbers), &
             dat%wvforeign(nWavenumbers), &
             dat%wvself_260(nWavenumbers),   &
             dat%wvself_296(nWavenumbers), STAT= stat)
    isError = stat /= 0
    if (isError) then
       print '("ERROR::read_file:: memory allocation ")'
       return
    endif

    ! read variables
    call readVarNC(ncid,"wavenumber",   dat%wavenumber)
    call readVarNC(ncid,"wvforeign",    dat%wvforeign)
    call readVarNC(ncid,"wvself_260",   dat%wvself_260)
    call readVarNC(ncid,"wvself_296",   dat%wvself_296)
    call check( nf_close(ncid) )

  end function getData


   subroutine check(status, varName, fatal)
     integer,                     intent(in) :: status
     character(len=*),  optional, intent(in) :: varName
     logical,           optional, intent(in) :: fatal

     logical             :: fatalLoc
     if (present(fatal)) then
        fatalLoc = fatal
      else
        fatalLoc = .true.
      end if

     if(status /= nf_noerr) then
       if (present(varName)) print *,'Processing: ', varName
       if (fatalLoc) then
         print *, 'netCDF error: ', status, ' : ', trim(nf_strerror(status))
        call exit(1)
      else
         print *, 'netCDF WARNING: ', status, ' : ', trim(nf_strerror(status))
         print *, 'status', status
    end if
     end if
   end subroutine check

   function inqDim(id, dimName, dimLen)
      integer(kind=4),          intent(in)    :: id
      character(len=*),         intent(in)    :: dimName
      integer(kind=4), optional, intent(inout):: dimLen
      integer(kind=4)                         :: dimId
      logical                                 :: inqDim

      if (dbg) print*, ' ncdfUtil::inqDim '
      inqDim = (nf_noerr == nf_inq_dimid(id, dimName, dimId))
      if (present(dimLen)) then
        if (inqDim) then
          call check( NF_INQ_DIMLEN(id, dimId,dimLen) )
        else
          dimLen = -1
        end if
      end if
   end function inqDim  

   subroutine readReal1D(id, varName, val, fatal)
      integer(kind=4),        intent(in)   :: id
      character(len=*),       intent(in)   :: varName
      real*4, dimension(:),   intent(inout):: val
      integer(kind=4)                      :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readReal1D '
      call check(nf_inq_varid(id, varName, varId), varName, fatal)
      call check(nf_get_var(id, varId, val), varName, fatal)
   end subroutine readReal1D

   subroutine readDouble1D(id, varName, val, fatal)
      integer(kind=4),        intent(in)   :: id
      character(len=*),       intent(in)   :: varName
      real*8, dimension(:),   intent(inout):: val
      integer(kind=4)                      :: varId
      logical,        optional, intent(in) :: fatal
      if (dbg) print*, ' ncdfUtil::readDouble1D '
      call check(nf_inq_varid(id, varName, varId), varName, fatal)
      call check(nf_get_var(id, varId, val), varName, fatal)
   end subroutine readDouble1D

end module read_file
