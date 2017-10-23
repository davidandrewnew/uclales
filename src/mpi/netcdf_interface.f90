!----------------------------------------------------------------------------
! This file is part of UMDLES.
!
! UMDLES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UMDLES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2015-2016, David A. New
!----------------------------------------------------------------------------
!

module netcdf_interface
  implicit none

  real    :: fill_value_real
  integer :: my_nf90_real, my_nf90_int

  interface write_var
     module procedure write_real_0d
     module procedure write_int_0d
     module procedure write_real_1d
     module procedure write_real_2d
     module procedure write_int_1d
     module procedure write_byte_1d
  end interface

  interface create_dim
     module procedure create_real_dim
     module procedure create_int_dim
  end interface

  contains

    !
    ! init_netcdf_interface
    !
    subroutine init_netcdf_interface
      use pnetcdf, only       : nf90_float, nf90_double, nf90_int
      use mpi_interface, only : appl_abort, myid
      implicit none

      real :: temp_real

      ! Create real fill values (NaN)
      temp_real       = 0.
      fill_value_real = 1./temp_real

      ! Get real type
      select case(kind(0.0))
      case(4)
         my_nf90_real = nf90_float
      case(8)
         my_nf90_real = nf90_double
      case default
         if (myid == 0) write(*,*) 'NetCDF real type not supported'
         call appl_abort(0)
      end select

      ! Get int type
      select case(kind(0))
      case(4)
         my_nf90_int = nf90_int
      case default
         if (myid == 0) write(*,*) 'NetCDF int type not supported'
         call appl_abort(0)
      end select

    end subroutine init_netcdf_interface

    !
    ! handle_error
    !
    subroutine handle_error(iret, message)
      use mpi_interface, only : appl_abort
      use pnetcdf, only       : nf90_noerr, nf90mpi_strerror

      implicit none

      integer, intent(in)          :: iret
      character(len=*), intent(in) :: message

      if (iret /= nf90_noerr) then
         write(6,*) trim(message), trim(nf90mpi_strerror(iret))
         call appl_abort(0)
      end if
    end subroutine handle_error

    !
    ! define_mode
    !
    subroutine define_mode(ncid)
      use pnetcdf, only : nf90mpi_redef, nf90_eindefine
      implicit none

      integer, intent(in) :: ncid

      integer :: iret

      ! Exit define mode
      iret = nf90mpi_redef(ncid)
      if (iret /= nf90_eindefine) call handle_error(iret, 'In nf90mpi_redef: ')

    end subroutine define_mode

    !
    ! independent_data_mode
    !
    subroutine independent_data_mode(ncid)
      use pnetcdf, only : nf90mpi_enddef, nf90_enotindefine, nf90_eindep, nf90mpi_begin_indep_data
      implicit none

      integer, intent(in) :: ncid

      integer :: iret

      iret = nf90mpi_enddef(ncid)
      if (iret == nf90_enotindefine) then
         iret = nf90mpi_begin_indep_data(ncid)
         if (iret /= nf90_eindep) call handle_error(iret, 'In nf90mpi_begin_indep_data: ')
      else
         call handle_error(iret, 'In nf90mpi_enddef: ')
      end if

    end subroutine independent_data_mode

    !
    ! collective_data_mode
    !
    subroutine collective_data_mode(ncid)
      use pnetcdf, only : nf90mpi_enddef, nf90_enotindefine, nf90_enotindep, nf90mpi_end_indep_data
      implicit none

      integer, intent(in) :: ncid

      integer :: iret

      iret = nf90mpi_enddef(ncid)
      if (iret == nf90_enotindefine) then
         iret = nf90mpi_end_indep_data(ncid)
         if (iret /= nf90_enotindep) call handle_error(iret, 'In nf90mpi_endd_indep_data: ')
      else
         call handle_error(iret, 'In nf90mpi_enddef: ')         
      end if

    end subroutine collective_data_mode

    !
    ! create_file
    !
    subroutine create_file(name, ncid)
      use pnetcdf, only : nf90_clobber, nf90_64bit_data, nf90mpi_create
      use mpi
      implicit none

      character(len=*), intent(in)   :: name
      integer, intent(out)           :: ncid

      integer :: iret, cmode

      ! Create file
      cmode = IOR(nf90_clobber, nf90_64bit_data)
      iret = nf90mpi_create(MPI_COMM_WORLD, trim(name), nf90_clobber, MPI_INFO_NULL, ncid)
      call handle_error(iret, 'In nf90mpi_create: ' // trim(name))

    end subroutine create_file

    !
    ! close_file
    !
    subroutine close_file(ncid)
      use pnetcdf, only : nf90mpi_close
      implicit none

      integer, intent(in) :: ncid

      integer :: iret

      ! Close file
      iret = nf90mpi_close(ncid)
      call handle_error(iret, 'In nf90mpi_close: ')

    end subroutine close_file

    !
    ! create_unlimited_dim
    !
    subroutine create_unlimited_dim(ncid, name, type)
      use pnetcdf, only : nf90mpi_def_dim, nf90mpi_unlimited, nf90_enameinuse
      implicit none

      integer, intent(in)      :: ncid, type
      character(*), intent(in) :: name

      integer :: iret, dimid

      !
      call define_mode(ncid)

      ! Create dimension
      iret = nf90mpi_def_dim(ncid, trim(name), nf90mpi_unlimited, dimid)
      if (iret /= nf90_enameinuse) then
         call handle_error(iret, 'In nf90mpi_def_dim: ' // trim(name))

         ! Create corresponding variable
         call create_var(ncid, trim(name), (/trim(name)/), type)
      end if

    end subroutine create_unlimited_dim

    !
    ! create_real_dim
    !
    subroutine create_real_dim(ncid, name, values)
      use pnetcdf, only : nf90mpi_def_dim, nf90mpi_unlimited, nf90_enameinuse
      use mpi
      implicit none

      integer, intent(in)           :: ncid
      character(*), intent(in)      :: name
      real, intent(inout), optional :: values(:)

      integer :: iret, dimid

      !
      call define_mode(ncid)

      ! Create dimension
      iret = nf90mpi_def_dim(ncid, trim(name), int(size(values), MPI_OFFSET_KIND), dimid)
      if (iret /= nf90_enameinuse) then
         call handle_error(iret, 'In nf90mpi_def_dim: ' // trim(name))

         ! Create corresponding variable
         call create_var(ncid, trim(name), (/trim(name)/), my_nf90_real)
         call write_var(ncid, trim(name), values)
      end if

    end subroutine create_real_dim

    !
    ! create_int_dim
    !
    subroutine create_int_dim(ncid, name, values)
      use pnetcdf, only : nf90mpi_def_dim, nf90_enameinuse
      use mpi
      implicit none

      integer, intent(in)      :: ncid
      character(*), intent(in) :: name
      integer, intent(inout)   :: values(:)

      integer :: iret, dimid

      ! 
      call define_mode(ncid)
      
      ! Create dimension
      iret = nf90mpi_def_dim(ncid, trim(name), int(size(values), MPI_OFFSET_KIND), dimid)
      if (iret /= nf90_enameinuse) then
         call handle_error(iret, 'In nf90mpi_def_dim: ' // trim(name))

         ! Create corresponding variable
         call create_var(ncid, trim(name), (/trim(name)/), my_nf90_int)
         call write_var(ncid, trim(name), values)
      end if

    end subroutine create_int_dim

    !
    ! create_var
    !
    subroutine create_var(ncid, var_name, dim_names, type, varid)
      use pnetcdf, only : nf90mpi_def_var, nf90mpi_inq_dimid
      implicit none

      integer, intent(in)            :: ncid, type
      character(*), intent(in)       :: var_name
      character(*), dimension(:)     :: dim_names
      integer, intent(out), optional :: varid

      integer              :: iret, ndims, idim, varid_temp
      integer, allocatable :: dimids(:)

      ! Ensure define mode
      call define_mode(ncid)

      ! Get number of dimensions
      ndims = size(dim_names)
      allocate(dimids(ndims))
      
      ! Get dimension IDs
      do idim = 1,ndims
         iret = nf90mpi_inq_dimid(ncid, trim(dim_names(idim)), dimids(idim))
         call handle_error(iret, 'In nf90mpi_inq_dimid: ' // trim(dim_names(idim)))
      end do

      ! Define variable
      iret = nf90mpi_def_var(ncid, trim(var_name), type, dimids, varid_temp)
      call handle_error(iret, 'In nf90mpi_def_var: ' // trim(var_name))

      ! Output varid
      if (present(varid)) varid = varid_temp

    end subroutine create_var

    !
    ! write_real_0d
    !
    subroutine write_real_0d(ncid, name, value, irec)
      use pnetcdf, only : nf90mpi_inq_varid, nf90mpi_iput_var, nf90mpi_put_var_all
      use mpi
      implicit none
      
      integer, intent(in)            :: ncid
      character(*), intent(in)       :: name
      real, intent(inout)            :: value
      integer, intent(in), optional  :: irec

      integer                       :: iret, varid
      integer(kind=MPI_OFFSET_KIND) :: start(1)

      ! Set start
      if (present(irec)) then
         start(1) = irec
      else
         start(1) = 1
      end if

      ! Get variable ID
      iret = nf90mpi_inq_varid(ncid, trim(name), varid)
      call handle_error(iret, 'In nf90mpi_inq_varid: ' // trim(name))

      ! Ensure collective data mode
      call collective_data_mode(ncid)

      ! Write data
      iret = nf90mpi_put_var_all(ncid, varid, value, start)
      call handle_error(iret, 'In nf90mpi_put_var: ' // trim(name))

    end subroutine write_real_0d

    !
    ! write_int_0d
    !
    subroutine write_int_0d(ncid, name, value, irec)
      use pnetcdf, only : nf90mpi_inq_varid, nf90mpi_iput_var, nf90mpi_put_var_all
      use mpi
      implicit none
      
      integer, intent(in)            :: ncid
      character(*), intent(in)       :: name
      integer, intent(inout)         :: value
      integer, intent(in), optional  :: irec

      integer                       :: iret, varid
      integer(kind=MPI_OFFSET_KIND) :: start(1)

      ! Set start
      if (present(irec)) then
         start(1) = irec
      else
         start(1) = 1
      end if

      ! Get variable ID
      iret = nf90mpi_inq_varid(ncid, trim(name), varid)
      call handle_error(iret, 'In nf90mpi_inq_varid: ' // trim(name))

      ! Ensure collective data mode
      call collective_data_mode(ncid)

      ! Write data
      iret = nf90mpi_put_var_all(ncid, varid, value, start)
      call handle_error(iret, 'In nf90mpi_put_var: ' // trim(name))

    end subroutine write_int_0d

    !
    ! write_real_1d
    !
    subroutine write_real_1d(ncid, name, values, irec, start_in, stride_in)
      use mpi
      use pnetcdf, only : nf90mpi_iput_var, nf90mpi_inq_varid, nf90mpi_put_var_all
      implicit none

      integer, intent(in)            :: ncid
      character(*), intent(in)       :: name
      real, intent(inout)            :: values(:)
      integer, intent(in), optional  :: irec, start_in(:), stride_in(:)

      integer                                    :: iret, varid, ndims, ndims_loop
      integer(kind=MPI_OFFSET_KIND), allocatable :: start(:), stride(:)
      
      ! Get variable ID
      iret = nf90mpi_inq_varid(ncid, trim(name), varid)
      call handle_error(iret, 'In nf90mpi_inq_varid: ' // trim(name))

      ! Get number of dimensions 
      if (present(start_in)) then
         ndims = size(start_in)
      else 
         ndims = size(shape(values))
      end if
      if (present(irec)) then
         ndims_loop = ndims
         ndims      = ndims + 1
      else
         ndims_loop = ndims
      end if

      ! Allocate arrays
      allocate(start(ndims))

      ! Get starts and strides
      if (present(start_in)) then
         start(1:ndims_loop) = int(start_in(1:ndims_loop), MPI_OFFSET_KIND)

         if (present(stride_in)) then
            allocate(stride(ndims))
            stride(1:ndims_loop) = int(stride_in(1:ndims_loop), MPI_OFFSET_KIND)
         end if
      else
         start(1:ndims_loop) = 1
      end if
      if (present(irec)) then
         start(ndims) = int(irec, MPI_OFFSET_KIND)
         if (present(start_in) .and. present(stride_in)) stride(ndims) = 1
      end if

      ! Ensure collective data mode
      call collective_data_mode(ncid)

      ! Write variable
      if (present(stride_in)) then
         iret = nf90mpi_put_var_all(ncid, varid, values, start, stride)
         call handle_error(iret, 'In nf90mpi_put_var: ' // trim(name))
      else
         iret = nf90mpi_put_var_all(ncid, varid, values, start)
         call handle_error(iret, 'In nf90mpi_put_var: ' // trim(name))
      end if

    end subroutine write_real_1d

    !
    ! write_real_2d
    !
    subroutine write_real_2d(ncid, name, values, irec, start_in, stride_in)
      use mpi
      use pnetcdf, only : nf90mpi_iput_var, nf90mpi_inq_varid, nf90mpi_put_var_all
      implicit none

      integer, intent(in)            :: ncid
      character(*), intent(in)       :: name
      real, intent(inout)            :: values(:,:)
      integer, intent(in), optional  :: irec, start_in(:), stride_in(:)

      integer                                    :: iret, varid, ndims, ndims_loop
      integer(kind=MPI_OFFSET_KIND), allocatable :: start(:), stride(:)
      
      ! Get variable ID
      iret = nf90mpi_inq_varid(ncid, trim(name), varid)
      call handle_error(iret, 'In nf90mpi_inq_varid: ' // trim(name))

      ! Get number of dimensions 
      if (present(start_in)) then
         ndims = size(start_in)
      else 
         ndims = size(shape(values))
      end if
      if (present(irec)) then
         ndims_loop = ndims
         ndims      = ndims + 1
      else
         ndims_loop = ndims
      end if

      ! Allocate arrays
      allocate(start(ndims))

      ! Get starts and strides
      if (present(start_in)) then
         start(1:ndims_loop) = int(start_in(1:ndims_loop), MPI_OFFSET_KIND)

         if (present(stride_in)) then
            allocate(stride(ndims))
            stride(1:ndims_loop) = int(stride_in(1:ndims_loop), MPI_OFFSET_KIND)
         end if
      else
         start(1:ndims_loop) = 1
      end if
      if (present(irec)) then
         start(ndims) = int(irec, MPI_OFFSET_KIND)
         if (present(start_in) .and. present(stride_in)) stride(ndims) = 1
      end if

      ! Ensure collective data mode
      call collective_data_mode(ncid)

      ! Write variable
      if (present(stride_in)) then
         iret = nf90mpi_put_var_all(ncid, varid, values, start, stride)
         call handle_error(iret, 'In nf90mpi_put_var: ' // trim(name))
      else
         iret = nf90mpi_put_var_all(ncid, varid, values, start)
         call handle_error(iret, 'In nf90mpi_put_var: ' // trim(name))
      end if

    end subroutine write_real_2d

    !
    ! write_int_1d
    !
    subroutine write_int_1d(ncid, name, values, irec, start_in, stride_in)
      use mpi
      use pnetcdf, only : nf90mpi_iput_var, nf90mpi_inq_varid, nf90mpi_put_var_all
      implicit none

      integer, intent(in)            :: ncid
      character(*), intent(in)       :: name
      integer, intent(inout)         :: values(:)
      integer, intent(in), optional  :: irec, start_in(:), stride_in(:)

      integer                                    :: iret, varid, ndims, ndims_loop
      integer(kind=MPI_OFFSET_KIND), allocatable :: start(:), stride(:)
      
      ! Get variable ID
      iret = nf90mpi_inq_varid(ncid, trim(name), varid)
      call handle_error(iret, 'In nf90mpi_inq_varid: ' // trim(name))

      ! Get number of dimensions 
      if (present(start_in)) then
         ndims = size(start_in)
      else 
         ndims = size(shape(values))
      end if
      if (present(irec)) then
         ndims_loop = ndims
         ndims      = ndims + 1
      else
         ndims_loop = ndims
      end if

      ! Allocate arrays
      allocate(start(ndims))

      ! Get starts and strides
      if (present(start_in)) then
         start(1:ndims_loop) = int(start_in(1:ndims_loop), MPI_OFFSET_KIND)

         if (present(stride_in)) then
            allocate(stride(ndims))
            stride(1:ndims_loop) = int(stride_in(1:ndims_loop), MPI_OFFSET_KIND)
         end if
      else
         start(1:ndims_loop) = 1
      end if
      if (present(irec)) then
         start(ndims) = int(irec, MPI_OFFSET_KIND)
         if (present(start_in) .and. present(stride_in)) stride(ndims) = 1
      end if

      ! Ensure collective data mode
      call collective_data_mode(ncid)

      ! Write variable
      if (present(stride_in)) then
         iret = nf90mpi_put_var_all(ncid, varid, values, start, stride)
         call handle_error(iret, 'In nf90mpi_put_var: ' // trim(name))
      else
         iret = nf90mpi_put_var_all(ncid, varid, values, start)
         call handle_error(iret, 'In nf90mpi_put_var: ' // trim(name))
      end if

    end subroutine write_int_1d

    !
    ! write_byte_1d
    !
    subroutine write_byte_1d(ncid, name, values, irec, start_in, stride_in)
      use mpi
      use pnetcdf, only : nf90mpi_iput_var, nf90mpi_inq_varid, nf90mpi_put_var_all
      implicit none

      integer, intent(in)            :: ncid
      character(*), intent(in)       :: name
      integer(kind=1), intent(inout) :: values(:)
      integer, intent(in), optional  :: irec, start_in(:), stride_in(:)

      integer                                    :: iret, varid, ndims, ndims_loop
      integer(kind=MPI_OFFSET_KIND), allocatable :: start(:), stride(:)
      
      ! Get variable ID
      iret = nf90mpi_inq_varid(ncid, trim(name), varid)
      call handle_error(iret, 'In nf90mpi_inq_varid: ' // trim(name))

      ! Get number of dimensions 
      if (present(start_in)) then
         ndims = size(start_in)
      else 
         ndims = size(shape(values))
      end if
      if (present(irec)) then
         ndims_loop = ndims
         ndims      = ndims + 1
      else
         ndims_loop = ndims
      end if

      ! Allocate arrays
      allocate(start(ndims))

      ! Get starts and strides
      if (present(start_in)) then
         start(1:ndims_loop) = int(start_in(1:ndims_loop), MPI_OFFSET_KIND)

         if (present(stride_in)) then
            allocate(stride(ndims))
            stride(1:ndims_loop) = int(stride_in(1:ndims_loop), MPI_OFFSET_KIND)
         end if
      else
         start(1:ndims_loop) = 1
      end if
      if (present(irec)) then
         start(ndims) = int(irec, MPI_OFFSET_KIND)
         if (present(start_in) .and. present(stride_in)) stride(ndims) = 1
      end if

      ! Ensure collective data mode
      call collective_data_mode(ncid)

      ! Write variable
      if (present(stride_in)) then
         iret = nf90mpi_put_var_all(ncid, varid, values, start, stride)
         call handle_error(iret, 'In nf90mpi_put_var: ' // trim(name))
      else
         iret = nf90mpi_put_var_all(ncid, varid, values, start)
         call handle_error(iret, 'In nf90mpi_put_var: ' // trim(name))
      end if

    end subroutine write_byte_1d

    !
    ! Graveyard?
    !

    !
    ! sync
    !
    subroutine sync(ncid)
      use pnetcdf, only : nf90mpi_sync
      implicit none

      integer, intent(in) :: ncid

      integer :: iret

      call collective_data_mode(ncid)
      iret = nf90mpi_sync(ncid)
      call handle_error(iret, 'In nf90mpi_sync: ')

    end subroutine sync

  end module netcdf_interface






