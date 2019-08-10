!----------------------------------------------------------------------------
! This file is part of UMDLES.
!
! UMDLES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UMDLES is distributsed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2015-2017, David A. New
!----------------------------------------------------------------------------

module modstat_io
  implicit none

  integer, parameter :: ntypes = 8

  character(*), parameter :: time_name = 'time'
  character(*), parameter :: rank_name = 'rank'
  character(*), parameter :: zt_name   = 'zt'
  character(*), parameter :: zm_name   = 'zm'
  character(*), parameter :: type_name = 'type'

  interface create_stat_var
     module procedure create_stat_var_0d
     module procedure create_stat_var_1d
  end interface create_stat_var

  interface write_stat_var
     module procedure write_stat_real_var_0d
     module procedure write_stat_real_var_1d
     module procedure write_stat_real_var_2d
     module procedure write_stat_int_var_1d
  end interface

contains

  !
  ! create_stat_file
  !
  subroutine create_stat_file(name, ncid, irec)
    use netcdf_interface, only : create_file, create_dim, create_unlimited_dim, create_var, &
                                 my_nf90_real, my_nf90_int, write_var
    use mpi_interface, only    : pecount 
    use grid, only             : zt, zm, dn0, dthldtls, dqtdtls, level, wfls, pi0, u0, v0, nzp
    use defs, only             : cp, R, alvl
    use step, only             : cntlat
    implicit none

    character(*), intent(in) :: name
    integer, intent(out)     :: ncid, irec

    integer              :: i
    integer, allocatable :: rank_values(:), type_values(:)
    real, allocatable    :: pi0_by_cp(:)
    real                 :: cp_copy, alvl_copy, R_copy, cntlat_copy

    ! Create file
    call create_file(name, ncid)

    ! Create time dimension (unlimited)
    call create_unlimited_dim(ncid, time_name, my_nf90_real)

    ! Create rank dimension if using MPI
    if (pecount > 1) then
       allocate(rank_values(pecount))
       do i = 1,pecount
          rank_values(i) = i
       end do
       call create_dim(ncid, rank_name, rank_values)
    end if

    ! Define altitude dimensions
    call create_dim(ncid, zt_name, zt)
    call create_dim(ncid, zm_name, zm)
    
    !
    allocate(type_values(ntypes))
    do i = 1,ntypes
       type_values(i) = i
    end do
    call create_dim(ncid, type_name, type_values)

    ! Create non-time-dependant parameters for simulation parameters
    call create_var(ncid=ncid, var_name='cp',   type=my_nf90_real)
    call create_var(ncid=ncid, var_name='alvl', type=my_nf90_real)
    call create_var(ncid=ncid, var_name='R',    type=my_nf90_real)
    call create_var(ncid=ncid, var_name='cntlat', type=my_nf90_real)

    call create_var(ncid, 'dn0',    my_nf90_real, (/zt_name/))
    call create_var(ncid, 'pi0',    my_nf90_real, (/zt_name/))

    ! Create variable for first and last sample times
    call create_var(ncid, 'nsmp',  my_nf90_int, (/time_name/))
    call create_var(ncid, 'fsttm', my_nf90_real, (/time_name/))
    call create_var(ncid, 'lsttm', my_nf90_real, (/time_name/))

    ! Needed for intent(inout) with parameter
    cp_copy     = cp
    alvl_copy   = alvl
    R_copy      = R
    cntlat_copy = cntlat
    allocate(pi0_by_cp(nzp))
    pi0_by_cp(:) = pi0(:)/cp

    ! Write parameters
    call write_var(ncid, 'cp',   cp_copy)
    call write_var(ncid, 'alvl', alvl_copy)
    call write_var(ncid, 'R',    R_copy)
    call write_var(ncid, 'cntlat', cntlat)

    call write_var(ncid, 'dn0', dn0)
    call write_var(ncid, 'pi0', pi0_by_cp)

    ! Initialize record index
    irec = 0

  end subroutine create_stat_file

  !
  ! create_stat_var_0d
  !
  subroutine create_stat_var_0d(ncid, var_name, data_type)
    use netcdf_interface, only : create_var, create_dim
    use mpi_interface, only    : pecount
    implicit none

    integer, intent(in)      :: ncid, data_type
    character(*), intent(in) :: var_name

    if (pecount == 1) then
       call create_var(ncid, var_name, data_type, (/time_name/))
    else
       call create_var(ncid, var_name, data_type, (/rank_name, time_name/))
    end if

  end subroutine create_stat_var_0d

  !
  ! create_stat_var_1d
  !
  subroutine create_stat_var_1d(ncid, var_name, z_dim_name, data_type)
    use netcdf_interface, only : create_var, create_dim
    use mpi_interface, only    : pecount
    implicit none

    integer, intent(in)      :: ncid, data_type
    character(*), intent(in) :: var_name, z_dim_name

    if (pecount == 1) then
       call create_var(ncid, var_name, data_type, (/z_dim_name, time_name/))
    else
       call create_var(ncid, var_name, data_type, (/z_dim_name, rank_name, time_name/))
    end if

  end subroutine create_stat_var_1d

  !
  ! create_stat_var_2d
  !
  subroutine create_stat_var_2d(ncid, var_name, z_dim_name, data_type)
    use netcdf_interface, only : create_var, create_dim
    use mpi_interface, only    : pecount
    implicit none

    integer, intent(in)      :: ncid, data_type
    character(*), intent(in) :: var_name, z_dim_name

    if (pecount == 1) then
       call create_var(ncid, var_name, data_type, (/z_dim_name, type_name, time_name/))
    else
       call create_var(ncid, var_name, data_type, (/z_dim_name, type_name, rank_name, time_name/))
    end if

  end subroutine create_stat_var_2d

  !
  ! advance_stat_time
  !
  subroutine advance_stat_time(time, fsttm, nsmp, ncid, irec)
    use netcdf_interface, only : write_var
    use grid, only             : pi1, nzp
    use defs, only             : cp
    implicit none

    real, intent(in)     :: time, fsttm
    integer, intent(in)  :: nsmp, ncid
    integer, intent(out) :: irec

    real    :: time_copy, fsttm_copy, stat_time
    integer :: nsmp_copy

    ! Because Parallel NetCDF does not allow intent(in) variables
    time_copy  = time
    nsmp_copy  = nsmp
    fsttm_copy = fsttm

    !
    stat_time = 0.5*(fsttm + time)

    ! Increment record counter
    irec = irec + 1

    ! Write to disk
    call write_var(ncid, time_name, stat_time, irec)
    call write_var(ncid, 'nsmp',  nsmp_copy,   irec)
    call write_var(ncid, 'fsttm', fsttm_copy,  irec)
    call write_var(ncid, 'lsttm', time_copy,   irec)

  end subroutine advance_stat_time

  !
  ! write_stat_real_var_0d
  !
  subroutine write_stat_real_var_0d(ncid, name, var, irec)
    use mpi_interface, only    : pecount, myid
    use grid, only             : nzp
    use netcdf_interface, only : write_var
    implicit none
    
    integer, intent(in)      :: ncid, irec
    character(*), intent(in) :: name
    real, intent(in)         :: var

    real :: var_copy(1)

    ! Because Parallel NetCDF does not allow intent(in) variables
    var_copy(1) = var

    if (pecount == 1) then
       call write_var(ncid, name, var_copy, irec)
    else
       call write_var(ncid, name, var_copy, irec, start_in=(/myid+1/), stride_in=(/1/))
    end if

  end subroutine write_stat_real_var_0d

  !
  ! write_stat_real_var_1d
  !
  subroutine write_stat_real_var_1d(ncid, name, var, irec)
    use mpi_interface, only    : pecount, myid
    use grid, only             : nzp
    use netcdf_interface, only : write_var
    implicit none
    
    integer, intent(in)      :: ncid, irec
    character(*), intent(in) :: name
    real, intent(in)         :: var(:)

    real, allocatable :: var_copy(:)

    ! Because Parallel NetCDF does not allow intent(in) variables
    allocate(var_copy(nzp))
    var_copy(:) = var(:)

    if (pecount == 1) then
       call write_var(ncid, name, var_copy, irec)
    else
       call write_var(ncid, name, var_copy, irec, start_in=(/1, myid+1/), stride_in=(/nzp, 1/))
    end if

  end subroutine write_stat_real_var_1d

  !
  ! write_stat_real_var_2d
  !
  subroutine write_stat_real_var_2d(ncid, name, var, irec)
    use mpi_interface, only    : pecount, myid
    use grid, only             : nzp
    use netcdf_interface, only : write_var
    implicit none
    
    integer, intent(in)      :: ncid, irec
    character(*), intent(in) :: name
    real, intent(in)         :: var(:,:)

    real, allocatable :: var_copy(:,:)

    ! Because Parallel NetCDF does not allow intent(in) variables
    allocate(var_copy(nzp,ntypes))
    var_copy(:,:) = var(:,:)

    if (pecount == 1) then
       call write_var(ncid, name, var_copy, irec)
    else
       call write_var(ncid, name, var_copy, irec, start_in=(/1, 1, myid+1/), stride_in=(/nzp, ntypes, 1/))
    end if

  end subroutine write_stat_real_var_2d

  !
  ! write_stat_int_var_1d
  !
  subroutine write_stat_int_var_1d(ncid, name, var, irec)
    use mpi_interface, only    : pecount, myid
    use grid, only             : nzp
    use netcdf_interface, only : write_var
    implicit none
    
    integer, intent(in)      :: ncid, irec
    character(*), intent(in) :: name
    integer, intent(in)      :: var(:)

    integer, allocatable :: var_copy(:)

    ! Because Parallel NetCDF does not allow intent(in) variables
    allocate(var_copy(nzp))
    var_copy(:) = var(:)

    if (pecount == 1) then
       call write_var(ncid, name, var_copy, irec)
    else
       call write_var(ncid, name, var_copy, irec, start_in=(/1, myid+1/), stride_in=(/nzp, 1/))
    end if

  end subroutine write_stat_int_var_1d
  
end module modstat_io




