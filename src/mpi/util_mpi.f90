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
module modutil_mpi
  implicit none

  ! par_max
  interface par_max
     module procedure real0D_par_max
     module procedure int0D_par_max
  end interface

  ! par_sum
  interface par_sum
    module procedure real0D_par_sum
    module procedure real1D_par_sum
    module procedure real2D_par_sum
    module procedure real3D_par_sum
    module procedure int0D_par_sum
    module procedure int1D_par_sum
  end interface par_sum

contains

  !
  ! par_max
  !

  !
  ! real0D_par_max
  !
  subroutine real0D_par_max(xxl, xxg)
    use mpi
    use mpi_interface, only : MY_SIZE
    implicit none

    real, intent(out) :: xxg
    real, intent(in)  :: xxl

    integer:: ierror

    call mpi_allreduce(xxl, xxg, 1, MY_SIZE, MPI_MAX, MPI_COMM_WORLD, ierror)

  end subroutine real0D_par_max

  !
  ! int0D_par_max
  !
  subroutine int0D_par_max(xxl, xxg)
    use mpi
    use mpi_interface, only : MY_INT
    implicit none

    integer, intent(out) :: xxg
    integer, intent(in)  :: xxl

    integer :: ierror

    call mpi_allreduce(xxl, xxg, 1, MY_INT, MPI_MAX, MPI_COMM_WORLD, ierror)

  end subroutine int0D_par_max

  !
  ! par_sum
  !

  !
  ! real0D_par_sum
  !
  subroutine real0D_par_sum(xxl, xxg)
    use mpi
    use mpi_interface, only : MY_SIZE
    implicit none

    real, intent(out) :: xxg
    real, intent(in)  :: xxl

    integer:: ierror

    call mpi_allreduce(xxl, xxg, 1, MY_SIZE, MPI_SUM, MPI_COMM_WORLD, ierror)

  end subroutine real0D_par_sum

  !
  ! real1D_par_sum
  !
  subroutine real1D_par_sum(xxl, xxg, n)
    use mpi
    use mpi_interface, only : MY_SIZE
    implicit none

    integer, intent(in) :: n
    real, intent(out)   :: xxg(:)
    real, intent(in)    :: xxl(:)

    integer             :: ierror

    call mpi_allreduce(xxl, xxg, n, MY_SIZE, MPI_SUM, MPI_COMM_WORLD, ierror)

  end subroutine real1D_par_sum

  !
  ! real2D_par_sum
  !
  subroutine real2D_par_sum(xxl, xxg, n)
    use mpi
    use mpi_interface, only : MY_SIZE
    implicit none

    integer, intent(in) :: n
    real, intent(out)   :: xxg(:,:)
    real, intent(in)    :: xxl(:,:)

    integer             :: ierror

    call mpi_allreduce(xxl, xxg, n, MY_SIZE, MPI_SUM, MPI_COMM_WORLD, ierror)

  end subroutine real2D_par_sum

  !
  ! real3D_par_sum
  !
  subroutine real3D_par_sum(xxl, xxg, n)
    use mpi
    use mpi_interface, only : MY_SIZE
    implicit none

    integer, intent(in) :: n
    real, intent(out)   :: xxg(:,:,:)
    real, intent(in)    :: xxl(:,:,:)

    integer             :: ierror

    call mpi_allreduce(xxl, xxg, n, MY_SIZE, MPI_SUM, MPI_COMM_WORLD, ierror)

  end subroutine real3D_par_sum

  !
  ! int_scalar_par_sum
  !
  subroutine int0D_par_sum(xxl, xxg)
    use mpi
    use mpi_interface, only : MY_INT
    implicit none

    integer, intent(out) :: xxg
    integer, intent(in)  :: xxl
    integer              :: ierror

    call mpi_allreduce(xxl, xxg, 1, MY_INT, MPI_SUM, MPI_COMM_WORLD, ierror)

  end subroutine int0D_par_sum

  !
  ! int1D_par_sum
  !
  subroutine int1D_par_sum(xxl, xxg, n)
    use mpi
    use mpi_interface, only : MY_INT
    implicit none

    integer, intent(in)  :: n
    integer, intent(out) :: xxg(:)
    integer, intent(in)  :: xxl(:)
    integer              :: ierror

    call mpi_allreduce(xxl, xxg, n, MY_INT, MPI_SUM, MPI_COMM_WORLD, ierror)

  end subroutine int1D_par_sum

end module modutil_mpi
