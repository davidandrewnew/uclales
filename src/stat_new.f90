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
! Copyright 2015-2018, David A. New
!----------------------------------------------------------------------------

module modstat
  implicit none

  ! Namelist parameters
  real    :: ssam_intvl = 30.   ! statistical sampling interval
  real    :: savg_intvl = 1800. ! statistical averaging interval

contains

  !
  ! init_stat
  !
  subroutine init_stat(time, runtype, cntlat)
    use grid, only         : nzp, nxp, nyp, nsmp
    use modstat_slab, only : init_stat_slab
    implicit none

    real, intent(in)              :: time, cntlat
    character (len=7), intent(in) :: runtype

    ! Initialize
    nsmp = 0

    ! Slab statitics
    call init_stat_slab(runtype, cntlat)

  end subroutine init_stat

  !
  ! sample_stat
  !
  subroutine sample_stat(time)
    use grid, only : nsmp, fsttm
    use modstat_slab, only : sample_stat_slab
    implicit none
    
    real, intent(in) :: time

    !
    if (nsmp == 0) fsttm = time

    ! Grid statistics
    call sample_stat_slab

  end subroutine sample_stat

  !
  ! stat_tendency
  !
  subroutine stat_tendency(itype, dt)
    use grid, only : a_up, a_vp, a_wp, a_ut, a_vt, a_wt, a_scr4, a_scr5, a_scr6, a_scr7, rkalpha, a_pexnr, nzp, nxp, nyp
    use util, only : velset
    use prss, only : poisson
    use modstat_slab, only : stat_slab_pressure, stat_slab_tendency
    use mpi_interface, only : myid

    integer, intent(in) :: itype
    real, intent(in)    :: dt

    a_scr4(:,:,:) = a_up(:,:,:) + rkalpha(1)*dt*a_ut(:,:,:)
    a_scr5(:,:,:) = a_vp(:,:,:) + rkalpha(1)*dt*a_vt(:,:,:)
    a_scr6(:,:,:) = a_wp(:,:,:) + rkalpha(1)*dt*a_wt(:,:,:)
    call velset(nzp, nxp, nyp, a_scr4, a_scr5, a_scr6)
    a_scr7(:,:,:) = a_pexnr(:,:,:)
    call poisson(dt, a_scr4, a_scr5, a_scr6, a_scr7)

    call stat_slab_pressure(itype, a_up, a_vp, a_wp, a_scr7)
    call stat_slab_tendency(itype)
  end subroutine stat_tendency

  !
  ! update_stat
  !
  subroutine update_stat
    use grid, only : nsmp
    use modstat_slab, only : update_stat_slab
    implicit none

    ! Grid statistics
    call update_stat_slab

    ! Increment sample count
    nsmp = nsmp + 1

  end subroutine update_stat

  !
  ! write_stat
  !
  subroutine write_stat(time)
    use grid, only          : nsmp, fsttm
    use mpi_interface, only : appl_abort
    use modstat_slab, only  : write_stat_slab
    implicit none

    real, intent(in) :: time 

    ! Grid statistics
    call write_stat_slab(time, fsttm, nsmp)
       
    ! Reset sample count
    fsttm = time
    nsmp  = 0

    ! Re-use last sample for last average as first sample for next average
    call update_stat
 
  end subroutine write_stat

  !
  ! exit_stat
  !
  subroutine exit_stat
    use modstat_slab, only : exit_stat_slab
    implicit none

    ! Grid statistics
    call exit_stat_slab

  end subroutine exit_stat


end module modstat




