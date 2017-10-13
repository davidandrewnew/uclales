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

module modstat_slab

implicit none

integer :: ncid, irec, nstats_1, nstats_2, n_local, n_global, n_s_local, n_s_global
real, dimension(:), pointer :: u1, v1, t1, q1, l1, th1, b1, qv1, &
                               u1_s, v1_s, t1_s, q1_s, l1_s, th1_s, b1_s, qv1_s, &
                               u2, v2, w2, t2, q2, l2, uw, vw, tw, qw, lw, tq, b2, bw, th2, qv2, thqv, &
                               u2_s, v2_s, w2_s, t2_s, q2_s, l2_s, uw_s, vw_s, tw_s, qw_s, lw_s, tq_s, b2_s, &
                               bw_s, th2_s, qv2_s, thqv_s  
real, dimension(:,:), allocatable, target :: x1, x2, x1_s, x2_s, x1_scratch

contains

  !
  ! init_stat_slab
  !
  subroutine init_stat_slab
    use grid, only             : filprf, level, nzp, nxp, nyp
    use modstat_io, only       : create_stat_file, create_stat_var, zt_name, zm_name
    use netcdf_interface, only : my_nf90_real, my_nf90_int, create_var
    use modutil_mpi, only      : par_sum
    implicit none

    ! Create statistics file
    call create_stat_file(trim(filprf)//'.stat.slab.nc', ncid, irec)

    ! Create first-order moment variables
    nstats_1 = 3
    call create_stat_var(ncid, 'u', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'v', zt_name, my_nf90_real)
    call create_stat_var(ncid, 't', zt_name, my_nf90_real)
    if (level >= 1) then
       nstats_1 = nstats_1 + 4
       call create_stat_var(ncid, 'q',  zt_name, my_nf90_real)
       call create_stat_var(ncid, 'b',  zt_name, my_nf90_real)
       call create_stat_var(ncid, 'th', zt_name, my_nf90_real)
       call create_stat_var(ncid, 'qv', zt_name, my_nf90_real)
    end if
    if (level >= 2) then
       nstats_1 = nstats_1 + 1
       call create_stat_var(ncid, 'l', zt_name, my_nf90_real)
    end if
    
    ! Create second-order moment variables
    nstats_2 = 7
    call create_stat_var(ncid, 'u2', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'v2', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'w2', zt_name, my_nf90_real)
    call create_stat_var(ncid, 't2', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'uw', zm_name, my_nf90_real)
    call create_stat_var(ncid, 'vw', zm_name, my_nf90_real)
    call create_stat_var(ncid, 'tw', zm_name, my_nf90_real)
    if (level >= 1) then
       nstats_2 = nstats_2 + 8
       call create_stat_var(ncid, 'q2',   zt_name, my_nf90_real)
       call create_stat_var(ncid, 'qw',   zm_name, my_nf90_real)
       call create_stat_var(ncid, 'tq',   zt_name, my_nf90_real)
       call create_stat_var(ncid, 'b2',   zt_name, my_nf90_real)
       call create_stat_var(ncid, 'bw',   zm_name, my_nf90_real)
       call create_stat_var(ncid, 'th2',  zt_name, my_nf90_real)
       call create_stat_var(ncid, 'qv2',  zt_name, my_nf90_real)
       call create_stat_var(ncid, 'thqv', zt_name, my_nf90_real)
    end if
    if (level >= 2) then
       nstats_2 = nstats_2 + 2
       call create_stat_var(ncid, 'l2', zt_name, my_nf90_real)
       call create_stat_var(ncid, 'lw', zm_name, my_nf90_real)
    end if

    ! Allocate memory for first-order moment statistics
    allocate(x1(nzp,nstats_1)) 
    x1(:,:) = 0.
    u1 => x1(:,1)
    v1 => x1(:,2)
    t1 => x1(:,3)
    if (level >= 1) then
       q1  => x1(:,4)
       b1  => x1(:,5)
       th1 => x1(:,6)
       qv1 => x1(:,7)
    end if
    if (level >= 2) l1 => x1(:,8)

    allocate(x1_s(nzp,nstats_1), x1_scratch(nzp,nstats_1))
    u1_s => x1_s(:,1)
    v1_s => x1_s(:,2)
    t1_s => x1_s(:,3)
    if (level >= 1) then
       q1_s  => x1_s(:,4)
       b1_s  => x1_s(:,5)
       th1_s => x1_s(:,6)
       qv1_s => x1_s(:,7)
    end if
    if (level >= 2) l1_s => x1_s(:,8)

    ! Allocate memory for second-order moment statistics
    allocate(x2(nzp,nstats_2))
    x2(:,:) = 0.
    u2 => x2(:,1)
    v2 => x2(:,2)
    w2 => x2(:,3)
    t2 => x2(:,4)
    uw => x2(:,5)
    vw => x2(:,6)
    tw => x2(:,7)
    if (level >= 1) then
       q2   => x2(:,8)
       qw   => x2(:,9)
       tq   => x2(:,10)
       b2   => x2(:,11)
       bw   => x2(:,12)
       th2  => x2(:,13)
       qv2  => x2(:,14)
       thqv => x2(:,15)
    end if
    if (level >= 2) then
       l2 => x2(:,16)
       lw => x2(:,17)
    end if

    allocate(x2_s(nzp,nstats_2))
    u2_s => x2_s(:,1)
    v2_s => x2_s(:,2)
    w2_s => x2_s(:,3)
    t2_s => x2_s(:,4)
    uw_s => x2_s(:,5)
    vw_s => x2_s(:,6)
    tw_s => x2_s(:,7)
    if (level >= 1) then
       q2_s   => x2_s(:,8)
       qw_s   => x2_s(:,9)
       tq_s   => x2_s(:,10)
       b2_s   => x2_s(:,11)
       bw_s   => x2_s(:,12)
       th2_s  => x2_s(:,13)
       qv2_s  => x2_s(:,14)
       thqv_s => x2_s(:,15)
    end if
    if (level >= 2) then
       l2_s => x2_s(:,16)
       lw_s => x2_s(:,17)
    end if

    ! Sample counts
    n_s_local = (nxp - 4)*(nyp - 4)
    call par_sum(n_s_local, n_s_global)

    n_local  = 0
    n_global = 0
    
  end subroutine init_stat_slab

  !
  ! stat_slab
  !
  subroutine stat_slab
    use grid, only        : a_up, a_vp, a_wp, a_tp, a_rp, liquid, nzp, nyp, nxp, level, a_b, a_theta, vapor
    use modutil_mpi, only : par_sum
    implicit none

    real :: f_local, f_global
    integer :: i, j, k

    ! Sample first-order statistics
    x1_s(:,:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 2,nzp-1
          u1_s(k) = u1_s(k) + a_up(k,i,j)
          v1_s(k) = v1_s(k) + a_vp(k,i,j)
          t1_s(k) = t1_s(k) + a_tp(k,i,j)
          if (level >= 1) then
             q1_s(k)  = q1_s(k)  + a_rp(k,i,j)
             b1_s(k)  = b1_s(k)  + a_b(k,i,j)
             th1_s(k) = th1_s(k) + a_theta(k,i,j)
             qv1_s(k) = qv1_s(k) + vapor(k,i,j)
          end if
          if (level >= 2) l1_s(k) = l1_s(k) + liquid(k,i,j)
       end do
    end do
    end do
    call par_sum(x1_s, x1_scratch, nzp*nstats_1)
    x1_s(:,:) = x1_scratch(:,:)/real(n_s_global)    

    ! Sample second-order statistics
    x2_s(:,:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 2,nzp-1
          u2_s(k) = u2_s(k) + (a_up(k,i,j) - u1_s(k))**2.
          v2_s(k) = v2_s(k) + (a_vp(k,i,j) - v1_s(k))**2.
          t2_s(k) = t2_s(k) + (a_tp(k,i,j) - t1_s(k))**2.
          uw_s(k) = uw_s(k) + a_wp(k,i,j)*(0.5*(a_up(k,i,j) + a_up(k+1,i,j)) - 0.5*(u1_s(k) + u1_s(k+1)))
          vw_s(k) = vw_s(k) + a_wp(k,i,j)*(0.5*(a_vp(k,i,j) + a_vp(k+1,i,j)) - 0.5*(v1_s(k) + v1_s(k+1)))
          tw_s(k) = tw_s(k) + a_wp(k,i,j)*(0.5*(a_tp(k,i,j) + a_tp(k+1,i,j)) - 0.5*(t1_s(k) + t1_s(k+1)))
          if (level >= 1) then
             q2_s(k)    = q2_s(k)   + (a_rp(k,i,j) - q1_s(k))**2.
             qw_s(k)    = qw_s(k)   + a_wp(k,i,j)*(0.5*(a_rp(k,i,j) + a_rp(k+1,i,j)) - 0.5*(q1_s(k) + q1_s(k+1)))
             tq_s(k)    = tq_s(k)   + (a_tp(k,i,j) - t1_s(k))*(a_rp(k,i,j) - q1_s(k))
             b2_s(k)    = b2_s(k)   + (a_b(k,i,j) - b1_s(k))**2.
             bw_s(k)    = bw_s(k)   + a_wp(k,i,j)*(0.5*(a_b(k,i,j) + a_b(k+1,i,j)) - 0.5*(b1_s(k) + b1_s(k+1)))
             th2_s(k)   = th2_s(k)  + (a_theta(k,i,j) - th1_s(k))**2.
             qv2_s(k)   = qv2_s(k)  + (vapor(k,i,j) - qv1_s(k))**2.
             thqv_s(k)  = thqv_s(k) + (a_theta(k,i,j) - th1_s(k))*(vapor(k,i,j) - qv1_s(k))
          end if
          if (level >= 2) then
             l2_s(k) = l2_s(k) + (liquid(k,i,j) - l1_s(k))**2.
             lw_s(k) = lw_s(k) + a_wp(k,i,j)*(0.5*(liquid(k,i,j) + liquid(k+1,i,j)) - 0.5*(l1_s(k) + l1_s(k+1)))
          end if
       end do
    end do
    end do
    x2_s(:,:) = x2_s(:,:)/real(n_s_local)

    ! Update second-order statistics
    f_local = real(n_s_local)/real(n_s_local + n_local)
    u2(:) = f_local*u2_s(:) + (1. - f_local)*u2(:) + f_local*(1. - f_local)*(u1_s(:) - u1(:))**2.
    v2(:) = f_local*v2_s(:) + (1. - f_local)*v2(:) + f_local*(1. - f_local)*(v1_s(:) - v1(:))**2.
    w2(:) = f_local*w2_s(:) + (1. - f_local)*w2(:)
    t2(:) = f_local*t2_s(:) + (1. - f_local)*t2(:) + f_local*(1. - f_local)*(t1_s(:) - t1(:))**2.
    uw(:) = f_local*uw_s(:) + (1. - f_local)*uw(:)
    vw(:) = f_local*vw_s(:) + (1. - f_local)*vw(:)
    tw(:) = f_local*tw_s(:) + (1. - f_local)*tw(:)
    if (level >= 1) then
       q2(:)   = f_local*q2_s(:)   + (1. - f_local)*q2(:)   + f_local*(1. - f_local)*(q1_s(:) - q1(:))**2.
       qw(:)   = f_local*qw_s(:)   + (1. - f_local)*qw(:)
       tq(:)   = f_local*tq_s(:)   + (1. - f_local)*tq(:)   + f_local*(1. - f_local)*(t1_s(:) - t1(:))*(q1_s(:) - q1(:))
       b2(:)   = f_local*b2_s(:)   + (1. - f_local)*b2(:)   + f_local*(1. - f_local)*(b1_s(:) - b1(:))**2.
       bw(:)   = f_local*bw_s(:)   + (1. - f_local)*bw(:)
       th2(:)  = f_local*th2_s(:)  + (1. - f_local)*th2(:)  + f_local*(1. - f_local)*(th1_s(:) - th1(:))**2.
       qv2(:)  = f_local*qv2_s(:)  + (1. - f_local)*qv2(:)  + f_local*(1. - f_local)*(qv1_s(:) - qv1(:))**2.
       thqv(:) = f_local*thqv_s(:) + (1. - f_local)*thqv(:) + f_local*(1. - f_local)*(th1_s(:) - th1(:))*(qv1_s(:) - qv1(:))
    end if
    if (level >= 2) then
       l2(:) = f_local*l2_s(:) + (1. - f_local)*l2(:) + f_local*(1. - f_local)*(l1_s(:) - l1(:))**2.
       lw(:) = f_local*lw_s(:) + (1. - f_local)*lw(:)
    end if
    n_local = n_local + n_s_local

    ! Update first-order statistics
    f_global = real(n_s_global)/real(n_s_global + n_global)
    x1(:,:)  = f_global*x1_s(:,:) + (1. - f_global)*x1(:,:)
    n_global = n_global + n_s_global

  end subroutine stat_slab

  !
  ! write_stat_slab
  !
  subroutine write_stat_slab(time, fsttm, nsmp)
    use modstat_io, only : write_stat_var, advance_stat_time
    use grid, only       : level, umean, vmean, th00
    implicit none

    real, intent(in)     :: time, fsttm
    integer, intent(in)  :: nsmp

    ! Add mean state
    u1(:) = u1(:) + umean
    v1(:) = v1(:) + vmean
    t1(:) = t1(:) + th00

    ! Write time
    call advance_stat_time(time, fsttm, nsmp, ncid, irec)

    ! Write statistics
    call write_stat_var(ncid, 'u', u1, irec)
    call write_stat_var(ncid, 'v', v1, irec)
    call write_stat_var(ncid, 't', t1, irec)
    if (level >= 1) then
       call write_stat_var(ncid, 'q',  q1,  irec)
       call write_stat_var(ncid, 'b',  b1,  irec)
       call write_stat_var(ncid, 'th', th1, irec)
       call write_stat_var(ncid, 'qv', qv1, irec)
    end if
    if (level >= 2) call write_stat_var(ncid, 'l', l1, irec)
    call write_stat_var(ncid, 'u2', u2, irec)
    call write_stat_var(ncid, 'v2', v2, irec)
    call write_stat_var(ncid, 't2', t2, irec)
    call write_stat_var(ncid, 'uw', uw, irec)
    call write_stat_var(ncid, 'vw', vw, irec)
    call write_stat_var(ncid, 'tw', tw, irec)
    if (level >= 1) then
       call write_stat_var(ncid, 'q2',   q2,   irec)
       call write_stat_var(ncid, 'qw',   qw,   irec)
       call write_stat_var(ncid, 'tq',   tq,   irec)
       call write_stat_var(ncid, 'b2',   b2,   irec)
       call write_stat_var(ncid, 'bw',   bw,   irec)
       call write_stat_var(ncid, 'th2',  th2,  irec)
       call write_stat_var(ncid, 'qv2',  qv2,  irec)
       call write_stat_var(ncid, 'thqv', thqv, irec)
    end if
    if (level >= 2) then
       call write_stat_var(ncid, 'l2', l2, irec)
       call write_stat_var(ncid, 'lw', lw, irec)
    end if

    ! Reset statistics computation
    n_local  = 0.
    n_global = 0.
    x1(:,:)  = 0.
    x2(:,:)  = 0.
  
  end subroutine write_stat_slab

  !
  ! exit_stat_slab
  !
  subroutine exit_stat_slab
    use netcdf_interface, only : close_file
    implicit none

    call close_file(ncid)

  end subroutine exit_stat_slab

end module modstat_slab

