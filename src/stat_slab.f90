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

integer :: ncid, irec, nstats_1, nstats_2, nstats_t, n_local, n_global, n_s_local, n_s_global
real, dimension(:), pointer :: u1, v1, t1, q1, l1, th1, b1, qv1, N2, ut, vt, tt, qt, &
                               u1_s, v1_s, t1_s, q1_s, l1_s, th1_s, b1_s, qv1_s, N2_s, ut_s, vt_s, tt_s, qt_s, &
                               u2, v2, w2, t2, q2, l2, uw, vw, tw, qw, lw, tq, b2, bw, th2, qv2, thqv, &
                               wDp, bDp, tDp, qDp, u2t, v2t, w2t, t2t, q2t, twt, qwt, tqt, &
                               u2_s, v2_s, w2_s, t2_s, q2_s, l2_s, uw_s, vw_s, tw_s, qw_s, lw_s, tq_s, b2_s, &
                               bw_s, th2_s, qv2_s, thqv_s, wDp_s, bDp_s, tDp_s, qDp_s, &
                               u2t_s, v2t_s, w2t_s, t2t_s, q2t_s, twt_s, qwt_s, tqt_s, &
                               w3_s, t2w_s, q2w_s, b2w_s, tw2_s, qw2_s, bw2_s, &
                               w3, t2w, q2w, b2w, tw2, qw2, bw2
real, dimension(:,:), allocatable, target :: x1, x2, xt, x1_s, x2_s, xt_s

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
    call create_stat_var(ncid, 'n', zt_name, my_nf90_int)

    ! Create first-order moment variables
    nstats_1 = 5
    call create_stat_var(ncid, 'u',  zt_name, my_nf90_real)
    call create_stat_var(ncid, 'v',  zt_name, my_nf90_real)
    call create_stat_var(ncid, 't',  zt_name, my_nf90_real)
    call create_stat_var(ncid, 'b',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 'N2', zt_name, my_nf90_real)
    if (level >= 1) then
       nstats_1 = nstats_1 + 3
       call create_stat_var(ncid, 'q',  zt_name, my_nf90_real)
       call create_stat_var(ncid, 'th', zt_name, my_nf90_real)
       call create_stat_var(ncid, 'qv', zt_name, my_nf90_real)
    end if
    if (level >= 2) then
       nstats_1 = nstats_1 + 1
       call create_stat_var(ncid, 'l', zt_name, my_nf90_real)
    end if

    nstats_t = 3
    call create_stat_var(ncid, 'ut', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'vt', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'tt', zt_name, my_nf90_real)
    if (level >= 1) then
       nstats_t = nstats_t + 1
       call create_stat_var(ncid, 'qt', zt_name, my_nf90_real)
    end if

    ! Create second-order moment variables
    nstats_2 = 22
    call create_stat_var(ncid, 'u2',  zt_name, my_nf90_real)
    call create_stat_var(ncid, 'v2',  zt_name, my_nf90_real)
    call create_stat_var(ncid, 'w2',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 'b2',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 't2',  zt_name, my_nf90_real)
    call create_stat_var(ncid, 'uw',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 'vw',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 'tw',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 'bw',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 'wDp', zm_name, my_nf90_real)
    call create_stat_var(ncid, 'bDp', zm_name, my_nf90_real)
    call create_stat_var(ncid, 'tDp', zm_name, my_nf90_real)
    call create_stat_var(ncid, 'w3',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 'b2w', zm_name, my_nf90_real)
    call create_stat_var(ncid, 't2w', zm_name, my_nf90_real)
    call create_stat_var(ncid, 'bw2', zm_name, my_nf90_real)
    call create_stat_var(ncid, 'tw2', zm_name, my_nf90_real)
    call create_stat_var(ncid, 'u2t', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'v2t', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'w2t', zm_name, my_nf90_real)
    call create_stat_var(ncid, 't2t', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'twt', zm_name, my_nf90_real)
    if (level >= 1) then
       nstats_2 = nstats_2 + 12
       call create_stat_var(ncid, 'q2',   zt_name, my_nf90_real)
       call create_stat_var(ncid, 'th2',  zt_name, my_nf90_real)
       call create_stat_var(ncid, 'qv2',  zt_name, my_nf90_real)
       call create_stat_var(ncid, 'tq',   zt_name, my_nf90_real)
       call create_stat_var(ncid, 'thqv', zt_name, my_nf90_real)
       call create_stat_var(ncid, 'qw',   zm_name, my_nf90_real)
       call create_stat_var(ncid, 'qDp',  zm_name, my_nf90_real)
       call create_stat_var(ncid, 'q2w',  zm_name, my_nf90_real)
       call create_stat_var(ncid, 'qw2',  zm_name, my_nf90_real)
       call create_stat_var(ncid, 'q2t',  zt_name, my_nf90_real)
       call create_stat_var(ncid, 'tqt',  zt_name, my_nf90_real)
       call create_stat_var(ncid, 'qwt',  zm_name, my_nf90_real)
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
    b1 => x1(:,3)
    N2 => x1(:,4)
    t1 => x1(:,5)
    if (level >= 1) then
       q1  => x1(:,6)
       th1 => x1(:,7)
       qv1 => x1(:,8)
    end if
    if (level >= 2) l1 => x1(:,9)

    allocate(x1_s(nzp,nstats_1))
    u1_s => x1_s(:,1)
    v1_s => x1_s(:,2)
    b1_s => x1_s(:,3)
    N2_s => x1_s(:,4)
    t1_s => x1_s(:,5)
    if (level >= 1) then
       q1_s  => x1_s(:,6)
       th1_s => x1_s(:,7)
       qv1_s => x1_s(:,8)
    end if
    if (level >= 2) l1_s => x1_s(:,9)

    allocate(xt(nzp,nstats_t)) 
    xt(:,:) = 0.
    ut => xt(:,1)
    vt => xt(:,2)
    tt => xt(:,3)
    if (level >= 1) qt => xt(:,4)

    allocate(xt_s(nzp,nstats_t)) 
    xt_s(:,:) = 0.
    ut_s => xt_s(:,1)
    vt_s => xt_s(:,2)
    tt_s => xt_s(:,3)
    if (level >= 1) qt_s => xt_s(:,4)

    ! Allocate memory for second-order moment statistics
    allocate(x2(nzp,nstats_2))
    x2(:,:) = 0.
    u2  => x2(:,1)
    v2  => x2(:,2)
    b2  => x2(:,3)
    w2  => x2(:,4)
    t2  => x2(:,5)
    uw  => x2(:,6)
    vw  => x2(:,7)
    bw  => x2(:,8)
    tw  => x2(:,9)
    wDp => x2(:,10)
    bDp => x2(:,11)
    tDp => x2(:,12)
    w3  => x2(:,13)
    b2w => x2(:,14)
    t2w => x2(:,15)
    bw2 => x2(:,16)
    tw2 => x2(:,17)
    u2t => x2(:,18)
    v2t => x2(:,19)
    w2t => x2(:,20)
    t2t => x2(:,21)
    twt => x2(:,22)
    if (level >= 1) then
       q2   => x2(:,23)
       qv2  => x2(:,24)
       th2  => x2(:,25)
       qw   => x2(:,26)
       qDp  => x2(:,27)
       tq   => x2(:,28)
       thqv => x2(:,29)
       q2w  => x2(:,30)       
       qw2  => x2(:,31)
       q2t  => x2(:,32)
       qwt  => x2(:,33)
       tqt  => x2(:,34)
    end if
    if (level >= 2) then
       l2 => x2(:,35)
       lw => x2(:,36)
    end if

    ! Allocate memory for second-order moment statistics
    allocate(x2_s(nzp,nstats_2))
    x2_s(:,:) = 0.
    u2_s  => x2_s(:,1)
    v2_s  => x2_s(:,2)
    w2_s  => x2_s(:,3)
    b2_s  => x2_s(:,4)
    t2_s  => x2_s(:,5)
    uw_s  => x2_s(:,6)
    vw_s  => x2_s(:,7)
    bw_s  => x2_s(:,8)
    tw_s  => x2_s(:,9)
    wDp_s => x2_s(:,10)
    bDp_s => x2_s(:,11)
    tDp_s => x2_s(:,12)
    w3_s  => x2_s(:,13)
    b2w_s => x2_s(:,14)
    t2w_s => x2_s(:,15)
    bw2_s => x2_s(:,16)
    tw2_s => x2_s(:,17)
    u2t_s => x2_s(:,18)
    v2t_s => x2_s(:,19)
    w2t_s => x2_s(:,20)
    t2t_s => x2_s(:,21)
    twt_s => x2_s(:,22)
    if (level >= 1) then
       q2_s   => x2_s(:,23)
       qv2_s  => x2_s(:,24)
       th2_s  => x2_s(:,25)
       qw_s   => x2_s(:,26)
       qDp_s  => x2_s(:,27)
       tq_s   => x2_s(:,28)
       thqv_s => x2_s(:,29)
       q2w_s  => x2_s(:,30)       
       qw2_s  => x2_s(:,31)
       q2t_s  => x2_s(:,32)
       qwt_s  => x2_s(:,33)
       tqt_s  => x2_s(:,34)
    end if
    if (level >= 2) then
       l2_s => x2_s(:,35)
       lw_s => x2_s(:,36)
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
    use grid, only        : a_up, a_vp, a_wp, a_tp, a_rp, liquid, nzp, nyp, nxp, level, a_b, a_theta, vapor, &
                            zm, a_pexnr, zt, a_ut, a_vt, a_wt, a_tt, a_rt
    use modutil_mpi, only : par_sum
    implicit none

    real :: f_local, f_global
    real, allocatable, dimension(:,:) :: x1_scratch, xt_scratch
    integer :: i, j, k, kp1

    allocate(x1_scratch(nzp,nstats_1), xt_scratch(nzp,nstats_t))
    
    ! Sample first-order statistics
    x1_s(:,:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 1,nzp
          kp1 = max(k,nzp-1)
          
          u1_s(k) = u1_s(k) + a_up(k,i,j)
          v1_s(k) = v1_s(k) + a_vp(k,i,j)
          b1_s(k) = b1_s(k) + 0.5*(a_b(k,i,j) + a_b(kp1,i,j))
          N2_s(k) = N2_s(k) + a_b(k,i,j) ! temporary
          t1_s(k) = t1_s(k) + a_tp(k,i,j)
          if (level >= 1) then
             q1_s(k)  = q1_s(k)  + a_rp(k,i,j)
             th1_s(k) = th1_s(k) + a_theta(k,i,j)
             qv1_s(k) = qv1_s(k) + vapor(k,i,j)
          end if
          if (level >= 2) l1_s(k) = l1_s(k) + liquid(k,i,j)
       end do
    end do
    end do
    call par_sum(x1_s, x1_scratch, nzp*nstats_1)
    x1_s(:,:) = x1_scratch(:,:)/real(n_s_global)

    ! Compute N**2. from buoyancy profile temporarily stored in N2
    N2_s(2:nzp-1) = (N2_s(3:nzp) - N2_s(2:nzp-1))/(zm(3:nzp) - zm(2:nzp-1))

    ! Sample first-order statistics
    xt_s(:,:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 1,nzp
          ut_s(k) = ut_s(k) + a_ut(k,i,j)
          vt_s(k) = vt_s(k) + a_vt(k,i,j)
          tt_s(k) = tt_s(k) + a_tt(k,i,j)
          if (level >= 1) qt_s(k) = qt_s(k) + a_rt(k,i,j)
       end do
    end do
    end do
    call par_sum(xt_s, xt_scratch, nzp*nstats_t)
    xt_s(:,:) = xt_scratch(:,:)/real(n_s_global)

    ! Sample high-order statistics
    x2_s(:,:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 2,nzp-1
          kp1 = max(k,nzp-1)

          u2_s(k)  = u2_s(k)  + (a_up(k,i,j) - u1_s(k))**2.
          v2_s(k)  = v2_s(k)  + (a_vp(k,i,j) - v1_s(k))**2.
          w2_s(k)  = w2_s(k)  + a_wp(k,i,j)**2.
          b2_s(k)  = b2_s(k)  + (0.5*(a_b(k,i,j) + a_b(kp1,i,j)) - b1_s(k))**2.
          t2_s(k)  = t2_s(k)  + (a_tp(k,i,j) - t1_s(k))**2.
          uw_s(k)  = uw_s(k)  + a_wp(k,i,j)*(0.5*(a_up(k,i,j) + a_up(kp1,i,j)) - 0.5*(u1_s(k) + u1_s(kp1)))
          vw_s(k)  = vw_s(k)  + a_wp(k,i,j)*(0.5*(a_vp(k,i,j) + a_vp(kp1,i,j)) - 0.5*(v1_s(k) + v1_s(kp1)))
          bw_s(k)  = bw_s(k)  + a_wp(k,i,j)*(0.5*(a_b(k,i,j) + a_b(kp1,i,j)) - b1_s(k))          
          tw_s(k)  = tw_s(k)  + a_wp(k,i,j)*(0.5*(a_tp(k,i,j) + a_tp(kp1,i,j)) - 0.5*(t1_s(k) + t1_s(kp1)))
          w3_s(k)  = w3_s(k)  + a_wp(k,i,j)**3.
          b2w_s(k) = b2w_s(k) + a_wp(k,i,j)*(0.5*(a_b(k,i,j) + a_b(kp1,i,j)) - b1_s(k))**2.
          t2w_s(k) = t2w_s(k) + a_wp(k,i,j)*(0.5*(a_tp(k,i,j) + a_tp(kp1,i,j)) - 0.5*(t1_s(k) + t1_s(kp1)))**2.
          bw2_s(k) = bw2_s(k) + a_wp(k,i,j)**2.*(0.5*(a_b(k,i,j) + a_b(kp1,i,j)) - b1_s(k))
          tw2_s(k) = tw2_s(k) + a_wp(k,i,j)**2.*(0.5*(a_tp(k,i,j) + a_tp(kp1,i,j)) - 0.5*(t1_s(k) + t1_s(kp1)))
          wDp_s(k) = wDp_s(k) + a_wp(k,i,j)*(a_pexnr(kp1,i,j) - a_pexnr(k,i,j))/(zt(kp1) - zt(k))
          bDp_s(k) = bDp_s(k) + 0.5*(a_b(k,i,j) + a_b(kp1,i,j))*(a_pexnr(kp1,i,j) - a_pexnr(k,i,j))/(zt(kp1) - zt(k))
          tDp_s(k) = tDp_s(k) + (0.5*(a_tp(k,i,j) + a_tp(kp1,i,j)) - 0.5*(t1_s(k) + t1_s(kp1)))*(a_pexnr(kp1,i,j) - a_pexnr(k,i,j))/(zt(kp1) - zt(k))
          u2t_s(k)  = u2t_s(k)  + 2.*(a_up(k,i,j) - u1_s(k))*(a_ut(k,i,j) - ut_s(k))
          v2t_s(k)  = v2t_s(k)  + 2.*(a_vp(k,i,j) - v1_s(k))*(a_vt(k,i,j) - vt_s(k))
          w2t_s(k)  = w2t_s(k)  + 2.*a_wp(k,i,j)*a_wt(k,i,j)
          t2t_s(k)  = t2t_s(k)  + 2.*(a_tp(k,i,j) - t1_s(k))*(a_tt(k,i,j) - tt_s(k))
          twt_s(k)  = twt_s(k)  + a_wp(k,i,j)*(0.5*(a_tp(k,i,j) + a_tp(kp1,i,j)) - 0.5*(t1_s(k) + t1_s(kp1))) &
                                * a_wt(k,i,j)*(0.5*(a_tt(k,i,j) + a_tt(kp1,i,j)) - 0.5*(tt_s(k) + tt_s(kp1)))
          if (level >= 1) then
             q2_s(k)   = q2_s(k)   + (a_rp(k,i,j) - q1_s(k))**2.
             th2_s(k)  = th2_s(k)  + (a_theta(k,i,j) - th1_s(k))**2.
             qv2_s(k)  = qv2_s(k)  + (vapor(k,i,j) - qv1_s(k))**2.
             tq_s(k)   = tq_s(k)   + (a_tp(k,i,j) - t1_s(k))*(a_rp(k,i,j) - q1_s(k))
             qw_s(k)   = qw_s(k)   + a_wp(k,i,j)*(0.5*(a_rp(k,i,j) + a_rp(kp1,i,j)) - 0.5*(q1_s(k) + q1_s(kp1)))
             thqv_s(k) = thqv_s(k) + (a_theta(k,i,j) - th1_s(k))*(vapor(k,i,j) - qv1_s(k))
             qDp_s(k)  = qDp_s(k)  + (0.5*(a_rp(k,i,j) + a_rp(kp1,i,j)) - 0.5*(q1_s(k) + q1_s(kp1)))*(a_pexnr(kp1,i,j) - a_pexnr(k,i,j))/(zt(kp1) - zt(k))
             q2w_s(k)  = q2w_s(k)  + a_wp(k,i,j)*(0.5*(a_rp(k,i,j) + a_rp(kp1,i,j)) - 0.5*(q1_s(k) + q1_s(kp1)))**2.
             qw2_s(k)  = qw2_s(k)  + a_wp(k,i,j)**2.*(0.5*(a_rp(k,i,j) + a_rp(kp1,i,j)) - 0.5*(q1_s(k) + q1_s(kp1)))
             q2t_s(k)  = q2t_s(k)  + 2.*(a_rp(k,i,j) - q1_s(k))*(a_rt(k,i,j) - qt_s(k))
             tqt_s(k)  = tqt_s(k)  + (a_tp(k,i,j) - t1_s(k))*(a_rp(k,i,j) - q1_s(k)) &
                                   + (a_rt(k,i,j) - qt_s(k))*(a_rt(k,i,j) - qt_s(k))
             qwt_s(k)  = qwt_s(k)  + a_wp(k,i,j)*(0.5*(a_rp(k,i,j) + a_rp(kp1,i,j)) - 0.5*(q1_s(k) + q1_s(kp1))) &
                                   * a_wt(k,i,j)*(0.5*(a_rt(k,i,j) + a_rt(kp1,i,j)) - 0.5*(qt_s(k) + qt_s(kp1)))
          end if
          if (level >= 2) then
             l2_s(k) = l2_s(k) + (liquid(k,i,j) - l1_s(k))**2.
             lw_s(k) = lw_s(k) + a_wp(k,i,j)*(0.5*(liquid(k,i,j) + liquid(kp1,i,j)) - 0.5*(l1_s(k) + l1_s(kp1)))
          end if
       end do
    end do
    end do
    x2_s(:,:) = x2_s(:,:)/real(n_s_local)

    ! Update third-order statistics
    f_local = real(n_s_local)/real(n_s_local + n_local)
    w3(:)  = f_local*w3_s(:)  + (1. - f_local)*w3(:)
    b2w(:) = f_local*b2w_s(:) + (1. - f_local)*b2w(:) + 2.*f_local*(1. - f_local)*(b1_s(:) - b1(:))*(bw_s(:) - bw(:))
    t2w(:) = f_local*t2w_s(:) + (1. - f_local)*t2w(:) + 2.*f_local*(1. - f_local)*(t1_s(:) - t1(:))*(tw_s(:) - tw(:))
    bw2(:) = f_local*bw2_s(:) + (1. - f_local)*bw2(:) + f_local*(1. - f_local)*(b1_s(:) - b1(:))*(w2_s(:) - w2(:))
    tw2(:) = f_local*tw2_s(:) + (1. - f_local)*tw2(:) + f_local*(1. - f_local)*(t1_s(:) - t1(:))*(w2_s(:) - w2(:))
    if (level >= 1) then
       q2w(:) = f_local*q2w_s(:) + (1. - f_local)*q2w(:) + 2.*f_local*(1. - f_local)*(q1_s(:) - q1(:))*(qw_s(:) - qw(:))
       qw2(:) = f_local*qw2_s(:) + (1. - f_local)*qw2(:) + f_local*(1. - f_local)*(q1_s(:) - q1(:))*(w2_s(:) - w2(:))
    end if

    ! Update second-order statistics
    u2(:)  = f_local*u2_s(:)  + (1. - f_local)*u2(:)  + f_local*(1. - f_local)*(u1_s(:) - u1(:))**2.
    v2(:)  = f_local*v2_s(:)  + (1. - f_local)*v2(:)  + f_local*(1. - f_local)*(v1_s(:) - v1(:))**2.
    w2(:)  = f_local*w2_s(:)  + (1. - f_local)*w2(:)
    b2(:)  = f_local*b2_s(:)  + (1. - f_local)*b2(:)  + f_local*(1. - f_local)*(b1_s(:) - b1(:))**2.
    t2(:)  = f_local*t2_s(:)  + (1. - f_local)*t2(:)  + f_local*(1. - f_local)*(t1_s(:) - t1(:))**2.
    uw(:)  = f_local*uw_s(:)  + (1. - f_local)*uw(:)
    vw(:)  = f_local*vw_s(:)  + (1. - f_local)*vw(:)
    bw(:)  = f_local*bw_s(:)  + (1. - f_local)*bw(:)
    tw(:)  = f_local*tw_s(:)  + (1. - f_local)*tw(:)
    wDp(:) = f_local*wDp_s(:) + (1. - f_local)*wDp(:)
    bDp(:) = f_local*bDp_s(:) + (1. - f_local)*bDp(:)
    tDp(:) = f_local*tDp_s(:) + (1. - f_local)*tDp(:)
    u2t(:) = f_local*u2t_s(:) + (1. - f_local)*u2t(:) + f_local*(1. - f_local)*(ut_s(:) - ut(:))**2.
    v2t(:) = f_local*v2t_s(:) + (1. - f_local)*v2t(:) + f_local*(1. - f_local)*(vt_s(:) - vt(:))**2.
    u2t(:) = f_local*w2t_s(:) + (1. - f_local)*w2t(:)
    t2t(:) = f_local*t2t_s(:) + (1. - f_local)*t2t(:) + f_local*(1. - f_local)*(tt_s(:) - tt(:))**2.
    twt(:) = f_local*twt_s(:) + (1. - f_local)*twt(:)
    if (level >= 1) then
       q2(:)   = f_local*q2_s(:)   + (1. - f_local)*q2(:)   + f_local*(1. - f_local)*(q1_s(:) - q1(:))**2.
       th2(:)  = f_local*th2_s(:)  + (1. - f_local)*th2(:)  + f_local*(1. - f_local)*(th1_s(:) - th1(:))**2.
       qv2(:)  = f_local*qv2_s(:)  + (1. - f_local)*qv2(:)  + f_local*(1. - f_local)*(qv1_s(:) - qv1(:))**2.
       qw(:)   = f_local*qw_s(:)   + (1. - f_local)*qw(:)
       tq(:)   = f_local*tq_s(:)   + (1. - f_local)*tq(:)   + f_local*(1. - f_local)*(t1_s(:) - t1(:))*(q1_s(:) - q1(:))
       thqv(:) = f_local*thqv_s(:) + (1. - f_local)*thqv(:) + f_local*(1. - f_local)*(th1_s(:) - th1(:))*(qv1_s(:) - qv1(:))
       qDp(:)  = f_local*qDp_s(:)  + (1. - f_local)*qDp(:)
       q2t(:)  = f_local*q2t_s(:)  + (1. - f_local)*q2t(:)  + f_local*(1. - f_local)*(qt_s(:) - qt(:))**2.
       tqt(:)  = f_local*tqt_s(:)  + (1. - f_local)*tqt(:)  + f_local*(1. - f_local)*(tt_s(:) - tt(:))*(qt_s(:) - qt(:))
       qwt(:)  = f_local*qwt_s(:)  + (1. - f_local)*qwt(:)
    end if
    if (level >= 2) then
       l2(:) = f_local*l2_s(:) + (1. - f_local)*l2(:) + f_local*(1. - f_local)*(l1_s(:) - l1(:))**2.
       lw(:) = f_local*lw_s(:) + (1. - f_local)*lw(:)
    end if

    ! Update first-order statistics
    f_global = real(n_s_global)/real(n_s_global + n_global)
    x1(:,:)  = f_global*x1_s(:,:) + (1. - f_global)*x1(:,:)
    xt(:,:)  = f_global*xt_s(:,:) + (1. - f_global)*xt(:,:)

    ! Update sample count
    n_local  = n_local + n_s_local
    n_global = n_global + n_s_global

  end subroutine stat_slab

  !
  ! write_stat_slab
  !
  subroutine write_stat_slab(time, fsttm, nsmp)
    use modstat_io, only : write_stat_var, advance_stat_time
    use grid, only       : level, umean, vmean, th00, nzp
    implicit none

    real, intent(in)     :: time, fsttm
    integer, intent(in)  :: nsmp

    integer, allocatable :: n_array(:)

    !
    allocate(n_array(nzp))
    n_array(:) = n_local

    ! Add mean state
    u1(:) = u1(:) + umean
    v1(:) = v1(:) + vmean
    t1(:) = t1(:) + th00

    ! Write time
    call advance_stat_time(time, fsttm, nsmp, ncid, irec)

    ! Write statistics
    call write_stat_var(ncid, 'n',   n_array, irec)
    call write_stat_var(ncid, 'u',   u1,      irec)
    call write_stat_var(ncid, 'v',   v1,      irec)
    call write_stat_var(ncid, 'b',   b1,      irec)
    call write_stat_var(ncid, 'N2',  N2,      irec)
    call write_stat_var(ncid, 't',   t1,      irec)

    call write_stat_var(ncid, 'ut',  ut,      irec)
    call write_stat_var(ncid, 'vt',  ut,      irec)
    call write_stat_var(ncid, 'tt',  tt,      irec)

    call write_stat_var(ncid, 'u2',  u2,      irec)
    call write_stat_var(ncid, 'v2',  v2,      irec)
    call write_stat_var(ncid, 'w2',  w2,      irec)
    call write_stat_var(ncid, 'b2',  b2,      irec)
    call write_stat_var(ncid, 't2',  t2,      irec)
    call write_stat_var(ncid, 'uw',  uw,      irec)
    call write_stat_var(ncid, 'vw',  vw,      irec)
    call write_stat_var(ncid, 'bw',  bw,      irec)
    call write_stat_var(ncid, 'tw',  tw,      irec)
    call write_stat_var(ncid, 'wDp', wDp,     irec)
    call write_stat_var(ncid, 'bDp', bDp,     irec)
    call write_stat_var(ncid, 'tDp', tDp,     irec)
    call write_stat_var(ncid, 'w3',  w3,      irec)
    call write_stat_var(ncid, 'b2w', b2w,     irec)
    call write_stat_var(ncid, 't2w', t2w,     irec)
    call write_stat_var(ncid, 'bw2', bw2,     irec)
    call write_stat_var(ncid, 'tw2', tw2,     irec)

    call write_stat_var(ncid, 'u2t', u2t,     irec)
    call write_stat_var(ncid, 'v2t', v2t,     irec)
    call write_stat_var(ncid, 'w2t', w2t,     irec)
    call write_stat_var(ncid, 't2t', t2t,     irec)
    call write_stat_var(ncid, 'twt', twt,     irec)

    if (level >= 1) then
       call write_stat_var(ncid, 'q',    q1,   irec)
       call write_stat_var(ncid, 'th',   th1,  irec)
       call write_stat_var(ncid, 'qv',   qv1,  irec)
       call write_stat_var(ncid, 'q2',   q2,   irec)
       call write_stat_var(ncid, 'th2',  th2,  irec)
       call write_stat_var(ncid, 'qv2',  qv2,  irec)
       call write_stat_var(ncid, 'tq',   tq,   irec)
       call write_stat_var(ncid, 'thqv', thqv, irec)
       call write_stat_var(ncid, 'qw',   qw,   irec)
       call write_stat_var(ncid, 'qDp',  qDp,  irec)
       call write_stat_var(ncid, 'q2w',  q2w,  irec)
       call write_stat_var(ncid, 'qw2',  qw2,  irec)

       call write_stat_var(ncid, 'q2t',  q2t,  irec)
       call write_stat_var(ncid, 'tqt',  tqt,  irec)
       call write_stat_var(ncid, 'qwt',  qwt,  irec)
    end if
    if (level >= 2) then
       call write_stat_var(ncid, 'l',  l1, irec)
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

