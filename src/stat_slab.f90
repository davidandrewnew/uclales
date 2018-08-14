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

integer :: ncid, irec, nstats_x1, nstats_x2, nstats_x1_misc, nstats_x2_misc, nstats_x_surf, &
           n_local, n_global, n_s_local, n_s_global

real, dimension(:), pointer :: u1, v1, w1, t1, q1, u1_s, v1_s, w1_s, t1_s, q1_s, &
                               u2, v2, w2, t2, q2, uw, vw, tw, qw, tq, &
                               u2_s, v2_s, w2_s, t2_s, q2_s, uw_s, vw_s, tw_s, qw_s, tq_s
real, dimension(:,:), pointer :: ut, vt, wt, tt, qt, ut_s, vt_s, wt_s, tt_s, qt_s, &
                                 ut_press, vt_press, wt_press, tt_press, qt_press, ut_press_s, vt_press_s, wt_press_s, tt_press_s, qt_press_s, &
                                 u2t, v2t, w2t, t2t, q2t, uwt, vwt, twt, qwt, tqt, &
                                 u2t_press, v2t_press, w2t_press, t2t_press, q2t_press, uwt_press, vwt_press, twt_press, qwt_press, tqt_press, &
                                 u2t_s, v2t_s, w2t_s, t2t_s, q2t_s, uwt_s, vwt_s, twt_s, qwt_s, tqt_s, &
                                 u2t_press_s, v2t_press_s, w2t_press_s, t2t_press_s, q2t_press_s, uwt_press_s, vwt_press_s, twt_press_s, qwt_press_s, tqt_press_s
real, dimension(:), pointer :: b1, N2, pi11, wfls1, dthldtls1, dqtdtls1, th1, l1, &
                               e, et_diss, et_shear, e_s, et_diss_s, et_shear_s, &
                               km, kh, km_s, kh_s, &
                               b1_s, N2_s, pi11_s, wfls1_s, dthldtls1_s, dqtdtls1_s, th1_s, l1_s, &
                               b2, bw, thb, lb, w3, tw2, qw2, bw2, t2w, q2w, b2w, &
                               b2_s, bw_s, thb_s, lb_s, w3_s, tw2_s, qw2_s, bw2_s, t2w_s, q2w_s, b2w_s
real, pointer :: uw_surf1, vw_surf1, tw_surf1, qw_surf1, t_surf1, q_surf1, uw_surf1_s, vw_surf1_s, tw_surf1_s, qw_surf1_s, t_surf1_s, q_surf1_s

real, dimension(:), allocatable, target     :: x_surf, x_surf_s
real, dimension(:,:), allocatable, target   :: x1, x2, x1_misc, x2_misc, x1_s, x2_s, x1_misc_s, x2_misc_s
real, dimension(:,:,:), allocatable, target :: xt, xt_press, x2t, x2t_press, xt_s, xt_press_s, x2t_s, x2t_press_s

contains

  !
  ! init_stat_slab
  !
  subroutine init_stat_slab
    use grid, only             : filprf, level, nzp, nxp, nyp, isfctyp
    use modstat_io, only       : create_stat_file, create_stat_var, create_stat_var_2d, zt_name, zm_name, ntypes
    use netcdf_interface, only : my_nf90_real, my_nf90_int, create_var
    use modutil_mpi, only      : par_sum
    implicit none

    !
    !
    !

    ! Initialize sample counts
    n_local  = 0
    n_global = 0

    ! Create statistics file
    call create_stat_file(trim(filprf)//'.stat.slab.nc', ncid, irec)
    call create_stat_var(ncid, 'n', zt_name, my_nf90_int)

    ! x_surf
    nstats_x_surf = 3
    if (level >= 1) nstats_x_surf = nstats_x_surf + 1
    if (isfctyp == 2 .or. isfctyp == 4) then
       nstats_x_surf = nstats_x_surf + 1
       if (level >= 1) nstats_x_surf = nstats_x_surf + 1
    end if
    allocate(x_surf(nstats_x_surf), x_surf_s(nstats_x_surf))
    x_surf(:) = 0.
    call create_stat_var(ncid, 'uw_surf', my_nf90_real)
    call create_stat_var(ncid, 'vw_surf', my_nf90_real)
    call create_stat_var(ncid, 'tw_surf', my_nf90_real)
    if (level >= 1) call create_stat_var(ncid, 'qw_surf', my_nf90_real)
    if (isfctyp == 2 .or. isfctyp == 4) then
       call create_stat_var(ncid, 't_surf', my_nf90_real)
       if (level >= 1) call create_stat_var(ncid, 'q_surf', my_nf90_real)
    end if

    ! x1
    nstats_x1 = 4
    if (level >= 1) nstats_x1 = nstats_x1 + 1
    allocate(x1(nzp,nstats_x1), x1_s(nzp,nstats_x1)) 
    x1(:,:) = 0.
    call create_stat_var(ncid, 'u', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'v', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'w', zm_name, my_nf90_real)
    call create_stat_var(ncid, 't', zt_name, my_nf90_real)
    if (level >= 1) call create_stat_var(ncid, 'q', zt_name, my_nf90_real)

    ! xt
    allocate(xt(nzp,ntypes,nstats_x1),xt_s(nzp,ntypes,nstats_x1)) 
    xt(:,:,:) = 0.
    call create_stat_var_2d(ncid, 'ut', zt_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'vt', zt_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'wt', zm_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'tt', zt_name, my_nf90_real)
    if (level >= 1) call create_stat_var_2d(ncid, 'qt', zt_name, my_nf90_real)

    ! xt_press
    allocate(xt_press(nzp,ntypes,nstats_x1),xt_press_s(nzp,ntypes,nstats_x1)) 
    xt_press(:,:,:) = 0.
    call create_stat_var_2d(ncid, 'ut_press', zt_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'vt_press', zt_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'wt_press', zm_name, my_nf90_real)

    ! x2
    nstats_x2 = 7
    if (level >= 1) nstats_x2 = nstats_x2 + 3
    allocate(x2(nzp,nstats_x2), x2_s(nzp,nstats_x2))
    x2(:,:) = 0.
    call create_stat_var(ncid, 'u2', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'v2', zt_name, my_nf90_real)
    call create_stat_var(ncid, 'w2', zm_name, my_nf90_real)
    call create_stat_var(ncid, 't2', zt_name, my_nf90_real)
    if (level >= 1) then
       call create_stat_var(ncid, 'q2', zt_name, my_nf90_real)
    end if
    call create_stat_var(ncid, 'uw',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 'vw',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 'tw',  zm_name, my_nf90_real)
    if (level >= 1) then
       call create_stat_var(ncid, 'qw', zm_name, my_nf90_real)
       call create_stat_var(ncid, 'tq', zt_name, my_nf90_real)
    end if

    ! x2t
    allocate(x2t(nzp,ntypes,nstats_x2), x2t_s(nzp,ntypes,nstats_x2))
    x2t(:,:,:) = 0.
    call create_stat_var_2d(ncid, 'u2t', zt_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'v2t', zt_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'w2t', zm_name, my_nf90_real)
    call create_stat_var_2d(ncid, 't2t', zt_name, my_nf90_real)
    if (level >= 1) then
       call create_stat_var_2d(ncid, 'q2t', zt_name, my_nf90_real)
    end if
    call create_stat_var_2d(ncid, 'uwt',  zm_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'vwt',  zm_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'twt',  zm_name, my_nf90_real)
    if (level >= 1) then
       call create_stat_var_2d(ncid, 'qwt', zm_name, my_nf90_real)
       call create_stat_var_2d(ncid, 'tqt', zt_name, my_nf90_real)
    end if

    ! x2t_press
    allocate(x2t_press(nzp,ntypes,nstats_x2), x2t_press_s(nzp,ntypes,nstats_x2))
    x2t_press(:,:,:) = 0.
    call create_stat_var_2d(ncid, 'u2t_press', zt_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'v2t_press', zt_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'w2t_press', zm_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'uwt_press',  zm_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'vwt_press',  zm_name, my_nf90_real)
    call create_stat_var_2d(ncid, 'twt_press',  zm_name, my_nf90_real)
    if (level >= 1) then
       call create_stat_var_2d(ncid, 'qwt_press', zm_name, my_nf90_real)
       call create_stat_var_2d(ncid, 'tqt_press', zt_name, my_nf90_real)
    end if

    ! x1_misc
    nstats_x1_misc = 10
    if (level >= 1) nstats_x1_misc = nstats_x1_misc + 3
    allocate(x1_misc(nzp,nstats_x1_misc), x1_misc_s(nzp,nstats_x1_misc)) 
    x1_misc(:,:) = 0.
    call create_stat_var(ncid, 'b',         zm_name, my_nf90_real)
    call create_stat_var(ncid, 'N2',        zm_name, my_nf90_real)
    call create_stat_var(ncid, 'pi1',       zt_name, my_nf90_real)
    call create_stat_var(ncid, 'wfls',      zt_name, my_nf90_real)
    call create_stat_var(ncid, 'dthldtls',  zt_name, my_nf90_real)
    call create_stat_var(ncid, 'e',         zm_name, my_nf90_real)
    call create_stat_var(ncid, 'et_shear',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 'et_diss',   zm_name, my_nf90_real)
    call create_stat_var(ncid, 'km',        zm_name, my_nf90_real)
    call create_stat_var(ncid, 'kh',        zm_name, my_nf90_real)
    if (level >= 2) then
       call create_stat_var(ncid, 'dqtdtls',  zt_name, my_nf90_real)
       call create_stat_var(ncid, 'th',       zt_name, my_nf90_real)
       call create_stat_var(ncid, 'l',        zt_name, my_nf90_real)
    end if

    ! x2_misc
    nstats_x2_misc = 7
    if (level >= 1) nstats_x2_misc = nstats_x2_misc + 2
    if (level >= 2) nstats_x2_misc = nstats_x2_misc + 2
    allocate(x2_misc(nzp,nstats_x2_misc), x2_misc_s(nzp,nstats_x2_misc))
    x2_misc(:,:) = 0.
    call create_stat_var(ncid, 'b2',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 'bw',  zm_name, my_nf90_real)
    if (level >= 2) then
       call create_stat_var(ncid, 'thb', zm_name, my_nf90_real)
       call create_stat_var(ncid, 'lb',  zm_name, my_nf90_real)
    end if
    call create_stat_var(ncid, 'w3',  zm_name, my_nf90_real)
    call create_stat_var(ncid, 'tw2', zm_name, my_nf90_real)
    if (level >= 1) call create_stat_var(ncid, 'qw2', zm_name, my_nf90_real)
    call create_stat_var(ncid, 'bw2', zm_name, my_nf90_real)
    call create_stat_var(ncid, 't2w', zm_name, my_nf90_real)
    if (level >= 1) call create_stat_var(ncid, 'q2w', zm_name, my_nf90_real)
    call create_stat_var(ncid, 'b2w', zm_name, my_nf90_real)

    !
    !
    !

    ! Sample counts
    n_s_local = (nxp - 4)*(nyp - 4)
    call par_sum(n_s_local, n_s_global)

    ! x_surf
    uw_surf1 => x_surf(1)
    vw_surf1 => x_surf(2)
    tw_surf1 => x_surf(3)
    if (level >= 1) qw_surf1 => x_surf(4)
    if (isfctyp == 2 .or. isfctyp == 4) then
       t_surf1 => x_surf(5)
       if (level >= 1) q_surf1 => x_surf(6)
    end if

    ! x_surf_s
    uw_surf1_s => x_surf_s(1)
    vw_surf1_s => x_surf_s(2)
    tw_surf1_s => x_surf_s(3)
    if (level >= 1) qw_surf1_s => x_surf_s(4)
    if (isfctyp == 2 .or. isfctyp == 4) then
       t_surf1_s => x_surf_s(5)
       if (level >= 1) q_surf1_s => x_surf_s(6)
    end if

    ! x1
    u1 => x1(:,1)
    v1 => x1(:,2)
    w1 => x1(:,3)
    t1 => x1(:,4)
    if (level >= 1) q1 => x1(:,5)

    ! x1_s
    u1_s => x1_s(:,1)
    v1_s => x1_s(:,2)
    w1_s => x1_s(:,3)
    t1_s => x1_s(:,4)
    if (level >= 1) q1_s => x1_s(:,5)

    ! xt_press
    ut_press => xt_press(:,:,1)
    vt_press => xt_press(:,:,2)
    wt_press => xt_press(:,:,3)
    tt_press => xt_press(:,:,4)
    if (level >= 1) qt_press => xt_press(:,:,5)

    ! xt_press_s
    ut_press_s => xt_press_s(:,:,1)
    vt_press_s => xt_press_s(:,:,2)
    wt_press_s => xt_press_s(:,:,3)
    tt_press_s => xt_press_s(:,:,4)
    if (level >= 1) qt_press_s => xt_press_s(:,:,5)

    ! xt
    ut => xt(:,:,1)
    vt => xt(:,:,2)
    wt => xt(:,:,3)
    tt => xt(:,:,4)
    if (level >= 1) qt => xt(:,:,5)

    ! xt_s
    ut_s => xt_s(:,:,1)
    vt_s => xt_s(:,:,2)
    wt_s => xt_s(:,:,3)
    tt_s => xt_s(:,:,4)
    if (level >= 1) qt_s => xt_s(:,:,5)

    ! x2
    u2 => x2(:,1)
    v2 => x2(:,2)
    w2 => x2(:,3)
    t2 => x2(:,4)
    uw => x2(:,5)
    vw => x2(:,6)
    tw => x2(:,7)
    if (level >= 1) then
       q2 => x2(:,8)
       qw => x2(:,9)
       tq => x2(:,10)
    end if

    ! x2_s
    u2_s => x2_s(:,1)
    v2_s => x2_s(:,2)
    w2_s => x2_s(:,3)
    t2_s => x2_s(:,4)
    uw_s => x2_s(:,5)
    vw_s => x2_s(:,6)
    tw_s => x2_s(:,7)
    if (level >= 1) then
       q2_s => x2_s(:,8)
       qw_s => x2_s(:,9)
       tq_s => x2_s(:,10)
    end if

    ! x2t
    u2t => x2t(:,:,1)
    v2t => x2t(:,:,2)
    w2t => x2t(:,:,3)
    t2t => x2t(:,:,4)
    uwt => x2t(:,:,5)
    vwt => x2t(:,:,6)
    twt => x2t(:,:,7)
    if (level >= 1) then
       q2t => x2t(:,:,8)
       qwt => x2t(:,:,9)
       tqt => x2t(:,:,10)
    end if

    ! x2t_s
    u2t_s => x2t_s(:,:,1)
    v2t_s => x2t_s(:,:,2)
    w2t_s => x2t_s(:,:,3)
    t2t_s => x2t_s(:,:,4)
    uwt_s => x2t_s(:,:,5)
    vwt_s => x2t_s(:,:,6)
    twt_s => x2t_s(:,:,7)
    if (level >= 1) then
       q2t_s => x2t_s(:,:,8)
       qwt_s => x2t_s(:,:,9)
       tqt_s => x2t_s(:,:,10)
    end if

    ! x2t_press
    u2t_press => x2t_press(:,:,1)
    v2t_press => x2t_press(:,:,2)
    w2t_press => x2t_press(:,:,3)
    t2t_press => x2t_press(:,:,4)
    uwt_press => x2t_press(:,:,5)
    vwt_press => x2t_press(:,:,6)
    twt_press => x2t_press(:,:,7)
    if (level >= 1) then
       q2t_press => x2t_press(:,:,8)
       qwt_press => x2t_press(:,:,9)
       tqt_press => x2t_press(:,:,10)
    end if

    ! x2t_press_s
    u2t_press_s => x2t_press_s(:,:,1)
    v2t_press_s => x2t_press_s(:,:,2)
    w2t_press_s => x2t_press_s(:,:,3)
    t2t_press_s => x2t_press_s(:,:,4)
    uwt_press_s => x2t_press_s(:,:,5)
    vwt_press_s => x2t_press_s(:,:,6)
    twt_press_s => x2t_press_s(:,:,7)
    if (level >= 1) then
       q2t_press_s => x2t_press_s(:,:,8)
       qwt_press_s => x2t_press_s(:,:,9)
       tqt_press_s => x2t_press_s(:,:,10)
    end if

    ! x1_misc
    b1        => x1_misc(:,1)
    N2        => x1_misc(:,2)
    pi11      => x1_misc(:,3)
    wfls1     => x1_misc(:,4)
    dthldtls1 => x1_misc(:,5)
    e         => x1_misc(:,6)
    et_shear  => x1_misc(:,7)
    et_diss   => x1_misc(:,8)
    km        => x1_misc(:,9)
    kh        => x1_misc(:,10)
    if (level >= 2) then
       dqtdtls1 => x1_misc(:,11)
       th1      => x1_misc(:,12)
       l1       => x1_misc(:,13)
    end if

    ! x1_misc_s
    b1_s        => x1_misc_s(:,1)
    N2_s        => x1_misc_s(:,2)
    pi11_s      => x1_misc_s(:,3)
    wfls1_s     => x1_misc_s(:,4)
    dthldtls1_s => x1_misc_s(:,5)
    e_s         => x1_misc_s(:,6)
    et_shear_s  => x1_misc_s(:,7)
    et_diss_s   => x1_misc_s(:,8)
    km_s        => x1_misc_s(:,9)
    kh_s        => x1_misc_s(:,10)
    if (level >= 2) then
       dqtdtls1_s => x1_misc_s(:,11)
       th1_s      => x1_misc_s(:,12)
       l1_s       => x1_misc_s(:,13) 
    end if

    ! x2_misc
    b2  => x2_misc(:,1)
    bw  => x2_misc(:,2)
    w3  => x2_misc(:,3)
    tw2 => x2_misc(:,4)
    bw2 => x2_misc(:,5)
    t2w => x2_misc(:,6)
    b2w => x2_misc(:,7)
    if (level >= 1) then
       qw2 => x2_misc(:,8)
       q2w => x2_misc(:,9)
    end if
    if (level >= 2) then
       thb => x2_misc(:,10)
       lb  => x2_misc(:,11)
    end if

    ! x2_misc_s
    b2_s  => x2_misc_s(:,1)
    bw_s  => x2_misc_s(:,2)
    w3_s  => x2_misc_s(:,3)
    tw2_s => x2_misc_s(:,4)
    bw2_s => x2_misc_s(:,5)
    t2w_s => x2_misc_s(:,6)
    b2w_s => x2_misc_s(:,7)
    if (level >= 1) then
       qw2_s => x2_misc_s(:,8)
       q2w_s => x2_misc_s(:,9)
    end if
    if (level >= 1) then
       thb_s => x2_misc_s(:,10)
       lb_s  => x2_misc_s(:,11)
    end if
    
  end subroutine init_stat_slab

  !
  ! sample_stat_slab_surface
  !
  subroutine sample_stat_slab_surface(sst)
    use grid, only        : a_up, a_vp, a_wp, a_tp, a_rp, liquid, nzp, nyp, nxp, level, a_theta, vapor, &
                            zm, a_pexnr, zt, a_ut, a_vt, a_wt, a_tt, a_rt, th00, dxi, dyi, &
                            uw_sfc, vw_sfc, wt_sfc, wq_sfc, psrf, isfctyp
    use modutil_mpi, only : par_sum
    use thrm, only        : rslf 
    use defs, only        : p00, rcp
    implicit none

    real, intent(in) :: sst

    real, allocatable, dimension(:), target :: x_surf_scratch
    integer :: i, j

    ! Sample surface fluxes
    x_surf_s(:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       uw_surf1_s = uw_surf1_s + uw_sfc(i,j)
       vw_surf1_s = vw_surf1_s + vw_sfc(i,j)
       tw_surf1_s = tw_surf1_s + wt_sfc(i,j)
    end do
    end do
    allocate(x_surf_scratch(nstats_x_surf))
    call par_sum(x_surf_s, x_surf_scratch, nstats_x_surf)
    x_surf_s(:) = x_surf_scratch(:)/real(n_s_global)

    ! Sample first-order statistics at the surface
    if (isfctyp == 2 .or. isfctyp == 4) then
       t_surf1_s = sst/(psrf/p00)**rcp
       if (level >= 1) q_surf1_s = rslf(psrf, sst)
    end if

  end subroutine sample_stat_slab_surface

  !
  ! sample_stat_slab
  !
  subroutine sample_stat_slab
    use grid, only        : a_up, a_vp, a_wp, a_tp, a_rp, liquid, nzp, nyp, nxp, level, a_theta, vapor, &
                            zm, a_pexnr, zt, a_ut, a_vt, a_wt, a_tt, a_rt, th00, dxi, dyi
    use modutil_mpi, only : par_sum
    implicit none

    real, pointer, dimension(:) :: u1_local, v1_local, w1_local, t1_local, q1_local
    real, allocatable, dimension(:,:), target :: x1_scratch
    integer :: i, j, k, kp1, ip1, jp1

    ! Sample first-order statistics
    x1_s(:,:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 1,nzp
          u1_s(k) = u1_s(k) + a_up(k,i,j)
          v1_s(k) = v1_s(k) + a_vp(k,i,j)
          w1_s(k) = w1_s(k) + a_wp(k,i,j)
          t1_s(k) = t1_s(k) + a_tp(k,i,j)
          if (level >= 1) q1_s(k) = q1_s(k) + a_rp(k,i,j)
       end do
    end do
    end do
    allocate(x1_scratch(nzp,nstats_x1))
    call par_sum(x1_s, x1_scratch, nzp*nstats_x1)
    x1_s(:,:) = x1_scratch(:,:)/real(n_s_global)

    ! Sample second-order statistics
    x2_s(:,:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 2,nzp-1
          kp1 = k+1
          ip1 = i+1
          jp1 = j+1
          
          u2_s(k) = u2_s(k) + (a_up(k,i,j) - u1_s(k))**2.
          v2_s(k) = v2_s(k) + (a_vp(k,i,j) - v1_s(k))**2.
          w2_s(k) = w2_s(k) + (a_wp(k,i,j) - w1_s(k))**2.
          t2_s(k) = t2_s(k) + (a_tp(k,i,j) - t1_s(k))**2.
          uw_s(k) = uw_s(k) + (a_wp(k,i,j) - w1_s(k))*(0.25*(a_up(k,i,j)+a_up(k,ip1,j)+a_up(kp1,i,j)+a_up(kp1,ip1,j)) &
                                                       - 0.5*(u1_s(k)+u1_s(kp1)))
          vw_s(k) = vw_s(k) + (a_wp(k,i,j) - w1_s(k))*(0.25*(a_vp(k,i,j)+a_vp(k,i,jp1)+a_vp(kp1,i,j)+a_vp(kp1,i,jp1)) &
                                                       - 0.5*(v1_s(k)+v1_s(kp1)))
          tw_s(k) = tw_s(k) + (a_wp(k,i,j) - w1_s(k))*(0.5*(a_tp(k,i,j)+a_tp(kp1,i,j)) - 0.5*(t1_s(k)+t1_s(kp1)))
          if (level >= 1) then
             q2_s(k) = q2_s(k) + (a_rp(k,i,j) - q1_s(k))**2.
             tq_s(k) = tq_s(k) + (a_tp(k,i,j) - t1_s(k))*(a_rp(k,i,j) - q1_s(k))
             qw_s(k) = qw_s(k) + (a_wp(k,i,j) - w1_s(k))*(0.5*(a_rp(k,i,j)+a_rp(kp1,i,j)) - 0.5*(q1_s(k)+q1_s(kp1)))
          end if
       end do
    end do
    end do
    x2_s(:,:) = x2_s(:,:)/real(n_s_local)

  end subroutine sample_stat_slab

  !
  ! stat_slab_tendency
  !
  subroutine stat_slab_tendency(itype)
    use grid, only        : a_up, a_vp, a_wp, a_tp, a_rp, liquid, nzp, nyp, nxp, level, a_theta, vapor, &
                            zm, a_pexnr, zt, a_ut, a_vt, a_wt, a_tt, a_rt, th00, dxi, dyi
    use modutil_mpi, only : par_sum
    use modstat_io, only  : ntypes

    implicit none

    integer, intent(in) :: itype

    real, pointer, dimension(:) :: ut_local, vt_local, wt_local, tt_local, qt_local
    real, allocatable, dimension(:,:), target :: xt_local, xt_scratch
    integer :: i, j, k, kp1, ip1, jp1

    allocate(xt_local(nzp,nstats_x1))
    ut_local => xt_local(:,1)
    vt_local => xt_local(:,2)
    wt_local => xt_local(:,3)
    tt_local => xt_local(:,4)
    if (level >= 1) qt_local => xt_local(:,5)
    
    ! Sample first-order statistics
    xt_local(:,:)  = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 1,nzp
          ut_local(k) = ut_local(k) + a_ut(k,i,j)
          vt_local(k) = vt_local(k) + a_vt(k,i,j)
          wt_local(k) = wt_local(k) + a_wt(k,i,j)
          tt_local(k) = tt_local(k) + a_tt(k,i,j)
          if (level >= 1) qt_local(k) = qt_local(k) + a_rt(k,i,j)
       end do
    end do
    end do
    allocate(xt_scratch(nzp,nstats_x1))
    call par_sum(xt_local, xt_scratch, nzp*nstats_x1)
    xt_local(:,:) = xt_scratch(:,:)/real(n_s_global)

    !
    if (itype == 1) xt_s(:,:,:) = 0.
    ut_s(:,itype) = ut_local(:)
    vt_s(:,itype) = vt_local(:)
    wt_s(:,itype) = wt_local(:)
    tt_s(:,itype) = tt_local(:)
    if (level >=1 ) qt_s(:,itype) = qt_local(:)

    ! Sample second-order statistics
    if (itype == 1) x2t_s(:,:,:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 2,nzp-1
          kp1 = k+1
          ip1 = i+1
          jp1 = j+1
          
          u2t_s(k,itype) = u2t_s(k,itype) + 2.*(a_up(k,i,j) - u1_s(k))*(a_ut(k,i,j) - ut_s(k,itype))
          v2t_s(k,itype) = v2t_s(k,itype) + 2.*(a_vp(k,i,j) - v1_s(k))*(a_vt(k,i,j) - vt_s(k,itype))
          w2t_s(k,itype) = w2t_s(k,itype) + 2.*(a_wp(k,i,j) - w1_s(k))*(a_wt(k,i,j) - wt_s(k,itype))
          t2t_s(k,itype) = t2t_s(k,itype) + 2.*(a_tp(k,i,j) - t1_s(k))*(a_tt(k,i,j) - tt_s(k,itype))
          uwt_s(k,itype) = uwt_s(k,itype) + (a_wp(k,i,j) - w1_s(k))*(0.25*(a_ut(k,i,j)+a_ut(k,ip1,j)+a_ut(kp1,i,j)+a_ut(kp1,ip1,j))     &
                                                                     - 0.5*(ut_s(k,itype)+ut_s(kp1,itype)))                             &
                                          + (a_wt(k,i,j) - wt_s(k,itype))*(0.25*(a_up(k,i,j)+a_up(k,ip1,j)+a_up(kp1,i,j)+a_up(kp1,ip1,j)) &
                                                                           - 0.5*(u1_s(k)+u1_s(kp1)))
          vwt_s(k,itype) = vwt_s(k,itype) + (a_wp(k,i,j) - w1_s(k))*(0.25*(a_vt(k,i,j)+a_vt(k,i,jp1)+a_vt(kp1,i,j)+a_vt(kp1,i,jp1))     &
                                                                     - 0.5*(vt_s(k,itype)+vt_s(kp1,itype)))                             &
                                          + (a_wt(k,i,j) - wt_s(k,itype))*(0.25*(a_vp(k,i,j)+a_vp(k,i,jp1)+a_vp(kp1,i,j)+a_vp(kp1,i,jp1)) &
                                                                           - 0.5*(v1_s(k)+v1_s(kp1)))                                
          twt_s(k,itype) = twt_s(k,itype) + (a_wp(k,i,j) - w1_s(k))*(0.5*(a_tt(k,i,j)+a_tt(kp1,i,j)) - 0.5*(tt_s(k,itype)+tt_s(kp1,itype))) &
                                          + (a_wt(k,i,j) - wt_s(k,itype))*(0.5*(a_tp(k,i,j)+a_tp(kp1,i,j)) - 0.5*(t1_s(k)+t1_s(kp1)))
          if (level >= 1) then
             q2t_s(k,itype) = q2t_s(k,itype) + 2.*(a_rp(k,i,j) - q1_s(k))*(a_rt(k,i,j) - qt_s(k,itype))
             tqt_s(k,itype) = tqt_s(k,itype) + (a_tp(k,i,j) - t1_s(k))*(a_rt(k,i,j) - qt_s(k,itype)) &
                                             + (a_tt(k,i,j) - tt_s(k,itype))*(a_rp(k,i,j) - q1_s(k))
             qwt_s(k,itype) = qwt_s(k,itype) + (a_wp(k,i,j) - w1_s(k))*(0.5*(a_rt(k,i,j)+a_rt(kp1,i,j))        &
                                                                        - 0.5*(qt_s(k,itype)+qt_s(kp1,itype))) &
                                             + (a_wt(k,i,j) - wt_s(k,itype))*(0.5*(a_rp(k,i,j)+a_rp(kp1,i,j)) - 0.5*(q1_s(k)+q1_s(kp1)))
          end if
       end do
    end do
    end do

    !
    x2t_s(:,itype,:) = x2t_s(:,itype,:)/real(n_s_local)

  end subroutine stat_slab_tendency

  !
  ! stat_slab_pressure
  !
  subroutine stat_slab_pressure(itype, u, v, w, pp)
    use grid, only        : a_ut, a_vt, a_wt, a_tp, a_rp, liquid, nzp, nyp, nxp, level, a_theta, vapor, &
                            zm, zt, th00, dxi, dyi, a_tt, a_rt
    use modutil_mpi, only : par_sum
    use modstat_io, only  : ntypes

    implicit none

    integer, intent(in) :: itype
    real, intent(in), dimension(:,:,:) :: u, v, w, pp

    real :: Dp_u, Dp_v, Dp_w
    real, pointer, dimension(:) :: ut_local, vt_local, wt_local
    real, allocatable, dimension(:,:), target :: xt_local, xt_scratch
    integer :: i, j, k, kp1, ip1, jp1

    allocate(xt_local(nzp,nstats_x1))
    ut_local => xt_local(:,1)
    vt_local => xt_local(:,2)
    wt_local => xt_local(:,3)
    
    ! Sample first-order statistics
    xt_local(:,:)  = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 2,nzp-1
          kp1 = k+1
          ip1 = i+1
          jp1 = j+1

          ut_local(k) = ut_local(k) - th00*(pp(k,ip1,j) - pp(k,i,j))*dxi
          vt_local(k) = vt_local(k) - th00*(pp(k,i,jp1) - pp(k,i,j))*dyi
          wt_local(k) = wt_local(k) - th00*(pp(kp1,i,j) - pp(k,i,j))/(zt(kp1) - zt(k))
       end do
    end do
    end do
    allocate(xt_scratch(nzp,nstats_x1))
    call par_sum(xt_local, xt_scratch, nzp*nstats_x1)
    xt_local(:,:) = xt_scratch(:,:)/real(n_s_global)

    !
    if (itype == 1) xt_press_s(:,:,:) = 0.
    ut_press_s(:,itype) = ut_local(:)
    vt_press_s(:,itype) = vt_local(:)
    wt_press_s(:,itype) = wt_local(:)
    tt_press_s(:,itype) = 0.
    if (level >= 1) qt_press_s(:,itype) = 0.

    ! Sample second-order statistics
    if (itype == 1) x2t_press_s(:,:,:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 2,nzp-1
          kp1 = k+1
          ip1 = i+1
          jp1 = j+1

          Dp_u = th00*(pp(k,ip1,j) - pp(k,i,j))*dxi
          Dp_v = th00*(pp(k,i,jp1) - pp(k,i,j))*dyi
          Dp_w = th00*(pp(kp1,i,j) - pp(k,i,j))/(zt(kp1) - zt(k))          

          u2t_press_s(k,itype) = u2t_press_s(k,itype) + 2.*(u(k,i,j) - u1_s(k))*(-Dp_u - ut_press_s(k,itype))
          v2t_press_s(k,itype) = v2t_press_s(k,itype) + 2.*(v(k,i,j) - v1_s(k))*(-Dp_v - vt_press_s(k,itype))
          w2t_press_s(k,itype) = w2t_press_s(k,itype) + 2.*(w(k,i,j) - w1_s(k))*(-Dp_w - wt_press_s(k,itype))
          twt_press_s(k,itype) = twt_press_s(k,itype) &
                               + (-Dp_w - wt_press_s(k,itype))*(0.5*(a_tp(k,i,j)+a_tp(kp1,i,j)) - 0.5*(t1_s(k)+t1_s(kp1)))
          if (level >= 1) then
             qwt_press_s(k,itype) = qwt_press_s(k,itype) + (-Dp_w - wt_press_s(k,itype))*(0.5*(a_rp(k,i,j)+a_rp(kp1,i,j)) - 0.5*(q1_s(k)+q1_s(kp1)))
          end if
       end do
    end do
    end do

    ! 
    x2t_press_s(:,itype,:) = x2t_press_s(:,itype,:)/real(n_s_local)

  end subroutine stat_slab_pressure

  !
  ! stat_slab_sfs
  !
  subroutine stat_slab_sfs(csx, delta, pr, km, kh, ri)
    use defs, only : pi
    use grid, only : nxp, nyp, nzp, dn0, a_tp, a_rp, zt, level
    implicit none

    real, intent(in)                   :: csx, delta, pr
    real, dimension(:,:,:), intent(in) :: km, kh, ri

    integer :: i, j, k
    real :: km_temp

    x1_misc_s(:,:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 2,nzp-1
          km_temp = km(k,i,j)/(0.5*(dn0(k)+dn0(k+1)))

          e_s(k)        = e_s(k)        + km_temp**2./(delta*pi*(csx**2))**2
          et_shear_s(k) = et_shear_s(k) + km_temp*kh(k,i,j)/(1. - ri(k,i,j)/pr)
          et_diss_s(k)  = et_diss_s(k)  + km_temp*kh(k,i,j)
          km_s(k)       = km_s(k)       + km_temp
          kh_s(k)       = kh_s(k)       + km_temp/pr
       end do
    end do
    end do

  end subroutine stat_slab_sfs

  !
  ! stat_slab_misc
  !
  subroutine stat_slab_misc
    use grid, only        : a_up, a_vp, a_wp, a_tp, a_rp, liquid, nzp, nyp, nxp, level, a_theta, vapor, &
                            zm, a_pexnr, zt, a_ut, a_vt, a_wt, a_tt, a_rt, th00, dxi, dyi, a_scr1, pi1, &
                            wfls, dthldtls, dqtdtls
    use modutil_mpi, only : par_sum
    use defs, only        : cp
    implicit none

    real, allocatable, dimension(:,:), target :: x1_misc_scratch
    integer :: i, j, k, kp1

    ! Sample first-order statistics
!    x1_misc_s(:,:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 1,nzp
          b1_s(k) = b1_s(k) + 0.5*(a_scr1(k,i,j)+a_scr1(k+1,i,j))
          N2_s(k) = N2_s(k) + a_scr1(k,i,j) ! temporary
          if (level >= 2) then
             th1_s(k) = th1_s(k) + a_theta(k,i,j)
             l1_s(k)  = l1_s(k)  + liquid(k,i,j)
          end if
       end do
    end do
    end do
    allocate(x1_misc_scratch(nzp,nstats_x1_misc))
    call par_sum(x1_misc_s, x1_misc_scratch, nzp*nstats_x1_misc)
    x1_misc_s(:,:) = x1_misc_scratch(:,:)/real(n_s_global)

    !
    pi11_s(:)      = pi1(:)/cp
    wfls1_s(:)     = wfls(:)
    dthldtls1_s(:) = dthldtls(:)
    dqtdtls1_s(:)  = dqtdtls(:)

    ! Compute N**2. from buoyancy profile temporarily stored in N2
    N2_s(2:nzp-1) = (N2_s(3:nzp) - N2_s(2:nzp-1))/(zm(3:nzp) - zm(2:nzp-1))

    ! Sample second-order statistics
    x2_misc_s(:,:) = 0.
    do j = 3,nyp-2
    do i = 3,nxp-2
       do k = 2,nzp-1
          kp1 = k+1
          
          b2_s(k)  = b2_s(k)  + (0.5*(a_scr1(k,i,j)+a_scr1(kp1,i,j))  - b1_s(k))**2.
          bw_s(k)  = bw_s(k)  + (0.5*(a_scr1(k,i,j)+a_scr1(kp1,i,j))  - b1_s(k))*(a_wp(k,i,j) - w1_s(k))
          w3_s(k)  = w3_s(k)  + (a_wp(k,i,j) - w1_s(k))**3.
          tw2_s(k) = tw2_s(k) + (a_wp(k,i,j) - w1_s(k))**2.*(0.5*(a_tp(k,i,j)+a_tp(kp1,i,j)) - 0.5*(t1_s(k)+t1_s(kp1)))
          bw2_s(k) = bw2_s(k) + (a_wp(k,i,j) - w1_s(k))**2.*(0.5*(a_scr1(k,i,j)+a_scr1(kp1,i,j)) - b1_s(k))
          t2w_s(k) = t2w_s(k) + (a_wp(k,i,j) - w1_s(k))*(0.5*(a_tp(k,i,j)+a_tp(kp1,i,j)) - 0.5*(t1_s(k)+t1_s(kp1)))**2.
          bw2_s(k) = bw2_s(k) + (a_wp(k,i,j) - w1_s(k))*(0.5*(a_scr1(k,i,j)+a_scr1(kp1,i,j)) - b1_s(k))**2.
          if (level >= 1) then
             qw2_s(k) = qw2_s(k) + (a_wp(k,i,j) - w1_s(k))**2.*(0.5*(a_rp(k,i,j)+a_rp(kp1,i,j)) - 0.5*(q1_s(k)+q1_s(kp1)))
             q2w_s(k) = q2w_s(k) + (a_wp(k,i,j) - w1_s(k))*(0.5*(a_rp(k,i,j)+a_rp(kp1,i,j)) - 0.5*(q1_s(k)+q1_s(kp1)))**2.
          end if
          if (level >= 2) then
             thb_s(k) = thb_s(k) + (0.5*(a_theta(k,i,j)+a_theta(kp1,i,j)) - 0.5*(th1_s(k)+th1_s(kp1)))&
                                 * (0.5*(a_scr1(k,i,j)+a_scr1(kp1,i,j)) - b1_s(k))
             lb_s(k)  = lb_s(k)  + (0.5*(liquid(k,i,j)+liquid(kp1,i,j)) - 0.5*(l1_s(k)+l1_s(kp1)))&
                                 * (0.5*(a_scr1(k,i,j)+a_scr1(kp1,i,j)) - b1_s(k))
          end if
       end do
    end do
    end do
    x2_misc_s(:,:) = x2_misc_s(:,:)/real(n_s_local)

  end subroutine stat_slab_misc

  !
  ! update_stat_slab
  !
  subroutine update_stat_slab
    use grid, only        : a_up, a_vp, a_wp, a_tp, a_rp, liquid, nzp, nyp, nxp, level, a_theta, vapor, &
                            zm, a_pexnr, zt, a_ut, a_vt, a_wt, a_tt, a_rt, th00, dxi, dyi, isfctyp
    use modutil_mpi, only : par_sum
    use modstat_io, only  : ntypes
    implicit none

    real :: f
    integer :: itype

    !
    f = real(n_s_local)/real(n_s_local + n_local)

    ! Update third-order statistics
    w3(:)  = f*w3_s(:)  + (1. - f)*w3(:)  
    b2w(:) = f*b2w_s(:) + (1. - f)*b2w(:) + 2.*f*(1. - f)*(b1_s(:) - b1(:))*(bw_s(:) - bw(:))
    t2w(2:nzp-1) = f*t2w_s(2:nzp-1) + (1. - f)*t2w(2:nzp-1) &
                 + 2.*f*(1. - f)*(0.5*(t1_s(2:nzp-1)+t1_s(3:nzp)) - 0.5*(t1(2:nzp-1)+t1(3:nzp)))*(tw_s(2:nzp-1) - tw(2:nzp-1))
    bw2(:) = f*bw2_s(:) + (1. - f)*bw2(:) + f*(1. - f)*(b1_s(:) - b1(:))*(w2_s(:) - w2(:))                                           
    tw2(2:nzp-1) = f*tw2_s(2:nzp-1) + (1. - f)*tw2(2:nzp-1) &
                 + f*(1. - f)*(0.5*(t1_s(2:nzp-1)+t1_s(3:nzp)) - 0.5*(t1(2:nzp-1)+t1(3:nzp)))*(w2_s(2:nzp-1) - w2(2:nzp-1))
    if (level >= 1) then
       q2w(2:nzp-1) = f*q2w_s(2:nzp-1) + (1. - f)*q2w(2:nzp-1) &
                    + 2.*f*(1. - f)*(0.5*(q1_s(2:nzp-1)+q1_s(3:nzp)) - 0.5*(q1(2:nzp-1)+q1(3:nzp)))*(qw_s(2:nzp-1) - qw(2:nzp-1))
       qw2(2:nzp-1) = f*qw2_s(2:nzp-1) + (1. - f)*qw2(2:nzp-1) &
                    + f*(1. - f)*(0.5*(q1_s(2:nzp-1)+q1_s(3:nzp)) - 0.5*(q1(2:nzp-1)+q1(3:nzp)))*(w2_s(2:nzp-1) - w2(2:nzp-1))
    end if

    ! Update second-order diagnostic variable statistics
    b2(:)  = f*b2_s(:)  + (1. - f)*b2(:)  + f*(1. - f)*(b1_s(:) - b1(:))**2.
    bw(:)  = f*bw_s(:)  + (1. - f)*bw(:)  
    thb(:) = f*thb_s(:) + (1. - f)*thb(:) + f*(1. - f)*(th1_s(:) - th1(:))*(b1_s(:) - b1(:))
    lb(:)  = f*lb_s(:)  + (1. - f)*lb(:)  + f*(1. - f)*(l1_s(:) - l1(:))*(b1_s(:) - b1(:))

    ! Update second-order prognostic variable statistics
    u2(:) = f*u2_s(:) + (1. - f)*u2(:) + f*(1. - f)*(u1_s(:) - u1(:))**2.
    v2(:) = f*v2_s(:) + (1. - f)*v2(:) + f*(1. - f)*(v1_s(:) - v1(:))**2.
    w2(:) = f*w2_s(:) + (1. - f)*w2(:) + f*(1. - f)*(w1_s(:) - w1(:))**2.
    t2(:) = f*t2_s(:) + (1. - f)*t2(:) + f*(1. - f)*(t1_s(:) - t1(:))**2.
    uw(:) = f*uw_s(:) + (1. - f)*uw(:) 
    vw(:) = f*vw_s(:) + (1. - f)*vw(:) 
    tw(:) = f*tw_s(:) + (1. - f)*tw(:) 
    if (level >= 1) then
       q2(:) = f*q2_s(:) + (1. - f)*q2(:) + f*(1. - f)*(q1_s(:) - q1(:))**2.
       qw(:) = f*qw_s(:) + (1. - f)*qw(:) 
       tq(:) = f*tq_s(:) + (1. - f)*tq(:) + f*(1. - f)*(t1_s(:) - t1(:))*(q1_s(:) - q1(:))
    end if

    ! Update second-order tendency statistics
    do itype = 1,ntypes
       u2t(:,itype) = f*u2t_s(:,itype) + (1. - f)*u2t(:,itype) + 2.*f*(1. - f)*(u1_s(:) - u1(:))*(ut_s(:,itype) - ut(:,itype))
       v2t(:,itype) = f*v2t_s(:,itype) + (1. - f)*v2t(:,itype) + 2.*f*(1. - f)*(v1_s(:) - v1(:))*(vt_s(:,itype) - vt(:,itype))
       w2t(:,itype) = f*w2t_s(:,itype) + (1. - f)*w2t(:,itype) + 2.*f*(1. - f)*(w1_s(:) - w1(:))*(wt_s(:,itype) - wt(:,itype))
       t2t(:,itype) = f*t2t_s(:,itype) + (1. - f)*t2t(:,itype) + 2.*f*(1. - f)*(t1_s(:) - t1(:))*(tt_s(:,itype) - tt(:,itype))
       uwt(2:nzp-1,itype) = f*uwt_s(2:nzp-2,itype) + (1. - f)*uwt(2:nzp-1,itype) & 
                          + f*(1. - f)*(wt_s(2:nzp-1,itype) - wt(2:nzp-1,itype))&
                            *(0.5*(u1_s(2:nzp-1)+u1_s(3:nzp)) - 0.5*(u1(2:nzp-1)+u1(3:nzp)))
       vwt(2:nzp-1,itype) = f*vwt_s(2:nzp-2,itype) + (1. - f)*vwt(2:nzp-1,itype) & 
                          + f*(1. - f)*(wt_s(2:nzp-1,itype) - wt(2:nzp-1,itype))&
                            *(0.5*(v1_s(2:nzp-1)+v1_s(3:nzp)) - 0.5*(v1(2:nzp-1)+v1(3:nzp)))
       twt(2:nzp-1,itype) = f*twt_s(2:nzp-2,itype) + (1. - f)*twt(2:nzp-1,itype) & 
                          + f*(1. - f)*(wt_s(2:nzp-1,itype) - wt(2:nzp-1,itype))&
                            *(0.5*(t1_s(2:nzp-1)+t1_s(3:nzp)) - 0.5*(t1(2:nzp-1)+t1(3:nzp)))
       if (level >= 1) then
          q2t(:,itype) = f*q2t_s(:,itype) + (1. - f)*q2t(:,itype) + 2.*f*(1. - f)*(q1_s(:) - q1(:))*(qt_s(:,itype) - qt(:,itype))
          qwt(2:nzp-1,itype) = f*qwt_s(2:nzp-2,itype) + (1. - f)*qwt(2:nzp-1,itype) & 
                             + f*(1. - f)*(wt_s(2:nzp-1,itype) - wt(2:nzp-1,itype))&
                               *(0.5*(q1_s(2:nzp-1)+q1_s(3:nzp)) - 0.5*(q1(2:nzp-1)+q1(3:nzp)))
          tqt(:,itype) = f*tqt_s(:,itype) + (1. - f)*tqt(:,itype) + f*(1. - f)*(t1_s(:) - t1(:))*(qt_s(:,itype) - qt(:,itype)) &
                                                                  + f*(1. - f)*(tt_s(:,itype) - tt(:,itype))*(q1_s(:) - q1(:))
       end if
    end do

    !
    do itype = 1,ntypes
       u2t_press(:,itype) = f*u2t_press_s(:,itype) + (1. - f)*u2t_press(:,itype) &
                          + 2.*f*(1. - f)*(u1_s(:) - u1(:))*(ut_press_s(:,itype) - ut_press(:,itype))
       v2t_press(:,itype) = f*v2t_press_s(:,itype) + (1. - f)*v2t_press(:,itype) &
                          + 2.*f*(1. - f)*(v1_s(:) - v1(:))*(vt_press_s(:,itype) - vt_press(:,itype))
       w2t_press(:,itype) = f*w2t_press_s(:,itype) + (1. - f)*w2t_press(:,itype) &
                          + 2.*f*(1. - f)*(w1_s(:) - w1(:))*(wt_press_s(:,itype) - wt_press(:,itype))
       uwt_press(2:nzp-1,itype) = f*uwt_press_s(2:nzp-1,itype) + (1. - f)*uwt_press(2:nzp-1,itype) &
                                + f*(1. - f)*(wt_press_s(2:nzp-1,itype) - wt_press(2:nzp-1,itype))&
                                  *(0.5*(u1_s(2:nzp-1)+u1_s(3:nzp)) - 0.5*(u1(2:nzp-1)+u1(3:nzp)))
       vwt_press(2:nzp-1,itype) = f*vwt_press_s(2:nzp-1,itype) + (1. - f)*vwt_press(2:nzp-1,itype) &
                                + f*(1. - f)*(wt_press_s(2:nzp-1,itype) - wt_press(2:nzp-1,itype))&
                                  *(0.5*(v1_s(2:nzp-1)+v1_s(3:nzp)) - 0.5*(v1(2:nzp-1)+v1(3:nzp)))
       twt_press(2:nzp-1,itype) = f*twt_press_s(2:nzp-1,itype) + (1. - f)*twt_press(2:nzp-1,itype) &
                                + f*(1. - f)*(wt_press_s(2:nzp-1,itype) - wt_press(2:nzp-1,itype))&
                                  *(0.5*(t1_s(2:nzp-1)+t1_s(3:nzp)) - 0.5*(t1(2:nzp-1)+t1(3:nzp)))
       if (level >= 1) then
          qwt_press(2:nzp-1,itype) = f*qwt_press_s(2:nzp-1,itype) + (1. - f)*qwt_press(2:nzp-1,itype) &
                                   + f*(1. - f)*(wt_press_s(2:nzp-1,itype) - wt_press(2:nzp-1,itype))&
                                     *(0.5*(q1_s(2:nzp-1)+q1_s(3:nzp)) - 0.5*(q1(2:nzp-1)+q1(3:nzp)))
          tqt_press(:,itype) = f*tqt_press_s(:,itype) + (1. - f)*tqt_press(:,itype) &
                             + f*(1. - f)*(t1_s(:) - t1(:))*(qt_press_s(:,itype) - qt_press(:,itype)) &
                             + f*(1. - f)*(tt_press_s(:,itype) - tt_press(:,itype))*(q1_s(:) - q1(:))
       end if
    end do
 
    ! Update first-order mean and tendency statistics
    x1(:,:)         = f*x1_s(:,:)         + (1. - f)*x1(:,:)
    xt(:,:,:)       = f*xt_s(:,:,:)       + (1. - f)*xt(:,:,:)
    xt_press(:,:,:) = f*xt_press_s(:,:,:) + (1. - f)*xt_press(:,:,:)
    x1_misc(:,:)    = f*x1_misc_s(:,:)    + (1. - f)*x1_misc(:,:)
    x_surf(:)       = f*x_surf_s(:)       + (1. - f)*x_surf(:)

    ! Update sample counts
    n_local  = n_local  + n_s_local
    n_global = n_global + n_s_global

  end subroutine update_stat_slab

  !
  ! write_stat_slab
  !
  subroutine write_stat_slab(time, fsttm, nsmp)
    use modstat_io, only : write_stat_var, advance_stat_time, ntypes
    use grid, only       : level, umean, vmean, th00, nzp, isfctyp
    use mpi_interface, only : myid
    implicit none

    real, intent(in)     :: time, fsttm
    integer, intent(in)  :: nsmp

    integer, allocatable :: n_array(:)
    integer :: istat, itype

    !
    allocate(n_array(nzp))
    n_array(:) = n_local

    ! Add mean state
    u1(:) = u1(:) + umean
    v1(:) = v1(:) + vmean
    t1(:) = t1(:) + th00

    ! Separate tendencies
    do istat = 1,nstats_x1
       do itype = ntypes,2,-1
          xt(:,itype,istat)       = xt(:,itype,istat)       - xt(:,itype-1,istat)
          xt_press(:,itype,istat) = xt_press(:,itype,istat) - xt_press(:,itype-1,istat)
       end do
    end do
    do istat = 1,nstats_x2
       do itype = ntypes,2,-1
          x2t(:,itype,istat)       = x2t(:,itype,istat)       - x2t(:,itype-1,istat)
          x2t_press(:,itype,istat) = x2t_press(:,itype,istat) - x2t_press(:,itype-1,istat)
       end do
    end do

    ! Write time
    call advance_stat_time(time, fsttm, nsmp, ncid, irec)

    ! Write statistics
    call write_stat_var(ncid, 'n', n_array, irec)

    call write_stat_var(ncid, 'u', u1, irec)
    call write_stat_var(ncid, 'v', v1, irec)
    call write_stat_var(ncid, 'w', w1, irec)
    call write_stat_var(ncid, 't', t1, irec)
    if (level >= 1) call write_stat_var(ncid, 'q', q1, irec)

    call write_stat_var(ncid, 'ut', ut, irec)
    call write_stat_var(ncid, 'vt', vt, irec)
    call write_stat_var(ncid, 'wt', wt, irec)
    call write_stat_var(ncid, 'tt', tt, irec)
    if (level >= 1) call write_stat_var(ncid, 'qt', qt, irec)

    call write_stat_var(ncid, 'ut_press', ut_press, irec)
    call write_stat_var(ncid, 'vt_press', vt_press, irec)
    call write_stat_var(ncid, 'wt_press', wt_press, irec)

    call write_stat_var(ncid, 'uw_surf', uw_surf1, irec)
    call write_stat_var(ncid, 'vw_surf', vw_surf1, irec)
    call write_stat_var(ncid, 'tw_surf', tw_surf1, irec)
    if (level >= 1) call write_stat_var(ncid, 'qw_surf', qw_surf1, irec)
    if (isfctyp == 2 .or. isfctyp == 4) then
       call write_stat_var(ncid, 't_surf', t_surf1, irec)
       if (level >= 1) call write_stat_var(ncid, 'q_surf', q_surf1, irec)
    end if

    write(*,*) time, myid, t_surf1

    call write_stat_var(ncid, 'u2', u2, irec)
    call write_stat_var(ncid, 'v2', v2, irec)
    call write_stat_var(ncid, 'w2', w2, irec)
    call write_stat_var(ncid, 't2', t2, irec)
    if (level >= 1) call write_stat_var(ncid, 'q2', q2, irec)
    call write_stat_var(ncid, 'uw', uw, irec)
    call write_stat_var(ncid, 'vw', vw, irec)
    call write_stat_var(ncid, 'tw', tw, irec)
    if (level >= 1) then
       call write_stat_var(ncid, 'qw', qw, irec)
       call write_stat_var(ncid, 'tq', tq, irec)
    end if

    call write_stat_var(ncid, 'u2t', u2t, irec)
    call write_stat_var(ncid, 'v2t', v2t, irec)
    call write_stat_var(ncid, 'w2t', w2t, irec)
    call write_stat_var(ncid, 't2t', t2t, irec)
    if (level >= 1) call write_stat_var(ncid, 'q2t', q2t, irec)
    call write_stat_var(ncid, 'uwt', uwt, irec)
    call write_stat_var(ncid, 'vwt', vwt, irec)
    call write_stat_var(ncid, 'twt', twt, irec)
    if (level >= 1) then
       call write_stat_var(ncid, 'qwt', qwt, irec)
       call write_stat_var(ncid, 'tqt', tqt, irec)
    end if

    call write_stat_var(ncid, 'u2t_press', u2t_press, irec)
    call write_stat_var(ncid, 'v2t_press', v2t_press, irec)
    call write_stat_var(ncid, 'w2t_press', w2t_press, irec)
    call write_stat_var(ncid, 'uwt_press', uwt_press, irec)    
    call write_stat_var(ncid, 'vwt_press', vwt_press, irec)
    call write_stat_var(ncid, 'twt_press', twt_press, irec)
    if (level >= 1) then
       call write_stat_var(ncid, 'qwt_press', qwt_press, irec)
       call write_stat_var(ncid, 'tqt_press', tqt_press, irec)
    end if

    call write_stat_var(ncid, 'b',         b1,        irec)
    call write_stat_var(ncid, 'N2',        N2,        irec)
    call write_stat_var(ncid, 'pi1',       pi11,      irec)
    call write_stat_var(ncid, 'wfls',      wfls1,     irec)
    call write_stat_var(ncid, 'dthldtls',  dthldtls1, irec)
    call write_stat_var(ncid, 'e',         e,         irec)
    call write_stat_var(ncid, 'et_shear',  et_shear,  irec)
    call write_stat_var(ncid, 'et_diss',   et_diss,   irec)
    call write_stat_var(ncid, 'km',        km,        irec)
    call write_stat_var(ncid, 'kh',        kh,        irec)
    if (level >= 2) then
       call write_stat_var(ncid, 'dqtdtls', dqtdtls1, irec)
       call write_stat_var(ncid, 'th',      th1,      irec)
       call write_stat_var(ncid, 'l',       l1,       irec)
    end if

    call write_stat_var(ncid, 'b2',  b2,  irec)
    call write_stat_var(ncid, 'bw',  bw,  irec)
    call write_stat_var(ncid, 'w3',  w3,  irec)
    call write_stat_var(ncid, 't2w', t2w, irec)
    if (level >= 1) call write_stat_var(ncid, 'q2w', q2w, irec)
    call write_stat_var(ncid, 'b2w', b2w, irec)
    call write_stat_var(ncid, 'tw2', tw2, irec)
    if (level >= 1) call write_stat_var(ncid, 'qw2', qw2, irec)
    call write_stat_var(ncid, 'bw2', bw2, irec)

    ! Reset statistics computation
    n_local  = 0.
    n_global = 0.

    x1(:,:)          = 0.
    xt(:,:,:)        = 0.
    xt_press(:,:,:)  = 0.
    x2(:,:)          = 0.
    x2t(:,:,:)       = 0.
    x2t_press(:,:,:) = 0.
    x1_misc(:,:)     = 0.
    x2_misc(:,:)     = 0.
    x_surf(:)        = 0.

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

