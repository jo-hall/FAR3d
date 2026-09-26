program test_radial_grid

! Unit tests for the radial discretization (src/grid.f90, src/findr.f90,
! src/dbyd.f90):
!   * the uniform grid (Auto_grid_on=1) and the packed "island" grid built by
!     findr (Auto_grid_on=0) have the documented shape;
!   * the 3-point nonuniform finite-difference weights dc1m/dc1p (first
!     derivative) and dc2m/dc2p (second derivative) are exact for quadratics,
!     as a second-order scheme must be, on both grids;
!   * the first-derivative error converges at second order under refinement;
!   * the cylindrical Laplacian del2c, (1/r)(r f')' - (m/r)^2 f, is exact on
!     r^m for m = 0, 1, 2.

  use param, only: IDP
  use cotrol, only: Auto_grid_on
  use domain
  use dbyd, only: dbydr0, d2bydr20, del2c
  use testing_mod

  implicit none

  real(IDP) :: e1, e2, e3

  ! ---- uniform grid ----
  call make_grid(100, auto=.true.)
  call check_close(r(0), 0.0_IDP, "uniform: r(0) = 0", atol=0.0_IDP)
  call check_close(r(mj), 1.0_IDP, "uniform: r(mj) = 1")
  call check_close(r(37), 0.37_IDP, "uniform: r(j) = j/mj")
  call check_fd_exact_on_quadratic("uniform")
  call check_laplacian_exact("uniform")

  ! ---- packed (findr) grid, parameters in the style of the DIII-D case ----
  call make_grid(100, auto=.false., ni_=40, nis_=30, ne_=30, delta_=0.25_IDP, rc_=0.625_IDP)
  call check_findr_shape(0.25_IDP, 0.625_IDP)
  call check_fd_exact_on_quadratic("findr")
  call check_laplacian_exact("findr")

  ! ---- second-order convergence of d/dr on uniform grids ----
  e1 = deriv_error(50); e2 = deriv_error(100); e3 = deriv_error(200)
  write(*,'("      d/dr sin(3r) max error: ",3es12.4)') e1, e2, e3
  call check_true(e1/e2 > 3.5_IDP .and. e1/e2 < 4.5_IDP, "d/dr error ratio 50->100 is ~4 (2nd order)")
  call check_true(e2/e3 > 3.5_IDP .and. e2/e3 < 4.5_IDP, "d/dr error ratio 100->200 is ~4 (2nd order)")

  call finish_tests("test_radial_grid")

contains

  ! Allocate the domain arrays exactly as far3d.f90 does, then call grid.
  subroutine make_grid(mj_in, auto, ni_, nis_, ne_, delta_, rc_)
    integer, intent(in) :: mj_in
    logical, intent(in) :: auto
    integer, intent(in), optional :: ni_, nis_, ne_
    real(IDP), intent(in), optional :: delta_, rc_
    external :: grid

    if (allocated(r)) deallocate(r, rinv, dc1m, dc1p, dc2m, dc2p, del2cm, del2cp)
    mj = mj_in
    allocate(r(0:mj), rinv(0:mj), dc1m(mj), dc1p(mj), dc2m(mj), dc2p(mj), del2cm(mj), del2cp(mj))
    r = 0.0_IDP; rinv = 0.0_IDP
    dc1m = 0.0_IDP; dc1p = 0.0_IDP; dc2m = 0.0_IDP; dc2p = 0.0_IDP; del2cm = 0.0_IDP; del2cp = 0.0_IDP

    if (auto) then
       Auto_grid_on = 1
    else
       Auto_grid_on = 0
       ni = ni_; nis = nis_; ne = ne_
       delta = delta_; rc = rc_
       fti = 0.95_IDP; fte = 0.95_IDP
    end if
    call grid
  end subroutine make_grid

  subroutine check_findr_shape(width, center)
    real(IDP), intent(in) :: width, center
    real(IDP) :: h, dx
    integer :: j, n_uniform
    logical :: monotone, in_band

    call check_close(r(0), 0.0_IDP, "findr: r(0) = 0", atol=0.0_IDP)
    call check_close(r(mj), 1.0_IDP, "findr: r(mj) = 1")
    monotone = all(r(1:mj) > r(0:mj-1))
    call check_true(monotone, "findr: r strictly increasing")

    ! The fine region is nis-1 uniform intervals of spacing delta/(nis-1)
    ! spanning [rc - delta/2, rc + delta/2]. Its ends are only located to the
    ! 1e-10 tolerance findr passes to its root finder (zeroin), hence 1e-8.
    dx = width/(nis-1)
    n_uniform = 0
    in_band = .true.
    do j = 1,mj
       h = r(j) - r(j-1)
       if (abs(h-dx) <= 1.0e-10_IDP*dx) then
          n_uniform = n_uniform + 1
          if (r(j-1) < center-0.5_IDP*width-1.0e-8_IDP .or. r(j) > center+0.5_IDP*width+1.0e-8_IDP) in_band = .false.
       end if
    end do
    call check_true(n_uniform == nis-1, "findr: fine region has nis-1 intervals of spacing delta/(nis-1)")
    call check_true(in_band, "findr: fine region lies inside [rc-delta/2, rc+delta/2]")
    call check_true(maxval(r(1:mj)-r(0:mj-1)) > 1.5_IDP*dx, "findr: grid coarsens outside the fine region")
  end subroutine check_findr_shape

  subroutine check_fd_exact_on_quadratic(label)
    character(len=*), intent(in) :: label
    real(IDP), allocatable :: f(:), d(:)
    allocate(f(0:mj), d(0:mj))

    f = 3.0_IDP*r**2 - 2.0_IDP*r + 0.5_IDP
    d = 0.0_IDP
    call dbydr0(d, f, 0.0_IDP, 1.0_IDP, 0)
    call check_close_array(d(1:mj-1), 6.0_IDP*r(1:mj-1) - 2.0_IDP, &
         label//": dc1m/dc1p exact first derivative of a quadratic", rtol=1.0e-10_IDP)

    d = 0.0_IDP
    call d2bydr20(d, f, 0.0_IDP, 1.0_IDP, 0)
    call check_close_array(d(1:mj-1), spread(6.0_IDP, 1, mj-1), &
         label//": dc2m/dc2p exact second derivative of a quadratic", rtol=1.0e-8_IDP)
  end subroutine check_fd_exact_on_quadratic

  subroutine check_laplacian_exact(label)
    ! del2c works on mode-indexed arrays a(0:mj,0:lmax); the mode's m enters
    ! through the -(m/r)^2 term. For f = r^m the exact result is 0 for m >= 1
    ! and 4 for f = r^2 with m = 0.
    character(len=*), intent(in) :: label
    real(IDP), allocatable :: a(:,:), d(:,:)

    if (allocated(mm)) deallocate(mm)
    lmax = 3
    allocate(mm(lmax), a(0:mj,0:lmax), d(0:mj,0:lmax))
    mm = (/ 0, 1, 2 /)
    mjm1 = mj-1
    a(:,1) = r**2          ! m = 0: laplacian = 4
    a(:,2) = r             ! m = 1: laplacian = 0
    a(:,3) = r**2          ! m = 2: laplacian = 0
    d = 0.0_IDP
    call del2c(d, a, 0.0_IDP, 1.0_IDP, 0)
    call check_close_array(d(1:mj-1,1), spread(4.0_IDP, 1, mj-1), label//": del2c(r^2), m=0 -> 4", rtol=1.0e-8_IDP)
    call check_close_array(d(1:mj-1,2), spread(0.0_IDP, 1, mj-1), label//": del2c(r), m=1 -> 0", atol=1.0e-8_IDP)
    call check_close_array(d(1:mj-1,3), spread(0.0_IDP, 1, mj-1), label//": del2c(r^2), m=2 -> 0", atol=1.0e-7_IDP)
  end subroutine check_laplacian_exact

  real(IDP) function deriv_error(mj_in)
    integer, intent(in) :: mj_in
    real(IDP), allocatable :: f(:), d(:)
    call make_grid(mj_in, auto=.true.)
    allocate(f(0:mj), d(0:mj))
    f = sin(3.0_IDP*r)
    d = 0.0_IDP
    call dbydr0(d, f, 0.0_IDP, 1.0_IDP, 0)
    deriv_error = maxval(abs(d(1:mj-1) - 3.0_IDP*cos(3.0_IDP*r(1:mj-1))))
  end function deriv_error

end program test_radial_grid
