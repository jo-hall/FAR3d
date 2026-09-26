program test_spectral_ops

! Unit tests for the Fourier-space operators behind the nonlinear terms
! (src/mult_mod.f90, src/dbyd.f90).
!
! FAR3d stores a real field f(theta,zeta) as coefficients on a list of modes
! (mm(l), nn(l)) that contains both (m,n) and its mirror (-m,-n). With
! phi = m*theta + n*zeta for the "positive" member of each pair (n > 0, or
! n = 0 and m > 0), a field of type +1 is
!     f = sum a*cos(phi) + b*sin(phi),  a = f(ll(m,n)), b = f(ll(-m,-n))
! and a field of type -1 swaps the roles:
!     f = sum a*sin(phi) + b*cos(phi).
! The (0,0) slot is the constant term for either type.
!
! Each test builds random coefficients, evaluates the fields on a
! (theta,zeta) grid, computes the operation pointwise, projects back onto the
! mode list, and compares with FAR3d's spectral result. The grid is fine
! enough that the projection is exact (no aliasing), so the comparison is
! exact to roundoff.

  use param, only: IDP
  use processor, only: myPE
  use var_para, only: mj_start, mj_end
  use domain
  use equil, only: qqinv
  use mult_mod, only: mult
  use dbyd, only: dbydth_par, grdpar
  use testing_mod

  implicit none

  integer, parameter :: MMODE = 4      ! |m| <= MMODE
  integer, parameter :: NMODE = 2      ! 0 <= |n| <= NMODE
  integer, parameter :: NTH = 48, NZT = 24
  real(IDP), parameter :: TWOPI = 8.0_IDP*atan(1.0_IDP)

  real(IDP), allocatable :: g(:,:), h(:,:), f(:,:), fref(:,:), fin(:,:)
  integer :: gt, ht

  myPE = 0
  call setup_modes()

  allocate(g(mj_start:mj_end,0:lmax), h(mj_start:mj_end,0:lmax), f(mj_start:mj_end,0:lmax), &
           fref(mj_start:mj_end,0:lmax), fin(mj_start:mj_end,0:lmax))

  ! ---- products for every combination of input parities ----
  do gt = -1, 1, 2
     do ht = -1, 1, 2
        call random_field(g); call random_field(h)
        f = 0.0_IDP
        call mult(f, g, gt, h, ht, 0.0_IDP, 1.0_IDP)
        call reference_product(g, gt, h, ht, fref)
        call check_close_array(pack(f(:,1:lmax), .true.), pack(fref(:,1:lmax), .true.), &
             "mult: types ("//tstr(gt)//","//tstr(ht)//") -> "//tstr(gt*ht)//" matches grid product", &
             rtol=1.0e-12_IDP)
     end do
  end do

  ! ---- accumulate semantics: f <- c1*f + c2*(g*h) ----
  call random_field(g); call random_field(h); call random_field(fin)
  f = fin
  call mult(f, g, 1, h, -1, 0.5_IDP, 3.0_IDP)
  call reference_product(g, 1, h, -1, fref)
  fref = 0.5_IDP*fin + 3.0_IDP*fref
  call check_close_array(pack(f(:,1:lmax), .true.), pack(fref(:,1:lmax), .true.), &
       "mult: f <- c1*f + c2*g*h", rtol=1.0e-12_IDP)

  ! ---- the product with the constant field 1 is the identity ----
  g = 0.0_IDP
  g(:,ll(0,0)) = 1.0_IDP
  call random_field(h)
  f = 0.0_IDP
  call mult(f, g, 1, h, -1, 0.0_IDP, 1.0_IDP)
  call check_close_array(pack(f(:,1:lmax), .true.), pack(h(:,1:lmax), .true.), &
       "mult: 1 * h = h", rtol=1.0e-14_IDP)

  ! ---- angular derivatives ----
  call check_dtheta(1)
  call check_dtheta(-1)
  call check_grdpar(1)
  call check_grdpar(-1)

  call finish_tests("test_spectral_ops")

contains

  function tstr(t) result(s)
    integer, intent(in) :: t
    character(len=2) :: s
    if (t > 0) then
       s = "+1"
    else
       s = "-1"
    end if
  end function tstr

  ! Mode list: n=0 with m = -MMODE..MMODE (symmetric, as FAR3d requires),
  ! and for n = 1..NMODE both (m,n) and (-m,-n) for m = -MMODE..MMODE.
  subroutine setup_modes()
    integer :: m, n, l, j

    lmax = (2*MMODE+1)*(1 + 2*NMODE)
    allocate(mm(lmax), nn(lmax), signl(lmax))
    l = 0
    do m = -MMODE, MMODE
       l = l+1; mm(l) = m; nn(l) = 0
    end do
    do n = 1, NMODE
       do m = -MMODE, MMODE
          l = l+1; mm(l) = m;  nn(l) = n
          l = l+1; mm(l) = -m; nn(l) = -n
       end do
    end do
    do l = 1, lmax
       signl(l) = 1
       if (nn(l) < 0 .or. (nn(l) == 0 .and. mm(l) < 0)) signl(l) = -1
       if (nn(l) == 0 .and. mm(l) == 0) signl(l) = 0
    end do

    mmin = minval(mm); mmax = maxval(mm)
    nmin = minval(nn); nmax = maxval(nn)
    mmaxx = max(mmax, abs(mmin))
    allocate(ll(mmin:mmax, nmin:nmax))
    ll = 0
    do l = 1, lmax
       ll(mm(l), nn(l)) = l
    end do

    call band_limits()

    ! band width of the complex representation, as setmod computes it
    mxmband = 0
    do n = 0, nmax
       mxmband = max(mxmband, mmend(n) - mmstart(n) + 1)
    end do

    ! a small radial grid; one MPI process owns all of it
    mj = 3
    mj_start = 0
    mj_end = mj
    allocate(r(0:mj), rinv(0:mj), qqinv(0:mj))
    r = (/ 0.0_IDP, 0.3_IDP, 0.6_IDP, 1.0_IDP /)
    rinv(0) = 0.0_IDP
    do j = 1, mj
       rinv(j) = 1.0_IDP/r(j)
    end do
    qqinv = (/ 1.1_IDP, 0.9_IDP, 0.6_IDP, 0.35_IDP /)
  end subroutine setup_modes

  ! m-range of each n row in the complex representation. This mirrors the
  ! perturbation-mode half of mmlims (src/initialize.f90), which is private
  ! to module initialize and so can't be called from here.
  subroutine band_limits()
    integer :: n, l
    allocate(mmstart(-nmax:nmax), mmend(-nmax:nmax))
    mmstart(0:nmax) = 1000
    mmend(0:nmax) = -1000
    do l = 1, lmax
       n = nn(l)
       if (n >= 0) then
          mmstart(n) = min(mmstart(n), mm(l))
          mmend(n) = max(mmend(n), mm(l))
       else
          mmstart(-n) = min(mmstart(-n), -mm(l))
          mmend(-n) = max(mmend(-n), -mm(l))
       end if
    end do
    do n = 0, nmax
       mmstart(-n) = -mmend(n)
       mmend(-n) = -mmstart(n)
    end do
  end subroutine band_limits

  subroutine random_field(z)
    real(IDP), intent(out) :: z(mj_start:,0:)
    call random_number(z)
    z = z - 0.5_IDP
    z(:,0) = 0.0_IDP
  end subroutine random_field

  ! basis function of slot l for a field of the given type, and its
  ! theta and zeta derivatives
  subroutine basis(l, ftype, th, zt, b, db_th, db_zt)
    integer, intent(in) :: l, ftype
    real(IDP), intent(in) :: th, zt
    real(IDP), intent(out) :: b, db_th, db_zt
    real(IDP) :: ph
    integer :: m, n
    logical :: use_cos

    ! work with the positive member (m,n) of the pair
    m = mm(l)*signl(l)
    n = nn(l)*signl(l)
    if (signl(l) == 0) then
       b = 1.0_IDP; db_th = 0.0_IDP; db_zt = 0.0_IDP
       return
    end if
    ph = m*th + n*zt
    ! type +1: slot of the positive member holds cos, mirror slot holds sin
    ! type -1: the reverse
    use_cos = (signl(l) > 0 .eqv. ftype > 0)
    if (use_cos) then
       b = cos(ph); db_th = -m*sin(ph); db_zt = -n*sin(ph)
    else
       b = sin(ph); db_th = m*cos(ph);  db_zt = n*cos(ph)
    end if
  end subroutine basis

  ! Project a function sampled on the grid onto the mode list for a field of
  ! type ftype: coefficient = <F*b>/<b*b>.
  subroutine project(fgrid, ftype, j, out)
    real(IDP), intent(in) :: fgrid(NTH,NZT)
    integer, intent(in) :: ftype, j
    real(IDP), intent(inout) :: out(mj_start:,0:)
    integer :: l, it, iz
    real(IDP) :: th, zt, b, d1, d2, num, den
    do l = 1, lmax
       num = 0.0_IDP; den = 0.0_IDP
       do iz = 1, NZT
          zt = TWOPI*(iz-1)/NZT
          do it = 1, NTH
             th = TWOPI*(it-1)/NTH
             call basis(l, ftype, th, zt, b, d1, d2)
             num = num + fgrid(it,iz)*b
             den = den + b*b
          end do
       end do
       out(j,l) = num/den
    end do
  end subroutine project

  subroutine evaluate(z, ztype, j, fgrid, dth, dzt)
    real(IDP), intent(in) :: z(mj_start:,0:)
    integer, intent(in) :: ztype, j
    real(IDP), intent(out) :: fgrid(NTH,NZT), dth(NTH,NZT), dzt(NTH,NZT)
    integer :: l, it, iz
    real(IDP) :: th, zt, b, d1, d2
    fgrid = 0.0_IDP; dth = 0.0_IDP; dzt = 0.0_IDP
    do iz = 1, NZT
       zt = TWOPI*(iz-1)/NZT
       do it = 1, NTH
          th = TWOPI*(it-1)/NTH
          do l = 1, lmax
             call basis(l, ztype, th, zt, b, d1, d2)
             fgrid(it,iz) = fgrid(it,iz) + z(j,l)*b
             dth(it,iz) = dth(it,iz) + z(j,l)*d1
             dzt(it,iz) = dzt(it,iz) + z(j,l)*d2
          end do
       end do
    end do
  end subroutine evaluate

  subroutine reference_product(g, gtype, h, htype, out)
    real(IDP), intent(in) :: g(mj_start:,0:), h(mj_start:,0:)
    integer, intent(in) :: gtype, htype
    real(IDP), intent(out) :: out(mj_start:,0:)
    real(IDP) :: gg(NTH,NZT), hh(NTH,NZT), w1(NTH,NZT), w2(NTH,NZT)
    integer :: j
    out = 0.0_IDP
    do j = mj_start, mj_end
       call evaluate(g, gtype, j, gg, w1, w2)
       call evaluate(h, htype, j, hh, w1, w2)
       call project(gg*hh, gtype*htype, j, out)
    end do
  end subroutine reference_product

  ! dbydth_par computes (1/r) d/dtheta; the result has the opposite type.
  subroutine check_dtheta(ftype)
    integer, intent(in) :: ftype
    real(IDP), allocatable :: a(:,:), d(:,:), dref(:,:)
    real(IDP) :: fg(NTH,NZT), dth(NTH,NZT), dzt(NTH,NZT)
    integer :: j
    allocate(a(mj_start:mj_end,0:lmax), d(mj_start:mj_end,0:lmax), dref(mj_start:mj_end,0:lmax))
    call random_field(a)
    d = 0.0_IDP; dref = 0.0_IDP
    call dbydth_par(d, a, ftype, 0.0_IDP, 1.0_IDP, 0)
    do j = 1, mj_end
       call evaluate(a, ftype, j, fg, dth, dzt)
       call project(dth*rinv(j), -ftype, j, dref)
    end do
    call check_close_array(pack(d(1:,1:lmax), .true.), pack(dref(1:,1:lmax), .true.), &
         "dbydth_par: (1/r) d/dtheta of a type "//tstr(ftype)//" field", rtol=1.0e-12_IDP)
  end subroutine check_dtheta

  ! grdpar computes (d/dzeta - iota*d/dtheta), i.e. multiplies each mode by
  ! -ltype*(n - m*iota); the result has the opposite type.
  subroutine check_grdpar(ftype)
    integer, intent(in) :: ftype
    real(IDP), allocatable :: a(:,:), d(:,:), dref(:,:)
    real(IDP) :: fg(NTH,NZT), dth(NTH,NZT), dzt(NTH,NZT)
    integer :: j
    allocate(a(mj_start:mj_end,0:lmax), d(mj_start:mj_end,0:lmax), dref(mj_start:mj_end,0:lmax))
    call random_field(a)
    d = 0.0_IDP; dref = 0.0_IDP
    call grdpar(d, a, ftype, 0.0_IDP, 1.0_IDP)
    do j = mj_start, mj_end
       call evaluate(a, ftype, j, fg, dth, dzt)
       call project(dzt - qqinv(j)*dth, -ftype, j, dref)
    end do
    call check_close_array(pack(d(:,1:lmax), .true.), pack(dref(:,1:lmax), .true.), &
         "grdpar: (d/dzeta - iota d/dtheta) of a type "//tstr(ftype)//" field", rtol=1.0e-12_IDP)
  end subroutine check_grdpar

end program test_spectral_ops
