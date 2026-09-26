program test_numerics_tools

! Unit tests for the numerical utilities in src/tools.f90 that the
! equilibrium and damping calculations rely on:
!   * spline/seval: the cubic interpolating spline used to map equilibrium
!     profiles onto FAR3d's radial grid reproduces cubic data exactly;
!   * zzdisp: the plasma dispersion function Z(zeta) used by the
!     electron-ion Landau damping terms, checked against exact identities;
!   * erf: FAR3d's own trapezoid-rule error function (it shadows the Fortran
!     intrinsic inside module tools) against the intrinsic.

  use param, only: IDP
  use tools, only: spline, seval, zzdisp, far3d_erf => erf
  use testing_mod

  implicit none

  real(IDP), parameter :: SQRTPI = 1.7724538509055160_IDP

  call test_spline_reproduces_cubic()
  call test_spline_interpolates()
  call test_plasma_dispersion()
  call test_erf()

  call finish_tests("test_numerics_tools")

contains

  real(IDP) elemental function cubic(x)
    real(IDP), intent(in) :: x
    cubic = 0.3_IDP - 1.2_IDP*x + 2.5_IDP*x**2 - 0.7_IDP*x**3
  end function cubic

  subroutine test_spline_reproduces_cubic()
    ! The end conditions of this spline (Forsythe, Malcolm & Moler) match the
    ! third derivative of cubics through the end points, so cubic data is
    ! reproduced exactly, on a nonuniform knot set too.
    integer, parameter :: n = 17
    real(IDP) :: x(n), y(n), b(n), c(n), d(n), u, err
    integer :: i

    do i = 1, n
       x(i) = ((i-1.0_IDP)/(n-1))**1.5_IDP
    end do
    y = cubic(x)
    call spline(n, x, y, b, c, d)

    err = 0.0_IDP
    do i = 0, 200
       u = i/200.0_IDP
       err = max(err, abs(seval(n, u, x, y, b, c, d) - cubic(u)))
    end do
    call check_true(err < 1.0e-12_IDP, "spline/seval reproduce a cubic exactly on nonuniform knots")
    call check_close_array(b, -1.2_IDP + 5.0_IDP*x - 2.1_IDP*x**2, &
         "spline: b(i) is the exact first derivative at the knots", rtol=1.0e-10_IDP)
  end subroutine test_spline_reproduces_cubic

  subroutine test_spline_interpolates()
    ! For smooth non-polynomial data: exact at knots, 4th-order convergence.
    real(IDP) :: e1, e2
    call check_close(spline_error(21, .true.), 0.0_IDP, "spline: passes through the knots", atol=1.0e-14_IDP)
    e1 = spline_error(21, .false.)
    e2 = spline_error(41, .false.)
    write(*,'("      spline error on sin(4x): ",2es12.4)') e1, e2
    call check_true(e1/e2 > 12.0_IDP, "spline: error drops ~16x when knot spacing halves (4th order)")
  end subroutine test_spline_interpolates

  real(IDP) function spline_error(n, at_knots)
    integer, intent(in) :: n
    logical, intent(in) :: at_knots
    real(IDP) :: x(n), y(n), b(n), c(n), d(n), u
    integer :: i
    do i = 1, n
       x(i) = (i-1.0_IDP)/(n-1)
    end do
    y = sin(4.0_IDP*x)
    call spline(n, x, y, b, c, d)
    spline_error = 0.0_IDP
    if (at_knots) then
       do i = 1, n
          spline_error = max(spline_error, abs(seval(n, x(i), x, y, b, c, d) - y(i)))
       end do
    else
       do i = 0, 999
          u = (i+0.5_IDP)/1000.0_IDP
          spline_error = max(spline_error, abs(seval(n, u, x, y, b, c, d) - sin(4.0_IDP*u)))
       end do
    end if
  end function spline_error

  subroutine test_plasma_dispersion()
    ! Z(zeta) = i*sqrt(pi)*w(zeta), w the Faddeeva function.
    !   Z(0)            = i*sqrt(pi)
    !   Im Z(x), x real = sqrt(pi)*exp(-x^2)
    !   Re Z(x), x real = -2*F(x), F the Dawson function
    !   Z(i*y), y > 0   = i*sqrt(pi)*exp(y^2)*erfc(y)
    !   Z(zeta)         ~ -1/zeta for |zeta| >> 1
    real(IDP) :: zr, zi, x, y
    integer :: k
    real(IDP), parameter :: xs(3) = (/ 0.5_IDP, 1.0_IDP, 2.0_IDP /)
    ! Dawson function F(x) at xs (Abramowitz & Stegun, Table 7.5)
    real(IDP), parameter :: dawson(3) = (/ 0.4244363835020223_IDP, 0.5380795069127684_IDP, 0.3013403889237920_IDP /)
    character(len=8) :: tag

    call zzdisp(0.0_IDP, 0.0_IDP, zr, zi)
    call check_close(zr, 0.0_IDP, "Z(0): real part 0", atol=1.0e-10_IDP)
    call check_close(zi, SQRTPI, "Z(0): imag part sqrt(pi)", rtol=1.0e-8_IDP)

    do k = 1, 3
       x = xs(k)
       write(tag,'(f4.1)') x
       call zzdisp(x, 0.0_IDP, zr, zi)
       call check_close(zi, SQRTPI*exp(-x*x), "Z(x) real x="//trim(tag)//": imag = sqrt(pi) exp(-x^2)", rtol=1.0e-8_IDP)
       call check_close(zr, -2.0_IDP*dawson(k), "Z(x) real x="//trim(tag)//": real = -2*Dawson(x)", rtol=1.0e-8_IDP)
       ! Z(-conj(zeta)) = -conj(Z(zeta)): the real part is odd in x. zzdisp
       ! only applies that reflection when x*y < 0, so for a purely real
       ! negative argument (y == 0 exactly) it returns Re Z with the wrong
       ! sign. In FAR3d, Re(zeta) is proportional to omegar, so this is hit
       ! only for omegar < 0 with zero collisionality.
       call zzdisp(-x, 0.0_IDP, zr, zi)
       call check_known_bug(abs(zr - 2.0_IDP*dawson(k)) < 1.0e-8_IDP, &
            "Z(-x) real x="//trim(tag)//": real part odd in x", "zzdisp sign error for real negative argument")
       call zzdisp(-x, 1.0e-12_IDP, zr, zi)
       call check_close(zr, 2.0_IDP*dawson(k), "Z(-x+0i) x="//trim(tag)//": real part odd in x (y > 0)", &
            rtol=1.0e-8_IDP)
    end do

    do k = 1, 3
       y = xs(k)
       write(tag,'(f4.1)') y
       call zzdisp(0.0_IDP, y, zr, zi)
       call check_close(zr, 0.0_IDP, "Z(iy) y="//trim(tag)//": real part 0", atol=1.0e-10_IDP)
       call check_close(zi, SQRTPI*erfc_scaled(y), "Z(iy) y="//trim(tag)//": imag = sqrt(pi) e^{y^2} erfc(y)", &
            rtol=1.0e-8_IDP)
    end do

    ! asymptotic regime: Z ~ -1/zeta - 1/(2 zeta^3), damped branch negligible for y > 0
    call zzdisp(30.0_IDP, 1.0_IDP, zr, zi)
    call check_close(zr, real(-1.0_IDP/cmplx(30.0_IDP,1.0_IDP,IDP) - 0.5_IDP/cmplx(30.0_IDP,1.0_IDP,IDP)**3), &
         "Z(30+i): matches asymptotic series (real)", rtol=1.0e-5_IDP)
    call check_close(zi, aimag(-1.0_IDP/cmplx(30.0_IDP,1.0_IDP,IDP) - 0.5_IDP/cmplx(30.0_IDP,1.0_IDP,IDP)**3), &
         "Z(30+i): matches asymptotic series (imag)", rtol=1.0e-3_IDP)
  end subroutine test_plasma_dispersion

  subroutine test_erf()
    ! FAR3d's erf integrates exp(-u^2) with a 1000-interval trapezoid rule,
    ! so it agrees with the intrinsic to about 1e-7.
    real(IDP) :: x, worst
    integer :: i
    worst = 0.0_IDP
    do i = -60, 60
       x = i/10.0_IDP
       worst = max(worst, abs(far3d_erf(x) - erf(x)))
    end do
    write(*,'("      max |erf_far3d - erf|: ",es12.4)') worst
    call check_true(worst < 1.0e-6_IDP, "tools erf agrees with the intrinsic erf to 1e-6 on [-6,6]")
    x = 0.7_IDP
    call check_close(far3d_erf(-x), -far3d_erf(x), "tools erf is odd", rtol=1.0e-15_IDP)
  end subroutine test_erf

end program test_numerics_tools
