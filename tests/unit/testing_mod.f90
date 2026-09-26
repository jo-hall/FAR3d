module testing_mod

! Minimal assertion helpers for the FAR3d unit tests. Each test program calls
! the check_* routines, which print one PASS/FAIL line per check, then calls
! finish_tests, which stops with a nonzero exit code if anything failed so
! CTest reports the test as failed.

  use param, only: IDP

  implicit none

  private
  public :: check_true, check_close, check_close_array, check_known_bug, finish_tests

  integer :: n_checks = 0
  integer :: n_failed = 0
  integer :: n_known = 0

contains

  subroutine record(ok, name, detail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: name, detail
    n_checks = n_checks + 1
    if (ok) then
       write(*,'("PASS  ",a)') name
    else
       n_failed = n_failed + 1
       write(*,'("FAIL  ",a,"  ",a)') name, detail
    end if
  end subroutine record

  subroutine check_true(cond, name)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: name
    call record(cond, name, "(condition false)")
  end subroutine check_true

  ! |actual - expected| <= atol + rtol*|expected|
  subroutine check_close(actual, expected, name, rtol, atol)
    real(IDP), intent(in) :: actual, expected
    character(len=*), intent(in) :: name
    real(IDP), intent(in), optional :: rtol, atol
    real(IDP) :: rt, at
    character(len=120) :: detail
    rt = 1.0e-12_IDP
    at = 0.0_IDP
    if (present(rtol)) rt = rtol
    if (present(atol)) at = atol
    write(detail,'("actual=",es23.15," expected=",es23.15)') actual, expected
    call record(abs(actual-expected) <= at + rt*abs(expected), name, trim(detail))
  end subroutine check_close

  ! Max-norm comparison: max|actual - expected| <= atol + rtol*max|expected|
  subroutine check_close_array(actual, expected, name, rtol, atol)
    real(IDP), dimension(:), intent(in) :: actual, expected
    character(len=*), intent(in) :: name
    real(IDP), intent(in), optional :: rtol, atol
    real(IDP) :: rt, at, err, scale
    character(len=120) :: detail
    rt = 1.0e-12_IDP
    at = 0.0_IDP
    if (present(rtol)) rt = rtol
    if (present(atol)) at = atol
    if (size(actual) /= size(expected)) then
       call record(.false., name, "(size mismatch)")
       return
    end if
    err = maxval(abs(actual-expected))
    scale = maxval(abs(expected))
    write(detail,'("max abs err=",es12.4," max|expected|=",es12.4)') err, scale
    call record(err <= at + rt*scale, name, trim(detail))
  end subroutine check_close_array

  ! A check that documents a known, not-yet-fixed defect in the code under
  ! test. cond is the CORRECT behavior. While the bug exists this prints XFAIL
  ! and does not fail the suite; once fixed it prints XPASS as a reminder to
  ! turn it into a regular check.
  subroutine check_known_bug(cond, name, bug)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: name, bug
    n_known = n_known + 1
    if (cond) then
       write(*,'("XPASS ",a,"  (known bug appears fixed: ",a,"; make this a regular check)")') name, bug
    else
       write(*,'("XFAIL ",a,"  (known bug: ",a,")")') name, bug
    end if
  end subroutine check_known_bug

  subroutine finish_tests(suite)
    character(len=*), intent(in) :: suite
    write(*,'(/,a,": ",i0," checks, ",i0," failed, ",i0," known-bug checks")') suite, n_checks, n_failed, n_known
    if (n_failed > 0) error stop 1
  end subroutine finish_tests

end module testing_mod
