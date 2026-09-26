program test_block_tridiag

! Unit test for the linear solver used by every time step (src/matrix.f90):
! decbt (block-tridiagonal LU) followed by solbt (forward/back substitution).
!
! For random, diagonally dominant block-tridiagonal systems -- including the
! two extra corner blocks decbt supports, (1,3) stored in c(:,:,1) and (n,n-2)
! stored in b(:,:,n) -- we pick a known solution x, build y = A*x with an
! independent dense matrix-vector product, then check that solbt recovers x.

  use param, only: IDP
  use matrix, only: decbt, solbt
  use testing_mod

  implicit none

  call run_case(m=1,  n=4,   corners=.false., name="scalar blocks, n=4")
  call run_case(m=3,  n=10,  corners=.false., name="3x3 blocks, n=10")
  call run_case(m=7,  n=50,  corners=.true.,  name="7x7 blocks, n=50, corner blocks")
  call run_case(m=56, n=40,  corners=.true.,  name="56x56 blocks (DIII-D size), n=40, corners")
  call check_rejects_small_n()

  call finish_tests("test_block_tridiag")

contains

  subroutine run_case(m, n, corners, name)
    integer, intent(in) :: m, n
    logical, intent(in) :: corners
    character(len=*), intent(in) :: name

    real(IDP), allocatable :: a(:,:,:), b(:,:,:), c(:,:,:), dense(:,:), x(:), y(:)
    integer, allocatable :: ip(:,:)
    integer :: k, ier, i, nn

    nn = m*n
    allocate(a(m,m,n), b(m,m,n), c(m,m,n), ip(m,n), dense(nn,nn), x(nn), y(nn))

    call random_matrix(a); call random_matrix(b); call random_matrix(c)
    ! diagonal dominance so the system is well conditioned (and so partial
    ! pivoting within blocks, which is all decbt does, is sufficient)
    do k = 1,n
       do i = 1,m
          a(i,i,k) = a(i,i,k) + 4.0_IDP*m
       end do
    end do
    if (.not. corners) then
       c(:,:,1) = 0.0_IDP
       b(:,:,n) = 0.0_IDP
    end if

    ! Dense reference: block row k has c(k) at column k-1, a(k) at k, b(k) at k+1,
    ! plus the corner blocks c(1) at (1,3) and b(n) at (n,n-2).
    dense = 0.0_IDP
    do k = 1,n
       call put(dense, k, k, a(:,:,k), m)
       if (k > 1) call put(dense, k, k-1, c(:,:,k), m)
       if (k < n) call put(dense, k, k+1, b(:,:,k), m)
    end do
    call put(dense, 1, 3, c(:,:,1), m)
    call put(dense, n, n-2, b(:,:,n), m)

    call random_number(x)
    x = x - 0.5_IDP
    y = matmul(dense, x)

    call decbt(m, n, a, b, c, ip, ier)
    call check_true(ier == 0, name//": decbt ier == 0")
    call solbt(m, n, a, b, c, y, ip)
    call check_close_array(y, x, name//": solbt recovers x", rtol=1.0e-11_IDP)
  end subroutine run_case

  subroutine check_rejects_small_n()
    ! decbt documents n >= 4 and returns ier = -1 otherwise
    real(IDP) :: a(2,2,3), b(2,2,3), c(2,2,3)
    integer :: ip(2,3), ier
    a = 1.0_IDP; b = 0.0_IDP; c = 0.0_IDP
    call decbt(2, 3, a, b, c, ip, ier)
    call check_true(ier == -1, "decbt rejects n < 4 with ier = -1")
  end subroutine check_rejects_small_n

  subroutine put(dense, bi, bj, blk, m)
    real(IDP), intent(inout) :: dense(:,:)
    integer, intent(in) :: bi, bj, m
    real(IDP), intent(in) :: blk(m,m)
    dense((bi-1)*m+1:bi*m, (bj-1)*m+1:bj*m) = dense((bi-1)*m+1:bi*m, (bj-1)*m+1:bj*m) + blk
  end subroutine put

  subroutine random_matrix(z)
    real(IDP), intent(out) :: z(:,:,:)
    call random_number(z)
    z = z - 0.5_IDP
  end subroutine random_matrix

end program test_block_tridiag
