module timers

! Per-time-step wall-clock timing, written to timing_<numrun> when the input
! flag timing_on is set (default off). One line per step gives the step's
! time split into exclusive categories:
!
!   linear     r.h.s. build and block-tridiagonal solves (solve.f90)
!   nonlinear  explicit nonlinear terms (solve.f90, nonlin /= 0 only)
!   comm       MPI layout transposes (trnsfr) and radial halo exchanges
!              (delstar, dlsq, dlstar)
!   gather     gathers to rank 0 (trnsfr0, trnsfr0e)
!   diag       diagnostics and output in the time loop (energy, hifreq, wrdump)
!   other      the rest of the step (total minus the categories above)
!
! Timers nest: timer_start pauses whichever category is running and
! timer_stop resumes it, so time is charged to the innermost category only
! (e.g. trnsfr called from inside a nonlinear block counts as comm) and the
! categories add up to the step total. Each column is the max over ranks,
! taken column by column. Timer calls outside the time loop (setup) are
! accumulated and then discarded; setup is reported as two wall-time spans
! in the file header instead.

  use mpi
  use param

  implicit none

  private
  public :: T_LINEAR, T_NONLIN, T_COMM, T_GATHER, T_DIAG
  public :: timers_init, timer_start, timer_stop, timers_mark_linstart, &
            timers_open, timers_step_begin, timers_step_end, timers_finalize

  integer, parameter :: T_LINEAR = 1, T_NONLIN = 2, T_COMM = 3, T_GATHER = 4, T_DIAG = 5
  integer, parameter :: NCAT = 5, MAXDEPTH = 16
  ! reduced/written per step: total, the NCAT categories, other
  integer, parameter :: NCOL = NCAT + 2

  logical :: enabled = .false.
  integer :: iout = -1, myrank = 0, nranks = 1, depth = 0, nsteps = 0
  integer, dimension(MAXDEPTH) :: stack
  real(IDP) :: t_start, t_last, t_step0, t_linstart0, t_linstart1, t_loop_end
  real(IDP), dimension(NCAT) :: acc
  real(IDP), dimension(NCOL) :: sums

contains

  ! Call right after MPI_INIT, before the input is read, so the setup time
  ! covers everything; whether timing is on is only known later (timers_open).
  subroutine timers_init

    integer :: ierr

    t_start = MPI_Wtime()
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)

  end subroutine timers_init

  subroutine timer_start(cat)

    integer, intent(in) :: cat
    real(IDP) :: now

    if (.not. enabled) return
    now = MPI_Wtime()
    if (depth > 0) acc(stack(depth)) = acc(stack(depth)) + (now - t_last)
    if (depth < MAXDEPTH) then
       depth = depth + 1
       stack(depth) = cat
    end if
    t_last = now

  end subroutine timer_start

  subroutine timer_stop(cat)

    integer, intent(in) :: cat
    real(IDP) :: now

    if (.not. enabled .or. depth == 0) return
    now = MPI_Wtime()
    acc(stack(depth)) = acc(stack(depth)) + (now - t_last)
    if (stack(depth) /= cat .and. myrank == 0) &
         write(0,'(" WARNING timers: stopped category ",i0," while ",i0," was running")') cat, stack(depth)
    depth = depth - 1
    t_last = now

  end subroutine timer_stop

  ! Bracket `call linstart` with mark_linstart(.true.) / mark_linstart(.false.).
  subroutine timers_mark_linstart(begin)

    logical, intent(in) :: begin

    if (begin) then
       t_linstart0 = MPI_Wtime()
    else
       t_linstart1 = MPI_Wtime()
    end if

  end subroutine timers_mark_linstart

  ! Called on every rank just before the time loop, after linstart. Turns
  ! timing on if requested and has rank 0 write the file header.
  subroutine timers_open(on, run_id, nonlin, mj, lmax, maxstp)

    logical, intent(in) :: on
    character(len=*), intent(in) :: run_id
    integer, intent(in) :: nonlin, mj, lmax, maxstp
    real(IDP), dimension(2) :: setup, setup_max
    integer :: ierr, nthreads
    character(len=3) :: build
!$  integer, external :: omp_get_max_threads

    enabled = on
    if (.not. enabled) return

    setup = [t_linstart0 - t_start, t_linstart1 - t_linstart0]
    call MPI_REDUCE(setup, setup_max, 2, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    nthreads = 1
!$  nthreads = omp_get_max_threads()
    build = 'cpu'
#ifdef _OPENACC
    build = 'gpu'
#endif

    if (myrank == 0) then
       open(newunit=iout, file='timing_'//trim(run_id), status='replace', form='formatted')
       write(iout,'("# FAR3d timing  ranks=",i0," omp_threads=",i0," build=",a," nonlin=",i0, &
            &" mj=",i0," lmax=",i0," maxstp=",i0)') nranks, nthreads, trim(build), nonlin, mj, lmax, maxstp
       write(iout,'("# setup_init_s=",es12.5,"  setup_linstart_s=",es12.5,"  (max over ranks)")') setup_max
       write(iout,'("# columns: seconds, max over ranks taken per column, so at ranks>1 total")')
       write(iout,'("# need not equal the sum of the other columns")')
       write(iout,'("#",a7,a12,7a12)') 'step', 'time', 'total', 'linear', 'nonlinear', 'comm', &
            'gather', 'diag', 'other'
       flush(iout)
    end if

    sums = 0.0_IDP
    nsteps = 0

  end subroutine timers_open

  subroutine timers_step_begin

    if (.not. enabled) return
    acc = 0.0_IDP
    depth = 0
    t_step0 = MPI_Wtime()
    t_last = t_step0

  end subroutine timers_step_begin

  subroutine timers_step_end(nstep, time)

    integer, intent(in) :: nstep
    real(IDP), intent(in) :: time
    real(IDP), dimension(NCOL) :: col, col_max
    integer :: ierr

    if (.not. enabled) return
    t_loop_end = MPI_Wtime()
    col(1) = t_loop_end - t_step0
    col(2:NCAT+1) = acc
    col(NCOL) = col(1) - sum(acc)
    sums = sums + col
    nsteps = nsteps + 1

    call MPI_REDUCE(col, col_max, NCOL, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (myrank == 0) then
       write(iout,'(i8,es12.4,7es12.4)') nstep, time, col_max
       flush(iout)
    end if

  end subroutine timers_step_end

  ! Called on every rank after the post-loop work (lincheck, final gathers
  ! and dump, endrun); writes the per-category totals and closes the file.
  subroutine timers_finalize

    real(IDP), dimension(NCOL+1) :: tot, tot_max
    integer :: ierr

    if (.not. enabled) return
    tot(1:NCOL) = sums
    tot(NCOL+1) = MPI_Wtime() - t_loop_end
    call MPI_REDUCE(tot, tot_max, NCOL+1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (myrank == 0) then
       write(iout,'("# summary: steps=",i0," sum_total=",es12.5," sum_linear=",es12.5, &
            &" sum_nonlinear=",es12.5," sum_comm=",es12.5," sum_gather=",es12.5," sum_diag=",es12.5, &
            &" sum_other=",es12.5," finalize_s=",es12.5)') nsteps, tot_max
       close(iout)
    end if
    enabled = .false.

  end subroutine timers_finalize

end module timers
