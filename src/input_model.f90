module input_model

! Reader for the human-readable "Input_Model" input format: one value per
! line, each preceded by a free-form "!!!!!!!!!!! varname: description"
! comment line. This is an additive alternative to the legacy "farin"
! raw-namelist format (nam_par/nam_arr in far3d.f90) -- it populates
! exactly the same cotrol/domain/equil/dynamo module variables that the
! legacy nam_par/nam_arr namelist reads fill, in the same relative order,
! so everything downstream of input parsing (dfault, inital, resume, etc.)
! is unaffected by which format was used.
!
! Field order below mirrors far3d.f90's nam_par (far3d.f90:34-43) then
! nam_arr (far3d.f90:44-45) declarations exactly, so every farin-settable
! variable has an Input_Model equivalent. Two deliberate differences from
! the field set found in a collaborator's older Input_Model-reading code
! (FAR3d_V2_2/Source/nonlinear/inputlist.f90): widthi/gammai are read as
! full per-mode arrays (dimension(lmax)) rather than a single scalar
! broadcast to every mode, matching what nam_arr already supports here;
! and the ~35 legacy analytic-equilibrium fields (ngeneq, q0, qlamb,
! chipr, ...) are omitted entirely, since this repo has exactly one
! equilibrium path (VMEC, via equilibrium.f90's seteq/vmec) and none of
! those fields are declared anywhere in this codebase.
!
! Optional trailing fields: entries after srcsinkEP2 may be omitted, so
! Input_Model files written before they existed still read unchanged.
! Currently just timing_on (per-step timing file, see timers.f90).

  use param
  use cotrol
  use domain
  use equil
  use dynamo

  implicit none

  private
  public :: peek_input_model_header, read_input_model, read_input_model_resume

contains

  ! Reads just the header fields (nstres, numrun, numruno, eq_name) and
  ! rewinds, without touching anything else or allocating. Mirrors the
  ! legacy farin path's throwaway header reads (far3d.f90:91-92, and again
  ! at far3d.f90:141 after `dfault` -- called twice for the same two
  ! reasons: once before `dfault` to name the farprt log file (needs
  ! numrun), and once after `dfault` to learn nstres/numruno for the
  ! continuation-run dump handling, since `dfault` may have reset them to
  ! defaults in between).
  subroutine peek_input_model_header(iunit)

    integer, intent(in) :: iunit
    character(len=1) :: cdum0

    read(iunit,'(a1)') cdum0
    read(iunit,*) nstres
    read(iunit,'(a1)') cdum0
    read(iunit,*) numrun(1),numrun(2)
    read(iunit,'(a1)') cdum0
    read(iunit,*) numruno(1),numruno(2),numruno(3)
    read(iunit,'(a1)') cdum0
    read(iunit,*) eq_name

    rewind(iunit)

  end subroutine peek_input_model_header

  ! Full read for a new run. Called once, after `dfault` has set every
  ! variable's default (same call site/order as the legacy path). Performs
  ! the same allocations far3d.f90:198 and far3d.f90:205 otherwise do
  ! inline for the legacy path, then reads every nam_par/nam_arr-equivalent
  ! field in nam_par's/nam_arr's own declared order.
  subroutine read_input_model(iunit)

    integer, intent(in) :: iunit

    call read_input_model_fields(iunit)

    allocate (mm(lmax),nn(lmax),mh(lmax),nh(lmax),mmeq(leqmax),nneq(leqmax),mheq(leqmax),nheq(leqmax))
    allocate (eta(0:mj),widthi(lmax),gammai(lmax))

    call read_input_model_arrays(iunit)

  end subroutine read_input_model

  ! Continuation-run counterpart to read_input_model: re-reads the same
  ! fields into already-allocated arrays, mirroring how initialize.f90's
  ! `resume` re-reads farin's nam_par/nam_arr itself (initialize.f90:136-139)
  ! rather than receiving parsed values as arguments.
  subroutine read_input_model_resume(iunit)

    integer, intent(in) :: iunit

    call read_input_model_fields(iunit)
    call read_input_model_arrays(iunit)

  end subroutine read_input_model_resume

  ! The nam_par-equivalent scalar fields, in far3d.f90:34-43's order.
  subroutine read_input_model_fields(iunit)

    integer, intent(in) :: iunit

    call skip_and_read_i(iunit,nstres)
    call skip_and_read_c2(iunit,numrun(1),numrun(2))
    call skip_and_read_c3(iunit,numruno(1),numruno(2),numruno(3))
    call skip_and_read_c(iunit,eq_name)

    call skip_and_read_i(iunit,maxstp)
    call skip_and_read_i(iunit,ndump)
    call skip_and_read_i(iunit,nprint)
    call skip_and_read_i(iunit,ndiag)
    call skip_and_read_i(iunit,lplots)
    call skip_and_read_i(iunit,itime)
    call skip_and_read_r(iunit,dt0)
    call skip_and_read_i(iunit,nonlin)
    call skip_and_read_i(iunit,mj)
    call skip_and_read_i(iunit,lmax)
    call skip_and_read_i(iunit,leqmax)
    call skip_and_read_i(iunit,ni)
    call skip_and_read_i(iunit,nis)
    call skip_and_read_i(iunit,ne)
    call skip_and_read_r(iunit,delta)
    call skip_and_read_r(iunit,rc)
    call skip_and_read_r(iunit,fti)
    call skip_and_read_r(iunit,fte)

    call skip_and_read_i(iunit,difnr_on)
    call skip_and_read_r(iunit,stdifp)
    call skip_and_read_r(iunit,stdifu)
    call skip_and_read_r(iunit,stdifv)
    call skip_and_read_r(iunit,eps)
    call skip_and_read_r(iunit,bet0)
    call skip_and_read_r(iunit,etascl)
    call skip_and_read_r(iunit,reta)
    call skip_and_read_r(iunit,eta0)
    call skip_and_read_r(iunit,etalmb)
    call skip_and_read_i(iunit,ietaeq)

    call skip_and_read_r(iunit,s)
    call skip_and_read_r(iunit,gamma)
    call skip_and_read_i(iunit,ipert)
    call skip_and_read_r(iunit,pertscl)
    call skip_and_read_i(iunit,m0dy)
    call skip_and_read_i(iunit,nocpl)
    call skip_and_read_r(iunit,xle)
    call skip_and_read_r(iunit,omcy)
    call skip_and_read_r(iunit,bet0_f)
    call skip_and_read_r(iunit,stdifnf)
    call skip_and_read_r(iunit,stdifvf)
    call skip_and_read_r(iunit,stdifnalp)
    call skip_and_read_r(iunit,stdifvalp)

    call skip_and_read_i(iunit,ext_prof)
    call skip_and_read_c(iunit,ext_prof_name)
    call skip_and_read_i(iunit,epflr_on)
    call skip_and_read_r(iunit,r_epflr)
    call skip_and_read_i(iunit,alpha_on)
    call skip_and_read_i(iunit,iflr_on)
    call skip_and_read_r(iunit,iflr)
    call skip_and_read_r(iunit,LcA0)
    call skip_and_read_r(iunit,LcA1)
    call skip_and_read_r(iunit,LcA2)
    call skip_and_read_r(iunit,LcA3)

    call skip_and_read_r(iunit,Adens)
    call skip_and_read_r(iunit,Bdens)
    call skip_and_read_i(iunit,twofl_on)
    call skip_and_read_r(iunit,dpres)
    call skip_and_read_r(iunit,bet0_alp)
    call skip_and_read_r(iunit,Adensalp)
    call skip_and_read_r(iunit,Bdensalp)
    call skip_and_read_r(iunit,LcA0alp)
    call skip_and_read_r(iunit,LcA1alp)
    call skip_and_read_r(iunit,LcA2alp)
    call skip_and_read_r(iunit,LcA3alp)
    call skip_and_read_r(iunit,omcyalp)
    call skip_and_read_r(iunit,r_epflralp)

    call skip_and_read_r(iunit,omegar)
    call skip_and_read_i(iunit,ieldamp_on)
    call skip_and_read_r(iunit,omcyb)
    call skip_and_read_r(iunit,rbound)
    call skip_and_read_i(iunit,trapped_on)
    call skip_and_read_r(iunit,betath_factor)
    call skip_and_read_i(iunit,spe1)
    call skip_and_read_i(iunit,spe2)
    call skip_and_read_i(iunit,EP_dens_on)

    call skip_and_read_i(iunit,EP_vel_on)
    call skip_and_read_i(iunit,Alpha_dens_on)
    call skip_and_read_i(iunit,Alpha_vel_on)
    call skip_and_read_i(iunit,DIIID_u)
    call skip_and_read_i(iunit,Eq_vel_on)
    call skip_and_read_i(iunit,Eq_velp_on)
    call skip_and_read_i(iunit,q_prof_on)
    call skip_and_read_r(iunit,deltaq)
    call skip_and_read_r(iunit,deltaiota)
    call skip_and_read_i(iunit,Eq_Presseq_on)

    call skip_and_read_i(iunit,Eq_Presstot_on)
    call skip_and_read_i(iunit,Edge_on)
    call skip_and_read_i(iunit,edge_p)
    call skip_and_read_i(iunit,Auto_grid_on)
    call skip_and_read_i(iunit,nopsievol_on)
    call skip_and_read_i(iunit,noprevol_on)
    call skip_and_read_i(iunit,nonfevol_on)
    call skip_and_read_i(iunit,nonalpevol_on)

    call skip_and_read_i(iunit,src_sink_th_on)
    call skip_and_read_i(iunit,src_sink_EP1_on)
    call skip_and_read_i(iunit,src_sink_EP2_on)
    call skip_and_read_i(iunit,src_sink_DIIID_on)
    call skip_and_read_i(iunit,src_sink_ITER_on)
    call skip_and_read_r(iunit,rsrc)
    call skip_and_read_r(iunit,wsrc)
    call skip_and_read_r(iunit,asrc)

    call skip_and_read_r(iunit,rsrc_EP1)
    call skip_and_read_r(iunit,wsrc_EP1)
    call skip_and_read_r(iunit,asrc_EP1)
    call skip_and_read_r(iunit,rsrc_EP2)
    call skip_and_read_r(iunit,wsrc_EP2)
    call skip_and_read_r(iunit,asrc_EP2)

    ! Meaning undocumented -- preserved for compatibility, see
    ! src/dfault.f90:151-155 for defaults. See a prior investigation
    ! (session history) before assuming these are safe to ignore.
    call skip_and_read_r(iunit,AWfctr)
    call skip_and_read_r(iunit,Nfctr)
    call skip_and_read_r(iunit,AWfctr_dif)
    call skip_and_read_r(iunit,Rfctr)
    call skip_and_read_r(iunit,Wfctr)

    call skip_and_read_i(iunit,B_par_on)
    call skip_and_read_l2(iunit,old_rd)

    ! matrix_out: declared (globals.f90:44) but not yet consumed anywhere
    ! in src/ (globals.f90:110 -- "not available yet"). Read here purely
    ! so the new format doesn't silently drop a field this repo already
    ! declares; harmless no-op today.
    call skip_and_read_l2(iunit,matrix_out)

  end subroutine read_input_model_fields

  ! The nam_arr-equivalent array fields, in far3d.f90:44-45's order.
  subroutine read_input_model_arrays(iunit)

    integer, intent(in) :: iunit
    character(len=1) :: cdum0

    read(iunit,'(a1)') cdum0
    read(iunit,*) mm(1:lmax)
    read(iunit,'(a1)') cdum0
    read(iunit,*) nn(1:lmax)
    read(iunit,'(a1)') cdum0
    read(iunit,*) mmeq(1:leqmax)
    read(iunit,'(a1)') cdum0
    read(iunit,*) nneq(1:leqmax)
    read(iunit,'(a1)') cdum0
    read(iunit,*) widthi(1:lmax)
    read(iunit,'(a1)') cdum0
    read(iunit,*) gammai(1:lmax)
    read(iunit,'(a1)') cdum0
    read(iunit,*) cnep(0:10)
    read(iunit,'(a1)') cdum0
    read(iunit,*) ctep(0:10)
    read(iunit,'(a1)') cdum0
    read(iunit,*) cvep(0:10)
    read(iunit,'(a1)') cdum0
    read(iunit,*) cnfp(0:10)
    read(iunit,'(a1)') cdum0
    read(iunit,*) cvfp(0:10)
    read(iunit,'(a1)') cdum0
    read(iunit,*) cnfpalp(0:10)
    read(iunit,'(a1)') cdum0
    read(iunit,*) cvfpalp(0:10)
    read(iunit,'(a1)') cdum0
    read(iunit,*) eqvt(0:10)
    read(iunit,'(a1)') cdum0
    read(iunit,*) eqvp(0:10)
    read(iunit,'(a1)') cdum0
    read(iunit,*) srcsinkth(0:10)
    read(iunit,'(a1)') cdum0
    read(iunit,*) srcsinkEP1(0:10)
    read(iunit,'(a1)') cdum0
    read(iunit,*) srcsinkEP2(0:10)

    ! Optional trailing fields (keep their dfault value when absent).
    call skip_and_read_i_opt(iunit,timing_on)

  end subroutine read_input_model_arrays

  subroutine skip_and_read_i(iunit,var)
    integer, intent(in) :: iunit
    integer, intent(out) :: var
    character(len=1) :: cdum0
    read(iunit,'(a1)') cdum0
    read(iunit,*) var
  end subroutine skip_and_read_i

  ! Like skip_and_read_i, but leaves var unchanged at end of file.
  subroutine skip_and_read_i_opt(iunit,var)
    integer, intent(in) :: iunit
    integer, intent(inout) :: var
    character(len=1) :: cdum0
    integer :: ios, tmp
    read(iunit,'(a1)',iostat=ios) cdum0
    if (ios /= 0) return
    read(iunit,*,iostat=ios) tmp
    if (ios == 0) var = tmp
  end subroutine skip_and_read_i_opt

  subroutine skip_and_read_r(iunit,var)
    integer, intent(in) :: iunit
    real(IDP), intent(out) :: var
    character(len=1) :: cdum0
    read(iunit,'(a1)') cdum0
    read(iunit,*) var
  end subroutine skip_and_read_r

  subroutine skip_and_read_l2(iunit,var)
    integer, intent(in) :: iunit
    logical, intent(out) :: var
    character(len=1) :: cdum0
    read(iunit,'(a1)') cdum0
    read(iunit,*) var
  end subroutine skip_and_read_l2

  subroutine skip_and_read_c(iunit,var)
    integer, intent(in) :: iunit
    character(len=*), intent(out) :: var
    character(len=1) :: cdum0
    read(iunit,'(a1)') cdum0
    read(iunit,*) var
  end subroutine skip_and_read_c

  subroutine skip_and_read_c2(iunit,var1,var2)
    integer, intent(in) :: iunit
    character(len=*), intent(out) :: var1,var2
    character(len=1) :: cdum0
    read(iunit,'(a1)') cdum0
    read(iunit,*) var1,var2
  end subroutine skip_and_read_c2

  subroutine skip_and_read_c3(iunit,var1,var2,var3)
    integer, intent(in) :: iunit
    character(len=*), intent(out) :: var1,var2,var3
    character(len=1) :: cdum0
    read(iunit,'(a1)') cdum0
    read(iunit,*) var1,var2,var3
  end subroutine skip_and_read_c3

end module input_model
