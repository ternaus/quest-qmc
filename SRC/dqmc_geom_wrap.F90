module DQMC_GEOM_WRAP
  use DQMC_GEOM_PARAM
  use DQMC_HAMILT
  use DQMC_SYMM
  use DQMC_LATT
  use DQMC_RECLATT
  use DQMC_BONDS
  use DQMC_Cfg
  use DQMC_STRUCT
  implicit none


  type GeomWrap
     type(lattice_t)       :: Lattice
     type(recip_lattice_t) :: RecipLattice
     type(recip_lattice_t) :: GammaLattice
     type(hamiltonian_t)   :: Hamilt
     type(bonds_t)         :: Bonds
     type(symm_operations) :: SymmOp
     type(pairing)         :: Pairs
  end type

  contains

  subroutine DQMC_Geom_Fill(gwrap, gfile, cfg, SOP)

    type(GeomWrap),    intent(inout)  :: gwrap
    type(config), intent(inout)       :: cfg
    character(len=60), intent(in)     :: gfile
    integer, intent(in)               :: SOP
    logical                           :: found, connected

    inquire(file=gfile,exist=found)
    if(found)then
     !Check whether file is already connected
     inquire(file=gfile,opened=connected)
     if(connected)then
      !If connected retrieve the unit
      inquire(file=gfile,number=inpunit)
      rewind(inpunit)
     else
      !If not connected, find a free unit
      do inpunit=7,99
        inquire(unit=inpunit,opened=connected)
        if(.not.connected)exit
      enddo
      !open file
      open(file=gfile,unit=inpunit)
     endif
    else 
     call DQMC_Error("cannot open geom def file "//gfile, 0)
    endif


    !Scan file to see which fields are specified
    call analyze_input
  
    !Initialize basic info about real space cluster
    call init_lattice(gwrap%Lattice, SOP)

    !Construct full lattice 
    call construct_lattice(gwrap%Lattice, SOP)
 
    !Initialize basic info about reciprocal lattice
    call init_recip_latt(gwrap%Lattice,gwrap%RecipLattice,.true.,cfg)
    call init_recip_latt(gwrap%Lattice,gwrap%GammaLattice,.false.,cfg)

    !construct reciprocal lattice shifted by k-point
    call construct_recip_lattice(gwrap%RecipLattice)

    !construct reciprocal lattice that includes Gamma
    call construct_recip_lattice(gwrap%GammaLattice)

    !Construct Hamiltonian
    call construct_hamilt(gwrap%Hamilt,gwrap%Lattice,gwrap%RecipLattice,cfg)
    
    !Read point-symmetry (optional)
    call read_symm(gwrap%SymmOp)

    !pair-to-pair map of action of each symmetry in real space (also translations)
    call map_symm_lattice(gwrap%SymmOp,gwrap%Lattice, gwrap%Hamilt, SOP)

    !point-to-point map of action of each symmetry in reciprocal space
    call map_symm_recip_lattice(gwrap%SymmOp,gwrap%RecipLattice,.true.)
    call map_symm_recip_lattice(gwrap%SymmOp,gwrap%GammaLattice,.false.)

    !Group pairs of lattice points into classes
    call construct_lattice_classes(gwrap%SymmOp,gwrap%Lattice)

    !Group k-points into classes
    call construct_recip_lattice_classes(gwrap%SymmOp,gwrap%RecipLattice,.true.)
    call construct_recip_lattice_classes(gwrap%SymmOp,gwrap%GammaLattice,.false.)

    ! Fill the matrix of fourier weights
    call DQMC_Fill_FourierC(Gwrap%RecipLattice, Gwrap%Lattice)
    call DQMC_Fill_FourierC(Gwrap%GammaLattice, Gwrap%Lattice)

    !Group hopping part of Hamiltonian in classes
    call count_hop_class(gwrap%Lattice,gwrap%Hamilt)

    !Group separately local classes for U and mu
    call count_local_classes(gwrap%Lattice,gwrap%Hamilt)

    !Assign phase for Green's function
    call assign_gf_phase(gwrap%Lattice,gwrap%RecipLattice%ktwist)

    !Read Bonds (optional input)   
    call read_bonds(gwrap%Bonds, SOP)

    if(gwrap%Bonds%initialized)then
       !bond-to-bond map of action of symmetry on bonds
       call map_symm_bonds(gwrap%Bonds,gwrap%SymmOp,gwrap%Lattice)

       !Group pairs of bonds into classes
       call construct_bond_classes(gwrap%Bonds,gwrap%SymmOp)

       !Map bonds throughout the entire lattice
       call construct_pairs(gwrap%Bonds,gwrap%Pairs,gwrap%Lattice, SOP)

       !pair-to-pair map of action of symmetry on bonds
       call map_symm_pairs(gwrap%Pairs, gwrap%SymmOp)

       !Group pairs of pairs into classes
       call construct_pair_classes(gwrap%Bonds,gwrap%Pairs)
    endif

    !Assign phase to each atom (optional input)
    call assign_phase(gwrap%Lattice)

    !Write some info
    !call write_files(gwrap%Lattice,gwrap%RecipLattice,gwrap%Hamilt)

  end subroutine DQMC_Geom_Fill

 !---------------------------------------------------------------------!

  subroutine DQMC_Geom_Init(gwrap, S, cfg)    

    type(GeomWrap), intent(in)    :: gwrap   
    type(Struct),   intent(inout) :: S       ! Struct
    type(config),   intent(inout) :: cfg

    ! ... local scalar ...
    integer  :: n                ! Order of matrix T and D 
    integer  :: i, j             ! Loop iterator
    integer  :: ic, ib
    real(wp), pointer  :: clab(:,:) => null()
    real(wp), pointer  :: tupvalue(:) => null()
    real(wp), pointer  :: tdnvalue(:) => null()
    integer, pointer   :: tmp(:,:) => null()

    ! ... Executable ...
    S%checklist=.false.
  
    n        = gwrap%Lattice%nsites
    S%nSite  = n
    S%nCell  = gwrap%Lattice%ncell
    S%nGroup = gwrap%Hamilt%nlocclass
    S%Name   = 'General Geometry - Free Format'
    allocate(tmp(n,n))

    !Fill T
    call group_hopping(gwrap%Hamilt, n, S%n_t, tmp, tupvalue, tdnvalue)
    call DQMC_CCS_Compress(n,-1, tmp, S%T)
    S%checklist(STRUCT_ADJ)=.true.

    !Fill ckb
    tmp = 0
    do ib = 0, size(gwrap%hamilt%tckb,2)-1
       i = gwrap%hamilt%tckb(1,ib) + 1
       j = gwrap%hamilt%tckb(2,ib) + 1
       tmp(i,j) = gwrap%hamilt%tckb(3,ib)
       tmp(j,i) = gwrap%hamilt%tckb(3,ib)
    enddo
    call DQMC_CCS_Compress(n,-1, tmp, S%ckb)

    !Symmetry Classes
    S%nClass = gwrap%Lattice%nclass
    allocate(S%D(n,n))
    allocate(S%F(S%nClass))
    allocate(S%clabel(S%nClass))

    S%D(1:n,1:n) =  gwrap%Lattice%myclass(0: n - 1, 0: n - 1)
    S%F(:)       =  gwrap%Lattice%class_size(:)
    clab         => gwrap%Lattice%class_label
    do ic = 1, S%nClass
       write(S%clabel(ic),'(2(i5),3(f10.4), i6.3)') (int(clab(ic,j)),j=4,5),(clab(ic,j),j=1,3), S%F(ic)
    enddo

    !store GF phase
    allocate(S%gf_phase(n,n))
    do i = 1, n
       do j = 1, n
          !if(abs(int(gwrap%lattice%gf_phase(i-1,j-1)))/=1)stop 'Problem with gf_phase'
          !S%gf_phase(i,j)=int(gwrap%Lattice%gf_phase(i-1,j-1))
          S%gf_phase(i, j) = gwrap%Lattice%gf_phase(i - 1, j - 1)
       enddo
    enddo
    allocate(S%chi_phase(n, n))
    S%chi_phase = 1
    S%checklist(STRUCT_CLASS)=.true.    

    allocate(S%map(n))
    S%map(1:n)=gwrap%Hamilt%mylocclass(0:n-1)
    
    !Fill B
    if (Found_Field(BONDS_F)) then
       tmp = 0
       S%n_b = gwrap%Pairs%nbond
       do ic = 0, size(gwrap%Pairs%nbondv)-1
          do ib = 1, gwrap%Pairs%nbondv(ic)
             i = gwrap%Pairs%bond_origin(ib,ic)
             j = gwrap%Pairs%bond_end(ib,ic)
             tmp(i+1,j+1) = gwrap%Pairs%bond_number(ib,ic)
          enddo
       enddo
       call DQMC_CCS_Compress(n,-1, tmp, S%B)

       !Store symmetry in S     
       allocate(S%class_b(S%n_b,S%n_b),S%size_b(gwrap%Pairs%nclass_p))
       S%class_b=gwrap%Pairs%myclass_p
       S%size_b=gwrap%Pairs%class_size_p  
       S%nClass_b=gwrap%Pairs%nclass_p
       S%checklist(STRUCT_BOND)=.true.    

       !Waves
       S%nWave=gwrap%Pairs%nWave
       allocate(S%wlabel(S%nWave))
       S%wlabel(:)=gwrap%Pairs%wave_label(:)
       if(Found_Field(PAIRS_F))then
          allocate(S%W(S%n_b,S%nWave))
          do i=1,S%nWave
             S%W(1:S%n_b,i)=gwrap%Pairs%bond_wgt(i,1:S%n_b)
          enddo
          S%checklist(STRUCT_WAVE)=.true.
       endif
    endif

    deallocate(tmp)
    
    if(Found_field(PHASE_F))then 
       allocate(S%P(n))
       S%P(1:n)=gwrap%Lattice%phase(0:n-1)
       S%checklist(STRUCT_PHASE)=.true.
    endif

    !Setting variables

    !Set variables that are otherwise read from main input
    call CFG_Set(cfg,"n",n)
    call CFG_Set(cfg,"t_up",S%n_t,tupvalue)
    call CFG_Set(cfg,"t_dn",S%n_t,tdnvalue)
    call CFG_Set(cfg,"U",S%nGroup,gwrap%Hamilt%Uvalue)
    call CFG_Set(cfg,"mu_up",S%nGroup,gwrap%Hamilt%muupvalue)
    call CFG_Set(cfg,"mu_dn",S%nGroup,gwrap%Hamilt%mudnvalue)

    deallocate(tupvalue)
    deallocate(tdnvalue)

    S%checklist(STRUCT_INIT)=.true.

  end subroutine DQMC_Geom_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_Print_HeaderFT(Gwrap, OPT, applytwist)
    use dqmc_mpi

    type(GeomWrap), intent(in)  :: Gwrap
    integer, intent(in)         :: OPT
    logical, intent(in) :: applytwist

    ! ... Local scalar ...
    integer :: na, nk, ii, i, ik, j, jj, nkpts
    integer, pointer :: myclass_k(:) 
    real*8, pointer  :: klist(:,:) 

    ! ... Executable ...

    if (qmc_sim%rank .ne. 0) return

    na   = Gwrap%Lattice%natom
    if(applytwist)then 
       nk         = Gwrap%RecipLattice%nclass_k
       nkpts      = Gwrap%RecipLattice%nkpts
       myclass_k => Gwrap%RecipLattice%myclass_k
       klist     => Gwrap%RecipLattice%klist
    else
       nk         = Gwrap%GammaLattice%nclass_k
       nkpts      = Gwrap%GammaLattice%nkpts
       myclass_k => Gwrap%GammaLattice%myclass_k
       klist     => Gwrap%GammaLattice%klist
    endif

    !Print general info about k-space
    if(applytwist)then
      write(OPT,'(A)')' Grid for Green''s function'
    else
      write(OPT,'(A)')' Grid for spin/charge correlations'
    endif

    write(OPT,'(A)')'  K-points'
    ii=0
    write(OPT,'(A)')'  Class'
    do i=1,nk
      jj=1
      do ik=1,nkpts
        if(myclass_k(ik)==i)then
          if(jj==1)then
            write(OPT,'(2x,i3,6x,3(f10.5))')i,(klist(ik,j),j=1,Gwrap%Lattice%ndim)
            jj=-1
          else
            write(OPT,'(11x,3(f10.5))')(klist(ik,j),j=1,Gwrap%Lattice%ndim)
          endif
        endif
      enddo
      write(OPT,*)
    enddo
    write(OPT,FMT_DBLINE)

  end subroutine DQMC_Print_HeaderFT

  !--------------------------------------------------------------------!

end module DQMC_GEOM_WRAP
