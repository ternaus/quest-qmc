module DQMC_Geom

  use DQMC_Struct

  implicit none 

  ! 
  ! This module defines subroutines to read geometry definition from files.
  
contains

  !---------------------------------------------------------------------!

  subroutine QUEST_GG_Init(S, prefix)
    !
    ! Purpose
    ! =======
    !    This subroutine reads in geometry data.
    !
    ! Arguments
    ! =========
    !
    character(*), intent(in)    :: prefix     ! input file name
    type(Struct), intent(inout) :: S           ! Struct

    ! ... Local Variable ...
    integer :: n, fh, fh2
    logical :: valid
    

    ! ... Executable ...

    S%checklist = .false.
    S%nSite = 0
    fh = 45
    
    ! neighbor table
    inquire(FILE=trim(prefix)//".neigh", EXIST=valid)
    if (.not. valid) then
       print *, prefix, ".neigh does not exist. "
       stop
    end if

    print *, "Open ", trim(prefix)//".neigh"
    open(fh, FILE = trim(prefix)//".neigh", form='formatted')
    call QUEST_Read_Adj(fh, S)
    close(fh)
    n = S%nSite
    
    ! class table 
    inquire(FILE=trim(prefix)//".class", EXIST=valid)
    if (.not. valid) then
       print *, prefix, ".class does not exist."
       stop
    end if

    ! label table
    fh2 = -1
    inquire(FILE=trim(prefix)//".label", EXIST=valid)
    if (.not. valid) then
       print *, prefix, ".label does not exist. Use defualt labels."
    else
       fh2 = 44
       print *, "Open ", trim(prefix)//".label"
       open(fh2, FILE = trim(prefix)//".label", form='formatted')
    end if

    print *, "Open ", trim(prefix)//".class"
    open(fh, FILE = trim(prefix)//".class", form='formatted')
    call QUEST_Read_Class(fh, fh2, S)
    close(fh)
    if (fh2 .gt. 0) then
       close(fh2)
    end if

    
    ! consistence check
    if (n .ne. S%nSite) then
       print *, "Number of sites are inconsistent for class files."
    end if

    ! phase table
    inquire(FILE=trim(prefix)//".phase", EXIST=valid)
    if (.not. valid) then
       print *, prefix, ".phase does not exist. "
       stop
    end if

    print *, "Open ", trim(prefix)//".phase"
    open(fh, FILE = trim(prefix)//".phase", form='formatted')
    call QUEST_Read_Phase(fh, S)
    close(fh)

    if (n .ne. S%nSite) then
       print *, "Number of sites are not inconsistent for phase files."
    end if

    ! other part of S
    allocate(S%Umap(n))
    allocate(S%mumap(n))
    S%Umap  = 1
    S%muMap = 1

    S%init  = .true. 
  end subroutine QUEST_GG_Init

  !---------------------------------------------------------------------!
  
  subroutine QUEST_Read_Adj(IPT, S)
    !
    ! Purpose
    ! =======
    !    This subroutine reads in adjacent tables.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)         :: IPT     ! input device
    type(Struct), intent(inout) :: S       ! Struct

    ! ... Local Variable ...
    integer  :: n, n_t, max_adj, line, cnt, ios
    integer  :: i, j, b
    integer, allocatable :: tbl(:,:), band_cnt(:)

    ! ... Executable ...

    ! read n, n_t, nmax
    read (IPT, *)  n, n_t
    
    ! initializes data 
    S%nSite = n
    S%n_t = n_t
    
    allocate(S%BS(n,n_t))
    S%BS = 0
    allocate(tbl(n,n))
    allocate(band_cnt(n_t))
    tbl  = 0

    ! read in data
    line = 1
    do 
       read (unit=IPT, iostat=ios, FMT=*)  i, j, b

       ! end of file
       if (ios .ne. 0) then
          exit
       end if

       ! validate input
       if (i.le.0.or.i.gt.n.or.j.le.0.or.j.gt.n.or.b.le.0.or.b.gt.n_t) then
          print *, "input data error", i,j,b, " in line ", line
       end if

       if (tbl(i,j) .ne. 0) then
          print *, "Neighbor ", i,j, " is redefined in line ", line
          stop
       end if

       tbl(i,j) = b

       line = line + 1
    end do

    ! symmetric check
    max_adj = 1
    do i = 1, n
       cnt = 0
       do j = 1, n
          if (tbl(i,j) .ne. tbl(j,i)) then
             print *, i,j,"=",tbl(i,j)," is not equal to ",j,i,"=",tbl(j,i)
             stop
          end if
          ! record band size
          if (tbl(i,j) .gt. 0) then
             S%BS(i, tbl(i,j)) =  S%BS(i, tbl(i,j)) + 1       
             cnt = cnt + 1
          end if
       end do

       ! count the maximum number of adjacency
       if (cnt .gt. max_adj) then
          max_adj = cnt
       end if
    end do
    
    ! build the T table
    allocate(S%T(n,max_adj))
    
    do i = 1, n
       band_cnt = 0
       cnt = 0
       do j = 1, n_t
          band_cnt(j) = cnt
          cnt = cnt + S%BS(i,j)
       end do

       do j = 1, n
          if (tbl(i,j) .ne. 0) then
             band_cnt(tbl(i,j)) = band_cnt(tbl(i,j)) + 1
             S%T(i, band_cnt(tbl(i,j))) = j
          end if
       end do
    end do

    ! checklist
    S%checklist(STRUCT_ADJ) = .true.

  end subroutine QUEST_Read_Adj

  !---------------------------------------------------------------------!

  subroutine QUEST_Read_Class(IPT, IPT2, S)
    !
    ! Purpose
    ! =======
    !    This subroutine reads in class tables.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)         :: IPT, IPT2     ! input device
    type(Struct), intent(inout) :: S             ! Struct

    ! ... Local Variable ...
    integer  :: i, j, c, nClass, n, line, ios, tmp

    ! ... Executable ...

    ! read n, n_t, nmax
    read (unit=IPT, FMT=*)  n, nClass
    
    ! initializes data 
    S%nSite  = n
    S%nClass = nClass
    
    allocate(S%D(n,n))
    allocate(S%F(nClass))
    allocate(S%label(nClass))
    nullify(S%cord)
    S%F = 0
    S%D = 0

    ! read in data
    line = 1
    do 
       read (unit=IPT, iostat=ios, FMT=*)  i, j, c

       ! end of file
       if (ios .ne. 0) then
          exit
       end if

       ! validate input
       if (i.le.0.or.i.gt.n.or.j.le.0.or.j.gt.n.or.c.le.0.or.c.gt.nClass) then
          print *, "input data error", i,j,c, " in line ", line
       end if
       
       if (S%D(i,j) .ne. 0) then
          print *, "Class ", i,j, " is redefined in line ", line
          stop
       end if

       S%D(i,j) = c
       S%F(c) = S%F(c) + 1
       line = line + 1
    end do

    ! check class table
    do i = 1, n
       do j = i, n
          if (S%D(i,j) .eq. 0) then
             print *, "Class of ", i, j, " is not initialized."
             stop
          end if
       end do
    end do

    do i = 1, nClass
       if (S%F(i) .eq. 0) then
          print *, "Class ", i, " = 0."
          stop
       end if
       
    end do

    if (IPT2 .gt. 0) then
       do i = 1, nClass
          read(IPT2, '(a)') S%label(i)
       end do
    else
       do i = 1, nClass
          write(S%label(i), "('class ',i5)") i
       end do
    end if

    ! checklist
    S%checklist(STRUCT_CLASS) = .true.

  end subroutine QUEST_Read_Class

  !---------------------------------------------------------------------!

  subroutine QUEST_Read_Phase(IPT, S)
    !
    ! Purpose
    ! =======
    !    This subroutine reads in phase tables.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)         :: IPT     ! input device
    type(Struct), intent(inout) :: S       ! Struct

    ! ... Local Variable ...
    integer  :: i, j, p, n, line, ios
    integer, allocatable :: phase(:)

    ! ... Executable ...

    ! read n, n_t, nmax
    read (unit=IPT, FMT=*)  n
    
    ! initializes data 
    S%nSite= n
    
    allocate(S%P(n,n))
    allocate(phase(n))
    phase = 0

    line = 1
    do 
       read (unit=IPT,iostat=ios, FMT=*)  i, p

       ! end of file
       if (ios .ne. 0) then
          exit
       end if
       
       ! validate input
       if (i.le.0.or.i.gt.n.or.(p.ne.1.and.p.ne.-1)) then
          print *, "input data error", i,p, " in line ", line
       end if
       
       if (phase(i) .ne. 0) then
          print *, "phase ", i, " is redefined in line ", line
          stop
       end if

       phase(i) = p
       line = line + 1
    end do
    
    ! check class table
    do i = 1, n
       do j = i, n
          S%P(i,j) = phase(i)*phase(j)
       end do
    end do
    
    ! checklist
    S%checklist(STRUCT_PHASE) = .true.
    
  end subroutine QUEST_Read_Phase

  !---------------------------------------------------------------------!

end module DQMC_Geom
