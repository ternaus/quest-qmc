module DQMC_STRUCT
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_CFG
  implicit none 

  ! 
  ! This module defines the data type for describing underline lattice
  ! struture used in DQMC. Basically, it contains six arrays. 
  ! Suppose the lattice has N sites. 
  !
  !   1. Adjacenty (T) is an NxK array, where K is the
  !      maximum number of adjacent sites (=max_adj). T(0,i) records
  !      the number of adjacency of site i. T(1:T(0,i),i) records
  !      the sites that are adjacent to site i. The order of 
  !      these adjacency information is important if the 'checkboard
  !      method' is used to construct Green's function. See DQMC_Gfun
  !      for more details.
  !
  !   2. Distance (D) is an NxN symmetric matrix. D(i,j) records the
  !      distance information between site i and site j. However, the
  !      data stored in D is not 'real distance' in any mathematical 
  !      sense. It is just an index of a distance table, which classifies
  !      all the possible distance in the lattice. In other words,
  !      if two elements in D have the same id, they have the same 
  !      distance. The total number of possible distance is nClass.
  !
  !   3. Frequency (F) is a length nClass vector, which stores the 
  !      counts for each unique distance. For example, F(i)=5 means
  !      there are 5 elements in D whose distance class is indexed i.
  !
  !   4. Phase (P) is also an NxN matrix, which is used in computing
  !      some phiscal measurements and Green's function. See DQMC_Phy0
  !      and DQMC_Phy2 for more details.
  !
  !   5. Neighborhood (N) is a NxB array, where B is the maximum
  !      number of neighbors (=max_nbr).  
  !      N(0,i) records the number of neighbors
  !      of site i. N(1:N(0,i),i) is a list of sites neighboring to
  !      site i. The order of neighbor should be consistant.
  !
  !   6. Wave (W) is a BxB matrix where B is the maximum number of
  !      neighbors. This function is used in computing pair measurements.
  !      See DQMC_Phy2 for more details.
  !
  !  Data Type
  !  =========
  !

  ! Compressed Column Storage for sparse matrices.
  ! This is used for hopping and neighbor matrices.
  ! The detail of CCS can be found on http://www.netlib.org/
  type CCS
     integer :: n                  ! dimension of matrix
     integer :: nnz                ! number of nonzeros
     integer, pointer :: A(:)      ! nonzeros elements, dim = nnz
     integer, pointer :: row(:)    ! row index, dim = nnz
     integer, pointer :: cstart(:) ! start position of columns, dim = n
  end type CCS

  ! The entire geometry related data structures.
  integer, parameter :: N_CHECKLIST   = 8
  integer, parameter :: STRUCT_INIT   = 1
  integer, parameter :: STRUCT_DIM    = 2
  integer, parameter :: STRUCT_ADJ    = 3
  integer, parameter :: STRUCT_CLASS  = 4
  integer, parameter :: STRUCT_WAVE   = 5
  integer, parameter :: STRUCT_BOND   = 6
  integer, parameter :: STRUCT_PHASE  = 7
  integer, parameter :: STRUCT_FT     = 8

  integer, parameter :: label_len     = 47
  integer, parameter :: gname_len     = 80

  type Struct
     integer           :: nSite         ! number of sites 
     integer           :: nCell
     character(gname_len):: name        ! Name of the structure
     integer, pointer  :: dim(:)        ! dim of the geometry

     integer           :: n_t           ! number of hopping types
     type(CCS)         :: T             ! hopping matrix

     integer           :: nClass        ! number of unique distance
     integer, pointer  :: D(:,:)        ! Distance 
     integer, pointer  :: F(:)          ! Frequency 
     integer           :: nGroup
     integer, pointer  :: map(:)        ! site classification
     integer, pointer  :: gf_phase(:,:) 
     integer, pointer  :: chi_phase(:,:) 
                              

     integer           :: n_b           ! number of neighbors types
     type(CCS)         :: B             ! Neighborhood matrix
     integer           :: nClass_b
     integer, pointer  :: class_b(:,:)  ! pair of bonds class
     integer, pointer  :: size_b(:)     ! size of bonds classes

     type(CCS)         :: ckb
     integer           :: nckb

     real(wp), pointer :: W(:,:)        ! Wave
     integer           :: nWave
     integer           :: nirrep        ! number of irreducible representations
     integer, pointer  :: wrepr(:)      ! symmetry group of wave  
     integer           :: nwclass       ! number of classes in the group
     integer, pointer  :: wclass(:)     ! symmetry group of wave  

     character(label_len), pointer :: clabel(:) ! Label for the distance table. 
     character(label_len), pointer :: wlabel(:) ! Label for the distance table. 

     real(wp), pointer :: P(:)          ! Phase

     real(wp), pointer :: FT(:,:)       ! FT matrix for Green
          
     logical :: checklist(N_CHECKLIST)

  end type Struct


  integer, parameter :: N_GEO_PARAM = 15
  integer, parameter :: GEMO_B      =  1
  integer, parameter :: GEMO_D      =  2
  integer, parameter :: GEMO_FT     =  3
  integer, parameter :: GEMO_P      =  4
  integer, parameter :: GEMO_T      =  5
  integer, parameter :: GEMO_W      =  6
  integer, parameter :: GEMO_CLabel =  7
  integer, parameter :: GEMO_dim    =  8
  integer, parameter :: GEMO_nClass =  9
  integer, parameter :: GEMO_nSite  = 10
  integer, parameter :: GEMO_nWave  = 11
  integer, parameter :: GEMO_n_b    = 12
  integer, parameter :: GEMO_n_t    = 13
  integer, parameter :: GEMO_Name   = 14
  integer, parameter :: GEMO_wLabel = 15

  character(len=*), parameter :: GEO_PARAM(N_GEO_Param) =  &
       &(/"B      ", &
       &  "D      ", &
       &  "FT     ", &
       &  "P      ", &
       &  "T      ", &
       &  "W      ", &
       &  "cLabel ", &
       &  "dim    ", &
       &  "nClass ", &       
       &  "nSite  ", &
       &  "nWave  ", &
       &  "n_b    ", &
       &  "n_t    ", &
       &  "name   ", &
       &  "wLabel "/)

contains
  
  !---------------------------------------------------------------------!

  subroutine DQMC_CCS_Free(sparA)
    !
    ! Purpose
    ! =======
    !    This subroutine deallocates arrays in CCS
    !
    ! Arguments
    ! =========
    !
    type(CCS), intent(inout) :: sparA

    ! ... Executable ...

    if (sparA%nnz .gt. 0) then
       deallocate(sparA%A)
       deallocate(sparA%row)
       deallocate(sparA%cstart)
    end if

    sparA%nnz = 0
    sparA%n   = 0

  end subroutine DQMC_CCS_Free

  !---------------------------------------------------------------------!

  subroutine DQMC_CCS_Print(sparA, OPT)
    !
    ! Purpose
    ! =======
    !    This subroutine deallocates arrays in CCS
    !
    ! Arguments
    ! =========
    !
    type(CCS), intent(in) :: sparA
    integer, intent(in)   :: OPT

    ! ... Local ...
    integer :: i, j

    ! ... Executable ...

    do i = 1, sparA%n
       do j = sparA%cstart(i), sparA%cstart(i+1)-1
          write(OPT, "(i5, i5, i5)") sparA%row(j), i, sparA%A(j)
       end do
    end do

  end subroutine DQMC_CCS_Print
  
  !---------------------------------------------------------------------!

  subroutine DQMC_CCS_Compress(n, nnz, A, sparA)
    !
    ! Purpose
    ! =======
    !    This subroutine converts a dense A to a sparse A in CCS format.
    !    If nnz is unknown, passed nnz=-1
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n, nnz
    integer, intent(in)      :: A(:,:)
    type(CCS), intent(inout) :: sparA

    ! ... Local varaibles ...
    integer :: i, j, k

    ! ... Executable ...

    ! two pass process
    if (nnz .lt. 0) then
       sparA%nnz = 0
       do i = 1, n
          do j = 1, n
             if (A(j,i) .ne. 0) then
                sparA%nnz = sparA%nnz + 1
             end if
          end do
       end do
    else
       sparA%nnz = nnz
    end if

    ! allocate space for sparse A
    sparA%n   = n
    allocate(sparA%A(sparA%nnz))
    allocate(sparA%row(sparA%nnz))
    allocate(sparA%cstart(n+1))

    ! fillin data
    k = 1
    do i = 1, n
       sparA%cstart(i) = k
       do j = 1, n
          if (A(j,i) .ne. 0) then
             sparA%row(k) = j
             sparA%A(k)   = A(j,i)
             k = k + 1
          end if
       end do
    end do
    
    ! validate 
    if (k .ne. sparA%nnz+1) then
       print *, k, sparA%nnz+1
       call DQMC_Error("The passing in nnz does not match actual &
            & number of nonzeros", 0)
    end if
    
    sparA%cstart(n+1) = k

  end subroutine DQMC_CCS_Compress

  !---------------------------------------------------------------------!

  subroutine DQMC_CCS_Fill(n, A, sparA)
    !
    ! Purpose
    ! =======
    !    This subroutine converts a sparse A to a dense A.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n
    integer, intent(inout)   :: A(:,:)
    type(CCS), intent(in)    :: sparA

    ! ... Local varaibles ...
    integer :: i, j

    ! ... Executable ...

    A = 0
    do i = 1, n
       do j = sparA%cstart(i), sparA%cstart(i+1)-1
          A(sparA%row(j),i) =  sparA%A(j)
       end do
    end do

  end subroutine DQMC_CCS_Fill

  !---------------------------------------------------------------------!

  subroutine DQMC_Struct_Free(S)
    !
    ! Purpose
    ! =======
    !    This subroutine deallocates arrays in S
    !
    ! Arguments
    ! =========
    !
    type(struct), intent(inout) :: S         ! Struct

    ! ... Executable ...

    if (S%checklist(STRUCT_INIT)) then
       if (S%checklist(STRUCT_DIM))   deallocate(S%dim)
       if (S%checklist(STRUCT_ADJ))   call DQMC_CCS_Free(S%T)
       if (S%checklist(STRUCT_CLASS)) deallocate(S%D, S%F, S%cLabel)
       if (S%checklist(STRUCT_BOND))  call DQMC_CCS_Free(S%B)
       if (S%checklist(STRUCT_WAVE))  deallocate(S%W, S%wLabel)
       if (S%checklist(STRUCT_PHASE)) deallocate(S%P)
       if (S%checklist(STRUCT_FT))    deallocate(S%FT)

    end if
    S%checklist = .false.

  end subroutine DQMC_Struct_Free

  !---------------------------------------------------------------------!

  subroutine DQMC_Geom_Read_Def(S, gfile, tableFormat)
    !
    ! Purpose
    ! =======
    !    This subrotine reads in parameters from a config file.
    !
    ! Arguments
    ! =========
    !
    type(struct), intent(inout)  :: S            ! geometry structure
    character(*), intent(in)     :: gfile        ! Input file handle
    logical, intent(out)         :: tableFormat

    ! ... Local Variable ...
    integer                :: stat, i, j, line, pos, idx, cnt, row, col, ielm
    real(wp)               :: relm
    character(len=llen)    :: str, attr, val
    logical                :: found, def_n, def_w, def_b, def_c
    integer, parameter     :: funit = 11
    integer, allocatable   :: tmp(:,:)
    ! ... Executable ...

    ! satinize
    S%checklist = .false.
    
    ! open file
    inquire(file=gfile, exist=found)
    if (.not. found) then
       call DQMC_Error("cannot open geom def file "//gfile, 0)
    end if
    open(unit = funit, file = gfile)

    !check whether this is in the "table" format
    tableFormat=.true.
    do
     read(funit,'(A)',iostat=stat)str
     if (stat .eq. 0) then
      pos=index(str,'#HAMILT')
      if(pos .ne. 0)then
       tableFormat=.false.
       write(*,'(A)') ' Detected Geometry Free Format.'
       exit
      endif
     else
      rewind(funit)
      exit
     endif
    enddo

    if(.not.tableFormat)return
    ! read geom def
    ! for fast access, sort records by name
    ! using insertion sort
    stat = STAT_COMMENT
    def_n = .false.
    def_b = .false.
    def_w = .false.
    def_c = .false.
    nullify(S%wlabel)
    nullify(S%clabel)
    line = 1
    S%name = "Not defined"
    do while (stat .ne. STAT_EOF)
       call DQMC_ReadLn(str, funit, stat)

       ! read in a geom definition 
       if (stat .eq. STAT_NORMAL) then
          pos = scan(str, SEPARAT, .false.)
          if (pos .ne. 0) then
             ! read name and data 
             attr = adjustl(str(1:pos-1))
             val  = adjustl(str(pos+1:llen))
             
             ! search the name (linear search)
             found = .false.
             do idx = 1, N_GEO_PARAM
                if (attr .eq. GEO_PARAM(idx)) then
                   found = .true.
                   exit
                end if
             end do
             
             ! found the name
             if (found) then
                select case(idx)
                case(GEMO_B     ) ! n*n
                   if (def_n) then
                      read(val,*) cnt
                      tmp = 0
                      do i = 1, cnt
                         read(funit, *) row, col, ielm
                         tmp(row, col) = ielm
                         line = line + 1
                      end do
                      ! compress
                      call DQMC_CCS_Compress(S%nSite, cnt, tmp, S%B)
                   else
                      call DQMC_Error("nSite &
                           & need be defined before Bond", 0)
                   end if
                   S%checklist(STRUCT_BOND) = .true.

                case(GEMO_D     ) ! n*n    
                   if (def_n) then
                      read(val,*) cnt
                      allocate(S%D(S%nSite, S%nSite))
                      S%D = 0
                      do i = 1, cnt
                         read(funit, *) row, col, ielm
                         S%D(row, col) = ielm
                         line = line + 1
                      end do
                   else
                      call DQMC_Error("nSite &
                           & need be defined before D-class", 0)
                   end if
                   S%checklist(STRUCT_CLASS) = .true.

                case(GEMO_FT    ) ! nclass*nclass   
                   if (def_c) then
                      read(val,*) cnt
                      allocate(S%FT(S%nClass, S%nClass))
                      do i = 1, cnt
                         read(funit, *) row, col, relm
                         S%FT(row, col) = relm
                         line = line + 1
                      end do
                   else
                      call DQMC_Error("nClass &
                           & need be defined before Fourier Trans", 0)
                   end if
                   S%checklist(STRUCT_FT) = .true.
                   
                case(GEMO_P     ) ! check n   
                   read(val,*) cnt
                   allocate(S%P(cnt))
                   do i = 1, cnt
                      read(funit, *) S%P(i)
                      line = line + 1
                   end do
                   S%checklist(STRUCT_PHASE) = .true.

                case(GEMO_T     ) ! n*n
                   if (def_n) then
                      read(val,*) cnt
                      tmp = 0
                      do i = 1, cnt
                         read(funit, *) row, col, ielm
                         tmp(row, col) = ielm
                         line = line + 1
                      end do
                      ! compress
                      call DQMC_CCS_Compress(S%nSite, cnt, tmp, S%T)
                   else
                      call DQMC_Error("nSite &
                           & need be defined before Hopping", 0)
                   end if
                   S%checklist(STRUCT_ADJ) = .true.

                case(GEMO_W     ) ! n_b*nwave
                   if (def_w .and. def_b) then
                      read(val,*) cnt
                      allocate(S%W(S%n_b, S%nWave))
                      S%W = ZERO
                      do i = 1, cnt
                         read(funit, *) row, col, relm
                         S%W(row, col) = relm
                         line = line + 1
                      end do
                   else
                      call DQMC_Error("n_b and nWave &
                           & need be defined before Wave funs", 0)
                   end if
                   S%checklist(STRUCT_WAVE) = .true.

                case(GEMO_CLabel)
                   read(val,*) cnt
                   allocate(S%clabel(cnt))
                   do i = 1, cnt
                      read(funit, '(A)') S%clabel(i)
                      line = line + 1
                   end do

                case(GEMO_dim   )
                   read(val, *) cnt
                   allocate(S%dim(cnt))
                   do i = 1, cnt
                      read(funit, *) S%dim(i)
                      line = line + 1
                   end do
                   S%checklist(STRUCT_DIM) = .true.

                case(GEMO_nClass)
                   read(val, *) S%nClass
                   def_c = .true.
                   
                case(GEMO_nSite )
                   read(val, *) S%nSite
                   def_n = .true.
                   allocate(tmp(S%nSite, S%nSite))

                case(GEMO_nWave )
                   read(val, *) S%nWave
                   def_w = .true.

                case(GEMO_n_b   )
                   read(val, *) S%n_b

                case(GEMO_n_t   )
                   read(val, *) S%n_t
                   def_b = .true.

                case(GEMO_Name  )
                   read(val, *) S%name

                case(GEMO_wLabel)
                   read(val,*) cnt
                   allocate(S%wlabel(cnt))
                   do i = 1, cnt
                      read(funit, '(A)') S%wlabel(i)
                      line = line + 1
                   end do
                end select

             else
                call DQMC_Warning("Warning: unknown geom input:"//trim(str), 1) 
             end if

          else
             call DQMC_Warning("cannot recog geom input line :", line)
          end if

       end if

       line = line + 1
    end do

    ! do some simple verification
    if (.not. def_n) then
       call DQMC_Error("nSite must be defined.", 0)
    end if

    if (S%checklist(STRUCT_PHASE)) then
       if (size(S%P) .ne. S%nSite) then
          call DQMC_Error("The length of phase assignment is not &
               & equal to nSite", 0)
       end if
    end if

    ! default class assignment
    if (.not. S%checklist(STRUCT_CLASS)) then
       allocate(S%D(S%nSite,S%nSite))
       cnt = 0
       do i = 1, S%nSite
          do j = i, S%nSite
             cnt = cnt + 1
             S%D(i,j) = cnt
             S%D(j,i) = cnt
          end do
       end do
       S%nClass = cnt
       def_c = .true.
       S%checklist(STRUCT_CLASS) = .true.
    end if

    ! nClass = max(D)
    if (.not. def_c) then
      cnt = 1
      do i = 1, S%nSite
         do j = 1, S%nSite
            if (cnt .lt. S%D(i,j)) then
               cnt = S%D(i,j)
            end if
         end do
      end do
      S%nClass = cnt
    end if

    ! discover F and Map
    allocate(S%F(S%nClass))
    allocate(S%map(S%nSite))
    call DQMC_Geom_Discover_Map(S%nSite, S%D, S%F, S%map, S%nGroup)
    call DQMC_Geom_Discover_F(S%nSite, S%D, S%F)

    if (.not. associated(S%clabel)) then
       allocate(S%clabel(S%nClass))
       do i = 1, S%nClass
          write(s%clabel(i), "('Class ',i3)") i
       end do
    end if

    ! wave label
    if (S%checklist(STRUCT_Wave)) then
       if (.not. associated(S%wlabel)) then
          allocate(S%wlabel(S%nWave))
          do i = 1, S%nWave
             write(s%wlabel(i), "('Wave ',i3)") i
          end do
       end if
    end if
    
    ! release temp variable
    deallocate(tmp)
    
    ! finish initialization
    S%checklist(STRUCT_INIT) = .true.
    close(funit)

  end subroutine DQMC_Geom_Read_Def

  !---------------------------------------------------------------------!
  
  subroutine DQMC_Geom_Discover_F(n, D, F)
    !
    ! Purpose
    ! =======
    !    This subroutine generate F
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)    :: n
    integer, intent(in)    :: D(:,:)
    integer, intent(inout) :: F(:)

    ! ... Local ...
    integer :: i, j, k

    ! ... Executable ...
    
    F = 0
    do i = 1, n
       do j = 1, n
          k = D(i,j)
          F(k) = F(k) + 1
       end do
    end do
  end subroutine DQMC_Geom_Discover_F

  !---------------------------------------------------------------------!
  
  subroutine DQMC_Geom_Discover_Map(n, D, F, Map, idx)
    !
    ! Purpose
    ! =======
    !    This subroutine generate F
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)    :: n
    integer, intent(in)    :: D(:,:)
    integer, intent(inout) :: F(:)        ! working space
    integer, intent(inout) :: Map(:)
    integer, intent(out)   :: idx

    ! ... Local ...
    integer :: i, k

    ! ... Executable ...
    
    F   = 0
    map = 0
    idx = 0
    do i = 1, n
       k = D(i,i)
       if (F(k) .eq. 0) then
          idx = idx + 1
          map(i) = idx
          F(k) = idx
       else
          map(i) = F(k)
       end if
    end do
  end subroutine DQMC_Geom_Discover_Map

  !---------------------------------------------------------------------!

  subroutine DQMC_Geom_Print(S, OPT)
    !
    ! Purpose
    ! =======
    !    This subroutine generate F
    !
    ! Arguments
    ! =========
    !    
    type(struct), intent(in)  :: S            ! geometry structure
    integer, intent(in)       :: OPT

    ! ... Local ...
    integer :: i, j, n
    
    ! ... Executable ...
    
    n = S%nSite
    write(OPT, "(a, a)") "name  = ", S%name
    write(OPT, "(a,i5)") "nSite = ", S%nSite

    write(OPT, "(a,i5)") "nWave = ", S%nWave

    if (S%checklist(STRUCT_DIM)) then
       write(OPT, "(a,i5)") "dim   = ", size(S%dim)
       do i = 1, size(S%dim)
          write(OPT, "(i5)") S%dim(i)
       end do
       write(OPT, *)
    end if

    if (S%checklist(STRUCT_ADJ)) then
       write(OPT, "(a,i5)") "n_t   = ", S%n_t
       write(OPT, "(a,i5)") "T     = ", S%T%nnz
       call DQMC_CCS_Print(S%T, OPT)
       write(OPT, *)
    end if

    if (S%checklist(STRUCT_BOND)) then
       write(OPT, "(a,i5)") "n_b   = ", S%n_b
       write(OPT, "(a,i5)") "B     = ", S%B%nnz
       call DQMC_CCS_Print(S%B, OPT)
       write(OPT, *)
    end if

    if (S%checklist(STRUCT_CLASS)) then
       write(OPT, "(a,i5)") "nClass= ", S%nClass
       write(OPT, "(a,i5)") "D     = ", S%nSite*S%nSite
       do i = 1, n
          do j = 1, n
             write(OPT, "(i5, i5, i5)") i, j, S%D(i, j)
          end do
       end do
       write(OPT, *)

       write(OPT, "(a,i5)") "cLabel = ", S%nClass
       do i = 1, S%nClass
          write(OPT, *) S%clabel(i)
       end do
       write(OPT, *)
    end if

    if (S%checklist(STRUCT_PHASE)) then
       write(OPT, "(a,i5)") "P      = ", n
       do i = 1, n
          write(OPT, *) S%P(i)
       end do
       write(OPT, *)
    end if
    
    if (S%checklist(STRUCT_WAVE)) then
       write(OPT, "(a,i5)") "nWave = ", S%nWave
       write(OPT, "(a,i5)") "nBond = ", S%n_b
       write(OPT, "(a,i5)") "W = ", S%nWave*S%n_b
       do i = 1, S%nWave
          do j = 1, S%n_b
             write(OPT, "(i5, i5, f15.8)") j, i, S%W(j, i)
          end do
       end do       
       write(OPT, *)

       write(OPT, "(a,i5)") "wLabel  = ", S%nWave
       do i = 1, S%nWave
          write(OPT, *) S%wlabel(i)
       end do
       write(OPT, *)
    end if
    
    if (S%checklist(STRUCT_FT)) then
       write(OPT, "(a,i5)") "FT = ", S%nClass*S%nClass
       do i = 1, S%nClass
          do j = 1, S%nClass
             write(OPT, *) i, j, S%FT(i, j)
          end do
       end do
    end if

  end subroutine DQMC_Geom_Print

end module DQMC_STRUCT
 
