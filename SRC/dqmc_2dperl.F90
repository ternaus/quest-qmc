module DQMC_2DPERL

  use DQMC_UTIL
  use DQMC_CFG
  use DQMC_STRUCT
  use DQMC_HUBBARD
  ! Chia-Chen: 09/06/2011
  ! added in order to access the number of processors
  use DQMC_MPI
  implicit none 

  ! 
  ! This module defines subroutines to initialize the data structure of 
  ! a two-dimensional periodic rectangular lattice (2DPerl).
  !     
  !  Subroutine List
  !  ===============
  !    DQMC_Read2PERL : read in data for 2D rectangular 
  !                                    lattice.
  !    DQMC_2DREC(nx, ny, S) : construct data for a 2D rectangular lattice.   
  !

  integer, parameter :: IMP_TRIANGLE  = 1
  integer, parameter :: IMP_RECTANGLE = 2
  
contains
  
  !---------------------------------------------------------------------!

  subroutine DQMC_Comp_2DPerl
    !
    ! Purpose
    ! =======
    !    This subroutine reads in data for a 2D rectangular lattice,
    !    and the calls DQMC_2DREC to construct the lattice structure.
    !    Basically, the only parameters needed are Nx and Ny, which
    !    the size of the lattice along the x-cord and the y-cord.
    !
    ! Arguments
    ! =========
    !
    integer :: OPT     ! Input/output handle

    ! ... Local scalar ...
    type(config)  :: cfg
    type(Hubbard) :: Hub          ! Hubbard model
    integer       :: nx, ny

    integer       :: i, j, k, nBin, nIter, slice
    real(wp)      :: randn(1)
    character(len=60) :: ofile

    ! ... Executable ...

    ! Read in run parameters from input file
    call DQMC_Read_Config(cfg)

    ! Initialize the geometry
    call CFG_Get(cfg, "nx", nx) 
    call CFG_Get(cfg, "ny", ny)
    call DQMC_Init_2DPerl(nx, ny, Hub%S, IMP_TRIANGLE)

    ! Fetch additional parameters
    call CFG_Get(cfg, "ofile", ofile)   ! output file name

    ! Initialize the Hubbard data structure
    call DQMC_Hub_Config(Hub, cfg)

    ! Execution MC loop
    ! Warmup sweep
    do i = 1, Hub%nWarm
       if (mod(i,10)==0) write(*,'(A,i6,1x,i3)')' Warmup Sweep, nwrap  : ', i, Hub%G_up%nwrap
       ! The second parameter means no measurement should be made.
       call DQMC_Hub_Sweep(Hub, NO_MEAS0)
       call DQMC_Hub_Sweep2(Hub, Hub%nTry)
    end do

    ! Measurement sweeps. Note that we divide up nPass measurement sweeps into nBins, 
    ! each having nPass/nBin sweeps. This is to reduce autocorrelation between measurements.
    nBin   = Hub%P0%nBin
    nIter  = Hub%nPass / Hub%tausk / nBin
    if (nIter > 0) then
       do i = 1, nBin
          do j = 1, nIter
             do k = 1, Hub%tausk
                call DQMC_Hub_Sweep(Hub, NO_MEAS0)
                call DQMC_Hub_Sweep2(Hub, Hub%nTry)
             enddo

             ! Fetch a random imaginary time slice for measurement 
             call ran0(1, randn, Hub%seed)
             slice = ceiling(randn(1)*Hub%L)

             write(*,'(a,3i6)') ' Measurement Sweep, bin, iter, slice : ', i, j, slice
             call DQMC_Hub_Meas(Hub, slice)
          end do

          ! Accumulate results for each bin
          call DQMC_Phy0_Avg(Hub%P0)

          if (Hub%meas2) then
             if(Hub%P2%diagonalize)then
               call DQMC_Phy2_Avg(Hub%P2, Hub%S)
             else
               call DQMC_Phy2_Avg(Hub%P2, Hub%S%W)
             endif
          end if
       end do
    else
      write(*,*) " Error : The number of measurement sweeps npass/(nbin*tausk) in each bin is less than zero!"
      write(*,*) "         Reset 'npass', 'nbin', and 'tausk' in the input and start over"
      stop
    endif



    !Compute average and error
    call DQMC_Phy0_GetErr(Hub%P0)
    if (Hub%meas2) then
       call DQMC_Phy2_GetErr(Hub%P2)
    end if

    ! Prepare output file
    call DQMC_open_file(adjustl(trim(ofile))//'.out', 'unknown', OPT)

    ! Print computed results
    call DQMC_Hub_Print(Hub, OPT)

    write(*,*)'Done Printing'

    ! Clean up the used storage
    call DQMC_Hub_Free(Hub)
    call DQMC_Config_Free(cfg)
    
  end subroutine DQMC_Comp_2DPerl
  
  !---------------------------------------------------------------------!

  subroutine DQMC_Init_2DPerl(nx, ny, S, IMP)    
    !
    ! Purpose
    ! =======
    !    This subroutine constuctures data structure for a 2D 
    !    periodic rectangular lattice.  
    !
    ! Details
    ! =======
    !    For a Nx*Ny 2D rectangular lattice.
    !
    !    1. The sites are numbered from 1, which is the site on the 
    !       south-west corner. The numbering is row major, which means
    !       it increases along x direction first (left to right) and 
    !       then y-direction (bottom up).
    !
    !    2. Adjacency (T) has exact 4 elements per site: left
    !       right, up and down. Since the adjacency is cyclic, 
    !       the code has had spacial treat for boundary sites.
    !       *** The way it computed here is not satisfied the 
    !           'checkboard' order. should change it later.
    !
    !    3. The number of unique distance is computed as follows.
    !       Let long = max(Nx, Ny) and short = min(Nx, Ny).
    !       The distinct distance sites form a trapezoid
    !       The bottom is (long/2+1), the top is ((long-short)/2+1)
    !       and the height is (short/2+1). Therefore, 
    !       
    !           nClass = (short/2+1)*(long-short/2+2)/2
    !       
    !    4. The phase (P) is computed in the rules that
    !       (a) adjacent sites have opposite phase.
    !       (b) site 1 is phased +.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)         :: nx, ny  ! dimension of the lattice
    type(Struct), intent(inout) :: S       ! Struct
    integer, intent(in)         :: IMP     

    ! ... local vars ...
    integer  :: n                ! Order of matrix T and D 
    integer  :: i, j, jx, jy     ! Loop iterator
    integer  :: up, dn, lt, rt   ! adj and neighbor
    integer  :: ne, nw, se, sw   ! 

    integer  :: idx              ! 
    integer  :: tmp(nx*ny, nx*ny)
    real(wp), pointer :: cord(:,:)
    
    ! ... parameters ...
    integer, parameter :: NADJ  = 4  ! Number of adjacencies
    integer, parameter :: NNBR  = 9  ! Number of neighbors
    
    real(wp), parameter  :: TWOPI     = 6.283185307179586

    ! ... Executable ...

    n   = nx*ny
    S%nSite = n
    write(S%name,'(A,I3,A,I3,A,I5)') &
         "2D Periodic Lattice; Nx=", &
         nx, "; Ny=", ny, "; total sites=", S%nSite
    
    ! Compute all distinct sites
    S%nWave = 9

    ! memory allocation
    allocate(S%P(n))
    allocate(S%W(NNBR, S%nWave))
    allocate(S%dim(2))
    allocate(S%map(n))

    S%dim(1) = nx
    S%dim(2) = ny
    S%Map    = 1
    S%nGroup = 1

    allocate(S%gf_phase(n,n))
    S%gf_phase=1

    ! Build adjacent matrix for checkerboard method
    S%n_t = 1
    tmp = 0
    do j = 1, ny
       do i = 1, nx
          idx = (j-1)*nx+i
          
          ! up link
          up  = idx + nx
          if (up .gt. n) then
             up = up - n
          end if
          
          ! down link
          dn = idx - nx
          if (dn .le. 0) then
             dn = dn + n
          end if
          
          ! left link
          if (i .eq. 1) then
             lt = idx + nx - 1
          else
             lt = idx - 1
          end if
             
          ! right link
          if (i .eq. nx) then
             rt = idx - nx + 1
          else
             rt = idx + 1
          end if

          tmp(idx, up) = 1
          tmp(idx, dn) = 1
          tmp(idx, lt) = 1
          tmp(idx, rt) = 1
       end do
    end do

    call DQMC_CCS_Compress(n, n*NADJ, tmp, S%T)

    ! Build neighboring matrix
    !        +---+---+---+
    !        | 7 | 4 | 1 |
    !        +---+---+---+
    !        | 8 | 5 | 2 |
    !        +---+---+---+
    !        | 9 | 6 | 3 |
    !        +---+---+---+
    !
    S%n_b = 9
    tmp = 0
    do j = 1, ny
       do i = 1, nx
          idx = (j-1)*nx+i
          
          ! up link
          up  = idx + nx
          if (up .gt. n) then
             up = up - n
          end if
          
          ! down link
          dn = idx - nx
          if (dn .le. 0) then
             dn = dn + n
          end if
          
          ! left link
          if (i .eq. 1) then
             lt = idx + nx - 1
          else
             lt = idx - 1
          end if
             
          ! right link
          if (i .eq. nx) then
             rt = idx - nx + 1
          else
             rt = idx + 1
          end if

          ! northeast link = up's right
          if (i .eq. nx) then
             ne = up + 1 - nx
          else
             ne = up + 1
          end if
          
          ! northwest link = up's left
          if (i .eq. 1) then
             nw = up - 1 + nx
          else
             nw = up - 1
          end if

          ! southeast link = down's right
          if (i .eq. nx) then
             se = dn + 1 - nx
          else
             se = dn + 1
          end if
          
          ! southwest link = down's left
          if (i .eq. 1) then
             sw = dn - 1 + nx
          else
             sw = dn - 1
          end if

          ! fill
          tmp(idx, ne) = 1
          tmp(idx, rt) = 2
          tmp(idx, se) = 3
          tmp(idx, up) = 4
          tmp(idx,idx) = 5
          tmp(idx, dn) = 6
          tmp(idx, nw) = 7
          tmp(idx, lt) = 8
          tmp(idx, sw) = 9
       end do
    end do

    call DQMC_CCS_Compress(n, n*NNBR, tmp, S%B)

    ! build up the distance matrix.
    if (IMP .eq. IMP_TRIANGLE) then
       call DQMC_2DPerl_DC_Imp1(n, nx, ny, S, cord)
    else
       call DQMC_2DPerl_DC_Imp2(n, nx, ny, S, cord)
    end if

    ! Initialize phase matrix
    S%P(1:nx:2) = -1
    S%P(2:nx:2) = 1
    do i = 2, ny, 2
       S%P((i-1)*nx+1:i*nx) =  -S%P(1:nx)
    end do

    do i = 3, ny, 2
       S%P((i-1)*nx+1:i*nx) = S%P(1:nx)
    end do

    ! Make wave matrix
    S%W(:,1) = (/ ZERO, ZERO, ZERO, ZERO,  ONE, ZERO, ZERO, ZERO, ZERO/)
    S%W(:,2) = (/ ZERO, HALF, ZERO, HALF, ZERO, HALF, ZERO, HALF, ZERO/)
    S%W(:,3) = (/ ZERO,-HALF, ZERO, HALF, ZERO, HALF, ZERO,-HALF, ZERO/)
    S%W(:,4) = (/ HALF, ZERO, HALF, ZERO, ZERO, ZERO, HALF, ZERO, HALF/)
    S%W(:,5) = (/-HALF, ZERO, HALF, ZERO, ZERO, ZERO, HALF, ZERO,-HALF/)
    S%W(:,6) = (/ ZERO, ZERO, ZERO,-HALF, ZERO, HALF, ZERO, ZERO, ZERO/)
    S%W(:,7) = (/ ZERO,-HALF, ZERO, ZERO, ZERO, ZERO, ZERO, HALF, ZERO/)
    S%W(:,8) = (/-HALF, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, HALF/)
    S%W(:,9) = (/ ZERO, ZERO,-HALF, ZERO, ZERO, ZERO, HALF, ZERO, ZERO/)

    ! label for each wave function
    allocate(S%wlabel(9))
    write(S%wlabel(1),"(a20)") "  S-Wave : "
    write(S%wlabel(2),"(a20)") " SX-Wave : "
    write(S%wlabel(3),"(a20)") "  D-Wave : "
    write(S%wlabel(4),"(a20)") "SXX-Wave : "
    write(S%wlabel(5),"(a20)") "DXX-Wave : "
    write(S%wlabel(6),"(a20)") " PX-Wave : "    
    write(S%wlabel(7),"(a20)") " PY-Wave : "    
    write(S%wlabel(8),"(a20)") "PXY-Wave : "    
    write(S%wlabel(9),"(a20)") "PYX-Wave : "


    ! Fourier Transformation
    allocate(S%FT(S%nClass, S%nClass))
    S%FT = ZERO
    
    do i = 1, S%nClass
       do jx = 0, nx - 1
          do jy = 0, ny - 1
             j = S%D(jy*nx + jx + 1, 1)
             S%FT(i,j) = S%FT(i,j) + &
                  cos(TWOPI*(cord(i,1)*jx/nx+cord(i,2)*jy/ny))

          end do
       end do
    end do

    ! enable the flag
    S%checklist= .true.
    deallocate(cord)

  end subroutine DQMC_INIT_2DPERL

  !---------------------------------------------------------------------!

  subroutine DQMC_2DPerl_DC_Imp1(n, nx, ny, S, cord)    
    !
    ! Purpose
    ! =======
    !    This subroutine implements distance classification, which
    !    use compact stragety
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)         :: n, nx, ny  ! dimension of the lattice
    type(Struct), intent(inout) :: S          ! Struct
    real(wp), pointer           :: cord(:,:) 

    ! ... local vars ...
    integer  :: L(nx*ny,nx*ny)   ! distance table
    integer  :: long, short      ! Used in computing nClass
    integer  :: i, j, ix, iy, jx, jy, dx, dy, idx

    ! site i to i is a special case, do not compute it
    long = max(nx, ny)
    short = min(nx, ny)
    S%nClass = (short/2+1)*(long-short/2+2)/2

    allocate(S%D(n,n))
    allocate(S%F(S%nClass))
    allocate(S%clabel(S%nClass))
    allocate(cord(S%nClass,2))

    S%D = 0
    S%F = 0

    S%F(1) = n
    write(S%clabel(1),FMT_POINT) 0, 0

    !! using lookup table (maybe hashing table )
    L = 0
    idx = 1 ! the first one is a special one
    do i = 1, n
       !! compute the index of i
       ix = mod(i-1, nx)+1
       iy = (i-1)/nx + 1

       !! initial the index of j
       do j = i+1, n
          !! compute the index of j
          jx = mod(j-1, nx)+1
          jy = (j-1)/nx + 1

          !! compute the distance
          dx = abs(ix-jx)
          dx = min(dx,nx-dx)
          dy = abs(iy-jy)
          dy = min(dy,ny-dy)
          long = max(dx, dy) + 1
          short = min(dx, dy) + 1
          
          ! not found
          if (L(long,short) .eq. 0) then

             idx = idx + 1             
             L(long, short) = idx
             S%D(i,j) = idx
             write(S%clabel(idx),FMT_POINT) long-1, short-1
             cord(idx, 1) = long  - 1
             cord(idx, 2) = short - 1

          else ! found

             S%D(i,j) = L(long, short)

          end if

          ! matrix D is symmetric
          S%D(j,i) = S%D(i,j)

          ! increase count by 2
          S%F(S%D(i,j)) = S%F(S%D(i,j)) + 2

       end do

       ! site i to i
       S%D(i,i) = 1
    end do

  end subroutine DQMC_2DPERL_DC_IMP1

  !---------------------------------------------------------------------!

  subroutine DQMC_2DPerl_DC_Imp2(n, nx, ny, S, cord)    
    !
    ! Purpose
    ! =======
    !    In this implementation of distance classication, the classes
    !    are in a rectangular form. 
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)         :: n, nx, ny  ! dimension of the lattice
    type(Struct), intent(inout) :: S          ! Struct
    real(wp), pointer           :: cord(:,:)

    ! ... parameters ...
    real(wp), parameter :: TWOPI = 6.283185307179586

    ! ... local vars ...
    integer  :: L(n,n)                        ! distance table
    integer  :: i, j                          ! loop iterators
    integer  :: ix, iy, jx, jy, dx, dy, idx   ! indices

    ! site i to i is a special case, do not compute it

    S%nClass = (nx/2+1)*(ny/2+1)

    allocate(S%D(n,n))
    allocate(S%F(S%nClass))
    allocate(S%clabel(S%nClass))
    allocate(cord(S%nClass,2))

    S%D  = 0
    S%F  = 0

    S%F(1) = n
    write(S%clabel(1),FMT_POINT) 0, 0

    !! using lookup table (maybe hashing table )
    L = 0
    idx = 1 ! the first one is a special one
    cord(1,1) = ZERO
    cord(1,2) = ZERO
    do i = 1, n
       !! compute the index of i
       ix = mod(i-1, nx)+1
       iy = (i-1)/nx + 1

       !! initial the index of j
       do j = i+1, n
          !! compute the index of j
          jx = mod(j-1, nx)+1
          jy = (j-1)/nx + 1

          !! compute the distance
          dx = abs(ix-jx)
          dx = min(dx,nx-dx)+1
          dy = abs(iy-jy)
          dy = min(dy,ny-dy)+1
          
          ! not found
          if (L(dx,dy) .eq. 0) then

             ! Creat a new node
             idx = idx + 1             
             L(dx, dy) = idx
             S%D(i,j) = idx
             write(S%clabel(idx),FMT_POINT) dx-1, dy-1

             ! Build a new row of COS table
             cord(idx, 1) = dx - 1
             cord(idx, 2) = dy - 1
          else ! found

             S%D(i,j) = L(dx, dy)

          end if

          ! matrix D is symmetric
          S%D(j,i) = S%D(i,j)

          ! increase count by 2
          S%F(S%D(i,j)) = S%F(S%D(i,j)) + 2

       end do

       ! site i to i
       S%D(i,i) = 1

    end do

  end subroutine DQMC_2DPERL_DC_IMP2

end module DQMC_2DPERL
 
