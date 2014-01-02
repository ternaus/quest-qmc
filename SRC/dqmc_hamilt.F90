module DQMC_HAMILT

use DQMC_GEOM_PARAM
use DQMC_LATT
use DQMC_RECLATT
use DQMC_Cfg

implicit none

 type :: Hamiltonian_t
 
   !number of sites having a non-zero J,U,t
   integer                :: nJsites, nUsites, ntsites            
 
   !maximum number of neighbors
   integer                :: maxneig                           
 
   !Number of neighbors of each site (0:nsites-1)
   integer, pointer       :: tnneig(:), Unneig(:), Jnneig(:)
 
   !Sites havind a non-zero J,U,t (0:nsites-1)
   integer, pointer       :: tsite(:), Usite(:), Jsite(:)
 
   !Sites neighboring each site   (0:nsites-1,nsites)
   integer, pointer       :: tneig(:,:), Uneig(:,:), Jneig(:,:)  

   !Pair of sites associated to each hopping throughout the system
   integer, pointer       :: tckb(:,:)
 
   !chemical potential
   real*8                 :: mu_up, mu_dn
 
   !value of U and J for each pair (0:nsites-1, 0:nsites-1)
   real*8, pointer        :: Uv(:,:), Jv(:,:)                   
 
   !values of U and mu inequivalent by symmetry (nlocclass)
   real*8, pointer        :: Uvalue(:), muupvalue(:), mudnvalue(:)
 
   !value of t for each pair (0:nsites-1, 0:nsites-1)
   complex*16, pointer    :: hopup(:,:), hopdn(:,:)
 
   !values of t inequivalent by symmetry (nhopclass)
   complex*16, pointer    :: tupvalue(:), tdnvalue(:)
 
   !wave function phase(not used in QMC)
   complex*16, pointer    :: phase(:)
 
   !number of different hoppings
   integer                :: nhopclass
 
   !number of sites with different U/mu
   integer                :: nlocclass
 
   !class for each site having U/mu (0:nsites-1)
   integer, pointer       :: mylocclass(:)
 
   !hopping class for each pair of neighboring site (0:nsites-1,maxval(tnneig))
   integer, pointer       :: myhopclass(:,:)

   !number of primitive links 
   integer                :: nplink
   !sites involved in primitive links (2, 0:nplink-1)
   !plink(1:2,il) stores the two orbitals involved in the link.
   !plink(1,il) is always contained in the original primitive cell
   integer, pointer       :: plink(:,:)
   !array (5,0:nplink-1) containing the following info:
   ! tlink(1, il) value of up-hopping   on primitive link il
   ! tlink(2, il) value of down-hopping on primitive link il
   ! tlink(3:5,il) = R_2 - R_1
   real(wp), pointer      :: tlink(:,:)
 
 
   logical                :: constructed
   logical                :: neig_found
   logical                :: analyzed
 
 end type

contains

!-----------------------------------------------------------------!

 subroutine free_hamilt(hamilt)
   !
   ! Free the space taken by pointers in Hamilt.
   !
    type(hamiltonian_t), intent(inout) :: hamilt

    if(associated(hamilt%tnneig))     deallocate(hamilt%tnneig)
    if(associated(hamilt%Unneig))     deallocate(hamilt%Unneig)
    if(associated(hamilt%Jnneig))     deallocate(hamilt%Jnneig) 
    if(associated(hamilt%tsite))      deallocate(hamilt%tsite)
    if(associated(hamilt%Usite))      deallocate(hamilt%Usite)
    if(associated(hamilt%Jsite))      deallocate(hamilt%Jsite)  
    if(associated(hamilt%Uvalue))     deallocate(hamilt%Uvalue) 
    if(associated(hamilt%muupvalue))  deallocate(hamilt%muupvalue) 
    if(associated(hamilt%mudnvalue))  deallocate(hamilt%mudnvalue) 
    if(associated(hamilt%tupvalue))   deallocate(hamilt%tupvalue)    
    if(associated(hamilt%tdnvalue))   deallocate(hamilt%tdnvalue)    
    if(associated(hamilt%phase))      deallocate(hamilt%phase)     
    if(associated(hamilt%mylocclass)) deallocate(hamilt%mylocclass) 
    if(associated(hamilt%tneig))      deallocate(hamilt%tneig)
    if(associated(hamilt%Uneig))      deallocate(hamilt%Uneig)
    if(associated(hamilt%Jneig))      deallocate(hamilt%Jneig) 
    if(associated(hamilt%Uv))         deallocate(hamilt%Uv)
    if(associated(hamilt%Jv))         deallocate(hamilt%Jv)   
    if(associated(hamilt%hopup))      deallocate(hamilt%hopup)  
    if(associated(hamilt%hopdn))      deallocate(hamilt%hopdn)  
    if(associated(hamilt%myhopclass)) deallocate(hamilt%myhopclass) 

 end subroutine free_hamilt

!-----------------------------------------------------------------!

 subroutine construct_hamilt(hamilt, lattice, recip_lattice, cfg)
   !
   ! Fill the hamiltonian : Fill all the variables making up 
   ! hamiltonian_t except those whose name end in "class'. Set 
   ! constructed and neig_found to true
   !
    type(lattice_t), target, intent(in) :: lattice
    type(recip_lattice_t),intent(in)    :: recip_lattice
    type(Hamiltonian_t),intent(out)     :: hamilt
    type(config), intent(inout)         :: cfg
 
    integer :: iat, jat, ios, ihop, jhop, il, it, is, js
    integer :: natom, nsites, ntcfg, nline, nhop
    real*8  :: ktwist(rdim), kpoint(rdim), hop3d(rdim)
    real*8  :: tijup, tijdn, U, twisthop
    logical :: ldum, doeshop
    character(len=50) :: string

    real*8, pointer :: pos(:,:), tcfg(:)


    character(len=*), parameter :: mu(2) = (/'mu_up','mu_dn'/)

    pos => lattice%cartpos
    nullify(tcfg)

    if(.not.lattice%constructed) &
       stop'Need to construct lattice before building Hamiltonian'
    if(.not.recip_lattice%initialized) &
       stop'Need to initialize recip_lattice before building Hamiltonian'

    !Set local alias
    natom  = lattice%natom
    nsites = lattice%nsites
    ktwist(1:rdim) = recip_lattice%ktwist(1:rdim)
    kpoint(1:rdim) = recip_lattice%kpoint(1:rdim)

    allocate(hamilt%hopup(0:nsites-1, 0:nsites-1))
    allocate(hamilt%hopdn(0:nsites-1, 0:nsites-1))
    allocate(hamilt%Uv(0:nsites-1, 0:nsites-1)) 
    allocate(hamilt%Jv(0:nsites-1, 0:nsites-1))
    allocate(hamilt%phase(0:nsites-1))

    !Read chemical potential
    do iat = 1, 2
       call CFG_Get(cfg, mu(iat), ntcfg, tcfg)
       if(ntcfg > 1)then
         write(*,'(A)')'WARNING: Only 1st entry for mu in input file is considered.'
       endif
       if (iat == 1) hamilt%mu_up = tcfg(1)
       if (iat == 2) hamilt%mu_dn = tcfg(1)
    enddo
    deallocate(tcfg)
    nullify(tcfg)

    !Find the hamiltonian field
    ldum = move_to_record(INPUT_FIELDS(HAMILT_F),inpunit)

    ! Count lines and hops (hops are off-diagonal matrix elements)
    nline = 0
    nhop  = 0
    do
       read(inpunit,'(A)') string
       read(string,*,iostat=ios) iat, jat, hop3d(1:3), tijup, tijdn, U
       if (ios .ne. 0) exit
       
       !Check label are in range
       if (iat>lattice%natom-1 .or. jat>lattice%natom-1) then
          write(*,*)'One of the atom in hopping is unspecified' 
          stop
       endif

       doeshop = hoptowho(iat, hop3d, jat, lattice) .ne. iat
       doeshop = doeshop .and. (abs(tijup).gt.1.d-6 .or. abs(tijdn).gt.1.d-6)
       if (doeshop) nhop = nhop + 1
       nline = nline + 1
    enddo

    hamilt%nplink = nhop
    allocate(hamilt%tlink(5,0:nhop-1))

    nhop = nhop*lattice%ncell
    allocate(hamilt%tckb(3,0:nhop-1))
    allocate(hamilt%plink(2,0:nhop-1))

    ihop = -1
    jhop = -1
    hamilt%tckb       = 0
    hamilt%hopup(:,:) = 0.d0 
    hamilt%hopdn(:,:) = 0.d0 
    hamilt%Jv(:,:)    = 0.d0 
    hamilt%Uv(:,:)    = 0.d0

    ldum = move_to_record(INPUT_FIELDS(HAMILT_F),inpunit)

    do il = 1, nline

       !Read next line in Hamilt
       read(inpunit,'(A)') string
       read(string,*) iat, jat, hop3d(1:3), tijup, tijdn, U

       doeshop = hoptowho(iat, hop3d, jat, lattice) .ne. iat
       doeshop = doeshop .and. (abs(tijup).gt.1.d-6 .or. abs(tijdn).gt.1.d-6)

       do it = 0, lattice%ncell-1

          is    = iat + it*natom
          js    = hoptowho(is, hop3d, jat, lattice)

          if (doeshop) then
             ihop = ihop + 1
             hamilt%tckb(1,ihop)  = min(is,js)
             hamilt%tckb(2,ihop)  = max(is,js)
             hamilt%plink(1,ihop) = is
             hamilt%plink(2,ihop) = js
             if (it .eq. 0) then
                jhop = jhop + 1
                hamilt%tlink(1,jhop)   = tijup
                hamilt%tlink(2,jhop)   = tijdn
                hamilt%tlink(3:5,jhop) = hop3d(1:3)
             endif
          endif

          twisthop = sum((pos(:,js) - pos(:,is) - hop3d(:))*ktwist(:))
          
          hamilt%hopup(is,js) = hamilt%hopup(is,js) + tijup * exp( im*twisthop)
          hamilt%hopdn(is,js) = hamilt%hopdn(is,js) + tijdn * exp( im*twisthop)
          hamilt%Uv(is,js)    = hamilt%Uv(is,js)    + U

          if (is .ne. js) then
             hamilt%hopup(js,is) = hamilt%hopup(js,is) + tijup * exp(-im*twisthop)
             hamilt%hopdn(js,is) = hamilt%hopdn(js,is) + tijdn * exp(-im*twisthop)
             hamilt%Uv(js,is)    = hamilt%Uv(js,is)    + U
          endif

       enddo

    enddo

    !Define phase on each atom compatibly with BC
    do iat = 0, nsites - 1
       jat = mod(iat,natom)
       hamilt%phase(iat) = exp(im*sum((kpoint(:)-ktwist(:))*(pos(:,iat)-pos(:,jat))))
       hamilt%phase(iat) = hamilt%phase(iat) * exp(im*(-sum(ktwist(:)*pos(:,iat))))
    enddo

    hamilt%constructed = .true.

    call find_neighbors(hamilt)
    hamilt%neig_found = .true.

    !if (sum(hamilt%tnneig) .ne. 2*nhop) then
    !   write(*,*) "neighbors and hops are incompatible. Bug"
    !   write(*,*) "nhop =", nhop
    !   write(*,*) "nneig =", sum(hamilt%tnneig)
    !   stop
    !endif
    call dqmc_hamilt_groupckb(hamilt)

 end subroutine 

!-------------------------------------------------------------------!

 subroutine find_neighbors(hamilt)
   !
   ! Given the matrices Uv,Jv and hop it finds which atoms 
   ! are neighbors and classify neighbors as with respect 
   ! to interaction or hopping. Returns tnneig(is): number
   ! of neighbors of is; tneig(is,j): j-th neighbor of "is".
   ! Analogous definition for Unneig,Uneig and Jnneig,Jneig.
   !
     type(Hamiltonian_t),intent(inout) :: hamilt

     integer is, js, ntsites, nUsites, nJsites, n
     integer, pointer :: tnneig(:),Unneig(:),Jnneig(:)
     integer, pointer :: tneig(:,:),Uneig(:,:),Jneig(:,:)
     integer, pointer :: tsite(:),Usite(:),Jsite(:)
     
     if(.not.hamilt%constructed) stop'Hamiltonian needs to be &
        &constructed before neig can be found'
    
     n=size(hamilt%hopup,1)

     nullify(tnneig, Unneig, Jnneig)
     nullify( tneig,  Uneig,  Jneig)
     nullify( tsite,  Usite,  Jsite)

     allocate(tnneig(0:n-1),Unneig(0:n-1),Jnneig(0:n-1))
     allocate(tneig(0:n-1,n),Uneig(0:n-1,n),Jneig(0:n-1,n))
     allocate(tsite(n),Usite(n),Jsite(n))
    
     !Count the number of sites connected to "is" by t, U and J and store site label
     tnneig = 0 
     Unneig = 0 
     Jnneig = 0
     do is = 0, n-1
        do js = 0, n-1
           if (abs(hamilt%hopup(is,js)).gt.1d-9 .and. is/=js) then
              tnneig(is) = tnneig(is) + 1
              tneig(is,tnneig(is)) = js
           elseif (abs(hamilt%hopdn(is,js)).gt.1d-9 .and. is/=js) then
              tnneig(is) = tnneig(is) + 1
              tneig(is,tnneig(is)) = js
           endif  
           if (abs(hamilt%Uv(is,js)).gt.1d-9) then
              Unneig(is) = Unneig(is) + 1
              Uneig(is,Unneig(is)) = js
           endif
           if (abs(hamilt%Jv(is,js)).gt.1d-9) then
              Jnneig(is) = Jnneig(is) + 1
              Jneig(is,Jnneig(is)) = js
           endif
        enddo
     enddo
    
     !Count the number of sites effectively having t,U or J on them
     ntsites = 0 
     nUsites = 0 
     nJsites = 0
     do is = 0, n-1
        if (tnneig(is).gt.0) then
           ntsites = ntsites + 1
           tsite(ntsites) = is
        endif 
        if(Unneig(is).gt.0)then
           nUsites = nUsites + 1
           Usite(nUsites) = is
        endif
        if(Jnneig(is).gt.0)then
           nJsites = nJsites + 1
           Jsite(nJsites) = is
        endif
     enddo
    
     !Fill hamiltonian variables
     hamilt%maxneig = max(maxval(tnneig), maxval(Unneig), maxval(Jnneig))
    
     !number of neighbors of each site
     hamilt%tnneig => tnneig
     hamilt%Unneig => Unneig
     hamilt%Jnneig => Jnneig
    
     !neighbors of each site
     hamilt%tneig => tneig
     hamilt%Uneig => Uneig
     hamilt%Jneig => Jneig
    
     !number of sites on which either t or U or J is different from 0
     hamilt%ntsites = ntsites
     hamilt%nUsites = nUsites
     hamilt%nJsites = nJsites
    
     !sites on which either t or U or J is different from 0
     hamilt%tsite => tsite
     hamilt%Usite => Usite
     hamilt%Jsite => Jsite
    
 end subroutine find_neighbors

!-------------------------------------------------------------------!

 subroutine count_hop_class(lattice, hamilt)
   !
   ! Construct a list of only hopping classes. These are a 
   ! subset of all the distance classes for which the hopping 
   ! is non zero. Find which hoppings are equivalent by 
   ! symmetry. nhopsite(iat,iclass) returns how many hoppings 
   ! in class iclass iat has. ordered_neig(iat,ineig) returns 
   ! the ineig-th neighbor of iat. Neighbors are here ordered 
   ! according to classes. So if nhopsite(2,1)=3 the first 
   ! three neigbors of 2 listed in ordered_neig are those
   ! belonging to class 1. tneig has a similar content but 
   ! neighbors are not ordered.
   !
    type(Hamiltonian_t)  :: hamilt
    type(lattice_t)      :: lattice

    integer              :: isite, ineig, jsite, newclass, ihc
    integer              :: maxtneig, nsites 
    integer              :: hopclass(lattice%nclass),nhopclass
    integer, allocatable :: pairhopclass(:,:)
    complex*16           :: tvaluetmpup(lattice%nclass)
    complex*16           :: tvaluetmpdn(lattice%nclass)
    
    nsites      = size(hamilt%tnneig)
    nhopclass   = 0
    hopclass(:) = 0
    
    allocate(pairhopclass(0:nsites-1,nsites))
    !class for the pair (isite, t-neighbour of site). 
    pairhopclass(:,:) = -1
    
    !Select classes for which hopping is different from 0
    !Uses symmetry info found in lattice.
    do isite = 0, nsites-1
       do ineig = 1, hamilt%tnneig(isite)
    
          jsite = hamilt%tneig(isite,ineig) 
          !(isite,jsite) is a pair of sites with non-zero hopping
          if(isite == jsite) cycle
          newclass = lattice%myclass(isite,jsite)
    
          !First check that this class has not been already found
          do ihc = 1, nhopclass
             if(newclass == hopclass(ihc))exit
          enddo
    
          !If hopclass is new, increase number of classes
          if(ihc == nhopclass+1)then 
             hopclass(ihc)  = newclass
             nhopclass      = ihc
             tvaluetmpup(ihc) = hamilt%hopup(isite,jsite)
             tvaluetmpdn(ihc) = hamilt%hopdn(isite,jsite)
          endif
    
          !Assign pair to the class
          pairhopclass(isite,ineig) = ihc
       enddo
    enddo
    hamilt%nhopclass = nhopclass
    maxtneig = maxval(hamilt%tnneig(:))
    
    allocate(hamilt%myhopclass(0:nsites-1,maxtneig))
    allocate(hamilt%tupvalue(nhopclass))
    allocate(hamilt%tdnvalue(nhopclass))

    !Assign t value to each class and class number to each t-pair
    hamilt%tupvalue(1:nhopclass) = tvaluetmpup(1:nhopclass)
    hamilt%tdnvalue(1:nhopclass) = tvaluetmpdn(1:nhopclass)

    do ineig = 1, maxtneig
       hamilt%myhopclass(0:nsites-1,ineig) = pairhopclass(0:nsites-1,ineig)
    enddo

    deallocate(pairhopclass)

 end subroutine count_hop_class

!-------------------------------------------------------------------!

 subroutine count_local_classes(lattice,hamilt)
  !
  ! Construct a list of only local classes i.e. the subset 
  ! of the distance classes defined on the same site. Note
  ! that we use the symmetry instead of the value to group
  ! sites together. This is because local classes refers to 
  ! both U and mu and two sites may have same value of U 
  ! but different mu. If equivalent by symmetry, the two 
  ! sites have, however, necessarily same u and mu.
  !
    type(Hamiltonian_t), intent(inout) :: hamilt
    type(lattice_t), intent(in)        :: lattice

    integer :: isite, iclass, jclass, nsites, nlocclass
    integer :: localtmp(lattice%nclass)
    real*8  :: Utmp(lattice%nclass) 
    real*8  :: muup(lattice%nclass)
    real*8  :: mudn(lattice%nclass)

    nlocclass = 0
    nsites    = size(hamilt%tnneig)
    allocate(hamilt%mylocclass(0:nsites-1))

    !Loop over all sites
    do isite = 0, nsites-1

       iclass = lattice%myclass(isite,isite)
    
       !Check whether class was already found
       do jclass = 1, nlocclass
          if(iclass.eq.localtmp(jclass))exit
       enddo
    
       !if not augment the number of local classes
       if(jclass>nlocclass)then
          nlocclass = nlocclass+1
          localtmp(jclass) = iclass
          !Save value of mu and U for this class
          Utmp(jclass) = hamilt%Uv(isite,isite)
          muup(jclass) = -dble(hamilt%hopup(isite,isite))+hamilt%mu_up
          mudn(jclass) = -dble(hamilt%hopdn(isite,isite))+hamilt%mu_dn
       endif
    
       !assign site to a class
       hamilt%mylocclass(isite) = jclass
    enddo
    
    !Store the number of different classes
    hamilt%nlocclass = nlocclass
    
    !Store the value of U and on-site energy (shifted by mu) for each class
    allocate(hamilt%Uvalue(nlocclass))
    allocate(hamilt%muupvalue(nlocclass))
    allocate(hamilt%mudnvalue(nlocclass))
    hamilt%Uvalue(1:nlocclass)  = Utmp(1:nlocclass)
    hamilt%muupvalue(1:nlocclass) = muup(1:nlocclass)
    hamilt%mudnvalue(1:nlocclass) = mudn(1:nlocclass)

 end subroutine

!-------------------------------------------------------------------!

 subroutine dqmc_hamilt_groupckb(hamilt)

    type(hamiltonian_t), intent(inout) :: hamilt

    integer :: nh, ickb, ih, jh, kh, ns, hop(3)
    logical :: skip 
    logical, pointer :: vs(:)

    ! Allocate array for visited sites (vs)
    ns = 1 + maxval(hamilt%tckb)
    allocate(vs(0:ns-1))

    ! Count the number of hoopings
    nh = size(hamilt%tckb,2)
    ickb = 0
    do ih = 0, nh-1
       ! Check whether hopping has been already assigned
       if (hamilt%tckb(3,ih) .gt. 0) cycle
       ! "ih" is the first of group ickb
       ickb = ickb + 1
       ! mark all sites as not visited ...
       vs(0:ns-1) = .false.
       ! ... except for those of "ih"
       vs(hamilt%tckb(1:2,ih)) = .true.
       ! save ckb group
       hamilt%tckb(3,ih) = ickb
       ! Start finding other hopping
       do jh = 0, nh-1
          ! If hopping already include or...
          skip = hamilt%tckb(3,jh) .gt. 0
          ! if one of the site was already visited then...
          skip = skip .or. vs(hamilt%tckb(1,jh))
          skip = skip .or. vs(hamilt%tckb(2,jh))
          ! do not include "jh" in ickb
          if (skip) cycle
          ! otherwise do it and...
          hamilt%tckb(3,jh) = ickb
          ! mark the sites as visited
          vs(hamilt%tckb(1:2,jh)) = .true.
       enddo
    enddo

    deallocate(vs)
   
    ! Sort them so that hops in same group are neighbors
    do ih = 1, nh-1
       hop(1:3) = hamilt%tckb(1:3,ih)
       do jh = 0, ih-1
          if (hop(3) < hamilt%tckb(3,jh)) exit
       enddo
       do kh = ih-1, jh, -1
          hamilt%tckb(1:3,kh+1) = hamilt%tckb(1:3,kh)
       enddo
       hamilt%tckb(1:3,jh) = hop(1:3)
    enddo

    !do ih = 0, nh-1
    !   write(*,*) hamilt%tckb(1:3, ih)
    !enddo

 end subroutine dqmc_hamilt_groupckb

!-------------------------------------------------------------------!

 subroutine group_hopping(hamilt, n, nt, tmap, tupvalue, tdnvalue)
  !
  !  Assign to every non-zero hopping a class based 
  !  on its value. Equal hopping matrix elements belong 
  !  to the same class : tmap(i,j)=it.
  !  tvalue(it) contains the value for a given class.
  !  Exclude local site energies.
  !
    type(Hamiltonian_t), intent(in) :: hamilt
    integer, intent(in)             :: n
    integer, intent(out)            :: tmap(n,n), nt
    real*8, pointer, intent(inout)  :: tupvalue(:)
    real*8, pointer, intent(inout)  :: tdnvalue(:)

    integer         :: is, js, it, jn
    real*8          :: tup, tdn

    !Initialize
    tmap = 0
    nt   = 0

    if (associated(tupvalue)) deallocate(tupvalue)
    if (associated(tdnvalue)) deallocate(tdnvalue)

    if(maxval(hamilt%tnneig) > 0) then
       allocate(tupvalue(hamilt%ntsites*hamilt%maxneig))
       allocate(tdnvalue(hamilt%ntsites*hamilt%maxneig))
       !Assign same class to identical matrix elements
       do is = 0, n-1
          do jn = 1, hamilt%tnneig(is)
             js = hamilt%tneig(is,jn)
             if (is .eq. js) cycle
             tup  = dble(hamilt%hopup(is,js))
             tdn  = dble(hamilt%hopdn(is,js))
             do it = 1, nt
                if (abs(tup-tupvalue(it)) < 1.d-6 .and. &
                 &  abs(tdn-tdnvalue(it)) < 1.d-6) exit
             enddo
             if (it .eq. nt+1) then
               nt = it
               tupvalue(it) = tup
               tdnvalue(it) = tdn
             endif
             tmap(is+1, js+1) = it
          enddo
       enddo
    else
       nt = 1
       allocate(tupvalue(1))
       allocate(tdnvalue(1))
       tupvalue = 0.d0
       tdnvalue = 0.d0
       do is = 1, n 
          tmap(is,is) = 1
       enddo
    endif

 end subroutine group_hopping

!-------------------------------------------------------------------!

end module DQMC_HAMILT
