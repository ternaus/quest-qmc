module DQMC_RECLATT

use DQMC_GEOM_PARAM
use DQMC_LATT
use DQMC_CFG

implicit none

type :: recip_lattice_t

 integer                :: ndim
 integer                :: nkpts              !number of k-points (equal to ncell)
 real*8, pointer        :: klist(:,:)         !list of k-points(nkpts,rdim)
 real*8                 :: kpoint(rdim)       !Input k-point 
 real*8                 :: ktwist(rdim)       !Twist vector
 real*8                 :: kcs(rdim,rdim)     !cartesian component of reciprocal superlattice***
 real*8                 :: ks(rdim,rdim)      !fractional components of reciprocal superlattice*** 
 real*8                 :: kc(rdim,rdim)      !cartesian components of reciprocal unit cell***
                                              !*** rows of these matrices are the vectors

 integer                :: nmomenta           !number of momenta for pair of particles
 real*8, pointer        :: ksum(:,:)          !list of momentum of pair of particles(nmomenta,rdim)
 integer, pointer       :: kmate(:,:)         !mate of each k-points that gives total momentum ksum
                                              !(nkpts, nmomenta)

 integer                :: nclass_k           !number of inequivalent k-points
 integer, pointer       :: myclass_k(:)       !class for each k-point (nkpts)
 integer, pointer       :: class_size_k(:)    !number of equivalent k-points in each class (nclass_k)
 integer, pointer       :: class_repr_k(:)    !representative k-point for each class (nclass_k)
   
 complex*16, pointer    :: FourierC(:,:) 

 logical                :: initialized
 logical                :: constructed
 logical                :: analyzed
end type recip_lattice_t


contains

!----------------------------------------------------------------------!

 subroutine free_reclatt(reclatt)

 type(recip_lattice_t), intent(inout) :: reclatt

 if(associated(reclatt%klist))         deallocate(reclatt%klist)
 if(associated(reclatt%FourierC))      deallocate(reclatt%FourierC)
 if(associated(reclatt%myclass_k))     deallocate(reclatt%myclass_k)
 if(associated(reclatt%class_size_k))  deallocate(reclatt%class_size_k)
 if(associated(reclatt%class_repr_k))  deallocate(reclatt%class_repr_k)
  
 end subroutine free_reclatt


!---------------------------------------------------------------------
! Read and fill most of the variables that define the lattices in
! real and reciprocal space.
!---------------------------------------------------------------------
subroutine init_recip_latt(lattice,recip_lattice,applytwist,cfg) 
 integer                                     :: ndim, i, j
 real*8                                      :: projk(rdim)
 real*8, pointer                             :: kc(:,:) 
 real*8, pointer                             :: kcs(:,:) 
 real*8, pointer                             :: ktwist(:) 
 real*8, pointer                             :: kpoint(:) 
 real*8, pointer                             :: ac(:,:) 
 real*8, pointer                             :: scc(:,:) 
 type(lattice_t),intent(in),target           :: lattice
 type(recip_lattice_t),intent(out),target    :: recip_lattice
 type(config),intent(in),target              :: cfg
 logical, intent(in)                         :: applytwist

 if(.not.lattice%initialized)stop'Need to initialize lattice before recip_lattice'

 !alias arrays
 scc=>lattice%scc
 ac=>lattice%ac
 kc=>recip_lattice%kc
 kcs=>recip_lattice%kcs
 ktwist=>recip_lattice%ktwist

 nullify(kpoint)
 ndim=lattice%ndim
 
 !twist of the boundary condition in untis of pi
 recip_lattice%kpoint=0.d0
 if(applytwist)then
    call cfg_get(cfg,'bcond',i,kpoint)
    if(i<ndim)stop'Number of boundary condition must at least equal ndim'
    recip_lattice%kpoint(1:ndim)=kpoint(1:ndim)
    !Since QUEST cannot handle complex hoppings we stop if the k-point is /=0.
    if(sum((kpoint(1:ndim)-nint(kpoint(1:ndim)))**2)>toll)then
     write(*,*)'Only integer twists are possible'
     stop
    endif
    deallocate(kpoint)
 endif
 kpoint=>recip_lattice%kpoint
  
 !Take care of k-space cells
 !Rows of kc are reciprocal lattice basis vectors in cartesian coordinates
 call get_inverse(ac,kc)
 kc(:,:)=kc*2.d0*pi

 !find reciprocal lattice vectors of the supercell
 call get_inverse(scc,kcs)
 kcs(:,:)=kcs(:,:)*2.d0*pi

 !convert twist in k-vector in cartesion coordinates
 do j=1,rdim 
  ktwist(j)=0.5*sum(kpoint(:)*kcs(j,:))
 enddo
 kpoint(:)=ktwist(:)

 !and project it inside the reciprocal supercell
 do i=1,rdim
  projk(i)=sum(kpoint(:)*scc(:,i))/(2.0*pi)
 enddo
 projk(:)=projk(:)-nint(projk(:))
 !then find its cartesian coordinates
 do i=1,rdim
  ktwist(i)=sum(projk(:)*kcs(:,i))
 enddo

 !find components of kcs in units of kc
 do i=1,3 
  do j=1,3
   recip_lattice%ks(i,j)=sum(kcs(i,:)*ac(:,j))/(2.d0*pi)
  enddo 
 enddo

 recip_lattice%ndim=ndim
 !number of k-points is equal to number of cell in supercell
 recip_lattice%nkpts=lattice%ncell

 recip_lattice%initialized=.true.

 !write to stdout
 ! if(applytwist)then
 !    write(*,*)
 !    write(*,*)' Reciprocal lattice basis vectors'
 !    write(*,'(3f14.7)')((kc(i,j),j=1,rdim),i=1,rdim)
 !    write(*,*)
 !    write(*,*)' Reciprocal super-lattice basis vectors'
 !    write(*,'(3f14.7)')((kcs(i,j),j=1,rdim),i=1,rdim)
 !    write(*,*)
 !    write(*,*)'Twist vector'
 !    write(*,'(3f14.7)')(ktwist(i),i=1,rdim)
 !    write(*,*)'================================================================'
 ! endif

end subroutine init_recip_latt




!----------------------------------------------------------------------------
! construct the list of k-vectors inside the 1st BZ (klist).
!----------------------------------------------------------------------------
subroutine construct_recip_lattice(recip_lattice)
 integer                :: nkx,nky,nkz,ndim,ikx,iky,ikz,ikv(rdim),ncubex,ncubey,ncubez, &
                           & ilist,nedge,nfound,nkpts,i,j
 integer*8, allocatable :: indedge(:)
 real*8                 :: invkc(rdim,rdim)
 real*8, allocatable    :: kset(:,:)
 real*8, pointer        :: klist(:,:) 
 type(recip_lattice_t), intent(inout)  :: recip_lattice

 if(.not.recip_lattice%initialized)stop'Need to initialize recip_lattice before construction'

 call get_inverse(recip_lattice%kc,invkc)

 ndim=recip_lattice%ndim
 allocate(klist(recip_lattice%nkpts,rdim))

 nkx=0; if(ndim>0)nkx=1 
 nky=0; if(ndim>1)nky=1 
 nkz=0; if(ndim>2)nkz=1
 !Define kset as the set of points surrouding the origin.
 !26 in 3D, 8 in 2D, 2 in 2D. This helps definining the
 !1st brillouin zone.
 allocate(kset(3**ndim-1,rdim),indedge(recip_lattice%nkpts))
 i=0
 do ikx=-nkx,nkx; ikv(1)=ikx
  do iky=-nky,nky; ikv(2)=iky
   do ikz=-nkz,nkz; ikv(3)=ikz
    if(sum(ikv**2)>0)then
     i=i+1 
     do j=1,rdim
      kset(i,j)=sum(ikv(:)*recip_lattice%kc(:,j))
     enddo
    endif
   enddo 
  enddo 
 enddo

 ncubex=0; ncubey=0; ncubez=0
 !Start by including the ktwist vector
 ilist=1; klist(ilist,:)=recip_lattice%ktwist(:); nedge=0
 if(ndim>0)then
  do 
   !Look at increasingly large "cubes" around the k-space origin.
   ncubex=ncubex+1
   if(ndim>1)ncubey=ncubey+1
   if(ndim>2)ncubez=ncubez+1
   nfound=0

   do ikz=-ncubez,ncubez
    ikv(3)=ikz

    !constant y edges : kz constant equal to ikz, ky constant equal to -ncubey
    !or ncubey, and ikx varying between +/-ncubex 
    if(ncubey/=0)then
     do iky=-ncubey,ncubey,2*ncubey
      ikv(2)=iky
      do ikx=-ncubex,ncubex
       ikv(1)=ikx
       !Include k-point associated to ikv(:) in klist if it belongs to 1st BZ
       call check_and_update
      enddo
     enddo
    endif

    !constant x edges : kz constant equal to ikz, kx constant equal to -ncubex
    !or ncubex, and ky varying in +/-(ncubey-1) (special care when ncubey=0 i.e. 1D)
    do ikx=-ncubex,ncubex,2*ncubex
     ikv(1)=ikx
     do iky=min(0,-ncubey+1),max(0,ncubey-1)
      ikv(2)=iky
      !Include k-point associated to ikv(:) in klist if it belongs to 1st BZ
      call check_and_update
     enddo
    enddo

   enddo

   !constant kz faces i.e. kz=-ncubez,ncubez
   if(ncubez>0)then
    do ikz=-ncubez,ncubez,2*ncubez
     ikv(3)=ikz
     do iky=-ncubey+1,ncubey-1
      ikv(2)=iky
      do ikx=-ncubex+1,ncubex-1
       ikv(1)=ikx
       !Include k-point associated to ikv(:) in klist if it belongs to 1st BZ
       call check_and_update
      enddo
     enddo
    enddo
   endif

   !no new k-points were found : exit
   if(nfound==0)exit
  enddo
 endif
 nkpts=ilist

 !Store k-vectors in asceding order of length
 call vsort(klist,nkpts)
 recip_lattice%klist=>klist

 
 ! write(*,*)'Reciprocal lattice (1st Brillouin zone)'
 ! write(*,'(/,A,1x,i4)')' Number of k-points found: ',nkpts
 ! if(nkpts/=recip_lattice%nkpts)then
 !  write(*,*)'This is different from', recip_lattice%nkpts
 !  stop 
 ! endif
 ! write(*,'(3f12.6)')((klist(i,j),j=1,rdim),i=1,recip_lattice%nkpts)
 ! write(*,*)'================================================================'

 recip_lattice%constructed=.true.

 contains

 subroutine  check_and_update
  integer i,ie
  real*8 :: kpt(rdim),projk(rdim),diff(rdim)
  logical on_edge,included

  !Get k-point associated to ikv(:)
  do i=1,3
    kpt(i)=sum(ikv(:)*recip_lattice%kcs(:,i))+recip_lattice%ktwist(i)
  enddo

  if(closer_to_zero(kpt,kset,on_edge,ndim))then
   !kpt is either inside or on edge of BZ
   included=.false.
   if(on_edge)then
    !Need to see whether equivalent was already included
    do ie=1,nedge
     !compute the distance with other points on edge.
     diff(:)=klist(indedge(ie),:)-kpt(:)
     !convert it in units of primitive reciprocal lattice
     do i=1,3
       projk(i)=sum(diff(:)*invkc(:,i))
     enddo
     !If all integers than point was already included
     projk(:)=projk(:)-nint(projk(:))
     if(sum(projk**2)<1.d-6)then
      included=.true.
     endif
    enddo
   endif

   if(.not.included)then
    !Include the point in klist 
    nfound=nfound+1
    ilist=ilist+1
    klist(ilist,:)=kpt(:)
    !save info if k-point was on edge
    if(on_edge)then
     nedge=nedge+1
     indedge(nedge)=ilist
    endif 
   endif
  endif
 end subroutine check_and_update

 subroutine vsort(vvec,n)
 integer n,i,j,k
 real*8 vvec(n,rdim),vlen(n),av(rdim),a
 do i=1,n
  vlen(i)=sum(vvec(i,:)**2)
 enddo
 do i=2,n
  a=vlen(i)
  av(:)=vvec(i,:)
  do j=1,i-1
   if(a+1.d-10<vlen(j))exit
  enddo
  do k=i-1,j,-1
   vlen(k+1)=vlen(k)
   vvec(k+1,:)=vvec(k,:)
  enddo
  vlen(j)=a
  vvec(j,:)=av(:)
 enddo
 end subroutine

end subroutine




!----------------------------------------------------------
!Determine whether a k-point (kpt) is closer to k=0 than
!to any of the points in kset. If half way set to true
!but also set on_edge to true.
!----------------------------------------------------------
logical function closer_to_zero(ktp,kset,on_edge,ndim)
integer, intent(in) :: ndim
real*8, intent(in)  :: ktp(rdim),kset(3**ndim-1,rdim)
real*8              :: diff(rdim),dist0,dist
logical             :: on_edge
integer             :: i
on_edge=.false.
diff(:)=ktp(:)**2
!Distance from origin
dist0=sum(diff)
!Loop over the points surrounding the origin
do i=1,3**ndim-1
 diff(:)=(ktp(:)-kset(i,:))**2
 !distance from kset(i,:)
 dist=sum(diff)
 if(dist+toll<dist0)then
  !ktp is closer to kset(i,:) than to k=0.0 : outside 1st BZ
  closer_to_zero=.false.
  return
 else
  !ktp is rougly equidistant from kset(i,:) and k=0.0 : on edge of 1st BZ
  if(dist-toll<dist0)on_edge=.true.
 endif
enddo 
closer_to_zero=.true.
end function closer_to_zero

  !--------------------------------------------------------------------!

  subroutine DQMC_init_kmate(reclatt, nmom, ksum)

  type(recip_lattice_t), intent(inout) :: reclatt
  integer, intent(in)  :: nmom
  real(wp), intent(in) :: ksum(nmom, rdim)

  integer :: ik, jk, i, im
  real(wp), dimension(rdim) :: kpt, kdiff, kproj
  real(wp) :: kcinv(rdim, rdim)

  reclatt%nmomenta = nmom
  allocate(reclatt%kmate(reclatt%nkpts, nmom))
  allocate(reclatt%ksum(nmom, rdim))

  reclatt%ksum = ksum

  call get_inverse(reclatt%kc, kcinv)
  !Loop over all k-points and find their partner
  do im = 1, reclatt%nmomenta
     do ik = 1, reclatt%nkpts
        !Compute k-point so that sum equals totk.
        kpt(:) = reclatt%ksum(im,:) - reclatt%klist(ik,:)
        !Find k-point in 1st BZ that corresponds to kpt
        do jk = 1, reclatt%nkpts
           !Compute difference in units of Primitive cell k-vectors
           kdiff(:) = kpt(:) - reclatt%klist(jk,:)
           do i = 1, rdim
              kproj(i) = sum(kdiff(:) * kcinv(i,:))
           enddo
           !If everything is an integer, we found the vector
           if(sum((kproj-nint(kproj))**2) < 1.d-6)then
              reclatt%kmate(ik,im) = jk
              exit
           endif
        enddo
     enddo
  enddo

  !Check
  do im = 1, reclatt%nmomenta
     do ik = 1, reclatt%nkpts
        jk=reclatt%kmate(ik,im)
        if(reclatt%kmate(jk,im) /= ik) then
           write(*,*)'K-mates do not match. Stop.'
           stop
        endif
     enddo
  enddo

  do im = 1, reclatt%nmomenta
     write(*,*)'Total momentum :', (ksum(im,i), i =1 ,rdim)
     write(*,*)'Pairs :'
     do ik = 1, reclatt%nkpts
        jk = reclatt%kmate(ik,im)
        write(*,'(3(i3, 3f10.5,10x))') ik, (reclatt%klist(ik,i), i=1,rdim),&
                   jk, (reclatt%klist(jk,i), i=1,rdim)
     enddo
     write(*,*)'==============================================================='
  enddo

  end subroutine 

  !--------------------------------------------------------------------!

  subroutine DQMC_Fill_FourierC(Reciplattice,lattice)
  ! Fill the matrix of Fourier coefficients
  type(lattice_t),intent(in),target           :: lattice
  type(recip_lattice_t),intent(inout),target    :: Reciplattice
  integer               :: nt,nk,i,ii,j
  real*8, pointer       :: tr(:,:)  
  real*8, pointer       :: kpts(:,:)
  integer, pointer      :: indx(:)  

  !initialize
  tr   => lattice%translation
  kpts => RecipLattice%klist
  indx => RecipLattice%class_repr_k
  nt   = lattice%ncell
  nk   = RecipLattice%nclass_k
  allocate(Reciplattice%FourierC(nt,nk))

  !compute
  do i = 0, nt-1
     ii = i + 1
     do j = 1, nk
        Reciplattice%FourierC(ii,j) = exp(im*sum(tr(:,i)*kpts(indx(j),:)))
     enddo
  enddo
  Reciplattice%FourierC = Reciplattice%FourierC / sqrt(dble(nt))

  end subroutine DQMC_Fill_FourierC


end module DQMC_RECLATT

