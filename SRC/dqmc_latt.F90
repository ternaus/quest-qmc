module DQMC_LATT

use DQMC_GEOM_PARAM
use DQMC_UTIL

implicit none

type :: lattice_t
 integer                :: nsites            !number of total sites 
 integer                :: natom             !number of sites inside the primitive cell
 integer                :: ncell             !number of cells inside supercell
 integer                :: ndim              !number of extended dimensions
 integer                :: sc(rdim,rdim)     !fractionary components of supercell**

 real*8                 :: ac(rdim,rdim)     !cartesian components of primitive cell**
 real*8                 :: scc(rdim,rdim)    !cartesian components of supercell** 
 real*8, pointer        :: pos(:,:)          !fractionary position of each site (rdim,0,nsites-1)
 real*8, pointer        :: cartpos(:,:)      !cartesian coordinates of each site (rdim,0,nsites-1)
 real*8, pointer        :: xat(:,:)          !fractional coordinate of site inside primitive cell(rdim,0:natom-1)
 real*8, pointer        :: phase(:)          !Phase for order parameter(0:nsites-1)
 real*8, pointer        :: translation(:,:)  !list of translation vectors(columns)(rdim,0:ncell-1)
                                                       !** columns of these matrices are the vectors

 integer                :: nclass            !number of classes for distance 
 integer, pointer       :: myclass(:,:)      !class for pair of sites (0:nsites-1, 0:nsites-1)
 integer, pointer       :: class_size(:)     !number of equivalent pairs in each class (nclass)
 real*8, pointer        :: class_label(:,:)  !label for pair classes
 integer, pointer       :: gf_phase(:,:) 

 character*3, pointer   :: olabel(:)         !label of each site in primitive cell

 logical                :: initialized
 logical                :: constructed
 logical                :: analyzed
end type lattice_t


contains

!--------------------------------------------------------------------!

subroutine free_lattice(latt)

type(lattice_t), intent(inout) :: latt

 if(associated(latt%pos))            deallocate(latt%pos)
 if(associated(latt%xat))            deallocate(latt%xat)
 if(associated(latt%phase))          deallocate(latt%phase)
 if(associated(latt%olabel))         deallocate(latt%olabel)
 if(associated(latt%cartpos))        deallocate(latt%cartpos)
 if(associated(latt%myclass))        deallocate(latt%myclass)
 if(associated(latt%gf_phase))       deallocate(latt%gf_phase)
 if(associated(latt%class_size))     deallocate(latt%class_size)
 if(associated(latt%translation))    deallocate(latt%translation)
 if(associated(latt%class_label))    deallocate(latt%class_label)

end subroutine free_lattice


!---------------------------------------------------------------------
! Read and fill most of the variables that define the lattices in
! real and reciprocal space.
!---------------------------------------------------------------------
subroutine init_lattice(lattice, SOP)
 integer, intent(in) :: SOP 
 integer          :: ndim,nsites,natom,ncell,ios,i,j
 real*8           :: ainv(rdim,rdim)
 real*8, pointer  :: ac(:,:) 
 real*8, pointer  :: scc(:,:) 
 integer, pointer :: sc(:,:) 
 character*3      :: olab1
 character*50     :: string
 logical          :: ldum
 type(lattice_t),target    :: lattice

 !alias arrays
 sc=>lattice%sc
 scc=>lattice%scc
 ac=>lattice%ac
 
 !Read number of dimensions
 ldum=move_to_record(INPUT_FIELDS(NDIM_F),inpunit)
 read(inpunit,*,iostat=ios)ndim
 if(ios/=0)stop ' Problem reading #NDIM field. Stop.'
  
 !read basis cell vectors(cartesian). Basis vectors are columns of ac.
 ldum=move_to_record(INPUT_FIELDS(PRIM_F),inpunit)
 ac(:,:)=0.d0
 do j=1,ndim
  read(inpunit,*,iostat=ios)(ac(i,j),i=1,ndim)
 enddo
 do j=ndim+1,rdim
   ac(j,j)=1.d3
 enddo
 !read(inpunit,*,iostat=ios)((ac(i,j),i=1,rdim),j=1,rdim)
 if(ios/=0)stop ' Problem reading #PRIM field. Stop.'
  
 !read supercell vectors in unit of the basis ones. Vectors are columns of sc.
 ldum=move_to_record(INPUT_FIELDS(SUPER_F),inpunit)
 sc(:,:)=0
 if(ndim>0)then
  do j=1,ndim
    read(inpunit,*,iostat=ios)(sc(i,j),i=1,ndim)
  enddo
  if(ios/=0)stop ' Problem reading #SUPER field. Stop.'
  do i=ndim+1,rdim 
   sc(i,i)=1 
  enddo
 endif

 !compute number of primitive cell inside the supercell
 !Use scc as temporary real*8 array
 scc(:,:)=dble(sc(:,:))
 ncell=nint(abs(get_det(scc)))
  
 !find cartesian component of the supercell
 do i=1,rdim 
  do j=1,rdim
   scc(i,j)=sum(ac(i,:)*sc(:,j))
  enddo
 enddo
  
 !cartesian coordinates for each orbital 
 natom=count_atom()
 nsites=ncell*natom
 ldum=move_to_record(INPUT_FIELDS(ORB_F),inpunit)
 allocate(lattice%olabel(0:natom-1),lattice%xat(rdim,0:natom-1))
 do i=0,natom-1
  read(inpunit,'(A)')string
  read(string,*,iostat=ios)olab1,lattice%xat(1:rdim,i)
  if(ios/=0)stop ' Problem reading #ORB field. Stop.'
  lattice%olabel(i)=olab1
 enddo
 
 !load values on lattice
 lattice%natom=natom
 lattice%nsites=nsites
 lattice%ndim=ndim
 lattice%ncell=ncell

 !convert cartesian atomic coordinates
 call get_inverse(ac,ainv)
 do i=0,natom-1
  call convert_to_fractional(lattice%xat(1:rdim,i),ainv)
 enddo

 lattice%initialized=.true.
 
 !write to stdout
 write(SOP,*)'================================================================'
 write(SOP,*)'Basic real space geometry info'
 write(SOP,*)
 write(SOP,*)'Crystal atomic basis'
 write(SOP,'(i3,3f14.7)')(j, lattice%xat(1:rdim,j),j=0,natom-1)
 write(SOP,*)
 write(SOP,*)'Basis cell vectors'
 write(SOP,'(3f14.7)')((ac(i,j),i=1,rdim),j=1,rdim)
 write(SOP,*)
 write(SOP,'(/,A)')' Supercell vectors (fractionary unit)'
 write(SOP,'(3i5)')((sc(i,j),i=1,rdim),j=1,rdim)
 write(SOP,*)
 write(SOP,*)'Super-Lattice vectors (cartesian)'
 write(SOP,'(3f14.7)')((scc(i,j),i=1,rdim),j=1,rdim)
 write(SOP,*)
 write(SOP,*)'================================================================'
 
end subroutine init_lattice




!-----------------------------------------------------------------------
! Construct the real space lattice.
!-----------------------------------------------------------------------
subroutine construct_lattice(lattice, SOP)
 type(lattice_t),intent(inout),target         :: lattice
 integer, intent(in)     :: SOP
 integer                 :: natom,nsites,iat,jat,xxmax(rdim),xxmin(rdim),ix,iy,iz,icount,j,jcount,ndim,it
 real*8, pointer         :: ac(:,:)  
 real*8, pointer         :: pos(:,:) 
 real*8, pointer         :: xat(:,:) 
 real*8, pointer         :: cartpos(:,:) 
 real*8, pointer         :: translation(:,:) 
 integer, pointer        :: sc(:,:) 
 real*8                  :: projk(rdim),xxat(rdim),xx(rdim),invscc(rdim,rdim),ainv(rdim,rdim),cmin,cmax
 
 if(.not.lattice%initialized)stop'Need to initialize lattice before construct_lattice'
 
 !Initialize local variables/pointers
 ndim   =  lattice%ndim
 natom  =  lattice%natom
 nsites =  lattice%nsites
 sc     => lattice%sc
 ac     => lattice%ac
 xat    => lattice%xat
 allocate(pos(rdim,0:nsites-1),cartpos(rdim,0:nsites-1),translation(rdim,0:lattice%ncell-1))

 !find a supercell that includes the previous one and which is an "easy" multiple of the primitive
 xxmax(:)=max(sc(:,1), sc(:,2), sc(:,3), sc(:,1)+sc(:,2), sc(:,1)+sc(:,3), sc(:,2)+sc(:,3), sc(:,1)+sc(:,2)+sc(:,3), 0)
 xxmin(:)=min(sc(:,1), sc(:,2), sc(:,3), sc(:,1)+sc(:,2), sc(:,1)+sc(:,3), sc(:,2)+sc(:,3), sc(:,1)+sc(:,2)+sc(:,3), 0)

 !Find lattice points inside the supercell without counting the
 !ones on edge twice
 icount=-1
 call get_inverse(lattice%scc,invscc)
 do iz=xxmin(3),xxmax(3)-1
  do iy=xxmin(2),xxmax(2)-1
   do ix=xxmin(1),xxmax(1)-1
    xx(1)=dble(ix); xx(2)=dble(iy); xx(3)=dble(iz)
    !cartesian coordinates
    do j=1,rdim
     xxat(j)=sum(xx(:)*ac(j,:))
    enddo
    !Project to see if inside supercell
    do j=1,rdim
     projk(j)=sum(xxat(:)*invscc(j,:))
    enddo
    cmin=minval(projk); cmax=maxval(projk)
    !If inside add it to the list of translations
    if(cmin>-toll .and. cmax<1.d0+toll)then
     jcount=icount+1
     !If on edge check whether translation has already  been added
     if(cmin<toll .and. cmax>1.d0-toll)then 
      !Ok, lattice point is on the edge of supercell
      do jcount=0,icount
       !See whether it differs by an already stored translation by a
       !vector belonging to the super-lattice and....
       do j=1,rdim
        projk(j)=sum((xxat(:)-translation(:,jcount))*invscc(j,:))
       enddo
       projk(1:ndim)=projk(1:ndim)-nint(projk(1:ndim))
       !...exit the loop if it does
       if(sum(projk(1:rdim)**2)<toll)exit
      enddo
     endif
     !Found new lattice point
     if(jcount>icount)then
      translation(:,jcount)=xxat(:)
      icount=jcount
     endif
    endif
   enddo
  enddo
 enddo

 !Deal with the  orbitals in such a way that the set [0,1...natom-1]
 !is translated into a set of the form [natom*it,natom*it+1,....,natom*(it+1)-1] 
 !for any arbitrary translation ("it" is an integer running from 0 to ncell-1)
 call get_inverse(ac,ainv)
 do iat=0,natom-1
  do j=1,rdim
   xxat(j)=sum(xat(:,iat)*ac(j,:))
  enddo
  do it=0,lattice%ncell-1
   jat=iat+natom*it
   cartpos(:,jat)=xxat(:)+translation(:,it)
   do j=1,rdim
    pos(j,jat)=sum(cartpos(:,jat)*ainv(j,:))
   enddo
  enddo
 enddo

 !Printing out
  write(SOP,*)'Real space lattice'
  write(SOP,*)
  write(SOP,*)'Number of orbitals in primitive cell: ',natom
  write(SOP,*)'Total number of orbitals:             ',nsites
  write(SOP,*)'index  label   type       X           Y         Z   '
  do iat=0,nsites-1
   icount=mod(iat,natom)
   write(SOP,'(i3,1x,A,1x,i3,3f14.5)')iat,lattice%olabel(icount),icount,(cartpos(j,iat),j=1,3)
  enddo
  write(SOP,*)'================================================================'

 lattice%pos => pos
 lattice%cartpos => cartpos
 lattice%translation => translation

 lattice%constructed=.true.

 !51 format(i4,5x,A3,3x,i4,5(1x,f10.5))
end subroutine



!----------------------------------------------------------
!Count the number of sites in the primitive cell
!----------------------------------------------------------
integer function count_atom()  result (natom)
 character*50  :: string
 character*3   :: olab1
 integer       :: ios
 real*8        :: x,y,z
 logical       :: ldum
 rewind(inpunit)
 ldum=move_to_record(INPUT_FIELDS(ORB_F),inpunit)
 natom=0
 do 
  read(inpunit,'(A)')string
  read(string,*,iostat=ios)olab1,x,y,z
  if(ios.ne.0)exit
  natom=natom+1
 enddo
 rewind(inpunit)
end function count_atom




!-------------------------------------------------------------
!Given cartesian coordinates in 3D space returns coordinates
!in units of primitive cell vectors. 
!-------------------------------------------------------------
subroutine convert_to_fractional(xat,ainv)
 integer::h
 real*8,intent(inout)::xat(rdim)
 real*8, intent(in) :: ainv(rdim,rdim)
 real*8:: xc(rdim)
 xc(:)=xat(:)
 do h=1,3
  xat(h)=sum(ainv(h,:)*xc(:))
 enddo
end subroutine convert_to_fractional
 



!------------------------------------------------------------------------------
!Returns phase(iat). Each orbital iat has now a phase (often +1 or -1).
!In input one specify a supercell (that must be contained in the
!bigger supercell used in simulation) and a phase for each of the orbitals
!inside it. Translational symmetry is then used to transfer the phase
!to the entire system.
!------------------------------------------------------------------------------
subroutine assign_phase(lattice)
type(lattice_t)     :: lattice
integer             :: pc(rdim,rdim),i,iat,jat,j,natom_ph,ios,ndim,nsites,natom
real*8              :: rpc(rdim,rdim),volume_ph,inv(rdim,rdim),projph(rdim,rdim),&
 &                     diff(rdim),projph2(rdim),ainv(rdim,rdim)
character*50        :: string
character*3         :: olab1
real*8, allocatable :: xat_ph(:,:),tmp_phase(:)
real*8, pointer     :: phase(:) 
logical             :: phase_assigned(0:lattice%nsites-1)

ndim=lattice%ndim
nsites=lattice%nsites
natom=lattice%natom

if(move_to_record(INPUT_FIELDS(PHASE_F),inpunit))then

 !Read the phase cell (pc)
 pc(:,:) = 0
 if(ndim>0)then
  read(inpunit,*)((pc(i,j),i=1,ndim),j=1,ndim)
  do i=ndim+1,rdim; pc(i,i)=1; enddo
 endif

 !Compute the volume of phase cell
 rpc(1:rdim,1:rdim)=dble(pc(1:rdim,1:rdim))
 volume_ph=abs(get_det(rpc))

 !check that the supercell is a multiple of the phase cell
 call get_inverse(rpc,inv)
 do i=1,rdim 
  do j=1,rdim
   projph(i,j)=sum(lattice%sc(:,j)*inv(i,:))
  enddo 
 enddo
 projph(:,:)=(projph(:,:)-nint(projph(:,:)))**2
 if(sqrt(sum(projph))>toll)stop 'Supercell frustrates the phase. Stop.'

 !Read atom positions and their phases
 natom_ph=nint(volume_ph*natom)
 allocate(phase(0:nsites-1),xat_ph(rdim,0:natom_ph-1),tmp_phase(0:natom_ph-1))
 do iat=0,natom_ph-1
  read(inpunit,'(A)')string
  read(string,*,iostat=ios)olab1,(xat_ph(i,iat),i=1,rdim),tmp_phase(iat)
  if(ios.ne.0)stop 'Problem with reading phases.'
 enddo

 !need to use inv corresponding to primitive cell
 call get_inverse(lattice%ac,ainv)
 do iat=0,natom_ph-1
  call convert_to_fractional(xat_ph(1,iat),ainv)
 enddo

 phase_assigned(:)=.false.
 !Assign the same phase to all sites equivalent by "phase-cell" translation 
 do iat=0,natom_ph-1
  do jat=0,nsites-1
   if(phase_assigned(jat))cycle
   !Determine whether site jat has the phase of iat
   diff(:)=xat_ph(:,iat)-lattice%pos(:,jat)
   do i=1,rdim
    projph2(i)=sum(diff(:)*inv(i,:))
   enddo
   projph2(:)=(projph2(:)-nint(projph2(:)))**2
   if(sqrt(sum(projph2))<toll)then
    !Assign to jat the phase of iat
    phase(jat)=tmp_phase(iat)
    phase_assigned(jat)=.true.
   endif
  enddo
 enddo

 !Store phase
 lattice%phase=>phase

 deallocate(xat_ph,tmp_phase)

endif
end subroutine assign_phase




!-----------------------------------------------------------------------------
! Given site and displacement in cartesian coordinates it returns the site
! to which the particle lands onto. The numerical label of the site is
! necessary since we are not excluding the case of two orbitals sitting
! at the same site (like Wannier function).
!----------------------------------------------------------------------------
integer function hoptowho(iat,delta,jat,lattice)
 type(lattice_t), intent(in) :: lattice
 integer,intent(in) :: iat,jat
 integer            :: i,j,ndim
 real*8, intent(in) :: delta(rdim)
 real*8             :: projk(rdim),xxat(rdim),invscc(rdim,rdim),xx(rdim)
 real*8, pointer    :: cartpos(:,:)

 if(.not.lattice%constructed)stop'Need to construct lattice before hoptowho'

 cartpos=>lattice%cartpos
 ndim=lattice%ndim
 call get_inverse(lattice%scc,invscc)
 !determine position after hopping
 xxat(:)=cartpos(:,iat)+delta(:)
 !Try to determine whether xxat corresponds to a site
 do j=jat,lattice%nsites-1,lattice%natom
  !compute distance vector in units of the supercell vectors
  xx(:)=cartpos(:,j)-xxat(:)
  do i=1,rdim
    projk(i)=sum(xx(:)*invscc(i,:))
  enddo
  projk(1:ndim)=projk(1:ndim)-nint(projk(1:ndim))
  if(sum(projk(:)**2)<toll)exit
 enddo
 if(j>=lattice%nsites)then
  !No site was found
  write(*,*)'Can''t find where',iat,' hops.'
  stop
 else
  hoptowho=j
 endif
end function



subroutine assign_gf_phase(lattice,twist)
type(lattice_t), intent(inout) :: lattice
real*8, intent(in) :: twist(3)
integer :: i, j, n, natom, ic, csize
real*8 :: d(3), d0(3), rphase, iphase
n=lattice%nsites
natom=lattice%natom
allocate(lattice%gf_phase(0:n-1,0:n-1))
do ic=1,lattice%nclass
   csize=0
   do j=0,n-1
      do i=0,n-1
         if(lattice%myclass(i,j)/=ic)cycle
         d=lattice%translation(:,j/natom)-lattice%translation(:,i/natom)
         if(csize==0)then
            d0=d
            lattice%gf_phase(j,i)=1
         else
            rphase=cos(sum(twist*(d0-d)))
            iphase=sin(sum(twist*(d0-d)))
            if(abs(iphase)<1.d-6)then
               if(abs(rphase-nint(rphase))<1.d-6)then
                  lattice%gf_phase(j,i)=nint(rphase)
               endif
            else
               rphase=cos(sum(twist*(d0+d)))
               iphase=sin(sum(twist*(d0+d)))
               if(abs(iphase)<1.d-6)then
                  if(abs(rphase-nint(rphase))<1.d-6)then
                     lattice%gf_phase(j,i)=nint(rphase)
                  else
                     stop'problem with phase'
                  endif
               endif
            endif
         endif
         csize=csize+1
      enddo
   enddo
   if(csize/=lattice%class_size(ic))stop'problem with classes'
enddo
end subroutine


!type(lattice_t), intent(inout) :: lattice
!real*8, intent(in) :: twist(3)
!integer :: i, j, k, n, natom, ic, jt, i0, j0, iat, jat
!real*8 :: R(3), delta(3)
!complex*16 :: phase1, phase2
!n=lattice%nsites
!natom=lattice%natom
!allocate(lattice%gf_phase(0:n-1,0:n-1))
!!Loop over cell
!do ic=0,lattice%ncell-1
!   !loop over translation
!   do jt=0,lattice%ncell-1
!      delta(:)=lattice%translation(:,jt)
!      !Find the two sites (orbital 0) connected by delta
!      i0=natom*ic
!      j0=hoptowho(i0,delta,0,lattice)
!      !Compute the winding vector
!      R=lattice%cartpos(:,j0)-lattice%cartpos(:,i0)-delta
!      phase=exp(-im*sum(R*twist))
!      !Loop over all orbitals contained in the two cells
!      do iat=0,natom-1
!         do jat=0,natom-1
!            i=i0+iat
!            j=j0+jat
!            lattice%gf_phase(j,i)=phase
!         enddo
!      enddo
!   enddo
!enddo   


end module DQMC_LATT

