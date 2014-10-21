module DQMC_SYMM

use DQMC_GEOM_PARAM
use DQMC_UTIL

implicit none

type :: symm_operations
 integer               :: ntotsymm              !number of symmetry read from input
 integer               :: nsymm                 !number of symmetries compatible with supercell
 integer               :: ntransl               !number of translations
 integer, pointer      :: map_symm(:,:)         !action of symmetry on site (0:nsites-1, nsymm)
                                                !map_symm(i,j) returns the site where site "i" is mapped 
                                                !by symmetry operation 'j"
 integer, pointer      :: map_symm_k(:,:)       !action of symmetry on k-point (nkpts, msymm)
                                                !map_symm_k(i,j) returns the k-point where k-point "i" is mapped 
                                                !by symmetry operation 'j"
 integer, pointer      :: map_symm_g(:,:)       !action of symmetry on k-point (nkpts, msymm)
                                                !map_symm_k(i,j) returns the k-point where k-point "i" is mapped 
                                                !by symmetry operation 'j"
                                                !msymm is equal to nsymm if addTimeRev==false. Otherwise msymm=nsymm+1
 integer, pointer      :: map_symm_b(:,:)       !action of symmetry on bond (ntotbond, nsymm)
                                                !A bond is a pair of sites where the first site is inside the unit cell
 integer, pointer      :: map_symm_p(:,:)       !action of symmetry on pair (nbond, nsymm)
                                                !Pairs are a subset of bonds used to describe pairing properties.
                                                !Their symmetry property are derived from those of bonds.
 integer, pointer      :: translback(:)         !number of the translation mapping a site into its untranslated 
                                                !image inside the primitive cell (0:nsites-1) 
 integer, pointer      :: translate(:,:)        !action of translations on site (0:nsites-1,0:ntransl-1)
 integer, pointer      :: valid_symm(:)         !label of valid symmetry operation(nsymm)
                                                !It returns a number between 1 and ntotsymm
 real*8, pointer       :: symmangle(:)          !Angle of rotation of symmetry operation (ntotsymm) in radiants
                                                !rotations
 real*8, pointer       :: symmpoint(:,:)        !point which symmetry operation goes through (rdim,ntotsymm)
                                                !rotations, reflections and inversions
 real*8, pointer       :: symmaxis(:,:)         !direction of the axis (rdim,ntotsymm)
                                                !rotations and reflections
 character*1, pointer  :: symmlabel(:)          !label of operation: C, D, I ,S (ntotsymm)
                                                !rotation, mirror plane, inversion center, rotoreflection
 
 logical               :: initialized
 logical               :: lattice_mapped
 logical               :: recip_lattice_mapped
 logical               :: bonds_mapped
 logical               :: addTimeRev

end type



contains

!------------------------------------------------------------------------------


 subroutine free_symm(symm)

  type(symm_operations), intent(inout) :: symm

    if(associated(symm%map_symm))    deallocate(symm%map_symm)
    if(associated(symm%symmaxis))    deallocate(symm%symmaxis)
    if(associated(symm%symmlabel))   deallocate(symm%symmlabel)
    if(associated(symm%translate))   deallocate(symm%translate)
    if(associated(symm%symmangle))   deallocate(symm%symmangle)
    if(associated(symm%symmpoint))   deallocate(symm%symmpoint)
    if(associated(symm%map_symm_k))  deallocate(symm%map_symm_k)
    if(associated(symm%map_symm_g))  deallocate(symm%map_symm_g)
    if(associated(symm%map_symm_b))  deallocate(symm%map_symm_b)
    if(associated(symm%map_symm_p))  deallocate(symm%map_symm_p)
    if(associated(symm%translback))  deallocate(symm%translback)
    if(associated(symm%valid_symm))  deallocate(symm%valid_symm)

 end subroutine free_symm


!------------------------------------------------------------------------------
! Read the point-symmetries from input.
! Rotation-axis : CN x y z x1 y1 z1
!  2\pi/N is the rotation angle
!  x y z specify a point belonging to the axis in cartesian coordinates
!  x1 y1 z1 specify the axis direction in cartesian coordinates
! Rotorefelction-axis : SN x y z x1 y1 z1
!  2\pi/N is the rotation angle
!  x y z specify a point belonging to the axis in cartesian coordinates
!  x1 y1 z1 specify the axis direction in cartesian coordinates
! Mirror plane : D x y z x1 y1 z1
!  x y z specify a point belonging to the plane in cartesian coordinate
!  x1 y1 z1 specify the direction normal to the plane in cartesian coordinates
! Inversion : I x y z
!  x y z specify position of inversion point in cartesian coordinates
!------------------------------------------------------------------
subroutine read_symm(SymmOp)
integer isymm,axis_order,ios1,ios2,i,nsymm
real*8 xpoint(3),xaxis(3),xnorm
character*50 string
character*1 label
logical dum
character*1, pointer :: symmlabel(:)  
real*8, pointer      :: symmangle(:)  
real*8, pointer      :: symmpoint(:,:)
real*8, pointer      :: symmaxis(:,:) 
type(symm_operations),intent(out) :: SymmOp

!Count the symmetry operation specified in input
nsymm=count_symmetry()
SymmOp%initialized=.false.
SymmOp%nsymm=nsymm
SymmOp%ntotsymm=nsymm
SymmOp%addTimeRev=.true.

if(nsymm/=0)then
 !allocate angle, points and axis for all symmetry operations
 allocate(symmangle(nsymm),symmpoint(3,nsymm),symmaxis(3,nsymm),symmlabel(nsymm))
 dum=move_to_record(INPUT_FIELDS(SYMM_F),inpunit)
 do isymm=1,nsymm
  !read entire line
  read(inpunit,'(A)')string
  string=adjustl(string)
  !read very first character
  read(string,'(A1)')label

  if(label=='C'.or.label=='c'.or.label=='S'.or.label=='s')then
   !Found a rotation or rotoreflection
   read(string,'(A1,i1)',iostat=ios1)label,axis_order
   read(string(4:50),*,iostat=ios2)(xpoint(i),i=1,3),(xaxis(i),i=1,3)
   if(ios1.ne.0.or.ios2.ne.0)then
    write(*,*)'Problem reading axis symmetry. Stop. Line:',isymm
    stop
   endif
   symmlabel(isymm)=label(1:1)
   symmpoint(:,isymm)=xpoint(:)
   !angle is 2pi/order
   symmangle(isymm)=4.d0*acos(0.d0)/axis_order
   !Save axis after normalizing it to 1.d0.
   xnorm=sqrt(sum(xaxis(:)**2))
   symmaxis(:,isymm)=xaxis(:)/xnorm

  elseif(label=='D'.or.label=='d')then
   !Found mirror plane
   read(string(3:50),*,iostat=ios2)(xpoint(i),i=1,3),(xaxis(i),i=1,3)
   if(ios2.ne.0)then
    write(*,*)'Problem reading plane symmetry. Stop. Line:',isymm
    stop
   endif
   symmlabel(isymm)=label(1:1)
   symmpoint(:,isymm)=xpoint(:)
   !Save axis after normalizing it to 1.d0.
   xnorm=sqrt(sum(xaxis(:)**2))
   symmaxis(:,isymm)=xaxis(:)/xnorm
   !Setting angle to 0.d0
   symmangle(isymm)=0.d0

  elseif(label=='I'.or.label=='i')then
   read(string(3:50),*,iostat=ios2)(xpoint(i),i=1,3)
   if(ios2.ne.0)then
    write(*,*)'Problem reading parity symmetry. Stop. Line:',isymm
    stop
   endif
   symmlabel(isymm)=label(1:1)
   symmpoint(:,isymm)=xpoint(:)
   !Setting angle and axis to 0.d0
   symmangle(isymm)=0.d0
   symmaxis(:,isymm)=0.d0
   SymmOp%addTimeRev=.false.

  else
   write(*,*)'Unknown label for symmetry. Stop'
   stop
  endif
 enddo

 !Fill SymmOp
 SymmOp%initialized=.true.
 SymmOp%nsymm=nsymm
 SymmOp%ntotsymm=nsymm
 !Save pointers
 SymmOp%symmlabel=>symmlabel
 SymmOp%symmaxis=>symmaxis
 SymmOp%symmpoint=>symmpoint
 SymmOp%symmangle=>symmangle

endif
end subroutine read_symm



!--------------------------------------------------------
!Count how many point-symmetries are specified in input
!--------------------------------------------------------
integer function count_symmetry()
integer msymm
character*1 label
character*30 string
logical ldum
msymm=0
if(Found_Field(SYMM_F))then
 ldum=move_to_record(INPUT_FIELDS(SYMM_F),inpunit)
 do
  read(inpunit,'(A)')string
  string=adjustl(string)
  read(string,'(A1)')label
  if(   label=='C'.or.label=='c' &
  & .or.label=='S'.or.label=='s' &
  & .or.label=='D'.or.label=='d' &
  & .or.label=='I'.or.label=='i' )then
   msymm=msymm+1
  else
   exit
  endif
 enddo
else
 write(*,*)'No symmetries have been specified'
endif
count_symmetry=msymm
end function count_symmetry




!----------------------------------------------------------------------
! Apply point-symmetry operation. label specifis the kind of
! symmetry (C, S, D or I). Point and axis specify the "position" of
! the symmetry operation (see comments to read_symmetry). set is a set
! of 3D vectors in cartesian coordinates upon which the symmetry acts.
! newset returns the transformed set.
!-----------------------------------------------------------------------
subroutine apply_point_symm(label,point,axis,theta,set,new_set,nset,reciprocal)
  ! label : C, S, D or I
  ! point, axis and theta specify the symmetry operation
  ! set is the set of point to be transformed (the lattice in real or k-space)
  ! nset is the number of points in set
  ! reciprocal is a flag specifying whether the set is made of k-vectors
integer,intent(in) :: nset
real*8,intent(in) :: point(3),axis(3),theta,set(3*nset)
logical, intent(in) :: reciprocal
real*8,intent(out) :: new_set(3,nset)
integer irefl,iset, alp, bet, gam, del, i, j
real*8 trans_set(3,nset),d,a1,a2,rot1(3,3),rot2(3,3),rot3(3,3),globrot(3,3), &
& cost,sint,dpoint(3),dset(3,nset),dset2(nset,3)
character*1 label

!If transformation is in k-space symmetry operation necessarily goes through 0.d0
!dpoint is the "origin" 
if(reciprocal)then
 dpoint(:)=0.d0
 dset=reshape(set,shape=(/3,nset/),order=(/2,1/))
else
 dpoint(:)=point(:)
 dset=reshape(set,shape=(/3,nset/),order=(/1,2/))
endif

if(label=='i'.or.label=='I')then
 !Deal with inversion around dpoint
 do iset=1,nset
  new_set(:,iset)=2*dpoint(:)-dset(:,iset)
 enddo
else
 !Deal with rotations, reflections or rotoreflections
 if(label=='c'.or.label=='C')then 
  irefl=1
 else
  irefl=-1
 endif

 !Note that cost and sint are 1 and 0 for mirror planes
 cost=cos(theta); sint=sin(theta)

 !translate the set
 do iset=1,nset
  trans_set(:,iset)=dset(:,iset)-dpoint(:)
 enddo

 !Define matrix that line up z-axis with rotation axis
 d=sqrt(axis(1)**2+axis(2)**2)
 if(d>1.d-6)then
  a1=axis(1)/d; a2=axis(2)/d
 else
  a1=1.d0; a2=0.d0
 endif

 !First rotate reference frame around z : R1
 rot1(1,1)= a1       ; rot1(1,2)= a2       ; rot1(1,3)=0.d0
 rot1(2,1)=-a2       ; rot1(2,2)= a1       ; rot1(2,3)=0.d0
 rot1(3,1)= 0.d0     ; rot1(3,2)= 0.d0     ; rot1(3,3)=1.d0
 !then Rotate reference frame around y : R2
 rot2(1,1)= axis(3)  ; rot2(1,2)= 0.d0     ; rot2(1,3)= -d
 rot2(2,1)= 0.d0     ; rot2(2,2)= 1.d0     ; rot2(2,3)= 0.d0
 rot2(3,1)= d        ; rot2(3,2)= 0.d0     ; rot2(3,3)= axis(3)
 !Finally rotate around z and/or reflect perpendicularly to xy (when irefl=-1): R3
 rot3(1,1)= cost      ; rot3(1,2)=-sint    ; rot3(1,3)=0.d0
 rot3(2,1)= sint      ; rot3(2,2)= cost    ; rot3(2,3)=0.d0
 rot3(3,1)= 0.d0      ; rot3(3,2)= 0.d0    ; rot3(3,3)=1.d0*irefl
 !compute globrot as R1^-1 R2^-1 R3 R2 R1 (undoing frame rotations)
 globrot(:,:)=0.d0
 do i=1,3; do j=1,3
  do alp=1,3; do bet=1,3; do gam=1,3; do del=1,3
   globrot(i,j)=globrot(i,j)+rot1(alp,i)*rot2(bet,alp)*rot3(bet,gam)*rot2(gam,del)*rot1(del,j)
  enddo; enddo; enddo; enddo
 enddo; enddo

 !Apply symmetry operation
 new_set(:,:)=0.d0
 do iset=1,nset
  do i=1,3
   do j=1,3
    new_set(i,iset)=new_set(i,iset)+globrot(i,j)*trans_set(j,iset)
   enddo
  enddo
 enddo

 !Translate back to the origin
 do iset=1,nset
  new_set(:,iset)=new_set(:,iset)+dpoint(:)
 enddo

endif

!If reciprocal transpose the matrix to follow the convention
if(reciprocal)then
 !Fancy way to get transpose
 dset2=reshape(new_set,shape=(/nset,3/),order=(/2,1/))
 new_set=reshape(dset2,shape=(/3,nset/))
endif

end subroutine apply_point_symm




!---------------------------------------------------------------------
! Map the action of the point-operation into the vector map_symm:
!map_symm(ifrom,isymm) returns the orbital on which orbital ifrom
!is transformed by isymm. Symmetries that are incompatible with the
!supercell are discarded. msymm is the number of compatible point
!symmetries.
! Construct translate(ifrom,itrans). Analogous to map_symm but returns
!the orbital to which ifrom is mapped by translation itrans. 
! Construct translback(ifrom). It returns the translation that applied
!to ifrom reduces the latter to the equivalent orbital inside the
!primitive cell.
!---------------------------------------------------------------------
subroutine map_symm_lattice(SymmOp,lattice, hamilt, SOP)
use DQMC_LATT
use DQMC_HAMILT

type(symm_operations)          :: SymmOp
type(lattice_t),intent(in)     :: lattice
type(hamiltonian_t),intent(in) :: hamilt
integer, intent(in)            :: SOP

integer :: i, istart, ipat, jpat, iat, jat, ii, it, itl, is
integer :: nsites, ntotsymm, ndim, msymm, nsymm, ntransl, natom
integer :: valid_symm(SymmOp%ntotsymm)
integer :: tmp_symm(0:lattice%nsites-1,SymmOp%ntotsymm) 

real*8  :: diff(rdim), projsc(rdim), invscc(rdim,rdim)
real*8  :: newpos(rdim,0:lattice%nsites-1)
real*8, pointer  :: symmangle(:)  
real*8, pointer  :: symmpoint(:,:)
real*8, pointer  :: symmaxis(:,:) 

logical :: mapped(0:lattice%nsites-1,SymmOp%ntotsymm)
logical :: mappedt(0:lattice%nsites-1,0:lattice%ncell-1)
logical :: equal

character*1, pointer :: symmlabel(:)

!Assign pointers
symmlabel => SymmOp%symmlabel 
symmangle => SymmOp%symmangle
symmpoint => SymmOp%symmpoint 
symmaxis  => SymmOp%symmaxis

!Initializa local variables
ntransl  = lattice%ncell
nsites   = lattice%nsites
natom    = lattice%natom
ntotsymm = SymmOp%ntotsymm
ndim     = lattice%ndim

call get_inverse(lattice%scc, invscc)

!Apply symmetry operations to sites...
mapped(:,:) = .false.
do i = 1, ntotsymm
   call apply_point_symm(symmlabel(i),symmpoint(:,i),symmaxis(:,i),symmangle(i),lattice%cartpos,newpos,nsites,.false.)
   !Find where "iat" was mapped by the symmetry operation
   do iat = 0, nsites-1
      ipat = mod(iat,natom)
      !Compare the new position (newpos) with the original one of all sites
      do jat = 0, nsites-1
         jpat = mod(iat,natom)
         !if jat was already "mapped onto" or label is different, skip jat
         if(mapped(jat,i) .or. lattice%olabel(jpat)/=lattice%olabel(ipat)) cycle
         !Compute distance in units of the supercell vectors
         diff(:) = newpos(:,iat) - lattice%cartpos(:,jat)
         do ii = 1, 3
            projsc(ii) = sum(diff(:)*invscc(ii,:))
         enddo
         !project inside the supercell - only along the "extended" dimensions
         projsc(1:ndim) = projsc(1:ndim) - nint(projsc(1:ndim))
         if(sqrt(sum(projsc(1:rdim)**2))<toll)then
            tmp_symm(iat,i) = jat
            mapped(jat,i)   = .true.
            exit
         endif
      enddo
   enddo
enddo

!Check whether symmetry is compatible with supercell 
msymm=0
do i = 1, ntotsymm
  do iat = 0, nsites-1
    if(.not.mapped(iat,i))exit
  enddo
  !if all components of mapped(:,i) are true we have a valid symmetry
  if(iat == nsites)then
    msymm = msymm + 1
    valid_symm(msymm) = i
  endif
enddo

!Check whether symmetry is compatible with hamiltonian
nsymm = msymm
do i = msymm, 1, -1
   ii = valid_symm(i)
   equal=.true.
   !Transform the 3 pieces of the hamiltonian
   do iat = 0, nsites-1
      ipat = tmp_symm(iat,ii)
      do jat = 0, nsites-1
         jpat = tmp_symm(jat,ii)
         equal = equal.and.( abs(hamilt%hopup(ipat,jpat)-hamilt%hopup(iat,jat))< 1.d-3 )
         equal = equal.and.( abs(hamilt%hopdn(ipat,jpat)-hamilt%hopdn(iat,jat))< 1.d-3 )
         equal = equal.and.( abs(hamilt%Uv(ipat,jpat)-hamilt%Uv(iat,jat))< 1.d-3 )
         equal = equal.and.( abs(hamilt%Jv(ipat,jpat)-hamilt%Jv(iat,jat))< 1.d-3 )
      enddo
   enddo
   !Exclude symmetry if incompatible with H
   if(.not.equal)then
    valid_symm(i) = valid_symm(nsymm)
    nsymm = nsymm - 1
   endif
enddo
msymm = nsymm

!Store which symmetries are valid and their action in the map_symm table
if(msymm > 0) &
  & allocate(SymmOp%map_symm(0:nsites-1,msymm),SymmOp%valid_symm(msymm))
do i = 1, msymm
  SymmOp%valid_symm(i) = valid_symm(i)
  SymmOp%map_symm(:,i) = tmp_symm(:,valid_symm(i))
enddo
!Save the number of valid symmetry operations
SymmOp%nsymm = msymm


!Deal with translation
SymmOp%ntransl = ntransl
if(ntransl>0)allocate(SymmOp%translate(0:nsites-1,0:ntransl-1),SymmOp%translback(0:nsites-1))
mappedt(:,:) = .false.
do it = 0, ntransl-1
  itl = it*natom
  !Apply translation to each site
  do iat = 0, nsites-1
    newpos(:,iat) = lattice%cartpos(:,iat) + lattice%translation(:,it)
  enddo
  !Find the site "jat" where "iat" is mapped by translation "it"
  do iat = 0, nsites-1
    istart = mod(iat,natom)
    do jat = istart, nsites-1, natom
      if(mappedt(jat,it))cycle
      diff(:) = newpos(:,iat) - lattice%cartpos(:,jat)
      do ii = 1, 3
        projsc(ii) = sum(diff(:) * invscc(ii,:))
      enddo
      projsc(1:ndim) = projsc(1:ndim) - nint(projsc(1:ndim))
      if(sum(projsc(1:rdim)**2) < toll)then
        SymmOp%translate(iat,it) = jat
        mappedt(jat,it)          = .true.
        exit
      endif
    enddo
  enddo
enddo

!Check correctness and construct "translback"
do is = 0, nsites-1
 do it = 0, ntransl-1
  !Check that mapping was done properly
  if(.not.mappedt(is,it))then
    write(*,'(A,i3,A)')'Translation ',it,' could not be mapped correctly'
    stop
  endif
  !Store the translation that maps "is" inside the primitive cell
  iat = SymmOp%translate(is,it)
  if(iat < natom)then
   SymmOp%translback(is) = it
   exit
  endif
 enddo
enddo

!Write out symmetry mappings
write(SOP,*)'Number of symmetry operations in input: ',ntotsymm
write(SOP,*)'Number of valid symmetry operations   : ',msymm
do i = 1, msymm
  write(SOP,*)'Mapping of Symmetry :',valid_symm(i),Symmlabel(valid_symm(i))
  do iat = 0, nsites-1
    write(SOP,*)iat,'->',SymmOp%map_symm(iat,i)
  enddo
enddo
do it = 0, ntransl-1
  write(SOP,*)'Mapping of Translation',it
  do iat = 0, nsites-1
    write(SOP,*)iat,'->',SymmOp%translate(iat,it)
  enddo
enddo
write(SOP,*)'Label of Translation mapping site in primitive cell'
do iat = 0, nsites-1
  write(SOP,*)iat,'->',SymmOp%translback(iat)
enddo

SymmOp%lattice_mapped = .true.

end subroutine map_symm_lattice




!--------------------------------------------------------------------------------
!k-space symmetry
!--------------------------------------------------------------------------------
subroutine map_symm_recip_lattice(SymmOp,recip_lattice,applytwist)
use DQMC_RECLATT
type(symm_operations)              :: SymmOp
type(recip_lattice_t),intent(in)   :: recip_lattice
logical, intent(in)                :: applytwist
integer                            :: i,j,ii,nsymm,ndim,ik,jk
real*8                             :: diff(rdim),projsc(rdim),invkc(rdim,rdim), &
                                      newklist(recip_lattice%nkpts,rdim),zerovec(rdim)
logical                            :: mapped(recip_lattice%nkpts), includesymm
character*1, pointer               :: symmlabel(:)  
real*8, pointer                    :: symmangle(:)  
real*8, pointer                    :: symmpoint(:,:)
real*8, pointer                    :: symmaxis(:,:) 
integer :: msymm, nsites, valid_symm(SymmOp%nsymm)
integer, pointer :: tmp_symm(:,:)  
integer, pointer :: tmp_symm_k(:,:)
integer, pointer :: tmp_valid(:)   

if(.not.SymmOp%lattice_mapped)stop'Need to map lattice symmetries before recip lattice ones'

!Assign pointers
symmlabel => SymmOp%symmlabel 
symmangle => SymmOp%symmangle
symmpoint => SymmOp%symmpoint
symmaxis  => SymmOp%symmaxis
!Only the "valid" symmetries will be considered
nsymm = SymmOp%nsymm

!Check whether we need to additionally apply time-reversal symmetry
!If there was no inversion add it for k-space
if(.not.SymmOp%addTimeRev)then
   do i=1,nsymm
      j=SymmOp%valid_symm(i)
      if(SymmOp%symmlabel(j)=='i'.or.SymmOp%symmlabel(j)=='I')exit
   enddo
   if(i>nsymm)SymmOp%addTimeRev=.true.
endif
if(SymmOp%addTimeRev) nsymm = nsymm + 1

ndim = recip_lattice%ndim
call get_inverse(recip_lattice%kc,invkc)
if(nsymm>0)then
   if(applytwist)then
      allocate(SymmOp%map_symm_k(recip_lattice%nkpts,nsymm))
   else
      allocate(SymmOp%map_symm_g(recip_lattice%nkpts,nsymm))
   endif
endif

!Loop over all valid symmetries
msymm = 0
do j = 1, SymmOp%nsymm
   i = SymmOp%valid_symm(j)
   mapped(:) = .false.
   includesymm = .true.
   !Apply point symmetry to all k-points
   call apply_point_symm(symmlabel(i),symmpoint(:,i),symmaxis(:,i),symmangle(i),&
    recip_lattice%klist,newklist,recip_lattice%nkpts,.true.)
   do ik = 1, recip_lattice%nkpts
      !Find which k-point, "jk", "ik" is mapped onto by symmetry "j" 
      do jk=1,recip_lattice%nkpts
         !For each k-point compute distance in unit of k-space unit cell
         if(mapped(jk))cycle
         diff(:) = newklist(ik,:) - recip_lattice%klist(jk,:)
         do ii = 1, rdim
            projsc(ii) = sum(diff(:)*invkc(:,ii))
         enddo
         !project inside the unit cell - only along the "extended" dimensions
         projsc(1:ndim) = projsc(1:ndim) - nint(projsc(1:ndim))
         if(sum(projsc(1:rdim)**2)<toll)then
            if(applytwist)then
               SymmOp%map_symm_k(ik,j) = jk
            else
               SymmOp%map_symm_g(ik,j) = jk
            endif
            mapped(jk) = .true.
            exit
         endif
      enddo
      if(jk>recip_lattice%nkpts)then
         if(applytwist)then
            includesymm = .false.
         else
            stop'Problem with symmetry in k-space'
         endif
      endif
   enddo
   if(includesymm)then
      msymm = msymm + 1
      valid_symm(msymm) = j
   endif
enddo

!Modify list of valid symmetry if twist spoils some of them
if(msymm /= SymmOp%nsymm.and.applytwist)then

  nsites = size(SymmOp%map_symm,1)

  allocate(tmp_symm(0:nsites-1,SymmOp%nsymm),tmp_symm_k(recip_lattice%nkpts,nsymm), &
           tmp_valid(SymmOp%nsymm))

  !Save symmetry mapping table
  tmp_symm   = SymmOp%map_symm
  tmp_symm_k = SymmOp%map_symm_k
  tmp_valid  = SymmOp%valid_symm

  !deallocate the old table
  deallocate(SymmOp%map_symm,SymmOp%map_symm_k,SymmOp%valid_symm)

  if(msymm>0)then

     !allocate new tables
     allocate(SymmOp%map_symm(0:nsites-1,msymm))
     allocate(SymmOp%valid_symm(msymm))
     nsymm = msymm
     if(SymmOp%addTimeRev) nsymm = nsymm + 1
     allocate(SymmOp%map_symm_k(recip_lattice%nkpts,nsymm))
   
     !copy valid symmetries
     do i = 1, msymm
        SymmOp%map_symm(:,i)   = tmp_symm(:,valid_symm(i))
        SymmOp%valid_symm(i)   = tmp_valid(valid_symm(i))
        SymmOp%map_symm_k(:,i) = tmp_symm_k(:,valid_symm(i))
     enddo

     !Save the number of valid symmetry operations
     SymmOp%nsymm = msymm
  endif

  deallocate(tmp_valid, tmp_symm, tmp_symm_k)

endif

!Add the additional symmetry time reversal symmetry when "I" is not a symmetry
!operation specified in input
if(SymmOp%addTimeRev)then
   mapped(:) = .false.
   zerovec(1:3)=(/0.d0, 0.d0, 0.d0/)
   call apply_point_symm('I', zerovec, zerovec, 0.d0,&
    recip_lattice%klist,newklist,recip_lattice%nkpts,.true.)
   do ik = 1, recip_lattice%nkpts
      !Find which k-point, "jk", "ik" is mapped onto by symmetry "j" 
      do jk = 1, recip_lattice%nkpts
         !For each k-point compute distance in unit of k-space unit cell
         if(mapped(jk))cycle
         diff(:) = newklist(ik,:)-recip_lattice%klist(jk,:)
         do ii = 1, rdim
            projsc(ii)=sum(diff(:)*invkc(:,ii))
         enddo
         !project inside the unit cell - only along the "extended" dimensions
         projsc(1:ndim) = projsc(1:ndim) - nint(projsc(1:ndim))
         if(sum(projsc(1:rdim)**2)<toll)then
          if(applytwist)then
             SymmOp%map_symm_k(ik,nsymm) = jk
          else
             SymmOp%map_symm_g(ik,nsymm) = jk
          endif
          mapped(jk) = .true.
          exit
         endif
      enddo
      if (jk > recip_lattice%nkpts) stop'Problem with symmetry in k-space'
   enddo
endif

SymmOp%recip_lattice_mapped = .true.

end subroutine map_symm_recip_lattice




!---------------------------------------------------------------------------------
! Given the list of bonds read from input, this routines does the following:
! 1) complete the list with the bonds which are equivalent, by symmetry,
!    to those specified in input.
! 2) Create a mapping (map_symm_b) that, given, bond "b" and a symmetry operation "s",
!    returns a the bond on which "b" is mapped by "s".
!---------------------------------------------------------------------------------
subroutine map_symm_bonds(Bonds,SymmOp,Lattice)
 use DQMC_BONDS
 type(symm_operations),intent(inout) :: SymmOp
 type(bonds_t),intent(inout) :: Bonds
 type(lattice_t),intent(in) :: Lattice
 integer :: nbclass,ibond,iat,jat,ntotpair,natom,nsites,bcl,newlabel,jbond,it, &
 & isymm,ndim
 integer, allocatable :: class(:),tag(:),pair_origin(:),pair_target(:),pair_label(:)
 integer, pointer :: bond_origin(:)
 integer, pointer :: bond_target(:)
 integer, pointer :: map_symm(:,:)
 real*8 :: invscc(rdim,rdim),proj(rdim)
 real*8, allocatable :: xxpair(:,:)
 logical, allocatable :: bond_on(:,:)

 if(.not.lattice%analyzed)stop'Need to analyze lattice before mapping bonds'
 if(Bonds%ntotbond==0)return

 natom=Lattice%natom; nsites=Lattice%nsites; ndim=Lattice%ndim
 call get_inverse(lattice%scc,invscc)
 nbclass=0
 newlabel=maxval(Bonds%bond_label)
 allocate(class(Bonds%ntotbond),tag(Bonds%ntotbond),bond_on(0:natom-1,0:nsites-1),bond_target(Bonds%ntotbond))
 bond_origin=>Bonds%bond_origin
 bond_on(:,:)=.false.

 !Assign to each bond a class based on lattice classes for distances
 do ibond=1,Bonds%ntotbond
  iat=bond_origin(ibond)
  jat=hoptowho(iat,Bonds%xxbond(1:rdim,ibond),Bonds%bond_target(ibond),Lattice)
  bond_target(ibond)=jat
  bond_on(iat,jat)=.true.
  class(ibond)=Lattice%myclass(iat,jat)
  !See if bond belongs to an already found class
  do jbond=1,ibond-1
   if(class(ibond)==class(jbond))exit
  enddo
  if(jbond==ibond)then
   !if not, create a new class
   nbclass=nbclass+1
   tag(nbclass)=class(ibond)
  endif
 enddo


 !Include all bonds which were left out but that are equivalent by symmetry...
 ntotpair=2*natom*nsites
 allocate(pair_label(ntotpair),pair_origin(ntotpair),pair_target(ntotpair),xxpair(3,ntotpair))
 !Save bond attribute in temporary "pair" variables 
 ntotpair=Bonds%ntotbond
 pair_origin(1:ntotpair)=bond_origin(1:ntotpair)
 pair_target(1:ntotpair)=bond_target(1:ntotpair)
 pair_label(1:ntotpair)=Bonds%bond_label(1:ntotpair)
 xxpair(:,1:ntotpair)=Bonds%xxbond(:,1:ntotpair)
 !... But only if there was no PAIR field specified in input
 if(.not.Found_Field(PAIRS_F))then
  !Loop over all pairs having the first atom inside the unit cell
  do iat=0,natom-1
   do jat=0,nsites-1
    !If bond is already "on" on this pair cycle
    if(bond_on(iat,jat))cycle
    !Otherwise loop over classes and found other bonds which belongs to same class
    do bcl=1,nbclass

     if(tag(bcl)==Lattice%myclass(iat,jat))then
      !We found a bond that needs to be included!
      !Assign a label to the bond
      newlabel=newlabel+1
      !Increase the number of total bonds
      ntotpair=ntotpair+1
      !Save its attribute
      pair_origin(ntotpair)=iat
      pair_target(ntotpair)=jat
      pair_label(ntotpair)=newlabel
      xxpair(:,ntotpair)=lattice%cartpos(:,jat)-lattice%cartpos(:,iat)
      !make xxpair as small as possible
      do it=1,rdim
       proj(it)=sum(xxpair(:,ntotpair)*invscc(it,:))
      enddo
      proj(1:ndim)=proj(1:ndim)-nint(proj(1:ndim))
      do it=1,rdim
       xxpair(it,ntotpair)=sum(lattice%scc(it,:)*proj(:))
      enddo
      !Switch on the bond flag
      bond_on(iat,jat)=.true.

      !Construct the opposite bond
      if(iat/=jat)then
       ntotpair=ntotpair+1
       pair_label(ntotpair)=-newlabel
       !Find the translation that maps iat inside the unit cell
       it=SymmOp%translback(jat) 
       !Translate jat and iat. jat returns the origin of the bond.
       pair_origin(ntotpair)=SymmOp%translate(jat,it)
       pair_target(ntotpair)=SymmOp%translate(iat,it)
       xxpair(:,ntotpair)=-xxpair(:,ntotpair-1)
       !Switch on the bond flag
       bond_on(SymmOp%translate(jat,it),SymmOp%translate(iat,it))=.true.
      endif
      !We switch this bons on. Stop looking over classes.
      exit
     endif

    enddo
   enddo
  enddo
 endif
 deallocate(class,tag,bond_target,bond_on)

 !Reload Bonds (new set completed with newly found bonds)
 if(Bonds%ntotbond/=ntotpair)then
  !First deallocate old stuff
  deallocate(Bonds%bond_origin,Bonds%xxbond,Bonds%bond_label,Bonds%bond_target)
  !Allocate and store temporary "pair" variable in "Bonds"
  allocate(Bonds%bond_origin(ntotpair),Bonds%bond_target(ntotpair), &
   Bonds%xxbond(3,ntotpair),Bonds%bond_label(ntotpair))
  Bonds%ntotbond=ntotpair
  Bonds%bond_origin(1:ntotpair)=pair_origin(1:ntotpair)
  Bonds%bond_target(1:ntotpair)=mod(pair_target(1:ntotpair),natom)
  Bonds%xxbond(:,1:ntotpair)=xxpair(:,1:ntotpair)
  Bonds%bond_label(1:ntotpair)=pair_label(1:ntotpair)
 endif

 !Store how a bond transforms under point symmetry
 allocate(map_symm(ntotpair,SymmOp%nsymm))
 do isymm=1,SymmOp%nsymm
  do ibond=1,ntotpair
   !Find the two sites where the bond origin and target are mapped into
   iat=SymmOp%map_symm(pair_origin(ibond),isymm)
   jat=SymmOp%map_symm(pair_target(ibond),isymm)
   !Translate them back do that origin is inside unit cell
   it=SymmOp%translback(iat) 
   iat=SymmOp%translate(iat,it)
   jat=SymmOp%translate(jat,it)
   !Find which bond is defined by (iat,jat) 
   do jbond=1,ntotpair
    if(pair_origin(jbond)==iat.and.pair_target(jbond)==jat)exit
   enddo
   if(jbond>ntotpair)stop 'Symmetry analysis : Cannot find equivalent bond'
   !Save the action of the symmetry operation
   map_symm(ibond,isymm)=jbond
  enddo
 enddo
 deallocate(pair_origin,pair_target,pair_label,xxpair)
 SymmOp%map_symm_b=>map_symm
 SymmOp%bonds_mapped=.true.

 !Write Info
 ! write(*,*)
 ! write(*,*)'Bonds (Set completed using symmetry)'
 ! do ibond=1,Bonds%ntotbond
 !  write(*,'(3i4,3f12.7)')ibond,Bonds%bond_label(ibond),Bonds%bond_origin(ibond),Bonds%xxbond(1:rdim,ibond)
 ! enddo
 !write(*,*)
 !write(*,*)'BOND MAPPING'
 !do isymm=1,SymmOp%nsymm
 ! write(*,*)'Symmetry',isymm
 ! do ibond=1,Bonds%ntotbond
 !  write(*,*)ibond,'-->',map_symm(ibond,isymm)
 ! enddo
 !enddo
 ! write(*,*)'====================================================================='

 Bonds%analyzed=.true.

end subroutine
 



!--------------------------------------------------------------------------------
! Map symmetry for pairs.
!--------------------------------------------------------------------------------
subroutine map_symm_pairs(Pairs, SymmOp)
 use DQMC_BONDS
 type(symm_operations),intent(inout) :: SymmOp
 type(pairing),intent(inout) :: Pairs
 integer :: isymm, ib, jb, newjb
 allocate(SymmOp%map_symm_p(Pairs%nbond,SymmOp%nsymm))
 do isymm=1,SymmOp%nsymm
   !loops over bonds defining the pairs
   do ib=1,Pairs%nbond
     !find the bond
     jb=Pairs%bond_map(ib)
     !map the bond
     newjb=SymmOp%map_symm_b(jb,isymm)
     !map the bond back into the pair list
     jb=Pairs%pair_map(newjb)
     if(jb==0)stop 'Symmetry analysis : Cannot find equivalent pair'
     !save it
     SymmOp%map_symm_p(ib,isymm)=jb
   enddo
 enddo
end subroutine




!------------------------------------------------------------------------------
!  Construct myclass(i,j). Given to sites i and j (not necessarily different)
! returns the class to which they belong. A class contains pairs of orbitals
! that, because of symmetry, are going to have identical pair-correlation
! functions.
!  Returns nclass, the number of classes, and class_size(iclass), the 
! number of pairs inside class iclass.
!  Returns class_label. This is a 4-components array. The first three are the
! cartesian separation of the two orbitals in the pair. The last component is
! the number of the atom inside the primitive cell that belongs to the pair.
!------------------------------------------------------------------------------
subroutine construct_lattice_classes(SymmOp,lattice)
use DQMC_LATT
integer               :: i,it,ip,is,j,isymm,iclass,istart,csize,csizenew,itransl,&
                         &jclass,idj,id,mclass, ip_transl,is_transl,ip2,is2,jstart,&
                         &nclass,nsites,natom,nsymm,ntransl
integer,allocatable   ::  patom(:,:),satom(:,:),csizev(:)
integer, pointer      :: myclass(:,:)
type(symm_operations) :: SymmOp
type(lattice_t)       :: lattice

if(.not.SymmOp%lattice_mapped)stop'Need to map symmetries over lattice before classes'

!initialize local variables
natom=lattice%natom
nsites=lattice%nsites
nsymm=SymmOp%nsymm
ntransl=SymmOp%ntransl

!allocate internal arrays
nclass=(natom*(natom+1))/2+natom*(nsites-natom)
allocate(patom(2,nclass),satom(2,nclass),csizev(nclass))
                                                
!At the beginning each distance is a separate class and only
!pairs with at least one atom in the primitive cell are considered.
!The pair (patom(ix,iclass) , satom(ix,iclas)) is the ix-th
!element of class "iclass"
allocate(myclass(0:nsites-1,0:nsites-1))
nclass=0
do ip=0,natom-1
 !first loop over sites inside primitive cell
 do is=ip,natom-1
  nclass=nclass+1
  if(is/=ip)then 
   csizev(nclass)=2
   patom(1,nclass)=ip; satom(1,nclass)=is
   patom(2,nclass)=is; satom(2,nclass)=ip
   myclass(ip,is)=nclass; myclass(is,ip)=nclass
  else 
   csizev(nclass)=1 
   patom(1,nclass)=ip; satom(1,nclass)=is
   myclass(ip,is)=nclass
  endif
 enddo
 do is=natom,nsites-1
  nclass=nclass+1
  patom(1,nclass)=ip; satom(1,nclass)=is
  csizev(nclass)=1
  myclass(ip,is)=nclass
 enddo
enddo

!Loop over symmetry operations. The "+1" symm op is pair permutation.
do isymm=1,nsymm+1
!do isymm=1,nsymm
 !Loop over classes (for the first symm op, classes are made
 !of all individual atom pairs in which the first atom lies
 !in the primitive cell and the second anywhere inside the supercell)
 do iclass=1,nclass
  istart=1
  do 
   csize=csizev(iclass)
   csizenew=csize
   do id=istart,csize
    !Map the atoms in the class under the symm operation
    if(isymm==nsymm+1)then
     ip=satom(id,iclass)
     is=patom(id,iclass)
    else
     ip=SymmOp%map_symm(patom(id,iclass),isymm)
     is=SymmOp%map_symm(satom(id,iclass),isymm)
    endif
    !Find the transformed pair translated such that the firts site is inside the
    !primitive cell
    itransl=SymmOp%translback(ip)
    ip2=SymmOp%translate(ip,itransl)
    is2=SymmOp%translate(is,itransl)
    !Find the class to which the pair belongs to
    jclass=myclass(ip2,is2)
    !if classes are different they need to be merged
    if(jclass/=iclass)then
     !all pairs of class jclass are transfered in class "iclass"
     jstart=csizenew
     csizenew=csizenew+csizev(jclass)
     call resize_class()
     do idj=1,csizev(jclass)
      ip=patom(idj,jclass); is=satom(idj,jclass)
      myclass(ip,is)=iclass
      patom(jstart+idj,iclass)=ip
      satom(jstart+idj,iclass)=is
     enddo
     !Size of jclass is nullified
     csizev(jclass)=0
    endif
   enddo
   !if class size did not change we have found all the classes equivalent 
   !to "iclass"
   if(csizenew==csize)exit
   !Update loop bounds to find new equivalence due to newly added elements
   istart=csizev(iclass)+1
   csizev(iclass)=csizenew
  enddo
 enddo
enddo

!Assign a class to all the remaining pair of atoms using translational symmetry
!redifine nclass as the number of final classes
mclass=0
do i=1,nclass
 if(csizev(i)>0)mclass=mclass+1
 do j=1,csizev(i)
  ip=patom(j,i)
  is=satom(j,i)
  myclass(ip,is)=mclass
  do it=1,ntransl-1
   ip_transl=SymmOp%translate(ip,it)
   is_transl=SymmOp%translate(is,it)
   myclass(ip_transl,is_transl)=mclass
  enddo
 enddo
enddo
deallocate(patom,satom,csizev)
nclass=mclass

!Define a class label using the distance between the sites
allocate(lattice%class_label(nclass,5))
do iclass=1,nclass
 prim:do i=0,natom-1
  super:do j=0,nsites-1
   if(myclass(i,j)==iclass)then
    lattice%class_label(iclass,1:3)=lattice%cartpos(1:3,j)-lattice%cartpos(1:3,i)
    ! 03/26/2013: added the second orbital index in lattice%class_label(iclass,5)
    lattice%class_label(iclass,4)=dble(i)
    lattice%class_label(iclass,5)=mod(j,natom)
    exit prim
   endif
  enddo super
 enddo prim
enddo

!Compute the size of each class
allocate(lattice%class_size(nclass))
lattice%class_size(:)=0
do i=0,nsites-1
 do j=0,nsites-1
  iclass=myclass(i,j)
  lattice%class_size(iclass)=lattice%class_size(iclass)+1
 enddo
enddo

!Save parameter in lattice
lattice%nclass=nclass
lattice%myclass=>myclass

lattice%analyzed=.true.

contains

  subroutine resize_class()
  implicit none
  integer :: curr_csize
  integer, allocatable :: tmpatom(:,:)
  curr_csize=size(patom,1)
  if(csizenew>curr_csize)then
    !Initially allocate temp array
    allocate(tmpatom(csizenew,nclass))
    !Update size of patom without loosing its content
    tmpatom(1:curr_csize,:)=patom(1:curr_csize,:)
    deallocate(patom); allocate(patom(csizenew,nclass))
    patom=tmpatom
    !Update size of satom without loosing its content
    tmpatom(1:curr_csize,:)=satom(1:curr_csize,:)
    deallocate(satom); allocate(satom(csizenew,nclass))
    satom=tmpatom
    !deallocate temp array
    deallocate(tmpatom)
  endif
  end subroutine resize_class

end subroutine construct_lattice_classes



!----------------------------------------------------------------------------------------
! Create classes of equivalent k-points.
!----------------------------------------------------------------------------------------
subroutine construct_recip_lattice_classes(SymmOp,recip_lattice,applytwist)
use DQMC_RECLATT
integer               :: i,ip,j,isymm,iclass,istart,csize,csizenew,&
                         &jclass,idj,id,mclass,ip2,jstart,&
                         &nsymm,nkpts
integer,allocatable   :: patom(:,:),csizev(:)
integer, pointer      :: myclass_k(:)
integer, pointer      :: map_symm(:,:)
type(symm_operations) :: SymmOp
type(recip_lattice_t) :: recip_lattice
logical,intent(in)    :: applytwist

if(.not.SymmOp%recip_lattice_mapped)stop'Need to map symmetries over lattice before classes (reciprocal)'

!initialize local variables
nsymm=SymmOp%nsymm
nkpts=recip_lattice%nkpts
if(SymmOp%addTimeRev)nsymm=nsymm+1

!Classes in k-space
allocate(csizev(nkpts),myclass_k(nkpts),patom(nkpts,nkpts))

!Each class initially contains one k-point
csizev(:)=1
do ip=1,nkpts
 csizev(ip)=1
 myclass_k(ip)=ip
 patom(1,ip)=ip
enddo

if(applytwist)then
   map_symm=>SymmOp%map_symm_k
else
   map_symm=>SymmOp%map_symm_g
endif

!loop over symmetry operations
do isymm=1,nsymm
 !loop over classes
 do iclass=1,nkpts
  istart=1
  do 
   csize=csizev(iclass)
   csizenew=csize
   do id=istart,csize
    ip2=map_symm(patom(id,iclass),isymm)
    jclass=myclass_k(ip2)
    if(jclass/=iclass)then
     !The two classes need to merged
     jstart=csizenew
     !increase class size for iclass
     csizenew=csizenew+csizev(jclass)
     !transfer jclass elements into iclass
     do idj=1,csizev(jclass)
      ip=patom(idj,jclass)
      myclass_k(ip)=iclass
      patom(jstart+idj,iclass)=ip
     enddo
     !annihilate jclass
     csizev(jclass)=0
    endif
   enddo
   !no new class was merged change iclass
   if(csizenew==csize)exit
   istart=csizev(iclass)+1
   csizev(iclass)=csizenew
  enddo
 enddo
enddo

!relabel classes to exclude classes which were annihilated
mclass=0
do i=1,nkpts
 if(csizev(i)>0)mclass=mclass+1
 do j=1,csizev(i)
  ip=patom(j,i)
  myclass_k(ip)=mclass
 enddo
enddo

!Associate to each class a representative and a size
allocate(recip_lattice%class_size_k(mclass),recip_lattice%class_repr_k(mclass))
j=0
do i=1,nkpts
 if(csizev(i)>0)then
  j=j+1
  recip_lattice%class_repr_k(j)=i
  recip_lattice%class_size_k(j)=csizev(i)
 endif
enddo
deallocate(patom,csizev)

!Save in recip_lattice
recip_lattice%nclass_k=mclass
recip_lattice%myclass_k=>myclass_k

recip_lattice%analyzed=.true.

end subroutine construct_recip_lattice_classes




!---------------------------------------------------------------------------------
! This routines construct my_class_b(ib,jb) where ib and jb are two bonds.
! my_class contains the symmetry class of the pair (ib,jb)
!---------------------------------------------------------------------------------
subroutine construct_bond_classes(Bonds,SymmOp)
use DQMC_BONDS
type(symm_operations),intent(in) :: SymmOp
type(bonds_t),intent(inout) :: Bonds
integer :: ib,nclass,ntotbond,ntotbondsq,isymm,iclass,istart,csize,csizenew,&
& id,bx,by,jclass,jstart,idj,mclass,jb,i,j
integer,pointer :: myclass(:,:) 
integer, allocatable :: bond1(:,:),bond2(:,:),csizev(:)

if(.not.SymmOp%bonds_mapped)stop'Need to map bonds before analyzing symmetry'

ntotbond=size(SymmOp%map_symm_b,1)
ntotbondsq=(ntotbond**2+ntotbond)/2

allocate(myclass(ntotbond,ntotbond),bond1(2,ntotbondsq), &
&        bond2(2,ntotbondsq),csizev(ntotbondsq))

!Initially Define classes as if all bonds were different
!the pair (bond1(ix,iclass) , bond2(ix,iclass)) is the ix-th member
!of class "iclass"
nclass=0
do ib=1,ntotbond
 nclass=nclass+1
 myclass(ib,ib)=nclass
 bond1(1,nclass)=ib
 bond2(1,nclass)=ib
 csizev(nclass)=1
enddo
do ib=1,ntotbond
 do jb=ib+1,ntotbond
  nclass=nclass+1
  myclass(ib,jb)=nclass
  myclass(jb,ib)=nclass
  bond1(1,nclass)=ib
  bond2(1,nclass)=jb
  bond1(2,nclass)=jb
  bond2(2,nclass)=ib
  csizev(nclass)=2
 enddo
enddo
!Try all symmetry operations
do isymm=1,SymmOp%nsymm
 !on all classes
 do iclass=1,nclass
  istart=1
  !we now loop over the elements of a class.
  !This number is increased as we found new equivalent
  !elements. That's why the loop is split in !1! and !2!
  do  !1!
   csize=csizev(iclass)
   csizenew=csize
   do id=istart,csize !2!
    !map the two bonds
    bx=SymmOp%map_symm_b(bond1(id,iclass),isymm)
    by=SymmOp%map_symm_b(bond2(id,iclass),isymm)
    !Find the new class
    jclass=myclass(bx,by)
    if(jclass/=iclass)then
     !Classes are different: merge them
     jstart=csizenew
     !Increase size of iclass
     csizenew=csizenew+csizev(jclass)
     call resize_class()
     !transfer jclass member to iclass
     do idj=1,csizev(jclass)
      bx=bond1(idj,jclass); by=bond2(idj,jclass)
      myclass(bx,by)=iclass
      bond1(jstart+idj,iclass)=bx
      bond2(jstart+idj,iclass)=by
     enddo
     !annihilate jclass
     csizev(jclass)=0
    endif
   enddo
   !Size has not changed so we cannot merge any other class into iclass
   if(csizenew==csize)exit
   !Update loop boundary to find new equivalence due to newly added classes.
   istart=csizev(iclass)+1
   csizev(iclass)=csizenew
  enddo
 enddo
enddo

!Exclude empty classes
mclass=0
do i=1,nclass
 if(csizev(i)>0)mclass=mclass+1
 do j=1,csizev(i)
  bx=bond1(j,i)
  by=bond2(j,i)
  myclass(bx,by)=mclass
 enddo
enddo
deallocate(bond1,bond2,csizev)

!Save classes in Bonds and determine class size
Bonds%nclass_b=mclass
Bonds%myclass_b=>myclass
allocate(Bonds%class_size_b(mclass))
Bonds%class_size_b(:) = 0
do ib=1,ntotbond
 do jb=1,ntotbond
  Bonds%class_size_b(myclass(ib,jb))=Bonds%class_size_b(myclass(ib,jb))+1
 enddo
enddo

Bonds%analyzed=.true.

contains

  subroutine resize_class()
    implicit none
    ! ... Local vars ...
    integer :: curr_csize
    integer, allocatable :: tmpbond(:,:)

    ! ... Executable ...
    curr_csize = size(bond1,1)
    if (csizenew > curr_csize) then
       !Initially allocate temp array
       allocate(tmpbond(csizenew,nclass))
       !Update size of patom without loosing its content
       tmpbond(1:curr_csize,:) = bond1(1:curr_csize,:)
       deallocate(bond1)
       allocate(bond1(csizenew,nclass))
       bond1 = tmpbond
       !Update size of satom without loosing its content
       tmpbond(1:curr_csize,:) = bond2(1:curr_csize,:)
       deallocate(bond2); allocate(bond2(csizenew,nclass))
       bond2=tmpbond
       !deallocate temp array
       deallocate(tmpbond)
    end if
  end subroutine resize_class

end subroutine




!---------------------------------------------------------------------------------
! Construct classes for pairs
!---------------------------------------------------------------------------------
subroutine construct_pair_classes(Bonds,Pairs)
  use DQMC_BONDS
  type(bonds_t),intent(in)    :: Bonds
  type(pairing),intent(inout) :: Pairs
  integer  :: nc, np, ip, jp, ib, jb, ic, jc
  integer, allocatable :: bclass(:)

  nc=0

  np=Pairs%nbond
  allocate(Pairs%myclass_p(np,np),bclass(Bonds%nclass_b))

  !bclass maps bond classes on pair classes. When 0 it means the
  !bond class has yet to be mapped
  bclass=0
  !loop over all pairs
  do ip=1,np
    !Find the corresponding one in Bonds
    ib=Pairs%bond_map(ip)
    do jp=ip,np
      jb=Pairs%bond_map(jp)
      !Find the class for pair (ib,jb)
      ic=Bonds%myclass_b(ib,jb)
      jc=bclass(ic)
      if(jc==0)then
        !we found a new class for pairs
        nc=nc+1
        !Save it
        Pairs%myclass_p(ip,jp)=nc
        Pairs%myclass_p(jp,ip)=nc
        bclass(ic)=nc
      else
        !assign existing class
        Pairs%myclass_p(ip,jp)=jc
        Pairs%myclass_p(jp,ip)=jc
      endif
    enddo
  enddo
  Pairs%nclass_p=nc
  deallocate(bclass)

  !Determine size of pair class
  allocate(Pairs%class_size_p(nc))
  Pairs%class_size_p=0
  do ip=1,np
    do jp=1,np
      jc = Pairs%myclass_p(ip,jp)
      Pairs%class_size_p(jc)=Pairs%class_size_p(jc)+1
    enddo
  enddo
end subroutine



end module DQMC_SYMM
