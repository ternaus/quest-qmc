module DQMC_BONDS

use DQMC_geom_param
use DQMC_LATT

implicit none

type :: bonds_t
 integer               :: ntotbond         !number of bonds read from input, possibly augmented by symmetry
 integer, pointer      :: bond_label(:)    !label of bond (ntotbond)
 integer, pointer      :: bond_origin(:)   !primitive cell site from which bond originates (ntotbond)
 integer, pointer      :: bond_target(:) 

 real(wp), pointer       :: xxbond(:,:)    !bond vector cartesian coordinates (rdim,ntotbond)

 integer               :: nclass_b         !number of inequivalent bond pairs
 integer, pointer      :: myclass_b(:,:)   !class for pair of bonds (ntotbond,ntotbond)
 integer, pointer      :: class_size_b(:)  !number of equivalent pairs of bond in each class (nclass_b)

 logical               :: initialized
 logical               :: analyzed

end type


type :: pairing
 integer               :: nwave                 !number of waves to analyze
 integer               :: nbond                 !number of bonds necessary to describe  waves
 integer, pointer      :: nbondv(:)             !number of bonds in each cell (0:ncell-1)(constant if no-dilution)
 integer, pointer      :: bond_origin(:,:)      !origin of bond in a given cell (nbond, 0:ncell-1)
 integer, pointer      :: bond_end(:,:)         !end of bond originating from cell (nbond, 0:ncell-1)
 integer, pointer      :: bond_map(:)           !mapping to the bonds in bonds_t (nbond)
 integer, pointer      :: pair_map(:)           !inverse of bond_map: Given a bond returns its label in pairing. 
                                                !If 0 bond does not enter pairing
 integer, pointer      :: bond_number(:,:)      !label of bonds starting from each cell (nbond,0:ncell-1)(needed when diluted)

 real(wp), pointer     :: bond_wgt(:,:)         !coefficient forming the wave (nwave, nbond)

 character*20, pointer :: wave_label(:)         !label of wave: s-wave, d-wave etc... (nwave)

 integer, pointer      :: myclass_p(:,:) 
 integer               :: nclass_p
 integer, pointer      :: class_size_p(:) 

 logical               :: initialized
end type


contains

!------------------------------------------------------!

subroutine free_bonds(bonds)
  
 type(bonds_t), intent(out) :: bonds

 if(associated(bonds%xxbond))       deallocate(bonds%xxbond)         
 if(associated(bonds%myclass_b))    deallocate(bonds%myclass_b)
 if(associated(bonds%bond_label))   deallocate(bonds%bond_label)
 if(associated(bonds%bond_origin))  deallocate(bonds%bond_origin)
 if(associated(bonds%bond_target))  deallocate(bonds%bond_target)
 if(associated(bonds%class_size_b)) deallocate(bonds%class_size_b)

end subroutine free_bonds

!------------------------------------------------------!

subroutine free_pairs(pairs)

 type(pairing), intent(out) :: pairs

 if(associated(pairs%nbondv))         deallocate(pairs%nbondv)
 if(associated(pairs%bond_end))       deallocate(pairs%bond_end)
 if(associated(pairs%bond_map))       deallocate(pairs%bond_map)
 if(associated(pairs%pair_map))       deallocate(pairs%pair_map)
 if(associated(pairs%bond_wgt))       deallocate(pairs%bond_wgt)
 if(associated(pairs%myclass_p))      deallocate(pairs%myclass_p)
 if(associated(pairs%wave_label))     deallocate(pairs%wave_label)
 if(associated(pairs%bond_origin))    deallocate(pairs%bond_origin)
 if(associated(pairs%bond_number))    deallocate(pairs%bond_number)
 if(associated(pairs%class_size_p))   deallocate(pairs%class_size_p)

end subroutine free_pairs


!-------------------------------------------------------------------------------------
! This subroutine defines 1)bond_label(ibond): the label of the ibond-th
! bond. This label is just the number of the record corresponding to that bond.
! Bonds between sites A and B (A/=B) have two labels: irec for A->B and -irec for
! B->A. 2)bond_origin(ibond): the site where the ibond-th bond originates.
! xxbond(1:3,ibond) are the cartesian coordinate describing the bond. When
! bond_label(jbond)=-bond_label(ibond) => xxbond(1:3,jbond)=-xxbond(1:3,ibond)
!-------------------------------------------------------------------------------------
subroutine read_bonds(Bonds, SOP)
type(bonds_t) :: Bonds
integer, intent(in) :: SOP
integer :: nrecord,iat,jat,ibond,irec,ios,i,ntotbond
integer, pointer :: bond_label(:) 
integer, pointer :: bond_origin(:) 
integer, pointer :: bond_target(:) 
real(wp) :: delta(rdim)
real(wp),pointer :: xxbond(:,:) 
character*100 string
logical ldum

ntotbond=0
nrecord=0
Bonds%initialized=.false.

if(move_to_record(INPUT_FIELDS(BONDS_F),inpunit))then

 !Count the bonds
 do 
  read(inpunit,'(A)')string
  read(string,*,iostat=ios)iat,jat,(delta(i),i=1,rdim)
  if(ios.ne.0)exit
  nrecord=nrecord+1
  ntotbond=ntotbond+1
  !Add the bond the goes in the opposite direction to the list
  if(.not.(iat==jat.and.sum(delta**2)<1.d-10))ntotbond=ntotbond+1
 enddo 

 allocate(bond_origin(ntotbond),bond_target(ntotbond),bond_label(ntotbond),xxbond(rdim,ntotbond))

 !rewind to the beginning of BONDS field
 ldum=move_to_record(INPUT_FIELDS(BONDS_F),inpunit)

 ibond=0
 !Read again and store bond information
 do irec=1,nrecord
  read(inpunit,'(A)')string
  read(string,*)iat,jat,(delta(i),i=1,rdim)
  ibond=ibond+1
  xxbond(:,ibond)=delta(:)
  bond_label(ibond)=irec
  bond_origin(ibond)=iat
  bond_target(ibond)=jat
  if(iat==jat.and.sum(delta**2)<1.d-10)cycle
  !Store info for "opposite" bond
  ibond=ibond+1
  xxbond(:,ibond)=-delta(:)
  bond_label(ibond)=-irec
  bond_origin(ibond)=jat
  bond_target(ibond)=iat
 enddo 

 !Print out
 write(SOP,*)'Bonds (from input)'
 do ibond=1,ntotbond
  write(SOP,'(3i4,3f12.7)')ibond,bond_label(ibond),bond_origin(ibond),xxbond(1:rdim,ibond)
 enddo

 !Save info in Bonds
 Bonds%ntotbond=ntotbond
 Bonds%bond_label=>bond_label
 Bonds%bond_origin=>bond_origin
 Bonds%bond_target=>bond_target
 Bonds%xxbond=>xxbond

 Bonds%initialized=.true.
 Found_Field(BONDS_F)=.true.

endif
end subroutine read_bonds




!------------------------------------------------------------------------------------------
! Construct all bonds i.e. pairs of atoms that enter pairing susceptibility.
! Each bond is characterized by pair_bond_origin(ibond,icell) which specify the
! site on cell icell where the ibond-th bond begins and pair_bond_end(ibond,icell) that
! specify the site on which it ends. Reads and store the names (wave_label) and 
! coeffiecients (pair_bond_wgt) with which the pairing susceptibility have to be
! combined to obtain pair state with a given angular momentum.
!------------------------------------------------------------------------------------------
subroutine construct_pairs(Bonds,Pairs,lattice, SOP)
type(bonds_t) :: Bonds
type(pairing) :: Pairs
type(lattice_t) :: lattice
integer, intent(in) :: SOP
integer :: ios,idum,ibond,iwave,iat,jat,icell,isite,jbond,nwave,npairbond,natom,nsites,ncell
integer, allocatable :: pair_bond_label(:)
integer, pointer     :: pair_bond_origin(:,:) 
integer, pointer     :: pair_bond_end(:,:) 
integer, pointer     :: pair_bond_map(:) 
integer, pointer     :: npairbondv(:) 
integer, pointer     :: pair_bond_number(:,:) 
real(wp)  :: rdum
real(wp), pointer     :: pair_bond_wgt(:,:) 
character*100 string
character*10 :: label
character*20,pointer :: wave_label(:) 
logical :: ldum, donullify

natom=lattice%natom
nsites=lattice%nsites
ncell=lattice%ncell
npairbond=0

if(.not.Bonds%initialized)stop'Need to initialize Bonds before calling make_pairs'
allocate(Pairs%pair_map(Bonds%ntotbond))

if(move_to_record(INPUT_FIELDS(PAIRS_F),inpunit))then

   read(inpunit,'(A)')string
   !count how many bonds are used to define pairing function
   do 
      if(npairbond>Bonds%ntotbond)stop 'Too many bonds in #PAIRING.'
      npairbond=npairbond+1
      read(string,*,iostat=ios)(idum,ibond=1,npairbond)
      if(ios/=0)exit
   enddo
   npairbond=npairbond-1
  
   !Count "waves" 
   nwave=0
   do 
      read(inpunit,*,iostat=ios)label,(rdum,ibond=1,npairbond)
      if(ios/=0)exit
      nwave=nwave+1
   enddo
  
   allocate(wave_label(nwave),pair_bond_wgt(nwave,npairbond),pair_bond_map(npairbond),pair_bond_label(npairbond))
   !rewind to beginning of #PAIR field
   ldum=move_to_record('#PAIR',inpunit)
   !Read and store the bonds defining the pairs
   read(inpunit,*)(pair_bond_label(ibond), ibond=1,npairbond)
   !Construct the mapping from pairs to bond (bond_map) and viceversa (pair_map)
   Pairs%pair_map(:) = 0
   do ibond = 1, npairbond
      do jbond = 1, Bonds%ntotbond
         if(pair_bond_label(ibond) == Bonds%bond_label(jbond))then
            pair_bond_map(ibond)  = jbond
            Pairs%pair_map(jbond) = ibond
            exit
         endif
      enddo
      if(jbond>Bonds%ntotbond)stop 'One or more labels in #PAIRING are wrong' 
   enddo
  
   do iwave = 1, nwave
      read(inpunit,*)wave_label(iwave),(pair_bond_wgt(iwave,ibond),ibond=1,npairbond)
      rdum=sqrt(sum(pair_bond_wgt(iwave,1:npairbond)**2))
      pair_bond_wgt(iwave,1:npairbond)=pair_bond_wgt(iwave,1:npairbond)/rdum
   enddo
  
   !Print
   write(SOP,'(1x,A)')'Wave coefficients'
   write(SOP,'(11x,20(5x,i2,3x))')(pair_bond_map(ibond),ibond=1,npairbond)
   do iwave = 1, nwave
      write(SOP,'(1x,A10,20f10.5)')wave_label(iwave),(pair_bond_wgt(iwave,ibond),ibond=1,npairbond)
   enddo
   write(SOP,'(76(''=''))')
   deallocate(pair_bond_label)
   donullify = .false.
else
   !The field #PAIR was not specified: Use as many waves as bonds
   npairbond = Bonds%ntotbond
   nwave = npairbond
   allocate(wave_label(nwave),pair_bond_map(npairbond))
   !label is temporarily given by an integer
   do iwave = 1, nwave
      write(wave_label(iwave),'(i2)')iwave
   enddo
   !mapping is trivial
   do ibond = 1, npairbond
      pair_bond_map(ibond)  = ibond
      Pairs%pair_map(ibond) = ibond
   enddo
   donullify = .true.
endif

!Define origin and target of all bonds in simulation cell.
allocate(pair_bond_origin(npairbond,0:ncell-1),pair_bond_end(npairbond,0:ncell-1),npairbondv(0:ncell-1), &
&  pair_bond_number(npairbond,0:ncell-1))
!store how many bonds are associated with a given cell
npairbondv(0:ncell-1) = npairbond !needed when system is site-diluted
!loop over all sites
do isite = 0, nsites-1
   !find atom type and cell where isite is contained
   iat = mod(isite,natom)
   icell = isite/natom
   do ibond = 1,npairbond
      jbond = pair_bond_map(ibond)
      !Check which bonds originate from that site
      if(Bonds%bond_origin(jbond) == iat)then
         !Save isite as the site in cell "icell" from where "ibond" originates
         pair_bond_origin(ibond,icell) = isite
         jat = Bonds%bond_target(jbond)
         !Save the site where ibond from icell is pointing
         pair_bond_end(ibond,icell) = hoptowho(isite,Bonds%xxbond(1:rdim,jbond),jat,lattice)
         pair_bond_number(ibond,icell) = ibond !needed when system is site-diluted
      endif
   enddo
enddo

!Fill members of Pairs
Pairs%nwave = nwave
Pairs%nbond = npairbond
Pairs%nbondv      =>  npairbondv
Pairs%bond_map    =>  pair_bond_map
Pairs%bond_number =>  pair_bond_number
Pairs%bond_end    =>  pair_bond_end
Pairs%bond_origin =>  pair_bond_origin
Pairs%wave_label  =>  wave_label

if(donullify)then
  nullify(Pairs%bond_wgt)
else
  Pairs%bond_wgt => pair_bond_wgt
endif

end subroutine construct_pairs


end module DQMC_BONDS

