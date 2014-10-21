module DQMC_KBONDS

use dqmc_reclatt
implicit none

type :: kbonds_t

   integer :: nak
   integer :: nbonds
   integer :: nmomenta

   integer,  pointer  :: map_symm_ak(:,:)    
   integer,  pointer  :: map_symm_bak(:,:,:) 
   integer,  pointer  :: class_size(:,:)     
   integer,  pointer  :: myclass(:,:,:) 
   integer,  pointer  :: nclass(:) 
   integer,  pointer  :: bond_origin(:, :) 
   integer,  pointer  :: bond_target(:, :) 
   integer,  pointer  :: bmap(:, :, :)    

   real(wp), pointer  :: ksum(:,:) 

end type

contains

!-----------------------------------------------------------------!

 subroutine init_kbonds(symm, lattice, reclatt, kbonds)

   use dqmc_symm
   use dqmc_latt

   ! This subroutine defines the action of symmetry operation on
   ! the pair (iat, k) where k is one of the k-points and iat is
   ! one of the orbital in the primitive cell. To each pair (iat,k)
   ! it is associated the integer iatk = (iat-1)*nk + k. The routine
   ! constructs the matrix map_symm_ak(iatk, isymm) that given the
   ! symmetry operation isymm returns the pair (jat, k') where (iat,k)
   ! is mapped.

   type(symm_operations), intent(in)    :: symm
   type(lattice_t),       intent(in)    :: lattice
   type(recip_lattice_t), intent(in)    :: reclatt
   type(kbonds_t),        intent(inout) :: kbonds
   integer :: iak, ia, ik, jk, ja, jak, na ,nk, ns, is, nm, nd, id
   character(len=1)  :: snd
   character(len=50) :: sfmt

   na = lattice%natom
   nk = reclatt%nkpts
   ns = Symm%nsymm
   nm = reclatt%nmomenta
   nd = reclatt%ndim

   kbonds%nak    = na * nk
   kbonds%nbonds = na * na * nk
   kbonds%nmomenta = nm
   allocate(kbonds%map_symm_ak(na*nk, ns))
   
   !Loop over atom types
   do ia = 0, na-1
      !Loop over k-points
      do ik = 1, nk
         !Construct state label
         iak = ia * nk + ik
         !Act with symmetry operation
         do is = 1, ns
            !Act on k-point
            jk = symm%map_symm_k(ik, is)
            !Act on site
            ja = symm%map_symm(ia, is)
            !Check this is an allowed mapping
            if(ja >= na) then
               write(*,*) 'Only symmetry operations mapping sites   &
               &  in the primitive cell to sites in the             &
               &  primitive cell are allowed.'
               stop
            endif
            !Get the new state
            jak = ja * nk + jk
            !Store the mapping
            kbonds%map_symm_ak(iak, is) = jak
         enddo
      enddo
   enddo

   write(snd,'(i1)')nd
   write(sfmt,'(5A)')'(i2,1x,',snd,'f9.5,1x,"->",1x,i2,1x,',snd,'f9.5)'
    
   open(unit=25, file='symmetry_map.info', status='unknown',position='append')
   do is = 1, ns
      id = symm%valid_symm(is)
      write(25,*)'Mapping of Symmetry :',id, symm%Symmlabel(id)
      do ia = 0, na-1
         do ik = 1, nk
            iak = ia * nk + ik
            jak = kbonds%map_symm_ak(iak, is)
            jk = mod(jak, nk)
            if(jk==0)jk=nk
            ja = (jak - jk)/nk
            write(25,sfmt) ia, (reclatt%klist(ik,id), id=1,nd),&
            &              ja, (reclatt%klist(jk,id), id=1,nd)
         enddo
      enddo
   enddo
   close(25)


   !Store total pair momenta in kbonds (handy).
   allocate(kbonds%ksum(nm, rdim))
   kbonds%ksum = reclatt%ksum

 end subroutine init_kbonds

!-----------------------------------------------------------------!

 subroutine construct_kbonds(reclatt, kbonds)

   ! Given the total number of momenta for the pair,
   ! construct bonds amongst states (iat,k) and (jat,k') 
   ! that have a total momentum k+k' that sums up to the 
   ! requested one. Note that a bond has a direction:
   ! it goes from (iat,k), where an up spin particle is
   ! created, to (jat,k), where a down one is created. 
   ! It is possible for some of the pairs to be identical
   ! (overcounting). This, however, does not pose any major
   ! problem apart from a small overhead during the computation
   ! of the 2-particle GF.

   type(recip_lattice_t), intent(in) :: reclatt
   type(kbonds_t), intent(inout)     :: kbonds

   integer :: nak, nmom, nb, nat, nk
   integer :: im, iak, jat, ik, jk, jak, ib
   integer :: nbond_from(kbonds%nak)
   logical :: assigned(kbonds%nak, kbonds%nak)

   nak    = kbonds%nak
   nb     = kbonds%nbonds
   nmom   = reclatt%nmomenta
   nk     = reclatt%nkpts
   nat    = nak / nk
   
 
   allocate(kbonds%bmap(nak,nak,nmom))
   allocate(kbonds%bond_origin(nb,nmom))
   allocate(kbonds%bond_target(nb,nmom))

   kbonds%bmap = -1

   !Loop over pair total momenta
   do im = 1, nmom
      ib = 0
      assigned = .false.
      nbond_from = 0
      !Loop over state
      do iak = 1, nak
         !Find k-point for the state
         ik = mod(iak-1, nk) + 1
         !and its mate jk i.e. the one that summed to ik gives
         ! the im-th total momentum.
         jk = reclatt%kmate(ik, im)
         !Loop over atom types
         do jat = 0, nat-1
            !and construct all states with momentum jk
            jak =  jat * nk + jk
            !assign a number to the pair iak, jak
            if(.not.assigned(iak,jak))then
               ib = ib + 1
               nbond_from(iak) = nbond_from(iak) + 1
               kbonds%bond_origin(ib, im) = iak
               kbonds%bond_target(ib, im) = jak
               kbonds%bmap(iak, jak, im)  = ib
               assigned(iak, jak)         = .true.
            endif
            !and to the pair jak, iak
            if(.not.assigned(jak,iak))then
               ib = ib + 1
               nbond_from(jak) = nbond_from(jak) + 1
               kbonds%bond_origin(ib, im) = jak
               kbonds%bond_target(ib, im) = iak
               kbonds%bmap(jak, iak, im)  = ib
               assigned(jak, iak)         = .true.
            endif
         enddo
      enddo
      !Check the number of pair equals the number of states
      if (ib /= nb) then
         write(*,*)'Pair momentum: ',(kbonds%ksum(im, ik), ik =1,3)
         write(*,*)'nb /= ib. Suspicious. Stop.'
         stop
      endif
   enddo
 
 end subroutine construct_kbonds

!-----------------------------------------------------------------!

 subroutine map_symm_kbond(kbonds)
 
 ! Given a pair (iat,k),(jat,k') this routines set up a matrix
 ! that map the action of all symmetry operations into the
 ! label of the final pair.
 
   type(kbonds_t), intent(inout) :: kbonds
   integer :: im, ib, jb, iak, jak, isymm, inew, jnew, nsymm, nb, nm
 
   nsymm = size(kbonds%map_symm_ak, 2)
   nb = kbonds%nbonds
   nm = kbonds%nmomenta
   allocate(kbonds%map_symm_bak(nb, nsymm, nm))
   
   !Loop over all pair momenta
   do im = 1, kbonds%nmomenta
   !Loop over all pairs with momenta im
      do ib = 1, nb
         !Get the two states
         iak = kbonds%bond_origin(ib, im)
         jak = kbonds%bond_target(ib, im)
         !Loop over all symmetries
         do isymm = 1, nsymm
            !Get the two new states
            inew = kbonds%map_symm_ak(iak, isymm)
            jnew = kbonds%map_symm_ak(jak, isymm)
            !and their pair number
            jb   = kbonds%bmap(inew, jnew, im)
            if (jb < 0) stop'Symmetry map is wrong'
            !Store the mapping
            kbonds%map_symm_bak(ib, isymm, im) = jb
         enddo
      enddo
   enddo
 
 end subroutine map_symm_kbond

!-----------------------------------------------------------------!

 subroutine construct_kbond_classes(kbonds)
 
 ! This routine construct my_class_b(ib,jb) where ib and jb are two bonds.
 ! my_class contains the symmetry class of the pair (ib,jb)
 
 
   type(kbonds_t),intent(inout)      :: kbonds
   integer :: nclass, ntotbond
   integer :: ib, isymm, iclass, istart, im, ak1 ,ak2
   integer :: csize, csizenew, maxclass, nsymm
   integer :: id, bx, by, jclass, jstart, idj, mclass, jb, i, j
   integer,pointer :: myclass(:,:) 
   integer, allocatable :: bond1(:,:), bond2(:,:), csizev(:)
   integer, pointer :: map_symm_b(:,:) 
   
   ntotbond   = kbonds%nbonds
   nsymm      = size(kbonds%map_symm_ak, 2)
   allocate(myclass(ntotbond,ntotbond))
   allocate(kbonds%myclass(ntotbond,ntotbond,kbonds%nmomenta))
   allocate(kbonds%nclass(kbonds%nmomenta))
      
   do im = 1, kbonds%nmomenta
   
      nclass = (ntotbond**2+ntotbond) / 2
      allocate(csizev(nclass))
      allocate(bond1(2,nclass))
      allocate(bond2(2,nclass))
      
      map_symm_b => kbonds%map_symm_bak(:,:,im)
      
      !Initially Define classes as if all bonds were different
      !the pair (bond1(ix,iclass) , bond2(ix,iclass)) is the ix-th member
      !of class "iclass"
      nclass = 0

      !Classes made up by a bond and itself
      do ib = 1, ntotbond
         nclass = nclass + 1
         myclass(ib,ib)  = nclass
         bond1(1,nclass) = ib
         bond2(1,nclass) = ib
         csizev(nclass)  = 1
      enddo
   
      !Classes made up by distinct bonds
      do ib = 1, ntotbond
         do jb = ib+1, ntotbond
            nclass = nclass + 1
            myclass(ib,jb)  = nclass
            myclass(jb,ib)  = nclass
            bond1(1,nclass) = ib
            bond2(1,nclass) = jb
            bond1(2,nclass) = jb
            bond2(2,nclass) = ib
            csizev(nclass)  = 2
         enddo
      enddo
      
      !Try all symmetry operations. The "+1" operation corresponds to
      !the symmetry between up-dn-dn-up and dn-up-up-dn.
      do isymm = 1, nsymm + 1
         !on all classes
         do iclass = 1, nclass
            istart = 1
            !we now loop over the elements of a class.
            !This number is increased as we found new equivalent
            !elements. 
            do  
               csize    = csizev(iclass)
               csizenew = csize
               do id = istart, csize 
                  !map the two bonds
                  if (isymm == nsymm + 1) then
                     ak1 = kbonds%bond_origin(bond1(id,iclass),im)
                     ak2 = kbonds%bond_target(bond1(id,iclass),im)
                     bx  = kbonds%bmap(ak2, ak1, im)
                     ak1 = kbonds%bond_origin(bond2(id,iclass),im)
                     ak2 = kbonds%bond_target(bond2(id,iclass),im)
                     by  = kbonds%bmap(ak2, ak1, im)
                  else
                     bx  = map_symm_b(bond1(id,iclass),isymm)
                     by  = map_symm_b(bond2(id,iclass),isymm)
                  endif
                  !Find the new class
                  jclass = myclass(bx,by)
                  if (jclass /= iclass) then
                     !Classes are different: merge them
                     jstart = csizenew
                     !Increase size of iclass
                     csizenew = csizenew + csizev(jclass)
                     call resize_class()
                     !transfer jclass member to iclass
                     do idj=1,csizev(jclass)
                        bx = bond1(idj, jclass) 
                        by = bond2(idj, jclass)
                        myclass(bx,by) = iclass
                        bond1(jstart+idj,iclass) = bx
                        bond2(jstart+idj,iclass) = by
                     enddo
                     !annihilate jclass
                     csizev(jclass) = 0
                  endif
               enddo
               !Size has not changed so we cannot merge any other class into iclass
               if (csizenew == csize) exit
               !Update loop boundary to find new equivalence due to newly added classes.
               istart = csizev(iclass) + 1
               csizev(iclass) = csizenew
            enddo
         enddo
      enddo
      
      !Exclude empty classes
      mclass = 0
      do i = 1, nclass
         if( csizev(i) > 0 ) mclass = mclass + 1
         do j = 1, csizev(i)
            bx = bond1(j,i)
            by = bond2(j,i)
            myclass(bx,by) = mclass
         enddo
      enddo
      
      !Save classes in Bonds and determine class size
      kBonds%nclass(im) = mclass
      kBonds%myclass(:,:,im) = myclass(:,:)

      deallocate(bond1, bond2, csizev)
   
   enddo ! npairk

   deallocate(myclass)

   maxclass = maxval(kbonds%nclass)
   allocate(kbonds%class_size(maxclass, kbonds%nmomenta))
   kbonds%class_size = 0
   do im = 1, kbonds%nmomenta
      do ib = 1, ntotbond
         do jb = 1, ntotbond
            id = kbonds%myclass(ib, jb, im)
            kbonds%class_size(id, im) = kbonds%class_size(id, im) + 1
         enddo
      enddo
   enddo

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
   
   end subroutine construct_kbond_classes

end module DQMC_KBONDS
