  subroutine DQMC_change_gtau_time(idir, tau, G_up, G_dn)
    !
    ! Purpose
    ! =======
    ! This subroutine computes a new Gtau which is adjecent 
    ! to the one stored in tau using the recursion relations. 
    ! idir specifies which one of the adjecent four G has to be 
    ! computed. 
    !
    ! tau contains (or may contain) two blocks : 
    ! G(i,j) and G(j,i) where i and j are time indices. tau%ii
    ! and tau%ib are assumed to contain the indices i and j
    ! (Note that somewhere else tau%ib contains the displacement
    ! from tau%ii instead). The variable tau%which says whether
    ! G(i,j) and/or G(j,i) are stored. 
    !
    ! This routine applies the following transformation:
    !   if (idir == 1) G(i,j)=>G(i+1,j) and/or G(j,i)=>G(j,i+1)
    !   if (idir == 2) G(i,j)=>G(i-1,j) and/or G(j,i)=>G(j,i-1)
    !   if (idir == 3) G(i,j)=>G(i,j+1) and/or G(j,i)=>G(j+1,i)
    !   if (idir == 4) G(i,j)=>G(i,j-1) and/or G(j,i)=>G(j-1,i)
    ! keeping correctly track of the case where i==j (either 
    ! initially or after the transformation). 
    ! i and j are always kept between 1 and L.
    !
    ! Arguments
    ! =========
    !
    type(Gtau),  intent(inout) :: tau
    type(G_fun), intent(in)    :: G_up, G_dn
    integer,     intent(in)    :: idir
  
    ! ... local ...
    integer :: i, j, id, n, L

    ! ... aliases ...
    real(wp), pointer :: up(:,:)          
    real(wp), pointer :: dn(:,:)
    real(wp), pointer :: W(:,:)
    type(MatB), pointer :: B_up, B_dn

    B_up => tau%B_up
    B_dn => tau%B_dn
    W => tau%W1

    n = tau%n
    L = tau%L

    if(tau%which .le. TAU_BOTH) then

       up => tau%upt0
       dn => tau%dnt0

       select case (idir)

       case (1) ! G(i,j)=> G(i+1,j) 

          i = tau%ii + 1
          if(i > L) i = 1

          !Multiply by B_{ii+1} 
          call DQMC_MultB_Left  (n, up, B_up, G_up%V(:,i), W)
          call DQMC_MultB_Left  (n, dn, B_dn, G_dn%V(:,i), W)

          !Time wrapped through beta. Need to change sign.
          if(i == 1)then
            up = -up
            dn = -dn
          endif

          !Final G is equal time. Handle G(i,j) properly.
          if (tau%ib == i) then
             do id = 1, n
                up(id,id) = 1.d0 + up(id,id)
                dn(id,id) = 1.d0 + dn(id,id)
             enddo
          endif

       case (2) ! G(i,j)=> G(i-1,j) 

          i = tau%ii

          !Initial G is equal time. Handle G(i,j) properly.
          if (tau%ib == tau%ii) then
             do id = 1, n
                up(id,id) = -1.d0 + up(id,id)
                dn(id,id) = -1.d0 + dn(id,id)
             enddo
          endif

          call DQMC_MultBi_Left(n, up, B_up, G_up%V(:,i), W)
          call DQMC_MultBi_Left(n, dn, B_dn, G_dn%V(:,i), W)

          !Time wrapped through zero. Need to change sign.
          if(i == 1)then
            up = -up
            dn = -dn
            i = L
          else
            i = i -1
          endif

       case (3) !G(i,j)=> G(i,j+1) 

          j = tau%ib + 1
          if (j > L) j = 1

          if (tau%ib == tau%ii) then
             do id = 1, n
                up(id,id) = -1.d0 + up(id,id)
                dn(id,id) = -1.d0 + dn(id,id)
             enddo
          endif

          call DQMC_MultBi_Right(n, up, B_up, G_up%V(:,j), W)
          call DQMC_MultBi_Right(n, dn, B_dn, G_dn%V(:,j), W)

          !Time wrapped through beta. Need to change sign.
          if(j == 1)then
            up = -up
            dn = -dn
          endif

       case(4) !G(i,j)=> G(i,j-1) 

          j = tau%ib

          call DQMC_MultB_Right(n, up, B_up, G_up%V(:,j), W)
          call DQMC_MultB_Right(n, dn, B_dn, G_dn%V(:,j), W)

          !Time wrapped through zero. Need to change sign.
          if(j == 1)then
            up = -up
            dn = -dn
            j = L
          else
            j = j - 1
          endif

          !Final G is equal time. Treat G(i,j) properly.
          if(tau%ii == j)then
             do id = 1, n
                up(id,id) = 1.d0 + up(id,id)
                dn(id,id) = 1.d0 + dn(id,id)
             enddo
          endif

       end select

    endif

    if(tau%which .ge. TAU_BOTH) then

       up => tau%up0t
       dn => tau%dn0t

       select case (idir)

       case (1) ! G(j,i)=>G(j,i+1)

          i = tau%ii+1
          if(i > L) i = 1

          !initial G is equal time. Handle G(j,i) properly.
          if (tau%ib == tau%ii) then
             do id = 1, n
                up(id,id) = -1.d0 + up(id,id)
                dn(id,id) = -1.d0 + dn(id,id)
             enddo
          endif

          !Multiply by B_{i+1} and its inverse
          call DQMC_MultBi_Right(n, up, B_up, G_up%V(:,i), W)
          call DQMC_MultBi_Right(n, dn, B_dn, G_dn%V(:,i), W)

          !Time wrapped through beta. Need to change sign.
          if(i == 1)then
            up = -up
            dn = -dn
          endif

       case(2) ! G(j,i)=>G(j,i-1)

          i = tau%ii

          call DQMC_MultB_Right(n, up, B_up, G_up%V(:,i), W)
          call DQMC_MultB_Right(n, dn, B_dn, G_dn%V(:,i), W)

          !Time wrapped through zero. Need to change sign.
          if(i == 1)then
            up = -up
            dn = -dn
            i = L
          else
            i = i -1
          endif

          !Final G is equal time. Handle G(j,i) properly.
          if(tau%ib == i)then
             do id=1,n
                up(id,id) = 1.d0 + up(id,id)
                dn(id,id) = 1.d0 + dn(id,id)
             enddo
          endif

       case(3) !G(j,i)=>G(j+1,i)

          j = tau%ib + 1
          if(j > L) j = 1

          call DQMC_MultB_Left  (n, up, B_up, G_up%V(:,j), W)
          call DQMC_MultB_Left  (n, dn, B_dn, G_dn%V(:,j), W)

          !Time wrapped through beta. Need to change sign.
          if(j == 1)then
            up = -up
            dn = -dn
          endif

          !Final G is equal time. Handle G(j,i) properly.
          if(tau%ii == j)then
             do id=1,n
                up(id,id) = 1.d0 + up(id,id)
                dn(id,id) = 1.d0 + dn(id,id)
             enddo
          endif

       case(4) ! G(j,i)=>G(j-1,i)

          j = tau%ib

          if (tau%ii == tau%ib) then
             do id = 1, n
                up(id,id) = -1.d0 + up(id,id)
                dn(id,id) = -1.d0 + dn(id,id)
             enddo
          endif

          call DQMC_MultBi_Left(n, up, B_up, G_up%V(:,j), W)
          call DQMC_MultBi_Left(n, dn, B_dn, G_dn%V(:,j), W)

          !Time wrapped through zero. Need to change sign.
          if(j == 1)then
             up = -up
             dn = -dn
             j = L
          else
             j = j - 1
          endif

       end select

    endif

    !Update block index
    select case (idir)
      case(1, 2)
         tau%ii = i
      case(3, 4)
         tau%ib = j
    end select 

  end subroutine DQMC_change_gtau_time
