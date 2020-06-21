c///////////////////////////////////////////////
c calculations.f - Radius of gyration, stresses,
c                  structure factors, etc.
c 
c Mike Hore, Hilton Beach, ON, July 2005
c//////////////////////////////////////////////

c
c bf_finalize_energy() - Finalizes the energy calculation by
c                     collecting the results at the master
c                     level.
c
      subroutine bf_finalize_energy()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      double precision utot_master

      call MPI_REDUCE(utot, utot_master, 1, MPI_DOUBLE_PRECISION,
     >                MPI_SUM, MASTER, MPI_COMM_WORLD, mpi_err)

      if (cluster_rank .eq. MASTER) then
         open (12, file="energy.data", access="append")
         write(12,*) itime, utot_master
         close(12)
      endif

      return
      end

c
c bf_find_clusters() - This is done at the master level
c                   whenever the checkpoint is written...
c                   this saves time in communication.
c
c (Mohamed Laradji)
      subroutine bf_find_clusters()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      double precision r_time1, r_time2
      double precision x_trans, y_trans, z_trans
  
      integer ix, iy, iz, jx, jy, jz, mxy, idx_i, idx_j
      integer icell, i, j, k, ii, jj, jcell, kcell
      integer iti, itj
      integer li,li0,lj,lj0
      integer i1,i2,j1,j2
      
      if (DISABLE_CLUSTERS .eq. 1) return
      if (cluster_rank .ne. MASTER) return

      ! Use the information in DBL_BUF to find clusters.
      r_time1 = MPI_WTIME()

      mxy = mx*my

      ! Clear out old data from list, head.
      do i=1, ncell
         head(i) = 0
      enddo

      do i=1, particle_count
         list(i) = 0
      enddo

      do i=1, particle_count
          clust_AA(i)=0
          size_AA(i)=0
          head_AA(i)=0
          list_AA(i)=0
          clust_BB(i)=0
          size_BB(i)=0
          head_BB(i)=0
          list_BB(i)=0
      enddo
      nbr_clust_AA=0
      nbr_clust_BB=0
      cl_AA=0
      cl_BB=0
 
      ! Consider new particle positions:
      do i=1, particle_count
         j = (i-1)*8
         
         if (DBL_BUF1(j+1) .le. 0.d0) then
             DBL_BUF1(j+1) = DBL_BUF1(j+1) + dim_x
         else if (DBL_BUF1(j+1) .gt. dim_x) then 
             DBL_BUF1(j+1) = DBL_BUF1(j+1) - dim_x
         endif

         if (DBL_BUF1(j+2) .le. 0.d0) then
             DBL_BUF1(j+2) = DBL_BUF1(j+2) + dim_y
         else if (DBL_BUF1(j+2) .gt. dim_y) then 
             DBL_BUF1(j+2) = DBL_BUF1(j+2) - dim_y
         endif

         if (DBL_BUF1(j+3) .le. 0.d0) then
             DBL_BUF1(j+3) = DBL_BUF1(j+3) + dim_z
         else if (DBL_BUF1(j+3) .gt. dim_z) then 
             DBL_BUF1(j+3) = DBL_BUF1(j+3) - dim_z
         endif

         icell = 1 + int(DBL_BUF1(j+1)*cellx_inv)    + 
     >               int(DBL_BUF1(j+2)*celly_inv)*mx +
     >               int(DBL_BUF1(j+3)*cellz_inv)*mxy

         ! Update table.
         j           = head(icell)
         head(icell) = i
         list(i)     = j
      enddo

      do iz=1, mz
        do iy=1, my
          do ix=1, mx
       
            ! Which cell are we in?
            icell = (iz-1)*mxy+(iy-1)*mx+ix

            ! Get head of this cell
            ii = head(icell)
              
            ! If i not equal to 0, the cell is not empty.
            if (ii .gt. 0) then
               do kcell=1, 14

                  i = ii ! Used in looping, so is not redundant

                  x_trans = 0.d0
                  y_trans = 0.d0
                  z_trans = 0.d0

                  ! Neighboring cell number:
                  jx = ix + nix(kcell)
                  jy = iy + niy(kcell)
                  jz = iz + niz(kcell)
                  
                  ! Periodic boundary:
                  if (ix .eq. mx .and. jx .gt. ix) then
                     jx = 1
                     x_trans = -dim_x
                  elseif (ix .eq. 1 .and. jx .lt. ix) then
                     jx = mx
                     x_trans = dim_x
                  endif

                  if (iy .eq. my .and. jy .gt. iy) then
                     jy = 1
                     y_trans = -dim_y
                  elseif (iy .eq. 1 .and. jy .lt. iy) then
                     jy = my
                     y_trans = dim_y 
                  endif

                  if (iz .eq. mz .and. jz .gt. iz) then
                     jz = 1 
                     z_trans = -dim_z
                  endif ! end PBC checks
                 
                  ! Now, calculate neighboring cell index:
                  jcell = (jz-1)*mxy+(jy-1)*mx+jx
              
                  ! Get head of neighboring cell:
                  j = head(jcell)
                  
                  ! If j not equal to 0, the cell is not empty.
                  if (j .gt. 0) then
                     ! If jcell = icell, we are processing this cell with itself, so
                     ! j will be the next particle in the cell, namely, list(i)
4174                 if (jcell .eq. icell)  j=list(i)
 
                     do while (j .gt. 0)
                       idx_i = (i-1)*8
                       idx_j = (j-1)*8

                       dx = DBL_BUF1(idx_i+1)-DBL_BUF1(idx_j+1)+ x_trans
                       dy = DBL_BUF1(idx_i+2)-DBL_BUF1(idx_j+2)+ y_trans
                       dz = DBL_BUF1(idx_i+3)-DBL_BUF1(idx_j+3)+ z_trans

                       dr = dx**2 + dy**2 + dz**2

                       ! Are the particles close enough to interact?
                       if (dr .lt. 1.d0) then
                         iti = int(DBL_BUF1(idx_i+4))
                         itj = int(DBL_BUF1(idx_j+4))
                         if (iti .eq. itj) then

c////////////////////////////////////////////////////////////////////////
C First, search for A- (polymer-) clusters
                       if(iti.eq.1)then
                        if(clust_AA(i).eq.0.and.clust_AA(j).eq.0)then
                          cl_AA=cl_AA+1
                          clust_AA(i)=cl_AA
                          clust_AA(j)=cl_AA
                          size_AA(cl_AA)=2
                          head_AA(cl_AA)=i
                          list_AA(i)=j
                          nbr_clust_AA=nbr_clust_AA+1
                        elseif(clust_AA(i).eq.0)then
                          lj=clust_AA(j)
                          clust_AA(i)=lj
                          size_AA(lj)=size_AA(lj)+1
                          j1=head_AA(lj)
                          head_AA(lj)=i
                          list_AA(i)=j1
                        elseif(clust_AA(j).eq.0)then
                          li=clust_AA(i)
                          clust_AA(j)=li
                          size_AA(li)=size_AA(li)+1
                          i1=head_AA(li)
                          head_AA(li)=j
                          list_AA(j)=i1
                        elseif(clust_AA(i).gt.clust_AA(j))then
                          li0=clust_AA(i)
                          lj0=clust_AA(j)
                          j1=head_AA(lj0)
                          i1=head_AA(li0)
                          head_AA(lj0)=i1
                          i2=i1
                          do while(i2.gt.0)
                            i1=i2
                            clust_AA(i1)=lj0
                            i2=list_AA(i1)
                          enddo
                          list_AA(i1)=j1
                          head_AA(li0)=0
                          size_AA(lj0)=size_AA(lj0)+size_AA(li0)
                          size_AA(li0)=0
                          nbr_clust_AA=nbr_clust_AA-1
                        elseif(clust_AA(i).lt.clust_AA(j))then
                          lj0=clust_AA(j)
                          li0=clust_AA(i)
                          i1=head_AA(li0)
                          j1=head_AA(lj0)
                          head_AA(li0)=j1
                          j2=j1
                          do while(j2.gt.0)
                            j1=j2
                            clust_AA(j1)=li0
                            j2=list_AA(j1)
                          enddo
                          list_AA(j1)=i1
                          head_AA(lj0)=0
                          size_AA(li0)=size_AA(li0)+size_AA(lj0)
                          size_AA(lj0)=0
                          nbr_clust_AA=nbr_clust_AA-1
                        endif
C Second, Search for B- (solvent-) clusters
                       elseif(iti.eq.-1)then
                        if(clust_BB(i).eq.0.and.clust_BB(j).eq.0)then
                          cl_BB=cl_BB+1
                          clust_BB(i)=cl_BB
                          clust_BB(j)=cl_BB
                          size_BB(cl_BB)=2
                          head_BB(cl_BB)=i
                          list_BB(i)=j
                          nbr_clust_BB=nbr_clust_BB+1
                        elseif(clust_BB(i).eq.0)then
                          lj=clust_BB(j)
                          clust_BB(i)=lj
                          size_BB(lj)=size_BB(lj)+1
                          j1=head_BB(lj)
                          head_BB(lj)=i
                          list_BB(i)=j1
                        elseif(clust_BB(j).eq.0)then
                          li=clust_BB(i)
                          clust_BB(j)=li
                          size_BB(li)=size_BB(li)+1
                          i1=head_BB(li)
                          head_BB(li)=j
                          list_BB(j)=i1
                        elseif(clust_BB(i).gt.clust_BB(j))then
                          li0=clust_BB(i)
                          lj0=clust_BB(j)
                          j1=head_BB(lj0)
                          i1=head_BB(li0)
                          head_BB(lj0)=i1
                          i2=i1
                          do while(i2.gt.0)
                            i1=i2
                            clust_BB(i1)=lj0
                            i2=list_BB(i1)
                          enddo
                          list_BB(i1)=j1
                          head_BB(li0)=0
                          size_BB(lj0)=size_BB(lj0)+size_BB(li0)
                          size_BB(li0)=0
                          nbr_clust_BB=nbr_clust_BB-1
                        elseif(clust_BB(i).lt.clust_BB(j))then
                          lj0=clust_BB(j)
                          li0=clust_BB(i)
                          i1=head_BB(li0)
                          j1=head_BB(lj0)
                          head_BB(li0)=j1
                          j2=j1
                          do while(j2.gt.0)
                            j1=j2
                            clust_BB(j1)=li0
                            j2=list_BB(j1)
                          enddo
                          list_BB(j1)=i1
                          head_BB(lj0)=0
                          size_BB(li0)=size_BB(li0)+size_BB(lj0)
                          size_BB(lj0)=0
                          nbr_clust_BB=nbr_clust_BB-1
                        endif
                       endif
c end clusters search
c////////////////////////////////////////////////////////////////////////
                         endif ! end if (iti .eq. itj)
                       endif ! end if (dr .lt. 1)
                       j = list(j)
                     enddo ! end do-while (j .gt. 0)
                     j = head(jcell)
                     i = list(i)

                     ! If there are more particles, proceed
                     if (i .gt. 0) goto 4174
                  endif ! end (j .gt. 0)                
               enddo ! end kcell loop
            endif ! end (i .gt. 0)
          enddo ! end ix loop
        enddo ! end iy loop
      enddo ! end iz loop

      open(41,file='A_clusters.data',access='APPEND')
      open(42,file='B_clusters.data',access='APPEND')
      write(41,*)itime,nbr_clust_AA
      write(42,*)itime,nbr_clust_BB
      do k=1,cl_AA
        if(size_AA(k).gt.0)then
          write(41,*)size_AA(k)
        endif
      enddo
      do k=1,cl_BB
        if(size_BB(k).gt.0)then
          write(42,*)size_BB(k)
        endif
      enddo
      close(41)
      close(42)

      r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then
      !   write(*,*) 'find_clusters() takes: ', r_time2-r_time1
      !endif
      return
      end

c
c bf_find_bridges() - Locates bridges of polymers through 
c                  successive dot products.
c
      subroutine bf_find_bridges()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      double precision r_time1, r_time2
      integer id_i, id_j
      integer i,j,k,l,m,n,idx
      integer bridge_count, procd, ttl_count,offset
      double precision v_dot, v1mag, v2mag
      double precision theta, net_theta

      if (DISABLE_BRIDGES .eq. YES) return
      !! Timing
      r_time1 = MPI_WTIME()

      ! Initialize
      bridge_count = 0
 
      ! If no polymers, we have no need to be here.
      if (chain_a .eq. 1 .and. chain_b .eq. 1) return

      do i=1, SUPP_SLAB_SIZE
         BRIDGE_LIST(i)=0
      enddo
      TOTAL_A_BRIDGE = 0
      TOTAL_B_BRIDGE = 0

      ! Let's loop and quantify how straight each polymer chain is by means of 
      ! dot products between adjacent orientation vectors between monomers.
      do i=1, cluster_size
        do j=1, a_poly

          ! First we find direction vectors
          do k=1, chain_a - 1
             idx  = (i-1)*local_no_orig + (j-1)*chain_a + k

             id_i = members(idx)
             id_j = members(idx+1)

             ! Get the components of the vector between adjacent monomers.
             if ( (id_i .ne. 0 .and. id_j .ne. 0) .and.
     >            (id_i .le. local_no .or. id_j .le. local_no)) then

                 l            = (idx - 1)*3
                 DBL_BUF(l+1) = poly_x(id_i) - poly_x(id_j)
                 DBL_BUF(l+2) = poly_y(id_i) - poly_y(id_j)
                 DBL_BUF(l+3) = poly_z(id_i) - poly_z(id_j)


                 v1mag = dsqrt(DBL_BUF(l+1)**2+DBL_BUF(l+2)**2+
     >                         DBL_BUF(l+3)**3)

                 DBL_BUF(l+1) = DBL_BUF(l+1)/v1mag
                 DBL_BUF(l+2) = DBL_BUF(l+2)/v1mag
                 DBL_BUF(l+3) = DBL_BUF(l+3)/v1mag
             endif
          enddo

          ! Now perform the dot products to determine if we have a candidate for a 
          ! bridge or not.
          do k=1, chain_a - BRIDGE_LENGTH_A
             net_theta = 0.d0
             idx       = (i-1)*local_no_orig + (j-1)*chain_a+k

             procd = 0
             do l=1, BRIDGE_LENGTH_A
               id_i = members(idx)
               id_j = members(idx+l)

               if ( (id_i .ne. 0 .and. id_j .ne. 0) .and.
     >              (id_i .le. local_no .or. id_j .le. local_no)) then
      
                   procd = procd + 1
                   m     = (idx - 1)*3
                   n     = (idx + l - 1)*3

                   v1mag = dsqrt(DBL_BUF(m+1)**2+DBL_BUF(m+2)**2+
     >                           DBL_BUF(m+3)**2)
                   v2mag = dsqrt(DBL_BUF(n+1)**2+DBL_BUF(n+2)**2+
     >                           DBL_BUF(n+3)**2)

                   v_dot = DBL_BUF(m+1)*DBL_BUF(n+1) +
     >                     DBL_BUF(m+2)*DBL_BUF(n+2) + 
     >                     DBL_BUF(m+3)*DBL_BUF(n+3)

                   if (v_dot .ge. 0.d0) then
                       theta     = dacos(v_dot/v1mag/v2mag)
                   else
                       theta     = pi/2.d0
                   endif
                   net_theta = net_theta + theta
               endif        
             enddo
            
             ! We have calculated the average orientation for the segment of length
             ! BRIDGE_LENGTH_A. Now... is it within our tolerance? 
             if ( (id_i .ne. 0 .and. id_j .ne. 0) .and.
     >            (id_i .le. local_no .or. id_j .le. local_no)) then
                  net_theta = net_theta/procd * 180/pi
                  if (net_theta .le. BRIDGE_EPSILON) then
                    if(p_flag(id_i) .ne.0) then
                     bridge_count           = bridge_count + 1
                     DBL_BUF1(bridge_count) = p_flag(id_i)
                    endif
                  endif
             endif
          enddo
        enddo
      enddo

      ! We need to collect the results at the master node, now.
      call MPI_REDUCE(bridge_count,ttl_count,1,MPI_INTEGER,MPI_SUM,
     >                MASTER,MPI_COMM_WORLD,mpi_err)
      if (cluster_rank .eq. MASTER) then
         do i=1, bridge_count
            BRIDGE_LIST(i)=int(DBL_BUF1(i))
         enddo
         ttl_count = bridge_count
      endif

      do i=MASTER+1, LAST_SLAVE
         ! Send the results to the master node if we're not the maser node.
         if (i .eq. cluster_rank .and. cluster_rank .ne. MASTER) then
            call MPI_SEND(bridge_count, 1, MPI_INTEGER, MASTER,
     >                    cluster_rank,MPI_COMM_WORLD,mpi_err)

            ! Only bother sending if there's something to send.
            if (bridge_count .gt. 0) then
               call MPI_SEND(DBL_BUF1,bridge_count,MPI_DOUBLE_PRECISION,
     >                       MASTER,cluster_rank,MPI_COMM_WORLD,mpi_err)
            endif 
         ! Otherwise, the master receives the information and writes it to disk.
         elseif (cluster_rank .eq. MASTER) then
            call MPI_RECV(bridge_count,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,
     >                    mpi_status(10),mpi_err)

            if (bridge_count .gt. 0) then
               call MPI_RECV(DBL_BUF1,bridge_count,MPI_DOUBLE_PRECISION,
     >                       i,i,MPI_COMM_WORLD,mpi_status(10),mpi_err)

              do j=1, bridge_count
                 BRIDGE_LIST(i+ttl_count)=int(DBL_BUF1(j))
              enddo
            endif
            ttl_count = ttl_count + bridge_count
         endif
      enddo

      offset = ttl_count
      if (cluster_rank .eq. MASTER) TOTAL_A_BRIDGE=offset

      ! Time to look for B bridges.
      if (chain_b .gt. 1) then
       bridge_count = 0
       ttl_count    = 0
       do i=1, cluster_size
        do j=1, b_poly

          ! First we find direction vectors
          do k=1, chain_b - 1
             idx  = (i-1)*local_no_orig + a_poly*chain_a + (j-1)*chain_b+k

             id_i = members(idx)
             id_j = members(idx+1)

             ! Get the components of the vector between adjacent monomers.
             if ( (id_i .ne. 0 .and. id_j .ne. 0) .and.
     >            (id_i .le. local_no .or. id_j .le. local_no)) then

                 l            = (idx - 1)*3
                 DBL_BUF(l+1) = poly_x(id_i) - poly_x(id_j)
                 DBL_BUF(l+2) = poly_y(id_i) - poly_y(id_j)
                 DBL_BUF(l+3) = poly_z(id_i) - poly_z(id_j)


                 v1mag = dsqrt(DBL_BUF(l+1)**2+DBL_BUF(l+2)**2+
     >                         DBL_BUF(l+3)**3)

                 DBL_BUF(l+1) = DBL_BUF(l+1)/v1mag
                 DBL_BUF(l+2) = DBL_BUF(l+2)/v1mag
                 DBL_BUF(l+3) = DBL_BUF(l+3)/v1mag
             endif
          enddo

          ! Now perform the dot products to determine if we have a candidate for a 
          ! bridge or not.
          do k=1, chain_b - BRIDGE_LENGTH_B
             net_theta = 0.d0
             idx       =(i-1)*local_no_orig+a_poly*chain_a+(j-1)*chain_b+k

             procd = 0
             do l=1, BRIDGE_LENGTH_B
               id_i = members(idx)
               id_j = members(idx+l)

               if ( (id_i .ne. 0 .and. id_j .ne. 0) .and.
     >              (id_i .le. local_no .or. id_j .le. local_no)) then
      
                   procd = procd + 1
                   m     = (idx - 1)*3
                   n     = (idx + l - 1)*3

                   v1mag = dsqrt(DBL_BUF(m+1)**2+DBL_BUF(m+2)**2+
     >                           DBL_BUF(m+3)**2)
                   v2mag = dsqrt(DBL_BUF(n+1)**2+DBL_BUF(n+2)**2+
     >                           DBL_BUF(n+3)**2)

                   v_dot = DBL_BUF(m+1)*DBL_BUF(n+1) +
     >                     DBL_BUF(m+2)*DBL_BUF(n+2) + 
     >                     DBL_BUF(m+3)*DBL_BUF(n+3)

                   if (v_dot .ge. 0.d0) then
                       theta     = dacos(v_dot/v1mag/v2mag)
                   else
                       theta     = pi/2.d0
                   endif
                   net_theta = net_theta + theta
               endif        
             enddo
            
             ! We have calculated the average orientation for the segment of length
             ! BRIDGE_LENGTH_B. Now... is it within our tolerance? 
             if ( (id_i .ne. 0 .and. id_j .ne. 0) .and.
     >            (id_i .le. local_no .or. id_j .le. local_no)) then
                  net_theta = net_theta/procd * 180/pi
                  if (net_theta .le. BRIDGE_EPSILON) then
                    if (p_flag(id_i) .ne. 0) then
                     bridge_count           = bridge_count + 1
                     DBL_BUF1(bridge_count) = p_flag(id_i)
                    endif
                  endif
             endif
          enddo
        enddo
       enddo

       ! We need to collect the results at the master node, now.
       call MPI_REDUCE(bridge_count,ttl_count,1,MPI_INTEGER,MPI_SUM,
     >                MASTER,MPI_COMM_WORLD,mpi_err)
       if (cluster_rank .eq. MASTER) then
         do i=1, bridge_count
            BRIDGE_LIST(offset+i)=int(DBL_BUF1(i))
         enddo
         ttl_count = bridge_count
         offset    = offset + bridge_count
       endif

       do i=MASTER+1, LAST_SLAVE
         ! Send the results to the master node if we're not the maser node.
         if (i .eq. cluster_rank .and. cluster_rank .ne. MASTER) then
            call MPI_SEND(bridge_count, 1, MPI_INTEGER, MASTER,
     >                    cluster_rank,MPI_COMM_WORLD,mpi_err)

            ! Only bother sending if there's something to send.
            if (bridge_count .gt. 0) then
               call MPI_SEND(DBL_BUF1,bridge_count,MPI_DOUBLE_PRECISION,
     >                       MASTER,cluster_rank,MPI_COMM_WORLD,mpi_err)
            endif 
         ! Otherwise, the master receives the information and writes it to disk.
         elseif (cluster_rank .eq. MASTER) then
            call MPI_RECV(bridge_count,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,
     >                    mpi_status(10),mpi_err)

            if (bridge_count .gt. 0) then
               call MPI_RECV(DBL_BUF1,bridge_count,MPI_DOUBLE_PRECISION,
     >                       i,i,MPI_COMM_WORLD,mpi_status(10),mpi_err)

               do j=1, bridge_count
                  BRIDGE_LIST(offset+j)=int(DBL_BUF1(j))
               enddo
            endif
            ttl_count = ttl_count + bridge_count
            offset    = offset + bridge_count
         endif
       enddo
      endif

      if (cluster_rank .eq. MASTER) TOTAL_A_BRIDGE=ttl_count

      !! Timing
      r_time2 = MPI_WTIME()
      if (cluster_rank .eq. MASTER) then
         !write(*,*) 'find_bridges() takes: ', r_time2-r_time1
      endif
      return
      end

c
c bf_finalize_bridges() - Finalizes the bridge detection by outputting the numbers and coordinates for
c                      the polymers in question. Called ONLY by MASTER.
c    
      subroutine bf_finalize_bridges()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      double precision r_time1, r_time2
      integer i, idx

      if (DISABLE_BRIDGES .eq. YES) return

      open(12, file="A_bridges.data", access="append")
      open(13, file="B_bridges.data", access="append")
   
      write(12,*) itime, TOTAL_A_BRIDGE
      write(13,*) itime, TOTAL_B_BRIDGE

      do i=1, TOTAL_A_BRIDGE+TOTAL_B_BRIDGE
         call bf_write_poly_master_level(BRIDGE_LIST(i))
      enddo

      close(12)
      close(13)
      return      
      end
c 
c bf_local_velocity() - Calculates the average local velocity
c                    of particles in the system.
c
      subroutine bf_local_velocity(npoints)
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      character*18 op_filename
      character*1  rank2
      character*1  rank1
      integer i,j,k, npoints
      double precision v_local
      double precision r_time1, r_time2

      !! Timing
      r_time1 = MPI_WTIME()

      ! Find the proper filename for recording.
      if (npoints .lt. 10) then
         rank1 = char(npoints+48)
         op_filename = "local_velo0" // rank1 // ".data"     
      else
         rank1 = char(int(npoints/10)+48)
         rank2 = char(mod(npoints,10)+48)
         op_filename = "local_velo" //rank1 // rank2 // ".data"     
      endif

      ! Local momentum calculation
      do i=1, npoints**3
         k = (i-1)*4
         DBL_BUF(k+1) = 0.d0
         DBL_BUF(k+2) = 0.d0
         DBL_BUF(k+3) = 0.d0
         DBL_BUF(k+4) = 0.d0
      enddo

      ! Add the velocities based on where the particle is.
      do i=1, local_no
         j = 1 + int(x(i)/npoints) + int(y(i)/npoints)*npoints + 
     >         int(z(i)/npoints)*npoints**2
 
         k = (j-1)*4
         DBL_BUF(k+1) = DBL_BUF(k+1) + vx(i)
         DBL_BUF(k+2) = DBL_BUF(k+2) + vy(i)
         DBL_BUF(k+3) = DBL_BUF(k+3) + vz(i)
         DBL_BUF(k+4) = DBL_BUF(k+4) + 1
      enddo

      ! Reduce the local velocities from other nodes.
      k = 4*npoints**3
      call MPI_REDUCE(DBL_BUF, DBL_BUF1,k,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
      
      ! Let the master combine and average.
      v_local = 0.d0
      if (cluster_rank .eq. MASTER) then
        do i=1, npoints**3
           k = (i-1)*4
           if (DBL_BUF1(k+4) .le. 1.d0) DBL_BUF1(k+4)=1.d0
           
           ! Average
           DBL_BUF1(k+1) = DBL_BUF1(k+1)/DBL_BUF1(k+4)
           DBL_BUF1(k+2) = DBL_BUF1(k+2)/DBL_BUF1(k+4)
           DBL_BUF1(k+3) = DBL_BUF1(k+3)/DBL_BUF1(k+4)
            
           v_local = v_local + dsqrt(DBL_BUF1(k+1)**2 +
     >                               DBL_BUF1(k+2)**2 + 
     >                               DBL_BUF1(k+3)**2)
        enddo

        ! Average out.
        v_local = v_local / npoints**3
        
        ! Write to file.
        open(17,file=op_filename,access='append')
        write(17,*) itime, v_local
        close(17)

        ! Timing
        r_time2 = MPI_WTIME()
        !write(*,*) 'local_velocity() takes: ', r_time2-r_time1
      endif
      return
      end

c
c bf_radius_gyration() - Calculates the radii of gyration for
c                     all polymers.
c
      subroutine bf_radius_gyration()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer poly_total, i, j, k, h,l
      double precision r_time1, r_time2
      double precision x_cm, y_cm, z_cm
      double precision rg_a, rg_b
      double precision tmp_a, tmp_b

      ! Exit if we have no purpose here, i.e. solvent only.
      if (chain_a .eq. 1 .and. chain_b .eq. 1) return
 
      tmp_a = 0.d0
      tmp_b = 0.d0

      r_time1 = MPI_WTIME()

      if (cluster_rank .eq. MASTER) then
         open (12, file='radius_gyr_a.data', access='append')
         open (13, file='radius_gyr_b.data', access='append')
         open (14, file='poly_bond_a.data', access='append')
         open (15, file='poly_bond_b.data', access='append')
      endif
     
      ! First, clear out array.
      do i=1, 3*particle_count
         DBL_BUF(i)  = 0.d0
         DBL_BUF1(i) = 0.d0
      enddo
    
      ! This is quick: calculate mean bond lengths.
      call MPI_REDUCE(poly_bond_a,tmp_a,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
      call MPI_REDUCE(poly_bond_b,tmp_b,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)

      if (cluster_rank .eq. MASTER) then
 
         ! Mean A bond length.
         write(14,*) itime, tmp_a/(chain_a-1)/cluster_size/a_poly

         ! Mean B bond length.
         write(15,*) itime, tmp_b/(chain_b-1)/cluster_size/b_poly

         close(14)
         close(15)
      endif
     
      ! We need to first calculate the center of mass for each chain.
      ! We can do this in a piece-wise fashion. 
      do i=1, local_no
         j = p_flag(i)
         k = (j-1)*3
         if (species(i) .eq. 1) l = chain_a
         if (species(i) .eq.-1) l = chain_b
         DBL_BUF(k+1)  = DBL_BUF(k+1)+poly_x(i)/l
         DBL_BUF(k+2)  = DBL_BUF(k+2)+poly_y(i)/l
         DBL_BUF(k+3)  = DBL_BUF(k+3)+poly_z(i)/l
      enddo

      ! Collect results at master level.
      if (chain_b .gt. 1) then
         k = 3*(cluster_size*n_poly)
      else
         k = 3*(cluster_size*n_poly)
      endif

      call MPI_REDUCE(DBL_BUF, DBL_BUF1, k,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
      
      
      ! Broadcast this information to all nodes.
      ! k = number of elements to send.
      if (chain_b .gt. 1) then
         k = 3*(cluster_size*n_poly)
      else
         k = 3*(cluster_size*n_poly)
      endif

      call MPI_BCAST(DBL_BUF1,k,MPI_DOUBLE_PRECISION,MASTER,
     >                   MPI_COMM_WORLD,mpi_err)
     
      ! We now need to calculate the radius of gyration from these values.
      rg_a = 0.d0
      rg_b = 0.d0
      
      do i=1, local_no
         j    = p_flag(i)
         k    = (j-1)*3
         x_cm = DBL_BUF1(k+1)
         y_cm = DBL_BUF1(k+2)
         z_cm = DBL_BUF1(k+3)
         if (species(i) .eq. 1) then
           rg_a = rg_a + (poly_x(i)-x_cm)**2 +
     >                   (poly_y(i)-y_cm)**2 +
     >                   (poly_z(i)-z_cm)**2
         else
           rg_b = rg_b + (poly_x(i)-x_cm)**2 +
     >                   (poly_y(i)-y_cm)**2 +
     >                   (poly_z(i)-z_cm)**2         
         endif
      enddo
      
      ! Sum across all nodes, collect at master node.
      call MPI_REDUCE(rg_a, tmp_a,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
      call MPI_REDUCE(rg_b, tmp_b,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)

      if (cluster_rank .eq. MASTER) then
       if(chain_a .gt. 1) then
         rg_a = tmp_a / (cluster_size * a_poly * chain_a)
         write(12,*) itime, dsqrt(rg_a) 
       endif

       if (chain_b .gt. 1) then
         write(13,*) itime, dsqrt(tmp_b/(chain_b * 
     >                      cluster_size * b_poly))
       endif
         close(12)
         close(13)
      endif
  
      r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then
      !   write(*,*) 'radius_gyration takes: ', r_time2-r_time1
      !endif
      return
      end


c
c bf_init_stress() - Sets the stresses to 0 -- I just wanted to get it out of the main code to keep the code
c                 tidy and organized.
c
      subroutine bf_init_stress()
      implicit none
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      
      sigma_xx = 0.d0
      sigma_xy = 0.d0
      sigma_xz = 0.d0
      sigma_yy = 0.d0
      sigma_yz = 0.d0
      sigma_zz = 0.d0

      return
      end

c
c bf_find_stress() - Takes the velocity, position, etc. of particle i and adds its contribution
c                 to the system's stresses. Stress calculations are finalized in finalize_stress().
c
      subroutine bf_add_stress(dx,dy,dz,force,intra)
      implicit none
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer intra
      double precision force

      ! Diagonal terms
      sigma_xx = sigma_xx - (dx**2)*force
      sigma_yy = sigma_yy - (dy**2)*force
      sigma_zz = sigma_zz - (dz**2)*force
      
      ! Off-diagonal terms -- if not an intrachain interaction.
      if (intra .eq. 0) then
        sigma_xy = sigma_xy - dx*dy*force
        sigma_xz = sigma_xz - dx*dz*force
        sigma_yz = sigma_yz - dy*dz*force
      endif

      return
      end

c
c bf_finalize_stress() - Take velocities into account, collect results at Master level
c                     and write to disk.
c
      subroutine bf_finalize_stress()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer i
      double precision s_xx,s_yy,s_zz,s_xy,s_xz,s_yz
      double precision r_time1, r_time2

      !! Timing
      r_time1 = MPI_WTIME()

      ! Add velocities:
      do i=1, local_no
        sigma_xx = sigma_xx + vx(i)**2
        sigma_xy = sigma_xy + vx(i)*vy(i)
        sigma_xz = sigma_xz + vx(i)*vz(i)
        sigma_yy = sigma_yy + vy(i)**2
        sigma_yz = sigma_yz + vy(i)*vz(i)
        sigma_zz = sigma_zz + vz(i)**2
      enddo 

      ! Collect results at master level.
      call MPI_REDUCE(sigma_xx,s_xx,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER,MPI_COMM_WORLD,mpi_err)
 
      call MPI_REDUCE(sigma_xy,s_xy,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER,MPI_COMM_WORLD,mpi_err)

      call MPI_REDUCE(sigma_xz,s_xz,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER,MPI_COMM_WORLD,mpi_err)

      call MPI_REDUCE(sigma_yy,s_yy,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER,MPI_COMM_WORLD,mpi_err)

      call MPI_REDUCE(sigma_yz,s_yz,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER,MPI_COMM_WORLD,mpi_err)

      call MPI_REDUCE(sigma_zz,s_zz,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER,MPI_COMM_WORLD,mpi_err)

      if (cluster_rank .eq. MASTER) then
        open(12, file='stress.data', access='append')
        write(12,34) itime,s_xx/particle_count, 
     >                    s_xy/particle_count,
     >                    s_xz/particle_count,
     >                    s_yy/particle_count,
     >                    s_yz/particle_count,
     >                    s_zz/particle_count
        close(12)

       !r_time2 = MPI_WTIME()
       !write(*,*) 'finalize_stress() takes: ', r_time2-r_time1
      endif
34    format(1I8,' ', E15.8, ' ', E15.8, ' ', E15.8, ' ', E15.8, ' ',  
     >       E15.8,' ', E15.8)
      return
      end
