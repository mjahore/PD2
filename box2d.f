c////////////////////////////////////////////////
c box2d.f - Handing off of particles to other
c         nodes, linked-list stuff, halo, etc.
c
c Mike Hore, University of Memphis, Summer 2004
c///////////////////////////////////////////////




c
c distribute2d() - Keep particles within system to satisfy the periodic boundary
c                condition in (x,y,z) and then distribute particles to all nodes
c                (that is, distribute particles which leave this CPU to the 
c                neighboring CPU to which they travel).
c
      subroutine distribute2d()
      implicit none
      include 'mpif.h'
      include 'common.f'

      integer i,j,k,l,m,n,dest,idx
      integer export_count, keep_count, recv_count 
      double precision r_time1, r_time2
      double precision slab_max, slab_min,y_tmp
      integer xfer_count
      integer sent_count
 
      if (cluster_size .eq. 1) return
     
      !! Timing
      r_time1 = MPI_WTIME()
      
      ! Slab boundaries:
      slab_min = (cluster_rank)   *slab_y
      slab_max = (cluster_rank+1) *slab_y


      ! Check periodic boundary:
      call periodic_boundary_processing2d(1, local_no)

      ! Here we send particles up, down, and keep ones which are
      ! still in our slab. The x,y,z and vx,vy,vz arrays are then
      ! restructured.
      export_count   = 0
      keep_count     = 0
      do i=1, local_no
         if (y(i) .gt. slab_max .or. y(i) .lt. slab_min) then
            k = export_count*13
            export_count  = export_count + 1
            DBL_BUF(k+1)  = x(i)
            DBL_BUF(k+2)  = y(i)
            DBL_BUF(k+3)  = vx(i)
            DBL_BUF(k+4)  = vy(i)
            DBL_BUF(k+5)  = dfloat(species(i))
            DBL_BUF(k+6)  = vx1(i)
            DBL_BUF(k+7)  = vy1(i)
            DBL_BUF(k+8) = fx1(i)
            DBL_BUF(k+9) = fy1(i)
            DBL_BUF(k+10) = dfloat(logical_id(i))
            DBL_BUF(k+11) = poly_x(i)
            DBL_BUF(k+12) = poly_y(i)
            DBL_BUF(k+13) = dfloat(p_flag(i))
          else
            k = keep_count*13
            keep_count     = keep_count + 1
            DBL_BUF1(k+1)  = x(i)
            DBL_BUF1(k+2)  = y(i)
            DBL_BUF1(k+3)  = dfloat(species(i))
            DBL_BUF1(k+4)  = dfloat(logical_id(i))
            DBL_BUF1(k+5)  = poly_x(i)
            DBL_BUF1(k+6)  = poly_y(i)
            DBL_BUF1(k+7)  = vx(i)
            DBL_BUF1(k+8) = vy(i)
            DBL_BUF1(k+9) = vx1(i)
            DBL_BUF1(k+10) = vy1(i)
            DBL_BUF1(k+11) = fx1(i)
            DBL_BUF1(k+12) = fy1(i)
            DBL_BUF1(k+13) = dfloat(p_flag(i))
          endif
      enddo

      
      ! Clear out members array
      do i=1, particle_count+1
         members(i) = 0
      enddo

      ! Put kept particles back into arrays.
      do i=1, keep_count
         k    = (i-1)*13

         x(i)       = DBL_BUF1(k+1)
         y(i)       = DBL_BUF1(k+2)
         species(i) = int(DBL_BUF1(k+3))
         idx        = int(DBL_BUF1(k+4))
         poly_x(i)  = DBL_BUF1(k+5)
         poly_y(i)  = DBL_BUF1(k+6)
         vx(i)      = DBL_BUF1(k+7)
         vy(i)      = DBL_BUF1(k+8)
         vx1(i)     = DBL_BUF1(k+9)
         vy1(i)     = DBL_BUF1(k+10) 
         fx1(i)     = DBL_BUF1(k+11)
         fy1(i)     = DBL_BUF1(k+12) 
         p_flag(i)  = int(DBL_BUF1(k+13))

         ! Take ownership of this particle
         logical_id(i) = idx
         members(idx)  = i
      enddo

      n = keep_count + 1 ! Index used for unpacking received particles.

      ! Send particles that have drifted out of this node's control to their new
      ! controller. All particles that need to be sent are stored in DBL_BUF. They've
      ! been separated from those particles that remain on this node.
      !
      ! This is somewhat inefficient, but saves on RAM (from previous algorithms) and
      ! should be fast enough that no real difference should be noticed in terms of 
      ! performance.
      recv_count = 0
      sent_count = 0
      do i=MASTER, LAST_SLAVE
         if (i .eq. cluster_rank) then

           ! Loop through all nodes and send the particles to each node.
           do j=MASTER, LAST_SLAVE
              xfer_count = 0

              if (j .ne. cluster_rank) then
                do k=1, export_count
                   l     = (k-1) * 13             
                   y_tmp = DBL_BUF(l+2) 
 
                   ! Determine where this particle goes.
                   dest = y_tmp/slab_y
         
                   ! Add the particle to the outgoing list.
                   if (dest .eq. j) then
                      m = xfer_count*13
                      DBL_BUF1(m+1)  = DBL_BUF(l+1)
                      DBL_BUF1(m+2)  = DBL_BUF(l+2)
                      DBL_BUF1(m+3)  = DBL_BUF(l+3)
                      DBL_BUF1(m+4)  = DBL_BUF(l+4)
                      DBL_BUF1(m+5)  = DBL_BUF(l+5)
                      DBL_BUF1(m+6)  = DBL_BUF(l+6)
                      DBL_BUF1(m+7)  = DBL_BUF(l+7)
                      DBL_BUF1(m+8)  = DBL_BUF(l+8)
                      DBL_BUF1(m+9)  = DBL_BUF(l+9)
                      DBL_BUF1(m+10) = DBL_BUF(l+10)
                      DBL_BUF1(m+11) = DBL_BUF(l+11)
                      DBL_BUF1(m+12) = DBL_BUF(l+12)
                      DBL_BUF1(m+13) = DBL_BUF(l+13)
                      xfer_count     = xfer_count + 1
                   endif
                enddo ! End loop through particles (DBL_BUF)

                ! Send the xfer_count, and if it is greater than 0,
                ! we send DBL_BUF1.
                sent_count = sent_count + xfer_count
                call MPI_SEND(xfer_count, 1, MPI_INTEGER, j, j, 
     >                        MPI_COMM_WORLD, mpi_err)

                if (xfer_count .gt. 0) then
                  m = xfer_count*13
                  call MPI_SEND(DBL_BUF1,m,MPI_DOUBLE_PRECISION, j, j,
     >                          MPI_COMM_WORLD, mpi_err)
                endif
              endif ! Endif (j .ne. cluster_rank)
           enddo ! End loop through all nodes except current node.
         else
           ! Else we will wait for a message from node i telling us whether there are particles
           ! here to receive.
           call MPI_RECV(xfer_count, 1, MPI_INTEGER, i, 
     >                   cluster_rank,MPI_COMM_WORLD, mpi_status(10),
     >                   mpi_err)

           if (xfer_count .gt. 0) then
               recv_count = recv_count + xfer_count
               m = xfer_count*13
               call MPI_RECV(DBL_BUF1,m,MPI_DOUBLE_PRECISION, i,
     >                       cluster_rank,MPI_COMM_WORLD,mpi_status(10),
     >                       mpi_err)

               ! And unpack:
               do j=1, xfer_count
                  m = (j-1)*13
                  
                  x(n)       = DBL_BUF1(m+1)
                  y(n)       = DBL_BUF1(m+2)
                  vx(n)      = DBL_BUF1(m+3)
                  vy(n)      = DBL_BUF1(m+4)
                  species(n) = int(DBL_BUF1(m+5))
                  vx1(n)     = DBL_BUF1(m+6)
                  vy1(n)     = DBL_BUF1(m+7)
                  fx1(n)     = DBL_BUF1(m+8)
                  fy1(n)     = DBL_BUF1(m+9)
                  idx        = int(DBL_BUF1(m+10))
                  poly_x(n)  = DBL_BUF1(m+11)
                  poly_y(n)  = DBL_BUF1(m+12)
                  p_flag(n)  = int(DBL_BUF1(m+13))
  
                  ! Take ownership of this particle
                  logical_id(n) = idx
                  members(idx)  = n

                  ! And update our index.
                  n = n + 1
               enddo
           endif
         endif
      enddo


      ! Reset count of total DPD particles.
      local_no = keep_count + recv_count

      r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then 
      !   write(*,*) 'distribute takes: ', r_time2-r_time1
      !endif
      return
      end

c
c distribute_halo_down2d() - Distributes particles from the lower section of current slab to
c                          CPU_down so that neighboring slabs interact with each other.
c
      subroutine distribute_halo_down2d()
      implicit none
      include 'mpif.h'
      include 'common.f'

      integer i, j, k, idx
      integer xfer_count, recv_count, data_recv_count, xfer_tag,recv_tag
      double precision slab_min
      double precision r_time1, r_time2
      double precision z_trans
 

      if (cluster_size .eq. 1) return

      xfer_count = 0
      slab_min   = cluster_rank*slab_y
    
      r_time1 = MPI_WTIME()

      ! Clear ext_list array
      do i=1, ext_send_count
         ext_list(i) = 0
      enddo

      ! We loop through particles and find ones which need
      ! to be sent down to CPU_down. They meet the condition
      ! that they're 1 unit or less from the slab boundary.
      do i=1, local_no
         ! If the particle falls within range, add it
         if (y(i) .le. (slab_min+sd)) then
             j = xfer_count*9
             xfer_count = xfer_count + 1
             DBL_BUF1(j+1) = x(i)
             DBL_BUF1(j+2) = y(i)
             DBL_BUF1(j+3) = vx(i)
             DBL_BUF1(j+4) = vy(i)
             DBL_BUF1(j+5) = dfloat(species(i))
             DBL_BUF1(j+6) = dfloat(logical_id(i))
             DBL_BUF1(j+7) = poly_x(i)
             DBL_BUF1(j+8)= poly_y(i)
             DBL_BUF1(j+9)= dfloat(p_flag(i))

             ! Record particles that are sent to neighbor so that
             ! we can process the interactions on THIS node as well.
             ext_list(xfer_count) = i
         endif
      enddo

      
      ! Now, send/receive the particles in a manner similar to distribute() above.
      do i=MASTER, LAST_SLAVE
         ! Do we send?
         if (i .eq. cluster_rank) then
           ! Send # of data, and then actual data.
           call MPI_SEND(xfer_count, 1, MPI_INTEGER, CPU_down,
     >                   CPU_down, MPI_COMM_WORLD, mpi_err)
          
           ! Define xfer_tag
           xfer_tag = CPU_down
         
           ! Define j as the number of elements to be sent.
           j = xfer_count * 9
  
           ! Send data.
           call MPI_SEND(DBL_BUF1, j, MPI_DOUBLE_PRECISION, CPU_down,
     >                   xfer_tag, MPI_COMM_WORLD, mpi_err)  
         ! Or do we receive?
         else if (i .eq. CPU_up) then
           ! Receive # of data, and then actual data.
           call MPI_RECV(recv_count, 1, MPI_INTEGER, CPU_up,
     >                   cluster_rank, MPI_COMM_WORLD, mpi_status(10),
     >                   mpi_err)

           ! Define recv_tag
           recv_tag = cluster_rank

           ! Define j as the number of elements to be received.
           j = recv_count * 9

           ! Receive data.
           call MPI_RECV(DBL_BUF, j, MPI_DOUBLE_PRECISION, CPU_up,
     >                   recv_tag,MPI_COMM_WORLD,mpi_status(10),mpi_err)

           ! That's it. Unpack below.
         endif
      enddo

      ! Set global variable corresponding to the # of particles from our
      ! neighbor that are being considered. This will save a TINY amount
      ! of time when we return the resulting forces in distribute_halo_up().
      ext_recv_count = recv_count
      ext_send_count = xfer_count

      ! Unpack and add to arrays.
      do i=1, ext_recv_count
         j    = local_no + i
         k    = (i-1)*9

         ! Unpack coordinates.
         x(j) = DBL_BUF(k+1)
         y(j) = DBL_BUF(k+2)

         ! Unpack velocities.
         vx(j) = DBL_BUF(k+3)
         vy(j) = DBL_BUF(k+4)
         
         ! Grab species.
         species(j) = int(DBL_BUF(k+5))
         idx        = int(DBL_BUF(k+6))
 
         ! Unpack polymer coordinates
         poly_x(j) = DBL_BUF(k+7)
         poly_y(j) = DBL_BUF(k+8)
 
         p_flag(j) = int(DBL_BUF(k+9))

         ! Take ownership of this particle
         logical_id(j) = idx
         members(idx)  = j
      enddo

      ! Periodic boundary processing
      call periodic_boundary_processing2d(local_no+1, 
     >                                    local_no+ext_recv_count)

      r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then 
      !   write(*,*) 'distribute_halo_down takes: ', r_time2-r_time1
      !endif
      return
      end

c
c distribute_halo_up2d() - Distributes the forces calculated between this slab
c                        and its upper neighbor BACK to the neighbor so that
c                        the particles' positions can be integrated.
c
      subroutine distribute_halo_up2d()
      implicit none
      include 'mpif.h'
      include 'common.f'

      integer i, j, k, recv_count, xfer_count
      integer xfer_tag, recv_tag
      double precision r_time1, r_time2


      if (cluster_size .eq. 1) return

      ! Timing.
      r_time1 = MPI_WTIME()

      ! Pack arrays for distribution to upper node.
      do i=1, ext_recv_count
         ! j is the index in the fx,fy,fz arrays.
         j = local_no + i
         k = (i-1)*2

         ! Let's pack the forces into the DBL_BUF array, now.        
         DBL_BUF(k+1) = fx(j)
         DBL_BUF(k+2) = fy(j)
      enddo

      ! Loop through all nodes, sending and receiving as required.
      do i=MASTER, LAST_SLAVE
         if (i .eq. cluster_rank) then
            ! Send forces up.
            xfer_count = ext_recv_count*2
            xfer_tag   = CPU_up
            
            call MPI_SEND(DBL_BUF, xfer_count, MPI_DOUBLE_PRECISION,
     >                    CPU_up, xfer_tag, MPI_COMM_WORLD, mpi_err)
         else if (i .eq. CPU_down) then
            ! Receive forces.
            recv_count = ext_send_count*2
            recv_tag   = cluster_rank

            call MPI_RECV(DBL_BUF1, recv_count, MPI_DOUBLE_PRECISION,
     >                    CPU_down, recv_tag, MPI_COMM_WORLD,
     >                    mpi_status(10), mpi_err)
         endif
      enddo

      ! Unpack received data.
      do i=1, ext_send_count
         j = (i-1)*2
         k = ext_list(i)
 
         fx(k) = fx(k) + DBL_BUF1(j+1)
         fy(k) = fy(k) + DBL_BUF1(j+2)
      enddo

      r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then
      !   write(*,*) 'distribute_halo_up() takes: ', r_time2 - r_time1
      !endif
      return
      end


c
c place2d() - Place all particles into a linked-list for processing pair-wise
c           interactions.
c
      subroutine place2d()
      implicit none
      include 'mpif.h'
      include 'common.f'

      integer i,j,k,icell,idx
      integer mxy

      mxy = mx*my

      ! PBC for 1 node only
      if (cluster_size .eq. 1) then
         ext_recv_count = 0
         call periodic_boundary_processing2d(1,local_no)
      endif

      ! Clear out old data from list, head.
      do i=1, ncell
         head(i) = 0
      enddo

      do i=1, particle_count
         list(i) = 0
      enddo

      ! Consider new particle positions:
      do i=1, local_no+ext_recv_count
         icell = 1 + int(x(i)*cellx_inv) + int(y(i)*celly_inv)*mx

         ! Update table.
         j           = head(icell)
         head(icell) = i
         list(i)     = j
      enddo

      return
      end

c
c periodic_boundary_processing2d() - Ensures that all particles within
c                                  lower_bound and upper_bound (inclusively)
c                                  adhere to the periodic boundary condition.
c
      subroutine periodic_boundary_processing2d(lower_bound,upper_bound)
      implicit none
      include 'mpif.h'
      include 'common.f'
      integer lower_bound, upper_bound, i

      ! Consider XYZ periodicity:
      do i=lower_bound, upper_bound
         if (x(i) .le. 0.d0) then
             x(i) = x(i) + dim_x
         else if (x(i) .gt. dim_x) then 
             x(i) = x(i) - dim_x
         endif

         if (y(i) .le. 0.d0) then
             y(i) = y(i) + dim_y
         else if (y(i) .gt. dim_y) then 
             y(i) = y(i) - dim_y
         endif
      enddo

      return
      end
