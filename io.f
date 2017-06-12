c////////////////////////////////////////////////
c io.f - Routines that read/write from/to
c        disk.
c
c Mike Hore, University of Memphis, Summer 2004
c///////////////////////////////////////////////

c
c collect_master_level() - Collects coordinates from all nodes and stores them in
c                          DBL_BUF1.
      subroutine collect_master_level3d()
      implicit none
      include 'mpif.h'
      include 'common.f'
      integer i,k,j, l,clust_idx
      integer recv_count, offset
      double precision r_time1, r_time2


      !! Timing
      r_time1 = MPI_WTIME()


      ! If master, open files and dump particle data to disk.
      if (cluster_rank .eq. MASTER) then
         clust_idx = local_no
         do i=1, local_no
            ! Store particle positions/species for
            ! cluster searches
            k = (i-1)*9
            DBL_BUF1(k+1)  = x(i)
            DBL_BUF1(k+2)  = y(i)
            DBL_BUF1(k+3)  = z(i)
            DBL_BUF1(k+4)  = species(i)
            DBL_BUF1(k+5)  = p_flag(i)
            DBL_BUF1(k+6)  = vx(i)
            DBL_BUF1(k+7)  = vy(i)
            DBL_BUF1(k+8)  = vz(i)
            DBL_BUF1(k+9) = logical_id(i)
         enddo
      endif


      ! Send information node-by-node to the master.
      do i=MASTER+1, LAST_SLAVE
         if (i .eq. cluster_rank) then
            do j=1, local_no
               k = (j-1)*9
               DBL_BUF(k+1)   = x(j)
               DBL_BUF(k+2)   = y(j)
               DBL_BUF(k+3)   = z(j)
               DBL_BUF(k+4)   = species(j)
               DBL_BUF(k+5)   = p_flag(j) 
               DBL_BUF(k+6)   = vx(j)
               DBL_BUF(k+7)   = vy(j)
               DBL_BUF(k+8)   = vz(j)
               DBL_BUF(k+9)  = logical_id(j)
            enddo

            ! Send transfer count to master
            call MPI_SEND(local_no, 1, MPI_INTEGER, MASTER, i, 
     >                    MPI_COMM_WORLD, mpi_err)
            ! Temporary
            k = local_no*9
            j = i + 1
            call MPI_SEND(DBL_BUF, k, MPI_DOUBLE_PRECISION, MASTER,
     >                    j, MPI_COMM_WORLD, mpi_err)
         else if (cluster_rank .eq. MASTER) then
            call MPI_RECV(recv_count, 1, MPI_INTEGER, i, i, 
     >                    MPI_COMM_WORLD, mpi_status(10), mpi_err)

            ! Temporary
            k = recv_count*9
            j = i + 1
            call MPI_RECV(DBL_BUF, k, MPI_DOUBLE_PRECISION, i, j,
     >                    MPI_COMM_WORLD, mpi_status(10), mpi_err)

            ! We have received all of the data from node i, let us now write
            ! it to disk.
            do j=1, recv_count
               k = (j-1)*9

               ! Store particle positions/species for
               ! cluster searches
               clust_idx = clust_idx + 1
               l = (clust_idx-1)*9
               DBL_BUF1(l+1) = DBL_BUF(k+1)
               DBL_BUF1(l+2) = DBL_BUF(k+2)
               DBL_BUF1(l+3) = DBL_BUF(k+3)
               DBL_BUF1(l+4) = DBL_BUF(k+4)
               DBL_BUF1(l+5) = DBL_BUF(k+5)
               DBL_BUF1(l+6) = DBL_BUF(k+6)
               DBL_BUF1(l+7) = DBL_BUF(k+7)
               DBL_BUF1(l+8) = DBL_BUF(k+8)
               DBL_BUF1(l+9) = DBL_BUF(k+9)
            enddo
         endif
      enddo

      r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then 
      ! write(*,*) 'collect_master_level() takes: ',
      !>             r_time2-r_time1
      !endif
      return
      end


c
c collect_master_level_poly() - Collects coordinates from all nodes and stores them in
c                               DBL_BUF1.
      subroutine collect_master_level3d_poly()
      implicit none
      include 'mpif.h'
      include 'common.f'
      integer i,k,j, l,clust_idx
      integer recv_count, offset
      double precision r_time1, r_time2


      !! Timing
      r_time1 = MPI_WTIME()


      ! If master, open files and dump particle data to disk.
      if (cluster_rank .eq. MASTER) then
         clust_idx = local_no
         do i=1, local_no
            ! Store particle positions/species for
            ! cluster searches
            k = (i-1)*5
            DBL_BUF1(k+1) = poly_x(i)
            DBL_BUF1(k+2) = poly_y(i)
            DBL_BUF1(k+3) = poly_z(i)
            DBL_BUF1(k+4) = species(i)
            DBL_BUF1(k+5) = p_flag(i)
         enddo
      endif


      ! Send information node-by-node to the master.
      do i=MASTER+1, LAST_SLAVE
         if (i .eq. cluster_rank) then
            do j=1, local_no
               k = (j-1)*5
               DBL_BUF(k+1)  = poly_x(j)
               DBL_BUF(k+2)  = poly_y(j)
               DBL_BUF(k+3)  = poly_z(j)
               DBL_BUF(k+4)  = species(j)
               DBL_BUF(k+5)  = p_flag(j)
            enddo

            ! Send transfer count to master
            call MPI_SEND(local_no, 1, MPI_INTEGER, MASTER, i, 
     >                    MPI_COMM_WORLD, mpi_err)
            ! Temporary
            k = local_no*5
            j = i + 1
            call MPI_SEND(DBL_BUF, k, MPI_DOUBLE_PRECISION, MASTER,
     >                    j, MPI_COMM_WORLD, mpi_err)
         else if (cluster_rank .eq. MASTER) then
            call MPI_RECV(recv_count, 1, MPI_INTEGER, i, i, 
     >                    MPI_COMM_WORLD, mpi_status(10), mpi_err)

            ! Temporary
            k = recv_count*5
            j = i + 1
            call MPI_RECV(DBL_BUF, k, MPI_DOUBLE_PRECISION, i, j,
     >                    MPI_COMM_WORLD, mpi_status(10), mpi_err)

            ! We have received all of the data from node i, let us now write
            ! it to disk.
            do j=1, recv_count
               k = (j-1)*5

               ! Store particle positions/species for
               ! cluster searches
               clust_idx = clust_idx + 1
               l = (clust_idx-1)*5
               DBL_BUF1(l+1)=DBL_BUF(k+1)
               DBL_BUF1(l+2)=DBL_BUF(k+2)
               DBL_BUF1(l+3)=DBL_BUF(k+3)
               DBL_BUF1(l+4)=DBL_BUF(k+4)
               DBL_BUF1(l+5)=DBL_BUF(k+5)
            enddo
         endif
      enddo

      return
      end

c
c collect_master_level2d() - Collects coordinates from all nodes and stores them in
c                          DBL_BUF1.
      subroutine collect_master_level2d()
      implicit none
      include 'mpif.h'
      include 'common.f'
      integer i,k,j, l,clust_idx
      integer recv_count, offset
      double precision r_time1, r_time2


      !! Timing
      r_time1 = MPI_WTIME()


      ! If master, open files and dump particle data to disk.
      if (cluster_rank .eq. MASTER) then
         clust_idx = local_no
         do i=1, local_no
            ! Store particle positions/species for
            ! cluster searches
            k = (i-1)*6
            DBL_BUF1(k+1) = x(i)
            DBL_BUF1(k+2) = y(i)
            DBL_BUF1(k+3) = species(i)
            DBL_BUF1(k+4) = p_flag(i)
            DBL_BUF1(k+5) = vx(i)
            DBL_BUF1(k+6) = vy(i)
         enddo
      endif


      ! Send information node-by-node to the master.
      do i=MASTER+1, LAST_SLAVE
         if (i .eq. cluster_rank) then
            do j=1, local_no
               k = (j-1)*6
               DBL_BUF(k+1)  = x(j)
               DBL_BUF(k+2)  = y(j)
               DBL_BUF(k+3)  = species(j)
               DBL_BUF(k+4)  = p_flag(j) 
               DBL_BUF(k+5)  = vx(i)
               DBL_BUF(k+6)  = vy(i)
            enddo

            ! Send transfer count to master
            call MPI_SEND(local_no, 1, MPI_INTEGER, MASTER, i, 
     >                    MPI_COMM_WORLD, mpi_err)
            ! Temporary
            k = local_no*6
            j = i + 1
            call MPI_SEND(DBL_BUF, k, MPI_DOUBLE_PRECISION, MASTER,
     >                    j, MPI_COMM_WORLD, mpi_err)
         else if (cluster_rank .eq. MASTER) then
            call MPI_RECV(recv_count, 1, MPI_INTEGER, i, i, 
     >                    MPI_COMM_WORLD, mpi_status(10), mpi_err)

            ! Temporary
            k = recv_count*6
            j = i + 1
            call MPI_RECV(DBL_BUF, k, MPI_DOUBLE_PRECISION, i, j,
     >                    MPI_COMM_WORLD, mpi_status(10), mpi_err)

            ! We have received all of the data from node i, let us now write
            ! it to disk.
            do j=1, recv_count
               k = (j-1)*6

               ! Store particle positions/species for
               ! cluster searches
               clust_idx = clust_idx + 1
               l = (clust_idx-1)*6
               DBL_BUF1(l+1)=DBL_BUF(k+1)
               DBL_BUF1(l+2)=DBL_BUF(k+2)
               DBL_BUF1(l+3)=DBL_BUF(k+3)
               DBL_BUF1(l+4)=DBL_BUF(k+4)
               DBL_BUF1(l+5)=DBL_BUF(k+5)
               DBL_BUF1(l+6)=DBL_BUF(k+6)
            enddo
         endif
      enddo

      r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then 
      ! write(*,*) 'collect_master_level() takes: ',
      !>             r_time2-r_time1
      !endif
      return
      end

