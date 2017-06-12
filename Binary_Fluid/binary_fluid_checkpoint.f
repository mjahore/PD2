c//////////////////////////////////////////////
c checkpoint.f - For resuming simulations.
c
c
c Mike Hore, University of Memphis Summer 2004
c//////////////////////////////////////////////


c
c bf_write_checkpoint_master_level() - Collects all information from all nodes
c                                   and writes it to disk. Very slow. Use carefully!
c
      subroutine bf_write_checkpoint_master_level()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer i,k,j, l,clust_idx
      integer recv_count
      double precision r_time1, r_time2

      if (DISABLE_CHECKPOINT .eq. YES) return

      !! Timing
      r_time1 = MPI_WTIME()


      ! If master, open files and dump particle data to disk.
      if (cluster_rank .eq. MASTER) then
         clust_idx = local_no
         write(*,*) ' '
         write(*,*) 'Writing checkpoint at master level...'
         open (12, file='checkpoint.data')
         do i=1, cluster_size*n_poly
            write(12,*) p_species(i)
         enddo
         write(12,*) itime
         write(12,*) cluster_rank, local_no
         do i=1, local_no
            ! Store particle positions/species for
            ! cluster searches
            k = (i-1)*4
            DBL_BUF1(k+1)=x(i)
            DBL_BUF1(k+2)=y(i)
            DBL_BUF1(k+3)=z(i)
            DBL_BUF1(k+4)=species(i)

            write(12,35) x(i), y(i), z(i)
            write(12,35) poly_x(i), poly_y(i), poly_z(i)
            write(12,35) fx(i), fy(i), fz(i)
            write(12,35) vx(i), vy(i), vz(i)
            write(12,*)  int(logical_id(i)), int(species(i)),
     >                   int(p_flag(i))
         enddo
      endif

      ! Send information node-by-node to the master.
      do i=MASTER+1, LAST_SLAVE
         if (i .eq. cluster_rank) then
            do j=1, local_no
               k = (j-1)*15
               DBL_BUF(k+1)  = x(j)
               DBL_BUF(k+2)  = y(j)
               DBL_BUF(k+3)  = z(j)
               DBL_BUF(k+4)  = poly_x(j)
               DBL_BUF(k+5)  = poly_y(j)
               DBL_BUF(k+6)  = poly_z(j)
               DBL_BUF(k+7)  = fx(j)
               DBL_BUF(k+8)  = fy(j)
               DBL_BUF(k+9)  = fz(j)
               DBL_BUF(k+10) = vx(j)
               DBL_BUF(k+11) = vy(j)
               DBL_BUF(k+12) = vz(j)
               DBL_BUF(k+13) = logical_id(j)
               DBL_BUF(k+14) = species(j)
               DBL_BUF(k+15) = p_flag(j)
            enddo

            ! Send transfer count to master
            call MPI_SEND(local_no, 1, MPI_INTEGER, MASTER, i, 
     >                    MPI_COMM_WORLD, mpi_err)
            ! Temporary
            k = local_no*15
            j = i + 1
            call MPI_SEND(DBL_BUF, k, MPI_DOUBLE_PRECISION, MASTER,
     >                    j, MPI_COMM_WORLD, mpi_err)
         else if (cluster_rank .eq. MASTER) then
            write(*,*) 'write_checkpoint_master_level(): MASTER <- NODE'
     >                 , ' ',i
            call MPI_RECV(recv_count, 1, MPI_INTEGER, i, i, 
     >                    MPI_COMM_WORLD, mpi_status(10), mpi_err)

            ! Temporary
            k = recv_count*15
            j = i + 1
            call MPI_RECV(DBL_BUF, k, MPI_DOUBLE_PRECISION, i, j,
     >                    MPI_COMM_WORLD, mpi_status(10), mpi_err)

            ! We have received all of the data from node i, let us now write
            ! it to disk.
            write(12,*) i, recv_count
            do j=1, recv_count
               k = (j-1)*15

               ! Store particle positions/species for
               ! cluster searches
               clust_idx = clust_idx + 1
               l = (clust_idx-1)*4
               DBL_BUF1(l+1)=DBL_BUF(k+1)
               DBL_BUF1(l+2)=DBL_BUF(k+2)
               DBL_BUF1(l+3)=DBL_BUF(k+3)
               DBL_BUF1(l+4)=DBL_BUF(k+14)

               write(12,35) DBL_BUF(k+1),  DBL_BUF(k+2),  DBL_BUF(k+3)
               write(12,35) DBL_BUF(k+4),  DBL_BUF(k+5),  DBL_BUF(k+6)
               write(12,35) DBL_BUF(k+7),  DBL_BUF(k+8),  DBL_BUF(k+9)
               write(12,35) DBL_BUF(k+10), DBL_BUF(k+11), DBL_BUF(k+12)
               write(12,*)  int(DBL_BUF(k+13)), int(DBL_BUF(k+14)), 
     >                      int(DBL_BUF(k+15))
            enddo
         endif
      enddo

      r_time2 = MPI_WTIME()
 
      if (cluster_rank .eq. MASTER) then 
         close(12)
      !   write(*,*) 'write_checkpoint_master_level() takes: ', 
      !>               r_time2-r_time1
      endif
35    format(E47.41, ' ', E47.41, ' ', E47.41)
      return
      end

c
c bf_read_checkpoint_master_level() - Restores a simulation from the master level.
c
      subroutine bf_read_checkpoint_master_level()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer i,k,j
      integer recv_count, local_no_tmp
      
      ! If master, open files and dump particle data to disk.
      if (cluster_rank .eq. MASTER) then
         write(*,*) ' '
         write(*,*) 'Restoring checkpoint at master level...'
         open (12, file='checkpoint.data')
         do i=1, cluster_size*n_poly
            read(12,*) p_species(i)
         enddo
         read(12,*) itime0
         read(12,*) k, local_no
         do i=1, local_no
            read(12,35) x(i), y(i), z(i)
            read(12,35) poly_x(i), poly_y(i), poly_z(i)
            read(12,35) fx(i), fy(i), fz(i)
            read(12,35) vx(i), vy(i), vz(i)
            read(12,*)  logical_id(i),species(i),p_flag(i)
            members(logical_id(i)) = i
         enddo
      endif

      ! Send information node-by-node to the master.
      do i=MASTER+1, LAST_SLAVE
         if (i .eq. cluster_rank) then
            call MPI_RECV(local_no, 1, MPI_INTEGER, MASTER, i, 
     >                    MPI_COMM_WORLD, mpi_status(10), mpi_err)

            call MPI_RECV(itime0, 1, MPI_INTEGER, MASTER, i, 
     >                    MPI_COMM_WORLD, mpi_status(10), mpi_err)

            ! Temporary
            k = local_no*15
            j = i + 1
            call MPI_RECV(DBL_BUF, k, MPI_DOUBLE_PRECISION, MASTER, j,
     >                    MPI_COMM_WORLD, mpi_status(10), mpi_err)

            ! Unpack data
            do j=1, local_no
               k = (j-1)*15
               x(j)          = DBL_BUF(k+1)
               y(j)          = DBL_BUF(k+2)
               z(j)          = DBL_BUF(k+3)
               poly_x(j)     = DBL_BUF(k+4)
               poly_y(j)     = DBL_BUF(k+5)
               poly_z(j)     = DBL_BUF(k+6)
               fx(j)         = DBL_BUF(k+7)
               fy(j)         = DBL_BUF(k+8)
               fz(j)         = DBL_BUF(k+9)
               vx(j)         = DBL_BUF(k+10)
               vy(j)         = DBL_BUF(k+11)
               vz(j)         = DBL_BUF(k+12)
               logical_id(j) = int(DBL_BUF(k+13))
               species(j)    = int(DBL_BUF(k+14))
               p_flag(j)     = int(DBL_BUF(k+15))

               members(logical_id(j)) = j
            enddo
         else if (cluster_rank .eq. MASTER) then
            write(*,*) 'read_checkpoint_master_level(): MASTER -> NODE '
     >                 ,i
            ! Read data, pack, and then send to node i.
            read(12,*) k, local_no_tmp
            do j=1, local_no_tmp
               k = (j-1)*15
               read(12,35) DBL_BUF(k+1),  DBL_BUF(k+2),  DBL_BUF(k+3) 
               read(12,35) DBL_BUF(k+4),  DBL_BUF(k+5),  DBL_BUF(k+6) 
               read(12,35) DBL_BUF(k+7),  DBL_BUF(k+8),  DBL_BUF(k+9) 
               read(12,35) DBL_BUF(k+10), DBL_BUF(k+11), DBL_BUF(k+12)
               read(12,*)  DBL_BUF(k+13), DBL_BUF(k+14), DBL_BUF(k+15)
            enddo

            ! Send transfer count to master
            call MPI_SEND(local_no_tmp, 1, MPI_INTEGER, i, i, 
     >                    MPI_COMM_WORLD, mpi_err)

            call MPI_SEND(itime0, 1, MPI_INTEGER, i, i, 
     >                    MPI_COMM_WORLD, mpi_err)
            ! Temporary
            k = local_no_tmp*15
            j = i + 1
            call MPI_SEND(DBL_BUF, k, MPI_DOUBLE_PRECISION, i,
     >                    j, MPI_COMM_WORLD, mpi_err)

         endif
      enddo

      if (cluster_rank .eq. MASTER) then
         close(12)
         write(*,*) 'Checkpoint restore complete.'
         write(*,*) ' '
      endif
35    format(E47.41, ' ', E47.41, ' ', E47.41)
      return
      end
     
