c///////////////////////////////////////////////
c dynamics.f - Rouse Mode analysis of dynamics
c 
c
c
c Mike Hore, hore@case.edu
c Case Western Reserve University
c April 2020
c//////////////////////////////////////////////


c
c gen_calc_rouse_modes() - Calculates the Rouse modes of each
c                          polymer chain. See the following refs:
c       Kalathi et al., Macromolecules 2014, 47, 6925
c       Jiang et al., J. Chem. Phys. 2007, 126, 044901.
c
      subroutine gen_calc_rouse_modes(p_type, t)
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer p_type, n_modes, idx, mode
      integer i,j,k,p,ierr
      integer N
      integer bf_monomer_id, bf_polymer_id
      integer t, offset
      double precision r_time1, r_time2
      double precision Xp(SUPP_SYS_SIZE)
      double precision Xp_tmp(SUPP_SYS_SIZE)

      !! Timing
      r_time1 = MPI_WTIME()

      ! There are at most (N - 1) modes:
      if (p_type .eq. 1) then
         N = chain_a
      else if (p_type .eq. -1) then
         N = chain_b
      endif
      n_modes = N 

      ! This subroutine does NOT need collect_master_level3d() to 
      ! be called first.
      do i=1, SUPP_SYS_SIZE
         Xp(i)     = 0.d0
         Xp_tmp(i) = 0.d0
      enddo

      ! Calculate all of the modes! 
      do i=1, local_no
        if (species(i) .eq. p_type) then
           ! Limit the size of the array by removing the component we do
           ! not care about.
           idx = (bf_polymer_id(i) - 1) * 3 * (n_modes+1)

           ! We need to know the position of this monomer in the
           ! chain.
          
           j = bf_monomer_id(i)

           ! This is the calculation of the Rouse modes Xp:
           k = idx
           do p=0, n_modes
              Xp_tmp(k+1) = Xp_tmp(k+1) + sqrt(2.d0/N) *
     >                                    dcos(p*PI*(j-0.5)/N)*poly_x(i)

              Xp_tmp(k+2) = Xp_tmp(k+2) + sqrt(2.d0/N) * 
     >                                    dcos(p*PI*(j-0.5)/N)*poly_y(i)

              Xp_tmp(k+3) = Xp_tmp(k+3) + sqrt(2.d0/N) * 
     >                                    dcos(p*PI*(j-0.5)/N)*poly_z(i)
              k = k + 3
           enddo
        endif
      enddo
     
      ! Collect the modes at the master level.
      if (p_type .eq. 1) then
         N = 3 * (n_modes+1) * cluster_size * a_poly
         call MPI_REDUCE(Xp_tmp,Xp,N,MPI_DOUBLE_PRECISION,
     >                   MPI_SUM, MASTER,MPI_COMM_WORLD,mpi_err)
         if (cluster_rank .eq. MASTER) then
            open(55, file="rouse_mode_a.data", access="stream")
         endif
      else if (p_type .eq. -1) then
         N = 3 * (n_modes+1) * cluster_size * b_poly
         call MPI_REDUCE(Xp_tmp,Xp,N,MPI_DOUBLE_PRECISION,
     >                   MPI_SUM, MASTER,MPI_COMM_WORLD,mpi_err)
         if (cluster_rank .eq. MASTER) then
            open(55, file="rouse_mode_b.data", access="stream")
         endif
      endif

      ! Write the modes to disk.
      if (cluster_rank .eq. MASTER) then

         ! Output length of the file.
         write(55) t

         ! Move to the end.
         call FSEEK(55, 0, 2, ierr)

         if (p_type .eq. 1) then
            N      = a_poly
         else 
            N      = b_poly
         endif

         do i=1, cluster_size*N
            idx = (i-1)*(n_modes+1)

            do p=0, n_modes
               k = idx + 3*p
               write(55) Xp(k+1), Xp(k+2), Xp(k+3)
            enddo
         enddo
         close(55)
      endif

      r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then 
      !   write(*,*) 'gen_calc_rouse_modes() takes: ', r_time2-r_time1
      !endif
      return
      end


c
c gen_calc_np_rouse_modes() - Calculates the Rouse modes of each
c                           polymer chain. See the following refs:
c       THIS IS FOR NANOPARTICLE SIMULATIONS
c       Kalathi et al., Macromolecules 2014, 47, 6925
c       Jiang et al., J. Chem. Phys. 2007, 126, 044901.
c
      subroutine gen_calc_np_rouse_modes(p_type, t)
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      include 'Nanoparticle/nanoparticle_var.f'
      integer p_type, n_modes, idx, mode
      integer i,j,k,p,ierr
      integer N
      integer np_monomer_id, np_polymer_id
      integer t, offset
      double precision r_time1, r_time2
      double precision Xp(SUPP_SYS_SIZE)
      double precision Xp_tmp(SUPP_SYS_SIZE)

      !! Timing
      r_time1 = MPI_WTIME()

      ! There are at most (N - 1) modes:
      if (p_type .eq. 1) then
         N = chain_a
      else if (p_type .eq. -1) then
         N = chain_b
      else if (p_type .eq. 5) then
         N = chain_g
      endif
      n_modes = N

      ! This subroutine does NOT need collect_master_level3d() to 
      ! be called first.
      do i=1, SUPP_SYS_SIZE
         Xp(i)     = 0.d0
         Xp_tmp(i) = 0.d0
      enddo

      ! Calculate all of the modes! 
      do i=1, local_no
        if (species(i) .eq. p_type) then
           ! Limit the size of the array by removing the component we do
           ! not care about.
           idx = (np_polymer_id(i) - 1) * 3 * (n_modes+1)

           ! We need to know the position of this monomer in the
           ! chain.
          
           j = np_monomer_id(i)

           ! This is the calculation of the Rouse modes Xp:
           k = idx
           do p=0, n_modes
              Xp_tmp(k+1) = Xp_tmp(k+1) + sqrt(2.d0/N) *
     >                                    dcos(p*PI*(j-0.5)/N)*poly_x(i)

              Xp_tmp(k+2) = Xp_tmp(k+2) + sqrt(2.d0/N) * 
     >                                    dcos(p*PI*(j-0.5)/N)*poly_y(i)

              Xp_tmp(k+3) = Xp_tmp(k+3) + sqrt(2.d0/N) * 
     >                                    dcos(p*PI*(j-0.5)/N)*poly_z(i)
              k = k + 3
           enddo
        endif
      enddo
     
      ! Collect the modes at the master level.
      if (p_type .eq. 1) then
         N = 3 * (n_modes+1) * (cluster_size * a_poly)
         call MPI_REDUCE(Xp_tmp,Xp,N,MPI_DOUBLE_PRECISION,
     >                   MPI_SUM, MASTER,MPI_COMM_WORLD,mpi_err)
         if (cluster_rank .eq. MASTER) then
            !open(55, file="rouse_mode_a.data", access="append")
            open(55, file="rouse_mode_a.data", access="stream")
         endif
      else if (p_type .eq. -1) then
         N = 3 * (n_modes+1) * cluster_size * b_poly
         call MPI_REDUCE(Xp_tmp,Xp,N,MPI_DOUBLE_PRECISION,
     >                   MPI_SUM, MASTER,MPI_COMM_WORLD,mpi_err)
         if (cluster_rank .eq. MASTER) then
            open(55, file="rouse_mode_b.data", access="stream")
         endif
      else if (p_type .eq. 5) then
         N = 3 * (n_modes+1) * cluster_size * n_local * (1 + n_graft)
         call MPI_REDUCE(Xp_tmp,Xp,N,MPI_DOUBLE_PRECISION,
     >                   MPI_SUM, MASTER,MPI_COMM_WORLD,mpi_err)
         if (cluster_rank .eq. MASTER) then
            !open(55, file="rouse_mode_g.data", access="append")
            open(55, file="rouse_mode_g.data", access="stream")
         endif
      endif

      ! Write the modes to disk.
      if (cluster_rank .eq. MASTER) then

         ! Output length of the file.
         write(55) t

         ! Move to the end.
         call FSEEK(55, 0, 2, ierr)

         if (p_type .eq. 1) then
            N      = a_poly 
         else if (p_type .eq. -1) then
            N      = b_poly
         else if (p_type .eq. 5) then
            N      = n_local*n_graft
         endif

         do i=1, cluster_size*N
            idx = (i-1)*(n_modes+1)

            do p=0, n_modes
               k = idx + 3*p
               !write(55,*) Xp(k+1), Xp(k+2), Xp(k+3)
               write(55) Xp(k+1), Xp(k+2), Xp(k+3)
            enddo
         enddo
         close(55)
      endif

      r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then 
      !   write(*,*) 'gen_calc_rouse_modes() takes: ', r_time2-r_time1
      !endif
      return
      end
