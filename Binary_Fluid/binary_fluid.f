c////////////////////////////////////////////////////////
c binary_fluid.f - The bulk of the simulation: handles
c                  forces and integration of eqs. of 
c                  motion.
c
c Mike Hore, University of Memphis, Summer 2004
c///////////////////////////////////////////////////////

c
c bf_read_dist_parameters() - Reads parameters from file and distributes them to
c                             all slave nodes.
      subroutine bf_read_dist_parameters()
      implicit none
      double precision sd_box
      parameter(sd_box=1.d0)
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'

c If we are the master node, read in the configuration
      if (cluster_rank .eq. MASTER) then
          ! Output the program configuration (bond type, etc.)
          write(*,*) 'Parallel/MPI code version: ', code_version
          write(*,*) 'Master node detects ', cluster_size, ' nodes.'
          write(*,*) 'Simulation type is Binary Fluid.'
          write(*,*) ' '
          if (CP_LOAD_TYPE .eq. LOAD_AT_RELAX) then
             write(*,*) 'Checkpoints load at relax.'
          else
             write(*,*) 'Checkpoints load at quench.'
          endif
          if (BOND_TYPE .eq. FENE) then
             write(*,*) 'Using FENE bonds.'
          else
             write(*,*) 'Using harmonic bonds.'
          endif
          write(*,*) ' '
 
          ! Read parameters and distribute:
          call bf_read_config()

          ! Put all integer variables into MPI_INT_BUF1 for
          ! broadcast.
          MPI_INT_BUF1(1) = load_cp
          MPI_INT_BUF1(2) = chain_a
          MPI_INT_BUF1(3) = chain_b
          MPI_INT_BUF1(4) = time
          MPI_INT_BUF1(5) = equil_time
          MPI_INT_BUF1(6) = iseed
          MPI_INT_BUF1(7) = check_freq
          MPI_INT_BUF1(8) = calc_freq

          ! Put all double precision variables into MPI_DBL_BUF1 
          ! for broadcast.
          MPI_DBL_BUF1(1)  = rho
          MPI_DBL_BUF1(2)  = temp
          MPI_DBL_BUF1(3)  = chi1
          MPI_DBL_BUF1(4)  = chi2
          MPI_DBL_BUF1(5)  = chi3
          MPI_DBL_BUF1(6)  = d_els
          MPI_DBL_BUF1(7)  = dt
          MPI_DBL_BUF1(8)  = dim_x
          MPI_DBL_BUF1(9)  = dim_y
          MPI_DBL_BUF1(10) = dim_z
          MPI_DBL_BUF1(11) = sigma
          MPI_DBL_BUF1(12) = comp
          MPI_DBL_BUF1(13) = sd
          MPI_DBL_BUF1(14) = max_bond
          MPI_DBL_BUF1(15) = vol_a
          MPI_DBL_BUF1(16) = vol_b
          MPI_DBL_BUF1(17) = phi_diff

          ! Broadcast MPI_INT_BUF1()
          call MPI_BCAST(MPI_INT_BUF1,8,MPI_INTEGER,MASTER,
     >                   MPI_COMM_WORLD,mpi_err)

          ! Broadcast MPI_DBL_BUF1()
          call MPI_BCAST(MPI_DBL_BUF1,17,MPI_DOUBLE_PRECISION,
     >                   MASTER,MPI_COMM_WORLD,mpi_err)
          write(*,*) 'Master node has distributed simulation parameters'
          write(*,*) ' '
c Else, wait to receive configuration from the master
      else
          call MPI_BCAST(MPI_INT_BUF1,8,MPI_INTEGER,MASTER,
     >                   MPI_COMM_WORLD,mpi_err)
          
          call MPI_BCAST(MPI_DBL_BUF1,17,MPI_DOUBLE_PRECISION,MASTER,
     >                   MPI_COMM_WORLD,mpi_err)

          ! Unpack receive buffer(MPI_INT_BUF1)
          load_cp    = MPI_INT_BUF1(1)
          chain_a    = MPI_INT_BUF1(2)
          chain_b    = MPI_INT_BUF1(3)
          time       = MPI_INT_BUF1(4)
          equil_time = MPI_INT_BUF1(5)
          iseed      = MPI_INT_BUF1(6)/cluster_rank + cluster_rank*1000
          check_freq = MPI_INT_BUF1(7)
          calc_freq  = MPI_INT_BUF1(8)

          ! Unpack receive buffer(MPI_DBL_BUF1)
          rho      = MPI_DBL_BUF1(1)
          temp     = MPI_DBL_BUF1(2)
          chi1     = MPI_DBL_BUF1(3)
          chi2     = MPI_DBL_BUF1(4)
          chi3     = MPI_DBL_BUF1(5)
          d_els    = MPI_DBL_BUF1(6)
          dt       = MPI_DBL_BUF1(7)
          dim_x    = MPI_DBL_BUF1(8)
          dim_y    = MPI_DBL_BUF1(9)
          dim_z    = MPI_DBL_BUF1(10)
          sigma    = MPI_DBL_BUF1(11)
          comp     = MPI_DBL_BUF1(12)
          sd       = MPI_DBL_BUF1(13)
          max_bond = MPI_DBL_BUF1(14)
          vol_a    = MPI_DBL_BUF1(15)
          vol_b    = MPI_DBL_BUF1(16)
          phi_diff = MPI_DBL_BUF1(17)
      endif

c Calculate some slab-related values:
      slabs   = cluster_size
      slab_z  = dim_z / float(slabs)
      slab_y  = dim_y
      slab_x  = dim_x
      
c Calculate some boxing-scheme/system partition information
      mx = int(dim_x/sd_box)
      my = int(dim_y/sd_box)
      mz = int(dim_z/sd_box)

      ncell = mx*my*mz
      cellx_inv = dfloat(mx)/dim_x
      celly_inv = dfloat(my)/dim_y
      cellz_inv = dfloat(mz)/dim_z

c Some DPD specfic values:
      gamma1 = sigma**2/(2.d0*temp)
      sigma  = sigma*dsqrt(3.d0/dt)

      return
      end

c
c bf_entry() - Entry point for binary fluid simulation.
c
      subroutine bf_entry()
      implicit none
      integer i, j, k, h, g
      double precision sd_box
      parameter(sd_box=1.d0)
      double precision px, py, pz
      double precision pxx, pyy, pzz
      double precision tenergy, tenergy_sys
      double precision stress_tmp
      double precision bf_sys_temp, t_temp
      common/temperature/t_temp
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'

c Create a configuration
      if (CONFIG_TYPE .eq. DISORDERED) then
         call bf_place_polymers()
      elseif (CONFIG_TYPE .eq. PHASE_SEPARATED) then
         call bf_place_polymers_slab()
      endif

c Transfer particles outside the slab to the appropriate neighbor:
       if (cluster_size .gt. 1) then
          call distribute3d()
          call distribute_halo_down3d()
       endif


c Calculate initial number of simulated particles.
       call bf_particle_check()

       ! Do we load a checkpoint?
       if (load_cp .eq. 1 .and. CP_LOAD_TYPE .eq. LOAD_AT_RELAX) then
          call bf_read_checkpoint_master_level_bin()
       else if (load_cp .eq. 1 .and. 
     >          CP_LOAD_TYPE .eq. LOAD_AT_QUENCH) then
           goto 200
       else
          itime0 = 0

          ! Calculate initial forces:
          call place3d()
          call bf_nhood_force_dpd(1)
       endif
       
       call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

c Equilibration loop:
       do itime=1+itime0, equil_time
 
         ! Process system.
         call bf_relax_dpd(1)

         if (mod(itime,calc_freq) .eq. 0) then
         
            ! Check to see that linear momentum is conserved.
            call bf_sys_momentum(px, py, pz)

            t_temp = bf_sys_temp()

            if (cluster_rank .eq. MASTER) then
               ! Screen output
               write(*,1258) itime, t_temp, px, py, pz

               open(12, file="temperature.data", access="append")
               write(12,*) itime, t_temp
               close(12)

               open(12, file="momentum.data", access="append")
               write(12,*) itime, px, py, pz
               close(12)
            endif
            
            ! Radii of gyration
            call bf_radius_gyration()
               
            ! Local velocities
            call bf_local_velocity(20)
         endif
      enddo

      ! Write checkpoint
      ! write(*,*) 'Turn checkpoint back on.'
      call bf_write_checkpoint_master_level_bin()

c Now that the system has relaxed, we proceed to process the system
c while considered particle species/spin/color
       if (cluster_rank .eq. MASTER) then
          write(*,*) ' '
          write(*,*) 'Quenching system...'
          write(*,*) ' '
       endif

c LOAD_AT_QUENCH entry point
200    continue

       ! Do we load a checkpoint?
       if (load_cp .eq. 1 .and. CP_LOAD_TYPE .eq. LOAD_AT_QUENCH) then
          call bf_read_checkpoint_master_level_bin()
       else
          itime0 = equil_time
       endif

       call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
       do itime=1+itime0, time+equil_time
         ! Process system.
         call bf_relax_dpd(0)
       
         ! Process stresses  
         call MPI_REDUCE(sigma_xy,stress_tmp,1,MPI_DOUBLE_PRECISION,
     >                   MPI_SUM, MASTER,MPI_COMM_WORLD,mpi_err)
         if (cluster_rank .eq. MASTER) then 
            sigma_xy = -stress_tmp/(dim_x*dim_y*dim_z)
         endif

         call MPI_REDUCE(sigma_xz,stress_tmp,1,MPI_DOUBLE_PRECISION,
     >                   MPI_SUM, MASTER,MPI_COMM_WORLD,mpi_err)
         if (cluster_rank .eq. MASTER) then 
            sigma_xz = -stress_tmp/(dim_x*dim_y*dim_z)
         endif

         call MPI_REDUCE(sigma_yz,stress_tmp,1,MPI_DOUBLE_PRECISION,
     >                   MPI_SUM, MASTER,MPI_COMM_WORLD,mpi_err)
         if (cluster_rank .eq. MASTER) then 
            sigma_yz = -stress_tmp/(dim_x*dim_y*dim_z)
         endif

         if (cluster_rank .eq. MASTER) then
             open(12, file="stress.data", access="append")
             write(12,*) (itime-itime0)*dt, sigma_xy, sigma_xz, sigma_yz
             close(12)
         endif

         if (mod(itime,calc_freq) .eq. 0) then

            ! Rouse Mode analysis
            call gen_calc_rouse_modes(1, (itime-itime0)/calc_freq)

            ! Collect particles at master level for snapshots and bridges and all of that.
            call collect_master_level3d()
            
            ! Check to see that linear momentum is conserved.
            call bf_sys_momentum(px, py, pz)

            ! Check system temperature.
            t_temp = bf_sys_temp()

            if (cluster_rank .eq. MASTER) then
               ! Screen output
               write(*,1258) itime, t_temp, px, py, pz
               
               ! System temperature
               open(12, file="temperature.data", access="append")
               write(12,*) itime, t_temp
               close(12)

               ! Momentum
               open(12, file="momentum.data", access="append")
               write(12,*) itime, px, py, pz
               close(12)
 
               ! Finalize bridge detection
               call bf_finalize_bridges()
            endif
        
            ! Find clusters.
            call bf_find_clusters()
   
            ! Total energy
            call bf_finalize_energy()

            ! Radii of gyration 
            call bf_radius_gyration()
            
            ! Structure factor
            call bf_comp_stfac()
            
            ! Local velocity
            call bf_local_velocity(20)

            ! Density profiles
            call bf_density_profiles()
         endif
         
         if (mod(itime, check_freq) .eq. 0) then
            ! Write checkpoint. 
            call collect_master_level3d()
            call bf_write_out_master_level()
            call bf_write_snapshot_master_level()
            call bf_write_checkpoint_master_level_bin()

            !call collect_master_level3d_poly()
            !call bf_write_out_master_level_poly()
         endif
         
      enddo

      ! Finalize and clean up.
      call MPI_FINALIZE(mpi_err)

1258  format('--> ', I7,' ', E20.9, ' ' E9.2, ' ', E9.2, ' ', E9.2) 
      stop
      end

c
c This subroutine is in the event something goes wrong and the program needs to 
c halt (properly)... such as bond lengths exploding and particles not being processed.
c
      subroutine bf_halt_on_error()
      include 'mpif.h'
      include 'common.f'

      ! Exit.
      if (cluster_rank .eq. MASTER) then
          write(*,*)'*** HALTING ON ERROR!!!! ***'
      endif

      ! Write checkpoint where we are.
      call bf_write_checkpoint_master_level_bin()

      ! Finalize and clean up.
      call MPI_FINALIZE(mpi_err)

      return
      end



c
c bf_intrachain_dpd() - Calculates intrachain interactions in polymers.
c                
      subroutine bf_intrachain_dpd(p_type, r_eq, l_bond)
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      double precision r_time1, r_time2, f_sum, f_comp, max_bond2
      double precision bond_strength, pref_bond, l_bond, r_eq
      integer p_type, i, j, k, l, chain_length
      integer id_i, id_j, lmn, last_particle
      integer species_ok

      ! These two indices are for the stress tensor calculation.
      integer idx,idy

      !! Timing
      r_time1 = MPI_WTIME()

      ! K is the number of interactions processed.
      k = 0

      ! d_els2 = max bond length squared.
      pref_bond     = l_bond
      max_bond2     = (l_bond - r_eq)**2
      bond_strength = chi3

      ! Loop through all potential members of this node (i.e. loop through
      ! all particles
      do i=1, local_no+ext_recv_count
         ! Match the species:
         if (species(i) .eq. p_type) then
            ! Get neighboring monomer, which is by definition the absolute index 
            ! of i + 1. Furthermore, get its LOCAL index.
            id_i = i
            id_j = members(logical_id(id_i)+1)

            ! Verify that these particles are part of the same chain.
            if (p_flag(id_i) .eq. p_flag(id_j) .and. id_j .gt. 0 .and.
     >         (id_i .le. local_no .or. id_j .le. local_no)) then
               ! We have an interacting pair.
               k = k + 1
               
               ! Separation distance:
               dx = poly_x(id_i) - poly_x(id_j)
               dy = poly_y(id_i) - poly_y(id_j)
               dz = poly_z(id_i) - poly_z(id_j)

               dr = dsqrt(dx**2 + dy**2 + dz**2)

               ! Force from FENE potential:
               if (BOND_TYPE .eq. FENE) then
                  f_sum = -bond_strength /(1.d0-(dr-r_eq)**2/max_bond2)
     >                     * (1.d0 - r_eq/dr)
               else if (BOND_TYPE .eq. HARMONIC) then
                  f_sum = -bond_strength * (1.d0 - pref_bond/dr)
                  utot  = utot + 0.5d0 * bond_strength * dr
               endif

               ! Add components to force arrays:
               f_comp   = f_sum*dx
               fx(id_i) = fx(id_i) + f_comp
               fx(id_j) = fx(id_j) - f_comp

               f_comp   = f_sum*dy
               fy(id_i) = fy(id_i) + f_comp
               fy(id_j) = fy(id_j) - f_comp

               f_comp   = f_sum*dz
               fz(id_i) = fz(id_i) + f_comp
               fz(id_j) = fz(id_j) - f_comp
 
               ! Stress calculation
               sigma_xy = sigma_xy + dx*f_sum*dy
               sigma_xz = sigma_xz + dx*f_sum*dz
               sigma_yz = sigma_yz + dy*f_sum*dz

               if (p_type .eq. 1) then 
                  poly_bond_a = poly_bond_a + dr
               else if (p_type .eq. -1) then 
                  poly_bond_b = poly_bond_b + dr
               else if (p_type .ge. 5) then
                  poly_bond_g = poly_bond_g + dr
               endif
            endif
         endif
         species_ok = 0
      enddo

      !! Timing
      !r_time2 = MPI_WTIME()

      ! Did we miss any particles?
      !call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
      !call MPI_REDUCE(k,j,1,MPI_INTEGER,MPI_SUM,MASTER,MPI_COMM_WORLD,
      !>     mpi_err)
      !if (cluster_rank .eq. MASTER) then
      !   if (p_type .eq. 1) then
      !     j = (chain_a-1)*cluster_size*a_poly-j
      !   else
      !     j = (chain_b-1)*cluster_size*b_poly-j
      !   endif
      !   if (j .ne. 0) stop
      !   !write(*,*) 'intrachain() takes: ', r_time2 - r_time1
      !endif
      return 
      end

c
c bf_entanglement() - Calculates entanglement interactions between
c                      macrocycles and their backbone chains (i.e.,
c                      does not allow bond crossing)
c                
c Mesoscopic simulation of entanglements using dissipative particle
c dynamics: Application to polymer brushes
c
c Florent Goujon, Patrice Malfreyt, and Dominic J. Tildesley
c
c Citation: The Journal of Chemical Physics 129, 034902 (2008); doi:
c           C10.1063/1.2954022
c
      subroutine bf_entanglement()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      double precision r_time1, r_time2, f_sum, f_comp, max_bond2
      double precision bond_strength, pref_bond
      integer p_type, i, j, k, l, chain_length
      integer id_i, id_j, lmn, last_particle
      integer id_k, id_l, id_start
      integer species_ok

      double precision Px_i, Py_i, Pz_i
      double precision Px_j, Py_j, Pz_j
      double precision dx1, dy1, dz1
      double precision t_i, t_j
      double precision dx_e, dy_e, dz_e, dr_e
      double precision denom
      double precision numer_x, numer_y, numer_z
      double precision Pi_Pj_diff_x
      double precision Pi_Pj_diff_y
      double precision Pi_Pj_diff_z
      double precision Rj, Ri
      double precision Rji

      ! Bond lengths and strengths.
      pref_bond     = 0.03
      max_bond2     = 0.03**2
      bond_strength = 80.d0

      ! Loop through all potential members of this node (i.e. loop through
      ! all particles
      do i=1, local_no+ext_recv_count
         ! Get neighboring monomer, which is by definition the absolute index 
         ! of i + 1. Furthermore, get its LOCAL index.
            id_i = i
            id_j = members(logical_id(id_i)+1)

            ! We need to enture id_k and id_l are on the 
            ! are on the same polymer backbone.
            if (p_flag(id_i) .ne. p_flag(id_j)) then 
                goto 100
            endif

            ! Separation distance.
            ! This defines the bond.
            dx = poly_x(id_i) - poly_x(id_j)
            dy = poly_y(id_i) - poly_y(id_j)
            dz = poly_z(id_i) - poly_z(id_j)

            ! This defines the halfway point of the bond.
            Px_i = 0.5d0*(x(id_i) + x(id_j))
            if (Px_i .gt. dim_x/2.d0) Px_i = Px_i - dim_x
            Py_i = 0.5d0*(y(id_i) + y(id_j))
            if (Py_i .gt. dim_y/2.d0) Py_i = Py_i - dim_y
            Pz_i = 0.5d0*(z(id_i) + z(id_j))
            if (Pz_i .gt. dim_z/2.d0) Pz_i = Pz_i - dim_z

            ! Loop through remaining particles on this code and
            ! compute entanglement force if necessary.
            do k=i, local_no+ext_recv_count
                  id_k     = k
                  id_l     = members(logical_id(id_k)+1)

                  ! We need to enture id_k and id_l are on the 
                  ! are on the same polymer backbone.
                  if (p_flag(id_k) .ne. p_flag(id_l)) then
                        goto 101
                  endif

                  Px_j = 0.5d0*(x(id_k) + x(id_l))
                  if (Px_j .gt. dim_x/2.d0) Px_j = Px_j - dim_x
                  Py_j = 0.5d0*(y(id_k) + y(id_l))
                  if (Py_j .gt. dim_y/2.d0) Py_j = Py_j - dim_y
                  Pz_j = 0.5d0*(z(id_k) + z(id_l))
                  if (Pz_j .gt. dim_z/2.d0) Pz_j = Pz_j - dim_z

                  dx1 = poly_x(id_k) - poly_x(id_l)
                  dy1 = poly_y(id_k) - poly_y(id_l)
                  dz1 = poly_z(id_k) - poly_z(id_l)

                  Rji = dx*dx1 + dy*dy1 + dz*dz1
                  Rj  = dx1**2 + dy1**2 + dz1**2
                  Ri  = dx**2  + dy**2  + dz**2

                  denom = Rji**2 - Ri*Rj
                  
                  Pi_Pj_diff_x = Px_i - Px_j
                  Pi_Pj_diff_y = Py_i - Py_j
                  Pi_Pj_diff_z = Pz_i - Pz_j

                  numer_x = Rj*dx - Rji*dx1
                  numer_y = Rj*dy - Rji*dy1
                  numer_z = Rj*dz - Rji*dz1

                  t_i = Pi_Pj_diff_x*numer_x +
     >                  Pi_Pj_diff_y*numer_y + 
     >                  Pi_Pj_diff_z*numer_z
                  t_i = t_i/denom

                  if (t_i .gt. 0.5) t_i = 0.5
                  if (t_i .lt. -0.5) t_i = -0.5

                  Pi_Pj_diff_x = -Pi_Pj_diff_x
                  Pi_Pj_diff_y = -Pi_Pj_diff_y
                  Pi_Pj_diff_z = -Pi_Pj_diff_z

                  numer_x = Ri*dx1 - Rji*dx
                  numer_y = Ri*dy1 - Rji*dy
                  numer_z = Ri*dz1 - Rji*dz

                  t_j = Pi_Pj_diff_x*numer_x +
     >                  Pi_Pj_diff_y*numer_y + 
     >                  Pi_Pj_diff_z*numer_z
                  t_j = t_j/denom

                  if (t_j .gt. 0.5) t_j = 0.5
                  if (t_j .lt. -0.5) t_j = -0.5

                  dx_e = Px_i + t_i*dx - (Px_j + t_j*dx1)
                  dy_e = Py_i + t_i*dy - (Py_j + t_j*dy1)
                  dz_e = Pz_i + t_i*dz - (Pz_j + t_j*dz1)
                  dr_e = dx_e**2 + dy_e**2 + dz_e**2

                  if (dr_e .lt. max_bond2) then
                      f_sum = bond_strength * (1.d0 - dr_e/pref_bond)

                      f_comp = f_sum * dx_e
                      fx(id_i)         = fx(id_i) + f_comp
                      fx(id_j)         = fx(id_j) + f_comp
                      fx(id_start+k)   = fx(id_start+k)   - f_comp 
                      fx(id_start+k-1) = fx(id_start+k-1) - f_comp 

                      f_comp = f_sum * dy_e
                      fy(id_i)         = fy(id_i) + f_comp
                      fy(id_j)         = fy(id_j) + f_comp
                      fy(id_start+k)   = fy(id_start+k)   - f_comp 
                      fy(id_start+k-1) = fy(id_start+k-1) - f_comp 

                      f_comp = f_sum * dz_e
                      fz(id_i)         = fz(id_i) + f_comp
                      fz(id_j)         = fz(id_j) + f_comp
                      fz(id_start+k)   = fz(id_start+k)   - f_comp 
                      fx(id_start+k-1) = fz(id_start+k-1) - f_comp 

                      ! Stress calculation
                      !sigma_xy = sigma_xy + dx*f_sum*dy
                      !sigma_xz = sigma_xz + dx*f_sum*dz
                      !sigma_yz = sigma_yz + dy*f_sum*dz
                  endif
101               continue
            enddo
100         continue
      enddo
      return
      end

c
c bf_torsion_dpd() - Calculates the three-body interaction between
c                    beads to enforce bond angles.
c
c Implicit-solvent mesoscale model based on soft-core potentials for
c self-assembled lipid membranes
c
c Joel D. Revalee, Mohamed Laradji, and P. B. Sunil Kumar
c
c Citation: The Journal of Chemical Physics 128, 035102 (2008); doi:
c C10.1063/1.2825300
c                
      subroutine bf_torsion_dpd(p_type, theta_i)
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      double precision r_time1, r_time2, f_sum, f_comp, max_bond2
      integer p_type, i, j, k, l, chain_length
      double precision k_bend, theta, theta_0, theta_i
      integer id_i, id_j, id_0, lmn, last_particle
      integer species_ok
      double precision rm1_x, rm1_y, rm1_z
      double precision rp1_x, rp1_y, rp1_z
      

      theta_0 = theta_i * pi/180.d0 ! in radians.
      k_bend = 20.d0

      write(*,*) 'This subroutine is broken. Please fix.'
      call bf_halt_on_error()

      ! Loop through all potential members of this node (i.e. loop through
      ! all particles
      do i=1, local_no+ext_recv_count
         ! Match the species:
         ! Note that the second clause (= 9) is to include the stopper
         ! group from the rotaxane module. It can be ignored for other
         ! simulation modules.
         if (species(i) .eq. p_type .or. abs(species(i)) .eq. p_type
     >       .and. abs(p_type) .eq. 2) then
            ! Get neighboring monomer, which is by definition the absolute index 
            ! of i + 1. Furthermore, get its LOCAL index.
            id_0 = i
            id_i = members(logical_id(id_i)-1)
            id_j = members(logical_id(id_i)+1)

            ! Verify that these particles are part of the same chain.
            if (p_flag(id_i) .eq. p_flag(id_j) .and. id_j .gt. 0 .and.
     >         (id_i .le. local_no .or. id_j .le. local_no) .and.
     >          species(id_i) .eq. species(id_j)) then

               ! We have an interacting pair.
               ! Separation distances:
               rm1_x = poly_x(id_0) - poly_x(id_i)
               rm1_y = poly_y(id_0) - poly_y(id_i)
               rm1_z = poly_z(id_0) - poly_z(id_i)

               rm1_x = poly_x(id_j) - poly_x(id_0)
               rm1_y = poly_y(id_j) - poly_y(id_0)
               rm1_z = poly_z(id_j) - poly_z(id_0)

               dx = poly_x(id_i) - poly_x(id_j)
               dy = poly_y(id_i) - poly_y(id_j)
               dz = poly_z(id_i) - poly_z(id_j)

               theta = (rm1_x*rp1_x + rm1_y*rp1_y + rm1_z*rp1_z) /
     >                 dsqrt(rm1_x**2+rm1_y**2+rm1_z**2)         /
     >                 dsqrt(rp1_x**2+rp1_y**2+rp1_z**2)

               f_sum = -k_bend * (cos(theta_0) - cos(theta))

               ! Add components to force arrays:
               f_comp   = f_sum*dx
               fx(id_i) = fx(id_i) + f_comp
               fx(id_j) = fx(id_j) - f_comp

               f_comp   = f_sum*dy
               fy(id_i) = fy(id_i) + f_comp
               fy(id_j) = fy(id_j) - f_comp

               f_comp   = f_sum*dz
               fz(id_i) = fz(id_i) + f_comp
               fz(id_j) = fz(id_j) - f_comp

               ! Stress calculation
               !sigma_xy = sigma_xy + dx*f_sum*dy
               !sigma_xz = sigma_xz + dx*f_sum*dz
               !sigma_yz = sigma_yz + dy*f_sum*dz
            endif
         endif
      enddo
      return 
      end
c
c bf_nhood_force_dpd() - Calculate pair-wise forces between all DPD particles.
c
      subroutine bf_nhood_force_dpd(relax)
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer ix, iy, iz, jx, jy, jz, mxy, idx_i, idx_j
      integer icell, i, j, k, ii, jj, jcell, kcell
      integer relax

      ! These indices for stress tensor
      integer idx, idy, interac_count, interactions

      ! Periodic boundary condition for XY
      double precision x_trans, y_trans, z_trans

      ! Temporary variables
      double precision f_sum, r_comp, f_comp
      double precision f_conservative, f_random, f_dissipative
      double precision dbl_tmp

      ! Timing
      double precision r_time1, r_time2


      ! Preliminary calculations:
      mxy  = mx*my
      utot = 0.d0
      !call bf_init_stress()


      sigma_xy = 0.d0
      sigma_xz = 0.d0  
      sigma_yz = 0.d0

      interac_count = 0
      interactions  = 0

      ! Empty force arrays
      do i=1, local_no+ext_recv_count
         fx(i) = 0.d0 
         fy(i) = 0.d0
         fz(i) = 0.d0
      enddo

      ! Initialize polymer bond values
      poly_bond_a = 0.d0
      poly_bond_b = 0.d0
      poly_bond_g = 0.d0

      !! Timing
      r_time1 = MPI_WTIME()
      dbl_tmp = 0.d0

      ! Do intrachain interactions:
      if (chain_a .gt. 1) then
         ! Do all intrachain interactions in all 'a' polymers:
         call bf_intrachain_dpd(1,d_els,max_bond)
      endif

      if (chain_b .gt. 1) then
         ! Do all intrachain interactions in all 'b' polymers:
         call bf_intrachain_dpd(-1,d_els,max_bond)
      endif
 
      ! Entanglement fores to prevent chain crossing.
      if (chain_a .gt. 1 .or. chain_b .gt. 1) then
        !call bf_entanglement()
      endif

      ! This is the main loop. Calculate all pair-wise forces for all DPD particles.
      k = 0
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
9999                 if (jcell .eq. icell)  j=list(i)
 
                     do while (j .gt. 0)
                       dx = x(i) - x(j) + x_trans
                       dy = y(i) - y(j) + y_trans
                       dz = z(i) - z(j) + z_trans

                       dr = dx**2 + dy**2 + dz**2

                       ! Are the particles close enough to interact?
                       if (dr .lt. 1.d0) then
                          ! Yes, they are.
                          dr = dsqrt(dr)
                          wr = 1.d0/dr - 1.d0

                          ! Choose interaction strength if we are not relaxing.
                          if (relax .eq. 1) then
                             a = chi1
                          else
                             if (species(i) .eq. species(j)) then
                                a = chi1
                             else
                                a = chi2
                             endif
                          endif

                          if (i .le. local_no .or. j .le. local_no) then
  
                            ! rdv = r . v (r-dot-v)
                            rdv = dx*(vx(i)-vx(j)) + 
     >                            dy*(vy(i)-vy(j)) +
     >                            dz*(vz(i)-vz(j))
                          
                            ! Get the random component used in the random force
                            r_comp = 2.d0*ran3() - 1.d0
                         
                            ! f_sum is the sum of the conservative, random, and dissipative forces
                            ! (per the DPD method)
                            f_conservative = a*wr
                            f_random       = sigma*wr*r_comp
                            f_dissipative  = -gamma1*wr**2*rdv

                            f_sum = f_conservative +
     >                              f_random       +
     >                              f_dissipative

                            ! Add f_sum to each component of the net force on each particle.
                            !! X component of force
                            f_comp = f_sum*dx
                            fx(i)       = fx(i) + f_comp
                            fx(j)       = fx(j) - f_comp

                            !! Y component of force 
                            f_comp = f_sum*dy
                            fy(i)       = fy(i) + f_comp
                            fy(j)       = fy(j) - f_comp

                            !! Z component of force
                            f_comp = f_sum*dz
                            fz(i)       = fz(i) + f_comp
                            fz(j)       = fz(j) - f_comp

                            utot = utot + 0.5d0*a*(dr - 1.d0)**2

                            ! Stress calculation
                            !sigma_xy = sigma_xy + dx*f_sum*dy
                            !sigma_xz = sigma_xz + dx*f_sum*dz
                            !sigma_yz = sigma_yz + dy*f_sum*dz
                          endif ! end if (i .le. local_no .or. j .le. local_no)
                       endif ! end if (dr .lt. 1)
1024                   j = list(j)
                     enddo ! end do-while (j .gt. 0)
                     j = head(jcell)
                     i = list(i)

                     ! If there are more particles, proceed
                     if (i .gt. 0) goto 9999
                  endif ! end (j .gt. 0)                
               enddo ! end kcell loop
            endif ! end (i .gt. 0)
          enddo ! end ix loop
        enddo ! end iy loop
      enddo ! end iz loop

      ! Here we must return the results from processing 'halo' 
      ! particles back up to our neighboring node.
      if (cluster_size .gt. 1) then
         ! Send forces up:
         call distribute_halo_up3d()
      endif

      ! Merge clusters together that cross node boundaries
      if (mod(itime, calc_freq) .eq. 0 .and. relax .eq. 0) then
          call bf_find_bridges()
          !call bf_finalize_stress()
      endif
      
      r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then
      !  write(*,*) 'nhood_force_dpd takes: ', r_time2-r_time1
      !endif

      return
      end


c
c bf_relax_dpd() - Integrates equations of motion according to the Velocity-Verlet algoritm.
c
c
      subroutine bf_relax_dpd(relax)
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer i,idx,j, relax
      double precision dt_2, fxt, fyt, fzt
      double precision del_lam
      double precision r_time1, r_time2

      !! timing
      r_time1 = MPI_WTIME()
 
      ! dt_2 = dt*dt
      dt_2 = dt**2

      ! del_lam = dt * lambda = 0.5*dt
      del_lam = LAMBDA*dt

      do i=1, local_no
        ! Store f_() values to scalars to cut down on array access.
        ! Store forces for correction portion of VV.
        fxt = fx(i)
        fyt = fy(i)
        fzt = fz(i)
       
        ! Store velocities for correction
        vx1(i) = vx(i)
        vy1(i) = vy(i)
        vz1(i) = vz(i)

        ! Update positions 
        x(i) = x(i) + dt*vx(i) + 0.5d0*dt_2*fxt
        y(i) = y(i) + dt*vy(i) + 0.5d0*dt_2*fyt
        z(i) = z(i) + dt*vz(i) + 0.5d0*dt_2*fzt
 
        ! Update positions (polymer)
        poly_x(i) = poly_x(i) + dt*vx(i) + 0.5d0*dt_2*fxt
        poly_y(i) = poly_y(i) + dt*vy(i) + 0.5d0*dt_2*fyt
        poly_z(i) = poly_z(i) + dt*vz(i) + 0.5d0*dt_2*fzt

        ! Predict velocities
        vx(i) = vx(i) + del_lam*fxt
        vy(i) = vy(i) + del_lam*fyt
        vz(i) = vz(i) + del_lam*fzt

        ! Store forces:
        fx1(i) = fxt
        fy1(i) = fyt
        fz1(i) = fzt
      enddo 

      ! Process particles again
      if (cluster_size .gt. 1) then
         call distribute3d()
         call distribute_halo_down3d()
      endif

      call place3d()
      call bf_nhood_force_dpd(relax)

      ! Update velocities
      sigma_xy = 0.d0
      sigma_xz = 0.d0
      sigma_yz = 0.d0
      do i=1, local_no
         vx(i) = vx1(i) + 0.5d0*dt*(fx(i)+fx1(i))
         vy(i) = vy1(i) + 0.5d0*dt*(fy(i)+fy1(i))
         vz(i) = vz1(i) + 0.5d0*dt*(fz(i)+fz1(i))        

         ! Stresses
         !sigma_xy = sigma_xy + vx(i)*vy(i) 
         !sigma_xz = sigma_xz + vx(i)*vz(i)
         !sigma_yz = sigma_yz + vy(i)*vz(i)
         sigma_xy = sigma_xy + x(i)*fy(i)
         sigma_xz = sigma_xz + x(i)*fz(i)
         sigma_yz = sigma_yz + y(i)*fz(i)
      enddo

      !r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then
      !   write(*,*) 'relax_dpd takes: ', r_time2-r_time1
      !endif
      return
      end

c
c bf_place_polymers(): Places polymers randomly into the system,
c                   taking into account system size and chain
c                   length for each type of polymer.

      subroutine bf_place_polymers()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer i, j, k, tmp_no,lmn
      double precision vx_total, vy_total, vz_total 
      double precision vx_tmp,   vy_tmp,   vz_tmp
      double precision dtmp, theta, phi
      double precision tmp1, tmp2, tmp3
      integer p_species_tmp(SUPP_SYS_SIZE), offset

      if (cluster_rank .eq. MASTER) then
          write(*,*) 'Creating disordered configuration...'
      endif

      ! Calculate a_no, the total number of dpd particles of
      ! type a.
      a_no   = int(slab_x*slab_y*slab_z*(phi_diff+1)/(2*vol_a))
      a_poly = int(a_no/chain_a)
      a_no   = a_poly * chain_a

      ! Similarly, calculate b_no.
      b_no   = int((slab_x*slab_y*slab_z-a_no*vol_a)/vol_b)
      b_poly = int(b_no/chain_b)
      b_no   = b_poly * chain_b

      ! Total number of polymers
      local_no = a_no + b_no
      n_poly = a_poly+b_poly

      ! IF we're master, write simulation details to disk.
      if (cluster_rank .eq. MASTER) then
         open (12, file='sim_details.data')
         write(12,*) 'Simulation details and parameters:'
         write(12,*) ' '
         write(12,*) 'DPD_NO: ', local_no*cluster_size
         write(12,*) 'A_NO  : ', a_no*cluster_size
         write(12,*) 'B_NO  : ', b_no*cluster_size
         write(12,*) 'A_POLY: ', a_poly*cluster_size
         write(12,*) 'B_POLY: ', b_poly*cluster_size
         write(12,*) '---------'
         write(12,*) ' '
         close(12)
      endif

      vx_total = 0.d0
      vy_total = 0.d0
      vz_total = 0.d0

      do i=1, cluster_size*local_no
         members(i) = 0
      enddo

      do i=1, cluster_size*n_poly
         p_species(i)     = 0
         p_species_tmp(i) = 0
      enddo

      ! Place all a particles into the system.
      ! (We must do this separately because of the varying 
      !  chain lengths between A and B particles)    
      do i=1, a_poly
          ! k is the monomer number.
          k = (i-1)*chain_a + 1         

          ! Start polymer at a random point in the system.
          x(k) = ran3()*slab_x
          y(k) = ran3()*slab_y
          z(k) = ran3()*slab_z + cluster_rank*slab_z

          ! Species
          species(k)    = 1

          ! Logical ID -- number of particle within entire system
          logical_id(k) = int(cluster_rank*local_no + k)

          ! Logical -> Physical Mapping
          members(logical_id(k)) = k

          ! Set poly arrays
          poly_x(k) = x(k)
          poly_y(k) = y(k)
          poly_z(k) = z(k)

          ! Group monomers into chains.
          p_flag(k) = i+cluster_rank*n_poly

          ! Make sure we know the particle type.
          tmp_no = i + cluster_rank*n_poly
          p_species(tmp_no) = 1

          ! Add rest of chain
          do j=1, chain_a-1
              k = k + 1

              phi   = ran3() * two_pi
              theta = ran3() * pi

              dx = d_els*dcos(phi)*dsin(theta)
              dy = d_els*dsin(phi)*dsin(theta)
              dz = d_els*dcos(theta)

              ! Place next monomer some short distance away 
              ! from previous monomer.
              x(k) = x(k-1) + dx
              y(k) = y(k-1) + dy
              z(k) = z(k-1) + dz

              ! Species
              species(k) = 1
            
              ! Logical ID -- number of particle within entire system
              logical_id(k) = cluster_rank*local_no + k

              ! Logical -> Physical Mapping
              members(logical_id(k)) = k

              ! Set poly arrays
              poly_x(k) = x(k)
              poly_y(k) = y(k)
              poly_z(k) = z(k)
 
              ! Group monomers into chains.
              p_flag(k) = i+cluster_rank*n_poly
         enddo
      enddo

      offset = k

      do i=1, b_poly
          ! k is the monomer number.
          k = (i-1)*chain_b + offset + 1       

          ! Start polymer at a random point in the system.
          x(k) = ran3()*slab_x
          y(k) = ran3()*slab_y
          z(k) = ran3()*slab_z + cluster_rank*slab_z

          ! Species
          species(k) = -1

          ! Logical ID -- number of particle within entire system
          logical_id(k) = cluster_rank*local_no + k

          ! Logical -> Physical Mapping
          members(logical_id(k)) = k

          ! Set poly arrays
          poly_x(k) = x(k)
          poly_y(k) = y(k)
          poly_z(k) = z(k)

          ! Group monomers into chains.
          p_flag(k) = i+cluster_rank*n_poly+a_poly

          ! Make sure we know the particle type.
          tmp_no = i + a_poly + cluster_rank*n_poly
          p_species(tmp_no) = -1

          ! Add rest of chain
          do j=1, chain_b-1
              k = k + 1

              phi   = ran3() * two_pi
              theta = ran3() * pi

              dx = d_els*dcos(phi)*dsin(theta)
              dy = d_els*dsin(phi)*dsin(theta)
              dz = d_els*dcos(theta)

              ! Place next monomer some short distance away 
              ! from previous monomer.
              x(k) = x(k-1) + dx
              y(k) = y(k-1) + dy
              z(k) = z(k-1) + dz

              ! Species
              species(k) = -1

              ! Logical ID -- number of particle within entire system
              logical_id(k) = cluster_rank*local_no + k

              ! Logical -> Physical Mapping
              members(logical_id(k)) = k

              ! Set poly arrays
              poly_x(k) = x(k)
              poly_y(k) = y(k)
              poly_z(k) = z(k)

              ! Group monomers into chains.
              p_flag(k) = i+cluster_rank*n_poly+a_poly
          enddo
      enddo

      ! Assign velocities to all particles.
      do i=1, local_no
          theta = pi*ran3()
          phi   = two_pi*ran3()

          ! Assign sqrt(temp) as a velocity in some random
          ! directon.
          vx(i) = dsqrt(temp) * dcos(phi) * dsin(theta)
          vy(i) = dsqrt(temp) * dsin(phi) * dsin(theta)
          vz(i) = dsqrt(temp) * dcos(theta)

          ! Add velocity to total velocities (for scaling below)
          vx_total = vx_total + vx(i)
          vy_total = vy_total + vy(i)
          vz_total = vz_total + vz(i)
      enddo

      call MPI_REDUCE(local_no,particle_count,1,MPI_INTEGER,MPI_SUM,
     >                MASTER,MPI_COMM_WORLD,mpi_err)
      call MPI_BCAST(particle_count,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,
     >               mpi_err)

      ! Let master recombine all velocities to get a net drift
      ! velocity.
      call MPI_REDUCE(vx_total, vx_tmp, 1, MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)

      call MPI_REDUCE(vy_total, vy_tmp, 1, MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)

      call MPI_REDUCE(vz_total, vz_tmp, 1, MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)

      ! Assign if MASTER node.
      if (cluster_rank .eq. MASTER) then
         vx_total = vx_tmp
         vy_total = vy_tmp
         vz_total = vz_tmp
      endif

      ! Broadcast new velocity totals:
      call MPI_BCAST(vx_total, 1, MPI_DOUBLE_PRECISION, MASTER,
     >               MPI_COMM_WORLD, mpi_err)

      call MPI_BCAST(vy_total, 1, MPI_DOUBLE_PRECISION, MASTER,
     >               MPI_COMM_WORLD, mpi_err)

      call MPI_BCAST(vz_total, 1, MPI_DOUBLE_PRECISION, MASTER,
     >               MPI_COMM_WORLD, mpi_err)

      ! Divide total velocities by number of particles.
      vx_total = vx_total / particle_count
      vy_total = vy_total / particle_count
      vz_total = vz_total / particle_count

      ! Scale velocities to eliminate a net-drift.
      do i=1, local_no
          vx(i) = vx(i) - vx_total
          vy(i) = vy(i) - vy_total
          vz(i) = vz(i) - vz_total

          ! Just some initialization stuff.
          fx1(i) = 0.d0
          fy1(i) = 0.d0
          fz1(i) = 0.d0
          vx1(i) = 0.d0
          vy1(i) = 0.d0
          vz1(i) = 0.d0
      enddo

      ! Set particle count.
      local_no_orig    = local_no
      ext_send_count = 0
      ext_recv_count = 0

      ! One last thing, we need a mapping of all polymers to a species.
      i = cluster_size*n_poly
      call MPI_REDUCE(p_species, p_species_tmp, i, MPI_INTEGER,
     >                MPI_SUM, MASTER, MPI_COMM_WORLD, mpi_err)

      ! For the master, re-assign to create map.
      if (cluster_rank .eq. MASTER) then
         do i=1, cluster_size*n_poly
            p_species(i) = p_species_tmp(i)
         enddo
      endif
      return
      end

c
c bf_place_polymers_slab(): Places polymers into a completely phase
c                        separated configuration.
c
      subroutine bf_place_polymers_slab()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer i, j, k, tmp_no,lmn
      double precision vx_total, vy_total, vz_total 
      double precision vx_tmp,   vy_tmp,   vz_tmp
      double precision dtmp, theta, phi
      double precision tmp1, tmp2, tmp3
      integer p_species_tmp(SUPP_SYS_SIZE), offset

      if (cluster_rank .eq. MASTER) then
          write(*,*) 'Creating phase separated configuration...'
      endif

      ! Calculate a_no, the total number of dpd particles of
      ! type a.
      a_no   = int(slab_x*slab_y*slab_z*(phi_diff+1)/(2*vol_a))
      a_poly = int(a_no/chain_a)
      a_no   = a_poly * chain_a

      ! Similarly, calculate b_no.
      b_no   = int((slab_x*slab_y*slab_z-a_no*vol_a)/vol_b)
      b_poly = int(b_no/chain_b)
      b_no   = b_poly * chain_b

      ! Total number of polymers
      local_no = a_no + b_no
      n_poly = a_poly+b_poly

      ! IF we're master, write simulation details to disk.
      if (cluster_rank .eq. MASTER) then
         open (12, file='sim_details.data')
         write(12,*) 'Simulation details and parameters:'
         write(12,*) ' '
         write(12,*) 'DPD_NO: ', local_no*cluster_size
         write(12,*) 'A_NO  : ', a_no*cluster_size
         write(12,*) 'B_NO  : ', b_no*cluster_size
         write(12,*) 'A_POLY: ', a_poly*cluster_size
         write(12,*) 'B_POLY: ', b_poly*cluster_size
         write(12,*) '---------'
         write(12,*) ' '
         close(12)
      endif

      vx_total = 0.d0
      vy_total = 0.d0
      vz_total = 0.d0

      do i=1, cluster_size*local_no
         members(i) = 0
      enddo

      do i=1, cluster_size*n_poly
         p_species(i)     = 0
         p_species_tmp(i) = 0
      enddo

      ! Place all a particles into the system.
      ! (We must do this separately because of the varying 
      !  chain lengths between A and B particles)    
      do i=1, a_poly
          ! k is the monomer number.
          k = (i-1)*chain_a + 1         

          ! Start polymer at a random point in the system.
          x(k) = ran3()*(slab_x/2.d0)
          y(k) = ran3()*slab_y
          z(k) = ran3()*slab_z + cluster_rank*slab_z

          ! Species
          species(k)    = 1

          ! Logical ID -- number of particle within entire system
          logical_id(k) = int(cluster_rank*local_no + k)

          ! Logical -> Physical Mapping
          members(logical_id(k)) = k

          ! Set poly arrays
          poly_x(k) = x(k)
          poly_y(k) = y(k)
          poly_z(k) = z(k)

          ! Group monomers into chains.
          if (chain_b .gt. 1) then
            p_flag(k) = i+cluster_rank*n_poly
          else
            p_flag(k) = 1+i+cluster_rank*a_poly
          endif

          ! Make sure we know the particle type.
          tmp_no = i + cluster_rank*n_poly
          p_species(tmp_no) = 1

          ! Add rest of chain
          do j=1, chain_a-1
              k = k + 1

              phi   = ran3() * two_pi
              theta = ran3() * pi

              dx = d_els*dcos(phi)*dsin(theta)
              dy = d_els*dsin(phi)*dsin(theta)
              dz = d_els*dcos(theta)

              ! Place next monomer some short distance away 
              ! from previous monomer.
              x(k) = x(k-1) + dx
              y(k) = y(k-1) + dy
              z(k) = z(k-1) + dz

              ! Species
              species(k) = 1
            
              ! Logical ID -- number of particle within entire system
              logical_id(k) = cluster_rank*local_no + k

              ! Logical -> Physical Mapping
              members(logical_id(k)) = k

              ! Set poly arrays
              poly_x(k) = x(k)
              poly_y(k) = y(k)
              poly_z(k) = z(k)
 
              ! Group monomers into chains.
              if (chain_b .gt. 1) then
                p_flag(k) = i+cluster_rank*n_poly
              else
                p_flag(k) = 1+i+cluster_rank*a_poly
              endif
         enddo
      enddo

      offset = k

      do i=1, b_poly
          ! k is the monomer number.
          k = (i-1)*chain_b + offset + 1       

          ! Start polymer at a random point in the system.
          x(k) = ran3()*(slab_x/2.d0) + slab_x/2.d0
          y(k) = ran3()*slab_y
          z(k) = ran3()*slab_z + cluster_rank*slab_z

          ! Species
          species(k) = -1

          ! Logical ID -- number of particle within entire system
          logical_id(k) = cluster_rank*local_no + k

          ! Logical -> Physical Mapping
          members(logical_id(k)) = k

          ! Set poly arrays
          poly_x(k) = x(k)
          poly_y(k) = y(k)
          poly_z(k) = z(k)

          ! Group monomers into chains.
          if (chain_b .gt. 1) then
             p_flag(k) = i+cluster_rank*n_poly+a_poly
          else
             p_flag(k) = 1
          endif

          ! Make sure we know the particle type.
          tmp_no = i + a_poly + cluster_rank*n_poly
          p_species(tmp_no) = -1

          ! Add rest of chain
          do j=1, chain_b-1
              k = k + 1

              phi   = ran3() * two_pi
              theta = ran3() * pi

              dx = d_els*dcos(phi)*dsin(theta)
              dy = d_els*dsin(phi)*dsin(theta)
              dz = d_els*dcos(theta)

              ! Place next monomer some short distance away 
              ! from previous monomer.
              x(k) = x(k-1) + dx
              y(k) = y(k-1) + dy
              z(k) = z(k-1) + dz

              ! Species
              species(k) = -1

              ! Logical ID -- number of particle within entire system
              logical_id(k) = cluster_rank*local_no + k

              ! Logical -> Physical Mapping
              members(logical_id(k)) = k

              ! Set poly arrays
              poly_x(k) = x(k)
              poly_y(k) = y(k)
              poly_z(k) = z(k)

              ! Group monomers into chains.
              if (chain_b .gt. 1) then
                 p_flag(k) = i+cluster_rank*n_poly+a_poly
              else
                 p_flag(k) = 1
              endif
          enddo
      enddo

      ! Assign velocities to all particles.
      do i=1, local_no
          theta = pi*ran3()
          phi   = two_pi*ran3()

          ! Assign sqrt(temp) as a velocity in some random
          ! directon.
          vx(i) = dsqrt(temp) * dcos(phi) * dsin(theta)
          vy(i) = dsqrt(temp) * dsin(phi) * dsin(theta)
          vz(i) = dsqrt(temp) * dcos(theta)

          ! Add velocity to total velocities (for scaling below)
          vx_total = vx_total + vx(i)
          vy_total = vy_total + vy(i)
          vz_total = vz_total + vz(i)
      enddo

      ! Translate the x-coordinates of all particles upward by 10 units.
      do i=1, local_no
         x(i)      = x(i)      + 10.d0
         poly_x(i) = poly_x(i) + 10.d0
      enddo

      ! Let master recombine all velocities to get a net drift
      ! velocity.
      call MPI_REDUCE(vx_total, vx_tmp, 1, MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)

      call MPI_REDUCE(vy_total, vy_tmp, 1, MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)

      call MPI_REDUCE(vz_total, vz_tmp, 1, MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)

      ! Assign if MASTER node.
      if (cluster_rank .eq. MASTER) then
         vx_total = vx_tmp
         vy_total = vy_tmp
         vz_total = vz_tmp
      endif

      ! Broadcast new velocity totals:
      call MPI_BCAST(vx_total, 1, MPI_DOUBLE_PRECISION, MASTER,
     >               MPI_COMM_WORLD, mpi_err)

      call MPI_BCAST(vy_total, 1, MPI_DOUBLE_PRECISION, MASTER,
     >               MPI_COMM_WORLD, mpi_err)

      call MPI_BCAST(vz_total, 1, MPI_DOUBLE_PRECISION, MASTER,
     >               MPI_COMM_WORLD, mpi_err)

      ! Divide total velocities by number of particles.
      vx_total = vx_total / (cluster_size*local_no)
      vy_total = vy_total / (cluster_size*local_no)
      vz_total = vz_total / (cluster_size*local_no)

      ! Scale velocities to eliminate a net-drift.
      do i=1, local_no
          vx(i) = vx(i) - vx_total
          vy(i) = vy(i) - vy_total
          vz(i) = vz(i) - vz_total

          ! Just some initialization stuff.
          fx1(i) = 0.d0
          fy1(i) = 0.d0
          fz1(i) = 0.d0
          vx1(i) = 0.d0
          vy1(i) = 0.d0
          vz1(i) = 0.d0
      enddo

      ! Set particle count.
      particle_count = cluster_size*local_no
      local_no_orig    = local_no
      ext_send_count = 0
      ext_recv_count = 0

      ! One last thing, we need a mapping of all polymers to a species.
      i = cluster_size*n_poly
      call MPI_REDUCE(p_species, p_species_tmp, i, MPI_INTEGER,
     >                MPI_SUM, MASTER, MPI_COMM_WORLD, mpi_err)

      ! For the master, re-assign to create map.
      if (cluster_rank .eq. MASTER) then
         do i=1, cluster_size*n_poly
            p_species(i) = p_species_tmp(i)
         enddo
      endif
      return
      end

