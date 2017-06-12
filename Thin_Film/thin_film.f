c////////////////////////////////////////////////////////////
c thin_film.f - Creation and initialization of slab
c               particles. Also handles introduction
c               of slab particles into main arrays
c               after each call of distribute().
c
c Mike Hore, University of Memphis, Summer 2006
c///////////////////////////////////////////////////////////

c
c tfilm_entry() - Entry point for the simulation.
c
      subroutine tfilm_entry()
      implicit none
      integer i, j, k, h, g
      double precision sd_box
      parameter(sd_box=1.d0)
      double precision px, py, pz
      double precision pxx, pyy, pzz
      double precision tenergy, tenergy_sys
      double precision bf_sys_temp, t_temp
      double precision r_time1, r_time2,ttl_time
      double precision MIN_T, MAX_T
      common/temperature/t_temp
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      include 'Thin_Film/thin_film_var.f'

      ttl_time = 0.d0

c Min/max thickness for film calculations. Not necessarily needed.
      if (USE_MIN_MAX .eq. YES) then
         MIN_T = film_min
         MAX_T = film_max
         if (cluster_rank .eq. MASTER) then
            write(*,*) 'Using MIN_T/MAX_T: ', MIN_T, MAX_T
         endif
      else
         MIN_T = wall_thick
         MAX_T = dim_y
         if (cluster_rank .eq. MASTER) then
            write(*,*) 'MIN_T/MAX_T disabled.'
         endif
      endif
   
c Create a configuration
      if (CONFIG_TYPE .eq. DISORDERED) then
         call tfilm_place_polymers()
      else if (CONFIG_TYPE .eq. PHASE_SEPARATED) then
         call tfilm_place_stripe()
      else
         call tfilm_place_wetting()
      endif
      
c Write initial configuration to disk.
      call collect_master_level3d()
      call bf_write_snapshot_master_level()

c Transfer particles outside the slab to the appropriate neighbor:
       if (cluster_size .gt. 1) then
          call distribute3d()
          call distribute_halo_down3d()
       endif

c Calculate initial number of simulated particles.
       call bf_particle_check()

       ! Do we load a checkpoint?
       if (load_cp .eq. 1 .and. CP_LOAD_TYPE .eq. LOAD_AT_RELAX) then
          call bf_read_checkpoint_master_level()
       else if (load_cp .eq. 1 .and. 
     >          CP_LOAD_TYPE .eq. LOAD_AT_QUENCH) then
           goto 200
       else
          itime0 = 0

          ! Calculate initial forces:
          call place3d()
          call tfilm_nhood_force_dpd(1)
       endif
       
       call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

c Equilibration loop:
       do itime=1+itime0, equil_time
 
         ! Process system.
         call tfilm_relax_dpd(1)


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
               
            ! Density profiles
            call tfilm_density_profiles()

            call collect_master_level3d()
            call bf_write_out_master_level()
         endif
      enddo

      ! Write checkpoint
      call bf_write_checkpoint_master_level()

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
          call bf_read_checkpoint_master_level()
       else
          itime0 = equil_time
       endif

       call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
       do itime=1+itime0, time+equil_time
         ! Process system.
         !r_time1 = MPI_WTIME()
         call tfilm_relax_dpd(0)
         !r_time2 = MPI_WTIME()

         !ttl_time = ttl_time + r_time2 - r_time1

         if (mod(itime,calc_freq) .eq. 0) then
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
            call tfilm_find_clusters(MIN_T,MAX_T)
   
            ! Total energy
            call bf_finalize_energy()

            ! Radii of gyration 
            call bf_radius_gyration()
            
            ! Structure factor
            call tfilm_slice_stfac(MIN_T,MAX_T)
            
            ! Density profiles
            call tfilm_density_profiles()
         endif
         
         if (mod(itime, check_freq) .eq. 0) then
            ! Write checkpoint. 
            call collect_master_level3d()
            call bf_write_out_master_level()
            call bf_write_snapshot_master_level()
            call bf_write_checkpoint_master_level()
         endif
         
      enddo

      !ttl_time = ttl_time / time
      !write(*,*) 'Average of ', ttl_time, ' seconds for ', time

      ! Finalize and clean up.
      call MPI_FINALIZE(mpi_err)

1258  format('--> ', I7,' ', E20.9, ' ' E9.2, ' ', E9.2, ' ', E9.2) 
      stop
      end

c
c tfilm_create_slab() - Creates and initializes slab particles.
c
c *The wall spans the XZ plane so that the wall is distributed among
c  all nodes. This is to ensure equal work load among nodes.
c
      subroutine tfilm_create_fslab()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Thin_Film/thin_film_var.f'
      double precision x0, y0, z0, slab_min, slab_max
      double precision xtemp, ytemp, ztemp
      integer i, j, k, n

      ! Lattice constant defined from fluid density.
      lattice_const = 0.2d0 * (4.d0/rho)**(1.d0/3.d0)

      ! Offset variables.
      x0 = 0.5d0 * dim_x
      y0 = 0.5d0 * wall_thick
      z0 = 0.5d0 * dim_z

      ! Number of particles in the fslab.
      n_fslab = 0

      ! Slab_min -- lower boundary of the SLAB.. too many slabs here. Grr!
      ! Good notation is good thinking. Bad notation is just annoying.
      slab_min = cluster_rank * slab_z
      slab_max = slab_min + slab_z

      do i=-400, 400
         do j=-400, 400
            do k=-400, 400

            ! Generate potential particles (FCC lattice)
            xtemp = x0 + lattice_const*(j+k)
            ytemp = y0 + lattice_const*(i+j)
            ztemp = z0 + lattice_const*(i+k)

            ! Does the point fall within the fslab dimensions? 
            ! (Fslab dimensions are dim_x by wall_thick by dim_x)
            if ( ((xtemp .le. dim_x) .and. (xtemp .gt. 0.d0)) .and.
     >           ((ytemp .le. wall_thick) .and. (ytemp .gt. 0.d0)) .and.
     >           ((ztemp .le. slab_max)   .and. (ztemp .gt. slab_min))) 
     >      then
                 n_fslab          = n_fslab + 1
                 fslab_x(n_fslab) = xtemp
                 fslab_y(n_fslab) = ytemp
                 fslab_z(n_fslab) = ztemp
            endif 
            enddo
         enddo
      enddo

      if (cluster_rank .eq. MASTER) then
         write(*,*) 'Created slab substrate.'
      endif
      return
      end


c
c tfilm_read_dist_parameters() - Reads parameters from file and distributes them to
c                             all slave nodes.
      subroutine tfilm_read_dist_parameters()
      implicit none
      double precision sd_box
      parameter(sd_box=1.d0)
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      include 'Thin_Film/thin_film_var.f'

c If we are the master node, read in the configuration
      if (cluster_rank .eq. MASTER) then
          ! Output the program configuration (bond type, etc.)
          write(*,*) 'Parallel/MPI code version: ', code_version
          write(*,*) 'Master node detects ', cluster_size, ' nodes.'
          write(*,*) 'Simulation type is Thin Film.'
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
          call tfilm_read_config()

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
          MPI_DBL_BUF1(18) = chi_wall_a
          MPI_DBL_BUF1(19) = chi_wall_b
          MPI_DBL_BUF1(20) = wall_thick
          MPI_DBL_BUF1(21) = film_min
          MPI_DBL_BUF1(22) = film_max

          ! Broadcast MPI_INT_BUF1()
          call MPI_BCAST(MPI_INT_BUF1,8,MPI_INTEGER,MASTER,
     >                   MPI_COMM_WORLD,mpi_err)

          ! Broadcast MPI_DBL_BUF1()
          call MPI_BCAST(MPI_DBL_BUF1,22,MPI_DOUBLE_PRECISION,
     >                   MASTER,MPI_COMM_WORLD,mpi_err)
          write(*,*) 'Master node has distributed simulation parameters'
          write(*,*) ' '
c Else, wait to receive configuration from the master
      else
          call MPI_BCAST(MPI_INT_BUF1,8,MPI_INTEGER,MASTER,
     >                   MPI_COMM_WORLD,mpi_err)
          
          call MPI_BCAST(MPI_DBL_BUF1,22,MPI_DOUBLE_PRECISION,MASTER,
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
          rho        = MPI_DBL_BUF1(1)
          temp       = MPI_DBL_BUF1(2)
          chi1       = MPI_DBL_BUF1(3)
          chi2       = MPI_DBL_BUF1(4)
          chi3       = MPI_DBL_BUF1(5)
          d_els      = MPI_DBL_BUF1(6)
          dt         = MPI_DBL_BUF1(7)
          dim_x      = MPI_DBL_BUF1(8)
          dim_y      = MPI_DBL_BUF1(9)
          dim_z      = MPI_DBL_BUF1(10)
          sigma      = MPI_DBL_BUF1(11)
          comp       = MPI_DBL_BUF1(12)
          sd         = MPI_DBL_BUF1(13)
          max_bond   = MPI_DBL_BUF1(14)
          vol_a      = MPI_DBL_BUF1(15)
          vol_b      = MPI_DBL_BUF1(16)
          phi_diff   = MPI_DBL_BUF1(17)
          chi_wall_a = MPI_DBL_BUF1(18)
          chi_wall_b = MPI_DBL_BUF1(19)
          wall_thick = MPI_DBL_BUF1(20)
          film_min   = MPI_DBL_BUF1(21)
          film_max   = MPI_DBL_BUF1(22)
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
c tfilm_read_config: Read simlation parameters from 'bfilm_config' in 
c              working directory.
c
      subroutine tfilm_read_config()
      include 'common.f'
      include 'Thin_Film/thin_film_var.f'
      include 'Binary_Fluid/binary_fluid_var.f'

      ! Open file
      open(12, file="tfilm_config")

      ! Read in parameters
      read(12,*) load_cp
      read(12,*) check_freq
      read(12,*) calc_freq
      read(12,*) rho
      read(12,*) chain_a
      read(12,*) chain_b
      read(12,*) time
      read(12,*) equil_time
      read(12,*) temp
      read(12,*) chi1
      read(12,*) chi2
      read(12,*) chi3
      read(12,*) chi_wall_a
      read(12,*) chi_wall_b
      read(12,*) d_els
      read(12,*) max_bond
      read(12,*) dt
      read(12,*) dim_x
      read(12,*) dim_y
      read(12,*) dim_z
      read(12,*) wall_thick
      read(12,*) iseed
      read(12,*) sigma
      read(12,*) comp      
      read(12,*) vol_a
      read(12,*) vol_b
      read(12,*) phi_diff
      read(12,*) sd
      read(12,*) film_min
      read(12,*) film_max

      ! Close file
      close(12)
      return
      end

c
c tfilm_nhood_force_dpd() - Calculate pair-wise forces between all DPD particles.
c
      subroutine tfilm_nhood_force_dpd(relax)
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      include 'Thin_Film/thin_film_var.f'
      integer ix, iy, iz, jx, jy, jz, mxy, idx_i, idx_j
      integer icell, i, j, k, ii, jj, jcell, kcell
      integer relax
      integer itij

      ! These indices for stress tensor
      integer idx, idy, interac_count, interactions

      ! Periodic boundary condition for XY
      double precision x_trans, y_trans, z_trans

      ! Temporary variables
      double precision f_sum, r_comp, f_comp
      double precision f_conservative, f_random, f_dissipative
      double precision dbl_tmp
      double precision wd

      ! Timing
      double precision r_time1, r_time2


      ! Preliminary calculations:
      mxy  = mx*my
      utot = 0.d0
      call bf_init_stress()

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

      !! Timing
      r_time1 = MPI_WTIME()
      dbl_tmp = 0.d0

      ! Do intrachain interactions:
      if (chain_a .gt. 1) then
         ! Do all intrachain interactions in all 'a' polymers:
         call bf_intrachain_dpd(1)
      endif

      if (chain_b .gt. 1) then
         ! Do all intrachain interactions in all 'b' polymers:
         call bf_intrachain_dpd(-1)
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
                          !wd = wr**(1.d0/2.d0)

                          if (species(i) .eq. 0 .and. species(j) .eq. 0)
     >                    then
                             goto 1024
                          else 
                             ! Choose interaction strength if we are not relaxing.
                             if (relax .eq. 1 .and. (species(i) .ne. 0
     >                           .and. species(j) .ne. 0)) then
                                a = chi1
                             else
                                itij = species(i) + species(j)
                                if (species(i) .eq. species(j)) then
                                   a = chi1
                                else if (itij .eq. 1) then
                                   if (relax .eq. 1) then 
                                      a = chi_wall_b
                                   else
                                      a = chi_wall_a
                                   endif
                                else if (itij .eq. -1) then
                                   a = chi_wall_b
                                else
                                   a = chi2
                                endif
                             endif

                             if ((i .le. local_no .or. j .le. local_no)
     >                        .or.
     >                        (species(i) .eq. 0.or. species(j) .eq. 0))
     >                       then
  
                               ! rdv = r . v (r-dot-v)
                               rdv = dx*(vx(i)-vx(j)) + 
     >                               dy*(vy(i)-vy(j)) +
     >                               dz*(vz(i)-vz(j))
                          
                               ! Get the random component used in the random force
                               r_comp = 2.d0*ran3() - 1.d0
                         
                               ! f_sum is the sum of the conservative, random, and dissipative forces
                               ! (per the DPD method)
                               f_conservative = a*wr
                               f_random       = sigma*wr*r_comp
                               f_dissipative  = -gamma1*wr**2*rdv

                               f_sum = f_conservative +
     >                                 f_random       +
     >                                 f_dissipative

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

                               ! Stresses
                               if (mod(itime, calc_freq) .eq. 0) then
                                 call bf_add_stress(dx,dy,dz,wr*a,0)
                               endif
                             endif ! end if (i .le. local_no .or. j .le. local_no)
                          endif ! end if species(i) .ne. 0 .and species(j) .ne. 0
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
          call bf_finalize_stress()
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
      subroutine tfilm_relax_dpd(relax)
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      include 'Thin_Film/thin_film_var.f'

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
        if (species(i) .eq. 0) goto 9979
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
9979    continue
      enddo 

      ! Process particles again
      if (cluster_size .gt. 1) then
         call distribute3d()
         call distribute_halo_down3d()
      endif

      call place3d()
      call tfilm_nhood_force_dpd(relax)

      ! Update velocities
      do i=1, local_no
         if (species(i) .eq. 0) goto 9980
         vx(i) = vx1(i) + 0.5d0*dt*(fx(i)+fx1(i))
         vy(i) = vy1(i) + 0.5d0*dt*(fy(i)+fy1(i))
         vz(i) = vz1(i) + 0.5d0*dt*(fz(i)+fz1(i))        
9980     continue
      enddo

      !r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then
      !   write(*,*) r_time2-r_time1, " seconds per timestep."
      !endif

      return
      end

c
c tfilm_place_polymers(): Places polymers randomly into the system,
c                   taking into account system size and chain
c                   length for each type of polymer.

      subroutine tfilm_place_polymers()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      include 'Thin_Film/thin_film_var.f'
      integer i, j, k, tmp_no,lmn
      double precision vx_total, vy_total, vz_total 
      double precision vx_tmp,   vy_tmp,   vz_tmp
      double precision dtmp, theta, phi
      double precision tmp1, tmp2, tmp3
      integer p_species_tmp(SUPP_SYS_SIZE), offset,id_offset
      integer local_count_tmp(256), local_count(256)
    
      do i=1, 256
         local_count(i)     = 0
         local_count_tmp(i) = 0
      enddo

      if (cluster_rank .eq. MASTER) then
          write(*,*) 'Creating disordered configuration...'
      endif

      ! Calculate a_no, the total number of dpd particles of
      ! type a.
      a_no   = int(slab_x*(slab_y-wall_thick)*slab_z*(phi_diff+1)/
     >             (2*vol_a))
      a_poly = int(a_no/chain_a)

      ! Similarly, calculate b_no.
      b_no   = int((slab_x*(slab_y-wall_thick)*slab_z-a_no*vol_a)/vol_b)
      b_poly = int(b_no/chain_b)

      ! Total number of polymers
      !local_no = a_no + b_no
      local_no  = a_poly*chain_a + b_poly*chain_b
      n_poly = a_poly+b_poly

      ! Find this node's offset, so that all particles have a UNIQUE id. (logical_id)
      local_count_tmp(cluster_rank+1) = local_no+n_fslab

      call MPI_REDUCE(local_count_tmp,local_count,256,MPI_INTEGER,
     >                MPI_SUM,MASTER,MPI_COMM_WORLD,mpi_err)
      call MPI_BCAST(local_count,256,MPI_INTEGER,MASTER,MPI_COMM_WORLD,
     >               mpi_err)
      id_offset = 0
      do i=MASTER,cluster_rank-1
         id_offset = id_offset + local_count(i+1)
      enddo

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
          y(k) = ran3()*(slab_y-wall_thick)+wall_thick
          z(k) = ran3()*slab_z + cluster_rank*slab_z

          ! Species
          species(k)    = 1

          ! Logical ID -- number of particle within entire system
          logical_id(k) = int(id_offset + k)

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
2219          continue
              phi   = ran3() * two_pi
              theta = ran3() * pi

              dx = d_els*dcos(phi)*dsin(theta)
              dy = d_els*dsin(phi)*dsin(theta)
              dz = d_els*dcos(theta)

              ! Place next monomer some short distance away 
              ! from previous monomer.
              x(k) = x(k-1) + dx
              y(k) = y(k-1) + dy
              if (y(k) .lt. wall_thick) goto 2219
              z(k) = z(k-1) + dz

              ! Species
              species(k) = 1
            
              ! Logical ID -- number of particle within entire system
              logical_id(k) = id_offset + k

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
          x(k) = ran3()*slab_x
          y(k) = ran3()*(slab_y-wall_thick)+wall_thick
          z(k) = ran3()*slab_z + cluster_rank*slab_z

          ! Species
          species(k) = -1

          ! Logical ID -- number of particle within entire system
          logical_id(k) = id_offset + k

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
2220          continue
              phi   = ran3() * two_pi
              theta = ran3() * pi

              dx = d_els*dcos(phi)*dsin(theta)
              dy = d_els*dsin(phi)*dsin(theta)
              dz = d_els*dcos(theta)

              ! Place next monomer some short distance away 
              ! from previous monomer.
              x(k) = x(k-1) + dx
              y(k) = y(k-1) + dy
              if (y(k) .lt. wall_thick) goto 2220
              z(k) = z(k-1) + dz

              ! Species
              species(k) = -1

              ! Logical ID -- number of particle within entire system
              logical_id(k) = id_offset + k

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


      ! Add slab particles into the arrays, setting species to 0 and giving them
      ! a negative polymer id.
      do i=1, n_fslab
         k = k + 1 
         x(k) = fslab_x(i)
         y(k) = fslab_y(i)
         z(k) = fslab_z(i)
         vx(k)= 0.d0
         vy(k)= 0.d0
         vz(k)= 0.d0
         poly_x(k) = fslab_x(i)
         poly_y(k) = fslab_y(i)
         poly_z(k) = fslab_z(i)
         species(k) = 0
         logical_id(k) = id_offset+k
         members(logical_id(k)) = k
         p_flag(k) = 0
      enddo

      ! Set particle count.
      local_no = local_no + n_fslab   ! Include slab particles.
      call MPI_REDUCE(local_no,particle_count,1,MPI_INTEGER,
     >                MPI_SUM,MASTER,MPI_COMM_WORLD,mpi_err)
      call MPI_BCAST(particle_count,1,MPI_INTEGER,MASTER,
     >               MPI_COMM_WORLD,mpi_err)


      ! Scale velocities to eliminate a net-drift.
      vx_total = vx_total / particle_count
      vy_total = vy_total / particle_count
      vz_total = vz_total / particle_count
      do i=1, local_no-n_fslab
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

      local_no_orig  = local_no
      ext_send_count = 0
      ext_recv_count = 0

      do i=1, cluster_size*local_no
         members(i) = 0
      enddo

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
c tfilm_place_stripe(): Places polymers randomly into the system with a stripe in the center,
c                       taking into account system size and chain length for each type of polymer.

      subroutine tfilm_place_stripe()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      include 'Thin_Film/thin_film_var.f'
      integer i, j, k, tmp_no,lmn
      double precision vx_total, vy_total, vz_total 
      double precision vx_tmp,   vy_tmp,   vz_tmp
      double precision dtmp, theta, phi
      double precision tmp1, tmp2, tmp3
      integer p_species_tmp(SUPP_SYS_SIZE), offset,id_offset
      integer local_count_tmp(256), local_count(256)
    
      do i=1, 256
         local_count(i)     = 0
         local_count_tmp(i) = 0
      enddo

      if (cluster_rank .eq. MASTER) then
          write(*,*) 'Creating striped configuration...'
      endif

      ! Calculate a_no, the total number of dpd particles of
      ! type a.
      a_no   = int(slab_x*(slab_y-wall_thick)*slab_z*(phi_diff+1)/
     >             (2*vol_a))
      a_poly = int(a_no/chain_a)

      ! Similarly, calculate b_no.
      b_no   = int((slab_x*(slab_y-wall_thick)*slab_z-a_no*vol_a)/vol_b)
      b_poly = int(b_no/chain_b)

      ! Total number of polymers
      !local_no = a_no + b_no
      local_no  = a_poly*chain_a + b_poly*chain_b
      n_poly = a_poly+b_poly

      ! Find this node's offset, so that all particles have a UNIQUE id. (logical_id)
      local_count_tmp(cluster_rank+1) = local_no+n_fslab

      call MPI_REDUCE(local_count_tmp,local_count,256,MPI_INTEGER,
     >                MPI_SUM,MASTER,MPI_COMM_WORLD,mpi_err)
      call MPI_BCAST(local_count,256,MPI_INTEGER,MASTER,MPI_COMM_WORLD,
     >               mpi_err)
      id_offset = 0
      do i=MASTER,cluster_rank-1
         id_offset = id_offset + local_count(i+1)
      enddo

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
          tmp1 = (slab_y-wall_thick)/2.d0
          x(k) = ran3()*slab_x
          y(k) = ran3()*tmp1+wall_thick+tmp1/2.d0
          z(k) = ran3()*slab_z + cluster_rank*slab_z

          ! Species
          species(k)    = 1

          ! Logical ID -- number of particle within entire system
          logical_id(k) = int(id_offset + k)

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
2219          continue
              phi   = ran3() * two_pi
              theta = ran3() * pi

              dx = d_els*dcos(phi)*dsin(theta)
              dy = d_els*dsin(phi)*dsin(theta)
              dz = d_els*dcos(theta)

              ! Place next monomer some short distance away 
              ! from previous monomer.
              x(k) = x(k-1) + dx
              y(k) = y(k-1) + dy
              if (y(k) .lt. wall_thick) goto 2219
              z(k) = z(k-1) + dz

              ! Species
              species(k) = 1
            
              ! Logical ID -- number of particle within entire system
              logical_id(k) = id_offset + k

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
          tmp1 = (slab_y-wall_thick)/4.d0 
          tmp2 = 2*ran3()-1.d0
          
          x(k) = ran3()*slab_x
          if (tmp2 .le. 0.d0) then
             y(k) = ran3()*tmp1+wall_thick
          else
             y(k) = ran3()*tmp1+wall_thick+3*tmp1
          endif
          z(k) = ran3()*slab_z + cluster_rank*slab_z

          ! Species
          species(k) = -1

          ! Logical ID -- number of particle within entire system
          logical_id(k) = id_offset + k

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
2220          continue
              phi   = ran3() * two_pi
              theta = ran3() * pi

              dx = d_els*dcos(phi)*dsin(theta)
              dy = d_els*dsin(phi)*dsin(theta)
              dz = d_els*dcos(theta)

              ! Place next monomer some short distance away 
              ! from previous monomer.
              x(k) = x(k-1) + dx
              y(k) = y(k-1) + dy
              if (y(k) .lt. wall_thick) goto 2220
              z(k) = z(k-1) + dz

              ! Species
              species(k) = -1

              ! Logical ID -- number of particle within entire system
              logical_id(k) = id_offset + k

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


      ! Add slab particles into the arrays, setting species to 0 and giving them
      ! a negative polymer id.
      do i=1, n_fslab
         k = k + 1 
         x(k) = fslab_x(i)
         y(k) = fslab_y(i)
         z(k) = fslab_z(i)
         vx(k)= 0.d0
         vy(k)= 0.d0
         vz(k)= 0.d0
         poly_x(k) = fslab_x(i)
         poly_y(k) = fslab_y(i)
         poly_z(k) = fslab_z(i)
         species(k) = 0
         logical_id(k) = id_offset+k
         members(logical_id(k)) = k
         p_flag(k) = 0
      enddo

      ! Set particle count.
      local_no = local_no + n_fslab   ! Include slab particles.
      call MPI_REDUCE(local_no,particle_count,1,MPI_INTEGER,
     >                MPI_SUM,MASTER,MPI_COMM_WORLD,mpi_err)
      call MPI_BCAST(particle_count,1,MPI_INTEGER,MASTER,
     >               MPI_COMM_WORLD,mpi_err)


      ! Scale velocities to eliminate a net-drift.
      vx_total = vx_total / particle_count
      vy_total = vy_total / particle_count
      vz_total = vz_total / particle_count
      do i=1, local_no-n_fslab
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

      local_no_orig  = local_no
      ext_send_count = 0
      ext_recv_count = 0

      do i=1, cluster_size*local_no
         members(i) = 0
      enddo

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
c tfilm_place_wetting(): Places polymers randomly into the system,
c                        with a slab parallel to the XZ plane. This
c                        is for a study of wetting/pre-wetting and
c                        contact angles.

      subroutine tfilm_place_wetting()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      include 'Thin_Film/thin_film_var.f'
      integer i, j, k, tmp_no,lmn
      double precision vx_total, vy_total, vz_total 
      double precision vx_tmp,   vy_tmp,   vz_tmp
      double precision dtmp, theta, phi
      double precision tmp1, tmp2, tmp3
      integer p_species_tmp(SUPP_SYS_SIZE), offset,id_offset
      integer local_count_tmp(256), local_count(256)
    
      do i=1, 256
         local_count(i)     = 0
         local_count_tmp(i) = 0
      enddo

      if (cluster_rank .eq. MASTER) then
          write(*,*) 'Creating disordered configuration...'
      endif

      ! Calculate a_no, the total number of dpd particles of
      ! type a.
      a_no   = int(slab_x*(slab_y-wall_thick)*slab_z*(phi_diff+1)/
     >             (2*vol_a))
      a_poly = int(a_no/chain_a)

      ! Similarly, calculate b_no.
      b_no   = int((slab_x*(slab_y-wall_thick)*slab_z-a_no*vol_a)/vol_b)
      b_poly = int(b_no/chain_b)

      ! Total number of polymers
      !local_no = a_no + b_no
      local_no  = a_poly*chain_a + b_poly*chain_b
      n_poly = a_poly+b_poly

      ! Find this node's offset, so that all particles have a UNIQUE id. (logical_id)
      local_count_tmp(cluster_rank+1) = local_no+n_fslab

      call MPI_REDUCE(local_count_tmp,local_count,256,MPI_INTEGER,
     >                MPI_SUM,MASTER,MPI_COMM_WORLD,mpi_err)
      call MPI_BCAST(local_count,256,MPI_INTEGER,MASTER,MPI_COMM_WORLD,
     >               mpi_err)
      id_offset = 0
      do i=MASTER,cluster_rank-1
         id_offset = id_offset + local_count(i+1)
      enddo

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
          x(k) = ran3()*slab_x/2.d0
          y(k) = ran3()*(slab_y-wall_thick)+wall_thick
          z(k) = ran3()*slab_z + cluster_rank*slab_z

          ! Species
          species(k)    = 1

          ! Logical ID -- number of particle within entire system
          logical_id(k) = int(id_offset + k)

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
2219          continue
              phi   = ran3() * two_pi
              theta = ran3() * pi

              dx = d_els*dcos(phi)*dsin(theta)
              dy = d_els*dsin(phi)*dsin(theta)
              dz = d_els*dcos(theta)

              ! Place next monomer some short distance away 
              ! from previous monomer.
              x(k) = x(k-1) + dx
              y(k) = y(k-1) + dy
              if (y(k) .lt. wall_thick) goto 2219
              z(k) = z(k-1) + dz

              ! Species
              species(k) = 1
            
              ! Logical ID -- number of particle within entire system
              logical_id(k) = id_offset + k

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
          x(k) = ran3()*slab_x/2.d0 + slab_x/2.d0
          y(k) = ran3()*(slab_y-wall_thick)+wall_thick
          z(k) = ran3()*slab_z + cluster_rank*slab_z

          ! Species
          species(k) = -1

          ! Logical ID -- number of particle within entire system
          logical_id(k) = id_offset + k

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
2220          continue
              phi   = ran3() * two_pi
              theta = ran3() * pi

              dx = d_els*dcos(phi)*dsin(theta)
              dy = d_els*dsin(phi)*dsin(theta)
              dz = d_els*dcos(theta)

              ! Place next monomer some short distance away 
              ! from previous monomer.
              x(k) = x(k-1) + dx
              y(k) = y(k-1) + dy
              if (y(k) .lt. wall_thick) goto 2220
              z(k) = z(k-1) + dz

              ! Species
              species(k) = -1

              ! Logical ID -- number of particle within entire system
              logical_id(k) = id_offset + k

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


      ! Add slab particles into the arrays, setting species to 0 and giving them
      ! a negative polymer id.
      do i=1, n_fslab
         k = k + 1 
         x(k) = fslab_x(i)
         y(k) = fslab_y(i)
         z(k) = fslab_z(i)
         vx(k)= 0.d0
         vy(k)= 0.d0
         vz(k)= 0.d0
         poly_x(k) = fslab_x(i)
         poly_y(k) = fslab_y(i)
         poly_z(k) = fslab_z(i)
         species(k) = 0
         logical_id(k) = id_offset+k
         members(logical_id(k)) = k
         p_flag(k) = 0
      enddo

      ! Set particle count.
      local_no = local_no + n_fslab   ! Include slab particles.
      call MPI_REDUCE(local_no,particle_count,1,MPI_INTEGER,
     >                MPI_SUM,MASTER,MPI_COMM_WORLD,mpi_err)
      call MPI_BCAST(particle_count,1,MPI_INTEGER,MASTER,
     >               MPI_COMM_WORLD,mpi_err)


      ! Scale velocities to eliminate a net-drift.
      vx_total = vx_total / particle_count
      vy_total = vy_total / particle_count
      vz_total = vz_total / particle_count
      do i=1, local_no-n_fslab
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

      local_no_orig  = local_no
      ext_send_count = 0
      ext_recv_count = 0

      do i=1, cluster_size*local_no
         members(i) = 0
      enddo

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
