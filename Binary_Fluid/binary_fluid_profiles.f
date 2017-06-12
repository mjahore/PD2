c///////////////////////////////////////////////
c profiles.f - Density/etc. profiles.
c 
c Mike Hore, Memphis, TN October 2005
c//////////////////////////////////////////////

c
c bf_density_profiles() - Calculates a density profile in the
c                      direction NORMAL to that in which a
c                      slab is created using CONFIG_TYPE=PHASE SEPARATED
c
      subroutine bf_density_profiles()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer i, slab, idx
      double precision slab_thickness, slab_volume
      double precision DENSITY(2*DENSITY_SLABS)
      double precision r_time1, r_time2

      if (DISABLE_DENSITY_PROFILES .eq. YES) return

      r_time1 = MPI_WTIME()

      ! Initialize
      slab_thickness = dim_x / dfloat(DENSITY_SLABS)

      do i=1, 2*DENSITY_SLABS
         DBL_BUF(i) = 0.d0
         DENSITY(i) = 0.d0
      enddo

      ! Place the particles into slabs, and talley the number of particles
      ! per slab (this node only).
      do i=1, local_no
         slab = int(x(i)/slab_thickness) + 1
         idx  = (slab - 1)*2
  
         if (species(i) .eq. 1) then
            DBL_BUF(idx+1) = DBL_BUF(idx+1) + 1 
         else
            DBL_BUF(idx+2) = DBL_BUF(idx+2) + 1
         endif
      enddo

      ! Combine results at master level. 
      i = 2*DENSITY_SLABS
      call MPI_REDUCE(DBL_BUF,DENSITY,i,MPI_DOUBLE_PRECISION,MPI_SUM,
     >                MASTER,MPI_COMM_WORLD,mpi_err)

      ! Let the master node normalize and then write the results to file.
      if (cluster_rank .eq. MASTER) then
         slab_volume = dim_y*dim_z*slab_thickness

         open(12, file="a_density_profile.data", access="append")
         open(13, file="b_density_profile.data", access="append")
         do i=1,DENSITY_SLABS
            idx = (i-1)*2
            write(12,*) i*slab_thickness, DENSITY(idx+1)/slab_volume,
     >                  itime
            write(13,*) i*slab_thickness, DENSITY(idx+2)/slab_volume,
     >                  itime
         enddo
         r_time2 = MPI_WTIME()

         !write(*,*) 'density_profiles() takes: ', r_time2-r_time1
      endif
      return
      end
