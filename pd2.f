c////////////////////////////////////////////////////////
c ParallelSim
c
c    DPD Simulation of large binary
c    fluids.
c
c Requires: MPI, Fortran77
c Integrator: Velocity-Verlet
c
c Mike Hore, The University of Memphis
c            Spring 2004 - Summer 2005
c
c
c////////////////////////////////////////////////////////
      program PD2
      implicit none
      integer i, j, k, h, g
   
c MPI functions/subroutines/defines
      include 'mpif.h'

c Variables/common blocks
      include 'common.f'

c Bring up MPI_COMM_WORLD communicator
      call MPI_INIT(mpi_err)
                               
c Get node position in cluster
      call MPI_COMM_RANK(MPI_COMM_WORLD, cluster_rank, mpi_err)
                                                   
c Get cluster size
      call MPI_COMM_SIZE(MPI_COMM_WORLD, cluster_size, mpi_err)
 
c A parameter we need
      LAST_SLAVE = (cluster_size) - 1
     
c Set values of nix, niy, niz arrays for the system cells (3D, pair-wise)
      nix(1)=0
      nix(2)=-1
      nix(3)=-1
      nix(4)=-1
      nix(5)=0
      nix(6)=0
      nix(7)=-1
      nix(8)=1
      nix(9)=-1
      nix(10)=0
      nix(11)=1
      nix(12)=-1
      nix(13)=0
      nix(14)=1

      niy(1)=0
      niy(2)=0
      niy(3)=-1
      niy(4)=1
      niy(5)=1
      niy(6)=0
      niy(7)=0
      niy(8)=0
      niy(9)=-1
      niy(10)=-1
      niy(11)=-1
      niy(12)=1
      niy(13)=1
      niy(14)=1

      niz(1)=0
      niz(2)=0
      niz(3)=0
      niz(4)=0
      niz(5)=0
      niz(6)=1
      niz(7)=1
      niz(8)=1
      niz(9)=1
      niz(10)=1
      niz(11)=1
      niz(12)=1
      niz(13)=1
      niz(14)=1

c Some node-specific values:
      CPU_up   = cluster_rank + 1
      CPU_down = cluster_rank - 1

c Periodic boundary conditions on CPUs.
      if (cluster_rank .eq. MASTER)           CPU_down = LAST_SLAVE
      if (cluster_rank .eq. LAST_SLAVE)       CPU_up   = MASTER

      if (cluster_rank .eq. MASTER) then
         write(*,*) 'PD2 by Michael J. A. Hore, Case Western Reserve
     > University'
         write(*,*) 'Copyright (c) 2004-2017 Michael J. A. Hore'
         write(*,*) ''
      endif

c Bootstrap the simulation depending on the type specified in
c common.f
      if (SIM_TYPE .eq. BINARY_FLUID) then
         call bf_read_dist_parameters()
         call bf_entry()
      else if (SIM_TYPE .eq. THIN_FILM) then
         call tfilm_read_dist_parameters()
         call tfilm_create_fslab()
         call tfilm_entry()
      endif

      stop
      end

c/////////////////////////////////////////////////////////////////////////////
c
c Subroutines for a binary polymer fluid. (DPD)
      include 'Binary_Fluid/binary_fluid.f'             ! VV-Integration, etc.
      include 'Binary_Fluid/binary_fluid_calc.f'        ! Calculations
      include 'Binary_Fluid/binary_fluid_checks.f'      ! Checks - temp, etc.
      include 'Binary_Fluid/binary_fluid_checkpoint.f'  ! Checkpointing
      include 'Binary_Fluid/binary_fluid_io.f'          ! File I/O
      include 'Binary_Fluid/binary_fluid_profiles.f'    ! Density profiles
      include 'Binary_Fluid/binary_fluid_stfac.f'       ! Structure Factor
c
c/////////////////////////////////////////////////////////////////////////////

c/////////////////////////////////////////////////////////////////////////////
c
c Subroutines for thin film support. (DPD)
      include 'Thin_Film/thin_film.f'                   ! Thin film, main code
      include 'Thin_Film/thin_film_calc.f'              ! Clusters, etc.
      include 'Thin_Film/thin_film_profiles.f'          ! Density profiles
      include 'Thin_Film/thin_film_stfac.f'             ! Stfac in XZ plane (2d)
c
c/////////////////////////////////////////////////////////////////////////////

c/////////////////////////////////////////////////////////////////////////////
c GENERAL MPI/PARALLELSIM STUFF:
c Subroutines for tracking particles within the nodes and for
c nearest-neighbor determination.
      include 'box2d.f'
      include 'box3d.f'

c Serial Fast-Fourier Transform (3D)
      include 'fourn.f'

c Subroutines for reading in data from disk.
      include 'io.f'

c Subroutine containing Quicksort algorithm.
      include 'qsort.f'

c ran3() Pseudo-random number generator, from Numerical Recipes.
      include 'ran3.f'

c Subroutines for output in standardized formats (e.g., XYZ)
      include 'Generic/gen_output.f'
      include 'Generic/dynamics.f'
c/////////////////////////////////////////////////////////////////////////////
