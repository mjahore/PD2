c////////////////////////////////////////////////
c common.f - Variable declarations
c
c
c Mike Hore, University of Memphis, Summer 2004
c///////////////////////////////////////////////
c
c The following lines are constants. Do not change
c them.
c
      character*8 code_version
      parameter(code_version='1.2.16')
      integer YES, NO
      parameter(YES=1, NO=0)
      integer FENE, HARMONIC
      parameter(FENE=0, HARMONIC=1)
      integer LOAD_AT_RELAX, LOAD_AT_QUENCH
      parameter(LOAD_AT_QUENCH=200, LOAD_AT_RELAX=100)
      integer PHASE_SEPARATED, DISORDERED,WETTING
      parameter(PHASE_SEPARATED=0, DISORDERED=1,WETTING=2)
      integer BINARY_FLUID, MEMBRANE, MC2D, THIN_FILM, ROTAXANE
      integer BOTTLEBRUSH, NANOPARTICLE,DIBLOCK
      parameter(BINARY_FLUID=0, MEMBRANE=1, MC2D=3, THIN_FILM=4,
     >          ROTAXANE=5,BOTTLEBRUSH=6,NANOPARTICLE=7,DIBLOCK=8)
     
      integer LINEAR,RING,POLYCATENANE
      parameter(LINEAR=0, RING=1, POLYCATENANE=2)


c
c SUPPORTED SYSTEM SIZE! VERY IMPORTANT PARAMETER!!!!
c
c The following parameters may be changed:
c
      integer SUPP_SYS_SIZE
      parameter(SUPP_SYS_SIZE=4000000)
      integer SUPP_SLAB_SIZE
      parameter(SUPP_SLAB_SIZE=1000000)
 
      integer DENSITY_SLABS
      parameter(DENSITY_SLABS=50)

      integer CP_LOAD_TYPE
      parameter(CP_LOAD_TYPE=LOAD_AT_QUENCH)

      integer BOND_TYPE
      parameter(BOND_TYPE=FENE)
 
      integer DISABLE_BRIDGES
      parameter(DISABLE_BRIDGES=YES)

      integer DISABLE_DENSITY_PROFILES
      parameter(DISABLE_DENSITY_PROFILES=YES)

      integer DISABLE_OUTPUT
      parameter(DISABLE_OUTPUT=NO)

      integer DISABLE_CHECKPOINT
      parameter(DISABLE_CHECKPOINT=NO)

      integer DISABLE_CLUSTERS
      parameter(DISABLE_CLUSTERS=YES)

      integer DISABLE_STFAC
      parameter(DISABLE_STFAC=NO)

      integer CONFIG_TYPE
      parameter(CONFIG_TYPE=DISORDERED)

      integer POLYMER_TYPE
      parameter(POLYMER_TYPE=LINEAR)

      double precision LAMBDA
      parameter(LAMBDA=0.5d0)

      double precision LATTICE_BOND, LATTICE_BOND_STRENGTH
      parameter(LATTICE_BOND=0.15d0, LATTICE_BOND_STRENGTH=100.d0)

      ! Valid choices for SIM_TYPE: ROTAXANE, BINARY_FLUID, MEMBRANE, BOTTLEBRUSH
      integer SIM_TYPE
      parameter(SIM_TYPE=NANOPARTICLE)

      integer USE_MIN_MAX
      parameter(USE_MIN_MAX=NO)
c
c This is the end of parameters that may be modified.
c

c//////////////////////////////////////////////////////
c
c MPI Variables
c
      integer          INT_BUF(SUPP_SYS_SIZE)     
      double precision DBL_BUF (20*SUPP_SYS_SIZE)
      double precision DBL_BUF1(20*SUPP_SYS_SIZE)
      common/buff/DBL_BUF,DBL_BUF1,INT_BUF

      integer mpi_err, cluster_size, cluster_rank
      integer mpi_status(10),MASTER, LAST_SLAVE
      integer CPU_up, CPU_down
      common/mpi/mpi_err, cluster_size, cluster_rank, CPU_up,
     >           CPU_down, LAST_SLAVE 
      parameter(MASTER=0)

c
c MPI Transfer Buffers
c
      integer           MPI_INT_BUF1(30)
      double precision  MPI_DBL_BUF1(30)
      common/buffer/MPI_DBL_BUF1,MPI_INT_BUF1

c
c DPD specific stuff.
c
      double precision wr      ! Weight factor
      double precision a       ! Interaction strength
      double precision rdv     ! (r . v)
      double precision gamma1  ! Dissipation term
      double precision sigma   ! Friction param
c
c Space/velocity/force
c
c (Change these array sizes for more particles)
c
      double precision  x(SUPP_SLAB_SIZE), y(SUPP_SLAB_SIZE), 
     >                  z(SUPP_SLAB_SIZE)
      double precision vx(SUPP_SLAB_SIZE),vy(SUPP_SLAB_SIZE),
     >                 vz(SUPP_SLAB_SIZE)
      double precision fx(SUPP_SLAB_SIZE),fy(SUPP_SLAB_SIZE),
     >                 fz(SUPP_SLAB_SIZE)

c
c Polymer coordinates
c
      double precision poly_x(SUPP_SLAB_SIZE), poly_y(SUPP_SLAB_SIZE), 
     >                 poly_z(SUPP_SLAB_SIZE)

c These arrays are needed for the velocity Verlet integration scheme.
c Note that this information needs to also be passed between nodes when
c particles are exchanged.
      double precision vx1(SUPP_SLAB_SIZE), vy1(SUPP_SLAB_SIZE), 
     >                 vz1(SUPP_SLAB_SIZE)
      double precision fx1(SUPP_SLAB_SIZE), fy1(SUPP_SLAB_SIZE), 
     >                 fz1(SUPP_SLAB_SIZE)

      ! Distances not needed as common:
      double precision dx, dy, dz, dr

      ! Integer array tracking which particles were send to CPU_down 
      ! to satisfy neighboring slab interactions. 
      integer ext_list(SUPP_SLAB_SIZE), ext_recv_count, ext_send_count
      integer particle_count, local_no_orig


      ! Particle Information
      common/coord/x,y,z
      common/veloc/vx,vy,vz
      common/force/fx,fy,fz
      common/vv_info/vx1,vy1,vz1
      common/ff_info/fx1,fy1,fz1
      common/poly/poly_x,poly_y,poly_z
      common/extrnl/ext_list, ext_recv_count, ext_send_count

c
c Species
c
c (Change these array sizes for more particles)
c
       integer species(0:SUPP_SLAB_SIZE),members(SUPP_SYS_SIZE)
       integer logical_id(SUPP_SLAB_SIZE)
       integer p_flag(0:SUPP_SLAB_SIZE), p_species(SUPP_SYS_SIZE)
       common/spin/species, logical_id, members, p_flag, p_species


c
c System cell variables
c
        integer head(SUPP_SLAB_SIZE), list(SUPP_SYS_SIZE)
        integer nix(14), niy(14), niz(14)
        integer mx, my, mz, ncell
        double precision cellx_inv, celly_inv, cellz_inv
        common/cell/head,list,nix,niy,niz
        common/cell2/cellx_inv,celly_inv,cellz_inv,mx,my,mz,ncell
c
c System slab variables
c
         integer slabs
         double precision slab_x, slab_y, slab_z
         common/slab/slab_x,slab_y,slab_z,slabs
c
c Pseudo-random Number Variables
c
        real*4 ran3
        integer a1,a2,a3,iseed,ma(55)
        common/rand/a1,a2,a3,iseed,ma

c
c Simulation parameters
c
        integer load_cp, time, equil_time, deg_frdm,
     >          check_freq, itime, itime0, calc_freq,
     >          local_no
 
        double precision rho, temp, chi1, chi2, chi3, dt, dim_x,
     >                   dim_y, dim_z, sd

        common/param/gamma1,rho,temp,chi1,chi2,chi3,dt,sd
     >        ,dim_x,dim_y,dim_z,sigma,load_cp ,time, local_no_orig 
     >        ,equil_time,deg_frdm, check_freq, itime,itime0,
     >        particle_count, calc_freq, local_no

c
c Energy
c
      double precision utot
      common/energ/utot

c
c Constants
c
        double precision pi, two_pi
        data     pi/3.141592654d0/
        data two_pi/6.283185308d0/
