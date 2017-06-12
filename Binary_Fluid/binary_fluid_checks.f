c///////////////////////////////////////////////////////////
c binary_fluid_checks.f - Temperature/momentum calcs, etc.
c
c
c Mike Hore, University of Memphis, Summer 2004
c///////////////////////////////////////////////////////////


c
c bf_momentum_check() - Check momentum of entire system.
c
      subroutine bf_momentum_check()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      double precision vx_t, vy_t, vz_t ! Velocity totals
      double precision vx1t, vy1t, vz1t ! Temporary
      integer i

      vx_t = 0.d0
      vy_t = 0.d0
      vz_t = 0.d0

      ! Calculate total velocities.
      do i=1, local_no
        vx_t = vx_t + vx(i)
        vy_t = vy_t + vy(i)
        vz_t = vz_t + vz(i)
      enddo


      ! Let the master node reduce all values.
      call MPI_REDUCE(vx_t, vx1t, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
      call MPI_REDUCE(vy_t, vy1t, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
      call MPI_REDUCE(vz_t, vz1t, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
      
      ! If master node, print:
      if (cluster_rank .eq. MASTER) then
        write(*,*) 'Momentum: ', vx1t, 
     >                           vy1t,
     >                           vz1t
      endif
      
      return
      end

c
c bf_temperature_check() - Calculate temperature of entire system.
c
      subroutine bf_temperature_check()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      double precision v_t ! Velocity total (squared)
      double precision v1t ! Temporary
      integer i, degrees_freedom

      v_t = 0.d0

      degrees_freedom = 3 * particle_count -3 

      ! Calculate total velocities.
      do i=1, local_no
        v_t = v_t + vx(i)**2 + vy(i)**2 + vz(i)**2
      enddo


      ! Let the master node reduce all values.
      call MPI_REDUCE(v_t, v1t, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
      
      ! If master node, print:
      if (cluster_rank .eq. MASTER) then
        write(*,*) 'Temperature: ', v1t/degrees_freedom
      endif
      
      return
      end


c
c bf_particle_check() - Gets number of DPD particles from each node and 
c                  displays it on the master node as a total of all
c                  particles simulated.
c
      subroutine bf_particle_check()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'

      integer total_local_no

      call MPI_REDUCE(local_no, total_local_no, 1, MPI_INTEGER, MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
 
      if (cluster_rank .eq. MASTER) then 
         write(*,*) 'Simulating ', total_local_no, ' particles.'
      endif
      return
      end

c
c bf_lmn() - Function to return the number of DPD particles being
c            simulated.
c
      integer function bf_lmn()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'

      call MPI_REDUCE(local_no, bf_lmn, 1, MPI_INTEGER, MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
      return
      end
 
c
c bf_sys_temp() - Returns the temperature of the TOTAL system. Note that this
c              function is valid only on the MASTER node.
c
      double precision function bf_sys_temp()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'

      double precision v_t ! Velocity total (squared)
      double precision v1t ! Temporary
      double precision degrees_freedom
      integer i, lmn

      v_t = 0.d0
      v1t = 0.d0

      degrees_freedom = 3 * (particle_count) - 3 
      !degrees_freedom = 3 * (particle_count)*cluster_size -3 

      ! Calculate total velocities.
      do i=1, local_no
        v_t = v_t + vx(i)**2 + vy(i)**2 + vz(i)**2
      enddo


      ! Let the master node reduce all values.
      call MPI_REDUCE(v_t, v1t, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)

      bf_sys_temp = v1t/degrees_freedom
      return
      end

c
c bf_sys_momentum() - Returns net momentum of system.
c
      subroutine bf_sys_momentum(px, py, pz)
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'

      double precision vx_t, vy_t, vz_t ! Velocity totals
      double precision vx1t, vy1t, vz1t ! Temporary
      double precision px, py, pz       ! Will be returned.
      double precision mass
      integer i

      vx_t = 0.d0
      vy_t = 0.d0
      vz_t = 0.d0

      ! Calculate total velocities.
      do i=1, local_no
        mass = 1.d0
        vx_t = vx_t + mass*vx(i)
        vy_t = vy_t + mass*vy(i)
        vz_t = vz_t + mass*vz(i)
      enddo


      ! Let the master node reduce all values.
      call MPI_REDUCE(vx_t, vx1t, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
      call MPI_REDUCE(vy_t, vy1t, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
      call MPI_REDUCE(vz_t, vz1t, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
      
      px = vx1t
      py = vy1t
      pz = vz1t

      return
      end
