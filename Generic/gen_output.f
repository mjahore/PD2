!
! Generic Output Routines for XYZ file formats, etc.
!
!
! Mike Hore, hore@case.edu
! Case Western Reserve University
! January 2017
!

c
c gen_write_xyz() - Write all particle coordinates into XYZ format for
c                   use with VMD, VisIt, etc.
c
c NOTE: You must call collect_master_level3d() first. The information
c       will be in DBL_BUF1
c
      subroutine gen_write_xyz(timestep)
      implicit none
      integer i,k,timestep
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'

      if (DISABLE_OUTPUT .eq. 1) return

      if (cluster_rank .eq. MASTER) then 
          write(*,*) ''
          write(*,*) 'Writing XYZ trajectory...'
          open(12, file="all_trajectory.xyz", access="append")
          open(13, file="trajectory.xyz")

          write(12,*) particle_count
          write(12,*) 'Particle Coordinates for timestep ', timestep
          write(13,*) particle_count
          write(13,*) 'Particle Coordinates for timestep ', timestep

          do i=1, particle_count
             k = (i-1)*11
             if (DBL_BUF1(k+4) .eq. 1) then
                write(12,*) 'C', DBL_BUF1(k+1), DBL_BUF1(k+2),
     >                           DBL_BUF1(k+3)
                write(13,*) 'C', DBL_BUF1(k+1), DBL_BUF1(k+2),
     >                           DBL_BUF1(k+3)
             elseif (DBL_BUF1(k+4) .eq. -1) then
                write(12,*) 'N', DBL_BUF1(k+1), DBL_BUF1(k+2),
     >                           DBL_BUF1(k+3)
                write(13,*) 'N', DBL_BUF1(k+1), DBL_BUF1(k+2),
     >                           DBL_BUF1(k+3)
             elseif (DBL_BUF1(k+4) .eq. 0) then
                write(12,*) 'H', DBL_BUF1(k+1), DBL_BUF1(k+2),
     >                           DBL_BUF1(k+3)
                write(13,*) 'H', DBL_BUF1(k+1), DBL_BUF1(k+2),
     >                           DBL_BUF1(k+3)
             elseif (DBL_BUF1(k+4) .eq. 2) then
                write(12,*) 'O', DBL_BUF1(k+1), DBL_BUF1(k+2),
     >                           DBL_BUF1(k+3)
                write(13,*) 'O', DBL_BUF1(k+1), DBL_BUF1(k+2),
     >                           DBL_BUF1(k+3)
             elseif (DBL_BUF1(k+4) .eq. 3) then
                write(12,*) 'S', DBL_BUF1(k+1), DBL_BUF1(k+2),
     >                           DBL_BUF1(k+3)
                write(13,*) 'S', DBL_BUF1(k+1), DBL_BUF1(k+2),
     >                           DBL_BUF1(k+3)
             endif
          enddo
          close(12)
          close(13)
      endif

      return
      end

