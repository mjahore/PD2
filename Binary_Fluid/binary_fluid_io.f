c//////////////////////////////////////////////////////////////////////
c binary_fluid_io.f - I/O routines, e.g. writing configurtions to disk
c                     and reading simulation parameters from disk.
c
c Mike Hore, January 2005, University of Memphis, Memphis, TN
c//////////////////////////////////////////////////////////////////////

c
c bf_read_config: Read simlation parameters from 'config' in 
c              working directory.
c
      subroutine bf_read_config()
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'

      ! Open file
      open(12, file="bf_config")

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
      read(12,*) d_els
      read(12,*) max_bond
      read(12,*) dt
      read(12,*) dim_x
      read(12,*) dim_y
      read(12,*) dim_z
      read(12,*) iseed
      read(12,*) sigma
      read(12,*) comp      
      read(12,*) vol_a
      read(12,*) vol_b
      read(12,*) phi_diff
      read(12,*) sd

      ! Close file
      close(12)
      return
      end

c
c bf_write_poly_master_level(): Writes all particles to disk
c                           at the master level. This method is much 
c                           slower than the write_out() subroutine. This requires
c                           that we call collect_master_level() first.
c
      subroutine bf_write_poly_master_level(p_id)
      implicit none
      include 'common.f'
      include 'mpif.h'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer i,j,k,recv_count,p_id

      if (DISABLE_OUTPUT .eq. 1) return

      ! Return if the code tries to write solvent particles.
      if (p_id .eq. 1 .and. chain_b .eq. 1) return

      if (cluster_rank .eq. MASTER) then
         open(12, file='polymer_by_id.data', access="append")
         write(12,*) itime, p_id
         do i=1, particle_count
            k = (i-1)*9
            if (int(DBL_BUF1(k+5)) .eq. p_id) then
               write(12,*) int(DBL_BUF1(k+5)),DBL_BUF1(k+1),
     >                     DBL_BUF1(k+2),DBL_BUF1(k+3),itime
            endif
         enddo
         write(12,*) ' '
         close(12)
      endif
      return
      end
c
c bf_write_out_master_level(): Writes all particles to disk
c                           at the master level. This method requires
c                           that we call collect_master_level() first.
c
      subroutine bf_write_out_master_level()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer i,j,k,recv_count

      if (DISABLE_OUTPUT .eq. 1) return

      if (cluster_rank .eq. MASTER) then
         open(12, file='species1.data')
         open(13, file='species-1.data')
         open(14, file="species0.data")
         do i=1, particle_count
            k = (i-1)*11
            if (int(DBL_BUF1(k+4)) .eq. 1) then
               write(12,*) DBL_BUF1(k+1),DBL_BUF1(k+2),DBL_BUF1(k+3),
     >                     DBL_BUF1(k+11)
            elseif (int(DBL_BUF1(k+4)) .eq. -1) then
               write(13,*) DBL_BUF1(k+1),DBL_BUF1(k+2),DBL_BUF1(k+3),
     >                     DBL_BUF1(k+11)
            elseif (int(DBL_BUF1(k+4)) .eq. 0) then
               write(14,*) DBL_BUF1(k+1),DBL_BUF1(k+2),DBL_BUF1(k+3),
     >                     DBL_BUF1(k+11)
            endif
         enddo
         close(12)
         close(13)
         close(14)
      endif
      return
      end

c
c bf_write_out_master_level_poly(): Writes all particles to disk
c                           at the master level. This method requires
c                           that we call collect_master_level() first.
c
      subroutine bf_write_out_master_level_poly()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer i,j,k,recv_count

      if (DISABLE_OUTPUT .eq. 1) return

      if (cluster_rank .eq. MASTER) then
         open(12, file='species1.data')
         open(13, file='species-1.data')
         open(14, file="species0.data")
         do i=1, particle_count
            k = (i-1)*6
            if (int(DBL_BUF1(k+4)) .eq. 1) then
               write(12,*) DBL_BUF1(k+1),DBL_BUF1(k+2),DBL_BUF1(k+3),
     >                     DBL_BUF1(k+5)
            elseif (int(DBL_BUF1(k+4)) .eq. -1) then
               write(13,*) DBL_BUF1(k+1),DBL_BUF1(k+2),DBL_BUF1(k+3),
     >                     DBL_BUF1(k+5)
            else
               write(14,*) DBL_BUF1(k+1),DBL_BUF1(k+2),DBL_BUF1(k+3),
     >                     DBL_BUF1(k+5)
            endif
         enddo
         close(12)
         close(13)
         close(14)
      endif
      return
      end

c
c bf_write_snapshot_master_level(): Writes all particles to disk
c                           at the master level. This method requires
c                           that we call collect_master_level() first.
c
      subroutine bf_write_snapshot_master_level()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      integer i,j,k,recv_count

      if (DISABLE_OUTPUT .eq. 1) return

      if (cluster_rank .eq. MASTER) then
         open(12, file='spec1_snapshot.data', access="append")
         open(13, file='spec-1_snapshot.data', access="append")
         open(14, file="spec1_velo_snapshot.data", access="append")
         open(15, file="spec-1_velo_snapshot.data", access="append")
         do i=1, particle_count
            k = (i-1)*11
            if (int(DBL_BUF1(k+4)) .eq. 1) then
               write(12,*) DBL_BUF1(k+1), DBL_BUF1(k+2), DBL_BUF1(k+3)
               write(14,*)int(DBL_BUF1(k+5)),DBL_BUF1(k+6),DBL_BUF1(k+7)
            elseif (int(DBL_BUF1(k+4)) .eq. -1) then
               write(13,*) DBL_BUF1(k+1),DBL_BUF1(k+2),DBL_BUF1(k+3)
               write(15,*)int(DBL_BUF1(k+5)),DBL_BUF1(k+6),DBL_BUF1(k+7)
            endif
         enddo
         close(12)
         close(13)
         close(14)
         close(15)
      endif
      return
      end
