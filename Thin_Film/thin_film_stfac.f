c//////////////////////////////////////////////
c stfac.f - Structure factor routines.
c      
c
c Mike Hore, University of Memphis, Summer 2006
c//////////////////////////////////////////////

c
c tfilm_slice_stfac() - Structure factor of the b-polymer, for R(t). Takes points along
c                       the center of the film instead of averaging over the entire thickness
c                       of the film.
c
      subroutine tfilm_slice_stfac(min_t, max_t)
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      include 'Thin_Film/thin_film_var.f'
      double precision r_time1, r_time2
      integer nn, l, nsite, number_q
      parameter(nn=2, number_q=256, l=2*number_q, nsite=l**2)
      real*4  data(2*nsite), data1(2*nsite), sf(nsite)
      real*4  sfunc(0:2*number_q)
      integer num(0:2*number_q)
      integer ndim(2)
      double precision xsize, ysize, zsize
      double precision param, qq
      double precision min_t, max_t
      integer i,j,k,kval
      integer ix,iy,iz,ii
 
      if (DISABLE_STFAC .eq. YES) return
     
      ! Timing.
      r_time1 = MPI_WTIME()
      
      xsize = dim_x/dfloat(l)
      ysize = (dim_y-wall_thick)/dfloat(l)
      zsize = dim_z/dfloat(l)
      param = xsize*ysize*zsize/dfloat(nsite)
      
      ndim(1) = l
      ndim(2) = l
      
      ! Initialize
      do i=1, 2*nsite
         data(i)  = 0.d0
         data1(i) = 0.d0
      enddo
      
      do i=0, l
         sfunc(i) = 0.d0
         num(i)   = 0.d0
      enddo
 
      ! Now we are going to look at the composition of the system, in parallel,
      ! and collect the results at the master level to FFT.    
      do i=1, local_no
         if (species(i) .eq. 0) goto 518
         
         ! Do a thickness check, make sure points are within the region
         ! of the film that we're interested in.
         if (y(i) .lt. min_t .or. y(i) .gt. max_t) goto 518

         ix = (x(i)/dim_x)*l + 1
         iz = (z(i)/dim_z)*l + 1
         j  = ix + (iz-1)*l
         k  = (j-1)*2
         
         ! Take species into account to get the structure.
         if (species(i) .eq. -1) then
            data1(k+1) = data1(k+1)+0.5
            data1(k+2) = 0.d0
         elseif (species(i) .eq. 1) then
            data1(k+1) = data1(k+1)-0.5
            data1(k+2) = 0.d0
         endif
518      continue
      enddo 
      
      ! Collect the composition information on the master node.
      k = 2*nsite
      call MPI_REDUCE(data1, data, k,MPI_REAL,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
     
      if (cluster_rank .eq. MASTER) then
         ! Fast Fourier Transform (fourn.f)
         call fourn(data,ndim,nn,1)
      
         ! |Q|**2 = re(Q)**2 + im(Q)**2
         do i=1, nsite
            k = (i-1)*2
            sf(i) = (data(k+1)**2 + data(k+2)**2)*param
         enddo
      
         ! Average spherically
         do i=1, nsite
            ix = mod(i-1, l)
            iz = int((i-1)/l)
         
            ! Fold back if we go over l/2
            if (ix .gt. l/2)ix=l-ix
            if (iz .gt. l/2)iz=l-iz
         
            qq   = dsqrt(dfloat(ix**2 + iz**2))
            kval = int(qq)
         
            sfunc(kval) = sfunc(kval) + sf(i)
            num(kval)   = num(kval) + 1
         enddo
      
         ! Finish average
         open(11, file="stfac.data", access="append")
         do i=0,l/2
            if (num(i) .gt. 0) then
               sfunc(i) = sfunc(i)/num(i)
            endif
         enddo
      
         do i=1, number_q
            qq = dfloat(i)*two_pi/dim_x
            write(11,*)qq, sfunc(i)
         enddo
         close(11)
      endif

      !r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then 
      !   write(*,*) 'comp_stfac() takes: ', r_time2-r_time1
      !endif
      return
      end

