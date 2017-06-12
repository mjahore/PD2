c//////////////////////////////////////////
c stfac.f - Structure factor routines.
c      
c
c Mike Hore, Hilton Beach, ON, July 2005
c/////////////////////////////////////////

c
c bf_comp_stfac() - Structure factor of the b-polymer, for R(t).
c                (parallel, of course with a serial FFT)
c
      subroutine bf_comp_stfac()
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      double precision r_time1, r_time2

      integer nn, l, nsite, number_q
      parameter(nn=3, number_q=64, l=2*number_q, nsite=l**3)
      real*4  data(2*nsite), data1(2*nsite), sf(nsite)
      real*4  sfunc(0:2*number_q)
      integer num(0:2*number_q)
      integer ndim(3)
      double precision xsize, ysize, zsize
      double precision param, qq
      integer i,j,k,kval
      integer ix,iy,iz,ii
 
      if (DISABLE_STFAC .eq. YES) return
     
      ! Timing.
      r_time1 = MPI_WTIME()
      
      xsize = dim_x/dfloat(l)
      ysize = dim_y/dfloat(l)
      zsize = dim_z/dfloat(l)
      param = xsize*ysize*zsize/dfloat(nsite)
      
      ndim(1) = l
      ndim(2) = l
      ndim(3) = l
      
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
         ix = (x(i)/dim_x)*l + 1
         iy = (y(i)/dim_y)*l
         iz = (z(i)/dim_z)*l
         j  = ix + iy*l + iz*l**2
         k  = (j-1)*2
         
         ! Take species into account to get the structure.
         if (species(i) .eq. -1) then
            data1(k+1) = data1(k+1)+0.5
            data1(k+2) = 0.d0
         else
            data1(k+1) = data1(k+1)-0.5
            data1(k+2) = 0.d0
         endif
      enddo 
      
      ! Collect the composition information on the master node.
      k = 2*nsite
      call MPI_REDUCE(data1, data, k,MPI_REAL,MPI_SUM,
     >                MASTER, MPI_COMM_WORLD, mpi_err)
     
      ! Fast Fourier Transform (fourn.f)
      call fourn(data,ndim,nn,1)
      
      ! |Q|**2 = re(Q)**2 + im(Q)**2
      do i=1, nsite
         k = (i-1)*2
         sf(i) = (data(k+1)**2 + data(k+2)**2)*param
      enddo
      
      ! Average spherically
      do i=1, nsite
         iz = int((i-1)/l**2)
         ii = mod(i-1, l**2) + 1
         ix = mod(ii-1, l)
         iy = int((ii-1)/l)
         
         ! Fold back if we go over l/2
         if (ix .gt. l/2)ix=l-ix
         if (iy .gt. l/2)iy=l-iy
         if (iz .gt. l/2)iz=l-iz
         
         qq   = dsqrt(dfloat(ix**2 + iy**2 + iz**2))
         kval = int(qq)
         
         sfunc(kval) = sfunc(kval) + sf(i)
         num(kval)   = num(kval) + 1
      enddo
      
      if (cluster_rank .eq. MASTER) then
         ! Finish average
         open(11, file="stfac.data", access="append")
         do i=0,l/2
            if (num(i) .gt. 0) then
               sfunc(i) = sfunc(i)/num(i)
            endif
         enddo
      
         do i=1, number_q
            qq = dfloat(i)*two_pi/dim_x
            write(11,9999)qq, sfunc(i)
         enddo
         close(11)
      endif

      r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then 
      !   write(*,*) 'comp_stfac() takes: ', r_time2-r_time1
      !endif
9999  format(10E15.7)
      return
      end

