c////////////////////////////////////////////////////////////
c thin_film_calc.f - Routines used for various calculations
c
c Mike Hore, Memphis, TN, Summer 2006
c////////////////////////////////////////////////////////////

c
c tfilm_find_clusters() - This is done at the master level
c                   whenever the checkpoint is written...
c                   this saves time in communication.
c
c (Mohamed Laradji)
      subroutine tfilm_find_clusters(min_t, max_t)
      implicit none
      include 'mpif.h'
      include 'common.f'
      include 'Binary_Fluid/binary_fluid_var.f'
      include 'Thin_Film/thin_film_var.f'
      double precision r_time1, r_time2
      double precision x_trans, y_trans, z_trans
      double precision min_t, max_t
      integer ix, iy, iz, jx, jy, jz, mxy, idx_i, idx_j
      integer icell, i, j, k, ii, jj, jcell, kcell
      integer iti, itj
      integer li,li0,lj,lj0
      integer i1,i2,j1,j2
      double precision ytmp_i, ytmp_j

      if (DISABLE_CLUSTERS .eq. 1) return
      if (cluster_rank .ne. MASTER) return

      ! Use the information in DBL_BUF to find clusters.
      r_time1 = MPI_WTIME()

      mxy = mx*my

      ! Clear out old data from list, head.
      do i=1, ncell
         head(i) = 0
      enddo

      do i=1, particle_count
         list(i) = 0
      enddo

      do i=1, particle_count
          clust_AA(i)=0
          size_AA(i)=0
          head_AA(i)=0
          list_AA(i)=0
          clust_BB(i)=0
          size_BB(i)=0
          head_BB(i)=0
          list_BB(i)=0
      enddo
      nbr_clust_AA=0
      nbr_clust_BB=0
      cl_AA=0
      cl_BB=0
 
      ! Consider new particle positions:
      do i=1, particle_count
         j = (i-1)*8
         
         if (DBL_BUF1(j+1) .le. 0.d0) then
             DBL_BUF1(j+1) = DBL_BUF1(j+1) + dim_x
         else if (DBL_BUF1(j+1) .gt. dim_x) then 
             DBL_BUF1(j+1) = DBL_BUF1(j+1) - dim_x
         endif

         if (DBL_BUF1(j+2) .le. 0.d0) then
             DBL_BUF1(j+2) = DBL_BUF1(j+2) + dim_y
         else if (DBL_BUF1(j+2) .gt. dim_y) then 
             DBL_BUF1(j+2) = DBL_BUF1(j+2) - dim_y
         endif

         if (DBL_BUF1(j+3) .le. 0.d0) then
             DBL_BUF1(j+3) = DBL_BUF1(j+3) + dim_z
         else if (DBL_BUF1(j+3) .gt. dim_z) then 
             DBL_BUF1(j+3) = DBL_BUF1(j+3) - dim_z
         endif

         icell = 1 + int(DBL_BUF1(j+1)*cellx_inv)    + 
     >               int(DBL_BUF1(j+2)*celly_inv)*mx +
     >               int(DBL_BUF1(j+3)*cellz_inv)*mxy

         ! Update table.
         j           = head(icell)
         head(icell) = i
         list(i)     = j
      enddo

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
4174                 if (jcell .eq. icell)  j=list(i)
 
                     do while (j .gt. 0)
                       idx_i = (i-1)*8
                       idx_j = (j-1)*8

                       ytmp_i = DBL_BUF1(idx_i+2)
                       ytmp_j = DBL_BUF1(idx_j+2)

                       dx = DBL_BUF1(idx_i+1)-DBL_BUF1(idx_j+1)+ x_trans
                       dy = DBL_BUF1(idx_i+2)-DBL_BUF1(idx_j+2)+ y_trans
                       dz = DBL_BUF1(idx_i+3)-DBL_BUF1(idx_j+3)+ z_trans

                       dr = dx**2 + dy**2 + dz**2

                       ! Are the particles between min_t and max_t?
                       if ((ytmp_i.lt.min_t .or. ytmp_i .gt. max_t).or.
     >                   (ytmp_j .lt. min_t .or. ytmp_j .gt. max_t))then
                         goto 109
                       endif
  
                       ! Are the particles close enough to interact?
                       if (dr .lt. 1.d0) then
                         iti = int(DBL_BUF1(idx_i+4))
                         itj = int(DBL_BUF1(idx_j+4))
                         if (iti .eq. itj) then

c////////////////////////////////////////////////////////////////////////
C First, search for A- (polymer-) clusters
                       if(iti.eq.1)then
                        if(clust_AA(i).eq.0.and.clust_AA(j).eq.0)then
                          cl_AA=cl_AA+1
                          clust_AA(i)=cl_AA
                          clust_AA(j)=cl_AA
                          size_AA(cl_AA)=2
                          head_AA(cl_AA)=i
                          list_AA(i)=j
                          nbr_clust_AA=nbr_clust_AA+1
                        elseif(clust_AA(i).eq.0)then
                          lj=clust_AA(j)
                          clust_AA(i)=lj
                          size_AA(lj)=size_AA(lj)+1
                          j1=head_AA(lj)
                          head_AA(lj)=i
                          list_AA(i)=j1
                        elseif(clust_AA(j).eq.0)then
                          li=clust_AA(i)
                          clust_AA(j)=li
                          size_AA(li)=size_AA(li)+1
                          i1=head_AA(li)
                          head_AA(li)=j
                          list_AA(j)=i1
                        elseif(clust_AA(i).gt.clust_AA(j))then
                          li0=clust_AA(i)
                          lj0=clust_AA(j)
                          j1=head_AA(lj0)
                          i1=head_AA(li0)
                          head_AA(lj0)=i1
                          i2=i1
                          do while(i2.gt.0)
                            i1=i2
                            clust_AA(i1)=lj0
                            i2=list_AA(i1)
                          enddo
                          list_AA(i1)=j1
                          head_AA(li0)=0
                          size_AA(lj0)=size_AA(lj0)+size_AA(li0)
                          size_AA(li0)=0
                          nbr_clust_AA=nbr_clust_AA-1
                        elseif(clust_AA(i).lt.clust_AA(j))then
                          lj0=clust_AA(j)
                          li0=clust_AA(i)
                          i1=head_AA(li0)
                          j1=head_AA(lj0)
                          head_AA(li0)=j1
                          j2=j1
                          do while(j2.gt.0)
                            j1=j2
                            clust_AA(j1)=li0
                            j2=list_AA(j1)
                          enddo
                          list_AA(j1)=i1
                          head_AA(lj0)=0
                          size_AA(li0)=size_AA(li0)+size_AA(lj0)
                          size_AA(lj0)=0
                          nbr_clust_AA=nbr_clust_AA-1
                        endif
C Second, Search for B- (solvent-) clusters
                       elseif(iti.eq.-1)then
                        if(clust_BB(i).eq.0.and.clust_BB(j).eq.0)then
                          cl_BB=cl_BB+1
                          clust_BB(i)=cl_BB
                          clust_BB(j)=cl_BB
                          size_BB(cl_BB)=2
                          head_BB(cl_BB)=i
                          list_BB(i)=j
                          nbr_clust_BB=nbr_clust_BB+1
                        elseif(clust_BB(i).eq.0)then
                          lj=clust_BB(j)
                          clust_BB(i)=lj
                          size_BB(lj)=size_BB(lj)+1
                          j1=head_BB(lj)
                          head_BB(lj)=i
                          list_BB(i)=j1
                        elseif(clust_BB(j).eq.0)then
                          li=clust_BB(i)
                          clust_BB(j)=li
                          size_BB(li)=size_BB(li)+1
                          i1=head_BB(li)
                          head_BB(li)=j
                          list_BB(j)=i1
                        elseif(clust_BB(i).gt.clust_BB(j))then
                          li0=clust_BB(i)
                          lj0=clust_BB(j)
                          j1=head_BB(lj0)
                          i1=head_BB(li0)
                          head_BB(lj0)=i1
                          i2=i1
                          do while(i2.gt.0)
                            i1=i2
                            clust_BB(i1)=lj0
                            i2=list_BB(i1)
                          enddo
                          list_BB(i1)=j1
                          head_BB(li0)=0
                          size_BB(lj0)=size_BB(lj0)+size_BB(li0)
                          size_BB(li0)=0
                          nbr_clust_BB=nbr_clust_BB-1
                        elseif(clust_BB(i).lt.clust_BB(j))then
                          lj0=clust_BB(j)
                          li0=clust_BB(i)
                          i1=head_BB(li0)
                          j1=head_BB(lj0)
                          head_BB(li0)=j1
                          j2=j1
                          do while(j2.gt.0)
                            j1=j2
                            clust_BB(j1)=li0
                            j2=list_BB(j1)
                          enddo
                          list_BB(j1)=i1
                          head_BB(lj0)=0
                          size_BB(li0)=size_BB(li0)+size_BB(lj0)
                          size_BB(lj0)=0
                          nbr_clust_BB=nbr_clust_BB-1
                        endif
                       endif
c end clusters search
c////////////////////////////////////////////////////////////////////////
                         endif ! end if (iti .eq. itj)
                       endif ! end if (dr .lt. 1)
109                    continue
                       j = list(j)
                     enddo ! end do-while (j .gt. 0)
                     j = head(jcell)
                     i = list(i)

                     ! If there are more particles, proceed
                     if (i .gt. 0) goto 4174
                  endif ! end (j .gt. 0)                
               enddo ! end kcell loop
            endif ! end (i .gt. 0)
          enddo ! end ix loop
        enddo ! end iy loop
      enddo ! end iz loop

      open(41,file='A_clusters.data',access='APPEND')
      open(42,file='B_clusters.data',access='APPEND')
      write(41,*)itime,nbr_clust_AA
      write(42,*)itime,nbr_clust_BB
      do k=1,cl_AA
        if(size_AA(k).gt.0)then
          write(41,*)size_AA(k)
        endif
      enddo
      do k=1,cl_BB
        if(size_BB(k).gt.0)then
          write(42,*)size_BB(k)
        endif
      enddo
      close(41)
      close(42)

      r_time2 = MPI_WTIME()
      !if (cluster_rank .eq. MASTER) then
      !   write(*,*) 'find_clusters() takes: ', r_time2-r_time1
      !endif
      return
      end
