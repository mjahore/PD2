c
c These variables hold the coordinates of the film-slab.
c That's about all that's needed.
      integer n_fslab
      double precision film_min, film_max
      common/film_reg/film_min,film_max
      double precision lattice_const,chi_wall_a,chi_wall_b,wall_thick
      double precision fslab_x(SUPP_SLAB_SIZE), fslab_y(SUPP_SLAB_SIZE),
     >                 fslab_z(SUPP_SLAB_SIZE)
      common/fslabs/fslab_x,fslab_y,fslab_z,chi_wall_a,chi_wall_b,
     >              lattice_const,wall_thick,n_fslab 
