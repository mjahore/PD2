c
c These are nanoparticle specific variables.
c
      ! Max number of nanoparticles
      integer MAX_NP
      parameter(MAX_NP=8192)

      ! Number of particles, subparticles, membership.
      integer n_count, n_graft
      integer n_local, chain_g
      common/np_int/n_count,n_graft,n_local,chain_g

      ! Core overlap
      double precision CHI_CORE
      parameter(CHI_CORE=500.d0)

      ! Radius, interaction parameter, etc.
      double precision n_radius, n_chi,n_mass, g_radius,g_mass,g_comp
      double precision np_bond
      common/np_dbl/n_radius, n_chi,n_mass,g_radius,g_mass,g_comp,
     >              np_bond
