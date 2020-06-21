c
c Bridge detection variables
c
c The following parameters may be changed:
c
c Bridge length --  the number of monomers to look at for
c                   potential bridge candidates.
c Bridge epsilon -- the number of degrees that the average
c                   orientation must be within to be considered 
c                   a bridge.
c
      integer          BRIDGE_LENGTH_A,BRIDGE_LENGTH_B
      integer          TOTAL_A_BRIDGE, TOTAL_B_BRIDGE
      integer          BRIDGE_LIST(SUPP_SLAB_SIZE)
      double precision BRIDGE_EPSILON
      parameter (BRIDGE_LENGTH_A=5, BRIDGE_LENGTH_B=5,
     >           BRIDGE_EPSILON=10.d0)
      common/bdglst/BRIDGE_LIST,TOTAL_A_BRIDGE,TOTAL_B_BRIDGE

c
c General variables specific to a binary fluid simulation.    
c
      integer chain_a, chain_b, a_no, b_no, b_poly, a_poly,
     >          n_poly

      double precision d_els, comp, max_bond, vol_a, vol_b, phi_diff

      common/bf_int/chain_a,chain_b,a_no,b_no,a_poly,b_poly,
     >                n_poly

      common/bf_dbl/d_els,comp,max_bond,vol_a,vol_b,phi_diff

c
c Polymer Bond Values
c 
      double precision poly_bond_a, poly_bond_b, poly_bond_g
      common/bond_lengths/poly_bond_a,poly_bond_b,poly_bond_g

c
c Clusters
c
      integer cl_AA,cl_BB
      integer nbr_clust_AA, nbr_clust_BB
      integer clust_BB(SUPP_SYS_SIZE),size_BB(SUPP_SYS_SIZE)
     >        ,head_BB(SUPP_SYS_SIZE),list_BB(SUPP_SYS_SIZE)
      integer clust_AA(SUPP_SYS_SIZE),size_AA(SUPP_SYS_SIZE)
     >        ,head_AA(SUPP_SYS_SIZE),list_AA(SUPP_SYS_SIZE)
      common/clust/size_AA,size_BB,head_AA,head_BB,list_AA,list_BB,
     >           clust_AA,clust_BB,cl_AA,cl_BB,nbr_clust_AA,nbr_clust_BB

c
c Stresses
c
      double precision sigma_xx, sigma_xy, sigma_xz
      double precision sigma_yy, sigma_yz
      double precision sigma_zz
      common/stress/sigma_xx,sigma_xy,sigma_xz,sigma_yy,sigma_yz,
     >              sigma_zz    

