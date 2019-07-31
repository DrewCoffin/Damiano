       MODULE parparam
#ifdef SINGLE
       use indeces
       integer :: ltheta(npt2),iprecip(npt2),iprecip_t
       real :: T_e,vth,mu_m(npt2),pa(npt2),vperp,mu_m2(npt2)
       real :: vpar0(npt2),vpar1(npt2),vpar2(npt2)
       real :: pa_check(npt2),penergy0(npt2)
       real :: Bo_p(npt2)
       real :: B1_p,B2_p,dB1(npt2),dB2(npt2)
       double precision :: r_p(npt2),ro_p(npt2)
       double precision :: theta_p(npt2)
       double precision :: dtheta_p,nu_p(npt2),h2p(npt2)
       double precision :: mu0_p(npt2),mu1_p(npt2),mu2_p(npt2)
       double precision :: h1p_1(npt2),h1p_2(npt2),e1p_2(npt2)
       double precision :: e1p_1(npt2)
       double precision :: pa_min,pa_max,mu_min,mu_max
       double precision :: den_grid(n21),j1p_grid(n21),sm_grid(n21)
       double precision :: den_temp(n21),j1p_temp(n21)
       double precision :: ratio_grid(n21),den_grid0(n21)
       double precision :: mu_grid(n21),mu_grid2(n21)
!      double precision :: pint_temp(n_int),x1_temp(n_int),F(n_int)
#else
       use indeces
       integer :: ltheta(npt2),iprecip(npt2),iprecip_t
       double precision :: T_e,vth,mu_m(npt2),pa(npt2),vperp,mu_m2(npt2)
       double precision :: vpar0(npt2),vpar1(npt2),vpar2(npt2)
       double precision :: r_p(npt2),ro_p(npt2),pa_check(npt2)
       double precision :: theta_p(npt2)
       double precision :: Bo_p(npt2),penergy0(npt2)
       double precision :: dtheta_p,nu_p(npt2),h2p(npt2)
       double precision :: mu0_p(npt2),mu1_p(npt2),mu2_p(npt2)
       double precision :: h1p_1(npt2),h1p_2(npt2),e1p_2(npt2)
       double precision :: B1_p,B2_p,dB1(npt2),dB2(npt2),e1p_1(npt2)
       double precision :: pa_min,pa_max,mu_min,mu_max
       double precision :: den_grid(n21),j1p_grid(n21),sm_grid(n21)
       double precision :: den_temp(n21),j1p_temp(n21)
       double precision :: ratio_grid(n21),den_grid0(n21)
       double precision :: mu_grid(n21),mu_grid2(n21)
!      double precision :: pint_temp(n_int),x1_temp(n_int),F(n_int)
#endif

       END MODULE parparam
