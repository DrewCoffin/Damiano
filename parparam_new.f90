       MODULE parparam
#ifdef SINGLE
       use indeces
       integer :: iprecip_t
       integer, dimension(:), allocatable :: ltheta,iprecip
       real :: ke,vth,vperp
       real, dimension(:), allocatable :: mu_m,pa,mu_m2
       real, dimension(:), allocatable :: vpar0,vpar1,vpar2
       real, dimension(:), allocatable :: pa_check,penergy0
       real, dimension(:), allocatable :: Bo_p
       real :: B1_p,B2_p
       real, dimension(:), allocatable :: dB1,dB2
       double precision, dimension(:), allocatable :: r_p,ro_p
       double precision, dimension(:), allocatable :: theta_p
       double precision :: dtheta_p
       double precision, dimension(:), allocatable :: nu_p,h2p
       double precision, dimension(:), allocatable :: mu0_p,mu1_p,mu2_p
       double precision, dimension(:), allocatable :: h1p_1,h1p_2,e1p_1,e1p_2
       double precision :: pa_min,pa_max,mu_min,mu_max
       double precision, dimension(:), allocatable :: den_grid,j1p_grid, &
                                       sm_grid,den_temp,j1p_temp, &
                                       ratio_grid,den_grid0, &
                                       mu_grid,mu_grid2
       double precision, dimension(:), allocatable :: pint_temp,x1_temp,F
#else
       use indeces
       integer :: iprecip_t
       integer, dimension(:), allocatable :: ltheta,iprecip
       double precision :: ke,vth,vperp
       double precision, dimension(:), allocatable :: mu_m,pa,mu_m2
       double precision, dimension(:), allocatable :: vpar0,vpar1,vpar2
       double precision, dimension(:), allocatable :: pa_check,penergy0
       double precision, dimension(:), allocatable :: Bo_p
       double precision :: B1_p,B2_p
       double precision, dimension(:), allocatable :: dB1,dB2
       double precision, dimension(:), allocatable :: r_p,ro_p
       double precision, dimension(:), allocatable :: theta_p
       double precision :: dtheta_p
       double precision, dimension(:), allocatable :: nu_p,h2p
       double precision, dimension(:), allocatable :: mu0_p,mu1_p,mu2_p
       double precision, dimension(:), allocatable :: h1p_1,h1p_2,e1p_1,e1p_2
       double precision :: pa_min,pa_max,mu_min,mu_max
       double precision, dimension(:), allocatable :: den_grid,j1p_grid, &
                                       sm_grid,den_temp,j1p_temp, &
                                       ratio_grid,den_grid0, &
                                       mu_grid,mu_grid2
       double precision, dimension(:), allocatable :: pint_temp,x1_temp,F
#endif

       END MODULE parparam
