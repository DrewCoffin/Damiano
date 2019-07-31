        MODULE param
        use indeces
        integer :: i_inertia,iparticle,nperiod,iprofile,seed,nbin,nt_p
        integer :: ii_vpamax,ipa_profile,icouple,n1p,n2p,ifilter
        integer :: np,index,nt,index_df,nt_df,index_fw,nt_fw,icoor
        integer :: index_p,ipwrite,n1p_b,n1p2,np2,ipress,imu,igyro
        integer :: iparticle_write
        logical :: ETEST,ETEST2,REFLECT,CENTRED,NEWBC,PRECIP,RESTART
        logical :: NO_TC
        double precision :: dn1p,dn1p_b
        double precision :: ro1,ro2,r_ionosphere,pi,t,dt,mphi,Re,M
        double precision :: mu_o,L,vn,Bn,rho_n,period,omega
        double precision :: tn,tmax,twrite,dtwrite,twrite_fw,dtwrite_fw
        double precision :: me,mp,c,cn,e,dro,dro_p,ro_d_min,ro_d_max
        double precision :: theta1,theta_n1,dtheta,rho_eq
        double precision :: twrite_p,dtwrite_p
        double precision :: constE_mu,const_ec,jn,en,kT
        double precision :: version,vmax,vmin,lambda_b,max_area,mu_mn
        double precision :: x1p1,x1pn1,kz,x2p1,x2pn2
        double precision :: const_jp,const_sm,const_mu,T_i
        double precision :: epsilon_o
        double precision :: fwhm_par,Amp,L_perp,fwhm
        double precision :: vp_min,vp_max
        END MODULE param
