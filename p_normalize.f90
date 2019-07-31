        SUBROUTINE p_normalize()
        use indeces
        use param
        use fields
        use mesh
        use parparam
        implicit none
        integer ii

        do ii=1,np
           mu0_p(ii)=mu0_p(ii)*L*L
           mu1_p(ii)=mu1_p(ii)*L*L
           nu_p(ii)=nu_p(ii)*L
           ro_p(ii)=ro_p(ii)/L
           r_p(ii)=r_p(ii)/L
           vpar0(ii)=vpar0(ii)/vn
           vpar1(ii)=vpar1(ii)/vn
           write(*,*) mu_m(ii)
           mu_m(ii)=mu_m(ii)/(me*vn*vn/Bn)
           write(*,*) mu_m(ii)
       enddo
       return
       END SUBROUTINE p_normalize
