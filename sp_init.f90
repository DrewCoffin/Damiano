       SUBROUTINE sp_init()
       use indeces
       use param
       use mesh
       use fields
       use parparam
       implicit none
       integer :: ii
       double precision :: temp,randf
       double precision :: v_t,theta_t,rho_t,ke_t(5)
       double precision :: Bfield,h_mu,vbar,mu,h1_temp,Bo_temp
       external randf,Bfield,h_mu,mu
       vbar=0
       np2=5
       np=np2
       ke_t(1)=1*1.609d-19
       ke_t(2)=5*1.609d-19
       ke_t(3)=10*1.609d-19
       ke_t(4)=100*1.609d-19
       ke_t(5)=1000*1.609d-19
       pa(1)=pi/6.
       pa(2)=pi/6.
       pa(3)=pi/6.
       pa(4)=pi/6.
       pa(5)=pi/6.
       theta_p(1)=pi/2.
       theta_p(2)=pi/2.
       theta_p(3)=pi/2.
       theta_p(4)=pi/2.
       theta_p(5)=pi/2.
       do ii=1,5
           v_t=sqrt(2.*ke_t(ii)/me)
           vperp=v_t*sin(pa(ii))
           vpar0(ii)=v_t*cos(pa(ii))
           vpar1(ii)=vpar0(ii)
           ro_p(ii)=10.0d0*Re
           r_p(ii)=ro_p(ii)*sin(theta_p(ii))*sin(theta_p(ii))
           mu0_p(ii)=mu(r_p(ii),theta_p(ii))
           nu_p(ii)=sin(theta_p(ii))*sin(theta_p(ii))/r_p(ii)
           dtheta_p=pi/2/100
           mu1_p(ii)=mu0_p(ii)
 
           h1_temp=h_mu(r_p(ii),theta_p(ii))
           Bo_temp=Bfield(M,h1_temp)
           mu_m(ii)=0.5*me*vperp**2/Bo_temp
      enddo
 55   format(5e13.5)
      return
      END SUBROUTINE sp_init
       

