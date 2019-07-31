          SUBROUTINE v_initialize(imode,PTEST1,PTEST3,DELTAF_DISTR,ro_r)

          use indeces
          use param
          use fields
          use mesh
          use parparam
          use e1terms
          use distgen
          implicit none
          integer :: i,imode,localseed
          integer :: localseed_vx,localseed_vy
          double precision :: temp_exp,h1_temp,Bo_temp,h_mu,Bfield,ro_r
          double precision :: vfil,B_eq,Btemp,const
          double precision :: localconst
          logical PTEST1,PTEST3,DELTAF_DISTR 

          external h_mu,Bfield,vfil
          localconst=Re*Re*Re
          seed=1000001
          do i=1,np2
             if(PTEST1)pa(i)=85.d0*pi/180.d0
             if(PTEST3)pa(i)=85.d0*pi/180.d0
             if(PTEST1)temp_exp=1.0
!            vx(i)=vx(i)*vth*vfil(ro_p(i))
!            vy(i)=vy(i)*vth*vfil(ro_p(i))
!            vz(i)=vz(i)*vth*vfil(ro_p(i))
             vx(i)=vx(i)*vth
             vy(i)=vy(i)*vth
             vz(i)=vz(i)*vth
             vperp=sqrt(vx(i)**2+vy(i)**2)/vn
 10          vpar0(i)=vz(i)/vn
             vpar1(i)=vpar0(i)
             pa(i)=atan(vperp/vpar0(i))
             B_eq=Bfield(M,ro_p(i)**3*localconst)/Bn
             Btemp=Bfield(M,h1p_1(i)*localconst)/Bn
             pa_check(i)=abs(asin(sqrt(B_eq/Btemp)*sin(pa(i))))
             if((pa_check(i)).lt.(3.0d0*pi/180.d0).and.(NO_TC))then
                 call distgen_extra(pi,localseed_vx,vx(i))
                 call distgen_extra(pi,localseed_vy,vy(i))
                 call distgen_extra(pi,localseed,vz(i))
                 vx(i)=vx(i)*vth
                 vy(i)=vy(i)*vth
                 vperp=sqrt(vx(i)**2+vy(i)**2)/vn
                 vz(i)=vz(i)*vth
                 goto 10
             endif
             if(imode.eq.2)theta_p(i)=pi/2.
             if(imode.eq.2)r_p(i)=10.
             h1_temp=h_mu(r_p(i)*L,theta_p(i))
             Bo_temp=Bfield(M,h1_temp)
             mu_m(i)=0.5*me*(vperp*vn)**2/Bo_temp/mu_mn
             mu_m2(i)=0.5*me*(vpar1(i)*vn)**2/Bo_temp/mu_mn
             IF(DELTAF_DISTR)THEN
                vpar0(i)=0d0
                vpar1(i)=0.d0
                mu_m(i)=0.d0
             ENDIF
!             penergy0(i)=sqrt(0.5*Bo_temp*mu_m(i)*mu_mn/me)/vn
          penergy0(i)=vperp*vperp+vpar0(i)*vpar0(i)
          enddo
          return
          END SUBROUTINE v_initialize

          FUNCTION vfil(ro_p)
          use param
          implicit none
          double precision :: ro_p,vfil,x2max,x2min
          double precision :: x2mid,sigma,x2temp
          x2max=1.d0/ro2
          x2min=1.d0/ro1
          x2mid=0.5d0*(x2min+x2max)
          sigma=0.85*(x2max-x2mid)
          x2temp=1.d0/ro_p
          vfil=exp(-0.5d0*((x2temp-x2mid)/sigma)**20)
          if(x2temp.gt.x2mid)vfil=1.d0
          return
          END FUNCTION vfil
          
