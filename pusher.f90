        SUBROUTINE push1(size2,rank,root)

        use indeces
        use param
        use mesh
        use fields
        use parparam
        use mpi
        implicit none
        integer :: ii,jj,j,lro
        integer :: k1,k2,k3,k4
        integer :: size2,rank,root,error,offset,k
        double precision :: h_mu2,length,h_nu2,l1,l2,w1,w2,dro_i
        double precision :: a1,a2,a3,a4,asum_i,c1,cost,sint,Ep_test
        double precision :: h1p_1i
        external h_mu2,length,h_nu2,Ep_test


        call MPI_BCAST(e1_1,n21,MPI_DOUBLE_PRECISION, &
                       root,MPI_COMM_WORLD,error)
                                                                                
        dro_i=1.d0/dro
        do 100 ii=1,np2
          if (iprecip(ii).eq.1)goto 100
          lro=(ro_p(ii)-ro1)*dro_i+1
          k1=lro+(ltheta(ii)-1)*n2 
          cost=cos(theta_p(ii))
          sint=sqrt(1.d0-cost*cost)
          h2p(ii)=h_nu2(r_p(ii),cost,sint)
          l2=abs(0.5*(h1(k1+n2)+h1p_1(ii))* &
                 (mu1_p(ii)-x1(ltheta(ii)+1)))
          l1=abs(0.5*(h1(k1)+h1p_1(ii))* &
                 (mu1_p(ii)-x1(ltheta(ii))))
!          l2=abs(length(ro_p(ii),pi/2-theta(lro,ltheta(ii)+1)
!     &              ,pi/2-theta1_p(ii)))
!          l1=abs(length(ro_p(ii),pi/2.-theta_p(ii),
!     &              pi/2.-theta(lro,ltheta(ii))))
          w1=abs(0.5*(h2(k1)+h2p(ii)) &
            *(x2(lro)-nu_p(ii)))
          w2=abs(0.5*(h2(k1+1)+h2p(ii)) &
            *(x2(lro+1)-nu_p(ii)))
!          c1=cosalpha(lro,ltheta(ii))
          a1=l1*w1
          a2=l1*w2
          a3=l2*w1
          a4=l2*w2
          asum_i=1.d0/(a1+a2+a3+a4)
          k1=lro+1+(ltheta(ii)+1-1)*n2
          k2=lro+(ltheta(ii)+1-1)*n2
          k3=lro+1+(ltheta(ii)-1)*n2
          k4=lro+(ltheta(ii)-1)*n2
          e1p_1(ii)=constE_mu*(a1*e1_1(k1)+ &
              a2*e1_1(k2)+ &
              a3*e1_1(k3)+ &
              a4*e1_1(k4))*asum_i
          dB1(ii)=(a1*dBo(k1)+ &
              a2*dBo(k2)+ &
              a3*dBo(k3)+ &
              a4*dBo(k4))*asum_i
          
!           call deriv(M,Bn,Re,ro_p(ii),theta_p(ii),dtheta_p,dB1(ii))
           h1p_1(ii)=h_mu2(r_p(ii),cost)
           h1p_1i=1.d0/h1p_1(ii)
!           if(ETEST2)e1p_1(ii)=constE_mu*Ep_test(theta_p(ii))
           vpar2(ii)=vpar0(ii)-2.d0*dt*(e1p_1(ii)+ &
                     mu_m(ii)*h1p_1i*dB1(ii))
           mu2_p(ii)=mu0_p(ii)+2.d0*dt*vpar1(ii)*h1p_1i
!           call rootfind(mu2_p(ii),nu_p(ii),r_p(ii),
!     &               theta_p(ii))
 100     enddo

        do ii=1,np2
          ltheta(ii)=0
          do j=1,n1-1
             if((mu2_p(ii).ge.x1(j)).and. &
                 (mu2_p(ii).le.x1(j+1)))then
                 ltheta(ii)=j
             endif
          enddo
        enddo

 55        format(5e13.5)
        return
        END SUBROUTINE push1

        SUBROUTINE push2(size2,rank,root)
        use indeces
        use param
        use mesh
        use fields
        use parparam
        use mpi
        implicit none
        integer :: ii,jj,j,lro
        integer :: k1,k2,k3,k4,ltheta_t
        integer :: size2,rank,root,error,k,offset
        double precision :: h_mu2,h_nu2,c1,cost,sint,dro_i
        double precision :: w1,w2,l1,l2,a1,a2,a3,a4,asum_i,length,Ep_test
        double precision :: h1p_1i,h1p_2i,j1p_fil(n21),x2_mid
        double precision :: j1delta_det
        double precision :: penergy
        external h_mu2,h_nu2,length,Ep_test,j1delta_det


        call MPI_BCAST(e1_2,n21,MPI_DOUBLE_PRECISION, &
                       root,MPI_COMM_WORLD,error)
        call MPI_BCAST(j1p_grid,n21,MPI_DOUBLE_PRECISION, &
                       root,MPI_COMM_WORLD,error)
                                                                                
        k=rank+1
        offset=np2*(k-1)

        dro_i=1.d0/dro
        do 100 ii=1,np2
           if (iprecip(ii).eq.1)goto 100
           lro=(ro_p(ii)-ro1)*dro_i+1
           k1=lro+(ltheta(ii)-1)*n2 
           cost=cos(theta_p(ii))
           sint=sqrt(1.d0-cost*cost)
           h2p(ii)=h_nu2(r_p(ii),cost,sint)
           l2=abs(0.5*(h1(k1+n2)+h1p_2(ii))* &
             (mu1_p(ii)-x1(ltheta(ii)+1)))
           l1=abs(0.5*(h1(k1)+h1p_2(ii))* &
             (mu1_p(ii)-x1(ltheta(ii))))
!           l2=abs(length(ro_p(ii),pi/2-theta(lro,ltheta(ii)+1)
!     &              ,pi/2-theta2_p(ii)))
!           l1=abs(length(ro_p(ii),pi/2.-theta2_p(ii),
!     &              pi/2.-theta(lro,ltheta(ii))))

           w1=abs(0.5*(h2(k1)+h2p(ii)) &
            *(x2(lro)-nu_p(ii)))
           w2=abs(0.5*(h2(k1+1)+h2p(ii)) &
            *(x2(lro+1)-nu_p(ii)))
!           c1=cosalpha(lro,ltheta(ii))
           a1=l1*w1
           a2=l1*w2
           a3=l2*w1
           a4=l2*w2
           asum_i=1.d0/(a1+a2+a3+a4)
           k1=lro+1+(ltheta(ii)+1-1)*n2
           k2=lro+(ltheta(ii)+1-1)*n2
           k3=lro+1+(ltheta(ii)-1)*n2
           k4=lro+(ltheta(ii)-1)*n2
           e1p_2(ii)=constE_mu*(a1*e1_2(k1)+ &
              a2*e1_2(k2)+ &
              a3*e1_2(k3)+ &
              a4*e1_2(k4))*asum_i
          dB2(ii)=(a1*dBo(k1)+ &
              a2*dBo(k2)+ &
              a3*dBo(k3)+ &
              a4*dBo(k4))*asum_i
          Bo_p(ii)=(a1*Bo(k1)+ &
              a2*Bo(k2)+ &
              a3*Bo(k3)+ &
              a4*Bo(k4))*asum_i
!           call deriv(M,Bn,Re,ro_p(ii),theta_p(ii),dtheta_p,dB2(ii))
!           if(ETEST2)e1p_2(ii)=constE_mu*Ep_test(theta_p(ii))
           h1p_2(ii)=h_mu2(r_p(ii),cost)
           h1p_1i=1.d0/h1p_1(ii)
           h1p_2i=1.d0/h1p_2(ii)
           vpar2(ii)=vpar1(ii)-0.5d0*dt*(e1p_1(ii)+e1p_2(ii)+ &
                 mu_m(ii)*(dB1(ii)*h1p_1i+dB2(ii)*h1p_2i))
           mu2_p(ii)=mu1_p(ii)+0.5d0*dt*(vpar1(ii)*h1p_1i+ &
                     vpar2(ii)*h1p_2i)
           call rootfind(mu2_p(ii),nu_p(ii),r_p(ii),theta_p(ii))
!        penergy=vpar2(ii)*vpar2(ii)+mu_m(ii)*Bo_p(ii)
!        penergy=vpar2(ii)*vpar2(ii)
!        penergy=sqrt(Bo_p(ii)*Bn*mu_m(ii)*mu_mn/me/0.5)/vn
!        penergy=sqrt(Bo_p(ii)*mu_m(ii)*2.d0)
         penergy=vpar2(ii)*vpar2(ii)+Bo_p(ii)*mu_m(ii)*2.d0
!           if(ii.eq.np2/2)write(*,*) penergy,penergy0(ii)
 100    enddo

        do ii=1,np2
           ltheta(ii)=0
           do j=1,n1-1
           if((mu2_p(ii).ge.x1(j)).and. &
                 (mu2_p(ii).le.x1(j+1)))then
                 ltheta(ii)=j
           endif
           enddo
        enddo


! particle boundaries

        do ii=1,np2
            IF((mu2_p(ii).lt.x1p1.or.mu2_p(ii).gt.x1pn1).and. &
            iprecip(ii).eq.0)THEN
            lro=(ro_p(ii)-ro1)*dro_i+1

            if(PRECIP)then
! forced precip of particles reaching boundary
                  iprecip(ii)=1
                  iprecip_t=iprecip_t+1
                  goto 999
            endif
! determine particle current total energy
            penergy=vpar2(ii)*vpar2(ii)+Bo_p(ii)*mu_m(ii)*2.d0

            if(mu2_p(ii).le.x1p1)then
               k=lro+(3-1)*n2
!               if(penergy.lt.1.2d0*penergy0(ii))then
!c reflection of particles
!                  mu2_p(ii)=x1p1+abs(mu2_p(ii)-x1p1)
!                  mu1_p(ii)=x1p1-abs(mu1_p(ii)-x1p1)
!                  call rootfind(mu2_p(ii),nu_p(ii),
!     &                  r_p(ii),theta_p(ii))
!                 vpar2(ii)=-vpar2(ii)
!                 vpar1(ii)=-vpar1(ii)
!               else
!c precip of particles
!             write(2,*) j1p_fil(k),j1delta1,pa_check(ii),mu2_p(ii),x1p1
                  iprecip(ii)=1
                  iprecip_t=iprecip_t+1
!               endif
            endif

            if(mu2_p(ii).ge.x1pn1)then
               k=lro+(n1-2-1)*n2
!               if(penergy.lt.1.2d0*penergy0(ii))then
!c reflection of particles
!                  mu2_p(ii)=x1pn1-abs(mu2_p(ii)-x1pn1)
!                  mu1_p(ii)=x1pn1+abs(mu1_p(ii)-x1pn1)
!                  call rootfind(mu2_p(ii),nu_p(ii),
!     &                  r_p(ii),theta_p(ii))
!                  vpar2(ii)=-vpar2(ii)
!                  vpar1(ii)=-vpar1(ii)
!               else
!c precip of particles
!              write(3,*) j1p_fil(k),j1delta1,pa_check(ii),mu2_p(ii),x1p1
                  iprecip(ii)=1
                  iprecip_t=iprecip_t+1
!               endif
            endif
            ENDIF
 999    enddo
 55     format(5e13.5)
        return
        END SUBROUTINE push2

        FUNCTION Bfield(M,h1)
          implicit none
          double precision :: Bfield,M,h1,pi,mu_o
          pi=4.d0*atan(1.d0)
          mu_o=4.d-7*pi
          Bfield=mu_o*M/(4.d0*pi*h1)
          return
        END FUNCTION Bfield


        SUBROUTINE deriv2(M,Bn,Re,ro,theta,dtheta,dB)
          implicit none
          double precision :: ro,theta,mu_p,mu_m,Re
          double precision :: theta1,theta2,r1,r2
          double precision :: mu,h1p,h1m,dtheta
          double precision :: Bp,Bm,dB,h_mu,Bfield,M,Bn
          external mu,h_mu,Bfield
          theta1=theta-dtheta
          theta2=theta+dtheta

          r1=ro*sin(theta1)**2
          r2=ro*sin(theta2)**2

          mu_p=mu(r1,theta1)
          mu_m=mu(r2,theta2)

          h1p=h_mu(r1,theta1)
          h1m=h_mu(r2,theta2)

          Bp=Bfield(M,h1p*Re**3)/Bn
          Bm=Bfield(M,h1m*Re**3)/Bn
          dB=(Bp-Bm)/(mu_p-mu_m)
          return
         END SUBROUTINE deriv2 

         SUBROUTINE deriv(M,Bn,Re,ro,theta,dtheta,dB)
          implicit none
          double precision :: c3,sint1,sint2,cost1,cost2
          double precision :: ro,theta,mu_p,mu_m,Re
          double precision :: theta1,theta2,r1,r2
          double precision :: mu2,h1p,h1m,dtheta
          double precision :: Bp,Bm,dB,h_mu2,Bfield,M,Bn
          external mu2,h_mu2,Bfield
          theta1=theta-dtheta
          theta2=theta+dtheta

          sint1=sin(theta1)
          sint2=sin(theta2)
          cost1=cos(theta1)
          cost2=cos(theta2)
          c3=Re*Re*Re
          
          r1=ro*sint1*sint1
          r2=ro*sint2*sint2
          
          mu_p=mu2(r1,cost1)
          mu_m=mu2(r2,cost2)

          h1p=h_mu2(r1,cost1)
          h1m=h_mu2(r2,cost2)

          Bp=Bfield(M,h1p*c3)/Bn
          Bm=Bfield(M,h1m*c3)/Bn

          dB=(Bp-Bm)/(mu_p-mu_m)
          return
         END SUBROUTINE deriv 

        SUBROUTINE root2(mu,nu,r_old,r,theta)
        implicit none
        double precision :: a,b,r,mu,nu
        double precision :: theta,r_old
        complex*16 :: cterm1,cterm2,cterm3
        complex*16 :: x(4),temp1,temp2,temp3
        a=mu**2
        b=nu
        cterm1=4.*2.**(1./3.)
        cterm2=(27.*a*b**2+sqrt(6912.*a**3+729.*a**2*b**4))**(1./3.)
        cterm3=(3.*a*2**(1./3.))
        temp1=-cterm1/cterm2+cterm2/cterm3
        temp2=cterm1/cterm2-cterm2/cterm3- &
           (2.*b)/(a*sqrt(temp1))
!        x(1)=sqrt(temp1)/2.-
!     &     sqrt(temp2)/2.
!        x(2)=sqrt(temp1)/2.+sqrt(temp2)/2.
        temp3=cterm1/cterm2-cterm2/cterm3+ &
           (2.*b)/(a*sqrt(temp1))
!        x(3)=-sqrt(temp1)/2.-sqrt(temp3)/2.
        x(4)=-sqrt(temp1)/2.+sqrt(temp3)/2.
!        do i=1,4
!           diff(i)=abs(real(x(i))-r_old)
!        enddo
!        diff_min=min(diff(1),diff(2),diff(3),diff(4))
!        do i=1,4
!           if (diff(i).eq.diff_min)then
!              r=real(x(i))
!              goto 10
!           endif
!        enddo
! 10     continue
        r=dreal(x(4))
        theta=acos(r*r*mu)
        return
       END SUBROUTINE root2
  
        SUBROUTINE rootfind(mu,nu,r,theta)
        use param
        implicit none
        double precision :: a,b,r,mu,nu
        double precision :: theta,c1,c2
        complex*16 :: cterm1,cterm2,cterm3
        complex*16 :: x(4),temp1,temp2,temp3
        if(abs(mu).le.1.d-50)then
          r=1.d0/nu
          theta=pi/2.d0
          goto 100
        endif
        c1=1.d0/3.d0
        a=mu*mu
        b=nu
        cterm1=4.d0*2.d0**c1
        cterm2=(27.d0*a*b*b+ &
               sqrt(6912.d0*a*a*a+729.d0*a*a*b*b*b*b))**c1
        cterm3=(3.d0*a*2**c1)
        temp1=-cterm1/cterm2+cterm2/cterm3
        c2=sqrt(temp1)
!        temp2=cterm1/cterm2-cterm2/cterm3-
!     &     (2.d0*b)/(a*c2)
        temp3=cterm1/cterm2-cterm2/cterm3+ &
           (2.d0*b)/(a*c2)
        x(4)=-0.5d0*c2+0.5d0*sqrt(temp3)
        r=dreal(x(4))
        theta=acos(r*r*mu)
 100    continue
        return
        END SUBROUTINE rootfind

        SUBROUTINE rootfind2(mu,nu,r,theta)
        implicit none
        double precision :: a,b,r,mu,nu,pi
        double precision :: theta,c1,c2
        complex*16 :: cterm1,cterm2,cterm3
        complex*16 :: x(4),temp1,temp2,temp3
        pi=4.d0*atan(1.d0)
        if(abs(mu).le.1.d-50)then
          r=1.d0/nu
          theta=pi/2.d0
          goto 100
        endif
        c1=1./3.
        a=mu*mu
        b=nu
        cterm1=4.*2.**c1
        cterm2=(27.*a*b*b+sqrt(6912.*a*a*a+729.*a*a*b*b*b*b))**c1
        cterm3=(3.*a*2**c1)
        temp1=-cterm1/cterm2+cterm2/cterm3
        c2=sqrt(temp1)
!        temp2=cterm1/cterm2-cterm2/cterm3-
!     &     (2.*b)/(a*c2)
        temp3=cterm1/cterm2-cterm2/cterm3+ &
           (2.*b)/(a*c2)
        x(4)=-0.5*c2+0.5*sqrt(temp3)
        r=dreal(x(4))
        theta=acos(r*r*mu)
 100    continue
        return
        END SUBROUTINE rootfind2

        FUNCTION ezt(theta)
        implicit none
        double precision :: theta,ezt
        ezt=1.e-19*(sin(theta))**14
        return
        END FUNCTION ezt
