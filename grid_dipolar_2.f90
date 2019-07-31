	  SUBROUTINE grid()
! set up grid equally spaced along field line
! change from grid_dipolar6.f -> cludge to make x1 points more symmetric about
! equator.
          use indeces
          use param
          use mesh
          use fields
          use parparam
          implicit none
          integer :: i,j,k
          double precision :: ro,rho_mu_nu,r_temp,h_nu,h_mu,h_phi
          double precision :: lambda_max1,dlambda1,lambda
          double precision :: mu,nu
          double precision :: length,ro_temp,h1_temp,length2
          double precision :: lambda1,lambda2,delta
          double precision :: Ax,Az,Bx,Bz,AdotB,absA,absB
          double precision :: beta,width_t,len_t,dx1
          double precision :: theta_i1,theta_in1,Ep_test
          double precision :: s(n1),s_temp,s_i1,s_total,ds
          double precision :: theta_temp(n1),theta_t,Bfield
          external length,length2,h_nu,Ep_test,Bfield
          ro=(ro1+ro2)/2
          x2(1)=1./ro
          theta_in1=asin((x2(1)*r_ionosphere)**0.5)
!         write(*,*) ro,ro1,ro2
!         write(*,*) theta_in1,r_ionosphere,x2(1)
          theta_i1=pi-theta_in1
          s_i1=length(ro,-1.0*pi/2.,pi/2.-theta_i1)
          s_total=length(ro,pi/2.-theta_i1,pi/2.-theta_in1)
          IF(CENTRED)THEN
            ds=s_total/float(n1-5)
            s(1)=s_i1-2.d0*ds
          ELSE
            ds=s_total/float(n1-4)
            s(1)=s_i1-1.5d0*ds
          ENDIF
          do i=2,n1
             s(i)=s(i-1)+ds
          enddo

          do i=1,n1
             do theta_t=pi,0.,-pi/float(100000)
                s_temp=length(ro,-1.0*pi/2.,pi/2.-theta_t)
                if(abs(s(i)-s_temp).lt.0.0002)then
                  theta_temp(i)=theta_t
!                   write(*,*) i,s(i),s_temp,theta(i)
                  goto 50
                endif
             enddo
 50          continue
          enddo


          do i=1,n1
             r_temp=ro*sin(theta_temp(i))**2
             x1(i)=mu(r_temp,theta_temp(i))
          enddo

          dro=(ro2-ro1)/float(n2-1)
          ro=ro1
          do i=1,n2
            x2(i)=1./ro
            ro=ro+dro
          enddo
! set particle boundaries as midway between first to mu values (and last two).
!
          x1p1=x1(2)
          x1pn1=x1(n1-1)
          delta=dro*0.1
          x2p1=x2(1)+0.0000001
          x2pn2=x2(n2)-0.0000001
          open(unit=2,file='grid.dat',status='unknown')
          open(unit=3,file='Bo.dat',status='unknown')
          open(unit=10,file='h1.dat',status='unknown')
          open(unit=11,file='h2.dat',status='unknown')
          open(unit=12,file='h3.dat',status='unknown')
          do i=1,n2
            do j=1,n1
                k=i+(j-1)*n2
                r(i,j)=rho_mu_nu(x1(j),x2(i))
                theta(i,j)=acos(r(i,j)**2*x1(j))
                lambda=pi/2.-theta(i,j)
                x(i,j)=r(i,j)*cos(lambda)
                z(i,j)=r(i,j)*sin(lambda)
                h1(k)=h_mu(r(i,j),theta(i,j))
                h2(k)=h_nu(r(i,j),theta(i,j))
                h3(k)=h_phi(r(i,j),theta(i,j))
                Bo(k)=Bfield(M,h1(k)*Re**3)
                if(ETEST)then
                  k=i+(j-1)*n2
                  e1_1(k)=Ep_test(theta(i,j))
                  e1_2(k)=e1_1(k)
                endif
                write(2,*) x(i,j),z(i,j)
!                write(3,55) x(i,j),z(i,j),Bo(k)
                write(3,55) x2(i),x1(j),Bo(k),theta(i,j)
                write(10,55) x(i,j),z(i,j),h1(k)
                write(11,55) x(i,j),z(i,j),h2(k)
                write(12,55) x(i,j),z(i,j),h3(k)
             enddo
             write(3,*)
             write(10,*)
             write(11,*)
             write(12,*)
          enddo
          close(2) 
          close(3)
          close(11)
          close(12)
  
! set normalization constants for B and rho
!
          Bn=10.d-9
          rho_n=1.d5*mp

! determine dBo (already normalized)
!
          open(unit=3,file='dBo.dat',status='unknown')
          do i=1,n2
            do j=2,n1-1
                k=i+(j-1)*n2
                dBo(k)=(Bo(k+n2)-Bo(k-n2))/(x1(j+1)-x1(j-1))/Bn
                write(3,55) x(i,j),z(i,j),dBo(k)
            enddo
            write(3,*)
          enddo
          close(3)

! set  equatorial value of rho
!
          rho_eq=1.0d6*mp

! vn to the Alfven velocity given Bn and rho_n
!
          vn=Bn/sqrt(mu_o*rho_n)
          tn=L/vn
          en=Bn*vn
          jn=Bn/(mu_o*L)
          constE_mu=abs(e*tn*Bn/me)
          const_ec=1.d0/(mu_o*epsilon_o*vn*vn)
          mu_mn=me*vn*vn/Bn
          const_jp=rho_n*e*vn/(mp*jn)
!          const_mu=L/en*(4.d0*pi)/(c*c)*(e/me)*Bn*mu_mn*rho_n/mp
!          const_sm=L*L/en*(4.d0*pi)/(c*c)*(e/L)*vn*vn*rho_n/mp
          const_mu=e*mu_o*L*vn*rho_n/(Bn*mp)
          const_sm=e*mu_o*L*vn*rho_n/(Bn*mp)
          open(unit=2,file='x2.dat',status='unknown')
          do i=1,n2
             write(2,*) i,x2(i)
          enddo
          close(2) 
          
          open(unit=2,file='x1.dat',status='unknown')
          do j=1,n1
             write(2,*) x1(j),theta(n2/2,j),z(n2/2,j)
          enddo
          close(2) 
          open(unit=3,file='length.dat',status='unknown')
          open(unit=4,file='angle.dat',status='unknown')
!          open(unit=10,file='area.dat',status='unknown')
          open(unit=11,file='width.dat',status='unknown')
          do i=1,n2-1
             do j=1,n1-1
                k=i+(j-1)*n2
                ro_temp=r(i,j)/(sin(theta(i,j)))**2
                lambda1=atan(z(i,j)/x(i,j))
                lambda2=atan(z(i,j+1)/x(i,j+1))
!                len(i,j)=length(ro_temp,lambda2,lambda1)
!                temp=abs(h1(i,j)+h1(i+1,j+1)
                width(i,j)=abs((h2(k)+h2(k+1))/2.*(x2(i)-x2(i+1)))
                 len_t=sqrt((x(i,j+1)-x(i,j))**2+(z(i,j+1)-z(i,j))**2)
!               len(i,j)=sqrt((x(i,j+1)-x(i,j))**2+(z(i,j+1)-z(i,j))**2)
                 len(i,j)=abs((h1(k)+h1(k+n2))/2.*(x1(j)-x1(j+1)))
!                 width_t=sqrt((x(i+1,j)-x(i,j))**2+(z(i+1,j)-z(i,j))**2)
                 Ax=x(i+1,j)-x(i,j)
                 Az=z(i+1,j)-z(i,j)
                 Bx=x(i,j+1)-x(i,j)
                 Bz=z(i,j+1)-z(i,j)
                 AdotB=Ax*Bx+Az*Bz
                 absA=sqrt(Ax*Ax+Az*Az)
                 absB=sqrt(Bx*Bx+Bz*Bz)
                 beta=acos(AdotB/(absA*absB))
                 alpha(i,j)=pi-beta
                 alpha(i,j)=0.d0
!                 area(i,j)=len(i,j)*width(i,j)*cos(alpha(i,j))
                 cosalpha(i,j)=cos(alpha(i,j))
                 max_area=max(max_area,area(i,j))
                 write(3,*) x(i,j),z(i,j),len(i,j),len_t
                 write(4,55) x(i,j),z(i,j),alpha(i,j)
                 write(10,55) x(i,j),z(i,j),area(i,j)
                 write(11,55) x(i,j),z(i,j),width(i,j)
             enddo
             write(3,*) 
             write(4,*) 
!             write(10,*) 
             write(11,*) 
          enddo
          close(3)
          close(4)
!          close(10)
          close(11)
          do i=1,n2
             do j=1,n1
!                area_ratio(i,j)=max_area/area(i,j)
                area_ratio(i,j)=1.d0
             enddo
          enddo
          if(ETEST)then
          open(unit=3,file='Ep_test.dat',status='unknown')
          do i=1,n2
             do j=1,n1
                k=i+(j-1)*n2
                write(3,55) x(i,j),z(i,j),e1_2(k)
             enddo
             write(3,*)
          enddo
          close(3)
          endif
          lambda_b=acos(sqrt(5.d0/ro))
 55       format(4e13.5)
          return
          END SUBROUTINE grid
 

          FUNCTION h_mu(r,theta)
          implicit none
          double precision :: r,theta,h_mu,c2
          c2=cos(theta)
          h_mu=r*r*r/sqrt((1+3.*c2*c2))
          return
          END FUNCTION h_mu
 
          FUNCTION h_mu2(r,cost)
          implicit none
          double precision :: r,cost,h_mu2
          h_mu2=r*r*r/sqrt((1+3.*cost*cost))
          return
          END FUNCTION h_mu2

          FUNCTION h_nu(r,theta)
          implicit none
          double precision :: r,theta,h_nu,c1
          c1=cos(theta)
          h_nu=r*r/(sin(theta)*sqrt(1.+3.*c1*c1))
          return
          END FUNCTION h_nu


          FUNCTION h_nu2(r,cost,sint)
          implicit none
          double precision :: r,theta,h_nu2,cost,sint
          h_nu2=r*r/(sint*sqrt(1.+3.*cost*cost))
          return
          END FUNCTION h_nu2

          FUNCTION h_phi(r,theta)
          implicit none
          double precision :: r,theta,h_phi
          h_phi=r*sin(theta)
          return
          END FUNCTION h_phi

          FUNCTION h_mu_xy(x,y,z)
          implicit none
          double precision :: x,y,z,h_mu_xy
          h_mu_xy=(x*x+y*y+z*z)**2/sqrt(x*x+y*y+4.*z*z)
          return
          END FUNCTION h_mu_xy
       
          FUNCTION h_nu_xy(x,y,z)
          implicit none
          double precision :: x,y,z,h_nu_xy
          h_nu_xy=(x*x+y*y+z*z)**2/(sqrt(x*x+y*y)*sqrt(x*x+y*y+4.*z*z))
          return
          END FUNCTION h_nu_xy

          FUNCTION h_phi_xy(x,y,z)
          implicit none
          double precision :: x,y,z,h_phi_xy
          h_phi_xy=sqrt(x*x+y*y)
          return
          END FUNCTION h_phi_xy
          
          FUNCTION rho_mu_nu(mu,nu)
          implicit none
          double precision :: mu,nu,term1,term2,term3,term4,term5
          double precision :: den1,rho_mu_nu,temp
          temp=(27.*mu**2*nu**2+ &
                sqrt(6912.*mu**6+729.*mu**4*nu**4))**(1./3.)              
          term1=-4.0*2**(1./3.)/temp
          term2=temp/(3.*2.**(1./3.)*mu**2)
          term3=4.0*2.**(1./3.)/temp
          term4=temp/(3.*2.**(1./3.)*mu**2)
          den1=-4.*2.**(1./3.)/temp+temp/(3.*2.**(1./3.)*mu**2)
          term5=2.*nu/(mu**2*sqrt(den1))         
          rho_mu_nu= &
              -0.5*sqrt(term1+term2)+0.5*sqrt(term3-term4+term5)
          return
          END FUNCTION rho_mu_nu      

          FUNCTION mu(r,theta)
          implicit none
          double precision :: r,theta,mu
          mu=cos(theta)/(r*r)
          return
          END FUNCTION mu

          FUNCTION mu2(r,cost)
          implicit none
          double precision :: r,cost,mu2
          mu2=cost/(r*r)
          return
          END FUNCTION mu2

          FUNCTION nu(r,theta)
          implicit none
          double precision :: r,theta,nu
          nu=sin(theta)*sin(theta)/r
          return
          END FUNCTION nu


        FUNCTION length2(ro,a,b)
        implicit none
        double precision :: ro,lambda,length2,asinh,a,b,pi
        pi=4.0*atan(1.0)
! asinh= ln(x+sqrt(x*x+1))
        length2=-(ro*asinh(sqrt(3.0)*sin(a)))/(2.0*sqrt(3.0))+ &
               ro*asinh(sqrt(3.0)*sin(b))/(2.0*sqrt(3.0))- &
               ro*sin(a)*sqrt(1+3.0*(sin(a))**2)/2.0+ &
               ro*sin(b)*sqrt(1+3.0*(sin(b))**2)/2.0
        return
        END FUNCTION length2

        FUNCTION length(ro,a,b)
        implicit none
        double precision :: ro,lambda,length,asinh,a,b,pi
        double precision :: c1,c2,c3
        c1=sqrt(3.0)
        c2=sin(a)
        c3=sin(b)
        pi=4.0*atan(1.0)
! asinh= ln(x+sqrt(x*x+1))
        length=-(ro*asinh(c1*c2))/(2.0*c1)+ &
               ro*asinh(c1*c3)/(2.0*c1)- &
               0.5*ro*c2*sqrt(1+3.0*c2*c2)+ &
               0.5*ro*c3*sqrt(1+3.0*c3*c3)
        return
        END FUNCTION length   

        FUNCTION asinh(x)
        implicit none
        double precision :: x,asinh
        asinh=log(x+sqrt(x*x+1))
        return
        END FUNCTION asinh

        FUNCTION Ep_test(theta)
        implicit none
        double precision :: Ep_test,theta
        Ep_test=1.e-2*(sin(theta))**14
        return
        END FUNCTION Ep_test
