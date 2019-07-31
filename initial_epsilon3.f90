          SUBROUTINE initial_epsilon3()
          use indeces
          use param
          use mesh
          use fields
          implicit none
          double precision :: temp,d1,d2,Amp,fwhm,e2_temp(n21)
          double precision :: f_e3,fmax
! pot=potential - only valid when rho_i=0
          double precision :: e3_temp(n21),pot(n21)
          double precision :: L_perp,f2(n1)
          integer :: i,j,k,k2,dummy
          open(unit=2,file='initial.dat',status='unknown')
! use d1=0.3 for r_ionosphere=5
          d1=0.45
          fwhm=0.1/L
          d2=fwhm/(2.d0*sqrt(-2.d0*log(1.d0/2.d0)))
          Amp=-3.0d5/vn
          temp=(ro1+ro2)/2.
          fmax=f_e3(theta(1,n1/2))
          L_perp=0.1
          do i=1,n2
             do j=3,n1-2
                k=i+(j-1)*n2
!                y(k+n21)=-Amp*
!                y(k+n21)=-Amp*f(theta(1,j))/fmax*
                y(k+n21)=-Amp*f_e3(theta(1,j))/fmax* &
!     &                exp(-((pi/2.-theta(n2/2,j)-0.0d0))**2/d1**2)*
                      exp(-((1./x2(i)-1./x2(n2/2))/L)**2/(2.*d2*d2))* &
                      sin(2*pi/x2(i)/L_perp)
                pot(k)=-Amp*Bo(k)/(2.*pi/L_perp)* &
                      exp(-((pi/2.-theta(n2/2,j)-0.0d0))**2/d1**2)* &
                      exp(-((1./x2(i)-1./x2(n2/2))/L)**2/(2.*d2*d2))* &
                      cos(2*pi/x2(i)/L_perp)
              if(j.eq.3.or.j.eq.n1-2)y(k+n21)=0.d0
              if(abs(y(k+n21)).lt.1.d-99) y(k+n21)=0.d0
             y(k+n21)=y(k+n21)/(1+ &
                      float(igyro)*rho_i(k)**2*(2.*pi/L_perp)**2)
             y0(k+n21)=y(k+n21)
             enddo
          enddo

          j=2
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j+1-1)*n2
             if(CENTRED)k2=i+(j+2-1)*n2
             y(k+n21)=-y(k2+n21)
          enddo

          j=1
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j+3-1)*n2
             if(CENTRED)k2=i+(j+4-1)*n2
             y(k+n21)=-y(k2+n21)
          enddo

          j=n1-1
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j-1-1)*n2
             if(CENTRED)k2=i+(j-2-1)*n2
             y(k+n21)=-y(k2+n21)
          enddo

          j=n1
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j-3-1)*n2
             if(CENTRED)k2=i+(j-4-1)*n2
             y(k+n21)=-y(k2+n21)
          enddo

          call eperp(y,e2_temp,e3_temp)

          do i=1,n2
            do j=1,n1
              k=i+(j-1)*n2
              write(2,55) x(i,j),z(i,j),y(k+n21),e2_temp(k)
            enddo
            write(2,*)
          enddo
          close(2)
          open(unit=2,file='initial_eq.dat',status='unknown')
          do i=1,n2
             k=i+(n1/2-1)*n2
             write(2,55) 1./x2(i),y(k+n21),e2_temp(k),pot(k), &
             exp(-(1./x2(i)-10.0)**18/0.00000000001) 
          enddo
          close(2)
  55      format(5e13.5)
          return
          END SUBROUTINE initial_epsilon3

          FUNCTION f_e3(theta)
          use indeces
          use param
          double precision :: theta,lambda,s,so,lambda_i
          double precision :: nu11,nu13,f_e3
          lambda=pi/2.d0-theta
!          lambda_i=65.8d0*pi/180.d0
!          lambda_i=56.0d0*pi/180.d0
          lambda_i=55.0d0*pi/180.d0
          s=sin(lambda)+sin(lambda_i)
          so=2.d0*sin(lambda_i)
!          nu11=0.393237d0
!          nu13=0.318194d0
          nu11=0.376853
          nu13=0.31261
!          f=cos(lambda)**3*(sin(pi*s/so)+(nu13/(9.d0*nu11-8.d0))*
!     &      sin(3.d0*pi*s/so))
          f_e3=cos(lambda)**3*( &
            0.87942*sin(1.91759*s) - &
            1.27822*10**(-10)*sin(3.83518*s) + &
            1.10489*sin(5.75276*s) + &
            6.08677*10**(-11)*sin(7.67035*s) - &
            0.376789*sin(9.58794*s) - &
            1.22111*10**(-13)*sin(11.5055*s) + &
            0.0252817*sin(13.4231*s) &
!    &       0.891869*sin(1.89472*s) - 
!    &       3.38772*10**(-10)*sin(3.78945*s) +
!    &       1.09828*sin(5.68417*s) + 
!    &       1.81485*10**(-10)*sin(7.57889*s)-
!    &       0.379119*sin(9.47361*s)-
!    &       1.51996*10**(-13)*sin(11.3683*s)+
!    &       0.0268327*sin(13.2631*s)
          )
          return
          END FUNCTION f_e3

