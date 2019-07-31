          SUBROUTINE initial_flr()
          use indeces
          use param
          use mesh
          use fields
          implicit none
          double precision :: temp,d2,f,fmax
          double precision :: profile_max
          integer :: i,j,k,k2
          open(unit=2,file='initial.dat',status='unknown')
!          d1=0.1
!          d1=0.10
! use d1=0.3 for r_ionosphere=5
!         d1=0.3
!         fwhm=0.25/L
          fwhm=fwhm/L
          d2=fwhm/(2.d0*sqrt(-2.d0*log(1.d0/2.d0)))
!          Amp=10.d4/vn
          Amp=Amp/vn
          temp=(ro1+ro2)/2.
          fmax=f(theta(1,n1/2))
          profile_max=0.1014
          L_perp=1.0
          do i=1,n2
             do j=3,n1-2
#ifdef POLOIDAL
                k=i+(j-1)*n2
#else
                k=i+(j-1)*n2+n21
#endif
                y(k)=-Amp*f(theta(1,j))/fmax* &
                      exp(-((1./x2(i)-1./x2(n2/2))/L)**2/(2.*d2*d2))* &
                      sin(2*pi/L_perp*(1./x2(i)-1./x2(n2/2)))
!                     exp(-((1./x2(i)-1./profile_max)/L)**2/(2.*d2*d2))
!             if(x2(i).lt.profile_max)then
!               y(k)=-Amp*f(theta(1,j))/fmax
!             endif
              if(j.eq.3.or.j.eq.n1-2)y(k)=0.d0
              y0(k)=y(k)
             enddo
          enddo

          j=2
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j+1-1)*n2
             if(CENTRED)k2=i+(j+2-1)*n2
#ifdef POLOIDAL
             y(k)=-y(k2)
#else
             y(k+n21)=-y(k2+n21)
#endif
          enddo

          j=1
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j+3-1)*n2
             if(CENTRED)k2=i+(j+4-1)*n2
#ifdef POLOIDAL
             y(k)=-y(k2)
#else
             y(k+n21)=-y(k2+n21)
#endif
          enddo

          j=n1-1
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j-1-1)*n2
             if(CENTRED)k2=i+(j-2-1)*n2
#ifdef POLOIDAL
             y(k)=-y(k2)
#else
             y(k+n21)=-y(k2+n21)
#endif
          enddo

          j=n1
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j-3-1)*n2
             if(CENTRED)k2=i+(j-4-1)*n2
#ifdef POLOIDAL
             y(k)=-y(k2)
#else
             y(k+n21)=-y(k2+n21)
#endif
          enddo

          do i=1,n2
            do j=1,n1
              k=i+(j-1)*n2
              write(2,*) x(i,j),z(i,j),y(k+n21),y(k)
            enddo
            write(2,*)
          enddo
          close(2)
          
          open(unit=2,file='initial_eq.dat',status='unknown')
          do i=1,n2
             k=i+(n1/2-1)*n2
             write(2,*) 1./x2(i),y(k+n21),y(k)
          enddo
          close(2)
          return
          END SUBROUTINE initial_flr

          FUNCTION f(theta)
          use indeces
          use param
          double precision :: theta,lambda,s,so,lambda_i 
          double precision :: nu11,nu13,f
          lambda=pi/2.d0-theta
!          lambda_i=71.54d0*pi/180.d0
          lambda_i=63.55d0*pi/180.d0
!          lambda_i=57.d0*pi/180.d0
          s=sin(lambda)+sin(lambda_i)
          so=2.d0*sin(lambda_i)
          nu11=0.393237d0
          nu13=0.318194d0
          f=cos(lambda)**3*(sin(pi*s/so)+(nu13/(9.d0*nu11-8.d0))* &
            sin(3.d0*pi*s/so))
          return
          END FUNCTION f

