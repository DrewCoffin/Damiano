          SUBROUTINE initial_pulse()
          use indeces
          use param
          use mesh
          use fields
          implicit none
          double precision :: temp,d1,d2,Amp,fwhm,e2_temp(n21)
! pot=potential - only valid when rho_i=0
          double precision :: e3_temp(n21),pot(n21)
          double precision :: L_perp
          integer :: i,j,k,k2
          open(unit=2,file='initial.dat',status='unknown')
! use d1=0.3 for r_ionosphere=5
          d1=0.05
          fwhm=0.05/L
          d2=fwhm/(2.d0*sqrt(-2.d0*log(1.d0/2.d0)))
          Amp=5.00d5/vn
          temp=(ro1+ro2)/2.
          L_perp=0.2
          do i=1,n2
             do j=3,n1-2
                k=i+(j-1)*n2
                y(k)=-Amp* &
                   exp(-((pi/2.-theta(n2/2,j)-0.0d0))**2/d1**2)* &
                   exp(-((1./x2(i)-1./x2(n2/2))/L)**2/(2.*d2*d2))* &
                   sin(2*pi/x2(i)/L_perp)
                pot(k)=-Amp*Bo(k)/(2.*pi/L_perp)* &
                      exp(-((pi/2.-theta(n2/2,j)-0.0d0))**2/d1**2)* &
                      exp(-((1./x2(i)-1./x2(n2/2))/L)**2/(2.*d2*d2))* &
                      cos(2*pi/x2(i)/L_perp)
              if(j.eq.3.or.j.eq.n1-2)y(k)=0.d0
              if(abs(y(k+n21)).lt.1.d-99) y(k+n21)=0.d0
!            y(k+n21)=y(k+n21)/(1+ &
!                     float(igyro)*rho_i(k)**2*(2.*pi/L_perp)**2)
             y0(k)=y(k)
             enddo
          enddo

          j=2
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j+1-1)*n2
             if(CENTRED)k2=i+(j+2-1)*n2
             y(k)=-y(k2)
          enddo

          j=1
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j+3-1)*n2
             if(CENTRED)k2=i+(j+4-1)*n2
             y(k)=-y(k2)
          enddo

          j=n1-1
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j-1-1)*n2
             if(CENTRED)k2=i+(j-2-1)*n2
             y(k)=-y(k2)
          enddo

          j=n1
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j-3-1)*n2
             if(CENTRED)k2=i+(j-4-1)*n2
             y(k)=-y(k2)
          enddo

          call eperp(y,e2_temp,e3_temp)

          do i=1,n2
            do j=1,n1
              k=i+(j-1)*n2
              write(2,55) x(i,j),z(i,j),y(k),e3_temp(k)
            enddo
            write(2,*)
          enddo
          close(2)
          open(unit=2,file='initial_eq.dat',status='unknown')
          do i=1,n2
             k=i+(n1/2-1)*n2
             write(2,55) 1./x2(i),y(k+n21),e3_temp(k),pot(k), &
             exp(-(1./x2(i)-10.0)**18/0.00000000001) 
          enddo
          close(2)
  55      format(5e13.5)
          return
          END SUBROUTINE initial_pulse
