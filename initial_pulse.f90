          SUBROUTINE initial_pulse()
          use indeces
          use param
          use mesh
          use fields
          implicit none
          double precision :: temp,sigma_par,sigma_perp,e2_temp(n21)
          double precision :: e3_temp(n21),pot(n21)
          integer :: i,j,k,k2
          open(unit=2,file='initial.dat',status='unknown')
          sigma_perp=fwhm/(2.d0*sqrt(-2.d0*log(1.d0/2.d0)))
          sigma_par=fwhm_par/(2.d0*sqrt(-2.d0*log(1.d0/2.d0)))
          Amp=Amp/vn
          temp=(ro1+ro2)/2.
          L_perp=L_perp/L
          do i=1,n2
             do j=3,n1-2
#ifdef POLOIDAL
                k=i+(j-1)*n2
#else
                k=i+(j-1)*n2+n21
#endif
                y(k)=-Amp* &
                   exp(-((pi/2.-theta(n2/2,j)-0.0d0))**2/ &
                   (2.*sigma_par**2))* &
                   exp(-((1./x2(i)-1./x2(n2/2)))**2/(2.*sigma_perp**2))* &
                   sin(2.*pi/x2(i)/L_perp)
                if(j.eq.3.or.j.eq.n1-2)y(k)=0.d0
                if(abs(y(k)).lt.1.d-99) y(k)=0.d0
                y(k+n21)=y(k+n21)/(1+ &
                       float(igyro)*rho_i(k)**2*(2.*pi/L_perp)**2)
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

          call eperp(y,e2_temp,e3_temp)

          do i=1,n2
            do j=1,n1
              k=i+(j-1)*n2
              write(2,55) x(i,j),z(i,j),y(k+n21),y(k)
            enddo
            write(2,*)
          enddo
          close(2)
          open(unit=2,file='initial_eq.dat',status='unknown')
          do i=1,n2
             k=i+(n1/2-1)*n2
             write(2,55) 1./x2(i),y(k+n21)
          enddo
          close(2)
  55      format(5e15.7)
          return
          END SUBROUTINE initial_pulse
