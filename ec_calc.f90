          SUBROUTINE ec_calc(Fc,ec)
          use indeces
          use param
          use mesh
          implicit none
          integer :: i,j,k
          double precision :: temp1(n2),temp2(n2)
          double precision :: Fc(*),ec(*)
          do j=1,n1
             do i=1,n2
                k=i+(j-1)*n2
                temp1(i)=h1(k)*h2(k)*h3(k)*Fc(k)
             enddo
             call RungeK(n2,x2,temp1,temp2)
             do i=1,n2
                k=i+(j-1)*n2
                ec(k)=temp2(i)/(h1(k)*h3(k))
             enddo
          enddo
          do i=1,n2
           IF(CENTRED)THEN
              ec(i+(1-1)*n2)=-ec(i+(5-1)*n2)
              ec(i+(2-1)*n2)=-ec(i+(4-1)*n2)
              ec(i+(n1-1)*n2)=-ec(i+(n1-4-1)*n2)
              ec(i+(n1-1-1)*n2)=-ec(i+(n1-3-1)*n2)
           ELSE
              ec(i+(1-1)*n2)=-ec(i+(4-1)*n2)
              ec(i+(2-1)*n2)=-ec(i+(3-1)*n2)
              ec(i+(n1-1)*n2)=-ec(i+(n1-3-1)*n2)
              ec(i+(n1-1-1)*n2)=-ec(i+(n1-2-1)*n2)
           ENDIF
          enddo
          return
          END SUBROUTINE ec_calc

          SUBROUTINE RungeK(n,x,F,y)
          integer :: n
          double precision :: x(n),F(n),y(n)
          double precision :: dydx(n)
          double precision :: k1,k2,k3,k4
!
! subroutine to solve dy/dx=F(x) using 
! a fifth order Runge-Kutte routine
!
          do i=1,n
             y(i)=0.d0
             dydx(i)=F(i)
          enddo
!          do i=2,n-1
!             dydx(i)=(F(i+1)-F(i-1))/(x(i+1)-x(i-1))
!          enddo
          do i=1,n-1
             k1=dydx(i)*(x(i+1)-x(i))
             k2=(dydx(i)+dydx(i+1))/2.*(x(i+1)-x(i))
             k3=k2
             k4=dydx(i+1)*(x(i+1)-x(i))
             y(i+1)=y(i)+1.0/6.0*(k1+2.*k2+2.*k3+k4)
          enddo
          return
          END SUBROUTINE RungeK
