          SUBROUTINE density1()
          use indeces
          use param
          use mesh
          implicit none
          double precision :: temp,ro_t,dp
          integer :: i,j,k
          external dp
          open(unit=2,file='density.dat',status='unknown')
          do i=1,n2
             ro_t=1/x2(i)
             do j=1,n1
                   k=i+(j-1)*n2
                   rho(k)=rho_eq*dp(x1(j))/dp(x1(3))
!                  write(2,*) x(i,j),z(i,j),rho(k)
                   write(2,*) x2(i),x1(j),rho(k)
             enddo
             write(2,*)
          enddo
          close(2)
          return
          END SUBROUTINE density1
