	  SUBROUTINE density2()
          use indeces
          use param
          use mesh
          use parparam
          implicit none
          double precision :: temp,norm
          character(len=7) :: temp2
          integer :: i,j,k,k_eq,k_ion
          k=1+(n1/2-1)*n2
          norm=den_grid0(k)
          write(*,*) 'norm',norm
!          call filter_r(n2,n1,m2,den_grid0,den_grid0)
!          call filter_p(n2,n1,m1,den_grid0,den_grid0)
          open(unit=2,file='density.dat',status='unknown')
          do i=1,n2
             do j=1,n1
                k=i+(j-1)*n2
                k_eq=n2/2+(n1/2-1)*n2
                k_ion=n2/2+(3-1)*n2
                rho(k)=rho_eq*den_grid0(k)/den_grid0(k_ion)
!               write(2,55) x(i,j),z(i,j),rho(k)
                write(2,55) x2(i),x1(j),rho(k)
             enddo
             write(2,*)
          enddo
          close(2) 
 10       format(a7)
 55       format(3e13.5)
          return
          END SUBROUTINE density2
