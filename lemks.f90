	  SUBROUTINE le_mks()
          use indeces
          use param
          use mesh
          implicit none
          double precision :: n_e
          integer :: i,j,k
          open(unit=2,file='le.dat',status='unknown')
          do i=1,n2
             do j=3,n1-2
                k=i+(j-1)*n2
                n_e=rho(k)/mp
                le2(k)=me/(mu_o*n_e*e*e)
                write(2,*) x(i,j),z(i,j),sqrt(le2(k)),sqrt(le2(k))/L
             enddo
             write(2,*)
          enddo
          close(2) 

          return
          END SUBROUTINE le_mks
