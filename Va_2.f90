	  SUBROUTINE Va_profile2()
          use indeces
          use param
          use mesh
          implicit none
          double precision :: temp
          integer :: i,j,k
          open(unit=2,file='va2.dat',status='unknown')
          do i=1,n2
             do j=1,n1
                k=i+(j-1)*n2
                temp=cos(theta(i,j))
                Bo(k)=Bo(k)*Bn
                Va(k)=Bo(k)/sqrt(mu_o*rho(k))
                write(2,55) x(i,j),z(i,j),Va(k)
             enddo
             write(2,*)
          enddo
          close(2) 
 55       format(3e13.5)
          return
          END SUBROUTINE Va_profile2
