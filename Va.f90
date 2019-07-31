	  SUBROUTINE Va_profile()
          use indeces
          use param
          use mesh
          implicit none
          double precision temp
          integer i,j,k
          open(unit=2,file='va1.dat',status='unknown')
          do i=1,n2
             do j=1,n1
                k=i+(j-1)*n2
                temp=cos(theta(i,j))
                Va(k)=Bo(k)/sqrt(mu_o*rho(k))
                gamma2(k)=1.d0
#ifdef JUPITER
                gamma2(k)=c*c/(c*c+Va(k)*Va(k))
#endif
             write(2,55) x2(i),x1(j),Va(k),gamma2(k)*Va(k)
             enddo
             write(2,*)
          enddo
          close(2) 
 55       format(4e13.5)
          return
          END SUBROUTINE Va_profile
