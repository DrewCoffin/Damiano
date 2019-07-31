	  SUBROUTINE rho_i_calc()
          use indeces
          use param
          use mesh

          implicit none
          double precision :: C_i
          integer :: i,j,k

          C_i=mp*sqrt(T_i/mp)/(e*Bn*L)
          open(unit=2,file='rho_i.dat',status='unknown')
! Bo normalized
          do i=1,n2
             do j=3,n1-2
                k=i+(j-1)*n2
                rho_i(k)=C_i/Bo(k)
                write(2,*) x2(i),x1(j),rho_i(k) 
             enddo
             write(2,*)
          enddo
          close(2) 
          return
          END SUBROUTINE rho_i_calc
