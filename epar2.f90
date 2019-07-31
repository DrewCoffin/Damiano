	  SUBROUTINE epar2()
          use indeces
          use param
          use mesh
          use fields
          implicit none
          integer :: k
          do k=1,n21
             e1_2(k)=e1_0(k)+2.*dt*e1_1(k)
          enddo
          return
          END SUBROUTINE epar2
