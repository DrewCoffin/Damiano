	  SUBROUTINE field_update()
          use indeces
          use param
          use mesh
          use fields
          implicit none
          integer :: k
          do k=1,n21
             j1_0(k)=j1_1(k)
             j1_1(k)=j1_2(k)
             e1_0(k)=e1_1(k)
             e1_1(k)=e1_2(k)
          enddo

          return
          END SUBROUTINE field_update
