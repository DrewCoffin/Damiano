          SUBROUTINE div_b(y,divb)
          use indeces
          use param
          use mesh
          implicit none
          integer :: i,j,k,kb1,kb2,kb3
          double precision :: y,divb
          dimension y(*),divb(*)

          do j=2,n1-2
          do i=2,n2-2
                k=i+n2*(j-1)
                kb1=k+n21*2
                kb2=k+n21*3
                kb3=k+n21*4
                divb(k)=1.d0/(h1(k)*h2(k)*h3(k))*( &
                (h2(k+n2)*h3(k+n2)*y(kb1+n2)- &
                 h2(k-n2)*h3(k-n2)*y(kb1-n2))/ &
                 (x1(j+1)-x1(j-1))+ &
                (h1(k+1)*h3(k+1)*y(kb2+1)- &
                 h1(k-1)*h3(k-1)*y(kb2-1))/ &
                 (x2(i+1)-x2(i-1))+ &
                 mphi*h1(k)*h2(k)*y(kb3))
             enddo
          enddo
          return
          END SUBROUTINE div_b
