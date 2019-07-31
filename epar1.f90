	  SUBROUTINE epar1()
          use indeces
          use param
          use mesh
          use fields
          implicit none
          integer :: i,j,k,k2
          do j=3,n1-2
          do i=1,n2
                k=i+(j-1)*n2
                e1_1(k)=le2(k)*(j1_2(k)-j1_0(k))/(2*dt)
             enddo
          enddo

       j=2
       do i=1,n2
          k=i+(j-1)*n2
          k2=i+(j+2-1)*n2
          e1_1(k)=e1_1(k2)
       enddo
 
       j=1
       do i=1,n2
          k=i+(j-1)*n2
          k2=i+(j+4-1)*n2
          e1_1(k)=e1_1(k2)
       enddo
 
 
       j=n1-1
       do i=1,n2
          k=i+(j-1)*n2
          k2=i+(j-2-1)*n2
          e1_1(k)=e1_1(k2)
       enddo
 
       j=n1
       do i=1,n2
          k=i+(j-1)*n2
          k2=i+(j-4-1)*n2
          e1_1(k)=e1_1(k2)
       enddo
       return
       END SUBROUTINE epar1
