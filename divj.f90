          SUBROUTINE div_j(j1,j2,j3,j1x1,j2x2,j3x3,divj)
          use indeces
          use param
          use mesh
          implicit none
          integer :: i,j,k
          double precision :: j1(*),j2(*),j3(*),divj(*)
          double precision :: j1x1(*),j2x2(*),j3x3(*)
          do i=1,n21
             divj(i)=0.d0
             j1x1(i)=0.d0
             j2x2(i)=0.d0
             j3x3(i)=0.d0
          enddo
          do j=3,n1-2
          do i=2,n2-1
                k=i+(j-1)*n2
                j1x1(k)=(h1(k+n2)*j1(k+n2)- &
                 h1(k-n2)*j1(k-n2))/ &
                 (x1(j+1)-x1(j-1))
!                if(.not.CENTRED.
!     &          and.(j.eq.3.or.j.eq.n1-2))j1x1=0.0
!                if(j.eq.3.or.j.eq.n1-2)j1x1=0.0
                 j2x2(k)=(h1(k+1)*h3(k+1)*j2(k+1)- &
                 h1(k-1)*h3(k-1)*j2(k-1))/ &
                 (x2(i+1)-x2(i-1)) 
                 j3x3(k)=mphi*h1(k)*h2(k)*j3(k)
                 divj(k)=1.d0/(h1(k)*h2(k)*h3(k))*( &
                         j1x1(k)+j2x2(k)+j3x3(k))
             enddo
          enddo           
!
          if(ifilter.eq.1)call filterd(n2,n1,divj,divj)
          return
          END SUBROUTINE div_j
