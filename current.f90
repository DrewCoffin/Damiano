          SUBROUTINE current(y,j1,j2,j3)
! subroutine to determine all component of current density. Note that this 
! routine only works for k_3=0.
! Peter Damiano, last modified Dec 6 2005.
!
          use indeces
          use param
          use mesh
          implicit none
          double precision :: term1,term2,y(*),j1(*)
          double precision :: j2(*),j3(*)
          integer :: i,j,k,kb1,kb2,kb3,k2
          do i=1,n21
             j1(k)=0.d0
             j2(k)=0.d0
             j3(k)=0.d0
          enddo
          do j=3,n1-2
             do i=2,n2-1
                   k=i+n2*(j-1)
                   kb1=k+n21*2
                   kb2=k+n21*3
                   kb3=k+n21*4
                   term1=(h3(k+1)*y(kb3+1)-h3(k-1)*y(kb3-1))/ &
                        (x2(i+1)-x2(i-1))
                   term2=mphi*h2(k)*y(kb2)
                   j1(k)=1.d0/(h2(k)*h3(k))*(term1-term2)
                   term1=mphi*h1(k)*y(kb1)
                   term2=(h3(k+n2)*y(kb3+n2)-h3(k-n2)*y(kb3-n2))/ &
                         (x1(j+1)-x1(j-1))
                   j2(k)=1.d0/(h1(k)*h3(k))*(term1-term2)
                   term1=(h2(k+n2)*y(kb2+n2)-h2(k-n2)*y(kb2-n2))/ &
                         (x1(j+1)-x1(j-1))
                   term2=(h1(k+1)*y(kb1+1)- &
                          h1(k-1)*y(kb1-1))/ &
                          (x2(i+1)-x2(i-1))
                   j3(k)=1.d0/(h1(k)*h2(k))*(term1-term2)
                   if(.not.CENTRED.and. &
                   (j.eq.3.or.j.eq.n1-2))j2(k)=0.d0
             enddo
! radial boundary conditions
!
             k=1+n2*(j-1)
             k2=3+n2*(j-1)
!            j1(k)=-j1(k2)
!            j2(k)=h1(k2)*h3(k2)*j2(k2)/(h1(k)*h3(k))
!            j3(k)=j3(k2)
             j1(k)=j1(k2)
             j2(k)=j2(k2)
             j3(k)=j3(k2)
             k=n2+n2*(j-1)
             k2=n2-2+n2*(j-1)
!            j1(k)=-j1(k2) 
!            j2(k)=h1(k2)*h3(k2)*j2(k2)/(h1(k)*h3(k))
!            j3(k)=j3(k2)
             j1(k)=j1(k2)
             j2(k)=j2(k2)
             j3(k)=j3(k2)
          enddo
! ionospheric boundary conditions
!
 100      continue
          j=2
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j+1-1)*n2
             if(CENTRED)k2=i+(j+2-1)*n2
             j1(k)=-j1(k2)
             if(NEWBC)j1(k)=j1(k2)*h1(k2)/h1(k)
             j2(k)=-j2(k2)
             j3(k)=-j3(k2)
          enddo

          j=1
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j+3-1)*n2
             if(CENTRED)k2=i+(j+4-1)*n2
             j1(k)=-j1(k2)
             if(NEWBC)j1(k)=j1(k2)*h1(k2)/h1(k)
             j2(k)=-j2(k2)
             j3(k)=-j3(k2)
          enddo


          j=n1-1
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j-1-1)*n2
             if(CENTRED)k2=i+(j-2-1)*n2
             j1(k)=-j1(k2)
             if(NEWBC)j1(k)=j1(k2)*h1(k2)/h1(k)
             j2(k)=-j2(k2)
          enddo

          j=n1
          do i=1,n2
             k=i+(j-1)*n2
             k2=i+(j-3-1)*n2
             if(CENTRED)k2=i+(j-4-1)*n2
             j1(k)=-j1(k2)
             if(NEWBC)j1(k)=j1(k2)*h1(k2)/h1(k)
             j2(k)=-j2(k2)
          enddo

          return
          END SUBROUTINE current

