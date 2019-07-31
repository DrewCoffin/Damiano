        SUBROUTINE filterd(n2,n1,sn,sn2)
        implicit none
        integer :: i,j,k,n1,n2
        double precision :: sn(*),sn2(*)
        double precision :: sn_t1(n2,n1),sn_t2(n2,n1)

! 1/4-1/2-1/4 digital filter as per Birdsall and Langdon Appendix C

        do j=1,n1
           do i=1,n2
              k=i+(j-1)*n2
              sn_t1(i,j)=sn(k)
           enddo
        enddo
       
! filter in x2 direction
        do j=1,n1
           do i=2,n2-1
              sn_t2(i,j)=0.25*sn_t1(i-1,j)+0.5*sn_t1(i,j) &
                           +0.25*sn_t1(i+1,j)
           enddo
        enddo          
        
! filter in x1 direction
        do j=3,n1-2
           do i=2,n2-1
              k=i+(j-1)*n2
              sn2(k)=0.25*sn_t2(i,j+1)+0.5*sn_t2(i,j) &
                       +0.25*sn_t2(i,j-1)
!              if(j.eq.3)
!     &          sn_t2(i,j)=0.5*sn_t1(i,j)+0.5*sn_t1(i,j+1)
!              if(j.eq.n2-2)
!     &          sn_t2(i,j)=0.5*sn_t1(i,j)+0.5*sn_t1(i,j-1)
           enddo
        enddo


        return
        END SUBROUTINE filterd
 


