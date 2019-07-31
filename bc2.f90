         subroutine bc2()
	 use indeces
         use param
         use mesh
         use fields
         implicit none
         integer :: i,j,k,ku2,ku3,kb1,kb2,kb3
         integer :: ku2_2,ku3_2,kb1_2,kb2_2,kb3_2
         double precision :: term1,term2,term3,term4,term5,term6,term7

! boundary conditions
!
          j=2
          do i=1,n2
            ku2=i+n2*(j-1)
            ku3=ku2+n21
            kb1=ku2+n21*2
            kb2=ku2+n21*3
            kb3=ku2+n21*4
            ku2_2=i+n2*(j+1-1)
            if(CENTRED)ku2_2=i+n2*(j+2-1)
            ku3_2=ku2_2+n21
            kb1_2=ku2_2+n21*2
            kb2_2=ku2_2+n21*3
            kb3_2=ku2_2+n21*4
            y(ku2)=-y(ku2_2)
            y(ku3)=-y(ku3_2)
            y(kb1)=-y(kb1_2)
            y(kb2)=y(kb2_2)*h2(ku2_2)/h2(ku2)
            y(kb3)=y(kb3_2)*h3(ku2_2)/h3(ku2)
          enddo


          j=1
          do i=1,n2
            ku2=i+n2*(j-1)
            ku3=ku2+n21
            kb1=ku2+n21*2
            kb2=ku2+n21*3
            kb3=ku2+n21*4
            ku2_2=i+n2*(j+3-1)
            if(CENTRED)ku2_2=i+n2*(j+4-1)
            ku3_2=ku2_2+n21
            kb1_2=ku2_2+n21*2
            kb2_2=ku2_2+n21*3
            kb3_2=ku2_2+n21*4
            y(ku2)=-y(ku2_2)
            y(ku3)=-y(ku3_2)
            y(kb1)=-y(kb1_2)
            y(kb2)=y(kb2_2)*h2(ku2_2)/h2(ku2)
            y(kb3)=y(kb3_2)*h3(ku2_2)/h3(ku2)
          enddo

          j=n1-1
          do i=1,n2
            ku2=i+n2*(j-1)
            ku3=ku2+n21
            kb1=ku2+n21*2
            kb2=ku2+n21*3
            kb3=ku2+n21*4
            ku2_2=i+n2*(j-1-1)
            if(CENTRED)ku2_2=i+n2*(j-2-1)
            ku3_2=ku2_2+n21
            kb1_2=ku2_2+n21*2
            kb2_2=ku2_2+n21*3
            kb3_2=ku2_2+n21*4
            y(ku2)=-y(ku2_2)
            y(ku3)=-y(ku3_2)
            y(kb1)=-y(kb1_2)
            y(kb2)=y(kb2_2)*h2(ku2_2)/h2(ku2)
            y(kb3)=y(kb3_2)*h3(ku2_2)/h3(ku2)
          enddo

          j=n1
          do i=1,n2
            ku2=i+n2*(j-1)
            ku3=ku2+n21
            kb1=ku2+n21*2
            kb2=ku2+n21*3
            kb3=ku2+n21*4
            ku2_2=i+n2*(j-3-1)
            if(CENTRED)ku2_2=i+n2*(j-4-1)
            ku3_2=ku2_2+n21
            kb1_2=ku2_2+n21*2
            kb2_2=ku2_2+n21*3
            kb3_2=ku2_2+n21*4
            y(ku2)=-y(ku2_2)
            y(ku3)=-y(ku3_2)
            y(kb1)=-y(kb1_2)
            y(kb2)=y(kb2_2)*h2(ku2_2)/h2(ku2)
            y(kb3)=y(kb3_2)*h3(ku2_2)/h3(ku2)
          enddo

!         do i=1,n2
!           j=3
!           ku2=i+n2*(j-1)
!           kb1=ku2+n21*2
!           y(kb1)=0.
!           j=n1-2
!           ku2=i+n2*(j-1)
!           kb1=ku2+n21*2
!           y(kb1)=0.
!         enddo

! radial boundary conditions
          i=1
          do j=1,n1
            ku2=i+n2*(j-1)
            ku3=ku2+n21
            kb1=ku2+n21*2
            kb2=ku2+n21*3
            kb3=ku2+n21*4
            ku2_2=i+2+n2*(j-1)
            ku3_2=ku2_2+n21
            kb1_2=ku2_2+n21*2
            kb2_2=ku2_2+n21*3
            kb3_2=ku2_2+n21*4
            y(ku2)=y(ku2_2)
            y(ku3)=y(ku3_2)
            y(kb1)=y(kb1_2)
            y(kb2)=y(kb2_2)
            y(kb3)=y(kb3_2)
          enddo

          i=n2
          do j=1,n1
            ku2=i+n2*(j-1)
            ku3=ku2+n21
            kb1=ku2+n21*2
            kb2=ku2+n21*3
            kb3=ku2+n21*4
            ku2_2=i-2+n2*(j-1)
            ku3_2=ku2_2+n21
            kb1_2=ku2_2+n21*2
            kb2_2=ku2_2+n21*3
            kb3_2=ku2_2+n21*4
            y(ku2)=y(ku2_2)
            y(ku3)=y(ku3_2)
            y(kb1)=y(kb1_2)
            y(kb2)=y(kb2_2)
            y(kb3)=y(kb3_2)
          enddo

          return
          END SUBROUTINE bc2
