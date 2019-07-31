         FUNCTION j1delta_det(j)
         use indeces
         use param
         use fields
         use mesh
         use parparam
         implicit none
         integer :: i,j,k
         double precision :: j1delta_det
         j1delta_det=0.d0 
         do i=n2-19,n2 
            k=i+(j-1)*n2
            j1delta_det=j1delta_det+abs(j1p_grid(k))
         enddo
         j1delta_det=j1delta_det/float(20)
         return
         END FUNCTION j1delta_det
