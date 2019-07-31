       MODULE fields
       use indeces
!      double precision :: y(n),y0(n),y2(n),f_1(n),f_2(n),y_sum(n)
       double precision, dimension(:), allocatable :: y,y0,y2,f_1,f_2,y_sum
!      double precision :: e2(n21),e2_sum(n21),e3(n21)
       double precision, dimension(:), allocatable :: e2,e2_sum,e3
!      double precision :: j1_0(n21),j1_1(n21),j1_2(n21)
       double precision, dimension(:), allocatable :: j1_0,j1_1,j1_2
!      double precision :: e1_0(n21),e1_1(n21),e1_2(n21)
       double precision, dimension(:), allocatable :: e1_0,e1_1,e1_2
!      double precision :: e1_t(n21),e1_sum(n21)
       double precision, dimension(:), allocatable :: e1_t,e1_sum
!      double precision :: j2(n21),j3(n21),divb(n21)
       double precision, dimension(:), allocatable :: j2,j3,divb
!      double precision :: ec0(n21),ec1(n21),divj_p(n21),divj_c(n21)
       double precision, dimension(:), allocatable :: ec0,ec1,divj_p,divj_c
!      double precision :: Fc0(n21),Fc1(n21)
       double precision, dimension(:), allocatable :: Fc0,Fc1
!      double precision :: j1x1(n21),j2x2(n21)
       double precision, dimension(:), allocatable :: j1x1,j2x2,j3x3
       END MODULE fields
