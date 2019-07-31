       MODULE mesh
       use param
!      double precision :: r(n2,n1),theta(n2,n1)
       double precision, dimension(:,:), allocatable :: r,theta
!      double precision :: x(n2,n1),z(n2,n1)
       double precision, dimension(:,:), allocatable :: x,z
!      double precision :: h2(n21),h1(n21),h3(n21)
       double precision, dimension(:), allocatable :: h2,h1,h3
!      double precision :: Bo(n21),Va(n21),rho(n21),dBo(n21)
       double precision, dimension(:), allocatable :: Bo,Va,rho,dBo
       double precision, dimension(:), allocatable :: gamma2
!      double precision :: x1(n1),x2(n2),le2(n21),rho_i(n21)
       double precision, dimension(:), allocatable :: x1,x2,le2,rho_i
!      double precision :: len(n2,n1),width(n2,n1)
       double precision, dimension(:,:), allocatable :: len,width
!      double precision :: alpha(n2,n1),area(n2,n1)
       double precision, dimension(:,:), allocatable :: alpha,area
!      double precision :: area_ratio(n2,n1),cosalpha(n2,n1)
       double precision, dimension(:,:), allocatable :: area_ratio,cosalpha
       double precision :: j1delta1,j1delta2
       END MODULE mesh
