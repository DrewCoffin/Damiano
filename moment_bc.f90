	SUBROUTINE moment_bc()
	use indeces
        use param
        use mesh
        use fields
        use parparam
        use mpi
	implicit none
	integer :: ii,jj,j,k1,k2,k3,k4,k5
        do ii=1,n2
           do jj=1,n1
              k1=ii+(jj-1)*n2
              den_grid(k1)=ratio_grid(k1)*den_grid(k1)/ &
                          h3(k1)
              j1p_grid(k1)=-const_jp*ratio_grid(k1)*j1p_grid(k1)/ &
                          h3(k1)
              sm_grid(k1)=const_sm*ratio_grid(k1)*sm_grid(k1)/ &
                          h3(k1)
              mu_grid(k1)=const_mu*ratio_grid(k1)*mu_grid(k1)/ &
                          h3(k1)
              mu_grid2(k1)=const_mu*ratio_grid(k1)*mu_grid2(k1)/ &
                          h3(k1)
           enddo
        enddo

        do ii=1,n2
            IF(CENTRED)THEN
               k1=ii+n2*(1-1)
               k2=ii+n2*(2-1)
               k3=ii+n2*(3-1)
               k4=ii+n2*(4-1)
               k5=ii+n2*(5-1)
            ELSE
               k1=ii+n2*(1-1)
               k2=ii+n2*(2-1)
               k4=ii+n2*(3-1)
               k5=ii+n2*(4-1)
            ENDIF
            j1p_grid(k1)=j1p_grid(k5)
            j1p_grid(k2)=j1p_grid(k4)
            IF(NEWBC)THEN
               j1p_grid(k1)=j1p_grid(k5)*h1(k5)/h1(k1)
               j1p_grid(k2)=j1p_grid(k4)*h1(k4)/h1(k2)
            ENDIF
            mu_grid(k1)=mu_grid(k5)
            mu_grid(k2)=mu_grid(k4)
            IF(CENTRED)THEN
               k1=ii+n2*(n1-1)
               k2=ii+n2*(n1-1-1)
               k3=ii+n2*(n1-2-1)
               k4=ii+n2*(n1-3-1)
               k5=ii+n2*(n1-4-1)
            ELSE
               k1=ii+n2*(n1-1)
               k2=ii+n2*(n1-1-1)
               k4=ii+n2*(n1-2-1)
               k5=ii+n2*(n1-3-1)
            ENDIF
            j1p_grid(k1)=j1p_grid(k5)
            j1p_grid(k2)=j1p_grid(k4)
            IF(NEWBC)THEN
               j1p_grid(k1)=j1p_grid(k5)*h1(k5)/h1(k1)
               j1p_grid(k2)=j1p_grid(k4)*h1(k4)/h1(k2)
            ENDIF
            mu_grid(k1)=mu_grid(k5)
            mu_grid(k2)=mu_grid(k4)
        enddo 
        do j=1,n1
           den_grid(1+(j-1)*n2)=den_grid(3+(j-1)*n2)
           den_grid(n2+(j-1)*n2)=den_grid(n2-2+(j-1)*n2)
           j1p_grid(1+(j-1)*n2)=-j1p_grid(3+(j-1)*n2)
           j1p_grid(n2+(j-1)*n2)=j1p_grid(n2-2+(j-1)*n2)
           sm_grid(1+(j-1)*n2)=sm_grid(3+(j-1)*n2)
           sm_grid(n2+(j-1)*n2)=sm_grid(n2-2+(j-1)*n2)
           mu_grid(1+(j-1)*n2)=mu_grid(3+(j-1)*n2)
           mu_grid(n2+(j-1)*n2)=mu_grid(n2-2+(j-1)*n2)
        enddo
	
        if(ifilter.eq.1)then
!           call filterd_r(n2,n1,j1p_grid,j1p_grid)
!           call filterd_r(n2,n1,sm_grid,sm_grid)
!           call filterd_r(n2,n1,mu_grid,mu_grid)
!           call filter_p(n2,n1,j1p_grid,j1p_grid)
        endif
	return
        END SUBROUTINE moment_bc
