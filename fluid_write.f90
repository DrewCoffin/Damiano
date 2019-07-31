	SUBROUTINE fluid_write(index_t)
	use indeces
	use mesh
	use fields
        use parparam
	implicit none
        character(len=7) :: fname
        integer :: index_t

        call getname('d_fl',index_t,fname)

        open(10,file=fname,form='unformatted',status='unknown')
        write(10) y
        write(10) j1delta1,j1delta2
        write(10) e1_2,e2,e3,j1_2,j2,j3,divj_c,j1x1,j2x2,divb,ec1, &
                  e1_sum,e2_sum,y_sum
        write(10) den_grid0,den_grid,j1p_grid,mu_grid,mu_grid2,sm_grid
        write(10) rho,va,le2
        close(10)
        return
        END SUBROUTINE fluid_write

