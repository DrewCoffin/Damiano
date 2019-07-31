	SUBROUTINE corr(dt)
	use indeces
	use fields
	implicit none
	integer :: ii
	double precision :: dt
        do ii=1,n
           y(ii)=y0(ii)+(f_1(ii)+f_2(ii))/2.*dt
	enddo
	return
	END SUBROUTINE corr

