	SUBROUTINE pred(dt)
	use indeces
	use fields
	implicit none
	integer :: ii
	double precision :: dt
        do 60 ii=1,n
              y2(ii)=y0(ii)+2.*f_1(ii)*dt
              y0(ii)=y(ii)
 60     continue
	return
	END SUBROUTINE pred
