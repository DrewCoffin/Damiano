	SUBROUTINE e1_term_write(index_t)
	use indeces
	use mesh
        use e1terms
	implicit none
        character(len=7) :: fname
        integer :: index_t,ii,jj,kk
        double precision :: beta_i,dummy1,dummy2

        call getname('d_tm',index_t,fname)
        open(unit=10,file=fname,status='unknown')
        do ii=1,n2
           do jj=1,n1
              kk=ii+n2*(jj-1)
              beta_i=1./(h2(kk)*h3(kk))
              dummy1=beta_i*(e2_term(ii+1,jj)-e2_term(ii-1,jj))/ &
                              (x2(ii+1)-x2(ii-1))
              dummy2=beta_i*(ec_term(ii+1,jj)-ec_term(ii-1,jj))/ &
                              (x2(ii+1)-x2(ii-1))
              write(10,55) x2(ii),x1(jj),dummy1,dummy2, &
                           mu_term(ii,jj),mu_term2(ii,jj)
           enddo
           write(10,*)
        enddo
 55     format(7e12.4)
        close(10)
        return
        END SUBROUTINE e1_term_write

