        SUBROUTINE fl_normalize()
        use indeces
        use param
        use fields
        use mesh
        use parparam
        implicit none
        integer ii,jj,kk
          
        open(unit=2,file='normal.dat',status='unknown')
        do ii=1,n2
            do jj=1,n1
                kk=ii+(jj-1)*n2
                Bo(kk)=Bo(kk)/Bn
                rho(kk)=rho(kk)/rho_n
                Va(kk)=Va(kk)/vn
                le2(kk)=le2(kk)/(L*L)
                write(2,55) x(ii,jj),z(ii,jj),Bo(kk), &
                         rho(kk),Va(kk),sqrt(le2(kk))
             enddo
             write(2,*)
          enddo
 55       format(6e13.5)
          close(2)
          return
          END SUBROUTINE fl_normalize
