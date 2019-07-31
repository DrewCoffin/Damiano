         SUBROUTINE scaling()
         use indeces
         use param
         use fields
         use mesh
         use parparam

         implicit none
         integer ii,jj,kk

! set up ratio_grid
!
         open(unit=2,file='ratio.dat',status='unknown')
         do ii=2,n2-1
            do jj=2,n1-1
               kk=ii+(jj-1)*n2
               ratio_grid(kk)=rho(kk)/den_grid0(kk)
               write(2,55) x(ii,jj),z(ii,jj),ratio_grid(kk)
            enddo
            write(2,*)
         enddo
         close(2)
 55      format(3e13.5)
!
         call moment_bc()
         return
         END SUBROUTINE scaling

