	  SUBROUTINE particle_update(size2,rank,root)
          use indeces
          use param
          use parparam
          implicit none
          integer :: k,ii,jj,size2,rank,root,offset
          do ii=1,np2
              mu0_p(ii)=mu1_p(ii)
              mu1_p(ii)=mu2_p(ii)
              vpar0(ii)=vpar1(ii)
              vpar1(ii)=vpar2(ii)
          enddo
          return
          END SUBROUTINE particle_update
