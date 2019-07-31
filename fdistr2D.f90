        SUBROUTINE fdistr2D(size2,rank,root)
        use indeces
        use param
        use mesh
        use fields
        use parparam
        use df
        use mpi
        implicit none
        integer :: ii,jj,i,j,kk,ibin,jbin,vv_bin(nbini,nbinj,n2,n1)
        integer :: lro,k1
        integer ::  size2,rank,root,error,offset,k
        integer :: ntp
        double precision :: vperp_min,vperp_max,vpar_min,vpar_max
        double precision :: vperp2,dvperp,dvpar,dro_i
        double precision :: h1_temp,h_mu,Bo_temp
        external h_mu

        vperp_min=vp_min/vn
        vperp_max=vp_max/vn
        vpar_min=vp_min/vn
        vpar_max=vp_max/vn

        dvperp=(vperp_max-vperp_min)/float(nbini)
        dvpar=(vpar_max-vpar_min)/float(nbinj)

        do ii=1,nbini
          do jj=1,nbinj
             do i=1,n2
               do j=1,n1
                vv_bin(ii,jj,i,j)=0
               enddo
             enddo
          enddo
        enddo
                                                                                
        do ii=1,np2
           ltheta(ii)=0
           do j=1,n1-1
           if((mu2_p(ii).ge.x1(j)).and. &
                 (mu2_p(ii).le.x1(j+1)))then
                 ltheta(ii)=j
           endif
           enddo
        enddo

        dro_i=1.d0/dro
        do 100 ii=1,np2
           if (iprecip(ii).eq.1)goto 100
             lro=(ro_p(ii)-ro1)*dro_i+1
             k1=lro+(ltheta(ii)-1)*n2
             h1_temp=h_mu(r_p(ii),theta_p(ii))
             Bo_temp=mu_o*M/(4.d0*pi*h1_temp*Re**3)/Bn
             vperp2=2.d0*Bo_p(ii)*mu_m(ii)
             ibin=int((sqrt(vperp2)-vperp_min)/dvperp)+1
             jbin=int((vpar2(ii)-vpar_min)/dvpar)+1
             if(ibin.lt.1.or.ibin.gt.nbini)goto 100
             if(jbin.lt.1.or.jbin.gt.nbinj)goto 100
            vv_bin(ibin,jbin,lro,ltheta(ii))= &
                          vv_bin(ibin,jbin,lro,ltheta(ii))+1
 100    continue


!      call MPI_REDUCE(bin_b(1,1,1),fbin(1,1,1),ntp,
!     &     MPI_DOUBLE_PRECISION,
!     &     MPI_SUM,root,MPI_COMM_WORLD,error)
        do 200 ii=1,n2
           do 250 jj=1,n1
           call MPI_REDUCE(vv_bin(1,1,ii,jj),df_array(1,1,ii,jj), &
              nbini*nbinj, &
              MPI_INTEGER, &
              MPI_SUM,root,MPI_COMM_WORLD,error)
 250       continue
 200    continue
        return
        END SUBROUTINE fdistr2D

