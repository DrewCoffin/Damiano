        SUBROUTINE pcurrent0(size2,rank,root)
! subroutine to determine initial particle density profile. 
        use indeces
        use param
        use mesh
        use fields
        use parparam
        use mpi
        implicit none
        integer :: i,ii,jj,k,j,lro
        integer :: k1,k2,k3,k4
        integer :: size2,rank,root,offset,error
        double precision :: h_mu,length,h_nu,l1,l2,w1,w2,atotal
        double precision :: c1,dro_i
        double precision :: a1,a2,a3,a4,asum_sq_i
        double precision :: den_b(n21),j1p_b(n21),mu_b(n21),sm_b(n21)
        double precision :: mu_b2(n21)
        external h_mu,length,h_nu

!        np2=np/size2
        k=rank+1
        offset=np2*(k-1)

        dro_i=1./dro

        do ii=1,n21
           den_b(ii)=0.d0
           den_grid0(ii)=0.d0
           j1p_b(ii)=0.d0
           mu_b(ii)=0.d0
           mu_b2(ii)=0.d0
           sm_b(ii)=0.d0
           den_grid(ii)=0.d0
           den_temp(ii)=0.d0
           j1p_grid(ii)=0.d0
           mu_grid(ii)=0.d0
           mu_grid2(ii)=0.d0
           sm_grid(ii)=0.d0
           j1p_temp(ii)=0.d0
        enddo


        do ii=1,np2
          do j=1,n1-1
             if((mu1_p(ii).ge.x1(j)).and. &
                 (mu1_p(ii).lt.x1(j+1)))then
                 ltheta(ii)=j
                 if((ltheta(ii).lt.0).or.(ltheta(ii).gt.n1))then
                    write(*,*) 'ltheta out of bounds'
                    stop
                 endif
             endif
          enddo
        enddo

        do ii=1,np2
          lro=(ro_p(ii)-ro1)*dro_i+1
          k1=lro+(ltheta(ii)-1)*n2
!          l2=abs(length(ro_p(ii),pi/2-theta(lro,ltheta(ii)+1)
!     &              ,pi/2-theta_p(ii)))
!          l1=abs(length(ro_p(ii),pi/2.-theta_p(ii),
!     &              pi/2.-theta(lro,ltheta(ii))))
          l2=abs(0.5*(h1(k1+n2)+h1p_1(ii))* &
             (mu0_p(ii)-x1(ltheta(ii)+1)))
          l1=abs(0.5*(h1(k1)+h1p_1(ii))* &
             (mu0_p(ii)-x1(ltheta(ii))))
          w1=abs((h2(k1)+h_nu(r_p(ii),theta_p(ii)))/2. &
            *(x2(lro)-nu_p(ii)))
          w2=abs((h2(k1+1)+h_nu(r_p(ii),theta_p(ii)))/2. &
            *(x2(lro+1)-nu_p(ii)))

          a1=l1*w1
          a2=l1*w2
          a3=l2*w1
          a4=l2*w2
          asum_sq_i=(1.d0/(a1+a2+a3+a4))**2
          k1=lro+1+(ltheta(ii)+1-1)*n2
          k2=lro+(ltheta(ii)+1-1)*n2
          k3=lro+1+(ltheta(ii)-1)*n2
          k4=lro+(ltheta(ii)-1)*n2
          den_b(k1)=den_b(k1)+ &
                       a1*asum_sq_i
          den_b(k2)=den_b(k2)+ &
                       a2*asum_sq_i
          den_b(k3)=den_b(k3)+ &
                       a3*asum_sq_i
          den_b(k4)=den_b(k4)+ &
                       a4*asum_sq_i
          j1p_b(k1)=j1p_b(k1)+ &
                vpar0(ii)*a1*asum_sq_i
          j1p_b(k2)=j1p_b(k2)+ &
                vpar0(ii)*a2*asum_sq_i
          j1p_b(k3)=j1p_b(k3)+ &
                vpar0(ii)*a3*asum_sq_i
          j1p_b(k4)=j1p_b(k4)+ &
                vpar0(ii)*a4*asum_sq_i
! magnetic moment moment
          mu_b(k1)=mu_b(k1)+ &
                mu_m(ii)*a1*asum_sq_i
          mu_b(k2)=mu_b(k2)+ &
                mu_m(ii)*a2*asum_sq_i
          mu_b(k3)=mu_b(k3)+ &
                mu_m(ii)*a3*asum_sq_i
          mu_b(k4)=mu_b(k4)+ &
                mu_m(ii)*a4*asum_sq_i
! vpar magnetic moment
!
! Bn is incorporated here because Bo is not normalized until the subroutine 
! fl_normalize

          mu_b2(k1)=mu_b2(k1)+ &
          0.5*vpar0(ii)*vpar0(ii)/(Bo(k1)/Bn) &
                *a1*asum_sq_i
          mu_b2(k2)=mu_b2(k2)+ &
          0.5*vpar0(ii)*vpar0(ii)/(Bo(k2)/Bn) &
                *a2*asum_sq_i
          mu_b2(k3)=mu_b2(k3)+ &
          0.5*vpar0(ii)*vpar0(ii)/(Bo(k3)/Bn) &
                *a3*asum_sq_i
          mu_b2(k4)=mu_b2(k4)+ &
          0.5*vpar0(ii)*vpar0(ii)/(Bo(k4)/Bn) &
                *a4*asum_sq_i
! second moment
          sm_b(k1)=sm_b(k1)+ &
             vpar0(ii)*vpar0(ii)*a1*asum_sq_i 
          sm_b(k2)=sm_b(k2)+ &
             vpar0(ii)*vpar0(ii)*a2*asum_sq_i
          sm_b(k3)=sm_b(k3)+ &
             vpar0(ii)*vpar0(ii)*a3*asum_sq_i
          sm_b(k4)=sm_b(k4)+ &
             vpar0(ii)*vpar0(ii)*a4*asum_sq_i
        enddo
      call MPI_REDUCE(den_b,den_grid0,n21,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,root,MPI_COMM_WORLD,error)
      call MPI_REDUCE(j1p_b,j1p_grid,n21,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,root,MPI_COMM_WORLD,error)
      call MPI_REDUCE(mu_b,mu_grid,n21,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,root,MPI_COMM_WORLD,error)
      call MPI_REDUCE(mu_b2,mu_grid2,n21,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,root,MPI_COMM_WORLD,error)
      call MPI_REDUCE(sm_b,sm_grid,n21,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,root,MPI_COMM_WORLD,error)
      call MPI_BCAST(den_grid,n21,MPI_DOUBLE_PRECISION, &
                       root,MPI_COMM_WORLD,error)
        IF(rank.eq.root)then
        do ii=2,n2-1
           do jj=2,n1-1
              k1=ii+(jj-1)*n2
              den_grid(k1)=den_grid0(k1)
              den_grid0(k1)=den_grid0(k1)/h3(k1)
           enddo
        enddo

! radial boundary conditions
!
        do jj=1,n1
           k1=1+(jj-1)*n2
           k2=2+(jj-1)*n2
           den_grid0(k1)=den_grid0(k2)
           k1=n2-1+(jj-1)*n2
           k2=n2+(jj-1)*n2
           den_grid0(k2)=den_grid0(k1)
        enddo
        ENDIF
        return
        END SUBROUTINE pcurrent0
