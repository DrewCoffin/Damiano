        SUBROUTINE p_initialize(size2,rank,root)
	use indeces
        use param
        use mesh
        use parparam
        use initial_arrays

        implicit none
        integer :: num,i,j,nsum1(1000),nsum2(1000),ip
        integer :: icount,ii,jj,nn,rank,root,size2,n2p2

        double precision :: v2,nmax
        double precision :: dv,temp1,temp2
        double precision :: drg0,sum
        double precision :: drgz0,D,rg0temp(n2p)
        double precision :: dmu0,dx2_temp
        double precision :: vb(nbin),rb(1000),drb,total,zb(1000),dzb

        double precision :: randf,temp
        double precision :: dp,mu_temp(n1p),r_temp,theta_temp,dth,mu,nu_t
        double precision :: vbar,theta_t1,theta_t2
	double precision :: Deltax1,Deltax2,dx1p,dx2p

	double precision :: px1(n1p/2),dpx1,A,rho_mu_nu
        double precision :: h_mu,h_nu,h_phi,r_t,th_t,h_mu2
        double precision :: x2max,x2min,dx2g
        external randf,dp,rho_mu_nu,mu,h_mu,h_nu,h_phi,h_mu2
        allocate(pint_temp(n_int))
        allocate(x1_temp(n_int))
        allocate(F(n_int))
        seed=10001

        n2p2=n2p/nproc
        np2=n2p2*n1p
        x2min=x2p1
        x2max=x2pn2
	Deltax2=(x2pn2-x2p1)
	Deltax1=(x1pn1-x1p1)
! pick probabiliities
!
        dpx1=1./float(n1p/2-1)
        do ii=2,n1p/2
           px1(ii)=px1(ii-1)+dpx1
        enddo
        dx1p=x1pn1/float(n_int)

! uniform density in x1
!
! ip is particle index here
! ii,jj index for n2p2,n1p respectively
        ip=1
        dx2_temp=(x2pn2-x2p1)/float(n2p-1)
        DO ii=1,n2p2
          nu_p(ip)=x2p1+(rank+float(ii-1)*size2)*dx2_temp
          x1_temp(1)=0.d0
          F(1)=0.d0
          do i=2,n_int
             x1_temp(i)=x1_temp(i-1)+dx1p
             r_t=rho_mu_nu(x1_temp(i),nu_p(ip))
             th_t=acos(r_t**2*x1_temp(i))
!            F(i)=h_mu(r_t,th_t)*h_nu(r_t,th_t)*h_phi(r_t,th_t)
             F(i)=dp(x1_temp(i))* &
!            F(i)=dp2(th_t)* &
                  h_mu(r_t,th_t)*h_nu(r_t,th_t)*h_phi(r_t,th_t)
 
          enddo
          do i=3,n_int
             pint_temp(i)=pint_temp(i-1)+dx1p*F(i)
          enddo
          nn=3
! loop goes from 1+np1/2 to n1p (in terms of mu_temp)
!
          do 20 i=1,n1p/2-1
             temp=pint_temp(n_int)*px1(i)
             do 10 j=nn,n_int
!                if(i.eq.n1p/2)write(*,55) pint_temp(j),temp,x1_temp(j),
!     &                        mu_temp(i+n1p/2-1),mu_temp(i+n1p/2-2)
                if(pint_temp(j).lt.temp) goto 10
!                if(i.eq.n1p/2)write(*,55) 1,pint_temp(j),x1_temp(j)
                nn=j-1
                mu_temp(i+n1p/2)=x1_temp(j)
                goto 20
 10          continue
 20       continue
          mu_temp(n1p)=x1pn1
 55       format(5e13.6)
          do i=1,n1p/2-1
             mu_temp(i)=-mu_temp(n1p-i+1)
          enddo
          do jj=1,n1p
                mu0_p(ip)=mu_temp(jj)
                mu1_p(ip)=mu0_p(ip)
                nu_p(ip)=x2p1+(rank+float(ii-1)*size2)*dx2_temp
                ro_p(ip)=1./nu_p(ip)
                ip=ip+1
           enddo
          ENDDO
          do i=1,np2
             call rootfind(mu0_p(i),nu_p(i),r_p(i), &
                     theta_p(i))
             h1p_1(i)=h_mu2(r_p(i),cos(theta_p(i)))
          enddo
! assign initial values to ltheta
!
        do i=1,np2
          do j=1,n1-1
             if((mu1_p(i).ge.x1(j)).and. &
                 (mu1_p(i).lt.x1(j+1)))then
                 ltheta(i)=j
                 if((ltheta(i).lt.0).or.(ltheta(i).gt.n1))then
                    write(*,*) 'ltheta out of bounds'
                    stop
                 endif
             endif
          enddo
        enddo

! parallel velocity distribution using Box-Muller
!

        vmin=-8.0
        vmax=8.0
        if(vmin.eq.0.) then
           vmin=-3.0
           vmax=3.0
        endif
        dv=(vmax-vmin)/nbin
        vb(1)=vmin
!*$*ASSERT DO(SERIAL)
        do i=2,nbin+1
           vb(i)=vmin+dv*float(i-1)
        enddo

        deallocate(pint_temp)
        deallocate(x1_temp)
        deallocate(F)

 999    return
        END SUBROUTINE p_initialize

        SUBROUTINE vrangefind(np,v,vmin,vmax)
        implicit none
        integer :: i,np
        double precision :: v(*),vmin,vmax
        vmin=0.d0
        vmax=0.d0
        do i=1,np
           vmax=max(vmax,v(i))
           vmin=min(vmin,v(i))
        enddo
        return
        END SUBROUTINE

        FUNCTION dp(x)
        double precision :: x,dp,no
!       if(x.gt.0)then
!          no=1.0
!       else
!          no=0.1
!       endif
        no=0.1
!       dp=(no+(0.1-no)*exp(-x**2/(0.015**2)))
        dp=(no+(4.0-no)*exp(-x**2/(0.006**2)))
!       dp=(no+(0.1-no)*exp(-x**2/(0.015**2)))
!       dp=exp(-x**2/(0.12**2))
        return
        END FUNCTION dp

