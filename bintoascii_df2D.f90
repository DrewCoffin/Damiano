         PROGRAM bintoascii_df2D
!
! program to construct 2D distribution functions from the data in the 
! d_paxxx_yy files. Max/Min vperp and vpar velocities must be specified
! as does the number of sample bins in each direction and the range in 
! x2 and x1 from which the particle data must be taken.
!
         implicit none
         integer :: i,ii,jj,istep,imin,imax,index,k
         integer :: ichoice,jchoice
         integer, parameter :: n1=512,n2=144,n21=n2*n1,n=n21*5,nbini=30, &
                               nbinj=30
         integer, parameter :: npt=144*2000000,nproc=144,np=npt/nproc
         integer, parameter :: npointsi=nbini+1,npointsj=nbinj+1
         integer :: iprecip(np),iparticle,ibin,jbin,bin(npointsi,npointsj)

! npoints=(vmax-vmin)/dv+1
! nbin=(vmax-vmin)/dv, dv=(vmax-vmin)/float(nbin)
!
         double precision :: y(n)
         double precision :: x(n2,n1),z(n2,n1)
         double precision :: x1(n1),x2(n2),dummy1,dummy2
         double precision :: rho(n21),va(n21),le2(n21)
         double precision :: Re,L,me,mp,c,e,Bn,rho_n,vn,en,pi
	 double precision :: mu_p(np),nu_p(np),ro_p(np),vpar(np),mu_m(np)
         double precision :: r_p,theta_p,h1_p,Bo_p,M,h_mu,vperp2,T_e_par
         double precision :: vperp_min,vperp_max,vpar_min,vpar_max
         double precision :: dvperp,dvpar,ro_d_min,ro_d_max
         double precision :: vbin_perp(npointsi),vbin_par(npointsj)
         double precision :: mu_o
         external h_mu

         character(len=7) stemp,fname
         character(len=11) fname2,rtemp
         vperp_min=-70.d0
         vperp_max=70.d0
         vpar_min=-70.d0
         vpar_max=70.d0

         dvperp=(vperp_max-vperp_min)/float(nbini)
         dvpar=(vpar_max-vpar_min)/float(nbinj)

         do ii=1,npointsi
            vbin_perp(ii)=vperp_min+dvperp*float(ii-1)
         enddo

         do jj=1,npointsj
            vbin_par(jj)=vpar_min+dvpar*float(jj-1)
         enddo

 
! zero bin array
!
        do ii=1,nbini
           do jj=1,nbinj
              bin(ii,jj)=0
           enddo
        enddo
! input indeces
!
         write(*,*) n1,n2,n21,n
         write(*,*) 'imin?'
         read(*,*) imin
         write(*,*) 'imax?'
         read(*,*) imax
         write(*,*) 'istep?'
         read(*,*) istep
         pi=4.0*atan(1.0)
         M=8.02e22
         Re=6.38e6
         L=Re
         me=9.109d-31
         mp=1.67d-27
         c=2.998d8
         e=1.609d-19
         Bn=10.d-9
         rho_n=1.d5*mp
         mu_o=4.d-7*pi 
         vn=Bn/sqrt(mu_o*rho_n)
         en=Bn*vn

! read in x1 from file
!
         open(2,file='x1.dat',status='old')
         do ii=1,n1
            read(2,*) x1(ii),dummy1,dummy2
         enddo
         close(2)
! read in x2 from file
!
         open(2,file='x2.dat',status='old')
         do ii=1,n2
            read(2,*) dummy1,x2(ii)
         enddo
         close(2)

! loop over index
!
         if(imin.eq.imax)then
           open(unit=12,file='electron.map',status='unknown')
         endif 
         do 444 index=imin,imax,istep 
            do 445 i=1,nproc
                call getname_mpi('d_pa',index,i-1,fname2)
                write(*,*) 'reading file ',fname2
                open(11,file=fname2,form='unformatted',status='old')
                read(11) mu_p,nu_p,vpar,mu_m,iprecip
                close(11)      
                open(unit=12,file='electron.map',status='unknown')
                jchoice=139
!                jchoice=149
!                jchoice=159
                ichoice=64
                ro_d_max=1./x2(ichoice-2)
                ro_d_min=1./x2(ichoice+2)
                do 446 ii=1,np
                   ro_p(ii)=1./nu_p(ii)
                   if(ro_p(ii).gt.ro_d_min.and. &
!                    ro_p(ii).lt.ro_d_max.and.mu_p(ii).gt.x1(jchoice-2). &
!                      and.mu_p(ii).le.x1(jchoice+2))then
                    ro_p(ii).lt.ro_d_max.and.mu_p(ii).gt.x1(jchoice-2) &
                    .and.mu_p(ii).le.x1(jchoice+2).and.vpar(ii).lt.5. &
                    .and.vpar(ii).gt.-5)then 
                      call rootfind(pi,mu_p(ii),nu_p(ii),r_p, &
                                theta_p)
                      h1_p=h_mu(r_p,theta_p)
                      Bo_p=mu_o*M/(4.d0*pi*h1_p*Re**3)/Bn
                      vperp2=2.d0*Bo_p*mu_m(ii)
                      ibin=int((sqrt(vperp2)-vperp_min)/dvperp)+1
                      jbin=int((vpar(ii)-vpar_min)/dvpar)+1
                      if(ibin.lt.1.or.ibin.gt.nbini)goto 50
                      if(jbin.lt.1.or.jbin.gt.nbinj)goto 50
                      bin(ibin,jbin)=bin(ibin,jbin)+1
!                      if(imin.eq.imax.and.vpar(ii).gt.0.d0)
!     &                write(12,*) i-1,ii,vpar(ii)
                      if(imin.eq.imax) &
                      write(12,*) i-1,ii,vpar(ii)
 50                   continue
                   endif
 446            continue
 445        continue
            call getname('d_vv',index,fname)
            open(10,file=fname,status='unknown')
            do ii=nbini/2+1,nbini
               do jj=1,nbinj
                  write(10,*) 0.5*(vbin_perp(ii)+vbin_perp(ii+1)), &
                              0.5*(vbin_par(jj)+vbin_par(jj+1)), &
                              bin(ii,jj)/ &
                           (0.5*(vbin_perp(ii)+vbin_perp(ii+1)))   
                  bin(ii,jj)=0
               enddo
!               write(10,*) 
            enddo
            close(10)
 444     continue
         if(imin.eq.imax)close(12)
 55      format(6e13.5)
         stop
         END PROGRAM bintoascii_df2D


          SUBROUTINE getname(temp,index,fname)
          implicit none
          integer one, ten, hundred,index
          character*4 temp
          character*7 fname
          hundred=index/100
          ten=(index-100*hundred)/10
          one=index - ten*10 - hundred*100
          fname(1:4)=temp(1:4)
          fname(5:5)=char(hundred+48)
          fname(6:6)=char(ten+48)
          fname(7:7)=char(one+48)
 10       format(a5)
          return
          END SUBROUTINE

          SUBROUTINE getname_mpi(temp,index,rank,fname)
          implicit none
          integer one, ten, hundred,index
          integer rank
          character*4 temp
          character*11 fname
          hundred=index/100
          ten=(index-100*hundred)/10
          one=index - ten*10 - hundred*100
          fname(1:4)=temp(1:4)
          fname(5:5)=char(hundred+48)
          fname(6:6)=char(ten+48)
          fname(7:7)=char(one+48)
          fname(8:8)='_'
          hundred=rank/100
          ten=(rank-100*hundred)/10
          one=rank - ten*10 - hundred*100
          fname(9:9)=char(hundred+48)
          fname(10:10)=char(ten+48)
          fname(11:11)=char(one+48)
 10       format(a5)
          return
          END SUBROUTINE

        SUBROUTINE rootfind(pi,mu,nu,r,theta)
        implicit none
        double precision a,b,r,mu,nu
        double precision theta,c1,c2,pi
        complex*16 cterm1,cterm2,cterm3
        complex*16 x(4),temp1,temp2,temp3
        if(mu.eq.0d0)then
          r=1.d0/nu
          theta=pi/2.d0
          goto 100
        endif
        c1=1.d0/3.d0
        a=mu*mu
        b=nu
        cterm1=4.d0*2.d0**c1
        cterm2=(27.d0*a*b*b+ &
               sqrt(6912.d0*a*a*a+729.d0*a*a*b*b*b*b))**c1
        cterm3=(3.d0*a*2**c1)
        temp1=-cterm1/cterm2+cterm2/cterm3
        c2=sqrt(temp1)
!        temp2=cterm1/cterm2-cterm2/cterm3-
!     &     (2.d0*b)/(a*c2)
        temp3=cterm1/cterm2-cterm2/cterm3+ &
           (2.d0*b)/(a*c2)
        x(4)=-0.5d0*c2+0.5d0*sqrt(temp3)
        r=dreal(x(4))
        theta=acos(r*r*mu)
 100    continue
        return
        END SUBROUTINE

         FUNCTION h_mu(r,theta)
         implicit none
         double precision r,theta,h_mu,c2
         c2=cos(theta)
         h_mu=r*r*r/sqrt((1+3.*c2*c2))
         return
         END FUNCTION

        SUBROUTINE fdistrb(index_df,np,ro_d_min,ro_d_max,vmin, &
        vmax,vpar2,ro_p)
        implicit none
        integer ii,jj,ic(500),np,index_df,nbin
        double precision dv,vbin(500),vpar2(*),ro_p(*)
        double precision vth,ro_d_min,ro_d_max,vmin,vmax
        character*7 fname
        nbin=100
        call getname('d_tf',index_df,fname)
        dv=(vmax-vmin)/float(nbin)
        do ii=1,nbin
           vbin(ii)=0.
           ic(ii)=0
        enddo
        vbin(1)=vmin
        do ii=2,nbin+1
            vbin(ii)=vmin+dv*float(ii-1)
        enddo

        do ii=1,np
           if(ro_p(ii).gt.ro_d_min.and. &
              ro_p(ii).lt.ro_d_max)then
           do jj=1,nbin
               if((vpar2(ii).gt.vbin(jj)).and. &
                  (vpar2(ii).le.vbin(jj+1)))then
                  ic(jj)=ic(jj)+1
               endif
           enddo
           endif
        enddo
!
        open(unit=9,file=fname,status='unknown')
        do  ii=1,nbin
           write(9,*) vbin(ii),ic(ii)
           write(9,*) vbin(ii+1),ic(ii)
        enddo
        close(9)

        return
        END SUBROUTINE


