         PROGRAM bintoascii_energy
! 
! Program to determine the components of energy for the dipolar model as well
! as the total. It computes energy density based on the pressure and mu moments 
! for vpar and vperp respectively. Total electron energy computed from these 
! moments (T_e0) agrees with that calculated from summing up the inividual 
! particles. All physical constants are in MKS units.
!
! d_aa - parallel energy density
! d_bb - perpendicular energy density
! re-investigate equality of const_sm and const_mu
!
! Peter Damiano, last modified Nov 26, 2007.
!
         implicit none
         integer :: i,ii,jj,istep,imin,imax,ichoice,index,k
         integer, parameter :: n1=256,n2=128,n21=n2*n1,n=n21*5
         integer, parameter :: npt=2**7*1000000,nproc=128,np=npt/nproc
         integer :: iprecip(np),iparticle,ix2
         double precision :: y(n),y_sum(n),dx2
         double precision :: e1_2(n21),e2(n21),e3(n21),e1_sum(n21)
         double precision :: e2_sum(n21)
         double precision :: j1_2(n21),j2(n21),j3(n21),divj_c(n21)
         double precision :: divb(n21),ec1(n21)
         double precision :: den_grid0(n21),den_grid(n21),j1p_grid(n21)
         double precision :: sm_grid(n21),mu_grid(n21),mu_grid2(n21)
         double precision :: j1x1(n21),j2x2(n21)
         double precision :: j1p_temp(n21),x(n2,n1),z(n2,n1)
         double precision :: h1(n21),h2(n21),h3(n21),Bo(n21)
         double precision :: x1(n1),x2(n2),dummy1,dummy2
         double precision :: rho(n21),va(n21),le2(n21)
         double precision :: elec_ed(n21),elec_ed2(n21)
         double precision :: T_u2,T_u3,T_b2,T_b3,T_e,T_e2,T_e3,TE,TE_i
         double precision :: Nprim,T_e0
         double precision :: T_eP,dummy,T_u3b30
         double precision :: a,b
         double precision :: Re,L,me,mp,c,e,Bn,rho_n,vn,en,const_sm,pi
         double precision :: const_mu,mu_o
         double precision :: mu_p(np),nu_p(np),ro_p(np)
! use double precision if -D SINGLE not used in compliation
!        double precision :: vpar(np),mu_m(np)
         real :: vpar(np),mu_m(np)
         double precision :: r_p,theta_p,h1_p,Bo_p,M,h_mu,vperp2,T_e_par
         double precision :: j1delta1,j1delta2
         integer :: index_temp(14),lll
         external :: h_mu

         character(len=7) :: stemp,fname
         character(len=11) :: fname2,rtemp


         write(*,*) n1,n2,n21,n
         write(*,*) 'imin?'
         read(*,*) imin
         write(*,*) 'imax?'
         read(*,*) imax
         write(*,*) 'istep?'
         read(*,*) istep
         ichoice=3
         iparticle=1
         pi=4.0*atan(1.0)
         M=8.02d22
         Re=6.38d6
         L=Re
         me=9.109d-31
         mp=1.67d-27
         c=2.998d8
         e=1.609d-19
         Bn=10.d-9
         rho_n=1.d5*mp
         en=Bn*vn
         mu_o=4.d-7*pi
         vn=Bn/sqrt(mu_o*rho_n)
!         const_sm=L*L/en*(4.d0*pi)/(c*c)*(e/L)*vn*vn*rho_n/mp
         const_mu=e*mu_o*L*vn*rho_n/(Bn*mp)
         const_sm=const_mu
         TE_i=0.d0 

         open(2,file='x1.dat',status='old')
         do ii=1,n1
            read(2,*) x1(ii),dummy1,dummy2
         enddo
         close(2)
         open(2,file='x2.dat',status='old')
         do ii=1,n2
            read(2,*) dummy1,x2(ii)
         enddo
         close(2)

         open(unit=10,file='h1.dat',status='old')
         open(unit=11,file='h2.dat',status='old')
         open(unit=12,file='h3.dat',status='old')
         open(unit=13,file='Bo.dat',status='old')
         do ii=1,n2
            do jj=1,n1
               k=ii+(jj-1)*n2
               read(10,*) x(ii,jj),z(ii,jj),h1(k)
               read(11,*) x(ii,jj),z(ii,jj),h2(k)
               read(12,*) x(ii,jj),z(ii,jj),h3(k)
               read(13,*) x(ii,jj),z(ii,jj),Bo(k)
               Bo(k)=Bo(k)/Bn 
            enddo
            read(10,12) stemp
            read(11,12) stemp
            read(12,12) stemp
            read(13,12) stemp
         enddo
         close(10)
         close(11)
         close(12)
         close(13)
 12      format(a7)
         open(unit=2,file='energy.dat',status='unknown')
!         do i=1,10000
!            read(2,*,end=13) rtemp
!         enddo
! 13      continue
         do 444 index=imin,imax,istep 
          write(*,*) index
!          do 444 lll=1,2,1
          index_temp(1)=1
          index_temp(2)=2
          index_temp(3)=3
          index_temp(4)=4
          index_temp(5)=45
          index_temp(6)=30
          index_temp(7)=35
          index_temp(8)=40
          index_temp(9)=46
          index_temp(10)=50
          index_temp(11)=55
          index_temp(12)=60
          index_temp(13)=65
          index_temp(14)=70
!          index=index_temp(lll)
           call getname('d_fl',index,fname)
!           if(index.gt.0)then
           open(10,file=fname,form='unformatted',status='old')
           read(10) y
           read(10) j1delta1,j1delta2
           read(10) e1_2,e2,e3,j1_2,j2,j3,divj_c,j1x1,j2x2,divb,ec1, &
                    e1_sum,e2_sum,y_sum
           read(10) den_grid0,den_grid,j1p_grid,mu_grid,mu_grid2,sm_grid
           read(10) rho,va,le2
           close(10)
           if(index.eq.1)call  energy(n1,n2,n21,const_sm,me,mp, &
                       x1,x2,h1,h2,h3,Bo,rho,sm_grid,mu_grid,y, &
                       T_u2,T_u3,T_b2,T_b3,T_e0,T_eP,TE,Nprim,elec_ed, &
                       elec_ed2)
           if(index.gt.1)call  energy(n1,n2,n21,const_sm,me,mp, &
                       x1,x2,h1,h2,h3,Bo,rho,sm_grid,mu_grid,y, &
                       T_u2,T_u3,T_b2,T_b3,dummy,T_eP,TE,Nprim, &
                       elec_ed,elec_ed2)
           if(index.eq.1)T_u3b30=T_u3+T_b3
           T_e2=0.d0
           T_e=0.d0
           T_e3=0.d0
           T_e_par=0.d0 
           if(iparticle.eq.1)then
	    do 445 i=1,nproc
                call getname_mpi('d_pa',index,i-1,fname2)
                write(*,*) 'reading file ',fname2
                open(11,file=fname2,form='unformatted',status='old')
                read(11) mu_p,nu_p,vpar,mu_m,iprecip
                close(11)      
		do ii=1,np
                  call rootfind(pi,mu_p(ii),nu_p(ii),r_p, &
                                theta_p)
                  h1_p=h_mu(r_p,theta_p)
                  Bo_p=mu_o*M/(4.d0*pi*h1_p*Re**3)/Bn
                  vperp2=2.d0*Bo_p*mu_m(ii)
                  if(iprecip(ii).eq.1)then
                    T_e2=T_e2+0.5*(vpar(ii)*vpar(ii)+vperp2)* &
                        me/mp*Nprim/float(npt)
                   endif
                   if(iprecip(ii).eq.0)then
                     T_e=T_e+0.5*(vpar(ii)*vpar(ii)+vperp2)* &
                         me/mp*Nprim/float(npt)
                     T_e_par=T_e_par+ &
                         0.5*(vpar(ii)*vpar(ii))* &
                         me/mp*Nprim/float(npt)
                   endif
                enddo
 445        continue
           endif   
!          TE=TE+T_e2+T_e
           TE=TE+T_e2
           if(index.eq.imin)TE_i=TE 
            write(2,55) 1.0*float(index),T_u2,T_b2,T_e, &
                       T_e2,TE,abs((TE_i-TE)/TE_i)
!           write(2,55) 2.0*float(index),T_u3,T_b3,T_e,
!     &                 T_e2,TE,T_e0
           a=T_u3b30
           b=T_u3+T_b3
!           write(2,55) 2.0*float(index),T_u3,
!     &                 T_e0-T_eP,a,b,b-a
!          write(2,55) 1.0*float(index),T_u3,T_b3,T_u3+T_b3,
!    &                 T_eP-T_e0,T_u3+T_b3+(T_eP-T_e0)
!          write(2,55) 2.0*float(index),T_u2,T_b2,T_u2+T_b2
           call flush(2)
           call getname('d_aa',index,fname)
           open(unit=10,file=fname,status='unknown')
           call getname('d_bb',index,fname)
           open(unit=11,file=fname,status='unknown')
           dx2=x2(2)-x2(1)
           ix2=(0.0984-x2(1))/dx2+1
           do ii=1,n2
              do jj=1,n1
                 k=ii+(jj-1)*n2
                if(jj.eq.1)elec_ed(k)=0.d0
                 write(10,75) x2(ii),x1(jj),elec_ed(k)
!                write(11,75) x2(ii),x1(jj),elec_ed2(k)
                 write(11,75) x2(ii),x1(jj),sqrt(elec_ed(k)**2+ &
                              elec_ed2(k)**2)
              enddo
              write(10,*)
              write(11,*)
           enddo
           close(10)
           close(11)
!           endif
 444     continue
         close(2)
 55      format(7es11.3)
 75      format(4e12.4)
         stop
         END PROGRAM

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
          END SUBROUTINE getname


         SUBROUTINE energy(n1,n2,n21,const_sm,me,mp, &
                           x1,x2,h1,h2,h3,Bo,rho,sm_grid,mu_grid,y, &
                           T_u2,T_u3,T_b2,T_b3,T_e0,T_eP,TE, &
                           Nprim,elec_ed,elec_ed2)
         implicit none
         integer :: i,j,k,n1,n2,n21
         double precision :: TE,T_u3,T_b3,T_e0,fact,dx1,dx2,const_sm
         double precision :: T_u2,T_b2
         double precision :: T_eP
         double precision :: me,mp 
         double precision :: y(*),x1(*),x2(*),h1(*),h2(*),h3(*),rho(*)
         double precision :: sm_grid(*),Nprim,Bo(*),mu_grid(*)
         double precision :: elec_ed(*),elec_ed2(*)
         TE=0.d0
         T_u2=0.d0
         T_u3=0.d0
         T_b2=0.d0
         T_b3=0.d0
         T_e0=0.d0
         T_eP=0.d0
         Nprim=0.d0
         dx2=x2(2)-x2(1)
         do i=2,n2-1
            do j=3,n1-2
               dx1=x1(j+1)-x1(j)
               fact=1.0
               k=i+(j-1)*n2
               if(i.eq.2.or.i.eq.n2-1)fact=0.5
               if(i.eq.2.and.j.eq.3)fact=0.25
               if(i.eq.2.and.j.eq.n1-2)fact=0.25
               if(i.eq.n2-1.and.j.eq.3)fact=0.25
               if(i.eq.n2-1.and.j.eq.n1-2)fact=0.25
!               if(x2(i).lt.0.10135)fact=0.0
               T_u2=T_u2+0.5d0*fact*dx2*h2(k)*dx1*h1(k)*h3(k)* &
                    rho(k)*y(k)*y(k)
               T_u3=T_u3+0.5d0*fact*dx2*h2(k)*dx1*h1(k)*h3(k)* &
                    rho(k)*y(k+n21)*y(k+n21)
               T_b2=T_b2+0.5d0*fact*dx2*h2(k)*dx1*h1(k)* &
                    h3(k)*y(k+n21*3)*y(k+n21*3)
               T_b3=T_b3+0.5d0*fact*dx2*h2(k)*dx1*h1(k)* &
                    h3(k)*y(k+n21*4)*y(k+n21*4)
               T_e0=T_e0+fact*dx2*h2(k)*dx1*h1(k)*h3(k)*( &
                    0.5*me/mp*sm_grid(k)/const_sm+ &
                    me/mp*mu_grid(k)*Bo(k)/const_sm)
               T_eP=T_eP+fact*dx2*h2(k)*dx1*h1(k)*h3(k)* &
                    me/mp*(0.5*sm_grid(k)+ &
                           mu_grid(k)*Bo(k))/const_sm



! elec_ed - electron energy density in v1
! elec_ed2 - electron energy density in v_perp
               elec_ed(k)=0.5*me/mp*sm_grid(k)/const_sm
               elec_ed2(k)=me/mp*mu_grid(k)*Bo(k)/const_sm
               Nprim=Nprim+fact*rho(k)*h1(k)*dx1*h2(k)*dx2*h3(k)               
            enddo
         enddo
         write(*,*) T_eP
         TE=T_u2+T_b2
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

