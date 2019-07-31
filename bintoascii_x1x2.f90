         PROGRAM bintoascii_x1x2

! Program to extract fluid and field variables from binary files and write
! them out in individual ascii files as a function of x1 and x2.

         implicit none
         integer :: ii,jj,istep,imin,imax,ichoice,index,k
         integer, parameter :: n1=512,n2=144,n21=n2*n1,n=n21*5
         double precision :: y(n),y_sum(n)
         double precision :: e1_2(n21),e2(n21),e3(n21),e1_sum(n21)
         double precision :: e2_sum(n21),Bo(n21)
         double precision :: j1_2(n21),j2(n21),j3(n21),divj_c(n21)
         double precision :: divb(n21),ec1(n21),j1x1(n21),j2x2(n21)
         double precision :: den_grid0(n21),den_grid(n21),j1p_grid(n21)
         double precision :: sm_grid(n21),mu_grid(n21),mu_grid2(n21)
         double precision :: j1p_temp(n21),x(n2,n1),z(n2,n1),h3(n2,n1)
         double precision :: x1(n1),x2(n2),dummy1,dummy2
         double precision :: rho(n21),va(n21),le2(n21)
         double precision :: j1delta1,j1delta2
         double precision :: x1dum,x2dum
! Poynting flux components
         double precision S1,S2

         character(len=7) :: fname1,fname2,fname3,fname4,fname5,fname6, &
         fname7,fname8,fname9,fname10,fname11,fname12,fname13,fname14, &
         fname15,fname16,fname17,fname18,fname19,fname,stemp, &
         fname20,fname21

         write(*,*) n1,n2,n21,n
         write(*,*) 'imin?'
         read(*,*) imin
         write(*,*) 'imax?'
         read(*,*) imax
         write(*,*) 'istep?'
         read(*,*) istep
         ichoice=3

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

         open(unit=10,file='h3.dat',status='unknown')
         do ii=1,n2
            do jj=1,n1
               read(10,*) x(ii,jj),z(ii,jj),h3(ii,jj)
            enddo
            read(10,12) stemp
         enddo
         close(10)

         open(unit=10,file='Bo.dat',status='unknown')
         do ii=1,n2
            do jj=1,n1
               k=ii+n2*(jj-1)
               read(10,*) x2dum,x1dum,Bo(k)
            enddo
!            read(10,12) stemp
         enddo
         close(10)


 12      format(a7)
         do 444 index=imin,imax,istep 
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
           write(*,*) index,j1delta1,j1delta2
!           endif

         call getname('d_u2',index,fname1)
         call getname('d_u3',index,fname2)
         call getname('d_b1',index,fname3)
         call getname('d_b2',index,fname4)
         call getname('d_b3',index,fname5)
         call getname('d_e1',index,fname6)
         call getname('d_e2',index,fname7)
         call getname('d_e3',index,fname8)
         call getname('d_j1',index,fname9)
         call getname('d_j2',index,fname10)
         call getname('d_j3',index,fname11)
         call getname('d_dj',index,fname12)
         call getname('d_db',index,fname13)
         call getname('d_et',index,fname14)
         call getname('d_dn',index,fname15)
         call getname('d_jp',index,fname16)
         call getname('d_mu',index,fname17)
         call getname('d_sm',index,fname18)
         call getname('d_ec',index,fname19)
         call getname('d_pf',index,fname20)
         call getname('d_jB',index,fname21)


         open(unit=10,file=fname1,status='unknown')
         open(unit=11,file=fname2,status='unknown')
         open(unit=12,file=fname3,status='unknown')
         open(unit=13,file=fname4,status='unknown')
         open(unit=14,file=fname5,status='unknown')
         open(unit=15,file=fname6,status='unknown')
         open(unit=16,file=fname7,status='unknown')
         open(unit=17,file=fname8,status='unknown')
         open(unit=18,file=fname9,status='unknown')
         open(unit=19,file=fname10,status='unknown')
         open(unit=20,file=fname11,status='unknown')
         open(unit=21,file=fname12,status='unknown')
         open(unit=22,file=fname13,status='unknown')
         open(unit=23,file=fname14,status='unknown')
         open(unit=24,file=fname15,status='unknown')
         open(unit=25,file=fname16,status='unknown')
         open(unit=26,file=fname17,status='unknown')
         open(unit=27,file=fname18,status='unknown')
         open(unit=28,file=fname19,status='unknown')
         open(unit=29,file=fname20,status='unknown')
         open(unit=30,file=fname21,status='unknown')

         if(imin.lt.0)goto 444
          do ii=1,n2
             do jj=1,n1
                k=ii+n2*(jj-1)
                write(10,55) x2(ii),x1(jj),y(k)
                write(11,55) x2(ii),x1(jj),y(k+n21)
                write(12,55) x2(ii),x1(jj),y(k+n21*2)
                write(13,55) x2(ii),x1(jj),y(k+n21*3)
!                write(14,55) x2(ii),x1(jj),y(k+n21*4)*h3(ii,jj)
! currently printing out MKS value of b3
!               write(14,55) x2(ii),x1(jj),y(k+n21*4)*10.d-9
                write(14,55) x2(ii),x1(jj),y(k+n21*4)
                write(15,55) x2(ii),x1(jj),e1_2(k),e1_sum(k)
                write(16,55) x2(ii),x1(jj),e2(k),e2_sum(k)
                write(17,55) x2(ii),x1(jj),e3(k)
                write(18,55) x2(ii),x1(jj),j1_2(k)
                write(19,55) x2(ii),x1(jj),j2(k)
                write(20,55) x2(ii),x1(jj),j3(k)
                write(21,55) x2(ii),x1(jj),divj_c(k)
                write(22,55) x2(ii),x1(jj),divb(k)
!                write(23,55) x2(ii),x1(jj),temp(k)
                write(24,55) x2(ii),x1(jj),den_grid(k),va(k),le2(k)
                write(25,55) x2(ii),x1(jj),j1p_grid(k)
                write(26,55) x2(ii),x1(jj),mu_grid(k),mu_grid2(k)
                write(27,55) x2(ii),x1(jj),sm_grid(k)
                write(28,55) x2(ii),x1(jj),ec1(k)
                S1=e2(k)*y(k+n21*4)
                S2=-e1_2(k)*y(k+n21*4)
                write(29,55) x2(ii),x1(jj),S1,S2
                write(30,55) x2(ii),x1(jj),j1p_grid(k)/Bo(k)
             enddo
             write(10,55)
             write(11,55)
             write(12,55)
             write(13,55)
             write(14,55)
             write(15,55)
             write(16,55)
             write(17,55)
             write(18,55)
             write(19,55)
             write(20,55)
             write(21,55)
             write(22,55)
             write(23,55)
             write(24,55)
             write(25,55)
             write(26,55)
             write(27,55)
             write(28,55)
             write(29,55)
             write(30,55)
          enddo
          
          close(10)
          close(11)
          close(12)
          close(13)
          close(14)
          close(15)
          close(16)
          close(17)
          close(18)
          close(19)
          close(20)
          close(21)
          close(22)
          close(23)
          close(24)
          close(25)
          close(26)
          close(27)
          close(28)
          close(29)
          close(30)
 444      continue
 55      format(5e15.7)
         stop
         END PROGRAM bintoascii_x1x2


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

