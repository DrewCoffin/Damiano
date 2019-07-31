        PROGRAM xslice
!
! program to take slice along a specified value of x2 from the ascii files 
! put out by bintoascii_x1x2.f. NOTE: By default the program averages the
! results for 1 adjacent field lines to either side of the one chosen.
!
! to do: move choice of field line to main program
!
! modified Feb 2013 to allow for input of both filtered and unfiltered data where
! filtering is being done by filter2d_dig.f
!
        implicit none
        integer, parameter :: n1=512,n2=144
        integer :: index1,index2,index,nof,ii,jj,istep,ifilter
        integer :: ix2
        real :: x1(n1),x2(n2),dummy1,dummy2
        real :: x(n2,n1),z(n2,n1),h3(n2,n1)
        real :: dx2
        character(len=7) :: fname,stemp
        character(len=4) :: temp
        character(len=2) :: field
        write(*,*) 'input data field (i.e. xx of d_xx)'
        read(*,5) field
 5      format(a2)
        write(*,*) 'filtered data? (0-no,1-yes)'
        read(*,*) ifilter
        if(ifilter.eq.1)then
           temp(1:2)='df'
        else
           temp(1:2)='d_'
        endif
        temp(3:4)=field
        write(*,*) 'input first'
        read(*,*) index1
        write(*,*) 'input last'
        read(*,*) index2
        write(*,*) 'istep'
        read(*,*) istep
        index=index1
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
         dx2=x2(2)-x2(1)
!        ix2=(0.103-x2(1))/dx2+1
         ix2=n2/2
! for e2
!       ix2=(0.1028-x2(1))/dx2+1
         open(unit=10,file='h3.dat',status='unknown')
         do ii=1,n2
            do jj=1,n1
               read(10,*) x(ii,jj),z(ii,jj),h3(ii,jj)
            enddo
            read(10,12) stemp
         enddo
         close(10)
 12      format(a7)


        do 100 ii=index1,index2,istep
           call getname(temp,index,fname)
           if(ifilter.eq.1)then
             call getdata_df(fname,index,ix2,x2,x1,x,z,n2,n1)
           else
             call getdata(fname,index,ix2,x2,x1,x,z,n2,n1)
           endif
           index=index+istep 
 100    continue
        stop
        END PROGRAM xslice


        SUBROUTINE getdata(fname,index,ix2,x2,x1,x,z,n2,n1)
        implicit none
        integer :: ii,jj,n2,n1,index,islice,ix2
        real :: x(n2,n1),z(n2,n1),val(n2,n1),dx,pi,aval(n1)
        real :: xm,zm,lambda1,lambda2,r,ro,length
        real :: x1(*),x2(*),dx2,dummy1,dummy2
        character(len=7) :: fname,fname2,temp
        character(len=2) :: field
        open(unit=2,file=fname,status='old')
        write(*,*) 'reading...',fname
        fname(1:1)='s'
        call getname(fname(1:4),index,fname2)
        open(unit=3,file=fname2,status='unknown')
        pi=4.0*atan(1.0)
        lambda1=0.d0
        do  ii=1,n2
           do jj=1,n1
              read(2,*) dummy1,dummy2,val(ii,jj)
           enddo
           read(2,10) temp
        enddo
 10     format(a7)
        ii=ix2
          do jj=n1/2,n1-2
              lambda2=atan(z(ii,jj)/x(ii,jj))
              r=z(ii,jj)/sin(lambda2)
              ro=r/(cos(lambda2))**2
              write(3,*) abs(length(ro,lambda1,lambda2)),val(ix2,jj), &
                         x1(jj),float(jj)
          enddo
        close(2)
        close(3)
        return
        END SUBROUTINE

        SUBROUTINE getdata_df(fname,index,ix2,x2,x1,x,z,n2,n1)
        implicit none
        integer :: ii,jj,n2,n1,index,islice,ix2
        real ::  x(n2,n1),z(n2,n1),val(n2,n1),dx,pi,aval(n1)
        real :: xm,zm,lambda1,lambda2,r,ro,length
        real :: x1(*),x2(*),dx2,dummy1,dummy2
        character(len=7) :: fname,fname2,temp
        character(len=2) :: field
        open(unit=2,file=fname,status='old')
        write(*,*) 'reading...',fname
        fname(1:1)='s'
        call getname(fname(1:4),index,fname2)
        open(unit=3,file=fname2,status='unknown')
        pi=4.0*atan(1.0)
        lambda1=0.d0
        do  ii=2,n2-1
           do jj=3,n1-2
              read(2,*) dummy1,dummy2,val(ii,jj)
           enddo
           read(2,10) temp
        enddo
 10     format(a7)
        ii=ix2
          do jj=n1/2,n1-2
              lambda2=atan(z(ii,jj)/x(ii,jj))
              r=z(ii,jj)/sin(lambda2)
              ro=r/(cos(lambda2))**2
              call average(ix2,n2,n1,val,aval)
              write(3,*) abs(length(ro,lambda1,lambda2)),val(ii,jj), &
                         x1(jj),float(jj)
          enddo
        close(2)
        close(3)
        return
        END SUBROUTINE getdata_df


        SUBROUTINE getname(temp,index,fname)
        implicit none
        integer :: one, ten, hundred,index
        character(len=4) :: temp
        character(len=7) fname
        hundred=index/100
        ten=(index-100*hundred)/10
        one=index - ten*10 - hundred*100
        fname(1:4)=temp(1:4)
        fname(5:5)=char(hundred+48)
        fname(6:6)=char(ten+48)
        fname(7:7)=char(one+48)
 10     format(a5)
        return
        END SUBROUTINE getname

        FUNCTION length(ro,a,b)
        implicit none
        real :: ro,lambda,length,pi,a,b,asinh
        pi=4.0*atan(1.0d0)
! asinh= ln(x+sqrt(x*x+1))
        length=-(ro*asinh(sqrt(3.0)*sin(a)))/(2.0*sqrt(3.0))+ &
               ro*asinh(sqrt(3.0)*sin(b))/(2.0*sqrt(3.0))- &
               ro*sin(a)*sqrt(1+3.0*(sin(a))**2)/2.0+ &
               ro*sin(b)*sqrt(1+3.0*(sin(b))**2)/2.0
        return
        END FUNCTION length

        FUNCTION asinh(x)
        implicit none
        real :: x,asinh
        asinh=log(x+sqrt(x*x+1.d0))
        return
        END FUNCTION asinh
        
        SUBROUTINE average(ix2,n2,n1,val,aval)
        integer :: ix2,n2,n1
        real :: val(n2,n1),aval(n1)
        do i=1,n2
           aval(i)=0.d0
        enddo

        do j=1,n1 
        do i=ix2-1,ix2+1
           aval(j)=aval(j)+val(i,j)
        enddo
        aval(j)=aval(j)/float(11)
        enddo       
        return
        END SUBROUTINE average
