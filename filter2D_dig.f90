        PROGRAM filter2d_dig
        implicit none
! Program filter2d_dig.f
        integer, parameter :: n2=144,n1=512
        integer :: ii,index1,index2,index,istep
        double precision :: x1(n2,n1),x2(n2,n1),data(n2,n1)
        character(len=7) :: fname
        character(len=4) :: temp
        character(len=2) :: field
        write(*,*) 'input data field (i.e. xx of d_xx)'
        read(*,5) field
 5      format(a2)
        temp(1:2)='d_' 
        temp(3:4)=field
        write(*,*) 'input first'
        read(*,*) index1
        write(*,*) 'input last'
        read(*,*) index2
        write(*,*) 'istep'
        read(*,*) istep
        index=index1
        do 100 ii=index1,index2,istep
           call getname(temp,index,fname)
           call getdata(fname,index,n2,n1,x2,x1,data)
           call filter(field,index,n2,n1,x2,x1,data)
           index=index+istep
 100    continue
        stop
        END PROGRAM filter2d_dig

        SUBROUTINE getdata(fname,index,n2,n1,x2,x1,data)
        implicit none
        integer :: ii,jj,n2,n1,index
        double precision :: x1(n2,n1),x2(n2,n1),data(n2,n1)
        character(len=7) fname,temp
        open(unit=2,file=fname,status='old')
        write(*,*) 'reading...',fname
        write(*,*) n2,n1
        do  ii=1,n2
           do jj=1,n1
              read(2,*) x2(ii,jj),x1(ii,jj),data(ii,jj)
           enddo
           read(2,10) temp
        enddo
 10     format(a7)
        close(2)
        return
        END SUBROUTINE getdata

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

        SUBROUTINE filter(field,index,n2,n1,x2,x1,sn)
        implicit none
        integer :: i,j,n2,n1,index,itemp
        double precision :: x1(n2,n1),x2(n2,n1)
        double precision sn(n2,n1),sn2(n2,n1),sn3(n2,n1)
        character(len=7) :: fname2,fname3
        character(len=4) :: temp2,temp3
        character(len=2) :: field
        write(*,*) index
        call filterd(n2,n1,sn,sn2)
        temp2(1:2)='df'
        temp2(3:4)=field
        temp3(1:2)='f_'
        temp3(3:4)=field
        call getname(temp2,index,fname2)
        call getname(temp3,index,fname3)
        open(unit=3,file=fname2,status='unknown')
        open(unit=4,file=fname3,status='unknown')
        do i=2,n2-1
        do j=3,n1-2
           write(3,*) x2(i,j),x1(i,j),sn2(i,j)
           if(j.eq.n1-3) write(4,*) x2(i,n1-2),sn2(i,j)
         enddo
        write(3,*)
        enddo
        close(3)
        close(4)
        return
        END SUBROUTINE filter

        SUBROUTINE filterd(n2,n1,sn,sn2)
        implicit none
        integer :: i,j,k,n1,n2,k2
        double precision :: sn(*),sn2(*)
        double precision :: sn_t1(n2,n1),sn_t2(n2,n1)

! 1/4-1/2-1/4 digital filter as per Birdsall and Langdon Appendix C

        do j=1,n1
           do i=1,n2
              k=i+(j-1)*n2
              sn_t1(i,j)=sn(k)
           enddo
        enddo
       
! filter in x2 direction
        do j=2,n1-1
           do i=2,n2-1
              sn_t2(i,j)=0.25*sn_t1(i-1,j)+0.5*sn_t1(i,j) &
                           +0.25*sn_t1(i+1,j)
           enddo
        enddo          
        
! filter in x1 direction
        do j=2,n1-1
           do i=2,n2-1
              k=i+(j-1)*n2
              sn2(k)=0.25*sn_t2(i,j+1)+0.5*sn_t2(i,j) &
                       +0.25*sn_t2(i,j-1)
           enddo
        enddo

        return
        END SUBROUTINE filterd
 


