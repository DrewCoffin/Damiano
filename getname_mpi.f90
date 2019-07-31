          SUBROUTINE getname_mpi(temp,index,rank,fname)
          implicit none
          integer :: one, ten, hundred,index
          integer :: rank
          character(len=4) :: temp
          character(len=11) :: fname
          hundred=index/100
          ten=(index-100*hundred)/10
          one=index - ten*10 - hundred*100
          fname(1:4)=temp(1:4)
          fname(5:5)=char(hundred+48)
          fname(6:6)=char(ten+48)
          fname(7:7)=char(one+48)
          fname(8:8)='_'
! currently will only work for less than 100 processors 
!
          hundred=rank/100
          ten=(rank-100*hundred)/10
          one=rank - ten*10 - hundred*100
          fname(9:9)=char(hundred+48)
          fname(10:10)=char(ten+48)
          fname(11:11)=char(one+48)
 10       format(a5)
          return
          END SUBROUTINE getname_mpi
