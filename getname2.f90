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
 10       format(a5)
          return
          END SUBROUTINE getname
