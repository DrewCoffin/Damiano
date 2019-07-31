       SUBROUTINE trid(n,a,d,c,b,x)
       implicit double precision(a-h,o-z)
       dimension a(n),b(n),c(n),d(n),x(n)
! forward elimination
       do 2 i=2,n
          xmult=a(i-1)/d(i-1)
          d(i)=d(i)-xmult*c(i-1)
          b(i)=b(i)-xmult*b(i-1)
 2     continue

! back substitution step

       x(n)=b(n)/d(n)
       
       do 3 i=n-1,1,-1
          x(i)=(b(i)-c(i)*x(i+1))/d(i)
 3     continue
       return
       END SUBROUTINE trid
