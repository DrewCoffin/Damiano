        SUBROUTINE eperp(y,e2,e3)
        use indeces
        use param
        use mesh
        implicit none
        integer :: ii,jj,k,k2
        double precision :: y,e2,e3
        dimension :: y(*),e2(*),e3(*)
	do jj=3,n1-2
           do ii=2,n2-1
              k=ii+(jj-1)*n2
              e2(k)=-y(k+n21)*Bo(k)
              e3(k)=y(k)*Bo(k)
	   enddo
           k=1+n2*(jj-1)
           k2=3+n2*(jj-1)
           e2(k)=h1(k2)*h3(k2)*e2(k2)/(h1(k)*h3(k))
           k=n2+n2*(jj-1)
           k2=n2-2+n2*(jj-1)
           e2(k)=h1(k2)*h3(k2)*e2(k2)/(h1(k)*h3(k))
        enddo
        do ii=1,n2
           IF(CENTRED)THEN
              e2(ii+(1-1)*n2)=-e2(ii+(5-1)*n2)
              e2(ii+(2-1)*n2)=-e2(ii+(4-1)*n2)
              e2(ii+(n1-1)*n2)=-e2(ii+(n1-4-1)*n2)
              e2(ii+(n1-1-1)*n2)=-e2(ii+(n1-3-1)*n2)
           ELSE
              e2(ii+(1-1)*n2)=-e2(ii+(4-1)*n2)
              e2(ii+(2-1)*n2)=-e2(ii+(3-1)*n2)
              e2(ii+(n1-1)*n2)=-e2(ii+(n1-3-1)*n2)
              e2(ii+(n1-1-1)*n2)=-e2(ii+(n1-2-1)*n2)
           ENDIF
        enddo
        return
        END SUBROUTINE eperp
