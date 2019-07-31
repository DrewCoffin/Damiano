          SUBROUTINE ec_calc_poloidal(Fc,ec)
          use indeces
          use param
          use mesh
          implicit none
          integer :: i,j,k
          double precision :: temp1(n2),temp2(n2)
          double precision :: Fc(*),ec(*)
          do j=1,n1
             do i=1,n2
                k=i+(j-1)*n2
                temp1(i)=h1(k)*h2(k)*h3(k)*Fc(k)
             enddo
             do i=1,n2
                k=i+(j-1)*n2
                ec(k)=temp1(i)/(mphi*h1(k)*h2(k))
             enddo
          enddo
          do i=1,n2
           IF(CENTRED)THEN
              ec(i+(1-1)*n2)=-ec(i+(5-1)*n2)
              ec(i+(2-1)*n2)=-ec(i+(4-1)*n2)
              ec(i+(n1-1)*n2)=-ec(i+(n1-4-1)*n2)
              ec(i+(n1-1-1)*n2)=-ec(i+(n1-3-1)*n2)
           ELSE
              ec(i+(1-1)*n2)=-ec(i+(4-1)*n2)
              ec(i+(2-1)*n2)=-ec(i+(3-1)*n2)
              ec(i+(n1-1)*n2)=-ec(i+(n1-3-1)*n2)
              ec(i+(n1-1-1)*n2)=-ec(i+(n1-2-1)*n2)
           ENDIF
          enddo
          return
          END SUBROUTINE ec_calc_poloidal
