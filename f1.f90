         SUBROUTINE F1()
	 use indeces
         use param
         use mesh
         use fields
         implicit none
         integer :: i,j,k,ku2,ku3,kb1,kb2,kb3
         integer :: ku2_2,ku3_2,kb1_2,kb2_2,kb3_2
         double precision :: term1,term2,term3,term4,term5,term6,term7
         do  j=3,n1-2
         do  i=2,n2-1
! equation 18

             ku2=i+n2*(j-1)
             k=ku2
             ku3=ku2+n21
             kb1=ku2+n21*2
             kb2=ku2+n21*3
             kb3=ku2+n21*4
             term1=(h2(k+n2)*y(kb2+n2)-h2(k-n2)*y(kb2-n2))/ &
                   (x1(j+1)-x1(j-1))
             term2=(h1(k+1)*y(kb1+1)-h1(k-1)*y(kb1-1))/ &
                   (x2(i+1)-x2(i-1))
             term3=(h3(k+n2)*y(kb3+n2)-h3(k-n2)*y(kb3-n2))/ &
                   (x1(j+1)-x1(j-1))
             term4=(h3(k+1)*e3(k+1)-h3(k-1)*e3(k-1))/ &
                   (x2(i+1)-x2(i-1))    
             term5=(h3(k+n2)*e3(k+n2)-h3(k-n2)*e3(k-n2))/ &
                   (x1(j+1)-x1(j-1))
             term6=(h2(k+n2)*e2(k+n2)-h2(k-n2)*e2(k-n2))/ &
                   (x1(j+1)-x1(j-1))
             term7=(h1(k+1)*e1_1(k+1)-h1(k-1)*e1_1(k-1))/ &
                   (x2(i+1)-x2(i-1))
! equation for du2/dt
#ifdef POLOIDAL
             f_1(ku2)=Bo(k)/(rho(k)*h1(k)*h2(k))* &
                   (term1-term2)
#endif
! equation for du3/dt
#ifdef TOROIDAL
             f_1(ku3)=-Bo(k)/(rho(k)*h1(k)*h3(k))*gamma2(k)* &
                   (mphi*h1(k)*y(kb1)-term3)
#endif
! equation for db1/dt
#ifdef B1
            f_1(kb1)=-1/(h2(k)*h3(k))* &
                  (term4+mphi*h2(k)*e2(k))
#endif
! equation for db2/dt
#ifdef POLOIDAL
             f_1(kb2)=1./(h1(k)*h3(k))* &
                   (mphi*h1(k)*e1_1(k)+term5) 
#endif
! equation for db3/dt
#ifdef TOROIDAL
             f_1(kb3)=-1./(h1(k)*h2(k))*(term6-term7)
#endif
           enddo
          enddo
          return
          END SUBROUTINE f1
