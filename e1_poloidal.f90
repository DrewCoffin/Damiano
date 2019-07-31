        SUBROUTINE epar_poloidal_1(e2,e3,e1_t)
! parallel electric field algorithm
!
        use indeces
        use param
        use mesh
        use e1terms
        use parparam
        implicit none
        integer :: i,j,k,k2
        double precision :: e2,e3,e1_t,divEp,dsm
        double precision :: dx_i,dx_ip,dx_im,beta_i,Fj
        double precision :: alpha_i,alpha_im,alpha_ip,dalpha
        double precision :: AA(n2),BB(n2),CC(n2),DD(n2),XX(n2)
        dimension e2(*),e3(*),e1_t(*)

! determine divergence
!
        do i=2,n2-1
	   do j=3,n1-2
              k=i+(j-1)*n2
              alpha_i=h2(k)/(h1(k)*h3(k))
              e3_term(i,j)=alpha_i*(h3(k+n2)*e3(k+n2)- &
                                  h3(k-n2)*e3(k-n2))/ &
                                  (x1(j+1)-x1(j-1)) 
              mu_term(i,j)=1.d0/h1(k)* &
                         (Bo(k+n2)-Bo(k-n2))/(x1(j+1)-x1(j-1))* &
                         mu_grid(k)
              mu_term2(i,j)=2.d0/h1(k)* &
                         (Bo(k+n2)-Bo(k-n2))/(x1(j+1)-x1(j-1))* &
                         mu_grid2(k)
              dsm=(sm_grid(k+n2)-sm_grid(k-n2))/(x1(j+1)-x1(j-1))
              if(j.eq.3) &
!     &           dsm=(sm_grid(k+n2)-sm_grid(k))/(x1(j+1)-x1(j))
                 dsm=(-sm_grid(k+n2*2)+4.d0*sm_grid(k+n2) &
                       -3.d0*sm_grid(k))/(x1(j+1)-x1(j-1))
              if(j.eq.n1-2) &
!     &           dsm=(sm_grid(k)-sm_grid(k-n2))/(x1(j)-x1(j-1))
                 dsm=(sm_grid(k-n2*2)-4.d0*sm_grid(k-n2) &
                      +3.d0*sm_grid(k))/(x1(j+1)-x1(j-1))
              sm_term(i,j)=1.d0/h1(k)*dsm
	   enddo
        enddo

        do j=3,n1-2
           do i=3,n2-1
              k=i+(j-1)*n2
              alpha_i=h2(k)/(h1(k)*h3(k))
              beta_i=mphi/(h2(k)*h3(k))
              sn(i,j)=beta_i*e3_term(i,j) &
                      +float(ipress)*float(icouple)*sm_term(i,j) &
                      +float(imu)*float(icouple)*mu_term(i,j) &
                      -float(imu)*float(icouple)*mu_term2(i,j) 
!              rhs_den(i,j)=h1(k)*mphi*beta_i* &
!                     (alpha_i+1.0/(h1(k)*le2(k)))
              rhs_den(i,j)=(1.0/le2(k)+mphi*alpha_i*beta_i*h1(k))
              e1_t(k)=-sn(i,j)/rhs_den(i,j)
           enddo
           k=1+(j-1)*n2
           k2=3+(j-1)*n2
           e1_t(k)=-e1_t(k2)
           k=n2+(j-1)*n2
           k2=n2-2+(j-1)*n2
           e1_t(k)=-e1_t(k2)
        enddo

       do j=1,n1
          k=n2-1+(j-1)*n2
          e1_t(k)=0.d0
       enddo

        if(ifilter.eq.1)then
           call filterd(n2,n1,e1_t,e1_t)
        endif       

!       j=2
!       do i=1,n2
!          k=i+(j-1)*n2
!          k2=i+(j+1-1)*n2
!          if(CENTRED)k2=i+(j+2-1)*n2
!          e1_t(k)=-e1_t(k2)
!          if(NEWBC)e1_t(k)=e1_t(k2)*h1(k2)/h1(k)
!       enddo

!       j=1
!       do i=1,n2
!          k=i+(j-1)*n2
!          k2=i+(j+3-1)*n2
!          if(CENTRED)k2=i+(j+4-1)*n2
!          e1_t(k)=-e1_t(k2)
!          if(NEWBC)e1_t(k)=e1_t(k2)*h1(k2)/h1(k)
!       enddo


!       j=n1-1
!       do i=1,n2
!          k=i+(j-1)*n2
!          k2=i+(j-1-1)*n2
!          if(CENTRED)k2=i+(j-2-1)*n2
!          e1_t(k)=e1_t(k2)
!          if(NEWBC)e1_t(k)=e1_t(k2)*h1(k2)/h1(k)
!       enddo

!       j=n1
!       do i=1,n2
!          k=i+(j-1)*n2
!          k2=i+(j-3-1)*n2
!          if(CENTRED)k2=i+(j-4-1)*n2
!          e1_t(k)=e1_t(k2)
!          if(NEWBC)e1_t(k)=e1_t(k2)*h1(k2)/h1(k)
!       enddo

 100   continue
       return
       END SUBROUTINE epar_poloidal_1


        SUBROUTINE epar_poloidal_2(e2,e3,ec1,e1_t)
        use indeces
        use param
        use mesh
        use e1terms
        use parparam
        implicit none
        integer :: i,j,k,k2
        double precision :: e2,e3,ec1,e1_t,divEp,dsm
        double precision :: dx_i,dx_ip,dx_im,beta_i,Fj
        double precision :: alpha_i,alpha_im,alpha_ip,dalpha
        double precision :: AA(n2),BB(n2),CC(n2),DD(n2),XX(n2)
        dimension e2(*),e3(*),ec1(*),e1_t(*)

! determine divergence
!
        do i=2,n2-1
	   do j=3,n1-2
              k=i+(j-1)*n2
              alpha_i=h2(k)/(h1(k)*h3(k))
              e3_term(i,j)=alpha_i* &
                     (h3(k+n2)*e3(k+n2)- &
                      h3(k-n2)*e3(k-n2))/ &
                         (x1(j+1)-x1(j-1)) 
              ec_term(i,j)=alpha_i* &
                     (h3(k+n2)*ec1(k+n2)- &
                      h3(k-n2)*ec1(k-n2))/ &
                         (x1(j+1)-x1(j-1))
              mu_term(i,j)=1.d0/h1(k)* &
                         (Bo(k+n2)-Bo(k-n2))/(x1(j+1)-x1(j-1))* &
                         mu_grid(k)
              mu_term2(i,j)=2.d0/h1(k)* &
                         (Bo(k+n2)-Bo(k-n2))/(x1(j+1)-x1(j-1))* &
                         mu_grid2(k)
              dsm=(sm_grid(k+n2)-sm_grid(k-n2))/(x1(j+1)-x1(j-1))
              if(j.eq.3) &
!     &           dsm=(sm_grid(k+n2)-sm_grid(k))/(x1(j+1)-x1(j))
                 dsm=(-sm_grid(k+n2*2)+4.d0*sm_grid(k+n2) &
                       -3.d0*sm_grid(k))/(x1(j+1)-x1(j-1))

              if(j.eq.n1-2) &
!     &           dsm=(sm_grid(k)-sm_grid(k-n2))/(x1(j)-x1(j-1))
                 dsm=(sm_grid(k-n2*2)-4.d0*sm_grid(k-n2) &
                      +3.d0*sm_grid(k))/(x1(j+1)-x1(j-1))
              sm_term(i,j)=1.d0/h1(k)*dsm
	   enddo
        enddo

        do j=3,n1-1
           do i=3,n2-2
              k=i+(j-1)*n2
              alpha_i=h2(k)/(h1(k)*h3(k))
              beta_i=mphi/(h2(k)*h3(k))
              sn(i,j)=beta_i*(e3_term(i,j)+ec_term(i,j)) &
                     +float(imu)*float(icouple)*mu_term(i,j) &
                     -float(imu)*float(icouple)*mu_term2(i,j) &
                     +float(ipress)*float(icouple)*sm_term(i,j)
! rhs denominotor
!             rhs_den(i,j)=h1(k)*mphi*beta_i* &
!                          (alpha_i+1.0/(h1(k)*le2(k)))
              rhs_den(i,j)=(1.0/le2(k)+mphi*alpha_i*beta_i*h1(k))
              e1_t(k)=-sn(i,j)/rhs_den(i,j)
           enddo
           k=1+(j-1)*n2
           k2=3+(j-1)*n2
           e1_t(k)=-e1_t(k2)
           k=n2+(j-1)*n2
           k2=n2-2+(j-1)*n2
           e1_t(k)=-e1_t(k2)
        enddo

        if(ifilter.eq.1)then
           call filterd(n2,n1,e1_t,e1_t)
        endif    

!       j=2
!       do i=1,n2
!          k=i+(j-1)*n2
!          k2=i+(j+1-1)*n2
!          if(CENTRED)k2=i+(j+2-1)*n2
!          e1_t(k)=-e1_t(k2)
!          if(NEWBC)e1_t(k)=e1_t(k2)*h1(k2)/h1(k)
!       enddo

!       j=1
!       do i=1,n2
!          k=i+(j-1)*n2
!          k2=i+(j+3-1)*n2
!          if(CENTRED)k2=i+(j+4-1)*n2
!          e1_t(k)=-e1_t(k2)
!          if(NEWBC)e1_t(k)=e1_t(k2)*h1(k2)/h1(k)
!       enddo


!       j=n1-1
!       do i=1,n2
!          k=i+(j-1)*n2
!          k2=i+(j-1-1)*n2
!          if(CENTRED)k2=i+(j-2-1)*n2
!          e1_t(k)=e1_t(k2)
!          if(NEWBC)e1_t(k)=e1_t(k2)*h1(k2)/h1(k)
!       enddo

!       j=n1
!       do i=1,n2
!          k=i+(j-1)*n2
!          k2=i+(j-3-1)*n2
!          if(CENTRED)k2=i+(j-3-1)*n2
!          e1_t(k)=e1_t(k2)
!          if(NEWBC)e1_t(k)=e1_t(k2)*h1(k2)/h1(k)
!       enddo


 100   continue
       return
       END SUBROUTINE epar_poloidal_2
