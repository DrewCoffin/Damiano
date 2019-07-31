        SUBROUTINE epar_1(e2,e3,e1_t)
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
              alpha_i=h3(k)/(h1(k)*h2(k))
              e2_term(i,j)=alpha_i*(h2(k+n2)*e2(k+n2)- &
                                  h2(k-n2)*e2(k-n2))/ &
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

        do i=1,n2
           do j=2,n1-1
              k=i+(j-1)*n2
              beta_i=1./(h2(k)*h3(k))
              sn(i,j)=beta_i*(e2_term(i+1,j)-e2_term(i-1,j))/ &
                              (x2(i+1)-x2(i-1)) &
                      +float(ipress)*float(icouple)*sm_term(i,j) &
                      +float(imu)*float(icouple)*mu_term(i,j) &
                      -float(imu)*float(icouple)*mu_term2(i,j) 
           enddo
        enddo

! AE_i-1+DE_i+CE_i+1=B_i
!
       do j=3,n1-2
          DD(2)=1.
          CC(2)=0.
!          BB(2)=sn(2,j)
          BB(2)=0.d0
          AA(n2-2)=0.
!          BB(n2-1)=sn(n2-2,j)
          BB(n2-1)=0.d0
          DD(n2-1)=1.
          CC(n2-1)=0.
          do i=3,n2-2             
             k=i+(j-1)*n2
             alpha_i=h3(k)/(h1(k)*h2(k))
             alpha_im=h3(k-1)/(h1(k-1)*h2(k-1))
             alpha_ip=h3(k+1)/(h1(k+1)*h2(k+1))
             dalpha=(alpha_ip-alpha_im)/(x2(i+1)-x2(i-1))
             beta_i=1./(h2(k)*h3(k))
             dx_i=x2(i+1)-x2(i-1)
             dx_ip=(x2(i+1)-x2(i))
             dx_im=(x2(i)-x2(i-1))
             Fj=dx_ip**2*dx_im+dx_im**2*dx_ip
             AA(i-1)=-beta_i*dalpha/dx_i+2.*beta_i*alpha_i*dx_ip/Fj
             DD(i)=-2.*beta_i*alpha_i*(dx_ip+dx_im)/Fj- &
                   1.d0/(h1(k)*le2(k))             
             CC(i)=beta_i*dalpha/dx_i+2*beta_i*alpha_i*dx_im/Fj
             BB(i)=sn(i,j)
          enddo

          call trid(n2-2,AA(2:n2-1),DD(2:n2-1),CC(2:n2-1), &
                    BB(2:n2-1),XX(2:n2-1))

          do i=1,n2
             k=i+(j-1)*n2
             e1_t(k)=XX(i)/(h1(k))
          enddo
          k=1+(j-1)*n2
          k2=3+(j-1)*n2
          e1_t(k)=-e1_t(k2)
          k=n2+(j-1)*n2
          k2=n2-2+(j-1)*n2
          e1_t(k)=-e1_t(k2)
       enddo
! 
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
       END SUBROUTINE epar_1


        SUBROUTINE epar_2(e2,e3,ec1,e1_t)
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
        do i=1,n2
	   do j=3,n1-2
              k=i+(j-1)*n2
              alpha_i=h3(k)/(h1(k)*h2(k))
              e2_term(i,j)=alpha_i* &
                     (h2(k+n2)*e2(k+n2)- &
                      h2(k-n2)*e2(k-n2))/ &
                         (x1(j+1)-x1(j-1)) 
              ec_term(i,j)=alpha_i* &
                     (h2(k+n2)*ec1(k+n2)- &
                      h2(k-n2)*ec1(k-n2))/ &
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

        do i=2,n2-1
           do j=3,n1-2
              k=i+(j-1)*n2
              beta_i=1./(h2(k)*h3(k))
              sn(i,j)=beta_i*((e2_term(i+1,j)-e2_term(i-1,j))/ &
                              (x2(i+1)-x2(i-1))+ &
                     float(icouple)*(ec_term(i+1,j)-ec_term(i-1,j))/ &
                              (x2(i+1)-x2(i-1))) &
                     +float(imu)*float(icouple)*mu_term(i,j) &
                     -float(imu)*float(icouple)*mu_term2(i,j) &
                     +float(ipress)*float(icouple)*sm_term(i,j)
           enddo
        enddo

! AE_i-1+DE_i+CE_i+1=B_i
!
       do j=3,n1-2
          DD(2)=1.
          CC(2)=0.
!          BB(2)=sn(2,j)
          BB(2)=0.d0
          AA(n2-1)=0.
!          BB(n2-1)=sn(n2-1,j)
          BB(n2-1)=0.d0
          DD(n2-1)=1.
          CC(n2-1)=0.
          do i=3,n2-2             
             k=i+(j-1)*n2
             alpha_i=h3(k)/(h1(k)*h2(k))
             alpha_im=h3(k-1)/(h1(k-1)*h2(k-1))
             alpha_ip=h3(k+1)/(h1(k+1)*h2(k+1))
             dalpha=(alpha_ip-alpha_im)/(x2(i+1)-x2(i-1))
             beta_i=1./(h2(k)*h3(k))
             dx_i=x2(i+1)-x2(i-1)
             dx_ip=(x2(i+1)-x2(i))
             dx_im=(x2(i)-x2(i-1))
             Fj=dx_ip**2*dx_im+dx_im**2*dx_ip
             AA(i-1)=-beta_i*dalpha/dx_i+2.*beta_i*alpha_i*dx_ip/Fj
             DD(i)=-2.*beta_i*alpha_i*(dx_ip+dx_im)/Fj- &
                   1.d0/(h1(k)*le2(k))             
             CC(i)=beta_i*dalpha/dx_i+2*beta_i*alpha_i*dx_im/Fj
             BB(i)=sn(i,j)
          enddo

!          call trid(n2,AA,DD,CC,BB,XX)

          call trid(n2-2,AA(2:n2-1),DD(2:n2-1),CC(2:n2-1), &
                    BB(2:n2-1),XX(2:n2-1))
          do i=1,n2
             k=i+(j-1)*n2
             e1_t(k)=XX(i)/(h1(k))
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
       END SUBROUTINE epar_2
