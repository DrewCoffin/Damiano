          PROGRAM main

! driver for program dipolar which is a 2D hybrid MHD-kinetic electron
! model of a standing shear Alfven wave in dipolar coordinates. The model 
! is 2D in the x1-x2 directions and assumes that k_x3=0.
! Peter Damiano, University of Alberta 2002.

          use indeces
          use param
          use fields
          use mesh
          use parparam
          use e1terms
          use distgen
          use df
          use mpi
          implicit none
          integer :: i,k,imin,ii,jj,kk,j,imode,i_acount,n_acount,iter
          integer :: rank,error,size2,root
          double precision :: h_mu,mu,Bfield,h_nu,randf,delta
          double precision :: h1_temp,Bo_temp,ro_r,temp,temp_exp
          double precision :: e1a(n21),e1a_f(n21),rtemp
          character(len=7) :: fname
          external h_mu,mu,Bfield,h_nu,randf
          logical :: NO_FLUID,LOAD_DISTR,PTEST1,PTEST2,PTEST3
          logical :: DELTAF_DISTR,PTEST4

          call MPI_INIT(error)
          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,error)
          call MPI_COMM_SIZE(MPI_COMM_WORLD,size2,error)
          root=0
          size2=nproc

          call field_arrays()
! parameters

          pi=4.0*atan(1.0)
!         write(*,*) size2,nproc
          if(size2.ne.nproc)stop
!
! n1 - number of points along field line
! n2 - number of points in the radial direction
!
! imode = modes
! imode = 1 run as with particles dictating density profile (for hybrid model)
! imode = 2 run with specified density profile 
          version=2019.1
          imode=1
          index_df=0
          index=1
          index_p=1
          i_inertia=0
#ifdef RUN_FLUID
          iparticle=0
          icouple=0
          imode=2
#else
          i_inertia=1
          iparticle=1
          icouple=1
#endif

#ifdef GYRO
          igyro=1
#else
          igyro=0
#endif

          ipress=1
          imu=1
          ipwrite=0
          iprofile=1
          ifilter=0
          i_acount=1
!          n_acount=20
          seed=10001
          nbin=100
!         RESTART=.FALSE.
          NO_FLUID=.FALSE.
          NO_TC=.FALSE.

          if(NO_FLUID)i_inertia=0
!*temporarily disconnected
! PTEST2 uses sp_init - use imode=2
! PTEST tracks first 5 particles - use imode=1
! PTEST3 rearranges n1p electrons form middle field line to first 
!        n1p in the array. Prints out 5. - use with imode=1
! PTEST4 tracks 5 electrons along middle field line. No readjustment 
!        of np value. Use with imode=1
!
!         PTEST1=.FALSE.
!         PTEST2=.FALSE.
!         PTEST3=.FALSE.
! originally true
! CHANGED
!         PTEST4=.TRUE.
!         if(PTEST1.and.PTEST2)then
!           write(*,*) '2 tests modes active...quitting'
!           stop
!         endif
!*temporarily disconected
          DELTAF_DISTR=.FALSE.
          ETEST=.FALSE.
          ETEST2=.FALSE.
          PRECIP=.FALSE.
          REFLECT=.FALSE.
          CENTRED=.TRUE.
          NEWBC=.TRUE.
          if((.NOT.PTEST1).and.(.NOT.PTEST2))ETEST=.FALSE.
          call input_read(size2,rank,root)
!         call MPI_BCAST(ro1,1,MPI_DOUBLE_PRECISION, &
!                        root,MPI_COMM_WORLD,error)

! make sure ro from sp_init is in range of grid
          if(PTEST2)then
            ro2=9.7
            ro1=10.3
          endif
          ro_r=(ro1+ro2)/2.d0
          ro_d_min=ro_r-0.04
          ro_d_max=ro_r+0.04
          n2p=n2pf
          n1p=n1pf
          np=n1p*n2p

!         pa_min=-pi*90./180.
!         pa_max=pi*180./180.
! ipa_profile=1 uniform pa profile
! ipa_profile=2 linear pa profile
! ipa_profile=3 weight as sin(alpha)
!         ipa_profile=3


! physical constants all in MKS units
!
          M=8.02d22
          Re=6.38d6
#ifdef JUPITER
          M=1.55e27
          Re=7.1492e7
#endif
          L=Re
          me=9.109d-31
          mp=1.67d-27*23.0
          c=2.998d8
          e=1.609d-19
!         mphi=250.
          mphi=0.d0
          mu_o=4.d-7*pi
          epsilon_o=8.85d-12
          T_e=T_e*1.609d-19
! in subroutine to generate electron distribution T_e is actually T_e/2 
          T_e=T_e/2.
          kT=1.0
          vth=sqrt(2.*T_e/me)
! convert eV to Joules
          T_i=T_i*1.609e-19


          t=0.

          call grid(size2,rank,root)

! normalized speed of light
         cn=c/vn

! define time step and writeout times
!
          tmax=tmax/tn
          dt=tmax/float(it)
          twrite=twrite/tn
          dtwrite=dtwrite/tn
          twrite_p=twrite_p/tn
          dtwrite_p=dtwrite_p/tn
          twrite_fw=twrite_fw/tn
          dtwrite_fw=dtwrite_fw/tn
          if(iparticle.eq.0)then
            twrite_fw=200.d0/tn
            dtwrite_fw=200.0d0/tn
          endif

          nt=it*(twrite/tmax)
          nt_p=it*(twrite_p/tmax)
          nt_fw=it*(twrite_fw/tmax)
          n_acount=5.0/tn*it/tmax


          if(imode.eq.1)then
             call p_initialize(size2,rank,root)
             call distgen_vx(pi,rank,np2,kT)
             call distgen_vy(pi,rank,np2,kT)
             call distgen_vz2(np2,kT)
             dtheta_p=pi/2/100
             call v_initialize(imode,PTEST1,PTEST3,DELTAF_DISTR,ro_r)
             call pcurrent0(size2,rank,root)
             IF(rank.eq.root)THEN
                call density2()
                call Va_profile()
                if(i_inertia.eq.1)call le_mks()
                call fl_normalize()
                call scaling()
                if(igyro.eq.1)call rho_i_calc()
                index=0
                call fluid_write(index)
                index=1
!                call j1delta_det()
              ENDIF
!              call MPI_BCAST(j1delta,1,MPI_DOUBLE_PRECISION,
!     &                  root,MPI_COMM_WORLD,error)
          elseif(imode.eq.2)then
             IF(rank.eq.root)THEN
                call density1()
                call Va_profile()
                if(i_inertia.eq.1)call le_mks()
                call fl_normalize()
                if(igyro.eq.1)call rho_i_calc()
             ENDIF
          else
             write(*,*) 'imode out of bounds...quitting.'
             stop
          endif
! send Bo info to all ranks
          call MPI_BCAST(Bo,n21,MPI_DOUBLE_PRECISION, &
                        root,MPI_COMM_WORLD,error)
          call MPI_BCAST(dBo,n21,MPI_DOUBLE_PRECISION, &
                        root,MPI_COMM_WORLD,error)


          open(19,file='par.dat',form='unformatted',status='unknown')
          write(19) n1,n2,n21,n
          close(19)

          if(rank.eq.root)then
              write(*,*) 'GKE - Parameters'
              write(*,*) '---------'
              write(*,*) 'n1p,n2,np',n1p,n2p,np
              write(*,*) 'n1,n2',n1,n2
              write(*,*) 'nproc=', nproc
              write(*,*) '---------'
          endif

       if(PTEST1.and.iparticle.eq.1)np=600
       if(PTEST2.and.iparticle.eq.1)then
             call sp_init()
             call p_normalize()
             do i=1,5
                write(*,*) vpar0(i)
             enddo
             n1p=np
       endif
       if(PTEST3.and.iparticle.eq.1)then
         do i=1,n1p
            k=(n2p/2-1)*n1p+i
            r_p(i)=r_p(k)
            theta_p(i)=theta_p(k)
            mu2_p(i)=mu2_p(k)
            mu1_p(i)=mu1_p(k)
            ro_p(i)=ro_p(k)
            nu_p(i)=nu_p(k)
            vpar0(i)=vpar0(k)
            vpar1(i)=vpar1(k)
            mu_m(i)=mu_m(k)
            np=n1p
         enddo
       endif
!
#ifdef NCAR

#else
       if(iparticle.eq.1)then
             vmax=0.
             mu_max=0.
             mu_min=0.
             do ii=1,np
                vmax=max(vmax,vpar0(ii))
                mu_max=max(mu_max,mu_m(ii))
             enddo
             vmax=2.*vmax
             vmin=-vmax
             if(DELTAF_DISTR)then
               vmax=1.0
               vmin=-1.0
             endif
             vmax=1.0
             vmin=-1.0
!             call fdistrb(vpar0,ro_p)
!             call fdistr_iN(vpar0,theta0_p)
!             call fdistr_iS(vpar0,theta0_p)
             call padistr(size2,rank,root)
             call mudistr(size2,rank,root)
       endif
#endif
!      call full_write(i+1)

! initialize arrays
!
          do i=1,n
             f_1(i)=0.
             f_2(i)=0.
             y(i)=0.
             y2(i)=0.
             y0(i)=0.
          enddo
          do i=1,n21
             e1_sum(i)=0.d0
             e2_sum(i)=0.d0
             e2(i)=0.
             e3(i)=0.
          enddo
          do i=1,n
             y_sum(i)=0.d0
          enddo

        do i=1,n21
           Fc1(i)=0.d0
           Fc0(i)=0.d0
        enddo
! initial perturbation
!
          if(rank.eq.root)then
#ifdef PULSE
          if(.NOT.NO_FLUID)call initial_pulse()
#else
          if(.NOT.NO_FLUID)call initial_flr()
#endif
#ifdef KFLR
          if(.NOT.NO_FLUID)call initial_epsilon5()
#endif
          endif

! call input_out
!
          if(rank.eq.root)call input_out()

        if(RESTART)then
           call full_read(size2,rank,root,imin)
           twrite_fw=twrite_fw + dtwrite_fw
           nt_fw=it*(twrite_fw/tmax)
        else
           imin=0
        endif

          do i=imin,it
! predictor step
!
             if(NO_FLUID) goto 99
             if(rank.eq.root)call F1()
             if(rank.eq.root)call pred(dt)
             if(rank.eq.root)call bc1()
 99          if(iparticle.eq.1) then
                call push1(size2,rank,root)
                call pcurrent_p(size2,rank,root)
             endif
             if(rank.eq.root)call current(y2,j1_2,j2,j3)
             if(icouple.eq.1)then
                if(rank.eq.root)call div_j(j1p_grid,j2,j3, &
                                           j1x1,j2x2,j3x3,divj_p)
             else
                if(rank.eq.root)call div_j(j1_2,j2,j3, &
                                           j1x1,j2x2,j3x3,divj_p)
             endif
#ifdef GYRO
             if(rank.eq.root)call eperp_gyro(y2,e2,e3)
#else
             if(rank.eq.root)call eperp(y2,e2,e3)
#endif
             if(i_inertia.eq.1) then
!                call epar1()
!                call epar2()
#ifdef POLOIDAL
               if(rank.eq.root)call epar_poloidal_1(e2,e3,e1_2)
#else
               if(rank.eq.root)call epar_1(e2,e3,e1_2)
#endif

             endif
! corrector step
!
             iter=0
 1000        continue
             if(NO_FLUID) goto 199
             if(rank.eq.root)call F2()
             if(rank.eq.root)call corr(dt)
             if(rank.eq.root)call bc2()
 199         if(iparticle.eq.1) then
                call push2(size2,rank,root)
                call pcurrent_c(size2,rank,root)
             endif
             if(rank.eq.root)call current(y,j1_2,j2,j3)
             if(icouple.eq.1)then
                if(rank.eq.root)call div_j(j1p_grid,j2,j3, &
                                           j1x1,j2x2,j3x3,divj_c)
             else
                if(rank.eq.root)call div_j(j1_2,j2,j3, &
                                           j1x1,j2x2,j3x3,divj_c)
             endif
             IF(rank.eq.root)then
             if(iparticle.eq.1.and.icouple.eq.1)then
             do ii=1,n2
                do jj=1,n1
                 kk=ii+(jj-1)*n2
                 Fc1(kk)=Fc0(kk)-dt*const_ec*0.5*(divj_c(kk)+divj_p(kk))
                 Fc0(kk)=Fc1(kk)
                enddo
             enddo
#ifdef POLOIDAL
             call ec_calc_poloidal(Fc1,ec1)
#else
             call ec_calc(Fc1,ec1)
#endif

!             do ii=1,n2
!                do jj=1,n1
!                 kk=ii+(jj-1)*n2
!                 rho(kk)=den_grid(kk)
!                 va(kk)=Bo(kk)/sqrt(den_grid(kk))
!                 le2(kk)=me*c*c/
!     &                  (4.*pi*den_grid(kk)*rho_n/mp*e*e*L*L)
!                enddo
!             enddo
             ENDIF
!	     if(ifilter.eq.1)then
!               call filterd(n2,n1,ec1,ec1)
!             endif
             endif
#ifdef GYRO
             if(rank.eq.root)call eperp_gyro(y,e2,e3)
#else
             if(rank.eq.root)call eperp(y,e2,e3)
#endif
!             if(iter.lt.3)then
!	       iter=iter+1
!               goto 1000
!             endif             

!              open(unit=10,file='test.dat',status='unknown')
!              do ii=1,n2
!                 do jj=1,n1
!                    kk=ii+(jj-1)*n2
!                    write(10,*) x(ii,jj),z(ii,jj),ec1(kk)
!                 enddo
!                 write(10,*)
!              enddo
!              close(10)

             if(i_inertia.eq.1) then
!                call epar1()
!                call epar2()
                 if(rank.eq.root)then
#ifdef POLOIDAL
                   call epar_poloidal_2(e2,e3,ec1,e1_2)
#else
                   call epar_2(e2,e3,ec1,e1_2)
#endif
                   do ii=1,n21
                      e1_sum(ii)=e1_sum(ii)+e1_2(ii)
                      e2_sum(ii)=e2_sum(ii)+e2(ii)
                   enddo
                   do ii=1,n
                      y_sum(ii)=y_sum(ii)+y(ii)
                   enddo
                 endif
             endif
            IF(rank.eq.root)then       
             if(i_inertia.eq.1.and.i_acount.le.n_acount-1)then
               do ii=1,n21
                  e1a(ii)=e1a(ii)+e1_2(ii)
               enddo
               i_acount=i_acount+1
             elseif(i_inertia.eq.1.and.i_acount.eq.n_acount)then
               do ii=1,n21
                  e1a_f(ii)=(e1a(ii)+e1_2(ii))/float(n_acount)
                  e1a(ii)=0.d0
               enddo
               i_acount=1
             endif
            ENDIF
               
               
               

! write out fluid fields at t=twrite
!
! do particles first before advancing index variable in next loop.
!
              if(i.eq.nt_p)then
               if(iparticle.eq.1)call fdistr2D(size2,rank,root)
               if(rank.eq.root)call df_write(index)
               if(iparticle.eq.1.and.iparticle_write.eq.1) &
                 call particle_write(size2,rank,root)
               twrite_p=twrite_p+dtwrite_p
               write(*,*) nt,nt_p,index
               nt_p=it*(twrite_p/tmax)
              endif

              if(i.eq.nt)then
!               if(rank.eq.root)write(*,*) index
               if(rank.eq.root)call div_b(y,divb)
               if(rank.eq.root)call fluid_write(index)
               if(rank.eq.root)call e1_term_write(index)
               index=index+1
               twrite=twrite+dtwrite
               nt=it*(twrite/tmax)
               index_df=index_df+1
!              if(rank.eq.root)call fdistrb(vpar2,ro_p)
              endif

! particle write
! commented out Tues Apr 17 2007 - program failed on new cluster because of a 
! writing problem associated with output1.dat
!
!           IF(rank.eq.root)then
!           if(iparticle.eq.1.and.ipwrite.eq.200)then
!             open(unit=21,file='output1.dat',status='unknown')
!             open(unit=22,file='output2.dat',status='unknown')
!             open(unit=23,file='output3.dat',status='unknown')
!             open(unit=24,file='output4.dat',status='unknown')
!             open(unit=25,file='output5.dat',status='unknown')
!             do ii=1,100000000
!                read(21,*,end=31) rtemp
!             enddo
! 31          do ii=1,100000000
!                read(22,*,end=32) rtemp
!             enddo
! 32          do ii=1,100000000
!                read(23,*,end=33) rtemp
!             enddo
! 
! 33          do ii=1,100000000
!                read(24,*,end=34) rtemp
!             enddo
! 34          do ii=1,100000000
!                read(25,*,end=35) rtemp
!             enddo
! 35          continue
!c              call particle_write()
!c              index_p=index_p+1
!c for np=1
!             if(np.eq.1)
!     &          write(21,55)t*tn, r_p(1),pi/2.-theta_p(1),
!     &                    r_p(1)*sin(theta_p(1)),
!     &                    r_p(1)*cos(theta_p(1)),vpar2(1)


!c for PTEST2
!            if(PTEST2.and.np.eq.5)then
!                do ii=1,5
!                  kk=ii
!                  write(20+ii,55)t*tn, r_p(kk),pi/2.-theta_p(kk),
!     &                    r_p(kk)*sin(theta_p(kk)),
!     &                    r_p(kk)*cos(theta_p(kk)),vpar2(kk)
!                enddo
!             endif
!
!
!
!c for PTEST1
!             if(PTEST1)then
!                do ii=1,5
!                  kk=ii*100
!                  write(20+ii,55)t*tn, r_p(kk),pi/2.-theta_p(kk),
!     &                    r_p(kk)*sin(theta_p(kk)),
!     &                    r_p(kk)*cos(theta_p(kk)),vpar2(kk)
!                enddo
!             endif
!
!             if(PTEST3)then
!                do ii=1,5
!                  kk=ii*n1p/5
!                  write(20+ii,55)t*tn, r_p(kk),pi/2.-theta_p(kk),
!     &                    r_p(kk)*sin(theta_p(kk)),
!     &                    r_p(kk)*cos(theta_p(kk)),vpar2(kk)
!                enddo
!             endif
!
!             if(PTEST4)then
!                do ii=1,5
!                  kk=n2p*n1p/2+ii*n1p/5
!                  write(20+ii,55)t*tn, r_p(kk),pi/2.-theta_p(kk),
!     &                    r_p(kk)*sin(theta_p(kk)),
!     &                    r_p(kk)*cos(theta_p(kk)),vpar2(kk)
!                enddo
!             endif
!
!c                kk=ii_vpamax
!                write(21,55)t*tn, r_p(kk),pi/2.-theta_p(kk),
!c     &                    r_p(kk)*sin(theta_p(kk)),
!c     &                    r_p(kk)*cos(theta_p(kk)),vpar2(kk)
!
!             ipwrite=0
!             close(21)
!             close(22)
!             close(23)
!             close(24)
!             close(25)
!           endif
!           ENDIF

! update particle fields

           if(rank.eq.root)call field_update()
           if(iparticle.eq.1)call particle_update(size2,rank,root)
           if (i.eq.nt_fw) then
              index_fw=index_fw+1
              call full_write(size2,rank,root,i+1)
              twrite_fw=twrite_fw + dtwrite_fw
              nt_fw=it*(twrite_fw/tmax)
            endif
           ipwrite=ipwrite+1
           t=t+dt
#ifdef PULSE
           if(t.gt.40.d0/tn)stop
#else
           if(t.gt.140.d0/tn)stop
#endif
          enddo
 55       format(6e13.5)
          call MPI_FINALIZE(error)
          stop
          END PROGRAM main
