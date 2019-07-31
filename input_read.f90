          SUBROUTINE input_read(size2,rank,root)
          use param
          use parparam
          use mpi
          implicit none
          integer :: size2,rank,root,error
          character(len=10) :: blank
          namelist /INPUT/  ro1,ro2,r_ionosphere,T_e,T_i,twrite,dtwrite,&
                            twrite_p,dtwrite_p,twrite_fw,dtwrite_fw,&
                            tmax,iparticle_write,RESTART,fwhm_par,fwhm, &
                            L_perp,Amp,vp_min,vp_max
          open(unit=10,file='INPUT')
          read(10,nml=INPUT)
          close(10)
!         call MPI_BCAST(ro1,1,MPI_DOUBLE_PRECISION, &
!                        root,MPI_COMM_WORLD,error)
!         call MPI_BCAST(ro2,1,MPI_DOUBLE_PRECISION, &
!                        root,MPI_COMM_WORLD,error)
!         call MPI_BCAST(r_ionosphere,1,MPI_DOUBLE_PRECISION, &
!                        root,MPI_COMM_WORLD,error)
!         call MPI_BCAST(ke,1,MPI_DOUBLE_PRECISION, &
!                        root,MPI_COMM_WORLD,error)
!         call MPI_BCAST(T_i,1,MPI_DOUBLE_PRECISION, &
!                        root,MPI_COMM_WORLD,error)
!         call MPI_BCAST(twrite,1,MPI_DOUBLE_PRECISION, &
!                        root,MPI_COMM_WORLD,error)
!         call MPI_BCAST(dtwrite,1,MPI_DOUBLE_PRECISION, &
!                        root,MPI_COMM_WORLD,error)
!         call MPI_BCAST(twrite_p,1,MPI_DOUBLE_PRECISION, &
!                        root,MPI_COMM_WORLD,error)
!         call MPI_BCAST(dtwrite_p,1,MPI_DOUBLE_PRECISION, &
!                        root,MPI_COMM_WORLD,error)
!         call MPI_BCAST(twrite_fw,1,MPI_DOUBLE_PRECISION, &
!                        root,MPI_COMM_WORLD,error)
!         call MPI_BCAST(dtwrite_fw,1,MPI_DOUBLE_PRECISION, &
!                        root,MPI_COMM_WORLD,error)
!         call MPI_BCAST(iparticle_write,1,MPI_INTEGER, &
!                        root,MPI_COMM_WORLD,error)
          END SUBROUTINE input_read
