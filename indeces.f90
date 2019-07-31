        MODULE indeces
!*****
! n1 - number of field aligned grid points
! n2 - number of perpendicular grid points
! n1pf - number of initial macro-particles along field line
! n2pf - number of perpendicular particle positions (nproc cannot exceed n2pf)
! nbini - number of perpendicular bins in determining distribution function.
! nbinj - number of parallel bins in determining distribution functions.
!*****
#ifdef RUN_FLUID
        integer, parameter :: n2=144,n1=512,n1pf=20000,n2pf=144
        integer, parameter :: n21=n2*n1,n=n21*5,npt=n1pf*n2pf,n_int=n1pf*10
        integer, parameter :: nbini=30,nbinj=30
        integer, parameter :: it=1000000,nproc=1,npt2=npt/nproc
#else
        integer, parameter :: n2=144,n1=512,n1pf=2000000,n2pf=144
        integer, parameter :: n21=n2*n1,n=n21*5,npt=n1pf*n2pf,n_int=n1pf*10
        integer, parameter :: nbini=30,nbinj=30
        integer, parameter :: it=1000000,nproc=144,npt2=npt/nproc
#endif
        END MODULE indeces
