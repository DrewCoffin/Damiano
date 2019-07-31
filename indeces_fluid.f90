        MODULE indeces_fluid
        integer, parameter :: m2=7,n2=2**m2,m1=6,n1=256,n1pf=30000,n2pf=2**7
        integer, parameter :: n21=n2*n1,n=n21*5,npt=n1pf*n2pf,n_int=n1pf*10
        integer, parameter :: it=400000,nproc=1,npt2=npt/nproc
        END MODULE indeces_fluid
