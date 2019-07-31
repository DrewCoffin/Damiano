        SUBROUTINE distgen_vz(pi,rank,no,kT)
        use indeces
        use distgen
        implicit none
        integer :: i,no,seed,rank
        double precision :: vpar0,kT,randf,temp,rho,theta
        double precision :: vbar,pi
        external randf
        seed=10000001*(rank+1)
        vbar=0.
        do i=1,no
           temp=randf(seed)
           theta=2.0*pi*randf(seed)
           rho=-kT*log(temp)
           vz(i)=sqrt(2.0*rho)*cos(theta)+vbar
        enddo
        return
        END SUBROUTINE distgen_vz
