        SUBROUTINE distgen_vy(pi,rank,no,kT)
        use indeces
        use distgen
        implicit none
        integer :: i,no,seed,rank
        double precision :: kT,randf,temp,rho,theta
        double precision :: vbar,pi
        external randf
        vbar=0.
        seed=100001*(rank+1)
        do i=1,no
           temp=randf(seed)
           theta=2.0*pi*randf(seed)
           rho=-kT*log(temp)
           vy(i)=sqrt(2.0*rho)*cos(theta)+vbar
        enddo
        return
        END SUBROUTINE distgen_vy
