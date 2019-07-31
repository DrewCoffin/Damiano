        SUBROUTINE distgen_vx(pi,rank,no,kT)
        use indeces
        use distgen
        implicit none
        integer :: i,no,seed,rank
        double precision :: kT,randf,temp,rho,theta
        double precision :: vbar,pi
        external randf
        seed=1001*(rank+1)
        vbar=0.
        do i=1,no
           temp=randf(seed)
           theta=2.0*pi*randf(seed)
           rho=-kT*log(temp)
           vx(i)=sqrt(2.0*rho)*cos(theta)+vbar
        enddo
        return
        END SUBROUTINE distgen_vx

        FUNCTION randf(iseed)
        implicit double precision(a-h,o-z)
! nov 4 notes
        data ia,ic,iq,ir /16807,2147483647,127773,2836/
! Note IQ=IC/IA, IR=MODULE(IC,IA)
!        common /sd/iseed
        ih=iseed
        il=mod(iseed,iq)
! mod(n,m)=n-[n/m]*m
        it=ia*il-ir*ih
        if(it.gt.0)then
          iseed=it
        else
          iseed=ic+it
        endif
        randf=iseed/float(ic)
        return
        END FUNCTION randf

