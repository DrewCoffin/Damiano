        SUBROUTINE distgen_extra(pi,localseed,vtemp)
! generate new vx in case pa is than specified limit
        implicit none
        integer :: i,no,localseed,rank
        double precision :: kT,randf,temp,rho,theta,vtemp
        double precision :: vbar,pi
        external randf
        vbar=0.
        kT=1.d0
        do i=1,1
           temp=randf(localseed)
           theta=2.0*pi*randf(localseed)
           rho=-kT*log(temp)
           vtemp=sqrt(2.0*rho)*cos(theta)+vbar
        enddo
        return
        END SUBROUTINE distgen_extra

