        SUBROUTINE pa_init(size2,rank,root)
        use indeces
        use param
        use parparam
        implicit none
        integer :: ii,size2,rank,root,localseed
        integer :: localseed_vx,localseed_vy
        double precision :: randf,slope,A,mt,b,P
        external randf,P
        localseed=(1001*(rank+1))
        localseed_vx=(100001*(rank+1))
        localseed_vy=(100000001*(rank+1))
        if(ipa_profile.eq.1)then
! for uniform pitch angle distribution
          do ii=1,np2
             pa(ii)=pa_min+randf(seed)*(pa_max-pa_min)
          enddo
! for linear pitch angle distribution (problem)
        elseif(ipa_profile.eq.2)then
           slope=(P(pa_max)-P(pa_min))/ &
                  (pa_max-pa_min)
           b=P(pa_min-pa_min)
           A=(slope/2.d0)*pa_max**2+b*pa_max- &
             ((slope/2.d0)*pa_min**2+b*pa_min)
           do  ii=1,np2
               mt=-((slope/2.d0)*pa_min**2 + b*pa_min + A*randf(seed))
!                mt=-((slope/2.d0)*pa_min**2 + b*pa_min + A)
                pa(ii)= (-b+sqrt(b*b-2.d0*slope*mt))/slope
                write(*,*) pa(ii),pa_min,pa_max
           enddo
        elseif(ipa_profile.eq.3)then
           do ii=1,np2
!              pa(ii)=acos(1.d0-randf(seed))
              pa(ii)=acos(1.d0-randf(localseed))
           enddo
        endif
        return
        END SUBROUTINE pa_init             

! FIX function

        FUNCTION P(x)
        double precision :: xmin,x,P
        P=-0.75*(x)+2.d0
        return
        END FUNCTION P

