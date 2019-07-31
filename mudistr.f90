        SUBROUTINE mudistr(size2,rank,root)
        use indeces
        use param
        use parparam
        implicit none
        integer :: ii,jj,ic(500),size2,rank,root
        double precision ::dv,mubin(500)
        character(len=11) :: fname
        call getname_mpi('d_mu',index_df,rank,fname)
        dv=(mu_max-mu_min)/float(nbin)
        do ii=1,nbin
           mubin(ii)=0.
           ic(ii)=0
        enddo
        mubin(1)=mu_min
        do ii=2,nbin+1
            mubin(ii)=mu_min+dv*float(ii-1)
        enddo

        do ii=1,np2
           if(ro_p(ii).gt.ro_d_min.and. &
              ro_p(ii).lt.ro_d_max)then
              do jj=1,nbin
                 if((mu_m(ii).gt.mubin(jj)).and. &
                     (mu_m(ii).le.mubin(jj+1)))then
                     ic(jj)=ic(jj)+1
                 endif
              enddo
           endif
        enddo
 
        open(unit=9,file=fname,status='unknown')
        do  ii=1,nbin
           write(9,*) mubin(ii),ic(ii)
           write(9,*) mubin(ii+1),ic(ii)
        enddo
        close(9)

        return
        END SUBROUTINE mudistr

