        SUBROUTINE padistr(size2,rank,root)
        use indeces
        use param
        use parparam
        implicit none
        integer :: ii,jj,kk,ic(500),ii_temp(10000)
        integer :: size2,rank,root
        double precision :: dv,pabin(500),vpamax
        character(len=11) :: fname
        call getname_mpi('d_ao',index_df,rank,fname)

        dv=(pa_max-pa_min)/float(nbin)
        do ii=1,nbin
           pabin(ii)=0.
           ic(ii)=0
        enddo
        pabin(1)=pa_min
        do ii=2,nbin+1
            pabin(ii)=pa_min+dv*float(ii-1)
        enddo
        kk=0
        vpamax=0.d0
        do ii=1,np2
           do jj=1,nbin
               if((pa(ii).gt.pabin(jj)).and. &
                  (pa(ii).le.pabin(jj+1)))then
                  ic(jj)=ic(jj)+1
               endif
           enddo
        enddo
 
        open(unit=9,file=fname,status='unknown')
        do  ii=1,nbin
           write(9,*) pabin(ii),ic(ii)
           write(9,*) pabin(ii+1),ic(ii)
        enddo
        close(9)

        return
        END SUBROUTINE padistr
