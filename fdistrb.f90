        SUBROUTINE fdistrb(vpar2,ro_p)
        use indeces
        use param
        implicit none
        integer :: ii,jj,ic(500)
        double precision :: dv,vbin(500),vpar2(*),ro_p(*)
        double precision :: vth
        character(len=7) :: fname
        call getname('d_df',index_df,fname)
        dv=(vmax-vmin)/float(nbin)
        do ii=1,nbin
           vbin(ii)=0.
           ic(ii)=0
        enddo
        vbin(1)=vmin
        do ii=2,nbin+1
            vbin(ii)=vmin+dv*float(ii-1)
        enddo

        do ii=1,np
           if(ro_p(ii).gt.ro_d_min.and. &
              ro_p(ii).lt.ro_d_max)then  
           do jj=1,nbin
               if((vpar2(ii).gt.vbin(jj)).and. &
                  (vpar2(ii).le.vbin(jj+1)))then
                  ic(jj)=ic(jj)+1
               endif
           enddo
           endif
        enddo

        open(unit=9,file=fname,status='unknown')
        do  ii=1,nbin
           write(9,*) vbin(ii),ic(ii)
           write(9,*) vbin(ii+1),ic(ii)
        enddo
        close(9)

        return
        END SUBROUTINE fdistrb

