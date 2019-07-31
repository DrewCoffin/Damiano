         PROGRAM df_extract
         use indeces
         use df
         implicit none
         integer,parameter :: npointsi=nbini+1,npointsj=nbinj+1
         integer :: i,j,ii,jj,k1,index
         integer :: df_sum(n2,n1)
         integer :: ichoice, jchoice
         character(len=7) :: fname, fname2
         integer :: df_array2(nbini,nbinj)
         double precision :: bin(npointsi,npointsj)
         double precision :: vbin_perp(npointsi),vbin_par(npointsj)
         double precision :: vperp_min,vperp_max,vpar_min,vpar_max
         double precision :: dvpar
         double precision :: dvperp

         vperp_min=vp_min/vn
         vperp_max=vp_max/vn
         vpar_min=vp_min/vn
         vpar_max=vp_max/vn
!        vperp_min=-70.d0
!        vperp_max=70.d0
!        vpar_min=-70.d0
!        vpar_max=70.d0

         dvperp=(vperp_max-vperp_min)/float(nbini)
         dvpar=(vpar_max-vpar_min)/float(nbinj)


         do 500 index=1,27,1
         do ii=1,npointsi
            vbin_perp(ii)=vperp_min+dvperp*float(ii-1)
         enddo

         do jj=1,npointsj
            vbin_par(jj)=vpar_min+dvpar*float(jj-1)
         enddo

         call getname('d_df',index,fname)
         open(10,file=fname,form='unformatted',status='old')
         read(10) df_array
         close(10)

        do ii=1,nbini
           do jj=1,nbinj
              df_array2(ii,jj)=0
           enddo
        enddo
        ichoice=72
        ichoice=61
! l_||=0.5
        jchoice=304
! l_||=1.37
!       jchoice=388
! l_||=0.5
!       jchoice=352
         do ii=nbini/2+1,nbini
            do jj=1,nbinj
               do i=ichoice-3,ichoice+3
                 do j=jchoice-2,jchoice+2
                    df_array2(ii,jj)=df_array2(ii,jj)+ &
                                     df_array(ii,jj,i,j)
                 enddo
               enddo
            enddo
         enddo

         call getname('d_vv',index,fname2)
         open(12,file=fname2,status='unknown')
         do ii=nbini/2+1,nbini
            do jj=1,nbinj
                  write(12,*) 0.5*(vbin_perp(ii)+vbin_perp(ii+1)), &
                              0.5*(vbin_par(jj)+vbin_par(jj+1)),&
                              float(df_array2(ii,jj)) &
                           /(0.5*(vbin_perp(ii)+vbin_perp(ii+1)))
                  df_sum(i,j)=df_sum(i,j)+df_array(ii,jj,i,j)
            enddo
         enddo
         close(12)
 500     continue

         do i=1,n2
         do j=1,n1
         df_sum(i,j)=0
         do ii=1,nbini
            do jj=1,nbinj
               df_sum(i,j)=df_sum(i,j)+df_array(ii,jj,i,j)
            enddo
         enddo
         enddo
         enddo


         do i=1,n2
            do j=1,n1
               write(99,*) i,j,df_sum(i,j)
            enddo
         enddo


 55      format(5e13.5)
         stop
         END PROGRAM df_extract

