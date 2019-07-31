        SUBROUTINE full_read(size2,rank,root,i)
        use indeces
        use param
        use fields
        use parparam
        implicit none
        integer :: rank,size2,root,i
        character(len=7) :: fname

        call getname('d_fs',rank,fname)

        open(12+rank,file=fname,form='unformatted',status='unknown')
        read(12+rank) i,t,index,nt,index_df,nt_df,index_fw,nt_fw,nt_p
        read(12+rank) ltheta,iprecip_t
        read(12+rank) twrite,dtwrite,twrite_fw,dtwrite_fw
        read(12+rank) twrite_p,dtwrite_p
        read(12+rank) ii_vpamax,vmin,vmax,index_p,ipwrite
        read(12+rank) y,y0
        read(12+rank) j1_0,j1_1,j1_2,j2,j3
        read(12+rank) e2,e3,Fc0
        read(12+rank) e1_0,e1_1,e1_2
        read(12+rank) vpar0,vpar1,vpar2
        read(12+rank) ro_p,r_p
        read(12+rank) theta_p
        read(12+rank) mu0_p,mu1_p,mu2_p
        read(12+rank) mu_m,pa,nu_p,iprecip
        read(12+rank) h1p_1
        read(12+rank) dB1
        read(12+rank) e1p_1
        read(12+rank) e1_sum,e2_sum,y_sum
        close(12+rank)
        return
        END SUBROUTINE full_read
