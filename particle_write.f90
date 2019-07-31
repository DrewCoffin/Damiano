	SUBROUTINE particle_write(size2,rank,root)
	use indeces
	use mesh
        use param
        use parparam
	implicit none
        integer :: size2,root,ii,rank,k,offset,jj
        character(len=11) fname
        call getname_mpi('d_pa',index,rank,fname)
        open(10+rank,file=fname,form='unformatted',status='unknown')
        write(10+rank) mu2_p,nu_p,vpar2,mu_m,iprecip
        close(10+rank)
        return
        END SUBROUTINE particle_write

