         SUBROUTINE input_out()
         use indeces
         use param
         use parparam
         implicit none

! write out all paramters to input.out
!
         open(unit=3,file='input.out',status='unknown')

         write(3,1) 'Program GKE - version:',version
         write(3,*) '-----------------------'
         write(3,*)
         write(3,10) 'tmax=       :',tmax
         write(3,10) 'mphi=       :',mphi
         write(3,10) 'Bn=         :',Bn
         write(3,10) 'L=          :',L
         write(3,10) 'tn=         :',tn
         write(3,20) 'np=         :',np
         write(3,10) 'r_ion=      :',r_ionosphere
         write(3,10) 'ro1=        :',ro1
         write(3,10) 'ro2=        :',ro2
         write(3,10) 'twrite=     :',twrite
         write(3,10) 'dtwrite=    :',dtwrite
         write(3,20) 'nt=         :',nt
         write(3,10) 'twrite_fw=  :',twrite_fw
         write(3,10) 'dtwrite_fw= :',dtwrite_fw
         write(3,20) 'nt_fw=      :',nt_fw
         write(3,10) 'dt=         :',dt
         write(3,20) 'it=         :',it
         write(3,20) 'n1=         :',n1
         write(3,20) 'n2=         :',n2
         write(3,20) 'n1p=        :',n1p
         write(3,20) 'n2p=        :',n2p
! convert back to eV (implied factor of 2 still as before)
         write(3,10) 'T_e=         :',T_e/1.609e-19
         write(3,20) 'nproc=      :',nproc
         
         write(3,20)
         write(3,20) 'Switches'
         write(3,20) '--------'
         write(3,20)
         write(3,20) 'i_inertia=  :',i_inertia
         write(3,20) 'iparticle=  :',iparticle
         write(3,20) 'iprofile =  :',iprofile
         write(3,20) 'icouple  =  :',icouple
         write(3,20) 'igyro    =  :',igyro
         write(3,*)


         write(3,*)
         write(3,*) 'Normalizations'
         write(3,*) '--------------'
         write(3,*) 
         write(3,10) 'vn=         :',vn
         write(3,10) 'Bn=         :',Bn
         write(3,10) 'rho_n=      :',rho_n
         write(3,10) 'jn=         :',jn
         write(3,10) 'tn=         :',tn
         write(3,10) 'en=         :',en
         write(3,10) 'mu_mn=      :',mu_mn
         write(3,10) 'constE_mu=  :',constE_mu
         write(3,*)

         close(3)
 1       format(a22,1x,f6.1)
 10      format(a12,1x,e10.4)
 20      format(a10,3x,i10)
         return
         END SUBROUTINE input_out
         
