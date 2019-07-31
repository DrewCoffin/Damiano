!
!**********************************************************************
!      subroutine  distgen_rick(vt,no,vth)
       SUBROUTINE distgen_vz2(no,vth)
!**********************************************************************
!..................R. Sydora (1/19/90)
      use indeces
      use distgen
      implicit double precision(a-h,o-z)
!      dimension xtemp(no),vpar(no)
      integer :: near2p,jj,jnext
      double precision :: powern,fract
!      real powern,fract
!**********************************************************************
!
!-----------
! Initialize base-2 bit-reversed number generation stuff (1/19/90):
!
!
! Find nearest power of 2 greater than no (assume no < 2**25):
      idone = 0
      do 125 i=0,25,1
        if (((2**i).ge.no) .and. (idone.eq.0)) then
          ipower = i
          idone = 1
        endif
 125  continue
      near2p = 2**ipower
! Counter for bit-reversed numbers:
      icount = 0
!
!-----------

!  1D parallel velocity distribution

      nv2 = 0
      vmax = 5.*vth
      dv = 2.*vmax/(no-1)
      vvnv2 = 1.
!   use xtemp as a dummy array for  cum. distribution function
      xtemp(1) = 0.
!   store the cumulative distribution in xtemp
      do  121 i = 2,no
      vv = ((i-1.5)*dv-vmax)/vth
      if(nv2.ne.0) vvnv2 = vv**nv2
      fv = vvnv2*exp(-.5*vv**2)
      xtemp(i) = xtemp(i-1) + max(fv,0.)
 121  continue
!
!
      df = xtemp(no)/no
      j = 1

      do 122 i = 1,no

         fv = (i-.5) * df
 123     if(fv.lt.xtemp(j+1)) go to 124
         j = j+1
         if (j .gt. (no-1)) then
           write(7,*) 'At i = ',i,' stop loop 122 in INTSTS; j = ',j
!           stop
            goto 199
         endif
         go to 123
 124     vv = dv*(j-1+(fv-xtemp(j))/(xtemp(j+1)-xtemp(j))) -vmax
!
!
!...find particle index from base-2 bit rev. numbers
!
 93   icount = icount + 1
      jj = icount
      fract = 0.0
      powern = 1.0
 91   jnext = jj/2
      powern = powern/2
      fract = fract + (jj - jnext*2)*powern
      if (jnext .eq. 0) go to 92
      jj = jnext
      go to 91
 92   iindex = fract*near2p + 1.0 + 1.0e-6
      if (iindex .gt. no) go to 93
! Avoid double hits:
      if (vz(iindex) .ne. 0.0) go to 93
!
      vz(iindex) =  vv
!
!
 122  continue
!
!
 199  return
      END SUBROUTINE distgen_vz2
