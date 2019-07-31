
program test
#ifdef DOUBLE
 double precision :: value=0.0D0
#else
 real :: value=0.0E0
#endif
write(*,*) value
end program test
