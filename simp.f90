!***************************************************************************
! simpson's rule integration routine                                       *
!***************************************************************************

      SUBROUTINE SIMP(Y,DX,INUM,VAL)
!     Simpson's rule routine
      IMPLICIT NONE
      DOUBLE PRECISION :: SUM,VAL,Y,DX,T
      DIMENSION Y(*)
      INTEGER :: I,INUM
      SUM=0.0
      VAL=0.0
      DO 50 I=1,INUM
         IF ((I.EQ.1).OR.(I.EQ.INUM)) THEN
            T=1.0
            GOTO 40
    
         ENDIF
         T=3.0 + (-1.0)**I
  40     SUM = SUM + T*Y(I)
  50  CONTINUE
      VAL = DX/3.0*SUM
      RETURN
      END SUBROUTINE SIMP
      

