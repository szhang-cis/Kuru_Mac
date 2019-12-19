      SUBROUTINE TRASTR(S,T,NS,ND,FLAG)
C******************************************************************
C
C**** ROUTINE TO TRANSFORM STRAIN TENSOR
C
C******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S(NS), T(ND,ND)
      DIMENSION C1(6), C2(6)
      CHARACTER*6 FLAG
C
      DATA C1/1.0D00,1.0D00,2.0D00,1.0D00,2.0D00,2.0D00/
      DATA C2/1.0D00,1.0D00,0.5D00,1.0D00,0.5D00,0.5D00/
C
      DO I=1,NS
        S(I)=S(I)*C2(I)
      ENDDO
C
      CALL TRASIG(S,T,NS,ND,FLAG)
C
      DO I=1,NS
        S(I)=S(I)*C1(I)
      ENDDO
C
      RETURN
      END
