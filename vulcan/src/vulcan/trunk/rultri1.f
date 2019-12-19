      SUBROUTINE RULTRI1(NDIME,NGAUS,POSGP,WEIGP)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP THE RADAU INTEGRATION CONSTANTS FOR
C     TRIANGLES AND TETRAHEDRA (OTHER OPTION)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION POSGP(NDIME,*), WEIGP(*)
C
C**** AREA INTEGRAL (TRIANGLES)
C
      IF(NDIME.EQ.2) THEN
       IF(NGAUS.EQ.1) THEN
c       call runendt('rultri1: not implemented')
        POSGP(1,1)= 1.D+00/3.D+00
        POSGP(2,1)= 1.D+00/3.D+00
C
        WEIGP(1)  = 1.D+00/2.D+00
        RETURN
       ENDIF
       IF(NGAUS.EQ.3) THEN
        POSGP(1,1)= 0.0D0
        POSGP(2,1)= 0.0D0
        POSGP(1,2)= 1.0D0
        POSGP(2,2)= 0.0D0
        POSGP(1,3)= 0.0D0
        POSGP(2,3)= 1.0D0
C
        WEIGP(1)  = 1.D+00/6.D+00
        WEIGP(2)  = 1.D+00/6.D+00
        WEIGP(3)  = 1.D+00/6.D+00
        RETURN
       ENDIF
       IF(NGAUS.EQ.4) THEN
c       call runendt('rultri1: not implemented')
        POSGP(1,1)= 1.D+00/3.D+00
        POSGP(2,1)= 1.D+00/3.D+00
        POSGP(1,2)= 1.D+00/5.D+00
        POSGP(2,2)= 1.D+00/5.D+00
        POSGP(1,3)= 3.D+00/5.D+00
        POSGP(2,3)= 1.D+00/5.D+00
        POSGP(1,4)= 1.D+00/5.D+00
        POSGP(2,4)= 3.D+00/5.D+00
C
        WEIGP(1)  =-27.D+00/96.D+00
        WEIGP(2)  = 25.D+00/96.D+00
        WEIGP(3)  = 25.D+00/96.D+00
        WEIGP(4)  = 25.D+00/96.D+00 
        RETURN 
       ENDIF
       IF(NGAUS.EQ.6) THEN
c       call runendt('rultri1: not implemented')
        EX1 = 0.81684 75729 80459
        ET1 = 0.09157 62135 09771
        EZ1 = 0.09157 62135 09771
        EX2 = 0.10810 30181 68070
        ET2 = 0.44594 84909 15965
        EZ2 = 0.44595 84909 15965
        POSGP(1,1)= EX1
        POSGP(2,1)= ET1
        POSGP(1,2)= ET1
        POSGP(2,2)= EZ1
        POSGP(1,3)= EZ1
        POSGP(2,3)= EX1
        POSGP(1,4)= EX2
        POSGP(2,4)= ET2
        POSGP(1,5)= ET2
        POSGP(2,5)= EZ2
        POSGP(1,6)= EZ2
        POSGP(2,6)= EX2
C
        A = 0.10995 17436 55322/2.0D00
        B = 0.22338 15896 78011/2.0D00
        WEIGP(1)  = A
        WEIGP(2)  = A
        WEIGP(3)  = A
        WEIGP(4)  = B
        WEIGP(5)  = B
        WEIGP(6)  = B
        RETURN 
       ENDIF
       IF(NGAUS.EQ.7) THEN
c       call runendt('rultri1: not implemented')
        EX1 = 0.33333 33333 33333
        ET1 = 0.33333 33333 33333
        EZ1 = 0.33333 33333 33333
        EX2 = 0.79742 69853 53087
        ET2 = 0.10128 65073 23456
        EZ2 = 0.10128 65073 23456
        EX3 = 0.47014 20641 05115
        ET3 = 0.47014 20641 05115
        EZ3 = 0.05971 58717 89770
        POSGP(1,1)= EX2
        POSGP(2,1)= ET2
        POSGP(1,2)= ET2
        POSGP(2,2)= EZ2
        POSGP(1,3)= EZ2
        POSGP(2,3)= EX2
        POSGP(1,4)= EX3
        POSGP(2,4)= ET3
        POSGP(1,5)= ET3
        POSGP(2,5)= EZ3
        POSGP(1,6)= EZ3
        POSGP(2,6)= EX3
        POSGP(1,6)= EX1
        POSGP(2,6)= ET1
C
        A = 0.12593 91805 44827/2.0D00
        B = 0.13239 41527 88506/2.0D00
        C = 0.22500 00000 00000/2.0D00
        WEIGP(1)  = A
        WEIGP(2)  = A
        WEIGP(3)  = A
        WEIGP(4)  = B
        WEIGP(5)  = B
        WEIGP(6)  = B
        WEIGP(7)  = C
        RETURN 
       ENDIF
C
      ENDIF
C
C**** VOLUME INTEGRAL (TETRAHEDRA)
C
      IF(NDIME.EQ.3) THEN
       IF(NGAUS.EQ.1) THEN
        POSGP(1,1)= 1.D+00/4.D+00
        POSGP(2,1)= 1.D+00/4.D+00
        POSGP(3,1)= 1.D+00/4.D+00
C
        WEIGP(1)  = 1.D+00/6.D+00
        RETURN
       ENDIF
C
       IF(NGAUS.EQ.4) THEN
        A=0.58541020D+00
        B=0.13819660D+00

c       A=0.333333333D0    ! to compute B-bar-vol for 5-noded tetrahedra
c       B=0.166666666D0

c       A=0.000000000D0
c       B=0.333333333D0

c       A=0.60541020D+00
c       B=0.110D+00
C
        POSGP(1,1)= B
        POSGP(2,1)= B
        POSGP(3,1)= B
        POSGP(1,2)= A
        POSGP(2,2)= B
        POSGP(3,2)= B
        POSGP(1,3)= B
        POSGP(2,3)= A
        POSGP(3,3)= B
        POSGP(1,4)= B
        POSGP(2,4)= B
        POSGP(3,4)= A
C
c       POSGP(1,1)= 0.0D0
c       POSGP(2,1)= 0.0D0
c       POSGP(3,1)= 0.0D0
c       POSGP(1,2)= 1.0D0
c       POSGP(2,2)= 0.0D0
c       POSGP(3,2)= 0.0D0
c       POSGP(1,3)= 0.0D0
c       POSGP(2,3)= 1.0D0
c       POSGP(3,3)= 0.0D0
c       POSGP(1,4)= 0.0D0
c       POSGP(2,4)= 0.0D0
c       POSGP(3,4)= 1.0D0
C
        WEIGP(1)  = 1.D+00/24.D+00
        WEIGP(2)  = 1.D+00/24.D+00
        WEIGP(3)  = 1.D+00/24.D+00
        WEIGP(4)  = 1.D+00/24.D+00
        RETURN
       ENDIF
       IF(NGAUS.EQ.5) THEN
c       call runendt('rultri1: not implemented')
        POSGP(1,1)= 1.D+00/4.D+00
        POSGP(2,1)= 1.D+00/4.D+00
        POSGP(3,1)= 1.D+00/4.D+00
        POSGP(1,2)= 1.D+00/6.D+00
        POSGP(2,2)= 1.D+00/6.D+00
        POSGP(3,2)= 1.D+00/6.D+00
        POSGP(1,3)= 1.D+00/3.D+00
        POSGP(2,3)= 1.D+00/6.D+00
        POSGP(3,3)= 1.D+00/6.D+00
        POSGP(1,4)= 1.D+00/6.D+00
        POSGP(2,4)= 1.D+00/3.D+00
        POSGP(3,4)= 1.D+00/6.D+00
        POSGP(1,5)= 1.D+00/6.D+00
        POSGP(2,5)= 1.D+00/6.D+00
        POSGP(3,5)= 1.D+00/3.D+00
C
        WEIGP(1)  =-4.D+00/30.D+00
        WEIGP(2)  = 3.D+00/40.D+00
        WEIGP(3)  = 3.D+00/40.D+00
        WEIGP(4)  = 3.D+00/40.D+00
        WEIGP(5)  = 3.D+00/40.D+00
        RETURN
       ENDIF
      ENDIF
C
      END
