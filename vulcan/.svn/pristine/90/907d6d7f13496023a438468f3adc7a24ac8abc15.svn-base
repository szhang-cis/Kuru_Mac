      SUBROUTINE RULTRI(NDIME,NGAUS,POSGP,WEIGP)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP THE RADAU INTEGRATION CONSTANTS FOR
C     TRIANGLES AND TETRAHEDRA
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
        POSGP(1,1)= 1.D+00/3.D+00
        POSGP(2,1)= 1.D+00/3.D+00
C
        WEIGP(1)  = 1.D+00/2.D+00
        RETURN
       ENDIF
       IF(NGAUS.EQ.3) THEN
        POSGP(1,1)= 1.D+00/6.D+00
        POSGP(2,1)= 1.D+00/6.D+00
        POSGP(1,2)= 2.D+00/3.D+00
        POSGP(2,2)= 1.D+00/6.D+00
        POSGP(1,3)= 1.D+00/6.D+00
        POSGP(2,3)= 2.D+00/3.D+00
C
        WEIGP(1)  = 1.D+00/6.D+00
        WEIGP(2)  = 1.D+00/6.D+00
        WEIGP(3)  = 1.D+00/6.D+00
        RETURN
       ENDIF
       IF(NGAUS.EQ.4) THEN
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
        EX1 = 0.81684 75729 80459D0
        ET1 = 0.09157 62135 09771D0
        EZ1 = 0.09157 62135 09771D0
        EX2 = 0.10810 30181 68070D0
        ET2 = 0.44594 84909 15965D0
        EZ2 = 0.44595 84909 15965D0
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
        A = 0.10995 17436 55322D0/2.0D00
        B = 0.22338 15896 78011D0/2.0D00
        WEIGP(1)  = A
        WEIGP(2)  = A
        WEIGP(3)  = A
        WEIGP(4)  = B
        WEIGP(5)  = B
        WEIGP(6)  = B
        RETURN 
       ENDIF
       IF(NGAUS.EQ.7) THEN
        EX1 = 0.33333 33333 33333D0 ! original: coord. seem to be wrong!
        ET1 = 0.33333 33333 33333D0
        EZ1 = 0.33333 33333 33333D0
        EX2 = 0.79742 69853 53087D0
        ET2 = 0.10128 65073 23456D0
        EZ2 = 0.10128 65073 23456D0
        EX3 = 0.47014 20641 05115D0
        ET3 = 0.47014 20641 05115D0
        EZ3 = 0.05971 58717 89770D0
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
c                                 ! Bathe (page 280) & Hughes (page 173)
        R1 = 0.10128 65073 23456D0
        R2 = 0.79742 69853 53087D0
        R3 = R1
        R4 = 0.47014 20641 05115D0
        R5 = R4
        R6 = 0.05971 58717 89770D0
        R7 = 0.33333 33333 33333D0

        S1 = R1
        S2 = R1
        S3 = R2
        S4 = R6
        S5 = R4
        S6 = R4
        S7 = R7

        POSGP(1,1)= R1
        POSGP(2,1)= S1
        POSGP(1,2)= R2
        POSGP(2,2)= S2
        POSGP(1,3)= R3
        POSGP(2,3)= S3
        POSGP(1,4)= R4
        POSGP(2,4)= S4
        POSGP(1,5)= R5
        POSGP(2,5)= S5
        POSGP(1,6)= R6
        POSGP(2,6)= S6
        POSGP(1,7)= R7
        POSGP(2,7)= S7
C
        A = 0.12593 91805 44827D0/2.0D00
        B = 0.13239 41527 88506D0/2.0D00
        C = 0.22500 00000 00000D0/2.0D00
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
       IF(NGAUS.EQ.13) THEN       ! Bathe (page 280) & Hughes (page 174)
        R1 = 0.06513 01029 02216D0
        R2 = 0.86973 97941 95568D0
        R3 = R1
        R4 = 0.31286 54960 04875D0
        R5 = 0.63844 41885 69809D0
        R6 = 0.04869 03154 25316D0
        R7 = R5
        R8 = R4
        R9 = R6
        R10= 0.26034 59660 79038D0
        R11= 0.47930 80678 41923D0
        R12= R10
        R13= 0.33333 33333 33333D0

        S1 = R1
        S2 = R1
        S3 = R2
        S4 = R6
        S5 = R4
        S6 = R5
        S7 = R6
        S8 = R5
        S9 = R4
        S10= R10
        S11= R10
        S12= R11
        S13= R13

        POSGP(1, 1)= R1
        POSGP(2, 1)= S1
        POSGP(1, 2)= R2
        POSGP(2, 2)= S2
        POSGP(1, 3)= R3
        POSGP(2, 3)= S3
        POSGP(1, 4)= R4
        POSGP(2, 4)= S4
        POSGP(1, 5)= R5
        POSGP(2, 5)= S5
        POSGP(1, 6)= R6
        POSGP(2, 6)= S6
        POSGP(1, 7)= R7
        POSGP(2, 7)= S7
        POSGP(1, 8)= R8
        POSGP(2, 8)= S8
        POSGP(1, 9)= R9
        POSGP(2, 9)= S9
        POSGP(1,10)= R10
        POSGP(2,10)= S10
        POSGP(1,11)= R11
        POSGP(2,11)= S11
        POSGP(1,12)= R12
        POSGP(2,12)= S12
        POSGP(1,13)= R13
        POSGP(2,13)= S13

        A = 0.05334 72356 08839D0/2.0D00
        B = 0.07711 37608 90257D0/2.0D00
        C = 0.17561 52574 33204D0/2.0D00
        D =-0.14957 00444 67670D0/2.0D00
        WEIGP( 1)  = A
        WEIGP( 2)  = A
        WEIGP( 3)  = A
        WEIGP( 4)  = B
        WEIGP( 5)  = B
        WEIGP( 6)  = B
        WEIGP( 7)  = B
        WEIGP( 8)  = B
        WEIGP( 9)  = B
        WEIGP(10)  = C
        WEIGP(11)  = C
        WEIGP(12)  = C
        WEIGP(13)  = D
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
        WEIGP(1)  = 1.D+00/24.D+00
        WEIGP(2)  = 1.D+00/24.D+00
        WEIGP(3)  = 1.D+00/24.D+00
        WEIGP(4)  = 1.D+00/24.D+00
        RETURN
       ENDIF
       IF(NGAUS.EQ.5) THEN
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
