      SUBROUTINE RULIUS(NGAUS,NDIME,POSGP,WEIGP)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP THE USER'S RULES INTEGRATION CONSTANTS
C
C
C     Notes:
C
C     The use of NRULE/NRULET=5 (integration points at nodes) in the
C     smoothing operation (with ISMO1/ISMO1T=1; see elm030.f and
C     elm005t.f) is not absolutely correct when NRULE/NRULET=1 (Gauss
C     integration points) is given in the mechanical/thermal input data
C     file because the location of nodes and Gauss points is different.
C
C     The use of NGAUS=1 and NRULE=5 for the standard computations is
C     not possible (see chek01.f, chek01t.f & chek01s.f).
C
C     For the smoothing operations (with ISMO1/ISMO1T=1; see elm030.f
C     and elm005t.f), the integration points and weights with NGAUS=1
C     are not used.
C
C     This routine must be improved! (March/2000) i.e.: NGAUS=1, etc
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION POSGP(NDIME,*), WEIGP(*)
C
C**** ONE DIMENSIONAL INTEGRAL
C
      IF(NDIME.EQ.1) THEN
       IF(NGAUS.EQ.2) THEN
        POSGP(1,1)=-1.0D0
	POSGP(1,2)= 1.0D0
C
        WEIGP(1)=1.0D0
	WEIGP(2)=1.0D0
       ENDIF
C
      ENDIF
C
C**** AREA INTEGRAL
C
      IF(NDIME.EQ.2) THEN
C
       IF(NGAUS.EQ.1) THEN     ! idem rultri.f (to avoid smooth. prob.)
        POSGP(1,1)= 1.D+00/3.D+00
        POSGP(2,1)= 1.D+00/3.D+00
C
        WEIGP(1)  = 1.D+00/2.D+00
        RETURN
       ENDIF
C
C**** TRIANGLES
C
       IF(NGAUS.EQ.3) THEN     ! Lobatto (almost comp. with rultri.f)
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
       ENDIF                   ! ngaus.eq.3
C
C**** QUADRILATERALS
C
       IF(NGAUS.EQ.4) THEN     ! Lobatto (compatible with rullob.f
C                                & almost compatible with rulgau.f)
        POSGP(1,1)=-1.D+00
        POSGP(2,1)=-1.D+00
        POSGP(1,2)=-1.D+00
        POSGP(2,2)= 1.D+00
        POSGP(1,3)= 1.D+00
        POSGP(2,3)=-1.D+00
        POSGP(1,4)= 1.D+00
        POSGP(2,4)= 1.D+00
C                              ! (compatible with LNODS numbering)
c       POSGP(1,1)=-1.D+00
c       POSGP(2,1)=-1.D+00
c       POSGP(1,2)= 1.D+00
c       POSGP(2,2)=-1.D+00
c       POSGP(1,3)= 1.D+00
c       POSGP(2,3)= 1.D+00
c       POSGP(1,4)=-1.D+00
c       POSGP(2,4)= 1.D+00
C
        WEIGP(1)  = 1.D+00
        WEIGP(2)  = 1.D+00
        WEIGP(3)  = 1.D+00
	WEIGP(4)  = 1.D+00
        RETURN
       ENDIF
C
       IF(NGAUS.EQ.9) THEN     ! Lobatto (compatible with rullob.f)
        POSGP(1,1)=-1.D+00
        POSGP(2,1)=-1.D+00
        POSGP(1,2)=-1.D+00
        POSGP(2,2)= 0.D+00
        POSGP(1,3)=-1.D+00
        POSGP(2,3)= 1.D+00
        POSGP(1,4)= 0.D+00
        POSGP(2,4)=-1.D+00
        POSGP(1,5)= 0.D+00
        POSGP(2,5)= 0.D+00
        POSGP(1,6)= 0.D+00
        POSGP(2,6)= 1.D+00
        POSGP(1,7)= 1.D+00
        POSGP(2,7)=-1.D+00
        POSGP(1,8)= 1.D+00
        POSGP(2,8)= 0.D+00
        POSGP(1,9)= 1.D+00
        POSGP(2,9)= 1.D+00
C
        WEIGP(1)  = 1.D+00/9.D+00
        WEIGP(2)  = 4.D+00/9.D+00
        WEIGP(3)  = 1.D+00/9.D+00
        WEIGP(4)  = 4.D+00/9.D+00
        WEIGP(5)  =16.D+00/9.D+00
        WEIGP(6)  = 4.D+00/9.D+00
        WEIGP(7)  = 1.D+00/9.D+00
        WEIGP(8)  = 4.D+00/9.D+00
        WEIGP(9)  = 1.D+00/9.D+00
        RETURN
       ENDIF
C
      ENDIF                    ! ndime.eq.2
C
C**** VOLUME INTEGRAL
C
      IF(NDIME.EQ.3) THEN
C
C**** TOBLERONES
C
       IF(NGAUS.EQ.2) THEN     ! Lobatto type in t direction
        POSGP(1,1)= 1.D+00/3.D+00
        POSGP(2,1)= 1.D+00/3.D+00
        POSGP(3,1)=-1.D+00
        POSGP(1,2)= 1.D+00/3.D+00
        POSGP(2,2)= 1.D+00/3.D+00
        POSGP(3,2)= 1.D+00
C
        WEIGP(1)  = 1.D+00/2.D+00
        WEIGP(2)  = 1.D+00/2.D+00
        RETURN
       ENDIF                   ! ngaus.eq.2
C
C**** TETRAHEDRA
C
       IF(NGAUS.EQ.4) THEN     ! Lobatto (almost comp. with rultri.f)
        POSGP(1,1)= 0.0D0
        POSGP(2,1)= 0.0D0
        POSGP(3,1)= 0.0D0
        POSGP(1,2)= 1.0D0
        POSGP(2,2)= 0.0D0
        POSGP(3,2)= 0.0D0
        POSGP(1,3)= 0.0D0
        POSGP(2,3)= 1.0D0
        POSGP(3,3)= 0.0D0
        POSGP(1,4)= 0.0D0
        POSGP(2,4)= 0.0D0
        POSGP(3,4)= 1.0D0
C
        WEIGP(1)  = 1.D+00/24.D+00
        WEIGP(2)  = 1.D+00/24.D+00
        WEIGP(3)  = 1.D+00/24.D+00
        WEIGP(4)  = 1.D+00/24.D+00
        RETURN
       ENDIF
C
C**** TOBLERONES
C
       IF(NGAUS.EQ.6) THEN     ! Lobatto type
        POSGP(1,1)= 0.D+00
        POSGP(2,1)= 0.D+00
        POSGP(3,1)=-1.D+00
        POSGP(1,2)= 1.D+00
        POSGP(2,2)= 0.D+00
        POSGP(3,2)=-1.D+00
        POSGP(1,3)= 0.D+00
        POSGP(2,3)= 1.D+00
        POSGP(3,3)=-1.D+00
        POSGP(1,4)= 0.D+00
        POSGP(2,4)= 0.D+00
        POSGP(3,4)= 1.D+00
        POSGP(1,5)= 1.D+00
        POSGP(2,5)= 0.D+00
        POSGP(3,5)= 1.D+00
        POSGP(1,6)= 0.D+00
        POSGP(2,6)= 1.D+00
        POSGP(3,6)= 1.D+00
C
        WEIGP(1)  = 1.D+00/6.D+00
        WEIGP(2)  = 1.D+00/6.D+00
        WEIGP(3)  = 1.D+00/6.D+00
        WEIGP(4)  = 1.D+00/6.D+00
        WEIGP(5)  = 1.D+00/6.D+00
        WEIGP(6)  = 1.D+00/6.D+00
       ENDIF                   ! ngaus.eq.6
C
C**** BRICKS
C
       IF(NGAUS.EQ.8) THEN     ! Lobatto (compatible with rullob.f)
        POSGP(1,1)=-1.D+00
        POSGP(2,1)=-1.D+00
        POSGP(3,1)=-1.D+00
        POSGP(1,2)=-1.D+00
        POSGP(2,2)=-1.D+00
        POSGP(3,2)= 1.D+00
        POSGP(1,3)=-1.D+00
        POSGP(2,3)= 1.D+00
        POSGP(3,3)=-1.D+00
        POSGP(1,4)=-1.D+00
        POSGP(2,4)= 1.D+00
        POSGP(3,4)= 1.D+00
        POSGP(1,5)= 1.D+00
        POSGP(2,5)=-1.D+00
        POSGP(3,5)=-1.D+00
        POSGP(1,6)= 1.D+00
        POSGP(2,6)=-1.D+00
        POSGP(3,6)= 1.D+00
        POSGP(1,7)= 1.D+00
        POSGP(2,7)= 1.D+00
        POSGP(3,7)=-1.D+00
        POSGP(1,8)= 1.D+00
        POSGP(2,8)= 1.D+00
        POSGP(3,8)= 1.D+00
C
        WEIGP(1)  = 1.D+00
        WEIGP(2)  = 1.D+00
        WEIGP(3)  = 1.D+00
        WEIGP(4)  = 1.D+00
        WEIGP(5)  = 1.D+00
        WEIGP(6)  = 1.D+00
        WEIGP(7)  = 1.D+00
        WEIGP(8)  = 1.D+00
        RETURN
       ENDIF
C
       IF(NGAUS.EQ.27) THEN    ! Lobatto (compatible with rullob.f)
        POSGP(1,1)= -1.D+00
        POSGP(2,1)= -1.D+00
        POSGP(3,1)= -1.D+00
        POSGP(1,2)= -1.D+00
        POSGP(2,2)= -1.D+00
        POSGP(3,2)=  0.D+00
        POSGP(1,3)= -1.D+00
        POSGP(2,3)= -1.D+00
        POSGP(3,3)=  1.D+00
        POSGP(1,4)= -1.D+00
        POSGP(2,4)=  0.D+00
        POSGP(3,4)= -1.D+00
        POSGP(1,5)= -1.D+00
        POSGP(2,5)=  0.D+00
        POSGP(3,5)=  0.D+00
        POSGP(1,6)= -1.D+00
        POSGP(2,6)=  0.D+00
        POSGP(3,6)=  1.D+00
        POSGP(1,7)= -1.D+00
        POSGP(2,7)=  1.D+00
        POSGP(3,7)= -1.D+00
        POSGP(1,8)= -1.D+00
        POSGP(2,8)=  1.D+00
        POSGP(3,8)=  0.D+00
        POSGP(1,9)= -1.D+00
        POSGP(2,9)=  1.D+00
        POSGP(3,9)=  1.D+00
        POSGP(1,10)= 0.D+00
        POSGP(2,10)=-1.D+00
        POSGP(3,10)=-1.D+00
        POSGP(1,11)= 0.D+00
        POSGP(2,11)=-1.D+00
        POSGP(3,11)= 0.D+00
        POSGP(1,12)= 0.D+00
        POSGP(2,12)=-1.D+00
        POSGP(3,12)= 1.D+00
        POSGP(1,13)= 0.D+00
        POSGP(2,13)= 0.D+00
        POSGP(3,13)=-1.D+00
        POSGP(1,14)= 0.D+00
        POSGP(2,14)= 0.D+00
        POSGP(3,14)= 0.D+00
        POSGP(1,15)= 0.D+00
        POSGP(2,15)= 0.D+00
        POSGP(3,15)= 1.D+00
        POSGP(1,16)= 0.D+00
        POSGP(2,16)= 1.D+00
        POSGP(3,16)=-1.D+00
        POSGP(1,17)= 0.D+00
        POSGP(2,17)= 1.D+00
        POSGP(3,17)= 0.D+00
        POSGP(1,18)= 0.D+00
        POSGP(2,18)= 1.D+00
        POSGP(3,18)= 1.D+00
        POSGP(1,19)= 1.D+00
        POSGP(2,19)=-1.D+00
        POSGP(3,19)=-1.D+00
        POSGP(1,20)= 1.D+00
        POSGP(2,20)=-1.D+00
        POSGP(3,20)= 0.D+00
        POSGP(1,21)= 1.D+00
        POSGP(2,21)=-1.D+00
        POSGP(3,21)= 1.D+00
        POSGP(1,22)= 1.D+00
        POSGP(2,22)= 0.D+00
        POSGP(3,22)=-1.D+00
        POSGP(1,23)= 1.D+00
        POSGP(2,23)= 0.D+00
        POSGP(3,23)= 0.D+00
        POSGP(1,24)= 1.D+00
        POSGP(2,24)= 0.D+00
        POSGP(3,24)= 1.D+00
        POSGP(1,25)= 1.D+00
        POSGP(2,25)= 1.D+00
        POSGP(3,25)=-1.D+00
        POSGP(1,26)= 1.D+00
        POSGP(2,26)= 1.D+00
        POSGP(3,26)= 0.D+00
        POSGP(1,27)= 1.D+00
        POSGP(2,27)= 1.D+00
        POSGP(3,27)= 1.D+00
C
        WEIGP(1)  = 1.D+00/27.D+00
        WEIGP(2)  = 4.D+00/27.D+00
        WEIGP(3)  = 1.D+00/27.D+00
        WEIGP(4)  = 4.D+00/27.D+00
        WEIGP(5)  =16.D+00/27.D+00
        WEIGP(6)  = 4.D+00/27.D+00
        WEIGP(7)  = 1.D+00/27.D+00
        WEIGP(8)  = 4.D+00/27.D+00
        WEIGP(9)  = 1.D+00/27.D+00
        WEIGP(10) = 4.D+00/27.D+00
        WEIGP(11) =16.D+00/27.D+00
        WEIGP(12) = 4.D+00/27.D+00
        WEIGP(13) =16.D+00/27.D+00
        WEIGP(14) =64.D+00/27.D+00
        WEIGP(15) =16.D+00/27.D+00
        WEIGP(16) = 4.D+00/27.D+00
        WEIGP(17) =16.D+00/27.D+00
        WEIGP(18) = 4.D+00/27.D+00
        WEIGP(19) = 1.D+00/27.D+00
        WEIGP(20) = 4.D+00/27.D+00
        WEIGP(21) = 1.D+00/27.D+00
        WEIGP(22) = 4.D+00/27.D+00
        WEIGP(23) =16.D+00/27.D+00
        WEIGP(24) = 4.D+00/27.D+00
        WEIGP(25) = 1.D+00/27.D+00
        WEIGP(26) = 4.D+00/27.D+00
        WEIGP(27) = 1.D+00/27.D+00
        RETURN
       ENDIF
C
      ENDIF                    ! ndime.eq.3
C
      END
