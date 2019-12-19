      SUBROUTINE ANGEUL( ROMTX, ANGLS, IFLAG )
C***********************************************************************
C
C     TRANSFORMATION BETWEEN EULER'S ANGLES AND ROTATION MATRIX
C
C     IFLAG=0 : ROTATION MATRIX FROM EULER'S ANGLES
C     IFLAG=1 : EULER'S ANGLES FROM ROTATION MATRIX
C
C              XG = R XN  : ANG GLOB -> NODL
C
C     ANGLS(1) : ROTATION ABOUT Z AXE
C     ANGLS(2) : ROTATION ABOUT X AXE
C     ANGLS(3) : ROTATION ABOUT Z AXE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER ( CUADR=90.0, SECIR=180.0, CUZER=1.0D-12 )
C
      DIMENSION ROMTX(3,3), ANGLS(3)
C                                                  ! ANGLES IN DEGREES
      IF (IFLAG.EQ.0) THEN     ! FROM EULER'S ANGLES TO ROTATION MATRIX
	 COSDA=DCOS(ANGLS(1))
	 SINDA=DSIN(ANGLS(1))
	 COSDB=DCOS(ANGLS(2))
	 SINDB=DSIN(ANGLS(2))
	 COSDG=DCOS(ANGLS(3))
	 SINDG=DSIN(ANGLS(3))
	 ROMTX(1,1)= COSDA*COSDG-SINDA*COSDB*SINDG
	 ROMTX(2,1)= SINDA*COSDG+COSDA*COSDB*SINDG
	 ROMTX(3,1)= SINDB*SINDG
	 ROMTX(1,2)=-COSDA*SINDG-SINDA*COSDB*COSDG
	 ROMTX(2,2)=-SINDA*SINDG+COSDA*COSDB*COSDG
	 ROMTX(3,2)= SINDB*COSDG
	 ROMTX(1,3)= SINDB*SINDA
	 ROMTX(2,3)=-SINDB*COSDA
	 ROMTX(3,3)= COSDB
      ELSE IF (IFLAG.EQ.1) THEN ! FROM ROTATION MATRIX TO EULER'S ANGLES
         IF ( ROMTX(3,3).GE.1.0 ) THEN
            ANGLS(2) = 0.0
         ELSEIF ( ROMTX(3,3).LE.-1.0) THEN
            ANGLS(2) = 180.0
         ELSE
   	    ANGLS(2)=DACOS( ROMTX(3,3) )
         ENDIF
         IF (ANGLS(2).EQ.0.0) THEN
            ANGLS(1)=0.0
         ELSE
            IF ( ABS(ROMTX(2,3)).LT.CUZER ) THEN
               IF (ROMTX(1,3).GT.0.0) THEN
                  ANGLS(1)= CUADR
               ELSE
                  ANGLS(1)=-CUADR
               ENDIF
            ELSE
               ANGLS(1)=DATAN( -ROMTX(1,3)/ROMTX(2,3) )
               IF (ROMTX(2,3).GT.0) ANGLS(1)=ANGLS(1)+SECIR
C                                    ! COSINE OF BETA NEGATIVE
            ENDIF
         END IF
         COSDA=DCOS(ANGLS(1))
         SINDA=DSIN(ANGLS(1))
         COSDG= ROMTX(1,1)*COSDA + ROMTX(2,1)*SINDA
         SINDG=-ROMTX(1,2)*COSDA - ROMTX(2,2)*SINDA
         IF (ABS(COSDG).LT.CUZER) THEN
            IF (SINDG.GT.0.0) THEN
               ANGLS(3)=CUADR
            ELSE
               ANGLS(3)=-CUADR
            ENDIF
         ELSE
            ANGLS(3)=DATAN( SINDG/COSDG )
            IF (COSDG.LT.0.0) ANGLS(3)=ANGLS(3)+SECIR
C                                    ! COSINE NEGATIVE
         ENDIF
      ENDIF
      RETURN
      END
