      SUBROUTINE QPARTT(PROPS,EFFPL,
     .                  QPART,FRACTF,FRACTA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PARTITION COEFFICIENT (EFFECTIVE
C     PLASTIC STRAIN DEPENDENT)
C     
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROPS(*)
C
      IPEP4=INT(PROPS(4))
C
      FRACTF=FRACT(1)
      FRACTA=FRACT(2)
C
      IPEP3=INT(PROPS(3))
      IF(IPEP3.EQ.1) THEN
       QPART=VPART(1,1)
       RETURN
      ENDIF
C
      IF(IPEP4.EQ.3) THEN                  ! non-standard: DP
       IF(EFFPL.LE.VPART(1,2)) QPART=VPART(1,1)
       DO IPART=2,NPART
        I1=IPART-1
        I2=IPART
        IF(EFFPL.GE.VPART(I1,2).AND.EFFPL.LE.VPART(I2,2))
     .   QPART=(VPART(I2,1)-VPART(I1,1))/(VPART(I2,2)-VPART(I1,2))*
     .         (EFFPL-VPART(I1,2))+VPART(I1,1)
       ENDDO
       IF(EFFPL.GE.VPART(NPART,2)) QPART=VPART(NPART,1)
      ENDIF               ! ipep4.eq.3
C
      RETURN
      END
