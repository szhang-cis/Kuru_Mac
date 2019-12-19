      SUBROUTINE CDAMAT(CDAMA,TEMPG,PROPS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE DAMAGE COEFFICIENT (TEMPERATURE
C     DEPENDENT)
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
C**** INITIALIZATION
C
      CDAMA=0.0D0
C
      IPEP4=INT(PROPS(4))
C
      IPEP3=INT(PROPS(3))
      IF(IPEP3.EQ.1) THEN
       IF(IDAMG.EQ.0) RETURN
       IF(IDAMG.GT.20) RETURN
       GO TO (21,21,22,21,21,21) IDAMG
   21  CDAMA=VDAMA(1,1)
       RETURN
   22  CONTINUE
       RETURN
      ENDIF
C
      IF(IPEP4.EQ.0) THEN
       IF(IDAMG.EQ.0) RETURN
       IF(IDAMG.GT.20) RETURN
C
       GO TO (31,31,32,31,31,31) IDAMG
   31  CONTINUE
C
C**** COEFFICIENT CDAMA
C
       IF(TEMPG.LE.VDAMA(1,2)) CDAMA=VDAMA(1,1)
       DO IDAMA=2,NDAMA
        I1=IDAMA-1
        I2=IDAMA
        IF(TEMPG.GE.VDAMA(I1,2).AND.TEMPG.LE.VDAMA(I2,2)) 
     .   CDAMA=(VDAMA(I2,1)-VDAMA(I1,1))/(VDAMA(I2,2)-VDAMA(I1,2))*
     .          (TEMPG-VDAMA(I1,2))+VDAMA(I1,1)
       ENDDO
       IF(TEMPG.GE.VDAMA(NDAMA,2)) CDAMA=VDAMA(NDAMA,1)
C
C**** TEMPERATURE DERIVATIVE OF CDAMA
C
       DCDAM=0.0D0
       RETURN
C
   32  CONTINUE
       RETURN
      ENDIF             ! ipep4.eq.0
C
      RETURN
      END
