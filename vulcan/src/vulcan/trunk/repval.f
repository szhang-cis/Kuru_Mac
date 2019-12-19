      SUBROUTINE REPVAL(IROTA,SGPRI,SIGMA)
C***********************************************************************
C
C**** THIS ROUTINE STORES THE REPRESENTATIVE VALUES OF NORMAL STRESSES
C     IN THE MATERIAL SYSTEM IN SGPRI
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SGPRI(*), SIGMA(*)
C
C**** THE MATERIAL SYSTEM IS ALREADY FIXED
C
      IF(IROTA.EQ.0) THEN
        SGPRI(1)=SIGMA(1)
        SGPRI(2)=SIGMA(2)
        SGPRI(3)=SIGMA(4)
C
      ELSE
C
C**** THE MATERIAL SYSTEM IS NOT FIXED YET     
C
C
C       LOOK FOR MAXIMUM & MINIMUM STRESSES PARALLEL TO CRACK    
C
        S11=SIGMA(2)
        S22=SIGMA(4)
        S12=SIGMA(6)
        SGPRI(1)=SIGMA(1)
        SGPRI(2)=(S11+S22)/2.0+SQRT(0.25*(S11-S22)**2+S12**2)
        SGPRI(3)=(S11+S22)/2.0-SQRT(0.25*(S11-S22)**2+S12**2)
C
      ENDIF
C
      RETURN
      END
