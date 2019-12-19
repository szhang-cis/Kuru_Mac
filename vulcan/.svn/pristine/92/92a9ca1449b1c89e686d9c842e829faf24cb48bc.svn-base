      SUBROUTINE FRDY30(PROPS,ELELM,ACELM,DVOLU,SHAPE)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE INTERNAL DYNAMIC FORCES
C    
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROPS(*), ELELM(*), ACELM(NDOFC,*),
     .          DVOLU(*), SHAPE(NNODL,*)
      DIMENSION ACELG(3)
C
      DENSE=PROPS(12)
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO 100 IGAUS=1,NGAUL
C
      DO IDOFN=1,NDOFN
       ACELG(IDOFN)=0.0
       DO INODL=1,NNODL
        ACELG(IDOFN)=ACELG(IDOFN)+SHAPE(INODL,IGAUS)*ACELM(IDOFN,INODL)
       END DO
      END DO
C
C**** INTEGRATE THE ACCELERATION INTO THE INTERNAL DYNAMIC FORCES
C
      DO INODL=1,NNODL
       IEVAB=(INODL-1)*NDOFN
       DO IDOFN=1,NDOFN
        IEVAB=IEVAB+1
        ELELM(IEVAB)=ELELM(IEVAB)+SHAPE(INODL,IGAUS)*DENSE*ACELG(IDOFN)*
     .               DVOLU(IGAUS)
       END DO 
      END DO 
C
  100 CONTINUE
C
      RETURN
      END
