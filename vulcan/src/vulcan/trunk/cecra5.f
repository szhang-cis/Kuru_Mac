      SUBROUTINE CECRA5(DVOLU,EHIST,ELCOD,KGAUS,NDIME,NNODE,NSTR1,
     .                  PROPS,SGTOT,SIGMA,STRAN,THICK)
C***********************************************************************
C
C****THIS ROUTINE EVALUATES TOTAL STRESSES FOR ELASTO/BRITTLE MATERIALS
C                     NCRIT= 5
C
C    Input variables:
C      
C      DVOLU               Volume of the integration point
C      ELCOD(NDIME,NNODE)  Coordinates of the element nodes
C      PROPS(NPROP)        Material properties
C      SIGMA(NSTR1)        'Trial' total stresses
C      STRAN(NSTR1)        Current total strains
C
C    Output variables:
C
C      SGTOT(NSTR1)        Current total stresses
C
C    Input/Output variables:
C
C      EHIST(1:9)          Transformation matrix
C      EHIST(10:12)        Current strains in cracks
C      EHIST(13:15)        Secant stiffness moduli for cracks
C      EHIST(16:18)        Characteristic lengths for cracks        
C      EHIST(19)           Time of onset of cracking
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'auxl_om.f'
      DIMENSION EHIST(*), ELCOD(NDIME,*), PROPS(*), 
     .          SGTOT(*), SIGMA(*),       STRAN(*)
C
C   CHECK IF THE POINT IS CRACKED AND HOW MANY CRACKS ARE PRESENT
C
C
      CALL SETCRK(EHIST(16),EHIST(1),ELCOD,DVOLU,NCRAK,NDIME,NNODE,
     .            NSTR1,PROPS,SIGMA,STRAN,EHIST(19),THICK)   
C
C   DEAL WITH EXISTING CRACKS
C
      CALL FIXCRK(EHIST(1),EHIST(10),EHIST(13),EHIST(16),
     .            NCRAK,NDIME,NSTR1,PROPS,SGTOT,SIGMA,STRAN)
C
      RETURN
      END
