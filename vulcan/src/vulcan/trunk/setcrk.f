      SUBROUTINE SETCRK(ALENG,BETAM,ELCOD,DVOLU,NCRAK,NDIME,NNODE,
     .                  NSTR1,PROPS,SIGMA,STRAN,TCRAK,THICK)
C***********************************************************************
C
C****THIS ROUTINE DETERMINES HOW MANY CRACKS ARE OPENED AND THE
C    CORRESPONDING MATERIAL SYSTEM
C
C    Input variables:
C      
C      DVOLU              Volume of the integration point
C      ELCOD(NDIME,NNODE) Coordinates of the element nodes
C      PROPS(NPROP)       Material properties 
C      SIGMA(NSTR1)       'Trial' total stresses
C      STRAN(NSTR1)       Current total strains
C
C    Output variables:
C
C      NCRAK              Number of active cracks
C      TCRAK              Time for onset of cracking
C
C    Input/Output variables:
C
C      ALENG(NCRAK)       Characteristic lengths for cracks        
C      BETAM(NDIME,NDIME) Matrix containing the material system
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION ALENG(*), BETAM(NDIME,*), ELCOD(NDIME,*), PROPS(*),
     .          SIGMA(*), STRAN(*)
      DIMENSION SGPRI(3), TRANS(6,6)
      DATA MCRAK,ZEROM/3,1.0E-10/
C
C***FIND THE NUMBER OF PREVIOUSLY EXISTING CRACKS
C
      NCRAK=0
CC$DIR SCALAR
      DO 10 ICRAK=1,MCRAK
   10 IF(ALENG(ICRAK).GT.ZEROM) NCRAK=NCRAK+1
C
C***IF THE POINT WAS PREVIOUSLY UNCRACKED CHECK FOR CRACKING
C
      IF(NCRAK.EQ.0) THEN
C
C       (a.1) find principal stress values 
C
        CALL PRIVAL(NSTR1,SIGMA,SGPRI)
C
C       (a.2) determine number of new cracks
C
        CALL NUMCRK(MCRAK,NCRAK,NEWCK,NEWSY,NTYPE,PROPS,SIGMA,
     .              SGPRI,    0)
C
C       If cracking occurs:
C
        IF(NCRAK.GT.0) THEN
C
C         (a.3) Set up material system,
C
          CALL DIRCOS(NDIME,NSTR1,SIGMA,SGPRI,BETAM)
C
C         (a.4) Find characteristic lengths
C 
          CALL LENCRK(ALENG,BETAM,ELCOD,DVOLU,NCRAK,NDIME,NNODE,
     .                NTYPE,THICK,IELEM)
C
C         (a.5) Rotate strains to material system
C
          CALL TRASTR(STRAN,BETAM,NSTRE,NDIME,'G_TO_L')
C
C         (a.6) Recompute stresses neglecting Poisson's effect
C
          CALL NOPOIS(NCRAK,NSTRS,PROPS,SIGMA,STRAN)
C
C         (a.7) Record time for onset of cracking
C
          TCRAK=TTIME
C
        ENDIF
C
      ELSE
C
C***IF THE POINT WAS PREVIOUSLY CRACKED, CHECK IF NEW CRACKS ARE OPENING
C
C       (b.1) Rotate strains to material system
C
        CALL TRASTR(STRAN,BETAM,NSTRE,NDIME,'G_TO_L')
C
C       (b.2) Recompute stresses neglecting Poisson's effect
C
        CALL NOPOIS(NCRAK,NSTRS,PROPS,SIGMA,STRAN)
C
C       (b.3) determine new number of cracks
C
        CALL NUMCRK(MCRAK,NCRAK,NEWCK,NEWSY,NTYPE,PROPS,SIGMA,
     .              SGPRI,    1)
C
C       (b.4) If necessary rotate to new material system
C
C        IF(NEWSY.EQ.1) 
C     .    CALL ROTSYS(BETAM,SIGMA,TRANS,NDIME,NSTRE,NSTRS)
C
C       (b.5) If necessary find characteristic lengths
C 
        IF(NEWCK.EQ.1)
     .    CALL LENCRK(ALENG,BETAM,ELCOD,DVOLU,NCRAK,NDIME,NNODE,
     .                NTYPE,THICK,IELEM)
C
      ENDIF
C
      RETURN
      END
