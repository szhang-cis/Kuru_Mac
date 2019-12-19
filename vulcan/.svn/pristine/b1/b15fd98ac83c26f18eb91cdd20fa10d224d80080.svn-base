      SUBROUTINE CONSET
C***********************************************************************
C
C**** THIS ROUTINE SETS-UP MOST OF THE REMAINING CONTROLLING PARAMETERS
C
C     Notes:
C
C     HOW TO COMPUTE NHIST (FOR KPROB=1):
C
C     NHIST=NSTR1+NSTR1*(NSTR1+1)/2+NBASE+NMECH   ! symm. const. tensor
C     NHIST=NSTR1+NSTR1*NSTR1+NBASE+NMECH         ! unsy. const. tensor
C
C     (NBASE defined in pointe.f; NMECH=1 for ITERME > 0)
C
C     MORE CORRECTLY:
C
C     NHIST=NSTR1+NSTRS*(NSTRS+1)/2+NBASE+NMECH   ! symm. const. tensor
C     NHIST=NSTR1+NSTRS*NSTRS+NBASE+NMECH         ! unsy. const. tensor
C
C     EXAMPLE FOR THE 2D CASE (ITERME>0):
C
C     NHIST=4+10+6+1=21   (symmetric const. tensor)
C     NHIST=4+16+6+1=31   (unsymmetric const. tensor)
C
C     EXAMPLE FOR THE 3D CASE (ITERME>0):
C
C     NHIST=6+21+8+1=36   (symmetric const. tensor)
C     NHIST=6+36+8+1=51   (unsymmetric const. tensor)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)                  ! not strictly necessary
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
C**** DEPENDING ON KPROB DEFINE SOME PARAMETERS
C
      IF(KPROB.EQ.0) THEN                   ! SEEPAGE
       NDOFN=1 
       NDOFC=NDOFN
       NSTR1=NDIME
       NHIST=NDIME
      ELSE IF(KPROB.EQ.1) THEN              ! STRUCTURAL,CONTINUUM
       NDOFN=NDIME
       NDOFC=NDOFN
       IF(KPORE.EQ.2)NDOFC=NDOFN+1
       IF(NDIME.EQ.1) NSTR1=1
       IF(NDIME.EQ.2) NSTR1=4 
       IF(NDIME.EQ.3) NSTR1=6
c      NHIST=36                     ! 3D symm. const. tensor, ITERME>0
c      NHIST=51                     ! 3D unsymm. const. tensor, ITERME>0
       NHIST=52                     ! 3D unsymm. const. tensor, NCHAI=3
c      NHIST=98                     ! 3D symm. const. tensor, dual phase
      ELSE IF(KPROB.EQ.2) THEN              ! STRUCTURAL,SHELLS
       IF(NDIME.EQ.2) THEN                  ! plates  
        NDOFN=3
        NSTR1=5
       ELSE                                 ! shells
        NDOFN=6
        NSTR1=8
       ENDIF   
       NDOFC=NDOFN
       NHIST=0
      ELSE IF(KPROB.EQ.3) THEN              ! SONIC WAVES
       NDOFN=1 
       NDOFC=NDOFN
       NSTR1=0
       NHIST=0
      ELSE IF(KPROB.EQ.4) THEN              ! THERMAL
       NDOFN=1 
       NDOFC=NDOFN
       NSTR1=NDIME
       NHIST=20
      ELSE IF(KPROB.EQ.5) THEN              ! INCOMPRESSIBILITY
       NDOFN=NDIME
       NDOFC=NDOFN
       NSTR1=1                        ! always
       NHIST=0
      ENDIF
C
C**** DEFINE SOME OTHER PARAMETERS 
C
      NEVAB=NNODE*NDOFN
      NEVAC=NNODE*NDOFC
      NTOTV=NPOIN*NDOFC
      NTOTG=NELEM*NGAUS
      NKOND=(NNODE+1)*NNODE/2    ! or NKOND=(NNODE-1)*NNODE/2+NNODE
      NMOVA=(NEVAB+1)*NEVAB/2    ! or NMOVA=(NEVAB-1)*NEVAB/2+NEVAB
      NKOVA=NMOVA
      IF(KSYMM.EQ.0) NKOVA=NEVAB*NEVAB   ! unsymmetric case
C
C**** WRITE SOME MAIN PROBLEM PARAMETERS 
C
      WRITE(LURES,905)     NPOIN,NELEM,NNODE,NGAUS,NGRUP,NMATS,NFUNC
      WRITE(LURES,910)     NDIME,NDOFC,NEVAC,NTOTV,NTOTG
      
#ifdef restricted
      WRITE(LURES,*) 'MAXELEM', MAXELEM, 'NELEM', NELEM
      if (NELEM.GT.MAXELEM)    ! defined at the preprocesor level
     . CALL RUNEND('RESTRICTED VERSION MAX ELEMENTS EXCEEDED')
      WRITE(LURES,*) 'MAXNODE', MAXNODE, 'NPOIN', NPOIN
      if (NPOIN.GT.MAXNODE)    ! defined at the preprocesor level
     . CALL RUNEND('RESTRICTED VERSION MAX NODES EXCEEDED')
#endif
C
      RETURN
  905 FORMAT(//5X,30HMAIN CONTROLLING PARAMETERS  :,//
     .        15X,30HNUMBER OF NODAL POINTS       =,I7/
     .        15X,30HNUMBER OF ELEMENTS           =,I7/
     .        15X,30HNODES PER ELEMENT (MAXIMUM)  =,2X,I5/
     .        15X,30HINT. PTS. PER ELEMENT (MAX.) =,2X,I5/
     .        15X,30HNUMBER OF ELEMENT  SETS      =,2X,I5/
     .        15X,30HNUMBER OF MATERIAL SETS      =,2X,I5/
     .        15X,30HNUMBER OF LOADING FUNCTIONS  =,2X,I5/)
  910 FORMAT(//5X,31HOTHER CONTROLLING PARAMETERS  :,//
     .        15X,30HNUMBER OF SPATIAL DIMENSIONS =,2X,I5/
     .        15X,30HNUMBER OF D.O.F PER NODE     =,2X,I5/
     .        15X,30HNUMBER OF D.O.F PER ELEMENT  =,2X,I5/
     .        15X,30HTOTAL NUMBER OF D.O.F.       =,I7/
     .        15X,30HTOTAL NUMBER OF INTEGT. PTS. =,I7/)
      END
