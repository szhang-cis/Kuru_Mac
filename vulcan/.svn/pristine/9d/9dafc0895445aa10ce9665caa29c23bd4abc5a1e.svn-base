      SUBROUTINE CONSETS
C***********************************************************************
C
C**** THIS ROUTINE SETS-UP MOST OF THE REMAINING CONTROLLING PARAMETERS
C
C***********************************************************************
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'     ! thermal-microstructural
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
C**** DEPENDING ON KPROB DEFINE SOME PARAMETERS
C
      IF(KPROBS.EQ.0) THEN                  ! SEEPAGE
       NDOFNS=1
       NDOFCS=NDOFNS
       NSTR1S=0
       NHISTS=0
      ELSE IF(KPROBS.EQ.1) THEN             ! STRUCTURAL,CONTINUUM
       NDOFNS=NDIMES
       NDOFCS=NDOFNS                        ! ndofc
       IF(KPORES.EQ.2)NDOFCS=NDOFNS+1
       IF(NDIMES.EQ.1) NSTR1S=1             ! nstr1
       IF(NDIMES.EQ.2) NSTR1S=4 
       IF(NDIMES.EQ.3) NSTR1S=6
       NHISTS=32                            ! nhist
      ELSE IF(KPROBS.EQ.2) THEN             ! STRUCTURAL,SHELLS
       NDOFNS=6
       NDOFCS=NDOFNS
       NSTR1S=8
       NHISTS=0
      ELSE IF(KPROBS.EQ.3) THEN             ! SONIC WAVES
       NDOFNS=1
       NDOFCS=NDOFNS
       NSTR1S=0
       NHISTS=0
      ELSE IF(KPROBS.EQ.4) THEN             ! THERMAL
       NDOFNS=1
       NDOFCS=NDOFNS
       NSTR1S=NDIMES
c      NHISTS=10+IMICR*50      ! 1-10: thermal var.; 10-60: microst.var.
       NHISTS=0                ! to be improved
      ENDIF
C
C**** DEFINE SOME OTHER PARAMETERS 
C
      NEVABS=NNODES*NDOFNS
      NEVACS=NNODES*NDOFCS
      NTOTVS=NPOINS*NDOFCS
      NTOTGS=NELEMS*NGAUSS
      NKONDS=(NNODES+1)*NNODES/2  ! or NKONDS=(NNODES-1)*NNODES/2+NNODES
      NMOVAS=(NEVABS+1)*NEVABS/2  ! or NMOVAS=(NEVABS-1)*NEVABS/2+NEVABS
      NKOVAS=NMOVAS
      IF(KSYMMS.EQ.0) NKOVAS=NEVABS*NEVABS   ! unsymmetric case
C
      NTOTVMS=NTOTVS*NDIMES    ! redefined in couinc.f for coupled prob.
      NDOFCMS=NDOFCS*NDIMES
C
C**** WRITE SOME MAIN PROBLEM PARAMETERS 
C
      WRITE(LURESS,905)      NPOINS,NELEMS,NNODES,NGAUSS,NGRUPS,NMATSS,
     .                       NFUNCS
      WRITE(LURESS,910)      NDIMES,NDOFCS,NEVACS,NTOTVS,NTOTGS
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
