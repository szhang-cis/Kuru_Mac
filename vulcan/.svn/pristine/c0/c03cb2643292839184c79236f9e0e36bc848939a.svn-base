      SUBROUTINE CONSETT
C***********************************************************************
C
C**** THIS ROUTINE SETS-UP MOST OF THE REMAINING CONTROLLING PARAMETERS
C
C***********************************************************************
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'     ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
C**** DEPENDING ON KPROB DEFINE SOME PARAMETERS
C
      IF(KPROBT.EQ.0) THEN                  ! SEEPAGE
       NDOFNT=1
       NDOFCT=NDOFNT
       NSTR1T=0
       NHISTT=0
      ELSE IF(KPROBT.EQ.1) THEN             ! STRUCTURAL,CONTINUUM
       NDOFNT=NDIMET
       NDOFCT=NDOFNT                        ! ndofc
       IF(KPORET.EQ.2)NDOFCT=NDOFNT+1
       IF(NDIMET.EQ.1) NSTR1T=1             ! nstr1
       IF(NDIMET.EQ.2) NSTR1T=4 
       IF(NDIMET.EQ.3) NSTR1T=6
       NHISTT=32                            ! nhist
      ELSE IF(KPROBT.EQ.2) THEN             ! STRUCTURAL,SHELLS
       NDOFNT=6
       NDOFCT=NDOFNT
       NSTR1T=8
       NHISTT=0
      ELSE IF(KPROBT.EQ.3) THEN             ! SONIC WAVES
       NDOFNT=1
       NDOFCT=NDOFNT
       NSTR1T=0
       NHISTT=0
      ELSE IF(KPROBT.EQ.4) THEN             ! THERMAL
       NDOFNT=1
       NDOFCT=NDOFNT
       NSTR1T=NDIMET
       NHIST1=10                     ! thermal var.; change if necessary
c      NHISTT=NHIST1*NMEMO3+NNUM4TMS ! + microst. var.; see setdatd.f
       NHISTT=NHIST1*NMEMO3
       IF(IMICR.EQ.1)
     .  NHISTT=NHISTT+NNUM4TMS       ! + microst. var.; see setdatd.f
      ENDIF
C
C**** DEFINE SOME OTHER PARAMETERS
C
      NEVABT=NNODET*NDOFNT
      NEVACT=NNODET*NDOFCT
      NTOTVT=NPOINT*NDOFCT
      NTOTGT=NELEMT*NGAUST
      NKONDT=(NNODET+1)*NNODET/2  ! or NKONDT=(NNODET-1)*NNODET/2+NNODET
      NMOVAT=(NEVABT+1)*NEVABT/2  ! or NMOVAT=(NEVABT-1)*NEVABT/2+NEVABT
      NKOVAT=NMOVAT
      IF(KSYMMT.EQ.0) NKOVAT=NEVABT*NEVABT   ! unsymmetric case
C
      NTOTVM=NTOTVT*NDIMET     ! redefined in couinc.f for coupled prob.
      NDOFCM=NDOFCT*NDIMET
C
C**** WRITE SOME MAIN PROBLEM PARAMETERS 
C
      WRITE(LUREST,905)      NPOINT,NELEMT,NNODET,NGAUST,NGRUPT,NMATST,
     .                       NFUNCT
      WRITE(LUREST,910)      NDIMET,NDOFCT,NEVACT,NTOTVT,NTOTGT
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
