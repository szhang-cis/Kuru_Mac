      SUBROUTINE ASEL04(CSTIF,ESTIF,WSTIF,PROPS,CSTI1,LNODS,
     .                  AUXMA,TRAMA,INFRI,COFRI)
C***********************************************************************
C
C**** THIS ROUTINE ASSEMBLES ELEMENTAL "EFFECTIVE" MATRIX  FOR
C     UNCOUPLED PROBLEMS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION CSTIF(NEVAB,NEVAB), ESTIF(NKOVA),
     .          WSTIF(NKOVA),       PROPS(NPROP,*)
      DIMENSION CSTI1(NEVAB,NEVAB)
      DIMENSION LNODS(NNODE)
C
      DIMENSION AUXMA(NEVAB,NEVAB), TRAMA(NEVAB,NEVAB),
     .          A1(3,3)
      DIMENSION INFRI(NPOIN),       COFRI(NSKEW,NDIME,*)
C
      IF(NMEMO6M.EQ.0) RETURN
C
C**** ASSEMBLE ELEMENTAL MATRIX
C
      IF(KDYNA.EQ.1) THEN          ! dynamic problems
       CTIM1=TALFA/(DTIME*DTIME)
       CTIM2=TBETA/ DTIME
       CTIM3=TGAMA
C
       CONS1=0.0D0
       CONS2=CTIM3
      END IF
C
C**** CHECK THE CORRECTNESS OF THE CONSTITUENT ELEMENTAL MATRICES
C
      ICHEK=0
      IF(ICHEK.EQ.1) THEN
       IF(KPROB.GE.1)
     .  CALL TESTM1(ESTIF,IELEM,NEVAB)
       IF(KDYNA.NE.0)
     .  CALL TESTM1(WSTIF,IELEM,NEVAB)
      ENDIF
C
      IF(NMEMO7M.EQ.0) THEN
C
C**** INITIALISE ASSEMBLED ELEMENTAL MATRIX
C
      DO IEVAB=1,NEVAB
       DO JEVAB=1,NEVAB
        CSTIF(IEVAB,JEVAB)=0.0D0
       ENDDO
      ENDDO
C
      IF(KSYMM.EQ.0) THEN              ! unsymmetric case
C
       IF(KDYNA.EQ.0) THEN
        IKOVA=0
        DO IEVAB=1,NEVAB
         DO JEVAB=1,NEVAB
          IKOVA=IKOVA+1
          CSTIF(IEVAB,JEVAB)=ESTIF(IKOVA)
         ENDDO
        ENDDO
       ELSE ! Dynamic problem
        IKOVA=0
        DO IEVAB=1,NEVAB
         DO JEVAB=1,NEVAB
          IKOVA=IKOVA+1
          CSTIF(IEVAB,JEVAB)=CONS1*WSTIF(IKOVA)+CONS2*ESTIF(IKOVA)
         ENDDO
        ENDDO
       ENDIF
C
      ELSE                             ! symmetric case
C
       IF(KDYNA.EQ.0) THEN
        IKOVA=0
        DO IEVAB=1,NEVAB
         DO JEVAB=IEVAB,NEVAB
          IKOVA=IKOVA+1
          CSTIF(IEVAB,JEVAB)=ESTIF(IKOVA)
          CSTIF(JEVAB,IEVAB)=CSTIF(IEVAB,JEVAB)
         ENDDO
        ENDDO
       ELSE ! Dynamic problem
        IKOVA=0
        DO IEVAB=1,NEVAB
         DO JEVAB=IEVAB,NEVAB
          IKOVA=IKOVA+1
          CSTIF(IEVAB,JEVAB)=CONS1*WSTIF(IKOVA)+CONS2*ESTIF(IKOVA)
          CSTIF(JEVAB,IEVAB)=CSTIF(IEVAB,JEVAB)
         ENDDO
        ENDDO
       ENDIF
C
      ENDIF
C
C**** PERFORMS LOCAL TRANSFORMATION (IF NECESSARY)
C
      IF(NSKEW.GT.0) THEN
C
       DO IEVAB=1,NEVAB
        DO JEVAB=1,NEVAB
         TRAMA(IEVAB,JEVAB)=0.0D0
        ENDDO
       ENDDO
C
       DO INODE=1,NNODL
        IPOIN=LNODS(INODE)
        IOEF1=INFRI(IPOIN)
C
        IF(IOEF1.EQ.0) THEN
         DO IDOFN=1,NDOFN
          DO JDOFN=1,NDOFN
           A1(IDOFN,JDOFN)=0.0D0
           IF(IDOFN.EQ.JDOFN) A1(IDOFN,JDOFN)=1.0D0
          ENDDO
         ENDDO
        ELSE
         DO IDOFN=1,NDOFN
          DO JDOFN=1,NDOFN
           A1(IDOFN,JDOFN)=COFRI(IOEF1,IDOFN,JDOFN)
          ENDDO
         ENDDO
        ENDIF
C
        DO IDOFN=1,NDOFN
         DO JDOFN=1,NDOFN
          IAUXM=NDOFN*(INODE-1)+IDOFN
          JAUXM=NDOFN*(INODE-1)+JDOFN
          TRAMA(IAUXM,JAUXM)=A1(IDOFN,JDOFN)
         ENDDO
        ENDDO
       ENDDO
C
C**** TRANSFORM GLOBAL JACOBIAN MATRIX TO LOCAL COORDINATE SYSTEM
C
       DO IEVAB=1,NEVAB
        DO JEVAB=1,NEVAB
         AUXMA(IEVAB,JEVAB)=0.0D0
         DO KEVAB=1,NEVAB
          DO LEVAB=1,NEVAB
            AUXMA(IEVAB,JEVAB)=AUXMA(IEVAB,JEVAB)+TRAMA(IEVAB,KEVAB)*
     .                         CSTIF(KEVAB,LEVAB)*TRAMA(JEVAB,LEVAB)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
C
       DO IEVAB=1,NEVAB
        DO JEVAB=1,NEVAB
         CSTIF(IEVAB,JEVAB)=AUXMA(IEVAB,JEVAB)
        ENDDO
       ENDDO
C
      ENDIF
C
      ELSE                  ! nmemo7m=1
C
C**** INITIALISE ASSEMBLED ELEMENTAL MATRIX
C
      DO IEVAB=1,NEVAB
       DO JEVAB=1,NEVAB
        CSTI1(IEVAB,JEVAB)=0.0D0
       ENDDO
      ENDDO
C
      IF(KSYMM.EQ.0) THEN          ! unsymmetric case
C
       IF(KDYNA.EQ.0) THEN         ! steady-state case
        IKOVA=0
        DO IEVAB=1,NEVAB
         DO JEVAB=1,NEVAB
          IKOVA=IKOVA+1
          CSTI1(IEVAB,JEVAB)=ESTIF(IKOVA)
         ENDDO
        ENDDO
       ELSE                        ! transient case
        IKOVA=0
        DO IEVAB=1,NEVAB
         DO JEVAB=1,NEVAB
          IKOVA=IKOVA+1 
          CSTI1(IEVAB,JEVAB)=CTIM2*WSTIF(IKOVA)+
     .                       CTIM3*ESTIF(IKOVA)
         ENDDO
        ENDDO
       ENDIF
C
      ELSE                         ! symmetric case
C
       IF(KDYNA.EQ.0) THEN         ! steady-state case
        IKOVA=0
        DO IEVAB=1,NEVAB
         DO JEVAB=IEVAB,NEVAB
          IKOVA=IKOVA+1
          CSTI1(IEVAB,JEVAB)=ESTIF(IKOVA)
          CSTI1(JEVAB,IEVAB)=CSTI1(IEVAB,JEVAB)
         ENDDO
        ENDDO
       ELSE                        ! transient case
        IKOVA=0
        DO IEVAB=1,NEVAB
         DO JEVAB=IEVAB,NEVAB
          IKOVA=IKOVA+1 
          CSTI1(IEVAB,JEVAB)=CTIM2*WSTIF(IKOVA)+
     .                       CTIM3*ESTIF(IKOVA)
          CSTI1(JEVAB,IEVAB)=CSTI1(IEVAB,JEVAB)
         ENDDO
        ENDDO
       ENDIF
C
      ENDIF
C
C**** PERFORMS LOCAL TRANSFORMATION (IF NECESSARY)
C
      IF(NSKEW.GT.0) THEN
C
       DO IEVAB=1,NEVAB
        DO JEVAB=1,NEVAB
         TRAMA(IEVAB,JEVAB)=0.0D0
        ENDDO
       ENDDO
C
       DO INODE=1,NNODL
        IPOIN=LNODS(INODE)
        IOEF1=INFRI(IPOIN)
C
        IF(IOEF1.EQ.0) THEN
         DO IDOFN=1,NDOFN
          DO JDOFN=1,NDOFN
           A1(IDOFN,JDOFN)=0.0D0
           IF(IDOFN.EQ.JDOFN) A1(IDOFN,JDOFN)=1.0D0
          ENDDO
         ENDDO
        ELSE
         DO IDOFN=1,NDOFN
          DO JDOFN=1,NDOFN
           A1(IDOFN,JDOFN)=COFRI(IOEF1,IDOFN,JDOFN)
          ENDDO
         ENDDO
        ENDIF
C
        DO IDOFN=1,NDOFN
         DO JDOFN=1,NDOFN
          IAUXM=NDOFN*(INODE-1)+IDOFN
          JAUXM=NDOFN*(INODE-1)+JDOFN
          TRAMA(IAUXM,JAUXM)=A1(IDOFN,JDOFN)
         ENDDO
        ENDDO
       ENDDO
C
C**** TRANSFORM GLOBAL JACOBIAN MATRIX TO LOCAL COORDINATE SYSTEM
C
       DO IEVAB=1,NEVAB
        DO JEVAB=1,NEVAB
         AUXMA(IEVAB,JEVAB)=0.0D0
         DO KEVAB=1,NEVAB
          DO LEVAB=1,NEVAB
            AUXMA(IEVAB,JEVAB)=AUXMA(IEVAB,JEVAB)+TRAMA(IEVAB,KEVAB)*
     .                         CSTI1(KEVAB,LEVAB)*TRAMA(JEVAB,LEVAB)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
C
       DO IEVAB=1,NEVAB
        DO JEVAB=1,NEVAB
         CSTI1(IEVAB,JEVAB)=AUXMA(IEVAB,JEVAB)
        ENDDO
       ENDDO
C
      ENDIF
C
      ENDIF             ! nmemo7m.eq.0
C
      RETURN
      END
