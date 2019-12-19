      SUBROUTINE ASELM1(MATNO,PROEL,PROPS,CSTIF,ESTIF,WSTIF,LNODS,
     .                  AUXMA,TRAMA,INFRI,COFRI)
C***********************************************************************
C
C**** THIS ROUTINE ASSEMBLES ELEMENTAL "EFFECTIVE" MATRIX  FOR
C     UNCOUPLED PROBLEMS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION MATNO(*),           PROEL(NPREL,*),
     .          PROPS(NPROP,*)
      DIMENSION CSTIF(NEVAB,NEVAB), ESTIF(NKOVA),
     .          WSTIF(NKOVA),       LNODS(NNODE,NELEM)
C
      DIMENSION AUXMA(NEVAB,NEVAB), TRAMA(NEVAB,NEVAB),
     .          A1(3,3)
      DIMENSION INFRI(NPOIN),       COFRI(NSKEW,NDIME,*)
C
C**** ASSEMBLE ELEMENTAL MATRIX
C
      IF(KDYNA.EQ.1) THEN
       CTIM1=TALFA/(DTIME*DTIME)
       CTIM2=TBETA/DTIME
       CTIM3=TGAMA
      ENDIF
C
C**** LOOP OVER ELEMENTS
C
      DO 1000 IELEM=1,NELEM
      LGRUP=MATNO(IELEM)
      LMATS=INT(PROEL(1,LGRUP))
      ITYPE=INT(PROEL(5,LGRUP))
      IF(KDYNA.EQ.1) THEN          ! dynamic problems
       IF(ITYPE.LT.100) THEN
        CONS1=CTIM1+CTIM2*(PROPS(21,LMATS))
        CONS2=CTIM3+CTIM2*(PROPS(22,LMATS))
       ELSE ! Boundary element
        CONS1=0.0D0
        CONS2=CTIM2
       END IF
       IF(ITYPE.EQ.4) THEN         ! contact element in dynamic problems
        CONS1=0.0D0
        CONS2=CTIM3
       ENDIF
      END IF
C
C**** READ MATRICES FROM DATA BASE
C 
                     CALL DATBAS(ESTIF,    7,    2)
      IF(KDYNA.EQ.1) CALL DATBAS(WSTIF,    8,    2)
C
C**** CHECK THE CORRECTNESS OF THE CONSTITUENT ELEMENTAL MATRICES
C
      ICHEK=0
      IF(ICHEK.EQ.1) THEN
        IF(KPROB.GE.1)
     .    CALL TESTM1(ESTIF,IELEM,NEVAB)
        IF(KDYNA.NE.0)
     .    CALL TESTM1(WSTIF,IELEM,NEVAB)
      ENDIF
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
       DO INODE=1,NNODE
        IPOIN=LNODS(INODE,IELEM)
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
C**** WRITE CSTIF TO DATA BASE
C 
      CALL DATBAS(CSTIF,    6,    1)
C
 1000 CONTINUE
C
      RETURN
      END
