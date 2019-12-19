      SUBROUTINE ASELM3(MATNO,PROEL,PROPS,CSTIF,ESTIF,WSTIF)
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
     .          WSTIF(NKOVA)
C
C**** ASSEMBLE ELEMENTAL MATRIX
C
      IF(KDYNA.EQ.1) THEN
ctm        CTIM1= TALFA/(DTIME**2)
        CTIM2= 1.0D00 / DTIME
        CTIM3= 1.0
      ENDIF
C
C**** LOOP OVER ELEMENTS
C
      DO 1000 IELEM=1,NELEM
      LGRUP=MATNO(IELEM)
      LMATS=INT(PROEL(1,LGRUP))
ctm      CODA1=CTIM2*(PROPS(21,LMATS))
ctm      CODA2=CTIM2*(PROPS(22,LMATS))
ctm poner coda1=0 (no se si son comentarios mios !!)
ctm       coda2=0
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
          CSTIF(IEVAB,JEVAB)=0.0
        ENDDO
      ENDDO
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
      ELSE
        IKOVA=0
        DO IEVAB=1,NEVAB
          DO JEVAB=IEVAB,NEVAB
            IKOVA=IKOVA+1 
            CSTIF(IEVAB,JEVAB)=CTIM2*WSTIF(IKOVA)+
     .                         CTIM3*ESTIF(IKOVA)
ctm              CSTIF(IEVAB,JEVAB)=(CTIM1+CODA1)*WSTIF(IKOVA)+
ctm     .                         (CTIM3+CODA2)*ESTIF(IKOVA)
            CSTIF(JEVAB,IEVAB)=CSTIF(IEVAB,JEVAB)
          ENDDO
        ENDDO
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
