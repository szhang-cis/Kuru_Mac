      SUBROUTINE MEMODT
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE MEMORY REQUIREMENTS & ADDRESS FOR 
C     DATABASE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_om.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
C**** TRY TO ALLOCATE DATA BASE IN VIRTUAL MEMORY
C
      ISTOP=NDISKDM    ! ndiskdm=0 out-of-core; ndiskdm=1 virtual memory
C
      NBLOC=11
      DO 10 IBLOC=NBLIM,1,-1
C
C**** CHECK POSSIBILITY TO ALLOCATE IBLOC ARRAYS IN VIRTUAL MEMORY
C
        IF((IBLOC.EQ.6).AND.(KDYNA.EQ.0).AND.(KPORE.NE.2)) GO TO 10
C
        NWORP=IDATP(IBLOC+1,3)
C
C**** MEMORY CONTROL
C
        IF(NWORP.GT.MWORP)
     .   CALL RUNEND('ERROR IN MEMODT: NWORP           ')
C
        LDABA=NWORP*NELEM+IDATP(12,3)
        LBYTB=LDABA*8
C
C**** SET APPROPRIATE FLAGS AND POINTERS FOR DATABASE
C                                 ISTOP=1 >> MEMORY HAS BEEN ALLOCATED
        IF(ISTOP.EQ.1) THEN
          DO ILOOP=1,IBLOC
            IDATP(ILOOP,5)=1
          ENDDO
          NSKIP=IDATP(IBLOC+1,4)
          NLENP=NLENP-NSKIP
          NRECP=NLENP*NELEM
          DO ILOOP=NBLOC+1,IBLOC+1,-1
            IDATP(ILOOP,4)=IDATP(ILOOP,4)-NSKIP
          ENDDO
          GO TO 100
        ENDIF
   10 CONTINUE
C
C**** ALLOCATE VM FOR DISTO
C
      IBLOC=0
      LDABA=IDATP(12,3)
      LBYTB=LDABA*8
C
  100 CONTINUE
C
      IDATP(12,5)=1
C
      RETURN
      END
