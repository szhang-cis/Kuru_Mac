      SUBROUTINE LIHEATS(CARTDS,DEPSVS,SGTOTS,DSTRAS,ELDISS,
     .                   GPCODS,LARGES,NDIMES,NDOFNS,NNODES,NSTR1S,
     .                   PROPSS,SHAPES,STRANS,STRA0S,TSTRAS,XJACMS,
     .                   NDOFCS,SIGMAS,VELCMS,DMATXS,KDYNAS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES TEMPERATURE GRADIENTS AND INCREMENTAL HEATS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inte_oms.f'
C
      DIMENSION CARTDS(NDIMES,*), SGTOTS(*),        DSTRAS(*),
     .          ELDISS(NDOFCS,*), PROPSS(*),
     .          GPCODS(*),        SHAPES(*),        STRANS(*),
     .          STRA0S(*),        TSTRAS(*),        XJACMS(NDIMES,*),
     .          SIGMAS(*),        VELCMS(*),        DMATXS(NSTR1S,*)
C
C**** CALCULATE THE TOTAL TEMPERATURE GRADIENT IN TSTRAT
C
      DO IDIMES=1,NDIMES
       TSTRAS(IDIMES)=0.0
       DO IEVABS=1,NNODES
        TSTRAS(IDIMES)=TSTRAS(IDIMES)+CARTDS(IDIMES,IEVABS)*
     .                 ELDISS(NDOFCS,IEVABS)
       END DO
      END DO
C
C**** CALCULATE THE INCREMENTAL TEMPERATURE GRADIENT IN DSTRAT
C
      DO IDIMES=1,NDIMES
       DSTRAS(IDIMES)=0.0
       IF(KDYNAS.EQ.1) THEN
        DO IEVABS=1,NNODES
         DSTRAS(IDIMES)=DSTRAS(IDIMES)+CARTDS(IDIMES,IEVABS)*
     .                  VELCMS(IEVABS)*DTIMES
        END DO
       ENDIF
      END DO
C
C**** CALCULATE THE EFFECTIVE HEATS 
C
      DO IDIMES=1,NDIMES
       SGTOTS(IDIMES)=0.0
       SIGMAS(IDIMES)=0.0
       DO JDIMES=1,NDIMES
        SGTOTS(IDIMES)=SGTOTS(IDIMES)+DMATXS(IDIMES,JDIMES)*
     .                 TSTRAS(JDIMES)
        SIGMAS(IDIMES)=SIGMAS(IDIMES)+DMATXS(IDIMES,JDIMES)*
     .                 DSTRAS(JDIMES)
       ENDDO
      ENDDO
C
      IF(KDYNAS.EQ.1) THEN         ! transient
       IF(KINTES.EQ.1) THEN        ! Euler's method
        DO IDIMES=1,NDIMES
         SGTOTS(IDIMES)=0.0
         SIGMAS(IDIMES)=0.0
         DO JDIMES=1,NDIMES
          SGTOTS(IDIMES)=SGTOTS(IDIMES)+DMATXS(IDIMES,JDIMES)*
     .                  (TALFAS*TSTRAS(JDIMES)+
     .                 (1.0D+00-TALFAS)*(TSTRAS(JDIMES)-DSTRAS(JDIMES)))
          SIGMAS(IDIMES)=SIGMAS(IDIMES)+DMATXS(IDIMES,JDIMES)*
     .                   DSTRAS(JDIMES)
         ENDDO
        ENDDO
       ENDIF
      ENDIF
C
C**** CALCULATES THE THERMAL DISIPATION
C
      auxtts=0.0
      auxtis=0.0
      do idimes=1,ndimes
      auxtts=auxtts+sgtots(idimes)*tstras(idimes)
      auxtis=auxtis+sigmas(idimes)*dstras(idimes)
      end do
C
      dicais=auxtis
      dicats=auxtts
C   
      RETURN
      END
