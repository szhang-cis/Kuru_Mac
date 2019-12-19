      SUBROUTINE MASS05S(DVOLUS,PROPSS,SHAPES,WSTIFS,EHISTS,WHAPES,
     .                   ELDISS,VELCMS,TEINIS,FPCHLS,DVOLIS)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE HEAT CAPACITY MATRIX ( ELEMENT NO. 5 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nued_om.f'   ! thermal-microstructural
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION DVOLUS(*),        PROPSS(*),
     .          SHAPES(NNODLS,*)
      DIMENSION WSTIFS(*),        EHISTS(NHISTS,*)
      DIMENSION WHAPES(NNODLS,*)
      DIMENSION ELDISS(NDOFCS,*), VELCMS(*),
     .          TEINIS(NDOFCS,*), FPCHLS(NFPCH,*),
     .          DVOLIS(*)
      DIMENSION EMATX(27,27)
C
C**** LOOP ON INTEGRATION POINTS
C
      DO 100 IGAUSS=1,NGAULS
C
C**** COMPUTE TEMPERATURE AT GAUSSIAN POINT
C
      IF(NMEMO3S.EQ.0) THEN
C
       DTEMPS=0.0
       TGAUSS=0.0
       TGAUIS=0.0
       PSEUDOS=0.0
       DO INODLS=1,NNODLS
        DTEMPS=DTEMPS+SHAPES(INODLS,IGAUSS)*VELCMS(INODLS)
        TGAUSS=TGAUSS+SHAPES(INODLS,IGAUSS)*ELDISS(1,INODLS)
        TGAUIS=TGAUIS+SHAPES(INODLS,IGAUSS)*TEINIS(1,INODLS)
        IF(IFILLS.EQ.1) THEN
         IF(IMICR.EQ.0) THEN
          IPSEUS=2*NNUPTS+1
         ELSE
          IPSEUS=2*NNUPTS+NNUPO+1
         ENDIF
         PSEUDOS=PSEUDOS+SHAPES(INODLS,IGAUSS)*FPCHLS(IPSEUS,INODLS)
        ENDIF
       END DO
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
       IF(KDYNAS.EQ.1) THEN
        IF(KINTES.EQ.1)              ! Euler method
     .   TGAUSS=TALFAS*TGAUSS+(1.0-TALFAS)*(TGAUSS-DTEMPS*DTIMES)
       ENDIF
C
       ILAHES=0
       ISINRS=1
C
       CALL CAPCOFS( BASMM, BASCC,PROPSS,TGAUSS,
     .              DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .              DSOURS,COUTDS, BASMI,TGAUIS,PSEUDOS)
      ELSE
       IX=NDIMETOS-1
       BASMM=EHISTS(1,IGAUSS)
       BASMI=EHISTS(6+IX,IGAUSS)
       BASCC=EHISTS(2,IGAUSS)
      ENDIF                          ! nmemo3s.eq.0
C
C**** DETERMINES CAPACITY
C
      BASMM=BASMI
C
      IF(ITERME.GT.0) THEN                      ! bidirectional coupling

       call runends('iterme > 0 not implemented in mass05s.f')

c      IF(ITERMD.GT.0) THEN                     ! deformed shape
c       IF(LARGET.EQ.1) THEN                    ! TLF
c        BASMM=BASMI
c       ENDIF
c       IF(LARGET.EQ.2) THEN                    ! ULF
c        detjj=1.0      ! determinant of F (from mechanical computation)
c        BASMM=BASMI/DETJJ
c       ENDIF
c       IF(LARGET.EQ.3) THEN                    ! Eulerian
c        BASMM=BASMI*DVOLIT(IGAUST)/DVOLUT(IGAUST)
c       ENDIF
c      ENDIF
      ENDIF
C
C**** ELEMENTAL HEAT CAPACITY MATRIX
C
      CALL WMHEATT(DVOLUS(IGAUSS),NDOFNS,NEVABS,NNODLS,PROPSS,
     .             SHAPES(1,IGAUSS),WSTIFS,BASCC,
     .             BASMM,COUTTS,
     .             WHAPES(1,IGAUSS),KSYMMS,EMATX)
C
  100 CONTINUE
C
      RETURN
      END
