      SUBROUTINE STIF05S(CARTDS,DVOLUS,EHISTS,ELDISS,EPMTXS,GPCODS,
     .                   PROPSS,SHAPES,STRSGS,ESTIFS,HSTIFS,ELCODS,
     .                   BMATXS,DMATXS,SIGMAS,XJACMS,VELCMS,
     .                   WARTDS,TEINIS,FPCHLS,DVOLIS,AUXS1S)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE CONDUCTIVITY MATRIX
C     ( ELEMENT TYPE NO. 5 )
C
C***********************************************************************
C
C     Index of variables:
C
C     IX=NDIMETO-1
C
C     EHISTT(   1) = Density
C     EHISTT(   2) = Specific Heat coefficient
C     EHISTT(   3) = Isotropic Conductivity or Conduct. x (orthot. mat.)
C     EHISTT(   4) = Conductivity y (orthotropic material)
C     EHISTT(   5) = Conductivity z (orthotropic material)
C     EHISTT(3:11) = Conductivity for the fully anisotropc mat. (3D)
C     EHISTT(4+IX) = L*Phase-change function
C     EHISTT(5+IX) = L*Phase-change function rate
C     EHISTT(6+IX) = Initial density
C     EHISTT(7+IX) = Coupling coefficient
C     EHISTT(8+IX) = Temperature derivative of phase-change function
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
      DIMENSION CARTDS(NDIMES,NNODLS,*), DVOLUS(*),
     .          EHISTS(NHISTS,*),        ELDISS(NDOFCS,*),
     .          EPMTXS(*),               GPCODS(NDIMES,*),
     .          PROPSS(*),               SHAPES(NNODLS,*),
     .          STRSGS(NSTR1S,*),        ESTIFS(*),
     .          HSTIFS(NEVABS,NNODES),   ELCODS(NDIMES,*)
      DIMENSION BMATXS(NSTR1S,*),        DMATXS(NSTR1S,*),
     .          SIGMAS(*),               XJACMS(NDIMES,*),
     .          VELCMS(*)
      DIMENSION WARTDS(NDIMES,NNODLS,*), TEINIS(NDOFCS,*),
     .          FPCHLS(NFPCH,*),         DVOLIS(*),
     .          AUXS1S(NEVABS,*)
      DIMENSION COEFKS(9)
C
C**** ENTER LOOPS FOR AREA NUMERICAL INTEGRATION
C
      DO IGAUSS=1,NGAULS
C
       IF(NMEMO3S.EQ.0) THEN
C
C**** COMPUTE TEMPERATURE AT GAUSSIAN POINT
C
        DTEMPS=0.0
        TGAUSS=0.0
        TGAUIS=0.0
        PSEUDOS=0.0
        DO INODLS=1,NNODLS
         IF(KDYNAS.EQ.1)
     .   DTEMPS=DTEMPS+SHAPES(INODLS,IGAUSS)*VELCMS(INODLS)
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
        ILAHES=-1
        ISINRS=1
C
        IDECIS=0
        IF(NMEMO10S.EQ.1) IDECIS=1               ! density changes
        IF(IDECIS.EQ.1)
     .   CALL CAPCOFS(BASMM,BASCC,PROPSS,TGAUSS,
     .                DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .                DSOURS,COUTDS, BASMI,TGAUIS,PSEUDOS)
        CALL CONCOFS(COEFKS,PROPSS,TGAUSS,PSEUDOS)
       ELSE
        IX=NDIMETOS-1
        BASMM= EHISTS(1,IGAUSS)
        BASMI= EHISTS(6+IX,IGAUSS)
        DO I=1,NDIMETOS
         COEFKS(I)=EHISTS(3+I-1,IGAUSS)
        ENDDO
       ENDIF                          ! nmemo3s.eq.0
C
C**** COMPUTES CONDUCTIVITY TENSOR    (only isotropic case available !!)
C
       IF(NMEMO10S.EQ.1) THEN
        IF(BASMM.GT.0.0) THEN
         DO I=1,NDIMETOS
          COEFKS(I)=COEFKS(I)*(BASMI/BASMM)
         ENDDO
        ENDIF
       ENDIF
C
       IF(ITERME.GT.0) THEN                     ! bidirectional coupling

        call runends('iterme > 0 not implemented in stif05s')

c       IF(ITERMD.GT.0) THEN                    ! deformed shape
c        IF(LARGET.EQ.1) THEN                   ! TLF
c         DO I=1,NDIMETO
c          detjj=1.0    ! determinant of F (from mechanical computation)
c          BASMM=BASMI/DETJJ
c          COEFKT(I)=COEFKT(I)*(BASMI/BASMM)
c         ENDDO
c        ENDIF
c        IF(LARGET.EQ.2) THEN                   ! ULF
c         DO I=1,NDIMETO
c          COEFKT(I)=COEFKT(I)
c         ENDDO
c        ENDIF
c        IF(LARGET.EQ.3) THEN                   ! Eulerian
c         DO I=1,NDIMETO
c          COEFKT(I)=COEFKT(I)
c          IF(IFREKT.EQ.2)
c    .      COEFKT(I)=COEFKT(I)*DVOLIT(IGAUST)/DVOLUT(IGAUST)
c         ENDDO
c        ENDIF
c       ENDIF
       ENDIF
C
       I=1
       DO IDIMES=1,NDIMES                ! nstr1s=ndimes to be improved
        DO JDIMES=1,NDIMES
         DMATXS(IDIMES,JDIMES)=0.0
         IF(ISOTRS.EQ.0.OR.ISOTRS.EQ.1) THEN
          IF(IDIMES.EQ.JDIMES) THEN
           IF(ISOTRS.EQ.1) I=IDIMES
           DMATXS(IDIMES,JDIMES)=COEFKS(I)
          ENDIF
          IF(ISOTRS.EQ.2) THEN
           CALL RUNENDS('ERROR: FULLY ANIST. MAT. NOT IMPLEMENTED')
          ENDIF
         ENDIF
        ENDDO
       ENDDO
C
C**** CALCULATES THE MATERIAL CONDUCTIVITY MATRIX
C
       NSTRES=NDIMES
       NSTRSS=NDIMES
       IF(LARGES.EQ.1.OR.LARGES.EQ.2)
     .  CALL RUNENDS('ERROR IN STIF05S: LARGE=1,2 NOT IMPLEMENTED')
C
       CALL KMATRIT(BMATXS,        CARTDS(1,1,IGAUSS),  DMATXS,
     .              DVOLUS(IGAUSS),ESTIFS,            GPCODS(1,IGAUSS),
     .              LARGES,KSYMMS, NDIMES,NDOFNS,NEVABS,NKOVAS, NNODLS,
     .              NSTRES,NSTRSS, NTYPES,SHAPES(1,IGAUSS),
     .              WARTDS(1,1,IGAUSS),AUXS1S)
C
      END DO         ! igauss=1,ngauls
C
      RETURN
      END
