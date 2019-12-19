      SUBROUTINE STIC05S(CARTDS,DVOLUS,EHISTS,ELDISS,EPMTXS,GPCODS,
     .                   PROPSS,SHAPES,STRSGS,ESTIFS,HSTIFS,BMATXS,
     .                   DMATXS,SIGMAS,XJACMS,WHAPES,ADVEMS,VELCMS,
     .                   TEINIS,FPCHLS)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE JACOBIAN MATRIX DUE TO CONVECTION
C     EFFECTS
C     ( ELEMENT TYPE NO. 5 )
C
C***********************************************************************
C
C     Index of variables:
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
     .          HSTIFS(NEVABS,NNODES)
      DIMENSION BMATXS(NSTR1S,*),        DMATXS(NSTR1S,*),
     .          SIGMAS(*),               XJACMS(NDIMES,*)
      DIMENSION WHAPES(NNODLS,*),        ADVEMS(NDIMES,*)
      DIMENSION VELCMS(*),               TEINIS(NDOFCS,*),
     .          FPCHLS(NFPCH,*)
      DIMENSION VELOCS(3)
C
C**** ENTER LOOPS FOR AREA NUMERICAL INTEGRATION
C
      DO IGAUSS=1,NGAULS
C
       DO IDIMES=1,NDIMES
        VELOCS(IDIMES)=0.0
        DO INODLS=1,NNODLS
         VELOCS(IDIMES)=VELOCS(IDIMES)+SHAPES(INODLS,IGAUSS)*
     .                  ADVEMS(IDIMES,INODLS)
        ENDDO
       END DO
C
C**** COMPUTE DENSITY*CAPACITY*VELOCITY
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
        ILAHES=0
        ISINRS=1
C
        CALL CAPCOFS(COEFMS,COEFCS,PROPSS,TGAUSS,
     .               DTEMPS,SOUR1S,COEFLS,ILAHES,ISINRS,
     .               DSOURS,COUTDS, BASMI,TGAUIS,PSEUDOS)
       ELSE
        IX=NDIMETOS-1
        COEFMS=EHISTS(1,IGAUSS)          ! COEFMS=Density
        BASMI= EHISTS(6+IX,IGAUSS)       ! BASMI=Initial density
        COEFCS=EHISTS(2,IGAUSS)          ! COEFCS=Capacity
       ENDIF                     ! nmemo3s.eq.0
C
C**** COMPUTES DENSITY
C
       COEFMS=BASMI
C
C**** COMPUTE DENSITY*CAPACITY*VELOCITY
C
       DO IDIMES=1,NDIMES               ! nstr1s=ndimes
        DMATXS(IDIMES,1)=COEFMS*COEFCS*VELOCS(IDIMES)
       ENDDO
C
C**** CALCULATES THE CONVECTIVE MATRIX
C
       IF(LARGES.EQ.1.OR.LARGES.EQ.2)
     .  CALL RUNENDS('ERROR IN STIF05S: LARGE=1,2 NOT IMPLEMENTED')
C
       CALL KMATRIC(BMATXS,        CARTDS(1,1,IGAUSS),  DMATXS,
     .              DVOLUS(IGAUSS),ESTIFS,            GPCODS(1,IGAUSS),
     .              LARGES,KSYMMS, NDIMES,NDOFNS,NEVABS,NKOVAS, NNODLS,
     .              NSTRES,NSTRSS, NTYPES,WHAPES(1,IGAUSS))
C
      END DO ! IGAUS=1,NGAUL
C
      RETURN
      END
