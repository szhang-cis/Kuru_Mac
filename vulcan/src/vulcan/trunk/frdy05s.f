      SUBROUTINE FRDY05S(PROPSS,VELCMS,ELELMS,DVOLUS,SHAPES,EHISTS,
     .                   ELDISS,DISIMS,ELELTS,WHAPES,TEINIS,FPCHLS,
     .                   DVOLIS)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "THERMAL DYNAMIC FORCES"
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
      DIMENSION PROPSS(*),        VELCMS(*),
     .          ELELMS(*),        DVOLUS(*),
     .          SHAPES(NNODLS,*), EHISTS(NHISTS,*),
     .          ELDISS(NDOFCS,*), DISIMS(*),
     .          ELELTS(*)
      DIMENSION WHAPES(NNODLS,*), TEINIS(NDOFCS,*),
     .          FPCHLS(NFPCH,*),  DVOLIS(*)
C
C**** INITIALIZE ELEMENT CONTRIBUTION
C
C     Note: this is only necessary if the transient terms are included
C           in the convergence criterion
C
      LHINCS=0                         ! better as input (see therdyt.f)
      IF(LHINCS.EQ.1) THEN
       DO IEVABS=1,NEVABS
        ELELTS(IEVABS)=0.0D00
       END DO
      ENDIF
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO 100 IGAUSS=1,NGAULS
C
      DDDTPS=0.0
      DTEMPS=0.0
      TGAUSS=0.0
      TGAUIS=0.0
      PSEUDOS=0.0
      DO INODLS=1,NNODLS
       DDDTPS=DDDTPS+SHAPES(INODLS,IGAUSS)*DISIMS(INODLS)
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
C**** UPDATES THE HEAT CAPACITY COEFFICIENT
C
      ILAHES=0
      ISINRS=2
C
      IF(NMEMO3S.EQ.0) THEN
       CALL CAPCOFS(BASMM,BASCC,PROPSS,TGAUSS,
     .              DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .              DSOURS,COUTDS, BASMI,TGAUIS,PSEUDOS)
      ELSE
       IX=NDIMETOS-1
       IF(IMICR.EQ.1) THEN
        BASMM=EHISTS(1,IGAUSS)                   ! computed in micr05s.f
        BASMI=EHISTS(6+IX,IGAUSS)                ! computed in micr05s.f
        BASCC=EHISTS(2,IGAUSS)                   ! computed in micr05s.f
       ELSE
        IF(ICONVS.EQ.0) THEN
         IDECIS=0
         IF(NMEMO10S.EQ.1) IDECIS=1              ! density changes
         IF(IDECIT.EQ.0) THEN
          CALL CAPCOFS(EHISTS(1,IGAUSS),EHISTS(2,IGAUSS),PROPSS,TGAUSS,
     .                 DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .                 DSOURS,EHISTS(7+IX,IGAUSS),EHISTS(6+IX,IGAUSS),
     .                                                   TGAUIS,PSEUDOS)
         ELSE
          CALL CAPCOFS(BASMM,EHISTS(2,IGAUSS),PROPSS,TGAUSS,
     .                 DTEMPS,SOUR1S,SOUR2S,ILAHES,ISINRS,
     .                 DSOURS,EHISTS(7+IX,IGAUSS),BASMI,TGAUIS,PSEUDOS)
         ENDIF
        ENDIF
        BASMM=EHISTS(1,IGAUSS)       ! computed in frin05s.f if IDECIS=1
C                                    ! or in fric05s.f if ICONVS=1
        BASMI=EHISTS(6+IX,IGAUSS)    ! computed in frin05s.f if IDECIS=1
C                                    ! or in fric05s.f if ICONVS=1
        BASCC=EHISTS(2,IGAUSS)       ! computed in fric05s.f if ICONVS=1
       ENDIF
      ENDIF
C
C**** EVALUATE EFFECTIVE HEATS
C
      DESIHS=BASMI*BASCC*DTEMPS
C
      IF(ITERME.GT.0) THEN                      ! bidirectional coupling

       call runends('iterme > 0 not implemented in frdy05s')

c      IF(ITERMD.GT.0) THEN                     ! deformed shape
c       IF(LARGET.EQ.1) THEN                    ! TLF
c        DESIHT=BASMI*BASCC*DTEMPT
c       ENDIF
c       IF(LARGET.EQ.2) THEN                    ! ULF
c        detjj=1.0      ! determinant of F (from mechanical computation)
c        BASMM=BASMI/DETJJ
c        DESIHT=BASMM*BASCC*DTEMPT
c       ENDIF
c       IF(LARGET.EQ.3) THEN                    ! Eulerian
c        BASMM=BASMI*DVOLIT(IGAUST)/DVOLUT(IGAUST)
c        DESIHT=BASMM*BASCC*DTEMPT
c       ENDIF
c      ENDIF
      ENDIF
C
C**** INTEGRATE THE HEATS INTO THE INTERNAL "HEATS FORCES" DUE TO
C     TRANSIENT EFFECTS
C
      DO IEVABS=1,NNODLS
       ELELMS(IEVABS)=ELELMS(IEVABS)+WHAPES(IEVABS,IGAUSS)*DESIHS*
     .                DVOLUS(IGAUSS)
       IF(LHINCS.EQ.1) THEN
        ELELTS(IEVABS)=ELELTS(IEVABS)+SHAPES(IEVABS,IGAUSS)*DESIHS*
     .                 DVOLUS(IGAUSS) 
       ENDIF
      END DO 
C
  100 CONTINUE
C
      RETURN
      END
