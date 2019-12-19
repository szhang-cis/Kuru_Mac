      SUBROUTINE FRIN04C(ELDIS,PROPS,BMSIG,TENOD,ELCOD,PWOEL,PREAL,
     .                   TGAPL,VNORL,DISIL)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESISTING FORCES OF THE NODAL 
C     GAP ELEMENT NO. 4
C
C     AUGMENTED LAGRANGIAN METHOD (ALM): 
C
C     VERSION 2 (modified Cante's thesis)       IAUGM=2
C     VERSION 3 (Cante's thesis)                IAUGM=3
C
C     Gap:      g_i=n_(1)_i*[u_(1)_i-u_(2)_i]  i=1,...,n_dim
C     Pressure: p_i=lambda_i                   i=1,...,n_dim
C     Force:    F_i=lambda_i*A_i    i=1,...,n_dim   (A_i=proyected area)
C
C     The pressure p_i is always negative in contact conditions.
C
C     The residual of the impenetrability condition can be written as:
C  
C     R=-g                              (IAUGM=2 & 3)
C
C
C     the following to be revised:
C
C     This formulation is almost equivalent to the classical theory:
C
C     Gap:      g_n=n_(1)*[u_(1)-u_(2)] = \sum g_i         i=1,...,n_dim
C     Pressure: p_n=H(g_n)*E_n*g_n                  (H=Heaviside funct.)
C     Force:    F_n=p_n*n_(1)
C            => F_i=p_n*n_(1)_i                            i=1,...,n_dim
C
C     neglecting the contribution of each g_j for all j ne i in the
C     force vector F_n and replacing H(g_n) by H(g_i) in F_i
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION ELDIS(NDOFN,NNODL), PROPS(NPROP),
     .          BMSIG(NDOFN*NNODL), TENOD(NNODL),
     .          ELCOD(NDIME,NNODL)
      DIMENSION PWOEL(NNODL),       PREAL(NNODL),
     .          TGAPL(NNODL),       VNORL(NDIME,NNODL)
      DIMENSION VERSO(3),           DISIL(NDOFN,*)
C
C**** PROPERTIES FOR THE CONTACT PROBLEM
C
C     PROPS(2)=RIGIN            : NORNAL STIFFNESS
C     PROPS(3)=IMODEL           : CONTACT MODEL
C     PROPS(5)=TEMPL            : LIQUIDUS TEMPERATURE
C     PROPS(6)=INOTE            : DENOTES WHICH NODE BELONGS TO PIECE
C
      RIGIN=PROPS(2)
      IMODEL=INT(PROPS(3))         ! not used
      TEMPL=PROPS(5)
      INOTE=INT(PROPS(6))
C
C**** OUTWARD UNIT NORMAL TO BODY 1
C
      TOLNOR=-1.0E-08
      VNORC=0.0
      DO IDIME=1,NDIME
       VNORA=VNORL(IDIME,1)*VNORL(IDIME,2)
       IF(VNORA.GT.TOLNOR) VNORL(IDIME,1)=0.0
       VNORA=VNORL(IDIME,1)
       VNORC=VNORC+VNORA*VNORA
      ENDDO
      VNORC=DSQRT(VNORC)
      DO IDIME=1,NDIME
       VNORA=VNORL(IDIME,1)
       VERSO(IDIME)=VNORA/VNORC
      ENDDO
C
C**** INITIALISE MECHANICAL COUPLING TERM
C
      IF(ITERME.GT.0) THEN       ! bidirectional coupling
       DO INODL=1,NNODL
        PWOEL(INODL)=0.0
       END DO
      ENDIF
C
C**** COMPUTE GAP FOR NDIME DIRECTIONS
C
C     DGAPN  > 0: CONTACT PROBLEM
C     DGAPN <= 0: NO CONTACT PROBLEM
C
      DO IDOFN=1,NDOFN
       PROYE=VERSO(IDOFN)
       ELDI1=ELDIS(IDOFN,1)
       ELDI2=ELDIS(IDOFN,2)
       DGAPN=PROYE*(ELDI1-ELDI2)
C
C**** COMPUTES CONTACT PRESSURE
C
       IF(IAUGM.EQ.2) THEN
        DISIL(IDOFN,3)=-RIGIN*DGAPN
       ENDIF
       IF(IAUGM.EQ.3) THEN
        IF(IAUG3.EQ.0) THEN
         DISIL(IDOFN,3)=-RIGIN*DGAPN
        ELSE
         DISIL(IDOFN,3)=0.0
        ENDIF
       ENDIF
C
       PRESN=ELDIS(IDOFN,3)+DISIL(IDOFN,3)
C
C**** KUHN-TUCKER CONDITION
C
       IF(TENOD(INOTE).LT.TEMPL) THEN
        IF(PRESN.GT.0.0) THEN                  ! positive pressure
         PRESN=0.0
         DISIL(IDOFN,3)=-ELDIS(IDOFN,3)        ! to obtain zero pressure
        ENDIF
       ENDIF
C
C**** EVALUATE ELEMENT CONTRIBUTION OF CONTACT FORCE
C     (AS A EXTERNAL FORCE >> -, AS A RESIDUAL FORCE >> +)
C
       BMSIG(IDOFN)        =-PROYE*PRESN
       BMSIG(IDOFN+NDOFN)  = PROYE*PRESN
      ENDDO                ! idofn=1,ndofn
C
C**** COMPUTE NORMAL GAP & PRESSURE TO CONVECTION-RADIATION COEFFICIENT
C
      DGAPN=0.0
      PRESN=0.0
      DO IDOFN=1,NDOFN
       PROYE=VERSO(IDOFN)
       ELDI1=ELDIS(IDOFN,1)
       ELDI2=ELDIS(IDOFN,2)
       DGAPN=DGAPN+PROYE*(ELDI1-ELDI2)
       PRESN=PRESN+ELDIS(IDOFN,3)-RIGIN*DGAPN       ! ?????
      END DO
C
      IF(PRESN.LE.0.0) THEN  ! contact
       DGAPN=0.0
       PRESN=-PRESN          ! absolute value
      ELSE                   ! a gap is produced
       PRESN=0.0 
       DGAPN=-DGAPN          ! absolute value
      ENDIF
C
      IF(ITERME.GE.0) THEN   ! uni or bidirectional coupling
       DO INODL=1,NNODL
        TGAPL(INODL)=DGAPN
        PREAL(INODL)=PRESN
       ENDDO
      ENDIF
C
      RETURN
      END
