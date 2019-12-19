      SUBROUTINE FRIN04(ELDIS,PROPS,BMSIG,TENOD,ELCOD,PWOEL,PREAL,TGAPL,
     .                  VNORL)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESISTING FORCES OF THE NODAL 
C     GAP ELEMENT NO. 4
C
C     Gap:      g_i=n_(1)_i*[u_(1)_i-u_(2)_i]  i=1,...,n_dim
C     Pressure: p_i=H(g_i)*E_n*g_i  i=1,...,n_dim   (H=Heaviside funct.)
C     Force:    F_i=p_i*A_i         i=1,...,n_dim   (A_i=proyected area)
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
      DIMENSION VERSO(3)
C
C**** PROPERTIES FOR THE CONTACT PROBLEM
C
C     PROPS(2)=RIGIN            : NORNAL STIFFNESS
C     PROPS(5)=TEMPL            : LIQUIDUS TEMPERATURE
C     PROPS(6)=INOTE            : DENOTES WHICH NODE BELONGS TO PIECE
C
      RIGIN=PROPS(2)
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
C**** COMPUTE CONTACT PRESSURE
C
       PRESN=RIGIN*DGAPN
       IF(TENOD(INOTE).LT.TEMPL) THEN
        IF(DGAPN.LE.0.0) THEN
         PRESN=0.0
        ENDIF
       ENDIF
C
C**** EVALUATE ELEMENT CONTRIBUTION OF CONTACT FORCE
C     (AS A EXTERNAL FORCE >> -, AS A RESIDUAL FORCE >> +)
C
       BMSIG(IDOFN)      = PROYE*PRESN
       BMSIG(IDOFN+NDOFN)=-PROYE*PRESN
      ENDDO                ! idofn=1,ndofn
C
C**** COMPUTE NORMAL GAP & PRESSURE TO CONVECTION-RADIATION COEFFICIENT
C
      DGAPN=0.0
      DO IDOFN=1,NDOFN
       PROYE=VERSO(IDOFN)
       ELDI1=ELDIS(IDOFN,1)
       ELDI2=ELDIS(IDOFN,2)
       DGAPN=DGAPN+PROYE*(ELDI1-ELDI2)
      END DO
C
      IF(DGAPN.GE.0.0) THEN      ! contact
       PRESN=RIGIN*DGAPN
       DGAPN=0.0
      ELSE                       ! a normal gap is produced
       PRESN=0.0
       DGAPN=-DGAPN
      ENDIF
C
      IF(ITERME.GE.0) THEN       ! uni or bidirectional coupling
       DO INODL=1,NNODL
        TGAPL(INODL)=DGAPN
        PREAL(INODL)=PRESN
       ENDDO
      ENDIF
C
      RETURN
      END
