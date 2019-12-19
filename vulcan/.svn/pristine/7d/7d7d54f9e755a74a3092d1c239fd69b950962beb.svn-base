      SUBROUTINE FRIN04Y(ELDIS,PROPS,BMSIG,TENOD,ELCOD,PWOEL,PREAL,
     .                   TGAPL,VNORL)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESISTING FORCES OF THE NODAL 
C     GAP ELEMENT NO. 4
C
C     AUGMENTED LAGRANGIAN METHOD (ALM)
C
C     VERSION 0: Agelet's thesis (IAUGM=0)
C     VERSION 1: modified Agelet's thesis (IAUGM=1)
C
C
C     Gap:      g_i=n_(1)_i*[u_(1)_i-u_(2)_i]  i=1,...,n_dim
C     Pressure: p_i=lambda_i                   i=1,...,n_dim
C     Force:    F_i=lambda_i*A_i    i=1,...,n_dim   (A_i=proyected area)
C
C     The pressure p_i is always negative in contact conditions.
C
C
C     VERSION 0:
C
C     The residual of the impenetrability condition can be written in
C     two forms:
C  
C     1) R=-g_i+p_i/E_n     (imodel=0)
C     2) R= g_i-p_i/E_n     (imodel=1)
C
C     The case 1) leads to a non definite positive symmetric jacobian
C     matrix. The case 2) gives a positive definite unsymmetric matrix. 
C
C
C     VERSION 1:
C
C     The residual of the impenetrability condition can be written in
C     two forms:
C
C     1) R=-g_i             (imodel=0)
C     2) R= g_i             (imodel=1)
C
C     The jacobian is assumed to be the same as the one for VERSION 0.
C
C
C     This formulation is almost equivalent to the classical theory:
C
C     Gap:      g_n=n_(1)*[u_(1)-u_(2)] = \sum g_i         i=1,...,n_dim
C     Pressure: p_n=lambda
C     Force:    F_n=H(g_n)*p_n*n_(1)               H: Heaviside function
C            => F_i=H(p_n)*p_n*n_(1)_i                     i=1,...,n_dim
C
C     neglecting the contribution of each g_j for all j ne i in the
C     force vector F_n and replacing H(g_n) by H(g_i) in F_i
C
C
C**** EXACT KUHN-TUCKER CONDITIONS  (they do not work for IAUGM=0)
C
C     IF(DGAPN.GT.0.0) THEN             ! gap closed (contact)
C      DGAPN=0.0
C     ENDIF
C
C     IF(PRESN.GT.0.0) THEN             ! positive pressure
C      PRESN=0.0
C     ENDIF
C
C     IF(DGAPN.LT.0.0) THEN             ! gap open
C      PRESN=0.0
C     ENDIF
C     IF(PRESN.LT.0.0) THEN             ! negative pressure (contact)
C      DGAPN=0.0
C     ENDIF
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
C     PROPS(3)=IMODEL           : CONTACT MODEL
C     PROPS(4)=IMODEX           : CONTACT MODEL FOR REGULARIZATION
C     PROPS(5)=TEMPL            : LIQUIDUS TEMPERATURE
C     PROPS(6)=INOTE            : DENOTES WHICH NODE BELONGS TO PIECE
C
      RIGIN=PROPS(2)
      IMODEL=INT(PROPS(3))  ! 0: symmetric matrix; 1: unsymmetric matrix
      IMODEX=INT(PROPS(4))  ! 0: regularization 1; 1: regularization 2
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
       PRESN=ELDIS(IDOFN,3)
       presz=presn
C
C**** APPROXIMATED KUHN-TUCKER CONDITIONS
C
       IF(TENOD(INOTE).LT.TEMPL) THEN
        if(imodex.eq.0) then
         IF(DGAPN.LE.0.0) THEN                  ! gap open
          DGAPN=0.0
         ENDIF
         IF(PRESN.GE.0.0) THEN                  ! positive pressure
          PRESN=0.0
          presz=0.0
         ENDIF
        else
c         IF(DGAPN.LE.0.0) THEN                  ! gap open (LT)
c          DGAPN=0.0
c          PRESN=0.0
c         ENDIF
         IF(PRESN.GE.0.0) THEN                  ! positive pressure (LT)
c          DGAPN=0.0
c          PRESN=0.0
          presz=0.0
         ENDIF
        endif
       ENDIF
C
C**** EVALUATE ELEMENT CONTRIBUTION OF CONTACT FORCE
C     (AS A EXTERNAL FORCE >> -, AS A RESIDUAL FORCE >> +)
C
c       BMSIG(IDOFN)        =-PROYE*PRESN
c       BMSIG(IDOFN+NDOFN)  = PROYE*PRESN
       BMSIG(IDOFN)        =-PROYE*PRESZ
       BMSIG(IDOFN+NDOFN)  = PROYE*PRESZ
C
       IF(IAUGM.EQ.0) THEN
        if(imodel.eq.0) then
         BMSIG(IDOFN+2*NDOFN)=-DGAPN*RIGIN-PRESN
        else
         BMSIG(IDOFN+2*NDOFN)= DGAPN*RIGIN+PRESN/RIGIN
        endif
       ELSE
        if(imodel.eq.0) then
         BMSIG(IDOFN+2*NDOFN)=-DGAPN
        else
         BMSIG(IDOFN+2*NDOFN)= DGAPN
        endif
       ENDIF
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
       PRESN=PRESN+PROYE*ELDIS(IDOFN,3)
      END DO
C
      IF(DGAPN.GE.0.0) THEN  ! contact
       DGAPN=0.0
       PRESN=-PRESN          ! absolute value (sign of gn ne sign of pn)
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
