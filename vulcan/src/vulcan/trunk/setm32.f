      SUBROUTINE SETM32(DVOLU,ELCOD,GPCOD,LNODS,PROPS,SHAPE,THICK,
     .                  DERIV,POSGP,WEIGP,XJACM,VNORL,STRSG,VTANL)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP SOME NEEDED CONSTANT MATRICES FOR
C     FUTURE USE ( ELEMENT NO. 32 )
C
C     Notes: VNORL is computed at each Gauss point (VNORL(NDIME,NGAUS))
C
C            SHAPE & DERIV must be dimensioned with
C            NNODN due to non-coincident contact mesh (NOCOL=1)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/TUNING4/RITEN,RITEF,TRUPL,TRUPM,TOLGA,TOLGAM
C
      COMMON/JACOBSA/IEROR,KEROR
C
      DIMENSION DERIV(NDIML,NNODN,*), DVOLU(*),
     .          ELCOD(NDIME,NNODL),   GPCOD(NDIME,NNODL),
     .          LNODS(*),             PROPS(*),
     .          SHAPE(NNODN,*)
      DIMENSION POSGP(NDIML,*),       WEIGP(*),
     .          XJACM(NDIME,NDIML),   VNORL(NDIME,*)
      DIMENSION STRSG(NSTR1,*),         ! useful for non-coincident mesh
     .          VTANL(NDIME,NDIME-1,*)
C
      DIMENSION EXIDI(3), ETADI(3), PVECT(3,9)
C
      TWOPI=6.283185307179586D0
C
C**** EVALUATE QUANTITIES ASSOCIATED WITH THE GAUSSIAN POINTS
C
      NNOBO=NNODL/2
      ICOLI=INT(PROPS(1))
      IF(ICOLI.EQ.1) THEN         ! contact
       IF(NPOIC.GT.0) THEN
        IF(NNODC.GT.0) NNOBO=NNODL/3
       ENDIF
      ENDIF
      IF(NOCOL.EQ.1) NNOBO=NNODL
C
C**** IDENTIFY SAMPLING POINTS AND WEIGHTS
C
      CALL RULEPW(NDIML,NGAUL,NRULE,POSGP,WEIGP)
C
C**** LOOP ON INTEGRATION POINTS
C
      DO IGAUS=1,NGAUL
C
C**** COMPUTE SHAPE FUNCTIONS AND DERIVATIVES
C
       EXISP=POSGP(1,IGAUS)
       ETASP=0.0D0
       IF(NDIML.EQ.2) ETASP=POSGP(2,IGAUS)
C
       CALL SHAFUN(DERIV(1,1,IGAUS),EXISP,ETASP,EZETA,NDIML,NNOBO,
     .                 1,    0,SHAPE(1,IGAUS))
C
       IF(NOCOL.EQ.1) THEN                         ! non-coincident mesh
        IAUXZ=1
        EXISP=STRSG(IAUXZ,IGAUS)               ! unnecessary for itask=1
        ETASP=0.0D0
        IF(NDIML.EQ.2) ETASP=STRSG(IAUXZ+1,IGAUS)
C
        CALL SHAFUN(DERIV(1,1,IGAUS),EXISP,ETASP,EZETA,NDIML,NNOBO,
     .                  1,NNOBO,SHAPE(1,IGAUS))
       ENDIF
C
C**** COMPUTE JACOBIAN ( NOTE CARTESIAN DERIVATIVES ARE NOT COMPUTED )
C
       CALL JACBOU(DERIV(1,1,IGAUS),DETJM,ELCOD,GPCOD(1,IGAUS),
     .             IELEM,NDIME,NDIML,NNOBO,SHAPE(1,IGAUS),XJACM,
     .             LURES,LUPRI)
       IF(IEROR.NE.0) RETURN
C
C**** INTEGRATION WEIGHTS
C
                        DVOLU(IGAUS)=WEIGP(IGAUS)*DETJM
       IF(NTYPE.EQ.3)   DVOLU(IGAUS)=DVOLU(IGAUS)*TWOPI*GPCOD(1,IGAUS)
       IF(NTYPE.EQ.5)   DVOLU(IGAUS)=DVOLU(IGAUS)*THICK ! 1-D
C
      END DO
C
C**** CHECKS GAUSSIAN VOLUME FOR AXISYMMETRIC PROBLEMS FOR OUTPUT OPER.
C
      IF(NTYPE.EQ.3) THEN
       CALL GAUCEK(DVOLU,   32,    1,NGAUL,NNODL,ICEKE)
       IF(ICEKE.NE.0) CALL RUNEND('ERROR: WRONG NGAUS IN gaucek.f-32')
      ENDIF
C
C**** DETERMINES BETWEEN CONTACT OR LINKING
C
      ICOLI=INT(PROPS(1))
      IF(ICOLI.EQ.1) THEN
C
C**** PROPERTIES FOR THE CONTACT PROBLEM
C
C     PROPS(6)=INOTE:DENOTES WHICH FACE BELONGS TO PIECE
C
       INOTE=INT(PROPS( 6))
       IFRIC=INT(PROPS(51))
C
       NEVBO=NNOBO*NDOFN
C
C**** LOOP OVER INTEGRATION POINTS
C
       DO 100 IGAUS=1,NGAUL
C
C**** COMPUTES OUTWARD UNIT NORMAL TO BODY 1 (should be INOTE ?)
C
C     Note: the outward unit normal vector of a body (1 in this case) is
C           obtained if a clockwise node connectivity is defined for the
C           contact element looking at it from the body (this is valid
C           for 2D & 3D problems). This is equivalent to an
C           anti-clockwise connectivity looking at the body from the
C           other contact body.
C           The normal vector is obtained as a cross product of the
C           tangential vectors, i.e., n_x = t_y x t_z
C           (in 2D, t_z=[0,0,1]).
C
       IF(NDIME.EQ.1) THEN                    ! 1D CASE
        EXIDI(1)=ELCOD(IDIME,1)*DERIV(1,1,IGAUS)
        EXINM=EXIDI(1)*EXIDI(1)
        PNORM=EXINM
        PVECT(1,IGAUS)=-EXIDI(1)                      ! normal to body 2
        PVECT(1,IGAUS)= EXIDI(1)                      ! normal to body 1
C
       ELSE IF(NDIME.EQ.2) THEN               ! 2D CASE
        EXINM=0.0D0
        DO IDIME=1,NDIME
         EXIDI(IDIME)=0.0D0
         DO INODE=1,NNOBO
          EXIDI(IDIME)=EXIDI(IDIME)+
     .                 ELCOD(IDIME,INODE)*DERIV(1,INODE,IGAUS)
         ENDDO
         EXINM=EXINM+EXIDI(IDIME)*EXIDI(IDIME)
        ENDDO
        PNORM=EXINM
        PVECT(1,IGAUS)= EXIDI(2)                      ! normal to body 2
        PVECT(2,IGAUS)=-EXIDI(1)
        PVECT(1,IGAUS)=-EXIDI(2)                      ! normal to body 1
        PVECT(2,IGAUS)= EXIDI(1)
C
        IF(IFRIC.EQ.1) THEN
         PNOR1=DSQRT(PNORM)
         VTANL(1,1,IGAUS)=-EXIDI(1)/PNOR1
         VTANL(2,1,IGAUS)=-EXIDI(2)/PNOR1
        ENDIF
C
       ELSE                                   ! 3D CASE
        DO IDIME=1,NDIME            ! compute two vectors on the surface
         EXIDI(IDIME)=0.0D0
         ETADI(IDIME)=0.0D0
         DO INODE=1,NNOBO
          EXIDI(IDIME)=EXIDI(IDIME)+DERIV(1,INODE,IGAUS)*
     .                 ELCOD(IDIME,INODE)
          ETADI(IDIME)=ETADI(IDIME)+DERIV(2,INODE,IGAUS)*
     .                 ELCOD(IDIME,INODE)
         ENDDO
        ENDDO
        PNORM=0.0D0
        DO IDIMA=1,NDIME            ! compute normal as cross product
         IDIMB=IDIMA+1-IDIMA/3*3
         IDIMC=IDIMB+1-IDIMB/3*3
         PVECT(IDIMA,IGAUS)=EXIDI(IDIMB)*ETADI(IDIMC)-
     .                      ETADI(IDIMB)*EXIDI(IDIMC) ! normal to body 1
         PNORM=PNORM+PVECT(IDIMA,IGAUS)*PVECT(IDIMA,IGAUS)
        ENDDO
C
        IF(IFRIC.EQ.1) THEN
         PNOR1=0.0D0
         PNOR2=0.0D0
         DO IDIME=1,NDIME
          PNOR1=PNOR1+EXIDI(IDIME)*EXIDI(IDIME)
          PNOR2=PNOR2+ETADI(IDIME)*ETADI(IDIME)
         ENDDO
         PNOR1=DSQRT(PNOR1)
         PNOR2=DSQRT(PNOR2)
         DO IDIME=1,NDIME
          VTANL(IDIME,1,IGAUS)=-EXIDI(IDIME)/PNOR1
          VTANL(IDIME,2,IGAUS)= ETADI(IDIME)/PNOR2
         ENDDO
        ENDIF
C
       ENDIF
C
       PNORM=DSQRT(PNORM)           ! compute normalized normal
       DO IDIME=1,NDIME
        PVECT(IDIME,IGAUS)=PVECT(IDIME,IGAUS)/PNORM
       ENDDO
C
  100  CONTINUE
C
C**** COMPUTES VNORL
C
       DO IDIME=1,NDIME
        DO IGAUS=1,NGAUL
         VNORL(IDIME,IGAUS)=PVECT(IDIME,IGAUS)
        ENDDO
       ENDDO
C
      ELSE
C
C**** LINKING PROBLEM
C
       ISIPO=0      ! better as input - 0: side to side; 1: node to node
       IF(ISIPO.EQ.1) THEN
        DO IGAUS=1,NGAUL
         DVOLU(IGAUS)=1.0D+00
        ENDDO
       ENDIF
C
      ENDIF             ! icoli.eq.1
C
      RETURN
      END
