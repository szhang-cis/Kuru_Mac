      SUBROUTINE STIF04Z(PROPS,ESTIF,TENOD,ELDIS,ELCOD,DISPL,VNORL)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE STIFFNESS MATRIX FOR THE NODAL GAP 
C     ELEMENT NO. 4
C
C     AUGMENTED LAGRANGIAN METHOD (ALM): VERSION 0 & 1
C
C**** PROPERTIES FOR THE CONTACT PROBLEM
C
C     PROPS(2)=RIGIN            : NORNAL STIFFNESS
C     PROPS(3)=IMODEL           : CONTACT MODEL
C     PROPS(4)=IMODEX           : CONTACT MODEL FOR REGULARIZATION
C     PROPS(5)=TEMPL            : LIQUIDUS TEMPERATURE
C     PROPS(6)=INOTE            : DENOTES WHICH NODE BELONGS TO PIECE
C
C**** STABILIZATION PARAMETERS (CONTROL AT MATERIAL PROPERTIES AND/OR
C     INTERVAL DATA LEVEL) COMING FROM A COMMON OF TUNING PARAMETERS
C
C     IITEN                     : INDEX FOR TUNING PARAMETERS
C     IITEF                     : ITERATION FROM WHICH TRUPL ACTS
C     TRUPL                     : FIRST STIFFNESS IN TENSION
C     TRUPM                     : SECOND STIFFNESS IN TENSION
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      COMMON/TUNING4/RITEN,RITEF,TRUPL,TRUPM,TOLGA,TOLGAM
C
      DIMENSION PROPS(NPROP),       ESTIF(NKOVA)
      DIMENSION TENOD(NNODL)
      DIMENSION ELDIS(NDOFN,NNODL), ELCOD(NDIME,NNODL)
      DIMENSION DISPL(NDOFN,*),     VNORL(NDIME,NNODL)
      DIMENSION AUXS1(3,3),         AUXS2(9,9), AUXS3(3,3),
     .          VERSO(3)
      DIMENSION RIGIN(3)
      DIMENSION RIGIX(3), RIGI3(3)
      DIMENSION RIGIY(3)
      DIMENSION RIGIZ(3)
C
      DO IDOFN=1,NDOFN
       RIGIN(IDOFN)=PROPS(2)
c       RIGIX(IDOFN)=1.0                    ! displacements-contact force
       RIGIX(IDOFN)=RIGIN(IDOFN)            ! displacements-contact force
       RIGI3(IDOFN)=1.0                    ! displacements-contact force
       RIGIY(IDOFN)=RIGIN(IDOFN)           ! displacements-displacements
c       RIGIZ(IDOFN)=RIGIN(IDOFN)           ! contact force-contact force
       RIGIZ(IDOFN)=1.0                     ! contact force-contact force
      ENDDO
      IMODEL=INT(PROPS(3))  ! 0: symmetric matrix; 1: unsymmetric matrix
      IMODEX=INT(PROPS(4))  ! 0: regularization 1; 1: regularization 2
      TEMPL=PROPS(5)
      INOTE=INT(PROPS(6))
C
      IITEF=INT(RITEF)
      IITEN=INT(RITEN)      ! 0: tun. not considered; 1: tun. considered
C
C**** REDEFINES OR CHECKS TRUPM FOR REGULARIZATION 1
C
      if(imodex.eq.0) then
       if(iiten.eq.0) then
        TRUPM=TRUPL         ! default value (10.0E-06)
       else
        IF(TRUPM.EQ.0.0) THEN
         TRUPM=10.E-06
         CALL RUNMEN('WARNING: TRUPM IS SET TO 10.E-06')
        ENDIF
       endif
      endif
C
C**** IF "JACOBIAN_REGULARIZATION_FACTOR" IS USED, REDEFINES TUNING PAR.
C
      NNN20=20
      INP10=INT(PROPS(10))
      IF(INP10.EQ.1) THEN
       RITEF=10.0
       TRUPL=PROPS(11)
       if(imodex.eq.0) then
        IF(TRUPL.LT.10.E-06) TRUPL=10.E-06
        TRUPM=TRUPL
       else
        TRUPM=0.0
       endif
       TOLGA=-1.0E-08
       NNN20=INT(PROPS(12))
      ENDIF
C
C**** OUTWARD UNIT NORMAL TO BODY 1
C
      TOLNOR=-1.0E-08       ! to avoid unrealistic situations
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
C**** INITIALISES AUXS1 & A PART OF AUXS2
C
      DO IDOFN=1,NDOFN
       DO JDOFN=1,NDOFN
        AUXS1(IDOFN,JDOFN)=0.0
        AUXS3(IDOFN,JDOFN)=0.0
       ENDDO
      ENDDO
      DO IDOFN=1,NDOFN*3
       DO JDOFN=1,NDOFN*3
        AUXS2(IDOFN,JDOFN)=0.0
       ENDDO
      ENDDO
C
C**** COMPUTES CONTACT RIGIDITY
C
      IF(TENOD(INOTE).LT.TEMPL) THEN
       DO IDOFN=1,NDOFN
C
C**** COMPUTES GAP FOR NDIME DIRECTIONS AT THE LAST CONVERGED STEP
C
        PROYE=VERSO(IDOFN)
        ELDI1=ELDIS(IDOFN,1)-DISPL(IDOFN,1)
        ELDI2=ELDIS(IDOFN,2)-DISPL(IDOFN,2)
        DGAPA=PROYE*(ELDI1-ELDI2)
C
C**** COMPUTES GAP FOR NDIME DIRECTIONS AT THE CURRENT STEP
C
        ELDI1=ELDIS(IDOFN,1)
        ELDI2=ELDIS(IDOFN,2)
        DGAPB=PROYE*(ELDI1-ELDI2)
C
C**** COMPUTES NORMAL PRESSURE
C
        PRESN=ELDIS(IDOFN,3)
        PRESA=ELDIS(IDOFN,3)-DISPL(IDOFN,3)
C
C**** "STABILIZED" CONTACT RIGIDITY
C
c        IF(IAUGM.EQ.0) THEN
cc         IF(DGAPA.LE.TOLGA.AND.DGAPB.LE.TOLGA) THEN
c         IF(DGAPA.LE.TOLGA.AND.DGAPB.LE.0.0) THEN
c          if(imodex.eq.0) then
c           RIGIZ(IDOFN)=RIGIZ(IDOFN)*TRUPM
c          else
c           RIGIX(IDOFN)=RIGIX(IDOFN)*TRUPM
c           RIGIY(IDOFN)=RIGIY(IDOFN)*TRUPM        ! not used
c          endif
c         ENDIF
c        ENDIF
C
C**** "IMPROVED" CONTACT RIGIDITY
C
c        IF(DGAPB.LE.TOLGA.OR.PRESN.GT.0.0) THEN
cc        IF(DGAPB.LE.0.0) THEN
c         DO III20=1,NNN20
c          IF(IITER.GT.(III20*IITEF)) THEN
c           if(imodex.eq.0) then
c            RIGIZ(IDOFN)=RIGIZ(IDOFN)*TRUPL
c           else
c            RIGIX(IDOFN)=RIGIX(IDOFN)*TRUPL
c            RIGIY(IDOFN)=RIGIY(IDOFN)*TRUPL       ! not used
c           endif
c          ENDIF
c         ENDDO
c        ENDIF



        idiego=1                  !!!!!! input
        if(idiego.eq.1) then

        prest=-rigin(idofn)*tolga
        if(prest.lt.0.0) prest=0.0
C
C**** "STABILIZED" CONTACT RIGIDITY
C
c         IF(PRESA.GT.0.0.AND.PRESN.GT.0.0) THEN
         IF(PRESA.GT.prest.AND.PRESN.GT.prest) THEN
          RIGI3(IDOFN)=RIGI3(IDOFN)*TRUPM
         ENDIF
C
C**** "IMPROVED" CONTACT RIGIDITY
C
c        IF(PRESN.GT.0.0) THEN
        IF(PRESN.GT.prest) THEN
         DO III20=1,NNN20
          IF(IITER.GT.(III20*IITEF)) THEN
           RIGI3(IDOFN)=RIGI3(IDOFN)*TRUPL
          ENDIF
         ENDDO
        ENDIF


c first approach
c        IF(PRESN.GT.0.0) THEN
c         RIGI3(IDOFN)=RIGI3(IDOFN)*TRUPL
c        ENDIF

       endif


       if(idiego.eq.2) then

        prest=rigin(idofn)*tolga
        IF(PRESA.LE.prest.AND.PRESN.LE.prest) THEN
         RIGI3(IDOFN)=RIGI3(IDOFN)
        ENDIF
        IF(PRESA.GT.prest.AND.PRESN.LE.prest) THEN
         RIGI3(IDOFN)=RIGI3(IDOFN)+PRESN/(PRESN-PRESA)
c         RIGI3(IDOFN)=PRESN/(PRESN-PRESA)
c         RIGI3(IDOFN)=RIGI3(IDOFN)
        ENDIF
        IF(PRESA.GT.prest.AND.PRESN.GT.prest) THEN
         RIGI3(IDOFN)=0.0
        ENDIF
        IF(PRESA.LE.prest.AND.PRESN.GT.prest) THEN
c         RIGI3(IDOFN)=PRESN/(PRESN-PRESA)
c         RIGI3(IDOFN)=-PRESA/(PRESN-PRESA)
         RIGI3(IDOFN)=RIGI3(IDOFN)-PRESA/(PRESN-PRESA)
        ENDIF


       endif


       if(idiego.eq.3) then

        prest=rigin(idofn)*tolga
C
C**** "STABILIZED" CONTACT RIGIDITY
C
         IF(PRESA.GT.0.0.AND.PRESN.GT.0.0) THEN
          RIGI3(IDOFN)=RIGI3(IDOFN)*TRUPM
         ENDIF
C
C**** "IMPROVED" CONTACT RIGIDITY
C
         DO III20=1,NNN20
          IF(PRESN.GT.(III20*PREST)) THEN
           RIGI3(IDOFN)=RIGI3(IDOFN)*TRUPL
          ENDIF
         ENDDO


       endif


       if(idiego.eq.4) then

        prest=rigin(idofn)*tolga
C
C**** "STABILIZED" CONTACT RIGIDITY
C
         IF(PRESN.GT.(IITEF*prest)) THEN
          RIGI3(IDOFN)=RIGI3(IDOFN)*TRUPL

c         write(7,*) 'presn > iitef*prest --- ielem',ielem
c         write(7,*) presn

         ENDIF
C
C**** "IMPROVED" CONTACT RIGIDITY
C
         IF(presn.gt.0.0.and.PRESN.LT.(IITEF*prest)) THEN
          RIGI3(IDOFN)=(-1.0/(iitef*prest))*PRESN+RIGI3(IDOFN)+
     .                 RIGI3(IDOFN)*TRUPL

c         write(7,*) 'presn > 0 < iitef*prest --- ielem',ielem
c         write(7,*) presn

         ENDIF

         IF(PRESA.GT.0.0.AND.PRESN.GT.0.0) THEN
          RIGI3(IDOFN)=RIGI3(IDOFN)*TRUPM
         ENDIF
C


       endif


        if(idiego.eq.5) then

        prest=rigin(idofn)*tolga
C
C**** "STABILIZED" CONTACT RIGIDITY
C
c         IF(PRESA.GT.0.0.AND.PRESN.GT.0.0) THEN
         IF(PRESA.GT.prest.AND.PRESN.GT.prest) THEN
          RIGI3(IDOFN)=RIGI3(IDOFN)*TRUPM
         ENDIF
C
C**** "IMPROVED" CONTACT RIGIDITY
C
c        IF(PRESN.GT.0.0) THEN
        IF(PRESN.GT.prest) THEN

         if(iiter.gt.iitef) rigi3(idofn)=rigi3(idofn)+
     .               (0.0-rigi3(idofn))/(miter-iitef)*(iiter-iitef)

         IF(IITER.GT.IITEF) THEN
          RIGI3(IDOFN)=RIGI3(IDOFN)*TRUPL
         ENDIF
        ENDIF


       endif


C
       END DO     !idofn=1,ndofn
      ENDIF
C
C**** COMPUTES ONLY A FORTH PART OF ELEMENTAL CONTACT MATRIX
C     (RELATES TO: DISPLACEMENTS-CONTACT FORCE)
C
      DO IDOFN=1,NDOFN
       PROYE=VERSO(IDOFN)
       DO JDOFN=1,NDOFN
        IF(IDOFN.EQ.JDOFN)
     .   AUXS1(IDOFN,JDOFN)=AUXS1(IDOFN,JDOFN)+PROYE*RIGIX(IDOFN)
       ENDDO
      ENDDO

      DO IDOFN=1,NDOFN
       PROYE=VERSO(IDOFN)
       DO JDOFN=1,NDOFN
        IF(IDOFN.EQ.JDOFN)
     .   AUXS3(IDOFN,JDOFN)=AUXS3(IDOFN,JDOFN)+PROYE*RIGI3(IDOFN)
       ENDDO
      ENDDO
C
C**** ELEMENTAL CONTACT MATRIX
C     (RELATES TO: DISPLACEMENTS-CONTACT FORCE)
C
      DO IDOFN=1,NDOFN
       DO JDOFN=1,NDOFN
        if(imodel.eq.0) then
         AUXS2(IDOFN,JDOFN+2*NDOFN)=-AUXS1(IDOFN,JDOFN)
         AUXS2(IDOFN+NDOFN,JDOFN+2*NDOFN)=AUXS1(IDOFN,JDOFN)
        else
c         AUXS2(IDOFN,JDOFN+2*NDOFN)=-AUXS1(IDOFN,JDOFN)
c         AUXS2(IDOFN+NDOFN,JDOFN+2*NDOFN)=AUXS1(IDOFN,JDOFN)
         AUXS2(IDOFN,JDOFN+2*NDOFN)=-AUXS3(IDOFN,JDOFN)
         AUXS2(IDOFN+NDOFN,JDOFN+2*NDOFN)=AUXS3(IDOFN,JDOFN)

         AUXS2(IDOFN+2*NDOFN,JDOFN)=AUXS1(IDOFN,JDOFN)
         AUXS2(IDOFN+2*NDOFN,JDOFN+NDOFN)=-AUXS1(IDOFN,JDOFN)

c  como si fuera imodel=0
c         AUXS2(IDOFN+2*NDOFN,JDOFN)=-AUXS1(IDOFN,JDOFN)
c         AUXS2(IDOFN+2*NDOFN,JDOFN+NDOFN)=AUXS1(IDOFN,JDOFN)

        endif
       ENDDO
      ENDDO
C
C**** ELEMENTAL CONTACT MATRIX
C     (RELATES TO: CONTACT FORCE-CONTACT FORCE)
C
      DO IDOFN=1,NDOFN
       DO JDOFN=IDOFN,NDOFN
        if(imodel.eq.0) then
         IF(IDOFN.EQ.JDOFN)
     .    AUXS2(IDOFN+2*NDOFN,JDOFN+2*NDOFN)=-1.0/RIGIZ(IDOFN)
        else
         IF(IDOFN.EQ.JDOFN)
     .    AUXS2(IDOFN+2*NDOFN,JDOFN+2*NDOFN)=1.0/RIGIZ(IDOFN)
        endif
       ENDDO
      ENDDO
C
C**** EVALUATE ELEMENT CONTRIBUTION
C
      if(imodel.eq.0) then
       DO IEVAB=1,NDOFN*3
        DO JEVAB=IEVAB,NDOFN*3
         ICONB=(2*NEVAB-IEVAB)*(IEVAB-1)/2+JEVAB
         IF(KSYMM.EQ.0) ICONB=(IEVAB-1)*NEVAB+JEVAB  ! unsymmetric
         ESTIF(ICONB)=AUXS2(IEVAB,JEVAB)
        ENDDO
       ENDDO
C
C**** LOAD ESTIF IN A SQUARE FORM FOR UNSYMMETRIC SOLVER
C
       IF(KSYMM.EQ.0) THEN
        DO IEVAB=1,NEVAB
         DO JEVAB=IEVAB,NEVAB
          KLOCS=(IEVAB-1)*NEVAB+JEVAB
          KLOCI=(JEVAB-1)*NEVAB+IEVAB
          ESTIF(KLOCI)=ESTIF(KLOCS)
         END DO
        END DO
       END IF
      else
       ICONB=0
c       DO IEVAB=1,NDOFN*3
c        DO JEVAB=1,NDOFN*3
       DO IEVAB=1,NEVAB
        DO JEVAB=1,NEVAB
         ICONB=ICONB+1
c         ESTIF(ICONB)=AUXS2(IEVAB,JEVAB)

         if(ievab.le.(3*ndofn).and.jevab.le.(3*ndofn)) then
          ESTIF(ICONB)=AUXS2(IEVAB,JEVAB)
         else
          ESTIF(ICONB)=0.0
         endif

        ENDDO
       ENDDO
      endif
C
      RETURN
      END
