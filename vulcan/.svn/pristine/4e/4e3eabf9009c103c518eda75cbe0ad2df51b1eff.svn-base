      SUBROUTINE POINTE
C***********************************************************************
C
C**** THIS ROUTINE DETERMINES POINTERS OF EHIST ARRAY FOR THE INTERNAL
C     VARIABLES
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
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
C**** INCOMPRESSIBILITY CONDITION
C
      IF(KPROB.EQ.5) RETURN
C
C**** DETERMINES (ACCORDING TO KPLAi):
C
C     1) NUMBER OF COMPONENTS OF THE INTERNAL VARIABLES WITHOUT PLASTIC
C        DEFORMATION (NBASE); see frin30.f
C     2) POINTERS OF INTERNAL VARIABLES (IPLAS); see frin30.f & outg30.f
C     3) POINTERS OF OUTPUT VECTORS OF INTERNAL VARIABLES
C        (IPLAO(I), I=1,NNUIN); see smoi30.f
C     4) POINTERS OF OUTPUT NODAL VECTORS OF INTERNAL VARIABLES
C        (from IPLAN(I) to IPLAM(I), I=1,NNUNO); see outnod.f
C        Note: IPLAN(I) is always 1 => it could be removed!
C     5) NUMBER OF COMPONENTS OF THE NUMBER OF INTERNAL VARIABLES TO
C        BE PRINTED IN THE POST FILE (NNUIN); see outnod.f
C     6) NUMBER OF INTERNAL VARIABLES TO BE SMOOTHED AND PRINTED IN
C        RESULTS FILE (NNUNO); see smoi30.f & outnod.f
C
C     NOTES:
C
C     a) THE COMPONENTS OF IPLAO AND IPLAN MUST BE EQUAL
C     b) NNUIN MUST BE EQUAL TO \sum_{i=1,nnuno} IPLAN(i) (IPLAM(i))
C     c) FOR SCALARS INTERNAL VARIABLES TO PRINT, IPLAN(I)=1 AND
C        IPLAM(I)=1 FOR ALL I (NNUIN=NNUNO)
C     d) PARAMETERS KPLAi DEFINES THE TYPE OF CONSTITUTIVE MODEL => 
C        THE FOLLOWING VARIABLES DEPEND ON KPLAI
C     e) IVLAS, IVLAO, IVLAN & IVLAM not used now
C
      NBASE=1                                     ! eff. plast. strain
      IF(KPLA1.EQ.1) NBASE=NBASE+1                ! isotropic hardening
      IF(KPLA2.EQ.1) NBASE=NBASE+NSTR1            ! kinematic hardening
C
      NNDAM=0
c     IF(KPLA3.EQ.1) NNDAM=1                      ! damage (standard)
      IF(KPLA3.EQ.1) NNDAM=2                      ! damage (concrete)
C
      IF(KPLA4.EQ.1) NBASE=NBASE+2                ! porosity+total hard.
C
      NSHRI=0
      IF(KPLA5.EQ.1) NSHRI=1                      ! shrinkage
C
      NFATI=0
      IF(KPLA6.EQ.1) NFATI=8                      ! fatigue
C
      IF(KPLA7.EQ.1) NBASE=NBASE+1                ! dot e^p
C
      IF(KPLA8.EQ.1) NBASE=NBASE+2                ! recrystall. (S & L)
C
      IF(KPLA9.EQ.1) NBASE=NBASE+1                ! old equiv. stress
C
      IF(KPLA10.EQ.1) NBASE=NBASE+(1+NCHAI)*NSTR1 ! viscous stress
C
C     IF(KPLA11.EQ.1) NBASE=...                   ! not used now!
C
      NDUAL=0
      NDUAX=2   ! int. variables; check dimensions of EBASEi in cepl53.f
      IF(KPLA12.EQ.1) NDUAL=2*NSTR1+2*NKOST+      ! e_p,Const tensor,
     .                      2*NSTR1+2*NDUAX+      ! sigma,e_eff,Cy,
     .                      2*NSTR1+1             ! e,strain-part-coeff
C
      NMECH=0
      IF(ITERME.GT.0) THEN                        ! bidirect. coupled
       IF(ITERMP.GT.0) THEN                       ! coupling term
        IF(NITERC.EQ.1.OR.NITERC.EQ.2) NMECH=1    ! improv. stagg. sche.
       ENDIF
      ENDIF
C
C**** PLASTIC MODELS
C
      IPLAS(1)=1                       ! STRAP(NSTR1)
      IPLAS(2)=IPLAS(1)+NSTR1          ! DMTEP(NKOST)
      IPLAS(3)=IPLAS(2)+NKOST          ! EBASE(NBASE)
      IPLAS(4)=IPLAS(3)+NBASE          ! DBASE(NNDAM)
      IPLAS(5)=IPLAS(4)+NNDAM          ! SHRIN
      IPLAS(6)=IPLAS(5)+NSHRI          ! COUTD
      IPLAS(7)=IPLAS(6)+NMECH          ! FATIG(NFATI)
      IPLAS(8)=IPLAS(7)+NFATI          ! VDUAL(NDUAL)
      IPLAS(9)=IPLAS(8)+NDUAL
      NPLAS   =IPLAS(9)-1
      IF(NPLAS.GT.NHIST)
     . CALL RUNEND('ERROR: NPLAS GT NHIST')
C
C**** OUTPUT POINTERS (only for scalar variables & second-order tensors
C                      to be printed)
C
      IPLAO(1)=IPLAS(3)                ! EBASE(1): effect. plast. strain
      IF(NBASE.GE.2) THEN              ! EBASE(>1)
       DO IBASE=2,NBASE
        IPLAO(IBASE)=IPLAO(1)+IBASE-1  ! isot. & kin. hard., etc.
       ENDDO
      ENDIF
      IF(NNDAM.GE.1) THEN              ! DBASE
       DO INDAM=1,NNDAM
        IPLAO(INDAM+NBASE)=IPLAO(NBASE)+INDAM
       ENDDO
      ENDIF
      IF(NSHRI.GE.1) THEN
       DO ISHRI=1,NSHRI                ! SHRIN
        IPLAO(ISHRI+NNDAM+NBASE)=IPLAO(NNDAM+NBASE)+ISHRI
       ENDDO
      ENDIF
      IF(NFATI.GE.1) THEN              ! FATIG
       DO IFATI=1,NFATI
        IPLAO(IFATI+NMECH+NSHRI+NNDAM+NBASE)=
     .        IPLAO(NMECH+NSHRI+NNDAM+NBASE)+IFATI
       ENDDO
      ENDIF
      IF(NDUAL.GE.1) THEN              ! VDUAL
       DO IDUAL=1,2*NSTR1+2+2          ! only sigma, e_eff & Cy
        IPLAO(IDUAL+NFATI+NMECH+NSHRI+NNDAM+NBASE)=
     .        IPLAO(NFATI+NMECH+NSHRI+NNDAM+NBASE)+2*NSTR1+2*NKOST+IDUAL
       ENDDO
      ENDIF
      ICX=2*NSTR1+2+2+NFATI+NMECH+NSHRI+NNDAM+NBASE
      IF(ICX.GT.50)                    ! see auxl_om.f & setdat.f
     . CALL RUNEND('ERROR IN pointe.f; INCREASE DIMENSION OF IPLAO') 
C
      IC=1                             ! effective plastic strain
      NC=1
      IPLAN(IC)=1
      IPLAM(IC)=IPLAN(IC)+NC-1
      IF(KPLA1.EQ.1) THEN              ! plastic isotropic hardening
       IC=IC+1
       NC=1
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
      ENDIF
      IF(KPLA2.EQ.1) THEN              ! kinematic hard. variable
       IC=IC+1
       NC=NSTR1
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
      ENDIF
      IF(KPLA4.EQ.1) THEN              ! porosity + total hard. function
       IC=IC+1
       NC=2
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
      ENDIF
      IF(KPLA7.EQ.1) THEN              ! dot e^p
       IC=IC+1
       NC=1
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
      ENDIF
      IF(KPLA8.EQ.1) THEN              ! deform. resistance & grain size
       IC=IC+1
       NC=2
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
      ENDIF
      IF(KPLA9.EQ.1) THEN              ! old equiv. stress
       IC=IC+1
       NC=1
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
      ENDIF
      IF(KPLA10.EQ.1) THEN             ! iso_inf & viscous stresses
       IC=IC+1
       NC=(1+NCHAI)*NSTR1
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
      ENDIF
      IF(KPLA11.EQ.1) THEN             ! not used now!

      ENDIF
C
      IF(KPLA3.EQ.1) THEN              ! damage or tau- & tau+
       IC=IC+1
       NC=NNDAM
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
      ENDIF
      IF(KPLA5.EQ.1) THEN              ! shrinkage
       IC=IC+1
       NC=NSHRI
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
      ENDIF
      IF(KPLA6.EQ.1) THEN              ! fatigue
       IC=IC+1
       NC=NFATI
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
      ENDIF
      IF(KPLA12.EQ.1) THEN             ! f & m variables
       IC=IC+1                         ! f stress
       NC=NSTR1
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
       IC=IC+1                         ! m stress
       NC=NSTR1
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
       IC=IC+1                         ! f effective plastic strain
       NC=1
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
       IC=IC+1                         ! f hardening function
       NC=1
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
       IC=IC+1                         ! m effective plastic strain
       NC=1
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
       IC=IC+1                         ! m effective plastic strain
       NC=1
       IPLAN(IC)=IPLAM(IC-1)+1
       IPLAM(IC)=IPLAN(IC)+NC-1
      ENDIF
      IF(IC.GT.50)                     ! see auxl_om.f & setdat.f
     . CALL RUNEND('ERROR IN pointe.f; INCREASE DIMENSION OF IPLAN/M') 
C
C**** NNUIN depending on KPLAi
C
      NNUIX=1                          ! check setdat.f
      IF(KPLA1.EQ.1) NNUIX=NNUIX+1                          ! EBASE
      IF(KPLA2.EQ.1) NNUIX=NNUIX+NSTR1
      IF(KPLA4.EQ.1) NNUIX=NNUIX+2
      IF(KPLA7.EQ.1) NNUIX=NNUIX+1
      IF(KPLA8.EQ.1) NNUIX=NNUIX+2
      IF(KPLA9.EQ.1) NNUIX=NNUIX+1
      IF(KPLA10.EQ.1) NNUIX=NNUIX+(1+NCHAI)*NSTR1
C
      IF(KPLA3.EQ.1) NNUIX=NNUIX+NNDAM                      ! DBASE
      IF(KPLA5.EQ.1) NNUIX=NNUIX+NSHRI                      ! SHRIN
      IF(KPLA6.EQ.1) NNUIX=NNUIX+NFATI                      ! FATIG
      IF(KPLA12.EQ.1) NNUIX=NNUIX+2*NSTR1+2+2               ! VDUAL
C
      IF(NNUIX.GT.NNUIN) THEN
C      write(7,*) 'nnuix,nnuin=',nnuix,nnuin
       CALL RUNEND('ERROR: INCREASE NNUIN IN setdat.f')
      ENDIF
      NNUIN=NNUIX
C
C**** NNUNO (maximum index of IPLAN or IPLAM=total number of variables)
C
      NNUNX=IC
      IF(NNUNX.GT.NNUNO) THEN
C      write(7,*) 'nnunx,nnuno=',nnunx,nnuno
       CALL RUNEND('ERROR: INCREASE NNUNO IN setdat.f')
      ENDIF
      NNUNO=NNUNX
C
C**** VISCOPLASTIC MODELS (not used now)
C
      IVLAS(1)=IPLAS(1)
      IVLAS(2)=IPLAS(2)
      IVLAS(3)=IPLAS(3)
      IVLAS(4)=IPLAS(4)
      IVLAS(5)=IPLAS(5)
      IVLAS(6)=IPLAS(6)
      IVLAS(7)=IPLAS(7)
      IVLAS(8)=IPLAS(8)
C
      RETURN
      END
