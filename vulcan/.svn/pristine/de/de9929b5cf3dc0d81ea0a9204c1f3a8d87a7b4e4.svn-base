      SUBROUTINE PROINP
C***********************************************************************
C
C**** THIS ROUTINE REDEFINES THE DEFAULT PROBLEM PARAMETERS 
C
C     PARAMETERS         CONTROLLING WORDS
C
C     KDYNA              'DYNAM'
C     KPORE              'COUPL'
C     KPROB              'SEEPA', 'STRUC', 'WAVES', 'THERM', 'INCOM'
C     LARGE              'LARGE'
C     others             'DIMEN'
C     NPOIC,NOCOI        'CONTA'
C     NSKEW              'LOCAL'
C     COMMUNICATION      'COMMU'
C     NACTI              'ACTIV'
C     NLDSF              'DEFOR'
C     NANIS              'ANISO'
C     IGALE              'WEAK_'
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      PARAMETER (MCOMM=17)
      CHARACTER*5 COMMD(MCOMM)
C
      DATA COMMD/'SEEPA','STRUC','WAVES','THERM','INCOM','DYNAM',
     .           'COUPL','LARGE','DIMEN','CONTA','LOCAL','COMMU',
     .           'ACTIV','DEFOR','ANISO','WEAK_','END_P'/
C
C**** READ PROBLEM PARAMETERS
C
      NPRIN=0
      ITAPE=LUDAT
      CALL LISTEN('PROINP',NPRIN,ITAPE)
      IF(WORDS(1).NE.'PROBL') GO TO 2000
 1000 CALL LISTEN('PROINP',NPRIN,ITAPE)
C
C**** IDENTIFY COMMAND
C
      DO ICOMM=1,MCOMM
       IF(WORDS(1).EQ.COMMD(ICOMM)) GO TO 100
      END DO
      GO TO 2000 
C
C**** EXECUTE APROPRIATE COMMAND
C
  100 CONTINUE
      GO TO (1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13), ICOMM
C
    1 IF(WORDS(1).EQ.'SEEPA')                         KPROB=0
      IF(WORDS(1).EQ.'STRUC')                         KPROB=1
      IF(WORDS(1).EQ.'STRUC'.AND.WORDS(2).EQ.'SHELL') KPROB=2
      IF(WORDS(1).EQ.'WAVES')                         KPROB=3
      IF(WORDS(1).EQ.'THERM')                         KPROB=4
      IF(WORDS(1).EQ.'INCOM')                         KPROB=5
      GO TO 1000
C
    2 KDYNA=1
      NDISR=3
      NDISO=5    ! maximum of the available KINTE options (see strinp.f)
      GO TO 1000
C
    3 KPORE=2
      IF(WORDS(2).EQ.'DISCR') KSMUS=-1
      IF(WORDS(2).EQ.'LOCAL') KSMUS=-2
      GO TO 1000
C
C**** LARGE DISPLACEMENTS/STRAINS ANALYSIS
C
C     Notes:
C
C     LARGE=1: Total Lagrangian Formulation
C     LARGE=2: Updated Lagrangian Formulation
C     LARGE=3: Eulerian Formulation
C
    4 LARGE=1                                   ! default
      IF(WORDS(2).NE.'     ') THEN
       IF(WORDS(2).NE.'TOTAL'.AND.WORDS(2).NE.'UPDAT'.AND.
     .    WORDS(2).NE.'EULER')
     .  CALL RUNEND('PROINP: WRONG LARGE STRAINS INPUT DATA')
       IF(WORDS(2).EQ.'TOTAL') LARGE=1
       IF(WORDS(2).EQ.'UPDAT') LARGE=2
c      IF(WORDS(2).EQ.'EULER') LARGE=3
       IF(WORDS(2).EQ.'EULER')
     .  CALL RUNEND('ERROR: EULERIAN FORMULATION NOT IMPLEMENTED')
      ENDIF
      GO TO 1000
C
    5 DWORD(1)='NPOIN'
      DWORD(2)='NELEM'
      DWORD(3)='NDIME'
      DWORD(4)='NNODE'
      DWORD(5)='NGAUS'
      DWORD(6)='NSETS'
      DWORD(7)='NMATS'
      DWORD(8)='NFUNC'
      DPARA(1)=0.
      DPARA(2)=0.
      DPARA(3)=2.
      DPARA(4)=4.
      DPARA(5)=4.
      DPARA(6)=1.
      DPARA(7)=1.
      DPARA(8)=1.
C
      CALL SORTER(1)
C
      NPOIN=INT(DPARA(1))
      NELEM=INT(DPARA(2))
      NDIME=INT(DPARA(3))
      NNODE=INT(DPARA(4))
      NGAUS=INT(DPARA(5))
      NGRUP=INT(DPARA(6))
      NMATS=INT(DPARA(7))
      NFUNC=INT(DPARA(8))
C
      IF(NGRUP.LT.NMATS)
     . CALL RUNEND('NGRUP CAN NOT BE LESS THAN NMATS')
      GO TO 1000
C
    6 CONTINUE
C
C**** MECHANICAL CONTACT
C
C     Notes:
C     For non-coincident mesh (NOCOI>0) and large displacements
C     (LARGC=1), the assumptions are (October-2003):
C     - The contact geometry (specifically, the three tasks according to
C     NDIME defined in setmtx.f for NOCOI>0) is incrementally (not
C     iteratively!) updated in the current version (i.e., LICOI must be
C     zero; see below). This is the reason why the routine setmtx.f is
C     only called by tmstep.f and not by iterat.f to compute the contact
C     forces.
C     - The geometric contact stiffness matrix is not considered (i.e.,
C     the derivative of the normal vector with respect to the displac.
C     vector is neglected).
C     - It is possible to use concident & non-coincident contact meshes
C     simultaneously (NOCOI=2). In this case, only the new contact
C     definition (i.e., with ISETC & IMAES defined at the set level)
C     must be used. Note that in this case the LARGE_DISPLACEMENTS is
C     only valid for the non-coincident contacts.
C
      IIX=2
      IIY=0
      IF(WORDS(IIX).EQ.'AUGME') THEN
       IIX=IIX+1
       IF(WORDS(IIX).EQ.'NPOIC') THEN
        IIX=IIX+1
        IIY=IIY+1
        NPOIC=INT(PARAM(IIY))
        IF(NPOIC.LT.0) THEN
         CALL RUNMEN('WARNING: CONTACT PROBLEM USING A.L. APPROACH WITH
     . NEGATIVE NUMBER OF CONTACT NODES')
         NPOIC=0
        ENDIF
        IF(WORDS(IIX).EQ.'VERSI') THEN
         IIX=IIX+1
         IIY=IIY+1
         IAUGM=INT(PARAM(IIY))
         IF(IAUGM.LT.0.OR.IAUGM.GT.5) THEN
          CALL RUNMEN('WARNING: WRONG NUMBER OF AUGM. LAGRAG. VERSION 
     .                                  ')
          IAUGM=0
         ENDIF
         IF(WORDS(IIX).EQ.'NNODC') THEN
          IIX=IIX+1
          IIY=IIY+1
          NNODC=INT(PARAM(IIY))
          IF(NNODC.LT.0) THEN
           CALL RUNMEN('WARNING: WRONG NUMBER OF AUGM. LAGRAG. NNODC
     .                                   ')
           NNODC=0
          ENDIF
         ENDIF
        ELSE
         IF(WORDS(IIX).EQ.'NNODC') THEN
          IIX=IIX+1
          IIY=IIY+1
          NNODC=INT(PARAM(IIY))
          IF(NNODC.LT.0) THEN
           CALL RUNMEN('WARNING: WRONG NUMBER OF AUGM. LAGRAG. NNODC
     .                                   ')
           NNODC=0
          ENDIF
         ENDIF
        ENDIF
       ELSE
        CALL RUNMEN('WARNING: CONTACT PROBLEM USING A.L. APPROACH WITH
     . ZERO NUMBER OF CONTACT NODES')
       ENDIF
      ENDIF
      IF(WORDS(IIX).EQ.'NON_C'.OR.WORDS(IIX).EQ.'BOTH_') THEN
       IF(WORDS(IIX).EQ.'NON_C') NOCOI=1 ! non-coincident mesh
       IF(WORDS(IIX).EQ.'BOTH_') NOCOI=2 ! both coinc. & non-coinc. mesh
       IIX=IIX+1
       IF(WORDS(IIX).EQ.'LARGE') THEN    ! large displac. for contact
        LARGC=1
        IIX=IIX+1
        IF(WORDS(IIX).EQ.'NON_L') THEN   ! non-linearized comput. of n
         LICOI=1
         IIX=IIX+1
c        call runend('ERROR: NON-LINEARIZED COMPUTATION OF N NOT IMPL.')
        ENDIF
        IF(WORDS(IIX).EQ.'MAXIM') THEN   ! maxim. gap possible contact
         IIX=IIX+1
         IIY=IIY+1
         GAPPC=PARAM(IIY)
        ENDIF
        IF(WORDS(IIX).EQ.'SKIPP') THEN   ! skipping steps for n comput.
         IIX=IIX+1
         IIY=IIY+1
         NSKIC=INT(PARAM(IIY))
        ENDIF
       ENDIF
      ENDIF
      GO TO 1000
C
    7 CONTINUE
C
C**** LOCAL SYSTEM OF COORDINATES (only useful for prescribed displac.)
C
      NSKEW=INT(PARAM(1))
      IF(NSKEW.LT.0)
     . CALL RUNEND('ERROR: WRONG NUMBER OF LOCAL SYSTEM OF COORD.')
      GO TO 1000
C
    8 CONTINUE
      GO TO 1000
C
    9 NACTI=1                         ! active elements
      GO TO 1000
C
   10 CONTINUE
C
C**** DEFORMATION-DEPENDENT FACE LOAD
C
C     Note:
C     For deformation-dependent loads (NLDSF=1), the assumptions are
C     (January-2004):
C     - only face loads can be deformation-dependent (i.e., points or
C     gravity loads cannot be deformation-dependent).
C     - The load surface is incrementally (not iteratively!) updated in 
C     the current version (i.e., LLDSF must be zero; see below).
C     - The deformation-dependent face load stiffness matrix is not
C     considered (i.e., the derivative of the def.-dep. face load with
C     respect to the displacement vector is neglected).
C
      NLDSF=1
      IF(WORDS(2).EQ.'NON_L') THEN     ! non-linearized computation
       LLDSF=1
       call runend('ERROR: NON-LIN. COMP. DEF.DEP.FACE LOAD NOT IMPL.')
      ENDIF
      GO TO 1000
C
   11 CONTINUE
C
C**** ANISOTROPIC CONSTITUTIVE MODELS
C
      NANIS=1
      NANIV=2                          ! default (anisotropic vectors)
      NANIC=3                          ! default (vector dimension)
      IF(WORDS(2).EQ.'FIBER') THEN
       NANIV=INT(PARAM(1))             ! number of fiber vectors (max.)
       NANIC=4                         ! dimension + volumetric fraction
      ENDIF
      GO TO 1000
C
   12 CONTINUE
C
C**** WEAK FORM
C
      IF(WORDS(2).EQ.'DISCO') THEN     ! Discontinuous Galerkin formul.
       IF(NMEMO7M.EQ.1)
     .  CALL RUNEND('ERROR: IGALE=1 & NMEMO7M=1 NOT IMPLEMENTED YET')
       IGALE=1
       IF(NDIME.EQ.0)
     .  CALL RUNEND('ERROR: PUT WEAK_FORM CARD AFTER DIMENSIONS CARD')
       IF(NDIME.EQ.1) NPRE4=1          ! NSTR1 defined in conset.f
       IF(NDIME.EQ.2) NPRE4=4
       IF(NDIME.EQ.3) NPRE4=6
       NPRE5=(NPRE4+1)*NPRE4/2         ! NKOST defined in addelm.f
       IF(KSYMM.EQ.0) NPRE5=NPRE4*NPRE4
       NPREA=NPRE1+NPRE2+NPRE3+NPRE4+NPRE5
      ENDIF
      GO TO 1000
C
   13 CONTINUE
      NPOIN=NPOIN+NPOIC
      NNODE=NNODE+NNODC
C
C**** ADDITIONAL CONTROLS
C
      IF(NPOIC.GT.NPOIN)
     . CALL RUNEND('ERROR: WRONG NUMBER OF CONTACT POINTS        ')
C
      IF(NSKEW.GT.NPOIN)
     . CALL RUNEND('ERROR: WRONG NUMBER OF LOCAL SYSTEM OF COORD.')
C
      IF(KDYNA.EQ.1.AND.NACTI.EQ.1)
     . CALL RUNEND('ERROR: KDYNA=NACTI=1 NOT IMPLEMENTED YET')
C
      IF(NOCOI.GT.0) THEN
       IF(LARGC.NE.0) THEN
        IF(NMEMO2M.EQ.0)
     .   CALL RUNEND('ERROR: NOCOI>0, LARGC=1 & NMEMO2M=0 ARE INCOMP.')
       ENDIF
      ENDIF
C
      IF(NLDSF.EQ.1.AND.LARGE.EQ.0)
     . CALL RUNEND('ERROR: LARGE MUST BE 1 FOR DEF.-DEP. LOADS')
C
      IF(NANIS.EQ.1.AND.NDIME.EQ.1)
     . CALL RUNEND('ERROR: 1D PROBLEMS ARE NOT POSSIBLE FOR ANISOTROPY')
C
      RETURN
 2000 CALL RUNEND('PROINP: ERROR IN PROBLEM DATA BLOCK')
      END
