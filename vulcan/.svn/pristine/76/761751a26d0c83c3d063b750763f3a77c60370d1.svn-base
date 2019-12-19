      SUBROUTINE POINTES
C***********************************************************************
C
C**** THIS ROUTINE DETERMINES POINTERS OF EHIST ARRAY FOR THE INTERNAL
C     VARIABLES FOR THE MICROSTRUCTURAL MODELS FOR THERMAL-MICROSTRUC-
C     TURAL PROBLEMS WITH MICRO EVOLUTION EQUATIONS (IEVFI=0)
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
      INCLUDE 'nuef_om.f'
C
C**** MICROSCOPICAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'auxl_oms.f'
C
C**** DETERMINES (ACCORDING TO KPLAiS):
C
C     1) NUMBER OF COMPONENTS OF THE INTERNAL VARIABLES (NBASES)
C     2) POINTERS OF INTERNAL VARIABLES (IPLASS)
C     3) POINTERS OF OUTPUT VECTORS OF INTERNAL VARIABLES
C        (IPLAOS(I), I=1,NNUINS); see smoi05t.f
C     4) POINTERS OF OUTPUT NODAL VECTORS OF INTERNAL VARIABLES
C        (from IPLANS(I) to IPLAMS(I), I=1,NNUNOS); see outnodt.f
C     5) NUMBER OF COMPONENTS OF THE NUMBER OF INTERNAL VARIABLES TO
C        BE PRINTED IN THE POST FILE (NNUINS); see outnodt.f
C     6) NUMBER OF INTERNAL VARIABLES TO BE SMOOTHED AND PRINTED IN
C        RESULTS FILE (NNUNOS); see smoi05t.f & outnodt.f
C
C     NOTES:
C
C     a) THE COMPONENTS OF IPLAOS, IPLANT AND IPLAMT MUST BE EQUAL
C     b) NNUINS MUST BE EQUAL TO \sum_{i=1,nnunos} IPLANS(i) (IPLAMS(i))
C     c) FOR SCALARS INTERNAL VARIABLES TO PRINT, IPLANS(I)=1 AND
C        IPLAMS(I)=1 FOR ALL I (NNUINS=NNUNOS)
C        (NNUINS=NNUNOS)
C     d) PARAMETERS KPLAiS DEFINES THE TYPE OF CONSTITUTIVE MODEL => 
C        THE FOLLOWING VARIABLES DEPEND ON KPLAiS
C
C     NHISTS: this variable is equal to NHISTT (see assprot.f).
C     NNUM4S: this variable is assigned to NNUM4T (see assprot.f).
C     NNUINS: this variable is assigned to NNUINT (see assprot.f).
C     NNUNOS: this variable is assigned to NNUNOT (see assprot.f).
C
C     Check (change if necessary) dimensions of:
C     NNUM4TS & NNUM4TMS in setdatd.f
C     MNU4SX & MNU4SMX in auxl_oms.f (ge than NNUM4TS & NNUM4TMS)
C     MNU4TX & MNU4TMX in auxl_omt.f (ge than NNUM4TS & NNUM4TMS)
C     MNUM4S & MNUM4SM in gend_oms.f (ge than NNUM4TS & NNUM4TMS)
C     MNUM4T & MNUM4TM in gend_omt.f (ge than NNUM4TS & NNUM4TMS)
C
C     The following indexes are automatically changed:
C     NHISTT          in consett.f
C     NNUINT & NNUNOT in setdatt.f
C     NNUINS & NNUNOS in setdats.f
C
C     The variable NNUM4TMS is up-to-date estimated in setdatd.f with
C     model KPLA9S=1 since it is the model that needs more variables
C     (see below; Jl/04)
C
C     List of subroutines to be changed (summary):
C     setdatd.f
C     auxl_oms.f, auxl_omt.f 
C     gend_oms.f, gend_omt.f
C
      NBASES=0
      IF(KPLA1S.EQ.1)  NBASES=NBASES+1+1   ! fpc1 & fpc2
      IF(KPLA2S.EQ.1)  NBASES=NBASES+14    ! fl,fa,fg,Ng,Rg,fc,Nc,Rc,
C                                          ! fgam,fgr,lg,Sg,ig,ic
      IF(KPLA3S.EQ.1)  NBASES=NBASES+14    ! fl,fd,fdg,fdi,Nd,Rd,
C                                          ! fe,feg,fie,fei,Ne,Re,id,ie
      IF(KPLA4S.EQ.1)  NBASES=NBASES+4+    ! fl,fa,fg,Si,
     .                        IMNMS+       ! Mn
     .                 (2+INNUM4S)*NNUM4S+ ! nnum4s*(Ng,Rg,Rgn),
     .                        3*IWEUS+2    ! fc,Nc,Rc,jo,ir
      IF(KPLA5S.EQ.1)  NBASES=NBASES+4+    ! fl,fa,fg,Si,
     .                        3*NNUM4S+2   ! nnum4s*(Ng,Rg,Ra),jo,ir
      IF(KPLA6S.EQ.1)  NBASES=NBASES+6     ! max. variables of AJM eq.
C                                          ! (iversi=1,3)
      IF(KPLA7S.EQ.1)  NBASES=NBASES+8     ! max. variables of hydration
C                                          ! (iversi=2)
      IF(KPLA8S.EQ.1)  NBASES=NBASES+4     ! AlN precip. & recrystal.
      IF(KPLA9S.EQ.1)  NBASES=NBASES+26+   ! fl,fa,fg,Si,...
     .                        5*NNUM4S+5   ! nnum4s*(Ng*3,Rg*2),jo,ir,
C                                          ! ix,irg,tf
      IF(KPLA10S.EQ.1) NBASES=NBASES+12    ! fl,fs,fe,Nc,Rc,ic,fid,fie,
C                                          ! fsopc,TSOCPC,TSOCCP,TSOCCE
      IF(KPLA11S.EQ.1) NBASES=NBASES+7     ! ff,Nf,Rf,fp,Np,Rp,Ts
      IF(KPLA12S.EQ.1) NBASES=NBASES+4+    ! ff,fp,cag,cprop,
     .                        6*NNUM4S+11  ! nnum4s*(lg,lf,an,Rf,Np,Rp),
C                                          ! cao,fsolaini,rfa,dperi,jp,
C                                          ! irec,tmin,tcaotic,caotoid,
C                                          ! ...
      IF(KPLA14S.EQ.1)  NBASES=NBASES+14   ! max. variables of AJM eq.
C                                          ! (iversi=1,3)
      IF(NBASES.GT.NNUM4SM)
     . CALL RUNENDS('POINTES: ERROR => NBASES GT NNUM4SM  ')
C
C**** MICROSTRUCTURAL MODELS (only for scalar variables)
C
      IPLASS(1)=NHISTS-NHISTM+1            ! Int. var. 1
      DO I=2,NBASES                        ! Int. var. 2-nbases
       IPLASS(I)=IPLASS(I-1)+1
      ENDDO
      NPLASS    =IPLASS(NBASES)
      IF(NPLASS.GT.NHISTS)
     . CALL RUNENDS('POINTES: ERROR => NPLASS GT NHISTS   ')
C
C**** OUTPUT POINTERS (only for scalar variables)
C
      DO I=1,NBASES
       IPLAOS(I)=IPLASS(I)
       IPLANS(I)=1
       IPLAMS(I)=1
      ENDDO
C
C**** NNUINS depending on KPLAiS
C     (FOR RATE PHASE-CHANGE FORMULATIONS (IPCFO=1))
C
C     For NPLAT phase-changes, the possible combinations are:
C
C     1) n (n ge 1) KPLA1S=1
C     2) n_1 (n_1 ge 1) KPLA1S=1 with n_2 (n_2=1) KPLA2S=1
C     3) n_1 (n_1 ge 1) KPLA1S=1 with n_2 (n_2=1) KPLA2S=1
C     4) n_1 (n_1=1) KPLA1S=1 with n_2 (n_2=1) KPLA2S=1
C
C     Several combinations are lacking !!
C
C     For one phase-change (NPLAT=1) for each material, the first term
C     is only valid in the computation of NNUINS.
C     For multiple phase-changes (NPLAT > 1) for each material,
C     additional terms have to be considered in the computation of
C     NNUINS.
C     i.e., for micro phase-changes corresponding to MODEL=1,
C           for NPLAT=1, => IF(KPLA1S.EQ.1)  NNUINS=NNUINS+1
C           for NPLAT=2, => IF(KPLA1S.EQ.1)  NNUINS=NNUINS+1+1
C
      NNUINS=0
      IF(KPLA1S.EQ.1)  NNUINS=NNUINS+1+1   ! check input data
      IF(KPLA2S.EQ.1)  NNUINS=NNUINS+12    ! only 12 to be printed
      IF(KPLA3S.EQ.1)  NNUINS=NNUINS+12    ! only 12 to be printed
      IF(KPLA4S.EQ.1)  NNUINS=NNUINS+4+    ! only those to be printed
     .                        IMNMS+(2+INNUM4S)*NNUM4S+3*IWEUS
      IF(KPLA5S.EQ.1)  NNUINS=NNUINS+4+    ! only those to be printed
     .                        3*NNUM4S
      IF(KPLA6S.EQ.1)  NNUINS=NNUINS+6     ! only those to be printed
      IF(KPLA7S.EQ.1)  NNUINS=NNUINS+8     ! max. variables of hydration
      IF(KPLA8S.EQ.1)  NNUINS=NNUINS+4     ! AlN precip. & recrystal.
      IF(KPLA9S.EQ.1)  NNUINS=NNUINS+26+   ! only those to be printed
     .                        5*NNUM4S
      IF(KPLA10S.EQ.1) NNUINS=NNUINS+12    ! only those to be printed
      IF(KPLA11S.EQ.1) NNUINS=NNUINS+7     ! ff,Nf,Rf,fp,Np,Rp,Ts
      IF(KPLA12S.EQ.1) NNUINS=NNUINS+4+    ! only those to be printed
     .                        6*NNUM4S
      IF(KPLA14S.EQ.1)  NNUINS=NNUINS+14   ! only those to be printed
C
C**** FOR SMOOTHING OPERATIONS TO PASS SMSTPT TO FPCHAT (see outsmot.f)
C
C     Notes:
C
C     IPLLLS is defined in propmic"i".f (i=1,10) according to INUPM and
C     gives the value "i".
C
C     INDEX1 is the number of microstructural variables to be smoothed
C     in the postprocess
C
C     INDEX2 is similar to KPLAiS but takes into account the order of
C     NNUPM (order given in the input data file).
C
C     INDEX3 is the number of variables to be transferred to the
C     mechanical problem or, for advective problems, it denotes the
C     number of variables defined in terms of evolution equations (by
C     means of a material derivative; see microst.f).
C     If INDEX3 > 0, add INDEX3 positions of array IPLXXS. This array
C     contains the number of the variable to be passed to FPCHAT in
C     outsmot.f (note that IPLXXS(I) should be > 1).
C
C     NNUPO: other microstructural variables to transfer to mechanical
C     problem by means of FPCHAT (excluding the phase-change functions).
C
      INDEX0=0
      NNUPO=0
      DO INUPM=1,NNUPM
       INDEX2=IPLLLS(INUPM)
       IF(INDEX2.EQ.1) THEN
        INDEX1=1               ! model variables to be smoothed
        INDEX3=0               ! other variables to fpchat
       ENDIF
       IF(INDEX2.EQ.2) THEN
        INDEX1=12              ! model variables to be smoothed
        INDEX3=0               ! be careful
c       IF(ICONVS.EQ.1) THEN
c        INDEX3=2
c        IPLXXS(1)=5           ! fifth variable: graphite radius
c        IPLXXS(2)=8           ! eighth variable: cementite radius
c       ENDIF
       ENDIF
       IF(INDEX2.EQ.3) THEN
        INDEX1=0
        INDEX3=0
       ENDIF
       IF(INDEX2.EQ.4) THEN
        INDEX1=4+              ! model variables to be smoothed
     .         IMNMS+(2+INNUM4S)*NNUM4S+3*IWEUS
        INDEX3=0               ! no additional variables to mech. prob.
       ENDIF
       IF(INDEX2.EQ.5) THEN
        INDEX1=4+3*NNUM4S      ! model variables to be smoothed
        INDEX3=0
       ENDIF
       IF(INDEX2.EQ.6) THEN
        INDEX1=6
        INDEX3=0               ! no additional variables to mech. prob.
        IPLXXS(1)=2            ! vfA
        IPLXXS(2)=3            ! cA
        IPLXXS(3)=4            ! cAR
       ENDIF
       IF(INDEX2.EQ.7) THEN
        INDEX1=7
        INDEX3=3               ! 3 var. to be transferred to mech. prob.
        IPLXXS(1)=2            ! second variable (compression strength)
        IPLXXS(2)=3            ! third variable  (tension strength)
        IPLXXS(3)=4            ! fourth variable (Young modulus)
       ENDIF
       IF(INDEX2.EQ.8) THEN
        INDEX1=4
        INDEX3=0
       ENDIF
       IF(INDEX2.EQ.9) THEN
        INDEX1=26+5*NNUM4S     ! model variables to be smoothed
        INDEX3=0
       ENDIF
       IF(INDEX2.EQ.10) THEN
        INDEX1=12              ! model variables to be smoothed
        INDEX3=0
       ENDIF
       IF(INDEX2.EQ.11) THEN   !
        INDEX1=7               ! model variables to be smoothed
        INDEX3=0
       ENDIF
       IF(INDEX2.EQ.12) THEN
        INDEX1=4+6*NNUM4S      ! model variables to be smoothed
        INDEX3=0
       ENDIF
       IF(INDEX2.EQ.14) THEN
        INDEX1=14              ! model variables to be smoothed
        INDEX3=0
       ENDIF
C
       INDEX0=INDEX0+1         ! phase-change function
       IPLUAS(INUPM)=INDEX0
C
       IF(INDEX3.GT.0) THEN
        DO INUPO=1,INDEX3
         IPLUOS(NNUPO+INUPO)=INDEX0-1+IPLXXS(INUPO)
        ENDDO
       ENDIF
C
       NNUPO=NNUPO+INDEX3
       INDEX0=INDEX0+INDEX1-1
      ENDDO
C
C**** NNUNOS (maximum index of IPLAN or IPLAM=total number of variables)
C
C     Note: taking into account that all the microstructural variables
C           of the models currently implemented are scalars  (Jan/97)
C           => NNUNOS=NNUINS
C
      NNUNOS=NNUINS
C
C**** CHECK
C
C     Note: change NFPCH in setdatt.f
C
      NLINEA=NFPCH
      NLINET=2*NNUPTS+NNUPO+NFILLS+IGALFAS
      IF(NLINET.GT.NLINEA)
     . CALL RUNENDS('POINTES: WRONG NUMBER OF PHASE_CHANGES')
C
      RETURN
      END
