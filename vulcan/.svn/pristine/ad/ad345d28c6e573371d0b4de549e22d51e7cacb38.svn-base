      SUBROUTINE ASSPROT(PROPST,PROPSS,INDEXS)
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS THE MICROSTRUCTURAL PROPERTIES TO THERMAL
C     PROPERTIES
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'                      ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C
C**** MICROSCOPICAL VARIABLES
C
      INCLUDE 'auxl_oms.f'
      INCLUDE 'prob_oms.f'
C
      DIMENSION PROPST(NPROPT,*), PROPSS(NPROPS,*)
C
C**** ASSIGNS MATERIALS PROPERTIES
C
C     NSKIPS=number of properties from thermal to micro
C
      IF(INDEXS.EQ.0) THEN            ! thermal to micro (see inppros.f)
       NSKIPS=1                       ! KPOROT
       DO IMATST=1,NMATST
        DO ISKIPS=1,NSKIPS
         PROPSS(ISKIPS,IMATST)=PROPST(NPROPM+ISKIPS,IMATST)
        ENDDO
       ENDDO
C
       IF(IEVFI.EQ.0) THEN
        NHISTS=NHISTT
        NPOROS=NPOROT
        NNUPCS=NNUPC
        NNUPTS=NNUPT
        NFILLS=NFILL
        ICONVS=ICONVT
        IGALFAS=IGALFA
       ENDIF
      ELSE                            ! micro to thermal
       NPROPX=NPROPT-NPROPM           ! number of microstructural prop.
       DO IMATST=1,NMATST
        DO IPROPT=1,NPROPX
         PROPST(IPROPT+NPROPM,IMATST)=PROPSS(IPROPT,IMATST)
        ENDDO
       ENDDO
C
       IF(IEVFI.EQ.0) THEN
        DO IBASES=1,NBASES
c        IPLAST(IBASES)=IPLASS(IBASES)                   ! not necessary
         IPLAOT(IBASES)=IPLAOS(IBASES)
         IPLANT(IBASES)=IPLANS(IBASES)
         IPLAMT(IBASES)=IPLAMS(IBASES)
        ENDDO
        IF(NNUPM.GT.0) THEN
         DO INUPM=1,NNUPM
          IPLUAT(INUPM)=IPLUAS(INUPM)
         ENDDO
        ENDIF
        IF(NNUPO.GT.0) THEN
         DO INUPO=1,NNUPO
          IPLUOT(INUPO)=IPLUOS(INUPO)
         ENDDO
        ENDIF
        IF(NNUINS.GT.NNUINT)          ! check setdatt.f
     .   CALL RUNENDT('ERROR IN ASSPROT: NNUINS GT NNUINT')
        NNUINT=NNUINS
        IF(NNUNOS.GT.NNUNOT)          ! check setdatt.f
     .   CALL RUNENDT('ERROR IN ASSPROT: NNUNOS GT NNUNOT')
        NNUNOT =NNUNOS
        NNUPC  =NNUPCS
        NNUPT  =NNUPTS
        NNUM4T =NNUM4S
        INNUM4T=INNUM4S
        IMNMT  =IMNMS
        IWEUT  =IWEUS
       ENDIF
      ENDIF                           ! indexs.eq.0
C
      RETURN
      END
