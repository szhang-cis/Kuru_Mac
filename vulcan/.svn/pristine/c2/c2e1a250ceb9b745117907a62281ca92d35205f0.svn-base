      SUBROUTINE ACTIVE(LACTI,IFFIX,LNODS,DISTO)
C***********************************************************************
C
C**** THIS ROUTINE READS AND SETS DATA FOR ACTIVE ELEMENTS
C
C
C     Notes: 
C
C     This routine must be entered after fixity.f because the
C     IFFIX array can be changed here for non active points belonging
C     to non active elements. Nevertheless, the defined boundary
C     conditions are stored and consequently used for active points.
C     If a non active point becomes active, the defined boundary
C     conditions are restored for it.
C
C     If the defined boundary conditions change, check the active
C     elements option.
C
C     LACTI(IELEM)=0 > non active element (NAE)
C     LACTI(IELEM)=1 > active element (AE)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'inpo_om.f'
      INCLUDE 'prob_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/LDFILE/ITAPE
C
      DIMENSION LACTI(NELEM)
      DIMENSION IFFIX(NTOTV,*), LNODS(NNODE,NELEM), DISTO(*)
C
      IERR1=1
      IF(ITAPE.NE.LUDAT) THEN
       ITAPE=LUACT      ! external .act file
       IF(IOACT.EQ.0) OPEN(UNIT=LUACT,FILE=CJ1,STATUS='OLD',ERR=1000)
       IOACT=1
       IERR1=0
C
 1000  IF(IERR1.NE.0) THEN
        IF(IERR1.EQ.1) WRITE(LURES,901)
        CALL RUNEND('ERROR IN OPENING FILE ')
       ENDIF
      ENDIF      ! itape.ne.ludat
C
      WRITE(LURES,900)
C
C**** RESTORES DEFINED BOUNDARY CONDITIONS
C
      DO ITOTV=1,NTOTV
       IF(KSMUS.EQ.0) THEN
        IFFIX(ITOTV,1)=IFFIX(ITOTV,2)
       ELSE
        IFFIX(ITOTV,1)=IFFIX(ITOTV,3)
        IFFIX(ITOTV,2)=IFFIX(ITOTV,4)
       ENDIF
      ENDDO
C
C**** READS ACTIVE ELEMENTS ARRAY
C
      NPRIN=1
      DO IELEM=1,NELEM
       CALL LISTEN('ACTIVE',NPRIN,ITAPE)
       IELEMX=DINT(PARAM(1))
       IF(IELEM.NE.IELEMX)
     .  CALL RUNEND('ERROR: NON COINCIDENT ELEMENT NUMBERING')
       LACTI(IELEM)=DINT(PARAM(2))
       WRITE(LURES,910) IELEM,LACTI(IELEM)
      ENDDO
C
      WRITE(LURES,925)
C
C**** CHANGES BOUNDARY CONDITIONS ACCORDING TO NON ACTIVE POINTS
C
      DO IELEM=1,NELEM
       IF(LACTI(IELEM).EQ.0) THEN             ! non active element (NAE)
        DO INODL=1,NNODL
         IPOIN=LNODS(INODL,IELEM)
         ISUMMX=0
C
         DO IELEMX=1,NELEM
          DO INODLX=1,NNODL
           IPOINX=LNODS(INODLX,IELEMX)
           IF(IPOIN.EQ.IPOINX) THEN
            IF(LACTI(IELEMX).EQ.1) ISUMMX=ISUMMX+1
           ENDIF             ! ipoin.eq.ipoinx
          ENDDO              ! inodlx=1,nnodl
         ENDDO               ! ielemx=1,nelem
C
         IF(ISUMMX.EQ.0) THEN                 ! NAE surrounded by NAEs
          DO IDOFC=1,NDOFC
           ITOTV=(IPOIN-1)*NDOFC+IDOFC
           IFFIX(ITOTV,1)=1
           IF(KSMUS.NE.0) IFFIX(ITOTV,2)=1
          ENDDO              ! idofc=1,ndofc
         ENDIF               ! isummx.eq.0
C
        ENDDO                ! inodl=1,nnodl
       ENDIF                 ! lacti(ielem).eq.0
      ENDDO                  ! ielem=1,nelem
C
C**** WRITES CHANGED BOUNDARY CONDITIONS
C
      WRITE(LURES,920)
C
      DO IPOIN=1,NPOIN
       ITOTV=(IPOIN-1)*NDOFC
       WRITE(LURES,921) IPOIN,(IFFIX(ITOTV+IDOFC,1),IDOFC=1,NDOFC)
      ENDDO
C
      WRITE(LURES,925)
C
      RETURN
  900 FORMAT(//10X,8H ELEMENT,10X,6H INDEX)
  901 FORMAT(' ERROR IN OPENING ACTIVE ELEMENT OPTION INPUT FILE 49')
  910 FORMAT(10X,I8,10X,I6)
  920 FORMAT(//10X,27HCHANGED BOUNDARY CONDITIONS,
     .        /10X,8H    NODE,10X,6H  CODE)
  921 FORMAT(10X,I8,13X,3I1)
  925 FORMAT(/)
      END
