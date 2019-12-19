      SUBROUTINE ACTIVET(LACTIT,IFFIXT,LNODST)
C***********************************************************************
C
C**** THIS ROUTINE READS AND SETS DATA FOR ACTIVE ELEMENTS
C
C
C     Notes: 
C
C     This routine must be entered after fixityt.f because the
C     IFFIXT array can be changed here for non active points belonging
C     to non active elements. Nevertheless, the defined boundary
C     conditions are stored and consequently used for active points.
C     If a non active point becomes active, the defined boundary
C     conditions are restored for it.
C
C     If the defined boundary conditions change, check the active
C     elements option.
C
C     LACTIT(IELEMT)=0 > non active element (NAE)
C     LACTIT(IELEMT)=1 > active element (AE)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'inpo_omt.f'
      INCLUDE 'prob_omt.f'
      INCLUDE 'auxl_omt.f'
C
      COMMON/LDFILET/ITAPET
C
      DIMENSION LACTIT(NELEMT)
      DIMENSION IFFIXT(NTOTVT,*), LNODST(NNODET,NELEMT)
C
      IERR1T=1
      IF(ITAPET.NE.LUDATT.AND.ITAPET.NE.LUACTT) GOTO 1000
      IF(ITAPET.EQ.LUACTT) THEN       ! external .act file
       IF(IOACTT.EQ.0) OPEN(UNIT=LUACTT,FILE=CH1T,STATUS='OLD',ERR=1000)
       IOACTT=1
       IERR1T=0
C
 1000  IF(IERR1T.NE.0) THEN
        IF(IERR1T.EQ.1) WRITE(LUREST,901)
        CALL RUNENDT('ERROR IN OPENING FILE ')
       ENDIF
      ENDIF      ! itapet.eq.luactt
C
      WRITE(LUREST,900)
C
C**** RESTORES DEFINED BOUNDARY CONDITIONS
C
      DO ITOTVT=1,NTOTVT
       IF(KSMUST.EQ.0) THEN
        IFFIXT(ITOTVT,1)=IFFIXT(ITOTVT,2)
       ELSE
        IFFIXT(ITOTVT,1)=IFFIXT(ITOTVT,3)
        IFFIXT(ITOTVT,2)=IFFIXT(ITOTVT,4)
       ENDIF
      ENDDO
C
C**** READS ACTIVE ELEMENTS ARRAY
C
      NPRINT=1
      DO IELEMT=1,NELEMT
       CALL LISTENT('ACTIVET',NPRINT,ITAPET)
       IELEMX=DINT(PARAMT(1))
       IF(IELEMT.NE.IELEMX)
     .  CALL RUNENDT('ERROR: NON COINCIDENT ELEMENT NUMBERING')
       LACTIT(IELEMT)=DINT(PARAMT(2))
       WRITE(LUREST,910) IELEMT,LACTIT(IELEMT)
      ENDDO
C
      WRITE(LUREST,925)
C
C**** CHANGES BOUNDARY CONDITIONS ACCORDING TO NON ACTIVE POINTS
C
      DO IELEMT=1,NELEMT
       IF(LACTIT(IELEMT).EQ.0) THEN           ! non active element (NAE)
        DO INODLT=1,NNODLT
         IPOINT=LNODST(INODLT,IELEMT)
         ISUMMX=0
C
         DO IELEMX=1,NELEMT
          DO INODLX=1,NNODLT
           IPOINX=LNODST(INODLX,IELEMX)
           IF(IPOINT.EQ.IPOINX) THEN
            IF(LACTIT(IELEMX).EQ.1) ISUMMX=ISUMMX+1
           ENDIF             ! ipoint.eq.ipoinx
          ENDDO              ! inodlx=1,nnodlt
         ENDDO               ! ielemx=1,nelemt
C
         IF(ISUMMX.EQ.0) THEN                 ! NAE surrounded by NAEs
          DO IDOFCT=1,NDOFCT
           ITOTVT=(IPOINT-1)*NDOFCT+IDOFCT
           IFFIXT(ITOTVT,1)=1
           IF(KSMUST.NE.0) IFFIXT(ITOTVT,2)=1
          ENDDO              ! idofct=1,ndofct
         ENDIF               ! isummx.eq.0
C
        ENDDO                ! inodlt=1,nnodlt
       ENDIF                 ! lactit(ielemt).eq.0
      ENDDO                  ! ielemt=1,nelemt
C
C**** WRITES CHANGED BOUNDARY CONDITIONS
C
      WRITE(LUREST,920)
C
      DO IPOINT=1,NPOINT
       ITOTVT=(IPOINT-1)*NDOFCT
       WRITE(LUREST,921) IPOINT,
     .                   (IFFIXT(ITOTVT+IDOFCT,1),IDOFCT=1,NDOFCT)
      ENDDO
C
      WRITE(LUREST,925)
C
      RETURN
  900 FORMAT(//10X,8H ELEMENT,10X,6H INDEX)
  901 FORMAT(' ERROR IN OPENING ACTIVE ELEMENT OPTION INPUT FILE 247',
     . '(59 LINUX)')
  910 FORMAT(10X,I8,10X,I6)
  920 FORMAT(//10X,27HCHANGED BOUNDARY CONDITIONS,
     .        /10X,8H    NODE,10X,6H  CODE)
  921 FORMAT(10X,I8,13X,3I1)
  925 FORMAT(/)
      END
