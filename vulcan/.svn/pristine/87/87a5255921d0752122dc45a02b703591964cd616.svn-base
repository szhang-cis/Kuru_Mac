      SUBROUTINE ADDSOL(NEQNS,NFRON,NLAST,NSTIF,WORK1)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE MEMORY REQUIREMENTS & ADDRESS FOR 
C     SOLVER
C
C***********************************************************************
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_om.f'
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
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1(:)
      REAL*8, ALLOCATABLE :: TEMP(:)               !Variable temporal
C
      NSIZEWORK1=SIZE(WORK1)                       !tamaño actual del arreglo
      ISOLA=1
      IF(NMEMO7M.EQ.1) ISOLA=1+LSTIF
C
C**** SKYLINE SOLVER
C
      IF(KSOLV.EQ.0) THEN
       IUNSY=0
       IF(KSYMM.EQ.0) IUNSY=1
       IFITE=0
       IF(MITCG.GT.0) IFITE=1
       ISOLV( 1)=ISOLA                             ! GSTDI(NEQNS)
       ISOLV( 2)=ISOLV( 1)+NEQNS                   ! GSTLO(NLAST*IUNSY)
       ISOLV( 3)=ISOLV( 2)+NLAST*IUNSY             ! GSTUP(NLAST)
       ISOLV( 4)=ISOLV( 3)+NLAST                   ! CSTIF(NEVAC,NEVAC)
       ISOLV( 5)=ISOLV( 4)+NEVAC*NEVAC             ! ELOAD(NEQNS)
       ISOLV( 6)=ISOLV( 5)+NEQNS                   ! LNUEQ(NTOTV)
       ISOLV( 7)=ISOLV( 6)+(NTOTV*ICHAL+4)/8       ! LPONT(NTOTV)
       ISOLV( 8)=ISOLV( 7)+(NTOTV*ICHAL+4)/8       ! DISIM(NEVAC)
       ISOLV( 9)=ISOLV( 8)+NEVAC                   ! FOREL(NEVAC)
       ISOLV(10)=ISOLV( 9)+NEVAC                   ! ALOAD(NEQNS*IFITE)
       ISOLV(11)=ISOLV(10)+NEQNS*IFITE             ! DELTA(NEQNS*IFITE)
       ISOLV(12)=ISOLV(11)+NEQNS*IFITE             ! DISIC(NEQNS*IFITE)
       ISOLV(13)=ISOLV(12)+NEQNS*IFITE             ! LOCAL(NEVAC*IFITE)
       LSOLV    =ISOLV(13)+(NEVAC*IFITE*ICHAL+4)/8
      ENDIF
C
C**** FRONTAL SOLVER
C
      IF(KSOLV.EQ.1) THEN
       IUNSY=0
       IF(KSYMM.EQ.0) IUNSY=1
       ISOLV(1) =ISOLA                           ! EQRHS(NBUFA)
       ISOLV(2) =ISOLV(1)+NBUFA                  ! EQUAT(NFRON,NBUFA)
       ISOLV(3) =ISOLV(2)+NFRON*NBUFA            ! GLOAD(NFRON)
       ISOLV(4) =ISOLV(3)+NFRON                  ! GSTIF(NSTIF)
       ISOLV(5) =ISOLV(4)+NSTIF                  ! LOCEL(NEVAC)
       ISOLV(6) =ISOLV(5)+(NEVAC*ICHAL+4)/8      ! NACVA(NFRON)
       ISOLV(7) =ISOLV(6)+(NFRON*ICHAL+4)/8      ! NAMEV(NBUFA)
       ISOLV(8) =ISOLV(7)+(NBUFA*ICHAL+4)/8      ! NDEST(NEVAC)
       ISOLV(9) =ISOLV(8)+(NEVAC*ICHAL+4)/8      ! NPIVO(NBUFA)
       ISOLV(10)=ISOLV(9)+(NBUFA*ICHAL+4)/8      ! VECRV(NFRON)
       ISOLV(11)=ISOLV(10)+NFRON                 ! CSTIF(NEVAC,NEVAC)
       ISOLV(12)=ISOLV(11)+NEVAC*NEVAC     ! EQCOL(IUNSY*NFRON*NBUFA)
       LSOLV    =ISOLV(12)+IUNSY*NFRON*NBUFA
      ENDIF
C
C**** PCG SOLVER
C
      IF(KSOLV.EQ.2) THEN
       NSIZE=NDOFC
       IF(NWIDT.LT.0) NSIZE=NDOFC*NDOFC
       ISOLV( 1)=ISOLA                             ! GSTDI(NPOIN*NSIZE)
       ISOLV( 2)=ISOLV( 1)+NPOIN*NSIZE             ! CSTIF(NEVAC,NEVAC)
       ISOLV( 3)=ISOLV( 2)+NEVAC*NEVAC             ! DISIM(NEVAC)
       ISOLV( 4)=ISOLV( 3)+NEVAC                   ! FOREL(NEVAC)
       ISOLV( 5)=ISOLV( 4)+NEVAC                   ! ALOAD(NTOTV)
       ISOLV( 6)=ISOLV( 5)+NTOTV                   ! DELTA(NTOTV)
       LSOLV    =ISOLV( 6)+NTOTV             
      ENDIF
C
C**** GMRES SOLVER
C
      IF(KSOLV.EQ.3) THEN
       IUNSY=0
       IFITE=0
       NLAST=0
       NWORKG=NEQNS*(MKRYL+1)+(MKRYL*(MKRYL+1))/2+4*MKRYL+2
       NCERO=0
C
       ISOLV( 1)=ISOLA                             ! GSTDI(NEQNS)
       ISOLV( 2)=ISOLV( 1)+NEQNS                   ! GSTLO(NLAST*IUNSY)
       ISOLV( 3)=ISOLV( 2)+NLAST*IUNSY             ! GSTUP(NLAST)
       ISOLV( 4)=ISOLV( 3)+NLAST                   ! CSTIF(NEVAC,NEVAC)
       ISOLV( 5)=ISOLV( 4)+NEVAC*NEVAC             ! ELOAD(NEQNS)
       ISOLV( 6)=ISOLV( 5)+NEQNS*NCERO             ! LNUEQ(NTOTV)
       ISOLV( 7)=ISOLV( 6)+(NTOTV*ICHAL+4)/8       ! LPONT(NTOTV)
       ISOLV( 8)=ISOLV( 7)+(NTOTV*ICHAL+4)/8       ! DISIM(NEVAC)
       ISOLV( 9)=ISOLV( 8)+NEVAC*NCERO             ! FOREL(NEVAC)
       ISOLV(10)=ISOLV( 9)+NEVAC*NCERO             ! ALOAD(NEQNS*IFITE)
       ISOLV(11)=ISOLV(10)+NEQNS*IFITE             ! DELTA(NEQNS*IFITE)
       ISOLV(12)=ISOLV(11)+NEQNS*IFITE             ! DISIC(NEQNS*IFITE)
       ISOLV(13)=ISOLV(12)+NEQNS*IFITE             ! LOCAL(NEVAC*IFITE)
       ISOLV(14)=ISOLV(13)+(NEVAC*IFITE*ICHAL+4)/8
       ISOLV(15)=ISOLV(14)                         ! BGMRE(NEQNS)
       ISOLV(16)=ISOLV(15)+NEQNS                   ! XGMRE(NEQNS)
       ISOLV(17)=ISOLV(16)+NEQNS                   ! WGMRE(NEQNS)
       LSOLV    =ISOLV(17)+NWORKG
      ENDIF
C
C**** PARDISO SOLVER
C
      IF(KSOLV.EQ.4) THEN
       IUNSY=0                                       !no cambio aparente
       IF(KSYMM.EQ.0) IUNSY=1
       IFITE=0
       IF(MITCG.GT.0) IFITE=1
       ISOLV( 1)=ISOLA                             ! GSTDI(NEQNS): Lo calcula en la misma rutina
       ISOLV( 2)=ISOLV( 1)+0                       ! GSTLO(NLAST*IUNSY): Lo calcula en la misma rutina
       ISOLV( 3)=ISOLV( 2)+0                       ! GSTUP(NLAST): Lo calcula en la misma rutina
       ISOLV( 4)=ISOLV( 3)+0                       ! CSTIF(NEVAC,NEVAC)
       ISOLV( 5)=ISOLV( 4)+NEVAC*NEVAC             ! ELOAD(NEQNS)
       ISOLV( 6)=ISOLV( 5)+NEQNS                   ! LNUEQ(NTOTV)
       ISOLV( 7)=ISOLV( 6)+(NTOTV*ICHAL+4)/8       ! LPONT(NTOTV)
       ISOLV( 8)=ISOLV( 7)+(NTOTV*ICHAL+4)/8       ! DISIM(NEVAC)
       ISOLV( 9)=ISOLV( 8)+NEVAC                   ! FOREL(NEVAC)
       ISOLV(10)=ISOLV( 9)+NEVAC                   ! ALOAD(NEQNS*IFITE)
       ISOLV(11)=ISOLV(10)+NEQNS*IFITE             ! DELTA(NEQNS*IFITE)
       ISOLV(12)=ISOLV(11)+NEQNS*IFITE             ! DISIC(NEQNS*IFITE)
       ISOLV(13)=ISOLV(12)+NEQNS*IFITE             ! LOCAL(NEVAC*IFITE)
       LSOLV    =ISOLV(13)+(NEVAC*IFITE*ICHAL+4)/8
      ENDIF
C
C**** MEMORY CONTROL
C
c      write(7,*) 'lsolv,mwork1=',lsolv,mwork1
      IF(LSOLV.GT.NSIZEWORK1) THEN                 !revisa el tamaño del arreglo
        WRITE(LURES,*) "IN ADDSOL REALLOCATE WORK1 FROM: ", NSIZEWORK1,
     .             "TO: ", LSOLV
       IERROR=0
       ALLOCATE(TEMP(LSOLV),STAT=IERROR)           !reserva el espacio nuevo
       IF (IERROR.NE.0)                            !Revisa si hay error para reservar
     .  CALL RUNEND("I CAN'T REALLOCATE WORK1 IN ADDSOL")
       TEMP(1:NSIZEWORK1)=WORK1(:)                 !Copia el contenido al nuevo espacio
       CALL MOVE_ALLOC(FROM=TEMP,TO=WORK1)         !Se asigna el nuevo espacio como WORK1
C       CALL RUNEND('ERROR IN DATBAS: MEMORY OUT OF SPACE ') 
      ENDIF
C
C**** DEFINE STARTING POSITION FOR SOLVER MEMORY
C
      LBYTS=LSOLV*8
      IRELE=1
C
      RETURN
      END
