      SUBROUTINE DATBAS(ARRAY,IFLAG,ITASK)
C***********************************************************************
C
C**** THIS ROUTINE INTERFACES VULCAN-DATA WITH DATABASE-FILE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'prob_om.f'
C
C**** SE DECLARAN BUFFER Y TEMP COMO VARIABLES DINAMICAS, INIT saber si se inicializó,
C
      INTEGER, SAVE :: INIT=0
      REAL*8, ALLOCATABLE :: BUFFER(:),TEMP(:) !TEMP es la variable temporal para redimensionar la variable buffer
C      COMMON/PRUEBA/BUFFER
      SAVE NBUFFER,BUFFER
      DIMENSION ARRAY(*)
C
C**** BEGIN
C
      NUNIT=LUDTS                                  ! .dts
C
      IF(IDATP(IFLAG,5).LT.0) THEN                 ! DATA IS OUT-CORE
C
       NRECD=(IELEM-1)*NLENP+IDATP(IFLAG,4)+1
       CALL RESTDK(NUNIT,ITASK,ARRAY,IDATP(IFLAG,1)*2,LENRC,NRECD)
C
      ELSE                                         ! DATA IS IN-CORE
C
       IF(IFLAG.EQ.12) THEN ! ACCESS TO DISTO
        NRECD=1
       ELSE 
        NRECD=NTOTV+(IELEM-1)*NWORP+IDATP(IFLAG,3)+1
       ENDIF
C*****************************************************************       
C
C**** INICIALIZA LA VARIABLE BUFFER
C
       IF (INIT.EQ.0) THEN
        NBUFFER=(NRECD+IDATP(IFLAG,1))*64 ! el tres es una aproximación muy subestimada de lo que va a necesitar más adelante
        IERROR=0
        ALLOCATE(BUFFER(NBUFFER),STAT=IERROR)
        WRITE(LURES,*) "INTIALIZE BUFFER IN datbas.f"
        IF (IERROR.NE.0) 
     .  CALL RUNEND("I CAN'T ALLOCATE BUFFER IN DATBAS")
        INIT=1
       ENDIF
C*****************************************************************       
C
C**** SE COMPRUEBA Y SE REDIMENSIONA SI ES NECESARIO LA VARIABLE BUFFER
C      

      NBUFFERTMP=(NRECD+IDATP(IFLAG,1))
      IF(NBUFFERTMP-1.GT.NBUFFER) THEN
        IERROR=0
        IF (1.5*NBUFFER.GT.NBUFFERTMP)
     .   NBUFFERTMP=1.5*NBUFFER
        
        WRITE(LURES,*) "REALLOCATE BUFFER (from, to, need)",
     .              NBUFFER,NBUFFERTMP,NRECD+IDATP(IFLAG,1)
        WRITE(LURES,*) 'REQUIRED DATABASE MEMORY (PERMANENT) ='
     .                 ,NBUFFERTMP/(1024.0*1024.0)*8.0,'MB'

        ALLOCATE(TEMP(NBUFFERTMP),STAT=IERROR)
        IF (IERROR.NE.0) 
     .   CALL RUNEND("I CAN'T REALLOCATE BUFFER IN DATBAS")
        TEMP(1:NBUFFER)=BUFFER(:)
        CALL MOVE_ALLOC(FROM=TEMP,TO=BUFFER)   
        NBUFFER=NBUFFERTMP
C        CALL RUNEND('ERROR IN DATBAS: MEMORY OUT OF SPACE ') 
      ENDIF   
C****************************************************************  
       CALL RESTVM(BUFFER,ITASK,ARRAY,IDATP(IFLAG,1),NRECD)
C
      ENDIF
C
      RETURN
      END
