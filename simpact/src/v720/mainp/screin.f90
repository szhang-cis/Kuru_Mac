SUBROUTINE screin(task,cpui)
!---------------------------------------
!  Write out headers
!---------------------------------------
USE param_db,ONLY: mich
USE name_db,ONLY : prognm
IMPLICIT NONE

  !--- Dummy variables
  REAL(kind=8),INTENT(OUT) :: cpui
  CHARACTER(len=*),INTENT(IN) :: task
  !--- Functions
  CHARACTER(len=mich):: inttoch
  !--- Local variables
  INTEGER(kind=4):: i, lng, fch(3), hra(3)

  INCLUDE 'prog_data.fpp'

  CALL timuse(cpui)
  lng = LEN_TRIM(prognm)
  WRITE(6,"(///////33X,<lng>(A,1X),/)") (prognm(i:i),i=1,lng)
  WRITE(55,"(33X,<lng>(A,1X),/)",ERR=9999) (prognm(i:i),i=1,lng)
  !PROGRAM version number and date
  WRITE( 6,"(33X,'version: ',A/33X,'date   : ',A,/)")TRIM(version),TRIM(verdate)
  WRITE(55,"(33X,'version: ',A/33X,'date   : ',A,/)",ERR=9999) TRIM(version),TRIM(verdate)
  DO i=1,7
    WRITE( 6,"(t15,a50)")titles(i)
    WRITE(55,"(t15,a50)",ERR=9999)titles(i)
  END DO
!ms$if (acad > 0)
  WRITE( 6,"(16X,'                Academic version         ',/)")
  WRITE(55,"(16X,'                Academic version         ',/)",ERR=9999)
!ms$endif
  WRITE (6,"(/)")
  WRITE (55,"(/)",ERR=9999)

  CALL idntdt(fch,hra)
  IF (TRIM(task) == 'BEGINP') THEN
    WRITE(6,101) TRIM(inttoch(hra(1),2)),TRIM(inttoch(hra(2),2)),TRIM(inttoch(hra(3),2)),  &
                 TRIM(inttoch(fch(3),2)),TRIM(inttoch(fch(2),2)),TRIM(inttoch(fch(1),2))
    WRITE(55,101,ERR=9999) TRIM(inttoch(hra(1),2)),TRIM(inttoch(hra(2),2)),TRIM(inttoch(hra(3),2)), &
                           TRIM(inttoch(fch(3),2)),TRIM(inttoch(fch(2),2)),TRIM(inttoch(fch(1),2))
  ELSE IF (TRIM(task) == 'RESTAR') THEN
    WRITE(6,102) TRIM(inttoch(hra(1),2)),TRIM(inttoch(hra(2),2)),TRIM(inttoch(hra(3),2)),  &
                 TRIM(inttoch(fch(3),2)),TRIM(inttoch(fch(2),2)),TRIM(inttoch(fch(1),2))
    WRITE(55,102,ERR=9999) TRIM(inttoch(hra(1),2)),TRIM(inttoch(hra(2),2)),TRIM(inttoch(hra(3),2)), &
                           TRIM(inttoch(fch(3),2)),TRIM(inttoch(fch(2),2)),TRIM(inttoch(fch(1),2))
  END IF

RETURN
  101 FORMAT(23X,'Started at  ',2(A,':'),A,'  on  ',2(A,'/'),A,//)
  102 FORMAT(23X,'RE-Started at  ',2(A,':'),A,'  on  ',2(A,'/'),A,//)
 9999 CALL runen2('Problems detected in "screin".')
END SUBROUTINE screin

SUBROUTINE screen(istep,dtime,ttime,cpui,endtm,tflag)
USE param_db,ONLY: mich
IMPLICIT NONE

  !--- Dummy arguments
  INTEGER(kind=4),INTENT(IN) :: istep
  REAL(kind=8),INTENT(IN) :: dtime, ttime, cpui, endtm
  LOGICAL,INTENT(INOUT):: tflag
  !Funcitons
  CHARACTER(len=mich):: inttoch
  !--- Local variables
  REAL(kind=8) :: cpuc, cpua, cucuc
  CHARACTER(len=120):: LineOUT
  INTEGER(kind=4) :: seg, iha, imi, isg, seg1, iha1, imi1, isg1, kstep
  INTEGER(kind=4),SAVE:: stepo=0 !Step old (previous screen print
  REAL (kind=8), SAVE :: cpuo !CPU old (previous screen print)

  IF (tflag) THEN
    tflag = .FALSE.
    stepo = istep
    CALL timuse(cpuo)
  ELSE
    CALL timuse(cpua)                 !present time
    cpuc = cpua - cpui                !time elapsed since start
    IF (cpuc < 0) cpuc=cpuc+86400     !if one day has passed
    seg = INT(cpuc)                   !round to seconds
    isg = MOD(seg,60)                 !seconds  
    imi = MOD(seg/60,60)              !minutes
    iha = seg/3600                    !hours
    kstep = istep - stepo             !steps since previous screen print
    cpuc = cpua - cpuo                !CPU since previous screen print
    cucuc = (endtm-ttime)/dtime*cpuc/kstep
    seg1 = INT(cucuc)        !round to seconds
    isg1 = MOD(seg1,60)      !seconds 
    imi1 = MOD(seg1/60,60)   !minutes
    iha1 = seg1/3600         !hours
    WRITE(LineOUT,"('Step ',A,2x,'DT=',e10.3,2x,'Time=',e10.4,2x,'CPU ',A,2(':',A),2x,         &
          'CPU to end st: ',A,2(':',A2))") TRIM(inttoch(istep,0)),dtime,ttime,      &
          TRIM(inttoch(iha,0)), TRIM(inttoch(imi,2)), TRIM(inttoch(isg,2)),                    &
          TRIM(inttoch(iha1,0)),TRIM(inttoch(imi1,2)),TRIM(inttoch(isg1,2))
    WRITE(6,"(A,A,$)") CHAR(13), TRIM(LineOUT)
    WRITE(55,"(' STEP ',A,3X,'DT= ',e11.4,3x,'Time= ',e11.5,3x,'CPU ',A,2(':',A),2x,      &
          'CPU to end stage (est.) ',A,2(':',A))",ERR=9) TRIM(inttoch(istep,0)),dtime,ttime, &
          TRIM(inttoch(iha,0)),TRIM(inttoch(imi,2)),TRIM(inttoch(isg,2)),                    &
          TRIM(inttoch(iha1,0)),TRIM(inttoch(imi1,2)),TRIM(inttoch(isg1,2))
    stepo = istep
    cpuo = cpua
!    !Variables para escribir el informe
!    CPU(1) = iha1
!    CPU(2) = imi1
!    CPU(3) = isg1
!    CALL wrtrpt('WTIME',newnum=istep,nrea=dtime,nint=3,vint=CPU)
  END IF

RETURN
 9 CALL runen2('ERROR WHILE WRITING "DEB" FILE.')
END SUBROUTINE screen

SUBROUTINE flushf(n)
INTEGER (kind = 4) n
!     dummy routine (in main-frame fortran it clears buffers for file n)
!      CALL flush(n)           !Not ALL FORTRANs
END SUBROUTINE flushf
