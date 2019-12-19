 SUBROUTINE dynamic(actio,ircon,nrest,kstep,istop,ecdis,actmsh,endmsh)
 !***********************************************************************
 !
 !                   EXPLICIT DYNAMIC strategy
 !
 !***********************************************************************
   ! Global data bases
   USE param_db,ONLY : milb                       !Global parameters
   USE ctrl_db, ONLY : ncdlt,istep,dtime,ttime,endtm,numct,neq,itemp,therm,npoin  !global control parameters
   USE kinc_db, ONLY : nvelr   !kinematic constraints
   USE loa_db,  ONLY : ifloa                                  !load sets
   USE npo_db,  ONLY : resid,fcont,velnp,tresi,tmass    !point information
   USE outp_db, ONLY : toutd,iwrit,cpui    !output asked, frequencies
   USE gvar_db, ONLY : updiv
   USE lispa0
   IMPLICIT NONE

   !Functions

   !Variables
   CHARACTER(len=milb) :: actio
   INTEGER  (kind=4) :: ircon,   &  ! counter of restart files
                        nrest,   &  ! restart dumping period
                        kstep,i
   INTEGER  (kind=4) :: istop, &  ! termination/error flag
                        ecdis     ! control DOF for output
   LOGICAL ::           actmsh   &  ! .TRUE. if mesh modification required in step
                       ,endmsh      ! .TRUE. if mesh modification required at end step

   LOGICAL ::           newid    &  ! .TRUE. if a new ID array needed
                       ,endst    &  ! .TRUE. if End Time of Strategy found
                       ,foutp    &  ! .TRUE. if Output Necessary
                       ,tflag    &  ! .TRUE. if reset time cpu forecast
                       ,actrst      ! .TRUE. if restart file must be written


   INTERFACE
     INCLUDE 'contac.h'
     INCLUDE 'elemnt.h'
     INCLUDE 'outdyn.h'
     INCLUDE 'screen.h'
   END INTERFACE

   newid = .TRUE.    ! compute ID in explit
   tflag = .TRUE.    ! Initialize CPU time forecast

   DO  ! loop over time steps

      !         Calculation of critical time step
      IF( ncdlt > 0 )THEN
        IF (MOD(istep,ncdlt) == 0 .OR. TRIM(actio) /= 'NSTRA0' ) THEN

          CALL timing(17,1)
          CALL incdlt ( )
          CALL timing(17,2)
        END IF
      END IF

      ! PREDICTION OF CPU TIME
      istep = istep + 1
      IF (MOD(istep,kstep)==0 .OR. tflag) CALL screen(istep,dtime,ttime,cpui,endtm,tflag)

      !           Calculation of contact forces
      IF(numct > 0)THEN

        CALL timing(16,1)
        DO i=1,npoin
          fcont(:,i) = 0d0
        END DO
        CALL contac('FORCES',iwrit,dtcal=dtime, ttime=ttime, velnp=velnp)

        CALL timing(16,2)
      END IF

      !           Nonconservative loading

      IF (ifloa > 0) CALL loadfl(ttime+dtime)
      !           Time integration

      CALL timing(13,1)
      CALL explit(newid)  ! standard central difference scheme
      !IF( itemp ) CALL texplit( ) !CALL heatstep ( dtime, ttime)
      ttime = ttime + dtime
      endst = (ttime+dtime/3d0 >= endtm)

      CALL timing(13,2)

      !           Calculation of internal nodal forces

      CALL timing(14,1)
      resid = 0d0  !; updiv = .FALSE.
      CALL elemnt('RESVPL', deltc=dtime, istop=istop, ttime= ttime)
      !IF( therm ) THEN
      !  tresi = 0d0
      !  CALL elemnt ('TRESID', deltc=dtime, ttime= ttime)
      !  ! frequency to recompute capacity matrix
      !  !CALL tlmass ( )              !compute capacity matrix
      !  !tmass = tmass*tmscal         !scale capacity matrix
      !END IF

      CALL timing(14,2)
      IF (istop == 1) THEN !istop=1 -> there was no convergence
         CALL outdyn(.FALSE.,.TRUE., ecdis)  ! final state
         WRITE(6 ,900)
         WRITE(55,900,err=9999)
      END IF

      !           Verify mesh modification criteria
      CALL mshcri(istep, ttime, dtime, actmsh, endmsh)
      IF (actmsh .OR. endmsh) endst=.TRUE.

      foutp = ( istop /= 0 .OR. endst )
      CALL timing(15,1)

      !           Results output

      CALL outdyn(.FALSE., foutp, ecdis)
      CALL timing(15,2)

      !           Restart file dumping
      IF (nrest > 0) THEN
        CALL timerst(.FALSE.,nrest,actrst)
      ELSE IF (nrest < 0) THEN
        actrst = MOD(istep,-nrest)==0
      ELSE
        actrst = .FALSE.
      END IF

      IF (actrst .OR. endst ) THEN
        CALL timing(15,1)
        CALL dumpin(ircon,endst,actmsh,endmsh)
        CALL timing(15,2)
        IF (endst.AND.(nrest > 0)) CALL timerst(.TRUE.,nrest,actrst)   !Initialize time for writting restart file
      END IF
      actio = 'NSTRA0'

      IF (istop /= 0 .OR.  endst ) EXIT

   END DO

   RETURN
 9999 CALL runen2('ERROR WHILE WRITING "DEB" FILE.')
 ! ++++FORMATS ++++++++++++++++++
  900 FORMAT(/,12x,' DIVERGENCE IN THE RETURN ALGORITHM : ', &
                   'PROGRAM STOPPED',/,17x,' FINAL GEOMETRY ', &
                   'SAVED FOR POSTPROCESSING!')
 END SUBROUTINE dynamic
