SUBROUTINE elemnt(TASK, name, deltc, ttime, istop, igrav, ivect, lnod,   &
                  flag1, flag2)
!********************************************************************
! Standard ELEMENT routine
!********************************************************************
USE ctrl_db, ONLY: npoin, top, bottom
USE esets_db,ONLY: esets, eset, nelms
USE npo_db
IMPLICIT NONE

  !--- Dummy variables
  CHARACTER(len=*),INTENT(IN):: TASK
    !Optional parameters
  LOGICAL,OPTIONAL:: flag1, flag2
  CHARACTER(len=*),OPTIONAL:: name
  INTEGER(kind=4),OPTIONAL:: istop, igrav, ivect(:)
  INTEGER (kind=4), POINTER, OPTIONAL :: lnod(:,:)
  REAL(kind=8),OPTIONAL:: deltc, ttime
  !--- Local variables
  INTEGER(kind=4):: iset, i

  INTERFACE
    INCLUDE 'elemt1.h'
    INCLUDE 'elemt2.h'
    INCLUDE 'elem04.h'
    INCLUDE 'elemt8.h'
    INCLUDE 'elem10.h'
    INCLUDE 'elem13.h'
    INCLUDE 'elem18.h'
  END INTERFACE

  IF (TRIM(task) == 'RESVPL') THEN
    IF (bottom) coorb=0d0
    IF (top) coort=0d0
    IF (top .OR. bottom) ifact=0
  END IF

  IF (PRESENT(flag2) .AND. TASK /= 'OUTDYN') flag2=.FALSE.

  DO iset=1,esets
    SELECT CASE(eset(iset))
    CASE (1)
      CALL elemt1(TASK, nelms, name, deltc, ttime, flag1, flag2)
    CASE (2)
      CALL elemt2(TASK, nelms, name, deltc, ttime, istop, flag1, flag2)
    CASE (4)
      CALL elem04(TASK, nelms, name, deltc, ttime, istop, &
                  lnod, flag1, flag2)
    CASE (8)
      CALL elemt8(TASK, nelms, name, deltc, ttime, istop, &
                  flag1, flag2)
    CASE (10)
      CALL elem10(TASK, nelms, name, istop, flag2)
    CASE (13)
      CALL elem13(TASK, nelms, name, deltc, ttime, istop, ivect, flag1, flag2)
    CASE (18)
      CALL elem18(TASK, nelms, name, deltc, ttime, istop, &
                  ivect, lnod, flag1, flag2)
    END SELECT

    IF (TRIM(task) == 'BOUNDA' .OR.         &
        TRIM(task) == 'DELETE' .OR.         &
        TRIM(task) == 'INIGAU' .OR.         &
        TRIM(task) == 'NODSET' .OR.         &
        TRIM(task) == 'RFSTBC' .OR.         &
        TRIM(task) == 'SEARCH' .OR.         &
        TRIM(task) == 'SMOOTH' .OR.         &
        TRIM(task) == 'SURFAC' ) THEN
      IF (flag2) THEN
        IF (TRIM(task) == 'SEARCH')  igrav=eset(iset)   !=etype (type of element is 'sname')
        EXIT        !set found, exit loop
      END IF
    END IF

  END DO
  IF (TRIM(task) == 'RESVPL') THEN
    IF (bottom) THEN
      DO i=1,npoin
        IF (ifact(i) > 0) THEN
          coorb(:,i) = coorb(:,i)/ifact(i) + coora(:,i)
        ELSE
          coorb(:,i) = coora(:,i)
        END IF
      END DO
    END IF
    IF (top) THEN
      DO i=1,npoin
        IF (ifact(i) > 0) THEN
          coort(:,i) = coort(:,i)/ifact(i) + coora(:,i)
        ELSE
          coort(:,i) = coora(:,i)
        END IF
      END DO
    END IF
  END IF

RETURN
END SUBROUTINE elemnt
