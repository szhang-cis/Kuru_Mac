SUBROUTINE elmdat(task,nelem,elsnam,itype)

! read element sets

  USE ctrl_db, ONLY: dtuser, ndime, ndofn, npoin, neulr
  USE outp_db, ONLY: iwrit
  USE esets_db, ONLY: nelms,nel
  USE npo_db
  IMPLICIT NONE

  CHARACTER(len=*),INTENT(IN) :: elsnam
  CHARACTER(len=*),INTENT(IN) :: task
  INTEGER (kind=4),INTENT(IN) :: itype
  INTEGER (kind=4),INTENT(IN OUT) :: nelem

  INTERFACE
    INCLUDE 'inpda1.h'
    INCLUDE 'inpda2.h'
    INCLUDE 'inpd04.h'
    INCLUDE 'inpda8.h'
    INCLUDE 'inpd10.h'
    INCLUDE 'inpd13.h'
    INCLUDE 'inpd18.h'
  END INTERFACE

  SELECT CASE (itype)
  CASE (1)
    CALL inpda1(task,ndime,neulr,nelem,iwrit,elsnam,nelms(1))
  CASE (2)
    CALL inpda2(task,nelem,iwrit,elsnam,nelms(2))
  CASE (4)
    CALL inpd04(task,iwrit,elsnam,nelms(4))
  CASE (8)
    CALL inpda8(task,nel(8),eule0,euler,coord,iwrit,elsnam,nelms(8))
  CASE (10)
    CALL inpd10(task,ndime,nelem,iwrit,elsnam,nelms(10))
  CASE (13)
    CALL inpd13(task,iwrit,elsnam,nelms(13))
  CASE (18)
    CALL inpd18(task,nel(18),iwrit,elsnam,nelms(18))
  CASE DEFAULT
    CALL runend('ELMDAT: ELEMENT_TYPE NOT EXISTENT  ')
  END SELECT

  RETURN
END SUBROUTINE elmdat
