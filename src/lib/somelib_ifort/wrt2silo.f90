MODULE wrt2silo
  !
  !------------------------------------------------------------------------
  !
  ! Module wrt2silo contains subroutines for writing Fortran data to Silo files
  ! The routines requires the Silo libraries to be installed.
  ! These can be obtained from:
  ! https://wci.llnl.gov/codes/visit/3rd_party
  !
  ! Included routines:
  !
  !      createSiloDb               Opens a Silo file and assigns a pointer
  !      closeSiloDb                Closes a Silo file
  !      mkQmeshCollinear           Creates a quad mesh collinear with the
  !                                 coordinate axes (2D or 3D)
  !
  !      mkQmeshNonCollinear        Creates a quad mesh non-collinear with
  !                                 the coordinate axes (2D or 3D)
  !
  !      setQmeshZoneVar            Assigns a zone variable to a quad mesh
  !                                 (2D or 3D, integer or double)
  !
  !      setCurve                   Creates a xy-curve (integer or double)
  !
  !   Last rev.:
  !      H. Hallberg, 2010-09-01
  !
  !------------------------------------------------------------------------
  !
  !USE precision

  IMPLICIT NONE
    INTEGER, PARAMETER::dp=SELECTED_REAL_KIND(15)
    !INCLUDE "silo.inc"
  INCLUDE "./silo/include/silo_f9x.inc"
  !
  ! Make module methods inaccessible to outside calls
  PRIVATE::mkQmesh2dUnstr, mkQmesh3dUnstr
  !
  ! Declare accessible module interfaces
  INTERFACE mkQmeshUnstr
     MODULE PROCEDURE mkQmesh2dUnstr
     MODULE PROCEDURE mkQmesh3dUnstr
  END INTERFACE
  !

  interface openSilo
     module procedure createSiloDb1
     !module procedure createSiloDb2
  end interface

  interface closeSilo
     module procedure closeSiloDb1
    ! module procedure closeSiloDb2
  end interface

CONTAINS
  !

  !------------------------------------------------------------------------
  !
  SUBROUTINE createSiloDb1(dbfile,siloname)
    !------------------------------------------------------------------------
    ! This routine closes an open Silo file.
    !
    !    Input:
    !          siloname          Name of the new Silo file
    !
    !    Output:
    !          dbfile            Pointer to the open Silo file
    !
    ! Last rev.:
    !    2010-09-01: H. Hallberg
    !------------------------------------------------------------------------

    IMPLICIT NONE
    CHARACTER(*),INTENT(IN)  :: siloname
    INTEGER, INTENT(OUT)     :: dbfile
    INTEGER                  :: ierr

    ! Create Silo file
    ierr = dbcreate(TRIM(siloname),LEN_TRIM(siloname),DB_CLOBBER,DB_LOCAL,&
         &DB_F77NULL,0,DB_PDB,dbfile)
    IF (ierr.EQ.-1) THEN
       WRITE(*,*) 'Could not create Silo file ', trim(siloname),&
            &'. Program execution stopped!'
       STOP
    END IF

    RETURN
  END SUBROUTINE createSiloDb1
  !
  !------------------------------------------------------------------------
  !
  SUBROUTINE closeSiloDb1(dbfile)
    ! =====================================================================
    ! This routine closes an open Silo file.
    !
    !    Input:
    !          dbfile            Pointer to the open Silo file
    !          outfileid         File-ID for program output
    !
    !    Output:
    !
    ! Last rev.:
    !    2010-09-01: H. Hallberg
    ! ======================================================================

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: dbfile
    INTEGER             :: ierr

    ierr = dbclose(dbfile)
    IF (ierr.eq.-1) THEN
       WRITE(*,*) 'Could not close Silo file. Program execution stopped!'
       STOP
    END IF

    RETURN

  END SUBROUTINE closeSiloDb1
  !
  !----------------------------------------------------------------------------
  !
  SUBROUTINE openExisting(siloname,dbfile)
    IMPLICIT NONE
    CHARACTER(*),INTENT(IN)::siloname
    INTEGER, INTENT(OUT)::dbfile
    INTEGER::ierr
!
    ierr=dbopen(TRIM(siloname),LEN_TRIM(siloname),DB_PDB,DB_APPEND,dbfile)
    IF (ierr.eq.-1) THEN
       WRITE(*,*) 'Could not open existing solo file.&
             & Program execution stopped!'
       STOP
    END IF
  END SUBROUTINE openExisting
  !
  !----------------------------------------------------------------------------
  !
  SUBROUTINE mkQmesh2dUnstr(dbfile,meshName,element,x,y,nels,nodelist,&
       &cycle,time)
    IMPLICIT NONE
    INTEGER, INTENT(IN)::dbfile,nodelist(:),nels,cycle
    REAL(dp),INTENT(IN)::x(:),y(:),time
    CHARACTER(*), INTENT(in):: meshname,element
    INTEGER::ierr,err,ndims,lnodelist,shapesize(1),shapecounts(1),&
         &nshapetypes,nnodes,optlistid,shapetype
    !
    ndims=2
    nnodes=UBOUND(x(:),1)
    lnodelist=UBOUND(nodelist(:),1)
    shapecounts=nels
    nshapetypes=1

    SELECT CASE(element)
    CASE('QUAD4')
       shapesize = 4
       shapetype = DB_ZONETYPE_QUAD
    CASE DEFAULT
       WRITE(*,*)'error'
    END SELECT
    !
    !----- Write the node connectivities ------
    !
    ierr=dbputzl2(dbfile,"zonelist",8,nels,ndims,nodelist,&
         &lnodelist,1,0,0,shapetype,shapesize,nels,1,DB_F77NULL,err)
    !
    IF (ierr.EQ.-1) then
       WRITE(12,*) 'Could not write unstructered 2D mesh ',&
            & TRIM(meshname), ' to Silo file. Program execution stopped!'
       STOP
    END IF

    !Option list, cykle and time
    ierr = dbmkoptlist(2, optlistid)
    ierr = dbaddiopt(optlistid, DBOPT_CYCLE, cycle)
    ierr = dbaddiopt(optlistid, DBOPT_DTIME, time)

    !Write the x and y coordinate values
    ierr=dbputum(dbfile,meshName,LEN_TRIM(meshName),ndims,x,y,DB_F77NULL,&
         &"x",1,"y",1,DB_F77NULL,0,DB_DOUBLE,nnodes,nels,"zonelist",8,&
         &DB_F77NULL,0,optlistid,err)
    IF (ierr.EQ.-1) then
       WRITE(*,*) 'Could not write unstructered 2D ',&
            & TRIM(meshname), ' to Silo file. Program execution stopped!'
       STOP
    END IF

  END SUBROUTINE mkQmesh2dUnstr
  !
  !------------------------------------------------------------------------
  !
  SUBROUTINE mkQmesh3dUnstr(dbfile,meshName,element,x,y,z,nels,nodelist,&
       &cycle,time)
    IMPLICIT NONE
    INTEGER, INTENT(IN)::dbfile,nodelist(:),nels,cycle
    REAL(dp),INTENT(IN)::x(:),y(:),z(:),time
    CHARACTER(*), INTENT(in):: meshname,element
    INTEGER::ierr,err,ndims,lnodelist,shapesize(1),shapecounts(1),&
         &nshapetypes,nnodes,optlistid,shapetype
    !
    ndims=3
    nnodes=UBOUND(x(:),1)
    lnodelist=UBOUND(nodelist(:),1)
    shapecounts=nels
    nshapetypes=1

    SELECT CASE(element)
    CASE('HEX8')
       shapesize = 8
       shapetype = DB_ZONETYPE_HEX
    CASE DEFAULT
       WRITE(*,*)'error'
    END SELECT
    !
    !----- Write the node connectivities ------
    !
    ierr=dbputzl2(dbfile,"zonelist",8,nels,ndims,nodelist,&
         &lnodelist,1,0,0,shapetype,shapesize,nels,1,DB_F77NULL,err)
    !
    IF (ierr.EQ.-1) then
       WRITE(12,*) 'Could not write unstructered 3D mesh ',&
            & TRIM(meshname), ' to Silo file. Program execution stopped!'
       STOP
    END IF

    !Option list, cykle and time
    ierr = dbmkoptlist(2, optlistid)
    ierr = dbaddiopt(optlistid, DBOPT_CYCLE, cycle)
    ierr = dbaddiopt(optlistid, DBOPT_DTIME, time)

    !Write the x and y coordinate values
    ierr=dbputum(dbfile,meshName,LEN_TRIM(meshName),ndims,x,y,z,&
         &"x",1,"y",1,"z",1,DB_DOUBLE,nnodes,nels,"zonelist",8,&
         &DB_F77NULL,0,optlistid,err)
    IF (ierr.EQ.-1) then
       WRITE(*,*) 'Could not write unstructered 3D ',&
            & TRIM(meshname), ' to Silo file. Program execution stopped!'
       STOP
    END IF

  END SUBROUTINE mkQmesh3dUnstr
  !
  !------------------------------------------------------------------------
  !
  SUBROUTINE putNodeVar(dbfile,meshname,plotVar,nodal)
    IMPLICIT NONE
    INTEGER, INTENT(IN)::dbfile
    REAL(8),INTENT(IN)::nodal(:)
    CHARACTER(*),INTENT(in):: meshname,plotVar
    INTEGER::err, ierr, nnodes

    nnodes=UBOUND(nodal(:),1)
    err = dbputuv1(dbFile,plotVar,LEN_TRIM(plotVar),meshName,&
          LEN_TRIM(meshName),nodal,nnodes,DB_F77NULL,0,DB_DOUBLE,&
          DB_NODECENT,DB_F77NULL,ierr)
    IF (ierr.EQ.-1) then
       WRITE(*,*) 'Could not write unstructered 2D ',&
            & TRIM(meshname), ' to Silo file. Program execution stopped!'
       STOP
    ENDIF

  END SUBROUTINE putNodeVar
  !
  !------------------------------------------------------------------------
  !
  SUBROUTINE defvars(dbfile,names,defs,lnames,ldefs)
    IMPLICIT NONE
    INTEGER, INTENT(IN)::dbfile,lnames(:),ldefs(:)
    CHARACTER(*), INTENT(IN)::names(:),defs(:)
    INTEGER::oldlen,err,ierr,types(3),optlist(3)

    types=(DB_VARTYPE_VECTOR)
    optlist=(DB_F77NULL)

    oldlen=dbget2dstrlen()
    err = dbset2dstrlen(60)
    err = dbputdefvars(dbfile,'Defvars',7,1,names,lnames,types,defs,&
         &ldefs,optlist,ierr)
!
    err = dbset2dstrlen(oldlen)
!
  END SUBROUTINE defvars
!

!
!--------------------plot the mesh of the file--------------------------
!
SUBROUTINE mkMesh(dbfile,inc,time,coord,enod)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::inc, enod(:,:), dbfile
  REAL(dp),INTENT(IN)::time, coord(:,:)
  REAL(dp),ALLOCATABLE::x(:),y(:)
  INTEGER,ALLOCATABLE::nodelist(:)
  INTEGER::i,j,k
  CHARACTER(72)::meshName
  !
  integer   ::nels, nod, nn

  nn   = maxval(enod)
  nels = size(enod,2)
  nod  = size(enod,1)

  ALLOCATE(nodelist(nels*nod),x(nn),y(nn))
  !
  meshName='mesh'

     x(:)=coord(1,:)
     y(:)=coord(2,:)
  !
  j=1
  DO i=1,nels !make the element conectivies to a vector
     k=i*nod
     nodelist(j:k)=enod(1:nod,i)
     j=j+nod
  ENDDO
  !
     CALL mkQmeshUnstr(dbfile,meshName,'QUAD4',x,y,nels,&
         &nodelist,inc,time)
!
END SUBROUTINE mkMesh


END MODULE wrt2silo
