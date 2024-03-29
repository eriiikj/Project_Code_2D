module wrt2vtk

! =======================================================================================================
!
!   Module wrt2vt contains subroutines for writing Fortran data to vtk files.
!   Only unstructured grid is implemented
!
!   Included routines:
!
!      VTKopen                    Open a VTK file and assign a pointer and a filename
!      VTKclose                   Close a VTK file
!      VTKPlotMesh                Creates an unstructured mesh (2D or 3D)
!      VTKPlotNodeVal             Write nodal values to an unstructured mesh
!      VTKPlotCellVal             Write cell values to an unstructured mesh
!
!   Last rev.:
!      M. Ristinmaa, 2020-02-25
!
! =======================================================================================================

implicit none

integer    :: filenumb=10
private filenumb

interface VTKOpen
   module procedure VTKopen1
   module procedure VTKopen2
end interface
private VTKopen1, VTKopen2

contains

! =======================================================================================================

subroutine VTKopen1(filename)
   ! ===========================================================
   ! This routine opens VTK file.
   !
   !    Input:
   !          filename          Name of the new VTK file
   !
   ! Last rev.:
   !    2020-02-25: M. Ristinmaa
   ! ===========================================================

   implicit none

   character(len=*)         :: filename

   open(unit=filenumb,file=filename,status='UNKNOWN')
   write(filenumb,"(A26)")'# vtk DataFile Version 2.0'
   write(filenumb,"(A14)")'SomeLib output'
   write(filenumb,"(A5)")'ASCII'
   write(filenumb,"(A25)")'DATASET UNSTRUCTURED_GRID'

   return

end subroutine VTKopen1

subroutine VTKopen2(filename,num)
   ! ===========================================================
   ! This routine opens VTK file.
   !
   !    Input:
   !          filename          Name of the new VTK file
   !          num               extension to VTK file filename.num
   !                            used for time series
   !
   ! Last rev.:
   !    2020-02-25: M. Ristinmaa
   ! ===========================================================

   implicit none

   character(len=*)         :: filename
   character(len=80)        :: filee
   character(len=30)        :: str
   integer                  :: num

!write(ci,'(I10)') step
!write(filename,'(A5,A10)') 'step_',adjustl(ci)
!filename = trim(filename) // '.vtk'
!open(unit=ivtk,file=filename,form='binary',convert='BIG_ENDIAN')

   write(str,*)num
   filee=trim(filename)//'.'//adjustl(str)
   open(unit=filenumb,file=filee,status='UNKNOWN')
   write(filenumb,"(A26)")'# vtk DataFile Version 2.0'
   write(filenumb,"(A14)")'SomeLib output'
   write(filenumb,"(A5)")'ASCII'
   write(filenumb,"(A25)")'DATASET UNSTRUCTURED_GRID'

   return

end subroutine VTKopen2

subroutine VTKclose()
   ! ===========================================================
   ! This routine closes VTK file.
   !
   ! Last rev.:
   !    2020-02-25: M. Ristinmaa
   ! ===========================================================

   implicit none

   close(unit=filenumb)

   return

end subroutine VTKclose

subroutine VTKplotMesh(eltype,enod,coord)

   implicit none

   character                           :: eltype*(*)
   integer, allocatable                :: enod(:,:)
   double precision, allocatable       :: coord(:,:)
   integer                             :: ie, im, nrelem, nrnode, nodeel, vtktype
   character(len=10)                   :: str

   nrelem=size(enod,2)
   nodeel=size(enod,1)
   nrnode=size(coord,2)

  ! write(str,*)nrnode
  !write(*,*)'Number of nods                :',maxval(enod), size(coord,1), size(coord,2)
  !write(*,*)'Number of elements            :',size(enod,2)

! write nodes
   write(filenumb,"(A7,X,I7,X,A6)")'POINTS ',nrnode,' float'
   if (size(coord,1).eq.3) then
     do im=1,nrnode
        write(filenumb,"(3f12.5)")coord(1,im),coord(2,im),coord(3,im)
     enddo
   else
     do im=1,nrnode
        write(filenumb,"(3f12.5)")coord(1,im),coord(2,im),0.d0
     enddo
   endif

! write cells or elements
   write(filenumb,"(A6,I7,X,I9)")'CELLS ',nrelem,nrelem*(nodeel+1)
   do ie=1,nrelem
       write(filenumb,"(40(I6,X))")nodeel,(enod(im,ie)-1,im=1,nodeel)
   enddo

! write vtk cell type
   if (eltype.eq.'tr3') vtktype=5
   if (eltype.eq.'qu4') vtktype=9
   if (eltype.eq.'brick8') vtktype=12
   write(filenumb,"(A11,X,I7)")'CELL_TYPES ',nrelem
   do ie=1,nrelem
       write(filenumb,"(I6)")vtktype
   enddo

end subroutine VTKPlotMesh

subroutine VTKPlotNodeVal(val,name)
   double precision                    :: val(:)
   character                           :: name*(*)
   integer                             :: nrnode, im
   CHARACTER(LEN=80) :: String

   nrnode=size(val)
   write(filenumb,"(A11,X,I7)")'POINT_DATA ',nrnode
   write(string,*)'SCALARS ',name,' float'
   !write(filenumb,"(A11,X,A6)")'SCALARS ',name,' float'
   write(filenumb,"(A80)")adjustl(string)
   write(filenumb,"(A20)")'LOOKUP_TABLE default'
   do im=1,nrnode
 !      write(filenumb,"(f12.5)")val(im)
       write(filenumb,*)val(im)
   enddo

end subroutine VTKPlotNodeVal

end module wrt2vtk
