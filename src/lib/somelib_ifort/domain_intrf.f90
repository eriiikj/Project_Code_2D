module domain_intf

use mater_intf
use elem_intf


! legru borde vara definerad så här
type domain_data
   integer              :: first_ele
   integer              :: last_ele
   integer              :: nr_nodes_ele
   character(len=80)    :: el_name
   character(len=80)    :: mat_name
end type

type(domain_data), allocatable    :: dom(:)



  interface domain_stiff
     module procedure domain_stiff1
     module procedure domain_stiff2
  end interface
  private domain_stiff1, domain_stiff2


!-----------------------------------------------------------------------------
contains

subroutine domain_init(legru,mp)
  implicit none
  integer                       :: legru(:,:)
  double precision              :: mp(:,:)
!  double precision              :: ep ! thickness for 2D

! same number of rows in legru and ep, i.e. one ep
! for each element group  
! initiate different material models

  call mater_init(legru,ep)
  call ele_init(legru,mp)

end subroutine ele_init

subroutine domain_stiff(ke,coord,enod,ed)
! calculate stiffness matrix for whole domain

end subroutine domain_stiff

subroutine domain_stiff_ele1(ke,coord,enod,ed,ie)
! stiffness matrix for one element
! 3D and axi-symmetry
  implicit none
  double precision               :: ke(:,:), coord(:,:),ed(:,:)
  integer                        :: enod(:,:), ie
  
! D, es should be found from material routine

  ielm=elGtoLm(ie)
  iele=elGtoLe(ie)
  call mater_get_D(D,ie)
  call mater_get_stress(es,ie)
  call ele_stiff(ke,coord,enod,ep,D,ed,es,ie) 
  
end subroutine domain_stiff_ele1

subroutine domain_stiff_ele2(ke,coord,ep,enod,ed,ie)
! stiffness matrix for one element
! 2D situation ep (thickness is present)
  implicit none
  double precision               :: ke(:,:), coord(:,:),ed(:,:)
  integer                        :: enod(:,:), ie
  
! D, es should be found from material routine

  ielm=elGtoLm(ie)
  iele=elGtoLe(ie)
  call mater_get_D(D,ie)
  call mater_get_stress(es,ie)
  call ele_stiff(ke,coord,ep,enod,ep,D,ed,es,ie) 
  
end subroutine domain_stiff_ele2

end module domain_intf
