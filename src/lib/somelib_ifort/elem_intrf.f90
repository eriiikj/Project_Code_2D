module elem_intf

! kanske skall det vara en ännu övergripande 
! module domian_intf

use domain_intf:: only domain_data


integer        :: elGtoL(:)    
private elGtoL

!  interface ele_stiff
!     module procedure name
!  end interface
!  private name

! LEGRU    DEFINING ELEMENT GROUP IGR                  
!    LEGRU(IGR,1) FIRST ELEMENT NUMBER          
!    LEGRU(IGR,2) LAST ELEMENT NUMBER           
!    LEGRU(IGR,3) NUMBER OF NODES IN AN ELEMENT 
!    LEGRU(IGR,4) ELEMENT IDENTIFIER
!    LEGRU(IGR,5) MATERIAL MODEL
!
! LEGRU(IGR,4)
!     3 : 3 node element
!     4 : 4 node element
!     8 : 8 node element
! vad med axi och 3D element
! och vad med diffusion och magnetism
!

!-----------------------------------------------------------------------------
contains

subroutine ele_init(legru,ep,mp)
  implicit none
  integer                       :: legru(:,:)
  double precision              :: ep   ! thickness if 2D

  nrgru=size(legru,1)
  


end subroutine ele_init


subroutine ele_stiff(ke,coord,enod,ed,ie)
  implicit none
  double precision               :: ke(:,:), coord(:,:),ed(:,:)
  integer                        :: enod(:,:), ie
  
! calls do differnt elements

! D, es should be found from material routine

  call mater_get_D(D,ie)
  call mater_get_stress(es,ie)
  call ele_stiff(ke,coord,enod,ep,D,ed,es,ie)
  
end subroutine ele_stiff

end module elem_intf
