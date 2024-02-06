module mater_magnlin

! Last modified
! M. Ristinmaa 2021-02-24
! M. Ristinmaa 2021-03-17
!  added new Bintr output to magmater_getlinH
!
!-----------------------------------------------------------------------------

!use memory_util

implicit none

 ! double precision, allocatable     :: mu(:,:), si(:,:)
  double precision, allocatable     :: propEM(:,:)
  double precision, allocatable     :: Hsave(:,:,:)
  double precision, allocatable     :: Bsave(:,:,:)

  integer, allocatable              :: matId(:)

  double precision, parameter       :: mu0=1.256637061435917d-06
  integer                           :: nrelem, nrgp
 
  private propEM, matID, Hsave, Bsave
  private mu0, nrelem, nrgp


interface magmater_tangentlin
   module procedure magmater_permelin_get1
   module procedure magmater_permelin_get2
end interface
private :: magmater_permelin_get1, magmater_permelin_get2

!-----------------------------------------------------------------------------
contains


subroutine magmater_initlin(propEMin,matIdin,nip)
implicit none

double precision                  :: propEMin(:,:)
integer                           :: matIdin(:)
integer                           :: nip            ! integration points

integer                           :: nels, ierr, iel, nmat, nprp

  nels=size(matIdin)
  nmat=size(propEmin,1)
  nprp=size(propEMin,2)
  
  nrelem=nels
  nrgp  =nip

!  write(*,*)'mater_magnlin'
!  write(*,*)'nrelem ',nels
!  write(*,*)'nrmat  ',size(propEM,1)
!  write(*,*)'nrprop ',size(propEM,2)
!  write(*,*)'prop   ',propEM(1,:)
!  write(*,*)'prop   ',propEM(2,:)
!  write(*,*)'matId  ',matId

  allocate(Hsave(3,nip,nels), stat=ierr)
  IF (ierr /= 0) STOP "*** Not enough memory ***"
  allocate(Bsave(3,nip,nels), stat=ierr)
  IF (ierr /= 0) STOP "*** Not enough memory ***"
  allocate(matId(nels), stat=ierr)
  allocate(propEM(nmat,nprp), stat=ierr)

  Hsave=0d0
  Bsave=0d0
  propEM=propEMin
  matId=matIdin


end subroutine magmater_initlin


subroutine magmater_magfield_get(H,B,ie) 
implicit none

double precision                  :: H(:,:), B(:,:)

double precision                  :: mu
integer                           :: ie, i

  mu=mu0*propEM(matId(ie),2)
  do i=1,nrgp
    H(:,i)=B(:,i)/mu
  end do
  
  if (size(H,1).eq.3) then
    Hsave(:,:,ie)=H
    Bsave(:,:,ie)=B
  else
    Hsave(1,:,ie)=H(1,:)
    Hsave(2,:,ie)=H(2,:)
    Hsave(3,:,ie)=0d0
    Bsave(1,:,ie)=B(1,:)
    Bsave(2,:,ie)=B(2,:)
    Bsave(3,:,ie)=0d0
  end if  

end subroutine magmater_magfield_get


subroutine magmater_permelin_get1(mu_out,si_out,iel)
implicit none

double precision                  :: mu_out(:,:,:), si_out(:,:,:)
integer                           :: iel, ig, dim, Id(3,3)
double precision                  :: mu, si

  dim=size(mu_out,1)
  Id(1,:)=(/1d0, 0d0, 0d0/)
  Id(2,:)=(/0d0, 1d0, 0d0/)
  Id(3,:)=(/0d0, 0d0, 1d0/)

  si=propEM(matId(iel),1)
  mu=mu0*propEM(matId(iel),2)
  do ig=1,nrgp
    mu_out(:,:,ig)=Id(1:dim,1:dim)/mu
    si_out(:,:,ig)=Id(1:dim,1:dim)*si
  end do  

  
  
end subroutine magmater_permelin_get1


subroutine magmater_getlinH(field,Hout,iel)
implicit none

character                         :: field*(*)
double precision                  :: Hout(:,:)
integer                           :: iel


  if (trim(field).eq.'Hfield') then
    if (size(Hout,1).eq.3) then
      Hout=Hsave(:,:,iel)
    else
      Hout=Hsave(1:2,:,iel)
    end if
  elseif (trim(field).eq.'Bfield') then
    if (size(Hout,1).eq.3) then
      Hout=Bsave(:,:,iel)
    else
      Hout=Bsave(1:2,:,iel)
    end if
  elseif (trim(field).eq.'Bintr') then
    Hout=0d0
  else
    write(*,*)'Field is note stored in mater_magnlin'
    write(*,*)field
    stop
  end if


end subroutine magmater_getlinH


subroutine magmater_permelin_get2(mu_out,si_out,iel)
implicit none

double precision                  :: mu_out(:), si_out(:)
integer                           :: iel

   si_out=propEM(matId(iel),1)
   mu_out=mu0*propEM(matId(iel),2)

end subroutine magmater_permelin_get2




end module mater_magnlin
