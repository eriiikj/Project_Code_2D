module mater_magnnonl1

! Last modified
! M. Ristinmaa 2021-02-24
! M. Ristinmaa 2021-03-17
!  added new Bintr output to magmater_getnonl1H
!-----------------------------------------------------------------------------
!
! nonlinear model, no Hysteresis
! B=mu0*H+2*Bs/pi*arctan(pi*mu0*(mur-1)/(2*Bs)*H)*(1-exp(T-Tc)/C)
!
use matrix_util, only: inv3 

implicit none

interface magmater_calcnonl1
   module procedure magmater_calcH1
   module procedure magmater_calcH2
end interface
private magmater_calcH1, magmater_calcH2

interface magmater_tangentnonl1
   module procedure magmater_tangent_getH1
   module procedure magmater_tangent_getH2
end interface
private magmater_tangent_getH1, magmater_tangent_getH2



  double precision, allocatable     :: propEM(:,:)
  double precision, allocatable     :: Hsave(:,:,:,:)
  double precision, allocatable     :: Bsave(:,:,:,:)

  integer, allocatable              :: matId(:)
  integer                           :: New, Old

  double precision, parameter       :: mu0=1.256637061435917d-06
  double precision, parameter       :: pi=3.141592653589793d0
  private propEM, matID, mu0, Hsave, Bsave,pi
  private New, Old

!-----------------------------------------------------------------------------
contains 


subroutine magmater_initnonl1(propEMin,matIdin,nip)
implicit none

double precision                  :: propEMin(:,:)
integer                           :: matIdin(:), nip

integer                           :: nels, nprp, ierr, nmat

  nels=size(matIdin)
  nmat=size(propEmin,1)
  nprp=size(propEMin,2)

!  write(*,*)'mater_magnlin'
!  write(*,*)'nrelem ',nels
!  write(*,*)'nrmat  ',size(propEMin,1)
!  write(*,*)'nrprop ',size(propEMin,2)
!  write(*,*)'prop   ',propEMin(1,:)
!  write(*,*)'prop   ',propEMin(2,:)
!  write(*,*)'matId  ',matIdin

  allocate(matId(nels), stat=ierr)
  allocate(propEM(nmat,5), stat=ierr)

  propEM=propEMin
  matId=matIdin

  allocate(Hsave(3,nip,nels,2), stat=ierr)
  IF (ierr /= 0) STOP "*** Not enough memory ***"
  allocate(Bsave(3,nip,nels,2), stat=ierr)
  IF (ierr /= 0) STOP "*** Not enough memory ***"

  Hsave=0d0
  Bsave=0d0
  New=2
  Old=1

end subroutine magmater_initnonl1


subroutine magmater_acceptnonl1
implicit none

!  if (Old.eq.1) then
!    New=1 
!    Old=2
!  else
!    New=2 
!    Old=1
!  end if

   Hsave(:,:,:,Old)=Hsave(:,:,:,New)


end subroutine magmater_acceptnonl1


subroutine magmater_tangent_getH1(mui,si,T,iel)
! routine for several gauss point
! provides the tangent dH=mui*dB
implicit none

double precision                 :: mui(:,:,:), si(:,:,:), T(:)
integer                          :: iel

integer                          :: ngp, i

  ngp=size(mui,3)
  mui=0d0
  si=0d0
  do i=1,ngp
    call magmater_tangent_getH2(mui(:,:,i),si(:,:,i),T(i),iel,i)
  end do
 

end subroutine magmater_tangent_getH1

 
subroutine magmater_tangent_getH2(mui,si,temper,iel,igp)
! routine for one gauss point
! provides the tangent dH=mui*dB
implicit none

double precision, intent(out)    :: mui(:,:), si(:,:)
double precision, intent(in)     :: temper
integer                          :: iel
integer, optional                :: igp

double precision                 :: H(3)
double precision                 :: A, mu0r, Bs, C, Tc
double precision                 :: Bee, Hee, dfdH, A1, A2, fH, delta, res, mu
double precision                 :: Id(3,3), mu3(3,3), mui3(3,3), tmp1, dBdHee, dmudHk(3,3)
integer                          :: gp, mat, ndim
  
! special for 3-node element with one gauss point, if call without igp
  if(present(igp))then
    gp=igp
  else
    gp=1  
  endif

! specific material
  mat=matId(iel)
  mu0r=propEM(mat,2)
  Bs=propEM(mat,3)
  C=propEM(mat,4)
  Tc=propEM(mat,5)
       
  H=Hsave(:,gp,iel,new)
  Hee=dsqrt(H(1)**2+H(2)**2+H(3)**2)  
   
  Id(1,:)=[1d0, 0d0, 0d0]
  Id(2,:)=[0d0, 1d0, 0d0]
  Id(3,:)=[0d0, 0d0, 1d0]

  A = pi*mu0*(mu0r - 1d0)/(2d0*Bs)
  A1 = 2d0*Bs/pi*datan(A*Hee)
  A2 = 1d0-dexp((temper-Tc)/C)
!  write(*,*)'tan Hee',Hee, temper, Tc,C
  if (Hee.ne.0d0) then
  
    Bee = mu0*Hee + A1*A2
    dBdHee = mu0 + 2d0*Bs/pi*(A/((A*Hee)**2 + 1d0))*(1d0-exp((temper-Tc)/C))

    tmp1=(dBdHee-Bee/Hee)/Hee**2
   
    dmudHk(1,:)=(/H(1)*H(1),H(1)*H(2),H(1)*H(3)/)
    dmudHk(2,:)=(/H(2)*H(1),H(2)*H(2),H(2)*H(3)/)
    dmudHk(3,:)=(/H(3)*H(1),H(3)*H(2),H(3)*H(3)/)
    dmudHk=dmudHk*tmp1
   
    mu=Bee/Hee;
    
    mu3=dmudHk+mu*Id
  else
 
    mu=mu0+2d0*Bs/pi*A*A2
    mu3=mu*Id
  
  end if
  
  call inv3(mui3,mu3)
  ndim=size(mui,2)
  if (ndim.eq.3) then
    mui(:,:)=mui3
    si=0d0
    si(1,1)=propEm(matId(iel),1)
    si(2,2)=propEm(matId(iel),2)
    si(3,3)=propEm(matId(iel),3)
  else
    mui(1:2,1:2)=mui3(1:2,1:2)
    si(1,1)=propEm(matId(iel),1)
    si(2,2)=propEm(matId(iel),2)
  end if 

end subroutine magmater_tangent_getH2

!subroutine magmater_calcH0(H,B,iel)

subroutine magmater_calcH1(H,B,T,iel)
implicit none

double precision                  :: H(:,:), B(:,:), T(:)
integer                           :: iel

integer                           :: igp, ngp

  ngp=size(B,2)
  do igp=1,ngp
    call magmater_calcH2(H(:,igp),B(:,igp),T(igp),iel,igp)
  end do
  
end subroutine magmater_calcH1


subroutine magmater_calcH2(H,B,temper,iel,igp)
implicit none

double precision                  :: H(:), B(:), temper
integer                           :: iel
integer, optional                 :: igp

double precision                  :: A, mu0r, Bs, C, Tc, Bsaves, DBee, Bold
double precision                  :: Bee, Hee, dfdH, A1, A2, fH, delta, res, mu
integer                           :: mat, j, i, steps, conv
integer, parameter                :: maxiter=8

    mat=matId(iel)
    mu0r=propEM(mat,2)
    Bs=propEM(mat,3)
    C=propEM(mat,4)
    Tc=propEM(mat,5)
  
    A = pi*mu0*(mu0r - 1d0)/(2d0*Bs)

    if (size(B).eq.2) then
       Bee=dsqrt(B(1)**2+B(2)**2)
    else
       Bee=dsqrt(B(1)**2+B(2)**2+B(3)**2)
    end if   
!
! Note that Bold is not Bold in the previous increment if 
! temperature changed. Bold is a quantity such that we start 
! from a state where the residual is zero, for this purpose 
! Bold is calculated from the residual equation
!
    A2 = 1d0-dexp((temper-Tc)/C)
    IF (Bee.ne.0d0) THEN
    
      Bsaves=Bee
      H=Hsave(:,igp,iel,old)
      Hee=dsqrt(H(1)**2+H(2)**2+H(3)**2)  ! start value
      A1 = 2d0*Bs/pi*datan(A*Hee)
      Bold=mu0*Hee + A1*A2        


! If no convergence we split the increment into steps
! until it will converge  
! tested, plot from simulation shows that we obtain Bee=f(Hee) curve
!       
      steps=1
      conv=0
      do while (conv.ne.1)
        Hee=dsqrt(H(1)**2+H(2)**2+H(3)**2)  ! start value
        DBee=(Bsaves-Bold)/dble(steps)
        Bee=Bold

        do i=1,steps
          Bee=Bee+DBee
    
          iters: DO j = 1,maxiter
            A1 = 2d0*Bs/pi*datan(A*Hee)
            fH = Bee - mu0*Hee - A1*A2
            dfdH = - mu0 - 2d0*Bs/pi*(A/((A*Hee)**2 + 1d0))*A2
            delta = fH/dfdH
            Hee = Hee - delta
            res = fH
 !           write(*,*)'iter'
            IF (abs(delta).lt.1.d-1.and.res.lt.1d-2) EXIT        
          ENDDO iters

          IF (j.gt.maxiter) then
            steps=steps+1
            conv=0
            exit
          else
            conv=1
          end if
                    
        end do
              
      end do

      mu=Bee/Hee
      H=B/mu
       
     else
       H=0d0  
     end if
      
   if (size(B).eq.3) then
     Hsave(:,igp,iel,new)=H
     Bsave(:,igp,iel,new)=B
   else
     ! 2D problem
     Hsave(1,igp,iel,new)=H(1)
     Hsave(2,igp,iel,new)=H(2)
     Hsave(3,igp,iel,new)=0d0
     Bsave(1,igp,iel,new)=B(1)
     Bsave(2,igp,iel,new)=B(2)
     Bsave(3,igp,iel,new)=0d0
   end if
     ! noitert 

end subroutine magmater_calcH2


subroutine magmater_getnonl1H(field,Hout,iel)
implicit none

character                         :: field*(*)
double precision                  :: Hout(:,:)
integer                           :: iel


  if (trim(field).eq.'Hfield') then
    if (size(Hout,1).eq.3) then
      Hout=Hsave(:,:,iel,new)
    else
      Hout=Hsave(1:2,:,iel,new)
    end if
  elseif (trim(field).eq.'Bfield') then
    if (size(Hout,1).eq.3) then
      Hout=Bsave(:,:,iel,new)
    else
      Hout=Bsave(1:2,:,iel,new)
    end if
  elseif (trim(field).eq.'Bintr') then
    Hout=0d0
  else
    write(*,*)'Field is note stored in mater_magnnonl1'
    write(*,*)field
    stop
  end if

end subroutine magmater_getnonl1H




end module  mater_magnnonl1
