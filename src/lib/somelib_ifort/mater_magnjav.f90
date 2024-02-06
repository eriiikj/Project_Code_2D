module mater_magnjav

! Last modified
! M. Ristinmaa 2021-02-13
!  initial version based on ode23
! M. Ristinmaa 2021-02-14
!  tangent is dB=mu*dH is implemented
! M. Ristinmaa 2021-02-14
!  changed order new,iel,gp to gp,iel,new
! M. Ristinmaa 2021-02-24
!  new call getH
! M. Ristinmaa 2021-03-17
!  N-R solution implemented
!  ode23 is used as failsafe if N-R solution fails
!  added new Bintr output to magmater_getjavH
!
!-----------------------------------------------------------------------------
!
! contains the Jiles-Atherton hysteresis model
!  D.C. Jiles and D.L. Atherton, L. Magn. Magn. Mater, vol 61, 48-60, 1986.
! based on the vector format by
!  A. Bergqvist IEEE transaction on magnetics, vol 32, 4213-4215, 1996
! and inverse format by
!  J.V. Leite et al IEEE transactions on magnetics, vol 10, 1769-1775, 2004

use matrix_util, only: inv3 
use fem_system, only: solveq

implicit none

interface magmater_calcjav
   module procedure magmater_javcalcH1
!   module procedure magmater_javcalcH2  ! ode23 replaced by magmater_javcalcNR
   module procedure magmater_javcalcNR
end interface
private magmater_javcalcH1, magmater_javcalcH2, magmater_javcalcNR
private mater_calcrev, mater_calcirr, mater_calcreveuler, mater_javnrstep

interface magmater_tangentjav
   module procedure magmater_javmater_getH1
   module procedure magmater_javmater_getH2
end interface
private magmater_javmater_getH1, magmater_javmater_getH2


private magmater_javtangent

  double precision, allocatable     :: propEM(:,:)
  double precision, allocatable     :: Bsave(:,:,:,:), Hsave(:,:,:,:)

  integer, allocatable              :: matId(:)
  integer                           :: New, Old

  double precision, parameter       :: mu0=1.256637061435917d-06

  private propEM, matID, mu0, Bsave, Hsave
  private New, Old

!-----------------------------------------------------------------------------
contains 


subroutine magmater_initjav(propEMin,matIdin,nip)
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
  allocate(propEM(nmat,18), stat=ierr)

  propEM=0d0
  if (nprp.eq.6) then 
! isotropy
    propEM(:,1)=propEmin(:,1)    !sx
    propEM(:,2)=propEmin(:,1)    !sy
    propEM(:,3)=propEmin(:,1)    !sz
    propEM(:,4)=propEmin(:,2)    !cx
    propEM(:,5)=propEmin(:,2)    !cy
    propEM(:,6)=propEmin(:,2)    !cz
    propEM(:,7)=propEmin(:,3)    !ax
    propEM(:,8)=propEmin(:,3)    !ay
    propEM(:,9)=propEmin(:,3)    !az
    propEM(:,10)=propEmin(:,4)    !alphax
    propEM(:,11)=propEmin(:,4)    !alphay
    propEM(:,12)=propEmin(:,4)    !alphaz
    propEM(:,13)=propEmin(:,5)   !Msx
    propEM(:,14)=propEmin(:,5)   !Msy
    propEM(:,15)=propEmin(:,5)   !Msz
    propEM(:,16)=propEmin(:,6)   !kx
    propEM(:,17)=propEmin(:,6)   !ky
    propEM(:,18)=propEmin(:,6)   !kz
  else
    propEM=propEMin
  end if

  matId=matIdin

  allocate(Bsave(3,nip,nels,2), stat=ierr)
  IF (ierr /= 0) STOP "*** Not enough memory ***"
  allocate(Hsave(3,nip,nels,2), stat=ierr)
  IF (ierr /= 0) STOP "*** Not enough memory ***"
!  allocate(Bsave(2,nels,nip,3), stat=ierr)
!  allocate(Hsave(2,nels,nip,3), stat=ierr)

  Bsave=0d0
  Hsave=0d0
  New=2
  Old=1

end subroutine magmater_initjav


subroutine magmater_acceptjav
implicit none

!  if (Old.eq.1) then
!    New=1 
!    Old=2
!  else
!    New=2 
!    Old=1
!  end if

! This is more correct 
! since we calculate the tangent for New
! in first iteration using the above 
! it will be the previous state
!

  Bsave(:,:,:,Old)=Bsave(:,:,:,New)
  Hsave(:,:,:,Old)=Hsave(:,:,:,New)

end subroutine magmater_acceptjav


subroutine magmater_javmater_getH1(mui,si,iel)
! routine for several gauss point
! provides the tangent dH=mui*dB
implicit none

double precision                 :: mui(:,:,:), si(:,:,:)
integer                          :: iel

integer                          :: ngp, i

  ngp=size(mui,3)
  mui=0d0
  si=0d0
  do i=1,ngp
!    call magmater_javmater_getH2(mui(:,:,i),si(:,:,i),iel,i)  !ode23
    call magmater_javmater_tangentNR(mui(:,:,i),si(:,:,i),iel,i)
  end do
 

end subroutine magmater_javmater_getH1

subroutine magmater_javcalcH1(H,B,iel)
implicit none

double precision                  :: H(:,:), B(:,:)
integer                           :: iel

integer                           :: igp, ngp
  
  ngp=size(B,2)
  do igp=1,ngp
!    call magmater_javcalcH2(H(:,igp),B(:,igp),iel,igp)  ! ode23
    call magmater_javcalcNR(H(:,igp),B(:,igp),iel,igp)
  end do
  
  
end subroutine magmater_javcalcH1
 

! ----- N-R solution

subroutine magmater_javcalcNR(H,B,iel,igp)
implicit none

double precision                  :: H(:), B(:)
integer                           :: iel, igp

double precision                  :: Hs(3), Hold(3), Hstart(3)
double precision                  :: Bstart(3), Bend(3), Bold(3), dB(3), dBstep(3)
double precision                  :: mpropEM(18)
integer                           :: conv, steps, i, flag, try_ode
  

! first test one step
  conv=0
  try_ode=0
  steps=1
  Bold=Bsave(:,igp,iel,old)
  if (size(B).eq.3) then
    dB=B-Bold
  else
    dB(1)=B(1)-Bold(1)
    dB(2)=B(2)-Bold(2)
    dB(3)=0d0
  end if
  
  mpropEM=propEm(matId(iel),:)

  
  do while (conv.eq.0.and.try_ode.eq.0)
  
    dBstep=dB/dble(steps)
    Hs=Hsave(:,igp,iel,old)
    conv=1
  
    i=1
    do while (i.le.steps.and.conv.eq.1) 
      
      Hstart=Hs
      Bstart=Bold+dBstep*dble(i-1)
      Bend=Bstart+dBstep
            
      call mater_javnrstep(Hs,flag,Bend,Bstart,Hstart,mpropEM)
     
      if (flag.eq.1) then
          conv=0
          !write(*,*)'Diverge'
          !pause
      end if
      i=i+1
      
    end do
  
    if (conv.eq.0) then
      if (steps.gt.500) then
        !write(*,*)'Problem convergence in mater_jav'
        !write(*,*)'Element and gauss point ',iel,igp
        conv=0
        try_ode=1
        !stop
        ! We should call ode23 as a backup solution
      end if
      steps=steps+4;
    end if
  
  end do

  if (try_ode.eq.1) then
    !write(*,*)'Newton-Raphson failed, failsafe switching to ode23'
    call magmater_javcalcH2(Hs,B,iel,igp)  ! ode23
 !   write(*,*)'H ',Hs
 !   stop
  end if

  if (size(B).eq.3) then
    H=Hs
    Hsave(:,igp,iel,new)=Hs
    Bsave(:,igp,iel,new)=B
  else
    ! 2D problem
    H(1)=Hs(1)
    H(2)=Hs(2)
    Hsave(1,igp,iel,new)=Hs(1)
    Hsave(2,igp,iel,new)=Hs(2)
    Hsave(3,igp,iel,new)=0d0
    Bsave(1,igp,iel,new)=B(1)
    Bsave(2,igp,iel,new)=B(2)
    Bsave(3,igp,iel,new)=0d0
  end if

end subroutine magmater_javcalcNR


subroutine mater_javnrstep(H,flag,B,Bold,Hold,mpropEm)
implicit none

double precision, intent(out)     :: H(:)
double precision, intent(in)      :: B(:), Bold(:), Hold(:), mpropEm(:)
integer                           :: flag

double precision                  :: val, dlambda, Mirr(3)


  call mater_calcreveuler(val,B,Bold,Hold,mpropEm)
 
  if (val.ge.0d0) then
    call mater_calcirr(H,Mirr,dlambda,flag,B,Bold,Hold,mpropEm)
    if (dlambda.lt.0.and.flag.eq.0) then
      call mater_calcrev(H,dlambda,flag,B,Bold,Hold,mpropEm)
    end if
  end if

  if (val.lt.0d0) then
    call mater_calcrev(H,dlambda,flag,B,Bold,Hold,mpropEm)
!    write(*,*)'Dlambda ',dlambda,flag
    if (dlambda.gt.0d0.and.flag.eq.0) then
      call mater_calcirr(H,Mirr,dlambda,flag,B,Bold,Hold,mpropEm) 
    end if
  end if

end subroutine mater_javnrstep

subroutine mater_calcreveuler(val,B,Bold,Hold,mpropEm)
!
! Predictor to find if rev or irr
! Euler step without irr-part
! 
implicit none

double precision, intent(out)     :: val
double precision, intent(in)      :: B(:), Bold(:), Hold(:), mpropEm(:)

double precision                  :: cx, cy, cz, cM(3,3)
double precision                  :: ax, ay, az, alphax, alphay, alphaz, alphaM(3,3)
double precision                  :: Msx, Msy, Msz, kx, ky, kz, kM(3,3), Id(3,3)

double precision                  :: dB(3), Mold(3), He(3), normHe, dH(3), Hguess(3)

double precision                  :: Manx, Many, Manz, dMan(3), dM(3), dHe(3)
double precision                  :: xsix, xsixy, xsixz, xsiy, xsiyx, xsiyz
double precision                  :: xsiz, xsizx, xsizy, xsiM(3,3)
double precision                  :: cothx, cothy, cothz
double precision                  :: chix, chiy, chiz, chiv(3), normchi
double precision                  :: tmpinv(3,3)

! parameters
  cx=mpropEM(4)
  cy=mpropEM(5)
  cz=mpropEM(6)
  cM(1,:)=[cx, 0d0, 0d0]
  cM(2,:)=[0d0, cy, 0d0]
  cM(3,:)=[0d0, 0d0, cz]
 

  ax=mpropEM(7)
  ay=mpropEM(8)
  az=mpropEM(9)

  alphax=mpropEM(10)
  alphay=mpropEM(11)
  alphaz=mpropEM(12)
  alphaM(1,:)=[alphax, 0d0, 0d0]
  alphaM(2,:)=[0d0, alphay, 0d0]
  alphaM(3,:)=[0d0, 0d0, alphaz]

  Msx=mpropEM(13)
  Msy=mpropEM(14)
  Msz=mpropEM(15)

  kx=mpropEM(16)
  ky=mpropEM(17)
  kz=mpropEM(18)
  kM(1,:)=[kx, 0d0, 0d0]
  kM(2,:)=[0d0, ky, 0d0]
  kM(3,:)=[0d0, 0d0, kz]

  Id(1,:)=[1d0, 0d0, 0d0]
  Id(2,:)=[0d0, 1d0, 0d0]
  Id(3,:)=[0d0, 0d0, 1d0]
  
! calculation
  dB=B-Bold
  Mold=B/mu0-Hold
  He=Hold+matmul(alphaM,Mold)
  normHe=dsqrt(He(1)**2+He(2)**2+He(3)**2)

  if (normHe.lt.1d-6) normHe=1d-6
        
  cothx=1d0/dtanh(normHe/ax)
  cothy=1d0/dtanh(normHe/ay)
  cothz=1d0/dtanh(normHe/az)
  Manx=Msx*(cothx-ax/normHe)*He(1)/normHe
  Many=Msy*(cothy-ay/normHe)*He(2)/normHe
  Manz=Msz*(cothz-az/normHe)*He(3)/normHe
  
  xsix=Msx/ax*(1d0-cothx**2+(ax/normHe)**2)*He(1)**2/normHe**2 &
     + Msx*(cothx-ax/normHe)*(1d0/normHe-He(1)**2/normHe**3)
  xsixy=Msx/ax*(1-cothx**2+(ax/normHe)**2)*He(1)*He(2)/normHe**2 &
     - Msx*(cothx-ax/normHe)*(He(1)*He(2)/normHe**3)
  xsixz=Msx/ax*(1-cothx**2+(ax/normHe)**2)*He(1)*He(3)/normHe**2 &
     - Msx*(cothx-ax/normHe)*(He(1)*He(3)/normHe**3)

  xsiy=Msy/ay*(1d0-cothy**2+(ay/normHe)**2)*He(2)**2/normHe**2 &
     + Msy*(cothy-ay/normHe)*(1d0/normHe-He(2)**2/normHe**3)
  xsiyx=Msy/ay*(1-cothy**2+(ay/normHe)**2)*He(2)*He(1)/normHe**2 &
     - Msy*(cothy-ay/normHe)*(He(2)*He(1)/normHe**3)
  xsiyz=Msy/ay*(1-cothy**2+(ay/normHe)**2)*He(2)*He(3)/normHe**2 &
     - Msy*(cothy-ay/normHe)*(He(2)*He(3)/normHe**3)

  xsiz=Msz/az*(1d0-cothz**2+(az/normHe)**2)*He(3)**2/normHe**2 &
     + Msz*(cothz-az/normHe)*(1d0/normHe-He(3)**2/normHe**3)
  xsizx=Msz/az*(1-cothy**2+(az/normHe)**2)*He(3)*He(1)/normHe**2 &
     - Msz*(cothy-az/normHe)*(He(3)*He(1)/normHe**3)
  xsizy=Msz/az*(1-cothy**2+(az/normHe)**2)*He(3)*He(2)/normHe**2 &
     - Msz*(cothy-az/normHe)*(He(3)*He(2)/normHe**3)

  xsiM(1,:)=[ xsix, xsixy, xsixz]
  xsiM(2,:)=[xsiyx,  xsiy, xsiyz]
  xsiM(3,:)=[xsizx, xsiyz,  xsiz]
  
  chix=1/kx*(Manx-Mold(1))
  chiy=1/ky*(Many-Mold(2))
  chiz=1/kz*(Manz-Mold(3))

  normchi=dsqrt(chix**2+chiy**2+chiz**2)

  chiv=[chix, chiy, chiz]

  Hguess=B/mu0-Mold
  dH=Hguess-Hold
  dMan=matmul(xsiM,dH)

  call inv3(tmpinv,Id+cM*xsiM*(Id-alphaM))
  dM=1d0/mu0*matmul(tmpinv,matmul((cM*xsiM),dB))
  dH=dB/mu0-dM
  dHe=dH+matmul(alphaM,dM)

!------------------------------------------
!
! 99% correct 
! can be used as a predictor to
! select rev or irr calculation
!
  val=dot_product(chiv,dHe+matmul(alphaM,matmul(cM,dMan)))

end subroutine mater_calcreveuler


subroutine mater_calcirr(H,Mirr,dlambda,flag,B,Bold,Hold,mpropEM)
! irreversible part is present
!
implicit none

double precision, intent(out)     :: H(3), Mirr(3), dlambda
integer, intent(out)              :: flag
double precision, intent(in)      :: B(3), Bold(:), Hold(:), mpropEm(:)


double precision                  :: cx, cy, cz, cM(3,3)
double precision                  :: ax, ay, az, alphax, alphay, alphaz, alphaM(3,3)
double precision                  :: Msx, Msy, Msz, kx, ky, kz, kM(3,3), Id(3,3), kin(3,3)

double precision                  :: Mold(3), He(3), Heold(3), DHe(3), normHe
double precision                  :: Mirrold(3), DMirr(3), M(3)

double precision                  :: Manx, Many, Manz, Man(3), Mann(3)
double precision                  :: xsix, xsixy, xsixz, xsiy, xsiyx, xsiyz
double precision                  :: xsiz, xsizx, xsizy, xsiM(3,3)
double precision                  :: cothx, cothy, cothz
double precision                  :: chi(3), normchi

double precision                  :: Res1(3), Res2(3), Res(6)
double precision                  :: tmpinv(3,3), tmp1(3,3), tmv1(3), beta(3,3)
double precision                  :: K(6,6)

integer                           :: iter, nriter
! parameters
  cx=mpropEM(4)
  cy=mpropEM(5)
  cz=mpropEM(6)
  cM(1,:)=[cx, 0d0, 0d0]
  cM(2,:)=[0d0, cy, 0d0]
  cM(3,:)=[0d0, 0d0, cz]
 

  ax=mpropEM(7)
  ay=mpropEM(8)
  az=mpropEM(9)

  alphax=mpropEM(10)
  alphay=mpropEM(11)
  alphaz=mpropEM(12)
  alphaM(1,:)=[alphax, 0d0, 0d0]
  alphaM(2,:)=[0d0, alphay, 0d0]
  alphaM(3,:)=[0d0, 0d0, alphaz]

  Msx=mpropEM(13)
  Msy=mpropEM(14)
  Msz=mpropEM(15)

  kx=mpropEM(16)
  ky=mpropEM(17)
  kz=mpropEM(18)
  kM(1,:)=[kx, 0d0, 0d0]
  kM(2,:)=[0d0, ky, 0d0]
  kM(3,:)=[0d0, 0d0, kz]
  kin(1,:)=[1d0/kx,0d0,0d0]
  kin(2,:)=[0d0,1d0/ky,0d0]
  kin(3,:)=[0d0,0d0,1d0/kz]
  

  Id(1,:)=[1d0, 0d0, 0d0]
  Id(2,:)=[0d0, 1d0, 0d0]
  Id(3,:)=[0d0, 0d0, 1d0]



! initiate
! to improve convergence
! Starting guess from last iterative state (not implemented)
    
  Mold=Bold/mu0-Hold
  He=Hold+matmul(alphaM,Mold)
  Heold=He 
  normHe=dsqrt(He(1)**2+He(2)**2+He(3)**2)

  if (normHe.le.1d-6) normHe=1d-6
 
  cothx=1d0/dtanh(normHe/ax)
  cothy=1d0/dtanh(normHe/ay)
  cothz=1d0/dtanh(normHe/az)
  Manx=Msx*(cothx-ax/normHe)*He(1)/normHe
  Many=Msy*(cothy-ay/normHe)*He(2)/normHe
  Manz=Msz*(cothz-az/normHe)*He(3)/normHe
  Mann=[Manx, Many, Manz]
  Man=Mann

  call inv3(tmpinv,Id-cM)
  Mirrold=matmul(tmpinv,Mold-matmul(cM,Mann))
  Mirr=Mirrold
  DMirr=0d0
  DHe=0d0
  
!  DB=B-Bold

  Dlambda=0d0

  chi=matmul(kin,Man-Mirr)
  normchi=dsqrt(chi(1)**2+chi(2)**2+chi(3)**2)
  if (normchi<1d-6) normchi=1d-6
   
  Res1=He+matmul(Id-alphaM,matmul(Id-cM,Mirr))+matmul((Id-alphaM),matmul(cM,Man))-B/mu0
  Res2=DMirr-Dlambda*chi/normchi;
  
  Res(1:3)=Res1
  Res(4:6)=Res2
  
  
! iterations

  iter=0
  nriter=6
  do while (norm2(Res).ge.1d-5.and.iter.lt.nriter)
   
    iter=iter+1
    
    xsix=Msx/ax*(1d0-cothx**2+(ax/normHe)**2)*He(1)**2/normHe**2 &
       + Msx*(cothx-ax/normHe)*(1d0/normHe-He(1)**2/normHe**3)
    xsixy=Msx/ax*(1-cothx**2+(ax/normHe)**2)*He(1)*He(2)/normHe**2 &
       - Msx*(cothx-ax/normHe)*(He(1)*He(2)/normHe**3)
    xsixz=Msx/ax*(1-cothx**2+(ax/normHe)**2)*He(1)*He(3)/normHe**2 &
       - Msx*(cothx-ax/normHe)*(He(1)*He(3)/normHe**3)

    xsiy=Msy/ay*(1d0-cothy**2+(ay/normHe)**2)*He(2)**2/normHe**2 &
       + Msy*(cothy-ay/normHe)*(1d0/normHe-He(2)**2/normHe**3)
    xsiyx=Msy/ay*(1-cothy**2+(ay/normHe)**2)*He(2)*He(1)/normHe**2 &
       - Msy*(cothy-ay/normHe)*(He(2)*He(1)/normHe**3)
    xsiyz=Msy/ay*(1-cothy**2+(ay/normHe)**2)*He(2)*He(3)/normHe**2 &
       - Msy*(cothy-ay/normHe)*(He(2)*He(3)/normHe**3)

    xsiz=Msz/az*(1d0-cothz**2+(az/normHe)**2)*He(3)**2/normHe**2 &
       + Msz*(cothz-az/normHe)*(1d0/normHe-He(3)**2/normHe**3)
    xsizx=Msz/az*(1-cothy**2+(az/normHe)**2)*He(3)*He(1)/normHe**2 &
       - Msz*(cothy-az/normHe)*(He(3)*He(1)/normHe**3)
    xsizy=Msz/az*(1-cothy**2+(az/normHe)**2)*He(3)*He(2)/normHe**2 &
       - Msz*(cothy-az/normHe)*(He(3)*He(2)/normHe**3)

    xsiM(1,:)=[ xsix, xsixy, xsixz]
    xsiM(2,:)=[xsiyx,  xsiy, xsiyz]
    xsiM(3,:)=[xsizx, xsiyz,  xsiz]
         
    tmp1=matmul(reshape(chi,(/3,1/)),reshape(chi,(/1,3/)))
    beta=Dlambda/normchi*Id-Dlambda*tmp1/normchi**3
     
    ! k11 
    k(1:3,1:3)=Id+matmul(Id-alphaM,matmul(cM,xsiM))
    ! k12
    k(1:3,4:6)=matmul(Id-alphaM,Id-cM)
   
   
    tmv1=chi+matmul(DHe,matmul(kin,xsiM))
    tmp1=matmul(reshape(chi,(/3,1/)),reshape(tmv1,(/1,3/)))
    ! k21
    k(4:6,1:3)=-tmp1/normchi-matmul(beta,matmul(kin,xsiM))
    
    tmv1=matmul(DHe,kin)
    tmp1=matmul(reshape(chi,(/3,1/)),reshape(tmv1,(/1,3/)))
    ! k22
    k(4:6,4:6)=Id+matmul(beta,kin)+tmp1/normchi

    call solveq(K,Res)
       
    He=He-Res(1:3);
    DHe=He-Heold;
    Mirr=Mirr-Res(4:6);
    DMirr=Mirr-Mirrold;
   
    normHe=dsqrt(He(1)**2+He(2)**2+He(3)**2)

    cothx=1d0/dtanh(normHe/ax)
    cothy=1d0/dtanh(normHe/ay)
    cothz=1d0/dtanh(normHe/az)
    Manx=Msx*(cothx-ax/normHe)*He(1)/normHe
    Many=Msy*(cothy-ay/normHe)*He(2)/normHe
    Manz=Msz*(cothz-az/normHe)*He(3)/normHe
    Man=[Manx, Many, Manz]
   
    chi=matmul(kin,Man-Mirr)
    normchi=sqrt(chi(1)**2+chi(2)**2+chi(3)**2)
    Dlambda=dot_product(chi,DHe)
   
    Res1=He+matmul(Id-alphaM,matmul(Id-cM,Mirr))+matmul((Id-alphaM),matmul(cM,Man))-B/mu0
    Res2=DMirr-Dlambda*chi/normchi;
  
    Res(1:3)=Res1
    Res(4:6)=Res2
        
         
  end do
  
  M=matmul(Id-cM,Mirr)+matmul(cM,Man)
  H=B/mu0-M
    
  if (iter.eq.nriter) then
      flag=1
  else
      flag=0
  end if
  
end subroutine mater_calcirr


subroutine mater_calcrev(H,dlambda,flag,B,Bold,Hold,mpropEm)
! irreversible part is present
!
implicit none

double precision, intent(out)     :: H(3), dlambda
integer, intent(out)              :: flag
double precision, intent(in)      :: B(:), Bold(:), Hold(:), mpropEm(:)


double precision                  :: cx, cy, cz, cM(3,3)
double precision                  :: ax, ay, az, alphax, alphay, alphaz, alphaM(3,3)
double precision                  :: Msx, Msy, Msz, kx, ky, kz, kM(3,3), Id(3,3), kin(3,3)

double precision                  :: Mold(3), He(3), Heold(3), DHe(3), normHe
double precision                  :: Mirrold(3), DMirr(3), M(3)

double precision                  :: Manx, Many, Manz, Man(3), Mann(3)
double precision                  :: xsix, xsixy, xsixz, xsiy, xsiyx, xsiyz
double precision                  :: xsiz, xsizx, xsizy, xsiM(3,3)
double precision                  :: cothx, cothy, cothz
double precision                  :: chi(3), normchi

double precision                  :: Res1(3), Res2(3), Res(3)
double precision                  :: tmpinv(3,3), tmp1(3,3), tmv1(3), beta(3,3)
double precision                  :: K(3,3)
integer                           :: iter, nriter

! parameters
  cx=mpropEM(4)
  cy=mpropEM(5)
  cz=mpropEM(6)
  cM(1,:)=[cx, 0d0, 0d0]
  cM(2,:)=[0d0, cy, 0d0]
  cM(3,:)=[0d0, 0d0, cz]
 

  ax=mpropEM(7)
  ay=mpropEM(8)
  az=mpropEM(9)

  alphax=mpropEM(10)
  alphay=mpropEM(11)
  alphaz=mpropEM(12)
  alphaM(1,:)=[alphax, 0d0, 0d0]
  alphaM(2,:)=[0d0, alphay, 0d0]
  alphaM(3,:)=[0d0, 0d0, alphaz]

  Msx=mpropEM(13)
  Msy=mpropEM(14)
  Msz=mpropEM(15)

  kx=mpropEM(16)
  ky=mpropEM(17)
  kz=mpropEM(18)
  kM(1,:)=[kx, 0d0, 0d0]
  kM(2,:)=[0d0, ky, 0d0]
  kM(3,:)=[0d0, 0d0, kz]
  kin(1,:)=[1d0/kx,0d0,0d0]
  kin(2,:)=[0d0,1d0/ky,0d0]
  kin(3,:)=[0d0,0d0,1d0/kz]
  

  Id(1,:)=[1d0, 0d0, 0d0]
  Id(2,:)=[0d0, 1d0, 0d0]
  Id(3,:)=[0d0, 0d0, 1d0]



  Mold=Bold/mu0-Hold
  He=Hold+matmul(alphaM,Mold)
  Heold=He

  normHe=dsqrt(He(1)**2+He(2)**2+He(3)**2)

  if (normHe.le.1d-6) normHe=1d-6
     
  cothx=1d0/dtanh(normHe/ax)
  cothy=1d0/dtanh(normHe/ay)
  cothz=1d0/dtanh(normHe/az)
  Manx=Msx*(cothx-ax/normHe)*He(1)/normHe
  Many=Msy*(cothy-ay/normHe)*He(2)/normHe
  Manz=Msz*(cothz-az/normHe)*He(3)/normHe
  Mann=[Manx, Many, Manz]
  Man=Mann

  call inv3(tmpinv,Id-cM)
  Mirrold=matmul(tmpinv,Mold-matmul(cM,Mann))

  DHe=0d0

! initial residual
  Res=He+matmul(Id-alphaM,matmul(Id-cM,Mirrold))+matmul(Id-alphaM,matmul(cM,Man))-B/mu0
  
  iter=0
  nriter=5
  do while (norm2(Res).ge.1d-5.and.iter.lt.nriter)
    
    iter=iter+1
   
    xsix=Msx/ax*(1d0-cothx**2+(ax/normHe)**2)*He(1)**2/normHe**2 &
       + Msx*(cothx-ax/normHe)*(1d0/normHe-He(1)**2/normHe**3)
    xsixy=Msx/ax*(1-cothx**2+(ax/normHe)**2)*He(1)*He(2)/normHe**2 &
       - Msx*(cothx-ax/normHe)*(He(1)*He(2)/normHe**3)
    xsixz=Msx/ax*(1-cothx**2+(ax/normHe)**2)*He(1)*He(3)/normHe**2 &
       - Msx*(cothx-ax/normHe)*(He(1)*He(3)/normHe**3)

    xsiy=Msy/ay*(1d0-cothy**2+(ay/normHe)**2)*He(2)**2/normHe**2 &
       + Msy*(cothy-ay/normHe)*(1d0/normHe-He(2)**2/normHe**3)
    xsiyx=Msy/ay*(1-cothy**2+(ay/normHe)**2)*He(2)*He(1)/normHe**2 &
       - Msy*(cothy-ay/normHe)*(He(2)*He(1)/normHe**3)
    xsiyz=Msy/ay*(1-cothy**2+(ay/normHe)**2)*He(2)*He(3)/normHe**2 &
       - Msy*(cothy-ay/normHe)*(He(2)*He(3)/normHe**3)

    xsiz=Msz/az*(1d0-cothz**2+(az/normHe)**2)*He(3)**2/normHe**2 &
       + Msz*(cothz-az/normHe)*(1d0/normHe-He(3)**2/normHe**3)
    xsizx=Msz/az*(1-cothy**2+(az/normHe)**2)*He(3)*He(1)/normHe**2 &
       - Msz*(cothy-az/normHe)*(He(3)*He(1)/normHe**3)
    xsizy=Msz/az*(1-cothy**2+(az/normHe)**2)*He(3)*He(2)/normHe**2 &
       - Msz*(cothy-az/normHe)*(He(3)*He(2)/normHe**3)

    xsiM(1,:)=[ xsix, xsixy, xsixz]
    xsiM(2,:)=[xsiyx,  xsiy, xsiyz]
    xsiM(3,:)=[xsizx, xsiyz,  xsiz]

    K=Id+matmul(Id-alphaM,matmul(cM,xsiM))

    call solveq(K,Res)
   
    He=He-Res
    DHe=DHe-Res

    normHe=dsqrt(He(1)**2+He(2)**2+He(3)**2)
    if (normHe.le.1d-6) normHe=1d-6
     
    cothx=1d0/dtanh(normHe/ax)
    cothy=1d0/dtanh(normHe/ay)
    cothz=1d0/dtanh(normHe/az)
    Manx=Msx*(cothx-ax/normHe)*He(1)/normHe
    Many=Msy*(cothy-ay/normHe)*He(2)/normHe
    Manz=Msz*(cothz-az/normHe)*He(3)/normHe
    Man=[Manx, Many, Manz]
   
    Res=He+matmul(Id-alphaM,matmul(Id-cM,Mirrold))+matmul(Id-alphaM,matmul(cM,Man))-B/mu0
    
  end do

! Check if guess is correct
  chi=matmul(kin,Man-Mirrold)
  dlambda=dot_product(chi,DHe)

  call inv3(tmpinv,Id-alphaM)
  H=Hold+matmul(tmpinv,DHe-1/mu0*matmul(alphaM,B-Bold))
  
  
  if (iter.eq.nriter) then
    flag=1
  else
    flag=0
  end if

end subroutine mater_calcrev


subroutine magmater_javmater_tangentNR(mui,si,iel,igp)
!
! algorithmic tangent for N-R solution
!
implicit none

double precision, intent(out)    :: mui(:,:), si(:,:)
integer                          :: iel
integer, optional                :: igp

double precision                  :: cx, cy, cz, cM(3,3)
double precision                  :: ax, ay, az, alphax, alphay, alphaz, alphaM(3,3)
double precision                  :: Msx, Msy, Msz, kx, ky, kz, kM(3,3), Id(3,3), kin(3,3)

double precision                  :: Manx, Many, Manz, Man(3), Mann(3)
double precision                  :: xsix, xsixy, xsixz, xsiy, xsiyx, xsiyz
double precision                  :: xsiz, xsizx, xsizy, xsiM(3,3)
double precision                  :: cothx, cothy, cothz
double precision                  :: chi(3), normchi

double precision                  :: M(3), B(3), H(3), He(3), Mirr(3), normHe
double precision                  :: Mold(3), Bold(3), Hold(3), Heold(3), DHe(3)
double precision                  :: mpropEm(18), Dlambda
double precision                  :: tmpinv(3,3), tmp1(3,3), tmv1(3), beta(3,3)
double precision                  :: K(6,6), ktr(3,3), k3(3,3)

! parameters

  mpropEM=propEm(matId(iel),:)

  cx=mpropEM(4)
  cy=mpropEM(5)
  cz=mpropEM(6)
  cM(1,:)=[cx, 0d0, 0d0]
  cM(2,:)=[0d0, cy, 0d0]
  cM(3,:)=[0d0, 0d0, cz]
 

  ax=mpropEM(7)
  ay=mpropEM(8)
  az=mpropEM(9)

  alphax=mpropEM(10)
  alphay=mpropEM(11)
  alphaz=mpropEM(12)
  alphaM(1,:)=[alphax, 0d0, 0d0]
  alphaM(2,:)=[0d0, alphay, 0d0]
  alphaM(3,:)=[0d0, 0d0, alphaz]

  Msx=mpropEM(13)
  Msy=mpropEM(14)
  Msz=mpropEM(15)

  kx=mpropEM(16)
  ky=mpropEM(17)
  kz=mpropEM(18)
  kM(1,:)=[kx, 0d0, 0d0]
  kM(2,:)=[0d0, ky, 0d0]
  kM(3,:)=[0d0, 0d0, kz]
  kin(1,:)=[1d0/kx,0d0,0d0]
  kin(2,:)=[0d0,1d0/ky,0d0]
  kin(3,:)=[0d0,0d0,1d0/kz]
  

  Id(1,:)=[1d0, 0d0, 0d0]
  Id(2,:)=[0d0, 1d0, 0d0]
  Id(3,:)=[0d0, 0d0, 1d0]


  Mold=Bsave(:,igp,iel,old)-mu0*Hsave(:,igp,iel,old)  
  Bold=Bsave(:,igp,iel,old)
  Hold=Hsave(:,igp,iel,old)
  M   =Bsave(:,igp,iel,new)-mu0*Hsave(:,igp,iel,new)  
  B   =Bsave(:,igp,iel,new)
  H   =Hsave(:,igp,iel,new)
  
! pre-calculations
  
  M=B/mu0-H
  He=H+matmul(alphaM,M)

  normHe=dsqrt(He(1)**2+He(2)**2+He(3)**2)

  if (normHe.le.1d-6) normHe=1d-6
     
  cothx=1d0/dtanh(normHe/ax)
  cothy=1d0/dtanh(normHe/ay)
  cothz=1d0/dtanh(normHe/az)
  Manx=Msx*(cothx-ax/normHe)*He(1)/normHe
  Many=Msy*(cothy-ay/normHe)*He(2)/normHe
  Manz=Msz*(cothz-az/normHe)*He(3)/normHe
  Man=[Manx, Many, Manz]
   
  xsix=Msx/ax*(1d0-cothx**2+(ax/normHe)**2)*He(1)**2/normHe**2 &
       + Msx*(cothx-ax/normHe)*(1d0/normHe-He(1)**2/normHe**3)
  xsixy=Msx/ax*(1-cothx**2+(ax/normHe)**2)*He(1)*He(2)/normHe**2 &
       - Msx*(cothx-ax/normHe)*(He(1)*He(2)/normHe**3)
  xsixz=Msx/ax*(1-cothx**2+(ax/normHe)**2)*He(1)*He(3)/normHe**2 &
       - Msx*(cothx-ax/normHe)*(He(1)*He(3)/normHe**3)

  xsiy=Msy/ay*(1d0-cothy**2+(ay/normHe)**2)*He(2)**2/normHe**2 &
       + Msy*(cothy-ay/normHe)*(1d0/normHe-He(2)**2/normHe**3)
  xsiyx=Msy/ay*(1-cothy**2+(ay/normHe)**2)*He(2)*He(1)/normHe**2 &
       - Msy*(cothy-ay/normHe)*(He(2)*He(1)/normHe**3)
  xsiyz=Msy/ay*(1-cothy**2+(ay/normHe)**2)*He(2)*He(3)/normHe**2 &
       - Msy*(cothy-ay/normHe)*(He(2)*He(3)/normHe**3)

  xsiz=Msz/az*(1d0-cothz**2+(az/normHe)**2)*He(3)**2/normHe**2 &
       + Msz*(cothz-az/normHe)*(1d0/normHe-He(3)**2/normHe**3)
  xsizx=Msz/az*(1-cothy**2+(az/normHe)**2)*He(3)*He(1)/normHe**2 &
       - Msz*(cothy-az/normHe)*(He(3)*He(1)/normHe**3)
  xsizy=Msz/az*(1-cothy**2+(az/normHe)**2)*He(3)*He(2)/normHe**2 &
       - Msz*(cothy-az/normHe)*(He(3)*He(2)/normHe**3)

  xsiM(1,:)=[ xsix, xsixy, xsixz]
  xsiM(2,:)=[xsiyx,  xsiy, xsiyz]
  xsiM(3,:)=[xsizx, xsiyz,  xsiz]


  if (norm2(M-Mold).lt.1e-1) then  

    ! rev  
    ktr=Id+matmul(Id-alphaM,matmul(cM,xsiM))
    call inv3(tmpinv,matmul(ktr,Id-cM))
    k3=1d0/mu0*matmul(tmpinv,Id-matmul(ktr,alphaM))
    
  else
 
    ! irr
    call inv3(tmpinv,Id-cM)
    Mirr=matmul(tmpinv,M-matmul(cM,Man))
    
    chi=matmul(kin,Man-Mirr)
    normchi=dsqrt(chi(1)**2+chi(2)**2+chi(3)**2)

    if (normchi.le.1d-6) normchi=1d-6

    Mold=Bold/mu0-Hold
    Heold=Hold+matmul(alphaM,Mold)
    DHe=He-Heold

    Dlambda=dot_product(chi,DHe)

    tmp1=matmul(reshape(chi,(/3,1/)),reshape(chi,(/1,3/)))
    beta=Dlambda/normchi*Id-Dlambda*tmp1/normchi**3
     
    ! k11 
    k(1:3,1:3)=Id+matmul(Id-alphaM,matmul(cM,xsiM))
    ! k12
    k(1:3,4:6)=matmul(Id-alphaM,Id-cM)
   
   
    tmv1=chi+matmul(DHe,matmul(kin,xsiM))
    tmp1=matmul(reshape(chi,(/3,1/)),reshape(tmv1,(/1,3/)))
    ! k21
    k(4:6,1:3)=-tmp1/normchi-matmul(beta,matmul(kin,xsiM))
    
    tmv1=matmul(DHe,kin)
    tmp1=matmul(reshape(chi,(/3,1/)),reshape(tmv1,(/1,3/)))
    ! k22
    k(4:6,4:6)=Id+matmul(beta,kin)+tmp1/normchi

    call inv3(tmpinv,K(4:6,4:6))
    ktr=K(1:3,1:3)-matmul(K(1:3,4:6),matmul(tmpinv,K(4:6,1:3)))

    call inv3(tmpinv,matmul(ktr,Id-cM))
    k3=1d0/mu0*matmul(tmpinv,Id-matmul(ktr,alphaM))
   
  end if

  if (size(mui,2).eq.3) then
    mui(:,:)=k3
    si=0d0
    si(1,1)=propEm(matId(iel),1)
    si(2,2)=propEm(matId(iel),2)
    si(3,3)=propEm(matId(iel),3)
  else
    mui(1:2,1:2)=k3(1:2,1:2)
    si(1,1)=propEm(matId(iel),1)
    si(2,2)=propEm(matId(iel),2)
  end if 


end subroutine magmater_javmater_tangentNR

! ----------------------------- ODE solver


subroutine magmater_javmater_getH2(mui,si,iel,igp)
!subroutine magmater_javmater_getH2(iel,igp)
!subroutine magmater_javmater_getH2(iel,igp)
! routine for one gauss point
! provides the tangent dH=mui*dB
implicit none

double precision, intent(out)    :: mui(:,:), si(:,:)
integer                          :: iel
integer, optional                :: igp

double precision                 :: K(3,3), Bcurr(3), Bcont(3), Mcurr(3), Hcurr(3)
double precision                 :: mpropEM(18), Id(3,3), mui3(3,3), dB(3), normdB
integer                          :: gp, ndim
  
! special for 3-node element with one gauss point, if call without igp
  if(present(igp))then
    gp=igp
  else
    gp=1  
  endif
  gp=igp

! specific material
  mpropEM=propEm(matId(iel),:)
    
! tangent should be calculated at current state Bcurr
! for this purpose we assume that the dB will continue in
! the same direction as in the previous increment providing Bcont
! this is only use for the selection revesible or irrevesible 
! We should scale dB such that the norm of it is not too big or small

  Bcurr=Bsave(:,gp,iel,new)
  Hcurr=Hsave(:,gp,iel,new)
  dB=(Bsave(:,gp,iel,new)-Bsave(:,gp,iel,old))
  normdB=dsqrt(dB(1)**2+dB(2)**2+dB(3)**2)
  Bcont=Bcurr+dB/normdB*1d-4
!  Bcont=Bcurr+(Bsave(:,gp,iel,new)-Bsave(:,gp,iel,old))*1d-3
  Mcurr=Bcurr/mu0-Hcurr
  
 ! write(*,*)'calculating tangent'
 
  call magmater_javtangent(K,Bcont,Bcurr,Mcurr,mpropEM)

  Id(1,:)=[1d0, 0d0, 0d0]
  Id(2,:)=[0d0, 1d0, 0d0]
  Id(3,:)=[0d0, 0d0, 1d0]
  
!  dM=K*dB
!  dH=dB/mu0-dM
  mui3=(Id/mu0-K)

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

end subroutine magmater_javmater_getH2




subroutine magmater_javcalcH2(H,B,iel,igp)
! Using ode23 solver for the J-A mode
! local error control is used to calculate 
! if sub-steps are needed.
! ode is  given by dM=f(M)=K*dB
implicit none

double precision                  :: H(:), B(:)
integer                           :: iel
integer, optional                 :: igp

double precision                  :: dB(3), dBstep(3), Bstart(3), Bend(3), Binit(3), K(3,3)
double precision                  :: Mstart(3), Msin(3), Bsin(3), H3(3), Mold(3), M(3)
double precision                  :: s1(3), s2(3), s3(3), s4(3), error(3), mpropEM(18)
double precision                  :: hs
integer                           :: conv, steps, gp, i
  
! special for 3-node element with one gauss point, if call without igp
  if(present(igp))then
    gp=igp
  else
    gp=1
  endif

! specific material
  mpropEM=propEm(matId(iel),:)
  
! We should have a general 3D format  
  if (size(B).eq.3) then
    dB=B-Bsave(:,gp,iel,old)
  else
    ! 2D problem
    dB(1)=B(1)-Bsave(1,gp,iel,old)
    dB(2)=B(2)-Bsave(2,gp,iel,old)
    dB(3)=0d0
  end if
  Mold=Bsave(:,gp,iel,old)/mu0-Hsave(:,gp,iel,old)
  Binit=Bsave(:,gp,iel,old)
  
  steps=1
  conv=0

  do while (conv.eq.0)

    dBstep=dB/dble(steps)
    M=Mold
    hs=1d0
    conv=1

    i=1
    do while (i.le.steps.and.conv.eq.1) 
      
      Mstart=M
      Bstart=Binit+dBstep*dble(i-1)
      Bend=Bstart+dBstep

      ! stage 1

      Msin=Mstart
      Bsin=Bstart
      call magmater_javtangent(K,Bend,Bsin,Msin,mpropEM)
      s1=matmul(K,Bend-Bstart)
      ! stage 2
 
      Msin=Mstart+hs/2d0*s1
      Bsin=Bstart+hs/2d0*dBstep
      call magmater_javtangent(K,Bend,Bsin,Msin,mpropEM)
      s2=matmul(K,Bend-Bstart)
      ! stage 3

      Msin=Mstart+3d0*hs/4d0*s2
      Bsin=Bstart+3d0*hs/4d0*dBstep
      call magmater_javtangent(K,Bend,Bsin,Msin,mpropEM)
      s3=matmul(K,Bend-Bstart)

      ! value

      M=Mstart+hs/9d0*(2d0*s1+3d0*s2+4d0*s3)

      ! error estimate

      call magmater_javtangent(K,Bend,Bend,M,mpropEM)
      s4=matmul(K,Bend-Bstart)

      error=hs/72d0*(-5d0*s1+6d0*s2+8d0*s3-9d0*s4)

      if (dsqrt(dot_product(error,error)).gt.20d0) then
        conv=0
      end if

      i=i+1

    end do

    if (conv.eq.0) then
      steps=steps+3
    end if

  end do

! store values

  if (size(B).eq.3) then
 !   H3=B/mu0-M
    H=B/mu0-M
    Hsave(:,gp,iel,new)=H
    Bsave(:,gp,iel,new)=B
  else
    ! 2D problem
    H(1)=B(1)/mu0-M(1)
    H(2)=B(2)/mu0-M(2)
    Hsave(1,gp,iel,new)=H(1)
    Hsave(2,gp,iel,new)=H(2)
    Hsave(3,gp,iel,new)=0d0
    Bsave(1,gp,iel,new)=B(1)
    Bsave(2,gp,iel,new)=B(2)
    Bsave(3,gp,iel,new)=0d0
  end if
  
end subroutine magmater_javcalcH2



subroutine magmater_javtangent(K,B,Bt,Mold,mpropEm)
! Calculate the tangent of the J-A model
! dM=K*dB where dB is a given steps, in ode-format
! we have that dM=f(M), i.e. input is M (Mold in call)
! Note that Bt and Mold must be calculated at the same 
! time instant, i.e. Bt(t_i) Mold(t_i)
! B the end value is only used for selection between
! reversible and irreversible
!  
implicit none

double precision, intent(in)                  :: B(3), Bt(3), Mold(3), mpropEm(18)
double precision, intent(out)                  :: K(3,3)

double precision                  :: Hold(3), He(3), Hguess(3), dH(3), dMan(3)
double precision                  :: normHe, val, Manx, Many, Manz

double precision                  :: ax, ay, az, Msx, Msy, Msz, xsix, xsiy, xsiz
double precision                  :: chix, chiy, chiz, kx, ky, kz
double precision                  :: cx, cy, cz, alphax, alphay, alphaz
double precision                  :: cM(3,3), alphaM(3,3), kM(3,3), xsiM(3,3), chim(3,1)

double precision                  :: tmp1(3,3), tmp2(3,3), tmpinv(3,3), normchi, Id(3,3)
double precision                  :: tmpv(3)


! define useful matrices and material parameters
  cx=mpropEM(4)
  cy=mpropEM(5)
  cz=mpropEM(6)
  cM(1,:)=[cx, 0d0, 0d0]
  cM(2,:)=[0d0, cy, 0d0]
  cM(3,:)=[0d0, 0d0, cz]
 

  ax=mpropEM(7)
  ay=mpropEM(8)
  az=mpropEM(9)

  alphax=mpropEM(10)
  alphay=mpropEM(11)
  alphaz=mpropEM(12)
  alphaM(1,:)=[alphax, 0d0, 0d0]
  alphaM(2,:)=[0d0, alphay, 0d0]
  alphaM(3,:)=[0d0, 0d0, alphaz]

  Msx=mpropEM(13)
  Msy=mpropEM(14)
  Msz=mpropEM(15)

  kx=mpropEM(16)
  ky=mpropEM(17)
  kz=mpropEM(18)
  kM(1,:)=[kx, 0d0, 0d0]
  kM(2,:)=[0d0, ky, 0d0]
  kM(3,:)=[0d0, 0d0, kz]

  Id(1,:)=[1d0, 0d0, 0d0]
  Id(2,:)=[0d0, 1d0, 0d0]
  Id(3,:)=[0d0, 0d0, 1d0]
  
! calculation
 ! dB=B-B0 
  Hold=Bt/mu0-Mold
  He=Hold+matmul(alphaM,Mold)
  normHe=dsqrt(dot_product(He,He))

  if (dabs(normHe).lt.1d-6) normHe=1d-6

  Manx=Msx*(1d0/dtanh(normHe/ax)-ax/normHe)*He(1)/normHe
  Many=Msy*(1d0/dtanh(normHe/ay)-ay/normHe)*He(2)/normHe
  Manz=Msz*(1d0/dtanh(normHe/az)-az/normHe)*He(3)/normHe

 
  xsix=Msx/ax*(1d0-1d0/dtanh(normHe/ax)**2+(ax/normHe)**2)*He(1)**2/normHe**2 &
     + Msx*(1d0/tanh(normHe/ax)-ax/normHe)*(1d0/normHe-He(1)**2/normHe**3)
  xsiy=Msy/ay*(1d0-1d0/dtanh(normHe/ay)**2+(ay/normHe)**2)*He(2)**2/normHe**2 &
     + Msy*(1d0/tanh(normHe/ay)-ay/normHe)*(1d0/normHe-He(2)**2/normHe**3)
  xsiz=Msz/az*(1d0-1d0/dtanh(normHe/az)**2+(az/normHe)**2)*He(3)**2/normHe**2 &
     + Msz*(1d0/dtanh(normHe/az)-az/normHe)*(1d0/normHe-He(3)**2/normHe**3)

  xsiM(1,:)=[xsix, 0d0, 0d0]
  xsiM(2,:)=[0d0, xsiy, 0d0]
  xsiM(3,:)=[0d0, 0d0, xsiz]

  chix=1/kx*(Manx-Mold(1))
  chiy=1/ky*(Many-Mold(2))
  chiz=1/kz*(Manz-Mold(3))
  chim(1,1)=chix
  chim(2,1)=chiy
  chim(3,1)=chiz
  normchi=dsqrt(dot_product(reshape(chim,(/3/)),reshape(chim,(/3/))))
  if (dabs(normchi).lt.1d-6) normchi=1d-6

! check 
  Hguess=B/mu0-Mold
  dH=Hguess-Hold
  dMan=matmul(xsiM,dH)

  val=dot_product(reshape(chim,(/3/)),dH+matmul(alphaM,matmul(cM,dMan)))

  if (val.gt.0d0) then
! irrevesible part exists

    tmp1=matmul(chim,transpose(chim))/normchi+matmul(cM,xsiM)
    tmp2=Id+matmul(matmul(chim,transpose(chim))/normchi,(Id-alphaM)) & 
        +matmul(cM,matmul(xsiM,(Id-alphaM)))
    call inv3(tmpinv,tmp2)
    K=1d0/mu0*matmul(tmpinv,tmp1)
    
  else
! reversible
! all matrices are diagonal

    tmpv(1)=1d0/(1d0+cx*xsix*(1d0-alphax))
    tmpv(2)=1d0/(1d0+cy*xsiy*(1d0-alphay))
    tmpv(3)=1d0/(1d0+cz*xsiz*(1d0-alphaz))
    K=0d0
    K(1,1)=tmpv(1)*cx*xsix/mu0
    K(2,2)=tmpv(2)*cy*xsiy/mu0
    K(3,3)=tmpv(3)*cz*xsiz/mu0

!    tmp2=Id+matmul(cM,matmul(xsiM,Id-alphaM))
!    call inv3(tmpinv,tmp2)
!    K=1d0/mu0*matmul(tmpinv,matmul(cM,xsiM))

  end if

end subroutine magmater_javtangent


subroutine magmater_getjavH(field,Hout,iel) 
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
    if (size(Hout,1).eq.3) then
      Hout=Bsave(:,:,iel,new)-mu0*Hsave(:,:,iel,new)
    else
      Hout=Bsave(1:2,:,iel,new)-mu0*Hsave(1:2,:,iel,new)
    end if
  else
    write(*,*)'Field is note stored in mater_magnjav'
    write(*,*)field
    stop
  end if

      
end subroutine magmater_getjavH


subroutine magmater_javpara_get(prop,iel)
implicit none

double precision                  :: prop(:)
integer                           :: iel

   prop=propEM(matId(iel),:)

end subroutine magmater_javpara_get


end module  mater_magnjav