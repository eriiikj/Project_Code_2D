module mater_hyperel

! last modified
! M. Ristinmaa 2011-05-06
!  - hyper elastic material implemented
! M. Ristinmaa 2011-05-06
!  - pushforward and pullback implemented
!  - upddefgr implemented
! M. Ristinmaa 2011-05-06
!  - pushforward, pullback and updefgr debugged for 2-d
!  - 3-d parts missing in pushforward and pullback for matrices
! M. Ristinmaa 2011-05-06
!  - stype introduced
! M. Ristinmaa 2011-05-16
!  - introduced new module hyper elastic material derived from mater_large
! M. Ristinmaa 2011-05-23
!  - Initiation of material parameters
! M. Ristinmaa 2011-09-07
!  - Introduced isochoric volumetric split for use with mixed elements
! M. Ristinmaa 2011-09-12
!  - Bugg dneohooke2 - missed pushforward in 3D
!  - Debugged neohooeXup and dneohookeXup
! M. Ristinmaa 2011-09-12
!  - Changed order stress and strain components
!    [11 22 33 12 23 13]
! M. Wallin 2014-04-07
!  - Added functions for calculating the strain energy
! M. Ristinmaa 2020-02-14
!  - Bugg size(f) should have been size(ef) when returning correct stress vetor
!    in stress calculation
! ------------------------------------------------------------------------------


use mater_large
use matrix_util
use some_constants

implicit none

double precision                  :: Kmp, Gmp
private Kmp, Gmp

interface dneohooke
  module procedure dneohooke1
  module procedure dneohooke2
  module procedure dneohooke1up
  module procedure dneohooke2up
  module procedure dneohooke3up
end interface

private dneohooke1, dneohooke2

interface neohookeEnergy
  module procedure neohookeEnergy1
  module procedure neohookeEnergy2
end interface


interface neohooke
  module procedure neohooke1
  module procedure neohooke2
  module procedure neohooke1up
  module procedure neohooke2up
  module procedure neohooke3up
end interface

private neohooke1, neohooke2

interface dineohooke
  module procedure dineohooke1
  module procedure dineohooke1a
  module procedure dineohooke2
end interface

private dineohooke1, dineohooke2


!private inv3, inv3s, det3

!-------------------------------------------------------------------------------

contains

subroutine neohooke_init(mp)
  implicit none
  double precision                :: mp(:)
  Kmp=mp(1)
  Gmp=mp(2)

  return
end subroutine neohooke_init

! Material routines for a pure displacement formulation

! Call for several gauss points
subroutine neohooke1(stype,stress,ef)
  implicit none
  double precision                :: stress(:,:), ef(:,:)
  character(len=*)                :: stype
  integer                         :: gp, nr_gp

  nr_gp=size(ef,2)

  do gp=1,nr_gp
    call neohooke2(stype,stress(:,gp), ef(:,gp))
  enddo

  return
end subroutine neohooke1

! Call for one gauss point
subroutine neohooke2(stype,stress,ef)
  implicit none
  double precision                :: stress(:), ef(:)
  character(len=*)                :: stype
  
  double precision                :: E, v
  double precision                :: F(3,3), C(3,3), iC(3,3), id(3,3)
  double precision                :: J, trC, S(3,3), invF(3,3)
  
  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif

! Check if a inverse deformation gradient is provided as input
  if (stype.eq.'Cauchy-inv') then
    call inv3(invF,F)
    F=invF
  end if
  
  C=matmul(transpose(F),F)
  call inv3s(iC,C)
  J=det3(F)
  Id=0d0
  Id(1,1)=1d0
  Id(2,2)=1d0
  Id(3,3)=1d0
  trC=C(1,1)+C(2,2)+C(3,3)

!write(*,*)'stype ',stype 
  S=Kmp/2d0*(J**2d0-1d0)*iC+Gmp*J**(-2d0/3d0)*(Id-trC/3d0*iC)

  if ((stype.eq.'Cauchy').or.(stype.eq.'Cauchy-inv')) then
    S=matmul(F,matmul(S,transpose(F)))/J
  elseif (stype.eq.'2ndPiola') then
  else
    stop 'stype not implemented'
  endif

  ! if (size(ef).eq.4) then
  !   stress=(/S(1,1),S(2,2),S(3,3),S(1,2)/)
  ! else
  !   stress=(/S(1,1),S(2,2),S(3,3),S(1,2),S(1,3),S(2,3)/)
  ! endif

  ! EJ overwrite, only 4 stress outputs
  stress=(/S(1,1),S(2,2),S(3,3),S(1,2)/)

 

  return
end subroutine neohooke2


! Routines for a mixed u-p-J formulation

! Call for several gauss points
subroutine neohooke1up(stype,stress,ef,th)
  implicit none
  double precision                :: stress(:,:), ef(:,:), th(:)
  character(len=*)                :: stype
  integer                         :: gp, nr_gp

  nr_gp=size(ef,2)

  do gp=1,nr_gp
    call neohooke2up(stype,stress(:,gp), ef(:,gp),th(gp))
  enddo

  return
end subroutine neohooke1up

! Call for one gauss point
subroutine neohooke2up(stype,stress,ef,th)
  implicit none
  double precision                :: stress(:), ef(:), th
  character(len=*)                :: stype
  
  double precision                :: E, v
  double precision                :: F(3,3), C(3,3), iC(3,3), id(3,3), invF(3,3)
  double precision                :: J, trC, S(3,3), b(3,3), trb, p
  
  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif

! Check if a inverse deformation gradient is provided as input
  if (stype.eq.'Cauchy-inv') then
    call inv3(invF,F)
    F=invF
  end if
  
  J=det3(F)
  Id=0d0
  Id(1,1)=1d0
  Id(2,2)=1d0
  Id(3,3)=1d0

  p=Kmp/2d0*(th-1d0/th)

  if ((stype.eq.'Cauchy').or.(stype.eq.'Cauchy-inv')) then
    b=J**(-2d0/3d0)*matmul(F,transpose(F))
    trb=b(1,1)+b(2,2)+b(3,3)
    b=b-trb/3d0*Id
    S=Gmp*b/J+p*id
  elseif (stype.eq.'2ndPiola') then
    C=matmul(transpose(F),F)
    call inv3s(iC,C)
    trC=C(1,1)+C(2,2)+C(3,3)
    S=p*J*iC+Gmp*J**(-2d0/3d0)*(Id-trC/3d0*iC)
  else
    stop 'stype not implemented'
  endif
  if (size(ef).eq.4) then
    stress=(/S(1,1),S(2,2),S(3,3),S(1,2)/)
  else
    stress=(/S(1,1),S(2,2),S(3,3),S(1,2),S(1,3),S(2,3)/)
  endif

  return
end subroutine neohooke2up

! Call for constant th
subroutine neohooke3up(stype,stress,ef,th)
  implicit none
  double precision                :: stress(:,:), ef(:,:), th
  character(len=*)                :: stype
  integer                         :: gp, nr_gp

  nr_gp=size(ef,2)

  do gp=1,nr_gp
    call neohooke2up(stype,stress(:,gp), ef(:,gp),th)
  enddo

  return
end subroutine neohooke3up



! Material tangent stiffness for a pure displacement formulation

subroutine dneohooke1(stype,D,ef)
  implicit none
  double precision                :: D(:,:,:), ef(:,:)
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(D,3)

  do gp=1,nr_gp
    call dneohooke2(stype,D(:,:,gp),ef(:,gp))
  enddo

  return
end subroutine dneohooke1


subroutine dneohooke2(stype,D,ef)
  implicit none
  double precision                :: D(:,:), ef(:)
  character(len=*)                :: stype
  
  double precision                :: E, v, K, G
  double precision                :: F(3,3), C(3,3), iC(3,3), Id(3,3)
  double precision                :: J, trC, a1, a2, a3, Ds(21)
  integer                         :: i, jj, kk, l, in(21,4), el
  
  K=Kmp
  G=Gmp

 
  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif
  
  C=matmul(transpose(F),F)
  call inv3s(iC,C)
  J=det3(F)
  Id=0d0
  Id(1,1)=1d0
  Id(2,2)=1d0
  Id(3,3)=1d0
  trC=C(1,1)+C(2,2)+C(3,3)

  a1=K*J**2d0+2d0/9d0*G*J**(-2d0/3d0)*trC
  a2=2d0*G/3d0*J**(-2d0/3d0)
  a3=G/3d0*J**(-2d0/3d0)*trC-K/2d0*(J**2d0-1d0)

  if (size(ef).eq.4) then
    in(1,:)=(/1, 1, 1, 1/)
    in(2,:)=(/1, 1, 2, 2/)
    in(3,:)=(/1, 1, 1, 2/)
    in(4,:)=(/2, 2, 2, 2/)
    in(5,:)=(/2, 2, 1, 2/)
    in(6,:)=(/1, 2, 1, 2/)

    do el=1,6
      i=in(el,1)
      jj=in(el,2)
      kk=in(el,3)
      l=in(el,4)

      Ds(el)=a1*iC(i,jj)*iC(kk,l)-a2*(Id(i,jj)*iC(kk,l)+iC(i,jj)*Id(kk,l)) &
       +a3*(iC(i,kk)*iC(jj,l)+iC(i,l)*iC(jj,kk))
    enddo

    D(1,:)=(/Ds(1), Ds(2), Ds(3)/)
    D(2,:)=(/Ds(2), Ds(4), Ds(5)/)
    D(3,:)=(/Ds(3), Ds(5), Ds(6)/)

  else
    in(1,:) =(/1, 1, 1, 1/)
    in(2,:) =(/1, 1, 2, 2/)
    in(3,:) =(/1, 1, 3, 3/)
    in(4,:) =(/1, 1, 1, 2/)
    in(5,:) =(/1, 1, 1, 3/)
    in(6,:) =(/1, 1, 2, 3/)    
    in(7,:) =(/2, 2, 2, 2/)    
    in(8,:) =(/2, 2, 3, 3/)    
    in(9,:) =(/2, 2, 1, 2/)    
    in(10,:)=(/2, 2, 1, 3/)    
    in(11,:)=(/2, 2, 2, 3/)    
    in(12,:)=(/3, 3, 3, 3/)    
    in(13,:)=(/3, 3, 1, 2/)    
    in(14,:)=(/3, 3, 1, 3/)    
    in(15,:)=(/3, 3, 2, 3/)    
    in(16,:)=(/1, 2, 1, 2/)    
    in(17,:)=(/1, 2, 1, 3/)    
    in(18,:)=(/1, 2, 2, 3/)    
    in(19,:)=(/1, 3, 1, 3/)    
    in(20,:)=(/1, 3, 2, 3/)    
    in(21,:)=(/2, 3, 2, 3/)    

    do el=1,21
      i=in(el,1)
      jj=in(el,2)
      kk=in(el,3)
      l=in(el,4)

      Ds(el)=a1*iC(i,jj)*iC(kk,l)-a2*(Id(i,jj)*iC(kk,l)+iC(i,jj)*Id(kk,l)) &
       +a3*(iC(i,kk)*iC(jj,l)+iC(i,l)*iC(jj,kk))
    enddo

    ! D(1,:)=(/Ds(1),  Ds(2),  Ds(3),  Ds(4),  Ds(5),  Ds(6)/)
    ! D(2,:)=(/Ds(2),  Ds(7),  Ds(8),  Ds(9), Ds(10), Ds(11)/)
    ! D(3,:)=(/Ds(3),  Ds(8), Ds(12), Ds(13), Ds(14), Ds(15)/)
    ! D(4,:)=(/Ds(4),  Ds(9), Ds(13), Ds(16), Ds(17), Ds(18)/)
    ! D(5,:)=(/Ds(5), Ds(10), Ds(14), Ds(17), Ds(19), Ds(20)/)
    ! D(6,:)=(/Ds(6), Ds(11), Ds(15), Ds(18), Ds(20), Ds(21)/)

    ! EJ, overwrite output to global plane strain case
    D(1,:)=(/Ds(1), Ds(2), Ds(4)/)
    D(2,:)=(/Ds(2), Ds(7), Ds(9)/)
    D(3,:)=(/Ds(4), Ds(9), Ds(16)/)
  endif 

  if (stype.eq.'ul') then
    call pushforward(D,D,ef,'j')
  elseif (stype.eq.'tl') then
  else
    stop 'stype not implemented'
  endif

  


  return
end subroutine dneohooke2


! Material tangent stiffness for a mixed displacement formulation

subroutine dneohooke1up(stype,D,Dth,ef,th)
  implicit none
  double precision                :: D(:,:,:), Dth(:), ef(:,:), th(:)
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(D,3)

  do gp=1,nr_gp
    call dneohooke2up(stype,D(:,:,gp), Dth(gp),ef(:,gp),th(gp))
  enddo

  return
end subroutine dneohooke1up


! Call when constant th in element
subroutine dneohooke3up(stype,D,Dth,ef,th)
  implicit none
  double precision                :: D(:,:,:), Dth(:), ef(:,:), th
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(D,3)

  do gp=1,nr_gp
    call dneohooke2up(stype,D(:,:,gp), Dth(gp),ef(:,gp),th)
  enddo

  return
end subroutine dneohooke3up


subroutine dneohooke2up(stype,D,Dth,ef,th)
  implicit none
  double precision                :: D(:,:), Dth, ef(:), th
  character(len=*)                :: stype
  
  double precision                :: F(3,3), b(3,3), bt(3,3), Id(3,3)
  double precision                :: J, trb, a1, a2, a3, Ds(21)
  double precision                :: tmp, tmp1, tmp2, tmp3, taud(3,3)
  integer                         :: i, jj, kk, l, in(21,4), el
  
  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif
  
  b=matmul(F,transpose(F))
  J=det3(F)
  Id=0d0
  Id(1,1)=1d0
  Id(2,2)=1d0
  Id(3,3)=1d0
  trb=b(1,1)+b(2,2)+b(3,3)
  
  taud=Gmp*J**(-2d0/3d0)*(b-1d0/3d0*trb*Id)

if (1.eq.1) then
  D(1,:)=(/2d0*taud(1,1),taud(1,1)+taud(2,2),taud(1,1)+taud(3,3), &
         taud(1,2),taud(1,3),taud(2,3)/)
  D(2,:)=(/taud(2,2)+taud(1,1),2d0*taud(2,2),taud(2,2)+taud(3,3), &
         taud(1,2),taud(1,3),taud(2,3)/)
  D(3,:)=(/taud(3,3)+taud(1,1),taud(3,3)+taud(2,2),2d0*taud(3,3), &
         taud(1,2),taud(1,3),taud(2,3)/)
  D(4,:)=(/taud(1,2),taud(1,2),taud(1,2),0d0,0d0,0d0/)
  D(5,:)=(/taud(1,3),taud(1,3),taud(1,3),0d0,0d0,0d0/)
  D(6,:)=(/taud(2,3),taud(2,3),taud(2,3),0d0,0d0,0d0/)
  D=D*(-1d0/J*2d0/3d0)

  tmp=-Gmp*trb*J**(-2d0/3d0)/J*2d0/3d0
  tmp1=-4d0/6d0*tmp
  tmp2=1d0/3d0*tmp
  tmp3=-0.5d0*tmp
  D(1,:)=D(1,:)+(/tmp1, tmp2, tmp2,  0d0,  0d0,  0d0/)
  D(2,:)=D(2,:)+(/tmp2, tmp1, tmp2,  0d0,  0d0,  0d0/)
  D(3,:)=D(3,:)+(/tmp2, tmp2, tmp1,  0d0,  0d0,  0d0/)
  D(4,:)=D(4,:)+(/ 0d0,  0d0,  0d0, tmp3,  0d0,  0d0/)
  D(5,:)=D(5,:)+(/ 0d0,  0d0,  0d0,  0d0, tmp3,  0d0/)
  D(6,:)=D(6,:)+(/ 0d0,  0d0,  0d0,  0d0,  0d0, tmp3/)
else
    in(1,:) =(/1, 1, 1, 1/)
    in(2,:) =(/1, 1, 2, 2/)
    in(3,:) =(/1, 1, 3, 3/)
    in(4,:) =(/1, 1, 1, 2/)
    in(5,:) =(/1, 1, 1, 3/)
    in(6,:) =(/1, 1, 2, 3/)    
    in(7,:) =(/2, 2, 2, 2/)    
    in(8,:) =(/2, 2, 3, 3/)    
    in(9,:) =(/2, 2, 1, 2/)    
    in(10,:)=(/2, 2, 1, 3/)    
    in(11,:)=(/2, 2, 2, 3/)    
    in(12,:)=(/3, 3, 3, 3/)    
    in(13,:)=(/3, 3, 1, 2/)    
    in(14,:)=(/3, 3, 1, 3/)    
    in(15,:)=(/3, 3, 2, 3/)    
    in(16,:)=(/1, 2, 1, 2/)    
    in(17,:)=(/1, 2, 1, 3/)    
    in(18,:)=(/1, 2, 2, 3/)    
    in(19,:)=(/1, 3, 1, 3/)    
    in(20,:)=(/1, 3, 2, 3/)    
    in(21,:)=(/2, 3, 2, 3/)    


    tmp1=-1d0/J*2d0/3d0
    tmp2=tmp1*Gmp*trb*J**(-2d0/3d0)
    do el=1,21
      i=in(el,1)
      jj=in(el,2)
      kk=in(el,3)
      l=in(el,4)

      Ds(el)=tmp1*(taud(i,jj)*Id(kk,l)+Id(i,jj)*taud(kk,l)) &
            +tmp2*(1d0/3d0*Id(i,jj)*Id(kk,l) &
       -1d0/2d0*(Id(i,kk)*Id(jj,l)+Id(i,l)*Id(jj,kk)))
    enddo

    D(1,:)=(/Ds(1),  Ds(2),  Ds(3),  Ds(4),  Ds(5),  Ds(6)/)
    D(2,:)=(/Ds(2),  Ds(7),  Ds(8),  Ds(9), Ds(10), Ds(11)/)
    D(3,:)=(/Ds(3),  Ds(8), Ds(12), Ds(13), Ds(14), Ds(15)/)
    D(4,:)=(/Ds(4),  Ds(9), Ds(13), Ds(16), Ds(17), Ds(18)/)
    D(5,:)=(/Ds(5), Ds(10), Ds(14), Ds(17), Ds(19), Ds(20)/)
    D(6,:)=(/Ds(6), Ds(11), Ds(15), Ds(18), Ds(20), Ds(21)/)

!    write(*,*)'D(1,:) ',D(1,:)
!    write(*,*)'D(2,:) ',D(2,:)
!    write(*,*)'D(3,:) ',D(3,:)
!    write(*,*)'D(4,:) ',D(4,:)
!    write(*,*)'D(5,:) ',D(5,:)
!    write(*,*)'D(6,:) ',D(6,:)
!stop
endif

  if (stype.eq.'ul') then
  elseif (stype.eq.'tl') then
    call pullback(D,D,ef,'j')
  else
    stop 'stype not implemented'
  endif

  Dth=1d0/2d0*Kmp*(1d0+1d0/th/th)

  return
end subroutine dneohooke2up


! Tangent for inverse motion problem

subroutine dineohooke1(stype,Duu,es,ef,scale)
  implicit none
  double precision                :: Duu(:,:,:)
  double precision                :: es(:,:)
  double precision                :: ef(:,:)
  double precision                :: scale(:)
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(Duu,3)

  do gp=1,nr_gp
    call dineohooke2(stype,Duu(:,:,gp),es(:,gp),ef(:,gp),scale(gp))
  enddo

  return
end subroutine dineohooke1


subroutine dineohooke1a(stype,Duu,es,ef)
  implicit none
  double precision                :: Duu(:,:,:)
  double precision                :: es(:,:)
  double precision                :: ef(:,:)
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(Duu,3)

  do gp=1,nr_gp
    call dineohooke2(stype,Duu(:,:,gp),es(:,gp),ef(:,gp),1d0)
  enddo

  return
end subroutine dineohooke1a


subroutine dineohooke2(stype,Duu,es,ef,scale)
  implicit none
  double precision                :: Duu(:,:)
  double precision                :: es(:), ef(:)
  double precision, optional      :: scale
  character(len=*)                :: stype
  
  double precision                :: f(3,3), fi(3,3)
  integer                         :: i, j, k, l, ii, jj
  
  double precision                :: duun(6,9), Duux(6,6)
  double precision                :: efn(3), efp(9), s(3,3)

  integer, parameter              :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                  (/1,2/),(/1,3/),(/2,3/)]
  integer, parameter              :: f_list(2,9)=[(/1,1/),(/1,2/),(/1,3/),&
                                                  (/2,1/),(/2,2/),(/2,3/),&
                                                  (/3,1/),(/3,2/),(/3,3/)]

  if (size(ef).eq.4) then
    f(1,:)=(/ef(1), ef(2), zero/)
    f(2,:)=(/ef(3), ef(4), zero/)
    f(3,:)=(/ zero,  zero,  one/)

    s(1,:)=(/es(1),es(3), zero/)
    s(2,:)=(/es(4),es(2), zero/)
    s(3,:)=(/zero , zero, zero/) ! 33 will not influence

  else

    f(1,:)=(/ef(1), ef(2), ef(3)/)
    f(2,:)=(/ef(4), ef(5), ef(6)/)
    f(3,:)=(/ef(7), ef(8), ef(9)/)

    s(1,:)=(/es(1),es(4),es(5)/)
    s(2,:)=(/es(4),es(2),es(6)/)
    s(3,:)=(/es(5),es(6),es(3)/)

  endif

  call inv3(fi,f)
  do ii=1,6
    i=d_list(1,ii)
    j=d_list(2,ii)
    do jj=1,9
      k=f_list(1,jj)
      l=f_list(2,jj)
      duun(ii,jj)=s(i,j)*fi(l,k)-fi(i,k)*s(l,j)-s(i,l)*fi(j,k)
    end do
  end do  

! Deformation gradient
  call inv3(fi,f)
  efp=(/fi(1,1),fi(1,2),fi(1,3),fi(2,1),fi(2,2),fi(2,3),fi(3,1),fi(3,2),fi(3,3)/)

! Tanget updated Lagrange formulation
  call dneohooke('ul',duux,efp)
  duux=duux*scale
! Note that the matrices dup and dpp found from the updated Lagrange
! scheme are the same as for the inverse motion problem, only duu and dpu need
! to be evaluated

! Duu matrix for inverse motion
! order 11 12 13 21 22 23 31 32 33
  duun(:,1)=duun(:,1)-duux(:,1)*fi(1,1)-duux(:,4)*fi(2,1)-duux(:,5)*fi(3,1)
  duun(:,2)=duun(:,2)-duux(:,4)*fi(1,1)-duux(:,2)*fi(2,1)-duux(:,6)*fi(3,1)
  duun(:,3)=duun(:,3)-duux(:,5)*fi(1,1)-duux(:,6)*fi(2,1)-duux(:,3)*fi(3,1)

  duun(:,4)=duun(:,4)-duux(:,1)*fi(1,2)-duux(:,4)*fi(2,2)-duux(:,5)*fi(3,2)
  duun(:,5)=duun(:,5)-duux(:,4)*fi(1,2)-duux(:,2)*fi(2,2)-duux(:,6)*fi(3,2)
  duun(:,6)=duun(:,6)-duux(:,5)*fi(1,2)-duux(:,6)*fi(2,2)-duux(:,3)*fi(3,2)

  duun(:,7)=duun(:,7)-duux(:,1)*fi(1,3)-duux(:,4)*fi(2,3)-duux(:,5)*fi(3,3)
  duun(:,8)=duun(:,8)-duux(:,4)*fi(1,3)-duux(:,2)*fi(2,3)-duux(:,6)*fi(3,3)
  duun(:,9)=duun(:,9)-duux(:,5)*fi(1,3)-duux(:,6)*fi(2,3)-duux(:,3)*fi(3,3)

  if (size(ef).eq.4) then
    duu(1,:)=duun(1,[1, 2, 4, 5])
    duu(2,:)=duun(2,[1, 2, 4, 5])
    duu(3,:)=duun(4,[1, 2, 4, 5])
  else 
    duu=duun
  end if

  return
end subroutine dineohooke2


subroutine neohookeEnergyInv(energy,ef)
  implicit none
  double precision                :: ef(:,:),energy(:),gpenergy
  integer                         :: gp,nr_gp
  
  double precision                :: efp(9), f(3,3), fi(3,3)

  nr_gp=size(ef,2)

  do gp=1,nr_gp
    if (size(ef).eq.4) then
      f(1,:)=(/ef(1,gp), ef(2,gp), zero/)
      f(2,:)=(/ef(3,gp), ef(4,gp), zero/)
      f(3,:)=(/    zero,     zero,  one/)
    else
      f(1,:)=(/ef(1,gp), ef(2,gp), ef(3,gp)/)
      f(2,:)=(/ef(4,gp), ef(5,gp), ef(6,gp)/)
      f(3,:)=(/ef(7,gp), ef(8,gp), ef(9,gp)/)
    endif
! Deformation gradient
    call inv3(fi,f)
    efp=(/fi(1,1),fi(1,2),fi(1,3),fi(2,1),fi(2,2),fi(2,3),fi(3,1),fi(3,2),fi(3,3)/)
    call neohookeEnergy1(gpenergy,efp,gp)
    energy(gp)=gpenergy
  enddo
  
  return
end subroutine neohookeEnergyInv



subroutine neohookeEnergy2(energy,ef)
  implicit none
  double precision                :: ef(:,:),energy(:),gpenergy
  integer                         :: gp,nr_gp

  nr_gp=size(ef,2)

  do gp=1,nr_gp
    call neohookeEnergy1(gpenergy,ef(:,gp),gp)
    energy(gp)=gpenergy
  enddo
  
  return
end subroutine neohookeEnergy2


subroutine neohookeEnergy1(energy,ef,gp)
  implicit none
  double precision                :: energy, ef(:)
  
  double precision                :: E, v, D1
  double precision                :: F(3,3), C(3,3)
  double precision                :: J, trC

  integer                         :: gp

  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif
  
  C=matmul(transpose(F),F)
  J=det3(F)

  trC=C(1,1)+C(2,2)+C(3,3)

  energy=Kmp/2d0*(0.5D0*(J**2d0-1d0)-dlog(J))+0.5D0*Gmp*(J**(-2d0/3d0)*trC-3d0)


  return
end subroutine neohookeEnergy1





end module mater_hyperel
