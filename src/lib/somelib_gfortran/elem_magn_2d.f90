module elem_magn_2d

! Last rev.:
!  M. Ristinmaa 2020-06-15
!   - initial version
!  M. Ristinmaa 2021-02-01
!   - Bugg in magnaxi4_B_rotA r=eci(gp,1->2) 

!
!------------------------------------------------------------------------------

! Module elem_magn_2d contains element subroutines for
! electormagnetics based on the magnetic vector potential A where
! B= \nabla x A
!
!  3-node element
!
!            3
!           / \
!          /   \
!         /     \
!        /       \
!       1---------2
!
!  4-node element
!
!      node location      gauss point location
!
!       4---------3                         eta (z-dir)
!       |         |            3----4        |
!       |         |            |    |        o-- xsi (r-dir)
!       |         |            1----2
!       1---------2
!
!------------------------------------------------------------------------------

implicit none

! This is an attempt to predefine some vectors for the 4-node element,
! to obtain speed. Note weights = 1d0
double precision  xsi4(4), eta4(4), G1, Ifm(9)
parameter        (G1=0.577350269189626D0)
parameter        (xsi4=(/-1D0,  1D0,-1D0,  1D0/)*G1 )
parameter        (eta4=(/-1D0, -1D0, 1D0,  1D0/)*G1 )
!parameter        (xsi4=(/ 1D0, 1D0,-1D0, -1D0/)*G1 ) ! Martins kod, gjort så att samma N_glo
!parameter        (eta4=(/-1D0, 1D0, 1D0, -1D0/)*G1 )
!parameter        (xsi4=(/-1D0,  1D0, 1D0, -1D0/)*G1 )  !test
!parameter        (eta4=(/-1D0, -1D0, 1D0,  1D0/)*G1 )

private xsi4, eta4, G1

double precision  NR4(4,4) ! 4 Gauss points and 4 nodes
data             (NR4(:,1)= (1D0-XSI4)*(1D0-ETA4)/4D0)
data             (NR4(:,2)= (1D0+XSI4)*(1D0-ETA4)/4D0)
data             (NR4(:,3)= (1D0+XSI4)*(1D0+ETA4)/4D0)
data             (NR4(:,4)= (1D0-XSI4)*(1D0+ETA4)/4D0)
private NR4

double precision  DNR4(8,4) ! Maybe one should use DNR(4,8) instead for speed
! derivate of shape functions with respect to xsi
data             (DNR4(1:8:2,1)= -(1D0-ETA4)/4D0)
data             (DNR4(1:8:2,2)=  (1D0-ETA4)/4D0)
data             (DNR4(1:8:2,3)=  (1D0+ETA4)/4D0)
data             (DNR4(1:8:2,4)= -(1D0+ETA4)/4D0)
! derivate of shape functions with respect to eta
data             (DNR4(2:8:2,1)=-(1D0-XSI4)/4D0)
data             (DNR4(2:8:2,2)=-(1D0+XSI4)/4D0)
data             (DNR4(2:8:2,3)= (1D0+XSI4)/4D0)
data             (DNR4(2:8:2,4)= (1D0-XSI4)/4D0)
private DNR4

! Interfaces
! ----------

interface magnaxi4_e
   module procedure magnaxi4_e1
   module procedure magnaxi4_e2
end interface
private :: magnaxi4_e1


interface magnaxi4_m
   module procedure magnaxi4_m1
   module procedure magnaxi4_m2
end interface
private :: magnaxi4_m1, magnaxi4_m2

private det2, inv2

!------------------------------------------------------------------------------
contains
!------ 4-node element axisymmetry

 subroutine magnaxi4_e1(Ke,Fe,coord,mu,Js)
! Symmetric around the x-axis
!
   implicit none

   double precision, intent(in)  :: mu(:), coord(:,:), Js(:)
   double precision, intent(out) :: Ke(:,:), Fe(:)

   double precision              :: JT(8,2), JTinv(2,2)
   double precision              :: Bm(2,4), dNdr(1,4), dNdz(1,4), N(1,4), Bsys(2,4)
   double precision              :: DetJ, eci(4,2), r
   integer                       :: gp, indx(2)
   integer, parameter            :: NGP=4

   eci = matmul(NR4,transpose(coord))  ! r-coord found in eci(:,2)
   JT=MATMUL(DNR4,transpose(coord))
   Ke=0d0
   Fe=0d0
   DO gp=1,NGP

     indx = (/ 2*gp-1, 2*gp /)
     detJ = det2(JT(indx,:))
     JTinv = JT(indx,:)
     call inv2(JTinv)
     Bm = matmul(JTinv,dNR4(indx,:))

     dNdz(1,:)=Bm(1,:)
     dNdr(1,:)=Bm(2,:)
     N(1,:)=NR4(gp,:)
     r=eci(gp,2)

     Bsys(1,:)=-dNdz(1,:)
     Bsys(2,:)=N(1,:)/r+dNdr(1,:)

     Ke=Ke+1d0/mu(gp)*MATMUL(TRANSPOSE(Bsys),Bsys)*detJ*r

!     Ke=Ke+1d0/mu(gp)*(MATMUL(TRANSPOSE(dNdr),dNdr) &
!            +MATMUL(TRANSPOSE(dNdz),dNdz) &
!            +1d0/r*(MATMUL(TRANSPOSE(dNdr),N) + MATMUL(TRANSPOSE(N),dNdr) &
!            +1d0/r*MATMUL(TRANSPOSE(N),N)))*detJ*r

     !write(*,*)'J ',Js(gp)
     Fe=Fe+Js(gp)*N(1,:)*detJ*r

   END DO

   RETURN
 end subroutine magnaxi4_e1


 subroutine magnaxi4_e2(Ke,Fe,coord,muinv,Js)
! Symmetric around the x-axis
!
   implicit none

   double precision, intent(in)  :: muinv(:,:,:), coord(:,:), Js(:)
   double precision, intent(out) :: Ke(:,:), Fe(:)

   double precision              :: JT(8,2), JTinv(2,2)
   double precision              :: Bm(2,4), dNdr(1,4), dNdz(1,4), N(1,4), Bsys(2,4)
   double precision              :: DetJ, eci(4,2), r
   integer                       :: gp, indx(2)
   integer, parameter            :: NGP=4

   eci = matmul(NR4,transpose(coord))  ! r-coord found in eci(:,2)
   JT=MATMUL(DNR4,transpose(coord))
   Ke=0d0
   Fe=0d0
   DO gp=1,NGP

     indx = (/ 2*gp-1, 2*gp /)
     detJ = det2(JT(indx,:))
     JTinv = JT(indx,:)
     call inv2(JTinv)
     Bm = matmul(JTinv,dNR4(indx,:))

     dNdz(1,:)=Bm(1,:)
     dNdr(1,:)=Bm(2,:)
     N(1,:)=NR4(gp,:)
     r=eci(gp,2)

     Bsys(1,:)=-dNdz(1,:)
     Bsys(2,:)=N(1,:)/r+dNdr(1,:)

     Ke=Ke+MATMUL(TRANSPOSE(Bsys),MATMUL(muinv(:,:,gp),Bsys))*detJ*r

!     Ke=Ke+1d0/mu(gp)*(MATMUL(TRANSPOSE(dNdr),dNdr) &
!            +MATMUL(TRANSPOSE(dNdz),dNdz) &
!            +1d0/r*(MATMUL(TRANSPOSE(dNdr),N) + MATMUL(TRANSPOSE(N),dNdr) &
!            +1d0/r*MATMUL(TRANSPOSE(N),N)))*detJ*r

     !write(*,*)'J ',Js(gp)
     Fe=Fe+Js(gp)*N(1,:)*detJ*r

   END DO

   RETURN
 end subroutine magnaxi4_e2


 subroutine magnaxi4_m1(Me,coord,sigma)
 ! Symmetric around the x-axis
 ! Note derivative d\mu/dA is not included in the tangent stiffness matrix
 !
    implicit none

    double precision, intent(in)  :: sigma(:,:,:), coord(:,:)
    double precision, intent(out) :: Me(:,:)

    double precision              :: JT(8,2)!, JTinv(2,2)
    double precision              :: N(1,4)
    double precision              :: DetJ, eci(4,2), r
    integer                       :: gp, indx(2)
    integer, parameter            :: NGP=4

    eci = matmul(NR4,transpose(coord))  ! r-coord found in eci(:,2)
    JT=MATMUL(DNR4,transpose(coord))

    Me=0d0

    DO gp=1,NGP

      indx = (/ 2*gp-1, 2*gp /)
      detJ = det2(JT(indx,:))
  !    JTinv = JT(indx,:)

      N(1,:)=NR4(gp,:)
      r=eci(gp,2)

      Me=Me+MATMUL(TRANSPOSE(N),MATMUL(sigma(:,:,gp),N))*detJ*r

    END DO

    RETURN
 end subroutine magnaxi4_m1


 subroutine magnaxi4_m2(Me,coord,sigma)
 ! Symmetric around the x-axis
 ! Note derivative d\mu/dA is not included in the tangent stiffness matrix
 !
    implicit none

    double precision, intent(in)  :: sigma(:), coord(:,:)
    double precision, intent(out) :: Me(:,:)

    double precision              :: JT(8,2)!, JTinv(2,2)
    double precision              :: N(1,4)
    double precision              :: DetJ, eci(4,2), r
    integer                       :: gp, indx(2)
    integer, parameter            :: NGP=4

    eci = matmul(NR4,transpose(coord))  ! r-coord found in eci(:,2)
    JT=MATMUL(DNR4,transpose(coord))

    Me=0d0

    DO gp=1,NGP

      indx = (/ 2*gp-1, 2*gp /)
      detJ = det2(JT(indx,:))
  !    JTinv = JT(indx,:)

      N(1,:)=NR4(gp,:)
      r=eci(gp,2)

      Me=Me+sigma(gp)*MATMUL(TRANSPOSE(N),N)*detJ*r

    END DO

    RETURN
 end subroutine magnaxi4_m2



 subroutine magnaxi4_fi(Fi,coord,H)
! Symmetric around the x-axis
! Note derivative d\mu/dA is not included in the tangent stiffness matrix
!
   implicit none

   double precision, intent(in)  :: coord(:,:), H(:,:)
   double precision, intent(out) :: Fi(:)

   double precision              :: JT(8,2), JTinv(2,2)
   double precision              :: Bm(2,4), dNdr(1,4), dNdz(1,4), N(1,4), Bsys(2,4)
   double precision              :: DetJ, eci(4,2), r
   integer                       :: gp, indx(2)
   integer, parameter            :: NGP=4

   eci = matmul(NR4,transpose(coord))  ! r-coord found in eci(:,2)
   JT=MATMUL(DNR4,transpose(coord))
   Fi=0d0
   DO gp=1,NGP

     indx = (/ 2*gp-1, 2*gp /)
     detJ = det2(JT(indx,:))
     JTinv = JT(indx,:)
     call inv2(JTinv)
     Bm = matmul(JTinv,dNR4(indx,:))

     dNdz(1,:)=Bm(1,:)
     dNdr(1,:)=Bm(2,:)
     N(1,:)=NR4(gp,:)
     r=eci(gp,2)
     
     Bsys(1,:)=-dNdz(1,:)
     Bsys(2,:)=N(1,:)/r+dNdr(1,:)

     Fi=Fi+matmul(transpose(Bsys),H(:,gp))*detJ*r

   END DO

   RETURN
 end subroutine magnaxi4_fi



  subroutine magnaxi4_B_rotA(B,coord,A)
  ! Symmetric around the x-axis
  ! Calculate magnetic induction B = rot(A)
  !
     implicit none

     double precision, intent(in)  :: A(:), coord(:,:)   ! nodal values
     double precision, intent(out) :: B(:,:) ! gauss point values

     double precision              :: JT(8,2), JTinv(2,2)
     double precision              :: Bm(2,4), dNdr(4), dNdz(4), N(4)
     double precision              :: DetJ, eci(4,2), r
     integer                       :: gp, indx(2)
     integer, parameter            :: NGP=4

     eci = matmul(NR4,transpose(coord))  ! r-coord found in eci(:,2)
     JT = matmul(DNR4,transpose(coord))


     DO gp=1,NGP

       indx = (/ 2*gp-1, 2*gp /)
       detJ = det2(JT(indx,:))
       JTinv = JT(indx,:)
       call inv2(JTinv)
       Bm = matmul(JTinv,dNR4(indx,:))

       dNdz(:)=Bm(1,:)
       dNdr(:)=Bm(2,:)
       N(:)=NR4(gp,:)
       r=eci(gp,2)

       B(1,gp) = -DOT_PRODUCT(dNdz,A)
       B(2,gp) = (DOT_PRODUCT(N,A)/r + DOT_PRODUCT(dNdr,A))

     END DO


  end subroutine magnaxi4_B_rotA


  subroutine magnaxi4_source(Q,cond,Ae)
     implicit none

     double precision, intent(in)  :: Ae(:), cond(:)   
     double precision, intent(out) :: Q(:) ! gauss point values

     double precision              :: N(4), gpval
     integer                       :: gp
     integer, parameter            :: NGP=4


     DO gp=1,NGP

       N(:)=NR4(gp,:)
       gpVal=DOT_PRODUCT(N,Ae) !Calculate the gauss point value
       Q(gp) = cond(gp)*(gpVal*gpVal)
      
     END DO

  end subroutine magnaxi4_source
  



! ----- help routines

function DET2(A)
  IMPLICIT NONE
  DOUBLE PRECISION                :: DET2
  DOUBLE PRECISION                :: A(:,:)

  DET2=A(1,1)*A(2,2)-A(2,1)*A(1,2)

RETURN
END function det2


subroutine inv2(A)
  implicit none
  double precision                :: A(:,:), Atmp(2,2)

  double precision                :: detA

  detA=A(1,1)*A(2,2)-A(1,2)*A(2,1)
  Atmp(1,:)=(/ A(2,2),-A(1,2)/)/detA
  Atmp(2,:)=(/-A(2,1), A(1,1)/)/detA
  A=Atmp

return
end subroutine inv2


end module elem_magn_2d
