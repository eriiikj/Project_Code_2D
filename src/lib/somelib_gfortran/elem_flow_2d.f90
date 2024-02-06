module elem_flow_2d

! Last rev.:
!  M. Ristinmaa 2020-02-28
!   - initial verison
!
!------------------------------------------------------------------------------

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
!       4---------3                         eta
!       |         |            3----4        |
!       |         |            |    |        o-- xsi
!       |         |            1----2
!       1---------2
!
!
!   Module elem_flow_2d contains element subroutines for
!   scalar problems
!
!------------------------------------------------------------------------------

implicit none

! This is an attempt to predefine some vectors for the 4-node element, to obtain speed
double precision  xsi4(4), eta4(4), G1, Ifm(9)
parameter        (G1=0.577350269189626D0)
parameter        (xsi4=(/-1D0,  1D0,-1D0,  1D0/)*G1 )
parameter        (eta4=(/-1D0, -1D0, 1D0,  1D0/)*G1 )

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

interface flow2d4_e
   module procedure flow2d4_e1
   module procedure flow2d4_e2
end interface
private :: flow2d4_e1, flow2d4_e2

interface flow2d4_s
!   module procedure flow2d4_s1
   module procedure flow2d4_s2
end interface
private :: flow2d4_s1, flow2d4_s2


!------------------------------------------------------------------------------
contains
!------ 4-node element axisymmetry

 subroutine flowaxi4_e(Ke,Fe,coord,kappa,Q)
! Symmetric around the x-axis
!
   implicit none

   double precision, intent(in)  :: kappa(:), coord(:,:), Q(:)
   double precision, intent(out) :: Ke(:,:), Fe(:)

   double precision              :: JT(8,2), JTinv(2,2)
   double precision              :: Bm(2,4), dNdr(1,4), dNdz(1,4), N(1,4)
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

     Ke=Ke+kappa(gp)*(MATMUL(TRANSPOSE(dNdr),dNdr) &
            +MATMUL(TRANSPOSE(dNdz),dNdz))*detJ*r

     Fe=Fe+Q(gp)*N(1,:)*detJ*r

   END DO

   RETURN
 end subroutine flowaxi4_e

 subroutine flowaxi4_m(Me,coord,para)
 ! Symmetric around the x-axis
 !
    implicit none

    double precision, intent(in)  :: para(:), coord(:,:)
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

      Me=Me+para(gp)*MATMUL(TRANSPOSE(N),N)*detJ*r

    END DO

    RETURN
  end subroutine flowaxi4_m



!------ 4-node element plane elements

subroutine flow2d4_e1(Ke,coord,ep,D)
   implicit none

   double precision, intent(in)  :: D(:,:,:), coord(:,:), ep
   double precision, intent(out) :: Ke(:,:)

   double precision              :: JT(8,2), JTinv(2,2), B(2,4), Dgp(2,2)
   double precision              :: DetJ
   integer                       :: GP_NR, indx(2)
   integer, parameter            :: NGP=4


   JT=MATMUL(DNR4,transpose(coord))
   Ke=0D0
   DO GP_NR=1,NGP

     indx = (/ 2*gp_nr-1, 2*gp_nr /)
     detJ = det2(JT(indx,:))
     JTinv = JT(indx,:)
     call inv2(JTinv)
     B = matmul(JTinv,dNR4(indx,:))

     Dgp=D(:,:,gp_nr)

     KE=KE+MATMUL(TRANSPOSE(B),MATMUL(Dgp,B))*DETJ*ep

   END DO

   RETURN
 end subroutine flow2d4_e1

subroutine flow2d4_e2(Ke,elcoord,ep,D)
   implicit none

   double precision, intent(in)  :: D(:,:), elcoord(:,:), ep
   double precision, intent(out) :: Ke(:,:)
   double precision              :: Dg(2,2,4)

   integer                       :: ii
   integer, parameter            :: NGP=4

   do ii=1,ngp
     Dg(:,:,ii)=D
   enddo
   call flow2d4_e1(Ke,elcoord,ep,Dg)

   return

end subroutine flow2d4_e2



subroutine flow2d4_s1(et,elcoord,ed)
   implicit none

   double precision, intent(in)    :: elcoord(:,:), ed(:)
   double precision, intent(inout) :: et(:,:)
   double precision                :: eci(4,2)

   call flow2d4_s2(et,elcoord,ed,eci)

   return

end subroutine flow2d4_s1


subroutine flow2d4_s2(et,coord,ed,eci)
   implicit none

   double precision, intent(in)  :: coord(:,:), ed(:)
   double precision, intent(out) :: eci(:,:), et(:,:)

   double precision              :: JT(8,2), JTinv(2,2), B(2,4)
   double precision              :: DetJ
   integer                       :: GP_NR, indx(2)
   integer, parameter            :: NGP=4


   eci = matmul(NR4,transpose(coord))
   JT=MATMUL(DNR4,transpose(coord))
   DO GP_NR=1,NGP

     indx = (/ 2*gp_nr-1, 2*gp_nr /)
     detJ = det2(JT(indx,:))
     JTinv = JT(indx,:)
     call inv2(JTinv)
     B = matmul(JTinv,dNR4(indx,:))

     et(:,gp_nr) = matmul(B,ed)

   END DO

   return
end subroutine flow2d4_s2


subroutine flow2d4_f(fint,coord,ep,es)
   implicit none

   double precision, intent(in)  :: coord(:,:), ep, es(:,:)
   double precision, intent(out) :: fint(:)


   double precision              :: JT(8,2), JTinv(2,2), B(2,4)
   double precision              :: DetJ, flux(2)
   integer                       :: GP_NR, indx(2)
   integer, parameter            :: NGP=4


   JT=MATMUL(DNR4,transpose(coord))
   fint=0D0
   DO GP_NR=1,NGP

     indx = (/ 2*gp_nr-1, 2*gp_nr /)
     detJ = det2(JT(indx,:))
     JTinv = JT(indx,:)
     call inv2(JTinv)
     B = matmul(JTinv,dNR4(indx,:))

     flux = es(:,gp_nr)
     fint = fint + matmul(transpose(B),flux)*detJ*ep

   END DO

   return
end subroutine flow2d4_f


!------ 3-node plane elements

subroutine flow2d3_e(Ke,elcoord,t,D)
  implicit none
  double precision        :: Ke(3,3), elcoord(:,:)
  double precision        :: t, D(2,2)

  double precision        :: detC, A, idetC, x(3), y(3)
  double precision        :: B(2,3)

  x=elcoord(1,:)
  y=elcoord(2,:)
  detC=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detC*0.5D0
  idetC=1d0/detC

  B=reshape((/y(2)-y(3), y(3)-y(1), y(1)-y(2), &
              x(3)-x(2), x(1)-x(3), x(2)-x(1)/), &
           [2,3],order=[2,1])*idetC

  Ke=matmul(matmul(transpose(B),D),B)*A*t

return
end subroutine flow2d3_e


subroutine flow2d3_m(Me,elcoord,t,rho)
  implicit none
  double precision        :: Me(3,3), elcoord(:,:)
  double precision        :: t, rho

  double precision        :: detC, A, x(3), y(3)

  x=elcoord(1,:)
  y=elcoord(2,:)
  detC=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detC*0.5D0

  Me(:,1)=(/2d0, 1d0, 1d0/)
  Me(:,2)=(/1d0, 2d0, 1d0/)
  Me(:,3)=(/1d0, 1d0, 2d0/)

  Me=rho*t*A/12D0*Me

return
end subroutine flow2d3_m


subroutine flow2d3_mlumped(Me,elcoord,t,rho)
  implicit none
  double precision        :: Me(3,3), elcoord(:,:)
  double precision        :: t, rho

  double precision        :: detC, A, x(3), y(3)

  x=elcoord(1,:)
  y=elcoord(2,:)
  detC=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detC*0.5D0

  Me(:,1)=(/1d0, 0d0, 0d0/)
  Me(:,2)=(/0d0, 1d0, 0d0/)
  Me(:,3)=(/0d0, 0d0, 1d0/)

  Me=rho*t*A/3D0*Me

return
end subroutine flow2d3_mlumped

subroutine flow2d3_f(fe,elcoord,t,val)
  implicit none
  double precision        :: fe(3), elcoord(:,:)
  double precision        :: t, val

  double precision        :: detC, A, x(3), y(3)

  x=elcoord(1,:)
  y=elcoord(2,:)
  detC=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detC*0.5D0

  fe=(/1d0, 1d0, 1d0/)*t*val*A/3d0

return
end subroutine flow2d3_f

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


end module elem_flow_2d
