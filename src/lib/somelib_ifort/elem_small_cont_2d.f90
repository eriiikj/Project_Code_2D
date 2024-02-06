module elem_small_cont_2d

! Last rev.:
!  H. Hallberg, 2011-04-26
!  M. Ristinmaa 2011-04-27
!   - Change order of variables in call
!   - s-command only calculates the strains
!  M. Ristinmaa 2011-04-28
!   - Renamed module, renamed elements inspired by Abaqus
!  M. Ristinmaa 2011-04-29
!   - included c2d3 element
!  M. Ristinmaa 2011-04-30
!   - renamed module lin to small
!  M. Ristinmaa 2012-06-11
!   - Bugg, Hooke plane strain change to 3x3 matrix according to manual
!  M. Ristinmaa 2020-03-12
!   - cleaning of 4-node element
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

interface c2d3_e
   module procedure plante
end interface

interface c2d3_s
   module procedure plants
end interface

interface c2d3_f
   module procedure plantf
end interface

interface c2d4_e
   module procedure plani4e_Ke
   module procedure plani4e_Ke1
end interface
private :: plani4e_Ke, plani4e_Ke1

interface c2d4_s
   module procedure plani4s_eset
   module procedure plani4s_eseteci
end interface
private :: plani4s_eset, plani4s_eseteci

private det2, inv2

!------------------------------------------------------------------------------
contains


!------ 3-node element

subroutine plante(Ke,elcoord,ep,D)
  implicit none
  double precision        :: Ke(6,6), elcoord(:,:)
  double precision        :: x(3), y(3), ep, D(3,3)

  double precision        :: detA, idetA, t, A
  double precision        :: B(3,6)

  x=elcoord(1,:)
  y=elcoord(2,:)
  detA=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detA*0.5D0
  idetA=1d0/detA

  t=ep

  B=reshape((/y(2)-y(3),0d0      ,y(3)-y(1),0d0      ,y(1)-y(2),0d0, &
            0d0      ,x(3)-x(2),0d0      ,x(1)-x(3),0d0      ,x(2)-x(1), &
            x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)/), &
           [3,6],order=[2,1])*idetA

  Ke=matmul(matmul(transpose(B),D),B)*A*t

return
end subroutine plante

subroutine plants(ee,elcoord,ed)
  implicit none
  double precision        :: elcoord(:,:)
  double precision        :: x(3), y(3), ee(3), ed(6)

  double precision        :: detA, idetA, A
  double precision        :: B(3,6)

  x=elcoord(1,:)
  y=elcoord(2,:)
  detA=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detA*0.5D0
  idetA=1d0/detA


  B=reshape((/y(2)-y(3),0d0      ,y(3)-y(1),0d0      ,y(1)-y(2),0d0, &
            0d0      ,x(3)-x(2),0d0      ,x(1)-x(3),0d0      ,x(2)-x(1), &
            x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)/), &
           [3,6],order=[2,1])*idetA!, magnaxid4_e2

  ee=matmul(B,ed)

return
end subroutine plants

subroutine plantf(fe,elcoord,ep,es)
  implicit none
  double precision        :: fe(6), elcoord(:,:), es(3)
  double precision        :: x(3), y(3), ep

  double precision        :: detA, idetA, t, A
  double precision        :: B(3,6)

  x=elcoord(1,:)
  y=elcoord(2,:)
  detA=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detA*0.5D0
  idetA=1d0/detA

  t=ep;

  B=reshape((/y(2)-y(3),0d0      ,y(3)-y(1),0d0      ,y(1)-y(2),0d0, &
            0d0      ,x(3)-x(2),0d0      ,x(1)-x(3),0d0      ,x(2)-x(1), &
            x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)/), &
           [3,6],order=[2,1])*idetA

  fe=matmul(transpose(B),es)*A*t

return
end subroutine plantf


!------ 4-node element

subroutine plani4e_Ke(Ke,coord,ep,D)
   implicit none

   double precision, intent(in)  :: D(:,:,:), coord(:,:), ep
   double precision, intent(out) :: Ke(:,:)

   double precision              :: JT(8,2), JTinv(2,2), DNX(2,4), B(3,8), Dgp(3,3)
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
     dNx = matmul(JTinv,dNR4(indx,:))

     B=0D0
     B(1,1:8:2)=dNx(1,:)
     B(2,2:8:2)=dNx(2,:)
     B(3,1:8:2)=dNx(2,:)
     B(3,2:8:2)=dNx(1,:)

     Dgp=D(:,:,gp_nr)

     KE=KE+MATMUL(TRANSPOSE(B),MATMUL(Dgp,B))*DETJ*ep

   END DO

   RETURN
 end subroutine plani4e_Ke

subroutine plani4e_Ke1(Ke,elcoord,ep,D)
   implicit none

   double precision, intent(in)  :: D(:,:), elcoord(:,:), ep
   double precision, intent(out) :: Ke(8,8)
   double precision              :: Dg(3,3,4)

   integer                       :: ii
   integer, parameter            :: NGP=4

   do ii=1,ngp
     Dg(:,:,ii)=D
   enddo
   call plani4e_Ke(Ke,elcoord,ep,Dg)

   return

end subroutine plani4e_Ke1



subroutine plani4s_eset(et,elcoord,ed)
   implicit none

   double precision, intent(in)    :: elcoord(:,:), ed(:)
   double precision, intent(inout) :: et(:,:)
   double precision                :: eci(4,2)

   call plani4s_eseteci(et,elcoord,ed,eci)

   return

end subroutine plani4s_eset


subroutine plani4s_eseteci(et,coord,ed,eci)
   implicit none

   double precision, intent(in)  :: coord(:,:), ed(:)
   double precision, intent(out) :: eci(:,:), et(:,:)

   double precision              :: JT(8,2), JTinv(2,2), DNX(2,4), B(3,8)
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
     dNx = matmul(JTinv,dNR4(indx,:))

     B=0D0
     B(1,1:8:2)=dNx(1,:)
     B(2,2:8:2)=dNx(2,:)
     B(3,1:8:2)=dNx(2,:)
     B(3,2:8:2)=dNx(1,:)

     et(:,gp_nr) = matmul(B,ed)

   END DO

   return
end subroutine plani4s_eseteci


subroutine c2d4_f(fint,coord,ep,es)
   implicit none

   double precision, intent(in)  :: coord(:,:), ep, es(:,:)
   double precision, intent(out) :: fint(:)


   double precision              :: JT(8,2), JTinv(2,2), DNX(2,4), B(3,8)
   double precision              :: DetJ, stress(3)
   integer                       :: GP_NR, indx(2)
   integer, parameter            :: NGP=4


   JT=MATMUL(DNR4,transpose(coord))
   fint=0D0
   DO GP_NR=1,NGP

     indx = (/ 2*gp_nr-1, 2*gp_nr /)
     detJ = det2(JT(indx,:))
     JTinv = JT(indx,:)
     call inv2(JTinv)
     dNx = matmul(JTinv,dNR4(indx,:))

     B=0D0
     B(1,1:8:2)=dNx(1,:)
     B(2,2:8:2)=dNx(2,:)
     B(3,1:8:2)=dNx(2,:)
     B(3,2:8:2)=dNx(1,:)

 !    stress = es((/1,2,4/),gp_nr)
     stress = es(:,gp_nr)
     fint = fint + matmul(transpose(B),stress)*detJ*ep

   END DO

   return
end subroutine c2d4_f

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


end module elem_small_cont_2d
