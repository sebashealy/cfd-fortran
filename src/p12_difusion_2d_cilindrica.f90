PROGRAM Difusion_2D_Cil !En la vertical o eje 'y' va r y en la horizontal o eje 'x' va theta
IMPLICIT NONE

!Declaración de variables
INTEGER i, j, n0, nr, max_iter !El 0 esta representando a theta
REAL*4 d0, dr, theta0, thetal, r0, rl, dv, Se, Sw, Sn, Ss, tolerance, residual, rp, rn, rs, pi
REAL*4, ALLOCATABLE, DIMENSION(:,:) :: aE, aW, aN, aS, aP, sP, phi
!phi es una variable que puede representar a cualquier otra variable como a la temperatura
REAL*4, ALLOCATABLE, DIMENSION(:) :: theta, thetac, r, rc

!Inicialización de variables
n0=40; nr=40

pi=ACOS(-1.0)

theta0=0.25*pi; thetal=0.75*pi !De la presentación, estamos definiendo los límites
r0=0.1; rl=1.0 !r0=0.1 porque cuando es 0, se complican mucho las cosas

max_iter=1000
tolerance=1e-5

d0=(thetal-theta0)/float(n0)
dr=(rl-r0)/float(nr)

Se=dr; sW=dr

!Alojamiento de la memoria dinámica
ALLOCATE(aP(n0,nr), aE(n0, nr), aW(n0, nr), aN(n0, nr), aS(n0, nr), sP(n0, nr), phi(0:n0+1, 0:nr+1))
ALLOCATE(r(0:nr), rc(0:nr+1), theta(0:n0), thetac(0:n0+1))

!Creación de la malla
CALL MESH_1D(n0, theta0, thetal, theta, thetac)
CALL MESH_1D(nr, r0, rl, r ,rc)

!Condiciones de frontera
phi=0.0
phi(:, nr+1)=1.0

!Cálculo de Coeficientes
do i=1, n0
    do j=1, nr
         !Radios, porque dependiendo de la posición de volumen de control el radio o tamaño de norte y sur va a cambiar
         rp=rc(j) !Este es el valor del radio en el centro del volumen de control
         rn=r(j)
         rs=r(j-1)

         Sn=rn*d0
         Ss=rs*d0
         dv=rp*dr*d0

         aE(i,j)=Se/(rp*d0)
         aW(i,j)=Sw/(rp*d0)
         aN(i,j)=Sn/dr
         aS(i,j)=Ss/dr
         aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)
         sP(i,j)=0.0
    end do
end do

!Correcciones por condiciones de frontera

!Este, no se modifica aunque se haya cambiado de lado el E y W porque también se cambia el origen de 'x' en este caso theta, en el lado derecho tenemos 1 y en el izquierdo theta. Es así por como se miden los ángulos
aP(n0,:)=aP(n0,:)+aE(n0,:)
sP(n0,:)=sP(n0,:)+2.0*aE(n0,:)*phi(n0+1,1:nr)
aE(n0,:)=0.0
!Oeste
aP(1,:)=aP(1,:)+aW(1,:)
sP(1,:)=sP(1,:)+2.0*aW(1,:)*phi(0,1:nr)
aW(1,:)=0.0
!Norte
aP(:,nr)=aP(:,nr)+aN(:,nr)
sP(:,nr)=sP(:,nr)+2.0*aN(:,nr)*phi(1:n0,nr+1)
aN(:,nr)=0.0
!Sur
aP(:,1)=aP(:,1)+aS(:,1)
sP(:,1)=sP(:,1)+2.0*aS(:,1)*phi(1:n0,0)
aS(:,1)=0.0

CALL Gauss_TDMA2D(phi, n0, nr,aP,aE,aW,aN,aS,sP,n0,nr,max_iter,tolerance,residual)

CALL WriteScalarField2D('temperature', 1, phi, thetac, rc, n0, nr)

END PROGRAM

!Generación de la malla en cada dirección
SUBROUTINE MESH_1D(nx,x0,xl,x,xc)
IMPLICIT NONE
INTEGER nx,i
REAL*4 dx,x0,xl
REAL*4 x(0:nx),xc(0:nx+1)

dx=1.0/FLOAT(nx)

DO i=0,nx
   x(i)=x0+FLOAT(i)*(xl-x0)*dx
END DO

xc(0)=x(0); xc(nx+1)=x(nx)

DO i=1,nx
   xc(i)=0.5*(x(i)+x(i-1))
END DO

END SUBROUTINE


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine TDMA(x,a,b,c,d,n)

implicit none
  integer n,k
  real a(n),b(n),c(n),d(n),x(n),m

    do k=2,N
     m=a(k)/b(k-1)
     b(k)=b(k)-m*c(k-1)
     d(k)=d(k)-m*d(k-1)
    end do

    x(n)=d(n)/b(n)

    do k=n-1,1,-1
     x(k)=(d(k)-c(k)*x(k+1))/b(k)
    end do

end subroutine

!***************************************************************

subroutine lineX_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)

implicit none
integer bi,ei,bj,ej,i,j,nx,ny
real phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real a(ei),b(ei),c(ei),d(ei)
bi=1; bj=1
do j=bj,ej
do i=bi,ei
     a(i)=-aW(i,j)
     b(i)=aP(i,j)
     c(i)=-aE(i,j)
     d(i)=sp(i,j) + aN(i,j) * phi(i,j+1) + aS(i,j) * phi(i,j-1)
end do
  call TDMA(phi(bi:ei,j), a, b, c ,d ,ei)
end do

end subroutine

!****************************************************************

subroutine lineY_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)

implicit none
integer bi,ei,bj,ej,i,j,nx,ny
real phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real a(ej),b(ej),c(ej),d(ej)
bi=1; bj=1

do i=bi,ei
do j=bj,ej
     a(j)=-aS(i,j)
     b(j)=aP(i,j)
     c(j)=-aN(i,j)
     d(j)=sp(i,j) + aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j)
end do
  call TDMA(phi(i,bj:ej), a, b, c ,d ,ej)
end do

end subroutine

!****************************
!****************************

real*4 function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
implicit none
integer bi,ei,bj,ej,i,j,nx,ny
real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real*4 acum(ei,ej),residual,NINV
bi=1; bj=1
acum=0
NINV = 1.0 / dfloat(ei*ej)
do i=bi,ei
do j=bj,ej
acum(i,j) = aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j) + aN(i,j) * phi(i,j+1) + &
aS(i,j) * phi(i,j-1) + sp(i,j) - aP(i,j) * phi(i,j)
end do
end do
residual = sqrt( NINV * sum(acum * acum) )
calcResidual=residual
end function

!************************************************************************************************************

!Solución usando el método de Gauss TDMA
subroutine Gauss_TDMA2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

implicit none

integer bi,ei,bj,ej,i,j,nx,ny,count_iter,max_iter
real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real*4 residual,tolerance

     interface
         real*4 function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
         implicit none
         integer bi,ei,bj,ej,i,j,nx,ny
         real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
         end function
     end interface

count_iter=0; residual=1.0

do while((count_iter <= max_iter).and.(residual > tolerance))
call lineX_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
call lineY_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
residual = calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
count_iter=count_iter+1
end do

end subroutine


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.
!Para que escriba los resultados que obtenemos
Subroutine WriteScalarField2D(Name,kx,T,xc,yc,nx,ny)
integer i,j,nx,ny,kx
real*4 T(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
character*(*)Name
character*50 txt,Filename
write(txt,'(i6)')kx
txt=ADJUSTL(txt)
Filename=name//txt(1:len_trim(txt))//".txt"

open(10,file=Filename(1:len_trim(Filename)))
    do j=0,ny+1
    do i=0,nx+1
    write(10,*)yc(j)*cos(xc(i)),yc(j)*sin(xc(i)),T(i,j)
    !AQuí modificamos el como se escribe el scalar field para que esté en coordenadas polares y no se vea un cuadrado en el gnuplot, lo que hacemos es escribir 'x' y 'y' pero utilizando r y theta usando senos y cosenos. El porque se ve en la foto pero es por lo de co/hip y ca/hipo
    end do
    write(10,*)''
    end do
close(10)
End Subroutine

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine interpolateToNodesUs(uc,us,nx,ny)

implicit none
integer nx,ny
real:: us(0:nx,0:ny+1),uc(0:nx+1,0:ny+1)
integer bi,ei,bj,ej,i,j
bi=1; bj=1;ej=ny; ei=nx

! Internal points
do i=bi,ei
do j=bj-1,ej+1
    uc(I,J) = ( us(I-1,J) + us(I,J) ) * 0.5
end do
end do

uc(bi-1,:)=us(0,:)
uc(ei+1,:)=us(nx,:)

end subroutine

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine interpolateToNodesVs(vc,vs,nx,ny)
implicit none
integer nx,ny
real:: vs(0:nx+1,0:ny),vc(0:nx+1,0:ny+1)
integer bi,ei,bj,ej,i,j
bi=1; bj=1;ej=ny; ei=nx

do i=bi-1,ei+1
do j=bj,ej
    vc(I,J) = ( vs(I,J) + vs(I,J-1) ) * 0.5
end do
end do

vc(:,bj-1)=vs(:,0)
vc(:,ej+1)=vs(:,ny)
end subroutine

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Subroutine WriteVectorField(Name,kx,uc,vc,xc,yc,nx,ny)
integer i,j,nx,ny,kx
real*4 uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
character*(*)Name
character*50 txt,Filename
write(txt,'(i6)')kx
txt=ADJUSTL(txt)
Filename=name//txt(1:len_trim(txt))//".txt"

open(11,file=Filename(1:len_trim(Filename)))
    do i=0,nx+1
    do j=0,ny+1
    write(11,*)xc(i),yc(j),uc(i,j),vc(i,j)
    end do
    write(11,*)''
    end do
close(11)
End Subroutine
