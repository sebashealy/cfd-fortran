PROGRAM Flujo_canal_obs_sol !Método Simplec
IMPLICIT NONE

!Declaración de variables
CHARACTER*50 Filename, itchar
CHARACTER*50 name
INTEGER i, j, nx, ny, maxit_SIMPLEC, itmax, ei, ej, it, max_iter, c
!el maxiter_S es para el método SIMPLEC (máximo de iteraciones), itmax esun contador del cíclo temporal, maxiter son las iteraciones para el método de solución TDMA
REAL*4 dx, dy, dv, x0, xl, y0, yl, Se, Sw, Sn, Ss, tolerance, residual, time, dt, ue, uw, vn, vs, Div, Re, &
    gamma, dbkx, dbky
!Div=Divergencia (criterio de diver o convergencia), dbkx y dbky son variables auxiliares para la solución del método
REAL*4, ALLOCATABLE, DIMENSION(:,:) :: aE, aW, aN, aS, aP, sP, u, u1, v, v1, P, Pp, d_h, d_v, uc, vc
!u y u1 son las velocidades actuales y la anterior, lo mismo con v. Pp es la presión corregida, de y dn son los coeficientes que se van a usar para conseguir las correcciones de la presión, P es la presión, la d_h y d_v son de y dn pero con otros nombres para seguir la solución de la presentación
REAL*4, ALLOCATABLE, DIMENSION(:) :: x, xc, y, yc
REAL*4 x0_obs, y0_obs, xl_obs, yl_obs, v_obs, u_obs
INTEGER, ALLOCATABLE :: mark_cells(:,:)

!inicialización de variables
nx=400; ny=40 !Número de volumenes de control
x0=0.0; xl=20.0 !Definir el dominio computacional en x
y0=-2.5; yl=2.5 !Definir el dominio computacional en y
Re=100.0

max_iter=1000
maxit_SIMPLEC=100
tolerance=1e-5
time=50.0 !Tiempo total de simulación
dt=0.001 !Paso temporal
it=0

itmax=int(time/dt)+1 !Número total de iteraciones

dx=(xl-x0)/float(nx)
dy=(yl-y0)/float(ny)

Se=dy; Sw=dy !El tamaño de las caras de los volúmenes de control
Sn=dx; Ss=dx

dv=dx*dy
gamma=1.0/Re

!Alojamiento de la memoria dinámica, le ponemos tamaño a los arreglos
ALLOCATE(aP(nx,ny), aE(nx,ny), aW(nx,ny), aN(nx,ny), aS(nx,ny), sP(nx,ny), uc(0:nx+1,0:ny+1), &
    vc(0:nx+1,0:ny+1))
!Los c's son para tener los valores en el centro de control
ALLOCATE(x(0:nx),xc(0:nx+1),y(0:ny),yc(0:ny+1),P(0:nx+1,0:ny+1),Pp(0:nx+1,0:ny+1))
!x solo toma las caras y xc toma los centros pero como los dos empiezan en el borde, xc tiene más valores. P es la presión normal y Pp la presión corregida o prima
ALLOCATE(u(0:nx,0:ny+1),u1(0:nx,0:ny+1),v(0:nx+1,0:ny), v1(0:nx+1,0:ny), d_h(0:nx+1,0:ny+1), &
    d_v(0:nx+1,0:ny+1))
!Se quita uno en x para las u y en y para las v por el desfasamiento de la malla. u es la actual y u1 es el valor anterior
ALLOCATE(mark_cells(0:nx+1,0:ny+1))

!Creación de malla
Call MESH_1D(nx,x0,xl,x,xc)
Call MESH_1D(ny,y0,yl,y,yc)

!Condiciones de fronteras e iniciales
u=1.0; v=0.0; P=0.0; Pp=0.0 !Iniciales, se pone u=1.0 para que no empiece en estado estacionario y ayudarle a la máquina.
d_h=0.0; d_v=0.0 !iniciales

!Flujo uniforme cara Oeste
u(0,:)=1.0 !Velocidad de entrada, los : son para decir que son en todo el dominio en 'y' de la variable u que es de 0:ny+1
!Cara norte y sur
u(:,0)=0.0; u(:,ny+1)=0.0

u1=u; v1=v

x0_obs=4.0; xl_obs=5.0
y0_obs=-0.5; yl_obs=0.5 !Localización del obstáculo sólido
mark_cells=0

OPEN(3, FILE='obstaculo.txt', STATUS='replace')
!Ciclo para marcar el obstáculo y anular velocidades dentro de él
Do i=1, nx
     Do j=1, ny
         if((xc(i) .gt. x0_obs .and. xc(i) .lt. xl_obs) .and. (yc(j) .gt. y0_obs) .and. (yc(j) .lt. yl_obs)) then
             mark_cells(i,j)=1; u(i,j)=0.0; v(i,j)=0.0;
             WRITE(3,*) xc(i), yc(j)
         end if
     end do
end do
CLOSE(3)

do i=1, nx
     do j=1, ny
         if(mark_cells(i,j) .eq. 0 .and. mark_cells(i+1,j) .eq. 1) then
             u(i,j)=0.0;
         end if
         if(mark_cells(i,j) .eq. 0 .and. mark_cells(i,j+1) .eq. 1) then
             v(i,j)=0.0;
         end if
     end do
end do


!Creación de una archivo de animación
OPEN(2, FILE='anim.gnp', STATUS='REPLACE')
WRITE(2,*)'set xrange [0:20]; set yrange[-2.5:2.5]; set size ratio 0.25'
WRITE(2,*)"set xlabel 'x'; set ylabel 'y' rotate by 0; set title 'Lid-Driven Re=100'"

!Inicia ciclo temporal
Do it=1, itmax

Div=1.0
!El Div lo inciamos en 1.0 para que entre el ciclo SIMPLEC aunque en verdad queremos que sea cercano a 0. La divergencia nos indica que tan bien está el programa a través de la convergecia
c=1

!Inicia Ciclo SIMPLEC
Do WHILE ((c .LT. maxit_SIMPLEC) .AND. (Div .GT. tolerance))

!Algoritmo de solución
!Cálculo de los coeficientes

!ECUACIÓN COMPONENTE U ########################### Vamos a resolver Navier-Stokes para u
aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; sP=0.0 !Se inicializan para que los coeficientes no tengan otros valores

ei=UBOUND(u,1)-1
!UPPER BOUND nos da el valor de la cara +1 o frontera que es i+1 pero le restamos 1 para que esté en el dominio computacional entonces es i, la nx
ej=UBOUND(u,2)-1 !UPPER BOUND, la ny. el 2 es porque es para el segundo índice, es decir j

!Cálculo de coeficientes
Do i=1, ei
    Do j=1, ej
    !Flujos en las caras
    ue=(u(i,j)+u(i+1,j))*0.5 !De la presentación
    uw=(u(i,j)+u(i-1,j))*0.5
    vn=(v(i,j)+v(i+1,j))*0.5
    vs=(v(i,j-1)+v(i+1,j-1))*0.5

    aE(i,j)=gamma*Se/dx-MIN(ue*Se, 0.0)
    !UPWIND lo que hace es darle más peso a los volúmenes de control que están del lado contrario al que el flujo está muy direccionado. Es decir, más peso flujo -> o <- flujo, más peso
    aW(i,j)=gamma*Sw/dx+MAX(uw*Sw, 0.0)
    aN(i,j)=gamma*Sn/dy-MIN(vn*Sn, 0.0)
    aS(i,j)=gamma*Ss/dy+MAX(vs*Ss, 0.0)

    aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)+dv/dt
    sP(i,j)=u(i,j)*dv/dt-(P(i+1,j)-P(i,j))*dv/dx

!Corrección de frontera inmersa: del obstáculo
!Este
if((mark_cells(i,j) .eq. 0) .and. (mark_cells(i+2,j) .eq. 1)) then
    sP(i,j)=sP(i,j)+aE(i,j)*u1(i+1,j)
    aE(i,j)=0.0
end if
!Oeste
if((mark_cells(i,j) .eq. 0) .and. (mark_cells(i-1,j) .eq. 1)) then
    sP(i,j)=sp(i,j)+aW(i,j)*u1(i-1,j)
    aW(i,j)=0.0
end if
!Norte
if((mark_cells(i,j) .eq. 0) .and. (mark_cells(i+1,j+1) .eq. 1) .or. mark_cells(i,j+1) .eq. 1) then
    aP(i,j)=aP(i,j)+aN(i,j)
    sP(i,j)=sP(i,j)+2.0*aN(i,j)*u1(i,j+1)
    aN(i,j)=0.0
end if
!Sur
if((mark_cells(i,j) .eq. 0) .and. (mark_cells(i+1,j-1) .eq. 1) .or. mark_cells(i,j-1) .eq. 1) then
    aP(i,j)=aP(i,j)+aS(i,j)
    sP(i,j)=sP(i,j)+2.0*aS(i,j)*u1(i,j-1)
    aS(i,j)=0.0
end if
!Obstáculo
if(mark_cells(i,j) .eq. 1 .or. (mark_cells(i,j) .eq. 0 .and. mark_cells(i+1,j) .eq. 1)) then
    aW(i,j)=0.0; aE(i,j)=0.0; aS(i,j)=0.0; aN(i,j)=0.0; u1(i,j)=0.0
    aP(i,j)=dv/dt ! también se puede hacer ap(i,j)=1.0 con sp(i,j)=u1(i,j)
    sP(i,j)=u1(i,j)*dv/dt
end if

end do
end do

!Corrección por condiciones de fronetera
dbkx=0.0; dbky=1.0
!dbky es 1 porque la malla no se desfasa en el eje y por lo que solo se modifican las condiciones de frontera para el eje x (dbkx=0). Se hace así para tener una solución generalizada que después vamos a poder reciclar

!Arreglo condiciones de frontera para malla Staggered
!Este
aP(ei,1:ej)=aP(ei,1:ej)+dbkx*aE(ei,1:ej)
!Se usa ei y ej porque los coeficientes que se están resolviendo estan dentro del dominio de control y ei y ej son los últimos volumenes de control en x y y respectivamente
sP(ei,1:ej)=sP(ei,1:ej)+aE(ei,1:ej)*0.0
aE(ei,1:ej)=0.0
!Oeste
aP(1,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
sP(1,1:ej)=sP(1,1:ej)+(1.0+dbkx)*aW(1,1:ej)*u(0,1:ej) !Aquí se observa el uso de dbkx porque s
aW(1,1:ej)=0.0
!Norte
aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
sP(1:ei,ej)=sP(1:ei,ej)+(1.0+dbky)*aN(1:ei,ej)*u(1:ei,ej+1)
aN(1:ei,ej)=0.0
!Sur
aP(1:ei,1)=aP(1:ei,1)+dbky*aS(1:ei,1)
sP(1:ei,1)=sP(1:ei,1)+(1.0+dbky)*aS(1:ei,1)*u(1:ei,0)
aS(1:ei,1)=0.0

!Cálculo de coeficientes para corrección de la presión de=Ae/Ap-sigma(a_nb)
d_h(1:ei,1:ej)=Se/(aP(1:ei,1:ej)-(aE(1:ei,1:ej)+aW(1:ei,1:ej)+aN(1:ei,1:ej)+aS(1:ei,1:ej)))

!Solución del sistema de ecuaciones
Call Gauss_TDMA2D(u1,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

!=================================================
!Ecuación componente V #####################################
aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; sP=0.0 !Se hacen cero para que no traigan los valores calculados anteriormente

ei=UBOUND(v,1)-1 !UPPER BOUND nos da el valor de la cara +1 como i+1 pero le restamos 1 entonces es i, la nx
ej=UBOUND(v,2)-1 !UPPER BOUND, la ny. el 2 es porque es para el segundo índice

!Cálculo de coeficientes
Do i=1, ei
    Do j=1, ej
    !Flujos en las caras
    ue=(u(i,j)+u(i,j+1))*0.5 !De la presentación
    uw=(u(i-1,j)+u(i-1,j+1))*0.5
    vn=(v(i,j)+v(i,j+1))*0.5
    vs=(v(i,j)+v(i,j-1))*0.5

    aE(i,j)=gamma*Se/dx-MIN(ue*Se,0.0)
    aW(i,j)=gamma*Sw/dx+MAX(uw*Sw,0.0)
    aN(i,j)=gamma*Sn/dy-MIN(vn*Sn,0.0)
    aS(i,j)=gamma*Ss/dy+MAX(vs*Ss,0.0)
    aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)+dv/dt
    sP(i,j)=v(i,j)*dv/dt-(P(i,j+1)-P(i,j))*dv/dy

!Corrección de frontera inmersa: del obstáculo
!Este
if((mark_cells(i,j) .eq. 0) .and. (mark_cells(i+1,j+1) .eq. 1) .or. mark_cells(i+1,j) .eq. 1) then
    aP(i,j)=aP(i,j)+aE(i,j)
    sP(i,j)=sP(i,j)+2.0*aE(i,j)*v1(i+1,j)
    aE(i,j)=0.0
end if
!Oeste
if((mark_cells(i,j) .eq. 0) .and. (mark_cells(i-1,j+1) .eq. 1) .or. mark_cells(i-1,j) .eq. 1) then
    aP(i,j)=aP(i,j)+aW(i,j)
    sP(i,j)=sP(i,j)+2.0*aW(i,j)*v1(i-1,j)
    aW(i,j)=0.0
end if
!Norte
if((mark_cells(i,j) .eq. 0) .and. (mark_cells(i,j+2) .eq. 1)) then
    sP(i,j)=sP(i,j)+aN(i,j)*v1(i,j+1)
    aN(i,j)=0.0
end if
!Sur
if((mark_cells(i,j) .eq. 0) .and. (mark_cells(i,j-1) .eq. 1)) then
    sP(i,j)=sP(i,j)+aS(i,j)*v1(i,j-1)
    aS(i,j)=0.0
end if
!Obstáculo
if(mark_cells(i,j) .eq. 1 .or. (mark_cells(i,j) .eq. 0 .and. mark_cells(i,j+1) .eq. 1)) then
    aW(i,j)=0.0; aE(i,j)=0.0; aS(i,j)=0.0; aN(i,j)=0.0; v1(i,j)=0.0
    aP(i,j)=dv/dt ! también se puede hacer ap(i,j)=1.0 con sp(i,j)=u1(i,j)
    sP(i,j)=v1(i,j)*dv/dt
end if

end do
end do

!Corrección por condiciones de fronetera
dbkx=1.0; dbky=0.0
!Ahora dbkx es 1 porque la malla no se desfasa en el eje x por lo que solo se modifican las condiciones de frontera para el eje y (dbky=0). Se hace así para tener una solución generalizada que después vamos a poder reciclar

!Arreglo condiciones de frontera para malla Staggered
!Este Aquí no se corrigió igual que en el caso de u porque el cambio es en x, no en y
aP(ei,1:ej)=aP(ei,1:ej)-dbkx*aE(ei,1:ej)
sP(ei,1:ej)=sP(ei,1:ej)+(dbkx)*aE(ei,1:ej)*0.0
aE(ei,1:ej)=0.0
!Oeste
aP(1,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
sP(1,1:ej)=sP(1,1:ej)+(1.0+dbkx)*aW(1,1:ej)*v(0,1:ej)
aW(1,1:ej)=0.0
!Norte
aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
sP(1:ei,ej)=sP(1:ei,ej)+(1.0+dbky)*aN(1:ei,ej)*v(1:ei,ej+1)
aN(1:ei,ej)=0.0
!Sur
aP(1:ei,1)=aP(1:ei,1)+dbky*aS(1:ei,1)
sP(1:ei,1)=sP(1:ei,1)+(1.0+dbky)*aS(1:ei,1)*v(1:ei,0)
aS(1:ei,1)=0.0

!Cálculo de coeficientes para corrección de la presión dn=An/Ap-sigma(a_nb)
d_v(1:ei,1:ej)=Sn/(aP(1:ei,1:ej)-(aE(1:ei,1:ej)+aW(1:ei,1:ej)+aN(1:ei,1:ej)+aS(1:ei,1:ej)))

!Solución del sistema de ecuaciones
Call Gauss_TDMA2D(v1,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

!=================================================
!Correcciones de la presión #####################################
aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; sP=0.0; Pp=0.0

ei=UBOUND(Pp,1)-1 !UPPER BOUND nos da el valor de la cara +1 como i+1 pero le restamos 1 entonces es i, la nx
ej=UBOUND(Pp,2)-1 !UPPER BOUND, la ny. el 2 es porque es para el segundo índice

!Cálculo de coeficientes
Do i=1, ei
    Do j=1, ej
        !Flujos en las caras obteniendolos de la gráfica de la diapositiva 92
        ue=u1(i,j)
        uw=u1(i-1,j)
        vn=v1(i,j)
        vs=v1(i,j-1)
        !Cálculo de los coeficientes a partir de las fórmulas obtenidas de la diapositiva 99
         aE(i,j)=d_h(i,j)*Se
         aW(i,j)=d_h(i-1,j)*Sw
         aN(i,j)=d_v(i,j)*Sn
         aS(i,j)=d_v(i,j-1)*Ss
         aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)
         sP(i,j)=-ue*Se+uw*Sw-vn*Sn+vs*Ss !Divergencia de u
    end do
end do

!Solución del sistema de ecuaciones para obtener la presión corregida (Pp=Presión prima)
Call Gauss_TDMA2D(Pp,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

!Correggimos la presión
P=P+Pp

!Corrección de velocidades, u y v
!COMPONENTE U. Se ponen así la nx y ny por el desfase que tiene la malla en el eje x
do i=1,nx-1
    do j=1,ny
         if(mark_cells(i,j) .eq. 0 .and. mark_cells(i+1,j) .ne. 1)then
             u1(i,j)=u1(i,j)+d_h(i,j)*(pp(i,j)-pp(i+1,j))
         end if
    end do
end do

!COMPONENTE V
do i=1,nx
    do j=1,ny-1
         if(mark_cells(i,j) .eq. 0 .and. mark_cells(i,j+1) .ne. 1)then
             v1(i,j)=v1(i,j)+d_v(i,j)*(pp(i,j)-pp(i,j+1))
         end if
    end do
end do

!Calcular el valor de la velocidad en la frontera Este por ser de tipo Neumann
u1(nx,:)=u1(nx-1,:)
v1(nx+1,:)=v1(nx,:)

!Verficiar la Div de U para cumnplir con el criterio del paso 6
Div=MAXVAL(ABS(sP(1:ei,1:ej))) !Divergencia de u=0, ec. de continuidad

c=c+1

END DO !Aquí termina el ciclo SIMPLEC

u=u1; v=v1
write(*,*)it,c, Div

!Escritra de resultados para la animación
IF (MOD(it,100) .EQ. 0) THEN !Lo que esto está haciendo es llamar a las subrutinas de abajo cada 100 iteraciones

    CALL interpolateToNodesUs(uc,u,nx,ny)
    CALL interpolateToNodesVs(vc,v,nx,ny)
    CALL WriteVectorField('vel',it,uc,vc,xc,yc,nx,ny)

    WRITE(itchar, '(i6)')it
    itchar=ADJUSTL(itchar) !El ADJUSTL es para acomodar la palabra todo a la izquierda y que no haya espacios vacios
    itchar=itchar(1:LEN_TRIM(itchar))
    !LEN_TRIM es para que corte el tamaño del character a sólo el número de characters que estamos usando en verdad

    WRITE(2,*)"p 'vel" //itchar(1:LEN_TRIM(itchar))//".txt' u 1:2:(0.23*3):(0.25*4) w vec lc 6"
   ! WRITE(2,*) 'pause 0.05'

end if

END DO !Términa el ciclo temporal

CALL interpolateToNodesUs(uc,u,nx,ny)
CALL interpolateToNodesVs(vc,v,nx,ny)
CALL WriteVectorField('vel',0,uc,vc,xc,yc,nx,ny)

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
write(10,*)xc(i),yc(j),T(i,j)
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
