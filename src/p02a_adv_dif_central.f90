PROGRAM adv_dif
IMPLICIT NONE
! COMENTARIO

!Declaración de variables
INTEGER nx, i, j
REAL*4 dx, dv, L, T0, Tl, Se, Sw, Pe, gamma, U
REAL*4, ALLOCATABLE :: T(:), aE(:), aW(:), aP(:), sP(:)

!Inicializar variables
nx=60; Se=1.0; Sw=1.0
T0=1.0; Tl=0.0
L=1.0
u=1.0
Pe=100.0

dx=L/float(nx)
gamma=1.0/Pe
dv=dx

!Alojamiento de memoria dinámica
ALLOCATE(T(0:nx+1), aE(nx), aW(nx), aP(nx), sP(nx))

!Condición inicial porque siempre que resuleves numericamente necesitamos un tiempo ficticio
T(1:nx)=0.5*(T(0)+T(nx+1))

!Condiciones de frontera
T(0)=T0; T(nx+1)=Tl

!Algoritmo de solución
!Calculo de Coeficientes
do i=1, nx
    aE(i)=gamma*Se/dx-0.5*U*Se
    aW(i)=gamma*Sw/dx+0.5*U*Sw
    aP(i)=aE(i)+aW(i)
    Sp(i)=0.0
end do

!Corrección por condiciones de frontera
!Oeste
aP(1)=aP(1)+aW(1)
sP(1)=sP(1)+2.0*aW(1)*T(0)
aW(1)=0.0

!Este
aP(nx)=aP(nx)+aE(nx)
sP(nx)=sP(nx)+2.0*aE(nx)*T(nx+1)
aE(nx)=0.0

!Solución del sistema de ecuaciones por el método de Jacobi
do j=1, 1000
    do i=1, nx
        t(i)=(aE(i)*T(i+1)+aW(i)*T(i-1)+sP(i))/aP(i)
    end do
end do


!Escritura de los resultados , ponemos 3 writes porque en ningun moment definimos una malla numérica entonces lo tenemos que hacer aquí
OPEN(1, FILE='Temperatura.txt', STATUS='REPLACE')
Write(1,*)2.0*0.0,T(0)
do i=1, nx
    WRITE(1,*)(float(i)-0.5)*dx, T(i)
end do
    WRITE(1, *)L,T(nx+1)
    CLOSE(1)

end program adv_dif
