PROGRAM adv_dif
IMPLICIT NONE
! COMENTARIO

!Declaración de variables
INTEGER nx, i, j, niter
REAL*8 dx, dv, L, T0, Tl, Se, Sw, Pe, gamma, U, err_targ, error
REAL*8, ALLOCATABLE :: T(:), aE(:), aW(:), aP(:), sP(:), Ta(:), x(:), T2(:)

!Inicializar variables
!Se le pone d0 a todos los que tienen decimales, en todo el codigo, para indicar que es una solución doble y no simple
nx=120; Se=1.0d0; Sw=1.0d0
T0=1.0d0; Tl=0.0d0
L=1.0d0
u=1.0d0
Pe=100.0d0

dx=L/dfloat(nx)
gamma=1.0d0/Pe
dv=dx

!Alojamiento de memoria dinámica
ALLOCATE(T(0:nx+1), aE(nx), aW(nx), aP(nx), sP(nx), Ta(0:nx+1), x(0:nx+1), T2(0:nx+1))

!Condición inicial porque siempre que resuelven numéricamente necesitamos un tiempo ficticio
T(1:nx)=0.5d0*(T(0)+T(nx+1))
Ta=0.0d0

!Condiciones de frontera
T(0)=T0; T(nx+1)=Tl

!Algoritmo de solución
!Calculo de Coeficientes
do i=1, nx
    !Esto era para Esquema Central
    !aE(i)=gamma*Se/dx-0.5*U*Se
    !aW(i)=gamma*Sw/dx+0.5*U*Sw

    !Esto es para el esquema UPWIND
    aE(i)=gamma*Se/dx-MIN(U*Se,0.0d0)
    aW(i)=gamma*Sw/dx+MAX(U*Sw,0.0d0)
    aP(i)=aE(i)+aW(i)
    Sp(i)=0.0d0
end do

!Corrección por condiciones de frontera
!Oeste
aP(1)=aP(1)+aW(1)
sP(1)=sP(1)+2.0d0*aW(1)*T(0)
aW(1)=0.0d0

!Este
aP(nx)=aP(nx)+aE(nx)
sP(nx)=sP(nx)+2.0d0*aE(nx)*T(nx+1)
aE(nx)=0.0d0


!Esto es para encontrar la convergencia de nuestro sistema y no estar probando diferentes iteraciones. Se calculan solas las iteraciones por lo que quitamos el número de iteraciones que ponemos manualmente
!Solución del sistema de ecuaciones por el método de Jacobi
err_targ=1.d-6 !Aquí fijamos el target o valor del error que queremos tener en nuestro problema
T2=0.0d0
niter=0.0d0
error=1.0d0

do while(error .gt. err_targ)
!gt es greater than. Entonces sería, si el error es más grande que el target error que queremos (err_targ) entonces entra al ciclo do
error=0.0d0
niter=niter+1
    do i=1, nx
        T(i)=(aE(i)*T(i+1)+aW(i)*T(i-1)+sP(i))/aP(i)
        error=error+abs(T(i)-T2(i))
    end do
    T2=T
end do
write(*,*)niter, error
!En la terminal va a escribir el número de iteraciones en las que convergieron las diferentes soluciones y la diferencia entre la solución actual y la vieja


! !Solución del sistema de ecuaciones por el método de Jacobi
! do j=1, 10000
!     do i=1, nx
!         t(i)=(aE(i)*T(i+1)+aW(i)*T(i-1)+sP(i))/aP(i)
!     end do
! end do

!Para crear nuestra malla
do i=1, nx
    x(i)=(dfloat(i)-0.5d0)*dx
end do

!Aquí ponemos nuestra solución análitica
do i=1, nx
    Ta(i)=(dEXP(Pe*U*x(i))-dEXP(Pe*U))/(1.0d0-dEXP(Pe*U))
end do


!Escritura de los resultados , ponemos 3 writes porque en ningun moment definimos una malla numérica entonces lo tenemos que hacer aquí
OPEN(1, FILE='TemperaturaUPWIND.txt', STATUS='REPLACE')
    Write(1,*)2.0d0*0.0d0,T(0), 1.0d0
do i=1, nx
    WRITE(1,*)x(i), T(i), Ta(i)
end do
    WRITE(1, *)L,T(nx+1),0.0d0
    CLOSE(1)

end program adv_dif
