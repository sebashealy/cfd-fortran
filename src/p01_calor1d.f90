PROGRAM Calor1D
IMPLICIT NONE
!Con IMPLICIT NONE borras todas las variables que ya podría haber

!Declaración de variables
INTEGER nx, i, j
REAL*4 Ta, Tb, L, q, k , Se, Sw, dx, dv
REAL*4, ALLOCATABLE :: T(:), Te(:), aE(:), aW(:), sP(:), x(:), ap(:)

!Inicializar variables
nx=30
Ta=100.0
Tb=200.0
k=0.5
q=1.0e6
L=0.02

Se=1.0
Sw=1.0
dx=L/float(nx)
dv=dx

!Alojamiento de memoria dinámica
ALLOCATE(T(0:nx+1), Te(0:nx+1), aP(nx), aE(nx), aW(nx), sP(nx), x(0:nx+1))

!Condición incial
T=0.0

!Condiciones de frontera
T(0)=Ta; T(nx+1)=Tb

!Algoritmo de solución
!Cálculo de los coeficientes
do i=1, nx
aE(i)=k*Se/dx
aW(i)=k*Sw/dx
aP(i)=aE(i)+aW(i)
sP(i)=q*dv
end do

!Corrección por condiciones de frontera
!oeste
aP(1)=aP(1)+aW(1)
sP(1)=sP(1)+2.0*aW(1)*T(0)
aW(1)=0.0

!este
aP(nx)=aP(nx)+aE(nx)
sP(nx)=sP(nx)+2.0*aE(nx)*T(nx+1)
aE(nx)=0.0

!Solución del sistema de ecuaciones por el método de Jacobi
do j=1,1000
do i =1,nx
      T(i)=(aE(i)*T(i+1)+aW(i)*T(i-1)+sP(i))/aP(i)
end do
end do

x(0)=0.0
do i=1, nx
x(i)=(float(i)-0.5)*dx
end do
x(nx+1)=L

!Solución analítica (De la presentación)
Te(0)=Ta
do i=1, nx
Te(i)=Ta+((Tb-Ta)/L)*x(i)-(q/(2.*k))*x(i)*(x(i)-L)
end do
Te(nx+1)=Tb

!Para escribir los resultados
OPEN(1, FILE='Temperatura.txt', STATUS='REPLACE')
DO i=0, nx+1
      WRITE(1, *)x(i), T(i), Te(i)
END DO

CLOSE(1)

end program Calor1D
