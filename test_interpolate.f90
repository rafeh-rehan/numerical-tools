program test_interpolate

use interpolate

implicit double precision (a-h,o-z)

double precision, allocatable :: xs(:), ys(:), zs(:,:)

pi = acos(-1.0)
n  = 100
allocate(xs(n),ys(n))

a = 50

!!!
! 
! 1D Interpolation (Barycentric)
! 
!!!

do i = 1, n
x = cos( (2*i-1)*pi / (2*n) ) * a ! chebyshev nodes on
                                  ! interval [a,-a]
!print*,x
xs(i) = -x
ys(i) = f(xs(i))
end do


do k = 1, 10**4
call random_number(dd)
z = -a + dd*(2*a)

call bary(n, xs, ys, z, y)
!print*,z,y,f(z)

rel   = abs(y - f(z))/abs(f(z))
error = max(error,rel)
end do

print*, 'Maximum relative error 1D barycentric interpolation:', error

!!!
! 
! 2D Interpolation (Bilinear)
! 
!!!

allocate(zs(n,n))

do i = 1, n
y = cos( (2*i-1)*pi / (2*n) ) * a
ys(i) = -y

do j = 1, n
zs(j,i) = g(xs(j),ys(i))
end do
end do

do k = 1, 15
call random_number(dd1)
call random_number(dd2)
x = -a + dd1 * (2*a)
y = -a + dd2 * (2*a)

!print*,x,y
!stop 

call bili(xs, n, ys, n, zs, x, y, z)

rel2   = abs(z - g(x,y))/abs(g(x,y))
error2 = max(error2,rel2)

end do

print*, 'Maximum relative error 2D bilinear interpolation:', error2


contains

function f(x)

f = 17*(x/3)**3 + 2.67
end function


function g(x,y)

g = 12*xy + 3*x + y
end function

end program
