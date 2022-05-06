module integrate

use tools


implicit double precision (a-h,o-z)


contains

! Reference : https://www.freecodecamp.org/news/simpsons-rule/
subroutine simp(f,a,b,n,y)
!
!  Computes the definite integral of a function
!  using Simpson's 3/8 rule.
!
!  f   -  function to be integrated
!  a,b -  domain of integration [a,b]
!  n   -  number of subintervals
!
!  y   -  value of the integral
!

double precision               :: x(n+1)
double precision, intent(out)  :: y

pi = acos(-1.0d0)

! check if number of subintervals is correct
if (mod(n,3) .ne. 0) then
print*, 'n must be a multiple of 3'
stop
end if

s  = 0.0d0
h  = (b-a)/n


x(1) = a
do i = 2,n+1
x(i) = a + h*i
end do

do i = 1, n/3
s = s + f(x(3*i-3+1)) + 3.0d0*f(x(3*i-2+1)) + 3.0d0*f(x(3*i-1+1)) + f(x(3*i+1))
end do

y = 3d0*h/8d0 * s

end subroutine


subroutine simp2(f,t,a,b,n,y)
!
!  Computes the definite integral of a function
!  using Simpson's 3/8 rule.
!
!  f   -  function to be integrated
!  t   -  additional parameter for f
!  a,b -  domain of integration [a,b]
!  n   -  number of subintervals
!
!  y   -  value of the integral
!

double precision               :: x(n+1)
double precision, intent(out)  :: y

pi = acos(-1.0d0)

! check if number of subintervals is correct
if (mod(n,3) .ne. 0) then
print*, 'n must be a multiple of 3'
stop
end if

s  = 0.0
h  = (b-a)/n

x(1) = a
do i = 2,n+1
x(i) = a + h*i
end do

do i = 1, n/3
s = s + f(x(3*i-3+1),t) + 3.0d0*f(x(3*i-2+1),t) + 3.0d0*f(x(3*i-1+1),t) + f(x(3*i+1),t)
end do

y = s*3d0*h/8.0d0

end subroutine


subroutine simp_lum(f,t,t2,a,b,n,y)
!
!  Computes the 1D definite integral of a function
!  using Simpson's 3/8 rule. For LUMINOSITY ONLY - accretion
!
!  f   -  function to be integrated
!  t   -  array, dependent on same var as f (chebyshev nodes)
!  t2   -  other parameters for f
!  a,b -  domain of integration [a,b]
!  n   -  number of subintervals
!
!  y   -  value of the integral
!

double precision, intent(out)  :: y
double precision               :: t2(n+1)
double precision               :: x(n+1)

pi = acos(-1.0d0)

! check if number of subintervals is correct
if (mod(n,3) .ne. 0) then
print*, 'n must be a multiple of 3'
stop
end if

sums = 0.0d0
h    = (b-a)/n

x(1) = a
do i = 2,n+1
x(i) = a + h*i
!print*,a+h*i
end do

do i = 1, n/3

sums = sums + f(x(3*i-3+1),t,t2(3*i-3+1)) + 3.0d0*f(x(3*i-2+1),t,t2(3*i-2+1)) &
       + 3.0d0*f(x(3*i-1+1),t,t2(3*i-1+1)) + f(x(3*i+1),t,t2(3*i+1))
end do

y = sums*3d0*h/8.0d0

end subroutine




! Gaussian Quadrature with n = 30
subroutine gauss(f,a,b,n_ints,f_integral)

!
!  Caluclates the integral
!
!  \int_a^b f(x)dx ~ \Sum^n_{i=1} w_i f(x_i)
!
!  By performing Legendre-Gauss quadrature.
!  Table of abscissa (x_i) and weigths (w_i)
!  are from:
!  https://pomax.github.io/bezierinfo/legendre-gauss.html#n30
!
!  n_ints: number of subintervals to break up the integral on,
!          each subinterval uses 30th order Gauss-Legendre
!          quadrature.


! Order of quadrature rule, array of abscissa and weights
integer, parameter :: n=30
double precision   :: xs(n),ws(n)
double precision   :: sub_ints(n_ints)

ws  = [0.1028526528935588d0, 0.1028526528935588d0, &
       0.1017623897484055d0, 0.1017623897484055d0, &
       0.0995934205867953d0, 0.0995934205867953d0, &
       0.0963687371746443d0, 0.0963687371746443d0, &
       0.0921225222377861d0, 0.0921225222377861d0, &
       0.0868997872010830d0, 0.0868997872010830d0, &
       0.0807558952294202d0, 0.0807558952294202d0, &
       0.0737559747377052d0, 0.0737559747377052d0, &
       0.0659742298821805d0, 0.0659742298821805d0, &
       0.0574931562176191d0, 0.0574931562176191d0, &
       0.0484026728305941d0, 0.0484026728305941d0, &
       0.0387991925696271d0, 0.0387991925696271d0, &
       0.0287847078833234d0, 0.0287847078833234d0, &
       0.0184664683110910d0, 0.0184664683110910d0, &
       0.0079681924961666d0, 0.0079681924961666d0]

xs  = [-0.0514718425553177d0, 0.0514718425553177d0, &
       -0.1538699136085835d0, 0.1538699136085835d0, &
       -0.2546369261678899d0, 0.2546369261678899d0, &
       -0.3527047255308781d0, 0.3527047255308781d0, &
       -0.4470337695380892d0, 0.4470337695380892d0, &
       -0.5366241481420199d0, 0.5366241481420199d0, &
       -0.6205261829892429d0, 0.6205261829892429d0, &
       -0.6978504947933158d0, 0.6978504947933158d0, &
       -0.7677774321048262d0, 0.7677774321048262d0, &
       -0.8295657623827684d0, 0.8295657623827684d0, &
       -0.8825605357920527d0, 0.8825605357920527d0, &
       -0.9262000474292743d0, 0.9262000474292743d0, &
       -0.9600218649683075d0, 0.9600218649683075d0, &
       -0.9836681232797472d0, 0.9836681232797472d0, &
       -0.9968934840746495d0, 0.9968934840746495d0]

! Create sub intervals
do i = 1,n_ints
zi = i
call remap(1d0,dble(n_ints),a,b,zi,z)
sub_ints(i) = z
end do


s = 0d0

!do i = 1,n_ints-1
do j = 1,n
!aa = sub_ints(i)
!bb = sub_ints(i+1)

!x = (bb-aa)/2d0*xs(j)+(bb+aa)/2d0
x = (b-a)/2d0*xs(j)+(b+a)/2d0
s =+ ws(j)*f(x)
end do
!end do

f_integral = (b-a)/2d0*s


end subroutine





end module
