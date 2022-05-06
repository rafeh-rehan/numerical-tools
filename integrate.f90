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



end module
