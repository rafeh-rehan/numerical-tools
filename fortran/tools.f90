module tools

implicit double precision (a-h,o-z)


contains


subroutine remap(a, b, c, d, x, y)

!
!   Maps an element in the interval [a,b]
!   on to the interval [c,d]
!

double precision, intent(out)  :: y

y = c + (d-c)/(b-a) * (x-a)

end subroutine

subroutine cart2sph(x, y, z, r, theta, phi)

!
!   Returns Spherical coordinates given
!   Cartesian Coordinates
!

double precision, intent(out) :: r, theta, phi
pi = acos(-1.0)

r    = sqrt(x**2 + y**2 + z**2)

theta = acos(z/r)

if (x .gt. 0.0) then
phi = atan(y/x)
else if (x .lt. 0.0) then
phi = atan(y/x) + pi
else
phi = pi/2.0
end if


end subroutine


subroutine linspace(a, b, n, ab)

!
! Creates an equispaced interval from [a,b]
! of size n+1
!

double precision, intent(out)  :: ab(n+1)

! Step size
h = (b-a)/n

do i = 1,n+1
ab(i) = a + (i-1)*h
end do


end subroutine


end module
