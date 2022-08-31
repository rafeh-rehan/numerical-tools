module interpolate

implicit double precision (a-h,o-z)

contains

! 1D Barycentric Interpolation. Reference:
! https://terpconnect.umd.edu/~petersd/666/BarycentricLagrange1.pdf
subroutine bary(n, xs, ys, x, y)
double precision, intent(in)  :: xs(n), ys(n)
double precision, intent(out) :: y

double precision :: ws(n)

eps = epsilon(1.0d0)


!size(xs) = n
if (size(xs) .ne. size(ys)) then
print *, "error: x array must be same size as y array"
stop
end if


! Return inputted node points
do i = 1, n
if (abs(x - xs(i)) .le. eps) then 
y = ys(i)
return
end if
end do

k  = n-1 ! degree of polynomial
pi = acos(-1.0d0)


! Compute weights

do i = 1, n 
ws(i) = 1.0d0

do j = 1, i-1 
ws(j) = ws(j)/(xs(j) - xs(i))
ws(i) = ws(i)/(xs(i) - xs(j))
end do
end do

do i =1,n
!print*,ws(i)
end do

! Interpolate using Barycentric Formula
s  = 0.0
ss = 0.0 

do i = 1,n
s  = s  + ws(i)*ys(i) / (x-xs(i)) 
ss = ss + ws(i) / (x-xs(i))
end do

y = s/ss

end subroutine




! 2D Billinear Interpolation. Reference: 
! https://en.wikipedia.org/wiki/Bilinear_interpolation
subroutine bili(xs, len_x, ys, len_y, zs, x, y, z)
double precision              :: xs(len_x), ys(len_y), zs(len_x,len_y)
double precision, intent(out) :: z

!print*, x,y

eps = epsilon(1.0)

!! Find closest nodes to inputted data
n1 = 0 ! x1
n2 = 0 ! x2

m1 = 0 ! y1
m2 = 0 ! y2

i  = 1
ii = 1


! Find nearest x nodes
do while (i .le. len_x-1)

if ( (xs(i) .le. x) .and. (x .le. xs(i+1)) ) then
n1 = i
n2 = i+1
exit
end if

i = i+1
if (i .eq. len_x) then 
print*, 'error: x is not in the data set', x
stop
endif 

end do


! Find nearest y nodes
do while (ii .le. len_y-1)

if ( (ys(ii) .le. y) .and. (y .le. ys(ii+1)) ) then
m1 = ii
m2 = ii+1
exit
end if

ii = ii+1
if (ii .eq. len_y) then 
print*, 'error: y is not in the data set', y
stop
end if

end do


! Interpolate
x1 = xs(n1)
x2 = xs(n2)

y1 = ys(m1)
y2 = ys(m2)

d = (x2 - x1)*(y2 - y1)
z = 1/d * (zs(n1,m1)*(x2-x)*(y2-y)+zs(n2,m1)*(x-x1)*(y2-y) &
    + zs(n1,m2)*(x2-x)*(y-y1) + zs(n2,m2)*(x-x1)*(y-y1))

end subroutine


end module
