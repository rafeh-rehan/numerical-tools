module roots

use tools


implicit real (a-h,o-z)


contains


subroutine newt(f,df,x0,x,ibug)
real, intent(out) :: x 
! 
!  Newton's method for finding roots of a function
! 
!  f  - function
!  df - derivative of function
!  x0 - initial guess
!  x  - found root
!  


eps = epsilon(1.)

iter = 10**3
x    = x0

do i = 1, iter

! Compute function and derivative at initial guess
fx  = f(x)
dfx = df(x)

! Update guess
h = fx/dfx
x = x - h

if (abs(fx) .lt. eps) exit
end do

if (ibug .eq. 1) print*,'At end, f(x) = ',fx,'iteration',i

end subroutine 


end module
