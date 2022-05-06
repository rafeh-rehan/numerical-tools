program test_integrate

use tools
use integrate

implicit double precision (a-h,o-z)

! Number of intervals (multiple of 3)
n = 450 


do i = 1, 10**4
call random_number(dd)
a = dd
call random_number(dd)
b = 10.0d0**5*dd

call simp(f,a,b,n,y)
z = f_int(a,b)


error = abs(y - z)/abs(z)
e_max = max(error,e_max)

end do
print*,'Max relative error simp:', e_max



contains


function f(x)
f = log(x)
end function


function f_int(a,b)
pi = acos(-1.0d0)
f_int = b*log(b)-b -(a*log(a)-a)
end function


end program
