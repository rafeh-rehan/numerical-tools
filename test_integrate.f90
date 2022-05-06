program test_integrate

use tools
use integrate

implicit double precision (a-h,o-z)

n = 450 ! take Big_O( (b-a)/2 * 10 )


do i = 1, 10**4
call random_number(dd)
a = dd
call random_number(dd)
b = 10.0d0**5*dd

call simp(f,a,b,n,y)
call gauss(f,a,b,n,f_integral)
z = f_int(a,b)


error = abs(y - z)/abs(z)
e_max = max(error,e_max)

gerror = abs(f_integral - z)/abs(z)
ge_max = max(gerror,ge_max)

end do
print*,'Max relative error simp:', e_max
print*,'Max relative error gauss:', ge_max



contains


function f(x)
f = log(x)
end function


function f_int(a,b)
pi = acos(-1.0d0)
f_int = b*log(b)-b -(a*log(a)-a)
end function


end program
