program test_roots


use tools
use roots

implicit real (a-h,o-z)
 
! Initialize error. Test algorithm 10**3 times 
error = 0.0
ibug = 0
do k = 1, 10**3

call random_number(dd)
dd = -25.0 + 50.0*dd ! Random number in [-25,25] as initial root guess

call newt(f,df,dd,root,ibug)
error = max(abs(f(root)), error)

end do 

print*, 'Maximum value of f(root)', error


contains 


! Polynomial test functions 
function f(x)

f = (x-1d0)*(x-3d0)*(x-10d0)

end function

function df(x)

df = (x-3d0)*(x-10d0)+(x-1d0)*(x-10d0)+(x-1d0)*(x-3d0)

end function



end program
