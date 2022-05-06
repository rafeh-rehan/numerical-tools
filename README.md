# numerical-tools
Some numerical analysis algorithms for integration and interpolation. 

integrate.f90
- Uses Simpson's 3/8 rule to compute definite integrals numerically
- Compiled using:
  gfortran tools.f90 integrate.f90 _module_.f90
- Using _module_ = test_integrate tests the accuracy of the algorithm

interpolate.f90
- Uses 1D barycentric interpolation and 2D bilinear interpolation 
  to numerically interpolate data
- Compiled using:
  gfortran interpolate.f90 _module_.f90
- Using _module_ = test_interpolate tests the accuracy of the algorithm
