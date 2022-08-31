import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as int

"""
This program numerically solves the simple harmonic oscillator's differential equation
in two ways and compares the results to the analytic solutions. The first numerical
method is to use scipy.integrate.odeint to numerically solve the ode. Second, using a
Forward Euler time stepping algorithm which is symplectic (energy conserving).

https://en.wikipedia.org/wiki/Harmonic_oscillator

"""


#Initial Conditions
x0 = 1
v0 = 0

#Time step for numerically solving ode
dt=0.01
t= np.arange(0.0,50.0,dt)
npts= len(t)

#Analytic solutions

x_a = np.cos(t)
v_a = -np.sin(t)

#Numerical solutions

#define a function for our system of first order odes and initial conditions
ics = [x0,v0]
def pend(y,t):
    x,v = y
    dydt = [v , -x]
    return dydt

# ODE solution using scipy's odeint
sol = int.odeint(pend, ics, t)
x_ode = sol[:,0]
v_ode = sol[:,1]

#Total energy of simple harmonic oscillator with k = m = 1
E = 0.5*(x_a**2 + v_a**2)

v_n = np.zeros(npts)
v_n[0] = v0 # set initial condition

x_n = np.zeros(npts)
x_n[0] = x0 # set initial condition

# Forward Euler time stepping algorithm
for i in range(1,npts):
    a_n= -x_n #using k = m = 1 and a Forward Euler time step

    v_n[i] = v_n[i-1] + a_n[i-1]*dt

    x_n[i] = x_n[i-1] + v_n[i]*dt

#Hard to see analytic soln b/c it overlaps almost entireley with numerical soln
plt.figure(1, figsize=(10,10))

#Velocity plot
plt.subplot(3,1,1)
plt.plot(t,v_a, t,v_n, t,v_ode)
plt.legend((r'Analytic',r'Numerical',r'ODE'), loc='lower right')
plt.ylabel('v(t)')
plt.grid('on')

#Position plot
plt.subplot(3,1,2)
plt.plot(t,x_a, t,x_n, t,x_ode)
plt.legend((r'Analytic',r'Numerical',r'ODE'), loc='lower right')
plt.ylabel('x(t)')
plt.xlabel('t')
plt.grid('on')

#Energy plot
plt.subplot(3,1,3)
plt.plot(t,E, label='Energy')
plt.legend(loc='lower right')
plt.ylabel('E(t)')
plt.xlabel('t')
plt.grid('on')

plt.suptitle('Numerical, Analytic, and ODE Solutions to SHO dt=' + ('%.3f' % dt))
