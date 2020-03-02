# -*- coding: utf-8 -*-
"""
Code for the Forward difference equation 
Partial Differential Equations
@author: flanagans

TRY with :
    Ut = 5*Uxx
    deltax = .01
    U(0,x) = x*(x-1)*(x+1)
    -1 <= x <=1
    t >=0
"""

import numpy as np
import matplotlib.pyplot as plt
xleft = -1
xright = 1
delta_t = 0.001
tfinal = 0.01
iterations = int(tfinal / delta_t)
delta_x = .1


gamma = 5
length = 1
n = int((xright - xleft)/delta_x) #this value comes from dividing l by delta_t
delta_x_squared = (delta_x)**2 #delta x = .1

m = (gamma*delta_t)/delta_x_squared #formula for calculating mu
#print(m)
#acceptable range of x values which will be evaluated by the piecewise function stored as array
x=np.linspace(-1, 1, n+1)
u_mastery = np.zeros(21)
#x2=np.zeros(11)

#print(x)
#print(x <= .2)
#print((x < 0.7) & (x >= 0.2))

#u_0 = np.piecewise(x, [x < 0.2, (x < 0.7) & (x >= 0.2), x >= 0.7], 
#                       [lambda x: -x, lambda x: x - 0.4, lambda x: 1-x ])
#print(u_0)


def mastery(x):
    "This is the function given for mastery check"
    y = x*(x-1)*(x+1)
    return y

u_0 = np.zeros(n+1)
        
for i in range(0, n):
    u_0[i] = mastery(x[i])
    

#Vector w alpha and beta values
b = np.zeros(n-1)
b[0] = u_0[0]
b[n-2] = u_0[n]
bt = np.matrix.transpose(b)
#Vector with the initial u conditions
u = np.zeros(n-1)
for i in range (0, (n-1)):
    u[i] = u_0[i+1]
    
ut = np.matrix.transpose(u)
#print(ut)

size = (n-1, n-1)
A = np.zeros(size)
#print(A)    

for i in range(n - 1):
    A[i, i] = 1 - 2*m
    
for i in range(n - 2):
    A[i, i + 1] = m
    
for i in range(n-2):
    A[i+1, i] = m
    
#print(A)

for i in range(iterations):
    ICs = np.matmul(A, ut)
    U_j = ICs + bt
    ut = U_j
    

Solution = np.zeros(n+1)
Solution[0] = b[0]
Solution[n] = b[n-2]
for i in range (0, n-1):
    Solution[i+1] = U_j[i]

#print(Solution)
#plt.plot(Solution)
plt.plot(x, Solution)
plt.ylabel('Temperature')
plt.xlabel('Position in X')
plt.title('Forward Difference Method')
plt.show()
