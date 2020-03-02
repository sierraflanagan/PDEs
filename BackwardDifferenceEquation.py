# -*- coding: utf-8 -*-
"""
Code for the Backward difference equation 
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
from numpy.linalg import inv

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


def mastery(x):
    "This is the function given for mastery check"
    y = x*(x-1)*(x+1)
    return y

u_1 = np.zeros(n+1)
        
for i in range(0, n):
    u_1[i] = mastery(x[i])
    

#Vector w alpha and beta values
b = np.zeros(n-1)
b[0] = u_1[0]
b[n-2] = u_1[n]
bt = np.matrix.transpose(b)
#Vector with the initial u conditions
u = np.zeros(n-1)
for i in range (0, (n-1)):
    u[i] = u_1[i+1]
    
ut = np.matrix.transpose(u)
#print(ut)

size = (n-1, n-1)
A_hat = np.zeros(size)
#print(A)    

for i in range(n - 1):
    A_hat[i, i] = 1 + 2*m
    
for i in range(n - 2):
    A_hat[i, i + 1] = -m
    
for i in range(n-2):
    A_hat[i+1, i] = -m
    

A_hat_inv = inv(A_hat)


for i in range(iterations):
    u_b = ut - bt
    u_j = np.matmul(A_hat_inv, u_b)
    u_t = u_j
        
Solution = np.zeros(n+1)
Solution[0] = b[0]
Solution[n] = b[n-2]
for i in range (0, n-1):
    Solution[i+1] = u_j[i]

print(Solution)
#plt.plot(Solution)
plt.plot(x, Solution)
plt.ylabel('Temperature')
plt.xlabel('Position In X')
plt.title('Backward Difference Method')
plt.show()
