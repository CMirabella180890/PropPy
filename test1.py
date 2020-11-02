from PropPy import Rotor 
import math 
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
# ========================================
R 	    = 7.60
B 	    = 3
V	    = 45
Clalfa	    = 2*math.pi
Cd	    = 0.01
rpm 	    = 1800
h	    = 3500
xx          = np.linspace(0,1,50)
index       = [0,1,2,3,4,5,6,7,8,9]
xx          = np.delete(xx, index)
mu          = 0
theta_hub   = 15 # [deg]
delta_theta = 4.953

x_Point_c = [0.150,0.385,0.485,0.660,0.846,0.950]
y_Point_c = [0.178,0.218,0.217,0.191,0.135,0.088]

x_Point_theta = [0.150,0.289,0.489,0.753,0.859,0.981] 
y_Point_theta = [15.7 ,10.5  ,7.1 ,5.3  ,3.020,1.611]
# ========================================
my_rotor = Rotor(R, B, V, Clalfa, Cd, rpm, h, xx, x_Point_c, y_Point_c, x_Point_theta, y_Point_theta, mu, theta_hub, delta_theta)
# ========================================
fig1 = plt.figure()
plt.plot(xx, my_rotor.theta)
plt.xlabel('$r$')              # x-label to the axes.
plt.ylabel(r'Theta - $\theta$')     # y-label to the axes.
plt.title(r'Theta diagram')    # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
# ========================================

print(my_rotor.mu)
print(my_rotor.xx)
print("Corde:\n", my_rotor.chord)
print("Solidità:\n", my_rotor.sigma)
print("Theta:\n", my_rotor.theta)
print("Lambda-i:\n", my_rotor.lambdai)
print("Phi-i:\n", my_rotor.phii)
print("Alpha-i:\n", my_rotor.alpha)
print("Cl:\n", my_rotor.Cl)
print("dTc:\n", my_rotor.dTc)
print("dQc:\n", my_rotor.dQc)
print("Tc:\n", my_rotor.Tc)
print("Qc:\n", my_rotor.Qc)
print("DATA:\n", my_rotor.data.to_string())
