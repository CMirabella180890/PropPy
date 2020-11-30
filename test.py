from PropPy import Propeller, Method_A, Method_B, Method_C
import numpy as np
import math
import matplotlib.pyplot as plt 
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# ========================================
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('multipage.pdf')
# ========================================
R 	= 1.587 
B 	= 3
V	= 35
Clalfa	= 2*math.pi
Cd	= 0.02
rpm 	= 1800
h	= 3500
xx      = np.linspace(0,1,100)
index   = [0,1,2,3,4,5,6,7,8,9]
xx      = np.delete(xx, index)
x_Point_c = [0.150,0.385,0.485,0.660,0.846,0.950]
y_Point_c = [0.178,0.218,0.217,0.191,0.135,0.088]
Cd        = 0.01
Cl        = 0.6

x_Point_theta = [0.150, 0.289, 0.489, 0.753, 0.859, 0.981] 
y_Point_theta = [17.2 , 10.35, 4.00 ,-1.08 ,-1.820,-2.950]

#x_Point_theta = [0.150, 0.385, 0.485,0.660,0.846,0.950] 
#y_Point_theta = [60   , 45   , 38   ,30   ,15   ,0.5]
# ========================================
# DEFINISCO UN ELICA
# ========================================
my_prop = Propeller(R, B, V, Clalfa, Cd, rpm, h, xx, x_Point_c, y_Point_c, x_Point_theta, y_Point_theta)
type(x_Point_c)
print("THETA:\n", my_prop.theta)
print("INFLOW:\n", my_prop.fie)
print("CHORD:\n", my_prop.chord)
print("ALFA:\n", my_prop.alpha)
print("STAZIONI LUNGO IL RAGGIO:\n", my_prop.xx)
print("Rapporto di avanzamento:\n", my_prop.lamb)
print("Velocità al tip:\n", my_prop.V_T)
print("Solidità:\n", my_prop.sigma)
print('Complessivo dati:\n', my_prop.data.to_string())
# ========================================
fig1 = plt.figure()
plt.plot(xx,my_prop.chord)
plt.xlabel('$r$')              # x-label to the axes.
plt.ylabel(r'Chord - $c$')     # y-label to the axes.
plt.title(r'Chord diagram')    # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
#pp.savefig(fig1)
# ========================================
fig2 = plt.figure()
plt.plot(xx,my_prop.theta)
plt.xlabel('$r$')                     # x-label to the axes.
plt.ylabel(r'Pitch - $\theta$ [deg]') # y-label to the axes.
plt.title(r'Pitch angle diagram')     # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
#pp.savefig(fig2)
# ========================================
fig3 = plt.figure()
plt.plot(xx,my_prop.sigma)
plt.xlabel('$r$')                  # x-label to the axes.
plt.ylabel(r'Solidity - $\sigma$') # y-label to the axes.
plt.title(r'Solidity diagram')     # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
#pp.savefig(fig3)
# ========================================
fig4 = plt.figure()
plt.plot(xx,my_prop.fie)
plt.xlabel('$r$')                    # x-label to the axes.
plt.ylabel(r'Inflow - $\phi$ [deg]') # y-label to the axes.
plt.title(r'Inflow angle diagram')   # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
#pp.savefig(fig4)
# ========================================
# METODO A
# ========================================
m1 = Method_A(my_prop.lamb, my_prop.V_T, my_prop.V_R, my_prop.sigma, my_prop.Clalfa, my_prop.xx, np.deg2rad(my_prop.fie), np.deg2rad(my_prop.theta), my_prop.B, my_prop.R)
# ========================================
fig5 = plt.figure()
plt.plot(xx,m1.wa)
plt.xlabel(r'$r$')                    # x-label to the axes.
plt.ylabel(r'$w_a$') # y-label to the axes.
plt.title(r'Axial induction - Method A')   # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
#pp.savefig(fig6)
#pp.close()
# ========================================
# ========================================
fig6 = plt.figure()
plt.plot(xx,m1.wt)
plt.xlabel(r'$r$')                    # x-label to the axes.
plt.ylabel(r'$w_t$') # y-label to the axes.
plt.title(r'Rotational induction - Method A')   # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
#pp.savefig(fig6)
#pp.close()
# ========================================
# ========================================
fig7 = plt.figure()
plt.plot(xx, np.rad2deg(m1.alfai))
plt.xlabel(r'$r$')                    # x-label to the axes.
plt.ylabel(r'$\alpha_i$ [deg]') # y-label to the axes.
plt.title(r'Induced angle of attack - Method A')   # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
#pp.savefig(fig6)
#pp.close()
# ========================================
# ========================================
fig8 = plt.figure()
plt.plot(xx,m1.dCT)
plt.xlabel(r'$r$')                    # x-label to the axes.
plt.ylabel(r'$\frac{dC_T}{dr}$') # y-label to the axes.
plt.title(r'Thrust coefficient - Method A')   # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
#pp.savefig(fig6)
#pp.close()
# ========================================
# ========================================
print('Complessivo dati:\n', m1.data.to_string())
print("alfai:\n", np.rad2deg(m1.alfai))
print("wa:\n", m1.wa)
print("dCT/dr:\n", m1.dCT)
print("CT:\n", m1.CT)
# ========================================
# METODO B
# ========================================
#m2 = Method_B(my_prop.lamb, my_prop.B, my_prop.V_T, my_prop.V_R, my_prop.sigma, my_prop.Clalfa, my_prop.xx, np.deg2rad(my_prop.fie), np.deg2rad(my_prop.theta))
# ========================================
#fig5 = plt.figure()
#plt.plot(xx,np.rad2deg(m2.alfai))
#plt.xlabel('$r$')                        # x-label to the axes.
#plt.ylabel(r'Inflow - $\alpha_i$ [deg]') # y-label to the axes.
#plt.title(r'Inflow - Method B')      # Title to the axes.
#plt.grid(True,linestyle='-.')
#plt.show()
#pp.savefig(fig5)
# ========================================
#fig6 = plt.figure()
#plt.plot(xx,m2.wa)
#plt.xlabel('$r$')                    # x-label to the axes.
#plt.ylabel('$w_a$') # y-label to the axes.
#plt.title(r'Axial induction - Method B')   # Title to the axes.
#plt.grid(True,linestyle='-.')
#plt.show()
#pp.savefig(fig6)
#pp.close()
# ========================================



