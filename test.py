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
V	= 45
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
# ========================================
my_prop = Propeller(R, B, V, Clalfa, Cd, rpm, h, xx, x_Point_c, y_Point_c, x_Point_theta, y_Point_theta)
type(x_Point_c)
print(my_prop.theta)
print(my_prop.lamb)
print(my_prop.V_T)
print("Angolo di inflow:\n", np.deg2rad(my_prop.fie))
print("Solidità:\n", my_prop.sigma)
fig1 = plt.figure()
plt.plot(xx,my_prop.chord)
plt.xlabel('$r$')              # x-label to the axes.
plt.ylabel(r'Chord - $c$')     # y-label to the axes.
plt.title(r'Chord diagram')    # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
pp.savefig(fig1)
# ========================================
fig2 = plt.figure()
plt.plot(xx,my_prop.theta)
plt.xlabel('$r$')                     # x-label to the axes.
plt.ylabel(r'Pitch - $\theta$ [deg]') # y-label to the axes.
plt.title(r'Pitch angle diagram')     # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
pp.savefig(fig2)
# ========================================
fig3 = plt.figure()
plt.plot(xx,my_prop.sigma)
plt.xlabel('$r$')                  # x-label to the axes.
plt.ylabel(r'Solidity - $\sigma$') # y-label to the axes.
plt.title(r'Solidity diagram')     # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
pp.savefig(fig3)
# ========================================
fig4 = plt.figure()
plt.plot(xx,my_prop.fie)
plt.xlabel('$r$')                    # x-label to the axes.
plt.ylabel(r'Inflow - $\phi$ [deg]') # y-label to the axes.
plt.title(r'Inflow angle diagram')   # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
pp.savefig(fig4)
# ========================================
m1 = Method_A(my_prop.lamb, my_prop.V_T, my_prop.V_R, my_prop.sigma, my_prop.Clalfa, my_prop.xx, np.deg2rad(my_prop.fie), np.deg2rad(my_prop.theta))
m2 = Method_B(my_prop.lamb, my_prop.B, my_prop.V_T, my_prop.V_R, my_prop.sigma, my_prop.Clalfa, my_prop.xx, np.deg2rad(my_prop.fie), np.deg2rad(my_prop.theta))
# ========================================
fff1 = plt.figure()
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
print("Correzione di Prandtl F:\n", m2.F_correct)
print("alfai:\n", np.rad2deg(m2.alfai))
print("wa:\n", m2.wa)
# ========================================
fig5 = plt.figure()
plt.plot(xx,np.rad2deg(m2.alfai))
plt.xlabel('$r$')                        # x-label to the axes.
plt.ylabel(r'Inflow - $\alpha_i$ [deg]') # y-label to the axes.
plt.title(r'Inflow - Method B')      # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
pp.savefig(fig5)
# ========================================
fig6 = plt.figure()
plt.plot(xx,m2.wa)
plt.xlabel('$r$')                    # x-label to the axes.
plt.ylabel('$w_a$') # y-label to the axes.
plt.title(r'Axial induction - Method B')   # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
pp.savefig(fig6)
pp.close()
# ========================================
m3 = Method_C(my_prop.alpha, my_prop.lamb, my_prop.B, my_prop.V_T, my_prop.V_R, my_prop.sigma, my_prop.Clalfa, my_prop.xx, np.deg2rad(my_prop.fie), np.deg2rad(my_prop.theta), my_prop.Omega, my_prop.R, my_prop.V, Cd, Cl)
print('Prova alfa:\n', my_prop.alpha)
print('Prova induzione a:\n', m3.a)
print('Prova induzione a_prime:\n', m3.a_prime)
print('Prova coefficienti di spinta dCT/dr:\n', m3.dCT)
print('Prova coefficienti di spinta dCP/dr:\n', m3.dCP)
print('Rapporto di funzionamento J:\n', m3.J)
print('Complessivo dati:\n', m3.data.to_string())
# ========================================
fig7 = plt.figure()
plt.plot(xx,my_prop.alpha)
plt.xlabel(r'$r$')                    # x-label to the axes.
plt.ylabel(r'$\alpha$') # y-label to the axes.
plt.title(r'Angle of attack')   # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
#pp.savefig(fig7)
#pp.close()
# ========================================
# ========================================
fig8 = plt.figure()
plt.plot(xx,m3.dCT)
plt.plot(xx,m3.dCP)
plt.xlabel(r'$r$')                    # x-label to the axes.
plt.ylabel(r'$\frac{dC_T}{dr}, \,\, \frac{dC_P}{dr}$') # y-label to the axes.
plt.title(r'Thrust and Power coefficient - Method C')   # Title to the axes.
plt.grid(True,linestyle='-.')
plt.show()
#pp.savefig(fig7)
#pp.close()
# ========================================



