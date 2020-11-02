# -*- coding: utf-8 -*-
#                                                                                                                          
#                                                                                                                          
# PPPPPPPPPPPPPPPPP                                                            PPPPPPPPPPPPPPPPP                            
# P::::::::::::::::P                                                           P::::::::::::::::P                           
# P::::::PPPPPP:::::P                                                          P::::::PPPPPP:::::P                          
# PP:::::P     P:::::P                                                         PP:::::P     P:::::P                         
#   P::::P     P:::::Prrrrr   rrrrrrrrr      ooooooooooo   ppppp   ppppppppp     P::::P     P:::::Pyyyyyyy           yyyyyyy
#   P::::P     P:::::Pr::::rrr:::::::::r   oo:::::::::::oo p::::ppp:::::::::p    P::::P     P:::::P y:::::y         y:::::y 
#   P::::PPPPPP:::::P r:::::::::::::::::r o:::::::::::::::op:::::::::::::::::p   P::::PPPPPP:::::P   y:::::y       y:::::y  
#   P:::::::::::::PP  rr::::::rrrrr::::::ro:::::ooooo:::::opp::::::ppppp::::::p  P:::::::::::::PP     y:::::y     y:::::y   
#   P::::PPPPPPPPP     r:::::r     r:::::ro::::o     o::::o p:::::p     p:::::p  P::::PPPPPPPPP        y:::::y   y:::::y    
#   P::::P             r:::::r     rrrrrrro::::o     o::::o p:::::p     p:::::p  P::::P                 y:::::y y:::::y     
#   P::::P             r:::::r            o::::o     o::::o p:::::p     p:::::p  P::::P                  y:::::y:::::y      
#   P::::P             r:::::r            o::::o     o::::o p:::::p    p::::::p  P::::P                   y:::::::::y       
# PP::::::PP           r:::::r            o:::::ooooo:::::o p:::::ppppp:::::::pPP::::::PP                  y:::::::y        
# P::::::::P           r:::::r            o:::::::::::::::o p::::::::::::::::p P::::::::P                   y:::::y         
# P::::::::P           r:::::r             oo:::::::::::oo  p::::::::::::::pp  P::::::::P                  y:::::y          
# PPPPPPPPPP           rrrrrrr               ooooooooooo    p::::::pppppppp    PPPPPPPPPP                 y:::::y           
#                                                           p:::::p                                      y:::::y            
#                                                           p:::::p                                     y:::::y             
#                                                          p:::::::p                                   y:::::y              
#                                                          p:::::::p                                  y:::::y  #                                                          p:::::::p                                 y:::::y 
#                                                          ppppppppp                                yyyyyyy                       
#                                                                                                                          
#
#
# ====================================================================== 
# An experimental Modules for propeller analysis 
# ====================================================================== 
"""
Created on Fri Oct 23 11:19:41 2020

@author: claum
"""
# =========================================
# Class begins
# ========================================= 
import numpy as np 
import math 
from scipy import integrate
from scipy.optimize import fsolve
import sympy as sp
import pandas as pd 
from sympy.core.symbol import symbols
from sympy.solvers.solveset import nonlinsolve
# =========================================
# Define a Propeller
# ========================================= 
class Propeller(object):
	"""
	Class for propeller input data and analysis
	INPUT 
	R       : Propeller radius (in meters)
	B       : Number of blades
	V       : Free stream velocity (in meters-per-seconds)
	Clalfa  : Blade-element lift curve slope
	Cd      : Blade-element drag coefficient
	rpm     : Rotation-per-minute
	h       : Flight level (in meters)
	D       : Propeller diameter (in meters)
	Omega   : Propeller angular velocity
	xx      : Stations along the propeller radius
	x_c     : Stations along propeller radius for chords (must be an array of six elements)
	c       : Chord evaluated at given stations along propeller radius (must be an array of six elements)
	x_theta : Stations along propeller radius for pitch angle (must be an array of six elements)
	y_theta : Pitch angle evaluate at given stations along propeller radius (must be an array of six elements)
	
	Other values stored inside Propeller:
	lamb    : Propeller advance ratio LAMBDA = V/(Omega*R)
	V_T     : Omega X xx
	V_R	: V_T X sqrt(lamb**2 + xx**2)
	sigma   : Propeller's solidity
	fie     : Propeller's inflow angle
	"""
	def __init__(self, R, B, V, Clalfa, Cd, rpm, h, xx, x_c, c, x_theta, y_theta):
		"""
		Inizialization of the object variable which contains propeller's input data
		Parameters
		----------
		R       --> Propeller radius
		B	--> Number of blades
		V       --> Working velocity
		Clalfa  --> Lift-curve slope of a blade element
		Cd 	--> Drag coefficient of a blade element
		rpm	--> Round per minutes
		h       --> Height from ground
		xx      --> Stations along the propeller radius
		x_c     --> Points along the propeller radius to interpolate c = c(r)
		c       --> Values of c = c(r) sampled along the propeller radius
		x_theta --> Points along the propeller radius to interpolate theta = theta(r)
		y_theta --> Values of theta = theta(r) sampled along the propeller radius
		Output
		----------
		Define a propeller as an object variable
		Example:
		my_prop = PropPy.Propeller(...)
		"""
		self.R, self.B, self.V, self.Clalfa, self.Cd, self.rpm, self.h, self.xx, self.x_c, self.c, self.x_theta, self.y_theta = R, B, V, 			Clalfa, Cd, rpm, h, xx, x_c, c, x_theta, y_theta
		self.n 	            = self.n(self.rpm)
		self.Omega 	    = self.Omega(self.n)
		self.D     	    = self.D(self.R)
		self.chord 	    = self.chord(self.x_c,self.c,self.xx)
		self.theta 	    = self.theta(self.x_theta,self.y_theta,self.xx)
		self.lamb           = self.advance_ratio(self.V, self.Omega, self.R)
		self.V_T	    = self.V_T(self.Omega, self.xx)
		self.V_R	    = self.V_R(self.V_T, self.lamb, self.xx)
		self.sigma          = self.sigma(self.B, self.chord, self.R)
		self.fie            = self.fie(self.lamb, self.xx)
# ========================================
# Define a Propeller
# ========================================
	def n(self, rpm):
		"""
		Function that converts round-per-minutes in round-per-seconds
		Parameters
		----------
		rpm --> Round-per-minutes

		Returns
		-------
		n   --> Round-per-seconds
		"""
		return self.rpm/60
# ========================================
	def Omega(self, n):
		"""
		Function that calculates Omega
		Parameters
		----------
		n     --> Round-per-minutes

		Returns
		-------
		Omega --> Propeller's angular velocity [radians/sec]
		"""
		return self.n*2*math.pi
# ========================================
	def D(self, R):
		"""
		Function that calculates propeller Diameter
		Parameters
		----------
		R --> Propeller's radius [m]

		Returns
		-------
		D --> Propeller's diameter [m]
		"""
		return self.R*2
# ========================================
	def advance_ratio(self, V, Omega, R):
		"""
		Function that calculates the advance ratio lambda
		Parameters
		----------
		V      --> Working velocity
		Omega  --> Propeller's angular velocity [radians/sec]
		R      --> Propeller's radius [m]

		Returns
		-------
		lambda --> Propeller's advance ratio
		"""
		return self.V/(self.Omega*self.R)
# ========================================
# Define chords and pitch angle distribution
# ========================================
	def chord(self, x_c, c, xx):
		"""
		Function that interpolates via Lagrange's polynomials the chord along the propeller's radius
		Parameters
		----------
		x_c --> Points along the propeller's radius in which chord has been sampled
		c   --> c = c(xx) sampled along the propeller's radius
		xx  --> Stations along the propeller's radius

		Returns
		-------
		Chord's values calculated through interpolation
		c_interpolated
		"""
		n = len(x_c)
		sum = 0
		for i in range(n):
			product = c[i]
			for j in range(n):
				if i != j:
					product = product*(xx - x_c[j])/(x_c[i]-x_c[j])
			sum = sum + product
		return sum
# ========================================
	def theta(self, x_theta, y_theta, xx):
		"""
		Function that interpolates via Lagrange's polynomials the pitch angles along the propeller's radius
		Parameters
		----------
		x_theta --> Points along the propeller's radius in which pitch angle's has been sampled
		y_theta --> theta = theta(xx) sampled along the propeller's radius
		xx      --> Stations along the propeller's radius

		Returns
		-------
		Pitch angle's values calculated through interpolation
		theta_interpolated
		"""
		n = len(x_theta)
		sum = 0
		for i in range(n):
			product = y_theta[i]
			for j in range(n):
				if i != j:
					product = product*(xx - x_theta[j])/(x_theta[i]-x_theta[j])
			sum = sum + product
		return sum
# ========================================
	def V_T(self, Omega, xx):
		"""
		Function that calculates V_T
		Parameters
		----------
		Omega --> Propeller's angular velocity
		xx    --> Stations along the propeller radius

		Returns
		-------
		See McCormick: Aerodynamics of VSTOL Flight (1967) - Chapter 4, pag 80
		Omega X xx
		"""
		return self.Omega*self.xx
# ========================================
	def V_R(self, V_T, lamb, xx):
		"""
		Function that calculates V_R
		Parameters
		----------
		V_T  --> Omega X xx
		lamb --> Propeller's advance ratio
		xx   --> Stations along the propeller's radius

		Returns
		-------
		See McCormick: Aerodynamics of VSTOL Flight (1967) - Chapter 4, pag 80
		V_R = V_T X sqrt(lamb**2 + xx**2)
		"""
		return self.V_T*np.sqrt(self.lamb**2 + self.xx**2)
# ========================================
	def sigma(self, B, chord, R):
		"""
		Propeller's solidity sigma
		Parameters
		----------
		B     --> Propeller's blades
		chord --> Propeller's chords along the radius
		R     --> Propeller's radius

		Returns
		-------
		Solidity sigma [Non-dimensional]
		"""
		return (self.B*self.chord)/(2*math.pi*self.R)
# ========================================
	def fie(self, lamb, xx):
		"""
		Propeller's inflow angles
		Parameters
		----------
		lamb --> Propeller's advance ratio
		xx   --> Stations along the propeller's radius

		Returns
		-------
		fie  --> Propeller's inflow angles
		"""
		return np.rad2deg(np.arctan(lamb/xx))
# ========================================
# Rotor class inheritance from Propeller
# ========================================	

class Rotor(Propeller):
	"""
	Class for rotor input data and analysis
	"""
	def __init__(self, R, B, V, Clalfa, Cd, rpm, h, xx, x_c, c, x_theta, y_theta, mu, theta_hub, delta_theta):
		"""
		Inizialization of the object variable which contains rotor's input data
		Parameters
		----------
		R           --> Rotor's radius
		B           --> Number of blades
		V           --> Velocity
		Clalfa      --> Lift slope of the blade element
		Cd          --> Drag coefficient for the blade elements
		rpm         --> Round-per-minutes
		h           --> Height above the ground
		xx          --> Stations along the rotor's radius
		x_c         --> Points along the rotor's radius where chord has been sampled (must be a six-points array)
		c           --> Chord sampled along the rotor's radius (must be a six-points array)
		x_theta	    --> Points along the rotor's radius where pitch angles has been sampled (must be a six-points array)
		y_theta     --> Pitch angles sampled along the rotor's radius (must be a six-points array)
		mu          --> Climb velocity V / (Omega X R)
		theta_hub   --> Pitch angle measured at the hub
		delta_theta --> Pitch angles' lapse-rate
		"""
		Propeller.__init__(self, R, B, V, Clalfa, Cd, rpm, h, xx, x_c, c, x_theta, y_theta)
		self.mu, self.theta_hub, self.delta_theta = mu, theta_hub, delta_theta
		self.chord    = Propeller.chord(self, self.x_c, self.c, self.xx)
		self.sigma    = 2*Propeller.sigma(self, self.B, self.chord, self.R)
		self.theta    = self.Theta(self.theta_hub, self.delta_theta, self.xx)
		self.lambdai  = self.Lambdai(self.mu, self.Clalfa, self.sigma, np.deg2rad(self.theta), self.xx)
		self.phii     = self.Phii(self.mu, self.lambdai, self.xx)
		self.alpha    = self.theta - self.phii
		self.Cl       = self.Clalfa*(np.deg2rad(self.alpha))
		self.dTc      = self.dTc(self.sigma, self.Cl, self.xx)
		self.dQc      = self.dQc(self.sigma, self.Cd, self.Cl, np.deg2rad(self.phii), self.xx)
		self.Tc       = self.TC(self.dTc, self.xx)
		self.Qc       = self.QC(self.dQc, self.xx)
		self.data     = self.my_rotor_table(self.xx, self.theta, self.sigma, self.lambdai, self.phii, self.Cl, self.Cd, self.dTc, self.dQc, 			self.Tc, self.Qc)
# ========================================	
	def Theta(self, theta_hub, delta_theta, xx):
		"""
		Defines pitch's angles variations with a linear law along the propeller's radius
		Parameters
		----------
		theta_hub   --> Pitch angles at the rotor's hub
		delta_theta --> Pitch lapse-rate along the rotor's radius
		xx          --> Stations along the rotor's radius

		Returns
		-------
		Linear variation of pitch angles along the rotor's radius
		theta = theta_hub - delta_theta X xx
		"""
		return self.theta_hub - self.delta_theta*self.xx
# ========================================
	def Lambdai(self, mu, Clalfa, sigma, theta, xx):
		"""
		Axial induction lambda_i
		Parameters
		----------
		mu     --> Climb velocity V / (Omega X R)
		Clalfa --> Lift slope of the blade element
		sigma  --> Rotor's solidity
		theta  --> Pitch angles along the rotor's radius
		xx     --> Stations along the rotor's radius

		Returns
		-------
		Axial induction on the rotor lambda_i, defined as positive root of the following equation

		lambda_i**2 + (mu + Clalfa * sigma/8) * lambda_i - xx * Clalfa * (sigma/8) * (theta - (mu/xx)) = 0

		"""
		bb    = self.mu + (self.Clalfa*self.sigma)/8
		cc    = ((self.xx*self.Clalfa*self.sigma)/8)*(self.theta - (self.mu)/self.xx)
		lambi = -0.5*bb + 0.5*((bb**2 + 4*cc)**0.5)

		return lambi
# ========================================
	def Phii(self, mu, lambdai, xx):
		"""
		Function that calculates inflow angles along the rotor's radius.
		Phi is expressed in [deg].
		Parameters
		----------
		mu      --> Climb velocity V / (Omega X R)
		lambdai --> Rotor's axial induction (see above)
		xx      --> Stations along the rotor's radius

		Returns
		-------
		fi_i [deg] calculated via the equation

	        (mu/xx) + (lambda_i/xx)

		"""
		return (self.mu)/self.xx + (self.lambdai)/self.xx
# ========================================
	def dTc(self, sigma, Cl, xx):
		"""
		Function that calculates dTc/dxx along the rotor's radius, where dTc is the Thrust non-dimensional coefficient
		of the single blade element at the xx station
		Parameters
		----------
		sigma --> Rotor's solidity
		Cl    --> Blade element's lift coefficient
		xx    --> Stations along the rotor's radius

		Returns
		-------
		dTc   --> Blade element's non-dimensional thrust coefficient
		"""
		return 0.5*self.sigma*self.Cl*(self.xx**2)
# ========================================
	def dQc(self, sigma, Cd, Cl, phii, xx):
		"""
		Function that calculates dQc/dxx along the rotor's radius, where dQc is the Torque non-dimensional coefficient
		of the single blade element at the xx station
		Parameters
		----------
		sigma --> Rotor's solidity
		Cd    --> Drag coefficient of the blade element
		Cl    --> Lift coefficient of the blade element
		phii  --> Inflow angle
		xx    --> Stations along the rotor's radius

		Returns
		-------
		dQc   --> Blade element non-dimensional torque coefficient
		"""
		return 0.5*self.sigma*(self.Cd + self.Cl*self.phii)*(self.xx**3)
# ========================================	
	def TC(self, dTc, xx):
		"""
		Rotor's global thrust coefficient evaluated with numerical integration along the rotor's ratio
		Parameters
		----------
		dTc --> Blade element's thrust coefficient
		xx  --> Stations along the rotor's radius

		Returns
		-------
		Tc  --> Rotor's thrust coefficient
		"""
		return integrate.simps(self.dTc, self.xx)
# ========================================	
	def QC(self, dQc, xx):
		"""
		Rotor's global torque coefficient evaluated with numerical integration along the rotor's ratio
		Parameters
		----------
		dQc --> Blade element's torque coefficien
		xx  --> Stations along the rotor's radius

		Returns
		-------
		Qc  --> Rotor's torque coefficient
		"""
		return integrate.simps(self.dQc, self.xx)
# ========================================	
	def my_rotor_table(self, xx, theta, sigma, lambdai, phii, Cl, Cd, dTc, dQc, Tc, Qc):
		"""
		Function that generates a comma separated file with outputs from the calculations above. Also, the function
		print to console all the values in columns format.
		Parameters
		----------
		xx      --> Stations along the rotor's radius
		theta   --> Pitch angles along the rotor's radius
		sigma   --> Rotor's solidity
		lambdai --> Rotor's axial induction
		phii    --> Inflow angles
		Cl      --> Blade element's lift coefficient
		Cd      --> Blade element's drag coefficient
		dTc     --> Blade element's thrust coefficient
		dQc     --> Blade element's torque coefficient
		Tc      --> Rotor's thrust coefficient
		Qc      --> Rotor's torque coefficient

		Returns
		-------
		Comma separated values (.csv) file
		"""
		xx      = self.xx
		theta   = self.theta
		sigma   = self.sigma
		lambdai = self.lambdai
		phii    = self.phii
		Cl      = self.Cl
		Cd      = self.Cd
		dTc     = self.dTc
		dQc     = self.dQc
		Tc      = self.Tc
		Qc      = self.Qc
		my_dict = {'r_bar': xx, 'theta': theta, 'sigma': sigma, 'lambdai': lambdai, 'phi': phii, 'Cl': Cl, 'Cd': Cd, 'dTc/dr_bar': dTc,
		'dQc/dr_bar': dQc}
		df1     = pd.DataFrame(my_dict)
		df2     = pd.DataFrame({'Tc': Tc, 'Qc': Qc}, index=np.arange(1))
		df      = pd.concat([df1, df2])
		df.to_csv("esempio")
		return df
			
# ++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++

# ========================================
# Momentum Theory
# ========================================
class Method_A(object):
	"""
	Class for propeller analysis with momentum theory
	INPUT 
	lamb     : Propeller advance ratio
	V_T      : Omega * xx
	V_R      : V_T * sqrt(lamb**2 + xx**2)
	Clalfa   : Lift slope coefficient 
	fie      : Inflow-angles distribution (must be in radians)
	theta    : Pitch-angles distribution (must be in radians)
	"""
	def __init__(self, lamb, V_T, V_R, sigma, Clalfa, xx, fie, theta):
		"""
		Inizialization of the object variable for Method A
		Parameters
		----------
		lamb     : Propeller advance ratio
		V_T      : Omega * xx
		V_R      : V_T * sqrt(lamb**2 + xx**2)
		sigma    : Propeller's solidity
		Clalfa   : Lift slope coefficient
		fie      : Inflow-angles distribution (must be in radians)
		theta    : Pitch-angles distribution (must be in radians)

		Results
		----------
		alfai    : Induced angle of attack distribution
		wa       : Axial induction
		wt       : Rotational induction
		"""
		self.lamb, self.V_T, self.V_R, self.sigma, self.Clalfa, self.xx, self.fie, self.theta = lamb, V_T, V_R, sigma, Clalfa, xx, fie, theta
		self.alfai = self.alfai(self.lamb, self.V_T, self.V_R, self.sigma, self.Clalfa, self.xx, self.fie, self.theta)
		self.wa    = self.wa(self.V_R, self.alfai, self.fie)
		self.wt    = self.wt(self.V_R, self.alfai, self.fie)
# ========================================
# Induced angle of attack
# ========================================
	def alfai(self, lamb, V_T, V_R, sigma, Clalfa, xx, fie, theta):
		"""
		Function that calculates the induced angle of attack distribution along the propeller's radius with Method A
		Parameters
		----------
		lamb   --> Propeller's advance ratio
		V_T    --> Omega * xx
		V_R    --> V_T * sqrt(lamb**2 + xx**2)
		sigma  --> Propeller's solidity
		Clalfa --> Lift slope coefficient
		xx     --> Stations along the propeller's radius
		fie    --> Inflow-angles distribution (must be in radians)
		theta  --> Pitch angles distribution (must be in radians)

		Returns
		-------
		alfai  --> Induced angle of attack distribution
		"""
		temp1 = (self.lamb)/(self.xx) + (self.sigma*self.Clalfa*self.V_R)/(8*self.V_T*self.xx**2)
		temp2 = temp1**2 + (self.sigma*self.Clalfa*self.V_R*(self.theta-self.fie))/(self.V_T*2*self.xx**2)
		return -0.5*(-temp1 + np.sqrt(temp2))
# ========================================
	def wa(self, V_R, alfai, fie):
		"""
		Function that calculates axial induction with Method A
		Parameters
		----------
		V_R   --> V_T * sqrt(lamb**2 + xx**2)
		alfai --> Induced angles of attack distribution
		fie   --> Inflow-angles distribution

		Returns
		-------
		wa    --> Axial induction
		"""
		return self.V_R*self.alfai*np.cos(self.fie + self.alfai)
# ========================================
	def wt(self, V_R, alfai, fie):
		"""
		Function that calculates rotational induction with Method A
		Parameters
		----------
		V_R   --> V_T * sqrt(lamb**2 + xx**2)
		alfai --> Induced angles of attack distribution
		fie   --> Inflow-angles distribution

		Returns
		-------
		wt    --> Rotational induction
		"""
		return self.V_R*self.alfai*np.sin(self.fie + self.alfai)
# ========================================	

# ++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++

# ========================================
# Small Disturbance Vortex Theory
# ========================================
class Method_B(object):
	"""
	Class for propeller analysis with small-disturbance-vortex-theory
	INPUT 
	lamb     : Propeller advance ratio
	V_T      : Omega * xx
	V_R      : V_T * sqrt(lamb**2 + xx**2)
	Clalfa   : Lift slope coefficient 
	fie      : Inflow-angles distribution (must be in radians)
	theta    : Pitch-angles distribution (must be in radians)
	"""
	def __init__(self, lamb, B, V_T, V_R, sigma, Clalfa, xx, fie, theta):
		"""
		Inizialization of the object variable for Method B
		Parameters
		----------
		lamb   --> Propeller's advance ratio
		B      --> Propeller's number of blades
		V_T    --> Omega * xx
		V_R    --> V_T * sqrt(lamb**2 + xx**2)
		sigma  --> Propeller's solidity
		Clalfa --> Lift slope coefficient
		xx     --> Stations along the propeller's radius
		fie    --> Inflow-angles distribution (must be in radians)
		theta  --> Pitch-angles distribution (must be in radians)
		Omega  --> Propeller's angular velocity
		R      --> Propeller's radius
		V      --> Working velocity
		"""
		self.lamb, self.B, self.V_T, self.V_R, self.sigma, self.Clalfa, self.xx, self.fie, self.theta = lamb, B, V_T, V_R, sigma, Clalfa, xx, 			fie, theta
		self.F_correct     = self.F_correct(self.B, self.xx, self.fie)
		self.alfai 	   = self.alfai(self.lamb, self.V_T, self.V_R, self.sigma, self.Clalfa, self.xx, self.fie, self.theta, self.F_correct)
		self.wa    	   = self.wa(self.V_R, self.alfai, self.fie)
		self.wt    	   = self.wt(self.V_R, self.alfai, self.fie)
# ========================================
# Induced angle of attack
# ========================================
	def F_correct(self, B, xx, fie):
		"""
		Prandtl's tip-loss factor (see for instance Aerodynamics of VSTO Flight, Chapter 4, pag 83, formula 4-29)
		Parameters
		----------
		B   --> Propeller's number of blades
		xx  --> Stations along the propeller's radius
		fie --> Inflow-angles distribution

		Returns
		-------
		F   --> Tip-loss factor
		"""
		return (2/math.pi)*np.arccos(np.exp(-(self.B*(1-self.xx))/(2*np.sin(self.fie))))
# ========================================
	def alfai(self, lamb, V_T, V_R, sigma, Clalfa, xx, fie, theta, F_correct):
		"""
		Induced angles of attack distribution along the propeller's radius with Method B
		Parameters
		----------
		lamb      --> Propeller's advance ratio
		V_T       --> Omega * xx
		V_R       --> V_T * sqrt(lamb**2 + xx**2)
		sigma     --> Propeller's solidity
		Clalfa    --> Lift slope coefficient
		xx        --> Stations along the propeller's radius
		fie       --> Inflow-angles distribution
		theta     --> Pitch-angles distribution
		F_correct --> Prandtl's tip-loss factor

		Returns
		-------
		alfai     --> Induced angles of attack distribution along the propeller's radius
		"""
		temp1 = (self.lamb)/(self.xx) + (self.sigma*self.Clalfa)/(8*(self.F_correct)*(self.xx)*np.cos(self.fie))
		temp2 = temp1**2 + (self.sigma*self.Clalfa*(self.theta - self.fie))/(2*(self.F_correct)*(self.xx)*np.cos(self.fie)) 
		return -0.5*(-temp1 + np.sqrt(temp2))
# ========================================
	def wa(self, V_R, alfai, fie):
		"""
		Function that calculates axial induction with Method B
		Parameters
		----------
		V_R   --> V_T * sqrt(lamb**2 + xx**2)
		alfai --> Induced angle of attack distribution
		fie   --> Inflow-angles distribution

		Returns
		-------
		wa    --> Axial induction
		"""
		return self.V_R*self.alfai*np.cos(self.fie + self.alfai)
# ========================================
	def wt(self, V_R, alfai, fie):
		"""
		Function that calculates rotational induction with Method B
		Parameters
		----------
		V_R   --> V_T * sqrt(lamb**2 + xx**2)
		alfai --> Induced angle of attack distribution
		fie   --> Inflow-angles distribution

		Returns
		-------
		wt    --> Rotational induction
		"""
		return self.V_R*self.alfai*np.sin(self.fie + self.alfai)
# ========================================	

# ========================================
# Full Vortex Theory
# ========================================
class Method_C(object):
	"""
	Class for propeller analysis with small-disturbance-vortex-theory
	INPUT 
	lamb     : Propeller advance ratio
	V_T      : Omega * xx
	V_R      : V_T * sqrt(lamb**2 + xx**2)
	Clalfa   : Lift slope coefficient 
	fie      : Inflow-angles distribution (must be in radians)
	theta    : Pitch-angles distribution (must be in radians)
	"""
	def __init__(self, lamb, B, V_T, V_R, sigma, Clalfa, xx, fie, theta, Omega, R, V):
		self.lamb, self.B, self.V_T, self.V_R, self.sigma, self.Clalfa, self.xx, self.fie, self.theta, self.Omega, self.R, self.V = lamb, B, 			V_T, V_R, sigma, Clalfa, xx, fie, theta, Omega, R, V
		self.F_correct     = self.F_correct(self.B, self.xx, self.fie)
		self.induction     = self.induction(self.lamb, self.V_T, self.V_R, self.sigma, self.Clalfa, self.xx, self.fie, self.theta, 			self.F_correct, self.Omega, self.R, self.V)
# ========================================
# Induced angle of attack
# ========================================
	def F_correct(self, B, xx, fie):
		return (2/math.pi)*np.arccos(np.exp(-(self.B*(1-self.xx))/(2*np.sin(self.fie))))
# ========================================
	def induction(self, lamb, V_T, V_R, sigma, Clalfa, xx, fie, theta, F_correct, Omega, R, V):
		#x = np.array(3)
		# ========================================
		# First approximation
		# ========================================
		temp1 = (self.lamb)/(self.xx) + (self.sigma*self.Clalfa)/(8*(self.F_correct)*(self.xx)*np.cos(self.fie))
		temp2 = temp1**2 + (self.sigma*self.Clalfa*(self.theta - self.fie))/(2*(self.F_correct)*(self.xx)*np.cos(self.fie)) 
		# ========================================
		alfai1 = -0.5*(-temp1 + np.sqrt(temp2))
		wa1    = self.V_R*alfai1*np.cos(self.fie + alfai1)
		wt1    = self.V_R*alfai1*np.sin(self.fie + alfai1)
		# ========================================
		x      = self.xx 
		lamb   = self.lamb
		Omega  = self.Omega
		vt     = self.V_T
		vr     = self.V_R
		r      = self.R
		v      = self.V	
		theta  = self.theta	
		f      = self.F_correct
		sigma  = self.sigma
		Clalfa = self.Clalfa
		y      = np.zeros(len(alfai1))
		dy     = np.zeros(len(alfai1))
		term1  = np.zeros(len(alfai1))
		term2  = np.zeros(len(alfai1)) 
		aa     = np.zeros(len(alfai1))
		bb     = np.zeros(len(alfai1))
		cc     = np.zeros(len(alfai1))
		da     = np.zeros(len(alfai1))
		db     = np.zeros(len(alfai1))
		dc     = np.zeros(len(alfai1))
		err    = 1E-9
		wt     = np.zeros(len(alfai1))
		wa     = np.zeros(len(alfai1))
		alfai  = np.zeros(len(alfai1))
		for i in range(len(alfai1)):
			wt1[i] = wt1[i] - 0.1*wt1[i]
			term1[i] = np.sqrt(v**2 + 4*wt1[i]*(Omega*x[i]*R - wt1[i])) - v
			term2[i] = lamb + np.sqrt(lamb**2 + 4*(wt1[i]/vt[i])*(x[i] - (wt1[i]/vt[i])))
			aa[i]    = theta[i] - np.arctan(2*wt1[i]/term1[i])
			bb[i]    = np.sqrt(0.25*(term2[i])**2 + (x[i] - (wt1[i]/vt[i]))**2)
			cc[i]    = 8*x[i]*f[i]*(wt1[i]/vt[i])
			da[i]    = ((-2)/((term1[i]**2)+(4*wt1[i]*wt1[i])))*(term1[i]-(wt1[i]*(((2*Omega*r*x[i])-(4*wt1[i]))/(term1[i]+v))))
			db[i]    = (((term2[i]*((x[i]/vt[i])-(2*wt1[i]/(vt[i]**2))))/(2*(term2[i]-lamb)))+(wt1[i]/(vt[i]**2))-(x[i]/vt[i]))/bb[i]
			dc[i]    = 8*x[i]*f[i]/vt[i]
			y[i]     = sigma[i]*Clalfa*aa[i]*bb[i] - cc[i]
			dy[i]    = sigma[i]*Clalfa*(aa[i]*db[i] + bb[i]*da[i]) - dc[i]
			for j in range(200):
				if abs(y[j]/dy[i]) < err:
					wt1[i] = wt1[i] - (y[i]/dy[i])
					wt[i]  = wt1[i]
			else:
				continue
				 
		return wt1, wt
			
		
			

			
	

	
 
		
