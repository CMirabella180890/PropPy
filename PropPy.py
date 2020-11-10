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
	V_T     : Omega X R
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
		self.V_T	    = self.V_T(self.Omega, self.R)
		self.V_R	    = self.V_R(self.V_T, self.lamb, self.xx)
		self.sigma          = self.sigma(self.B, self.chord, self.R)
		self.fie            = self.fie(self.lamb, self.xx)
		self.alpha          = self.theta - self.fie 
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
	def V_T(self, Omega, R):
		"""
		Function that calculates V_T
		Parameters
		----------
		Omega --> Propeller's angular velocity
		R     --> Propeller's radius

		Returns
		-------
		See McCormick: Aerodynamics of VSTOL Flight (1967) - Chapter 4, pag 80
		V_tip = Omega X R
		"""
		return self.Omega*self.R
# ========================================
	def V_R(self, V_T, lamb, xx):
		"""
		Function that calculates V_R
		Parameters
		----------
		V_T  --> Omega X R
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
	V_T      : Omega * R
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
		V_T      : Omega * R
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
		V_T    --> Omega * R
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
		temp2 = temp1**2 + (self.sigma*self.Clalfa*self.V_R*(self.theta - self.fie))/(self.V_T*2*self.xx**2)
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
	V_T      : Omega * R
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
		V_T    --> Omega * R
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
		V_T       --> Omega * R
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
	alfa     : Angle of attack
	lamb     : Propeller advance ratio
	V_T      : Omega * R
	V_R      : V_T * sqrt(lamb**2 + xx**2)
	Clalfa   : Lift slope coefficient 
	fie      : Inflow-angles distribution (must be in radians)
	theta    : Pitch-angles distribution (must be in radians)
	"""
	def __init__(self, alpha, lamb, B, V_T, V_R, sigma, Clalfa, xx, fie, theta, Omega, R, V, Cd, Cl):
		self.alpha, self.lamb, self.B, self.V_T, self.V_R, self.sigma, self.Clalfa, self.xx, self.fie, self.theta, self.Omega, self.R, self.V, 			self.Cd, self.Cl   = alpha, lamb, B, V_T, V_R, sigma, Clalfa, xx, fie, theta, Omega, R, V, Cd, Cl
		self.F_correct     = self.F_correct(self.B, self.xx, self.fie)
		self.lamb1         = self.lamb1(self.Cl, self.Cd, self.fie)
		self.lamb2         = self.lamb2(self.Cl, self.Cd, self.fie)
		self.a             = self.a(self.sigma, self.lamb1, self.fie)
		self.a_prime       = self.a_prime(self.sigma, self.lamb2, self.fie)
		self.dCT           = self.F_correct*self.dCTdr(self.sigma, self.lamb1, self.xx, self.a_prime, self.fie)
		self.dCQ           = self.dCQdr(self.sigma, self.lamb2, self.xx, self.a_prime, self.fie)
		self.dCP           = self.dCQ*2*math.pi
		self.J             = self.advance_ratio(self.xx, self.a_prime, self.a, self.fie)
		self.data          = self.my_prop_table(self.xx, self.alpha, self.Cl, self.Cd, self.fie, self.lamb1, self.lamb2, self.a, self.a_prime, 			self.dCT, self.dCQ, self.J)
# ========================================
# Induced angle of attack
# ========================================
	def F_correct(self, B, xx, fie):
		return (2/math.pi)*np.arccos(np.exp(-(self.B*(1-self.xx))/(2*np.sin(self.fie))))
# ========================================
	def lamb1(self, Cl, Cd, fie):
		return self.Cl*np.cos(self.fie) - self.Cd*np.sin(self.fie)
	def lamb2(self, Cl, Cd, fie):
		return self.Cl*np.sin(self.fie) - self.Cd*np.cos(self.fie)
# ========================================
	def a(self, sigma, lamb1, fie):
		return  ((self.sigma*self.lamb1)/(2*(1-np.cos(2*self.fie))))/(1 - 
			(self.sigma*self.lamb1)/(2*(1-np.cos(2*self.fie))))
	def a_prime(self, sigma, lamb2, fie):
		return ((self.sigma*self.lamb2)/(2*np.sin(2*self.fie)))/(1 + 
			(self.sigma*self.lamb2)/(2*np.sin(2*self.fie)))
	def V_e(self, V, a, fie):
		"""
		Returns V equivalent squared!
		"""
		return ((self.V**2)*(1+self.a)**2)/((np.sin(self.fie))**2)
# ========================================
	def dCTdr(self, sigma, lamb1, xx, a_prime, fie):
		return ((math.pi**3)/4)*self.sigma*self.lamb1*(self.xx**3)*((1 - self.a_prime)**2/((np.cos(self.fie))**2))
	def dCQdr(self, sigma, lamb2, xx, a_prime, fie): 
		return ((1-self.a_prime)**2/(np.cos(self.fie))**2)*((math.pi**3)/8)*self.sigma*self.lamb2*(self.xx**4)
# ========================================
	def advance_ratio(self, xx, a_prime, a, fie):
		return ((1-self.a_prime)/(1+self.a))*math.pi*self.xx*np.tan(self.fie)
# ========================================	
	def my_prop_table(self, xx, alpha, Cl, Cd, fie, lamb1, lamb2, a, a_prime, dCT, dCQ, J):
		"""
		Function that generates a comma separated file with outputs from the calculations above. Also, the function
		print to console all the values in columns format.
		Parameters
		----------

		Returns
		-------
		Comma separated values (.csv) file
		"""
		xx      = self.xx
		alpha   = self.alpha
		Cl      = self.Cl
		Cd      = self.Cd
		fie     = self.fie
		lamb1   = self.lamb1
		lamb2   = self.lamb2
		a       = self.a
		a_prime = self.a_prime
		dCT     = self.dCT
		dCQ     = self.dCQ
		J       = self.J
		my_dict = {'xx': xx, 'alpha': alpha, 'Cl': Cl, 'Cd': Cd, 'fie': fie, 'lamb1': lamb1, 'lamb2': lamb2, 'a': a, 'a_prime': a_prime,
		'dCT/dr_bar': dCT, 'dCQ/dr_bar': dCQ, 'J': J}
		df1     = pd.DataFrame(my_dict)
		#df2     = pd.DataFrame({'Tc': Tc, 'Qc': Qc}, index=np.arange(1))
		#df      = pd.concat([df1, df2])
		#df1.to_csv("esempio")
		return df1
			
# ++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++

	
		
			
		
			

			
	

	
 
		
