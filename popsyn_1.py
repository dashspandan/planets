#simulate planet distribution given stellar and disk parameters iteration 2 with realistic disk parameters from Lupus survey

#Star mass range is in 0.08 to 0.2 solar masses only

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import random as rn
from coeffsstandard import getcoeffdisk
from buildplanet import buildplanet
from hillradius import findorbits
import massselection as ms

a = 1000 # Number of simulations required

#s_axis = [] # List of semi major axis of planets
#mass_planet = [] #List of simulated planet masses

star_mass_list, disk_mass_list = ms.mass_selection_pop_syn_lupus(a)  # For disks from Lupus survel relation
#star_mass_list, disk_mass_list = ms.mass_selection_pop_syn_ormel(a) # For disks with Ormel's criterion

age_list = np.linspace(1e6, 1e7, 10000).tolist()

sim = 0

while sim < len(star_mass_list):

	p = 1. #surface density profile factor, must be less than 2
	z_o = 0.02 # metallicity of disk
	z = 0.5   # silicate to pebble ratio (From Ormel and Yamila, personal communication)
	m_star = star_mass_list[sim] # star mass in solar masses 	
	m_disk = disk_mass_list[sim] #disk mass in stellar masses
	M_dot_accre = 1e-10 # in solar masses per year
	r_out = 200. #outer radius of disk (To make sufficient pebbles)
	h_a = 0.03 #disk aspect ratio, assumed constant for the inner disk (Ormel 2017)
	disk_edge = 0.0102*np.power(0.08/m_star, 1./7.) # in au (from Ormel 2017)
	m_planet_isolation = m_star*np.power(h_a, 3.)*333060.402 #convert from solar to earth masses
	iceline = 0.125*(m_star/0.1) #From Ormel 2017

	orbits = findorbits(m_star, m_planet_isolation, disk_edge, iceline)
	m_dot_2D_coeff, m_dot_3D_coeff, y, a_dot_coeff = getcoeffdisk(p, z_o, z, m_star, m_disk, M_dot_accre, r_out, h_a)

	count = 1 #initiate counter
	o = rn.choice(age_list) #lifetime of disk
	buildplanet(0., count, o, m_dot_2D_coeff, m_dot_3D_coeff, y, a_dot_coeff, m_planet_isolation, orbits, iceline)

	print 'sim %d completed' %sim
	sim = sim+1


#TRAPPIST system parameters

trm = [1.02, 1.16, 0.3, 0.77, 0.93, 1.15, 0.33]
trsa = [0.012, 0.016, 0.022, 0.029, 0.038, 0.047, 0.062]
trmuerr = [0.154, 0.142, 0.039, 0.079, 0.080, 0.098, 0.056]
trmlerr = [0.143, 0.131, 0.035, 0.075, 0.078, 0.095, 0.049]
trmerr = [trmlerr, trmuerr]



plt.errorbar(trsa, trm, yerr = trmerr, color = 'red', fmt='o', label = 'TRAPPIST system')
plt.legend()
plt.show()
