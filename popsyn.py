#simulate planet distribution given stellar and disk parameters iteration 1 with just random selections

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import random as rn
from coeffsstandard import getcoeffdisk
from buildplanet import buildplanet
from hillradius import findorbits


mass_list = np.linspace(0.08, 0.2, 1000).tolist()
age_list = np.linspace(1e6, 1e7, 10000).tolist()
sim = 1
while sim <= 1000:

	p = 1. #surface density profile factor, must be less than 2
	z_o = 0.02 # metallicity of disk
	z = 0.5   # silicate to pebble ratio (From Ormel and Yamila, personal communication)
	m_star = rn.choice(mass_list) # star mass in solar masses 	
	m_disk = 0.04*m_star #disk mass in stellar masses, correlation of total disk to star is about 2-6% (Andrews 2013)
	#m_disk = np.exp(1.8*np.log(m_star)-2.6) #disk mass in stellar masses (Lupus 2016)
	M_dot_accre = 1e-10 # in solar masses per year
	r_out = 200. #outer radius of disk (To make sufficient pebbles)
	h_a = 0.03 #disk aspect ratio, assumed constant for the inner disk (Ormel 2017)
	disk_edge = 0.0102*np.power(0.08/m_star, 1./7.) # in au (from Ormel 2017)
	m_planet_isolation = m_star*np.power(h_a, 3.)*333060.402 #convert from solar to earth masses
	iceline = 0.125*(m_star/0.1)

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
