# Planet mass-distance profile
from coeffsstandard import getcoeffdisk
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

list_m = [0.08, 0.1, 0.2] #list of masses to plot profiles

num = 0

while num < len(list_m):
	p = 1. #surface density profile factor, must be less than 2
	z_o = 0.02 # metallicity of disk
	z = 0.5   # silicate to pebble ratio
	m_star = list_m[num] # star mass in solar masses	
	m_disk = 0.04*m_star #disk mass in stellar masses
	M_dot_accre = 1e-10 # in solar masses per year
	r_out = 200. #outer radius of disk
	h_a = 0.03 #disk aspect ratio
	f = 0. #Start of internal time
	o = 3e6 #Lifetime of disk
	e = 0.0102*np.power(0.08/m_star, 1./7.) #inner disk radius
	iceline = 0.125*(m_star/0.1) #From Ormel 2017
	m_planet_isolation = m_star*np.power(h_a, 3.)*333060.402 #convert from solar to earth masses

	m_dot_2D_coeff, m_dot_3D_coeff, y, a_dot_coeff = getcoeffdisk(p, z_o, z, m_star, m_disk, M_dot_accre, r_out, h_a)

	a = []
	b = []
	c = []


	w = 0. #initiate internal time (time in planet's frame of reference)
	r = iceline + 0.02 #Outer radius with midplane pebble-gas ratio 1
	
	#Exterior to iceline

	j = lambda w: a_dot_coeff*np.power(m_dot_2D_coeff*(np.power(w,0.67))+np.power(1e-4,0.33),3.) #For migration
	while r > iceline:
		h,k = integrate.quad(j, 0., w)
		r = np.power(np.add(h,np.power(iceline + 0.02,y)), 1./y)
		c.append(w)
		a.append(w+f)
		b.append(r)
		w = w+1e2

	M_p = np.power(m_dot_2D_coeff*(np.power(c,0.67))+np.power(1e-4,0.33),3.)  #Mass-behaviour
	Mp = M_p.tolist()
	tt1 = w
	m1 = Mp[-1]

	print 'Planetesimal crosses iceline with mass %0.3f M_earth after %d years.'%(m1,a[-1])
	
	#Interior to iceline	
	z = iceline
	i = lambda w: a_dot_coeff*m1*np.exp(m_dot_3D_coeff*(np.power(w,0.67)-np.power(tt1,0.67)))

	m2 = Mp[-1]
	while m2 < m_planet_isolation and a[-1] < o:
		m = m1*np.exp(m_dot_3D_coeff*(np.power(w,0.67)-np.power(tt1,0.67)))
		a.append(w+f)
		c.append(w)
		if z >= e:
			g,s = integrate.quad(i, tt1, w)
			z = np.power(np.power(iceline,y)+g, 1./y)
			b.append(z)
			Mp.append(m)
			m2 = m
			w = w+1e2
		else:
			b.append(e)
			Mp.append(m)
			m2 = m
			w = w+1e2

	print 'Simulation stops after %d years with planet mass %0.3f M_earth.'%(a[-1],Mp[-1])


	plt.plot(b, Mp, label = str(m_star)+' $M_{\odot}$')
	plt.title('Growth-Migration behaviour of a sample planet around stars', size = 28)
	plt.xlabel('Semi-major axis in AU', size = 22)
	plt.ylabel('Mass of planet in $M_{\oplus}$', size = 22)
	plt.xscale('log')
	plt.yscale('log')
	#plt.grid()
	plt.tick_params(labelsize = 18)
	plt.legend(prop={'size': 18})

	num = num + 1


plt.show()	
