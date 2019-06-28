# Grow the planet

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

def buildplanet(f, count, o, m_dot_2D_coeff, m_dot_3D_coeff, y, a_dot_coeff, m_planet_isolation, orbits, iceline): #f is total time, count is the number of planets, o is disk lifetime
	
	a = []
	b = []
	c = []
	q = count
	d = np.linspace(1,100,100).tolist()  # a sample list of number of planets that can be formed
#	x = np.linspace(0.01,0.03,100) # sample spacing between planets
	x = orbits
	e = x[d.index(q)]
	file1 = open("simulplanetormelnoimf2.txt","a")
	
	w = 0. #initiate internal time (time in planet's frame of reference)
#	r = 0.12 #Outer radius with midplane pebble-gas ratio 1

	r = iceline + 0.02 #assume width of iceline is 0.02 au 
	
#	j = lambda t: (a_dot_coeff)*np.power(1.76e-4*(np.power(t,0.67))+np.power(1e-4,0.33),3.) #Just the integration to find mass
	j = lambda w: a_dot_coeff*np.power(m_dot_2D_coeff*(np.power(w,0.67))+np.power(1e-4,0.33),3.)
	while r>iceline:
		h,k = integrate.quad(j, 0., w)
#		r = np.power(np.add(h,np.power(0.12,y)), 1./y)
		r = np.power(np.add(h,np.power(iceline+0.02,y)), 1./y)
		c.append(w)
		a.append(w+f)
		b.append(r)
		w = w+1e2

	M_p = np.power(m_dot_2D_coeff*(np.power(c,0.67))+np.power(1e-4,0.33),3.)
	Mp = M_p.tolist()
	tt1 = w
	m1 = Mp[-1]

	print 'Planetesimal crosses iceline with mass %0.3f M_earth after %d years.'%(m1,a[-1])
	
	
	z = iceline
	i = lambda w: a_dot_coeff*m1*np.exp(m_dot_3D_coeff*(np.power(w,0.67)-np.power(tt1,0.67)))

	m2 = Mp[-1]
	while m2 < m_planet_isolation and a[-1] < o:
		m = m1*np.exp(m_dot_3D_coeff*(np.power(w,0.67)-np.power(tt1,0.67)))
		a.append(w+f)
		c.append(w)
		if z >= e:
			g,s = integrate.quad(i, tt1, w)
#			z = np.power(np.power(0.1,y)+g, 1./y)
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
#	else:
#		print 'Disk dispersed before isolation mass. Final planet mass is %0.3f M_earth'%(Mp[-1])

	print 'Simulation stops after %d years with planet mass %0.3f M_earth.'%(a[-1],Mp[-1])
#	plt.scatter(b[-1],Mp[-1], label = str(q)+'th planet')
#	plt.figure()
#	ax = plt.gca()
#	return b[-1], Mp[-1] # Return planet semi-major axis and mass values to save them separately
	if Mp[-1] < m_planet_isolation:
		file1.write(str(m1)+','+str(b[-1])+','+str(Mp[-1])+'\n')
	
		plt.scatter(b[-1],Mp[-1],color = 'blue')
		plt.xscale('log')
		plt.yscale('log')
		plt.title('Distribution of planets', size = 22)
		plt.xlabel('Semi-major axis in AU', size = 16)
		plt.ylabel('Mass of planet in $M_{\oplus}$', size = 16)
		plt.xlim(0.005, 1.)
#	plt.grid()
		plt.tick_params(labelsize=14)
#	plt.legend()


#	ax.scatter(b[-1],Mp[-1])
#	ax.set_title('Distribution of planets', size = 22)
#	ax.set_xlabel('Semi-major axis in AU', size = 16)
#	ax.set_ylabel('Mass of planet in $M_{\oplus}$', size = 16)
#	ax.set_xlim(0.005, 0.1)
#	plt.grid()
#	plt.tick_params(labelsize=14)
#	plt.legend()


#	q = count
	if a[-1] < o:
		q = q + 1
		buildplanet(a[-1], q, o, m_dot_2D_coeff, m_dot_3D_coeff, y, a_dot_coeff, m_planet_isolation, orbits, iceline)
	else:
		print 'No more planets can be formed as disk dispersed.'

