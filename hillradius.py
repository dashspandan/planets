# Calculate Hill radius and consequently the orbits of planets and returen list

import numpy as np

def findorbits(m_star, m_planet_isolation, disk_edge, iceline):

	m_star = m_star #solar masses
	m_star_kg = np.multiply(m_star, 1.989e30)
	m_planet_isolation = m_planet_isolation # earth masses
	m_planet_kg = np.multiply(m_planet_isolation, 5.972e24)
	disk_edge = disk_edge # in au
	factor = np.power(np.divide(m_planet_kg, 3.*m_star_kg), 1./3.)

	b = [disk_edge]
	d = disk_edge
	while d <= iceline:
		c = 5.*d*factor
		d = d+c
		b.append(d)

	return b

