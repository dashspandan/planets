# Set coeffs for accretion interior and exterior to iceline for given disk parameters

import numpy as np

def getcoeffdisk(p, z_o, z, m_star, m_disk, M_dot_accre, r_out, h_a):

	p = p #surface density profile factor, must be less than 2
	z_o = z_o # metallicity of disk
	z = z   # silicate to pebble ratio
	m_star = m_star # star mass in solar masses
	m_star_kg = np.multiply(m_star, 1.989e30) # mass of star
	m_disk = m_disk #in solar masses
	m_disk_kg = np.multiply(m_disk, 1.989e30) # in kg
	M_dot_accre = M_dot_accre # in solar masses per year
	M_dot_accre_kg = np.multiply(M_dot_accre, 1.989e30) # in kg per yr
	r_out = r_out # in AU
	r_out_m = np.multiply(r_out, 1.496e11) # in metres
	G = np.divide(6.67e-11, np.power(3.17e-8,2.)) #gravitational constant in m3kg-1yr-2
	G1 = np.multiply(G, np.power(6.68e-12, 3.)) #gravitational constant in AU3kg-1yr-2
	h = h_a #disk aspect ratio

	# Pebble to gas flux in form of f_p/g = f_coeff_norm(t/1e5)^(-1/3)

	f_coeff = np.divide(2.*m_disk*np.power(z_o, 5./3.)*np.power(np.multiply(G,m_star_kg),1./3.), 3.*M_dot_accre*r_out_m*np.power(10., 2./3.))
	f_coeff_norm = np.divide(f_coeff,np.power(1e5, 1./3.))

	print 'Pebble to gas flux for the disk is f_p/g = %.2f*(t/1e5 yr)^(-1/3)' %(f_coeff_norm) 


	#Exterior to iceline

	#The accretion is of form m_pl_dot_2D = e_pa_2D*f_p/g*M_dot_accre_kg = m_pl*(2D_coeff*t^(2/3) + M_o^(1/3))^3
	#e_pa_2D = e_pa_2D_coeff*M_pl^(2/3)

	e_pa_2D_coeff = 184.61*np.power(m_star_kg, -2./3.)
	m_dot_2D_coeff = np.divide(0.5*np.multiply(np.multiply(e_pa_2D_coeff, f_coeff), M_dot_accre_kg), np.power(6.e24, 1./3.)) #in earth masses

	print 'm_pl_2D = (%.2E*t^(2/3) + (1e^-4)^(1/3))^3' %m_dot_2D_coeff 


	#Interior to iceline

	#The accretion is of form m_pl_dot_3D = e_pa_3D*f_s/g*M_dot_accre_kg
	#e_pa_3D = e_pa_3D_coeff*M_pl

	e_pa_3D_coeff = np.divide(0.07, np.multiply(m_star_kg, 1e-5))
	m_dot_3D_coeff = 1.5*np.multiply(np.multiply(np.multiply(f_coeff, np.power(z, 5./3.)), e_pa_3D_coeff), M_dot_accre_kg)

	print 'm_pl_3D = m_pl_core_ice*exp(%.2E*(t^(2/3) - t_icecross^(2/3)))' %m_dot_3D_coeff


	#Gas surface density of form of sigma_g = sigma_g_coeff*r^(-p)

	sigma_g_coeff = (2.-p)*np.divide(m_disk_kg, 2.*np.pi*np.power(r_out, 2.-p))

	print 'sigma_g = %.2E*a^(-%d) where a in AU, others in SI' %(sigma_g_coeff, p)

	#Type I migration a_dot_coeff

	y = (1.5+p-2.) #exp factor for migration integration, needed in the code for growth
	a_dot_coeff = -y*np.divide(6.e24*sigma_g_coeff*2.8*np.power(G1, 0.5), np.power(m_star_kg, 1.5)*np.power(h, 2.))

	print 'a^(%.2f)-a_o^(%.2f) = %.2E*integration of M(in earth masses)dt' %(y, y, a_dot_coeff)

	return m_dot_2D_coeff, m_dot_3D_coeff, y, a_dot_coeff


