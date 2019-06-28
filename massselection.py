import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import expon
import csv

def mass_selection_ormel():

	m_star_list = np.linspace(0.08, 0.2, 1000)

	dm = m_star_list[1] - m_star_list[0]

	IMF = np.power(m_star_list, -1.3) #From Chandler 2003

	#n = (1./0.3)*((1./(0.08**0.3)) - (1./(m_star**0.3)))

	number = IMF*dm


	loc, scale = expon.fit(number.tolist())

	t = expon.rvs(loc, scale, size = 1)

	if t[0] < number[0] and t[0] > number[len(number)-1]:
		return 0.04*np.power(np.divide(t[0], dm), -1./1.3)
	else:
		pass

def mass_selection_lupus():

	m_star_list = np.linspace(0.08, 0.2, 1000)

	dm = m_star_list[1] - m_star_list[0]

	IMF = np.power(m_star_list, -1.3) #From Chandler 2003

	#n = (1./0.3)*((1./(0.08**0.3)) - (1./(m_star**0.3)))

	number = IMF*dm


	loc, scale = expon.fit(number.tolist())

	t = expon.rvs(loc, scale, size = 1)

	if t[0] < number[0] and t[0] > number[len(number)-1]:
		return np.exp(1.8*np.log(np.power(np.divide(t[0], dm), -1./1.3))-2.6)
	else:
		pass

def mass_selection_pop_syn_lupus(a):
#m_disk_ormel = []
	m_disk_lupus = []


#while len(m_disk_ormel) < 40:
#	m = mass_selection_ormel()
#	if m is not None:
#		m_disk_ormel.append(m)
#	else:
#		pass

#l = [i for i in m_disk_ormel if i is not None]

#print len(l)

	while len(m_disk_lupus) < a:
		m = mass_selection_lupus()
		if m is not None:
			m_disk_lupus.append(m)
		else:
			pass

	star_mass = np.exp((np.log(m_disk_lupus) + 2.6)/1.8)

	return star_mass.tolist(), m_disk_lupus


def mass_selection_pop_syn_ormel(a):
#m_disk_ormel = []
	m_disk_ormel = []


#while len(m_disk_ormel) < 40:
#	m = mass_selection_ormel()
#	if m is not None:
#		m_disk_ormel.append(m)
#	else:
#		pass

#l = [i for i in m_disk_ormel if i is not None]

#print len(l)

	while len(m_disk_ormel) < a:
		m = mass_selection_ormel()
		if m is not None:
			m_disk_ormel.append(m)
		else:
			pass

	star_mass = np.divide(m_disk_ormel, 0.04)

	return star_mass.tolist(), m_disk_ormel

#k = [i for i in m_disk_lupus if i is not None]

#print len(k)

#m_disk_actual = []

#with open('disks.csv','r') as csvfile:
#    plots = csv.reader(csvfile, delimiter=',')
#    for row in plots:
#		if float(row[0]) > 0.08 and float(row[0]) < 0.2:
#			m_disk_actual.append(float(row[1]))
#		else:
#			pass

#bins = np.logspace(np.log10(1.e-6),np.log10(1.e-2), 10)

#plt.hist(l, bins = bins, label = 'Ormel', alpha = 0.4) # Ormel 2017
#plt.hist(k, bins = bins, label = 'Lupus', alpha = 0.6) # Lupus survey
#plt.hist(m_disk_actual, bins = bins, label = 'Data from Yamila', alpha = 0.7) # Actual masses 
#plt.xlabel('Mass of Disks in M$_{\odot}$', size = 18)
#plt.ylabel('Freq.', size = 18)
#plt.title('Number distribution of Disks')
#plt.xscale('log')
#plt.legend(loc = 'upper left')
#plt.show()
