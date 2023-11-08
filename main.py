import numpy as np
import matplotlib.pyplot as plt

gama = 3
theta_obs = np.radians(90) #check in 90 and 0 (supposed to come out blank)
theta_gama = np.arcsin(1/gama) #regular arcsin formula
opening_angle = 90 #in degrees, up for change
dot_num_phi = 10000 #number of dotes planted in phi and theta
dot_num_theta = 300

azimuthal_start = np.radians(0) #phi start
azimuthal_end = np.radians(360) #phi end
polar_start = np.radians(90 - opening_angle) #theta start
polar_end = np.radians(90 + opening_angle) #theta end

phi = np.linspace(azimuthal_start, azimuthal_end, dot_num_phi)
theta = np.linspace(polar_start, polar_end, dot_num_theta)

shape = (dot_num_phi, dot_num_theta)
zeta = np.zeros(shape)
zeta_beaming = np.ndarray.copy(zeta)
for i in range(len(phi)):
  for k in range(len(theta)):
    zeta[i,k] = np.arccos(np.cos(theta_obs)*np.cos(theta[k])+np.sin(theta_obs)*np.sin(theta[k])*np.sin(phi[i]))
    zeta_beaming[i,k] = zeta[i,k]
    if (zeta_beaming[i,k] > theta_gama):
      zeta_beaming[i,k] = 0

'''
#zeta = np.arccos(np.cos(theta_obs)*np.cos(theta)+np.sin(theta_obs)*np.sin(theta)*np.sin(phi)) #formula num 1 (govreen-segal) - changed cos(phi) to sin(phi) + why?
zeta_beaming = np.copy(zeta)

beaming = zeta > theta_gama
zeta_beaming[beaming] = 0 #zeroing each item smaller then arcsin.1/gama
'''

print(zeta)