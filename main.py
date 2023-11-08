import numpy as np
import matplotlib.pyplot as plt

gama = 3
theta_obs = np.radians(90) #check in 90 and 0 (supposed to come out blank)
theta_gama = np.arcsin(1/gama) #regular arcsin formula
opening_angle = 20 #in degrees, up for change
dot_num = 300 #number of dotes planted in phi and theta

azimuthal_start = np.radians(0) #phi start
azimuthal_end = np.radians(360) #phi end
polar_start = np.radians(90 - opening_angle) #theta start
polar_end = np.radians(90 + opening_angle) #theta end

phi = np.linspace(azimuthal_start, azimuthal_end, dot_num)
theta = np.linspace(polar_start, polar_end, dot_num)

shape = (dot_num, dot_num)
zeta = np.zeros(shape)
zeta_beaming = np.copy(zeta)
for i in range(len(phi)):
  for k in range(len(theta)):
    zeta[i,k] = np.arccos(np.cos(theta_obs)*np.cos(theta[k])+np.sin(theta_obs)*np.sin(theta[k])*np.sin(phi[i]))
    zeta_beaming[i,k] = zeta[i,k]
    if (zeta_beaming[i,k] > theta_gama):
      zeta_beaming[i,k] = 0
    #if (zeta_beaming[i,k] != 0):
      #print(zeta_beaming[i,k], i, k, phi[i], theta[k])

print(zeta_beaming)