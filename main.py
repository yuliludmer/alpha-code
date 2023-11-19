import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

gama = 3
theta_obs = np.radians(90) #check in 90 and 0 (supposed to come out blank)
theta_gama = np.arcsin(1/gama) #regular arcsin formula
opening_angle = 20 #in degrees, up for change
dot_num = 300 #number of dotes planted in phi and theta
t_0 = 500*24*60*60 #500 days
c = constants.c #speed of light

azimuthal_start = np.radians(0) #phi start
azimuthal_end = np.radians(360) #phi end
polar_start = np.radians(90 - opening_angle) #theta start
polar_end = np.radians(90 + opening_angle) #theta end

phi = np.linspace(azimuthal_start, azimuthal_end, dot_num) #making the dotes
theta = np.linspace(polar_start, polar_end, dot_num) #making the dotes

shape = (dot_num, dot_num) #size for the 2 dimensional arrays
psi = np.zeros(shape) #base values psi
zeta = np.zeros(shape) #base values zeta
zeta_beaming = np.copy(zeta) #same as zeta
for i in range(len(phi)): #run on every phi dot
  for k in range(len(theta)): #run on evry theta dot
    zeta[i,k] = np.arccos(np.cos(theta_obs)*np.cos(theta[k])+np.sin(theta_obs)*np.sin(theta[k])*np.sin(phi[i])) #zeta formula (changed )
    zeta_beaming[i,k] = zeta[i,k] #copy zeta value
    if (zeta_beaming[i,k] > theta_gama): #filter each zeta value (beaming)
      zeta_beaming[i,k] = 0 #zero any zeta_beaming that is bigger then theta_gama

    psi[i,k] = np.arcsin(-(np.cos(theta[k])-np.cos(theta_obs)*np.cos(zeta[i,k]))/(np.sin(zeta[i,k])*np.sin(theta_obs))) #formula for psi (used zeta and not zeta_beaming)

    #if (zeta_beaming[i,k] != 0):
      #print(zeta_beaming[i,k], i, k, phi[i], theta[k])

R = np.sqrt(1-(1/gama**2))*c*t_0 #R formula
r = R*np.sin(zeta_beaming)

x = np.zeros(shape)
y = np.zeros(shape)
for i in range(len(phi)): #run on every phi dot
  for k in range(len(theta)): #run on evry theta dot
    x[i,k] = r[i,k]*np.cos(psi[i,k])
    y[i,k] = r[i,k]*np.sin(psi[i,k])

plt.scatter(x, y)
plt.show()