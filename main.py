import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
import math 

gamma = 3
theta_obs = np.radians(90) #check in 90 and 0 (supposed to come out blank)
theta_gamma = np.arcsin(1/gamma) #regular arcsin formula
opening_angle = 25 #in degrees, up for change
dot_num = 799 #number of dotes planted in phi and theta
t_0 = 500*24*60*60 #500 days
c = constants.c #speed of light

azimuthal_start = np.radians(0) #phi start
azimuthal_end = np.radians(360) #phi end
polar_start = np.radians(90 - opening_angle) #theta start
polar_end = np.radians(90 + opening_angle) #theta end

phi = np.linspace(azimuthal_start, azimuthal_end, dot_num) #making the dotes
theta = np.linspace(polar_start, polar_end, dot_num) #making the dotes

shape = (dot_num, dot_num) #size for the 2 dimensional arrays
cos_zeta = np.zeros(shape) #base values zeta
sin_zeta = np.zeros(shape)
#zeta_beaming = np.copy(zeta) #same as zeta
for i in range(len(phi)): #run on every phi dot
  for k in range(len(theta)): #run on evry theta dot
    cos_zeta[i,k] = np.cos(theta_obs)*np.cos(theta[k])+np.sin(theta_obs)*np.sin(theta[k])*np.sin(phi[i]) #zeta formula (changed)
    sin_zeta[i,k] = np.sqrt(1-(cos_zeta[i, k])**2)
    #if (zeta[i,k]<0):
     #print(zeta[i,k])

'''
sin_mask = sin_zeta >= np.sin(theta_gamma) #ask gilad
cos_mask = cos_zeta >= np.cos(theta_gamma)
sin_zeta_beaming = np.copy(sin_zeta)
cos_zeta_beaming = np.copy(cos_zeta)
sin_zeta_beaming = np.ma.masked_array(sin_zeta, sin_mask)
cos_zeta_beaming = np.ma.masked_array(cos_zeta, cos_mask)
#print(zeta_beaming<0)
'''
R = np.sqrt(1-(1/gamma**2))*c*t_0 #R formula
r = R*sin_zeta

sin_psi = np.zeros(np.shape(cos_zeta))
cos_psi = np.zeros(np.shape(cos_zeta))

#highest = -999
#lowest = 999

x = np.zeros(shape)
y = np.zeros(shape)
for i in range(len(x)): #run on every phi dot
  for k in range(len(y)): #run on evry theta dot
    sin_psi[i,k] = -(np.cos(theta[k])-np.cos(theta_obs)*cos_zeta[i,k])/(sin_zeta[i,k])*np.sin(theta_obs) #formula for psi 
    cos_psi[i,k] = np.sqrt(1-(sin_psi[i,k])**2)
    x[i,k] = r[i,k]*cos_psi[i,k]
    y[i,k] = r[i,k]*sin_psi[i,k]
    '''
    if(x[i,k]<0):
      print(x[i,k])
    if(math.isnan(psi[i,k])):
      continue
    else:
      print(psi[i,k])
      if(psi[i,k]>highest):
        highest = psi[i,k]
      if(psi[i,k]<lowest):
        lowest = psi[i,k]

print(highest, lowest)
'''
'''
x1 = R*np.sin(zeta)*np.cos(phi)*np.sin(theta)
y1 = R*np.sin(zeta)*np.sin(phi)*np.sin(theta)
plt.scatter(x1.flatten(),y1.flatten())
plt.show()
'''
plt.scatter(x,y)
plt.show()