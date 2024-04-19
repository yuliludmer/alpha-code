import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
import math 

gamma = 3
theta_obs = np.radians(60) #check in 90 and 0 (supposed to come out blank)
theta_gamma = np.arcsin(1/gamma) #regular arcsin formula
opening_angle = 30 #in degrees, up for change
dot_num = 800 #number of dotes planted in phi and theta
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
zeta = np.zeros(shape)
for i in range(len(phi)): #run on every phi dot
  for k in range(len(theta)): #run on every theta dot
    cos_zeta[i,k] = (np.sin(-theta_obs)*np.sin(theta[k])*np.sin(phi[i]))+(np.cos(-theta_obs)*np.cos(theta[k]))
    sin_zeta[i,k] = np.sqrt(1-(cos_zeta[i, k])**2)
    if ((i+1)>=(dot_num/4) and (i+1)<=(dot_num*1.5/2)):
      sin_zeta[i,k] = -sin_zeta[i,k]
    if (i+1)>=(dot_num/2) and ((i+1)<=(dot_num)):
      cos_zeta[i,k] = -cos_zeta[i,k]
        
mask = sin_zeta >= np.sin(theta_gamma)

sin_zeta_beaming = np.copy(sin_zeta)
cos_zeta_beaming = np.copy(cos_zeta)

sin_zeta_beaming = np.ma.masked_array(sin_zeta, mask)
cos_zeta_beaming = np.ma.masked_array(cos_zeta, mask)

sin_mask = sin_zeta <= -np.sin(theta_gamma)
sin_zeta_beaming = np.ma.masked_array(sin_zeta, sin_mask)

f = np.zeros(shape)
f = np.where(sin_zeta < np.sin(theta_gamma), 1, 0)

R = np.sqrt(1-(1/gamma**2))*c*t_0 #R formula
r = R*sin_zeta_beaming
unfilterd_r = R*sin_zeta

sin_psi = np.zeros(np.shape(cos_zeta))
cos_psi = np.zeros(np.shape(cos_zeta))

x = np.zeros(shape)
y = np.zeros(shape)

for i in range(len(x)): #run on every phi dot
  for k in range(len(y)): #run on evry theta dot
    sin_psi[i,k] = -(np.cos(theta[k])-np.cos(-theta_obs)*(cos_zeta_beaming[i,k]))/((sin_zeta_beaming[i,k])*np.sin(-theta_obs))  #formula for psi
    cos_psi[i,k] = np.sqrt(1-(sin_psi[i,k])**2)
    #if(~np.isnan(cos_psi[i,k]) and ~np.isnan(sin_psi[i,k])):
    x[i,k] = r[i,k]*cos_psi[i,k]
    y[i,k] = r[i,k]*sin_psi[i,k]
'''
unfilterd_sin_psi = np.zeros(shape)
unfilterd_cos_psi = np.zeros(shape)
unfilterd_x = np.zeros(shape)
unfilterd_y = np.zeros(shape)

for i in range(len(x)): #run on every phi dot
  for k in range(len(y)): #run on evry theta dot
    unfilterd_sin_psi[i,k] = -(np.cos(theta[k])-np.cos(-theta_obs)*np.abs(cos_zeta[i,k]))/((sin_zeta[i,k])*np.sin(-theta_obs))  #formula for psi
    if(unfilterd_sin_psi[i,k]>1 or unfilterd_sin_psi[i,k]<-1):
      print("sin_psi:",unfilterd_sin_psi[i,k], "sin_zeta:",sin_zeta[i,k], "cos_zeta",cos_zeta[i,k],"theta:", theta[k],"phi:", phi[i])
    
    unfilterd_cos_psi[i,k] = np.sqrt(1-(unfilterd_sin_psi[i,k])**2)

    unfilterd_x[i,k] = unfilterd_r[i,k]*unfilterd_cos_psi[i,k]
    unfilterd_y[i,k] = unfilterd_r[i,k]*unfilterd_sin_psi[i,k]
print(unfilterd_cos_psi)

integral = np.trapz(np.trapz(f, unfilterd_y, axis=0), unfilterd_x, axis=0)

print(integral)
'''
plt.scatter(x,y)
plt.gca().set_aspect('equal')
plt.title("gamma = 1.5 ,theta_obs = 45")#+ str(int(np.degrees(theta_obs))))

plt.show()