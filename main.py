import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from scipy import interpolate
import math 

gamma = 3
theta_obs = np.radians(135)
theta_gamma = np.arcsin(1/gamma)
opening_angle = 60
dot_num = 800 
t_0 = 500*24*60*60 #500 days
c = constants.c 

print(theta_gamma)

azimuthal_start = np.radians(0) #phi start
azimuthal_end = np.radians(360) #phi end
polar_start = np.radians(90 - opening_angle) #theta start
polar_end = np.radians(90 + opening_angle) #theta end

phi = np.linspace(azimuthal_start, azimuthal_end, dot_num) 
theta = np.linspace(polar_start, polar_end, dot_num) 

shape = (dot_num, dot_num) 
cos_zeta = np.zeros(shape) 
sin_zeta = np.zeros(shape)
zeta = np.zeros(shape)

for i in range(len(phi)): #run on every phi dot
  for k in range(len(theta)): #run on every theta dot
    cos_zeta[i,k] = (np.sin(-theta_obs)*np.sin(theta[k])*np.sin(phi[i]))+(np.cos(-theta_obs)*np.cos(theta[k]))
    sin_zeta[i,k] = np.sqrt(1-(cos_zeta[i, k])**2)
    if ((i+1)>=(dot_num/4) and (i+1)<=(dot_num*1.5/2)):
      sin_zeta[i,k] = -sin_zeta[i,k]
    if (cos_zeta[i,k]>=0):
     zeta[i,k] = np.arccos(cos_zeta[i,k])
    else:
      zeta[i,k] = (np.pi-np.arccos(np.abs(cos_zeta[i,k])))

mask = np.where(zeta >= theta_gamma, 1, 0)

sin_zeta_beaming = np.copy(sin_zeta)
cos_zeta_beaming = np.copy(cos_zeta)
sin_zeta_beaming = np.where((mask==1), np.nan, sin_zeta)
cos_zeta_beaming = np.where((mask==1), np.nan, cos_zeta)

f = np.zeros(shape)
f = np.where(zeta < theta_gamma, 1, 0)

print(np.sum(f))

R = np.sqrt(1-(1/gamma**2))*c*t_0 
r = R*sin_zeta_beaming
unfilterd_r = R*sin_zeta

sin_psi = np.zeros(shape)
cos_psi = np.zeros(shape)
x = np.zeros(shape)
y = np.zeros(shape)

for i in range(len(x)): 
  for k in range(len(y)): 
    sin_psi[i,k] = (-(np.cos(theta[k])-np.cos(-theta_obs)*(cos_zeta_beaming[i,k])))/((sin_zeta_beaming[i,k])*np.sin(-theta_obs))
    cos_psi[i,k] = np.sqrt(1-(sin_psi[i,k])**2)
    x[i,k] = r[i,k]*cos_psi[i,k]
    y[i,k] = r[i,k]*sin_psi[i,k]

unfilterd_sin_psi = np.zeros(shape)
unfilterd_cos_psi = np.zeros(shape)
unfilterd_x = np.zeros(shape)
unfilterd_y = np.zeros(shape)

for i in range(len(x)): 
  for k in range(len(y)):
    unfilterd_sin_psi[i,k] = -(np.cos(theta[k])-np.cos(-theta_obs)*cos_zeta[i,k])/(sin_zeta[i,k]*np.sin(-theta_obs))  
    unfilterd_cos_psi[i,k] = np.sqrt(1-(unfilterd_sin_psi[i,k])**2)
    unfilterd_x[i,k] = unfilterd_r[i,k]*unfilterd_cos_psi[i,k]
    unfilterd_y[i,k] = unfilterd_r[i,k]*unfilterd_sin_psi[i,k]
'''
min_nan_x = np.max(unfilterd_x)
max_nan_x = np.min(unfilterd_x)
min_nan_y = np.max(unfilterd_y)
max_nan_y = np.min(unfilterd_y)

for i in range(len(x)): 
  for k in range(len(y)):
    if(np.isnan(x[i,k]) and unfilterd_x[i,k]<min_nan_x):
      min_nan_x = unfilterd_x[i,k]
    if(np.isnan(x[i,k]) and unfilterd_x[i,k]>max_nan_x):
      min_nan_x = unfilterd_x[i,k]
    if(np.isnan(y[i,k]) and unfilterd_y[i,k]<min_nan_y):
      min_nan_y = unfilterd_y[i,k]
    if(np.isnan(y[i,k]) and unfilterd_y[i,k]>max_nan_y):
      min_nan_y = unfilterd_y[i,k]
'''
x_vec = np.linspace(np.min(unfilterd_x), np.max(unfilterd_x), dot_num)
y_vec = np.linspace(np.min(unfilterd_y), np.max(unfilterd_y), dot_num)

points = np.array((unfilterd_x.flatten(), unfilterd_y.flatten())).T
values = np.ones(dot_num**2) #f.flatten()
(x_vec_mat,y_vec_mat) = np.meshgrid(x_vec,y_vec)

z_mat = interpolate.griddata((unfilterd_x.flatten(), unfilterd_y.flatten()),values,(x_vec_mat,y_vec_mat))
zmat_clean = np.nan_to_num(z_mat, copy=True, nan=0.0, posinf=None, neginf=None)
'''
for i in range(dot_num):
  for k in range(dot_num):
    if(x_vec_mat[i,k]>min_nan_x and x_vec_mat[i,k]<max_nan_x and y_vec_mat[i,k]>min_nan_y and y_vec_mat[i,k]<max_nan_y):
      zmat_clean[i,k] = 0
'''
for i in range(dot_num):
  for k in range(dot_num):
    if(f[i,k] == 0):
      zmat_clean[i,k] = 0

integral_y = np.trapz(zmat_clean, y_vec, axis=0)
integral_x = np.trapz(integral_y, x_vec, axis=0)
integral = np.sum(integral_x)

print(integral)

plt.scatter(x,y)
plt.gca().set_aspect('equal')

plt.show()