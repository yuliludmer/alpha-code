import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from scipy import interpolate

gamma = 1.5
theta_obs = np.radians(150)
theta_gamma = np.arcsin(1/gamma)
opening_angle = 30
dot_num = 800 

t_0 = 500*24*60*60 #500 days
c = constants.c 

azimuthal_start = np.radians(0)
azimuthal_end = np.radians(360)
polar_start = np.radians(90 - opening_angle) 
polar_end = np.radians(90 + opening_angle) 

phi = np.linspace(azimuthal_start, azimuthal_end, dot_num) 
theta = np.linspace(polar_start, polar_end, dot_num) 

shape = (dot_num, dot_num) 
cos_zeta = np.zeros(shape) 
sin_zeta = np.zeros(shape)
zeta = np.zeros(shape)

for i in range(len(phi)): 
  for k in range(len(theta)): 
    cos_zeta[i,k] = (np.sin(-theta_obs)*np.sin(theta[k])*np.sin(phi[i]))+(np.cos(-theta_obs)*np.cos(theta[k]))
    sin_zeta[i,k] = np.sqrt(1-(cos_zeta[i, k])**2)
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

R = np.sqrt(1-(1/gamma**2))*c*t_0 
r = R*sin_zeta_beaming
unfilterd_r = R*sin_zeta

sin_psi = np.zeros(shape)
cos_psi = np.zeros(shape)
x = np.zeros(shape)
y = np.zeros(shape)

for i in range(len(x)): 
  for k in range(len(y)): 
    sin_psi[i,k] = np.round(-((np.cos(theta[k])-np.cos(-theta_obs)*(cos_zeta_beaming[i,k])))/((sin_zeta_beaming[i,k])*np.sin(-theta_obs)),10)
    cos_psi[i,k] = np.sqrt(1-(sin_psi[i,k])**2)
    x[i,k] = r[i,k]*cos_psi[i,k]
    y[i,k] = r[i,k]*sin_psi[i,k]

unfilterd_sin_psi = np.zeros(shape)
unfilterd_cos_psi = np.zeros(shape)
unfilterd_x = np.zeros(shape)
unfilterd_y = np.zeros(shape)
x_nonan = np.zeros(shape)
y_nonan = np.zeros(shape)

for i in range(len(x)): 
  for k in range(len(y)):
    unfilterd_sin_psi[i,k] = np.round(-(np.cos(theta[k])-np.cos(-theta_obs)*cos_zeta[i,k])/(sin_zeta[i,k]*np.sin(-theta_obs)),10)  
    unfilterd_cos_psi[i,k] = np.sqrt(1-(unfilterd_sin_psi[i,k])**2)
    unfilterd_x[i,k] = unfilterd_r[i,k]*unfilterd_cos_psi[i,k]
    unfilterd_y[i,k] = unfilterd_r[i,k]*unfilterd_sin_psi[i,k]
    if(~np.isnan(x[i,k])):
      x_nonan[i,k] = x[i,k]
      y_nonan[i,k] = y[i,k]

x_vec = np.linspace(np.min(x_nonan), np.max(x_nonan), dot_num)
y_vec = np.linspace(np.min(y_nonan), np.max(y_nonan), dot_num)

points = np.array((unfilterd_x.flatten(), unfilterd_y.flatten())).T
values = f.flatten().T.astype(float) 
(x_vec_mat,y_vec_mat) = np.meshgrid(x_vec,y_vec)

z_mat = interpolate.griddata(points,values,(x_vec_mat,y_vec_mat))

z_mat[np.isnan(z_mat)] = 0
z_mat[z_mat > 0] = 1

integral_y = np.trapz(z_mat, y_vec, axis=0)
integral_x = np.trapz(integral_y, x_vec, axis=0)
integral = integral_x*2

plt.scatter(x, y, color = 'blue')
plt.scatter(-x, y, color = 'blue')
plt.title(f"2d integral of ejected mass: {integral}", fontsize = 11)
plt.suptitle(f"gamma - {gamma} theta obs - 150 opening angle - {opening_angle}", fontsize = 9)
plt.xlabel("x - axis")
plt.ylabel("y - axis")
plt.gca().set_aspect('equal')

plt.show()