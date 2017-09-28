import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as lin
from scipy.integrate import quad 

### Original motivation was from ECEN451 class - Dr.Huff 


### constant

epsilon = 8.854e-12
mu = 4*np.pi*1e-7
c = 3e8
f = 300e6
wavelength = c/f
N = 101 
a = 0.005 * wavelength 
k = 2 * np.pi / wavelength 
L = 0.47 * wavelength 
dz = L/N
feedpt = int(round(N/2))
omega = 2*np.pi*f 
zscale = 1j*dz/(omega * epsilon)
waveimp = np.sqrt(mu/epsilon)
zscale2 = 1j * dz * waveimp /k


### initialize matrices
v = np.zeros(N)
i = np.zeros(N)
zmn = np.zeros((N,N),complex)
pos = np.zeros(N)

### define position matrix
for i in range(N):
	pos[i] = (dz/2) + ((i-1)*dz) - (L/2)

### define integral 

for m in range(N):
	for n in range(N):
		zm = pos[m]
		zn = pos[n]
		
		def G_Re(z):
			R = np.sqrt(a**2 + (zm-z)**2)
			return np.real(np.exp(-1j * k * R)/(4*np.pi*R))
		
		def G_Im(z):
			R = np.sqrt(a**2 + (zm-z)**2)
			return np.imag(np.exp(-1j*k*R)/(4*np.pi*R))

		z_Re = quad(G_Re, zn - (dz/2) , zn + (dz/2))[0] 
		z_Im = quad(G_Im, zn - (dz/2) , zn + (dz/2))[0]

		R_upp = np.sqrt(a**2 + (zm - (zn + (dz/2)))**2)
		z_upp = (zm - (zn+ (dz/2))) * np.exp(-1j*k*R_upp) / (4*np.pi*R_upp **3) * ((1j*k*R_upp) + 1)

		R_low = np.sqrt(a**2 + (zm - (zn - (dz/2)))**2)
		z_low = (zm - (zn - (dz/2))) * np.exp(-1j*k*R_low) / (4*np.pi*R_low**3) * ((1j*k*R_low) + 1)

		zmn[m,n] = (k**2 * (z_Re + (1j*z_Im))) + (z_upp - z_low)

### scale impedance
Zmn = np.dot(zmn, zscale)
#print('Zmn = ',Zmn)


### define voltage matrix - Delta-gap source model
v[feedpt] = 1.
#print('voltage = ',v)


### define current matrix
ymn = lin.inv(Zmn)
i = np.dot(ymn, v)

### Find the input impedance
Zin = v[feedpt] / i[feedpt]
print('Zin = ',Zin)

### plot current distribution 
plt.figure()
plt.title('Current Distribution along thin wire dipole antenna')
plt.xlabel('zm [m]')
plt.ylabel('Current [A]')
plt.plot(pos,abs(i))
plt.show()	


### initialize electric field matrix
theta = np.linspace( -np.pi ,np.pi , 100)
phi = np.linspace(-np.pi,np.pi,100)
E_theta_real = np.zeros(np.size(theta))
E_theta_imag = np.zeros(np.size(theta))
E_theta_realXY = np.zeros(np.size(phi))
E_theta_imagXY = np.zeros(np.size(phi))

### define space factor for XZ and YZ planes 
for m in range(np.size(theta)):
	for n in range(np.size(pos)):
		E_theta_real[m] += np.real(i[n]*np.sin(theta[m])*np.exp(1j*k*pos[n]*np.cos(theta[m]))*dz)
		E_theta_imag[m] += np.imag(i[n]*np.sin(theta[m])*np.exp(1j*k*pos[n]*np.cos(theta[m]))*dz)

### define space factor for XY plane 
for m in range(np.size(phi)):
	for n in range(np.size(pos)):
		E_theta_realXY[m] += np.real(i[n]*np.sin(k*L/2) * np.exp(1j*k*pos[n]*np.cos(k*L/2)) * dz)
		E_theta_imagXY[m] += np.imag(i[n]*np.sin(k*L/2) * np.exp(1j*k*pos[n]*np.cos(k*L/2)) * dz)


### plane information
### XZ plane = (phi=0 plane)
### YZ plane = (phi=90 plane)
### XY plane = (theta=90 plane)
### Special note: "phi" angle can be found on the polar coordinate.

### Magnitude of field for XZ and YZ planes  
E_theta_mag = (E_theta_real**2) + (E_theta_imag**2) 
E_theta_norm = E_theta_mag / np.max(E_theta_mag)


### Magnitude of field for XY plane
E_theta_magXY = (E_theta_realXY**2) + (E_theta_imagXY**2)
E_theta_normXY = (E_theta_magXY) / np.max(E_theta_magXY)


### dB version - recall the argument should be dimensionless for dB 
E_theta_dB = 10*np.log10(E_theta_norm)
E_thetaXY_dB = 10*np.log10(E_theta_normXY)

### magnitude plots
plt.figure(1)
fig1 = plt.subplot(221, polar = True)
plt.polar(theta , E_theta_mag)
plt.title(r'XZ plane - phi=0: E$\theta$')
fig1.set_theta_zero_location('N')

fig2 = plt.subplot(222, polar = True)
plt.polar(phi, E_theta_magXY)
plt.title(r'XY plane - theta=90: E$\theta$')
fig2.set_theta_zero_location('N')

fig3 = plt.subplot(223,polar = True)
plt.polar(theta, E_theta_mag)
plt.title(r'YZ plane - phi=90: E$\theta$')
fig3.set_theta_zero_location('N')

fig4 = plt.subplot(224)
plt.plot(pos,abs(i))
plt.title('zm vs. current: z-directed dipole')
plt.xlabel('zm')
plt.ylabel('Ampere')


### dB plots 
plt.figure(2)
fig1 = plt.subplot(221, polar = True)
plt.polar(theta, E_theta_dB)
plt.title(r'XZ plane - phi=0: E$\theta$ [dB]')
fig1.set_ylim(-40,10)
fig1.set_yticks([-40, -30, -20, -10, 0 ,10])
fig1.set_theta_zero_location('N')

fig2 = plt.subplot(222, polar = True)
plt.polar(phi, E_thetaXY_dB)
plt.title(r'XY plane - theta=90: E$\theta$ [dB]')
fig2.set_ylim(-40,10)
fig2.set_yticks([-40, -30, -20, -10, 0 ,10])
fig2.set_theta_zero_location('N')

fig3 = plt.subplot(223, polar = True)
plt.polar(theta, E_theta_dB)
plt.title(r'YZ plane - phi=90: E$\theta$ [dB]')
fig3.set_ylim(-40,10)
fig3.set_yticks([-40, -30, -20, -10, 0 ,10])
fig3.set_theta_zero_location('N')

fig4 = plt.subplot(224)
plt.plot(pos, abs(i))
plt.title('zm vs. current: z-directed dipole')
plt.xlabel('zm')
plt.ylabel('Ampere')

plt.show()


