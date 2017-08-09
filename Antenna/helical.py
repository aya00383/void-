import numpy as np
import matplotlib.pyplot as plt

c = 3e8
L = input('Length of Helical? [m]') # meter
N = float(input('Number of turns?'))
S = float(L/N) 
print('Spacing between each turn = ',S,' [meter]')
alpha = 13 * np.pi / 180.0  
D = S / np.pi / np.tan(alpha)
print('Diameter of Helix = ',D,' [meter]')
C = np.pi * D
print('C = ',C)
lamda = C*1
print('Wavelength = ',lamda,' [meter]')

f = c/lamda
print('center frequency of Helix = ',f/1.0e6,' [MHz]')
gnd_plane = 3 * lamda / 4.0 
print('Ground plane = ',gnd_plane, ' [meter]')

Rin = 140.0 * C /lamda
Rin_max = Rin + (Rin * 20/100)
Rin_min = Rin - (Rin * 20/100)
print('Rin = ',Rin)
print('Max Rin = ',Rin_max)
print('Min Rin = ',Rin_min)

D0 = 15 * N * C**2 * S / lamda**3
print('Directivity = ',D0)

L0 = np.sqrt(S**2 + C**2)

# for ordinary end-fire
p1 = L0/C/ (S/C + 1)

# for hansen-woodyard end-fire
p2 = L0 / C / ((S/C) + ((2*N) +1) / (2*N))

k = 2 * np.pi / C

def psy(p):
	return k * (S * np.cos(theta) - (L0/p))

def E_normalized(p):
	return np.sin(np.pi/2/N) * np.cos(np.sin(N/2 * psy(p)) / np.sin(psy(p)/2))

theta = np.linspace(-np.pi,np.pi,100)

G = 10.8 + 10*np.log10( (C/lamda)**2 * N * (S/lamda) )
print('Gain = ',G)

print('Far-Field Distance = ',2*L**2/lamda,' [m]')
'''
plt.figure()
ax1 = plt.subplot(211, polar = True)
ax1.polar(theta, E_normalized(p1)
plt.title('Normalized far-field pattern - Ordinary end-fire mode')

ax2 = plt.subplot(212, polar = True)
plt.polar(theta, E_normalized(p2))
plt.title('Normalized far-field pattern - Hansen-Woodyard end-fire mode')
plt.show()

'''
