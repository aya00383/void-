import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as lin
from scipy import integrate


#--- varibles
epsilon = 8.854e-12     # F/m
N = input('Number of elements that you want to use in approximation?')
a = float(input('radius of wire? [m]'))
l = float(input('length of wire? [m]'))
dyn = l/N              # each element yn
dym = l/N		# each element ym


#--- initializing the matrix
V = np.zeros(N)
I = np.zeros(N)
ym = np.zeros(N)
yn = np.zeros(N+1)
Z = np.zeros((N,N))

#--- define voltage 
for m in range(N):
	V[m]= 4*np.pi*epsilon
	ym[m]= dym + (dym*m) -(dym/2) 
	yn[m]= dyn + (dyn*m) 
#--- equation
for m in range(N):
	for n in range(N):
		Z[m,n]= integrate.quad(lambda y: 1./(np.sqrt((a**2)+(ym[m]-y)**2)) ,yn[n-1],yn[n])[0]		

Z_inv = lin.inv(Z)
I = np.dot(Z_inv, V)

def pulse(x, pos, h, w):
	return h * (x>= pos) * (x <= (pos + w))
x = np.linspace(0,l,1000)

print('I = ',I)
print('Z = ',Z)
plt.figure()
plt.title('Charge Distribution')
for t in range(N):
	plt.plot(x,pulse(x,yn[t-1],I[t],dyn),color='b')
plt.xlabel('y-axis')
plt.ylabel('Amplitude')
plt.show()

