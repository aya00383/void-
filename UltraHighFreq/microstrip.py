import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt

print("Welcome to Microstrip Line Calculator!")
q1 = input("(1) W/d is given or (2) Z0 is given? (3) Z0 and either W or d is given. Type the corresponding number.  ")

e_r = float(input('The relative dielectric constant?  '))


if q1 == 1:
	d = float(input('The thickness of microstrip line? [meters]  '))
	W = float(input('The width of microstrip line? [meters]  '))

	e_eff = ((e_r+1)/2.) + ((e_r-1)/2*(1/np.sqrt(1+(12*d/W))))
	print('W/d ratio = ',W/d)

	if W/d <1:
		Z0 = 60/np.sqrt(e_eff)*np.log((8*d/W)+(W/d/4))
		print('The characteristic impedance Z0 = ',Z0,' [ohm]')
	elif W/d >=1:
		Z0 = 120*np.pi/np.sqrt(e_eff)/( (W/d) + 1.393 + (0.667*np.log(W/d + 1.444)))

	print('The calculated characteristic impedance Z0 = ',Z0,' ohm')

elif q1 == 2:
	Z0 = float(input('The characteristic impedance? [ohm]  '))
	A = Z0/60 * np.sqrt((e_r+1)/2) + (e_r-1)/(e_r+1)*(0.23+ (0.11/e_r))
	B = 377*np.pi/2./Z0/np.sqrt(e_r)
	# the variable "ratio" refers to "W/d ratio"
	q2 = input("Would you want the W/d ratio to be bigger than 2? Type \'1'\ if yes. Otherwise, type \'2\'  ") 

	if q2 == 1:
		ratio = 2/np.pi*(B-1-np.log(2*B -1) + ((e_r-1)/2/e_r *(np.log(B-1)+0.39-(0.61/e_r))))
	elif q2 == 2:
		ratio = 8*np.exp(A)/(np.exp(2*A)-2)
	else:
		print('You typed something other than the given options. Program terminated.')

	print('The calculated W/d ratio = ',ratio)

elif q1 == 3:
	Z0 = float(input('The characteristic impedance? [ohm]  '))
	W = float(input('The width of microstrip line? [meters] (Type \'0\' if this value tbd)  '))
	d = float(input('The thickness of microstrip line? [meters] (Type\'0\' if this value tbd.)  '))
	if W == 0 and d == 0:
		print('You type zeros for both \'W\' and \'d\'. This should have been option #1.')
	#	break
	else:
		A = Z0/60 * np.sqrt((e_r+1)/2) + (e_r-1)/(e_r+1)*(0.23+ (0.11/e_r))
		B = 377*np.pi/2./Z0/np.sqrt(e_r)
		
		if d == 0:
			ratio1 = 2/np.pi*(B-1-np.log(2*B -1) + ((e_r-1)/2/e_r *(np.log(B-1)+0.39-(0.61/e_r))))
			ratio2 = 8*np.exp(A)/(np.exp(2*A)-2)
			if ratio1 < 2:
			#	print('The ratio is less than 2.')
				#print('This is ratio2 value (this should be greater than 2) = ', ratio2)
				d = W/ratio1 
				print('d = ',d,' [meters]')

			elif ratio2 > 2:
			#	print('The ratio is greater than 2.')
				#print('This is ratio1 value (this should be less than 2) = ', ratio1)
				d = W/ratio2
				print('d = ',d,' [meters]')

			else:
				print('Problem occurs. Try Again.')
		elif W == 0:
			ratio1 = 2/np.pi*(B-1-np.log(2*B -1) + ((e_r-1)/2/e_r *(np.log(B-1)+0.39-(0.61/e_r))))
			ratio2 = 8*np.exp(A)/(np.exp(2*A)-2)
			if ratio1 < 2:
			#	print('The ratio is less than 2.')
				#print('This is ratio2 value = ', ratio2)
				W = d*ratio1 
				print('W = ',W, ' [meters]')

			elif ratio2 > 2:
			#	print('The ratio is greater than 2.')
				#print('This is ratio1 value = ', ratio1)
				W = d * ratio2
				print('W = ',W,' [meters]')

			else:
				print('Problem occurs. Try Again.')
	e_eff = ((e_r+1)/2.) + ((e_r-1)/2*(1/np.sqrt(1+(12*d/W))))
	print('Calculated Effective Dielectric constant = ',e_eff)
					
else:
	print('You typed something other than the given options. Restart the program.')


