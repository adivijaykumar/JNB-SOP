#######################################################################################################################################
########This function calculates the entanglement of the system specified by the variables declared below #############################
#######################################################################################################################################

import numpy as np
from numpy.linalg import *
import cmath
from sympy.physics.quantum import Dagger
from DOP_Heff_Brillouin import *
import matplotlib.pyplot as plt

def Brillouin(J,p1,p2,K1,K2,epsilon,T):
	im = cmath.sqrt(-1)

#	J=10.0			#Eigenvalue of J^2 operator, enter value as a float
#	p1=1.0			#Moment of Inertia of first top
#	p2=1.0			#Moment of Inertia of first top
#	K1=1.0			#Spring constant of Top 1
#	K2=1.0			#Spring constant of Top 2
#	epsilon=1.0		#Coupling strength of the two tops

	Hamilt = Heff_Brillouin(J,p1,p2,K1,K2,epsilon,T)
	eigenvalues,eigenvectors=eigh(Hamilt)	#Calculates eigenvectors and eigenvalues of Hamiltonian using routines in numpy.linalg
	eigenvectors=np.array(eigenvectors)
	eigenvectors=np.transpose(np.array(eigenvectors))
	N=int(2*J+1)
	N_sq=N*N
	return eigenvalues, eigenvectors