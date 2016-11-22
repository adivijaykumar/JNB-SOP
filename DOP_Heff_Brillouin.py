###################################################################################################################################################
########This function calculates the Effective Hamiltonian of the Coupled Kicked Rotor (Dirac Detla) Potential function using Brillouin Wigner
###################################################################################################################################################

import numpy as np
from numpy.linalg import *
import cmath
from sympy.physics.quantum import Dagger
from DOP_Potential import *
from DOP_H0 import *

im = cmath.sqrt(-1)

def Heff_Brillouin(J,p1,p2,K1,K2,epsilon,T):

	Omega=2*np.pi
	H0=KineticEnergy(p1,p2,J)
	V=Potential(K1,K2,J,epsilon)
	Hamilt = H0 + V/T + 1/12.0*(V*H0*V-V*V*V/T)
	return Hamilt