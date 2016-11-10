import numpy as np
from numpy.linalg import *
import scipy.linalg as sp
import cmath
from sympy.physics.quantum import Dagger, TensorProduct
from DOP_Heff import *
import matplotlib.pyplot as plt

def Floquet(J,p1,p2,K1,K2,epsilon):

	im = cmath.sqrt(-1)	

	#J=10.0
	#p1=1.0
	#p2=1.0
	#K1=1.0
	#K2=1.0
	#epsilon = 1.0

	N = int(2*J+1)
	N_sq = N*N
	m_arr = np.linspace(J,-J,num=N)	
	Jy= [[0 for x in range(N)] for x in range(N)]	#To calculate the representation of Jy
	Identity = np.identity(N)
 
#To calculate the representation of Jy using ladder operators
	for index in range(N-1):
		Jy[index+1][index] = -np.sqrt((J+m_arr[index])*(J-m_arr[index]+1))
    	Jy[index][index+1] = np.sqrt((J-m_arr[index+1])*(J+m_arr[index+1]+1))
	Jy = np.matrix(Jy)
	Jy=(-0.5*im)*Jy


#To calculate the representation of Jz
	Jz=[[0 for x in range(N)] for y in range(N)]
	for x in range(N):
		m=-x+J
		Jz[x][x]=m
	Jz=np.matrix(Jz)
	Jz_sq=Jz*Jz
	J1=(Jz+0.5*Identity)
	J_sq=J1*J1

	JyEigenvalues, JyBasis_inv = eigh(Jy)
	JyBasis = Dagger(JyBasis_inv)
	Jy_prime=JyBasis_inv*Jy*JyBasis

	Uf_prime=np.zeros((N,N),dtype=np.complex64)
	for x in range(N):
		Uf_prime[x,x]=np.exp(-im*0.5*np.pi*Jy_prime[x,x])

	Uf_orig=JyBasis*Uf_prime*JyBasis_inv

	Uk_1=np.zeros((N,N),dtype=np.complex64)
	for x in range(N):
		Uk_1[x,x]=np.exp(-im*K1*Jz_sq[x,x]/(2*J))

	Uk_2=np.zeros((N,N),dtype=np.complex64)
	for x in range(N):
		Uk_2[x,x]=np.exp(-im*K2*J_sq[x,x]/(2*J))

	U_12=np.zeros((N_sq,N_sq),dtype=np.complex64)
	JzJ1 = TensorProduct(Jz,J1)
	for x in range(N_sq):
		U_12[x,x]=np.exp(-im*epsilon*JzJ1[x,x]/J)	#J_sq could be Jz*(Jz+0.5*Identity)

	U1=Uk_1*Uf_orig
	U2=Uk_2*Uf_orig
	UT=U_12*(TensorProduct(U1,U2))

	eigval, eigvec = eig(UT)

	eigvec=np.array(eigvec)
	eigvec=np.transpose(eigvec)
	return eigval, eigvec
	
#	Sv_array=[]		#To store the entanglement values for different eigenstates
##To calculate the reduced density matrix of the system (C_Dagger * C)
#	for vector in eigvec:
#		C=np.zeros((N,N),dtype=np.complex64)
#		dummy_index=0
#
#		for i in range(N_sq):
#			C[dummy_index][i%N]=vector[i]
#			if i%N == N-1:
#				dummy_index=dummy_index+1
#
#		C=np.matrix(C)
#		DensityMatrix1 = C*Dagger(C)
#		Lambda=eigvalsh(DensityMatrix1)
#		Sv=0.0
#		for x in Lambda:
#			if x > 0:
#				Sv=Sv-x*np.log(x)
#	
#		Sv_array.append(Sv)
#
#	plt.axis([0, N_sq, -0.5, 4])
#	plt.xlabel('Eigenvalues')
#	plt.ylabel('Entanglement')
#	plt.plot(Sv_array,'bo')
#	plt.show()#