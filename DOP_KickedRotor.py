from DOP_Dalibard import *
from DOP_Floquet import *
from DOP_Brillouin import *
import matplotlib.pyplot as plt
import numpy as np
import cmath

im = cmath.sqrt(-1)	
J= 10.0	
p1=1/J
p2=1/J
K1=10.0
K2=15.0
T=1.0

for epsilon in [0.0]:#,2.0,5.0,10.0,50.0,100.0,500.0,1000.0,10000.0]:
	N=int(2*J+1)
	N_sq=int(N*N)
	

	Energy,Eigenvectors1=Dalibard(J,p1,p2,K1,K2,epsilon)
	QuasiEnergy,Eigenvectors2=Floquet(J,p1,p2,K1,K2,epsilon)
	Energy_B, Eigenvectors_B = Brillouin(J,p1,p2,K1,K2,epsilon,T)
	QuasiEnergy=-im*np.log(QuasiEnergy) + np.pi
	Energy=Energy%(2*np.pi)
	
	Sv_array1=[]
	Sv_array2=[]		#To store the entanglement values for different eigenstates
	Sv_array3=[]

	#To calculate the reduced density matrix of the system (C_Dagger * C)
	for vector in Eigenvectors1:
		C=np.zeros((N,N),dtype=np.complex64)
		dummy_index=0
		for i in range(N_sq):
			C[dummy_index][i%N]=vector[i]
			if i%N == N-1:
				dummy_index=dummy_index+1
		C=np.matrix(C)
		DensityMatrix1 = C*Dagger(C)
		Lambda1=eigvalsh(DensityMatrix1)
		Sv=0.0
		for x in Lambda1:
			if x > 0:
				Sv=Sv-x*np.log(x)
		Sv_array1.append(Sv)
	
	for vector in Eigenvectors2:
		C=np.zeros((N,N),dtype=np.complex64)
		dummy_index=0
		for i in range(N_sq):
			C[dummy_index][i%N]=vector[i]
			if i%N == N-1:
				dummy_index=dummy_index+1
		C=np.matrix(C)
		DensityMatrix2 = C*Dagger(C)
		Lambda2=eigvalsh(DensityMatrix2)
		Sv=0.0
		for x in Lambda2:
			if x > 0:
				Sv=Sv-x*np.log(x)
		Sv_array2.append(Sv)

	for vector in Eigenvectors_B:
		C=np.zeros((N,N),dtype=np.complex64)
		dummy_index=0
		for i in range(N_sq):
			C[dummy_index][i%N]=vector[i]
			if i%N == N-1:
				dummy_index=dummy_index+1
		C=np.matrix(C)
		DensityMatrix3 = C*Dagger(C)
		Lambda1=eigvalsh(DensityMatrix3)
		Sv=0.0
		for x in Lambda1:
			if x > 0:
				Sv=Sv-x*np.log(x)
		Sv_array3.append(Sv)
	
	plt.axis([0, N_sq, -0.5, 4])
	plt.title('J=%s p1=%s p2=%s K1=%s K2=%s epsilon=%s' %(J,p1,p2,K1,K2,epsilon))
	plt.xlabel('Eigenvectors')
	plt.ylabel('Entanglement')
	plt.plot(Sv_array1,'bo',label='From Dalibard')
	plt.plot(Sv_array2,'ro',label='From Floquet')
	plt.plot(Sv_array3,'go',label='From Brillouin')
	plt.legend()
	plt.savefig('2_epsilon_%s.png'%epsilon)
	plt.close()
	
	print 'epsilon %s done' %epsilon


	#plt.show()

	#QE = np.real(QuasiEnergy)
	#plt.figure()
	#plt.axis([0, 7, 0, 75])
	#
	#plt.hist(QE)
	#plt.savefig('QuasiEnergy.png')
	#plt.figure()
	#plt.axis([0, 7, 0, 75])
#
	#plt.hist(Energy)
	#plt.savefig('Energy.png')
#
	#plt.show()#