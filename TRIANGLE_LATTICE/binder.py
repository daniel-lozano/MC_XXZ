import numpy as np
import matplotlib.pyplot as plt
s=input("J_perp (format \pm x.xx)=")
FILE=np.loadtxt("Magnetization_Jperp"+s+".txt")

N=50
T=FILE[:,0]
M22=FILE[:,1]
M4=FILE[:,2]
lens=len(T)//N
g=0.5*(1-(1./3)*M4/M22)
labels=["L=10","L=15","L=12"]
for i in range(lens):
	plt.plot(T[N*i:(N)*(i+1)],g[N*i:(N)*(i+1)],"o-",label=labels[i])
plt.legend()
plt.xlabel("$T$")
plt.ylabel("$ G $")
plt.title("$J_{\perp}=$"+s)
plt.show()

for i in range(lens):
    plt.plot(T[N*i:(N)*(i+1)],M4[N*i:(N)*(i+1)],"o-",label=labels[i])
    plt.errorbar(T[N*i:(N)*(i+1)],M4[N*i:(N)*(i+1)],xerr=None,yerr=FILE[:,4])
plt.legend()
plt.xlabel("$T$")
plt.ylabel("$ \\langle M^4 \\rangle $")
plt.title("$J_{\perp}=$"+s)
plt.show()

#plt.legend()
#plt.xlim(2.15,2.35)
#plt.show()

#
#N=50
#trials=FILE.shape[0]//50
#print(FILE.shape,N,trials)
#
#for i in range(trials):
#    T=FILE[i*50:(i+1)*50,0]
#    M22=FILE[i*50:(i+1)*50,1]
#    M4=FILE[i*50:(i+1)*50,2]
#
#    plt.plot(T,M22,"--",label=str(FILE[i*50:(i+1)*50,3][0]))
#    plt.plot(T,M4,"--",label=str(FILE[i*50:(i+1)*50,3][0]))
#plt.legend()
#plt.show()
#
#for i in range(trials):
#    T=FILE[i*50:(i+1)*50,0]
#    M22=FILE[i*50:(i+1)*50,1]
#    M4=FILE[i*50:(i+1)*50,2]
#
#    plt.plot(T,1-0.33*M4/M22,"--",label=str(FILE[i*50:(i+1)*50,3][0]))
#plt.legend()
#plt.xlim(2.15,2.35)
#plt.show()
