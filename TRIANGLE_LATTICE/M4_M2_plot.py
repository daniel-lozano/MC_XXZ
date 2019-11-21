import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

s=input("J_perp (format \pm x.xx)=")
FILE=np.loadtxt("Magnetization_Jperp"+s+".txt")

N=50
T=FILE[:,0]
M22=FILE[:,1]
M4=FILE[:,2]
Error_M2=FILE[:,3]
Error_M4=FILE[:,4]
lens=len(T)//N
g=0.5*(3-M4/M22)
G_func=np.zeros(N)

labels=["L=10","L=12","L=15"]
plt.subplot(121)
for i in range(lens):

    plt.plot(T[N*i:(N)*(i+1)],M4[N*i:(N)*(i+1)])
    plt.errorbar(T[N*i:(N)*(i+1)],M4[N*i:(N)*(i+1)],xerr=None,yerr=Error_M4[N*i:(N)*(i+1)],label=labels[i])

plt.legend()
plt.ylim(0.0,max(M4)*1.05)
plt.xlabel("$T$")
plt.ylabel("$ \\langle M^4 \\rangle $")
plt.title("PMC, $J_{\perp}=$"+s)
#plt.show()
plt.subplot(122)
for i in range(lens):
    
    plt.plot(T[N*i:(N)*(i+1)],np.sqrt(M22[N*i:(N)*(i+1)]))
    
    plt.errorbar(T[N*i:(N)*(i+1)],np.sqrt(M22[N*i:(N)*(i+1)]),xerr=None,yerr=Error_M2[N*i:(N)*(i+1)],label=labels[i])
    
plt.legend()
plt.ylim(0.0,max(np.sqrt(M22))*1.05)
plt.xlabel("$T$")
plt.ylabel("$ \\langle M^2 \\rangle $")
plt.title("PMC, $J_{\perp}=$"+s)

plt.show()

