import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

s=input("J_perp (format \pm x.xx)=")
FILE=np.loadtxt("Magnetization_Jperp"+s+".txt")

N=50
T=FILE[:,0]
M22=FILE[:,1]
M4=FILE[:,2]
lens=len(T)//N
g=0.5*(3-M4/M22)
G_func=np.zeros(N)

labels=["L=10","L=15","L=12"]
for i in range(lens):
	plt.plot(T[N*i:(N)*(i+1)],g[N*i:(N)*(i+1)],"-",label=labels[i])
plt.legend()
plt.ylim(0.6,1)
plt.xlabel("$T$")
plt.ylabel("$ G $")
plt.title("PMC, $J_{\perp}=$"+s)
plt.show()

G10=g[0:50]
G12=g[50:100]
G15=g[100:150]
Temp=T[0:50]
F=G10+G12-2*G15

func=np.zeros(len(Temp))
parameters=np.polyfit(Temp,F,10)
func=np.zeros(len(Temp[20:]))
#print(parameters)

for i in range(len(parameters)):
    func+=parameters[i]*Temp[20:]**(len(parameters)-i-1)

plt.plot(Temp,F,"o")
plt.plot(Temp[20:],func,"--")
plt.grid()
plt.show()

def function(x):
    total_func=0
    for i in range(len(parameters)):
        total_func+=parameters[i]*x**(len(parameters)-i-1)
    return total_func

print(brentq(function,2.2,2.5))


