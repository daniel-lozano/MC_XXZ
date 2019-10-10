import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


set=int(input("number 3n+1="))
keyword1="Magnetization Density^2"
keyword2="Magnetization Density^4"

File10=np.load("run"+str(set)+".npy")
File12=np.load("run"+str(set+1)+".npy")
File15=np.load("run"+str(set+2)+".npy")

print("Jx=-"+str(File10[0]['Jx']))
##--------------------------------------------------------------------
s=input("J_perp (format \pm x.xx)=")
FILE=np.loadtxt("Magnetization_Jperp"+s+".txt")

N=50
Temp=FILE[:,0]
M22=FILE[:,1]
M4=FILE[:,2]
lens=len(Temp)//N
g=0.5*(3-M4/M22)
G_func=np.zeros(N)

labels=["L=10 PMC" ,"L=15 PMC","L=12 PMC"]
for i in range(lens):
#    plt.plot(T[N*i:(N)*(i+1)],M4[N*i:(N)*(i+1)],"o",label=labels[i])
    plt.plot(Temp[N*i:(N)*(i+1)],g[N*i:(N)*(i+1)],"o",label=labels[i])
plt.legend()
#plt.ylim(0.6,1)
plt.xlabel("$T$")
plt.ylabel("$ G $")
plt.title("PMC, $J_{\perp}=$"+s)

##--------------------------------------------------------------------

#set=int(input("number 3n+1="))
#keyword1="Magnetization Density^2"
#keyword2="Magnetization Density^4"
#
#File10=np.load("run"+str(set)+".npy")
#File12=np.load("run"+str(set+1)+".npy")
#File15=np.load("run"+str(set+2)+".npy")
#
#print("Jx="+str(File10[0]['Jx']))

M4_10=np.zeros(len(File10))
M2_10=np.zeros(len(File10))

M4_12=np.zeros(len(File10))
M2_12=np.zeros(len(File10))

M4_15=np.zeros(len(File10))
M2_15=np.zeros(len(File10))
T=np.zeros(len(File10))
for i in range(len(File10)):
    M4_10[i]=File10[i][keyword2]
    M2_10[i]=File10[i][keyword1]
    
    M4_12[i]=File12[i][keyword2]
    M2_12[i]=File12[i][keyword1]
    
    M4_15[i]=File15[i][keyword2]
    M2_15[i]=File15[i][keyword1]
    
    T[i]=(File10[i]["T"])

G10=0.5*(3-M4_10/M2_10**2)
G12=0.5*(3-M4_12/M2_12**2)
G15=0.5*(3-M4_15/M2_15**2)

plt.plot(T,G10,label="10 QMC")
plt.plot(T,G12,label="12 QMC")
plt.plot(T,G15,label="15 QMC")

plt.legend()
plt.ylim(0.6,1)
plt.title("QMC, $J_\perp$="+str(File10[0]['Jx']))
plt.show()

F=G10+G12-2*G15
parameters=np.polyfit(T[10:],F[10:],6)
func=np.zeros(len(T[10:]))
#print(parameters)

for i in range(len(parameters)):
    func+=parameters[i]*T[10:]**(len(parameters)-i-1)

plt.plot(T,F,"o")
plt.plot(T[10:],func,"--")
plt.grid()
plt.show()

def function(x):
    total_func=0
    for i in range(len(parameters)):
        total_func+=parameters[i]*x**(len(parameters)-i-1)
    return total_func

print(brentq(function,2.2,2.5))



for i in range(lens):
    plt.plot(Temp[N*i:(N)*(i+1)],M4[N*i:(N)*(i+1)],"o",label=labels[i]+" M4")
plt.legend()
#plt.ylim(0.6,1)
plt.xlabel("$T$")
plt.ylabel("$ \\langle M^4\\rangle $")
plt.title(" Get rid of the factor  @ $ J_{\perp}=$"+s)

factor=M4[0]/M4_10[-1]
print(factor)
plt.plot(T,M4_10*factor,label="10 QMC M4")
plt.plot(T,M4_12*factor,label="12 QMC M4")
plt.plot(T,M4_15*factor,label="15 QMC M4")
plt.show()


for i in range(lens):
    plt.plot(Temp[N*i:(N)*(i+1)],M22[N*i:(N)*(i+1)],"o",label=labels[i]+" M22")
plt.legend()
#plt.ylim(0.6,1)
plt.xlabel("$T$")
plt.ylabel("$ \\langle M^2\\rangle^2  $")
plt.title(" Get rid of the factor $J_{\perp}=$"+s)

factor=M4[0]/M4_10[-1]
print(factor)
plt.plot(T,factor*M2_10**2,label="10 QMC")
plt.plot(T,factor*M2_12**2,label="12 QMC")
plt.plot(T,factor*M2_15**2,label="15 QMC")
plt.legend()


plt.show()
