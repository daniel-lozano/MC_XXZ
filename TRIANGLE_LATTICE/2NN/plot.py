import numpy as np 
import matplotlib.pyplot as plt
from sys import  argv

size=argv[1]#input("Enter the size of the lattice=")

FILE=np.loadtxt("Energy_magnetization_L"+size+".txt")




print(FILE.shape)

T=FILE[:,0]
E=FILE[:,1]
E_error=FILE[:,9]
M=FILE[:,2]
Cv=FILE[:,3]
Dd=FILE[:,4]

plt.figure(figsize=(10,7))
plt.subplot(411)
plt.plot(T,E,"o-")
plt.errorbar(T,E,xerr=None,yerr=E_error)
plt.ylabel("$ E/J $")
#plt.xlabel("$ T $")
plt.subplot(412)
plt.plot(T,M,"o-")
plt.ylabel(" $ M_{q=0} $")
#plt.xlabel("$ T $")
plt.subplot(413)
plt.plot(T,Cv,"o-")
#plt.plot(T[:-1],(E[1:]-E[:-1])/(T[1]-T[0]),".-")
plt.ylabel(" $ C $")
plt.xlabel("$ T $")
plt.subplot(414)
plt.plot(T,Dd,"o-")
#plt.plot(T[:-1],(E[1:]-E[:-1])/(T[1]-T[0]),".-")
plt.ylabel(" $ \\rho_{\\mathrm{droplet}} $")
plt.xlabel("$ T $")
plt.show()
print(T[:-1][np.argmax((E[1:]-E[:-1])/(T[1]-T[0]))])


if(FILE.shape[1]>5):
    M2A=FILE[:,5]
    M4A=FILE[:,6]
    X1=FILE[:,7]
    X3=FILE[:,8]



    plt.plot(T,0.5*(3-M4A/M2A**2),"o-",label=str(size)+" Binder ratio")

    if(len(argv)>2 ):
        size2=argv[2]
        FILE2=np.loadtxt("Energy_magnetization_L"+size2+".txt")
        M2B=FILE2[:,5]
        M4B=FILE2[:,6]
        plt.plot(T,0.5*(3-M4B/M2B**2),"o-",label=str(size2))
        
    if(len(argv)>3 and len(argv)<5):
        size3=argv[3]
        FILE3=np.loadtxt("Energy_magnetization_L"+size3+".txt")
        M2C=FILE3[:,5]
        M4C=FILE3[:,6]
        plt.plot(T,0.5*(3-M4C/M2C**2),"o-",label=str(size3))

    plt.legend()
    plt.show()
    
    plt.plot(T,X1,label="X1")
    plt.plot(T,X3,label="X3")
    plt.legend()
    plt.show()

