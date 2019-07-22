import numpy as np 
import matplotlib.pyplot as plt
from sys import  argv

size=argv[1]#input("Enter the size of the lattice=")
FILE=np.loadtxt("Energy_magnetization_L"+size+".txt")
print(FILE.shape)

T=FILE[:,0]
E=FILE[:,1]
Cv=FILE[:,2]
M=FILE[:,3]
MA=FILE[:,4]
MB=FILE[:,5]


print(T[np.argmax(Cv)])
plt.figure(figsize=(10,7))
plt.subplot(311)
plt.plot(T,E,"o-")
plt.ylabel("$ Energy $")
#plt.xlabel("$ T $")
plt.subplot(312)
plt.plot(T,M,"o-")
plt.ylabel(" $ Magnetization $")
#plt.xlabel("$ T $")
plt.subplot(313)
plt.plot(T,Cv,"o-",label="$ T_c= $"+str(T[np.argmax(Cv)]))
plt.ylabel(" $ C $")
plt.xlabel("$ T $")
plt.legend()
plt.show()


plt.subplot(311)
plt.plot(T,2*MA,label="Magnetization A")
plt.legend()
plt.subplot(312)
plt.plot(T,2*MB,label="Magnetization B")
plt.legend()
plt.subplot(313)
plt.plot(T,MA+MB,label="Total magnetization sum")
plt.plot(T,M,"--",label="Total magnetization")
plt.legend()
plt.show()
