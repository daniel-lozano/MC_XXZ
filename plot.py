import numpy as np 
import matplotlib.pyplot as plt
from sys import  argv

size=argv[1]#input("Enter the size of the lattice=")
FILE=np.loadtxt("Energy_magnetization_L"+size+".txt")
print(FILE.shape)

T=FILE[:,0]
E=FILE[:,1]
M=FILE[:,2]
Cv=FILE[:,3]
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
plt.plot(T,Cv,"o-")
plt.plot(T,Cv,"o-")
plt.ylabel(" $ C $")
plt.xlabel("$ T $")
plt.show()



