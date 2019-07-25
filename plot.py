import numpy as np 
import matplotlib.pyplot as plt
from sys import  argv

size=argv[1]#input("Enter the size of the lattice=")
FILE=np.loadtxt("Energy_magnetization_L"+size+".txt")
N=50
FILE_O=np.loadtxt("Energy_magnetization_L"+str(N)+"_no_Jperp.txt")

print(FILE.shape)

T=FILE[:,0]
E=FILE[:,1]
Cv=FILE[:,2]
M=FILE[:,3]
MA=FILE[:,4]
MB=FILE[:,5]


T_O=FILE_O[:,0]
E_O=FILE_O[:,1]
Cv_O=FILE_O[:,2]
M_O=FILE_O[:,3]


print(T[np.argmax(Cv)])
plt.figure(figsize=(10,7))
plt.subplot(311)
plt.plot(T,E,"^-",label=size+"$\\times $"+size)
plt.plot(T_O,E_O,"o-",label=str(N)+"$\\times $"+str(N)+"$\ J_\perp=0 $")
plt.ylabel("$ Energy $")
plt.legend()
plt.ylim(-2.5,0.0)
plt.subplot(312)
plt.plot(T,M,"^-",label=size+"$\\times $"+size)
plt.plot(T_O,M_O,"o-",label=str(N)+"$\\times $"+str(N)+"$\ J_\perp=0 $")
plt.ylabel(" $ Magnetization $")
plt.legend()

#plt.xlabel("$ T $")
plt.subplot(313)
plt.plot(T,Cv,"^-",label="$ T_c= $"+str(T[np.argmax(Cv)]))
plt.plot(T_O,Cv_O,"o-",label="$ T_c= $"+str(T_O[np.argmax(Cv_O)]))
plt.ylabel(" $ C $")
plt.xlabel("$ T $")
plt.ylim(0,2.3)
plt.xlim(1.5,3)
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
