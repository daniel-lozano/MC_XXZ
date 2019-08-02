import numpy as np 
import matplotlib.pyplot as plt
from sys import  argv

size=argv[1]#input("Enter the size of the lattice=")
Jperp=input("Value of Jperp=")
FILE=np.loadtxt("Energy_magnetization_L"+size+"_Jperp"+Jperp +".txt")
compare=input("Comparing the methods? yes (1) or no (0)=")

if(compare=="1"):
    FILE1=np.loadtxt("Energy_magnetization_J123_L"+size+"_Jperp"+Jperp +".txt")
    T1=FILE1[:,0]
    E1=FILE1[:,1]
    Cv1=FILE1[:,2]
    M1=FILE1[:,3]
    MA1=FILE1[:,4]
    MB1=FILE1[:,5]

N=50
FILE_O=np.loadtxt("Data_"+str(N)+"_no_Jperpp.txt")


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
plt.title("$ J_{\\perp}=$" + Jperp+"$J$")
plt.plot(T,E,"k.-",label=size+"$\\times $"+size)
plt.plot(T_O,E_O,"r--",label=str(N)+"$\\times $"+str(N)+"$\ J_\perp=0 $")
if(compare=="1"):
    plt.plot(T1,E1,"b.-",label=size+"$\\times $"+size+" J123")

plt.ylabel("$ Energy $")
plt.legend()
plt.ylim(-2.5,0.0)
plt.subplot(312)
plt.plot(T,M,"k^-",label=size+"$\\times $"+size)
plt.plot(T_O,M_O,"r.-",label=str(N)+"$\\times $"+str(N)+"$\ J_\perp=0 $")
if(compare=="1"):
    plt.plot(T1,M1,"b.-",label=size+"$\\times $"+size+" J123")
plt.ylabel(" $ Magnetization $")
plt.legend()

#plt.xlabel("$ T $")
plt.subplot(313)
plt.plot(T,Cv,"k^-",label="$ T_c= $"+str(T[np.argmax(Cv)]))
plt.plot(T_O,Cv_O,"r.-",label="$ T_c= $"+str(T_O[np.argmax(Cv_O)]))
if(compare=="1"):
    plt.plot(T1,Cv1,"b.-",label=size+"$\\times $"+size+" J123")

plt.ylabel(" $ C $")
plt.xlabel("$ T $")
plt.ylim(0,3)
#plt.xlim(1.5,max(T))
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
