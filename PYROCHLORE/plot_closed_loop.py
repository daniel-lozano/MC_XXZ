import numpy as np 
import matplotlib.pyplot as plt


sizex=input("Enter size of x=")#input("Enter the size of the lattice=")
sizey=input("Enter size of y=")
sizez=input("Enter size of z=")
print("Energy_"+sizex+"_"+sizey+"_"+sizez+"_H.txt")
FILE=np.loadtxt("Energy_"+sizex+"_"+sizey+"_"+sizez+"_H.txt")
T=FILE[:,0]
E_file=FILE[:,1]
C_v=FILE[:,2]
#S2S3=FILE[:,3]


Cv=np.array((E_file[1:]-E_file[0:-1])/(T[1]-T[0]))
T_red=np.array(T[1:])

N=50
print(len(E_file)//N)

label=["-0.001","-0.01","-0.02","-0.05","-0.1"]
col=["k","r","g","b","violet"]

#for i in range(len(E_file)//N):
#
#    T=FILE[i*N:(i+1)*N,0]
#
#
#    S2S3=FILE[i*N:(i+1)*N,3]
#    plt.plot(T,S2S3,label="s2s3")
#    plt.legend()
#    plt.xlabel("$ T $")
#    plt.ylabel("$ \\langle s_2s_2 -s_2s_3 \\rangle  $ ")
#
#plt.show()

plt.figure(figsize=(12,5))
for i in range(len(E_file)//N):
    
    T=FILE[i*N:(i+1)*N,0]
    E=FILE[i*N:(i+1)*N,1]
    C_v=FILE[i*N:(i+1)*N,2]##/(4*int(sizex)**3)
    
    plt.subplot(121)
    plt.semilogx(T,E,"^-")
    plt.xlabel("$ T/J $",size=15)
    plt.ylabel("Energy",size=15)
    plt.title("Energy per spin, "+sizex+"x"+sizey+"x"+sizez)
    #plt.yticks(np.linspace(-1,0,9))
    plt.grid()
    plt.subplot(122)
    plt.semilogx(T,C_v,"^-")
    plt.xlabel("$ T/J $",size=15)
    plt.ylabel("$ C_v $",size=15)
    plt.title("MC Specific heat, "+sizex+"x"+sizey+"x"+sizez)
    plt.legend()
    plt.grid()

plt.show()

KB=1 #8.6E-2

plt.figure(figsize=(12,5))
for i in range(len(E_file)//N):
    
    T=FILE[i*N:(i+1)*N,0]
    E=FILE[i*N:(i+1)*N,1]
    C_v=FILE[i*N:(i+1)*N,2]##/(4*int(sizex)**3)

    plt.subplot(121)
    plt.plot(T*KB,E,"^-")
    plt.xlabel("$ T/J $",size=15)
    plt.ylabel("Energy",size=15)
    plt.title("Energy per spin, "+sizex+"x"+sizey+"x"+sizez)
    plt.grid()
    plt.subplot(122)
    plt.plot(T,C_v,"^-")
    plt.xlabel("$ T/J $",size=15)
    plt.ylabel("$ C_v/T $",size=15)
    plt.title("Numerical Specific heat, "+sizex+"x"+sizey+"x"+sizez)
    plt.grid()

plt.show()

plt.figure(figsize=(5,5))
for i in range(len(E_file)//N):
    T=FILE[i*N:(i+1)*N,0]
    E=FILE[i*N:(i+1)*N,1]
    C_v=FILE[i*N:(i+1)*N,2]##/(4*int(sizex)**3)
    plt.grid()
    plt.semilogx(T,C_v,"^-",label="$ J_{3a}= $"+label[i],color=col[i])

#plt.xlim(5E-3,10)
plt.legend()
plt.xlabel("$ T/J $",size=20)
plt.ylabel(" $ C_v(k_B) $ ",size=20)
plt.title("MC Specific heat, "+sizex+"x"+sizey+"x"+sizez)
plt.show()



#    plt.subplot(211)
#    plt.grid()
#    plt.semilogx(T_red*KB,(Cv/T_red),"^-",label="Numerical derivative")
#    plt.ylabel("$ C_v/T $",size=15)
#    plt.title("Numerical Specific heat, "+sizex+"x"+sizey+"x"+sizez)
#    plt.legend()
