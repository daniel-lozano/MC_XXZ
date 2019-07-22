import numpy as np
import matplotlib.pyplot as plt
FILE=np.loadtxt("Magnetization.txt")
N=50
trials=FILE.shape[0]//50
print(FILE.shape,N,trials)

for i in range(trials):
    T=FILE[i*50:(i+1)*50,0]
    M22=FILE[i*50:(i+1)*50,1]
    M4=FILE[i*50:(i+1)*50,2]
    
    plt.plot(T,M22,"--",label=str(FILE[i*50:(i+1)*50,3][0]))
    plt.plot(T,M4,"--",label=str(FILE[i*50:(i+1)*50,3][0]))
plt.legend()
plt.show()

for i in range(trials):
    T=FILE[i*50:(i+1)*50,0]
    M22=FILE[i*50:(i+1)*50,1]
    M4=FILE[i*50:(i+1)*50,2]
    
    plt.plot(T,1-0.33*M4/M22,"--",label=str(FILE[i*50:(i+1)*50,3][0]))
plt.legend()
plt.xlim(2.15,2.35)
plt.show()
