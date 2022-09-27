import numpy as np
import matplotlib.pyplot as plt

f = np.loadtxt('iterations.txt')
N = f[:,0]
i = f[:,1]
s = i/N
ss = i/N**2


plt.figure()
plt.subplot(311)
plt.plot(N,i, 'o')
plt.xlabel('N'); plt.ylabel('i')
plt.text(1,70000, 'Iterations')
plt.subplot(312)
plt.plot(N,s, 'o')
plt.text(1,350, 'Iterations divided by N.')
plt.xlabel('N'); plt.ylabel('i/N')
plt.subplot(313)
plt.plot(N,ss, 'o')
plt.text(75,1.6, 'Iterations divided by N^2.')
plt.xlabel('N'); plt.ylabel('i/N^2')



i_N = 2*N**2
plt.figure()
plt.plot(N, i_N, label='$i(N) = 2N^2$')
plt.plot(N,i, 'o', label='iterations')
plt.xlabel('N')
plt.ylabel('i')
plt.legend()