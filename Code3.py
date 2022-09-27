import numpy as np 
import matplotlib.pyplot as plt

vfile = np.loadtxt('v.txt')                     # Analytic eigenvectors
eigvecfile = np.loadtxt('eigenvectors.txt')     # Numeric eigenvectors

x            = np.arange(0,11)  # x-hat. 
v_analytic   = vfile[0:3,:]                         # Three first eigenvectors. 
v_numeric    = np.transpose(eigvecfile[:, 0:3])     # Three first eigenvectors. 

v_an  = np.zeros([3,11])
v_num = np.zeros([3,11])
for i in range(0,3):
    v_an[i]   = np.insert(v_analytic[i], [0,9], [0,0]) # Adds boundary
    v_num[i]  = np.insert(v_numeric[i],  [0,9], [0,0]) # conditions. 

plt.figure()
plt.subplot(311)
plt.plot(x, v_an[0], 'o', label='analytic')
plt.plot(x, v_num[0], 'o', label='numeric')
plt.xlabel('$\hat{x}$'); plt.ylabel('$v_1$')
plt.legend()
plt.subplot(312)
plt.plot(x, v_an[1], 'o', label='analytic')
plt.plot(x, -v_num[1], 'o', label='numeric')
plt.xlabel('$\hat{x}$'); plt.ylabel('$v_2$')
plt.legend()
plt.subplot(313)
plt.plot(x, v_an[2], 'o', label='analytic')
plt.plot(x, -v_num[2], 'o', label='numeric')
plt.xlabel('$\hat{x}$'); plt.ylabel('$v_3$')
plt.legend()