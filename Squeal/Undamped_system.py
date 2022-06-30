# -*- coding: utf-8 -*-
"""
Created on Thu May  7 14:12:57 2020

@author: SaiKiran
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from numpy.linalg import inv
import math
########################inputs###########################
m=1;delta=0.01
K1=0.18;alpha1=150;K2=2.5;alpha2=30;K3=1.33
K11=K1*math.cos(math.pi*alpha1/180)**2+K2*math.cos(math.pi*alpha2/180)**2
K12=K1*math.cos(math.pi*alpha1/180)*math.sin(math.pi*alpha1/180)+K2*math.cos(math.pi*alpha2/180)*math.sin(math.pi*alpha2/180)
K22=K1*math.sin(math.pi*alpha1/180)**2+K2*math.sin(math.pi*alpha2/180)**2+K3
######################################################################
####
M = np.array([[m, 0],[0,  m]])
f1=[];f2=[];r1=[];r2=[]
for u in np.arange(0,1,delta): 
    v=np.zeros((2,1),dtype=complex)
    K=np.array([[K11, (K12-u*K3)],[K12, K22-0*K3]])
    values=linalg.eigvals(K,M)
    for i in range(2):
        v[i]=values[i]
    frequencies=np.sqrt(v)
    actualfre=complex(0,1)*frequencies
    realpart=actualfre.real
    frequency=actualfre.imag
    r1.extend(realpart[0])
    r2.extend(realpart[1])
    f1.extend(frequency[0])
    f2.extend(frequency[1])
coef=[round(u,2) for u in np.arange(0,1,delta)]
plt.plot(coef,f1)
plt.plot(coef,f2)
plt.xlabel('Friction coefficient')
plt.ylabel('ImaginaryPart')
plt.show()
plt.plot(coef,r1)
plt.plot(coef,r2)
plt.xlabel('Friction coefficient')
plt.ylabel('Real Part')
plt.show()
plt.figure(figsize=(10,10))
plt.scatter(r1,f1)
plt.scatter(r2,f2)
plt.xlabel('Real part-fre')
plt.ylabel('frequency(rad/s)')
plt.show()
