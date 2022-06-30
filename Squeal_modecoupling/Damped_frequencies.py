"""
Created on Sun May 10 22:15:38 2020

@author: Sai Kiran
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from numpy.linalg import inv
import math
############################interference fit#############
########################inputs###########################
m=1;delta=0.02;c1=0.3;c2=0.1;c3=0
K1=0.18;alpha1=150;K2=2.5;alpha2=30;K3=1.33
K11=K1*math.cos(math.pi*alpha1/180)**2+K2*math.cos(math.pi*alpha2/180)**2
K12=K1*math.cos(math.pi*alpha1/180)*math.sin(math.pi*alpha1/180)+K2*math.cos(math.pi*alpha2/180)*math.sin(math.pi*alpha2/180)
K22=K1*math.sin(math.pi*alpha1/180)**2+K2*math.sin(math.pi*alpha2/180)**2+K3
######################################################################
M = np.array([[m, 0],[0,  m]])
C=np.array([[c1,c3],[c3,c2]])
f1=[];f2=[];r1=[];r2=[];df1=[];df2=[];dr1=[];dr2=[]
for u in np.arange(0,1,delta):
    v=np.zeros((2,1),dtype=complex);
    vd=np.zeros((2,1),dtype=complex);
    K=np.array([[K11, (K12-u*K3)],[K12, K22-0*K3]])
    values=linalg.eigvals(K,M)
    for i in range(2):
        v[i]=values[i]
    frequencies=np.sqrt(v)
    actfreq=complex(0,1)*frequencies
    a1=1;a2=c1+c2;a3=K11+K22+(c1*c2);a4=K11*c2+c1*K22;a5=K11*K22-(K12*(K12-u*K3))
    dampedfreq=np.roots([a1,a2,a3,a4,a5])
    vd[0]=dampedfreq[0];vd[1]=dampedfreq[2]
    f1.extend(actfreq.imag[0])
    f2.extend(actfreq.imag[1])
    r1.extend(actfreq.real[0])
    r2.extend(actfreq.real[1])
    df1.extend(vd.imag[0])
    df2.extend(vd.imag[1])
    dr1.extend(vd.real[0])
    dr2.extend(vd.real[1])
#######################Frequency vs coeff plot################################
coef=[round(u,2) for u in np.arange(0,1,delta)]
plt.plot(coef,f1)
plt.plot(coef,f2)
plt.scatter(coef,df1,color='k')
plt.scatter(coef,df2,color='k')
plt.xlabel('coefficient')
plt.ylabel('frequency (rad/s)')
plt.show()
#######################Realpart vs coeff plot#####
plt.plot(coef,r1)
plt.plot(coef,r2)
plt.scatter(coef,dr1,color='k')
plt.scatter(coef,dr2,color='k')
plt.xlabel('coefficient')
plt.ylabel('Real part')
plt.show()
######################frequency vs Real part plot###
plt.scatter(r1,f1,color='g')
plt.scatter(r2,f2,color='g')
plt.scatter(dr1,df1,color='k')
plt.scatter(dr2,df2,color='k')
plt.xlabel('Real part')
plt.ylabel('Frequency (rad/s')
plt.show()
