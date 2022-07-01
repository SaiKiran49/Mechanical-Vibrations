import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from numpy.linalg import inv
import math
############################interference fit#############
########################inputs###########################
m=1;u=0.7;c1=0.2;c2=0.0;c3=0
K1=0.18;alpha1=150;K2=2.5;alpha2=30;K3=1.33
K11=K1*math.cos(math.pi*alpha1/180)**2+K2*math.cos(math.pi*alpha2/180)**2
K12=K1*math.cos(math.pi*alpha1/180)*math.sin(math.pi*alpha1/180)+K2*math.cos(math.pi*alpha2/180)*math.sin(math.pi*alpha2/180)
K22=K1*math.sin(math.pi*alpha1/180)**2+K2*math.sin(math.pi*alpha2/180)**2+K3
######################################################################
M = np.array([[m, 0],[0,  m]])
C=np.array([[c1,c3],[c3,c2]])
K=np.array([[K11, (K12-u*K3)],[K12, K22-0*K3]])
values=linalg.eigvals(K,M)
frequencies=np.sqrt(values)
a1=1;a2=c1+c2;a3=K11+K22+(c1*c2);a4=K11*c2+c1*K22;a5=K11*K22-(K12*(K12-u*K3))
dampedfreq=np.roots([a1,a2,a3,a4,a5])
###################################################################
#####################state space form ##############
#####################################################################
delta=1.0e-3
endtime=40
F=np.zeros((4,1))
A=np.zeros((4,4))
B=np.zeros((4,4))
Y=np.zeros((4,1))
Y[2]=1;
F=np.zeros((4,1))
I=np.identity(2)
A[0:2,0:2]=M
A[2:4,2:4]=I
B[0:2,2:4]=K
B[2:4,0:2]=-I
B[0:2,0:2]=C
A_inv=inv(A)
x=[];y=[];xdot=[];ydot=[]
###################integration State space solution###############
for t in np.arange(0,endtime,delta):
    Ynew=Y+delta*A_inv.dot(F-B.dot(Y))
    Y=Ynew
    x.extend(Y[2])  
    y.extend(Y[3])
    xdot.extend(Y[0])
    ydot.extend(Y[1])
time=[round(t,5) for t in np.arange(0,endtime,delta)]
plt.plot(time,x)
plt.plot(time,y)
plt.xlabel('time (s)')
plt.ylabel('displacement')
plt.title('response curve at u=%.2f' %u)
plt.legend(['in plane','out of plane'],loc='upper right')
plt.show()
############################plotting energy, Power##############
Pfric=u*K3*np.array(y)*np.array(xdot)
efric=u*K3*np.array(y)*np.array(x)
#Etotal=0.5*(K11*np.array(x)**2+(K12+K12-(u*K3))*np.array(x)*np.array(y)+K22*np.array(y)**2)+0.5*(np.array(xdot)**2+np.array(ydot)**2)
plt.plot(time,Pfric)
#plt.plot(time,efric)
plt.xlabel('time (s)')
plt.ylabel('Frictional power')
plt.show()
plt.plot(time,efric)
plt.xlabel('time (s)')
plt.ylabel('Frictional energy')
plt.show()
####################phase plots###########################
#plt.plot(x,xdot)
#plt.xlabel('x')
#plt.ylabel('xdot')
####################################################
#plt.show()
#plt.plot(y,ydot)
#plt.xlabel('y')
#plt.ylabel('ydot')
#plt.show()
