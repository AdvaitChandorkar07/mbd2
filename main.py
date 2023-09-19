import math
import matplotlib.pyplot as plt
import numpy as np



def graph_plotter(x,z,color,labe):
    '''Function used to plot the graph to
    display the solution of euler method '''
    plt.title('MBD assignment question 2')
    #naming the x-axis and y-axis
    plt.xlabel('time')
    plt.ylabel('x(t)')
    #plt.plot(x, z)
    plt.plot(x, z, color=color, linewidth=1,  marker='o', markerfacecolor='blue', markersize=2,label=labe)
    #plt.legend(labe)



t=float(input("input tau"))
tau=np.array([0,0,t])
r=float(input("import radius"))
l=float(input("import slant height"))
m=float(input("input mass m"))

h=math.sqrt(l**2-r**2)
sinb=math.sin(r/l)
cosb=math.sqrt(1-sinb**2)


delta=float(input("time step is ="))
ti=float(input("Initial time ="))
tf=float(input("Final time ="))

i=ti
time=[]
while i<tf+delta:
    time.append((i))
    #time=[i for i in range(ti,tf+delta,delta)]
    i+=delta

#initializing the lists
theta =[0 for i in range(len(time))]
omega = [np.array([0,0,0]) for i in range(len(time))]
alpha = [np.array([0,0,0]) for i in range(len(time))]
omega[0]=np.array([0,0,0])

I_old=np.array([[3/10*m*r*r,0,0],[0,3/20*m*r*r + 3/80*m*h*h,0],[0,0,3/20*m*r*r + 3/80*m*h*h]])
A=np.array([[1,0,0],[0,cosb,sinb],[0,-sinb,cosb]])
I_og=np.dot(np.dot(np.linalg.inv(A),I_old),A)
I_inst=I_og

omegaz=[]
omegan=[]
omegacheck=[]
omegana=[]

for i in range(1,len(time)):
    n_cap=np.array([cosb*math.cos(theta[i-1]),cosb*math.sin(theta[i-1]),sinb])
    rcmc=3*h/4*n_cap
    f_term= tau-(np.cross(omega[i-1],np.dot(I_inst,omega[i-1])))-(-np.cross(rcmc,np.cross(omega[i-1],(np.cross(omega[i-1],rcmc)))))
    b_term= I_inst + np.dot(rcmc,rcmc)*np.identity(3)
    omega[i]=omega[i-1]+delta*np.dot(np.linalg.inv(b_term),f_term)
    alpha[i]=(omega[i]-omega[i-1])/delta
    theta[i]=theta[i-1]+delta*omega[i][2]
    R = np.array(
        [[math.cos(theta[i]), -math.sin(theta[i]), 0], [math.sin(theta[i]), math.cos(theta[i]), 0],
         [0, 0, 1]])
    I_inst = np.array(np.dot(np.dot(np.linalg.inv(R), I_og), R))
alpha[0]=(omega[1]-omega[0])/delta
alphamod=[]
for i in range(len(time)):

    omegaz+=[omega[i][2]]
    alphamod+=[np.linalg.norm(alpha[i])]
    omegan+=[-r/l*omegaz[i]]

    #omegana+=[omega[i]-omega[i][2]*n_cap]
    #omegacheck+=[-r/l*omegaz[i]]
    #print(time[i], "  ",omegana[i],np.linalg.norm(omegana[i]),omegacheck[i])

graph_plotter(time,omegaz,'red','1')
graph_plotter(time,omegan,'blue','2')
graph_plotter(time,alphamod,'green','3')
#plt.legend('1','2','3')
plt.legend(bbox_to_anchor =(0.75, 1.15), ncol = 3)
plt.show()