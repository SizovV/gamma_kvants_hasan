import numpy as np
import subprocess
import numpy.linalg as linalg
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pandas as pd
import seaborn as sea
import subprocess
sea.set()
import matplotlib.pyplot as plt
<<<<<<< HEAD

def main():
    subprocess.call("FirstWork.exe", shell=True)

    d11 = open("data_elips1.txt").read()
    d11 = d11.split("\n")
    r=0
    dn=[]
    print(len(d11))
    while(r<len(d11)-1):
        i=0
        d1_new= d11[r].split("\t")
        d2 = []
        for i in range(len(d1_new)-1):
            if (float(d1_new[i])!=0.0):
                d2.append(float(d1_new[i]))
        dn.append(d2)
        r+=1
    print("Fuuhhh")
    dn = pd.DataFrame(dn)
    dn = dn.apply(lambda x: pd.Series(x.dropna().values))
    dn = dn.sort_values(by=[0, 1])
    fig = plt.figure()
    print("Fuuuuuuuhhh-n-da")

    ax1 = fig.add_subplot(111)
    #ax2 = fig.add_subplot(2, 2, 2)
    #ax3 = fig.add_subplot(2, 2, 3)
    #ax4 = fig.add_subplot(2, 2, 4)
    ax1.set(ylabel="\\nu", xlabel="E, MeV")
    #ax2.set(title="h=7.5", ylabel="\\nu")
    #ax3.set(title="h=15", ylabel="\\nu", xlabel="E, MeV")
    #ax4.set(title="h=30", ylabel="\\nu", xlabel="E, MeV")
    num = 10000
    ax1.plot(dn[0].tolist(), dn[1].rolling(num, min_periods=1).mean().tolist(), "-", color='g',  label = 'h=0 cm')
    ax1.plot(dn[0].tolist(), dn[2].rolling(num, min_periods=1).mean().tolist(), "-", color='r',  label = 'h=7.5 cm')
    ax1.plot(dn[0].tolist(), dn[3].rolling(num, min_periods=1).mean().tolist(), "-", color='b',  label = 'h=15 cm')
    ax1.plot(dn[0].tolist(), dn[5].rolling(num, min_periods=1).mean().tolist(), "-", color='y',  label = 'h=30 cm')
    ax1.legend()
    #plt.suptitle('symmetry axes: side of the cylinder')
    plt.suptitle('symmetry axes: center of the cylinder')

    plt.show()

if __name__ == "__main__":
    main()
=======
#matplotlib.use('Qt5Agg')

#def main():
subprocess.call("FirstWork.exe", shell=True)

#    d11 = open("data_elips1.txt").read()
#    d11 = d11.split("\n")
#    r=0
#    dn=[]
#    print(len(d11))
#    while(r<len(d11)-1):
#        i=0
#        d1_new= d11[r].split("\t")
#        d2 = []
#        for i in range(len(d1_new)-1):
#            if (float(d1_new[i])!=0.0):
#                d2.append(float(d1_new[i]))
#        dn.append(d2)
#        r+=1
#    print("Fuuhhh")
#    dn = pd.DataFrame(dn)
#    dn = dn.apply(lambda x: pd.Series(x.dropna().values))
#    dn = dn.sort_values(by=[0, 1])
#    fig = plt.figure()
#    print("Fuuuuuuuhhh-n-da")
#    ax1 = fig.add_subplot(2, 2, 1)
#    ax2 = fig.add_subplot(2, 2, 2)
#    ax3 = fig.add_subplot(2, 2, 3)
#    ax4 = fig.add_subplot(2, 2, 4)
#    ax1.set(title="h=0", ylabel="\\nu")
#    ax2.set(title="h=7.5", ylabel="\\nu")
#    ax3.set(title="h=15", ylabel="\\nu", xlabel="E, MeV")
#    ax4.set(title="h=30", ylabel="\\nu", xlabel="E, MeV")
#    num = 5000
#    ax1.plot(dn[0].tolist(), dn[1].rolling(num, min_periods=1).mean().tolist(), "-", color='g')
#    ax2.plot(dn[0].tolist(), dn[2].rolling(num, min_periods=1).mean().tolist(), "-", color='r')
#    ax3.plot(dn[0].tolist(), dn[3].rolling(num, min_periods=1).mean().tolist(), "-", color='b')
#    ax4.plot(dn[0].tolist(), dn[5].rolling(num, min_periods=1).mean().tolist(), "-", color='y')

    #plt.suptitle('symmetry axes: side of the cylinder')
    #plt.suptitle('symmetry axes: center of the cylinder')

    #plt.show()

#if __name__ == "__main__":
    #main()
>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33





<<<<<<< HEAD
#d1 = open("data_elips.txt").read()
#d1 = d1.split("|")
#r=0
#d=[]
#
#while(r<len(d1)-1):
#    i=0
#    d1_new= d1[r].split()
#    d2 = []
#    while(i<len(d1_new)):
#        d3 = []
#        for k in range(3):
#            d3.append(float(d1_new[i+k]))
#        d2.append(d3)
#        i+=3
#    d.append(d2)
#    r+=1

#fig = plt.figure()
#ax = fig.add_subplot(111, projection="3d")
#ax = Axes3D(fig)


#A = np.array([[1/(10**2),0,0],[0,1/(30**2),0],[0,0,1/1e-20]])
#center = [0,0,0]

# find the rotation matrix and radii of the axes
#U, s, rotation = linalg.svd(A)
#radii = 1.0/np.sqrt(s)

# now carry on with EOL's answer
#u = np.linspace(0.0, 2.0 * np.pi, 100)
#v = np.linspace(0.0, np.pi, 100)
#x = radii[0] * np.outer(np.cos(u), np.sin(v))
#y = radii[1] * np.outer(np.sin(u), np.sin(v))
#z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
#for i in range(len(x)):
#    for j in range(len(x)):
#        [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation) + center

# plot
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color='b', alpha=0.2)
#ax.set_zlim(0, 20)
#ax.set_ylim(-30, 30)
#ax.set_xlim(-30, 30)

#ax.set_xlabel("x")
#ax.set_ylabel("y")
#ax.set_zlabel("z")
#lns_new = []
#for i in range(len(d)):
#    colors = cm.rainbow(np.linspace(0, 1, len(d[i])))
#    d_new = np.transpose(d[i]).tolist()
#    ls1 = ax.plot3D(d_new[0], d_new[1], d_new[2], color='r', linewidth=0.5)
#    for j, c in zip(range(2), colors):
#        ls = ax.scatter(d_new[0][j], d_new[1][j], d_new[2][j], s=15, color=c)
#        lns_new.append(ls)

#plt.show()
=======
d1 = open("data_elips.txt").read()
d1 = d1.split("|")
r=0
d=[]

while(r<len(d1)-1):
    i=0
    d1_new= d1[r].split()
    d2 = []
    while(i<len(d1_new)):
        d3 = []
        for k in range(3):
            d3.append(float(d1_new[i+k]))
        d2.append(d3)
        i+=3
    d.append(d2)
    r+=1

fig = plt.figure()
#ax = fig.add_subplot(111, projection="3d")
ax = Axes3D(fig)


A = np.array([[1/(10**2),0,0],[0,1/(30**2),0],[0,0,1/1e-20]])
center = [0,0,0]

# find the rotation matrix and radii of the axes
U, s, rotation = linalg.svd(A)
radii = 1.0/np.sqrt(s)

# now carry on with EOL's answer
u = np.linspace(0.0, 2.0 * np.pi, 100)
v = np.linspace(0.0, np.pi, 100)
x = radii[0] * np.outer(np.cos(u), np.sin(v))
y = radii[1] * np.outer(np.sin(u), np.sin(v))
z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
for i in range(len(x)):
    for j in range(len(x)):
        [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation) + center

# plot
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color='b', alpha=0.2)
ax.set_zlim(0, 20)
ax.set_ylim(-30, 30)
ax.set_xlim(-30, 30)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
lns_new = []
for i in range(len(d)):
    colors = cm.rainbow(np.linspace(0, 1, len(d[i])))
    d_new = np.transpose(d[i]).tolist()
    ls1 = ax.plot3D(d_new[0], d_new[1], d_new[2], color='r', linewidth=0.5)
    for j, c in zip(range(2), colors):
        ls = ax.scatter(d_new[0][j], d_new[1][j], d_new[2][j], s=15, color=c)
        lns_new.append(ls)

plt.show()

>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33
