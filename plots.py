'''
Uses the csv file created by nbody.c to perform analysis
Plots the position data for x, y, and z as a function of time
Renders 2D animations in the xy, xz, and yz planes as well as a 3D animation
'''

import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import animation

writer_vid = animation.FFMpegWriter(fps=60) #FFMpeg writer to save animations
data_csv = 'Data\\pos_data.csv' #Data file from nbody.c

#Constant values taken from the corresponding run of nbody.c
tf = 10000
dt = 0.1
dim_n = 3
part_n = 50
part_num = range(part_n)

#Extracting number of time steps from csv file
with open(data_csv) as f:
    reader = csv.reader(f, delimiter=',')
    time_n = int((sum(1 for row in reader) - 1)/part_n)

#3D array containing position data (x,y,z) for each particle at each time step
pos_data = np.empty([time_n, part_n, dim_n])
time = np.linspace(0, tf, time_n)

#Filling pos_data array with data from the csv
print('Reading csv')
with open(data_csv) as f:
    reader = csv.reader(f, delimiter=',') 
    line = 0
    t = 0 #time tracker
    p = 0 #particle tracker
    for i in reader:
        if line == 0:
            line += 1
            continue
  
        pos_data[t][p][0] = i[2]
        pos_data[t][p][1] = i[3]
        pos_data[t][p][2] = i[4]
        p += 1
        
        if p == part_n:
            p = 0
            t += 1

#List of particle numbers for the legend of the plots
label_list = []
for i in range(part_n):
    label_list.append(str(i))

#Plotting each particle's x position as a function of time
print('Plotting x postion')
for i in range(part_n):
    plt.plot(time, pos_data[:,i,0])
plt.ylim((-100,100))
plt.xlabel('Time')
plt.ylabel('px')
plt.savefig('Figures\\px.png')
plt.clf()

#Plotting each particle's y position as a function of time
print('Plotting y postion')
for i in range(part_n):
    plt.plot(time, pos_data[:,i,1])
plt.ylim((-100,100))
plt.xlabel('Time')
plt.ylabel('py')
plt.savefig('Figures\\py.png')
plt.clf()

#Plotting each particle's z position as a function of time
print('Plotting z postion')
for i in range(part_n):
    plt.plot(time, pos_data[:,i,2])
plt.ylim((-100,100))
plt.xlabel('Time')
plt.ylabel('pz')
plt.savefig('Figures\\pz.png')
plt.clf()

print('Plots saved')

#Lists to hold the points rendered in the 2D animations
pts_xy = []
pts_xz = []
pts_yz = []

#Setting up 2D animation figure
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.set_size_inches(18, 5)

#Adjustable axis limits
ax1.set_xlim((-100,100))
ax1.set_ylim((-100,100))
ax2.set_xlim((-100,100))
ax2.set_ylim((-100,100))
ax3.set_xlim((-100,100))
ax3.set_ylim((-100,100))

ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax2.set_xlabel('x')
ax2.set_ylabel('z')
ax3.set_xlabel('y')
ax3.set_ylabel('z')

#Appends the appropriate amount of point objects to our lists
for i in range(part_n):
    pts_xy.append(ax1.plot([], [], 'o', ms=5)[0])
    pts_xz.append(ax2.plot([], [], 'o', ms=5)[0])
    pts_yz.append(ax3.plot([], [], 'o', ms=5)[0])

ax1.legend(label_list, loc='upper right', bbox_to_anchor=(-0.2, 1))

#Draw function for animation.FuncAnimation, fills data for each point at given time step
def drawframe(t):
    for i, pt in enumerate(pts_xy):
        pt.set_data(pos_data[t,i,0], pos_data[t,i,1])
        
    for i, pt in enumerate(pts_xz):
        pt.set_data(pos_data[t,i,0], pos_data[t,i,2])
        
    for i, pt in enumerate(pts_yz):
        pt.set_data(pos_data[t,i,1], pos_data[t,i,2])
    
    return

#Rendering animation (using 1 out of every 200 time steps)
print('Creating 2D Animation')
anim_2d = animation.FuncAnimation(fig, drawframe, frames=range(0, len(pos_data), 200), interval=10, blit=False)

#Saving animation to mp4 file
print('Saving 2D Animation')
anim_2d.save('Figures\\2d.mp4', writer=writer_vid)
plt.close()

print('2D animation saved')

#Setting up 3D animation figure
fig_3d = plt.figure(figsize=(10,10))
ax = plt.axes(projection='3d')

#Points will be plotted in the figure as a scatter plot
scat = ax.scatter([], [], [], s=50)

#Adjustable axis limits
ax.set_xlim((-100,100))
ax.set_ylim((-100,100))
ax.set_zlim((-100,100))

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

#Draw function for animation.FuncAnimation, fills the scatter plot with x, y, and z positions at the given time step
def drawframe_3d(t):
    scat._offsets3d = (pos_data[t,:,0], pos_data[t,:,1], pos_data[t,:,2])
    return

#Rendering animation (using 1 out of every 200 time steps)
print('Creating 3D Animation')
anim_3d = animation.FuncAnimation(fig_3d, drawframe_3d, frames=range(0, len(pos_data), 200), interval=10, blit=False)

#Saving animation to mp4 file
print('Saving 3D Animation')
anim_3d.save('Figures\\3d.mp4', writer=writer_vid)
plt.close()

print('3D animation saved')