import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.colors import LightSource
import numpy as np

file = 'curvature.dat'

# Open the input file in read mode
with open(file, 'r') as input_file:
    # Read all lines from the file
    lines = input_file.readlines()

kmax = 0.005
kres = 10

nkp=2*kres
mid = kres+1
node = 1


offset = []
offset.append((-0.017659606952654991,0.046513917396043679,0.43965460613976798))#+ve
offset.append((0.017665681958398235,0.046638430945586576,0.47514974714462382)) #-ve

k_x = []
k_y = []
k_z = []
c_x = []
c_y = []
c_z = []

magnitude= []
#tolerance = (1/(kmax*0.2))**2


#       2 000 000 at 0.005 
# 200 000 000 000 at 0.00005

for line in lines:
    split_line=line.split()
    for val in split_line:
        val = float(val)

    length = (float(split_line[3])**2+float(split_line[4])**2+float(split_line[5])**2)**0.5
    magnitude.append(length)
    if length > 1*10**3:
        k_x.append(float(split_line[0]))
        k_y.append(float(split_line[1]))
        k_z.append(float(split_line[2]))

        c_x.append(float(split_line[3]))
        c_y.append(float(split_line[4]))
        c_z.append(float(split_line[5]))
    else:
        continue



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_proj_type('ortho')

max_mag = max(magnitude)
min_mag = min(magnitude)

print(max_mag)

for i in range(len(c_x)):
    c_x[i] /= max_mag
    c_y[i] /= max_mag
    c_z[i] /= max_mag


ax.quiver(k_x, k_y, k_z, c_x, c_y, c_z, length=0.15*kmax)
ax.plot(offset[1][0],offset[1][1],offset[1][2],color ='red', marker='x')


plt.show()