import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.colors import LightSource
import numpy as np

file = 'miller.dat'


# Open the input file in read mode
with open(file, 'r') as input_file:
    # Read all lines from the file
    lines = input_file.readlines()

data = []
bins=180
histrange=[0,10]
x=[]
y=[]

for line in lines:
    x.append(float(line.split()[0]))
    y.append(float(line.split()[1]))

#
fig = plt.figure()
plt.scatter(x,y,)

plt.show()