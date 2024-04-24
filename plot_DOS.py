import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.colors import LightSource
import numpy as np

file = 'DOS.dat'

# Open the input file in read mode
with open(file, 'r') as input_file:
    # Read all lines from the file
    lines = input_file.readlines()

data = []
col = []
bins=180
histrange=[0,10]

for line in lines:
    data.append(float(line.split()[0]))
    # col.append(float(line.split()[1]))
    # for val in line.split():
    #     val = float(val)
    #     data.append(val)

fig = plt.figure()
plt.hist(data,bins,histrange,density=True)

plt.show()