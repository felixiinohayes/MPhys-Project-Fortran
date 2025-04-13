import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.colors import LightSource
import numpy as np

file = '../supercomputer/eigenvalues.dat'

# Open the input file in read mode
with open(file, 'r') as input_file:
    # Read all lines from the file
    lines = input_file.readlines()

data = []
bins=1000
histrange=[-1,6]
count = 0
for line in lines:
    value = float(line.split()[0])
    if value != 0:
        data.append(value)
        if abs(value) <0.3:
            count += 1
    # col.append(float(line.split()[1]))
    # for val in line.split():
    #     val = float(val)
    #     data.append(val)

print(count)

fig = plt.figure()
plt.xticks(np.arange(histrange[0], histrange[1], 0.1))
plt.xlim(min(data),max(data))
plt.hist(data,bins,histrange,density=False)

plt.show()
