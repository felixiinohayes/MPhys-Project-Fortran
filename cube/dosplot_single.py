import matplotlib.pyplot as plt
import numpy as np

file = 'data/cube_9_triv_bulk.dat'

# Open the input file in read mode
with open(file, 'r') as input_file:
    # Read all lines from the file
    lines = input_file.readlines()

# Initialize an empty list to hold the data
data = []

# Loop through each line in the file and split into x and y
for line in lines:
    # Split each line by space and convert to float
    x, y = map(float, line.split())
    data.append([x, y])

# Convert the list to a NumPy array
data = np.array(data)

# Create the plot
fig = plt.figure()
plt.plot(data[:, 1], data[:, 0])  # Plot x (column 0) vs y (column 1)
plt.title("Topological Surface eta=0.001")

# Show the plot
plt.show()