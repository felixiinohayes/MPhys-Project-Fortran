import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

# Load the data from the CSV file
df = pd.read_csv('runtime.csv')
df = df.drop(df.index[-1])

# Extract columns as numpy arrays
nev = df['nev'].to_numpy()
nblocks = df['nblocks'].to_numpy()
nprocs = df['nprocs'].to_numpy()
runtime = df['runtime'].to_numpy()

nev = nev / 100
nblocks = nblocks / 10

def model(variables, a, b, c, d):
    nblocks = variables[0]
    nev = variables[1]
    nprocs = variables[2]
    return a*nev + b*(nblocks**c)/nprocs + d

# Use curve fitting to find the best a and b
params, covariance = curve_fit(model, (nblocks, nev, nprocs), runtime)
a, b, c, d = params

# Print the best-fit parameters
print(f"a = {a:.4f}")
print(f"b = {b:.4f}")
print(f"c = {c:.4f}")
print(f"d = {d:.4f}")
# print(f"e = {e:.4f}")
# print(f"f = {f:.4f}")

# Calculate model predictions using the fitted parameters
predicted_times = model((nblocks, nev, nprocs), a, b, c, d)

# Calculate Mean Square Error (MSE)
mse = np.mean((runtime - predicted_times) ** 2)
print(f"Mean Square Error = {mse:.8f}")