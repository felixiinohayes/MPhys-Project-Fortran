import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

blocks = np.array([7, 10, 15, 20])
nev = np.array([120, 100, 300, 700])

log_nev = np.log(nev)
print(log_nev)
slope, intercept, r_value, p_value, std_err = linregress(blocks, log_nev)
print("Slope:", slope)
print("Intercept:", intercept)

blocks_extrap = np.linspace(blocks.min(), 100, 200)
log_nev_extrap = slope * blocks_extrap + intercept

print(np.exp(slope * 30 + intercept))

# Create the plot
plt.figure(figsize=(8, 6))
plt.scatter(blocks, log_nev, color='blue', label='Data points')
plt.plot(blocks_extrap, log_nev_extrap, color='red', label='Linear fit & Extrapolation')
plt.xlabel("Blocks")
plt.ylabel("log(nev)")
plt.title("Plot of log(nev) vs Blocks with Extrapolation")
plt.legend()
plt.grid(True)
plt.show()