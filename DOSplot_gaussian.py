import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# Load the eigenvalues from the generated file
eigenvalues = np.loadtxt('cube_eigenvalues_7.dat')

# Set up the energy bins for the DOS calculation
emin, emax, eres = 3.6, 4.8, 70
energy_bins = np.linspace(emin, emax, eres)

# Initialize DOS
dos = np.zeros(len(energy_bins))

# Gaussian broadening parameters
eta1 = 0.2  # Width of the Gaussian, similar to eta1 in the Fortran code

# Calculate the DOS using Gaussian broadening
for eig in eigenvalues:
    factor = (energy_bins - eig) / eta1
    dos += np.exp(-0.5 * factor**2) / (np.sqrt(2 * np.pi) * eta1)

# Optionally, apply Gaussian smoothing to the DOS for a smoother output
# dos_smoothed = gaussian_filter1d(dos, sigma=1.0)

# Plot the Density of States (DOS)
plt.figure(figsize=(8, 6))
plt.plot(energy_bins, dos, label="Gaussian Blurred DOS", color='b')
plt.xlabel('Energy (eV)')
plt.ylabel('Density of States')
plt.title('Density of States with Gaussian Broadening')
plt.legend()
plt.grid(True)
plt.show()