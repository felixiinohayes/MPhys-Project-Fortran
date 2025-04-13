import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

nb = []
put = []
pupp = []

with open('./mem_usage.csv', 'r') as f:
    lines = f.readlines()
    # Skip header row if it exists
    for line in lines:
        parts = line.strip().split(',')
        try:
            nb.append(int(parts[0]))
            put.append(float(parts[1]))
            pupp.append(float(parts[2]))
        except ValueError:
            # Skip lines that can't be converted (like headers)
            continue

# Convert to numpy arrays after collecting all data
nb = np.array(nb)
put = np.array(put)
pupp = np.array(pupp)

# Define inverse cubic function with non-negative constraint: f(x) = a/(x^3) + b, where b >= 0
def inverse_cubic(x, a, b):
    # Ensure b is non-negative to prevent values below 0
    return a / (x**3) + abs(b)

# Fit the model to the data
put_params, put_cov = curve_fit(inverse_cubic, nb, put)
pupp_params, pupp_cov = curve_fit(inverse_cubic, nb, pupp)

# Create a range of nb values for extrapolation
max_nb = np.max(nb)
extrapolated_nb = np.linspace(max_nb, max_nb*5, 100)  # Extrapolating to 5x the max value

# Calculate extrapolated values
extrapolated_put = inverse_cubic(extrapolated_nb, *put_params)
extrapolated_pupp = inverse_cubic(extrapolated_nb, *pupp_params)

# Print model parameters
print("PUT model: a = {:.4f}, b = {:.4f}".format(put_params[0], put_params[1]))
print("PUPP model: a = {:.4f}, b = {:.4f}".format(pupp_params[0], pupp_params[1]))

# Plot original data
plt.figure(figsize=(10, 6))
plt.scatter(nb, put, label='put (original)', marker='o')
plt.scatter(nb, pupp, label='pupp (original)', marker='s')

# Plot fitted curves for the original range
x_fit = np.linspace(np.min(nb), np.max(nb), 100)
plt.plot(x_fit, inverse_cubic(x_fit, *put_params), 'b-', label='put fit')
plt.plot(x_fit, inverse_cubic(x_fit, *pupp_params), 'r-', label='pupp fit')

# Plot extrapolated values
plt.plot(extrapolated_nb, extrapolated_put, 'b--', label='put extrapolated')
plt.plot(extrapolated_nb, extrapolated_pupp, 'r--', label='pupp extrapolated')

plt.xlabel('nb')
plt.ylabel('Memory Usage (%)')
plt.ylim(bottom=0,top=100)  # Ensure y-axis starts at 0
plt.title('Memory Usage Comparison with Inverse Cubic Extrapolation')
plt.legend()
plt.grid(True)
plt.show()

# Print extrapolated values at specific points
print("\nExtrapolated values at key points:")
test_points = [max_nb, max_nb*2, max_nb*3, max_nb*4, max_nb*5]
for point in test_points:
    put_val = inverse_cubic(point, *put_params)
    pupp_val = inverse_cubic(point, *pupp_params)
    print(f"nb = {point:.0f}: PUT = {put_val:.4f}%, PUPP = {pupp_val:.4f}%")





