import numpy as np
import matplotlib.pyplot as plt
from scipy.special import eval_genlaguerre
from scipy.integrate import simpson

def R_ns(n, Z, r):
    # Calculates the radial part of the wavefunction for a hydrogen-like atom
    return (1.0 / (np.sqrt(2) * n)) * ((2 * Z / n)**1.5) * eval_genlaguerre(n - 1, 1, 2 * Z * r / n) * np.exp(-Z * r / n)

# Create a spatial range from the nucleus out to 8 Bohr radii
r = np.linspace(0, 8, 1000)
Z = 3  # The atomic number of Lithium

# Calculate the probability density for the 1s and 2s orbitals 
P1s_01 = (r**2) * (R_ns(1, Z, r)**2)
P2s_01 = (r**2) * (R_ns(2, Z, r)**2)

# These coefficients determine how much the higher energy virtual states (3s, 4s, etc.)
c_coeffs = [(2*4.125 - 0.916) / -17.004, (2*2.344 - 0.582) / -22.957, 
            (2*1.590 - 0.413) / -25.712, (2*1.176 - 0.312) / -27.209]

# Create the new perturbed 2s orbital by mixing the baseline with virtual states
R2s_2nd = R_ns(2, Z, r) + sum(c * R_ns(n, Z, r) for c, n in zip(c_coeffs, [3,4,5,6]))

# Re-normalize the new function so the total probability remains equal to 1
R2s_2nd /= np.sqrt(simpson(y=(r**2)*(R2s_2nd**2), x=r))
P2s_2nd = (r**2) * (R2s_2nd**2)

# Calculate the densities using the Zeff found in the Variational Method section
P1s_lit = (r**2) * (R_ns(1, 2.68, r)**2)
P2s_lit = (r**2) * (R_ns(2, 1.87, r)**2)


plt.figure(figsize=(12, 6))

plt.xlabel("Distance (a0)", fontweight='bold')
plt.ylabel("P(r)", fontweight='bold')
plt.title("Radial Probability Density", fontweight='bold')

plt.plot(r, P1s_01, label="Zeroth/First Order 1s", color="purple", linestyle="--")
plt.plot(r, P1s_lit, label="Variational Parameter Fit of 1s (Zeff=2.68)", color="orange", linewidth=2)

plt.plot(r, P2s_01, label="Zeroth/First Order 2s", color="cyan", linestyle="-.")
plt.plot(r, P2s_2nd, label="Second Order 2s", color="green", linestyle=":")
plt.plot(r, P2s_lit, label="Variational Parameter Fit of 2s (Zeff=1.87)", color="red", linewidth=2)

plt.legend(prop={'weight': 'bold'}) 
plt.grid(True, alpha=0.3)
plt.show()
