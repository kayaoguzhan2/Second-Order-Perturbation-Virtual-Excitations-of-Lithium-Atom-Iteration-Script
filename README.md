In this repository, `Jprime_and_Kprime_Iterations.py` script is designed to calculate the dynamical electron correlation energy of the Lithium atom ($Z=3$) ground state, utilizing Second-Order Perturbation Theory. Specifically, it evaluates the spatial transition integrals for single-electron virtual excitations ($2s \to 3s, 4s, 5s$, and $6s$) by applying the `scipy.integrate.dblquad` adaptive quadrature algorithm. The `Radial_Probability_Densities.py` script is designed to plot radial probability densities for each perturbation correction and for the Variational Method. 

### J' and K' Iteration Script

### Adaptive Quadrature Algorithm

The scipy.integrate.dblquad function utilizes an adaptive quadrature algorithm that dynamically refines its integration step size in regions where the mathematical function changes rapidly. This method was specifically chosen because it evaluates the entire spatial domain from zero to infinity, completely eliminating the truncation and discretization errors inherent to fixed-grid approaches. As a result, it seamlessly captures both the ultra-sharp probability density peak of the $1s$ core electron near the nucleus and the highly diffuse exponential tails of upper virtual orbitals simultaneously.

### Computational Approach

This computational approach executes a double numerical integration across the entire spatial domain (from zero to infinity, `np.inf`), dynamically adjusting its mathematical step size to ensure a better precision. Consequently, it achieves high precision and accurately models complex quantum mechanical behaviors. These include effectively capturing the diffuse exponential tails of higher virtual states (such as the $6s$ orbital) alongside the highly concentrated probability density peak of the $1s$ core electron.

### Mathematical Variables

* **$r, r_1, r_2$**: Radial distances of the interacting electrons from the nucleus.
* **$n$**: Principal quantum number ($3, 4, 5, 6$) that determines the energy level and spatial spread of the virtual orbital.
* **$Z$**: Atomic number (Nuclear charge).
* **$J'$** (Coulomb Integral): Represents the electrostatic repulsion between the $1s$ core and the virtual transition density.
* **$K'$** (Exchange Integral): Denotes the quantum mechanical spatial exchange interaction resulting from the Pauli Exclusion Principle.
* **$\Delta E$** (Energy Denominator): The zeroth-order energy difference between the ground and excited configurations.
* **$E^{(2)}$**: The final second-order energy contribution, computed via the Rayleigh-Schrödinger formula.
  
#### Wavefunction Visualization and Perturbation Analysis

This part of the script shows the radial probability density of the Lithium atom's ground and perturbed states. It visually explains how electron correlation and screening effects change the spatial distribution of the wavefunctions.

## Purpose and Physical Mechanism

The main goal is to visually compare the independent-particle model (zeroth-order baseline) with the correlated model (second-order perturbation theory) and the optimized variational model. 

In the independent-particle baseline, we assume the electrons feel the full nuclear charge ($Z=3$) and we completely ignore the repulsion between electrons. To include this dynamical correlation in our physics model, the second-order perturbation approach mixes higher-energy virtual orbitals ($3s, 4s, 5s, 6s$) into the baseline $2s$ state. This mixing changes the shape of the wavefunction to reduce the electron-electron repulsion. We then compare this new distribution to the variational method, which uses effective nuclear charges to simulate how inner-shell electrons shield the nucleus.

## Computational Approach

The script creates a spatial grid starting from the nucleus up to $8$ Bohr radii. It calculates the analytical hydrogenic wavefunctions, creates a linear combination using the calculated first-order perturbation coefficients, and uses numerical integration to re-normalize the new state. Finally, it plots all probability curves together to show how the orbitals relax and contract.

## Mathematical Variables and Functions

* **$R_{ns}(n, Z, r)$**: The radial part of the hydrogenic wavefunction. It shows the probability amplitude of finding an electron at a radial distance $r$. Here, $n$ is the principal quantum number and $Z$ is the nuclear charge.
* **$P(r)$**: The radial probability density. It describes the real probability of finding the electron inside a spherical shell of radius $r$ and thickness $dr$. The formula is:
  $$P(r) = r^2 |R_{ns}(r)|^2$$
* **$c$** (`c_coeffs`): The first-order correction coefficients. They show how much each virtual $ns$ state mixes into the $2s$ state to lower the Coulomb repulsion.
* **$Z_{eff}$**: The effective nuclear charge. The variational approach uses it to model how inner-shell electrons shield the electrostatic pull of the nucleus ($Z_{eff} = 2.68$ for the $1s$ orbital, $Z_{eff} = 1.87$ for the $2s$ orbital).

