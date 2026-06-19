This repository features a Python script designed to calculate the dynamical electron correlation energy of the Lithium atom ($Z=3$) ground state, utilizing Second-Order Perturbation Theory. Specifically, it evaluates the spatial transition integrals for single-electron virtual excitations ($2s \to 3s, 4s, 5s$, and $6s$) by applying the `scipy.integrate.dblquad` adaptive quadrature algorithm.

## Adaptive Quadrature Algorithm

The scipy.integrate.dblquad function utilizes an adaptive quadrature algorithm that dynamically refines its integration step size in regions where the mathematical function changes rapidly. This method was specifically chosen because it evaluates the entire spatial domain from zero to infinity, completely eliminating the truncation and discretization errors inherent to fixed-grid approaches. As a result, it seamlessly captures both the ultra-sharp probability density peak of the $1s$ core electron near the nucleus and the highly diffuse exponential tails of upper virtual orbitals simultaneously.

## Computational Approach

This computational approach executes a double numerical integration across the entire spatial domain (from zero to infinity, `np.inf`), dynamically adjusting its mathematical step size to ensure a better precision. Consequently, it achieves high precision and accurately models complex quantum mechanical behaviors. These include effectively capturing the diffuse exponential tails of higher virtual states (such as the $6s$ orbital) alongside the highly concentrated probability density peak of the $1s$ core electron.

## Mathematical Variables

* **$r, r_1, r_2$**: Radial distances of the interacting electrons from the nucleus.
* **$n$**: Principal quantum number ($3, 4, 5, 6$) that determines the energy level and spatial spread of the virtual orbital.
* **$Z$**: Atomic number (Nuclear charge).
* **$J'$** (Coulomb Integral): Represents the electrostatic repulsion between the $1s$ core and the virtual transition density.
* **$K'$** (Exchange Integral): Denotes the quantum mechanical spatial exchange interaction resulting from the Pauli Exclusion Principle.
* **$\Delta E$** (Energy Denominator): The zeroth-order energy difference between the ground and excited configurations.
* **$E^{(2)}$**: The final second-order energy contribution, computed via the Rayleigh-Schrödinger formula.
## Wavefunction Visualization and Perturbation Analysis

This section of the script visualizes the radial probability density of the Lithium atom's ground and perturbed states, illustrating how electron correlation and screening modify the spatial distribution of the wavefunctions.

### Purpose and Physical Mechanism

The primary objective is to visually contrast the independent-particle model (zeroth-order baseline) against the correlated model (second-order perturbation theory) and the optimized variational model. In the independent-particle baseline, electrons are assumed to experience an unscreened nuclear charge ($Z=3$), completely neglecting inter-electronic repulsion. 

To physically account for this dynamical correlation, the second-order perturbation approach admixes higher-energy virtual orbital configurations ($3s, 4s, 5s, 6s$) into the baseline $2s$ state. This expansion deforms the wavefunction to minimize electron-electron repulsion. This perturbed distribution is graphically compared against the variational method, which uses effective nuclear charges to simulate static inner-shell screening effects.

### Computational Approach

The script generates a precise radial grid extending from the nucleus outwards to $8$ Bohr radii. It evaluates the analytic hydrogenic wavefunctions, executes a linear combination weighted by the calculated first-order perturbation coefficients, and applies numerical integration to re-normalize the perturbed state. Finally, it maps all probability profiles simultaneously to demonstrate orbital relaxation and contraction.

### Mathematical Variables and Functions

* **$R_{ns}(n, Z, r)$**: The radial part of the hydrogenic wavefunction. It represents the probability amplitude of finding an electron at a radial distance $r$ for a given principal quantum number $n$ and nuclear charge $Z$.
* **$P(r)$**: The radial probability density, defined as $P(r) = r^2 |R_{ns}(r)|^2$. It describes the actual probability of finding the electron within a spherical shell of radius $r$ and thickness $dr$.
* **$c$** (`c_coeffs`): The first-order correction coefficients to the wavefunction. They dictate the amplitude weight (admixture) of each virtual $ns$ state mixing into the perturbed $2s$ state to minimize Coulomb repulsion.
* **$Z_{eff}$**: The effective nuclear charge. Used in the variational approach to model the electrostatic shielding of the nucleus caused by inner-shell electrons ($Z_{eff} = 2.68$ for $1s$, $Z_{eff} = 1.87$ for $2s$).
