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
