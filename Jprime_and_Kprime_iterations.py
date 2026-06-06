import numpy as np
import pandas as pd
from scipy.integrate import dblquad
from scipy.special import eval_genlaguerre

def R_ns(r, n, Z=3):
    """
    Calculates the radial wavefunction.
    r: Radial distance from the nucleus.
    n: Principal quantum number.
    Z: Atomic number.
    eval_genlaguerre: Generates Generalized Laguerre Polynomials to map the radial nodes mathematically.
    """
    if n == 1:
        return 2 * (Z**1.5) * np.exp(-Z * r)
    elif n == 2:
        return (1 / np.sqrt(2)) * (Z**1.5) * (1 - Z * r / 2) * np.exp(-Z * r / 2)
    else:
        rho = 2 * Z * r / n
        norm = (2 * (Z**1.5)) / (n**2.5)
        poly = eval_genlaguerre(n - 1, 1, rho)
        return norm * poly * np.exp(-Z * r / n)

def multipole_term(r1, r2):
    """
    Evaluates the multipole expansion term for spherical symmetry (1 / r_greater).
    r1, r2: Radial coordinates of the first and second interacting electrons.
    """
    r_less, r_greater = min(r1, r2), max(r1, r2)
    return 1 / r_greater if r_greater > 0 else 0

def compute_integral(n):
    """
    Calculates the spatial transition integrals J_prime and K_prime.
    dblquad: SciPy's adaptive quadrature algorithm for double numerical integration.
    lambda: Defines anonymous, inline functions for the integrand and spatial boundaries.
    np.inf: Represents mathematical infinity to integrate over the entire space.
    epsabs, epsrel: Parameters setting the absolute and relative error tolerances (1e-5) for algorithmic convergence.
    """
    J_val, _ = dblquad(
        lambda r2, r1: R_ns(r1, n) * R_ns(r1, 2) * multipole_term(r1, r2) * R_ns(r2, 1) * R_ns(r2, 1) * (r1**2) * (r2**2),
        0, np.inf, lambda r1: 0, lambda r1: np.inf, epsabs=1e-5, epsrel=1e-5
    )
    
    K_val, _ = dblquad(
        lambda r2, r1: R_ns(r1, n) * R_ns(r1, 1) * multipole_term(r1, r2) * R_ns(r2, 1) * R_ns(r2, 2) * (r1**2) * (r2**2),
        0, np.inf, lambda r1: 0, lambda r1: np.inf, epsabs=1e-5, epsrel=1e-5
    )
    return J_val, K_val

hartree_to_ev = 27.2114
E_0 = -275.51 

print(f"{'Channel':<10} | {'J_prime (eV)':<12} | {'K_prime (eV)':<12} | {'Delta_E (eV)':<12} | {'E_2_Contrib (eV)':<15}")
print("-" * 75)

total_E2 = 0

for n in range(3, 7):
    """
    Iterates through the virtual state channels (3s, 4s, 5s, 6s).
    delta_E: Calculates the energy denominator required for the perturbation expansion.
    E2: Computes the specific correlation energy contribution for the given orbital.
    """
    J, K = compute_integral(n)
    
    E_uncorrected = 2 * (-122.45) - 13.6057 * (9 / (n**2))
    delta_E = E_0 - E_uncorrected
    
    E2 = ((2 * J - K)**2 * hartree_to_ev**2) / delta_E
    total_E2 += E2
    
    print(f"2s -> {n}s  | {J * hartree_to_ev:<12.3f} | {K * hartree_to_ev:<12.3f} | {delta_E:<12.3f} | {E2:<15.3f}")

print("-" * 75)
print(f"Total Second-Order Energy Contribution: {total_E2:.3f} eV")
