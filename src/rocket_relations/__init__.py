import numpy as np

""" 
Imports two main functions from the ideal rocket model: chraacteristic velocity, and thrust coefficient. The two make the following assumptions:
    - Chamber conditions created by constant-pressure heating process
    - Working fluid is a non-reacting, calorically perfect gas
    -Nozzle flow is isentropic
    - Flow is choked at nozzle throat
    - Flow is steady and quasi-1D
    - Ignore thermal conductivity and viscous or other body forces in the flow

Both print results to the console and return the calculated arrays.

Characteristic Velocity:
    Parameters
    ----------
    gamma: Specific heat ratio (Cp/Cv), dimensionless
    R: Specific gas constant (J/kg-K)
    T0: Stagnation temperature (K)
    array_choice (only change if input arrays are the same size!): 'y' to treat inputs as paired arrays, 'n' to calculate all variations

    Returns
    -------
    char_vel: Characteristic velocity (m/s)


Thrust Coefficient:

    Parameters
    ----------
    gamma: Specific heat ratio (Cp/Cv), dimensionless
    A_ratio: Nozzle area ratio (Ae/A*), dimensionless
    pe_p0: Nozzle exit pressure ratio (pe/p0), dimensionless
    pa_p0: Ambient pressure ratio (pa/p0), dimensionless
    array_choice (only change if input arrays are the same size!): 'y' to treat inputs as paired arrays, 'n' to calculate all variations

    Returns
    -------
    t_coeff: Thrust coefficient, dimensionless

"""

from .ideal import characteristic_velocity, thrust_coefficient

__all__ = ["characteristic_velocity", "thrust_coefficient"]