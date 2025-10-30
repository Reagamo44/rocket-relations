
import numpy as np

def characteristic_velocity(gamma, R, T0, array_choice):
    """ Calculate the characteristice velocity of a rocket with the following ideal rocket assumptions
    - Chamber conditions created by constant-pressure heating process
    - Working fluid is a non-reacting, calorically perfect gas
    -Nozzle flow is isentropic
    - Flow is choked at nozzle throat
    - Flow is steady and quasi-1D
    - Ignore thermal conductivity and iscous or other body forces in the flow

    Intended to intake numpy arrays and calculate either all variations or treat inputs as paired arrays.

    Parameters
    ----------
    gamma: Specific heat ratio (Cp/Cv), dimensionless
    R: Specific gas constant (J/kg-K)
    T0: Stagnation temperature (K)
    array_choice: 'y' to treat inputs as paired arrays, 'n' to calculate all variations

    Returns
    -------
    char_vel: Characteristic velocity (m/s)
      """


# Calculate exponent to save space
    exponent = (gamma + 1)/(gamma - 1)
    

    if array_choice == "y": # paired arrays function
        char_vel = np.sqrt((1/gamma)*(np.power((gamma + 1)/2, exponent)*R*T0))

# Print results
        for i, c in enumerate(char_vel):
            print(f"γ={gamma[i]}, R={R[i]}, T0={T0[i]} → characteristic velocity = {char_vel[i]:.4f}")

    else:

# Variables to hold expanded arrays to save space
        G = gamma[:, None, None]
        R_1 = R[None, :, None]
        T0_1 = T0[None, None, :]

        char_vel = np.sqrt((1/G)*(np.power((G + 1)/2, exponent[:, None, None])*R_1*T0_1))

# Print results
        for i, g in enumerate(gamma):
            for j, r in enumerate(R):
                for k, t in enumerate(T0):
                    print(f"γ={g}, R={r}, T0={t} → characteristic velocity ={char_vel[i,j,k]:.4f}")

    return char_vel

def thrust_coefficient(gamma, A_ratio, pe_p0, pa_p0, array_choice):
    """ Calculate the thrust coefficient of a rocket with the following ideal rocket assumptions
    - Chamber conditions created by constant-pressure heating process
    - Working fluid is a non-reacting, calorically perfect gas
    -Nozzle flow is isentropic
    - Flow is choked at nozzle throat
    - Flow is steady and quasi-1D
    - Ignore thermal conductivity and iscous or other body forces in the flow

    Parameters
    ----------
    gamma: Specific heat ratio (Cp/Cv), dimensionless
    A_ratio: Nozzle area ratio (Ae/A*), dimensionless
    pe_p0: Nozzle exit pressure ratio (pe/p0), dimensionless
    pa_p0: Ambient pressure ratio (pa/p0), dimensionless
    array_choice: 'y' to treat inputs as paired arrays, 'n' to calculate all variations

    Returns
    -------
    t_coeff: Thrust coefficient, dimensionless
      """

# Calculate exponents to save space
    exponent_1 = (gamma + 1)/(gamma - 1)
    exponent_2 = (gamma - 1)/gamma


    if array_choice == "y": # paired arrays function
        t_coeff = np.sqrt(((2*np.power(gamma, 2))/2)*np.power(2/(gamma + 1), exponent_1)*(1 - np.power(pe_p0, exponent_2))) + (pe_p0 - pa_p0)*A_ratio
        for i, t in enumerate(t_coeff):
            print(f"γ={gamma[i]}, Ae/A* = {A_ratio[i]}, pe/p0 = {pe_p0[i]}, pa/p0 = {pa_p0[i]} → thrust coefficient = {t_coeff[i]:.4f}")
    else:
        # Variables to hold expanded arrays to save space
        G = gamma[:, None, None, None]
        A = A_ratio[None, :, None, None]
        Pe = pe_p0[None, None, :, None] #type error!!!!
        Pa = pa_p0[None, None, None, :]
        E1 = exponent_1[:, None, None, None]
        E2 = exponent_2[:, None, None, None]

        t_coeff = np.sqrt(((2*np.power(G, 2))/2)*np.power(2/(G + 1), E1)*(1 - np.power(Pe, E2))) + (Pe - Pa)*A

# Print results
        for i, g in enumerate(gamma):
            for j, a in enumerate(A_ratio):
                for k, pe in enumerate(pe_p0):
                    for l, pa in enumerate(pa_p0):
                        print(f"γ={g}, Ae/A* = {a}, pe/p0 = {pe}, pa/p0 = {pa} → thrust coefficient ={t_coeff[i,j,k,l]:.4f}")

    return t_coeff
    

def check_input(
    raw_input,
    name,
    min_val=None,
    max_val=None,
    min_inclusive=True,
    max_inclusive=True
):
    """
    Parse and validate comma separated numeric input into a 1-D NumPy array.

    Parameters
    ----------
    raw_input : str or np.ndarray
        Input string (comma separated) or existing numpy array.
    name : str
        Variable name for error messages.
    min_val, max_val : float, optional
        Lower/upper bounds for allowed values.
    min_inclusive, max_inclusive : bool, optional
        Whether the bounds are inclusive.

    Returns
    -------
    np.ndarray
        Cleaned, validated 1-D NumPy array of floats.
    """

    # Convert to np.ndarray
    if isinstance(raw_input, np.ndarray):
        arr = raw_input.astype(float, copy=False)
    elif isinstance(raw_input, str):
        try:
            arr = np.fromstring(raw_input, sep=',', dtype=float)
        except Exception:
            raise TypeError(f"{name} must be comma separated numeric values.")
    else:
        raise TypeError(f"{name} must be a comma separated string or NumPy array.")

    # Basic structure checks
    if arr.size == 0:
        raise ValueError(f"{name} must contain at least one numeric value.")
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1-D; got shape {arr.shape}.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains NaN or infinite values.")

    # Range checks
    if min_val is not None:
        bad = ~(arr >= min_val) if min_inclusive else ~(arr > min_val)
        if np.any(bad):
            cmp = ">=" if min_inclusive else ">"
            raise ValueError(f"{name} must be {cmp} {min_val}. Offending values: {arr[bad]}")
    if max_val is not None:
        bad = ~(arr <= max_val) if max_inclusive else ~(arr < max_val)
        if np.any(bad):
            cmp = "<=" if max_inclusive else "<"
            raise ValueError(f"{name} must be {cmp} {max_val}. Offending values: {arr[bad]}")

    return arr


def safe_input(prompt, **kwargs):
    """ If valid input is not provided, reprompt the user until valid input is given. 
    Parameters
    ----------
    prompt : str
        The input prompt to display to the user.
    **kwargs : dict
        Additional keyword arguments to pass to check_input.
    
    """
    while True:
        try:
            return check_input(input(prompt), **kwargs)
        except ValueError as e:
            print(f"Error: {e}")
        except TypeError as e:
            print(f"Error: {e}")


# User chooses which value to calculate, in future iterations, calculate both at once
calc_choice = input("Would you like to calculate thrust coefficient ('tc') or characteristic velocity (anything else)?: ").strip().lower()

# gamma defenition
gamma = safe_input("Please enter list of specific heat ratios (comma separated): ",
                    name="Specific heat ratio",
                    min_val=1.0,
                    max_val=1.8, min_inclusive=False)
        

if calc_choice != 'tc': # Thrust coefficient

# R definition
    R = safe_input("Please enter list of specific gas constants (comma separated): ",
                    name="Specific gas constant",
                    min_val = 1)

# T0 definition
    T0 = safe_input("Please enter list of stagnation temperatures (comma separated): ",
                    name="Stagnation temperature",
                    min_val = 0,
                    min_inclusive = False)

# Check if arrays are same size and ask user if they want to treat as paired inputs
    if gamma.size == R.size == T0.size:
        array_choice = input("Your arrays are the same size! Would like to treat them as paired inputs rather than calculate all variations? (y/n): ").strip().lower()
        if array_choice not in ['y', 'n']:
            array_choice = input("Invalid input. Please enter 'y' or 'n': ").strip().lower()
    else:
        array_choice = "n"

    characteristic_velocity(gamma, R, T0, array_choice)

else: # Characteristic velocity

# A_ratio definition
    A_ratio = safe_input("Please enter list of nozzle area ratios (Ae/A*), comma separated: ",
                    name="Nozzle area ratio",
                    min_val = 0,
                    min_inclusive = False)

# pe_p0 definition
    pe_p0 = safe_input("Please enter list of nozzle exit pressure ratios (pe/p0), comma separated: ",
                    name="Nozzle exit pressure ratio",
                    min_val = 0,
                    max_val = 1,
                    max_inclusive = False)
# pa_p0 definition
    pa_p0 = safe_input("Please enter list of ambient pressure ratios (pa/p0), comma separated: ",
                    name="Ambient pressure ratio",
                    min_val = 0,
                    max_val = 1,
                    max_inclusive= False)

# Check if arrays are same size and ask user if they want to treat as paired inputs
    if gamma.size == A_ratio.size == pe_p0.size == pa_p0.size:
        array_choice = input("Your arrays are the same size! Would like to treat them as paired inputs rather than calculate all variations? (y/n): ").strip().lower()
    else:
        array_choice = "n"

    thrust_coefficient(gamma, A_ratio, pe_p0, pa_p0, array_choice)