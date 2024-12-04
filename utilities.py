import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import sympy as sp
import cantera as ct


# Constants
R = 8.314  # Universal gas constant, J/(mol·K)
T0 = 298.15  # Standard temperature, K
P = 101325  # Standard pressure (1 atm in Pa)
df = pd.read_excel('SP_298K.xlsx')  # Load data from an Excel file containing thermodynamic properties

def standard_properties(reactants, products):
    # Create explicit copies of reactants and products' properties from the data frame
    reactants_properties_298K = df[df['Species'].isin(reactants[:, 0])].copy()
    products_properties_298K = df[df['Species'].isin(products[:, 0])].copy()

    # Add columns for coefficients and weighted properties
    reactants_properties_298K["Coefficient"] = reactants[:, 1].astype(float)
    products_properties_298K["Coefficient"] = products[:, 1].astype(float)

    # Calculate weighted properties for reactants
    reactants_properties_298K["Weighted_DHf"] = reactants_properties_298K["DHf°[kJ/mol]"] * reactants_properties_298K["Coefficient"]
    reactants_properties_298K["Weighted_DGf"] = reactants_properties_298K["DGf°[kJ/mol]"] * reactants_properties_298K["Coefficient"]
    reactants_properties_298K["Weighted_S"] = reactants_properties_298K["S°[J/K·mol]"] * reactants_properties_298K["Coefficient"]

    # Calculate weighted properties for products
    products_properties_298K["Weighted_DHf"] = products_properties_298K["DHf°[kJ/mol]"] * products_properties_298K["Coefficient"]
    products_properties_298K["Weighted_DGf"] = products_properties_298K["DGf°[kJ/mol]"] * products_properties_298K["Coefficient"]
    products_properties_298K["Weighted_S"] = products_properties_298K["S°[J/K·mol]"] * products_properties_298K["Coefficient"]

    # Compute the weighted sums for products and reactants
    dH_f0 = products_properties_298K["Weighted_DHf"].sum() - reactants_properties_298K["Weighted_DHf"].sum()
    dG_f0 = products_properties_298K["Weighted_DGf"].sum() - reactants_properties_298K["Weighted_DGf"].sum()

    # Option to calculate ΔSº directly from the table
    # dS0 = products_properties_298K["Weighted_S"].sum() - reactants_properties_298K["Weighted_S"].sum()

    # Assuming ΔHº and ΔSº are temperature-independent, calculate ΔSº using the relation ΔGº = ΔHº - TΔSº
    dS0 = (dG_f0 - dH_f0) / T0

    # Print the filtered rows for debug or validation
    print("Reactants Properties")
    print(reactants_properties_298K)
    print("\n")
    print("Products Properties")
    print(products_properties_298K)
    print("\n")

    # Output results with a qualitative description
    print(f"ΔGº: {dG_f0} [kJ/mol], and it {'is' if dG_f0 > 0 else 'is not'} spontaneous.")
    print(f"ΔHº: {dH_f0} [kJ/mol], and it {'is endothermic.' if dH_f0 > 0 else 'is exothermic.'}")
    print(f"ΔSº: {dS0} [J/K·mol], and the disorder {'increases.' if dH_f0 > 0 else 'decreases.'}")

    # Return calculated values and data frames for further use
    return dG_f0, dH_f0, dS0, reactants_properties_298K, products_properties_298K

def properties_temperature(temperatures, dG_f0, dS0):
    # Calculation of ΔGº for each temperature
    # Using the relation: ΔGº(T2) = ΔGº(T1) + ΔSº(T2 - T1)
    dG_values = dG_f0 + dS0 * (temperatures - T0)

    # Calculation of K_a and ln(K_a)
    # ln(K_a) = -ΔGº / (R * T) and K_a = exp(ln(K_a))
    ln_Ka_values = -dG_values / (R * temperatures)
    Ka_values = np.exp(ln_Ka_values)

    # Create a DataFrame with the results
    results = pd.DataFrame({
        'T (K)': temperatures,             # Temperature in Kelvin
        'ΔGº (kJ/mol)': dG_values,        # Gibbs free energy change at each temperature
        'ln(Ka)': ln_Ka_values,           # Natural logarithm of the equilibrium constant
        'Ka': Ka_values                   # Equilibrium constant
        # 'Xeq': Xeq_values                # Uncomment if equilibrium composition is added
    })
    
    # Print the results table
    print(results)
    
    # Return the results as a DataFrame
    return results

def convert_equilibrium_constants(K_type, K_value, T=None, delta_n=0, P_total=None, C_total=None):
    """
    Convierte entre diferentes constantes de equilibrio: Kc, Kp, Ky, Kx, Ka.
    
    Parámetros:
        K_type (str): El tipo de constante dada ('c', 'p', 'y', 'x', 'a').
        K_value (float): El valor de la constante.
        T (float): Temperatura en Kelvin (obligatoria para Kp y Ky).
        delta_n (float): Cambio en moles gaseosos (productos - reactivos).
        P_total (float): Presión total (obligatoria para Ky, Ka y Kp si es relevante).
        C_total (float): Concentración total en el sistema (obligatoria para Ky y Kx si es relevante).
    
    Retorna:
        dict: Equivalencias entre las constantes de equilibrio.
    """
    R = 0.08206  # Constante de los gases en L·atm·mol^-1·K^-1
    
    Kc = Kp = Ky = Kx = Ka = None

    if K_type.lower() == 'c':  # Kc está dado
        Kc = K_value
        if T is not None:
            Kp = Kc * (R * T) ** delta_n
        if P_total is not None and C_total is not None:
            Ky = Kc * (C_total / P_total) ** delta_n
            Kx = Ky  # Si no hay diferencias específicas para Kx
        if Ky is not None and P_total is not None:
            Ka = Ky * P_total

    elif K_type.lower() == 'p':  # Kp está dado
        Kp = K_value
        if T is not None:
            Kc = Kp / (R * T) ** delta_n
        if P_total is not None:
            Ky = Kp / (P_total ** delta_n)
            Kx = Ky  # Si no hay diferencias específicas para Kx
        if Ky is not None and P_total is not None:
            Ka = Ky * P_total

    elif K_type.lower() == 'y':  # Ky está dado
        Ky = K_value
        if P_total is not None:
            Kp = Ky * P_total ** delta_n
        if C_total is not None and P_total is not None:
            Kc = Ky * (P_total / C_total) ** delta_n
        if P_total is not None:
            Ka = Ky * P_total
        Kx = Ky  # Si no hay diferencias específicas para Kx

    elif K_type.lower() == 'x':  # Kx está dado
        Kx = K_value
        Ky = Kx  # Asumiendo equivalencia entre Kx y Ky
        if P_total is not None:
            Kp = Ky * P_total ** delta_n
        if C_total is not None and P_total is not None:
            Kc = Ky * (P_total / C_total) ** delta_n
        if P_total is not None:
            Ka = Ky * P_total

    elif K_type.lower() == 'a':  # Ka está dado
        Ka = K_value
        if P_total is not None:
            Ky = Ka / P_total
            Kx = Ky  # Si no hay diferencias específicas para Kx
            if T is not None:
                Kp = Ky * P_total ** delta_n
            if C_total is not None:
                Kc = Ky * (P_total / C_total) ** delta_n

    else:
        raise ValueError("K_type debe ser 'c' (Kc), 'p' (Kp), 'y' (Ky), 'x' (Kx), o 'a' (Ka).")
    
    return {"Kc": Kc, "Kp": Kp, "Ky": Ky, "Kx": Kx, "Ka": Ka}



