import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import sympy as sp
import cantera as ct
from sympy import symbols, Eq, solve
import math

# Constants
R = 8.314  # Universal gas constant in J/(mol·K)
T0 = 298.15  # Standard temperature in K
P = 1  # Standard atmospheric pressure in atm
# data source: https://www2.chem.wisc.edu/deptfiles/genchem/netorial/modules/thermodynamics/table.htm
df = pd.read_excel('SP_298K.xlsx')  # Load thermodynamic properties from an Excel file
x = symbols('x')

def standard_properties(reactants, products, P = 1, P0 = 1):
    """
    Calculate standard thermodynamic properties (ΔH°, ΔG°, ΔS°) at 298 K for a reaction.

    Parameters:
        reactants (list of tuples): Reactant species and their stoichiometric coefficients.
        products (list of tuples): Product species and their stoichiometric coefficients.

    Returns:
        tuple: ΔG°, ΔH°, ΔS°, overall stoichiometric coefficient, and data frames for reactants and products.
    """
    # Map reactants and products to their stoichiometric coefficients
    r = {reactant[0]: float(reactant[1]) for reactant in reactants}
    p = {product[0]: float(product[1]) for product in products}

    # Filter properties for reactants and products based on the species names
    r_298K = df[df['Species'].isin(r.keys())].copy()
    p_298K = df[df['Species'].isin(p.keys())].copy()

    # Assign stoichiometric coefficients to the respective species
    r_298K["Coefficient"] = r_298K["Species"].map(r)
    p_298K["Coefficient"] = p_298K["Species"].map(p)

    # Compute the total stoichiometric coefficient
    vi = p_298K["Coefficient"].sum() - r_298K["Coefficient"].sum()

    # Calculate weighted properties for reactants
    r_298K["Weighted_DHf"] = r_298K["DHf°[kJ/mol]"] * r_298K["Coefficient"]
    r_298K["Weighted_DGf"] = r_298K["DGf°[kJ/mol]"] * r_298K["Coefficient"]
    r_298K["Weighted_S"] = r_298K["S°[J/K·mol]"] * r_298K["Coefficient"]

    # Calculate weighted properties for products
    p_298K["Weighted_DHf"] = p_298K["DHf°[kJ/mol]"] * p_298K["Coefficient"]
    p_298K["Weighted_DGf"] = p_298K["DGf°[kJ/mol]"] * p_298K["Coefficient"]
    p_298K["Weighted_S"] = p_298K["S°[J/K·mol]"] * p_298K["Coefficient"]

    # Compute overall property changes (ΔH°, ΔG°, ΔS°)
    dH_f0 = p_298K["Weighted_DHf"].sum() - r_298K["Weighted_DHf"].sum()
    dG_f0 = p_298K["Weighted_DGf"].sum() - r_298K["Weighted_DGf"].sum()

    # Calculate ΔS° using the relationship ΔG° = ΔH° - TΔS°
    dS0 = (dH_f0 - dG_f0) / T0

    ln_Ka_0 = -dG_f0 * 1000 / (R * T0)  # Convert ΔG to J/mol for consistency
    Ka_0 = np.exp(ln_Ka_0)

    # Mass balance using conversion
    # Moles per component
    r_298K["n [mol]"] = r_298K["Coefficient"] * (1-x)
    p_298K["n [mol]"] = p_298K["Coefficient"] * x
    # Total number of moles
    n_total = r_298K["n [mol]"].sum() + p_298K["n [mol]"].sum()
    # Mass fraction of each component
    r_298K["y_i"] = r_298K["n [mol]"]/n_total
    p_298K["y_i"] = p_298K["n [mol]"]/n_total
    # Calculate Ky
    K_y = ((r_298K["y_i"]**r_298K["Coefficient"]).prod())/((p_298K["y_i"]**p_298K["Coefficient"]).prod()) 
    K_y = K_y *(P/P0)**-vi

    # Display detailed reactant and product properties for debugging
    print("Reactant Properties:\n", r_298K, "\n")
    print("Product Properties:\n", p_298K, "\n")

    # Print thermodynamic results with descriptions
    print(f"ΔG°: {dG_f0} kJ/mol, {'non-spontaneous' if dG_f0 > 0 else 'spontaneous'}.")
    print(f"ΔH°: {dH_f0} kJ/mol, {'endothermic' if dH_f0 > 0 else 'exothermic'}.")
    print(f"ΔS°: {dS0} kJ/(K·mol), {'increased disorder' if dS0 > 0 else 'decreased disorder'}.")
    print(f"ln(K_a): {ln_Ka_0}")
    print(f"K_a: {Ka_0}")
    print(f"K_y: {K_y}")
    print(f"Net Stoichiometric Coefficient (v_i): {vi}")
    print(f"Total number of moles (n_t): {n_total}")

    return dG_f0, dH_f0, dS0, vi, K_y, r_298K, p_298K

def properties_temperature(temperatures, dH_f0, dS0, K_y = None, Convertion = 0):
    """
    Compute temperature-dependent properties (ΔG° and equilibrium constants).

    Parameters:
        temperatures (array): Array of temperatures in Kelvin.
        dH_f0 (float): Standard enthalpy change at 298 K in kJ/mol.
        dS0 (float): Standard entropy change at 298 K in J/(K·mol).

    Returns:
        DataFrame: Temperature, ΔG°, ln(K_a), and K_a values.
    """
    # Calculate ΔG° at each temperature
    dG_values = dH_f0 - dS0 * temperatures

    # Compute ln(K_a) and equilibrium constants K_a
    ln_Ka_values = -dG_values * 1000 / (R * temperatures)  # Convert ΔG to J/mol for consistency
    Ka_values = np.exp(ln_Ka_values)

    # Create a results DataFrame
    results = pd.DataFrame({
        'Temperature (K)': temperatures,
        'ΔG° (kJ/mol)': dG_values,
        'ln(K_a)': ln_Ka_values,
        'K_a': Ka_values,
    })

    if K_y and Convertion == 0:
        Convertion = []
        for K_a_i in results["K_a"]:
            if K_a_i >  1e5:
                Convertion.append(1)
            else:
                equation = Eq(K_y, K_a_i)
                solutions = solve(equation, x)
                solutions = [sol for sol in solutions if sol.is_real and sol > 0]
                Convertion.append(solutions[0])
        results['X'] = Convertion
    elif K_y and Convertion > 0:
        K_y = K_y.subs(x,0.8).evalf()
        print(f"K_y = {K_y} for a convertion of {Convertion}")
        Tx = (dH_f0*1000) / (dS0*1000 - R*math.log(K_y))
        print(f"For a conversion {Convertion} a temperature of {Tx} K is required")
    

    # Display results
    print(results)        
    return results

def graphs_T(results):
    """
    Generate plots for thermodynamic relationships.

    Parameters:
        results (DataFrame): DataFrame containing temperature, ΔG°, ln(K_a), and K_a values.
    """
    # Plot ln(K_a) vs 1/T
    plt.figure()
    plt.plot(1 / np.array(results["Temperature (K)"]), results["ln(K_a)"], marker='o', label='ln(K_a) vs 1/T')
    plt.xlabel('1/T (1/K)')
    plt.ylabel('ln(K_a)')
    plt.title('ln(K_a) vs 1/T')
    plt.grid()
    plt.show()

    # Plot ln(K_a) vs T
    plt.figure()
    plt.plot(results["Temperature (K)"], results["ln(K_a)"], marker='o', color='orange', label='ln(K_a) vs T')
    plt.xlabel('Temperature (K)')
    plt.ylabel('ln(K_a)')
    plt.title('ln(K_a) vs Temperature')
    plt.grid()
    plt.show()

    # Plot ΔG° vs T
    plt.figure()
    plt.plot(results["Temperature (K)"], results["ΔG° (kJ/mol)"], marker='o', color='orange', label='ΔG° vs T')
    plt.xlabel('Temperature (K)')
    plt.ylabel('ΔG° (kJ/mol)')
    plt.title('ΔG° vs Temperature')
    plt.grid()
    plt.show()

    if "X" in results.columns:
        # Plot X vs T
        plt.figure()
        plt.plot(results["Temperature (K)"], results["X"], marker='o', color='blue', label='X vs T')
        plt.xlabel('Temperature (K)')
        plt.ylabel('X')
        plt.title('X vs Temperature')
        plt.grid()
        plt.show()

def analysis_pressures(results, dG_f0, pressures, temperature, n):
    """
    Perform analysis of thermodynamic properties at varying pressures.
    """
    temperatures = np.full(n, temperature)
    results["Temperature (K)"] = temperatures
    results["ΔG° (kJ/mol)"] = dG_f0 + 8.314 / 1000 * temperatures * np.log(pressures)
    results["ln(K_a)"] = -np.array(results["ΔG° (kJ/mol)"]) * 1000 / (8.314 * temperatures)
    results["K_a"] = np.exp(results["ln(K_a)"])
    results["Pressure (atm)"] = pressures
    # Move the "Pressure (atm)" column to the second position
    cols = list(results.columns)  # Get all column names as a list
    cols.insert(1, cols.pop(cols.index("Pressure (atm)")))  # Insert the column at the second position
    results = results[cols]  # Reorganize the DataFrame with the new column order
    print(results)
    return results

def graphs_P(results):
    """
    Generate plots for thermodynamic relationships involving pressure (P).

    Parameters:
        results (DataFrame): DataFrame containing Pressure (atm), ΔG°, ln(K_a), and K_a values.
    """
    # Plot ln(K_a) vs Pressure
    plt.figure()
    plt.plot(1/ np.array(results["Pressure (atm)"]), results["ln(K_a)"], marker='o', label='ln(K_a) vs Pressure')
    plt.xlabel('1/P (1/atm)')
    plt.ylabel('ln(K_a)')
    plt.title('ln(K_a) vs Pressure')
    plt.grid()
    plt.show()

    # Plot ln(K_a) vs Pressure
    plt.figure()
    plt.plot(results["Pressure (atm)"], results["ln(K_a)"], marker='o', color='green', label='K_a vs Pressure')
    plt.xlabel('Pressure (atm)')
    plt.ylabel('ln(K_a)')
    plt.title('ln(K_a) vs Pressure')
    plt.grid()
    plt.show()

    # Plot ΔG° vs Pressure
    plt.figure()
    plt.plot(results["Pressure (atm)"], results["ΔG° (kJ/mol)"], marker='o', color='orange', label='ΔG° vs Pressure')
    plt.xlabel('Pressure (atm)')
    plt.ylabel('ΔG° (kJ/mol)')
    plt.title('ΔG° vs Pressure')
    plt.grid()
    plt.show()
