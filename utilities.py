import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import sympy as sp
import cantera as ct
from sympy import symbols, Eq, solve

# Constants
R = 8.314  # Universal gas constant in J/(mol·K)
T0 = 298.15  # Standard temperature in K
P = 1  # Standard atmospheric pressure in atm
df = pd.read_excel('SP_298K.xlsx')  # Load thermodynamic properties from an Excel file

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
    reactants_dict = {reactant[0]: float(reactant[1]) for reactant in reactants}
    products_dict = {product[0]: float(product[1]) for product in products}

    # Filter properties for reactants and products based on the species names
    reactants_properties_298K = df[df['Species'].isin(reactants_dict.keys())].copy()
    products_properties_298K = df[df['Species'].isin(products_dict.keys())].copy()

    # Assign stoichiometric coefficients to the respective species
    reactants_properties_298K["Coefficient"] = reactants_properties_298K["Species"].map(reactants_dict)
    products_properties_298K["Coefficient"] = products_properties_298K["Species"].map(products_dict)

    # Compute the total stoichiometric coefficient
    vi = products_properties_298K["Coefficient"].sum() - reactants_properties_298K["Coefficient"].sum()

    # Calculate weighted properties for reactants
    reactants_properties_298K["Weighted_DHf"] = reactants_properties_298K["DHf°[kJ/mol]"] * reactants_properties_298K["Coefficient"]
    reactants_properties_298K["Weighted_DGf"] = reactants_properties_298K["DGf°[kJ/mol]"] * reactants_properties_298K["Coefficient"]
    reactants_properties_298K["Weighted_S"] = reactants_properties_298K["S°[J/K·mol]"] * reactants_properties_298K["Coefficient"]

    # Calculate weighted properties for products
    products_properties_298K["Weighted_DHf"] = products_properties_298K["DHf°[kJ/mol]"] * products_properties_298K["Coefficient"]
    products_properties_298K["Weighted_DGf"] = products_properties_298K["DGf°[kJ/mol]"] * products_properties_298K["Coefficient"]
    products_properties_298K["Weighted_S"] = products_properties_298K["S°[J/K·mol]"] * products_properties_298K["Coefficient"]

    # Compute overall property changes (ΔH°, ΔG°, ΔS°)
    dH_f0 = products_properties_298K["Weighted_DHf"].sum() - reactants_properties_298K["Weighted_DHf"].sum()
    dG_f0 = products_properties_298K["Weighted_DGf"].sum() - reactants_properties_298K["Weighted_DGf"].sum()

    # Calculate ΔS° using the relationship ΔG° = ΔH° - TΔS°
    dS0 = (dH_f0 - dG_f0) / T0

    ln_Ka_0 = -dG_f0 * 1000 / (R * T0)  # Convert ΔG to J/mol for consistency
    Ka_0 = np.exp(ln_Ka_0)

    # Mass balance using conversion
    x = symbols('x')
    # Moles per component
    reactants_properties_298K["n [mol]"] = reactants_properties_298K["Coefficient"] * (1-x)
    products_properties_298K["n [mol]"] = products_properties_298K["Coefficient"] * x
    # Total number of moles
    n_total = reactants_properties_298K["n [mol]"].sum() + products_properties_298K["n [mol]"].sum()
    # Mass fraction of each component
    reactants_properties_298K["y_i"] = reactants_properties_298K["n [mol]"]/n_total
    products_properties_298K["y_i"] = products_properties_298K["n [mol]"]/n_total
    # Calculate Ky
    K_y = ((reactants_properties_298K["y_i"]**reactants_properties_298K["Coefficient"]).prod())/((products_properties_298K["y_i"]**products_properties_298K["Coefficient"]).prod()) 
    K_y = K_y *(P/P0)**vi

    # Display detailed reactant and product properties for debugging
    print("Reactant Properties:\n", reactants_properties_298K, "\n")
    print("Product Properties:\n", products_properties_298K, "\n")

    # Print thermodynamic results with descriptions
    print(f"ΔG°: {dG_f0} kJ/mol, {'non-spontaneous' if dG_f0 > 0 else 'spontaneous'}.")
    print(f"ΔH°: {dH_f0} kJ/mol, {'endothermic' if dH_f0 > 0 else 'exothermic'}.")
    print(f"ΔS°: {dS0} kJ/(K·mol), {'increased disorder' if dS0 > 0 else 'decreased disorder'}.")
    print(f"ln(K_a): {ln_Ka_0}")
    print(f"K_a: {Ka_0}")
    print(f"K_y: {K_y}")
    print(f"Net Stoichiometric Coefficient (v_i): {vi}")
    print(f"Total number of moles (n_t): {n_total}")

    return dG_f0, dH_f0, dS0, vi, K_y, reactants_properties_298K, products_properties_298K

def properties_temperature(temperatures, dH_f0, dS0, K_y = None, Convertion = 0.8):
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

    if K_y:
        x = symbols('x')
        if isinstance(K_y, Eq):
            Convertion = []
            for K_a_i in results["K_a"]:
                equation = Eq(K_y, K_a_i)
                solutions = solve(equation, x)
                Convertion.append(solutions[0])
        else:
            Convertion = K_y

    
        results['X'] = Convertion

    # Display results
    print(results)
    return results

def graphs(results):
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



