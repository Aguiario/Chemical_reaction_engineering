import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import sympy as sp
import cantera as ct


# Constantes
R = 8.314  # Constante universal de gases, J/(mol·K)
T0 = 298.15  # Temperatura estándar, K
P = 101325  # Presión estándar (1 atm en Pa)
df = pd.read_excel('SP_298K.xlsx')

def standar_propieties(reactants, products):
    # Crear copias explícitas
    reactants_properties_298K = df[df['Species'].isin(reactants[:, 0])].copy()
    products_properties_298K = df[df['Species'].isin(products[:, 0])].copy()

    # Añadir columnas para coeficientes y propiedades ponderadas
    reactants_properties_298K["Coefficient"] = reactants[:, 1].astype(float)
    products_properties_298K["Coefficient"] = products[:, 1].astype(float)

    reactants_properties_298K["Weighted_DHf"] = reactants_properties_298K["DHf°[kJ/mol]"] * reactants_properties_298K["Coefficient"]
    reactants_properties_298K["Weighted_DGf"] = reactants_properties_298K["DGf°[kJ/mol]"] * reactants_properties_298K["Coefficient"]
    reactants_properties_298K["Weighted_S"] = reactants_properties_298K["S°[J/K·mol]"] * reactants_properties_298K["Coefficient"]

    products_properties_298K["Weighted_DHf"] = products_properties_298K["DHf°[kJ/mol]"] * products_properties_298K["Coefficient"]
    products_properties_298K["Weighted_DGf"] = products_properties_298K["DGf°[kJ/mol]"] * products_properties_298K["Coefficient"]
    products_properties_298K["Weighted_S"] = products_properties_298K["S°[J/K·mol]"] * products_properties_298K["Coefficient"]

    # Suma ponderada de productos y reactivos
    dH_f0 = products_properties_298K["Weighted_DHf"].sum() - reactants_properties_298K["Weighted_DHf"].sum()
    dG_f0 = products_properties_298K["Weighted_DGf"].sum() - reactants_properties_298K["Weighted_DGf"].sum()
    # Si se hace con la informacion de la tabla
    #dS0 = products_properties_298K["Weighted_S"].sum() - reactants_properties_298K["Weighted_S"].sum()
    # Suponiendo que ΔHº y ΔSº son independientes de temperatura se calcula ΔSº con los valores de la última tabla (T=298,15K) a partir de la relación ΔGº = ΔHº- T ΔSº
    dS0 = (dG_f0 - dH_f0)/T0

    # Imprime las filas filtradas

    print("Reactants Properties")
    print(reactants_properties_298K)
    print("\n")
    print("Products Properties")
    print(products_properties_298K)
    print("\n")

    print(f"ΔGº: {dG_f0} [kJ/mol], y {'es' if dG_f0 > 0 else 'no es'} espontáneo.")
    print(f"ΔHº: {dH_f0} [kJ/mol], y {'es endotérmico.' if dH_f0 > 0 else 'es exotérmico.'}")
    print(f"Sº: {dS0} [J/K·mol], y el desorden {'aumenta.' if dH_f0 > 0 else 'disminuye.'}")

    return dG_f0, dH_f0, dS0, reactants_properties_298K, products_properties_298K