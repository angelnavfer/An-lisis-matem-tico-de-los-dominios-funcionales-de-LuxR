import numpy as np
import matplotlib.pyplot as plt

# Parámetros base
K1_base = 0.2  # μM^-1 (afinidad proteína-ligando)
K3_base = 1.0  # μM^-1 (afinidad proteína-ADN)
R2 = 100  # Cantidad de receptor dimerizado (μM)
L = 1  # Concentración fija de inductor (μM)
Phi = 10  # Un único valor de Phi

# Definir el rango de valores para K1 y K3
K1_values = np.logspace(-1, 3, 250)  # De 0,1 a 1000 μM^-1
K3_values = np.logspace(-1, 3, 250)  # De 0,1 a 1000 μM^-1

# Función para calcular G*
def G_star(L, K1, K3, R2, Phi):
    Kapp = (1 / K1) * ((1 / (K3 * (R2 / Phi))) / (Phi + (1 / (K3 * (R2 / Phi)))))  # Afinidad efectiva
    Gm = 100 * (Phi / (Phi + (1 / (K3 * (R2 / Phi)))))  # Amplitud máxima
    return Gm * (L / (Kapp + L))  # Expresión de G*

# Calcular G* para cada valor de K1 y K3 con el Phi fijo
G1_values = np.array([G_star(L, K1, K3_base, R2, Phi) for K1 in K1_values])
G3_values = np.array([G_star(L, K1_base, K3, R2, Phi) for K3 in K3_values])

# Calcular la primera derivada usando diferencias finitas
dG1_dK1 = np.gradient(G1_values, K1_values)  # dG*/dK1
dG3_dK3 = np.gradient(G3_values, K3_values)  # dG*/dK3

# Calcular la elasticidad para K1 y K3
elasticidad_K1 = (dG1_dK1 * K1_values) / G1_values
elasticidad_K3 = (dG3_dK3 * K3_values) / G3_values

# Crear figura con curvas solapadas de sensibilidad
plt.figure(figsize=(10, 6))

plt.semilogx(K1_values, dG1_dK1, label="d[G*]/dK1 (Sensibilidad)", color="blue", linestyle="dashed")
plt.semilogx(K3_values, dG3_dK3, label="d[G*]/dK3 (Sensibilidad)", color="red", linestyle="solid")

plt.xlabel("Valor de K1 o K3 [μM^-1]")
plt.ylabel("Derivada de [G*]")
plt.title(f"Sensibilidad de [G*] a K1 y K3")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)

plt.show()

# Crear figura para la elasticidad
plt.figure(figsize=(10, 6))

plt.semilogx(K1_values, elasticidad_K1, label="Elasticidad respecto a K1", color="orange", linestyle="dashed")
plt.semilogx(K3_values, elasticidad_K3, label="Elasticidad respecto a K3", color="green", linestyle="solid")

plt.xlabel("Valor de K1 o K3 [μM^-1]")
plt.ylabel("Elasticidad de [G*]")
plt.title("Elasticidad de [G*] respecto a K1 y K3")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)

plt.show()