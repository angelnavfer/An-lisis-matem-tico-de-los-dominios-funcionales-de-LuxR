import numpy as np
import matplotlib.pyplot as plt

# Parámetros base
K1_base = 0.2  # μM^-1 (afinidad proteína-ligando)
K3_base = 1.0  # μM^-1 (afinidad proteína-ADN)
R2 = 100       # Cantidad de receptor dimerizado (μM)

# Rango de concentraciones de inductor [L]
L_values = np.logspace(-4, 2, 250)  # Desde 0.001 hasta 100 μM

# Rango de valores para Phi
Phi_values = [0.1, 1, 10]  # Diferentes intensidades del RBS

# Rango de valores para K1 y K3
K1_values = np.logspace(-1, 3, 5)  # Variación de K1 entre 1 y 10^3 μM^-1
K3_values = np.logspace(-1, 3, 5)  # Variación de K3 entre 1 y 10^3 μM^-1

# Función para calcular G*
def G_star(L, K1, K3, R2, Phi):
    Kapp = (1 / K1) * ((1 / (K3 * (R2 / Phi))) / (Phi + (1 / (K3 * (R2 / Phi)))))  # Afinidad efectiva
    Gm = 100 * (Phi / (Phi + (1 / (K3 * (R2 / Phi)))))  # Amplitud máxima de la señal
    return Gm * (L / (Kapp + L))  # Expresión de G*

for Phi in Phi_values:
    for K1 in K1_values:
        G_values = [G_star(L, K1, K3_base, R2, Phi) for L in L_values]
        print(f"K1={K1}, Phi={Phi}, G_values[:5]={G_values[:5]}")

for Phi in Phi_values:
    for K3 in K3_values:
        G_values = [G_star(L, K1_base, K3, R2, Phi) for L in L_values]
        print(f"K3={K3}, Phi={Phi}, G_values[:5]={G_values[:5]}")

# Gráficos para variaciones en K1 con diferentes Phi
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Gráfico 1: Variando K1
for Phi in Phi_values:
    for K1 in K1_values:
        G_values = [G_star(L, K1, K3_base, R2, Phi) for L in L_values]
        axes[0].semilogx(L_values, G_values, label=f"K1={K1}, Phi={Phi}")

axes[0].set_xlabel("Concentración de inductor [L] [μM]")
axes[0].set_ylabel("[G*] [μM]")
axes[0].set_title("Variación de K1 en G* para distintos Phi")
axes[0].legend()
axes[0].legend(loc="upper left", bbox_to_anchor=(1.05, 1))
axes[1].legend(loc="upper left", bbox_to_anchor=(1.05, 1))
axes[0].grid(True, which="both", linestyle="--", linewidth=0.5)

# Gráfico 2: Variando K3
for Phi in Phi_values:
    for K3 in K3_values:
        G_values = [G_star(L, K1_base, K3, R2, Phi) for L in L_values]
        axes[1].semilogx(L_values, G_values, label=f"K3={K3}, Phi={Phi}")

axes[1].set_xlabel("Concentración de inductor [L] [μM]")
axes[1].set_ylabel("[G*] [μM]")
axes[1].set_title("Variación de K3 en G* para distintos Phi")
axes[1].legend()
axes[0].legend(loc="upper left", bbox_to_anchor=(1.05, 1))
axes[1].legend(loc="upper left", bbox_to_anchor=(1.05, 1))
axes[1].grid(True, which="both", linestyle="--", linewidth=0.5)

# Ajustar disposición y mostrar gráficos
plt.tight_layout()
plt.show()





