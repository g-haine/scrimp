import numpy as np
import matplotlib.pyplot as plt

# Définition de la fonction u(t)
def u(t):
    max_val = np.maximum(90, t)
    num = (max_val - 100)**2
    denom = (500 - max_val)**2
    exp_part = np.exp(-num / denom)
    return 0.2 * np.minimum(1, exp_part)

# Génération des valeurs de t et calcul de u(t)
t_values = np.linspace(0, 500, 1000)
u_values = u(t_values)

# Tracé de la courbe
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 32
fig=plt.figure(figsize=(20, 2), layout="constrained")
plt.plot(t_values, u_values, label=r'$u_{\bf t}(t)$', color='blue')
plt.title('Time evolution of the tangential control')
plt.xlabel(r'Time $t$ (s)')
plt.ylabel(r'$u_{\bf t}(t)$')
plt.grid(True)
fig.legend([r'$u_{\bf t}(t)$'],loc='outside center right')
plt.show()
