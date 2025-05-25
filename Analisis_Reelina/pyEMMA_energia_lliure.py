#!/home/marct/miniforge3/envs/pyemma_env/bin/python


import pyemma
import matplotlib.pyplot as plt
import numpy as np
import glob

# Llista de trajectòries a analitzar
trajectories = glob.glob("sequencia_0*.nc")

# Emmagatzemar energies lliures
free_energies = []

for traj in trajectories:
    # Càrrega de dades (ajusta segons les teves variables)
    data = np.loadtxt(traj)  # Exemple, substitueix amb la teva extracció de dades
    tica=pyemma.coordinates.tica(data, lag=10)
    reduced_data = tica.get_output()
    
    # Estimar energia lliure
    fe = pyemma.plots.plot_free_energy(*reduced_data.T)
    free_energies.append(fe)

# Comparació de les energies lliures de diverses simulacions
plt.figure(figsize=(8, 6))

for i, fe in enumerate(free_energies):
    plt.plot(fe[0], fe[1], label=f"Simulació {i+1}")

plt.xlabel("Coordenada de reacció")
plt.ylabel("Energia lliure (kcal/mol)")
plt.title("Comparació de les energies lliures entre simulacions")
plt.legend()
plt.show()
