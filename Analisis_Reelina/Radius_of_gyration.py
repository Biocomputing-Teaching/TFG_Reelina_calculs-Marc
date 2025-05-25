import sys
sys.path.append("/home/10033944/paquets")
import MDAnalysis as mda
import pandas as pd
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Load your trajectory and topology files
import os
INPUT_NC="Sortides_MD_reelina"
INPUT_PARM7="seq_parm7"

fitxers_nc = [f for f in os.listdir(INPUT_NC) if f.endswith(".nc")]

# Generar els noms dels fitxers .parm7 corresponents
fitxers_parm7 = [f.replace(".nc", ".parm7") for f in fitxers_nc]

# Iterar sobre les parelles de fitxers .nc i .parm7
for conf_nc, conf_parm7 in zip(fitxers_nc, fitxers_parm7):
    u = mda.Universe(
        os.path.expanduser(f"{INPUT_PARM7}/{conf_parm7}"),
        os.path.expanduser(f"{INPUT_NC}/{conf_nc}"),
        topology_format="PARM7"
    )

    # Seleccionar els àtoms de la proteïna
    protein = u.select_atoms("protein")

    # Calcular el radi de giració
    time = []
    rgyr = []

    for frame in u.trajectory:
        time.append(u.trajectory.time)
        rgyr.append(protein.radius_of_gyration())

    # Crear un CSV amb els resultats
    rgyr_df = pd.DataFrame(list(zip(time, rgyr)),
                           columns=['Time (ps)', 'Radius of gyration (A)'],
                           index=None)
    rgyr_df.to_csv(f"RG_{conf_nc.replace('.nc', '')}.csv", index=False)


import os
import pandas as pd
import matplotlib.pyplot as plt

# Llista tots els fitxers CSV creats (que comencen amb "RG_")
csv_files = [f for f in os.listdir() if f.startswith("RG_") and f.endswith(".csv")]

plt.figure(figsize=(12, 8))
for csv_file in csv_files:
    # Carrega el fitxer CSV
    df = pd.read_csv(csv_file)
    # Ploteja la corba (Temps vs. Radi de giració)
    plt.plot(df["Time (ps)"], df["Radius of gyration (A)"], label=csv_file)

plt.xlabel("Time (ps)")
plt.ylabel("Radius of gyration (A)")
plt.title("Superposició del radi de giració per a cada simulació")
plt.legend(loc="best")
plt.show()
