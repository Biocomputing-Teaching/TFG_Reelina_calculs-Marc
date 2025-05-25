#El RMSF (Root Mean Square Fluctuation) és una mesura de la 
#flexibilitat d'una proteïna, que calcula les fluctuacions 
#mitjanes de les posicions dels àtoms al llarg del temps. 
#Això et permet identificar quines parts de la proteïna són 
#més flexibles o rígides durant una simulació de dinàmica molecular.

import sys
sys.path.append("/home/10033944/paquets")
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSF
import matplotlib.pyplot as plt
import warnings
import pandas as pd
warnings.filterwarnings('ignore')

# Load your trajectory and topology files
import os

# Directori dels fitxers de trajectòria (.nc) i topologia (.parm7)
INPUT_NC = "Sortides_MD_reelina"
INPUT_PARM7 = "seq_parm7"

# Llista de fitxers .nc
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
    selection = u.select_atoms("protein")

    # Calcular RMSF
    rmsf = RMSF(selection).run()

    # Convertir les dades a un DataFrame de pandas
    rmsf_df = pd.DataFrame({
        "Residue ID": selection.resids,
        "RMSF (A)": rmsf.rmsf
    })

    # Guardar el DataFrame en un CSV amb el nom del fitxer original
    rmsf_df.to_csv(f"RMSF_{conf_nc.replace('.nc', '')}.csv", index=False)

# Després de finalitzar el bucle que processa cada trajectòria,
# llegim els CSV generats i plotejem els RMSF per cada simulació

plt.figure(figsize=(12, 8))
for fitxer in fitxers_nc:
    # Generem el nom del CSV a partir del nom del fitxer .nc
    nom_csv = f"RMSF_{fitxer.replace('.nc', '')}.csv"
    # Carreguem el CSV en un DataFrame
    df = pd.read_csv(nom_csv)
    # Plotejem la línia corresponent: l'eix X és el "Residue ID" i l'eix Y el "RMSF (A)"
    plt.plot(df["Residue ID"], df["RMSF (A)"], label=f"Simulació {fitxer}")

plt.xlabel("Residue ID")
plt.ylabel("RMSF (Å)")
plt.title("Comparació del RMSF per cada simulació")
plt.legend(loc="best")
plt.show()
