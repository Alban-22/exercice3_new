import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import re


a=20

plt.rcParams.update({
    'figure.figsize': (7,5),      # taille des figures (largeur, hauteur)
    'font.size': a,               # taille globale de police
    'axes.labelsize': a,          # taille des labels axes
    'axes.titlesize': a,          # taille des titres
    'axes.titleweight': 'bold',    # titre en gras
    'legend.fontsize': a,         # taille légende
    'legend.title_fontsize': a,   # titre légende si tu en as
    'xtick.labelsize': a,         # taille ticks x
    'ytick.labelsize': a,         # taille ticks y
    'xtick.direction': 'in',       # ticks vers l'intérieur
    'ytick.direction': 'in',
    'lines.linewidth': 2,          # épaisseur des courbes
    'lines.markersize': 8,         # taille des points
})

# ============================================================
# USER SETTINGS
# ============================================================

folder = rf"C:/EPFL/Semestre_4/Physique_numérique/exercice3_new/Scan_epsilon_Gravit_s_0.9_epsilon_1e+03/"

# ============================================================
# Output folder
# ============================================================


fig_dir = os.path.join(folder, "figures")
os.makedirs(fig_dir, exist_ok=True)

# ============================================================
# Scan files
# ============================================================

files = sorted(glob.glob(os.path.join(folder, "*.txt")))

files = [
    f for f in files
    if re.match(r".*_[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?\.txt$", os.path.basename(f))
]

datasets = []
epsilon_values = []

for f in files:
    name = os.path.basename(f)[:-4]
    parts = name.split("_")

    param_name = parts[-2]
    value = float(parts[-1])

    if param_name != "epsilon":
        continue

    data = np.loadtxt(f)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    datasets.append(data)
    epsilon_values.append(value)

print(f"Found {len(datasets)} datasets.")

if len(datasets) == 0:
    raise ValueError("Aucun dataset trouvé. Vérifie le dossier ou le nom des fichiers.")

# ============================================================
# Accès données
# ============================================================

#Pour n=0 et N=2.

i = 0
data = datasets[i]

t = data[:, 0]

# particule 0
x0  = data[:, 1]
y0  = data[:, 2]
vx0 = data[:, 3]
vy0 = data[:, 4]

# particule 1
x1  = data[:, 5]
y1  = data[:, 6]
vx1 = data[:, 7]
vy1 = data[:, 8]

# énergies globales totale.
Ec = data[:, -5]
Ep = data[:, -4]
Em = data[:, -3]
P  = data[:, -2]
dt = data[:, -1]

# ============================================================
#Energies totales
# ============================================================

plt.figure()

plt.plot(t, Em, color="black", label=rf"$E_m \,\,[\mathrm{{J}}]$")
plt.plot(t, Ec, color="blue", label=rf"$E_c \,\,[\mathrm{{J}}]$")
plt.plot(t, Ep, color="red", label=rf"$E_p \,\,[\mathrm{{J}}]$")

plt.xlabel(rf"$t \,\,[\mathrm{{s}}]$")
plt.ylabel(rf"$E \,\,[\mathrm{{J}}]$")
#plt.title(rf"Evolution des énergies du système Terre–Lune")

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(os.path.join(fig_dir, "energies_evolution.png"), dpi=300)


# ============================================================
#distance Terre-Lune
# ============================================================
d= 384748*1000 # distance terre-lune en m

plt.figure()
distance = np.sqrt((x1 - x0)**2 + (y1 - y0)**2)
plt.plot(t, distance-d, color="purple", label=rf"Distance Terre–Lune")
plt.hlines(0, t[0], t[-1], color="black", linestyle="--", label=rf"Distance cible $d$")

plt.xlabel(rf"$t \,\,[\mathrm{{s}}]$")
plt.ylabel(rf"$||\vec{{r_T}} - \vec{{r_L}}|| - d \,\,[\mathrm{{m}}]$")

plt.grid(True, alpha=0.3)   
plt.legend()
plt.tight_layout()

plt.savefig(os.path.join(fig_dir, "distance_terre_lune.png"), dpi=300)

# ============================================================
#distance Terre-Lune 2
# ============================================================
d= 384748*1000 # distance terre-lune en m

plt.figure()
distance = np.sqrt((x1 - x0)**2 + (y1 - y0)**2)
plt.plot(t, np.abs((distance-d)/d), color="purple", label=rf"Distance Terre–Lune")
plt.hlines(0, t[0], t[-1], color="black", linestyle="--", label=rf"Distance cible $d$")

plt.xlabel(rf"$t \,\,[\mathrm{{s}}]$")
plt.ylabel(rf"$\frac{{| ||\vec{{r_T}} - \vec{{r_L}}|| - d|}}{{d}} $")

plt.grid(True, alpha=0.3)   
plt.legend()
plt.tight_layout()

plt.savefig(os.path.join(fig_dir, "distance_terre_lune2.png"), dpi=300)


# ============================================================
# conservation de la quantité de mouvement
# ============================================================

plt.figure()
plt.plot(t, P, color="red", label=rf"Quantité de mouvement $P$")

plt.xlabel(rf"$t \,\,[\mathrm{{s}}]$")
plt.ylabel(rf"$P \,\,[\mathrm{{kg\cdot m\cdot s^{{-1}}}}]$")

plt.grid(True, alpha=0.3)   
plt.legend()
plt.tight_layout()

plt.savefig(os.path.join(fig_dir, "quantite_mouvement.png"), dpi=300)