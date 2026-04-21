import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import re
import math


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

folder = rf"C:/EPFL/Semestre_4/Physique_numérique/exercice3_new/Scan_alphadeg0_Gravit_s_0.9_epsilon_1e-08/"

# ============================================================
# Output folder
# ============================================================

fig_dir = os.path.join(folder, "figures")
os.makedirs(fig_dir, exist_ok=True)

# Récupérer tous les fichiers .txt
files = sorted(glob.glob(os.path.join(folder, "*.txt")))



datasets = []
parameterscan = [] # On commence avec une liste VIDE
parameterscan_name = "" 

for f in files:
    name = os.path.basename(f)[:-4]
    parts = name.split("_")

    param_name = parts[-2]
    value = float(parts[-1])

    # On vérifie qu'on lit bien le bon paramètre
    if "alphadeg" not in param_name: 
        continue
    
    parameterscan_name = param_name 

    data = np.loadtxt(f)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    datasets.append(data)
    parameterscan.append(value) # Ici, on remplit dynamiquement

# Une fois la boucle finie, on peut trier les deux listes ensemble 
# pour que les angles soient dans l'ordre croissant sur les graphes
indices = np.argsort(parameterscan)
parameterscan = [parameterscan[i] for i in indices]
datasets = [datasets[i] for i in indices]

print(f"Found {len(datasets)} datasets for parameter: {parameterscan_name}")
# ============================================================
# Accès données
# ============================================================

n=1
N=2

t=[]
x=[]
y=[]
vx=[]
vy=[]
acc=[]
Ec=[]
Ep=[]
Em=[]
P=[]
dt=[]


for i in range(len(datasets)):
    
    data = datasets[i]
    #les variables t_i x_i etc de chaque simulation coorspondent a la ieme simulation,
    # ce sont des vecteurs contenant toutes les valeurs de la simulations
    t_ = data[:, 0]
    #sur cette boucle, on range les x suivant le numéro: x--> x0, x1,.xj.. 
    # ou c est des coordonnees de xj pendant la simulatioon i. 
    x_ = []
    y_ = []  
    vx_ = []
    vy_ = []
    acc_ = []
    #x[j] donne la simulation de xj pour la i eme simulation
    for j in range(n+N):
        x_.append(data[:, 1+5*j])
        y_.append(data[:, 2+5*j])
        vx_.append(data[:, 3+5*j])
        vy_.append(data[:, 4+5*j])
        acc_.append(data[:, 5+5*j])
    Ec_ = data[:, -5]
    Ep_ = data[:, -4]
    Em_ = data[:, -3]
    P_  = data[:, -2]
    dt_ = data[:, -1]
    
    t.append(t_)
    x.append(x_)
    y.append(y_)
    vx.append(vx_)
    vy.append(vy_)
    acc.append(acc_)
    Ec.append(Ec_)
    Ep.append(Ep_)
    Em.append(Em_)
    P.append(P_)
    dt.append(dt_)
    
t=np.array(t)
x=np.array(x)
y=np.array(y)
vx=np.array(vx)
vy=np.array(vy)
acc=np.array(acc)
Ec=np.array(Ec)
Ep=np.array(Ep)
Em=np.array(Em)
P=np.array(P)
dt=np.array(dt)


# ============================================================
# Bilan : t[i] donne le vecteur temps de la i eme simulation,
#
# acc[i][j] donne le vecteur des acc de la j eme particule de la i eme simulation, etc
# vy[i][j] donne le vecteur des vy de la j eme particule de la i eme simulation, etc
# vx[i][j] donne le vecteur des vx de la j eme particule de la i eme simulation, etc
# y[i][j] donne le vecteur des y de la j eme particule de la i eme simulation, etc
# x[i][j] donne le vecteur des x de la j eme particule de la i eme simulation, etc
#acc[i][j] donne le vecteur des acc de la j eme particule de la i eme simulation, etc
# Ec[i] donne le vecteur de l'energie cinetique de la i eme simulation, etc
# Ep[i] donne le vecteur de l'energie potentielle de la i eme simulation, etc
# Em[i] donne le vecteur de l'energie mecanique de la i eme simulation, etc
# P[i] donne le vecteur de la quantite de mouvement de la i eme simulation, etc
# dt[i] donne le vecteur des pas de temps de la i eme simulation, etc
# ============================================================

#============================================================
R_T= 6378.1e3
parameterscan=  np.linspace(120, 240, 6)
#============================================================



# ============================================================
# Trajectoires toutes les simulations
# ============================================================

nb_simu = len(x)
nb_cols = 3
nb_rows = math.ceil(nb_simu / nb_cols)

fig, axes = plt.subplots(nb_rows, nb_cols, figsize=(15, 5 * nb_rows))

if nb_rows == 1:
    axes = [axes]

nb_particules = len(x[0])

# Couleurs fixes pour chaque particule
colors = plt.cm.tab10(range(nb_particules))

xmin = min(min(min(x[i][j]) for j in range(len(x[i]))) for i in range(nb_simu))
xmax = max(max(max(x[i][j]) for j in range(len(x[i]))) for i in range(nb_simu))
ymin = min(min(min(y[i][j]) for j in range(len(y[i]))) for i in range(nb_simu))
ymax = max(max(max(y[i][j]) for j in range(len(y[i]))) for i in range(nb_simu))

for i in range(nb_simu):
    row = i // nb_cols
    col = i % nb_cols

    ax = axes[row][col] if nb_rows > 1 else axes[col]

    for j in range(len(x[i])):
        ax.plot(x[i][j], y[i][j], color=colors[j])

    ax.set_title(rf"$\alpha_0={parameterscan[i]:.3f}^\circ$")
    ax.set_xlabel(rf"$x \,\,[\mathrm{{m}}]$")
    ax.set_ylabel(rf"$y \,\,[\mathrm{{m}}]$")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.grid(True)

# Supprimer les cases vides
for i in range(nb_simu, nb_rows * nb_cols):
    row = i // nb_cols
    col = i % nb_cols
    fig.delaxes(axes[row][col] if nb_rows > 1 else axes[col])

fig.suptitle("Trajectoires des particules pour chaque simulation", fontsize=16)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "trajectoires.png"), dpi=300)    

# ============================================================
# Trajectoires une simulation:  i
# ============================================================

def plot_simulation(i):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(6,6))

    nb_particules = len(x[i])
    colors = plt.cm.tab10(range(nb_particules))

    xmin = min(min(x[i][j]) for j in range(nb_particules))
    xmax = max(max(x[i][j]) for j in range(nb_particules))
    ymin = min(min(y[i][j]) for j in range(nb_particules))
    ymax = max(max(y[i][j]) for j in range(nb_particules))

    for j in range(nb_particules):
        plt.plot(x[i][j], y[i][j], color=colors[j])

    plt.xlabel(rf"$x \,\,[\mathrm{{m}}]$")
    plt.ylabel(rf"$y \,\,[\mathrm{{m}}]$")
    plt.title(rf"$\alpha_0={parameterscan[i]:.3f}^\circ$")

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.grid(True)
    plt.gca().set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"trajectoire_simulation_{i}.png"), dpi=300)    

plot_simulation(0)

# ============================================================
# hauteur 
# ============================================================

r = []

for i in range(len(x)):

    r0 = np.array([x[i][0], y[i][0]])   # (2, N)
    r1 = np.array([x[i][1], y[i][1]])   # (2, N)

    ri = np.linalg.norm(r0 - r1, axis=0)  # norme pour chaque k

    r.append(ri)

r = np.array(r)
h = r - R_T
h_min = np.min(h, axis=1)

plt.plot(parameterscan, h_min, 'o-')

plt.xlabel(rf"$\alpha_0 \,\,[\mathrm{{°}}]$")
plt.ylabel(rf"$h_{{min}} \,\,[\mathrm{{m}}]$")

plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "hauteur_min.png"), dpi=300)  


# ============================================================
#acceleration
# ============================================================

acc_max = np.max(acc[:, 0, :], axis=1)
plt.plot(parameterscan, acc_max, 'o-')
plt.xlabel(rf"$\alpha_0 \,\,[\mathrm{{°}}]$")
plt.ylabel(rf"$a_{{max}} \,\,[\mathrm{{m \cdot s^{{-2}} }}]$")
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "acceleration_max.png"), dpi=300)

