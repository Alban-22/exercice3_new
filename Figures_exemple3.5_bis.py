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


colors = [
    '#F4D03F',  # sonde (jaune doux)
    '#1B4F72',  # Terre (bleu nuit)
    '#D5D8DC'   # Lune (gris lunaire)
]

colors = [
    '#FFD700',  # Sonde → or (visible, attire l'œil)
    '#2E86C1',  # Terre → bleu profond
    '#BFC9CA'   # Lune → gris clair
]

plt.rcParams.update({
    'figure.figsize': (7,5),

    # Texte
    'font.size': a,
    'axes.labelsize': a,
    'axes.titlesize': a+2,
    'axes.titleweight': 'bold',

    # Axes
    'axes.edgecolor': '#333333',
    'axes.linewidth': 1.2,

    # Grille
    'axes.grid': True,
    'grid.color': "#0B0C60",
    'grid.linestyle': '--',
    'grid.alpha': 0.5,

    # Ticks
    'xtick.labelsize': a-1,
    'ytick.labelsize': a-1,
    'xtick.direction': 'in',
    'ytick.direction': 'in',

    # Courbes
    'lines.linewidth': 2.2,

    # Fond
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
})


#===========================================================================
#option espace
#==========================================================================

# ============================================================
# USER SETTINGS MODIFER A CHAQUE NOUVELLE SIMULATION
# ============================================================

folder = rf"C:/EPFL/Semestre_4/Physique_numérique/exercice3_new/Scan_alphadeg0_Gravit_s_0.9_epsilon_1e-06/"

 # une taille par particule
R_T= 6378.1e3
R_L =1737.4e3
R_s= 2.51

R = [R_s, R_T, R_L]

scale = 300 

point_sizes = [(r / R_T)**2 * scale for r in R]

parameterscan=  variable_array =np.linspace(180 - 2/60,180+2/60, 5)

#============================================================

#rf"C:/EPFL/Semestre_4/Physique_numérique/exercice3_new/Scan_alphadeg0_Gravit_s_0.9_epsilon_1e-08/"

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
#colors = plt.cm.tab10(range(nb_particules))


for i in range(nb_simu):
    row = i // nb_cols
    col = i % nb_cols

    ax = axes[row][col] if nb_rows > 1 else axes[col]

    
    for j in range(len(x[i])):
        ax.plot(x[i][j], y[i][j], color=colors[j], alpha=0.9)
        
        ax.scatter(x[i][j][-1], y[i][j][-1],
           color=colors[j],
           s=point_sizes[j],
           edgecolors='black',
           linewidths=0.6,
           zorder=3)
        
        x_terre = x[i][1][-1]
        y_terre = y[i][1][-1]

        ax.scatter(x_terre, y_terre,
               color='red',
               marker='+',
               s=scale/10,linewidths=2,
           zorder=4)
           
        
            
    ax.set_title(rf"$\alpha_0={parameterscan[i]:.3f}^\circ$")
    ax.set_xlabel(rf"$x \,\,[\mathrm{{m}}]$")
    ax.set_ylabel(rf"$y \,\,[\mathrm{{m}}]$")
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.5)

# Supprimer les cases vides
for i in range(nb_simu, nb_rows * nb_cols):
    row = i // nb_cols
    col = i % nb_cols
    fig.delaxes(axes[row][col] if nb_rows > 1 else axes[col])


plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "trajectoires.png"), dpi=300)   
plt.close() 

# ============================================================
# Trajectoires une simulation:  i
# ============================================================

def plot_simulation(i):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(6,6))

    nb_particules = len(x[i])
    colors = plt.cm.tab10(range(nb_particules))


    for j in range(nb_particules):
        plt.plot(x[i][j], y[i][j], color=colors[j])
        
        plt.scatter(x[i][j][-1], y[i][j][-1],
               color=colors[j],
               s=point_sizes[j])
    


    plt.xlabel(rf"$x \,\,[\mathrm{{m}}]$")
    plt.ylabel(rf"$y \,\,[\mathrm{{m}}]$")
    plt.title(rf"$\alpha_0={parameterscan[i]:.3f}^\circ$")

    plt.grid(True, alpha=0.5)
    plt.gca().set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"trajectoire_simulation_{i}.png"), dpi=300)    
    plt.close()

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

print("Hauteurs minimales:", h_min)

plt.plot(parameterscan, h_min, 'o-', color='blue')
plt.xlabel(rf"$\alpha_0 \,\,[\mathrm{{°}}]$")
plt.ylabel(rf"$h_{{min}} \,\,[\mathrm{{m}}]$")

plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "hauteur_min.png"), dpi=300)  
plt.close()


# ============================================================
#acceleration
# ============================================================

acc_max = np.max(acc[:, 0, :], axis=1)
plt.plot(parameterscan, acc_max, 'o-', color= 'red')
plt.xlabel(rf"$\alpha_0 \,\,[\mathrm{{°}}]$")
plt.ylabel(rf"$a_{{max}} \,\,[\mathrm{{m \cdot s^{{-2}} }}]$")
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "acceleration_max.png"), dpi=300)
plt.close()
