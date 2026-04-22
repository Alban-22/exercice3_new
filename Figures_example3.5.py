from matplotlib.patches import Circle
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
    "#190CA6",  # Sonde → or (visible, attire l'œil)
    "red",  
    'green'   # Lune → gris clair
]

u= 3
u_L=7
couleur = "#F98A02" 
couleur2 = "#6EF082"
# ============================================================
# USER SETTINGS MODIFER A CHAQUE NOUVELLE SIMULATION
# ============================================================

folder = rf"C:/EPFL/Semestre_4/Physique_numérique/exercice3_new/Scan_alphadeg0_ATM3_s_0.9_epsilon_1e-05/"
 
#folder = rf"C:/EPFL/Semestre_4/Physique_numérique/exercice3_new/Scan_alphadeg0_Gravit_s_0.9_epsilon_1e-15/"
 # une taille par particule
R_T= 6378.1e3
R_L =1737.4e3
R_s= 2.51

R = [R_s, R_T, R_L]

scale = 300 

point_sizes = [(r / R_T)**2 * scale for r in R]

parameterscan=variable_arrayvariable_array =np.linspace(165,175,10)

# 1: np.linspace(150,210,30)
# 2: np.linspace(168,193,30)
# 3: np.linspace(165,175,10)



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
    
t=np.array(t, dtype=object)
x=np.array(x, dtype=object)
y=np.array(y, dtype=object)
vx=np.array(vx, dtype=object)
vy=np.array(vy, dtype=object)
#acc=np.array(acc, dtype=object)
Ec=np.array(Ec, dtype=object)
Ep=np.array(Ep, dtype=object)
Em=np.array(Em, dtype=object)
P=np.array(P, dtype=object)
dt=np.array(dt, dtype=object)


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




# Couleurs fixes pour chaque particule
#colors = plt.cm.tab10(range(nb_particules))


# Nouvelle fonction pour dessiner un astre : point noir (taille réelle) + disque visible
def draw_body(ax, x, y, R_reel, color_visible, R_visible=None, color_point='black', zorder=2):
    # Disque visible (plus grand, couleur flashy)
    if R_visible is not None:
        circle_visible = Circle(
            (x, y),
            R_visible,
            facecolor=color_visible,
            edgecolor=None,
            alpha=1.0,
            zorder=zorder
        )
        ax.add_patch(circle_visible)
    # Disque noir (taille réelle)
    circle_reel = Circle(
        (x, y),
        R_reel,
        facecolor=color_point,
        edgecolor='black',
        linewidth=1.2,
        alpha=1.0,
        zorder=zorder+1
    )
    ax.add_patch(circle_reel)
    
nb_simu = len(x)
nb_cols = 3
nb_rows = math.ceil(nb_simu / nb_cols)


fig, axes = plt.subplots(nb_rows, nb_cols, figsize=(15, 5 * nb_rows))

if nb_rows == 1:
    axes = [axes]

nb_particules = len(x[0])


for i in range(nb_simu):
    row = i // nb_cols
    col = i % nb_cols
    ax = axes[row][col] if nb_rows > 1 else axes[col]

    # Trajectoires
    ax.plot(x[i][0], y[i][0], color='blue', zorder=3, label='Sonde')
    ax.plot(x[i][1], y[i][1], color='red', zorder=3, label='Terre')
    ax.plot(x[i][2], y[i][2], color='green', zorder=3, label='Lune')

    # Positions finales
    x_terre = x[i][1][-1]
    y_terre = y[i][1][-1]
    x_lune = x[i][2][-1]
    y_lune = y[i][2][-1]

    # Terre : disque orange visible + point noir taille réelle
    draw_body(ax, x_terre, y_terre, R_T, color_visible=couleur, R_visible=R_T*u, color_point='black', zorder=1)
    # Lune : disque vert visible + point noir taille réelle
    draw_body(ax, x_lune, y_lune, R_L, color_visible=couleur2, R_visible=R_L*u_L, color_point='black', zorder=1)

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

    fig, ax = plt.subplots(figsize=(6,6))

    nb_particules = len(x[i])

    for j in range(nb_particules):
        ax.plot(x[i][j], y[i][j], color=colors[j], zorder=3)

    # positions finales
    x_terre = x[i][1][-1]
    y_terre = y[i][1][-1]

    x_lune = x[i][2][-1]
    y_lune = y[i][2][-1]

    # Terre : disque orange visible + point noir taille réelle
    draw_body(ax, x_terre, y_terre, R_T, color_visible=couleur, R_visible=R_T*u, color_point='black', zorder=1)
    # Lune : disque vert visible + point noir taille réelle
    draw_body(ax, x_lune, y_lune, R_L, color_visible=couleur2, R_visible=R_L*u_L, color_point='black', zorder=1)

    ax.set_xlabel(rf"$x \,\,[\mathrm{{m}}]$")
    ax.set_ylabel(rf"$y \,\,[\mathrm{{m}}]$")
    ax.set_title(rf"$\alpha_0={parameterscan[i]:.3f}^\circ$")

    ax.grid(True, alpha=0.5)
    ax.set_aspect('equal')
    ax.get_legend()

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"trajectoire_simulation_{i}.png"), dpi=300)
    plt.close()
    
    
plot_simulation(0)

# ============================================================
# hauteur 
# ============================================================
'''
def min_list(x):
    return np.array([np.min(xi) for xi in x])

def max_list(x):
    return np.array([np.max(xi) for xi in x])

'''

r = []

for i in range(len(x)):

    r0 = np.array([x[i][0], y[i][0]])   # (2, N)
    r1 = np.array([x[i][1], y[i][1]])   # (2, N)

    ri = np.linalg.norm(r0 - r1, axis=0)  # norme pour chaque k

    r.append(ri)

#r = np.array(r, dtype=object)  # (len(x), N)
'''
#h_min = np.min(h, axis=1)
h_min = min_list(h)
#h_min = min_list([hi for hi in h])
'''
r = np.array(r, dtype=object)
h = r-R_T
h_min = []
for i in range(len(h)):
    hmin = np.min(h[i])
    h_min.append(hmin)
    
h_min = np.array(h_min)


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
#acc_max = max_list(acc[:, 0])
#acc_max = max_list([acc_i[0] for acc_i in acc])
#acc_max = np.max(acc[:, 0, :], axis=1)

acc_max = []
for i in range(len(acc)):
    acc_max_i = np.max(acc[i][0])  # max de cette simulation 
    acc_max.append(acc_max_i)  # on ajoute à la liste globale

acc_max = np.array(acc_max) # acceleration max de la sonde sur simulation i 

plt.plot(parameterscan, acc_max, 'o-', color= 'red')
plt.xlabel(rf"$\alpha_0 \,\,[\mathrm{{°}}]$")
plt.ylabel(rf"$a_{{max}} \,\,[\mathrm{{m \cdot s^{{-2}} }}]$")
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "acceleration_max.png"), dpi=300)
plt.close()
