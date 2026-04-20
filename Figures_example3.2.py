import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import re

# ============================================================
# USER SETTINGS
# ============================================================

folder = r"C:/EPFL/Semestre_4/Physique_numérique/Exercice_3/Partie/Exercice3/Scan_epsilon_Gravit_s_0.9_epsilon_1e+03"
R_T = 6.3781e6   # rayon de la Terre (m)

# Références analytiques
h_ref = 1.0e4
vmax_ref = 1.1122e4

# ============================================================
# Output folder
# ============================================================

fig_dir = os.path.join(folder, "figures_c")
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

# Tri croissant en epsilon
order = np.argsort(epsilon_values)
epsilon_values = np.array(epsilon_values)[order]
datasets = [datasets[i] for i in order]

# ============================================================
# Compute h_min, v_max, number of saved points
# ============================================================

hmins = []
vmaxs = []
npoints = []

for data in datasets:
    t  = data[:, 0]
    x  = data[:, 1]
    y  = data[:, 2]
    vx = data[:, 3]
    vy = data[:, 4]

    r = np.sqrt(x**2 + y**2)
    h = r - R_T
    v = np.sqrt(vx**2 + vy**2)

    hmins.append(np.min(h))
    vmaxs.append(np.max(v))
    npoints.append(len(t))

hmins = np.array(hmins)
vmaxs = np.array(vmaxs)
npoints = np.array(npoints)

# ============================================================
# Errors
# ============================================================

err_h = np.abs(hmins - h_ref)
err_v = np.abs(vmaxs - vmax_ref)

rel_err_h = err_h / abs(h_ref)
rel_err_v = err_v / abs(vmax_ref)

# ============================================================
# Print numerical values
# ============================================================

print("\nResults:")
for eps, hmin, vmax, eh, ev, n in zip(epsilon_values, hmins, vmaxs, rel_err_h, rel_err_v, npoints):
    print(f"epsilon = {eps:12.6g}   h_min = {hmin:12.6f} m   "
          f"vmax = {vmax:12.6f} m/s   rel_err_h = {eh:.3e}   "
          f"rel_err_v = {ev:.3e}   saved_points = {n}")

# ============================================================
# Terre pour la trajectoire
# ============================================================

theta = np.linspace(0, 2*np.pi, 500)
xT = R_T * np.cos(theta)
yT = R_T * np.sin(theta)

# ============================================================
# Plot 1 : toutes les trajectoires
# ============================================================

plt.figure()

for data, eps in zip(datasets, epsilon_values):
    x = data[:, 1]
    y = data[:, 2]
    plt.plot(x, y, label=f"eps={eps:.0e}")

plt.plot(xT, yT, 'k', label="Terre")
plt.axis("equal")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Trajectoire de la sonde (pas adaptatif)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "trajectory_all.png"), dpi=300)

# ============================================================
# Plot 2 : altitude vs time pour toutes les epsilon
# ============================================================

plt.figure()

for data, eps in zip(datasets, epsilon_values):
    t = data[:, 0]
    x = data[:, 1]
    y = data[:, 2]
    h = np.sqrt(x**2 + y**2) - R_T

    plt.plot(t/3600, h, label=f"eps={eps:.0e}")

plt.xlabel("t (h)")
plt.ylabel("h (m)")
plt.title("Altitude de la sonde (pas adaptatif)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "altitude_vs_time_all.png"), dpi=300)

# ============================================================
# Plot 3 : vitesse vs time pour toutes les epsilon
# ============================================================

plt.figure()

for data, eps in zip(datasets, epsilon_values):
    t  = data[:, 0]
    vx = data[:, 3]
    vy = data[:, 4]
    v = np.sqrt(vx**2 + vy**2)

    plt.plot(t/3600, v, label=f"eps={eps:.0e}")

plt.xlabel("t (h)")
plt.ylabel("v (m/s)")
plt.title("Vitesse de la sonde (pas adaptatif)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "speed_vs_time_all.png"), dpi=300)

# ============================================================
# Plot 4 : convergence de h_min en fonction de epsilon
# ============================================================

plt.figure()
plt.plot(epsilon_values, hmins, 'o-', label=r"$h_{\min}$ numérique")
plt.axhline(h_ref, linestyle='--', label=r"$h_{\min}$ analytique")
plt.xscale("log")
plt.xlabel(r"$\epsilon$")
plt.ylabel(r"$h_{\min}$ (m)")
plt.title(r"Convergence de $h_{\min}$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "hmin_vs_epsilon.png"), dpi=300)

# ============================================================
# Plot 5 : convergence de vmax en fonction de epsilon
# ============================================================

plt.figure()
plt.plot(epsilon_values, vmaxs, 'o-', label=r"$v_{\max}$ numérique")
plt.axhline(vmax_ref, linestyle='--', label=r"$v_{\max}$ analytique")
plt.xscale("log")
plt.xlabel(r"$\epsilon$")
plt.ylabel(r"$v_{\max}$ (m/s)")
plt.title(r"Convergence de $v_{\max}$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "vmax_vs_epsilon.png"), dpi=300)

# ============================================================
# Plot 6 : erreur relative sur h_min
# ============================================================

plt.figure()
plt.plot(epsilon_values, rel_err_h, 'o-', label=r"Erreur relative sur $h_{\min}$")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\epsilon$")
plt.ylabel("Erreur relative")
plt.title(r"Erreur relative sur $h_{\min}$")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "relerr_hmin_vs_epsilon.png"), dpi=300)

# ============================================================
# Plot 7 : erreur relative sur vmax
# ============================================================

plt.figure()
plt.plot(epsilon_values, rel_err_v, 'o-', label=r"Erreur relative sur $v_{\max}$")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\epsilon$")
plt.ylabel("Erreur relative")
plt.title(r"Erreur relative sur $v_{\max}$")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "relerr_vmax_vs_epsilon.png"), dpi=300)

# ============================================================
# Plot 8 : nombre de points sauvegardés vs epsilon
# Attention : ce n'est un bon proxy du nombre de pas
# que si sampling = 1
# ============================================================

plt.figure()
plt.plot(epsilon_values, npoints, 'o-')
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\epsilon$")
plt.ylabel("Nombre de points sauvegardés")
plt.title("Coût numérique en fonction de $\epsilon$")
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "npoints_vs_epsilon.png"), dpi=300)

# ============================================================
# Plot 9 : dt(t) pour toutes les epsilon
# On suppose que la colonne 7 est dt_current
# ============================================================

plt.figure()

for data, eps in zip(datasets, epsilon_values):
    t = data[:, 0]
    dt_current = data[:, 7]

    plt.plot(t/3600, dt_current, label=f"eps={eps:.0e}")

plt.xlabel("t (h)")
plt.ylabel(r"$\Delta t$ (s)")
plt.title(r"Pas de temps adaptatif $\Delta t(t)$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "dt_vs_time_all.png"), dpi=300)

# ============================================================
# Plot 10 : comparaison distance Terre-sonde et dt(t)
# Ici on prend le run le plus précis = plus petit epsilon
# ============================================================

best_index = 0
best_data = datasets[best_index]
best_eps = epsilon_values[best_index]

t = best_data[:, 0]
x = best_data[:, 1]
y = best_data[:, 2]
r = np.sqrt(x**2 + y**2)
h = r - R_T
dt_current = best_data[:, 7]

fig, ax1 = plt.subplots()

ax1.plot(t/3600, h, label="Altitude", linewidth=2)
ax1.set_xlabel("t (h)")
ax1.set_ylabel("h (m)")
ax1.grid(True)

ax2 = ax1.twinx()
ax2.plot(t/3600, dt_current, '--', label=r"$\Delta t$", linewidth=2)
ax2.set_ylabel(r"$\Delta t$ (s)")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")

plt.title(f"Altitude et pas adaptatif (eps={best_eps:.0e})")
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "altitude_and_dt_vs_time_best.png"), dpi=300)

# ============================================================
# Plot 11 : zoom près du périgée pour le run le plus précis
# ============================================================

imin = np.argmin(h)
tmin = t[imin] / 3600

plt.figure()
plt.plot(t/3600, h, label="Altitude")
plt.axvline(tmin, linestyle='--', label="Périgée")
plt.xlim(max(0, tmin - 2), tmin + 2)
plt.ylim(bottom=0)
plt.xlabel("t (h)")
plt.ylabel("h (m)")
plt.title(f"Zoom sur le périgée (eps={best_eps:.0e})")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "altitude_zoom_perigee_best.png"), dpi=300)

print("\nFigures saved in:", fig_dir)