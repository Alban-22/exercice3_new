import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import re

# ============================================================
# USER SETTINGS
# ============================================================

folder = r"C:/EPFL/Semestre_4/Physique_numérique/Exercice_3/Partie/Exercice3/Scan_epsilon_Gravit3_s_0.9_epsilon_1e+03"
R_T = 6.3781e6          # rayon Terre (m)
rho_0 = 1.2             # densité au niveau de la mer (kg/m^3)
lamb = 7238.2           # hauteur caractéristique (m)
Cx = 0.3                # coefficient de traînée
mA = 8500.0             # masse Artemis (kg)
dA = 5.02               # diamètre Artemis (m)

G = 6.674e-11
M_T = 5.972e24
mu = G * M_T

S = np.pi * (dA/2.0)**2

# Références éventuelles
h_ref = 1.0e4
vmax_ref = 1.1122e4

# ============================================================
# Output folder
# ============================================================

fig_dir = os.path.join(folder, "figures_33a")
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
# Fonctions physiques
# ============================================================

def rho_atm(r):
    return rho_0 * np.exp(-(r - R_T) / lamb)

def compute_quantities(data):
    t  = data[:, 0]
    x  = data[:, 1]
    y  = data[:, 2]
    vx = data[:, 3]
    vy = data[:, 4]

    r = np.sqrt(x**2 + y**2)
    h = r - R_T
    v = np.sqrt(vx**2 + vy**2)

    r_safe = np.maximum(r, 1e-12)
    v_safe = np.maximum(v, 1e-12)

    rho = rho_atm(r)

    # Accélération gravitationnelle
    ax_g = -mu * x / r_safe**3
    ay_g = -mu * y / r_safe**3

    # Accélération de traînée
    factor_drag = -0.5 * rho * S * Cx * v_safe / mA
    ax_d = factor_drag * vx
    ay_d = factor_drag * vy

    # Accélération totale
    ax_tot = ax_g + ax_d
    ay_tot = ay_g + ay_d

    a_tot = np.sqrt(ax_tot**2 + ay_tot**2)
    a_drag = np.sqrt(ax_d**2 + ay_d**2)

    # Force de traînée
    Fx_d = mA * ax_d
    Fy_d = mA * ay_d

    # Puissance de traînée : P = F.v
    P_drag = Fx_d * vx + Fy_d * vy

    return {
        "t": t,
        "x": x,
        "y": y,
        "vx": vx,
        "vy": vy,
        "r": r,
        "h": h,
        "v": v,
        "rho": rho,
        "a_tot": a_tot,
        "a_drag": a_drag,
        "P_drag": P_drag
    }

# ============================================================
# Compute global quantities
# ============================================================

all_q = []

hmins = []
vmaxs = []
amaxs = []
pmaxs = []
npoints = []

for data in datasets:
    q = compute_quantities(data)
    all_q.append(q)

    hmins.append(np.min(q["h"]))
    vmaxs.append(np.max(q["v"]))
    amaxs.append(np.max(q["a_tot"]))
    pmaxs.append(np.max(np.abs(q["P_drag"])))
    npoints.append(len(q["t"]))

hmins = np.array(hmins)
vmaxs = np.array(vmaxs)
amaxs = np.array(amaxs)
pmaxs = np.array(pmaxs)
npoints = np.array(npoints)

# ============================================================
# Errors
# ============================================================

rel_err_h = np.abs(hmins - h_ref) / abs(h_ref)
rel_err_v = np.abs(vmaxs - vmax_ref) / abs(vmax_ref)

# ============================================================
# Print numerical values
# ============================================================

print("\nResults:")
for eps, hmin, vmax, amax, pmax, eh, ev, n in zip(
    epsilon_values, hmins, vmaxs, amaxs, pmaxs, rel_err_h, rel_err_v, npoints
):
    print(
        f"epsilon = {eps:12.6g}   "
        f"h_min = {hmin:12.6f} m   "
        f"vmax = {vmax:12.6f} m/s   "
        f"a_max = {amax:12.6f} m/s^2   "
        f"|P_drag|max = {pmax:12.6e} W   "
        f"rel_err_h = {eh:.3e}   rel_err_v = {ev:.3e}   "
        f"saved_points = {n}"
    )

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

for q, eps in zip(all_q, epsilon_values):
    plt.plot(q["x"], q["y"], label=f"eps={eps:.0e}")

plt.plot(xT, yT, 'k', label="Terre")
plt.axis("equal")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Trajectoire de la sonde (3.3.a)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "trajectory_all.png"), dpi=300)

# ============================================================
# Plot 2 : altitude vs time
# ============================================================

plt.figure()

for q, eps in zip(all_q, epsilon_values):
    plt.plot(q["t"]/3600, q["h"], label=f"eps={eps:.0e}")

plt.xlabel("t (h)")
plt.ylabel("h (m)")
plt.title("Altitude de la sonde")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "altitude_vs_time_all.png"), dpi=300)

# ============================================================
# Plot 3 : vitesse vs time
# ============================================================

plt.figure()

for q, eps in zip(all_q, epsilon_values):
    plt.plot(q["t"]/3600, q["v"], label=f"eps={eps:.0e}")

plt.xlabel("t (h)")
plt.ylabel("v (m/s)")
plt.title("Vitesse de la sonde")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "speed_vs_time_all.png"), dpi=300)

# ============================================================
# Plot 4 : accélération totale vs time
# ============================================================

plt.figure()

for q, eps in zip(all_q, epsilon_values):
    plt.plot(q["t"]/3600, q["a_tot"], label=f"eps={eps:.0e}")

plt.xlabel("t (h)")
plt.ylabel(r"$a(t)$ (m/s$^2$)")
plt.title("Accélération totale de la sonde")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "acceleration_vs_time_all.png"), dpi=300)

# ============================================================
# Plot 5 : puissance de traînée vs time
# ============================================================

plt.figure()

for q, eps in zip(all_q, epsilon_values):
    plt.plot(q["t"]/3600, q["P_drag"], label=f"eps={eps:.0e}")

plt.xlabel("t (h)")
plt.ylabel(r"$P_{\mathrm{drag}}(t)$ (W)")
plt.title("Puissance instantanée de la traînée")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "power_drag_vs_time_all.png"), dpi=300)

# ============================================================
# Plot 6 : h_min vs epsilon
# ============================================================

plt.figure()
plt.plot(epsilon_values, hmins, 'o-', label=r"$h_{\min}$ numérique")
plt.axhline(h_ref, linestyle='--', label=r"$h_{\min}$ référence")
plt.xscale("log")
plt.xlabel(r"$\epsilon$")
plt.ylabel(r"$h_{\min}$ (m)")
plt.title(r"Convergence de $h_{\min}$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "hmin_vs_epsilon.png"), dpi=300)

# ============================================================
# Plot 7 : vmax vs epsilon
# ============================================================

plt.figure()
plt.plot(epsilon_values, vmaxs, 'o-', label=r"$v_{\max}$ numérique")
plt.axhline(vmax_ref, linestyle='--', label=r"$v_{\max}$ référence")
plt.xscale("log")
plt.xlabel(r"$\epsilon$")
plt.ylabel(r"$v_{\max}$ (m/s)")
plt.title(r"Convergence de $v_{\max}$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "vmax_vs_epsilon.png"), dpi=300)

# ============================================================
# Plot 8 : accélération maximale vs epsilon
# ============================================================

plt.figure()
plt.plot(epsilon_values, amaxs, 'o-')
plt.xscale("log")
plt.xlabel(r"$\epsilon$")
plt.ylabel(r"$a_{\max}$ (m/s$^2$)")
plt.title(r"Convergence de l'accélération maximale")
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "amax_vs_epsilon.png"), dpi=300)

# ============================================================
# Plot 9 : puissance de traînée max vs epsilon
# ============================================================

plt.figure()
plt.plot(epsilon_values, pmaxs, 'o-')
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\epsilon$")
plt.ylabel(r"$|P_{\mathrm{drag}}|_{\max}$ (W)")
plt.title(r"Convergence de la puissance maximale de traînée")
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "pdragmax_vs_epsilon.png"), dpi=300)

# ============================================================
# Plot 10 : nombre de points sauvegardés vs epsilon
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
# Plot 11 : run le plus précis
# ============================================================

best_index = 0
best_q = all_q[best_index]
best_eps = epsilon_values[best_index]

t = best_q["t"]
h = best_q["h"]
a_tot = best_q["a_tot"]
P_drag = best_q["P_drag"]

fig, ax1 = plt.subplots()

ax1.plot(t/3600, h, label="Altitude", linewidth=2)
ax1.set_xlabel("t (h)")
ax1.set_ylabel("h (m)")
ax1.grid(True)

ax2 = ax1.twinx()
ax2.plot(t/3600, a_tot, '--', label="Accélération", linewidth=2)
ax2.set_ylabel(r"$a$ (m/s$^2$)")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")

plt.title(f"Altitude et accélération (eps={best_eps:.0e})")
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "altitude_and_acceleration_best.png"), dpi=300)

plt.figure()
plt.plot(t/3600, P_drag, label="Puissance de traînée")
plt.xlabel("t (h)")
plt.ylabel(r"$P_{\mathrm{drag}}(t)$ (W)")
plt.title(f"Puissance de traînée (eps={best_eps:.0e})")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "power_drag_best.png"), dpi=300)

print("\nFigures saved in:", fig_dir)