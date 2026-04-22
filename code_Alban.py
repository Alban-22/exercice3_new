import glob
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# Settings
# ============================================================

base_dir = os.path.dirname(os.path.abspath(__file__))
folder = os.path.join(base_dir, "Scan_alphadeg0_ATM3_s_0.9_epsilon_1e-05")



if not glob.glob(os.path.join(folder, "*.txt")):
    scan_candidates = sorted(
        glob.glob(os.path.join(base_dir, "Scan_alphadeg0_ATM3_s_0.9_epsilon_1e-05")),
        key=os.path.getmtime,
        reverse=True,
    )
    if scan_candidates:
        folder = scan_candidates[0]

R_T = 6.3781e6

fig_dir = os.path.join(folder, "figures_33b")
os.makedirs(fig_dir, exist_ok=True)


def label_for_param(param_name):
    labels = {
        "epsilon": r"$\epsilon$",
        "dt": r"$\Delta t$ (s)",
        "angle_deg": "Angle de lancement (deg)",
        "rho_0": r"$\rho_0$ (kg/m$^3$)",
        "s": "s",
    }
    return labels.get(param_name, param_name)


def parse_scan_file(path):
    name = os.path.basename(path)[:-4]
    parts = name.split("_")

    for split_idx in range(len(parts) - 1, 0, -1):
        value_str = "_".join(parts[split_idx:])
        try:
            value = float(value_str)
        except ValueError:
            continue

        param_name = parts[split_idx - 1]
        if split_idx >= 2 and parts[split_idx - 2] == "angle":
            param_name = "angle_deg"

        return param_name, value

    raise ValueError(f"Cannot parse scan parameter from file name: {name}")


# ============================================================
# Read scan files
# ============================================================

files = sorted(glob.glob(os.path.join(folder, "*.txt")))
files = [
    f
    for f in files
    if re.match(r".*_[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?\.txt$", os.path.basename(f))
]

datasets = []
param_values = []
param_names = []

for f in files:
    param_name, value = parse_scan_file(f)
    data = np.loadtxt(f)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    datasets.append(data)
    param_values.append(value)
    param_names.append(param_name)

if len(datasets) == 0:
    raise ValueError(f"No dataset found in {folder}. Run parameterscan3.3.B.py first.")

unique_param_names = sorted(set(param_names))
if len(unique_param_names) != 1:
    raise ValueError(f"Several scanned parameters found: {unique_param_names}")

param_name = unique_param_names[0]
param_label = label_for_param(param_name)

order = np.argsort(param_values)
param_values = np.array(param_values)[order]
datasets = [datasets[i] for i in order]

print(f"Found {len(datasets)} datasets.")
print(f"Scanned parameter: {param_name}")

# ============================================================
# Quantities (compatible with engine2.cpp output)
# ============================================================

N = 2  # nombre de planètes
n = 1  # nombre de sondes
nb_corps = N + n  # 3 (sonde, Terre, Lune)

all_q = []

for data in datasets:
    t = data[:, 0]
    # Pour chaque corps, extraire x, y, vx, vy, acc
    x = [data[:, 1 + 5 * i] for i in range(nb_corps)]
    y = [data[:, 2 + 5 * i] for i in range(nb_corps)]
    vx = [data[:, 3 + 5 * i] for i in range(nb_corps)]
    vy = [data[:, 4 + 5 * i] for i in range(nb_corps)]
    acc = [data[:, 5 + 5 * i] for i in range(nb_corps)]

    # Pour la sonde (indice 0)
    r_s = np.sqrt(x[0] ** 2 + y[0] ** 2)
    h_s = r_s - R_T
    v_s = np.sqrt(vx[0] ** 2 + vy[0] ** 2)

    # Energies (fin de ligne)
    Ec = data[:, -5]
    Ep = data[:, -4]
    Em = data[:, -3]
    P = data[:, -2]
    dt = data[:, -1]

    all_q.append(
        {
            "t": t,
            "x": x,
            "y": y,
            "vx": vx,
            "vy": vy,
            "acc": acc,
            "h": h_s,
            "v": v_s,
            "Ec": Ec,
            "Ep": Ep,
            "Em": Em,
            "P": P,
            "dt": dt,
        }
    )

hmins = np.array([np.min(q["h"]) for q in all_q])
vmaxs = np.array([np.max(q["v"]) for q in all_q])
amaxs = np.array([np.nanmax(q["acc"][0]) for q in all_q])  # max acc sonde
fmaxs = np.full_like(hmins, np.nan)  # pas de F_frott dans engine2.cpp
npoints = np.array([len(q["t"]) for q in all_q])
final_times = np.array([q["t"][-1] for q in all_q])

print("\nResults:")
for value, hmin, vmax, amax, fmax, n, tfinal in zip(
    param_values, hmins, vmaxs, amaxs, fmaxs, npoints, final_times
):
    print(
        f"{param_name} = {value:12.6g}   "
        f"h_min = {hmin:12.6f} m   "
        f"vmax = {vmax:12.6f} m/s   "
        f"a_max = {amax:12.6f} m/s^2   "
        f"|F_frott|max = {fmax:12.6e} N   "
        f"saved_points = {n}   t_final = {tfinal:.6g} s"
    )

# ============================================================
# Plots versus scanned parameter
# ============================================================

collision_mask = hmins <= 0.0

tf_ref = 2 * 24 * 3600

# collision si la simu s'est arrêtée avant le temps final
collision_mask = final_times < (tf_ref - 1e-9)
nocollision_mask = ~collision_mask

if np.any(collision_mask):
    first_collision_angle = np.min(param_values[collision_mask])
    print(
        f"\nPlus petit angle avec collision: {first_collision_angle:.15f} deg "
        f"({param_name})"
    )
else:
    first_collision_angle = None
    print("\nAucune collision detectee dans ce scan.")

plt.figure()
plt.scatter(param_values[nocollision_mask], hmins[nocollision_mask],
            label="Sans collision")
plt.scatter(param_values[collision_mask], hmins[collision_mask],
            label="Collision")
plt.axhline(0.0, linestyle="--", color="k", label="Surface Terre")
if first_collision_angle is not None:
    plt.axvline(
        first_collision_angle,
        linestyle=":",
        color="tab:orange",
        label=rf"Premier angle collision = {first_collision_angle:.6f}$^\circ$",
    )
plt.xlabel(param_label)
plt.ylabel(r"$h_{\min}$ (m)")
plt.title("Altitude minimale en fonction de l'angle")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, f"hmin_vs_{param_name}.png"), dpi=300)



plt.figure()
for q, value in zip(all_q, param_values):
    # Sonde
    plt.plot(q["x"][0], q["y"][0], label=f"Sonde, {param_name}={value:g}", color="tab:blue")
    # Terre
    plt.plot(q["x"][1], q["y"][1], label=f"Terre, {param_name}={value:g}", color="tab:red", linestyle='--', alpha=0.7)
    # Lune
    plt.plot(q["x"][2], q["y"][2], label=f"Lune, {param_name}={value:g}", color="tab:green", linestyle='-.', alpha=0.7)

# Surface Terre
theta = np.linspace(0, 2 * np.pi, 500)
plt.plot(R_T * np.cos(theta), R_T * np.sin(theta), "k", label="Terre (surface)")
plt.axis("equal")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Trajectoires (sonde, Terre, Lune)")
plt.grid(True)
#plt.legend(fontsize=9, loc='best', ncol=2)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "trajectory_all.png"), dpi=300)

# Zoom sur la Terre
plt.figure()
for q, value in zip(all_q, param_values):
    plt.plot(q["x"][0], q["y"][0], label=f"Sonde, {param_name}={value:g}", color="tab:blue")
    plt.plot(q["x"][1], q["y"][1], label=f"Terre, {param_name}={value:g}", color="tab:red", linestyle='--', alpha=0.7)
    plt.plot(q["x"][2], q["y"][2], label=f"Lune, {param_name}={value:g}", color="tab:green", linestyle='-.', alpha=0.7)
plt.plot(R_T * np.cos(theta), R_T * np.sin(theta), "k", label="Terre (surface)")
earth_zoom_margin = 2.0e6
earth_zoom_lim = R_T + earth_zoom_margin
plt.axis("equal")
plt.xlim(-earth_zoom_lim, earth_zoom_lim)
plt.ylim(-earth_zoom_lim, earth_zoom_lim)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Trajectoires - zoom autour de la Terre")
plt.grid(True)
plt.legend(fontsize=8, loc='best', ncol=2)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "trajectory_zoom_earth.png"), dpi=300)

print("\nFigures saved in:", fig_dir)