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
    raise ValueError(f"No dataset found in {folder}")

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
# NEW FORMAT PARSING
# ============================================================

all_q = []

for data in datasets:
    t = data[:, 0]
    ncols = data.shape[1]

    # last 5 columns: energies + momentum + dt
    E_c = data[:, -5]
    E_p = data[:, -4]
    E_meca = data[:, -3]
    p_norm = data[:, -2]
    dt_current = data[:, -1]

    # particle block
    remaining = ncols - 1 - 5
    if remaining % 5 != 0:
        raise ValueError("Invalid file format: particle block not divisible by 5")

    n_particles = remaining // 5

    # take particle 0 as main object
    i = 0
    base = 1 + 5 * i

    x = data[:, base + 0]
    y = data[:, base + 1]
    vx = data[:, base + 2]
    vy = data[:, base + 3]
    a_tot = data[:, base + 4]

    r = np.sqrt(x**2 + y**2)
    h = r - R_T
    v = np.sqrt(vx**2 + vy**2)

    all_q.append(
        {
            "t": t,
            "x": x,
            "y": y,
            "h": h,
            "v": v,
            "E_meca": E_meca,
            "a_tot": a_tot,
            "F_frott": a_tot,
        }
    )

# ============================================================
# Quantities
# ============================================================

hmins = np.array([np.min(q["h"]) for q in all_q])
vmaxs = np.array([np.max(q["v"]) for q in all_q])
amaxs = np.array([np.nanmax(q["a_tot"]) for q in all_q])
fmaxs = np.array([np.nanmax(q["F_frott"]) for q in all_q])
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
        f"|a|max = {fmax:12.6e}   "
        f"saved_points = {n}   t_final = {tfinal:.6g} s"
    )

# ============================================================
# Collision detection
# ============================================================

tf_ref = 2 * 24 * 3600
collision_mask = final_times < tf_ref
nocollision_mask = ~collision_mask

if np.any(collision_mask):
    first_collision = np.min(param_values[collision_mask])
    print(f"\nPlus petit paramètre avec collision: {first_collision}")
else:
    first_collision = None
    print("\nAucune collision detectee.")

# ============================================================
# PLOT 1: hmin
# ============================================================

plt.figure()
plt.scatter(param_values[nocollision_mask], hmins[nocollision_mask],
            label="Sans collision")
plt.scatter(param_values[collision_mask], hmins[collision_mask],
            label="Collision")

plt.axhline(0.0, linestyle="--", color="k", label="Surface Terre")

if first_collision is not None:
    plt.axvline(first_collision, linestyle=":", color="orange",
                label=f"Collision = {first_collision}")

plt.xlabel(param_label)
plt.ylabel(r"$h_{\min}$ (m)")
plt.title("Altitude minimale")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, f"hmin_vs_{param_name}.png"), dpi=300)

# ============================================================
# PLOT 2: Trajectoires
# ============================================================

plt.figure()
theta = np.linspace(0, 2 * np.pi, 500)

for q, value in zip(all_q, param_values):
    plt.plot(q["x"], q["y"], label=f"{value:g}")

plt.plot(R_T * np.cos(theta), R_T * np.sin(theta), "k")

plt.axis("equal")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Trajectoires")
plt.grid()
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "trajectory_all.png"), dpi=300)

# ============================================================
# ZOOM EARTH
# ============================================================

plt.figure()

for q in all_q:
    plt.plot(q["x"], q["y"])

plt.plot(R_T * np.cos(theta), R_T * np.sin(theta), "k")

lim = R_T + 2e6
plt.axis("equal")
plt.xlim(-lim, lim)
plt.ylim(-lim, lim)

plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Zoom Terre")
plt.grid()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "trajectory_zoom_earth.png"), dpi=300)

# ============================================================
# PLOT 3: amax
# ============================================================

plt.figure()

plt.plot(param_values, amaxs, color="0.7")

plt.scatter(param_values[~collision_mask], amaxs[~collision_mask],
            label="Sans collision")

plt.scatter(param_values[collision_mask], amaxs[collision_mask],
            label="Collision")

if first_collision is not None:
    plt.axvline(first_collision, linestyle=":", color="orange")

plt.xlabel(param_label)
plt.ylabel(r"$a_{max}$")
plt.yscale("log")
plt.title("Acceleration max")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, f"amax_vs_{param_name}.png"), dpi=300)

print("\nFigures saved in:", fig_dir)