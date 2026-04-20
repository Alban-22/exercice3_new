import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import re

# ============================================================
# USER SETTINGS
# ============================================================

folder = r"C:/EPFL/Semestre_4/Physique_numérique/Exercice_3/Partie/Exercice3/Scan_dt_Gravit_s_0.9_epsilon_0_rk4_order_fine/"
R_T = 6.3781e6   # rayon de la Terre (m)

use_interpolated_reference = True
h_ref_manual = 1.0e4
vmax_ref_manual = 1.1122e4

# ============================================================
# Output folder
# ============================================================

fig_dir = os.path.join(folder, "figures_b")
os.makedirs(fig_dir, exist_ok=True)
zoom_dir = os.path.join(fig_dir, "trajectory_zoom_earth_by_dt")
os.makedirs(zoom_dir, exist_ok=True)

# ============================================================
# Scan files
# ============================================================

files = sorted(glob.glob(os.path.join(folder, "*.txt")))

files = [
    f for f in files
    if re.match(r".*_[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?\.txt$", os.path.basename(f))
]

datasets = []
dt_values = []

for f in files:
    name = os.path.basename(f)[:-4]
    parts = name.split("_")

    param_name = parts[-2]
    value = float(parts[-1])

    if param_name != "dt":
        continue

    data = np.loadtxt(f)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    datasets.append(data)
    dt_values.append(value)

print(f"Found {len(datasets)} datasets.")

if len(datasets) == 0:
    raise ValueError("Aucun dataset trouvé.")

# tri croissant
order = np.argsort(dt_values)
dt_values = np.array(dt_values)[order]
datasets = [datasets[i] for i in order]

# ============================================================
# Compute h_min and v_max
# ============================================================

def quadratic_extremum(t, q, mode="min"):
    """Return an interpolated extremum from a local parabola through 3 points."""
    if len(q) < 3:
        idx = np.argmin(q) if mode == "min" else np.argmax(q)
        return q[idx], t[idx]

    idx = np.argmin(q) if mode == "min" else np.argmax(q)

    if idx == 0 or idx == len(q) - 1:
        return q[idx], t[idx]

    t_loc = t[idx-1:idx+2]
    q_loc = q[idx-1:idx+2]

    a, b, c = np.polyfit(t_loc, q_loc, 2)

    if abs(a) < 1e-30:
        return q[idx], t[idx]

    t_star = -b / (2.0 * a)

    # Use the interpolated value only if the vertex is inside the bracket
    # and the parabola has the expected curvature.
    good_min = (mode == "min" and a > 0.0)
    good_max = (mode == "max" and a < 0.0)

    if t_loc[0] <= t_star <= t_loc[-1] and (good_min or good_max):
        q_star = a * t_star**2 + b * t_star + c
        return q_star, t_star

    return q[idx], t[idx]

hmins = []
vmaxs = []
t_hmins = []
t_vmaxs = []
final_states = []

for data in datasets:
    t = data[:,0]
    x = data[:,1]
    y = data[:,2]
    vx = data[:,3]
    vy = data[:,4]

    r = np.sqrt(x**2 + y**2)
    h = r - R_T
    v = np.sqrt(vx**2 + vy**2)

    hmin_interp, t_hmin_interp = quadratic_extremum(t, h, mode="min")
    vmax_interp, t_vmax_interp = quadratic_extremum(t, v, mode="max")

    hmins.append(hmin_interp)
    vmaxs.append(vmax_interp)
    t_hmins.append(t_hmin_interp)
    t_vmaxs.append(t_vmax_interp)
    final_states.append(data[-1, 1:5])

hmins = np.array(hmins)
vmaxs = np.array(vmaxs)
t_hmins = np.array(t_hmins)
t_vmaxs = np.array(t_vmaxs)
final_states = np.array(final_states)

if use_interpolated_reference:
    h_ref = hmins[0]
    vmax_ref = vmaxs[0]
else:
    h_ref = h_ref_manual
    vmax_ref = vmax_ref_manual

# ============================================================
# Errors
# ============================================================

rel_err_h = np.abs(hmins - h_ref) / abs(h_ref)
rel_err_v = np.abs(vmaxs - vmax_ref) / abs(vmax_ref)
rel_err_h_plot = np.where(rel_err_h > 0.0, rel_err_h, np.nan)
rel_err_v_plot = np.where(rel_err_v > 0.0, rel_err_v, np.nan)
ref_label = "reference interpolee" if use_interpolated_reference else "reference analytique"

# Pour tester l'ordre de RK4, il vaut mieux comparer l'etat au meme temps final.
# On prend le plus petit dt comme reference numerique et on ne le trace pas.
ref_state = final_states[0]
ref_dt = dt_values[0]

pos_scale = max(np.linalg.norm(ref_state[0:2]), 1.0)
vel_scale = max(np.linalg.norm(ref_state[2:4]), 1.0)

scaled_state_errors = np.sqrt(
    np.sum(((final_states[:, 0:2] - ref_state[0:2]) / pos_scale)**2, axis=1)
    + np.sum(((final_states[:, 2:4] - ref_state[2:4]) / vel_scale)**2, axis=1)
)

position_errors = (
    np.linalg.norm(final_states[:, 0:2] - ref_state[0:2], axis=1) / pos_scale
)
velocity_errors = (
    np.linalg.norm(final_states[:, 2:4] - ref_state[2:4], axis=1) / vel_scale
)

conv_mask = dt_values > ref_dt
roundoff_floor = 1e-12
fit_mask = conv_mask & (scaled_state_errors > roundoff_floor)
if np.count_nonzero(fit_mask) < 2:
    fit_mask = conv_mask.copy()

def fitted_order(dt, err):
    mask = (err > 0) & np.isfinite(err)
    if np.count_nonzero(mask) < 2:
        return np.nan
    return np.polyfit(np.log(dt[mask]), np.log(err[mask]), 1)[0]

order_state = fitted_order(dt_values[fit_mask], scaled_state_errors[fit_mask])
order_pos = fitted_order(dt_values[fit_mask], position_errors[fit_mask])
order_vel = fitted_order(dt_values[fit_mask], velocity_errors[fit_mask])

# ============================================================
# Print results
# ============================================================

print("\nResults:")
print(f"Reference h_min = {h_ref:.12g} m")
print(f"Reference vmax  = {vmax_ref:.12g} m/s")
for dt, hmin, vmax, th, tv, eh, ev in zip(
    dt_values, hmins, vmaxs, t_hmins, t_vmaxs, rel_err_h, rel_err_v
):
    print(f"dt = {dt:12.6g}   h_min = {hmin:12.6f} m at {th/3600:9.4f} h   "
          f"vmax = {vmax:12.6f} m/s at {tv/3600:9.4f} h   "
          f"rel_err_h = {eh:.3e}   rel_err_v = {ev:.3e}")

print(f"\nReference numerique pour l'ordre RK4: dt_ref = {ref_dt:g} s")
print(f"Pente estimee avec erreur finale > {roundoff_floor:.0e}.")
print(f"Ordre estime sur l'etat final normalise : {order_state:.3f}")
print(f"Ordre estime sur la position finale      : {order_pos:.3f}")
print(f"Ordre estime sur la vitesse finale       : {order_vel:.3f}")

# ============================================================
# Trajectoires (TOUS les dt)
# ============================================================

theta = np.linspace(0, 2*np.pi, 500)
xT = R_T * np.cos(theta)
yT = R_T * np.sin(theta)

plt.figure()

for data, dt in zip(datasets, dt_values):
    x = data[:,1]
    y = data[:,2]
    plt.plot(x, y, label=f"dt={dt:g}")

plt.plot(xT, yT, 'k', label="Terre")
plt.axis("equal")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Trajectoire de la sonde")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "trajectory_all.png"), dpi=300)

# ============================================================
# Trajectoire zoomee autour de la Terre, une figure par dt
# ============================================================

earth_zoom_margin = 1.5e7
earth_zoom_lim = R_T + earth_zoom_margin

for data, dt in zip(datasets, dt_values):
    x = data[:,1]
    y = data[:,2]

    fig, ax = plt.subplots()
    ax.plot(x, y, linewidth=1.5, label=f"dt={dt:g}")
    ax.plot(xT, yT, 'k', linewidth=1.5, label="Terre")

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-earth_zoom_lim, earth_zoom_lim)
    ax.set_ylim(-earth_zoom_lim, earth_zoom_lim)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Trajectoire autour de la Terre, dt={dt:g} s")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()

    dt_label = f"{dt:g}".replace(".", "p").replace("-", "m")
    fig.savefig(os.path.join(zoom_dir, f"trajectory_zoom_earth_dt_{dt_label}.png"), dpi=300)
    plt.close(fig)

# ============================================================
# Altitude vs time (TOUS les dt)
# ============================================================

plt.figure()

for data, dt in zip(datasets, dt_values):
    t = data[:,0]
    x = data[:,1]
    y = data[:,2]

    h = np.sqrt(x**2 + y**2) - R_T
    plt.plot(t/3600, h, label=f"dt={dt:g}")

plt.xlabel("t (h)")
plt.ylabel("h (m)")
plt.title("Altitude de la sonde")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "altitude_vs_time_all.png"), dpi=300)

# ============================================================
# Vitesse vs time (TOUS les dt)
# ============================================================

plt.figure()

for data, dt in zip(datasets, dt_values):
    t = data[:,0]
    vx = data[:,3]
    vy = data[:,4]

    v = np.sqrt(vx**2 + vy**2)
    plt.plot(t/3600, v, label=f"dt={dt:g}")

plt.xlabel("t (h)")
plt.ylabel("v (m/s)")
plt.title("Vitesse de la sonde")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "speed_vs_time_all.png"), dpi=300)

# ============================================================
# Convergence plots
# ============================================================

plt.figure()
plt.plot(dt_values, hmins, 'o-', label=r"$h_{\min}$ numérique")
plt.axhline(h_ref, linestyle='--', label=rf"$h_{{\min}}$ {ref_label}")
plt.xscale("log")
plt.xlabel(r"$\Delta t$ (s)")
plt.ylabel(r"$h_{\min}$ (m)")
plt.title(r"Convergence de $h_{\min}$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "hmin_vs_dt.png"), dpi=300)

plt.figure()
plt.plot(dt_values, vmaxs, 'o-', label=r"$v_{\max}$ numérique")
plt.axhline(vmax_ref, linestyle='--', label=rf"$v_{{\max}}$ {ref_label}")
plt.xscale("log")
plt.xlabel(r"$\Delta t$ (s)")
plt.ylabel(r"$v_{\max}$ (m/s)")
plt.title(r"Convergence de $v_{\max}$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "vmax_vs_dt.png"), dpi=300)

plt.figure()
plt.plot(dt_values, rel_err_h_plot, 'o-', label=r"Erreur relative sur $h_{\min}$")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\Delta t$ (s)")
plt.ylabel("Erreur relative")
plt.title(r"Erreur relative sur $h_{\min}$")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "relerr_hmin_vs_dt.png"), dpi=300)

plt.figure()
plt.plot(dt_values, rel_err_v_plot, 'o-', label=r"Erreur relative sur $v_{\max}$")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\Delta t$ (s)")
plt.ylabel("Erreur relative")
plt.title(r"Erreur relative sur $v_{\max}$")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "relerr_vmax_vs_dt.png"), dpi=300)

rel_err_h = np.abs(hmins - h_ref) / abs(h_ref)
rel_err_v = np.abs(vmaxs - vmax_ref) / abs(vmax_ref)
rel_err_h_plot = np.where(rel_err_h > 0.0, rel_err_h, np.nan)
rel_err_v_plot = np.where(rel_err_v > 0.0, rel_err_v, np.nan)

def compare_err(x, p, C=1):
    return C * x**p

dt_continu = np.logspace(np.log10(dt_values.min()), np.log10(dt_values.max()), 1000)

def reference_power_law(dt_grid, dt_data, err_data, power):
    mask = (err_data > roundoff_floor) & np.isfinite(err_data)
    if not np.any(mask):
        return np.full_like(dt_grid, np.nan)
    anchor_dt = dt_data[mask][0]
    anchor_err = err_data[mask][0]
    return anchor_err * (dt_grid / anchor_dt)**power

plt.figure()
plt.plot(dt_values, rel_err_h_plot, 'o-', label="Erreur relative sur h_min")
plt.plot(dt_continu, compare_err(dt_continu, p=4, C=1e1), '-.', label=r"$O(dt^4)$")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\Delta t$")
plt.ylabel("Erreur relative")
plt.legend()
plt.grid(True, which="both")
plt.savefig(os.path.join(fig_dir, "relerr_hmin_convergence.png"), dpi=300)

plt.figure()
plt.plot(dt_values, rel_err_v_plot, 'o-', label="Erreur relative sur v_max")
plt.plot(dt_continu, compare_err(dt_continu, p=4, C=1e3), '-.', label=r"$O(dt^4)$")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\Delta t$")
plt.ylabel("Erreur relative")
plt.legend()
plt.grid(True, which="both")
plt.savefig(os.path.join(fig_dir, "relerr_vmax_convergence.png"), dpi=300)

dt_conv = dt_values[conv_mask]
state_err_conv = scaled_state_errors[conv_mask]
pos_err_conv = position_errors[conv_mask]
vel_err_conv = velocity_errors[conv_mask]
dt_continu_conv = np.logspace(np.log10(dt_conv.min()), np.log10(dt_conv.max()), 1000)

plt.figure()
plt.plot(dt_conv, state_err_conv, 'o-', label="Etat final normalise")
plt.plot(dt_conv, pos_err_conv, 's-', label="Position finale")
plt.plot(dt_conv, vel_err_conv, '^-', label="Vitesse finale")
plt.plot(
    dt_continu_conv,
    reference_power_law(dt_continu_conv, dt_conv, state_err_conv, power=4),
    '--',
    label=r"$O(dt^4)$"
)
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\Delta t$")
plt.ylabel("Erreur relative finale")
plt.title(r"Test direct de l'ordre RK4 au temps final")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "final_state_convergence.png"), dpi=300)
