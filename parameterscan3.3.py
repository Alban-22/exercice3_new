import numpy as np
import subprocess
import os

G = 6.67430e-11
M_T = 5.972e24

def compute_vx0_vy0_from_h(
    h_target,
    r0=314159e3,
    v0=1200.0,
    R_T=6.3781e6,
    clockwise=False
):
    mu = G * M_T
    rp = R_T + h_target

    if rp <= 0.0:
        raise ValueError("rp <= 0")
    if rp >= r0:
        raise ValueError("rp >= r0")

    vp2 = v0**2 + 2.0 * mu * (1.0 / rp - 1.0 / r0)
    vp = np.sqrt(vp2)

    vt0 = (rp / r0) * vp
    vr0 = -np.sqrt(max(v0**2 - vt0**2, 0.0))

    vx0 = vr0
    vy0 = vt0 if not clockwise else -vt0
    return vx0, vy0

# Parameters
repertoire = 'C:/EPFL/Semestre_4/Physique_numérique/Exercice_3/Partie/Exercice3/'
executable = 'engine.exe'
input_filename = 'configuration.in.example'

input_parameters = {
    'N': 2,
    'n': 1,
    'tf': 2 * 24 * 3600,
    'dt': 10,
    'adaptive': 1,
    'm2': 7.3477e8,
    'epsilon': 1e4,
    's': 0.9,
    'sampling': 10,
    'rho_0': 1.2,          # atmosphère active pour 3.3.b
    'R_T': 6.3781e6,
    'lambda': 7238.2,
    'Cx0': 0,
    'Cx1': 0.3,
    'Cx2': 0
}

# ici on scanne h_target, pas epsilon
h_values = np.array([1e4, 2e4, 5e4, 1e5, 2e5, 5e5])

outdir = "Scan_htarget_atmosphere"
os.makedirs(outdir, exist_ok=True)

for h_target in h_values:
    params = input_parameters.copy()

    vx0, vy0 = compute_vx0_vy0_from_h(
        h_target=h_target,
        r0=314159e3,
        v0=1200.0,
        R_T=params['R_T'],
        clockwise=False
    )

    params['vx0'] = vx0
    params['vy0'] = vy0

    output_file = f"run_htarget_{h_target:.6g}.txt"
    output_path = os.path.join(outdir, output_file)

    param_string = " ".join(f"{k}={v:.15g}" for k, v in params.items())

    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{param_string} output={output_path}"
    )

    print(f"h_target = {h_target:.6g} m, vx0 = {vx0:.6f}, vy0 = {vy0:.6f}")
    subprocess.run(cmd, shell=True)
    print("Done.")