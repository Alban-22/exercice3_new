import numpy as np
import subprocess
import os

# Parameters
repertoire = 'C:/EPFL/Semestre_4/Physique_numérique/exercice3_new/'
executable = 'engine2.exe'
input_filename = 'configuration.in.example' # Strictly no longer needed, but we keep it for now to avoid having to change the code in engine.cpp



input_parameters = {
    'N': 2,
    'n': 1,
    'tf': 1 * 24 * 3600,
    'dt': 10,
    'adaptive': 1,
    'epsilon': 1e-6,
    's': 0.9,
    'sampling': 1,
    'rho_0': 0.0,
    'alphadeg0': 5,
    'alphadeg1': 0,
    'alphadeg2': 0,
    'bool_alpha0': True,    
    'bool_alpha1': False,
    'bool_alpha2': False,
    'vnorm0': 1200
}

# -------------------------------------------------

# Updated from last time, the code below can now be used to scan any parameter, just make sure to update the paramstr and the variable_array accordingly

paramstr = 'alphadeg0' # The parameter to scan, must be one of the keys in input_parameters

variable_array =np.linspace(180 - 2/60,180+2/60, 5)


outstr = f"Gravit_s_{input_parameters['s']:.2g}_epsilon_{input_parameters['epsilon']:.2g}"

# -------------------------------------------------
# Create output directory (2 significant digits)
# -------------------------------------------------
outdir = f"Scan_{paramstr}_{outstr}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)


for i in range(len(variable_array)):

    # Copy parameters and overwrite scanned one
    params = input_parameters.copy()
    params[paramstr] = variable_array[i]

    output_file = f"{outstr}_{paramstr}_{variable_array[i]}.txt"
    output_path = os.path.join(outdir, output_file)

    # Build parameter string
    param_string = " ".join(f"{k}={v:.15g}" for k, v in params.items())

    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{param_string} output={output_path}"
    )

    print(cmd)
    subprocess.run(cmd, shell=True)
    print("Done.")


