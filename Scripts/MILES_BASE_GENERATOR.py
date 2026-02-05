import numpy as np
import astropy.io.fits as fits
from astropy.wcs import WCS
import os
import re

# Pasta com os FITS originais
FITS_PATH = "../Starlightv04/MILES_PADOVA00_KU_baseFe/"
# Pasta onde os .spec serão salvos
SPEC_DIR = "../Starlightv04/BasesMiles/"
# Arquivo de massas
MASS_FILE = "../Starlightv04/MILES_PADOVA00_KU_baseFe/out_mass_KU_PADOVA00"
# Arquivo final de
OUTPUT_BASE_FILE = "../Starlightv04/BaseFiles/Base.Miles.Total"

Z_SUN = 0.019  # Referência ???


def get_data_from_fits(filepath):
    with fits.open(filepath) as file:
        flux = file[0].data
        header = file[0].header
        wcs = WCS(header)
        pixel_indices = np.arange(len(flux))
        lambdas = wcs.pixel_to_world_values(pixel_indices)
        if isinstance(lambdas, tuple):
            lambdas = lambdas[0]
    return lambdas, flux


def load_mass_map(filepath):
    mass_map = {}
    if not os.path.exists(filepath):
        print(f"AVISO: Arquivo de massas {filepath} não encontrado. Usando Mstar=1.0")
        return {}
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            mh = round(float(cols[2]), 2)
            age = round(float(cols[3]), 4)
            mstar = float(cols[5])
            mass_map[(mh, age)] = mstar
    return mass_map


def main():
    # Fixar diretório de trabalho no local do script
    script_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_path)
    os.chdir(script_dir)

    if not os.path.exists(SPEC_DIR):
        os.makedirs(SPEC_DIR)

    # 1. Conversão FITS -> SPEC
    print("--- Passo 1: Convertendo FITS para SPEC ---")
    if os.path.exists(FITS_PATH):
        fits_files = [f for f in os.listdir(FITS_PATH) if f.endswith(".fits")]
        for filename in fits_files:
            filepath = os.path.join(FITS_PATH, filename)
            output_filename = filename.replace(".fits", ".spec")
            output_filepath = os.path.join(SPEC_DIR, output_filename)

            if not os.path.exists(output_filepath):
                lambdas, flux = get_data_from_fits(filepath)
                np.savetxt(output_filepath, np.column_stack((lambdas, flux)))
        print(f"Processamento de FITS concluído.")
    else:
        print(f"ERRO: Pasta FITS {FITS_PATH} não encontrada.")
        return

    # 2. Geração do Arquivo de Base Total
    print("\n--- Passo 2: Gerando arquivo de Base Total (Starlight) ---")
    mass_map = load_mass_map(MASS_FILE)
    # Padrão MILES: Mku1.30Z[m/p][M/H]T[Idade]...
    # ai made
    pattern = re.compile(r"Mku1.30Z([mp])(\d+\.\d+)T(\d+\.\d+)")

    spec_files = sorted([f for f in os.listdir(SPEC_DIR) if f.endswith(".spec")])
    base_elements = []

    for filename in spec_files:
        match = pattern.search(filename)
        if match:
            sign, z_val, age_gyr_str = match.groups()
            mh_val = float(z_val) * (-1 if sign == "m" else 1)
            age_gyr = float(age_gyr_str)

            # Busca Mstar
            mstar = mass_map.get((round(mh_val, 2), round(age_gyr, 4)), 1.0000)

            age_yr = age_gyr * 1e9
            z_mass = Z_SUN * (10**mh_val)
            code = f"t{float(age_gyr):.4f}_z{sign}{z_val}"

            base_elements.append(
                {
                    "file": filename,
                    "age": f"{age_yr:.5e}",
                    "z": f"{z_mass:.5f}",
                    "code": code[:15],
                    "mstar": f"{mstar:.4f}",
                    "yav": "0",
                    "afe": "0.0000",
                }
            )

    # Escrita final no formato do Starlight
    # ai made
    with open(OUTPUT_BASE_FILE, "w") as f:
        f.write(str(len(base_elements)) + "            [N_base]\n")
        for e in base_elements:
            line = f"{e['file']:45s} {e['age']:14s} {e['z']:11s} {e['code']:15s} {e['mstar']:10s} {e['yav']:5s} {e['afe']:s}\n"
            f.write(line)
        f.write(
            "# spec-file                                   age [yr]       Z           code            Mstar      YAV?  a/Fe\n"
        )
        f.write("# Total MILES base generated with real Mstar values from IAC tables\n")

    print(
        f"Sucesso! Base total com {len(base_elements)} componentes gerada em {os.path.abspath(OUTPUT_BASE_FILE)}"
    )


main()
