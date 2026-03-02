import os
import re

import astropy.io.fits as fits
import numpy as np
from astropy.wcs import WCS

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = SCRIPT_DIR

FITS_PATH = os.path.join(PROJECT_ROOT, "Starlightv04", "MILES_PADOVA00_KU_baseFe")
SPEC_DIR = os.path.join(PROJECT_ROOT, "Starlightv04", "BasesMiles")
MASS_FILE = os.path.join(
    PROJECT_ROOT, "Starlightv04", "MILES_PADOVA00_KU_baseFe", "out_mass_KU_PADOVA00"
)
OUTPUT_BASE_DIR = os.path.join(PROJECT_ROOT, "Starlightv04", "BaseFiles")

Z_SUN = 0.019  # Referência ???

pattern = re.compile(r"Mku1.30Z([mp])(\d+\.\d+)T(\d+\.\d+)")


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


def convert_fits_to_spec(fits_path, spec_dir):
    print("Verificando/Convertendo FITS para SPEC")
    fits_files = []
    all_dotfits = os.listdir(FITS_PATH)
    for f in all_dotfits:
        if f.endswith(".fits"):
            fits_files.append(f)
    converted_count = 0

    for filename in fits_files:
        filepath = os.path.join(fits_path, filename)
        output_filename = filename.replace(".fits", ".spec")
        output_filepath = os.path.join(spec_dir, output_filename)

        if not os.path.exists(output_filepath):
            lambdas, flux = get_data_from_fits(filepath)
            np.savetxt(output_filepath, np.column_stack((lambdas, flux)))
            converted_count += 1

    if converted_count > 0:
        print(
            f"Processamento concluído: {converted_count} novos arquivos .spec gerados."
        )
    else:
        print("Todos os arquivos .spec já existem. Pulando conversão.")

    return True


def select_ages_and_metallicities():
    all_specs = os.listdir(SPEC_DIR)
    ages_list = set()
    mh_list = set()
    for f in all_specs:
        if f.endswith(".spec"):
            match = pattern.search(f)
            if match:
                sign, z_val, age_gyr_str = match.groups()
                mh_val = float(z_val) * (-1 if sign == "m" else 1)
                age_gyr = float(age_gyr_str)
                ages_list.add(age_gyr)
                mh_list.add(mh_val)
    ages_list = sorted(list(ages_list))
    mh_list = sorted(list(mh_list))

    print("\nIdades disponíveis (Gyr):")
    for i, age in enumerate(ages_list):
        print(f"{i}: {age:.4f} Gyr")

    selected_ages = input(
        "\nDigite os índices das idades desejadas (ex: 0,2,4) ou 'all' para todas: "
    )

    if selected_ages.strip().lower() == "all":
        selected_ages_gyr = ages_list
    else:
        selected_ages_gyr = []
        for i in selected_ages.split(","):
            stripped = i.strip()
            if stripped.isdigit() and int(stripped) < len(ages_list):
                selection = ages_list[int(stripped)]
                selected_ages_gyr.append(selection)

    print("\nMetallicidades disponíveis ([M/H]):")
    for i, mh in enumerate(mh_list):
        print(f"{i}: {mh:.4f} [M/H]")

    selected_mh = input(
        "Digite os índices das metallicidades desejadas (ex: 0,1) ou 'all' para todas: "
    )

    if selected_mh.strip().lower() == "all":
        selected_mh_list = mh_list
    else:
        selected_mh_list = []
        for i in selected_mh.split(","):
            stripped = i.strip()
            if stripped.isdigit() and int(stripped) < len(mh_list):
                selection = mh_list[int(stripped)]
                selected_mh_list.append(selection)

    return selected_ages_gyr, selected_mh_list


def generate_filtered_base(spec_dir, mass_file):
    output_base_name = input(
        f"\nDigite o nome do arquivo de saída para a base filtrada (ex: Base.Miles.X): "
    )

    output_base_file = os.path.join(OUTPUT_BASE_DIR, output_base_name)

    selected_ages_list, selected_mh_list = select_ages_and_metallicities()
    mass_map = load_mass_map(mass_file)

    spec_files = []
    for f in os.listdir(spec_dir):
        spec_files.append(f)

    base_elements = []

    for filename in spec_files:
        match = pattern.search(filename)
        if match:
            sign, z_val, age_gyr_str = match.groups()
            mh_val = float(z_val) * (-1 if sign == "m" else 1)
            age_gyr = float(age_gyr_str)

            # Filtro
            if mh_val not in selected_mh_list:
                continue

            if age_gyr not in selected_ages_list:
                continue

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

    with open(output_base_file, "w") as f:
        f.write(str(len(base_elements)) + "            [N_base]\n")
        for e in base_elements:
            line = f"{e['file']:45s} {e['age']:14s} {e['z']:11s} {e['code']:15s} {e['mstar']:10s} {e['yav']:5s} {e['afe']:s}\n"
            f.write(line)
        f.write(
            "# spec-file                                   age [yr]       Z           code            Mstar      YAV?  a/Fe\n"
        )
        f.write(
            "# Base MILES filtrada gerada no arquivo Base.Miles.Total. Idades e metallicidades selecionadas:\n"
        )
    return output_base_name


def run():
    print("\n--- Gerador de Base MILES ---")
    convert_fits_to_spec(FITS_PATH, SPEC_DIR)
    base_name = generate_filtered_base(SPEC_DIR, MASS_FILE)
    
    # Criar link simbólico para a pasta atual (STARLIGHT raiz)
    link_dir = PROJECT_ROOT
    os.makedirs(link_dir, exist_ok=True)
    
    base_original_path = os.path.join(OUTPUT_BASE_DIR, base_name)
    link_path = os.path.join(link_dir, base_name)
    
    if os.path.exists(link_path) or os.path.islink(link_path):
        os.remove(link_path)
    os.symlink(base_original_path, link_path)
    print(f"  [LINK CRIADO] {link_path} -> {base_original_path}")
    
    return base_name


if __name__ == "__main__":
    run()
