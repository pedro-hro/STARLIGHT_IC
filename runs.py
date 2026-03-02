import glob
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import config
import pandas as pd

base_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(base_dir, "Scripts"))

import starlight_output_analysis as sl_analysis


def select_targets(input_dir):
    all_files = sorted(glob.glob(os.path.join(input_dir, "*.in")))

    if not all_files:
        print(
            "Nenhum arquivo .in encontrado! Rode os scripts de pré-processamento primeiro."
        )
        return []

    targets = [os.path.basename(f).replace(".in", "") for f in all_files]

    print("\nAlvos disponíveis:")
    for i, target in enumerate(targets):
        print(f"{i}: {target}")

    selected_indices = input('\nDigite os índices ou "all" para rodar todos: ')

    if selected_indices.strip().lower() == "all":
        return targets
    else:
        indices = [
            int(idx.strip())
            for idx in selected_indices.split(",")
            if idx.strip().isdigit()
        ]
        selected_targets = [targets[idx] for idx in indices if idx < len(targets)]
        return selected_targets


def main():

    work_dirs = ["grids", "outputs", "logs"]
    for d in work_dirs:
        full_path = os.path.join(config.PIPELINE_DIR, d)
        if os.path.exists(full_path):
            shutil.rmtree(full_path)
        os.makedirs(full_path)

    selected_targets = select_targets(config.INPUTS_DIR)

    s = config.STARLIGHT_PARAMS
    base_name = s["base"]

    chunks = [
        selected_targets[i : i + config.CHUNK_SIZE]
        for i in range(0, len(selected_targets), config.CHUNK_SIZE)
    ]

    rel_inputs_dir = os.path.relpath(config.INPUTS_DIR, config.PIPELINE_DIR)

    if not rel_inputs_dir.endswith(os.sep):
        rel_inputs_dir += os.sep

    for i, chunk in enumerate(chunks):
        grid_filename = os.path.join(config.PIPELINE_DIR, f"grids/grid_{i + 1}.in")
        n = len(chunk)

        header = f"""{n}    [Number of fits to run]
{s["base_dir"]}    [base_dir]
{rel_inputs_dir}    [obs_dir]
./             [mask_dir]
outputs/             [out_dir]
{s["seed"]}      [random seed]
{s["lllow_SN"]}         [llow_SN] lower-lambda of S/N window
{s["llup_SN"]}         [lupp_SN] upper-lambda of S/N window
{s["olsyn_ini"]}         [Olsyn_ini] lower-lambda for fit
{s["olsyn_fin"]}         [Olsyn_fin] upper-lambda for fit
{s["delta_lambda"]}            [Odlsyn] delta-lambda for fit
{s["fscale_chi2"]}            [fscale_chi2] fudge-factor for chi2
{s["kine"]}            [FIT/FXK] Fit or Fix kinematics
{s["is_err"]}              [IsErrSpecAvailable] 1/0 = Yes/No
{s["is_flag"]}              [IsFlagSpecAvailable] 1/0 = Yes/No
"""
        with open(grid_filename, "w") as f:
            f.write(header)
            for infile in chunk:
                # Formato do grid: spectro.in config_file base_file mask extinction v0 vd spectro.out
                line = f"{infile}.in   {s['template_config']}   {base_name}   {s['mask']}   {s['extinction_law']}   {s['v0']}   {s['vd']}   {infile}.out\n"
                f.write(line)

    # Paralelização: roda os grids gerados usando subprocess e ThreadPoolExecutor
    grids = sorted(glob.glob(os.path.join(config.PIPELINE_DIR, "grids/grid_*.in")))

    def executar_starlight(i, grid_path):
        log_file = os.path.join(config.PIPELINE_DIR, "logs", f"grid_{i + 1}.log")
        cmd = f"{config.STARLIGHT_EXE} < {grid_path} > {log_file} 2>&1"
        subprocess.run(cmd, shell=True, cwd=config.PIPELINE_DIR)
        print(f"  [CONCLUÍDO] {os.path.basename(grid_path)}")

    total_cores = os.cpu_count() or 1
    max_workers = int(total_cores * (5 / 6)) or 1

    print(
        f"  > Hardware detectado: {total_cores} núcleos. Alocando {max_workers} threads simultâneas.\n"
    )

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for i, grid_path in enumerate(grids):
            executor.submit(executar_starlight, i, grid_path)

    print(f"   Outputs salvos em: {os.path.join(config.PIPELINE_DIR, 'outputs/')}")

    # Fazer analise do output
    for res in sorted(glob.glob(os.path.join(config.PIPELINE_DIR, "outputs/*.out"))):
        starlight_output = sl_analysis.StarlightOutput(res)
        props = starlight_output.calculate_mean_properties()

        df = pd.DataFrame(
            {
                "Target": os.path.basename(res).replace(".out", ""),
                "Mean Age(by light)": props["mean_age_light_gyr"],
                "Mean Z (by light)": props["mean_Z_light"],
                "Mean Age (by mass)": props["mean_age_mass_gyr"],
                "Mean Z (by mass)": props["mean_Z_mass"],
                "A_V": starlight_output.av,
                "chi2": starlight_output.chi2,
                "adev": starlight_output.adev,
                "Clip %": (starlight_output.nclip / starlight_output.n0) * 100,
            },
            index=[0],
        )

        df.to_csv(
            os.path.join(config.PIPELINE_DIR, "outputs", "summary.csv"),
            mode="a",
            header=not os.path.exists(
                os.path.join(config.PIPELINE_DIR, "outputs", "summary.csv")
            ),
            index=False,
        )


if __name__ == "__main__":
    main()
