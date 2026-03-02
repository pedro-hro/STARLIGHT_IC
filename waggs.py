import glob
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d

err_s = 0.0245
fwhm_target = 2.51
lambda_min = 3200
lambda_max = 9100
step = 1.0
R_dados = 6800
plots = True

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = SCRIPT_DIR
WAGGS_DIR = os.path.join(PROJECT_ROOT, "Dados", "WAGGS-DR1")

OUTPUT_IN_DIR = os.path.join(PROJECT_ROOT, "inputs_waggs")
OUTPUT_PLOT_DIR = os.path.join(PROJECT_ROOT, "input_plots_waggs")

os.makedirs(OUTPUT_IN_DIR, exist_ok=True)
os.makedirs(OUTPUT_PLOT_DIR, exist_ok=True)


def get_data_from_fits(filepath):
    """
    Lê fluxo, header e erro do arquivo FITS e retorna lambda, fluxo, erro.
    """
    with fits.open(filepath) as file:
        flux = file[0].data
        header = file[0].header
        try:
            error = file[1].data
        except IndexError:
            error = np.full_like(flux, np.nan)

        wcs = WCS(header)
        pixel_indices = np.arange(len(flux))
        lambdas = wcs.pixel_to_world_values(pixel_indices)
        lambdas = lambdas * 1e10

    return lambdas, flux, error


def processar_espectros(
    files_dict,
    lambda_min,
    lambda_max,
    step,
    output_name,
    fwhm_target,
    err_sistematico,
    R_dados,  # Resolução espectral dos dados de entrada (WAGGS: 6800),
):
    """
    Combina várias bandas num único arquivo de entrada considerando a resolução da base.
    Recebe um dicionário com as bandas e caminhos dos arquivos FITS, processa cada banda para igualar a resolução (se necessário), interpola para um grid comum, e combina os resultados priorizando os dados com menor erro.

    Args:
        files_dict [dict]: Dicionário {'Banda': 'caminho/arquivo.fits'}
        lambda_min [float]: Lambda inicial
        lambda_max [float]: Lambda final
        step [float]: Passo do grid
        output_name [str]: Nome do arquivo de saída
        fwhm_target [float]: FWHM alvo em Angstroms (ex: BC03=3.0, MILES=2.51)
        err_sistematico [float]: Fração de erro sistemático adicionado em quadratura
        R_dados [float]: Resolução espectral (R = λ/Δλ) dos dados de entrada (WAGGS: 6800)
    """

    master_lambda = np.arange(lambda_min, lambda_max + 1, step)
    master_flux = np.full_like(master_lambda, np.nan)
    master_error = np.full_like(master_lambda, np.nan)

    for band, path in files_dict.items():
        try:
            l_orig, f_orig, e_orig = get_data_from_fits(path)
        except FileNotFoundError:
            print(f"    [ERRO] Arquivo não encontrado: {path}")
            continue

        l_central = np.nanmedian(l_orig)
        fwhm_dados = l_central / R_dados
        fwhm_kernel_sq = fwhm_target**2 - fwhm_dados**2

        delta_l = np.nanmedian(np.diff(l_orig))

        if fwhm_kernel_sq > 0:
            sigma_angstroms = np.sqrt(fwhm_kernel_sq) / 2.355
            sigma_pix = sigma_angstroms / delta_l

            if sigma_pix > 0.5:
                f_orig = gaussian_filter1d(f_orig, sigma_pix)
                e_orig = np.sqrt(gaussian_filter1d(e_orig**2, sigma_pix))

        # Interpolação
        interp_flux = interp1d(l_orig, f_orig, bounds_error=False, fill_value=np.nan)
        interp_error = interp1d(l_orig, e_orig, bounds_error=False, fill_value=np.nan)

        f_resampled = interp_flux(master_lambda)
        e_resampled = interp_error(master_lambda)

        # Overlap e shift
        has_new_data = ~np.isnan(f_resampled)
        is_master_empty = np.isnan(master_flux)
        mask_overlap = has_new_data & (~is_master_empty)

        if np.any(mask_overlap):
            with np.errstate(divide="ignore", invalid="ignore"):
                ratios = master_flux[mask_overlap] / f_resampled[mask_overlap]
                ratio = np.nanmedian(ratios)

            if not np.isnan(ratio) and ratio > 0:
                f_resampled *= ratio
                e_resampled *= ratio

        mask_fill = has_new_data & is_master_empty
        master_flux[mask_fill] = f_resampled[mask_fill]
        master_error[mask_fill] = e_resampled[mask_fill]

        safe_master_err = np.nan_to_num(master_error, nan=999999)
        safe_new_err = np.nan_to_num(e_resampled, nan=999999)

        better_error = safe_new_err < safe_master_err
        mask_update = mask_overlap & better_error

        master_flux[mask_update] = f_resampled[mask_update]
        master_error[mask_update] = e_resampled[mask_update]

    # Flags
    flag = np.zeros_like(master_flux)
    bad_pixels = (
        np.isnan(master_flux)
        | (master_flux <= 0)
        | np.isnan(master_error)
        | (master_error <= 0)
    )

    flag[bad_pixels] = 2
    master_flux[bad_pixels] = 0.0
    master_error[bad_pixels] = 1e19
    master_error = np.nan_to_num(master_error, nan=1e19)

    # Erro sistemático
    good_pixels = ~bad_pixels
    master_error[good_pixels] = np.sqrt(
        master_error[good_pixels] ** 2
        + (err_sistematico * master_flux[good_pixels]) ** 2
    )

    output_data = np.column_stack((master_lambda, master_flux, master_error, flag))

    np.savetxt(output_name, output_data, fmt=["%.1f", "%.16e", "%.16e", "%d"])

    return master_lambda, master_flux, master_error


print(f"Procurando arquivos em: {os.path.abspath(WAGGS_DIR)}")
all_files = glob.glob(os.path.join(WAGGS_DIR, "norm_*.fits"))
print(f"Encontrados {len(all_files)} arquivos FITS brutos.")

targets = {}
for file_path in all_files:
    filename = os.path.basename(file_path)
    parts = filename.split("_")
    target_name = parts[1]
    band = parts[2][0]
    if target_name not in targets:
        targets[target_name] = {}
        print(f"  Alvo: {target_name}")
    targets[target_name][band] = file_path

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for target_name, bands_dict in targets.items():
        output_file = os.path.join(OUTPUT_IN_DIR, f"{target_name}.in")
        plot_file = os.path.join(OUTPUT_PLOT_DIR, f"{target_name}_processed.png")

        try:
            l_ambda, fluxo, erro = processar_espectros(
                bands_dict,
                lambda_min=lambda_min,
                lambda_max=lambda_max,
                step=step,
                output_name=output_file,
                err_sistematico=err_s,
                fwhm_target=fwhm_target,
                R_dados=R_dados,
            )

            mask_valid = erro < 70
            if plots:
                plt.figure(figsize=(12, 5))
                plt.plot(l_ambda, fluxo, color="black", lw=0.5, label="Fluxo")
                plt.fill_between(
                    l_ambda[mask_valid],
                    (fluxo - erro)[mask_valid],
                    (fluxo + erro)[mask_valid],
                    color="gray",
                    alpha=0.3,
                    label="Erro",
                )
                plt.title(f"Dados WAGGS: {target_name}")
                plt.xlabel(r"Comprimento de Onda $[\AA]$")
                plt.ylabel("Fluxo Normalizado")
                plt.xlim(lambda_min, lambda_max)
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                plt.savefig(plot_file, dpi=150)
                plt.close()
        except Exception as e:
            print(f"  [ERRO] Falha: {target_name}: {e}")
