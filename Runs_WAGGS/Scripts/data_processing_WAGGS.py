import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d


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


def processar_espectros(files_dict, lambda_min=3200, lambda_max=9100, step=1.0, output_name='NGC0104_v2_orig_err.in',
                        fwhm_target=3.0, err_sistematico=0.05):
    """
    Combina várias bandas num único arquivo de entrada considerando a resolução da base.

    Args:
        files_dict [dict]: Dicionário {'Banda': 'caminho/arquivo.fits'}
        lambda_min [float]: Lambda inicial
        lambda_max [float]: Lambda final
        step [float]: Passo do grid
        output_name [str]: Nome do arquivo de saída
        fwhm_target: 3 # FWHM alvo em Angstroms (BC03: 3A - Depois botar referência)
        err_sistematico
    """

    master_lambda = np.arange(lambda_min, lambda_max + 1, step)
    master_flux = np.full_like(master_lambda, np.nan)
    master_error = np.full_like(master_lambda, np.nan)

    sigma_angstroms = fwhm_target / 2.355

    for banda, path in files_dict.items():
        try:
            l_orig, f_orig, e_orig = get_data_from_fits(path)
        except FileNotFoundError:
            print(f"    [ERRO] Arquivo não encontrado: {path}")
            continue

        # Ajuste de Resolução
        delta_l = np.nanmedian(np.diff(l_orig))
        if delta_l <= 0 or np.isnan(delta_l): delta_l = 0.4

        sigma_pix = sigma_angstroms / delta_l

        if sigma_pix > 0.5:
            f_orig = gaussian_filter1d(f_orig, sigma_pix)
            e_orig = gaussian_filter1d(e_orig, sigma_pix)  # O erro original também é suavizado

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
            with np.errstate(divide='ignore', invalid='ignore'):
                ratios = master_flux[mask_overlap] / f_resampled[mask_overlap]
                ratio = np.nanmedian(ratios)

            if not np.isnan(ratio) and ratio > 0:
                f_resampled *= ratio
                e_resampled *= ratio

        # 1. Preenche onde está vazio
        mask_fill = has_new_data & is_master_empty
        master_flux[mask_fill] = f_resampled[mask_fill]
        master_error[mask_fill] = e_resampled[mask_fill]

        # 2. Atualiza onde ambos têm dados, priorizando menor erro
        safe_master_err = np.nan_to_num(master_error, nan=999999)
        safe_new_err = np.nan_to_num(e_resampled, nan=999999)

        better_error = safe_new_err < safe_master_err
        mask_update = mask_overlap & better_error

        master_flux[mask_update] = f_resampled[mask_update]
        master_error[mask_update] = e_resampled[mask_update]

    # Flags
    flag = np.zeros_like(master_flux)
    bad_pixels = np.isnan(master_flux) | (master_flux <= 0) | np.isnan(master_error) | (master_error <= 0)

    flag[bad_pixels] = 2
    master_flux[bad_pixels] = 0.0
    master_error[bad_pixels] = 1e19
    master_error = np.nan_to_num(master_error, nan=1e19)

    # Adiciona erro sistemático
    good_pixels = ~bad_pixels
    master_error[good_pixels] = np.sqrt(
        master_error[good_pixels] ** 2 + (err_sistematico * master_flux[good_pixels]) ** 2)

    output_data = np.column_stack((master_lambda, master_flux, master_error, flag))

    np.savetxt(output_name, output_data, fmt=['%.1f', '%.16e', '%.16e', '%d'])
    print(f"Arquivo gerado: {output_name}")

    return master_lambda, master_flux, master_error
