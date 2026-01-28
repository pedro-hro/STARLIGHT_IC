import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.interpolate import interp1d

def get_data_from_fits(filepath):
    """
    Lê fluxo, header e erro do arquivo FITS e retorna lambda, fluxo, erro.
    """
    with fits.open(filepath) as file:
        flux = file[0].data
        header = file[0].header
        # Tenta ler erro da extensão 1, se não existir cria array de zeros ou nan
        try:
            error = file[1].data
        except IndexError:
            error = np.full_like(flux, np.nan)

        wcs = WCS(header)
        pixel_indices = np.arange(len(flux))
        lambdas = wcs.pixel_to_world_values(pixel_indices)
        
        # Ajuste de unidade
        lambdas = lambdas * 1e10

    return lambdas, flux, error

def processar_espectros(files_dict, lambda_min=3200, lambda_max=9100, step=1.0, output_name='NGC0104.in', err_sistematico=0.05):
    """
    Combina várias bandas num único arquivo de entrada.
    
    Args:
        files_dict [dict]: Dicionário {'Banda': 'caminho/arquivo.fits'}
        lambda_min [float]: Lambda inicial
        lambda_max [float]: Lambda final
        step [float]: Passo do grid
        output_name [str]: Nome do arquivo de saída
    """

    master_lambda = np.arange(lambda_min, lambda_max + 1, step)
    master_flux = np.full_like(master_lambda, np.nan)
    master_error = np.full_like(master_lambda, np.nan)

    for banda, path in files_dict.items():
        try:
            l_orig, f_orig, e_orig = get_data_from_fits(path)
        except FileNotFoundError:
            print(f"    [ERRO] Arquivo não encontrado: {path}")
            continue

        # Interpolação
        interp_flux = interp1d(l_orig, f_orig, bounds_error=False, fill_value=np.nan)
        interp_error = interp1d(l_orig, e_orig, bounds_error=False, fill_value=np.nan)

        f_resampled = interp_flux(master_lambda)
        e_resampled = interp_error(master_lambda)

        # Lógica de combinação (prioriza quem tem menor erro nas pontas)
        has_new_data = ~np.isnan(f_resampled)
        is_master_empty = np.isnan(master_flux)

        # 1. Preenche onde está vazio
        mask_fill = has_new_data & is_master_empty
        master_flux[mask_fill] = f_resampled[mask_fill]
        master_error[mask_fill] = e_resampled[mask_fill]

        # 2. Onde sobrepõe, usa o menor erro
        mask_overlap = has_new_data & (~is_master_empty)

        safe_master_err = np.nan_to_num(master_error, nan=999999)
        safe_new_err = np.nan_to_num(e_resampled, nan=999999)
        
        better_error = safe_new_err < safe_master_err
        mask_update = mask_overlap & better_error

        master_flux[mask_update] = f_resampled[mask_update]
        master_error[mask_update] = e_resampled[mask_update]

    flag = np.zeros_like(master_flux)
    bad_pixels = np.isnan(master_flux)

    flag[bad_pixels] = 2
    master_flux[bad_pixels] = 0.0
    master_error[bad_pixels] = 0.0

    # Adiciona erro
    master_error = np.sqrt(master_error**2 + (err_sistematico * master_flux)**2)

    # Substitui nan por zero para evitar problemas
    master_error = np.nan_to_num(master_error, nan=0.0)

    output_data = np.column_stack((master_lambda, master_flux, master_error, flag))
    
    # Salva (botei 16 casas porque é float64, verificar no header)
    np.savetxt(output_name, output_data, fmt=['%.1f', '%.16e', '%.16e', '%d'])
    print(f"Arquivo gerado: {output_name}")
    
    return master_lambda, master_flux, master_error