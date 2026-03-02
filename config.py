import os

# CONFIGURAÇÕES DA PIPELINE
PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = PIPELINE_DIR

STARLIGHT_EXE = os.path.join(PROJECT_ROOT, "Starlightv04", "StarlightChains_v04.exe")
INPUTS_DIR = os.path.join(PROJECT_ROOT, "inputs_waggs")

# PARÂMETROS RUN
STARLIGHT_PARAMS = {
    # Arquivos Essenciais
    "base": "Base.Miles.MaxCut",
    "base_dir": "BasesMiles/",
    "mask": "Masks.EmLines.SDSS.gm",
    "template_config": "StCv04.C99-A.config",
    # Parâmetros Físicos e Numéricos
    "seed": 112017,
    "lllow_SN": 4900,
    "llup_SN": 5100,
    "olsyn_ini": 3600,
    "olsyn_fin": 7200,
    "delta_lambda": 1.0,  # Passo de onda (odlsyn)
    "fscale_chi2": 1.0,
    # Controle e Setup
    "kine": "FXK",  # 'FIT' ou 'FXK'
    "is_err": 1,  # 1 = Yes, 0 = No
    "is_flag": 1,  # 1 = Yes, 0 = No
    "extinction_law": "CCM",
    "v0": 0.0,
    "vd": 150.0,
}

# Configuração de Paralelização
CHUNK_SIZE = 1
