import numpy as np
import matplotlib.pyplot as plt
import os

class StarlightOutput:
    """
    Classe para ler arquivos de output do STARLIGHT.
    """
    def __init__(this, filepath):
        this.filepath = filepath
        this.filename = os.path.basename(filepath)

        # Extraído do cabeçalho
        this.chi2 = None   # Qui^2
        this.adev = None   # Desvio
        this.av = None     # Extinção (A_V) [magnitudes]
        this.v0 = None     # V [km/s]
        this.vd = None     # Dispersão de velocidade [km/s]
        this.n0 = None     # Número de pontos usados no ajuste
        this.nl = None     # Número de pontos efetivos
        this.nclip = None  # Número de pontos clipados
        this.base = None  # Base usada (não implementado)

        # Tabela principal/população (começa com "# j x_j[%] Mini_j[%] ...")
        this.population = {
            'j': [],     # Índice do componente da base
            'x_j': [],   # Fração de LUZ [%]
            'm_ini': [], # Fração de MASSA INICIAL [%]
            'm_cor': [], # Fração de MASSA ATUAL [%]
            'age': [],   # Idade [anos]
            'Z': []      # Metalicidade
        }

        # Tabela final "Synthetic spectrum"
        this.spectrum = {
            'l_obs': [], # Comprimento de onda (Angstrom)
            'f_obs': [], # Fluxo Observado
            'f_syn': [], # Fluxo Modelo
            'wei': []    # Peso: >0 usado, <=0 mascarado/clipado
        }

        this.min_lambda = None
        this.max_lambda = None

        # Executa a leitura imediata ao instanciar a classe
        this.read_file()

    def read_file(this):
        """
        Lê o arquivo .out linha por linha e importa os atributos da classe.
        Usa flags de estado (reading_population, reading_spectrum) para navegar pelas seções.
        """
        if not os.path.exists(this.filepath):
            raise FileNotFoundError(f"Arquivo não encontrado: {this.filepath}")

        with open(this.filepath, 'r') as f:
            lines = f.readlines()

        reading_population = False
        reading_spectrum = False

        for i, line in enumerate(lines):
            line = line.strip()

            if '[chi2/Nl_eff]' in line:
                try: this.chi2 = float(line.split()[0])
                except: this.chi2 = np.nan
            elif '[adev (%)]' in line:
                try: this.adev = float(line.split()[0])
                except: this.adev = np.nan
            elif '[AV_min' in line:
                try: this.av = float(line.split()[0])
                except: this.av = np.nan
            elif '[v0_min' in line:
                try: this.v0 = float(line.split()[0])
                except: this.v0 = np.nan
            elif '[vd_min' in line:
                try: this.vd = float(line.split()[0])
                except: this.vd = np.nan
            elif '[NOl_eff]' in line:
                try: this.n0 = int(line.split()[0])
                except: this.n0 = np.nan
            elif '[Nl_eff]' in line:
                try: this.nl = int(line.split()[0])
                except: this.nl = np.nan
            elif '[Ntot_cliped & clip_method]' in line:
                try: this.nclip = int(line.split()[0])
                except: this.nclip = np.nan
            elif '[arq_base]' in line:
                try: this.base = line.split()[0]
                except: this.base = None

            # Identifica o cabeçalho: "# j x_j(%)..."
            if line.startswith('# j') and 'x_j(%)' in line:
                reading_population = True
                continue
            
            if reading_population:
                # O fim da tabela é marcado por uma linha vazia ou início de nova seção com "##"
                if len(line) == 0 or line.startswith('##'):
                    reading_population = False
                else:
                    cols = line.split()
                    try:
                        this.population['j'].append(int(cols[0]))
                        this.population['x_j'].append(float(cols[1])) # Fração de luz [%]
                        this.population['m_ini'].append(float(cols[2])) # Fração de massa inicial [%]
                        this.population['m_cor'].append(float(cols[3])) # Fração de massa atual [%]
                        this.population['age'].append(float(cols[4])) # [anos]
                        this.population['Z'].append(float(cols[5]))
                    except (ValueError, IndexError):
                        pass 

            if '## Synthetic spectrum (Best Model) ##' in line:
                reading_spectrum = True
                continue 

            if reading_spectrum:
                cols = line.split()
                # pula o número de pontos na primeira linha
                if len(cols) == 1: 
                    continue
                
                #lambda obs_flux syn_flux weight
                if len(cols) >= 3:
                    try:
                        this.spectrum['l_obs'].append(float(cols[0]))
                        this.spectrum['f_obs'].append(float(cols[1]))
                        this.spectrum['f_syn'].append(float(cols[2]))
                        this.spectrum['wei'].append(float(cols[3]))
                    except ValueError:
                        pass

        # Converte listas python para numpy
        for key in this.population:
            this.population[key] = np.array(this.population[key])
        for key in this.spectrum:
            this.spectrum[key] = np.array(this.spectrum[key])

        this.min_lambda = np.min(this.spectrum['l_obs']) if len(this.spectrum['l_obs']) > 0 else None
        this.max_lambda = np.max(this.spectrum['l_obs']) if len(this.spectrum['l_obs']) > 0 else None

    def calculate_mean_properties(this):
        """
        Calcula as propriedades médias da população (Idade e Metalicidade).
        Realiza médias ponderadas tanto por LUZ (x_j) quanto por MASSA (Mini_j).
        
        Returns:
            dict: Dicionário com médias logarítmicas e lineares.
        """
        x = this.population['x_j']
        m = this.population['m_ini']
        age = this.population['age']
        Z = this.population['Z']

        sum_x = np.sum(x)
        sum_m = np.sum(m)

        # Não sei se precisa:
        if sum_x == 0 or sum_m == 0:
            return None

        # Média Ponderada por LUZ
        mean_log_age_light = np.sum(x * np.log10(age)) / sum_x
        mean_Z_light = np.sum(x * Z) / sum_x

        # Média Ponderada por MASSA
        mean_log_age_mass = np.sum(m * np.log10(age)) / sum_m
        mean_Z_mass = np.sum(m * Z) / sum_m

        return {
            'mean_log_age_light': mean_log_age_light,
            'mean_age_light_gyr': (10**mean_log_age_light)/1e9, # Converte log(anos) para Giga-anos
            'mean_Z_light': mean_Z_light,
            
            'mean_log_age_mass': mean_log_age_mass,
            'mean_age_mass_gyr': (10**mean_log_age_mass)/1e9,
            'mean_Z_mass': mean_Z_mass,
        }

    def plot_fit(this):
        """
        Plota o ajuste espectral (Observado vs Modelo com as bases)
        Destaca pontos clipados e qualidade (Chi2, Adev).
        """
        l = this.spectrum['l_obs']
        fo = this.spectrum['f_obs']
        fs = this.spectrum['f_syn']
        wei = this.spectrum['wei']

        plt.figure(figsize=(12, 6))

        plt.plot(l, fo, 'k', label='Observado', lw=0.5)
        plt.plot(l, fs, 'r', label='Modelo', lw=1.0, alpha=0.8)
        
        # plota pontos ignorados
        masked = wei <= 0
        if np.any(masked):
            plt.scatter(l[masked], fo[masked], color='blue', marker='x', s=15, label='Clipado', zorder=5)

        plt.xlabel('Comprimento de Onda ($\AA$)')
        plt.ylabel('Fluxo Normalizado')

        title_str = f'Ajuste: {this.filename}\n'
        if this.chi2 is not None: title_str += f'$\chi^2/N_{{eff}}={this.chi2:.2f}$  '
        if this.adev is not None: title_str += f'Adev={this.adev:.2f}%  '
        if this.av is not None:   title_str += f'$A_V={this.av:.2f}$ mag'
            
        plt.title(title_str)
        plt.legend(frameon=True)
        plt.grid(True, alpha=0.3)
        plt.savefig(this.filename.replace('.out', '_fit.png'))
        plt.show()