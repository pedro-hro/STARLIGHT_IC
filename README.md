# Análise de Populações Estelares: NGC 104 (47 Tucanae)

Este repositório contém os scripts de processamento, configurações de grid e notebooks de análise desenvolvidos na minha Iniciação Científica. Inicialmente para o estudo do Aglomerado Globular **NGC 104 (47 Tuc)**.

---

## 📂 Estrutura do Repositório

- **`STARLIGHT/Runs_WAGGS/`**:
  - **`NGC0104/`**: Contém os arquivos específicos do alvo.
    - `grid_NGC0104.in`: Arquivo de definição do *grid* para o STARLIGHT.
    - `StCv04.C99.config`: Arquivo de configuração (parâmetros de MCMC, *clipping*, etc.).
    - `ngc0104.ipynb`: Jupyter Notebook principal com a análise dos resultados (`.out`) e geração de plots.
  - **`Scripts/`**: Rotinas auxiliares.
    - `data_processing_WAGGS.py`: Processamento dos fits brutos, fusão de bandas e tratamento de erros.
    - `starlight_analysis_WAGGS.py`: Classes para leitura e interpretação das saídas do STARLIGHT.

> **Nota:** Os dados brutos (`.fits`), executáveis do STARLIGHT e arquivos de base (`Base.BC03.N`) **não** estão incluídos neste repositório para manter a leveza e respeitar direitos de uso.

---

## 🔭 Física

### 1. Dados (WAGGS)
Espectros integrados do *WiFeS Atlas of Galactic Globular Clusters* (WAGGS DR1).
- **Cobertura:** 3270 Å a 9050 Å.
- **Pré-processamento:** Os dados foram reamostrados para $\Delta\lambda = 1.0$ Å e receberam um erro sistemático de 5% para estabilidade numérica.

### 2. Configuração do STARLIGHT
*  **Versão do STARLIGHT:** v04 (Cid Fernandes 2007).
*   **Base Estelar:** `Base.BC03.N`.
    *   45 SSPs.
    *   3 Metalicidades ($Z = 0.004, 0.02, 0.05$) e 15 Idades.
*   **Configuração:** `StCv04.C99.config`.
*   **Extinção:** Lei de avermelhamento CCM.

---

## 📊 Resultados Preliminares

| Parâmetro | Valor Obtido 
| :--- | :--- |
| **Idade (Luz)** | ~4.9 Gyr | 
| **Idade (Massa)** | ~8.4 Gyr | 
| **Metalicidade** | $Z \approx 0.019$  
| **Extinção($A_V$)** | ~0.0 mag 

---

## ⚙️ Instalação e Reprodução

### 1. Pré-requisitos
*   **Python 3.8+**
*   Bibliotecas: `numpy`, `matplotlib`, `astropy`, `scipy`.

### 2. Configurando o STARLIGHT
1.  Baixe o **STARLIGHT v04** no [site oficial](http://www.starlight.ufsc.br).
2. Baixe os dados **WAGGS** do [repositório oficial](https://researchdata.edu.au/wiggs-wifes-atlas-galactic-globular-clusters/165145).
3.  Coloque o executável na pasta `STARLIGHT/Runs_WAGGS/NGC0104/` (ou aponte o caminho no *grid*).
4.  Certifique-se de ter a pasta `BasesDir/` contendo os arquivos da base.

### 3. Execução
```bash
cd Runs_WAGGS/NGC0104/
./StarlightChains_v04.exe < grid_NGC0104.in
```
