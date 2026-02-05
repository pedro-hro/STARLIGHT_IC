import numpy as np
import os
import re

def get_base_info(master_base):
    """
    Retorna as idades (Gyr) e metalicidades ([M/H]) únicas presentes na base mestre.
    """
    if not os.path.exists(master_base):
        print(f"Erro: Base {master_base} não encontrada.")
        return [], []
    
    with open(master_base, 'r') as f:
        lines = f.readlines()
    
    pattern = re.compile(r"Mku1.30Z([mp])(\d+\.\d+)T(\d+\.\d+)")
    ages = set()
    mhs = set()
    
    for line in lines[1:]:
        if line.startswith('#') or not line.strip(): continue
        filename = line.split()[0]
        match = pattern.search(filename)
        if match:
            sign, z_val, age_gyr_str = match.groups()
            mh = float(z_val) * (-1 if sign == 'm' else 1)
            ages.add(float(age_gyr_str))
            mhs.add(mh)
            
    return sorted(list(ages)), sorted(list(mhs))

def create_custom_base(master_base, output_name, ages_gyr=None, mh_list=None):
    """
    Cria uma sub-base filtrada extraindo infos do nome do arquivo (evita Z_SUN).
    """
    if not os.path.exists(master_base):
        print(f"Erro: Base mestre não encontrada.")
        return

    with open(master_base, 'r') as f:
        lines = f.readlines()

    data_lines = [l for l in lines[1:] if not l.startswith('#') and l.strip()]
    header_comments = [l for l in lines if l.startswith('#')]
    pattern = re.compile(r"Mku1.30Z([mp])(\d+\.\d+)T(\d+\.\d+)")

    selected_lines = []
    for line in data_lines:
        filename = line.split()[0]
        match = pattern.search(filename)
        if not match: continue
        
        sign, z_val, age_gyr_str = match.groups()
        mh_val = float(z_val) * (-1 if sign == 'm' else 1)
        age_val = float(age_gyr_str)

        match_mh = True
        if mh_list is not None:
            match_mh = any(np.isclose(mh_val, target, atol=0.01) for target in mh_list)

        match_age = True
        if ages_gyr is not None:
            match_age = any(np.isclose(age_val, target, rtol=0.01) for target in ages_gyr)

        if match_mh and match_age:
            selected_lines.append(line)

    with open(output_name, 'w') as f:
        f.write(f"{len(selected_lines)}            [N_base]\n")
        f.writelines(selected_lines)
        f.writelines(header_comments)

    print(f"Sucesso: Base '{output_name}' gerada com {len(selected_lines)} componentes.")

def show_selection_ui(master_base, output_name):
    """
    Gera widgets de seleção para serem usados no Jupyter Notebook.
    """
    import ipywidgets as widgets
    from IPython.display import display

    ages, mhs = get_base_info(master_base)
    
    age_sel = widgets.SelectMultiple(options=[(f"{a:.4f} Gyr", a) for a in ages], description='Idades', rows=10)
    mh_sel = widgets.SelectMultiple(options=[(f"{m:+.2f}", m) for m in mhs], description='[M/H]', rows=7)
    button = widgets.Button(description="Gerar Base Customizada", button_style='success')

    def on_button_clicked(b):
        create_custom_base(master_base, output_name, list(age_sel.value), list(mh_sel.value))

    button.on_click(on_button_clicked)
    display(widgets.HBox([age_sel, mh_sel]), button)