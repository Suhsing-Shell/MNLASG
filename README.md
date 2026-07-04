# MNLASG

Este repositório contém a implementação em Python de um modelo numérico transiente de parâmetros concentrados para a simulação de sistemas de refrigeração termoelétrica baseados em módulos Peltier (TEC1-12706). O modelo descreve o comportamento térmico dinâmico do sistema considerando o acoplamento eletrotérmico, a dependência das propriedades do módulo com a temperatura, os mecanismos de transferência de calor e diferentes configurações de dissipação térmica. O algoritmo foi desenvolvido para apoiar o dimensionamento de refrigeradores portáteis, permitindo análises de sensibilidade, otimização paramétrica e avaliação do desempenho térmico sob diferentes condições ambientais. Os resultados obtidos foram comparados com dados experimentais disponíveis na literatura para verificar a capacidade preditiva do modelo.

- Modelagem transiente por parâmetros concentrados
- Simulação de módulos Peltier TEC1-12706
- Propriedades termoelétricas dependentes da temperatura
- Dissipação por ar e por watercooler
- Análise de sensibilidade via matriz jacobiana
- Simulações em Python

# Parâmentros

        # --- GEOMETRIA E ISOLAMENTO ---
        'dim_bloco': (0.19, 0.11, 0.14),           # Dimensões externas da câmara/bloco (Comprimento, Largura, Altura) em [m]
        'L_ins': 0.020,                             # Espessura da camada de isolamento térmico em [m] (Ex: 20 mm)
        'k_ins': 0.035,                             # Condutividade térmica do material isolante em [W/m·K] (Ex: Poliuretano/EPS)
        
        # --- COEFICIENTES DE CONVECÇÃO ---
        'h_int': 8.0,                               # Coeficiente de película de convecção interna em [W/m²·K]
        'h_ext': 22.0,                              # Coeficiente de película de convecção externa (ar ambiente) em [W/m²·K]
        
        # --- CONDIÇÕES AMBIENTAIS E CONFIGURAÇÃO ---
        'T_amb': 28.2,                              # Temperatura do ar ambiente externo em [°C]
        'n_TEC': 1,                                 # Quantidade de pastilhas Peltier ligadas no sistema
        
        # --- PROPRIEDADES NOMINAIS DA PELTIER (Na temperatura de referência) ---
        'alpha0': 0.050,                            # Coeficiente Seebeck nominal do módulo em [V/K]
        'R0': 2.0,                                  # Resistência elétrica nominal do módulo em [Ω]
        'K0': 0.8,                                  # Condutância térmica nominal do módulo em [W/K]
        
        # --- COEFICIENTES DE DEPENDÊNCIA TÉRMICA (Efeito da temperatura nas propriedades) ---
        'beta_alpha': -0.0015,                      # Coeficiente de variação do Coeficiente Seebeck com a temperatura [1/K]
        'beta_R': 0.004,                            # Coeficiente de variação da Resistência Elétrica com a temperatura [1/K]
        'beta_K': 0.0015,                           # Coeficiente de variação da Condutância Térmica com a temperatura [1/K]
        
        # --- SISTERMA DE RESFRIAMENTO DO LADO QUENTE ---
        'cooling_type': 'air',                      # Tipo de dissipação do lado quente ('air' para dissipador/ventoinha ou 'water' para waterblock)
        'R_hs': 0.32164,                            # Resistência térmica total do dissipador de calor em [K/W] (Calibrado com contato)
        
        # --- PARÂMETROS SE RESFRIAMENTO POR ÁGUA (Ignorados se cooling_type='air') ---
        'm_dot_w': 0.05,                            # Vazão mássica de água no bloco de resfriamento em [kg/s] # type: ignore
        'cp_w': 4180.0,                             # Calor específico da água em [J/kg·K]
        'T_w_in': 28.2,                             # Temperatura da água na entrada do waterblock em [°C]
        
        # --- PARÂMETROS DE CONTATO E SIMULAÇÃO DINÂMICA ---
        'R_block': 0.02,                            # Resistência térmica de contato entre a pastilha e o bloco de acoplamento em [K/W]
        'T_init': 28.2,                             # Temperatura inicial de todo o sistema no instante t=0 em [°C]
        
        # --- CRITÉRIOS DE CONVERGÊNCIA / ESTABILIZAÇÃO ---
        'stabilization_tolerance': 0.5,             # Janela de tolerância absoluta para considerar a temperatura estável em [°C]
        'stabilization_rate_tolerance': 0.01,        # Taxa de variação máxima (dT/dt) permitida para definir regime permanente em [K/s]
        'stabilization_max_time': 7200.0            # Tempo limite máximo de execução da simulação em [s] (Equivale a 2 horas)

    # --- PARÂMETROS DE ENTRADA DO ENSAIO ---
    I_test = 5.9                                    # Corrente elétrica contínua aplicada para testar a pastilha Peltier em [A]
