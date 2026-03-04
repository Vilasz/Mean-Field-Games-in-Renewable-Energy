# Mean Field Games em Sistemas de Energia Renovável

Este projeto investiga a aplicação da teoria de **Mean Field Games** a sistemas elétricos com grande participação de fontes renováveis — solar, eólica e hidrelétrica — tendo como caso de estudo o **Sistema Interligado Nacional (SIN)** do Brasil em 2025.

O SIN é a espinha dorsal do setor elétrico brasileiro, interconectando quatro subsistemas regionais (Norte, Nordeste, Sudeste e Sul) e coordenando o despacho de uma das matrizes energéticas mais diversificadas do mundo: hidrelétricas de grande porte, extensos parques eólicos no Nordeste, expansão acelerada de solar fotovoltaica, usinas nucleares e um parque térmico de backup. Entender como essas fontes interagem e como o custo marginal de operação responde às condições do sistema é o problema central que este repositório aborda.

## Motivação

O operador do sistema elétrico enfrenta diariamente o desafio de equilibrar oferta e demanda em tempo real, minimizando custos e garantindo segurança de suprimento. Fontes renováveis intermitentes (solar e eólica) trazem variabilidade que precisa ser compensada por fontes controláveis (hidrelétrica e térmica). Ao mesmo tempo, o custo marginal de operação (CMO) reflete a escassez ou abundância de recursos em cada momento e região, servindo como sinal de preço para o mercado.

Este projeto constrói um modelo de **despacho social ótimo** que, a partir de dados operacionais reais do ONS, determina o nível ideal de geração térmica/nuclear hora a hora, respeitando restrições de capacidade, rampa e reserva. Os resultados são validados contra o despacho efetivamente observado, e o CMO calculado pelo modelo DECOMP é analisado para revelar padrões sazonais e relações com o mix de geração.

---

## Formulação Matemática

### Variáveis e Notação

Seja $s \in \mathcal{S} = \{N, NE, SE, S\}$ o índice de subsistema e $t \in \{0, 1, \ldots, T-1\}$ o índice temporal (horas). Definimos:

| Símbolo | Descrição | Tipo |
|---------|-----------|------|
| $g^s_{t,s}$ | Geração solar | Exógena (observada) |
| $g^r_{t,s}$ | Geração eólica | Exógena (observada) |
| $g^h_{t,s}$ | Geração hidrelétrica | Exógena (observada) |
| $g^n_{t,s}$ | Geração controlável (térmica + nuclear) | **Variável de decisão** |
| $D_{t,s}$ | Demanda efetiva (carga medida) | Exógena |
| $x_{int,t,s}$ | Importação líquida via intercâmbio | Exógena |
| $c_{t,s}$ | Curtailment (desperdício renovável) | **Variável de decisão** |
| $u_{t,s}$ | Déficit (corte de carga) | **Variável de decisão** |
| $R_{t,s}$ | Reserva operacional alocada | **Variável de decisão** |

A demanda líquida vista pelos geradores locais incorpora os fluxos de interconexão:

$$
D_{net,t,s} = D_{t,s} - x_{int,t,s}
$$

onde $x_{int,t,s} > 0$ indica importação (reduz demanda local) e $x_{int,t,s} < 0$ indica exportação.

### Problema do Planejador Social

Para cada subsistema $s$, o planejador minimiza o custo operacional total:

$$
\min_{g^n, c, u, R} \; J_s = \sum_{t=0}^{T-1} \left[ c_1 g^n_{t,s} + \frac{c_2}{2}\left(g^n_{t,s}\right)^2 + \frac{\gamma}{2}\left(g^n_{t,s} - g^n_{t-1,s}\right)^2 + \pi_u u_{t,s} + \pi_c c_{t,s} + \kappa R_{t,s} \right]
$$

sujeito às restrições:

$$
\begin{aligned}
g^s_{t,s} + g^r_{t,s} + g^h_{t,s} + g^n_{t,s} &= D_{net,t,s} + c_{t,s} - u_{t,s} && \text{(balanço de potência)} \quad [\lambda_{t,s}] \\[4pt]
0 \leq g^n_{t,s} &\leq K_{n,s} && \text{(limites de capacidade)} \quad [\underline{\mu}_{t,s}, \overline{\mu}_{t,s}] \\[4pt]
\left|g^n_{t,s} - g^n_{t-1,s}\right| &\leq \rho_s && \text{(restrição de rampa)} \\[4pt]
R_{t,s} + sR_{t,s} &\geq z_\alpha \cdot \sigma_{\epsilon}(t,s) && \text{(reserva mínima)} \\[4pt]
g^n_{t,s} + R_{t,s} &\leq K_{n,s} && \text{(viabilidade de reserva)} \\[4pt]
c_{t,s}, \; u_{t,s}, \; R_{t,s}, \; sR_{t,s} &\geq 0 && \text{(não-negatividade)}
\end{aligned}
$$

onde os termos entre colchetes $[\cdot]$ indicam os multiplicadores duais (preços-sombra) associados.

### Interpretação dos Termos da Função Objetivo

- $c_1 g^n_{t,s}$: custo variável linear de geração térmica (combustível)
- $\frac{c_2}{2}(g^n_{t,s})^2$: componente quadrática que captura custos crescentes de despacho (heat-rate curves)
- $\frac{\gamma}{2}(g^n_{t,s} - g^n_{t-1,s})^2$: penalização de rampa que incentiva suavidade no despacho, refletindo custos de start-up/shutdown e estresse térmico
- $\pi_u u_{t,s}$: custo de déficit, calibrado no Value of Lost Load (VOLL), tipicamente $\pi_u \sim 10^4$ R\$/MWh
- $\pi_c c_{t,s}$: custo de curtailment (desperdício de renovável), tipicamente $\pi_c \sim 10$ R\$/MWh
- $\kappa R_{t,s}$: custo de manter reserva operacional disponível

A variável $sR_{t,s} \geq 0$ é um slack de reserva penalizado na função objetivo com peso $\Pi_{sR} = 10^6$ para garantir viabilidade mesmo quando $K_n$ é insuficiente para cobrir simultaneamente despacho e reserva.

### Calibração de Parâmetros a Partir de Dados

Os parâmetros $K_{n,s}$ e $\rho_s$ são estimados a partir da série histórica observada de geração controlável:

$$
K_{n,s} = Q_{0.995}\left(\{g^n_{obs,t,s}\}_{t=0}^{T-1}\right), \qquad
\rho_s = Q_{0.995}\left(\{|g^n_{obs,t,s} - g^n_{obs,t-1,s}|\}_{t=1}^{T-1}\right)
$$

O uso do percentil $99.5\%$ em vez do máximo garante robustez a outliers operacionais (emergências, erros de medição) que não representam a operação normal.

### Incerteza Renovável e Requisito de Reserva

A incerteza $\sigma_\epsilon(t,s)$ é estimada a partir do desvio-padrão do erro de previsão renovável, estratificado por mês $m$ e hora do dia $h$:

$$
\epsilon_{t,s} = (g^s_{t,s} + g^r_{t,s}) - \mu_{m(t), h(t), s}
$$

$$
\sigma_\epsilon(m, h, s) = \sqrt{\frac{1}{N_{m,h}-1} \sum_{t \in \mathcal{T}_{m,h}} \epsilon_{t,s}^2}
$$

onde $\mu_{m,h,s} = \mathbb{E}[g^s + g^r \mid \text{mês}=m, \text{hora}=h, \text{subsistema}=s]$ e $\mathcal{T}_{m,h}$ é o conjunto de horas com mês $m$ e hora $h$.

O requisito de reserva garante cobertura probabilística:

$$
R_{t,s} \geq z_\alpha \cdot \sigma_\epsilon(m(t), h(t), s), \qquad z_{0.05} = 1.645
$$

assegurando que $95\%$ das flutuações renováveis sejam cobertas pela reserva alocada.

### Condições KKT e Preço-Sombra

O Lagrangiano do problema é:

$$
\mathcal{L} = J_s + \sum_t \lambda_{t,s}\left[D_{net,t,s} + c_{t,s} - u_{t,s} - g^s_{t,s} - g^r_{t,s} - g^h_{t,s} - g^n_{t,s}\right] + \ldots
$$

As condições de Karush-Kuhn-Tucker (KKT) para a restrição de balanço produzem o **preço-sombra** $\lambda_{t,s}$, que é o custo marginal de servir uma unidade adicional de demanda:

$$
\lambda_{t,s} = \frac{\partial J_s^*}{\partial D_{net,t,s}}
$$

A análise de complementaridade revela três regimes operacionais:

**Regime 1 — Déficit** ($u_{t,s} > 0$, $c_{t,s} = 0$): O sistema está em escassez; mesmo com geração no máximo, a demanda não é atendida. A complementaridade de $u_{t,s}$ impõe $\lambda_{t,s} = \pi_u$. O preço-sombra bate no teto do custo de déficit.

**Regime 2 — Curtailment** ($c_{t,s} > 0$, $u_{t,s} = 0$): Há excesso de geração renovável que precisa ser descartada. A complementaridade de $c_{t,s}$ impõe $\lambda_{t,s} = -\pi_c$. O preço se torna negativo, sinalizando abundância.

**Regime 3 — Interior** ($u_{t,s} = 0$, $c_{t,s} = 0$): Operação normal sem slacks ativos. A condição de primeira ordem para $g^n_{t,s}$ fornece:

$$
\lambda_{t,s} = c_1 + c_2 g^n_{t,s} + \gamma\left(2g^n_{t,s} - g^n_{t-1,s} - g^n_{t+1,s}\right) + \overline{\mu}_{t,s} - \underline{\mu}_{t,s}
$$

Quando as restrições de capacidade não estão ativas ($\overline{\mu} = \underline{\mu} = 0$) e os efeitos de rampa são pequenos, obtém-se a relação clássica de custo marginal:

$$
\lambda_{t,s} \approx c_1 + c_2 g^n_{t,s}
$$

### Balanço Energético e Demanda Residual

A equação de balanço de potência pode ser reescrita para isolar a demanda residual — a quantidade que as fontes controláveis devem suprir:

$$
g^n_{t,s} = \underbrace{D_{net,t,s}}_{\text{demanda líquida}} - \underbrace{(g^s_{t,s} + g^r_{t,s} + g^h_{t,s})}_{\text{renováveis + hidro}} + c_{t,s} - u_{t,s}
$$

Definindo a demanda residual:

$$
D_{res,t,s} \equiv D_{net,t,s} - (g^s_{t,s} + g^r_{t,s} + g^h_{t,s})
$$

temos $g^n_{t,s} = D_{res,t,s} + c_{t,s} - u_{t,s}$. Nos períodos sem curtailment nem déficit, $g^n_{t,s} = D_{res,t,s}$.

### Interconexão entre Subsistemas

O fluxo de interconexão $x_{int,t,s}$ é calculado a partir dos dados de intercâmbio bilateral. Para um fluxo registrado do subsistema de origem $o$ para o destino $d$ com magnitude $f_{t,o \to d}$:

$$
x_{int,t,s} = \sum_{j: d_j = s} f_{t,j} - \sum_{i: o_i = s} f_{t,i}
$$

Isto é, o intercâmbio líquido de $s$ é a soma das importações menos a soma das exportações. A demanda líquida $D_{net} = D - x_{int}$ é a demanda "vista" pelos geradores locais.

### Custo Marginal de Operação (CMO)

O CMO oficial do SIN é calculado pelo modelo DECOMP do ONS e publicado com duas granularidades:

**CMO Semanal** (por patamar de carga): Para cada semana operativa $w$, subsistema $s$ e patamar $p \in \{\text{leve}, \text{médio}, \text{pesado}\}$:

$$
\text{CMO}_{w,s,p} \in \mathbb{R}_+ \quad \text{(R\$/MWh)}
$$

O spread entre patamares $\Delta_{w,s} = \text{CMO}_{w,s,\text{pesada}} - \text{CMO}_{w,s,\text{leve}}$ reflete a inflexibilidade do sistema para responder a variações de carga.

**CMO Semi-horário**: Com resolução $\Delta t = 30$ min, reflete o preço de curto prazo:

$$
\text{CMO}_{t,s} \in \mathbb{R}_+ \quad \text{para } t = 0, \tfrac{1}{2}, 1, \tfrac{3}{2}, \ldots \text{ (horas)}
$$

Para alinhamento com o painel horário, agregamos por média aritmética:

$$
\text{CMO}^h_{t,s} = \frac{1}{2}\left(\text{CMO}_{t,s} + \text{CMO}_{t+\frac{1}{2},s}\right)
$$

### Análise de Estacionariedade

A estacionariedade do CMO é avaliada pelo teste Augmented Dickey-Fuller (ADF). Sob $H_0$ a série possui raiz unitária (é não-estacionária):

$$
\Delta y_t = \alpha + \beta t + \phi y_{t-1} + \sum_{j=1}^{p} \delta_j \Delta y_{t-j} + \varepsilon_t
$$

Rejeitamos $H_0$ quando a estatística ADF é suficientemente negativa (p-valor < 0.05), indicando que a série de preços é estacionária — propriedade relevante para modelagem e previsão.

### Classificação de Regimes de Preço

O CMO horário é classificado em regimes utilizando os percentis $Q_{0.25}$ e $Q_{0.75}$ da distribuição empírica:

$$
\text{regime}(t,s) = \begin{cases}
\text{baixo} & \text{se } \text{CMO}^h_{t,s} \leq Q_{0.25} \\
\text{alto} & \text{se } \text{CMO}^h_{t,s} \geq Q_{0.75} \\
\text{normal} & \text{caso contrário}
\end{cases}
$$

### Métricas de Validação

A qualidade do modelo é avaliada comparando o despacho controlável otimizado $\hat{g}^n_{t,s}$ com o observado $g^n_{obs,t,s}$:

$$
\text{MAE}_s = \frac{1}{T}\sum_{t} \left|\hat{g}^n_{t,s} - g^n_{obs,t,s}\right|
$$

$$
\text{RMSE}_s = \sqrt{\frac{1}{T}\sum_{t} \left(\hat{g}^n_{t,s} - g^n_{obs,t,s}\right)^2}
$$

$$
R^2_s = 1 - \frac{\sum_t (\hat{g}^n_{t,s} - g^n_{obs,t,s})^2}{\sum_t (g^n_{obs,t,s} - \bar{g}^n_{obs,s})^2}
$$

A comparação entre os cenários com e sem intercâmbio permite isolar o efeito causal da interconexão: qualquer diferença $\Delta\text{RMSE}_s = \text{RMSE}_s^{(\text{sem interc})} - \text{RMSE}_s^{(\text{com interc})}$ é inteiramente atribuível aos fluxos $x_{int}$.

---

## Estrutura do Projeto

```
Mean-Field-Games-in-Renewable-Energy/
│
├── README.md
├── LICENSE
│
├── model1/                                    # Estudos com dados sintéticos
│   ├── initial_study.ipynb                   # Exploração inicial
│   ├── base_model.ipynb                      # Modelo sintético completo
│   └── data/                                  # Resultados de sweeps e Monte Carlo
│
├── validate_model/
│   ├── __init__.py
│   ├── pipeline.py                            # Módulo compartilhado (loaders, painel, visualizações)
│   │
│   ├── data/
│   │   ├── demanda_efetiva/                  # Curvas de carga reais por subsistema
│   │   │   ├── CURVA_CARGA_NORTE_2025.csv
│   │   │   ├── CURVA_CARGA_NORDESTE_2025.csv
│   │   │   ├── CURVA_CARGA_SUDESTE_2025.csv
│   │   │   └── CURVA_CARGA_SUL_2025.csv
│   │   │
│   │   ├── demanda_esperada/                 # Previsões dia-seguinte por patamar
│   │   │   ├── DEMANDA_N_2025.csv
│   │   │   ├── DEMANDA_NE_2025.csv
│   │   │   ├── DEMANDA_SE_2025.csv
│   │   │   └── DEMANDA_S_2025.csv
│   │   │
│   │   ├── producao_solar/                   # Geração solar fotovoltaica
│   │   │   └── fotovoltaicas_2025.csv
│   │   │
│   │   ├── producao_eolica/                  # Geração eólica
│   │   │   └── eolicas_2025.csv
│   │   │
│   │   ├── producao_hidroelétrica/           # Geração hidrelétrica (~60% da matriz)
│   │   │   └── hidroeletricas_2025.csv
│   │   │
│   │   ├── producao_non_renewable/           # Geração controlável
│   │   │   ├── nuclear_2025.csv
│   │   │   └── TERMICAS_2025.csv
│   │   │
│   │   ├── intercambio/                      # Fluxos de potência
│   │   │   ├── Intercambio_do_SIN_2025.csv          # Externo (internacional)
│   │   │   └── intercambio_interno_2025.csv         # Interno (entre subsistemas)
│   │   │
│   │   └── precos/                           # Custo Marginal de Operação (CMO)
│   │       ├── cmo_semanal2025.csv                  # DECOMP semanal (por patamar)
│   │       └── cmo_semihorario2025.csv              # Semi-horário (30 min)
│   │
│   └── outputs/                              # Predições do modelo
│       ├── predictions_social_dispatch_2025.csv
│       └── predictions_social_dispatch_2025_with_intercambio.csv
│
├── validation.ipynb                           # Notebook original de validação
├── 01_analise_exploratoria.ipynb              # Análise exploratória de dados
├── 02_modelos_despacho.ipynb                  # Modelos de despacho social
└── 03_analise_precos.ipynb                    # Análise de preços (CMO)
```

## Notebooks

### `01_analise_exploratoria.ipynb` — Análise Exploratória

Primeiro contato com os dados operacionais do SIN. O notebook carrega todas as fontes de dados através do módulo `pipeline.py` e constrói o painel unificado (hora × subsistema). Em seguida, realiza a exploração visual e estatística: estatísticas descritivas por subsistema, composição percentual da matriz energética, heatmaps dia × hora revelando padrões diurnos e sazonais de cada variável, análise detalhada da geração hidrelétrica (participação, perfil horário, sazonalidade, variabilidade), curvas de duração com geração empilhada, perfis diários médios, matrizes de correlação entre fontes, distribuição da demanda residual $D_{res} = D_{net} - (g^s + g^r + g^h)$ e balanço energético agregado.

### `02_modelos_despacho.ipynb` — Modelos de Despacho Social

Implementação e avaliação do modelo de otimização. Define a formulação do problema (parâmetros, estimação de $K_n$ e $\rho$ via quantis, incerteza renovável $\sigma_\epsilon$ por mês/hora), implementa o solver via CVXPY com fallbacks (OSQP → SCS → ECOS → heurística) e executa dois cenários: despacho com intercâmbio ($D_{net}$) e sem intercâmbio ($D$). Compara quantitativamente o impacto da interconexão através de métricas (MAE, RMSE, $R^2$), analisa preços-sombra $\lambda_{t,s}$ por regime operacional (normal, curtailment, déficit), visualiza séries temporais e produz scatter plots de observado vs. predito.

### `03_analise_precos.ipynb` — Análise de Preços (CMO)

Análise focada no Custo Marginal de Operação calculado pelo modelo DECOMP do ONS. Utiliza dois datasets — semanal (com discriminação por patamar de carga: leve, médio, pesado) e semi-horário (30 min) — para explorar: estatísticas descritivas e distribuições, evolução temporal e perfil intra-diário, testes de estacionariedade (ADF) sob $H_0$ de raiz unitária, classificação de regimes de preço via percentis, padrões mensais típicos de custo (boxplots sazonais, spreads $\Delta_{w,s}$ entre patamares) e correlação do CMO com variáveis físicas do sistema (demanda, geração térmica, demanda residual, mix de geração).

### `validation.ipynb` — Notebook Original

Notebook monolítico original que combina carregamento de dados, análise exploratória e modelagem. Mantido por referência — as funcionalidades foram reorganizadas nos três notebooks acima com melhorias e correções.

## Módulo Compartilhado: `validate_model/pipeline.py`

O módulo `pipeline.py` centraliza toda a infraestrutura de dados compartilhada pelos notebooks:

- **Leitura robusta de CSV**: detecção automática de separador e encoding, tratamento de formato numérico brasileiro (vírgula como decimal), correção de artefatos de encoding em strings
- **Normalização**: nomes de colunas padronizados, mapeamento de subsistemas (NORTE → N, NORDESTE → NE, etc.)
- **Loaders especializados**: funções para cada tipo de dado (demanda efetiva, demanda prevista, geração por fonte, intercâmbio interno/externo, CMO semanal e semi-horário)
- **Classe `SINPaths`**: encapsula todos os caminhos de dados para um dado ano, facilitando reutilização
- **Função `build_panel()`**: constrói o painel unificado (hora × subsistema) com uma única chamada
- **Visualizações**: heatmaps, curvas de duração e janelas temporais reutilizáveis

## Como Começar

### Pré-requisitos

```bash
pip install numpy pandas matplotlib seaborn scipy cvxpy scikit-learn statsmodels
```

### Instalação

```bash
git clone https://github.com/yourusername/Mean-Field-Games-in-Renewable-Energy.git
cd Mean-Field-Games-in-Renewable-Energy
```

### Dados

Os dados operacionais do SIN podem ser obtidos em:
- **Link**: https://drive.google.com/drive/folders/1mhDPyeKm5SD1Ba0SKNl9hDBhWuHTnN68?usp=drive_link
- Coloque os CSVs nos subdiretórios correspondentes dentro de `validate_model/data/`

### Execução

Os notebooks devem ser executados na seguinte ordem:

1. `01_analise_exploratoria.ipynb` — entendimento dos dados
2. `02_modelos_despacho.ipynb` — otimização e validação
3. `03_analise_precos.ipynb` — análise de custos marginais
