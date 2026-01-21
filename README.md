# Mean Field Games in Renewable Energy Systems

Here is an overview for an initial analysis. Modeling Brazil's energy consumption 

## Overview

This project ends is to explore the application of **Mean Field Game theory** to energy systems with significant renewable energy sources (solar, wind, and hydroelectric). The main objectives are as far as we have in these notebooks are:

1. **Model Development**: Create a deterministic dispatch optimization model that balances:
   - Intermittent renewable generation (solar \(g^s\), wind \(g^r\))
   - **Hydroelectric generation (\(g^h\))** - controllable renewable baseload
   - Thermal/nuclear generation (\(g^n\)) - fully controllable
   - Demand (\(D\))
   - Curtailment (\(c\)) and deficit (\(u\))

2. **Theoretical Analysis**: Study the system behavior through:
   - Parameter sensitivity analysis
   - Marginal value curves for capacity expansion
   - Shadow price (λ) dynamics and regime classification
   - Ramping constraints and reserve requirements

3. **Real-World Validation**: Test the model against actual data from Brazil's SIN, including:
   - 4 subsystems (North, Northeast, Southeast, South)
   - Hourly data for 2025 across **all major generation sources**
   - Internal interconnection flows
   - **Hydroelectric dispatch** (Brazil's dominant source, ~60% of electricity)
   - Observed vs. predicted thermal/nuclear dispatch comparison

## Mathematical Framework

### Social Planner Problem

The system minimizes total operational cost over horizon \(T\):

\[
\min_{g^n, c, u} \sum_{t=0}^{T-1} \left[ c_1 g^n_t + \frac{1}{2}c_2 (g^n_t)^2 + \pi_u u_t + \pi_c c_t + \frac{\gamma}{2}(g^n_t - g^n_{t-1})^2 \right]
\]

**Subject to:**
- **Energy balance**: \(g^s_t + g^r_t + g^h_t + g^n_t = D_t + c_t - u_t\)
- **Capacity limits**: \(0 \leq g^n_t \leq K_n\)
- **Ramping constraint**: \(|g^n_t - g^n_{t-1}| \leq \rho\)
- **Reserve requirement**: \(g^n_t + R_t \leq K_n\), where \(R_t = z_\alpha \sigma_\epsilon(t)\)
- **Non-negativity**: \(c_t, u_t \geq 0\)
- **Hydroelectric**: \(g^h_t\) is **exogenous** (observed, not optimized) in validation

**Parameters:**
- \(c_1, c_2\): thermal generation cost (linear and quadratic)
- \(\pi_u\): deficit penalty (Value of Lost Load)
- \(\pi_c\): curtailment penalty
- \(\gamma\): ramping cost coefficient
- \(K_n, K_s, K_r\): installed capacities (controllable, solar, wind)
- \(\rho\): maximum ramp rate

### Shadow Price (λ)

The marginal cost of serving demand at time \(t\):

\[
\lambda_t = \frac{\partial C^*}{\partial D_t}
\]

**Regime classification:**
- **Deficit**: \(u_t > 0 \Rightarrow \lambda_t \approx \pi_u\)
- **Curtailment**: \(c_t > 0 \Rightarrow \lambda_t \approx \pi_c\)
- **Interior**: Otherwise, \(\lambda_t \approx c_1 + c_2 g^n_t\) (+ ramping effects if \(\gamma > 0\))

## Project Structure

```
Mean-Field-Games-in-Renewable-Energy/
│
├── README.md                      # This file
├── LICENSE
│
├── model1/                        # Synthetic model studies
│   ├── initial_study.ipynb       # First exploration with synthetic data
│   ├── base_model.ipynb          # Main synthetic model with full analysis
│   └── data/                      # Results from parameter sweeps
│       ├── baseline_summary.csv
│       ├── sweep_*.csv           # Parameter sensitivity results
│       ├── mc_100_summary.csv    # Monte Carlo simulation results
│       ├── lhs_300_summary.csv   # Latin Hypercube Sampling results
│       └── mv_*.csv              # Marginal value curves
│
├── validation.ipynb               # Real-world data validation (SIN)
│
└── validate_model/
    ├── data/                      # Real data from Brazilian power system
    │   ├── demanda_efetiva/      # Actual demand (load curves)
    │   │   ├── CURVA_CARGA_NORDESTE_2025.csv
    │   │   ├── CURVA_CARGA_NORTE_2025.csv
    │   │   ├── CURVA_CARGA_SUDESTE_2025.csv
    │   │   └── CURVA_CARGA_SUL_2025.csv
    │   │
    │   ├── demanda_esperada/     # Expected demand (day-ahead forecasts)
    │   │   ├── DEMANDA_N_2025.csv
    │   │   ├── DEMANDA_NE_2025.csv
    │   │   ├── DEMANDA_S_2025.csv
    │   │   └── DEMANDA_SE_2025.csv
    │   │
    │   ├── producao_eolica/      # Wind generation
    │   │   └── eolicas_2025.csv
    │   │
    │   ├── producao_solar/       # Solar generation
    │   │   └── fotovoltaicas_2025.csv
    │   │
    │   ├── producao_hidroelétrica/ # Hydroelectric generation
    │   │   └── hidroeletricas_2025.csv   # Brazil's largest source (~60%)
    │   │
    │   ├── producao_non_renewable/ # Controllable generation
    │   │   ├── nuclear_2025.csv
    │   │   └── TERMICAS_2025.csv
    │   │
    │   └── intercambio/          # Interconnection flows
    │       ├── Intercambio_do_SIN_2025.csv       # External (international)
    │       └── intercambio_interno_2025.csv      # Internal (inter-subsystem, between Regions)
    │
    └── outputs/                   # Model predictions
        ├── predictions_social_dispatch_2025.csv
        └── predictions_social_dispatch_2025_with_intercambio.csv
```

## Getting Started

### Prerequisites

```bash
pip install numpy pandas matplotlib scipy cvxpy scikit-learn
```

### Quick Start

1. **Clone the repository**
```bash
git clone https://github.com/yourusername/Mean-Field-Games-in-Renewable-Energy.git
cd Mean-Field-Games-in-Renewable-Energy
```

2. **Download the data**
   - Access: https://drive.google.com/drive/folders/1mhDPyeKm5SD1Ba0SKNl9hDBhWuHTnN68?usp=drive_link
   - Place CSV files in the appropriate `validate_model/data/` subdirectories

3. **Run the notebooks**
   - Finally `validation.ipynb` for model analysis with data

## Notebooks Description

### 1. `validation.ipynb` - Real-World Application

**Purpose**: Validate the model using actual data from Brazil's SIN (2025).

**Workflow:**
1. **Data Loading & Cleaning**
   - Robust CSV parser (handles encoding, separators)
   - Normalization of subsystem names (N, NE, SE, S)
   - Datetime alignment and aggregation
   - **Load hydroelectric data** (Brazil's dominant source)

2. **Data Preparation**
   - Merge demand, generation (solar, wind, **hydro**, nuclear, thermal), and interconnection data
   - Compute net demand: \(D_{net} = D - x_{int}\)
   - Estimate renewable uncertainty: \(\sigma_\epsilon(t,s)\) by hour (solar + wind only)

3. **Hydroelectric Treatment**
   - Treat \(g^h_t\) as **exogenous controllable** (observed, not optimized)
   - Accounts for ~60% of Brazil's electricity
   - Complex reservoir dynamics beyond hourly dispatch scope
   - **Residual demand** for thermal: \(D_{net} - (g^s + g^r + g^h)\)

4. **Parameter Calibration**
   - Estimate \(K_n\) from 99.5th percentile of observed \(g^n\) (thermal/nuclear only)
   - Estimate \(\rho\) from 99.5th percentile of observed ramps
   - Grid search over (γ, c_2) to minimize RMSE

5. **Dispatch Optimization**
   - Solve social planner problem with CVXPY (OSQP/ECOS solvers)
   - **Balance**: \(g^s + g^r + g^h + g^n = D_{net} + c - u\) where \(g^h\) is fixed (observed)
   - Incorporate reserve constraints: \(g^n_t + R_t \leq K_n\)
   - Extract dual variables (shadow prices λ)

6. **Validation & Diagnostics**
   - **Metrics**: MAE, RMSE, R² for \(g^n_{hat}\) vs \(g^n_{obs}\)
   - **Heatmaps**: Day×Hour patterns for demand, solar, wind, **hydro**, interconnection
   - **Duration curves**: Load ordering showing **hydro baseload + intermittent renewables + thermal residual**
   - **Time series**: Model vs. observed dispatch windows with all generation sources
   - **Reserve slack**: Check for constraint violations
   - **Shadow price vs. marginal cost**: Verify KKT conditions

