# Proposta de Modelagem para Geração de Energia

Este documento descreve uma proposta de modelagem matemática para a geração de energia, dividida por fontes (solar, outras renováveis e não renováveis). O objetivo é estabelecer uma base conceitual para análises e simulações de sistemas de energia.

---

## 1. Estrutura Fundamental do Modelo

A geração total de energia do sistema, $G_{total}(t)$, em um determinado instante $t$, pode ser expressa como a soma das contribuições das principais categorias de fontes:

$$G_{total}(t) = g_s(t) + g_r(t) + g_n(t)$$

Onde:
- **$g_s(t)$**: Geração de energia solar.
- **$g_r(t)$**: Geração de outras fontes renováveis (eólica, hídrica, etc.).
- **$g_n(t)$**: Geração de fontes não renováveis (termelétricas a gás, carvão, etc.).

A seguir, propomos uma modelagem básica para cada uma dessas funções.

---

## 2. Modelagem dos Componentes de Geração

### 2.1. Geração Solar: $g_s(t)$

A geração solar é altamente dependente de fatores astronômicos e atmosféricos.

**Parâmetros Chave:**
- `$A$`: Área total de superfície dos painéis solares [m²].
- `$\eta$`: Eficiência média de conversão dos painéis (e.g., 0.15-0.22).
- `$I(t, \lambda, \phi)$`: Irradiância solar [W/m²], que depende do tempo $t$ e da localização (latitude $\lambda$, longitude $\phi$).
- `$c(t)$`: Fator de nebulosidade (0 para céu claro, 1 para encoberto).
- `$P_{temp}(t)$`: Fator de perda por temperatura.

**Função de Geração Solar Proposta:**
$$g_s(t) = A \cdot \eta \cdot I(t, \lambda, \phi) \cdot (1 - c(t)) \cdot P_{temp}(t)$$

### 2.2. Geração de Outras Renováveis: $g_r(t)$

Esta categoria agrega fontes como eólica e hídrica.

#### Geração Eólica ($g_w(t)$)

A geração eólica é uma função cúbica da velocidade do vento $v(t)$, operando dentro de limites específicos.

**Função de Geração Eólica Proposta (simplificada):**
$$
g_w(t) \propto
\begin{cases}
    0 & \text{se } v(t) < v_{\text{cut-in}} \text{ ou } v(t) > v_{\text{cut-out}} \\
    P_{\text{rated}} \left( \frac{v(t)^3 - v_{\text{cut-in}}^3}{v_{\text{rated}}^3 - v_{\text{cut-in}}^3} \right) & \text{se } v_{\text{cut-in}} \leq v(t) < v_{\text{rated}} \\
    P_{\text{rated}} & \text{se } v_{\text{rated}} \leq v(t) \leq v_{\text{cut-out}}
\end{cases}
$$

#### Geração Hídrica ($g_h(t)$)

A geração hídrica depende do fluxo de água $Q(t)$ e da altura da queda $H$.

**Função de Geração Hídrica Proposta:**
$$g_h(t) = \eta_h \cdot \rho \cdot g \cdot H \cdot Q(t)$$

### 2.3. Geração Não Renovável: $g_n(t)$

A geração não renovável é **despachável**, ou seja, pode ser controlada ativamente.

**Parâmetros Chave:**
- `$C_n$`: Capacidade máxima de geração da planta [MW].
- `$\delta_n(t)$`: Fator de despacho/operação (parâmetro de controle entre 0 e 1).

**Função de Geração Não Renovável Proposta:**
$$g_n(t) = C_n \cdot \delta_n(t)$$

---

## 3. A Condição de Equilíbrio e a Dinâmica do Sistema

Para manter a estabilidade da rede, a geração total deve atender à demanda $D(t)$ a cada instante. A análise mais rica vem da dinâmica do sistema (suas taxas de variação).

**Condição de Equilíbrio Dinâmico:**
Derivando a equação de equilíbrio $G_{total}(t) = D(t)$, obtemos:
$$\frac{dG_{total}}{dt} = \frac{dD}{dt}$$

O que se expande para:
$$g'_s(t) + g'_r(t) + g'_n(t) = D'(t)$$

**Interpretação:**
A capacidade de rampa das fontes controláveis ($g'_n(t)$) deve ser, a todo momento, suficiente para compensar a variação da demanda ($D'(t)$) e a volatilidade das fontes renováveis intermitentes ($g'_s(t)$ e $g'_r(t)$). A principal equação operacional é:

$$g'_n(t) = D'(t) - g'_s(t) - g'_r(t)$$
