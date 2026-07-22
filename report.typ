
#set page(
  paper: "a4",
  margin: (x: 2.05cm, y: 2.00cm),
  numbering: "1",
)
#set text(lang: "pt", region: "BR", size: 10.05pt)
#set par(justify: true, leading: 0.585em)
#set heading(numbering: "1.")
#set enum(numbering: "(i)")
#set list(indent: 1.15em, body-indent: 0.55em)

#show heading.where(level: 1): it => block(above: 1.30em, below: 0.55em)[#it]
#show heading.where(level: 2): it => block(above: 1.00em, below: 0.36em)[#it]
#show heading.where(level: 3): it => block(above: 0.80em, below: 0.26em)[#it]

#let argmin = math.op("argmin")
#let argmax = math.op("argmax")
#let diag = math.op("diag")
#let grav = math.op("grav")
#let Up = math.op("Up")
#let Rev = math.op("Rev")
#let CMO = math.op("CMO")
#let PLD = math.op("PLD")
#let BR = math.op("BR")
#let pos = math.op("pos")
#let softmax = math.op("softmax")
#let KL = math.op("KL")
#let EE = math.op("E")
#let ind = math.op("1")

#let keybox(title, body) = block(
  inset: 8pt,
  stroke: 0.62pt + gray,
  radius: 5pt,
  fill: luma(249),
)[
  #strong(title)
  #v(0.26em)
  #body
]

#let warnbox(title, body) = block(
  inset: 8pt,
  stroke: 0.72pt + rgb("8a5a00"),
  radius: 5pt,
  fill: rgb("fff8e6"),
)[
  #strong(title)
  #v(0.26em)
  #body
]

#let resultbox(title, body) = block(
  inset: 8pt,
  stroke: 0.68pt + rgb("2f5f7f"),
  radius: 5pt,
  fill: rgb("f4f9fc"),
)[
  #strong(title)
  #v(0.26em)
  #body
]

#let imgfig(title, path, note: none) = figure(
  block(width: 100%)[
    #image(path, width: 100%)
    #if note != none [
      #v(0.18em)
      #text(size: 8.6pt, fill: gray.darken(30%))[#note]
    ]
  ],
  caption: title,
)

#align(center)[
  #text(17.2pt, weight: "bold")[MFG - Renewable Energy]

  #v(0.28em)


  #v(0.62em)
  #text(10.1pt)[João Felipe Vilas Boas]

  #v(0.18em)

]

#v(0.80em)

#block(inset: 9pt, stroke: 0.65pt + gray, radius: 5pt, fill: luma(250))[
  *Resumo.* Este documento tenta desenvolver uma estrutura empírico-matemática para estudar como a expansão de fontes renováveis intermitentes altera o equilíbrio de mercado, o despacho hidrotérmico, o valor econômico da água, o curtailment, as rampas líquidas e os incentivos privados de investimento no setor elétrico brasileiro. A implementação computacional reconstrói dados horários de 2025 por subsistema, com cobertura comum de 1º de janeiro a 31 de outubro, separando carga observada da Rede Básica, carga global, geração solar centralizada, MMGD, solar total, geração eólica, CMO observado, proxy de PLD e preço-sombra endógeno do modelo. A metodologia possui três camadas. A primeira é exploratória e documenta fatos estilizados: curva do pato, fatores de capacidade, rampas, captura de preço, diferenças sazonais e heterogeneidade por subsistema. A segunda é um modelo situacional finito com dois produtores solares, uma térmica marginal e uma hidrelétrica com reservatório proxy, usado para isolar canibalização de receita, valor da água, déficit horário, curtailment e entrada de agente solar. A terceira é um Mean Field Game locacional discreto, no qual agentes heterogêneos escolhem capacidade e localização tomando preços como dados, enquanto a distribuição populacional determina oferta agregada, preços locacionais, déficit, curtailment e incentivos de entrada. Os resultados preliminares indicam forte canibalização solar nos fatores de captura, heterogeneidade regional relevante e divergência expressiva entre capacidade privada e capacidade social ótima. A etapa MFG, porém, ainda deve ser interpretada como aproximação iterativa, pois a tolerância de convergência não foi atingida até o momento.
]

#v(0.42em)

= Introdução

A expansão de fontes renováveis intermitentes muda a pergunta econômica central da operação elétrica. Em um sistema com baixa penetração renovável, a questão natural é quanto de geração convencional cara pode ser substituída por geração solar ou eólica de baixo custo marginal. Em um sistema com alta penetração renovável, a pergunta fica mais sutil: não basta saber se existe energia renovável ao longo do dia, é preciso saber se essa energia aparece no local correto, na hora correta e com flexibilidade suficiente para acompanhar a carga residual.

Esse ponto é crucial para o Brasil. O sistema brasileiro combina uma base hidrelétrica relevante, uma expansão rápida da solar centralizada e da MMGD, grande concentração eólica em alguns subsistemas, transmissão limitada entre regiões e formação de preço de curto prazo que não pode ser confundida com o preço-sombra de um modelo acadêmico. A hidrelétrica, principal fonte de carga no sistema brasileiro, é um ativo intertemporal: turbinar água hoje pode reduzir custo hoje, mas reduz flexibilidade amanhã. Da mesma forma, a solar não deve ser tratada apenas como energia barata. Quando muitos produtores solares geram nas mesmas horas, a própria solar reduz a carga líquida exatamente quando vende energia, pressionando para baixo o preço capturado.

A pergunta principal que tentamos responder é:

#align(center)[
  #box(inset: 8pt, stroke: 0.65pt + black, radius: 4pt)[
    *Como o aumento da participação de fontes renováveis intermitentes altera o despacho hidrotérmico, o preço de curto prazo, e os incentivos privados de investimento em um mercado locacional de energia?*
  ]
]

Primeiro, organizamos os dados horários. Depois, transformamos esses dados em diagnósticos operacionais: carga líquida, carga residual, rampas, preço médio, preço capturado e fatores de captura. Em seguida, resolvemos um despacho centralizado para entender o papel do preço-sombra. Depois, usamos um laboratório situacional finito para explicar o mecanismo econômico de forma transparente: duas solares, uma térmica e uma hidrelétrica. Por fim, passamos para um Mean Field Game locacional, no qual uma população grande de agentes pequenos escolhe capacidade e localização respondendo a preços endógenos.

== Contribuições

Para iniciar, formalizamos uma ponte entre dados horários, despacho hidrotérmico e modelos de equilíbrio com muitos agentes. Segundo, explicitamos a diferença entre objetos que frequentemente são misturados: carga e demanda, CMO e PLD, preço observado e preço-sombra, solar centralizada e MMGD, carga da Rede Básica e carga global. Terceiro, propõe uma implementação discreta de MFG locacional: estado individual de localização, capacidade e tipo; ação de ajuste de capacidade e realocação; distribuição populacional; best response; atualização forward; clearing locacional; e iteração de ponto fixo.

== Resultado quantitativo em uma frase

#resultbox([Síntese dos resultados preliminares])[
  Os dados cobrem o período comum de 2025-01-01 a 2025-10-31. A capacidade solar efetiva total proxy é 37,2 GW, distribuída em 2,30 GW no Norte, 11,40 GW no Nordeste, 17,97 GW no Sudeste/Centro-Oeste e 5,54 GW no Sul. Os fatores de captura solares agregados usando CMO observado são inferiores a um em todos os subsistemas: 0,671 no Norte, 0,575 no Nordeste, 0,706 no Sudeste/Centro-Oeste e 0,703 no Sul. Na comparação de investimento, encontramos capacidade social ótima de 11,2 GW, capacidade observada/proxy de 37,2 GW e maior capacidade privada com margem não negativa de 59,5 GW. A iteração MFG termina com 51,4 GW de capacidade agregada aproximada, mas sem convergência plena: erro final de preço 2176,537, erro final de distribuição 0,078, tolerância 0,005 e 25 iterações.
]

= Notação, dados e objetos empíricos

A primeira exigência de rigor é separar notação física, notação de mercado e notação de modelo. Nesta versão, a letra $L$ é usada para carga física observada ou reconstruída. A letra $D$ fica reservada para demanda econômica, quando houver uma curva de demanda. 

== Índices e conjuntos


#keybox("Dicionário operacional das variáveis principais")[
  A notação abaixo será usada de forma consistente no texto. Ela é importante porque evita confundir objetos fisicamente parecidos, mas economicamente distintos.

  #table(
    columns: (0.70fr, 1.65fr, 2.25fr),
    inset: 4pt,
    align: left,
    [Símbolo], [Objeto de referência (notebook)], [Interpretação],
    [$L_(ell,t)^"RB"$], [`L_RB`], [carga observada pela Rede Básica no subsistema $ell$ e hora $t$],
    [$G_(ell,t)^"MMGD"$], [`G_s_MMGD`], [geração distribuída solar estimada, majoritariamente não despachável pelo operador],
    [$L_(ell,t)^"global"$], [`L_global`], [carga física recomposta antes do abatimento da MMGD],
    [$G_(ell,t)^"s,cen"$], [`G_s_centralizada`], [geração solar centralizada observada ou inferida],
    [$G_(ell,t)^"s,all"$], [`G_s_all`], [solar centralizada mais MMGD],
    [$G_(ell,t)^"w"$], [`G_w_total`], [geração eólica agregada],
    [$L_(ell,t)^"res"$], [`L_res_global` ou `L_res_grid`], [carga residual usada para diagnóstico de flexibilidade],
    [$CMO_(ell,t)^"obs"$], [`CMO_obs`], [custo marginal de operação observado/divulgado],
    [$hat(PLD)_(ell,t)$], [`PLD_hat`], [proxy de PLD construída a partir do CMO e limites regulatórios],
    [$lambda_(ell,t)^"model"$], [`lambda_model`], [preço-sombra endógeno do balanço no problema de otimização],
  )
]


O tempo é discreto e horário:

$ t = 1, dots, T. $

Os subsistemas são

$ cal(L) = ("N", "NE", "SE", "S"). $

O índice $ell in cal(L)$ representa localização. O índice $i$ representa produtor solar no modelo situacional. No MFG, o estado individual é $x$, a ação é $a$, a distribuição populacional é $m$ e o preço-sombra locacional é $lambda_(ell,t)$.

== Carga da Rede Básica, carga global e MMGD

Usamos a seguinte decomposição:

#table(
  columns: (1.22fr, 1.02fr, 2.80fr),
  inset: 4pt,
  align: left,
  [`L_RB`], [$L_(ell,t)^"RB"$], [Carga vista pela Rede Básica. A MMGD aparece, em grande parte, como abatimento dessa carga.],
  [`L_global`], [$L_(ell,t)^"global"$], [Carga global reconstruída. No notebook, $L^"global" = L^"RB" + G^"MMGD"$.],
  [`G_s_centralizada`], [$G_(ell,t)^"s,cen"$], [Geração solar centralizada conectada ao sistema.],
  [`G_s_MMGD`], [$G_(ell,t)^"MMGD"$], [Geração solar distribuída estimada atrás do medidor.],
  [`G_s_all`], [$G_(ell,t)^"s,all"$], [Solar total: centralizada mais MMGD.],
  [`L_net_s_grid`], [$L_(ell,t)^"net,grid"$], [Carga da Rede Básica líquida de solar centralizada.],
  [`L_res_grid`], [$L_(ell,t)^"res,grid"$], [Carga residual da Rede Básica após solar centralizada e eólica.],
  [`L_net_s_global`], [$L_(ell,t)^"net,global"$], [Carga global líquida de solar total.],
  [`L_res_global`], [$L_(ell,t)^"res,global"$], [Carga global residual após solar total e eólica.],
)

A ótica `grid` é mais próxima do problema operacional centralizado da Rede Básica. A ótica `global` é mais adequada para entender a física econômica ampla, pois recompõe a carga atendida pela MMGD. Manteremos como default a ótica global:

$ L^"default" = L^"global", quad G^"s,default" = G^"s,all", quad L^"res,default" = L^"res,global". $

Essa escolha é útil para estudar o impacto agregado da solar, mas há uma ressalva: a MMGD não é um recurso plenamente despachável pelo operador. Ela afeta a carga líquida e a formação de cenários, mas não será tratada como uma usina centralizada controlável.

== Cobertura e grandezas principais

#table(
  columns: (1.3fr, 1.25fr, 1.25fr, 0.80fr, 0.95fr),
  inset: 4pt,
  align: left,
  [Base], [Início], [Fim], [Obs.], [Missing],
  [Painel SIN], [2025-01-01 00:00], [2025-10-31 23:00], [29.184], [0,934%],
  [CMO horário], [2025-01-01 00:00], [2025-10-31 23:00], [29.088], [0,000%],
  [Hidrologia diária], [2025-01-01 00:00], [2025-10-31 00:00], [304], [0,000%],
  [Solar bruta], [2025-01-01 00:00], [2025-10-31 23:00], [841.440], [0,034%],
)

O período comum usado nos dados de referência é de 2025-01-01 00:00 a 2025-10-31 00:00. Como a base ainda não cobre todo o ano de 2025, as estatísticas devem ser interpretadas como resultados do período disponível, não como estatísticas anuais fechadas.

#table(
  columns: (0.62fr, 0.90fr, 0.90fr, 0.90fr, 0.90fr, 0.90fr, 0.90fr),
  inset: 4pt,
  align: right,
  [Subs.], [$L^"global"$ médio], [$L^"global"$ máx.], [$L^"res"$ médio], [$G^"s,cen"$], [$G^"MMGD"$], [$G^"s,all"$],
  [N], [8.789 MW], [12.204 MW], [8.019 MW], [8.742 MWh], [4.220.831 MWh], [4.229.573 MWh],
  [NE], [14.464 MW], [20.156 MW], [-1.181 MW], [14.130.877 MWh], [9.444.872 MWh], [23.575.749 MWh],
  [SE], [47.490 MW], [73.061 MW], [42.425 MW], [12.271.938 MWh], [24.740.992 MWh], [37.012.930 MWh],
  [S], [14.978 MW], [27.444 MW], [13.068 MW], [4.840 MWh], [8.809.161 MWh], [8.814.001 MWh],
)

No Nordeste, a carga residual média global aparece negativa no período disponível, porque a geração eólica e solar total excede a carga global média em várias horas. 

== CMO, PLD e preço-sombra

Há três preços diferentes. O primeiro é o CMO observado, denotado $CMO_(ell,t)^"obs"$. O segundo é o PLD, usado na liquidação de curto prazo; no momento ainda não puxamos dados de PLD observado, então seguimos com a implementação de `PLD_hat` como proxy. O terceiro é o preço-sombra do modelo, $lambda_(ell,t)^"model"$, multiplicador do balanço no problema de despacho.

#table(
  columns: (0.60fr, 0.95fr, 0.95fr, 0.95fr, 0.95fr, 0.95fr),
  inset: 4pt,
  align: right,
  [Subs.], [$CMO$ médio], [$CMO$ p95], [$CMO$ máx.], [`PLD_hat` médio], [`PLD_hat` máx.],
  [N], [139,84], [340,35], [668,1], [166,19], [668,1],
  [NE], [126,96], [333,64], [668,0], [156,53], [668,0],
  [SE], [202,75], [365,44], [1652,1], [215,29], [940,0],
  [S], [215,04], [368,74], [1653,2], [226,15], [940,0],
)

A diferença entre CMO e PLD proxy é fundamental. A proxy de PLD foi truncada por limites regulatórios ou operacionais; por isso, no Sudeste e no Sul, máximos de CMO acima de 1600 R\$/MWh aparecem como máximo de `PLD_hat` de 940 R\$/MWh.

= Diagnóstico operacional: carga residual, curva do pato e captura de preço

A carga líquida solar global é

$ L_(ell,t)^"net,global" = L_(ell,t)^"global" - G_(ell,t)^"s,all". $

A carga residual global é

$ L_(ell,t)^"res,global" = L_(ell,t)^"global" - G_(ell,t)^"s,all" - G_(ell,t)^"w". $

A rampa residual é

$ Delta L_(ell,t)^"res" = L_(ell,t+1)^"res" - L_(ell,t)^"res". $

Quando a solar cresce, o vale diurno da carga líquida aprofunda. A questão não é só a queda do nível médio, mas o aumento da inclinação de fim de tarde. Esse fenômeno é conhecido como curva do pato. Economicamente, a curva do pato significa que o sistema precisa reduzir geração controlável no meio do dia e subir rapidamente ao entardecer. Essa transição tem valor econômico mesmo quando a energia renovável total do dia é abundante.

#imgfig(
  [Figura 1 — Perfis horários por subsistema.],
  "outputs/figures/fig01_perfis_horarios.png",
  note: [A figura reconstrói, por subsistema, as variáveis que sustentam o restante do artigo: carga, solar, eólica, preços e hidrologia.]
)

#imgfig(
  [Figura 2 — Diagnóstico operacional de carga líquida e rampas.],
  "outputs/figures/fig02_diagnostico_operacional.png",
  note: [A figura evidencia que os problemas relevantes não são apenas de energia diária, mas de horário de mínima carga residual e magnitude da rampa.]
)

#imgfig(
  [Figura 3 — Curva do pato contrafactual por escala solar.],
  "outputs/figures/fig04_curva_pato_contrafactual.png",
  note: [A curva contrafactual mostra a passagem de um sistema em que a solar reduz carga líquida para outro em que a flexibilidade passa a ser a restrição marginal.]
)

#resultbox([Principais rampas observadas no diagnóstico consolidado])[
  No consolidado, a maior rampa residual positiva aparece no Sudeste/Centro-Oeste, com 8066,55 MW/h. O Nordeste também apresenta rampas relevantes, com rampa residual positiva máxima de 5154,81 MW/h e rampa líquida solar positiva máxima de 5630,36 MW/h. Esses números ajudam a explicar por que o problema não é apenas adicionar MWh renovável, mas garantir que recursos flexíveis acompanhem a transição horária da carga residual.
]

== Fator de captura e canibalização

Suponha que o preço horário seja uma função crescente da carga residual:

$ lambda_t = Lambda(L_t - G_t^"s" - G_t^"w"), quad Lambda' > 0. $

Para um produtor solar $i$,

$ g_(i,t)^"s" = K_i^"s" a_(i,t)^"s". $

A receita modelada é

$ Rev_i^"model" = sum_(t=1)^T lambda_t g_(i,t)^"s". $

A receita capturada por MWh é

$ R_i^"cap,model" = frac(sum_(t=1)^T lambda_t K_i^"s" a_(i,t)^"s", sum_(t=1)^T K_i^"s" a_(i,t)^"s") = frac(sum_(t=1)^T lambda_t a_(i,t)^"s", sum_(t=1)^T a_(i,t)^"s"). $

Se a capacidade agregada solar é $K^"agg"$ e a geração agregada é $G_t^"s" = K^"agg" a_t^"agg"$, então

$ frac(partial lambda_t, partial K^"agg") = - Lambda_t^' a_t^"agg". $

Logo,

$ frac(partial R_i^"cap,model", partial K^"agg") = - frac(sum_t a_(i,t)^"s" a_t^"agg" Lambda_t^', sum_t a_(i,t)^"s"). $

Se $Lambda_t^' > 0$ e há coincidência temporal entre o produtor $i$ e a solar agregada, o numerador é positivo e a derivada é negativa. Essa é a demonstração algébrica da canibalização, reduz-se o preço nas horas em que ela própria gera.

Empiricamente, quando se usa CMO ou PLD proxy, o fator de captura é

$ F^"capture" = frac(P^"solar", P^"avg"), $

com

$ P^"solar" = frac(sum_t P_t G_t^"s", sum_t G_t^"s"), quad P^"avg" = frac(1,T) sum_t P_t. $

#table(
  columns: (0.65fr, 0.90fr, 0.90fr, 0.90fr, 0.90fr),
  inset: 4pt,
  align: right,
  [Subs.], [$F^"solar"$ CMO], [$F^"wind"$ CMO], [$F^"solar"$ PLD_hat], [$F^"wind"$ PLD_hat],
  [N], [0,671], [1,356], [0,783], [1,245],
  [NE], [0,575], [1,220], [0,727], [1,146],
  [SE], [0,706], [1,011], [0,775], [0,998],
  [S], [0,703], [1,065], [0,772], [1,052],
)

#imgfig(
  [Figura 4 — Fatores de captura agregados por fonte e subsistema.],
  "outputs/figures/fig13_fator_captura_agregado.png",
  note: [A solar captura preço abaixo da média simples em todos os subsistemas. A eólica, por outro lado, captura preço acima ou próximo da média na maior parte dos subsistemas.]
)

#imgfig(
  [Figura 5 — Fator de captura por subsistema, estação, tipo de dia e métrica de preço.],
  "outputs/figures/fig05_fatores_captura_detalhado.png",
  note: [A figura mostra que a canibalização solar não é homogênea: ela muda por estação e por tipo de dia.]
)

#resultbox([Observado aqui])[
  A evidência preliminar é consistente com canibalização da receita solar. Usando CMO observado, o fator de captura solar fica entre 0,575 e 0,706, enquanto o fator de captura eólico fica entre 1,011 e 1,356. A leitura econômica é clara: a solar gera mais justamente quando o preço relativo tende a ser menor, enquanto a eólica, no período disponível, está mais alinhada a horas de maior valor ou menos sujeita à coincidência homogênea da solar.
]

= Hidrologia proxy e valor econômico da água

A hidrologia é tratada como proxy agregada por subsistema. Construímos um reservatório equivalente em energia para capturar o principal mecanismo econômico: água hoje versus água amanhã.

O estado hídrico agregado é $S_(ell,t)$. A dinâmica reduzida é

$ S_(ell,t+1) = S_(ell,t) + A_(ell,t) - alpha_ell^"h" h_(ell,t) - v_(ell,t). $

As restrições são

$ 0 <= S_(ell,t) <= overline(S)_ell, quad 0 <= h_(ell,t) <= overline(H)_ell, quad v_(ell,t) >= 0. $

O valor terminal da água é representado por

$ Phi(S_(ell,T+1)) = - theta_ell^"S" S_(ell,T+1). $

Como o problema é de minimização, o sinal negativo significa que terminar o horizonte com mais água reduz o valor objetivo, pois aumenta a flexibilidade futura.

#imgfig(
  [Figura 6 — Séries hidrológicas proxy por subsistema.],
  "outputs/figures/fig06_hidrologia_proxy.png",
  note: [A figura deve ser lida como diagnóstico econômico da água, não como simulação hidráulica física por usina.]
)

A versão física completa exigiria, para cada usina $j$,

$ P_(j,t)^"h" = eta_(j,t) rho^"w" grav H_(j,t)^"net" Q_(j,t)^"tur", $

com

$ H_(j,t)^"net" = Z_j^"mon"(V_(j,t)) - Z_j^"jus"(Q_(j,t)^"tur" + Q_(j,t)^"spill") - H_(j,t)^"loss". $

Essa extensão é relevante, mas estaríamos correndo o risco de fugir um pouco do escopo. Para o objetivo atual, a proxy é suficiente para testar se a solar preserva água, desloca térmica e altera preço-sombra.

= Despacho centralizado e preço-sombra

O despacho centralizado é o benchmark social de curto prazo. Ele responde à pergunta: dada a carga, a geração renovável disponível, a hidrologia proxy e a capacidade térmica, qual combinação de recursos minimiza o custo de atender o sistema?

Para um subsistema, o balanço é escrito como

$ L_t + c_t - u_t - G_t^"s" - G_t^"w" - h_t - n_t = 0. $

A função objetivo simplificada é

$ min_(n,h,c,u,S,v) sum_(t=1)^T [C^"T"(n_t) + pi_c c_t + pi_u u_t] - theta^"S" S_(T+1). $

O custo térmico é convexo:

$ C^"T"(n_t) = c_1 n_t + frac(c_2,2)n_t^2. $

Calibramos os principais parâmetros como:

#table(
  columns: (1.55fr, 1.00fr),
  inset: 4pt,
  align: left,
  [Parâmetro], [Valor],
  [$c_1$], [123,48 R\$/MWh],
  [$c_2$], [0,023424 R\$/MWh²],
  [$K^"TH"$ SIN], [13.676 MW],
  [$K^"H"$ SIN], [81.102 MW],
  [$K^"NUC"$ SIN], [2.016 MW],
  [Rampa térmica], [683 MW/h],
  [$pi_c$], [30 R\$/MWh],
  [$pi_u$ e VOLL], [3500 R\$/MWh],
  [Custo de oportunidade da água], [80 R\$/MWh],
  [Valor-sombra terminal da água proxy], [60 R\$/MWh],
)

== Lagrangiana e KKT

Introduza $lambda_t$ para o balanço e $mu_t$ para a dinâmica hídrica. A Lagrangiana, omitindo temporariamente limites de capacidade, é

$ cal(L) = sum_t [C^"T"(n_t) + pi_c c_t + pi_u u_t] - theta^"S" S_(T+1) + sum_t lambda_t (L_t + c_t - u_t - G_t^"s" - G_t^"w" - h_t - n_t) + sum_t mu_t (S_(t+1) - S_t - A_t + alpha^"h" h_t + v_t). $

Pelo teorema do envelope,

$ frac(partial Q^*, partial L_t) = lambda_t. $

Logo, $lambda_t$ é o custo marginal modelado de atender uma unidade adicional de carga naquela hora.

Para a térmica interior,

$ frac(partial cal(L), partial n_t) = C_T^"prime"(n_t) - lambda_t = 0, $

portanto

$ lambda_t = c_1 + c_2 n_t. $

Se houver penalidade quadrática de rampa (nossa assumption do modelo anterior)

$ frac(gamma,2)(n_t - n_(t-1))^2 + frac(gamma,2)(n_(t+1)-n_t)^2, $

sua derivada em $n_t$ é

$ gamma(n_t - n_(t-1)) - gamma(n_(t+1)-n_t) = gamma(2n_t - n_(t-1) - n_(t+1)). $

Então, em ponto interior,

$ lambda_t = c_1 + c_2 n_t + gamma(2n_t - n_(t-1) - n_(t+1)). $

Para curtailment,

$ frac(partial cal(L), partial c_t) = pi_c + lambda_t. $

Como $c_t >= 0$, se $c_t > 0$ temos $lambda_t = -pi_c$. Para déficit,

$ frac(partial cal(L), partial u_t) = pi_u - lambda_t. $

Como $u_t >= 0$, se $u_t > 0$ temos $lambda_t = pi_u$. Portanto,

$ -pi_c <= lambda_t <= pi_u. $

Essa banda é a ligação entre penalidades e preço-sombra. Em excesso renovável, o preço-sombra tende ao piso. Em escassez, tende ao teto.

Para a água, a derivada em relação a $h_t$ dá

$ -lambda_t + alpha^"h" mu_t = 0, $

ou

$ lambda_t = alpha^"h" mu_t. $

Portanto, a hidro gera até o ponto em que o valor da energia no balanço iguala o valor de oportunidade da água consumida.

#imgfig(
  [Figura 7 — Despacho com reservatório proxy: exemplo de dia típico no Sudeste/Centro-Oeste.],
  "outputs/figures/fig07_despacho_dia_tipico_SE.png",
  note: [A figura ilustra como o despacho usa térmica, hidro e armazenamento proxy para fechar a carga residual.]
)

#warnbox([Interpretação de custo negativo])[
  Quando o valor objetivo do despacho fica negativo, isso decorre do crédito terminal da água, $-theta^S S_(T+1)$. O resultado significa que o modelo atribui valor econômico a terminar o horizonte com água preservada.
]

= Modelo situacional finito

O modelo situacional é o ponto em que deixamos de ser apenas diagnóstico e passa a estudar interação estratégica simples. O objetivo é explicar, em um ambiente pequeno, os mecanismos que depois aparecem no MFG.

Há dois produtores solares, $i=1,2$, uma térmica e uma hidrelétrica. A geração solar é

$ g_(i,t)^"s" = K_i^"s" a_(i,t)^"s". $

A solar total é

$ G_t^"s" = g_(1,t)^"s" + g_(2,t)^"s". $

A escala do segundo produtor é

$ kappa = frac(K_2^"s", K_1^"s"). $

Aqui usamos um horizonte situacional de 192 horas, isto é, 8 blocos de 24 horas. A capacidade solar total de referência é 37,2 GW, a capacidade térmica situacional é 16,5 GW, a capacidade hidrelétrica é 85,7 GW e a nuclear inflexível é aproximadamente 2,0 GW. O fator solar médio no bloco situacional é 0,270 e o fator eólico médio é 0,544.

O lucro do agente solar $i$ é

$ Pi_i = sum_t lambda_t g_(i,t)^"s" - I_i(K_i^"s"), $

com custo de investimento

$ I_i(K_i^"s") = F_i + q_i K_i^"s" + frac(eta_i,2)(K_i^"s")^2. $

Quando $kappa=0$, o segundo agente não existe. Por isso, $Rev_2=0$ e $R_2^"cap"$ deve ser indefinido, não uma receita capturada zero. 

#imgfig(
  [Figura 8 — Sensibilidade situacional ao tamanho do segundo produtor solar.],
  "outputs/figures/fig09_situacional_sensibilidade.png",
  note: [A figura mostra como receita, preço capturado, déficit e curtailment respondem ao aumento de $kappa$.]
)

#resultbox([Modelo situacional])[
  O modelo situacional tem o objetivo de exemplificar um possível comportamento esperado, apesar de interpretação limitada dado que ainda estamos trabalhando com um modelo teórico e abstrato. Ele mostra que a entrada solar não é apenas adição de energia barata. Ela muda o perfil de carga residual, desloca térmica, pode preservar água, mas também reduz o preço-sombra em horas solares. Assim, a entrada do segundo agente pode reduzir a receita capturada do incumbente. O mecanismo é exatamente o que o MFG generaliza: cada agente é individualmente pequeno ou limitado, mas a massa agregada de decisões muda preços e restrições.
]

== Déficit horário versus energia agregada

Um resultado importante é que pode haver energia agregada suficiente e, ainda assim, déficit em algumas horas. Formalmente, a condição

$ sum_t (G_t^"s" + G_t^"w" + overline(N) + E^"hydro") >= sum_t L_t $

não implica

$ G_t^"s" + G_t^"w" + overline(N)_t + overline(H)_t >= L_t quad "para todo" t. $

A primeira é uma condição de energia agregada. A segunda é uma condição de potência e flexibilidade hora a hora. A diferença entre elas é exatamente a razão pela qual rampa, armazenamento e transmissão entram no problema.

= Sensibilidade hidrológica

A sensibilidade hidrológica compara cenários seco, base e úmido. O objetivo é saber se o valor da solar depende do estado da água. Os resultados agregados foram:

#table(
  columns: (1.00fr, 1.05fr, 1.05fr, 1.05fr, 0.90fr, 1.05fr, 1.05fr, 1.05fr),
  inset: 4pt,
  align: right,
  [Cenário], [Custo R\$ mi], [Térmica GWh], [Hidro GWh], [Spill GWh], [Curt. solar GWh], [Déficit MWh], [$lambda$ médio],
  [Seca -30%], [4899,51], [2667,59], [6110,35], [773,17], [0,14], [1.205.553], [2315,86],
  [Base], [949,75], [1427,06], [8388,35], [0,00], [0,00], [167.858], [671,45],
  [Úmida +30%], [152,55], [287,46], [9695,81], [1706,15], [0,01], [0], [142,76],
)

#imgfig(
  [Figura 9 — Sensibilidade hidrológica: seca, base e úmida.],
  "outputs/figures/fig18_valor_agua.png",
  note: [A figura completa os valores da tabela para base e úmida. O ponto conceitual é que o mesmo MWh solar tem valor diferente dependendo do estoque e da afluência.]
)

#resultbox([])[
  Em cenário seco, o preço-sombra médio e o déficit aumentam fortemente, e a energia solar passa a ter um valor sistêmico adicional: preservar água em horas de geração solar pode reduzir escassez em horas futuras. Em cenário úmido, o valor marginal da água tende a cair e pode aparecer maior risco de vertimento ou curtailment. Portanto, o valor da solar não é constante; ele depende do estado hidrológico.
]

= Investimento: ótimo social, orçamento fixo e livre entrada

A etapa de investimento conecta o despacho ao problema de entrada. Há três objetos diferentes: capacidade observada/proxy, capacidade privada e capacidade social.

No ótimo social, o planejador escolhe $K$ para minimizar custo total anualizado:

$ W(K) = I(K) + 365 C^"op,day"(K). $

Se houver vários dias típicos $omega$ com pesos $p_omega$,

$ W(K) = I(K) + 365 sum_(omega in Omega) p_omega C_omega^"op,day"(K). $

A condição de primeira ordem interior é

$ frac(partial I, partial K_ell^"s") + 365 sum_omega p_omega frac(partial C_omega^"op", partial K_ell^"s") = 0. $

Pelo envelope theorem,

$ frac(partial C_omega^"op", partial K_ell^"s") = - sum_t lambda_(ell,t,omega) a_(ell,t,omega)^"s". $

Logo,

$ frac(partial I, partial K_ell^"s") = 365 sum_omega p_omega sum_t lambda_(ell,t,omega) a_(ell,t,omega)^"s". $

O custo marginal anualizado de nova capacidade deve igualar o benefício marginal sistêmico esperado da geração solar.

No problema privado, o agente entra se

$ Pi_ell(K) = 365 sum_omega p_omega sum_t P_(ell,t,omega) K a_(ell,t,omega)^"s" - I_ell(K) >= 0. $

Se $P$ for PLD/proxy, a condição é uma regra de mercado. Se $P$ for $lambda^"model"$, é uma receita teórica do modelo.

#table(
  columns: (1.40fr, 1.00fr, 2.90fr),
  inset: 4pt,
  align: left,
  [Regime], [Capacidade], [Interpretação],
  [Orçamento fixo], [37,21 GW], [Capacidade total igual à observada/proxy, redistribuída por margem locacional.],
  [Livre entrada privada], [59,54 GW], [Maior capacidade na grade com margem anual não negativa.],
  [Ótimo social], [11,16 GW], [Capacidade que minimiza custo operacional mais investimento na grade.],
  [Observada/proxy], [37,21 GW], [Capacidade solar efetiva total por p99 da geração observada.],
)

A alocação com orçamento fixo encontrada foi: 8,95 GW no Norte, 9,98 GW no Nordeste, 10,42 GW no Sudeste/Centro-Oeste e 7,86 GW no Sul. Comparando com a capacidade observada/proxy, o modelo desloca capacidade para o Norte e para o Sul em relação à observação, e reduz a concentração relativa no Sudeste/Centro-Oeste.

#imgfig(
  [Figura 10 — Comparação entre capacidade privada, social e observada/proxy.],
  "outputs/figures/fig11_investimento_privado_social.png",
  note: [A diferença entre 59,54 GW privado e 11,16 GW social deve ser lida com cautela, pois depende da escala temporal, da proxy de preço e da penalização de curtailment.]
)

#resultbox([])[
  A divergência entre capacidade privada e social é o resultado mais importante da etapa de investimento. No nosso modelo, maior capacidade privada viável econtrada foi de 59,54 GW, enquanto o ótimo social na grade foi 11,16 GW. A interpretação é que a receita privada, sob a proxy adotada, ainda permite entrada em níveis nos quais o benefício sistêmico marginal já caiu. Isso é coerente com canibalização, curtailment, custo de conexão e externalidades locacionais. Como a versão atual mantém simplificações, o número absoluto deve ser tratado como preliminar; o mecanismo, porém, é central para o artigo.
]

= Mean Field Game locacional discreto


A razão para introduzir um MFG é que, em mercados com muitos projetos renováveis potenciais, a decisão relevante não é apenas a melhor resposta de um produtor isolado. Cada produtor é pequeno individualmente, mas a massa agregada de produtores altera os preços, a carga residual, o congestionamento e o curtailment. 

Em um jogo finito com $N$ agentes solares, o agente $i$ escolheria capacidade e localização resolvendo

$ max_(a_i) J_i(a_i, a_(-i)), $

onde $a_(-i)$ representa as ações de todos os outros agentes. Esse problema cresce rapidamente com $N$, porque cada agente precisaria antecipar a reação estratégica dos demais. A aproximação de campo médio substitui o vetor enorme $a_(-i)$ por uma distribuição agregada $m_t$ sobre tipos, capacidades e localizações. Assim, o problema individual passa a ser

$ max_(a) J(x,a; m, lambda), $

e a consistência exige que a distribuição gerada pelas políticas ótimas seja a própria distribuição usada para formar preços:

$ lambda_0 arrow.r pi_0^*(m,lambda_0) arrow.r G_0^"s"(m) arrow.r lambda_1(m; G_0) arrow.r pi_1^*(m,lambda_1) arrow.r dots $

Essa cadeia é o núcleo do MFG implementado que tentamos implementar. Criamos um laboratório matemático em que a intermitência renovável, a heterogeneidade locacional e os custos de entrada podem ser analisados em equilíbrio.

Do ponto de vista teórico, um MFG discreto possui quatro blocos:

#enum[
  um problema individual, no qual cada tipo escolhe ação ótima dado o preço esperado;
][
  uma equação forward, que transporta a distribuição de agentes conforme as ações escolhidas;
][
  um clearing de mercado, que transforma geração agregada em preços-sombra locacionais;
][
  uma condição de ponto fixo, exigindo consistência entre distribuição, política e preço.
]

No caso contínuo clássico, esses blocos aparecem como uma equação de Hamilton-Jacobi-Bellman para o valor individual e uma equação de Kolmogorov-Fokker-Planck para a distribuição. Usamos a versão discreta desses objetos: `best_response` substitui a HJB, `update_distribution` substitui a equação forward, e o despacho locacional substitui a condição de clearing, para essa versão locacional inicial.

== Estado, ação e distribuição

O estado individual é

$ x = (ell, b, theta), $

onde $ell$ é subsistema, $b$ é bin de capacidade e $theta$ é tipo. O conjunto de estados é

$ cal(X) = cal(L) times cal(B) times cal(Theta). $

A ação é

$ a = (Delta b, ell^"new"), $

isto é, ajuste discreto de capacidade e escolha de localização futura. A transição simplificada é

$ b^"new" = min(max(b + Delta b, 0), B), quad ell^"next" = ell^"new". $

A distribuição populacional é

$ m_t(x) = m_t(ell,b,theta), quad sum_(x in cal(X)) m_t(x)=1. $

A capacidade agregada em $ell$ é

$ K_(ell,t)^"agg"(m_t) = M sum_(b,theta) K_b m_t(ell,b,theta). $

Definimos a unidade de capacidade por bin é 93,0 MW por agente quando `N_AGENTS=100`. Os custos anualizados por região são: 514.922 R\$/MW/ano no Norte, 437.347 R\$/MW/ano no Nordeste, 446.031 R\$/MW/ano no Sudeste/Centro-Oeste e 475.324 R\$/MW/ano no Sul. Os limites máximos locacionais usados na regularização são 4,59 GW no Norte, 22,80 GW no Nordeste, 35,95 GW no Sudeste/Centro-Oeste e 11,08 GW no Sul.

== Problema individual

Dada uma trajetória de preços $lambda_(ell,t)$, cada agente resolve

$ V_t(x) = max_(a in cal(A)(x)) ( r_t(x,a;lambda,m) + beta EE( V_(t+1)(X_(t+1)) | x,a ) ). $

A recompensa é decomposta em receita, custo de investimento, custo de realocação e penalização de concentração:

$ r_t(x,a;lambda,m) = Rev_t(x,a;lambda) - I(x,a) - C^"move"(x,a) - C^"risk"(x,a,m). $

Se o agente escolhe localização $ell^"new"$ e capacidade $K_(b^"new")$,

$ Rev_t(x,a;lambda) = lambda_(ell^"new",t) K_(b^"new") a_(ell^"new",theta,t)^"s". $

O custo de ajuste é

$ I(x,a) = q_(ell^"new",theta) pos(K_(b^"new") - K_b)
  + frac(eta_(ell^"new",theta),2) pos(K_(b^"new") - K_b)^2
  + F_(ell^"new",theta) ind_("entrada"). $

O custo de realocação é

$ C^"move"(x,a) = rho_(ell,ell^"new") ind_(ell^"new" != ell). $

Usamos regularização com custo de ajuste `CHI_ADJ = 2000`, custo de realocação para outro subsistema de 500.000.000 e penalidade locacional `GAMMA_LOC = 5.000.000`. O objetivo dessa regularização é evitar deslocamentos explosivos da massa para uma única região por diferenças pequenas de preço.

== Equação forward e clearing

Dada a política ótima $pi_t(a|x)$,

$ m_(t+1)(x^"new") = sum_(x in cal(X)) sum_(a in cal(A)(x)) P(x^"new" | x,a) pi_t(a|x) m_t(x). $

A geração solar agregada é

$ G_(ell,t)^"s"(m_t) = M sum_(b,theta) K_b a_(ell,theta,t)^"s" m_t(ell,b,theta). $

O clearing locacional resolve

$ min_(n,h,c,u,F) sum_(ell in cal(L)) [C_ell^"T"(n_(ell,t)) + pi_c c_(ell,t) + pi_u u_(ell,t)] $

sujeito ao balanço

$ L_(ell,t) + c_(ell,t) - u_(ell,t) - G_(ell,t)^"s"(m_t) - G_(ell,t)^"w" - h_(ell,t) - n_(ell,t) - sum_(k in cal(L)) F_(k ell,t) + sum_(k in cal(L)) F_(ell k,t) = 0, $

e aos limites de transmissão

$ -overline(F)_(ell k) <= F_(ell k,t) <= overline(F)_(ell k). $

O multiplicador desse balanço é o preço locacional $lambda_(ell,t)$.

== Definição de equilíbrio

Um equilíbrio MFG é uma coleção

$ (pi_t, m_t, lambda_t, F_t, c_t, u_t)_(t=1)^T $

que satisfaz:

#enum[
  dado $lambda$, a política $pi$ resolve o problema individual;
  dado $pi$, a distribuição $m$ evolui pela equação forward;
  dado $m$, o clearing locacional fecha o balanço e gera $lambda$;
  fluxos, geração, curtailment, déficit e estoques respeitam restrições físicas;
  a trajetória de preços usada pelos agentes coincide com a trajetória gerada pelo clearing.
]

Em forma de ponto fixo,

$ lambda = cal(C)(cal(F)(BR(lambda))). $

O algoritmo usa relaxação:

$ lambda^(k+1) = (1 - omega_lambda) lambda^k + omega_lambda hat(lambda)^(k+1), $

$ m^(k+1) = (1 - omega_m) m^k + omega_m hat(m)^(k+1). $

Os erros são

$ e_lambda^k = max_(ell,t) |lambda_(ell,t)^(k+1) - lambda_(ell,t)^k|, $

$ e_m^k = max_x |m^(k+1)(x) - m^k(x)|. $

#imgfig(
  [Figura 11 — Convergência da iteração MFG.],
  "outputs/figures/fig12_mfg_convergencia.png",
  note: [O erro de preço não atinge a tolerância; os resultados do MFG devem ser interpretados como aproximações iterativas.]
)

#table(
  columns: (0.55fr, 0.85fr, 0.85fr, 0.85fr, 0.85fr, 0.95fr, 0.95fr, 0.95fr),
  inset: 4pt,
  align: right,
  [Subs.], [$K^"obs"$ MW], [$K^"MFG"$ MW], [$Delta K$], [$lambda$ médio], [$R^"cap"$], [Curt. MWh], [Déficit MWh],
  [N], [2297], [3261], [41,96%], [550,50], [1107,35], [0], [0],
  [NE], [11402], [15782], [38,41%], [550,50], [1042,62], [0], [0],
  [SE], [17975], [25534], [42,06%], [550,50], [1078,80], [0], [0],
  [S], [5539], [6833], [23,37%], [550,50], [1114,22], [0], [0],
)

#imgfig(
  [Figura 12 — Resultados agregados do MFG.],
  "outputs/figures/fig13_mfg_resultados_1.png",
  note: [A figura resume capacidade, preço, curtailment e déficit na aproximação MFG.]
)

#imgfig(
  [Figura 13 — Preços locacionais e distribuição do MFG.],
  "outputs/figures/fig13_mfg_precos.png",
  note: [A igualdade de preço médio reportada em 550,50 R\$/MWh sugere que a versão atual ainda não gerou separação locacional forte no resumo final.]
)

#imgfig(
  [Figura 14 — Distribuição populacional do MFG por localização, capacidade e tipo.],
  "outputs/figures/fig13_mfg_distribuicao.png",
  note: [A distribuição populacional é o objeto característico do MFG: não acompanhamos cada firma, mas a massa de agentes em cada estado.]
)

#warnbox([Status do MFG])[
  A iteração MFG discreta não atingiu a tolerância estabelecida. O erro final de preço foi 2176,537 (muito acima do aceitável), o erro final da distribuição foi 0,078, a tolerância era 0,005 e o algoritmo parou em 25 iterações. Portanto, os números acima não devem ser chamados de equilíbrio MFG plenamente convergido. Eles são uma aproximação iterativa útil para diagnóstico de direção econômica.
]

#bibliography("references.bib", title: "Referências")
